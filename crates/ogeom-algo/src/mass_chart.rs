//! Mass properties of a trimmed face, integrated in its own chart.
//!
//! A face is a region of its surface's `(u, v)` chart bounded by its
//! pcurves. Green's theorem turns the integral over that region into one
//! round its boundary:
//!
//! ```text
//! ∬ f du dv = ∮ F dv,   F(u, v) = ∫ f(s, v) ds from u_ref to u
//! ```
//!
//! so each quadrature node on a boundary loop carries an inner integral
//! along `u` from a fixed `u_ref` to the node. Every sample of that inner
//! integral is a point of the surface with its `n dA` and a weight, which
//! is what the mass accumulators take: the region is never meshed and its
//! boundary is the exact pcurves, whatever surface and whatever trim.
//!
//! The surfaces taken are those whose integrands are trigonometric
//! polynomials or polynomials in the chart: the analytic ones, and lines
//! and conics swept or revolved. Panels break at a pcurve's knots and at
//! every quarter turn, where the ten-point Gauss rule integrates such an
//! integrand to rounding, and the whole is run again on panels twice as
//! fine until two runs agree to a part in ten billion. A short boundary
//! panel on an analytic surface takes as few points as its error bound
//! allows, and the inner integrals there, exact on quarter turns, are not
//! refined between runs.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{
    Curve, Curve2d as _, Curve3d as _, PlanarCurve, Surface as _, SurfaceGeometry,
    Transformable as _,
};
use ogeom_math::{Point, Point2, Vector, Vector2, gauss_legendre_rule};
use ogeom_topo::{EdgeRepr, Model, NodeData, Orientation, Shape};

/// One pcurve piece of a boundary loop, walked from `t0` to `t1` and moved
/// by `shift` in the chart: a whole number of periods where the walk
/// crossed a periodic surface's join.
struct Segment {
    edge: Shape,
    curve: PlanarCurve,
    t0: f64,
    t1: f64,
    shift: Vector2,
}

impl Segment {
    fn at(&self, t: f64, tol: Tolerances) -> OgeomResult<(Point2, Vector2)> {
        let d = self.curve.derivatives_at(t, 1, tol)?;
        Ok((Point2::ORIGIN + d[0] + self.shift, d[1]))
    }
}

/// A face as its surface and the closed chart loops that bound it, each
/// with the sign that makes its enclosed region count positive for the
/// face's outer boundary and negative for a hole.
pub(crate) struct ChartFace {
    surface: SurfaceGeometry,
    loops: Vec<(Vec<Segment>, f64)>,
    /// The face's orientation: `-1` where it uses its surface reversed.
    sign: f64,
    /// Where the inner integrals start.
    u_ref: f64,
    /// The chart's size, for weighing a miss against it.
    scale: f64,
    /// The surface's own knot lines in `u` and in `v`, where panels break:
    /// read once, not once per inner integral.
    knot_lines: (Vec<f64>, Vec<f64>),
    /// The straight chart steps closing each junction where one piece ends
    /// short of where the next begins, as their start, their step and their
    /// loop's sign. A fitted pcurve may end anywhere within its edge's
    /// tolerance of the vertex, and an open loop's integral depends on
    /// where the inner integrals start: a step of `dv` costs the strip from
    /// `u_ref` across, the chart's whole width on a face wrapping round its
    /// axis. Closed by these steps, every loop measures the region it
    /// bounds whatever `u_ref` is.
    bridges: Vec<(Point2, Vector2, f64)>,
}

/// A face's chart loops, or `None` where they cannot be had exactly: a
/// surface whose integrand is not a trigonometric polynomial, an edge
/// without a pcurve on the face or whose pcurve strays from it, a placement
/// that scales, or a loop whose pieces do not meet.
pub(crate) fn chart_face(model: &Model, face: &Shape, tol: Tolerances) -> Option<ChartFace> {
    loops_of(model, face, tol).ok().flatten()
}

/// Whether a surface's integrands are trigonometric polynomials in its
/// chart, or polynomials: the analytic surfaces, and a line or conic swept
/// or revolved. A spline's are neither, and a rule that integrates them to
/// rounding costs more than the mesh.
fn integrable(surface: &SurfaceGeometry) -> bool {
    fn conic(curve: &Curve) -> bool {
        match curve {
            Curve::Line(_)
            | Curve::Circle(_)
            | Curve::Ellipse(_)
            | Curve::Hyperbola(_)
            | Curve::Parabola(_) => true,
            Curve::Trimmed(trimmed) => conic(trimmed.basis()),
            _ => false,
        }
    }
    match surface {
        SurfaceGeometry::Plane(_)
        | SurfaceGeometry::Cylinder(_)
        | SurfaceGeometry::Cone(_)
        | SurfaceGeometry::Sphere(_)
        | SurfaceGeometry::Torus(_) => true,
        // A spline, and a spline curve swept or revolved, is a polynomial
        // (or rational) piece by knot span: taken with the panels broken
        // at its knots.
        SurfaceGeometry::BSpline(_) => true,
        SurfaceGeometry::Extrusion(e) => conic(e.curve()) || spline_knots(e.curve()).is_some(),
        SurfaceGeometry::Revolution(r) => conic(r.curve()) || spline_knots(r.curve()).is_some(),
        _ => false,
    }
}

/// A spline curve's distinct knots, through a trim.
fn spline_knots(curve: &Curve) -> Option<Vec<f64>> {
    match curve {
        Curve::BSpline(b) => Some(b.knots().distinct().into_iter().map(|(k, _)| k).collect()),
        Curve::Trimmed(t) => spline_knots(t.basis()),
        _ => None,
    }
}

/// Where a surface's integrand changes piece: its knot lines in `u` and
/// in `v`.
fn knot_lines(surface: &SurfaceGeometry) -> (Vec<f64>, Vec<f64>) {
    let distinct = |k: &ogeom_math::KnotVector| -> Vec<f64> {
        k.distinct().into_iter().map(|(x, _)| x).collect()
    };
    match surface {
        SurfaceGeometry::BSpline(b) => (distinct(b.u_knots()), distinct(b.v_knots())),
        SurfaceGeometry::Extrusion(e) => (spline_knots(e.curve()).unwrap_or_default(), Vec::new()),
        SurfaceGeometry::Revolution(r) => (Vec::new(), spline_knots(r.curve()).unwrap_or_default()),
        _ => (Vec::new(), Vec::new()),
    }
}

fn loops_of(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Option<ChartFace>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    if !integrable(surface) {
        return Ok(None);
    }
    let Some(walked) = walked(model, face, true, tol)? else {
        return Ok(None);
    };
    let ((u0, u1), (v0, v1)) = surface.domain();
    let u_ref = if surface.is_periodic_u() {
        walked.lo.x
    } else {
        walked.lo.x.clamp(u0, u1)
    };
    let knot_lines = knot_lines(&walked.placed);
    let period = Vector2::new(
        if surface.is_periodic_u() {
            u1 - u0
        } else {
            0.0
        },
        if surface.is_periodic_v() {
            v1 - v0
        } else {
            0.0
        },
    );
    let mut bridges = Vec::new();
    for (segments, region) in &walked.loops {
        for (k, segment) in segments.iter().enumerate() {
            let next = &segments[(k + 1) % segments.len()];
            let (end, _) = segment.at(segment.t1, tol)?;
            let (start, _) = next.at(next.t0, tol)?;
            // A loop wrapping round a periodic chart closes a whole period
            // from where it began; only what is left over is a step.
            let gap = start - end;
            let step = gap - fold(gap, period);
            if step.y != 0.0 {
                bridges.push((end, step, *region));
            }
        }
    }
    Ok(Some(ChartFace {
        knot_lines,
        bridges,
        surface: walked.placed,
        loops: walked.loops,
        sign: walked.handedness
            * if face.orientation() == Orientation::Reversed {
                -1.0
            } else {
                1.0
            },
        u_ref,
        scale: walked.scale,
    }))
}

/// `1` for a placement that keeps a chart's metric and handedness, `-1` for
/// one that keeps its metric and reflects it, `None` for one that scales:
/// the pcurves are the unplaced surface's, and a scale changes the metric
/// under them.
pub(crate) fn rigid_handedness(placement: &ogeom_math::Transform) -> Option<f64> {
    match placement.kind() {
        ogeom_math::TransformKind::Identity
        | ogeom_math::TransformKind::Translation
        | ogeom_math::TransformKind::Rotation => Some(1.0),
        ogeom_math::TransformKind::PlaneMirror | ogeom_math::TransformKind::PointMirror => {
            Some(-1.0)
        }
        _ => None,
    }
}

/// A face's boundary walked into closed chart loops.
struct Walked {
    /// The surface, placed as the face is.
    placed: SurfaceGeometry,
    /// `-1` where the placement reflects: the placed chart's own normal
    /// then points against the face's.
    handedness: f64,
    /// Each loop, with the sign that turns it round its region: `1` where
    /// the face lies to the left of the walk, `-1` where to the right.
    loops: Vec<(Vec<Segment>, f64)>,
    /// The chart's lower corner and its size.
    lo: Point2,
    scale: f64,
}

/// Points along a face's boundary, each with its edge and the chart
/// direction the face lies in from there, read off the walked loops'
/// windings. `None` where the boundary cannot be walked into closed loops.
pub(crate) fn material_sides(
    model: &Model,
    face: &Shape,
    tol: Tolerances,
) -> Option<Vec<(Shape, Point2, Vector2)>> {
    let walked = walked(model, face, false, tol).ok()??;
    let mut out = Vec::new();
    for (segments, region) in &walked.loops {
        for segment in segments {
            for k in 1..=4 {
                let t = segment.t0 + (segment.t1 - segment.t0) * f64::from(k) / 5.0;
                let (at, d) = segment.at(t, tol).ok()?;
                let heading = if segment.t1 < segment.t0 { -d } else { d };
                out.push((segment.edge.clone(), at, heading.perpendicular() * *region));
            }
        }
    }
    Some(out)
}

/// A face's loops in its chart, or `None` where they cannot be had: an
/// edge without a pcurve on the face, a placement that scales, or a loop
/// whose pieces do not meet; and where `strict`, a pcurve straying from
/// its edge.
fn walked(
    model: &Model,
    face: &Shape,
    strict: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Walked>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    // The pcurves are the unplaced surface's; a rigid placement keeps its
    // chart, a scaling one would change the metric under them.
    let placement = face.transform(model.datums())?;
    let Some(handedness) = rigid_handedness(&placement) else {
        return Ok(None);
    };
    let placed = surface.clone().transformed(&placement, tol)?;
    let ((u0, u1), (v0, v1)) = surface.domain();
    let period = Vector2::new(
        if surface.is_periodic_u() {
            u1 - u0
        } else {
            0.0
        },
        if surface.is_periodic_v() {
            v1 - v0
        } else {
            0.0
        },
    );

    let mut loops = Vec::new();
    for wire in model.ordered_children_of(face)? {
        let Some(segments) = walk(model, data.surface, &placed, &wire, period, strict, tol)? else {
            return Ok(None);
        };
        loops.push(segments);
    }
    if loops.is_empty() {
        return Ok(None);
    }

    // The chart's extent, from a few points of every piece.
    let (mut lo, mut hi) = (
        Point2::new(f64::INFINITY, f64::INFINITY),
        Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for segments in &loops {
        for segment in segments {
            for k in 0..=4 {
                let t = segment.t0 + (segment.t1 - segment.t0) * f64::from(k) / 4.0;
                let (p, _) = segment.at(t, tol)?;
                lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
                hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
            }
        }
    }
    let scale = (hi.x - lo.x).max(hi.y - lo.y);
    if !scale.is_finite() || scale <= 0.0 {
        return Ok(None);
    }

    // Each loop must close, piece to piece, to a millionth of the chart, or
    // to what the two edges meeting there own: fitted sections each miss
    // their junction by up to their stated tolerance, possibly opposite
    // ways, and the vertex there records any wider gap it absorbed. The
    // chart reads that slack through the surface's stretch there.
    let reach = scale * 1e-6 + tol.parametric();
    let owned = |edge: &Shape| {
        let own = model
            .node(edge)
            .and_then(|n| n.data().as_edge())
            .map_or(0.0, |d| d.tolerance.get());
        let ends = crate::edge_vertices(model, edge)
            .ok()
            .flatten()
            .map_or(0.0, |(a, b)| {
                [a, b]
                    .iter()
                    .filter_map(|v| model.node(v).and_then(|n| n.data().as_vertex()))
                    .map(|d| d.tolerance.get())
                    .fold(0.0, f64::max)
            });
        (own, ends)
    };
    for segments in &loops {
        for (k, segment) in segments.iter().enumerate() {
            let next = &segments[(k + 1) % segments.len()];
            let (end, _) = segment.at(segment.t1, tol)?;
            let (start, _) = next.at(next.t0, tol)?;
            let (a, b) = (owned(&segment.edge), owned(&next.edge));
            let slack = (a.0 + b.0).max(a.1).max(b.1);
            let stretch = placed.d1_at(end.x, end.y, tol).map_or(1.0, |(du, dv)| {
                du.magnitude().min(dv.magnitude()).max(tol.confusion())
            });
            if end.distance(start) > reach + slack / stretch {
                return Ok(None);
            }
        }
    }

    // Which loop is the boundary and which the holes: the boundary encloses
    // the most chart, whichever way each was wound. A spline is integrated
    // knot span by knot span: one rule across a closed spline of many spans
    // can miss its area entirely and give the wrong winding.
    let mut areas = Vec::with_capacity(loops.len());
    for segments in &loops {
        let mut area = 0.0;
        for segment in segments {
            let (lo, hi) = (segment.t0.min(segment.t1), segment.t0.max(segment.t1));
            let mut breaks = vec![segment.t0, segment.t1];
            if let PlanarCurve::BSpline(spline) = &segment.curve {
                breaks.extend(
                    spline
                        .knots()
                        .distinct()
                        .into_iter()
                        .map(|(k, _)| k)
                        .filter(|k| *k > lo && *k < hi),
                );
            }
            breaks.sort_by(f64::total_cmp);
            if segment.t1 < segment.t0 {
                breaks.reverse();
            }
            for pair in breaks.windows(2) {
                for (t, w) in gauss_legendre_rule(pair[0], pair[1]) {
                    let (p, d) = segment.at(t, tol)?;
                    area += p.x * d.y * w;
                }
            }
        }
        areas.push(area);
    }
    let Some(outer) = (0..areas.len()).max_by(|a, b| areas[*a].abs().total_cmp(&areas[*b].abs()))
    else {
        return Ok(None);
    };
    let loops = loops
        .into_iter()
        .zip(&areas)
        .enumerate()
        .map(|(k, (segments, area))| {
            let region = if k == outer { 1.0 } else { -1.0 };
            (segments, region * area.signum())
        })
        .collect();
    Ok(Some(Walked {
        placed,
        handedness,
        loops,
        lo,
        scale,
    }))
}

/// A wire's pcurve pieces in walking order, or `None` where one is missing.
///
/// A seam bounds its face twice, once down each side of the chart; which
/// side an occurrence takes is the one continuing the point already walked
/// to, so the walk starts off a seam where it can. A piece that starts a
/// whole period from where the last one ended is moved by that period.
fn walk(
    model: &Model,
    surface: ogeom_topo::SurfaceId,
    placed: &SurfaceGeometry,
    wire: &Shape,
    period: Vector2,
    strict: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Segment>>> {
    let mut edges = model.ordered_children_of(wire)?;
    let is_seam = |e: &Shape| {
        model
            .node(e)
            .and_then(|n| n.data().as_edge())
            .and_then(|d| d.pcurve_for(surface, e.location()))
            .is_some_and(|r| matches!(r, EdgeRepr::Seam { .. }))
    };
    if let Some(start) = edges.iter().position(|e| !is_seam(e)) {
        edges.rotate_left(start);
    }
    let mut segments: Vec<Segment> = Vec::with_capacity(edges.len());
    let mut last: Option<Point2> = None;
    for edge in &edges {
        let Some(repr) = model
            .node(edge)
            .and_then(|n| n.data().as_edge())
            .and_then(|d| d.pcurve_for(surface, edge.location()))
        else {
            return Ok(None);
        };
        let reversed = edge.orientation() == Orientation::Reversed;
        let (ids, range) = match repr {
            EdgeRepr::PCurve { curve, range, .. } => (vec![*curve], *range),
            EdgeRepr::Seam {
                forward,
                reversed: back,
                range,
                ..
            } => {
                let preferred = if reversed {
                    [*back, *forward]
                } else {
                    [*forward, *back]
                };
                (preferred.to_vec(), *range)
            }
            _ => return Ok(None),
        };
        let (t0, t1) = if reversed { (range.1, range.0) } else { range };
        // Of the candidate sides, the one whose start lies nearest the walk,
        // after folding by whole periods.
        let mut best: Option<(f64, Segment)> = None;
        for id in ids {
            let Some(curve) = model.geometry().pcurve(id) else {
                return Ok(None);
            };
            let start = Point2::ORIGIN + curve.derivatives_at(t0, 0, tol)?[0];
            let shift = match last {
                Some(at) => fold(at - start, period),
                None => Vector2::new(0.0, 0.0),
            };
            let miss = last.map_or(0.0, |at| (start + shift).distance(at));
            if best.as_ref().is_none_or(|(held, _)| miss < *held) {
                best = Some((
                    miss,
                    Segment {
                        edge: edge.clone(),
                        curve: curve.clone(),
                        t0,
                        t1,
                        shift,
                    },
                ));
            }
        }
        let Some((_, segment)) = best else {
            return Ok(None);
        };
        if strict && !lies_on_edge(model, edge, placed, &segment, tol)? {
            return Ok(None);
        }
        last = Some(segment.at(segment.t1, tol)?.0);
        segments.push(segment);
    }
    if segments.is_empty() {
        return Ok(None);
    }
    Ok(Some(segments))
}

/// Whether a piece's pcurve, lifted through the surface, runs along its
/// edge's own curve: to a hundred times the confusion distance, or to the
/// edge's own stated tolerance where that is looser. A fitted pcurve can
/// stray from its edge by the edge's tolerance, and the region it bounds is
/// measured that closely; one straying further is left to the mesh, which
/// takes its boundary from the edge. An edge with no curve of its own (a
/// pole) has nothing to stray from.
///
/// The two curves need not share a parameter, so each lifted point is
/// measured against the nearest point of the edge's curve over its range.
fn lies_on_edge(
    model: &Model,
    edge: &Shape,
    placed: &SurfaceGeometry,
    segment: &Segment,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .and_then(|d| d.curve3d())
    else {
        return Ok(true);
    };
    let Some(curve) = model.geometry().curve(*curve) else {
        return Ok(false);
    };
    let curve = curve
        .clone()
        .transformed(&edge.transform(model.datums())?, tol)?;
    // The edge's own tolerance is the slop it states: a pcurve within it
    // bounds the region that closely, far closer than a mesh would.
    let stated = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .map_or(0.0, |d| d.tolerance.get());
    let reach = (tol.confusion() * 100.0).max(stated);
    for k in 1..=5 {
        let t = segment.t0 + (segment.t1 - segment.t0) * f64::from(k) / 6.0;
        let (at, _) = segment.at(t, tol)?;
        let at = into_domain(placed, at);
        let lifted = placed.point_at(at.x, at.y, tol)?;
        if nearest(&curve, *range, lifted, tol)? > reach {
            return Ok(false);
        }
    }
    Ok(true)
}

/// A chart point brought into the surface's domain where a pcurve fitted
/// along its border strays past it by rounding, a millionth of the span or
/// less; any further is left for the evaluation to refuse.
pub(crate) fn into_domain(surface: &SurfaceGeometry, at: Point2) -> Point2 {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let onto = |x: f64, lo: f64, hi: f64, periodic: bool| {
        let slack = (hi - lo) * 1e-6;
        if periodic || !(lo..=hi).contains(&x) && (x < lo - slack || x > hi + slack) {
            x
        } else {
            x.clamp(lo, hi)
        }
    };
    Point2::new(
        onto(at.x, u0, u1, surface.is_periodic_u()),
        onto(at.y, v0, v1, surface.is_periodic_v()),
    )
}

/// The distance from `target` to `curve` over `range`: the best of a
/// sampling, narrowed by golden sections about it.
fn nearest(curve: &Curve, range: (f64, f64), target: Point, tol: Tolerances) -> OgeomResult<f64> {
    const SAMPLES: u32 = 32;
    let at = |t: f64| -> OgeomResult<f64> { Ok(curve.point_at(t, tol)?.distance(target)) };
    let step = (range.1 - range.0) / f64::from(SAMPLES);
    let mut best = (range.0, at(range.0)?);
    for k in 1..=SAMPLES {
        let t = range.0 + step * f64::from(k);
        let d = at(t)?;
        if d < best.1 {
            best = (t, d);
        }
    }
    let (mut a, mut b) = (
        (best.0 - step).max(range.0.min(range.1)),
        (best.0 + step).min(range.0.max(range.1)),
    );
    let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
    for _ in 0..60 {
        let (c, d) = (b - (b - a) * ratio, a + (b - a) * ratio);
        if at(c)? < at(d)? {
            b = d;
        } else {
            a = c;
        }
    }
    Ok(best.1.min(at(f64::midpoint(a, b))?))
}

/// The whole number of periods nearest to `gap`, along each periodic
/// direction.
fn fold(gap: Vector2, period: Vector2) -> Vector2 {
    let along = |g: f64, p: f64| if p > 0.0 { (g / p).round() * p } else { 0.0 };
    Vector2::new(along(gap.x, period.x), along(gap.y, period.y))
}

/// The most times the panels are doubled before the face is left to the
/// mesh.
const DOUBLINGS: u32 = 5;

/// Agreement asked of two runs, the second on panels twice as fine,
/// relative to the size of what they integrate.
const AGREE: f64 = 1e-10;

/// A quarter turn: the widest panel on an angular parameter, over which a
/// trigonometric polynomial integrates to rounding under the ten-point
/// rule.
const QUARTER: f64 = core::f64::consts::FRAC_PI_2;

/// The boundary rules a short panel may take, by their number of points.
const ORDERS: [usize; 5] = [3, 4, 5, 6, 7];

/// The Gauss-Legendre rules of [`ORDERS`] on `[-1, 1]`: nodes and weights.
const RULES: [(&[f64], &[f64]); 5] = [
    (
        &[-0.7745966692414834, 0.0, 0.7745966692414834],
        &[0.5555555555555557, 0.8888888888888888, 0.5555555555555557],
    ),
    (
        &[
            -0.8611363115940526,
            -0.33998104358485626,
            0.33998104358485626,
            0.8611363115940526,
        ],
        &[
            0.34785484513745357,
            0.6521451548625464,
            0.6521451548625464,
            0.34785484513745357,
        ],
    ),
    (
        &[
            -0.906179845938664,
            -0.5384693101056831,
            0.0,
            0.5384693101056831,
            0.906179845938664,
        ],
        &[
            0.23692688505618928,
            0.4786286704993663,
            0.5688888888888887,
            0.4786286704993663,
            0.23692688505618928,
        ],
    ),
    (
        &[
            -0.9324695142031519,
            -0.6612093864662645,
            -0.2386191860831969,
            0.2386191860831969,
            0.6612093864662645,
            0.9324695142031519,
        ],
        &[
            0.17132449237917027,
            0.3607615730481387,
            0.46791393457269104,
            0.46791393457269104,
            0.3607615730481387,
            0.17132449237917027,
        ],
    ),
    (
        &[
            -0.9491079123427586,
            -0.7415311855993945,
            -0.4058451513773972,
            0.0,
            0.4058451513773972,
            0.7415311855993945,
            0.9491079123427586,
        ],
        &[
            0.12948496616886973,
            0.27970539148927687,
            0.3818300505051187,
            0.4179591836734693,
            0.3818300505051187,
            0.27970539148927687,
            0.12948496616886973,
        ],
    ),
];

/// The highest frequency taken for a boundary integrand in an angle: the
/// second moments' on a sphere or a torus, a cubic in the point times the
/// `n dA`.
const FREQUENCY: f64 = 6.0;

/// What a boundary panel's rule may miss by, against the panel's own size:
/// below what ten points over a quarter turn are bounded by.
const PANEL_MISS: f64 = 1e-16;

/// The `order`-point Gauss-Legendre rule on `[a, b]`: one of [`ORDERS`],
/// or ten points.
fn rule(order: usize, a: f64, b: f64) -> Vec<(f64, f64)> {
    let (half, middle) = ((b - a) * 0.5, f64::midpoint(a, b));
    match ORDERS.iter().position(|&n| n == order) {
        Some(i) => {
            let (nodes, weights) = RULES[i];
            nodes
                .iter()
                .zip(weights)
                .map(|(x, w)| (middle + half * x, w * half))
                .collect()
        }
        None => gauss_legendre_rule(a, b).to_vec(),
    }
}

/// One sample of a face's integral: the surface point, its outward
/// `Su x Sv`, and the quadrature weight.
type Sample = (Point, Vector, f64);

impl ChartFace {
    /// A point of the face's surface, to take moments from.
    pub(crate) fn anchor(&self, tol: Tolerances) -> OgeomResult<Point> {
        let segment = &self.loops[0].0[0];
        let (at, _) = segment.at(segment.t0, tol)?;
        let at = into_domain(&self.surface, at);
        self.surface.point_at(at.x, at.y, tol)
    }

    /// The face's integral summed into an accumulator from `fresh`: every
    /// sample is handed to `contribute` as the surface point, its `n dA`
    /// with the weight's magnitude folded in, and the weight's sign (an area
    /// takes `|n dA|` times that sign, a volume the product). Each run sums
    /// into an accumulator of its own, and the one that settles is
    /// returned; `None` where none did.
    pub(crate) fn integrate<A>(
        &self,
        reference: Point,
        tol: Tolerances,
        fresh: impl Fn() -> A,
        contribute: impl Fn(&mut A, Point, Vector, f64),
    ) -> Option<A> {
        let mut held = self.run(1, reference, tol, &mut |_| {}).ok()?;
        for doubling in 1..=DOUBLINGS {
            // Each doubling costs twice the last: a cancelled watch is
            // honoured between them, and the caller's own checkpoint then
            // reports it.
            if ogeom_core::progress::checkpoint().is_err() {
                return None;
            }
            let mut sum = fresh();
            let proxy = self
                .run(1 << doubling, reference, tol, &mut |(p, n, w)| {
                    contribute(&mut sum, p, n * w.abs(), w.signum());
                })
                .ok()?;
            if settled(held, proxy) {
                return Some(sum);
            }
            held = proxy;
        }
        None
    }

    /// The face's samples with every panel split `fine` ways, each handed to
    /// `sink` in a fixed order, and the integrals of a few measures over
    /// them to compare runs by.
    fn run(
        &self,
        fine: u32,
        reference: Point,
        tol: Tolerances,
        sink: &mut dyn FnMut(Sample),
    ) -> OgeomResult<[f64; 5]> {
        let mut proxy = [0.0; 5];
        let size = self.scale.max(1.0);
        let mut take = |(p, n, w): Sample| {
            let e = p - reference;
            let row = [n.magnitude(), n.x, n.y, n.z, e.dot(n) / size];
            for (acc, x) in proxy.iter_mut().zip(row) {
                *acc += x * w;
            }
            sink((p, n, w));
        };
        for (segments, region) in &self.loops {
            for segment in segments {
                // A straight piece along which `v` does not move adds
                // nothing.
                if let PlanarCurve::Line(_) = segment.curve {
                    let (_, d) = segment.at(segment.t0, tol)?;
                    if d.y.abs() <= 1e-14 * d.x.abs() {
                        continue;
                    }
                }
                let breaks = self.outer_breaks(segment, tol)?;
                for pair in breaks.windows(2) {
                    let order = self.outer_order(segment, pair[0], pair[1], tol)?;
                    for k in 0..fine {
                        let a = pair[0] + (pair[1] - pair[0]) * f64::from(k) / f64::from(fine);
                        let b = pair[0] + (pair[1] - pair[0]) * f64::from(k + 1) / f64::from(fine);
                        for (t, wt) in rule(order, a, b) {
                            let (at, d) = segment.at(t, tol)?;
                            self.inner(at, region * wt * d.y, fine, tol, &mut take)?;
                        }
                    }
                }
            }
        }
        for &(start, step, region) in &self.bridges {
            for (t, wt) in gauss_legendre_rule(0.0, 1.0) {
                self.inner(start + step * t, region * wt * step.y, fine, tol, &mut take)?;
            }
        }
        Ok(proxy)
    }

    /// Where a boundary piece's panels break: its own ends, its knots, and
    /// every quarter turn it makes in an angular chart.
    fn outer_breaks(&self, segment: &Segment, tol: Tolerances) -> OgeomResult<Vec<f64>> {
        let (t0, t1) = (segment.t0, segment.t1);
        let mut breaks = vec![t0, t1];
        let (lo, hi) = (t0.min(t1), t0.max(t1));
        let knots = match &segment.curve {
            PlanarCurve::BSpline(spline) => spline
                .knots()
                .distinct()
                .into_iter()
                .map(|(k, _)| k)
                .collect(),
            _ => Vec::new(),
        };
        breaks.extend(knots.into_iter().filter(|k| *k > lo && *k < hi));
        // Where the piece crosses the surface's own knot lines, found on a
        // sampling and settled by bisection.
        let (u_knots, v_knots) = &self.knot_lines;
        if !u_knots.is_empty() || !v_knots.is_empty() {
            const N: usize = 32;
            let mut prev: Option<(f64, Point2)> = None;
            for k in 0..=N {
                #[allow(clippy::cast_precision_loss)]
                let t = t0 + (t1 - t0) * k as f64 / N as f64;
                let (p, _) = segment.at(t, tol)?;
                if let Some((tp, pp)) = prev {
                    for (lines, read) in [(&u_knots, 0usize), (&v_knots, 1usize)] {
                        let coord = |q: Point2| if read == 0 { q.x } else { q.y };
                        for line in lines.iter() {
                            let (a, b) = (coord(pp) - line, coord(p) - line);
                            if a * b >= 0.0 {
                                continue;
                            }
                            let (mut lo_t, mut hi_t, mut f_lo) = (tp, t, a);
                            for _ in 0..60 {
                                let mid = 0.5 * (lo_t + hi_t);
                                let f = coord(segment.at(mid, tol)?.0) - line;
                                if f.signum() == f_lo.signum() {
                                    (lo_t, f_lo) = (mid, f);
                                } else {
                                    hi_t = mid;
                                }
                            }
                            breaks.push(0.5 * (lo_t + hi_t));
                        }
                    }
                }
                prev = Some((t, p));
            }
        }
        // How far the piece reaches across the chart, in quarter turns of
        // whichever parameters are angles.
        let mut pieces = match segment.curve {
            PlanarCurve::Line(_) => 1.0,
            _ => ((t1 - t0).abs() / QUARTER).ceil(),
        };
        if !matches!(self.surface, SurfaceGeometry::Plane(_)) {
            let (mut u, mut v) = (
                (f64::INFINITY, f64::NEG_INFINITY),
                (f64::INFINITY, f64::NEG_INFINITY),
            );
            for k in 0..=8 {
                let t = t0 + (t1 - t0) * f64::from(k) / 8.0;
                let (p, _) = segment.at(t, tol)?;
                u = (u.0.min(p.x), u.1.max(p.x));
                v = (v.0.min(p.y), v.1.max(p.y));
            }
            pieces = pieces
                .max(((u.1 - u.0) / QUARTER).ceil())
                .max(((v.1 - v.0) / QUARTER).ceil());
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let pieces = pieces.clamp(1.0, 64.0) as u32;
        for k in 1..pieces {
            breaks.push(t0 + (t1 - t0) * f64::from(k) / f64::from(pieces));
        }
        breaks.sort_by(f64::total_cmp);
        if t1 < t0 {
            breaks.reverse();
        }
        breaks.dedup();
        Ok(breaks)
    }

    /// How many Gauss points a boundary panel from `a` to `b` takes: the
    /// fewest whose error bound is below [`PANEL_MISS`] of the panel's own
    /// size, or ten.
    ///
    /// The integrand along a panel is a trigonometric polynomial of the
    /// chart point, of frequency at most [`FREQUENCY`] in an angle (a
    /// polynomial in a length, whose chart size stands in for a radian),
    /// times the pcurve's `v` speed. On a polynomial piece of degree three
    /// or less, written `c(s) = c0 + A1 s + A2 s^2 + A3 s^3` over the panel
    /// as `s` runs over `[-1, 1]`, the integrand is analytic in the
    /// Bernstein ellipse of every `rho`, bounded there through the bounds
    /// of `|Im c|` and `|c'|`, and the `n`-point rule misses by at most
    /// `64/15 M rho^(2 - 2n) / (rho^2 - 1)`. Any other piece, or a surface
    /// whose integrand is of another kind, takes ten.
    fn outer_order(
        &self,
        segment: &Segment,
        a: f64,
        b: f64,
        tol: Tolerances,
    ) -> OgeomResult<usize> {
        let cubic = match &segment.curve {
            PlanarCurve::Line(_) => true,
            PlanarCurve::BSpline(spline) => !spline.is_rational() && spline.knots().degree() <= 3,
            _ => false,
        };
        let length = 1.0 / self.scale.max(f64::MIN_POSITIVE);
        let (wu, wv) = match self.surface {
            SurfaceGeometry::Sphere(_) | SurfaceGeometry::Torus(_) => (1.0, 1.0),
            SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_) => (1.0, length),
            SurfaceGeometry::Plane(_) => (length, length),
            _ => return Ok(10),
        };
        if !cubic {
            return Ok(10);
        }
        // The cubic's coefficients from its points at s = -1, -1/3, 1/3, 1.
        let (middle, half) = (f64::midpoint(a, b), (b - a) * 0.5);
        let mut y = [Point2::ORIGIN; 4];
        for (slot, s) in y.iter_mut().zip([-1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0]) {
            *slot = segment.at(middle + half * s, tol)?.0;
        }
        let odd = (y[3] - y[0]) * 0.5;
        let odd_inner = (y[2] - y[1]) * 0.5;
        let a3 = (odd - odd_inner * 3.0) * (9.0 / 8.0);
        let a1 = odd - a3;
        let a2 = ((y[0].to_vector() + y[3].to_vector()) - (y[1].to_vector() + y[2].to_vector()))
            * (9.0 / 16.0);
        let size = a1.y.abs() + a2.y.abs() + a3.y.abs();
        if size == 0.0 {
            return Ok(ORDERS[0]);
        }
        let reach = |c: Vector2| c.x.abs() * wu + c.y.abs() * wv;
        let (r1, r2, r3) = (reach(a1), reach(a2), reach(a3));
        let mut best = [f64::INFINITY; ORDERS.len()];
        for k in 0..48 {
            let rho = 1.05 * 1.25_f64.powi(k);
            let (big, small) = (0.5 * (rho + 1.0 / rho), 0.5 * (rho - 1.0 / rho));
            // |Im c| on the ellipse, and |dv/ds| there against the panel's
            // own size.
            let im = r1 * small
                + r2 * 2.0 * big * small
                + r3 * (3.0 * big * big + small * small) * small;
            let exponent = FREQUENCY * im;
            if exponent > 600.0 {
                break;
            }
            let speed = (a1.y.abs() + 2.0 * big * a2.y.abs() + 3.0 * big * big * a3.y.abs()) / size;
            let common = 64.0 / 15.0 * exponent.exp() * speed * rho * rho / (rho * rho - 1.0);
            for (slot, &n) in best.iter_mut().zip(&ORDERS) {
                #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
                let miss = common * rho.powi(-2 * n as i32);
                *slot = slot.min(miss);
            }
        }
        Ok(ORDERS
            .iter()
            .zip(best)
            .find(|&(_, miss)| miss <= PANEL_MISS)
            .map_or(10, |(&n, _)| n))
    }

    /// The inner integral from `u_ref` to the boundary point `at`, along
    /// `u` at its `v`, its samples weighted by `outer`.
    fn inner(
        &self,
        at: Point2,
        outer: f64,
        fine: u32,
        tol: Tolerances,
        sink: &mut dyn FnMut(Sample),
    ) -> OgeomResult<()> {
        let at = into_domain(&self.surface, at);
        let (ua, ub) = (self.u_ref, at.x);
        if ua == ub || outer == 0.0 {
            return Ok(());
        }
        // Along `u` at a fixed `v`, an analytic surface's point, `n dA` and
        // its length are trigonometric polynomials in `u` (polynomials on a
        // plane), which a quarter-turn panel takes to rounding: only the
        // boundary's own panels are refined between runs there.
        let refined = match self.surface {
            SurfaceGeometry::Plane(_)
            | SurfaceGeometry::Cylinder(_)
            | SurfaceGeometry::Cone(_)
            | SurfaceGeometry::Sphere(_)
            | SurfaceGeometry::Torus(_) => 1,
            _ => fine,
        };
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let pieces = if matches!(self.surface, SurfaceGeometry::Plane(_)) {
            1
        } else {
            ((ub - ua).abs() / QUARTER).ceil().clamp(1.0, 64.0) as u32
        } * refined;
        let u_knots = &self.knot_lines.0;
        let (lo, hi) = (ua.min(ub), ua.max(ub));
        let mut cuts: Vec<f64> = (0..=pieces)
            .map(|k| ua + (ub - ua) * f64::from(k) / f64::from(pieces))
            .collect();
        cuts.extend(u_knots.iter().copied().filter(|k| *k > lo && *k < hi));
        cuts.sort_by(f64::total_cmp);
        if ub < ua {
            cuts.reverse();
        }
        cuts.dedup();
        let isoline = Isoline::of(&self.surface, at.y, tol);
        for pair in cuts.windows(2) {
            let (a, b) = (pair[0], pair[1]);
            for (u, wu) in gauss_legendre_rule(a, b) {
                let (p, du, dv) = match isoline.as_ref().and_then(|line| line.at(u, tol)) {
                    Some(found) => found,
                    None => self.surface.point_d1_at(u, at.y, tol)?,
                };
                sink((p, du.cross(dv) * self.sign, outer * wu));
            }
        }
        Ok(())
    }
}

/// A polynomial spline patch read along one `v`: its control net summed
/// across `v` once, into the control points of the row at `v` and of its
/// `v` derivative. Every sample of an inner integral shares that `v`, so a
/// sample costs one row's worth of basis functions instead of the grid's.
struct Isoline<'a> {
    knots: &'a ogeom_math::KnotVector,
    row: Vec<Point>,
    across: Vec<Vector>,
}

impl<'a> Isoline<'a> {
    fn of(surface: &'a SurfaceGeometry, v: f64, tol: Tolerances) -> Option<Self> {
        let SurfaceGeometry::BSpline(patch) = surface else {
            return None;
        };
        if patch.is_rational() {
            return None;
        }
        let v_knots = patch.v_knots();
        let span = v_knots.span(v, tol).ok()?;
        let q = v_knots.degree();
        let basis = v_knots.basis_derivatives(span, v, 1);
        let grid = patch.grid();
        let mut row = Vec::with_capacity(grid.u_count());
        let mut across = Vec::with_capacity(grid.u_count());
        for i in 0..grid.u_count() {
            let (mut p, mut d) = (Vector::ZERO, Vector::ZERO);
            for k in 0..=q {
                let c = grid.get(i, span - q + k)?.scaled.to_vector();
                p += c * basis[0][k];
                d += c * basis[1][k];
            }
            row.push(Point::ORIGIN + p);
            across.push(d);
        }
        Some(Self {
            knots: patch.u_knots(),
            row,
            across,
        })
    }

    /// The point, `du` and `dv` at `u`; `None` off the knots' domain.
    fn at(&self, u: f64, tol: Tolerances) -> Option<(Point, Vector, Vector)> {
        let span = self.knots.span(u, tol).ok()?;
        let p = self.knots.degree();
        let basis = self.knots.basis_derivatives(span, u, 1);
        let (mut point, mut du, mut dv) = (Vector::ZERO, Vector::ZERO, Vector::ZERO);
        for k in 0..=p {
            let c = self.row[span - p + k].to_vector();
            point += c * basis[0][k];
            du += c * basis[1][k];
            dv += self.across[span - p + k] * basis[0][k];
        }
        Some((Point::ORIGIN + point, du, dv))
    }
}

/// Whether two runs agree, against the size of what they integrate.
fn settled(a: [f64; 5], b: [f64; 5]) -> bool {
    let size: f64 = b.iter().map(|x| x.abs()).sum();
    let miss = a
        .iter()
        .zip(&b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max);
    miss <= AGREE * size || size == 0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Each short rule integrates every polynomial of degree below twice
    /// its points exactly, on a panel away from the origin.
    #[test]
    fn the_short_rules_are_exact_to_their_degree() {
        let (a, b) = (0.75, 2.0);
        for n in ORDERS {
            let points = rule(n, a, b);
            assert_eq!(points.len(), n);
            for degree in 0..2 * n {
                #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
                let k = degree as i32;
                let exact = (b.powi(k + 1) - a.powi(k + 1)) / f64::from(k + 1);
                let sum: f64 = points.iter().map(|&(x, w)| w * x.powi(k)).sum();
                assert!(
                    (sum - exact).abs() <= 1e-14 * exact.abs(),
                    "{n} points, degree {degree}: {sum} against {exact}"
                );
            }
        }
    }
}
