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
//! fine until two runs agree to a part in ten billion.

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
        SurfaceGeometry::Extrusion(e) => conic(e.curve()),
        SurfaceGeometry::Revolution(r) => conic(r.curve()),
        _ => false,
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
    // The pcurves are the unplaced surface's; a rigid placement keeps its
    // chart, a scaling one would change the metric under them.
    let placement = face.transform(model.datums())?;
    if !matches!(
        placement.kind(),
        ogeom_math::TransformKind::Identity
            | ogeom_math::TransformKind::Translation
            | ogeom_math::TransformKind::Rotation
    ) {
        return Ok(None);
    }
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
        let Some(segments) = walk(model, data.surface, &placed, &wire, period, tol)? else {
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

    // Each loop must close, piece to piece, to a millionth of the chart.
    let reach = scale * 1e-6 + tol.parametric();
    for segments in &loops {
        for (k, segment) in segments.iter().enumerate() {
            let next = &segments[(k + 1) % segments.len()];
            let (end, _) = segment.at(segment.t1, tol)?;
            let (start, _) = next.at(next.t0, tol)?;
            if end.distance(start) > reach {
                return Ok(None);
            }
        }
    }

    // Which loop is the boundary and which the holes: the boundary encloses
    // the most chart, whichever way each was wound.
    let mut areas = Vec::with_capacity(loops.len());
    for segments in &loops {
        let mut area = 0.0;
        for segment in segments {
            for (t, w) in gauss_legendre_rule(segment.t0, segment.t1) {
                let (p, d) = segment.at(t, tol)?;
                area += p.x * d.y * w;
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
    let u_ref = if surface.is_periodic_u() {
        lo.x
    } else {
        lo.x.clamp(u0, u1)
    };
    Ok(Some(ChartFace {
        surface: placed,
        loops,
        sign: if face.orientation() == Orientation::Reversed {
            -1.0
        } else {
            1.0
        },
        u_ref,
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
        if !lies_on_edge(model, edge, placed, &segment, tol)? {
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
/// edge's own curve to a hundred times the confusion distance. A fitted
/// pcurve can stray from its edge by the edge's tolerance, and the region
/// it bounds would be measured that far off; the mesh, which takes its
/// boundary from the edge, is asked instead. An edge with no curve of its
/// own (a pole) has nothing to stray from.
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
    let reach = tol.confusion() * 100.0;
    for k in 1..=5 {
        let t = segment.t0 + (segment.t1 - segment.t0) * f64::from(k) / 6.0;
        let (at, _) = segment.at(t, tol)?;
        let lifted = placed.point_at(at.x, at.y, tol)?;
        if nearest(&curve, *range, lifted, tol)? > reach {
            return Ok(false);
        }
    }
    Ok(true)
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

/// One sample of a face's integral: the surface point, its outward
/// `Su x Sv`, and the quadrature weight.
type Sample = (Point, Vector, f64);

impl ChartFace {
    /// A point of the face's surface, to take moments from.
    pub(crate) fn anchor(&self, tol: Tolerances) -> OgeomResult<Point> {
        let segment = &self.loops[0].0[0];
        let (at, _) = segment.at(segment.t0, tol)?;
        self.surface.point_at(at.x, at.y, tol)
    }

    /// Feed every sample of the face's integral to `contribute`, as the
    /// surface point, its `n dA` with the weight's magnitude folded in, and
    /// the weight's sign: an area takes `|n dA|` times that sign, a volume
    /// the product. `false` where the rule did not settle, and nothing was
    /// fed.
    pub(crate) fn integrate(
        &self,
        reference: Point,
        tol: Tolerances,
        contribute: &mut dyn FnMut(Point, Vector, f64),
    ) -> bool {
        let Ok((mut held, _)) = self.run(1, reference, tol) else {
            return false;
        };
        for doubling in 1..=DOUBLINGS {
            let Ok((proxy, samples)) = self.run(1 << doubling, reference, tol) else {
                return false;
            };
            if settled(held, proxy) {
                for (p, n, w) in samples {
                    contribute(p, n * w.abs(), w.signum());
                }
                return true;
            }
            held = proxy;
        }
        false
    }

    /// The face's samples with every panel split `fine` ways, and the
    /// integrals of a few measures over them to compare runs by.
    fn run(
        &self,
        fine: u32,
        reference: Point,
        tol: Tolerances,
    ) -> OgeomResult<([f64; 5], Vec<Sample>)> {
        let mut samples = Vec::new();
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
                    for k in 0..fine {
                        let a = pair[0] + (pair[1] - pair[0]) * f64::from(k) / f64::from(fine);
                        let b = pair[0] + (pair[1] - pair[0]) * f64::from(k + 1) / f64::from(fine);
                        for (t, wt) in gauss_legendre_rule(a, b) {
                            let (at, d) = segment.at(t, tol)?;
                            self.inner(at, region * wt * d.y, fine, tol, &mut samples)?;
                        }
                    }
                }
            }
        }
        let mut proxy = [0.0; 5];
        let size = self.scale.max(1.0);
        for (p, n, w) in &samples {
            let e = *p - reference;
            let row = [n.magnitude(), n.x, n.y, n.z, e.dot(*n) / size];
            for (acc, x) in proxy.iter_mut().zip(row) {
                *acc += x * w;
            }
        }
        Ok((proxy, samples))
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

    /// The inner integral from `u_ref` to the boundary point `at`, along
    /// `u` at its `v`, its samples weighted by `outer`.
    fn inner(
        &self,
        at: Point2,
        outer: f64,
        fine: u32,
        tol: Tolerances,
        samples: &mut Vec<Sample>,
    ) -> OgeomResult<()> {
        let (ua, ub) = (self.u_ref, at.x);
        if ua == ub || outer == 0.0 {
            return Ok(());
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let pieces = if matches!(self.surface, SurfaceGeometry::Plane(_)) {
            1
        } else {
            ((ub - ua).abs() / QUARTER).ceil().clamp(1.0, 64.0) as u32
        } * fine;
        for k in 0..pieces {
            let a = ua + (ub - ua) * f64::from(k) / f64::from(pieces);
            let b = ua + (ub - ua) * f64::from(k + 1) / f64::from(pieces);
            for (u, wu) in gauss_legendre_rule(a, b) {
                let p = self.surface.point_at(u, at.y, tol)?;
                let (du, dv) = self.surface.d1_at(u, at.y, tol)?;
                samples.push((p, du.cross(dv) * self.sign, outer * wu));
            }
        }
        Ok(())
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
