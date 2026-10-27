//! Where a face's boundary crosses itself in its surface's chart.
//!
//! A face is the region of its surface its wires enclose, read in the
//! surface's chart, the material on the left of every edge as the face
//! walks it. Where two of its edges' pcurves cross, the boundary walks a
//! figure eight: the region on one side of the crossing is wound the other
//! way round from the region on the other, and the face encloses no single
//! region. The mesher, the classifier and a boolean's arrangement all read
//! the face in its chart and each make something different of it.
//!
//! Two edges meet at the vertex they share, and fitted pcurves wander about
//! it by the tolerance their edges state: a crossing that lands within that
//! reach of a vertex of each edge is where the two meet, not a crossing.
//! Nor is a crossing of two consecutive edges whose loop, from it to the
//! vertex they meet at, the edges' own tolerances account for. Each
//! pcurve's image stands within its edge's tolerance of the edge's curve,
//! so where the images cross the curves stand within the two tolerances
//! of each other, and toward the vertex, where the curves meet, they draw
//! together: the images of two curves that do not cross close a loop no
//! wider than twice the two tolerances. A sharp corner's fitted pcurves
//! close such a sliver as readily as they miss each other. A touch, where
//! the two pcurves meet without crossing, is not reported either; nor is
//! anything about a degenerate edge, which stands for one point however
//! it runs in the chart.
//!
//! Each pcurve is drawn as a polyline, the walk carried across a closed
//! surface's seam so that consecutive edges meet in one copy of the chart,
//! and the segments of two edges that cross seed Newton on the pcurves
//! themselves. On a closed surface each pair is also asked a period apart
//! either way. Two pcurves that cross within the sag of their polylines
//! may go unseen; what is reported is a crossing of the pcurves
//! themselves, found to rounding.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve2d as _, PlanarCurve, Surface as _};
use ogeom_math::{Point, Point2, Vector2, solve};
use ogeom_topo::{EdgeRepr, Model, Orientation, Shape, TShapeId};

/// How many segments a curved pcurve is drawn with. A straight one is
/// drawn as the one segment it is.
const SEGMENTS: u32 = 32;

/// A point where a face's boundary crosses itself in its surface's chart.
#[derive(Debug, Clone, PartialEq)]
pub struct ChartCrossing {
    /// The two edges crossing, each as (wire, edge): the index of the wire
    /// among the face's wires as stored, and of the edge in that wire's
    /// order. The first comes first in that order.
    pub edges: [(usize, usize); 2],
    /// Where on each edge's pcurve, in that pcurve's own parameter.
    pub parameters: [f64; 2],
    /// Whole periods of the chart the second edge's pcurve is moved by,
    /// against the first's, to meet it there, beyond those the walk
    /// carries it by across the seam. Zero on a surface that does not
    /// close.
    pub apart: [i64; 2],
    /// Where, on the surface, placed as the face is.
    pub at: Point,
}

/// One edge as the face walks it.
struct Walked<'m> {
    wire: usize,
    edge: usize,
    node: TShapeId,
    pcurve: &'m PlanarCurve,
    range: (f64, f64),
    /// Added to the pcurve's points to carry it to where the walk has got
    /// to across the seam.
    shift: Vector2,
    /// The edge's vertices, placed, each with its tolerance.
    vertices: Vec<(Point, f64)>,
    tolerance: f64,
    /// Drawn and asked: not a degenerate edge.
    drawn: bool,
    /// Walked from its range's end to its start.
    reversed: bool,
    /// The vertices it is walked from and to.
    ends: Option<(TShapeId, TShapeId)>,
    /// The polyline: parameters and points, shifted.
    samples: Vec<(f64, Point2)>,
}

/// Every place where the boundary of `face` crosses itself in its
/// surface's chart, as the module documentation describes, in the order of
/// the edges crossing. Nothing where an edge of the face has no pcurve on
/// its surface: what cannot be drawn is not read.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the face
/// or its geometry is not in the model.
pub fn chart_crossings(
    model: &Model,
    face: &Shape,
    tol: Tolerances,
) -> OgeomResult<Vec<ChartCrossing>> {
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "face names a surface not in this model");
    };
    let ((u0, u1), (v0, v1)) = surface.domain();
    let periods = [
        (surface.is_periodic_u() || surface.is_closed_u(tol))
            .then_some(u1 - u0)
            .filter(|p| *p > 0.0 && p.is_finite()),
        (surface.is_periodic_v() || surface.is_closed_v(tol))
            .then_some(v1 - v0)
            .filter(|p| *p > 0.0 && p.is_finite()),
    ];
    let placement = face.transform(model.datums())?;
    let stored = face.oriented(Orientation::Forward);

    let mut walked: Vec<Walked<'_>> = Vec::new();
    for (w, wire) in model.ordered_children_of(&stored)?.into_iter().enumerate() {
        let mut last: Option<Point2> = None;
        for (i, edge) in model.ordered_children_of(&wire)?.into_iter().enumerate() {
            let Some(edge_data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                return Ok(Vec::new());
            };
            let reversed = edge.orientation() == Orientation::Reversed;
            let (id, range) = match edge_data.pcurve_for(data.surface, edge.location()) {
                Some(EdgeRepr::PCurve { curve, range, .. }) => (*curve, *range),
                Some(EdgeRepr::Seam {
                    forward,
                    reversed: back,
                    range,
                    ..
                }) => (if reversed { *back } else { *forward }, *range),
                _ => return Ok(Vec::new()),
            };
            let Some(pcurve) = model.geometry().pcurve(id) else {
                ogeom_bail!(Dangling, "an edge names geometry not in this model");
            };
            let (first, end) = {
                let (a, b) = (
                    pcurve.point_at(range.0, tol)?,
                    pcurve.point_at(range.1, tol)?,
                );
                if reversed { (b, a) } else { (a, b) }
            };
            // Whole periods that bring this edge's start to where the last
            // one ended.
            let mut shift = Vector2::new(0.0, 0.0);
            if let Some(last) = last {
                for (axis, period) in periods.iter().enumerate() {
                    let Some(period) = period else { continue };
                    let (to, from) = if axis == 0 {
                        (last.x, first.x)
                    } else {
                        (last.y, first.y)
                    };
                    let whole = ((to - from) / period).round() * period;
                    if axis == 0 {
                        shift.x = whole;
                    } else {
                        shift.y = whole;
                    }
                }
            }
            last = Some(end + shift);
            let mut vertices = Vec::with_capacity(2);
            for v in model.children_of(&edge)? {
                let Some(vd) = model.node(&v).and_then(|n| n.data().as_vertex()) else {
                    continue;
                };
                vertices.push((
                    v.transform(model.datums())?.apply(vd.point),
                    vd.tolerance.get(),
                ));
            }
            let pieces = if matches!(pcurve, PlanarCurve::Line(_)) {
                1
            } else {
                SEGMENTS
            };
            let samples = if edge_data.degenerate {
                Vec::new()
            } else {
                (0..=pieces)
                    .map(|k| {
                        let t = range.0 + (range.1 - range.0) * f64::from(k) / f64::from(pieces);
                        pcurve.point_at(t, tol).map(|p| (t, p + shift))
                    })
                    .collect::<OgeomResult<_>>()?
            };
            let ends =
                crate::build::edge_vertices(model, &edge)?.map(|(a, b)| (a.node(), b.node()));
            walked.push(Walked {
                wire: w,
                edge: i,
                node: edge.node(),
                reversed,
                ends,
                pcurve,
                range,
                shift,
                vertices,
                tolerance: edge_data.tolerance.get(),
                drawn: !edge_data.degenerate,
                samples,
            });
        }
    }

    // The period moves each pair is asked at.
    let mut moves: Vec<[i64; 2]> = vec![[0, 0]];
    for (axis, period) in periods.iter().enumerate() {
        if period.is_some() {
            let mut more = Vec::new();
            for m in &moves {
                for k in [-1, 1] {
                    let mut moved = *m;
                    moved[axis] = k;
                    more.push(moved);
                }
            }
            moves.extend(more);
        }
    }
    let offset = |m: [i64; 2]| {
        #[allow(clippy::cast_precision_loss, reason = "one period either way")]
        Vector2::new(
            periods[0].map_or(0.0, |p| p * m[0] as f64),
            periods[1].map_or(0.0, |p| p * m[1] as f64),
        )
    };
    let boxes: Vec<(Point2, Point2)> = walked.iter().map(|e| bounds(&e.samples)).collect();

    let mut found: Vec<ChartCrossing> = Vec::new();
    for a in 0..walked.len() {
        if !walked[a].drawn {
            continue;
        }
        for b in a + 1..walked.len() {
            if !walked[b].drawn || walked[a].node == walked[b].node {
                continue;
            }
            for &m in &moves {
                let by = offset(m);
                let (lo, hi) = (boxes[b].0 + by, boxes[b].1 + by);
                if lo.x > boxes[a].1.x
                    || hi.x < boxes[a].0.x
                    || lo.y > boxes[a].1.y
                    || hi.y < boxes[a].0.y
                {
                    continue;
                }
                found.extend(pair_crossings(
                    &walked[a], &walked[b], by, m, surface, &placement, tol,
                )?);
            }
        }
    }
    Ok(found)
}

/// The box of a polyline's points.
fn bounds(samples: &[(f64, Point2)]) -> (Point2, Point2) {
    let mut lo = Point2::new(f64::INFINITY, f64::INFINITY);
    let mut hi = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    for (_, p) in samples {
        lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
        hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
    }
    (lo, hi)
}

/// The crossings of two walked edges, the second moved by `by`, that are
/// no meeting at a vertex and no touch, each once, in order along the
/// first.
fn pair_crossings(
    a: &Walked<'_>,
    b: &Walked<'_>,
    by: Vector2,
    apart: [i64; 2],
    surface: &ogeom_geom::SurfaceGeometry,
    placement: &ogeom_math::Transform,
    tol: Tolerances,
) -> OgeomResult<Vec<ChartCrossing>> {
    let mut found: Vec<(Point2, ChartCrossing)> = Vec::new();
    for i in 1..a.samples.len() {
        let (pa, qa) = (a.samples[i - 1], a.samples[i]);
        for j in 1..b.samples.len() {
            let (pb, qb) = (b.samples[j - 1], b.samples[j]);
            let Some((sa, sb)) = segments_cross(pa.1, qa.1, pb.1 + by, qb.1 + by) else {
                continue;
            };
            let seed = [pa.0 + (qa.0 - pa.0) * sa, pb.0 + (qb.0 - pb.0) * sb];
            let Some([ta, tb]) = polish(a, b, by, seed, tol) else {
                continue;
            };
            let uv = a.pcurve.point_at(ta, tol)?;
            // Two segments either side of a sample seed the same crossing.
            if found
                .iter()
                .any(|(at, _)| at.distance(uv) <= tol.confusion())
            {
                continue;
            }
            let (Ok((point, du, dv)), Ok(da), Ok(db)) = (
                surface.point_d1_at(uv.x, uv.y, tol),
                a.pcurve.d1_at(ta, tol),
                b.pcurve.d1_at(tb, tol),
            ) else {
                continue;
            };
            let at = placement.apply(point);
            // A touch is no crossing: the two must meet at an angle on the
            // surface.
            let (ta3, tb3) = (du * da.x + dv * da.y, du * db.x + dv * db.y);
            let (la, lb) = (ta3.magnitude(), tb3.magnitude());
            if la <= 0.0
                || lb <= 0.0
                || ta3.cross(tb3).magnitude() <= crate::self_cross::CROSSING_SINE * la * lb
            {
                continue;
            }
            // Where the two meet at a vertex of each, they meet.
            let widest = a.tolerance.max(b.tolerance).max(tol.approximation());
            let near = |e: &Walked<'_>| {
                e.vertices
                    .iter()
                    .any(|(v, t)| at.distance(*v) <= crate::check::CHART_GAP_SLACK * widest.max(*t))
            };
            if near(a) && near(b) {
                continue;
            }
            if apart == [0, 0] && sliver(a, ta, b, tb, surface, tol)? {
                continue;
            }
            found.push((
                uv,
                ChartCrossing {
                    edges: [(a.wire, a.edge), (b.wire, b.edge)],
                    parameters: [ta, tb],
                    apart,
                    at,
                },
            ));
        }
    }
    Ok(found.into_iter().map(|(_, crossing)| crossing).collect())
}

/// Whether `a` and `b`, crossing at the given parameters of their
/// pcurves, are consecutive in their wire and run within twice the sum of
/// their tolerances of each other, on the surface, from the crossing to
/// the vertex they meet at (see the module documentation).
fn sliver(
    a: &Walked<'_>,
    ta: f64,
    b: &Walked<'_>,
    tb: f64,
    surface: &ogeom_geom::SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<bool> {
    const SAMPLES: u32 = 32;
    let (Some(ea), Some(eb)) = (a.ends, b.ends) else {
        return Ok(false);
    };
    // The edge walked into the vertex they meet at, and the one walked
    // out of it.
    let ((into, t_into), (out, t_out)) = if ea.1 == eb.0 {
        ((a, ta), (b, tb))
    } else if eb.1 == ea.0 {
        ((b, tb), (a, ta))
    } else {
        return Ok(false);
    };
    if into.wire != out.wire {
        return Ok(false);
    }
    let walked_end = |e: &Walked<'_>, last: bool| {
        if last != e.reversed {
            e.range.1
        } else {
            e.range.0
        }
    };
    let stretch = |e: &Walked<'_>, from: f64, to: f64| -> OgeomResult<Vec<Point>> {
        (0..=SAMPLES)
            .map(|k| {
                let uv = e
                    .pcurve
                    .point_at(from + (to - from) * f64::from(k) / f64::from(SAMPLES), tol)?;
                surface.point_at(uv.x, uv.y, tol)
            })
            .collect()
    };
    let first = stretch(into, t_into, walked_end(into, true))?;
    let second = stretch(out, walked_end(out, false), t_out)?;
    let width = 2.0 * (into.tolerance + out.tolerance);
    let within = |points: &[Point], line: &[Point]| {
        points.iter().all(|p| {
            line.windows(2)
                .any(|s| distance_to_segment(*p, s[0], s[1]) <= width)
        })
    };
    Ok(within(&first, &second) && within(&second, &first))
}

/// How far `p` stands from the segment from `a` to `b`.
fn distance_to_segment(p: Point, a: Point, b: Point) -> f64 {
    let along = b - a;
    let length = along.dot(along);
    if length <= 0.0 {
        return p.distance(a);
    }
    let t = ((p - a).dot(along) / length).clamp(0.0, 1.0);
    p.distance(a + along * t)
}

/// Where two chart segments cross, as the fraction along each; ends
/// included, parallel segments never.
fn segments_cross(p: Point2, q: Point2, r: Point2, s: Point2) -> Option<(f64, f64)> {
    let (d, e) = (q - p, s - r);
    let denominator = d.cross(e);
    if denominator == 0.0 {
        return None;
    }
    let w = r - p;
    let (t, u) = (w.cross(e) / denominator, w.cross(d) / denominator);
    ((0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u)).then_some((t, u))
}

/// Newton from `seed` onto where the two pcurves, the second moved by
/// `by`, meet, inside both edges' ranges; `None` where it does not land
/// there.
fn polish(
    a: &Walked<'_>,
    b: &Walked<'_>,
    by: Vector2,
    seed: [f64; 2],
    tol: Tolerances,
) -> Option<[f64; 2]> {
    let clamp = |t: f64, (lo, hi): (f64, f64)| t.max(lo.min(hi)).min(lo.max(hi));
    let system = |x: &[f64; 2]| {
        let (ta, tb) = (clamp(x[0], a.range), clamp(x[1], b.range));
        let pa = a.pcurve.point_at(ta, tol).unwrap_or(Point2::ORIGIN) + a.shift;
        let pb = b.pcurve.point_at(tb, tol).unwrap_or(Point2::ORIGIN) + b.shift + by;
        let da = a.pcurve.d1_at(ta, tol).unwrap_or(Vector2::new(0.0, 0.0));
        let db = b.pcurve.d1_at(tb, tol).unwrap_or(Vector2::new(0.0, 0.0));
        ([pa.x - pb.x, pa.y - pb.y], [[da.x, -db.x], [da.y, -db.y]])
    };
    let criteria = solve::Criteria {
        residual: tol.confusion() * 0.01,
        step: tol.parametric(),
        max_iterations: 40,
    };
    let (x, _, _, _) = solve::newton_system_fixed(system, seed, criteria).ok()?;
    let (ta, tb) = (clamp(x[0], a.range), clamp(x[1], b.range));
    let pa = a.pcurve.point_at(ta, tol).ok()? + a.shift;
    let pb = b.pcurve.point_at(tb, tol).ok()? + b.shift + by;
    (pa.distance(pb) <= tol.confusion()).then_some([ta, tb])
}
