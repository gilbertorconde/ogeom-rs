//! The medial axis of a planar region: the locus of centres of maximal
//! inscribed circles — what tool-path generation and midline extraction are
//! built on.
//!
//! For a **convex** polygon the medial axis coincides with the straight
//! skeleton, and the shrinking-polygon construction computes it exactly:
//! every edge moves inward at unit speed, every vertex rides its angular
//! bisector, and each event — two neighbouring bisectors meeting — retires
//! an edge and starts a new skeleton branch. Convexity is what makes this
//! exact: no reflex vertex, so no split events, so every branch is a
//! straight segment between circumcentre-like meets.
//!
//! Everything else is refused by name: a face with holes, a reflex corner,
//! a curved boundary. Each of those changes the mathematics — holes and
//! reflex corners introduce split events, arcs introduce parabolic
//! bisectors — and a wrong axis is worse than a named refusal, because tool
//! paths gouge quietly.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Point, Point2, Vector2};
use ogeom_topo::{EdgeRepr, Model, NodeData, Shape};

/// The medial axis of a face, as straight segments between branch points.
#[derive(Debug, Clone)]
pub struct MedialAxis {
    /// The skeleton's segments, each an inward branch: from a boundary
    /// vertex or an earlier meet, to a meet or the centre.
    pub segments: Vec<(Point, Point)>,
    /// The inscribed-circle radius at each segment's inner end — the
    /// clearance a tool of that radius has there.
    pub clearance: Vec<f64>,
}

/// The medial axis of a convex planar polygonal face.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction), each
/// by name: a face that is not planar, carries inner wires, has a curved
/// edge, or turns a reflex corner. Those need split events and parabolic
/// bisectors this construction does not have; refusing is the honest answer
/// until it does.
pub fn medial_axis(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<MedialAxis> {
    let Some(data) = model.node(face).and_then(|n| match n.data() {
        NodeData::Face(d) => Some(d.clone()),
        _ => None,
    }) else {
        ogeom_bail!(Construction, "the shape is not a face");
    };
    let Some(ogeom_geom::SurfaceGeometry::Plane(plane)) = model.geometry().surface(data.surface)
    else {
        ogeom_bail!(
            Construction,
            "the medial axis is computed for planar faces; this face's \
             surface is not a plane"
        );
    };
    let frame = plane.plane().frame();
    let placement = face.transform(model.datums())?;

    let wires = model.ordered_children_of(face)?;
    if wires.len() != 1 {
        ogeom_bail!(
            Construction,
            "a face with inner wires needs split events; the medial axis of \
             a region with holes is not constructed yet"
        );
    }
    let mut ring: Vec<Point2> = Vec::new();
    for edge in model.ordered_children_of(&wires[0])? {
        let Some(edge_data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = edge_data.curve3d() else {
            ogeom_bail!(Construction, "a boundary edge carries no curve");
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            ogeom_bail!(Construction, "a boundary curve is not in this model");
        };
        if !matches!(geometry, ogeom_geom::Curve::Line(_)) {
            ogeom_bail!(
                Construction,
                "a curved boundary bisects along parabolas; the medial axis \
                 of arcs is not constructed yet"
            );
        }
        let (t0, t1) = if edge.orientation() == ogeom_topo::Orientation::Reversed {
            (range.1, range.0)
        } else {
            (range.0, range.1)
        };
        let start = placement.apply(geometry.point_at(t0, tol)?);
        let _ = t1;
        let local = frame.to_local(start);
        ring.push(Point2::new(local.x, local.y));
    }
    if ring.len() < 3 {
        ogeom_bail!(Construction, "a polygon needs three corners");
    }
    // Wind counter-clockwise, so inward is to the left.
    if signed_area(&ring) < 0.0 {
        ring.reverse();
    }
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        let c = ring[(i + 2) % ring.len()];
        let cross = (b - a).cross(c - b);
        if cross < -tol.confusion() {
            ogeom_bail!(
                Construction,
                "a reflex corner needs split events; the medial axis of a \
                 non-convex region is not constructed yet"
            );
        }
    }

    // The shrink loop. Each active vertex rides the bisector of its two
    // neighbouring edges; the earliest meeting of adjacent riders retires
    // the edge between them.
    let lift = |p: Point2, out: &mut Vec<(Point, Point)>, q: Point2| {
        let a = frame.origin() + frame.x().vector() * p.x + frame.y().vector() * p.y;
        let b = frame.origin() + frame.x().vector() * q.x + frame.y().vector() * q.y;
        out.push((a, b));
    };
    let mut segments: Vec<(Point, Point)> = Vec::new();
    let mut clearance: Vec<f64> = Vec::new();
    // Active loop: (position, the two edge directions meeting there).
    let n = ring.len();
    let mut active: Vec<(Point2, Vector2, Vector2, f64)> = Vec::with_capacity(n);
    for i in 0..n {
        let prev = ring[(i + n - 1) % n];
        let here = ring[i];
        let next = ring[(i + 1) % n];
        let e_in = (here - prev).normalized(tol)?;
        let e_out = (next - here).normalized(tol)?;
        active.push((here, e_in, e_out, 0.0));
    }

    while active.len() > 2 {
        // The earliest adjacent meet.
        let m = active.len();
        let mut best: Option<(usize, Point2, f64)> = None;
        for i in 0..m {
            let j = (i + 1) % m;
            let (pi, ei_in, ei_out, ti) = active[i];
            let (pj, ej_in, ej_out, tj) = active[j];
            let bi = bisector_dir(ei_in, ei_out);
            let bj = bisector_dir(ej_in, ej_out);
            let Some(meet) = ray_meet(pi, bi, pj, bj) else {
                continue;
            };
            // The event time is the inward distance of the shared edge —
            // both riders reach the meet as that edge's offset sweeps it.
            let speed_i = rider_speed(ei_in, ei_out);
            let t = ti + (meet - pi).magnitude() / speed_i;
            let _ = tj;
            let _ = ej_in;
            let _ = ej_out;
            if best.as_ref().is_none_or(|(_, _, held)| t < *held) {
                best = Some((i, meet, t));
            }
        }
        let Some((i, meet, t)) = best else {
            ogeom_bail!(
                Construction,
                "the shrink found no event; the ring is degenerate"
            );
        };
        let j = (i + 1) % active.len();
        let (pi, ei_in, _, _) = active[i];
        let (pj, _, ej_out, _) = active[j];
        lift(pi, &mut segments, meet);
        lift(pj, &mut segments, meet);
        clearance.push(t);
        clearance.push(t);
        // The two riders merge into one on the bisector of the surviving
        // outer edges.
        let merged = (meet, ei_in, ej_out, t);
        if i < j {
            active[i] = merged;
            active.remove(j);
        } else {
            active[j] = merged;
            active.remove(i);
        }
    }
    // Two riders left: they close on one another along the shared axis.
    if let [a, b] = active.as_slice() {
        lift(a.0, &mut segments, b.0);
        clearance.push(a.3.max(b.3));
    }
    Ok(MedialAxis {
        segments,
        clearance,
    })
}

/// The inward bisector direction at a corner between edge directions.
fn bisector_dir(e_in: Vector2, e_out: Vector2) -> Vector2 {
    let n_in = Vector2::new(-e_in.y, e_in.x);
    let n_out = Vector2::new(-e_out.y, e_out.x);
    (n_in + n_out)
        .normalized(Tolerances::millimetres())
        .unwrap_or(n_in)
}

/// How fast a rider moves along its bisector per unit of inward offset.
fn rider_speed(e_in: Vector2, e_out: Vector2) -> f64 {
    // The bisector makes angle θ/2 with each edge normal, where θ is the
    // turn; unit inward speed of the edges means 1/cos(θ/2) along it — but
    // 1/sin(half interior angle) in edge terms. Derived from the offset of
    // both edges staying on the rider.
    let n_in = Vector2::new(-e_in.y, e_in.x);
    let b = bisector_dir(e_in, e_out);
    let denom = b.dot(n_in).max(1e-12);
    1.0 / denom
}

/// Where two rays meet, `None` when parallel or behind either origin.
fn ray_meet(p: Point2, d: Point2Dir, q: Point2, e: Point2Dir) -> Option<Point2> {
    let denom = d.cross(e);
    if denom.abs() < 1e-14 {
        return None;
    }
    let w = q - p;
    let t = w.cross(e) / denom;
    let s = w.cross(d) / denom;
    if t < -1e-9 || s < -1e-9 {
        return None;
    }
    Some(p + d * t)
}

type Point2Dir = Vector2;

fn signed_area(ring: &[Point2]) -> f64 {
    let mut sum = 0.0;
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        sum += a.x * b.y - b.x * a.y;
    }
    sum / 2.0
}
