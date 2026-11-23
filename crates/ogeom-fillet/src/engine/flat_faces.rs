//! Closed-form reads of a face that the rolled, carried and stepped rounds
//! share, with the setbacks and seats that ask the same questions: a planar
//! face's plane and straight sides, whether a point lies inside it and how
//! far one is from it, whether a segment enters a polygon, and a face's
//! placed surface.

use crate::engine::consumed_seat::segment;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::SurfaceGeometry;
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, NodeData, Orientation, Shape, ShapeType, explore_unique};

/// A planar face's plane as a point and its outward unit normal, and its
/// boundary as straight segments; `None` where it is curved or bounded by
/// anything else.
#[derive(Clone)]
pub(crate) struct Flat {
    pub(crate) origin: Point,
    pub(crate) normal: Vector,
    pub(crate) sides: Vec<[Point; 2]>,
}

pub(crate) fn flat(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Option<Flat>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    let Some(SurfaceGeometry::Plane(plane)) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let placement = face.transform(model.datums())?;
    let mut normal = placement.apply_vector(plane.plane().normal().vector());
    if face.orientation() == Orientation::Reversed {
        normal = -normal;
    }
    let origin = placement.apply(plane.plane().origin());
    let mut sides = Vec::new();
    for e in explore_unique(model, face, ShapeType::Edge)? {
        let Some(points) = segment(model, &e, tol)? else {
            return Ok(None);
        };
        sides.push(points);
    }
    Ok(Some(Flat {
        origin,
        normal,
        sides,
    }))
}

/// A unit basis of the plane square to `normal`.
pub(crate) fn basis(normal: Vector) -> (Vector, Vector) {
    let helper = if normal.x.abs() < 0.9 {
        Vector::new(1.0, 0.0, 0.0)
    } else {
        Vector::new(0.0, 1.0, 0.0)
    };
    let e1 = normal.cross(helper);
    let e1 = e1 / e1.magnitude();
    (e1, normal.cross(e1))
}

/// Whether `p`, in the plane of `sides`, lies inside the region they
/// bound, read by the parity of the sides a ray along `e1` crosses.
pub(crate) fn inside(sides: &[[Point; 2]], p: Point, (e1, e2): (Vector, Vector)) -> bool {
    let at = |q: Point| ((q - p).dot(e1), (q - p).dot(e2));
    let mut odd = false;
    for [a, b] in sides {
        let ((ax, ay), (bx, by)) = (at(*a), at(*b));
        if (ay > 0.0) != (by > 0.0) {
            let x = ax + (bx - ax) * (0.0 - ay) / (by - ay);
            if x > 0.0 {
                odd = !odd;
            }
        }
    }
    odd
}

pub(crate) fn segment_distance(p: Point, [a, b]: [Point; 2]) -> f64 {
    let ab = b - a;
    let l2 = ab.dot(ab);
    let t = if l2 > 0.0 {
        ((p - a).dot(ab) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    p.distance(a + ab * t)
}

/// The distance from `p` to a planar face.
pub(crate) fn face_distance(face: &Flat, p: Point) -> f64 {
    let height = (p - face.origin).dot(face.normal);
    let foot = p - face.normal * height;
    if inside(&face.sides, foot, basis(face.normal)) {
        return height.abs();
    }
    face.sides
        .iter()
        .map(|s| segment_distance(p, *s))
        .fold(f64::INFINITY, f64::min)
}

/// The greatest of `a cos s + b sin s` over `[lo, hi]`.
pub(crate) fn trig_max(a: f64, b: f64, lo: f64, hi: f64) -> f64 {
    let (lo, hi) = (lo.min(hi), lo.max(hi));
    let m = a.hypot(b);
    let phase = b.atan2(a);
    let tau = core::f64::consts::TAU;
    let peak = phase + tau * ((lo - phase) / tau).ceil();
    if peak <= hi {
        return m;
    }
    let f = |s: f64| a * s.cos() + b * s.sin();
    f(lo).max(f(hi))
}

/// Whether some stretch of the segment `side` lies inside the simple
/// polygon `ring` in the plane square to `normal`, more than `margin` from
/// its boundary.
pub(crate) fn enters(ring: &[Point], normal: Vector, [p, q]: [Point; 2], margin: f64) -> bool {
    let frame = basis(normal);
    let sides: Vec<[Point; 2]> = (0..ring.len())
        .map(|i| [ring[i], ring[(i + 1) % ring.len()]])
        .collect();
    let (e1, e2) = frame;
    let flat = |x: Point| ((x - p).dot(e1), (x - p).dot(e2));
    let (qx, qy) = flat(q);
    let mut cuts = vec![0.0, 1.0];
    for [a, b] in &sides {
        let ((ax, ay), (bx, by)) = (flat(*a), flat(*b));
        let (dx, dy) = (bx - ax, by - ay);
        let den = qx * dy - qy * dx;
        if den == 0.0 {
            continue;
        }
        let t = (ax * dy - ay * dx) / den;
        let u = (ax * qy - ay * qx) / den;
        if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
            cuts.push(t);
        }
    }
    cuts.sort_by(f64::total_cmp);
    for pair in cuts.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        if b - a <= 0.0 {
            continue;
        }
        for k in 1..8 {
            let x = p + (q - p) * (a + (b - a) * f64::from(k) / 8.0);
            let clear = sides.iter().all(|s| segment_distance(x, *s) > margin);
            if clear && inside(&sides, x, frame) {
                return true;
            }
        }
    }
    false
}

/// A face's placed surface, with its outward sign.
pub(crate) fn host(
    model: &Model,
    face: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<(SurfaceGeometry, f64)>> {
    use ogeom_geom::Transformable as _;
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let placed = surface
        .clone()
        .transformed(&face.transform(model.datums())?, tol)?;
    let sign = if face.orientation() == Orientation::Reversed {
        -1.0
    } else {
        1.0
    };
    Ok(Some((placed, sign)))
}
