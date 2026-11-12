//! P5: guards, independent of how the result was built.
//!
//! Before the solid is made: no ring crosses itself, no face is turned
//! inside out, no two rings cross, and every ring is planar (a vertex
//! facet that is not is refused by name: its corner closes with a patch).
//! After: the kernel's `check`, and no two faces folded flat along an
//! edge (one plane, facing opposite ways), which neither a volume nor a
//! check of each face alone would see. A guard that fails is a refusal,
//! never a search for another construction.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::Origin;
use super::assemble::Ring;

/// A ring's normal by Newell's sum, twice its area long.
fn newell(points: &[Point], ring: &[usize]) -> Vector {
    let mut n = Vector::ZERO;
    let a = points[ring[0]];
    for k in 1..ring.len().saturating_sub(1) {
        n += (points[ring[k]] - a).cross(points[ring[k + 1]] - a);
    }
    n
}

/// Refuse a ring that crosses itself, a face turned inside out, or two
/// rings crossing.
pub(crate) fn rings_fit(points: &[Point], rings: &[Ring], tol: Tolerances) -> OgeomResult<()> {
    for (ring, origin, _) in rings {
        let n = newell(points, ring);
        let n = n / n.magnitude();
        let a = points[ring[0]];
        let off = ring
            .iter()
            .map(|i| (points[*i] - a).dot(n).abs())
            .fold(0.0, f64::max);
        if off > tol.confusion() * 10.0 {
            if origin.iter().all(|o| matches!(o, Origin::Vertex(_))) {
                ogeom_bail!(
                    Construction,
                    "a chamfer corner of {} bevels whose contact lines do not cross \
                         on one plane closes with a patch through them, which is not \
                         built yet; see docs/CAPABILITIES.md, fillet.chamfers",
                    ring.len()
                );
            }
            ogeom_bail!(
                Construction,
                "a chamfer face is not planar; a corner sits {off} off"
            );
        }
    }
    for (ring, _, facing) in rings {
        if let Some(facing) = facing
            && newell(points, ring).dot(*facing) <= 0.0
        {
            ogeom_bail!(
                Construction,
                "a chamfer turns a face inside out: the distance is past what \
                     the face holds"
            );
        }
    }
    let frames: Vec<(Vector, Vector, Vector, Point)> = rings
        .iter()
        .map(|(ring, _, _)| {
            let n = newell(points, ring);
            let n = n / n.magnitude();
            let a = points[ring[0]];
            let x = (points[ring[1]] - a).cross(n).cross(n);
            let x = -x / x.magnitude();
            (n, x, n.cross(x), a)
        })
        .collect();
    let flat = |ri: usize, p: Point| -> (f64, f64) {
        let (_, x, y, a) = frames[ri];
        ((p - a).dot(x), (p - a).dot(y))
    };
    for (ri, (ring, _, _)) in rings.iter().enumerate() {
        let pts: Vec<(f64, f64)> = ring.iter().map(|i| flat(ri, points[*i])).collect();
        if crosses_itself(&pts, tol.confusion()) {
            ogeom_bail!(
                Construction,
                "a chamfer's setbacks cross over on a face: the distance is past \
                     what the face holds"
            );
        }
    }
    // A segment of one ring passing through another's inside.
    for (ri, (ring, _, _)) in rings.iter().enumerate() {
        let (n, _, _, a) = frames[ri];
        let pts: Vec<(f64, f64)> = ring.iter().map(|i| flat(ri, points[*i])).collect();
        for (rj, (other, _, _)) in rings.iter().enumerate() {
            if ri == rj {
                continue;
            }
            for k in 0..other.len() {
                let (p, q) = (points[other[k]], points[other[(k + 1) % other.len()]]);
                let (sp, sq) = ((p - a).dot(n), (q - a).dot(n));
                if sp * sq >= 0.0 || sp.abs() <= tol.confusion() || sq.abs() <= tol.confusion() {
                    continue;
                }
                let x = p + (q - p) * (sp / (sp - sq));
                if strictly_inside(&pts, flat(ri, x), tol.confusion()) {
                    ogeom_bail!(
                        Construction,
                        "two faces of the chamfered solid cross: the distance is past \
                             what the solid holds"
                    );
                }
            }
        }
    }
    Ok(())
}

/// The built solid's own guards: the kernel's check, and no two faces folded
/// flat along an edge.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// either fails.
pub(crate) fn solid_fits(model: &Model, solid: &Shape, tol: Tolerances) -> OgeomResult<()> {
    let found = ogeom_algo::check(model, solid, tol)?;
    if !found.is_valid() {
        ogeom_bail!(
            Construction,
            "the blended solid fails its check: {:?}",
            found.problems
        );
    }
    for e in explore_unique(model, solid, ShapeType::Edge)? {
        let faces = ogeom_topo::ancestors_of(model, solid, &e, ShapeType::Face)?;
        if let [a, b] = faces.as_slice() {
            let (_, na) = ogeom_algo::face_normal(model, a, tol)?;
            let (_, nb) = ogeom_algo::face_normal(model, b, tol)?;
            if na.dot(nb) < -1.0 + tol.angular() {
                ogeom_bail!(
                    Construction,
                    "two faces of the blended solid lie folded flat on one another"
                );
            }
        }
    }
    Ok(())
}

/// Whether the polygon `pts` crosses itself: two edges not next to each
/// other meeting, or touching within `eps`.
fn crosses_itself(pts: &[(f64, f64)], eps: f64) -> bool {
    let n = pts.len();
    for i in 0..n {
        for j in (i + 1)..n {
            if j == i + 1 || (i == 0 && j == n - 1) {
                continue;
            }
            if segments_meet(pts[i], pts[(i + 1) % n], pts[j], pts[(j + 1) % n], eps) {
                return true;
            }
        }
    }
    false
}

fn segments_meet(a: (f64, f64), b: (f64, f64), c: (f64, f64), d: (f64, f64), eps: f64) -> bool {
    let cross = |o: (f64, f64), p: (f64, f64), q: (f64, f64)| {
        (p.0 - o.0) * (q.1 - o.1) - (p.1 - o.1) * (q.0 - o.0)
    };
    let (d1, d2) = (cross(c, d, a), cross(c, d, b));
    let (d3, d4) = (cross(a, b, c), cross(a, b, d));
    let scale = eps * ((b.0 - a.0).hypot(b.1 - a.1) + (d.0 - c.0).hypot(d.1 - c.1)).max(1.0);
    (d1 * d2 < -scale * scale || d1.abs() <= scale || d2.abs() <= scale)
        && (d3 * d4 < -scale * scale || d3.abs() <= scale || d4.abs() <= scale)
        && {
            // The two boxes overlap.
            a.0.min(b.0) <= c.0.max(d.0) + eps
                && c.0.min(d.0) <= a.0.max(b.0) + eps
                && a.1.min(b.1) <= c.1.max(d.1) + eps
                && c.1.min(d.1) <= a.1.max(b.1) + eps
        }
}

/// Whether `p` lies inside the polygon `pts`, farther than `eps` from its
/// edges.
fn strictly_inside(pts: &[(f64, f64)], p: (f64, f64), eps: f64) -> bool {
    let n = pts.len();
    let mut inside = false;
    for i in 0..n {
        let (a, b) = (pts[i], pts[(i + 1) % n]);
        let (ab, ap) = ((b.0 - a.0, b.1 - a.1), (p.0 - a.0, p.1 - a.1));
        let len2 = ab.0 * ab.0 + ab.1 * ab.1;
        let t = ((ap.0 * ab.0 + ap.1 * ab.1) / len2).clamp(0.0, 1.0);
        let near = (a.0 + ab.0 * t - p.0).hypot(a.1 + ab.1 * t - p.1);
        if near <= eps {
            return false;
        }
        if (a.1 > p.1) != (b.1 > p.1) {
            let x = a.0 + (p.1 - a.1) / (b.1 - a.1) * (b.0 - a.0);
            if p.0 < x {
                inside = !inside;
            }
        }
    }
    inside
}
