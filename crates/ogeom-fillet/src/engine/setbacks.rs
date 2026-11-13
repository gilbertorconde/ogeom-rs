//! P2: where blends end, by intersection, not by case.
//!
//! On each face, the contact lines of the blends on its two edges at a
//! vertex meet at the face's setback corner there: three planes, the face's
//! and the two lines' (each square to the face). An edge not blended keeps
//! its own line, so a corner of no setback stays where it is.

use ogeom_core::{FastMap, OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Point, Vector};

use super::classify::Planar;
use super::{Id, id};

/// Every face corner's setback corner (holes' corners too), keyed by face
/// index and vertex, or `None` where two edges of a face run on in one line
/// at a vertex and are set back by different distances there.
pub(crate) fn setback_corners(
    planar: &Planar,
    setback: &FastMap<(Id, Id), f64>,
    tol: Tolerances,
) -> OgeomResult<Option<FastMap<(usize, Id), Point>>> {
    let mut out = FastMap::default();
    for (fi, face) in planar.faces.iter().enumerate() {
        for ring in face.loops() {
            let n = ring.len();
            for k in 0..n {
                let (prev, v, next) = (ring[(k + n - 1) % n], ring[k], ring[(k + 1) % n]);
                // An edge's direction along the ring, the unit step into the
                // face square to it, and the face's setback from it.
                let leg = |a: Id, b: Id| -> OgeomResult<(Vector, Vector, f64)> {
                    let along = planar.points[&b] - planar.points[&a];
                    let length = along.magnitude();
                    if length <= tol.confusion() {
                        ogeom_bail!(Construction, "a face has an edge of no length");
                    }
                    let along = along / length;
                    let edge = planar.edges[&(a.min(b), a.max(b))];
                    let d = setback
                        .get(&(id(&face.shape), edge))
                        .copied()
                        .unwrap_or(0.0);
                    Ok((along, face.normal.cross(along), d))
                };
                let (u1, m1, d1) = leg(prev, v)?;
                let (u2, m2, d2) = leg(v, next)?;
                let at = planar.points[&v];
                let Some(offset) = corner_offset((u1, m1, d1), (u2, m2, d2), tol) else {
                    return Ok(None);
                };
                out.insert((fi, v), at + offset);
            }
        }
    }
    Ok(Some(out))
}

/// The step from a face's corner to its setback corner, in the face's
/// plane: `d1` in from the edge arriving along `u1` (stepping in along
/// `m1`) and `d2` in from the edge leaving along `u2` (along `m2`). An edge
/// not set back keeps the corner on its own line, so a corner of no setback
/// stays exactly where it is. `None` where the two edges run on in one line
/// with different setbacks.
fn corner_offset(
    (u1, m1, d1): (Vector, Vector, f64),
    (u2, m2, d2): (Vector, Vector, f64),
    tol: Tolerances,
) -> Option<Vector> {
    let turn = u1.cross(u2).magnitude();
    if turn <= tol.angular() {
        // Straight on: one line, one setback, or none to give.
        return ((d1 - d2).abs() <= tol.confusion()).then(|| m1 * d1);
    }
    if d1 == 0.0 && d2 == 0.0 {
        return Some(Vector::ZERO);
    }
    if d1 == 0.0 {
        // Along the arriving edge, back from the corner, to where the
        // leaving edge's setback line crosses it.
        return Some(u1 * (d2 / m2.dot(u1)));
    }
    if d2 == 0.0 {
        return Some(u2 * (d1 / m1.dot(u2)));
    }
    // In the plane: x = a u2 + b m2 with m1·x = d1 and m2·x = d2.
    let b = d2;
    let a = (d1 - b * m1.dot(m2)) / m1.dot(u2);
    Some(u2 * a + m2 * b)
}
