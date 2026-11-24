//! The convex hull of a few points in space, built incrementally with
//! exact orientation tests.
//!
//! Each point is tested against every face. A point strictly outside some
//! face is outside the hull; it then sees every face it is not strictly
//! inside, those in its plane included, and they are removed and the
//! horizon they leave closed with new faces through the point. A point in
//! a face's plane and outside the hull is outside that face, so removing
//! the face keeps the new hull convex, and no new face has its point on
//! the line of its horizon edge: both faces along that edge share the
//! line's plane, and both are seen. The hull may keep neighbouring faces
//! in one plane. Quadratic in the points, for the few a working subset
//! holds.

use ogeom_core::predicates::{Exact, Predicates, Sign};
use ogeom_math::Vector;

/// The hull's triangles, by index into the points, counter-clockwise seen
/// from outside; `None` where the points lie in one plane.
pub(super) fn hull(points: &[Vector]) -> Option<Vec<[usize; 3]>> {
    let at = |i: usize| points[i].to_array();
    let orient =
        |a: usize, b: usize, c: usize, d: usize| Exact::orient3d(at(a), at(b), at(c), at(d));
    let first = start(points)?;
    let mut faces: Vec<[usize; 3]> = Vec::new();
    let [a, b, c, d] = first;
    for (face, opposite) in [
        ([a, b, c], d),
        ([a, b, d], c),
        ([a, c, d], b),
        ([b, c, d], a),
    ] {
        // The opposite corner is inside: below, where the face winds
        // clockwise, so the face winds counter-clockwise from outside.
        let [p, q, r] = face;
        if orient(p, q, r, opposite) == Sign::Positive {
            faces.push([p, q, r]);
        } else {
            faces.push([p, r, q]);
        }
    }
    for p in 0..points.len() {
        if first.contains(&p) {
            continue;
        }
        let sides: Vec<Sign> = faces.iter().map(|&[a, b, c]| orient(a, b, c, p)).collect();
        if !sides.contains(&Sign::Negative) {
            continue;
        }
        let seen: Vec<bool> = sides.iter().map(|s| *s != Sign::Positive).collect();
        let mut rim: Vec<(usize, usize)> = Vec::new();
        for (face, _) in faces.iter().zip(&seen).filter(|(_, s)| **s) {
            rim.extend([(face[0], face[1]), (face[1], face[2]), (face[2], face[0])]);
        }
        let horizon: Vec<(usize, usize)> = rim
            .iter()
            .copied()
            .filter(|&(u, v)| !rim.contains(&(v, u)))
            .collect();
        let mut keep = seen.iter();
        faces.retain(|_| keep.next().is_some_and(|s| !s));
        faces.extend(horizon.into_iter().map(|(u, v)| [u, v, p]));
    }
    Some(faces)
}

/// Four points not in one plane, the first of them the first point.
fn start(points: &[Vector]) -> Option<[usize; 4]> {
    let at = |i: usize| points[i].to_array();
    let a = 0;
    let b = (1..points.len()).find(|&i| points[i] != points[a])?;
    for c in 0..points.len() {
        for d in 0..points.len() {
            if Exact::orient3d(at(a), at(b), at(c), at(d)) != Sign::Zero {
                return Some([a, b, c, d]);
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A cube's corners and points on its faces, edges and inside: twelve
    /// triangles' worth of faces, every corner a vertex and every face's
    /// normal outward.
    #[test]
    fn a_cubes_hull_has_its_corners() {
        let mut points = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    points.push(Vector::new(f64::from(i), f64::from(j), f64::from(k)));
                }
            }
        }
        let faces = hull(&points).unwrap_or_else(|| panic!("a hull"));
        let centre = Vector::new(1.0, 1.0, 1.0);
        let mut area = 0.0;
        for &[a, b, c] in &faces {
            let n = (points[b] - points[a]).cross(points[c] - points[a]);
            assert!(n.dot(points[a] - centre) > 0.0, "{a} {b} {c} faces inward");
            area += n.magnitude() / 2.0;
        }
        assert!((area - 24.0).abs() < 1e-12, "{area}");
        let mut used: Vec<usize> = faces.iter().flatten().copied().collect();
        used.sort_unstable();
        used.dedup();
        for corner in [0, 2, 6, 8, 18, 20, 24, 26] {
            assert!(used.contains(&corner), "corner {corner} missing");
        }
        assert!(!used.contains(&13), "the centre is not on the hull");
    }

    #[test]
    fn points_in_one_plane_have_no_hull() {
        let points: Vec<Vector> = (0..10)
            .map(|i| Vector::new(f64::from(i), f64::from(i * i % 7), 0.0))
            .collect();
        assert!(hull(&points).is_none());
    }
}
