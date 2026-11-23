//! The convex hull of a point cloud in space, by quickhull with exact
//! orientation tests.
//!
//! Each face keeps the points strictly outside it that no earlier face
//! took. The face with points left takes its furthest one as the eye; the
//! eye sees every face it is not strictly inside, found by walking from
//! that face across its neighbours, those in its plane included (the
//! argument in `fit/zone/hull.rs` for why that keeps the hull convex holds
//! here unchanged). The horizon the seen faces leave is closed with new
//! faces through the eye, and the points the seen faces kept go to the
//! first new face they are strictly outside, or are dropped as inside.
//! Neighbouring faces may lie in one plane.
//!
//! The walk and the hand-over of points keep the work near `n log n` for
//! the clouds a mesh gives, where the subset hull in `fit/zone` tests every
//! point against every face.

use ogeom_core::predicates::{Exact, Predicates, Sign};
use ogeom_core::{FastMap, FastSet};
use ogeom_math::{Point, Vector};

/// A closed convex hull: triangles by index into the points it was built
/// from, counter-clockwise seen from outside, and across each triangle's
/// edge `k` (corner `k` to corner `k + 1`) the triangle on the other side.
pub(crate) struct Hull {
    pub(crate) faces: Vec<[u32; 3]>,
    pub(crate) across: Vec<[usize; 3]>,
}

struct Face {
    corners: [u32; 3],
    across: [usize; 3],
    alive: bool,
    outside: Vec<u32>,
    normal: Vector,
}

/// The hull of `points`, or `None` where they lie in one plane, or where
/// the horizon the eye sees is not one loop (it always is for points in
/// general position; the walk only fails on a degenerate input the exact
/// tests cannot settle into a convex shape).
pub(crate) fn convex_hull(points: &[Point]) -> Option<Hull> {
    let at = |i: u32| points[i as usize].to_array();
    let orient = |f: &[u32; 3], p: u32| Exact::orient3d(at(f[0]), at(f[1]), at(f[2]), at(p));
    let mut seen_point: FastSet<[u64; 3]> = FastSet::default();
    let unique: Vec<u32> = (0..u32::try_from(points.len()).ok()?)
        .filter(|&i| {
            let p = points[i as usize];
            p.x.is_finite()
                && p.y.is_finite()
                && p.z.is_finite()
                && seen_point.insert([p.x.to_bits(), p.y.to_bits(), p.z.to_bits()])
        })
        .collect();
    let [a, b, c, d] = start(points, &unique)?;

    let plane = |corners: [u32; 3]| {
        let [p, q, r] = corners.map(|i| points[i as usize]);
        let n = (q - p).cross(r - p);
        let length = n.magnitude();
        if length > 0.0 { n / length } else { n }
    };
    let mut faces: Vec<Face> = Vec::new();
    for (corners, opposite) in [
        ([a, b, c], d),
        ([a, b, d], c),
        ([a, c, d], b),
        ([b, c, d], a),
    ] {
        let corners = if orient(&corners, opposite) == Sign::Positive {
            corners
        } else {
            [corners[0], corners[2], corners[1]]
        };
        let normal = plane(corners);
        faces.push(Face {
            corners,
            across: [usize::MAX; 3],
            alive: true,
            outside: Vec::new(),
            normal,
        });
    }
    let mut by_edge: FastMap<(u32, u32), usize> = FastMap::default();
    for (f, face) in faces.iter().enumerate() {
        for k in 0..3 {
            by_edge.insert((face.corners[k], face.corners[(k + 1) % 3]), f);
        }
    }
    for face in &mut faces {
        for k in 0..3 {
            face.across[k] = by_edge[&(face.corners[(k + 1) % 3], face.corners[k])];
        }
    }
    for &p in &unique {
        if [a, b, c, d].contains(&p) {
            continue;
        }
        if let Some(face) = faces
            .iter_mut()
            .find(|f| orient(&f.corners, p) == Sign::Negative)
        {
            face.outside.push(p);
        }
    }

    let mut pending: Vec<usize> = (0..faces.len()).collect();
    let mut mark: Vec<u32> = vec![0; faces.len()];
    let mut epoch = 0_u32;
    while let Some(f) = pending.pop() {
        if !faces[f].alive || faces[f].outside.is_empty() {
            continue;
        }
        let height = |p: u32| faces[f].normal.dot(points[p as usize] - Point::ORIGIN);
        let eye = faces[f]
            .outside
            .iter()
            .copied()
            .max_by(|&p, &q| height(p).total_cmp(&height(q)).then(q.cmp(&p)))?;

        epoch += 1;
        let mut seen: Vec<usize> = vec![f];
        let mut walk: Vec<usize> = vec![f];
        mark[f] = epoch;
        // Each horizon edge, from a seen face: its two corners in the seen
        // face's winding and the unseen face across it.
        let mut horizon: Vec<(u32, u32, usize)> = Vec::new();
        while let Some(g) = walk.pop() {
            for k in 0..3 {
                let h = faces[g].across[k];
                if mark[h] == epoch {
                    continue;
                }
                if orient(&faces[h].corners, eye) == Sign::Positive {
                    let corners = faces[g].corners;
                    horizon.push((corners[k], corners[(k + 1) % 3], h));
                } else {
                    mark[h] = epoch;
                    seen.push(h);
                    walk.push(h);
                }
            }
        }
        // A face first found behind a horizon edge may be seen later from
        // another side: that edge is then inside the seen region.
        horizon.retain(|&(_, _, h)| mark[h] != epoch);

        let first_new = faces.len();
        let mut from: FastMap<u32, usize> = FastMap::default();
        let mut to: FastMap<u32, usize> = FastMap::default();
        for (i, &(u, v, h)) in horizon.iter().enumerate() {
            let n = first_new + i;
            if from.insert(u, n).is_some() || to.insert(v, n).is_some() {
                return None;
            }
            let back = faces[h].corners.iter().position(|&x| x == v)?;
            if faces[h].corners[(back + 1) % 3] != u {
                return None;
            }
            faces[h].across[back] = n;
            let corners = [u, v, eye];
            let normal = plane(corners);
            faces.push(Face {
                corners,
                across: [h, usize::MAX, usize::MAX],
                alive: true,
                outside: Vec::new(),
                normal,
            });
            mark.push(0);
        }
        for face in &mut faces[first_new..] {
            let [u, v, _] = face.corners;
            face.across[1] = *from.get(&v)?;
            face.across[2] = *to.get(&u)?;
        }

        let mut orphans: Vec<u32> = Vec::new();
        for &g in &seen {
            faces[g].alive = false;
            orphans.append(&mut faces[g].outside);
        }
        for p in orphans {
            if p == eye {
                continue;
            }
            if let Some(n) =
                (first_new..faces.len()).find(|&n| orient(&faces[n].corners, p) == Sign::Negative)
            {
                faces[n].outside.push(p);
            }
        }
        pending.extend(first_new..faces.len());
    }

    let mut index: Vec<usize> = vec![usize::MAX; faces.len()];
    let mut kept = 0;
    for (f, face) in faces.iter().enumerate() {
        if face.alive {
            index[f] = kept;
            kept += 1;
        }
    }
    let live = faces.iter().filter(|f| f.alive);
    Some(Hull {
        faces: live.clone().map(|f| f.corners).collect(),
        across: live.map(|f| f.across.map(|g| index[g])).collect(),
    })
}

/// Four of the points not in one plane: the lowest in `x`, the one
/// furthest from it, the one furthest from their line and the one furthest
/// from the plane of those three, checked exactly, with a search where the
/// furthest are still coplanar.
fn start(points: &[Point], unique: &[u32]) -> Option<[u32; 4]> {
    let p = |i: u32| points[i as usize];
    let furthest = |score: &dyn Fn(Point) -> f64| {
        unique
            .iter()
            .copied()
            .max_by(|&i, &j| score(p(i)).total_cmp(&score(p(j))).then(j.cmp(&i)))
    };
    let a = unique
        .iter()
        .copied()
        .min_by(|&i, &j| p(i).x.total_cmp(&p(j).x).then(i.cmp(&j)))?;
    let b = furthest(&|q| (q - p(a)).magnitude())?;
    let along = p(b) - p(a);
    let c = furthest(&|q| along.cross(q - p(a)).magnitude())?;
    let normal = along.cross(p(c) - p(a));
    let d = furthest(&|q| normal.dot(q - p(a)).abs())?;
    let orient = |i: u32, j: u32, k: u32, l: u32| {
        Exact::orient3d(
            p(i).to_array(),
            p(j).to_array(),
            p(k).to_array(),
            p(l).to_array(),
        )
    };
    if orient(a, b, c, d) != Sign::Zero {
        return Some([a, b, c, d]);
    }
    for &c in unique {
        for &d in unique {
            if orient(a, b, c, d) != Sign::Zero {
                return Some([a, b, c, d]);
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::cast_precision_loss, reason = "test code")]
    use super::*;

    fn check_closed_and_convex(points: &[Point], hull: &Hull) {
        for (f, corners) in hull.faces.iter().enumerate() {
            for k in 0..3 {
                let g = hull.across[f][k];
                let (u, v) = (corners[k], corners[(k + 1) % 3]);
                let back = hull.faces[g].iter().position(|&x| x == v).unwrap();
                assert_eq!(
                    hull.faces[g][(back + 1) % 3],
                    u,
                    "{f} and {g} do not share {u}-{v}"
                );
                assert_eq!(hull.across[g][back], f);
            }
            let at = |i: u32| points[i as usize].to_array();
            for i in 0..u32::try_from(points.len()).unwrap() {
                let side = Exact::orient3d(at(corners[0]), at(corners[1]), at(corners[2]), at(i));
                assert_ne!(side, Sign::Negative, "point {i} outside face {f}");
            }
        }
    }

    /// A 3 × 3 × 3 grid of a cube: every corner a vertex, the points on
    /// its faces, edges and inside none, and the faces' area the cube's.
    #[test]
    fn a_grid_cubes_hull_is_its_six_squares() {
        let mut points = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    points.push(Point::new(f64::from(i), f64::from(j), f64::from(k)));
                }
            }
        }
        let hull = convex_hull(&points).unwrap();
        check_closed_and_convex(&points, &hull);
        let area: f64 = hull
            .faces
            .iter()
            .map(|&[a, b, c]| {
                let [a, b, c] = [a, b, c].map(|i| points[i as usize]);
                (b - a).cross(c - a).magnitude() / 2.0
            })
            .sum();
        assert!((area - 24.0).abs() < 1e-12, "{area}");
        let used: FastSet<u32> = hull.faces.iter().flatten().copied().collect();
        assert!(!used.contains(&13), "the centre is not on the hull");
    }

    /// Points on a sphere, every one on the hull: Euler's count of faces
    /// for a triangulated sphere, `2n − 4`.
    #[test]
    fn points_on_a_sphere_are_all_corners() {
        let n = 600;
        let golden = core::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
        let points: Vec<Point> = (0..n)
            .map(|i| {
                let z = 1.0 - 2.0 * (f64::from(i) + 0.5) / f64::from(n);
                let r = (1.0 - z * z).sqrt();
                let t = golden * f64::from(i);
                Point::new(r * t.cos(), r * t.sin(), z)
            })
            .collect();
        let hull = convex_hull(&points).unwrap();
        check_closed_and_convex(&points, &hull);
        assert_eq!(hull.faces.len(), 2 * n as usize - 4);
    }

    #[test]
    fn points_in_one_plane_have_no_hull() {
        let points: Vec<Point> = (0..10)
            .map(|i| Point::new(f64::from(i), f64::from(i * i % 7), 0.0))
            .collect();
        assert!(convex_hull(&points).is_none());
    }
}
