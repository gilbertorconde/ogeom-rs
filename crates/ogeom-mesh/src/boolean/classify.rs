//! The winding number of a closed mesh about a point, exactly.
//!
//! A ray from the point along a coordinate axis crosses the mesh's
//! triangles; each crossing counts `+1` where the triangle faces along the
//! ray and `-1` where it faces back, and the sum is the winding number: 1
//! inside a closed, outward-wound mesh and 0 outside it.
//!
//! A ray through an edge or a vertex of the mesh is resolved by simulation
//! of simplicity: the point is moved by `(ε, ε²)` in the two coordinates
//! across the ray, `ε` infinitesimal. The moved point lies on no edge, so
//! every crossing is a proper one, and the sign of each test on it is read
//! from the leading term that does not vanish. Every test is exact.

use ogeom_core::{OgeomResult, Sign, ogeom_bail};
use ogeom_math::{Aabb, Point};

use super::exact::times;
use super::pairs::Tri;
use super::points::{Def, Hom, Pt, axes, normal_sign, orient2d, side};
use crate::Bvh;

/// The sign of `a - b`.
fn sign_of_difference(a: f64, b: f64) -> Sign {
    if a > b {
        Sign::Positive
    } else if a < b {
        Sign::Negative
    } else {
        Sign::Zero
    }
}

/// Where the ray from `p` along `axis`, the way `dir` points, crosses
/// triangle `t`: `None` when it does not, else the crossing's count.
///
/// # Errors
///
/// [`ogeom_core::OgeomError::Invariant`] if `p` lies on `t`.
fn crossing(t: &Tri, p: &impl Hom, axis: usize, dir: Sign) -> OgeomResult<Option<Sign>> {
    let facing = normal_sign(&t.p, axis);
    if facing.is_zero() {
        return Ok(None);
    }
    let (i, j) = axes(axis);
    let corner = t.p.map(|q| Pt::new(Def::Vertex(q)));
    for k in 0..3 {
        let (a, b) = (&t.p[k], &t.p[(k + 1) % 3]);
        let mut o = orient2d(axis, &corner[k], &corner[(k + 1) % 3], p);
        if o.is_zero() {
            o = sign_of_difference(a[j], b[j]);
        }
        if o.is_zero() {
            o = sign_of_difference(b[i], a[i]);
        }
        if o != facing {
            return Ok(None);
        }
    }
    let s = side(&t.p, p);
    if s.is_zero() {
        ogeom_bail!(Invariant, "a point to classify lies on the other mesh");
    }
    // Ahead along the ray when the point is behind the plane as the ray
    // sees it.
    let along = times(dir, facing);
    Ok((s == along.reversed()).then_some(along))
}

/// The winding number of the mesh whose triangles are `tris`, indexed by
/// `bvh`, about point `p`, counted along `axis` the way `dir` points,
/// leaving out triangle `skip`.
///
/// # Errors
///
/// [`ogeom_core::OgeomError::Invariant`] if `p` lies on a triangle
/// counted.
pub(super) fn winding(
    tris: &[Tri],
    bvh: &Bvh,
    p: &impl Hom,
    axis: usize,
    dir: Sign,
    skip: Option<u32>,
) -> OgeomResult<i64> {
    let h = p.iv();
    let bounds = [h[0].over(&h[3]), h[1].over(&h[3]), h[2].over(&h[3])];
    let candidates = match bounds {
        [Some(x), Some(y), Some(z)] => {
            let (lo, hi) = ([x.lo, y.lo, z.lo], [x.hi, y.hi, z.hi]);
            bvh.triangles_where(|b: &Aabb| {
                let (Some(bl), Some(bh)) = (b.low(), b.high()) else {
                    return false;
                };
                let (bl, bh) = (coords(bl), coords(bh));
                let (i, j) = axes(axis);
                let across = bl[i] <= hi[i] && lo[i] <= bh[i] && bl[j] <= hi[j] && lo[j] <= bh[j];
                let ahead = if dir == Sign::Positive {
                    bh[axis] >= lo[axis]
                } else {
                    bl[axis] <= hi[axis]
                };
                across && ahead
            })
        }
        _ => (0..u32::try_from(tris.len()).unwrap_or(u32::MAX)).collect(),
    };
    let mut total = 0_i64;
    for t in candidates {
        if Some(t) == skip {
            continue;
        }
        match crossing(&tris[t as usize], p, axis, dir)? {
            Some(Sign::Positive) => total += 1,
            Some(Sign::Negative) => total -= 1,
            _ => {}
        }
    }
    Ok(total)
}

const fn coords(p: Point) -> [f64; 3] {
    [p.x, p.y, p.z]
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::super::points::Mid;
    use super::*;
    use ogeom_topo::Triangulation;

    fn cube() -> (Vec<Tri>, Bvh) {
        let positions: Vec<Point> = (0..8)
            .map(|k| {
                let bit = |b: u32| f64::from((k >> b) & 1);
                Point::new(bit(0), bit(1), bit(2))
            })
            .collect();
        let quads = [
            [0, 2, 3, 1],
            [4, 5, 7, 6],
            [0, 1, 5, 4],
            [2, 6, 7, 3],
            [0, 4, 6, 2],
            [1, 3, 7, 5],
        ];
        let mut mesh = Triangulation::new();
        mesh.positions = positions;
        for [a, b, c, d] in quads {
            mesh.triangles.push([a, b, c]);
            mesh.triangles.push([a, c, d]);
        }
        let tris = mesh
            .triangles
            .iter()
            .enumerate()
            .map(|(id, &v)| Tri {
                side: 0,
                id: u32::try_from(id).unwrap(),
                v,
                p: v.map(|k| coords(mesh.positions[k as usize])),
            })
            .collect();
        (tris, Bvh::build(&mesh).unwrap())
    }

    #[test]
    fn winding_counts_through_edges_and_corners() {
        let (tris, bvh) = cube();
        let at = |p: [f64; 3]| Pt::new(Def::Vertex(p));
        // Rays along the cube's diagonals of faces and through its corners.
        for (p, want) in [
            ([0.5, 0.5, 0.5], 1),
            ([0.25, 0.25, 0.25], 1),
            ([0.5, 0.0 + 0.5, 0.75], 1),
            ([-1.0, 0.0, 0.0], 0),
            ([-1.0, 1.0, 1.0], 0),
            ([2.0, 0.5, 0.5], 0),
        ] {
            for axis in 0..3 {
                for dir in [Sign::Positive, Sign::Negative] {
                    let w = winding(&tris, &bvh, &at(p), axis, dir, None).unwrap();
                    assert_eq!(w, want, "{p:?} axis {axis} {dir:?}");
                }
            }
        }
        // A point on the cube is refused, not counted.
        assert!(winding(&tris, &bvh, &at([0.0, 0.5, 0.5]), 0, Sign::Positive, None).is_err());
        // From a face's centre, outward with that face left out, nothing.
        let c = tris[0].p.map(|q| Pt::new(Def::Vertex(q)));
        let mid = Mid([&c[0], &c[1], &c[2]]);
        let w = winding(&tris, &bvh, &mid, 2, Sign::Negative, Some(0)).unwrap();
        assert_eq!(w, 0);
    }
}
