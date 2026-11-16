//! Holding a B-spline patch within a convex region of planes by its control
//! net.
//!
//! A B-spline whose weights are all positive lies in the convex hull of its
//! control points, and a plane's height is affine, so a net within every
//! plane of a convex region puts the whole surface within it.
//! [`net_held_within`] keeps the net's border rows as they are (they carry
//! the patch's sides and its cross-boundary derivatives, so its tangency
//! there) and moves each other point standing past a plane to the nearest
//! point of the region: the one minimum of the squared distance under the
//! planes' inequalities, found by trying the planes' active sets in a fixed
//! order. The planes are sorted first, so the result does not depend on the
//! order they are handed in.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::BSplineSurface;
use ogeom_math::{ControlGrid, Point, Vector, Weighted};

/// A fixed control point standing past a plane by more than the fit: no net
/// with that point kept lies within the region.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NetOutside {
    /// The plane's index in the list handed in.
    pub plane: usize,
    /// The point's indices in the net, along `u` and along `v`.
    pub at: (usize, usize),
    /// How far past the plane it stands.
    pub height: f64,
}

/// `surface` with its net held within the planes: every plane `(normal,
/// at)` bounds the points `p` with `normal · (p - at) <= 0`, `normal` of
/// unit length. The `fixed` outermost rows and columns of the net are kept,
/// each required within `fit` of every plane; every other point standing
/// past a plane moves to the nearest point within all of them. The surface
/// then lies within `fit` of the region (and within it where the kept rows
/// are). `Err(NetOutside)` names a kept point standing past a plane by more
/// than `fit`, and the region's nearest point where none is found (planes
/// with no common side).
///
/// # Errors
///
/// [`ogeom_core::OgeomError::Construction`] where a weight is not positive
/// (the hull bounds no such surface) or a normal is not of unit length.
pub fn net_held_within(
    surface: &BSplineSurface,
    planes: &[(Vector, Point)],
    fixed: usize,
    fit: f64,
) -> OgeomResult<Result<BSplineSurface, NetOutside>> {
    for (normal, _) in planes {
        if (normal.magnitude() - 1.0).abs() > 1e-9 {
            ogeom_bail!(Construction, "a plane's normal is not of unit length");
        }
    }
    // Sorted by normal, then offset: a fixed order whatever order the
    // planes come in, so the active sets are tried alike.
    let mut order: Vec<usize> = (0..planes.len()).collect();
    let key = |k: usize| {
        let (n, at) = planes[k];
        [n.x, n.y, n.z, n.dot(at.to_vector())]
    };
    order.sort_by(|&a, &b| {
        key(a)
            .iter()
            .zip(key(b).iter())
            .map(|(x, y)| x.total_cmp(y))
            .find(|o| o.is_ne())
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    let sorted: Vec<(Vector, Point)> = order.iter().map(|&k| planes[k]).collect();
    let height = |p: Point| -> (usize, f64) {
        sorted
            .iter()
            .enumerate()
            .map(|(k, (n, at))| (k, n.dot(p - *at)))
            .fold((0, f64::NEG_INFINITY), |a, b| if b.1 > a.1 { b } else { a })
    };

    let grid = surface.grid();
    let (nu, nv) = (grid.u_count(), grid.v_count());
    let mut points = grid.points().to_vec();
    for w in &points {
        if w.weight.is_nan() || w.weight <= 0.0 {
            ogeom_bail!(Construction, "a control weight is not positive");
        }
    }
    let scale = points
        .iter()
        .map(|w| w.point().to_vector().magnitude())
        .fold(1.0_f64, f64::max);
    for i in 0..nu {
        for j in 0..nv {
            let k = i * nv + j;
            let p = points[k].point();
            let (plane, h) = height(p);
            let kept = i < fixed || j < fixed || i + fixed >= nu || j + fixed >= nv;
            if kept {
                if h > fit {
                    return Ok(Err(NetOutside {
                        plane: order[plane],
                        at: (i, j),
                        height: h,
                    }));
                }
                continue;
            }
            if h <= 0.0 {
                continue;
            }
            let Some(q) = nearest_within(&sorted, p, scale) else {
                return Ok(Err(NetOutside {
                    plane: order[plane],
                    at: (i, j),
                    height: h,
                }));
            };
            points[k] = Weighted::new(q, points[k].weight, Tolerances::default())?;
        }
    }
    Ok(Ok(BSplineSurface::rational(
        surface.u_knots().clone(),
        surface.v_knots().clone(),
        ControlGrid::new(points, nu, nv)?,
    )?))
}

/// The point within every plane nearest `p`: the squared distance is
/// strictly convex and the region convex, so the minimum is one point, the
/// one whose active planes' multipliers are all non-negative and which no
/// other plane stands under. Active sets of up to three planes (all a point
/// of space can rest on independently) are tried in order of size, then of
/// index. `None` where no set of up to three gives a point within all.
fn nearest_within(planes: &[(Vector, Point)], p: Point, scale: f64) -> Option<Point> {
    let slack = scale * 1e-12;
    let within = |q: Point| planes.iter().all(|(n, at)| n.dot(q - *at) <= slack);
    let m = planes.len();
    let mut sets: Vec<Vec<usize>> = Vec::new();
    for a in 0..m {
        sets.push(vec![a]);
    }
    for a in 0..m {
        for b in a + 1..m {
            sets.push(vec![a, b]);
        }
    }
    for a in 0..m {
        for b in a + 1..m {
            for c in b + 1..m {
                sets.push(vec![a, b, c]);
            }
        }
    }
    for set in sets {
        // On every plane of the set: p - sum(l_s n_s), with the Gram system
        // G l = h, h the heights at p.
        let n: Vec<Vector> = set.iter().map(|&s| planes[s].0).collect();
        let h: Vec<f64> = set
            .iter()
            .map(|&s| planes[s].0.dot(p - planes[s].1))
            .collect();
        let Some(l) = solve_gram(&n, &h) else {
            continue;
        };
        if l.iter().any(|x| *x < 0.0) {
            continue;
        }
        let q = n.iter().zip(&l).fold(p, |q, (n, l)| q - *n * *l);
        if within(q) {
            return Some(q);
        }
    }
    None
}

/// The multipliers `l` with `sum_t (n_s · n_t) l_t = h_s`, by Cramer's rule
/// on the one to three planes; `None` where their normals are dependent.
fn solve_gram(n: &[Vector], h: &[f64]) -> Option<Vec<f64>> {
    let g = |a: usize, b: usize| n[a].dot(n[b]);
    match n.len() {
        1 => Some(vec![h[0] / g(0, 0)]),
        2 => {
            let det = g(0, 0) * g(1, 1) - g(0, 1) * g(1, 0);
            if det.abs() <= 1e-12 {
                return None;
            }
            Some(vec![
                (h[0] * g(1, 1) - g(0, 1) * h[1]) / det,
                (g(0, 0) * h[1] - h[0] * g(1, 0)) / det,
            ])
        }
        3 => {
            let m = [
                [g(0, 0), g(0, 1), g(0, 2)],
                [g(1, 0), g(1, 1), g(1, 2)],
                [g(2, 0), g(2, 1), g(2, 2)],
            ];
            let det3 = |m: [[f64; 3]; 3]| {
                m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
            };
            let det = det3(m);
            if det.abs() <= 1e-12 {
                return None;
            }
            let mut l = Vec::with_capacity(3);
            for c in 0..3 {
                let mut mc = m;
                for (r, row) in mc.iter_mut().enumerate() {
                    row[c] = h[r];
                }
                l.push(det3(mc) / det);
            }
            Some(l)
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;
    use ogeom_geom::Surface as _;
    use ogeom_math::KnotVector;

    const T: Tolerances = Tolerances::millimetres();

    /// A 6 x 6 net over the unit square at height `z`, its middle points
    /// raised to `peak`.
    fn bump(peak: f64) -> BSplineSurface {
        let mut points = Vec::new();
        for i in 0..6 {
            for j in 0..6 {
                let inner = (2..4).contains(&i) && (2..4).contains(&j);
                let z = if inner { peak } else { 0.0 };
                points.push(Point::new(f64::from(i) / 5.0, f64::from(j) / 5.0, z));
            }
        }
        let knots = KnotVector::clamped_uniform(3, 6).unwrap();
        BSplineSurface::new(
            knots.clone(),
            knots,
            &ControlGrid::new(points, 6, 6).unwrap(),
            T,
        )
        .unwrap()
    }

    fn surface_heights(s: &BSplineSurface, planes: &[(Vector, Point)]) -> f64 {
        let mut worst = f64::NEG_INFINITY;
        for a in 0..=40 {
            for b in 0..=40 {
                let p = s
                    .point_at(f64::from(a) / 40.0, f64::from(b) / 40.0, T)
                    .unwrap();
                for (n, at) in planes {
                    worst = worst.max(n.dot(p - *at));
                }
            }
        }
        worst
    }

    #[test]
    fn a_bump_past_a_roof_is_held_under_it_its_border_kept() {
        // Two roof planes meeting over the middle at 0.6, the bump's net
        // raised to 2 through them.
        let s = 0.5_f64.sqrt();
        let planes = [
            (Vector::new(s, 0.0, s), Point::new(0.5, 0.0, 0.6)),
            (Vector::new(-s, 0.0, s), Point::new(0.5, 0.0, 0.6)),
        ];
        let before = bump(2.0);
        assert!(surface_heights(&before, &planes) > 0.05);
        let held = net_held_within(&before, &planes, 2, 1e-9).unwrap().unwrap();
        assert!(surface_heights(&held, &planes) <= 1e-12);
        for i in 0..6 {
            for j in 0..6 {
                if !((2..4).contains(&i) && (2..4).contains(&j)) {
                    let k = i * 6 + j;
                    assert_eq!(held.grid().points()[k], before.grid().points()[k]);
                }
            }
        }
        // Either order of the planes, the same net to the bit.
        let turned = net_held_within(&before, &[planes[1], planes[0]], 2, 1e-9)
            .unwrap()
            .unwrap();
        assert_eq!(held.grid().points(), turned.grid().points());
    }

    #[test]
    fn a_held_point_moves_to_the_region_s_nearest_point() {
        // The middle points raised to 0.4 under the planes z <= 0.1 and
        // x <= 0.5: the nearest point within both of (0.6, y, 0.4) is
        // (0.5, y, 0.1).
        let planes = [
            (Vector::new(0.0, 0.0, 1.0), Point::new(0.0, 0.0, 0.1)),
            (Vector::new(1.0, 0.0, 0.0), Point::new(0.5, 0.0, 0.0)),
        ];
        let held = net_held_within(&bump(0.4), &planes, 2, 1.0)
            .unwrap()
            .unwrap();
        let p = held.grid().points()[3 * 6 + 2].point();
        assert!(
            (p.x - 0.5).abs() < 1e-15 && (p.z - 0.1).abs() < 1e-15,
            "{p:?}"
        );
        let q = held.grid().points()[2 * 6 + 2].point();
        assert!(
            (q.x - 0.4).abs() < 1e-15 && (q.z - 0.1).abs() < 1e-15,
            "{q:?}"
        );
    }

    #[test]
    fn a_kept_point_past_a_plane_is_named() {
        let planes = [(Vector::new(0.0, 0.0, 1.0), Point::new(0.0, 0.0, -0.1))];
        let out = net_held_within(&bump(0.4), &planes, 2, 1e-3)
            .unwrap()
            .unwrap_err();
        assert_eq!((out.plane, out.at), (0, (0, 0)));
        assert!((out.height - 0.1).abs() < 1e-12);
    }
}
