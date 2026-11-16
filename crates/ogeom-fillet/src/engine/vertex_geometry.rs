//! Half-space geometry shared by vertex-only tools and their adapters.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape};

/// A convex polytope from the half-spaces that bound it, each given by
/// its inward normal and a point on its plane.
///
/// Its corners are the feasible meetings of three planes; each plane's
/// face is those of its corners that lie on it, walked round the face's
/// centroid. The polyhedron builder checks what this hands it (planarity,
/// every edge shared by two faces) so a set of half-spaces that bounds
/// nothing, or bounds a sliver, is refused rather than built.
pub(crate) fn convex_block(
    model: &mut Model,
    walls: &[(Vector, Point)],
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let slack = tol.confusion() * 100.0;
    let n = walls.len();
    let mut points: Vec<Point> = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                let Some(p) = planes_meet(&walls[i], &walls[j], &walls[k]) else {
                    continue;
                };
                if walls
                    .iter()
                    .any(|(normal, on)| normal.dot(p - *on) < -slack)
                {
                    continue;
                }
                if points.iter().any(|q| q.distance(p) <= slack) {
                    continue;
                }
                points.push(p);
            }
        }
    }
    let mut rings: Vec<Vec<usize>> = Vec::new();
    for (normal, on) in walls {
        let mine: Vec<usize> = (0..points.len())
            .filter(|&i| normal.dot(points[i] - *on).abs() <= slack)
            .collect();
        if mine.len() < 3 {
            continue;
        }
        let centroid = mine
            .iter()
            .fold(Vector::ZERO, |acc, &i| acc + (points[i] - Point::ORIGIN))
            * (1.0 / f64::from(u32::try_from(mine.len()).unwrap_or(u32::MAX)));
        let centroid = Point::ORIGIN + centroid;
        let axis = normal.normalized(tol)?;
        let first = (points[mine[0]] - centroid).normalized(tol)?;
        let second = axis.cross(first);
        let mut ordered: Vec<(f64, usize)> = mine
            .iter()
            .map(|&i| {
                let v = points[i] - centroid;
                (v.dot(second).atan2(v.dot(first)), i)
            })
            .collect();
        ordered.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(core::cmp::Ordering::Equal));
        rings.push(ordered.into_iter().map(|(_, i)| i).collect());
    }
    Ok(ogeom_algo::make_polyhedron(model, &points, &rings, tol)?.shape)
}

/// Where three planes meet, or `None` where they do not span.
pub(crate) fn planes_meet(
    a: &(Vector, Point),
    b: &(Vector, Point),
    c: &(Vector, Point),
) -> Option<Point> {
    let m = [a.0, b.0, c.0];
    let rhs = [
        a.0.dot(a.1 - Point::ORIGIN),
        b.0.dot(b.1 - Point::ORIGIN),
        c.0.dot(c.1 - Point::ORIGIN),
    ];
    let x = solve3(&m, rhs)?;
    Some(Point::ORIGIN + x)
}

/// `m_k · x = rhs_k` for three rows, or `None` where they do not span.
pub(crate) fn solve3(m: &[Vector; 3], rhs: [f64; 3]) -> Option<Vector> {
    let a: [[f64; 3]; 3] = std::array::from_fn(|i| [m[i].x, m[i].y, m[i].z]);
    let det3 = |a: &[[f64; 3]; 3]| -> f64 {
        a[0][0].mul_add(
            a[1][1].mul_add(a[2][2], -(a[1][2] * a[2][1])),
            -a[0][1].mul_add(
                a[1][0].mul_add(a[2][2], -(a[1][2] * a[2][0])),
                -(a[0][2] * a[1][0].mul_add(a[2][1], -(a[1][1] * a[2][0]))),
            ),
        )
    };
    let det = det3(&a);
    if !det.is_finite() || det.abs() <= 1e-12 {
        return None;
    }
    let mut x = [0.0_f64; 3];
    for (k, xk) in x.iter_mut().enumerate() {
        let mut ak = a;
        for (row, r) in ak.iter_mut().zip(rhs) {
            row[k] = r;
        }
        *xk = det3(&ak) / det;
    }
    Some(Vector::new(x[0], x[1], x[2]))
}
