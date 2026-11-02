//! The curvature of a triangle mesh at its vertices, read from a local
//! quadric fitted to the vertices about each one.
//!
//! # The local quadric
//!
//! About a vertex `p`, the neighbourhood is the [`NEIGHBOURS`] vertices
//! nearest it within `size`, `p` itself among them. Its radius `r` is the
//! distance to the next vertex out, the first left out, or `size` where
//! fewer are that near. Each neighbour at distance `d` weighs
//! `(1 − (d/r)²)²`: one at `p`, falling smoothly to zero at `r`, so a
//! vertex entering or leaving the neighbourhood as the mesh moves changes
//! nothing abruptly.
//!
//! Each vertex has a normal: the sum of the area normals of the faces at
//! it, as they wind. The neighbourhood's normal `n` is the weighted sum of
//! its vertices' normals. Over the plane through `p` normal to `n`, with
//! coordinates `x`, `y` along two directions in it, the height of each
//! neighbour along `n` is fitted by a quadric,
//! `h(x, y) = a + b x + c y + d x² + e xy + f y²`, by weighted least
//! squares. The fit has linear terms, so a neighbourhood normal off the
//! surface's by a small angle tilts the frame and not the surface it finds.
//!
//! The quadric's surface above `p` is the point `p + a n`. Its curvature
//! there is that of the graph of `h` at the origin: the shape operator from
//! the first and second fundamental forms of `h`, exact for the quadric
//! however the frame is tilted. Its eigenvalues are the principal
//! curvatures, and its eigenvectors, lifted onto the graph, the principal
//! directions.
//!
//! A neighbourhood whose vertices do not fix the six coefficients, fewer
//! than six with weight or all on one conic of the plane, has no quadric:
//! the normal equations are taken as singular where a Cholesky pivot is no
//! more than the rounding of the largest diagonal entry.
//!
//! # Sign
//!
//! Curvature is positive where the surface bends away from its outward
//! normal, as everywhere on a sphere wound outward: a sphere of radius `R`
//! reads `1/R` for both principal curvatures.

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::{Point, Vector};
use ogeom_topo::Triangulation;

use crate::point_index::{Neighbour, PointIndex};

/// How many vertices nearest each one its local quadric is fitted to.
///
/// The quadric has six coefficients; its height at the vertex, which
/// smoothing moves the vertex to, is drawn from about this many
/// independent heights, so the noise in them falls about as the square
/// root of their number.
pub const NEIGHBOURS: usize = 100;

/// The curvature of a mesh's surface at one vertex.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Curvature {
    /// The mean of the principal curvatures.
    pub mean: f64,
    /// The product of the principal curvatures.
    pub gaussian: f64,
    /// The larger principal curvature.
    pub k1: f64,
    /// The smaller principal curvature.
    pub k2: f64,
    /// The unit direction along which the surface bends by `k1`.
    pub dir1: Vector,
    /// The unit direction along which the surface bends by `k2`: the
    /// surface normal crossed with `dir1`.
    pub dir2: Vector,
}

/// The curvature of `mesh` at each of its vertices, from the local quadric
/// over the vertices within `size` of it, as the [module
/// documentation](self) describes. Computed in parallel; the answers are
/// the same at any thread count.
///
/// `None` for a vertex no triangle uses, and for one whose neighbourhood
/// fixes no quadric. The mesh's own normals, if it has any, are not read:
/// the faces' winding gives the outward side.
///
/// # Errors
///
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if `size` is
/// not positive and finite or a position is not finite;
/// [`OgeomError::Range`](ogeom_core::OgeomError::Range) if a triangle names
/// a vertex the mesh does not have.
pub fn curvature(mesh: &Triangulation, size: f64) -> OgeomResult<Vec<Option<Curvature>>> {
    check_size(size)?;
    let count = mesh.positions.len();
    if let Some(&k) = mesh
        .triangles
        .iter()
        .flatten()
        .find(|&&k| k as usize >= count)
    {
        ogeom_bail!(Range, "a triangle names vertex {k} of {count}");
    }
    let normals = vertex_normals(&mesh.positions, mesh.triangles.iter().copied());
    let members: Vec<u32> = (0..count)
        .filter(|&v| normals[v].is_some())
        .map(|v| u32::try_from(v).unwrap_or(u32::MAX))
        .collect();
    let cloud = Cloud::build(&mesh.positions, &normals, members)?;
    let vertices: Vec<u32> = (0..count)
        .map(|v| u32::try_from(v).unwrap_or(u32::MAX))
        .collect();
    Ok(in_parallel(&vertices, |v| {
        normals[v as usize]?;
        let p = mesh.positions[v as usize];
        let near = cloud.near(p, size);
        cloud.fit(p, &near, |_| true)?.curvature()
    }))
}

/// Refuse a neighbourhood size that is not a positive finite length.
pub(crate) fn check_size(size: f64) -> OgeomResult<()> {
    if !size.is_finite() || size <= 0.0 {
        ogeom_bail!(Numeric, "a neighbourhood size of {size} is not a length");
    }
    Ok(())
}

/// `f` over `items` in parallel, in batches, answers in item order.
pub(crate) fn in_parallel<R: Send>(items: &[u32], f: impl Fn(u32) -> R + Sync) -> Vec<R> {
    let batches: Vec<&[u32]> = items.chunks(crate::bvh::BATCH).collect();
    ogeom_core::parallel::map_ordered(&batches, |_, batch| {
        batch.iter().map(|&v| f(v)).collect::<Vec<R>>()
    })
    .into_iter()
    .flatten()
    .collect()
}

/// The unit normal at each vertex: the sum of the area normals of the
/// faces at it, as they wind; `None` where no face, or none with area,
/// meets it.
pub(crate) fn vertex_normals(
    positions: &[Point],
    faces: impl Iterator<Item = [u32; 3]>,
) -> Vec<Option<Vector>> {
    let mut sum = vec![Vector::new(0.0, 0.0, 0.0); positions.len()];
    for [a, b, c] in faces {
        let (pa, pb, pc) = (
            positions[a as usize],
            positions[b as usize],
            positions[c as usize],
        );
        let n = (pb - pa).cross(pc - pa);
        for v in [a, b, c] {
            sum[v as usize] += n;
        }
    }
    sum.into_iter()
        .map(|n| {
            let length = n.magnitude();
            (length > 0.0 && length.is_finite()).then(|| n / length)
        })
        .collect()
}

/// The vertices a neighbourhood is drawn from, with their normals.
pub(crate) struct Cloud<'a> {
    index: PointIndex,
    /// The vertex at each point of the index.
    members: Vec<u32>,
    positions: &'a [Point],
    normals: &'a [Option<Vector>],
}

/// A vertex's neighbourhood: the vertices in it, by their slot in the
/// cloud, and its radius.
pub(crate) struct Near {
    pub(crate) found: Vec<Neighbour>,
    pub(crate) radius: f64,
}

impl<'a> Cloud<'a> {
    /// The cloud of `members`, vertices of `positions` with a normal.
    pub(crate) fn build(
        positions: &'a [Point],
        normals: &'a [Option<Vector>],
        members: Vec<u32>,
    ) -> OgeomResult<Self> {
        let points: Vec<Point> = members.iter().map(|&v| positions[v as usize]).collect();
        Ok(Self {
            index: PointIndex::build(&points)?,
            members,
            positions,
            normals,
        })
    }

    /// The vertex at a slot.
    pub(crate) fn vertex(&self, slot: u32) -> u32 {
        self.members[slot as usize]
    }

    /// The neighbourhood of `p`: the [`NEIGHBOURS`] vertices nearest it
    /// within `size`, its radius the distance to the next one out, or
    /// `size` where there is none.
    pub(crate) fn near(&self, p: Point, size: f64) -> Near {
        let mut found = self.index.k_nearest(p, NEIGHBOURS + 1, size);
        let radius = if found.len() > NEIGHBOURS {
            found.pop().map_or(size, |n| n.distance)
        } else {
            size
        };
        Near { found, radius }
    }

    /// The local quadric about `p` over the vertices of `near` that `keep`
    /// accepts, or `None` where they fix none.
    pub(crate) fn fit(&self, p: Point, near: &Near, keep: impl Fn(u32) -> bool) -> Option<Local> {
        let r = near.radius;
        if r.is_nan() || r <= 0.0 {
            return None;
        }
        let taken: Vec<(Point, Vector, f64)> = near
            .found
            .iter()
            .filter_map(|n| {
                let v = self.vertex(n.index);
                let q = n.distance / r;
                let w = (1.0 - q * q).powi(2);
                let normal = self.normals[v as usize]?;
                (q < 1.0 && keep(v)).then(|| (self.positions[v as usize], normal, w))
            })
            .collect();
        let normal = taken
            .iter()
            .fold(Vector::new(0.0, 0.0, 0.0), |acc, &(_, n, w)| acc + n * w);
        let length = normal.magnitude();
        if length.is_nan() || length <= 0.0 {
            return None;
        }
        let n = normal / length;
        let (t1, t2) = tangents(n);
        let mut m = [[0.0; 6]; 6];
        let mut rhs = [0.0; 6];
        for &(q, _, w) in &taken {
            let d = q - p;
            let (u, v, h) = (d.dot(t1) / r, d.dot(t2) / r, d.dot(n) / r);
            let row = [1.0, u, v, u * u, u * v, v * v];
            for i in 0..6 {
                rhs[i] += w * row[i] * h;
                for j in 0..=i {
                    m[i][j] += w * row[i] * row[j];
                }
            }
        }
        let c = cholesky_solve(&mut m, rhs)?;
        Some(Local {
            origin: p,
            frame: [t1, t2, n],
            height: c[0] * r,
            gradient: [c[1], c[2]],
            hessian: [2.0 * c[3] / r, c[4] / r, 2.0 * c[5] / r],
        })
    }
}

/// Two unit directions making a right-handed frame with the unit `n`.
fn tangents(n: Vector) -> (Vector, Vector) {
    let (ax, ay, az) = (n.x.abs(), n.y.abs(), n.z.abs());
    let axis = if ax <= ay && ax <= az {
        Vector::new(1.0, 0.0, 0.0)
    } else if ay <= az {
        Vector::new(0.0, 1.0, 0.0)
    } else {
        Vector::new(0.0, 0.0, 1.0)
    };
    let t1 = axis - n * axis.dot(n);
    let t1 = t1 / t1.magnitude();
    (t1, n.cross(t1))
}

/// Solve `m c = rhs` for the symmetric `m`, given by its lower triangle,
/// by Cholesky; `None` where a pivot is no more than the rounding of the
/// largest diagonal entry, as for a singular `m`.
fn cholesky_solve(m: &mut [[f64; 6]; 6], rhs: [f64; 6]) -> Option<[f64; 6]> {
    let largest = (0..6).map(|i| m[i][i]).fold(0.0_f64, f64::max);
    let floor = 6.0 * f64::EPSILON * largest;
    for j in 0..6 {
        let row = m[j];
        let pivot = row[..j].iter().fold(row[j], |acc, &x| acc - x * x);
        if pivot.is_nan() || pivot <= floor {
            return None;
        }
        let root = pivot.sqrt();
        m[j][j] = root;
        let above = m[j];
        for below in &mut m[j + 1..] {
            let s = below[..j]
                .iter()
                .zip(&above[..j])
                .fold(below[j], |acc, (&a, &b)| acc - a * b);
            below[j] = s / root;
        }
    }
    let mut y = rhs;
    for i in 0..6 {
        for k in 0..i {
            y[i] -= m[i][k] * y[k];
        }
        y[i] /= m[i][i];
    }
    for i in (0..6).rev() {
        for k in i + 1..6 {
            y[i] -= m[k][i] * y[k];
        }
        y[i] /= m[i][i];
    }
    y.iter().all(|c| c.is_finite()).then_some(y)
}

/// A local quadric height field over a vertex's tangent plane.
pub(crate) struct Local {
    origin: Point,
    /// The frame's two tangent directions and its normal.
    frame: [Vector; 3],
    /// The quadric's height above the origin.
    height: f64,
    /// Its slope there, along the two tangent directions.
    gradient: [f64; 2],
    /// Its second derivatives there: `h_xx`, `h_xy`, `h_yy`.
    hessian: [f64; 3],
}

impl Local {
    /// The point of the quadric's surface above the origin, along the
    /// frame's normal.
    pub(crate) fn projected(&self) -> Point {
        self.origin + self.frame[2] * self.height
    }

    /// The curvature of the quadric's surface above the origin, or `None`
    /// where it does not come out finite.
    pub(crate) fn curvature(&self) -> Option<Curvature> {
        let [t1, t2, n] = self.frame;
        let [fx, fy] = self.gradient;
        let [fxx, fxy, fyy] = self.hessian;
        let w = (1.0 + fx * fx + fy * fy).sqrt();
        // First and second fundamental forms of the graph, the second
        // taken against the normal facing away from the frame's, so that a
        // surface bending away from its normal reads positive.
        let (e, f, g) = (1.0 + fx * fx, fx * fy, 1.0 + fy * fy);
        let (l, m, nn) = (-fxx / w, -fxy / w, -fyy / w);
        let det = e * g - f * f;
        // The shape operator: the first form's inverse times the second.
        let a = (g * l - f * m) / det;
        let b = (g * m - f * nn) / det;
        let c = (e * m - f * l) / det;
        let d = (e * nn - f * m) / det;
        let mean = 0.5 * (a + d);
        let gaussian = (l * nn - m * m) / det;
        // Self-adjoint under the first form, the operator has real
        // eigenvalues; rounding can leave the discriminant just below zero.
        let spread = (mean * mean - gaussian).max(0.0).sqrt();
        let (k1, k2) = (mean + spread, mean - spread);
        // An eigenvector for `k1` in the plane's coordinates: the larger of
        // the two the operator's rows give.
        let (p, q) = ((b, k1 - a), (k1 - d, c));
        let (u, v) = if p.0.hypot(p.1) >= q.0.hypot(q.1) {
            p
        } else {
            q
        };
        let (u, v) = if u == 0.0 && v == 0.0 {
            (1.0, 0.0)
        } else {
            (u, v)
        };
        let lifted = t1 * u + t2 * v + n * (fx * u + fy * v);
        let dir1 = lifted / lifted.magnitude();
        let normal = (n - t1 * fx - t2 * fy) / w;
        let dir2 = normal.cross(dir1);
        let dir2 = dir2 / dir2.magnitude();
        let all = [
            mean, gaussian, k1, k2, dir1.x, dir1.y, dir1.z, dir2.x, dir2.y, dir2.z,
        ];
        all.iter().all(|x| x.is_finite()).then_some(Curvature {
            mean,
            gaussian,
            k1,
            k2,
            dir1,
            dir2,
        })
    }
}
