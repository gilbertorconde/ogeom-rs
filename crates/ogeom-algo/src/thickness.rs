//! Wall thickness at many points of a mesh or a shape, by rays or by
//! inscribed spheres.
//!
//! [`thickness`] measures against a triangle mesh through its [`Bvh`], and
//! [`thickness_of_shape`] against a B-rep's own surfaces, never a
//! tessellation of them. Each point carries its outward normal, and the
//! wall is measured into the material from the point, against the normal.
//! Both answer every point independently, in parallel, in the points'
//! order, and the same points give the same answers bit for bit on any
//! thread count.
//!
//! # Rays
//!
//! [`Method::Ray`] casts from the point along its inward normal and reads
//! the distance to the first surface the ray crosses, on either side of it:
//! a triangle's winding, or a face's orientation, does not matter. With a
//! cone, `rays` rays are cast, the first along the normal and the others
//! spread over the cone by a golden-angle spiral that covers its solid
//! angle evenly, and the shortest reading wins.
//!
//! # Spheres
//!
//! [`Method::Sphere`] reads the diameter of the largest ball that touches
//! the point, has its centre on the inward normal, and holds no surface: the
//! medial ball through the point, whose diameter is the wall measured
//! square across it. Every such ball is tangent to the point's tangent
//! plane, so the smaller ones nest in the larger, and the search starts
//! from the ball centred halfway along the inward ray to its first
//! crossing, which the ray's far end bounds: the ball is never longer than
//! the ray. A ball tangent at `p` along the inward unit normal `u` passes
//! through `x` at radius `|x − p|² / (2 (x − p)·u)`, the circumradius of
//! the shrinking-ball iteration.
//!
//! On a mesh the medial ball is found in one pass, exactly: a triangle
//! outside the starting ball is outside every ball nested in it, so the
//! medial radius is the least, over the triangles the starting ball holds,
//! of the radius at which a ball growing from `p` first touches each one
//! (inside its face where the ball is tangent to its plane, on an edge
//! where it is tangent to the edge's line, or at a corner, whichever is
//! least). The ball is the medial ball of the mesh as faceted: where it
//! would touch a convex curved wall along a whole circle, as a rod's
//! across its diameter, the facets turned in round the point it touches
//! cut it, and it reads less than the smooth wall would.
//!
//! On a shape the iteration is the shrinking ball's: the surface nearest
//! the centre is found by projection, and while it lies inside the ball by
//! more than the confusion distance the ball shrinks to the one through
//! it. A step from a point `δ` inside a ball of radius `r` shrinks it by
//! `δ (2r − δ) / (2 (x − p)·u)`, at least `δ / 4`, so the steps end, with a
//! ball that holds no surface deeper inside it than the confusion distance.
//!
//! A wedge between two planes at an angle `α` therefore reads `2 s tan(α/2)`
//! at a distance `s` from its edge along one face, where a ray along the
//! normal reads `s tan α`.
//!
//! # What is measured to
//!
//! A surface within [`Thickness::facing_limit`] of edge-on to the ray (for
//! a sphere, to the point's normal) is passed over: the side walls of a
//! plate do not cut off the plate's thickness near its rim. On a mesh, the
//! triangles within the confusion distance of the point are its own and
//! are passed over too: a ball that touches a corner of a faceted surface
//! would otherwise be cut by the facets meeting there. On a shape the faces
//! the point lies on are passed over at the ray's start only, so the far
//! side of a curved face that folds back is still met.

mod shape;

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail, parallel};
use ogeom_math::{Point, Vector};
use ogeom_mesh::Bvh;

pub use shape::thickness_of_shape;

/// How [`thickness`] and [`thickness_of_shape`] read a wall.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Method {
    /// Along rays from the point into the material: the shortest of `rays`
    /// rays, the first along the inward normal and the rest spread within
    /// `cone` radians of it, in `[0, π/2)`. One ray, or no cone, casts
    /// along the normal alone.
    Ray {
        /// The half-angle of the cone the rays are spread over.
        cone: f64,
        /// How many rays to cast, at least one.
        rays: u32,
    },
    /// The diameter of the largest ball touching the point with its centre
    /// on the inward normal and no surface inside.
    Sphere,
}

/// How [`thickness`] and [`thickness_of_shape`] measure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Thickness {
    /// Rays or spheres.
    pub method: Method,
    /// The angle, in radians in `[0, π/2]`, from edge-on to the ray within
    /// which a surface is passed over. Zero passes over only a surface
    /// exactly edge-on.
    pub facing_limit: f64,
    /// The thickest a wall may read; a point whose wall is thicker answers
    /// `None`. May be infinite.
    pub max: f64,
    /// The tolerances: on a mesh, the triangles within their confusion
    /// distance of a point are its own; on a shape, the precision faces
    /// are met and projected to.
    pub tol: Tolerances,
}

impl Default for Thickness {
    /// Along the normal alone, passing over surfaces exactly edge-on, at
    /// any thickness, in millimetres.
    fn default() -> Self {
        Self {
            method: Method::Ray { cone: 0.0, rays: 1 },
            facing_limit: 0.0,
            max: f64::INFINITY,
            tol: Tolerances::millimetres(),
        }
    }
}

impl Thickness {
    /// The options checked: `max` not negative or undefined, the facing
    /// limit in `[0, π/2]`, and for rays at least one ray and a cone in
    /// `[0, π/2)`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for any of
    /// them out of its range.
    pub fn check(&self) -> OgeomResult<()> {
        if self.max.is_nan() || self.max < 0.0 {
            ogeom_bail!(
                Domain,
                "the thickest wall {} must not be negative",
                self.max
            );
        }
        if !(0.0..=core::f64::consts::FRAC_PI_2).contains(&self.facing_limit) {
            ogeom_bail!(
                Domain,
                "the facing limit {} must be in [0, π/2]",
                self.facing_limit
            );
        }
        if let Method::Ray { cone, rays } = self.method {
            if rays == 0 {
                ogeom_bail!(Domain, "at least one ray must be cast");
            }
            if !(0.0..core::f64::consts::FRAC_PI_2).contains(&cone) {
                ogeom_bail!(Domain, "the cone {cone} must be in [0, π/2)");
            }
        }
        Ok(())
    }

    /// The options checked against the points and their normals, and the
    /// sine of the facing limit: a surface is met where its unit normal's
    /// dot with the unit ray is larger than it in magnitude.
    fn checked(&self, points: usize, normals: usize) -> OgeomResult<f64> {
        self.check()?;
        if normals != points {
            ogeom_bail!(Dimension, "{normals} normals for {points} points");
        }
        Ok(self.facing_limit.sin())
    }
}

/// The wall thickness at each point of the mesh `target` was built over,
/// measured into the material against the point's outward normal; in the
/// points' order, `None` where it is thicker than [`Thickness::max`] or the
/// inward ray leaves the mesh without crossing it.
///
/// See the [module](self) for what each method reads and what it measures
/// to. The points are meant to lie on the mesh (a scan's vertices, with
/// their normals): a point off it by more than the confusion distance has
/// no triangles of its own, and reads the surface it stands off as its
/// wall. A point that is not finite, or whose normal has no direction,
/// answers `None`.
///
/// # Errors
///
/// As [`Thickness::check`];
/// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) if there is
/// not one normal per point.
pub fn thickness(
    target: &Bvh,
    points: &[Point],
    normals: &[Vector],
    options: &Thickness,
) -> OgeomResult<Vec<Option<f64>>> {
    let sin = options.checked(points.len(), normals.len())?;
    let faces = face_normals(target);
    let mesh = Mesh {
        target,
        faces: &faces,
        sin,
        options,
    };
    Ok(parallel::map_ordered(points, |i, &p| {
        let u = inward(p, normals[i])?;
        mesh.at(p, u)
    }))
}

/// The unit inward normal for a point and its outward normal, or `None`
/// where the point is not finite or the normal has no direction.
fn inward(p: Point, normal: Vector) -> Option<Vector> {
    if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
        return None;
    }
    let length = normal.magnitude();
    (length > 0.0 && length.is_finite()).then(|| -normal / length)
}

/// Each triangle's unit normal, zero for a degenerate one.
fn face_normals(target: &Bvh) -> Vec<Vector> {
    (0..target.triangle_count())
        .map(|t| {
            u32::try_from(t)
                .ok()
                .and_then(|t| target.triangle(t))
                .map_or(Vector::ZERO, |[a, b, c]| {
                    let n = (b - a).cross(c - a);
                    let length = n.magnitude();
                    if length > 0.0 && length.is_finite() {
                        n / length
                    } else {
                        Vector::ZERO
                    }
                })
        })
        .collect()
}

/// Whether a surface with unit normal `m` is met by a ray along the unit
/// `v`: turned further than the facing limit, whose sine is `sin`, from
/// edge-on. A surface with no normal is met by nothing.
fn meets(m: Vector, v: Vector, sin: f64) -> bool {
    m.dot(v).abs() > sin
}

/// The radius at which a ball tangent at `p` along the unit `u` passes
/// through `x`, or `None` where no ball tangent there does.
fn through(p: Point, u: Vector, x: Point) -> Option<f64> {
    let v = x - p;
    let height = v.dot(u);
    (height > 0.0).then(|| v.dot(v) / (2.0 * height))
}

/// The `k`th of `rays` directions spread within `cone` of the unit `axis`:
/// the axis first, then a golden-angle spiral whose rings are equally
/// spaced in the cosine of their angle from the axis, out to the cone's
/// rim, so each ray stands for an equal share of the cone's solid angle.
fn spread(axis: Vector, cone: f64, rays: u32, k: u32) -> Vector {
    if k == 0 || rays < 2 {
        return axis;
    }
    let (b1, b2) = basis(axis);
    let share = f64::from(k) / f64::from(rays - 1);
    let cos = 1.0 - share * (1.0 - cone.cos());
    let sin = (1.0 - cos * cos).max(0.0).sqrt();
    let turn = f64::from(k) * core::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    let v = axis * cos + (b1 * turn.cos() + b2 * turn.sin()) * sin;
    v / v.magnitude()
}

/// Two unit vectors completing the unit `n` to a right-handed orthonormal
/// frame, continuous everywhere but where `n.z` changes sign.
fn basis(n: Vector) -> (Vector, Vector) {
    let sign = 1.0_f64.copysign(n.z);
    let a = -1.0 / (sign + n.z);
    let b = n.x * n.y * a;
    (
        Vector::new(1.0 + sign * n.x * n.x * a, sign * b, -sign * n.x),
        Vector::new(b, sign + n.y * n.y * a, -n.y),
    )
}

/// A mesh measured against, with what every point asks of it.
struct Mesh<'a> {
    target: &'a Bvh,
    faces: &'a [Vector],
    sin: f64,
    options: &'a Thickness,
}

impl Mesh<'_> {
    /// The wall at `p`, inward along the unit `u`.
    fn at(&self, p: Point, u: Vector) -> Option<f64> {
        let own = self.target.within_ball(p, self.options.tol.confusion());
        let max = self.options.max;
        match self.options.method {
            Method::Ray { cone, rays } => (0..rays)
                .filter_map(|k| self.ray(p, spread(u, cone, rays, k), max, &own))
                .min_by(f64::total_cmp),
            Method::Sphere => {
                let ray = self.ray(p, u, max, &own);
                self.sphere(p, u, ray, &own)
            }
        }
    }

    /// The distance along the unit `v` from `p` to the first triangle met,
    /// not `p`'s own, within `max`.
    fn ray(&self, p: Point, v: Vector, max: f64, own: &[u32]) -> Option<f64> {
        self.target
            .hits_along(p, v, max, true)
            .into_iter()
            .find(|hit| {
                own.binary_search(&hit.triangle).is_err()
                    && meets(self.faces[hit.triangle as usize], v, self.sin)
            })
            .map(|hit| hit.distance)
    }

    /// The diameter of the medial ball at `p` along the unit `u`, shrunk
    /// from half the inward ray's reading `ray`, or from half the thickest
    /// wall where the ray met nothing within it.
    fn sphere(&self, p: Point, u: Vector, ray: Option<f64>, own: &[u32]) -> Option<f64> {
        let mut radius = match ray {
            Some(length) => length / 2.0,
            None if self.options.max.is_finite() => self.options.max / 2.0,
            None => return None,
        };
        // Every ball smaller than this one nests in it, so a triangle
        // outside it is outside them all: the least touching radius over
        // the triangles it holds is the medial ball's.
        let mut bounded = ray.is_some();
        let start = radius;
        for t in self.target.within_ball(p + u * start, start) {
            if own.binary_search(&t).is_ok() || !meets(self.faces[t as usize], u, self.sin) {
                continue;
            }
            bounded = true;
            if let Some(r) = touch(p, u, self.target.triangle(t)?)
                && r < radius
            {
                radius = r;
            }
        }
        bounded.then_some(2.0 * radius)
    }
}

/// The radius at which a ball tangent at `p` along the unit `u` first
/// touches the triangle `t`, growing from nothing, or `None` where none
/// does: the least of the radii through each point of it, found inside
/// the face where the ball is tangent to its plane, on an edge where it is
/// tangent to the edge's line, or at a corner.
fn touch(p: Point, u: Vector, t: [Point; 3]) -> Option<f64> {
    let mut least: Option<f64> = None;
    let mut offer = |r: f64| {
        if r.is_finite() && least.is_none_or(|held| r < held) {
            least = Some(r);
        }
    };
    for &x in &t {
        if let Some(r) = through(p, u, x) {
            offer(r);
        }
    }
    for k in 0..3 {
        let (a, b) = (t[k], t[(k + 1) % 3]);
        for s in edge_critical(p, u, a, b) {
            if let Some(r) = through(p, u, a + (b - a) * s) {
                offer(r);
            }
        }
    }
    let normal = (t[1] - t[0]).cross(t[2] - t[0]);
    let length = normal.magnitude();
    if length > 0.0 && length.is_finite() {
        let mut m = normal / length;
        let mut offset = (t[0] - p).dot(m);
        if offset < 0.0 {
            m = -m;
            offset = -offset;
        }
        let lean = 1.0 + u.dot(m);
        if offset > 0.0 && lean > 0.0 {
            let r = offset / lean;
            // The ball's centre is `r` short of the plane along `m`, so it
            // touches the plane at the centre moved `r` along `m`.
            let x = p + u * r + m * r;
            if inside(t, normal, x) {
                offer(r);
            }
        }
    }
    least
}

/// The parameters in `(0, 1)` along the edge from `a` to `b` where the
/// radius of the ball tangent at `p` along `u` through the edge's point is
/// stationary: the roots of `Aβ s² + 2Aα s + (2Bα − Cβ)`, with `f = a − p`,
/// `e = b − a`, `A = e·e`, `B = f·e`, `C = f·f`, `α = f·u` and `β = e·u`,
/// the derivative's numerator.
fn edge_critical(p: Point, u: Vector, a: Point, b: Point) -> impl Iterator<Item = f64> {
    let (f, e) = (a - p, b - a);
    let (aa, bb, cc) = (e.dot(e), f.dot(e), f.dot(f));
    let (alpha, beta) = (f.dot(u), e.dot(u));
    let quadratic = aa * beta;
    let linear = 2.0 * aa * alpha;
    let constant = 2.0 * bb * alpha - cc * beta;
    let mut roots = [f64::NAN; 2];
    if quadratic == 0.0 {
        if linear != 0.0 {
            roots[0] = -constant / linear;
        }
    } else {
        let discriminant = linear * linear - 4.0 * quadratic * constant;
        if discriminant >= 0.0 {
            let q = -0.5 * (linear + discriminant.sqrt().copysign(linear));
            roots[0] = q / quadratic;
            if q != 0.0 {
                roots[1] = constant / q;
            }
        }
    }
    roots.into_iter().filter(|s| *s > 0.0 && *s < 1.0)
}

/// Whether `x`, in the triangle's plane, lies inside the triangle or on its
/// edges, against the triangle's own normal `normal`.
fn inside(t: [Point; 3], normal: Vector, x: Point) -> bool {
    (0..3).all(|k| {
        let (a, b) = (t[k], t[(k + 1) % 3]);
        (b - a).cross(x - a).dot(normal) >= 0.0
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    /// A ball tangent at the origin along `+z` touches a triangle far below
    /// a tilted plane where the closed form for the plane says, a corner
    /// where only the corner reaches, and an edge at the edge's own
    /// tangency, each checked against a dense scan of the triangle.
    #[test]
    fn the_touching_radius_is_the_least_over_the_triangle() {
        let p = Point::ORIGIN;
        let u = Vector::Z;
        let cases = [
            [
                Point::new(-5.0, -5.0, 2.0),
                Point::new(5.0, -5.0, 3.0),
                Point::new(0.0, 6.0, 2.5),
            ],
            [
                Point::new(3.0, 0.0, 1.0),
                Point::new(6.0, 1.0, 1.5),
                Point::new(5.0, 4.0, 2.0),
            ],
            [
                Point::new(-4.0, 2.0, 1.0),
                Point::new(4.0, 2.0, 1.0),
                Point::new(0.0, 6.0, 9.0),
            ],
        ];
        for t in cases {
            let found = touch(p, u, t).unwrap();
            let mut least = f64::INFINITY;
            let n = 600;
            for i in 0..=n {
                for j in 0..=(n - i) {
                    let (a, b) = (f64::from(i) / f64::from(n), f64::from(j) / f64::from(n));
                    let x = t[0] + (t[1] - t[0]) * a + (t[2] - t[0]) * b;
                    if let Some(r) = through(p, u, x) {
                        least = least.min(r);
                    }
                }
            }
            assert!(found <= least + 1e-12, "{found} above the scan's {least}");
            assert!(
                found >= least - 1e-4 * least,
                "{found} below the scan's {least}"
            );
        }
    }

    /// The spread rays start on the axis, stay within the cone, and reach
    /// its rim.
    #[test]
    fn spread_rays_cover_the_cone_out_to_its_rim() {
        let axis = Vector::new(0.3, -0.5, -0.8);
        let axis = axis / axis.magnitude();
        let cone = 0.4;
        let rays = 50;
        let angles: Vec<f64> = (0..rays)
            .map(|k| {
                spread(axis, cone, rays, k)
                    .dot(axis)
                    .clamp(-1.0, 1.0)
                    .acos()
            })
            .collect();
        assert!(angles[0].abs() < 1e-12);
        assert!(angles.iter().all(|&a| a <= cone + 1e-12));
        assert!((angles[rays as usize - 1] - cone).abs() < 1e-9);
    }
}
