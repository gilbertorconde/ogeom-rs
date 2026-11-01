//! Signed distances from many points to a mesh or a shape.
//!
//! [`to_mesh`] measures to a triangle mesh through its [`Bvh`], and
//! [`to_shape`] to a B-rep's own surfaces, curves and points, never to a
//! tessellation of them. Both answer every point independently, in
//! parallel, in the points' order, and the same points give the same
//! answers bit for bit on any thread count.
//!
//! # The sign
//!
//! A distance is positive outside and negative inside. On a mesh the side
//! is read from the closest feature's angle-weighted pseudonormal: the
//! triangle's normal where the closest point is inside a triangle, the sum
//! of the two triangles' normals where it is on an edge, and the sum of the
//! triangles' normals at a vertex weighted by their corner angles there.
//! For a closed, consistently wound mesh the displacement from the closest
//! point to the query lies in that feature's normal cone, and its dot with
//! the pseudonormal has the sign of the side, so a point beside an edge or
//! a corner is never read on the wrong side. Corners that share a position
//! are one vertex, so a mesh written with its vertices repeated per face
//! reads as the surface it draws.
//!
//! On a shape the side is the face's own: the oriented normal at the foot
//! where the closest point is inside a face, the sum of the faces' normals
//! where it is on an edge (the displacement then lies in the plane across
//! the edge, between the two normals outside or their opposites inside),
//! and the solid's own classification of the point where it is at a vertex.
//! A shape with no closed boundary, or a point the classification puts on
//! it, reads a vertex by the sum of its faces' normals.
//!
//! # Facing
//!
//! With [`DistanceOptions::facing_limit`], each point carries a normal of
//! its own (a scan's) and only surface turned within the limit of it is
//! measured to: on a mesh, the triangles whose normals are; on a shape,
//! the face feet whose normals are, and the edges and vertices where one of
//! their faces' normals is. The nearest of what is left is the answer, so
//! a point scanned on one side of a thin wall reads that side however near
//! the other is. The test is on cosines, `n·m ≥ cos(limit)` for unit
//! normals, so a surface exactly at the limit is in when the rounding of
//! the cosine allows. The sign is still the side of the feature the
//! answer lies on, read as above.

mod mesh;
mod shape;

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::{Point, Vector};
use ogeom_topo::Shape;

pub use mesh::to_mesh;
pub(crate) use mesh::{Pseudonormals, to_mesh_with};
pub use shape::to_shape;

#[cfg(doc)]
use ogeom_mesh::Bvh;

/// How [`to_mesh`] and [`to_shape`] measure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistanceOptions {
    /// The farthest a closest point may be; a point with nothing within it
    /// answers `None`. May be infinite.
    pub max: f64,
    /// The greatest angle, in radians, between a point's own normal and
    /// the surface it is measured to, in `[0, π]`; `None` measures to all
    /// of it.
    pub facing_limit: Option<f64>,
    /// Whether a point inside reads negative. Unsigned distances are never
    /// negative.
    pub signed: bool,
}

impl Default for DistanceOptions {
    /// Everything, at any distance, signed.
    fn default() -> Self {
        Self {
            max: f64::INFINITY,
            facing_limit: None,
            signed: true,
        }
    }
}

impl DistanceOptions {
    /// The options checked: `max` not negative or undefined, and the facing
    /// limit, where there is one, in `[0, π]`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for either out
    /// of its range.
    pub fn check(&self) -> OgeomResult<()> {
        if self.max.is_nan() || self.max < 0.0 {
            ogeom_bail!(
                Domain,
                "the farthest distance {} must not be negative",
                self.max
            );
        }
        if let Some(limit) = self.facing_limit
            && !(0.0..=core::f64::consts::PI).contains(&limit)
        {
            ogeom_bail!(Domain, "the facing limit {limit} must be in [0, π]");
        }
        Ok(())
    }

    /// The options checked against the points and their normals: normals
    /// where there is a facing limit, as many of them as points, and the
    /// cosine of the limit where there is one.
    fn checked(&self, points: usize, normals: Option<&[Vector]>) -> OgeomResult<Option<f64>> {
        self.check()?;
        if let Some(normals) = normals
            && normals.len() != points
        {
            ogeom_bail!(Dimension, "{} normals for {points} points", normals.len());
        }
        let Some(limit) = self.facing_limit else {
            return Ok(None);
        };
        if normals.is_none() {
            ogeom_bail!(Domain, "a facing limit needs the points' normals");
        }
        Ok(Some(limit.cos()))
    }
}

/// What a point is nearest to.
#[derive(Debug, Clone, PartialEq)]
pub struct Closest<F> {
    /// The distance to [`Closest::point`]: negative inside where the
    /// distance is signed.
    pub distance: f64,
    /// The closest point.
    pub point: Point,
    /// The face, edge or vertex it lies on.
    pub feature: F,
    /// The unit normal the sign is read from: the face's, or at an edge or
    /// a vertex the pseudonormal there. Zero where the faces' normals
    /// cancel.
    pub normal: Vector,
}

/// A feature of a triangle mesh, by a triangle that has it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MeshFeature {
    /// Inside the triangle.
    Face(u32),
    /// On the triangle's edge `k`, from its corner `k` to its corner
    /// `k + 1` (mod 3).
    Edge(u32, u8),
    /// At the triangle's corner `k`.
    Vertex(u32, u8),
}

/// A feature of a shape: the face, edge or vertex the closest point lies
/// on, as the shape holds it.
#[derive(Debug, Clone, PartialEq)]
pub enum ShapeFeature {
    /// Inside a face.
    Face(Shape),
    /// On an edge.
    Edge(Shape),
    /// At a vertex.
    Vertex(Shape),
}

/// A point's own normal made unit, or `None` where it has no direction.
fn unit(v: Vector) -> Option<Vector> {
    let length = v.magnitude();
    (length > 0.0 && length.is_finite()).then(|| v / length)
}

/// `v` made unit, or zero where it has no direction.
fn unit_or_zero(v: Vector) -> Vector {
    unit(v).unwrap_or(Vector::ZERO)
}

/// The distance signed by the side `normal` puts `p` on, seen from the
/// closest point `at`: negative where the displacement runs against it.
fn signed(distance: f64, p: Point, at: Point, normal: Vector, sign: bool) -> f64 {
    if sign && (p - at).dot(normal) < 0.0 {
        -distance
    } else {
        distance
    }
}
