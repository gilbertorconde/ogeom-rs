//! Rigid and similarity transforms, and general affine transforms.
//!
//! [`Transform`] is a similarity: an orthonormal linear part, a uniform scale
//! and a translation. That covers everything a solid modeller applies to a
//! shape — placement, rotation, mirroring, uniform scaling — while preserving
//! the two properties the geometry depends on: angles are unchanged, and an
//! analytic surface stays the same *kind* of analytic surface. A cylinder
//! remains a cylinder.
//!
//! [`GeneralTransform`] drops both guarantees, allowing non-uniform scaling and
//! shear. It is a separate type on purpose: applying one turns a circle into an
//! ellipse and a cylinder into something with no analytic form at all, so it
//! cannot be used interchangeably.
//!
//! # Form classification
//!
//! Every [`Transform`] carries a [`TransformKind`], and applying one dispatches
//! on it: a translation adds a vector, the identity does nothing at all. That
//! matters because transforms are applied to every control point of every
//! curve, every vertex of every tessellation, over an entire model — the
//! difference between a branch and nine multiplies, repeated a hundred million
//! times, is real.
//!
//! The kind is *derived from* the data rather than asserted alongside it, so it
//! cannot drift out of agreement with the matrix it describes. Every
//! constructor routes through one private classifier; there is no way to build a
//! transform that claims more structure than it has.

use core::ops::Mul;

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};

use crate::{
    Axis, Direction, Direction2, Frame, Frame2, Matrix2, Matrix3, Point, Point2, Quaternion,
    Vector, Vector2,
};

/// How much structure a [`Transform`] has, for dispatch.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TransformKind {
    /// Does nothing.
    #[default]
    Identity,
    /// Translation only.
    Translation,
    /// Rotation about an axis through the origin, possibly with a translation.
    Rotation,
    /// Reflection through a point — equivalently, a scale of `-1`.
    PointMirror,
    /// Reflection in a plane.
    PlaneMirror,
    /// Uniform scaling, possibly with a translation.
    Scale,
    /// Anything else: a combination with no simpler description.
    Compound,
}

/// A similarity transform: orthonormal rotation or reflection, uniform scale,
/// translation.
///
/// Applied as `p -> linear * (scale * p) + translation`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Transform {
    linear: Matrix3,
    scale: f64,
    translation: Vector,
    kind: TransformKind,
}

/// A similarity transform in the plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Transform2 {
    linear: Matrix2,
    scale: f64,
    translation: Vector2,
    kind: TransformKind,
}

/// A general affine transform: any linear part, plus a translation.
///
/// Non-uniform scaling and shear are allowed, so angles are not preserved and
/// analytic geometry does not survive intact. Kept distinct from [`Transform`]
/// so that cannot happen by accident.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeneralTransform {
    /// The linear part.
    pub linear: Matrix3,
    /// The translation.
    pub translation: Vector,
}

/// How far from `1` a scale factor may sit and still count as unit, and how far
/// a matrix may stray from a canonical form and still be classified as it.
/// Dimensionless; classification is a fast-path hint, and a value that just
/// misses the threshold is merely applied by the general path.
const CLASSIFY_EPS: f64 = 1e-12;

impl Default for Transform {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Transform {
    /// The identity.
    pub const IDENTITY: Self = Self {
        linear: Matrix3::IDENTITY,
        scale: 1.0,
        translation: Vector::ZERO,
        kind: TransformKind::Identity,
    };

    /// Derive the kind from the parts, and build the transform.
    ///
    /// The single constructor everything else routes through, so a transform's
    /// kind can never disagree with what it actually does.
    fn build(linear: Matrix3, scale: f64, translation: Vector) -> Self {
        let kind = Self::classify(&linear, scale, translation);
        Self {
            linear,
            scale,
            translation,
            kind,
        }
    }

    /// A transform from the parts it is stored as.
    ///
    /// `linear` is the orthonormal part alone and `scale` the uniform factor
    /// beside it, which is how a [`Transform`] holds them. The kind is
    /// re-derived rather than taken on trust, since it is a function of the
    /// other three.
    ///
    /// For reading a document back. Going the long way round — multiplying the
    /// scale into the matrix and asking
    /// [`GeneralTransform::to_similarity`](crate::GeneralTransform::to_similarity)
    /// to factor it out again — recovers a transform that is *close*, not the
    /// one that was written, and a round trip that drifts a little each time is
    /// not a round trip.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `linear` is
    /// not orthonormal within `eps`, or `scale` is not finite and non-zero — a
    /// placement that squashes space is not a placement.
    pub fn from_parts(
        linear: Matrix3,
        scale: f64,
        translation: Vector,
        eps: f64,
    ) -> OgeomResult<Self> {
        if !scale.is_finite() || scale == 0.0 {
            ogeom_bail!(
                Construction,
                "a placement's scale must be finite and non-zero; got {scale}"
            );
        }
        if !linear.is_orthonormal(eps) {
            ogeom_bail!(
                Construction,
                "a placement's linear part must be orthonormal; this one shears                  or scales unevenly"
            );
        }
        Ok(Self::build(linear, scale, translation))
    }

    /// Classify a similarity from its parts.
    fn classify(linear: &Matrix3, scale: f64, translation: Vector) -> TransformKind {
        let is_identity_linear = linear.is_equal(&Matrix3::IDENTITY, CLASSIFY_EPS);
        let unit_scale = (scale - 1.0).abs() <= CLASSIFY_EPS;
        let negative_unit_scale = (scale + 1.0).abs() <= CLASSIFY_EPS;
        let no_translation = translation.square_magnitude() == 0.0;

        if is_identity_linear && unit_scale {
            return if no_translation {
                TransformKind::Identity
            } else {
                TransformKind::Translation
            };
        }
        if is_identity_linear && negative_unit_scale {
            return TransformKind::PointMirror;
        }
        if is_identity_linear {
            return TransformKind::Scale;
        }
        if unit_scale && linear.is_orthonormal(CLASSIFY_EPS) {
            // Determinant separates a rotation from a reflection; both are
            // orthonormal, and confusing them flips the sense of every face.
            return if linear.determinant() > 0.0 {
                TransformKind::Rotation
            } else {
                TransformKind::PlaneMirror
            };
        }
        TransformKind::Compound
    }

    /// A translation.
    #[must_use]
    pub fn translation(v: Vector) -> Self {
        Self::build(Matrix3::IDENTITY, 1.0, v)
    }

    /// A rotation of `angle` radians about `axis`.
    #[must_use]
    pub fn rotation(axis: Axis, angle: f64) -> Self {
        let linear = Matrix3::rotation(axis.direction, angle);
        // Rotating about an axis that misses the origin: move the axis point to
        // the origin, rotate, move back. Written as one translation so the
        // result stays a single transform.
        let p = axis.location.to_vector();
        Self::build(linear, 1.0, p - linear * p)
    }

    /// A rotation given as a quaternion, about an axis through the origin.
    #[must_use]
    pub fn from_quaternion(q: Quaternion) -> Self {
        Self::build(q.to_matrix(), 1.0, Vector::ZERO)
    }

    /// A uniform scaling about `centre`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `factor` is
    /// zero or non-finite. A zero scale collapses every shape to a point and is
    /// not invertible, so it is refused rather than allowed to produce
    /// degenerate geometry.
    pub fn scaling(centre: Point, factor: f64, tol: Tolerances) -> OgeomResult<Self> {
        if !factor.is_finite() || factor.abs() <= tol.confusion() {
            ogeom_bail!(Construction, "scale factor {factor} is degenerate");
        }
        let c = centre.to_vector();
        Ok(Self::build(Matrix3::IDENTITY, factor, c - c * factor))
    }

    /// Reflection through a point.
    #[must_use]
    pub fn point_mirror(centre: Point) -> Self {
        let c = centre.to_vector();
        Self::build(Matrix3::IDENTITY, -1.0, c + c)
    }

    /// Reflection in the plane through `origin` with the given `normal`.
    #[must_use]
    pub fn plane_mirror(origin: Point, normal: Direction) -> Self {
        let linear = Matrix3::reflection(normal);
        let p = origin.to_vector();
        Self::build(linear, 1.0, p - linear * p)
    }

    /// Reflection in a line — a half turn about it.
    #[must_use]
    pub fn axis_mirror(axis: Axis) -> Self {
        Self::rotation(axis, core::f64::consts::PI)
    }

    /// The transform taking world coordinates into `frame`'s local coordinates.
    #[must_use]
    pub fn to_frame(frame: &Frame) -> Self {
        let linear = frame.to_matrix().transposed();
        Self::build(linear, 1.0, -(linear * frame.origin().to_vector()))
    }

    /// The transform taking `frame`'s local coordinates into world coordinates.
    #[must_use]
    pub fn from_frame(frame: &Frame) -> Self {
        Self::build(frame.to_matrix(), 1.0, frame.origin().to_vector())
    }

    /// The transform taking `from`'s local coordinates into `to`'s.
    #[must_use]
    pub fn between_frames(from: &Frame, to: &Frame) -> Self {
        Self::to_frame(to) * Self::from_frame(from)
    }

    /// This transform's classification.
    #[must_use]
    pub const fn kind(&self) -> TransformKind {
        self.kind
    }

    /// The orthonormal part.
    #[must_use]
    pub const fn linear(&self) -> Matrix3 {
        self.linear
    }

    /// The uniform scale factor. Negative for a point mirror.
    #[must_use]
    pub const fn scale_factor(&self) -> f64 {
        self.scale
    }

    /// The translation.
    #[must_use]
    pub const fn translation_vector(&self) -> Vector {
        self.translation
    }

    /// Whether this transform preserves handedness.
    ///
    /// A shape transformed by a transform that does not must have its
    /// orientation flipped to stay consistent — otherwise a mirrored solid ends
    /// up inside out.
    #[must_use]
    pub fn preserves_handedness(&self) -> bool {
        self.linear.determinant() * self.scale.signum() > 0.0
    }

    /// Apply to a point.
    #[must_use]
    pub fn apply(&self, p: Point) -> Point {
        match self.kind {
            TransformKind::Identity => p,
            TransformKind::Translation => p + self.translation,
            TransformKind::PointMirror | TransformKind::Scale => {
                Point::from_vector(p.to_vector() * self.scale + self.translation)
            }
            TransformKind::Rotation | TransformKind::PlaneMirror => {
                Point::from_vector(self.linear * p.to_vector() + self.translation)
            }
            TransformKind::Compound => {
                Point::from_vector(self.linear * (p.to_vector() * self.scale) + self.translation)
            }
        }
    }

    /// Apply to a free vector. Translation does not affect it.
    #[must_use]
    pub fn apply_vector(&self, v: Vector) -> Vector {
        match self.kind {
            TransformKind::Identity | TransformKind::Translation => v,
            TransformKind::PointMirror | TransformKind::Scale => v * self.scale,
            TransformKind::Rotation | TransformKind::PlaneMirror => self.linear * v,
            TransformKind::Compound => self.linear * (v * self.scale),
        }
    }

    /// Apply to a direction, renormalizing.
    ///
    /// A similarity maps unit vectors to vectors of length `|scale|`, so the
    /// result is rescaled. Under a negative scale the direction reverses, which
    /// is the correct behaviour and not a sign error.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the result
    /// cannot be normalized, which a valid similarity never produces.
    pub fn apply_direction(&self, d: Direction, tol: Tolerances) -> OgeomResult<Direction> {
        match self.kind {
            TransformKind::Identity | TransformKind::Translation | TransformKind::Scale => Ok(d),
            TransformKind::PointMirror => Ok(d.reversed()),
            _ => Direction::new(self.apply_vector(d.vector()), tol),
        }
    }

    /// Apply to a frame, transforming origin and all three axes.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// transformed axes cannot be renormalized.
    pub fn apply_frame(&self, f: &Frame, tol: Tolerances) -> OgeomResult<Frame> {
        Frame::from_axes(
            self.apply(f.origin()),
            self.apply_direction(f.x(), tol)?,
            self.apply_direction(f.y(), tol)?,
            self.apply_direction(f.z(), tol)?,
            tol,
        )
    }

    /// The inverse.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if the linear part is
    /// singular, which a valid similarity never is.
    pub fn inverse(&self) -> OgeomResult<Self> {
        if self.kind == TransformKind::Identity {
            return Ok(Self::IDENTITY);
        }
        if self.kind == TransformKind::Translation {
            return Ok(Self::translation(-self.translation));
        }
        if self.scale == 0.0 {
            ogeom_bail!(Numeric, "transform has a zero scale and no inverse");
        }
        // The linear part is orthonormal, so its inverse is its transpose — no
        // need to go through a general inversion, and no rounding beyond the
        // transpose itself.
        let inv_linear = self.linear.transposed();
        let inv_scale = 1.0 / self.scale;
        Ok(Self::build(
            inv_linear,
            inv_scale,
            -(inv_linear * self.translation) * inv_scale,
        ))
    }

    /// Whether two transforms agree in effect.
    #[must_use]
    pub fn is_equal(&self, other: &Self, tol: Tolerances) -> bool {
        (self.scale - other.scale).abs() <= CLASSIFY_EPS
            && self.linear.is_equal(&other.linear, CLASSIFY_EPS)
            && self.translation.is_equal(other.translation, tol)
    }

    /// This transform as a general affine one.
    #[must_use]
    pub fn to_general(&self) -> GeneralTransform {
        GeneralTransform {
            linear: self.linear * self.scale,
            translation: self.translation,
        }
    }
}

impl Mul for Transform {
    type Output = Self;
    /// Composition. `(a * b)` applies `b` first, then `a`.
    fn mul(self, b: Self) -> Self {
        if self.kind == TransformKind::Identity {
            return b;
        }
        if b.kind == TransformKind::Identity {
            return self;
        }
        Self::build(
            self.linear * b.linear,
            self.scale * b.scale,
            self.linear * (b.translation * self.scale) + self.translation,
        )
    }
}

impl Default for Transform2 {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Transform2 {
    /// The identity.
    pub const IDENTITY: Self = Self {
        linear: Matrix2::IDENTITY,
        scale: 1.0,
        translation: Vector2::ZERO,
        kind: TransformKind::Identity,
    };

    fn build(linear: Matrix2, scale: f64, translation: Vector2) -> Self {
        let is_identity_linear = linear.is_equal(&Matrix2::IDENTITY, CLASSIFY_EPS);
        let unit_scale = (scale - 1.0).abs() <= CLASSIFY_EPS;
        let kind = if is_identity_linear && unit_scale {
            if translation.square_magnitude() == 0.0 {
                TransformKind::Identity
            } else {
                TransformKind::Translation
            }
        } else if is_identity_linear && (scale + 1.0).abs() <= CLASSIFY_EPS {
            TransformKind::PointMirror
        } else if is_identity_linear {
            TransformKind::Scale
        } else if unit_scale && (linear.determinant().abs() - 1.0).abs() <= CLASSIFY_EPS {
            if linear.determinant() > 0.0 {
                TransformKind::Rotation
            } else {
                TransformKind::PlaneMirror
            }
        } else {
            TransformKind::Compound
        };
        Self {
            linear,
            scale,
            translation,
            kind,
        }
    }

    /// A translation.
    #[must_use]
    pub fn translation(v: Vector2) -> Self {
        Self::build(Matrix2::IDENTITY, 1.0, v)
    }

    /// A rotation about `centre`.
    #[must_use]
    pub fn rotation(centre: Point2, angle: f64) -> Self {
        let linear = Matrix2::rotation(angle);
        let c = centre.to_vector();
        Self::build(linear, 1.0, c - linear * c)
    }

    /// A uniform scaling about `centre`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `factor` is
    /// degenerate.
    pub fn scaling(centre: Point2, factor: f64, tol: Tolerances) -> OgeomResult<Self> {
        if !factor.is_finite() || factor.abs() <= tol.confusion() {
            ogeom_bail!(Construction, "scale factor {factor} is degenerate");
        }
        let c = centre.to_vector();
        Ok(Self::build(Matrix2::IDENTITY, factor, c - c * factor))
    }

    /// Reflection in the line through `origin` with the given `normal`.
    #[must_use]
    pub fn line_mirror(origin: Point2, normal: Direction2) -> Self {
        let (x, y) = (normal.x(), normal.y());
        let linear = Matrix2::new([
            [(-2.0f64).mul_add(x * x, 1.0), -2.0 * x * y],
            [-2.0 * x * y, (-2.0f64).mul_add(y * y, 1.0)],
        ]);
        let p = origin.to_vector();
        Self::build(linear, 1.0, p - linear * p)
    }

    /// This transform's classification.
    #[must_use]
    pub const fn kind(&self) -> TransformKind {
        self.kind
    }

    /// The orthonormal part.
    #[must_use]
    pub const fn linear(&self) -> Matrix2 {
        self.linear
    }

    /// The uniform scale factor. Negative for a point mirror.
    #[must_use]
    pub const fn scale_factor(&self) -> f64 {
        self.scale
    }

    /// The translation.
    #[must_use]
    pub const fn translation_vector(&self) -> Vector2 {
        self.translation
    }

    /// Whether this transform preserves handedness.
    #[must_use]
    pub fn preserves_handedness(&self) -> bool {
        self.linear.determinant() * self.scale.signum() > 0.0
    }

    /// Apply to a point.
    #[must_use]
    pub fn apply(&self, p: Point2) -> Point2 {
        match self.kind {
            TransformKind::Identity => p,
            TransformKind::Translation => p + self.translation,
            TransformKind::PointMirror | TransformKind::Scale => {
                Point2::from_vector(p.to_vector() * self.scale + self.translation)
            }
            TransformKind::Rotation | TransformKind::PlaneMirror => {
                Point2::from_vector(self.linear * p.to_vector() + self.translation)
            }
            TransformKind::Compound => {
                Point2::from_vector(self.linear * (p.to_vector() * self.scale) + self.translation)
            }
        }
    }

    /// Apply to a direction, renormalizing.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the result
    /// cannot be normalized, which a valid similarity never produces.
    pub fn apply_direction(&self, d: Direction2, tol: Tolerances) -> OgeomResult<Direction2> {
        match self.kind {
            TransformKind::Identity | TransformKind::Translation | TransformKind::Scale => Ok(d),
            TransformKind::PointMirror => Ok(d.reversed()),
            _ => Direction2::new(self.apply_vector(d.vector()), tol),
        }
    }

    /// Apply to a frame, transforming the origin and both axes.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// transformed axes cannot be renormalized or are no longer perpendicular.
    pub fn apply_frame(&self, f: &Frame2, tol: Tolerances) -> OgeomResult<Frame2> {
        Frame2::from_axes(
            self.apply(f.origin()),
            self.apply_direction(f.x(), tol)?,
            self.apply_direction(f.y(), tol)?,
            tol,
        )
    }

    /// Apply to a free vector.
    #[must_use]
    pub fn apply_vector(&self, v: Vector2) -> Vector2 {
        match self.kind {
            TransformKind::Identity | TransformKind::Translation => v,
            TransformKind::PointMirror | TransformKind::Scale => v * self.scale,
            TransformKind::Rotation | TransformKind::PlaneMirror => self.linear * v,
            TransformKind::Compound => self.linear * (v * self.scale),
        }
    }

    /// The inverse.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if the transform is
    /// degenerate.
    pub fn inverse(&self) -> OgeomResult<Self> {
        if self.kind == TransformKind::Identity {
            return Ok(Self::IDENTITY);
        }
        if self.scale == 0.0 {
            ogeom_bail!(Numeric, "transform has a zero scale and no inverse");
        }
        let inv_linear = self.linear.transposed();
        let inv_scale = 1.0 / self.scale;
        Ok(Self::build(
            inv_linear,
            inv_scale,
            -(inv_linear * self.translation) * inv_scale,
        ))
    }

    /// Whether two transforms agree in effect.
    #[must_use]
    pub fn is_equal(&self, other: &Self, tol: Tolerances) -> bool {
        (self.scale - other.scale).abs() <= CLASSIFY_EPS
            && self.linear.is_equal(&other.linear, CLASSIFY_EPS)
            && self.translation.is_equal(other.translation, tol)
    }
}

impl Mul for Transform2 {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        if self.kind == TransformKind::Identity {
            return b;
        }
        if b.kind == TransformKind::Identity {
            return self;
        }
        Self::build(
            self.linear * b.linear,
            self.scale * b.scale,
            self.linear * (b.translation * self.scale) + self.translation,
        )
    }
}

impl Default for GeneralTransform {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl GeneralTransform {
    /// The identity.
    pub const IDENTITY: Self = Self {
        linear: Matrix3::IDENTITY,
        translation: Vector::ZERO,
    };

    /// From a linear part and a translation.
    #[must_use]
    pub const fn new(linear: Matrix3, translation: Vector) -> Self {
        Self {
            linear,
            translation,
        }
    }

    /// Non-uniform scaling about the origin.
    #[must_use]
    pub const fn scaling_xyz(x: f64, y: f64, z: f64) -> Self {
        Self::new(Matrix3::scaling_xyz(x, y, z), Vector::ZERO)
    }

    /// Apply to a point.
    #[must_use]
    pub fn apply(&self, p: Point) -> Point {
        Point::from_vector(self.linear * p.to_vector() + self.translation)
    }

    /// Apply to a free vector.
    #[must_use]
    pub fn apply_vector(&self, v: Vector) -> Vector {
        self.linear * v
    }

    /// Apply to a normal vector.
    ///
    /// Normals transform by the inverse transpose, not by the linear part
    /// itself. Using the linear part directly is only correct for a similarity;
    /// under any shear or non-uniform scale it tilts normals off the surface
    /// they belong to, which then breaks every orientation test downstream.
    ///
    /// The result is not renormalized — it is a direction, not a length.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if the linear part is
    /// singular.
    pub fn apply_normal(&self, n: Vector) -> OgeomResult<Vector> {
        Ok(self.linear.inverse()?.transposed() * n)
    }

    /// Whether this transform preserves handedness.
    #[must_use]
    pub fn preserves_handedness(&self) -> bool {
        self.linear.determinant() > 0.0
    }

    /// The factor by which volumes are multiplied. Negative if handedness flips.
    #[must_use]
    pub fn volume_ratio(&self) -> f64 {
        self.linear.determinant()
    }

    /// The inverse.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if the linear part is
    /// singular.
    pub fn inverse(&self) -> OgeomResult<Self> {
        let inv = self.linear.inverse()?;
        Ok(Self::new(inv, -(inv * self.translation)))
    }

    /// Whether this is a similarity, and so can be narrowed to a [`Transform`].
    #[must_use]
    pub fn is_similarity(&self, eps: f64) -> bool {
        self.to_similarity(eps).is_some()
    }

    /// This transform as a [`Transform`], if it is in fact a similarity.
    ///
    /// Returns `None` when the linear part contains shear or non-uniform
    /// scaling, since no similarity describes it.
    #[must_use]
    pub fn to_similarity(&self, eps: f64) -> Option<Transform> {
        // A similarity's linear part is `s * R` with `R` orthonormal, so its
        // columns are mutually orthogonal and all of length |s|.
        let det = self.linear.determinant();
        if det == 0.0 {
            return None;
        }
        let scale = det.abs().cbrt() * det.signum();
        let rotation = self.linear * (1.0 / scale);
        if !rotation.is_orthonormal(eps) {
            return None;
        }
        Some(Transform::build(rotation, scale, self.translation))
    }
}

impl Mul for GeneralTransform {
    type Output = Self;
    /// Composition. `(a * b)` applies `b` first, then `a`.
    fn mul(self, b: Self) -> Self {
        Self::new(
            self.linear * b.linear,
            self.linear * b.translation + self.translation,
        )
    }
}

impl From<Transform> for GeneralTransform {
    fn from(t: Transform) -> Self {
        t.to_general()
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

    fn sample_points() -> [Point; 4] {
        [
            Point::ORIGIN,
            Point::new(1.0, 0.0, 0.0),
            Point::new(-3.0, 7.5, 2.25),
            Point::new(1e3, -1e3, 0.5),
        ]
    }

    #[test]
    fn classification_matches_what_the_transform_does() {
        assert_eq!(Transform::IDENTITY.kind(), TransformKind::Identity);
        assert_eq!(
            Transform::translation(Vector::X).kind(),
            TransformKind::Translation
        );
        assert_eq!(
            Transform::rotation(Axis::Z, 0.5).kind(),
            TransformKind::Rotation
        );
        assert_eq!(
            Transform::point_mirror(Point::ORIGIN).kind(),
            TransformKind::PointMirror
        );
        assert_eq!(
            Transform::plane_mirror(Point::ORIGIN, Direction::Z).kind(),
            TransformKind::PlaneMirror
        );
        assert_eq!(
            Transform::scaling(Point::ORIGIN, 3.0, T).unwrap().kind(),
            TransformKind::Scale
        );
        let compound =
            Transform::rotation(Axis::Z, 0.5) * Transform::scaling(Point::ORIGIN, 3.0, T).unwrap();
        assert_eq!(compound.kind(), TransformKind::Compound);
    }

    #[test]
    fn a_zero_rotation_classifies_as_identity_not_rotation() {
        // The classification is derived from the data, so it cannot claim more
        // structure than the transform has — or less.
        assert_eq!(
            Transform::rotation(Axis::Z, 0.0).kind(),
            TransformKind::Identity
        );
        assert_eq!(
            Transform::translation(Vector::ZERO).kind(),
            TransformKind::Identity
        );
        assert_eq!(
            Transform::scaling(Point::ORIGIN, 1.0, T).unwrap().kind(),
            TransformKind::Identity
        );
    }

    #[test]
    fn every_dispatch_path_gives_the_same_answer_as_the_general_one() {
        // The whole point of classification is speed, so each fast path must
        // agree exactly with the general formula it replaces.
        let cases = [
            Transform::IDENTITY,
            Transform::translation(Vector::new(1.0, -2.0, 3.0)),
            Transform::rotation(Axis::new(Point::new(1.0, 0.0, 0.0), Direction::Z), 0.7),
            Transform::point_mirror(Point::new(2.0, 0.0, -1.0)),
            Transform::plane_mirror(Point::new(0.0, 1.0, 0.0), Direction::Y),
            Transform::scaling(Point::new(1.0, 1.0, 1.0), 2.5, T).unwrap(),
        ];
        for t in cases {
            for p in sample_points() {
                let general = Point::from_vector(
                    t.linear() * (p.to_vector() * t.scale_factor()) + t.translation_vector(),
                );
                assert!(
                    t.apply(p).is_equal(general, T),
                    "fast path for {:?} disagrees",
                    t.kind()
                );
            }
        }
    }

    #[test]
    fn rotation_about_an_off_origin_axis_leaves_the_axis_fixed() {
        let axis = Axis::new(Point::new(5.0, 3.0, 0.0), Direction::Z);
        let t = Transform::rotation(axis, 1.234);
        assert!(t.apply(axis.location).is_equal(axis.location, T));
        assert!(
            t.apply(axis.point_at(10.0))
                .is_equal(axis.point_at(10.0), T)
        );
        // A point off the axis moves, staying at the same radius.
        let p = Point::new(6.0, 3.0, 0.0);
        assert_relative_eq!(axis.distance_to(t.apply(p)), 1.0, epsilon = 1e-14);
    }

    #[test]
    fn scaling_about_a_centre_leaves_the_centre_fixed() {
        let c = Point::new(3.0, -1.0, 2.0);
        let t = Transform::scaling(c, 4.0, T).unwrap();
        assert!(t.apply(c).is_equal(c, T));
        let p = c + Vector::new(1.0, 0.0, 0.0);
        assert!(t.apply(p).is_equal(c + Vector::new(4.0, 0.0, 0.0), T));
    }

    #[test]
    fn degenerate_scales_are_refused() {
        assert!(Transform::scaling(Point::ORIGIN, 0.0, T).is_err());
        assert!(Transform::scaling(Point::ORIGIN, f64::NAN, T).is_err());
        assert!(Transform::scaling(Point::ORIGIN, f64::INFINITY, T).is_err());
        assert!(
            Transform::scaling(Point::ORIGIN, -2.0, T).is_ok(),
            "negative is fine"
        );
    }

    #[test]
    fn handedness_tracks_mirroring() {
        assert!(Transform::rotation(Axis::Z, 1.0).preserves_handedness());
        assert!(Transform::translation(Vector::X).preserves_handedness());
        assert!(
            Transform::scaling(Point::ORIGIN, 3.0, T)
                .unwrap()
                .preserves_handedness()
        );
        assert!(!Transform::plane_mirror(Point::ORIGIN, Direction::Z).preserves_handedness());
        assert!(!Transform::point_mirror(Point::ORIGIN).preserves_handedness());
        // Two mirrors make a rotation.
        let twice = Transform::plane_mirror(Point::ORIGIN, Direction::Z)
            * Transform::plane_mirror(Point::ORIGIN, Direction::X);
        assert!(twice.preserves_handedness());
    }

    #[test]
    fn axis_mirror_is_a_half_turn() {
        let t = Transform::axis_mirror(Axis::Z);
        assert!(
            t.apply(Point::new(1.0, 0.0, 5.0))
                .is_equal(Point::new(-1.0, 0.0, 5.0), T)
        );
        assert!(t.preserves_handedness(), "a half turn is a rotation");
    }

    #[test]
    fn inverse_round_trips_for_every_kind() {
        let cases = [
            Transform::IDENTITY,
            Transform::translation(Vector::new(1.0, -2.0, 3.0)),
            Transform::rotation(Axis::new(Point::new(1.0, 2.0, 3.0), Direction::Y), 2.1),
            Transform::point_mirror(Point::new(1.0, 1.0, 1.0)),
            Transform::plane_mirror(Point::new(0.0, 0.0, 4.0), Direction::Z),
            Transform::scaling(Point::new(-1.0, 0.0, 0.0), 0.25, T).unwrap(),
        ];
        for t in cases {
            let inv = t.inverse().unwrap();
            for p in sample_points() {
                assert!(inv.apply(t.apply(p)).is_equal(p, T), "{:?}", t.kind());
                assert!(t.apply(inv.apply(p)).is_equal(p, T), "{:?}", t.kind());
            }
        }
    }

    #[test]
    fn composition_applies_right_to_left() {
        let a = Transform::translation(Vector::new(10.0, 0.0, 0.0));
        let b = Transform::rotation(Axis::Z, core::f64::consts::FRAC_PI_2);
        let p = Point::new(1.0, 0.0, 0.0);
        assert!((a * b).apply(p).is_equal(a.apply(b.apply(p)), T));
        assert!((b * a).apply(p).is_equal(b.apply(a.apply(p)), T));
        // And the two orders genuinely differ.
        assert!(!(a * b).is_equal(&(b * a), T));
    }

    #[test]
    fn composition_is_associative() {
        let a = Transform::rotation(Axis::X, 0.3);
        let b = Transform::scaling(Point::new(1.0, 0.0, 0.0), 2.0, T).unwrap();
        let c = Transform::translation(Vector::new(0.0, 5.0, 0.0));
        assert!(((a * b) * c).is_equal(&(a * (b * c)), T));
    }

    #[test]
    fn vectors_ignore_translation_and_directions_stay_unit() {
        let t = Transform::translation(Vector::new(100.0, 0.0, 0.0))
            * Transform::rotation(Axis::Z, 0.9);
        let v = Vector::new(1.0, 2.0, 3.0);
        assert!(
            t.apply_vector(v)
                .is_equal(Transform::rotation(Axis::Z, 0.9).apply_vector(v), T)
        );
        let d = t.apply_direction(Direction::X, T).unwrap();
        assert_relative_eq!(d.vector().magnitude(), 1.0, epsilon = 1e-15);
    }

    #[test]
    fn a_point_mirror_reverses_directions() {
        let t = Transform::point_mirror(Point::new(5.0, 5.0, 5.0));
        assert!(
            t.apply_direction(Direction::X, T)
                .unwrap()
                .is_equal(-Direction::X, T)
        );
        // A positive scale does not.
        let s = Transform::scaling(Point::ORIGIN, 3.0, T).unwrap();
        assert!(
            s.apply_direction(Direction::X, T)
                .unwrap()
                .is_equal(Direction::X, T)
        );
    }

    #[test]
    fn frame_transforms_round_trip_through_world() {
        let f = Frame::new(
            Point::new(1.0, 2.0, 3.0),
            Direction::from_coords(1.0, 1.0, 0.0, T).unwrap(),
            Direction::Z,
            T,
        )
        .unwrap();
        let to = Transform::to_frame(&f);
        let from = Transform::from_frame(&f);
        for p in sample_points() {
            assert!(from.apply(to.apply(p)).is_equal(p, T));
            // And they agree with the frame's own conversion.
            assert!(to.apply(p).is_equal(f.to_local(p), T));
            assert!(from.apply(f.to_local(p)).is_equal(p, T));
        }
    }

    #[test]
    fn between_frames_composes_correctly() {
        let a = Frame::new(Point::new(1.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
        let b = Frame::new(Point::new(0.0, 5.0, 0.0), Direction::X, Direction::Y, T).unwrap();
        let t = Transform::between_frames(&a, &b);
        // A point at local (1,2,3) in `a` must land at the same world position
        // when read back out of `b`.
        let local = Point::new(1.0, 2.0, 3.0);
        assert!(b.to_world(t.apply(local)).is_equal(a.to_world(local), T));
    }

    #[test]
    fn general_transform_normals_use_the_inverse_transpose() {
        // Non-uniform scaling: the plane z = x has normal (1, 0, -1) up to
        // scale. Scale x by 2 and the plane becomes z = x/2, whose normal is
        // (1, 0, -2) up to scale — not (2, 0, -1), which is what applying the
        // linear part directly would give.
        let g = GeneralTransform::scaling_xyz(2.0, 1.0, 1.0);
        let n = Vector::new(1.0, 0.0, -1.0);
        let transformed = g.apply_normal(n).unwrap();

        let on_plane = Vector::new(1.0, 0.0, 1.0);
        assert_relative_eq!(n.dot(on_plane), 0.0, epsilon = 1e-15);
        assert_relative_eq!(
            transformed.dot(g.apply_vector(on_plane)),
            0.0,
            epsilon = 1e-14,
            max_relative = 1e-14
        );
        // The naive answer does not stay perpendicular.
        assert!(g.apply_vector(n).dot(g.apply_vector(on_plane)).abs() > 1e-6);
    }

    #[test]
    fn general_transform_recognizes_similarities() {
        let similar: GeneralTransform = Transform::rotation(Axis::Z, 0.4).into();
        assert!(similar.is_similarity(1e-12));
        let narrowed = similar.to_similarity(1e-12).unwrap();
        assert_eq!(narrowed.kind(), TransformKind::Rotation);

        let scaled: GeneralTransform = Transform::scaling(Point::ORIGIN, 3.0, T).unwrap().into();
        assert!(scaled.is_similarity(1e-12));
        assert_relative_eq!(
            scaled.to_similarity(1e-12).unwrap().scale_factor(),
            3.0,
            epsilon = 1e-12
        );

        assert!(!GeneralTransform::scaling_xyz(1.0, 2.0, 3.0).is_similarity(1e-12));
        let shear = GeneralTransform::new(
            Matrix3::new([[1.0, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            Vector::ZERO,
        );
        assert!(!shear.is_similarity(1e-12));
    }

    #[test]
    fn general_transform_volume_ratio_and_inverse() {
        let g = GeneralTransform::scaling_xyz(2.0, 3.0, 4.0);
        assert_relative_eq!(g.volume_ratio(), 24.0);
        assert!(g.preserves_handedness());
        let inv = g.inverse().unwrap();
        for p in sample_points() {
            assert!(inv.apply(g.apply(p)).is_equal(p, T));
        }
        let flip = GeneralTransform::scaling_xyz(-1.0, 1.0, 1.0);
        assert!(!flip.preserves_handedness());
        assert!(
            GeneralTransform::scaling_xyz(0.0, 1.0, 1.0)
                .inverse()
                .is_err()
        );
    }

    #[test]
    fn transform2_behaves_like_its_3d_counterpart() {
        let r = Transform2::rotation(Point2::new(1.0, 1.0), core::f64::consts::FRAC_PI_2);
        assert_eq!(r.kind(), TransformKind::Rotation);
        assert!(
            r.apply(Point2::new(1.0, 1.0))
                .is_equal(Point2::new(1.0, 1.0), T)
        );
        assert!(
            r.apply(Point2::new(2.0, 1.0))
                .is_equal(Point2::new(1.0, 2.0), T)
        );
        assert!(
            r.inverse()
                .unwrap()
                .apply(r.apply(Point2::ORIGIN))
                .is_equal(Point2::ORIGIN, T)
        );

        let m = Transform2::line_mirror(Point2::ORIGIN, Direction2::Y);
        assert_eq!(m.kind(), TransformKind::PlaneMirror);
        assert!(!m.preserves_handedness());
        assert!(
            m.apply(Point2::new(3.0, 2.0))
                .is_equal(Point2::new(3.0, -2.0), T)
        );

        assert!(Transform2::scaling(Point2::ORIGIN, 0.0, T).is_err());
    }
}
