//! Unit vectors, with the invariant enforced by the type.
//!
//! A [`Direction`] is always unit length. Every constructor normalizes and can
//! fail; there is no way to build one from components without that check.
//!
//! This matters more than it looks. Surface normals, axis directions and
//! parameterization references are all directions, and an algorithm that
//! assumes unit length — as almost all of them do, implicitly, when they skip a
//! division — silently produces scaled results when handed a vector that is not.
//! Making the invariant unrepresentable-if-false removes the whole class.

use core::ops::{Mul, Neg};

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};

use crate::{Vector, Vector2};

/// A unit vector in space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Direction(Vector);

/// A unit vector in the plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Direction2(Vector2);

impl Direction {
    /// +X.
    pub const X: Self = Self(Vector::X);
    /// +Y.
    pub const Y: Self = Self(Vector::Y);
    /// +Z.
    pub const Z: Self = Self(Vector::Z);

    /// Normalize `v` into a direction.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `v` is
    /// non-finite or shorter than `tol.confusion()`.
    pub fn new(v: Vector, tol: Tolerances) -> OgeomResult<Self> {
        Ok(Self(v.normalized(tol)?))
    }

    /// Normalize components into a direction.
    ///
    /// # Errors
    ///
    /// As [`Direction::new`].
    pub fn from_coords(x: f64, y: f64, z: f64, tol: Tolerances) -> OgeomResult<Self> {
        Self::new(Vector::new(x, y, z), tol)
    }

    /// A direction from a vector that is *already* a unit vector.
    ///
    /// Checks rather than normalizes, and the distinction is the whole reason
    /// it exists: dividing a unit vector by its own magnitude does not give it
    /// back, it gives something a bit or two away. That is invisible until
    /// something has to reproduce a direction exactly — reading a document back
    /// from a file, above all, where the drift turns a round trip that should
    /// be the identity into one that changes the model a little every time.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `v` is
    /// non-finite, or its length differs from one by more than
    /// `tol.confusion()`.
    pub fn unit(v: Vector, tol: Tolerances) -> OgeomResult<Self> {
        if !v.is_finite() {
            ogeom_bail!(Construction, "a direction must be finite; got {v:?}");
        }
        let length = v.magnitude();
        if (length - 1.0).abs() > tol.confusion() {
            ogeom_bail!(
                Construction,
                "expected a unit vector, got one of length {length}"
            );
        }
        Ok(Self(v))
    }

    /// The underlying unit vector.
    #[must_use]
    pub const fn vector(self) -> Vector {
        self.0
    }

    /// X component.
    #[must_use]
    pub const fn x(self) -> f64 {
        self.0.x
    }

    /// Y component.
    #[must_use]
    pub const fn y(self) -> f64 {
        self.0.y
    }

    /// Z component.
    #[must_use]
    pub const fn z(self) -> f64 {
        self.0.z
    }

    /// Components as an array.
    #[must_use]
    pub const fn to_array(self) -> [f64; 3] {
        self.0.to_array()
    }

    /// Dot product with another direction — the cosine of the angle between
    /// them, in `[-1, 1]` up to rounding.
    #[must_use]
    pub fn dot(self, other: Self) -> f64 {
        self.0.dot(other.0)
    }

    /// Dot product with a free vector.
    #[must_use]
    pub fn dot_vector(self, v: Vector) -> f64 {
        self.0.dot(v)
    }

    /// Cross product, as a free vector. Its magnitude is the sine of the angle
    /// between the two directions, so it is *not* itself a direction — for
    /// nearly parallel inputs it is nearly null.
    #[must_use]
    pub fn cross_vector(self, other: Self) -> Vector {
        self.0.cross(other.0)
    }

    /// Cross product with a free vector.
    ///
    /// Its magnitude is the component of `v` perpendicular to this direction,
    /// which makes it the accurate way to get a perpendicular distance:
    /// subtracting the parallel component instead cancels catastrophically for
    /// a point far along the direction.
    #[must_use]
    pub fn cross_with(self, v: Vector) -> Vector {
        self.0.cross(v)
    }

    /// Cross product, renormalized into a direction.
    ///
    /// Collinearity is judged against the *angular* tolerance, not the linear
    /// one: for unit inputs the cross product's magnitude is the sine of the
    /// angle between them, a dimensionless quantity that a length tolerance
    /// does not describe.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the two
    /// directions are collinear.
    pub fn cross(self, other: Self, tol: Tolerances) -> OgeomResult<Self> {
        let v = self.cross_vector(other);
        let m = v.magnitude();
        if m <= tol.angular() {
            ogeom_bail!(Construction, "cross product of collinear directions");
        }
        Ok(Self(v / m))
    }

    /// The unit normal to two free vectors.
    ///
    /// The right way to build a normal from two edges of a triangle. Naively
    /// normalizing `a.cross(b)` compares its magnitude — which is twice the
    /// triangle's area, and so scales as the *square* of the size — against a
    /// length tolerance. A triangle a micron across then looks degenerate even
    /// though its normal is perfectly well determined. The test here is
    /// relative: `|a x b| > tol.angular() * |a| * |b|`, which asks the question
    /// that actually matters, whether the two vectors are collinear, and gives
    /// the same answer at every scale.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `a` and `b`
    /// are collinear, or either is null.
    pub fn from_cross(a: Vector, b: Vector, tol: Tolerances) -> OgeomResult<Self> {
        let v = a.cross(b);
        let m = v.magnitude();
        if m <= tol.angular() * a.magnitude() * b.magnitude() {
            ogeom_bail!(Construction, "cannot take a normal to collinear vectors");
        }
        Ok(Self(v / m))
    }

    /// Angle to `other`, in `[0, π]`.
    #[must_use]
    pub fn angle(self, other: Self) -> f64 {
        // atan2 rather than acos of the dot product: acos loses roughly half its
        // significant digits near 0 and π, which is where these tests matter.
        self.cross_vector(other).magnitude().atan2(self.dot(other))
    }

    /// Whether the two point the same way, within `tol.angular()`.
    #[must_use]
    pub fn is_equal(self, other: Self, tol: Tolerances) -> bool {
        self.angle(other) <= tol.angular()
    }

    /// Whether the two point opposite ways, within `tol.angular()`.
    #[must_use]
    pub fn is_opposite(self, other: Self, tol: Tolerances) -> bool {
        core::f64::consts::PI - self.angle(other) <= tol.angular()
    }

    /// Whether the two are parallel, ignoring sense.
    #[must_use]
    pub fn is_parallel(self, other: Self, tol: Tolerances) -> bool {
        self.is_equal(other, tol) || self.is_opposite(other, tol)
    }

    /// Whether the two are perpendicular, within `tol.angular()`.
    #[must_use]
    pub fn is_normal(self, other: Self, tol: Tolerances) -> bool {
        (core::f64::consts::FRAC_PI_2 - self.angle(other)).abs() <= tol.angular()
    }

    /// Some direction perpendicular to this one.
    ///
    /// Which one is unspecified but deterministic. Chosen by crossing with
    /// whichever axis this direction is least aligned with, so the cross product
    /// is never near-degenerate and the result is numerically sound for every
    /// input.
    #[must_use]
    pub fn any_perpendicular(self) -> Self {
        let [ax, ay, az] = [self.x().abs(), self.y().abs(), self.z().abs()];
        let axis = if ax <= ay && ax <= az {
            Vector::X
        } else if ay <= az {
            Vector::Y
        } else {
            Vector::Z
        };
        let v = self.0.cross(axis);
        // Guaranteed non-degenerate: `axis` is the least-aligned unit axis, so
        // the angle between them is at least acos(1/sqrt(3)) ~= 54.7 degrees.
        Self(v / v.magnitude())
    }

    /// This direction reflected through the origin.
    #[must_use]
    pub const fn reversed(self) -> Self {
        Self(Vector::new(-self.0.x, -self.0.y, -self.0.z))
    }

    /// This direction with the Z component dropped, renormalized.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if this
    /// direction is parallel to Z, leaving nothing to project.
    pub fn to_2d(self, tol: Tolerances) -> OgeomResult<Direction2> {
        Direction2::new(self.0.xy(), tol)
    }
}

impl Direction2 {
    /// +X.
    pub const X: Self = Self(Vector2::X);
    /// +Y.
    pub const Y: Self = Self(Vector2::Y);

    /// Normalize `v` into a direction.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `v` is
    /// non-finite or shorter than `tol.confusion()`.
    pub fn new(v: Vector2, tol: Tolerances) -> OgeomResult<Self> {
        Ok(Self(v.normalized(tol)?))
    }

    /// Normalize components into a direction.
    ///
    /// # Errors
    ///
    /// As [`Direction2::new`].
    pub fn from_coords(x: f64, y: f64, tol: Tolerances) -> OgeomResult<Self> {
        Self::new(Vector2::new(x, y), tol)
    }

    /// A direction from a vector that is *already* a unit vector.
    ///
    /// As [`Direction::unit`]: it checks rather than normalizes, so a direction
    /// read back from a document is the one that was written and not something
    /// a bit or two away from it.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `v` is
    /// non-finite, or its length differs from one by more than
    /// `tol.confusion()`.
    pub fn unit(v: Vector2, tol: Tolerances) -> OgeomResult<Self> {
        if !v.is_finite() {
            ogeom_bail!(Construction, "a direction must be finite; got {v:?}");
        }
        let length = v.magnitude();
        if (length - 1.0).abs() > tol.confusion() {
            ogeom_bail!(
                Construction,
                "expected a unit vector, got one of length {length}"
            );
        }
        Ok(Self(v))
    }

    /// The direction at `angle` radians counter-clockwise from +X.
    #[must_use]
    pub fn from_angle(angle: f64) -> Self {
        let (sin, cos) = angle.sin_cos();
        Self(Vector2::new(cos, sin))
    }

    /// The underlying unit vector.
    #[must_use]
    pub const fn vector(self) -> Vector2 {
        self.0
    }

    /// X component.
    #[must_use]
    pub const fn x(self) -> f64 {
        self.0.x
    }

    /// Y component.
    #[must_use]
    pub const fn y(self) -> f64 {
        self.0.y
    }

    /// Components as an array.
    #[must_use]
    pub const fn to_array(self) -> [f64; 2] {
        self.0.to_array()
    }

    /// Dot product — the cosine of the angle between the two.
    #[must_use]
    pub fn dot(self, other: Self) -> f64 {
        self.0.dot(other.0)
    }

    /// Scalar cross product — the sine of the signed angle from `self` to
    /// `other`.
    #[must_use]
    pub fn cross(self, other: Self) -> f64 {
        self.0.cross(other.0)
    }

    /// This direction rotated a quarter turn counter-clockwise. Exact.
    #[must_use]
    pub const fn perpendicular(self) -> Self {
        Self(self.0.perpendicular())
    }

    /// Angle from +X, in `(-π, π]`.
    #[must_use]
    pub fn to_angle(self) -> f64 {
        self.0.y.atan2(self.0.x)
    }

    /// Signed angle to `other`, in `(-π, π]`, positive counter-clockwise.
    #[must_use]
    pub fn angle(self, other: Self) -> f64 {
        self.cross(other).atan2(self.dot(other))
    }

    /// Whether the two point the same way, within `tol.angular()`.
    #[must_use]
    pub fn is_equal(self, other: Self, tol: Tolerances) -> bool {
        self.angle(other).abs() <= tol.angular()
    }

    /// Whether the two point opposite ways, within `tol.angular()`.
    #[must_use]
    pub fn is_opposite(self, other: Self, tol: Tolerances) -> bool {
        core::f64::consts::PI - self.angle(other).abs() <= tol.angular()
    }

    /// Whether the two are parallel, ignoring sense.
    #[must_use]
    pub fn is_parallel(self, other: Self, tol: Tolerances) -> bool {
        self.is_equal(other, tol) || self.is_opposite(other, tol)
    }

    /// Whether the two are perpendicular, within `tol.angular()`.
    #[must_use]
    pub fn is_normal(self, other: Self, tol: Tolerances) -> bool {
        (core::f64::consts::FRAC_PI_2 - self.angle(other).abs()).abs() <= tol.angular()
    }

    /// This direction reflected through the origin.
    #[must_use]
    pub const fn reversed(self) -> Self {
        Self(Vector2::new(-self.0.x, -self.0.y))
    }

    /// This direction embedded in the XY plane.
    #[must_use]
    pub const fn to_3d(self) -> Direction {
        Direction(Vector::new(self.0.x, self.0.y, 0.0))
    }
}

impl Neg for Direction {
    type Output = Self;
    fn neg(self) -> Self {
        self.reversed()
    }
}

impl Neg for Direction2 {
    type Output = Self;
    fn neg(self) -> Self {
        self.reversed()
    }
}

impl Mul<f64> for Direction {
    type Output = Vector;
    /// Scaling a direction yields a free vector: the result is no longer unit
    /// length, so it is no longer a direction.
    fn mul(self, s: f64) -> Vector {
        self.0 * s
    }
}

impl Mul<Direction> for f64 {
    type Output = Vector;
    fn mul(self, d: Direction) -> Vector {
        d.0 * self
    }
}

impl Mul<f64> for Direction2 {
    type Output = Vector2;
    fn mul(self, s: f64) -> Vector2 {
        self.0 * s
    }
}

impl Mul<Direction2> for f64 {
    type Output = Vector2;
    fn mul(self, d: Direction2) -> Vector2 {
        d.0 * self
    }
}

impl From<Direction> for Vector {
    fn from(d: Direction) -> Self {
        d.0
    }
}

impl From<Direction2> for Vector2 {
    fn from(d: Direction2) -> Self {
        d.0
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn every_construction_path_yields_unit_length() {
        let cases = [
            Direction::new(Vector::new(3.0, 4.0, 12.0), T).unwrap(),
            Direction::from_coords(-1.0, 2.0, -0.5, T).unwrap(),
            Direction::X.any_perpendicular(),
            Direction::new(Vector::new(1.0, 1.0, 1.0), T)
                .unwrap()
                .reversed(),
            Direction::X.cross(Direction::Y, T).unwrap(),
        ];
        for d in cases {
            assert_relative_eq!(d.vector().magnitude(), 1.0, epsilon = 1e-15);
        }
    }

    #[test]
    fn degenerate_input_is_refused() {
        assert!(Direction::new(Vector::ZERO, T).is_err());
        assert!(Direction::from_coords(f64::NAN, 0.0, 0.0, T).is_err());
        assert!(Direction2::new(Vector2::ZERO, T).is_err());
        // Collinear directions have a null cross product.
        assert!(Direction::X.cross(Direction::X, T).is_err());
        assert!(Direction::X.cross(-Direction::X, T).is_err());
        assert!(Direction::from_cross(Vector::X, Vector::X * 3.0, T).is_err());
        assert!(Direction::from_cross(Vector::ZERO, Vector::Y, T).is_err());
    }

    #[test]
    fn a_normal_to_a_tiny_triangle_is_still_well_defined() {
        // The trap: |a x b| is twice the triangle's area, so it scales as the
        // square of the size. Comparing it against a length tolerance rejects
        // small-but-perfectly-valid triangles.
        for scale in [1e-6_f64, 1e-3, 1.0, 1e3] {
            let a = Vector::new(scale, 0.0, 0.0);
            let b = Vector::new(0.0, scale, 0.0);
            let n = Direction::from_cross(a, b, T).unwrap();
            assert!(n.is_equal(Direction::Z, T), "failed at scale {scale}");
        }
        // Whereas the naive route does reject them, which is why it is not used.
        assert!(
            Direction::new(
                Vector::new(1e-6, 0.0, 0.0).cross(Vector::new(0.0, 1e-6, 0.0)),
                T
            )
            .is_err()
        );
    }

    #[test]
    fn any_perpendicular_is_sound_for_every_axis_alignment() {
        // The failure mode this guards: crossing with a fixed axis gives a
        // near-null result when the input happens to be parallel to that axis.
        let cases = [
            Direction::X,
            Direction::Y,
            Direction::Z,
            -Direction::X,
            -Direction::Z,
            Direction::from_coords(1.0, 1.0, 1.0, T).unwrap(),
            Direction::from_coords(1.0, 1e-14, 1e-14, T).unwrap(),
            Direction::from_coords(1e-14, 1e-14, 1.0, T).unwrap(),
        ];
        for d in cases {
            let p = d.any_perpendicular();
            assert_relative_eq!(p.vector().magnitude(), 1.0, epsilon = 1e-14);
            assert_relative_eq!(d.dot(p), 0.0, epsilon = 1e-14);
        }
    }

    #[test]
    fn angle_relations() {
        assert_relative_eq!(Direction::X.angle(Direction::X), 0.0);
        assert_relative_eq!(Direction::X.angle(-Direction::X), core::f64::consts::PI);
        assert_relative_eq!(
            Direction::X.angle(Direction::Y),
            core::f64::consts::FRAC_PI_2
        );
        assert!(Direction::X.is_equal(Direction::X, T));
        assert!(Direction::X.is_opposite(-Direction::X, T));
        assert!(Direction::X.is_parallel(-Direction::X, T));
        assert!(!Direction::X.is_equal(-Direction::X, T));
        assert!(Direction::X.is_normal(Direction::Y, T));
    }

    #[test]
    fn scaling_a_direction_gives_a_free_vector() {
        // The type change is the point: the result is not unit length, so it
        // must not keep claiming to be a direction.
        let v: Vector = Direction::X * 5.0;
        assert_eq!(v, Vector::new(5.0, 0.0, 0.0));
        assert_eq!(5.0 * Direction::X, v);
    }

    #[test]
    fn direction2_angle_round_trips() {
        // Compare directions rather than angles. Angles are only defined modulo
        // 2*pi and `to_angle` has a branch cut, so a direct comparison fails at
        // the cut for reasons that say nothing about correctness: `from_angle`
        // of exactly -pi produces a tiny negative y, which `to_angle` maps back
        // to -pi rather than +pi. Both name the same direction.
        for turns in 0..32 {
            let a = f64::from(turns) * core::f64::consts::PI / 16.0 - core::f64::consts::PI;
            let d = Direction2::from_angle(a);
            assert_relative_eq!(d.vector().magnitude(), 1.0, epsilon = 1e-15);
            assert!(
                Direction2::from_angle(d.to_angle()).is_equal(d, T),
                "round trip failed at {a}"
            );
        }
    }

    #[test]
    fn direction2_perpendicular_is_exact_and_has_period_four() {
        let d = Direction2::from_angle(0.37);
        assert_eq!(
            d.perpendicular()
                .perpendicular()
                .perpendicular()
                .perpendicular(),
            d
        );
        assert_eq!(d.perpendicular().dot(d), 0.0, "exactly zero");
    }

    #[test]
    fn direction2_signed_angle() {
        let quarter = core::f64::consts::FRAC_PI_2;
        assert_relative_eq!(Direction2::X.angle(Direction2::Y), quarter);
        assert_relative_eq!(Direction2::Y.angle(Direction2::X), -quarter);
    }

    #[test]
    fn dimension_round_trip() {
        let d = Direction2::from_angle(0.9);
        let up = d.to_3d();
        assert_relative_eq!(up.z(), 0.0);
        assert!(up.to_2d(T).unwrap().is_equal(d, T));
        // A direction with nothing in the XY plane cannot be projected into it.
        assert!(Direction::Z.to_2d(T).is_err());
    }
}
