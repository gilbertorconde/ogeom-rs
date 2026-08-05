//! Free vectors in 2D and 3D.
//!
//! A vector has magnitude and direction and is unaffected by translation. See
//! [`Direction`](crate::Direction) for the unit-length variant, whose invariant
//! the type system enforces, and [`Point`](crate::Point) for positions.

use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};

/// A free vector in space.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vector {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

/// A free vector in the plane.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vector2 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
}

impl Vector {
    /// The zero vector.
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0);
    /// Unit vector along +X.
    pub const X: Self = Self::new(1.0, 0.0, 0.0);
    /// Unit vector along +Y.
    pub const Y: Self = Self::new(0.0, 1.0, 0.0);
    /// Unit vector along +Z.
    pub const Z: Self = Self::new(0.0, 0.0, 1.0);

    /// A vector from components.
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// A vector with all components equal.
    #[must_use]
    pub const fn splat(v: f64) -> Self {
        Self::new(v, v, v)
    }

    /// Components as an array.
    #[must_use]
    pub const fn to_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// A vector from an array.
    #[must_use]
    pub const fn from_array([x, y, z]: [f64; 3]) -> Self {
        Self::new(x, y, z)
    }

    /// Component by index, `0..3`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Range`](ogeom_core::OgeomError::Range) if `index >= 3`.
    pub fn coord(self, index: usize) -> OgeomResult<f64> {
        match index {
            0 => Ok(self.x),
            1 => Ok(self.y),
            2 => Ok(self.z),
            _ => ogeom_bail!(Range, "vector component {index} of 3"),
        }
    }

    /// Dot product.
    #[must_use]
    pub fn dot(self, other: Self) -> f64 {
        self.x
            .mul_add(other.x, self.y.mul_add(other.y, self.z * other.z))
    }

    /// Cross product. Right-handed: `X.cross(Y) == Z`.
    ///
    /// Each component is a two-term difference `ab - cd`, evaluated without a
    /// fused multiply-add on purpose. An FMA rounds one product and not the
    /// other, which destroys the exact cancellation that makes the cross
    /// product of collinear vectors come out at exactly zero.
    #[must_use]
    pub fn cross(self, other: Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    /// Scalar triple product `self · (a × b)` — the signed volume of the
    /// parallelepiped the three vectors span.
    #[must_use]
    pub fn triple(self, a: Self, b: Self) -> f64 {
        self.dot(a.cross(b))
    }

    /// Squared magnitude. Prefer this to [`Vector::magnitude`] when comparing
    /// lengths: it avoids a square root and the rounding that comes with it.
    #[must_use]
    pub fn square_magnitude(self) -> f64 {
        self.dot(self)
    }

    /// Magnitude.
    #[must_use]
    pub fn magnitude(self) -> f64 {
        self.square_magnitude().sqrt()
    }

    /// This vector scaled to unit length.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// magnitude is below `tol.confusion()`, or if any component is non-finite.
    /// Normalizing a near-zero vector amplifies whatever noise it holds into an
    /// arbitrary direction, so it is refused rather than approximated.
    pub fn normalized(self, tol: Tolerances) -> OgeomResult<Self> {
        if !self.is_finite() {
            ogeom_bail!(Construction, "cannot normalize a non-finite vector");
        }
        let m = self.magnitude();
        if m <= tol.confusion() {
            ogeom_bail!(Construction, "cannot normalize a vector of magnitude {m}");
        }
        Ok(self / m)
    }

    /// Whether every component is finite.
    #[must_use]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    /// Whether the magnitude is within `tol.confusion()` of zero.
    #[must_use]
    pub fn is_zero(self, tol: Tolerances) -> bool {
        self.magnitude() <= tol.confusion()
    }

    /// Whether two vectors agree component-wise within `tol.confusion()`.
    #[must_use]
    pub fn is_equal(self, other: Self, tol: Tolerances) -> bool {
        (self - other).magnitude() <= tol.confusion()
    }

    /// Angle to `other`, in `[0, π]`.
    ///
    /// Uses `atan2` of the cross and dot products rather than `acos` of the
    /// normalized dot product: `acos` loses most of its precision for nearly
    /// parallel or nearly antiparallel vectors, which is exactly where angle
    /// tests matter.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if either
    /// vector is degenerate.
    pub fn angle(self, other: Self, tol: Tolerances) -> OgeomResult<f64> {
        if self.is_zero(tol) || other.is_zero(tol) {
            ogeom_bail!(Construction, "angle is undefined for a null vector");
        }
        Ok(self.cross(other).magnitude().atan2(self.dot(other)))
    }

    /// Angle to `other` measured about `reference`, in `(-π, π]`.
    ///
    /// Positive when `self` turns towards `other` counter-clockwise as seen from
    /// the tip of `reference`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if any vector
    /// is degenerate, or if `reference` is not perpendicular to both within
    /// `tol.angular()`.
    pub fn signed_angle(self, other: Self, reference: Self, tol: Tolerances) -> OgeomResult<f64> {
        let unsigned = self.angle(other, tol)?;
        let normal = self.cross(other);
        if normal.is_zero(tol) {
            // Parallel or antiparallel: the sign is meaningless, and `unsigned`
            // is already exactly 0 or π.
            return Ok(unsigned);
        }
        if reference.is_zero(tol) {
            ogeom_bail!(Construction, "reference vector is null");
        }
        Ok(if normal.dot(reference) < 0.0 {
            -unsigned
        } else {
            unsigned
        })
    }

    /// Whether `self` and `other` point the same way, within `tol.angular()`.
    #[must_use]
    pub fn is_parallel(self, other: Self, tol: Tolerances) -> bool {
        self.angle(other, tol).is_ok_and(|a| a <= tol.angular())
    }

    /// Whether `self` and `other` are parallel, ignoring sense.
    #[must_use]
    pub fn is_collinear(self, other: Self, tol: Tolerances) -> bool {
        self.angle(other, tol)
            .is_ok_and(|a| a <= tol.angular() || (core::f64::consts::PI - a) <= tol.angular())
    }

    /// Whether `self` and `other` are perpendicular, within `tol.angular()`.
    #[must_use]
    pub fn is_normal(self, other: Self, tol: Tolerances) -> bool {
        self.angle(other, tol)
            .is_ok_and(|a| (core::f64::consts::FRAC_PI_2 - a).abs() <= tol.angular())
    }

    /// Component of `self` along `other`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `other` is
    /// degenerate.
    pub fn projected_onto(self, other: Self, tol: Tolerances) -> OgeomResult<Self> {
        let square = other.square_magnitude();
        if square <= tol.confusion() * tol.confusion() {
            ogeom_bail!(Construction, "cannot project onto a null vector");
        }
        Ok(other * (self.dot(other) / square))
    }

    /// Linear interpolation, `t = 0` giving `self`.
    #[must_use]
    pub fn lerp(self, other: Self, t: f64) -> Self {
        self + (other - self) * t
    }

    /// Component-wise minimum.
    #[must_use]
    pub fn min(self, other: Self) -> Self {
        Self::new(
            self.x.min(other.x),
            self.y.min(other.y),
            self.z.min(other.z),
        )
    }

    /// Component-wise maximum.
    #[must_use]
    pub fn max(self, other: Self) -> Self {
        Self::new(
            self.x.max(other.x),
            self.y.max(other.y),
            self.z.max(other.z),
        )
    }

    /// The 2D vector formed by dropping the Z component.
    #[must_use]
    pub const fn xy(self) -> Vector2 {
        Vector2::new(self.x, self.y)
    }
}

impl Vector2 {
    /// The zero vector.
    pub const ZERO: Self = Self::new(0.0, 0.0);
    /// Unit vector along +X.
    pub const X: Self = Self::new(1.0, 0.0);
    /// Unit vector along +Y.
    pub const Y: Self = Self::new(0.0, 1.0);

    /// A vector from components.
    #[must_use]
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Components as an array.
    #[must_use]
    pub const fn to_array(self) -> [f64; 2] {
        [self.x, self.y]
    }

    /// A vector from an array.
    #[must_use]
    pub const fn from_array([x, y]: [f64; 2]) -> Self {
        Self::new(x, y)
    }

    /// Dot product.
    ///
    /// Two terms, so no fused multiply-add — see [`Vector::cross`] for why.
    #[must_use]
    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y
    }

    /// The scalar cross product — the signed area of the parallelogram the two
    /// vectors span. Positive when `other` is counter-clockwise from `self`.
    #[must_use]
    pub fn cross(self, other: Self) -> f64 {
        self.x * other.y - self.y * other.x
    }

    /// Squared magnitude.
    #[must_use]
    pub fn square_magnitude(self) -> f64 {
        self.dot(self)
    }

    /// Magnitude.
    #[must_use]
    pub fn magnitude(self) -> f64 {
        self.square_magnitude().sqrt()
    }

    /// This vector rotated a quarter turn counter-clockwise. Exact — no
    /// trigonometry involved.
    #[must_use]
    pub const fn perpendicular(self) -> Self {
        Self::new(-self.y, self.x)
    }

    /// This vector scaled to unit length.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// magnitude is below `tol.confusion()`, or if any component is non-finite.
    pub fn normalized(self, tol: Tolerances) -> OgeomResult<Self> {
        if !self.is_finite() {
            ogeom_bail!(Construction, "cannot normalize a non-finite vector");
        }
        let m = self.magnitude();
        if m <= tol.confusion() {
            ogeom_bail!(Construction, "cannot normalize a vector of magnitude {m}");
        }
        Ok(self / m)
    }

    /// Whether both components are finite.
    #[must_use]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite()
    }

    /// Whether the magnitude is within `tol.confusion()` of zero.
    #[must_use]
    pub fn is_zero(self, tol: Tolerances) -> bool {
        self.magnitude() <= tol.confusion()
    }

    /// Whether two vectors agree component-wise within `tol.confusion()`.
    #[must_use]
    pub fn is_equal(self, other: Self, tol: Tolerances) -> bool {
        (self - other).magnitude() <= tol.confusion()
    }

    /// Angle to `other`, in `(-π, π]`, positive counter-clockwise.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if either
    /// vector is degenerate.
    pub fn angle(self, other: Self, tol: Tolerances) -> OgeomResult<f64> {
        if self.is_zero(tol) || other.is_zero(tol) {
            ogeom_bail!(Construction, "angle is undefined for a null vector");
        }
        Ok(self.cross(other).atan2(self.dot(other)))
    }

    /// Linear interpolation, `t = 0` giving `self`.
    #[must_use]
    pub fn lerp(self, other: Self, t: f64) -> Self {
        self + (other - self) * t
    }

    /// This vector with a Z component of zero.
    #[must_use]
    pub const fn to_3d(self) -> Vector {
        Vector::new(self.x, self.y, 0.0)
    }
}

macro_rules! impl_vector_ops {
    ($t:ty, $($f:ident),+) => {
        impl Add for $t {
            type Output = Self;
            fn add(self, o: Self) -> Self { Self { $($f: self.$f + o.$f),+ } }
        }
        impl Sub for $t {
            type Output = Self;
            fn sub(self, o: Self) -> Self { Self { $($f: self.$f - o.$f),+ } }
        }
        impl Neg for $t {
            type Output = Self;
            fn neg(self) -> Self { Self { $($f: -self.$f),+ } }
        }
        impl Mul<f64> for $t {
            type Output = Self;
            fn mul(self, s: f64) -> Self { Self { $($f: self.$f * s),+ } }
        }
        impl Mul<$t> for f64 {
            type Output = $t;
            fn mul(self, v: $t) -> $t { v * self }
        }
        impl Div<f64> for $t {
            type Output = Self;
            fn div(self, s: f64) -> Self { Self { $($f: self.$f / s),+ } }
        }
        impl AddAssign for $t {
            fn add_assign(&mut self, o: Self) { *self = *self + o; }
        }
        impl SubAssign for $t {
            fn sub_assign(&mut self, o: Self) { *self = *self - o; }
        }
        impl MulAssign<f64> for $t {
            fn mul_assign(&mut self, s: f64) { *self = *self * s; }
        }
        impl DivAssign<f64> for $t {
            fn div_assign(&mut self, s: f64) { *self = *self / s; }
        }
    };
}

impl_vector_ops!(Vector, x, y, z);
impl_vector_ops!(Vector2, x, y);

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn cross_product_is_right_handed() {
        assert_eq!(Vector::X.cross(Vector::Y), Vector::Z);
        assert_eq!(Vector::Y.cross(Vector::Z), Vector::X);
        assert_eq!(Vector::Z.cross(Vector::X), Vector::Y);
        assert_eq!(Vector::Y.cross(Vector::X), -Vector::Z);
    }

    #[test]
    fn normalizing_a_null_vector_is_refused_not_approximated() {
        assert!(Vector::ZERO.normalized(T).is_err());
        assert!(Vector::new(1e-12, 0.0, 0.0).normalized(T).is_err());
        assert!(Vector::new(f64::NAN, 0.0, 0.0).normalized(T).is_err());
        assert!(Vector::new(f64::INFINITY, 0.0, 0.0).normalized(T).is_err());
        assert!(Vector::new(3.0, 4.0, 0.0).normalized(T).is_ok());
    }

    #[test]
    fn normalized_has_unit_magnitude() {
        let v = Vector::new(3.0, 4.0, 12.0).normalized(T).unwrap();
        assert_relative_eq!(v.magnitude(), 1.0, epsilon = 1e-15);
    }

    #[test]
    fn angle_is_accurate_for_nearly_parallel_vectors() {
        // The case acos-based formulations get wrong: for a tiny angle, the dot
        // product is 1 - O(angle^2), so acos loses half the mantissa. atan2 of
        // cross over dot keeps full precision.
        let tiny: f64 = 1e-9;
        let a = Vector::X;
        let b = Vector::new(tiny.cos(), tiny.sin(), 0.0);
        assert_relative_eq!(a.angle(b, T).unwrap(), tiny, max_relative = 1e-9);
    }

    #[test]
    fn angle_endpoints() {
        assert_relative_eq!(Vector::X.angle(Vector::X, T).unwrap(), 0.0);
        assert_relative_eq!(
            Vector::X.angle(-Vector::X, T).unwrap(),
            core::f64::consts::PI
        );
        assert_relative_eq!(
            Vector::X.angle(Vector::Y, T).unwrap(),
            core::f64::consts::FRAC_PI_2
        );
        assert!(Vector::X.angle(Vector::ZERO, T).is_err());
    }

    #[test]
    fn signed_angle_respects_the_reference_direction() {
        let a = Vector::X;
        let b = Vector::Y;
        let quarter = core::f64::consts::FRAC_PI_2;
        assert_relative_eq!(a.signed_angle(b, Vector::Z, T).unwrap(), quarter);
        assert_relative_eq!(a.signed_angle(b, -Vector::Z, T).unwrap(), -quarter);
        // Antiparallel: no meaningful sign, and pi either way.
        assert_relative_eq!(
            a.signed_angle(-a, Vector::Z, T).unwrap(),
            core::f64::consts::PI
        );
    }

    #[test]
    fn parallel_collinear_and_normal() {
        let a = Vector::new(1.0, 2.0, 3.0);
        assert!(a.is_parallel(a * 5.0, T));
        assert!(!a.is_parallel(a * -5.0, T), "antiparallel is not parallel");
        assert!(a.is_collinear(a * -5.0, T), "but it is collinear");
        assert!(Vector::X.is_normal(Vector::Y, T));
        assert!(!Vector::X.is_normal(Vector::X, T));
    }

    #[test]
    fn projection_onto_an_axis() {
        let v = Vector::new(3.0, 4.0, 5.0);
        let p = v.projected_onto(Vector::X, T).unwrap();
        assert_eq!(p, Vector::new(3.0, 0.0, 0.0));
        // The residual is perpendicular to the axis, by construction.
        assert_relative_eq!((v - p).dot(Vector::X), 0.0, epsilon = 1e-15);
        assert!(v.projected_onto(Vector::ZERO, T).is_err());
    }

    #[test]
    fn triple_product_is_the_signed_volume() {
        assert_relative_eq!(Vector::X.triple(Vector::Y, Vector::Z), 1.0);
        assert_relative_eq!(Vector::X.triple(Vector::Z, Vector::Y), -1.0);
        // Coplanar vectors span no volume.
        assert_relative_eq!(Vector::X.triple(Vector::Y, Vector::new(1.0, 1.0, 0.0)), 0.0);
    }

    #[test]
    fn component_access_is_bounds_checked() {
        let v = Vector::new(1.0, 2.0, 3.0);
        assert_eq!(v.coord(0).unwrap(), 1.0);
        assert_eq!(v.coord(2).unwrap(), 3.0);
        assert!(v.coord(3).is_err());
    }

    #[test]
    fn vector2_cross_is_the_signed_area() {
        assert_relative_eq!(Vector2::X.cross(Vector2::Y), 1.0);
        assert_relative_eq!(Vector2::Y.cross(Vector2::X), -1.0);
        assert_relative_eq!(Vector2::X.cross(Vector2::X), 0.0);
    }

    #[test]
    fn vector2_perpendicular_is_an_exact_quarter_turn() {
        let v = Vector2::new(0.1, 0.7);
        let p = v.perpendicular();
        assert_eq!(p, Vector2::new(-0.7, 0.1));
        assert_eq!(p.dot(v), 0.0, "exactly zero, not merely small");
        assert_eq!(p.perpendicular().perpendicular().perpendicular(), v);
    }

    #[test]
    fn vector2_angle_is_signed() {
        let quarter = core::f64::consts::FRAC_PI_2;
        assert_relative_eq!(Vector2::X.angle(Vector2::Y, T).unwrap(), quarter);
        assert_relative_eq!(Vector2::Y.angle(Vector2::X, T).unwrap(), -quarter);
    }

    #[test]
    fn arithmetic_operators() {
        let a = Vector::new(1.0, 2.0, 3.0);
        let b = Vector::new(4.0, 5.0, 6.0);
        assert_eq!(a + b, Vector::new(5.0, 7.0, 9.0));
        assert_eq!(b - a, Vector::splat(3.0));
        assert_eq!(a * 2.0, Vector::new(2.0, 4.0, 6.0));
        assert_eq!(2.0 * a, a * 2.0);
        assert_eq!(a / 2.0, Vector::new(0.5, 1.0, 1.5));
        let mut c = a;
        c += b;
        c -= b;
        assert_eq!(c, a);
    }

    #[test]
    fn lerp_hits_both_endpoints() {
        let a = Vector::new(1.0, 0.0, 0.0);
        let b = Vector::new(3.0, 4.0, 0.0);
        assert_eq!(a.lerp(b, 0.0), a);
        assert_eq!(a.lerp(b, 1.0), b);
        assert_eq!(a.lerp(b, 0.5), Vector::new(2.0, 2.0, 0.0));
    }
}
