//! Positions in 2D and 3D.
//!
//! A point is affected by translation; a [`Vector`] is not. Keeping them
//! distinct types means the compiler rejects the classic errors — translating a
//! normal, adding two positions — rather than letting them produce plausible
//! nonsense.
//!
//! The algebra that results is the affine one: point − point is a vector, point
//! + vector is a point, and point + point is not defined.

use core::ops::{Add, AddAssign, Sub, SubAssign};

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};

use crate::{Vector, Vector2};

/// A position in space.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
    /// Z coordinate.
    pub z: f64,
}

/// A position in the plane.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point2 {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
}

impl Point {
    /// The origin.
    pub const ORIGIN: Self = Self::new(0.0, 0.0, 0.0);

    /// A point from coordinates.
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Coordinates as an array.
    #[must_use]
    pub const fn to_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// A point from an array.
    #[must_use]
    pub const fn from_array([x, y, z]: [f64; 3]) -> Self {
        Self::new(x, y, z)
    }

    /// The position vector from the origin.
    #[must_use]
    pub const fn to_vector(self) -> Vector {
        Vector::new(self.x, self.y, self.z)
    }

    /// The point at the tip of `v` placed at the origin.
    #[must_use]
    pub const fn from_vector(v: Vector) -> Self {
        Self::new(v.x, v.y, v.z)
    }

    /// Coordinate by index, `0..3`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Range`](ogeom_core::OgeomError::Range) if `index >= 3`.
    pub fn coord(self, index: usize) -> OgeomResult<f64> {
        match index {
            0 => Ok(self.x),
            1 => Ok(self.y),
            2 => Ok(self.z),
            _ => ogeom_bail!(Range, "point coordinate {index} of 3"),
        }
    }

    /// Squared distance to `other`. Prefer this when comparing distances.
    #[must_use]
    pub fn square_distance(self, other: Self) -> f64 {
        (other - self).square_magnitude()
    }

    /// Distance to `other`.
    #[must_use]
    pub fn distance(self, other: Self) -> f64 {
        self.square_distance(other).sqrt()
    }

    /// Whether two points coincide within `tol.confusion()`.
    #[must_use]
    pub fn is_equal(self, other: Self, tol: Tolerances) -> bool {
        self.square_distance(other) <= tol.confusion() * tol.confusion()
    }

    /// Whether two points coincide within an explicit distance.
    #[must_use]
    pub fn is_within(self, other: Self, distance: f64) -> bool {
        self.square_distance(other) <= distance * distance
    }

    /// Whether every coordinate is finite.
    #[must_use]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    /// The midpoint of `self` and `other`.
    #[must_use]
    pub fn midpoint(self, other: Self) -> Self {
        Self::new(
            f64::midpoint(self.x, other.x),
            f64::midpoint(self.y, other.y),
            f64::midpoint(self.z, other.z),
        )
    }

    /// Linear interpolation, `t = 0` giving `self`.
    #[must_use]
    pub fn lerp(self, other: Self, t: f64) -> Self {
        self + (other - self) * t
    }

    /// The centroid of a set of points.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the set is
    /// empty.
    pub fn centroid(points: &[Self]) -> OgeomResult<Self> {
        let Some((first, rest)) = points.split_first() else {
            ogeom_bail!(Construction, "centroid of an empty point set");
        };
        // Accumulate offsets from the first point rather than absolute
        // coordinates: for a cluster far from the origin the absolute sum loses
        // precision to cancellation, while the offsets stay small.
        let mut sum = Vector::ZERO;
        for p in rest {
            sum += *p - *first;
        }
        #[allow(clippy::cast_precision_loss)]
        Ok(*first + sum / points.len() as f64)
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

    /// This point with the Z coordinate dropped.
    #[must_use]
    pub const fn xy(self) -> Point2 {
        Point2::new(self.x, self.y)
    }
}

impl Point2 {
    /// The origin.
    pub const ORIGIN: Self = Self::new(0.0, 0.0);

    /// A point from coordinates.
    #[must_use]
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Coordinates as an array.
    #[must_use]
    pub const fn to_array(self) -> [f64; 2] {
        [self.x, self.y]
    }

    /// A point from an array.
    #[must_use]
    pub const fn from_array([x, y]: [f64; 2]) -> Self {
        Self::new(x, y)
    }

    /// The position vector from the origin.
    #[must_use]
    pub const fn to_vector(self) -> Vector2 {
        Vector2::new(self.x, self.y)
    }

    /// The point at the tip of `v` placed at the origin.
    #[must_use]
    pub const fn from_vector(v: Vector2) -> Self {
        Self::new(v.x, v.y)
    }

    /// Squared distance to `other`.
    #[must_use]
    pub fn square_distance(self, other: Self) -> f64 {
        (other - self).square_magnitude()
    }

    /// Distance to `other`.
    #[must_use]
    pub fn distance(self, other: Self) -> f64 {
        self.square_distance(other).sqrt()
    }

    /// Whether two points coincide within `tol.confusion()`.
    #[must_use]
    pub fn is_equal(self, other: Self, tol: Tolerances) -> bool {
        self.square_distance(other) <= tol.confusion() * tol.confusion()
    }

    /// Whether both coordinates are finite.
    #[must_use]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite()
    }

    /// The midpoint of `self` and `other`.
    #[must_use]
    pub fn midpoint(self, other: Self) -> Self {
        Self::new(
            f64::midpoint(self.x, other.x),
            f64::midpoint(self.y, other.y),
        )
    }

    /// Linear interpolation, `t = 0` giving `self`.
    #[must_use]
    pub fn lerp(self, other: Self, t: f64) -> Self {
        self + (other - self) * t
    }

    /// This point in the XY plane of space.
    #[must_use]
    pub const fn to_3d(self) -> Point {
        Point::new(self.x, self.y, 0.0)
    }
}

impl Add<Vector> for Point {
    type Output = Self;
    fn add(self, v: Vector) -> Self {
        Self::new(self.x + v.x, self.y + v.y, self.z + v.z)
    }
}

impl Sub<Vector> for Point {
    type Output = Self;
    fn sub(self, v: Vector) -> Self {
        Self::new(self.x - v.x, self.y - v.y, self.z - v.z)
    }
}

impl Sub for Point {
    type Output = Vector;
    /// The displacement from `other` to `self`.
    fn sub(self, other: Self) -> Vector {
        Vector::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl AddAssign<Vector> for Point {
    fn add_assign(&mut self, v: Vector) {
        *self = *self + v;
    }
}

impl SubAssign<Vector> for Point {
    fn sub_assign(&mut self, v: Vector) {
        *self = *self - v;
    }
}

impl Add<Vector2> for Point2 {
    type Output = Self;
    fn add(self, v: Vector2) -> Self {
        Self::new(self.x + v.x, self.y + v.y)
    }
}

impl Sub<Vector2> for Point2 {
    type Output = Self;
    fn sub(self, v: Vector2) -> Self {
        Self::new(self.x - v.x, self.y - v.y)
    }
}

impl Sub for Point2 {
    type Output = Vector2;
    fn sub(self, other: Self) -> Vector2 {
        Vector2::new(self.x - other.x, self.y - other.y)
    }
}

impl AddAssign<Vector2> for Point2 {
    fn add_assign(&mut self, v: Vector2) {
        *self = *self + v;
    }
}

impl SubAssign<Vector2> for Point2 {
    fn sub_assign(&mut self, v: Vector2) {
        *self = *self - v;
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn affine_algebra() {
        let a = Point::new(1.0, 2.0, 3.0);
        let b = Point::new(4.0, 6.0, 3.0);
        let d: Vector = b - a;
        assert_eq!(d, Vector::new(3.0, 4.0, 0.0));
        assert_eq!(a + d, b);
        assert_eq!(b - d, a);
        assert_relative_eq!(a.distance(b), 5.0);
        assert_relative_eq!(a.square_distance(b), 25.0);
    }

    #[test]
    fn midpoint_does_not_overflow_for_extreme_coordinates() {
        // The naive (a + b) / 2 overflows to infinity here; f64::midpoint does
        // not. Coordinates this large are pathological, but a kernel that
        // produces infinities on them is worse than one that does not.
        let a = Point::new(f64::MAX, 0.0, 0.0);
        let b = Point::new(f64::MAX, 0.0, 0.0);
        assert!(a.midpoint(b).is_finite());
        assert_eq!(a.midpoint(b).x, f64::MAX);
    }

    #[test]
    fn midpoint_and_lerp_agree() {
        let a = Point::new(-1.0, 5.0, 2.0);
        let b = Point::new(3.0, 1.0, -4.0);
        assert_eq!(a.midpoint(b), a.lerp(b, 0.5));
        assert_eq!(a.lerp(b, 0.0), a);
        assert_eq!(a.lerp(b, 1.0), b);
    }

    #[test]
    fn centroid_keeps_precision_far_from_the_origin() {
        // Summing absolute coordinates around 1e9 and dividing loses the
        // millimetre-scale detail we care about. Summing offsets does not.
        let base = 1.0e9;
        let pts = [
            Point::new(base, base, base),
            Point::new(base + 2.0, base, base),
            Point::new(base + 1.0, base + 3.0, base),
        ];
        let c = Point::centroid(&pts).unwrap();
        assert_relative_eq!(c.x, base + 1.0, epsilon = 1e-9);
        assert_relative_eq!(c.y, base + 1.0, epsilon = 1e-9);
        assert_relative_eq!(c.z, base, epsilon = 1e-9);
    }

    #[test]
    fn centroid_of_nothing_is_an_error_not_the_origin() {
        assert!(Point::centroid(&[]).is_err());
        let single = Point::new(1.0, 2.0, 3.0);
        assert_eq!(Point::centroid(&[single]).unwrap(), single);
    }

    #[test]
    fn equality_uses_tolerance() {
        let a = Point::new(1.0, 1.0, 1.0);
        let near = Point::new(1.0 + 1e-9, 1.0, 1.0);
        let far = Point::new(1.0 + 1e-3, 1.0, 1.0);
        assert!(a.is_equal(near, T));
        assert!(!a.is_equal(far, T));
        assert!(a.is_within(far, 1e-2));
    }

    #[test]
    fn coordinate_access_is_bounds_checked() {
        let p = Point::new(7.0, 8.0, 9.0);
        assert_eq!(p.coord(1).unwrap(), 8.0);
        assert!(p.coord(3).is_err());
    }

    #[test]
    fn dimension_round_trip() {
        let p = Point::new(1.0, 2.0, 3.0);
        assert_eq!(p.xy(), Point2::new(1.0, 2.0));
        assert_eq!(p.xy().to_3d(), Point::new(1.0, 2.0, 0.0));
    }

    #[test]
    fn point2_affine_algebra() {
        let a = Point2::new(1.0, 2.0);
        let b = Point2::new(4.0, 6.0);
        assert_eq!(b - a, Vector2::new(3.0, 4.0));
        assert_relative_eq!(a.distance(b), 5.0);
        assert_eq!(a + (b - a), b);
    }
}
