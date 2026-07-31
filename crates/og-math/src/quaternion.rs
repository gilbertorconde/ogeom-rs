//! Unit quaternions for rotation.
//!
//! Preferred over matrices wherever rotations are *composed* or *interpolated*:
//! repeated matrix products drift away from orthonormality and have to be
//! re-orthonormalized, while a quaternion only needs renormalizing, and there is
//! no meaningful way to interpolate two rotation matrices directly.
//!
//! [`Quaternion`] is not constrained to unit length by construction — the
//! arithmetic needs unnormalized intermediates — but every rotation operation
//! either requires or restores unit length, and says which in its
//! documentation.

use core::ops::{Add, Mul, Neg, Sub};

use og_core::{OgResult, Tolerances, og_bail};

use crate::{Direction, Matrix3, Vector};

/// How far a matrix may stray from orthonormal and still be accepted as a
/// rotation. Dimensionless, and generous enough to admit a matrix assembled
/// from a chain of rotations without admitting a scaling.
const ORTHONORMAL_EPS: f64 = 1e-10;

/// A quaternion `w + xi + yj + zk`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quaternion {
    /// Scalar part.
    pub w: f64,
    /// `i` coefficient.
    pub x: f64,
    /// `j` coefficient.
    pub y: f64,
    /// `k` coefficient.
    pub z: f64,
}

impl Default for Quaternion {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Quaternion {
    /// The identity rotation.
    pub const IDENTITY: Self = Self::new(1.0, 0.0, 0.0, 0.0);

    /// From components.
    #[must_use]
    pub const fn new(w: f64, x: f64, y: f64, z: f64) -> Self {
        Self { w, x, y, z }
    }

    /// The rotation of `angle` radians about `axis`, right-handed.
    #[must_use]
    pub fn from_axis_angle(axis: Direction, angle: f64) -> Self {
        let (s, c) = (angle * 0.5).sin_cos();
        Self::new(c, axis.x() * s, axis.y() * s, axis.z() * s)
    }

    /// The shortest rotation taking `from` to `to`.
    ///
    /// # Errors
    ///
    /// [`OgError::Construction`](og_core::OgError::Construction) if the two are
    /// antiparallel: infinitely many shortest rotations exist and picking one
    /// arbitrarily would make the result depend on unobservable rounding.
    pub fn between(from: Direction, to: Direction, tol: Tolerances) -> OgResult<Self> {
        let d = from.dot(to);
        if d < -1.0 + tol.confusion() {
            og_bail!(
                Construction,
                "rotation between antiparallel directions is not unique"
            );
        }
        let axis = from.cross_vector(to);
        // w = 1 + cos(theta), (x,y,z) = sin(theta) * axis. Normalizing this
        // halves the angle, which is what the quaternion needs.
        Self::new(1.0 + d, axis.x, axis.y, axis.z).normalized(tol)
    }

    /// From a rotation matrix.
    ///
    /// Uses Shepperd's method: pick the largest of the four possible divisors so
    /// the division is never by something near zero. The naive `w`-first
    /// formulation loses precision for rotations near π, where `w → 0`.
    ///
    /// # Errors
    ///
    /// [`OgError::Construction`](og_core::OgError::Construction) if `m` is not
    /// orthonormal with determinant `+1`.
    pub fn from_matrix(m: &Matrix3, tol: Tolerances) -> OgResult<Self> {
        if !m.is_orthonormal(ORTHONORMAL_EPS) {
            og_bail!(Construction, "rotation matrix is not orthonormal");
        }
        if m.determinant() < 0.0 {
            og_bail!(Construction, "matrix is a reflection, not a rotation");
        }
        let r = &m.rows;
        let trace = m.trace();
        let q = if trace > 0.0 {
            let s = (trace + 1.0).sqrt() * 2.0;
            Self::new(
                0.25 * s,
                (r[2][1] - r[1][2]) / s,
                (r[0][2] - r[2][0]) / s,
                (r[1][0] - r[0][1]) / s,
            )
        } else if r[0][0] > r[1][1] && r[0][0] > r[2][2] {
            let s = (1.0 + r[0][0] - r[1][1] - r[2][2]).sqrt() * 2.0;
            Self::new(
                (r[2][1] - r[1][2]) / s,
                0.25 * s,
                (r[0][1] + r[1][0]) / s,
                (r[0][2] + r[2][0]) / s,
            )
        } else if r[1][1] > r[2][2] {
            let s = (1.0 - r[0][0] + r[1][1] - r[2][2]).sqrt() * 2.0;
            Self::new(
                (r[0][2] - r[2][0]) / s,
                (r[0][1] + r[1][0]) / s,
                0.25 * s,
                (r[1][2] + r[2][1]) / s,
            )
        } else {
            let s = (1.0 - r[0][0] - r[1][1] + r[2][2]).sqrt() * 2.0;
            Self::new(
                (r[1][0] - r[0][1]) / s,
                (r[0][2] + r[2][0]) / s,
                (r[1][2] + r[2][1]) / s,
                0.25 * s,
            )
        };
        q.normalized(tol)
    }

    /// Squared norm.
    #[must_use]
    pub fn square_norm(self) -> f64 {
        self.w.mul_add(
            self.w,
            self.x
                .mul_add(self.x, self.y.mul_add(self.y, self.z * self.z)),
        )
    }

    /// Norm.
    #[must_use]
    pub fn norm(self) -> f64 {
        self.square_norm().sqrt()
    }

    /// This quaternion scaled to unit norm.
    ///
    /// # Errors
    ///
    /// [`OgError::Construction`](og_core::OgError::Construction) if the norm is
    /// below `tol.confusion()` or any component is non-finite.
    pub fn normalized(self, tol: Tolerances) -> OgResult<Self> {
        if !self.is_finite() {
            og_bail!(Construction, "cannot normalize a non-finite quaternion");
        }
        let n = self.norm();
        if n <= tol.confusion() {
            og_bail!(Construction, "cannot normalize a quaternion of norm {n}");
        }
        Ok(Self::new(self.w / n, self.x / n, self.y / n, self.z / n))
    }

    /// Whether every component is finite.
    #[must_use]
    pub fn is_finite(self) -> bool {
        self.w.is_finite() && self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    /// The conjugate. For a unit quaternion this is the inverse rotation.
    #[must_use]
    pub const fn conjugate(self) -> Self {
        Self::new(self.w, -self.x, -self.y, -self.z)
    }

    /// The multiplicative inverse.
    ///
    /// # Errors
    ///
    /// [`OgError::Numeric`](og_core::OgError::Numeric) if the norm is
    /// degenerate.
    pub fn inverse(self, tol: Tolerances) -> OgResult<Self> {
        let n2 = self.square_norm();
        if n2 <= tol.confusion() * tol.confusion() {
            og_bail!(Numeric, "quaternion of norm {} has no inverse", n2.sqrt());
        }
        let c = self.conjugate();
        Ok(Self::new(c.w / n2, c.x / n2, c.y / n2, c.z / n2))
    }

    /// Apply this rotation to a vector. Assumes unit norm.
    ///
    /// Evaluated as `v + 2w(u × v) + 2(u × (u × v))` with `u` the vector part,
    /// which costs fewer operations than building the matrix and avoids the
    /// intermediate rounding of a full quaternion sandwich product.
    #[must_use]
    pub fn rotate(self, v: Vector) -> Vector {
        let u = Vector::new(self.x, self.y, self.z);
        let uv = u.cross(v);
        v + uv * (2.0 * self.w) + u.cross(uv) * 2.0
    }

    /// The equivalent rotation matrix. Assumes unit norm.
    #[must_use]
    pub fn to_matrix(self) -> Matrix3 {
        let (w, x, y, z) = (self.w, self.x, self.y, self.z);
        let (xx, yy, zz) = (x * x, y * y, z * z);
        let (xy, xz, yz) = (x * y, x * z, y * z);
        let (wx, wy, wz) = (w * x, w * y, w * z);
        Matrix3::new([
            [
                (-2.0f64).mul_add(yy + zz, 1.0),
                2.0 * (xy - wz),
                2.0 * (xz + wy),
            ],
            [
                2.0 * (xy + wz),
                (-2.0f64).mul_add(xx + zz, 1.0),
                2.0 * (yz - wx),
            ],
            [
                2.0 * (xz - wy),
                2.0 * (yz + wx),
                (-2.0f64).mul_add(xx + yy, 1.0),
            ],
        ])
    }

    /// The rotation axis and angle. Assumes unit norm; angle is in `[0, π]`.
    ///
    /// # Errors
    ///
    /// [`OgError::Construction`](og_core::OgError::Construction) if the rotation
    /// is the identity, where the axis is undefined.
    pub fn to_axis_angle(self, tol: Tolerances) -> OgResult<(Direction, f64)> {
        let u = Vector::new(self.x, self.y, self.z);
        let sin_half = u.magnitude();
        if sin_half <= tol.confusion() {
            og_bail!(Construction, "identity rotation has no defined axis");
        }
        // atan2 rather than acos(w) or asin: accurate across the whole range,
        // including angles near 0 and near pi.
        Ok((Direction::new(u, tol)?, 2.0 * sin_half.atan2(self.w)))
    }

    /// Dot product, as 4-vectors.
    #[must_use]
    pub fn dot(self, o: Self) -> f64 {
        self.w
            .mul_add(o.w, self.x.mul_add(o.x, self.y.mul_add(o.y, self.z * o.z)))
    }

    /// Spherical linear interpolation, `t = 0` giving `self`.
    ///
    /// Takes the shorter of the two arcs, and falls back to normalized linear
    /// interpolation when the two are nearly coincident, where the `sin` in the
    /// denominator of the spherical form goes to zero.
    ///
    /// # Errors
    ///
    /// [`OgError::Construction`](og_core::OgError::Construction) if either
    /// operand cannot be normalized.
    pub fn slerp(self, other: Self, t: f64, tol: Tolerances) -> OgResult<Self> {
        let a = self.normalized(tol)?;
        let mut b = other.normalized(tol)?;
        let mut cos = a.dot(b);
        // q and -q are the same rotation; pick the representative that gives the
        // shorter path.
        if cos < 0.0 {
            b = -b;
            cos = -cos;
        }
        if cos > 1.0 - 1e-9 {
            return (a * (1.0 - t) + b * t).normalized(tol);
        }
        let theta = cos.clamp(-1.0, 1.0).acos();
        let sin = theta.sin();
        let wa = ((1.0 - t) * theta).sin() / sin;
        let wb = (t * theta).sin() / sin;
        (a * wa + b * wb).normalized(tol)
    }

    /// The rotation angle to `other`, in `[0, π]`. Assumes unit norms.
    #[must_use]
    pub fn angle_to(self, other: Self) -> f64 {
        let cos = self.dot(other).abs().clamp(-1.0, 1.0);
        2.0 * cos.acos()
    }
}

impl Mul for Quaternion {
    type Output = Self;
    /// Hamilton product. `(a * b)` rotates by `b` first, then by `a`.
    fn mul(self, o: Self) -> Self {
        Self::new(
            self.w
                .mul_add(o.w, -self.x.mul_add(o.x, self.y.mul_add(o.y, self.z * o.z))),
            self.w.mul_add(
                o.x,
                self.x.mul_add(o.w, self.y.mul_add(o.z, -(self.z * o.y))),
            ),
            self.w.mul_add(
                o.y,
                self.y.mul_add(o.w, self.z.mul_add(o.x, -(self.x * o.z))),
            ),
            self.w.mul_add(
                o.z,
                self.z.mul_add(o.w, self.x.mul_add(o.y, -(self.y * o.x))),
            ),
        )
    }
}

impl Mul<f64> for Quaternion {
    type Output = Self;
    fn mul(self, s: f64) -> Self {
        Self::new(self.w * s, self.x * s, self.y * s, self.z * s)
    }
}

impl Add for Quaternion {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.w + o.w, self.x + o.x, self.y + o.y, self.z + o.z)
    }
}

impl Sub for Quaternion {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.w - o.w, self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

impl Neg for Quaternion {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.w, -self.x, -self.y, -self.z)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

    fn axis() -> Direction {
        Direction::from_coords(1.0, 2.0, -1.0, T).unwrap()
    }

    #[test]
    fn rotate_agrees_with_the_matrix_form() {
        let q = Quaternion::from_axis_angle(axis(), 1.234);
        let m = q.to_matrix();
        for v in [
            Vector::X,
            Vector::new(3.0, -2.0, 5.0),
            Vector::new(-1e3, 1e-3, 7.0),
        ] {
            assert!(
                q.rotate(v).is_equal(m * v, T),
                "quaternion and matrix disagree"
            );
        }
    }

    #[test]
    fn rotation_preserves_length_and_is_invertible() {
        let q = Quaternion::from_axis_angle(axis(), 2.1);
        let v = Vector::new(1.0, 2.0, 3.0);
        assert_relative_eq!(q.rotate(v).magnitude(), v.magnitude(), epsilon = 1e-14);
        assert!(q.conjugate().rotate(q.rotate(v)).is_equal(v, T));
        assert!(q.inverse(T).unwrap().rotate(q.rotate(v)).is_equal(v, T));
    }

    #[test]
    fn composition_applies_right_to_left() {
        let a = Quaternion::from_axis_angle(Direction::Z, core::f64::consts::FRAC_PI_2);
        let b = Quaternion::from_axis_angle(Direction::X, core::f64::consts::FRAC_PI_2);
        let v = Vector::Y;
        assert!((a * b).rotate(v).is_equal(a.rotate(b.rotate(v)), T));
    }

    #[test]
    fn matrix_round_trip_is_exact_near_pi() {
        // The naive w-first extraction divides by something approaching zero
        // here. Shepperd's largest-divisor choice does not.
        // Including angles right at pi, where w -> 0 and the naive extraction
        // divides by something vanishing.
        let near_pi = core::f64::consts::PI - 1e-8;
        for angle in [0.0_f64, 0.1, 1.0, 3.0, near_pi, core::f64::consts::PI] {
            let q = Quaternion::from_axis_angle(axis(), angle);
            let back = Quaternion::from_matrix(&q.to_matrix(), T).unwrap();
            // q and -q are the same rotation, so compare the rotations.
            assert!(
                back.to_matrix().is_equal(&q.to_matrix(), 1e-12),
                "round trip failed at angle {angle}"
            );
        }
    }

    #[test]
    fn reflections_are_rejected_as_rotations() {
        let m = Matrix3::reflection(Direction::Z);
        assert!(Quaternion::from_matrix(&m, T).is_err());
        assert!(Quaternion::from_matrix(&Matrix3::scaling(2.0), T).is_err());
    }

    #[test]
    fn axis_angle_round_trip() {
        let a = axis();
        for angle in [0.01_f64, 0.5, 1.5, 3.0] {
            let q = Quaternion::from_axis_angle(a, angle);
            let (back_axis, back_angle) = q.to_axis_angle(T).unwrap();
            assert_relative_eq!(back_angle, angle, epsilon = 1e-12);
            assert!(back_axis.is_equal(a, T));
        }
        assert!(Quaternion::IDENTITY.to_axis_angle(T).is_err());
    }

    #[test]
    fn between_gives_the_shortest_rotation() {
        let from = Direction::X;
        let to = Direction::from_coords(1.0, 1.0, 0.0, T).unwrap();
        let q = Quaternion::between(from, to, T).unwrap();
        assert!(
            Direction::new(q.rotate(from.vector()), T)
                .unwrap()
                .is_equal(to, T)
        );
        let (_, angle) = q.to_axis_angle(T).unwrap();
        assert_relative_eq!(angle, core::f64::consts::FRAC_PI_4, epsilon = 1e-12);
    }

    #[test]
    fn between_refuses_the_ambiguous_antiparallel_case() {
        assert!(Quaternion::between(Direction::X, -Direction::X, T).is_err());
        // Identical directions are fine: the answer is the identity.
        let q = Quaternion::between(Direction::X, Direction::X, T).unwrap();
        assert!(q.rotate(Vector::Y).is_equal(Vector::Y, T));
    }

    #[test]
    fn slerp_hits_the_endpoints_and_stays_unit() {
        let a = Quaternion::from_axis_angle(Direction::Z, 0.2);
        let b = Quaternion::from_axis_angle(Direction::X, 1.9);
        assert!(
            a.slerp(b, 0.0, T)
                .unwrap()
                .to_matrix()
                .is_equal(&a.to_matrix(), 1e-12)
        );
        assert!(
            a.slerp(b, 1.0, T)
                .unwrap()
                .to_matrix()
                .is_equal(&b.to_matrix(), 1e-12)
        );
        for i in 0..=10 {
            let q = a.slerp(b, f64::from(i) / 10.0, T).unwrap();
            assert_relative_eq!(q.norm(), 1.0, epsilon = 1e-14);
        }
    }

    #[test]
    fn slerp_takes_the_short_way_round() {
        let a = Quaternion::from_axis_angle(Direction::Z, 0.0);
        // Same rotation, opposite representative. Naive slerp would sweep the
        // long way; the sign correction must prevent that.
        let b = -Quaternion::from_axis_angle(Direction::Z, 0.4);
        let mid = a.slerp(b, 0.5, T).unwrap();
        let expected = Quaternion::from_axis_angle(Direction::Z, 0.2);
        assert!(mid.to_matrix().is_equal(&expected.to_matrix(), 1e-12));
    }

    #[test]
    fn slerp_survives_nearly_identical_inputs() {
        let a = Quaternion::from_axis_angle(Direction::Z, 1.0);
        let b = Quaternion::from_axis_angle(Direction::Z, 1.0 + 1e-12);
        let mid = a.slerp(b, 0.5, T).unwrap();
        assert!(mid.is_finite());
        assert_relative_eq!(mid.norm(), 1.0, epsilon = 1e-14);
    }

    #[test]
    fn degenerate_quaternions_are_refused() {
        let zero = Quaternion::new(0.0, 0.0, 0.0, 0.0);
        assert!(zero.normalized(T).is_err());
        assert!(zero.inverse(T).is_err());
        assert!(
            Quaternion::new(f64::NAN, 0.0, 0.0, 1.0)
                .normalized(T)
                .is_err()
        );
    }

    #[test]
    fn angle_to_ignores_the_sign_representative() {
        let a = Quaternion::from_axis_angle(Direction::Z, 0.0);
        let b = Quaternion::from_axis_angle(Direction::Z, 1.0);
        assert_relative_eq!(a.angle_to(b), 1.0, epsilon = 1e-12);
        assert_relative_eq!(a.angle_to(-b), 1.0, epsilon = 1e-12);
    }
}
