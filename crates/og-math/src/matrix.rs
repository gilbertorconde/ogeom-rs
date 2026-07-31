//! Small dense matrices — the linear part of a transform.
//!
//! Fixed 3×3 and 2×2, row-major, stored inline. Deliberately not general: these
//! sit on the hot path of every transform application in the kernel, and the
//! sizes are known. General linear algebra lives in [`crate::solve`].

use core::ops::{Add, Mul, Neg, Sub};

use og_core::{OgResult, og_bail};

use crate::{Direction, Vector, Vector2};

/// A 3×3 matrix, row-major.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix3 {
    /// Rows, each `[m[i][0], m[i][1], m[i][2]]`.
    pub rows: [[f64; 3]; 3],
}

/// A 2×2 matrix, row-major.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix2 {
    /// Rows.
    pub rows: [[f64; 2]; 2],
}

impl Default for Matrix3 {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Default for Matrix2 {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Matrix3 {
    /// The identity.
    pub const IDENTITY: Self = Self::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    /// All zeros.
    pub const ZERO: Self = Self::new([[0.0; 3]; 3]);

    /// From rows.
    #[must_use]
    pub const fn new(rows: [[f64; 3]; 3]) -> Self {
        Self { rows }
    }

    /// From column vectors.
    #[must_use]
    pub const fn from_columns(a: Vector, b: Vector, c: Vector) -> Self {
        Self::new([[a.x, b.x, c.x], [a.y, b.y, c.y], [a.z, b.z, c.z]])
    }

    /// From row vectors.
    #[must_use]
    pub const fn from_rows(a: Vector, b: Vector, c: Vector) -> Self {
        Self::new([[a.x, a.y, a.z], [b.x, b.y, b.z], [c.x, c.y, c.z]])
    }

    /// A uniform scaling.
    #[must_use]
    pub const fn scaling(s: f64) -> Self {
        Self::new([[s, 0.0, 0.0], [0.0, s, 0.0], [0.0, 0.0, s]])
    }

    /// A non-uniform scaling along the coordinate axes.
    #[must_use]
    pub const fn scaling_xyz(x: f64, y: f64, z: f64) -> Self {
        Self::new([[x, 0.0, 0.0], [0.0, y, 0.0], [0.0, 0.0, z]])
    }

    /// A right-handed rotation of `angle` radians about `axis`.
    ///
    /// Rodrigues' formula. Orthonormal to within rounding for any unit `axis`,
    /// which [`Direction`] guarantees.
    #[must_use]
    pub fn rotation(axis: Direction, angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        let t = 1.0 - c;
        let (x, y, z) = (axis.x(), axis.y(), axis.z());
        Self::new([
            [
                t.mul_add(x * x, c),
                t.mul_add(x * y, -(s * z)),
                t.mul_add(x * z, s * y),
            ],
            [
                t.mul_add(x * y, s * z),
                t.mul_add(y * y, c),
                t.mul_add(y * z, -(s * x)),
            ],
            [
                t.mul_add(x * z, -(s * y)),
                t.mul_add(y * z, s * x),
                t.mul_add(z * z, c),
            ],
        ])
    }

    /// Reflection in the plane through the origin with normal `n`.
    #[must_use]
    pub fn reflection(n: Direction) -> Self {
        let (x, y, z) = (n.x(), n.y(), n.z());
        Self::new([
            [(-2.0f64).mul_add(x * x, 1.0), -2.0 * x * y, -2.0 * x * z],
            [-2.0 * x * y, (-2.0f64).mul_add(y * y, 1.0), -2.0 * y * z],
            [-2.0 * x * z, -2.0 * y * z, (-2.0f64).mul_add(z * z, 1.0)],
        ])
    }

    /// Element at `(row, col)`.
    ///
    /// # Errors
    ///
    /// [`OgError::Range`](og_core::OgError::Range) if either index exceeds 2.
    pub fn get(&self, row: usize, col: usize) -> OgResult<f64> {
        if row > 2 || col > 2 {
            og_bail!(Range, "matrix index ({row}, {col}) of 3x3");
        }
        Ok(self.rows[row][col])
    }

    /// Column `i` as a vector.
    ///
    /// # Errors
    ///
    /// [`OgError::Range`](og_core::OgError::Range) if `i > 2`.
    pub fn column(&self, i: usize) -> OgResult<Vector> {
        if i > 2 {
            og_bail!(Range, "matrix column {i} of 3");
        }
        Ok(Vector::new(
            self.rows[0][i],
            self.rows[1][i],
            self.rows[2][i],
        ))
    }

    /// Row `i` as a vector.
    ///
    /// # Errors
    ///
    /// [`OgError::Range`](og_core::OgError::Range) if `i > 2`.
    pub fn row(&self, i: usize) -> OgResult<Vector> {
        if i > 2 {
            og_bail!(Range, "matrix row {i} of 3");
        }
        Ok(Vector::from_array(self.rows[i]))
    }

    /// The transpose.
    #[must_use]
    pub const fn transposed(&self) -> Self {
        let m = &self.rows;
        Self::new([
            [m[0][0], m[1][0], m[2][0]],
            [m[0][1], m[1][1], m[2][1]],
            [m[0][2], m[1][2], m[2][2]],
        ])
    }

    /// The determinant.
    #[must_use]
    pub fn determinant(&self) -> f64 {
        let m = &self.rows;
        m[0][0].mul_add(
            m[1][1].mul_add(m[2][2], -(m[1][2] * m[2][1])),
            m[0][1].mul_add(
                -m[1][0].mul_add(m[2][2], -(m[1][2] * m[2][0])),
                m[0][2] * m[1][0].mul_add(m[2][1], -(m[1][1] * m[2][0])),
            ),
        )
    }

    /// The trace.
    #[must_use]
    pub fn trace(&self) -> f64 {
        self.rows[0][0] + self.rows[1][1] + self.rows[2][2]
    }

    /// The inverse.
    ///
    /// # Errors
    ///
    /// [`OgError::Numeric`](og_core::OgError::Numeric) if the matrix is
    /// singular.
    ///
    /// Singularity is a question of numerical conditioning, not of modelling
    /// tolerance, so this takes no [`Tolerances`](og_core::Tolerances). The determinant is compared
    /// against the rounding error incurred computing it: for an `n x n` matrix
    /// with entries bounded by `s`, that is on the order of `n * n! * eps *
    /// s^n`. The bound scales with the entries and has no absolute floor,
    /// because a matrix scaled by `1e-3` has a determinant scaled by `1e-9` and
    /// is no less invertible for it.
    pub fn inverse(&self) -> OgResult<Self> {
        let d = self.determinant();
        let scale = self
            .rows
            .iter()
            .flatten()
            .fold(0.0_f64, |acc, v| acc.max(v.abs()));
        if d.abs() <= 18.0 * f64::EPSILON * scale * scale * scale {
            og_bail!(Numeric, "matrix is singular (determinant {d})");
        }
        let m = &self.rows;
        let cof =
            |a: usize, b: usize, c: usize, e: usize| m[a][b].mul_add(m[c][e], -(m[a][e] * m[c][b]));
        // Adjugate (transposed cofactor matrix), divided by the determinant.
        Ok(Self::new([
            [
                cof(1, 1, 2, 2) / d,
                -cof(0, 1, 2, 2) / d,
                cof(0, 1, 1, 2) / d,
            ],
            [
                -cof(1, 0, 2, 2) / d,
                cof(0, 0, 2, 2) / d,
                -cof(0, 0, 1, 2) / d,
            ],
            [
                cof(1, 0, 2, 1) / d,
                -cof(0, 0, 2, 1) / d,
                cof(0, 0, 1, 1) / d,
            ],
        ]))
    }

    /// Whether the matrix is orthonormal — its rows form an orthonormal basis,
    /// so it represents a rotation or a reflection and its inverse is its
    /// transpose.
    #[must_use]
    pub fn is_orthonormal(&self, eps: f64) -> bool {
        let p = *self * self.transposed();
        p.is_equal(&Self::IDENTITY, eps)
    }

    /// Whether every element agrees with `other` within `eps`.
    ///
    /// Takes a bare epsilon rather than [`Tolerances`](og_core::Tolerances): matrix entries are
    /// dimensionless ratios, and a length tolerance does not apply to them.
    #[must_use]
    pub fn is_equal(&self, other: &Self, eps: f64) -> bool {
        self.rows
            .iter()
            .flatten()
            .zip(other.rows.iter().flatten())
            .all(|(a, b)| (a - b).abs() <= eps)
    }

    /// Whether every element is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.rows.iter().flatten().all(|v| v.is_finite())
    }
}

impl Matrix2 {
    /// The identity.
    pub const IDENTITY: Self = Self::new([[1.0, 0.0], [0.0, 1.0]]);
    /// All zeros.
    pub const ZERO: Self = Self::new([[0.0; 2]; 2]);

    /// From rows.
    #[must_use]
    pub const fn new(rows: [[f64; 2]; 2]) -> Self {
        Self { rows }
    }

    /// A counter-clockwise rotation of `angle` radians.
    #[must_use]
    pub fn rotation(angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        Self::new([[c, -s], [s, c]])
    }

    /// A uniform scaling.
    #[must_use]
    pub const fn scaling(s: f64) -> Self {
        Self::new([[s, 0.0], [0.0, s]])
    }

    /// The transpose.
    #[must_use]
    pub const fn transposed(&self) -> Self {
        let m = &self.rows;
        Self::new([[m[0][0], m[1][0]], [m[0][1], m[1][1]]])
    }

    /// The determinant.
    #[must_use]
    pub fn determinant(&self) -> f64 {
        let m = &self.rows;
        m[0][0].mul_add(m[1][1], -(m[0][1] * m[1][0]))
    }

    /// The inverse.
    ///
    /// # Errors
    ///
    /// [`OgError::Numeric`](og_core::OgError::Numeric) if the matrix is
    /// singular. See [`Matrix3::inverse`] for why this takes no tolerance.
    pub fn inverse(&self) -> OgResult<Self> {
        let d = self.determinant();
        let scale = self
            .rows
            .iter()
            .flatten()
            .fold(0.0_f64, |acc, v| acc.max(v.abs()));
        if d.abs() <= 4.0 * f64::EPSILON * scale * scale {
            og_bail!(Numeric, "matrix is singular (determinant {d})");
        }
        let m = &self.rows;
        Ok(Self::new([
            [m[1][1] / d, -m[0][1] / d],
            [-m[1][0] / d, m[0][0] / d],
        ]))
    }

    /// Whether every element agrees with `other` within `eps`.
    #[must_use]
    pub fn is_equal(&self, other: &Self, eps: f64) -> bool {
        self.rows
            .iter()
            .flatten()
            .zip(other.rows.iter().flatten())
            .all(|(a, b)| (a - b).abs() <= eps)
    }
}

impl Mul<Vector> for Matrix3 {
    type Output = Vector;
    fn mul(self, v: Vector) -> Vector {
        let m = &self.rows;
        Vector::new(
            m[0][0].mul_add(v.x, m[0][1].mul_add(v.y, m[0][2] * v.z)),
            m[1][0].mul_add(v.x, m[1][1].mul_add(v.y, m[1][2] * v.z)),
            m[2][0].mul_add(v.x, m[2][1].mul_add(v.y, m[2][2] * v.z)),
        )
    }
}

impl Mul for Matrix3 {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        let mut out = [[0.0_f64; 3]; 3];
        for (i, row) in out.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = self.rows[i][0].mul_add(
                    o.rows[0][j],
                    self.rows[i][1].mul_add(o.rows[1][j], self.rows[i][2] * o.rows[2][j]),
                );
            }
        }
        Self::new(out)
    }
}

impl Mul<f64> for Matrix3 {
    type Output = Self;
    fn mul(self, s: f64) -> Self {
        let mut out = self.rows;
        for cell in out.iter_mut().flatten() {
            *cell *= s;
        }
        Self::new(out)
    }
}

impl Add for Matrix3 {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        let mut out = self.rows;
        for (i, row) in out.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell += o.rows[i][j];
            }
        }
        Self::new(out)
    }
}

impl Sub for Matrix3 {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        self + (-o)
    }
}

impl Neg for Matrix3 {
    type Output = Self;
    fn neg(self) -> Self {
        self * -1.0
    }
}

impl Mul<Vector2> for Matrix2 {
    type Output = Vector2;
    fn mul(self, v: Vector2) -> Vector2 {
        let m = &self.rows;
        Vector2::new(
            m[0][0].mul_add(v.x, m[0][1] * v.y),
            m[1][0].mul_add(v.x, m[1][1] * v.y),
        )
    }
}

impl Mul for Matrix2 {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        let mut out = [[0.0_f64; 2]; 2];
        for (i, row) in out.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = self.rows[i][0].mul_add(o.rows[0][j], self.rows[i][1] * o.rows[1][j]);
            }
        }
        Self::new(out)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use og_core::Tolerances;

    /// Only needed where a [`Direction`] has to be constructed.
    const T: Tolerances = Tolerances::millimetres();
    /// Matrix entries are dimensionless, so comparisons take a bare epsilon.
    /// A few hundred ulps covers accumulated rounding in a product of rotations.
    const EPS: f64 = 1e-13;

    #[test]
    fn identity_is_neutral() {
        let v = Vector::new(1.0, 2.0, 3.0);
        assert_eq!(Matrix3::IDENTITY * v, v);
        let m = Matrix3::rotation(Direction::Z, 0.7);
        assert!((m * Matrix3::IDENTITY).is_equal(&m, EPS));
        assert!((Matrix3::IDENTITY * m).is_equal(&m, EPS));
    }

    #[test]
    fn rotation_is_orthonormal_and_preserves_length() {
        let axis = Direction::from_coords(1.0, 2.0, 3.0, T).unwrap();
        for k in 0..8 {
            let m = Matrix3::rotation(axis, f64::from(k) * 0.7);
            assert!(m.is_orthonormal(EPS));
            assert_relative_eq!(m.determinant(), 1.0, epsilon = 1e-14);
            let v = Vector::new(3.0, -1.0, 2.0);
            assert_relative_eq!((m * v).magnitude(), v.magnitude(), epsilon = 1e-13);
        }
    }

    #[test]
    fn rotation_about_z_matches_the_hand_computation() {
        let m = Matrix3::rotation(Direction::Z, core::f64::consts::FRAC_PI_2);
        let v = m * Vector::X;
        assert_relative_eq!(v.x, 0.0, epsilon = 1e-15);
        assert_relative_eq!(v.y, 1.0, epsilon = 1e-15);
        assert_relative_eq!(v.z, 0.0, epsilon = 1e-15);
    }

    #[test]
    fn rotation_composes_additively_in_angle() {
        let axis = Direction::from_coords(0.0, 1.0, 1.0, T).unwrap();
        let a = Matrix3::rotation(axis, 0.3);
        let b = Matrix3::rotation(axis, 0.4);
        let ab = Matrix3::rotation(axis, 0.7);
        assert!((a * b).is_equal(&ab, EPS));
    }

    #[test]
    fn reflection_is_an_involution_with_negative_determinant() {
        let n = Direction::from_coords(1.0, 1.0, 0.0, T).unwrap();
        let m = Matrix3::reflection(n);
        assert_relative_eq!(m.determinant(), -1.0, epsilon = 1e-14);
        assert!((m * m).is_equal(&Matrix3::IDENTITY, EPS));
        // A vector in the mirror plane is unmoved; the normal is negated.
        let in_plane = Vector::new(1.0, -1.0, 0.0);
        assert!((m * in_plane).is_equal(in_plane, T));
        assert!((m * n.vector()).is_equal(-n.vector(), T));
    }

    #[test]
    fn inverse_round_trips() {
        let m = Matrix3::new([[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 4.0]]);
        let inv = m.inverse().unwrap();
        assert!((m * inv).is_equal(&Matrix3::IDENTITY, EPS));
        assert!((inv * m).is_equal(&Matrix3::IDENTITY, EPS));
    }

    #[test]
    fn singular_matrices_are_refused() {
        // Third row is the sum of the first two.
        let m = Matrix3::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [5.0, 7.0, 9.0]]);
        assert!(m.inverse().is_err());
        assert!(Matrix3::ZERO.inverse().is_err());
    }

    #[test]
    fn singularity_threshold_is_relative_not_absolute() {
        // A well-conditioned matrix scaled down by 1e-3 has a determinant of
        // 1e-9 — below an absolute confusion threshold, but perfectly
        // invertible. An absolute test would reject it.
        let m = Matrix3::IDENTITY * 1e-3;
        assert_relative_eq!(m.determinant(), 1e-9, max_relative = 1e-12);
        let inv = m.inverse().unwrap();
        assert!((m * inv).is_equal(&Matrix3::IDENTITY, EPS));
        // And the same matrix scaled up stays invertible too.
        let big = Matrix3::IDENTITY * 1e6;
        assert!((big * big.inverse().unwrap()).is_equal(&Matrix3::IDENTITY, EPS));
    }

    #[test]
    fn orthonormal_inverse_equals_transpose() {
        let axis = Direction::from_coords(2.0, -1.0, 0.5, T).unwrap();
        let m = Matrix3::rotation(axis, 1.1);
        assert!(m.inverse().unwrap().is_equal(&m.transposed(), EPS));
        assert!(!Matrix3::scaling(2.0).is_orthonormal(EPS));
    }

    #[test]
    fn determinant_and_trace() {
        let m = Matrix3::scaling_xyz(2.0, 3.0, 4.0);
        assert_relative_eq!(m.determinant(), 24.0);
        assert_relative_eq!(m.trace(), 9.0);
    }

    #[test]
    fn columns_and_rows_are_bounds_checked() {
        let m = Matrix3::from_columns(Vector::X, Vector::Y, Vector::Z);
        assert_eq!(m.column(0).unwrap(), Vector::X);
        assert_eq!(m.row(1).unwrap(), Vector::Y);
        assert!(m.column(3).is_err());
        assert!(m.row(3).is_err());
        assert!(m.get(0, 3).is_err());
        assert!(m.is_equal(&Matrix3::IDENTITY, EPS));
    }

    #[test]
    fn matrix2_rotation_and_inverse() {
        let m = Matrix2::rotation(core::f64::consts::FRAC_PI_2);
        let v = m * Vector2::X;
        assert_relative_eq!(v.x, 0.0, epsilon = 1e-15);
        assert_relative_eq!(v.y, 1.0, epsilon = 1e-15);
        assert_relative_eq!(m.determinant(), 1.0, epsilon = 1e-15);
        assert!((m * m.inverse().unwrap()).is_equal(&Matrix2::IDENTITY, EPS));
        assert!(Matrix2::ZERO.inverse().is_err());
    }
}
