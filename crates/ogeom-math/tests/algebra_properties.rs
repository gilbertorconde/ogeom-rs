//! Property tests for the core algebra.
//!
//! Identities rather than examples: Lagrange's identity, anticommutativity of
//! the cross product, invariance of length and angle under rotation, quaternion
//! and matrix agreement. These are the facts every algorithm downstream will
//! silently assume.

#![allow(clippy::unwrap_used, clippy::expect_used)]

use ogeom_core::Tolerances;
use ogeom_math::{Direction, Matrix3, Point, Quaternion, Vector, Vector2};
use proptest::prelude::*;

const T: Tolerances = Tolerances::millimetres();

/// Coordinates spanning ten orders of magnitude, so the identities are checked
/// where cancellation actually bites and not only near unity.
fn coord() -> impl Strategy<Value = f64> {
    prop_oneof![-1.0f64..1.0, -1e3f64..1e3, -1e-3f64..1e-3]
}

fn vector() -> impl Strategy<Value = Vector> {
    (coord(), coord(), coord()).prop_map(|(x, y, z)| Vector::new(x, y, z))
}

fn vector2() -> impl Strategy<Value = Vector2> {
    (coord(), coord()).prop_map(|(x, y)| Vector2::new(x, y))
}

fn point() -> impl Strategy<Value = Point> {
    (coord(), coord(), coord()).prop_map(|(x, y, z)| Point::new(x, y, z))
}

/// A direction, by rejecting vectors too short to normalize.
fn direction() -> impl Strategy<Value = Direction> {
    vector().prop_filter_map("degenerate", |v| Direction::new(v, T).ok())
}

/// Relative comparison, so the assertion means the same thing at every scale.
fn close(a: f64, b: f64, rel: f64) -> bool {
    let scale = a.abs().max(b.abs()).max(1.0);
    (a - b).abs() <= rel * scale
}

proptest! {
    #[test]
    fn cross_is_anticommutative_and_orthogonal(a in vector(), b in vector()) {
        let c = a.cross(b);
        prop_assert!(c.is_equal(-b.cross(a), T));
        // Orthogonality holds to the precision of the products involved, which
        // is set by the magnitudes, not by an absolute distance.
        let scale = a.magnitude() * b.magnitude();
        prop_assert!(c.dot(a).abs() <= 1e-12 * scale * a.magnitude().max(1.0));
        prop_assert!(c.dot(b).abs() <= 1e-12 * scale * b.magnitude().max(1.0));
    }

    /// Lagrange's identity: |a x b|^2 + (a . b)^2 = |a|^2 |b|^2.
    #[test]
    fn lagrange_identity(a in vector(), b in vector()) {
        let lhs = a.cross(b).square_magnitude() + a.dot(b) * a.dot(b);
        let rhs = a.square_magnitude() * b.square_magnitude();
        prop_assert!(close(lhs, rhs, 1e-12), "{lhs} != {rhs}");
    }

    #[test]
    fn cross_of_a_vector_with_itself_is_exactly_zero(a in vector()) {
        // Exactly, not approximately: each component is `xy - yx` with both
        // products rounding identically. This is the property the deliberate
        // absence of a fused multiply-add in `cross` buys.
        prop_assert_eq!(a.cross(a), Vector::ZERO);
    }

    #[test]
    fn triple_product_is_cyclic(a in vector(), b in vector(), c in vector()) {
        let scale = (a.magnitude() * b.magnitude() * c.magnitude()).max(1.0);
        prop_assert!((a.triple(b, c) - b.triple(c, a)).abs() <= 1e-11 * scale);
        prop_assert!((a.triple(b, c) + a.triple(c, b)).abs() <= 1e-11 * scale);
    }

    #[test]
    fn normalization_is_idempotent_and_yields_unit_length(v in vector()) {
        if let Ok(n) = v.normalized(T) {
            prop_assert!(close(n.magnitude(), 1.0, 1e-14));
            let again = n.normalized(T).unwrap();
            prop_assert!(again.is_equal(n, T));
        }
    }

    #[test]
    fn rotation_preserves_length_and_angle(
        axis in direction(),
        angle in -10.0f64..10.0,
        a in vector(),
        b in vector(),
    ) {
        let m = Matrix3::rotation(axis, angle);
        let (ra, rb) = (m * a, m * b);
        prop_assert!(close(ra.magnitude(), a.magnitude(), 1e-12));
        prop_assert!(close(rb.magnitude(), b.magnitude(), 1e-12));
        let scale = a.magnitude() * b.magnitude();
        prop_assert!((ra.dot(rb) - a.dot(b)).abs() <= 1e-11 * scale.max(1.0));
        prop_assert!(m.is_orthonormal(1e-13));
    }

    #[test]
    fn rotation_by_an_angle_and_its_negation_cancel(
        axis in direction(),
        angle in -10.0f64..10.0,
        v in vector(),
    ) {
        let forward = Matrix3::rotation(axis, angle);
        let back = Matrix3::rotation(axis, -angle);
        prop_assert!((back * (forward * v)).is_equal(v, T));
    }

    #[test]
    fn a_rotation_leaves_its_own_axis_fixed(axis in direction(), angle in -10.0f64..10.0) {
        let m = Matrix3::rotation(axis, angle);
        prop_assert!((m * axis.vector()).is_equal(axis.vector(), T));
    }

    #[test]
    fn quaternion_and_matrix_rotations_agree(
        axis in direction(),
        angle in -10.0f64..10.0,
        v in vector(),
    ) {
        let q = Quaternion::from_axis_angle(axis, angle);
        prop_assert!(q.rotate(v).is_equal(q.to_matrix() * v, T));
    }

    #[test]
    fn quaternion_survives_a_matrix_round_trip(axis in direction(), angle in -3.0f64..3.0) {
        let q = Quaternion::from_axis_angle(axis, angle);
        let back = Quaternion::from_matrix(&q.to_matrix(), T).unwrap();
        // q and -q name the same rotation, so compare what they do, not what
        // they are.
        prop_assert!(back.to_matrix().is_equal(&q.to_matrix(), 1e-11));
    }

    #[test]
    fn quaternion_composition_matches_matrix_composition(
        a1 in direction(), t1 in -3.0f64..3.0,
        a2 in direction(), t2 in -3.0f64..3.0,
    ) {
        let (p, q) = (
            Quaternion::from_axis_angle(a1, t1),
            Quaternion::from_axis_angle(a2, t2),
        );
        prop_assert!((p * q).to_matrix().is_equal(&(p.to_matrix() * q.to_matrix()), 1e-11));
    }

    #[test]
    fn conjugation_undoes_rotation(axis in direction(), angle in -10.0f64..10.0, v in vector()) {
        let q = Quaternion::from_axis_angle(axis, angle);
        prop_assert!(q.conjugate().rotate(q.rotate(v)).is_equal(v, T));
    }

    #[test]
    fn matrix_inverse_round_trips_when_invertible(
        a in vector(), b in vector(), c in vector(), v in vector(),
    ) {
        let m = Matrix3::from_columns(a, b, c);
        let Ok(inv) = m.inverse() else { return Ok(()) };

        // How much accuracy survives a round trip is set by the conditioning of
        // the matrix, not by any fixed distance. Columns spanning six orders of
        // magnitude give a condition number around 1e6, and the result is then
        // accurate to about 1e-16 * 1e6 relative — which can exceed an absolute
        // millimetre-scale tolerance while being perfectly correct. Asserting
        // `is_equal` here would be asserting something false about floating
        // point, not something true about the inverse.
        let norm = |x: &Matrix3| {
            x.rows.iter().map(|r| r.iter().map(|e| e.abs()).sum::<f64>()).fold(0.0, f64::max)
        };
        let condition = norm(&m) * norm(&inv);
        let round_trip = inv * (m * v);
        let error = (round_trip - v).magnitude();
        let allowed = 1e-13 * condition * v.magnitude().max(1.0);

        prop_assert!(
            error <= allowed,
            "error {error} exceeds {allowed} at condition {condition}"
        );
    }

    #[test]
    fn point_and_vector_form_an_affine_space(p in point(), q in point(), v in vector()) {
        prop_assert!(((p + v) - p).is_equal(v, T));
        prop_assert!((p + (q - p)).is_equal(q, T));
        prop_assert!((q - p).is_equal(-(p - q), T));
        prop_assert!(close(p.distance(q), q.distance(p), 1e-15));
    }

    #[test]
    fn midpoint_is_equidistant(p in point(), q in point()) {
        let m = p.midpoint(q);
        prop_assert!(close(m.distance(p), m.distance(q), 1e-12));
        prop_assert!(close(m.distance(p) * 2.0, p.distance(q), 1e-12));
    }

    #[test]
    fn centroid_of_a_translated_set_translates_with_it(
        ps in prop::collection::vec(point(), 1..12),
        v in vector(),
    ) {
        let c = Point::centroid(&ps).unwrap();
        let moved: Vec<_> = ps.iter().map(|p| *p + v).collect();
        prop_assert!(Point::centroid(&moved).unwrap().is_equal(c + v, T));
    }

    #[test]
    fn vector2_perpendicular_is_exact(v in vector2()) {
        // Exact zero, again by construction: the two products are the same
        // magnitudes with opposite signs.
        prop_assert_eq!(v.perpendicular().dot(v), 0.0);
        prop_assert_eq!(v.perpendicular().perpendicular(), -v);
        prop_assert!(close(v.perpendicular().magnitude(), v.magnitude(), 1e-15));
    }

    #[test]
    fn vector2_cross_is_antisymmetric(a in vector2(), b in vector2()) {
        prop_assert_eq!(a.cross(b), -b.cross(a));
        prop_assert_eq!(a.cross(a), 0.0);
    }

    #[test]
    fn any_perpendicular_is_always_unit_and_orthogonal(d in direction()) {
        let p = d.any_perpendicular();
        prop_assert!(close(p.vector().magnitude(), 1.0, 1e-14));
        prop_assert!(d.dot(p).abs() <= 1e-14);
    }

    #[test]
    fn direction_angle_is_symmetric_and_bounded(a in direction(), b in direction()) {
        prop_assert!(close(a.angle(b), b.angle(a), 1e-14));
        prop_assert!((0.0..=core::f64::consts::PI + 1e-12).contains(&a.angle(b)));
        prop_assert!(a.angle(a) <= 1e-12);
    }
}
