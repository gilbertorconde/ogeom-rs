//! Property tests for frames and transforms.
//!
//! The laws that matter here: composition is associative, inversion undoes
//! application, a similarity preserves ratios of distances and all angles, and
//! the classification used for dispatch never disagrees with the general
//! formula it replaces.

#![allow(clippy::unwrap_used, clippy::expect_used)]

use ogeom_core::Tolerances;
use ogeom_math::{
    Axis, Direction, Frame, GeneralTransform, Handedness, Matrix3, Point, Transform, TransformKind,
    Vector,
};
use proptest::prelude::*;

const T: Tolerances = Tolerances::millimetres();

fn coord() -> impl Strategy<Value = f64> {
    prop_oneof![-10.0f64..10.0, -1e3f64..1e3]
}

fn point() -> impl Strategy<Value = Point> {
    (coord(), coord(), coord()).prop_map(|(x, y, z)| Point::new(x, y, z))
}

fn vector() -> impl Strategy<Value = Vector> {
    (coord(), coord(), coord()).prop_map(|(x, y, z)| Vector::new(x, y, z))
}

fn direction() -> impl Strategy<Value = Direction> {
    vector().prop_filter_map("degenerate", |v| Direction::new(v, T).ok())
}

fn axis() -> impl Strategy<Value = Axis> {
    (point(), direction()).prop_map(|(p, d)| Axis::new(p, d))
}

fn frame() -> impl Strategy<Value = Frame> {
    (point(), direction()).prop_map(|(p, d)| Frame::about(p, d))
}

/// Every kind of similarity, so the laws are checked across all dispatch paths.
fn transform() -> impl Strategy<Value = Transform> {
    prop_oneof![
        Just(Transform::IDENTITY),
        vector().prop_map(Transform::translation),
        (axis(), -10.0f64..10.0).prop_map(|(a, t)| Transform::rotation(a, t)),
        point().prop_map(Transform::point_mirror),
        (point(), direction()).prop_map(|(p, d)| Transform::plane_mirror(p, d)),
        (point(), prop_oneof![-10.0f64..-0.1, 0.1f64..10.0])
            .prop_filter_map("degenerate scale", |(p, s)| Transform::scaling(p, s, T)
                .ok()),
        (axis(), -3.0f64..3.0, point(), 0.2f64..5.0).prop_filter_map("compound", |(a, t, p, s)| {
            Some(Transform::rotation(a, t) * Transform::scaling(p, s, T).ok()?)
        }),
    ]
}

/// Relative comparison of points, so the assertion means the same thing at
/// every distance from the origin. Transforms with a large translation move
/// coordinates far enough that an absolute millimetre tolerance stops being the
/// right question.
fn points_close(a: Point, b: Point, rel: f64) -> bool {
    let scale = a
        .to_vector()
        .magnitude()
        .max(b.to_vector().magnitude())
        .max(1.0);
    a.distance(b) <= rel * scale
}

proptest! {
    #[test]
    fn classification_never_lies(t in transform()) {
        // Whatever fast path `apply` takes must land in the same place as the
        // general formula. This is the property that makes classification a
        // pure optimisation rather than a second implementation.
        let general = |p: Point| {
            Point::from_vector(
                t.linear() * (p.to_vector() * t.scale_factor()) + t.translation_vector(),
            )
        };
        for p in [Point::ORIGIN, Point::new(1.0, -2.0, 3.5), Point::new(1e3, 1e3, -1e3)] {
            prop_assert!(points_close(t.apply(p), general(p), 1e-12));
        }
    }

    #[test]
    fn inversion_undoes_application(t in transform(), p in point()) {
        let inv = t.inverse().unwrap();
        prop_assert!(points_close(inv.apply(t.apply(p)), p, 1e-11));
        prop_assert!(points_close(t.apply(inv.apply(p)), p, 1e-11));
    }

    #[test]
    fn composition_is_associative(a in transform(), b in transform(), c in transform(), p in point()) {
        let left = ((a * b) * c).apply(p);
        let right = (a * (b * c)).apply(p);
        prop_assert!(points_close(left, right, 1e-11));
    }

    #[test]
    fn composition_applies_right_to_left(a in transform(), b in transform(), p in point()) {
        prop_assert!(points_close((a * b).apply(p), a.apply(b.apply(p)), 1e-11));
    }

    #[test]
    fn identity_is_neutral(t in transform(), p in point()) {
        prop_assert!(points_close((t * Transform::IDENTITY).apply(p), t.apply(p), 1e-12));
        prop_assert!(points_close((Transform::IDENTITY * t).apply(p), t.apply(p), 1e-12));
    }

    /// The defining property of a similarity: all distances scale by the same
    /// factor. If this held only approximately, analytic surfaces would not
    /// survive being transformed.
    #[test]
    fn a_similarity_scales_every_distance_equally(
        t in transform(), a in point(), b in point(), c in point(),
    ) {
        let ratio = |p: Point, q: Point| {
            let before = p.distance(q);
            let after = t.apply(p).distance(t.apply(q));
            (before, after)
        };
        let (b1, a1) = ratio(a, b);
        let (b2, a2) = ratio(b, c);
        prop_assume!(b1 > 1e-6 && b2 > 1e-6);
        // a1/b1 must equal a2/b2, tested without dividing.
        let scale = (a1 * b2).abs().max(a2 * b1).max(1.0);
        prop_assert!((a1 * b2 - a2 * b1).abs() <= 1e-9 * scale);
    }

    #[test]
    fn a_similarity_preserves_angles(t in transform(), a in vector(), b in vector()) {
        let (ta, tb) = (t.apply_vector(a), t.apply_vector(b));
        let before = a.angle(b, T);
        let after = ta.angle(tb, T);
        if let (Ok(before), Ok(after)) = (before, after) {
            // A handedness-reversing transform maps an angle to itself, since
            // the unsigned angle is what is preserved.
            prop_assert!((before - after).abs() <= 1e-9, "{before} vs {after}");
        }
    }

    #[test]
    fn handedness_composes_multiplicatively(a in transform(), b in transform()) {
        prop_assert_eq!(
            (a * b).preserves_handedness(),
            a.preserves_handedness() == b.preserves_handedness()
        );
    }

    #[test]
    fn directions_stay_unit_under_any_similarity(t in transform(), d in direction()) {
        let out = t.apply_direction(d, T).unwrap();
        prop_assert!((out.vector().magnitude() - 1.0).abs() <= 1e-12);
    }

    #[test]
    fn translation_never_moves_a_free_vector(v in vector(), w in vector()) {
        prop_assert!(Transform::translation(w).apply_vector(v).is_equal(v, T));
    }

    #[test]
    fn a_rotation_leaves_its_axis_pointwise_fixed(a in axis(), angle in -10.0f64..10.0, t in -1e3f64..1e3) {
        let r = Transform::rotation(a, angle);
        let on_axis = a.point_at(t);
        prop_assert!(points_close(r.apply(on_axis), on_axis, 1e-11));
    }

    #[test]
    fn a_mirror_is_its_own_inverse(p in point(), n in direction(), q in point()) {
        let m = Transform::plane_mirror(p, n);
        prop_assert!(points_close((m * m).apply(q), q, 1e-11));
        prop_assert!(!m.preserves_handedness());
    }

    #[test]
    fn scaling_multiplies_distances_by_its_factor(
        c in point(), f in prop_oneof![-10.0f64..-0.1, 0.1f64..10.0], p in point(),
    ) {
        let t = Transform::scaling(c, f, T).unwrap();
        let before = c.distance(p);
        let after = c.distance(t.apply(p));
        let scale = before.max(1.0) * f.abs().max(1.0);
        prop_assert!((after - before * f.abs()).abs() <= 1e-10 * scale);
    }

    #[test]
    fn frames_round_trip_between_local_and_world(f in frame(), p in point()) {
        prop_assert!(points_close(f.to_world(f.to_local(p)), p, 1e-11));
        prop_assert!(points_close(f.to_local(f.to_world(p)), p, 1e-11));
    }

    #[test]
    fn a_frame_is_always_orthonormal_and_right_handed(f in frame()) {
        prop_assert!(f.to_matrix().is_orthonormal(1e-13));
        prop_assert_eq!(f.handedness(), Handedness::Right);
        let triple = f.x().vector().triple(f.y().vector(), f.z().vector());
        prop_assert!((triple - 1.0).abs() <= 1e-12);
    }

    #[test]
    fn frame_transforms_agree_with_frame_conversion(f in frame(), p in point()) {
        prop_assert!(points_close(Transform::to_frame(&f).apply(p), f.to_local(p), 1e-11));
        prop_assert!(points_close(Transform::from_frame(&f).apply(p), f.to_world(p), 1e-11));
    }

    #[test]
    fn between_frames_preserves_world_position(a in frame(), b in frame(), local in point()) {
        let t = Transform::between_frames(&a, &b);
        prop_assert!(points_close(b.to_world(t.apply(local)), a.to_world(local), 1e-10));
    }

    #[test]
    fn reversing_a_frames_primary_direction_preserves_handedness(f in frame()) {
        let r = f.with_z_reversed();
        prop_assert_eq!(r.handedness(), Handedness::Right);
        prop_assert!(r.z().is_equal(-f.z(), T));
        prop_assert!(r.to_matrix().is_orthonormal(1e-13));
    }

    #[test]
    fn axis_projection_is_idempotent_and_minimal(a in axis(), p in point(), t in -100.0f64..100.0) {
        let proj = a.project(p);
        prop_assert!(points_close(a.project(proj), proj, 1e-10));
        // No other point on the axis is closer.
        prop_assert!(p.distance(proj) <= p.distance(a.point_at(a.parameter_of(p) + t)) + 1e-9);
    }

    #[test]
    fn general_transform_normals_stay_perpendicular_to_the_surface(
        sx in 0.1f64..10.0, sy in 0.1f64..10.0, sz in 0.1f64..10.0,
        n in vector(), tangent in vector(),
    ) {
        // Take any vector perpendicular to `n`; after transforming, the
        // transformed normal must still be perpendicular to it. That is the
        // whole reason normals use the inverse transpose.
        prop_assume!(n.magnitude() > 1e-3 && tangent.magnitude() > 1e-3);
        let t = tangent - n * (tangent.dot(n) / n.square_magnitude());
        prop_assume!(t.magnitude() > 1e-3);

        let g = GeneralTransform::scaling_xyz(sx, sy, sz);
        let tn = g.apply_normal(n).unwrap();
        let tt = g.apply_vector(t);
        let scale = tn.magnitude() * tt.magnitude();
        prop_assert!(tn.dot(tt).abs() <= 1e-10 * scale.max(1.0));
    }

    #[test]
    fn general_transform_recognizes_exactly_the_similarities(t in transform()) {
        let g: GeneralTransform = t.into();
        let narrowed = g.to_similarity(1e-10);
        prop_assert!(narrowed.is_some(), "a similarity must be recognized as one");
        let n = narrowed.unwrap();
        for p in [Point::ORIGIN, Point::new(2.0, -1.0, 0.5)] {
            prop_assert!(points_close(n.apply(p), t.apply(p), 1e-9));
        }
    }

    #[test]
    fn non_uniform_scaling_is_not_a_similarity(sx in 0.1f64..10.0, sy in 0.1f64..10.0) {
        prop_assume!((sx - sy).abs() > 0.1);
        prop_assert!(!GeneralTransform::scaling_xyz(sx, sy, 1.0).is_similarity(1e-10));
    }

    #[test]
    fn a_transform_built_from_a_rotation_matrix_is_classified_as_one(
        d in direction(), angle in 0.1f64..3.0,
    ) {
        let t = Transform::from_quaternion(ogeom_math::Quaternion::from_axis_angle(d, angle));
        prop_assert_eq!(t.kind(), TransformKind::Rotation);
        prop_assert!(t.linear().is_equal(&Matrix3::rotation(d, angle), 1e-13));
    }
}
