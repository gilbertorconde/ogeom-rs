//! Property tests for the analytic primitives.
//!
//! The laws here are the defining ones — the properties that make a shape that
//! shape, checked independently of how it is stored: the focal sum of an
//! ellipse, the constant distance of a sphere, the invariance of every shape
//! under a similarity.

#![allow(clippy::unwrap_used, clippy::expect_used)]

use ogeom_core::Tolerances;
use ogeom_math::{
    Axis, Circle, Cone, Cylinder, Direction, Ellipse, Frame, Hyperbola, Plane, Point, Sphere,
    Torus, TorusKind, Transform, Vector, conic::complete_elliptic_e,
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

fn frame() -> impl Strategy<Value = Frame> {
    (point(), direction()).prop_map(|(p, d)| Frame::about(p, d))
}

/// A radius large enough that millimetre-scale tolerances are meaningful.
fn radius() -> impl Strategy<Value = f64> {
    0.01f64..100.0
}

/// Similarities only — the transforms under which analytic shapes survive.
fn similarity() -> impl Strategy<Value = Transform> {
    prop_oneof![
        vector().prop_map(Transform::translation),
        (point(), direction(), -10.0f64..10.0)
            .prop_map(|(p, d, a)| Transform::rotation(Axis::new(p, d), a)),
        (point(), 0.05f64..20.0)
            .prop_filter_map("scale", |(p, s)| Transform::scaling(p, s, T).ok()),
        (point(), direction()).prop_map(|(p, d)| Transform::plane_mirror(p, d)),
    ]
}

fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * a.abs().max(b.abs()).max(1.0)
}

proptest! {
    /// The focal definition of an ellipse: every point on it has the same sum
    /// of distances to the two foci, namely twice the major radius.
    #[test]
    fn ellipse_focal_sum_is_constant(
        f in frame(), a in radius(), ratio in 0.01f64..1.0, angle in -6.0f64..6.0,
    ) {
        let b = a * ratio;
        let e = Ellipse::new(f, a, b, T).unwrap();
        let (f1, f2) = e.foci();
        // A point on the curve, from the parametric form.
        let p = e.centre() + f.x() * (a * angle.cos()) + f.y() * (b * angle.sin());
        prop_assert!(close(p.distance(f1) + p.distance(f2), 2.0 * a, 1e-9));
    }

    /// And the hyperbolic one: the *difference* of distances is constant.
    #[test]
    fn hyperbola_focal_difference_is_constant(
        f in frame(), a in radius(), b in radius(), t in -2.0f64..2.0,
    ) {
        let h = Hyperbola::new(f, a, b, T).unwrap();
        let (f1, f2) = h.foci();
        let p = h.centre() + f.x() * (a * t.cosh()) + f.y() * (b * t.sinh());
        prop_assert!(close((p.distance(f2) - p.distance(f1)).abs(), 2.0 * a, 1e-8));
    }

    #[test]
    fn circle_points_are_equidistant_from_the_centre(
        f in frame(), r in radius(), angle in -6.0f64..6.0,
    ) {
        let c = Circle::new(f, r, T).unwrap();
        let p = c.centre() + f.x() * (r * angle.cos()) + f.y() * (r * angle.sin());
        prop_assert!(close(c.centre().distance(p), r, 1e-12));
        prop_assert!(c.distance_to(p) <= 1e-9 * r.max(1.0));
    }

    #[test]
    fn circle_through_three_points_passes_through_all_three(
        a in point(), b in point(), c in point(),
    ) {
        let Ok(circle) = Circle::through(a, b, c, T) else { return Ok(()) };
        for p in [a, b, c] {
            prop_assert!(close(circle.centre().distance(p), circle.radius(), 1e-8));
        }
        // And its plane contains them.
        let plane = Plane::new(circle.frame());
        for p in [a, b, c] {
            prop_assert!(plane.distance_to(p) <= 1e-8 * circle.radius().max(1.0));
        }
    }

    #[test]
    fn sphere_surface_points_are_at_exactly_the_radius(
        c in point(), r in radius(), d in direction(),
    ) {
        let s = Sphere::centred(c, r, T).unwrap();
        let on = c + d * r;
        prop_assert!(close(c.distance(on), r, 1e-12));
        prop_assert!(s.distance_to(on) <= 1e-9 * r.max(1.0));
        prop_assert!(s.encloses(c + d * (r * 0.5), T));
        prop_assert!(!s.encloses(c + d * (r * 1.5), T));

        // The normal is recovered by subtracting the centre from a surface
        // point. For a small sphere far from the origin that subtraction
        // cancels: a centre at z = 400 and a radius of 0.01 leaves an absolute
        // error around 4e-14 in a vector of length 0.01, so the recovered
        // direction is off by some 4e-12 radians. `is_equal` compares against
        // `angular` at 1e-12, which is a threshold for *stored* directions, and
        // asserting it here would be asserting something false about floating
        // point rather than something true about spheres.
        let normal = s.normal_at(on, T).unwrap();
        let coordinate_scale = c.to_vector().magnitude().max(r);
        let allowed = 1e-13 * (coordinate_scale / r).max(1.0);
        prop_assert!(
            normal.angle(d) <= allowed,
            "normal off by {} rad, allowed {allowed}",
            normal.angle(d)
        );
    }

    #[test]
    fn cylinder_surface_points_are_at_exactly_the_radius(
        f in frame(), r in radius(), angle in -6.0f64..6.0, h in -100.0f64..100.0,
    ) {
        let c = Cylinder::new(f, r, T).unwrap();
        let p = f.origin() + f.x() * (r * angle.cos()) + f.y() * (r * angle.sin()) + f.z() * h;
        prop_assert!(c.distance_to(p) <= 1e-9 * r.max(h.abs()).max(1.0));
        // The normal is radial, so perpendicular to the axis.
        let n = c.normal_at(p, T).unwrap();
        prop_assert!(n.dot(f.z()).abs() <= 1e-9);
    }

    #[test]
    fn cone_surface_points_lie_on_the_cone(
        f in frame(), r0 in 0.0f64..50.0, half in 0.05f64..1.5,
        z in -50.0f64..50.0, angle in -6.0f64..6.0,
    ) {
        let c = Cone::new(f, r0, half, T).unwrap();
        let r = c.radius_at(z);
        prop_assume!(r >= 0.0);
        let p = f.origin() + f.x() * (r * angle.cos()) + f.y() * (r * angle.sin()) + f.z() * z;
        let scale = r.max(z.abs()).max(1.0);
        prop_assert!(
            c.distance_to(p) <= 1e-9 * scale,
            "distance {} at r = {r}, z = {z}",
            c.distance_to(p)
        );
    }

    #[test]
    fn cone_apex_is_where_the_radius_vanishes(f in frame(), r0 in 0.0f64..50.0, half in 0.05f64..1.5) {
        let c = Cone::new(f, r0, half, T).unwrap();
        let apex = c.apex();
        // The apex is on the axis, and the radius there is zero.
        prop_assert!(c.axis().distance_to(apex) <= 1e-9 * r0.max(1.0));
        let z = f.z().dot_vector(apex - f.origin());
        prop_assert!(c.radius_at(z).abs() <= 1e-9 * r0.max(1.0));
        prop_assert!(c.distance_to(apex) <= 1e-9 * r0.max(1.0));
    }

    #[test]
    fn torus_surface_points_lie_on_the_torus(
        f in frame(), major in radius(), minor in radius(),
        u in -6.0f64..6.0, v in -6.0f64..6.0,
    ) {
        let t = Torus::new(f, major, minor, T).unwrap();
        // The standard parameterization: go out to the tube centre, then round
        // the tube.
        let out = f.x() * u.cos() + f.y() * u.sin();
        let p = f.origin() + out * major.mul_add(1.0, minor * v.cos()) + f.z() * (minor * v.sin());
        prop_assert!(t.distance_to(p) <= 1e-9 * (major + minor));
    }

    #[test]
    fn torus_kind_is_decided_by_the_two_radii(major in radius(), minor in radius()) {
        let t = Torus::new(Frame::WORLD, major, minor, T).unwrap();
        let expected = if (major - minor).abs() <= T.confusion() {
            TorusKind::Horn
        } else if major > minor {
            TorusKind::Ring
        } else {
            TorusKind::Spindle
        };
        prop_assert_eq!(t.kind(T), expected);
    }

    #[test]
    fn plane_projection_is_idempotent_and_lands_on_the_plane(
        f in frame(), p in point(),
    ) {
        let plane = Plane::new(f);
        let projected = plane.project(p);
        prop_assert!(plane.distance_to(projected) <= 1e-9 * p.to_vector().magnitude().max(1.0));
        prop_assert!(plane.project(projected).is_equal(projected, T));
        // The displacement is along the normal. Stated as a relative bound on
        // the cross product rather than through `is_collinear`, whose angular
        // tolerance of 1e-12 is a comparison for exact directions: the offset
        // here is a difference of coordinates that may be a thousand units
        // apart, so its own rounding already exceeds that.
        let offset = p - projected;
        let cross = offset.cross(plane.normal().vector()).magnitude();
        prop_assert!(cross <= 1e-9 * offset.magnitude().max(1.0));
    }

    /// The reason `Transform` is restricted to similarities: an analytic shape
    /// stays the same kind of shape, and a point on it stays on it.
    #[test]
    fn similarities_map_spheres_to_spheres(
        c in point(), r in radius(), d in direction(), t in similarity(),
    ) {
        let s = Sphere::centred(c, r, T).unwrap();
        let moved = s.transformed(&t, T).unwrap();
        prop_assert!(close(moved.radius(), r * t.scale_factor().abs(), 1e-10));
        let on = c + d * r;
        let scale = moved.radius().max(moved.centre().to_vector().magnitude()).max(1.0);
        prop_assert!(moved.distance_to(t.apply(on)) <= 1e-9 * scale);
    }

    #[test]
    fn similarities_map_cylinders_to_cylinders(
        f in frame(), r in radius(), angle in -6.0f64..6.0, h in -50.0f64..50.0,
        t in similarity(),
    ) {
        let c = Cylinder::new(f, r, T).unwrap();
        let moved = c.transformed(&t, T).unwrap();
        prop_assert!(close(moved.radius(), r * t.scale_factor().abs(), 1e-10));
        let p = f.origin() + f.x() * (r * angle.cos()) + f.y() * (r * angle.sin()) + f.z() * h;
        let scale = moved.radius().max(t.apply(p).to_vector().magnitude()).max(1.0);
        prop_assert!(moved.distance_to(t.apply(p)) <= 1e-8 * scale);
    }

    #[test]
    fn similarities_preserve_a_cones_half_angle(
        f in frame(), r0 in 0.0f64..50.0, half in 0.05f64..1.5, t in similarity(),
    ) {
        let c = Cone::new(f, r0, half, T).unwrap();
        let moved = c.transformed(&t, T).unwrap();
        prop_assert_eq!(moved.half_angle(), half);
        prop_assert!(close(moved.reference_radius(), r0 * t.scale_factor().abs(), 1e-10));
    }

    #[test]
    fn similarities_scale_a_planes_distances_uniformly(
        f in frame(), p in point(), t in similarity(),
    ) {
        let plane = Plane::new(f);
        let moved = plane.transformed(&t, T).unwrap();
        let before = plane.distance_to(p);
        let after = moved.distance_to(t.apply(p));
        let scale = before.max(after).max(1.0) * t.scale_factor().abs().max(1.0);
        prop_assert!((after - before * t.scale_factor().abs()).abs() <= 1e-8 * scale);
    }

    #[test]
    fn radii_are_never_negative_after_any_similarity(
        f in frame(), r in radius(), t in similarity(),
    ) {
        // A mirror has a negative scale factor; the radius must not follow it.
        prop_assert!(Circle::new(f, r, T).unwrap().transformed(&t, T).unwrap().radius() > 0.0);
        prop_assert!(Cylinder::new(f, r, T).unwrap().transformed(&t, T).unwrap().radius() > 0.0);
        prop_assert!(Sphere::new(f, r, T).unwrap().transformed(&t, T).unwrap().radius() > 0.0);
    }

    #[test]
    fn complete_elliptic_integral_is_bounded_and_monotone(m in 0.0f64..1.0) {
        let e = complete_elliptic_e(m);
        prop_assert!((1.0..=core::f64::consts::FRAC_PI_2 + 1e-15).contains(&e));
        // Decreasing in m.
        if m < 0.99 {
            prop_assert!(complete_elliptic_e(m + 0.01) < e);
        }
    }

    /// The perimeter must sit between the two elementary bounds it is squeezed
    /// by: no shorter than the inscribed circle, no longer than the
    /// circumscribed one.
    #[test]
    fn ellipse_perimeter_is_bracketed(a in radius(), ratio in 0.01f64..1.0) {
        let b = a * ratio;
        let e = Ellipse::new(Frame::WORLD, a, b, T).unwrap();
        let length = e.length();
        prop_assert!(length >= core::f64::consts::TAU * b - 1e-9);
        prop_assert!(length <= core::f64::consts::TAU * a + 1e-9);
        // And it exceeds the perimeter of the inscribed rhombus.
        prop_assert!(length >= 4.0 * a.hypot(b) - 1e-9);
    }
}
