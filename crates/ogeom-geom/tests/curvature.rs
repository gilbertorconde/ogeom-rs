//! Principal curvatures where the closed forms are known.
#![allow(clippy::unwrap_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::{CylinderSurface, PlaneSurface, SphereSurface, Surface, TorusSurface};
use ogeom_math::{Cylinder, Frame, Plane, Sphere, Torus};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_plane_has_no_curvature_anywhere() {
    let plane = PlaneSurface::new(Plane::new(Frame::WORLD));
    let c = plane.curvature_at(3.0, -2.0, T).unwrap();
    assert!(c.max.abs() < 1e-12 && c.min.abs() < 1e-12, "{c:?}");
    assert!(c.is_umbilic(T));
    assert!(c.max_direction.vector().dot(c.min_direction.vector()).abs() < 1e-12);
}

#[test]
fn a_sphere_curves_the_same_every_way() {
    let sphere = SphereSurface::new(Sphere::new(Frame::WORLD, 2.5, T).unwrap());
    for (u, v) in [(0.3, 0.2), (2.0, -0.9), (4.5, 1.1)] {
        let c = sphere.curvature_at(u, v, T).unwrap();
        // Signed against the outward normal: the sphere curves away from it.
        assert!((c.max + 1.0 / 2.5).abs() < 1e-9, "max {}", c.max);
        assert!((c.min + 1.0 / 2.5).abs() < 1e-9, "min {}", c.min);
        assert!(c.is_umbilic(T));
        assert!((c.gaussian() - 1.0 / (2.5 * 2.5)).abs() < 1e-9);
        assert!(c.max_direction.vector().dot(c.normal.vector()).abs() < 1e-9);
        assert!(c.min_direction.vector().dot(c.normal.vector()).abs() < 1e-9);
        assert!(c.max_direction.vector().dot(c.min_direction.vector()).abs() < 1e-9);
    }
}

#[test]
fn a_cylinder_curves_round_and_not_along() {
    let cylinder =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, 4.0, T).unwrap(), (-1.0, 1.0)).unwrap();
    let c = cylinder.curvature_at(1.0, 0.3, T).unwrap();
    assert!(!c.is_umbilic(T));
    // The flat direction is the ruling, the curved one runs round.
    assert!(c.max.abs() < 1e-9, "along the ruling: {}", c.max);
    assert!((c.min + 0.25).abs() < 1e-9, "round: {}", c.min);
    assert!(
        c.max_direction.vector().z.abs() > 1.0 - 1e-9,
        "{:?}",
        c.max_direction
    );
    assert!(
        c.min_direction.vector().z.abs() < 1e-9,
        "{:?}",
        c.min_direction
    );
    assert!((c.mean() + 0.125).abs() < 1e-9);
    assert!(c.gaussian().abs() < 1e-12);
}

#[test]
fn a_torus_is_saddle_inside_and_dome_outside() {
    let torus = TorusSurface::new(Torus::new(Frame::WORLD, 5.0, 1.0, T).unwrap());
    // On the outer equator both principal curvatures are negative against
    // the outward normal: 1/minor round the tube and 1/(major + minor)
    // round the ring; on the inner equator the ring's turns positive.
    let outer = torus.curvature_at(0.7, 0.0, T).unwrap();
    let outer_pair = [outer.max, outer.min];
    assert!(
        outer_pair.iter().any(|k| (k + 1.0).abs() < 1e-9),
        "{outer:?}"
    );
    assert!(
        outer_pair.iter().any(|k| (k + 1.0 / 6.0).abs() < 1e-9),
        "{outer:?}"
    );
    assert!(outer.gaussian() > 0.0);
    let inner = torus.curvature_at(0.7, std::f64::consts::PI, T).unwrap();
    let inner_pair = [inner.max, inner.min];
    assert!(
        inner_pair.iter().any(|k| (k + 1.0).abs() < 1e-9),
        "{inner:?}"
    );
    assert!(
        inner_pair.iter().any(|k| (k - 1.0 / 4.0).abs() < 1e-9),
        "{inner:?}"
    );
    assert!(inner.gaussian() < 0.0, "a saddle inside");
}
