//! A parameter that is not a number, or is infinite, is refused by every
//! curve and surface, periodic ones included: wrapping it would answer a
//! point that is not a number.
#![allow(clippy::unwrap_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::{CircleCurve, Curve3d as _, CylinderSurface, Surface as _, SurfaceGeometry};
use ogeom_math::{Circle, Cylinder, Frame};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_circle_refuses_a_parameter_that_is_not_finite() {
    let circle = CircleCurve::new(Circle::new(Frame::WORLD, 2.0, T).unwrap());
    for t in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(circle.point_at(t, T).is_err(), "{t}");
    }
    // A finite parameter a few turns out still wraps.
    let p = circle
        .point_at(4.0 * core::f64::consts::TAU + 0.5, T)
        .unwrap();
    let q = circle.point_at(0.5, T).unwrap();
    assert!(p.distance(q) < 1e-9);
}

#[test]
fn a_cylinder_refuses_a_parameter_that_is_not_finite() {
    let cylinder: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, 1.0, T).unwrap(), (0.0, 5.0))
            .unwrap()
            .into();
    for u in [f64::NAN, f64::INFINITY] {
        assert!(cylinder.point_at(u, 1.0, T).is_err(), "{u}");
    }
    let p = cylinder
        .point_at(-core::f64::consts::TAU + 0.25, 1.0, T)
        .unwrap();
    let q = cylinder.point_at(0.25, 1.0, T).unwrap();
    assert!(p.distance(q) < 1e-9);
}
