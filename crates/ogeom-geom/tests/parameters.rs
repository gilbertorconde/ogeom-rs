//! Parameters at the edges of what a curve or surface accepts: one that is
//! not a number or is infinite is refused, periodic geometry included, and
//! a trim's parameters move with its basis's under a transform.
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

/// A trim moves with its basis's parameters under a scaling: a cylinder
/// stored over heights from zero, trimmed to heights 1 to 3, scaled by two
/// about the origin, is trimmed to heights 2 to 6; and a trim of a trim of
/// a line keeps its place along the line.
#[test]
fn a_trim_scales_with_its_basis() {
    use ogeom_geom::{Curve, LineCurve, Transformable as _, TrimmedCurve, TrimmedSurface};
    use ogeom_math::{Axis, Direction, Point, Transform};
    let twice = Transform::scaling(Point::ORIGIN, 2.0, T).unwrap();
    let cylinder: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, 1.0, T).unwrap(), (0.0, 10.0))
            .unwrap()
            .into();
    let trimmed: SurfaceGeometry =
        TrimmedSurface::new(cylinder, (0.0, core::f64::consts::PI), (1.0, 3.0), T)
            .unwrap()
            .into();
    let scaled = trimmed.transformed(&twice, T).unwrap();
    let (_, (v0, v1)) = scaled.domain();
    assert!(
        (v0 - 2.0).abs() < 1e-12 && (v1 - 6.0).abs() < 1e-12,
        "{v0} {v1}"
    );
    let bottom = scaled.point_at(0.0, v0, T).unwrap();
    assert!((bottom.z - 2.0).abs() < 1e-12, "{bottom:?}");

    let line: Curve = LineCurve::over(Axis::new(Point::ORIGIN, Direction::X), 0.0, 10.0)
        .unwrap()
        .into();
    let once: Curve = TrimmedCurve::new(line, 2.0, 8.0, T).unwrap().into();
    let twice_trimmed: Curve = TrimmedCurve::new(once, 3.0, 5.0, T).unwrap().into();
    let moved = twice_trimmed.transformed(&twice, T).unwrap();
    let (a, b) = moved.domain();
    let (pa, pb) = (moved.point_at(a, T).unwrap(), moved.point_at(b, T).unwrap());
    assert!(
        (pa.x - 6.0).abs() < 1e-9 && (pb.x - 10.0).abs() < 1e-9,
        "{pa:?} {pb:?}"
    );
}
