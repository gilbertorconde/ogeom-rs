//! An offset of an analytic surface is the analytic surface it names: the
//! generic offset's points stand on it, grown whichever way the basis
//! normal points.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::{
    ConeSurface, CylinderSurface, OffsetSurface, PlaneSurface, SphereSurface, Surface,
    SurfaceGeometry, TorusSurface,
};
use ogeom_math::{Cone, Cylinder, Frame, Plane, Point, Sphere, Torus};

const T: Tolerances = Tolerances::millimetres();

/// Every generic offset point, over a grid inside the basis's window, stands
/// on the analytic surface: its distance from it by `on`, the analytic
/// surface's own equation.
fn agrees(basis: SurfaceGeometry, distance: f64, on: impl Fn(&SurfaceGeometry, Point) -> f64) {
    let offset = OffsetSurface::new(basis, distance).unwrap();
    let analytic = offset.analytic(T).unwrap().expect("an analytic basis");
    let ((u0, u1), (v0, v1)) = offset.domain();
    let (u0, u1) = (u0.max(-5.0), u1.min(5.0));
    let (v0, v1) = (v0.max(-5.0), v1.min(5.0));
    // Inside the window: a ball's poles have no normal to offset along.
    for i in 1..6 {
        for j in 1..6 {
            let u = u0 + (u1 - u0) * f64::from(i) / 6.0;
            let v = v0 + (v1 - v0) * f64::from(j) / 6.0;
            let p = offset.point_at(u, v, T).unwrap();
            let off = on(&analytic, p);
            assert!(off.abs() < 1e-9, "{p:?} stands {off} off");
        }
    }
}

#[test]
fn offsets_of_the_analytics_are_analytics() {
    let axis_distance = |p: Point| p.x.hypot(p.y);
    agrees(
        PlaneSurface::over(Plane::new(Frame::WORLD), (-5.0, 5.0), (-5.0, 5.0))
            .unwrap()
            .into(),
        1.5,
        |s, p| match s {
            SurfaceGeometry::Plane(q) => q.plane().signed_distance_to(p),
            _ => f64::NAN,
        },
    );
    for d in [0.5, -0.5] {
        agrees(
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 2.0, T).unwrap(), (-3.0, 3.0))
                .unwrap()
                .into(),
            d,
            |s, p| match s {
                SurfaceGeometry::Cylinder(c) => axis_distance(p) - c.cylinder().radius(),
                _ => f64::NAN,
            },
        );
        agrees(
            SphereSurface::new(Sphere::new(Frame::WORLD, 2.0, T).unwrap()).into(),
            d,
            |s, p| match s {
                SurfaceGeometry::Sphere(b) => p.distance(Point::ORIGIN) - b.sphere().radius(),
                _ => f64::NAN,
            },
        );
        agrees(
            TorusSurface::new(Torus::new(Frame::WORLD, 5.0, 1.0, T).unwrap()).into(),
            d,
            |s, p| match s {
                SurfaceGeometry::Torus(t) => {
                    (axis_distance(p) - 5.0).hypot(p.z) - t.torus().minor_radius()
                }
                _ => f64::NAN,
            },
        );
    }
    agrees(
        ConeSurface::new(Cone::new(Frame::WORLD, 2.0, 0.3, T).unwrap(), (0.0, 4.0))
            .unwrap()
            .into(),
        0.25,
        |s, p| match s {
            SurfaceGeometry::Cone(c) => {
                // Distance from the cone's surface, square to its rulings.
                let cone = c.cone();
                (axis_distance(p) - cone.radius_at(p.z)) * cone.half_angle().cos()
            }
            _ => f64::NAN,
        },
    );
}
