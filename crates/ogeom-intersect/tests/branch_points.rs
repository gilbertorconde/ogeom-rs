//! Branch-point stitching, pinned on the classical configuration: two
//! equal perpendicular cylinders meet in two ellipses that cross at two
//! branch points, and the marcher's honest fragments reassemble into
//! exactly those ellipses.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::{CylinderSurface, SurfaceGeometry};
use ogeom_intersect::{Marching, Stopped, branches};
use ogeom_math::{Cylinder, Direction, Frame, Point};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn equal_cylinders_come_back_as_two_closed_ellipses() {
    let a: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, 1.0, T).unwrap(), (-3.0, 3.0))
            .unwrap()
            .into();
    let along_x = Frame::new(Point::ORIGIN, Direction::X, Direction::Y, T).unwrap();
    let b: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(along_x, 1.0, T).unwrap(), (-3.0, 3.0))
            .unwrap()
            .into();

    let found = branches(&a, &b, Marching::default(), T).unwrap();
    assert_eq!(found.len(), 2, "two ellipses, not a pile of fragments");
    for branch in &found {
        assert_eq!(
            branch.stopped,
            Stopped::Closed,
            "each closes through the branch points"
        );
        let length: f64 = branch
            .points
            .windows(2)
            .map(|pair| pair[0].distance(pair[1]))
            .sum();
        // The √2-by-1 ellipse's perimeter, which has no closed form but
        // has a value: 7.6404 by series. The polyline inscribes it.
        assert!(
            (length - 7.6404).abs() < 5e-3,
            "an ellipse's worth of curve, measured: {length}"
        );
        // Every point lies on both cylinders.
        for p in &branch.points {
            let on_a = (p.x.hypot(p.y) - 1.0).abs();
            let on_b = (p.z.hypot(p.y) - 1.0).abs();
            assert!(on_a < 1e-5 && on_b < 1e-5, "off-surface point {p:?}");
        }
    }
}

#[test]
fn tangential_contact_traces_as_the_circle_it_is() {
    use ogeom_geom::SphereSurface;
    use ogeom_intersect::{Contact, trace_tangential};
    use ogeom_math::Sphere;

    // A unit sphere seated inside a unit cylinder: tangent along the whole
    // equator, no transversal crossing anywhere — the case the crossing
    // walker honestly refuses and the valley walker owns.
    let drum: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, 1.0, T).unwrap(), (-2.0, 2.0))
            .unwrap()
            .into();
    let ball: SurfaceGeometry =
        SphereSurface::new(Sphere::centred(Point::ORIGIN, 1.0, T).unwrap()).into();

    // Seed at the contact: cylinder (u=0, v=0) is (1,0,0); the sphere's
    // equator holds the same point.
    use ogeom_geom::Surface as _;
    let seed_b = {
        // Find the sphere parameters of (1,0,0) by projection: equator at
        // whatever convention the chart uses.
        let ((ua, ub), (va, vb)) = ball.domain();
        let mut best = ((ua, va), f64::INFINITY);
        for i in 0..=64 {
            for j in 0..=64 {
                let u = ua + (ub - ua) * f64::from(i) / 64.0;
                let v = va + (vb - va) * f64::from(j) / 64.0;
                if let Ok(p) = ball.point_at(u, v, T) {
                    let d = p.distance(Point::new(1.0, 0.0, 0.0));
                    if d < best.1 {
                        best = ((u, v), d);
                    }
                }
            }
        }
        best.0
    };
    let contact = Contact {
        on_a: (0.0, 0.0),
        on_b: seed_b,
        point: Point::new(1.0, 0.0, 0.0),
    };

    let traced = trace_tangential(&drum, &ball, contact, Marching::default(), T).unwrap();
    assert_eq!(traced.stopped, Stopped::Closed, "the equator closes");
    let length: f64 = traced
        .points
        .windows(2)
        .map(|pair| pair[0].distance(pair[1]))
        .sum();
    assert!(
        (length - core::f64::consts::TAU).abs() < 1e-2,
        "one full equator: {length}"
    );
    for p in &traced.points {
        assert!(
            (p.x.hypot(p.y) - 1.0).abs() < 1e-5 && p.z.abs() < 1e-5,
            "on the equator: {p:?}"
        );
    }
}
