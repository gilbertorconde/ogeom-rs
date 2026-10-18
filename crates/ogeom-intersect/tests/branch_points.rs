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
    // equator, no transversal crossing anywhere: the case the crossing
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

/// Two walls swept square to each other from quadratic profiles, each
/// profile ending on its apex: one along x, its apex line y = 5.64 at
/// z = -15.68, one down z from z = -11.76, its apex line y = 5.64 at
/// x = -11.21. They touch where the apex lines cross, normals parallel,
/// and the section runs into that touch ever more shallowly until the walk
/// gives out short of it. One profile's apex is drawn a nanometre past its
/// end, as rounding leaves it in a part read from a file: the touch lies a
/// hair outside that wall's chart, where the search for it must reach onto
/// the bound and back. The section ends on the touch, to the confusion
/// distance either way the apex is pushed.
#[test]
fn a_section_into_a_touch_where_both_walls_end_ends_on_it() {
    use ogeom_geom::{BSplineCurve, Curve, ExtrusionSurface};
    use ogeom_math::{KnotVector, Vector};
    let quadratic = |control: [Point; 3]| -> Curve {
        BSplineCurve::new(
            KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap(),
            control.to_vec(),
            T,
        )
        .unwrap()
        .into()
    };
    let touch = Point::new(-11.21, 5.64, -15.68);
    // A nanometre off the middle control point tilts the profile's end
    // off level by two in a billion, and the apex moves just past that
    // end: lowered, on the wall along x, whose profile ends on its apex,
    // and on the wall down z, whose profile starts on it.
    for (along_x, down_z) in [(-1e-9, 0.0), (0.0, 1e-9)] {
        let a: SurfaceGeometry = ExtrusionSurface::new(
            quadratic([
                Point::new(-13.13, -3.29, 0.0),
                Point::new(-13.13, 5.64 + along_x, 0.0),
                Point::new(-13.13, 5.64, -15.68),
            ]),
            Direction::X,
            22.96,
        )
        .unwrap()
        .into();
        let b: SurfaceGeometry = ExtrusionSurface::new(
            quadratic([
                Point::new(-11.21, 5.64, -11.76),
                Point::new(-9.21, 5.64 - down_z, -11.76),
                Point::new(-9.21, 3.64, -11.76),
            ]),
            Direction::new(Vector::new(0.0, 0.0, -1.0), T).unwrap(),
            3.92,
        )
        .unwrap()
        .into();
        // A chord a boolean marches a part of this size at.
        let marching = Marching {
            chord: 1e-6,
            ..Marching::default()
        };
        let found = branches(&a, &b, marching, T).unwrap();
        assert_eq!(found.len(), 1, "one section from the touch to the top");
        let ends = [
            found[0].points[0],
            found[0].points[found[0].points.len() - 1],
        ];
        let miss = ends
            .iter()
            .map(|p| p.distance(touch))
            .fold(f64::INFINITY, f64::min);
        assert!(
            miss <= T.confusion(),
            "the section stops {miss} short of the touch ({along_x}, {down_z})"
        );
    }
}
