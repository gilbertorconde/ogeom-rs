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
