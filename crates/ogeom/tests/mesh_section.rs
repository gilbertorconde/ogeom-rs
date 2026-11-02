//! A scan-sized mesh sectioned by planes: one plane and a stack of them.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

#[path = "support/cpu_time.rs"]
mod cpu_time;
use cpu_time::cpu_time;

use core::f64::consts::TAU;
use std::time::Duration;

use ogeom::math::{Direction, Plane, Point, Vector};
use ogeom::mesh::{section, sections};
use ogeom::topo::Triangulation;

/// A torus about the `z` axis, `rings` by `sides` quads, two triangles each.
fn torus(major: f64, minor: f64, rings: u32, sides: u32) -> Triangulation {
    let mut mesh = Triangulation::new();
    for i in 0..rings {
        let u = f64::from(i) / f64::from(rings) * TAU;
        for j in 0..sides {
            let v = f64::from(j) / f64::from(sides) * TAU;
            let r = major + minor * v.cos();
            mesh.positions
                .push(Point::new(r * u.cos(), r * u.sin(), minor * v.sin()));
        }
    }
    let at = |i: u32, j: u32| (i % rings) * sides + (j % sides);
    for i in 0..rings {
        for j in 0..sides {
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            mesh.triangles.push([a, b, c]);
            mesh.triangles.push([a, c, d]);
        }
    }
    mesh
}

/// A million triangles, cut across the tube: two closed loops, the outer
/// facing away from the axis and the inner towards it, in under
/// 100 ms of processor time on one thread, in the test profile as well as
/// an optimised build.
#[test]
fn a_million_triangles_section_within_a_tenth_of_a_second() {
    let mesh = torus(10.0, 3.0, 1000, 500);
    assert_eq!(mesh.triangles.len(), 1_000_000);
    let plane = Plane::through(
        Point::new(0.0, 0.0, 0.37),
        Direction::new(
            Vector::new(0.01, 0.02, 1.0),
            ogeom::core::Tolerances::millimetres(),
        )
        .unwrap(),
    );
    let (s, took) = cpu_time(|| section(&mesh, &plane).unwrap());
    assert!(s.unlinked.is_empty());
    assert_eq!(s.chains.len(), 2);
    for c in &s.chains {
        assert!(c.closed);
        let r = c.points[0].x.hypot(c.points[0].y);
        let facing = c.outward[0].x * c.points[0].x + c.outward[0].y * c.points[0].y;
        assert_eq!(facing > 0.0, r > 10.0, "{r} {facing}");
    }
    eprintln!("one section of 1M triangles: {took:?}");
    assert!(
        took < Duration::from_millis(100),
        "{took:?} of processor time"
    );

    // Clear of the tube both ends, and through it between.
    let mut offsets = vec![-4.5];
    offsets.extend((0..14).map(|k| -2.0 + 0.25 * f64::from(k)));
    offsets.push(4.5);
    let (stack, took) = cpu_time(|| sections(&mesh, &plane, &offsets).unwrap());
    eprintln!("a stack of 16 sections of 1M triangles: {took:?}");
    assert_eq!(stack.len(), 16);
    // Below the tube and above it the planes miss; across it, two loops.
    assert!(stack[0].chains.is_empty() && stack[15].chains.is_empty());
    assert!(stack[1..15].iter().all(|s| s.chains.len() == 2));
}
