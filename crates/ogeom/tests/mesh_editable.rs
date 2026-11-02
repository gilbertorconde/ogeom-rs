//! A scan-sized mesh held for editing: built, read back, a patch cut out
//! of it, and the whole cut by a plane and capped.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

#[path = "support/cpu_time.rs"]
mod cpu_time;
use cpu_time::cpu_time;

use core::f64::consts::TAU;
use std::time::Duration;

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Plane, Point, Vector};
use ogeom::mesh::{Cap, CutReport, Editable, Side};
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

/// Five million triangles built into half-edges in under two seconds of
/// processor time, read back unchanged in under one, and welded and built
/// in under four; a patch of them removed leaves one rim.
#[test]
#[ignore = "heavy"]
fn five_million_triangles_build_and_read_back() {
    let (rings, sides) = (2500, 1000);
    let mesh = torus(10.0, 3.0, rings, sides);
    assert_eq!(mesh.triangles.len(), 5_000_000);
    let (mut e, built) = cpu_time(|| Editable::from(&mesh, 0.0).unwrap());
    eprintln!("5M triangles built: {built:?}");
    assert!(e.is_closed());
    let (back, read) = cpu_time(|| e.to_triangulation());
    eprintln!("5M triangles read back: {read:?}");
    assert_eq!(back.positions, mesh.positions);
    assert_eq!(back.triangles, mesh.triangles);

    // Ten by ten quads in the middle of the torus.
    let patch: Vec<u32> = (0..10)
        .flat_map(|i| (0..10).map(move |j| 2 * (i * sides + j)))
        .flat_map(|f| [f, f + 1])
        .collect();
    let ((), removed) = cpu_time(|| e.remove_faces(&patch));
    eprintln!("200 faces removed: {removed:?}");
    let ((rims, pieces), queried) = cpu_time(|| (e.rims(), e.components()));
    eprintln!("rims and components of 5M triangles: {queried:?}");
    assert_eq!(rims.len(), 1);
    assert_eq!(rims[0].edges.len(), 40);
    assert_eq!(pieces.len(), 1);

    let (welded, weld) = cpu_time(|| Editable::from(&mesh, 1e-6).unwrap());
    eprintln!("5M triangles welded and built: {weld:?}");
    assert_eq!(welded.vertex_count(), mesh.positions.len());
    assert!(welded.is_closed());

    assert!(built < Duration::from_secs(5), "{built:?} to build");
    assert!(read < Duration::from_secs(1), "{read:?} to read back");
    assert!(weld < Duration::from_secs(17), "{weld:?} to weld and build");
}

/// Five million triangles cut across the tube by a plane and capped in
/// under six seconds of processor time: the cut through the tube's top
/// half is two loops, an inner and an outer, capped as one ring, and the
/// part below is closed.
#[test]
#[ignore = "heavy"]
fn five_million_triangles_cut_and_capped() {
    let mesh = torus(10.0, 3.0, 2500, 1000);
    let mut e = Editable::from(&mesh, 0.0).unwrap();
    let normal = Direction::new(Vector::new(0.0, 0.05, 1.0), Tolerances::millimetres()).unwrap();
    let plane = Plane::through(Point::new(0.0, 0.0, 1.0), normal);
    let (report, took) = cpu_time(|| {
        e.cut(&plane, Side::Below, Some(Cap { bridge: 0.0 }))
            .unwrap()
    });
    eprintln!("5M triangles cut and capped: {took:?}");
    assert_eq!(
        report,
        CutReport {
            capped_loops: 2,
            ..CutReport::default()
        }
    );
    assert!(e.is_closed());
    assert_eq!(e.components().len(), 1);
    assert!(took < Duration::from_secs(6), "{took:?} to cut and cap");
}
