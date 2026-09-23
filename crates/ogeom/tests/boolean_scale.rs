//! A boolean's cost follows what the tool touches, not the size of the
//! solid it cuts: a solid of thousands of small faces takes a thin drill in
//! about the time its faces take to copy.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use std::time::{Duration, Instant};

use ogeom::algo::{MeshSolidOptions, check, solid_from_mesh, volume_properties};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, ShapeType, Triangulation, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A closed slab over an `n` by `n` grid of unit cells: flat underneath and
/// on its four sides, and on top a flat disc of radius `flat` round the
/// middle, rising in a rippled bowl outside it, so every top triangle off
/// the disc is its own face.
fn rippled_slab(n: u32, height: f64, flat: f64) -> Triangulation {
    let centre = f64::from(n) / 2.0;
    let top_z = |x: f64, y: f64| {
        let d = (x - centre).hypot(y - centre);
        let rise = (d - flat).max(0.0);
        height + 0.05 * rise * rise * (1.0 + 0.3 * (1.7 * x).sin() * (1.3 * y).cos())
    };
    let mut mesh = Triangulation::new();
    let side = n + 1;
    for layer in 0..2 {
        for j in 0..side {
            for i in 0..side {
                let (x, y) = (f64::from(i), f64::from(j));
                let z = if layer == 0 { top_z(x, y) } else { 0.0 };
                mesh.positions.push(Point::new(x, y, z));
            }
        }
    }
    let top = |i: u32, j: u32| j * side + i;
    let bottom = |i: u32, j: u32| side * side + j * side + i;
    for j in 0..n {
        for i in 0..n {
            let (a, b, c, d) = (top(i, j), top(i + 1, j), top(i + 1, j + 1), top(i, j + 1));
            mesh.triangles.push([a, b, c]);
            mesh.triangles.push([a, c, d]);
            let (a, b, c, d) = (
                bottom(i, j),
                bottom(i + 1, j),
                bottom(i + 1, j + 1),
                bottom(i, j + 1),
            );
            mesh.triangles.push([a, c, b]);
            mesh.triangles.push([a, d, c]);
        }
    }
    // The rim, walked anticlockwise from above: each step's wall faces out.
    let mut rim: Vec<(u32, u32)> = Vec::new();
    rim.extend((0..n).map(|i| (i, 0)));
    rim.extend((0..n).map(|j| (n, j)));
    rim.extend((1..=n).rev().map(|i| (i, n)));
    rim.extend((1..=n).rev().map(|j| (0, j)));
    for k in 0..rim.len() {
        let (p, q) = (rim[k], rim[(k + 1) % rim.len()]);
        let (bp, bq, tq, tp) = (
            bottom(p.0, p.1),
            bottom(q.0, q.1),
            top(q.0, q.1),
            top(p.0, p.1),
        );
        mesh.triangles.push([bp, bq, tq]);
        mesh.triangles.push([bp, tq, tp]);
    }
    mesh
}

/// A converted mesh of a couple of thousand faces is drilled through its
/// flat middle in seconds, valid, and short by exactly the bore.
#[test]
fn a_thin_drill_through_a_many_faced_solid_is_quick() {
    let (n, height, flat, radius) = (32, 10.0, 6.0, 2.5);
    let mut model = Model::new();
    let out = solid_from_mesh(
        &mut model,
        &rippled_slab(n, height, flat),
        &MeshSolidOptions::default(),
        T,
    )
    .unwrap();
    assert!(out.closed, "{:?}", out.report);
    let faces = explore_unique(&model, &out.shape, ShapeType::Face)
        .unwrap()
        .len();
    assert!(faces > 1500, "{faces}");
    let before = volume_properties(&model, &out.shape, Deflection::default(), T)
        .unwrap()
        .mass;

    let middle = f64::from(n) / 2.0;
    let frame = Frame::new(
        Point::new(middle, middle, -500.0),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, frame, radius, 1000.0, T)
        .unwrap()
        .shape;
    let started = Instant::now();
    let cut = ogeom::boolean::cut(&mut model, &out.shape, &drill, T)
        .unwrap()
        .shape;
    let took = started.elapsed();
    let diagnosis = check(&model, &cut, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let after = volume_properties(&model, &cut, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass;
    let bore = core::f64::consts::PI * radius * radius * height;
    assert!(
        ((before - after) - bore).abs() / bore < 1e-3,
        "removed {} against {bore}",
        before - after
    );
    assert!(
        took < Duration::from_secs(10),
        "{faces} faces took {took:?}"
    );
}
