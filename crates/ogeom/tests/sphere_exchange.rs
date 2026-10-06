//! Spheres through STEP and IGES: a sphere face bounded by its seam walked
//! both ways (with no edges at its poles, as an exchange file states it)
//! reads back whole, and a solid's spherical void stays a void.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

#[path = "support/walks.rs"]
mod walks;

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

/// The volume of the shape's mesh at a fine chord: the same measure on
/// both sides of the trip, whichever pcurves the reader fitted.
fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::mesh::triangulate(model, shape, Deflection::with_chord(0.01).unwrap(), T)
        .unwrap()
        .volume()
}

fn cases() -> Vec<(&'static str, Model, Shape)> {
    let mut out = Vec::new();
    let mut m = Model::new();
    let s = ogeom::algo::make_sphere(&mut m, at(0.0, 0.0, 0.0), 10.0, T)
        .unwrap()
        .shape;
    out.push(("whole sphere", m, s));

    let mut m = Model::new();
    let s = ogeom::algo::make_sphere(&mut m, at(0.0, 0.0, 0.0), 10.0, T)
        .unwrap()
        .shape;
    let b = ogeom::algo::make_box(&mut m, at(5.0, -20.0, -20.0), (40.0, 40.0, 40.0), T)
        .unwrap()
        .shape;
    let s = ogeom::boolean::cut(&mut m, &s, &b, T).unwrap().shape;
    out.push(("sphere cut on its side", m, s));

    let mut m = Model::new();
    let b = ogeom::algo::make_box(&mut m, at(-20.0, -20.0, -20.0), (40.0, 40.0, 40.0), T)
        .unwrap()
        .shape;
    let s = ogeom::algo::make_sphere(&mut m, at(0.0, 0.0, 0.0), 10.0, T)
        .unwrap()
        .shape;
    let s = ogeom::boolean::cut(&mut m, &b, &s, T).unwrap().shape;
    out.push(("spherical void", m, s));

    let mut m = Model::new();
    let a = ogeom::algo::make_sphere(&mut m, at(0.0, 0.0, 0.0), 10.0, T)
        .unwrap()
        .shape;
    let c = ogeom::algo::make_sphere(&mut m, at(12.0, 0.0, 0.0), 8.0, T)
        .unwrap()
        .shape;
    let s = ogeom::boolean::fuse(&mut m, &a, &c, T).unwrap().shape;
    out.push(("two spheres fused", m, s));
    out
}

fn holds(name: &str, format: &str, original: (&Model, &Shape), back: (&Model, &Shape)) {
    let (want, got) = (volume(original.0, original.1), volume(back.0, back.1));
    assert!(
        (got - want).abs() < want * 1e-5,
        "{name} through {format}: {got} against {want}"
    );
    assert!(
        ogeom::algo::check(back.0, back.1, T).unwrap().is_valid(),
        "{name} through {format}"
    );
    let shells = |m: &Model, s: &Shape| explore_unique(m, s, ShapeType::Shell).unwrap().len();
    assert_eq!(
        shells(back.0, back.1),
        shells(original.0, original.1),
        "{name} through {format}"
    );
}

#[test]
fn spheres_round_trip_through_step() {
    for (name, model, solid) in cases() {
        let mut document = ogeom::doc::Document::over(model.clone());
        document.add_part("part", solid.clone());
        let text = ogeom::io::write_step(&document, T).unwrap();
        let import = ogeom::io::read_step(&text, T).unwrap();
        assert_eq!(import.solids.len(), 1, "{name}");
        holds(
            name,
            "STEP",
            (&model, &solid),
            (import.document.model(), &import.solids[0]),
        );
    }
}

#[test]
fn spheres_round_trip_through_iges() {
    for (name, model, solid) in cases() {
        let mut document = ogeom::doc::Document::over(model.clone());
        document.add_part("part", solid.clone());
        let text = ogeom::io::write_iges(&document, T).unwrap();
        let import = ogeom::io::read_iges(&text, T).unwrap();
        assert_eq!(import.solids.len(), 1, "{name}");
        assert_eq!(
            walks::edges_walked_one_way(import.document.model(), &import.solids[0]),
            0,
            "{name}"
        );
        holds(
            name,
            "IGES",
            (&model, &solid),
            (import.document.model(), &import.solids[0]),
        );
    }
}
