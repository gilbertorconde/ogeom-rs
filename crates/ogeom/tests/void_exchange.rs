//! A solid's cavity through STEP: written as a void, and read as one
//! whichever way round the file's flags put it.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Filter, Model, Shape, ShapeType, explore, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A 10 mm cube with a 2 mm cubic void at its middle.
fn hollow_cube() -> (Model, Shape) {
    let mut model = Model::new();
    let outer = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(4.0, 4.0, 4.0), Direction::Z, Direction::X, T).unwrap();
    let inner = ogeom::algo::make_box(&mut model, at, (2.0, 2.0, 2.0), T)
        .unwrap()
        .shape;
    let o = explore(&model, &outer, Filter::OfType(ShapeType::Shell)).unwrap()[0].clone();
    let i = explore(&model, &inner, Filter::OfType(ShapeType::Shell)).unwrap()[0].clone();
    let solid = ogeom::algo::make_solid(&mut model, &[o, i.reversed()])
        .unwrap()
        .shape;
    (model, solid)
}

fn written(model: &Model, solid: &Shape) -> String {
    let mut document = ogeom::doc::Document::over(model.clone());
    document.add_part("cube", solid.clone());
    ogeom::io::write_step(&document, T).unwrap()
}

/// One solid of two shells, 992 mm³ about the cube's centre, valid.
fn reads_as_the_hollow_cube(text: &str) {
    let import = ogeom::io::read_step(text, T).unwrap();
    assert_eq!(import.solids.len(), 1);
    let (model, solid) = (import.document.model(), &import.solids[0]);
    assert_eq!(
        explore_unique(model, solid, ShapeType::Shell)
            .unwrap()
            .len(),
        2
    );
    let properties =
        ogeom::algo::volume_properties(model, solid, Deflection::default(), T).unwrap();
    assert!(
        (properties.mass - 992.0).abs() < 1e-6,
        "{}",
        properties.mass
    );
    assert!(properties.centre.distance(Point::new(5.0, 5.0, 5.0)) < 1e-9);
    let diagnosis = ogeom::algo::check(model, solid, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
}

#[test]
fn a_solid_with_a_void_is_written_with_it() {
    let (model, solid) = hollow_cube();
    let text = written(&model, &solid);
    assert!(text.contains("BREP_WITH_VOIDS"));
    assert_eq!(text.matches("ORIENTED_CLOSED_SHELL").count(), 1);
    assert!(text.contains(".F.);") || text.contains(".F.)"));
    reads_as_the_hollow_cube(&text);
}

/// A void written the wrong way out (its flags turn it to add material)
/// is still the cavity the solid's void list names.
#[test]
fn a_void_written_the_wrong_way_out_reads_as_a_void() {
    let (model, solid) = hollow_cube();
    let text = written(&model, &solid);
    let line = text
        .lines()
        .find(|l| l.contains("ORIENTED_CLOSED_SHELL"))
        .unwrap()
        .to_string();
    let flipped = text.replace(&line, &line.replace(".F.", ".T."));
    assert_ne!(flipped, text);
    reads_as_the_hollow_cube(&flipped);
}
