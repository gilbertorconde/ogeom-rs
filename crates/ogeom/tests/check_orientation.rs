//! What `check` sees that topology alone cannot: a face turned inside out.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_box, make_shell, make_solid};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn cube(model: &mut Model, at: (f64, f64, f64)) -> Shape {
    let frame = Frame::new(Point::new(at.0, at.1, at.2), Direction::Z, Direction::X, T).unwrap();
    make_box(model, frame, (2.0, 3.0, 4.0), T).unwrap().shape
}

#[test]
fn a_box_is_valid_and_one_with_a_face_turned_inside_out_is_not() {
    let mut model = Model::new();
    let block = cube(&mut model, (0.0, 0.0, 0.0));
    assert!(check(&model, &block, T).unwrap().is_valid());

    let mut faces = explore_unique(&model, &block, ShapeType::Face).unwrap();
    faces[2] = faces[2].reversed();
    let shell = make_shell(&mut model, &faces).unwrap().shape;
    let flipped = make_solid(&mut model, &[shell]).unwrap().shape;
    let diagnosis = check(&model, &flipped, T).unwrap();
    assert!(!diagnosis.is_valid());
    let inward: Vec<_> = diagnosis
        .problems
        .iter()
        .filter(|p| p.kind == ShapeType::Face && p.what.contains("orientation is reversed"))
        .collect();
    assert_eq!(inward.len(), 1, "{diagnosis}");
    assert_eq!(inward[0].at.node(), faces[2].node());
}
