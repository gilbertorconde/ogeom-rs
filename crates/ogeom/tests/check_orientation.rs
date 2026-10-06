//! What `check` sees that topology alone cannot: a face turned inside out.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use std::f64::consts::TAU;

use ogeom::algo::{
    check, make_box, make_face, make_polygon, make_revolution, make_shell, make_solid,
};
use ogeom::core::Tolerances;
use ogeom::geom::PlaneSurface;
use ogeom::math::{Axis, Direction, Frame, Plane, Point};
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

/// A drum of radius 5 and `height` along `y`, its top a cone rising
/// `rise` to the axis: a cap all but flat and not planar.
fn domed_drum(model: &mut Model, height: f64, rise: f64) -> Shape {
    let pts = [(0.0, 0.0), (5.0, 0.0), (5.0, height), (0.0, height + rise)]
        .map(|(x, y)| Point::new(x, y, 0.0));
    let wire = make_polygon(model, &pts, true, T).unwrap().shape;
    let plane = PlaneSurface::new(Plane::new(Frame::WORLD)).into();
    let profile = make_face(model, plane, &[wire], T).unwrap().shape;
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Y,
    };
    make_revolution(model, &profile, axis, TAU, T)
        .unwrap()
        .shape
}

/// Each face of a drum with a nearly flat domed cap turned to face into
/// the drum is named, and only that face: the cap, the wall and the base,
/// on a tall drum and on one thinner than the probe's first step.
#[test]
fn a_face_of_a_drum_with_a_nearly_flat_cap_turned_inside_out_is_named() {
    for (height, rise) in [(10.0, 0.05), (10.0, 1.0), (0.05, 0.02)] {
        let mut model = Model::new();
        let drum = domed_drum(&mut model, height, rise);
        let diagnosis = check(&model, &drum, T).unwrap();
        assert!(diagnosis.is_valid(), "{height} {rise}: {diagnosis}");
        let faces = explore_unique(&model, &drum, ShapeType::Face).unwrap();
        assert_eq!(faces.len(), 3);
        for turned in 0..faces.len() {
            let mut held = faces.clone();
            held[turned] = held[turned].reversed();
            let shell = make_shell(&mut model, &held).unwrap().shape;
            let flipped = make_solid(&mut model, &[shell]).unwrap().shape;
            let diagnosis = check(&model, &flipped, T).unwrap();
            let named: Vec<_> = diagnosis
                .problems
                .iter()
                .filter(|p| p.what.contains("orientation is reversed"))
                .map(|p| p.at.node())
                .collect();
            assert_eq!(
                named,
                [faces[turned].node()],
                "{height} {rise}, face {turned}: {diagnosis}"
            );
        }
    }
}
