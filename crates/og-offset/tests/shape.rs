//! The solid offset and the shells cut from it.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use og_core::Tolerances;
use og_math::{Frame, Point};
use og_topo::{Filter, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

fn fine() -> og_mesh::Deflection {
    og_mesh::Deflection {
        chord: 1e-4,
        ..og_mesh::Deflection::default()
    }
}

fn volume(model: &og_topo::Model, shape: &og_topo::Shape) -> f64 {
    og_algo::volume_properties(model, shape, fine(), T)
        .unwrap()
        .mass
}

#[test]
fn a_box_offset_outward_is_the_bigger_box() {
    let mut model = og_topo::Model::new();
    let block = og_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let result = og_offset::offset_shape(&mut model, &block.shape, 0.5, T).unwrap();

    let diagnosis = og_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - 27.0).abs() < 1e-9,
        "outward box offset volume {measured} against 27"
    );
    // Topology preserved one-for-one.
    let faces = explore(&model, &result.shape, Filter::OfType(ShapeType::Face)).unwrap();
    assert_eq!(faces.len(), 6);
}

#[test]
fn a_box_offset_inward_is_the_smaller_box() {
    let mut model = og_topo::Model::new();
    let block = og_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let result = og_offset::offset_shape(&mut model, &block.shape, -0.5, T).unwrap();
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - 1.0).abs() < 1e-9,
        "inward box offset volume {measured} against 1"
    );
}

#[test]
fn a_cylinder_offset_grows_radius_and_caps_alike() {
    let mut model = og_topo::Model::new();
    let drum = og_algo::make_cylinder(&mut model, Frame::WORLD, 1.0, 2.0, T).unwrap();
    let result = og_offset::offset_shape(&mut model, &drum.shape, 0.5, T).unwrap();

    let diagnosis = og_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let expected = core::f64::consts::PI * 1.5 * 1.5 * 3.0;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 5e-3,
        "cylinder offset volume {measured} against {expected}"
    );
}

#[test]
fn an_offset_that_swallows_the_box_is_refused() {
    let mut model = og_topo::Model::new();
    let block = og_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    assert!(og_offset::offset_shape(&mut model, &block.shape, -1.5, T).is_err());
}

/// The face of `solid` whose interior contains `probe`.
fn face_at(model: &og_topo::Model, solid: &og_topo::Shape, probe: Point) -> og_topo::Shape {
    explore(model, solid, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find(|f| {
            og_algo::classify_on_face(model, f, probe, fine(), T)
                .map(|c| c == og_algo::Containment::In)
                .unwrap_or(false)
        })
        .expect("the solid has a face there")
}

#[test]
fn a_shelled_box_keeps_its_walls_and_opens_its_top() {
    let mut model = og_topo::Model::new();
    let block = og_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let top = face_at(&model, &block.shape, Point::new(1.0, 1.0, 2.0));

    let t = 0.2;
    let result =
        og_offset::make_thick_solid(&mut model, &block.shape, std::slice::from_ref(&top), t, T)
            .unwrap();

    let diagnosis = og_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let expected = 8.0 - 1.6 * 1.6 * 1.8;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 1e-9,
        "box shell volume {measured} against {expected}"
    );

    // Five outer walls, the top ring, and five cavity walls.
    let faces = explore(&model, &result.shape, Filter::OfType(ShapeType::Face)).unwrap();
    assert_eq!(faces.len(), 11);
    assert!(result.history.is_deleted(&top));

    // The cavity is real: its centre is outside the material.
    let inside =
        og_algo::classify_in_solid_exact(&model, &result.shape, Point::new(1.0, 1.0, 1.0), T)
            .unwrap();
    assert_eq!(inside, og_algo::Containment::Out);
    // And the wall is material.
    let wall =
        og_algo::classify_in_solid_exact(&model, &result.shape, Point::new(0.1, 1.0, 1.0), T)
            .unwrap();
    assert_eq!(wall, og_algo::Containment::In);
}

#[test]
fn a_shelled_cylinder_becomes_a_cup() {
    let mut model = og_topo::Model::new();
    let drum = og_algo::make_cylinder(&mut model, Frame::WORLD, 1.0, 2.0, T).unwrap();
    let top = face_at(&model, &drum.shape, Point::new(0.0, 0.0, 2.0));

    let t = 0.2;
    let result =
        og_offset::make_thick_solid(&mut model, &drum.shape, std::slice::from_ref(&top), t, T)
            .unwrap();

    let diagnosis = og_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let pi = core::f64::consts::PI;
    let expected = pi * 2.0 - pi * 0.8 * 0.8 * 1.8;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 1e-3,
        "cup volume {measured} against {expected}"
    );
    assert!(result.history.is_deleted(&top));
}
