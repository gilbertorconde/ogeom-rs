//! Bodies under reflecting placements: the chart's natural normal flips
//! against every orientation flag, and each consumer must fold the
//! handedness back in.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Frame, Point};

const T: Tolerances = Tolerances::millimetres();

fn vol(model: &ogeom_topo::Model, s: &ogeom_topo::Shape) -> f64 {
    ogeom_algo::volume_properties(
        model,
        s,
        ogeom_mesh::Deflection {
            chord: 1e-3,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

fn block_and_mirror(model: &mut ogeom_topo::Model) -> (ogeom_topo::Shape, ogeom_topo::Shape) {
    let block = ogeom_algo::make_box(
        model,
        Frame::new(
            Point::new(1.0, 0.0, 0.0),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap(),
        (10.0, 10.0, 10.0),
        T,
    )
    .unwrap()
    .shape;
    let mirror =
        ogeom_math::Transform::plane_mirror(Point::new(1.0, 0.0, 0.0), ogeom_math::Direction::X);
    let mirrored = model.placed(&block, mirror);
    (block, mirrored)
}

#[test]
fn a_mirrored_body_measures_right_side_out() {
    let mut model = ogeom_topo::Model::new();
    let (_, mirrored) = block_and_mirror(&mut model);
    let properties = ogeom_algo::volume_properties(
        &model,
        &mirrored,
        ogeom_mesh::Deflection {
            chord: 1e-3,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap();
    assert!((properties.mass - 1000.0).abs() < 1e-6);
    assert!(
        (properties.centre.x + 4.0).abs() < 1e-6,
        "the centre mirrored across x = 1: {:?}",
        properties.centre
    );

    // And the bake restates it right side out.
    let baked = ogeom_algo::baked_shape(&mut model, &mirrored, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &baked.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!((vol(&model, &baked.shape) - 1000.0).abs() < 1e-6);
}

#[test]
fn a_body_fuses_with_its_mirror_across_the_shared_face() {
    let mut model = ogeom_topo::Model::new();
    let (block, mirrored) = block_and_mirror(&mut model);
    let fused = ogeom_bool::fuse(&mut model, &block, &mirrored, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &fused.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!(
        (vol(&model, &fused.shape) - 2000.0).abs() < 1e-6,
        "the halves joined across their shared face"
    );
}

#[test]
fn a_mirrored_drum_fuses_through_its_fitted_images() {
    // A curved body under a reflection: the boolean restates it in world
    // coordinates first, and whatever pcurves that restatement had to fit
    // carry their slop on the record — the melt's snap reaches it, or the
    // contact dangles a hair from the boundary it paved.
    let mut model = ogeom_topo::Model::new();
    let drum = ogeom_algo::make_cylinder(
        &mut model,
        Frame::new(
            Point::new(3.0, 0.0, 0.0),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap(),
        2.0,
        6.0,
        T,
    )
    .unwrap()
    .shape;
    let mirror = ogeom_math::Transform::plane_mirror(Point::ORIGIN, ogeom_math::Direction::X);
    let mirrored = model.placed(&drum, mirror);
    let block = ogeom_algo::make_box(
        &mut model,
        Frame::new(
            Point::new(-2.0, -5.0, 0.0),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap(),
        (2.0, 10.0, 6.0),
        T,
    )
    .unwrap()
    .shape;
    let fused = ogeom_bool::fuse(&mut model, &block, &mirrored, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &fused.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    // Block plus drum less the circular segment they share.
    let pi = core::f64::consts::PI;
    let theta = 2.0 * core::f64::consts::FRAC_PI_3;
    let segment = 4.0 * (theta - theta.sin()) / 2.0 * 6.0;
    let expected = 6.0_f64.mul_add(2.0 * 10.0, pi * 4.0 * 6.0) - segment;
    let measured = vol(&model, &fused.shape);
    assert!(
        (measured - expected).abs() < expected * 1e-3,
        "fused volume {measured} against {expected}"
    );
}
