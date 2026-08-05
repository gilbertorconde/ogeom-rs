//! The boolean's API completions: scaled placements bake automatically and
//! run the analytic pipeline.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Transform};
use ogeom::mesh::Deflection;
use ogeom::topo::Model;

const T: Tolerances = Tolerances::millimetres();

/// A placement with scale no longer refuses: the boolean bakes it — every
/// surface restated exactly in its own analytic vocabulary, pcurves
/// re-derived — and the drilled result measures.
#[test]
fn a_scaled_placement_bakes_and_the_boolean_runs() {
    let mut model = Model::new();
    let unit = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let scaled = ogeom::algo::transformed(
        &mut model,
        &unit,
        Transform::scaling(Point::ORIGIN, 2.0, T).unwrap(),
    )
    .unwrap()
    .shape;

    // The bake alone is exact: the doubled box measures eight thousand.
    let baked = ogeom::algo::baked_shape(&mut model, &scaled, T)
        .unwrap()
        .shape;
    let v = ogeom::algo::volume_properties(&model, &baked, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((v - 8000.0).abs() < 1e-6, "the bake is exact: {v}");

    // And the boolean takes the scaled shape directly.
    let f = Frame::new(Point::new(10.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, f, 3.0, 22.0, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(&mut model, &scaled, &drill, T)
        .unwrap()
        .shape;
    let v = ogeom::algo::volume_properties(&model, &cut, Deflection::with_chord(1e-2).unwrap(), T)
        .unwrap()
        .mass;
    let expected = 8000.0 - core::f64::consts::PI * 9.0 * 20.0;
    assert!(
        (v - expected).abs() / expected < 1e-3,
        "the scaled cut measures: {v} vs {expected}"
    );
}
