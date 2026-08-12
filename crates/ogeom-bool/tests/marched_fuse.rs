//! Booleans whose sections only the marching machinery can trace.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Frame, Point};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &ogeom_topo::Model, shape: &ogeom_topo::Shape, chord: f64) -> f64 {
    ogeom_algo::volume_properties(
        model,
        shape,
        ogeom_mesh::Deflection {
            chord,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

#[test]
fn a_branch_standing_on_a_drum_fuses_through_the_marched_seam() {
    // The seam winds the branch's whole chart and is metres of marching at
    // the old fixed chord: the walk has to close within its budget, the
    // fitted section has to meet the tolerance it states, and the branch
    // wall has to split along it.
    let mut model = ogeom_topo::Model::new();
    let main = Frame::new(
        Point::new(-20.0, 0.0, 0.0),
        ogeom_math::Direction::X,
        ogeom_math::Direction::Y,
        T,
    )
    .unwrap();
    let drum = ogeom_algo::make_cylinder(&mut model, main, 10.0, 40.0, T).unwrap();
    let branch = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T).unwrap();
    let joined = ogeom_bool::fuse(&mut model, &drum.shape, &branch.shape, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &joined.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // The drum plus the branch's protruding stub: the stub is the branch
    // above the saddle, and the saddle dips below the drum's crown by at
    // most the sagitta of a five-wide chord on a ten-radius circle.
    let pi = core::f64::consts::PI;
    let drum_volume = pi * 100.0 * 40.0;
    let stub_low = pi * 25.0 * 10.0;
    let stub_high = pi * 25.0 * (10.0 + 1.34);
    let measured = volume(&model, &joined.shape, 1e-3);
    assert!(
        measured > drum_volume + stub_low && measured < drum_volume + stub_high,
        "fused volume {measured} outside [{}, {}]",
        drum_volume + stub_low,
        drum_volume + stub_high
    );
}
