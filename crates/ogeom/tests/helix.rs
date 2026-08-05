//! A thread from the new §3 helix: the exact curve swept into a pipe by the
//! skinned sweep, its volume answering to the closed-form arc length.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::{Curve, Curve3d, HelixCurve};
use ogeom::math::Frame;
use ogeom::mesh::Deflection;
use ogeom::topo::Model;

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_helical_pipe_is_a_spring_that_measures() {
    let mut model = Model::new();
    // Two turns, radius 10, pitch 8 — wide enough that the tube nowhere
    // approaches itself and the pipe volume is the tube formula.
    let helix = HelixCurve::new(Frame::WORLD, 10.0, 8.0, 2.0).unwrap();
    let length = helix.arc_length(0.0, 2.0 * core::f64::consts::TAU);
    let domain = Curve3d::domain(&helix);
    let spine = ogeom::algo::make_edge(&mut model, Curve::Helix(helix), domain, T)
        .unwrap()
        .shape;
    let spring = ogeom::offset::make_pipe_skinned(&mut model, &spine, 1.5, 5e-3, T)
        .unwrap()
        .shape;

    // The tube is thin against the coil: measure at a chord an order under
    // the tube radius so the inscribed mesh's undercut stays inside the
    // volume check's own budget.
    let fine = Deflection::with_chord(1e-2).unwrap();
    let volume = ogeom::algo::volume_properties(&model, &spring, fine, T)
        .unwrap()
        .mass;
    // A tube around a constant-speed spine encloses length times section
    // area, exactly, when it stays clear of itself; the fit and the mesh
    // spend the rest of the budget.
    let expected = length * core::f64::consts::PI * 1.5 * 1.5;
    assert!(
        (volume - expected).abs() / expected < 0.01,
        "spring volume {volume} against tube formula {expected}"
    );
}
