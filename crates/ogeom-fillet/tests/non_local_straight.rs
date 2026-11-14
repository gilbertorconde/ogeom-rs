//! A straight round along a cylinder's ruling reaching a face it does not
//! end at, through the blend engine's one boolean (ADR 0001, M5a and M5c):
//! the strip closed into its tool, capped on its end faces, and cut.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

/// A drum of radius 5, 10 tall, cut flat at `x = 3` (its part `x <= 3`
/// kept), with a hole 0.1 by 0.15 drilled 3 down from its top within where
/// the round of its straight edge at `(3, 4)` at 2 takes: the round takes
/// the hole with it. Through the one boolean, valid, the edge generating
/// its strip, the volume the flat drum's less the spandrel's section along
/// the 10: the corner `P = (3, 4)`, the round's foot on the flat `L = (3, 2
/// sqrt 2)` and on the wall `C = 5/3 (1, 2 sqrt 2)`, the triangle `PLC` with
/// the wall's segment over `CP`, less the ball's over `LC`.
#[test]
fn a_flat_drum_s_straight_edge_rounds_past_a_hole_within_its_spandrel() {
    let mut m = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut m, at(0.0, 0.0, 0.0), 5.0, 10.0, T)
        .unwrap()
        .shape;
    let beyond = ogeom_algo::make_box(&mut m, at(3.0, -6.0, -1.0), (3.0, 12.0, 12.0), T)
        .unwrap()
        .shape;
    let flat = ogeom_bool::cut(&mut m, &drum, &beyond, T).unwrap().shape;
    let hole = ogeom_algo::make_box(&mut m, at(2.85, 3.6, 7.0), (0.1, 0.15, 4.0), T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(&mut m, &flat, &hole, T).unwrap().shape;
    let edge = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let bb = ogeom_algo::shape_bounds(&m, e, T).unwrap();
            let (l, h) = (bb.low().unwrap(), bb.high().unwrap());
            (l.x - 3.0).abs() < 1e-6
                && (h.x - 3.0).abs() < 1e-6
                && (l.y - 4.0).abs() < 1e-6
                && (h.y - 4.0).abs() < 1e-6
                && (h.z - l.z - 10.0).abs() < 1e-5
        })
        .expect("the straight edge at (3, 4)");
    let built =
        ogeom_fillet::fillet_edges(&mut m, &s, std::slice::from_ref(&edge), 2.0, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    let r2 = 2.0_f64.sqrt();
    let (p, l, c) = ((3.0, 4.0), (3.0, 2.0 * r2), (5.0 / 3.0, 10.0 * r2 / 3.0));
    let triangle = 0.5 * ((l.0 - p.0) * (c.1 - p.1) - (c.0 - p.0) * (l.1 - p.1)).abs();
    let segment = |radius: f64, cos: f64| {
        let phi = cos.clamp(-1.0, 1.0).acos();
        0.5 * radius * radius * (phi - phi.sin())
    };
    let wall = segment(5.0, (c.0 * p.0 + c.1 * p.1) / 25.0);
    let ball = segment(2.0, 1.0 / 3.0);
    let spandrel = triangle + wall - ball;
    let flat_area = 25.0 * core::f64::consts::PI - 25.0 * 0.6_f64.acos() + 12.0;
    let want = (flat_area - spandrel) * 10.0;
    let got = ogeom_algo::volume_properties(&m, &built.shape, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((got - want).abs() < 1e-6, "{got} against {want}");
    assert!(!built.history.generated(&edge).is_empty(), "the strip");
}
