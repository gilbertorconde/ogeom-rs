//! A rim's round taking a narrow ring beside it off whole.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{chamfer_edges, fillet_edges};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn vol(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(
        m,
        s,
        ogeom_mesh::Deflection {
            chord: 1e-4,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

fn drum(m: &mut Model) -> Shape {
    let frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
    ogeom_algo::make_cylinder(m, frame, 10.0, 10.0, T)
        .unwrap()
        .shape
}

/// The circular edge of `s` of radius `r` at height `z`.
fn circle(m: &Model, s: &Shape, r: f64, z: f64) -> Shape {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let Ok(Some((c, (a, b)))) = ogeom_algo::placed_edge_curve(m, e, T) else {
                return false;
            };
            [a, 0.5 * (a + b)].iter().all(|&t| {
                let p = c.point_at(t, T).unwrap();
                (p.z - z).abs() < 1e-9 && (p.x.hypot(p.y) - r).abs() < 1e-9
            })
        })
        .expect("the circle")
}

/// The drum's top rim rounded at 3 on its own: what the ring's round is.
fn plain() -> f64 {
    let mut m = Model::new();
    let d = drum(&mut m);
    let rim = circle(&m, &d, 10.0, 10.0);
    let built = fillet_edges(&mut m, &d, &[rim], 3.0, T).unwrap();
    vol(&m, &built.shape)
}

/// Round `edge` of `ringed` at 3: valid, no face of the ring left, and the
/// plain drum's round.
fn takes_the_ring(m: &mut Model, ringed: &Shape, edge: &Shape, ring_faces: usize) {
    let faces = |m: &Model, s: &Shape| explore_unique(m, s, ShapeType::Face).unwrap().len();
    assert_eq!(faces(m, ringed), 3 + ring_faces);
    let built = fillet_edges(m, ringed, std::slice::from_ref(edge), 3.0, T).unwrap();
    let diagnosis = ogeom_algo::check(m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    // Cap, wall, floor and the round: the ring is gone.
    assert_eq!(faces(m, &built.shape), 4);
    let (v, w) = (vol(m, &built.shape), plain());
    assert!((v - w).abs() <= 1e-6 * w, "{v} against {w}");
}

#[test]
fn a_chamfered_rim_rounds_larger_taking_the_chamfer() {
    let mut m = Model::new();
    let d = drum(&mut m);
    let rim = circle(&m, &d, 10.0, 10.0);
    let chamfered = chamfer_edges(&mut m, &d, &[rim], 0.5, T).unwrap().shape;
    // The edge between the cap and the chamfer, and the one between the
    // chamfer and the wall: either rounds the same.
    for (r, z) in [(9.5, 10.0), (10.0, 9.5)] {
        let mut mm = m.clone();
        let edge = circle(&mm, &chamfered, r, z);
        takes_the_ring(&mut mm, &chamfered, &edge, 1);
    }
}

#[test]
fn a_small_round_rounds_larger_taking_the_small_one() {
    let mut m = Model::new();
    let d = drum(&mut m);
    let rim = circle(&m, &d, 10.0, 10.0);
    let rounded = fillet_edges(&mut m, &d, &[rim], 0.5, T).unwrap().shape;
    let edge = circle(&m, &rounded, 9.5, 10.0);
    takes_the_ring(&mut m, &rounded, &edge, 1);
}
