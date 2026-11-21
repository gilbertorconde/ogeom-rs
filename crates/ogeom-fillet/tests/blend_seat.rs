//! The seat reading on its own: a face's reach for a round of two
//! distances, the other distance small enough that only that face stops
//! it, read without building.
#![allow(clippy::unwrap_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{BlendSizing, RoundProfile, blend_seat};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_plate_s_top_and_front_reaches_are_its_depth_and_height() {
    let mut m = ogeom_topo::Model::new();
    let frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
    let s = ogeom_algo::make_box(&mut m, frame, (20.0, 4.0, 20.0), T)
        .unwrap()
        .shape;
    let edge = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap();
            [a, b].iter().all(|&t| {
                let p = c.point_at(t, T).unwrap();
                p.y.abs() < 1e-9 && (p.z - 20.0).abs() < 1e-9
            })
        })
        .unwrap();
    let face_at = |on: fn(Point) -> bool| {
        explore_unique(&m, &s, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                let b = ogeom_algo::shape_bounds(&m, f, T).unwrap();
                on(b.low().unwrap()) && on(b.high().unwrap())
            })
            .unwrap()
    };
    let top = face_at(|p| (p.z - 20.0).abs() < 1e-6);
    let front = face_at(|p| p.y.abs() < 1e-6);
    for (face, reach) in [(&top, 4.0), (&front, 20.0)] {
        let sizing = BlendSizing::RoundDistances {
            edge: &edge,
            face,
            on_face: 1.0,
            on_other: 1e-5,
            profile: RoundProfile::Elliptical,
        };
        let seat = blend_seat(&m, &s, &sizing, 40.0, T).unwrap().unwrap();
        assert!((seat - reach).abs() < 1e-6, "{seat} against {reach}");
    }
}
