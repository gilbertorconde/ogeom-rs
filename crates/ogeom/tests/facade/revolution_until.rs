//! A revolution that stops on a surface: each point of the profile turns
//! until its circle first meets the limit, at a different angle for each.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::PlaneSurface;
use ogeom::math::{Axis, Direction, Frame, Plane, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The rectangle x in [5, 10], z in [0, 10] on the plane y = 0 turned about
/// Z until it meets the face x = -2 of a block beside the axis: a point at
/// radius r turns through acos(-2 / r), so the solid is
/// 10 * integral over r in [5, 10] of acos(-2 / r) r dr.
#[test]
fn a_revolution_stops_on_a_wall_beside_its_axis() {
    let mut model = Model::new();
    let pts =
        [(5.0, 0.0), (10.0, 0.0), (10.0, 10.0), (5.0, 10.0)].map(|(x, z)| Point::new(x, 0.0, z));
    let wire = ogeom::algo::make_polygon(&mut model, &pts, true, T)
        .unwrap()
        .shape;
    let plane = Plane::new(Frame::new(Point::ORIGIN, -Direction::Y, Direction::X, T).unwrap());
    let profile = ogeom::algo::make_face(&mut model, PlaneSurface::new(plane).into(), &[wire], T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(-20.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let block = ogeom::algo::make_box(&mut model, frame, (18.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let wall = explore_unique(&model, &block, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let b = ogeom::algo::shape_bounds(&model, f, T).unwrap();
            (b.low().unwrap().x + 2.0).abs() < 1e-6 && (b.high().unwrap().x + 2.0).abs() < 1e-6
        })
        .expect("the wall at x = -2");
    let axis = Axis::new(Point::ORIGIN, Direction::Z);
    let turned = ogeom::offset::make_revolution_until(&mut model, &profile, axis, &wall, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &turned, T).unwrap().is_valid());
    let fine = Deflection::with_chord(1e-3).unwrap();
    let volume = |model: &Model, shape| {
        ogeom::algo::volume_properties(model, shape, fine, T)
            .unwrap()
            .mass
    };
    let n = 20_000;
    let want: f64 = (0..n)
        .map(|i| {
            let r = 5.0 + 5.0 * (f64::from(i) + 0.5) / f64::from(n);
            (-2.0 / r).acos() * r
        })
        .sum::<f64>()
        * 5.0
        / f64::from(n)
        * 10.0;
    let v = volume(&model, &turned);
    assert!((v - want).abs() < want * 2e-3, "{v} against {want}");
    let joined = ogeom::boolean::fuse(&mut model, &block, &turned, T)
        .unwrap()
        .shape;
    let both = volume(&model, &joined);
    assert!((both - (3600.0 + want)).abs() < want * 2e-3, "{both}");
}
