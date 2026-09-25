//! How big a body is: the smallest box holding it, not the carriers'.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{make_box, make_face, make_polygon, make_revolution, make_sphere, tight_bounds};
use ogeom::core::Tolerances;
use ogeom::geom::PlaneSurface;
use ogeom::math::{Axis, Direction, Frame, Plane, Point};
use ogeom::topo::Model;

const T: Tolerances = Tolerances::millimetres();

fn near(a: Point, b: Point) {
    assert!(a.distance(b) < 1e-6, "{a:?} against {b:?}");
}

#[test]
fn a_revolved_tube_reaches_its_outer_radius() {
    let mut model = Model::new();
    let pts =
        [(4.0, 0.0), (6.0, 0.0), (6.0, 10.0), (4.0, 10.0)].map(|(x, y)| Point::new(x, y, 0.0));
    let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
    let face = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Y,
    };
    let tube = make_revolution(&mut model, &face, axis, core::f64::consts::TAU, T)
        .unwrap()
        .shape;
    let b = tight_bounds(&model, &tube, T).unwrap();
    near(b.low().unwrap(), Point::new(-6.0, 0.0, -6.0));
    near(b.high().unwrap(), Point::new(6.0, 10.0, 6.0));
}

#[test]
fn a_ball_and_a_box_are_as_big_as_they_are() {
    let mut model = Model::new();
    let ball = make_sphere(&mut model, Frame::WORLD, 5.0, T).unwrap().shape;
    let b = tight_bounds(&model, &ball, T).unwrap();
    near(b.low().unwrap(), Point::new(-5.0, -5.0, -5.0));
    near(b.high().unwrap(), Point::new(5.0, 5.0, 5.0));
    let block = make_box(&mut model, Frame::WORLD, (3.0, 4.0, 5.0), T)
        .unwrap()
        .shape;
    let b = tight_bounds(&model, &block, T).unwrap();
    near(b.low().unwrap(), Point::ORIGIN);
    near(b.high().unwrap(), Point::new(3.0, 4.0, 5.0));
}
