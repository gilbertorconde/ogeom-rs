//! The divide family: faces cut along their own parameter lines, the
//! solid unchanged in shape and still closed.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::{Continuity, SurfaceGeometry};
use ogeom::heal::IsoLine;
use ogeom::math::Frame;
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn faces(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Face).unwrap()
}

/// The solid is valid, closed and encloses `volume`.
fn holds(model: &Model, shape: &Shape, volume: f64, band: f64) {
    let diagnosis = ogeom::algo::check(model, shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    // Fine enough that a spline wall's mesh is within the band.
    let fine = Deflection {
        chord: 1e-3,
        angular: 0.02,
        ..Deflection::default()
    };
    let measured = ogeom::algo::volume_properties(model, shape, fine, T)
        .unwrap()
        .mass;
    assert!(
        (measured - volume).abs() <= volume * band,
        "{measured} against {volume}"
    );
}

fn surface_of<'a>(model: &'a Model, face: &Shape) -> &'a SurfaceGeometry {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        panic!("not a face");
    };
    model.geometry().surface(data.surface).unwrap()
}

#[test]
fn a_box_face_cut_across_leaves_the_box() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let face = faces(&model, &block)[0].clone();
    // The middle of the face's own chart, whichever way it runs.
    let divided = ogeom::heal::divide_face(&mut model, &block, &face, IsoLine::U(5.0), T)
        .or_else(|_| ogeom::heal::divide_face(&mut model, &block, &face, IsoLine::U(-5.0), T))
        .unwrap()
        .shape;
    assert_eq!(faces(&model, &divided).len(), 7);
    holds(&model, &divided, 1000.0, 1e-9);
}

#[test]
fn a_drum_divides_into_quarter_turns() {
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T)
        .unwrap()
        .shape;
    let divided = ogeom::heal::divide_by_angle(&mut model, &drum, core::f64::consts::FRAC_PI_2, T)
        .unwrap()
        .shape;
    let walls = faces(&model, &divided)
        .into_iter()
        .filter(|f| matches!(surface_of(&model, f), SurfaceGeometry::Cylinder(_)))
        .count();
    assert_eq!(walls, 4, "four quarter walls");
    holds(&model, &divided, core::f64::consts::PI * 4.0 * 5.0, 1e-3);
}

#[test]
fn a_ball_divides_by_angle_both_ways() {
    let mut model = Model::new();
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 3.0, T)
        .unwrap()
        .shape;
    let divided = ogeom::heal::divide_by_angle(&mut model, &ball, core::f64::consts::FRAC_PI_2, T)
        .unwrap()
        .shape;
    assert_eq!(
        faces(&model, &divided).len(),
        8,
        "four lunes, halved at the equator"
    );
    holds(
        &model,
        &divided,
        4.0 / 3.0 * core::f64::consts::PI * 27.0,
        1e-3,
    );
}

#[test]
fn a_converted_drum_divides_at_its_weak_knots() {
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T)
        .unwrap()
        .shape;
    let nurbs = ogeom::algo::to_nurbs(&mut model, &drum, T).unwrap().shape;
    let before = faces(&model, &nurbs).len();
    holds(&model, &nurbs, core::f64::consts::PI * 4.0 * 5.0, 1e-3);
    let divided = ogeom::heal::divide_by_continuity(&mut model, &nurbs, Continuity::C1, T)
        .unwrap()
        .shape;
    assert!(
        faces(&model, &divided).len() > before,
        "the circle's double knots cut the wall"
    );
    holds(&model, &divided, core::f64::consts::PI * 4.0 * 5.0, 1e-3);
    // Every piece is now smooth throughout.
    let again = ogeom::heal::divide_by_continuity(&mut model, &divided, Continuity::C1, T)
        .unwrap()
        .shape;
    assert_eq!(faces(&model, &again).len(), faces(&model, &divided).len());
}

#[test]
fn every_piece_of_a_bezier_drum_is_one_span() {
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T)
        .unwrap()
        .shape;
    let bezier = ogeom::heal::to_bezier(&mut model, &drum, T).unwrap().shape;
    for face in faces(&model, &bezier) {
        let SurfaceGeometry::BSpline(s) = surface_of(&model, &face) else {
            panic!("a converted face is a spline");
        };
        assert_eq!(s.u_knots().distinct().len(), 2, "one span in u");
        assert_eq!(s.v_knots().distinct().len(), 2, "one span in v");
    }
    holds(&model, &bezier, core::f64::consts::PI * 4.0 * 5.0, 1e-3);
}

#[test]
fn a_box_divides_until_no_face_is_larger_than_asked() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let divided = ogeom::heal::divide_by_area(&mut model, &block, 30.0, T)
        .unwrap()
        .shape;
    let all = faces(&model, &divided);
    assert_eq!(all.len(), 24, "each side in four");
    for face in &all {
        let area = ogeom::algo::surface_properties(&model, face, Deflection::default(), T)
            .unwrap()
            .mass;
        assert!(area <= 30.0, "{area}");
    }
    holds(&model, &divided, 1000.0, 1e-9);
}

#[test]
fn a_plate_with_a_hole_divides_through_and_beside_it() {
    let mut model = Model::new();
    let plate = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 10.0, 2.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(
        ogeom::math::Point::new(6.0, 5.0, -1.0),
        ogeom::math::Direction::Z,
        ogeom::math::Direction::X,
        T,
    )
    .unwrap();
    let post = ogeom::algo::make_cylinder(&mut model, at, 2.0, 4.0, T)
        .unwrap()
        .shape;
    let holed = ogeom::boolean::cut(&mut model, &plate, &post, T)
        .unwrap()
        .shape;
    let volume = 400.0 - core::f64::consts::PI * 4.0 * 2.0;
    let before = faces(&model, &holed).len();
    // Every face quartered by area crosses the hole or passes beside it.
    let divided = ogeom::heal::divide_by_area(&mut model, &holed, 60.0, T)
        .unwrap()
        .shape;
    assert!(faces(&model, &divided).len() > before);
    holds(&model, &divided, volume, 1e-3);
    // And the hole's wall by angle.
    let turned = ogeom::heal::divide_by_angle(&mut model, &divided, 1.0, T)
        .unwrap()
        .shape;
    holds(&model, &turned, volume, 1e-3);
}

#[test]
fn a_torus_divides_by_angle_both_ways() {
    let mut model = Model::new();
    let ring = ogeom::algo::make_torus(&mut model, Frame::WORLD, 5.0, 1.0, T)
        .unwrap()
        .shape;
    let divided = ogeom::heal::divide_by_angle(&mut model, &ring, core::f64::consts::PI, T)
        .unwrap()
        .shape;
    assert_eq!(faces(&model, &divided).len(), 4);
    holds(
        &model,
        &divided,
        2.0 * core::f64::consts::PI.powi(2) * 5.0,
        1e-3,
    );
}
