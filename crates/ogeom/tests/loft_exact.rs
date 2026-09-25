//! A smooth loft holds its sections exactly where they are exactly
//! representable: a corner stays a corner and a circle a circle.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_edge, make_polygon, make_wire, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::CircleCurve;
use ogeom::math::{Circle, Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn square(model: &mut Model, half: f64, z: f64) -> Shape {
    let pts = [(-half, -half), (half, -half), (half, half), (-half, half)]
        .map(|(x, y)| Point::new(x, y, z));
    make_polygon(model, &pts, true, T).unwrap().shape
}

fn circle(model: &mut Model, r: f64, z: f64) -> Shape {
    let frame = Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X, T).unwrap();
    let c = Circle::new(frame, r, T).unwrap();
    let e = make_edge(
        model,
        CircleCurve::new(c).into(),
        (0.0, core::f64::consts::TAU),
        T,
    )
    .unwrap()
    .shape;
    make_wire(model, &[e], T).unwrap().shape
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    assert!(check(model, shape, T).unwrap().is_valid());
    let fine = Deflection {
        chord: 1e-3,
        angular: 0.02,
        ..Deflection::default()
    };
    volume_properties(model, shape, fine, T).unwrap().mass
}

fn loft(model: &mut Model, sections: &[Shape]) -> Shape {
    ogeom::offset::make_loft_skinned(model, sections, 1e-3, T)
        .unwrap()
        .shape
}

#[test]
fn two_squares_loft_to_a_prism() {
    let mut model = Model::new();
    let (a, b) = (square(&mut model, 5.0, 0.0), square(&mut model, 5.0, 15.0));
    let solid = loft(&mut model, &[a, b]);
    let v = volume(&model, &solid);
    assert!((v - 1500.0).abs() < 1500.0 * 1e-4, "{v}");
}

#[test]
fn two_circles_loft_to_a_drum() {
    let mut model = Model::new();
    let (a, b) = (circle(&mut model, 5.0, 0.0), circle(&mut model, 5.0, 15.0));
    let want = core::f64::consts::PI * 25.0 * 15.0;
    let solid = loft(&mut model, &[a, b]);
    let v = volume(&model, &solid);
    assert!((v - want).abs() < want * 1e-4, "{v} against {want}");
}

#[test]
fn three_squares_loft_to_a_prism() {
    let mut model = Model::new();
    let sections: Vec<Shape> = [0.0, 7.0, 15.0]
        .iter()
        .map(|z| square(&mut model, 5.0, *z))
        .collect();
    let solid = loft(&mut model, &sections);
    let v = volume(&model, &solid);
    assert!((v - 1500.0).abs() < 1500.0 * 1e-4, "{v}");
}

#[test]
fn three_equal_circles_loft_to_a_drum() {
    let mut model = Model::new();
    let sections: Vec<Shape> = [0.0, 6.0, 15.0]
        .iter()
        .map(|z| circle(&mut model, 5.0, *z))
        .collect();
    let want = core::f64::consts::PI * 25.0 * 15.0;
    let solid = loft(&mut model, &sections);
    let v = volume(&model, &solid);
    assert!((v - want).abs() < want * 1e-4, "{v} against {want}");
}

#[test]
fn circles_of_three_radii_loft_round() {
    let mut model = Model::new();
    let sections = vec![
        circle(&mut model, 5.0, 0.0),
        circle(&mut model, 7.0, 5.0),
        circle(&mut model, 4.0, 10.0),
    ];
    let solid = loft(&mut model, &sections);
    let v = volume(&model, &solid);
    // Between the cones through the radii and the drum of the largest.
    assert!(
        v > core::f64::consts::PI * 16.0 * 10.0 && v < core::f64::consts::PI * 49.0 * 10.0,
        "{v}"
    );
}
