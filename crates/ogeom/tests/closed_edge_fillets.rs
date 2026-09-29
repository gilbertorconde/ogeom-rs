//! Fillet chains holding closed edges: both rims of a drum, and a rim the
//! solid keeps as two halves.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{make_box, make_cylinder, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::Curve;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{EdgeRepr, Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// The edges on circles, with the height of each circle's centre.
fn circular_edges(model: &Model, shape: &Shape) -> Vec<(Shape, f64)> {
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter_map(|edge| {
            let data = model.node(&edge)?.data().as_edge()?;
            let Some(EdgeRepr::Curve3d { curve, .. }) = data.curve3d() else {
                return None;
            };
            match model.geometry().curve(*curve)? {
                Curve::Circle(c) => Some((edge.clone(), c.circle().centre().z)),
                _ => None,
            }
        })
        .collect()
}

/// A convex rim of radius `r` on a drum's end, rounded by `f`: the ring the
/// ball leaves behind, by Pappus.
fn rim_removed(r: f64, f: f64) -> f64 {
    // The corner's cross-section: an f x f square less a quarter disc.
    let area = f * f * (1.0 - PI / 4.0);
    // Its centroid, measured in from the square's outer corner.
    let inset = f * (10.0 - 3.0 * PI) / (3.0 * (4.0 - PI));
    2.0 * PI * (r - inset) * area
}

/// A drum with both rims in one chain: each rounds.
#[test]
fn both_rims_of_a_drum_round_in_one_chain() {
    let mut model = Model::new();
    let drum = make_cylinder(&mut model, Frame::WORLD, 10.0, 20.0, T)
        .unwrap()
        .shape;
    let rims: Vec<Shape> = circular_edges(&model, &drum)
        .into_iter()
        .map(|(e, _)| e)
        .collect();
    assert_eq!(rims.len(), 2);
    let rounded = ogeom::fillet::fillet_edges(&mut model, &drum, &rims, 1.0, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &rounded, T).unwrap().is_valid());
    let want = PI * 100.0 * 20.0 - 2.0 * rim_removed(10.0, 1.0);
    let v = volume(&model, &rounded);
    assert!((v - want).abs() < want * 1e-4, "{v} against {want}");
}

/// A blind hole's top rim, held as two half circles: both halves in the
/// chain round the whole rim once, and either half alone rounds it the
/// same way.
#[test]
fn a_rim_in_two_halves_rounds_as_one_loop() {
    let mut model = Model::new();
    let frame = Frame::new(Point::new(-20.0, -20.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let block = make_box(&mut model, frame, (40.0, 40.0, 10.0), T)
        .unwrap()
        .shape;
    let seat = Frame::new(Point::new(0.0, 0.0, 5.0), Direction::Z, Direction::X, T).unwrap();
    let bore = make_cylinder(&mut model, seat, 5.0, 10.0, T).unwrap().shape;
    let holed = ogeom::boolean::cut(&mut model, &block, &bore, T)
        .unwrap()
        .shape;
    let rim: Vec<Shape> = circular_edges(&model, &holed)
        .into_iter()
        .filter(|(_, z)| (z - 10.0).abs() < 1e-9)
        .map(|(e, _)| e)
        .collect();
    assert_eq!(rim.len(), 2, "the rim is two halves");
    let both = ogeom::fillet::fillet_edges(&mut model, &holed, &rim, 1.0, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &both, T).unwrap().is_valid());
    let one = ogeom::fillet::fillet_edges(&mut model, &holed, &rim[..1], 1.0, T)
        .unwrap()
        .shape;
    let (a, b) = (volume(&model, &both), volume(&model, &one));
    assert!((a - b).abs() < a * 1e-9, "{a} against {b}");
    // The rim is convex: the ball takes off the corner's ring, whose
    // centroid stands outside the bore by the corner's inset.
    let before = volume(&model, &holed);
    let area = 1.0 - PI / 4.0;
    let inset = (10.0 - 3.0 * PI) / (3.0 * (4.0 - PI));
    let removed = 2.0 * PI * (5.0 + inset) * area;
    assert!(
        (before - a - removed).abs() < removed * 1e-3,
        "{a} from {before}"
    );
}
