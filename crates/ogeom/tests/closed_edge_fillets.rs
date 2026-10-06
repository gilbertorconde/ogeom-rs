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

/// Two crossing bores: one of radius 4 along x at y = z = 10, one of
/// radius 3 along y at x = 10, z = 12, through a 20 mm cube. The narrower
/// rises out of the wider's top, so the two walls meet in one closed loop,
/// a convex edge the solid keeps in pieces; the narrower bore's wall is
/// the second cut's tool face, presented the other way round.
fn crossed_bores(model: &mut Model) -> Shape {
    let block = make_box(model, Frame::WORLD, (20.0, 20.0, 20.0), T)
        .unwrap()
        .shape;
    let along_x = Frame::new(Point::new(-1.0, 10.0, 10.0), Direction::X, Direction::Y, T).unwrap();
    let along_y = Frame::new(Point::new(10.0, -1.0, 12.0), Direction::Y, Direction::Z, T).unwrap();
    let wide = make_cylinder(model, along_x, 4.0, 22.0, T).unwrap().shape;
    let narrow = make_cylinder(model, along_y, 3.0, 22.0, T).unwrap().shape;
    let once = ogeom::boolean::cut(model, &block, &wide, T).unwrap().shape;
    ogeom::boolean::cut(model, &once, &narrow, T).unwrap().shape
}

/// What a small ball takes off the crossed bores' loop, over its radius
/// squared, in the limit: the loop's integral of `tan(t/2) - t/2`, `t`
/// the turn between the walls' normals, the cross-section a ball leaves
/// in a corner of that turn. The loop is traced in closed form, by the
/// narrower bore's angle, on each side of the wider's axis.
fn crossed_loop_share() -> f64 {
    let top = (2.0_f64 / 3.0).asin();
    let (from, to) = (PI - top, top + 2.0 * PI);
    let steps = 20_000;
    let mut share = 0.0;
    for side in [1.0, -1.0] {
        let at = |a: f64| {
            let z = 3.0f64.mul_add(a.sin(), 12.0);
            let y = side * (16.0 - (z - 10.0).powi(2)).max(0.0).sqrt() + 10.0;
            Point::new(3.0f64.mul_add(a.cos(), 10.0), y, z)
        };
        for i in 0..steps {
            let a0 = from + (to - from) * f64::from(i) / f64::from(steps);
            let a1 = from + (to - from) * f64::from(i + 1) / f64::from(steps);
            let p = at(f64::midpoint(a0, a1));
            let wide = Point::new(p.x, 10.0, 10.0) - p;
            let narrow = Point::new(10.0, p.y, 12.0) - p;
            let turn = (wide.dot(narrow) / 12.0).clamp(-1.0, 1.0).acos();
            share += ((turn / 2.0).tan() - turn / 2.0) * at(a0).distance(at(a1));
        }
    }
    share
}

/// The loop where the crossed bores meet rounds from any piece of it:
/// the ball rolls round the whole tangent loop, and what it takes off at
/// radii 0.4 and 0.2, extrapolated to the ball's vanishing, is the loop's
/// share in closed form. Rounding starts from the short piece by the
/// wider bore's seam, where the narrower wall's trim turns back.
#[test]
fn the_loop_of_two_crossed_bores_rounds_to_its_share() {
    let mut model = Model::new();
    let part = crossed_bores(&mut model);
    let seam = Point::new(7.763_929_322, 14.0, 10.0);
    let piece = explore_unique(&model, &part, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let Some((a, b)) = ogeom::algo::edge_vertices(&model, e).unwrap() else {
                return false;
            };
            let at = |v: &Shape| model.node(v).unwrap().data().as_vertex().unwrap().point;
            let ends = [at(&a), at(&b)];
            ends.iter().any(|p| p.distance(seam) < 1e-6)
                && ends.iter().any(|p| p.y > 13.0 && p.z > 10.5)
        })
        .expect("the piece of the loop at the wider bore's seam");
    let before = volume(&model, &part);
    let mut per_square = Vec::new();
    for r in [0.4, 0.2] {
        let mut copy = model.clone();
        let rounded =
            ogeom::fillet::fillet_edges(&mut copy, &part, std::slice::from_ref(&piece), r, T)
                .unwrap()
                .shape;
        let diagnosis = ogeom::algo::check(&copy, &rounded, T).unwrap();
        assert!(diagnosis.is_valid(), "r {r}: {:?}", diagnosis.problems);
        per_square.push((before - volume(&copy, &rounded)) / (r * r));
    }
    // What a ball takes off over its radius squared is the share plus a
    // term in the radius and one in its square: twice the second less the
    // first leaves the share less half the square's term, about two
    // thousandths of the share at these radii.
    let extrapolated = 2.0f64.mul_add(per_square[1], -per_square[0]);
    let want = crossed_loop_share();
    assert!(
        (extrapolated - want).abs() < want * 5e-3,
        "{per_square:?} extrapolate to {extrapolated}, the loop's share is {want}"
    );
}
