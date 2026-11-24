//! A solid whose boundary passes through itself, resolved by its union.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A box by its low corner and its size.
type Placed = ((f64, f64, f64), (f64, f64, f64));

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom_algo::volume_properties(model, shape, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// Two boxes, `low` corner and size each, joined as one solid of two
/// shells: where they overlap, the boundary encloses the material twice.
fn two_boxes_as_one(model: &mut Model, boxes: [Placed; 2]) -> Shape {
    let shells: Vec<Shape> = boxes
        .iter()
        .map(|&((x, y, z), size)| {
            let frame = Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap();
            let solid = ogeom_algo::make_box(model, frame, size, T).unwrap().shape;
            explore_unique(model, &solid, ShapeType::Shell).unwrap()[0].clone()
        })
        .collect();
    model.add_solid(&shells).unwrap()
}

fn resolved(boxes: [Placed; 2], union: f64) {
    let mut model = Model::new();
    let solid = two_boxes_as_one(&mut model, boxes);
    let before = ogeom_algo::check(&model, &solid, T).unwrap();
    assert!(!before.is_valid(), "the two boxes cross each other");
    let built = ogeom_bool::self_union(&mut model, &solid, T).unwrap();
    let after = ogeom_algo::check(&model, &built.shape, T).unwrap();
    assert!(after.is_valid(), "{:?}", after.problems);
    assert_eq!(model.kind_of(&built.shape).unwrap(), ShapeType::Solid);
    assert_eq!(
        explore_unique(&model, &built.shape, ShapeType::Shell)
            .unwrap()
            .len(),
        1
    );
    let measured = volume(&model, &built.shape);
    assert!(
        (measured - union).abs() <= 1e-9 * union,
        "volume {measured}, the union's {union}"
    );
}

#[test]
fn two_boxes_crossing_in_one_solid_resolve_to_their_union() {
    // 8 + 4 - the 1 x 1 x 1.5 overlap.
    resolved(
        [
            ((0.0, 0.0, 0.0), (2.0, 2.0, 2.0)),
            ((1.0, 0.5, 0.5), (2.0, 1.0, 2.0)),
        ],
        10.5,
    );
}

#[test]
fn two_boxes_sharing_a_floor_and_a_roof_resolve_to_their_union() {
    // The second box's floor and roof lie in the first one's, facing the
    // same way, over the overlap: 8 + 4 - the 1 x 1 x 2 overlap.
    resolved(
        [
            ((0.0, 0.0, 0.0), (2.0, 2.0, 2.0)),
            ((1.0, 0.5, 0.0), (2.0, 1.0, 2.0)),
        ],
        10.0,
    );
}

/// The prism over a polygon in the plane z = 0, one high.
fn prism(model: &mut Model, corners: &[(f64, f64)]) -> Shape {
    let points: Vec<Point> = corners
        .iter()
        .map(|&(x, y)| Point::new(x, y, 0.0))
        .collect();
    let wire = ogeom_algo::make_polygon(model, &points, true, T)
        .unwrap()
        .shape;
    let plane = ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(Frame::WORLD));
    let face = ogeom_algo::make_face(model, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    ogeom_algo::make_prism(model, &face, ogeom_math::Vector::new(0.0, 0.0, 1.0), T)
        .unwrap()
        .shape
}

/// A strip one wide wound round three sides of a square and down its
/// fourth, its last arm running back into its first: the arms
/// [0, 4] x [0, 1] and [0.2, 1.2] x [0.5, 2] overlap over a half. Built as
/// two prisms sewn where the last arm leaves the top one: one closed
/// shell, its walls crossing and its floor and roof lying twice over the
/// overlap, as a file's offset wall thicker than its offset is. The wall
/// the two prisms share is taken out of both, or, where `joint`, kept in
/// both, facing either way.
fn wound_strip_resolves(joint: bool) {
    let mut model = Model::new();
    let bends = prism(
        &mut model,
        &[
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 3.0),
            (0.2, 3.0),
            (0.2, 2.0),
            (1.2, 2.0),
            (3.0, 2.0),
            (3.0, 1.0),
            (0.0, 1.0),
        ],
    );
    let last = prism(
        &mut model,
        &[(0.2, 0.5), (1.2, 0.5), (1.2, 2.0), (0.2, 2.0)],
    );
    let mut faces = Vec::new();
    for part in [&bends, &last] {
        for face in explore_unique(&model, part, ShapeType::Face).unwrap() {
            let bound = ogeom_algo::face_bounds(&model, &face).unwrap();
            let (low, high) = (bound.low().unwrap(), bound.high().unwrap());
            // The bound is widened alike each way, so its middle is exact.
            let shared = high.y - low.y < 0.5
                && ((low.y + high.y) / 2.0 - 2.0).abs() < 1e-9
                && ((low.x + high.x) / 2.0 - 0.7).abs() < 1e-9;
            if joint || !shared {
                faces.push(face);
            }
        }
    }
    let sewn = ogeom_algo::sew(&mut model, &faces, T).unwrap();
    assert_eq!(sewn.shells.len(), 1);
    assert!(sewn.free_edges.is_empty());
    let solid = model.add_solid(&sewn.shells).unwrap();
    // The strip's area counts the overlap twice: 4 + 2 + 2.8 + 1.5.
    assert!((volume(&model, &solid) - 10.3).abs() <= 1e-9);
    assert!(
        !ogeom_algo::self_crossings(&model, &solid, T)
            .unwrap()
            .is_empty()
    );

    let built = ogeom_bool::self_union(&mut model, &solid, T).unwrap();
    let after = ogeom_algo::check(&model, &built.shape, T).unwrap();
    assert!(after.is_valid(), "{:?}", after.problems);
    assert_eq!(
        explore_unique(&model, &built.shape, ShapeType::Shell)
            .unwrap()
            .len(),
        1
    );
    let measured = volume(&model, &built.shape);
    assert!(
        (measured - 9.8).abs() <= 1e-9,
        "volume {measured}, the union's 9.8"
    );
}

#[test]
fn a_strip_wound_back_through_itself_in_one_shell_resolves_to_its_union() {
    wound_strip_resolves(false);
}

#[test]
fn a_wound_strip_with_a_wall_inside_twice_resolves_to_its_union() {
    // The wall lying twice, facing either way, has the material on both
    // sides of it and goes, with the faces either side of it it shares
    // edges with kept.
    wound_strip_resolves(true);
}

#[test]
fn a_solid_whose_boundary_does_not_cross_itself_comes_back_as_it_is() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T)
        .unwrap()
        .shape;
    let built = ogeom_bool::self_union(&mut model, &solid, T).unwrap();
    assert!(built.shape.is_same(&solid));
    assert!(built.history.is_empty());
}
