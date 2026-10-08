//! A small cut into a large solid costs what the cut touches, not what the
//! solid holds.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{make_box, make_cylinder};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A pin of radius 1.5 standing 20 tall from `at`.
fn pin(model: &mut Model, at: Point) -> Shape {
    let frame = Frame::new(at, Direction::Z, Direction::X, T).unwrap();
    make_cylinder(model, frame, 1.5, 20.0, T).unwrap().shape
}

/// A plate 10 thick drilled through with a `rows` by `rows` grid of holes
/// of radius 1.5, 10 apart and 10 from its sides, in one cut, with the pin
/// drilling one more hole at its corner.
fn plate_and_corner(rows: u32, faces: usize) -> (Model, Shape, Shape) {
    let mut model = Model::new();
    let size = 10.0 * f64::from(rows) + 10.0;
    let block = make_box(&mut model, Frame::WORLD, (size, size, 10.0), T)
        .unwrap()
        .shape;
    let mut pins = Vec::new();
    for i in 0..rows {
        for j in 0..rows {
            let at = Point::new(10.0 + 10.0 * f64::from(i), 10.0 + 10.0 * f64::from(j), -5.0);
            pins.push(pin(&mut model, at));
        }
    }
    let pins = model.add_compound(&pins).unwrap();
    let plate = ogeom::boolean::cut(&mut model, &block, &pins, T)
        .unwrap()
        .shape;
    assert_eq!(
        explore_unique(&model, &plate, ShapeType::Face)
            .unwrap()
            .len(),
        faces
    );
    let corner = pin(&mut model, Point::new(5.0, 5.0, -5.0));
    (model, plate, corner)
}

/// The corner hole into a plate of 582 faces costs within twice what the
/// same hole costs in a plate of 22: the faces it does not reach are set
/// aside, not split, classified or rebuilt, and what a boolean reads of the
/// whole plate (each face's edges and boxes, the classifier's preparation)
/// is kept in the model and read back while the plate stands unchanged.
/// The bound covers what the large plate still pays for every face it
/// passes through (`side`, the history copies, the shells put back
/// together), linear and memory-bound, on top of the work the cut reaches.
#[test]
#[ignore = "heavy"]
fn a_corner_hole_costs_what_it_touches_not_what_the_plate_holds() {
    let mut small = plate_and_corner(4, 22);
    let mut large = plate_and_corner(24, 582);
    // The least of many runs each, taken in turns of a few runs, the cut
    // made again each time in the model holding its plate.
    let least = |(model, plate, corner): &mut (Model, Shape, Shape)| {
        let before = explore_unique(model, plate, ShapeType::Face).unwrap().len();
        (0..6).fold(f64::INFINITY, |least, _| {
            let started = std::time::Instant::now();
            let cut = ogeom::boolean::cut(model, plate, corner, T).unwrap();
            let elapsed = started.elapsed().as_secs_f64();
            let after = explore_unique(model, &cut.shape, ShapeType::Face)
                .unwrap()
                .len();
            assert_eq!(after, before + 1);
            least.min(elapsed)
        })
    };
    // A machine busy with other work slows some runs of either side; the
    // least of each side over a round reads past that, and a round is
    // measured again, up to four, while the bound is not met.
    let mut ratios = Vec::new();
    for _ in 0..4 {
        let (mut least_small, mut least_large) = (f64::INFINITY, f64::INFINITY);
        for _ in 0..12 {
            least_small = least_small.min(least(&mut small));
            least_large = least_large.min(least(&mut large));
        }
        ratios.push((least_large, least_small));
        if least_large <= 2.0 * least_small {
            return;
        }
    }
    panic!("582 faces against 22, (large s, small s) per round: {ratios:?}");
}

/// How many faces of `result` are faces of `plate` passed through as they
/// stand.
fn passed_through(model: &Model, plate: &Shape, result: &Shape) -> usize {
    let own: ogeom::core::FastSet<_> = explore_unique(model, plate, ShapeType::Face)
        .unwrap()
        .iter()
        .map(Shape::node)
        .collect();
    explore_unique(model, result, ShapeType::Face)
        .unwrap()
        .iter()
        .filter(|face| own.contains(&face.node()))
        .count()
}

/// The vertices of the hole nearest the corner pin.
fn nearest_hole_vertices(model: &Model, plate: &Shape) -> Vec<Shape> {
    explore_unique(model, plate, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .filter(|v| {
            let at = model.node(v).unwrap().data().as_vertex().unwrap().point;
            at.x < 12.0 && at.y < 12.0 && at.x > 8.0 && at.y > 8.0
        })
        .collect()
}

/// Widen `vertex` through [`Model::widen`].
fn widen_vertex(model: &mut Model, vertex: &Shape) {
    model
        .widen(vertex, ogeom::core::Tolerance::new(1.1).unwrap())
        .unwrap();
}

/// Widen `vertex` through the node handed out for editing.
fn widen_in_node(model: &mut Model, vertex: &Shape) {
    if let Some(node) = model.node_mut(vertex)
        && let ogeom::topo::NodeData::Vertex(data) = node.data_mut()
    {
        data.tolerance = ogeom::core::Tolerance::new(1.1).unwrap();
    }
}

/// A boolean on a solid a boolean has read before reads it back as it was
/// read, and reads it as it stands after an edit in place: a vertex of the hole beside the corner pin,
/// widened until that hole's wall reaches the pin, takes the wall into the
/// cut, as the same cut into the plate edited before any boolean does.
#[test]
fn a_solid_edited_in_place_is_read_as_it_stands_by_the_next_boolean() {
    let editors: [fn(&mut Model, &Shape); 2] = [widen_vertex, widen_in_node];
    let (model, plate, corner) = plate_and_corner(8, 70);
    let mut unedited = model.clone();
    let cut = ogeom::boolean::cut(&mut unedited, &plate, &corner, T).unwrap();
    let untouched = passed_through(&unedited, &plate, &cut.shape);
    for edit in editors {
        // Read by a first cut and read back by a second, then edited.
        let mut read = model.clone();
        ogeom::boolean::cut(&mut read, &plate, &corner, T).unwrap();
        let back = ogeom::boolean::cut(&mut read, &plate, &corner, T).unwrap();
        assert_eq!(passed_through(&read, &plate, &back.shape), untouched);
        for vertex in nearest_hole_vertices(&read, &plate) {
            edit(&mut read, &vertex);
        }
        let again = ogeom::boolean::cut(&mut read, &plate, &corner, T).unwrap();
        // Edited before any boolean read it.
        let mut fresh = model.clone();
        for vertex in nearest_hole_vertices(&fresh, &plate) {
            edit(&mut fresh, &vertex);
        }
        let first = ogeom::boolean::cut(&mut fresh, &plate, &corner, T).unwrap();
        let expected = passed_through(&fresh, &plate, &first.shape);
        assert!(expected < untouched, "{expected} of {untouched}");
        assert_eq!(passed_through(&read, &plate, &again.shape), expected);
    }
}
