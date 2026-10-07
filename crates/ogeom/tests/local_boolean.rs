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
/// aside, not split, classified or rebuilt.
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
    let (mut least_small, mut least_large) = (f64::INFINITY, f64::INFINITY);
    for _ in 0..12 {
        least_small = least_small.min(least(&mut small));
        least_large = least_large.min(least(&mut large));
    }
    assert!(
        least_large <= 2.0 * least_small,
        "582 faces {least_large:.5} s, 22 faces {least_small:.5} s"
    );
}
