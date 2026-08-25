//! A written subset is a subset: `write(model, roots)` carries the closure
//! of `roots`, not the model (issue #16).
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::io::native::{WriteOptions, read, write};
use ogeom::math::{Direction, Frame, Point};
use ogeom::topo::Model;

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &ogeom::topo::Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, ogeom::mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// One box's snapshot out of a two-box model holds one box, and reads back
/// as that box — while both roots still round-trip whole.
#[test]
fn a_snapshot_of_one_root_is_a_fraction_and_round_trips() {
    let mut model = Model::new();
    let b1 = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let far = Frame::new(Point::new(100.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let b2 = ogeom::algo::make_box(&mut model, far, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;

    let one = write(&model, std::slice::from_ref(&b1), WriteOptions::default()).unwrap();
    let all = write(&model, &[b1, b2], WriteOptions::default()).unwrap();

    // The subset subsets. Half the nodes, half the surfaces: comfortably
    // under two thirds of the whole, not within a rounding error of it.
    assert!(
        one.len() * 3 < all.len() * 2,
        "one box wrote {} bytes against {} for both",
        one.len(),
        all.len()
    );

    // And it is still that box.
    let (m1, roots) = read(&one).unwrap();
    assert_eq!(roots.len(), 1);
    assert!((volume(&m1, &roots[0]) - 1000.0).abs() < 1e-6);

    // The whole still round-trips whole, both roots alive.
    let (m2, roots) = read(&all).unwrap();
    assert_eq!(roots.len(), 2);
    assert!((volume(&m2, &roots[0]) - 1000.0).abs() < 1e-6);
    assert!((volume(&m2, &roots[1]) - 1000.0).abs() < 1e-6);
}

/// Subsetting is a property of the *roots*, not a new format: asking for
/// every root writes the model as it always did, handles unrenumbered.
#[test]
fn asking_for_every_root_writes_the_model_unchanged() {
    let mut model = Model::new();
    let b1 = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let whole = write(&model, std::slice::from_ref(&b1), WriteOptions::default()).unwrap();
    let (m, roots) = read(&whole).unwrap();
    assert_eq!(roots.len(), 1);
    assert!((volume(&m, &roots[0]) - 1000.0).abs() < 1e-6);
}
