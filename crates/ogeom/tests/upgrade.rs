//! The upgrade family, end to end: a fused pair of half-boxes carries the
//! splits its construction left, and the upgrades take them back out
//! without changing what the solid is.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_fused_pair_of_halves_unifies_back_to_a_box() {
    let mut model = Model::new();
    let left = ogeom::algo::make_box(&mut model, Frame::WORLD, (5.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let f = Frame::new(Point::new(5.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let right = ogeom::algo::make_box(&mut model, f, (5.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let fused = ogeom::boolean::fuse(&mut model, &left, &right, T)
        .unwrap()
        .shape;
    assert_eq!(
        explore_unique(&model, &fused, ShapeType::Face)
            .unwrap()
            .len(),
        10,
        "the fuse leaves four walls split"
    );

    let unified = ogeom::heal::unify_same_domain(&mut model, &fused, T)
        .unwrap()
        .0
        .shape;
    assert_eq!(
        explore_unique(&model, &unified, ShapeType::Face)
            .unwrap()
            .len(),
        6,
        "one face per box side again"
    );
    let volume = ogeom::algo::volume_properties(&model, &unified, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!(
        (volume - 1000.0).abs() < 1e-6,
        "nothing changed shape: {volume}"
    );

    // The unified walls still carry split collinear boundary edges. The
    // edge merge takes those out too.
    let merged = ogeom::heal::merge_edges(&mut model, &unified, T)
        .unwrap()
        .0
        .shape;
    let edges = explore_unique(&model, &merged, ShapeType::Edge)
        .unwrap()
        .len();
    assert!(
        edges <= 12,
        "a box has twelve edges; the merge got to {edges}"
    );
    let volume = ogeom::algo::volume_properties(&model, &merged, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!(
        (volume - 1000.0).abs() < 1e-6,
        "still the same box: {volume}"
    );
}

#[test]
fn tolerances_shrink_back_to_what_is_measured() {
    let mut model = Model::new();
    let solid = ogeom::algo::make_box(&mut model, Frame::WORLD, (4.0, 5.0, 6.0), T)
        .unwrap()
        .shape;
    // Widen one edge artificially, as an import might have.
    let edge = explore_unique(&model, &solid, ShapeType::Edge).unwrap()[0].clone();
    if let Some(node) = model.node_mut(&edge)
        && let ogeom::topo::NodeData::Edge(data) = node.data_mut()
    {
        data.tolerance = data.tolerance.widen_to(0.5);
    }
    let shrunk = ogeom::heal::reduce_tolerances(&mut model, &solid, T).unwrap();
    assert!(shrunk >= 1, "the widened claim shrank");
    if let Some(node) = model.node(&edge) {
        let claimed = node.data().tolerance().unwrap().get();
        assert!(claimed < 1e-3, "back to the measured agreement: {claimed}");
    }
}

/// Four boxes in a row, fused: the top, the bottom and both long sides
/// are each split in four. Unified around one piece of the top, the top
/// becomes one face and every other split is left alone; unified whole,
/// all of them merge.
#[test]
fn unifying_around_a_face_merges_its_carrier_and_nothing_else() {
    let face_count = |model: &Model, s: &ogeom::topo::Shape| {
        explore_unique(model, s, ShapeType::Face).unwrap().len()
    };
    let volume = |model: &Model, s: &ogeom::topo::Shape| {
        ogeom::algo::volume_properties(model, s, Deflection::default(), T)
            .unwrap()
            .mass
    };
    let mut model = Model::new();
    let mut row = ogeom::algo::make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
        .unwrap()
        .shape;
    for k in 1..4 {
        let at = Frame::new(
            Point::new(f64::from(k), 0.0, 0.0),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        let next = ogeom::algo::make_box(&mut model, at, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        row = ogeom::boolean::fuse(&mut model, &row, &next, T)
            .unwrap()
            .shape;
    }
    let faces = explore_unique(&model, &row, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 18);
    let top = faces
        .iter()
        .find(|f| {
            let (p, n) = ogeom::algo::face_normal(&model, f, T).unwrap();
            n.z > 0.999 && p.x < 1.0
        })
        .unwrap()
        .clone();

    let (around, _) =
        ogeom::heal::unify_same_domain_around(&mut model, &row, std::slice::from_ref(&top), T)
            .unwrap();
    assert_eq!(face_count(&model, &around.shape), 15);
    assert!(
        ogeom::algo::check(&model, &around.shape, T)
            .unwrap()
            .is_valid()
    );
    assert!((volume(&model, &around.shape) - 4.0).abs() < 1e-9);

    let (whole, _) = ogeom::heal::unify_same_domain(&mut model, &row, T).unwrap();
    assert_eq!(face_count(&model, &whole.shape), 6);
}
