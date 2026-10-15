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

/// Equal inputs in fresh models unify to equal results: the merged faces
/// enter the arena in one order, so node ids and face order repeat from
/// run to run.
#[test]
fn unify_same_domain_is_deterministic() {
    let run = || {
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
        let (whole, _) = ogeom::heal::unify_same_domain(&mut model, &row, T).unwrap();
        explore_unique(&model, &whole.shape, ShapeType::Face)
            .unwrap()
            .iter()
            .map(|f| f.node().index())
            .collect::<Vec<_>>()
    };
    let first = run();
    assert_eq!(first.len(), 6);
    for _ in 0..4 {
        assert_eq!(run(), first, "the merged faces keep their node ids");
    }
}

fn planar_faces(model: &Model, shape: &ogeom::topo::Shape) -> usize {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .iter()
        .filter(|f| {
            let Some(ogeom::topo::NodeData::Face(data)) = model.node(f).map(|n| n.data()) else {
                return false;
            };
            matches!(
                model.geometry().surface(data.surface),
                Some(ogeom::geom::SurfaceGeometry::Plane(_))
            )
        })
        .count()
}

fn exact_volume(model: &Model, shape: &ogeom::topo::Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// A flat face of `outer` with `inner` cut out of it, at height `z`, swept
/// up by `h`. The rings are given as walked; the face winds them itself.
fn pad(
    model: &mut Model,
    outer: &[(f64, f64)],
    inner: &[(f64, f64)],
    z: f64,
    h: f64,
) -> ogeom::topo::Shape {
    let mut wires = Vec::new();
    for ring in [outer, inner] {
        if ring.is_empty() {
            continue;
        }
        let points: Vec<Point> = ring.iter().map(|&(x, y)| Point::new(x, y, z)).collect();
        wires.push(
            ogeom::algo::make_polygon(model, &points, true, T)
                .unwrap()
                .shape,
        );
    }
    let at = Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X, T).unwrap();
    let plane = ogeom::geom::PlaneSurface::new(ogeom::math::Plane::new(at));
    let face = ogeom::algo::make_face(model, plane.into(), &wires, T)
        .unwrap()
        .shape;
    ogeom::algo::make_prism(model, &face, ogeom::math::Vector::new(0.0, 0.0, h), T)
        .unwrap()
        .shape
}

/// A plate with a rim padded round its edge, the rim's profile walked the
/// other way round from the plate's: the rim's outer walls store their
/// planes facing into the part and are used reversed, the plate's sides
/// face out. Each pair is one wall, facing one way, and merges: the open
/// box has eleven faces, as it was drawn.
#[test]
fn walls_storing_their_planes_opposite_ways_still_merge() {
    let mut model = Model::new();
    let plate = pad(
        &mut model,
        &[(0.0, 0.0), (60.0, 0.0), (60.0, 180.0), (0.0, 180.0)],
        &[],
        0.0,
        2.5,
    );
    let rim = pad(
        &mut model,
        &[(60.0, 0.0), (0.0, 0.0), (0.0, 180.0), (60.0, 180.0)],
        &[(2.5, 2.5), (57.5, 2.5), (57.5, 177.5), (2.5, 177.5)],
        2.5,
        28.5,
    );
    let open_box = ogeom::boolean::fuse(&mut model, &plate, &rim, T)
        .unwrap()
        .shape;
    let before = exact_volume(&model, &open_box);

    let refined = ogeom::heal::unify_same_domain(&mut model, &open_box, T)
        .unwrap()
        .0
        .shape;
    assert_eq!(
        explore_unique(&model, &refined, ShapeType::Face)
            .unwrap()
            .len(),
        11
    );
    let diagnosis = ogeom::algo::check(&model, &refined, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let after = exact_volume(&model, &refined);
    assert!(
        (after - before).abs() < before * 1e-12,
        "{after} against {before}"
    );
}

/// Two plates side by side, fused, then drilled through on the line they
/// meet along: the top halves each hold half the hole, so the merged top
/// is bounded by two rings, the plate's outline and the hole's rim. So is
/// the bottom. Each merges into one face with a hole.
#[test]
fn a_merged_face_keeps_the_hole_its_pieces_shared() {
    let mut model = Model::new();
    let left = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 20.0, 5.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(10.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let right = ogeom::algo::make_box(&mut model, at, (10.0, 20.0, 5.0), T)
        .unwrap()
        .shape;
    let plate = ogeom::boolean::fuse(&mut model, &left, &right, T)
        .unwrap()
        .shape;
    let axis = Frame::new(Point::new(10.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, axis, 3.0, 7.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom::boolean::cut(&mut model, &plate, &drill, T)
        .unwrap()
        .shape;
    assert_eq!(
        planar_faces(&model, &drilled),
        10,
        "the fuse left the plate split"
    );
    let before = exact_volume(&model, &drilled);

    let refined = ogeom::heal::unify_same_domain(&mut model, &drilled, T)
        .unwrap()
        .0
        .shape;
    assert_eq!(planar_faces(&model, &refined), 6);
    let diagnosis = ogeom::algo::check(&model, &refined, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let after = exact_volume(&model, &refined);
    assert!(
        (after - before).abs() < before * 1e-12,
        "{after} against {before}"
    );
}
