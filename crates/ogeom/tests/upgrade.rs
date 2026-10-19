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

/// How many faces of `shape` lie on a surface of `kind`.
fn faces_on(model: &Model, shape: &ogeom::topo::Shape, kind: ogeom::geom::SurfaceKind) -> usize {
    use ogeom::geom::Surface as _;
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .iter()
        .filter(|f| {
            let Some(ogeom::topo::NodeData::Face(data)) = model.node(f).map(|n| n.data()) else {
                return false;
            };
            model
                .geometry()
                .surface(data.surface)
                .is_some_and(|s| s.kind() == kind)
        })
        .count()
}

/// A block carrying a round boss, a tapered boss and a dome, drilled
/// through, all centred on the plane x = 0, cut down to its half on one
/// side and fused with its own mirror image: the fuse leaves every curved
/// face the plane crosses in two, the image's half on a frame turned or
/// run the other way. Unified, each is one face again, the bore a whole
/// tube closed by a seam, the solid the whole block it was cut from.
#[test]
fn a_half_fused_with_its_mirror_image_unifies_its_curved_faces() {
    use ogeom::geom::SurfaceKind;
    let mut model = Model::new();
    let along_y = |x: f64, y: f64, z: f64| {
        Frame::new(Point::new(x, y, z), Direction::Y, Direction::Z, T).unwrap()
    };
    let corner = Frame::new(Point::new(-10.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let mut whole = ogeom::algo::make_box(&mut model, corner, (20.0, 30.0, 5.0), T)
        .unwrap()
        .shape;
    let boss = ogeom::algo::make_cylinder(&mut model, along_y(0.0, 2.0, 5.0), 3.0, 7.0, T)
        .unwrap()
        .shape;
    let taper = ogeom::algo::make_cone(&mut model, along_y(0.0, 12.0, 5.0), 3.0, 2.0, 6.0, T)
        .unwrap()
        .shape;
    let polar = Frame::new(Point::new(0.0, 24.0, 4.0), Direction::X, Direction::Y, T).unwrap();
    let dome = ogeom::algo::make_sphere(&mut model, polar, 3.0, T)
        .unwrap()
        .shape;
    for part in [boss, taper, dome] {
        whole = ogeom::boolean::fuse(&mut model, &whole, &part, T)
            .unwrap()
            .shape;
    }
    let drill = ogeom::algo::make_cylinder(&mut model, along_y(0.0, -1.0, 2.5), 1.0, 32.0, T)
        .unwrap()
        .shape;
    let whole = ogeom::boolean::cut(&mut model, &whole, &drill, T)
        .unwrap()
        .shape;
    let side = Frame::new(Point::new(0.0, -2.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let keep = ogeom::algo::make_box(&mut model, side, (12.0, 34.0, 12.0), T)
        .unwrap()
        .shape;
    let half = ogeom::boolean::common(&mut model, &whole, &keep, T)
        .unwrap()
        .shape;
    let copy = ogeom::algo::copied(&mut model, &half).unwrap().shape;
    let mirror = ogeom::math::Transform::plane_mirror(Point::ORIGIN, Direction::X);
    let image = ogeom::algo::transformed(&mut model, &copy, mirror)
        .unwrap()
        .shape;
    let fused = ogeom::boolean::fuse(&mut model, &half, &image, T)
        .unwrap()
        .shape;
    // The boss and the bore, the tapered boss, the dome.
    let kinds = [
        (SurfaceKind::Cylinder, 2),
        (SurfaceKind::Cone, 1),
        (SurfaceKind::Sphere, 1),
    ];
    for (kind, whole_count) in kinds {
        assert_eq!(faces_on(&model, &whole, kind), whole_count, "{kind:?}");
        assert_eq!(
            faces_on(&model, &fused, kind),
            2 * whole_count,
            "the fuse leaves each {kind:?} face split"
        );
    }

    let refined = ogeom::heal::unify_same_domain(&mut model, &fused, T)
        .unwrap()
        .0
        .shape;
    for (kind, whole_count) in kinds {
        assert_eq!(faces_on(&model, &refined, kind), whole_count, "{kind:?}");
    }
    let diagnosis = ogeom::algo::check(&model, &refined, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let mesh = ogeom::mesh::triangulate(&model, &refined, Deflection::default(), T).unwrap();
    assert!(
        mesh.is_closed(),
        "the refined solid draws closed, the bore's seam and all"
    );
    let (before, after) = (exact_volume(&model, &whole), exact_volume(&model, &refined));
    assert!(
        (after - before).abs() < before * 1e-9,
        "{after} against {before}"
    );
}

/// Unify `shape`, and check what unification must keep: the result is
/// valid, its mesh closes and its volume is `shape`'s. The volumes are
/// exact integrals over the faces, so the two agree to rounding; 1e-9
/// relative leaves room for the merged faces' boundaries being summed in
/// another order.
fn unified_keeping_the_solid(model: &mut Model, shape: &ogeom::topo::Shape) -> ogeom::topo::Shape {
    let before = exact_volume(model, shape);
    let refined = ogeom::heal::unify_same_domain(model, shape, T)
        .unwrap()
        .0
        .shape;
    let diagnosis = ogeom::algo::check(model, &refined, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let mesh = ogeom::mesh::triangulate(model, &refined, Deflection::default(), T).unwrap();
    assert!(mesh.is_closed(), "the unified solid draws closed");
    let after = exact_volume(model, &refined);
    assert!(
        (after - before).abs() < before * 1e-9,
        "{after} against {before}"
    );
    refined
}

/// Two coaxial cylinders fused end to end, on one frame: each side is a
/// whole band with its own seam, the two seams running up one line, and
/// the bands meet round the circle between them. Unified, they are one
/// band whose seam runs on from one band's into the other's.
#[test]
fn stacked_cylinder_bands_unify_into_one() {
    use ogeom::geom::SurfaceKind;
    let (r, h) = (3.0, 5.0);
    let mut model = Model::new();
    let lower = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, r, h, T)
        .unwrap()
        .shape;
    let second = Frame::new(Point::new(0.0, 0.0, h), Direction::Z, Direction::X, T).unwrap();
    let upper = ogeom::algo::make_cylinder(&mut model, second, r, h, T)
        .unwrap()
        .shape;
    let stacked = ogeom::boolean::fuse(&mut model, &lower, &upper, T)
        .unwrap()
        .shape;
    assert_eq!(faces_on(&model, &stacked, SurfaceKind::Cylinder), 2);

    let refined = unified_keeping_the_solid(&mut model, &stacked);
    assert_eq!(faces_on(&model, &refined, SurfaceKind::Cylinder), 1);
    assert_eq!(planar_faces(&model, &refined), 2);
    let drum = core::f64::consts::PI * r * r * 2.0 * h;
    let volume = exact_volume(&model, &refined);
    assert!(
        (volume - drum).abs() < drum * 1e-9,
        "{volume} against {drum}"
    );
}

/// Bands whose seams meet end to end but run opposite ways from the vertex
/// between them: two cylinders stacked, the second pointing down from the
/// top; and a sphere (axis along x) cut at its equator and fused with its
/// mirror image. One face would walk a seam line made of two edges, one
/// forward and one back, which offsets do not read as one seam. So each
/// pair stays split, and the solid still draws closed.
#[test]
fn seams_running_both_ways_from_a_vertex_stay_split() {
    use ogeom::geom::SurfaceKind;
    for kind in [SurfaceKind::Cylinder, SurfaceKind::Sphere] {
        let mut model = Model::new();
        let split = if kind == SurfaceKind::Cylinder {
            let lower = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 3.0, 5.0, T)
                .unwrap()
                .shape;
            let down = Frame::new(
                Point::new(0.0, 0.0, 10.0),
                Direction::Z.reversed(),
                Direction::X,
                T,
            )
            .unwrap();
            let upper = ogeom::algo::make_cylinder(&mut model, down, 3.0, 5.0, T)
                .unwrap()
                .shape;
            ogeom::boolean::fuse(&mut model, &lower, &upper, T)
                .unwrap()
                .shape
        } else {
            let polar = Frame::new(Point::ORIGIN, Direction::X, Direction::Y, T).unwrap();
            let ball = ogeom::algo::make_sphere(&mut model, polar, 5.0, T)
                .unwrap()
                .shape;
            let corner =
                Frame::new(Point::new(0.0, -6.0, -6.0), Direction::Z, Direction::X, T).unwrap();
            let keep = ogeom::algo::make_box(&mut model, corner, (6.0, 12.0, 12.0), T)
                .unwrap()
                .shape;
            let half = ogeom::boolean::common(&mut model, &ball, &keep, T)
                .unwrap()
                .shape;
            let copy = ogeom::algo::copied(&mut model, &half).unwrap().shape;
            let mirror = ogeom::math::Transform::plane_mirror(Point::ORIGIN, Direction::X);
            let image = ogeom::algo::transformed(&mut model, &copy, mirror)
                .unwrap()
                .shape;
            ogeom::boolean::fuse(&mut model, &half, &image, T)
                .unwrap()
                .shape
        };
        assert_eq!(faces_on(&model, &split, kind), 2, "{kind:?}");

        let refined = unified_keeping_the_solid(&mut model, &split);
        assert_eq!(faces_on(&model, &refined, kind), 2, "{kind:?}");
    }
}

/// Bands whose seams start at different angles round the circle they
/// share: two cylinders stacked, the second turned a quarter round; the
/// upper half of a torus fused with its mirror image turned 0.7 rad; and a
/// sphere's half fused with its mirror image turned 0.7 rad about the
/// axis. One face would need a seam running along part of the shared
/// circle between the two seams, which meshes, booleans and offsets do not
/// read as a seam. So each pair stays split, and the solid still draws
/// closed and takes a box cut through it.
#[test]
fn bands_whose_seams_start_apart_stay_split() {
    use ogeom::geom::SurfaceKind;
    use ogeom::math::{Axis, Transform};
    let corner = |x: f64, y: f64, z: f64| {
        Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
    };
    let mirrored_and_turned = |model: &mut Model,
                               solid: &ogeom::topo::Shape,
                               keep: (Frame, (f64, f64, f64)),
                               axis: Direction| {
        let keep = ogeom::algo::make_box(model, keep.0, keep.1, T)
            .unwrap()
            .shape;
        let half = ogeom::boolean::common(model, solid, &keep, T)
            .unwrap()
            .shape;
        let copy = ogeom::algo::copied(model, &half).unwrap().shape;
        let turn = Transform::rotation(Axis::new(Point::ORIGIN, axis), 0.7)
            * Transform::plane_mirror(Point::ORIGIN, axis);
        let image = ogeom::algo::transformed(model, &copy, turn).unwrap().shape;
        ogeom::boolean::fuse(model, &half, &image, T).unwrap().shape
    };
    for kind in [
        SurfaceKind::Cylinder,
        SurfaceKind::Torus,
        SurfaceKind::Sphere,
    ] {
        let mut model = Model::new();
        let (split, cutter) = match kind {
            SurfaceKind::Cylinder => {
                let lower = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 3.0, 5.0, T)
                    .unwrap()
                    .shape;
                let second =
                    Frame::new(Point::new(0.0, 0.0, 5.0), Direction::Z, Direction::Y, T).unwrap();
                let upper = ogeom::algo::make_cylinder(&mut model, second, 3.0, 5.0, T)
                    .unwrap()
                    .shape;
                let stacked = ogeom::boolean::fuse(&mut model, &lower, &upper, T)
                    .unwrap()
                    .shape;
                (stacked, (corner(1.0, -4.0, 2.0), (5.0, 8.0, 6.0)))
            }
            SurfaceKind::Torus => {
                let torus = ogeom::algo::make_torus(&mut model, Frame::WORLD, 10.0, 3.0, T)
                    .unwrap()
                    .shape;
                let keep = (corner(-15.0, -15.0, 0.0), (30.0, 30.0, 5.0));
                let fused = mirrored_and_turned(&mut model, &torus, keep, Direction::Z);
                (fused, (corner(2.0, -20.0, -5.0), (20.0, 40.0, 10.0)))
            }
            _ => {
                let polar = Frame::new(Point::ORIGIN, Direction::X, Direction::Y, T).unwrap();
                let ball = ogeom::algo::make_sphere(&mut model, polar, 5.0, T)
                    .unwrap()
                    .shape;
                let keep = (corner(0.0, -6.0, -6.0), (6.0, 12.0, 12.0));
                let fused = mirrored_and_turned(&mut model, &ball, keep, Direction::X);
                (fused, (corner(-1.0, -1.0, -6.0), (8.0, 8.0, 12.0)))
            }
        };
        assert_eq!(faces_on(&model, &split, kind), 2, "{kind:?}");

        let refined = unified_keeping_the_solid(&mut model, &split);
        assert_eq!(faces_on(&model, &refined, kind), 2, "{kind:?}");
        let block = ogeom::algo::make_box(&mut model, cutter.0, cutter.1, T)
            .unwrap()
            .shape;
        let cut = ogeom::boolean::cut(&mut model, &refined, &block, T)
            .unwrap()
            .shape;
        let mesh = ogeom::mesh::triangulate(&model, &cut, Deflection::default(), T).unwrap();
        assert!(mesh.is_closed(), "{kind:?}: the cut draws closed");
    }
}

/// A prism on an ellipse centred on the plane x = 0, cut to its half on
/// one side and fused with its mirror image. The image's wall is swept
/// from the mirrored ellipse, whose eccentric angle runs the other way.
/// Unified, the wall is one face again, a whole tube closed by a seam.
#[test]
fn an_elliptic_wall_fused_with_its_mirror_image_unifies() {
    use ogeom::geom::SurfaceKind;
    let (major, minor, h) = (4.0, 2.0, 5.0);
    let mut model = Model::new();
    // The ellipse's `x` runs along y and its `y` along -x, so the mirror
    // in x = 0 runs its angle backwards.
    let frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::Y, T).unwrap();
    let ellipse = ogeom::math::Ellipse::new(frame, major, minor, T).unwrap();
    let rim = ogeom::algo::make_edge(
        &mut model,
        ogeom::geom::Curve::Ellipse(ogeom::geom::EllipseCurve::new(ellipse)),
        (0.0, core::f64::consts::TAU),
        T,
    )
    .unwrap()
    .shape;
    let wire = ogeom::algo::make_wire(&mut model, &[rim], T).unwrap().shape;
    let plane = ogeom::geom::PlaneSurface::new(ogeom::math::Plane::new(Frame::WORLD));
    let base = ogeom::algo::make_face(&mut model, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let whole =
        ogeom::algo::make_prism(&mut model, &base, ogeom::math::Vector::new(0.0, 0.0, h), T)
            .unwrap()
            .shape;
    assert_eq!(faces_on(&model, &whole, SurfaceKind::Extrusion), 1);
    let side = Frame::new(Point::new(0.0, -6.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let keep = ogeom::algo::make_box(&mut model, side, (6.0, 12.0, h + 2.0), T)
        .unwrap()
        .shape;
    let half = ogeom::boolean::common(&mut model, &whole, &keep, T)
        .unwrap()
        .shape;
    let copy = ogeom::algo::copied(&mut model, &half).unwrap().shape;
    let mirror = ogeom::math::Transform::plane_mirror(Point::ORIGIN, Direction::X);
    let image = ogeom::algo::transformed(&mut model, &copy, mirror)
        .unwrap()
        .shape;
    let fused = ogeom::boolean::fuse(&mut model, &half, &image, T)
        .unwrap()
        .shape;
    assert_eq!(
        faces_on(&model, &fused, SurfaceKind::Extrusion),
        2,
        "the fuse leaves the wall split"
    );

    let refined = unified_keeping_the_solid(&mut model, &fused);
    assert_eq!(faces_on(&model, &refined, SurfaceKind::Extrusion), 1);
    assert_eq!(planar_faces(&model, &refined), 2);
    let prism = core::f64::consts::PI * major * minor * h;
    let volume = exact_volume(&model, &refined);
    assert!(
        (volume - prism).abs() < prism * 1e-9,
        "{volume} against {prism}"
    );
}

/// A half torus swept from a disk, each half covering the whole tube and
/// so closed on itself by the tube's seam, fused with a copy turned half
/// round the axis: the two meet along both their end circles. Unified,
/// they are the whole torus: one face, closed both ways round.
#[test]
fn torus_halves_round_the_whole_tube_unify() {
    use ogeom::geom::SurfaceKind;
    let (major, minor) = (10.0, 3.0);
    let mut model = Model::new();
    let at = Frame::new(Point::new(major, 0.0, 0.0), Direction::Y, Direction::X, T).unwrap();
    let rim = ogeom::algo::make_edge(
        &mut model,
        ogeom::geom::Curve::Circle(ogeom::geom::CircleCurve::new(
            ogeom::math::Circle::new(at, minor, T).unwrap(),
        )),
        (0.0, core::f64::consts::TAU),
        T,
    )
    .unwrap()
    .shape;
    let wire = ogeom::algo::make_wire(&mut model, &[rim], T).unwrap().shape;
    let disk = ogeom::geom::PlaneSurface::new(ogeom::math::Plane::new(at));
    let disk = ogeom::algo::make_face(&mut model, disk.into(), &[wire], T)
        .unwrap()
        .shape;
    let axis = ogeom::math::Axis::new(Point::ORIGIN, Direction::Z);
    let half = ogeom::algo::make_revolution(&mut model, &disk, axis, core::f64::consts::PI, T)
        .unwrap()
        .shape;
    let copy = ogeom::algo::copied(&mut model, &half).unwrap().shape;
    let turn = ogeom::math::Transform::rotation(axis, core::f64::consts::PI);
    let other = ogeom::algo::transformed(&mut model, &copy, turn)
        .unwrap()
        .shape;
    let whole = ogeom::boolean::fuse(&mut model, &half, &other, T)
        .unwrap()
        .shape;
    assert_eq!(faces_on(&model, &whole, SurfaceKind::Torus), 2);

    let refined = unified_keeping_the_solid(&mut model, &whole);
    assert_eq!(
        explore_unique(&model, &refined, ShapeType::Face)
            .unwrap()
            .len(),
        1
    );
    // Pappus: the disk's area swept round the circle its centre runs.
    let torus = core::f64::consts::PI * minor * minor * core::f64::consts::TAU * major;
    let volume = exact_volume(&model, &refined);
    assert!(
        (volume - torus).abs() < torus * 1e-9,
        "{volume} against {torus}"
    );
}

/// `shape` with one of its curved faces cut along each line in turn, each
/// time the first face the line crosses.
fn cut_along(
    model: &mut Model,
    shape: &ogeom::topo::Shape,
    lines: &[ogeom::heal::IsoLine],
) -> ogeom::topo::Shape {
    let mut shape = shape.clone();
    for &line in lines {
        let face = explore_unique(model, &shape, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                planar_faces(model, f) == 0
                    && ogeom::heal::divide_face(&mut model.clone(), &shape, f, line, T).is_ok()
            })
            .unwrap();
        shape = ogeom::heal::divide_face(model, &shape, &face, line, T)
            .unwrap()
            .shape;
    }
    shape
}

/// A drum's side and a torus cut into pieces along their own parameters,
/// so every piece lies on the one surface: cut in two round the axis, the
/// halves share two edges, one inside the old face and the old seam,
/// which each half reads on its own side; cut again, three or four pieces
/// meet round the axis and across it. Unified, the pieces are one face
/// again, the old seam its seam. So are the halves where the old seam
/// carries one pcurve for both, as a file can store it, and one half reads
/// it a period away from the rest of its boundary.
#[test]
fn pieces_of_one_surface_unify_round_the_axis() {
    use core::f64::consts::PI;
    use ogeom::geom::SurfaceKind;
    use ogeom::heal::IsoLine;
    use ogeom::topo::{EdgeRepr, NodeData};
    let (r, h) = (3.0, 5.0);
    let (major, minor) = (10.0, 3.0);
    let drum = PI * r * r * h;
    // Pappus: the tube's disk swept round the circle its centre runs.
    let ring = PI * minor * minor * 2.0 * PI * major;
    let cases: [(SurfaceKind, &[IsoLine], bool); 5] = [
        (SurfaceKind::Cylinder, &[IsoLine::U(PI)], false),
        (SurfaceKind::Cylinder, &[IsoLine::U(PI)], true),
        (
            SurfaceKind::Cylinder,
            &[IsoLine::U(PI), IsoLine::U(1.0)],
            false,
        ),
        (
            SurfaceKind::Cylinder,
            &[IsoLine::U(PI), IsoLine::U(1.0), IsoLine::V(2.0)],
            false,
        ),
        (SurfaceKind::Torus, &[IsoLine::U(PI), IsoLine::V(PI)], false),
    ];
    for (kind, lines, one_pcurve) in cases {
        let mut model = Model::new();
        let (solid, closed_form) = if kind == SurfaceKind::Torus {
            let torus = ogeom::algo::make_torus(&mut model, Frame::WORLD, major, minor, T)
                .unwrap()
                .shape;
            (torus, ring)
        } else {
            let cylinder = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, r, h, T)
                .unwrap()
                .shape;
            (cylinder, drum)
        };
        let pieces = cut_along(&mut model, &solid, lines);
        assert_eq!(faces_on(&model, &pieces, kind), lines.len() + 1);
        if one_pcurve {
            for edge in explore_unique(&model, &pieces, ShapeType::Edge).unwrap() {
                let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(|n| n.data_mut()) else {
                    continue;
                };
                for repr in &mut data.representations {
                    if let EdgeRepr::Seam {
                        forward,
                        surface,
                        location,
                        range,
                        ..
                    } = repr.clone()
                    {
                        *repr = EdgeRepr::PCurve {
                            curve: forward,
                            surface,
                            location,
                            range,
                        };
                    }
                }
            }
            let diagnosis = ogeom::algo::check(&model, &pieces, T).unwrap();
            assert!(diagnosis.is_valid(), "{diagnosis}");
        }

        let refined = unified_keeping_the_solid(&mut model, &pieces);
        assert_eq!(
            faces_on(&model, &refined, kind),
            1,
            "{kind:?} cut along {lines:?}"
        );
        let volume = exact_volume(&model, &refined);
        assert!(
            (volume - closed_form).abs() < closed_form * 1e-9,
            "{volume} against {closed_form}"
        );
    }
}

/// A cylinder cut to its half on one side of x = 0 and fused with its
/// mirror image, the image's pcurves on the two lines the halves share
/// restated as linear splines over twice their range: the same chart
/// lines, read at another pace. Unified, the halves are one tube again,
/// the line whose sides land a period apart its seam, both its sides the
/// one curve a period across.
#[test]
fn half_tubes_whose_shared_lines_read_at_another_pace_unify() {
    use ogeom::geom::{Curve2d as _, SurfaceGeometry, SurfaceKind};
    use ogeom::topo::{EdgeRepr, NodeData};
    let (r, h) = (3.0, 5.0);
    let mut model = Model::new();
    let whole = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, r, h, T)
        .unwrap()
        .shape;
    let side = Frame::new(Point::new(0.0, -5.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let keep = ogeom::algo::make_box(&mut model, side, (5.0, 10.0, h + 2.0), T)
        .unwrap()
        .shape;
    let half = ogeom::boolean::common(&mut model, &whole, &keep, T)
        .unwrap()
        .shape;
    let copy = ogeom::algo::copied(&mut model, &half).unwrap().shape;
    let mirror = ogeom::math::Transform::plane_mirror(Point::ORIGIN, Direction::X);
    let image = ogeom::algo::transformed(&mut model, &copy, mirror)
        .unwrap()
        .shape;
    let fused = ogeom::boolean::fuse(&mut model, &half, &image, T)
        .unwrap()
        .shape;
    let tubes: Vec<_> = explore_unique(&model, &fused, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| faces_on(&model, f, SurfaceKind::Cylinder) == 1)
        .collect();
    assert_eq!(tubes.len(), 2, "the fuse leaves the tube split");

    let first: Vec<_> = explore_unique(&model, &tubes[0], ShapeType::Edge)
        .unwrap()
        .iter()
        .map(|e| e.node())
        .collect();
    let Some(NodeData::Face(data)) = model.node(&tubes[1]).map(|n| n.data().clone()) else {
        panic!("a face holds face data");
    };
    assert!(matches!(
        model.geometry().surface(data.surface),
        Some(SurfaceGeometry::Cylinder(_))
    ));
    let mut restated = 0;
    for edge in explore_unique(&model, &tubes[1], ShapeType::Edge).unwrap() {
        if !first.contains(&edge.node()) {
            continue;
        }
        let edge_data = model.node(&edge).unwrap().data().as_edge().unwrap();
        let Some(EdgeRepr::PCurve { curve, range, .. }) =
            edge_data.pcurve_for(data.surface, edge.location()).cloned()
        else {
            panic!("a shared line carries a pcurve on each half");
        };
        let line = model.geometry().pcurve(curve).unwrap().clone();
        let ends = vec![
            line.point_at(range.0, T).unwrap(),
            line.point_at(range.1, T).unwrap(),
        ];
        let span = 2.0 * (range.1 - range.0);
        let knots = ogeom::math::KnotVector::new(vec![0.0, 0.0, span, span], 1).unwrap();
        let paced = ogeom::geom::BSpline2d::new(knots, ends, T).unwrap();
        if let Some(NodeData::Edge(held)) = model.node_mut(&edge).map(|n| n.data_mut()) {
            held.representations
                .retain(|r| r.surface() != Some(data.surface));
        }
        ogeom::algo::attach_pcurve(
            &mut model,
            &edge,
            paced.into(),
            data.surface,
            edge.location().clone(),
            (0.0, span),
        )
        .unwrap();
        restated += 1;
    }
    assert_eq!(restated, 2, "the halves share two lines");
    let diagnosis = ogeom::algo::check(&model, &fused, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");

    let refined = unified_keeping_the_solid(&mut model, &fused);
    assert_eq!(faces_on(&model, &refined, SurfaceKind::Cylinder), 1);
    let drum = core::f64::consts::PI * r * r * h;
    let volume = exact_volume(&model, &refined);
    assert!(
        (volume - drum).abs() < drum * 1e-9,
        "{volume} against {drum}"
    );
}

/// A whole drum with a tube stacked on it, the tube made of two halves
/// fused across x = 0, and the drum's seam damaged as a file can carry it:
/// both its sides read on one sheet, so the drum closes into no band and
/// the three walls make no one face. The drum is left as it is, and the
/// tube's halves still merge into one tube, as with no drum beneath them.
#[test]
fn a_band_whose_seam_closes_nothing_leaves_its_neighbours_merging() {
    use ogeom::geom::SurfaceKind;
    use ogeom::topo::{EdgeRepr, NodeData};
    let (r, h) = (3.0, 5.0);
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, r, h, T)
        .unwrap()
        .shape;
    let up = Frame::new(Point::new(0.0, 0.0, h), Direction::Z, Direction::X, T).unwrap();
    let tube = ogeom::algo::make_cylinder(&mut model, up, r, h, T)
        .unwrap()
        .shape;
    let side = Frame::new(
        Point::new(0.0, -5.0, h - 1.0),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    let keep = ogeom::algo::make_box(&mut model, side, (5.0, 10.0, h + 2.0), T)
        .unwrap()
        .shape;
    let half = ogeom::boolean::common(&mut model, &tube, &keep, T)
        .unwrap()
        .shape;
    let copy = ogeom::algo::copied(&mut model, &half).unwrap().shape;
    let mirror = ogeom::math::Transform::plane_mirror(Point::ORIGIN, Direction::X);
    let image = ogeom::algo::transformed(&mut model, &copy, mirror)
        .unwrap()
        .shape;
    let halves = ogeom::boolean::fuse(&mut model, &half, &image, T)
        .unwrap()
        .shape;
    let stacked = ogeom::boolean::fuse(&mut model, &drum, &halves, T)
        .unwrap()
        .shape;
    assert_eq!(faces_on(&model, &stacked, SurfaceKind::Cylinder), 3);
    let mut damaged = 0;
    for edge in explore_unique(&model, &stacked, ShapeType::Edge).unwrap() {
        let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(|n| n.data_mut()) else {
            continue;
        };
        for repr in &mut data.representations {
            if let EdgeRepr::Seam {
                forward, reversed, ..
            } = repr
            {
                *reversed = *forward;
                damaged += 1;
            }
        }
    }
    assert_eq!(damaged, 1, "only the drum walks a seam");

    let refined = ogeom::heal::unify_same_domain(&mut model, &stacked, T)
        .unwrap()
        .0
        .shape;
    assert_eq!(faces_on(&model, &refined, SurfaceKind::Cylinder), 2);
}
