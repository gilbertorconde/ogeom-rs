//! A model compacted in place keeps what its roots reach under the handles
//! it had, drops the rest, and goes on working as the uncompacted model
//! does.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::FastSet;

use ogeom::algo::{make_box, make_cylinder};
use ogeom::core::Tolerances;
use ogeom::math::{Aabb, Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Filter, Model, Shape, ShapeType, SurfaceId, explore, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A pin of radius 1.5 standing 20 tall from `at`.
fn pin(model: &mut Model, at: Point) -> Shape {
    let frame = Frame::new(at, Direction::Z, Direction::X, T).unwrap();
    make_cylinder(model, frame, 1.5, 20.0, T).unwrap().shape
}

/// A plate 10 thick with a `rows` by `rows` grid of holes of radius 1.5,
/// 10 apart and 10 from its sides, each hole its own cut, all in one
/// model. Also returns the first pin, which no later result holds.
fn holed_plate(rows: u32) -> (Model, Shape, Shape) {
    let mut model = Model::new();
    let size = 10.0 * f64::from(rows) + 10.0;
    let mut plate = make_box(&mut model, Frame::WORLD, (size, size, 10.0), T)
        .unwrap()
        .shape;
    let mut first = None;
    for i in 0..rows {
        for j in 0..rows {
            let at = Point::new(10.0 + 10.0 * f64::from(i), 10.0 + 10.0 * f64::from(j), -5.0);
            let tool = pin(&mut model, at);
            first.get_or_insert_with(|| tool.clone());
            plate = ogeom::boolean::cut(&mut model, &plate, &tool, T)
                .unwrap()
                .shape;
        }
    }
    (model, plate, first.unwrap())
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// Each face of `shape`, with its surface handle and its bounds.
fn faces(model: &Model, shape: &Shape) -> Vec<(Shape, SurfaceId, Aabb)> {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .map(|face| {
            let surface = model.node(&face).unwrap().data().as_face().unwrap().surface;
            let bounds = ogeom::algo::face_bounds(model, &face).unwrap();
            (face, surface, bounds)
        })
        .collect()
}

/// The distinct nodes below `shape`, itself included.
fn reached_nodes(model: &Model, shape: &Shape) -> usize {
    explore(model, shape, Filter::All)
        .unwrap()
        .iter()
        .map(Shape::node)
        .collect::<FastSet<_>>()
        .len()
}

/// The plate after one more hole, at its corner: (volume, faces, valid,
/// nodes an editor of the result copies because the plate holds them too).
fn corner_hole(model: &mut Model, plate: &Shape) -> (f64, usize, bool, usize) {
    let corner = pin(model, Point::new(5.0, 5.0, -5.0));
    let cut = ogeom::boolean::cut(model, plate, &corner, T).unwrap().shape;
    let faces = explore_unique(model, &cut, ShapeType::Face).unwrap().len();
    let valid = ogeom::algo::check(model, &cut, T).unwrap().is_valid();
    let measured = volume(model, &cut);
    let copied = model.unshare(&cut).unwrap().len();
    (measured, faces, valid, copied)
}

/// Every face handle of the plate resolves after the compaction, to the
/// same surface and the same bounds; the nodes left are the plate's own;
/// a dropped pin and its geometry and lineage no longer resolve; and a
/// cut on the compacted plate builds the solid the same cut builds on the
/// uncompacted model.
#[test]
fn a_retained_model_keeps_its_handles_and_drops_the_rest() {
    let (mut model, plate, first) = holed_plate(3);
    let mut whole = model.clone();
    let before = faces(&model, &plate);
    let plate_volume = volume(&model, &plate);
    let lineage: Vec<_> = before
        .iter()
        .map(|(face, _, _)| model.roots_of(face))
        .collect();
    // The first pin's end caps lie clear of the plate, so nothing the plate
    // is made of shares their surfaces or derives from them.
    let cap = explore_unique(&model, &first, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|face| {
            ogeom::algo::face_bounds(&model, face)
                .unwrap()
                .high()
                .is_some_and(|high| high.z < 0.0)
        })
        .expect("the pin has a cap below the plate");
    let cap_surface = model.node(&cap).unwrap().data().as_face().unwrap().surface;
    let cap_entity = model.identity_of(&cap).unwrap();
    let (nodes, entries) = (model.node_count(), model.provenance().iter().count());
    let (curves, pcurves, surfaces) = model.geometry().counts();

    model
        .retain_reachable(std::slice::from_ref(&plate))
        .unwrap();

    assert_eq!(model.node_count(), reached_nodes(&model, &plate));
    assert!(
        model.node_count() < nodes,
        "{} of {nodes}",
        model.node_count()
    );
    let (c, p, s) = model.geometry().counts();
    assert!(c < curves && p < pcurves && s < surfaces);
    assert!(model.provenance().iter().count() < entries);

    assert_eq!(faces(&model, &plate), before);
    for ((face, _, _), roots) in before.iter().zip(&lineage) {
        assert_eq!(&model.roots_of(face), roots);
    }
    assert!(model.node(&first).is_none());
    assert!(model.kind_of(&cap).is_err());
    assert!(model.geometry().surface(cap_surface).is_none());
    assert!(model.provenance().get(cap_entity).is_none());
    assert_eq!(model.shape_of(cap_entity), None);
    // Measured on the same faces under the same deflection, so the same
    // sum; the bound is rounding.
    assert!((volume(&model, &plate) - plate_volume).abs() <= plate_volume * 1e-12);

    let after = corner_hole(&mut model, &plate);
    let want = corner_hole(&mut whole, &plate);
    assert!(after.2 && want.2, "both cuts are valid solids");
    assert_eq!(after.1, want.1);
    assert!(want.3 > 0, "the cut passes faces of the plate through");
    assert_eq!(after.3, want.3, "the result holds what it did");
    // The same cut on the same faces; the bound is rounding.
    assert!(
        (after.0 - want.0).abs() <= want.0 * 1e-9,
        "{after:?} {want:?}"
    );
    // A handle made after the compaction does not land on a dropped one.
    assert!(model.node(&first).is_none());
}

/// A compacted model writes out and reads back, whole and by root, and
/// takes a document in, as a model that never dropped anything does.
#[test]
fn a_retained_model_writes_reads_and_absorbs() {
    let (mut model, plate, _) = holed_plate(2);
    model
        .retain_reachable(std::slice::from_ref(&plate))
        .unwrap();
    let want = volume(&model, &plate);
    let options = ogeom::io::native::WriteOptions::default();

    let (small, roots) =
        ogeom::io::native::compacted(&model, std::slice::from_ref(&plate)).unwrap();
    assert!(ogeom::algo::check(&small, &roots[0], T).unwrap().is_valid());
    assert!((volume(&small, &roots[0]) - want).abs() <= want * 1e-12);

    let text = ogeom::io::native::write(&model, &[], options).unwrap();
    let (read, _) = ogeom::io::native::read(&text).unwrap();
    assert_eq!(read.node_count(), model.node_count());

    let mut other = Model::new();
    let block = make_box(&mut other, Frame::WORLD, (1.0, 2.0, 3.0), T)
        .unwrap()
        .shape;
    let text = ogeom::io::native::write(&other, std::slice::from_ref(&block), options).unwrap();
    let absorbed = ogeom::io::native::read_into(&mut model, &text).unwrap();
    // A box 1 by 2 by 3; the bound is rounding.
    assert!((volume(&model, &absorbed.shapes[0]) - 6.0).abs() < 1e-9);
    assert!((volume(&model, &plate) - want).abs() <= want * 1e-12);
}

/// The plate of 582 faces, built one cut per hole: its model compacts to
/// what the plate reaches, clones in a few milliseconds, and cuts as the
/// uncompacted model does.
#[test]
#[ignore = "heavy"]
fn a_582_face_plate_compacts_to_what_it_reaches() {
    let (mut model, plate, _) = holed_plate(24);
    let mut whole = model.clone();
    let before = faces(&model, &plate);
    assert_eq!(before.len(), 582);

    model
        .retain_reachable(std::slice::from_ref(&plate))
        .unwrap();
    assert!(model.node_count() < 20_000, "{} nodes", model.node_count());
    let least = (0..5)
        .map(|_| {
            let started = std::time::Instant::now();
            let copy = model.clone();
            let elapsed = started.elapsed().as_secs_f64();
            drop(copy);
            elapsed
        })
        .fold(f64::INFINITY, f64::min);
    assert!(least < 0.020, "a clone takes {:.1} ms", least * 1e3);
    assert_eq!(faces(&model, &plate), before);

    let after = corner_hole(&mut model, &plate);
    let want = corner_hole(&mut whole, &plate);
    assert!(after.2 && want.2, "both cuts are valid solids");
    assert_eq!(after.1, want.1);
    assert!(want.3 > 0, "the cut passes faces of the plate through");
    assert_eq!(after.3, want.3, "the result holds what it did");
    // The same cut on the same faces; the bound is rounding.
    assert!(
        (after.0 - want.0).abs() <= want.0 * 1e-9,
        "{after:?} {want:?}"
    );
}
