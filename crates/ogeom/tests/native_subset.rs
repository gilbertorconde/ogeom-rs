//! A written subset is a subset: `write(model, roots)` carries the closure
//! of `roots`, not the model.
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
/// as that box, while both roots still round-trip whole.
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
/// every root writes the whole model, handles unrenumbered.
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

/// A box cut by an `n` x `n` grid of through holes, one cut each, in a model
/// that keeps every intermediate plate's provenance.
fn holed_plate(n: u32) -> (Model, ogeom::topo::Shape) {
    let mut model = Model::new();
    let size = 10.0f64.mul_add(f64::from(n), 10.0);
    let mut plate = ogeom::algo::make_box(&mut model, Frame::WORLD, (size, size, 10.0), T)
        .unwrap()
        .shape;
    for i in 0..n {
        for j in 0..n {
            let at = Point::new(
                10.0f64.mul_add(f64::from(i), 10.0),
                10.0f64.mul_add(f64::from(j), 10.0),
                -5.0,
            );
            let hole =
                ogeom::algo::make_cylinder(&mut model, Frame::WORLD.with_origin(at), 1.5, 20.0, T)
                    .unwrap()
                    .shape;
            plate = ogeom::boolean::cut(&mut model, &plate, &hole, T)
                .unwrap()
                .shape;
        }
    }
    (model, plate)
}

/// Every sub-shape of `shape`, itself included.
fn every_sub_shape(model: &Model, shape: &ogeom::topo::Shape) -> Vec<ogeom::topo::Shape> {
    use ogeom::topo::ShapeType;
    let mut all = vec![shape.clone()];
    for kind in [
        ShapeType::Shell,
        ShapeType::Face,
        ShapeType::Wire,
        ShapeType::Edge,
        ShapeType::Vertex,
    ] {
        all.extend(ogeom::topo::explore_unique(model, shape, kind).unwrap());
    }
    all
}

fn faces(model: &Model, shape: &ogeom::topo::Shape) -> usize {
    ogeom::topo::explore_unique(model, shape, ogeom::topo::ShapeType::Face)
        .unwrap()
        .len()
}

/// Read back, the plate is the same solid and every identity its shapes
/// carry still resolves its entry.
fn assert_read_back_whole(model: &Model, plate: &ogeom::topo::Shape, text: &str) {
    let (restored, roots) = read(text).unwrap();
    let back = &roots[0];
    assert!(ogeom::algo::check(&restored, back, T).unwrap().is_valid());
    assert_eq!(faces(&restored, back), faces(model, plate));
    let (v0, v1) = (volume(model, plate), volume(&restored, back));
    assert!((v0 - v1).abs() <= 1e-9 * v0, "volume {v1} against {v0}");
    let ids = identified(&restored, back);
    assert_eq!(ids, identified(model, plate));
    for id in ids {
        assert!(
            restored.provenance().get(id).is_some(),
            "{id:?} does not resolve"
        );
    }
}

/// The identities the nodes under `shape` carry, ascending.
fn identified(model: &Model, shape: &ogeom::topo::Shape) -> Vec<ogeom::core::EntityId> {
    let mut ids: Vec<_> = every_sub_shape(model, shape)
        .iter()
        .filter_map(|s| model.identity_of(s))
        .collect();
    ids.sort_unstable();
    ids
}

/// A snapshot keeps the provenance its shapes reach and nothing else: the
/// entries of the identities it carries and the ancestry those name. Every
/// lineage query on the written shapes answers as it did in the model.
#[test]
fn a_snapshot_keeps_only_the_provenance_its_shapes_reach() {
    let (model, plate) = holed_plate(3);
    let text = write(
        &model,
        std::slice::from_ref(&plate),
        WriteOptions::default(),
    )
    .unwrap();

    // What a reader of the plate can reach, walked here from the model.
    let mut reach = ogeom::core::FastSet::default();
    let mut stack: Vec<_> = every_sub_shape(&model, &plate)
        .iter()
        .filter_map(|s| model.identity_of(s))
        .collect();
    while let Some(id) = stack.pop() {
        if reach.insert(id) {
            stack.extend(model.provenance().get(id).unwrap().inputs().iter().copied());
        }
    }
    let written = text.lines().filter(|l| l.starts_with("entity ")).count();
    assert_eq!(written, reach.len());
    assert!(
        written * 3 < model.provenance().len(),
        "{written} entries written of {}",
        model.provenance().len()
    );
    // Ids keep their numbers, so the table has gaps, which version 3 reads.
    assert!(text.starts_with("ogeom 3\n"));

    assert_read_back_whole(&model, &plate, &text);
    let (restored, roots) = read(&text).unwrap();
    for shape in every_sub_shape(&model, &plate) {
        let Some(id) = model.identity_of(&shape) else {
            continue;
        };
        let there = restored.shape_of(id).unwrap();
        assert_eq!(restored.roots_of(&there), model.roots_of(&shape));
    }
    // An entry nothing reaches is gone from the read model, its id still
    // issued.
    let dropped = (1..=model.provenance().len() as u64)
        .filter_map(ogeom::core::EntityId::from_raw)
        .find(|id| !reach.contains(id) && id.get() < restored.provenance().len() as u64)
        .unwrap();
    assert!(restored.provenance().get(dropped).is_none());
    // And a read model writes the same bytes again.
    assert_eq!(
        write(&restored, &roots, WriteOptions::default()).unwrap(),
        text
    );

    // Absorbed into a live model, the kept entries land under their new
    // ids and the dropped ones stay dropped.
    let mut live = Model::new();
    ogeom::algo::make_box(&mut live, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
    let absorbed = ogeom::io::native::read_into(&mut live, &text).unwrap();
    let ids = identified(&live, &absorbed.shapes[0]);
    assert!(!ids.is_empty());
    for id in ids {
        assert!(live.provenance().get(id).is_some());
    }
    assert!(live.provenance().get(absorbed.entities[&dropped]).is_none());
}

/// A whole model's table has no gaps, and is written at the version every
/// earlier reader follows.
#[test]
fn a_gapless_table_is_written_at_version_two() {
    let mut model = Model::new();
    ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T).unwrap();
    let text = write(&model, &[], WriteOptions::default()).unwrap();
    assert!(text.starts_with("ogeom 2\n"));
}

/// The plate of a 24 x 24 grid of holes, each cut on its own, carries the
/// provenance of all 576 intermediate plates. Its snapshot carries its own.
#[test]
#[ignore = "heavy"]
fn a_582_face_plate_snapshot_is_small_and_quick() {
    let (model, plate) = holed_plate(24);
    assert_eq!(faces(&model, &plate), 582);
    let options = WriteOptions::default();
    // The fastest of ten, so a busy machine does not decide it.
    let mut fastest = std::time::Duration::MAX;
    let mut text = String::new();
    for _ in 0..10 {
        let started = ogeom::core::clock::Instant::now();
        text = write(&model, std::slice::from_ref(&plate), options).unwrap();
        fastest = fastest.min(started.elapsed());
    }
    assert!(text.len() < 2_000_000, "{} bytes", text.len());
    assert!(fastest.as_millis() < 50, "written in {fastest:?}");
    assert_read_back_whole(&model, &plate, &text);
}
