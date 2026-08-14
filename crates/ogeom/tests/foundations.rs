//! The §1 foundations M7 owed: cancellation that stops a long operation at
//! its next checkpoint, and parallelism whose answer is bit-identical at any
//! thread count.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::{OgeomError, Tolerances, parallel, progress};
use ogeom::math::{Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

/// One deterministic part with enough faces to be worth threads: a box, a
/// drum and a torus, tessellated together under one compound-free model.
fn build(model: &mut Model) -> Vec<Shape> {
    let a = ogeom::algo::make_box(model, Frame::WORLD, (30.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let b = ogeom::algo::make_cylinder(model, Frame::WORLD, 6.0, 25.0, T)
        .unwrap()
        .shape;
    let c = ogeom::algo::make_torus(model, Frame::WORLD, 12.0, 3.0, T)
        .unwrap()
        .shape;
    vec![a, b, c]
}

/// The same model tessellated at one thread and at four serializes to the
/// same bytes — scheduling cannot reach the answer.
#[test]
fn parallel_tessellation_is_bit_identical() {
    let write = |threads: usize| {
        parallel::set_threads(threads);
        let mut model = Model::new();
        let shapes = build(&mut model);
        for shape in &shapes {
            ogeom::mesh::tessellate(&mut model, shape, Deflection::default(), T).unwrap();
        }
        parallel::set_threads(0);
        ogeom::io::native::write(&model, &shapes, ogeom::io::native::WriteOptions::default())
            .unwrap()
    };
    let serial = write(1);
    let threaded = write(4);
    assert_eq!(
        serial, threaded,
        "the native bytes must not depend on the thread count"
    );
}

/// A pre-cancelled watch stops the tessellation at its first checkpoint,
/// and the error says cancelled — not a partial result, not a stall.
#[test]
fn a_cancelled_watch_stops_tessellation() {
    let mut model = Model::new();
    let shapes = build(&mut model);
    let watch = progress::Watch::new();
    watch.canceller().cancel();
    let outcome = progress::watched(&watch, || {
        ogeom::mesh::tessellate(&mut model, &shapes[0], Deflection::default(), T)
    });
    assert!(matches!(outcome, Err(OgeomError::Cancelled)));
}

/// The stage names of a boolean reach the watch's sink in pipeline order,
/// and the cut still measures — a sink is a listener, not a participant.
#[test]
fn a_boolean_reports_its_stages() {
    use std::sync::{Arc, Mutex};
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let drill_frame = Frame::new(
        ogeom::math::Point::new(10.0, 10.0, -1.0),
        ogeom::math::Direction::Z,
        ogeom::math::Direction::X,
        T,
    )
    .unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, drill_frame, 3.0, 12.0, T)
        .unwrap()
        .shape;

    let heard: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
    let record = Arc::clone(&heard);
    let watch = progress::Watch::with_sink(move |name| {
        record.lock().unwrap().push(name.to_owned());
    });
    let cut = progress::watched(&watch, || {
        ogeom::boolean::cut(&mut model, &block, &drill, T)
    })
    .unwrap()
    .shape;

    let stages = heard.lock().unwrap().clone();
    assert!(
        stages.iter().any(|s| s == "boolean: gather")
            && stages.iter().any(|s| s == "boolean: split"),
        "the pipeline names its stages: {stages:?}"
    );
    let volume = ogeom::algo::volume_properties(&model, &cut, Deflection::default(), T)
        .unwrap()
        .mass;
    let expected = 20.0f64 * 20.0 * 10.0 - core::f64::consts::PI * 9.0 * 10.0;
    assert!(
        (volume - expected).abs() / expected < 0.01,
        "{volume} vs {expected}"
    );
}

/// The whole-shape mesh is the same at one thread and at four, to the bit.
///
/// `tessellate` writes its meshes into the model, so the test above can
/// compare serialized bytes. `triangulate` hands its mesh straight back, and
/// what has to be identical is that mesh: every position, every index, in
/// order. A face landing in the wrong slot, or a piece appended out of turn,
/// would show here and nowhere else.
#[test]
fn parallel_triangulation_returns_the_same_mesh() {
    let mesh_at = |threads: usize| {
        parallel::set_threads(threads);
        let mut model = Model::new();
        let shapes = build(&mut model);
        let meshes: Vec<_> = shapes
            .iter()
            .map(|shape| ogeom::mesh::triangulate(&model, shape, Deflection::default(), T).unwrap())
            .collect();
        parallel::set_threads(0);
        meshes
    };
    let serial = mesh_at(1);
    let threaded = mesh_at(4);

    assert_eq!(serial.len(), threaded.len());
    for (one, many) in serial.iter().zip(&threaded) {
        assert_eq!(
            one.triangles, many.triangles,
            "the triangles must not depend on the thread count"
        );
        assert_eq!(
            one.positions.len(),
            many.positions.len(),
            "the vertex count must not depend on the thread count"
        );
        for (a, b) in one.positions.iter().zip(&many.positions) {
            assert!(
                a.x.to_bits() == b.x.to_bits()
                    && a.y.to_bits() == b.y.to_bits()
                    && a.z.to_bits() == b.z.to_bits(),
                "a vertex moved between thread counts: {a:?} against {b:?}"
            );
        }
    }
}

/// A boolean's answer does not depend on how many threads paved its sections.
///
/// The filler measures each section against the boundaries it crosses, and
/// that measuring now runs in parallel. What must not move is where the
/// crossings land: a pave arriving in a different order would split an edge
/// at the same places in a different sequence, and the rebuilt solid would
/// differ in its vertex numbering even where its geometry agreed. Serialized
/// bytes catch exactly that.
#[test]
fn a_parallel_boolean_is_bit_identical() {
    let cut_at = |threads: usize| {
        parallel::set_threads(threads);
        let mut model = Model::new();
        let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let mut solid = block;
        for (x, y) in [(5.0, 5.0), (15.0, 5.0), (5.0, 15.0)] {
            let frame = Frame::new(
                Point::new(x, y, -1.0),
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
                T,
            )
            .unwrap();
            let drill = ogeom::algo::make_cylinder(&mut model, frame, 2.0, 12.0, T)
                .unwrap()
                .shape;
            solid = ogeom::boolean::cut(&mut model, &solid, &drill, T)
                .unwrap()
                .shape;
        }
        parallel::set_threads(0);
        let shapes = vec![solid];
        ogeom::io::native::write(&model, &shapes, ogeom::io::native::WriteOptions::default())
            .unwrap()
    };
    let serial = cut_at(1);
    let threaded = cut_at(4);
    assert_eq!(
        serial, threaded,
        "the cut solid's bytes must not depend on the thread count"
    );
}
