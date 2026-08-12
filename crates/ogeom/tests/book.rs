//! The guide's code, run. Every example in `docs/book/` is included from
//! this file by anchor, so a snippet that stops compiling — or stops being
//! true — fails the build instead of lying in the documentation.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The edge of `shape` whose midpoint is nearest `near`.
fn edge_near(model: &Model, shape: &Shape, near: Point) -> Shape {
    use ogeom::geom::Curve3d as _;
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| {
            let mid = |e: &Shape| {
                let data = model.node(e).unwrap().data().as_edge().unwrap();
                let ogeom::topo::EdgeRepr::Curve3d { curve, range, .. } = data.curve3d().unwrap()
                else {
                    unreachable!()
                };
                model
                    .geometry()
                    .curve(*curve)
                    .unwrap()
                    .point_at(f64::midpoint(range.0, range.1), T)
                    .unwrap()
                    .distance(near)
            };
            mid(a)
                .partial_cmp(&mid(b))
                .unwrap_or(core::cmp::Ordering::Equal)
        })
        .expect("some edge")
}

#[test]
fn first_solid() {
    // ANCHOR: first_solid
    let mut model = Model::new();
    let tol = Tolerances::millimetres();

    // A 20×20×10 block, and a Ø8 hole through its middle.
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), tol)
        .unwrap()
        .shape;
    let axis = Frame::new(Point::new(10.0, 10.0, 0.0), Direction::Z, Direction::X, tol).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, axis, 4.0, 10.0, tol)
        .unwrap()
        .shape;
    let part = ogeom::boolean::cut(&mut model, &block, &drill, tol)
        .unwrap()
        .shape;

    // The result is measured, not assumed: volume against the closed form.
    let volume = ogeom::algo::volume_properties(&model, &part, Deflection::default(), tol)
        .unwrap()
        .mass;
    let exact = 20.0 * 20.0 * 10.0 - core::f64::consts::PI * 4.0 * 4.0 * 10.0;
    assert!((volume - exact).abs() / exact < 0.01);
    // ANCHOR_END: first_solid
}

#[test]
fn fillet_an_edge() {
    let mut model = Model::new();
    // ANCHOR: fillet_an_edge
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 12.0), T)
        .unwrap()
        .shape;

    // Pick the top edge along y = 0 and round it at radius 2.
    let edge = edge_near(&model, &block, Point::new(20.0, 0.0, 12.0));
    let rounded = ogeom::fillet::fillet_edge(&mut model, &block, &edge, 2.0, T)
        .unwrap()
        .shape;

    // A fillet removes the square corner and leaves the quarter cylinder:
    // ΔV = (1 − π/4)·r²·length, exactly.
    let volume = ogeom::algo::volume_properties(&model, &rounded, Deflection::default(), T)
        .unwrap()
        .mass;
    let exact = 40.0 * 30.0 * 12.0 - (1.0 - core::f64::consts::FRAC_PI_4) * 4.0 * 40.0;
    assert!((volume - exact).abs() / exact < 0.01);
    // ANCHOR_END: fillet_an_edge
}

#[test]
fn a_part_survives_step() {
    let mut model = Model::new();
    // ANCHOR: step_roundtrip
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 10.0, 5.0), T)
        .unwrap()
        .shape;

    // Exchange works on documents: a model plus product structure, colours,
    // PMI, views. A bare part is a document with one product.
    let mut document = ogeom::doc::Document::over(model);
    document.add_part("block", block);

    let text = ogeom::io::write_step(&document, T).unwrap();
    let import = ogeom::io::read_step(&text, T).unwrap();

    // What came back is the same solid, measured.
    let back = &import.document;
    let root = back.roots()[0];
    let occurrence = &back.occurrences_of(root).unwrap()[0];
    let volume =
        ogeom::algo::volume_properties(back.model(), &occurrence.shape, Deflection::default(), T)
            .unwrap()
            .mass;
    assert!((volume - 1000.0).abs() / 1000.0 < 0.01);
    // ANCHOR_END: step_roundtrip
}

#[test]
fn refused_by_name() {
    let mut model = Model::new();
    // ANCHOR: refused_by_name
    // A cylinder of zero radius is not a small cylinder; it is a mistake,
    // and the kernel says which one instead of producing garbage geometry.
    let err = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 0.0, 10.0, T).unwrap_err();
    assert!(err.to_string().contains("cylinder radius"));
    // ANCHOR_END: refused_by_name
}
