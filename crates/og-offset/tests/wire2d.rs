//! The 2D wire offset, measured against Minkowski's arithmetic.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use og_core::Tolerances;
use og_geom::Curve3d as _;
use og_math::Point;
use og_offset::Join;
use og_topo::{EdgeRepr, Filter, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

fn rectangle(model: &mut og_topo::Model) -> og_topo::Shape {
    let corners = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(4.0, 0.0, 0.0),
        Point::new(4.0, 3.0, 0.0),
        Point::new(0.0, 3.0, 0.0),
    ];
    og_algo::make_polygon(model, &corners, true, T)
        .unwrap()
        .shape
}

/// The wire's enclosed area, by shoelace over a fine sampling.
fn area(model: &og_topo::Model, wire: &og_topo::Shape) -> f64 {
    let mut samples: Vec<(f64, f64, f64)> = Vec::new();
    for edge in explore(model, wire, Filter::OfType(ShapeType::Edge)).unwrap() {
        let node = model.node(&edge).unwrap();
        let data = node.data().as_edge().unwrap();
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            panic!("an edge without a curve");
        };
        let geometry = model.geometry().curve(*curve).unwrap().clone();
        let (a, b) = *range;
        let reversed = edge.orientation() == og_topo::Orientation::Reversed;
        for i in 0..256 {
            let f = f64::from(i) / 256.0;
            let t = if reversed {
                b - (b - a) * f
            } else {
                a + (b - a) * f
            };
            let p = geometry.point_at(t, T).unwrap();
            samples.push((p.x, p.y, p.z));
        }
    }
    // The vector area: half the norm of the cross-product sum, which works
    // in whatever plane the wire lies.
    let mut sum = og_math::Vector::new(0.0, 0.0, 0.0);
    for i in 0..samples.len() {
        let (p, q) = (samples[i], samples[(i + 1) % samples.len()]);
        sum += og_math::Vector::new(
            p.1 * q.2 - p.2 * q.1,
            p.2 * q.0 - p.0 * q.2,
            p.0 * q.1 - p.1 * q.0,
        );
    }
    sum.magnitude() / 2.0
}

#[test]
fn an_outward_offset_with_arc_joins_grows_by_minkowski() {
    let mut model = og_topo::Model::new();
    let wire = rectangle(&mut model);
    let w = 0.5;
    let result = og_offset::offset_wire(&mut model, &wire, w, Join::Arc, T).unwrap();

    assert!(og_algo::is_wire_closed(&model, &result.shape, T).unwrap());
    let edges = explore(&model, &result.shape, Filter::OfType(ShapeType::Edge)).unwrap();
    assert_eq!(edges.len(), 8, "four sides and four corner arcs");

    let expected = 12.0 + 14.0 * w + core::f64::consts::PI * w * w;
    let measured = area(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 1e-3,
        "area {measured} against {expected}"
    );
}

#[test]
fn an_outward_offset_with_intersection_joins_is_the_bigger_rectangle() {
    let mut model = og_topo::Model::new();
    let wire = rectangle(&mut model);
    let w = 0.5;
    let result = og_offset::offset_wire(&mut model, &wire, w, Join::Intersection, T).unwrap();

    assert!(og_algo::is_wire_closed(&model, &result.shape, T).unwrap());
    let edges = explore(&model, &result.shape, Filter::OfType(ShapeType::Edge)).unwrap();
    assert_eq!(edges.len(), 4, "the corners extend to sharp meetings");

    let expected = 5.0 * 4.0;
    let measured = area(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 1e-9,
        "area {measured} against {expected}"
    );
}

#[test]
fn an_inward_offset_trims_to_the_smaller_rectangle() {
    let mut model = og_topo::Model::new();
    let wire = rectangle(&mut model);
    let result = og_offset::offset_wire(&mut model, &wire, -0.5, Join::Arc, T).unwrap();

    assert!(og_algo::is_wire_closed(&model, &result.shape, T).unwrap());
    let edges = explore(&model, &result.shape, Filter::OfType(ShapeType::Edge)).unwrap();
    assert_eq!(edges.len(), 4, "inward corners trim rather than join");

    let expected = 3.0 * 2.0;
    let measured = area(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 1e-9,
        "area {measured} against {expected}"
    );
}

#[test]
fn offsetting_the_offset_composes_like_one_bigger_offset() {
    // The second pass exercises arc pieces and tangent-continuous corners:
    // a rounded rectangle's sides meet its corner arcs with no gap at all.
    let mut model = og_topo::Model::new();
    let wire = rectangle(&mut model);
    let first = og_offset::offset_wire(&mut model, &wire, 0.5, Join::Arc, T).unwrap();
    let second = og_offset::offset_wire(&mut model, &first.shape, 0.5, Join::Arc, T).unwrap();

    assert!(og_algo::is_wire_closed(&model, &second.shape, T).unwrap());
    let expected = 12.0 + 14.0 * 1.0 + core::f64::consts::PI;
    let measured = area(&model, &second.shape);
    assert!(
        (measured - expected).abs() < 1e-3,
        "area {measured} against {expected}"
    );
}

#[test]
fn an_offset_that_consumes_the_wire_is_refused() {
    let mut model = og_topo::Model::new();
    let wire = rectangle(&mut model);
    assert!(og_offset::offset_wire(&mut model, &wire, -1.6, Join::Arc, T).is_err());
}

#[test]
fn history_reports_the_edit_edge_by_edge() {
    let mut model = og_topo::Model::new();
    let wire = rectangle(&mut model);
    let inputs = explore(&model, &wire, Filter::OfType(ShapeType::Edge)).unwrap();
    let result = og_offset::offset_wire(&mut model, &wire, 0.5, Join::Arc, T).unwrap();
    for edge in &inputs {
        assert_eq!(
            result.history.modified(edge).len(),
            1,
            "every side moved to exactly one offset side"
        );
    }
    assert_eq!(result.history.modified(&wire).len(), 1);
}

/// The wire's total enclosed area works for compounds too: sum per wire.
fn area_of_all(model: &og_topo::Model, shape: &og_topo::Shape) -> (usize, f64) {
    let wires = if model.kind_of(shape).unwrap() == og_topo::ShapeType::Wire {
        vec![shape.clone()]
    } else {
        explore(model, shape, Filter::OfType(og_topo::ShapeType::Wire)).unwrap()
    };
    let total = wires.iter().map(|w| area(model, w)).sum();
    (wires.len(), total)
}

#[test]
fn an_open_path_offsets_into_its_rounded_outline() {
    let mut model = og_topo::Model::new();
    let path = og_algo::make_polygon(
        &mut model,
        &[Point::new(0.0, 0.0, 0.0), Point::new(6.0, 0.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let w = 0.5;
    let result = og_offset::offset_wire(&mut model, &path, w, Join::Arc, T).unwrap();
    assert!(og_algo::is_wire_closed(&model, &result.shape, T).unwrap());
    let expected = 2.0 * w * 6.0 + core::f64::consts::PI * w * w;
    let measured = area(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 1e-3,
        "rounded outline area {measured} against {expected}"
    );
}

#[test]
fn an_open_path_offsets_into_its_square_outline() {
    let mut model = og_topo::Model::new();
    let path = og_algo::make_polygon(
        &mut model,
        &[Point::new(0.0, 0.0, 0.0), Point::new(6.0, 0.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let w = 0.5;
    let result = og_offset::offset_wire(&mut model, &path, w, Join::Intersection, T).unwrap();
    assert!(og_algo::is_wire_closed(&model, &result.shape, T).unwrap());
    let expected = 2.0 * w * (6.0 + 2.0 * w);
    let measured = area(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 1e-9,
        "square outline area {measured} against {expected}"
    );
}

#[test]
fn a_collapsing_inward_offset_resolves_into_its_islands() {
    // A dumbbell: two 4-squares joined by a thin neck. Offsetting inward
    // past the neck's half-width severs it, and the survivors are two loops.
    let mut model = og_topo::Model::new();
    let outline = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(4.0, 0.0, 0.0),
        Point::new(4.0, 1.5, 0.0),
        Point::new(5.0, 1.5, 0.0),
        Point::new(5.0, 0.0, 0.0),
        Point::new(9.0, 0.0, 0.0),
        Point::new(9.0, 4.0, 0.0),
        Point::new(5.0, 4.0, 0.0),
        Point::new(5.0, 2.5, 0.0),
        Point::new(4.0, 2.5, 0.0),
        Point::new(4.0, 4.0, 0.0),
        Point::new(0.0, 4.0, 0.0),
    ];
    let wire = og_algo::make_polygon(&mut model, &outline, true, T)
        .unwrap()
        .shape;
    let result = og_offset::offset_wire(&mut model, &wire, -0.8, Join::Arc, T).unwrap();
    let (loops, total) = area_of_all(&model, &result.shape);
    assert_eq!(loops, 2, "the severed neck leaves two islands");
    assert!(
        total > 9.0 && total < 13.0,
        "island area sum {total} within the plausible band"
    );
}
