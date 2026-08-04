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
    let mut samples: Vec<(f64, f64)> = Vec::new();
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
            samples.push((p.x, p.y));
        }
    }
    let mut doubled = 0.0;
    for i in 0..samples.len() {
        let (p, q) = (samples[i], samples[(i + 1) % samples.len()]);
        doubled += p.0.mul_add(q.1, -(q.0 * p.1));
    }
    (doubled / 2.0).abs()
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
