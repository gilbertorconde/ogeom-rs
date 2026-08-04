//! The chamfer's other spellings, and the corner blends of the sketch plane.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use og_core::Tolerances;
use og_geom::Curve3d as _;
use og_math::{Frame, Point};
use og_topo::{Filter, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

/// The top edge of the box along y at x = 2, z = 2, and the top face z = 2.
fn top_edge_and_face(
    model: &og_topo::Model,
    block: &og_topo::Shape,
) -> (og_topo::Shape, og_topo::Shape) {
    let edge = explore(model, block, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            og_algo::edge_vertices(model, e)
                .unwrap()
                .is_some_and(|(a, b)| {
                    let p = |v: &og_topo::Shape| {
                        model
                            .node(v)
                            .and_then(|n| n.data().as_vertex().map(|d| d.point))
                            .unwrap()
                    };
                    let (pa, pb) = (p(&a), p(&b));
                    (pa.x - 2.0).abs() < 1e-9
                        && (pa.z - 2.0).abs() < 1e-9
                        && (pb.x - 2.0).abs() < 1e-9
                        && (pb.z - 2.0).abs() < 1e-9
                })
        })
        .expect("the box has that edge");
    let face = explore(model, block, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find(|f| {
            og_algo::classify_on_face(model, f, Point::new(1.0, 1.0, 2.0), fine(), T)
                .map(|c| c == og_algo::Containment::In)
                .unwrap_or(false)
        })
        .expect("the box has its top face");
    (edge, face)
}

fn fine() -> og_mesh::Deflection {
    og_mesh::Deflection {
        chord: 1e-3,
        ..og_mesh::Deflection::default()
    }
}

#[test]
fn an_asymmetric_chamfer_cuts_each_face_its_own_distance() {
    let mut model = og_topo::Model::new();
    let block = og_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let (edge, top) = top_edge_and_face(&model, &block.shape);

    let (d_top, d_side) = (0.5, 0.25);
    let result =
        og_fillet::chamfer_edge_distances(&mut model, &block.shape, &edge, &top, d_top, d_side, T)
            .unwrap();

    let diagnosis = og_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let props = og_algo::volume_properties(&model, &result.shape, fine(), T).unwrap();
    let exact = 8.0 - d_top * d_side / 2.0 * 2.0;
    assert!(
        (props.mass - exact).abs() < 1e-9,
        "asymmetric chamfer volume {} against {exact}",
        props.mass
    );

    // The named face carries the first distance: the bevel runs from
    // x = 2 − d_top on the top face to z = 2 − d_side on the side.
    let mid = Point::new(2.0 - d_top / 2.0, 1.0, 2.0 - d_side / 2.0);
    let faces = explore(&model, &result.shape, Filter::OfType(ShapeType::Face)).unwrap();
    let on_bevel = faces.iter().any(|f| {
        og_algo::classify_on_face(&model, f, mid, fine(), T)
            .map(|c| c == og_algo::Containment::In)
            .unwrap_or(false)
    });
    assert!(on_bevel, "the bevel runs between the two unequal contacts");
    assert!(result.history.is_deleted(&edge));
}

#[test]
fn a_distance_angle_chamfer_at_forty_five_degrees_is_the_symmetric_one() {
    let mut model = og_topo::Model::new();
    let block = og_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let (edge, top) = top_edge_and_face(&model, &block.shape);

    let result = og_fillet::chamfer_edge_angle(
        &mut model,
        &block.shape,
        &edge,
        &top,
        0.5,
        core::f64::consts::FRAC_PI_4,
        T,
    )
    .unwrap();
    let props = og_algo::volume_properties(&model, &result.shape, fine(), T).unwrap();
    let exact = 8.0 - 0.5 * 0.5 / 2.0 * 2.0;
    assert!(
        (props.mass - exact).abs() < 1e-9,
        "45-degree chamfer volume {} against {exact}",
        props.mass
    );
}

#[test]
fn a_distance_angle_chamfer_that_never_meets_the_other_face_is_refused() {
    let mut model = og_topo::Model::new();
    let block = og_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let (edge, top) = top_edge_and_face(&model, &block.shape);
    // Parallel to the reference face: the bevel never reaches the other.
    let refused = og_fillet::chamfer_edge_angle(
        &mut model,
        &block.shape,
        &edge,
        &top,
        0.5,
        core::f64::consts::PI - 1e-6,
        T,
    );
    assert!(refused.is_err());
}

/// A rectangle wire in the xy plane, and its corner vertex at (4, 3).
fn rectangle_with_corner(model: &mut og_topo::Model) -> (og_topo::Shape, og_topo::Shape) {
    let corners = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(4.0, 0.0, 0.0),
        Point::new(4.0, 3.0, 0.0),
        Point::new(0.0, 3.0, 0.0),
    ];
    let wire = og_algo::make_polygon(model, &corners, true, T)
        .unwrap()
        .shape;
    let vertex = explore(model, &wire, Filter::OfType(ShapeType::Vertex))
        .unwrap()
        .into_iter()
        .find(|v| {
            model
                .node(v)
                .and_then(|n| n.data().as_vertex().map(|d| d.point))
                .is_some_and(|p| (p.x - 4.0).abs() < 1e-9 && (p.y - 3.0).abs() < 1e-9)
        })
        .expect("the rectangle has that corner");
    (wire, vertex)
}

#[test]
fn a_filleted_corner_becomes_a_tangent_arc() {
    let mut model = og_topo::Model::new();
    let (wire, vertex) = rectangle_with_corner(&mut model);

    let radius = 0.5;
    let result = og_fillet::fillet_corner_2d(&mut model, &wire, &vertex, radius, T).unwrap();

    assert!(og_algo::is_wire_closed(&model, &result.shape, T).unwrap());
    let edges = explore(&model, &result.shape, Filter::OfType(ShapeType::Edge)).unwrap();
    assert_eq!(edges.len(), 5);

    // The arc's midpoint bulges toward the old corner from the centre at
    // (3.5, 2.5): tangency is what places it exactly.
    let arc_mid = Point::new(
        3.5 + radius / core::f64::consts::SQRT_2,
        2.5 + radius / core::f64::consts::SQRT_2,
        0.0,
    );
    let reached = edges.iter().any(|e| {
        let Some(node) = model.node(e) else {
            return false;
        };
        let Some(data) = node.data().as_edge() else {
            return false;
        };
        let Some(og_topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            return false;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            return false;
        };
        let mid = f64::midpoint(range.0, range.1);
        geometry
            .point_at(mid, T)
            .map(|p| p.is_within(arc_mid, 1e-9))
            .unwrap_or(false)
    });
    assert!(reached, "the arc passes through its tangent midpoint");

    // History: the corner is gone, the arc came from it, both edges trimmed.
    assert!(result.history.is_deleted(&vertex));
    assert_eq!(result.history.generated(&vertex).len(), 1);
}

#[test]
fn a_chamfered_corner_becomes_a_straight_cut() {
    let mut model = og_topo::Model::new();
    let (wire, vertex) = rectangle_with_corner(&mut model);

    let result = og_fillet::chamfer_corner_2d(&mut model, &wire, &vertex, 0.5, 1.0, T).unwrap();
    assert!(og_algo::is_wire_closed(&model, &result.shape, T).unwrap());
    let edges = explore(&model, &result.shape, Filter::OfType(ShapeType::Edge)).unwrap();
    assert_eq!(edges.len(), 5);

    // The cut runs from 0.5 back along the earlier edge to 1.0 along the
    // later; its midpoint says which side got which.
    let cut_mid = Point::new(f64::midpoint(4.0, 3.0), f64::midpoint(2.5, 3.0), 0.0);
    let reached = edges.iter().any(|e| {
        let Some(node) = model.node(e) else {
            return false;
        };
        let Some(data) = node.data().as_edge() else {
            return false;
        };
        let Some(og_topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            return false;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            return false;
        };
        let mid = f64::midpoint(range.0, range.1);
        geometry
            .point_at(mid, T)
            .map(|p| p.is_within(cut_mid, 1e-9))
            .unwrap_or(false)
    });
    assert!(reached, "the cut runs between its two unequal trims");
}

#[test]
fn a_corner_blend_that_consumes_a_whole_edge_is_refused() {
    let mut model = og_topo::Model::new();
    let (wire, vertex) = rectangle_with_corner(&mut model);
    assert!(og_fillet::fillet_corner_2d(&mut model, &wire, &vertex, 5.0, T).is_err());
}
