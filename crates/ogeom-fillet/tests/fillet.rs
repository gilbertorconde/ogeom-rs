//! The rolling-ball fillet, where the ball's envelope is a cylinder.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Frame, Point};
use ogeom_topo::{Filter, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

/// The top edge of the box along y at x = 2, z = 2.
fn top_edge(model: &ogeom_topo::Model, block: &ogeom_topo::Shape) -> ogeom_topo::Shape {
    explore(model, block, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(model, e)
                .unwrap()
                .is_some_and(|(a, b)| {
                    let p = |v: &ogeom_topo::Shape| {
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
        .expect("the box has that edge")
}

#[test]
fn a_filleted_box_edge_gains_a_tangent_cylinder() {
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let edge = top_edge(&model, &block.shape);

    let radius = 0.5;
    let result = ogeom_fillet::fillet_edge(&mut model, &block.shape, &edge, radius, T).unwrap();

    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // The removed sliver is the corner square minus the quarter disc, run
    // along the edge. The mesh volume converges to it with the chord.
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    let props = ogeom_algo::volume_properties(&model, &result.shape, fine, T).unwrap();
    let exact = 8.0 - (radius * radius - core::f64::consts::PI * radius * radius / 4.0) * 2.0;
    assert!(
        (props.mass - exact).abs() < 1e-3,
        "fillet volume {} against {exact}",
        props.mass
    );

    // Seven faces: six of the box (two now trimmed back to the tangency
    // lines) plus the blend cylinder.
    let faces = explore(&model, &result.shape, Filter::OfType(ShapeType::Face)).unwrap();
    assert_eq!(faces.len(), 7);

    // The history knows the rounded edge is gone.
    assert!(result.history.is_deleted(&edge));

    // The blend passes through the arc midpoint: the ball's centre line runs
    // at (1.5, y, 1.5), and the surface lies half a radius further out along
    // the diagonal.
    let mid = Point::new(
        1.5 + radius / core::f64::consts::SQRT_2,
        1.0,
        1.5 + radius / core::f64::consts::SQRT_2,
    );
    let on_blend = faces.iter().any(|f| {
        ogeom_algo::classify_on_face(&model, f, mid, fine, T)
            .map(|c| c == ogeom_algo::Containment::In)
            .unwrap_or(false)
    });
    assert!(on_blend, "the blend face passes through the arc midpoint");

    // Tangency is the point of a fillet: where the blend meets the top face
    // there is no crease. The two surfaces share their normal along the
    // tangency line at (1.5, y, 2.0) — the cylinder's radial direction there
    // is +z, the plane's normal exactly +z.
    let tangent = Point::new(1.5, 1.0, 2.0);
    let touching = faces
        .iter()
        .filter(|f| {
            ogeom_algo::classify_on_face(&model, f, tangent, fine, T)
                .map(|c| c != ogeom_algo::Containment::Out)
                .unwrap_or(false)
        })
        .count();
    assert!(
        touching >= 2,
        "the tangency line lies on both the top face and the blend, found {touching}"
    );
}

#[test]
fn a_cylinder_rim_gains_a_toroidal_blend() {
    let mut model = ogeom_topo::Model::new();
    let (radius, height, blend) = (1.0, 2.0, 0.3);
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, radius, height, T).unwrap();

    // The top rim: the circular edge at z = height.
    let edge = explore(&model, &drum.shape, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(&model, e)
                .unwrap()
                .is_some_and(|(a, _)| {
                    model
                        .node(&a)
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                        .is_some_and(|p| (p.z - height).abs() < 1e-9 && p.x.hypot(p.y) > 0.5)
                })
        })
        .expect("the cylinder has its top rim");

    let result = ogeom_fillet::fillet_edge(&mut model, &drum.shape, &edge, blend, T).unwrap();

    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // Pappus over the meridian cusp: the removed ring is the corner square
    // minus the quarter disc, each weighted by its distance from the axis.
    let pi = core::f64::consts::PI;
    let removed = 2.0
        * pi
        * (radius * blend * blend
            - blend * blend * blend / 2.0
            - (radius - blend) * pi * blend * blend / 4.0
            - blend * blend * blend / 3.0);
    let exact = pi * radius * radius * height - removed;
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    let props = ogeom_algo::volume_properties(&model, &result.shape, fine, T).unwrap();
    assert!(
        (props.mass - exact).abs() < 2e-3,
        "rim fillet volume {} against {exact}",
        props.mass
    );

    // Four faces: the wall and both caps trimmed back, plus the torus.
    let faces = explore(&model, &result.shape, Filter::OfType(ShapeType::Face)).unwrap();
    assert_eq!(faces.len(), 4);

    assert!(result.history.is_deleted(&edge));

    // The blend passes through the meridian midpoint of the quarter tube.
    let mid = Point::new(
        (radius - blend) + blend / core::f64::consts::SQRT_2,
        0.0,
        (height - blend) + blend / core::f64::consts::SQRT_2,
    );
    let on_blend = faces.iter().any(|f| {
        ogeom_algo::classify_on_face(&model, f, mid, fine, T)
            .map(|c| c != ogeom_algo::Containment::Out)
            .unwrap_or(false)
    });
    assert!(
        on_blend,
        "the blend face passes through the tube's meridian midpoint"
    );
}

#[test]
fn a_rim_fillet_that_swallows_the_axis_is_refused() {
    let mut model = ogeom_topo::Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 1.0, 2.0, T).unwrap();
    let edge = explore(&model, &drum.shape, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(&model, e)
                .unwrap()
                .is_some_and(|(a, _)| {
                    model
                        .node(&a)
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                        .is_some_and(|p| (p.z - 2.0).abs() < 1e-9 && p.x.hypot(p.y) > 0.5)
                })
        })
        .expect("the cylinder has its top rim");
    assert!(ogeom_fillet::fillet_edge(&mut model, &drum.shape, &edge, 1.0, T).is_err());
}

#[test]
fn a_variable_radius_fillet_widens_along_its_edge() {
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let edge = top_edge(&model, &block.shape);

    let (r0, r1) = (0.3, 0.6);
    let result =
        ogeom_fillet::fillet_edge_variable(&mut model, &block.shape, &edge, r0, r1, T).unwrap();

    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // The removed sliver integrates the section area along the edge: for the
    // linear law, L (1 - pi/4) times the integral of r(u) squared.
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    let integral = (r1 * r1 * r1 - r0 * r0 * r0) / (3.0 * (r1 - r0));
    let exact = 8.0 - 2.0 * (1.0 - core::f64::consts::PI / 4.0) * integral;
    let props = ogeom_algo::volume_properties(&model, &result.shape, fine, T).unwrap();
    assert!(
        (props.mass - exact).abs() < 2e-3,
        "variable fillet volume {} against {exact}",
        props.mass
    );

    assert_eq!(
        explore(&model, &result.shape, Filter::OfType(ShapeType::Face))
            .unwrap()
            .len(),
        7
    );
    assert!(result.history.is_deleted(&edge));

    // The blend passes through the mid-section's arc midpoint, where the
    // radius is the mean of the two ends.
    let r = f64::midpoint(r0, r1);
    let mid = Point::new(
        (2.0 - r) + r / core::f64::consts::SQRT_2,
        1.0,
        (2.0 - r) + r / core::f64::consts::SQRT_2,
    );
    let faces = explore(&model, &result.shape, Filter::OfType(ShapeType::Face)).unwrap();
    let on_blend = faces.iter().any(|f| {
        ogeom_algo::classify_on_face(&model, f, mid, fine, T)
            .map(|c| c != ogeom_algo::Containment::Out)
            .unwrap_or(false)
    });
    assert!(
        on_blend,
        "the blend passes through the mean-radius midpoint"
    );
}

#[test]
fn a_concave_edge_gains_an_additive_fillet() {
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let seat = Frame::new(
        Point::new(1.0, -0.5, 1.0),
        ogeom_math::Direction::Z,
        ogeom_math::Direction::X,
        T,
    )
    .unwrap();
    let notch = ogeom_algo::make_box(&mut model, seat, (2.0, 3.0, 2.0), T).unwrap();
    let cut = ogeom_bool::cut(&mut model, &block.shape, &notch.shape, T).unwrap();

    // The re-entrant edge of the L: along y at x = 1, z = 1.
    let edge = explore(&model, &cut.shape, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(&model, e)
                .unwrap()
                .is_some_and(|(a, b)| {
                    let p = |v: &ogeom_topo::Shape| {
                        model
                            .node(v)
                            .and_then(|n| n.data().as_vertex().map(|d| d.point))
                            .unwrap()
                    };
                    let (pa, pb) = (p(&a), p(&b));
                    (pa.x - 1.0).abs() < 1e-9
                        && (pa.z - 1.0).abs() < 1e-9
                        && (pb.x - 1.0).abs() < 1e-9
                        && (pb.z - 1.0).abs() < 1e-9
                })
        })
        .expect("the L has its re-entrant edge");

    let radius = 0.25;
    let result = ogeom_fillet::fillet_edge(&mut model, &cut.shape, &edge, radius, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // The L keeps its volume plus the added sliver along the re-entrant
    // edge: the corner square minus the quarter disc, run along its length.
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    let added = radius * radius * (1.0 - core::f64::consts::PI / 4.0) * 2.0;
    let exact = 6.0 + added;
    let props = ogeom_algo::volume_properties(&model, &result.shape, fine, T).unwrap();
    assert!(
        (props.mass - exact).abs() < 1e-3,
        "concave fillet volume {} against {exact}",
        props.mass
    );
    assert!(result.history.is_deleted(&edge));
}
