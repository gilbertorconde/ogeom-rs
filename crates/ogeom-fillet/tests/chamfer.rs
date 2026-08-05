//! The chamfer: P4's opening stone, standing on M3's booleans.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Frame, Point};
use ogeom_topo::{Filter, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_chamfered_box_edge_loses_exactly_the_wedge() {
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();

    // The top edge along y at x = 2, z = 2.
    let edge = explore(&model, &block.shape, Filter::OfType(ShapeType::Edge))
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
                    (pa.x - 2.0).abs() < 1e-9
                        && (pa.z - 2.0).abs() < 1e-9
                        && (pb.x - 2.0).abs() < 1e-9
                        && (pb.z - 2.0).abs() < 1e-9
                })
        })
        .expect("the box has that edge");

    let distance = 0.5;
    let result = ogeom_fillet::chamfer_edge(&mut model, &block.shape, &edge, distance, T).unwrap();

    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let fine = ogeom_mesh::Deflection {
        chord: 1e-3,
        ..ogeom_mesh::Deflection::default()
    };
    let props = ogeom_algo::volume_properties(&model, &result.shape, fine, T).unwrap();
    let exact = 8.0 - distance * distance / 2.0 * 2.0;
    assert!(
        (props.mass - exact).abs() < 1e-9,
        "chamfer volume {} against {exact}",
        props.mass
    );

    // Seven faces: six of the box (two now notched, two trimmed) plus the
    // bevel.
    let faces = explore(&model, &result.shape, Filter::OfType(ShapeType::Face)).unwrap();
    assert_eq!(faces.len(), 7);

    // The history knows the beveled edge is gone.
    assert!(result.history.is_deleted(&edge));

    // The bevel face lies at 45 degrees: its centroid sits where the wedge's
    // hypotenuse ran.
    let mid = Point::new(2.0 - distance / 2.0, 1.0, 2.0 - distance / 2.0);
    let on_bevel = faces.iter().any(|f| {
        ogeom_algo::classify_on_face(&model, f, mid, fine, T)
            .map(|c| c == ogeom_algo::Containment::In)
            .unwrap_or(false)
    });
    assert!(on_bevel, "the bevel face passes through the wedge diagonal");
}
