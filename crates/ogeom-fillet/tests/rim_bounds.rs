//! A carried round can rest only on the bounded circular edge.
#![allow(clippy::unwrap_used, reason = "test constructions")]

use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_shelf_round_refuses_contact_beyond_the_circular_edges_end() {
    let mut model = Model::new();
    let cylinder = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 10.0, 26.0, T)
        .unwrap()
        .shape;
    // Above the shelf the boss ends at y=1. Its circular foot consequently
    // ends there too, while the planar shelf continues past it to y=5.
    let frame = Frame::new(Point::new(-11.0, 1.0, 21.0), Direction::Z, Direction::X, T).unwrap();
    let tool = ogeom_algo::make_box(&mut model, frame, (22.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let cylinder = ogeom_bool::cut(&mut model, &cylinder, &tool, T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(0.0, -5.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let block = ogeom_algo::make_box(&mut model, frame, (12.6, 10.0, 21.0), T)
        .unwrap()
        .shape;
    let solid = ogeom_bool::fuse(&mut model, &cylinder, &block, T)
        .unwrap()
        .shape;
    let (unified, _) = ogeom_heal::unify_same_domain(&mut model, &solid, T).unwrap();
    let solid = unified.shape;
    let edges: Vec<_> = explore_unique(&model, &solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|edge| {
            ogeom_algo::edge_vertices(&model, edge)
                .unwrap()
                .is_some_and(|(a, b)| {
                    [a, b].iter().all(|v| {
                        let p = v
                            .transform(model.datums())
                            .unwrap()
                            .apply(model.node(v).unwrap().data().as_vertex().unwrap().point);
                        (p.x - 12.6).abs() < T.confusion() && (p.z - 21.0).abs() < T.confusion()
                    })
                })
        })
        .collect();
    assert_eq!(edges.len(), 1);
    // Where the two planes are free at y=2, their R=3 round has centre
    // (9.6,2,18). A transfer onto the continuation of the trimmed-away
    // circular rim would raise it and miss that exact round by 0.00456 mm.
    // The carried-rim family must refuse this extra support change.
    let refusal = ogeom_fillet::fillet_edges(&mut model, &solid, &edges, 3.0, T).unwrap_err();
    assert!(
        refusal
            .to_string()
            .contains("a carried shelf contact leaves its circular rim"),
        "{refusal}"
    );
}
