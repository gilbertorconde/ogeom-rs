//! The additive blends: concave edges, and the revolved seats in every sign.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Frame, Point};
use ogeom_topo::{Filter, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

fn fine() -> ogeom_mesh::Deflection {
    ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    }
}

fn volume(model: &ogeom_topo::Model, shape: &ogeom_topo::Shape) -> f64 {
    ogeom_algo::volume_properties(model, shape, fine(), T)
        .unwrap()
        .mass
}

/// An L: a box with a notch, and its re-entrant edge along y at (1, 1).
fn l_bracket(model: &mut ogeom_topo::Model) -> (ogeom_topo::Shape, ogeom_topo::Shape) {
    let block = ogeom_algo::make_box(model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let seat = Frame::new(
        Point::new(1.0, -0.5, 1.0),
        ogeom_math::Direction::Z,
        ogeom_math::Direction::X,
        T,
    )
    .unwrap();
    let notch = ogeom_algo::make_box(model, seat, (2.0, 3.0, 2.0), T).unwrap();
    let cut = ogeom_bool::cut(model, &block.shape, &notch.shape, T).unwrap();
    let edge = explore(model, &cut.shape, Filter::OfType(ShapeType::Edge))
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
                    (pa.x - 1.0).abs() < 1e-9
                        && (pa.z - 1.0).abs() < 1e-9
                        && (pb.x - 1.0).abs() < 1e-9
                        && (pb.z - 1.0).abs() < 1e-9
                })
        })
        .expect("the L has its re-entrant edge");
    (cut.shape, edge)
}

#[test]
fn a_concave_chamfer_fills_the_corner_with_its_bevel() {
    let mut model = ogeom_topo::Model::new();
    let (bracket, edge) = l_bracket(&mut model);
    let d = 0.25;
    let result = ogeom_fillet::chamfer_edge(&mut model, &bracket, &edge, d, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let exact = 6.0 + d * d / 2.0 * 2.0;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - exact).abs() < 1e-9,
        "concave chamfer volume {measured} against {exact}"
    );
    assert!(result.history.is_deleted(&edge));
}

/// The rim of a hole: the top edge of a tube's bore.
#[test]
fn a_hole_rim_gains_the_mirrored_toroidal_blend() {
    let mut model = ogeom_topo::Model::new();
    let outer = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 2.0, 2.0, T).unwrap();
    let bore = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 1.0, 2.0, T).unwrap();
    let tube = ogeom_bool::cut(&mut model, &outer.shape, &bore.shape, T).unwrap();

    // The bore's top rim: radius 1 at z = 2.
    let edge = explore(&model, &tube.shape, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(&model, e)
                .unwrap()
                .is_some_and(|(a, _)| {
                    model
                        .node(&a)
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                        .is_some_and(|p| {
                            (p.z - 2.0).abs() < 1e-9 && (p.x.hypot(p.y) - 1.0).abs() < 1e-6
                        })
                })
        })
        .expect("the tube has its bore rim");

    let r = 0.3;
    let result = ogeom_fillet::fillet_edge(&mut model, &tube.shape, &edge, r, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // Pappus over the mirrored meridian cusp at the bore radius.
    let pi = core::f64::consts::PI;
    let hole = 1.0;
    let removed = 2.0
        * pi
        * (hole * r * r + r * r * r / 2.0 - (hole + r) * pi * r * r / 4.0 + r * r * r / 3.0);
    let exact = pi * (4.0 - 1.0) * 2.0 - removed;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - exact).abs() < 2e-3,
        "hole rim fillet volume {measured} against {exact}"
    );
    assert!(result.history.is_deleted(&edge));
}

/// The base of a boss: the concave circular seat where a cylinder stands on
/// a plate.
#[test]
fn a_boss_base_gains_an_additive_toroidal_blend() {
    let mut model = ogeom_topo::Model::new();
    let plate = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 2.0), T).unwrap();
    let seat = Frame::new(
        Point::new(5.0, 5.0, 2.0),
        ogeom_math::Direction::Z,
        ogeom_math::Direction::X,
        T,
    )
    .unwrap();
    let boss = ogeom_algo::make_cylinder(&mut model, seat, 1.0, 3.0, T).unwrap();
    let joined = ogeom_bool::fuse(&mut model, &plate.shape, &boss.shape, T).unwrap();

    // The boss's base circle: radius 1 at z = 2 about (5, 5).
    let edge = explore(&model, &joined.shape, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(&model, e)
                .unwrap()
                .is_some_and(|(a, _)| {
                    model
                        .node(&a)
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                        .is_some_and(|p| {
                            (p.z - 2.0).abs() < 1e-9
                                && ((p.x - 5.0).hypot(p.y - 5.0) - 1.0).abs() < 1e-6
                        })
                })
        })
        .expect("the joined part has the boss base circle");

    let r = 0.3;
    let result = ogeom_fillet::fillet_edge(&mut model, &joined.shape, &edge, r, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let pi = core::f64::consts::PI;
    let boss_r = 1.0;
    let added = 2.0
        * pi
        * (boss_r * r * r + r * r * r / 2.0 - (boss_r + r) * pi * r * r / 4.0 + r * r * r / 3.0);
    let exact = 200.0 + pi * 3.0 + added;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - exact).abs() < 2e-3,
        "boss base fillet volume {measured} against {exact}"
    );
    assert!(result.history.is_deleted(&edge));
}
