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

/// The same mouth in a plate, drilled right through: the bore's rim rounded
/// over where the drill was longer than the plate is thick.
///
/// The drill's length is the whole point. A bore exactly as deep as the
/// plate leaves its mouth on the drill cylinder's own rim, and that circle
/// starts where the wall's chart does; a deeper drill leaves it as a
/// section, starting wherever the cut split it. The blend's wedge shares
/// the wall with the plate, and the shared rim's image on that chart is a
/// straight pcurve running from a quarter turn round to a quarter turn
/// past the seam, which the boolean's tear remover, seeing two samples
/// three quarters of a turn apart, read as a wrap and turned back on
/// itself. The wall then never split, its strip inside the wedge survived,
/// and the shell would not close.
#[test]
fn a_bore_mouth_in_a_plate_rounds_over_however_deep_the_drill() {
    let pi = core::f64::consts::PI;
    let (plate, bore, r) = (10.0_f64, 4.0_f64, 1.0_f64);
    for deeper in [false, true] {
        let mut model = ogeom_topo::Model::new();
        let block = ogeom_algo::make_box(
            &mut model,
            Frame::new(
                Point::new(-plate, -plate, 0.0),
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
                T,
            )
            .unwrap(),
            (plate * 2.0, plate * 2.0, 10.0),
            T,
        )
        .unwrap();
        let (low, tall) = if deeper { (-1.0, 12.0) } else { (0.0, 10.0) };
        let seat = Frame::new(
            Point::new(0.0, 0.0, low),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap();
        let drill = ogeom_algo::make_cylinder(&mut model, seat, bore, tall, T).unwrap();
        let holed = ogeom_bool::cut(&mut model, &block.shape, &drill.shape, T).unwrap();

        // The mouth, asked for away from the wall's seam so the seam's own
        // edge cannot answer instead.
        let mouth = explore(&model, &holed.shape, Filter::OfType(ShapeType::Edge))
            .unwrap()
            .into_iter()
            .find(|e| {
                use ogeom_geom::Curve3d as _;
                let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = model
                    .node(e)
                    .and_then(|n| n.data().as_edge())
                    .and_then(ogeom_topo::EdgeData::curve3d)
                else {
                    return false;
                };
                let Some(c) = model.geometry().curve(*curve) else {
                    return false;
                };
                (0..=8).all(|i| {
                    let t = range.0 + (range.1 - range.0) * f64::from(i) / 8.0;
                    c.point_at(t, T).is_ok_and(|p| {
                        (p.z - 10.0).abs() < 1e-9 && (p.x.hypot(p.y) - bore).abs() < 1e-6
                    })
                })
            })
            .expect("the plate has its bore mouth");

        let built = ogeom_fillet::fillet_edge(&mut model, &holed.shape, &mouth, r, T)
            .unwrap_or_else(|e| panic!("drill deeper {deeper}: {e}"));
        let diagnosis = ogeom_algo::check(&model, &built.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        assert_eq!(
            explore(&model, &built.shape, Filter::OfType(ShapeType::Face))
                .unwrap()
                .len(),
            8,
            "drill deeper {deeper}: six walls, the bore and the band"
        );
        // Pappus over the mirrored meridian cusp at the bore's radius, the
        // same figure the tube's rim answers to.
        let removed = 2.0
            * pi
            * (bore * r * r + r * r * r / 2.0 - (bore + r) * pi * r * r / 4.0 + r * r * r / 3.0);
        let want = plate * plate * 4.0 * 10.0 - pi * bore * bore * 10.0 - removed;
        let got = volume(&model, &built.shape);
        assert!(
            (got - want).abs() < want * 1e-5,
            "drill deeper {deeper}: {got} against {want}"
        );
    }
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
