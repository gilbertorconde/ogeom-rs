//! M4's closing bar: a mechanical part, end to end.
//!
//! One part through the whole manufacturing vocabulary: a block shelled into
//! an open tray, its corner edges chamfered, filleted at constant radius and
//! filleted at a running radius, a lofted boss and a swept nozzle fused onto
//! its floor. Every operation contributes a closed-form volume change, so
//! the final part's volume is an exact sum — the end-to-end answer is
//! checked against arithmetic, not against itself.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use og_core::Tolerances;
use og_math::{Frame, Point};
use og_topo::{Filter, Model, Shape, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

fn fine() -> og_mesh::Deflection {
    og_mesh::Deflection {
        chord: 1e-3,
        ..og_mesh::Deflection::default()
    }
}

/// The face of `solid` whose interior contains `probe`.
fn face_at(model: &Model, solid: &Shape, probe: Point) -> Shape {
    explore(model, solid, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find(|f| {
            og_algo::classify_on_face(model, f, probe, fine(), T)
                .map(|c| c == og_algo::Containment::In)
                .unwrap_or(false)
        })
        .expect("the solid has a face there")
}

/// The vertical edge of `solid` running the full height at `(x, y)`.
fn vertical_edge_at(model: &Model, solid: &Shape, x: f64, y: f64) -> Shape {
    explore(model, solid, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            og_algo::edge_vertices(model, e)
                .unwrap()
                .is_some_and(|(a, b)| {
                    let p = |v: &Shape| {
                        model
                            .node(v)
                            .and_then(|n| n.data().as_vertex().map(|d| d.point))
                            .unwrap()
                    };
                    let (pa, pb) = (p(&a), p(&b));
                    (pa.x - x).abs() < 1e-9
                        && (pa.y - y).abs() < 1e-9
                        && (pb.x - x).abs() < 1e-9
                        && (pb.y - y).abs() < 1e-9
                        && (pa.z - pb.z).abs() > 5.0
                })
        })
        .expect("the solid has that vertical edge")
}

#[test]
fn a_mechanical_part_survives_the_whole_vocabulary() {
    let mut model = Model::new();
    let pi = core::f64::consts::PI;
    let mut expected = 0.0_f64;

    // A 40 x 30 x 10 block, shelled into an open tray with 2 mm walls.
    let block = og_algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 10.0), T).unwrap();
    let top = face_at(&model, &block.shape, Point::new(20.0, 15.0, 10.0));
    let tray =
        og_offset::make_thick_solid(&mut model, &block.shape, std::slice::from_ref(&top), 2.0, T)
            .unwrap();
    expected += 40.0 * 30.0 * 10.0 - 36.0 * 26.0 * 8.0;

    // Chamfer one outer corner edge.
    let edge = vertical_edge_at(&model, &tray.shape, 0.0, 0.0);
    let chamfered = og_fillet::chamfer_edge(&mut model, &tray.shape, &edge, 1.5, T).unwrap();
    expected -= 1.5 * 1.5 / 2.0 * 10.0;

    // Fillet the next at constant radius.
    let edge = vertical_edge_at(&model, &chamfered.shape, 40.0, 0.0);
    let filleted = og_fillet::fillet_edge(&mut model, &chamfered.shape, &edge, 1.5, T).unwrap();
    expected -= 1.5 * 1.5 * (1.0 - pi / 4.0) * 10.0;

    // And a third at a radius that runs from 0.8 to 1.8 along the edge.
    let edge = vertical_edge_at(&model, &filleted.shape, 40.0, 30.0);
    let varied =
        og_fillet::fillet_edge_variable(&mut model, &filleted.shape, &edge, 0.8, 1.8, T).unwrap();
    // The law's own integral of r^2 over the unit parameter, times height.
    let integral = (1.8_f64.powi(3) - 0.8_f64.powi(3)) / (3.0 * (1.8 - 0.8));
    expected -= (1.0 - pi / 4.0) * integral * 10.0;

    // A lofted boss on the tray floor: 8 x 8 down to 4 x 4, six tall.
    let square = |model: &mut Model, half: f64, z: f64, cx: f64, cy: f64| {
        let corners = [
            Point::new(cx - half, cy - half, z),
            Point::new(cx + half, cy - half, z),
            Point::new(cx + half, cy + half, z),
            Point::new(cx - half, cy + half, z),
        ];
        og_algo::make_polygon(model, &corners, true, T)
            .unwrap()
            .shape
    };
    let bottom = square(&mut model, 4.0, 2.0, 12.0, 15.0);
    let upper = square(&mut model, 2.0, 8.0, 12.0, 15.0);
    let boss = og_offset::make_loft(&mut model, &bottom, &upper, T).unwrap();
    let with_boss = og_bool::fuse(&mut model, &varied.shape, &boss.shape, T).unwrap();
    expected += 6.0 / 3.0 * (64.0 + 16.0 + 32.0);

    // A swept nozzle standing on the floor: a straight pipe, fused on.
    let spine = {
        let line = og_geom::LineCurve::segment(
            Point::new(28.0, 15.0, 2.0),
            Point::new(28.0, 15.0, 8.0),
            T,
        )
        .unwrap();
        let curve = og_geom::Curve::Line(line);
        let domain = og_geom::Curve3d::domain(&curve);
        og_algo::make_edge(&mut model, curve, domain, T)
            .unwrap()
            .shape
    };
    let nozzle = og_offset::make_pipe(&mut model, &spine, 1.5, T).unwrap();
    let part = og_bool::fuse(&mut model, &with_boss.shape, &nozzle.shape, T).unwrap();
    expected += pi * 1.5 * 1.5 * 6.0;

    // The part holds together: closed, valid, and its volume is the sum of
    // everything done to it.
    let diagnosis = og_algo::check(&model, &part.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let measured = og_algo::volume_properties(&model, &part.shape, fine(), T)
        .unwrap()
        .mass;
    assert!(
        (measured - expected).abs() < 0.05,
        "part volume {measured} against the operations' own sum {expected}"
    );

    // The cavity is open above and hollow within.
    let hollow =
        og_algo::classify_in_solid_exact(&model, &part.shape, Point::new(30.0, 25.0, 5.0), T)
            .unwrap();
    assert_eq!(hollow, og_algo::Containment::Out);
    // The boss and nozzle are material.
    for probe in [Point::new(12.0, 15.0, 5.0), Point::new(28.0, 15.0, 5.0)] {
        let inside = og_algo::classify_in_solid_exact(&model, &part.shape, probe, T).unwrap();
        assert_eq!(inside, og_algo::Containment::In, "at {probe:?}");
    }
    // And the walls are walls.
    let wall = og_algo::classify_in_solid_exact(&model, &part.shape, Point::new(1.0, 15.0, 5.0), T)
        .unwrap();
    assert_eq!(wall, og_algo::Containment::In);
}
