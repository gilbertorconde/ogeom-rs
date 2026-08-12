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

/// The top rim of a cylinder, found by a vertex at the rim's height and
/// radius.
fn rim_edge(
    model: &ogeom_topo::Model,
    solid: &ogeom_topo::Shape,
    height: f64,
    radius: f64,
) -> ogeom_topo::Shape {
    explore(model, solid, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(model, e)
                .unwrap()
                .is_some_and(|(a, _)| {
                    model
                        .node(&a)
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                        .is_some_and(|p| {
                            (p.z - height).abs() < 1e-9 && (p.x.hypot(p.y) - radius).abs() < 1e-6
                        })
                })
        })
        .expect("the solid has the rim")
}

#[test]
fn a_cylinder_rim_chamfers_to_a_cone_with_the_distance_angle_form() {
    let mut model = ogeom_topo::Model::new();
    let (radius, height) = (5.0, 10.0);
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, radius, height, T).unwrap();
    let edge = rim_edge(&model, &drum.shape, height, radius);

    // The named face is the wall; at forty-five degrees the derived distance
    // along the cap equals the axial one.
    let wall = explore(&model, &drum.shape, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find(|f| {
            model
                .node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Cylinder(_)))
        })
        .expect("the drum has its wall");
    let (distance, angle) = (1.0, core::f64::consts::FRAC_PI_4);
    let result =
        ogeom_fillet::chamfer_edge_angle(&mut model, &drum.shape, &edge, &wall, distance, angle, T)
            .unwrap();

    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // The removed ring is a cylinder slab minus the cone frustum it leaves:
    // slab pi R^2 w, frustum pi w/3 (R^2 + R(R-c) + (R-c)^2), legs w = c = 1.
    let pi = core::f64::consts::PI;
    let (w, c) = (distance, distance * angle.tan());
    let inner = radius - c;
    let removed = pi * radius * radius * w
        - pi * w / 3.0 * (radius * radius + radius * inner + inner * inner);
    let exact = pi * radius * radius * height - removed;
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    let measured = ogeom_algo::volume_properties(&model, &result.shape, fine, T)
        .unwrap()
        .mass;
    // The budget is the mesh's, not the chamfer's: inscribed chords on a
    // radius-five drum undercut the true volume by a few hundredths at this
    // deflection.
    assert!(
        (measured - exact).abs() < 5e-2,
        "rim chamfer volume {measured} against {exact}"
    );

    // Four faces — both caps and the wall trimmed back, plus the bevel —
    // and the bevel is a genuine cone.
    let faces = explore(&model, &result.shape, Filter::OfType(ShapeType::Face)).unwrap();
    assert_eq!(faces.len(), 4);
    let cones = faces
        .iter()
        .filter(|f| {
            model
                .node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Cone(_)))
        })
        .count();
    assert_eq!(cones, 1, "the bevel face is a cone");
    assert!(result.history.is_deleted(&edge));

    // At forty-five degrees the distance-angle form is the symmetric one.
    let mut second = ogeom_topo::Model::new();
    let again = ogeom_algo::make_cylinder(&mut second, Frame::WORLD, radius, height, T).unwrap();
    let rim = rim_edge(&second, &again.shape, height, radius);
    let symmetric =
        ogeom_fillet::chamfer_edge(&mut second, &again.shape, &rim, distance, T).unwrap();
    let symmetric_volume = ogeom_algo::volume_properties(&second, &symmetric.shape, fine, T)
        .unwrap()
        .mass;
    assert!(
        (symmetric_volume - measured).abs() < 1e-9,
        "the two spellings should build the same solid"
    );
}

#[test]
fn a_boss_base_chamfer_adds_the_bevel_ring() {
    // The concave revolved seat: the wedge sits in the open corner where the
    // boss meets the plate, and the boolean fuses it.
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

    let d = 0.3;
    let result = ogeom_fillet::chamfer_edge(&mut model, &joined.shape, &edge, d, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // Pappus over the corner triangle: legs d up the boss and d out along
    // the plate, centroid a third of the way out.
    let pi = core::f64::consts::PI;
    let boss_r = 1.0;
    let added = 2.0 * pi * (boss_r + d / 3.0) * (d * d / 2.0);
    let exact = 200.0 + pi * 3.0 + added;
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    let measured = ogeom_algo::volume_properties(&model, &result.shape, fine, T)
        .unwrap()
        .mass;
    assert!(
        (measured - exact).abs() < 2e-3,
        "boss base chamfer volume {measured} against {exact}"
    );
    assert!(result.history.is_deleted(&edge));
}

#[test]
fn a_rim_chamfer_that_swallows_the_axis_is_refused() {
    let mut model = ogeom_topo::Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 1.0, 2.0, T).unwrap();
    let edge = rim_edge(&model, &drum.shape, 2.0, 1.0);
    assert!(
        ogeom_fillet::chamfer_edge(&mut model, &drum.shape, &edge, 1.5, T).is_err(),
        "a cap leg past the axis has nothing to stand on"
    );
}
