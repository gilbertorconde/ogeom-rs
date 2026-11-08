//! The chamfer: a wedge subtracted through the boolean.
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

    // Four faces (both caps and the wall trimmed back, plus the bevel)
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

/// A box up to 2e7 from the origin on a tilted frame, chamfered along one
/// edge by a tenth: the wedge's faces are wound by the turn of their own
/// corners, so the box loses the wedge there as at the origin.
#[test]
fn a_chamfer_far_from_the_origin_loses_exactly_the_wedge() {
    let n = ogeom_math::Direction::new(ogeom_math::Vector::new(1.0, 1.0, 1.0), T).unwrap();
    for o in [0.0, 1e6, 1e7, 2e7] {
        let mut model = ogeom_topo::Model::new();
        let frame = Frame::new(Point::new(o, -o, o), n, ogeom_math::Direction::X, T).unwrap();
        let block = ogeom_algo::make_box(&mut model, frame, (2.0, 2.0, 2.0), T).unwrap();
        // The edge along the frame's y at local x = 2, z = 2.
        let local = |v: &ogeom_topo::Shape| {
            let p = model
                .node(v)
                .and_then(|n| n.data().as_vertex().map(|d| d.point))
                .unwrap();
            frame.to_local(p)
        };
        let edge = explore(&model, &block.shape, Filter::OfType(ShapeType::Edge))
            .unwrap()
            .into_iter()
            .find(|e| {
                ogeom_algo::edge_vertices(&model, e)
                    .unwrap()
                    .is_some_and(|(a, b)| {
                        [local(&a), local(&b)]
                            .iter()
                            .all(|p| (p.x - 2.0).abs() < 1e-6 && (p.z - 2.0).abs() < 1e-6)
                    })
            })
            .expect("the box has that edge");
        let distance = 0.1;
        let result = ogeom_fillet::chamfer_edge(&mut model, &block.shape, &edge, distance, T)
            .unwrap_or_else(|e| panic!("o {o:e}: {e}"));
        let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "o {o:e}: {:?}", diagnosis.problems);
        let props = ogeom_algo::volume_properties(
            &model,
            &result.shape,
            ogeom_mesh::Deflection::default(),
            T,
        )
        .unwrap();
        let exact = 8.0 - distance * distance / 2.0 * 2.0;
        assert!(
            (props.mass - exact).abs() < 1e-6,
            "o {o:e}: chamfer volume {} against {exact}",
            props.mass
        );
    }
}

/// A chamfer of unequal distances names the face its first distance runs
/// along, and that face must meet the edge. On a 20 x 20 x 2 plate built
/// as a prism of a square, the top face is the bottom face placed up the
/// sweep, one node in two places; a top edge with its first distance
/// named along the bottom face was taken as named along the top, and
/// bevelled there. It is refused, as on a box; named along the top it
/// loses the wedge of legs 0.5 and 1 along its length of 20.
#[test]
fn a_chamfer_named_on_a_prism_s_far_end_face_is_refused() {
    use ogeom_math::{Direction, Vector};
    let mut model = ogeom_topo::Model::new();
    let wire = ogeom_algo::make_polygon(
        &mut model,
        &[
            Point::new(0.0, 0.0, 0.0),
            Point::new(20.0, 0.0, 0.0),
            Point::new(20.0, 20.0, 0.0),
            Point::new(0.0, 20.0, 0.0),
        ],
        true,
        T,
    )
    .unwrap()
    .shape;
    let plane = ogeom_geom::PlaneSurface::over(
        ogeom_math::Plane::through(Point::ORIGIN, Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(&mut model, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let plate = ogeom_algo::make_prism(&mut model, &face, Vector::new(0.0, 0.0, 2.0), T)
        .unwrap()
        .shape;
    let at_height = |model: &ogeom_topo::Model, z: f64| {
        ogeom_topo::explore_unique(model, &plate, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                let b = ogeom_algo::shape_bounds(model, f, T).unwrap();
                (b.centre().unwrap().z - z).abs() < 1e-6 && b.size().z < 1e-3
            })
            .unwrap()
    };
    let (bottom, top) = (at_height(&model, 0.0), at_height(&model, 2.0));
    assert!(bottom.is_partner(&top), "the prism's ends share one node");
    let edge = ogeom_topo::explore_unique(&model, &plate, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let b = ogeom_algo::shape_bounds(&model, e, T).unwrap();
            let c = b.centre().unwrap();
            (c.y).abs() < 1e-6 && (c.z - 2.0).abs() < 1e-6 && b.size().x > 19.0
        })
        .unwrap();
    let spec = |face: &ogeom_topo::Shape| {
        vec![(
            edge.clone(),
            ogeom_fillet::Chamfer::Distances {
                face: face.clone(),
                on_face: 0.5,
                on_other: 1.0,
            },
        )]
    };
    let Err(ogeom_core::OgeomError::Construction(why)) =
        ogeom_fillet::chamfer_edges_with(&mut model, &plate, &spec(&bottom), T)
    else {
        panic!("a chamfer named along a face away from its edge was built");
    };
    assert!(
        why.to_string()
            .contains("the named face does not meet the edge being chamfered"),
        "{why}"
    );
    let built = ogeom_fillet::chamfer_edges_with(&mut model, &plate, &spec(&top), T).unwrap();
    let fine = ogeom_mesh::Deflection {
        chord: 1e-3,
        ..ogeom_mesh::Deflection::default()
    };
    let lost = 800.0
        - ogeom_algo::volume_properties(&model, &built.shape, fine, T)
            .unwrap()
            .mass;
    assert!((lost - 0.5 * 1.0 / 2.0 * 20.0).abs() < 1e-9, "lost {lost}");
}

/// The faces of `shape` and whether each one's own check passes.
fn faces_sound(model: &ogeom_topo::Model, shape: &ogeom_topo::Shape) {
    let diagnosis = ogeom_algo::check(model, shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    for face in ogeom_topo::explore_unique(model, shape, ShapeType::Face).unwrap() {
        let diagnosis = ogeom_algo::check(model, &face, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    }
}

/// Two bevels of 1.9 on the top and bottom edges of one 2 mm wall of a
/// 20 x 20 x 2 plate cross on that wall: the wall is cut away and the two
/// bevel planes meet in a ridge. The result is the plate less both
/// wedges, whose cross-sections, the triangles of legs 1.9 at the top and
/// bottom corners, share the triangle of base 1.8 on the wall and height
/// 0.9 to the ridge: the plate loses 20 x (2 x 1.805 - 0.81). The plate
/// built as a box or as a prism of a square answers alike, valid.
#[test]
fn bevels_crossing_on_a_thin_wall_cut_both_wedges() {
    use ogeom_math::{Direction, Vector};
    let fine = ogeom_mesh::Deflection {
        chord: 1e-3,
        ..ogeom_mesh::Deflection::default()
    };
    let plate = |model: &mut ogeom_topo::Model, prism: bool| {
        if !prism {
            return ogeom_algo::make_box(model, Frame::WORLD, (20.0, 20.0, 2.0), T)
                .unwrap()
                .shape;
        }
        let wire = ogeom_algo::make_polygon(
            model,
            &[
                Point::new(0.0, 0.0, 0.0),
                Point::new(20.0, 0.0, 0.0),
                Point::new(20.0, 20.0, 0.0),
                Point::new(0.0, 20.0, 0.0),
            ],
            true,
            T,
        )
        .unwrap()
        .shape;
        let plane = ogeom_geom::PlaneSurface::over(
            ogeom_math::Plane::through(Point::ORIGIN, Direction::Z),
            (-100.0, 100.0),
            (-100.0, 100.0),
        )
        .unwrap();
        let face = ogeom_algo::make_face(model, plane.into(), &[wire], T)
            .unwrap()
            .shape;
        ogeom_algo::make_prism(model, &face, Vector::new(0.0, 0.0, 2.0), T)
            .unwrap()
            .shape
    };
    for prism in [false, true] {
        let mut model = ogeom_topo::Model::new();
        let part = plate(&mut model, prism);
        let specs: Vec<_> = [2.0, 0.0]
            .into_iter()
            .map(|z: f64| {
                let edge = ogeom_topo::explore_unique(&model, &part, ShapeType::Edge)
                    .unwrap()
                    .into_iter()
                    .find(|e| {
                        let b = ogeom_algo::shape_bounds(&model, e, T).unwrap();
                        let c = b.centre().unwrap();
                        c.y.abs() < 1e-6 && (c.z - z).abs() < 1e-6 && b.size().x > 19.0
                    })
                    .unwrap();
                (edge, ogeom_fillet::Chamfer::Symmetric(1.9))
            })
            .collect();
        let built = ogeom_fillet::chamfer_edges_with(&mut model, &part, &specs, T).unwrap();
        faces_sound(&model, &built.shape);
        let lost = 800.0
            - ogeom_algo::volume_properties(&model, &built.shape, fine, T)
                .unwrap()
                .mass;
        let exact = 20.0 * (2.0 * 1.9 * 1.9 / 2.0 - 1.8 * 0.9 / 2.0);
        assert!(
            (lost - exact).abs() < 1e-9,
            "lost {lost} against {exact} (prism {prism})"
        );
    }
}

/// The two rims of a disc of radius 5 and thickness 2, both bevelled by
/// 1.5, cross on its wall: each wedge's cross-section is the triangle of
/// legs 1.5 at its corner, area 1.125 about radius 4.5, and the two share
/// the triangle (4.5, 1), (5, 0.5), (5, 1.5), area 0.25 about radius
/// 29/6. By Pappus the disc loses 2 pi (2 x 1.125 x 4.5 - 0.25 x 29/6),
/// and the result is valid.
#[test]
fn rim_bevels_crossing_on_a_thin_drum_wall_cut_both_rings() {
    use core::f64::consts::{PI, TAU};
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    let mut model = ogeom_topo::Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 2.0, T)
        .unwrap()
        .shape;
    let specs: Vec<_> = ogeom_topo::explore_unique(&model, &drum, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| ogeom_algo::shape_bounds(&model, e, T).unwrap().size().z < 1e-6)
        .map(|e| (e, ogeom_fillet::Chamfer::Symmetric(1.5)))
        .collect();
    assert_eq!(specs.len(), 2);
    let built = ogeom_fillet::chamfer_edges_with(&mut model, &drum, &specs, T).unwrap();
    faces_sound(&model, &built.shape);
    let lost = PI * 25.0 * 2.0
        - ogeom_algo::volume_properties(&model, &built.shape, fine, T)
            .unwrap()
            .mass;
    let exact = TAU * (2.0 * 1.125 * 4.5 - 0.25 * 29.0 / 6.0);
    assert!(
        (lost - exact).abs() < 1e-6 * exact,
        "lost {lost} against {exact}"
    );
}
