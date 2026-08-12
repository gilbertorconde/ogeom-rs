//! The solid offset and the shells cut from it.
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

#[test]
fn a_box_offset_outward_is_the_bigger_box() {
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let result = ogeom_offset::offset_shape(&mut model, &block.shape, 0.5, T).unwrap();

    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - 27.0).abs() < 1e-9,
        "outward box offset volume {measured} against 27"
    );
    // Topology preserved one-for-one.
    let faces = explore(&model, &result.shape, Filter::OfType(ShapeType::Face)).unwrap();
    assert_eq!(faces.len(), 6);
}

#[test]
fn a_box_offset_inward_is_the_smaller_box() {
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let result = ogeom_offset::offset_shape(&mut model, &block.shape, -0.5, T).unwrap();
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - 1.0).abs() < 1e-9,
        "inward box offset volume {measured} against 1"
    );
}

#[test]
fn a_cylinder_offset_grows_radius_and_caps_alike() {
    let mut model = ogeom_topo::Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 1.0, 2.0, T).unwrap();
    let result = ogeom_offset::offset_shape(&mut model, &drum.shape, 0.5, T).unwrap();

    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let expected = core::f64::consts::PI * 1.5 * 1.5 * 3.0;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 5e-3,
        "cylinder offset volume {measured} against {expected}"
    );
}

#[test]
fn an_offset_that_swallows_the_box_is_refused() {
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    assert!(ogeom_offset::offset_shape(&mut model, &block.shape, -1.5, T).is_err());
}

/// The face of `solid` whose interior contains `probe`.
fn face_at(
    model: &ogeom_topo::Model,
    solid: &ogeom_topo::Shape,
    probe: Point,
) -> ogeom_topo::Shape {
    explore(model, solid, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find(|f| {
            ogeom_algo::classify_on_face(model, f, probe, fine(), T)
                .map(|c| c == ogeom_algo::Containment::In)
                .unwrap_or(false)
        })
        .expect("the solid has a face there")
}

#[test]
fn a_shelled_box_keeps_its_walls_and_opens_its_top() {
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let top = face_at(&model, &block.shape, Point::new(1.0, 1.0, 2.0));

    let t = 0.2;
    let result =
        ogeom_offset::make_thick_solid(&mut model, &block.shape, std::slice::from_ref(&top), t, T)
            .unwrap();

    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let expected = 8.0 - 1.6 * 1.6 * 1.8;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 1e-9,
        "box shell volume {measured} against {expected}"
    );

    // Five outer walls, the top ring, and five cavity walls.
    let faces = explore(&model, &result.shape, Filter::OfType(ShapeType::Face)).unwrap();
    assert_eq!(faces.len(), 11);
    assert!(result.history.is_deleted(&top));

    // The cavity is real: its centre is outside the material.
    let inside =
        ogeom_algo::classify_in_solid_exact(&model, &result.shape, Point::new(1.0, 1.0, 1.0), T)
            .unwrap();
    assert_eq!(inside, ogeom_algo::Containment::Out);
    // And the wall is material.
    let wall =
        ogeom_algo::classify_in_solid_exact(&model, &result.shape, Point::new(0.1, 1.0, 1.0), T)
            .unwrap();
    assert_eq!(wall, ogeom_algo::Containment::In);
}

#[test]
fn a_shelled_cylinder_becomes_a_cup() {
    let mut model = ogeom_topo::Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 1.0, 2.0, T).unwrap();
    let top = face_at(&model, &drum.shape, Point::new(0.0, 0.0, 2.0));

    let t = 0.2;
    let result =
        ogeom_offset::make_thick_solid(&mut model, &drum.shape, std::slice::from_ref(&top), t, T)
            .unwrap();

    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let pi = core::f64::consts::PI;
    let expected = pi * 2.0 - pi * 0.8 * 0.8 * 1.8;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 1e-3,
        "cup volume {measured} against {expected}"
    );
    assert!(result.history.is_deleted(&top));
}

/// A box with its top edge filleted: the part whose offset meets a partial
/// cylinder, arcs with vertices, and vertices seated on curved faces.
fn filleted_box(model: &mut ogeom_topo::Model) -> ogeom_topo::Shape {
    let block = ogeom_algo::make_box(model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let edge = explore(model, &block.shape, Filter::OfType(ShapeType::Edge))
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
        .expect("the box has that edge");
    ogeom_fillet::fillet_edge(model, &block.shape, &edge, 0.5, T)
        .unwrap()
        .shape
}

#[test]
fn a_filleted_box_offsets_with_its_blend() {
    let pi = core::f64::consts::PI;
    for (w, side, blend) in [(0.2, 2.4, 0.7), (-0.2, 1.6, 0.3)] {
        let mut model = ogeom_topo::Model::new();
        let part = filleted_box(&mut model);
        let result = ogeom_offset::offset_shape(&mut model, &part, w, T).unwrap();
        let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        let expected = side * side * side - (1.0 - pi / 4.0) * blend * blend * side;
        let measured = volume(&model, &result.shape);
        assert!(
            (measured - expected).abs() < 2e-3,
            "offset {w}: volume {measured} against {expected}"
        );
        // Topology preserved: still seven faces.
        assert_eq!(
            explore(&model, &result.shape, Filter::OfType(ShapeType::Face))
                .unwrap()
                .len(),
            7
        );
    }
}

#[test]
fn a_filleted_box_shells_blend_and_all() {
    let mut model = ogeom_topo::Model::new();
    let part = filleted_box(&mut model);
    let bottom = face_at(&model, &part, Point::new(1.0, 1.0, 0.0));
    let t = 0.2;
    let result =
        ogeom_offset::make_thick_solid(&mut model, &part, std::slice::from_ref(&bottom), t, T)
            .unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let pi = core::f64::consts::PI;
    let outer = 8.0 - (1.0 - pi / 4.0) * 0.25 * 2.0;
    let cavity = 1.6 * 1.6 * 1.8 - (1.0 - pi / 4.0) * 0.09 * 1.6;
    let expected = outer - cavity;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 2e-3,
        "filleted shell volume {measured} against {expected}"
    );
    assert!(result.history.is_deleted(&bottom));
}

#[test]
fn a_cone_offsets_to_the_parallel_cone_apex_and_all() {
    // Base radius 3, height 4: the 3-4-5 cone. Offset by 0.5, the lateral
    // surface becomes the parallel cone, the base drops by the amount, and
    // the apex — a vertex with no normal to offer — lands on the parallel
    // cone's own apex, 0.5 / sin over the old one.
    let mut model = ogeom_topo::Model::new();
    let cone = ogeom_algo::make_cone(&mut model, Frame::WORLD, 3.0, 0.0, 4.0, T).unwrap();
    let result = ogeom_offset::offset_shape(&mut model, &cone.shape, 0.5, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // Radius 4 at the moved base, apex height 32/6 above it.
    let pi = core::f64::consts::PI;
    let expected = pi / 3.0 * 16.0 * (16.0 / 3.0);
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 1e-2,
        "offset cone volume {measured} against {expected}"
    );
    assert_eq!(
        explore(&model, &result.shape, Filter::OfType(ShapeType::Face))
            .unwrap()
            .len(),
        2
    );
    // The apex vertex itself: at 29/6 above the original base plane.
    let apex = explore(&model, &result.shape, Filter::OfType(ShapeType::Vertex))
        .unwrap()
        .into_iter()
        .any(|v| {
            model
                .node(&v)
                .and_then(|n| n.data().as_vertex().map(|d| d.point))
                .is_some_and(|p| (p.z - 29.0 / 6.0).abs() < 1e-6 && p.x.hypot(p.y) < 1e-6)
        });
    assert!(apex, "the apex rode the parallel cone");
}

#[test]
fn a_rim_filleted_cylinder_offsets_with_its_toroidal_band() {
    // The fillet's torus keeps its centre circle and grows its tube, which
    // preserves the tangencies — so the offset part is the rim-filleted
    // version of the offset cylinder, and the same closed form prices it.
    let mut model = ogeom_topo::Model::new();
    let (radius, height, blend, w) = (2.0, 4.0, 0.5, 0.3);
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, radius, height, T).unwrap();
    let rim = explore(&model, &drum.shape, Filter::OfType(ShapeType::Edge))
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
                            (p.z - height).abs() < 1e-9 && (p.x.hypot(p.y) - radius).abs() < 1e-6
                        })
                })
        })
        .expect("the drum has its top rim");
    let part = ogeom_fillet::fillet_edge(&mut model, &drum.shape, &rim, blend, T)
        .unwrap()
        .shape;

    let result = ogeom_offset::offset_shape(&mut model, &part, w, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let pi = core::f64::consts::PI;
    let (r, h, b) = (radius + w, height + 2.0 * w, blend + w);
    let removed =
        2.0 * pi * (r * b * b - b * b * b / 2.0 - (r - b) * pi * b * b / 4.0 - b * b * b / 3.0);
    let expected = pi * r * r * h - removed;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 5e-3,
        "offset rim-filleted drum volume {measured} against {expected}"
    );
    // Topology preserved: both caps, the wall, and the band.
    assert_eq!(
        explore(&model, &result.shape, Filter::OfType(ShapeType::Face))
            .unwrap()
            .len(),
        4
    );
}
