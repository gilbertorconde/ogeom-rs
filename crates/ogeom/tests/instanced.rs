//! Solids whose topology is instanced: a prism's far cap reuses the
//! profile's nodes under the travel, and every algorithm that resolves by
//! node alone reads one name as two places.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

#[path = "support/walks.rs"]
mod walks;

use ogeom_core::Tolerances;
use ogeom_math::{Plane, Point, Vector};
use ogeom_topo::{Filter, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &ogeom_topo::Model, shape: &ogeom_topo::Shape) -> f64 {
    ogeom_algo::volume_properties(
        model,
        shape,
        ogeom_mesh::Deflection {
            chord: 1e-3,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

/// An 8 by 6 rectangle extruded 4 up: every top node is a placed reuse of a
/// bottom one.
fn rect_prism(model: &mut ogeom_topo::Model) -> ogeom_topo::Shape {
    let corners = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(8.0, 0.0, 0.0),
        Point::new(8.0, 6.0, 0.0),
        Point::new(0.0, 6.0, 0.0),
    ];
    let verts: Vec<ogeom_topo::Shape> = corners
        .iter()
        .map(|p| ogeom_algo::make_vertex(model, *p).shape)
        .collect();
    let mut edges = Vec::new();
    for i in 0..corners.len() {
        let j = (i + 1) % corners.len();
        let line = ogeom_geom::LineCurve::segment(corners[i], corners[j], T).unwrap();
        let curve = ogeom_geom::Curve::Line(line);
        let domain = ogeom_geom::Curve3d::domain(&curve);
        edges.push(
            ogeom_algo::make_edge_between(model, curve, domain, &verts[i], &verts[j], T)
                .unwrap()
                .shape,
        );
    }
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(corners[0], ogeom_math::Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face_with_pcurves(model, plane.into(), &[edges], T)
        .unwrap()
        .shape;
    ogeom_algo::make_prism(model, &face, Vector::new(0.0, 0.0, 4.0), T)
        .unwrap()
        .shape
}

fn face_where(
    model: &ogeom_topo::Model,
    solid: &ogeom_topo::Shape,
    pick: impl Fn(Point, Vector) -> bool,
) -> ogeom_topo::Shape {
    explore(model, solid, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find(|f| {
            let d = model.node(f).and_then(|n| n.data().as_face());
            let Some(d) = d else { return false };
            let Some(ogeom_geom::SurfaceGeometry::Plane(p)) = model.geometry().surface(d.surface)
            else {
                return false;
            };
            let placed = f.transform(model.datums()).unwrap();
            let origin = placed.apply(p.plane().frame().origin());
            let n = placed.apply_vector(p.plane().normal().vector());
            pick(origin, n)
        })
        .expect("the picked face")
}

#[test]
fn a_prism_bakes_with_its_far_cap_where_the_travel_put_it() {
    let mut model = ogeom_topo::Model::new();
    let solid = rect_prism(&mut model);
    let baked = ogeom_algo::baked_shape(&mut model, &solid, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &baked.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!((volume(&model, &baked.shape) - 192.0).abs() < 1e-6);
}

#[test]
fn a_prism_takes_a_general_transform_whole() {
    let mut model = ogeom_topo::Model::new();
    let solid = rect_prism(&mut model);
    let transform = ogeom_math::Transform::translation(Vector::new(1.0, 2.0, 3.0)).to_general();
    let moved = ogeom_algo::general_transformed_shape(&mut model, &solid, &transform, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &moved.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!((volume(&model, &moved.shape) - 192.0).abs() < 1e-6);
}

#[test]
fn a_prism_shells_open_at_the_top() {
    let mut model = ogeom_topo::Model::new();
    let solid = rect_prism(&mut model);
    let top = face_where(&model, &solid, |origin, n| {
        (n.z.abs() - 1.0).abs() < 1e-9 && (origin.z - 4.0).abs() < 1e-9
    });
    let shelled =
        ogeom_offset::make_thick_solid(&mut model, &solid, std::slice::from_ref(&top), 0.5, T)
            .unwrap();
    let diagnosis = ogeom_algo::check(&model, &shelled.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert_eq!(walks::edges_walked_one_way(&model, &shelled.shape), 0);
    // The walls hold half a unit each way. The open top keeps the cavity's
    // full height: 192 less 7 by 5 by 3.5.
    assert!((volume(&model, &shelled.shape) - 69.5).abs() < 1e-6);
}

#[test]
fn a_prism_wall_takes_a_draft() {
    let mut model = ogeom_topo::Model::new();
    let solid = rect_prism(&mut model);
    let wall = face_where(&model, &solid, |_, n| {
        n.z.abs() < 1e-9 && (n.y.abs() - 1.0).abs() < 1e-9
    });
    let neutral = Plane::through(Point::new(0.0, 0.0, 0.0), ogeom_math::Direction::Z);
    let angle = 2.0_f64.to_radians();
    let drafted = ogeom_offset::apply_draft(
        &mut model,
        &solid,
        std::slice::from_ref(&wall),
        neutral,
        ogeom_math::Direction::Z,
        angle,
        T,
    )
    .unwrap();
    let diagnosis = ogeom_algo::check(&model, &drafted.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert_eq!(walks::edges_walked_one_way(&model, &drafted.shape), 0);
    // One wall leans in from the neutral plane at the base: the removed
    // wedge is the wall's length times the leaning triangle.
    let expected = 8.0_f64.mul_add(-(4.0 * 4.0 / 2.0 * angle.tan()), 192.0);
    let measured = volume(&model, &drafted.shape);
    assert!(
        (measured - expected).abs() < 1e-6,
        "drafted prism volume {measured} against {expected}"
    );
}

/// The prism's wall drafted 2 degrees, then offset out and in: every face
/// a plane moved by the offset, the drafted wall's plane along its own
/// leaned normal, so the solid is the box grown by `d` on five sides with
/// the leaned wall `d` out from where it stood. The drafted wall's window
/// stands where its face does in the leaned plane's own chart; read in the
/// old plane's numbers it missed the face, and the offset's walls met the
/// leaned one beside its trim, 0.5 off their curves.
#[test]
fn a_drafted_prism_offsets_to_its_grown_planes() {
    use ogeom_geom::{Curve2d as _, Surface as _};
    let mut model = ogeom_topo::Model::new();
    let solid = rect_prism(&mut model);
    let wall = face_where(&model, &solid, |_, n| {
        n.z.abs() < 1e-9 && (n.y.abs() - 1.0).abs() < 1e-9
    });
    let neutral = Plane::through(Point::new(0.0, 0.0, 0.0), ogeom_math::Direction::Z);
    let angle = 2.0_f64.to_radians();
    let drafted = ogeom_offset::apply_draft(
        &mut model,
        &solid,
        std::slice::from_ref(&wall),
        neutral,
        ogeom_math::Direction::Z,
        angle,
        T,
    )
    .unwrap()
    .shape;
    // Every pcurve of every face inside its surface's window.
    for face in ogeom_topo::explore_unique(&model, &drafted, ShapeType::Face).unwrap() {
        let data = model.node(&face).and_then(|n| n.data().as_face()).unwrap();
        let surface = model.geometry().surface(data.surface).unwrap();
        let ((u0, u1), (v0, v1)) = surface.domain();
        for edge in ogeom_topo::explore_unique(&model, &face, ShapeType::Edge).unwrap() {
            let edge_data = model.node(&edge).and_then(|n| n.data().as_edge()).unwrap();
            let Some(ogeom_topo::EdgeRepr::PCurve { curve, range, .. }) =
                edge_data.pcurve_for(data.surface, edge.location())
            else {
                panic!("a prism edge has a pcurve on each face");
            };
            let pcurve = model.geometry().pcurve(*curve).unwrap();
            for t in [range.0, range.1] {
                let p = pcurve.point_at(t, T).unwrap();
                assert!(
                    (u0..=u1).contains(&p.x) && (v0..=v1).contains(&p.y),
                    "{p:?} outside {:?}",
                    surface.domain()
                );
            }
        }
    }
    let (cos, tan) = (angle.cos(), angle.tan());
    for d in [0.5, -0.5] {
        let grown = ogeom_offset::offset_shape(&mut model, &drafted, d, T)
            .unwrap()
            .shape;
        let diagnosis = ogeom_algo::check(&model, &grown, T).unwrap();
        assert!(diagnosis.is_valid(), "offset {d}: {:?}", diagnosis.problems);
        // x over 8 + 2d, z over 4 + 2d, y from the leaned wall (moved d
        // along its normal, d / cos further along y) to 6 + d: the mean
        // width at the mean height z = 2.
        let height = 2.0f64.mul_add(d, 4.0);
        let expected = 2.0f64.mul_add(d, 8.0) * height * (6.0 + d + d / cos - 2.0 * tan);
        let measured = volume(&model, &grown);
        assert!(
            (measured - expected).abs() < 1e-6,
            "offset {d}: {measured} against {expected}"
        );
    }
}
