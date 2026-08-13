//! The blind ring groove: a full-turn revolved tool whose own top face lies
//! in the face it cuts into.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Frame, Plane, Point, Vector};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &ogeom_topo::Model, s: &ogeom_topo::Shape) -> f64 {
    ogeom_algo::volume_properties(
        model,
        s,
        ogeom_mesh::Deflection {
            chord: 1e-3,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

/// A square ring profile in the meridian half-plane `seam_off` radians round,
/// wound `flip`-wise.
fn ring(model: &mut ogeom_topo::Model, seam_off: f64, flip: bool) -> ogeom_topo::Shape {
    let (sn, cs) = seam_off.sin_cos();
    let dirp = Vector::new(cs, sn, 0.0);
    let side = ogeom_math::Direction::new(Vector::new(-sn, cs, 0.0), T).unwrap();
    let mut corners = [
        Point::ORIGIN + dirp * 7.0 + Vector::new(0.0, 0.0, 4.0),
        Point::ORIGIN + dirp * 9.0 + Vector::new(0.0, 0.0, 4.0),
        Point::ORIGIN + dirp * 9.0 + Vector::new(0.0, 0.0, 6.0),
        Point::ORIGIN + dirp * 7.0 + Vector::new(0.0, 0.0, 6.0),
    ];
    if flip {
        corners.reverse();
    }
    let verts: Vec<ogeom_topo::Shape> = corners
        .iter()
        .map(|p| ogeom_algo::make_vertex(model, *p).shape)
        .collect();
    let mut edges = Vec::new();
    for i in 0..4 {
        let j = (i + 1) % 4;
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
        Plane::through(corners[0], side),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face_with_pcurves(model, plane.into(), &[edges], T)
        .unwrap()
        .shape;
    ogeom_algo::make_revolution(
        model,
        &face,
        ogeom_math::Axis {
            location: Point::ORIGIN,
            direction: ogeom_math::Direction::Z,
        },
        core::f64::consts::TAU,
        T,
    )
    .unwrap()
    .shape
}

#[test]
fn a_revolution_stands_right_side_out_whichever_way_its_profile_winds() {
    // The material side follows the wire's own walk, measured from the
    // traversal itself — both windings of one profile revolve into the same
    // ring, right side out.
    for flip in [false, true] {
        let mut model = ogeom_topo::Model::new();
        let solid = ring(&mut model, 0.0, flip);
        let diagnosis = ogeom_algo::check(&model, &solid, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        let expected = core::f64::consts::PI * (81.0 - 49.0) * 2.0;
        let measured = volume(&model, &solid);
        assert!(
            (measured - expected).abs() < expected * 1e-3,
            "flip {flip}: ring volume {measured} against {expected}"
        );
    }
}

#[test]
fn a_blind_groove_cuts_into_the_top_wherever_the_seam_parks() {
    // The ring's own top face lies in the block's top plane: the coplanar
    // pair melts, the groove opens, and the seam's meridian may sit
    // anywhere.
    for seam_off in [0.0_f64, 1.0] {
        let mut model = ogeom_topo::Model::new();
        let block = ogeom_algo::make_box(
            &mut model,
            Frame::new(
                Point::new(-15.0, -15.0, 0.0),
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
                T,
            )
            .unwrap(),
            (30.0, 30.0, 6.0),
            T,
        )
        .unwrap()
        .shape;
        let tool = ring(&mut model, seam_off, false);
        let grooved = ogeom_bool::cut(&mut model, &block, &tool, T).unwrap();
        let diagnosis = ogeom_algo::check(&model, &grooved.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        let expected = 30.0_f64.mul_add(30.0 * 6.0, -(core::f64::consts::PI * 32.0 * 2.0));
        let measured = volume(&model, &grooved.shape);
        assert!(
            (measured - expected).abs() < expected * 1e-3,
            "seam {seam_off}: grooved volume {measured} against {expected}"
        );
    }
}
