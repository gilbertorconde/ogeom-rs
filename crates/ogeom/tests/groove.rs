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

/// The ring of the through groove: radii 6..8, `y` in [2, 4], turned about
/// the y axis with its seam parked `seam_off` radians round.
fn through_ring(model: &mut ogeom_topo::Model, seam_off: f64) -> ogeom_topo::Shape {
    let (sn, cs) = seam_off.sin_cos();
    let dirp = Vector::new(cs, 0.0, -sn);
    let side = ogeom_math::Direction::new(Vector::new(sn, 0.0, cs), T).unwrap();
    let corners = [
        Point::ORIGIN + dirp * 6.0 + Vector::new(0.0, 2.0, 0.0),
        Point::ORIGIN + dirp * 8.0 + Vector::new(0.0, 2.0, 0.0),
        Point::ORIGIN + dirp * 8.0 + Vector::new(0.0, 4.0, 0.0),
        Point::ORIGIN + dirp * 6.0 + Vector::new(0.0, 4.0, 0.0),
    ];
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
            direction: ogeom_math::Direction::Y,
        },
        core::f64::consts::TAU,
        T,
    )
    .unwrap()
    .shape
}

#[test]
fn a_through_groove_cuts_out_of_both_faces_wherever_the_seam_parks() {
    // The blind groove's ring stops inside the block. This one passes clean
    // through it and out of both z faces, in four lobes, and its inner wall
    // stands at the block's own height — so the block's top plane touches it
    // along a single line. The walls are cylinders about the turn's axis, and
    // naming them as such is what lets the plane meet them in exact rulings
    // rather than in a marched curve carrying a fit error.
    for seam_off in [0.0_f64, 1.0] {
        let mut model = ogeom_topo::Model::new();
        let block = ogeom_algo::make_box(
            &mut model,
            Frame::new(
                Point::new(-10.0, -10.0, 0.0),
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
                T,
            )
            .unwrap(),
            (20.0, 20.0, 6.0),
            T,
        )
        .unwrap()
        .shape;
        let tool = through_ring(&mut model, seam_off);

        let grooved = ogeom_bool::cut(&mut model, &block, &tool, T).unwrap();

        let diagnosis = ogeom_algo::check(&model, &grooved.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

        // What the ring takes is its annulus clipped to z in [0, 6], twice its
        // 2 mm thickness: the half annulus, less what stands above z = 6 —
        // which is the outer circle's segment there, the inner circle being
        // tangent to that plane and contributing none.
        let pi = core::f64::consts::PI;
        let above = 64.0_f64.mul_add(0.75_f64.acos(), -(6.0 * 28.0_f64.sqrt()));
        let expected = 2400.0 - 2.0 * (14.0 * pi - above);
        let measured = volume(&model, &grooved.shape);
        assert!(
            (measured - expected).abs() < expected * 1e-3,
            "seam {seam_off}: grooved volume {measured} against {expected}"
        );

        // The groove is interior to the block's footprint, so a cut that only
        // ever removes material leaves the bounds exactly the block's.
        let mesh = ogeom_mesh::triangulate(
            &model,
            &grooved.shape,
            ogeom_mesh::Deflection {
                chord: 1e-3,
                ..ogeom_mesh::Deflection::default()
            },
            T,
        )
        .unwrap();
        for (axis, (low, high)) in [(-10.0, 10.0), (-10.0, 10.0), (0.0, 6.0)]
            .into_iter()
            .enumerate()
        {
            let measured_low = mesh
                .positions
                .iter()
                .map(|p| [p.x, p.y, p.z][axis])
                .fold(f64::INFINITY, f64::min);
            let measured_high = mesh
                .positions
                .iter()
                .map(|p| [p.x, p.y, p.z][axis])
                .fold(f64::NEG_INFINITY, f64::max);
            assert!(
                measured_low >= low - 1e-6 && measured_high <= high + 1e-6,
                "seam {seam_off}: axis {axis} spans {measured_low}..{measured_high}, \
                 outside the block's {low}..{high}"
            );
        }
    }
}
