//! Booleans whose sections only the marching machinery can trace.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Frame, Point};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &ogeom_topo::Model, shape: &ogeom_topo::Shape, chord: f64) -> f64 {
    ogeom_algo::volume_properties(
        model,
        shape,
        ogeom_mesh::Deflection {
            chord,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

#[test]
fn a_branch_standing_on_a_drum_fuses_through_the_marched_seam() {
    // The seam winds the branch's whole chart and is metres of marching:
    // the walk has to close within its budget, the fitted section has to
    // meet the tolerance it states, and the branch wall has to split along
    // it.
    let mut model = ogeom_topo::Model::new();
    let main = Frame::new(
        Point::new(-20.0, 0.0, 0.0),
        ogeom_math::Direction::X,
        ogeom_math::Direction::Y,
        T,
    )
    .unwrap();
    let drum = ogeom_algo::make_cylinder(&mut model, main, 10.0, 40.0, T).unwrap();
    let branch = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T).unwrap();
    let joined = ogeom_bool::fuse(&mut model, &drum.shape, &branch.shape, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &joined.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // The drum plus the branch's protruding stub: the stub is the branch
    // above the saddle, and the saddle dips below the drum's crown by at
    // most the sagitta of a five-wide chord on a ten-radius circle.
    let pi = core::f64::consts::PI;
    let drum_volume = pi * 100.0 * 40.0;
    let stub_low = pi * 25.0 * 10.0;
    let stub_high = pi * 25.0 * (10.0 + 1.34);
    let measured = volume(&model, &joined.shape, 1e-3);
    assert!(
        measured > drum_volume + stub_low && measured < drum_volume + stub_high,
        "fused volume {measured} outside [{}, {}]",
        drum_volume + stub_low,
        drum_volume + stub_high
    );
}

#[test]
fn two_equal_drums_cross_in_their_two_exact_ellipses() {
    // Equal radii are the one crossing whose quartic factors: two ellipses
    // in the axes' bisector planes, spoken exactly, so the fuse holds
    // together through the tangent poles no march can walk. The volume is
    // the classic identity: two drums less their bicylinder.
    let mut model = ogeom_topo::Model::new();
    let upright = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T).unwrap();
    let across_frame = Frame::new(
        Point::new(-20.0, 0.0, 10.0),
        ogeom_math::Direction::X,
        ogeom_math::Direction::Y,
        T,
    )
    .unwrap();
    let across = ogeom_algo::make_cylinder(&mut model, across_frame, 5.0, 40.0, T).unwrap();
    let joined = ogeom_bool::fuse(&mut model, &upright.shape, &across.shape, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &joined.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let ellipses = ogeom_topo::explore_unique(&model, &joined.shape, ogeom_topo::ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            model
                .node(e)
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.curve3d())
                .and_then(|r| {
                    let ogeom_topo::EdgeRepr::Curve3d { curve, .. } = r else {
                        return None;
                    };
                    model.geometry().curve(*curve)
                })
                .is_some_and(|c| matches!(c, ogeom_geom::Curve::Ellipse(_)))
        })
        .count();
    assert!(
        ellipses >= 4,
        "the seam rebuilt as exact ellipse arcs: {ellipses}"
    );

    let pi = core::f64::consts::PI;
    let bicylinder = 16.0 * 125.0 / 3.0;
    let expected = pi * 25.0 * 20.0 + pi * 25.0 * 40.0 - bicylinder;
    let measured = volume(&model, &joined.shape, 1e-3);
    assert!(
        (measured - expected).abs() < expected * 1e-3,
        "Steinmetz union {measured} against {expected}"
    );
}

/// A torus elbow fused onto a drum, cut at a mitre plane, and a second
/// drum cut at the same plane, fused face to face. Each mitre cap is
/// bounded by an exact ellipse on one side and the torus's fitted section
/// on the other, and the two leave each other tangentially, parting into
/// lunes from about a micron to a millimetre wide, the fitted one crossing
/// the ellipse at a shallow turn in some directions. The two pieces only
/// touch, so the fuse holds their volumes' sum. In all but one direction
/// the fuse refused.
#[test]
fn mitred_pieces_whose_caps_part_tangentially_fuse() {
    use ogeom_geom::{CircleCurve, PlaneSurface};
    use ogeom_math::{Axis, Circle, Direction, Plane, Vector};
    let disc = |model: &mut ogeom_topo::Model, frame: Frame, r: f64| {
        let circle = Circle::new(frame, r, T).unwrap();
        let edge = ogeom_algo::make_edge(
            model,
            CircleCurve::new(circle).into(),
            (0.0, 2.0 * core::f64::consts::PI),
            T,
        )
        .unwrap()
        .shape;
        let wire = ogeom_algo::make_wire(model, &[edge], T).unwrap().shape;
        ogeom_algo::make_face(
            model,
            PlaneSurface::new(Plane::new(frame)).into(),
            &[wire],
            T,
        )
        .unwrap()
        .shape
    };
    let (big, r, half) = (10.0_f64, 1.0_f64, 0.5_f64);
    for (sx, sz) in [
        (0.0, 0.6),
        (0.0, 0.3),
        (0.0, -0.6),
        (0.6, 0.0),
        (-0.6, 0.0),
        (0.4, 0.4),
        (0.0, 0.05),
    ] {
        let mut model = ogeom_topo::Model::new();
        let corner = Point::new(big, 0.0, 0.0);
        let c0 = Point::new(big * half.cos(), -big * half.sin(), 0.0);
        let tangent = Vector::new(half.sin(), half.cos(), 0.0);
        let radial = Vector::new(half.cos(), -half.sin(), 0.0);
        let f0 = Frame::new(
            c0,
            Direction::new(tangent, T).unwrap(),
            Direction::new(radial, T).unwrap(),
            T,
        )
        .unwrap();
        let d0 = disc(&mut model, f0, r);
        let axis = Axis {
            location: Point::ORIGIN,
            direction: Direction::Z,
        };
        let arc = ogeom_algo::make_revolution(&mut model, &d0, axis, half, T)
            .unwrap()
            .shape;
        let fy = Frame::new(corner, Direction::Y, Direction::X, T).unwrap();
        let on = ogeom_algo::make_cylinder(&mut model, fy, r, 3.0, T)
            .unwrap()
            .shape;
        let torus = ogeom_bool::fuse(&mut model, &arc, &on, T).unwrap().shape;
        let d = Vector::new(sx, 1.0, sz).normalized(T).unwrap();
        let side = Vector::new(0.0, 1.0, 0.0).cross(d).normalized(T).unwrap();
        let fc = Frame::new(
            corner - d * 3.0,
            Direction::new(d, T).unwrap(),
            Direction::new(side, T).unwrap(),
            T,
        )
        .unwrap();
        let drum = ogeom_algo::make_cylinder(&mut model, fc, r, 6.0, T)
            .unwrap()
            .shape;
        let n = (Vector::new(0.0, 1.0, 0.0) + d).normalized(T).unwrap();
        let pf = Frame::new(
            corner,
            Direction::new(n, T).unwrap(),
            Direction::new(side, T).unwrap(),
            T,
        )
        .unwrap();
        let plane = PlaneSurface::over(Plane::new(pf), (-50.0, 50.0), (-50.0, 50.0)).unwrap();
        let corners: Vec<Point> = [(-20.0, -20.0), (20.0, -20.0), (20.0, 20.0), (-20.0, 20.0)]
            .iter()
            .map(|&(a, b)| pf.to_world(Point::new(a, b, 0.0)))
            .collect();
        let wire = ogeom_algo::make_polygon(&mut model, &corners, true, T)
            .unwrap()
            .shape;
        let face = ogeom_algo::make_face(&mut model, plane.into(), &[wire], T)
            .unwrap()
            .shape;
        let after = ogeom_algo::make_half_space(&mut model, &face, corner + n * 5.0, T)
            .unwrap()
            .shape;
        let before = ogeom_bool::cut(&mut model, &torus, &after, T)
            .unwrap()
            .shape;
        let beyond = ogeom_bool::common(&mut model, &drum, &after, T)
            .unwrap()
            .shape;
        let pieces = volume(&model, &before, 1e-3) + volume(&model, &beyond, 1e-3);
        let fused = ogeom_bool::fuse(&mut model, &before, &beyond, T)
            .unwrap_or_else(|e| panic!("direction ({sx}, {sz}): {e}"))
            .shape;
        assert!(
            ogeom_algo::check(&model, &fused, T).unwrap().is_valid(),
            "direction ({sx}, {sz})"
        );
        let got = volume(&model, &fused, 1e-3);
        assert!(
            (got - pieces).abs() < 1e-8,
            "direction ({sx}, {sz}): {got} against the pieces' {pieces}"
        );
    }
}
