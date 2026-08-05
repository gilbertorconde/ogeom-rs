//! Pipes and lofts, against Pappus and the frustum formulae.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::Curve3d as _;
use ogeom_math::{Circle, Frame, Point};
use ogeom_topo::Filter;

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
fn a_straight_pipe_is_a_cylinder() {
    let mut model = ogeom_topo::Model::new();
    let line = ogeom_geom::LineCurve::segment(Point::ORIGIN, Point::new(0.0, 0.0, 2.0), T).unwrap();
    let curve = ogeom_geom::Curve::Line(line);
    let domain = curve.domain();
    let spine = ogeom_algo::make_edge(&mut model, curve, domain, T)
        .unwrap()
        .shape;

    let result = ogeom_offset::make_pipe(&mut model, &spine, 0.3, T).unwrap();
    let expected = core::f64::consts::PI * 0.09 * 2.0;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 5e-4,
        "straight pipe volume {measured} against {expected}"
    );
    assert!(!result.history.generated(&spine).is_empty());
}

#[test]
fn a_quarter_arc_pipe_is_a_torus_segment() {
    let mut model = ogeom_topo::Model::new();
    let circle = Circle::new(Frame::WORLD, 2.0, T).unwrap();
    let curve = ogeom_geom::Curve::Circle(ogeom_geom::CircleCurve::new(circle));
    let spine = ogeom_algo::make_edge(&mut model, curve, (0.0, core::f64::consts::FRAC_PI_2), T)
        .unwrap()
        .shape;

    let result = ogeom_offset::make_pipe(&mut model, &spine, 0.3, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // Pappus: the tube's area rides the spine's length.
    let expected = core::f64::consts::PI * 0.09 * 2.0 * core::f64::consts::FRAC_PI_2;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 1e-3,
        "arc pipe volume {measured} against {expected}"
    );

    // Four faces: two half tubes and two meridian caps.
    let faces = ogeom_topo::explore(
        &model,
        &result.shape,
        Filter::OfType(ogeom_topo::ShapeType::Face),
    )
    .unwrap();
    assert_eq!(faces.len(), 4);
}

#[test]
fn a_closed_circular_pipe_is_the_whole_torus() {
    let mut model = ogeom_topo::Model::new();
    let circle = Circle::new(Frame::WORLD, 2.0, T).unwrap();
    let curve = ogeom_geom::Curve::Circle(ogeom_geom::CircleCurve::new(circle));
    let domain = curve.domain();
    let spine = ogeom_algo::make_edge(&mut model, curve, domain, T)
        .unwrap()
        .shape;

    let result = ogeom_offset::make_pipe(&mut model, &spine, 0.3, T).unwrap();
    let expected = 2.0 * core::f64::consts::PI * core::f64::consts::PI * 2.0 * 0.09;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 2e-3,
        "torus pipe volume {measured} against {expected}"
    );
}

#[test]
fn a_pipe_that_swallows_its_spine_is_refused() {
    let mut model = ogeom_topo::Model::new();
    let circle = Circle::new(Frame::WORLD, 0.5, T).unwrap();
    let curve = ogeom_geom::Curve::Circle(ogeom_geom::CircleCurve::new(circle));
    let domain = curve.domain();
    let spine = ogeom_algo::make_edge(&mut model, curve, domain, T)
        .unwrap()
        .shape;
    assert!(ogeom_offset::make_pipe(&mut model, &spine, 0.6, T).is_err());
}

fn square(model: &mut ogeom_topo::Model, half: f64, z: f64) -> ogeom_topo::Shape {
    let corners = [
        Point::new(-half, -half, z),
        Point::new(half, -half, z),
        Point::new(half, half, z),
        Point::new(-half, half, z),
    ];
    ogeom_algo::make_polygon(model, &corners, true, T)
        .unwrap()
        .shape
}

#[test]
fn a_polygonal_loft_is_the_frustum_pyramid() {
    let mut model = ogeom_topo::Model::new();
    let bottom = square(&mut model, 1.0, 0.0);
    let top = square(&mut model, 0.5, 2.0);

    let result = ogeom_offset::make_loft(&mut model, &bottom, &top, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // The pyramidal frustum: h/3 (A1 + A2 + sqrt(A1 A2)).
    let expected = 2.0 / 3.0 * (4.0 + 1.0 + 2.0);
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 1e-9,
        "polygonal loft volume {measured} against {expected}"
    );
    assert!(!result.history.generated(&bottom).is_empty());
}

#[test]
fn a_circular_loft_is_the_cone_frustum() {
    let mut model = ogeom_topo::Model::new();
    let ring = |model: &mut ogeom_topo::Model, r: f64, z: f64| {
        let frame = Frame::new(
            Point::new(0.0, 0.0, z),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap();
        let circle = Circle::new(frame, r, T).unwrap();
        let curve = ogeom_geom::Curve::Circle(ogeom_geom::CircleCurve::new(circle));
        let domain = curve.domain();
        let edge = ogeom_algo::make_edge(model, curve, domain, T)
            .unwrap()
            .shape;
        ogeom_algo::make_wire(model, std::slice::from_ref(&edge), T)
            .unwrap()
            .shape
    };
    let bottom = ring(&mut model, 1.0, 0.0);
    let top = ring(&mut model, 0.5, 2.0);

    let result = ogeom_offset::make_loft(&mut model, &bottom, &top, T).unwrap();
    let pi = core::f64::consts::PI;
    let expected = pi * 2.0 / 3.0 * 0.5_f64.mul_add(0.5, 1.0f64.mul_add(1.0, 1.0 * 0.5));
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 1e-3,
        "circular loft volume {measured} against {expected}"
    );
}

#[test]
fn a_twisted_loft_is_refused_as_skew() {
    let mut model = ogeom_topo::Model::new();
    let bottom = square(&mut model, 1.0, 0.0);
    // The top square rotated 45 degrees: every wall would be skew.
    let corners = [
        Point::new(0.0, -0.7, 2.0),
        Point::new(0.7, 0.0, 2.0),
        Point::new(0.0, 0.7, 2.0),
        Point::new(-0.7, 0.0, 2.0),
    ];
    let top = ogeom_algo::make_polygon(&mut model, &corners, true, T)
        .unwrap()
        .shape;
    assert!(ogeom_offset::make_loft(&mut model, &bottom, &top, T).is_err());
}

#[test]
fn a_skinned_loft_through_cone_sections_measures_as_the_frustum() {
    let mut model = ogeom_topo::Model::new();
    let ring = |model: &mut ogeom_topo::Model, r: f64, z: f64| {
        let frame = Frame::new(
            Point::new(0.0, 0.0, z),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap();
        let circle = Circle::new(frame, r, T).unwrap();
        let curve = ogeom_geom::Curve::Circle(ogeom_geom::CircleCurve::new(circle));
        let domain = curve.domain();
        let edge = ogeom_algo::make_edge(model, curve, domain, T)
            .unwrap()
            .shape;
        ogeom_algo::make_wire(model, std::slice::from_ref(&edge), T)
            .unwrap()
            .shape
    };
    let sections = [
        ring(&mut model, 1.0, 0.0),
        ring(&mut model, 0.75, 1.0),
        ring(&mut model, 0.5, 2.0),
    ];
    let result = ogeom_offset::make_loft_skinned(&mut model, &sections, 1e-3, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // Linear radii through cone sections: the skin reproduces the frustum.
    let pi = core::f64::consts::PI;
    let expected = pi * 2.0 / 3.0 * (1.0 + 0.5 + 0.25);
    let measured = volume(&model, &result.shape);
    // A skin at its own stated error: the volume deficit is the fitted
    // sections riding just inside their circles.
    assert!(
        (measured - expected).abs() < 1e-2,
        "skinned frustum volume {measured} against {expected}"
    );
}

#[test]
fn a_pipe_along_a_free_form_spine_holds_pappus() {
    let mut model = ogeom_topo::Model::new();
    // A gentle S in the xz plane.
    let spine_curve = ogeom_geom::Curve::BSpline(
        ogeom_geom::BSplineCurve::new(
            ogeom_math::KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 3).unwrap(),
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(2.0, 0.0, 1.0),
                Point::new(4.0, 0.0, -1.0),
                Point::new(6.0, 0.0, 0.0),
            ],
            T,
        )
        .unwrap(),
    );
    let domain = ogeom_geom::Curve3d::domain(&spine_curve);
    let spine = ogeom_algo::make_edge(&mut model, spine_curve.clone(), domain, T)
        .unwrap()
        .shape;
    let r = 0.2;
    let result = ogeom_offset::make_pipe_skinned(&mut model, &spine, r, 1e-4, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // Pappus for a rotation-minimizing tube: area times spine length, to
    // second order in curvature times radius.
    let length = ogeom_algo::curve_length(&spine_curve, domain, T).unwrap();
    let pi = core::f64::consts::PI;
    let expected = pi * r * r * length;
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < expected * 0.01,
        "free-form pipe volume {measured} against {expected}"
    );
    assert!(!result.history.generated(&spine).is_empty());
}
