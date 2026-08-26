//! Pipes and lofts, against Pappus and the frustum formulae.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::Curve3d as _;
use ogeom_math::{Circle, Frame, Point};
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

/// A square profile face of side `side`, square to `tangent` at `centre`.
fn square_profile(
    model: &mut ogeom_topo::Model,
    centre: ogeom_math::Point,
    tangent: ogeom_math::Vector,
    side: f64,
) -> ogeom_topo::Shape {
    use ogeom_math::Direction;
    let normal = Direction::new(tangent, T).unwrap();
    let plane = ogeom_math::Plane::through(centre, normal);
    let frame = plane.frame();
    let h = side / 2.0;
    let corners: Vec<ogeom_math::Point> = [(-h, -h), (h, -h), (h, h), (-h, h)]
        .iter()
        .map(|(a, b)| centre + frame.x().vector() * *a + frame.y().vector() * *b)
        .collect();
    let wire = ogeom_algo::make_polygon(model, &corners, true, T)
        .unwrap()
        .shape;
    let surface: ogeom_geom::SurfaceGeometry =
        ogeom_geom::PlaneSurface::over(plane, (-side * 2.0, side * 2.0), (-side * 2.0, side * 2.0))
            .unwrap()
            .into();
    ogeom_algo::make_face(model, surface, std::slice::from_ref(&wire), T)
        .unwrap()
        .shape
}

/// The quarter arc of radius `r` about the origin in the xy plane, starting
/// at `(r, 0, 0)` heading `+y`.
fn quarter_arc(model: &mut ogeom_topo::Model, r: f64) -> ogeom_topo::Shape {
    let circle = ogeom_math::Circle::new(Frame::WORLD, r, T).unwrap();
    let curve: ogeom_geom::Curve = ogeom_geom::CircleCurve::new(circle).into();
    ogeom_algo::make_edge(model, curve, (0.0, core::f64::consts::FRAC_PI_2), T)
        .unwrap()
        .shape
}

#[test]
fn a_square_face_along_an_arc_sweeps_the_volume_pappus_names() {
    let mut model = ogeom_topo::Model::new();
    let r = 20.0;
    let spine = quarter_arc(&mut model, r);
    let start = Point::new(r, 0.0, 0.0);
    let profile = square_profile(&mut model, start, ogeom_math::Vector::Y, 4.0);

    let result =
        ogeom_offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // Pappus: the centroid rides the spine, so the volume is area times the
    // arc's own length.
    let expected = 16.0 * (core::f64::consts::FRAC_PI_2 * r);
    let measured = ogeom_algo::volume_properties(
        &model,
        &result.shape,
        ogeom_mesh::Deflection::with_chord(1e-3).unwrap(),
        T,
    )
    .unwrap()
    .mass;
    assert!(
        (measured - expected).abs() / expected < 0.01,
        "pipe shell volume {measured} against {expected}"
    );
    assert!(!result.history.generated(&spine).is_empty());
    assert!(!result.history.generated(&profile).is_empty());
}

#[test]
fn a_triangle_along_a_helix_makes_a_thread() {
    use ogeom_geom::Curve3d as _;
    let mut model = ogeom_topo::Model::new();
    let helix = ogeom_geom::HelixCurve::new(Frame::WORLD, 5.0, 4.0, 2.0).unwrap();
    let curve: ogeom_geom::Curve = helix.into();
    let domain = curve.domain();
    let length = {
        // Chord-sum over a fine sampling: the closed form is the hypotenuse
        // law, but measuring it keeps the test honest about the curve.
        let mut sum = 0.0;
        let mut last = curve.point_at(domain.0, T).unwrap();
        for i in 1..=512 {
            let t = domain.0 + (domain.1 - domain.0) * f64::from(i) / 512.0;
            let p = curve.point_at(t, T).unwrap();
            sum += last.distance(p);
            last = p;
        }
        sum
    };
    let start = curve.point_at(domain.0, T).unwrap();
    let tangent = curve.d1_at(domain.0, T).unwrap();
    let spine = ogeom_algo::make_edge(&mut model, curve, domain, T)
        .unwrap()
        .shape;

    // A triangular wire centred on the spine, square to it: the thread form.
    let normal = ogeom_math::Direction::new(tangent, T).unwrap();
    let plane = ogeom_math::Plane::through(start, normal);
    let frame = plane.frame();
    let corners: Vec<Point> = [(0.6, 0.0), (-0.3, 0.45), (-0.3, -0.45)]
        .iter()
        .map(|(a, b)| start + frame.x().vector() * *a + frame.y().vector() * *b)
        .collect();
    let profile = ogeom_algo::make_polygon(&mut model, &corners, true, T)
        .unwrap()
        .shape;

    let result =
        ogeom_offset::make_pipe_shell(&mut model, &profile, &spine, true, 1e-3, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // The triangle's area times the helix length brackets the thread: the
    // Frenet frames turn the profile with the spine, and the centroid rides
    // it, so the volume sits near the Pappus figure.
    let area = 0.9 * 0.45; // base 0.9, height 0.9, halved
    let expected = area * length;
    // Coarse deflection on purpose: a hundred-station helical skin at a
    // fine chord costs minutes and the claim here is a ten-percent band.
    let measured =
        ogeom_algo::volume_properties(&model, &result.shape, ogeom_mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    assert!(
        (measured - expected).abs() / expected < 0.1,
        "thread volume {measured} against {expected}"
    );
}

#[test]
fn a_profile_with_a_hole_sweeps_the_hole() {
    let mut model = ogeom_topo::Model::new();
    let r = 20.0;
    let spine = quarter_arc(&mut model, r);
    let start = Point::new(r, 0.0, 0.0);

    // A square face with a round hole: the hole must ride the sweep.
    let normal = ogeom_math::Direction::new(ogeom_math::Vector::Y, T).unwrap();
    let plane = ogeom_math::Plane::through(start, normal);
    let frame = plane.frame();
    let h = 2.0;
    let corners: Vec<Point> = [(-h, -h), (h, -h), (h, h), (-h, h)]
        .iter()
        .map(|(a, b)| start + frame.x().vector() * *a + frame.y().vector() * *b)
        .collect();
    let outer = ogeom_algo::make_polygon(&mut model, &corners, true, T)
        .unwrap()
        .shape;
    let hole_circle =
        ogeom_math::Circle::new(ogeom_math::Frame::about(start, normal), 1.0, T).unwrap();
    let hole_curve: ogeom_geom::Curve = ogeom_geom::CircleCurve::new(hole_circle).into();
    let hole_domain = ogeom_geom::Curve3d::domain(&hole_curve);
    let hole_edge = ogeom_algo::make_edge(&mut model, hole_curve, hole_domain, T)
        .unwrap()
        .shape;
    let hole = ogeom_algo::make_wire(&mut model, std::slice::from_ref(&hole_edge), T)
        .unwrap()
        .shape;
    let surface: ogeom_geom::SurfaceGeometry =
        ogeom_geom::PlaneSurface::over(plane, (-8.0, 8.0), (-8.0, 8.0))
            .unwrap()
            .into();
    let profile = ogeom_algo::make_face(&mut model, surface, &[outer, hole], T)
        .unwrap()
        .shape;

    let result =
        ogeom_offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let pi = core::f64::consts::PI;
    let expected = (16.0 - pi) * (core::f64::consts::FRAC_PI_2 * r);
    let measured = ogeom_algo::volume_properties(
        &model,
        &result.shape,
        ogeom_mesh::Deflection::with_chord(1e-3).unwrap(),
        T,
    )
    .unwrap()
    .mass;
    assert!(
        (measured - expected).abs() / expected < 0.01,
        "holed pipe shell volume {measured} against {expected}"
    );
}

#[test]
fn a_frenet_frame_on_a_straight_spine_is_refused_by_name() {
    let mut model = ogeom_topo::Model::new();
    let line =
        ogeom_geom::LineCurve::segment(Point::new(0.0, 0.0, 0.0), Point::new(0.0, 10.0, 0.0), T)
            .unwrap();
    let curve: ogeom_geom::Curve = line.into();
    let domain = ogeom_geom::Curve3d::domain(&curve);
    let spine = ogeom_algo::make_edge(&mut model, curve, domain, T)
        .unwrap()
        .shape;
    let profile = square_profile(
        &mut model,
        Point::new(0.0, 0.0, 0.0),
        ogeom_math::Vector::Y,
        2.0,
    );
    let err =
        ogeom_offset::make_pipe_shell(&mut model, &profile, &spine, true, 1e-3, T).unwrap_err();
    assert!(err.to_string().contains("Frenet"), "{err}");
}

#[test]
fn a_leaning_pipe_shell_profile_is_refused_by_name() {
    let mut model = ogeom_topo::Model::new();
    let spine = quarter_arc(&mut model, 20.0);
    // The profile's plane contains the start tangent instead of crossing it.
    let profile = square_profile(
        &mut model,
        Point::new(20.0, 0.0, 0.0),
        ogeom_math::Vector::Z,
        4.0,
    );
    let err =
        ogeom_offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T).unwrap_err();
    assert!(err.to_string().contains("leans"), "{err}");
}

#[test]
fn a_closed_loft_through_four_sections_is_a_watertight_ring() {
    // Four circles stood around a ring, each square to the ring's own
    // tangent, lofted closed: one face bounding itself both ways round, no
    // caps anywhere.
    let mut model = ogeom_topo::Model::new();
    let ring_r = 10.0;
    let mut sections = Vec::new();
    for i in 0..16 {
        let angle = core::f64::consts::TAU / 16.0 * f64::from(i);
        let centre = Point::new(ring_r * angle.cos(), ring_r * angle.sin(), 0.0);
        let tangent = ogeom_math::Vector::new(-angle.sin(), angle.cos(), 0.0);
        let normal = ogeom_math::Direction::new(tangent, T).unwrap();
        // Alignment is the caller's authorship: every section shares its
        // local axes, so the loop does not twist.
        let frame = Frame::new(centre, normal, ogeom_math::Direction::Z, T).unwrap();
        let circle = Circle::new(frame, 1.0, T).unwrap();
        let curve: ogeom_geom::Curve = ogeom_geom::CircleCurve::new(circle).into();
        let domain = curve.domain();
        let edge = ogeom_algo::make_edge(&mut model, curve, domain, T)
            .unwrap()
            .shape;
        sections.push(
            ogeom_algo::make_wire(&mut model, std::slice::from_ref(&edge), T)
                .unwrap()
                .shape,
        );
    }
    let result = ogeom_offset::make_loft_skinned_closed(&mut model, &sections, 2e-2, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // One face, closed both ways round.
    assert_eq!(
        explore(&model, &result.shape, Filter::OfType(ShapeType::Face))
            .unwrap()
            .len(),
        1
    );
    // Sixteen sections make the loop dense enough for the closed C1 solve;
    // the volume then sits close to the torus the ring approximates.
    let expected = core::f64::consts::PI * (core::f64::consts::TAU * ring_r);
    let measured = ogeom_algo::volume_properties(
        &model,
        &result.shape,
        ogeom_mesh::Deflection::with_chord(1e-3).unwrap(),
        T,
    )
    .unwrap()
    .mass;
    assert!(
        (measured - expected).abs() / expected < 0.05,
        "closed loft volume {measured} against {expected}"
    );
    for section in &sections {
        assert!(!result.history.generated(section).is_empty());
    }
}

#[test]
fn a_rectangle_lofted_to_a_point_is_the_pyramid_the_closed_form_names() {
    let mut model = ogeom_topo::Model::new();
    let corners = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(4.0, 0.0, 0.0),
        Point::new(4.0, 3.0, 0.0),
        Point::new(0.0, 3.0, 0.0),
    ];
    let base = ogeom_algo::make_polygon(&mut model, &corners, true, T)
        .unwrap()
        .shape;
    // A skew apex on purpose: pyramid walls are triangles wherever it sits.
    let apex = ogeom_algo::make_vertex(&mut model, Point::new(1.0, 1.0, 6.0)).shape;
    let result = ogeom_offset::make_loft(&mut model, &base, &apex, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let expected = 4.0 * 3.0 * 6.0 / 3.0;
    let measured =
        ogeom_algo::volume_properties(&model, &result.shape, ogeom_mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    assert!(
        (measured - expected).abs() < 1e-6,
        "pyramid volume {measured} against {expected}"
    );
    assert!(!result.history.generated(&base).is_empty());
    assert!(!result.history.generated(&apex).is_empty());
}

#[test]
fn a_circle_lofted_to_a_point_on_its_axis_is_a_cone() {
    let mut model = ogeom_topo::Model::new();
    let circle = Circle::new(Frame::WORLD, 2.0, T).unwrap();
    let curve: ogeom_geom::Curve = ogeom_geom::CircleCurve::new(circle).into();
    let domain = curve.domain();
    let ring = ogeom_algo::make_edge(&mut model, curve, domain, T)
        .unwrap()
        .shape;
    let base = ogeom_algo::make_wire(&mut model, std::slice::from_ref(&ring), T)
        .unwrap()
        .shape;
    let apex = ogeom_algo::make_vertex(&mut model, Point::new(0.0, 0.0, 5.0)).shape;
    let result = ogeom_offset::make_loft(&mut model, &base, &apex, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let pi = core::f64::consts::PI;
    let expected = pi * 4.0 * 5.0 / 3.0;
    let measured = ogeom_algo::volume_properties(
        &model,
        &result.shape,
        ogeom_mesh::Deflection::with_chord(1e-4).unwrap(),
        T,
    )
    .unwrap()
    .mass;
    assert!(
        (measured - expected).abs() < 1e-3,
        "cone volume {measured} against {expected}"
    );

    // Off the axis the cone is oblique, which is the skinned machinery's.
    let mut second = ogeom_topo::Model::new();
    let circle = Circle::new(Frame::WORLD, 2.0, T).unwrap();
    let curve: ogeom_geom::Curve = ogeom_geom::CircleCurve::new(circle).into();
    let domain = curve.domain();
    let ring = ogeom_algo::make_edge(&mut second, curve, domain, T)
        .unwrap()
        .shape;
    let base = ogeom_algo::make_wire(&mut second, std::slice::from_ref(&ring), T)
        .unwrap()
        .shape;
    let leaning = ogeom_algo::make_vertex(&mut second, Point::new(1.0, 0.0, 5.0)).shape;
    assert!(ogeom_offset::make_loft(&mut second, &base, &leaning, T).is_err());
}

#[test]
fn an_alignment_hint_untwists_a_loft() {
    // Two squares whose traversals start a corner apart: left to their own
    // starts the skin shears, with the hints it is the prism it should be.
    let build = |model: &mut ogeom_topo::Model, rotate: usize, z: f64| -> ogeom_topo::Shape {
        let corners = [
            Point::new(0.0, 0.0, z),
            Point::new(2.0, 0.0, z),
            Point::new(2.0, 2.0, z),
            Point::new(0.0, 2.0, z),
        ];
        let rotated: Vec<Point> = (0..4).map(|i| corners[(i + rotate) % 4]).collect();
        ogeom_algo::make_polygon(model, &rotated, true, T)
            .unwrap()
            .shape
    };
    let mut model = ogeom_topo::Model::new();
    let bottom = build(&mut model, 0, 0.0);
    let top = build(&mut model, 1, 3.0);
    let hints = [Point::new(0.0, 0.0, 0.0), Point::new(0.0, 0.0, 3.0)];
    let aligned = ogeom_offset::make_loft_skinned_aligned(
        &mut model,
        &[bottom.clone(), top.clone()],
        &hints,
        5e-2,
        T,
    )
    .unwrap();
    let volume_of = |model: &ogeom_topo::Model, shape: &ogeom_topo::Shape| {
        ogeom_algo::volume_properties(
            model,
            shape,
            ogeom_mesh::Deflection::with_chord(1e-3).unwrap(),
            T,
        )
        .unwrap()
        .mass
    };
    let straight = volume_of(&model, &aligned.shape);
    assert!(
        (straight - 12.0).abs() / 12.0 < 0.05,
        "aligned loft volume {straight} against 12"
    );

    // Without the hints the rows pair a corner apart and the skin twists —
    // visibly less volume, which is the defect the hint exists to fix.
    let twisted = ogeom_offset::make_loft_skinned(&mut model, &[bottom, top], 5e-2, T).unwrap();
    let sheared = volume_of(&model, &twisted.shape);
    assert!(
        sheared < straight * 0.95,
        "the twist should cost volume: {sheared} vs {straight}"
    );
}

#[test]
fn a_ruled_loft_between_tilted_polygons_still_builds() {
    // Non-parallel sections were never the refusal — only skew walls are.
    // A top square turned about the x axis keeps every wall planar.
    let mut model = ogeom_topo::Model::new();
    let bottom = ogeom_algo::make_polygon(
        &mut model,
        &[
            Point::new(0.0, 0.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
            Point::new(2.0, 2.0, 0.0),
            Point::new(0.0, 2.0, 0.0),
        ],
        true,
        T,
    )
    .unwrap()
    .shape;
    let angle = 25.0_f64.to_radians();
    let turn = |y: f64, z: f64| -> (f64, f64) {
        let (dy, dz) = (y - 1.0, z - 3.0);
        (
            angle.cos().mul_add(dy, -(angle.sin() * dz)) + 1.0,
            angle.sin().mul_add(dy, angle.cos() * dz) + 3.0,
        )
    };
    let corners: Vec<Point> = [(0.0_f64, 0.0_f64), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]
        .iter()
        .map(|(x, y)| {
            let (ty, tz) = turn(*y, 3.0);
            Point::new(*x, ty, tz)
        })
        .collect();
    let top = ogeom_algo::make_polygon(&mut model, &corners, true, T)
        .unwrap()
        .shape;
    let result = ogeom_offset::make_loft(&mut model, &bottom, &top, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let measured =
        ogeom_algo::volume_properties(&model, &result.shape, ogeom_mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    assert!(
        measured > 1.0,
        "the tilted loft encloses volume: {measured}"
    );
}

#[test]
fn a_pipe_shell_round_a_closed_circle_matches_the_torus() {
    // The case with an exact answer, which is what pins the holonomy
    // correction: a circular profile round a circular spine is a torus.
    let mut model = ogeom_topo::Model::new();
    let ring_r = 10.0;
    let circle = Circle::new(Frame::WORLD, ring_r, T).unwrap();
    let curve: ogeom_geom::Curve = ogeom_geom::CircleCurve::new(circle).into();
    let domain = curve.domain();
    let spine = ogeom_algo::make_edge(&mut model, curve, domain, T)
        .unwrap()
        .shape;

    let start = Point::new(ring_r, 0.0, 0.0);
    let normal = ogeom_math::Direction::new(ogeom_math::Vector::Y, T).unwrap();
    let section = Circle::new(
        Frame::new(start, normal, ogeom_math::Direction::Z, T).unwrap(),
        1.0,
        T,
    )
    .unwrap();
    let scurve: ogeom_geom::Curve = ogeom_geom::CircleCurve::new(section).into();
    let sdomain = scurve.domain();
    let sedge = ogeom_algo::make_edge(&mut model, scurve, sdomain, T)
        .unwrap()
        .shape;
    let profile = ogeom_algo::make_wire(&mut model, std::slice::from_ref(&sedge), T)
        .unwrap()
        .shape;

    let result =
        ogeom_offset::make_pipe_shell(&mut model, &profile, &spine, false, 5e-3, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let expected = 2.0 * core::f64::consts::PI * core::f64::consts::PI * ring_r;
    let measured = ogeom_algo::volume_properties(
        &model,
        &result.shape,
        ogeom_mesh::Deflection::with_chord(1e-3).unwrap(),
        T,
    )
    .unwrap()
    .mass;
    assert!(
        (measured - expected).abs() / expected < 0.01,
        "closed pipe shell volume {measured} against {expected}"
    );
    assert!(!result.history.generated(&spine).is_empty());
}

#[test]
fn a_round_profile_along_a_closed_square_spine_is_a_ring() {
    // The closed square spine, its corners rounded the way a real ring's
    // are: four straights and four quarter arcs, one G1 loop. The sharp
    // corner is refused by name — no skin can turn a section through a
    // finite angle over no arc — and the refusal test below pins that.
    let mut model = ogeom_topo::Model::new();
    let (half, r) = (8.0, 2.0);
    let flat = half - r;
    // Tangent points and corner arcs, walked counter-clockwise from the
    // middle of the +x side.
    let mut edges: Vec<ogeom_topo::Shape> = Vec::new();
    let mut vertices: Vec<(ogeom_topo::Shape, Point)> = Vec::new();
    let corner_centres = [
        Point::new(flat, flat, 0.0),
        Point::new(-flat, flat, 0.0),
        Point::new(-flat, -flat, 0.0),
        Point::new(flat, -flat, 0.0),
    ];
    // Each side's straight run, then the arc at its far corner.
    for (i, _) in corner_centres.iter().enumerate() {
        let angle = core::f64::consts::FRAC_PI_2 * f64::from(u8::try_from(i).unwrap());
        let (c, s_) = (angle.cos(), angle.sin());
        // Outward side direction and travel direction for side i.
        let out = ogeom_math::Vector::new(c, s_, 0.0);
        let along = ogeom_math::Vector::new(-s_, c, 0.0);
        let from = Point::new(0.0, 0.0, 0.0) + out * half - along * flat;
        let to = Point::new(0.0, 0.0, 0.0) + out * half + along * flat;
        vertices.push((ogeom_algo::make_vertex(&mut model, from).shape, from));
        vertices.push((ogeom_algo::make_vertex(&mut model, to).shape, to));
        let _ = &corner_centres[i];
    }
    for i in 0..4 {
        let (vf, pf) = vertices[2 * i].clone();
        let (vt, pt) = vertices[2 * i + 1].clone();
        let line = ogeom_geom::LineCurve::segment(pf, pt, T).unwrap();
        let curve: ogeom_geom::Curve = line.into();
        let domain = curve.domain();
        edges.push(
            ogeom_algo::make_edge_between(&mut model, curve, domain, &vf, &vt, T)
                .unwrap()
                .shape,
        );
        // The arc from this side's end to the next side's start, about the
        // shared corner centre.
        let (vn, _) = vertices[(2 * i + 2) % 8].clone();
        let centre = corner_centres[i];
        let frame = Frame::new(
            centre,
            ogeom_math::Direction::Z,
            ogeom_math::Direction::new(pt - centre, T).unwrap(),
            T,
        )
        .unwrap();
        let circle = Circle::new(frame, r, T).unwrap();
        let curve: ogeom_geom::Curve = ogeom_geom::CircleCurve::new(circle).into();
        edges.push(
            ogeom_algo::make_edge_between(
                &mut model,
                curve,
                (0.0, core::f64::consts::FRAC_PI_2),
                &vt,
                &vn,
                T,
            )
            .unwrap()
            .shape,
        );
    }
    let spine = ogeom_algo::make_wire(&mut model, &edges, T).unwrap().shape;

    // The profile at the first edge's own start, square to it.
    let start = vertices[0].1;
    let tangent = vertices[1].1 - vertices[0].1;
    let normal = ogeom_math::Direction::new(tangent, T).unwrap();
    let section = Circle::new(
        Frame::new(start, normal, ogeom_math::Direction::Z, T).unwrap(),
        1.0,
        T,
    )
    .unwrap();
    let scurve: ogeom_geom::Curve = ogeom_geom::CircleCurve::new(section).into();
    let sdomain = scurve.domain();
    let sedge = ogeom_algo::make_edge(&mut model, scurve, sdomain, T)
        .unwrap()
        .shape;
    let profile = ogeom_algo::make_wire(&mut model, std::slice::from_ref(&sedge), T)
        .unwrap()
        .shape;

    let result =
        ogeom_offset::make_pipe_shell(&mut model, &profile, &spine, false, 5e-2, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // One face bounding itself both ways round, and a volume in a coarse
    // band of Pappus: the corners are smoothed by the skin.
    assert_eq!(
        explore(&model, &result.shape, Filter::OfType(ShapeType::Face))
            .unwrap()
            .len(),
        1
    );
    // Pappus round the rounded square: perimeter = four flats and a full
    // circle of corner arcs.
    let perimeter = (2.0 * flat).mul_add(4.0, core::f64::consts::TAU * r);
    let expected = core::f64::consts::PI * perimeter;
    let measured = ogeom_algo::volume_properties(
        &model,
        &result.shape,
        ogeom_mesh::Deflection::with_chord(1e-3).unwrap(),
        T,
    )
    .unwrap()
    .mass;
    assert!(
        (measured - expected).abs() / expected < 0.02,
        "square ring volume {measured} against {expected}"
    );
}

#[test]
fn a_sharp_cornered_closed_spine_is_refused_by_name() {
    let mut model = ogeom_topo::Model::new();
    let corners = [
        Point::new(8.0, -8.0, 0.0),
        Point::new(8.0, 8.0, 0.0),
        Point::new(-8.0, 8.0, 0.0),
        Point::new(-8.0, -8.0, 0.0),
    ];
    let spine = ogeom_algo::make_polygon(&mut model, &corners, true, T)
        .unwrap()
        .shape;
    let start = corners[0];
    let normal = ogeom_math::Direction::new(corners[1] - corners[0], T).unwrap();
    let section = Circle::new(
        Frame::new(start, normal, ogeom_math::Direction::Z, T).unwrap(),
        1.0,
        T,
    )
    .unwrap();
    let scurve: ogeom_geom::Curve = ogeom_geom::CircleCurve::new(section).into();
    let sdomain = scurve.domain();
    let sedge = ogeom_algo::make_edge(&mut model, scurve, sdomain, T)
        .unwrap()
        .shape;
    let profile = ogeom_algo::make_wire(&mut model, std::slice::from_ref(&sedge), T)
        .unwrap()
        .shape;
    let err =
        ogeom_offset::make_pipe_shell(&mut model, &profile, &spine, false, 5e-2, T).unwrap_err();
    assert!(err.to_string().contains("sharp corner"), "{err}");
}

#[test]
fn an_l_spine_mitres_its_corner_and_the_runs_share_the_ring() {
    // An open cornered spine: two straight legs at a right angle. Each leg
    // sweeps as its own ruled wall between end rings, the corner's twin
    // stations throw both boundary rings onto the bisector plane, and the
    // sew joins the runs along that one mitred ring. The mitre passes
    // through the centreline corner, so Pappus prices the whole elbow at
    // area times the legs' summed length — exactly.
    let mut model = ogeom_topo::Model::new();
    let a = Point::new(0.0, 0.0, 0.0);
    let b = Point::new(20.0, 0.0, 0.0);
    let c = Point::new(20.0, 20.0, 0.0);
    let va = ogeom_algo::make_vertex(&mut model, a).shape;
    let vb = ogeom_algo::make_vertex(&mut model, b).shape;
    let vc = ogeom_algo::make_vertex(&mut model, c).shape;
    let seg = |model: &mut ogeom_topo::Model,
               f: (&ogeom_topo::Shape, Point),
               t: (&ogeom_topo::Shape, Point)|
     -> ogeom_topo::Shape {
        let line = ogeom_geom::LineCurve::segment(f.1, t.1, T).unwrap();
        let curve = ogeom_geom::Curve::Line(line);
        let d = ogeom_geom::Curve3d::domain(&curve);
        ogeom_algo::make_edge_between(model, curve, d, f.0, t.0, T)
            .unwrap()
            .shape
    };
    let e1 = seg(&mut model, (&va, a), (&vb, b));
    let e2 = seg(&mut model, (&vb, b), (&vc, c));
    let spine = ogeom_algo::make_wire(&mut model, &[e1, e2], T)
        .unwrap()
        .shape;
    let profile = square_profile(&mut model, a, ogeom_math::Vector::X, 4.0);

    let result =
        ogeom_offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - 640.0).abs() < 640.0 * 1e-3,
        "mitred elbow volume {measured} against 640"
    );
}

#[test]
fn a_corner_against_a_curved_leg_is_refused_by_name() {
    // The mitre is exact only where the leg is straight: a ruled wall
    // between end rings *is* the trimmed prism. Against an arc it is not,
    // and the refusal says so.
    let mut model = ogeom_topo::Model::new();
    let r = 20.0;
    // A quarter arc ending at (0, 20), then a straight leg heading +x —
    // a genuine corner between a curved run and a straight one.
    let arc_curve: ogeom_geom::Curve =
        ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(Frame::WORLD, r, T).unwrap()).into();
    let a = Point::new(r, 0.0, 0.0);
    let b = Point::new(0.0, r, 0.0);
    let c = Point::new(0.0, r + 20.0, 0.0);
    let va = ogeom_algo::make_vertex(&mut model, a).shape;
    let vb = ogeom_algo::make_vertex(&mut model, b).shape;
    let vc = ogeom_algo::make_vertex(&mut model, c).shape;
    let arc = ogeom_algo::make_edge_between(
        &mut model,
        arc_curve,
        (0.0, core::f64::consts::FRAC_PI_2),
        &va,
        &vb,
        T,
    )
    .unwrap()
    .shape;
    let line = ogeom_geom::LineCurve::segment(b, c, T).unwrap();
    // The arc leaves its end heading -x; the leg turns square up +y.
    let lcurve = ogeom_geom::Curve::Line(line);
    let ldomain = ogeom_geom::Curve3d::domain(&lcurve);
    let leg = ogeom_algo::make_edge_between(&mut model, lcurve, ldomain, &vb, &vc, T)
        .unwrap()
        .shape;
    let spine = ogeom_algo::make_wire(&mut model, &[arc, leg], T)
        .unwrap()
        .shape;
    let profile = square_profile(&mut model, a, ogeom_math::Vector::Y, 4.0);
    let err =
        ogeom_offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T).unwrap_err();
    assert!(err.to_string().contains("mitre"), "{err}");
}

/// A skinned loft ends at a point: circles narrowing to an apex, the apex
/// deliberately off every axis so no exact cone could stand in.
///
/// The reference is piecewise: Pappus's frustum from the two rings, plus the
/// cone from the top ring to the apex — whose shear off the axis changes
/// nothing, volume being shear-invariant. The skin smooths the crease where
/// the pieces meet, so the assertion carries the fit's honesty, not the
/// mesher's.
#[test]
fn a_skinned_loft_to_an_offset_point_measures_as_frustum_plus_cone() {
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
    let s0 = ring(&mut model, 10.0, 0.0);
    let s1 = ring(&mut model, 6.0, 6.0);
    let apex = ogeom_algo::make_vertex(&mut model, Point::new(2.0, 1.0, 15.0)).shape;
    let built = ogeom_offset::make_loft_skinned(&mut model, &[s0, s1, apex], 1e-2, T).unwrap();
    assert!(
        ogeom_algo::check(&model, &built.shape, T)
            .unwrap()
            .is_valid(),
        "the apex loft is a valid solid"
    );
    let frustum = core::f64::consts::PI * 6.0 / 3.0 * (100.0 + 60.0 + 36.0);
    let cone = core::f64::consts::PI * 36.0 * 9.0 / 3.0;
    let measured = volume(&model, &built.shape);
    let reference = frustum + cone;
    assert!(
        (measured - reference).abs() / reference < 0.01,
        "apex loft volume {measured} against {reference}"
    );
}

/// A wavy middle section skins: planarity is the caps' requirement, and only
/// the end sections carry caps.
#[test]
fn a_non_planar_middle_section_lofts() {
    let mut model = ogeom_topo::Model::new();
    let flat = |model: &mut ogeom_topo::Model, r: f64, z: f64| {
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
    // The middle ring waves out of plane: a trigonometric-spline circle
    // whose z oscillates, fitted closed.
    let wavy = {
        let n = 64_i32;
        let pts: Vec<Point> = (0..=n)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let a = core::f64::consts::TAU * f64::from(i % n) / f64::from(n);
                Point::new(8.0 * a.cos(), 8.0 * a.sin(), 5.0 + 0.5 * (3.0 * a).sin())
            })
            .collect();
        let fitted = ogeom_geom::fit::fit_points_closed(&pts, 3, 1e-3, T).unwrap();
        let curve = ogeom_geom::Curve::BSpline(fitted.curve);
        let domain = curve.domain();
        let edge = ogeom_algo::make_edge(&mut model, curve, domain, T)
            .unwrap()
            .shape;
        ogeom_algo::make_wire(&mut model, std::slice::from_ref(&edge), T)
            .unwrap()
            .shape
    };
    let s0 = flat(&mut model, 10.0, 0.0);
    let s2 = flat(&mut model, 9.0, 10.0);
    let built = ogeom_offset::make_loft_skinned(&mut model, &[s0, wavy, s2], 5e-2, T).unwrap();
    assert!(
        ogeom_algo::check(&model, &built.shape, T)
            .unwrap()
            .is_valid(),
        "the wavy loft is a valid solid"
    );
    // Coarse expectation only: between the r=10 and r=9 caps through an
    // r=8 waist, the volume sits between the two bounding cylinders.
    let v = volume(&model, &built.shape);
    let lo = core::f64::consts::PI * 64.0 * 10.0;
    let hi = core::f64::consts::PI * 100.0 * 10.0;
    assert!(
        lo < v && v < hi,
        "wavy loft volume {v} outside ({lo}, {hi})"
    );
}
