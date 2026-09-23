//! Draft on walls of revolution: the round-boss case.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Frame, Plane, Point};
use ogeom_topo::{Filter, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &ogeom_topo::Model, shape: &ogeom_topo::Shape) -> f64 {
    // Chord 1e-3 on purpose: the cone tessellation converges here, while
    // finer chords currently mis-sample the slant.
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

/// The face of `solid` whose surface is a cylinder.
fn wall_of(model: &ogeom_topo::Model, solid: &ogeom_topo::Shape) -> ogeom_topo::Shape {
    explore(model, solid, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find(|f| {
            model
                .node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Cylinder(_)))
        })
        .expect("the solid has a cylindrical wall")
}

#[test]
fn a_drafted_boss_wall_is_a_cone_holding_its_neutral_circle() {
    let mut model = ogeom_topo::Model::new();
    let plate = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 2.0), T).unwrap();
    let seat = Frame::new(
        Point::new(10.0, 10.0, 2.0),
        ogeom_math::Direction::Z,
        ogeom_math::Direction::X,
        T,
    )
    .unwrap();
    let boss = ogeom_algo::make_cylinder(&mut model, seat, 5.0, 10.0, T).unwrap();
    let joined = ogeom_bool::fuse(&mut model, &plate.shape, &boss.shape, T).unwrap();
    let wall = wall_of(&model, &joined.shape);

    let angle = 2.0_f64.to_radians();
    let neutral = Plane::through(Point::new(0.0, 0.0, 2.0), ogeom_math::Direction::Z);
    let result = ogeom_offset::apply_draft(
        &mut model,
        &joined.shape,
        std::slice::from_ref(&wall),
        neutral,
        ogeom_math::Direction::Z,
        angle,
        T,
    )
    .unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // The wall came back a cone at exactly the draft's half-angle, holding
    // radius five on the neutral plane.
    let cone = explore(&model, &result.shape, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find_map(|f| {
            model
                .node(&f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .and_then(|s| match s {
                    ogeom_geom::SurfaceGeometry::Cone(c) => Some(c.cone()),
                    _ => None,
                })
        })
        .expect("the drafted wall is a cone");
    assert!(
        (cone.half_angle().abs() - angle).abs() < 1e-9,
        "half angle {} against {angle}",
        cone.half_angle()
    );
    let neutral_height =
        (Point::new(10.0, 10.0, 2.0) - cone.frame().origin()).dot(cone.frame().z().vector());
    assert!(
        (cone.radius_at(neutral_height) - 5.0).abs() < 1e-9,
        "the base circle moved off five"
    );

    // Plate plus the frustum the leaning boss became.
    let pi = core::f64::consts::PI;
    let r1 = 10.0_f64.mul_add(-angle.tan(), 5.0);
    let expected = 800.0 + pi * 10.0 / 3.0 * (5.0_f64.mul_add(5.0 + r1, r1 * r1));
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 0.15,
        "drafted boss volume {measured} against {expected}"
    );
}

#[test]
fn a_bare_drum_drafts_into_the_frustum_band_and_all() {
    // The seamed wall takes the wholesale band rebuild on the turned cone.
    let mut model = ogeom_topo::Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T).unwrap();
    let wall = wall_of(&model, &drum.shape);
    let angle = 2.0_f64.to_radians();
    let neutral = Plane::through(Point::new(0.0, 0.0, 0.0), ogeom_math::Direction::Z);
    let result = ogeom_offset::apply_draft(
        &mut model,
        &drum.shape,
        std::slice::from_ref(&wall),
        neutral,
        ogeom_math::Direction::Z,
        angle,
        T,
    )
    .unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    let pi = core::f64::consts::PI;
    let r1 = 10.0_f64.mul_add(-angle.tan(), 5.0);
    let expected = pi * 10.0 / 3.0 * (5.0_f64.mul_add(5.0 + r1, r1 * r1));
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 0.15,
        "drafted drum volume {measured} against {expected}"
    );
}

#[test]
fn a_neutral_plane_along_the_axis_is_refused_by_name() {
    let mut model = ogeom_topo::Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T).unwrap();
    let wall = wall_of(&model, &drum.shape);
    let neutral = Plane::through(Point::new(0.0, 0.0, 0.0), ogeom_math::Direction::X);
    let err = ogeom_offset::apply_draft(
        &mut model,
        &drum.shape,
        std::slice::from_ref(&wall),
        neutral,
        ogeom_math::Direction::Z,
        2.0_f64.to_radians(),
        T,
    )
    .unwrap_err();
    // A plane through the axis cuts the wall in two lines: two hinges,
    // and no one draft.
    assert!(err.to_string().contains("more than once"), "{err}");
}

#[test]
fn a_draft_that_swallows_the_apex_is_refused_by_name() {
    let mut model = ogeom_topo::Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T).unwrap();
    let wall = wall_of(&model, &drum.shape);
    let neutral = Plane::through(Point::new(0.0, 0.0, 0.0), ogeom_math::Direction::Z);
    let err = ogeom_offset::apply_draft(
        &mut model,
        &drum.shape,
        std::slice::from_ref(&wall),
        neutral,
        ogeom_math::Direction::Z,
        40.0_f64.to_radians(),
        T,
    )
    .unwrap_err();
    assert!(err.to_string().contains("apex"), "{err}");
}

#[test]
fn the_angles_sign_picks_inward_or_outward() {
    // Positive narrows the solid in the pull direction, negative widens it;
    // the two wedges mirror about the undrafted volume.
    let mut measured = Vec::new();
    for angle in [10.0_f64, -10.0] {
        let mut model = ogeom_topo::Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T).unwrap();
        let wall = explore(&model, &block.shape, Filter::OfType(ShapeType::Face))
            .unwrap()
            .into_iter()
            .find(|f| {
                let d = model.node(f).and_then(|n| n.data().as_face());
                let Some(d) = d else { return false };
                let Some(ogeom_geom::SurfaceGeometry::Plane(p)) =
                    model.geometry().surface(d.surface)
                else {
                    return false;
                };
                let placed = f.transform(model.datums()).unwrap();
                let n = placed.apply_vector(p.plane().normal().vector());
                n.z.abs() < 1e-9 && (n.y - 1.0).abs() < 1e-9
            })
            .expect("the +y wall");
        let neutral = Plane::through(Point::new(0.0, 0.0, 0.0), ogeom_math::Direction::Z);
        let drafted = ogeom_offset::apply_draft(
            &mut model,
            &block.shape,
            std::slice::from_ref(&wall),
            neutral,
            ogeom_math::Direction::Z,
            angle.to_radians(),
            T,
        )
        .unwrap();
        measured.push(volume(&model, &drafted.shape));
    }
    let wedge = 10.0 * (10.0 * 10.0 / 2.0) * 10.0_f64.to_radians().tan();
    assert!(
        (measured[0] - (1000.0 - wedge)).abs() < 1e-3,
        "inward volume {} against {}",
        measured[0],
        1000.0 - wedge
    );
    assert!(
        (measured[1] - (1000.0 + wedge)).abs() < 1e-3,
        "outward volume {} against {}",
        measured[1],
        1000.0 + wedge
    );
}

#[test]
fn a_negative_draft_widens_the_drum() {
    // The revolved path honours the sign the same way: the wall flares
    // outward into the frustum whose base sits on the neutral circle.
    let mut model = ogeom_topo::Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T).unwrap();
    let wall = wall_of(&model, &drum.shape);
    let angle = (-2.0_f64).to_radians();
    let neutral = Plane::through(Point::new(0.0, 0.0, 0.0), ogeom_math::Direction::Z);
    let result = ogeom_offset::apply_draft(
        &mut model,
        &drum.shape,
        std::slice::from_ref(&wall),
        neutral,
        ogeom_math::Direction::Z,
        angle,
        T,
    )
    .unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let pi = core::f64::consts::PI;
    let r1 = 10.0_f64.mul_add(2.0_f64.to_radians().tan(), 5.0);
    let expected = pi * 10.0 / 3.0 * 5.0_f64.mul_add(5.0 + r1, r1 * r1);
    let measured = volume(&model, &result.shape);
    assert!(
        (measured - expected).abs() < 0.15,
        "widened drum volume {measured} against {expected}"
    );
}

/// A wavy spline profile closed into a slab
/// footprint, extruded along `z`, whose front wall is an extruded-spline
/// surface. `amplitude` and `frequency` set how tightly the wall curls.
fn spline_prism(
    model: &mut ogeom_topo::Model,
    amplitude: f64,
    frequency: f64,
) -> ogeom_topo::Shape {
    use ogeom_geom::Curve;
    let points: Vec<Point> = (0..=16)
        .map(|i| {
            let x = 20.0 * f64::from(i) / 16.0;
            Point::new(x, (x * frequency).sin() * amplitude, 0.0)
        })
        .collect();
    let spline = ogeom_geom::fit::fit_points(&points, 3, 1e-6, T)
        .unwrap()
        .curve;
    let curve: Curve = Curve::BSpline(spline);
    let dom = {
        use ogeom_geom::Curve3d as _;
        curve.domain()
    };
    let a = ogeom_algo::make_vertex(model, points[0]).shape;
    let b = ogeom_algo::make_vertex(model, points[16]).shape;
    let c = ogeom_algo::make_vertex(model, Point::new(20.0, -8.0, 0.0)).shape;
    let d = ogeom_algo::make_vertex(model, Point::new(0.0, -8.0, 0.0)).shape;
    let e_spline = ogeom_algo::make_edge_between(model, curve, dom, &a, &b, T)
        .unwrap()
        .shape;
    let seg = |m: &mut ogeom_topo::Model, p: Point, q: Point, vp, vq| {
        let line: Curve = Curve::Line(ogeom_geom::LineCurve::segment(p, q, T).unwrap());
        let ld = {
            use ogeom_geom::Curve3d as _;
            line.domain()
        };
        ogeom_algo::make_edge_between(m, line, ld, vp, vq, T)
            .unwrap()
            .shape
    };
    let e1 = seg(model, points[16], Point::new(20.0, -8.0, 0.0), &b, &c);
    let e2 = seg(
        model,
        Point::new(20.0, -8.0, 0.0),
        Point::new(0.0, -8.0, 0.0),
        &c,
        &d,
    );
    let e3 = seg(model, Point::new(0.0, -8.0, 0.0), points[0], &d, &a);
    let plane = Plane::through(Point::ORIGIN, ogeom_math::Direction::Z);
    // Counter-clockwise about +z, so the face's material is the footprint.
    let face = ogeom_algo::make_face_with_pcurves(
        model,
        ogeom_geom::PlaneSurface::over(plane, (-40.0, 40.0), (-40.0, 40.0))
            .unwrap()
            .into(),
        &[vec![
            e3.reversed(),
            e2.reversed(),
            e1.reversed(),
            e_spline.reversed(),
        ]],
        T,
    )
    .unwrap()
    .shape;
    ogeom_algo::make_prism(model, &face, ogeom_math::Vector::new(0.0, 0.0, 10.0), T)
        .unwrap()
        .shape
}

/// The face of `solid` on an extrusion surface.
fn extruded_wall_of(model: &ogeom_topo::Model, solid: &ogeom_topo::Shape) -> ogeom_topo::Shape {
    explore(model, solid, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find(|f| {
            model
                .node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Extrusion(_)))
        })
        .expect("the solid has an extruded wall")
}

#[test]
fn an_extruded_spline_wall_drafts_to_the_requested_angle() {
    use ogeom_geom::Surface as _;
    let mut model = ogeom_topo::Model::new();
    let solid = spline_prism(&mut model, 1.5, 0.4);
    let before = volume(&model, &solid);
    let wall = extruded_wall_of(&model, &solid);
    let angle = 0.1_f64;
    let drafted = ogeom_offset::apply_draft(
        &mut model,
        &solid,
        std::slice::from_ref(&wall),
        Plane::through(Point::ORIGIN, ogeom_math::Direction::Z),
        ogeom_math::Direction::Z,
        angle,
        T,
    )
    .unwrap();
    let after = volume(&model, &drafted.shape);
    assert!(
        after < before && before - after < before * 0.2,
        "a draft shaves a wedge: {before} -> {after}"
    );

    // The drafted wall is the fitted face; its normal leans off the pull by
    // exactly the requested angle, at three sampled heights.
    let fitted = explore(&model, &drafted.shape, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find(|f| {
            model
                .node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::BSpline(_)))
        })
        .expect("the drafted wall is fitted");
    let surface = {
        let d = model
            .node(&fitted)
            .unwrap()
            .data()
            .as_face()
            .unwrap()
            .clone();
        model.geometry().surface(d.surface).unwrap().clone()
    };
    let ((u0, u1), (v0, v1)) = surface.domain();
    for frac in [0.25, 0.5, 0.75] {
        let (u, v) = (f64::midpoint(u0, u1), v0 + (v1 - v0) * frac);
        let (du, dv) = surface.d1_at(u, v, T).unwrap();
        let n = du.cross(dv);
        let lean = (n / n.magnitude())
            .dot(ogeom_math::Vector::new(0.0, 0.0, 1.0))
            .asin()
            .abs();
        assert!(
            (lean - angle).abs() < 1e-4,
            "the wall leans {lean} at height {frac}, wanted {angle}"
        );
    }
}

#[test]
fn a_draft_that_folds_the_wall_refuses_by_name() {
    // A profile curled tighter than the draft's reach: the turned rulings
    // cross inside the drafted window, and the fold is refused before
    // anything is fitted.
    let mut model = ogeom_topo::Model::new();
    let solid = spline_prism(&mut model, 2.0, 0.7);
    let wall = extruded_wall_of(&model, &solid);
    let err = ogeom_offset::apply_draft(
        &mut model,
        &solid,
        std::slice::from_ref(&wall),
        Plane::through(Point::ORIGIN, ogeom_math::Direction::Z),
        ogeom_math::Direction::Z,
        0.1,
        T,
    )
    .unwrap_err()
    .to_string();
    assert!(
        err.contains("folds the wall"),
        "the fold names itself: {err}"
    );
}

/// The angle a fitted wall's rulings make with `pull`, sampled at three
/// heights up the middle of its chart: the draft angle, by definition: a
/// drafted wall is ruled along the pull turned by the draft, whatever its
/// hinge does.
fn leans_of(
    model: &ogeom_topo::Model,
    solid: &ogeom_topo::Shape,
    pull: ogeom_math::Vector,
) -> Vec<f64> {
    use ogeom_geom::Surface as _;
    let fitted = explore(model, solid, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find(|f| {
            model
                .node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::BSpline(_)))
        })
        .expect("the drafted wall is fitted");
    let surface = {
        let d = model
            .node(&fitted)
            .unwrap()
            .data()
            .as_face()
            .unwrap()
            .clone();
        model.geometry().surface(d.surface).unwrap().clone()
    };
    let ((u0, u1), (v0, v1)) = surface.domain();
    [0.25, 0.5, 0.75]
        .into_iter()
        .map(|frac| {
            let (u, v) = (f64::midpoint(u0, u1), v0 + (v1 - v0) * frac);
            let (_, dv) = surface.d1_at(u, v, T).unwrap();
            (dv / dv.magnitude()).dot(pull).abs().acos()
        })
        .collect()
}

/// A drum drafted about a neutral plane *tilted* against its axis: the
/// hinge is an ellipse, the drafted wall the ruled surface through it with
/// every ruling at the draft angle from the pull.
#[test]
fn a_drum_drafts_about_an_oblique_neutral() {
    let mut model = ogeom_topo::Model::new();
    let solid = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 10.0, 20.0, T)
        .unwrap()
        .shape;
    let before = volume(&model, &solid);
    let wall = wall_of(&model, &solid);
    let tilt = 0.35_f64;
    let neutral = Plane::through(
        Point::new(0.0, 0.0, 10.0),
        ogeom_math::Direction::new(ogeom_math::Vector::new(tilt.sin(), 0.0, tilt.cos()), T)
            .unwrap(),
    );
    let angle = 0.1_f64;
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
    // The hinge sits mid-height: the wall narrows above it and widens
    // below, and the two nearly cancel. What a draft is measured by is its
    // angle.
    let after = volume(&model, &drafted.shape);
    assert!(
        (after - before).abs() < before * 0.02,
        "a mid-height draft keeps the volume: {before} -> {after}"
    );
    for lean in leans_of(
        &model,
        &drafted.shape,
        ogeom_math::Vector::new(0.0, 0.0, 1.0),
    ) {
        assert!(
            (lean - angle).abs() < 2e-3,
            "the wall leans {lean} off the pull, wanted {angle}"
        );
    }
}

/// A wall on a raw fitted patch (a skinned loft's, with no ruling to turn)
/// drafts the same way, and comes out the frustum a drafted cylinder is.
#[test]
fn a_fitted_patch_wall_drafts_to_the_requested_angle() {
    let mut model = ogeom_topo::Model::new();
    let ring = |model: &mut ogeom_topo::Model, z: f64| {
        let frame = Frame::new(
            Point::new(0.0, 0.0, z),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap();
        let circle = ogeom_math::Circle::new(frame, 10.0, T).unwrap();
        let curve = ogeom_geom::Curve::Circle(ogeom_geom::CircleCurve::new(circle));
        let domain = {
            use ogeom_geom::Curve3d as _;
            curve.domain()
        };
        let edge = ogeom_algo::make_edge(model, curve, domain, T)
            .unwrap()
            .shape;
        ogeom_algo::make_wire(model, std::slice::from_ref(&edge), T)
            .unwrap()
            .shape
    };
    let sections = [
        ring(&mut model, 0.0),
        ring(&mut model, 5.0),
        ring(&mut model, 10.0),
    ];
    // A cubic through forty-eight samples of a circle of radius ten sits
    // seven microns off it; the skin's target says so.
    let solid = ogeom_offset::make_loft_skinned(&mut model, &sections, 1e-2, T)
        .unwrap()
        .shape;
    let wall = explore(&model, &solid, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find(|f| {
            model
                .node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::BSpline(_)))
        })
        .expect("the skinned wall");
    let angle = 0.1_f64;
    let drafted = ogeom_offset::apply_draft(
        &mut model,
        &solid,
        std::slice::from_ref(&wall),
        Plane::through(Point::ORIGIN, ogeom_math::Direction::Z),
        ogeom_math::Direction::Z,
        angle,
        T,
    )
    .unwrap();
    let diagnosis = ogeom_algo::check(&model, &drafted.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    // The base circle is the hinge, so the top narrows to 10 − 10 tan(0.1):
    // the frustum's volume, to what the fit resolves.
    let top = 10.0 - 10.0 * angle.tan();
    let frustum = core::f64::consts::PI * 10.0 / 3.0 * (100.0 + 10.0 * top + top * top);
    let after = volume(&model, &drafted.shape);
    assert!(
        (after - frustum).abs() < frustum * 5e-3,
        "the drafted skin measures {after} against the frustum's {frustum}"
    );
    for lean in leans_of(
        &model,
        &drafted.shape,
        ogeom_math::Vector::new(0.0, 0.0, 1.0),
    ) {
        assert!(
            (lean - angle).abs() < 2e-3,
            "the wall leans {lean} off the pull, wanted {angle}"
        );
    }
}
