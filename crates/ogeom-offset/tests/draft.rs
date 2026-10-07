//! Draft on walls of revolution: the round-boss case.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Frame, Plane, Point};
use ogeom_topo::{Filter, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &ogeom_topo::Model, shape: &ogeom_topo::Shape) -> f64 {
    // Chord 1e-3: where the cone tessellation converges. A finer chord
    // mis-samples the slant, which this test does not pin.
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
    let (curve, points) = wavy_profile(amplitude, frequency);
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

/// The spline through a sine of `amplitude` and `frequency` along twenty
/// units of x in the XY plane, and the points it passes through.
fn wavy_profile(amplitude: f64, frequency: f64) -> (ogeom_geom::Curve, Vec<Point>) {
    let points: Vec<Point> = (0..=16)
        .map(|i| {
            let x = 20.0 * f64::from(i) / 16.0;
            Point::new(x, (x * frequency).sin() * amplitude, 0.0)
        })
        .collect();
    let spline = ogeom_geom::fit::fit_points(&points, 3, 1e-6, T)
        .unwrap()
        .curve;
    (ogeom_geom::Curve::BSpline(spline), points)
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

/// A drafted spline wall is the wall the draft names: every point of the
/// fitted face lies on a ruling through the profile, the pull turned
/// inwards about the profile's tangent by the angle, measured densely
/// against the rulings rather than at the samples the fit was made from.
#[test]
fn a_drafted_spline_wall_lies_on_its_turned_rulings() {
    use ogeom_geom::{Curve3d as _, Surface as _};
    let (amplitude, frequency, angle) = (1.0, 0.9, 0.1_f64);
    let mut model = ogeom_topo::Model::new();
    let solid = spline_prism(&mut model, amplitude, frequency);
    let wall = extruded_wall_of(&model, &solid);
    let drafted = ogeom_offset::apply_draft(
        &mut model,
        &solid,
        std::slice::from_ref(&wall),
        Plane::through(Point::ORIGIN, ogeom_math::Direction::Z),
        ogeom_math::Direction::Z,
        angle,
        T,
    )
    .unwrap()
    .shape;
    let diagnosis = ogeom_algo::check(&model, &drafted, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    // The rulings: through the profile point, the pull leaned towards the
    // material (away from the footprint's outside, +y here).
    let (profile, _) = wavy_profile(amplitude, frequency);
    let (t0, t1) = profile.domain();
    let ruling_at = |t: f64| -> (Point, ogeom_math::Vector) {
        let c = profile.point_at(t, T).unwrap();
        let d = profile.d1_at(t, T).unwrap();
        let out = ogeom_math::Vector::new(-d.y, d.x, 0.0);
        let out = out / out.magnitude();
        (
            c,
            ogeom_math::Vector::new(0.0, 0.0, angle.cos()) - out * angle.sin(),
        )
    };
    let off_ruling = |p: Point, t: f64| -> f64 {
        let (c, r) = ruling_at(t);
        let w = p - c;
        (w - r * w.dot(r)).magnitude()
    };
    let at_k = |k: usize| t0 + (t1 - t0) * f64::from(u32::try_from(k).unwrap()) / 8000.0;
    let rulings: Vec<(Point, ogeom_math::Vector)> =
        (0..=8000).map(|k| ruling_at(at_k(k))).collect();
    // The nearest ruling: the closest of the sampled ones, then the
    // distance minimized between its neighbours.
    let to_rulings = |p: Point| -> (f64, usize) {
        let mut best = (f64::INFINITY, 0);
        for (k, (c, r)) in rulings.iter().enumerate() {
            if (c.x - p.x).abs() > 2.0 {
                continue;
            }
            let w = p - *c;
            let off = (w - *r * w.dot(*r)).magnitude();
            if off < best.0 {
                best = (off, k);
            }
        }
        let k = best.1;
        let (mut a, mut b) = (at_k(k.saturating_sub(1)), at_k((k + 1).min(8000)));
        for _ in 0..60 {
            let (m1, m2) = (a + (b - a) / 3.0, b - (b - a) / 3.0);
            if off_ruling(p, m1) < off_ruling(p, m2) {
                b = m2;
            } else {
                a = m1;
            }
        }
        (off_ruling(p, f64::midpoint(a, b)).min(best.0), k)
    };
    let fitted = explore(&model, &drafted, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find_map(|f| {
            let d = model.node(&f)?.data().as_face()?.clone();
            match model.geometry().surface(d.surface)? {
                ogeom_geom::SurfaceGeometry::BSpline(b) => Some(b.clone()),
                _ => None,
            }
        })
        .expect("the drafted wall is fitted");
    // Away from the profile's ends, where the wall runs on past them to be
    // trimmed: the fitted wall against the rulings.
    let inner = 800..=rulings.len() - 801;
    let ((u0, u1), (v0, v1)) = fitted.domain();
    // Densely along each parameter in turn, whichever runs along the wall.
    let mut worst = 0.0_f64;
    for i in 0..=2000 {
        for j in 0..=6 {
            let (dense, sparse) = (f64::from(i) / 2000.0, f64::from(j) / 6.0);
            for (a, b) in [(dense, sparse), (sparse, dense)] {
                let (u, v) = (u0 + (u1 - u0) * a, v0 + (v1 - v0) * b);
                let (off, k) = to_rulings(fitted.point_at(u, v, T).unwrap());
                if inner.contains(&k) {
                    worst = worst.max(off);
                }
            }
        }
    }
    eprintln!("drafted wall off its rulings by {worst}");
    assert!(
        worst <= 1e-4,
        "the drafted wall strays {worst} from its rulings"
    );
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
    // seven microns off it; the skin's target says so. The aligned loft
    // skins by fitting where the plain one would build the drum exactly.
    let hints: Vec<Point> = [0.0, 5.0, 10.0]
        .iter()
        .map(|z| Point::new(10.0, 0.0, *z))
        .collect();
    let solid = ogeom_offset::make_loft_skinned_aligned(&mut model, &sections, &hints, 1e-2, T)
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

/// A drum drafted about a neutral plane tilted well off its axis: the
/// drafted wall is the ruled surface through the ellipse where the plane
/// cuts the drum, each ruling the pull turned by the draft about the
/// ellipse's tangent. Sampled densely along the wall, its hinge stands on
/// that ellipse and its rulings run that way out to the window's far
/// rows, within the draft's fit target of 1e-4.
#[test]
fn an_oblique_drafted_drum_holds_its_rulings_between_stations() {
    use ogeom_geom::Surface as _;
    let mut model = ogeom_topo::Model::new();
    let solid = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 10.0, 20.0, T)
        .unwrap()
        .shape;
    let wall = wall_of(&model, &solid);
    let tilt = 0.7_f64;
    let up = ogeom_math::Vector::new(tilt.sin(), 0.0, tilt.cos());
    let neutral = Plane::through(
        Point::new(0.0, 0.0, 10.0),
        ogeom_math::Direction::new(up, T).unwrap(),
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
    .unwrap()
    .shape;
    let diagnosis = ogeom_algo::check(&model, &drafted, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let surface = explore(&model, &drafted, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find_map(|f| {
            let d = model.node(&f)?.data().as_face()?.clone();
            match model.geometry().surface(d.surface)? {
                s @ ogeom_geom::SurfaceGeometry::BSpline(_) => Some(s.clone()),
                _ => None,
            }
        })
        .expect("the drafted wall is fitted");
    let ((u0, u1), (v0, v1)) = surface.domain();
    // The chart's `v` is the height along the ruling, the hinge at naught.
    let reach = v0.abs().max(v1.abs());
    let mut worst = 0.0_f64;
    for k in 0..=4000 {
        let u = u0 + (u1 - u0) * f64::from(k) / 4000.0;
        let hinge = surface.point_at(u, 0.0, T).unwrap();
        let off_drum = (hinge.x.hypot(hinge.y) - 10.0).abs();
        let off_plane = neutral.distance_to(hinge);
        // The exact ruling there: the pull turned by the draft about the
        // hinge's tangent, square to the drum's normal and the plane's.
        let radial = ogeom_math::Vector::new(hinge.x, hinge.y, 0.0);
        let tangent = radial.cross(up);
        let (_, along) = surface.d1_at(u, 0.0, T).unwrap();
        let along = along / along.magnitude();
        let mut lean = f64::INFINITY;
        for sense in [1.0, -1.0] {
            let turn = ogeom_math::Transform::rotation(
                ogeom_math::Axis::new(hinge, ogeom_math::Direction::new(tangent, T).unwrap()),
                angle * sense,
            );
            let ruling = turn.apply_vector(ogeom_math::Vector::Z);
            lean = lean.min(ruling.cross(along).magnitude());
        }
        worst = worst.max(off_drum.hypot(off_plane) + lean * reach);
    }
    eprintln!("oblique drafted drum off its rulings by {worst} out to {reach}");
    assert!(worst <= 1e-4, "the wall strays {worst} from its rulings");
}

/// The fitted surface of `solid`'s one B-spline face, with the face.
fn bspline_surface_of(
    model: &ogeom_topo::Model,
    solid: &ogeom_topo::Shape,
) -> (ogeom_topo::Shape, ogeom_geom::SurfaceGeometry) {
    explore(model, solid, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find_map(|f| {
            let d = model.node(&f)?.data().as_face()?.clone();
            match model.geometry().surface(d.surface)? {
                s @ ogeom_geom::SurfaceGeometry::BSpline(_) => Some((f, s.clone())),
                _ => None,
            }
        })
        .expect("a fitted face")
}

/// A skinned loft's wall closes on its seam to position only: its normal
/// turns across the seam, and the exact rulings of its draft with it, one
/// set for each side. The exact drafted wall there is the two sides' ruled
/// surfaces, each continued past the seam on its own tangent plane until
/// they meet. With the seam twisted up the wall and the neutral plane
/// tilted across it, the drafted wall near the seam lies on one side or
/// the other within the draft's fit target, though at the seam it stands
/// half the turn times the reach (7e-4) off either side's last ruling.
#[test]
fn a_draft_across_a_seam_closed_only_to_position_meets_both_sides() {
    use ogeom_geom::Surface as _;
    use ogeom_math::Vector;
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
    // The seam turns 0.02 about the axis per unit of height.
    let hints: Vec<Point> = [0.0_f64, 5.0, 10.0]
        .iter()
        .map(|z| Point::new(10.0 * (z * 0.02).cos(), 10.0 * (z * 0.02).sin(), *z))
        .collect();
    // A cubic through the skin's samples of a circle kinks where it
    // closes: the normal turns by 1.4e-3 across the seam.
    let solid = ogeom_offset::make_loft_skinned_aligned(&mut model, &sections, &hints, 1e-3, T)
        .unwrap()
        .shape;
    let (wall, skin) = bspline_surface_of(&model, &solid);
    // Tilted about the x axis, so the hinge crosses the seam climbing and
    // its tangent there is not square to the pull.
    let tilt = 0.15_f64;
    let up = Vector::new(0.0, -tilt.sin(), tilt.cos());
    let neutral = Plane::through(
        Point::new(0.0, 0.0, 3.0),
        ogeom_math::Direction::new(up, T).unwrap(),
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
    .unwrap()
    .shape;
    let diagnosis = ogeom_algo::check(&model, &drafted, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let (_, wall) = bspline_surface_of(&model, &drafted);

    // The exact hinge, ruling and hinge tangent (the way `u` runs) at the
    // skin's `u`: where its `u` column crosses the plane, the pull turned
    // about the crossing's tangent by the draft, in the sense that leans
    // the outward normal towards the pull.
    let ((su0, su1), (sv0, sv1)) = skin.domain();
    let exact = |u: f64| -> (Point, Vector, Vector) {
        let mut v = f64::midpoint(sv0, sv1);
        for _ in 0..12 {
            let p = skin.point_at(u, v, T).unwrap();
            let (_, dv) = skin.d1_at(u, v, T).unwrap();
            v -= neutral.signed_distance_to(p) / up.dot(dv);
        }
        let hinge = skin.point_at(u, v, T).unwrap();
        let (du, _) = skin.d1_at(u, v, T).unwrap();
        let mut normal = skin.normal_at(u, v, T).unwrap().vector();
        if normal.dot(Vector::new(hinge.x, hinge.y, 0.0)) < 0.0 {
            normal = -normal;
        }
        let mut tangent = normal.cross(up);
        tangent = tangent / tangent.magnitude();
        if tangent.dot(du) < 0.0 {
            tangent = -tangent;
        }
        let axis = ogeom_math::Axis::new(hinge, ogeom_math::Direction::new(tangent, T).unwrap());
        let turn = [angle, -angle]
            .into_iter()
            .map(|a| ogeom_math::Transform::rotation(axis, a))
            .max_by(|a, b| {
                let lean = |t: &ogeom_math::Transform| t.apply_vector(normal).z;
                lean(a).total_cmp(&lean(b))
            })
            .unwrap();
        (hinge, turn.apply_vector(Vector::Z), tangent)
    };
    let off_line = |p: Point, h: Point, r: Vector| (p - h).cross(r).magnitude();
    // Rulings sampled over a sixteenth of the skin either side of its
    // seam, alternately from each end.
    let near = (su1 - su0) / 16.0;
    let samples = 512;
    let step = near / f64::from(samples);
    let lines: Vec<(Point, Vector)> = (0..=samples)
        .flat_map(|k| [su0 + step * f64::from(k), su1 - step * f64::from(k)])
        .map(|u| {
            let (h, r, _) = exact(u);
            (h, r)
        })
        .collect();
    // Each side's continuation: its last ruling and hinge tangent, the
    // tangent pointing out past the seam.
    let ends = [(su0, -1.0), (su1, 1.0)].map(|(u, out)| {
        let (h, r, t) = exact(u);
        (h, r, t * out)
    });
    // The wedge the turn opens between the two last rulings, at the
    // window's furthest row.
    let ((u0, u1), (v0, v1)) = wall.domain();
    let reach = v0.abs().max(v1.abs());
    let gap = ends[0].1.cross(ends[1].1).magnitude() * reach;
    let off_exact = |p: Point| -> f64 {
        let (mut best, mut at) = (f64::INFINITY, 0usize);
        for (k, (h, r)) in lines.iter().enumerate() {
            let d = off_line(p, *h, *r);
            if d < best {
                (best, at) = (d, k);
            }
        }
        // Refined about the nearest sample, on its side.
        #[allow(clippy::cast_precision_loss)]
        let (u, out) = if at % 2 == 0 {
            (su0 + step * (at / 2) as f64, 1.0)
        } else {
            (su1 - step * (at / 2) as f64, -1.0)
        };
        let (mut lo, mut hi) = (u - step * out, u + step * out);
        let clamp = |x: f64| x.clamp(su0, su1);
        for _ in 0..30 {
            let (a, b) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
            let (ea, eb) = (exact(clamp(a)), exact(clamp(b)));
            if off_line(p, ea.0, ea.1) < off_line(p, eb.0, eb.1) {
                hi = b;
            } else {
                lo = a;
            }
        }
        let e = exact(clamp(f64::midpoint(lo, hi)));
        best = best.min(off_line(p, e.0, e.1));
        // Past either end, on that side's continuation: within its tangent
        // plane, out to twice the gap the turn opens at the window's rows.
        for (h, r, t) in ends {
            let normal = t.cross(r);
            let normal = normal / normal.magnitude();
            let out = r.cross(normal);
            if (0.0..=gap * 2.0).contains(&(p - h).dot(out)) {
                best = best.min((p - h).dot(normal).abs());
            }
        }
        best
    };
    let band = (u1 - u0) / 24.0;
    let (mut worst, mut own) = (0.0_f64, 0.0_f64);
    for k in 0..=100 {
        let f = band * f64::from(k) / 100.0;
        for u in [u0 + f, u1 - f] {
            for frac in [0.0, 0.5, 1.0] {
                let p = wall.point_at(u, v0 + (v1 - v0) * frac, T).unwrap();
                worst = worst.max(off_exact(p));
                if k == 0 {
                    own = own.max(off_line(p, ends[0].0, ends[0].1));
                }
            }
        }
    }
    eprintln!(
        "drafted skin off its sides by {worst} near the seam, {own} off one side's last ruling, \
         out to {reach}, gap {gap}"
    );
    assert!(worst <= 1e-4, "the wall strays {worst} from both sides");
}
