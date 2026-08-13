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
    assert!(err.to_string().contains("square to"), "{err}");
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
    // Positive narrows the solid in the pull direction, negative widens it
    // — the two wedges mirror about the undrafted volume.
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
