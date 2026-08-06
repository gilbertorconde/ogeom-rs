//! §C1 of `docs/PLAN.md`: a profile swept along a planar spine, the way a
//! moulding runs round a frame.
//!
//! Every claim is a closed form. A square spine with a square profile is
//! prisms along the sides and quarter-turn wedges at the corners, and its
//! volume is the sum of four boxes and four quarter-cylinders — which is a
//! volume nobody has to approximate.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::with_chord(1e-4).unwrap(), T)
        .unwrap()
        .mass
}

/// A closed square spine in `z = 0`, running counter-clockwise from the
/// corner nearest the origin.
fn square_spine(model: &mut Model, side: f64) -> Shape {
    ogeom::algo::make_polygon(
        model,
        &[
            Point::new(0.0, 0.0, 0.0),
            Point::new(side, 0.0, 0.0),
            Point::new(side, side, 0.0),
            Point::new(0.0, side, 0.0),
        ],
        true,
        T,
    )
    .unwrap()
    .shape
}

/// A closed rectangular profile standing on the spine's start, in the plane
/// `x = 0` — square to the spine's first run along `+x`, and containing the
/// spine's normal `z`, which is what the sweep requires. It stands *outside*
/// the spine, from `y = -outer` to `y = -inner`, and rises to `height`.
fn upright_profile(model: &mut Model, inner: f64, outer: f64, height: f64) -> Shape {
    ogeom::algo::make_polygon(
        model,
        &[
            Point::new(0.0, -outer, 0.0),
            Point::new(0.0, -inner, 0.0),
            Point::new(0.0, -inner, height),
            Point::new(0.0, -outer, height),
        ],
        true,
        T,
    )
    .unwrap()
    .shape
}

/// A square spine and an upright rectangular profile: four prisms and four
/// quarter-turn wedges, and the volume says so.
///
/// The spine runs counter-clockwise about `+z`, so every corner turns left by
/// a right angle and the wedge is a quarter of an annulus of the profile's own
/// cross-section — the join the 2D offset makes, for the same reason.
#[test]
fn a_square_spine_sweeps_prisms_and_quarter_turn_corners() {
    let mut model = Model::new();
    let side = 20.0;
    let (inner, outer, height) = (1.0, 3.0, 2.0);
    let spine = square_spine(&mut model, side);
    let profile = upright_profile(&mut model, inner, outer, height);

    let built = ogeom::offset::make_evolved(&mut model, &spine, &profile, T).unwrap();
    let measured = volume(&model, &built.shape);

    // Four straight runs of the profile's area over the spine's own length,
    // and four quarter-annulus wedges of that same area.
    let area = (outer - inner) * height;
    let pi = core::f64::consts::PI;
    let straight = 4.0 * area * side;
    let wedge = pi * outer.mul_add(outer, -(inner * inner)) / 4.0 * height;
    let want = straight + 4.0 * wedge;
    assert!(
        (measured - want).abs() < want * 1e-3,
        "four runs and four quarter turns: {measured} against {want}"
    );
    assert!(
        ogeom::algo::check(&model, &built.shape, T)
            .unwrap()
            .is_valid(),
        "and the moulding is a valid solid"
    );
}

/// The same spine given as a *face*. The profile is now an open wire that
/// starts and ends on the spine's own plane, so the sweep leaves that plane
/// open — and the face says to close it there. What comes back is a volume
/// where the wire spine would leave a shell.
#[test]
fn a_face_spine_closes_the_sweep_into_a_volume() {
    let mut model = Model::new();
    let side = 20.0;
    let (inner, outer, height) = (2.0, 5.0, 3.0);

    let outline = square_spine(&mut model, side);
    let plane = ogeom::geom::PlaneSurface::over(
        ogeom::math::Plane::through(Point::ORIGIN, Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom::algo::make_face(&mut model, plane.into(), &[outline], T)
        .unwrap()
        .shape;

    // An open profile: up the outside, across the top, and back down to the
    // plane on the inside. Its two ends are on `z = 0`.
    let profile = ogeom::algo::make_polygon(
        &mut model,
        &[
            Point::new(0.0, -outer, 0.0),
            Point::new(0.0, -outer, height),
            Point::new(0.0, -inner, height),
            Point::new(0.0, -inner, 0.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;

    let built = ogeom::offset::make_evolved(&mut model, &face, &profile, T).unwrap();
    assert_eq!(
        ogeom::topo::explore_unique(&model, &built.shape, ShapeType::Solid)
            .unwrap()
            .len(),
        1,
        "a face spine closes the sweep into a solid"
    );
    let measured = volume(&model, &built.shape);
    let area = (outer - inner) * height;
    let pi = core::f64::consts::PI;
    let want = 4.0f64.mul_add(
        pi * outer.mul_add(outer, -(inner * inner)) / 4.0 * height,
        4.0 * area * side,
    );
    assert!(
        (measured - want).abs() < want * 1e-3,
        "the same material, now closed: {measured} against {want}"
    );
}

/// An arc in the spine turns the profile about that arc's own axis, so the
/// swept piece is exact rather than chorded. A quarter-round corner spine —
/// two straight runs joined by an arc — is measured against the closed form
/// for its own annulus.
#[test]
fn an_arc_in_the_spine_turns_the_profile_about_its_own_axis() {
    let mut model = Model::new();
    let (inner, outer, height) = (1.0, 2.0, 1.5);
    let bend = 5.0;
    let run = 10.0;

    // Along +x to the start of the bend, a quarter arc of radius `bend`
    // centred at (run, bend, 0), then along +y.
    let start = Point::new(0.0, 0.0, 0.0);
    let arc_start = Point::new(run, 0.0, 0.0);
    let arc_end = Point::new(run + bend, bend, 0.0);
    let end = Point::new(run + bend, bend + run, 0.0);

    // One vertex per junction, so the three edges are a wire rather than
    // three edges that happen to touch.
    let a = ogeom::algo::make_vertex(&mut model, start).shape;
    let b = ogeom::algo::make_vertex(&mut model, arc_start).shape;
    let c = ogeom::algo::make_vertex(&mut model, arc_end).shape;
    let d = ogeom::algo::make_vertex(&mut model, end).shape;

    let straight = |model: &mut Model, from: Point, to: Point, v0: &Shape, v1: &Shape| {
        let direction = Direction::new(to - from, T).unwrap();
        let line = ogeom::geom::LineCurve::new(ogeom::math::Axis {
            location: from,
            direction,
        });
        ogeom::algo::make_edge_between(model, line.into(), (0.0, from.distance(to)), v0, v1, T)
            .unwrap()
            .shape
    };
    let first = straight(&mut model, start, arc_start, &a, &b);
    let centre = Point::new(run, bend, 0.0);
    let circle = ogeom::math::Circle::new(
        Frame::new(centre, Direction::Z, -Direction::Y, T).unwrap(),
        bend,
        T,
    )
    .unwrap();
    let quarter = core::f64::consts::FRAC_PI_2;
    let arc = ogeom::algo::make_edge_between(
        &mut model,
        ogeom::geom::CircleCurve::new(circle).into(),
        (0.0, quarter),
        &b,
        &c,
        T,
    )
    .unwrap()
    .shape;
    let last = straight(&mut model, arc_end, end, &c, &d);

    let spine = ogeom::algo::make_wire(&mut model, &[first, arc, last], T)
        .unwrap()
        .shape;

    let profile = upright_profile(&mut model, inner, outer, height);
    let built = ogeom::offset::make_evolved(&mut model, &spine, &profile, T).unwrap();
    let measured = volume(&model, &built.shape);

    // Two straight runs, plus a quarter turn of the profile about the bend's
    // axis. The turned piece is an annular quarter whose radii are the bend
    // offset by the profile's own reach — Pappus, exactly.
    let area = (outer - inner) * height;
    let pi = core::f64::consts::PI;
    let straight = area * (run + run);
    let turned = pi / 2.0 * ((bend + outer).powi(2) - (bend + inner).powi(2)) / 2.0 * height;
    let want = straight + turned;
    assert!(
        (measured - want).abs() < want * 1e-3,
        "two runs and one bend: {measured} against {want}"
    );
}

/// The refusals, by name. A profile that leans out of the plane containing
/// the spine's normal is not square to the spine, and a spine that leaves its
/// own plane has no square profile to carry at all.
#[test]
fn a_leaning_profile_and_a_twisted_spine_are_refused_by_name() {
    let mut model = Model::new();
    let spine = square_spine(&mut model, 10.0);
    // A profile lying *in* the spine's plane: its own normal is the spine's,
    // so the spine's normal is not in it.
    let flat = ogeom::algo::make_polygon(
        &mut model,
        &[
            Point::new(0.0, -3.0, 0.0),
            Point::new(0.0, -1.0, 0.0),
            Point::new(2.0, -1.0, 0.0),
            Point::new(2.0, -3.0, 0.0),
        ],
        true,
        T,
    )
    .unwrap()
    .shape;
    let refused = ogeom::offset::make_evolved(&mut model, &spine, &flat, T);
    let message = format!("{}", refused.unwrap_err());
    assert!(
        message.contains("square to the spine"),
        "the refusal names what is wrong: {message}"
    );

    // A spine that leaves its plane.
    let twisted = ogeom::algo::make_polygon(
        &mut model,
        &[
            Point::new(0.0, 0.0, 0.0),
            Point::new(10.0, 0.0, 0.0),
            Point::new(10.0, 10.0, 4.0),
            Point::new(0.0, 10.0, 0.0),
        ],
        true,
        T,
    )
    .unwrap()
    .shape;
    let profile = upright_profile(&mut model, 1.0, 3.0, 2.0);
    let refused = ogeom::offset::make_evolved(&mut model, &twisted, &profile, T);
    let message = format!("{}", refused.unwrap_err());
    assert!(
        message.contains("leaves its own plane"),
        "and so does this one: {message}"
    );

    // And an open profile along a *wire* spine: there is no plane to close it
    // against, so what it sweeps is a shell and the refusal says which spine
    // would make it a volume.
    let open = ogeom::algo::make_polygon(
        &mut model,
        &[
            Point::new(0.0, -3.0, 0.0),
            Point::new(0.0, -3.0, 2.0),
            Point::new(0.0, -1.0, 2.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    let refused = ogeom::offset::make_evolved(&mut model, &spine, &open, T);
    let message = format!("{}", refused.unwrap_err());
    assert!(
        message.contains("sweeps a shell, not a volume"),
        "and the third names the spine that would close it: {message}"
    );
}
