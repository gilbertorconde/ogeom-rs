//! The screw sweep: a planar profile, every point of it running its own
//! helix. Its volume is what Pappus says (the profile's area times the path
//! of its centroid round the axis) and it reaches exactly as far as the
//! profile does.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_face, make_polygon, shape_bounds, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::PlaneSurface;
use ogeom::math::{Axis, Direction, Frame, Plane, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

/// The rectangle x in [10, 12], z in [0, 2] on the XZ plane.
fn rectangle(model: &mut Model) -> Shape {
    let pts =
        [(10.0, 0.0), (12.0, 0.0), (12.0, 2.0), (10.0, 2.0)].map(|(x, z)| Point::new(x, 0.0, z));
    let wire = make_polygon(model, &pts, true, T).unwrap().shape;
    let plane = Plane::new(Frame::new(Point::ORIGIN, -Direction::Y, Direction::X, T).unwrap());
    make_face(model, PlaneSurface::new(plane).into(), &[wire], T)
        .unwrap()
        .shape
}

fn z_axis() -> Axis {
    Axis {
        location: Point::ORIGIN,
        direction: Direction::Z,
    }
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    let fine = Deflection {
        chord: 1e-3,
        angular: 0.02,
        ..Deflection::default()
    };
    volume_properties(model, shape, fine, T).unwrap().mass
}

#[test]
fn a_square_thread_is_as_big_as_pappus_says() {
    let mut model = Model::new();
    let profile = rectangle(&mut model);
    let thread =
        ogeom::offset::make_helical_sweep(&mut model, &profile, z_axis(), 5.0, 4.0, false, 0.0, T)
            .unwrap()
            .shape;
    assert!(check(&model, &thread, T).unwrap().is_valid());
    let want = 4.0 * 4.0 * core::f64::consts::TAU * 11.0;
    let v = volume(&model, &thread);
    assert!((v - want).abs() < want * 1e-4, "{v} against {want}");
    let bound = shape_bounds(&model, &thread, T).unwrap();
    let (lo, hi) = (bound.low().unwrap(), bound.high().unwrap());
    assert!(
        (hi.z - 22.0).abs() < 1e-3 && lo.z.abs() < 1e-3,
        "{lo:?} {hi:?}"
    );
}

#[test]
fn a_fine_pitch_over_ten_turns_measures() {
    let mut model = Model::new();
    let profile = rectangle(&mut model);
    let thread =
        ogeom::offset::make_helical_sweep(&mut model, &profile, z_axis(), 3.0, 10.0, true, 0.0, T)
            .unwrap()
            .shape;
    let want = 10.0 * 4.0 * core::f64::consts::TAU * 11.0;
    let v = volume(&model, &thread);
    assert!((v - want).abs() < want * 1e-4, "{v} against {want}");
}

#[test]
fn a_profile_as_tall_as_the_pitch_is_refused() {
    let mut model = Model::new();
    let profile = rectangle(&mut model);
    assert!(
        ogeom::offset::make_helical_sweep(&mut model, &profile, z_axis(), 2.0, 3.0, false, 0.0, T)
            .is_err()
    );
}

#[test]
fn a_round_wire_spring_is_as_big_as_pappus_says() {
    let mut model = Model::new();
    // A circle of radius 1 about (8, 0, 0) in the XZ plane.
    let frame = Frame::new(Point::new(8.0, 0.0, 0.0), -Direction::Y, Direction::X, T).unwrap();
    let circle = ogeom::math::Circle::new(frame, 1.0, T).unwrap();
    let edge = ogeom::algo::make_edge(
        &mut model,
        ogeom::geom::CircleCurve::new(circle).into(),
        (0.0, core::f64::consts::TAU),
        T,
    )
    .unwrap()
    .shape;
    let wire = ogeom::algo::make_wire(&mut model, &[edge], T)
        .unwrap()
        .shape;
    let plane = Plane::new(frame);
    let profile = make_face(&mut model, PlaneSurface::new(plane).into(), &[wire], T)
        .unwrap()
        .shape;
    let spring =
        ogeom::offset::make_helical_sweep(&mut model, &profile, z_axis(), 3.0, 3.0, false, 0.0, T)
            .unwrap()
            .shape;
    assert!(check(&model, &spring, T).unwrap().is_valid());
    let want = core::f64::consts::PI * 3.0 * core::f64::consts::TAU * 8.0;
    let v = volume(&model, &spring);
    assert!((v - want).abs() < want * 1e-3, "{v} against {want}");
}

#[test]
fn a_tapered_thread_builds_and_grows() {
    let mut model = Model::new();
    let profile = rectangle(&mut model);
    let thread =
        ogeom::offset::make_helical_sweep(&mut model, &profile, z_axis(), 5.0, 2.0, false, 0.5, T)
            .unwrap()
            .shape;
    assert!(check(&model, &thread, T).unwrap().is_valid());
    let bound = shape_bounds(&model, &thread, T).unwrap();
    assert!(bound.high().unwrap().x > 12.5);
}

#[test]
fn a_thread_measures_to_pappus_at_the_default_deflection() {
    let mut model = Model::new();
    let profile = rectangle(&mut model);
    let thread =
        ogeom::offset::make_helical_sweep(&mut model, &profile, z_axis(), 5.0, 4.0, false, 0.0, T)
            .unwrap()
            .shape;
    let want = 4.0 * 4.0 * core::f64::consts::TAU * 11.0;
    let v = volume_properties(&model, &thread, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((v - want).abs() < want * 1e-4, "{v} against {want}");
}

/// Under the Frenet frame a profile square to a helix rides the helix's own
/// screw: the solid is as big as the profile's area times its centroid's
/// path, and it reaches as far out as the profile does.
#[test]
fn a_frenet_pipe_on_a_helix_is_a_screw() {
    use ogeom::geom::{Curve3d as _, HelixCurve};
    let mut model = Model::new();
    let helix = HelixCurve::new(Frame::WORLD, 11.0, 5.0, 4.0).unwrap();
    let range = helix.domain();
    let start = helix.point_at(range.0, T).unwrap();
    let tangent = helix.d1_at(range.0, T).unwrap();
    let length = {
        let per_turn = (core::f64::consts::TAU * 11.0).hypot(5.0);
        per_turn * 4.0
    };
    let edge = ogeom::algo::make_edge(&mut model, helix.into(), range, T)
        .unwrap()
        .shape;
    let spine = ogeom::algo::make_wire(&mut model, &[edge], T)
        .unwrap()
        .shape;
    let normal = Direction::new(tangent, T).unwrap();
    // The square's sides: one radial, one square to it in the section.
    let radial = Direction::new(start - Point::new(0.0, 0.0, start.z), T).unwrap();
    let frame = Frame::new(start, normal, radial, T).unwrap();
    let (x, y) = (frame.x().vector(), frame.y().vector());
    let pts: Vec<Point> = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
        .iter()
        .map(|(a, b)| start + x * *a + y * *b)
        .collect();
    let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
    let face = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let pipe = ogeom::offset::make_pipe_shell(&mut model, &face, &spine, true, 1e-3, T)
        .unwrap()
        .shape;
    let v = volume(&model, &pipe);
    let want = 4.0 * length;
    assert!((v - want).abs() < want * 1e-3, "{v} against {want}");
}

/// A coil kept inside a block that cuts across its turns measures what it
/// is at the default deflection: the fitted sections the boolean leaves on
/// its walls meet within their stated tolerances, and the integral takes
/// them there rather than a coarse mesh. Each point of the section, `r` out
/// from the axis, lies inside for `acos(-x0 / r)` of every half turn. A
/// block whose face holds the axis (`x0 = 0`) keeps half the coil.
#[test]
fn a_coil_trimmed_by_a_block_measures_its_closed_form() {
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Y,
    };
    for x0 in [-3.0_f64, 0.0, 0.5, 2.0] {
        let mut model = Model::new();
        let pts =
            [(5.0, 0.0), (6.0, 0.0), (6.0, 1.0), (5.0, 1.0)].map(|(x, y)| Point::new(x, y, 0.0));
        let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
        let square = make_face(
            &mut model,
            PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
            &[wire],
            T,
        )
        .unwrap()
        .shape;
        let coil =
            ogeom::offset::make_helical_sweep(&mut model, &square, axis, 3.0, 4.0, false, 0.0, T)
                .unwrap()
                .shape;
        let frame = Frame::new(Point::new(x0, -5.0, -20.0), Direction::Z, Direction::X, T).unwrap();
        let size = (20.0 - x0, 30.0, 40.0);
        let block = ogeom::algo::make_box(&mut model, frame, size, T)
            .unwrap()
            .shape;
        let n = 20_000;
        let want: f64 = (0..n)
            .map(|i| {
                let r = 5.0 + (f64::from(i) + 0.5) / f64::from(n);
                8.0 * r * (x0 / r).clamp(-1.0, 1.0).acos()
            })
            .sum::<f64>()
            / f64::from(n);
        let kept = ogeom::boolean::common(&mut model, &block, &coil, T)
            .unwrap()
            .shape;
        let cut = ogeom::boolean::cut(&mut model, &block, &coil, T)
            .unwrap()
            .shape;
        let default = |shape: &Shape| {
            volume_properties(&model, shape, Deflection::default(), T)
                .unwrap()
                .mass
        };
        let (common, rest) = (default(&kept), default(&cut));
        assert!(
            (common - want).abs() < want * 2e-3,
            "x0 {x0}: common {common} against {want}"
        );
        let block_volume = size.0 * size.1 * size.2;
        assert!(
            (rest - (block_volume - want)).abs() < want * 2e-3,
            "x0 {x0}: cut {rest} against {}",
            block_volume - want
        );
    }
}

/// The block on the other side of the axis keeps the other half: the coil
/// is symmetric about the plane through its axis, turn for turn.
#[test]
fn a_coil_halved_by_a_block_on_its_axis_keeps_half_either_side() {
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Y,
    };
    let mut model = Model::new();
    let pts = [(5.0, 0.0), (6.0, 0.0), (6.0, 1.0), (5.0, 1.0)].map(|(x, y)| Point::new(x, y, 0.0));
    let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
    let square = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let coil =
        ogeom::offset::make_helical_sweep(&mut model, &square, axis, 3.0, 4.0, false, 0.0, T)
            .unwrap()
            .shape;
    let frame = Frame::new(
        Point::new(-20.0, -5.0, -20.0),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    let block = ogeom::algo::make_box(&mut model, frame, (20.0, 30.0, 40.0), T)
        .unwrap()
        .shape;
    let kept = ogeom::boolean::common(&mut model, &block, &coil, T)
        .unwrap()
        .shape;
    let v = volume_properties(&model, &kept, Deflection::default(), T)
        .unwrap()
        .mass;
    let want = core::f64::consts::PI * 4.0 * 5.5;
    assert!((v - want).abs() < want * 2e-3, "{v} against {want}");
}

/// No pitch and a taper is a flat spiral: the profile turns in its plane
/// square to the axis while moving out, three turns of a unit square two
/// out per turn. Its volume is the area times the arc its centroid runs,
/// `2 pi * turns * r0 + pi * taper * turns^2`, and it stays as high as the
/// profile.
#[test]
fn no_pitch_and_a_taper_is_a_flat_spiral() {
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Y,
    };
    let mut model = Model::new();
    let pts = [(5.0, 0.0), (6.0, 0.0), (6.0, 1.0), (5.0, 1.0)].map(|(x, y)| Point::new(x, y, 0.0));
    let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
    let square = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let spiral =
        ogeom::offset::make_helical_sweep(&mut model, &square, axis, 0.0, 3.0, false, 2.0, T)
            .unwrap()
            .shape;
    assert!(check(&model, &spiral, T).unwrap().is_valid());
    let v = volume(&model, &spiral);
    let pi = core::f64::consts::PI;
    let want = 2.0 * pi * 3.0 * 5.5 + pi * 2.0 * 9.0;
    assert!((v - want).abs() < want * 5e-3, "{v} against {want}");
    let b = ogeom::algo::tight_bounds(&model, &spiral, T).unwrap();
    assert!(b.low().unwrap().y > -1e-6 && b.high().unwrap().y < 1.0 + 1e-6);
    // A taper that does not clear the profile's width meets the last turn.
    assert!(
        ogeom::offset::make_helical_sweep(&mut model, &square, axis, 0.0, 3.0, false, 0.5, T)
            .is_err()
    );
}

/// A groove swept two turns down into a blind bore that a cylinder
/// primitive drilled. The primitive's wall is stored over exactly its own
/// height, and the groove's sections with it must still run out across the
/// bore's rim. The same bore made by extruding a circle is the reference.
#[test]
fn a_groove_cuts_into_a_bore_drilled_by_a_cylinder_primitive() {
    let cut = |primitive: bool| {
        let mut m = Model::new();
        let block = ogeom::algo::make_box(&mut m, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let frame = Frame::new(Point::new(10.0, 10.0, 2.0), Direction::Z, Direction::X, T).unwrap();
        let bore = if primitive {
            ogeom::algo::make_cylinder(&mut m, frame, 2.5, 8.0, T)
                .unwrap()
                .shape
        } else {
            let circle = ogeom::math::Circle::new(frame, 2.5, T).unwrap();
            let edge = ogeom::algo::make_edge(
                &mut m,
                ogeom::geom::CircleCurve::new(circle).into(),
                (0.0, core::f64::consts::TAU),
                T,
            )
            .unwrap()
            .shape;
            let wire = ogeom::algo::make_wire(&mut m, &[edge], T).unwrap().shape;
            let disc = make_face(
                &mut m,
                PlaneSurface::new(Plane::new(frame)).into(),
                &[wire],
                T,
            )
            .unwrap()
            .shape;
            ogeom::algo::make_prism(&mut m, &disc, ogeom::math::Vector::new(0.0, 0.0, 8.0), T)
                .unwrap()
                .shape
        };
        let drilled = ogeom::boolean::cut(&mut m, &block, &bore, T).unwrap().shape;
        let pts = [(2.4, 11.4), (3.0, 11.0625), (3.0, 10.9375), (2.4, 10.6)]
            .map(|(d, z)| Point::new(10.0 + d, 10.0, z));
        let wire = make_polygon(&mut m, &pts, true, T).unwrap().shape;
        let plane = Plane::new(
            Frame::new(Point::new(10.0, 10.0, 10.0), -Direction::Y, Direction::X, T).unwrap(),
        );
        let profile = make_face(&mut m, PlaneSurface::new(plane).into(), &[wire], T)
            .unwrap()
            .shape;
        let axis = Axis {
            location: Point::new(10.0, 10.0, 10.0),
            direction: -Direction::Z,
        };
        let groove =
            ogeom::offset::make_helical_sweep(&mut m, &profile, axis, 1.0, 2.0, false, 0.0, T)
                .unwrap()
                .shape;
        let result = ogeom::boolean::cut(&mut m, &drilled, &groove, T)
            .unwrap_or_else(|e| panic!("primitive {primitive}: {e}"))
            .shape;
        assert!(check(&m, &result, T).unwrap().is_valid());
        volume(&m, &result)
    };
    let (a, b) = (cut(true), cut(false));
    assert!((a - b).abs() < 1e-3, "{a} against {b}");
}

/// A disc lying level, square to the axis, climbing one pitch as it turns:
/// every level cut is the disc, so the solid holds the disc's area times
/// the rise, runs from 0 to the pitch along the axis, and reaches as far
/// out as the disc does.
#[test]
fn a_level_profile_climbs_the_helix() {
    use ogeom::algo::{make_edge, make_wire};
    use ogeom::geom::CircleCurve;
    use ogeom::math::Circle;
    let mut model = Model::new();
    let frame = Frame::new(Point::new(10.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let circle = CircleCurve::new(Circle::new(frame, 1.0, T).unwrap());
    let edge = make_edge(
        &mut model,
        circle.into(),
        (0.0, 2.0 * core::f64::consts::PI),
        T,
    )
    .unwrap()
    .shape;
    let wire = make_wire(&mut model, &[edge], T).unwrap().shape;
    let disc = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let sweep =
        ogeom::offset::make_helical_sweep(&mut model, &disc, z_axis(), 10.0, 1.0, false, 0.0, T)
            .unwrap()
            .shape;
    let diagnosis = check(&model, &sweep, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let want = core::f64::consts::PI * 10.0;
    let v = volume(&model, &sweep);
    assert!((v - want).abs() < want * 1e-3, "{v} against {want}");
    let mesh =
        ogeom::mesh::triangulate(&model, &sweep, Deflection::with_chord(1e-3).unwrap(), T).unwrap();
    let (low, high, reach) = mesh.positions.iter().fold(
        (f64::INFINITY, f64::NEG_INFINITY, 0.0_f64),
        |(lo, hi, r), p| (lo.min(p.z), hi.max(p.z), r.max(p.x.hypot(p.y))),
    );
    assert!(
        low.abs() < 1e-6 && (high - 10.0).abs() < 1e-6,
        "z {low} .. {high}"
    );
    assert!((reach - 11.0).abs() < 1e-3, "reaches {reach}");
}

/// A level disc the axis runs through has points with no helix to follow,
/// and is refused rather than swept through itself.
#[test]
fn a_level_profile_on_the_axis_is_refused() {
    use ogeom::algo::{make_edge, make_wire};
    use ogeom::geom::CircleCurve;
    use ogeom::math::Circle;
    let mut model = Model::new();
    let frame = Frame::new(Point::new(0.5, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let circle = CircleCurve::new(Circle::new(frame, 1.0, T).unwrap());
    let edge = make_edge(
        &mut model,
        circle.into(),
        (0.0, 2.0 * core::f64::consts::PI),
        T,
    )
    .unwrap()
    .shape;
    let wire = make_wire(&mut model, &[edge], T).unwrap().shape;
    let disc = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let refused =
        ogeom::offset::make_helical_sweep(&mut model, &disc, z_axis(), 10.0, 1.0, false, 0.0, T);
    assert!(
        matches!(refused, Err(ogeom::core::OgeomError::Construction(_))),
        "{refused:?}"
    );
}

/// A planar face through the Z axis, turned `theta` about it, from
/// (distance from the axis, height) pairs.
fn axial_face(model: &mut Model, pts: &[(f64, f64)], theta: f64) -> Shape {
    let (c, s) = (theta.cos(), theta.sin());
    let p: Vec<Point> = pts
        .iter()
        .map(|&(r, z)| Point::new(r * c, r * s, z))
        .collect();
    let wire = make_polygon(model, &p, true, T).unwrap().shape;
    let normal = Direction::new(ogeom::math::Vector::new(s, -c, 0.0), T).unwrap();
    let x = Direction::new(ogeom::math::Vector::new(c, s, 0.0), T).unwrap();
    let plane = Plane::new(Frame::new(Point::ORIGIN, normal, x, T).unwrap());
    make_face(model, PlaneSurface::new(plane).into(), &[wire], T)
        .unwrap()
        .shape
}

/// A thread groove cut into a chamfered rod is the same cut whatever angle
/// about the axis it starts at. The groove's sections run from one patch of
/// the sweep to the next, and each ends on the patch edge exactly, where
/// the next one starts; a section stopping short of it left the wire open
/// at some start angles.
#[test]
fn a_thread_cut_into_a_rod_is_the_same_at_any_start_angle() {
    let mut volumes = Vec::new();
    for degrees in [0.0_f64, 37.0, 91.0, 195.0, 286.0, 351.0] {
        let mut model = Model::new();
        let section = axial_face(
            &mut model,
            &[
                (0.0, 0.0),
                (2.2, 0.0),
                (2.5, 0.3),
                (2.5, 11.7),
                (2.2, 12.0),
                (0.0, 12.0),
            ],
            0.0,
        );
        let rod =
            ogeom::algo::make_revolution(&mut model, &section, z_axis(), core::f64::consts::TAU, T)
                .unwrap()
                .shape;
        let s = 1.2;
        let profile = axial_face(
            &mut model,
            &[
                (2.9, s - 0.36),
                (2.01, s - 0.05),
                (2.01, s + 0.05),
                (2.9, s + 0.36),
            ],
            degrees.to_radians(),
        );
        let groove = ogeom::offset::make_helical_sweep(
            &mut model,
            &profile,
            z_axis(),
            0.8,
            12.0,
            false,
            0.0,
            T,
        )
        .unwrap()
        .shape;
        let cut = ogeom::boolean::cut(&mut model, &rod, &groove, T)
            .unwrap_or_else(|e| panic!("at {degrees} degrees: {e}"))
            .shape;
        let diagnosis = check(&model, &cut, T).unwrap();
        assert!(diagnosis.is_valid(), "at {degrees} degrees: {diagnosis}");
        volumes.push((degrees, volume(&model, &cut)));
    }
    let first = volumes[0].1;
    for (degrees, v) in &volumes {
        assert!(
            (v - first).abs() < first * 1e-5,
            "{v} at {degrees} degrees against {first}"
        );
    }
}
