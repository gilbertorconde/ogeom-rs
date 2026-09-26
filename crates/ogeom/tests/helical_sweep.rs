//! The screw sweep: a profile in a plane through the axis, every point of
//! it running its own helix. Its volume is what Pappus says (the profile's
//! area times the path of its centroid round the axis) and it reaches
//! exactly as far as the profile does.
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
/// from the axis, lies inside for `acos(-x0 / r)` of every half turn.
#[test]
fn a_coil_trimmed_by_a_block_measures_its_closed_form() {
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Y,
    };
    for x0 in [-3.0_f64, 0.5, 2.0] {
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
