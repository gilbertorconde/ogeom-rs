//! Rotation-minimizing frames along curves and wires, against closed forms:
//! a helix turns against its Frenet frame at minus its torsion, a circle's
//! frame is its Frenet frame turned once, a line's is constant, a planar
//! curve's out-of-plane normal stays put through an inflection, a
//! polygon's turns are the least rotations at its corners and its
//! holonomy their product, and a skew spline's frames and holonomy match
//! the frame equation integrated on its own.
#![allow(clippy::unwrap_used, clippy::cast_precision_loss, reason = "test code")]

use core::f64::consts::{PI, TAU};

use ogeom_algo::{CurveFrame, make_edge, make_polygon, rotation_minimizing_frames};
use ogeom_core::Tolerances;
use ogeom_geom::fit::{Spacing, interpolate};
use ogeom_geom::{CircleCurve, Curve, Curve3d as _, HelixCurve, LineCurve};
use ogeom_math::{Circle, Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn edge_of(model: &mut Model, curve: Curve) -> Shape {
    let domain = curve.domain();
    make_edge(model, curve, domain, T).unwrap().shape
}

fn shares(n: usize) -> Vec<f64> {
    (0..=n).map(|k| k as f64 / n as f64).collect()
}

/// The signed angle from `a` to `b` about `axis`, all square to `axis`.
fn angle_about(a: Vector, b: Vector, axis: Vector) -> f64 {
    a.cross(b).dot(axis).atan2(a.dot(b))
}

/// The difference of two angles, wrapped into `(-π, π]`.
fn wrapped(a: f64) -> f64 {
    (a + PI).rem_euclid(TAU) - PI
}

/// Each frame unit, square and right-handed to `1e-12`.
fn assert_orthonormal(frames: &[CurveFrame]) {
    for f in frames {
        let (t, n, b) = (f.tangent.vector(), f.normal.vector(), f.binormal());
        assert!((t.magnitude() - 1.0).abs() <= 1e-12, "{f:?}");
        assert!((n.magnitude() - 1.0).abs() <= 1e-12, "{f:?}");
        assert!(t.dot(n).abs() <= 1e-12, "tangent·normal {}", t.dot(n));
        assert!((b.magnitude() - 1.0).abs() <= 1e-12, "{f:?}");
    }
}

/// Rodrigues: `v` turned about the unit `axis` by `angle`.
fn turned(v: Vector, axis: Vector, angle: f64) -> Vector {
    v * angle.cos() + axis.cross(v) * angle.sin() + axis * (axis.dot(v) * (1.0 - angle.cos()))
}

/// The least rotation taking unit `from` to unit `to`, applied to `v`.
fn least_turn(v: Vector, from: Vector, to: Vector) -> Vector {
    let axis = from.cross(to);
    let s = axis.magnitude();
    if s <= 1e-15 {
        return v;
    }
    turned(v, axis / s, s.atan2(from.dot(to)))
}

#[test]
fn a_helix_turns_against_its_frenet_frame_at_minus_its_torsion() {
    // Radius 10, pitch 5, three turns: torsion c/(r²+c²) with c the rise
    // per radian, so the frame turns against the Frenet frame by minus the
    // torsion times the length, -3·2π·c/√(r²+c²) in all.
    let (r, pitch, turns) = (10.0, 5.0, 3.0);
    let c = pitch / TAU;
    let mut model = Model::new();
    let helix = edge_of(
        &mut model,
        HelixCurve::new(Frame::WORLD, r, pitch, turns)
            .unwrap()
            .into(),
    );
    let total = -turns * TAU * c / r.hypot(c);
    let frenet_normal = |t: f64| Vector::new(-t.cos(), -t.sin(), 0.0);
    let at = shares(240);
    let frames =
        rotation_minimizing_frames(&model, &helix, &at, Some(frenet_normal(0.0)), T).unwrap();
    assert!(frames.closing_twist.is_none());
    assert_orthonormal(&frames.frames);
    let (mut worst_twist, mut worst_tangent, mut worst_point) = (0.0_f64, 0.0_f64, 0.0_f64);
    for f in &frames.frames {
        // Length runs evenly with the turn angle on a helix.
        let t = f.share * turns * TAU;
        let on = Point::new(r * t.cos(), r * t.sin(), c * t);
        let tangent = Vector::new(-r * t.sin(), r * t.cos(), c) / r.hypot(c);
        worst_point = worst_point.max(f.at.distance(on));
        worst_tangent = worst_tangent.max((f.tangent.vector() - tangent).magnitude());
        let twist = angle_about(frenet_normal(t), f.normal.vector(), tangent);
        worst_twist = worst_twist.max(wrapped(twist - total * f.share).abs());
    }
    eprintln!(
        "helix: twist off {worst_twist:e} rad (end {total}), tangent off {worst_tangent:e}, \
         point off {worst_point:e}"
    );
    assert!(worst_twist <= 1e-6, "twist off by {worst_twist}");
    assert!(worst_tangent <= 1e-9, "tangent off by {worst_tangent}");
    assert!(worst_point <= 1e-6, "point off by {worst_point}");
}

#[test]
fn a_circles_frame_is_its_frenet_frame_turned_once_and_comes_home_untwisted() {
    let r = 7.0;
    let frame = Frame::new(
        Point::new(1.0, -2.0, 3.0),
        Direction::from_coords(1.0, 2.0, 2.0, T).unwrap(),
        Direction::X,
        T,
    )
    .unwrap();
    let mut model = Model::new();
    let circle = edge_of(
        &mut model,
        CircleCurve::new(Circle::new(frame, r, T).unwrap()).into(),
    );
    // A start normal leaning out of the plane: the angle to the Frenet
    // normal, which points at the centre, holds all the way round.
    let at = shares(97);
    let first = rotation_minimizing_frames(&model, &circle, &[0.0], None, T)
        .unwrap()
        .frames[0];
    let lean = first.normal.vector() * 0.6 + first.binormal() * 0.8;
    let frames = rotation_minimizing_frames(&model, &circle, &at, Some(lean), T).unwrap();
    assert_orthonormal(&frames.frames);
    let centre = frame.origin();
    let frenet = |f: &CurveFrame| (centre - f.at) / (centre - f.at).magnitude();
    let start = angle_about(
        frenet(&frames.frames[0]),
        frames.frames[0].normal.vector(),
        frames.frames[0].tangent.vector(),
    );
    let mut worst = 0.0_f64;
    for f in &frames.frames {
        let a = angle_about(frenet(f), f.normal.vector(), f.tangent.vector());
        worst = worst.max(wrapped(a - start).abs());
    }
    let twist = frames.closing_twist.unwrap();
    eprintln!("circle: angle to Frenet drifts {worst:e}, closing twist {twist:e}");
    assert!(
        worst <= 1e-12,
        "the angle to the Frenet frame drifts {worst}"
    );
    assert!(
        twist.abs() <= 1e-12,
        "a circle comes home turned by {twist}"
    );
}

#[test]
fn a_lines_frame_is_constant_and_the_start_normal_is_made_square() {
    let mut model = Model::new();
    let (a, b) = (Point::new(1.0, 2.0, 3.0), Point::new(4.0, -2.0, 15.0));
    let line = edge_of(&mut model, LineCurve::segment(a, b, T).unwrap().into());
    let asked = Vector::new(0.3, 1.0, -0.2);
    let frames = rotation_minimizing_frames(&model, &line, &shares(10), Some(asked), T).unwrap();
    assert_orthonormal(&frames.frames);
    let t = (b - a) / (b - a).magnitude();
    let square = asked - t * asked.dot(t);
    let square = square / square.magnitude();
    let worst = frames
        .frames
        .iter()
        .map(|f| (f.normal.vector() - square).magnitude())
        .fold(0.0, f64::max);
    eprintln!("line: normal moves {worst:e}");
    assert!(worst <= 1e-12, "the line's normal moves {worst}");
    // A normal along the line sets no frame.
    assert!(rotation_minimizing_frames(&model, &line, &[0.5], Some(t), T).is_err());
}

/// An S-curve in the plane through the origin square to `(1, 1, 1)`, with
/// one inflection where its bend changes side.
fn s_curve(model: &mut Model) -> (Shape, Vector) {
    let normal = Vector::new(1.0, 1.0, 1.0) / 3.0_f64.sqrt();
    let u = Vector::new(1.0, -1.0, 0.0) / 2.0_f64.sqrt();
    let v = normal.cross(u);
    let points: Vec<Point> = (0..=40)
        .map(|k| {
            let x = f64::from(k) / 40.0 * 30.0;
            Point::ORIGIN + u * x + v * (6.0 * (x / 30.0 * TAU).sin())
        })
        .collect();
    let curve = interpolate(&points, 3, Spacing::Chordal, T).unwrap();
    (edge_of(model, Curve::BSpline(curve)), normal)
}

#[test]
fn a_planar_curve_keeps_its_normal_through_an_inflection() {
    let mut model = Model::new();
    let (s, plane) = s_curve(&mut model);
    let at = shares(300);
    // Out of the plane: the normal is the plane's normal everywhere.
    let out = rotation_minimizing_frames(&model, &s, &at, Some(plane), T).unwrap();
    assert_orthonormal(&out.frames);
    let off = out
        .frames
        .iter()
        .map(|f| (f.normal.vector() - plane).magnitude())
        .fold(0.0, f64::max);
    // In the plane: it stays in the plane and never flips, where a Frenet
    // normal would turn over at the inflection.
    let inside = rotation_minimizing_frames(&model, &s, &at, Some(plane.cross(Vector::X)), T)
        .unwrap()
        .frames;
    let leaves = inside
        .iter()
        .map(|f| f.normal.vector().dot(plane).abs())
        .fold(0.0, f64::max);
    let least_step = inside
        .windows(2)
        .map(|w| w[0].normal.vector().dot(w[1].normal.vector()))
        .fold(1.0, f64::min);
    eprintln!("planar: normal off the plane's {off:e}, in-plane normal leaves by {leaves:e}");
    assert!(off <= 1e-12, "the out-of-plane normal moves {off}");
    assert!(
        leaves <= 1e-12,
        "the in-plane normal leaves the plane by {leaves}"
    );
    assert!(least_step > 0.99, "the normal flips: {least_step}");
}

#[test]
fn a_skew_polygons_frames_turn_least_at_its_corners_and_its_holonomy_is_their_product() {
    let corners = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(10.0, 0.0, 0.0),
        Point::new(10.0, 8.0, 3.0),
        Point::new(2.0, 9.0, 7.0),
    ];
    let mut model = Model::new();
    let wire = make_polygon(&mut model, &corners, true, T).unwrap().shape;
    let legs: Vec<Vector> = (0..4)
        .map(|i| {
            let d = corners[(i + 1) % 4] - corners[i];
            d / d.magnitude()
        })
        .collect();
    let lengths: Vec<f64> = (0..4)
        .map(|i| corners[i].distance(corners[(i + 1) % 4]))
        .collect();
    let total: f64 = lengths.iter().sum();
    // A station in the middle of each leg and one on each inner corner.
    let mut at = Vec::new();
    let mut run = 0.0;
    for l in &lengths {
        at.push(run / total);
        at.push((run + l / 2.0) / total);
        run += l;
    }
    at.push(1.0);
    let n0 = Vector::new(0.0, 0.6, 0.8);
    let frames = rotation_minimizing_frames(&model, &wire, &at, Some(n0), T).unwrap();
    assert_orthonormal(&frames.frames);
    // The oracle: the start normal turned by the least rotation at each
    // corner in turn, constant along each leg.
    let mut expected = n0;
    let mut worst = 0.0_f64;
    for (i, leg) in legs.iter().enumerate() {
        if i > 0 {
            expected = least_turn(expected, legs[i - 1], *leg);
        }
        for f in &frames.frames[2 * i..2 * i + 2] {
            worst = worst.max((f.normal.vector() - expected).magnitude());
            assert!((f.tangent.vector() - *leg).magnitude() <= 1e-12);
        }
    }
    // The end of the path takes the frame arriving.
    worst = worst.max((frames.frames[8].normal.vector() - expected).magnitude());
    let home = least_turn(expected, legs[3], legs[0]);
    let holonomy = angle_about(n0, home, legs[0]);
    let twist = frames.closing_twist.unwrap();
    eprintln!("polygon: frames off {worst:e}, holonomy {twist} vs {holonomy}");
    assert!(worst <= 1e-12, "the frames are off by {worst}");
    assert!(wrapped(twist - holonomy).abs() <= 1e-12);
    // Paid off, the end frame is the start frame.
    let closed = frames.closed_up(T).unwrap();
    assert!((closed[0].normal.vector() - n0).magnitude() <= 1e-12);
    let end = least_turn(closed[8].normal.vector(), legs[3], legs[0]);
    assert!((end - n0).magnitude() <= 1e-12, "{end:?}");
}

/// The frame equation `n' = -(n·t')t` integrated by fourth-order
/// Runge-Kutta in `steps` steps of `curve`'s parameter from `n` at `u0`
/// to `u1`.
fn integrated(curve: &Curve, (u0, u1): (f64, f64), n: Vector, steps: usize) -> Vector {
    let rate = |u: f64, n: Vector| -> Vector {
        let d = curve.derivatives_at(u, 2, T).unwrap();
        let speed = d[1].magnitude();
        let t = d[1] / speed;
        let dt = (d[2] - t * t.dot(d[2])) / speed;
        t * -n.dot(dt)
    };
    let h = (u1 - u0) / steps as f64;
    let mut n = n;
    for k in 0..steps {
        let u = u0 + h * k as f64;
        let k1 = rate(u, n);
        let k2 = rate(u + h / 2.0, n + k1 * (h / 2.0));
        let k3 = rate(u + h / 2.0, n + k2 * (h / 2.0));
        let k4 = rate(u + h, n + k3 * h);
        n += (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0);
        n /= n.magnitude();
    }
    n
}

/// A trefoil (the (2, 3) torus knot) interpolated through 120 points, its
/// first point repeated last: a closed skew curve whose ends meet at a
/// slight corner.
fn trefoil() -> Curve {
    let points: Vec<Point> = (0..=120)
        .map(|k| {
            let s = TAU * f64::from(k % 120) / 120.0;
            let rr = 8.0 + 3.0 * (3.0 * s).cos();
            Point::new(
                rr * (2.0 * s).cos(),
                rr * (2.0 * s).sin(),
                3.0 * (3.0 * s).sin(),
            )
        })
        .collect();
    Curve::BSpline(interpolate(&points, 3, Spacing::Centripetal, T).unwrap())
}

#[test]
fn a_skew_splines_frames_and_holonomy_match_the_frame_equation() {
    let curve = trefoil();
    let (lo, hi) = curve.domain();
    let mut model = Model::new();
    let knot = edge_of(&mut model, curve.clone());
    let at = shares(16);
    let frames = rotation_minimizing_frames(&model, &knot, &at, None, T).unwrap();
    assert_orthonormal(&frames.frames);
    let n0 = frames.frames[0].normal.vector();
    let mut worst = 0.0_f64;
    for f in &frames.frames {
        let expected = integrated(&curve, (lo, f.parameter), n0, 20_000);
        worst = worst.max((f.normal.vector() - expected).magnitude());
        let tangent = curve.d1_at(f.parameter, T).unwrap();
        assert!((f.tangent.vector() - tangent / tangent.magnitude()).magnitude() <= 1e-12);
    }
    // Home: the frame equation to the end, then the least turn across the
    // join where the spline's ends meet.
    let t_end = curve.d1_at(hi, T).unwrap();
    let t_start = curve.d1_at(lo, T).unwrap();
    let (t_end, t_start) = (t_end / t_end.magnitude(), t_start / t_start.magnitude());
    let end = integrated(&curve, (lo, hi), n0, 40_000);
    let home = least_turn(end, t_end, t_start);
    let holonomy = angle_about(n0, home, t_start);
    let twist = frames.closing_twist.unwrap();
    let off = wrapped(twist - holonomy).abs();
    eprintln!("trefoil: frames off {worst:e}, holonomy {twist} vs {holonomy} (off {off:e})");
    // The stations hold each step to about a sixty-fourth of a turn (149
    // of them here, the worst step 7°) and the spline is only C2 at its
    // knots, so double reflection lands within about 1e-5 of the
    // integrated frame (measured 9.1e-6).
    assert!(worst <= 2e-5, "the frames are off by {worst}");
    assert!(off <= 2e-5, "the holonomy is off by {off}");
}

#[test]
fn stations_outside_the_path_are_refused() {
    let mut model = Model::new();
    let line = edge_of(
        &mut model,
        LineCurve::segment(Point::ORIGIN, Point::new(1.0, 0.0, 0.0), T)
            .unwrap()
            .into(),
    );
    for bad in [-0.1, 1.1, f64::NAN] {
        assert!(rotation_minimizing_frames(&model, &line, &[bad], None, T).is_err());
    }
}
