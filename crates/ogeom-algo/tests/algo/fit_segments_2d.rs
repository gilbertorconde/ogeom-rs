//! Lines and arcs fitted to 2D polylines, against the outlines the points
//! were drawn from: a slot, a rectangle with rounded or filleted corners, a
//! hook and a shallow arc, with and without seeded Gaussian noise. Every
//! test checks the hard bound itself: each input point lies within the
//! tolerance of the piece whose range holds it, measured to the segment or
//! the arc, and consecutive pieces meet end to end.
#![allow(
    clippy::unwrap_used,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "test code"
)]

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use ogeom_algo::{Piece2, Segment2, Segments, fit_segments_2d};
use ogeom_core::OgeomError;
use ogeom_math::{Point2, Vector2};

/// Gaussian noise from a fixed seed.
struct Noise(u64);

impl Noise {
    fn uniform(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn gauss(&mut self, sigma: f64) -> f64 {
        let (u, v) = (self.uniform().max(1e-300), self.uniform());
        sigma * (-2.0 * u.ln()).sqrt() * (TAU * v).cos()
    }

    fn shake(&mut self, points: &mut [Point2], sigma: f64) {
        for p in points {
            *p = Point2::new(p.x + self.gauss(sigma), p.y + self.gauss(sigma));
        }
    }
}

/// Points along a line, every `step` or a little less, without its end.
fn line(points: &mut Vec<Point2>, from: Point2, to: Point2, step: f64) {
    let n = (from.distance(to) / step).ceil() as usize;
    points.extend((0..n).map(|k| from.lerp(to, k as f64 / n as f64)));
}

/// Points along an arc, every `step` or a little less, without its end.
fn arc(points: &mut Vec<Point2>, centre: Point2, radius: f64, start: f64, sweep: f64, step: f64) {
    let n = (radius * sweep.abs() / step).ceil() as usize;
    points.extend((0..n).map(|k| {
        let a = sweep.mul_add(k as f64 / n as f64, start);
        Point2::new(
            radius.mul_add(a.cos(), centre.x),
            radius.mul_add(a.sin(), centre.y),
        )
    }));
}

/// A rectangle from the origin to `(w, h)`, counter-clockwise, its corners
/// rounded with radius `r`, or sharp for `r = 0`.
fn rounded_rectangle(w: f64, h: f64, r: f64, step: f64) -> Vec<Point2> {
    let mut p = Vec::new();
    let corners = [
        (Point2::new(w - r, r), -FRAC_PI_2),
        (Point2::new(w - r, h - r), 0.0),
        (Point2::new(r, h - r), FRAC_PI_2),
        (Point2::new(r, r), PI),
    ];
    for k in 0..4 {
        let (c, a) = corners[k];
        let (next, b) = corners[(k + 1) % 4];
        if r > 0.0 {
            arc(&mut p, c, r, a, FRAC_PI_2, step);
        }
        let from = Point2::new(r.mul_add(b.cos(), c.x), r.mul_add(b.sin(), c.y));
        let to = Point2::new(r.mul_add(b.cos(), next.x), r.mul_add(b.sin(), next.y));
        line(&mut p, from, to, step);
    }
    // Start on the bottom side, away from any corner.
    let shift = p.len() - p.len() / 8;
    p.rotate_left(shift);
    p
}

fn options(tolerance: f64, closed: bool) -> Segments {
    Segments {
        tolerance,
        corners: true,
        axis_snap: None,
        min_bow: 0.0,
        closed,
    }
}

fn distance(p: Point2, piece: &Piece2) -> f64 {
    match *piece {
        Piece2::Line { from, to } => {
            let d = to - from;
            let t = ((p - from).dot(d) / d.dot(d)).clamp(0.0, 1.0);
            p.distance(from + d * t)
        }
        Piece2::Arc {
            centre,
            radius,
            start,
            sweep,
        } => {
            let d = p - centre;
            let mut a = d.y.atan2(d.x) - start;
            if sweep < 0.0 {
                a = -a;
            }
            if a.rem_euclid(TAU) <= sweep.abs() {
                (d.magnitude() - radius).abs()
            } else {
                p.distance(piece.start()).min(p.distance(piece.end()))
            }
        }
        Piece2::Circle { centre, radius, .. } => (p.distance(centre) - radius).abs(),
    }
}

/// The unit direction of travel at the start or the end of a piece.
fn heading(piece: &Piece2, at_end: bool) -> Vector2 {
    match *piece {
        Piece2::Line { from, to } => (to - from) / from.distance(to),
        Piece2::Arc { start, sweep, .. } => {
            let a = if at_end { start + sweep } else { start };
            Vector2::new(-a.sin(), a.cos()) * sweep.signum()
        }
        Piece2::Circle { .. } => unreachable!("not in a chain"),
    }
}

/// Every point within the tolerance of its piece, every point in some
/// piece's range, the ranges running on from one to the next, and the
/// pieces meeting end to end. The turn at each joint, in radians.
fn check(points: &[Point2], segments: &[Segment2], tolerance: f64, closed: bool) -> Vec<f64> {
    let n = points.len();
    let mut covered = vec![false; n];
    for s in segments {
        let (first, last) = s.range;
        let count = if last >= first {
            last - first
        } else {
            last + n - first
        };
        let mut sum = 0.0;
        for k in 0..=count {
            let i = (first + k) % n;
            let d = distance(points[i], &s.piece);
            assert!(d <= tolerance, "point {i} is {d} from {s:?}");
            sum += d * d;
            covered[i] = true;
        }
        assert!((s.sigma - (sum / (count + 1) as f64).sqrt()).abs() < 1e-12);
    }
    assert!(covered.iter().all(|c| *c), "a point is in no range");
    let m = segments.len();
    let joints = if closed { m } else { m - 1 };
    (0..joints)
        .map(|k| {
            let (a, b) = (&segments[k], &segments[(k + 1) % m]);
            assert_eq!(a.range.1, b.range.0, "{a:?} then {b:?}");
            let gap = a.piece.end().distance(b.piece.start());
            assert!(gap < 1e-9, "a gap of {gap} between {a:?} and {b:?}");
            let (u, v) = (heading(&a.piece, true), heading(&b.piece, false));
            u.cross(v).atan2(u.dot(v)).abs()
        })
        .collect()
}

fn lines(segments: &[Segment2]) -> usize {
    segments
        .iter()
        .filter(|s| matches!(s.piece, Piece2::Line { .. }))
        .count()
}

fn radii(segments: &[Segment2]) -> Vec<f64> {
    segments
        .iter()
        .filter_map(|s| match s.piece {
            Piece2::Arc { radius, .. } => Some(radius),
            _ => None,
        })
        .collect()
}

/// A slot of two straights 20 long and two half circles of radius 5,
/// scanned with noise of σ = 0.02: two lines and two arcs, tangent, the
/// radii within 1 %.
#[test]
fn a_noisy_slot_gives_two_lines_and_two_tangent_arcs() {
    let mut p = Vec::new();
    line(&mut p, Point2::new(0.0, -5.0), Point2::new(20.0, -5.0), 0.1);
    arc(&mut p, Point2::new(20.0, 0.0), 5.0, -FRAC_PI_2, PI, 0.1);
    line(&mut p, Point2::new(20.0, 5.0), Point2::new(0.0, 5.0), 0.1);
    arc(&mut p, Point2::new(0.0, 0.0), 5.0, FRAC_PI_2, PI, 0.1);
    Noise(0x9e37_79b9_7f4a_7c15).shake(&mut p, 0.02);
    let tolerance = 0.1;
    let s = fit_segments_2d(&p, &options(tolerance, true)).unwrap();
    assert_eq!(s.len(), 4, "{s:#?}");
    assert_eq!(lines(&s), 2, "{s:#?}");
    for r in radii(&s) {
        assert!((r - 5.0).abs() < 0.05, "radius {r}");
    }
    for turn in check(&p, &s, tolerance, true) {
        assert!(turn < 1e-9, "a joint turns by {turn}");
    }
}

/// A rectangle 40 by 20 whose corners a scan has rounded off with radius
/// 1, noise σ = 0.02, fitted within 0.4: the rounding lies within 0.3 of
/// the corner's two sides, so it gives way to four lines meeting at the
/// sharp corners.
#[test]
fn a_rectangle_with_corner_blobs_gives_four_lines_meeting_at_its_corners() {
    let mut p = rounded_rectangle(40.0, 20.0, 1.0, 0.1);
    Noise(0x2545_f491_4f6c_dd1d).shake(&mut p, 0.02);
    let tolerance = 0.4;
    let s = fit_segments_2d(&p, &options(tolerance, true)).unwrap();
    assert_eq!(s.len(), 4, "{s:#?}");
    assert_eq!(lines(&s), 4, "{s:#?}");
    let turns = check(&p, &s, tolerance, true);
    for turn in turns {
        assert!((turn - FRAC_PI_2).abs() < 0.01, "a corner turns by {turn}");
    }
    let corners = [(0.0, 0.0), (40.0, 0.0), (40.0, 20.0), (0.0, 20.0)];
    for seg in &s {
        let at = seg.piece.end();
        let nearest = corners
            .iter()
            .map(|&(x, y)| at.distance(Point2::new(x, y)))
            .fold(f64::INFINITY, f64::min);
        assert!(nearest < 0.05, "a corner at {at:?}");
    }
    // Without `corners` the lines meet on the rounding, where the points
    // are, a chamfer the tolerance allows rather than the corner.
    let chamfered = fit_segments_2d(
        &p,
        &Segments {
            corners: false,
            ..options(tolerance, true)
        },
    )
    .unwrap();
    assert_eq!(lines(&chamfered), 4, "{chamfered:#?}");
    check(&p, &chamfered, tolerance, true);
    for seg in &chamfered {
        let at = seg.piece.end();
        let corner = corners
            .iter()
            .map(|&(x, y)| at.distance(Point2::new(x, y)))
            .fold(f64::INFINITY, f64::min);
        let input = p
            .iter()
            .map(|q| at.distance(*q))
            .fold(f64::INFINITY, f64::min);
        assert!(corner > 0.1 && input < tolerance, "a joint at {at:?}");
    }
}

/// The same rectangle with fillets of radius 1 a tolerance of 0.1 cannot
/// hide: each fillet stays an arc of radius 1, tangent to its sides.
#[test]
fn a_fillet_wider_than_the_tolerance_is_kept() {
    let mut p = rounded_rectangle(40.0, 20.0, 1.0, 0.05);
    Noise(0x6a09_e667_f3bc_c909).shake(&mut p, 0.01);
    let tolerance = 0.1;
    let s = fit_segments_2d(&p, &options(tolerance, true)).unwrap();
    assert_eq!(s.len(), 8, "{s:#?}");
    assert_eq!(lines(&s), 4, "{s:#?}");
    for r in radii(&s) {
        assert!((r - 1.0).abs() < 0.02, "radius {r}");
    }
    for turn in check(&p, &s, tolerance, true) {
        assert!(turn < 1e-9, "a joint turns by {turn}");
    }
}

/// An open hook, a line, a quarter circle and a line, noise free: three
/// pieces from its first point to its last, tangent where they meet.
#[test]
fn an_open_hook_runs_from_its_first_point_to_its_last() {
    let mut p = Vec::new();
    line(&mut p, Point2::new(0.0, 0.0), Point2::new(10.0, 0.0), 0.1);
    arc(
        &mut p,
        Point2::new(10.0, 3.0),
        3.0,
        -FRAC_PI_2,
        FRAC_PI_2,
        0.1,
    );
    line(&mut p, Point2::new(13.0, 3.0), Point2::new(13.0, 10.0), 0.1);
    p.push(Point2::new(13.0, 10.0));
    let tolerance = 0.01;
    let s = fit_segments_2d(&p, &options(tolerance, false)).unwrap();
    assert_eq!(s.len(), 3, "{s:#?}");
    assert_eq!(lines(&s), 2, "{s:#?}");
    assert_eq!(s[0].range.0, 0);
    assert_eq!(s[2].range.1, p.len() - 1);
    assert!(s[0].piece.start().distance(p[0]) < 1e-9);
    assert!(s[2].piece.end().distance(p[p.len() - 1]) < 1e-9);
    assert!((radii(&s)[0] - 3.0).abs() < 1e-3);
    for turn in check(&p, &s, tolerance, false) {
        assert!(turn < 1e-9, "a joint turns by {turn}");
    }
}

/// A rectangle turned 0.05° off the axes: with `axis_snap` its sides turn
/// back onto them, and without it they stay as drawn.
#[test]
fn axis_snap_turns_near_axis_lines_onto_the_axes() {
    let turn = 0.05f64.to_radians();
    let (sin, cos) = turn.sin_cos();
    let p: Vec<Point2> = rounded_rectangle(40.0, 20.0, 0.0, 0.5)
        .into_iter()
        .map(|q| {
            let (x, y) = (q.x - 20.0, q.y - 10.0);
            Point2::new(x.mul_add(cos, -y * sin), x.mul_add(sin, y * cos))
        })
        .collect();
    let tolerance = 0.1;
    let drawn = fit_segments_2d(&p, &options(tolerance, true)).unwrap();
    let snapped = fit_segments_2d(
        &p,
        &Segments {
            axis_snap: Some(1f64.to_radians()),
            ..options(tolerance, true)
        },
    )
    .unwrap();
    let off_axis = |s: &Segment2| {
        let (a, b) = (s.piece.start(), s.piece.end());
        (a.x - b.x).abs().min((a.y - b.y).abs())
    };
    assert_eq!(drawn.len(), 4, "{drawn:#?}");
    assert_eq!(snapped.len(), 4, "{snapped:#?}");
    check(&p, &drawn, tolerance, true);
    check(&p, &snapped, tolerance, true);
    assert!(drawn.iter().all(|s| off_axis(s) > 1e-3), "{drawn:#?}");
    assert!(snapped.iter().all(|s| off_axis(s) < 1e-12), "{snapped:#?}");
}

/// An arc of radius 100 over a chord of 20 bows 0.5 from it: one arc
/// within 0.1, and with `min_bow` of 1 the fewest lines, two.
#[test]
fn an_arc_bowing_less_than_min_bow_is_taken_as_lines() {
    let mut p = Vec::new();
    let half = (10.0f64 / 100.0).asin();
    arc(
        &mut p,
        Point2::new(0.0, -100.0),
        100.0,
        FRAC_PI_2 + half,
        -2.0 * half,
        0.1,
    );
    let tolerance = 0.1;
    let one = fit_segments_2d(&p, &options(tolerance, false)).unwrap();
    assert_eq!(one.len(), 1, "{one:#?}");
    assert_eq!(radii(&one).len(), 1);
    check(&p, &one, tolerance, false);
    let straight = fit_segments_2d(
        &p,
        &Segments {
            min_bow: 1.0,
            ..options(tolerance, false)
        },
    )
    .unwrap();
    assert_eq!(straight.len(), 2, "{straight:#?}");
    assert_eq!(lines(&straight), 2);
    check(&p, &straight, tolerance, false);
}

/// A closed loop one circle holds is that circle.
#[test]
fn a_closed_circle_is_one_circle() {
    let mut p = Vec::new();
    arc(&mut p, Point2::new(3.0, 4.0), 2.0, 0.0, TAU, 0.05);
    let s = fit_segments_2d(&p, &options(1e-3, true)).unwrap();
    assert_eq!(s.len(), 1);
    assert!(matches!(
        s[0].piece,
        Piece2::Circle {
            counter_clockwise: true,
            ..
        }
    ));
    assert_eq!(s[0].range, (0, 0));
}

#[test]
fn bad_input_is_refused() {
    let p = [Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)];
    let refused = |p: &[Point2], o: &Segments| fit_segments_2d(p, o).unwrap_err();
    assert!(matches!(
        refused(&p, &options(0.0, false)),
        OgeomError::Construction(_)
    ));
    assert!(matches!(
        refused(&p, &options(0.1, true)),
        OgeomError::Construction(_)
    ));
    assert!(matches!(
        refused(&p[..1], &options(0.1, false)),
        OgeomError::Construction(_)
    ));
    assert!(matches!(
        refused(&[p[0], Point2::new(f64::NAN, 0.0)], &options(0.1, false)),
        OgeomError::Numeric(_)
    ));
    assert!(matches!(
        refused(
            &p,
            &Segments {
                axis_snap: Some(1.0),
                ..options(0.1, false)
            }
        ),
        OgeomError::Construction(_)
    ));
    let one = fit_segments_2d(&p, &options(0.1, false)).unwrap();
    assert_eq!(one.len(), 1);
    assert_eq!(one[0].range, (0, 1));
}

/// Random chains of lines and arcs, sharp and smooth turns, noise from
/// none to the tolerance, open and closed: every fit holds every point
/// within the tolerance, the fallbacks included.
#[test]
fn random_chains_hold_the_bound() {
    let mut noise = Noise(0x0123_4567_89ab_cdef);
    for case in 0..300 {
        let mut p = Vec::new();
        let (mut at, mut heading) = (Point2::new(0.0, 0.0), 0.0f64);
        for _ in 0..1 + case % 7 {
            let step = 0.05 + 0.3 * noise.uniform();
            if noise.uniform() < 0.5 {
                let length = 0.5 + 10.0 * noise.uniform();
                let to = at + Vector2::new(heading.cos(), heading.sin()) * length;
                line(&mut p, at, to, step);
                at = to;
            } else {
                let radius = 0.2 + 5.0 * noise.uniform();
                let sweep = (noise.uniform() - 0.5) * 2.0 * PI;
                let side = sweep.signum() * FRAC_PI_2;
                let centre =
                    at + Vector2::new((heading + side).cos(), (heading + side).sin()) * radius;
                let start = (at.y - centre.y).atan2(at.x - centre.x);
                arc(&mut p, centre, radius, start, sweep, step);
                heading += sweep;
                at = Point2::new(
                    radius.mul_add((start + sweep).cos(), centre.x),
                    radius.mul_add((start + sweep).sin(), centre.y),
                );
            }
            if noise.uniform() < 0.5 {
                heading += (noise.uniform() - 0.5) * 3.0;
            }
        }
        p.push(at);
        let tolerance = 0.01 + 0.2 * noise.uniform();
        let sigma = tolerance * noise.uniform() / 3.0;
        noise.shake(&mut p, sigma);
        let closed = case % 2 == 1 && p.len() > 3;
        let o = Segments {
            corners: case % 3 != 0,
            axis_snap: (case % 5 == 0).then_some(0.05),
            min_bow: if case % 4 == 0 { 2.0 * tolerance } else { 0.0 },
            ..options(tolerance, closed)
        };
        let s = fit_segments_2d(&p, &o).unwrap();
        if !matches!(s[0].piece, Piece2::Circle { .. }) {
            check(&p, &s, tolerance, closed);
        }
    }
}

/// A "D" of radius 10, every `step`: the half circle over the top, then
/// the diameter back along y = 0, started `shift` points on.
fn d_shape(step: f64, shift: usize) -> Vec<Point2> {
    let mut p = Vec::new();
    arc(&mut p, Point2::new(0.0, 0.0), 10.0, 0.0, PI, step);
    line(
        &mut p,
        Point2::new(-10.0, 0.0),
        Point2::new(10.0, 0.0),
        step,
    );
    p.rotate_left(shift);
    p
}

/// The D's arc and line: the arc's radius within `radius` of 10 and its
/// sweep π, the line 20 long along y = 0, both ends of each within
/// `within` of the corners (10, 0) and (-10, 0).
fn assert_d(s: &[Segment2], radius: f64, within: f64) {
    assert_eq!(s.len(), 2, "{s:#?}");
    assert_eq!(lines(s), 1, "{s:#?}");
    for seg in s {
        match seg.piece {
            Piece2::Arc {
                radius: r, sweep, ..
            } => {
                assert!((r - 10.0).abs() <= radius, "radius {r}");
                assert!((sweep.abs() - PI).abs() <= within / 10.0, "sweep {sweep}");
            }
            Piece2::Line { from, to } => {
                assert!(from.y.abs().max(to.y.abs()) <= within, "{seg:?}");
                assert!((from.distance(to) - 20.0).abs() <= 2.0 * within, "{seg:?}");
            }
            Piece2::Circle { .. } => panic!("{seg:?}"),
        }
        for corner in [Point2::new(10.0, 0.0), Point2::new(-10.0, 0.0)] {
            let near = [seg.piece.start(), seg.piece.end()]
                .iter()
                .map(|e| e.distance(corner))
                .fold(f64::INFINITY, f64::min);
            assert!(near <= within, "no end of {seg:?} at {corner:?}");
        }
    }
}

/// A closed "D", a half circle of radius 10 closed along its diameter,
/// gives the arc and the line its open chain does, meeting at its
/// corners, with and without `corners`, with `axis_snap` and at a looser
/// tolerance. A loop of two runs has its one neighbour on both sides.
#[test]
fn a_closed_d_gives_one_arc_and_one_line() {
    let on = |t: f64| Point2::new(10.0 * t.cos(), 10.0 * t.sin());
    let mut d: Vec<Point2> = (0..=40).map(|k| on(k as f64 * PI / 40.0)).collect();
    d.extend((1..20).map(|k| Point2::new(-10.0 + k as f64, 0.0)));
    let o = |closed| Segments {
        tolerance: 0.05,
        corners: true,
        axis_snap: None,
        min_bow: 0.1,
        closed,
    };
    let open = fit_segments_2d(&d, &o(false)).unwrap();
    assert_eq!(open.len(), 2, "{open:#?}");
    for options in [
        o(true),
        Segments {
            corners: false,
            ..o(true)
        },
        Segments {
            axis_snap: Some(0.01),
            ..o(true)
        },
        Segments {
            tolerance: 0.2,
            ..o(true)
        },
    ] {
        let s = fit_segments_2d(&d, &options).unwrap();
        check(&d, &s, options.tolerance, true);
        assert_d(&s, 1e-9, 1e-9);
    }
}

/// D shapes started on the arc and mid-line, noise free and with noise of
/// σ = 0.02, with and without `corners`: the radius within 1e-9 noise free
/// and within σ with noise, the corners within the tolerance.
#[test]
fn closed_d_shapes_meet_at_their_corners() {
    let tolerance = 0.1;
    for (sigma, shift) in [(0.0, 0), (0.0, 400), (0.02, 0), (0.02, 400)] {
        let mut p = d_shape(0.05, shift);
        Noise(0x9e37_79b9_7f4a_7c15).shake(&mut p, sigma);
        for corners in [true, false] {
            let s = fit_segments_2d(
                &p,
                &Segments {
                    corners,
                    ..options(tolerance, true)
                },
            )
            .unwrap();
            check(&p, &s, tolerance, true);
            if sigma == 0.0 {
                assert_d(&s, 1e-9, 1e-9);
            } else {
                assert_d(&s, sigma, tolerance);
            }
        }
    }
}

/// A lens, two arcs of radius 10 about (0, ±6) crossing at (±8, 0): two
/// arcs, each sweeping 2 asin(0.8), meeting at the crossings.
#[test]
fn a_closed_lens_gives_two_arcs_meeting_at_its_corners() {
    let half = 0.8f64.asin();
    let mut p = Vec::new();
    arc(
        &mut p,
        Point2::new(0.0, -6.0),
        10.0,
        FRAC_PI_2 - half,
        2.0 * half,
        0.05,
    );
    arc(
        &mut p,
        Point2::new(0.0, 6.0),
        10.0,
        -FRAC_PI_2 - half,
        2.0 * half,
        0.05,
    );
    let tolerance = 0.01;
    let s = fit_segments_2d(&p, &options(tolerance, true)).unwrap();
    assert_eq!(s.len(), 2, "{s:#?}");
    check(&p, &s, tolerance, true);
    for seg in &s {
        let Piece2::Arc { radius, sweep, .. } = seg.piece else {
            panic!("{seg:?}");
        };
        assert!((radius - 10.0).abs() < 1e-9, "radius {radius}");
        assert!((sweep - 2.0 * half).abs() < 1e-9, "sweep {sweep}");
        let at = seg.piece.end();
        assert!((at.x.abs() - 8.0).abs().max(at.y.abs()) < 1e-9, "{at:?}");
    }
}

/// Closed slots and rounded rectangles, loops of more than two runs,
/// noise free and with noise of σ = 0.02, with and without `corners`:
/// their lines and their tangent arcs at the drawn lengths and radii.
#[test]
fn closed_slots_and_rounded_rectangles_keep_their_pieces() {
    let tolerance = 0.1;
    for sigma in [0.0, 0.02] {
        let mut slot = Vec::new();
        line(
            &mut slot,
            Point2::new(0.0, -5.0),
            Point2::new(20.0, -5.0),
            0.1,
        );
        arc(&mut slot, Point2::new(20.0, 0.0), 5.0, -FRAC_PI_2, PI, 0.1);
        line(
            &mut slot,
            Point2::new(20.0, 5.0),
            Point2::new(0.0, 5.0),
            0.1,
        );
        arc(&mut slot, Point2::new(0.0, 0.0), 5.0, FRAC_PI_2, PI, 0.1);
        let mut rect = rounded_rectangle(40.0, 20.0, 3.0, 0.1);
        let mut noise = Noise(0x2545_f491_4f6c_dd1d);
        noise.shake(&mut slot, sigma);
        noise.shake(&mut rect, sigma);
        let slack = if sigma == 0.0 { 1e-9 } else { sigma };
        for corners in [true, false] {
            let o = Segments {
                corners,
                ..options(tolerance, true)
            };
            for (p, count, radius, sides) in [
                (&slot, 4, 5.0, &[20.0][..]),
                (&rect, 8, 3.0, &[34.0, 14.0][..]),
            ] {
                let s = fit_segments_2d(p, &o).unwrap();
                assert_eq!(s.len(), count, "{s:#?}");
                assert_eq!(lines(&s), count / 2, "{s:#?}");
                for r in radii(&s) {
                    assert!((r - radius).abs() <= slack, "radius {r}");
                }
                for seg in &s {
                    if let Piece2::Line { from, to } = seg.piece {
                        let length = from.distance(to);
                        let nearest = sides
                            .iter()
                            .map(|l| (length - l).abs())
                            .fold(f64::INFINITY, f64::min);
                        assert!(nearest <= 2.0 * tolerance, "length {length}");
                    }
                }
                for turn in check(p, &s, tolerance, true) {
                    assert!(turn < 1e-9, "a joint turns by {turn}");
                }
            }
        }
    }
}
