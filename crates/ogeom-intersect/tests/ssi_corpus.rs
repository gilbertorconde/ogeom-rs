//! A small surface/surface intersection corpus, each case scored against
//! its exact section: the two-sided Hausdorff distance between the traced
//! polylines and the exact curves, and how many section ends lie on neither
//! a touch nor the edge of a chart. Each case holds both distances within a
//! few chords and the free ends within a stated count.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::{
    BSplineSurface, CylinderSurface, PlaneSurface, SphereSurface, Surface as _, SurfaceGeometry,
};
use ogeom_intersect::{Marching, Traced, branches};
use ogeom_math::{ControlGrid, Cylinder, Direction, Frame, KnotVector, Plane, Point, Sphere};

const T: Tolerances = Tolerances::millimetres();

/// A bicubic patch `z = Σ c·xᵖ·yᵠ` over a box: each monomial's control
/// height the product of its blossoms.
fn polynomial(
    terms: &[(i32, i32, f64)],
    (x0, x1): (f64, f64),
    (y0, y1): (f64, f64),
) -> SurfaceGeometry {
    let blossom = |power: i32, at: [f64; 3]| match power {
        0 => 1.0,
        1 => (at[0] + at[1] + at[2]) / 3.0,
        2 => (at[0] * at[1] + at[1] * at[2] + at[0] * at[2]) / 3.0,
        _ => at[0] * at[1] * at[2],
    };
    let args = |i: usize, lo: f64, hi: f64| {
        let mut at = [lo; 3];
        at[..i].fill(hi);
        at
    };
    let mut net = Vec::new();
    for i in 0..=3_usize {
        for j in 0..=3_usize {
            let (ax, ay) = (args(i, x0, x1), args(j, y0, y1));
            let z: f64 = terms
                .iter()
                .map(|&(p, q, c)| c * blossom(p, ax) * blossom(q, ay))
                .sum();
            let (fi, fj) = (
                f64::from(u8::try_from(i).unwrap()),
                f64::from(u8::try_from(j).unwrap()),
            );
            net.push(Point::new(
                x0 + (x1 - x0) * fi / 3.0,
                y0 + (y1 - y0) * fj / 3.0,
                z,
            ));
        }
    }
    let cubic = || KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 3).unwrap();
    BSplineSurface::new(cubic(), cubic(), &ControlGrid::new(net, 4, 4).unwrap(), T)
        .unwrap()
        .into()
}

fn ground() -> SurfaceGeometry {
    PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-2.0, 2.0),
        (-2.0, 2.0),
    )
    .unwrap()
    .into()
}

/// The exact section as polylines dense enough to stand within a chord of
/// it, and the points it may end at.
struct Truth {
    curves: Vec<Vec<Point>>,
    ends: Vec<Point>,
}

/// Samples of the segment from `p` to `q`.
fn segment(p: Point, q: Point, count: u32) -> Vec<Point> {
    (0..=count)
        .map(|k| p + (q - p) * (f64::from(k) / f64::from(count)))
        .collect()
}

/// Samples of a closed parametric curve over `0..=turn`.
fn closed(count: u32, turn: f64, at: impl Fn(f64) -> Point) -> Vec<Point> {
    (0..=count)
        .map(|k| at(turn * f64::from(k) / f64::from(count)))
        .collect()
}

/// The distance from a point to a segment.
fn to_segment(p: Point, from: Point, to: Point) -> f64 {
    let along = to - from;
    let length = along.square_magnitude();
    let t = if length > 0.0 {
        ((p - from).dot(along) / length).clamp(0.0, 1.0)
    } else {
        0.0
    };
    p.distance(from + along * t)
}

/// The distance from a point to a polyline, looked for first among the
/// segments beside `last`, the one the point before was nearest, and over
/// the whole polyline where those stand further off than a tenth of
/// `within`: never less than the true distance.
fn to_polyline_from(p: Point, polyline: &[Point], last: &mut usize, within: f64) -> f64 {
    let segments = polyline.len().saturating_sub(1);
    if segments == 0 {
        return f64::INFINITY;
    }
    let best = |range: core::ops::Range<usize>| {
        range
            .map(|k| (to_segment(p, polyline[k], polyline[k + 1]), k))
            .fold((f64::INFINITY, 0), |x, y| if y.0 < x.0 { y } else { x })
    };
    let from = last.saturating_sub(64).min(segments);
    let (d, k) = best(from..(*last + 64).min(segments));
    let (d, k) = if d <= within {
        (d, k)
    } else {
        best(0..segments)
    };
    *last = k;
    d
}

/// Both directed Hausdorff distances, traced to exact and exact to traced,
/// and the ends lying on neither the truth's ends nor a chart's edge.
fn score(
    found: &[Traced],
    truth: &Truth,
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    within: f64,
) -> (f64, f64, usize) {
    let mut lasts = vec![0; truth.curves.len()];
    let traced_to_truth = found
        .iter()
        .flat_map(|s| s.points.iter())
        .map(|p| {
            truth
                .curves
                .iter()
                .zip(lasts.iter_mut())
                .map(|(c, last)| to_polyline_from(*p, c, last, within))
                .fold(f64::INFINITY, f64::min)
        })
        .fold(0.0, f64::max);
    // Every exact sample's distance to the nearest traced chord, over all
    // of them: a gap anywhere shows.
    let truth_to_traced = truth
        .curves
        .iter()
        .flat_map(|c| c.iter().step_by(10))
        .map(|q| {
            found
                .iter()
                .flat_map(|s| s.points.windows(2))
                .map(|w| to_segment(*q, w[0], w[1]))
                .fold(f64::INFINITY, f64::min)
        })
        .fold(0.0, f64::max);
    let on_edge = |surface: &SurfaceGeometry, (u, v): (f64, f64)| {
        let ((u0, u1), (v0, v1)) = surface.domain();
        let near = |x: f64, lo: f64, hi: f64| {
            let span = (hi - lo).abs().max(1.0) * 1e-7;
            (x - lo).abs() <= span || (x - hi).abs() <= span
        };
        near(u, u0, u1) || near(v, v0, v1)
    };
    let mut free = 0;
    for s in found {
        if s.closed() {
            continue;
        }
        let n = s.points.len();
        for k in [0, n - 1] {
            let at_end = truth
                .ends
                .iter()
                .any(|e| e.distance(s.points[k]) <= T.confusion());
            if !at_end && !on_edge(a, s.on_a[k]) && !on_edge(b, s.on_b[k]) {
                free += 1;
            }
        }
    }
    (traced_to_truth, truth_to_traced, free)
}

fn check(
    name: &str,
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    chord: f64,
    truth: &Truth,
    (within, free_ends): (f64, usize),
) -> Vec<Traced> {
    let marching = Marching {
        chord,
        ..Marching::default()
    };
    let found = branches(a, b, marching, T).unwrap();
    let (out, missed, free) = score(&found, truth, (a, b), within);
    println!("{name}: traced to exact {out:.2e}, exact to traced {missed:.2e}, free ends {free}");
    assert!(
        out <= within,
        "{name}: a section stands {out:e} off the exact one"
    );
    assert!(
        missed <= within,
        "{name}: the exact section is missed by {missed:e}"
    );
    assert!(
        free <= free_ends,
        "{name}: {free} ends lie on neither a touch nor an edge"
    );
    found
}

/// `z = K·(x³ − s²x)(y³ − s²y)` over a box, cut by `z = 0`: three lines
/// each way crossing at nine touches. At `s = 0.1` a straight stride of the
/// march spans two crossings. Where the stretches between crossings are
/// joined, each is traced once: the length traced exceeds the exact one by
/// no more than a hundred chords. At `K = 1`, `s = 0.03` the plane and the
/// patch part by less than the confusion distance well away from the lines,
/// and only coverage is asked: every line is traced, wherever between the
/// crossings it is placed. Where the patch crosses the plane at a slope
/// `t`, a section is placed only to within the confusion distance over `t`:
/// near a crossing `t` is `4K·s⁴·r` at `r` from it, about 1e-4 a few chords
/// out at `K = 10`, `s = 0.1`, and between the crossings at `s = 0.03` it
/// stays under `K·s⁵`.
#[test]
fn a_grid_of_lines_crossing_at_nine_touches() {
    for (s, k, chord, within, placed) in [
        (0.3, 100.0, 1e-4, 1e-3, true),
        (0.1, 100.0, 1e-4, 1e-3, true),
        (0.1, 100.0, 1e-7, 1e-6, true),
        (0.1, 10.0, 1e-6, 2e-4, true),
        (0.03, 100.0, 1e-7, 1e-4, true),
        (0.03, 1.0, 1e-7, 1e-4, false),
    ] {
        grid_of_lines(s, k, chord, within, placed);
    }
}

fn grid_of_lines(s: f64, k: f64, chord: f64, within: f64, placed: bool) {
    let ((x0, x1), (y0, y1)) = ((-0.5, 0.6), (-0.45, 0.55));
    let s2 = s * s;
    let patch = polynomial(
        &[
            (3, 3, k),
            (3, 1, -s2 * k),
            (1, 3, -s2 * k),
            (1, 1, s2 * s2 * k),
        ],
        (x0, x1),
        (y0, y1),
    );
    let mut curves = Vec::new();
    let mut ends = Vec::new();
    for c in [-s, 0.0, s] {
        curves.push(segment(
            Point::new(c, y0, 0.0),
            Point::new(c, y1, 0.0),
            2000,
        ));
        curves.push(segment(
            Point::new(x0, c, 0.0),
            Point::new(x1, c, 0.0),
            2000,
        ));
        for d in [-s, 0.0, s] {
            ends.push(Point::new(c, d, 0.0));
        }
    }
    let truth = Truth { curves, ends };
    let name = format!("line grid s={s} K={k} chord={chord:e}");
    let plane = ground();
    if !placed {
        let marching = Marching {
            chord,
            ..Marching::default()
        };
        let found = branches(&patch, &plane, marching, T).unwrap();
        let (_, missed, _) = score(&found, &truth, (&patch, &plane), within);
        assert!(
            missed <= within,
            "{name}: the exact section is missed by {missed:e}"
        );
        return;
    }
    let found = check(&name, (&patch, &plane), chord, &truth, (within, 0));
    let traced: f64 = found
        .iter()
        .flat_map(|f| f.points.windows(2).map(|w| w[0].distance(w[1])))
        .sum();
    let exact = 3.0 * (x1 - x0) + 3.0 * (y1 - y0);
    assert!(
        traced - exact <= 100.0 * chord,
        "{name}: {:e} of the section is traced twice",
        traced - exact
    );
}

/// `z = x·y·(y − d)` over a box, a patch linear in `x` and quadratic in
/// `y`, cut by `z = 0`: the lines `x = 0`, `y = 0` and `y = d`, crossing
/// at two touches closer together than the march strides along `x = 0`.
/// With `d` a few thousandths, the patch stays within `d²` of the plane
/// across the strip between `y = 0` and `y = d`, and a march found one
/// crossing of the two and lost the other line whole.
#[test]
fn a_saddle_crossing_a_plane_along_three_lines() {
    for (d, (x0, x1), (y0, y1), chord) in [
        (0.05, (-1.0, 1.0), (-1.1, 0.9), 1e-4),
        (0.002, (-1.0, 1.3), (-0.6, 1.2), 1e-4),
        (0.001, (-1.0, 1.0), (-1.1, 1.2), 1e-5),
        (0.003, (-0.7, 1.0), (-1.1, 0.9), 1e-5),
    ] {
        saddle_along_three_lines(d, (x0, x1), (y0, y1), chord);
    }
}

fn saddle_along_three_lines(d: f64, (x0, x1): (f64, f64), (y0, y1): (f64, f64), chord: f64) {
    let blossom = |s: f64, t: f64| s * t - d * (s + t) / 2.0;
    let mut net = Vec::new();
    for x in [x0, x1] {
        for (y, (s, t)) in [
            (y0, (y0, y0)),
            (f64::midpoint(y0, y1), (y0, y1)),
            (y1, (y1, y1)),
        ] {
            net.push(Point::new(x, y, x * blossom(s, t)));
        }
    }
    let patch: SurfaceGeometry = BSplineSurface::new(
        KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap(),
        KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap(),
        &ControlGrid::new(net, 2, 3).unwrap(),
        T,
    )
    .unwrap()
    .into();
    let plane: SurfaceGeometry = PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap()
    .into();
    let truth = Truth {
        curves: vec![
            segment(Point::new(0.0, y0, 0.0), Point::new(0.0, y1, 0.0), 4000),
            segment(Point::new(x0, 0.0, 0.0), Point::new(x1, 0.0, 0.0), 4000),
            segment(Point::new(x0, d, 0.0), Point::new(x1, d, 0.0), 4000),
        ],
        ends: vec![Point::ORIGIN, Point::new(0.0, d, 0.0)],
    };
    check(
        &format!("saddle d={d} chord={chord:e}"),
        (&patch, &plane),
        chord,
        &truth,
        (10.0 * chord, 0),
    );
}

/// `z = x·y·((x − a)² + (y − a)² − ρ²)` cut by `z = 0`: the axes, crossing
/// at the origin, and a circle passing between them without meeting them.
#[test]
fn two_axes_and_a_circle_beside_them() {
    let (a, rho, k, chord) = (0.2, 0.15, 10.0, 1e-4);
    let ((x0, x1), (y0, y1)) = ((-0.5, 0.6), (-0.4, 0.7));
    let c0 = 2.0 * a * a - rho * rho;
    let patch = polynomial(
        &[
            (3, 1, k),
            (1, 3, k),
            (2, 1, -2.0 * a * k),
            (1, 2, -2.0 * a * k),
            (1, 1, c0 * k),
        ],
        (x0, x1),
        (y0, y1),
    );
    let circle = closed(4000, core::f64::consts::TAU, |t| {
        Point::new(a + rho * t.cos(), a + rho * t.sin(), 0.0)
    });
    let truth = Truth {
        curves: vec![
            segment(Point::new(x0, 0.0, 0.0), Point::new(x1, 0.0, 0.0), 4000),
            segment(Point::new(0.0, y0, 0.0), Point::new(0.0, y1, 0.0), 4000),
            circle,
        ],
        ends: vec![Point::ORIGIN],
    };
    check(
        "axes and circle",
        (&patch, &ground()),
        chord,
        &truth,
        (10.0 * chord, 0),
    );
}

/// Two equal cylinders crossed square: two ellipses crossing where the
/// cylinders touch.
#[test]
fn equal_cylinders_crossed_square() {
    let (r, chord) = (30.0, 1e-6);
    let a: SurfaceGeometry = CylinderSurface::new(
        Cylinder::new(Frame::WORLD, r, T).unwrap(),
        (-3.0 * r, 3.0 * r),
    )
    .unwrap()
    .into();
    let along_x = Frame::new(Point::ORIGIN, Direction::X, Direction::Y, T).unwrap();
    let b: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(along_x, r, T).unwrap(), (-3.0 * r, 3.0 * r))
            .unwrap()
            .into();
    let ellipse = |sign: f64| {
        closed(20_000, core::f64::consts::TAU, |t| {
            Point::new(r * t.cos(), r * t.sin(), sign * r * t.cos())
        })
    };
    let truth = Truth {
        curves: vec![ellipse(1.0), ellipse(-1.0)],
        ends: vec![Point::new(0.0, r, 0.0), Point::new(0.0, -r, 0.0)],
    };
    check(
        "crossed cylinders",
        (&a, &b),
        chord,
        &truth,
        (10.0 * chord, 0),
    );
}

/// A sphere of radius `2r` cut by a cylinder of radius `r` through its
/// centre: Viviani's curve, a figure eight crossing itself where the
/// cylinder touches the sphere.
#[test]
fn vivianis_curve() {
    let (r, chord) = (10.0, 1e-6);
    let sphere: SurfaceGeometry =
        SphereSurface::new(Sphere::centred(Point::ORIGIN, 2.0 * r, T).unwrap()).into();
    let frame = Frame::new(Point::new(r, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let drum: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(frame, r, T).unwrap(), (-3.0 * r, 3.0 * r))
            .unwrap()
            .into();
    let curve = closed(20_000, 2.0 * core::f64::consts::TAU, |t| {
        Point::new(r * (1.0 + t.cos()), r * t.sin(), 2.0 * r * (t / 2.0).sin())
    });
    let truth = Truth {
        curves: vec![curve],
        ends: vec![Point::new(2.0 * r, 0.0, 0.0)],
    };
    // Two walks stall where they meet halfway round a lobe, each end
    // beside the other.
    check(
        "Viviani",
        (&sphere, &drum),
        chord,
        &truth,
        (10.0 * chord, 2),
    );
}
