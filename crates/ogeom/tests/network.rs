//! A network of patches fitted to a scan at once: a cylinder band in
//! sectors and rows, a box's edge rounded by a fillet beside a sharp one,
//! the corner of a box rounded on all three edges, and a three-sided
//! patch, each against the shape it samples, known in closed form.
#![allow(clippy::unwrap_used, reason = "test code")]
#![allow(clippy::cast_precision_loss, reason = "small counts")]

use core::f64::consts::{FRAC_PI_2, PI, TAU};
use std::time::Duration;

use ogeom::algo::{
    ChainedPatch, CurveDef, Network, NetworkOptions, PatchDef, PatchSide, check, fit_network,
    fit_network_chained,
};
use ogeom::core::clock::Instant;
use ogeom::core::{OgeomError, Tolerances};
use ogeom::geom::fit::Ended;
use ogeom::geom::{BSplineSurface, CircleCurve, Curve, Surface as _};
use ogeom::math::{Circle, Direction, Frame, Point};
use ogeom::mesh::PointIndex;
use ogeom::topo::Model;

const T: Tolerances = Tolerances::millimetres();

/// The noise on the noisy clouds, standard deviation in millimetres.
const NOISE: f64 = 0.02;

/// The tolerance a noisy cloud is held to: five deviations of the noise,
/// past what it reaches over twenty thousand points (four and a half).
const NOISY_TOLERANCE: f64 = 5.0 * NOISE;

/// The tolerance a clean cloud is held to.
const CLEAN_TOLERANCE: f64 = 1e-3;

/// The widest a gap along a shared curve may be: the edge and both
/// patches' boundaries are evaluated from the same poles, so they differ
/// by the rounding of the evaluation, a few units in the last place of
/// coordinates about ten.
const ROUNDING_GAP: f64 = 1e-12;

/// The widest angle across a curve held tangent: its pole triples are
/// solved for exactly, so the normals differ by rounding.
const ROUNDING_ANGLE: f64 = 1e-9;

/// A fixed sequence of draws (xorshift), so every run sees the same points.
struct Draws(u64);

impl Draws {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// Uniform in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1_u64 << 53) as f64
    }

    /// Standard normal, by Box-Muller.
    fn gaussian(&mut self) -> f64 {
        let u = 1.0 - self.unit();
        let v = self.unit();
        (-2.0 * u.ln()).sqrt() * (TAU * v).cos()
    }
}

fn options(tolerance: f64, creases: Vec<usize>) -> NetworkOptions {
    NetworkOptions {
        tolerance,
        smoothing: 1e-12,
        creases,
        spans: 1,
        most: 16,
    }
}

fn fit(
    curves: &[CurveDef],
    patches: &[PatchDef],
    cloud: &[Point],
    options: &NetworkOptions,
) -> (Model, Network) {
    let index = PointIndex::build(cloud).unwrap();
    let mut model = Model::new();
    let network = fit_network(&mut model, patches, curves, &index, options, T).unwrap();
    (model, network)
}

/// The widest distance from the network's surfaces, read on a 41 by 41
/// grid of each, to the true shape as `off` measures it.
fn widest_off(network: &Network, off: impl Fn(Point) -> f64) -> f64 {
    let mut widest = 0.0_f64;
    for surface in &network.surfaces {
        for i in 0..=40 {
            for j in 0..=40 {
                let p = surface
                    .point_at(f64::from(i) / 40.0, f64::from(j) / 40.0, T)
                    .unwrap();
                widest = widest.max(off(p));
            }
        }
    }
    widest
}

/// The shell is sound: nothing but its open border is reported, and that
/// as `open` free edges.
fn assert_sound_shell(model: &Model, network: &Network, open: usize) {
    let diagnosis = check(model, &network.shell, T).unwrap();
    assert!(diagnosis.is_usable(), "{diagnosis:?}");
    assert_eq!(diagnosis.problems.len(), 1, "{diagnosis:?}");
    assert!(
        diagnosis.problems[0]
            .what
            .starts_with(&format!("{open} edge(s) used an odd number of times")),
        "{diagnosis:?}"
    );
    assert_eq!(network.faces.len(), network.surfaces.len());
}

const RADIUS: f64 = 10.0;
const HEIGHT: f64 = 8.0;

fn on_cylinder(angle: f64, z: f64) -> Point {
    Point::new(RADIUS * angle.cos(), RADIUS * angle.sin(), z)
}

/// A band of the cylinder of radius 10 round `z`, 8 high, in `sectors`
/// round it and `rows` up it: its curves (the seams as segments, the arcs
/// as circles) and patches.
fn band(sectors: usize, rows: usize) -> (Vec<CurveDef>, Vec<PatchDef>) {
    let angle = |k: usize| TAU * (k % sectors) as f64 / sectors as f64;
    let height = |l: usize| HEIGHT * l as f64 / rows as f64;
    let mut curves = Vec::new();
    for k in 0..sectors {
        for l in 0..rows {
            curves.push(CurveDef::Polyline(vec![
                on_cylinder(angle(k), height(l)),
                on_cylinder(angle(k), height(l + 1)),
            ]));
        }
    }
    for k in 0..sectors {
        for l in 0..=rows {
            let frame = Frame::new(
                Point::new(0.0, 0.0, height(l)),
                Direction::Z,
                Direction::X,
                T,
            )
            .unwrap();
            curves.push(CurveDef::Curve {
                curve: Curve::Circle(CircleCurve::new(Circle::new(frame, RADIUS, T).unwrap())),
                range: (angle(k), angle(k) + TAU / sectors as f64),
            });
        }
    }
    let seam = |k: usize, l: usize| (k % sectors) * rows + l;
    let arc = |k: usize, l: usize| sectors * rows + k * (rows + 1) + l;
    let mut patches = Vec::new();
    for k in 0..sectors {
        for l in 0..rows {
            patches.push(PatchDef::four([
                (arc(k, l), false),
                (seam(k + 1, l), false),
                (arc(k, l + 1), true),
                (seam(k, l), true),
            ]));
        }
    }
    (curves, patches)
}

/// `count` points on the band, off it radially by noise of deviation
/// `noise`.
fn band_points(count: usize, noise: f64, seed: u64) -> Vec<Point> {
    let mut draws = Draws(seed);
    (0..count)
        .map(|_| {
            let angle = TAU * draws.unit();
            let z = HEIGHT * draws.unit();
            let r = noise.mul_add(draws.gaussian(), RADIUS);
            Point::new(r * angle.cos(), r * angle.sin(), z)
        })
        .collect()
}

fn off_cylinder(p: Point) -> f64 {
    (p.x.hypot(p.y) - RADIUS).abs()
}

#[test]
fn a_band_in_two_to_four_sectors_fits_the_cylinder_tangent_and_shut() {
    let cloud = band_points(20_000, 0.0, 7);
    for sectors in 2..=4 {
        let (curves, patches) = band(sectors, 1);
        let (model, network) = fit(&curves, &patches, &cloud, &options(CLEAN_TOLERANCE, vec![]));
        assert_eq!(network.ended, Ended::Met, "{network:?}");
        let off = widest_off(&network, off_cylinder);
        assert!(
            off <= CLEAN_TOLERANCE,
            "{sectors} sectors: {off} off the cylinder"
        );
        assert!(network.worst_gap <= ROUNDING_GAP, "{}", network.worst_gap);
        assert!(
            network.worst_tangent <= ROUNDING_ANGLE,
            "{}",
            network.worst_tangent
        );
        assert!(network.corners.is_empty());
        assert_eq!(network.outside, 0);
        assert_eq!(
            network.per_patch.iter().map(|p| p.points).sum::<usize>(),
            cloud.len()
        );
        // Two free arcs per sector, top and bottom.
        assert_sound_shell(&model, &network, 2 * sectors);
    }
}

#[test]
fn corners_of_four_patches_hold_the_band_tangent_through_them() {
    // Three sectors in two rows: the middle ring's three corners each join
    // four patches, where the seams and the middle arcs cross.
    let cloud = band_points(20_000, 0.0, 7);
    let (curves, patches) = band(3, 2);
    let (model, network) = fit(&curves, &patches, &cloud, &options(CLEAN_TOLERANCE, vec![]));
    assert_eq!(network.ended, Ended::Met, "{network:?}");
    assert!(widest_off(&network, off_cylinder) <= CLEAN_TOLERANCE);
    assert!(network.worst_gap <= ROUNDING_GAP, "{}", network.worst_gap);
    assert!(
        network.worst_tangent <= ROUNDING_ANGLE,
        "{}",
        network.worst_tangent
    );
    assert!(network.corners.is_empty());
    assert_sound_shell(&model, &network, 6);
}

#[test]
fn a_noisy_band_fits_within_its_tolerance_and_not_the_noise() {
    let cloud = band_points(20_000, NOISE, 11);
    let (curves, patches) = band(3, 2);
    let (_, network) = fit(&curves, &patches, &cloud, &options(NOISY_TOLERANCE, vec![]));
    let off = widest_off(&network, off_cylinder);
    assert!(off <= NOISY_TOLERANCE / 2.0, "{off} off the cylinder");
    assert!(
        network.worst_tangent <= ROUNDING_ANGLE,
        "{}",
        network.worst_tangent
    );
    for patch in &network.per_patch {
        // A fit that followed the noise would come nearer the points than
        // the noise stands off the shape.
        assert!(patch.sigma >= 0.9 * NOISE, "{patch:?}");
        assert!(patch.worst <= NOISY_TOLERANCE, "{patch:?}");
    }
}

#[test]
fn the_same_points_fit_the_same_network_to_the_bit() {
    let cloud = band_points(20_000, NOISE, 17);
    let (curves, patches) = band(3, 2);
    let (_, first) = fit(&curves, &patches, &cloud, &options(NOISY_TOLERANCE, vec![]));
    let (_, second) = fit(&curves, &patches, &cloud, &options(NOISY_TOLERANCE, vec![]));
    assert_eq!(first.surfaces, second.surfaces);
    assert_eq!(first.per_patch, second.per_patch);
    assert_eq!(
        first.worst_tangent.to_bits(),
        second.worst_tangent.to_bits()
    );
}

/// The fillet's radius, and the box's length along it.
const FILLET: f64 = 3.0;
const LENGTH: f64 = 20.0;

/// The box's profile across its rounded edge, in `(y, z)`, segment by
/// segment at `s` from nought to one: the front face up `y = 0` to
/// `z = 7`, the fillet round `(3, 7)`, the top at `z = 10` to the sharp
/// edge at `y = 15`, and the back face down from it.
fn profile(segment: usize, s: f64) -> (f64, f64) {
    match segment {
        0 => (0.0, 7.0 * s),
        1 => {
            let a = FRAC_PI_2.mul_add(-s, PI);
            (
                FILLET.mul_add(a.cos(), FILLET),
                FILLET.mul_add(a.sin(), 7.0),
            )
        }
        2 => (12.0f64.mul_add(s, 3.0), 10.0),
        _ => (15.0, 10.0 * (1.0 - s)),
    }
}

fn profile_length(segment: usize) -> f64 {
    [7.0, FRAC_PI_2 * FILLET, 12.0, 10.0][segment]
}

/// The distance from `(y, z)` to the profile.
fn off_profile(y: f64, z: f64) -> f64 {
    let segment = |a: (f64, f64), b: (f64, f64)| {
        let (dy, dz) = (b.0 - a.0, b.1 - a.1);
        let t = ((y - a.0).mul_add(dy, (z - a.1) * dz) / dy.mul_add(dy, dz * dz)).clamp(0.0, 1.0);
        (t.mul_add(-dy, y - a.0)).hypot(t.mul_add(-dz, z - a.1))
    };
    let (cy, cz) = (y - FILLET, z - 7.0);
    let arc = if cy <= 0.0 && cz >= 0.0 {
        (cy.hypot(cz) - FILLET).abs()
    } else {
        f64::INFINITY
    };
    segment((0.0, 0.0), (0.0, 7.0))
        .min(segment((3.0, 10.0), (15.0, 10.0)))
        .min(segment((15.0, 10.0), (15.0, 0.0)))
        .min(arc)
}

/// The box edge's network: a patch per profile segment, `u` along `x`.
/// Curves 0 to 4 run along `x` at the profile's breaks (1 and 2 the
/// fillet's borders, 3 the sharp edge); 5 to 8 are the profile at `x = 0`,
/// 9 to 12 at the far end.
fn box_edge() -> (Vec<CurveDef>, Vec<PatchDef>) {
    let at = |segment: usize, s: f64, x: f64| {
        let (y, z) = profile(segment, s);
        Point::new(x, y, z)
    };
    let mut curves = Vec::new();
    for b in 0..5 {
        let (segment, s) = if b < 4 { (b, 0.0) } else { (3, 1.0) };
        curves.push(CurveDef::Polyline(vec![
            at(segment, s, 0.0),
            at(segment, s, LENGTH),
        ]));
    }
    for x in [0.0, LENGTH] {
        for segment in 0..4 {
            curves.push(CurveDef::Polyline(
                (0..=64)
                    .map(|i| at(segment, f64::from(i) / 64.0, x))
                    .collect(),
            ));
        }
    }
    let patches = (0..4)
        .map(|i| PatchDef::four([(i, false), (9 + i, false), (i + 1, true), (5 + i, true)]))
        .collect();
    (curves, patches)
}

/// `count` points over the box edge, off it along its normal by noise of
/// deviation `noise`.
fn box_points(count: usize, noise: f64, seed: u64) -> Vec<Point> {
    let total: f64 = (0..4).map(profile_length).sum();
    let mut draws = Draws(seed);
    (0..count)
        .map(|_| {
            let mut along = total * draws.unit();
            let mut segment = 0;
            while segment < 3 && along > profile_length(segment) {
                along -= profile_length(segment);
                segment += 1;
            }
            let s = along / profile_length(segment);
            let (y, z) = profile(segment, s);
            let (ny, nz) = match segment {
                0 => (-1.0, 0.0),
                1 => {
                    let a = FRAC_PI_2.mul_add(-s, PI);
                    (a.cos(), a.sin())
                }
                2 => (0.0, 1.0),
                _ => (1.0, 0.0),
            };
            let e = noise * draws.gaussian();
            Point::new(LENGTH * draws.unit(), ny.mul_add(e, y), nz.mul_add(e, z))
        })
        .collect()
}

#[test]
fn a_box_edge_keeps_its_crease_and_rounds_its_fillet_smoothly() {
    let (curves, patches) = box_edge();
    for (noise, tolerance) in [(0.0, CLEAN_TOLERANCE), (NOISE, NOISY_TOLERANCE)] {
        let cloud = box_points(20_000, noise, 5);
        let (model, network) = fit(&curves, &patches, &cloud, &options(tolerance, vec![3]));
        let off = widest_off(&network, |p| off_profile(p.y, p.z));
        assert!(off <= tolerance, "noise {noise}: {off} off the box");
        assert!(network.worst_gap <= ROUNDING_GAP, "{}", network.worst_gap);
        // The fillet's borders are held tangent; the sharp edge is not a
        // tangent join and is not measured as one.
        assert!(
            network.worst_tangent <= ROUNDING_ANGLE,
            "{}",
            network.worst_tangent
        );
        // The sharp edge keeps its right angle. A face within the
        // tolerance of a plane tilts from it by at most twice the
        // tolerance over its width, the front face's 7 the narrowest.
        let tilt = 2.0 * tolerance / 7.0;
        for i in 0..=10 {
            let u = f64::from(i) / 10.0;
            let top = network.surfaces[2].normal_at(u, 1.0, T).unwrap().vector();
            let back = network.surfaces[3].normal_at(u, 0.0, T).unwrap().vector();
            let angle = top.cross(back).magnitude().atan2(top.dot(back));
            assert!(
                (angle - FRAC_PI_2).abs() <= tilt,
                "noise {noise}: the edge turns {angle} at {u}"
            );
        }
        // The two ends' profiles and the outer lines along `x`.
        assert_sound_shell(&model, &network, 10);
    }
}

const BALL: f64 = 5.0;

fn on_ball(x: f64, y: f64, z: f64) -> Point {
    let n = (x * x + y * y + z * z).sqrt();
    Point::new(BALL * x / n, BALL * y / n, BALL * z / n)
}

/// The great circle arc from `a` to `b` on the ball, as points.
fn ball_arc(a: Point, b: Point) -> CurveDef {
    let (a, b) = (a.to_vector(), b.to_vector());
    CurveDef::Polyline(
        (0..=64)
            .map(|i| {
                let t = f64::from(i) / 64.0;
                let v = a * (1.0 - t) + b * t;
                Point::ORIGIN + v * (BALL / v.magnitude())
            })
            .collect(),
    )
}

/// `count` points over the octant of the ball in `x, y, z ≥ 0`, off it
/// radially by noise of deviation `noise`.
fn octant_points(count: usize, noise: f64, seed: u64) -> Vec<Point> {
    let mut draws = Draws(seed);
    (0..count)
        .map(|_| {
            let (a, b, c) = (
                draws.gaussian().abs(),
                draws.gaussian().abs(),
                draws.gaussian().abs(),
            );
            let n = (a * a + b * b + c * c).sqrt();
            let r = noise.mul_add(draws.gaussian(), BALL);
            Point::new(r * a / n, r * b / n, r * c / n)
        })
        .collect()
}

fn off_ball(p: Point) -> f64 {
    (p.to_vector().magnitude() - BALL).abs()
}

/// The corner of a box rounded on its three edges: the octant of a ball
/// of radius 5, in three patches meeting at its middle. Curves 0 to 5 run
/// round its border, 6 to 8 in from the middles of its sides to its
/// middle.
fn rounded_corner() -> (Vec<CurveDef>, Vec<PatchDef>) {
    let (x, y, z) = (
        on_ball(1.0, 0.0, 0.0),
        on_ball(0.0, 1.0, 0.0),
        on_ball(0.0, 0.0, 1.0),
    );
    let (xy, yz, zx, middle) = (
        on_ball(1.0, 1.0, 0.0),
        on_ball(0.0, 1.0, 1.0),
        on_ball(1.0, 0.0, 1.0),
        on_ball(1.0, 1.0, 1.0),
    );
    let curves = vec![
        ball_arc(x, xy),
        ball_arc(x, zx),
        ball_arc(y, yz),
        ball_arc(y, xy),
        ball_arc(z, zx),
        ball_arc(z, yz),
        ball_arc(xy, middle),
        ball_arc(yz, middle),
        ball_arc(zx, middle),
    ];
    let patches = vec![
        PatchDef::four([(0, false), (6, false), (8, true), (1, true)]),
        PatchDef::four([(2, false), (7, false), (6, true), (3, true)]),
        PatchDef::four([(4, false), (8, false), (7, true), (5, true)]),
    ];
    (curves, patches)
}

/// The widest angle a curve held tangent by its spans may show, read off
/// the patches: the rows are solved for exactly, so the normals differ by
/// the rounding of the elimination and of evaluating unit normals of
/// patches with coordinates about five, a few hundred units in the last
/// place at most.
const EXACT_ANGLE: f64 = 1e-12;

/// The pole `w` along side `k` of a patch's net, walked as [`PatchDef`]'s
/// sides are.
fn side_pole(surface: &BSplineSurface, k: usize, w: usize) -> Point {
    let (nu, nv) = (surface.grid().u_count(), surface.grid().v_count());
    let (i, j) = match k {
        0 => (w, 0),
        1 => (nu - 1, w),
        2 => (nu - 1 - w, nv - 1),
        _ => (0, nv - 1 - w),
    };
    surface.grid().get(i, j).unwrap().point()
}

/// The chart point at share `f` of the walk along side `k`.
fn on_side(k: usize, f: f64) -> (f64, f64) {
    match k {
        0 => (f, 0.0),
        1 => (1.0, f),
        2 => (1.0 - f, 1.0),
        _ => (0.0, 1.0 - f),
    }
}

/// Across every curve two of `patches` share, read off the fitted patches
/// at 201 points along it, its ends included: the widest angle between
/// the two patches' normals and the widest distance between their points.
/// Each such curve's poles must be the same in both patches' nets, to the
/// bit.
fn across_shared(network: &Network, patches: &[PatchDef]) -> (f64, f64) {
    let mut users: Vec<Vec<(usize, usize, bool)>> = Vec::new();
    for (p, patch) in patches.iter().enumerate() {
        for (k, side) in patch.sides.iter().enumerate() {
            let Some(side) = side else { continue };
            if users.len() <= side.curve {
                users.resize(side.curve + 1, Vec::new());
            }
            users[side.curve].push((p, k, side.reversed));
        }
    }
    let (mut angle, mut gap) = (0.0_f64, 0.0_f64);
    for pair in users.iter().filter(|u| u.len() == 2) {
        let [(p, kp, rp), (q, kq, rq)] = [pair[0], pair[1]];
        let (a, b) = (&network.surfaces[p], &network.surfaces[q]);
        let along = |k: usize, s: &BSplineSurface| {
            if k.is_multiple_of(2) {
                s.grid().u_count()
            } else {
                s.grid().v_count()
            }
        };
        let n = along(kp, a);
        assert_eq!(n, along(kq, b));
        for w in 0..n {
            // The curve's pole `w` in its own order, along either walk.
            let (wp, wq) = (
                if rp { n - 1 - w } else { w },
                if rq { n - 1 - w } else { w },
            );
            assert_eq!(
                side_pole(a, kp, wp),
                side_pole(b, kq, wq),
                "curve shared by {p}, {q}"
            );
        }
        for i in 0..=200 {
            let t = f64::from(i) / 200.0;
            let (ca, cb) = (
                on_side(kp, if rp { 1.0 - t } else { t }),
                on_side(kq, if rq { 1.0 - t } else { t }),
            );
            let (na, nb) = (
                a.normal_at(ca.0, ca.1, T).unwrap().vector(),
                b.normal_at(cb.0, cb.1, T).unwrap().vector(),
            );
            angle = angle.max(na.cross(nb).magnitude().atan2(na.dot(nb)));
            let (pa, pb) = (
                a.point_at(ca.0, ca.1, T).unwrap(),
                b.point_at(cb.0, cb.1, T).unwrap(),
            );
            gap = gap.max(pa.distance(pb));
        }
    }
    (angle, gap)
}

#[test]
fn a_clean_rounded_corner_is_held_exactly_tangent_by_its_spans() {
    // The ball's octant from a clean scan, to a thousandth: holding its
    // corner with `γ` linear along each whole curve into it binds each
    // curve to a single cubic and misses the points by 0.00107 at any
    // spans, so the corner is held by the span of each curve next to it.
    let (curves, patches) = rounded_corner();
    let cloud = octant_points(20_000, 0.0, 3);
    let (model, network) = fit(&curves, &patches, &cloud, &options(CLEAN_TOLERANCE, vec![]));
    assert_eq!(network.ended, Ended::Met, "{network:?}");
    assert!(network.per_patch.iter().all(|p| p.worst <= CLEAN_TOLERANCE));
    let off = widest_off(&network, off_ball);
    assert!(off <= CLEAN_TOLERANCE, "{off} off the ball");
    let (angle, gap) = across_shared(&network, &patches);
    assert!(angle <= EXACT_ANGLE, "{angle}");
    assert!(gap <= ROUNDING_GAP, "{gap}");
    assert!(
        network.worst_tangent <= EXACT_ANGLE,
        "{}",
        network.worst_tangent
    );
    let [corner] = network.corners.as_slice() else {
        panic!("{:?}", network.corners);
    };
    assert_eq!(corner.valence, 3);
    assert!(corner.angle <= EXACT_ANGLE, "{corner:?}");
    // The corner may slide along the ball from where the curves were drawn
    // to meet, but stays on it.
    assert!(off_ball(corner.point) <= CLEAN_TOLERANCE, "{corner:?}");
    // Refined past one span each way, so each curve into the corner runs
    // on past its span next to the corner, held there by the triples.
    assert!(
        network
            .per_patch
            .iter()
            .all(|p| p.spans.0 >= 2 && p.spans.1 >= 2),
        "{:?}",
        network.per_patch
    );
    assert_sound_shell(&model, &network, 6);
}

#[test]
fn a_corner_left_near_tangent_is_reported_with_its_angle() {
    // A cap of six patches from a clean scan, to a thousandth: an even
    // corner held by its spans misses the points at their seed's
    // parameters, and its gluing is not spread, so the curves into it are
    // fitted near tangent and the corner's widest angle reported, as read
    // off the patches.
    let (curves, patches) = cap(6);
    let cloud = cap_points(20_000, 0.0, 29);
    let (model, network) = fit(&curves, &patches, &cloud, &options(CLEAN_TOLERANCE, vec![]));
    assert_eq!(network.ended, Ended::Met, "{network:?}");
    assert!(widest_off(&network, off_ball) <= CLEAN_TOLERANCE);
    let [corner] = network.corners.as_slice() else {
        panic!("{:?}", network.corners);
    };
    assert_eq!(corner.valence, 6);
    let normals: Vec<_> = network
        .surfaces
        .iter()
        .map(|s| s.normal_at(0.0, 0.0, T).unwrap().vector())
        .collect();
    let mut widest = 0.0_f64;
    for (i, a) in normals.iter().enumerate() {
        for b in &normals[i + 1..] {
            widest = widest.max(a.cross(*b).magnitude().atan2(a.dot(*b)));
        }
    }
    // The same normals, evaluated the same way: equal to rounding.
    assert!(
        (corner.angle - widest).abs() <= EXACT_ANGLE,
        "{corner:?} against {widest}"
    );
    // Near tangent, not held: well past what rounding leaves.
    assert!(corner.angle > EXACT_ANGLE, "{corner:?}");
    assert!(network.worst_tangent >= corner.angle);
    assert_sound_shell(&model, &network, 12);
}

#[test]
fn clean_caps_of_five_and_seven_patches_are_held_exactly_tangent() {
    // Clean scans, to a thousandth. Holding the corner with `γ` linear
    // along each whole curve binds the curve to a single cubic, and by its
    // spans asks `γ` to fall to nought within a span, which the points at
    // their seed's parameters do not admit: the seed's gluing runs
    // linearly along the whole curve. Spread along the curve, every knot
    // inside its strip doubled, it refines as far as the points need.
    for n in [5, 7] {
        let (curves, patches) = cap(n);
        let cloud = cap_points(20_000, 0.0, 29);
        let (model, network) = fit(&curves, &patches, &cloud, &options(CLEAN_TOLERANCE, vec![]));
        assert_eq!(network.ended, Ended::Met, "{n}: {network:?}");
        assert!(network.per_patch.iter().all(|p| p.worst <= CLEAN_TOLERANCE));
        let off = widest_off(&network, off_ball);
        assert!(off <= CLEAN_TOLERANCE, "{n}: {off} off the ball");
        let (angle, gap) = across_shared(&network, &patches);
        assert!(angle <= EXACT_ANGLE, "{n}: {angle}");
        assert!(gap <= ROUNDING_GAP, "{n}: {gap}");
        assert!(
            network.worst_tangent <= EXACT_ANGLE,
            "{n}: {}",
            network.worst_tangent
        );
        let [corner] = network.corners.as_slice() else {
            panic!("{:?}", network.corners);
        };
        assert_eq!(corner.valence, n);
        assert!(corner.angle <= EXACT_ANGLE, "{corner:?}");
        assert!(off_ball(corner.point) <= CLEAN_TOLERANCE, "{corner:?}");
        // Refined past one span each way, so each curve into the corner
        // has a doubled knot inside it.
        assert!(
            network
                .per_patch
                .iter()
                .all(|p| p.spans.0 >= 2 && p.spans.1 >= 2),
            "{n}: {:?}",
            network.per_patch
        );
        assert_sound_shell(&model, &network, 2 * n);
    }
}

#[test]
fn a_three_sided_patch_shrinks_its_missing_side_to_a_point() {
    // The octant as one patch: along the equator from `x` to `y`, up to the
    // pole, the third side shrunk to the pole, and down to `x`.
    let curves = vec![
        ball_arc(on_ball(1.0, 0.0, 0.0), on_ball(0.0, 1.0, 0.0)),
        ball_arc(on_ball(0.0, 1.0, 0.0), on_ball(0.0, 0.0, 1.0)),
        ball_arc(on_ball(0.0, 0.0, 1.0), on_ball(1.0, 0.0, 0.0)),
    ];
    let side = |curve: usize| {
        Some(PatchSide {
            curve,
            reversed: false,
        })
    };
    let patches = vec![PatchDef {
        sides: [side(0), side(1), None, side(2)],
    }];
    let cloud = octant_points(20_000, 0.0, 3);
    let (model, network) = fit(&curves, &patches, &cloud, &options(CLEAN_TOLERANCE, vec![]));
    assert_eq!(network.ended, Ended::Met, "{network:?}");
    assert!(widest_off(&network, off_ball) <= CLEAN_TOLERANCE);
    let pole = network.surfaces[0].point_at(0.0, 1.0, T).unwrap();
    for i in 1..=10 {
        let at = network.surfaces[0]
            .point_at(f64::from(i) / 10.0, 1.0, T)
            .unwrap();
        assert!(at.distance(pole) <= ROUNDING_GAP, "{at:?}");
    }
    assert_sound_shell(&model, &network, 3);
}

#[test]
fn a_network_that_does_not_close_is_refused() {
    let (curves, patches) = band(3, 1);
    let cloud = band_points(1000, 0.0, 7);
    let index = PointIndex::build(&cloud).unwrap();
    let refused = |patches: &[PatchDef], curves: &[CurveDef]| {
        matches!(
            fit_network(
                &mut Model::new(),
                patches,
                curves,
                &index,
                &options(CLEAN_TOLERANCE, vec![]),
                T
            ),
            Err(OgeomError::Construction(_))
        )
    };
    // A shared seam walked the same way by both its patches.
    let mut same_way = patches.clone();
    same_way[0].sides[1] = Some(PatchSide {
        curve: 1,
        reversed: true,
    });
    assert!(refused(&same_way, &curves));
    // A seam bounding three patches.
    let mut third = patches.clone();
    third.push(patches[0]);
    assert!(refused(&third, &curves));
    // A seam standing off the corner its neighbours meet at.
    let mut apart = curves.clone();
    apart[1] = CurveDef::Polyline(vec![
        on_cylinder(TAU / 3.0 + 0.01, 0.0),
        on_cylinder(TAU / 3.0, HEIGHT),
    ]);
    assert!(refused(&patches, &apart));
    // A patch with two sides missing.
    let mut two_missing = patches.clone();
    two_missing[0].sides[1] = None;
    two_missing[0].sides[3] = None;
    assert!(refused(&two_missing, &curves));
}

#[test]
#[ignore = "heavy"]
fn six_patches_on_a_hundred_thousand_points_in_bounded_time() {
    let cloud = band_points(100_000, NOISE, 13);
    let (curves, patches) = band(3, 2);
    let start = Instant::now();
    let (_, network) = fit(&curves, &patches, &cloud, &options(NOISY_TOLERANCE, vec![]));
    let took = start.elapsed();
    eprintln!("six patches on 100k points: {took:?}, {:?}", network.ended);
    assert!(widest_off(&network, off_cylinder) <= NOISY_TOLERANCE / 2.0);
    assert!(took < Duration::from_secs(120), "{took:?}");
}

#[test]
fn the_network_reads_back_from_step_as_written() {
    let cloud = band_points(20_000, NOISE, 19);
    let (curves, patches) = band(3, 2);
    let (model, network) = fit(&curves, &patches, &cloud, &options(NOISY_TOLERANCE, vec![]));
    let mut document = ogeom::doc::Document::over(model);
    document.add_part("network", network.shell.clone());
    let text = ogeom::io::write_step(&document, T).unwrap();
    let import = ogeom::io::read_step(&text, T).unwrap();
    let [shell] = import.shells.as_slice() else {
        panic!("{} shells", import.shells.len());
    };
    let model = import.document.model();
    let diagnosis = check(model, shell, T).unwrap();
    assert!(diagnosis.is_usable(), "{diagnosis:?}");
    let faces = ogeom::topo::explore_unique(model, shell, ogeom::topo::ShapeType::Face).unwrap();
    assert_eq!(faces.len(), network.surfaces.len());
    // Each face's patch reads back as one of the fitted patches, pole for
    // pole to the precision the file prints, so the joins read back as
    // they were measured.
    for face in &faces {
        let data = model.node(face).unwrap().data().as_face().unwrap();
        let Some(ogeom::geom::SurfaceGeometry::BSpline(read)) =
            model.geometry().surface(data.surface)
        else {
            panic!("a face lost its patch");
        };
        let matches = network.surfaces.iter().any(|fitted| {
            let (a, b) = (read.grid().points(), fitted.grid().points());
            a.len() == b.len()
                && a.iter()
                    .zip(b)
                    .all(|(p, q)| p.point().distance(q.point()) <= ROUNDING_GAP)
        });
        assert!(matches, "a read patch matches none fitted");
    }
}

/// The band's sweep a T-junction network covers, and where along it the
/// junction stands: a third of the way, so the curves either side of it
/// differ in length.
const T_SWEEP: f64 = TAU / 3.0;
const T_AT: f64 = T_SWEEP / 3.0;

/// A circle arc of the band at height `z` from angle `a` to `b`.
fn band_arc(z: f64, a: f64, b: f64) -> CurveDef {
    let frame = Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X, T).unwrap();
    CurveDef::Curve {
        curve: Curve::Circle(CircleCurve::new(Circle::new(frame, RADIUS, T).unwrap())),
        range: (a, b),
    }
}

/// A seam of the band at angle `a` from height `z0` to `z1`.
fn band_seam(a: f64, z0: f64, z1: f64) -> CurveDef {
    CurveDef::Polyline(vec![on_cylinder(a, z0), on_cylinder(a, z1)])
}

/// A third of the cylinder band in three rows: a patch below, one above it
/// across the whole sweep, and two side by side on top, the seam between
/// them (curve 11) ending on the middle row's top side at a T-junction,
/// a third of the way along. That side is curves 2 and 3; its opposite,
/// curve 1, is one curve and carries the junction on to the bottom patch.
fn t_junction() -> (Vec<CurveDef>, Vec<ChainedPatch>) {
    let (h1, h2) = (HEIGHT / 3.0, 2.0 * HEIGHT / 3.0);
    let curves = vec![
        band_arc(0.0, 0.0, T_SWEEP),
        band_arc(h1, 0.0, T_SWEEP),
        band_arc(h2, 0.0, T_AT),
        band_arc(h2, T_AT, T_SWEEP),
        band_arc(HEIGHT, 0.0, T_AT),
        band_arc(HEIGHT, T_AT, T_SWEEP),
        band_seam(0.0, 0.0, h1),
        band_seam(T_SWEEP, 0.0, h1),
        band_seam(0.0, h1, h2),
        band_seam(T_SWEEP, h1, h2),
        band_seam(0.0, h2, HEIGHT),
        band_seam(T_AT, h2, HEIGHT),
        band_seam(T_SWEEP, h2, HEIGHT),
    ];
    let patches = vec![
        ChainedPatch::four([
            vec![(0, false)],
            vec![(7, false)],
            vec![(1, true)],
            vec![(6, true)],
        ]),
        ChainedPatch::four([
            vec![(1, false)],
            vec![(9, false)],
            vec![(3, true), (2, true)],
            vec![(8, true)],
        ]),
        ChainedPatch::four([
            vec![(2, false)],
            vec![(11, false)],
            vec![(4, true)],
            vec![(10, true)],
        ]),
        ChainedPatch::four([
            vec![(3, false)],
            vec![(12, false)],
            vec![(5, true)],
            vec![(11, true)],
        ]),
    ];
    (curves, patches)
}

/// `count` points on the band's sweep from angle nought to `T_SWEEP`.
fn sweep_points(count: usize, seed: u64) -> Vec<Point> {
    let mut draws = Draws(seed);
    (0..count)
        .map(|_| on_cylinder(T_SWEEP * draws.unit(), HEIGHT * draws.unit()))
        .collect()
}

fn fit_chained(
    curves: &[CurveDef],
    patches: &[ChainedPatch],
    cloud: &[Point],
    corrections: usize,
) -> (Model, Network) {
    let index = PointIndex::build(cloud).unwrap();
    let mut model = Model::new();
    let network = fit_network_chained(
        &mut model,
        patches,
        curves,
        &index,
        &options(CLEAN_TOLERANCE, vec![]),
        corrections,
        T,
    )
    .unwrap();
    (model, network)
}

#[test]
fn a_t_junction_shares_its_poles_and_holds_tangent() {
    let (curves, patches) = t_junction();
    let cloud = sweep_points(20_000, 23);
    let (model, network) = fit_chained(&curves, &patches, &cloud, 0);
    assert_eq!(network.ended, Ended::Met, "{network:?}");
    let off = widest_off(&network, off_cylinder);
    assert!(off <= CLEAN_TOLERANCE, "{off} off the cylinder");
    // Every curve meets its patches' boundaries to rounding: the halves of
    // the split side, the seam ending at the junction, and the curve below
    // that carries the junction through to the bottom patch.
    assert!(network.worst_gap <= ROUNDING_GAP, "{}", network.worst_gap);
    // Tangent across every curve, the seam ending at the junction too.
    assert!(
        network.worst_tangent <= ROUNDING_ANGLE,
        "{}",
        network.worst_tangent
    );
    // A patch counts twice at a junction its side runs straight through:
    // the junction is a regular corner, not an irregular one.
    assert!(network.corners.is_empty(), "{:?}", network.corners);
    // The middle and bottom patches take a knot of multiplicity three at
    // the junction, and are smooth across it. The junction's share of the
    // way along is the drawn arcs' lengths' (a third of the sweep against
    // two), each drawn to chords within the tolerance of the circle, which
    // shortens it by a far smaller fraction than the bound.
    for p in [0, 1] {
        let surface = &network.surfaces[p];
        let knots = surface.u_knots().knots();
        let junction = knots
            .iter()
            .copied()
            .find(|&k| k > 0.0 && k < 1.0 && knots.iter().filter(|&&x| x == k).count() == 3)
            .unwrap();
        assert!((junction - 1.0 / 3.0).abs() <= 1e-3, "{junction}");
        for i in 0..=10 {
            let v = f64::from(i) / 10.0;
            let before = surface.normal_at(junction.next_down(), v, T).unwrap();
            let after = surface.normal_at(junction, v, T).unwrap();
            let (a, b) = (before.vector(), after.vector());
            let angle = a.cross(b).magnitude().atan2(a.dot(b));
            assert!(angle <= ROUNDING_ANGLE, "patch {p} turns {angle} at {v}");
        }
    }
    // The border: the arcs below and above and the seams at either end.
    assert_sound_shell(&model, &network, 9);
}

/// The band's sweep in three rows with nested junctions: a patch below
/// across the whole sweep, its top split at a third of the way along
/// (curves 1 and 2) by the seam between the two patches of the middle row
/// (curve 12), and the left of those split again at half its width
/// (curves 3 and 4) by the seam between two of the three top patches
/// (curve 15). Curve 1 is a curve of the bottom patch's split side and
/// carries that second junction down into it, which crosses three strips.
fn nested_junctions() -> (Vec<CurveDef>, Vec<ChainedPatch>) {
    let (h1, h2) = (HEIGHT / 3.0, 2.0 * HEIGHT / 3.0);
    let half = T_AT / 2.0;
    let curves = vec![
        band_arc(0.0, 0.0, T_SWEEP),
        band_arc(h1, 0.0, T_AT),
        band_arc(h1, T_AT, T_SWEEP),
        band_arc(h2, 0.0, half),
        band_arc(h2, half, T_AT),
        band_arc(h2, T_AT, T_SWEEP),
        band_arc(HEIGHT, 0.0, half),
        band_arc(HEIGHT, half, T_AT),
        band_arc(HEIGHT, T_AT, T_SWEEP),
        band_seam(0.0, 0.0, h1),
        band_seam(T_SWEEP, 0.0, h1),
        band_seam(0.0, h1, h2),
        band_seam(T_AT, h1, h2),
        band_seam(T_SWEEP, h1, h2),
        band_seam(0.0, h2, HEIGHT),
        band_seam(half, h2, HEIGHT),
        band_seam(T_AT, h2, HEIGHT),
        band_seam(T_SWEEP, h2, HEIGHT),
    ];
    let patch = |bottom: Vec<(usize, bool)>, right: usize, top: Vec<(usize, bool)>, left: usize| {
        ChainedPatch::four([bottom, vec![(right, false)], top, vec![(left, true)]])
    };
    let patches = vec![
        patch(vec![(0, false)], 10, vec![(2, true), (1, true)], 9),
        patch(vec![(1, false)], 12, vec![(4, true), (3, true)], 11),
        patch(vec![(2, false)], 13, vec![(5, true)], 12),
        patch(vec![(3, false)], 15, vec![(6, true)], 14),
        patch(vec![(4, false)], 16, vec![(7, true)], 15),
        patch(vec![(5, false)], 17, vec![(8, true)], 16),
    ];
    (curves, patches)
}

#[test]
fn nested_junctions_share_their_poles_and_hold_tangent() {
    let (curves, patches) = nested_junctions();
    let cloud = sweep_points(20_000, 23);
    let (model, network) = fit_chained(&curves, &patches, &cloud, 0);
    assert_eq!(network.ended, Ended::Met, "{network:?}");
    let off = widest_off(&network, off_cylinder);
    assert!(off <= CLEAN_TOLERANCE, "{off} off the cylinder");
    // Curve 1 crosses two strips, and the bottom patch's side along the
    // arcs below it three; every curve the same poles in both patches.
    let mut parts = vec![1; curves.len()];
    parts[1] = 2;
    parts[0] = 3;
    let (angle, gap) = across_chained(&network, &patches, &parts);
    assert!(angle <= EXACT_ANGLE, "{angle}");
    assert!(gap <= ROUNDING_GAP, "{gap}");
    assert!(network.worst_gap <= ROUNDING_GAP, "{}", network.worst_gap);
    assert!(
        network.worst_tangent <= EXACT_ANGLE,
        "{}",
        network.worst_tangent
    );
    assert!(network.corners.is_empty(), "{:?}", network.corners);
    // The bottom patch takes a knot of multiplicity three at both
    // junctions, a sixth and a third of the way along, and the middle
    // left one at the second: each is held smooth across them.
    let shares = side_shares(&network.surfaces[0], 0);
    assert_eq!(shares.len(), 4, "{shares:?}");
    assert!((shares[1] - 1.0 / 6.0).abs() <= 1e-3, "{shares:?}");
    assert!((shares[2] - 1.0 / 3.0).abs() <= 1e-3, "{shares:?}");
    for p in [0, 1] {
        let line = across_junction_line(&network, p);
        assert!(
            line <= EXACT_ANGLE,
            "patch {p} turns {line} across a junction"
        );
    }
    // The border: the arcs below and above and the seams at either end.
    assert_sound_shell(&model, &network, 10);
}

#[test]
fn sides_split_unlike_are_refused() {
    let (curves, patches) = t_junction();
    let cloud = sweep_points(1000, 23);
    let index = PointIndex::build(&cloud).unwrap();
    let refused = |patches: &[ChainedPatch]| {
        matches!(
            fit_network_chained(
                &mut Model::new(),
                patches,
                &curves,
                &index,
                &options(CLEAN_TOLERANCE, vec![]),
                0,
                T
            ),
            Err(OgeomError::Construction(_))
        )
    };
    // The bottom side of the middle patch split in three while its top is
    // split in two.
    let mut unlike = patches.clone();
    unlike[1].sides[0] = vec![
        PatchSide {
            curve: 1,
            reversed: false,
        },
        PatchSide {
            curve: 6,
            reversed: false,
        },
        PatchSide {
            curve: 7,
            reversed: false,
        },
    ];
    assert!(refused(&unlike));
}

#[test]
fn parameter_correction_brings_a_curved_band_nearer_its_points() {
    let cloud = band_points(20_000, 0.0, 7);
    let (curves, patches) = band(3, 2);
    let patches: Vec<ChainedPatch> = patches.into_iter().map(ChainedPatch::from).collect();
    let sigma = |network: &Network| {
        let (sum, count) = network.per_patch.iter().fold((0.0, 0), |(s, n), p| {
            (p.sigma.mul_add(p.sigma * p.points as f64, s), n + p.points)
        });
        (sum / count as f64).sqrt()
    };
    let (_, seeded) = fit_chained(&curves, &patches, &cloud, 0);
    let (_, corrected) = fit_chained(&curves, &patches, &cloud, 3);
    assert!(sigma(&corrected) < sigma(&seeded));
    // The same spans, still within the tolerance, still tangent.
    assert_eq!(corrected.ended, Ended::Met);
    let spans = |n: &Network| n.per_patch.iter().map(|p| p.spans).collect::<Vec<_>>();
    assert_eq!(spans(&corrected), spans(&seeded));
    assert!(
        corrected
            .per_patch
            .iter()
            .all(|p| p.worst <= CLEAN_TOLERANCE)
    );
    assert!(corrected.worst_tangent <= ROUNDING_ANGLE);
}

#[test]
fn a_rounded_corner_of_three_patches_is_held_exactly_tangent() {
    // The ball's octant from a scan whose noise its tolerance allows for:
    // each curve into the corner runs on to the border, so the patches are
    // held tangent across it and at the corner, to rounding.
    let (curves, patches) = rounded_corner();
    let cloud = octant_points(20_000, NOISE, 3);
    let (model, network) = fit(&curves, &patches, &cloud, &options(NOISY_TOLERANCE, vec![]));
    assert_eq!(network.ended, Ended::Met, "{network:?}");
    let off = widest_off(&network, off_ball);
    assert!(off <= NOISY_TOLERANCE / 2.0, "{off} off the ball");
    assert!(network.worst_gap <= ROUNDING_GAP, "{}", network.worst_gap);
    assert!(
        network.worst_tangent <= ROUNDING_ANGLE,
        "{}",
        network.worst_tangent
    );
    let [corner] = network.corners.as_slice() else {
        panic!("{:?}", network.corners);
    };
    assert_eq!(corner.valence, 3);
    assert!(corner.angle <= ROUNDING_ANGLE, "{corner:?}");
    assert_sound_shell(&model, &network, 6);
}

/// The polar angle of the ball's cap a five-patch network covers.
const CAP: f64 = PI / 3.0;

/// The point of the ball at polar angle `polar` and azimuth `azimuth`.
fn on_ball_at(polar: f64, azimuth: f64) -> Point {
    Point::new(
        BALL * polar.sin() * azimuth.cos(),
        BALL * polar.sin() * azimuth.sin(),
        BALL * polar.cos(),
    )
}

/// A cap of the ball in five patches round its top, a corner of five: the
/// spokes (curves 0 to 4) run from the top to the middles of the five
/// arcs of the cap's rim, each arc two curves (5 to 14).
fn five_cap() -> (Vec<CurveDef>, Vec<PatchDef>) {
    cap(5)
}

/// A cap of the ball in `n` patches round its top, a corner of `n`, laid
/// out as [`five_cap`]'s.
fn cap(n: usize) -> (Vec<CurveDef>, Vec<PatchDef>) {
    let step = TAU / n as f64;
    let mut curves: Vec<CurveDef> = (0..n)
        .map(|k| ball_arc(on_ball(0.0, 0.0, 1.0), on_ball_at(CAP, step * k as f64)))
        .collect();
    let rim = |a: f64, b: f64| {
        CurveDef::Polyline(
            (0..=64)
                .map(|i| on_ball_at(CAP, (b - a).mul_add(f64::from(i) / 64.0, a)))
                .collect(),
        )
    };
    for k in 0..n {
        let a = step * k as f64;
        curves.push(rim(a, a + step / 2.0));
        curves.push(rim(a + step / 2.0, a + step));
    }
    let patches = (0..n)
        .map(|k| {
            PatchDef::four([
                (k, false),
                (n + 2 * k, false),
                (n + 1 + 2 * k, false),
                ((k + 1) % n, true),
            ])
        })
        .collect();
    (curves, patches)
}

/// `count` points over the cap, off the ball radially by noise of
/// deviation `noise`.
fn cap_points(count: usize, noise: f64, seed: u64) -> Vec<Point> {
    let mut draws = Draws(seed);
    (0..count)
        .map(|_| {
            // Even over the cap's area: its height is even.
            let z = (1.0 - CAP.cos()).mul_add(-draws.unit(), 1.0);
            let polar = z.acos();
            let azimuth = TAU * draws.unit();
            let r = noise.mul_add(draws.gaussian(), BALL);
            on_ball_at(polar, azimuth).to_vector() * (r / BALL)
        })
        .map(|v| Point::ORIGIN + v)
        .collect()
}

#[test]
fn a_corner_of_five_patches_is_held_exactly_tangent() {
    let (curves, patches) = five_cap();
    let cloud = cap_points(20_000, NOISE, 29);
    let (model, network) = fit(&curves, &patches, &cloud, &options(NOISY_TOLERANCE, vec![]));
    assert_eq!(network.ended, Ended::Met, "{network:?}");
    let off = widest_off(&network, off_ball);
    assert!(off <= NOISY_TOLERANCE / 2.0, "{off} off the ball");
    assert!(network.worst_gap <= ROUNDING_GAP, "{}", network.worst_gap);
    assert!(
        network.worst_tangent <= ROUNDING_ANGLE,
        "{}",
        network.worst_tangent
    );
    let [corner] = network.corners.as_slice() else {
        panic!("{:?}", network.corners);
    };
    assert_eq!(corner.valence, 5);
    assert!(corner.point.distance(on_ball(0.0, 0.0, 1.0)) <= NOISY_TOLERANCE);
    assert_sound_shell(&model, &network, 10);
}

/// The curved cover's height over `(x, y)`: curved both ways.
fn cover(x: f64, y: f64) -> f64 {
    8.0 - x * x / 100.0 - y * y / 60.0
}

/// The cover over `x` in ±30 and `y` in ±25, in six patches along `x` by
/// five along `y`: the lines along `x` first, then those along `y`, each
/// drawn as points on the cover.
fn cover_network() -> (Vec<CurveDef>, Vec<PatchDef>) {
    let xs: Vec<f64> = (0..=6).map(|i| -30.0 + 10.0 * f64::from(i)).collect();
    let ys: Vec<f64> = (0..=5).map(|j| -25.0 + 10.0 * f64::from(j)).collect();
    let line = |a: (f64, f64), b: (f64, f64)| {
        CurveDef::Polyline(
            (0..=32)
                .map(|k| {
                    let t = f64::from(k) / 32.0;
                    let (x, y) = ((b.0 - a.0).mul_add(t, a.0), (b.1 - a.1).mul_add(t, a.1));
                    Point::new(x, y, cover(x, y))
                })
                .collect(),
        )
    };
    let mut curves = Vec::new();
    // Along `x`: segment `i` of line `j`.
    for &y in &ys {
        for i in 0..6 {
            curves.push(line((xs[i], y), (xs[i + 1], y)));
        }
    }
    // Along `y`: segment `j` of line `i`.
    let base = curves.len();
    for &x in &xs {
        for j in 0..5 {
            curves.push(line((x, ys[j]), (x, ys[j + 1])));
        }
    }
    let along_x = |i: usize, j: usize| j * 6 + i;
    let along_y = |i: usize, j: usize| base + i * 5 + j;
    let mut patches = Vec::new();
    for j in 0..5 {
        for i in 0..6 {
            patches.push(PatchDef::four([
                (along_x(i, j), false),
                (along_y(i + 1, j), false),
                (along_x(i, j + 1), true),
                (along_y(i, j), true),
            ]));
        }
    }
    (curves, patches)
}

#[test]
fn a_thirty_patch_cover_reads_back_from_step_tangent_and_shut() {
    // Thirty patches over a cover curved both ways, tangent within a
    // thousandth of a degree and shut within a micrometre, read back from
    // STEP.
    let mut draws = Draws(31);
    let cloud: Vec<Point> = (0..60_000)
        .map(|_| {
            let (x, y) = (
                60.0f64.mul_add(draws.unit(), -30.0),
                50.0f64.mul_add(draws.unit(), -25.0),
            );
            Point::new(x, y, NOISE.mul_add(draws.gaussian(), cover(x, y)))
        })
        .collect();
    let (curves, patches) = cover_network();
    let (model, network) = fit(&curves, &patches, &cloud, &options(NOISY_TOLERANCE, vec![]));
    assert_eq!(network.surfaces.len(), 30);
    let off = widest_off(&network, |p| (p.z - cover(p.x, p.y)).abs());
    assert!(off <= NOISY_TOLERANCE / 2.0, "{off} off the cover");

    let mut document = ogeom::doc::Document::over(model);
    document.add_part("cover", network.shell.clone());
    let text = ogeom::io::write_step(&document, T).unwrap();
    let import = ogeom::io::read_step(&text, T).unwrap();
    let [shell] = import.shells.as_slice() else {
        panic!("{} shells", import.shells.len());
    };
    let model = import.document.model();
    let diagnosis = check(model, shell, T).unwrap();
    assert!(diagnosis.is_usable(), "{diagnosis:?}");
    let faces = ogeom::topo::explore_unique(model, shell, ogeom::topo::ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 30);
    // Each fitted patch's read-back patch, matched pole for pole.
    let read: Vec<ogeom::geom::BSplineSurface> = network
        .surfaces
        .iter()
        .map(|fitted| {
            faces
                .iter()
                .find_map(|face| {
                    let data = model.node(face).unwrap().data().as_face().unwrap();
                    let Some(ogeom::geom::SurfaceGeometry::BSpline(read)) =
                        model.geometry().surface(data.surface)
                    else {
                        return None;
                    };
                    let (a, b) = (read.grid().points(), fitted.grid().points());
                    (a.len() == b.len()
                        && a.iter()
                            .zip(b)
                            .all(|(p, q)| p.point().distance(q.point()) <= ROUNDING_GAP))
                    .then(|| read.clone())
                })
                .unwrap()
        })
        .collect();
    // Across every inner line, read off the patches read back: the
    // widest angle between their normals and the widest gap.
    let (mut angle, mut gap) = (0.0_f64, 0.0_f64);
    let mut across = |a: &ogeom::geom::BSplineSurface,
                      at_a: (f64, f64),
                      b: &ogeom::geom::BSplineSurface,
                      at_b: (f64, f64)| {
        let (na, nb) = (
            a.normal_at(at_a.0, at_a.1, T).unwrap().vector(),
            b.normal_at(at_b.0, at_b.1, T).unwrap().vector(),
        );
        angle = angle.max(na.cross(nb).magnitude().atan2(na.dot(nb)));
        let (pa, pb) = (
            a.point_at(at_a.0, at_a.1, T).unwrap(),
            b.point_at(at_b.0, at_b.1, T).unwrap(),
        );
        gap = gap.max(pa.distance(pb));
    };
    for j in 0..5 {
        for i in 0..6 {
            let here = &read[j * 6 + i];
            for k in 0..=40 {
                let s = f64::from(k) / 40.0;
                if i < 5 {
                    across(here, (1.0, s), &read[j * 6 + i + 1], (0.0, s));
                }
                if j < 4 {
                    across(here, (s, 1.0), &read[(j + 1) * 6 + i], (s, 0.0));
                }
            }
        }
    }
    assert!(angle <= 0.001_f64.to_radians(), "{angle}");
    assert!(gap <= 1e-3, "{gap}");
}

/// How far the box's edge fillets run back from the rounded corner.
const FILLET_RUN: f64 = 4.0;

/// `p` turned a third of a turn about the diagonal `k` times: `x` to `y`,
/// `y` to `z`, `z` to `x`.
fn turned(p: Point, k: usize) -> Point {
    (0..k % 3).fold(p, |p, _| Point::new(p.z, p.x, p.y))
}

/// The corner of a box rounded on its three edges with the edges' fillets
/// attached: the ball's octant of [`rounded_corner`] (curves 0 to 8,
/// patches 0 to 2), and along each edge, back from the octant, a quarter
/// of the cylinder of the ball's radius round that edge's axis, in two
/// patches either side of the line that runs straight on from the
/// octant's curve into the corner (curves `9 + 5k` to `13 + 5k` for the
/// fillet `k`). The ends of the octant's curves into its corner are now
/// corners of four patches.
fn filleted_corner() -> (Vec<CurveDef>, Vec<PatchDef>) {
    let (mut curves, mut patches) = rounded_corner();
    let line = |a: Point, b: Point| CurveDef::Polyline(vec![a, b]);
    let ring = |x: f64, from: f64, to: f64| {
        CurveDef::Polyline(
            (0..=64)
                .map(|i| {
                    let a = (to - from).mul_add(f64::from(i) / 64.0, from);
                    Point::new(x, BALL * a.cos(), BALL * a.sin())
                })
                .collect(),
        )
    };
    // The fillet round `x`, from the octant's side in `x = 0` back to
    // `x = −FILLET_RUN`, then turned onto the other two edges.
    let back = -FILLET_RUN;
    let (y, yz, z) = (
        Point::new(0.0, BALL, 0.0),
        on_ball(0.0, 1.0, 1.0),
        Point::new(0.0, 0.0, BALL),
    );
    let at_back = |p: Point| Point::new(back, p.y, p.z);
    let fillet = [
        line(yz, at_back(yz)),
        line(y, at_back(y)),
        line(z, at_back(z)),
        ring(back, 0.0, FRAC_PI_2 / 2.0),
        ring(back, FRAC_PI_2, FRAC_PI_2 / 2.0),
    ];
    // The octant's curves along its side in `x = 0`, from `y` and from `z`
    // to the side's middle, and those the turns take them to.
    let sides = [(2, 5), (4, 1), (0, 3)];
    for (k, &(from_y, from_z)) in sides.iter().enumerate() {
        let base = curves.len();
        for curve in &fillet {
            let CurveDef::Polyline(points) = curve else {
                unreachable!()
            };
            curves.push(CurveDef::Polyline(
                points.iter().map(|&p| turned(p, k)).collect(),
            ));
        }
        patches.push(PatchDef::four([
            (from_y, true),
            (base + 1, false),
            (base + 3, false),
            (base, true),
        ]));
        patches.push(PatchDef::four([
            (from_z, false),
            (base, false),
            (base + 4, true),
            (base + 2, true),
        ]));
    }
    (curves, patches)
}

/// `count` points over the filleted corner, even over its area, off it
/// along its normal by noise of deviation `noise`.
fn filleted_points(count: usize, noise: f64, seed: u64) -> Vec<Point> {
    let octant = PI * BALL * BALL / 2.0;
    let fillet = FRAC_PI_2 * BALL * FILLET_RUN;
    let mut draws = Draws(seed);
    let ball = octant_points(count, noise, seed ^ 0x9e37);
    let mut at = 0;
    (0..count)
        .map(|_| {
            let pick = (octant + 3.0 * fillet) * draws.unit();
            if pick < octant {
                at += 1;
                return ball[at - 1];
            }
            let k = (1..3)
                .filter(|&k| pick - octant > fillet * f64::from(k))
                .count();
            let a = FRAC_PI_2 * draws.unit();
            let r = noise.mul_add(draws.gaussian(), BALL);
            let p = Point::new(-FILLET_RUN * draws.unit(), r * a.cos(), r * a.sin());
            turned(p, k)
        })
        .collect()
}

/// The distance from `p` to the filleted corner: to the nearest of the
/// octant and the three quarter cylinders, each point of them taken at
/// `p`'s nearest point there.
fn off_filleted(p: Point) -> f64 {
    let ball = {
        let q = Point::new(p.x.max(0.0), p.y.max(0.0), p.z.max(0.0)).to_vector();
        p.distance(Point::ORIGIN + q * (BALL / q.magnitude()))
    };
    let fillet = |k: usize| {
        // Turned back so the fillet runs along `x`.
        let q = turned(p, 3 - k % 3);
        let (y, z) = (q.y.max(0.0), q.z.max(0.0));
        let r = y.hypot(z);
        let on = Point::new(q.x.clamp(-FILLET_RUN, 0.0), BALL * y / r, BALL * z / r);
        q.distance(on)
    };
    ball.min(fillet(0)).min(fillet(1)).min(fillet(2))
}

#[test]
fn a_filleted_box_corner_is_held_exactly_tangent() {
    // The rounded corner with its edges' fillets attached: each curve into
    // the corner runs on to a corner of four patches, where it carries on
    // straight across a fillet, so the corner is held by its spans and
    // every curve, past the span next to it, by the pole triples. At the
    // finer noise the octant's patches take two spans each way, so the
    // knot past each span next to the corner is doubled.
    let (curves, patches) = filleted_corner();
    for noise in [NOISE, NOISE / 4.0] {
        // Five deviations of the noise, as the other noisy clouds.
        let tolerance = 5.0 * noise;
        let cloud = filleted_points(40_000, noise, 37);
        let (model, network) = fit(&curves, &patches, &cloud, &options(tolerance, vec![]));
        assert_eq!(network.ended, Ended::Met, "noise {noise}: {network:?}");
        let off = widest_off(&network, off_filleted);
        assert!(
            off <= tolerance / 2.0,
            "noise {noise}: {off} off the corner"
        );
        let (angle, gap) = across_shared(&network, &patches);
        assert!(angle <= EXACT_ANGLE, "noise {noise}: {angle}");
        assert!(gap <= ROUNDING_GAP, "noise {noise}: {gap}");
        assert!(
            network.worst_tangent <= EXACT_ANGLE,
            "noise {noise}: {}",
            network.worst_tangent
        );
        let [corner] = network.corners.as_slice() else {
            panic!("{:?}", network.corners);
        };
        assert_eq!(corner.valence, 3);
        assert!(corner.angle <= EXACT_ANGLE, "{corner:?}");
        if noise < NOISE {
            assert!(
                network.per_patch[..3]
                    .iter()
                    .all(|p| p.spans.0 >= 2 && p.spans.1 >= 2),
                "{:?}",
                network.per_patch
            );
        }
        // The far ends of the fillets: two curves along each and the arcs
        // across their ends.
        assert_sound_shell(&model, &network, 12);
    }
}

#[test]
fn a_corner_of_six_patches_is_held_exactly_tangent() {
    // An even corner: the twists round it close in a cycle that is always
    // singular, and the rows holding the corner by its spans leave that
    // cycle's one condition on the curves' own poles.
    let (curves, patches) = cap(6);
    let cloud = cap_points(20_000, NOISE, 29);
    let (model, network) = fit(&curves, &patches, &cloud, &options(NOISY_TOLERANCE, vec![]));
    assert_eq!(network.ended, Ended::Met, "{network:?}");
    let off = widest_off(&network, off_ball);
    assert!(off <= NOISY_TOLERANCE / 2.0, "{off} off the ball");
    let (angle, gap) = across_shared(&network, &patches);
    assert!(angle <= EXACT_ANGLE, "{angle}");
    assert!(gap <= ROUNDING_GAP, "{gap}");
    let [corner] = network.corners.as_slice() else {
        panic!("{:?}", network.corners);
    };
    assert_eq!(corner.valence, 6);
    assert!(corner.angle <= EXACT_ANGLE, "{corner:?}");
    assert_sound_shell(&model, &network, 12);
}

/// The rounded box corner of [`filleted_corner`], each edge's fillet one
/// patch: its side along the octant runs from `z` to the side's middle and
/// on to `y`, two curves meeting at a junction, where the octant's curve
/// into its corner ends (curves `9 + 3k` to `11 + 3k` for the fillet `k`:
/// its edges back from `y` and from `z`, and the arc across its far end,
/// which carries the junction on).
fn filleted_corner_in_one() -> (Vec<CurveDef>, Vec<ChainedPatch>) {
    let (mut curves, octant) = rounded_corner();
    let mut patches: Vec<ChainedPatch> = octant.into_iter().map(ChainedPatch::from).collect();
    let back = -FILLET_RUN;
    let (y, z) = (Point::new(0.0, BALL, 0.0), Point::new(0.0, 0.0, BALL));
    let at_back = |p: Point| Point::new(back, p.y, p.z);
    let fillet = [
        vec![y, at_back(y)],
        vec![z, at_back(z)],
        (0..=64)
            .map(|i| {
                let a = FRAC_PI_2 * f64::from(i) / 64.0;
                Point::new(back, BALL * a.cos(), BALL * a.sin())
            })
            .collect::<Vec<_>>(),
    ];
    for (k, &(from_y, from_z)) in [(2, 5), (4, 1), (0, 3)].iter().enumerate() {
        let base = curves.len();
        for points in &fillet {
            curves.push(CurveDef::Polyline(
                points.iter().map(|&p| turned(p, k)).collect(),
            ));
        }
        patches.push(ChainedPatch::four([
            vec![(from_z, false), (from_y, true)],
            vec![(base, false)],
            vec![(base + 2, false)],
            vec![(base + 1, true)],
        ]));
    }
    (curves, patches)
}

/// The walk shares at which side `k` of `surface` is split into its
/// curves: the knots of multiplicity three inside the side's direction,
/// in walking order, with the ends.
fn side_shares(surface: &BSplineSurface, k: usize) -> Vec<f64> {
    let knots = if k.is_multiple_of(2) {
        surface.u_knots().knots()
    } else {
        surface.v_knots().knots()
    };
    let mut inner: Vec<f64> = knots
        .iter()
        .copied()
        .filter(|&x| x > 0.0 && x < 1.0 && knots.iter().filter(|&&y| y == x).count() == 3)
        .collect();
    inner.dedup();
    let mut shares = vec![0.0];
    if k < 2 {
        shares.extend(inner);
    } else {
        shares.extend(inner.iter().rev().map(|x| 1.0 - x));
    }
    shares.push(1.0);
    shares
}

/// A patch a curve bounds: the patch, its side, the share of the side's
/// walk the curve covers, and whether it runs against the walk.
type ChainedUser = (usize, usize, (f64, f64), bool);

/// As [`across_shared`], over patches whose sides may be chains of
/// curves: each curve read off either patch over its share of the side,
/// `parts` how many of the side's strips each curve crosses where more
/// than one (nought or missing for one).
fn across_chained(network: &Network, patches: &[ChainedPatch], parts: &[usize]) -> (f64, f64) {
    let crosses = |c: usize| parts.get(c).copied().unwrap_or(1).max(1);
    let mut users: Vec<Vec<ChainedUser>> = Vec::new();
    for (p, patch) in patches.iter().enumerate() {
        for (k, side) in patch.sides.iter().enumerate() {
            let shares = if side.len() > 1 {
                side_shares(&network.surfaces[p], k)
            } else {
                vec![0.0, 1.0]
            };
            let strips: usize = side.iter().map(|s| crosses(s.curve)).sum();
            if side.len() > 1 {
                assert_eq!(shares.len(), strips + 1, "patch {p} side {k}");
            }
            let mut at = 0;
            for piece in side {
                if users.len() <= piece.curve {
                    users.resize(piece.curve + 1, Vec::new());
                }
                let end = if side.len() > 1 {
                    at + crosses(piece.curve)
                } else {
                    1
                };
                users[piece.curve].push((p, k, (shares[at], shares[end]), piece.reversed));
                at = end;
            }
        }
    }
    let (mut angle, mut gap) = (0.0_f64, 0.0_f64);
    for pair in users.iter().filter(|u| u.len() == 2) {
        let at = |(p, k, (a, b), reversed): ChainedUser, t: f64| {
            let f = (b - a).mul_add(if reversed { 1.0 - t } else { t }, a);
            let (u, v) = on_side(k, f);
            let surface = &network.surfaces[p];
            (
                surface.normal_at(u, v, T).unwrap().vector(),
                surface.point_at(u, v, T).unwrap(),
            )
        };
        for i in 0..=200 {
            let t = f64::from(i) / 200.0;
            let ((na, pa), (nb, pb)) = (at(pair[0], t), at(pair[1], t));
            angle = angle.max(na.cross(nb).magnitude().atan2(na.dot(nb)));
            gap = gap.max(pa.distance(pb));
        }
    }
    (angle, gap)
}

/// The widest angle between the normals either side of the junction line
/// across patch `p`'s `u` direction, read at 201 points along it.
fn across_junction_line(network: &Network, p: usize) -> f64 {
    let surface = &network.surfaces[p];
    let shares = side_shares(surface, 0);
    let mut widest = 0.0_f64;
    for &junction in &shares[1..shares.len() - 1] {
        for i in 0..=200 {
            let v = f64::from(i) / 200.0;
            let before = surface.normal_at(junction.next_down(), v, T).unwrap();
            let after = surface.normal_at(junction, v, T).unwrap();
            let (a, b) = (before.vector(), after.vector());
            widest = widest.max(a.cross(b).magnitude().atan2(a.dot(b)));
        }
    }
    widest
}

#[test]
fn a_box_corner_with_its_fillets_one_patch_each_is_held_exactly_tangent() {
    // Each curve into the rounded corner ends at a junction on a fillet's
    // side, which pins the curve's ratio to the fillet's knot spacings
    // either side of it, and the corner pins it to its own: the line across
    // the fillet through the junction is held at the corner's ratio
    // instead, tangent though not `C¹`.
    let (curves, patches) = filleted_corner_in_one();
    for noise in [NOISE, NOISE / 4.0] {
        let tolerance = 5.0 * noise;
        let cloud = filleted_points(40_000, noise, 37);
        let index = PointIndex::build(&cloud).unwrap();
        let mut model = Model::new();
        let network = fit_network_chained(
            &mut model,
            &patches,
            &curves,
            &index,
            &options(tolerance, vec![]),
            0,
            T,
        )
        .unwrap();
        assert_eq!(network.ended, Ended::Met, "noise {noise}: {network:?}");
        let off = widest_off(&network, off_filleted);
        assert!(
            off <= tolerance / 2.0,
            "noise {noise}: {off} off the corner"
        );
        let (angle, gap) = across_chained(&network, &patches, &[]);
        assert!(angle <= EXACT_ANGLE, "noise {noise}: {angle}");
        assert!(gap <= ROUNDING_GAP, "noise {noise}: {gap}");
        for p in 3..6 {
            let line = across_junction_line(&network, p);
            assert!(
                line <= EXACT_ANGLE,
                "noise {noise}: fillet {p} turns {line}"
            );
        }
        assert!(
            network.worst_tangent <= EXACT_ANGLE,
            "noise {noise}: {}",
            network.worst_tangent
        );
        let [corner] = network.corners.as_slice() else {
            panic!("{:?}", network.corners);
        };
        assert_eq!(corner.valence, 3);
        assert!(corner.angle <= EXACT_ANGLE, "{corner:?}");
        // The far ends of the fillets: two edges along each and the arc
        // across its end.
        assert_sound_shell(&model, &network, 9);
    }
}

/// A cap of the ball in four patches with two corners of three, `A` and
/// `B` either side of the top, joined by curve 0 over it: two more curves
/// from each run out to the border (1 and 2 from `B`, 3 and 4 from `A`),
/// a hexagon of great circle arcs with its corners on the cap's rim (5 to
/// 10), so every patch turns at each of its corners.
fn two_corner_cap() -> (Vec<CurveDef>, Vec<PatchDef>) {
    let (a, b) = (on_ball_at(0.3, PI), on_ball_at(0.3, 0.0));
    let rim_at = |azimuth: f64| on_ball_at(CAP, azimuth);
    let rim = |from: f64, to: f64| ball_arc(rim_at(from), rim_at(to));
    let q = PI / 4.0;
    let curves = vec![
        ball_arc(a, b),
        ball_arc(b, rim_at(q)),
        ball_arc(b, rim_at(-q)),
        ball_arc(a, rim_at(3.0 * q)),
        ball_arc(a, rim_at(5.0 * q)),
        rim(q, 3.0 * q),
        rim(5.0 * q, 7.0 * q),
        rim(3.0 * q, PI),
        rim(PI, 5.0 * q),
        rim(-q, 0.0),
        rim(0.0, q),
    ];
    let patches = vec![
        PatchDef::four([(0, false), (1, false), (5, false), (3, true)]),
        PatchDef::four([(0, true), (4, false), (6, false), (2, true)]),
        PatchDef::four([(3, false), (7, false), (8, false), (4, true)]),
        PatchDef::four([(2, false), (9, false), (10, false), (1, true)]),
    ];
    (curves, patches)
}

#[test]
fn two_corners_joined_by_a_curve_are_held_exactly_tangent() {
    // The curve between the corners is pinned to a ratio at either end:
    // both corners are rescaled to the geometric mean of the two, each
    // through its curves that run on to the border. At the lower noise the
    // corners held by their spans miss the points beside them, and the
    // gluing is spread along the curves instead, `γ` running along the
    // curve between them from one corner's value to the other's.
    let (curves, patches) = two_corner_cap();
    for noise in [NOISE, NOISE / 4.0] {
        let tolerance = 5.0 * noise;
        let cloud = cap_points(20_000, noise, 29);
        let (model, network) = fit(&curves, &patches, &cloud, &options(tolerance, vec![]));
        assert_eq!(network.ended, Ended::Met, "{noise}: {network:?}");
        let off = widest_off(&network, off_ball);
        assert!(off <= tolerance / 2.0, "{noise}: {off} off the ball");
        let (angle, gap) = across_shared(&network, &patches);
        assert!(angle <= EXACT_ANGLE, "{noise}: {angle}");
        assert!(gap <= ROUNDING_GAP, "{noise}: {gap}");
        assert!(
            network.worst_tangent <= EXACT_ANGLE,
            "{noise}: {}",
            network.worst_tangent
        );
        assert_eq!(network.corners.len(), 2, "{:?}", network.corners);
        for corner in &network.corners {
            assert_eq!(corner.valence, 3);
            assert!(corner.angle <= EXACT_ANGLE, "{corner:?}");
            // The corners may slide along the ball from where the curves
            // were drawn to meet, but stay on it.
            assert!(off_ball(corner.point) <= tolerance, "{corner:?}");
        }
        assert_sound_shell(&model, &network, 6);
    }
}

/// The ball covered as a cube: a patch per face of the cube in `x, y, z`
/// from −1 to 1, pushed out on to the ball, its curves the great circle
/// arcs over the cube's edges. Every corner is one of three patches, and
/// every curve runs between two of them. Each patch's parameters run up
/// the axes along its face, and each curve up its axis, so no strip turns.
fn ball_cube() -> (Vec<CurveDef>, Vec<PatchDef>) {
    // A cube corner's index: a bit per axis, set where it is at +1.
    let corner = |i: usize| {
        let at = |bit: usize| if i & (1 << bit) == 0 { -1.0 } else { 1.0 };
        on_ball(at(0), at(1), at(2))
    };
    let mut ends: Vec<(usize, usize)> = Vec::new();
    let mut curves = Vec::new();
    let mut patches = Vec::new();
    for axis in 0..3 {
        for high in [false, true] {
            // The face's other two axes, in the order that turns its walk
            // counter-clockwise seen from outside: `e_b × e_c` is `e_axis`
            // taken cyclically.
            let (mut b, mut c) = ((axis + 1) % 3, (axis + 2) % 3);
            if !high {
                (b, c) = (c, b);
            }
            let base = if high { 1 << axis } else { 0 };
            let walk = [0, 1 << b, (1 << b) | (1 << c), 1 << c].map(|x| x | base);
            let mut sides = [(0, false); 4];
            for (k, side) in sides.iter_mut().enumerate() {
                let (from, to) = (walk[k], walk[(k + 1) % 4]);
                let key = (from.min(to), from.max(to));
                let c = ends.iter().position(|&e| e == key).unwrap_or_else(|| {
                    ends.push(key);
                    curves.push(ball_arc(corner(key.0), corner(key.1)));
                    ends.len() - 1
                });
                *side = (c, from > to);
            }
            patches.push(PatchDef::four(sides));
        }
    }
    (curves, patches)
}

/// `count` points over the whole ball, off it radially by noise of
/// deviation `noise`.
fn ball_points(count: usize, noise: f64, seed: u64) -> Vec<Point> {
    let mut draws = Draws(seed);
    (0..count)
        .map(|_| {
            let (a, b, c) = (draws.gaussian(), draws.gaussian(), draws.gaussian());
            let n = (a * a + b * b + c * c).sqrt();
            let r = noise.mul_add(draws.gaussian(), BALL);
            Point::new(r * a / n, r * b / n, r * c / n)
        })
        .collect()
}

#[test]
fn a_ball_covered_as_a_cube_is_held_exactly_tangent_at_every_corner() {
    // Every curve is pinned to a ratio at both its ends, and every curve of
    // every corner is asked one: the corners' scales and the curves'
    // ratios are settled together.
    let (curves, patches) = ball_cube();
    for noise in [NOISE, NOISE / 4.0] {
        let tolerance = 5.0 * noise;
        let cloud = ball_points(40_000, noise, 41);
        let (model, network) = fit(&curves, &patches, &cloud, &options(tolerance, vec![]));
        assert_eq!(network.ended, Ended::Met, "{network:?}");
        let off = widest_off(&network, off_ball);
        assert!(off <= tolerance / 2.0, "{off} off the ball");
        let (angle, gap) = across_shared(&network, &patches);
        assert!(angle <= EXACT_ANGLE, "{angle}");
        assert!(gap <= ROUNDING_GAP, "{gap}");
        assert!(
            network.worst_tangent <= EXACT_ANGLE,
            "{}",
            network.worst_tangent
        );
        assert_eq!(network.corners.len(), 8, "{:?}", network.corners);
        for corner in &network.corners {
            assert_eq!(corner.valence, 3);
            assert!(corner.angle <= EXACT_ANGLE, "{corner:?}");
            assert!(off_ball(corner.point) <= tolerance, "{corner:?}");
        }
        // Closed: every edge used twice.
        let diagnosis = check(&model, &network.shell, T).unwrap();
        assert!(diagnosis.problems.is_empty(), "{diagnosis:?}");
    }
}

#[test]
fn a_clean_ball_covered_as_a_cube_meets_its_tolerance_exactly_tangent() {
    // The seeds take a few points by the patch beside the face of the cube
    // they lie over, up to 0.0022 past the great circle between the two,
    // and every fit at the seed's parameters leaves them about as far past
    // that patch's side: no try meets 0.001 until a held fit is finished,
    // those points taken by the patch across the side.
    let (curves, patches) = ball_cube();
    let cloud = ball_points(40_000, 0.0, 41);
    let (model, network) = fit(&curves, &patches, &cloud, &options(CLEAN_TOLERANCE, vec![]));
    assert_eq!(network.ended, Ended::Met, "{network:?}");
    assert_eq!(network.outside, 0);
    assert_eq!(
        network.per_patch.iter().map(|p| p.points).sum::<usize>(),
        cloud.len()
    );
    assert!(
        network.per_patch.iter().all(|p| p.worst <= CLEAN_TOLERANCE),
        "{:?}",
        network.per_patch
    );
    let off = widest_off(&network, off_ball);
    assert!(off <= CLEAN_TOLERANCE, "{off} off the ball");
    let (angle, gap) = across_shared(&network, &patches);
    assert!(angle <= EXACT_ANGLE, "{angle}");
    assert!(gap <= ROUNDING_GAP, "{gap}");
    assert!(
        network.worst_tangent <= EXACT_ANGLE,
        "{}",
        network.worst_tangent
    );
    assert_eq!(network.corners.len(), 8, "{:?}", network.corners);
    for corner in &network.corners {
        assert!(corner.angle <= EXACT_ANGLE, "{corner:?}");
        assert!(off_ball(corner.point) <= CLEAN_TOLERANCE, "{corner:?}");
    }
    let diagnosis = check(&model, &network.shell, T).unwrap();
    assert!(diagnosis.problems.is_empty(), "{diagnosis:?}");
}

/// A cap with a corner of six patches at its pole, one of its curves
/// running on to a corner of three between the pole and the rim, at
/// azimuth nought: curve 0 joins the two, curves 1 and 2 run from the
/// second corner to the rim at azimuths `±π/6`, and curves 3 to 7 from the
/// pole to the rim at azimuths `π/3` to `5π/3`.
fn even_joined_cap() -> (Vec<CurveDef>, Vec<PatchDef>) {
    let (a, b) = (on_ball(0.0, 0.0, 1.0), on_ball_at(0.5, 0.0));
    let rim_at = |azimuth: f64| on_ball_at(CAP, azimuth);
    let sixth = PI / 3.0;
    let mut curves = vec![
        ball_arc(a, b),
        ball_arc(b, rim_at(sixth / 2.0)),
        ball_arc(b, rim_at(-sixth / 2.0)),
    ];
    curves.extend((1..6).map(|j| ball_arc(a, rim_at(sixth * j as f64))));
    let mut rim = |from: f64, to: f64| {
        curves.push(ball_arc(rim_at(from), rim_at(to)));
        curves.len() - 1
    };
    let mut patches = vec![
        PatchDef::four([
            (0, false),
            (1, false),
            (rim(sixth / 2.0, sixth), false),
            (3, true),
        ]),
        PatchDef::four([
            (0, true),
            (7, false),
            (rim(5.0 * sixth, 5.5 * sixth), false),
            (2, true),
        ]),
    ];
    for j in 1..5 {
        let from = sixth * j as f64;
        patches.push(PatchDef::four([
            (2 + j, false),
            (rim(from, from + sixth / 2.0), false),
            (rim(from + sixth / 2.0, from + sixth), false),
            (3 + j, true),
        ]));
    }
    patches.push(PatchDef::four([
        (2, false),
        (rim(-sixth / 2.0, 0.0), false),
        (rim(0.0, sixth / 2.0), false),
        (1, true),
    ]));
    (curves, patches)
}

#[test]
fn an_even_corner_joined_to_an_odd_one_is_held_exactly_tangent() {
    // A corner of six patches joined by a curve to one of three: the curve
    // is pinned to a ratio at either end, and the even corner is rescaled
    // to their mean through the curves of its own cycle of scales. Held by
    // its spans, an even corner leaves points about 0.09 off at the seed's
    // parameters on this ball (a cap of six patches 0.097), so the noise is
    // taken a half wider than elsewhere, the tolerance five deviations.
    let (curves, patches) = even_joined_cap();
    let noise = 1.5 * NOISE;
    let tolerance = 5.0 * noise;
    let cloud = cap_points(20_000, noise, 29);
    let (model, network) = fit(&curves, &patches, &cloud, &options(tolerance, vec![]));
    assert_eq!(network.ended, Ended::Met, "{network:?}");
    let off = widest_off(&network, off_ball);
    assert!(off <= tolerance / 2.0, "{off} off the ball");
    let (angle, gap) = across_shared(&network, &patches);
    assert!(angle <= EXACT_ANGLE, "{angle}");
    assert!(gap <= ROUNDING_GAP, "{gap}");
    assert!(
        network.worst_tangent <= EXACT_ANGLE,
        "{}",
        network.worst_tangent
    );
    let mut valences: Vec<usize> = network.corners.iter().map(|c| c.valence).collect();
    valences.sort_unstable();
    assert_eq!(valences, [3, 6], "{:?}", network.corners);
    for corner in &network.corners {
        assert!(corner.angle <= EXACT_ANGLE, "{corner:?}");
        assert!(off_ball(corner.point) <= tolerance, "{corner:?}");
    }
    assert_sound_shell(&model, &network, 2 * patches.len() - 2);
}
