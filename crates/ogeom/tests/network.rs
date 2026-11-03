//! A network of patches fitted to a scan at once: a cylinder band in
//! sectors and rows, a box's edge rounded by a fillet beside a sharp one,
//! the corner of a box rounded on all three edges, and a three-sided
//! patch, each against the shape it samples, known in closed form.
#![allow(clippy::unwrap_used, reason = "test code")]
#![allow(clippy::cast_precision_loss, reason = "small counts")]

use core::f64::consts::{FRAC_PI_2, PI, TAU};
use std::time::Duration;

use ogeom::algo::{CurveDef, Network, NetworkOptions, PatchDef, PatchSide, check, fit_network};
use ogeom::core::clock::Instant;
use ogeom::core::{OgeomError, Tolerances};
use ogeom::geom::fit::Ended;
use ogeom::geom::{CircleCurve, Curve, Surface as _};
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

#[test]
fn a_corner_of_three_patches_is_reported_with_its_angle() {
    let (curves, patches) = rounded_corner();
    let middle = on_ball(1.0, 1.0, 1.0);
    for (noise, tolerance) in [(0.0, CLEAN_TOLERANCE), (NOISE, NOISY_TOLERANCE)] {
        let cloud = octant_points(20_000, noise, 3);
        let (model, network) = fit(&curves, &patches, &cloud, &options(tolerance, vec![]));
        assert!(widest_off(&network, off_ball) <= tolerance);
        assert!(network.worst_gap <= ROUNDING_GAP, "{}", network.worst_gap);
        let [corner] = network.corners.as_slice() else {
            panic!("{:?}", network.corners);
        };
        assert_eq!(corner.valence, 3);
        assert!(corner.point.distance(middle) <= tolerance, "{corner:?}");
        // Near tangent: a patch's normal near a curve is that of its first
        // poles in from the curve, which a fit within the tolerance tilts
        // by at most twice the tolerance over their spacing, a third of a
        // span (the sides are arcs of about 0.6155 radians, the arc from
        // an octant's side's middle to its middle); the two patches either
        // side tilt opposite ways at worst.
        let spans = network
            .per_patch
            .iter()
            .map(|p| p.spans.0.max(p.spans.1))
            .max()
            .unwrap();
        let spacing = BALL * 0.6155 / (3.0 * spans as f64);
        let bound = 4.0 * tolerance / spacing;
        assert!(corner.angle > 0.0, "{corner:?}");
        assert!(corner.angle <= bound, "{corner:?} beyond {bound}");
        assert!(network.worst_tangent >= corner.angle);
        assert!(
            network.worst_tangent <= bound,
            "{} beyond {bound}",
            network.worst_tangent
        );
        assert_sound_shell(&model, &network, 6);
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
