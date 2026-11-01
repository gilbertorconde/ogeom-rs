//! Best fits of primitives to noisy points, against the shapes the points
//! were drawn from: ten thousand points with seeded Gaussian noise of
//! 0.02 mm land every parameter within four standard errors (4σ/√n) of the
//! truth and the reported sigma within 5 % of the noise, with or without
//! 2 % of the points thrown a millimetre off and a 3σ cutoff.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::cast_precision_loss,
    reason = "test code"
)]

use ogeom_algo::Cutoff;
use ogeom_algo::fit::{self, Fit, FitFixed, FitOptions, Side};
use ogeom_core::{OgeomError, Tolerances};
use ogeom_math::{Axis, Direction, Frame, Point, Point2, Vector, Vector2};

const T: Tolerances = Tolerances::millimetres();
const SIGMA: f64 = 0.02;
const COUNT: usize = 10_000;
/// One point in fifty is an outlier.
const OUTLIER_EVERY: usize = 50;
const OUTLIER: f64 = 1.0;

/// A fixed sequence of draws (xorshift), so every run sees the same points.
struct Seeded(u64);

impl Seeded {
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

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    /// A standard normal draw, by Box-Muller.
    fn gaussian(&mut self) -> f64 {
        let u = 1.0 - self.unit();
        let v = self.unit();
        (-2.0 * u.ln()).sqrt() * (core::f64::consts::TAU * v).cos()
    }
}

/// A sample on the true shape and the unit directions across it the noise
/// moves it along: one for a surface, two for a curve in space.
struct Sample {
    at: Point,
    across: Vector,
    also: Option<Vector>,
}

/// The samples moved by Gaussian noise whose orthogonal RMS is `SIGMA`, and
/// with `outliers` every fiftieth also thrown `OUTLIER` along `across`.
fn noisy(samples: &[Sample], seed: u64, outliers: bool) -> Vec<Point> {
    let mut draw = Seeded(seed);
    samples
        .iter()
        .enumerate()
        .map(|(i, s)| {
            let mut p = match s.also {
                None => s.at + s.across * (SIGMA * draw.gaussian()),
                Some(v) => {
                    let each = SIGMA / 2.0_f64.sqrt();
                    s.at + s.across * (each * draw.gaussian()) + v * (each * draw.gaussian())
                }
            };
            if outliers && i % OUTLIER_EVERY == 0 {
                p += s.across * OUTLIER;
            }
            p
        })
        .collect()
}

fn tilted() -> Frame {
    Frame::new(
        Point::new(12.0, -7.0, 31.0),
        Direction::new(Vector::new(0.3, -0.4, 1.0), T).unwrap(),
        Direction::new(Vector::new(1.0, 0.5, 0.05), T).unwrap(),
        T,
    )
    .unwrap()
}

/// 4σ/√n: four standard errors of a parameter that moves each point's
/// distance by as much as it moves itself.
fn bound() -> f64 {
    4.0 * SIGMA / (COUNT as f64).sqrt()
}

fn near(what: &str, error: f64, variant: &str) {
    near_by(what, error, 1.0, variant);
}

/// Within four standard errors of a parameter that moves a point's
/// distance by, on average and squared, `share` of its own change: a
/// sphere's centre moves a point's distance by its component along the
/// point's direction, a third of the whole on average over the sphere, and
/// is determined to `σ/√(n/3)`. Where the share is one, the bound is
/// 4σ/√n.
fn near_by(what: &str, error: f64, share: f64, variant: &str) {
    let within = bound() / share.sqrt();
    assert!(
        error.abs() <= within,
        "{variant}: {what} off by {error:e}, more than 4σ/√(n·{share}) = {within:e}"
    );
}

/// The fit's sigma within 5 % of the noise, and with outliers, none of them
/// kept.
fn fits_well<S>(fit: &Fit<S>, variant: &str, outliers: bool) {
    assert!(
        (fit.sigma / SIGMA - 1.0).abs() <= 0.05,
        "{variant}: sigma {} against noise {SIGMA}",
        fit.sigma
    );
    assert!(fit.converged, "{variant}: not converged");
    assert_eq!(fit.residuals.len(), COUNT);
    if outliers {
        for i in (0..COUNT).step_by(OUTLIER_EVERY) {
            assert!(!fit.used[i], "{variant}: outlier {i} kept");
        }
    }
}

/// The three ways each kind is fitted: clean points, outliers with a 3σ
/// cutoff, and outliers with a 3σ cutoff from a consensus start.
fn variants<S>() -> [(&'static str, bool, FitOptions<'static, S>); 3] {
    [
        ("clean", false, FitOptions::default()),
        (
            "outliers",
            true,
            FitOptions {
                cutoff: Cutoff::Sigma(3.0),
                ..FitOptions::default()
            },
        ),
        (
            "outliers, robust start",
            true,
            FitOptions {
                cutoff: Cutoff::Sigma(3.0),
                robust_start: true,
                ..FitOptions::default()
            },
        ),
    ]
}

/// The standard deviation of the samples' coordinate along `d` about `o`:
/// the lever an axis tilt toward `d` moves the points by.
fn lever(samples: &[Sample], o: Point, d: Vector) -> f64 {
    let n = samples.len() as f64;
    let mean = samples.iter().map(|s| (s.at - o).dot(d)).sum::<f64>() / n;
    (samples
        .iter()
        .map(|s| ((s.at - o).dot(d) - mean).powi(2))
        .sum::<f64>()
        / n)
        .sqrt()
}

/// How far a fitted axis direction tilts from the true one toward each of
/// the frame's `x` and `y`, the sense of the fitted one free.
fn tilts(frame: Frame, fitted: Direction) -> (f64, f64) {
    let v = if fitted.vector().dot(frame.z().vector()) < 0.0 {
        -fitted.vector()
    } else {
        fitted.vector()
    };
    (v.dot(frame.x().vector()), v.dot(frame.y().vector()))
}

/// Where a fitted axis crosses the true axis's plane through `at`, in the
/// true frame's `x` and `y`: how far across the true axis it passes.
fn crossing(frame: Frame, at: Point, location: Point, direction: Direction) -> (f64, f64) {
    let z = frame.z().vector();
    let d = direction.vector();
    let t = (at - location).dot(z) / d.dot(z);
    let hit = location + d * t - at;
    (hit.dot(frame.x().vector()), hit.dot(frame.y().vector()))
}

#[test]
fn a_plane_lands_within_four_standard_errors() {
    let frame = tilted();
    let mut draw = Seeded(11);
    let samples: Vec<Sample> = (0..COUNT)
        .map(|_| {
            let (x, y) = (draw.range(-20.0, 20.0), draw.range(-15.0, 15.0));
            Sample {
                at: frame.origin() + frame.x().vector() * x + frame.y().vector() * y,
                across: frame.z().vector(),
                also: None,
            }
        })
        .collect();
    let lx = lever(&samples, frame.origin(), frame.x().vector());
    let ly = lever(&samples, frame.origin(), frame.y().vector());
    for (variant, outliers, options) in variants() {
        let points = noisy(&samples, 12, outliers);
        let fit = fit::plane(&points, &options, T).unwrap();
        fits_well(&fit, variant, outliers);
        let (tx, ty) = tilts(frame, fit.shape.normal());
        near("tilt toward x", tx * lx, variant);
        near("tilt toward y", ty * ly, variant);
        near(
            "offset",
            fit.shape.signed_distance_to(frame.origin()),
            variant,
        );
    }
}

#[test]
fn a_sphere_lands_within_four_standard_errors() {
    let (centre, radius) = (Point::new(3.0, -2.0, 5.0), 12.0);
    let mut draw = Seeded(21);
    let samples: Vec<Sample> = (0..COUNT)
        .map(|_| {
            let z = draw.range(-1.0, 1.0);
            let a = draw.range(0.0, core::f64::consts::TAU);
            let s = (1.0 - z * z).sqrt();
            let n = Vector::new(s * a.cos(), s * a.sin(), z);
            Sample {
                at: centre + n * radius,
                across: n,
                also: None,
            }
        })
        .collect();
    for (variant, outliers, options) in variants() {
        let points = noisy(&samples, 22, outliers);
        let fit = fit::sphere(&points, &options, T).unwrap();
        fits_well(&fit, variant, outliers);
        let c = fit.shape.centre();
        near_by("centre x", c.x - centre.x, 1.0 / 3.0, variant);
        near_by("centre y", c.y - centre.y, 1.0 / 3.0, variant);
        near_by("centre z", c.z - centre.z, 1.0 / 3.0, variant);
        near("radius", fit.shape.radius() - radius, variant);
    }
}

fn cylinder_samples(frame: Frame, radius: f64, arc: f64, length: f64, seed: u64) -> Vec<Sample> {
    let mut draw = Seeded(seed);
    (0..COUNT)
        .map(|_| {
            let u = draw.range(0.0, arc);
            let h = draw.range(-length / 2.0, length / 2.0);
            let n = frame.x().vector() * u.cos() + frame.y().vector() * u.sin();
            Sample {
                at: frame.origin() + n * radius + frame.z().vector() * h,
                across: n,
                also: None,
            }
        })
        .collect()
}

#[test]
fn a_cylinder_lands_within_four_standard_errors() {
    let (frame, radius) = (tilted(), 8.0);
    let samples = cylinder_samples(frame, radius, core::f64::consts::TAU, 30.0, 31);
    let lx = lever(&samples, frame.origin(), frame.z().vector());
    for (variant, outliers, options) in variants() {
        let points = noisy(&samples, 32, outliers);
        let fit = fit::cylinder(&points, &options, T).unwrap();
        fits_well(&fit, variant, outliers);
        let f = fit.shape.frame();
        near("radius", fit.shape.radius() - radius, variant);
        let (tx, ty) = tilts(frame, f.z());
        // Each moves a point's distance by its share along the point's own
        // direction round the axis, half on average.
        near_by("tilt toward x", tx * lx, 0.5, variant);
        near_by("tilt toward y", ty * lx, 0.5, variant);
        let (cx, cy) = crossing(frame, frame.origin(), f.origin(), f.z());
        near_by("axis across x", cx, 0.5, variant);
        near_by("axis across y", cy, 0.5, variant);
    }
}

#[test]
fn a_cone_lands_within_four_standard_errors() {
    let frame = tilted();
    let (reference, angle) = (6.0, 20.0_f64.to_radians());
    let (s, c) = angle.sin_cos();
    let mut draw = Seeded(41);
    let samples: Vec<Sample> = (0..COUNT)
        .map(|_| {
            let u = draw.range(0.0, core::f64::consts::TAU);
            let h = draw.range(-10.0, 10.0);
            let radial = frame.x().vector() * u.cos() + frame.y().vector() * u.sin();
            Sample {
                at: frame.origin()
                    + radial * angle.tan().mul_add(h, reference)
                    + frame.z().vector() * h,
                across: radial * c - frame.z().vector() * s,
                also: None,
            }
        })
        .collect();
    let lz = lever(&samples, frame.origin(), frame.z().vector());
    for (variant, outliers, options) in variants() {
        let points = noisy(&samples, 42, outliers);
        let fit = fit::cone(&points, &options, T).unwrap();
        fits_well(&fit, variant, outliers);
        let f = fit.shape.frame();
        let (tx, ty) = tilts(frame, f.z());
        near("tilt toward x", tx * lz, variant);
        near("tilt toward y", ty * lz, variant);
        let (cx, cy) = crossing(frame, frame.origin(), f.origin(), f.z());
        near("axis across x", cx, variant);
        near("axis across y", cy, variant);
        let z = (frame.origin() - f.origin()).dot(f.z().vector());
        near(
            "radius at the true reference",
            (fit.shape.radius_at(z) - reference) * c,
            variant,
        );
        near("half angle", (fit.shape.half_angle() - angle) * lz, variant);
    }
}

#[test]
fn a_torus_lands_within_four_standard_errors() {
    let (frame, major, minor) = (tilted(), 20.0, 5.0);
    let mut draw = Seeded(51);
    let samples: Vec<Sample> = (0..COUNT)
        .map(|_| {
            let u = draw.range(0.0, core::f64::consts::TAU);
            let v = draw.range(0.0, core::f64::consts::TAU);
            let radial = frame.x().vector() * u.cos() + frame.y().vector() * u.sin();
            let n = radial * v.cos() + frame.z().vector() * v.sin();
            Sample {
                at: frame.origin() + radial * major + n * minor,
                across: n,
                also: None,
            }
        })
        .collect();
    let lever_tilt = lever(&samples, frame.origin(), frame.x().vector());
    for (variant, outliers, options) in variants() {
        let points = noisy(&samples, 52, outliers);
        let fit = fit::torus(&points, &options, T).unwrap();
        fits_well(&fit, variant, outliers);
        let t = &fit.shape;
        near("minor radius", t.minor_radius() - minor, variant);
        // The major radius and the centre's height move a point's distance
        // by the cosine and the sine of its angle round the tube, half on
        // average; the centre across the axis by those times the cosine
        // round the axis, a quarter.
        near_by("major radius", t.major_radius() - major, 0.5, variant);
        let off = t.centre() - frame.origin();
        near_by("centre x", off.dot(frame.x().vector()), 0.25, variant);
        near_by("centre y", off.dot(frame.y().vector()), 0.25, variant);
        near_by("centre z", off.dot(frame.z().vector()), 0.5, variant);
        let (tx, ty) = tilts(frame, t.frame().z());
        near("tilt toward x", tx * lever_tilt, variant);
        near("tilt toward y", ty * lever_tilt, variant);
    }
}

#[test]
fn a_circle_in_space_lands_within_four_standard_errors() {
    let (frame, radius) = (tilted(), 15.0);
    let mut draw = Seeded(61);
    let samples: Vec<Sample> = (0..COUNT)
        .map(|_| {
            let u = draw.range(0.0, core::f64::consts::TAU);
            let n = frame.x().vector() * u.cos() + frame.y().vector() * u.sin();
            Sample {
                at: frame.origin() + n * radius,
                across: n,
                also: Some(frame.z().vector()),
            }
        })
        .collect();
    let l = lever(&samples, frame.origin(), frame.x().vector());
    for (variant, outliers, options) in variants() {
        let points = noisy(&samples, 62, outliers);
        let fit = fit::circle3(&points, &options, T).unwrap();
        fits_well(&fit, variant, outliers);
        near("radius", fit.shape.radius() - radius, variant);
        let off = fit.shape.centre() - frame.origin();
        near("centre x", off.dot(frame.x().vector()), variant);
        near("centre y", off.dot(frame.y().vector()), variant);
        near("centre z", off.dot(frame.z().vector()), variant);
        let (tx, ty) = tilts(frame, fit.shape.frame().z());
        near("tilt toward x", tx * l, variant);
        near("tilt toward y", ty * l, variant);
    }
}

#[test]
fn a_line_in_space_lands_within_four_standard_errors() {
    let frame = tilted();
    let mut draw = Seeded(71);
    let samples: Vec<Sample> = (0..COUNT)
        .map(|_| Sample {
            at: frame.origin() + frame.z().vector() * draw.range(-25.0, 25.0),
            across: frame.x().vector(),
            also: Some(frame.y().vector()),
        })
        .collect();
    let l = lever(&samples, frame.origin(), frame.z().vector());
    for (variant, outliers, options) in variants() {
        let points = noisy(&samples, 72, outliers);
        let fit = fit::line3(&points, &options, T).unwrap();
        fits_well(&fit, variant, outliers);
        let (tx, ty) = tilts(frame, fit.shape.direction);
        near("tilt toward x", tx * l, variant);
        near("tilt toward y", ty * l, variant);
        let (cx, cy) = crossing(
            frame,
            frame.origin(),
            fit.shape.location,
            fit.shape.direction,
        );
        near("across x", cx, variant);
        near("across y", cy, variant);
    }
}

#[test]
fn a_circle_in_the_plane_lands_within_four_standard_errors() {
    let (centre, radius) = (Point2::new(1.0, 2.0), 10.0);
    let mut draw = Seeded(81);
    let samples: Vec<Sample> = (0..COUNT)
        .map(|_| {
            let u = draw.range(0.0, core::f64::consts::TAU);
            let n = Vector::new(u.cos(), u.sin(), 0.0);
            Sample {
                at: centre.to_3d() + n * radius,
                across: n,
                also: None,
            }
        })
        .collect();
    for (variant, outliers, options) in variants() {
        let points: Vec<Point2> = noisy(&samples, 82, outliers)
            .iter()
            .map(|p| p.xy())
            .collect();
        let fit = fit::circle2(&points, &options, T).unwrap();
        fits_well(&fit, variant, outliers);
        near("radius", fit.shape.radius() - radius, variant);
        // Each moves a point's distance by the cosine or the sine of its
        // angle, half on average.
        near_by("centre x", fit.shape.centre().x - centre.x, 0.5, variant);
        near_by("centre y", fit.shape.centre().y - centre.y, 0.5, variant);
    }
}

#[test]
fn a_line_in_the_plane_lands_within_four_standard_errors() {
    let (origin, along) = (Point2::new(-3.0, 4.0), Vector2::new(0.8, 0.6));
    let across = Vector2::new(-0.6, 0.8);
    let mut draw = Seeded(91);
    let samples: Vec<Sample> = (0..COUNT)
        .map(|_| Sample {
            at: (origin.to_vector() + along * draw.range(-30.0, 30.0))
                .to_3d()
                .into_point(),
            across: across.to_3d(),
            also: None,
        })
        .collect();
    let l = lever(&samples, origin.to_3d(), along.to_3d());
    for (variant, outliers, options) in variants() {
        let points: Vec<Point2> = noisy(&samples, 92, outliers)
            .iter()
            .map(|p| p.xy())
            .collect();
        let fit = fit::line2(&points, &options, T).unwrap();
        fits_well(&fit, variant, outliers);
        let d = fit.shape.direction.vector();
        let d = if d.dot(along) < 0.0 { d * -1.0 } else { d };
        near("tilt", d.dot(across) * l, variant);
        near(
            "offset",
            (origin - fit.shape.location).to_3d().dot(across.to_3d()),
            variant,
        );
    }
}

trait IntoPoint {
    fn into_point(self) -> Point;
}

impl IntoPoint for Vector {
    fn into_point(self) -> Point {
        Point::from_vector(self)
    }
}

/// A 60° arc of a cylinder, with no normals: its radius comes back within
/// half a percent.
#[test]
fn a_sixty_degree_arc_of_a_cylinder_recovers_its_radius() {
    let (frame, radius) = (tilted(), 25.0);
    let samples = cylinder_samples(frame, radius, 60.0_f64.to_radians(), 20.0, 101);
    let points = noisy(&samples, 102, false);
    let fit = fit::cylinder(&points, &FitOptions::default(), T).unwrap();
    let error = (fit.shape.radius() / radius - 1.0).abs();
    assert!(
        error <= 0.005,
        "radius {} off by {error}",
        fit.shape.radius()
    );
    assert!((fit.sigma / SIGMA - 1.0).abs() <= 0.05, "{}", fit.sigma);
}

/// A precision ball scanned over its upper half, fitted with its radius
/// held a little short: the radius stays exactly as held, and the centre
/// moves up toward the scanned cap by what the linearized least squares
/// says, `δ·Σcos θ / Σcos² θ`, to second order in the shift.
#[test]
fn a_held_radius_holds_exactly_and_the_centre_moves_for_it() {
    let (radius, short) = (10.0, 0.05);
    let mut draw = Seeded(111);
    let mut points = Vec::new();
    let (mut sc, mut scc) = (0.0, 0.0);
    for _ in 0..COUNT {
        let z = draw.range(0.0, 1.0);
        let a = draw.range(0.0, core::f64::consts::TAU);
        let s = (1.0 - z * z).sqrt();
        points.push(Point::new(
            radius * s * a.cos(),
            radius * s * a.sin(),
            radius * z,
        ));
        sc += z;
        scc += z * z;
    }
    let held = radius - short;
    let options = FitOptions {
        fixed: FitFixed {
            radius: Some(held),
            ..FitFixed::default()
        },
        ..FitOptions::default()
    };
    let fit = fit::sphere(&points, &options, T).unwrap();
    assert_eq!(fit.shape.radius().to_bits(), held.to_bits());
    let shift = short * sc / scc;
    let c = fit.shape.centre();
    assert!(
        (c.z - shift).abs() <= shift * shift / radius,
        "centre moved {} against {shift}",
        c.z
    );
    assert!(c.x.abs() < 1e-3 && c.y.abs() < 1e-3, "{c:?}");
    // Free, the same points give the true ball.
    let free = fit::sphere(&points, &FitOptions::default(), T).unwrap();
    assert!((free.shape.radius() - radius).abs() < 1e-9);
    assert!(free.shape.centre().distance(Point::ORIGIN) < 1e-9);
}

/// What is held stays exactly as given: a plane's normal, a cylinder's
/// axis, a cone's half angle.
#[test]
fn what_is_held_stays_as_given() {
    let frame = tilted();
    let samples = cylinder_samples(frame, 8.0, core::f64::consts::TAU, 30.0, 121);
    let points = noisy(&samples, 122, false);
    let axis = Axis::new(frame.origin(), frame.z());
    let fit = fit::cylinder(
        &points,
        &FitOptions {
            fixed: FitFixed {
                axis: Some(axis),
                ..FitFixed::default()
            },
            ..FitOptions::default()
        },
        T,
    )
    .unwrap();
    let z = fit.shape.frame().z();
    assert!(z == frame.z() || z == -frame.z(), "{z:?}");
    assert!(
        Axis::new(fit.shape.frame().origin(), z).location.distance(
            frame.origin()
                + frame.z().vector()
                    * (fit.shape.frame().origin() - frame.origin()).dot(frame.z().vector())
        ) < 1e-9
    );
    near(
        "radius about the held axis",
        fit.shape.radius() - 8.0,
        "held axis",
    );

    let tipped = Direction::new(Vector::new(0.0, 0.01, 1.0), T).unwrap();
    let plane = fit::plane(
        &points,
        &FitOptions {
            fixed: FitFixed {
                direction: Some(tipped),
                ..FitFixed::default()
            },
            ..FitOptions::default()
        },
        T,
    )
    .unwrap();
    assert!(plane.shape.normal() == tipped || plane.shape.normal() == -tipped);

    let half = 0.3_f64;
    let cone_points: Vec<Point> = (0..400)
        .map(|i| {
            let (u, h) = (f64::from(i % 20) * 0.31, f64::from(i / 20) * 0.5);
            let r = half.tan().mul_add(h, 4.0) + 0.01 * f64::from(i % 3);
            Point::new(r * u.cos(), r * u.sin(), h)
        })
        .collect();
    let cone = fit::cone(
        &cone_points,
        &FitOptions {
            fixed: FitFixed {
                half_angle: Some(half),
                ..FitFixed::default()
            },
            ..FitOptions::default()
        },
        T,
    )
    .unwrap();
    assert_eq!(cone.shape.half_angle().to_bits(), half.to_bits());
}

/// A bore's normals face the axis: the side read from them is the
/// material outside, and a point pushed into the material (away from the
/// axis) has a negative residual. Told the side, the fit signs the same
/// way without normals.
#[test]
fn a_bore_signs_its_residuals_into_the_material() {
    let frame = tilted();
    let samples = cylinder_samples(frame, 8.0, core::f64::consts::TAU, 30.0, 131);
    let mut points: Vec<Point> = samples.iter().map(|s| s.at).collect();
    points[0] += samples[0].across * 0.5;
    let inward: Vec<Vector> = samples.iter().map(|s| s.across * -1.0).collect();
    let read = fit::cylinder(
        &points,
        &FitOptions {
            normals: Some(&inward),
            ..FitOptions::default()
        },
        T,
    )
    .unwrap();
    assert_eq!(read.side, Some(Side::Outside));
    assert!(read.residuals[0] < -0.4, "{}", read.residuals[0]);
    let told = fit::cylinder(
        &points,
        &FitOptions {
            side: Some(Side::Outside),
            ..FitOptions::default()
        },
        T,
    )
    .unwrap();
    assert_eq!(told.side, Some(Side::Outside));
    assert!(told.residuals[0] < -0.4, "{}", told.residuals[0]);
    let boss = fit::cylinder(&points, &FitOptions::default(), T).unwrap();
    assert_eq!(boss.side, None);
    assert!(boss.residuals[0] > 0.4, "{}", boss.residuals[0]);
}

/// Refusals name their reason: too few points, a degenerate spread, a
/// radius outside the range asked for, and an option the kind has no use
/// for.
#[test]
fn refusals_name_their_reason() {
    let message = |e: OgeomError| e.to_string();
    let few = [Point::ORIGIN, Point::new(1.0, 0.0, 0.0)];
    let e = fit::plane(&few, &FitOptions::default(), T).unwrap_err();
    assert!(message(e).contains("needs at least 3 points"));

    let collinear: Vec<Point> = (0..10)
        .map(|i| Point::new(f64::from(i), 0.0, 0.0))
        .collect();
    let e = fit::plane(&collinear, &FitOptions::default(), T).unwrap_err();
    assert!(message(e).contains("degenerate"), "plane");
    let e = fit::circle3(&collinear, &FitOptions::default(), T).unwrap_err();
    assert!(message(e).contains("degenerate"), "circle");
    let flat: Vec<Point> = (0..100)
        .map(|i| Point::new(f64::from(i % 10), f64::from(i / 10), 0.0))
        .collect();
    let e = fit::sphere(&flat, &FitOptions::default(), T).unwrap_err();
    assert!(message(e).contains("degenerate"), "sphere");

    let ring: Vec<Point2> = (0..12)
        .map(|i| {
            let a = f64::from(i) * 0.5;
            Point2::new(5.0 * a.cos(), 5.0 * a.sin())
        })
        .collect();
    let e = fit::circle2(
        &ring,
        &FitOptions {
            radius_range: Some((6.0, 8.0)),
            ..FitOptions::default()
        },
        T,
    )
    .unwrap_err();
    assert!(matches!(e, OgeomError::Range(_)), "{e}");
    assert!(
        fit::circle2(
            &ring,
            &FitOptions {
                radius_range: Some((4.0, 6.0)),
                ..FitOptions::default()
            },
            T,
        )
        .is_ok()
    );

    let e = fit::plane(
        &flat,
        &FitOptions {
            side: Some(Side::Inside),
            ..FitOptions::default()
        },
        T,
    )
    .unwrap_err();
    assert!(message(e).contains("does not apply"));
    let e = fit::cone(
        &flat,
        &FitOptions {
            fixed: FitFixed {
                radius: Some(1.0),
                ..FitFixed::default()
            },
            ..FitOptions::default()
        },
        T,
    )
    .unwrap_err();
    assert!(message(e).contains("does not apply"));
}

/// The same points give the same fit, bit for bit.
#[test]
fn a_fit_is_the_same_every_time() {
    let frame = tilted();
    let samples = cylinder_samples(frame, 8.0, core::f64::consts::TAU, 30.0, 141);
    let points = noisy(&samples, 142, true);
    let options = FitOptions {
        cutoff: Cutoff::Sigma(3.0),
        robust_start: true,
        ..FitOptions::default()
    };
    let first = fit::cylinder(&points, &options, T).unwrap();
    let again = fit::cylinder(&points, &options, T).unwrap();
    assert_eq!(first, again);
}

fn holding<S>(fixed: FitFixed) -> FitOptions<'static, S> {
    FitOptions {
        fixed,
        ..FitOptions::default()
    }
}

/// Exact points on each shape come back as that shape to well within a
/// micron, through every start: normals given, a direction or an axis
/// held, a guess.
#[test]
fn every_start_lands_on_exact_points() {
    let frame = tilted();
    let axis = Axis::new(frame.origin(), frame.z());
    let held_axis = FitFixed {
        axis: Some(axis),
        ..FitFixed::default()
    };
    let held_direction = FitFixed {
        direction: Some(frame.z()),
        ..FitFixed::default()
    };

    let cylinder = cylinder_samples(frame, 8.0, 2.0, 12.0, 151);
    let points: Vec<Point> = cylinder.iter().map(|s| s.at).collect();
    let normals: Vec<Vector> = cylinder.iter().map(|s| s.across).collect();
    for (what, options) in [
        (
            "normals",
            FitOptions {
                normals: Some(&normals),
                ..FitOptions::default()
            },
        ),
        ("direction", holding(held_direction)),
        ("axis", holding(held_axis)),
        (
            "guess",
            FitOptions {
                guess: Some(
                    ogeom_math::Cylinder::new(
                        Frame::about(
                            frame.origin() + Vector::new(0.3, -0.2, 0.1),
                            Direction::new(frame.z().vector() + Vector::new(0.05, 0.0, 0.0), T)
                                .unwrap(),
                        ),
                        7.5,
                        T,
                    )
                    .unwrap(),
                ),
                ..FitOptions::default()
            },
        ),
    ] {
        let fit = fit::cylinder(&points, &options, T).unwrap();
        assert!(
            (fit.shape.radius() - 8.0).abs() < 1e-6,
            "cylinder from {what}: {:?}",
            fit.shape
        );
        assert!(fit.sigma < 1e-6, "cylinder from {what}: {}", fit.sigma);
    }

    let (reference, angle) = (6.0, 0.35_f64);
    let mut draw = Seeded(161);
    let (mut points, mut normals) = (Vec::new(), Vec::new());
    for _ in 0..2000 {
        let (u, h) = (draw.range(0.0, 2.5), draw.range(-5.0, 5.0));
        let radial = frame.x().vector() * u.cos() + frame.y().vector() * u.sin();
        points.push(
            frame.origin() + radial * angle.tan().mul_add(h, reference) + frame.z().vector() * h,
        );
        normals.push(radial * angle.cos() - frame.z().vector() * angle.sin());
    }
    for (what, options) in [
        (
            "normals",
            FitOptions {
                normals: Some(&normals),
                ..FitOptions::default()
            },
        ),
        ("direction", holding(held_direction)),
        ("axis", holding(held_axis)),
    ] {
        let fit = fit::cone(&points, &options, T).unwrap();
        assert!(
            (fit.shape.half_angle() - angle).abs() < 1e-7,
            "cone from {what}: {:?}",
            fit.shape
        );
        assert!(fit.sigma < 1e-6, "cone from {what}: {}", fit.sigma);
    }

    // A small patch of a thick torus, its tube many times the patch.
    let (major, minor) = (40.0, 10.0);
    let (mut points, mut normals) = (Vec::new(), Vec::new());
    for _ in 0..2000 {
        let (u, v) = (draw.range(0.3, 0.4), draw.range(1.0, 1.3));
        let radial = frame.x().vector() * u.cos() + frame.y().vector() * u.sin();
        let n = radial * v.cos() + frame.z().vector() * v.sin();
        points.push(frame.origin() + radial * major + n * minor);
        normals.push(n);
    }
    for (what, options) in [
        (
            "normals",
            FitOptions {
                normals: Some(&normals),
                ..FitOptions::default()
            },
        ),
        (
            "direction",
            FitOptions {
                normals: Some(&normals),
                fixed: held_direction,
                ..FitOptions::default()
            },
        ),
        ("axis", holding(held_axis)),
    ] {
        let fit = fit::torus(&points, &options, T).unwrap();
        let t = fit.shape;
        assert!(
            (t.minor_radius() - minor).abs() < 1e-6 && (t.major_radius() - major).abs() < 1e-6,
            "torus from {what}: {t:?}"
        );
    }

    let mut ring = Vec::new();
    for _ in 0..500 {
        let u = draw.range(0.0, 3.0);
        ring.push(
            frame.origin() + (frame.x().vector() * u.cos() + frame.y().vector() * u.sin()) * 15.0,
        );
    }
    for (what, options) in [
        ("free", FitOptions::default()),
        ("direction", holding(held_direction)),
        ("axis", holding(held_axis)),
    ] {
        let fit = fit::circle3(&ring, &options, T).unwrap();
        assert!(
            (fit.shape.radius() - 15.0).abs() < 1e-7
                && fit.shape.centre().distance(frame.origin()) < 1e-6,
            "circle from {what}: {:?}",
            fit.shape
        );
    }
}
