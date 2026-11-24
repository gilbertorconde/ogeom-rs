//! Scattered points fitted over a seed surface: noisy samples of a cylinder
//! band, a spherical cap and shapes that are no height field over any
//! plane, each fitted over a seed and measured against the true shape,
//! which is known in closed form.
#![allow(clippy::unwrap_used, reason = "test code")]
#![allow(clippy::cast_precision_loss, reason = "small counts")]

use core::f64::consts::{PI, TAU};

use ogeom_core::Tolerances;
use ogeom_geom::fit::{Ended, ScatteredFit, ScatteredOptions, Target, fit_surface_scattered_over};
use ogeom_geom::{
    CylinderSurface, PlaneSurface, SphereSurface, Surface as _, SurfaceGeometry, project_on_surface,
};
use ogeom_math::{Cylinder, Direction, Frame, Plane, Point, Sphere, Vector};

const T: Tolerances = Tolerances::millimetres();

/// The noise on every cloud, standard deviation in millimetres.
const NOISE: f64 = 0.02;

/// The tolerance a point is held to: five deviations of the noise, past
/// what the noise reaches over a few thousand points (under four and a
/// half), so a point beyond it is a miss of the shape.
const TOLERANCE: f64 = 5.0 * NOISE;

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

fn options(spans: (usize, usize), most: (usize, usize), target: Target) -> ScatteredOptions {
    ScatteredOptions {
        degree: 3,
        spans,
        most,
        fairing: 1e-12,
        target,
    }
}

/// The widest distance from the fitted surface, read on a 61 by 61 grid
/// of its parameters, to the true shape as `off` measures it.
fn widest_off(fit: &ScatteredFit, off: impl Fn(Point) -> f64) -> f64 {
    let mut widest = 0.0_f64;
    for i in 0..=60 {
        for j in 0..=60 {
            let p = fit
                .surface
                .point_at(f64::from(i) / 60.0, f64::from(j) / 60.0, T)
                .unwrap();
            widest = widest.max(off(p));
        }
    }
    widest
}

/// `count` points on a radius-`r(angle, z)` surface round the `z` axis,
/// angles in `angles`, heights in `heights`, moved off it radially by the
/// noise.
fn radial_cloud(
    count: usize,
    angles: (f64, f64),
    heights: (f64, f64),
    r: impl Fn(f64, f64) -> f64,
    draws: &mut Draws,
) -> Vec<Point> {
    (0..count)
        .map(|_| {
            let a = angles.0 + (angles.1 - angles.0) * draws.unit();
            let z = heights.0 + (heights.1 - heights.0) * draws.unit();
            let rho = r(a, z) + NOISE * draws.gaussian();
            Point::new(rho * a.cos(), rho * a.sin(), z)
        })
        .collect()
}

/// A seed cylinder off the true one: a millimetre smaller, its axis half
/// a millimetre aside and tilted a little.
fn cylinder_seed() -> SurfaceGeometry {
    let frame = Frame::new(
        Point::new(0.5, -0.3, 0.0),
        Direction::new(Vector::new(0.01, 0.0, 1.0), T).unwrap(),
        Direction::X,
        T,
    )
    .unwrap();
    SurfaceGeometry::Cylinder(
        CylinderSurface::new(Cylinder::new(frame, 19.0, T).unwrap(), (0.0, 1.0)).unwrap(),
    )
}

/// A 120 degree band of a radius-20 cylinder, straddling the seed's seam at
/// angle nought, recovered within the tolerance everywhere over the fit,
/// its residual the noise.
#[test]
fn a_cylinder_band_across_the_seam_fits_over_a_cylinder_seed() {
    let mut draws = Draws(0x5eed_0001);
    let band = radial_cloud(
        3000,
        (-PI / 3.0, PI / 3.0),
        (0.0, 30.0),
        |_, _| 20.0,
        &mut draws,
    );
    let tolerance = TOLERANCE;
    let fit = fit_surface_scattered_over(
        &band,
        &cylinder_seed(),
        &options((2, 2), (16, 16), Target::Tolerance(tolerance)),
        T,
    )
    .unwrap();
    assert!(fit.met(), "{:?} at {:?} spans", fit.ended, fit.spans);
    let off = widest_off(&fit, |p| (p.x.hypot(p.y) - 20.0).abs());
    assert!(off < tolerance, "the fit strays {off} from the cylinder");
    assert!(
        (fit.sigma / NOISE - 1.0).abs() < 0.1,
        "the residual {} is not the noise",
        fit.sigma
    );
}

/// A cap of a radius-25 sphere, 40 degrees about a point of the seed's
/// equator on its seam, fitted over a smaller sphere seed off centre.
#[test]
fn a_spherical_cap_fits_over_a_sphere_seed() {
    let mut draws = Draws(0x5eed_0002);
    let reach = 40.0_f64.to_radians();
    let cap: Vec<Point> = (0..3000)
        .map(|_| {
            // Uniform over the cap's area: the cosine of the angle from
            // its centre uniform between the rim's and one.
            let c = 1.0 - (1.0 - reach.cos()) * draws.unit();
            let s = (1.0 - c * c).sqrt();
            let a = TAU * draws.unit();
            let along = Vector::new(c, s * a.cos(), s * a.sin());
            let rho = 25.0 + NOISE * draws.gaussian();
            Point::from_vector(along * rho)
        })
        .collect();
    let seed = SurfaceGeometry::Sphere(SphereSurface::new(
        Sphere::new(
            Frame::new(Point::new(1.0, 0.5, -0.5), Direction::Z, Direction::X, T).unwrap(),
            23.0,
            T,
        )
        .unwrap(),
    ));
    let tolerance = TOLERANCE;
    let fit = fit_surface_scattered_over(
        &cap,
        &seed,
        &options((2, 2), (16, 16), Target::Tolerance(tolerance)),
        T,
    )
    .unwrap();
    assert!(fit.met(), "{:?} at {:?} spans", fit.ended, fit.spans);
    // The fit runs over the square the cap's parameters span, past the
    // cap's round rim at the corners, where no point holds it: it is read
    // where the points are, at each point's foot.
    let mut off = 0.0_f64;
    for &(u, v) in &fit.parameters {
        let p = fit.surface.point_at(u, v, T).unwrap();
        off = off.max((p.to_vector().magnitude() - 25.0).abs());
    }
    assert!(off < tolerance, "the fit strays {off} from the sphere");
    assert!((fit.sigma / NOISE - 1.0).abs() < 0.1, "{}", fit.sigma);
}

/// A wavy band three quarters round a cylinder, folding over every plane:
/// the principal-plane fit cannot unfold it, the seeded one does, and
/// refines along the waves only.
#[test]
fn a_wavy_band_no_plane_unfolds_fits_over_a_cylinder_seed() {
    let mut draws = Draws(0x5eed_0003);
    let r = |a: f64, _: f64| 0.8_f64.mul_add((4.0 * a).sin(), 20.0);
    let band = radial_cloud(4000, (0.0, 1.5 * PI), (0.0, 20.0), r, &mut draws);
    let tolerance = TOLERANCE;
    let fit = fit_surface_scattered_over(
        &band,
        &cylinder_seed(),
        &options((2, 2), (64, 64), Target::Tolerance(tolerance)),
        T,
    )
    .unwrap();
    assert!(fit.met(), "{:?} at {:?} spans", fit.ended, fit.spans);
    // The radial gap bounds the distance to the true surface from above.
    let off = widest_off(&fit, |p| {
        let a = p.y.atan2(p.x).rem_euclid(TAU);
        (p.x.hypot(p.y) - r(a, p.z)).abs()
    });
    assert!(off < tolerance, "the fit strays {off} from the wavy band");
    // The height runs straight, so its spans stay as they started.
    assert_eq!(fit.spans.1, 2, "{:?}", fit.spans);
    assert!(fit.spans.0 > 2, "{:?}", fit.spans);

    let flat = ogeom_geom::fit::fit_surface_scattered(&band, 3, (12, 12), 1e-3, T).unwrap();
    assert!(
        flat.error > 1.0,
        "the principal plane fit was expected to fold, missing by {}",
        flat.error
    );
}

/// A plane with one bump near a corner: refinement splits only the spans
/// over the bump. A point planted off the surface is reported as the
/// worst, at the distance the projector finds to the whole surface, and
/// the refinement does not bend the fit to it.
#[test]
fn refinement_splits_only_where_the_fit_misses_and_names_the_worst_point() {
    let mut draws = Draws(0x5eed_0004);
    let bump = |x: f64, y: f64| 1.5 * (-((x - 70.0).powi(2) + (y - 75.0).powi(2)) / 50.0).exp();
    let mut cloud: Vec<Point> = (0..6000)
        .map(|_| {
            let (x, y) = (100.0 * draws.unit(), 100.0 * draws.unit());
            Point::new(x, y, bump(x, y) + NOISE * draws.gaussian())
        })
        .collect();
    let seed = SurfaceGeometry::Plane(PlaneSurface::new(Plane::XY));
    let opts = options((4, 4), (64, 64), Target::Tolerance(TOLERANCE));
    let fit = fit_surface_scattered_over(&cloud, &seed, &opts, T).unwrap();
    assert!(fit.met(), "{:?} at {:?} spans", fit.ended, fit.spans);
    assert!(fit.worst.1 <= TOLERANCE && fit.sigma <= TOLERANCE);
    let off = widest_off(&fit, |p| (p.z - bump(p.x, p.y)).abs());
    assert!(off < TOLERANCE, "the fit strays {off} from the bump");

    // Every knot added lies over the bump: x and y from 55 to 90 of the
    // 0 to 100 the points span, give or take the half span a split puts
    // its knot inside.
    let added = |knots: &ogeom_math::KnotVector| -> Vec<f64> {
        knots
            .distinct()
            .iter()
            .map(|k| k.0)
            .filter(|k| ![0.0, 0.25, 0.5, 0.75, 1.0].contains(k))
            .collect()
    };
    let (au, av) = (added(fit.surface.u_knots()), added(fit.surface.v_knots()));
    assert!(!au.is_empty() && !av.is_empty(), "{au:?} {av:?}");
    for k in au.iter().chain(&av) {
        assert!((0.5..=0.9).contains(k), "a knot at {k} is off the bump");
    }
    // Knots as close as the closest added, all over, would take 32 spans
    // a side; the bump takes a handful.
    let finest = |knots: &ogeom_math::KnotVector| {
        let d = knots.distinct();
        d.windows(2).map(|w| w[1].0 - w[0].0).fold(1.0, f64::min)
    };
    assert!(finest(fit.surface.u_knots()) <= 1.0 / 32.0);
    assert!(au.len() + av.len() <= 12, "{au:?} {av:?}");
    assert_eq!((fit.spans.0 - 4, fit.spans.1 - 4), (au.len(), av.len()));

    let planted = 1234;
    cloud[planted].z += 0.4;
    let marred = fit_surface_scattered_over(&cloud, &seed, &opts, T).unwrap();
    assert_eq!(marred.ended, Ended::Strays);
    assert_eq!(marred.worst.0, planted);
    let surface = SurfaceGeometry::BSpline(marred.surface.clone());
    let mut widest = 0.0_f64;
    for p in &cloud {
        widest = widest.max(project_on_surface(&surface, *p, 16, T).unwrap().distance);
    }
    assert!(
        (marred.worst.1 - widest).abs() < 1e-9,
        "{} {widest}",
        marred.worst.1
    );
    assert!(
        (marred.worst.1 - 0.4).abs() < 4.0 * NOISE,
        "{}",
        marred.worst.1
    );
    assert_eq!(marred.spans, fit.spans);
}

/// A saddle over a plane seed fits within the noise.
#[test]
fn a_saddle_fits_over_a_plane_seed() {
    let mut draws = Draws(0x5eed_0005);
    let z = |x: f64, y: f64| (x * x - y * y) / 40.0;
    let cloud: Vec<Point> = (0..3000)
        .map(|_| {
            let (x, y) = (40.0 * draws.unit() - 20.0, 40.0 * draws.unit() - 20.0);
            Point::new(x, y, z(x, y) + NOISE * draws.gaussian())
        })
        .collect();
    let seed = SurfaceGeometry::Plane(PlaneSurface::new(Plane::XY));
    let fit = fit_surface_scattered_over(
        &cloud,
        &seed,
        &options((2, 2), (16, 16), Target::Tolerance(TOLERANCE)),
        T,
    )
    .unwrap();
    assert!(fit.met(), "{:?}", fit.ended);
    let off = widest_off(&fit, |p| (p.z - z(p.x, p.y)).abs());
    assert!(off < TOLERANCE, "{off}");
}

/// Fits seeded each with the one before converge. On a cylinder seed off
/// the true axis the parameters run unevenly over the shape; on a fit they
/// are its nearest points, so over the same spans each fit comes closer
/// than the one before, by less each time. With a tolerance, the fit seeded
/// with the first meets it too.
#[test]
fn a_previous_fit_seeds_the_next() {
    let mut draws = Draws(0x5eed_0006);
    let r = |a: f64, z: f64| 0.5_f64.mul_add((3.0 * a).sin() * (z / 5.0).cos(), 20.0);
    let band = radial_cloud(4000, (0.0, PI), (0.0, 20.0), r, &mut draws);
    let reseeded = |fit: &ScatteredFit, opts: &ScatteredOptions| {
        let seed = SurfaceGeometry::BSpline(fit.surface.clone());
        fit_surface_scattered_over(&band, &seed, opts, T).unwrap()
    };

    let fixed = options((4, 3), (4, 3), Target::Fixed);
    let mut fit = fit_surface_scattered_over(&band, &cylinder_seed(), &fixed, T).unwrap();
    let mut step = f64::INFINITY;
    for _ in 0..4 {
        let next = reseeded(&fit, &fixed);
        let fell = fit.sigma - next.sigma;
        assert!(fell > 0.0 && fell < step, "{} -> {}", fit.sigma, next.sigma);
        (fit, step) = (next, fell);
    }

    let opts = options((2, 2), (64, 64), Target::Tolerance(TOLERANCE));
    let first = fit_surface_scattered_over(&band, &cylinder_seed(), &opts, T).unwrap();
    assert!(first.met(), "{:?}", first.ended);
    let second = reseeded(&first, &opts);
    assert!(second.met(), "{:?}", second.ended);
    let off = widest_off(&second, |p| {
        let a = p.y.atan2(p.x).rem_euclid(TAU);
        (p.x.hypot(p.y) - r(a, p.z)).abs()
    });
    assert!(off < TOLERANCE, "the second fit strays {off} from the band");
}

/// The residual left is the noise. A tolerance inside the noise's reach
/// refines a wavy band until a split would fit noise rather than shape, and
/// stops there: the residual is the noise's deviation, and points the fit
/// did not see stand as close. A net fine enough to follow the noise shows
/// it: its residual falls under the noise while the points it did not see
/// stand farther off.
#[test]
fn noise_is_fitted_to_its_level_and_overfitting_shows() {
    let mut draws = Draws(0x5eed_0007);
    let r = |a: f64, _: f64| 0.8_f64.mul_add((4.0 * a).sin(), 20.0);
    let band = radial_cloud(3000, (0.0, PI), (0.0, 20.0), r, &mut draws);
    let seed = cylinder_seed();
    let under = fit_surface_scattered_over(
        &band,
        &seed,
        &options((2, 2), (32, 32), Target::Tolerance(2.0 * NOISE)),
        T,
    )
    .unwrap();
    assert_eq!(under.ended, Ended::Noise, "{:?}", under.spans);
    assert!((under.sigma / NOISE - 1.0).abs() < 0.05, "{}", under.sigma);
    assert!(
        under.cross_sigma < 1.05 * under.sigma,
        "{} {}",
        under.cross_sigma,
        under.sigma
    );
    assert!(
        under.spans.0 < 32 && under.spans.1 < 32,
        "{:?}",
        under.spans
    );

    let over =
        fit_surface_scattered_over(&band, &seed, &options((24, 24), (24, 24), Target::Fixed), T)
            .unwrap();
    assert_eq!(over.ended, Ended::Fixed);
    assert!(over.sigma < 0.9 * NOISE, "{}", over.sigma);
    assert!(
        over.cross_sigma > 1.5 * over.sigma,
        "{} {}",
        over.cross_sigma,
        over.sigma
    );
}

/// The same points and options give the same fit, to the bit.
#[test]
fn a_seeded_fit_is_deterministic() {
    let mut draws = Draws(0x5eed_0008);
    let r = |a: f64, _: f64| 0.8_f64.mul_add((4.0 * a).sin(), 20.0);
    let band = radial_cloud(2000, (0.0, PI), (0.0, 20.0), r, &mut draws);
    let fit = || {
        fit_surface_scattered_over(
            &band,
            &cylinder_seed(),
            &options((2, 2), (32, 32), Target::Tolerance(TOLERANCE)),
            T,
        )
        .unwrap()
    };
    assert_eq!(fit(), fit());
}

/// What the fit cannot answer is refused by name: a tolerance that is no
/// distance, a span budget below its start, and a seed that does not unfold
/// the points (a line along a cylinder seed's axis has one angle).
#[test]
fn a_seeded_fit_refuses_what_it_cannot_answer() {
    let mut draws = Draws(0x5eed_0009);
    let band = radial_cloud(200, (0.0, PI), (0.0, 20.0), |_, _| 20.0, &mut draws);
    let seed = cylinder_seed();
    for opts in [
        options((2, 2), (4, 4), Target::Tolerance(0.0)),
        options((2, 2), (4, 4), Target::Tolerance(f64::NAN)),
        options((4, 4), (2, 2), Target::Fixed),
        options((0, 2), (4, 4), Target::Fixed),
    ] {
        assert!(
            fit_surface_scattered_over(&band, &seed, &opts, T).is_err(),
            "{opts:?}"
        );
    }
    let line: Vec<Point> = (0..50)
        .map(|i| Point::new(20.0, 0.0, f64::from(i)))
        .collect();
    let upright = SurfaceGeometry::Cylinder(
        CylinderSurface::new(Cylinder::new(Frame::WORLD, 19.0, T).unwrap(), (0.0, 1.0)).unwrap(),
    );
    let fixed = options((2, 2), (4, 4), Target::Fixed);
    let error = fit_surface_scattered_over(&line, &upright, &fixed, T).unwrap_err();
    assert!(error.to_string().contains("do not spread"), "{error}");
}
