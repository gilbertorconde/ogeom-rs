//! Measuring an intersection against ground truth.
//!
//! `docs/SCOPE.md` §7 gates the boolean pipeline on this: surface/surface
//! quality is *measured* before anything is built on top of it, and if it does
//! not clear the bar the project ships as a geometry library rather than
//! spending years on a boolean over an intersector that cannot carry one.
//!
//! This is the instrument. It is not a test — it is the thing a test asserts
//! about, and the thing a report quotes.
//!
//! # What ground truth is here
//!
//! Not another kernel's answer. The two surfaces themselves.
//!
//! An intersection curve has exactly one defining property: every point of it
//! lies on **both** surfaces. That is checkable without reference to anything
//! else, it is checkable to machine precision, and it does not care how the
//! curve was arrived at. So the measure is the largest distance from any sampled
//! point of the reported curve to either surface — which is zero for a correct
//! result and says how wrong an incorrect one is.
//!
//! # What it deliberately does not measure
//!
//! *Completeness.* Every point being on both surfaces says the answer is not
//! wrong; it does not say the answer is all of it. An intersector that returned
//! one of two circles would score perfectly here. Completeness needs a second
//! instrument — sampling one surface and asking whether points near the other
//! are covered — and the general intersector is where that starts to matter,
//! since these closed forms return the whole answer by construction.
//!
//! Saying so is the point. A benchmark that quietly measured one thing while
//! being read as measuring another would be worse than none, because the number
//! would be trusted.

use og_core::{OgResult, Tolerances};
use og_geom::{Curve, Curve3d, Surface, SurfaceGeometry};
use og_math::Point;

use crate::surface::{Meeting, surface_surface};

/// How well one intersection was solved.
#[derive(Debug, Clone, PartialEq)]
pub struct Measured {
    /// What the intersector said.
    pub meeting: Meeting,
    /// The largest distance from any sampled point to the surfaces it should
    /// lie on.
    ///
    /// `None` when there was nothing to sample — the surfaces are apart, or the
    /// same — because a deviation of zero and no measurement at all are
    /// different things and averaging them together would flatter the result.
    pub deviation: Option<f64>,
    /// How many points were sampled.
    pub samples: usize,
}

/// A summary over many cases.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Report {
    /// Cases attempted.
    pub cases: usize,
    /// Cases the closed forms answered.
    pub solved: usize,
    /// Cases reported as needing the general intersector.
    pub deferred: usize,
    /// The worst deviation seen across every case that produced a curve.
    pub worst: f64,
    /// Which case that was.
    pub worst_case: Option<String>,
}

impl Report {
    /// Whether every solved case met a bar.
    ///
    /// Deliberately not a stored field: the bar belongs to whoever is asking,
    /// and burying one in the report would make every reader inherit a
    /// threshold they did not choose.
    #[must_use]
    pub fn within(&self, bar: f64) -> bool {
        self.worst <= bar
    }
}

/// Measure one pair.
///
/// # Errors
///
/// As [`surface_surface`] — a pair with no closed form is reported rather than
/// measured, and the caller decides whether that counts against it.
pub fn measure(a: &SurfaceGeometry, b: &SurfaceGeometry, tol: Tolerances) -> OgResult<Measured> {
    /// Enough to catch a curve that is right at its ends and wrong in the
    /// middle, which is the shape of a wrong answer that endpoints miss.
    const SAMPLES: usize = 64;

    let meeting = surface_surface(a, b, tol)?;
    let mut worst: Option<f64> = None;
    let mut samples = 0;

    let mut check = |p: Point| {
        let off = distance_to(a, p, tol).max(distance_to(b, p, tol));
        worst = Some(worst.map_or(off, |w: f64| w.max(off)));
        samples += 1;
    };

    match &meeting {
        Meeting::Along(curves) => {
            for curve in curves {
                let (lo, hi) = sampling_range(curve);
                for i in 0..=SAMPLES {
                    #[allow(clippy::cast_precision_loss)]
                    let u = lo + (hi - lo) * i as f64 / SAMPLES as f64;
                    if let Ok(p) = curve.point_at(u, tol) {
                        check(p);
                    }
                }
            }
        }
        Meeting::Touching(points) => {
            for p in points {
                check(*p);
            }
        }
        // Nothing to sample: see `Measured::deviation`.
        Meeting::Apart | Meeting::Same => {}
    }

    Ok(Measured {
        meeting,
        deviation: worst,
        samples,
    })
}

/// Measure a whole set of named cases.
///
/// # Errors
///
/// Never for a case with no closed form — those are counted as deferred, which
/// is the honest reading: the closed forms are not claiming to answer them.
pub fn measure_all(
    cases: &[(String, SurfaceGeometry, SurfaceGeometry)],
    tol: Tolerances,
) -> Report {
    let mut report = Report {
        cases: cases.len(),
        ..Report::default()
    };
    for (name, a, b) in cases {
        match measure(a, b, tol) {
            Ok(found) => {
                report.solved += 1;
                if let Some(off) = found.deviation
                    && off > report.worst
                {
                    report.worst = off;
                    report.worst_case = Some(name.clone());
                }
            }
            Err(_) => report.deferred += 1,
        }
    }
    report
}

/// The range to sample a curve over.
///
/// An unbounded line's own domain reaches a billion units either way, and
/// sampling that says nothing useful about an intersection near the origin — the
/// interesting part is where the surfaces actually are. A bounded curve is
/// sampled over all of itself.
fn sampling_range(curve: &Curve) -> (f64, f64) {
    let (lo, hi) = curve.domain();
    if matches!(curve, Curve::Line(_)) {
        return (-100.0, 100.0);
    }
    (lo, hi)
}

/// How far a point is from a surface.
///
/// Uses each quadric's own closed-form distance where it has one, which is what
/// makes this an independent check rather than a restatement of the
/// intersector's arithmetic.
fn distance_to(surface: &SurfaceGeometry, p: Point, tol: Tolerances) -> f64 {
    use SurfaceGeometry as S;
    match surface {
        S::Plane(x) => x.plane().distance_to(p),
        S::Sphere(x) => x.sphere().distance_to(p),
        S::Cylinder(x) => x.cylinder().distance_to(p),
        S::Cone(x) => x.cone().distance_to(p),
        S::Torus(x) => x.torus().distance_to(p),
        // No closed form: fall back to the nearest sampled point, which is a
        // weaker check and is only reached for surfaces the closed forms do not
        // handle anyway.
        other => nearest_on(other, p, tol),
    }
}

/// The distance to the nearest sampled point of a surface.
fn nearest_on(surface: &SurfaceGeometry, p: Point, tol: Tolerances) -> f64 {
    const STEPS: usize = 128;
    let ((ua, ub), (va, vb)) = surface.domain();
    let mut best = f64::MAX;
    for i in 0..=STEPS {
        for j in 0..=STEPS {
            #[allow(clippy::cast_precision_loss)]
            let (s, t) = (i as f64 / STEPS as f64, j as f64 / STEPS as f64);
            if let Ok(q) = surface.point_at(ua + (ub - ua) * s, va + (vb - va) * t, tol) {
                best = best.min(p.distance(q));
            }
        }
    }
    best
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::print_stdout)]
mod tests {
    use super::*;
    use og_geom::{CylinderSurface, PlaneSurface, SphereSurface};
    use og_math::{Cylinder, Direction, Frame, Plane, Point, Sphere, Vector};

    const T: Tolerances = Tolerances::millimetres();

    fn plane(origin: Point, normal: Vector) -> SurfaceGeometry {
        PlaneSurface::new(Plane::through(origin, Direction::new(normal, T).unwrap())).into()
    }

    fn sphere(centre: Point, radius: f64) -> SurfaceGeometry {
        SphereSurface::new(Sphere::centred(centre, radius, T).unwrap()).into()
    }

    fn cylinder(origin: Point, axis: Vector, radius: f64) -> SurfaceGeometry {
        let frame = Frame::new(
            origin,
            Direction::new(axis, T).unwrap(),
            Direction::from_cross(axis, Vector::new(0.3, 0.5, 0.9), T).unwrap(),
            T,
        )
        .unwrap();
        CylinderSurface::new(Cylinder::new(frame, radius, T).unwrap(), (-10.0, 10.0))
            .unwrap()
            .into()
    }

    /// Every case with a closed form, named.
    fn corpus() -> Vec<(String, SurfaceGeometry, SurfaceGeometry)> {
        let mut out = Vec::new();
        let mut add = |name: &str, a: SurfaceGeometry, b: SurfaceGeometry| {
            out.push((name.to_string(), a, b));
        };

        add(
            "plane/plane crossing",
            plane(Point::ORIGIN, Vector::Z),
            plane(Point::ORIGIN, Vector::X),
        );
        add(
            "plane/plane oblique",
            plane(Point::new(1.0, 2.0, 3.0), Vector::new(1.0, 1.0, 1.0)),
            plane(Point::new(-2.0, 0.5, 1.0), Vector::new(0.2, -1.0, 0.7)),
        );
        add(
            "plane/sphere through the centre",
            plane(Point::ORIGIN, Vector::Z),
            sphere(Point::ORIGIN, 3.0),
        );
        add(
            "plane/sphere off centre",
            plane(Point::new(0.0, 0.0, 1.5), Vector::Z),
            sphere(Point::ORIGIN, 3.0),
        );
        add(
            "plane/sphere oblique",
            plane(Point::new(0.4, -0.2, 0.9), Vector::new(1.0, 2.0, 3.0)),
            sphere(Point::new(1.0, 1.0, 1.0), 4.0),
        );
        add(
            "plane/cylinder perpendicular",
            plane(Point::new(0.0, 0.0, 2.0), Vector::Z),
            cylinder(Point::ORIGIN, Vector::Z, 2.0),
        );
        add(
            "plane/cylinder oblique",
            plane(Point::ORIGIN, Vector::new(0.0, 1.0, 1.0)),
            cylinder(Point::ORIGIN, Vector::Z, 2.0),
        );
        add(
            "plane/cylinder along the axis",
            plane(Point::new(0.5, 0.0, 0.0), Vector::X),
            cylinder(Point::ORIGIN, Vector::Z, 2.0),
        );
        add(
            "plane/cylinder tangent",
            plane(Point::new(2.0, 0.0, 0.0), Vector::X),
            cylinder(Point::ORIGIN, Vector::Z, 2.0),
        );
        add(
            "sphere/sphere crossing",
            sphere(Point::ORIGIN, 3.0),
            sphere(Point::new(4.0, 0.0, 0.0), 2.0),
        );
        add(
            "sphere/sphere oblique",
            sphere(Point::new(1.0, -2.0, 0.5), 5.0),
            sphere(Point::new(-3.0, 1.0, 2.0), 4.0),
        );
        add(
            "cylinder/sphere coaxial",
            cylinder(Point::ORIGIN, Vector::Z, 1.5),
            sphere(Point::ORIGIN, 3.0),
        );
        out
    }

    #[test]
    fn every_closed_form_lands_on_both_surfaces() {
        // The gate's own measurement, run as an assertion. Every point of every
        // reported curve is on both surfaces to machine precision — which is
        // the defining property of an intersection curve and the only one that
        // can be checked without a second implementation to compare against.
        let report = measure_all(&corpus(), T);
        println!(
            "intersection benchmark: {} cases, {} solved, {} deferred, worst \
             deviation {:e}{}",
            report.cases,
            report.solved,
            report.deferred,
            report.worst,
            report
                .worst_case
                .as_ref()
                .map_or(String::new(), |c| format!(" ({c})"))
        );

        assert_eq!(
            report.deferred, 0,
            "every case here should have a closed form"
        );
        assert_eq!(report.solved, report.cases);
        assert!(
            report.within(1e-12),
            "worst deviation {:e} at {:?}",
            report.worst,
            report.worst_case
        );
    }

    #[test]
    fn the_kind_of_answer_is_right_and_not_only_its_accuracy() {
        // Landing on both surfaces is necessary and not sufficient: an
        // intersector returning one circle of two would score perfectly. These
        // pin the *shape* of each answer.
        let case = |a: SurfaceGeometry, b: SurfaceGeometry| surface_surface(&a, &b, T).unwrap();

        assert!(matches!(
            case(plane(Point::ORIGIN, Vector::Z), plane(Point::ORIGIN, Vector::X)),
            Meeting::Along(ref c) if c.len() == 1
        ));
        assert_eq!(
            case(
                plane(Point::ORIGIN, Vector::Z),
                plane(Point::new(0.0, 0.0, 1.0), Vector::Z)
            ),
            Meeting::Apart
        );
        assert_eq!(
            case(
                plane(Point::ORIGIN, Vector::Z),
                plane(Point::new(3.0, 4.0, 0.0), Vector::Z)
            ),
            Meeting::Same
        );

        // A sphere resting on a plane touches; it does not meet along anything.
        assert!(matches!(
            case(
                plane(Point::ORIGIN, Vector::Z),
                sphere(Point::new(0.0, 0.0, 2.0), 2.0)
            ),
            Meeting::Touching(ref p) if p.len() == 1
        ));

        // A plane through a cylinder's axis cuts two lines, not one.
        assert!(matches!(
            case(
                plane(Point::ORIGIN, Vector::X),
                cylinder(Point::ORIGIN, Vector::Z, 2.0)
            ),
            Meeting::Along(ref c) if c.len() == 2
        ));

        // A sphere larger than a coaxial cylinder cuts it in two circles.
        assert!(matches!(
            case(
                cylinder(Point::ORIGIN, Vector::Z, 1.0),
                sphere(Point::ORIGIN, 3.0)
            ),
            Meeting::Along(ref c) if c.len() == 2
        ));
    }

    #[test]
    fn an_oblique_plane_cuts_a_cylinder_in_an_ellipse_of_the_right_size() {
        // The closed form's whole value: not a fitted curve that is nearly an
        // ellipse, but the ellipse, with the radii geometry says it has.
        let angle = core::f64::consts::FRAC_PI_3;
        let radius = 2.0;
        let cut = plane(Point::ORIGIN, Vector::new(0.0, angle.sin(), angle.cos()));
        let drum = cylinder(Point::ORIGIN, Vector::Z, radius);
        let Meeting::Along(curves) = surface_surface(&cut, &drum, T).unwrap() else {
            panic!("an oblique cut should meet along a curve");
        };
        assert_eq!(curves.len(), 1);
        let Curve::Ellipse(e) = &curves[0] else {
            panic!(
                "an oblique cut of a cylinder is an ellipse, got {:?}",
                curves[0]
            );
        };
        approx::assert_relative_eq!(e.ellipse().minor_radius(), radius, max_relative = 1e-12);
        approx::assert_relative_eq!(
            e.ellipse().major_radius(),
            radius / angle.cos(),
            max_relative = 1e-12
        );
    }

    #[test]
    fn a_pair_with_no_closed_form_is_deferred_rather_than_guessed() {
        // The honest half. Two cylinders on skew axes meet in a quartic space
        // curve; returning something plausible would be the single worst thing
        // this module could do, because the boolean above it would trust it.
        let a = cylinder(Point::ORIGIN, Vector::Z, 1.0);
        let b = cylinder(Point::ORIGIN, Vector::X, 1.0);
        let err = surface_surface(&a, &b, T).unwrap_err();
        assert!(
            err.to_string().contains("marching"),
            "unexpected message: {err}"
        );

        let deferred = measure_all(&[("cylinder/cylinder skew".to_string(), a, b)], T);
        assert_eq!(deferred.deferred, 1);
        assert_eq!(deferred.solved, 0);
        assert_eq!(deferred.worst, 0.0, "a deferred case scores nothing");
    }

    #[test]
    fn nothing_to_sample_is_reported_as_nothing_rather_than_as_zero_error() {
        // A deviation of zero and no measurement at all are different, and
        // averaging them together would flatter every report containing a pair
        // that simply misses.
        let found = measure(
            &plane(Point::ORIGIN, Vector::Z),
            &plane(Point::new(0.0, 0.0, 5.0), Vector::Z),
            T,
        )
        .unwrap();
        assert_eq!(found.meeting, Meeting::Apart);
        assert_eq!(found.deviation, None);
        assert_eq!(found.samples, 0);
    }
}
