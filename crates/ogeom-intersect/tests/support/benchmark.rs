//! Measuring an intersection against ground truth.
//!
//! `docs/PLAN.md` gates the boolean pipeline on this: surface/surface
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

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve3d, Surface, SurfaceGeometry};
use ogeom_math::Point;

use ogeom_intersect::{Meeting, surface_surface};

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
pub fn measure(a: &SurfaceGeometry, b: &SurfaceGeometry, tol: Tolerances) -> OgeomResult<Measured> {
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
