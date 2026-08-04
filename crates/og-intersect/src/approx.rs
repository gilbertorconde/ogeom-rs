//! The approximation stage: a traced branch becomes curves.
//!
//! A traced branch is a polyline with a stated chord tolerance — honest, and
//! not what anything downstream wants to hold. An edge wants a curve in space;
//! a face wants that curve in its *own parameter space*, because splitting a
//! face happens there and a curve the face cannot express is a curve it cannot
//! be split along (`docs/DATA_MODEL.md` §6).
//!
//! So one branch becomes three fits sharing one tolerance: the 3D curve, and
//! one pcurve per surface, each fitted from the samples the tracer already
//! recorded. The tracer kept the parameters on both surfaces at every point
//! precisely for this moment — re-deriving them here would be a projection per
//! point, solving again what the marcher already solved.
//!
//! # The tolerance story, stated once
//!
//! The result's tolerance is a *sum of stated parts*, not a hope: the trace
//! sits within its chord tolerance of the true intersection, and the fit sits
//! within its own reported error of the trace. Both numbers are carried, and
//! the total is what an edge built on this curve must widen its tolerance to.
//! Nothing here rounds a miss up to a hit — a fit that could not reach its
//! target says so, and the caller decides whether the looser curve is usable.
//!
//! # Seams
//!
//! A branch crossing a periodic surface's seam has parameter samples that jump
//! by a period — the pcurve polyline tears even though the curve in space is
//! smooth. The samples are unwrapped before fitting: each step is folded to
//! the nearest image, so the pcurve runs continuously past the seam and may
//! legitimately leave `[0, 2π)`. That is what a pcurve on a periodic surface
//! is; folding it back would re-tear it.

use og_core::{OgResult, Tolerances, og_bail};
use og_geom::{BSpline2d, BSplineCurve, Surface, SurfaceGeometry};
use og_math::Point2;

use crate::march::Traced;

/// A branch of an intersection, as curves.
#[derive(Debug, Clone, PartialEq)]
pub struct IntersectionCurve {
    /// The curve in space.
    pub curve: BSplineCurve,
    /// The same curve in the first surface's parameter space.
    pub on_a: BSpline2d,
    /// And in the second's.
    pub on_b: BSpline2d,
    /// How far the *fits* may sit from the traced polyline.
    ///
    /// The worst of the three fits' reported errors. The distance to the true
    /// intersection adds the trace's own chord tolerance on top; both are
    /// stated so an edge built on this knows what to carry.
    pub fit_error: f64,
    /// Whether every fit met the tolerance it was asked for.
    pub met: bool,
    /// Whether the branch is a closed loop.
    pub closed: bool,
}

/// Fit one traced branch to curves, within `tolerance`.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the branch has
/// fewer than two points or the tolerance is not a positive distance.
pub fn approximate_branch(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    branch: &Traced,
    tolerance: f64,
    tol: Tolerances,
) -> OgResult<IntersectionCurve> {
    if branch.points.len() < 2 {
        og_bail!(
            Construction,
            "a branch of {} points is not a curve",
            branch.points.len()
        );
    }

    // One fit in seven dimensions — the curve and both parameter images
    // together. Fitted separately, each fit's parameter correction drifts
    // its parameterization independently and the three results silently stop
    // being same-parameter: the boolean found pcurves claiming 1e-7 that
    // evaluated millimetres from their own curve. Jointly, one
    // parameterization and one knot vector serve all three, and the reported
    // error bounds every coordinate.
    let unwrapped_a = unwrap_periodic(a, &branch.on_a);
    let unwrapped_b = unwrap_periodic(b, &branch.on_b);
    let (space, on_a, on_b) = og_geom::fit::fit_points_joint(
        &branch.points,
        &unwrapped_a,
        &unwrapped_b,
        3,
        tolerance,
        tol,
    )?;

    Ok(IntersectionCurve {
        fit_error: space
            .error
            .max(space_error(a, &(on_a.clone(), space.met, space.error), tol))
            .max(space_error(b, &(on_b.clone(), space.met, space.error), tol)),
        met: space.met,
        curve: space.curve,
        on_a,
        on_b,
        closed: branch.closed(),
    })
}

/// The fitted pcurve's error, converted back into space.
///
/// The pcurve was fitted in parameter units, against a scale estimated from
/// the whole branch — but the surface's stretch varies along the curve, so an
/// error acceptable in parameter units may be worse in millimetres where the
/// surface stretches hardest. This converts the fit's parameter-space error
/// through the local stretch at samples along the pcurve and reports the
/// worst, so the number the caller reads is in the units the caller measures
/// everything else in.
fn space_error(surface: &SurfaceGeometry, fitted: &(BSpline2d, bool, f64), tol: Tolerances) -> f64 {
    use og_geom::Curve2d;
    let (pcurve, _, parameter_error) = fitted;
    // Convert the parameter-space error back through the surface's local
    // stretch at a few places; take the worst.
    let (lo, hi) = pcurve.domain();
    let mut worst = 0.0_f64;
    for i in 0..=16 {
        #[allow(clippy::cast_precision_loss)]
        let u = lo + (hi - lo) * f64::from(i) / 16.0;
        let Ok(at) = pcurve.point_at(u, tol) else {
            continue;
        };
        let Ok((du, dv)) = surface.d1_at(at.x, at.y, tol) else {
            continue;
        };
        let stretch = du.magnitude().max(dv.magnitude());
        worst = worst.max(parameter_error * stretch);
    }
    worst
}

/// Unfold parameter samples across a periodic surface's seam.
///
/// Each step is folded to the nearest image of the next sample, so a branch
/// crossing `u = 0` continues to `-0.1` rather than tearing to `2π - 0.1`. The
/// result may leave the surface's stated domain, which is what a pcurve
/// crossing a seam *is*.
fn unwrap_periodic(surface: &SurfaceGeometry, samples: &[(f64, f64)]) -> Vec<Point2> {
    let ((ua, ub), (va, vb)) = surface.domain();
    let u_period = if surface.is_periodic_u() {
        Some(ub - ua)
    } else {
        None
    };
    let v_period = if surface.is_periodic_v() {
        Some(vb - va)
    } else {
        None
    };
    let fold = |previous: f64, next: f64, period: Option<f64>| match period {
        None => next,
        Some(period) => {
            let mut candidate = next;
            while candidate - previous > period * 0.5 {
                candidate -= period;
            }
            while previous - candidate > period * 0.5 {
                candidate += period;
            }
            candidate
        }
    };

    let mut out = Vec::with_capacity(samples.len());
    let mut at = Point2::new(samples[0].0, samples[0].1);
    out.push(at);
    for sample in &samples[1..] {
        at = Point2::new(
            fold(at.x, sample.0, u_period),
            fold(at.y, sample.1, v_period),
        );
        out.push(at);
    }
    out
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::march::{Marching, branches};
    use og_geom::{Curve2d, Curve3d, CylinderSurface, PlaneSurface, SphereSurface};
    use og_math::{Cylinder, Direction, Frame, Plane, Point, Sphere, Vector};

    const T: Tolerances = Tolerances::millimetres();

    fn sphere(radius: f64) -> SurfaceGeometry {
        SphereSurface::new(Sphere::centred(Point::ORIGIN, radius, T).unwrap()).into()
    }

    fn cylinder(radius: f64) -> SurfaceGeometry {
        CylinderSurface::new(Cylinder::new(Frame::WORLD, radius, T).unwrap(), (-4.0, 4.0))
            .unwrap()
            .into()
    }

    fn plane(origin: Point, normal: Vector) -> SurfaceGeometry {
        PlaneSurface::over(
            Plane::through(origin, Direction::new(normal, T).unwrap()),
            (-6.0, 6.0),
            (-6.0, 6.0),
        )
        .unwrap()
        .into()
    }

    fn options() -> Marching {
        Marching {
            chord: 1e-5,
            ..Marching::default()
        }
    }

    /// The distance of a fitted curve from both surfaces, sampled densely.
    ///
    /// This is the measure the whole stage exists for: the *fit* — not the
    /// polyline it came from — is what downstream code holds, so the fit is
    /// what must lie on both surfaces.
    fn fitted_deviation(a: &SurfaceGeometry, b: &SurfaceGeometry, curve: &BSplineCurve) -> f64 {
        let off = |surface: &SurfaceGeometry, p: Point| match surface {
            SurfaceGeometry::Plane(x) => x.plane().distance_to(p),
            SurfaceGeometry::Sphere(x) => x.sphere().distance_to(p),
            SurfaceGeometry::Cylinder(x) => x.cylinder().distance_to(p),
            _ => 0.0,
        };
        let (lo, hi) = curve.knots().domain();
        let mut worst = 0.0_f64;
        for i in 0..=800 {
            #[allow(clippy::cast_precision_loss)]
            let u = lo + (hi - lo) * f64::from(i) / 800.0;
            if let Ok(p) = curve.point_at(u, T) {
                worst = worst.max(off(a, p).abs().max(off(b, p).abs()));
            }
        }
        worst
    }

    #[test]
    fn a_fitted_branch_lies_on_both_surfaces_to_the_stated_total() {
        // The tolerance story end to end: trace within 1e-5, fit within 1e-4,
        // so the fitted curve is within the sum of the two of the true
        // intersection — measured against the surfaces, not the polyline.
        let a = sphere(3.0);
        let b = cylinder(1.5);
        let found = branches(&a, &b, options(), T).unwrap();
        assert_eq!(found.len(), 2);

        for branch in &found {
            let fitted = approximate_branch(&a, &b, branch, 1e-4, T).unwrap();
            assert!(fitted.met, "fit error {:e}", fitted.fit_error);
            assert!(fitted.closed);
            let off = fitted_deviation(&a, &b, &fitted.curve);
            assert!(
                off <= 1e-4 + 1e-5,
                "the fitted curve is {off:e} off the surfaces"
            );
            // And it is compact: a curve, not a decorated polyline.
            assert!(
                fitted.curve.control_points().len() * 4 < branch.points.len(),
                "{} control points for {} samples",
                fitted.curve.control_points().len(),
                branch.points.len()
            );
        }
    }

    #[test]
    fn the_pcurves_lift_back_onto_the_curve() {
        // A pcurve is only worth having if evaluating it and lifting through
        // its surface lands on the intersection. Checked through both
        // surfaces at matched ends and sampled interiors.
        let a = sphere(3.0);
        let b = cylinder(1.5);
        let found = branches(&a, &b, options(), T).unwrap();
        let branch = &found[0];
        let fitted = approximate_branch(&a, &b, branch, 1e-4, T).unwrap();

        for (surface, pcurve) in [(&a, &fitted.on_a), (&b, &fitted.on_b)] {
            let (lo, hi) = pcurve.domain();
            for i in 0..=200 {
                #[allow(clippy::cast_precision_loss)]
                let u = lo + (hi - lo) * f64::from(i) / 200.0;
                let at = pcurve.point_at(u, T).unwrap();
                let lifted = surface.point_at(at.x, at.y, T).unwrap();
                // The lifted point is on its own surface by construction; what
                // matters is that it is on the *other* one too, i.e. on the
                // intersection.
                let off = match (surface as &SurfaceGeometry, &a, &b) {
                    _ if core::ptr::eq(surface, &a) => match &b {
                        SurfaceGeometry::Cylinder(c) => c.cylinder().distance_to(lifted),
                        _ => 0.0,
                    },
                    _ => match &a {
                        SurfaceGeometry::Sphere(s) => s.sphere().distance_to(lifted),
                        _ => 0.0,
                    },
                };
                assert!(
                    off.abs() < 5e-4,
                    "a lifted pcurve point is {off:e} off the intersection"
                );
            }
        }
    }

    #[test]
    fn a_branch_across_the_seam_gets_a_continuous_pcurve() {
        // A plane through a cylinder's axis at an angle produces an ellipse
        // whose pcurve crosses the cylinder's u = 0 seam. Folded naively the
        // pcurve tears by 2π; unwrapped it runs smoothly and leaves the stated
        // domain, which is what crossing a seam means.
        let a = cylinder(2.0);
        let b = plane(Point::ORIGIN, Vector::new(0.0, 0.4, 1.0));
        let found = branches(&a, &b, options(), T).unwrap();
        assert_eq!(found.len(), 1, "an oblique plane cuts one ellipse");
        let fitted = approximate_branch(&a, &b, &found[0], 1e-4, T).unwrap();

        // Continuity: no two adjacent samples of the fitted pcurve jump by
        // anything near a period.
        let (lo, hi) = fitted.on_a.domain();
        let mut previous = fitted.on_a.point_at(lo, T).unwrap();
        for i in 1..=400 {
            #[allow(clippy::cast_precision_loss)]
            let u = lo + (hi - lo) * f64::from(i) / 400.0;
            let at = fitted.on_a.point_at(u, T).unwrap();
            assert!(
                (at.x - previous.x).abs() < 1.0,
                "the pcurve tears at the seam: {} to {}",
                previous.x,
                at.x
            );
            previous = at;
        }
    }

    #[test]
    fn what_cannot_be_fitted_is_refused() {
        let a = sphere(1.0);
        let b = plane(Point::ORIGIN, Vector::Z);
        let found = branches(&a, &b, options(), T).unwrap();
        assert!(approximate_branch(&a, &b, &found[0], 0.0, T).is_err());
        assert!(approximate_branch(&a, &b, &found[0], -1.0, T).is_err());

        let empty = Traced {
            points: vec![],
            on_a: vec![],
            on_b: vec![],
            stopped: crate::march::Stopped::Stalled,
        };
        assert!(approximate_branch(&a, &b, &empty, 1e-4, T).is_err());
    }
}
