//! Exact B-spline forms for the analytic curves.
//!
//! Every curve here has a *rational* B-spline form that is exact: not a fit, not
//! an approximation to a tolerance. A circle is a piecewise rational quadratic
//! and lands on the circle at every parameter, which is the whole reason
//! rational weights exist and why `docs/PLAN.md` calls them load-bearing rather
//! than an optional extra.
//!
//! # What conversion is for
//!
//! Three things need it. Exchange formats describe free-form geometry and
//! nothing else, so an exact circle has to become a NURBS to be written at all.
//! A general affine transform (a shear, a non-uniform scale) carries a circle
//! to an ellipse and an ellipse to something with no analytic name, but carries
//! a NURBS to a NURBS by moving its control points, exactly. And an algorithm
//! that only knows one representation can be given every shape in it.
//!
//! # The parameter does not survive, and cannot
//!
//! A circle's parameter is its angle. Its rational quadratic form's is not, and
//! no reparameterization of a rational quadratic makes it one; the two are
//! related by an arctangent. So conversion preserves the *curve* and not the
//! parameterization, and every converted curve is handed back on `[0, 1]`.
//!
//! That is why this is a geometry operation rather than a topology one. An edge
//! converted this way needs its range restated and each of its pcurves re-derived
//! against a surface whose parameterization has also moved, and re-deriving a
//! pcurve is a fit rather than a construction. See `docs/PLAN.md`.
//!
//! Because the parameterization moves, an *arc* is built as an arc rather than
//! built whole and trimmed: the span is what is converted, so the result covers
//! exactly it.

use core::f64::consts::TAU;

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Frame, KnotVector, Point, Weighted};

use crate::curve::{BSplineCurve, Curve};
use crate::fit::{self, Fitted};
use crate::traits::{Curve3d, Surface};
use ogeom_core::ogeom_err;

/// The widest span one rational quadratic Bézier is allowed to cover.
///
/// A quarter turn. The construction degrades as the span approaches half a
/// turn (the tangents meet further and further away and the weight falls to
/// zero), so the arc is split until every span is comfortably inside that.
const MAX_SPAN: f64 = core::f64::consts::FRAC_PI_2;

impl Curve {
    /// This curve as a B-spline, exactly, over `[0, 1]`.
    ///
    /// Exact rather than fitted: the result passes through the same points as
    /// the original at corresponding parameters, to rounding. The
    /// *correspondence* is not the identity; see the module documentation for
    /// why it cannot be.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the curve's
    /// range is degenerate or its geometry cannot be evaluated.
    pub fn to_bspline(&self, tol: Tolerances) -> OgeomResult<BSplineCurve> {
        let (lo, hi) = self.domain();
        self.to_bspline_over((lo, hi), tol)
    }

    /// This curve as a B-spline over part of its range.
    ///
    /// The span is what gets converted, rather than the whole curve being
    /// converted and then trimmed. For a conic those differ: trimming a
    /// rational quadratic needs knot insertion at a parameter that has to be
    /// solved for, while building the arc directly is a closed form.
    ///
    /// # Errors
    ///
    /// As [`Curve::to_bspline`].
    pub fn to_bspline_over(&self, range: (f64, f64), tol: Tolerances) -> OgeomResult<BSplineCurve> {
        let (lo, hi) = range;
        if !lo.is_finite() || !hi.is_finite() || hi <= lo + tol.parametric() {
            ogeom_bail!(Construction, "cannot convert an empty range [{lo}, {hi}]");
        }
        match self {
            // A segment is a degree-one B-spline with two control points, and
            // that is not an approximation of a line, it is a line.
            Self::Line(_) => segment(self.point_at(lo, tol)?, self.point_at(hi, tol)?),

            Self::Circle(c) => {
                let circle = c.circle();
                let (a, b) = oriented(lo, hi, c.is_reversed());
                conic_arc(circle.frame(), circle.radius(), circle.radius(), a, b, tol)
            }
            Self::Ellipse(e) => {
                let ellipse = e.ellipse();
                let (a, b) = oriented(lo, hi, e.is_reversed());
                conic_arc(
                    ellipse.frame(),
                    ellipse.major_radius(),
                    ellipse.minor_radius(),
                    a,
                    b,
                    tol,
                )
            }

            // A parabola is a quadratic, so one *polynomial* Bézier covers any
            // span of it exactly, no weights needed. The middle control point
            // is where the tangents at the ends meet.
            Self::Parabola(_) | Self::Hyperbola(_) => tangent_quadratic(self, lo, hi, tol),

            // Already one. Trimmed to the span, which for a spline is knot
            // insertion and therefore exact.
            Self::BSpline(s) => {
                let (a, b) = s.knots().domain();
                let mut out = s.clone();
                if hi < b - tol.parametric() {
                    out = out.split_at(hi, tol)?.0;
                }
                if lo > a + tol.parametric() {
                    out = out.split_at(lo, tol)?.1;
                }
                normalized(out)
            }

            // A helix is transcendental: no rational B-spline states it
            // exactly, and this function's contract is exactness. Fit it
            // through the fitting machinery at a stated tolerance instead.
            Self::Helix(_) => ogeom_bail!(
                Construction,
                "a helix has no exact B-spline form; fit it at a stated tolerance instead"
            ),

            // The same contract refuses the derived types whose spelling is
            // a quotient or a composition: nothing rational states them.
            Self::Offset(_) => ogeom_bail!(
                Construction,
                "an offset curve has no exact B-spline form; fit it at a stated tolerance instead"
            ),
            Self::OnSurface(_) => ogeom_bail!(
                Construction,
                "a surface curve has no exact B-spline form in general; fit it at a stated \
                 tolerance instead"
            ),

            Self::Trimmed(t) => {
                // A trim does not renumber anything: its parameter *is* its
                // basis's, restricted to a sub-interval. Only a reversed trim
                // moves one, and it mirrors within the trim's own range.
                let (ta, tb) = t.domain();
                if !t.is_reversed() {
                    return t.basis().to_bspline_over((lo, hi), tol);
                }
                let at = |u: f64| ta + tb - u;
                // Mirroring swaps the ends, so the basis is converted forwards
                // and the result turned round.
                let forwards = t.basis().to_bspline_over((at(hi), at(lo)), tol)?;
                reverse(&forwards)
            }
        }
    }
}

/// A degree-one B-spline between two points.
fn segment(from: Point, to: Point) -> OgeomResult<BSplineCurve> {
    let knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1)?;
    BSplineCurve::rational(
        knots,
        vec![
            Weighted {
                scaled: from,
                weight: 1.0,
            },
            Weighted {
                scaled: to,
                weight: 1.0,
            },
        ],
    )
}

/// The angular range to build, accounting for a curve that runs backwards.
///
/// A reversed conic evaluates at the *negated* angle (not at a mirrored one
/// within its range), so converting the span `[lo, hi]` of it means walking the
/// underlying conic from `-lo` to `-hi`, which runs the other way round.
const fn oriented(lo: f64, hi: f64, reversed: bool) -> (f64, f64) {
    if reversed { (-lo, -hi) } else { (lo, hi) }
}

/// A circular or elliptical arc as a piecewise rational quadratic.
///
/// One construction serves both, because an ellipse is the affine image of a
/// circle and the rational quadratic form is preserved by an affine map: the
/// control points move with it and the *weights do not change at all*. Writing
/// the ellipse case separately would be writing the same thing twice with two
/// chances to get it wrong.
fn conic_arc(
    frame: Frame,
    major: f64,
    minor: f64,
    from: f64,
    to: f64,
    tol: Tolerances,
) -> OgeomResult<BSplineCurve> {
    let sweep = to - from;
    if sweep.abs() > TAU + tol.parametric() {
        ogeom_bail!(
            Construction,
            "an arc of {sweep} radians covers the conic more than once"
        );
    }
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss,
        reason = "a span count bounded by four; the ceiling below is exact"
    )]
    let spans = ((sweep.abs() / MAX_SPAN).ceil() as usize).max(1);
    #[allow(clippy::cast_precision_loss)]
    let step = sweep / spans as f64;
    // Half the span, which is the angle between a chord and the tangent at its
    // end. Its cosine is the middle control point's weight, and the reciprocal
    // is how far along the bisector that point sits.
    let half = step * 0.5;
    let (cos_half, reach) = (half.cos(), 1.0 / half.cos());
    if cos_half <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a span of {step} radians is too wide for one rational quadratic"
        );
    }

    let at = |angle: f64| frame.to_world(Point::new(major * angle.cos(), minor * angle.sin(), 0.0));
    let shoulder = |angle: f64| {
        frame.to_world(Point::new(
            major * reach * angle.cos(),
            minor * reach * angle.sin(),
            0.0,
        ))
    };

    let mut control = Vec::with_capacity(2 * spans + 1);
    control.push(Weighted {
        scaled: at(from),
        weight: 1.0,
    });
    for k in 0..spans {
        #[allow(clippy::cast_precision_loss)]
        let start = from + step * k as f64;
        let middle = start + half;
        let end = start + step;
        // Stored homogeneous (the point already multiplied by its weight)
        // because that is the form evaluation wants and converting on the way
        // in and out again would only add rounding.
        control.push(Weighted {
            scaled: Point::from_vector(shoulder(middle).to_vector() * cos_half),
            weight: cos_half,
        });
        control.push(Weighted {
            scaled: at(end),
            weight: 1.0,
        });
    }

    let mut knots = vec![0.0, 0.0, 0.0];
    for k in 1..spans {
        #[allow(clippy::cast_precision_loss)]
        let at_knot = k as f64 / spans as f64;
        knots.push(at_knot);
        knots.push(at_knot);
    }
    knots.extend([1.0, 1.0, 1.0]);
    BSplineCurve::rational(KnotVector::new(knots, 2)?, control)
}

/// A span of a curve whose second derivative is constant in its own frame, as
/// one quadratic Bézier through the meeting point of its end tangents.
///
/// Exact for a parabola, which *is* a quadratic. For a hyperbola the same
/// construction is exact with a weight on the middle point, and the weight
/// falls out of requiring the curve to pass through its own midpoint, which is
/// what is solved for here rather than quoted from a table, so it stays right
/// for any span.
fn tangent_quadratic(
    curve: &Curve,
    lo: f64,
    hi: f64,
    tol: Tolerances,
) -> OgeomResult<BSplineCurve> {
    let (start, end) = (curve.point_at(lo, tol)?, curve.point_at(hi, tol)?);
    let (ta, tb) = (curve.d1_at(lo, tol)?, curve.d1_at(hi, tol)?);

    // Where the two end tangents meet. For a conic that point is the middle
    // control point of its quadratic form.
    let Some(shoulder) = meet(start, ta, end, tb, tol) else {
        ogeom_bail!(
            Construction,
            "the end tangents of this span are parallel, so it has no quadratic \
             form; split the range"
        );
    };

    // The weight that makes the Bézier pass through the curve's own midpoint.
    // At the Bézier's centre the three basis values are 1/4, 1/2, 1/4, so the
    // point there is (A + 2wS + B) / (2 + 2w), and solving for w against the
    // real midpoint is one division.
    let middle = curve.point_at(f64::midpoint(lo, hi), tol)?;
    let numerator = (start.to_vector() + end.to_vector()) * 0.5 - middle.to_vector();
    let denominator = middle.to_vector() - shoulder.to_vector();
    let weight = if denominator.magnitude() <= tol.confusion() {
        1.0
    } else {
        let w = numerator.dot(denominator) / denominator.square_magnitude();
        if !w.is_finite() || w <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "this span has no rational quadratic form; split the range"
            );
        }
        w
    };

    let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2)?;
    BSplineCurve::rational(
        knots,
        vec![
            Weighted {
                scaled: start,
                weight: 1.0,
            },
            Weighted {
                scaled: Point::from_vector(shoulder.to_vector() * weight),
                weight,
            },
            Weighted {
                scaled: end,
                weight: 1.0,
            },
        ],
    )
}

/// Where two lines meet, or `None` if they are parallel or skew.
fn meet(
    a: Point,
    along_a: ogeom_math::Vector,
    b: Point,
    along_b: ogeom_math::Vector,
    tol: Tolerances,
) -> Option<Point> {
    let between = b - a;
    let cross = along_a.cross(along_b);
    let denominator = cross.square_magnitude();
    if denominator <= tol.confusion() * tol.confusion() {
        return None;
    }
    let t = between.cross(along_b).dot(cross) / denominator;
    let found = a + along_a * t;
    // Skew lines have a nearest approach rather than a meeting; only a real
    // intersection is a control point.
    let s = between.cross(along_a).dot(cross) / denominator;
    if found.distance(b + along_b * s) > tol.confusion() {
        return None;
    }
    Some(found)
}

/// A spline traced the other way.
///
/// The control points reverse and the knots mirror within their own span. No
/// geometry moves; this is the same curve, walked backwards.
fn reverse(curve: &BSplineCurve) -> OgeomResult<BSplineCurve> {
    let (a, b) = curve.knots().domain();
    let mut knots: Vec<f64> = curve.knots().knots().iter().map(|k| a + b - k).collect();
    knots.reverse();
    let mut control = curve.control_points().to_vec();
    control.reverse();
    BSplineCurve::rational(KnotVector::new(knots, curve.degree())?, control)
}

/// A spline over `[0, 1]`, whatever it was over.
impl Curve {
    /// This curve as a B-spline *fitted* over `range` to a stated
    /// tolerance, at the curve's own parameters.
    ///
    /// The approximation [`Curve::to_bspline`] refuses to make silently:
    /// a helix, an offset curve, a curve on a surface (anything with no
    /// exact rational form) is sampled at its own parameters and fitted
    /// through them, the fit measured against the curve *between* the
    /// samples as well as at them, and the sampling doubled until the
    /// tolerance is met or the budget runs out. The result is
    /// same-parameter with the curve, so every chart already speaking the
    /// curve's parameter still does, and its `error` is what was measured,
    /// not what was asked.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
    /// the range is empty or the tolerance is not a distance; as
    /// [`fit::fit_points_at`].
    pub fn fitted_bspline_over(
        &self,
        range: (f64, f64),
        tolerance: f64,
        tol: Tolerances,
    ) -> OgeomResult<Fitted<BSplineCurve>> {
        same_parameter_fit(self, range, 3, tolerance, tol)
    }
}

impl BSplineCurve {
    /// This curve at a degree no higher than `max_degree`, to a stated
    /// tolerance: itself where it already is, and otherwise a fit at its
    /// own parameters: what an exchange format with a degree limit needs
    /// written. Same-parameter with the original, error measured.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
    /// `max_degree` is zero, the curve is periodic (a ring has no ends to
    /// fit between), or the tolerance is not a distance.
    pub fn restricted_to_degree(
        &self,
        max_degree: usize,
        tolerance: f64,
        tol: Tolerances,
    ) -> OgeomResult<Fitted<Self>> {
        if max_degree == 0 {
            ogeom_bail!(Construction, "a curve needs a degree of at least one");
        }
        if self.degree() <= max_degree {
            return Ok(Fitted {
                curve: self.clone(),
                error: 0.0,
                met: true,
            });
        }
        if self.is_periodic() {
            ogeom_bail!(
                Construction,
                "a periodic curve has no ends to fit between; reseam it first"
            );
        }
        let domain = self.domain();
        same_parameter_fit(
            &Curve::BSpline(self.clone()),
            domain,
            max_degree,
            tolerance,
            tol,
        )
    }
}

impl crate::surface::BSplineSurface {
    /// This patch at degrees no higher than `max_degree` in either
    /// direction, to a stated tolerance: itself where it already is, and
    /// otherwise a fit through a grid of its own points, the grid doubled
    /// until the fit holds every sample to the tolerance or the budget
    /// runs out. The fit's parameterization is its own (chord-length
    /// through the grid, not the patch's), so a pcurve spoken against the
    /// patch must be re-derived against the result.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
    /// `max_degree` is zero or the tolerance is not a distance; as
    /// [`fit::fit_surface_grid`].
    pub fn restricted_to_degree(
        &self,
        max_degree: usize,
        tolerance: f64,
        tol: Tolerances,
    ) -> OgeomResult<Fitted<Self>> {
        if max_degree == 0 {
            ogeom_bail!(Construction, "a patch needs a degree of at least one");
        }
        if self.u_knots().degree() <= max_degree && self.v_knots().degree() <= max_degree {
            return Ok(Fitted {
                curve: self.clone(),
                error: 0.0,
                met: true,
            });
        }
        grid_fitted(
            |u, v| self.point_at(u, v, tol),
            self.domain(),
            max_degree,
            tolerance,
            tol,
        )
    }
}

/// A patch fitted at `degree` through a grid of `point` over `domain`, the
/// grid doubled until the fit holds every sample to `tolerance` or the
/// budget runs out; the best fit either way, its error measured at the
/// samples.
fn grid_fitted(
    point: impl Fn(f64, f64) -> OgeomResult<Point>,
    domain: ((f64, f64), (f64, f64)),
    degree: usize,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Fitted<crate::surface::BSplineSurface>> {
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        ogeom_bail!(Construction, "a tolerance of {tolerance} is not a distance");
    }
    let ((ua, ub), (va, vb)) = domain;
    if ![ua, ub, va, vb].iter().all(|x| x.is_finite()) {
        ogeom_bail!(Construction, "an unbounded surface cannot be fitted");
    }
    let mut samples = 16usize;
    let mut best: Option<Fitted<crate::surface::BSplineSurface>> = None;
    for _ in 0..4 {
        let mut rows: Vec<Vec<Point>> = Vec::with_capacity(samples + 1);
        for j in 0..=samples {
            #[allow(clippy::cast_precision_loss)]
            let v = va + (vb - va) * j as f64 / samples as f64;
            let mut row = Vec::with_capacity(samples + 1);
            for i in 0..=samples {
                #[allow(clippy::cast_precision_loss)]
                let u = ua + (ub - ua) * i as f64 / samples as f64;
                row.push(point(u, v)?);
            }
            rows.push(row);
        }
        let fitted = fit::fit_surface_grid(&rows, degree, tolerance, tol)?;
        if fitted.met {
            return Ok(fitted);
        }
        if best.as_ref().is_none_or(|b| fitted.error < b.error) {
            best = Some(fitted);
        }
        samples *= 2;
    }
    best.ok_or_else(|| ogeom_err!(Construction, "the patch could not be sampled"))
}

impl crate::surface::SurfaceGeometry {
    /// This surface as a B-spline *fitted* over its own domain to a stated
    /// tolerance: the approximation [`to_bspline`](Self::to_bspline)
    /// refuses to make silently, for an offset surface or anything else
    /// with no exact rational form. A grid of the surface's points is
    /// fitted at degree three, the grid doubled until every sample is
    /// within the tolerance or the budget runs out, and `error` is what
    /// was measured. The fit's parameterization is its own, so a pcurve
    /// spoken against the surface must be re-derived against the result.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
    /// the domain is unbounded or the tolerance is not a distance; as
    /// [`fit::fit_surface_grid`].
    pub fn fitted_bspline(
        &self,
        tolerance: f64,
        tol: Tolerances,
    ) -> OgeomResult<Fitted<crate::surface::BSplineSurface>> {
        use crate::traits::Surface as _;
        grid_fitted(
            |u, v| self.point_at(u, v, tol),
            self.domain(),
            3,
            tolerance,
            tol,
        )
    }
}

/// A fit of `curve` over `range` at the curve's own parameters, measured
/// between the samples as well as at them, the sampling doubled until the
/// tolerance is met or the budget runs out.
fn same_parameter_fit(
    curve: &Curve,
    range: (f64, f64),
    degree: usize,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Fitted<BSplineCurve>> {
    let (lo, hi) = range;
    if !lo.is_finite() || !hi.is_finite() || hi <= lo + tol.parametric() {
        ogeom_bail!(Construction, "cannot fit over an empty range [{lo}, {hi}]");
    }
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        ogeom_bail!(Construction, "a tolerance of {tolerance} is not a distance");
    }
    let mut samples = 32usize;
    let mut best: Option<Fitted<BSplineCurve>> = None;
    for _ in 0..8 {
        #[allow(clippy::cast_precision_loss)]
        let params: Vec<f64> = (0..=samples)
            .map(|i| lo + (hi - lo) * i as f64 / samples as f64)
            .collect();
        let points = params
            .iter()
            .map(|t| curve.point_at(*t, tol))
            .collect::<OgeomResult<Vec<Point>>>()?;
        let fitted = fit::fit_points_at(&params, &points, degree, tolerance, tol)?;
        let mut error = fitted.error;
        for pair in params.windows(2) {
            let t = f64::midpoint(pair[0], pair[1]);
            error = error.max(
                curve
                    .point_at(t, tol)?
                    .distance(fitted.curve.point_at(t, tol)?),
            );
        }
        let candidate = Fitted {
            curve: fitted.curve,
            error,
            met: error <= tolerance,
        };
        if candidate.met {
            return Ok(candidate);
        }
        if best.as_ref().is_none_or(|b| error < b.error) {
            best = Some(candidate);
        }
        samples *= 2;
    }
    best.ok_or_else(|| ogeom_err!(Construction, "the curve could not be sampled"))
}

fn normalized(curve: BSplineCurve) -> OgeomResult<BSplineCurve> {
    let knots = curve.knots().reparameterized(0.0, 1.0)?;
    BSplineCurve::rational(knots, curve.control_points().to_vec())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use core::f64::consts::{FRAC_PI_2, PI};
    use ogeom_math::{Circle, Direction, Ellipse, Hyperbola, Parabola, Vector};

    use crate::curve::{
        CircleCurve, EllipseCurve, HyperbolaCurve, LineCurve, ParabolaCurve, TrimmedCurve,
    };

    const T: Tolerances = Tolerances::millimetres();

    /// A helix has no exact form and the exact conversion still says so;
    /// asked for a fit at a stated tolerance, it comes back same-parameter
    /// with the helix, within the tolerance everywhere, the error reported
    /// as measured.
    #[test]
    fn a_helix_fits_to_a_stated_tolerance_at_its_own_parameters() {
        let helix: Curve = crate::curve::HelixCurve::new(ogeom_math::Frame::WORLD, 5.0, 4.0, 2.0)
            .unwrap()
            .into();
        assert!(helix.to_bspline(T).is_err(), "the exact conversion refuses");
        let domain = helix.domain();
        let fitted = helix.fitted_bspline_over(domain, 1e-3, T).unwrap();
        assert!(fitted.met && fitted.error <= 1e-3, "error {}", fitted.error);
        let (lo, hi) = fitted.curve.domain();
        assert!((lo - domain.0).abs() < 1e-12 && (hi - domain.1).abs() < 1e-12);
        for i in 0..=200 {
            let t = domain.0 + (domain.1 - domain.0) * f64::from(i) / 200.0;
            let gap = helix
                .point_at(t, T)
                .unwrap()
                .distance(fitted.curve.point_at(t, T).unwrap());
            assert!(gap <= 1e-3, "same-parameter within tolerance at {t}: {gap}");
        }
    }

    /// A cubic raised to a quintic restricted back to degree three is the
    /// cubic again, to rounding; one already low enough is itself.
    #[test]
    fn a_curve_restricted_in_degree_holds_its_tolerance() {
        let cubic = BSplineCurve::new(
            ogeom_math::KnotVector::clamped_uniform(3, 6).unwrap(),
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.0, 2.0, 0.5),
                Point::new(2.5, 1.0, -0.5),
                Point::new(4.0, 3.0, 1.0),
                Point::new(5.0, 0.5, 0.0),
                Point::new(6.0, 2.0, 2.0),
            ],
            T,
        )
        .unwrap();
        let quintic = cubic.elevated(T).unwrap().elevated(T).unwrap();
        assert_eq!(quintic.degree(), 5);
        let same = quintic.restricted_to_degree(5, 1e-6, T).unwrap();
        assert!(same.met && same.error == 0.0 && same.curve.degree() == 5);
        let back = quintic.restricted_to_degree(3, 1e-6, T).unwrap();
        assert!(back.met && back.curve.degree() == 3, "error {}", back.error);
        let (lo, hi) = cubic.domain();
        for i in 0..=50 {
            let t = lo + (hi - lo) * f64::from(i) / 50.0;
            let gap = cubic
                .point_at(t, T)
                .unwrap()
                .distance(back.curve.point_at(t, T).unwrap());
            assert!(gap < 1e-6, "the cubic again at {t}: {gap}");
        }
    }

    /// A bicubic patch restricted to degree two holds every sample of
    /// itself to the tolerance asked, at the lower degree both ways.
    #[test]
    fn a_patch_restricted_in_degree_holds_its_samples() {
        let (nu, nv) = (5, 4);
        let mut points = Vec::with_capacity(nu * nv);
        for i in 0..nu {
            for j in 0..nv {
                #[allow(clippy::cast_precision_loss)]
                let (x, y) = (i as f64, j as f64);
                points.push(Point::new(x, y, (x * 0.7).sin() * (y * 0.5).cos() * 0.3));
            }
        }
        let patch = crate::surface::BSplineSurface::new(
            ogeom_math::KnotVector::clamped_uniform(3, nu).unwrap(),
            ogeom_math::KnotVector::clamped_uniform(3, nv).unwrap(),
            &ogeom_math::ControlGrid::new(points, nu, nv).unwrap(),
            T,
        )
        .unwrap();
        let lower = patch.restricted_to_degree(2, 1e-2, T).unwrap();
        assert!(lower.met && lower.error <= 1e-2, "error {}", lower.error);
        assert!(lower.curve.u_knots().degree() == 2 && lower.curve.v_knots().degree() == 2);
        let same = patch.restricted_to_degree(3, 1e-2, T).unwrap();
        assert!(same.met && same.error == 0.0);
    }

    /// The greatest distance from any point of the converted curve to the
    /// original, sampled densely.
    ///
    /// The parameterizations differ, so the comparison is *geometric*: for each
    /// sample of the conversion, how far is it from the curve it came from. For
    /// a conic that distance has a closed form, so no search is involved.
    fn deviation(original: &Curve, converted: &BSplineCurve, samples: usize) -> f64 {
        let (a, b) = converted.knots().domain();
        let mut worst = 0.0_f64;
        for i in 0..=samples {
            #[allow(clippy::cast_precision_loss)]
            let u = a + (b - a) * i as f64 / samples as f64;
            let p = converted.point_at(u, T).unwrap();
            worst = worst.max(distance_to(original, p));
        }
        worst
    }

    /// Distance from a point to an analytic curve, in closed form.
    fn distance_to(curve: &Curve, p: Point) -> f64 {
        match curve {
            Curve::Circle(c) => c.circle().distance_to(p),
            Curve::Ellipse(e) => {
                // Not closed form, so fall back to a dense parameter search.
                nearest(curve, p, e.ellipse().major_radius())
            }
            _ => nearest(curve, p, 1.0),
        }
    }

    /// Distance from a point to a curve, by a coarse scan and then a local
    /// refinement.
    ///
    /// The scan alone is not good enough to measure a conversion by: its error
    /// is set by the step it takes, so a test built on it reports the sampling
    /// resolution rather than the conversion's accuracy. Refining around the
    /// best sample takes it to where the answer is about the curve again.
    fn nearest(curve: &Curve, p: Point, _scale: f64) -> f64 {
        const SCAN: usize = 2_000;
        let (a, b) = curve.domain();
        let at = |u: f64| curve.point_at(u, T).map_or(f64::MAX, |q| p.distance(q));

        let mut best = (a, f64::MAX);
        for i in 0..=SCAN {
            #[allow(clippy::cast_precision_loss)]
            let u = a + (b - a) * i as f64 / SCAN as f64;
            let d = at(u);
            if d < best.1 {
                best = (u, d);
            }
        }
        // Ternary search on the bracket either side of the best sample. The
        // distance is unimodal there for every curve this is used on.
        #[allow(clippy::cast_precision_loss)]
        let step = (b - a) / SCAN as f64;
        let (mut lo, mut hi) = (best.0 - step, best.0 + step);
        for _ in 0..200 {
            let one = lo + (hi - lo) / 3.0;
            let two = hi - (hi - lo) / 3.0;
            if at(one) < at(two) {
                hi = two;
            } else {
                lo = one;
            }
        }
        at(f64::midpoint(lo, hi)).min(best.1)
    }

    #[test]
    fn a_line_becomes_a_degree_one_spline_through_its_own_ends() {
        let from = Point::new(1.0, 2.0, 3.0);
        let to = Point::new(4.0, -1.0, 0.5);
        let line: Curve = LineCurve::segment(from, to, T).unwrap().into();
        let spline = line.to_bspline(T).unwrap();

        assert_eq!(spline.degree(), 1);
        assert!(!spline.is_rational(), "a line needs no weights");
        assert!(spline.point_at(0.0, T).unwrap().is_equal(from, T));
        assert!(spline.point_at(1.0, T).unwrap().is_equal(to, T));
        assert!(spline.point_at(0.5, T).unwrap().is_equal(
            Point::from_vector((from.to_vector() + to.to_vector()) * 0.5),
            T
        ));
    }

    #[test]
    fn a_full_circle_becomes_an_exact_rational_quadratic() {
        // Exact, not fitted: every sample lands on the circle to rounding. This
        // is the property that makes rational weights load-bearing: a
        // polynomial spline cannot represent a circle at all, only approach it.
        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 2.5, T).unwrap()).into();
        let spline = circle.to_bspline(T).unwrap();

        assert_eq!(spline.degree(), 2);
        assert!(spline.is_rational(), "a circle needs its weights");
        assert!(
            deviation(&circle, &spline, 500) < 1e-12,
            "off the circle by {}",
            deviation(&circle, &spline, 500)
        );
        // And it closes.
        assert!(
            spline
                .point_at(0.0, T)
                .unwrap()
                .is_equal(spline.point_at(1.0, T).unwrap(), T)
        );
    }

    #[test]
    fn an_arc_covers_its_own_span_and_no_more() {
        // Built as an arc rather than built whole and trimmed, so the ends are
        // the arc's ends exactly.
        let circle = Circle::new(Frame::WORLD, 3.0, T).unwrap();
        let curve: Curve = CircleCurve::new(circle).into();
        for (from, to) in [
            (0.0, FRAC_PI_2),
            (0.3, 1.9),
            (PI, PI * 1.5),
            (0.0, PI * 1.75),
        ] {
            let spline = curve.to_bspline_over((from, to), T).unwrap();
            assert!(
                spline
                    .point_at(0.0, T)
                    .unwrap()
                    .is_equal(curve.point_at(from, T).unwrap(), T),
                "arc [{from}, {to}] starts in the wrong place"
            );
            assert!(
                spline
                    .point_at(1.0, T)
                    .unwrap()
                    .is_equal(curve.point_at(to, T).unwrap(), T),
                "arc [{from}, {to}] ends in the wrong place"
            );
            assert!(deviation(&curve, &spline, 200) < 1e-12);
        }
    }

    #[test]
    fn an_ellipse_uses_the_same_construction_and_the_same_weights() {
        // An ellipse is the affine image of a circle, and an affine map carries
        // a rational quadratic to a rational quadratic with the weights
        // untouched. One routine, not two.
        let ellipse: Curve =
            EllipseCurve::new(Ellipse::new(Frame::WORLD, 5.0, 2.0, T).unwrap()).into();
        let spline = ellipse.to_bspline(T).unwrap();
        assert_eq!(spline.degree(), 2);
        assert!(deviation(&ellipse, &spline, 400) < 1e-9);

        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 5.0, T).unwrap()).into();
        let round = circle.to_bspline(T).unwrap();
        let weights: Vec<f64> = spline.control_points().iter().map(|c| c.weight).collect();
        let same: Vec<f64> = round.control_points().iter().map(|c| c.weight).collect();
        assert_eq!(weights, same, "the weights should not depend on the radii");
    }

    #[test]
    fn a_conic_in_a_tilted_frame_converts_where_it_actually_is() {
        let frame = Frame::new(
            Point::new(3.0, -2.0, 7.0),
            Direction::from_coords(1.0, 1.0, 1.0, T).unwrap(),
            Direction::from_coords(1.0, -1.0, 0.0, T).unwrap(),
            T,
        )
        .unwrap();
        let circle: Curve = CircleCurve::new(Circle::new(frame, 4.0, T).unwrap()).into();
        let spline = circle.to_bspline(T).unwrap();
        assert!(deviation(&circle, &spline, 300) < 1e-12);
    }

    #[test]
    fn a_parabola_is_a_polynomial_quadratic_exactly() {
        // Degree two and no weights: a parabola *is* a quadratic, so the
        // conversion is not even rational.
        let parabola: Curve = ParabolaCurve::new(Parabola::new(Frame::WORLD, 1.5, T).unwrap(), 4.0)
            .unwrap()
            .into();
        let spline = parabola.to_bspline(T).unwrap();
        assert_eq!(spline.degree(), 2);
        assert!(deviation(&parabola, &spline, 300) < 1e-9);
        for c in spline.control_points() {
            assert!(
                (c.weight - 1.0).abs() < 1e-9,
                "a parabola should need no weights, got {}",
                c.weight
            );
        }
    }

    #[test]
    fn a_hyperbola_becomes_a_rational_quadratic() {
        let hyperbola: Curve =
            HyperbolaCurve::new(Hyperbola::new(Frame::WORLD, 2.0, 1.0, T).unwrap(), 1.0)
                .unwrap()
                .into();
        let spline = hyperbola.to_bspline(T).unwrap();
        assert_eq!(spline.degree(), 2);
        assert!(
            deviation(&hyperbola, &spline, 300) < 1e-9,
            "off by {}",
            deviation(&hyperbola, &spline, 300)
        );
    }

    #[test]
    fn a_spline_converts_to_itself_and_a_trimmed_one_to_its_piece() {
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0], 2).unwrap();
        let control = vec![
            Point::ORIGIN,
            Point::new(1.0, 2.0, 0.0),
            Point::new(3.0, 2.0, 1.0),
            Point::new(4.0, 0.0, 0.0),
        ];
        let spline: Curve = BSplineCurve::new(knots, control, T).unwrap().into();
        let same = spline.to_bspline(T).unwrap();
        assert!(deviation(&spline, &same, 200) < 1e-12);

        let trimmed: Curve = Curve::Trimmed(Box::new(
            TrimmedCurve::new(spline.clone(), 0.25, 0.75, T).unwrap(),
        ));
        let piece = trimmed.to_bspline(T).unwrap();
        assert!(
            piece
                .point_at(0.0, T)
                .unwrap()
                .is_equal(spline.point_at(0.25, T).unwrap(), T)
        );
        assert!(
            piece
                .point_at(1.0, T)
                .unwrap()
                .is_equal(spline.point_at(0.75, T).unwrap(), T)
        );
        assert!(deviation(&spline, &piece, 200) < 1e-12);
    }

    #[test]
    fn a_reversed_conic_converts_to_the_curve_it_actually_traces() {
        use crate::traits::Reversible;
        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 1.0, T).unwrap()).into();
        let backwards = circle.reversed();
        let spline = backwards.to_bspline_over((0.0, FRAC_PI_2), T).unwrap();

        assert!(
            spline
                .point_at(0.0, T)
                .unwrap()
                .is_equal(backwards.point_at(0.0, T).unwrap(), T),
            "a reversed arc should start where the reversed curve does"
        );
        assert!(
            spline
                .point_at(1.0, T)
                .unwrap()
                .is_equal(backwards.point_at(FRAC_PI_2, T).unwrap(), T)
        );
    }

    #[test]
    fn an_empty_range_is_refused() {
        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 1.0, T).unwrap()).into();
        assert!(circle.to_bspline_over((1.0, 1.0), T).is_err());
        assert!(circle.to_bspline_over((1.0, 0.0), T).is_err());
        assert!(circle.to_bspline_over((0.0, f64::NAN), T).is_err());
        let _ = Vector::ZERO;
    }
}

// --- surfaces ----------------------------------------------------------------

impl crate::surface::SurfaceGeometry {
    /// This surface as a B-spline patch, exactly, over `[0, 1]` in both
    /// directions.
    ///
    /// Every analytic surface here is a surface of revolution or a ruled one,
    /// so its exact form falls out of the curve conversion above rather than
    /// needing a construction of its own: revolve an exactly-converted profile
    /// and the patch is exact wherever the profile was.
    ///
    /// As for curves, the *parameterization* does not survive: a cylinder's
    /// `u` is an angle and the patch's is not. See the module documentation.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// surface's own extents are degenerate, or its geometry cannot be
    /// evaluated.
    pub fn to_bspline(&self, tol: Tolerances) -> OgeomResult<crate::surface::BSplineSurface> {
        self.to_bspline_over(self.domain(), tol)
    }

    /// This surface over `window` of its own parameters as a B-spline patch,
    /// exactly: [`to_bspline`](Self::to_bspline) for a window other than the
    /// surface's own extents. A trimmed surface converts as its basis over
    /// the trim's window.
    ///
    /// # Errors
    ///
    /// As [`to_bspline`](Self::to_bspline).
    fn to_bspline_over(
        &self,
        window: ((f64, f64), (f64, f64)),
        tol: Tolerances,
    ) -> OgeomResult<crate::surface::BSplineSurface> {
        use crate::surface::SurfaceGeometry as S;
        let ((ua, ub), (va, vb)) = window;
        match self {
            // An offset of a free-form basis has no exact spline form
            // (the unit normal is a quotient), and this function's contract
            // is exactness.
            S::Offset(_) => ogeom_bail!(
                Construction,
                "an offset surface has no exact B-spline form; offset the basis analytically \
                 or fit at a stated tolerance instead"
            ),

            // Ruled both ways: four corners and a bilinear patch, which is a
            // degree-one B-spline in each direction and exact for a plane.
            S::Plane(_) => {
                let corner = |u: f64, v: f64| self.point_at(u, v, tol);
                bilinear(
                    corner(ua, va)?,
                    corner(ub, va)?,
                    corner(ua, vb)?,
                    corner(ub, vb)?,
                )
            }

            // Ruled along `v`: the profile at each end of the height, lofted.
            // Exact because a cylinder and a cone are straight in `v`.
            S::Cylinder(_) | S::Cone(_) => {
                let profile = |v: f64| -> OgeomResult<Vec<Weighted<Point>>> {
                    let mut out = Vec::new();
                    let ring = self.section_at(v, (ua, ub), tol)?;
                    out.extend_from_slice(ring.control_points());
                    Ok(out)
                };
                let (knots, _) = {
                    let ring = self.section_at(va, (ua, ub), tol)?;
                    (ring.knots().clone(), ())
                };
                loft(&knots, &profile(va)?, &profile(vb)?)
            }

            // Closed in both directions and revolved: build the profile once
            // and turn it, which is the same construction the surface is.
            S::Sphere(sp) => {
                let sphere = sp.sphere();
                let frame = sphere.frame();
                // The meridian at u = 0, framed so its own angle is the
                // latitude exactly: x one radius out along the equator, y
                // toward the north pole.
                let meridian = Frame::new(frame.origin(), -frame.y(), frame.x(), tol)?;
                let circle = ogeom_math::Circle::new(meridian, sphere.radius(), tol)?;
                let profile: Curve = crate::curve::CircleCurve::new(circle).into();
                revolved_patch(
                    &profile,
                    (va, vb),
                    ogeom_math::Axis {
                        location: frame.origin(),
                        direction: frame.z(),
                    },
                    (ua, ub),
                    tol,
                )
            }
            S::Torus(t) => {
                let torus = t.torus();
                let frame = torus.frame();
                let centre = frame.origin() + frame.x().vector() * torus.major_radius();
                let tube = Frame::new(centre, -frame.y(), frame.x(), tol)?;
                let circle = ogeom_math::Circle::new(tube, torus.minor_radius(), tol)?;
                let profile: Curve = crate::curve::CircleCurve::new(circle).into();
                revolved_patch(
                    &profile,
                    (va, vb),
                    ogeom_math::Axis {
                        location: frame.origin(),
                        direction: frame.z(),
                    },
                    (ua, ub),
                    tol,
                )
            }
            S::Revolution(r) => revolved_patch(r.curve(), (va, vb), r.axis(), (ua, ub), tol),

            S::Extrusion(e) => {
                let base = e.curve().to_bspline_over((ua, ub), tol)?;
                let along = e.direction().vector() * (vb - va);
                let start: Vec<Weighted<Point>> = base
                    .control_points()
                    .iter()
                    .map(|c| Weighted {
                        scaled: Point::from_vector(
                            c.scaled.to_vector() + e.direction().vector() * va * c.weight,
                        ),
                        weight: c.weight,
                    })
                    .collect();
                let end: Vec<Weighted<Point>> = start
                    .iter()
                    .map(|c| Weighted {
                        scaled: Point::from_vector(c.scaled.to_vector() + along * c.weight),
                        weight: c.weight,
                    })
                    .collect();
                loft(base.knots(), &start, &end)
            }

            S::BSpline(s) => {
                let whole = s.domain();
                let near = |a: (f64, f64), b: (f64, f64)| {
                    (a.0 - b.0).abs() <= tol.parametric() && (a.1 - b.1).abs() <= tol.parametric()
                };
                if near(whole.0, (ua, ub)) && near(whole.1, (va, vb)) {
                    Ok(s.clone())
                } else {
                    s.segment((ua, ub), (va, vb), tol)
                }
            }

            S::Trimmed(t) => t.basis().to_bspline_over(window, tol),
        }
    }

    /// The profile of a surface of revolution at one height, as a spline.
    fn section_at(
        &self,
        v: f64,
        u_range: (f64, f64),
        tol: Tolerances,
    ) -> OgeomResult<crate::curve::BSplineCurve> {
        use crate::surface::SurfaceGeometry as S;
        let (frame, radius) = match self {
            S::Cylinder(c) => (c.cylinder().frame(), c.cylinder().radius()),
            S::Cone(c) => (c.cone().frame(), c.cone().radius_at(v).abs()),
            _ => ogeom_bail!(Construction, "this surface has no circular section"),
        };
        let at = Frame::new(frame.origin() + frame.z() * v, frame.z(), frame.x(), tol)?;
        conic_arc(at, radius, radius, u_range.0, u_range.1, tol)
    }
}

/// A rational profile revolved about an axis, exactly.
///
/// The construction every surface of revolution shares. The profile converts
/// to its exact rational form; the turn is the exact rational unit circle
/// over the swept angle; and each surface control point is a profile control
/// point *rotated by an arc control point*, with the weights multiplied. The
/// identity behind it: rotation about the axis is linear in `(cos u, sin u)`,
/// the tensor product factors, and the patch evaluates to precisely
/// `Rot_u(profile(v))` wherever both conversions were exact.
fn revolved_patch(
    profile: &Curve,
    v_range: (f64, f64),
    axis: ogeom_math::Axis,
    u_range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<crate::surface::BSplineSurface> {
    let frame = Frame::about(axis.location, axis.direction);
    let turn: Curve =
        crate::curve::CircleCurve::new(ogeom_math::Circle::new(frame, 1.0, tol)?).into();
    let arc = turn.to_bspline_over(u_range, tol)?;
    let pro = profile.to_bspline_over(v_range, tol)?;

    let locals: Vec<(f64, f64, f64, f64)> = pro
        .control_points()
        .iter()
        .map(|c| {
            let l = frame.to_local(Point::from_vector(c.scaled.to_vector() / c.weight));
            (l.x, l.y, l.z, c.weight)
        })
        .collect();
    let mut points = Vec::with_capacity(arc.control_points().len() * locals.len());
    for ci in arc.control_points() {
        let l = frame.to_local(Point::from_vector(ci.scaled.to_vector() / ci.weight));
        let (a, b) = (l.x, l.y);
        for &(x, y, z, wj) in &locals {
            let weight = ci.weight * wj;
            let rotated =
                frame.to_world(Point::new(a.mul_add(x, -(b * y)), b.mul_add(x, a * y), z));
            points.push(Weighted {
                scaled: Point::from_vector(rotated.to_vector() * weight),
                weight,
            });
        }
    }
    let grid = ogeom_math::ControlGrid::new(points, arc.control_points().len(), locals.len())?;
    crate::surface::BSplineSurface::rational(arc.knots().clone(), pro.knots().clone(), grid)
}

/// A bilinear patch through four corners.
fn bilinear(a: Point, b: Point, c: Point, d: Point) -> OgeomResult<crate::surface::BSplineSurface> {
    let line = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1)?;
    let grid = ogeom_math::ControlGrid::new(
        vec![
            Weighted {
                scaled: a,
                weight: 1.0,
            },
            Weighted {
                scaled: c,
                weight: 1.0,
            },
            Weighted {
                scaled: b,
                weight: 1.0,
            },
            Weighted {
                scaled: d,
                weight: 1.0,
            },
        ],
        2,
        2,
    )?;
    crate::surface::BSplineSurface::rational(line.clone(), line, grid)
}

/// A patch ruled between two rows of control points sharing one knot vector.
fn loft(
    across: &KnotVector,
    start: &[Weighted<Point>],
    end: &[Weighted<Point>],
) -> OgeomResult<crate::surface::BSplineSurface> {
    if start.len() != end.len() {
        ogeom_bail!(
            Dimension,
            "a ruled patch needs the same control points at each end, got {} \
             and {}",
            start.len(),
            end.len()
        );
    }
    let mut points = Vec::with_capacity(start.len() * 2);
    for (a, b) in start.iter().zip(end) {
        points.push(*a);
        points.push(*b);
    }
    let grid = ogeom_math::ControlGrid::new(points, start.len(), 2)?;
    let along = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1)?;
    crate::surface::BSplineSurface::rational(across.clone(), along, grid)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod surface_tests {
    use super::*;
    use crate::surface::{
        ConeSurface, CylinderSurface, ExtrusionSurface, PlaneSurface, SphereSurface,
        SurfaceGeometry, TorusSurface, TrimmedSurface,
    };
    use ogeom_math::{Circle, Cone, Cylinder, Direction, Plane, Sphere, Torus};

    const T: Tolerances = Tolerances::millimetres();

    /// How far the converted patch strays from the surface it came from.
    ///
    /// Measured *implicitly* (the distance from each sampled point of the
    /// patch to the analytic surface) rather than by comparing the two at
    /// proportional parameters. Comparing parameters would be measuring the
    /// wrong thing: a rational quadratic's parameter is not proportional to the
    /// angle it sweeps, so even an exact conversion disagrees pointwise, and a
    /// test built that way reports the reparameterization instead of the error.
    fn deviation(distance: impl Fn(Point) -> f64, patch: &crate::surface::BSplineSurface) -> f64 {
        let ((pa, pb), (qa, qb)) = patch.domain();
        let mut worst = 0.0_f64;
        for i in 0..=60 {
            for j in 0..=60 {
                #[allow(clippy::cast_precision_loss)]
                let (s, t) = (i as f64 / 60.0, j as f64 / 60.0);
                if let Ok(p) = patch.point_at(pa + (pb - pa) * s, qa + (qb - qa) * t, T) {
                    worst = worst.max(distance(p).abs());
                }
            }
        }
        worst
    }

    /// Whether a patch reaches the same corners as the surface it converted.
    ///
    /// The implicit test says the patch lies *on* the surface; this says it
    /// covers the same piece of it, which the implicit test alone cannot.
    fn spans_the_same(original: &SurfaceGeometry, patch: &crate::surface::BSplineSurface) -> bool {
        let ((ua, ub), (va, vb)) = original.domain();
        let ((pa, pb), (qa, qb)) = patch.domain();
        [
            (ua, va, pa, qa),
            (ua, vb, pa, qb),
            (ub, va, pb, qa),
            (ub, vb, pb, qb),
        ]
        .iter()
        .all(|(u, v, p, q)| {
            match (original.point_at(*u, *v, T), patch.point_at(*p, *q, T)) {
                (Ok(a), Ok(b)) => a.is_equal(b, T),
                _ => false,
            }
        })
    }

    #[test]
    fn a_plane_becomes_a_bilinear_patch() {
        let plane: SurfaceGeometry =
            PlaneSurface::over(Plane::new(Frame::WORLD), (-2.0, 5.0), (-1.0, 3.0))
                .unwrap()
                .into();
        let patch = plane.to_bspline(T).unwrap();
        assert!(!patch.is_rational(), "a plane needs no weights");
        let flat = Plane::new(Frame::WORLD);
        assert!(deviation(|p| flat.distance_to(p), &patch) < 1e-12);
        assert!(spans_the_same(&plane, &patch));
    }

    #[test]
    fn a_cylinder_becomes_an_exact_rational_patch() {
        // Circular in `u` and straight in `v`, so the exact patch is the exact
        // circle lofted, and it lands on the cylinder everywhere, not near it.
        let cylinder: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 2.0, T).unwrap(), (0.0, 5.0))
                .unwrap()
                .into();
        let patch = cylinder.to_bspline(T).unwrap();
        assert!(patch.is_rational());
        let exact = Cylinder::new(Frame::WORLD, 2.0, T).unwrap();
        let off = deviation(|p| exact.distance_to(p), &patch);
        assert!(off < 1e-12, "off the cylinder by {off}");
        assert!(spans_the_same(&cylinder, &patch));
    }

    #[test]
    fn a_cone_becomes_an_exact_rational_patch() {
        let cone: SurfaceGeometry = ConeSurface::new(
            Cone::new(Frame::WORLD, 1.0, 0.4_f64.atan(), T).unwrap(),
            (0.0, 4.0),
        )
        .unwrap()
        .into();
        let patch = cone.to_bspline(T).unwrap();
        let exact = Cone::new(Frame::WORLD, 1.0, 0.4_f64.atan(), T).unwrap();
        let off = deviation(|p| exact.distance_to(p), &patch);
        assert!(off < 1e-12, "off the cone by {off}");
        assert!(spans_the_same(&cone, &patch));
    }

    #[test]
    fn an_extrusion_becomes_its_profile_lofted() {
        let circle = crate::curve::CircleCurve::new(Circle::new(Frame::WORLD, 3.0, T).unwrap());
        let extrusion: SurfaceGeometry = ExtrusionSurface::new(circle.into(), Direction::Z, 6.0)
            .unwrap()
            .into();
        let patch = extrusion.to_bspline(T).unwrap();
        // A circle swept along its own axis is a cylinder, so the implicit test
        // is the cylinder's.
        let exact = Cylinder::new(Frame::WORLD, 3.0, T).unwrap();
        let off = deviation(|p| exact.distance_to(p), &patch);
        assert!(off < 1e-12, "off the swept circle by {off}");
        assert!(spans_the_same(&extrusion, &patch));
    }

    #[test]
    fn a_patch_converts_to_itself() {
        let cylinder: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 1.0, T).unwrap(), (0.0, 1.0))
                .unwrap()
                .into();
        let patch = cylinder.to_bspline(T).unwrap();
        let again: SurfaceGeometry = patch.clone().into();
        let twice = again.to_bspline(T).unwrap();
        assert_eq!(patch, twice);
    }

    #[test]
    fn a_sphere_becomes_an_exact_rational_patch() {
        let sphere = Sphere::new(Frame::WORLD, 2.5, T).unwrap();
        let surface: SurfaceGeometry = SphereSurface::new(sphere).into();
        let patch = surface.to_bspline(T).unwrap();
        assert!(patch.is_rational(), "a sphere needs weights");
        assert!(deviation(|p| sphere.distance_to(p), &patch) < 1e-12);
        assert!(spans_the_same(&surface, &patch));
    }

    #[test]
    fn a_torus_becomes_an_exact_rational_patch() {
        let torus = Torus::new(Frame::WORLD, 3.0, 1.0, T).unwrap();
        let surface: SurfaceGeometry = TorusSurface::new(torus).into();
        let patch = surface.to_bspline(T).unwrap();
        assert!(patch.is_rational(), "a torus needs weights");
        assert!(deviation(|p| torus.distance_to(p), &patch) < 1e-12);
        assert!(spans_the_same(&surface, &patch));
    }

    #[test]
    fn a_revolution_becomes_the_exact_patch_its_own_construction_is() {
        // A line parallel to the axis, revolved three quarters of a turn: the
        // surface is a cylinder wall, so the patch can be measured against
        // the cylinder's own signed distance, an independent authority, not
        // the revolution evaluating itself.
        use crate::curve::LineCurve;
        let line =
            LineCurve::segment(Point::new(2.0, 0.0, 0.0), Point::new(2.0, 0.0, 5.0), T).unwrap();
        let surface: SurfaceGeometry = crate::surface::RevolutionSurface::new(
            line.into(),
            ogeom_math::Axis {
                location: Point::new(0.0, 0.0, 0.0),
                direction: ogeom_math::Direction::Z,
            },
            1.5 * core::f64::consts::PI,
        )
        .unwrap()
        .into();
        let cylinder = ogeom_math::Cylinder::new(Frame::WORLD, 2.0, T).unwrap();
        let patch = surface.to_bspline(T).unwrap();
        assert!(deviation(|p| cylinder.distance_to(p), &patch) < 1e-12);
        assert!(spans_the_same(&surface, &patch));
    }

    #[test]
    fn a_trimmed_surface_converts_as_its_basis_over_the_window() {
        let plane: SurfaceGeometry = PlaneSurface::new(Plane::new(Frame::WORLD)).into();
        let trimmed: SurfaceGeometry = SurfaceGeometry::Trimmed(Box::new(
            TrimmedSurface::new(plane, (0.0, 1.0), (2.0, 5.0), T).unwrap(),
        ));
        let patch = trimmed.to_bspline(T).unwrap();
        let flat = Plane::new(Frame::WORLD);
        assert!(deviation(|p| flat.distance_to(p), &patch) < 1e-12);
        assert!(spans_the_same(&trimmed, &patch));
    }

    #[test]
    fn a_trimmed_spline_converts_to_its_piece() {
        let plane: SurfaceGeometry =
            PlaneSurface::over(Plane::new(Frame::WORLD), (0.0, 4.0), (0.0, 4.0))
                .unwrap()
                .into();
        let spline: SurfaceGeometry = plane.to_bspline(T).unwrap().into();
        let trimmed: SurfaceGeometry = SurfaceGeometry::Trimmed(Box::new(
            TrimmedSurface::new(spline, (0.25, 0.5), (0.0, 0.75), T).unwrap(),
        ));
        let patch = trimmed.to_bspline(T).unwrap();
        assert!(spans_the_same(&trimmed, &patch));
    }

    #[test]
    fn a_spline_segment_keeps_its_parameters() {
        let cylinder: SurfaceGeometry = CylinderSurface::new(
            ogeom_math::Cylinder::new(Frame::WORLD, 2.0, T).unwrap(),
            (0.0, 3.0),
        )
        .unwrap()
        .into();
        let whole = cylinder.to_bspline(T).unwrap();
        let ((ua, ub), (va, vb)) = whole.domain();
        let u = (ua + 0.2 * (ub - ua), ua + 0.7 * (ub - ua));
        let v = (va + 0.1 * (vb - va), va + 0.6 * (vb - va));
        let piece = whole.segment(u, v, T).unwrap();
        assert_eq!(piece.domain(), (u, v));
        for i in 0..=8 {
            for j in 0..=8 {
                let s = u.0 + (u.1 - u.0) * f64::from(i) / 8.0;
                let t = v.0 + (v.1 - v.0) * f64::from(j) / 8.0;
                let off = piece
                    .point_at(s, t, T)
                    .unwrap()
                    .distance(whole.point_at(s, t, T).unwrap());
                assert!(off < 1e-12, "{off} off at ({s}, {t})");
            }
        }
    }
}

// --- general affine transforms ------------------------------------------------

impl Curve {
    /// This curve carried through a general affine transform.
    ///
    /// A shear or an uneven scale is not a placement: it carries a circle to an
    /// ellipse, and an ellipse to a conic with no analytic name here. So the
    /// curve is converted to its exact B-spline form first and the *control
    /// points* are moved, which an affine map does exactly: a B-spline is an
    /// affine combination of its control points, so transforming them and
    /// transforming every point of the curve are the same thing.
    ///
    /// The result is a B-spline whatever went in, and its parameterization is
    /// the converted one rather than the original. That is the price of a
    /// transform the analytic types cannot express, and it is why
    /// [`transformed`](crate::traits::Transformable::transformed) takes only a
    /// [`Transform`](ogeom_math::Transform): a placement keeps the type, and only
    /// this does not.
    ///
    /// # Errors
    ///
    /// As [`Curve::to_bspline`], plus
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// transform collapses the curve to a point.
    pub fn general_transformed(
        &self,
        t: &ogeom_math::GeneralTransform,
        tol: Tolerances,
    ) -> OgeomResult<BSplineCurve> {
        let spline = self.to_bspline(tol)?;
        let moved: Vec<Weighted<Point>> = spline
            .control_points()
            .iter()
            .map(|c| {
                // Stored homogeneous, so the point has already been multiplied
                // by its weight. An affine map is not linear (it has a
                // translation), so the translation has to be scaled by the
                // weight too, or a rational curve's control points drift apart
                // from its weights and the curve leaves the shape entirely.
                let position = c.point();
                Weighted {
                    scaled: Point::from_vector(t.apply(position).to_vector() * c.weight),
                    weight: c.weight,
                }
            })
            .collect();
        if moved
            .iter()
            .all(|c| c.point().is_equal(moved[0].point(), tol))
        {
            ogeom_bail!(
                Construction,
                "this transform collapses the curve to a point; it is singular \
                 in the curve's own directions"
            );
        }
        BSplineCurve::rational(spline.knots().clone(), moved)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod affine_tests {
    use super::*;
    use crate::curve::{CircleCurve, LineCurve};
    use ogeom_math::{Circle, GeneralTransform, Matrix3, Vector};

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn an_uneven_scale_carries_a_circle_onto_the_ellipse_it_should() {
        // The case a placement cannot express. A circle of radius one scaled by
        // three in x and one in y is the ellipse with those radii, and the
        // converted curve lands on it exactly rather than near it.
        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 1.0, T).unwrap()).into();
        let stretch = GeneralTransform::new(
            Matrix3::new([[3.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            Vector::ZERO,
        );
        let moved = circle.general_transformed(&stretch, T).unwrap();

        let (a, b) = moved.knots().domain();
        for i in 0..=400 {
            #[allow(clippy::cast_precision_loss)]
            let u = a + (b - a) * i as f64 / 400.0;
            let p = moved.point_at(u, T).unwrap();
            let on_ellipse = (p.x / 3.0).powi(2) + p.y.powi(2);
            assert!(
                (on_ellipse - 1.0).abs() < 1e-12,
                "at {u} the point {p:?} is not on the ellipse: {on_ellipse}"
            );
        }
    }

    #[test]
    fn a_shear_is_exact_because_a_spline_is_an_affine_combination() {
        // Transforming the control points and transforming every point of the
        // curve are the same operation, which is what makes this exact rather
        // than fitted.
        let line: Curve = LineCurve::segment(Point::ORIGIN, Point::new(2.0, 3.0, 0.0), T)
            .unwrap()
            .into();
        let shear = GeneralTransform::new(
            Matrix3::new([[1.0, 0.7, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            Vector::new(1.0, 2.0, 3.0),
        );
        let moved = line.general_transformed(&shear, T).unwrap();

        for i in 0..=50 {
            #[allow(clippy::cast_precision_loss)]
            let u = i as f64 / 50.0;
            let before =
                line.point_at(line.domain().0 + (line.domain().1 - line.domain().0) * u, T);
            let after = moved.point_at(u, T).unwrap();
            assert!(shear.apply(before.unwrap()).is_equal(after, T));
        }
    }

    #[test]
    fn a_rational_curves_weights_survive_the_move() {
        // The trap: control points are stored already multiplied by their
        // weight, so a transform with a translation has to scale the
        // translation by the weight too. Getting that wrong leaves a circle's
        // control points and weights inconsistent, and the curve wanders off
        // the shape entirely, most visibly under a pure translation, where
        // nothing should change but the position.
        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 2.0, T).unwrap()).into();
        let shift = GeneralTransform::new(Matrix3::IDENTITY, Vector::new(10.0, -4.0, 6.0));
        let moved = circle.general_transformed(&shift, T).unwrap();

        let centre = Point::new(10.0, -4.0, 6.0);
        let (a, b) = moved.knots().domain();
        for i in 0..=300 {
            #[allow(clippy::cast_precision_loss)]
            let u = a + (b - a) * i as f64 / 300.0;
            let p = moved.point_at(u, T).unwrap();
            assert!(
                (p.distance(centre) - 2.0).abs() < 1e-12,
                "at {u} the radius is {}",
                p.distance(centre)
            );
        }
    }

    #[test]
    fn a_transform_that_collapses_the_curve_is_refused() {
        // Projecting a circle in the xy plane onto the z axis leaves a point.
        // A "curve" that is one point is not a curve, and returning it would
        // hand back something every later algorithm divides by the length of.
        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 1.0, T).unwrap()).into();
        let squash = GeneralTransform::new(Matrix3::new([[0.0; 3]; 3]), Vector::ZERO);
        assert!(circle.general_transformed(&squash, T).is_err());
    }
}
