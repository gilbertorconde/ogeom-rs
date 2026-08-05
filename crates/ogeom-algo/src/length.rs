//! Arc length, and distributing points along a curve by it.
//!
//! A curve's parameter is not its length. A B-spline traverses its own
//! parameter at whatever speed its knots imply; a cone's slant runs at a rate
//! set by its half angle. So "a point every two millimetres" and "twenty evenly
//! spaced points" are questions about *length*, and answering them from
//! parameter values is answering a different question.
//!
//! # Two operations, one of them an inversion
//!
//! Length is an integral: the speed `|c'(u)|` integrated over the range. That
//! is [`curve_length`], and it is exact to a stated tolerance rather than
//! summed from a polyline.
//!
//! Placing a point *at* a length is the inverse of that integral, which has no
//! closed form for anything but a line and a circle. It is solved rather than
//! approximated — the length from the start is strictly increasing wherever the
//! parameterization is regular, so a bracketed root find always converges, and
//! there is no risk of the multiple-root trouble a general solve would have.
//!
//! # This is not tessellation
//!
//! [`ogeom_mesh::discretize()`] places points where the curve *bends*,
//! which is what a mesh wants and what a drawing wants. These place points
//! where the caller asked, evenly along the curve, which is what a toolpath, a
//! dimension chain or a sampling pattern wants. Neither substitutes for the
//! other: an evenly spaced polyline through a tight corner misses it, and a
//! deflection-driven one has no even spacing to speak of.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d};
use ogeom_math::{Point, integrate, solve};

/// The arc length of a curve over `range`.
///
/// Integrated from the curve's own speed, so it is the length of the *curve*
/// rather than of a polyline that approximates it. A tessellated length is
/// always short — every chord cuts a corner — and the shortfall is exactly what
/// a deflection tolerance permits, which is far larger than this.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if `range` is not finite;
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the integral does not
/// converge, which means the parameterization is singular somewhere in the
/// range rather than that the curve is long.
pub fn curve_length(curve: &Curve, range: (f64, f64), tol: Tolerances) -> OgeomResult<f64> {
    let (lo, hi) = range;
    if !lo.is_finite() || !hi.is_finite() {
        ogeom_bail!(Domain, "cannot measure the length of [{lo}, {hi}]");
    }
    if (hi - lo).abs() <= tol.parametric() {
        return Ok(0.0);
    }
    // The magnitude of the derivative is the speed along the curve, and its
    // integral is the distance covered. A failure to evaluate is a zero
    // contribution rather than a panic: the integrator samples inside the
    // range, and a curve that cannot be differentiated there has a singular
    // parameterization, which is what the integrator will then report.
    let speed = |u: f64| curve.d1_at(u, tol).map_or(0.0, |d| d.magnitude());
    let length = integrate(speed, lo, hi, tol.confusion())?;
    Ok(length.abs())
}

/// The parameter at which a given arc length from the start of `range` is
/// reached.
///
/// The inverse of [`curve_length`]. `target` is measured from `range.0`, and
/// must lie between zero and the curve's total length over the range — asking
/// for a point beyond the end is refused rather than clamped, because a
/// clamped answer is indistinguishable from a correct one at the end.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if `target` is negative or
/// past the end; [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the solve
/// does not converge.
pub fn parameter_at_length(
    curve: &Curve,
    range: (f64, f64),
    target: f64,
    tol: Tolerances,
) -> OgeomResult<f64> {
    let total = curve_length(curve, range, tol)?;
    if !target.is_finite() || target < -tol.confusion() {
        ogeom_bail!(
            Domain,
            "arc length {target} is not a distance along a curve"
        );
    }
    if target > total + tol.confusion() {
        ogeom_bail!(
            Domain,
            "asked for the point {target} along a curve {total} long; clamping \
             it would give an answer indistinguishable from a correct one at \
             the end"
        );
    }
    if target <= tol.confusion() {
        return Ok(range.0);
    }
    if target >= total - tol.confusion() {
        return Ok(range.1);
    }

    // Length from the start is strictly increasing wherever the speed is
    // non-zero, so this has exactly one root in the range and a bracketed
    // method cannot land on the wrong one.
    let residual = |u: f64| curve_length(curve, (range.0, u), tol).unwrap_or(0.0) - target;
    let criteria = solve::Criteria {
        // The residual is a *length*, so it is measured against a spatial
        // tolerance; the step is a parameter and is measured against a
        // parametric one.
        residual: tol.confusion(),
        step: tol.parametric(),
        ..solve::Criteria::default()
    };
    Ok(solve::brent(residual, range.0, range.1, criteria)?.value)
}

/// `count` points evenly spaced *by arc length* along a curve, ends included.
///
/// `count` is the number of points, so two gives the ends and nothing between.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if fewer than two
/// points are asked for; otherwise as [`parameter_at_length`].
pub fn points_by_count(
    curve: &Curve,
    range: (f64, f64),
    count: usize,
    tol: Tolerances,
) -> OgeomResult<Vec<(f64, Point)>> {
    if count < 2 {
        ogeom_bail!(
            Construction,
            "a distribution along a curve needs at least its two ends, got \
             {count}"
        );
    }
    let total = curve_length(curve, range, tol)?;
    #[allow(clippy::cast_precision_loss)]
    let step = total / (count - 1) as f64;
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        #[allow(clippy::cast_precision_loss)]
        let at = parameter_at_length(curve, range, step * i as f64, tol)?;
        out.push((at, curve.point_at(at, tol)?));
    }
    Ok(out)
}

/// Points along a curve at a fixed arc-length `spacing`.
///
/// The first point is at the start. The last is at the end *whatever the
/// spacing divides to*, so the final gap is short rather than the curve being
/// left unfinished — a distribution that stops before the end is a different
/// answer from the one asked for, and silently dropping the tail is how a
/// toolpath ends up not reaching the edge of the material.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `spacing` is
/// not finite and positive; otherwise as [`parameter_at_length`].
pub fn points_by_spacing(
    curve: &Curve,
    range: (f64, f64),
    spacing: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<(f64, Point)>> {
    if !spacing.is_finite() || spacing <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "spacing {spacing} must be finite and positive"
        );
    }
    let total = curve_length(curve, range, tol)?;
    let mut out = Vec::new();
    let mut at_length = 0.0;
    while at_length < total - tol.confusion() {
        let at = parameter_at_length(curve, range, at_length, tol)?;
        out.push((at, curve.point_at(at, tol)?));
        at_length += spacing;
    }
    out.push((range.1, curve.point_at(range.1, tol)?));
    Ok(out)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use core::f64::consts::{PI, TAU};
    use ogeom_geom::{BSplineCurve, CircleCurve, LineCurve};
    use ogeom_math::{Circle, Frame, KnotVector};

    const T: Tolerances = Tolerances::millimetres();

    fn circle(radius: f64) -> Curve {
        CircleCurve::new(Circle::new(Frame::WORLD, radius, T).unwrap()).into()
    }

    #[test]
    fn a_lines_length_is_the_distance_between_its_ends() {
        let line: Curve = LineCurve::segment(Point::ORIGIN, Point::new(3.0, 4.0, 0.0), T)
            .unwrap()
            .into();
        assert_relative_eq!(
            curve_length(&line, (0.0, 5.0), T).unwrap(),
            5.0,
            epsilon = 1e-12
        );
        // Its parameter is already arc length, so half the length is half way.
        assert_relative_eq!(
            parameter_at_length(&line, (0.0, 5.0), 2.5, T).unwrap(),
            2.5,
            epsilon = 1e-9
        );
    }

    #[test]
    fn a_circles_length_is_its_circumference_and_an_arcs_is_the_fraction() {
        let c = circle(2.0);
        assert_relative_eq!(
            curve_length(&c, (0.0, TAU), T).unwrap(),
            TAU * 2.0,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            curve_length(&c, (0.0, PI), T).unwrap(),
            PI * 2.0,
            epsilon = 1e-9
        );
    }

    #[test]
    fn length_is_measured_on_the_curve_not_on_a_polyline_through_it() {
        // The distinction that makes this worth having. A chord always cuts a
        // corner, so a tessellated length is short — and short by whatever the
        // deflection permits, which is far more than this integral's error.
        let c = circle(1.0);
        let exact = TAU;
        let integrated = curve_length(&c, (0.0, TAU), T).unwrap();
        assert!(
            (integrated - exact).abs() < 1e-9,
            "got {integrated} against {exact}"
        );

        let mesh =
            ogeom_mesh::discretize(&c, (0.0, TAU), ogeom_mesh::Deflection::default(), T).unwrap();
        assert!(
            mesh.length() < exact - 1e-4,
            "a coarse polyline should be visibly short, got {}",
            mesh.length()
        );
    }

    #[test]
    fn points_by_count_are_evenly_spaced_along_the_curve() {
        // On a circle, even in length is even in angle — which is the check
        // that the inversion is doing its job rather than returning parameters.
        let c = circle(3.0);
        let points = points_by_count(&c, (0.0, TAU), 9, T).unwrap();
        assert_eq!(points.len(), 9);

        let step = TAU / 8.0;
        for (i, (at, _)) in points.iter().enumerate() {
            #[allow(clippy::cast_precision_loss)]
            let want = step * i as f64;
            assert!((at - want).abs() < 1e-6, "point {i} at {at}, wanted {want}");
        }
        // And the chords between consecutive points are all the same length.
        let first = points[0].1.distance(points[1].1);
        for pair in points.windows(2) {
            assert_relative_eq!(pair[0].1.distance(pair[1].1), first, max_relative = 1e-6);
        }
    }

    #[test]
    fn an_unevenly_parameterized_curve_is_still_evenly_divided() {
        // The case a parameter-space distribution gets wrong. This spline's
        // knots make it cover ground at very different speeds, so equal
        // parameter steps are not equal distances and equal distances are not
        // equal parameter steps.
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 0.2, 1.0, 1.0, 1.0, 1.0], 3).unwrap();
        let control = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(0.5, 4.0, 0.0),
            Point::new(3.0, 4.0, 0.0),
            Point::new(9.0, 0.5, 0.0),
            Point::new(10.0, 0.0, 0.0),
        ];
        let spline: Curve = BSplineCurve::new(knots, control, T).unwrap().into();
        let range = spline.domain();

        let points = points_by_count(&spline, range, 12, T).unwrap();
        // Equal *arc* length between consecutive points, which is what was
        // asked for. Not equal chords: a chord cuts the corner, so where this
        // curve bends hardest its chord is several percent shorter than the arc
        // it spans, and asserting on chords would be asserting the curve is
        // straight.
        let step = curve_length(&spline, (points[0].0, points[1].0), T).unwrap();
        for pair in points.windows(2) {
            let along = curve_length(&spline, (pair[0].0, pair[1].0), T).unwrap();
            assert_relative_eq!(along, step, max_relative = 1e-6);
        }

        // And the parameters are *not* evenly spaced, which is the whole point.
        let steps: Vec<f64> = points.windows(2).map(|w| w[1].0 - w[0].0).collect();
        let spread = steps.iter().fold(0.0_f64, |a, b| a.max(*b))
            / steps.iter().fold(f64::MAX, |a, b| a.min(*b));
        assert!(
            spread > 1.5,
            "this curve's parameter should be visibly uneven, spread {spread}"
        );
    }

    #[test]
    fn spacing_always_reaches_the_end_even_when_it_does_not_divide() {
        let c = circle(1.0);
        let total = TAU;
        // Deliberately does not divide the circumference.
        let points = points_by_spacing(&c, (0.0, total), 1.0, T).unwrap();
        assert!(points.len() >= 7);
        assert_relative_eq!(points[0].0, 0.0, epsilon = 1e-12);
        assert_relative_eq!(points[points.len() - 1].0, total, epsilon = 1e-9);

        // Every gap but the last is the spacing asked for; the last is short.
        for pair in points[..points.len() - 1].windows(2) {
            let along = curve_length(&c, (pair[0].0, pair[1].0), T).unwrap();
            assert_relative_eq!(along, 1.0, max_relative = 1e-6);
        }
        let tail = curve_length(
            &c,
            (points[points.len() - 2].0, points[points.len() - 1].0),
            T,
        )
        .unwrap();
        assert!(
            tail <= 1.0 + 1e-9,
            "the last gap should be short, got {tail}"
        );
    }

    #[test]
    fn asking_beyond_the_end_is_refused_rather_than_clamped() {
        // A clamped answer sits exactly where a correct one at the end would,
        // so the caller cannot tell the difference.
        let c = circle(1.0);
        assert!(parameter_at_length(&c, (0.0, PI), PI * 2.0, T).is_err());
        assert!(parameter_at_length(&c, (0.0, PI), -1.0, T).is_err());
        assert!(parameter_at_length(&c, (0.0, PI), PI, T).is_ok());
    }

    #[test]
    fn distributions_that_describe_nothing_are_refused() {
        let c = circle(1.0);
        assert!(points_by_count(&c, (0.0, TAU), 1, T).is_err());
        assert!(points_by_count(&c, (0.0, TAU), 0, T).is_err());
        assert!(points_by_spacing(&c, (0.0, TAU), 0.0, T).is_err());
        assert!(points_by_spacing(&c, (0.0, TAU), -1.0, T).is_err());
        assert!(points_by_spacing(&c, (0.0, TAU), f64::NAN, T).is_err());
    }
}
