//! Curves fitted through a trace at its own parameters and held to their
//! tolerance between the samples.
//!
//! A fit's own error is measured at the points it was fitted through. A
//! trace with more detail than the samples (a many-span spline, a polyline
//! with more corners than samples) leaves a fit that meets its target at
//! every sample and misses the trace between them. These fits sample a
//! source curve at its own breaks, check each interval between samples
//! at its eighths, split the intervals that miss, and report the worst
//! distance measured anywhere, not the worst at the samples. The fits are
//! [`ogeom_geom::fit::fit_curve_sampled`] and its planar twin; this module
//! chooses where they start.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::fit::Fitted;
use ogeom_geom::{BSpline2d, BSplineCurve, Curve, Curve2d as _, Curve3d as _, PlanarCurve};
use ogeom_math::{KnotVector, Point, Point2};

/// The most intervals a trace is sampled at, as many as a sampled fit
/// refines to.
const MOST: usize = ogeom_geom::fit::SAMPLED_SPANS;

/// The parameters to start a fit of `curve` over `range` from: each span of
/// the curve (between its distinct knots) cut in `degree + 1` pieces, and
/// at least `at_least` pieces overall.
#[must_use]
pub fn stations(curve: &Curve, range: (f64, f64), at_least: usize) -> Vec<f64> {
    let (breaks, degree) = breaks_3d(curve);
    spread(range, &breaks, degree, at_least)
}

/// As [`stations`], for a planar curve.
#[must_use]
pub fn stations_2d(curve: &PlanarCurve, range: (f64, f64), at_least: usize) -> Vec<f64> {
    let (breaks, degree) = breaks_2d(curve);
    spread(range, &breaks, degree, at_least)
}

/// As [`stations`], for a range cut at `breaks` into pieces of `degree`:
/// a patch's direction between its distinct knots, say.
#[must_use]
pub fn stations_between(
    range: (f64, f64),
    breaks: &[f64],
    degree: usize,
    at_least: usize,
) -> Vec<f64> {
    spread(range, breaks, degree, at_least)
}

/// Where a space curve's pieces join, and the degree of its pieces.
fn breaks_3d(curve: &Curve) -> (Vec<f64>, usize) {
    match curve {
        Curve::BSpline(spline) => (
            spline
                .knots()
                .distinct()
                .into_iter()
                .map(|(value, _)| value)
                .collect(),
            spline.knots().degree(),
        ),
        Curve::Trimmed(trim) => {
            let (lo, hi) = trim.domain();
            let (inner, degree) = breaks_3d(trim.basis());
            if trim.is_reversed() {
                (inner.into_iter().map(|b| lo + hi - b).collect(), degree)
            } else {
                (inner, degree)
            }
        }
        Curve::Offset(offset) => breaks_3d(offset.basis()),
        Curve::OnSurface(on) => breaks_2d(on.pcurve()),
        _ => (Vec::new(), 3),
    }
}

/// Where a planar curve's pieces join, and the degree of its pieces.
fn breaks_2d(curve: &PlanarCurve) -> (Vec<f64>, usize) {
    match curve {
        PlanarCurve::BSpline(spline) => (
            spline
                .knots()
                .distinct()
                .into_iter()
                .map(|(value, _)| value)
                .collect(),
            spline.knots().degree(),
        ),
        PlanarCurve::Trimmed(trim) => {
            let (lo, hi) = trim.domain();
            let (inner, degree) = breaks_2d(trim.basis());
            if trim.is_reversed() {
                (inner.into_iter().map(|b| lo + hi - b).collect(), degree)
            } else {
                (inner, degree)
            }
        }
        PlanarCurve::Offset(offset) => breaks_2d(offset.basis()),
        _ => (Vec::new(), 3),
    }
}

/// `range` cut at the `breaks` inside it, each piece in `degree + 1`, and
/// evened out to at least `at_least` pieces; never more than [`MOST`].
fn spread(range: (f64, f64), breaks: &[f64], degree: usize, at_least: usize) -> Vec<f64> {
    let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
    let uniform = |n: usize| -> Vec<f64> {
        (0..=n)
            .map(|i| {
                #[allow(clippy::cast_precision_loss, reason = "a sample index")]
                let f = i as f64 / n as f64;
                if i == n { hi } else { lo + (hi - lo) * f }
            })
            .collect()
    };
    let at_least = at_least.clamp(1, MOST);
    if hi - lo <= 0.0 {
        return vec![lo, hi];
    }
    let mut cuts = vec![lo];
    cuts.extend(
        breaks
            .iter()
            .copied()
            .filter(|&b| b > lo + (hi - lo) * 1e-9 && b < hi - (hi - lo) * 1e-9),
    );
    cuts.push(hi);
    let pieces = (degree + 1).max(2);
    if (cuts.len() - 1) * pieces > MOST {
        return uniform(MOST);
    }
    let mut out = Vec::with_capacity((cuts.len() - 1) * pieces + 1);
    for pair in cuts.windows(2) {
        for k in 0..pieces {
            #[allow(clippy::cast_precision_loss, reason = "a small count")]
            let f = k as f64 / pieces as f64;
            out.push(pair[0] + (pair[1] - pair[0]) * f);
        }
    }
    out.push(hi);
    if out.len() > at_least {
        return out;
    }
    // Too few: every piece cut evenly until there are enough.
    let per = (at_least).div_ceil(out.len() - 1);
    let mut finer = Vec::with_capacity((out.len() - 1) * per + 1);
    for pair in out.windows(2) {
        for k in 0..per {
            #[allow(clippy::cast_precision_loss, reason = "a small count")]
            let f = k as f64 / per as f64;
            finer.push(pair[0] + (pair[1] - pair[0]) * f);
        }
    }
    finer.push(hi);
    finer
}

/// Where a planar spline over `range` turns a corner: its interior knots
/// of multiplicity at least its degree, inside `range`.
#[must_use]
pub fn corners_2d(curve: &BSpline2d, range: (f64, f64)) -> Vec<f64> {
    let knots = curve.knots();
    let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
    knots
        .distinct()
        .into_iter()
        .filter(|&(value, mult)| mult >= knots.degree() && value > lo && value < hi)
        .map(|(value, _)| value)
        .collect()
}

/// Polynomial splines of one degree, each clamped and each starting where
/// the one before ends in parameter, joined end to end into one spline
/// that is continuous at the joins and no smoother. Each join takes the
/// middle of the two ends that meet there; the second value is how far
/// that moved either end.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// there are no pieces, they differ in degree, one is rational, or their
/// parameter ranges do not abut.
pub fn joined(pieces: &[BSplineCurve], tol: Tolerances) -> OgeomResult<(BSplineCurve, f64)> {
    let Some(first) = pieces.first() else {
        ogeom_bail!(Construction, "nothing to join");
    };
    let degree = first.degree();
    let mut knots: Vec<f64> = Vec::new();
    let mut control: Vec<Point> = Vec::new();
    let mut moved = 0.0_f64;
    for (i, piece) in pieces.iter().enumerate() {
        if piece.degree() != degree || piece.is_rational() {
            ogeom_bail!(Construction, "only polynomial pieces of one degree join");
        }
        let own = piece.knots().knots();
        let points: Vec<Point> = piece.control_points().iter().map(|c| c.point()).collect();
        if own.len() < 2 * (degree + 1) || points.is_empty() {
            ogeom_bail!(Construction, "a piece too short to join");
        }
        if i == 0 {
            knots.extend_from_slice(&own[..own.len() - 1]);
            control.extend(points);
            continue;
        }
        let at = own[0];
        let held = knots[knots.len() - 1];
        if (at - held).abs() > tol.parametric().max(held.abs() * 1e-12) {
            ogeom_bail!(Construction, "pieces to join must abut in parameter");
        }
        let last = control.len() - 1;
        let gap = control[last].distance(points[0]);
        moved = moved.max(gap * 0.5);
        control[last] = control[last].midpoint(points[0]);
        // The join keeps `degree` copies of its knot; the piece's first
        // knot copies past the clamp's are dropped with its first point.
        knots.extend_from_slice(&own[degree + 1..own.len() - 1]);
        control.extend(points.into_iter().skip(1));
    }
    let end = pieces[pieces.len() - 1].knots().knots();
    knots.push(end[end.len() - 1]);
    let curve = BSplineCurve::new(KnotVector::new(knots, degree)?, control, tol)?;
    Ok((curve, moved))
}

/// A space curve fitted to `trace` at its own parameters, starting from the
/// samples at `stations` (strictly increasing): the open case of
/// [`ogeom_geom::fit::fit_curve_sampled`]. Each interval is measured at its
/// eighths as well as at the samples, and the intervals that miss
/// `tolerance` are split until it holds everywhere measured or the samples
/// reach [`ogeom_geom::fit::SAMPLED_SPANS`] intervals. The best fit either
/// way, its `error` the worst distance measured and `met` whether that is
/// within `tolerance`.
///
/// # Errors
///
/// As `trace`, and as [`ogeom_geom::fit::fit_curve_sampled`].
pub fn fit_traced(
    trace: impl FnMut(f64) -> OgeomResult<Point>,
    stations: &[f64],
    degree: usize,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Fitted<BSplineCurve>> {
    ogeom_geom::fit::fit_curve_sampled(trace, stations, false, degree, tolerance, tol)
}

/// As [`fit_traced`], in the plane: a chart image fitted at the parameters
/// of the curve it annotates. The distances are in the trace's own units.
///
/// # Errors
///
/// As `trace`, and as [`ogeom_geom::fit::fit_curve_2d_sampled`].
pub fn fit_traced_2d(
    trace: impl FnMut(f64) -> OgeomResult<Point2>,
    stations: &[f64],
    degree: usize,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Fitted<BSpline2d>> {
    ogeom_geom::fit::fit_curve_2d_sampled(trace, stations, false, degree, tolerance, tol)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;

    const T: Tolerances = Tolerances::millimetres();

    /// A trace weaving faster than the starting samples is fitted to its
    /// tolerance everywhere, and a fit that cannot reach it says by how
    /// much at the worst point measured, not at its samples.
    #[test]
    fn a_traced_fit_holds_between_its_samples() {
        let trace =
            |t: f64| -> OgeomResult<Point> { Ok(Point::new(t, 0.05 * (t * 173.0).sin(), 0.0)) };
        let start: Vec<f64> = (0..=8).map(|i| f64::from(i) / 8.0).collect();
        let fitted = fit_traced(trace, &start, 3, 1e-4, T).unwrap();
        assert!(fitted.met, "{}", fitted.error);
        let mut worst = 0.0_f64;
        for i in 0..=5000 {
            let t = f64::from(i) / 5000.0;
            let p = fitted.curve.point_at(t, T).unwrap();
            worst = worst.max(p.distance(trace(t).unwrap()));
        }
        assert!(worst <= 1e-4, "{worst}");
    }

    /// A many-span spline is sampled at every span, not at a fixed count.
    #[test]
    fn stations_follow_the_spans() {
        let control: Vec<Point> = (0..67)
            .map(|i| Point::new(f64::from(i), f64::from(i % 2), 0.0))
            .collect();
        let knots = KnotVector::clamped_uniform(3, control.len()).unwrap();
        let curve = Curve::BSpline(BSplineCurve::new(knots, control, T).unwrap());
        let s = stations(&curve, (0.0, 1.0), 33);
        assert_eq!(s.len(), 64 * 4 + 1);
        assert!(s.windows(2).all(|w| w[1] > w[0]));
        let few = stations(&curve, (0.25, 0.5), 33);
        assert!(few.len() >= 34 && few[0] == 0.25 && few[few.len() - 1] == 0.5);
    }
}
