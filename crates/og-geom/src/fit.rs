//! Fitting a B-spline to points, to a stated error target.
//!
//! The missing half of fitting. Interpolation and fixed-count approximation
//! exist in `og-algo`; what they cannot do is *choose* — the caller names a
//! control-point count and hopes. This module is the loop that closes that:
//! fit, measure where the fit is worst, refine the knots exactly there, and
//! repeat until the error target is met.
//!
//! `docs/SCOPE.md` deferred this with a warning worth repeating: a fit that
//! silently picks its own resolution and reports success is the shape of answer
//! that gets trusted. So the result here carries the error actually reached and
//! whether the target was met, and a fit that ran out of room says so rather
//! than rounding "close" up to "done".
//!
//! # Where the knots go
//!
//! Refinement is *where the error is*, not everywhere. Splitting every span
//! doubles the control points per round and most of them buy nothing — a curve
//! that is straight for most of its length and tight in one corner needs its
//! knots in the corner. Each round measures the error per span and splits only
//! the spans that exceed the target, so the knot density ends up tracking the
//! curvature, which is where it belongs.
//!
//! # Who this is for
//!
//! The marching intersector, first: a traced branch is a polyline with a stated
//! chord tolerance, and downstream code wants a curve, so the polyline is fitted
//! to the same tolerance and the result is as good as the trace. But nothing
//! here knows about intersections — it fits points, in three dimensions or two,
//! which is also what a digitized profile or an imported polyline needs.

use og_core::{OgResult, Tolerances, og_bail};
use og_math::{KnotVector, Point, Point2};

use crate::curve::BSplineCurve;
use crate::curve2d::BSpline2d;

/// What a fit produced.
#[derive(Debug, Clone, PartialEq)]
pub struct Fitted<C> {
    /// The curve.
    pub curve: C,
    /// The largest distance from any input point to the curve at its
    /// parameter.
    pub error: f64,
    /// Whether the error target was met.
    ///
    /// A fit can run out of room — more control points than points to fit
    /// solves nothing, since at that ratio least squares *is* interpolation —
    /// and then this is `false` and `error` says how close it got. Reported
    /// rather than rounded up to success.
    pub met: bool,
}

/// Fit a spline through 3D points, refining until `tolerance` is met.
///
/// The first and last points are honoured exactly: they are where the curve
/// joins whatever comes next, and a fit that drifts at its ends produces gaps
/// at every junction built on it. A closed polyline — first point repeated at
/// the end — therefore comes back closed.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if there are
/// fewer than two distinct points or the tolerance is not a positive distance.
pub fn fit_points(
    points: &[Point],
    degree: usize,
    tolerance: f64,
    tol: Tolerances,
) -> OgResult<Fitted<BSplineCurve>> {
    let (knots, control, error, met) = fit::<3>(
        &points.iter().map(|p| [p.x, p.y, p.z]).collect::<Vec<_>>(),
        degree,
        tolerance,
        tol,
    )?;
    let curve = BSplineCurve::new(
        knots,
        control
            .into_iter()
            .map(|c| Point::new(c[0], c[1], c[2]))
            .collect(),
        tol,
    )?;
    Ok(Fitted { curve, error, met })
}

/// Fit a spline through 2D points, refining until `tolerance` is met.
///
/// The planar twin of [`fit_points`], for pcurves: an intersection curve lives
/// on both surfaces, and each face needs it in its own parameter space.
///
/// # Errors
///
/// As [`fit_points`].
pub fn fit_points_2d(
    points: &[Point2],
    degree: usize,
    tolerance: f64,
    tol: Tolerances,
) -> OgResult<Fitted<BSpline2d>> {
    let (knots, control, error, met) = fit::<2>(
        &points.iter().map(|p| [p.x, p.y]).collect::<Vec<_>>(),
        degree,
        tolerance,
        tol,
    )?;
    let curve = BSpline2d::new(
        knots,
        control
            .into_iter()
            .map(|c| Point2::new(c[0], c[1]))
            .collect(),
        tol,
    )?;
    Ok(Fitted { curve, error, met })
}

/// Fit one curve living in three spaces at once: a 3D curve and its two
/// parameter-space images, as a single seven-dimensional fit.
///
/// One parameterization, one knot vector, one correction: the three results
/// are same-parameter *by construction*, which separate fits cannot promise —
/// each fit's parameter correction drifts its parameterization independently,
/// and the drift is invisible to every per-fit error measure. The boolean
/// found that: pcurves claiming 1e-7 evaluated millimetres from their own
/// curve. The reported error bounds the worst deviation across all seven
/// coordinates, so it bounds each space's deviation too.
///
/// # Errors
///
/// As [`fit_points`], and the three inputs must be equally long.
#[allow(clippy::type_complexity)]
pub fn fit_points_joint(
    points: &[Point],
    on_a: &[Point2],
    on_b: &[Point2],
    degree: usize,
    tolerance: f64,
    tol: Tolerances,
) -> OgResult<(Fitted<BSplineCurve>, BSpline2d, BSpline2d)> {
    if points.len() != on_a.len() || points.len() != on_b.len() {
        og_bail!(
            Construction,
            "a joint fit needs the same trace seen in every space"
        );
    }
    let joined: Vec<[f64; 7]> = points
        .iter()
        .zip(on_a)
        .zip(on_b)
        .map(|((p, a), b)| [p.x, p.y, p.z, a.x, a.y, b.x, b.y])
        .collect();
    let (knots, control, error, met) = fit::<7>(&joined, degree, tolerance, tol)?;
    let curve = BSplineCurve::new(
        knots.clone(),
        control
            .iter()
            .map(|c| Point::new(c[0], c[1], c[2]))
            .collect(),
        tol,
    )?;
    let pa = BSpline2d::new(
        knots.clone(),
        control.iter().map(|c| Point2::new(c[3], c[4])).collect(),
        tol,
    )?;
    let pb = BSpline2d::new(
        knots,
        control.iter().map(|c| Point2::new(c[5], c[6])).collect(),
        tol,
    )?;
    Ok((Fitted { curve, error, met }, pa, pb))
}

/// Fit a pcurve at *fixed* parameters — the source curve's own.
///
/// The fixed-parameter twin of [`fit_points_2d`], and the difference is the
/// contract: parameter correction is what makes a free fit's residual honest,
/// and it is exactly what a *same-parameter* fit must never do, because the
/// parameters are not a guess to be improved — they are the 3D curve's own,
/// and drifting them is how a pcurve ends up evaluating away from the curve
/// it annotates. Here the parameters stay put, refinement adds knots where
/// the error says, and the reported error is the true same-parameter
/// deviation in the chart.
///
/// # Errors
///
/// As [`fit_points_2d`], and the parameters must be strictly increasing and
/// as many as the points.
pub fn fit_points_2d_at(
    parameters: &[f64],
    points: &[Point2],
    degree: usize,
    tolerance: f64,
    tol: Tolerances,
) -> OgResult<Fitted<BSpline2d>> {
    if parameters.len() != points.len() {
        og_bail!(Construction, "one parameter per point, or the fit is a lie");
    }
    if parameters.windows(2).any(|w| w[1] <= w[0]) {
        og_bail!(Construction, "fixed parameters must strictly increase");
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        og_bail!(Construction, "a tolerance of {tolerance} is not a distance");
    }
    if degree == 0 {
        og_bail!(Construction, "a fit needs a degree of at least one");
    }
    let data: Vec<[f64; 2]> = points.iter().map(|p| [p.x, p.y]).collect();
    if data.len() < 2 {
        og_bail!(Construction, "a fit needs at least two points");
    }
    let degree = degree.min(data.len() - 1);
    let (a, b) = (parameters[0], parameters[parameters.len() - 1]);
    let mut knots = single_span(degree, a, b)?;

    const ROUNDS: usize = 32;
    let mut best: Option<(KnotVector, Vec<[f64; 2]>, f64)> = None;
    for _ in 0..ROUNDS {
        let control = least_squares::<2>(&knots, &data, parameters)?;
        let errors = residuals::<2>(&knots, &control, &data, parameters);
        let worst = errors.iter().fold(0.0_f64, |acc, e| acc.max(e.1));
        if best.as_ref().is_none_or(|(_, _, held)| worst < *held) {
            best = Some((knots.clone(), control.clone(), worst));
        }
        if worst <= tolerance {
            let curve = BSpline2d::new(
                knots,
                control
                    .into_iter()
                    .map(|c| Point2::new(c[0], c[1]))
                    .collect(),
                tol,
            )?;
            return Ok(Fitted {
                curve,
                error: worst,
                met: true,
            });
        }
        if knots.control_point_count() >= data.len() {
            break;
        }
        let Some(refined) = refined_where_bad(&knots, &errors, tolerance)? else {
            break;
        };
        knots = refined;
    }
    let (knots, control, error) =
        best.ok_or_else(|| og_core::og_err!(Construction, "the fixed-parameter fit never solved"))?;
    let curve = BSpline2d::new(
        knots,
        control
            .into_iter()
            .map(|c| Point2::new(c[0], c[1]))
            .collect(),
        tol,
    )?;
    Ok(Fitted {
        curve,
        error,
        met: false,
    })
}

/// The dimension-generic core.
///
/// Least squares is solved coordinate by coordinate: the collocation matrix
/// depends only on the parameters and the knots, so the expensive part is
/// shared and each coordinate is one more right-hand side.
#[allow(clippy::type_complexity)]
fn fit<const D: usize>(
    points: &[[f64; D]],
    degree: usize,
    tolerance: f64,
    tol: Tolerances,
) -> OgResult<(KnotVector, Vec<[f64; D]>, f64, bool)> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        og_bail!(Construction, "a tolerance of {tolerance} is not a distance");
    }
    if degree == 0 {
        og_bail!(Construction, "a fit needs a degree of at least one");
    }
    let points = collapse::<D>(points, tol);
    if points.len() < 2 {
        og_bail!(
            Construction,
            "a fit needs at least two distinct points, got {}",
            points.len()
        );
    }
    let degree = degree.min(points.len() - 1);
    let parameters = centripetal::<D>(&points);

    // Start with the fewest control points a clamped curve of this degree can
    // have: one Bézier span. Refinement adds knots only where the error says.
    let (a, b) = (parameters[0], parameters[parameters.len() - 1]);
    let mut knots = single_span(degree, a, b)?;

    // Each round may split every offending span, so the count grows by at most
    // a factor of two a round; a cap keeps a pathological input from spinning.
    const ROUNDS: usize = 32;
    let mut parameters = parameters;
    let mut best: Option<(KnotVector, Vec<[f64; D]>, f64)> = None;
    for _ in 0..ROUNDS {
        // Parameter correction in an earlier round can slide the data out of
        // a span that refinement checked *before* the slide, and the solve
        // reports itself singular. The rounds before it stand: keep the best
        // of them rather than promoting a bookkeeping casualty to an error.
        let control = match least_squares::<D>(&knots, &points, &parameters) {
            Ok(control) => control,
            Err(e) => {
                if let Some((knots, control, worst)) = best {
                    return Ok((knots, control, worst, false));
                }
                return Err(e);
            }
        };
        // Parameter correction, and it is not a refinement — it is what makes
        // the residual mean anything. The residual is measured at each point's
        // assigned parameter, and the centripetal assignment is a guess: where
        // it drifts from the curve's own flow, a *perfect* curve still shows an
        // error at the assigned spot, the loop reads that as the curve's
        // fault, and refinement adds knots forever against a floor it can
        // never get under. Projecting each point onto the current curve —
        // Newton on the foot of the perpendicular — removes the
        // parameterization's share of the error and leaves the curve's.
        for _ in 0..2 {
            correct_parameters::<D>(&knots, &control, &points, &mut parameters);
        }
        let errors = residuals::<D>(&knots, &control, &points, &parameters);
        let worst = errors.iter().fold(0.0_f64, |acc, e| acc.max(e.1));
        if best.as_ref().is_none_or(|(_, _, held)| worst < *held) {
            best = Some((knots.clone(), control.clone(), worst));
        }
        if worst <= tolerance {
            return Ok((knots, control, worst, true));
        }

        // More control points than data points is interpolation wearing a
        // different name, and past that adding knots buys nothing.
        if knots.control_point_count() >= points.len() {
            break;
        }
        let Some(refined) = refined_where_bad(&knots, &errors, tolerance)? else {
            break;
        };
        knots = refined;
    }

    #[allow(clippy::unwrap_used, reason = "at least one round always runs")]
    let (knots, control, worst) = best.unwrap();
    Ok((knots, control, worst, false))
}

/// Drop consecutive duplicates, which contribute a zero-length chord and make
/// the parameterization stall.
fn collapse<const D: usize>(points: &[[f64; D]], tol: Tolerances) -> Vec<[f64; D]> {
    let mut out: Vec<[f64; D]> = Vec::with_capacity(points.len());
    for p in points {
        if out
            .last()
            .is_some_and(|q| distance::<D>(p, q) <= tol.confusion() * 0.01)
        {
            continue;
        }
        out.push(*p);
    }
    out
}

fn distance<const D: usize>(a: &[f64; D], b: &[f64; D]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt()
}

/// Centripetal parameters over the points, on `[0, 1]`.
/// Fit a rectangular grid of points with a tensor-product B-spline surface.
///
/// The deferred grid fit, kept honest the way the curve fit is: rows first,
/// columns second, both passes at *fixed* parameters — a grid's
/// parameterization is shared property, and correcting it per row is how a
/// grid stops being one. Each pass adapts one shared knot vector against the
/// worst residual across its whole family, so every row rides the same
/// basis, which is what makes the second pass a fit over control points
/// rather than a guess. The reported error is measured at the end, surface
/// against every input point, and `met` does not round up.
///
/// `rows[j][i]` runs `i` along `u` and `j` along `v`.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the grid is
/// not rectangular or is smaller than two by two, or the tolerance is not a
/// distance.
pub fn fit_surface_grid(
    rows: &[Vec<Point>],
    degree: usize,
    tolerance: f64,
    tol: Tolerances,
) -> OgResult<Fitted<crate::BSplineSurface>> {
    use crate::traits::Surface as _;
    if !tolerance.is_finite() || tolerance <= 0.0 {
        og_bail!(Construction, "a tolerance of {tolerance} is not a distance");
    }
    let nv = rows.len();
    if nv < 2 {
        og_bail!(Construction, "a surface fit needs at least two rows");
    }
    let nu = rows[0].len();
    if nu < 2 || rows.iter().any(|r| r.len() != nu) {
        og_bail!(Construction, "a surface fit needs a rectangular grid");
    }
    let raw: Vec<Vec<[f64; 3]>> = rows
        .iter()
        .map(|r| r.iter().map(|p| [p.x, p.y, p.z]).collect())
        .collect();

    // Averaged centripetal parameters: one shared assignment per direction.
    let average = |families: &[Vec<[f64; 3]>]| -> Vec<f64> {
        let mut sums = vec![0.0; families[0].len()];
        for family in families {
            for (s, p) in sums.iter_mut().zip(centripetal::<3>(family)) {
                *s += p;
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let n = families.len() as f64;
        sums.iter().map(|s| s / n).collect()
    };
    let u_params = average(&raw);
    let columns: Vec<Vec<[f64; 3]>> = (0..nu)
        .map(|i| raw.iter().map(|r| r[i]).collect())
        .collect();
    let v_params = average(&columns);

    // Pass one: every row on one shared knot vector.
    let (u_knots, row_controls) = fit_family::<3>(&raw, &u_params, degree, tolerance * 0.5)?;
    // Pass two: the columns of control points, against the v parameters.
    let k = u_knots.control_point_count();
    let control_columns: Vec<Vec<[f64; 3]>> = (0..k)
        .map(|i| row_controls.iter().map(|r| r[i]).collect())
        .collect();
    let (v_knots, column_controls) =
        fit_family::<3>(&control_columns, &v_params, degree, tolerance * 0.5)?;
    let l = v_knots.control_point_count();

    // Assemble: `column_controls[i][j]` is the control at u-index i,
    // v-index j; the grid is row-major in u.
    let mut net: Vec<Point> = Vec::with_capacity(k * l);
    for column in &column_controls {
        for c in column {
            net.push(Point::new(c[0], c[1], c[2]));
        }
    }
    let grid = og_math::ControlGrid::new(net, k, l)?;
    let surface = crate::BSplineSurface::new(u_knots, v_knots, &grid, tol)?;

    // The honest error: the surface against every input point, at the
    // grid's own parameters.
    let mut worst = 0.0_f64;
    for (j, row) in rows.iter().enumerate() {
        for (i, p) in row.iter().enumerate() {
            let at = surface.point_at(u_params[i], v_params[j], tol)?;
            worst = worst.max(at.distance(*p));
        }
    }
    Ok(Fitted {
        curve: surface,
        error: worst,
        met: worst <= tolerance,
    })
}

/// Fit a family of point rows sharing parameters onto one knot vector,
/// refined against the worst residual across the whole family, parameters
/// held fixed.
/// The best round a family fit reached: knots, one control row per member,
/// and the worst residual.
type FamilyRound<const D: usize> = (KnotVector, Vec<Vec<[f64; D]>>, f64);

fn fit_family<const D: usize>(
    family: &[Vec<[f64; D]>],
    parameters: &[f64],
    degree: usize,
    tolerance: f64,
) -> OgResult<(KnotVector, Vec<Vec<[f64; D]>>)> {
    let degree = degree.min(parameters.len() - 1).max(1);
    let (a, b) = (parameters[0], parameters[parameters.len() - 1]);
    let mut knots = single_span(degree, a, b)?;
    const ROUNDS: usize = 24;
    let mut best: Option<FamilyRound<D>> = None;
    for _ in 0..ROUNDS {
        let mut controls = Vec::with_capacity(family.len());
        let mut merged: Vec<(f64, f64)> = parameters.iter().map(|u| (*u, 0.0)).collect();
        let mut solvable = true;
        for row in family {
            match least_squares::<D>(&knots, row, parameters) {
                Ok(control) => {
                    for (slot, entry) in residuals::<D>(&knots, &control, row, parameters)
                        .iter()
                        .zip(merged.iter_mut())
                    {
                        entry.1 = entry.1.max(slot.1);
                    }
                    controls.push(control);
                }
                Err(_) => {
                    solvable = false;
                    break;
                }
            }
        }
        if !solvable {
            break;
        }
        let worst = merged.iter().fold(0.0_f64, |acc, e| acc.max(e.1));
        if best.as_ref().is_none_or(|(_, _, held)| worst < *held) {
            best = Some((knots.clone(), controls, worst));
        }
        if worst <= tolerance || knots.control_point_count() >= parameters.len() {
            break;
        }
        let Some(refined) = refined_where_bad(&knots, &merged, tolerance)? else {
            break;
        };
        knots = refined;
    }
    let Some((knots, controls, _)) = best else {
        og_bail!(NotDone, "the family fit solved no round at all");
    };
    Ok((knots, controls))
}

fn centripetal<const D: usize>(points: &[[f64; D]]) -> Vec<f64> {
    let mut out = Vec::with_capacity(points.len());
    out.push(0.0);
    let mut total = 0.0;
    for pair in points.windows(2) {
        total += distance::<D>(&pair[0], &pair[1]).sqrt();
        out.push(total);
    }
    if total > 0.0 {
        for u in &mut out {
            *u /= total;
        }
    }
    // Exactness at the far end matters: the last parameter must be the domain
    // end, not a rounding neighbour of it.
    if let Some(last) = out.last_mut() {
        *last = 1.0;
    }
    out
}

/// A clamped knot vector with a single span.
fn single_span(degree: usize, a: f64, b: f64) -> OgResult<KnotVector> {
    let mut knots = Vec::with_capacity(2 * (degree + 1));
    knots.extend(core::iter::repeat_n(a, degree + 1));
    knots.extend(core::iter::repeat_n(b, degree + 1));
    KnotVector::new(knots, degree)
}

/// Least-squares control points for a fixed knot vector.
///
/// The ends are pinned to the first and last data points — they are where the
/// curve joins its neighbours — and the interior is solved. With no interior
/// there is nothing to solve and the pinned Bézier is the answer.
fn least_squares<const D: usize>(
    knots: &KnotVector,
    points: &[[f64; D]],
    parameters: &[f64],
) -> OgResult<Vec<[f64; D]>> {
    let n = knots.control_point_count();
    let m = points.len();
    let degree = knots.degree();
    let interior = n.saturating_sub(2);

    let mut control = vec![[0.0; D]; n];
    control[0] = points[0];
    control[n - 1] = points[m - 1];
    if interior == 0 {
        return Ok(control);
    }

    // The collocation rows, with the pinned ends moved to the right-hand side.
    let mut normal = nalgebra::DMatrix::<f64>::zeros(interior, interior);
    let mut rhs = vec![nalgebra::DVector::<f64>::zeros(interior); D];

    let mut rows: Vec<(usize, Vec<(usize, f64)>)> = Vec::with_capacity(m);
    for (k, &u) in parameters.iter().enumerate() {
        let span = knots.span_unchecked(u);
        let basis = knots.basis(span, u);
        let first = span - degree;
        rows.push((
            k,
            basis
                .iter()
                .enumerate()
                .map(|(j, b)| (first + j, *b))
                .collect(),
        ));
    }

    for (k, row) in &rows {
        // The residual this row wants to explain, after the pinned ends.
        let mut target = points[*k];
        for (index, b) in row {
            if *index == 0 {
                for d in 0..D {
                    target[d] -= b * points[0][d];
                }
            } else if *index == n - 1 {
                for d in 0..D {
                    target[d] -= b * points[m - 1][d];
                }
            }
        }
        for (i, bi) in row {
            if *i == 0 || *i == n - 1 {
                continue;
            }
            for (j, bj) in row {
                if *j == 0 || *j == n - 1 {
                    continue;
                }
                normal[(i - 1, j - 1)] += bi * bj;
            }
            for d in 0..D {
                rhs[d][i - 1] += bi * target[d];
            }
        }
    }

    // The normal matrix can be singular when a span has no parameter in it —
    // a knot was placed where there is no data to say where the curve goes.
    let Some(inverted) = normal.clone().try_inverse() else {
        og_bail!(
            NotDone,
            "the fitting system is singular: a knot span contains no data"
        );
    };
    for d in 0..D {
        let solved = &inverted * &rhs[d];
        for i in 0..interior {
            control[i + 1][d] = solved[i];
        }
    }
    Ok(control)
}

/// Move each parameter to the foot of the perpendicular from its point.
///
/// One Newton step per call on `g(u) = (C(u) - Q) · C'(u) = 0`. The ends stay
/// pinned: they are where the curve joins its neighbours, and letting them
/// slide would trade end accuracy for interior accuracy silently.
fn correct_parameters<const D: usize>(
    knots: &KnotVector,
    control: &[[f64; D]],
    points: &[[f64; D]],
    parameters: &mut [f64],
) {
    let (lo, hi) = knots.domain();
    let last = parameters.len() - 1;
    for (k, u) in parameters.iter_mut().enumerate() {
        if k == 0 || k == last {
            continue;
        }
        let (at, d1, d2) = evaluate::<D>(knots, control, *u);
        let gap: [f64; D] = core::array::from_fn(|d| at[d] - points[k][d]);
        let dot = |a: &[f64; D], b: &[f64; D]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
        let numerator = dot(&gap, &d1);
        let denominator = dot(&d1, &d1) + dot(&gap, &d2);
        if denominator.abs() <= f64::MIN_POSITIVE {
            continue;
        }
        let stepped = *u - numerator / denominator;
        if stepped.is_finite() {
            *u = stepped.clamp(lo, hi);
        }
    }
    // Projection can reorder neighbours near a tight turn; the fit assumes the
    // parameters walk forward with the points.
    for k in 1..parameters.len() {
        if parameters[k] < parameters[k - 1] {
            parameters[k] = parameters[k - 1];
        }
    }
}

/// A curve point and its first two derivatives, from raw knots and control.
fn evaluate<const D: usize>(
    knots: &KnotVector,
    control: &[[f64; D]],
    u: f64,
) -> ([f64; D], [f64; D], [f64; D]) {
    let degree = knots.degree();
    let span = knots.span_unchecked(u);
    let table = knots.basis_derivatives(span, u, 2);
    let first = span - degree;
    let mut out = [[0.0; D]; 3];
    for (order, row) in table.iter().enumerate().take(3) {
        for (j, b) in row.iter().enumerate() {
            for d in 0..D {
                out[order][d] += b * control[first + j][d];
            }
        }
    }
    (out[0], out[1], out[2])
}

/// The distance from each point to the curve at its parameter.
fn residuals<const D: usize>(
    knots: &KnotVector,
    control: &[[f64; D]],
    points: &[[f64; D]],
    parameters: &[f64],
) -> Vec<(f64, f64)> {
    let degree = knots.degree();
    parameters
        .iter()
        .zip(points)
        .map(|(&u, p)| {
            let span = knots.span_unchecked(u);
            let basis = knots.basis(span, u);
            let first = span - degree;
            let mut at = [0.0; D];
            for (j, b) in basis.iter().enumerate() {
                for d in 0..D {
                    at[d] += b * control[first + j][d];
                }
            }
            (u, distance::<D>(&at, p))
        })
        .collect()
}

/// The knot vector with every offending span split.
///
/// Split at the *median parameter* inside the span, not its geometric middle.
/// The middle is where a textbook puts it and it stalls in practice: when the
/// data crowds into one half of a bad span, a middle knot leaves the other
/// half empty, an empty half makes the least-squares system singular, and the
/// span can never be refined again — the fit then converges to just above the
/// target and sticks there. The median always leaves data on both sides.
///
/// `None` when no span can be split further — every bad span holds fewer than
/// two parameters, and a knot needs data on both sides to be supported.
fn refined_where_bad(
    knots: &KnotVector,
    errors: &[(f64, f64)],
    tolerance: f64,
) -> OgResult<Option<KnotVector>> {
    let distinct = knots.distinct();
    let mut refined = knots.clone();
    let mut changed = false;
    for window in distinct.windows(2) {
        let (lo, hi) = (window[0].0, window[1].0);
        let inside: Vec<f64> = errors
            .iter()
            .filter(|(u, _)| *u >= lo && *u < hi)
            .map(|(u, _)| *u)
            .collect();
        let bad = errors
            .iter()
            .any(|(u, e)| *u >= lo && *u < hi && *e > tolerance);
        if !bad || inside.len() < 2 {
            continue;
        }
        // The knot between the two middle parameters, kept strictly interior.
        let at = f64::midpoint(inside[inside.len() / 2 - 1], inside[inside.len() / 2])
            .clamp(lo + (hi - lo) * 1e-6, hi - (hi - lo) * 1e-6);
        // Both sides must keep data, or the new span is unsupported.
        let left = inside.iter().any(|u| *u < at);
        let right = inside.iter().any(|u| *u >= at);
        if left && right {
            refined = refined.with_knot_inserted(at, 1)?;
            changed = true;
        }
    }
    Ok(if changed { Some(refined) } else { None })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod grid_tests {
    use super::*;
    use crate::traits::Surface as _;
    use og_core::Tolerances;

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn a_torus_patch_grid_fits_to_tolerance_on_and_off_the_grid() {
        let torus = og_math::Torus::new(og_math::Frame::WORLD, 2.0, 0.5, T).unwrap();
        let surface = crate::TorusSurface::new(torus);
        let (nu, nv) = (25, 17);
        let span_u = 1.2_f64;
        let span_v = 0.9_f64;
        let sample = |fu: f64, fv: f64| surface.point_at(span_u * fu, span_v * fv, T).unwrap();
        let mut rows = Vec::new();
        for j in 0..nv {
            let mut row = Vec::new();
            for i in 0..nu {
                row.push(sample(
                    f64::from(i) / f64::from(nu - 1),
                    f64::from(j) / f64::from(nv - 1),
                ));
            }
            rows.push(row);
        }
        let fitted = fit_surface_grid(&rows, 3, 1e-4, T).unwrap();
        assert!(fitted.met, "error {} above the target", fitted.error);

        // Off the grid too: the fit describes the surface, not just its
        // samples. The fitted chart and the torus's differ, so compare by
        // distance to the true surface rather than at matched parameters.
        let (ud, vd) = fitted.curve.domain();
        for i in 0..8 {
            for j in 0..8 {
                let u = ud.0 + (ud.1 - ud.0) * (0.07 + 0.9 * f64::from(i) / 7.0);
                let v = vd.0 + (vd.1 - vd.0) * (0.07 + 0.9 * f64::from(j) / 7.0);
                let p = fitted.curve.point_at(u, v, T).unwrap();
                let d = torus.distance_to(p);
                assert!(d < 5e-4, "off-grid deviation {d} at ({u}, {v})");
            }
        }
    }

    #[test]
    fn a_grid_the_basis_can_represent_fits_to_rounding() {
        // Points from a bilinear patch: degree one in both directions.
        let corner = |x: f64, y: f64| Point::new(x, y, 0.3 * x - 0.2 * y);
        let mut rows = Vec::new();
        for j in 0..6 {
            let mut row = Vec::new();
            for i in 0..6 {
                row.push(corner(f64::from(i), 2.0 * f64::from(j)));
            }
            rows.push(row);
        }
        let fitted = fit_surface_grid(&rows, 1, 1e-9, T).unwrap();
        assert!(fitted.met, "error {} above rounding", fitted.error);
    }

    #[test]
    fn a_ragged_grid_is_refused() {
        let rows = vec![
            vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0)],
            vec![Point::new(0.0, 1.0, 0.0)],
        ];
        assert!(fit_surface_grid(&rows, 2, 1e-6, T).is_err());
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::traits::{Curve2d as _, Curve3d as _};
    use core::f64::consts::TAU;

    const T: Tolerances = Tolerances::millimetres();

    /// Distance from a point to a curve: scan, then ternary-refine.
    fn nearest(curve: &BSplineCurve, p: Point) -> f64 {
        let at = |u: f64| curve.point_at(u, T).map_or(f64::MAX, |q| p.distance(q));
        let mut best = (0.0, f64::MAX);
        for i in 0..=4000 {
            #[allow(clippy::cast_precision_loss)]
            let u = i as f64 / 4000.0;
            let d = at(u);
            if d < best.1 {
                best = (u, d);
            }
        }
        let (mut lo, mut hi) = ((best.0 - 5e-4).max(0.0), (best.0 + 5e-4).min(1.0));
        for _ in 0..100 {
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

    /// Points along a circle, the standard curve a spline cannot represent
    /// exactly and can approach as closely as asked.
    fn circle_points(n: usize, radius: f64) -> Vec<Point> {
        (0..=n)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let a = TAU * i as f64 / n as f64;
                Point::new(radius * a.cos(), radius * a.sin(), 0.0)
            })
            .collect()
    }

    #[test]
    fn a_fit_meets_the_tolerance_it_was_asked_for_and_says_what_it_reached() {
        let points = circle_points(200, 5.0);
        for tolerance in [1e-2, 1e-4, 1e-6] {
            let fitted = fit_points(&points, 3, tolerance, T).unwrap();
            assert!(fitted.met, "target {tolerance:e} not met");
            assert!(
                fitted.error <= tolerance,
                "reported {:e} over target {tolerance:e}",
                fitted.error
            );
            // And the report is honest: measure independently, point to curve.
            // A bare scan reports its own step size — 2000 samples over a
            // circumference of thirty is an 8e-3 grid, and comparing that to a
            // 1e-6 fit measures the scan — so the scan brackets and a ternary
            // search finishes.
            let mut worst = 0.0_f64;
            for p in &points {
                worst = worst.max(nearest(&fitted.curve, *p));
            }
            assert!(
                worst <= tolerance * 1.5,
                "independent measurement found {worst:e} against {tolerance:e}"
            );
        }
    }

    #[test]
    fn a_tighter_tolerance_never_uses_fewer_control_points() {
        let points = circle_points(300, 3.0);
        let coarse = fit_points(&points, 3, 1e-2, T).unwrap();
        let fine = fit_points(&points, 3, 1e-6, T).unwrap();
        assert!(
            fine.curve.control_points().len() > coarse.curve.control_points().len(),
            "{} then {}",
            coarse.curve.control_points().len(),
            fine.curve.control_points().len()
        );
        // And the coarse one is genuinely coarse: far fewer control points
        // than input points, or the fit is interpolation in disguise.
        assert!(coarse.curve.control_points().len() < 30);
    }

    #[test]
    fn knots_go_where_the_error_is() {
        // A straight run with one tight corner. Uniform refinement would
        // spread knots evenly; adaptive refinement must put them in the
        // corner.
        let mut points = Vec::new();
        for i in 0..=100 {
            points.push(Point::new(f64::from(i) * 0.1, 0.0, 0.0));
        }
        for i in 1..=50 {
            let a = f64::from(i) / 50.0 * core::f64::consts::FRAC_PI_2;
            points.push(Point::new(10.0 + a.sin() * 0.5, (1.0 - a.cos()) * 0.5, 0.0));
        }
        for i in 1..=100 {
            points.push(Point::new(10.5, 0.5 + f64::from(i) * 0.1, 0.0));
        }

        let fitted = fit_points(&points, 3, 1e-4, T).unwrap();
        assert!(fitted.met);

        // Knot parameters cluster around the corner, which sits at roughly
        // half way through the arc length.
        let distinct = fitted.curve.knots().distinct();
        let interior: Vec<f64> = distinct[1..distinct.len() - 1]
            .iter()
            .map(|(u, _)| *u)
            .collect();
        let near_corner = interior
            .iter()
            .filter(|u| (0.40..0.60).contains(*u))
            .count();
        assert!(
            near_corner * 2 > interior.len(),
            "only {near_corner} of {} interior knots are near the corner",
            interior.len()
        );
    }

    #[test]
    fn the_ends_are_honoured_exactly_and_a_closed_loop_stays_closed() {
        let points = circle_points(64, 2.0);
        let fitted = fit_points(&points, 3, 1e-3, T).unwrap();
        let (a, b) = fitted.curve.knots().domain();
        let start = fitted.curve.point_at(a, T).unwrap();
        let end = fitted.curve.point_at(b, T).unwrap();
        assert!(start.is_equal(points[0], T), "the start drifted");
        assert!(end.is_equal(*points.last().unwrap(), T), "the end drifted");
        assert!(start.is_equal(end, T), "the loop opened");
    }

    #[test]
    fn an_impossible_target_is_reported_not_rounded_up_to_success() {
        // Five points cannot pin a curve to a picometre unless the curve
        // interpolates them, and past that adding knots buys nothing. The fit
        // must say it fell short and how far.
        let points = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(2.0, -1.0, 0.0),
            Point::new(3.0, 1.0, 0.0),
            Point::new(4.0, 0.0, 0.0),
        ];
        let fitted = fit_points(&points, 3, 1e-15, T).unwrap();
        // With as many control points as points it may interpolate and land
        // at rounding; either way the flags must be consistent.
        assert_eq!(fitted.met, fitted.error <= 1e-15);
    }

    #[test]
    fn the_2d_fit_is_the_same_machinery() {
        let points: Vec<Point2> = (0..=100)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let a = TAU * f64::from(i) / 100.0;
                Point2::new(3.0 * a.cos(), 3.0 * a.sin())
            })
            .collect();
        let fitted = fit_points_2d(&points, 3, 1e-4, T).unwrap();
        assert!(fitted.met);
        assert!(fitted.error <= 1e-4);
        // Spot-check a few points independently, with refinement past the
        // scan's own resolution.
        for p in points.iter().step_by(17) {
            let scan = |u: f64| {
                fitted
                    .curve
                    .point_at(u, T)
                    .map_or(f64::MAX, |q| p.distance(q))
            };
            let mut best = (0.0, f64::MAX);
            for i in 0..=2000 {
                #[allow(clippy::cast_precision_loss)]
                let u = i as f64 / 2000.0;
                let d = scan(u);
                if d < best.1 {
                    best = (u, d);
                }
            }
            let (mut lo, mut hi) = ((best.0 - 1e-3).max(0.0), (best.0 + 1e-3).min(1.0));
            for _ in 0..100 {
                let one = lo + (hi - lo) / 3.0;
                let two = hi - (hi - lo) / 3.0;
                if scan(one) < scan(two) {
                    hi = two;
                } else {
                    lo = one;
                }
            }
            let found = scan(f64::midpoint(lo, hi)).min(best.1);
            assert!(found < 2e-4, "a 2d point is {found:e} off the fit");
        }
    }

    #[test]
    fn inputs_that_describe_nothing_are_refused() {
        let p = Point::ORIGIN;
        assert!(fit_points(&[], 3, 1e-3, T).is_err());
        assert!(fit_points(&[p], 3, 1e-3, T).is_err());
        assert!(
            fit_points(&[p, p, p], 3, 1e-3, T).is_err(),
            "all duplicates"
        );
        let two = [p, Point::new(1.0, 0.0, 0.0)];
        assert!(fit_points(&two, 3, 0.0, T).is_err());
        assert!(fit_points(&two, 3, -1.0, T).is_err());
        assert!(fit_points(&two, 3, f64::NAN, T).is_err());
        assert!(fit_points(&two, 0, 1e-3, T).is_err());
        // Two points always fit: the segment between them.
        assert!(fit_points(&two, 3, 1e-9, T).unwrap().met);
    }
}
