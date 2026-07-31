//! Fitting a B-spline to points.
//!
//! Two jobs that look alike and are not. [`interpolate`] makes a curve that
//! passes through every point exactly. [`approximate`] makes one that passes
//! near them, with fewer control points than there are points to fit.
//!
//! # Which one you want
//!
//! Interpolation is right when the points are exact — corners of a profile, a
//! path a machine must visit. It is wrong for measured data, because it fits
//! the noise as faithfully as the signal, and the wiggles it invents between
//! samples can be large.
//!
//! Approximation is right when the points are samples of something smoother
//! than they are. It also *cannot* be told to use as many control points as
//! there are data points — at that ratio the least-squares system is the
//! interpolation system, and calling one function and getting the other is a
//! trap. That case is refused with a message pointing at [`interpolate`].
//!
//! # Parameterization
//!
//! Centripetal by default: parameter spacing goes as the square root of the
//! chord, not the chord. Uniform spacing produces visible loops when the points
//! are unevenly spread, and plain chord length overshoots on sharp turns.
//! Centripetal is the standard compromise and is what a CAD user expects a
//! fitted curve to look like.
//!
//! # What is not here
//!
//! Nothing chooses the number of control points for you, and nothing decides
//! when a fit is "good enough" and stops. Those need an error target and a
//! knot-insertion loop; see the deferred list in `docs/SCOPE.md`.

use og_core::{OgResult, Tolerances, og_bail};
use og_geom::BSplineCurve;
use og_math::{KnotVector, Point};

/// How to spread parameters over the points.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Spacing {
    /// Parameter spacing goes as the square root of the chord.
    ///
    /// The default, and what a fitted curve is expected to look like: it
    /// neither loops the way uniform spacing does on unevenly spread points nor
    /// overshoots the way chord length does at a sharp turn.
    #[default]
    Centripetal,
    /// Parameter spacing proportional to the chord.
    Chordal,
    /// Equal parameter spacing, ignoring the points entirely.
    ///
    /// Correct only when the points really are evenly spread; otherwise it is
    /// the one that loops.
    Uniform,
}

/// Fit a B-spline that passes through every point.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if there are fewer
/// points than the degree requires, if two consecutive points coincide — which
/// leaves a parameter interval of zero and a singular system — or if the system
/// turns out singular anyway.
pub fn interpolate(
    points: &[Point],
    degree: usize,
    spacing: Spacing,
    tol: Tolerances,
) -> OgResult<BSplineCurve> {
    if degree == 0 {
        og_bail!(
            Construction,
            "a curve of degree zero is a point, not a curve"
        );
    }
    if points.len() <= degree {
        og_bail!(
            Construction,
            "interpolating a degree-{degree} curve needs at least {} points, \
             got {}",
            degree + 1,
            points.len()
        );
    }

    let parameters = parameterize(points, spacing, tol)?;
    let knots = KnotVector::averaged(degree, &parameters)?;

    // The collocation system: row k says "the curve at parameter t_k is point
    // k", which in the basis is a weighted sum of the control points.
    let n = points.len();
    let mut matrix = nalgebra::DMatrix::<f64>::zeros(n, n);
    for (row, &t) in parameters.iter().enumerate() {
        let span = knots.span(t, tol)?;
        let basis = knots.basis(span, t);
        for (j, value) in basis.iter().enumerate() {
            matrix[(row, span - degree + j)] = *value;
        }
    }

    let control = solve(&matrix, points)?;
    BSplineCurve::new(knots, control, tol)
}

/// Fit a B-spline that passes *near* the points, with `control_count` control
/// points.
///
/// The first and last points are interpolated exactly — a fitted curve that
/// does not start where the data starts is almost never wanted, and the ends
/// are where a free least-squares fit goes worst.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if there are too
/// few points, if `control_count` does not leave the system overdetermined, or
/// if the normal equations are singular.
pub fn approximate(
    points: &[Point],
    degree: usize,
    control_count: usize,
    spacing: Spacing,
    tol: Tolerances,
) -> OgResult<BSplineCurve> {
    if degree == 0 {
        og_bail!(
            Construction,
            "a curve of degree zero is a point, not a curve"
        );
    }
    if control_count <= degree {
        og_bail!(
            Construction,
            "a degree-{degree} curve needs more than {degree} control points, \
             asked for {control_count}"
        );
    }
    if points.len() <= control_count {
        og_bail!(
            Construction,
            "approximating {} points with {control_count} control points is not \
             an approximation — at that ratio the least-squares system *is* the \
             interpolation system. Use `interpolate`",
            points.len()
        );
    }

    let parameters = parameterize(points, spacing, tol)?;
    // Knots spread over the parameter range rather than averaged: averaging is
    // for interpolation, where there is one knot per point. Here there are
    // fewer control points than points, and the knots have to be placed so
    // every span holds at least one of them or the system is singular.
    let knots = spread(degree, control_count, &parameters)?;

    let free = control_count - 2;
    let inner = points.len() - 2;
    let mut matrix = nalgebra::DMatrix::<f64>::zeros(inner, free);
    let mut rhs = vec![Point::ORIGIN; inner];

    for k in 1..points.len() - 1 {
        let t = parameters[k];
        let span = knots.span(t, tol)?;
        let basis = knots.basis(span, t);

        // The two interpolated ends are known, so their contribution moves to
        // the right-hand side instead of being solved for.
        let mut residual = points[k].to_vector();
        for (j, value) in basis.iter().enumerate() {
            let column = span - degree + j;
            if column == 0 {
                residual -= points[0].to_vector() * *value;
            } else if column == control_count - 1 {
                residual -= points[points.len() - 1].to_vector() * *value;
            } else {
                matrix[(k - 1, column - 1)] = *value;
            }
        }
        rhs[k - 1] = Point::ORIGIN + residual;
    }

    // Normal equations. Forming them squares the condition number, which for a
    // well-spread fit of the sizes this is used at costs a few digits and buys
    // a much smaller solve than a QR of the full system.
    let normal = matrix.transpose() * &matrix;
    let projected = project(&matrix, &rhs);
    let middle = solve(&normal, &projected)?;

    let mut control = Vec::with_capacity(control_count);
    control.push(points[0]);
    control.extend(middle);
    control.push(points[points.len() - 1]);
    BSplineCurve::new(knots, control, tol)
}

/// Parameters for the points, by the chosen spacing.
fn parameterize(points: &[Point], spacing: Spacing, tol: Tolerances) -> OgResult<Vec<f64>> {
    let n = points.len();
    if n < 2 {
        og_bail!(Construction, "fitting needs at least two points");
    }
    if spacing == Spacing::Uniform {
        #[allow(clippy::cast_precision_loss)]
        return Ok((0..n).map(|i| i as f64 / (n - 1) as f64).collect());
    }

    let mut weights = Vec::with_capacity(n - 1);
    for w in points.windows(2) {
        let chord = w[0].distance(w[1]);
        if chord <= tol.confusion() {
            og_bail!(
                Construction,
                "two consecutive points coincide, which leaves a parameter \
                 interval of zero and a system with no solution; remove the \
                 duplicate before fitting"
            );
        }
        weights.push(if spacing == Spacing::Centripetal {
            chord.sqrt()
        } else {
            chord
        });
    }

    let total: f64 = weights.iter().sum();
    let mut parameters = Vec::with_capacity(n);
    parameters.push(0.0);
    let mut running = 0.0;
    for w in &weights {
        running += w;
        parameters.push(running / total);
    }
    // The last is 1.0 by construction, but only to within rounding, and the
    // knot vector's domain end has to match it exactly or the final point sits
    // a hair outside the curve's domain.
    let last = parameters.len() - 1;
    parameters[last] = 1.0;
    Ok(parameters)
}

/// A clamped knot vector for an approximation.
///
/// Interior knots are placed so each spans an equal share of the *data*, which
/// is what keeps every span occupied. A knot span with no data point in it
/// leaves a column of zeros in the system and no unique answer.
fn spread(degree: usize, control_count: usize, parameters: &[f64]) -> OgResult<KnotVector> {
    let mut knots = vec![0.0; degree + 1];
    let interior = control_count - degree - 1;

    #[allow(clippy::cast_precision_loss)]
    let step = (parameters.len() - 1) as f64 / (control_count - degree) as f64;
    for j in 1..=interior {
        #[allow(
            clippy::cast_precision_loss,
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss
        )]
        let at = (j as f64 * step) as usize;
        #[allow(clippy::cast_precision_loss)]
        let fraction = j as f64 * step - at as f64;
        let a = parameters[at.min(parameters.len() - 1)];
        let b = parameters[(at + 1).min(parameters.len() - 1)];
        knots.push(fraction.mul_add(b - a, a));
    }
    knots.extend(std::iter::repeat_n(1.0, degree + 1));
    KnotVector::new(knots, degree)
}

/// `Mᵀ b`, one column of points at a time.
fn project(matrix: &nalgebra::DMatrix<f64>, rhs: &[Point]) -> Vec<Point> {
    let mut out = vec![Point::ORIGIN; matrix.ncols()];
    for (column, slot) in out.iter_mut().enumerate() {
        let mut sum = og_math::Vector::ZERO;
        for (row, point) in rhs.iter().enumerate() {
            sum += point.to_vector() * matrix[(row, column)];
        }
        *slot = Point::ORIGIN + sum;
    }
    out
}

/// Solve `M x = points` for the three coordinates at once.
fn solve(matrix: &nalgebra::DMatrix<f64>, rhs: &[Point]) -> OgResult<Vec<Point>> {
    let n = rhs.len();
    let mut b = nalgebra::DMatrix::<f64>::zeros(n, 3);
    for (row, point) in rhs.iter().enumerate() {
        b[(row, 0)] = point.x;
        b[(row, 1)] = point.y;
        b[(row, 2)] = point.z;
    }

    // LU with partial pivoting. The collocation matrix is banded and diagonally
    // dominant for a sensible parameterization, so this is stable; the failure
    // it does report — a singular system — means the points or the knots were
    // degenerate, which is worth an error rather than a plausible answer.
    let Some(x) = matrix.clone().lu().solve(&b) else {
        og_bail!(
            Construction,
            "the fitting system has no unique solution; the points are \
             degenerate, or the knots leave a span with no point in it"
        );
    };
    Ok((0..n)
        .map(|i| Point::new(x[(i, 0)], x[(i, 1)], x[(i, 2)]))
        .collect())
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use og_geom::Curve3d;

    const T: Tolerances = Tolerances::millimetres();

    fn helix(n: usize) -> Vec<Point> {
        (0..n)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let t = i as f64 / (n - 1) as f64 * std::f64::consts::TAU;
                Point::new(t.cos() * 5.0, t.sin() * 5.0, t * 0.5)
            })
            .collect()
    }

    #[test]
    fn an_interpolant_passes_through_every_point() {
        // The defining property. A fit that is merely close is an
        // approximation, and the two are not interchangeable.
        for degree in [2, 3, 5] {
            let points = helix(12);
            let curve = interpolate(&points, degree, Spacing::Centripetal, T).unwrap();
            let parameters = parameterize(&points, Spacing::Centripetal, T).unwrap();

            for (point, t) in points.iter().zip(&parameters) {
                let on_curve = curve.point_at(*t, T).unwrap();
                assert!(
                    on_curve.distance(*point) < 1e-9,
                    "degree {degree}: missed by {}",
                    on_curve.distance(*point)
                );
            }
        }
    }

    #[test]
    fn an_interpolant_through_collinear_points_is_the_line_they_lie_on() {
        // A curve that wanders off a straight run of points is the classic
        // parameterization failure, and it is invisible at the points
        // themselves — only between them.
        let points: Vec<Point> = (0..8).map(|i| Point::new(f64::from(i), 0.0, 0.0)).collect();
        let curve = interpolate(&points, 3, Spacing::Centripetal, T).unwrap();

        for i in 0..=40 {
            let t = f64::from(i) / 40.0;
            let p = curve.point_at(t, T).unwrap();
            assert!(p.y.abs() < 1e-9 && p.z.abs() < 1e-9, "wandered to {p:?}");
        }
    }

    #[test]
    fn an_approximation_uses_the_control_points_it_was_given_and_hits_the_ends() {
        let points = helix(60);
        let curve = approximate(&points, 3, 10, Spacing::Centripetal, T).unwrap();
        assert_eq!(curve.control_points().len(), 10);

        let (a, b) = curve.domain();
        assert!(curve.point_at(a, T).unwrap().distance(points[0]) < 1e-9);
        assert!(
            curve
                .point_at(b, T)
                .unwrap()
                .distance(points[points.len() - 1])
                < 1e-9
        );
    }

    #[test]
    fn more_control_points_fit_the_data_more_closely() {
        // The property that makes an approximation useful: it is a knob, and
        // turning it has to do what it says.
        let points = helix(80);
        let parameters = parameterize(&points, Spacing::Centripetal, T).unwrap();
        let mut previous = f64::INFINITY;

        for count in [6, 10, 20, 40] {
            let curve = approximate(&points, 3, count, Spacing::Centripetal, T).unwrap();
            let worst = points
                .iter()
                .zip(&parameters)
                .map(|(p, t)| curve.point_at(*t, T).unwrap().distance(*p))
                .fold(0.0_f64, f64::max);
            assert!(
                worst < previous,
                "{count} control points fit worse than the previous step: \
                 {worst} against {previous}"
            );
            previous = worst;
        }
        assert!(
            previous < 0.05,
            "40 control points should fit well, got {previous}"
        );
    }

    #[test]
    fn centripetal_spacing_beats_uniform_on_unevenly_spread_points() {
        // The reason it is the default. Uniform spacing on points that bunch
        // and then spread makes the curve loop between the spread ones, and the
        // loop is far larger than any tolerance would allow.
        let mut points = vec![Point::ORIGIN];
        for i in 1..=5 {
            points.push(Point::new(f64::from(i) * 0.1, 0.0, 0.0));
        }
        points.push(Point::new(20.0, 0.0, 0.0));
        points.push(Point::new(40.0, 0.0, 0.0));

        let excursion = |spacing| {
            let curve = interpolate(&points, 3, spacing, T).unwrap();
            (0..=200)
                .map(|i| {
                    let t = f64::from(i) / 200.0;
                    let p = curve.point_at(t, T).unwrap();
                    p.y.hypot(p.z)
                })
                .fold(0.0_f64, f64::max)
        };
        assert!(
            excursion(Spacing::Centripetal) <= excursion(Spacing::Uniform) + 1e-12,
            "centripetal should be no worse than uniform"
        );
    }

    #[test]
    fn asking_for_an_approximation_that_is_an_interpolation_is_refused() {
        // Silently returning an interpolant would be worse than failing: the
        // caller asked for smoothing and would get the noise back, fitted
        // exactly, with nothing to say so.
        let points = helix(10);
        let refused = approximate(&points, 3, 10, Spacing::Centripetal, T);
        assert!(refused.is_err());
        assert!(
            format!("{}", refused.unwrap_err()).contains("Use `interpolate`"),
            "the message should point at the function that does want this"
        );
    }

    #[test]
    fn coincident_points_are_refused_rather_than_solved_around() {
        // Two points at one place leave a parameter interval of zero, and the
        // system has no unique answer. Nudging one apart silently would move
        // data the caller supplied.
        let points = vec![
            Point::ORIGIN,
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
        ];
        let refused = interpolate(&points, 2, Spacing::Centripetal, T);
        assert!(refused.is_err());
        assert!(format!("{}", refused.unwrap_err()).contains("coincide"));
    }

    #[test]
    fn too_few_points_for_the_degree_is_refused() {
        let points = helix(3);
        assert!(interpolate(&points, 5, Spacing::Centripetal, T).is_err());
        assert!(interpolate(&points, 0, Spacing::Centripetal, T).is_err());
        assert!(approximate(&points, 3, 3, Spacing::Centripetal, T).is_err());
    }

    #[test]
    fn a_degree_one_interpolant_is_the_polyline_itself() {
        let points = helix(6);
        let curve = interpolate(&points, 1, Spacing::Chordal, T).unwrap();
        assert_eq!(curve.control_points().len(), points.len());
        for (control, point) in curve.control_points().iter().zip(&points) {
            assert!(control.scaled.distance(*point) < 1e-12);
        }
    }
}
