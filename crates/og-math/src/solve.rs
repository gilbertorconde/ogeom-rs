//! Root finding, polynomial roots and minimization.
//!
//! The numerical substrate every intersection, projection and extrema algorithm
//! in the kernel sits on. Nothing here is geometric; it is deliberately kept
//! separate so those algorithms are about geometry rather than about
//! convergence.
//!
//! # What to reach for
//!
//! - A root inside a known bracket: [`brent`]. Guaranteed to converge, and
//!   nearly as fast as Newton in practice.
//! - A root with a known derivative and a good starting point: [`newton`],
//!   which falls back to bisection whenever a step would leave the bracket.
//!   Unsafeguarded Newton diverges on the configurations that matter — a
//!   tangential intersection is exactly where the derivative vanishes.
//! - Roots of a polynomial up to quartic: [`roots`]. Closed form, and the
//!   quadratic is written to avoid the cancellation the schoolbook formula
//!   suffers.
//! - A system of equations: [`newton_system`]. Surface projection is two
//!   equations in two unknowns; intersection marching is much the same.
//! - A minimum without derivatives: [`minimize`].

use nalgebra::{DMatrix, DVector};
use og_core::{OgResult, og_bail};

/// How a solver finished.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Convergence {
    /// The residual fell below the requested tolerance.
    Residual,
    /// The step size fell below the requested tolerance.
    Step,
    /// The iteration limit was reached first. The result is the best estimate
    /// found, and is *not* to be treated as a root.
    Exhausted,
}

impl Convergence {
    /// Whether the solver actually converged.
    #[must_use]
    pub const fn is_converged(self) -> bool {
        !matches!(self, Self::Exhausted)
    }
}

/// A solver result: the estimate, the residual there, and how it finished.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Solution {
    /// The parameter value.
    pub value: f64,
    /// The function's value there.
    pub residual: f64,
    /// How the iteration ended.
    pub convergence: Convergence,
    /// Iterations taken.
    pub iterations: usize,
}

/// Stopping criteria for an iterative solver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Criteria {
    /// Stop when `|f(x)|` falls to this.
    pub residual: f64,
    /// Stop when the step falls to this.
    pub step: f64,
    /// Give up after this many iterations.
    pub max_iterations: usize,
}

impl Default for Criteria {
    fn default() -> Self {
        Self {
            // Tight enough for geometric work in f64 without chasing the last
            // couple of bits, which costs iterations and buys nothing.
            residual: 1e-13,
            step: 1e-14,
            max_iterations: 100,
        }
    }
}

impl Criteria {
    /// Criteria with a given residual tolerance and the default step limit.
    #[must_use]
    pub fn with_residual(residual: f64) -> Self {
        Self {
            residual,
            ..Self::default()
        }
    }
}

/// Find a root of `f` in `[a, b]` by Brent's method.
///
/// Combines bisection, the secant method and inverse quadratic interpolation,
/// taking whichever step is both safe and fast. Guaranteed to converge for a
/// continuous function that changes sign across the bracket, and superlinear in
/// practice — the right default when a bracket is available.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the bracket is
/// malformed, or if `f` does not change sign across it, which means there is no
/// root to find by this method.
pub fn brent<F>(mut f: F, a: f64, b: f64, criteria: Criteria) -> OgResult<Solution>
where
    F: FnMut(f64) -> f64,
{
    if !a.is_finite() || !b.is_finite() || a >= b {
        og_bail!(Construction, "bracket [{a}, {b}] is empty or non-finite");
    }
    let (mut fa, mut fb) = (f(a), f(b));
    if fa == 0.0 {
        return Ok(Solution {
            value: a,
            residual: 0.0,
            convergence: Convergence::Residual,
            iterations: 0,
        });
    }
    if fb == 0.0 {
        return Ok(Solution {
            value: b,
            residual: 0.0,
            convergence: Convergence::Residual,
            iterations: 0,
        });
    }
    if fa * fb > 0.0 {
        og_bail!(
            Construction,
            "f does not change sign across [{a}, {b}]: f(a) = {fa}, f(b) = {fb}"
        );
    }

    let (mut a, mut b) = (a, b);
    // `b` is kept as the better estimate throughout.
    if fa.abs() < fb.abs() {
        core::mem::swap(&mut a, &mut b);
        core::mem::swap(&mut fa, &mut fb);
    }
    let mut c = a;
    let mut fc = fa;
    let mut previous_step = b - a;
    let mut used_bisection = true;

    for iteration in 1..=criteria.max_iterations {
        let mut s = if fa != fc && fb != fc {
            // Inverse quadratic interpolation, when three distinct values allow
            // it.
            a * fb * fc / ((fa - fb) * (fa - fc))
                + b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
        } else {
            b - fb * (b - a) / (fb - fa)
        };

        // Bisect instead whenever the interpolated step is outside the bracket
        // or is not shrinking fast enough. This is what turns a fast but
        // unreliable method into a guaranteed one.
        let bounds = ((3.0 * a + b) / 4.0, b);
        let outside = if bounds.0 < bounds.1 {
            s < bounds.0 || s > bounds.1
        } else {
            s < bounds.1 || s > bounds.0
        };
        let step = (s - b).abs();
        let stalled = if used_bisection {
            step >= (b - c).abs() / 2.0
        } else {
            step >= previous_step.abs() / 2.0
        };
        if outside || stalled || previous_step.abs() < criteria.step {
            s = f64::midpoint(a, b);
            used_bisection = true;
        } else {
            used_bisection = false;
        }

        let fs = f(s);
        previous_step = b - c;
        c = b;
        fc = fb;
        if fa * fs < 0.0 {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        if fa.abs() < fb.abs() {
            core::mem::swap(&mut a, &mut b);
            core::mem::swap(&mut fa, &mut fb);
        }

        if fb.abs() <= criteria.residual {
            return Ok(Solution {
                value: b,
                residual: fb,
                convergence: Convergence::Residual,
                iterations: iteration,
            });
        }
        if (b - a).abs() <= criteria.step {
            return Ok(Solution {
                value: b,
                residual: fb,
                convergence: Convergence::Step,
                iterations: iteration,
            });
        }
    }
    Ok(Solution {
        value: b,
        residual: fb,
        convergence: Convergence::Exhausted,
        iterations: criteria.max_iterations,
    })
}

/// Find a root of `f` near `start`, using its derivative, safeguarded by a
/// bracket.
///
/// Takes a Newton step when that lands inside `[a, b]` and reduces the residual,
/// and bisects otherwise. Plain Newton is not usable here: it diverges wherever
/// the derivative is small, and small derivatives are precisely the tangential
/// configurations a geometry kernel spends its time on.
///
/// `f` returns the value and the derivative together, since evaluating them
/// separately usually repeats most of the work.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the bracket is
/// malformed or `f` does not change sign across it.
pub fn newton<F>(mut f: F, a: f64, b: f64, start: f64, criteria: Criteria) -> OgResult<Solution>
where
    F: FnMut(f64) -> (f64, f64),
{
    if !a.is_finite() || !b.is_finite() || a >= b {
        og_bail!(Construction, "bracket [{a}, {b}] is empty or non-finite");
    }
    let (mut low, mut high) = (a, b);
    let (fa, _) = f(low);
    let (fb, _) = f(high);
    if fa == 0.0 {
        return Ok(Solution {
            value: low,
            residual: 0.0,
            convergence: Convergence::Residual,
            iterations: 0,
        });
    }
    if fb == 0.0 {
        return Ok(Solution {
            value: high,
            residual: 0.0,
            convergence: Convergence::Residual,
            iterations: 0,
        });
    }
    if fa * fb > 0.0 {
        og_bail!(Construction, "f does not change sign across [{a}, {b}]");
    }
    // Orient so that f(low) < 0 < f(high); the bracket update is then a single
    // comparison rather than a sign product.
    if fa > 0.0 {
        core::mem::swap(&mut low, &mut high);
    }

    let mut x = start.clamp(a, b);
    let mut previous_step = (b - a).abs();

    for iteration in 1..=criteria.max_iterations {
        let (value, slope) = f(x);
        if value.abs() <= criteria.residual {
            return Ok(Solution {
                value: x,
                residual: value,
                convergence: Convergence::Residual,
                iterations: iteration,
            });
        }
        if value < 0.0 {
            low = x;
        } else {
            high = x;
        }

        let newton_step = if slope == 0.0 {
            f64::INFINITY
        } else {
            value / slope
        };
        let candidate = x - newton_step;
        let out_of_bracket = (candidate - low) * (candidate - high) > 0.0;
        // A step that has not at least halved is a sign Newton is not making
        // progress here, so fall back rather than grind.
        let too_slow = (2.0 * newton_step).abs() > previous_step;

        let next = if out_of_bracket || too_slow || !candidate.is_finite() {
            f64::midpoint(low, high)
        } else {
            candidate
        };
        previous_step = (next - x).abs();
        x = next;

        if previous_step <= criteria.step {
            let (residual, _) = f(x);
            return Ok(Solution {
                value: x,
                residual,
                convergence: Convergence::Step,
                iterations: iteration,
            });
        }
    }
    let (residual, _) = f(x);
    Ok(Solution {
        value: x,
        residual,
        convergence: Convergence::Exhausted,
        iterations: criteria.max_iterations,
    })
}

/// The real roots of a polynomial, in increasing order.
///
/// `coefficients` are in ascending power order: `c[0] + c[1] x + c[2] x^2 ...`.
/// Degrees up to four are solved in closed form; above that the polynomial is
/// deflated by companion-matrix eigenvalues.
///
/// Repeated roots are returned once each, since a geometry caller wants the
/// distinct parameter values.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if every
/// coefficient is zero, where every value is a root.
pub fn roots(coefficients: &[f64], tolerance: f64) -> OgResult<Vec<f64>> {
    // Drop leading zeros so the true degree drives the choice of method; a
    // "cubic" whose cubic term is zero is a quadratic and must be solved as
    // one, or the leading division blows up.
    let mut c = coefficients;
    while let Some((&last, rest)) = c.split_last() {
        if last.abs() <= tolerance * c.iter().fold(0.0_f64, |m, v| m.max(v.abs())).max(1.0) {
            c = rest;
        } else {
            break;
        }
    }

    let mut out = match c.len() {
        0 => og_bail!(
            Construction,
            "the zero polynomial has every value as a root"
        ),
        1 => Vec::new(),
        2 => vec![-c[0] / c[1]],
        3 => quadratic_roots(c[2], c[1], c[0]),
        4 => cubic_roots(c[3], c[2], c[1], c[0]),
        _ => companion_roots(c, tolerance),
    };

    out.retain(|r| r.is_finite());
    out.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
    out.dedup_by(|a, b| (*a - *b).abs() <= tolerance * a.abs().max(1.0));
    Ok(out)
}

/// Real roots of `a x^2 + b x + c`.
///
/// Uses the citardauq form for whichever root would otherwise be computed as a
/// difference of nearly equal numbers. The schoolbook formula loses most of its
/// precision for the smaller root when `b^2 >> 4ac` — which is the common case
/// for a ray grazing a sphere.
#[must_use]
pub fn quadratic_roots(a: f64, b: f64, c: f64) -> Vec<f64> {
    if a == 0.0 {
        return if b == 0.0 { Vec::new() } else { vec![-c / b] };
    }
    let discriminant = b.mul_add(b, -(4.0 * a * c));
    if discriminant < 0.0 {
        return Vec::new();
    }
    if discriminant == 0.0 {
        return vec![-b / (2.0 * a)];
    }
    let sqrt = discriminant.sqrt();
    // Add magnitudes rather than subtract them, then get the other root from
    // the product relation.
    let q = -0.5 * (b + b.signum() * sqrt);
    let (r1, r2) = (q / a, if q == 0.0 { 0.0 } else { c / q });
    if r1 <= r2 { vec![r1, r2] } else { vec![r2, r1] }
}

/// Real roots of `a x^3 + b x^2 + c x + d`.
///
/// Depressed cubic plus the trigonometric solution in the three-real-roots
/// case, which avoids the complex arithmetic Cardano's formula would otherwise
/// need there.
#[must_use]
pub fn cubic_roots(a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
    if a == 0.0 {
        return quadratic_roots(b, c, d);
    }
    let (b, c, d) = (b / a, c / a, d / a);
    let shift = b / 3.0;
    // x = t - b/3 removes the quadratic term.
    let p = shift.mul_add(-b, c);
    let q = (2.0 / 27.0 * b * b).mul_add(b, shift.mul_add(-c, d));

    let half_q = q / 2.0;
    let third_p = p / 3.0;
    let discriminant = half_q.mul_add(half_q, third_p * third_p * third_p);

    if discriminant > 0.0 {
        let sqrt = discriminant.sqrt();
        let u = (-half_q + sqrt).cbrt();
        let v = (-half_q - sqrt).cbrt();
        vec![u + v - shift]
    } else if discriminant == 0.0 {
        if p == 0.0 {
            vec![-shift]
        } else {
            let u = (-half_q).cbrt();
            let mut r = vec![2.0 * u - shift, -u - shift];
            r.sort_by(|x, y| x.partial_cmp(y).unwrap_or(core::cmp::Ordering::Equal));
            r
        }
    } else {
        // Three distinct real roots, via trigonometry.
        let radius = (-third_p).sqrt();
        let cos = (-half_q / (radius * radius * radius)).clamp(-1.0, 1.0);
        let angle = cos.acos() / 3.0;
        let scale = 2.0 * radius;
        let tau_third = core::f64::consts::TAU / 3.0;
        let mut r = vec![
            scale.mul_add(angle.cos(), -shift),
            scale.mul_add((angle - tau_third).cos(), -shift),
            scale.mul_add((angle + tau_third).cos(), -shift),
        ];
        r.sort_by(|x, y| x.partial_cmp(y).unwrap_or(core::cmp::Ordering::Equal));
        r
    }
}

/// Real roots of a polynomial of any degree, via companion-matrix eigenvalues.
fn companion_roots(c: &[f64], tolerance: f64) -> Vec<f64> {
    let n = c.len() - 1;
    let lead = c[n];
    let mut m = DMatrix::<f64>::zeros(n, n);
    for i in 0..n {
        m[(i, n - 1)] = -c[i] / lead;
        if i + 1 < n {
            m[(i + 1, i)] = 1.0;
        }
    }
    // Only the real eigenvalues are roots; complex conjugate pairs are not.
    m.complex_eigenvalues()
        .iter()
        .filter(|e| e.im.abs() <= tolerance.max(1e-9) * e.re.abs().max(1.0))
        .map(|e| e.re)
        .collect()
}

/// Minimize a scalar function on `[a, b]` without derivatives.
///
/// Brent's method again: golden-section search with parabolic interpolation
/// wherever the parabola is well behaved. Converges for any continuous function
/// and is not fooled by the flat regions near a minimum, where a derivative
/// method has nothing to work with.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the bracket is
/// malformed.
pub fn minimize<F>(mut f: F, a: f64, b: f64, criteria: Criteria) -> OgResult<Solution>
where
    F: FnMut(f64) -> f64,
{
    if !a.is_finite() || !b.is_finite() || a >= b {
        og_bail!(Construction, "bracket [{a}, {b}] is empty or non-finite");
    }
    // 1 - 1/phi: the golden-section step.
    const GOLDEN: f64 = 0.381_966_011_250_105_15;

    let (mut low, mut high) = (a, b);
    let mut x = GOLDEN.mul_add(b - a, a);
    let (mut w, mut v) = (x, x);
    let mut fx = f(x);
    let (mut fw, mut fv) = (fx, fx);
    let mut step = 0.0_f64;
    let mut previous_step = 0.0_f64;

    for iteration in 1..=criteria.max_iterations {
        let middle = f64::midpoint(low, high);
        let tolerance = criteria.step.mul_add(x.abs(), criteria.step);
        if (x - middle).abs() <= 2.0f64.mul_add(tolerance, -((high - low) / 2.0)) {
            return Ok(Solution {
                value: x,
                residual: fx,
                convergence: Convergence::Step,
                iterations: iteration,
            });
        }

        let mut use_golden = true;
        if previous_step.abs() > tolerance {
            // Fit a parabola through the three best points so far.
            let r = (x - w) * (fx - fv);
            let q = (x - v) * (fx - fw);
            let mut p = (x - v) * q - (x - w) * r;
            let mut q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            }
            q = q.abs();
            // Accept the parabolic step only if it stays inside the bracket and
            // is smaller than half the step before last.
            if p.abs() < (0.5 * q * previous_step).abs() && p > q * (low - x) && p < q * (high - x)
            {
                step = p / q;
                let candidate = x + step;
                if candidate - low < 2.0 * tolerance || high - candidate < 2.0 * tolerance {
                    step = if x < middle { tolerance } else { -tolerance };
                }
                use_golden = false;
            }
        }
        if use_golden {
            previous_step = if x < middle { high - x } else { low - x };
            step = GOLDEN * previous_step;
        }

        let next = if step.abs() >= tolerance {
            x + step
        } else if step > 0.0 {
            x + tolerance
        } else {
            x - tolerance
        };
        let fnext = f(next);

        if fnext <= fx {
            if next < x {
                high = x;
            } else {
                low = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = next;
            fx = fnext;
        } else {
            if next < x {
                low = next;
            } else {
                high = next;
            }
            if fnext <= fw || w == x {
                v = w;
                fv = fw;
                w = next;
                fw = fnext;
            } else if fnext <= fv || v == x || v == w {
                v = next;
                fv = fnext;
            }
        }
        previous_step = step;
    }
    Ok(Solution {
        value: x,
        residual: fx,
        convergence: Convergence::Exhausted,
        iterations: criteria.max_iterations,
    })
}

/// The outcome of solving a system of equations.
#[derive(Debug, Clone, PartialEq)]
pub struct SystemSolution {
    /// The estimate.
    pub value: Vec<f64>,
    /// The norm of the residual vector there.
    pub residual: f64,
    /// How the iteration ended.
    pub convergence: Convergence,
    /// Iterations taken.
    pub iterations: usize,
}

/// Solve `f(x) = 0` for a vector `x`, by damped Newton.
///
/// `f` returns the residual vector and the Jacobian, row-major. The step is
/// halved until it actually reduces the residual — undamped Newton overshoots
/// badly from a poor start, and a geometry caller's start is often only a rough
/// guess from a coarse sampling.
///
/// Surface projection is this with two equations in two unknowns; so is a step
/// of a surface/surface intersection march.
///
/// Where no root exists the residual has a positive minimum, and no damping
/// finds a downhill step from it. That is reported as
/// [`Convergence::Exhausted`] with the best estimate attached — "no root here"
/// is a useful answer, and far better than iterating to the limit.
///
/// # Errors
///
/// [`OgError::Dimension`](og_core::OgError::Dimension) if the Jacobian's shape
/// disagrees with the residual, and
/// [`OgError::Numeric`](og_core::OgError::Numeric) if the Jacobian is singular
/// at the starting point.
pub fn newton_system<F>(mut f: F, start: &[f64], criteria: Criteria) -> OgResult<SystemSolution>
where
    F: FnMut(&[f64]) -> (Vec<f64>, Vec<Vec<f64>>),
{
    let n = start.len();
    let mut x = DVector::from_row_slice(start);

    let evaluate = |x: &DVector<f64>, f: &mut F| {
        let (r, j) = f(x.as_slice());
        (r, j)
    };

    let (mut residual, mut jacobian) = evaluate(&x, &mut f);
    if residual.len() != n || jacobian.len() != n || jacobian.iter().any(|row| row.len() != n) {
        og_bail!(
            Dimension,
            "expected a {n}-vector residual and {n}x{n} Jacobian"
        );
    }
    let mut norm = residual.iter().map(|v| v * v).sum::<f64>().sqrt();

    for iteration in 1..=criteria.max_iterations {
        if norm <= criteria.residual {
            return Ok(SystemSolution {
                value: x.as_slice().to_vec(),
                residual: norm,
                convergence: Convergence::Residual,
                iterations: iteration - 1,
            });
        }

        let j = DMatrix::from_fn(n, n, |r, c| jacobian[r][c]);
        let rhs = DVector::from_row_slice(&residual);
        let Some(delta) = j.lu().solve(&rhs) else {
            og_bail!(Numeric, "Jacobian is singular after {iteration} iterations");
        };

        // Damping: keep halving until the residual actually falls. Without it,
        // Newton happily steps past the solution and never comes back.
        let mut scale = 1.0;
        let mut accepted = None;
        for _ in 0..30 {
            let candidate = &x - &delta * scale;
            let (r, jj) = evaluate(&candidate, &mut f);
            let candidate_norm = r.iter().map(|v| v * v).sum::<f64>().sqrt();
            if candidate_norm < norm || candidate_norm <= criteria.residual {
                accepted = Some((candidate, r, jj, candidate_norm));
                break;
            }
            scale *= 0.5;
        }

        let Some((next, r, jj, next_norm)) = accepted else {
            // No downhill step exists: this is a local minimum of the residual,
            // not a root, and saying so is more useful than iterating forever.
            return Ok(SystemSolution {
                value: x.as_slice().to_vec(),
                residual: norm,
                convergence: Convergence::Exhausted,
                iterations: iteration,
            });
        };

        let step = (&next - &x).norm();
        x = next;
        residual = r;
        jacobian = jj;
        norm = next_norm;

        if norm <= criteria.residual {
            return Ok(SystemSolution {
                value: x.as_slice().to_vec(),
                residual: norm,
                convergence: Convergence::Residual,
                iterations: iteration,
            });
        }
        if step <= criteria.step {
            return Ok(SystemSolution {
                value: x.as_slice().to_vec(),
                residual: norm,
                convergence: Convergence::Step,
                iterations: iteration,
            });
        }
    }
    Ok(SystemSolution {
        value: x.as_slice().to_vec(),
        residual: norm,
        convergence: Convergence::Exhausted,
        iterations: criteria.max_iterations,
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const C: Criteria = Criteria {
        residual: 1e-13,
        step: 1e-14,
        max_iterations: 100,
    };

    #[test]
    fn brent_finds_a_simple_root() {
        let s = brent(|x| x * x - 2.0, 0.0, 2.0, C).unwrap();
        assert!(s.convergence.is_converged());
        assert_relative_eq!(s.value, core::f64::consts::SQRT_2, epsilon = 1e-12);
    }

    #[test]
    fn brent_handles_a_root_at_a_bracket_end() {
        let s = brent(|x| x, -1.0, 0.0, C).unwrap();
        assert_relative_eq!(s.value, 0.0);
        assert_eq!(s.iterations, 0);
    }

    #[test]
    fn brent_refuses_a_bracket_without_a_sign_change() {
        assert!(brent(|x| x * x + 1.0, -1.0, 1.0, C).is_err());
        assert!(brent(|x| x, 1.0, 0.0, C).is_err(), "reversed bracket");
        assert!(brent(|x| x, 0.0, f64::NAN, C).is_err());
    }

    #[test]
    fn brent_converges_on_a_function_that_defeats_the_secant_method() {
        // Very flat near the root, then steep: pure secant crawls, bisection
        // alone is slow, and the hybrid must do better than either.
        let s = brent(|x| x.powi(15) - 0.5, 0.0, 2.0, C).unwrap();
        assert!(s.convergence.is_converged());
        assert!(s.residual.abs() < 1e-12);
        assert!(s.iterations < 60, "took {} iterations", s.iterations);
    }

    #[test]
    fn newton_converges_faster_than_bisection_when_it_can() {
        let s = newton(|x| (x * x - 2.0, 2.0 * x), 0.5, 2.0, 1.0, C).unwrap();
        assert!(s.convergence.is_converged());
        assert_relative_eq!(s.value, core::f64::consts::SQRT_2, epsilon = 1e-12);
        assert!(s.iterations < 12, "took {} iterations", s.iterations);
    }

    #[test]
    fn newton_survives_a_vanishing_derivative() {
        // f(x) = x^3 has f'(0) = 0. Unsafeguarded Newton stalls; the bisection
        // fallback must carry it through.
        let s = newton(|x| (x * x * x, 3.0 * x * x), -1.0, 2.0, 1.9, C).unwrap();
        assert!(s.value.abs() < 1e-4, "landed at {}", s.value);
    }

    #[test]
    fn newton_survives_a_terrible_starting_point() {
        for start in [-0.999_f64, 0.0, 1.999, 1.0] {
            let s = newton(|x| (x * x - 2.0, 2.0 * x), -1.0, 2.0, start, C).unwrap();
            assert!(
                (s.value - core::f64::consts::SQRT_2).abs() < 1e-9,
                "start {start} gave {}",
                s.value
            );
        }
    }

    #[test]
    fn quadratic_roots_stay_accurate_when_the_roots_are_far_apart() {
        // x^2 - (1e8 + 1e-8) x + 1 has roots 1e8 and 1e-8. The schoolbook
        // formula computes the small one as a difference of nearly equal
        // numbers and loses almost all of it.
        let r = quadratic_roots(1.0, -(1e8 + 1e-8), 1.0);
        assert_eq!(r.len(), 2);
        assert_relative_eq!(r[0], 1e-8, max_relative = 1e-10);
        assert_relative_eq!(r[1], 1e8, max_relative = 1e-14);
    }

    #[test]
    fn quadratic_edge_cases() {
        assert_eq!(
            quadratic_roots(1.0, 0.0, 1.0),
            Vec::<f64>::new(),
            "no real roots"
        );
        assert_eq!(quadratic_roots(1.0, -2.0, 1.0), vec![1.0], "double root");
        assert_eq!(
            quadratic_roots(0.0, 2.0, -4.0),
            vec![2.0],
            "degenerates to linear"
        );
        assert_eq!(quadratic_roots(0.0, 0.0, 1.0), Vec::<f64>::new());
        let r = quadratic_roots(1.0, 0.0, -4.0);
        assert_relative_eq!(r[0], -2.0);
        assert_relative_eq!(r[1], 2.0);
    }

    #[test]
    fn cubic_with_three_real_roots() {
        // (x + 3)(x - 1)(x - 2) = x^3 - 7x + 6
        let r = cubic_roots(1.0, 0.0, -7.0, 6.0);
        assert_eq!(r.len(), 3);
        assert_relative_eq!(r[0], -3.0, epsilon = 1e-12);
        assert_relative_eq!(r[1], 1.0, epsilon = 1e-12);
        assert_relative_eq!(r[2], 2.0, epsilon = 1e-12);
    }

    #[test]
    fn cubic_with_one_real_root() {
        // x^3 + x + 1 has a single real root near -0.6823
        let r = cubic_roots(1.0, 0.0, 1.0, 1.0);
        assert_eq!(r.len(), 1);
        assert_relative_eq!(r[0], -0.682_327_803_828_019_3, epsilon = 1e-12);
    }

    #[test]
    fn cubic_with_repeated_roots() {
        // (x - 2)^2 (x + 1) = x^3 - 3x^2 + 4
        let r = cubic_roots(1.0, -3.0, 0.0, 4.0);
        assert_eq!(r.len(), 2, "a repeated root is reported once");
        assert_relative_eq!(r[0], -1.0, epsilon = 1e-9);
        assert_relative_eq!(r[1], 2.0, epsilon = 1e-9);
        // A triple root.
        let t = cubic_roots(1.0, 0.0, 0.0, 0.0);
        assert_eq!(t, vec![0.0]);
    }

    #[test]
    fn roots_strips_leading_zeros_before_choosing_a_method() {
        // Presented as a cubic, but with a zero cubic term: solving it as one
        // would divide by zero.
        let r = roots(&[-4.0, 0.0, 1.0, 0.0], 1e-12).unwrap();
        assert_eq!(r.len(), 2);
        assert_relative_eq!(r[0], -2.0, epsilon = 1e-12);
        assert_relative_eq!(r[1], 2.0, epsilon = 1e-12);
    }

    #[test]
    fn roots_of_a_quartic() {
        // (x-1)(x-2)(x-3)(x-4) = x^4 - 10x^3 + 35x^2 - 50x + 24
        let r = roots(&[24.0, -50.0, 35.0, -10.0, 1.0], 1e-9).unwrap();
        assert_eq!(r.len(), 4);
        for (got, want) in r.iter().zip([1.0, 2.0, 3.0, 4.0]) {
            assert_relative_eq!(got, &want, epsilon = 1e-7);
        }
    }

    #[test]
    fn roots_of_a_high_degree_polynomial() {
        // (x-1)(x-2)(x-3)(x-4)(x-5)
        let r = roots(&[-120.0, 274.0, -225.0, 85.0, -15.0, 1.0], 1e-9).unwrap();
        assert_eq!(r.len(), 5);
        for (got, want) in r.iter().zip([1.0, 2.0, 3.0, 4.0, 5.0]) {
            assert_relative_eq!(got, &want, epsilon = 1e-6);
        }
    }

    #[test]
    fn roots_degenerate_cases() {
        assert!(roots(&[], 1e-12).is_err());
        assert!(roots(&[0.0, 0.0], 1e-12).is_err());
        assert_eq!(
            roots(&[5.0], 1e-12).unwrap(),
            Vec::<f64>::new(),
            "a nonzero constant"
        );
        assert_eq!(roots(&[0.0, 1.0], 1e-12).unwrap(), vec![0.0]);
    }

    #[test]
    fn minimize_finds_a_smooth_minimum() {
        let s = minimize(|x| (x - 0.3) * (x - 0.3) + 1.0, -2.0, 2.0, C).unwrap();
        assert_relative_eq!(s.value, 0.3, epsilon = 1e-7);
        assert_relative_eq!(s.residual, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn minimize_handles_a_flat_minimum() {
        // Quartic: the gradient vanishes to third order at the minimum, so a
        // derivative-based method has almost nothing to follow.
        let s = minimize(|x: f64| (x - 0.5).powi(4), -1.0, 2.0, C).unwrap();
        assert!((s.value - 0.5).abs() < 1e-3, "landed at {}", s.value);
        assert!(s.residual < 1e-12);
    }

    #[test]
    fn minimize_refuses_a_malformed_bracket() {
        assert!(minimize(|x| x, 1.0, 0.0, C).is_err());
        assert!(minimize(|x| x, 0.0, f64::INFINITY, C).is_err());
    }

    #[test]
    fn newton_system_solves_a_two_by_two() {
        // x^2 + y^2 = 25, x - y = 1  ->  (4, 3)
        let s = newton_system(
            |v| {
                let (x, y) = (v[0], v[1]);
                (
                    vec![x.mul_add(x, y * y) - 25.0, x - y - 1.0],
                    vec![vec![2.0 * x, 2.0 * y], vec![1.0, -1.0]],
                )
            },
            &[5.0, 1.0],
            C,
        )
        .unwrap();
        assert!(s.convergence.is_converged());
        assert_relative_eq!(s.value[0], 4.0, epsilon = 1e-10);
        assert_relative_eq!(s.value[1], 3.0, epsilon = 1e-10);
    }

    #[test]
    fn newton_system_damping_survives_a_start_where_plain_newton_diverges() {
        // arctan is the classic case: an undamped Newton step from |x| > 1.4
        // overshoots to a *larger* residual, and each subsequent step is worse,
        // so the iteration runs away. From (5, 5) the first full step lands
        // near -30. Halving until the residual actually falls is what recovers
        // it.
        let s = newton_system(
            |v| {
                let (x, y) = (v[0], v[1]);
                (
                    vec![x.atan(), y.atan()],
                    vec![
                        vec![x.mul_add(x, 1.0).recip(), 0.0],
                        vec![0.0, y.mul_add(y, 1.0).recip()],
                    ],
                )
            },
            &[5.0, 5.0],
            C,
        )
        .unwrap();
        assert!(s.convergence.is_converged(), "{s:?}");
        assert!(s.value[0].abs() < 1e-9 && s.value[1].abs() < 1e-9, "{s:?}");
    }

    #[test]
    fn newton_system_reports_a_residual_minimum_rather_than_looping() {
        // No root exists: x^2 + 1 is never zero. The solver must stop at the
        // residual's minimum and say it did not converge, rather than iterate
        // to its limit or present the estimate as a solution.
        let s = newton_system(
            |v| {
                let (x, y) = (v[0], v[1]);
                (
                    vec![x.mul_add(x, 1.0), y],
                    vec![vec![2.0 * x, 0.0], vec![0.0, 1.0]],
                )
            },
            &[2.0, 2.0],
            C,
        )
        .unwrap();
        assert!(!s.convergence.is_converged());
        assert!(s.residual >= 1.0, "the residual cannot go below 1 here");
    }

    #[test]
    fn newton_system_reports_a_singular_jacobian_rather_than_looping() {
        let s = newton_system(
            |v| {
                (
                    vec![v[0] * v[0], v[1]],
                    vec![vec![2.0 * v[0], 0.0], vec![0.0, 0.0]],
                )
            },
            &[1.0, 1.0],
            C,
        );
        assert!(s.is_err());
    }

    #[test]
    fn newton_system_checks_its_shapes() {
        let s = newton_system(|_| (vec![1.0], vec![vec![1.0, 2.0]]), &[0.0, 0.0], C);
        assert!(s.is_err());
    }

    #[test]
    fn exhausted_is_reported_not_hidden() {
        // One iteration cannot possibly converge; the result must say so rather
        // than present the first guess as an answer.
        let s = brent(
            |x| x * x - 2.0,
            0.0,
            2.0,
            Criteria {
                max_iterations: 1,
                ..C
            },
        )
        .unwrap();
        assert_eq!(s.convergence, Convergence::Exhausted);
        assert!(!s.convergence.is_converged());
    }
}
