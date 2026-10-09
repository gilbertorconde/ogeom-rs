//! Root finding, polynomial roots and minimization.
//!
//! The numerical substrate every intersection, projection and extrema algorithm
//! in the kernel sits on. Nothing here is geometric. It is kept separate so
//! those algorithms are about geometry rather than about convergence.
//!
//! # What to reach for
//!
//! - A root inside a known bracket: [`brent`]. Guaranteed to converge, and
//!   nearly as fast as Newton in practice.
//! - A root with a known derivative and a good starting point: [`newton`],
//!   which falls back to bisection whenever a step would leave the bracket.
//!   Unsafeguarded Newton diverges on the configurations that matter: a
//!   tangential intersection is exactly where the derivative vanishes.
//! - Roots of a polynomial: [`roots`]. Closed form up to the cubic, the
//!   quadratic written to avoid the cancellation the schoolbook formula
//!   suffers. Above that, or with no allocation, [`real_roots`]: each root
//!   bracketed between its derivative's roots.
//! - A system of equations: [`newton_system`]. Surface projection is two
//!   equations in two unknowns. Intersection marching is much the same.
//! - A minimum without derivatives: [`minimize`].

use nalgebra::{DMatrix, DVector};
use ogeom_core::{OgeomResult, ogeom_bail};

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
/// practice, the right default when a bracket is available.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the bracket is
/// malformed, or if `f` does not change sign across it, which means there is no
/// root to find by this method.
pub fn brent<F>(mut f: F, a: f64, b: f64, criteria: Criteria) -> OgeomResult<Solution>
where
    F: FnMut(f64) -> f64,
{
    if !a.is_finite() || !b.is_finite() || a >= b {
        ogeom_bail!(Construction, "bracket [{a}, {b}] is empty or non-finite");
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
        ogeom_bail!(
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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the bracket is
/// malformed or `f` does not change sign across it.
pub fn newton<F>(mut f: F, a: f64, b: f64, start: f64, criteria: Criteria) -> OgeomResult<Solution>
where
    F: FnMut(f64) -> (f64, f64),
{
    if !a.is_finite() || !b.is_finite() || a >= b {
        ogeom_bail!(Construction, "bracket [{a}, {b}] is empty or non-finite");
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
        ogeom_bail!(Construction, "f does not change sign across [{a}, {b}]");
    }
    // Orient so that f(low) < 0 < f(high). The bracket update is then a single
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
/// Degrees up to three are solved in closed form. Above that each root is
/// bracketed between the derivative's roots, as [`real_roots`] does.
///
/// Repeated roots are returned once each, since a geometry caller wants the
/// distinct parameter values, and a double root (a tangency) is found
/// although rounding leaves the polynomial a hair either side of zero there.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if every
/// coefficient is zero, where every value is a root.
pub fn roots(coefficients: &[f64], tolerance: f64) -> OgeomResult<Vec<f64>> {
    // Drop leading zeros so the true degree drives the choice of method. A
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
        0 => ogeom_bail!(
            Construction,
            "the zero polynomial has every value as a root"
        ),
        1 => Vec::new(),
        2 => vec![-c[0] / c[1]],
        3 => quadratic_roots(c[2], c[1], c[0]),
        4 => cubic_roots(c[3], c[2], c[1], c[0]),
        n => {
            let mut found = vec![0.0; n - 1];
            let count = real_roots(c, 0.0, &mut found);
            found.truncate(count);
            found
        }
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
/// precision for the smaller root when `b^2 >> 4ac`, which is the common case
/// for a ray grazing a sphere.
#[must_use]
pub fn quadratic_roots(a: f64, b: f64, c: f64) -> Vec<f64> {
    if a == 0.0 {
        return if b == 0.0 { Vec::new() } else { vec![-c / b] };
    }
    let discriminant = b.mul_add(b, -(4.0 * a * c));
    // A double root's discriminant is zero only up to the rounding of the
    // coefficients, which can leave it a hair below.
    let rounding = 8.0 * f64::EPSILON * b.mul_add(b, (4.0 * a * c).abs());
    if discriminant.abs() <= rounding {
        return vec![-b / (2.0 * a)];
    }
    if discriminant < 0.0 {
        return Vec::new();
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
    // What rounding leaves in `p`, `q` and so in the discriminant: a double
    // root's is zero only to within it.
    let p_rounding = 8.0 * f64::EPSILON * ((shift * b).abs() + c.abs());
    let q_rounding =
        8.0 * f64::EPSILON * ((2.0 / 27.0 * b * b * b).abs() + (shift * c).abs() + d.abs());
    let rounding = half_q.abs() * q_rounding
        + third_p * third_p * p_rounding
        + 8.0 * f64::EPSILON * (half_q * half_q + (third_p * third_p * third_p).abs());

    if discriminant.abs() <= rounding {
        if p.abs() <= p_rounding {
            vec![-shift]
        } else {
            // t^3 + p t + q with a double root: 3q/p once, -3q/(2p) twice.
            let mut r = vec![3.0 * q / p - shift, -1.5 * q / p - shift];
            r.sort_by(|x, y| x.partial_cmp(y).unwrap_or(core::cmp::Ordering::Equal));
            r
        }
    } else if discriminant > 0.0 {
        let sqrt = discriminant.sqrt();
        let u = (-half_q + sqrt).cbrt();
        let v = (-half_q - sqrt).cbrt();
        vec![u + v - shift]
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

/// Polynomials up to this degree are solved with their working on the stack.
const STACK_DEGREE: usize = 16;

/// The real roots of a polynomial, in increasing order, written to the front
/// of `out`; returns how many there are.
///
/// `coefficients` are in ascending power order, as for [`roots`]; leading
/// coefficients that are exactly zero are dropped. Nothing is allocated up
/// to degree 16.
///
/// The derivative's real roots, found the same way down to a line, split
/// the real line into intervals on which the polynomial is monotone. Each
/// interval whose ends differ in sign holds one root, which safeguarded
/// Newton finds with bisection to fall back on (Yuksel, "High-Performance
/// Polynomial Root Finding for Graphics", HPG 2022). A root of even
/// multiplicity, a tangency, has no sign change around it: it is a critical
/// point where the polynomial's value is zero to rounding, or no larger
/// than `touch` with no root on either side. `touch` is in the polynomial's
/// own units, since only the caller knows what rounding its coefficients
/// carry. Such a root comes back once.
///
/// A polynomial whose every coefficient is zero, or one that is not finite,
/// has no roots returned.
///
/// # Panics
///
/// If `out` has room for fewer values than the polynomial's degree.
pub fn real_roots(coefficients: &[f64], touch: f64, out: &mut [f64]) -> usize {
    let mut c = coefficients;
    while let Some((&0.0, rest)) = c.split_last() {
        c = rest;
    }
    if c.len() < 2 || c.iter().any(|v| !v.is_finite()) {
        return 0;
    }
    let n = c.len() - 1;
    assert!(out.len() >= n, "room for {} of {n} roots", out.len());
    // Fujiwara's bound on the roots' moduli, doubled so that no root stands
    // near an end where rounding could flip the value's sign. By Gauss-Lucas
    // every derivative's roots lie within it too.
    let lead = c[n];
    let mut bound = 0.0_f64;
    for k in 1..=n {
        let ratio = (c[n - k] / lead).abs();
        let ratio = if k == n { ratio / 2.0 } else { ratio };
        #[allow(clippy::cast_precision_loss, reason = "a degree")]
        let term = if k == 1 {
            ratio
        } else {
            ratio.powf(1.0 / k as f64)
        };
        bound = bound.max(term);
    }
    let bound = if bound > 0.0 { 4.0 * bound } else { 1.0 };
    if !bound.is_finite() {
        return 0;
    }
    monotone_roots(c, touch, bound, out)
}

/// The roots of `c` (degree at least one, nonzero lead) within `bound`.
fn monotone_roots(c: &[f64], touch: f64, bound: f64, out: &mut [f64]) -> usize {
    let n = c.len() - 1;
    if n == 1 {
        out[0] = -c[0] / c[1];
        return 1;
    }
    if n <= STACK_DEGREE {
        let mut slope = [0.0; STACK_DEGREE];
        let mut critical = [0.0; STACK_DEGREE];
        split_roots(
            c,
            touch,
            bound,
            out,
            &mut slope[..n],
            &mut critical[..n - 1],
        )
    } else {
        let (mut slope, mut critical) = (vec![0.0; n], vec![0.0; n - 1]);
        split_roots(c, touch, bound, out, &mut slope, &mut critical)
    }
}

/// The roots of `c`, from the roots of its derivative, with `slope` and
/// `critical` as working space for the derivative and its roots.
fn split_roots(
    c: &[f64],
    touch: f64,
    bound: f64,
    out: &mut [f64],
    slope: &mut [f64],
    critical: &mut [f64],
) -> usize {
    let n = c.len() - 1;
    for (k, s) in slope.iter_mut().enumerate() {
        #[allow(clippy::cast_precision_loss, reason = "a degree")]
        let power = (k + 1) as f64;
        *s = power * c[k + 1];
    }
    let found = monotone_roots(slope, 0.0, bound, critical);
    // What Horner's rule and the coefficients' own rounding can leave in a
    // value, relative to its terms' size, with room to spare: a value no
    // larger is zero.
    #[allow(clippy::cast_precision_loss, reason = "a degree")]
    let rounding = 32.0 * n as f64 * f64::EPSILON;

    let mut count = 0;
    let mut push = |x: f64| {
        if count < out.len() {
            out[count] = x;
            count += 1;
        }
    };
    // Walk the monotone intervals left to right. A critical point within
    // `touch` of zero is a root only once the intervals either side of it
    // are known not to cross: if they do, its value is on the far side of
    // zero and the crossings are the roots.
    let mut a = -bound;
    let (mut fa, _, _) = value_and_slope(c, a);
    let mut a_touches = false;
    let mut crossed_before = false;
    let ends = critical[..found].iter().map(|&x| (x, true));
    for (b, is_critical) in ends.chain(core::iter::once((bound, false))) {
        if b <= a {
            continue;
        }
        let (mut fb, _, size) = value_and_slope(c, b);
        let b_touches = is_critical && fb.abs() <= touch;
        if is_critical && fb.abs() <= rounding * size {
            fb = 0.0;
        }
        let crosses = fa * fb < 0.0;
        if fa == 0.0 || (a_touches && !crossed_before && !crosses) {
            push(a);
        }
        if crosses {
            push(monotone_root(c, a, b, fa, fb));
        }
        crossed_before = crosses;
        a_touches = b_touches;
        (a, fa) = (b, fb);
    }
    if fa == 0.0 {
        push(a);
    }
    count
}

/// A polynomial's value and slope at `x`, and the size of its terms there,
/// `sum |c_k| |x|^k`, which bounds the rounding in the value.
fn value_and_slope(c: &[f64], x: f64) -> (f64, f64, f64) {
    let (mut p, mut dp, mut size) = (0.0_f64, 0.0_f64, 0.0_f64);
    let magnitude = x.abs();
    for &coefficient in c.iter().rev() {
        dp = dp.mul_add(x, p);
        p = p.mul_add(x, coefficient);
        size = size.mul_add(magnitude, coefficient.abs());
    }
    (p, dp, size)
}

/// The one root of `c` in `[a, b]`, where `fa` and `fb` differ in sign:
/// Newton from the secant's point, bisecting whenever a step would leave
/// the bracket or does not halve the step before last. The bracket keeps
/// the sign change, so the result is a root even where `c` is not monotone.
fn monotone_root(c: &[f64], mut a: f64, mut b: f64, fa: f64, fb: f64) -> f64 {
    let rising = fb > 0.0;
    let mut x = a - fa * (b - a) / (fb - fa);
    if !(x > a && x < b) {
        x = f64::midpoint(a, b);
    }
    let (mut step, mut previous) = (b - a, b - a);
    for _ in 0..128 {
        let (f, df, _) = value_and_slope(c, x);
        if f == 0.0 {
            return x;
        }
        if (f > 0.0) == rising {
            b = x;
        } else {
            a = x;
        }
        let delta = f / df;
        let newton = x - delta;
        let next = if newton > a && newton < b && delta.abs() * 2.0 <= previous.abs() {
            previous = step;
            step = delta;
            newton
        } else {
            previous = step;
            step = 0.5 * (b - a);
            f64::midpoint(a, b)
        };
        // A step below rounding has converged; so has a bracket with no
        // float left inside it.
        if next == x || step.abs() <= 2.0 * f64::EPSILON * next.abs() || next == a || next == b {
            return next;
        }
        x = next;
    }
    x
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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the bracket is
/// malformed.
pub fn minimize<F>(mut f: F, a: f64, b: f64, criteria: Criteria) -> OgeomResult<Solution>
where
    F: FnMut(f64) -> f64,
{
    if !a.is_finite() || !b.is_finite() || a >= b {
        ogeom_bail!(Construction, "bracket [{a}, {b}] is empty or non-finite");
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
/// halved until it actually reduces the residual. Undamped Newton overshoots
/// badly from a poor start, and a geometry caller's start is often only a rough
/// guess from a coarse sampling.
///
/// Surface projection is this with two equations in two unknowns, and so is a
/// step of a surface/surface intersection march.
///
/// Where no root exists the residual has a positive minimum, and no damping
/// finds a downhill step from it. That is reported as
/// [`Convergence::Exhausted`] with the best estimate attached. "No root here"
/// is a useful answer, and far better than iterating to the limit.
///
/// # Errors
///
/// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) if the Jacobian's shape
/// disagrees with the residual, and
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if the Jacobian is singular
/// at the starting point.
pub fn newton_system<F>(mut f: F, start: &[f64], criteria: Criteria) -> OgeomResult<SystemSolution>
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
        ogeom_bail!(
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
            ogeom_bail!(Numeric, "Jacobian is singular after {iteration} iterations");
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

/// A fixed-size [`newton_system`] for `N` unknowns, allocation-free.
///
/// The same damped iteration (halving until the residual falls, the same
/// three verdicts), on stack arrays and a fixed-size LU: the intersectors'
/// three- and four-unknown systems run this millions of times per model,
/// and the general path allocates for every residual, Jacobian and step.
///
/// # Errors
///
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if the
/// Jacobian is singular.
pub fn newton_system_fixed<const N: usize, F>(
    mut f: F,
    start: [f64; N],
    criteria: Criteria,
) -> OgeomResult<([f64; N], f64, Convergence, usize)>
where
    F: FnMut(&[f64; N]) -> ([f64; N], [[f64; N]; N]),
{
    // The Jacobian of the point last measured: the lazy iteration asks for
    // it only right after measuring the same point.
    let last = std::cell::Cell::new([[0.0; N]; N]);
    newton_system_fixed_lazy(
        |x| {
            let (residual, jacobian) = f(x);
            last.set(jacobian);
            residual
        },
        |_| Some(last.get()),
        start,
        criteria,
    )
}

/// [`newton_system_fixed`] with the residual and the Jacobian measured
/// apart, the Jacobian only where a step is taken.
///
/// The damped search measures every trial point and keeps few, and where
/// the Jacobian costs more than the residual (a spline surface's
/// derivatives against its point) most of that goes on Jacobians never
/// used. The iteration is the one [`newton_system_fixed`] runs on
/// `|x| (residual(x), jacobian(x))`, a `None` Jacobian standing for a
/// point that cannot be measured: an infinite residual and a zero
/// Jacobian. So it answers the same to the bit.
///
/// # Errors
///
/// As [`newton_system_fixed`].
pub fn newton_system_fixed_lazy<const N: usize, R, J>(
    mut residual_at: R,
    mut jacobian_at: J,
    start: [f64; N],
    criteria: Criteria,
) -> OgeomResult<([f64; N], f64, Convergence, usize)>
where
    R: FnMut(&[f64; N]) -> [f64; N],
    J: FnMut(&[f64; N]) -> Option<[[f64; N]; N]>,
{
    let norm_of = |r: &[f64; N]| r.iter().map(|v| v * v).sum::<f64>().sqrt();
    let mut x = start;
    let first = residual_at(&x);
    let (mut residual, mut jacobian) =
        jacobian_at(&x).map_or(([f64::INFINITY; N], [[0.0; N]; N]), |j| (first, j));
    let mut norm = norm_of(&residual);
    for iteration in 1..=criteria.max_iterations {
        if norm <= criteria.residual {
            return Ok((x, norm, Convergence::Residual, iteration - 1));
        }
        let Some(delta) = solve_fixed(jacobian, residual) else {
            ogeom_bail!(Numeric, "Jacobian is singular after {iteration} iterations");
        };
        let mut scale = 1.0;
        let mut accepted = None;
        for _ in 0..30 {
            let mut candidate = x;
            for (value, d) in candidate.iter_mut().zip(delta.iter()) {
                *value -= d * scale;
            }
            let r = residual_at(&candidate);
            let candidate_norm = norm_of(&r);
            // A point whose Jacobian cannot be measured is unmeasured, and
            // no step is accepted onto it.
            if (candidate_norm < norm || candidate_norm <= criteria.residual)
                && let Some(jj) = jacobian_at(&candidate)
            {
                accepted = Some((candidate, r, jj, candidate_norm));
                break;
            }
            scale *= 0.5;
        }
        let Some((next, r, jj, next_norm)) = accepted else {
            return Ok((x, norm, Convergence::Exhausted, iteration));
        };
        let step = next
            .iter()
            .zip(&x)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt();
        x = next;
        residual = r;
        jacobian = jj;
        norm = next_norm;
        if norm <= criteria.residual {
            return Ok((x, norm, Convergence::Residual, iteration));
        }
        if step <= criteria.step {
            return Ok((x, norm, Convergence::Step, iteration));
        }
    }
    Ok((x, norm, Convergence::Exhausted, criteria.max_iterations))
}

/// `A x = b` by LU with partial pivoting, `None` at a zero pivot.
///
/// The arithmetic is the general solver's LU step for step: the first
/// largest magnitude pivots, the multipliers are the column times the
/// pivot's reciprocal, and the back substitution runs column by column from
/// the last. So a system moved from [`newton_system`] to
/// [`newton_system_fixed`] solves to the same bits.
fn solve_fixed<const N: usize>(mut a: [[f64; N]; N], mut b: [f64; N]) -> Option<[f64; N]> {
    for col in 0..N {
        let mut pivot = col;
        for row in col + 1..N {
            if a[row][col].abs() > a[pivot][col].abs() {
                pivot = row;
            }
        }
        let diag = a[pivot][col];
        if diag == 0.0 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        let inverse = 1.0 / diag;
        let head = a[col];
        for row in col + 1..N {
            let factor = a[row][col] * inverse;
            for (entry, above) in a[row].iter_mut().zip(&head).skip(col + 1) {
                *entry -= factor * above;
            }
            b[row] -= factor * b[col];
        }
    }
    for col in (0..N).rev() {
        b[col] /= a[col][col];
        let known = b[col];
        for row in 0..col {
            b[row] -= a[row][col] * known;
        }
    }
    Some(b)
}

/// A two-unknown [`newton_system`], allocation-free.
///
/// The foot-point projection runs this system millions of times per real
/// model, and the general path pays a heap allocation for every residual,
/// Jacobian, vector and factorization of every damped step. The algorithm
/// here is the same (damped Newton, halving until the residual falls, the
/// same three convergence verdicts), with the two-by-two solve written out:
/// partial pivoting is one comparison, and singularity is a vanishing
/// pivot.
///
/// # Errors
///
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if the Jacobian
/// is singular.
pub fn newton_system_2<F>(
    mut f: F,
    start: [f64; 2],
    criteria: Criteria,
) -> OgeomResult<([f64; 2], f64, Convergence, usize)>
where
    F: FnMut([f64; 2]) -> ([f64; 2], [[f64; 2]; 2]),
{
    let mut x = start;
    let (mut residual, mut jacobian) = f(x);
    let mut norm = residual[0].hypot(residual[1]);

    for iteration in 1..=criteria.max_iterations {
        if norm <= criteria.residual {
            return Ok((x, norm, Convergence::Residual, iteration - 1));
        }

        // Solve J * delta = residual, partial pivoting on the first column.
        let (row0, row1, rhs0, rhs1) = if jacobian[0][0].abs() >= jacobian[1][0].abs() {
            (jacobian[0], jacobian[1], residual[0], residual[1])
        } else {
            (jacobian[1], jacobian[0], residual[1], residual[0])
        };
        if row0[0].abs() <= f64::EPSILON * (row1[0].abs() + row0[1].abs()).max(1.0) {
            ogeom_bail!(Numeric, "Jacobian is singular after {iteration} iterations");
        }
        let factor = row1[0] / row0[0];
        let denom = factor.mul_add(-row0[1], row1[1]);
        if denom.abs() <= f64::EPSILON * row0[1].abs().max(1.0) {
            ogeom_bail!(Numeric, "Jacobian is singular after {iteration} iterations");
        }
        let d1 = factor.mul_add(-rhs0, rhs1) / denom;
        let d0 = d1.mul_add(-row0[1], rhs0) / row0[0];
        let delta = [d0, d1];

        // Damping: keep halving until the residual actually falls.
        let mut scale = 1.0;
        let mut accepted = None;
        for _ in 0..30 {
            let candidate = [
                delta[0].mul_add(-scale, x[0]),
                delta[1].mul_add(-scale, x[1]),
            ];
            let (r, jj) = f(candidate);
            let candidate_norm = r[0].hypot(r[1]);
            if candidate_norm < norm || candidate_norm <= criteria.residual {
                accepted = Some((candidate, r, jj, candidate_norm));
                break;
            }
            scale *= 0.5;
        }
        let Some((next, r, jj, next_norm)) = accepted else {
            return Ok((x, norm, Convergence::Exhausted, iteration));
        };

        let step = (next[0] - x[0]).hypot(next[1] - x[1]);
        x = next;
        residual = r;
        jacobian = jj;
        norm = next_norm;

        if norm <= criteria.residual {
            return Ok((x, norm, Convergence::Residual, iteration));
        }
        if step <= criteria.step {
            return Ok((x, norm, Convergence::Step, iteration));
        }
    }
    Ok((x, norm, Convergence::Exhausted, criteria.max_iterations))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    /// Every value in `expected` is among `found`, to the precision a double
    /// root keeps (half the digits).
    fn has_each(found: &[f64], expected: &[f64]) {
        for e in expected {
            assert!(
                found
                    .iter()
                    .any(|f| (f - e).abs() <= 1e-6 * e.abs().max(1.0)),
                "{e} missing from {found:?}"
            );
        }
    }

    /// A quartic whose roots pair off by sign (a ray through a torus's
    /// middle) comes back, where the Schur iteration alone ran on for ever.
    #[test]
    fn a_quartic_with_roots_paired_by_sign_comes_back() {
        // (t^2 - 4)(t^2 - 16) = t^4 - 20 t^2 + 64.
        let found = roots(&[64.0, 0.0, -20.0, 0.0, 1.0], 1e-12).unwrap();
        assert_eq!(found.len(), 4, "{found:?}");
        for (got, want) in found.iter().zip([-4.0, -2.0, 2.0, 4.0]) {
            assert!((got - want).abs() < 1e-9, "{found:?}");
        }
        // And one with no real root at all.
        assert!(
            roots(&[64.0, 0.0, 4.0, 0.0, 1.0], 1e-12)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn double_roots_survive_rounding() {
        for r in [0.1, 0.3, 0.7, 1.1, 3.3, 1000.0 / 3.0, -2.9, 1e-3 / 7.0] {
            // (x - r)^2
            has_each(&quadratic_roots(1.0, -2.0 * r, r * r), &[r]);
            // (x - 1)(x - r)^2 = x^3 - (1 + 2r) x^2 + (2r + r^2) x - r^2
            has_each(
                &cubic_roots(1.0, -(1.0 + 2.0 * r), 2.0f64.mul_add(r, r * r), -(r * r)),
                &[1.0, r],
            );
            // (x - 2)(x + 1)(x - r)^2, through the general solver.
            let quadratic = [r * r, -2.0 * r, 1.0];
            let pair = [-2.0, -1.0, 1.0];
            let mut quartic = [0.0; 5];
            for (i, a) in quadratic.iter().enumerate() {
                for (j, b) in pair.iter().enumerate() {
                    quartic[i + j] += a * b;
                }
            }
            has_each(&roots(&quartic, 1e-12).unwrap(), &[-1.0, 2.0, r]);
        }
    }
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
        // f(x) = x^3 has f'(0) = 0. Unsafeguarded Newton stalls. The bisection
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
        // One iteration cannot possibly converge. The result must say so rather
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

    /// The fixed-size solver lands on the general one's bits: the same
    /// iterates, verdict and count, on nonlinear systems of three, four and
    /// five unknowns whose Jacobians need pivoting and have ties for it.
    #[test]
    fn the_fixed_solver_matches_the_general_one_to_the_bit() {
        fn check<const N: usize>(start: [f64; N]) {
            let system = |x: &[f64]| {
                let mut r = vec![0.0; N];
                let mut j = vec![vec![0.0; N]; N];
                for i in 0..N {
                    let k = (i + 1) % N;
                    #[allow(clippy::cast_precision_loss)]
                    let weight = 1.0 + i as f64 * 0.37;
                    r[i] = x[i].mul_add(x[k], -weight) + x[k].sin() * 0.3;
                    j[i][i] += x[k];
                    j[i][k] += x[i] + x[k].cos() * 0.3;
                }
                (r, j)
            };
            let criteria = Criteria {
                residual: 1e-14,
                step: 1e-15,
                max_iterations: 60,
            };
            let general = newton_system(system, &start, criteria).unwrap();
            let fixed = newton_system_fixed(
                |x: &[f64; N]| {
                    let (r, j) = system(x);
                    let mut rows = [[0.0; N]; N];
                    for (to, from) in rows.iter_mut().zip(&j) {
                        to.copy_from_slice(from);
                    }
                    (r.try_into().unwrap(), rows)
                },
                start,
                criteria,
            )
            .unwrap();
            assert_eq!(general.value, fixed.0.to_vec());
            assert_eq!(general.residual.to_bits(), fixed.1.to_bits());
            assert_eq!(general.convergence, fixed.2);
            assert_eq!(general.iterations, fixed.3);
        }
        check([1.0, 1.0, 1.0]);
        check([0.5, 2.0, -1.0, 1.5]);
        check([1.0, -1.0, 1.0, -1.0, 2.0]);
        check([3.0, 0.2, 0.7, 1.1, 0.4]);
    }

    /// The bracketing solver against two oracles: roots built from known
    /// factors, and the companion matrix's real eigenvalues.
    mod bracketing {
        use super::*;
        use proptest::prelude::*;

        /// The real eigenvalues of the companion matrix, with near-real
        /// pairs kept where the polynomial all but vanishes, and Durand-Kerner
        /// where the Schur iteration does not settle.
        fn companion(c: &[f64]) -> Vec<f64> {
            let n = c.len() - 1;
            let lead = c[n];
            let mut m = DMatrix::<f64>::zeros(n, n);
            for i in 0..n {
                m[(i, n - 1)] = -c[i] / lead;
                if i + 1 < n {
                    m[(i + 1, i)] = 1.0;
                }
            }
            let eigenvalues: Vec<nalgebra::Complex<f64>> =
                match nalgebra::linalg::Schur::try_new(m, f64::EPSILON, 1000) {
                    Some(schur) => schur.complex_eigenvalues().iter().copied().collect(),
                    None => durand_kerner(c),
                };
            let mut out: Vec<f64> = eigenvalues
                .iter()
                .filter_map(|e| {
                    let scale = e.re.abs().max(1.0);
                    if e.im.abs() <= 1e-9 * scale {
                        return Some(e.re);
                    }
                    if e.im.abs() > 1e-7 * scale {
                        return None;
                    }
                    let (p, _, size) = value_and_slope(c, e.re);
                    (p.abs() <= 1e-10 * size).then_some(e.re)
                })
                .collect();
            out.sort_by(f64::total_cmp);
            out
        }

        fn durand_kerner(c: &[f64]) -> Vec<nalgebra::Complex<f64>> {
            use nalgebra::Complex;
            let n = c.len() - 1;
            let monic: Vec<f64> = c.iter().map(|x| x / c[n]).collect();
            let radius = 1.0 + monic[..n].iter().fold(0.0_f64, |m, x| m.max(x.abs()));
            #[allow(clippy::cast_precision_loss)]
            let mut z: Vec<Complex<f64>> = (0..n)
                .map(|k| {
                    Complex::from_polar(radius, 0.4 + core::f64::consts::TAU * k as f64 / n as f64)
                })
                .collect();
            let value = |x: Complex<f64>| {
                let mut p = Complex::new(1.0, 0.0);
                for &coefficient in monic[..n].iter().rev() {
                    p = p * x + coefficient;
                }
                p
            };
            for _ in 0..500 {
                let mut largest = 0.0_f64;
                for i in 0..n {
                    let mut denominator = Complex::new(1.0, 0.0);
                    for j in 0..n {
                        if i != j {
                            denominator *= z[i] - z[j];
                        }
                    }
                    if denominator.norm() == 0.0 {
                        continue;
                    }
                    let step = value(z[i]) / denominator;
                    z[i] -= step;
                    largest = largest.max(step.norm() / z[i].norm().max(1.0));
                }
                if largest <= f64::EPSILON * 4.0 {
                    break;
                }
            }
            z
        }

        /// The ascending coefficients of `lead * prod (x - r)`.
        fn from_roots(lead: f64, roots: &[f64]) -> Vec<f64> {
            let mut c = vec![lead];
            for &r in roots {
                let mut next = vec![0.0; c.len() + 1];
                for (k, a) in c.iter().enumerate() {
                    next[k + 1] += a;
                    next[k] -= a * r;
                }
                c = next;
            }
            c
        }

        /// `from_roots` times `x^2 - 2 re x + re^2 + im^2`, a pair with no
        /// real root.
        fn with_pair(c: &[f64], re: f64, im: f64) -> Vec<f64> {
            let quadratic = [re.mul_add(re, im * im), -2.0 * re, 1.0];
            let mut out = vec![0.0; c.len() + 2];
            for (i, a) in c.iter().enumerate() {
                for (j, b) in quadratic.iter().enumerate() {
                    out[i + j] += a * b;
                }
            }
            out
        }

        /// The roots of a line-torus quartic, touching within the
        /// intersection's threshold.
        fn on_torus(c: &[f64; 5]) -> Vec<f64> {
            solve(c, 1e-9 * c.iter().map(|x| x.abs()).sum::<f64>())
        }

        fn solve(c: &[f64], touch: f64) -> Vec<f64> {
            let mut out = vec![0.0; c.len()];
            let n = real_roots(c, touch, &mut out);
            out.truncate(n);
            out
        }

        /// `found` is `expected` one for one, in order, each within `relative`
        /// of the root's magnitude (and of `scale` near zero).
        fn matches(found: &[f64], expected: &[f64], relative: f64, scale: f64) -> bool {
            found.len() == expected.len()
                && found
                    .iter()
                    .zip(expected)
                    .all(|(f, e)| (f - e).abs() <= relative * e.abs().max(scale))
        }

        fn sorted(mut v: Vec<f64>) -> Vec<f64> {
            v.sort_by(f64::total_cmp);
            v
        }

        /// The line-torus quartic, in the frame and units the intersection
        /// solves it in: offset `m` from the centre at the line's nearest
        /// approach, unit direction `d`, radii over their sum.
        fn torus_quartic(m: [f64; 3], d: [f64; 3], big: f64) -> [f64; 5] {
            let small = 1.0 - big;
            let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
            let a = dot(d, d);
            let b = 2.0 * dot(m, d);
            let c = dot(m, m) + big * big - small * small;
            let p = d[0] * d[0] + d[1] * d[1];
            let q = 2.0 * (m[0] * d[0] + m[1] * d[1]);
            let s = m[0] * m[0] + m[1] * m[1];
            let four = 4.0 * big * big;
            [
                c * c - four * s,
                2.0 * b * c - four * q,
                b * b + 2.0 * a * c - four * p,
                2.0 * a * b,
                a * a,
            ]
        }

        #[test]
        fn simple_and_repeated_roots_come_back_once_each() {
            // Four simple roots, both signs.
            let c = from_roots(1.0, &[-3.0, -0.5, 1.0, 7.0]);
            assert!(matches(
                &solve(&c, 0.0),
                &[-3.0, -0.5, 1.0, 7.0],
                1e-12,
                1.0
            ));
            // A double root between two simple ones, and at either end of
            // the real roots.
            for (roots, distinct) in [
                ([-1.0, 2.0, 2.0, 5.0], [-1.0, 2.0, 5.0]),
                ([-4.0, -4.0, 0.5, 3.0], [-4.0, 0.5, 3.0]),
                ([-2.0, 1.5, 6.0, 6.0], [-2.0, 1.5, 6.0]),
            ] {
                let c = from_roots(2.0, &roots);
                let found = solve(&c, 0.0);
                assert!(
                    matches(&found, &distinct, 1e-7, 1.0),
                    "{roots:?}: {found:?}"
                );
            }
            // Two double roots, and a triple and a quadruple root.
            let found = solve(&from_roots(1.0, &[-1.0, -1.0, 3.0, 3.0]), 0.0);
            assert!(matches(&found, &[-1.0, 3.0], 1e-7, 1.0), "{found:?}");
            let found = solve(&from_roots(1.0, &[0.5, 2.0, 2.0, 2.0]), 0.0);
            assert!(matches(&found, &[0.5, 2.0], 1e-5, 1.0), "{found:?}");
            let found = solve(&from_roots(1.0, &[1.25; 4]), 0.0);
            assert!(matches(&found, &[1.25], 1e-3, 1.0), "{found:?}");
        }

        #[test]
        fn a_root_at_zero_or_at_a_critical_point_is_found_once() {
            // x (x - 1)(x + 2)(x - 3): a root exactly at zero.
            let found = solve(&from_roots(1.0, &[0.0, 1.0, -2.0, 3.0]), 0.0);
            assert!(
                matches(&found, &[-2.0, 0.0, 1.0, 3.0], 1e-12, 1.0),
                "{found:?}"
            );
            // x^4 - 2 x^2 = x^2 (x^2 - 2): the double root at zero is a
            // critical point whose value is exactly zero.
            let found = solve(&[0.0, 0.0, -2.0, 0.0, 1.0], 0.0);
            let r = 2.0_f64.sqrt();
            assert!(matches(&found, &[-r, 0.0, r], 1e-12, 1.0), "{found:?}");
            // x^4 alone, and x^4 - 1 with roots at the bound's scale.
            assert_eq!(solve(&[0.0, 0.0, 0.0, 0.0, 1.0], 0.0), vec![0.0]);
            let found = solve(&[-1.0, 0.0, 0.0, 0.0, 1.0], 0.0);
            assert!(matches(&found, &[-1.0, 1.0], 1e-14, 1.0), "{found:?}");
        }

        #[test]
        fn degenerate_input_has_no_roots_and_high_degree_works() {
            assert_eq!(real_roots(&[0.0; 5], 0.0, &mut [0.0; 4]), 0);
            assert_eq!(real_roots(&[1.0, f64::NAN, 1.0], 0.0, &mut [0.0; 2]), 0);
            assert_eq!(real_roots(&[3.0, 0.0, 0.0], 0.0, &mut []), 0);
            // Leading exact zeros drop the degree.
            let found = solve(&[-4.0, 0.0, 1.0, 0.0, 0.0], 0.0);
            assert!(matches(&found, &[-2.0, 2.0], 1e-14, 1.0), "{found:?}");
            // Degree twenty, past the stack's working space.
            let roots: Vec<f64> = (0..20).map(|k| f64::from(k) * 0.25 - 2.0).collect();
            let found = solve(&from_roots(1.0, &roots), 0.0);
            assert!(matches(&found, &roots, 1e-6, 1.0), "{found:?}");
        }

        #[test]
        fn a_near_tangency_is_a_touch_only_within_the_threshold() {
            // ((x - 1)^2 + e)((x + 3)^2 + 1): no real root, the minimum 17 e
            // above zero at one. Within rounding it is a root whatever the
            // threshold; past rounding, only within the threshold.
            for (e, touch, touches) in [
                (1e-14, 0.0, true),
                (1e-11, 0.0, false),
                (1e-11, 1e-9, true),
                (1e-8, 1e-9, false),
            ] {
                let c = [1.0 + e, -2.0, 1.0];
                let c = with_pair(&c, -3.0, 1.0);
                let found = solve(&c, touch);
                if touches {
                    assert!(matches(&found, &[1.0], 1e-9, 1.0), "{e}: {found:?}");
                } else {
                    assert!(found.is_empty(), "{e}: {found:?}");
                }
            }
            // ((x - 1)^2 - e)((x + 3)^2 + 1): two real roots, or one where
            // rounding cannot tell them apart, never a third.
            for e in [1e-8, 1e-11, 1e-14] {
                let c = with_pair(&[1.0 - e, -2.0, 1.0], -3.0, 1.0);
                let found = solve(&c, 1e-9);
                let (low, high) = (1.0 - e.sqrt(), 1.0 + e.sqrt());
                assert!(
                    matches(&found, &[low, high], 1e-8, 1.0) || matches(&found, &[1.0], 1e-6, 1.0),
                    "{e}: {found:?}"
                );
            }
        }

        fn real_root() -> impl Strategy<Value = f64> {
            -1.0..1.0f64
        }

        proptest! {
            /// Distinct real roots at any scale, with or without a complex
            /// pair beside them, are all found to within rounding of their
            /// conditioning.
            #[test]
            fn known_simple_roots_are_found_at_any_scale(
                first in real_root(),
                gaps in proptest::collection::vec(0.05..0.7f64, 1..=3),
                exponent in -6i32..=6,
                lead in prop_oneof![-1e3..-1e-3f64, 1e-3..1e3f64],
                pair in proptest::option::of((real_root(), 0.1..1.0f64)),
            ) {
                let mut roots = vec![first];
                for gap in gaps {
                    roots.push(roots[roots.len() - 1] + gap);
                }
                let scale = 10f64.powi(exponent);
                let roots: Vec<f64> = roots.iter().map(|r| r * scale).collect();
                let mut c = from_roots(lead, &roots);
                if roots.len() <= 2 && let Some((re, im)) = pair {
                    c = with_pair(&c, re * scale, im * scale);
                }
                let found = solve(&c, 0.0);
                prop_assert!(matches(&found, &roots, 1e-9, scale), "{roots:?}: {found:?}");
            }

            /// A double root among simple ones comes back once, to half the
            /// digits; a pair a little apart comes back as two.
            #[test]
            fn double_and_near_double_roots(
                double in real_root(),
                apart in (0.1..1.0f64, 0.1..1.0f64),
                sides in (any::<bool>(), any::<bool>()),
                split in prop_oneof![Just(0.0), 1e-4..1e-2f64],
                exponent in -4i32..=4,
            ) {
                let scale = 10f64.powi(exponent);
                // One root either side, or both on one side a step apart.
                let side = |left: bool, d: f64| if left { double - d } else { double + d };
                let second = if sides.0 == sides.1 { apart.0 + apart.1 } else { apart.1 };
                let others = [side(sides.0, apart.0), side(sides.1, second)];
                let pair = [double - split, double + split];
                let all = [pair[0], pair[1], others[0], others[1]].map(|r| r * scale);
                let c = from_roots(1.0, &all);
                let found = solve(&c, 0.0);
                let expected = if split == 0.0 {
                    sorted(vec![double, others[0], others[1]])
                } else {
                    sorted(vec![pair[0], pair[1], others[0], others[1]])
                };
                let expected: Vec<f64> = expected.iter().map(|r| r * scale).collect();
                let precision = if split == 0.0 { 1e-6 } else { 1e-7 / split };
                prop_assert!(
                    matches(&found, &expected, precision, scale),
                    "{expected:?}: {found:?}"
                );
            }

            /// Any quartic: every root the companion matrix finds where the
            /// quartic crosses cleanly is found here, and nothing else.
            #[test]
            fn agrees_with_the_companion_matrix(
                mut c in proptest::collection::vec(-10.0..10.0f64, 4),
                lead in prop_oneof![-10.0..-0.1f64, 0.1..10.0f64],
            ) {
                c.push(lead);
                let new = solve(&c, 1e-10);
                let old = companion(&c);
                let clean = |r: &f64| {
                    let (_, slope, size) = value_and_slope(&c, *r);
                    slope.abs() * (1.0 + r.abs()) > 1e-3 * size
                };
                let near = |set: &[f64], r: f64| {
                    set.iter().any(|s| (s - r).abs() <= 1e-8 * r.abs().max(1.0))
                };
                for r in old.iter().filter(|r| clean(r)) {
                    prop_assert!(near(&new, *r), "{r} of {old:?} missing from {new:?}");
                }
                for r in new.iter().filter(|r| clean(r)) {
                    prop_assert!(near(&old, *r), "{r} of {new:?} not in {old:?}");
                }
            }

            /// Line-torus quartics: every crossing the companion matrix finds
            /// is found here, every root found makes the quartic vanish, and a
            /// tangency comes back once.
            #[test]
            fn torus_quartics_agree_with_the_companion_matrix(
                offset in (-1.5..1.5f64, -1.5..1.5f64, -0.5..0.5f64),
                direction in (-1.0..1.0f64, -1.0..1.0f64, -1.0..1.0f64),
                big in 0.55..0.95f64,
            ) {
                let d = [direction.0, direction.1, direction.2];
                let length = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                prop_assume!(length > 0.1);
                let d = d.map(|x| x / length);
                let m = [offset.0, offset.1, offset.2];
                let along = m[0] * d[0] + m[1] * d[1] + m[2] * d[2];
                let m = [0, 1, 2].map(|k| m[k] - along * d[k]);
                let c = torus_quartic(m, d, big);
                let new = on_torus(&c);
                let old = companion(&c);
                for r in &new {
                    let (p, _, _) = value_and_slope(&c, *r);
                    let size = c.iter().map(|x| x.abs()).sum::<f64>() * r.abs().max(1.0).powi(4);
                    prop_assert!(p.abs() <= 1e-9 * size, "{r}: {p} of {size}");
                }
                prop_assert!(new.windows(2).all(|w| w[1] > w[0]), "{new:?}");
                for r in &old {
                    prop_assert!(
                        new.iter().any(|s| (s - r).abs() <= 1e-6),
                        "{r} of {old:?} missing from {new:?}"
                    );
                }
            }
        }

        /// Lines against a torus: through the hole, along the axis, grazing
        /// the tube's top and the outer equator, and across the middle.
        #[test]
        fn line_torus_quartics() {
            let big = 0.75;
            let small = 0.25;
            // Across the middle in the torus's plane: four crossings.
            let found = on_torus(&torus_quartic([0.0; 3], [1.0, 0.0, 0.0], big));
            assert!(
                matches(&found, &[-1.0, -0.5, 0.5, 1.0], 1e-12, 1.0),
                "{found:?}"
            );
            // Along the axis, and steeply through the hole: no crossing.
            assert!(on_torus(&torus_quartic([0.0; 3], [0.0, 0.0, 1.0], big)).is_empty());
            let found = on_torus(&torus_quartic([0.1, 0.05, 0.0], [0.1, 0.0, 0.995], big));
            assert!(found.is_empty(), "{found:?}");
            // Along the tube's top: tangent twice.
            let found = on_torus(&torus_quartic([0.0, 0.0, small], [1.0, 0.0, 0.0], big));
            assert!(matches(&found, &[-big, big], 1e-7, 1.0), "{found:?}");
            // Grazing the outer equator: one tangency.
            let found = on_torus(&torus_quartic([0.0, 1.0, 0.0], [1.0, 0.0, 0.0], big));
            assert!(matches(&found, &[0.0], 1e-7, 1.0), "{found:?}");
            // Grazing the inner equator: one tangency between two crossings.
            let found = on_torus(&torus_quartic([0.0, 0.5, 0.0], [1.0, 0.0, 0.0], big));
            let outer = 0.75_f64.sqrt();
            assert!(
                matches(&found, &[-outer, 0.0, outer], 1e-7, 1.0),
                "{found:?}"
            );
        }
    }
}
