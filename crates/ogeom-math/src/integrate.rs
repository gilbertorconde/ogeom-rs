//! Numerical integration.
//!
//! Gauss–Legendre quadrature, applied adaptively. A kernel integrates for arc
//! length, for area and volume over a parametric patch, and for the moments
//! that follow from those — all integrands that are smooth almost everywhere
//! and awkward exactly where a feature is.
//!
//! # Why Gauss rather than Simpson
//!
//! An `n`-point Gauss rule is exact for polynomials of degree `2n - 1`, against
//! Simpson's `3` for the same three evaluations. Since a B-spline *is* a
//! piecewise polynomial and the speed of a curve along one is a square root of
//! a polynomial, the integrands here are close enough to polynomial that the
//! difference is large. Ten points integrate most spans to machine precision in
//! one go.
//!
//! # Why adaptive on top
//!
//! A fixed rule cannot report its own error. Halving the interval and comparing
//! can: if the whole and the two halves agree, the estimate has converged
//! *there*, and if they do not, only that part is subdivided. So a curve that
//! is straight over most of its length and sharp in one place costs what the
//! sharp place costs, not what the sharp place would cost applied everywhere.
//!
//! The recursion is bounded, and a result that hit the bound says so rather
//! than being returned as though it converged.

use ogeom_core::{OgeomResult, ogeom_bail};

/// Nodes of the ten-point Gauss–Legendre rule on `[-1, 1]`, positive half.
///
/// The rule is symmetric, so the negative nodes are these negated and the
/// weights are shared. Values are the standard ones — roots of the degree-ten
/// Legendre polynomial — quoted to full `f64` precision.
const NODES: [f64; 5] = [
    0.148_874_338_981_631_21,
    0.433_395_394_129_247_2,
    0.679_409_568_299_024_4,
    0.865_063_366_688_984_5,
    0.973_906_528_517_171_7,
];

/// Weights matching [`NODES`].
const WEIGHTS: [f64; 5] = [
    0.295_524_224_714_752_87,
    0.269_266_719_309_996_35,
    0.219_086_362_515_982_04,
    0.149_451_349_150_580_6,
    0.066_671_344_308_688_14,
];

/// The most times [`integrate`] will subdivide one interval.
///
/// An integrand that has not converged by here has a singularity rather than a
/// resolution problem, and the depth limit turns that into a reported failure
/// instead of an exhausted stack.
const MAX_DEPTH: u32 = 24;

/// Integrate `f` over `[a, b]` with the fixed ten-point rule.
///
/// Exact for polynomials up to degree nineteen. No error estimate — for that,
/// use [`integrate`], which is this applied adaptively.
///
/// A reversed interval integrates to the negative, as it should: the rule
/// carries the sign of `b - a` rather than quietly sorting its arguments.
pub fn gauss_legendre<F: FnMut(f64) -> f64>(mut f: F, a: f64, b: f64) -> f64 {
    let half = (b - a) * 0.5;
    let middle = f64::midpoint(a, b);
    let mut total = 0.0;
    for (node, weight) in NODES.iter().zip(&WEIGHTS) {
        let offset = half * node;
        total += weight * (f(middle - offset) + f(middle + offset));
    }
    total * half
}

/// Integrate `f` over `[a, b]` to an absolute tolerance.
///
/// Subdivides where — and only where — the estimate has not settled, so a
/// mostly-smooth integrand costs about what the smooth part costs.
///
/// # What it will not do
///
/// Each half is given half its parent's budget, so the budgets sum to the one
/// asked for and the result is bounded by it. The cost is that an integrand
/// with an *infinite derivative* at an endpoint — `sqrt(1 - x^2)` at `x = 1`,
/// which is a circle's own equation — has a budget shrinking faster than its
/// error does, and cannot be squeezed arbitrarily. In practice it manages
/// about `1e-7` on that shape, and lands within `1e-14` when it does; asked for
/// `1e-8` it reports that it could not rather than returning the number it
/// reached.
///
/// This does not affect arc length, which is what the routine is mostly for:
/// the speed along a curve is `|c'(u)|`, smooth and positive wherever the
/// parameterization is regular. A singularity here means a genuinely singular
/// parameterization, which is worth being told about.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if the interval is not finite;
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if some part of it did not
/// converge within the depth limit. That is reported rather than returned as a
/// number, because an integral that silently stopped improving is the shape of
/// answer that gets trusted.
pub fn integrate<F: FnMut(f64) -> f64>(
    mut f: F,
    a: f64,
    b: f64,
    tolerance: f64,
) -> OgeomResult<f64> {
    if !a.is_finite() || !b.is_finite() {
        ogeom_bail!(Domain, "cannot integrate over [{a}, {b}]");
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        ogeom_bail!(Domain, "integration tolerance {tolerance} must be positive");
    }
    if a == b {
        return Ok(0.0);
    }
    let whole = gauss_legendre(&mut f, a, b);
    refine(&mut f, a, b, tolerance, whole, 0)
}

/// One step of the adaptive halving.
fn refine<F: FnMut(f64) -> f64>(
    f: &mut F,
    a: f64,
    b: f64,
    tolerance: f64,
    whole: f64,
    depth: u32,
) -> OgeomResult<f64> {
    let middle = f64::midpoint(a, b);
    let left = gauss_legendre(&mut *f, a, middle);
    let right = gauss_legendre(&mut *f, middle, b);
    let split = left + right;

    if (split - whole).abs() <= tolerance {
        // Richardson: the halved estimate is the better one, and the difference
        // is a fair estimate of what is still missing from it.
        return Ok(split + (split - whole) / 1023.0);
    }
    // Nothing left here worth resolving. A Gauss rule does not converge in
    // *relative* terms against a square-root singularity — the error stays a
    // roughly fixed fraction of the contribution — so an interval containing
    // one can fail the comparison above at every depth, while the quantity it
    // is failing about shrinks to nothing. Once the total magnitude on this
    // interval is inside the budget, no amount of refining it can move the
    // answer by more than the budget, so refining it is not worth doing.
    //
    // Compared against `|left| + |right|` rather than `|split|`: two halves
    // that cancel would otherwise look negligible while each is large.
    if left.abs() + right.abs() <= tolerance {
        return Ok(split);
    }
    if depth >= MAX_DEPTH {
        ogeom_bail!(
            NotDone,
            "the integral over [{a}, {b}] did not converge to {tolerance} \
             within {MAX_DEPTH} subdivisions; the integrand has a singularity \
             there rather than a resolution problem"
        );
    }
    // Half the tolerance to each half, so the halves' errors sum to the whole's
    // rather than each being allowed the whole budget.
    let half = tolerance * 0.5;
    Ok(
        refine(f, a, middle, half, left, depth + 1)?
            + refine(f, middle, b, half, right, depth + 1)?,
    )
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use core::f64::consts::PI;

    #[test]
    fn a_polynomial_within_the_rules_degree_is_exact_in_one_go() {
        // Degree nineteen is what a ten-point rule integrates exactly, and
        // "exactly" here should mean to rounding, not to a tolerance.
        let f = |x: f64| x.powi(19) + 3.0 * x.powi(4) - 7.0 * x + 2.0;
        let exact = 1.0 / 20.0 + 3.0 / 5.0 - 7.0 / 2.0 + 2.0;
        assert_relative_eq!(gauss_legendre(f, 0.0, 1.0), exact, epsilon = 1e-14);
    }

    #[test]
    fn transcendental_integrands_converge() {
        assert_relative_eq!(
            integrate(f64::sin, 0.0, PI, 1e-12).unwrap(),
            2.0,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            integrate(|x| 1.0 / x, 1.0, core::f64::consts::E, 1e-12).unwrap(),
            1.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn an_infinite_derivative_at_an_endpoint_is_handled_to_a_stated_limit() {
        // The quarter circle. Its integrand's derivative blows up at x = 1, so
        // the halved budget shrinks faster than the error there does and the
        // method has a floor. Where it converges it is far better than asked —
        // and where it does not, it says so instead of returning what it
        // reached, which is the whole difference between a limit and a bug.
        let quarter = |x: f64| (1.0 - x * x).max(0.0).sqrt();
        let found = integrate(quarter, 0.0, 1.0, 1e-7).unwrap();
        assert_relative_eq!(found, PI / 4.0, epsilon = 1e-12);
        assert!(
            integrate(quarter, 0.0, 1.0, 1e-8).is_err(),
            "asked for more than the method can give, it should say so"
        );
    }

    #[test]
    fn a_reversed_interval_integrates_to_the_negative() {
        // Rather than being quietly sorted, which would make an arc length
        // computed backwards come out positive and hide the caller's mistake.
        let forward = integrate(f64::sin, 0.0, PI, 1e-12).unwrap();
        let backward = integrate(f64::sin, PI, 0.0, 1e-12).unwrap();
        assert_relative_eq!(forward, -backward, epsilon = 1e-12);
    }

    #[test]
    fn an_empty_interval_integrates_to_nothing() {
        assert_eq!(integrate(f64::sin, 1.0, 1.0, 1e-12).unwrap(), 0.0);
    }

    #[test]
    fn an_integrand_that_will_not_converge_says_so() {
        // 1/x towards zero has no finite integral. Returning a large number
        // would be worse than failing, because a caller cannot tell it apart
        // from a genuinely large answer.
        assert!(integrate(|x| 1.0 / x, 0.0, 1.0, 1e-12).is_err());
    }

    #[test]
    fn non_finite_bounds_and_tolerances_are_refused() {
        assert!(integrate(f64::sin, 0.0, f64::NAN, 1e-9).is_err());
        assert!(integrate(f64::sin, f64::NEG_INFINITY, 0.0, 1e-9).is_err());
        assert!(integrate(f64::sin, 0.0, 1.0, 0.0).is_err());
        assert!(integrate(f64::sin, 0.0, 1.0, -1.0).is_err());
    }
}
