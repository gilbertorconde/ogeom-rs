//! Numerical integration.
//!
//! Gauss–Legendre quadrature, applied adaptively. A kernel integrates for arc
//! length, for area and volume over a parametric patch, and for the moments
//! that follow from those: all integrands that are smooth almost everywhere
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
//! A fixed rule cannot report its own error. A Gauss–Kronrod pair can: the
//! fifteen-point Kronrod rule shares the seven-point Gauss rule's nodes, so
//! one pass over an interval yields two estimates, and their difference is
//! a fair measure of what the better one still misses. Where that is inside
//! the budget the estimate has converged *there*; where it is not, only that
//! part is subdivided. So a curve that is straight over most of its length
//! and sharp in one place costs what the sharp place costs, not what the
//! sharp place would cost applied everywhere. The same pair, in tensor
//! form, integrates over a rectangle of parameters the same way.
//!
//! The recursion is bounded, and a result that hit the bound says so rather
//! than being returned as though it converged.

use ogeom_core::{OgeomResult, ogeom_bail};

/// Nodes of the ten-point Gauss–Legendre rule on `[-1, 1]`, positive half.
///
/// The rule is symmetric, so the negative nodes are these negated and the
/// weights are shared. Values are the standard ones (roots of the degree-ten
/// Legendre polynomial), quoted to full `f64` precision.
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

/// Nodes of the fifteen-point Kronrod rule on `[-1, 1]`, positive half,
/// outermost first; every other one from the second is a node of the
/// seven-point Gauss rule it extends.
const KRONROD_NODES: [f64; 8] = [
    0.991_455_371_120_812_6,
    0.949_107_912_342_758_5,
    0.864_864_423_359_769_1,
    0.741_531_185_599_394_5,
    0.586_087_235_467_691_1,
    0.405_845_151_377_397_2,
    0.207_784_955_007_898_48,
    0.0,
];

/// Weights matching [`KRONROD_NODES`].
const KRONROD_WEIGHTS: [f64; 8] = [
    0.022_935_322_010_529_224,
    0.063_092_092_629_978_56,
    0.104_790_010_322_250_19,
    0.140_653_259_715_525_92,
    0.169_004_726_639_267_9,
    0.190_350_578_064_785_42,
    0.204_432_940_075_298_89,
    0.209_482_141_084_727_82,
];

/// Weights of the seven-point Gauss rule, at the Kronrod nodes of odd
/// index (the second, fourth, sixth and the centre).
const GAUSS7_WEIGHTS: [f64; 4] = [
    0.129_484_966_168_869_7,
    0.279_705_391_489_276_64,
    0.381_830_050_505_118_9,
    0.417_959_183_673_469_4,
];

/// The most times [`integrate`] will subdivide one interval.
///
/// An integrand that has not converged by here has a singularity rather than a
/// resolution problem, and the depth limit turns that into a reported failure
/// instead of an exhausted stack.
const MAX_DEPTH: u32 = 24;

/// The most times [`integrate_2d`] will halve one cell, either way: a
/// singular integrand is told about at the bottom rather than pursued.
const MAX_DEPTH_2D: u32 = 24;

/// Integrate `f` over `[a, b]` with the fixed ten-point rule.
///
/// Exact for polynomials up to degree nineteen. No error estimate; for that,
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

/// Integrate `f` over `[a, b]` with the seven-point Gauss and fifteen-point
/// Kronrod pair: the Kronrod value, and the magnitude of its difference
/// from the Gauss value as the estimate of what it still misses.
///
/// Fifteen evaluations, shared. The Kronrod rule is exact for polynomials
/// up to degree twenty-two, the Gauss rule up to thirteen; where the two
/// agree the integrand is polynomial enough that both are right, and the
/// gap is a fair measure of the error where they are not. A reversed
/// interval integrates to the negative, as [`gauss_legendre`] does.
pub fn gauss_kronrod<F: FnMut(f64) -> f64>(mut f: F, a: f64, b: f64) -> (f64, f64) {
    let half = (b - a) * 0.5;
    let middle = f64::midpoint(a, b);
    let mut kronrod = 0.0;
    let mut gauss = 0.0;
    for (i, (node, weight)) in KRONROD_NODES.iter().zip(&KRONROD_WEIGHTS).enumerate() {
        let offset = half * node;
        let pair = if *node == 0.0 {
            f(middle)
        } else {
            f(middle - offset) + f(middle + offset)
        };
        kronrod += weight * pair;
        if i % 2 == 1 {
            gauss += GAUSS7_WEIGHTS[i / 2] * pair;
        }
    }
    (kronrod * half, ((kronrod - gauss) * half).abs())
}

/// Integrate `f` over `[a, b]` to an absolute tolerance.
///
/// Subdivides where, and only where, the estimate has not settled, so a
/// mostly-smooth integrand costs about what the smooth part costs.
///
/// # What it will not do
///
/// Each half is given half its parent's budget, so the budgets sum to the one
/// asked for and the result is bounded by it. The cost is that an integrand
/// with an *infinite derivative* at an endpoint (`sqrt(1 - x^2)` at `x = 1`,
/// which is a circle's own equation) has a budget shrinking faster than its
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
    refine(&mut f, a, b, tolerance, 0)
}

/// One step of the adaptive halving.
fn refine<F: FnMut(f64) -> f64>(
    f: &mut F,
    a: f64,
    b: f64,
    tolerance: f64,
    depth: u32,
) -> OgeomResult<f64> {
    let (value, error) = gauss_kronrod(&mut *f, a, b);
    if error <= tolerance {
        return Ok(value);
    }
    // Nothing left here worth resolving. A Gauss rule does not converge in
    // *relative* terms against a square-root singularity (the error stays a
    // roughly fixed fraction of the contribution), so an interval containing
    // one can fail the comparison above at every depth, while the quantity it
    // is failing about shrinks to nothing. Once the total magnitude on this
    // interval is inside the budget, no amount of refining it can move the
    // answer by more than the budget, so refining it is not worth doing.
    if value.abs() + error <= tolerance {
        return Ok(value);
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
    let middle = f64::midpoint(a, b);
    let half = tolerance * 0.5;
    Ok(refine(f, a, middle, half, depth + 1)? + refine(f, middle, b, half, depth + 1)?)
}

/// Integrate `f(u, v)` over the rectangle `[a, b] x [c, d]` to an absolute
/// tolerance: a patch's area, its moments, anything spread over a chart.
///
/// The Gauss–Kronrod pair in tensor form: one pass over a cell is fifteen
/// by fifteen evaluations and yields the Kronrod estimate and, from the
/// same values, the estimate with the Gauss rule in `u` and the one with
/// it in `v`, each gap the error owed to that direction. A cell whose
/// worse gap is inside its budget is done; one whose is not is halved
/// *along the rougher direction*, each half given half the budget, so the
/// cells' errors sum to the whole's and a ridge running across the chart
/// (a crease in an integrand, a seam) costs a line of cells rather than a
/// field of them.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if a bound is not
/// finite or the tolerance is not positive;
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if some cell did
/// not converge within the depth limit.
pub fn integrate_2d<F: FnMut(f64, f64) -> f64>(
    mut f: F,
    u: (f64, f64),
    v: (f64, f64),
    tolerance: f64,
) -> OgeomResult<f64> {
    let (a, b) = u;
    let (c, d) = v;
    if !a.is_finite() || !b.is_finite() || !c.is_finite() || !d.is_finite() {
        ogeom_bail!(Domain, "cannot integrate over [{a}, {b}] x [{c}, {d}]");
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        ogeom_bail!(Domain, "integration tolerance {tolerance} must be positive");
    }
    if a == b || c == d {
        return Ok(0.0);
    }
    refine_2d(&mut f, (a, b), (c, d), tolerance, 0)
}

/// The tensor Gauss–Kronrod pair over one cell: the Kronrod value, and the
/// gaps to the estimates with the Gauss rule in `u` and in `v`.
fn gauss_kronrod_2d<F: FnMut(f64, f64) -> f64>(
    f: &mut F,
    (a, b): (f64, f64),
    (c, d): (f64, f64),
) -> (f64, [f64; 2]) {
    let (half_u, middle_u) = ((b - a) * 0.5, f64::midpoint(a, b));
    let (half_v, middle_v) = ((d - c) * 0.5, f64::midpoint(c, d));
    // Each node's two mirror images, or the centre once.
    let stations = |half: f64, middle: f64, node: f64| -> [Option<f64>; 2] {
        if node == 0.0 {
            [Some(middle), None]
        } else {
            [Some(middle - half * node), Some(middle + half * node)]
        }
    };
    let mut kronrod = 0.0;
    let mut gauss_u = 0.0;
    let mut gauss_v = 0.0;
    for (i, (node_u, weight_u)) in KRONROD_NODES.iter().zip(&KRONROD_WEIGHTS).enumerate() {
        for (j, (node_v, weight_v)) in KRONROD_NODES.iter().zip(&KRONROD_WEIGHTS).enumerate() {
            let mut cell = 0.0;
            for uu in stations(half_u, middle_u, *node_u).into_iter().flatten() {
                for vv in stations(half_v, middle_v, *node_v).into_iter().flatten() {
                    cell += f(uu, vv);
                }
            }
            kronrod += weight_u * weight_v * cell;
            if i % 2 == 1 {
                gauss_u += GAUSS7_WEIGHTS[i / 2] * weight_v * cell;
            }
            if j % 2 == 1 {
                gauss_v += weight_u * GAUSS7_WEIGHTS[j / 2] * cell;
            }
        }
    }
    let scale = half_u * half_v;
    (
        kronrod * scale,
        [
            ((kronrod - gauss_u) * scale).abs(),
            ((kronrod - gauss_v) * scale).abs(),
        ],
    )
}

/// One step of the adaptive quartering.
fn refine_2d<F: FnMut(f64, f64) -> f64>(
    f: &mut F,
    (a, b): (f64, f64),
    (c, d): (f64, f64),
    tolerance: f64,
    depth: u32,
) -> OgeomResult<f64> {
    let (value, [error_u, error_v]) = gauss_kronrod_2d(f, (a, b), (c, d));
    let error = error_u.max(error_v);
    if error <= tolerance || value.abs() + error <= tolerance {
        return Ok(value);
    }
    if depth >= MAX_DEPTH_2D {
        ogeom_bail!(
            NotDone,
            "the integral over [{a}, {b}] x [{c}, {d}] did not converge to \
             {tolerance} within {MAX_DEPTH_2D} halvings; the integrand has a \
             singularity there rather than a resolution problem"
        );
    }
    let half = tolerance * 0.5;
    if error_u >= error_v {
        let mu = f64::midpoint(a, b);
        Ok(refine_2d(f, (a, mu), (c, d), half, depth + 1)?
            + refine_2d(f, (mu, b), (c, d), half, depth + 1)?)
    } else {
        let mv = f64::midpoint(c, d);
        Ok(refine_2d(f, (a, b), (c, mv), half, depth + 1)?
            + refine_2d(f, (a, b), (mv, d), half, depth + 1)?)
    }
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
        // method has a floor. Where it converges it is far better than asked;
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

    /// The pair's constants, checked by what they must do: the weights sum
    /// to the interval, the Gauss rule is exact to degree thirteen and the
    /// Kronrod rule to twenty-two, and the gap between them is zero where
    /// both are exact.
    #[test]
    fn the_gauss_kronrod_pair_is_exact_to_its_degrees() {
        let (weight_sum, _) = gauss_kronrod(|_| 1.0, -1.0, 1.0);
        assert_relative_eq!(weight_sum, 2.0, epsilon = 1e-15);
        let thirteen = |x: f64| x.powi(13) + 2.0 * x.powi(8) - x;
        let (value, error) = gauss_kronrod(thirteen, 0.0, 1.0);
        assert_relative_eq!(value, 1.0 / 14.0 + 2.0 / 9.0 - 0.5, epsilon = 1e-14);
        assert!(error < 1e-14, "both rules exact, no gap: {error}");
        let twenty_two = |x: f64| x.powi(22) - 3.0 * x.powi(17);
        let (value, error) = gauss_kronrod(twenty_two, 0.0, 1.0);
        assert_relative_eq!(value, 1.0 / 23.0 - 3.0 / 18.0, epsilon = 1e-14);
        assert!(
            error > 1e-6,
            "the Gauss rule is not exact here, and the gap says so: {error}"
        );
    }

    /// A rectangle of parameters: exact for a product of polynomials in one
    /// pass, the area of a sphere from its own chart, and a ridge that
    /// forces quartering on one side of the cell only.
    #[test]
    fn a_rectangle_integrates_to_a_stated_tolerance() {
        let product =
            integrate_2d(|u, v| u * u * v * v * v, (0.0, 1.0), (0.0, 1.0), 1e-12).unwrap();
        assert_relative_eq!(product, 1.0 / 12.0, epsilon = 1e-13);
        let sphere = integrate_2d(|_, v| v.sin(), (0.0, 2.0 * PI), (0.0, PI), 1e-10).unwrap();
        assert_relative_eq!(sphere, 4.0 * PI, epsilon = 1e-9);
        let ridge = integrate_2d(|u, v| (u - 0.3).abs() + v, (0.0, 1.0), (0.0, 1.0), 1e-9).unwrap();
        // ∫|u − 0.3| du over [0, 1] = 0.045 + 0.245 = 0.29; ∫ v dv = 0.5.
        assert_relative_eq!(ridge, 0.29 + 0.5, epsilon = 1e-8);
        assert!(integrate_2d(|u, _| 1.0 / u, (0.0, 1.0), (0.0, 1.0), 1e-9).is_err());
        assert_eq!(
            integrate_2d(|u, v| u + v, (1.0, 1.0), (0.0, 1.0), 1e-9).unwrap(),
            0.0
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
