//! B-spline algorithms over control points: evaluation, refinement, elevation.
//!
//! Everything here is generic over [`Blend`], the affine structure a control
//! point needs. That is what lets one implementation serve curves and surfaces,
//! 2D and 3D, and — through the homogeneous trick — rational and non-rational
//! alike, instead of four near-copies that drift apart.
//!
//! # Rational curves
//!
//! A rational B-spline is a non-rational one in one higher dimension: weight
//! each control point, carry the weight as an extra coordinate, evaluate as
//! usual, then divide through. Every algorithm here therefore applies unchanged
//! to rational geometry via [`Weighted`], which matters because exact circles,
//! cylinders and spheres are *only* representable rationally.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};

use crate::{KnotVector, Point, Point2, Vector, Vector2};

/// The affine structure a control point needs: scaling and addition.
///
/// Implemented for vectors, points and scalars. de Boor and the refinement
/// algorithms take only affine combinations — coefficients summing to one — so
/// applying them to positions is meaningful even though positions have no
/// meaningful sum on their own.
pub trait Blend: Copy {
    /// The additive identity.
    fn zero() -> Self;
    /// Scale by a factor.
    fn scale(self, k: f64) -> Self;
    /// Add another value.
    fn add(self, other: Self) -> Self;

    /// `self * (1 - t) + other * t`.
    #[must_use]
    fn lerp(self, other: Self, t: f64) -> Self {
        self.scale(1.0 - t).add(other.scale(t))
    }

    /// Subtract, via scaling by `-1`.
    #[must_use]
    fn sub(self, other: Self) -> Self {
        self.add(other.scale(-1.0))
    }
}

impl Blend for f64 {
    fn zero() -> Self {
        0.0
    }
    fn scale(self, k: f64) -> Self {
        self * k
    }
    fn add(self, other: Self) -> Self {
        self + other
    }
}

impl Blend for Vector {
    fn zero() -> Self {
        Self::ZERO
    }
    fn scale(self, k: f64) -> Self {
        self * k
    }
    fn add(self, other: Self) -> Self {
        self + other
    }
}

impl Blend for Vector2 {
    fn zero() -> Self {
        Self::ZERO
    }
    fn scale(self, k: f64) -> Self {
        self * k
    }
    fn add(self, other: Self) -> Self {
        self + other
    }
}

impl Blend for Point {
    fn zero() -> Self {
        Self::ORIGIN
    }
    fn scale(self, k: f64) -> Self {
        Self::from_vector(self.to_vector() * k)
    }
    fn add(self, other: Self) -> Self {
        Self::from_vector(self.to_vector() + other.to_vector())
    }
}

impl Blend for Point2 {
    fn zero() -> Self {
        Self::ORIGIN
    }
    fn scale(self, k: f64) -> Self {
        Self::from_vector(self.to_vector() * k)
    }
    fn add(self, other: Self) -> Self {
        Self::from_vector(self.to_vector() + other.to_vector())
    }
}

/// A control point carrying a weight, for rational geometry.
///
/// Stored in *homogeneous* form — the point is already multiplied through by
/// the weight — because that is the form every algorithm needs, and converting
/// on each access would be both slower and a source of drift.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Weighted<P> {
    /// The point scaled by the weight.
    pub scaled: P,
    /// The weight.
    pub weight: f64,
}

impl<P: Blend> Weighted<P> {
    /// A weighted control point from a position and a weight.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the weight
    /// is not finite and positive. A zero weight makes the projection undefined
    /// and a negative one makes the curve leave its control polygon's convex
    /// hull, so neither is admitted.
    pub fn new(point: P, weight: f64, tol: Tolerances) -> OgeomResult<Self> {
        if !weight.is_finite() || weight <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "control point weight {weight} must be finite and positive"
            );
        }
        Ok(Self {
            scaled: point.scale(weight),
            weight,
        })
    }

    /// The unweighted position.
    #[must_use]
    pub fn point(self) -> P {
        self.scaled.scale(1.0 / self.weight)
    }
}

impl<P: Blend> Blend for Weighted<P> {
    fn zero() -> Self {
        Self {
            scaled: P::zero(),
            weight: 0.0,
        }
    }
    fn scale(self, k: f64) -> Self {
        Self {
            scaled: self.scaled.scale(k),
            weight: self.weight * k,
        }
    }
    fn add(self, other: Self) -> Self {
        Self {
            scaled: self.scaled.add(other.scaled),
            weight: self.weight + other.weight,
        }
    }
}

/// Check that a control point count matches a knot vector.
fn check_shape<P>(knots: &KnotVector, control: &[P]) -> OgeomResult<()> {
    if control.len() != knots.control_point_count() {
        ogeom_bail!(
            Dimension,
            "knot vector describes {} control points, got {}",
            knots.control_point_count(),
            control.len()
        );
    }
    Ok(())
}

/// Evaluate a B-spline at `u` by de Boor's algorithm.
///
/// Numerically the right way to do it: a sequence of convex combinations of
/// control points, so the result stays inside their hull and no intermediate
/// can blow up. Evaluating the basis functions and taking a weighted sum gives
/// the same answer in exact arithmetic but is less stable, and expanding the
/// polynomial in monomials is far worse.
///
/// # Errors
///
/// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) if the control point
/// count disagrees with the knot vector; [`OgeomError::Domain`](ogeom_core::OgeomError::Domain)
/// if `u` is outside the domain.
pub fn evaluate<P: Blend>(
    knots: &KnotVector,
    control: &[P],
    u: f64,
    tol: Tolerances,
) -> OgeomResult<P> {
    check_shape(knots, control)?;
    let span = knots.span(u, tol)?;
    let p = knots.degree();

    let mut d: Vec<P> = (0..=p).map(|i| control[span - p + i]).collect();
    let k = knots.knots();
    for r in 1..=p {
        for j in (r..=p).rev() {
            let left = k[span + j - p];
            let right = k[span + j + 1 - r];
            // The span lookup guarantees this width is positive: a zero would
            // mean a knot of multiplicity above the degree, which the knot
            // vector's own validation rejects.
            let alpha = (u - left) / (right - left);
            d[j] = d[j - 1].lerp(d[j], alpha);
        }
    }
    Ok(d[p])
}

/// Evaluate a B-spline and its derivatives up to order `n`.
///
/// `result[0]` is the point; `result[k]` is the `k`th derivative. Orders above
/// the degree are zero.
///
/// # Errors
///
/// As [`evaluate`].
pub fn derivatives<P: Blend>(
    knots: &KnotVector,
    control: &[P],
    u: f64,
    n: usize,
    tol: Tolerances,
) -> OgeomResult<Vec<P>> {
    check_shape(knots, control)?;
    let span = knots.span(u, tol)?;
    let p = knots.degree();
    let basis = knots.basis_derivatives(span, u, n);

    Ok((0..=n)
        .map(|order| {
            let mut sum = P::zero();
            for i in 0..=p {
                sum = sum.add(control[span - p + i].scale(basis[order][i]));
            }
            sum
        })
        .collect())
}

/// Insert `value` into the knot vector `count` times, adjusting control points
/// so the curve is unchanged.
///
/// Boehm's algorithm. The foundation of nearly everything else: splitting a
/// curve, converting to Bézier form, and raising continuity constraints all
/// reduce to knot insertion.
///
/// # Errors
///
/// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) on a shape mismatch,
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if `value` is outside the
/// domain, and [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// insertion would push a multiplicity above the degree.
pub fn insert_knot<P: Blend>(
    knots: &KnotVector,
    control: &[P],
    value: f64,
    count: usize,
    tol: Tolerances,
) -> OgeomResult<Spline<P>> {
    check_shape(knots, control)?;
    if count == 0 {
        return Ok((knots.clone(), control.to_vec()));
    }
    let span = knots.span(value, tol)?;
    let p = knots.degree();
    let existing = knots.multiplicity_of(value);
    if existing + count > p {
        ogeom_bail!(
            Construction,
            "inserting {count} copies of {value} would reach multiplicity {}, above degree {p}",
            existing + count
        );
    }

    let new_knots = knots.with_knot_inserted(value, count)?;
    let last = control.len() - 1;
    let k = knots.knots();
    let (s, r) = (existing, count);

    let mut points: Vec<P> = vec![P::zero(); control.len() + r];
    // Control points outside the affected window are unchanged; those before it
    // keep their index, those after it shift right by the number inserted.
    points[..=span - p].copy_from_slice(&control[..=span - p]);
    points[span - s + r..=last + r].copy_from_slice(&control[span - s..=last]);

    // The window that the insertion actually reworks, refined in place. Each
    // pass is a set of convex combinations, so the points stay in the hull.
    let mut window: Vec<P> = (0..=p - s).map(|i| control[span - p + i]).collect();
    let mut window_start = span - p;
    for j in 1..=r {
        window_start = span - p + j;
        for i in 0..=p - j - s {
            let left = k[window_start + i];
            let right = k[i + span + 1];
            let alpha = (value - left) / (right - left);
            window[i] = window[i].lerp(window[i + 1], alpha);
        }
        points[window_start] = window[0];
        points[span + r - j - s] = window[p - j - s];
    }

    // Whatever the passes left in the middle of the window.
    if window_start + 1 < span - s {
        let width = (span - s) - (window_start + 1);
        points[window_start + 1..span - s].copy_from_slice(&window[1..=width]);
    }

    Ok((new_knots, points))
}

/// A B-spline: a knot vector paired with its control points.
pub type Spline<P> = (KnotVector, Vec<P>);

/// A Bézier segment: the parameter interval it covers, and its control points.
pub type BezierSegment<P> = ((f64, f64), Vec<P>);

/// Split a B-spline at `u` into two, each with its own clamped knot vector.
///
/// Works by raising the multiplicity at `u` to the degree, at which point the
/// control points either side are already independent.
///
/// # Errors
///
/// As [`insert_knot`], plus [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if
/// `u` is at either end of the domain, where one half would be empty.
pub fn split<P: Blend>(
    knots: &KnotVector,
    control: &[P],
    u: f64,
    tol: Tolerances,
) -> OgeomResult<(Spline<P>, Spline<P>)> {
    check_shape(knots, control)?;
    let (start, end) = knots.domain();
    if u <= start + tol.parametric() || u >= end - tol.parametric() {
        ogeom_bail!(
            Domain,
            "cannot split at {u}, an end of the domain [{start}, {end}]"
        );
    }
    let p = knots.degree();
    let existing = knots.multiplicity_of(u);
    let (refined, points) = insert_knot(knots, control, u, p - existing, tol)?;

    // After refinement the two halves meet at a control point they share.
    let cut = refined.knots().partition_point(|k| *k < u);
    let left_points = points[..cut].to_vec();
    let right_points = points[cut - 1..].to_vec();

    let mut left_knots = refined.knots()[..cut + p].to_vec();
    left_knots.push(u);
    let mut right_knots = vec![u];
    right_knots.extend_from_slice(&refined.knots()[cut..]);

    Ok((
        (KnotVector::new(left_knots, p)?, left_points),
        (KnotVector::new(right_knots, p)?, right_points),
    ))
}

/// Decompose a B-spline into its Bézier segments.
///
/// Returns one control-point array per segment, each of `degree + 1` points,
/// together with the parameter interval it covers. Many algorithms — plotting,
/// intersection, conversion to exchange formats — are far simpler on Bézier
/// pieces than on the whole spline.
///
/// # Errors
///
/// As [`insert_knot`].
pub fn to_bezier_segments<P: Blend>(
    knots: &KnotVector,
    control: &[P],
    tol: Tolerances,
) -> OgeomResult<Vec<BezierSegment<P>>> {
    check_shape(knots, control)?;
    let p = knots.degree();
    let (start, end) = knots.domain();

    // Raise every interior knot to full multiplicity; the control points then
    // partition directly into segments.
    let mut current_knots = knots.clone();
    let mut current_points = control.to_vec();
    for (value, multiplicity) in knots.distinct() {
        if value <= start || value >= end {
            continue;
        }
        let needed = p - multiplicity;
        if needed > 0 {
            let (k, c) = insert_knot(&current_knots, &current_points, value, needed, tol)?;
            current_knots = k;
            current_points = c;
        }
    }

    let breaks: Vec<f64> = core::iter::once(start)
        .chain(
            current_knots
                .distinct()
                .into_iter()
                .filter(|(v, _)| *v > start && *v < end)
                .map(|(v, _)| v),
        )
        .chain(core::iter::once(end))
        .collect();

    Ok(breaks
        .windows(2)
        .enumerate()
        .map(|(i, w)| ((w[0], w[1]), current_points[i * p..i * p + p + 1].to_vec()))
        .collect())
}

/// Raise the degree by one, leaving the curve unchanged.
///
/// Works segment by segment on the Bézier decomposition, where degree elevation
/// is the exact closed form `Q[i] = (i/(p+1)) P[i-1] + (1 - i/(p+1)) P[i]`, and
/// reassembles by removing the knots that were introduced.
///
/// # Errors
///
/// As [`to_bezier_segments`].
pub fn elevate_degree<P: Blend>(
    knots: &KnotVector,
    control: &[P],
    tol: Tolerances,
) -> OgeomResult<Spline<P>> {
    check_shape(knots, control)?;
    let p = knots.degree();
    let segments = to_bezier_segments(knots, control, tol)?;

    let mut points: Vec<P> = Vec::with_capacity(segments.len() * (p + 1) + 1);
    let mut new_knots: Vec<f64> = Vec::new();

    for (index, ((a, b), segment)) in segments.iter().enumerate() {
        // The elevated Bezier segment has p + 2 control points.
        let mut elevated: Vec<P> = Vec::with_capacity(p + 2);
        elevated.push(segment[0]);
        #[allow(clippy::cast_precision_loss)]
        for i in 1..=p {
            let t = i as f64 / (p + 1) as f64;
            elevated.push(segment[i - 1].lerp(segment[i], 1.0 - t));
        }
        elevated.push(segment[p]);

        if index == 0 {
            points.extend_from_slice(&elevated);
            new_knots.extend(core::iter::repeat_n(*a, p + 2));
        } else {
            // The shared endpoint is already present.
            points.extend_from_slice(&elevated[1..]);
            new_knots.extend(core::iter::repeat_n(*a, p + 1));
        }
        if index == segments.len() - 1 {
            new_knots.extend(core::iter::repeat_n(*b, p + 2));
        }
    }

    Ok((KnotVector::new(new_knots, p + 1)?, points))
}

/// Reverse the parameter direction, leaving the curve's shape unchanged.
#[must_use]
pub fn reverse<P: Blend>(knots: &KnotVector, control: &[P]) -> Spline<P> {
    let mut points = control.to_vec();
    points.reverse();
    (knots.reversed(), points)
}

/// Evaluate a rational B-spline: de Boor in homogeneous coordinates, then
/// divide through by the weight.
///
/// # Errors
///
/// As [`evaluate`], plus [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if the
/// accumulated weight vanishes, which positive input weights make impossible.
pub fn evaluate_rational<P: Blend>(
    knots: &KnotVector,
    control: &[Weighted<P>],
    u: f64,
    tol: Tolerances,
) -> OgeomResult<P> {
    let h = evaluate(knots, control, u, tol)?;
    if h.weight.abs() <= tol.confusion() {
        ogeom_bail!(Numeric, "rational evaluation produced a vanishing weight");
    }
    Ok(h.point())
}

/// Evaluate a rational B-spline and its derivatives up to order `n`.
///
/// The quotient rule applied to the homogeneous form. Differentiating the
/// projected curve directly is not an option: the projection is a quotient, so
/// its derivatives mix all lower orders.
///
/// # Errors
///
/// As [`evaluate_rational`].
pub fn rational_derivatives<P: Blend>(
    knots: &KnotVector,
    control: &[Weighted<P>],
    u: f64,
    n: usize,
    tol: Tolerances,
) -> OgeomResult<Vec<P>> {
    let homogeneous = derivatives(knots, control, u, n, tol)?;
    if homogeneous[0].weight.abs() <= tol.confusion() {
        ogeom_bail!(Numeric, "rational evaluation produced a vanishing weight");
    }

    // C^(k) = ( A^(k) - sum_{i=1..k} C(k,i) w^(i) C^(k-i) ) / w
    let mut out: Vec<P> = Vec::with_capacity(n + 1);
    for (order, term) in homogeneous.iter().enumerate() {
        let mut value = term.scaled;
        for i in 1..=order {
            #[allow(clippy::cast_precision_loss)]
            let binomial = binomial_coefficient(order, i) as f64;
            value = value.sub(out[order - i].scale(binomial * homogeneous[i].weight));
        }
        out.push(value.scale(1.0 / homogeneous[0].weight));
    }
    Ok(out)
}

/// `n choose k`, computed multiplicatively so it stays exact for the small
/// values derivative formulas need.
#[must_use]
pub fn binomial_coefficient(n: usize, k: usize) -> u64 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut result = 1_u64;
    for i in 0..k {
        result = result * (n - i) as u64 / (i as u64 + 1);
    }
    result
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

    fn cubic_curve() -> (KnotVector, Vec<Point>) {
        let control = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 2.0, 0.0),
            Point::new(3.0, 3.0, 1.0),
            Point::new(5.0, 1.0, 2.0),
            Point::new(6.0, -1.0, 1.0),
            Point::new(8.0, 0.0, 0.0),
        ];
        (
            KnotVector::clamped_uniform(3, control.len()).unwrap(),
            control,
        )
    }

    fn sample(knots: &KnotVector, control: &[Point], n: usize) -> Vec<Point> {
        let (a, b) = knots.domain();
        (0..=n)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let u = a + (b - a) * (i as f64 / n as f64);
                evaluate(knots, control, u, T).unwrap()
            })
            .collect()
    }

    #[test]
    fn a_clamped_curve_interpolates_its_end_points() {
        let (k, c) = cubic_curve();
        let (a, b) = k.domain();
        assert!(evaluate(&k, &c, a, T).unwrap().is_equal(c[0], T));
        assert!(evaluate(&k, &c, b, T).unwrap().is_equal(c[c.len() - 1], T));
    }

    #[test]
    fn de_boor_agrees_with_the_basis_function_sum() {
        // Two independent routes to the same value; they must agree.
        let (k, c) = cubic_curve();
        for i in 0..=50 {
            let u = f64::from(i) / 50.0;
            let span = k.span(u, T).unwrap();
            let basis = k.basis(span, u);
            let mut sum = Vector::ZERO;
            for j in 0..=k.degree() {
                sum += c[span - k.degree() + j].to_vector() * basis[j];
            }
            let de_boor = evaluate(&k, &c, u, T).unwrap();
            assert!(de_boor.is_equal(Point::from_vector(sum), T), "at u = {u}");
        }
    }

    #[test]
    fn shape_mismatches_and_out_of_domain_parameters_are_refused() {
        let (k, c) = cubic_curve();
        assert!(
            evaluate(&k, &c[..3], 0.5, T).is_err(),
            "too few control points"
        );
        assert!(evaluate(&k, &c, -0.1, T).is_err());
        assert!(evaluate(&k, &c, 1.1, T).is_err());
    }

    #[test]
    fn derivatives_agree_with_finite_differences() {
        let (k, c) = cubic_curve();
        let h = 1e-6;
        for i in 1..20 {
            let u = f64::from(i) / 20.0;
            let d = derivatives(&k, &c, u, 2, T).unwrap();
            assert!(d[0].is_equal(evaluate(&k, &c, u, T).unwrap(), T));

            let ahead = evaluate(&k, &c, u + h, T).unwrap();
            let behind = evaluate(&k, &c, u - h, T).unwrap();
            let numeric = (ahead - behind) * (1.0 / (2.0 * h));
            assert!(
                (d[1].to_vector() - numeric).magnitude() < 1e-5,
                "first derivative disagrees at {u}"
            );
        }
    }

    #[test]
    fn knot_insertion_does_not_move_the_curve() {
        let (k, c) = cubic_curve();
        let before = sample(&k, &c, 100);
        for (value, count) in [(0.25, 1), (0.5, 2), (0.75, 3), (0.1, 1)] {
            let (k2, c2) = insert_knot(&k, &c, value, count, T).unwrap();
            assert_eq!(c2.len(), c.len() + count);
            assert_eq!(k2.multiplicity_of(value), k.multiplicity_of(value) + count);
            let after = sample(&k2, &c2, 100);
            for (a, b) in before.iter().zip(&after) {
                assert!(
                    a.is_equal(*b, T),
                    "inserting {count} at {value} moved the curve"
                );
            }
        }
    }

    #[test]
    fn repeated_insertion_matches_a_single_multiple_insertion() {
        let (k, c) = cubic_curve();
        let (ka, ca) = insert_knot(&k, &c, 0.4, 3, T).unwrap();

        let (k1, c1) = insert_knot(&k, &c, 0.4, 1, T).unwrap();
        let (k2, c2) = insert_knot(&k1, &c1, 0.4, 1, T).unwrap();
        let (kb, cb) = insert_knot(&k2, &c2, 0.4, 1, T).unwrap();

        assert_eq!(ka.knots(), kb.knots());
        for (a, b) in ca.iter().zip(&cb) {
            assert!(a.is_equal(*b, T));
        }
    }

    #[test]
    fn insertion_beyond_the_degree_is_refused() {
        let (k, c) = cubic_curve();
        assert!(insert_knot(&k, &c, 0.5, 4, T).is_err());
        assert!(insert_knot(&k, &c, 0.5, 3, T).is_ok());
        assert!(
            insert_knot(&k, &c, 2.0, 1, T).is_err(),
            "outside the domain"
        );
    }

    #[test]
    fn splitting_reproduces_both_halves_of_the_original() {
        let (k, c) = cubic_curve();
        let cut = 0.4;
        let ((lk, lc), (rk, rc)) = split(&k, &c, cut, T).unwrap();

        assert_relative_eq!(lk.domain().1, cut, epsilon = 1e-15);
        assert_relative_eq!(rk.domain().0, cut, epsilon = 1e-15);
        assert!(lk.is_clamped() && rk.is_clamped());

        for i in 0..=40 {
            let t = f64::from(i) / 40.0;
            let left_u = lk.domain().0 + (cut - lk.domain().0) * t;
            let right_u = cut + (rk.domain().1 - cut) * t;
            assert!(
                evaluate(&lk, &lc, left_u, T)
                    .unwrap()
                    .is_equal(evaluate(&k, &c, left_u, T).unwrap(), T),
                "left half diverges at {left_u}"
            );
            assert!(
                evaluate(&rk, &rc, right_u, T)
                    .unwrap()
                    .is_equal(evaluate(&k, &c, right_u, T).unwrap(), T),
                "right half diverges at {right_u}"
            );
        }
    }

    #[test]
    fn splitting_at_an_end_of_the_domain_is_refused() {
        let (k, c) = cubic_curve();
        assert!(split(&k, &c, 0.0, T).is_err());
        assert!(split(&k, &c, 1.0, T).is_err());
    }

    #[test]
    fn bezier_decomposition_covers_the_curve_exactly() {
        let (k, c) = cubic_curve();
        let segments = to_bezier_segments(&k, &c, T).unwrap();
        // Two interior knots means three segments.
        assert_eq!(segments.len(), 3);
        for (_, points) in &segments {
            assert_eq!(points.len(), k.degree() + 1);
        }

        // Each segment, evaluated as a Bezier, must match the original curve
        // over its own interval.
        for ((a, b), points) in &segments {
            let bezier = KnotVector::clamped_uniform(k.degree(), points.len())
                .unwrap()
                .reparameterized(*a, *b)
                .unwrap();
            for i in 0..=20 {
                let u = a + (b - a) * (f64::from(i) / 20.0);
                assert!(
                    evaluate(&bezier, points, u, T)
                        .unwrap()
                        .is_equal(evaluate(&k, &c, u, T).unwrap(), T),
                    "segment [{a}, {b}] diverges at {u}"
                );
            }
        }
    }

    #[test]
    fn degree_elevation_does_not_move_the_curve() {
        let (k, c) = cubic_curve();
        let before = sample(&k, &c, 100);
        let (k2, c2) = elevate_degree(&k, &c, T).unwrap();
        assert_eq!(k2.degree(), k.degree() + 1);
        assert_eq!(k2.domain(), k.domain());

        let after = sample(&k2, &c2, 100);
        for (a, b) in before.iter().zip(&after) {
            assert!(a.is_equal(*b, T), "elevation moved the curve");
        }
    }

    #[test]
    fn elevation_twice_is_still_the_same_curve() {
        let (k, c) = cubic_curve();
        let before = sample(&k, &c, 60);
        let (k1, c1) = elevate_degree(&k, &c, T).unwrap();
        let (k2, c2) = elevate_degree(&k1, &c1, T).unwrap();
        assert_eq!(k2.degree(), 5);
        for (a, b) in before.iter().zip(&sample(&k2, &c2, 60)) {
            assert!(a.is_equal(*b, T));
        }
    }

    #[test]
    fn reversal_traverses_the_same_points_backwards() {
        let (k, c) = cubic_curve();
        let (rk, rc) = reverse(&k, &c);
        let (a, b) = k.domain();
        for i in 0..=40 {
            let t = f64::from(i) / 40.0;
            let forward = evaluate(&k, &c, a + (b - a) * t, T).unwrap();
            let backward = evaluate(&rk, &rc, a + (b - a) * (1.0 - t), T).unwrap();
            assert!(forward.is_equal(backward, T), "at t = {t}");
        }
    }

    /// A quarter circle, exactly, as a rational quadratic. This is the reason
    /// rational geometry exists: no polynomial curve is a circular arc.
    fn quarter_circle() -> (KnotVector, Vec<Weighted<Point>>) {
        let w = core::f64::consts::FRAC_1_SQRT_2;
        let control = vec![
            Weighted::new(Point::new(1.0, 0.0, 0.0), 1.0, T).unwrap(),
            Weighted::new(Point::new(1.0, 1.0, 0.0), w, T).unwrap(),
            Weighted::new(Point::new(0.0, 1.0, 0.0), 1.0, T).unwrap(),
        ];
        (KnotVector::clamped_uniform(2, 3).unwrap(), control)
    }

    #[test]
    fn a_rational_quadratic_traces_an_exact_circular_arc() {
        let (k, c) = quarter_circle();
        for i in 0..=100 {
            let u = f64::from(i) / 100.0;
            let p = evaluate_rational(&k, &c, u, T).unwrap();
            // Every point is at exactly unit distance from the origin — which
            // no non-rational B-spline can achieve.
            assert_relative_eq!(p.to_vector().magnitude(), 1.0, epsilon = 1e-14);
            assert_relative_eq!(p.z, 0.0, epsilon = 1e-15);
        }
        assert!(
            evaluate_rational(&k, &c, 0.0, T)
                .unwrap()
                .is_equal(Point::new(1.0, 0.0, 0.0), T)
        );
        assert!(
            evaluate_rational(&k, &c, 1.0, T)
                .unwrap()
                .is_equal(Point::new(0.0, 1.0, 0.0), T)
        );
    }

    #[test]
    fn rational_derivatives_agree_with_finite_differences() {
        let (k, c) = quarter_circle();
        let h = 1e-6;
        for i in 1..20 {
            let u = f64::from(i) / 20.0;
            let d = rational_derivatives(&k, &c, u, 2, T).unwrap();
            assert!(d[0].is_equal(evaluate_rational(&k, &c, u, T).unwrap(), T));

            let ahead = evaluate_rational(&k, &c, u + h, T).unwrap();
            let behind = evaluate_rational(&k, &c, u - h, T).unwrap();
            let numeric = (ahead - behind) * (1.0 / (2.0 * h));
            assert!(
                (d[1].to_vector() - numeric).magnitude() < 1e-5,
                "at u = {u}: {:?} vs {numeric:?}",
                d[1]
            );
        }
    }

    #[test]
    fn the_tangent_of_a_circular_arc_is_perpendicular_to_its_radius() {
        let (k, c) = quarter_circle();
        for i in 0..=20 {
            let u = f64::from(i) / 20.0;
            let d = rational_derivatives(&k, &c, u, 1, T).unwrap();
            let radius = d[0].to_vector();
            let tangent = d[1].to_vector();
            assert!(
                radius.dot(tangent).abs() < 1e-12,
                "not perpendicular at {u}: {}",
                radius.dot(tangent)
            );
        }
    }

    #[test]
    fn knot_insertion_preserves_a_rational_curve_too() {
        let (k, c) = quarter_circle();
        let (k2, c2) = insert_knot(&k, &c, 0.5, 1, T).unwrap();
        for i in 0..=50 {
            let u = f64::from(i) / 50.0;
            let a = evaluate_rational(&k, &c, u, T).unwrap();
            let b = evaluate_rational(&k2, &c2, u, T).unwrap();
            assert!(a.is_equal(b, T), "at {u}");
            assert_relative_eq!(b.to_vector().magnitude(), 1.0, epsilon = 1e-14);
        }
    }

    #[test]
    fn degenerate_weights_are_refused() {
        assert!(Weighted::new(Point::ORIGIN, 0.0, T).is_err());
        assert!(Weighted::new(Point::ORIGIN, -1.0, T).is_err());
        assert!(Weighted::new(Point::ORIGIN, f64::NAN, T).is_err());
        assert!(Weighted::new(Point::ORIGIN, f64::INFINITY, T).is_err());
        assert!(Weighted::new(Point::ORIGIN, 2.0, T).is_ok());
    }

    #[test]
    fn weighted_round_trips_through_its_homogeneous_form() {
        let p = Point::new(3.0, -1.0, 2.0);
        let w = Weighted::new(p, 2.5, T).unwrap();
        assert!(w.point().is_equal(p, T));
        assert!(w.scaled.is_equal(Point::new(7.5, -2.5, 5.0), T));
    }

    #[test]
    fn binomial_coefficients() {
        assert_eq!(binomial_coefficient(0, 0), 1);
        assert_eq!(binomial_coefficient(5, 0), 1);
        assert_eq!(binomial_coefficient(5, 5), 1);
        assert_eq!(binomial_coefficient(5, 2), 10);
        assert_eq!(binomial_coefficient(10, 5), 252);
        assert_eq!(binomial_coefficient(3, 4), 0);
    }

    #[test]
    fn scalar_and_planar_control_points_work_too() {
        // The Blend abstraction has to serve every control point type, not just
        // 3D positions.
        let k = KnotVector::clamped_uniform(2, 4).unwrap();
        let scalars = vec![0.0_f64, 1.0, 3.0, 2.0];
        assert_relative_eq!(evaluate(&k, &scalars, 0.0, T).unwrap(), 0.0);
        assert_relative_eq!(evaluate(&k, &scalars, 1.0, T).unwrap(), 2.0);

        let planar = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 2.0),
            Point2::new(3.0, 1.0),
            Point2::new(4.0, 0.0),
        ];
        assert!(
            evaluate(&k, &planar, 0.0, T)
                .unwrap()
                .is_equal(planar[0], T)
        );
        assert!(
            evaluate(&k, &planar, 1.0, T)
                .unwrap()
                .is_equal(planar[3], T)
        );
    }
}

/// A rectangular grid of control points for a tensor-product surface.
///
/// Stored row-major: `points[i * v_count + j]` is the point at `u` index `i` and
/// `v` index `j`. Carrying the shape with the data means the surface functions
/// cannot be handed a grid with the wrong stride, which is the mistake that
/// otherwise produces a plausible but transposed surface.
#[derive(Debug, Clone, PartialEq)]
pub struct ControlGrid<P> {
    points: Vec<P>,
    u_count: usize,
    v_count: usize,
}

impl<P: Blend> ControlGrid<P> {
    /// A grid from row-major points.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) if the point count
    /// is not `u_count * v_count`, or either count is zero.
    pub fn new(points: Vec<P>, u_count: usize, v_count: usize) -> OgeomResult<Self> {
        if u_count == 0 || v_count == 0 {
            ogeom_bail!(Dimension, "control grid must be at least 1x1");
        }
        if points.len() != u_count * v_count {
            ogeom_bail!(
                Dimension,
                "a {u_count}x{v_count} grid needs {} points, got {}",
                u_count * v_count,
                points.len()
            );
        }
        Ok(Self {
            points,
            u_count,
            v_count,
        })
    }

    /// Number of control points along `u`.
    #[must_use]
    pub const fn u_count(&self) -> usize {
        self.u_count
    }

    /// Number of control points along `v`.
    #[must_use]
    pub const fn v_count(&self) -> usize {
        self.v_count
    }

    /// The point at `(i, j)`, or `None` if either index is out of range.
    #[must_use]
    pub fn get(&self, i: usize, j: usize) -> Option<P> {
        if i >= self.u_count || j >= self.v_count {
            return None;
        }
        self.points.get(i * self.v_count + j).copied()
    }

    /// All points, row-major.
    #[must_use]
    pub fn points(&self) -> &[P] {
        &self.points
    }

    /// This grid with `u` and `v` exchanged.
    #[must_use]
    pub fn transposed(&self) -> Self {
        let mut points = Vec::with_capacity(self.points.len());
        for j in 0..self.v_count {
            for i in 0..self.u_count {
                points.push(self.points[i * self.v_count + j]);
            }
        }
        Self {
            points,
            u_count: self.v_count,
            v_count: self.u_count,
        }
    }

    /// Apply `f` to every point.
    #[must_use]
    pub fn map<Q: Blend>(&self, f: impl Fn(P) -> Q) -> ControlGrid<Q> {
        ControlGrid {
            points: self.points.iter().map(|p| f(*p)).collect(),
            u_count: self.u_count,
            v_count: self.v_count,
        }
    }
}

/// Check that a grid's shape matches its two knot vectors.
fn check_grid_shape<P>(ku: &KnotVector, kv: &KnotVector, grid: &ControlGrid<P>) -> OgeomResult<()> {
    if grid.u_count != ku.control_point_count() || grid.v_count != kv.control_point_count() {
        ogeom_bail!(
            Dimension,
            "knot vectors describe a {}x{} grid, got {}x{}",
            ku.control_point_count(),
            kv.control_point_count(),
            grid.u_count,
            grid.v_count
        );
    }
    Ok(())
}

/// Evaluate a tensor-product B-spline surface at `(u, v)`.
///
/// Sums the `(p+1) x (q+1)` non-zero basis products over the control window.
/// Only that window contributes — the basis has local support — so cost depends
/// on the degrees, not on the size of the surface.
///
/// # Errors
///
/// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) on a shape mismatch, and
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if a parameter is outside its
/// knot vector's domain.
pub fn evaluate_surface<P: Blend>(
    ku: &KnotVector,
    kv: &KnotVector,
    grid: &ControlGrid<P>,
    u: f64,
    v: f64,
    tol: Tolerances,
) -> OgeomResult<P> {
    check_grid_shape(ku, kv, grid)?;
    let (p, q) = (ku.degree(), kv.degree());
    let (su, sv) = (ku.span(u, tol)?, kv.span(v, tol)?);
    let (nu, nv) = (ku.basis(su, u), kv.basis(sv, v));

    let mut total = P::zero();
    for (i, &weight_u) in nu.iter().enumerate() {
        // Accumulate along v first, then weight the row: one multiply per row
        // instead of one per point.
        let mut row = P::zero();
        for (j, &weight_v) in nv.iter().enumerate() {
            let Some(point) = grid.get(su - p + i, sv - q + j) else {
                ogeom_bail!(Dimension, "control grid index out of range");
            };
            row = row.add(point.scale(weight_v));
        }
        total = total.add(row.scale(weight_u));
    }
    Ok(total)
}

/// Evaluate a surface and its partial derivatives up to total order `order`.
///
/// `result[k][l]` is the derivative taken `k` times in `u` and `l` times in `v`,
/// so `result[0][0]` is the point itself.
///
/// # Errors
///
/// As [`evaluate_surface`].
pub fn surface_derivatives<P: Blend>(
    ku: &KnotVector,
    kv: &KnotVector,
    grid: &ControlGrid<P>,
    u: f64,
    v: f64,
    order: usize,
    tol: Tolerances,
) -> OgeomResult<Vec<Vec<P>>> {
    check_grid_shape(ku, kv, grid)?;
    let (p, q) = (ku.degree(), kv.degree());
    let (su, sv) = (ku.span(u, tol)?, kv.span(v, tol)?);
    let du = ku.basis_derivatives(su, u, order);
    let dv = kv.basis_derivatives(sv, v, order);

    let mut out = vec![vec![P::zero(); order + 1]; order + 1];
    for (k, row) in out.iter_mut().enumerate() {
        for (l, cell) in row.iter_mut().enumerate() {
            // Derivatives past the degree in either direction vanish, and the
            // basis returns them as exact zeros, so this sums to zero without
            // needing a special case.
            let mut total = P::zero();
            for (i, &weight_u) in du[k].iter().enumerate() {
                let mut inner = P::zero();
                for (j, &weight_v) in dv[l].iter().enumerate() {
                    let Some(point) = grid.get(su - p + i, sv - q + j) else {
                        ogeom_bail!(Dimension, "control grid index out of range");
                    };
                    inner = inner.add(point.scale(weight_v));
                }
                total = total.add(inner.scale(weight_u));
            }
            *cell = total;
        }
    }
    Ok(out)
}

/// Evaluate a rational tensor-product surface: homogeneous evaluation, then
/// divide through.
///
/// # Errors
///
/// As [`evaluate_surface`], plus
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if the accumulated weight
/// vanishes, which positive input weights make impossible.
pub fn evaluate_rational_surface<P: Blend>(
    ku: &KnotVector,
    kv: &KnotVector,
    grid: &ControlGrid<Weighted<P>>,
    u: f64,
    v: f64,
    tol: Tolerances,
) -> OgeomResult<P> {
    let h = evaluate_surface(ku, kv, grid, u, v, tol)?;
    if h.weight.abs() <= tol.confusion() {
        ogeom_bail!(
            Numeric,
            "rational surface evaluation produced a vanishing weight"
        );
    }
    Ok(h.point())
}

/// Evaluate a rational surface and its partial derivatives up to total order
/// `order`.
///
/// The two-parameter quotient rule. Each mixed partial subtracts the weight's
/// influence in `u`, in `v`, and in both together; dropping the last of those
/// three sums is the classic error, and it only shows up on genuinely rational
/// surfaces with mixed derivatives — which is to say, on exactly the spheres and
/// tori where the answer matters.
///
/// # Errors
///
/// As [`evaluate_rational_surface`].
pub fn rational_surface_derivatives<P: Blend>(
    ku: &KnotVector,
    kv: &KnotVector,
    grid: &ControlGrid<Weighted<P>>,
    u: f64,
    v: f64,
    order: usize,
    tol: Tolerances,
) -> OgeomResult<Vec<Vec<P>>> {
    let h = surface_derivatives(ku, kv, grid, u, v, order, tol)?;
    let w0 = h[0][0].weight;
    if w0.abs() <= tol.confusion() {
        ogeom_bail!(
            Numeric,
            "rational surface evaluation produced a vanishing weight"
        );
    }

    let mut s = vec![vec![P::zero(); order + 1]; order + 1];
    for k in 0..=order {
        for l in 0..=order {
            let mut value = h[k][l].scaled;
            #[allow(clippy::cast_precision_loss)]
            for i in 1..=k {
                let c = binomial_coefficient(k, i) as f64;
                value = value.sub(s[k - i][l].scale(c * h[i][0].weight));
            }
            #[allow(clippy::cast_precision_loss)]
            for j in 1..=l {
                let c = binomial_coefficient(l, j) as f64;
                value = value.sub(s[k][l - j].scale(c * h[0][j].weight));
            }
            #[allow(clippy::cast_precision_loss)]
            for i in 1..=k {
                let ci = binomial_coefficient(k, i) as f64;
                for j in 1..=l {
                    let cj = binomial_coefficient(l, j) as f64;
                    value = value.sub(s[k - i][l - j].scale(ci * cj * h[i][j].weight));
                }
            }
            s[k][l] = value.scale(1.0 / w0);
        }
    }
    Ok(s)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod surface_tests {
    use super::*;
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

    /// A bicubic patch with some genuine curvature.
    fn patch() -> (KnotVector, KnotVector, ControlGrid<Point>) {
        let (nu, nv) = (5, 4);
        let mut points = Vec::with_capacity(nu * nv);
        for i in 0..nu {
            for j in 0..nv {
                #[allow(clippy::cast_precision_loss)]
                let (x, y) = (i as f64, j as f64);
                points.push(Point::new(x, y, (x * 0.7).sin() * (y * 0.5).cos()));
            }
        }
        (
            KnotVector::clamped_uniform(3, nu).unwrap(),
            KnotVector::clamped_uniform(2, nv).unwrap(),
            ControlGrid::new(points, nu, nv).unwrap(),
        )
    }

    #[test]
    fn grid_shape_is_checked_on_construction() {
        assert!(ControlGrid::new(vec![Point::ORIGIN; 6], 2, 3).is_ok());
        assert!(ControlGrid::new(vec![Point::ORIGIN; 6], 3, 3).is_err());
        assert!(ControlGrid::new(Vec::<Point>::new(), 0, 3).is_err());
    }

    #[test]
    fn grid_indexing_is_row_major_and_bounds_checked() {
        let g = ControlGrid::new(
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
                Point::new(0.0, 2.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
                Point::new(1.0, 1.0, 0.0),
                Point::new(1.0, 2.0, 0.0),
            ],
            2,
            3,
        )
        .unwrap();
        assert_eq!(g.get(1, 2), Some(Point::new(1.0, 2.0, 0.0)));
        assert_eq!(g.get(0, 1), Some(Point::new(0.0, 1.0, 0.0)));
        assert_eq!(g.get(2, 0), None);
        assert_eq!(g.get(0, 3), None);
    }

    #[test]
    fn transposing_twice_is_the_identity() {
        let (_, _, g) = patch();
        let t = g.transposed();
        assert_eq!(t.u_count(), g.v_count());
        assert_eq!(t.v_count(), g.u_count());
        for i in 0..g.u_count() {
            for j in 0..g.v_count() {
                assert_eq!(t.get(j, i), g.get(i, j));
            }
        }
        assert_eq!(t.transposed(), g);
    }

    #[test]
    fn a_clamped_patch_interpolates_its_corner_control_points() {
        let (ku, kv, g) = patch();
        let ((u0, u1), (v0, v1)) = (ku.domain(), kv.domain());
        let corners = [
            (u0, v0, g.get(0, 0).unwrap()),
            (u0, v1, g.get(0, g.v_count() - 1).unwrap()),
            (u1, v0, g.get(g.u_count() - 1, 0).unwrap()),
            (u1, v1, g.get(g.u_count() - 1, g.v_count() - 1).unwrap()),
        ];
        for (u, v, expected) in corners {
            assert!(
                evaluate_surface(&ku, &kv, &g, u, v, T)
                    .unwrap()
                    .is_equal(expected, T),
                "corner ({u}, {v})"
            );
        }
    }

    #[test]
    fn surface_shape_mismatches_are_refused() {
        let (ku, kv, g) = patch();
        let wrong = ControlGrid::new(g.points().to_vec(), 4, 5).unwrap();
        assert!(evaluate_surface(&ku, &kv, &wrong, 0.5, 0.5, T).is_err());
        assert!(evaluate_surface(&ku, &kv, &g, 1.5, 0.5, T).is_err());
        assert!(evaluate_surface(&ku, &kv, &g, 0.5, -0.5, T).is_err());
    }

    #[test]
    fn surface_partials_agree_with_finite_differences() {
        let (ku, kv, g) = patch();
        let h = 1e-6;
        for iu in 1..6 {
            for iv in 1..6 {
                let (u, v) = (f64::from(iu) / 6.0, f64::from(iv) / 6.0);
                let d = surface_derivatives(&ku, &kv, &g, u, v, 2, T).unwrap();
                assert!(d[0][0].is_equal(evaluate_surface(&ku, &kv, &g, u, v, T).unwrap(), T));

                let du = (evaluate_surface(&ku, &kv, &g, u + h, v, T).unwrap()
                    - evaluate_surface(&ku, &kv, &g, u - h, v, T).unwrap())
                    * (1.0 / (2.0 * h));
                let dv = (evaluate_surface(&ku, &kv, &g, u, v + h, T).unwrap()
                    - evaluate_surface(&ku, &kv, &g, u, v - h, T).unwrap())
                    * (1.0 / (2.0 * h));
                assert!((d[1][0].to_vector() - du).magnitude() < 1e-5 * du.magnitude().max(1.0));
                assert!((d[0][1].to_vector() - dv).magnitude() < 1e-5 * dv.magnitude().max(1.0));

                // The mixed partial, which the naive quotient rule drops.
                let mixed = (evaluate_surface(&ku, &kv, &g, u + h, v + h, T).unwrap()
                    - evaluate_surface(&ku, &kv, &g, u + h, v - h, T).unwrap()
                    - (evaluate_surface(&ku, &kv, &g, u - h, v + h, T).unwrap()
                        - evaluate_surface(&ku, &kv, &g, u - h, v - h, T).unwrap()))
                    * (1.0 / (4.0 * h * h));
                assert!(
                    (d[1][1].to_vector() - mixed).magnitude() < 1e-3 * mixed.magnitude().max(1.0),
                    "mixed partial wrong at ({u}, {v})"
                );
            }
        }
    }

    /// A hemisphere, exactly, as a rational biquadratic. Only a rational
    /// surface can be one.
    fn rational_hemisphere() -> (KnotVector, KnotVector, ControlGrid<Weighted<Point>>) {
        let w = core::f64::consts::FRAC_1_SQRT_2;
        // A quarter arc in u, swept through a quarter turn in v.
        let rows: [[(Point, f64); 3]; 3] = [
            [
                (Point::new(1.0, 0.0, 0.0), 1.0),
                (Point::new(1.0, 1.0, 0.0), w),
                (Point::new(0.0, 1.0, 0.0), 1.0),
            ],
            [
                (Point::new(1.0, 0.0, 1.0), w),
                (Point::new(1.0, 1.0, 1.0), w * w),
                (Point::new(0.0, 1.0, 1.0), w),
            ],
            [
                (Point::new(0.0, 0.0, 1.0), 1.0),
                (Point::new(0.0, 0.0, 1.0), w),
                (Point::new(0.0, 0.0, 1.0), 1.0),
            ],
        ];
        let points: Vec<_> = rows
            .iter()
            .flatten()
            .map(|(p, w)| Weighted::new(*p, *w, T).unwrap())
            .collect();
        (
            KnotVector::clamped_uniform(2, 3).unwrap(),
            KnotVector::clamped_uniform(2, 3).unwrap(),
            ControlGrid::new(points, 3, 3).unwrap(),
        )
    }

    #[test]
    fn a_rational_biquadratic_traces_an_exact_sphere() {
        let (ku, kv, g) = rational_hemisphere();
        for iu in 0..=10 {
            for iv in 0..=10 {
                let (u, v) = (f64::from(iu) / 10.0, f64::from(iv) / 10.0);
                let p = evaluate_rational_surface(&ku, &kv, &g, u, v, T).unwrap();
                assert_relative_eq!(
                    p.to_vector().magnitude(),
                    1.0,
                    epsilon = 1e-13,
                    max_relative = 1e-13
                );
            }
        }
    }

    #[test]
    fn rational_surface_partials_agree_with_finite_differences() {
        let (ku, kv, g) = rational_hemisphere();
        let h = 1e-6;
        let at = |u: f64, v: f64| evaluate_rational_surface(&ku, &kv, &g, u, v, T).unwrap();
        for iu in 1..6 {
            for iv in 1..6 {
                let (u, v) = (f64::from(iu) / 6.0, f64::from(iv) / 6.0);
                let d = rational_surface_derivatives(&ku, &kv, &g, u, v, 2, T).unwrap();
                assert!(d[0][0].is_equal(at(u, v), T));

                let du = (at(u + h, v) - at(u - h, v)) * (1.0 / (2.0 * h));
                let dv = (at(u, v + h) - at(u, v - h)) * (1.0 / (2.0 * h));
                assert!(
                    (d[1][0].to_vector() - du).magnitude() < 1e-5 * du.magnitude().max(1.0),
                    "du wrong at ({u}, {v})"
                );
                assert!(
                    (d[0][1].to_vector() - dv).magnitude() < 1e-5 * dv.magnitude().max(1.0),
                    "dv wrong at ({u}, {v})"
                );

                // The mixed partial is where the cross term in the two-parameter
                // quotient rule matters; without it this is visibly wrong.
                let mixed =
                    (at(u + h, v + h) - at(u + h, v - h) - (at(u - h, v + h) - at(u - h, v - h)))
                        * (1.0 / (4.0 * h * h));
                assert!(
                    (d[1][1].to_vector() - mixed).magnitude() < 1e-2 * mixed.magnitude().max(1.0),
                    "mixed partial wrong at ({u}, {v}): {:?} vs {mixed:?}",
                    d[1][1]
                );
            }
        }
    }

    #[test]
    fn a_spheres_normal_is_radial() {
        // Independent of the derivative formulas: on a unit sphere centred at
        // the origin, du x dv must be parallel to the position vector.
        let (ku, kv, g) = rational_hemisphere();
        for iu in 1..8 {
            for iv in 1..8 {
                let (u, v) = (f64::from(iu) / 8.0, f64::from(iv) / 8.0);
                let d = rational_surface_derivatives(&ku, &kv, &g, u, v, 1, T).unwrap();
                let radius = d[0][0].to_vector();
                let normal = d[1][0].to_vector().cross(d[0][1].to_vector());
                assert!(
                    normal.magnitude() > 1e-6,
                    "degenerate tangents at ({u}, {v})"
                );
                let sine =
                    radius.cross(normal).magnitude() / (radius.magnitude() * normal.magnitude());
                assert!(sine < 1e-9, "normal not radial at ({u}, {v}): sine {sine}");
            }
        }
    }

    #[test]
    fn uniform_weights_reduce_to_the_polynomial_surface() {
        let (ku, kv, g) = patch();
        let weighted = g.map(|p| Weighted {
            scaled: p.scale(2.0),
            weight: 2.0,
        });
        for iu in 0..=6 {
            for iv in 0..=6 {
                let (u, v) = (f64::from(iu) / 6.0, f64::from(iv) / 6.0);
                let plain = evaluate_surface(&ku, &kv, &g, u, v, T).unwrap();
                let rational = evaluate_rational_surface(&ku, &kv, &weighted, u, v, T).unwrap();
                assert!(plain.is_equal(rational, T));
            }
        }
    }
}
