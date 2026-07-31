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

use og_core::{OgResult, Tolerances, og_bail};

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
    /// [`OgError::Construction`](og_core::OgError::Construction) if the weight
    /// is not finite and positive. A zero weight makes the projection undefined
    /// and a negative one makes the curve leave its control polygon's convex
    /// hull, so neither is admitted.
    pub fn new(point: P, weight: f64, tol: Tolerances) -> OgResult<Self> {
        if !weight.is_finite() || weight <= tol.confusion() {
            og_bail!(
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
fn check_shape<P>(knots: &KnotVector, control: &[P]) -> OgResult<()> {
    if control.len() != knots.control_point_count() {
        og_bail!(
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
/// [`OgError::Dimension`](og_core::OgError::Dimension) if the control point
/// count disagrees with the knot vector; [`OgError::Domain`](og_core::OgError::Domain)
/// if `u` is outside the domain.
pub fn evaluate<P: Blend>(
    knots: &KnotVector,
    control: &[P],
    u: f64,
    tol: Tolerances,
) -> OgResult<P> {
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
) -> OgResult<Vec<P>> {
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
/// [`OgError::Dimension`](og_core::OgError::Dimension) on a shape mismatch,
/// [`OgError::Domain`](og_core::OgError::Domain) if `value` is outside the
/// domain, and [`OgError::Construction`](og_core::OgError::Construction) if the
/// insertion would push a multiplicity above the degree.
pub fn insert_knot<P: Blend>(
    knots: &KnotVector,
    control: &[P],
    value: f64,
    count: usize,
    tol: Tolerances,
) -> OgResult<Spline<P>> {
    check_shape(knots, control)?;
    if count == 0 {
        return Ok((knots.clone(), control.to_vec()));
    }
    let span = knots.span(value, tol)?;
    let p = knots.degree();
    let existing = knots.multiplicity_of(value);
    if existing + count > p {
        og_bail!(
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
/// As [`insert_knot`], plus [`OgError::Domain`](og_core::OgError::Domain) if
/// `u` is at either end of the domain, where one half would be empty.
pub fn split<P: Blend>(
    knots: &KnotVector,
    control: &[P],
    u: f64,
    tol: Tolerances,
) -> OgResult<(Spline<P>, Spline<P>)> {
    check_shape(knots, control)?;
    let (start, end) = knots.domain();
    if u <= start + tol.parametric() || u >= end - tol.parametric() {
        og_bail!(
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
) -> OgResult<Vec<BezierSegment<P>>> {
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
) -> OgResult<Spline<P>> {
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
/// As [`evaluate`], plus [`OgError::Numeric`](og_core::OgError::Numeric) if the
/// accumulated weight vanishes, which positive input weights make impossible.
pub fn evaluate_rational<P: Blend>(
    knots: &KnotVector,
    control: &[Weighted<P>],
    u: f64,
    tol: Tolerances,
) -> OgResult<P> {
    let h = evaluate(knots, control, u, tol)?;
    if h.weight.abs() <= tol.confusion() {
        og_bail!(Numeric, "rational evaluation produced a vanishing weight");
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
) -> OgResult<Vec<P>> {
    let homogeneous = derivatives(knots, control, u, n, tol)?;
    if homogeneous[0].weight.abs() <= tol.confusion() {
        og_bail!(Numeric, "rational evaluation produced a vanishing weight");
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
