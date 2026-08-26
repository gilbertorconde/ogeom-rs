//! Knot vectors and B-spline basis functions.
//!
//! The basis is the foundation of every free-form curve and surface in the
//! kernel. Everything else — de Boor evaluation, knot insertion, degree
//! elevation, Bézier decomposition — is built on the functions here.
//!
//! # Representation
//!
//! A [`KnotVector`] stores the *flat* non-decreasing sequence, with repeated
//! knots written out. That is what every algorithm wants, and deriving it from
//! a distinct-knots-plus-multiplicities form on each call would cost an
//! allocation in the hottest loop in the crate.
//!
//! Repeated knots must be bit-identical, and every operation here preserves
//! that: knot insertion copies the inserted value rather than recomputing it.
//! Multiplicity is therefore an exact question, not a tolerance one, which
//! matters because multiplicity determines continuity — a knot of multiplicity
//! `p` in a degree-`p` curve is a corner, and "nearly a corner" is not a thing.
//!
//! # Conventions
//!
//! For degree `p` and `n` control points the flat vector has `n + p + 1`
//! entries. A *clamped* vector repeats its first and last knots `p + 1` times,
//! so the curve passes through its first and last control points; that is the
//! usual form for a bounded curve and the one [`KnotVector::clamped_uniform`]
//! produces.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use smallvec::SmallVec;

/// Basis values for one span, sized to avoid allocating for typical degrees.
pub type BasisValues = SmallVec<[f64; 8]>;

/// One row of basis values per derivative order, inline up to the jet
/// orders the kernel asks for.
pub type DerivativeRows = SmallVec<[BasisValues; 4]>;

/// A non-decreasing knot sequence with an associated degree.
#[derive(Debug, Clone, PartialEq)]
pub struct KnotVector {
    knots: Vec<f64>,
    degree: usize,
}

impl KnotVector {
    /// A knot vector from a flat non-decreasing sequence.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// sequence is too short for the degree, is not non-decreasing, contains a
    /// non-finite value, or has an interior knot of multiplicity greater than
    /// the degree — which would disconnect the curve rather than merely make it
    /// sharp.
    pub fn new(knots: Vec<f64>, degree: usize) -> OgeomResult<Self> {
        if degree == 0 {
            ogeom_bail!(Construction, "degree must be at least 1");
        }
        // n + p + 1 entries for n control points, and n >= p + 1 for the basis
        // to be well defined over a non-empty domain.
        let minimum = 2 * (degree + 1);
        if knots.len() < minimum {
            ogeom_bail!(
                Construction,
                "degree {degree} needs at least {minimum} knots, got {}",
                knots.len()
            );
        }
        if !knots.iter().all(|k| k.is_finite()) {
            ogeom_bail!(Construction, "knot vector contains a non-finite value");
        }
        if knots.windows(2).any(|w| w[1] < w[0]) {
            ogeom_bail!(Construction, "knot vector is not non-decreasing");
        }

        let this = Self { knots, degree };
        if this.domain_start() >= this.domain_end() {
            ogeom_bail!(Construction, "knot vector spans an empty domain");
        }

        // Interior multiplicity above the degree splits the curve in two.
        let (first, last) = (this.degree, this.knots.len() - this.degree - 1);
        let mut index = first;
        while index < last {
            let value = this.knots[index];
            let mut count = 0;
            while index < last && this.knots[index] == value {
                count += 1;
                index += 1;
            }
            // The two clamp knots at either end of the domain are allowed their
            // full multiplicity; only strictly interior ones are constrained.
            if value > this.domain_start() && value < this.domain_end() && count > this.degree {
                ogeom_bail!(
                    Construction,
                    "interior knot {value} has multiplicity {count}, above degree {}",
                    this.degree
                );
            }
        }
        Ok(this)
    }

    /// A clamped uniform knot vector for `control_points` control points.
    ///
    /// The domain is `[0, 1]`, the ends are clamped, and the interior knots are
    /// evenly spaced.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if there are
    /// too few control points for the degree.
    pub fn clamped_uniform(degree: usize, control_points: usize) -> OgeomResult<Self> {
        if control_points < degree + 1 {
            ogeom_bail!(
                Construction,
                "degree {degree} needs at least {} control points, got {control_points}",
                degree + 1
            );
        }
        let interior = control_points - degree - 1;
        let mut knots = Vec::with_capacity(control_points + degree + 1);
        knots.extend(core::iter::repeat_n(0.0, degree + 1));
        for i in 1..=interior {
            #[allow(clippy::cast_precision_loss)]
            knots.push(i as f64 / (interior + 1) as f64);
        }
        knots.extend(core::iter::repeat_n(1.0, degree + 1));
        Self::new(knots, degree)
    }

    /// A clamped knot vector from parameter values, for interpolation.
    ///
    /// Uses the averaging rule, which places interior knots so that the
    /// resulting interpolation system is well conditioned — a uniform vector
    /// over unevenly spaced parameters gives a nearly singular one.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if there are
    /// too few parameters, or they are not strictly increasing.
    pub fn averaged(degree: usize, parameters: &[f64]) -> OgeomResult<Self> {
        if parameters.len() < degree + 1 {
            ogeom_bail!(
                Construction,
                "degree {degree} needs at least {} parameters",
                degree + 1
            );
        }
        if parameters.windows(2).any(|w| w[1] <= w[0]) {
            ogeom_bail!(Construction, "parameters must be strictly increasing");
        }
        let n = parameters.len();
        let mut knots = Vec::with_capacity(n + degree + 1);
        knots.extend(core::iter::repeat_n(parameters[0], degree + 1));
        #[allow(clippy::cast_precision_loss)]
        for j in 1..n - degree {
            let mean: f64 = parameters[j..j + degree].iter().sum::<f64>() / degree as f64;
            knots.push(mean);
        }
        knots.extend(core::iter::repeat_n(parameters[n - 1], degree + 1));
        Self::new(knots, degree)
    }

    /// The degree.
    #[must_use]
    pub const fn degree(&self) -> usize {
        self.degree
    }

    /// The flat knot sequence.
    #[must_use]
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }

    /// The number of control points this vector describes.
    #[must_use]
    pub const fn control_point_count(&self) -> usize {
        self.knots.len() - self.degree - 1
    }

    /// The first parameter of the usable domain.
    #[must_use]
    pub fn domain_start(&self) -> f64 {
        self.knots[self.degree]
    }

    /// The last parameter of the usable domain.
    #[must_use]
    pub fn domain_end(&self) -> f64 {
        self.knots[self.knots.len() - self.degree - 1]
    }

    /// The usable domain.
    #[must_use]
    pub fn domain(&self) -> (f64, f64) {
        (self.domain_start(), self.domain_end())
    }

    /// Whether the ends are clamped, so the curve meets its first and last
    /// control points.
    #[must_use]
    pub fn is_clamped(&self) -> bool {
        let last = self.knots.len() - 1;
        self.knots[..=self.degree]
            .iter()
            .all(|k| *k == self.knots[0])
            && self.knots[last - self.degree..]
                .iter()
                .all(|k| *k == self.knots[last])
    }

    /// The multiplicity of the knot value at `index`.
    ///
    /// Exact: repeated knots are bit-identical by construction.
    #[must_use]
    pub fn multiplicity_at(&self, index: usize) -> usize {
        let Some(&value) = self.knots.get(index) else {
            return 0;
        };
        self.knots.iter().filter(|k| **k == value).count()
    }

    /// The multiplicity of `value`, or zero if it is not a knot.
    #[must_use]
    pub fn multiplicity_of(&self, value: f64) -> usize {
        self.knots.iter().filter(|k| **k == value).count()
    }

    /// The distinct knot values with their multiplicities, in order.
    #[must_use]
    pub fn distinct(&self) -> Vec<(f64, usize)> {
        let mut out: Vec<(f64, usize)> = Vec::new();
        for &k in &self.knots {
            match out.last_mut() {
                Some((value, count)) if *value == k => *count += 1,
                _ => out.push((k, 1)),
            }
        }
        out
    }

    /// Whether `u` lies in the usable domain, within `tol.parametric()`.
    #[must_use]
    pub fn contains(&self, u: f64, tol: Tolerances) -> bool {
        let (start, end) = self.domain();
        u >= start - tol.parametric() && u <= end + tol.parametric()
    }

    /// The index of the knot span containing `u`.
    ///
    /// Returns `i` with `knots[i] <= u < knots[i+1]`, clamped so that the end of
    /// the domain resolves to the last non-empty span rather than falling off
    /// it. Binary search, so cost is logarithmic in the knot count.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if `u` is outside the
    /// domain by more than `tol.parametric()`.
    pub fn span(&self, u: f64, tol: Tolerances) -> OgeomResult<usize> {
        let (start, end) = self.domain();
        if !u.is_finite() || u < start - tol.parametric() || u > end + tol.parametric() {
            ogeom_bail!(Domain, "parameter {u} outside knot domain [{start}, {end}]");
        }
        Ok(self.span_unchecked(u))
    }

    /// The knot span containing `u`, clamping out-of-range values into the
    /// domain rather than reporting them.
    #[must_use]
    pub fn span_unchecked(&self, u: f64) -> usize {
        let last = self.control_point_count() - 1;
        // The end of the domain belongs to the last span; without this it would
        // land one past it, since the search looks for `knots[i] <= u`.
        if u >= self.knots[last + 1] {
            return last;
        }
        if u <= self.knots[self.degree] {
            return self.degree;
        }
        let mut low = self.degree;
        let mut high = last + 1;
        while high - low > 1 {
            let mid = usize::midpoint(low, high);
            if u < self.knots[mid] {
                high = mid;
            } else {
                low = mid;
            }
        }
        low
    }

    /// The `degree + 1` non-zero basis functions at `u`.
    ///
    /// Entry `i` is the value of basis function `span - degree + i`. They are
    /// non-negative and sum to exactly one up to rounding — the partition of
    /// unity, which is what makes a B-spline curve lie in the convex hull of its
    /// control points.
    ///
    /// Cox-de Boor, in the triangular form that avoids evaluating the zero
    /// functions and never divides by a zero knot difference.
    #[must_use]
    pub fn basis(&self, span: usize, u: f64) -> BasisValues {
        let p = self.degree;
        let mut n = BasisValues::with_capacity(p + 1);
        n.push(1.0);
        let mut left = BasisValues::with_capacity(p + 1);
        let mut right = BasisValues::with_capacity(p + 1);
        left.push(0.0);
        right.push(0.0);

        for j in 1..=p {
            left.push(u - self.knots[span + 1 - j]);
            right.push(self.knots[span + j] - u);
            let mut saved = 0.0;
            n.push(0.0);
            for r in 0..j {
                // `right[r + 1] + left[j - r]` is the width of the union of two
                // adjacent supports, which is positive whenever the basis
                // function is, so this cannot divide by zero for a valid vector.
                let denominator = right[r + 1] + left[j - r];
                let temp = n[r] / denominator;
                n[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            n[j] = saved;
        }
        n
    }

    /// The non-zero basis functions and their derivatives up to order `n`.
    ///
    /// `result[k][i]` is the `k`th derivative of basis function
    /// `span - degree + i`. Orders above the degree are identically zero and
    /// are returned as such rather than as noise.
    #[must_use]
    pub fn basis_derivatives(&self, span: usize, u: f64, n: usize) -> DerivativeRows {
        let p = self.degree;
        let order = n.min(p);

        // `ndu` holds the basis values and the knot differences from the
        // triangular recurrence; both halves are needed to build derivatives.
        // Every scratch row lives inline for the degrees the kernel actually
        // meets: this is the innermost loop of every spline evaluation, and
        // it used to be the kernel's single largest allocation source.
        let mut ndu: SmallVec<[BasisValues; 8]> =
            core::iter::repeat_with(|| BasisValues::from_elem(0.0, p + 1))
                .take(p + 1)
                .collect();
        ndu[0][0] = 1.0;
        let mut left = BasisValues::from_elem(0.0, p + 1);
        let mut right = BasisValues::from_elem(0.0, p + 1);

        for j in 1..=p {
            left[j] = u - self.knots[span + 1 - j];
            right[j] = self.knots[span + j] - u;
            let mut saved = 0.0;
            for r in 0..j {
                ndu[j][r] = right[r + 1] + left[j - r];
                let temp = ndu[r][j - 1] / ndu[j][r];
                ndu[r][j] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            ndu[j][j] = saved;
        }

        let mut derivatives: DerivativeRows =
            core::iter::repeat_with(|| BasisValues::from_elem(0.0, p + 1))
                .take(n + 1)
                .collect();
        for (j, slot) in derivatives[0].iter_mut().enumerate() {
            *slot = ndu[j][p];
        }
        // The rows above `order` stay zero: a derivative past the degree of a
        // piecewise polynomial is identically zero, not merely small.

        // Two alternating rows of coefficients, per the standard algorithm.
        let mut a = [
            BasisValues::from_elem(0.0, p + 1),
            BasisValues::from_elem(0.0, p + 1),
        ];
        for r in 0..=p {
            let (mut s1, mut s2) = (0_usize, 1_usize);
            a[0][0] = 1.0;
            for k in 1..=order {
                let mut d = 0.0;
                let rk = r as isize - k as isize;
                let pk = p - k;
                if r >= k {
                    a[s2][0] = a[s1][0] / ndu[pk + 1][rk as usize];
                    d = a[s2][0] * ndu[rk as usize][pk];
                }
                let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
                let j2 = if r as isize - 1 <= pk as isize {
                    k - 1
                } else {
                    p - r
                };
                for j in j1..=j2 {
                    let index = (rk + j as isize) as usize;
                    a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][index];
                    d += a[s2][j] * ndu[index][pk];
                }
                if r <= pk {
                    a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                    d += a[s2][k] * ndu[r][pk];
                }
                derivatives[k][r] = d;
                core::mem::swap(&mut s1, &mut s2);
            }
        }

        // Multiply through by the falling factorial the recurrence omits.
        let mut factor = p;
        for (k, row) in derivatives.iter_mut().enumerate().take(order + 1).skip(1) {
            #[allow(clippy::cast_precision_loss)]
            let scale = factor as f64;
            for value in row.iter_mut() {
                *value *= scale;
            }
            factor = factor.saturating_mul(p.saturating_sub(k));
        }
        derivatives
    }

    /// Insert `value` into the sequence, `count` times.
    ///
    /// Only the knots change; adjusting control points to keep the shape is
    /// [`crate::bspline::insert_knot`].
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the result
    /// would push a knot's multiplicity above the degree.
    pub fn with_knot_inserted(&self, value: f64, count: usize) -> OgeomResult<Self> {
        let mut knots = self.knots.clone();
        let position = knots.partition_point(|k| *k <= value);
        for _ in 0..count {
            knots.insert(position, value);
        }
        Self::new(knots, self.degree)
    }

    /// This vector with its domain mapped onto `[start, end]`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the target
    /// range is empty or non-finite.
    pub fn reparameterized(&self, start: f64, end: f64) -> OgeomResult<Self> {
        if !start.is_finite() || !end.is_finite() || end <= start {
            ogeom_bail!(Construction, "target range [{start}, {end}] is empty");
        }
        let (a, b) = self.domain();
        let scale = (end - start) / (b - a);
        // Map, then overwrite the repeated end knots with the exact endpoints:
        // the arithmetic would otherwise give each copy a slightly different
        // value and silently destroy the clamping.
        let mut knots: Vec<f64> = self
            .knots
            .iter()
            .map(|k| (k - a).mul_add(scale, start))
            .collect();
        for k in &mut knots {
            if *k <= start {
                *k = start;
            } else if *k >= end {
                *k = end;
            }
        }
        let last = knots.len() - 1;
        for i in 0..self.knots.len() {
            if self.knots[i] == a {
                knots[i] = start;
            }
            if self.knots[last - i] == b {
                knots[last - i] = end;
            }
        }
        Self::new(knots, self.degree)
    }

    /// This vector with the parameter direction reversed.
    ///
    /// The domain is preserved and the sequence of interior spacings is
    /// mirrored. Reversing a curve reverses its knots and its control points
    /// together.
    ///
    /// Multiplicity is preserved *exactly* — equal knots map through the same
    /// arithmetic and so stay equal — which is what continuity depends on. The
    /// interior knot *values* are not bit-exactly restored by reversing twice,
    /// since `a + b - k` is not an exact involution in floating point; they
    /// return to within one ulp.
    #[must_use]
    pub fn reversed(&self) -> Self {
        let (a, b) = self.domain();
        let sum = a + b;
        let mut knots: Vec<f64> = self.knots.iter().rev().map(|k| sum - k).collect();
        // Same reasoning as `reparameterized`: restore the endpoints exactly.
        let last = knots.len() - 1;
        for i in 0..knots.len() {
            if knots[i] <= a {
                knots[i] = a;
            }
            if knots[last - i] >= b {
                knots[last - i] = b;
            }
        }
        Self {
            knots,
            degree: self.degree,
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

    fn cubic() -> KnotVector {
        // Degree 3, 7 control points, two interior knots.
        KnotVector::new(
            vec![0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0],
            3,
        )
        .unwrap()
    }

    #[test]
    fn malformed_vectors_are_refused() {
        assert!(KnotVector::new(vec![0.0, 1.0], 3).is_err(), "too short");
        assert!(
            KnotVector::new(vec![0.0, 0.0, 1.0, 0.5, 1.0, 1.0], 2).is_err(),
            "not non-decreasing"
        );
        assert!(
            KnotVector::new(vec![0.0, 0.0, f64::NAN, 1.0, 1.0, 1.0], 2).is_err(),
            "non-finite"
        );
        assert!(KnotVector::new(vec![0.0; 8], 3).is_err(), "empty domain");
        assert!(KnotVector::new(vec![], 0).is_err(), "degree zero");
        // Interior multiplicity above the degree disconnects the curve.
        assert!(KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0], 2).is_err());
        // At the degree it is merely a corner, which is legitimate.
        assert!(KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0], 2).is_ok());
    }

    #[test]
    fn clamped_uniform_has_the_expected_shape() {
        let k = KnotVector::clamped_uniform(3, 7).unwrap();
        assert_eq!(k.knots().len(), 11);
        assert_eq!(k.control_point_count(), 7);
        assert_eq!(k.domain(), (0.0, 1.0));
        assert!(k.is_clamped());
        assert_eq!(k.multiplicity_of(0.0), 4);
        assert_eq!(k.multiplicity_of(1.0), 4);
        assert_relative_eq!(k.knots()[4], 1.0 / 4.0);

        // A Bezier: no interior knots at all.
        let b = KnotVector::clamped_uniform(3, 4).unwrap();
        assert_eq!(b.knots(), &[0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]);
        assert!(KnotVector::clamped_uniform(3, 3).is_err());
    }

    #[test]
    fn span_lookup_handles_both_ends_of_the_domain() {
        let k = cubic();
        assert_eq!(k.span(0.0, T).unwrap(), 3, "domain start");
        assert_eq!(k.span(0.1, T).unwrap(), 3);
        assert_eq!(k.span(0.25, T).unwrap(), 4, "on a knot, span to its right");
        assert_eq!(k.span(0.3, T).unwrap(), 4);
        assert_eq!(k.span(0.9, T).unwrap(), 6);
        // The domain end must resolve to the last span, not one past it.
        assert_eq!(k.span(1.0, T).unwrap(), 6);
        assert!(k.span(-0.1, T).is_err());
        assert!(k.span(1.1, T).is_err());
        assert!(k.span(f64::NAN, T).is_err());
    }

    #[test]
    fn basis_functions_form_a_partition_of_unity() {
        let k = cubic();
        for i in 0..=200 {
            let u = f64::from(i) / 200.0;
            let span = k.span(u, T).unwrap();
            let n = k.basis(span, u);
            assert_eq!(n.len(), 4);
            let sum: f64 = n.iter().sum();
            assert_relative_eq!(sum, 1.0, epsilon = 1e-14);
            assert!(n.iter().all(|v| *v >= -1e-15), "basis must be non-negative");
        }
    }

    #[test]
    fn basis_matches_the_bernstein_polynomials_for_a_bezier() {
        // With no interior knots the B-spline basis is exactly Bernstein.
        let k = KnotVector::clamped_uniform(3, 4).unwrap();
        for i in 0..=20 {
            let u = f64::from(i) / 20.0;
            let span = k.span(u, T).unwrap();
            let n = k.basis(span, u);
            let v = 1.0 - u;
            let expected = [v * v * v, 3.0 * u * v * v, 3.0 * u * u * v, u * u * u];
            for (got, want) in n.iter().zip(expected) {
                assert_relative_eq!(*got, want, epsilon = 1e-14);
            }
        }
    }

    #[test]
    fn basis_is_an_interpolant_at_a_clamped_end() {
        let k = cubic();
        let n = k.basis(k.span(0.0, T).unwrap(), 0.0);
        assert_relative_eq!(n[0], 1.0, epsilon = 1e-15);
        assert!(n[1..].iter().all(|v| v.abs() < 1e-15));

        let n = k.basis(k.span(1.0, T).unwrap(), 1.0);
        assert_relative_eq!(n[3], 1.0, epsilon = 1e-15);
        assert!(n[..3].iter().all(|v| v.abs() < 1e-15));
    }

    #[test]
    fn basis_derivatives_agree_with_finite_differences() {
        let k = cubic();
        let h = 1e-6;
        for i in 1..20 {
            let u = f64::from(i) / 20.0;
            let span = k.span(u, T).unwrap();
            let d = k.basis_derivatives(span, u, 2);

            // Zeroth order must reproduce the plain basis.
            let plain = k.basis(span, u);
            for (a, b) in d[0].iter().zip(plain.iter()) {
                assert_relative_eq!(a, b, epsilon = 1e-14);
            }

            // First order against a central difference, evaluated in the same
            // span so the basis indices line up.
            let ahead = k.basis(span, u + h);
            let behind = k.basis(span, u - h);
            for j in 0..=k.degree() {
                let numeric = (ahead[j] - behind[j]) / (2.0 * h);
                assert_relative_eq!(d[1][j], numeric, epsilon = 1e-5);
            }
        }
    }

    #[test]
    fn basis_derivatives_sum_to_zero() {
        // The basis sums to one everywhere, so every derivative of that sum is
        // identically zero — a strong check on the whole recurrence.
        let k = cubic();
        for i in 0..=50 {
            let u = f64::from(i) / 50.0;
            let span = k.span(u, T).unwrap();
            let d = k.basis_derivatives(span, u, 3);
            assert_relative_eq!(d[0].iter().sum::<f64>(), 1.0, epsilon = 1e-13);
            for (order, row) in d.iter().enumerate().skip(1) {
                let sum: f64 = row.iter().sum();
                assert!(sum.abs() < 1e-8, "order {order} sums to {sum}");
            }
        }
    }

    #[test]
    fn derivatives_above_the_degree_are_zero() {
        let k = cubic();
        let span = k.span(0.4, T).unwrap();
        let d = k.basis_derivatives(span, 0.4, 5);
        assert_eq!(d.len(), 6);
        for (order, row) in d.iter().enumerate().skip(k.degree() + 1) {
            assert!(row.iter().all(|v| *v == 0.0), "order {order} is not zero");
        }
    }

    #[test]
    fn multiplicity_and_distinct_knots() {
        let k = cubic();
        assert_eq!(k.multiplicity_of(0.0), 4);
        assert_eq!(k.multiplicity_of(0.5), 1);
        assert_eq!(k.multiplicity_of(0.6), 0);
        assert_eq!(k.multiplicity_at(0), 4);
        assert_eq!(
            k.distinct(),
            vec![(0.0, 4), (0.25, 1), (0.5, 1), (0.75, 1), (1.0, 4)]
        );
    }

    #[test]
    fn knot_insertion_raises_multiplicity() {
        let k = cubic();
        let inserted = k.with_knot_inserted(0.5, 2).unwrap();
        assert_eq!(inserted.multiplicity_of(0.5), 3);
        assert_eq!(inserted.knots().len(), k.knots().len() + 2);
        assert_eq!(inserted.domain(), k.domain());
        // Beyond the degree it would disconnect the curve.
        assert!(k.with_knot_inserted(0.5, 3).is_err());
    }

    #[test]
    fn reparameterization_preserves_clamping_exactly() {
        let k = cubic();
        let r = k.reparameterized(-2.0, 6.0).unwrap();
        assert_eq!(r.domain(), (-2.0, 6.0));
        assert!(r.is_clamped(), "the repeated end knots must stay identical");
        assert_eq!(r.multiplicity_of(-2.0), 4);
        assert_eq!(r.multiplicity_of(6.0), 4);
        // Interior knots map proportionally.
        assert_relative_eq!(r.knots()[4], 0.0, epsilon = 1e-12);
        assert!(k.reparameterized(1.0, 1.0).is_err());
        assert!(k.reparameterized(1.0, f64::NAN).is_err());
    }

    #[test]
    fn reversal_mirrors_the_spacing_and_keeps_the_domain() {
        // Deliberately uneven interior spacing, so a mirror is observable.
        let k = KnotVector::new(vec![0.0, 0.0, 0.0, 0.1, 0.8, 1.0, 1.0, 1.0], 2).unwrap();
        let r = k.reversed();
        assert_eq!(r.domain(), k.domain());
        assert!(r.is_clamped());
        assert_relative_eq!(r.knots()[3], 0.2, epsilon = 1e-15);
        assert_relative_eq!(r.knots()[4], 0.9, epsilon = 1e-15);
        // Reversing twice restores the vector to within rounding. Not exactly:
        // `a + b - k` is not an exact involution in floating point.
        for (got, want) in r.reversed().knots().iter().zip(k.knots()) {
            assert_relative_eq!(got, want, epsilon = 1e-15);
        }
        // What must hold exactly is multiplicity, since continuity depends on
        // it: two knots that were equal stay equal through any number of
        // reversals.
        let multiplicities: Vec<usize> = r.distinct().iter().map(|(_, m)| *m).collect();
        let original: Vec<usize> = k.distinct().iter().map(|(_, m)| *m).collect();
        assert_eq!(multiplicities, original);
    }

    #[test]
    fn averaged_knots_follow_the_parameters() {
        let params = [0.0, 0.1, 0.4, 0.9, 1.0];
        let k = KnotVector::averaged(3, &params).unwrap();
        assert_eq!(k.control_point_count(), 5);
        assert_eq!(k.domain(), (0.0, 1.0));
        assert!(k.is_clamped());
        // One interior knot, the mean of parameters 1..4.
        assert_relative_eq!(k.knots()[4], (0.1 + 0.4 + 0.9) / 3.0, epsilon = 1e-15);
        assert!(KnotVector::averaged(3, &[0.0, 1.0]).is_err());
        assert!(
            KnotVector::averaged(2, &[0.0, 0.5, 0.5, 1.0]).is_err(),
            "not increasing"
        );
    }

    #[test]
    fn basis_at_a_repeated_interior_knot_is_still_a_partition_of_unity() {
        // Multiplicity equal to the degree: a corner, where the recurrence has
        // the most opportunity to divide by something vanishing.
        let k = KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0], 2).unwrap();
        for u in [0.0_f64, 0.25, 0.499_999, 0.5, 0.500_001, 0.75, 1.0] {
            let span = k.span(u, T).unwrap();
            let n = k.basis(span, u);
            assert_relative_eq!(n.iter().sum::<f64>(), 1.0, epsilon = 1e-14);
            assert!(n.iter().all(|v| v.is_finite()), "non-finite basis at {u}");
        }
    }
}
