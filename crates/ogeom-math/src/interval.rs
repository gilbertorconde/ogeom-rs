//! Interval arithmetic for filtered predicates.
//!
//! An [`Interval`] encloses a real number the program cannot represent: every
//! operation returns bounds that certainly contain the true result, obtained
//! by computing with the IEEE operations — which are correctly rounded, so
//! the true value lies within one ulp of the computed one — and then widening
//! each bound one step outward. The enclosure is conservative, never wrong.
//!
//! The point of carrying bounds is [`Interval::certain_sign`]: a sign
//! decision made through an interval is either *certain*, because zero lies
//! outside the bounds, or honestly undecided, because it does not. That is
//! the filter a predicate wants — answer fast when floating point can, and
//! say so when it cannot, instead of reading rounding noise as a direction.
//!
//! The vocabulary is the arithmetic predicates need: add, subtract, multiply,
//! negate, square, absolute value, square root, and division away from zero.
//! Transcendentals are deliberately absent — the standard library does not
//! state error bounds for them, and an enclosure that might not enclose is
//! worse than none.

use ogeom_core::predicates::Sign;

/// A closed interval `[lo, hi]` certainly containing a true real value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    lo: f64,
    hi: f64,
}

impl Interval {
    /// The degenerate interval holding one exactly-represented value.
    #[must_use]
    pub const fn point(value: f64) -> Self {
        Self {
            lo: value,
            hi: value,
        }
    }

    /// An interval from stated bounds. Panics if `lo > hi` or either is NaN,
    /// because such an "enclosure" encloses nothing.
    #[must_use]
    pub fn new(lo: f64, hi: f64) -> Self {
        assert!(lo <= hi, "an interval needs ordered, comparable bounds");
        Self { lo, hi }
    }

    /// A value with a symmetric absolute uncertainty.
    #[must_use]
    pub fn about(value: f64, radius: f64) -> Self {
        assert!(radius >= 0.0, "an uncertainty is a magnitude");
        Self::new(value - radius, value + radius)
    }

    /// The lower bound.
    #[must_use]
    pub const fn lo(&self) -> f64 {
        self.lo
    }

    /// The upper bound.
    #[must_use]
    pub const fn hi(&self) -> f64 {
        self.hi
    }

    /// The width of the enclosure.
    #[must_use]
    pub fn width(&self) -> f64 {
        self.hi - self.lo
    }

    /// Whether the enclosure contains `value`.
    #[must_use]
    pub fn contains(&self, value: f64) -> bool {
        self.lo <= value && value <= self.hi
    }

    /// The sign of the true value, where the bounds decide it: `None` means
    /// zero lies inside the enclosure and floating point genuinely cannot
    /// tell — which is an answer, not a failure.
    #[must_use]
    pub fn certain_sign(&self) -> Option<Sign> {
        if self.lo > 0.0 {
            Some(Sign::Positive)
        } else if self.hi < 0.0 {
            Some(Sign::Negative)
        } else if self.lo == 0.0 && self.hi == 0.0 {
            Some(Sign::Zero)
        } else {
            None
        }
    }

    /// The negation, exact — negation never rounds.
    #[must_use]
    pub const fn neg(&self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }

    /// The sum.
    #[must_use]
    pub fn add(&self, other: &Self) -> Self {
        Self {
            lo: (self.lo + other.lo).next_down(),
            hi: (self.hi + other.hi).next_up(),
        }
    }

    /// The difference.
    #[must_use]
    pub fn sub(&self, other: &Self) -> Self {
        Self {
            lo: (self.lo - other.hi).next_down(),
            hi: (self.hi - other.lo).next_up(),
        }
    }

    /// The product: the extremes over the four bound products, widened.
    #[must_use]
    pub fn mul(&self, other: &Self) -> Self {
        let products = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ];
        let mut lo = products[0];
        let mut hi = products[0];
        for p in &products[1..] {
            lo = lo.min(*p);
            hi = hi.max(*p);
        }
        Self {
            lo: lo.next_down(),
            hi: hi.next_up(),
        }
    }

    /// The square — tighter than `mul` with itself, because a square cannot
    /// be negative even when the interval straddles zero.
    #[must_use]
    pub fn square(&self) -> Self {
        let (a, b) = (self.lo * self.lo, self.hi * self.hi);
        if self.lo <= 0.0 && self.hi >= 0.0 {
            Self {
                lo: 0.0,
                hi: a.max(b).next_up(),
            }
        } else {
            Self {
                lo: a.min(b).next_down().max(0.0),
                hi: a.max(b).next_up(),
            }
        }
    }

    /// The absolute value, exact.
    #[must_use]
    pub fn abs(&self) -> Self {
        if self.lo >= 0.0 {
            *self
        } else if self.hi <= 0.0 {
            self.neg()
        } else {
            Self {
                lo: 0.0,
                hi: self.hi.max(-self.lo),
            }
        }
    }

    /// The square root, for enclosures of non-negative values. A lower bound
    /// pushed below zero by widening is clamped — the true value it encloses
    /// was non-negative. An interval entirely below zero has no real root
    /// and returns `None`.
    #[must_use]
    pub fn sqrt(&self) -> Option<Self> {
        if self.hi < 0.0 {
            return None;
        }
        let lo = if self.lo <= 0.0 {
            0.0
        } else {
            self.lo.sqrt().next_down().max(0.0)
        };
        Some(Self {
            lo,
            hi: self.hi.sqrt().next_up(),
        })
    }

    /// The quotient, defined only when the divisor certainly excludes zero.
    #[must_use]
    pub fn checked_div(&self, other: &Self) -> Option<Self> {
        if other.lo <= 0.0 && other.hi >= 0.0 {
            return None;
        }
        let quotients = [
            self.lo / other.lo,
            self.lo / other.hi,
            self.hi / other.lo,
            self.hi / other.hi,
        ];
        let mut lo = quotients[0];
        let mut hi = quotients[0];
        for q in &quotients[1..] {
            lo = lo.min(*q);
            hi = hi.max(*q);
        }
        Some(Self {
            lo: lo.next_down(),
            hi: hi.next_up(),
        })
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn finite() -> impl Strategy<Value = f64> {
        // Magnitudes a kernel actually computes with, both signs, zero included.
        prop_oneof![
            Just(0.0),
            -1e12..1e12f64,
            (-1.0..1.0f64).prop_map(|x| x * 1e-12),
        ]
    }

    proptest! {
        /// The defining law: an operation on point intervals encloses the
        /// floating-point result, which is within half an ulp of the truth.
        #[test]
        fn point_operations_enclose_their_own_result(a in finite(), b in finite()) {
            let (x, y) = (Interval::point(a), Interval::point(b));
            prop_assert!(x.add(&y).contains(a + b));
            prop_assert!(x.sub(&y).contains(a - b));
            prop_assert!(x.mul(&y).contains(a * b));
            prop_assert!(x.square().contains(a * a));
            prop_assert!(x.abs().contains(a.abs()));
            if a >= 0.0 {
                prop_assert!(x.sqrt().unwrap().contains(a.sqrt()));
            }
            if b != 0.0 {
                prop_assert!(x.checked_div(&y).unwrap().contains(a / b));
            }
        }

        /// Containment is monotone: whatever holds a value before an
        /// operation holds the operated value after.
        #[test]
        fn enclosures_stay_enclosures(a in finite(), b in finite(), r in 0.0..1e-6f64) {
            let x = Interval::about(a, r);
            let y = Interval::about(b, r);
            // The true values are a and b themselves; every combination of
            // them must land inside.
            prop_assert!(x.add(&y).contains(a + b));
            prop_assert!(x.mul(&y).contains(a * b));
            prop_assert!(x.sub(&y).contains(a - b));
        }
    }

    #[test]
    fn signs_are_certain_only_away_from_zero() {
        assert_eq!(
            Interval::new(1e-300, 2e-300).certain_sign(),
            Some(Sign::Positive)
        );
        assert_eq!(
            Interval::new(-2.0, -1e-300).certain_sign(),
            Some(Sign::Negative)
        );
        assert_eq!(Interval::point(0.0).certain_sign(), Some(Sign::Zero));
        assert_eq!(Interval::new(-1e-300, 1e-300).certain_sign(), None);
    }

    #[test]
    fn the_classic_rounding_case_is_enclosed() {
        // 0.1 + 0.2 in doubles is famously not 0.3; the enclosure holds the
        // computed sum and stays within a couple of ulps.
        let z = Interval::point(0.1).add(&Interval::point(0.2));
        assert!(z.contains(0.1 + 0.2));
        assert!(z.width() <= 4.0 * f64::EPSILON);
    }

    #[test]
    fn squares_of_straddling_intervals_start_at_zero() {
        let s = Interval::new(-2.0, 3.0).square();
        assert_eq!(s.lo(), 0.0);
        assert!(s.contains(9.0) && s.contains(0.25));
    }

    #[test]
    fn division_through_zero_refuses() {
        assert!(
            Interval::point(1.0)
                .checked_div(&Interval::new(-1.0, 1.0))
                .is_none()
        );
    }

    /// The use the type exists for: a cross-product magnitude computed from
    /// normals carrying a stated residual either certainly clears a floor or
    /// the interval says the arithmetic cannot decide.
    #[test]
    fn a_filtered_gate_decision() {
        let residual = 1e-12;
        let (ax, ay) = (
            Interval::about(1.0, residual),
            Interval::about(0.0, residual),
        );
        let (bx, by) = (
            Interval::about(1.0, residual),
            Interval::about(1e-6, residual),
        );
        // The 2D cross product ax*by - ay*bx encloses the true sine-scale.
        let cross = ax.mul(&by).sub(&ay.mul(&bx));
        assert_eq!(cross.certain_sign(), Some(Sign::Positive));
        // Shrink the angle to the residual's own scale and the sign is
        // honestly undecided.
        let by = Interval::about(1e-12, residual);
        let cross = ax.mul(&by).sub(&ay.mul(&bx));
        assert_eq!(cross.certain_sign(), None);
    }
}
