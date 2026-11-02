//! The two number types every predicate is evaluated in: an interval that
//! decides most signs in floating point, and an exact dyadic number that
//! decides the rest.
//!
//! A predicate is one polynomial, written once over [`Ring`]. It is first
//! evaluated over [`Iv`]; when zero lies inside the enclosure, the sign is
//! not yet known and the same polynomial is evaluated over [`Big`], whose
//! arithmetic never rounds. Every input coordinate is an `f64`, which is a
//! dyadic rational, so sums and products of them are dyadic rationals too
//! and [`Big`] holds each one exactly.

use ogeom_core::Sign;

/// The arithmetic a predicate's polynomial needs.
pub(super) trait Ring: Clone {
    /// An input coordinate, exactly.
    fn of(value: f64) -> Self;
    fn add(&self, other: &Self) -> Self;
    fn sub(&self, other: &Self) -> Self;
    fn mul(&self, other: &Self) -> Self;
    fn neg(&self) -> Self;
    /// The sign, or `None` while the value is not known closely enough.
    fn sign(&self) -> Option<Sign>;
}

/// A closed interval certainly holding a real value.
///
/// Rounded outward on every operation: a bound is moved one step away only
/// when the operation rounded, so values computed without rounding stay
/// points, and an exact zero stays a certain zero. A bound that overflows
/// is clamped to the largest finite value on its side, and a bound that
/// cannot be told at all opens the interval to the whole line.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Iv {
    pub lo: f64,
    pub hi: f64,
}

impl Iv {
    const WHOLE: Self = Self {
        lo: f64::NEG_INFINITY,
        hi: f64::INFINITY,
    };

    fn bounded(lo: f64, hi: f64) -> Self {
        if lo.is_nan() || hi.is_nan() {
            return Self::WHOLE;
        }
        Self {
            lo: if lo == f64::INFINITY { f64::MAX } else { lo },
            hi: if hi == f64::NEG_INFINITY {
                -f64::MAX
            } else {
                hi
            },
        }
    }

    const fn is_zero(&self) -> bool {
        self.lo == 0.0 && self.hi == 0.0
    }

    /// The quotient `self / den`, enclosed, when `den` is certainly not
    /// zero.
    pub fn over(&self, den: &Self) -> Option<Self> {
        if !(den.lo > 0.0 || den.hi < 0.0) {
            return None;
        }
        if self.lo == self.hi && den.lo == den.hi && self.lo.abs() >= EXACT_ERROR_FLOOR {
            // A zero residual, exact by the fused multiply-add this far
            // above underflow, says the quotient did not round.
            let q = self.lo / den.lo;
            if q.is_finite() && q.mul_add(den.lo, -self.lo) == 0.0 {
                return Some(Self::of(q));
            }
        }
        let q = [
            self.lo / den.lo,
            self.lo / den.hi,
            self.hi / den.lo,
            self.hi / den.hi,
        ];
        if q.iter().any(|v| !v.is_finite()) {
            return None;
        }
        let lo = q.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = q.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Some(Self::bounded(lo.next_down(), hi.next_up()))
    }

    /// The value nearest the middle, and the farthest the true value can lie
    /// from it.
    pub fn middle(&self) -> (f64, f64) {
        if self.lo == self.hi {
            return (self.lo, 0.0);
        }
        let mid = self.lo / 2.0 + self.hi / 2.0;
        (mid, (mid - self.lo).max(self.hi - mid).next_up())
    }
}

/// The least product magnitude whose rounding error is itself a normal
/// `f64`: the error of a product of `2^e` is a multiple of `2^(e - 105)`,
/// which stays above the subnormal range while `e >= -900`.
const EXACT_ERROR_FLOOR: f64 = f64::from_bits((1023 - 900) << 52);

/// `x + y` and the rounding error's sign: negative when the true sum lies
/// below the rounded one.
fn sum(x: f64, y: f64) -> (f64, f64) {
    let s = x + y;
    let bv = s - x;
    let av = s - bv;
    (s, (x - av) + (y - bv))
}

fn sum_down(x: f64, y: f64) -> f64 {
    let (s, err) = sum(x, y);
    if err < 0.0 { s.next_down() } else { s }
}

fn sum_up(x: f64, y: f64) -> f64 {
    let (s, err) = sum(x, y);
    if err > 0.0 { s.next_up() } else { s }
}

impl Ring for Iv {
    fn of(value: f64) -> Self {
        Self {
            lo: value,
            hi: value,
        }
    }

    fn add(&self, other: &Self) -> Self {
        Self::bounded(sum_down(self.lo, other.lo), sum_up(self.hi, other.hi))
    }

    fn sub(&self, other: &Self) -> Self {
        Self::bounded(sum_down(self.lo, -other.hi), sum_up(self.hi, -other.lo))
    }

    fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::of(0.0);
        }
        if self.lo == self.hi && other.lo == other.hi {
            let p = self.lo * other.lo;
            let err = self.lo.mul_add(other.lo, -p);
            if p.is_finite() && p.abs() >= EXACT_ERROR_FLOOR {
                // The rounding error of a product, exact by the fused
                // multiply-add while the product is this far above
                // underflow, says which way the true product lies.
                let lo = if err < 0.0 { p.next_down() } else { p };
                let hi = if err > 0.0 { p.next_up() } else { p };
                return Self::bounded(lo, hi);
            }
        }
        let p = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ];
        if p.iter().any(|v| v.is_nan()) {
            return Self::WHOLE;
        }
        let lo = p.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = p.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Self::bounded(lo.next_down(), hi.next_up())
    }

    fn neg(&self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }

    fn sign(&self) -> Option<Sign> {
        if self.lo > 0.0 {
            Some(Sign::Positive)
        } else if self.hi < 0.0 {
            Some(Sign::Negative)
        } else if self.is_zero() {
            Some(Sign::Zero)
        } else {
            None
        }
    }
}

/// An exact dyadic rational: `±mag · 2^exp`, the magnitude in 64-bit limbs,
/// least significant first, with no zero limb at either end. Zero has no
/// limbs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Big {
    neg: bool,
    mag: Vec<u64>,
    exp: i64,
}

impl Big {
    fn zero() -> Self {
        Self {
            neg: false,
            mag: Vec::new(),
            exp: 0,
        }
    }

    fn trimmed(neg: bool, mut mag: Vec<u64>, mut exp: i64) -> Self {
        while mag.last() == Some(&0) {
            mag.pop();
        }
        let low = mag.iter().take_while(|&&l| l == 0).count();
        if low == mag.len() {
            return Self::zero();
        }
        if low > 0 {
            mag.drain(..low);
            exp += 64 * i64::try_from(low).unwrap_or(i64::MAX / 128);
        }
        Self { neg, mag, exp }
    }

    /// The magnitude shifted left by `bits`.
    fn shifted(mag: &[u64], bits: u64) -> Vec<u64> {
        let limbs = usize::try_from(bits / 64).unwrap_or(usize::MAX);
        let rest = bits % 64;
        let mut out = vec![0_u64; limbs];
        out.reserve(mag.len() + 1);
        if rest == 0 {
            out.extend_from_slice(mag);
        } else {
            let mut carry = 0_u64;
            for &l in mag {
                out.push((l << rest) | carry);
                carry = l >> (64 - rest);
            }
            out.push(carry);
        }
        out
    }

    fn cmp_mag(a: &[u64], b: &[u64]) -> std::cmp::Ordering {
        let (a, b) = (strip(a), strip(b));
        a.len()
            .cmp(&b.len())
            .then_with(|| a.iter().rev().cmp(b.iter().rev()))
    }

    fn add_mag(a: &[u64], b: &[u64]) -> Vec<u64> {
        let (long, short) = if a.len() >= b.len() { (a, b) } else { (b, a) };
        let mut out = Vec::with_capacity(long.len() + 1);
        let mut carry = false;
        for (i, &l) in long.iter().enumerate() {
            let (s, c1) = l.overflowing_add(short.get(i).copied().unwrap_or(0));
            let (s, c2) = s.overflowing_add(u64::from(carry));
            out.push(s);
            carry = c1 || c2;
        }
        out.push(u64::from(carry));
        out
    }

    /// `a - b` for `a >= b`.
    fn sub_mag(a: &[u64], b: &[u64]) -> Vec<u64> {
        let mut out = Vec::with_capacity(a.len());
        let mut borrow = false;
        for (i, &l) in a.iter().enumerate() {
            let (s, b1) = l.overflowing_sub(b.get(i).copied().unwrap_or(0));
            let (s, b2) = s.overflowing_sub(u64::from(borrow));
            out.push(s);
            borrow = b1 || b2;
        }
        out
    }

    /// A tight enclosure of the value in floating point.
    pub fn enclose(&self) -> Iv {
        let n = self.mag.len();
        if n == 0 {
            return Iv::of(0.0);
        }
        let take = n.min(2);
        let mut top = 0_u128;
        for &l in self.mag[n - take..].iter().rev() {
            top = (top << 64) | u128::from(l);
        }
        let shift = self.exp + 64 * i64::try_from(n - take).unwrap_or(i64::MAX / 128);
        #[allow(clippy::cast_precision_loss)]
        let approx = scale(top as f64, shift);
        let (lo, hi) = (approx.next_down().next_down(), approx.next_up().next_up());
        let lo = if lo.is_nan() { 0.0 } else { lo.max(0.0) };
        let hi = if hi.is_nan() { f64::INFINITY } else { hi };
        if self.neg {
            Iv::bounded(-hi, -lo)
        } else {
            Iv::bounded(lo, hi)
        }
    }
}

fn strip(a: &[u64]) -> &[u64] {
    let mut n = a.len();
    while n > 0 && a[n - 1] == 0 {
        n -= 1;
    }
    &a[..n]
}

/// `value · 2^shift`, exact while the result stays normal.
fn scale(mut value: f64, mut shift: i64) -> f64 {
    let step = |s: i64| f64::from_bits(u64::try_from(1023 + s).unwrap_or(0) << 52);
    while shift > 1000 {
        value *= step(1000);
        shift -= 1000;
    }
    while shift < -1000 {
        value *= step(-1000);
        shift += 1000;
    }
    value * step(shift)
}

impl Ring for Big {
    fn of(value: f64) -> Self {
        if value == 0.0 || !value.is_finite() {
            return Self::zero();
        }
        let bits = value.to_bits();
        let biased = i64::try_from((bits >> 52) & 0x7ff).unwrap_or(0);
        let frac = bits & ((1_u64 << 52) - 1);
        let (mut mant, mut exp) = if biased == 0 {
            (frac, -1074)
        } else {
            (frac | (1_u64 << 52), biased - 1075)
        };
        let tz = mant.trailing_zeros();
        mant >>= tz;
        exp += i64::from(tz);
        Self {
            neg: value < 0.0,
            mag: vec![mant],
            exp,
        }
    }

    fn add(&self, other: &Self) -> Self {
        if self.mag.is_empty() {
            return other.clone();
        }
        if other.mag.is_empty() {
            return self.clone();
        }
        let exp = self.exp.min(other.exp);
        let a = Self::shifted(&self.mag, self.exp.abs_diff(exp));
        let b = Self::shifted(&other.mag, other.exp.abs_diff(exp));
        if self.neg == other.neg {
            return Self::trimmed(self.neg, Self::add_mag(&a, &b), exp);
        }
        match Self::cmp_mag(&a, &b) {
            std::cmp::Ordering::Equal => Self::zero(),
            std::cmp::Ordering::Greater => Self::trimmed(self.neg, Self::sub_mag(&a, &b), exp),
            std::cmp::Ordering::Less => Self::trimmed(other.neg, Self::sub_mag(&b, &a), exp),
        }
    }

    fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    fn mul(&self, other: &Self) -> Self {
        if self.mag.is_empty() || other.mag.is_empty() {
            return Self::zero();
        }
        let mut out = vec![0_u64; self.mag.len() + other.mag.len()];
        for (i, &x) in self.mag.iter().enumerate() {
            let mut carry = 0_u128;
            for (j, &y) in other.mag.iter().enumerate() {
                let t = u128::from(x) * u128::from(y) + u128::from(out[i + j]) + carry;
                #[allow(clippy::cast_possible_truncation)]
                {
                    out[i + j] = t as u64;
                }
                carry = t >> 64;
            }
            let mut k = i + other.mag.len();
            while carry > 0 {
                let t = u128::from(out[k]) + carry;
                #[allow(clippy::cast_possible_truncation)]
                {
                    out[k] = t as u64;
                }
                carry = t >> 64;
                k += 1;
            }
        }
        Self::trimmed(self.neg != other.neg, out, self.exp + other.exp)
    }

    fn neg(&self) -> Self {
        if self.mag.is_empty() {
            return self.clone();
        }
        Self {
            neg: !self.neg,
            mag: self.mag.clone(),
            exp: self.exp,
        }
    }

    fn sign(&self) -> Option<Sign> {
        Some(if self.mag.is_empty() {
            Sign::Zero
        } else if self.neg {
            Sign::Negative
        } else {
            Sign::Positive
        })
    }
}

/// The product of two signs.
pub(super) const fn times(a: Sign, b: Sign) -> Sign {
    match (a, b) {
        (Sign::Zero, _) | (_, Sign::Zero) => Sign::Zero,
        (Sign::Positive, s) => s,
        (Sign::Negative, s) => s.reversed(),
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn exact_sum(values: &[f64]) -> Big {
        values
            .iter()
            .fold(Big::zero(), |acc, &v| acc.add(&Big::of(v)))
    }

    #[test]
    fn big_cancels_what_floating_point_loses() {
        // 1e30 + 1 - 1e30 is 1 exactly, and 0 in floating point.
        let s = exact_sum(&[1e30, 1.0, -1e30]);
        assert_eq!(s, Big::of(1.0));
        let p = Big::of(0.1).mul(&Big::of(3.0)).sub(&Big::of(0.3));
        // 0.1 and 0.3 are not tenths; 3 · fl(0.1) - fl(0.3) is a tiny
        // positive dyadic, and the interval cannot tell its sign.
        assert_eq!(p.sign(), Some(Sign::Positive));
        let iv = Iv::of(0.1).mul(&Iv::of(3.0)).sub(&Iv::of(0.3));
        assert_eq!(iv.sign(), None);
    }

    #[test]
    fn interval_keeps_exact_zeros_and_encloses_the_rest() {
        let z = Iv::of(2.5).sub(&Iv::of(2.5));
        assert_eq!(z.sign(), Some(Sign::Zero));
        let big = Big::of(1e-300).mul(&Big::of(1e-300)).add(&Big::of(-1e-310));
        let iv = Iv::of(1e-300).mul(&Iv::of(1e-300)).add(&Iv::of(-1e-310));
        let e = big.enclose();
        assert!(iv.lo <= e.hi && e.lo <= iv.hi);
        assert_eq!(big.sign(), Some(Sign::Negative));
        let huge = Iv::of(f64::MAX).add(&Iv::of(f64::MAX));
        assert_eq!(huge.sign(), Some(Sign::Positive));
        let nan = Iv::of(f64::MAX).mul(&Iv::of(f64::MAX)).sub(&huge);
        assert_eq!(nan.sign(), None);
    }

    #[test]
    fn enclosure_holds_the_value() {
        for &(a, b) in &[(1.5, 3.25), (-7.0, 1e-20), (1e100, 3e150), (0.1, 0.7)] {
            let v = Big::of(a).mul(&Big::of(b)).add(&Big::of(a));
            let e = v.enclose();
            let f = a.mul_add(b, a);
            assert!(e.lo <= f.next_up() && f.next_down() <= e.hi, "{a} {b}");
            assert!(e.hi - e.lo <= 8.0 * f.abs() * f64::EPSILON);
        }
    }
}
