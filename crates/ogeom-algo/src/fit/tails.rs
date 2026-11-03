//! The two-sided tails of the normal and the Student t distributions, and
//! the t quantile at a normal tail: how far a residual must lie, in units
//! of a noise estimated from a few points, to be as rare as `k` standard
//! deviations of a noise that is known. Likewise the F quantile at a normal
//! tail, for a ratio of two such estimates.

/// The most terms a continued fraction or a series sums; each converges
/// geometrically well before this for every argument used here.
const MOST_TERMS: usize = 100_000;

/// A continued fraction's terms are summed until one changes it by less
/// than this, relatively.
const SETTLED: f64 = 1e-15;

/// The floor a continued fraction's denominators are kept from zero by.
const TINY: f64 = 1e-300;

/// `ln Γ(x)` for `x > 0`, by Lanczos's approximation with `g = 7` and
/// nine coefficients, accurate to about fifteen digits.
fn ln_gamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // The reflection formula.
        let pi = core::f64::consts::PI;
        return (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let z = x - 1.0;
    let mut sum = C[0];
    for (i, c) in C.iter().enumerate().skip(1) {
        #[allow(clippy::cast_precision_loss, reason = "nine terms")]
        let i = i as f64;
        sum += c / (z + i);
    }
    let t = z + G + 0.5;
    0.5f64.mul_add((2.0 * core::f64::consts::PI).ln(), (z + 0.5) * t.ln()) - t + sum.ln()
}

/// The regularised upper incomplete gamma function `Q(a, x)`.
fn gamma_q(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    let front = (a * x.ln() - x - ln_gamma(a)).exp();
    if x < a + 1.0 {
        // The series for P(a, x).
        let (mut term, mut sum, mut n) = (1.0 / a, 1.0 / a, a);
        for _ in 0..MOST_TERMS {
            n += 1.0;
            term *= x / n;
            sum += term;
            if term.abs() < sum.abs() * SETTLED {
                break;
            }
        }
        return 1.0 - front * sum;
    }
    // The continued fraction for Q(a, x), by Lentz's method.
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / TINY;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..MOST_TERMS {
        #[allow(clippy::cast_precision_loss, reason = "far below 2^52 terms")]
        let i = i as f64;
        let an = -i * (i - a);
        b += 2.0;
        d = an.mul_add(d, b);
        if d.abs() < TINY {
            d = TINY;
        }
        c = b + an / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let step = d * c;
        h *= step;
        if (step - 1.0).abs() < SETTLED {
            break;
        }
    }
    front * h
}

/// The continued fraction of the incomplete beta function, by Lentz's
/// method.
fn beta_fraction(a: f64, b: f64, x: f64) -> f64 {
    let guard = |v: f64| if v.abs() < TINY { TINY } else { v };
    let mut c = 1.0;
    let mut d = 1.0 / guard(1.0 - (a + b) * x / (a + 1.0));
    let mut h = d;
    for m in 1..MOST_TERMS {
        #[allow(clippy::cast_precision_loss, reason = "far below 2^52 terms")]
        let m = m as f64;
        let m2 = 2.0 * m;
        let even = m * (b - m) * x / ((a + m2 - 1.0) * (a + m2));
        d = 1.0 / guard(even.mul_add(d, 1.0));
        c = guard(1.0 + even / c);
        h *= d * c;
        let odd = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1.0));
        d = 1.0 / guard(odd.mul_add(d, 1.0));
        c = guard(1.0 + odd / c);
        let step = d * c;
        h *= step;
        if (step - 1.0).abs() < SETTLED {
            break;
        }
    }
    h
}

/// The regularised incomplete beta function `I_x(a, b)`.
fn beta_i(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let front = (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (-x).ln_1p()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        front * beta_fraction(a, b, x) / a
    } else {
        1.0 - front * beta_fraction(b, a, 1.0 - x) / b
    }
}

/// `P(|Z| > k)` for a standard normal `Z`.
pub(super) fn normal_tail(k: f64) -> f64 {
    gamma_q(0.5, 0.5 * k * k)
}

/// The `k` at which the normal's two-sided tail is `tail`: zero for a tail
/// of one or more, infinite for none.
pub(super) fn normal_at_tail(tail: f64) -> f64 {
    if tail.is_nan() || tail <= 0.0 {
        return f64::INFINITY;
    }
    if tail >= 1.0 {
        return 0.0;
    }
    // The tail falls as k grows: double past it, then halve the bracket
    // until it cannot narrow.
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    while normal_tail(hi) > tail {
        lo = hi;
        hi *= 2.0;
        if !hi.is_finite() {
            return f64::INFINITY;
        }
    }
    loop {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            return hi;
        }
        if normal_tail(mid) > tail {
            lo = mid;
        } else {
            hi = mid;
        }
    }
}

/// `P(|T| > t)` for Student's `T` with `dof` degrees of freedom.
fn t_tail(t: f64, dof: f64) -> f64 {
    beta_i(0.5 * dof, 0.5, dof / t.mul_add(t, dof))
}

/// The `t` at which Student's two-sided tail with `dof` degrees of freedom
/// is the normal's at `k`: `k` itself in the limit of many, and larger the
/// fewer there are. Infinite with none, and for a `k` with no tail.
pub(super) fn t_at_normal(k: f64, dof: usize) -> f64 {
    if dof == 0 || k.is_infinite() {
        return f64::INFINITY;
    }
    let tail = normal_tail(k);
    if tail.is_nan() || tail <= 0.0 {
        return f64::INFINITY;
    }
    #[allow(clippy::cast_precision_loss, reason = "point counts")]
    let dof = dof as f64;
    // The t tail falls as t grows: double past it, then halve the bracket
    // until it cannot narrow.
    let (mut lo, mut hi) = (0.0_f64, k.max(1.0));
    while t_tail(hi, dof) > tail {
        lo = hi;
        hi *= 2.0;
        if !hi.is_finite() {
            return f64::INFINITY;
        }
    }
    loop {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            return hi;
        }
        if t_tail(mid, dof) > tail {
            lo = mid;
        } else {
            hi = mid;
        }
    }
}

/// `P(F > f)` for Snedecor's `F` with `d1` and `d2` degrees of freedom.
fn f_tail(f: f64, d1: f64, d2: f64) -> f64 {
    beta_i(0.5 * d2, 0.5 * d1, d2 / d1.mul_add(f, d2))
}

/// The `f` at which the upper tail of Snedecor's `F` with `d1` and `d2`
/// degrees of freedom is the normal's two-sided tail at `k`. Infinite with
/// either count zero, and for a `k` with no tail.
pub(super) fn f_at_normal(k: f64, d1: usize, d2: usize) -> f64 {
    if d1 == 0 || d2 == 0 || k.is_infinite() {
        return f64::INFINITY;
    }
    let tail = normal_tail(k);
    if tail.is_nan() || tail <= 0.0 {
        return f64::INFINITY;
    }
    #[allow(clippy::cast_precision_loss, reason = "point counts")]
    let (d1, d2) = (d1 as f64, d2 as f64);
    // The tail falls as f grows: double past it, then halve the bracket
    // until it cannot narrow.
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    while f_tail(hi, d1, d2) > tail {
        lo = hi;
        hi *= 2.0;
        if !hi.is_finite() {
            return f64::INFINITY;
        }
    }
    loop {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            return hi;
        }
        if f_tail(mid, d1, d2) > tail {
            lo = mid;
        } else {
            hi = mid;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_tails_match_their_closed_forms_and_tables() {
        // P(|Z| > 3), and P(|Z| > 1.959963984540054) = 0.05.
        assert!((normal_tail(3.0) - 0.002_699_796_063_260_186_6).abs() < 1e-15);
        assert!((normal_tail(1.959_963_984_540_054) - 0.05).abs() < 1e-14);
        // One degree of freedom is the Cauchy: P(|T| > t) = 1 - 2 atan(t) / π.
        let p = normal_tail(3.0);
        let cauchy = (core::f64::consts::FRAC_PI_2 * (1.0 - p)).tan();
        assert!((t_at_normal(3.0, 1) / cauchy - 1.0).abs() < 1e-10);
        // Two: P(|T| > t) = 1 - t / √(2 + t²).
        let q = 1.0 - p;
        let two = q * (2.0 / (1.0 - q * q)).sqrt();
        assert!((t_at_normal(3.0, 2) / two - 1.0).abs() < 1e-10);
        // The 5 % two-sided point of t with ten degrees of freedom.
        assert!((t_at_normal(1.959_963_984_540_054, 10) - 2.228_138_851_986_274).abs() < 1e-9);
        // Many: k (1 + (k² + 1) / 4ν), the expansion's first term.
        let many = t_at_normal(3.0, 1_000_000);
        assert!((many - 3.0 * (1.0 + 10.0 / 4e6)).abs() < 1e-9, "{many}");
        assert!(t_at_normal(3.0, 0).is_infinite());
        // The 5 % and 1 % two-sided points of the normal, and back.
        assert!((normal_at_tail(0.05) - 1.959_963_984_540_054).abs() < 1e-12);
        assert!((normal_at_tail(0.01) - 2.575_829_303_548_901).abs() < 1e-12);
        assert!((normal_at_tail(normal_tail(3.0)) - 3.0).abs() < 1e-12);
        assert_eq!(normal_at_tail(1.0), 0.0);
        assert!(normal_at_tail(0.0).is_infinite());
        assert!(t_at_normal(f64::INFINITY, 5).is_infinite());
    }

    #[test]
    fn the_f_quantile_matches_its_closed_forms_and_tables() {
        // F with one degree of freedom on top is T squared.
        for dof in [1, 4, 30] {
            let t = t_at_normal(3.0, dof);
            assert!((f_at_normal(3.0, 1, dof) / (t * t) - 1.0).abs() < 1e-10);
        }
        // With two on top, P(F > f) = (1 + 2 f / d2)^(-d2 / 2).
        let p = normal_tail(3.0);
        let two = 5.0 * (p.powf(-2.0 / 10.0) - 1.0);
        assert!((f_at_normal(3.0, 2, 10) / two - 1.0).abs() < 1e-10);
        // The 5 % points of F(6, 10) and F(3, 20).
        let five = normal_at_tail(0.05);
        assert!((f_at_normal(five, 6, 10) - 3.217_174_547_398_7).abs() < 1e-9);
        assert!((f_at_normal(five, 3, 20) - 3.098_391_212_074_7).abs() < 1e-9);
        assert!(f_at_normal(3.0, 0, 5).is_infinite());
        assert!(f_at_normal(3.0, 5, 0).is_infinite());
        assert!(f_at_normal(f64::INFINITY, 2, 5).is_infinite());
    }
}
