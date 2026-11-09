//! Radius laws: how a variable blend's radius runs along its chain.
//!
//! A law is read by the fraction of the chain's arc length from its start.
//! Between given radii it is a cubic Hermite interpolant whose tangents
//! follow Fritsch and Butland (a weighted harmonic mean of the neighbouring
//! secants, zero where they change sign), so each piece stays between the
//! radii at its two ends: the law never overshoots its data, and a law of
//! positive radii is positive everywhere.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d as _};

/// How a variable-radius blend's radius runs along its chain of edges, by
/// the fraction of the chain's arc length from the chain's start.
#[derive(Debug, Clone, PartialEq)]
pub enum RadiusLaw {
    /// From the chain's start to its end, linearly in arc length.
    Linear {
        /// The radius at the chain's start.
        start: f64,
        /// The radius at its end.
        end: f64,
    },
    /// Radii at fractions of the chain's arc length (`0..=1`, increasing),
    /// interpolated smoothly (C1) between them without overshooting them;
    /// held level before the first and after the last on an open chain,
    /// and wrapped round on a closed one.
    Points(Vec<(f64, f64)>),
}

/// A validated law, as cubic Hermite pieces over the fractions.
#[derive(Debug, Clone)]
pub(crate) struct Profile {
    /// The fractions the radii stand at, increasing.
    at: Vec<f64>,
    /// The radii.
    radius: Vec<f64>,
    /// The law's slope at each, per unit fraction.
    slope: Vec<f64>,
    /// Whether the chain closes: the law runs on past the last fraction
    /// round to the first.
    closed: bool,
}

impl Profile {
    /// The law on a chain that is `closed` or not, refused where it states
    /// no usable radius.
    pub(crate) fn new(law: &RadiusLaw, closed: bool, tol: Tolerances) -> OgeomResult<Self> {
        let usable = |r: f64| r.is_finite() && r > tol.confusion();
        let mut points: Vec<(f64, f64)> = match law {
            RadiusLaw::Linear { start, end } => {
                if closed && (start - end).abs() > tol.confusion() {
                    ogeom_bail!(
                        Construction,
                        "a linear law from {start} to {end} jumps where a closed \
                         chain meets its start; give its radii as points"
                    );
                }
                vec![(0.0, *start), (1.0, *end)]
            }
            RadiusLaw::Points(points) => points.clone(),
        };
        if points.is_empty() {
            ogeom_bail!(Construction, "a radius law of no points states no radius");
        }
        for &(s, r) in &points {
            if !s.is_finite() || !(0.0..=1.0).contains(&s) {
                ogeom_bail!(
                    Construction,
                    "a radius law's fraction {s} lies off its chain (0 to 1)"
                );
            }
            if !usable(r) {
                ogeom_bail!(Construction, "a fillet of radius {r} rounds nothing");
            }
        }
        if points.windows(2).any(|w| w[1].0 <= w[0].0) {
            ogeom_bail!(
                Construction,
                "a radius law's fractions must increase along the chain"
            );
        }
        // On a closed chain the end is the start: a radius given at both
        // is one radius, given twice.
        if closed && points.len() > 1 {
            let (first, last) = (points[0], points[points.len() - 1]);
            if first.0 <= tol.parametric() && last.0 >= 1.0 - tol.parametric() {
                if (first.1 - last.1).abs() > tol.confusion() {
                    ogeom_bail!(
                        Construction,
                        "a closed chain's law gives {} at its start and {} at \
                         its end, which are one point",
                        first.1,
                        last.1
                    );
                }
                points.pop();
            }
        }
        let at: Vec<f64> = points.iter().map(|p| p.0).collect();
        let radius: Vec<f64> = points.iter().map(|p| p.1).collect();
        let n = at.len();
        let mut slope = vec![0.0; n];
        if n > 1 {
            // The secant and span from point `i` to the next, wrapping
            // round a closed chain.
            let span = |i: usize| -> f64 {
                if i + 1 < n {
                    at[i + 1] - at[i]
                } else {
                    at[0] + 1.0 - at[i]
                }
            };
            let secant = |i: usize| -> f64 { (radius[(i + 1) % n] - radius[i]) / span(i) };
            for (i, m) in slope.iter_mut().enumerate() {
                let inner = closed || (i > 0 && i + 1 < n);
                if inner {
                    let before = if i == 0 { n - 1 } else { i - 1 };
                    let (h0, h1) = (span(before), span(i));
                    let (d0, d1) = (secant(before), secant(i));
                    *m = if d0 * d1 > 0.0 {
                        let (w0, w1) = (2.0f64.mul_add(h1, h0), 2.0f64.mul_add(h0, h1));
                        (w0 + w1) / (w0 / d0 + w1 / d1)
                    } else {
                        0.0
                    };
                } else {
                    // An open chain's first or last point: the secant to
                    // its neighbour where the point is the chain's end, and
                    // level where the law is held level beyond it.
                    let (end, d) = if i == 0 {
                        (at[0] <= tol.parametric(), secant(0))
                    } else {
                        (at[i] >= 1.0 - tol.parametric(), secant(i - 1))
                    };
                    *m = if end { d } else { 0.0 };
                }
            }
        }
        Ok(Self {
            at,
            radius,
            slope,
            closed,
        })
    }

    /// The radius and its slope per unit fraction at fraction `s`.
    pub(crate) fn at(&self, s: f64) -> (f64, f64) {
        let n = self.at.len();
        if n == 1 {
            return (self.radius[0], 0.0);
        }
        let hermite = |i: usize, j: usize, x0: f64, h: f64, s: f64| -> (f64, f64) {
            let t = (s - x0) / h;
            let (y0, y1) = (self.radius[i], self.radius[j]);
            let (m0, m1) = (self.slope[i] * h, self.slope[j] * h);
            let t2 = t * t;
            let t3 = t2 * t;
            let value = (2.0 * t3 - 3.0 * t2 + 1.0) * y0
                + (t3 - 2.0 * t2 + t) * m0
                + (3.0 * t2 - 2.0 * t3) * y1
                + (t3 - t2) * m1;
            let rate = (6.0 * t2 - 6.0 * t) * (y0 - y1)
                + (3.0 * t2 - 4.0 * t + 1.0) * m0
                + (3.0 * t2 - 2.0 * t) * m1;
            (value, rate / h)
        };
        let s = if self.closed { s.rem_euclid(1.0) } else { s };
        if self.closed {
            let (first, last) = (self.at[0], self.at[n - 1]);
            if s < first || s >= last {
                let x0 = last;
                let h = first + 1.0 - last;
                let s = if s < first { s + 1.0 } else { s };
                return hermite(n - 1, 0, x0, h, s);
            }
        } else {
            if s < self.at[0] {
                return (self.radius[0], 0.0);
            }
            if s > self.at[n - 1] {
                return (self.radius[n - 1], 0.0);
            }
        }
        let i = self
            .at
            .partition_point(|&x| x <= s)
            .saturating_sub(1)
            .min(n - 2);
        hermite(i, i + 1, self.at[i], self.at[i + 1] - self.at[i], s)
    }

    /// The law over one stretch of its chain, `from` to `to` (the
    /// stretch's start and end by fraction, `to` below `from` where the
    /// stretch runs against the chain), read by the fraction of the
    /// stretch: an open law of the same radii. Each Hermite piece is
    /// restricted with its value and slope at the stretch's ends, so the
    /// window is the law itself there, not a fresh interpolation.
    pub(crate) fn window(&self, from: f64, to: f64) -> Self {
        let span = to - from;
        let (lo, hi) = (from.min(to), from.max(to));
        let mut nodes: Vec<(f64, f64, f64)> = Vec::new();
        let mut push = |s: f64| {
            let (r, m) = self.at(s);
            nodes.push(((s - from) / span, r, m * span));
        };
        push(from);
        push(to);
        for &s in &self.at {
            if s > lo && s < hi {
                push(s);
            }
        }
        nodes.sort_by(|a, b| a.0.total_cmp(&b.0));
        Self {
            at: nodes.iter().map(|n| n.0).collect(),
            radius: nodes.iter().map(|n| n.1).collect(),
            slope: nodes.iter().map(|n| n.2).collect(),
            closed: false,
        }
    }

    /// The largest radius the law reaches: at a given point, since no
    /// piece overshoots its ends.
    pub(crate) fn largest(&self) -> f64 {
        self.radius
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// The smallest.
    pub(crate) fn least(&self) -> f64 {
        self.radius.iter().copied().fold(f64::INFINITY, f64::min)
    }

    /// Whether the law was laid on a chain that closes.
    pub(crate) const fn closed(&self) -> bool {
        self.closed
    }

    /// Whether the law holds one radius all along.
    pub(crate) fn level(&self, tol: Tolerances) -> bool {
        self.largest() - self.least() <= tol.confusion()
    }

    /// The open law as a cubic B-spline over the fraction: the knots,
    /// doubled inside (C1), and the radii's control values. Each piece is
    /// its Hermite cubic's own Bezier form, and the C1 join is the point
    /// knot removal drops.
    pub(crate) fn spline(&self, tol: Tolerances) -> (Vec<f64>, Vec<f64>) {
        // The pieces over the whole chain: the level runs before the first
        // point and after the last, and the Hermite pieces between.
        let mut pieces: Vec<(f64, f64, [f64; 4])> = Vec::new();
        let n = self.at.len();
        let level = |x0: f64, x1: f64, r: f64| (x0, x1, [r, r, r, r]);
        if self.at[0] > tol.parametric() {
            pieces.push(level(0.0, self.at[0], self.radius[0]));
        }
        for i in 0..n.saturating_sub(1) {
            let (x0, x1) = (self.at[i], self.at[i + 1]);
            let h = x1 - x0;
            let (y0, y1) = (self.radius[i], self.radius[i + 1]);
            pieces.push((
                x0,
                x1,
                [
                    y0,
                    self.slope[i].mul_add(h / 3.0, y0),
                    (-self.slope[i + 1]).mul_add(h / 3.0, y1),
                    y1,
                ],
            ));
        }
        if self.at[n - 1] < 1.0 - tol.parametric() {
            pieces.push(level(self.at[n - 1], 1.0, self.radius[n - 1]));
        }
        if pieces.is_empty() {
            pieces.push(level(0.0, 1.0, self.radius[0]));
        }
        // The chain's own ends stand at 0 and 1 exactly.
        let last = pieces.len() - 1;
        pieces[0].0 = 0.0;
        pieces[last].1 = 1.0;
        let mut knots = vec![0.0; 4];
        let mut values = vec![pieces[0].2[0]];
        for (k, (_, x1, bezier)) in pieces.iter().enumerate() {
            values.push(bezier[1]);
            values.push(bezier[2]);
            if k == last {
                values.push(bezier[3]);
                knots.extend([1.0; 4]);
            } else {
                knots.extend([*x1, *x1]);
            }
        }
        (knots, values)
    }
}

/// How many pieces a curve's arc length is tabulated over.
const RULER_PIECES: u32 = 256;

/// Gauss-Legendre nodes and weights on `[-1, 1]`, five points: exact for
/// polynomials to degree nine, so a piece of a smooth curve's speed is
/// integrated to rounding.
const GAUSS: [(f64, f64); 5] = [
    (0.0, 0.568_888_888_888_888_9),
    (-0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
    (0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
    (-0.906_179_845_938_664, 0.236_926_885_056_189_1),
    (0.906_179_845_938_664, 0.236_926_885_056_189_1),
];

/// A curve's arc length from the start of a range, tabulated, read back
/// as the fraction of the whole run.
#[derive(Debug, Clone)]
pub(crate) struct Ruler {
    curve: Curve,
    range: (f64, f64),
    /// The length up to each piece's start, and the whole at the end.
    cumulative: Vec<f64>,
    closed: bool,
}

impl Ruler {
    pub(crate) fn new(
        curve: Curve,
        range: (f64, f64),
        closed: bool,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        let mut cumulative: Vec<f64> = Vec::with_capacity(RULER_PIECES as usize + 1);
        cumulative.push(0.0);
        let step = (range.1 - range.0) / f64::from(RULER_PIECES);
        let mut whole = 0.0;
        for k in 0..RULER_PIECES {
            let a = f64::from(k).mul_add(step, range.0);
            whole += speed_integral(&curve, a, a + step, tol)?;
            cumulative.push(whole);
        }
        if whole <= tol.confusion() {
            ogeom_bail!(Construction, "a chain of no length carries no radius law");
        }
        Ok(Self {
            curve,
            range,
            cumulative,
            closed,
        })
    }

    /// The whole run's length.
    pub(crate) fn length(&self) -> f64 {
        self.cumulative[self.cumulative.len() - 1]
    }

    /// The fraction of the run's length up to parameter `t`, and its rate
    /// per unit parameter.
    pub(crate) fn fraction(&self, t: f64, tol: Tolerances) -> OgeomResult<(f64, f64)> {
        let (lo, hi) = self.range;
        let span = hi - lo;
        let (t, turns) = if self.closed {
            let turns = ((t - lo) / span).floor();
            (turns.mul_add(-span, t), turns)
        } else {
            (t.clamp(lo, hi), 0.0)
        };
        let step = span / f64::from(RULER_PIECES);
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "a clamped piece index"
        )]
        let k = (((t - lo) / step).floor().max(0.0) as u32).min(RULER_PIECES - 1);
        let a = f64::from(k).mul_add(step, lo);
        let length = self.cumulative[k as usize] + speed_integral(&self.curve, a, t, tol)?;
        let speed = self.curve.d1_at(t, tol)?.magnitude();
        let whole = self.length();
        Ok((length / whole + turns, speed / whole))
    }
}

/// The curve's length from `a` to `b`, by five-point Gauss-Legendre.
fn speed_integral(curve: &Curve, a: f64, b: f64, tol: Tolerances) -> OgeomResult<f64> {
    let (mid, half) = (f64::midpoint(a, b), (b - a) / 2.0);
    let mut sum = 0.0;
    for (x, w) in GAUSS {
        sum += w * curve.d1_at(half.mul_add(x, mid), tol)?.magnitude();
    }
    Ok(sum * half)
}

/// A law laid along the curve a march walks: the radius at the curve's
/// parameter, and its rate per unit parameter, for the marcher's
/// Jacobian.
#[derive(Debug, Clone)]
pub(crate) struct RadiusAlong {
    profile: Profile,
    ruler: Ruler,
    /// The fraction of the curve's run at the chain's start.
    origin: f64,
    /// Whether the chain runs against the curve's parameter.
    reversed: bool,
}

impl RadiusAlong {
    pub(crate) const fn new(profile: Profile, ruler: Ruler, origin: f64, reversed: bool) -> Self {
        Self {
            profile,
            ruler,
            origin,
            reversed,
        }
    }

    /// The radius at the curve's parameter `t`, and its rate per unit `t`.
    pub(crate) fn at(&self, t: f64, tol: Tolerances) -> OgeomResult<(f64, f64)> {
        let (s, ds) = self.ruler.fraction(t, tol)?;
        let s = s - self.origin;
        let (s, ds) = if self.reversed { (-s, -ds) } else { (s, ds) };
        let (r, dr) = self.profile.at(s);
        Ok((r, dr * ds))
    }

    /// The law's own radii, smallest and largest.
    pub(crate) fn least(&self) -> f64 {
        self.profile.least()
    }

    pub(crate) fn largest(&self) -> f64 {
        self.profile.largest()
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::{Profile, RadiusLaw};
    use ogeom_core::Tolerances;

    /// A window is the law itself over its stretch, read by the stretch's
    /// own fraction: forward, against the chain, and up to a closed chain's
    /// end, where the law wraps.
    #[test]
    fn a_window_is_the_law_over_its_stretch() {
        let tol = Tolerances::millimetres();
        let law = RadiusLaw::Points(vec![(0.1, 0.5), (0.4, 1.5), (0.7, 0.8)]);
        let closed = Profile::new(&law, true, tol).unwrap();
        for (from, to) in [(0.05, 0.55), (0.9, 0.3), (0.65, 1.0), (0.2, 0.25)] {
            let window = closed.window(from, to);
            for i in 0..=100 {
                let u = f64::from(i) / 100.0;
                let (r, m) = window.at(u);
                let (want, rate) = closed.at(from + (to - from) * u);
                assert!(
                    (r - want).abs() < 1e-12,
                    "{from}..{to} at {u}: {r} against {want}"
                );
                assert!(
                    (m - rate * (to - from)).abs() < 1e-9,
                    "{from}..{to} at {u}: slope {m} against {}",
                    rate * (to - from)
                );
            }
        }
    }
}
