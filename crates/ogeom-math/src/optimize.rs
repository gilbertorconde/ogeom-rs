//! Minima of functions of several variables over a box: local, and global.
//!
//! [`minimize_local`] walks downhill from a start by the Nelder-Mead simplex,
//! which asks for values only. [`global_minimum`] finds the least value over
//! the whole box by branch and bound: each sub-box's floor is its centre's
//! value less a Lipschitz constant times its half-diagonal, the box with the
//! lowest floor is split, and the search stops when no floor lies more than
//! the tolerance below the best value found. The constant is estimated from
//! the function's own slopes, so the certificate holds as far as that
//! estimate does. [`swarm_minimum`] is the particle swarm, for a function
//! too rough for a Lipschitz bound to say much.

use std::collections::BinaryHeap;

use ogeom_core::{OgeomResult, ogeom_bail};

/// A minimum found over a box.
#[derive(Debug, Clone, PartialEq)]
pub struct Minimum {
    /// Where.
    pub point: Vec<f64>,
    /// The function's value there.
    pub value: f64,
    /// Whether the search proved no point of the box lies lower by more
    /// than its tolerance (under the slope bound it estimated), rather than
    /// stopping at its evaluation budget.
    pub certified: bool,
    /// Function evaluations spent.
    pub evaluations: usize,
}

fn check_box(lower: &[f64], upper: &[f64]) -> OgeomResult<()> {
    if lower.is_empty() || lower.len() != upper.len() {
        ogeom_bail!(
            Construction,
            "a box needs matching, non-empty bounds; got {} and {}",
            lower.len(),
            upper.len()
        );
    }
    if lower
        .iter()
        .zip(upper)
        .any(|(a, b)| !a.is_finite() || !b.is_finite() || b <= a)
    {
        ogeom_bail!(Construction, "a box's bounds must be finite and increasing");
    }
    Ok(())
}

/// One evaluation, counted, a NaN read as no minimum.
fn call<F: FnMut(&[f64]) -> f64>(f: &mut F, x: &[f64], evaluations: &mut usize) -> f64 {
    *evaluations += 1;
    let v = f(x);
    if v.is_nan() { f64::INFINITY } else { v }
}

fn clamp_into(x: &mut [f64], lower: &[f64], upper: &[f64]) {
    for ((v, a), b) in x.iter_mut().zip(lower).zip(upper) {
        *v = v.clamp(*a, *b);
    }
}

/// A local minimum of `f` near `start` within the box, by the Nelder-Mead
/// simplex from a first simplex of `step` along each axis. Points are held
/// to the box.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// box is malformed or `start` does not match it.
pub fn minimize_local<F: FnMut(&[f64]) -> f64>(
    mut f: F,
    start: &[f64],
    lower: &[f64],
    upper: &[f64],
    step: f64,
    tolerance: f64,
    max_evaluations: usize,
) -> OgeomResult<Minimum> {
    check_box(lower, upper)?;
    if start.len() != lower.len() {
        ogeom_bail!(
            Construction,
            "the start has {} coordinates, the box {}",
            start.len(),
            lower.len()
        );
    }
    let n = start.len();
    let mut evaluations = 0usize;
    let mut eval = |x: &mut Vec<f64>, evaluations: &mut usize| -> f64 {
        clamp_into(x, lower, upper);
        call(&mut f, x, evaluations)
    };
    let mut simplex: Vec<(Vec<f64>, f64)> = Vec::with_capacity(n + 1);
    let mut first = start.to_vec();
    let v = eval(&mut first, &mut evaluations);
    simplex.push((first, v));
    for i in 0..n {
        let mut x = start.to_vec();
        x[i] += if x[i] + step <= upper[i] { step } else { -step };
        let v = eval(&mut x, &mut evaluations);
        simplex.push((x, v));
    }
    while evaluations < max_evaluations {
        simplex.sort_by(|a, b| a.1.total_cmp(&b.1));
        let (best, worst) = (simplex[0].1, simplex[n].1);
        let size = simplex[1..]
            .iter()
            .map(|(x, _)| {
                x.iter()
                    .zip(&simplex[0].0)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0_f64, f64::max)
            })
            .fold(0.0_f64, f64::max);
        if (worst - best).abs() <= tolerance && size <= tolerance.sqrt() * 1e-3 + 1e-14 {
            break;
        }
        if size <= 1e-15 {
            break;
        }
        #[allow(clippy::cast_precision_loss)]
        let centroid: Vec<f64> = (0..n)
            .map(|k| simplex[..n].iter().map(|(x, _)| x[k]).sum::<f64>() / n as f64)
            .collect();
        let toward = |t: f64| -> Vec<f64> {
            centroid
                .iter()
                .zip(&simplex[n].0)
                .map(|(c, w)| c + t * (c - w))
                .collect()
        };
        let mut reflected = toward(1.0);
        let fr = eval(&mut reflected, &mut evaluations);
        if fr < simplex[0].1 {
            let mut expanded = toward(2.0);
            let fe = eval(&mut expanded, &mut evaluations);
            simplex[n] = if fe < fr {
                (expanded, fe)
            } else {
                (reflected, fr)
            };
        } else if fr < simplex[n - 1].1 {
            simplex[n] = (reflected, fr);
        } else {
            let mut contracted = if fr < simplex[n].1 {
                toward(0.5)
            } else {
                toward(-0.5)
            };
            let fc = eval(&mut contracted, &mut evaluations);
            if fc < simplex[n].1.min(fr) {
                simplex[n] = (contracted, fc);
            } else {
                // Shrink toward the best.
                let anchor = simplex[0].0.clone();
                for entry in simplex.iter_mut().skip(1) {
                    let mut x: Vec<f64> = entry
                        .0
                        .iter()
                        .zip(&anchor)
                        .map(|(p, a)| a + 0.5 * (p - a))
                        .collect();
                    let v = eval(&mut x, &mut evaluations);
                    *entry = (x, v);
                }
            }
        }
    }
    simplex.sort_by(|a, b| a.1.total_cmp(&b.1));
    let (point, value) = simplex.swap_remove(0);
    Ok(Minimum {
        point,
        value,
        certified: false,
        evaluations,
    })
}

/// One box of the search, ordered by its floor, lowest first.
struct Cell {
    floor: f64,
    lower: Vec<f64>,
    upper: Vec<f64>,
}

impl PartialEq for Cell {
    fn eq(&self, other: &Self) -> bool {
        self.floor.total_cmp(&other.floor).is_eq()
    }
}
impl Eq for Cell {}
impl PartialOrd for Cell {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Cell {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        other.floor.total_cmp(&self.floor)
    }
}

/// The bounds a box's floor is read from.
struct Bounds {
    lipschitz: f64,
    curvature: f64,
}

impl Bounds {
    /// A box's centre value and floor: the tighter of the slope bound from
    /// its centre, and the centre's gradient with the curvature bound,
    /// which closes on a smooth minimum as the box shrinks.
    fn probe<F: FnMut(&[f64]) -> f64>(
        &self,
        f: &mut F,
        c: &[f64],
        r: f64,
        evaluations: &mut usize,
    ) -> (f64, f64) {
        let fc = call(f, c, evaluations);
        let mut gradient = 0.0_f64;
        let step = (r * 1e-3).max(1e-9);
        for k in 0..c.len() {
            let mut a = c.to_vec();
            let mut b = c.to_vec();
            a[k] += step;
            b[k] -= step;
            let d = (call(f, &a, evaluations) - call(f, &b, evaluations)) / (2.0 * step);
            gradient += d * d;
        }
        let first = fc - self.lipschitz * r;
        let second = fc - gradient.sqrt() * r - 0.5 * self.curvature * r * r;
        (fc, first.max(second))
    }
}

/// The least value of `f` over the box, to within `tolerance`, by
/// Lipschitz branch and bound with local polishing.
///
/// A box's floor is the tighter of two: its centre's value less a slope
/// bound times its half-diagonal, and less the centre's gradient times the
/// half-diagonal and half a curvature bound times its square, which closes
/// on a smooth minimum as the box shrinks. Both bounds are estimated from
/// finite differences over a sampling of the box and doubled; a function
/// steeper or more sharply bent somewhere than any sample shows can hide a
/// minimum from the certificate. `max_evaluations`
/// caps the work; a search stopped by it reports `certified: false`.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// box is malformed or the tolerance is not positive.
pub fn global_minimum<F: FnMut(&[f64]) -> f64>(
    mut f: F,
    lower: &[f64],
    upper: &[f64],
    tolerance: f64,
    max_evaluations: usize,
) -> OgeomResult<Minimum> {
    check_box(lower, upper)?;
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        ogeom_bail!(Construction, "a tolerance of {tolerance} is not positive");
    }
    let n = lower.len();
    let mut evaluations = 0usize;

    // The slope and curvature bounds: finite differences at a
    // low-discrepancy sampling, along the axes and the diagonals.
    let widths: Vec<f64> = lower.iter().zip(upper).map(|(a, b)| b - a).collect();
    let mut best_point: Vec<f64> = lower
        .iter()
        .zip(upper)
        .map(|(a, b)| 0.5 * (a + b))
        .collect();
    let mut best = call(&mut f, &best_point, &mut evaluations);
    let (mut slope, mut bend) = (0.0_f64, 0.0_f64);
    let directions: Vec<Vec<f64>> = (0..n)
        .map(|k| (0..n).map(|j| if j == k { 1.0 } else { 0.0 }).collect())
        .chain((0..n).map(|k| {
            #[allow(clippy::cast_precision_loss)]
            let norm = (n as f64).sqrt();
            (0..n)
                .map(|j| if j < k { -1.0 / norm } else { 1.0 / norm })
                .collect()
        }))
        .collect();
    let h = widths.iter().copied().fold(f64::INFINITY, f64::min) * 1e-3;
    for i in 1..=32 * n {
        let x: Vec<f64> = (0..n)
            .map(|k| {
                let t = halton(i, PRIMES[k % PRIMES.len()]);
                lower[k] + h + (widths[k] - 2.0 * h) * t
            })
            .collect();
        let fx = call(&mut f, &x, &mut evaluations);
        if fx < best {
            (best, best_point) = (fx, x.clone());
        }
        for d in &directions {
            let ahead: Vec<f64> = x.iter().zip(d).map(|(v, e)| v + h * e).collect();
            let behind: Vec<f64> = x.iter().zip(d).map(|(v, e)| v - h * e).collect();
            let (fa, fb) = (
                call(&mut f, &ahead, &mut evaluations),
                call(&mut f, &behind, &mut evaluations),
            );
            if fx.is_finite() && fa.is_finite() && fb.is_finite() {
                slope = slope.max((fa - fb).abs() / (2.0 * h));
                bend = bend.max((fa - 2.0 * fx + fb).abs() / (h * h));
            }
        }
    }
    // Doubled for safety, and the slope bound widened by the dimension:
    // the samples read it along a few directions only.
    #[allow(clippy::cast_precision_loss)]
    let lipschitz = (2.0 * slope * (n as f64).sqrt()).max(tolerance);
    #[allow(clippy::cast_precision_loss)]
    let curvature = (2.0 * bend * n as f64).max(tolerance);

    let half_diagonal = |lo: &[f64], hi: &[f64]| -> f64 {
        lo.iter()
            .zip(hi)
            .map(|(a, b)| (0.5 * (b - a)).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    let bounds = Bounds {
        lipschitz,
        curvature,
    };
    let mut heap = BinaryHeap::new();
    let centre: Vec<f64> = lower
        .iter()
        .zip(upper)
        .map(|(a, b)| 0.5 * (a + b))
        .collect();
    let (_, floor) = bounds.probe(
        &mut f,
        &centre,
        half_diagonal(lower, upper),
        &mut evaluations,
    );
    heap.push(Cell {
        floor,
        lower: lower.to_vec(),
        upper: upper.to_vec(),
    });
    let mut certified = false;
    let step = widths.iter().copied().fold(f64::INFINITY, f64::min) * 1e-2;
    while let Some(cell) = heap.pop() {
        if cell.floor >= best - tolerance {
            certified = true;
            break;
        }
        if evaluations >= max_evaluations {
            heap.push(cell);
            break;
        }
        // Split across the longest side.
        let k = (0..n)
            .max_by(|a, b| {
                (cell.upper[*a] - cell.lower[*a]).total_cmp(&(cell.upper[*b] - cell.lower[*b]))
            })
            .unwrap_or(0);
        let mid = 0.5 * (cell.lower[k] + cell.upper[k]);
        for half in 0..2 {
            let (mut lo, mut hi) = (cell.lower.clone(), cell.upper.clone());
            if half == 0 {
                hi[k] = mid;
            } else {
                lo[k] = mid;
            }
            let c: Vec<f64> = lo.iter().zip(&hi).map(|(a, b)| 0.5 * (a + b)).collect();
            let (fc, floor) = bounds.probe(&mut f, &c, half_diagonal(&lo, &hi), &mut evaluations);
            if fc < best - tolerance {
                // A new basin: polished before it is trusted.
                let polished =
                    minimize_local(&mut f, &c, lower, upper, step, tolerance * 1e-3, 400)?;
                evaluations += polished.evaluations;
                if polished.value < fc {
                    (best, best_point) = (polished.value, polished.point);
                } else {
                    (best, best_point) = (fc, c.clone());
                }
            } else if fc < best {
                (best, best_point) = (fc, c.clone());
            }
            heap.push(Cell {
                floor,
                lower: lo,
                upper: hi,
            });
        }
    }
    if heap.is_empty() {
        certified = true;
    }
    Ok(Minimum {
        point: best_point,
        value: best,
        certified,
        evaluations,
    })
}

/// A minimum of `f` over the box by particle swarm: `particles` points
/// moving under their own best and the swarm's best for `iterations`
/// rounds, from a deterministic scattering, then polished locally. No
/// certificate: the swarm reports the best it saw.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// box is malformed or there are no particles.
pub fn swarm_minimum<F: FnMut(&[f64]) -> f64>(
    mut f: F,
    lower: &[f64],
    upper: &[f64],
    particles: usize,
    iterations: usize,
) -> OgeomResult<Minimum> {
    check_box(lower, upper)?;
    if particles == 0 {
        ogeom_bail!(Construction, "a swarm needs particles");
    }
    let n = lower.len();
    let mut evaluations = 0usize;
    let widths: Vec<f64> = lower.iter().zip(upper).map(|(a, b)| b - a).collect();
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let mut random = move || -> f64 {
        // xorshift64*, deterministic so a result reproduces.
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        #[allow(clippy::cast_precision_loss)]
        let r = (state.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64 / (1u64 << 53) as f64;
        r
    };
    let mut position: Vec<Vec<f64>> = (1..=particles)
        .map(|i| {
            (0..n)
                .map(|k| lower[k] + widths[k] * halton(i, PRIMES[k % PRIMES.len()]))
                .collect()
        })
        .collect();
    let mut velocity: Vec<Vec<f64>> = vec![vec![0.0; n]; particles];
    let mut own_best: Vec<(Vec<f64>, f64)> = position
        .iter()
        .map(|x| (x.clone(), call(&mut f, x, &mut evaluations)))
        .collect();
    let mut swarm_best = own_best
        .iter()
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .cloned()
        .unwrap_or_else(|| (position[0].clone(), f64::INFINITY));
    let (inertia, own, social) = (0.72, 1.49, 1.49);
    for _ in 0..iterations {
        for p in 0..particles {
            for k in 0..n {
                let (r1, r2) = (random(), random());
                velocity[p][k] = inertia * velocity[p][k]
                    + own * r1 * (own_best[p].0[k] - position[p][k])
                    + social * r2 * (swarm_best.0[k] - position[p][k]);
                velocity[p][k] = velocity[p][k].clamp(-widths[k], widths[k]);
                position[p][k] = (position[p][k] + velocity[p][k]).clamp(lower[k], upper[k]);
            }
            let v = call(&mut f, &position[p], &mut evaluations);
            if v < own_best[p].1 {
                own_best[p] = (position[p].clone(), v);
                if v < swarm_best.1 {
                    swarm_best = (position[p].clone(), v);
                }
            }
        }
    }
    let step = widths.iter().copied().fold(f64::INFINITY, f64::min) * 1e-3;
    let polished = minimize_local(&mut f, &swarm_best.0, lower, upper, step, 1e-14, 400)?;
    let (point, best) = if polished.value < swarm_best.1 {
        (polished.point, polished.value)
    } else {
        swarm_best
    };
    Ok(Minimum {
        point,
        value: best,
        certified: false,
        evaluations: evaluations + polished.evaluations,
    })
}

const PRIMES: [usize; 8] = [2, 3, 5, 7, 11, 13, 17, 19];

/// The `i`th element of the van der Corput sequence in `base`.
fn halton(mut i: usize, base: usize) -> f64 {
    let mut result = 0.0;
    #[allow(clippy::cast_precision_loss)]
    let mut fraction = 1.0 / base as f64;
    while i > 0 {
        #[allow(clippy::cast_precision_loss)]
        {
            result += fraction * (i % base) as f64;
            fraction /= base as f64;
        }
        i /= base;
    }
    result
}
