//! Chebyshev fits by sequential linear programming in a trust region.
//!
//! Each step linearises every point's distance about the current shape
//! and solves the linear program for the step, within a box of the trust
//! radius, that makes the linearised zone narrowest; the step is taken if
//! the true objective falls by a fair share of what the model predicted,
//! and the radius grows or shrinks by how well the model predicted it
//! (Nocedal and Wright, *Numerical Optimization*, algorithm 4.1, applied to
//! a polyhedral function of smooth distances as in Fletcher's model
//! algorithm for composite nondifferentiable optimisation). From any start
//! the steps converge to a point where no first-order change of the
//! parameters narrows the zone.
//!
//! The steps run until the model predicts no fall beyond rounding. The
//! fit is then optimal if the linear model, solved over a box as wide as
//! the points themselves, predicts no narrowing beyond the confusion
//! distance: the Chebyshev conditions, that the contacts' gradients
//! balance so no change of the parameters moves every contact inward at
//! once, hold within that distance.
//!
//! Each linear program has one row per point and side; it is solved over
//! a working set of rows, seeded with the extreme points and grown by the
//! rows the working optimum breaks until it breaks none, which is the
//! optimum over every row.

use ogeom_math::{Point, Vector};

use super::lp::{self, ROUNDING, Solved};

/// The most parameters a shape has: an axis's four.
pub(super) const MOST_PARAMETERS: usize = 4;

/// The most steps a fit takes. Near a solution with as many contacts as
/// the zone has freedoms the steps converge as Newton's method does, in a
/// handful; a fit still moving after this many is reported not optimal
/// rather than left to run on.
const MOST_STEPS: usize = 200;

/// A step is taken when the objective falls by at least this share of the
/// predicted fall (any value in `[0, 1/4)` keeps the convergence proof).
const ACCEPT: f64 = 0.1;
/// Below this share the model was poor and the radius shrinks to a
/// quarter of the step; above [`GROW`], it doubles.
const SHRINK: f64 = 0.25;
/// See [`SHRINK`].
const GROW: f64 = 0.75;

/// What the zone makes least.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Goal {
    /// The greatest distance less the least: a minimum zone.
    Zone,
    /// The greatest distance: a minimum circumscribed shape.
    Outer,
    /// The least distance, made greatest: a maximum inscribed shape.
    Inner,
}

/// A shape fitted by its distances to points.
pub(super) trait Shape: Sized {
    /// How many parameters it has, at most [`MOST_PARAMETERS`].
    fn parameters(&self) -> usize;
    /// The point's distance from the shape, and its gradient in the
    /// parameters, zero where the distance is.
    fn measure(&self, p: Point, gradient: &mut [f64; MOST_PARAMETERS]) -> f64;
    /// The shape moved by a step in its parameters.
    fn moved(&self, step: &[f64]) -> Option<Self>;
}

/// A fit's end.
#[derive(Debug, Clone)]
pub(super) struct Settled<S> {
    pub shape: S,
    /// Each point's distance from the shape.
    pub distances: Vec<f64>,
    /// Whether the Chebyshev conditions hold.
    pub optimal: bool,
}

struct Measured {
    d: Vec<f64>,
    g: Vec<[f64; MOST_PARAMETERS]>,
    least: f64,
    most: f64,
    objective: f64,
}

fn measure<S: Shape>(shape: &S, points: &[Point], goal: Goal) -> Measured {
    let mut d = Vec::with_capacity(points.len());
    let mut g = Vec::with_capacity(points.len());
    let (mut least, mut most) = (f64::INFINITY, f64::NEG_INFINITY);
    for p in points {
        let mut gradient = [0.0; MOST_PARAMETERS];
        let distance = shape.measure(*p, &mut gradient);
        least = least.min(distance);
        most = most.max(distance);
        d.push(distance);
        g.push(gradient);
    }
    let objective = match goal {
        Goal::Zone => most - least,
        Goal::Outer => most,
        Goal::Inner => -least,
    };
    Measured {
        d,
        g,
        least,
        most,
        objective,
    }
}

/// Settle `start` to the Chebyshev optimum of `goal` nearest it. `reach` is
/// the points' size, the widest step the optimality test allows; `slack`
/// the narrowing below which the zone counts as least.
pub(super) fn settle<S: Shape>(
    start: S,
    points: &[Point],
    goal: Goal,
    reach: f64,
    slack: f64,
) -> Settled<S> {
    let k = start.parameters();
    let mut shape = start;
    let mut at = measure(&shape, points, goal);
    let mut radius = (at.most - at.least).clamp(slack, reach.max(slack));
    for _ in 0..MOST_STEPS {
        let Some((step, value)) = model(&at, k, goal, radius) else {
            break;
        };
        // A fall within rounding of the distances is no fall.
        let floor = ROUNDING * (at.most.abs() + reach);
        let predicted = at.objective - value;
        if predicted <= floor {
            if radius >= reach {
                break;
            }
            match model(&at, k, goal, reach) {
                Some((_, far)) if at.objective - far > floor => {
                    radius = (radius * 2.0).min(reach);
                    continue;
                }
                _ => break,
            }
        }
        let length = step.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let Some(trial) = shape.moved(&step) else {
            radius = length / 4.0;
            continue;
        };
        let next = measure(&trial, points, goal);
        let ratio = (at.objective - next.objective) / predicted;
        if ratio >= ACCEPT {
            shape = trial;
            at = next;
        }
        if ratio < SHRINK {
            radius = length / 4.0;
        } else if ratio > GROW {
            radius = (radius * 2.0).min(reach);
        }
        if radius <= reach * f64::EPSILON {
            break;
        }
    }
    let optimal = model(&at, k, goal, reach).is_some_and(|(_, far)| at.objective - far <= slack)
        // A point at distance zero sits on the distance's cone point,
        // where no gradient speaks for it and the test says nothing.
        && (goal == Goal::Outer || at.least > 0.0);
    Settled {
        shape,
        distances: at.d,
        optimal,
    }
}

/// Which rows of a point are in the working set.
const UPPER: u8 = 1;
const LOWER: u8 = 2;

/// The step within `radius` that makes the linearised objective least, and
/// that least value; `None` where the program could not be solved.
fn model(at: &Measured, k: usize, goal: Goal, radius: f64) -> Option<(Vec<f64>, f64)> {
    // Variables: the step, then the upper bound, then the lower.
    let (upper, lower) = match goal {
        Goal::Zone => (Some(k), Some(k + 1)),
        Goal::Outer => (Some(k), None),
        Goal::Inner => (None, Some(k)),
    };
    let n = k + usize::from(upper.is_some()) + usize::from(lower.is_some());
    let mut c = vec![0.0; n];
    if let Some(u) = upper {
        c[u] = 1.0;
    }
    if let Some(l) = lower {
        c[l] = -1.0;
    }

    let mut rows: Vec<(Vec<f64>, f64)> = Vec::new();
    for j in 0..k {
        for sign in [1.0, -1.0] {
            let mut a = vec![0.0; n];
            a[j] = sign;
            rows.push((a, radius));
        }
    }
    let mut held = vec![0_u8; at.d.len()];
    let row = |i: usize, side: u8| -> (Vec<f64>, f64) {
        let mut a = vec![0.0; n];
        let g = &at.g[i];
        if side == UPPER {
            // d + g·s ≤ t_upper
            a[..k].copy_from_slice(&g[..k]);
            if let Some(u) = upper {
                a[u] = -1.0;
            }
            (a, -at.d[i])
        } else {
            // d + g·s ≥ t_lower
            for j in 0..k {
                a[j] = -g[j];
            }
            if let Some(l) = lower {
                a[l] = 1.0;
            }
            (a, at.d[i])
        }
    };
    let seeds = k + 2;
    if upper.is_some() {
        for i in largest(at.d.len(), seeds, |i| Some(at.d[i])) {
            held[i] |= UPPER;
            rows.push(row(i, UPPER));
        }
    }
    if lower.is_some() {
        for i in largest(at.d.len(), seeds, |i| Some(-at.d[i])) {
            held[i] |= LOWER;
            rows.push(row(i, LOWER));
        }
    }

    loop {
        let borrowed: Vec<(&[f64], f64)> = rows.iter().map(|(a, b)| (a.as_slice(), *b)).collect();
        let Solved::Optimal(z) = lp::minimise(&c, &borrowed) else {
            return None;
        };
        let s = &z[..k];
        let moved = |i: usize| {
            let g = &at.g[i];
            let (mut sum, mut size) = (at.d[i], at.d[i].abs());
            for j in 0..k {
                sum += g[j] * s[j];
                size += (g[j] * s[j]).abs();
            }
            (sum, size)
        };
        let mut added = false;
        if let Some(u) = upper {
            let t = z[u];
            for i in largest(at.d.len(), seeds, |i| {
                if held[i] & UPPER != 0 {
                    return None;
                }
                let (value, size) = moved(i);
                let excess = value - t;
                (excess > ROUNDING * (size + t.abs())).then_some(excess)
            }) {
                held[i] |= UPPER;
                rows.push(row(i, UPPER));
                added = true;
            }
        }
        if let Some(l) = lower {
            let t = z[l];
            for i in largest(at.d.len(), seeds, |i| {
                if held[i] & LOWER != 0 {
                    return None;
                }
                let (value, size) = moved(i);
                let excess = t - value;
                (excess > ROUNDING * (size + t.abs())).then_some(excess)
            }) {
                held[i] |= LOWER;
                rows.push(row(i, LOWER));
                added = true;
            }
        }
        if !added {
            let value = c.iter().zip(&z).map(|(a, b)| a * b).sum();
            return Some((s.to_vec(), value));
        }
    }
}

/// The indices of the `count` greatest keys, greatest first, ties to the
/// lower index; an index whose key is `None` is passed over.
pub(super) fn largest(len: usize, count: usize, key: impl Fn(usize) -> Option<f64>) -> Vec<usize> {
    let mut best: Vec<(f64, usize)> = Vec::with_capacity(count + 1);
    for i in 0..len {
        let Some(v) = key(i) else {
            continue;
        };
        if best.len() == count && best.last().is_some_and(|&(w, _)| v <= w) {
            continue;
        }
        let at = best.partition_point(|&(w, _)| w >= v);
        best.insert(at, (v, i));
        best.truncate(count);
    }
    best.into_iter().map(|(_, i)| i).collect()
}

/// A sphere's centre, its distances those to the centre.
#[derive(Debug, Clone, Copy)]
pub(super) struct Centre(pub Point);

impl Shape for Centre {
    fn parameters(&self) -> usize {
        3
    }

    fn measure(&self, p: Point, gradient: &mut [f64; MOST_PARAMETERS]) -> f64 {
        let v = p - self.0;
        let d = v.magnitude();
        if d > 0.0 {
            gradient[0] = -v.x / d;
            gradient[1] = -v.y / d;
            gradient[2] = -v.z / d;
        }
        d
    }

    fn moved(&self, step: &[f64]) -> Option<Self> {
        Some(Self(self.0 + Vector::new(step[0], step[1], step[2])))
    }
}

/// A line in space, its distances those to the line.
///
/// Its parameters are lengths: the shift of the axis across itself at the
/// origin, and the shift its tilt makes at `half` along it, the points'
/// greatest reach from the origin along the axis. A unit step in any of
/// them moves no point's distance by more than a unit.
#[derive(Debug, Clone, Copy)]
pub(super) struct Line {
    pub origin: Point,
    pub x: Vector,
    pub y: Vector,
    pub z: Vector,
    pub half: f64,
}

impl Line {
    /// The line through `origin` along unit `z`, scaled to `points`.
    pub(super) fn new(origin: Point, z: Vector, points: &[Point]) -> Option<Self> {
        let z = unit(z)?;
        let seed = if z.x.abs() <= z.y.abs() && z.x.abs() <= z.z.abs() {
            Vector::new(1.0, 0.0, 0.0)
        } else if z.y.abs() <= z.z.abs() {
            Vector::new(0.0, 1.0, 0.0)
        } else {
            Vector::new(0.0, 0.0, 1.0)
        };
        let x = unit(seed - z * seed.dot(z))?;
        let y = z.cross(x);
        let half = points
            .iter()
            .fold(0.0_f64, |m, p| m.max((*p - origin).dot(z).abs()));
        Some(Self {
            origin,
            x,
            y,
            z,
            half,
        })
    }
}

fn unit(v: Vector) -> Option<Vector> {
    let m = v.magnitude();
    (m > 0.0 && m.is_finite()).then(|| v / m)
}

impl Shape for Line {
    fn parameters(&self) -> usize {
        4
    }

    fn measure(&self, p: Point, gradient: &mut [f64; MOST_PARAMETERS]) -> f64 {
        let v = p - self.origin;
        let (px, py, pz) = (v.dot(self.x), v.dot(self.y), v.dot(self.z));
        let d = px.hypot(py);
        if d > 0.0 {
            let along = if self.half > 0.0 { pz / self.half } else { 0.0 };
            gradient[0] = -px / d;
            gradient[1] = -py / d;
            gradient[2] = -px / d * along;
            gradient[3] = -py / d * along;
        }
        d
    }

    fn moved(&self, step: &[f64]) -> Option<Self> {
        let origin = self.origin + self.x * step[0] + self.y * step[1];
        let tilt = if self.half > 0.0 {
            self.x * (step[2] / self.half) + self.y * (step[3] / self.half)
        } else {
            Vector::new(0.0, 0.0, 0.0)
        };
        let z = unit(self.z + tilt)?;
        let x = unit(self.x - z * self.x.dot(z))?;
        let y = z.cross(x);
        Some(Self {
            origin,
            x,
            y,
            z,
            half: self.half,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Each shape's gradient is its distance's, against central
    /// differences through `moved`.
    #[test]
    fn the_gradients_are_the_distances_own() {
        let points = [
            Point::new(3.0, 1.0, 2.0),
            Point::new(-2.0, 4.0, -1.5),
            Point::new(0.5, -3.0, 4.0),
        ];
        let line = Line::new(
            Point::new(0.2, -0.1, 0.3),
            Vector::new(0.1, 0.2, 1.0),
            &points,
        )
        .unwrap_or_else(|| panic!("a line"));
        check(&line, &points);
        check(&Centre(Point::new(0.3, 0.2, -0.4)), &points);
    }

    fn check<S: Shape>(shape: &S, points: &[Point]) {
        let h = 1e-6;
        for p in points {
            let mut g = [0.0; MOST_PARAMETERS];
            shape.measure(*p, &mut g);
            for j in 0..shape.parameters() {
                let at = |sign: f64| {
                    let mut step = vec![0.0; shape.parameters()];
                    step[j] = sign * h;
                    let moved = shape.moved(&step).unwrap_or_else(|| panic!("moved"));
                    moved.measure(*p, &mut [0.0; MOST_PARAMETERS])
                };
                let numeric = (at(1.0) - at(-1.0)) / (2.0 * h);
                assert!(
                    (numeric - g[j]).abs() < 1e-7,
                    "parameter {j}: {} against {numeric}",
                    g[j]
                );
            }
        }
    }

    #[test]
    fn largest_keeps_the_greatest_in_order() {
        let keys = [3.0, 9.0, 1.0, 9.0, 4.0];
        assert_eq!(largest(keys.len(), 3, |i| Some(keys[i])), vec![1, 3, 4]);
        assert_eq!(
            largest(keys.len(), 2, |i| (i != 1).then_some(keys[i])),
            vec![3, 4]
        );
    }
}
