//! One walker, several conditions.
//!
//! Following a curve nobody can write down is the same problem every time.
//! Surface intersection tracks *on both surfaces*; a silhouette tracks *the
//! normal is square to the view*; a rolling-ball blend tracks *the ball
//! touches both supports and its section stands where the guide says*. The
//! conditions are different and the geometry is different, but the walk is
//! not: take a step along the curve's own direction, correct back onto the
//! condition, measure how far the chord sagged, and set the next step from
//! that.
//!
//! So the walk lives here once, over a [`Condition`], and what changes per
//! problem is the condition's own residual and derivatives. The step control,
//! the stall reporting and the closure test are written once and inherited —
//! which matters, because they are the parts that took the longest to get
//! right and would be the easiest to get subtly wrong a second time.
//!
//! # What a condition owes the walker
//!
//! `n` unknowns and `n − 1` equations. That shortfall is not an oversight: the
//! solution set of `n − 1` equations in `n` unknowns *is* a curve, which is
//! what there is to follow. The walker supplies the missing equation itself —
//! a plane across the direction of travel, saying how far along to land — and
//! that is what turns "somewhere on the curve" into "the next point".
//!
//! The direction of travel comes free. The curve's tangent in parameter space
//! is the null vector of the condition's own Jacobian, and a condition that
//! has a cheaper or more careful formula for it — the intersector does, and
//! uses it to refuse a crossing too shallow to trust — says so by overriding
//! [`Condition::tangent`].

use crate::march::{Marching, Stopped};
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_math::{Point, Vector, solve};

/// A curve stated as what it satisfies, and everything needed to follow it.
///
/// The parameter vector is whatever the condition is posed in: four numbers
/// for a surface pair, two for a silhouette, five for a blend section marching
/// a guide. The walker never interprets them.
pub trait Condition {
    /// How many unknowns the condition is posed in.
    fn unknowns(&self) -> usize;

    /// Where a parameter vector puts the curve in space.
    fn position(&self, x: &[f64], tol: Tolerances) -> Option<Point>;

    /// How the position moves with each unknown — one vector per unknown.
    ///
    /// The walker needs this to write its own travel equation, which is a
    /// statement about where the *point* goes rather than about the
    /// parameters.
    fn position_gradient(&self, x: &[f64], tol: Tolerances) -> Option<Vec<Vector>>;

    /// The condition itself: `n − 1` residuals, and the Jacobian of them.
    ///
    /// `None` where the condition cannot be evaluated there at all, which the
    /// walker reads as a stall rather than as a zero.
    fn system(&self, x: &[f64], tol: Tolerances) -> Option<(Vec<f64>, Vec<Vec<f64>>)>;

    /// Bring a parameter vector back into the region the condition is posed
    /// on. Called before every evaluation, so a condition may assume it.
    fn clamp(&self, x: &mut [f64]);

    /// Whether a parameter vector has left that region.
    fn outside(&self, x: &[f64], tol: Tolerances) -> bool;

    /// Whether it is at the edge of it — which is how a stall at a boundary is
    /// told apart from a stall at a singularity.
    fn near_edge(&self, x: &[f64]) -> bool;

    /// A length scale for the step control: how big the thing being walked is.
    fn extent(&self) -> f64;

    /// Whether [`Condition::tangent`]'s *sign* is its own, continuous along
    /// the curve, or arbitrary from point to point.
    ///
    /// A null vector's sign is whatever the arithmetic gave it, so the default
    /// answer is no and the walker keeps its own heading. Saying yes is a
    /// claim, and a load-bearing one: where two surfaces touch, the cross
    /// product of their normals swaps sides, and a walker that quietly turned
    /// it back round would march from one branch onto the other straight
    /// through the tangency — two thin curves through two touching points
    /// coming back as one confident loop that is on neither of them. The flip
    /// is the signal, not noise.
    fn tangent_is_oriented(&self) -> bool {
        false
    }

    /// The direction the curve runs, as a unit vector in space.
    ///
    /// The default derives it from the condition's own Jacobian: the tangent
    /// in parameter space is that matrix's null vector, and the space tangent
    /// is the position gradient applied to it. A condition with a cheaper or
    /// more careful formula overrides this — and "more careful" is not
    /// hypothetical, since the null vector says nothing about whether the
    /// direction it found is real or is the residual's own noise.
    fn tangent(&self, x: &[f64], tol: Tolerances) -> Option<Vector> {
        let (_, jacobian) = self.system(x, tol)?;
        let null = null_vector(&jacobian, self.unknowns())?;
        let gradient = self.position_gradient(x, tol)?;
        let mut out = Vector::ZERO;
        for (g, n) in gradient.iter().zip(&null) {
            out += *g * *n;
        }
        let length = out.magnitude();
        if length <= tol.confusion() {
            return None;
        }
        Some(out / length)
    }
}

/// One walked curve.
#[derive(Debug, Clone, PartialEq)]
pub struct Walked {
    /// The parameter vector at each point, in order.
    pub states: Vec<Vec<f64>>,
    /// Where each is in space.
    pub points: Vec<Point>,
    /// Why it stopped.
    pub stopped: Stopped,
}

/// Follow a condition's curve both ways from a starting point.
///
/// Forwards first; if that closes, the curve is a loop and there is nothing
/// behind. Otherwise the backward half is walked and the two are joined.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// settings are unusable, or the start does not have the condition's own
/// number of unknowns.
pub fn follow<C: Condition + ?Sized>(
    condition: &C,
    start: &[f64],
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Walked> {
    options.validate()?;
    if start.len() != condition.unknowns() {
        ogeom_core::ogeom_bail!(
            Construction,
            "the condition is posed in {} unknowns and the start has {}",
            condition.unknowns(),
            start.len()
        );
    }
    let ahead = walk_one_way(condition, start, 1.0, options, tol)?;
    if ahead.stopped == Stopped::Closed {
        return Ok(ahead);
    }
    let behind = walk_one_way(condition, start, -1.0, options, tol)?;

    let mut states = behind.states;
    let mut points = behind.points;
    states.reverse();
    points.reverse();
    states.pop();
    points.pop();
    states.extend(ahead.states);
    points.extend(ahead.points);

    // The worse of the two reasons: a curve truncated at either end is
    // truncated.
    let stopped = if ahead.stopped == Stopped::RanOut || behind.stopped == Stopped::RanOut {
        Stopped::RanOut
    } else if ahead.stopped == Stopped::Stalled || behind.stopped == Stopped::Stalled {
        Stopped::Stalled
    } else {
        Stopped::LeftTheDomain
    };
    Ok(Walked {
        states,
        points,
        stopped,
    })
}

/// Walk one way from a start.
///
/// # Errors
///
/// Only through the progress sink; a walk that goes nowhere reports why in
/// [`Walked::stopped`] rather than failing.
pub fn walk_one_way<C: Condition + ?Sized>(
    condition: &C,
    start: &[f64],
    sense: f64,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Walked> {
    let mut at: Vec<f64> = start.to_vec();
    condition.clamp(&mut at);
    let Some(from) = condition.position(&at, tol) else {
        return Ok(Walked {
            states: vec![at],
            points: Vec::new(),
            stopped: Stopped::Stalled,
        });
    };
    let mut states = vec![at.clone()];
    let mut points = vec![from];
    let mut stopped = Stopped::RanOut;

    // The step is set by how far the chord may sag from the arc, and the sag
    // is measured rather than assumed: about `h · turn / 8`, where `turn` is
    // the angle between successive tangents. So the step that just meets the
    // tolerance is found by control rather than by a constant.
    let reach = condition.extent();
    let ceiling = reach / 8.0;
    let mut step = (options.chord * reach)
        .sqrt()
        .clamp(tol.confusion(), ceiling);
    // The null-space tangent's sign is arbitrary from point to point, so the
    // walk carries the direction it is going and keeps to it.
    let mut heading: Option<Vector> = None;

    while points.len() < options.max_points {
        ogeom_core::progress::checkpoint()?;
        let Some(direction) = oriented(condition, &at, heading, sense, tol) else {
            stopped = Stopped::Stalled;
            break;
        };
        let here = points[points.len() - 1];

        let mut taken = None;
        for _ in 0..40 {
            let Some(next) = correct(condition, &at, (here, direction, step), tol) else {
                step *= 0.5;
                if step <= tol.confusion() {
                    break;
                }
                continue;
            };
            // Like against like: the *travel* direction at the next point,
            // sensed the same way, or a backward walk would read every step
            // as a half turn and crawl to a halt.
            let turn = oriented(condition, &next.0, Some(direction), sense, tol)
                .map_or(0.0, |t| direction.dot(t).clamp(-1.0, 1.0).acos());
            let sag = step * turn / 8.0;
            if sag <= options.chord || step <= tol.confusion() * 8.0 {
                // Aim the next step at exactly the tolerance. Sag grows with
                // the square of the step, so the correction is a square root,
                // damped so one tight corner does not make the rest of the
                // curve expensive nor one straight stretch overshoot.
                let scale = if sag > 0.0 {
                    (options.chord / sag).sqrt().clamp(0.5, 2.0)
                } else {
                    2.0
                };
                taken = Some((next, (step * scale).clamp(tol.confusion(), ceiling)));
                break;
            }
            step *= (options.chord / sag).sqrt().clamp(0.25, 0.9);
        }
        let Some(((next_state, next_point), following)) = taken else {
            // A stall right at a domain edge is the edge, not a singularity:
            // the walk converges on the boundary from inside and the
            // correction starts failing when the step would cross it, so the
            // last accepted point sits a fraction of a step short.
            stopped = if condition.near_edge(&at) {
                Stopped::LeftTheDomain
            } else {
                Stopped::Stalled
            };
            break;
        };

        // Back where we started: a closed loop. Only checked once the walk has
        // gone far enough to have left, or every curve would close at once.
        if points.len() > 3 && next_point.distance(from) <= step {
            states.push(states[0].clone());
            points.push(from);
            stopped = Stopped::Closed;
            break;
        }
        if condition.outside(&next_state, tol) {
            stopped = Stopped::LeftTheDomain;
            break;
        }

        heading = Some(direction);
        states.push(next_state.clone());
        points.push(next_point);
        at = next_state;
        step = following;
    }

    Ok(Walked {
        states,
        points,
        stopped,
    })
}

/// The tangent, turned to keep going the way the walk is going.
fn oriented<C: Condition + ?Sized>(
    condition: &C,
    at: &[f64],
    heading: Option<Vector>,
    sense: f64,
    tol: Tolerances,
) -> Option<Vector> {
    let direction = condition.tangent(at, tol)?;
    if condition.tangent_is_oriented() {
        // The condition's own sign, kept exactly — including where it flips.
        return Some(direction * sense);
    }
    let along = match heading {
        // A null vector's sign is whatever the arithmetic gave it; what the
        // walk means by "onward" is the way it was already going.
        Some(previous) if direction.dot(previous) < 0.0 => -direction,
        _ => direction,
    };
    Some(if heading.is_none() {
        along * sense
    } else {
        along
    })
}

/// Bring a guess onto the condition, landing a stated distance along.
///
/// The condition's own `n − 1` equations say *on the curve*; the walker's one
/// more says *this far along it*. Without that row the system would be
/// underdetermined and Newton would wander along the curve instead of
/// converging to a point on it.
fn correct<C: Condition + ?Sized>(
    condition: &C,
    from: &[f64],
    (anchor, along, reach): (Point, Vector, f64),
    tol: Tolerances,
) -> Option<(Vec<f64>, Point)> {
    let n = condition.unknowns();
    let system = |x: &[f64]| {
        let mut at = x.to_vec();
        condition.clamp(&mut at);
        let (mut residual, mut jacobian) = condition
            .system(&at, tol)
            .unwrap_or_else(|| (vec![0.0; n - 1], vec![vec![0.0; n]; n - 1]));
        let point = condition.position(&at, tol).unwrap_or(Point::ORIGIN);
        let gradient = condition
            .position_gradient(&at, tol)
            .unwrap_or_else(|| vec![Vector::ZERO; n]);
        residual.push((point - anchor).dot(along) - reach);
        jacobian.push(gradient.iter().map(|g| g.dot(along)).collect());
        (residual, jacobian)
    };
    let criteria = solve::Criteria {
        residual: tol.confusion() * 0.01,
        step: tol.parametric(),
        max_iterations: 40,
    };
    let found = solve::newton_system(system, from, criteria).ok()?;
    if found.residual > tol.confusion() {
        return None;
    }
    let mut at = found.value;
    condition.clamp(&mut at);
    let point = condition.position(&at, tol)?;
    Some((at, point))
}

/// The null vector of an `(n − 1) × n` matrix: the generalized cross product.
///
/// Component `i` is the determinant of the matrix with column `i` struck out,
/// signed by `(−1)^i` — which is exactly the cross product for `n = 3` and the
/// perpendicular for `n = 2`, and is the direction the curve runs for any `n`.
/// `None` where the matrix has full rank, which means the "curve" is a point
/// and there is nothing to follow.
fn null_vector(jacobian: &[Vec<f64>], n: usize) -> Option<Vec<f64>> {
    if n == 0 || jacobian.len() + 1 != n {
        return None;
    }
    let mut out = Vec::with_capacity(n);
    for column in 0..n {
        let minor: Vec<Vec<f64>> = jacobian
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .filter(|(k, _)| *k != column)
                    .map(|(_, v)| *v)
                    .collect()
            })
            .collect();
        let sign = if column % 2 == 0 { 1.0 } else { -1.0 };
        out.push(sign * determinant(&minor));
    }
    let length = out.iter().map(|v| v * v).sum::<f64>().sqrt();
    if length <= f64::MIN_POSITIVE {
        return None;
    }
    for v in &mut out {
        *v /= length;
    }
    Some(out)
}

/// A small square determinant, by expansion. The sizes here are at most four.
fn determinant(matrix: &[Vec<f64>]) -> f64 {
    match matrix.len() {
        0 => 1.0,
        1 => matrix[0][0],
        2 => matrix[0][0].mul_add(matrix[1][1], -(matrix[0][1] * matrix[1][0])),
        n => {
            let mut total = 0.0;
            for column in 0..n {
                let minor: Vec<Vec<f64>> = matrix[1..]
                    .iter()
                    .map(|row| {
                        row.iter()
                            .enumerate()
                            .filter(|(k, _)| *k != column)
                            .map(|(_, v)| *v)
                            .collect()
                    })
                    .collect();
                let sign = if column % 2 == 0 { 1.0 } else { -1.0 };
                total += sign * matrix[0][column] * determinant(&minor);
            }
            total
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;

    const T: Tolerances = Tolerances::millimetres();

    /// A circle of radius `r` about the origin in the `z = h` plane, posed in
    /// three unknowns — the point's own coordinates — with two equations. A
    /// deliberately silly condition, chosen because its answer is known
    /// exactly and its Jacobian has nothing in common with a surface pair's.
    struct CircleAt {
        radius: f64,
        height: f64,
    }

    impl Condition for CircleAt {
        fn unknowns(&self) -> usize {
            3
        }
        fn position(&self, x: &[f64], _tol: Tolerances) -> Option<Point> {
            Some(Point::new(x[0], x[1], x[2]))
        }
        fn position_gradient(&self, _x: &[f64], _tol: Tolerances) -> Option<Vec<Vector>> {
            Some(vec![Vector::X, Vector::Y, Vector::Z])
        }
        fn system(&self, x: &[f64], _tol: Tolerances) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
            Some((
                vec![
                    x[0].mul_add(x[0], x[1] * x[1]) - self.radius * self.radius,
                    x[2] - self.height,
                ],
                vec![vec![2.0 * x[0], 2.0 * x[1], 0.0], vec![0.0, 0.0, 1.0]],
            ))
        }
        fn clamp(&self, _x: &mut [f64]) {}
        fn outside(&self, _x: &[f64], _tol: Tolerances) -> bool {
            false
        }
        fn near_edge(&self, _x: &[f64]) -> bool {
            false
        }
        fn extent(&self) -> f64 {
            self.radius * 4.0
        }
    }

    /// The walker follows a condition it has never heard of, closes the loop,
    /// and lands on the circle to the chord it was given — the tangent coming
    /// from the null space alone, since this condition supplies no formula.
    #[test]
    fn a_condition_the_walker_knows_nothing_about_is_followed_to_its_chord() {
        let circle = CircleAt {
            radius: 3.0,
            height: 1.5,
        };
        let options = Marching {
            chord: 1e-5,
            ..Marching::default()
        };
        let walked = follow(&circle, &[3.0, 0.0, 1.5], options, T).unwrap();
        assert_eq!(walked.stopped, Stopped::Closed, "a circle closes");
        assert!(walked.points.len() > 20, "{} points", walked.points.len());

        for p in &walked.points {
            assert!((p.x.hypot(p.y) - 3.0).abs() < 1e-9, "on the circle: {p:?}");
            assert!((p.z - 1.5).abs() < 1e-9, "in its plane: {p:?}");
        }
        // The polyline's length is the circumference, to the chord's own sag.
        let length: f64 = walked.points.windows(2).map(|w| w[0].distance(w[1])).sum();
        let circumference = 2.0 * core::f64::consts::PI * 3.0;
        assert!(
            length <= circumference && length > circumference * (1.0 - 1e-4),
            "the inscribed polygon: {length} against {circumference}"
        );
    }

    /// The null vector is the direction the curve runs, for the shapes a
    /// condition actually has.
    #[test]
    fn the_null_vector_is_the_generalized_cross_product() {
        // Two unknowns, one equation: the perpendicular.
        let null = null_vector(&[vec![3.0, 4.0]], 2).unwrap();
        assert!((null[0] - 0.8).abs() < 1e-12 && (null[1] + 0.6).abs() < 1e-12);
        // Three unknowns, two equations: the cross product of the rows.
        let null = null_vector(&[vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]], 3).unwrap();
        assert!(null[0].abs() < 1e-12 && null[1].abs() < 1e-12 && null[2].abs() - 1.0 < 1e-12);
        // A matrix whose rows are dependent has no curve to follow.
        assert!(null_vector(&[vec![1.0, 2.0, 3.0], vec![2.0, 4.0, 6.0]], 3).is_none());
    }
}
