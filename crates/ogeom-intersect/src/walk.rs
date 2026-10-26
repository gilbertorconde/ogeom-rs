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
//! the stall reporting and the closure test are written once and inherited.
//!
//! # What a condition owes the walker
//!
//! `n` unknowns and `n − 1` equations. That shortfall is not an oversight: the
//! solution set of `n − 1` equations in `n` unknowns *is* a curve, which is
//! what there is to follow. The walker supplies the missing equation itself
//! (a plane across the direction of travel, saying how far along to land), and
//! that is what turns "somewhere on the curve" into "the next point".
//!
//! The direction of travel comes free. The curve's tangent in parameter space
//! is the null vector of the condition's own Jacobian, and a condition that
//! has a cheaper or more careful formula for it (the intersector does, and
//! uses it to refuse a crossing too shallow to trust) says so by overriding
//! [`Condition::tangent`].

use crate::march::{Marching, Stopped};
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_math::{Point, Vector, solve};
use smallvec::SmallVec;

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

    /// How the position moves with each unknown: one vector per unknown.
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

    /// [`Condition::system`], [`Condition::position`] and
    /// [`Condition::position_gradient`] at one parameter vector, as the
    /// walker's correction asks for all three at every step. The default asks
    /// each; a condition whose position is one of the points its system
    /// already evaluates overrides this to evaluate it once.
    #[allow(clippy::type_complexity, reason = "the three answers, together")]
    fn system_at(
        &self,
        x: &[f64],
        tol: Tolerances,
    ) -> Option<((Vec<f64>, Vec<Vec<f64>>), Point, Vec<Vector>)> {
        Some((
            self.system(x, tol)?,
            self.position(x, tol)?,
            self.position_gradient(x, tol)?,
        ))
    }

    /// Bring a parameter vector back into the region the condition is posed
    /// on. Called before every evaluation, so a condition may assume it.
    fn clamp(&self, x: &mut [f64]);

    /// Whether a parameter vector has left that region.
    fn outside(&self, x: &[f64], tol: Tolerances) -> bool;

    /// Whether it is at the edge of it, which is how a stall at a boundary is
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
    /// through the tangency: two thin curves through two touching points
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
    /// more careful formula overrides this, and "more careful" is not
    /// hypothetical, since the null vector says nothing about whether the
    /// direction it found is real or is the residual's own noise.
    fn tangent(&self, x: &[f64], tol: Tolerances) -> Option<Vector> {
        let (_, jacobian) = self.system(x, tol)?;
        let gradient = self.position_gradient(x, tol)?;
        null_tangent(&jacobian, &gradient, tol)
    }

    /// [`Condition::tangent`] at `x`, given the system's Jacobian and the
    /// position gradient already evaluated there: the walker's correction
    /// ends on them, so the tangent at the point it lands on need not
    /// evaluate the condition again.
    ///
    /// The default asks [`Condition::tangent`], which is right for any
    /// condition with its own formula; one that keeps the null-space default
    /// answers with [`null_tangent`] over what it is given.
    fn tangent_from(
        &self,
        x: &[f64],
        jacobian: &[Vec<f64>],
        gradient: &[Vector],
        tol: Tolerances,
    ) -> Option<Vector> {
        let _ = (jacobian, gradient);
        self.tangent(x, tol)
    }
}

/// The unit space tangent of a condition's curve from its Jacobian (the
/// `n − 1` rows of its equations) and its position gradient: the Jacobian's
/// null vector carried into space. What [`Condition::tangent`] answers by
/// default. `None` where the Jacobian has full rank or the tangent vanishes.
#[must_use]
pub fn null_tangent(jacobian: &[Vec<f64>], gradient: &[Vector], tol: Tolerances) -> Option<Vector> {
    let null = null_vector(jacobian, gradient.len())?;
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
/// behind. Otherwise the backward half is walked: if that closes, it is the
/// loop, and otherwise the two are joined.
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
    // A backward half that closed is the whole loop as well, solved step by
    // step the other way round: the forward half stopped short of closing
    // it (out of points, or stalled where the other way passed), and joined
    // onto it would run past the start a second time. It is turned to run
    // the way the forward half set off.
    if behind.stopped == Stopped::Closed {
        let mut whole = behind;
        whole.states.reverse();
        whole.points.reverse();
        return Ok(whole);
    }

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
    // The tangent at the point just accepted, measured there to judge the
    // step's turn: the same question the next step opens with, asked with
    // the same heading, so it is answered once.
    let mut ahead: Option<Option<Vector>> = None;

    while points.len() < options.max_points {
        ogeom_core::progress::checkpoint()?;
        let direction = match ahead.take() {
            Some(known) => known,
            None => oriented(condition, &at, heading, sense, tol),
        };
        let Some(direction) = direction else {
            stopped = Stopped::Stalled;
            break;
        };
        let here = points[points.len() - 1];

        let mut taken = None;
        for _ in 0..40 {
            let Some((state, point, evaluated)) =
                correct(condition, &at, (here, direction, step), options.chord, tol)
            else {
                step *= 0.5;
                if step <= tol.confusion() {
                    break;
                }
                continue;
            };
            let next = (state, point);
            // Like against like: the *travel* direction at the next point,
            // sensed the same way, or a backward walk would read every step
            // as a half turn and crawl to a halt.
            let there = match evaluated {
                Some((jacobian, gradient)) => orient(
                    condition,
                    condition.tangent_from(&next.0, &jacobian, &gradient, tol),
                    Some(direction),
                    sense,
                ),
                None => oriented(condition, &next.0, Some(direction), sense, tol),
            };
            let turn = there.map_or(0.0, |t| direction.dot(t).clamp(-1.0, 1.0).acos());
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
                taken = Some((next, (step * scale).clamp(tol.confusion(), ceiling), there));
                break;
            }
            step *= (options.chord / sag).sqrt().clamp(0.25, 0.9);
        }
        let Some(((next_state, next_point), following, there)) = taken else {
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

        // Back at the start: a closed loop, if the step onto the start is
        // one the walk can take. Only asked once the walk has gone far
        // enough to have left, or every curve would close at once.
        if points.len() > 3
            && here.distance(from) <= step * 2.0
            && closes(
                condition,
                &at,
                (here, direction, step),
                from,
                options.chord,
                tol,
            )
        {
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
        states.push(next_state);
        points.push(next_point);
        at.clone_from(&states[states.len() - 1]);
        step = following;
        ahead = Some(there);
    }

    Ok(Walked {
        states,
        points,
        stopped,
    })
}

/// Whether the walk at `at` closes on its start: the step from `here` onto
/// `from` is solved on the condition and held to the chord like any other.
///
/// Coming back near the start is not enough. Across the seam of a closed
/// but non-periodic surface the curve can break (a spline sphere's
/// silhouette does by about 0.004, as its normals do not agree there), and
/// a walk arriving at the seam from one side finds its start a step away
/// on the other: a loop joined across that gap claims a closed curve that
/// is not there. So the start must lie within the step the walk is taking,
/// the condition solved on the plane through it must land on it, and the
/// condition solved halfway must keep to the chord between them. The step
/// is never longer than one the walk already accepted, and a start that is
/// further on is closed from the next point instead.
fn closes<C: Condition + ?Sized>(
    condition: &C,
    at: &[f64],
    (here, direction, step): (Point, Vector, f64),
    from: Point,
    chord: f64,
    tol: Tolerances,
) -> bool {
    let reach = (from - here).dot(direction);
    if reach <= 0.0 || reach > step {
        return false;
    }
    let near = chord.max(tol.confusion());
    let Some((_, landed, _)) = correct(condition, at, (here, direction, reach), chord, tol) else {
        return false;
    };
    if landed.distance(from) > near {
        return false;
    }
    let Some((_, middle, _)) = correct(condition, at, (here, direction, reach * 0.5), chord, tol)
    else {
        return false;
    };
    distance_to_segment(middle, here, from) <= near
}

/// Distance from a point to a segment.
fn distance_to_segment(p: Point, a: Point, b: Point) -> f64 {
    let along = b - a;
    let length = along.square_magnitude();
    if length <= f64::MIN_POSITIVE {
        return p.distance(a);
    }
    let t = ((p - a).dot(along) / length).clamp(0.0, 1.0);
    p.distance(a + along * t)
}

/// The tangent, turned to keep going the way the walk is going.
fn oriented<C: Condition + ?Sized>(
    condition: &C,
    at: &[f64],
    heading: Option<Vector>,
    sense: f64,
    tol: Tolerances,
) -> Option<Vector> {
    orient(condition, condition.tangent(at, tol), heading, sense)
}

/// A tangent turned to keep going the way the walk is going.
fn orient<C: Condition + ?Sized>(
    condition: &C,
    direction: Option<Vector>,
    heading: Option<Vector>,
    sense: f64,
) -> Option<Vector> {
    let direction = direction?;
    if condition.tangent_is_oriented() {
        // The condition's own sign, kept exactly, including where it flips.
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

/// Where a correction landed: the state, its point, and the condition's
/// own Jacobian rows and position gradient there when the correction's last
/// evaluation was at that state.
type Landed = (Vec<f64>, Point, Option<(Vec<Vec<f64>>, Vec<Vector>)>);

/// Bring a guess onto the condition, landing a stated distance along.
///
/// The condition's own `n − 1` equations say *on the curve*; the walker's one
/// more says *this far along it*. Without that row the system would be
/// underdetermined and Newton would wander along the curve instead of
/// converging to a point on it. The point it lands on stands within
/// `within` of the curve as far as the solve can tell: see
/// [`settling_residual`].
///
/// The sizes walked here (two to five unknowns) solve on the stack.
fn correct<C: Condition + ?Sized>(
    condition: &C,
    from: &[f64],
    travel: (Point, Vector, f64),
    within: f64,
    tol: Tolerances,
) -> Option<Landed> {
    match condition.unknowns() {
        2 => correct_fixed::<C, 2>(condition, from, travel, within, tol),
        3 => correct_fixed::<C, 3>(condition, from, travel, within, tol),
        4 => correct_fixed::<C, 4>(condition, from, travel, within, tol),
        5 => correct_fixed::<C, 5>(condition, from, travel, within, tol),
        _ => correct_any(condition, from, travel, tol),
    }
}

/// The correction's stopping rule: on the curve to a hundredth of the
/// confusion, or steps under the parametric tolerance.
fn correction_criteria(tol: Tolerances) -> solve::Criteria {
    solve::Criteria {
        residual: tol.confusion() * 0.01,
        step: tol.parametric(),
        max_iterations: 40,
    }
}

/// The residual a correction landing at a state settles to: the stopping
/// rule's own, or less where that would leave the point further than
/// `within` from the curve.
///
/// With the travel row held, a residual `r` in the condition's rows is
/// what one more Newton step would take out, and that step moves the point
/// by `G·J⁻¹·r`: `J` the system's Jacobian with the travel row, `G` the
/// position gradient. So a residual of at most `within` over the Frobenius
/// norm of `G·J⁻¹`, restricted to the condition's rows, leaves the point
/// within `within` of where the curve crosses the travel plane. Where the
/// condition's rows are well conditioned that norm is about one and the
/// stopping rule's residual is far tighter already. Where they are nearly
/// dependent it is not: two surfaces crossing at an angle `θ` part by
/// `d·sin θ` at `d` across the section, so the norm is about `1 / sin θ`,
/// and a residual of a hundredth of the confusion leaves the point
/// `ε / (100·sin θ)` across, far outside the chord at a near tangency.
pub(crate) fn settling_residual<const N: usize>(
    jacobian: &[[f64; N]; N],
    gradient: &[Vector],
    within: f64,
    stop: f64,
) -> f64 {
    let mut spread = 0.0;
    for row in 0..N - 1 {
        let mut unit = [0.0; N];
        unit[row] = 1.0;
        let Some(column) = solve::linear_fixed(*jacobian, unit) else {
            return stop;
        };
        let moved = gradient
            .iter()
            .zip(&column)
            .fold(Vector::ZERO, |sum, (g, c)| sum + *g * *c);
        spread += moved.square_magnitude();
    }
    let spread = spread.sqrt();
    if spread.is_finite() && spread * stop > within {
        within / spread
    } else {
        stop
    }
}

/// A landing solved further, down to its [`settling_residual`]: the state
/// the solve settles on, or `None` to keep the landing as it is. A landing
/// whose next Newton step, `J⁻¹·r` for its residual `r`, would move the
/// point by no more than `within` is kept without measuring the rest.
///
/// The settling residual reads the distance left to the curve off Newton's
/// linear model at the landing, and that reading holds only as far as the
/// model does. Kantorovich's condition for Newton's iteration from the
/// landing to converge to the one root beside it asks that the Jacobian's
/// change over the move, seen through the landing's inverse, stay within a
/// half. So the further solve is kept only where `‖J₀⁻¹·(J₁ − J₀)‖ ≤ ½`
/// over the move it made, in the Frobenius norm, which is never less than
/// the operator norm. Beside a touch, where two stretches of the curve
/// cross within the distance a further solve would move, the Jacobian
/// turns faster than that and the landing is kept.
pub(crate) fn settle<const N: usize>(
    mut system: impl FnMut(&[f64; N]) -> ([f64; N], [[f64; N]; N]),
    (found, residual): ([f64; N], &[f64; N]),
    (jacobian, gradient): (&[[f64; N]; N], &[Vector]),
    within: f64,
    criteria: solve::Criteria,
) -> Option<[f64; N]> {
    let step = solve::linear_fixed(*jacobian, *residual)?;
    let moved = gradient
        .iter()
        .zip(&step)
        .fold(Vector::ZERO, |sum, (g, c)| sum + *g * *c);
    if moved.magnitude() <= within {
        return None;
    }
    let norm = residual.iter().map(|v| v * v).sum::<f64>().sqrt();
    let residual = settling_residual(jacobian, gradient, within, criteria.residual);
    if norm <= residual {
        return None;
    }
    let tighter = solve::Criteria {
        residual,
        ..criteria
    };
    let (settled, left, _, _) = solve::newton_system_fixed(&mut system, found, tighter).ok()?;
    if left > norm {
        return None;
    }
    let (_, after) = system(&settled);
    let mut drift = 0.0;
    for column in 0..N {
        let change: [f64; N] =
            core::array::from_fn(|row| after[row][column] - jacobian[row][column]);
        let moved = solve::linear_fixed(*jacobian, change)?;
        drift += moved.iter().map(|v| v * v).sum::<f64>();
    }
    (drift.sqrt() <= 0.5).then_some(settled)
}

/// One evaluation of a correction: the clamped state, its point, the
/// condition's Jacobian rows, the position gradient, and the whole
/// system's Jacobian and residual.
type Evaluated<const N: usize> = (
    [f64; N],
    Point,
    Vec<Vec<f64>>,
    Vec<Vector>,
    [[f64; N]; N],
    [f64; N],
);

/// [`correct`] for `N` unknowns, allocation-free in the solve, keeping the
/// last evaluation so the landing hands back its Jacobian.
fn correct_fixed<C: Condition + ?Sized, const N: usize>(
    condition: &C,
    from: &[f64],
    (anchor, along, reach): (Point, Vector, f64),
    within: f64,
    tol: Tolerances,
) -> Option<Landed> {
    let start: [f64; N] = from.try_into().ok()?;
    let last: core::cell::Cell<Option<Evaluated<N>>> = core::cell::Cell::new(None);
    let system = |x: &[f64; N]| {
        let mut at = *x;
        condition.clamp(&mut at);
        last.set(None);
        // Where the condition cannot be evaluated the residual is infinite,
        // so the damped step backs off; a zero there would read as a root.
        let mut residual = [f64::INFINITY; N];
        let mut jacobian = [[0.0; N]; N];
        let Some(((rows, matrix), point, gradient)) = condition.system_at(&at, tol) else {
            return (residual, jacobian);
        };
        if rows.len() + 1 != N
            || matrix.len() + 1 != N
            || gradient.len() != N
            || matrix.iter().any(|row| row.len() != N)
        {
            return (residual, jacobian);
        }
        residual[..N - 1].copy_from_slice(&rows);
        for (to, row) in jacobian.iter_mut().zip(&matrix) {
            to.copy_from_slice(row);
        }
        residual[N - 1] = (point - anchor).dot(along) - reach;
        for (entry, g) in jacobian[N - 1].iter_mut().zip(&gradient) {
            *entry = g.dot(along);
        }
        last.set(Some((at, point, matrix, gradient, jacobian, residual)));
        (residual, jacobian)
    };
    let criteria = correction_criteria(tol);
    let (mut found, norm, _, _) = solve::newton_system_fixed(system, start, criteria).ok()?;
    if norm > tol.confusion() {
        return None;
    }
    let landed_at = |found: &[f64; N]| {
        let mut at = *found;
        condition.clamp(&mut at);
        at
    };
    let mut evaluated = last.take();
    if evaluated.as_ref().is_none_or(|e| e.0 != landed_at(&found)) {
        system(&found);
        evaluated = last.take();
    }
    if let Some((_, _, _, gradient, jacobian, residual)) = &evaluated
        && let Some(settled) = settle(
            system,
            (found, residual),
            (jacobian, gradient),
            within,
            criteria,
        )
    {
        found = settled;
        evaluated = last.take();
    }
    let at = landed_at(&found);
    if let Some((evaluated, point, matrix, gradient, _, _)) = evaluated
        && evaluated == at
    {
        return Some((at.to_vec(), point, Some((matrix, gradient))));
    }
    let point = condition.position(&at, tol)?;
    Some((at.to_vec(), point, None))
}

/// [`correct`] for any number of unknowns, on the general solver, to the
/// stopping rule alone: the conditions walked here are posed in two to five
/// unknowns and settle on the fixed path.
fn correct_any<C: Condition + ?Sized>(
    condition: &C,
    from: &[f64],
    (anchor, along, reach): (Point, Vector, f64),
    tol: Tolerances,
) -> Option<Landed> {
    let n = condition.unknowns();
    let system = |x: &[f64]| {
        let mut at = x.to_vec();
        condition.clamp(&mut at);
        let Some(((mut residual, mut jacobian), point, gradient)) = condition.system_at(&at, tol)
        else {
            return (vec![f64::INFINITY; n], vec![vec![0.0; n]; n]);
        };
        residual.push((point - anchor).dot(along) - reach);
        jacobian.push(gradient.iter().map(|g| g.dot(along)).collect());
        (residual, jacobian)
    };
    let found = solve::newton_system(system, from, correction_criteria(tol)).ok()?;
    if found.residual > tol.confusion() {
        return None;
    }
    let mut at = found.value;
    condition.clamp(&mut at);
    let point = condition.position(&at, tol)?;
    Some((at, point, None))
}

/// The null vector of an `(n − 1) × n` matrix: the generalized cross product.
///
/// Component `i` is the determinant of the matrix with column `i` struck out,
/// signed by `(−1)^i`, which is exactly the cross product for `n = 3` and the
/// perpendicular for `n = 2`, and is the direction the curve runs for any `n`.
/// `None` where the matrix has full rank, which means the "curve" is a point
/// and there is nothing to follow.
fn null_vector(jacobian: &[Vec<f64>], n: usize) -> Option<Vec<f64>> {
    if n == 0 || jacobian.len() + 1 != n || jacobian.iter().any(|row| row.len() < n) {
        return None;
    }
    let all: Columns = (0..n).collect();
    let mut out = Vec::with_capacity(n);
    for column in 0..n {
        let kept = without(&all, column);
        let sign = if column % 2 == 0 { 1.0 } else { -1.0 };
        out.push(sign * determinant(jacobian, &kept));
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

/// Column indices of a minor, on the stack for the sizes walked here.
type Columns = SmallVec<[usize; 8]>;

/// `columns` with the entry at `position` struck out.
fn without(columns: &[usize], position: usize) -> Columns {
    columns
        .iter()
        .enumerate()
        .filter(|(k, _)| *k != position)
        .map(|(_, c)| *c)
        .collect()
}

/// The determinant of the square minor of `rows` (from the last
/// `columns.len()` rows up) on `columns`, by expansion along its first row.
/// The sizes here are at most five.
fn determinant(rows: &[Vec<f64>], columns: &[usize]) -> f64 {
    let rows = &rows[rows.len() - columns.len()..];
    match columns.len() {
        0 => 1.0,
        1 => rows[0][columns[0]],
        2 => {
            let (a, b) = (&rows[0], &rows[1]);
            a[columns[0]].mul_add(b[columns[1]], -(a[columns[1]] * b[columns[0]]))
        }
        n => {
            let mut total = 0.0;
            for position in 0..n {
                let sign = if position % 2 == 0 { 1.0 } else { -1.0 };
                total += sign
                    * rows[0][columns[position]]
                    * determinant(&rows[1..], &without(columns, position));
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

    /// The null vector by cofactors, computed on owned minors: what the
    /// stack version must reproduce bit for bit.
    fn cofactor_null(jacobian: &[Vec<f64>], n: usize) -> Vec<f64> {
        fn det(m: &[Vec<f64>]) -> f64 {
            match m.len() {
                0 => 1.0,
                1 => m[0][0],
                2 => m[0][0].mul_add(m[1][1], -(m[0][1] * m[1][0])),
                n => (0..n)
                    .map(|c| {
                        let minor: Vec<Vec<f64>> = m[1..]
                            .iter()
                            .map(|r| (0..n).filter(|&k| k != c).map(|k| r[k]).collect())
                            .collect();
                        let sign = if c % 2 == 0 { 1.0 } else { -1.0 };
                        sign * m[0][c] * det(&minor)
                    })
                    .fold(0.0, |a, b| a + b),
            }
        }
        let mut out: Vec<f64> = (0..n)
            .map(|c| {
                let minor: Vec<Vec<f64>> = jacobian
                    .iter()
                    .map(|r| (0..n).filter(|&k| k != c).map(|k| r[k]).collect())
                    .collect();
                let sign = if c % 2 == 0 { 1.0 } else { -1.0 };
                sign * det(&minor)
            })
            .collect();
        let length = out.iter().map(|v| v * v).sum::<f64>().sqrt();
        for v in &mut out {
            *v /= length;
        }
        out
    }

    #[test]
    fn the_null_vector_matches_plain_cofactors_exactly() {
        let mut seed = 0x2545_f491_4f6c_dd1d_u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            // The top 32 bits, exact in an `f64`, in [-0.5, 0.5).
            f64::from((seed >> 32) as u32) / f64::from(u32::MAX) - 0.5
        };
        for n in 2..=5 {
            for _ in 0..50 {
                let jacobian: Vec<Vec<f64>> = (0..n - 1)
                    .map(|_| (0..n).map(|_| next()).collect())
                    .collect();
                let got = null_vector(&jacobian, n).unwrap();
                let want = cofactor_null(&jacobian, n);
                assert_eq!(
                    got.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                    want.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
                );
            }
        }
    }

    /// A circle of radius `r` about the origin in the `z = h` plane, posed in
    /// three unknowns (the point's own coordinates) with two equations. A
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
    /// and lands on the circle to the chord it was given, the tangent coming
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

    /// An ellipse with semi-axes `a` and `b` about the origin in the
    /// `z = 0` plane, posed like [`CircleAt`]: its curvature changes along
    /// it, so a walk round it takes a different number of steps each way.
    struct Ellipse {
        a: f64,
        b: f64,
    }

    impl Condition for Ellipse {
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
            let (aa, bb) = (self.a * self.a, self.b * self.b);
            Some((
                vec![
                    (x[0] * x[0]).mul_add(1.0 / aa, x[1] * x[1] / bb) - 1.0,
                    x[2],
                ],
                vec![
                    vec![2.0 * x[0] / aa, 2.0 * x[1] / bb, 0.0],
                    vec![0.0, 0.0, 1.0],
                ],
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
            self.a * 4.0
        }
    }

    /// A loop the backward half closes is the loop, whatever stopped the
    /// forward half. From a start where going round backwards takes fewer
    /// points than forwards, with points enough only for the backward way,
    /// the forward half runs out and the backward one closes: joined, the
    /// two ran past the start a second time and were reported truncated.
    /// What comes back is the backward loop, turned to run forwards.
    #[test]
    fn a_loop_closed_going_backwards_is_the_loop() {
        let ellipse = Ellipse { a: 3.0, b: 1.0 };
        let options = Marching {
            chord: 1e-5,
            ..Marching::default()
        };
        let (start, ahead, behind) = (0..64)
            .map(|k| {
                let t = f64::from(k) * 0.1;
                let start = [3.0 * t.cos(), t.sin(), 0.0];
                let ahead = walk_one_way(&ellipse, &start, 1.0, options, T).unwrap();
                let behind = walk_one_way(&ellipse, &start, -1.0, options, T).unwrap();
                (start, ahead, behind)
            })
            .find(|(_, ahead, behind)| {
                ahead.stopped == Stopped::Closed
                    && behind.stopped == Stopped::Closed
                    && behind.points.len() < ahead.points.len()
            })
            .expect("a start the backward walk closes from in fewer points");

        let options = Marching {
            max_points: behind.points.len(),
            ..options
        };
        let walked = follow(&ellipse, &start, options, T).unwrap();
        assert_eq!(walked.stopped, Stopped::Closed, "the backward half closed");
        assert_eq!(walked.points.len(), behind.points.len(), "once round");
        let (first, last) = (walked.points[0], walked.points[walked.points.len() - 1]);
        assert_eq!(first, last, "closed on its start");
        assert!(
            (walked.points[1] - first).dot(ahead.points[1] - first) > 0.0,
            "it runs the way the forward half set off"
        );
        for p in &walked.points {
            let on = (p.x / 3.0).hypot(p.y) - 1.0;
            assert!(
                on.abs() < 1e-9 && p.z.abs() < 1e-12,
                "on the ellipse: {p:?}"
            );
        }
    }

    /// A curve on a closed but non-periodic chart, the cylinder
    /// `(5 cos 2πu, 5 sin 2πu, v)` over `u` in `[0, 1]`, climbing `rise` as
    /// it goes round: it meets itself across the seam only when `rise` is
    /// zero, and otherwise breaks there by `rise`, as a spline sphere's
    /// silhouette does where the sphere's normals disagree across its seam.
    struct Climbing {
        rise: f64,
    }

    impl Climbing {
        /// How far `p` stands from the curve, at its own angle round the
        /// axis, the seam's point counted on both sides.
        fn miss(&self, p: Point) -> f64 {
            let u = p.y.atan2(p.x).rem_euclid(core::f64::consts::TAU) / core::f64::consts::TAU;
            [u, u + 1.0]
                .into_iter()
                .filter(|u| *u <= 1.0)
                .map(|u| {
                    let a = u * core::f64::consts::TAU;
                    p.distance(Point::new(5.0 * a.cos(), 5.0 * a.sin(), self.rise * u))
                })
                .fold(f64::INFINITY, f64::min)
        }
    }

    impl Condition for Climbing {
        fn unknowns(&self) -> usize {
            2
        }
        fn position(&self, x: &[f64], _tol: Tolerances) -> Option<Point> {
            let a = x[0] * core::f64::consts::TAU;
            Some(Point::new(5.0 * a.cos(), 5.0 * a.sin(), x[1]))
        }
        fn position_gradient(&self, x: &[f64], _tol: Tolerances) -> Option<Vec<Vector>> {
            let a = x[0] * core::f64::consts::TAU;
            Some(vec![
                Vector::new(-a.sin(), a.cos(), 0.0) * (5.0 * core::f64::consts::TAU),
                Vector::Z,
            ])
        }
        fn system(&self, x: &[f64], _tol: Tolerances) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
            Some((
                vec![self.rise.mul_add(-x[0], x[1])],
                vec![vec![-self.rise, 1.0]],
            ))
        }
        fn clamp(&self, x: &mut [f64]) {
            x[0] = x[0].clamp(0.0, 1.0);
            x[1] = x[1].clamp(-1.0, 1.0);
        }
        fn outside(&self, x: &[f64], tol: Tolerances) -> bool {
            let band = tol.parametric();
            x[0] < -band || x[0] > 1.0 + band || x[1] < -1.0 - band || x[1] > 1.0 + band
        }
        fn near_edge(&self, x: &[f64]) -> bool {
            x[0] <= 1e-6 || x[0] >= 1.0 - 1e-6
        }
        fn extent(&self) -> f64 {
            14.0
        }
    }

    /// A walk closes only through a step it solved. One that comes back
    /// beside its start across a seam the curve breaks at, 0.004 here as a
    /// spline sphere's silhouette breaks by about that, ends open at the
    /// seam rather than claiming a loop: joined straight onto the start,
    /// the last chord missed the curve by half the break. Every chord of
    /// what comes back keeps to the curve, from a start on either side of
    /// the seam or away from it; and where the curve does meet itself
    /// across the seam, the loop still closes.
    #[test]
    fn a_walk_closes_only_through_a_step_it_solved() {
        let chord = 1e-5;
        let options = Marching {
            chord,
            ..Marching::default()
        };
        let broken = Climbing { rise: 0.004 };
        for u in [0.0, 1e-4, 0.3, 0.9999, 1.0] {
            let walked = follow(&broken, &[u, broken.rise * u], options, T).unwrap();
            assert_ne!(
                walked.stopped,
                Stopped::Closed,
                "from {u}: a curve broken at the seam is not a loop"
            );
            for pair in walked.points.windows(2) {
                let middle = pair[0] + (pair[1] - pair[0]) * 0.5;
                let miss = broken.miss(middle);
                assert!(
                    miss <= chord,
                    "from {u}: the chord {:?} to {:?} misses the curve by {miss:e}",
                    pair[0],
                    pair[1]
                );
            }
            let (first, last) = (&walked.states[0], &walked.states[walked.states.len() - 1]);
            assert!(
                first[0].min(last[0]) <= 1e-6 && first[0].max(last[0]) >= 1.0 - 1e-6,
                "from {u}: the whole way round, seam to seam: {first:?} to {last:?}"
            );
        }

        let whole = Climbing { rise: 0.0 };
        let walked = follow(&whole, &[0.0, 0.0], options, T).unwrap();
        assert_eq!(
            walked.stopped,
            Stopped::Closed,
            "a circle met across the seam"
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
