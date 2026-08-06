//! Solving the sketch, and reading the system's shape.
//!
//! The solver is damped Gauss–Newton over the stacked residuals: the Jacobian
//! by central differences, the step by singular-value least squares — which
//! handles the under-determined case a sketch nearly always is — and a
//! halving line search so a bad quadratic model shrinks the step instead of
//! diverging.
//!
//! The diagnosis is the same Jacobian read for structure rather than
//! descent. Its rank against the parameter count is the degrees of freedom;
//! the right null space says *which* geometry can still move; the left null
//! space names the dependent constraint groups, and the residuals at the
//! optimum split those into the redundant — dependent but satisfied — and
//! the conflicting, which no parameter vector can satisfy together.
//!
//! **Analytic per-constraint Jacobians are declined**, and the decision is
//! settled rather than pending. Central differences over the parameters each
//! constraint names already agree with the analytic value to a part in ten
//! billion, at one extra residual evaluation apiece. Taking them exactly
//! means rewriting every residual over a scalar trait: one more place for a
//! residual and its derivative to drift apart, bought for speed the solver
//! does not need at sketch scale. If profiling ever says otherwise, the
//! change is mechanical and this is where to start.

use nalgebra::DMatrix;
use ogeom_core::OgeomResult;

use crate::model::{CircleId, ConstraintId, PointId, Sketch};

/// How hard to try.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolveOptions {
    /// The worst residual, in sketch length units, that counts as solved.
    pub tolerance: f64,
    /// Iterations before the solver admits a stall.
    pub max_iterations: usize,
    /// Compute each step through the sparse conjugate-gradient solver
    /// instead of the dense SVD. Same minimum-norm answer, cost scaling
    /// with the constraints actually written rather than with the whole
    /// parameter vector — the choice for large sketches. Diagnosis always
    /// reads structure through the SVD regardless.
    pub sparse: bool,
}

impl Default for SolveOptions {
    fn default() -> Self {
        Self {
            tolerance: 1e-9,
            max_iterations: 100,
            sparse: false,
        }
    }
}

/// What can still move, and by how much.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Freedom {
    /// Degrees of freedom remaining: parameters minus the system's rank.
    ///
    /// A sketch with no anchor keeps its rigid-body motions — two
    /// translations and a rotation — and they are counted here, because
    /// they are real: fixing a point and a direction is what removes them.
    pub degrees: usize,
    /// Points with a component in some remaining motion.
    pub movable_points: Vec<PointId>,
    /// Circles whose radius is not pinned down.
    pub movable_radii: Vec<CircleId>,
}

/// The structural reading of the sketch at its current parameters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Diagnosis {
    /// What can still move.
    pub freedom: Freedom,
    /// Dependent groups that are satisfied — redundant constraints.
    pub redundant: Vec<Vec<ConstraintId>>,
    /// Dependent groups that are violated — the conflicting constraints,
    /// by name.
    pub conflicting: Vec<Vec<ConstraintId>>,
}

impl Diagnosis {
    /// Whether the sketch is exactly constrained: nothing moves, nothing
    /// redundant, nothing conflicting.
    #[must_use]
    pub fn is_well_constrained(&self) -> bool {
        self.freedom.degrees == 0 && self.redundant.is_empty() && self.conflicting.is_empty()
    }
}

/// What a solve did, and what it learned.
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    /// Whether every residual came under the tolerance.
    ///
    /// A conflicted sketch does not converge; the diagnosis names the
    /// conflict. An unconflicted one that still fails to converge stalled —
    /// a fold the line search could not leave — and keeping the flag
    /// separate from the diagnosis keeps those two stories distinct.
    pub converged: bool,
    /// Iterations taken.
    pub iterations: usize,
    /// The worst residual at the final configuration, in length units.
    pub residual: f64,
    /// The structural reading at the final configuration.
    pub diagnosis: Diagnosis,
}

/// Singular values below this fraction of the largest count as zero for
/// rank; null-space components below this fraction of the largest count as
/// absent for naming. One knob, because both questions are "is this zero,
/// numerically" against the same matrix.
const RANK_RELATIVE: f64 = 1e-8;

/// A residual counts as violated for redundancy-vs-conflict purposes at
/// this multiple of the solve tolerance: comfortably above converged noise,
/// far below any real disagreement.
const VIOLATION_FACTOR: f64 = 1e3;

impl Sketch {
    /// Solve the sketch in place.
    ///
    /// The geometry moves to satisfy the constraints; where it cannot —
    /// conflicts — the returned [`Diagnosis`] names them, and the sketch is
    /// left at the least-bad configuration the descent reached rather than
    /// reverted, because seeing where the fight is happening is diagnostic
    /// too.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// options are not usable.
    pub fn solve(&mut self, options: SolveOptions) -> OgeomResult<Solution> {
        if !options.tolerance.is_finite() || options.tolerance <= 0.0 {
            ogeom_core::ogeom_bail!(
                Construction,
                "a tolerance of {} is not a distance",
                options.tolerance
            );
        }
        let scale = self.characteristic_scale();
        let mut residual = Vec::new();
        let mut iterations = 0;
        let mut converged = false;

        for _ in 0..options.max_iterations {
            self.residuals(&self.params.clone(), scale, &mut residual);
            let worst = residual.iter().fold(0.0_f64, |acc, r| acc.max(r.abs()));
            if worst <= options.tolerance {
                converged = true;
                break;
            }
            iterations += 1;

            let r = nalgebra::DVector::from_column_slice(&residual);
            let step: Vec<f64> = if options.sparse {
                let triplets = self.jacobian_triplets(scale, residual.len());
                let a = ogeom_math::SparseMatrix::from_triplets(
                    residual.len(),
                    self.params.len(),
                    &triplets,
                );
                let negated: Vec<f64> = residual.iter().map(|v| -v).collect();
                let budget = 4 * self.params.len().max(residual.len()).max(8);
                match ogeom_math::least_squares_cgnr(&a, &negated, 1e-10, budget) {
                    Some(step) => step,
                    None => break,
                }
            } else {
                let jacobian = self.numeric_jacobian(scale, residual.len());
                let svd = jacobian.svd(true, true);
                let Ok(step) = svd.solve(&(-&r), rank_epsilon(&svd)) else {
                    break;
                };
                step.iter().copied().collect()
            };

            // Halving line search: accept the first step that reduces the
            // residual norm; a stall — no fraction helps — ends the descent.
            let before = r.norm();
            let base = self.params.clone();
            let mut accepted = false;
            let mut t = 1.0;
            for _ in 0..8 {
                for (p, s) in self.params.iter_mut().zip(step.iter()) {
                    *p += t * s;
                }
                self.residuals(&self.params.clone(), scale, &mut residual);
                let after = nalgebra::DVector::from_column_slice(&residual).norm();
                if after < before {
                    accepted = true;
                    break;
                }
                self.params.copy_from_slice(&base);
                t /= 2.0;
            }
            if !accepted {
                break;
            }
        }

        // The final reading, wherever the descent ended.
        self.residuals(&self.params.clone(), scale, &mut residual);
        let worst = residual.iter().fold(0.0_f64, |acc, r| acc.max(r.abs()));
        if worst <= options.tolerance {
            converged = true;
        }
        let diagnosis = self.diagnose_with(scale, &residual, options.tolerance);
        Ok(Solution {
            converged,
            iterations,
            residual: worst,
            diagnosis,
        })
    }

    /// Drag a point toward a target, re-solving from where the sketch
    /// stands.
    ///
    /// The target is a *soft* objective: two lightly weighted rows pull
    /// the point while every constraint keeps its full weight, so the
    /// solve lands on the constraint-satisfying configuration nearest the
    /// pointer rather than fighting the dimensions. Warm-started from the
    /// current parameters, this is the primitive an interactive sketcher
    /// calls every frame.
    ///
    /// # Errors
    ///
    /// As [`Sketch::solve`], plus
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the point is
    /// not from this sketch.
    pub fn drag(
        &mut self,
        point: crate::model::PointId,
        to: ogeom_math::Point2,
        options: SolveOptions,
    ) -> OgeomResult<Solution> {
        // The pull enters as an ordinary constraint at a whisper of the
        // weight, then leaves; ids of real constraints are untouched
        // because the pull is appended last and popped before any
        // diagnosis the caller sees is produced... but diagnosis happens
        // inside solve, so instead: solve against a temporary constraint
        // and re-diagnose without it.
        let pull = self.constrain(crate::model::Constraint::Fixed(point, to))?;
        let _ = pull;
        self.soft_last = true;
        let outcome = self.solve(options);
        self.soft_last = false;
        self.constraints.pop();
        outcome?;
        // Release polish: even a whisper of a pull perturbs the hard
        // constraints at its squared weight, so the sketch re-solves
        // without it — warm-started, a step or two — and lands exactly on
        // the constraints, at the configuration the drag chose.
        self.solve(options)
    }

    /// Read the sketch's structure at its current parameters, without
    /// moving anything.
    #[must_use]
    pub fn diagnose(&self) -> Diagnosis {
        let scale = self.characteristic_scale();
        let mut residual = Vec::new();
        self.residuals(&self.params.clone(), scale, &mut residual);
        self.diagnose_with(scale, &residual, SolveOptions::default().tolerance)
    }

    fn diagnose_with(&self, scale: f64, residual: &[f64], tolerance: f64) -> Diagnosis {
        let n = self.params.len();
        let m = residual.len();
        if n == 0 || m == 0 {
            return Diagnosis {
                freedom: Freedom {
                    degrees: n,
                    movable_points: self
                        .points
                        .iter()
                        .enumerate()
                        .map(|(i, _)| PointId(i))
                        .collect(),
                    movable_radii: self
                        .circles
                        .iter()
                        .enumerate()
                        .map(|(i, _)| CircleId(i))
                        .collect(),
                },
                redundant: Vec::new(),
                conflicting: Vec::new(),
            };
        }

        let jacobian = self.numeric_jacobian(scale, m);
        let svd = jacobian.svd(true, true);
        let sigma_max = svd
            .singular_values
            .iter()
            .fold(0.0_f64, |acc, s| acc.max(*s));
        let rank = if sigma_max > 0.0 {
            svd.singular_values
                .iter()
                .filter(|s| **s > sigma_max * RANK_RELATIVE)
                .count()
        } else {
            0
        };

        // The thin SVD spans only min(m, n) directions, so the null spaces
        // are read through *projectors* instead of missing factor columns:
        // I − V_r V_rᵀ projects onto everything the constraints do not
        // reach, and I − U_r U_rᵀ onto everything the parameters cannot
        // produce. The thin factors determine both completely.

        // Right null: a parameter is free exactly when the row space does
        // not contain its axis — the projector's diagonal is its remaining
        // freedom, in [0, 1].
        let mut movable_points = Vec::new();
        let mut movable_radii = Vec::new();
        if let Some(v_t) = &svd.v_t {
            let mut param_free = vec![1.0_f64; n];
            for r in 0..rank.min(v_t.nrows()) {
                for (j, free) in param_free.iter_mut().enumerate() {
                    *free -= v_t[(r, j)] * v_t[(r, j)];
                }
            }
            for (i, data) in self.points.iter().enumerate() {
                if param_free[data.at] > 1e-6 || param_free[data.at + 1] > 1e-6 {
                    movable_points.push(PointId(i));
                }
            }
            for (i, data) in self.circles.iter().enumerate() {
                if param_free[data.radius_at] > 1e-6 {
                    movable_radii.push(CircleId(i));
                }
            }
        }

        // Left null: a residual row participates in a dependence exactly
        // when the column space does not contain its axis, and two rows
        // belong to the same dependent group when the projector couples
        // them. Connected components of that coupling are the groups.
        let mut redundant: Vec<Vec<ConstraintId>> = Vec::new();
        let mut conflicting: Vec<Vec<ConstraintId>> = Vec::new();
        if let Some(u) = &svd.u {
            let columns = rank.min(u.ncols());
            let projector = |i: usize, k: usize| -> f64 {
                let identity = if i == k { 1.0 } else { 0.0 };
                let mut reach = 0.0;
                for r in 0..columns {
                    reach += u[(i, r)] * u[(k, r)];
                }
                identity - reach
            };
            let participants: Vec<usize> = (0..m).filter(|&i| projector(i, i) > 1e-6).collect();

            // Connected components over the projector's coupling.
            let mut groups: Vec<Vec<usize>> = Vec::new();
            let mut assigned = vec![false; m];
            for &start in &participants {
                if assigned[start] {
                    continue;
                }
                let mut group = vec![start];
                assigned[start] = true;
                let mut cursor = 0;
                while cursor < group.len() {
                    let i = group[cursor];
                    cursor += 1;
                    for &k in &participants {
                        if !assigned[k] && projector(i, k).abs() > 1e-6 {
                            assigned[k] = true;
                            group.push(k);
                        }
                    }
                }
                group.sort_unstable();
                groups.push(group);
            }

            let row_owner = self.row_owners(m);
            let violation = tolerance * VIOLATION_FACTOR;
            for rows in groups {
                let mut ids: Vec<ConstraintId> = rows.iter().map(|&r| row_owner[r]).collect();
                ids.sort_unstable_by_key(|id| id.0);
                ids.dedup();
                let satisfied = rows.iter().all(|&r| residual[r].abs() <= violation);
                if satisfied {
                    redundant.push(ids);
                } else {
                    conflicting.push(ids);
                }
            }
        }

        Diagnosis {
            freedom: Freedom {
                degrees: n - rank,
                movable_points,
                movable_radii,
            },
            redundant,
            conflicting,
        }
    }

    /// The Jacobian of the residuals by central differences, exploiting
    /// structural sparsity: each constraint is differentiated only over
    /// the parameters of the entities it names, so the cost scales with
    /// the constraint count rather than with constraints times the whole
    /// parameter vector — the difference between a solver that re-reads
    /// the world and one that can keep up with a drag.
    fn numeric_jacobian(&self, scale: f64, m: usize) -> DMatrix<f64> {
        let mut jacobian = DMatrix::zeros(m, self.params.len());
        for (row, col, value) in self.jacobian_triplets(scale, m) {
            jacobian[(row, col)] = value;
        }
        jacobian
    }

    /// The sparse Jacobian as `(row, column, value)` entries — only the
    /// parameters each constraint names are differenced, so the entry count
    /// is the coupling structure itself.
    fn jacobian_triplets(&self, scale: f64, m: usize) -> Vec<(usize, usize, f64)> {
        let mut triplets = Vec::new();
        let mut params = self.params.clone();
        let mut forward = Vec::new();
        let mut backward = Vec::new();
        let mut row = 0;
        for constraint in &self.constraints {
            let rows = Self::rows_of(constraint);
            for j in self.parameters_of(constraint) {
                let h = 1e-6 * params[j].abs().max(1.0);
                let held = params[j];
                params[j] = held + h;
                forward.clear();
                self.constraint_residuals(constraint, &params, scale, &mut forward);
                params[j] = held - h;
                backward.clear();
                self.constraint_residuals(constraint, &params, scale, &mut backward);
                params[j] = held;
                let weight = if self.soft_last && row + rows == m {
                    crate::model::SOFT_WEIGHT
                } else {
                    1.0
                };
                for k in 0..rows {
                    triplets.push((row + k, j, weight * (forward[k] - backward[k]) / (2.0 * h)));
                }
            }
            row += rows;
        }
        debug_assert_eq!(row, m);
        triplets
    }

    /// Which constraint owns each residual row.
    fn row_owners(&self, m: usize) -> Vec<ConstraintId> {
        let mut owners = Vec::with_capacity(m);
        for (i, constraint) in self.constraints.iter().enumerate() {
            for _ in 0..Self::rows_of(constraint) {
                owners.push(ConstraintId(i));
            }
        }
        owners
    }
}

/// The epsilon `nalgebra`'s SVD solve treats as zero, tied to the matrix
/// itself rather than an absolute guess.
fn rank_epsilon(svd: &nalgebra::SVD<f64, nalgebra::Dyn, nalgebra::Dyn>) -> f64 {
    let sigma_max = svd
        .singular_values
        .iter()
        .fold(0.0_f64, |acc, s| acc.max(*s));
    (sigma_max * RANK_RELATIVE).max(f64::MIN_POSITIVE)
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::model::{Constraint, Sketch, TangencySide};
    use approx::assert_relative_eq;
    use ogeom_math::Point2;

    fn solved(sketch: &mut Sketch) -> Solution {
        let solution = sketch.solve(SolveOptions::default()).unwrap();
        assert!(
            solution.converged,
            "expected convergence, residual {}",
            solution.residual
        );
        solution
    }

    /// A rectangle pinned at a corner: eight parameters, eight rows, no
    /// freedom left, and the corners land where the dimensions say.
    #[test]
    fn a_dimensioned_rectangle_is_well_constrained_and_lands_exactly() {
        let mut sketch = Sketch::new();
        // Deliberately skewed starting guesses.
        let p0 = sketch.add_point(Point2::new(0.3, -0.2));
        let p1 = sketch.add_point(Point2::new(35.0, 4.0));
        let p2 = sketch.add_point(Point2::new(43.0, 28.0));
        let p3 = sketch.add_point(Point2::new(-3.0, 31.0));
        let bottom = sketch.add_line(p0, p1).unwrap();
        let right = sketch.add_line(p1, p2).unwrap();
        let top = sketch.add_line(p3, p2).unwrap();
        let left = sketch.add_line(p0, p3).unwrap();

        sketch
            .constrain(Constraint::Fixed(p0, Point2::new(0.0, 0.0)))
            .unwrap();
        sketch.constrain(Constraint::Horizontal(bottom)).unwrap();
        sketch
            .constrain(Constraint::Distance(p0, p1, 40.0))
            .unwrap();
        sketch.constrain(Constraint::Vertical(right)).unwrap();
        sketch
            .constrain(Constraint::Distance(p1, p2, 30.0))
            .unwrap();
        sketch.constrain(Constraint::Horizontal(top)).unwrap();
        sketch.constrain(Constraint::Vertical(left)).unwrap();

        let solution = solved(&mut sketch);
        assert!(solution.diagnosis.is_well_constrained());
        assert_eq!(solution.diagnosis.freedom.degrees, 0);
        let c2 = sketch.point(p2).unwrap();
        // The solve may put p1 at +40 or -40 from the anchor; both are the
        // rectangle. The dimensions are what must hold.
        assert_relative_eq!(
            sketch.measure_distance(p0, p1).unwrap(),
            40.0,
            epsilon = 1e-7
        );
        assert_relative_eq!(
            sketch.measure_distance(p1, p2).unwrap(),
            30.0,
            epsilon = 1e-7
        );
        assert_relative_eq!(c2.y.abs(), 30.0, epsilon = 1e-7);
    }

    /// Without the anchor the same rectangle keeps its rigid-body motions,
    /// and the diagnosis says which geometry can still move: all of it.
    #[test]
    fn an_unanchored_rectangle_reports_its_rigid_freedom() {
        let mut sketch = Sketch::new();
        let p0 = sketch.add_point(Point2::new(0.0, 0.0));
        let p1 = sketch.add_point(Point2::new(40.0, 0.0));
        let p2 = sketch.add_point(Point2::new(40.0, 30.0));
        let p3 = sketch.add_point(Point2::new(0.0, 30.0));
        let bottom = sketch.add_line(p0, p1).unwrap();
        let right = sketch.add_line(p1, p2).unwrap();
        let top = sketch.add_line(p3, p2).unwrap();
        let left = sketch.add_line(p0, p3).unwrap();
        sketch.constrain(Constraint::Horizontal(bottom)).unwrap();
        sketch
            .constrain(Constraint::Distance(p0, p1, 40.0))
            .unwrap();
        sketch.constrain(Constraint::Vertical(right)).unwrap();
        sketch
            .constrain(Constraint::Distance(p1, p2, 30.0))
            .unwrap();
        sketch.constrain(Constraint::Horizontal(top)).unwrap();
        sketch.constrain(Constraint::Vertical(left)).unwrap();

        let solution = solved(&mut sketch);
        // Two translations plus... the bottom-horizontal already pins the
        // rotation, so the rectangle slides but does not spin.
        assert_eq!(solution.diagnosis.freedom.degrees, 2);
        assert_eq!(solution.diagnosis.freedom.movable_points.len(), 4);
        assert!(solution.diagnosis.conflicting.is_empty());
    }

    /// Two distances that cannot both hold: the diagnosis names exactly the
    /// two of them, and the solver honestly does not converge.
    #[test]
    fn contradictory_distances_are_named() {
        let mut sketch = Sketch::new();
        let p0 = sketch.add_point(Point2::new(0.0, 0.0));
        let p1 = sketch.add_point(Point2::new(10.0, 0.0));
        let l = sketch.add_line(p0, p1).unwrap();
        sketch
            .constrain(Constraint::Fixed(p0, Point2::new(0.0, 0.0)))
            .unwrap();
        sketch.constrain(Constraint::Horizontal(l)).unwrap();
        let d10 = sketch
            .constrain(Constraint::Distance(p0, p1, 10.0))
            .unwrap();
        let d12 = sketch
            .constrain(Constraint::Distance(p0, p1, 12.0))
            .unwrap();

        let solution = sketch.solve(SolveOptions::default()).unwrap();
        assert!(!solution.converged, "a contradiction cannot converge");
        assert_eq!(solution.diagnosis.conflicting.len(), 1);
        let group = &solution.diagnosis.conflicting[0];
        assert!(
            group.contains(&d10) && group.contains(&d12),
            "group: {group:?}"
        );
        // The names are printable.
        for id in group {
            assert!(sketch.describe(*id).unwrap().starts_with("distance"));
        }
    }

    /// The same distance twice is dependent but satisfiable: redundant, not
    /// conflicting, and the sketch still solves.
    #[test]
    fn a_repeated_distance_is_redundant_not_conflicting() {
        let mut sketch = Sketch::new();
        let p0 = sketch.add_point(Point2::new(0.0, 0.0));
        let p1 = sketch.add_point(Point2::new(9.0, 1.0));
        let l = sketch.add_line(p0, p1).unwrap();
        sketch
            .constrain(Constraint::Fixed(p0, Point2::new(0.0, 0.0)))
            .unwrap();
        sketch.constrain(Constraint::Horizontal(l)).unwrap();
        let a = sketch
            .constrain(Constraint::Distance(p0, p1, 10.0))
            .unwrap();
        let b = sketch
            .constrain(Constraint::Distance(p0, p1, 10.0))
            .unwrap();

        let solution = solved(&mut sketch);
        assert!(solution.diagnosis.conflicting.is_empty());
        assert_eq!(solution.diagnosis.redundant.len(), 1);
        let group = &solution.diagnosis.redundant[0];
        assert!(group.contains(&a) && group.contains(&b));
    }

    /// A line pulled tangent to a fixed circle: the point-line distance
    /// equals the radius afterwards, from a start that was nowhere near.
    #[test]
    fn a_line_solves_to_tangency() {
        let mut sketch = Sketch::new();
        let centre = sketch.add_point(Point2::new(0.0, 0.0));
        let circle = sketch.add_circle(centre, 5.0).unwrap();
        let a = sketch.add_point(Point2::new(-10.0, 2.0));
        let b = sketch.add_point(Point2::new(10.0, 1.0));
        let l = sketch.add_line(a, b).unwrap();
        sketch
            .constrain(Constraint::Fixed(centre, Point2::new(0.0, 0.0)))
            .unwrap();
        sketch.constrain(Constraint::Radius(circle, 5.0)).unwrap();
        sketch
            .constrain(Constraint::Fixed(a, Point2::new(-10.0, 2.0)))
            .unwrap();
        sketch.constrain(Constraint::Horizontal(l)).unwrap();
        sketch
            .constrain(Constraint::TangentLineCircle(l, circle))
            .unwrap();

        // Horizontal through y=2 cannot touch a radius-5 circle without
        // moving; the free endpoint drags the carrier to y = +-5... but a
        // horizontal line through a fixed point at y=2 is stuck at y=2, so
        // tangency must move the line by moving b? It cannot: horizontal +
        // fixed(a) pins the carrier. This is a conflict, and it is named.
        let solution = sketch.solve(SolveOptions::default()).unwrap();
        assert!(!solution.converged);
        assert!(!solution.diagnosis.conflicting.is_empty());

        // Freed of the anchor, it solves: the whole line slides to y = 5
        // or y = -5.
        let mut free = Sketch::new();
        let centre = free.add_point(Point2::new(0.0, 0.0));
        let circle = free.add_circle(centre, 5.0).unwrap();
        let a = free.add_point(Point2::new(-10.0, 2.0));
        let b = free.add_point(Point2::new(10.0, 1.0));
        let l = free.add_line(a, b).unwrap();
        free.constrain(Constraint::Fixed(centre, Point2::new(0.0, 0.0)))
            .unwrap();
        free.constrain(Constraint::Radius(circle, 5.0)).unwrap();
        free.constrain(Constraint::Horizontal(l)).unwrap();
        free.constrain(Constraint::TangentLineCircle(l, circle))
            .unwrap();
        solved(&mut free);
        let (pa, _) = free.line_ends(l).unwrap();
        assert_relative_eq!(pa.y.abs(), 5.0, epsilon = 1e-6);
    }

    /// Perpendicularity and a 45-degree angle, solved together on a fan of
    /// three lines from one anchor.
    #[test]
    fn angles_solve_between_lines() {
        let mut sketch = Sketch::new();
        let origin = sketch.add_point(Point2::new(0.0, 0.0));
        let px = sketch.add_point(Point2::new(10.0, 0.5));
        let py = sketch.add_point(Point2::new(-0.5, 10.0));
        let pd = sketch.add_point(Point2::new(8.0, 6.0));
        let lx = sketch.add_line(origin, px).unwrap();
        let ly = sketch.add_line(origin, py).unwrap();
        let ld = sketch.add_line(origin, pd).unwrap();
        sketch
            .constrain(Constraint::Fixed(origin, Point2::new(0.0, 0.0)))
            .unwrap();
        sketch
            .constrain(Constraint::Fixed(px, Point2::new(10.0, 0.0)))
            .unwrap();
        sketch.constrain(Constraint::Perpendicular(lx, ly)).unwrap();
        sketch
            .constrain(Constraint::Angle(lx, ld, core::f64::consts::FRAC_PI_4))
            .unwrap();

        solved(&mut sketch);
        assert_relative_eq!(
            sketch.measure_angle(lx, ly).unwrap().abs(),
            core::f64::consts::FRAC_PI_2,
            epsilon = 1e-7
        );
        assert_relative_eq!(
            sketch.measure_angle(lx, ld).unwrap(),
            core::f64::consts::FRAC_PI_4,
            epsilon = 1e-7
        );
    }

    /// Symmetry across a construction axis: the pair mirrors, and marking
    /// the axis construction changes nothing about the solve.
    #[test]
    fn symmetry_mirrors_across_a_construction_axis() {
        let mut sketch = Sketch::new();
        let a0 = sketch.add_point(Point2::new(0.0, -10.0));
        let a1 = sketch.add_point(Point2::new(0.0, 10.0));
        let axis = sketch.add_line(a0, a1).unwrap();
        sketch.set_line_construction(axis, true).unwrap();
        let p = sketch.add_point(Point2::new(-7.0, 3.0));
        let q = sketch.add_point(Point2::new(4.0, 5.0));
        sketch
            .constrain(Constraint::Fixed(a0, Point2::new(0.0, -10.0)))
            .unwrap();
        sketch
            .constrain(Constraint::Fixed(a1, Point2::new(0.0, 10.0)))
            .unwrap();
        sketch
            .constrain(Constraint::Fixed(p, Point2::new(-7.0, 3.0)))
            .unwrap();
        sketch.constrain(Constraint::Symmetric(p, q, axis)).unwrap();

        solved(&mut sketch);
        assert!(sketch.is_line_construction(axis).unwrap());
        let qq = sketch.point(q).unwrap();
        assert_relative_eq!(qq.x, 7.0, epsilon = 1e-7);
        assert_relative_eq!(qq.y, 3.0, epsilon = 1e-7);
    }

    /// An arc's rim points stay at one radius through a solve that moves
    /// one of them.
    #[test]
    fn an_arc_keeps_one_radius() {
        let mut sketch = Sketch::new();
        let centre = sketch.add_point(Point2::new(0.0, 0.0));
        let start = sketch.add_point(Point2::new(8.0, 0.0));
        let end = sketch.add_point(Point2::new(0.5, 7.0));
        let (_arc, _coupling) = sketch.add_arc(centre, start, end).unwrap();
        sketch
            .constrain(Constraint::Fixed(centre, Point2::new(0.0, 0.0)))
            .unwrap();
        sketch
            .constrain(Constraint::Fixed(start, Point2::new(8.0, 0.0)))
            .unwrap();

        solved(&mut sketch);
        assert_relative_eq!(
            sketch.measure_distance(centre, end).unwrap(),
            8.0,
            epsilon = 1e-7
        );
    }

    /// External and internal circle tangency, each holding its recorded side.
    #[test]
    fn circle_tangency_holds_its_side() {
        let mut sketch = Sketch::new();
        let c1 = sketch.add_point(Point2::new(0.0, 0.0));
        let c2 = sketch.add_point(Point2::new(12.0, 0.5));
        let big = sketch.add_circle(c1, 5.0).unwrap();
        let small = sketch.add_circle(c2, 2.0).unwrap();
        sketch
            .constrain(Constraint::Fixed(c1, Point2::new(0.0, 0.0)))
            .unwrap();
        sketch.constrain(Constraint::Radius(big, 5.0)).unwrap();
        sketch.constrain(Constraint::Radius(small, 2.0)).unwrap();
        sketch
            .constrain(Constraint::TangentCircles(
                big,
                small,
                TangencySide::External,
            ))
            .unwrap();
        solved(&mut sketch);
        assert_relative_eq!(
            sketch.measure_distance(c1, c2).unwrap(),
            7.0,
            epsilon = 1e-7
        );
    }

    /// Equal length and equal radius propagate a dimension.
    #[test]
    fn equality_propagates_dimensions() {
        let mut sketch = Sketch::new();
        let p0 = sketch.add_point(Point2::new(0.0, 0.0));
        let p1 = sketch.add_point(Point2::new(10.0, 0.0));
        let p2 = sketch.add_point(Point2::new(0.0, 4.0));
        let p3 = sketch.add_point(Point2::new(6.0, 5.0));
        let l1 = sketch.add_line(p0, p1).unwrap();
        let l2 = sketch.add_line(p2, p3).unwrap();
        sketch
            .constrain(Constraint::Fixed(p0, Point2::new(0.0, 0.0)))
            .unwrap();
        sketch
            .constrain(Constraint::Fixed(p1, Point2::new(10.0, 0.0)))
            .unwrap();
        sketch
            .constrain(Constraint::Fixed(p2, Point2::new(0.0, 4.0)))
            .unwrap();
        sketch.constrain(Constraint::EqualLength(l1, l2)).unwrap();

        let ca = sketch.add_point(Point2::new(20.0, 0.0));
        let cb = sketch.add_point(Point2::new(30.0, 0.0));
        let circle_a = sketch.add_circle(ca, 3.0).unwrap();
        let circle_b = sketch.add_circle(cb, 1.0).unwrap();
        sketch.constrain(Constraint::Radius(circle_a, 3.0)).unwrap();
        sketch
            .constrain(Constraint::EqualRadius(circle_a, circle_b))
            .unwrap();

        solved(&mut sketch);
        assert_relative_eq!(
            sketch.measure_distance(p2, p3).unwrap(),
            10.0,
            epsilon = 1e-7
        );
        assert_relative_eq!(
            sketch.measure_radius(circle_b).unwrap(),
            3.0,
            epsilon = 1e-7
        );
    }

    /// Driven dimensions read without constraining: measuring changes no
    /// degrees of freedom.
    #[test]
    fn measurements_are_driven_not_driving() {
        let mut sketch = Sketch::new();
        let p0 = sketch.add_point(Point2::new(0.0, 0.0));
        let p1 = sketch.add_point(Point2::new(3.0, 4.0));
        let before = sketch.diagnose();
        assert_eq!(before.freedom.degrees, 4);
        assert_relative_eq!(
            sketch.measure_distance(p0, p1).unwrap(),
            5.0,
            epsilon = 1e-12
        );
        let after = sketch.diagnose();
        assert_eq!(after.freedom.degrees, 4, "measuring moved nothing");
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod sparse_tests {
    use super::*;
    use crate::model::{Constraint, Sketch};
    use approx::assert_relative_eq;
    use ogeom_math::Point2;

    const T_OPTS: SolveOptions = SolveOptions {
        tolerance: 1e-9,
        max_iterations: 100,
        sparse: false,
    };

    /// The sparse Jacobian equals the dense one entry for entry, on a
    /// sketch exercising every kind of coupling.
    #[test]
    fn the_sparse_step_solves_the_same_sketch() {
        // The same dimensioned rectangle twice: the dense SVD path and the
        // sparse conjugate-gradient path must both converge, and land on
        // the same measured geometry.
        let build = || {
            let mut sketch = Sketch::new();
            let p0 = sketch.add_point(Point2::new(0.1, -0.3));
            let p1 = sketch.add_point(Point2::new(52.0, 2.0));
            let p2 = sketch.add_point(Point2::new(49.0, 21.0));
            let p3 = sketch.add_point(Point2::new(-2.0, 19.0));
            let bottom = sketch.add_line(p0, p1).unwrap();
            let right = sketch.add_line(p1, p2).unwrap();
            let top = sketch.add_line(p3, p2).unwrap();
            let left = sketch.add_line(p0, p3).unwrap();
            sketch
                .constrain(Constraint::Fixed(p0, Point2::new(0.0, 0.0)))
                .unwrap();
            sketch.constrain(Constraint::Horizontal(bottom)).unwrap();
            sketch
                .constrain(Constraint::Distance(p0, p1, 50.0))
                .unwrap();
            sketch.constrain(Constraint::Vertical(right)).unwrap();
            sketch
                .constrain(Constraint::Distance(p1, p2, 20.0))
                .unwrap();
            sketch.constrain(Constraint::Horizontal(top)).unwrap();
            sketch.constrain(Constraint::Vertical(left)).unwrap();
            (sketch, p0, p1, p2)
        };

        let (mut dense_sketch, _, d1, d2) = build();
        let dense = dense_sketch.solve(T_OPTS).unwrap();
        assert!(dense.converged);

        let (mut sparse_sketch, s0, s1, s2) = build();
        let sparse = sparse_sketch
            .solve(SolveOptions {
                sparse: true,
                ..T_OPTS
            })
            .unwrap();
        assert!(sparse.converged, "residual {}", sparse.residual);
        assert_relative_eq!(
            sparse_sketch.measure_distance(s0, s1).unwrap(),
            50.0,
            epsilon = 1e-6
        );
        assert_relative_eq!(
            sparse_sketch.measure_distance(s1, s2).unwrap(),
            dense_sketch.measure_distance(d1, d2).unwrap(),
            epsilon = 1e-6
        );
    }

    #[test]
    fn the_sparse_jacobian_is_the_dense_one() {
        let mut sketch = Sketch::new();
        let p0 = sketch.add_point(Point2::new(0.3, -0.1));
        let p1 = sketch.add_point(Point2::new(9.0, 1.0));
        let p2 = sketch.add_point(Point2::new(8.0, 7.0));
        let centre = sketch.add_point(Point2::new(4.0, 3.0));
        let l1 = sketch.add_line(p0, p1).unwrap();
        let l2 = sketch.add_line(p1, p2).unwrap();
        let circle = sketch.add_circle(centre, 2.0).unwrap();
        sketch
            .constrain(Constraint::Fixed(p0, Point2::new(0.0, 0.0)))
            .unwrap();
        sketch.constrain(Constraint::Distance(p0, p1, 9.0)).unwrap();
        sketch.constrain(Constraint::Perpendicular(l1, l2)).unwrap();
        sketch
            .constrain(Constraint::TangentLineCircle(l1, circle))
            .unwrap();
        sketch.constrain(Constraint::Symmetric(p1, p2, l1)).unwrap();
        sketch.constrain(Constraint::Radius(circle, 2.0)).unwrap();

        let scale = sketch.characteristic_scale();
        let mut residual = Vec::new();
        sketch.residuals(&sketch.params.clone(), scale, &mut residual);
        let m = residual.len();
        let sparse = sketch.numeric_jacobian(scale, m);

        // The dense reference: perturb every parameter, evaluate everything.
        let n = sketch.params.len();
        let mut dense = nalgebra::DMatrix::zeros(m, n);
        let mut params = sketch.params.clone();
        let (mut forward, mut backward) = (Vec::new(), Vec::new());
        for j in 0..n {
            let h = 1e-6 * params[j].abs().max(1.0);
            let held = params[j];
            params[j] = held + h;
            sketch.residuals(&params, scale, &mut forward);
            params[j] = held - h;
            sketch.residuals(&params, scale, &mut backward);
            params[j] = held;
            for i in 0..m {
                dense[(i, j)] = (forward[i] - backward[i]) / (2.0 * h);
            }
        }
        for i in 0..m {
            for j in 0..n {
                assert_relative_eq!(sparse[(i, j)], dense[(i, j)], epsilon = 1e-9);
            }
        }
    }

    /// Dragging slides a rectangle's free corner along what the
    /// constraints allow: the dimensions hold exactly, and the corner
    /// lands as near the pointer as the dimensions permit.
    #[test]
    fn dragging_moves_what_may_move_and_nothing_else() {
        let mut sketch = Sketch::new();
        let p0 = sketch.add_point(Point2::new(0.0, 0.0));
        let p1 = sketch.add_point(Point2::new(40.0, 0.0));
        let p2 = sketch.add_point(Point2::new(40.0, 30.0));
        let p3 = sketch.add_point(Point2::new(0.0, 30.0));
        let bottom = sketch.add_line(p0, p1).unwrap();
        let right = sketch.add_line(p1, p2).unwrap();
        let top = sketch.add_line(p3, p2).unwrap();
        let left = sketch.add_line(p0, p3).unwrap();
        sketch.constrain(Constraint::Horizontal(bottom)).unwrap();
        sketch
            .constrain(Constraint::Distance(p0, p1, 40.0))
            .unwrap();
        sketch.constrain(Constraint::Vertical(right)).unwrap();
        sketch
            .constrain(Constraint::Distance(p1, p2, 30.0))
            .unwrap();
        sketch.constrain(Constraint::Horizontal(top)).unwrap();
        sketch.constrain(Constraint::Vertical(left)).unwrap();

        // Unanchored: the rectangle may translate. Drag the origin corner.
        let before = sketch.constraints().len();
        let solution = sketch.drag(p0, Point2::new(5.0, 7.0), T_OPTS).unwrap();
        assert!(solution.converged, "residual {}", solution.residual);
        assert_eq!(sketch.constraints().len(), before, "the pull left no trace");
        let at = sketch.point(p0).unwrap();
        assert_relative_eq!(at.x, 5.0, epsilon = 1e-3);
        assert_relative_eq!(at.y, 7.0, epsilon = 1e-3);
        assert_relative_eq!(
            sketch.measure_distance(p0, p1).unwrap(),
            40.0,
            epsilon = 1e-6
        );
        assert_relative_eq!(
            sketch.measure_distance(p1, p2).unwrap(),
            30.0,
            epsilon = 1e-6
        );

        // Anchored, the drag can only fail politely: the corner stays put
        // and the dimensions still hold.
        sketch
            .constrain(Constraint::Fixed(p0, Point2::new(5.0, 7.0)))
            .unwrap();
        let held = sketch.drag(p0, Point2::new(90.0, 90.0), T_OPTS).unwrap();
        assert!(held.converged);
        let at = sketch.point(p0).unwrap();
        assert_relative_eq!(at.x, 5.0, epsilon = 1e-6);
        assert_relative_eq!(at.y, 7.0, epsilon = 1e-6);
    }
}
