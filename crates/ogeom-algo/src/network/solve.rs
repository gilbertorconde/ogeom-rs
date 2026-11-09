//! The linear algebra of a network fit: the unknowns left free by linear
//! constraints on the poles, and a sparse symmetric positive definite
//! solve.

use std::collections::BTreeMap;

use ogeom_core::FastMap;

/// The unit roundoff.
const EPS: f64 = f64::EPSILON * 0.5;

/// How far below the largest diagonal entry a Cholesky pivot may fall
/// before the normal matrix counts as singular, as the single-patch fits
/// hold their own systems to.
const SINGULAR: f64 = 1e-12;

/// A coefficient and a bound on its rounding error.
#[derive(Debug, Clone, Copy)]
struct Term {
    value: f64,
    error: f64,
}

/// Each slot as a combination of the free unknowns.
#[derive(Debug, Clone)]
pub(super) struct Basis {
    /// The number of free unknowns.
    pub free: usize,
    /// Each slot's combination: `(unknown, coefficient)`.
    pub of_slot: Vec<Vec<(usize, f64)>>,
    /// Each slot's constant part, where some slots are prescribed: the
    /// slot is its combination plus this, coordinate by coordinate.
    pub offset: Option<Vec<[f64; 3]>>,
}

impl Basis {
    /// No constraint: each slot its own unknown.
    pub fn identity(slots: usize) -> Self {
        Self {
            free: slots,
            of_slot: (0..slots).map(|s| vec![(s, 1.0)]).collect(),
            offset: None,
        }
    }
}

/// An expression in the free slots plus a constant, each coordinate of the
/// constant with a bound on its rounding error.
type Affine = (BTreeMap<usize, Term>, [Term; 3]);

/// [`null_space`] where some slots are prescribed: each `(slot, value)` of
/// `known` holds that slot at that value, which the rows then carry as
/// constants, coordinate by coordinate.
///
/// The same elimination, each expression carrying its constant beside its
/// terms. A row left with no term but a constant beyond its rounding bound
/// asks the prescribed slots for what they do not hold: `None`.
pub(super) fn null_space_known(
    slots: usize,
    rows: &[Vec<(usize, f64)>],
    known: &[(usize, [f64; 3])],
) -> Option<Basis> {
    let zero = Term {
        value: 0.0,
        error: 0.0,
    };
    let mut solved: Vec<Option<Affine>> = vec![None; slots];
    for &(s, value) in known {
        let constant = value.map(|v| Term {
            value: v,
            error: 0.0,
        });
        solved[s] = Some((BTreeMap::new(), constant));
    }
    let mut order: Vec<usize> = Vec::new();
    let shift = |c: &mut [Term; 3], by: &[Term; 3], a: f64, a_error: f64| {
        for (to, from) in c.iter_mut().zip(by) {
            let value = a * from.value;
            to.value += value;
            to.error += a.abs() * from.error
                + a_error * from.value.abs()
                + EPS * (value.abs() + to.value.abs());
        }
    };
    for row in rows {
        let mut reduced: BTreeMap<usize, Term> = BTreeMap::new();
        let mut constant = [zero; 3];
        for &(s, a) in row {
            match &solved[s] {
                Some((expr, c)) => {
                    for (&t, b) in expr {
                        add(&mut reduced, t, a * b.value, a.abs() * b.error);
                    }
                    shift(&mut constant, c, a, 0.0);
                }
                None => add(&mut reduced, s, a, 0.0),
            }
        }
        reduced.retain(|_, t| t.value.abs() > t.error);
        let mut largest: Option<(usize, Term)> = None;
        for (&s, &t) in &reduced {
            if largest.is_none_or(|(_, b)| t.value.abs() > b.value.abs()) {
                largest = Some((s, t));
            }
        }
        let Some((pivot, head)) = largest else {
            if constant.iter().any(|c| c.value.abs() > c.error) {
                return None;
            }
            continue;
        };
        let mut expr: BTreeMap<usize, Term> = BTreeMap::new();
        for (&t, term) in &reduced {
            if t == pivot {
                continue;
            }
            let value = -term.value / head.value;
            let error =
                (term.error + value.abs() * head.error) / head.value.abs() + EPS * value.abs();
            expr.insert(t, Term { value, error });
        }
        let pivot_constant = constant.map(|c| {
            let value = -c.value / head.value;
            Term {
                value,
                error: (c.error + value.abs() * head.error) / head.value.abs() + EPS * value.abs(),
            }
        });
        for &s in &order {
            let Some((other, other_constant)) = solved[s].as_mut() else {
                continue;
            };
            let Some(at) = other.remove(&pivot) else {
                continue;
            };
            for (&t, term) in &expr {
                add(
                    other,
                    t,
                    at.value * term.value,
                    at.value.abs() * term.error + term.value.abs() * at.error,
                );
            }
            other.retain(|_, t| t.value.abs() > t.error);
            shift(other_constant, &pivot_constant, at.value, at.error);
        }
        solved[pivot] = Some((expr, pivot_constant));
        order.push(pivot);
    }
    let mut unknown = vec![usize::MAX; slots];
    let mut free = 0;
    for (s, done) in solved.iter().enumerate() {
        if done.is_none() {
            unknown[s] = free;
            free += 1;
        }
    }
    let of_slot = (0..slots)
        .map(|s| match &solved[s] {
            None => vec![(unknown[s], 1.0)],
            Some((expr, _)) => expr.iter().map(|(&t, c)| (unknown[t], c.value)).collect(),
        })
        .collect();
    let offset = (0..slots)
        .map(|s| match &solved[s] {
            None => [0.0; 3],
            Some((_, c)) => c.map(|t| t.value),
        })
        .collect();
    Some(Basis {
        free,
        of_slot,
        offset: Some(offset),
    })
}

/// The slots left free by homogeneous linear constraints, each row a list
/// of `(slot, coefficient)` summing to nought, and every slot as a
/// combination of them.
///
/// Gaussian elimination row by row, each row reduced by the slots already
/// solved for and solved for its largest coefficient. Each coefficient
/// carries a bound on its rounding error, and one no larger than its bound
/// is nought: a row with nothing left repeats rows before it (the rows
/// round a corner of four patches close on themselves), and adds nothing.
pub(super) fn null_space(slots: usize, rows: &[Vec<(usize, f64)>]) -> Basis {
    let mut solved: Vec<Option<BTreeMap<usize, Term>>> = vec![None; slots];
    let mut order: Vec<usize> = Vec::new();
    for row in rows {
        let mut reduced: BTreeMap<usize, Term> = BTreeMap::new();
        for &(s, a) in row {
            match &solved[s] {
                Some(expr) => {
                    for (&t, b) in expr {
                        add(&mut reduced, t, a * b.value, a.abs() * b.error);
                    }
                }
                None => add(&mut reduced, s, a, 0.0),
            }
        }
        reduced.retain(|_, t| t.value.abs() > t.error);
        let mut largest: Option<(usize, Term)> = None;
        for (&s, &t) in &reduced {
            if largest.is_none_or(|(_, b)| t.value.abs() > b.value.abs()) {
                largest = Some((s, t));
            }
        }
        let Some((pivot, head)) = largest else {
            continue;
        };
        let mut expr: BTreeMap<usize, Term> = BTreeMap::new();
        for (&t, term) in &reduced {
            if t == pivot {
                continue;
            }
            let value = -term.value / head.value;
            let error =
                (term.error + value.abs() * head.error) / head.value.abs() + EPS * value.abs();
            expr.insert(t, Term { value, error });
        }
        for &s in &order {
            let Some(other) = solved[s].as_mut() else {
                continue;
            };
            let Some(at) = other.remove(&pivot) else {
                continue;
            };
            for (&t, term) in &expr {
                add(
                    other,
                    t,
                    at.value * term.value,
                    at.value.abs() * term.error + term.value.abs() * at.error,
                );
            }
            other.retain(|_, t| t.value.abs() > t.error);
        }
        solved[pivot] = Some(expr);
        order.push(pivot);
    }
    let mut unknown = vec![usize::MAX; slots];
    let mut free = 0;
    for (s, done) in solved.iter().enumerate() {
        if done.is_none() {
            unknown[s] = free;
            free += 1;
        }
    }
    let of_slot = (0..slots)
        .map(|s| match &solved[s] {
            None => vec![(unknown[s], 1.0)],
            Some(expr) => expr.iter().map(|(&t, c)| (unknown[t], c.value)).collect(),
        })
        .collect();
    Basis {
        free,
        of_slot,
        offset: None,
    }
}

/// Add `value`, itself in error by up to `error`, to the entry at `slot`.
fn add(into: &mut BTreeMap<usize, Term>, slot: usize, value: f64, error: f64) {
    let entry = into.entry(slot).or_insert(Term {
        value: 0.0,
        error: 0.0,
    });
    entry.value += value;
    entry.error += error + EPS * (value.abs() + entry.value.abs());
}

/// A symmetric matrix by its lower triangle, entries summed as they come.
#[derive(Debug, Default)]
pub(super) struct Lower {
    entries: FastMap<(usize, usize), f64>,
}

impl Lower {
    /// Add `value` at `(row, column)` of the full matrix; only the lower
    /// triangle is kept, the upper being its mirror.
    pub fn add(&mut self, row: usize, column: usize, value: f64) {
        if row >= column {
            *self.entries.entry((row, column)).or_insert(0.0) += value;
        }
    }
}

/// Solve `A x = b` for each right-hand side, `A` symmetric positive
/// definite of order `n`, by Cholesky over the envelope of `A` reordered
/// by reverse Cuthill-McKee. `None` where a pivot falls to the singular
/// floor below the largest diagonal entry.
pub(super) fn solve_spd(n: usize, a: &Lower, rhs: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(r, c) in a.entries.keys() {
        if r != c {
            neighbours[r].push(c);
            neighbours[c].push(r);
        }
    }
    for list in &mut neighbours {
        list.sort_unstable();
        list.dedup();
    }
    let order = reverse_cuthill_mckee(&neighbours);
    let mut position = vec![0; n];
    for (k, &i) in order.iter().enumerate() {
        position[i] = k;
    }
    // Each permuted row's first column inside the envelope.
    let mut first: Vec<usize> = (0..n).collect();
    for &(r, c) in a.entries.keys() {
        let (pr, pc) = (position[r], position[c]);
        let (hi, lo) = if pr >= pc { (pr, pc) } else { (pc, pr) };
        first[hi] = first[hi].min(lo);
    }
    let mut rows: Vec<Vec<f64>> = (0..n).map(|i| vec![0.0; i - first[i] + 1]).collect();
    // Summed in a fixed order, so the factor does not hang on the map's.
    let mut placed: Vec<(usize, usize, f64)> = a
        .entries
        .iter()
        .map(|(&(r, c), &v)| {
            let (pr, pc) = (position[r], position[c]);
            if pr >= pc { (pr, pc, v) } else { (pc, pr, v) }
        })
        .collect();
    placed.sort_by_key(|x| (x.0, x.1));
    for (r, c, v) in placed {
        rows[r][c - first[r]] += v;
    }
    let largest = (0..n).fold(0.0_f64, |m, i| m.max(rows[i][i - first[i]].abs()));
    let floor = largest * SINGULAR;
    for i in 0..n {
        for j in first[i]..=i {
            let start = first[i].max(first[j]);
            let mut sum = rows[i][j - first[i]];
            for k in start..j {
                sum -= rows[i][k - first[i]] * rows[j][k - first[j]];
            }
            if i == j {
                if !(sum > floor && sum.is_finite()) {
                    return None;
                }
                rows[i][i - first[i]] = sum.sqrt();
            } else {
                rows[i][j - first[i]] = sum / rows[j][j - first[j]];
            }
        }
    }
    let solved = rhs
        .iter()
        .map(|b| {
            let mut y: Vec<f64> = order.iter().map(|&i| b[i]).collect();
            for i in 0..n {
                let mut sum = y[i];
                for k in first[i]..i {
                    sum -= rows[i][k - first[i]] * y[k];
                }
                y[i] = sum / rows[i][i - first[i]];
            }
            for i in (0..n).rev() {
                y[i] /= rows[i][i - first[i]];
                let yi = y[i];
                for k in first[i]..i {
                    y[k] -= rows[i][k - first[i]] * yi;
                }
            }
            let mut x = vec![0.0; n];
            for (k, &i) in order.iter().enumerate() {
                x[i] = y[k];
            }
            x
        })
        .collect();
    Some(solved)
}

/// An order of the nodes that keeps neighbours close: breadth first from
/// a node of least degree in each component, neighbours by rising degree,
/// the whole reversed. Ties go to the lower index.
fn reverse_cuthill_mckee(neighbours: &[Vec<usize>]) -> Vec<usize> {
    let n = neighbours.len();
    let mut seen = vec![false; n];
    let mut order = Vec::with_capacity(n);
    let mut by_degree: Vec<usize> = (0..n).collect();
    by_degree.sort_by_key(|&i| (neighbours[i].len(), i));
    for &start in &by_degree {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        let mut head = order.len();
        order.push(start);
        while head < order.len() {
            let at = order[head];
            head += 1;
            let mut next: Vec<usize> = neighbours[at]
                .iter()
                .copied()
                .filter(|&m| !seen[m])
                .collect();
            next.sort_by_key(|&m| (neighbours[m].len(), m));
            for m in next {
                seen[m] = true;
                order.push(m);
            }
        }
    }
    order.reverse();
    order
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn a_cycle_of_constraints_leaves_its_repeat_out() {
        // x0 = x1, x1 = x2, x2 = x0: the third repeats the first two.
        let rows = vec![
            vec![(0, 1.0), (1, -1.0)],
            vec![(1, 1.0), (2, -1.0)],
            vec![(2, 1.0), (0, -1.0)],
        ];
        let basis = null_space(3, &rows);
        assert_eq!(basis.free, 1);
        for combo in &basis.of_slot {
            assert_eq!(combo.len(), 1);
            assert!((combo[0].1 - 1.0).abs() < 1e-15);
        }
    }

    #[test]
    fn the_solve_matches_a_dense_one() {
        // A path graph's Laplacian plus the identity, against its known
        // product.
        let n = 7;
        let mut a = Lower::default();
        for i in 0..n {
            a.add(i, i, 3.0);
            if i + 1 < n {
                a.add(i + 1, i, -1.0);
            }
        }
        let x: Vec<f64> = (0..n)
            .map(|i| f64::from(u32::try_from(i).unwrap()).sin())
            .collect();
        let b: Vec<f64> = (0..n)
            .map(|i| {
                let mut s = 3.0 * x[i];
                if i > 0 {
                    s -= x[i - 1];
                }
                if i + 1 < n {
                    s -= x[i + 1];
                }
                s
            })
            .collect();
        let got = solve_spd(n, &a, &[b]).unwrap();
        for (g, w) in got[0].iter().zip(&x) {
            assert!((g - w).abs() < 1e-14, "{g} against {w}");
        }
    }
}
