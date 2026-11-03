//! A linear program in a few free variables: minimise `c·z` subject to
//! rows `a·z ≤ b`.
//!
//! Solved as its dual, `min b·y` subject to `Σ y_r a_r = −c` and `y ≥ 0`:
//! a tableau with one row per variable however many constraint rows there
//! are, so a pivot costs the variables times the rows. Phase one drives
//! out one artificial column per variable; phase two minimises. Pivots
//! follow Bland's rule, the least column that improves and the least basic
//! row among tied ratios, under which the simplex method cannot cycle. The
//! variables are the dual's simplex multipliers, read off the artificial
//! columns' reduced costs; the rows a dual basis holds are the constraints
//! the optimum stands on.

/// A row's entries are within a few hundred roundings of the largest
/// entry or cost they are compared with; below that, an entry is zero.
pub(super) const ROUNDING: f64 = 256.0 * f64::EPSILON;

/// What a linear program came to.
#[derive(Debug, Clone, PartialEq)]
pub(super) enum Solved {
    /// The optimum `z`.
    Optimal(Vec<f64>),
    /// The objective falls without bound over the rows.
    Unbounded,
    /// The tableau could not be brought to an optimum: the rows hold no
    /// point at all, or the pivots ran past their bound.
    Failed,
}

/// Minimise `c·z` subject to `a·z ≤ b` for each row `(a, b)`, every `a` as
/// long as `c`.
pub(super) fn minimise(c: &[f64], rows: &[(&[f64], f64)]) -> Solved {
    let n = c.len();
    let m = rows.len();
    let width = m + n;
    // The dual: row i reads Σ_r a_r[i] y_r = −c[i], made non-negative on
    // the right by its sign, with artificial column m + i.
    let sign: Vec<f64> = c
        .iter()
        .map(|ci| if -ci >= 0.0 { 1.0 } else { -1.0 })
        .collect();
    let mut table = vec![0.0; n * width];
    let mut rhs = vec![0.0; n];
    for i in 0..n {
        for (r, (a, _)) in rows.iter().enumerate() {
            table[i * width + r] = sign[i] * a[i];
        }
        table[i * width + m + i] = 1.0;
        rhs[i] = -sign[i] * c[i];
    }
    let mut basis: Vec<usize> = (m..width).collect();

    let entry = rows
        .iter()
        .flat_map(|(a, _)| a.iter())
        .chain(c)
        .fold(1.0_f64, |big, v| big.max(v.abs()))
        * ROUNDING;
    let cost_scale = rows.iter().fold(0.0_f64, |big, (_, b)| big.max(b.abs()));

    // Phase one: the artificials' sum made least.
    let artificial = |j: usize| if j >= m { 1.0 } else { 0.0 };
    let phase_one = pivots(
        &mut table,
        &mut rhs,
        &mut basis,
        width,
        &artificial,
        width,
        entry,
        entry,
    );
    if phase_one != Pivoted::Optimal {
        return Solved::Failed;
    }
    let left: f64 = (0..n).filter(|&i| basis[i] >= m).map(|i| rhs[i]).sum();
    if left > entry * c.iter().fold(1.0_f64, |big, v| big.max(v.abs())) {
        // No multipliers make the rows' combination −c: some direction
        // lowers c·z and no row stops it.
        return Solved::Unbounded;
    }
    // An artificial still basic at zero is pivoted out on any structural
    // entry its row has; a row with none is a combination of the others.
    for i in 0..n {
        if basis[i] < m {
            continue;
        }
        let pick = (0..m)
            .filter(|&j| table[i * width + j].abs() > entry)
            .max_by(|&a, &b| {
                table[i * width + a]
                    .abs()
                    .total_cmp(&table[i * width + b].abs())
            });
        if let Some(j) = pick {
            pivot(&mut table, &mut rhs, &mut basis, width, i, j);
        }
    }

    // Phase two: the constraint rows' bounds, over the structural columns.
    let cost = |j: usize| if j < m { rows[j].1 } else { 0.0 };
    let tolerance = cost_scale.max(f64::MIN_POSITIVE) * ROUNDING;
    match pivots(
        &mut table, &mut rhs, &mut basis, width, &cost, m, entry, tolerance,
    ) {
        Pivoted::Optimal => {}
        Pivoted::Unbounded => return Solved::Failed,
        Pivoted::Stalled => return Solved::Failed,
    }

    // z_i = σ_i π_i, and an artificial's reduced cost is −π_i.
    let z = (0..n)
        .map(|i| {
            let column = m + i;
            let reduced = cost(column)
                - (0..n)
                    .map(|k| cost(basis[k]) * table[k * width + column])
                    .sum::<f64>();
            -sign[i] * reduced
        })
        .collect();
    Solved::Optimal(z)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Pivoted {
    Optimal,
    Unbounded,
    Stalled,
}

/// Bland's-rule pivots until no column below `enter_below` improves the
/// objective `cost` by more than `tolerance`.
#[allow(clippy::too_many_arguments, reason = "one tableau, passed apart")]
fn pivots(
    table: &mut [f64],
    rhs: &mut [f64],
    basis: &mut [usize],
    width: usize,
    cost: &dyn Fn(usize) -> f64,
    enter_below: usize,
    entry: f64,
    tolerance: f64,
) -> Pivoted {
    let n = basis.len();
    // Bland's rule ends in at most the number of bases, but each pivot
    // here is cheap and a tableau this small settles in tens of pivots;
    // the bound only stops a tableau that rounding has set cycling.
    let most = 64 * (width + n) * n;
    for _ in 0..most {
        let entering = (0..enter_below).find(|&j| {
            if basis.contains(&j) {
                return false;
            }
            let reduced = cost(j)
                - (0..n)
                    .map(|k| cost(basis[k]) * table[k * width + j])
                    .sum::<f64>();
            reduced < -tolerance
        });
        let Some(j) = entering else {
            return Pivoted::Optimal;
        };
        let mut leave: Option<(usize, f64)> = None;
        for i in 0..n {
            let a = table[i * width + j];
            if a <= entry {
                continue;
            }
            let ratio = rhs[i].max(0.0) / a;
            leave = match leave {
                Some((k, best)) if ratio > best || (ratio == best && basis[k] < basis[i]) => {
                    Some((k, best))
                }
                _ => Some((i, ratio)),
            };
        }
        let Some((i, _)) = leave else {
            return Pivoted::Unbounded;
        };
        pivot(table, rhs, basis, width, i, j);
    }
    Pivoted::Stalled
}

fn pivot(
    table: &mut [f64],
    rhs: &mut [f64],
    basis: &mut [usize],
    width: usize,
    i: usize,
    j: usize,
) {
    let n = basis.len();
    let p = table[i * width + j];
    for v in &mut table[i * width..(i + 1) * width] {
        *v /= p;
    }
    rhs[i] /= p;
    for k in 0..n {
        if k == i {
            continue;
        }
        let f = table[k * width + j];
        if f == 0.0 {
            continue;
        }
        for col in 0..width {
            table[k * width + col] -= f * table[i * width + col];
        }
        rhs[k] -= f * rhs[i];
    }
    basis[i] = j;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn solve(c: &[f64], rows: &[(Vec<f64>, f64)]) -> Solved {
        let borrowed: Vec<(&[f64], f64)> = rows.iter().map(|(a, b)| (a.as_slice(), *b)).collect();
        minimise(c, &borrowed)
    }

    fn optimum(solved: Solved) -> Vec<f64> {
        match solved {
            Solved::Optimal(z) => z,
            other => panic!("expected an optimum, got {other:?}"),
        }
    }

    /// A box's corner: min −x − 2y over 0 ≤ x ≤ 3, 0 ≤ y ≤ 2, x + y ≤ 4
    /// stands at (2, 2).
    #[test]
    fn a_polygon_gives_its_vertex() {
        let rows = vec![
            (vec![1.0, 0.0], 3.0),
            (vec![-1.0, 0.0], 0.0),
            (vec![0.0, 1.0], 2.0),
            (vec![0.0, -1.0], 0.0),
            (vec![1.0, 1.0], 4.0),
        ];
        let z = optimum(solve(&[-1.0, -2.0], &rows));
        assert!(
            (z[0] - 2.0).abs() < 1e-12 && (z[1] - 2.0).abs() < 1e-12,
            "{z:?}"
        );
    }

    /// The Chebyshev line through (0, 0), (1, 1), (2, 0): y = 0.5 with a
    /// worst error of 0.5, every point touching.
    #[test]
    fn a_chebyshev_line_alternates() {
        // Variables (a, b, t): |a + b x − y| ≤ t, minimise t.
        let data = [(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)];
        let mut rows = Vec::new();
        for (x, y) in data {
            rows.push((vec![1.0, x, -1.0], y));
            rows.push((vec![-1.0, -x, -1.0], -y));
        }
        let z = optimum(solve(&[0.0, 0.0, 1.0], &rows));
        assert!((z[0] - 0.5).abs() < 1e-12, "{z:?}");
        assert!(z[1].abs() < 1e-12, "{z:?}");
        assert!((z[2] - 0.5).abs() < 1e-12, "{z:?}");
    }

    /// x has only a lower bound, and minimising −x runs off.
    #[test]
    fn an_open_direction_is_unbounded() {
        let rows = vec![
            (vec![-1.0, 0.0], 0.0),
            (vec![0.0, 1.0], 1.0),
            (vec![0.0, -1.0], 1.0),
        ];
        assert_eq!(solve(&[-1.0, 0.0], &rows), Solved::Unbounded);
    }

    /// Degenerate rows, many through one vertex, settle without cycling.
    #[test]
    fn many_rows_through_one_vertex_settle() {
        let mut rows = Vec::new();
        for k in 0..40 {
            let a = f64::from(k) * 0.05;
            rows.push((vec![a.cos(), a.sin()], 0.0));
        }
        let z = optimum(solve(&[-1.0, -1.0], &rows));
        assert!(z[0].abs() < 1e-12 && z[1].abs() < 1e-12, "{z:?}");
        for (a, b) in &rows {
            assert!(a[0] * z[0] + a[1] * z[1] <= b + 1e-12, "{z:?} breaks {a:?}");
        }
    }
}
