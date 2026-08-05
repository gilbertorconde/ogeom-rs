//! Sparse linear algebra: a compressed-sparse-row matrix and the conjugate
//! gradient least-squares solver built on it.
//!
//! The consumer this exists for is a constraint solver whose Jacobian has a
//! row per residual and entries only under the parameters that residual
//! names — a matrix that is nearly all zeros and whose dense factorization
//! cost grows with the whole parameter vector rather than with the
//! constraints. CGNR works entirely through products with the matrix and
//! its transpose, so its cost per iteration is the number of stored entries.
//!
//! Starting from zero, CGNR converges to the least-squares solution lying in
//! the row space — the *minimum-norm* solution — which is the same answer a
//! pseudo-inverse gives, so a caller can swap it for an SVD solve and expect
//! agreement, not merely feasibility.

/// A read-only sparse matrix in compressed-sparse-row form.
#[derive(Debug, Clone)]
pub struct SparseMatrix {
    rows: usize,
    cols: usize,
    row_starts: Vec<usize>,
    col_indices: Vec<usize>,
    values: Vec<f64>,
}

impl SparseMatrix {
    /// Build from `(row, column, value)` triplets. Duplicate positions sum,
    /// which is what an assembly loop wants; explicit zeros are kept out.
    ///
    /// # Panics
    /// If a triplet indexes outside `rows × cols`.
    #[must_use]
    pub fn from_triplets(rows: usize, cols: usize, triplets: &[(usize, usize, f64)]) -> Self {
        let mut sorted: Vec<(usize, usize, f64)> = triplets
            .iter()
            .inspect(|(r, c, _)| {
                assert!(*r < rows && *c < cols, "triplet outside the matrix");
            })
            .copied()
            .collect();
        sorted.sort_by_key(|&(r, c, _)| (r, c));

        let mut row_starts = vec![0usize; rows + 1];
        let mut col_indices = Vec::with_capacity(sorted.len());
        let mut values: Vec<f64> = Vec::with_capacity(sorted.len());
        let mut next = sorted.into_iter().peekable();
        for (r, start) in row_starts.iter_mut().enumerate().take(rows) {
            let row_begin = col_indices.len();
            *start = row_begin;
            while let Some(&(tr, c, v)) = next.peek() {
                if tr != r {
                    break;
                }
                next.next();
                if col_indices.len() > row_begin && col_indices.last() == Some(&c) {
                    let last = values.len() - 1;
                    values[last] += v;
                } else {
                    col_indices.push(c);
                    values.push(v);
                }
            }
        }
        row_starts[rows] = col_indices.len();
        Self {
            rows,
            cols,
            row_starts,
            col_indices,
            values,
        }
    }

    /// The shape, `(rows, cols)`.
    #[must_use]
    pub fn shape(&self) -> (usize, usize) {
        (self.rows, self.cols)
    }

    /// The number of stored entries.
    #[must_use]
    pub fn stored(&self) -> usize {
        self.values.len()
    }

    /// `y = A·x`.
    ///
    /// # Panics
    /// If `x` is not `cols` long.
    #[must_use]
    pub fn multiply(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.cols, "vector length must match columns");
        let mut y = vec![0.0; self.rows];
        for (y_r, window) in y.iter_mut().zip(self.row_starts.windows(2)) {
            let mut sum = 0.0;
            for i in window[0]..window[1] {
                sum = self.values[i].mul_add(x[self.col_indices[i]], sum);
            }
            *y_r = sum;
        }
        y
    }

    /// `y = Aᵀ·x`.
    ///
    /// # Panics
    /// If `x` is not `rows` long.
    #[must_use]
    pub fn transpose_multiply(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.rows, "vector length must match rows");
        let mut y = vec![0.0; self.cols];
        for (window, xr) in self.row_starts.windows(2).zip(x) {
            for i in window[0]..window[1] {
                y[self.col_indices[i]] = self.values[i].mul_add(*xr, y[self.col_indices[i]]);
            }
        }
        y
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).fold(0.0, |acc, (x, y)| x.mul_add(*y, acc))
}

/// Minimize `‖A·x − b‖` by conjugate gradient on the normal equations,
/// from zero, so a rank-deficient system yields the minimum-norm
/// least-squares solution. Convergence is declared when the gradient
/// `‖Aᵀ(b − A·x)‖` falls to `tolerance` relative to `‖Aᵀb‖`; `None` means
/// the iteration budget ran out before that happened.
///
/// # Panics
/// If `b` is not `rows` long.
#[must_use]
pub fn least_squares_cgnr(
    a: &SparseMatrix,
    b: &[f64],
    tolerance: f64,
    max_iterations: usize,
) -> Option<Vec<f64>> {
    let (rows, cols) = a.shape();
    assert_eq!(b.len(), rows, "right-hand side must match rows");
    let mut x = vec![0.0; cols];
    let mut r = b.to_vec();
    let mut z = a.transpose_multiply(&r);
    let target = tolerance * dot(&z, &z).sqrt().max(f64::MIN_POSITIVE);
    let mut p = z.clone();
    let mut zz = dot(&z, &z);
    if zz.sqrt() <= target {
        return Some(x);
    }
    for _ in 0..max_iterations {
        let w = a.multiply(&p);
        let ww = dot(&w, &w);
        if ww <= 0.0 {
            // A conjugate direction the matrix annihilates: the gradient is
            // as reduced as this arithmetic can make it.
            return Some(x);
        }
        let alpha = zz / ww;
        for (xi, pi) in x.iter_mut().zip(&p) {
            *xi = alpha.mul_add(*pi, *xi);
        }
        for (ri, wi) in r.iter_mut().zip(&w) {
            *ri = alpha.mul_add(-wi, *ri);
        }
        z = a.transpose_multiply(&r);
        let zz_next = dot(&z, &z);
        if zz_next.sqrt() <= target {
            return Some(x);
        }
        let beta = zz_next / zz;
        zz = zz_next;
        for (pi, zi) in p.iter_mut().zip(&z) {
            *pi = beta.mul_add(*pi, *zi);
        }
    }
    None
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn triplets_assemble_sum_and_multiply() {
        // [ 2 0 1 ]      duplicates at (0,0) sum to 2
        // [ 0 3 0 ]
        let a = SparseMatrix::from_triplets(
            2,
            3,
            &[(0, 0, 1.0), (0, 2, 1.0), (1, 1, 3.0), (0, 0, 1.0)],
        );
        assert_eq!(a.stored(), 3);
        assert_eq!(a.multiply(&[1.0, 2.0, 3.0]), vec![5.0, 6.0]);
        assert_eq!(a.transpose_multiply(&[1.0, 1.0]), vec![2.0, 3.0, 1.0]);
    }

    #[test]
    fn an_overdetermined_system_lands_on_the_normal_equation_answer() {
        // Fit y = c0 + c1 t through (0,1), (1,2), (2,4): the closed-form
        // least-squares line is c = (5/6, 3/2).
        let a = SparseMatrix::from_triplets(
            3,
            2,
            &[
                (0, 0, 1.0),
                (0, 1, 0.0),
                (1, 0, 1.0),
                (1, 1, 1.0),
                (2, 0, 1.0),
                (2, 1, 2.0),
            ],
        );
        let x = least_squares_cgnr(&a, &[1.0, 2.0, 4.0], 1e-14, 100).unwrap();
        assert!((x[0] - 5.0 / 6.0).abs() < 1e-10, "{x:?}");
        assert!((x[1] - 1.5).abs() < 1e-10, "{x:?}");
    }

    #[test]
    fn a_rank_deficient_system_returns_the_minimum_norm_solution() {
        // One equation, two unknowns: x + y = 2. The minimum-norm answer is
        // (1, 1) — the same point a pseudo-inverse names.
        let a = SparseMatrix::from_triplets(1, 2, &[(0, 0, 1.0), (0, 1, 1.0)]);
        let x = least_squares_cgnr(&a, &[2.0], 1e-14, 50).unwrap();
        assert!(
            (x[0] - 1.0).abs() < 1e-12 && (x[1] - 1.0).abs() < 1e-12,
            "{x:?}"
        );
    }

    #[test]
    fn a_zero_matrix_answers_zero_rather_than_spinning() {
        let a = SparseMatrix::from_triplets(2, 2, &[]);
        let x = least_squares_cgnr(&a, &[1.0, 1.0], 1e-12, 10).unwrap();
        assert_eq!(x, vec![0.0, 0.0]);
    }
}
