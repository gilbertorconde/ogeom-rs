//! The height patch an N-sided filling is fitted as.
//!
//! The patch stands over a plane: `S(u, v) = O + u·e1 + v·e2 + h(u, v)·n`,
//! with `(u, v)` the plane coordinates and `h` a scalar cubic B-spline over
//! a rectangle covering the hole. With the frame fixed, every condition a
//! filling asks for is linear in the control heights: a point the surface
//! passes through fixes `h`, a tangent plane fixes `h_u` and `h_v`, and a
//! second fundamental form (once the tangent plane is fixed) fixes the
//! three second derivatives. The heights minimise the conditions' squared
//! residuals plus a weighted thin-plate bending energy over the whole
//! rectangle, which settles every control the conditions leave free.
//!
//! The surface itself is an ordinary B-spline patch: the plane part is
//! carried by control points at the Greville abscissae, which reproduce a
//! linear function exactly, so the patch's chart is the plane's.
//!
//! A loop whose corner is seen smooth along the plane's normal is no
//! height field over it: the two sides meeting there ask for two slopes at
//! one point. Nor is a hole whose sides' supports stand square to the
//! plane, as a tube's wall does to its rim's. The free patch fits all three
//! coordinates of its control points over a chart the caller draws from
//! the loop itself, each condition fixing a derivative of the point
//! `S(u, v)`: a point condition fixes `S`, and a tangent ribbon fixes the
//! derivative of `S` across the side, a vector that may run square to the
//! plane. Every condition is still linear, and the three coordinates share
//! one matrix.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::BSplineSurface;
use ogeom_math::{ControlGrid, KnotVector, Point, Point2, Vector, Vector2};

/// The degree of the height patch in both directions: cubic, so the
/// curvature a G2 side asks for is continuous across the patch.
pub(crate) const DEGREE: usize = 3;

/// The plane the patch stands over: an origin, two in-plane unit axes and
/// the unit normal `e1 × e2`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct PlaneFrame {
    pub origin: Point,
    pub e1: Vector,
    pub e2: Vector,
    pub n: Vector,
}

impl PlaneFrame {
    /// A point's plane coordinates.
    pub fn chart(&self, p: Point) -> Point2 {
        let d = p - self.origin;
        Point2::new(d.dot(self.e1), d.dot(self.e2))
    }

    /// A point's height above the plane.
    pub fn height(&self, p: Point) -> f64 {
        (p - self.origin).dot(self.n)
    }
}

/// One linear condition on the patch's `D` fitted components (the height,
/// or a point's three coordinates): a derivative of them at `at` equals
/// `target`, the row scaled by `weight` before squaring. The derivative is
/// a sum of partial derivatives, each `(order_u, order_v, coefficient)`,
/// orders at most two.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Condition<const D: usize> {
    pub at: Point2,
    pub terms: [(usize, usize, f64); 3],
    pub target: [f64; D],
    pub weight: f64,
}

impl<const D: usize> Condition<D> {
    /// The partial derivative of order `(order_u, order_v)`.
    pub fn partial(at: Point2, order: (usize, usize), target: [f64; D], weight: f64) -> Self {
        Self {
            at,
            terms: [(order.0, order.1, 1.0), (0, 0, 0.0), (0, 0, 0.0)],
            target,
            weight,
        }
    }

    /// The first (`second` false) or second derivative along the unit
    /// chart direction `d`.
    pub fn along(at: Point2, d: Vector2, second: bool, target: [f64; D], weight: f64) -> Self {
        let terms = if second {
            [
                (2, 0, d.x * d.x),
                (1, 1, 2.0 * d.x * d.y),
                (0, 2, d.y * d.y),
            ]
        } else {
            [(1, 0, d.x), (0, 1, d.y), (0, 0, 0.0)]
        };
        Self {
            at,
            terms,
            target,
            weight,
        }
    }
}

/// A symmetric positive definite matrix held as its lower band.
struct Banded {
    size: usize,
    band: usize,
    /// `rows[k][d]` is entry `(k, k - d)`.
    rows: Vec<Vec<f64>>,
}

impl Banded {
    fn new(size: usize, band: usize) -> Self {
        Self {
            size,
            band,
            rows: vec![vec![0.0; band + 1]; size],
        }
    }

    fn add(&mut self, i: usize, j: usize, value: f64) {
        let (hi, lo) = if i >= j { (i, j) } else { (j, i) };
        let d = hi - lo;
        debug_assert!(d <= self.band, "an entry outside the band");
        if d <= self.band {
            self.rows[hi][d] += value;
        }
    }

    fn get(&self, i: usize, j: usize) -> f64 {
        let (hi, lo) = if i >= j { (i, j) } else { (j, i) };
        let d = hi - lo;
        if d > self.band { 0.0 } else { self.rows[hi][d] }
    }

    /// The banded Cholesky factor, `factor[k][d]` being `L(k, k - d)`;
    /// `None` where the matrix is not positive definite to working
    /// precision.
    fn factor(&self) -> Option<Vec<Vec<f64>>> {
        let n = self.size;
        let b = self.band;
        let mut factor = vec![vec![0.0f64; b + 1]; n];
        for k in 0..n {
            let first = k.saturating_sub(b);
            for j in first..=k {
                let mut sum = self.get(k, j);
                let from = first.max(j.saturating_sub(b));
                for m in from..j {
                    sum -= factor[k][k - m] * factor[j][j - m];
                }
                if j == k {
                    if sum <= 0.0 || !sum.is_finite() {
                        return None;
                    }
                    factor[k][0] = sum.sqrt();
                } else {
                    factor[k][k - j] = sum / factor[j][0];
                }
            }
        }
        Some(factor)
    }

    /// Solve with the factor [`Banded::factor`] gave.
    fn solve(&self, factor: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
        let n = self.size;
        let b = self.band;
        let mut y = vec![0.0f64; n];
        for k in 0..n {
            let mut sum = rhs[k];
            for m in k.saturating_sub(b)..k {
                sum -= factor[k][k - m] * y[m];
            }
            y[k] = sum / factor[k][0];
        }
        let mut x = vec![0.0f64; n];
        for k in (0..n).rev() {
            let mut sum = y[k];
            for m in (k + 1)..n.min(k + b + 1) {
                sum -= factor[m][m - k] * x[m];
            }
            x[k] = sum / factor[k][0];
        }
        x
    }
}

/// The non-zero basis functions and their first two derivatives at `t`:
/// the first index, and the values by order.
fn basis_of(knots: &KnotVector, t: f64) -> (usize, [Vec<f64>; 3]) {
    let span = knots.span_unchecked(t);
    let rows = knots.basis_derivatives(span, t, 2);
    (
        span - DEGREE,
        [rows[0].to_vec(), rows[1].to_vec(), rows[2].to_vec()],
    )
}

/// `∫ N_i^(r) N_j^(r)` over the knot domain, for `r` = 0, 1, 2: the Gram
/// matrices the bending energy is built from, exact by four-point Gauss on
/// every knot span (the products are polynomials of degree at most six).
fn gram(knots: &KnotVector, count: usize) -> [Vec<Vec<f64>>; 3] {
    const NODES: [(f64, f64); 4] = [
        (-0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
        (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
    ];
    let mut out = [
        vec![vec![0.0; count]; count],
        vec![vec![0.0; count]; count],
        vec![vec![0.0; count]; count],
    ];
    let distinct: Vec<f64> = knots.distinct().iter().map(|(k, _)| *k).collect();
    for w in distinct.windows(2) {
        let (a, b) = (w[0], w[1]);
        let half = 0.5 * (b - a);
        let mid = 0.5 * (a + b);
        for (x, weight) in NODES {
            let t = half.mul_add(x, mid);
            let span = knots.span_unchecked(t);
            let rows = knots.basis_derivatives(span, t, 2);
            let first = span - DEGREE;
            for (r, gram_r) in out.iter_mut().enumerate() {
                for (i, bi) in rows[r].iter().enumerate() {
                    for (j, bj) in rows[r].iter().enumerate() {
                        gram_r[first + i][first + j] += weight * half * bi * bj;
                    }
                }
            }
        }
    }
    out
}

/// A fitted patch and how far it misses its conditions.
pub(crate) struct PatchFit {
    pub surface: BSplineSurface,
    /// The conditions' weighted squared residuals, summed.
    pub residual: f64,
}

/// The controls of one fit: the knots, the `D` components of every
/// control (row-major, `v` fastest) and the residual.
struct Net<const D: usize> {
    u_knots: KnotVector,
    v_knots: KnotVector,
    values: Vec<[f64; D]>,
    residual: f64,
}

/// Fit `D` components over `domain` with `controls` controls per
/// direction to `conditions`, with the thin-plate energy at `smoothing`
/// on each component. The components share one matrix and are solved
/// apart.
fn solve_net<const D: usize>(
    domain: ((f64, f64), (f64, f64)),
    controls: (usize, usize),
    conditions: &[Condition<D>],
    smoothing: f64,
) -> OgeomResult<Net<D>> {
    let (nu, nv) = controls;
    let ((ua, ub), (va, vb)) = domain;
    let u_knots = KnotVector::clamped_uniform(DEGREE, nu)?.reparameterized(ua, ub)?;
    let v_knots = KnotVector::clamped_uniform(DEGREE, nv)?.reparameterized(va, vb)?;
    let size = nu * nv;
    let band = DEGREE * nv + DEGREE;
    let mut matrix = Banded::new(size, band);
    let mut rhs = vec![[0.0f64; D]; size];

    // The bending energy ∫∫ h_uu² + 2 h_uv² + h_vv².
    let [gu0, gu1, gu2] = gram(&u_knots, nu);
    let [gv0, gv1, gv2] = gram(&v_knots, nv);
    for i in 0..nu {
        for i2 in i.saturating_sub(DEGREE)..(i + DEGREE + 1).min(nu) {
            for j in 0..nv {
                for j2 in j.saturating_sub(DEGREE)..(j + DEGREE + 1).min(nv) {
                    let (k, k2) = (i * nv + j, i2 * nv + j2);
                    if k2 > k {
                        continue;
                    }
                    let e = gu2[i][i2] * gv0[j][j2]
                        + 2.0 * gu1[i][i2] * gv1[j][j2]
                        + gu0[i][i2] * gv2[j][j2];
                    matrix.add(k, k2, smoothing * e);
                }
            }
        }
    }

    // The conditions, as least-squares rows.
    let mut rows = Vec::with_capacity(conditions.len());
    for c in conditions {
        let (fu, bu) = basis_of(&u_knots, c.at.x.clamp(ua, ub));
        let (fv, bv) = basis_of(&v_knots, c.at.y.clamp(va, vb));
        let mut local = [[0.0f64; DEGREE + 1]; DEGREE + 1];
        for &(ou, ov, coefficient) in &c.terms {
            for (line, a) in local.iter_mut().zip(&bu[ou]) {
                for (slot, b) in line.iter_mut().zip(&bv[ov]) {
                    *slot += coefficient * a * b;
                }
            }
        }
        let mut row: Vec<(usize, f64)> = Vec::with_capacity((DEGREE + 1) * (DEGREE + 1));
        for (i, line) in local.iter().enumerate() {
            for (j, value) in line.iter().enumerate() {
                row.push(((fu + i) * nv + (fv + j), c.weight * value));
            }
        }
        let target = c.target.map(|t| c.weight * t);
        for &(k, a) in &row {
            for (slot, t) in rhs[k].iter_mut().zip(target) {
                *slot += a * t;
            }
            for &(k2, b) in &row {
                if k2 <= k {
                    matrix.add(k, k2, a * b);
                }
            }
        }
        rows.push((row, target));
    }

    let Some(factor) = matrix.factor() else {
        ogeom_bail!(
            Numeric,
            "the filling's system is singular at {nu}x{nv} controls"
        );
    };
    let mut values = vec![[0.0f64; D]; size];
    for d in 0..D {
        let column: Vec<f64> = rhs.iter().map(|r| r[d]).collect();
        for (slot, x) in values.iter_mut().zip(matrix.solve(&factor, &column)) {
            slot[d] = x;
        }
    }
    let mut residual = 0.0f64;
    for (row, target) in &rows {
        for (d, t) in target.iter().enumerate() {
            let value: f64 = row.iter().map(|&(k, a)| a * values[k][d]).sum();
            residual += (value - t).powi(2);
        }
    }
    Ok(Net {
        u_knots,
        v_knots,
        values,
        residual,
    })
}

/// Fit the height over `domain` with `controls` control heights per
/// direction to `conditions`, with the thin-plate energy at `smoothing`.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// the control counts cannot carry the degree, and
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if the system
/// is singular.
pub(crate) fn fit_height(
    frame: &PlaneFrame,
    domain: ((f64, f64), (f64, f64)),
    controls: (usize, usize),
    conditions: &[Condition<1>],
    smoothing: f64,
    tol: Tolerances,
) -> OgeomResult<PatchFit> {
    let (nu, nv) = controls;
    let net = solve_net(domain, controls, conditions, smoothing)?;
    // Greville abscissae carry the plane part exactly.
    let greville = |knots: &KnotVector, i: usize| -> f64 {
        let k = knots.knots();
        let window = &k[i + 1..=i + DEGREE];
        #[expect(clippy::cast_precision_loss, reason = "the degree, three")]
        let count = window.len() as f64;
        window.iter().sum::<f64>() / count
    };
    let mut control = Vec::with_capacity(nu * nv);
    for i in 0..nu {
        let u = greville(&net.u_knots, i);
        for j in 0..nv {
            let v = greville(&net.v_knots, j);
            let h = net.values[i * nv + j][0];
            control.push(frame.origin + frame.e1 * u + frame.e2 * v + frame.n * h);
        }
    }
    let grid = ControlGrid::new(control, nu, nv)?;
    Ok(PatchFit {
        surface: BSplineSurface::new(net.u_knots, net.v_knots, &grid, tol)?,
        residual: net.residual,
    })
}

/// Fit a patch whose control points' three coordinates are all free over
/// `domain`, with `controls` controls per direction, to `conditions` on
/// the point `S(u, v)`, with the thin-plate energy at `smoothing` on each
/// coordinate. Where the height patch keeps `(u, v)` as the plane's
/// coordinates, this one leaves the chart to the conditions: `(u, v)` is
/// whatever chart they are placed on.
///
/// # Errors
///
/// As [`fit_height`].
pub(crate) fn fit_free(
    domain: ((f64, f64), (f64, f64)),
    controls: (usize, usize),
    conditions: &[Condition<3>],
    smoothing: f64,
    tol: Tolerances,
) -> OgeomResult<PatchFit> {
    let (nu, nv) = controls;
    let net = solve_net(domain, controls, conditions, smoothing)?;
    let control: Vec<Point> = net
        .values
        .iter()
        .map(|[x, y, z]| Point::new(*x, *y, *z))
        .collect();
    let grid = ControlGrid::new(control, nu, nv)?;
    Ok(PatchFit {
        surface: BSplineSurface::new(net.u_knots, net.v_knots, &grid, tol)?,
        residual: net.residual,
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code")]
mod tests {
    use super::*;
    use ogeom_geom::Surface as _;

    const T: Tolerances = Tolerances::millimetres();

    fn flat() -> PlaneFrame {
        PlaneFrame {
            origin: Point::ORIGIN,
            e1: Vector::new(1.0, 0.0, 0.0),
            e2: Vector::new(0.0, 1.0, 0.0),
            n: Vector::new(0.0, 0.0, 1.0),
        }
    }

    #[test]
    fn a_quadratic_height_fixed_on_a_dense_grid_is_reproduced() {
        // h = x² - x·y + 0.5 y² lies in the cubic spline space, so values on
        // a grid dense enough to pin every control give it back to rounding.
        let h = |x: f64, y: f64| x * x - x * y + 0.5 * y * y;
        let mut conditions = Vec::new();
        for i in 0..=20 {
            for j in 0..=20 {
                let (x, y) = (f64::from(i) / 10.0 - 1.0, f64::from(j) / 10.0 - 1.0);
                conditions.push(Condition::partial(
                    Point2::new(x, y),
                    (0, 0),
                    [h(x, y)],
                    1.0,
                ));
            }
        }
        let patch = fit_height(
            &flat(),
            ((-1.0, 1.0), (-1.0, 1.0)),
            (7, 6),
            &conditions,
            1e-12,
            T,
        )
        .unwrap()
        .surface;
        for (x, y) in [(0.13, -0.71), (-0.9, 0.4), (0.55, 0.55)] {
            let p = patch.point_at(x, y, T).unwrap();
            assert!((p.x - x).abs() < 1e-12 && (p.y - y).abs() < 1e-12, "{p:?}");
            assert!((p.z - h(x, y)).abs() < 1e-9, "{p:?} against {}", h(x, y));
        }
    }

    #[test]
    fn slope_and_bend_conditions_set_the_derivatives_they_name() {
        // A plane's own value at one point, slopes everywhere along a line,
        // and a bend: the fit honours the derivative orders, read back from
        // the surface's own derivatives.
        let mut conditions = vec![Condition::partial(
            Point2::new(0.0, 0.0),
            (0, 0),
            [0.25],
            1.0,
        )];
        for k in 0..=10 {
            let y = f64::from(k) / 5.0 - 1.0;
            for (order, target) in [((1, 0), 0.5), ((0, 1), -0.25), ((2, 0), 0.0)] {
                conditions.push(Condition::partial(
                    Point2::new(0.0, y),
                    order,
                    [target],
                    1.0,
                ));
            }
        }
        let patch = fit_height(
            &flat(),
            ((-1.0, 1.0), (-1.0, 1.0)),
            (6, 6),
            &conditions,
            1e-9,
            T,
        )
        .unwrap()
        .surface;
        let (du, dv) = patch.d1_at(0.0, 0.3, T).unwrap();
        assert!((du.z - 0.5).abs() < 1e-6, "{du:?}");
        assert!((dv.z + 0.25).abs() < 1e-6, "{dv:?}");
        let p = patch.point_at(0.7, -0.4, T).unwrap();
        let plane = 0.25 + 0.5 * 0.7 + 0.25 * 0.4;
        assert!(
            (p.z - plane).abs() < 1e-6,
            "the energy's minimiser is the plane: {p:?}"
        );
    }
}
