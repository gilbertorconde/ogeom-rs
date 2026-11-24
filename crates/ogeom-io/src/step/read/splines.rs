//! Spline knots and degrees as the file states or implies them, and spline patches joined.

use super::*;

/// The knots a STEP B-spline form implies when it states none: a Bézier
/// form's pieces joined at whole numbers, a uniform form's knots one apart
/// from `-degree`, a quasi-uniform form's the same clamped at both ends.
pub(super) fn implied_knots(form: &str, degree: usize, count: usize) -> Vec<f64> {
    let p = degree.max(1);
    #[allow(clippy::cast_precision_loss)]
    let at = |i: usize| i as f64;
    if form.starts_with("BEZIER") {
        let pieces = (count.saturating_sub(1) / p).max(1);
        let mut out = vec![0.0; p + 1];
        for s in 1..pieces {
            out.extend(std::iter::repeat_n(at(s), p));
        }
        out.extend(std::iter::repeat_n(at(pieces), p + 1));
        out
    } else if form.starts_with("QUASI_UNIFORM") {
        let interior = count.saturating_sub(p + 1);
        let mut out = vec![0.0; p + 1];
        out.extend((1..=interior).map(at));
        out.extend(std::iter::repeat_n(at(interior + 1), p + 1));
        out
    } else {
        #[allow(clippy::cast_precision_loss)]
        (0..count + p + 1).map(|i| i as f64 - p as f64).collect()
    }
}

/// Two splines joined end to start: the lower degree raised to the
/// higher, the second's weights scaled to meet the first's at the join (a
/// rational curve is unchanged by scaling all its weights), and the join
/// knot left at full multiplicity.
pub(super) fn join_splines(
    held: (KnotVector, Vec<Weighted<Point>>),
    next: &BSplineCurve,
    tol: Tolerances,
) -> OgeomResult<(KnotVector, Vec<Weighted<Point>>)> {
    let (mut knots, mut control) = held;
    let mut next = next.clone();
    while next.degree() < knots.degree() {
        next = next.elevated(tol)?;
    }
    while knots.degree() < next.degree() {
        (knots, control) = ogeom_math::bspline::elevate_degree(&knots, &control, tol)?;
    }
    let (Some(end), Some(start)) = (control.last(), next.control_points().first()) else {
        ogeom_bail!(Construction, "a composite segment has no control points");
    };
    let factor = end.weight / start.weight;
    let scaled: Vec<Weighted<Point>> = next
        .control_points()
        .iter()
        .map(|w| Weighted {
            scaled: w.scaled.scale(factor),
            weight: w.weight * factor,
        })
        .collect();
    ogeom_math::bspline::join(&(knots, control), &(next.knots().clone(), scaled))
}

/// One patch of a composite surface as a tensor spline: its knots each way,
/// on the unit interval, and its net indexed `[u][v]`.
#[derive(Debug, Clone)]
pub(super) struct Patch {
    u: KnotVector,
    v: KnotVector,
    net: Vec<Vec<Weighted<Point>>>,
}

impl Patch {
    pub(super) fn of(surface: &ogeom_geom::BSplineSurface) -> OgeomResult<Self> {
        let grid = surface.grid();
        let (nu, nv) = (grid.u_count(), grid.v_count());
        let points = grid.points();
        let net = (0..nu)
            .map(|i| (0..nv).map(|j| points[i * nv + j]).collect())
            .collect();
        Ok(Self {
            u: surface.u_knots().reparameterized(0.0, 1.0)?,
            v: surface.v_knots().reparameterized(0.0, 1.0)?,
            net,
        })
    }

    pub(super) fn transposed(self) -> Self {
        let (nu, nv) = (self.net.len(), self.net[0].len());
        let net = (0..nv)
            .map(|j| (0..nu).map(|i| self.net[i][j]).collect())
            .collect();
        Self {
            u: self.v,
            v: self.u,
            net,
        }
    }

    pub(super) fn reversed_u(mut self) -> Self {
        self.net.reverse();
        self.u = self.u.reversed();
        self
    }

    /// Apply a curve operation along `u` to every column of the net.
    fn along_u(
        self,
        op: impl Fn(&KnotVector, &[Weighted<Point>]) -> OgeomResult<(KnotVector, Vec<Weighted<Point>>)>,
    ) -> OgeomResult<Self> {
        let nv = self.net[0].len();
        let mut knots = self.u.clone();
        let mut columns = Vec::with_capacity(nv);
        for j in 0..nv {
            let column: Vec<Weighted<Point>> = self.net.iter().map(|row| row[j]).collect();
            let (k, c) = op(&self.u, &column)?;
            knots = k;
            columns.push(c);
        }
        let nu = columns[0].len();
        let net = (0..nu)
            .map(|i| (0..nv).map(|j| columns[j][i]).collect())
            .collect();
        Ok(Self {
            u: knots,
            v: self.v,
            net,
        })
    }
}

/// A grid of patches (`grid[i][j]`, `i` along `u`) joined into one surface.
pub(super) fn join_patches(
    mut grid: Vec<Vec<Patch>>,
    tol: Tolerances,
) -> OgeomResult<SurfaceGeometry> {
    let (m, n) = (grid.len(), grid[0].len());
    // Degrees raised to the grid's highest each way.
    let pu = grid
        .iter()
        .flatten()
        .map(|p| p.u.degree())
        .max()
        .unwrap_or(1);
    let pv = grid
        .iter()
        .flatten()
        .map(|p| p.v.degree())
        .max()
        .unwrap_or(1);
    for patch in grid.iter_mut().flatten() {
        let mut p = patch.clone();
        while p.u.degree() < pu {
            p = p.along_u(|k, c| ogeom_math::bspline::elevate_degree(k, c, tol))?;
        }
        p = p.transposed();
        while p.u.degree() < pv {
            p = p.along_u(|k, c| ogeom_math::bspline::elevate_degree(k, c, tol))?;
        }
        *patch = p.transposed();
    }
    // Knots unified: along u within each grid row of patches sharing an
    // index `i` (their u knots must agree to share v boundaries), along v
    // within each index `j`.
    let unify = |patches: Vec<Patch>, along_v: bool| -> OgeomResult<Vec<Patch>> {
        let turned: Vec<Patch> = if along_v {
            patches.into_iter().map(Patch::transposed).collect()
        } else {
            patches
        };
        let mut wanted: Vec<(f64, usize)> = Vec::new();
        for p in &turned {
            for (value, count) in p.u.distinct() {
                if value <= 1e-12 || value >= 1.0 - 1e-12 {
                    continue;
                }
                match wanted.iter_mut().find(|(v, _)| (v - value).abs() <= 1e-12) {
                    Some(slot) => slot.1 = slot.1.max(count),
                    None => wanted.push((value, count)),
                }
            }
        }
        let mut out = Vec::with_capacity(turned.len());
        for mut p in turned {
            for &(value, count) in &wanted {
                let have = p.u.multiplicity_of(value);
                if have < count {
                    p = p.along_u(|k, c| {
                        ogeom_math::bspline::insert_knot(k, c, value, count - have, tol)
                    })?;
                }
            }
            out.push(if along_v { p.transposed() } else { p });
        }
        Ok(out)
    };
    for row in &mut grid {
        *row = unify(std::mem::take(row), true)?;
    }
    for j in 0..n {
        let column: Vec<Patch> = grid.iter().map(|row| row[j].clone()).collect();
        for (row, p) in grid.iter_mut().zip(unify(column, false)?) {
            row[j] = p;
        }
    }
    // Joined: shared boundaries once, each join a knot of full multiplicity.
    let near = |a: &Weighted<Point>, b: &Weighted<Point>| {
        a.point().distance(b.point()) <= tol.confusion() * 100.0
            && (a.weight - b.weight).abs() <= 1e-9 * a.weight.abs().max(1.0)
    };
    let mut rows: Vec<Vec<Weighted<Point>>> = Vec::new();
    for (i, grid_row) in grid.iter().enumerate() {
        let nu = grid_row[0].net.len();
        for a in 0..nu {
            if i > 0 && a == 0 {
                // The row this patch shares with the one before it.
                let prev = rows.len() - 1;
                let mut shared = 0;
                for (j, patch) in grid_row.iter().enumerate() {
                    for (b, cell) in patch.net[0].iter().enumerate() {
                        if j > 0 && b == 0 {
                            continue;
                        }
                        if !near(&rows[prev][shared], cell) {
                            ogeom_bail!(Construction, "neighbouring patches part along u");
                        }
                        shared += 1;
                    }
                }
                continue;
            }
            let mut row = Vec::new();
            for (j, patch) in grid_row.iter().enumerate() {
                for (b, cell) in patch.net[a].iter().enumerate() {
                    if j > 0 && b == 0 {
                        if !near(row.last().unwrap_or(cell), cell) {
                            ogeom_bail!(Construction, "neighbouring patches part along v");
                        }
                        continue;
                    }
                    row.push(*cell);
                }
            }
            rows.push(row);
        }
    }
    // Each patch's knots shifted to its own unit of the whole; at a join
    // one copy of the shared end goes, leaving the multiplicity at the
    // degree, and the next patch's run follows past its own clamped start.
    let joined_knots = |parts: Vec<&KnotVector>, degree: usize| -> OgeomResult<KnotVector> {
        let mut out: Vec<f64> = Vec::new();
        for (k, knots) in parts.iter().enumerate() {
            #[allow(clippy::cast_precision_loss)]
            let shift = k as f64;
            let values: Vec<f64> = knots.knots().iter().map(|x| x + shift).collect();
            if k == 0 {
                out.extend(values);
            } else {
                out.pop();
                out.extend_from_slice(&values[degree + 1..]);
            }
        }
        KnotVector::new(out, degree)
    };
    let u = joined_knots((0..m).map(|i| &grid[i][0].u).collect(), pu)?;
    let v = joined_knots((0..n).map(|j| &grid[0][j].v).collect(), pv)?;
    let (nu, nv) = (rows.len(), rows[0].len());
    let cells: Vec<Weighted<Point>> = rows.into_iter().flatten().collect();
    Ok(
        ogeom_geom::BSplineSurface::rational(u, v, ogeom_math::ControlGrid::new(cells, nu, nv)?)?
            .into(),
    )
}

/// The highest spline degree a file may state. The standard sets no limit;
/// exporters stay far below this, and a degree past it would only size
/// allocations from a number in the file.
const MOST_DEGREE: usize = 64;

/// A spline degree read from the file, `missing` where the argument is
/// absent, refused past [`MOST_DEGREE`] or where it is not a count.
pub(super) fn file_degree(arg: Option<&Arg>, missing: f64) -> OgeomResult<usize> {
    let degree = arg.and_then(Arg::number).unwrap_or(missing);
    #[allow(clippy::cast_precision_loss)]
    if !degree.is_finite() || degree < 0.0 || degree > MOST_DEGREE as f64 {
        ogeom_bail!(Construction, "a spline of degree {degree}");
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    Ok(degree as usize)
}

/// Knots written once each with their multiplicities, spelt out. A
/// multiplicity past one more than [`MOST_DEGREE`] is refused: no valid
/// spline repeats a knot more often than its order.
pub(super) fn expanded_knots(mults: &[Arg], knots: &[Arg]) -> OgeomResult<Vec<f64>> {
    let mut out = Vec::new();
    for (m, k) in mults.iter().zip(knots) {
        let count = m.number().unwrap_or(1.0);
        #[allow(clippy::cast_precision_loss)]
        if !count.is_finite() || count < 0.0 || count > (MOST_DEGREE + 1) as f64 {
            ogeom_bail!(Construction, "a knot of multiplicity {count}");
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        out.extend(std::iter::repeat_n(
            k.number().unwrap_or(0.0),
            count as usize,
        ));
    }
    Ok(out)
}
