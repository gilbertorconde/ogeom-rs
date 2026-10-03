//! The nearest point on a surface.
//!
//! A foot point is where the displacement from the surface to a target is
//! perpendicular to both tangents. It is found by Newton from a seed, and
//! the seed is what decides which foot: a grid over the surface's spans,
//! a stored copy of that grid for many targets, or a caller's own guess.

use crate::{Surface, SurfaceGeometry, SurfaceJet};
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_math::Point;

/// Where a point projects onto a surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceProjection {
    /// The parameters of the nearest point found.
    pub parameters: (f64, f64),
    /// The nearest point.
    pub point: Point,
    /// The distance to it.
    pub distance: f64,
}

/// The nearest point on a surface to `target`.
///
/// A coarse grid to bracket, then Newton on the two conditions that define a
/// foot point: the displacement from the surface to the target is perpendicular
/// to both tangents. Grid resolution is `samples` per direction: the distance
/// to a surface is generally multi-modal, and the grid sets how fine a basin
/// it can tell apart.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if the surface cannot be
/// evaluated over its own domain.
pub fn project_on_surface(
    surface: &SurfaceGeometry,
    target: Point,
    samples: usize,
    tol: Tolerances,
) -> OgeomResult<SurfaceProjection> {
    // A plane, a cylinder and a sphere have their nearest point in closed
    // form, unique off the axis or the centre: where it lies inside the
    // surface's window, the grid below would find the same basin by a
    // thousand evaluations.
    if let Some(seed) = closed_form_foot(surface, target, tol)
        && let Ok(found) = refine_foot(surface, target, seed, tol)
    {
        return Ok(found);
    }
    let (us, vs) = seed_lines(surface, samples);

    let mut scan = Scan::default();
    for &u in &us {
        let mut row = Row::with_capacity(vs.len());
        for &v in &vs {
            let d = surface
                .point_at(u, v, tol)
                .map_or(f64::INFINITY, |p| p.square_distance(target));
            row.push((u, v, d));
        }
        scan.push_row(row);
    }

    scan.finish().refine(surface, target, tol)
}

/// The parameters of the nearest point where a closed form gives it and it
/// needs no clamping into the window: a plane's orthogonal foot, a
/// cylinder's for a point off its axis, a sphere's for a point off its
/// centre. `None` for any other surface, or a foot outside the window.
fn closed_form_foot(
    surface: &SurfaceGeometry,
    target: Point,
    tol: Tolerances,
) -> Option<(f64, f64)> {
    use ogeom_math::elementary;
    let (u, v) = match surface {
        SurfaceGeometry::Plane(p) => elementary::plane_parameters(&p.plane(), target),
        SurfaceGeometry::Cylinder(c) => {
            elementary::cylinder_parameters(&c.cylinder(), target, tol).ok()?
        }
        SurfaceGeometry::Sphere(s) => {
            elementary::sphere_parameters(&s.sphere(), target, tol).ok()?
        }
        _ => return None,
    };
    let ((u0, u1), (v0, v1)) = surface.domain();
    (u >= u0 && u <= u1 && v >= v0 && v <= v1).then_some((u, v))
}

/// One row of a seed scan: `(u, v, square distance)` per cell, a gap where
/// the surface would not evaluate.
type Row = smallvec::SmallVec<[(f64, f64, f64); 64]>;

/// A seed scan that keeps the grid's local minima, row by row.
///
/// Three rows are enough to know whether a cell of the middle one beats
/// its eight neighbours, so the scan never holds the grid: a projection
/// onto a thread flank seeds thousands of cells, and a caller projecting
/// every sample of an edge would pay that grid each time.
#[derive(Debug, Default)]
struct Scan {
    before: Option<Row>,
    last: Option<Row>,
    rows: usize,
    starts: Starts,
}

impl Scan {
    /// Take the next row; the previous row's minima are now decidable.
    fn push_row(&mut self, row: Row) {
        if let Some(last) = self.last.take() {
            self.starts
                .minima(self.rows - 1, self.before.as_ref(), &last, Some(&row));
            self.before = Some(last);
        }
        self.last = Some(row);
        self.rows += 1;
    }

    /// The last row's minima, then the picks.
    fn finish(mut self) -> Starts {
        if let Some(last) = self.last.take() {
            self.starts
                .minima(self.rows - 1, self.before.as_ref(), &last, None);
        }
        self.starts
    }
}

/// The few best basins of a grid scan.
///
/// The nearest seed is not always in the right basin. A helical flank
/// stacks its turns a pitch apart, and where the pitch is shorter than the
/// chord between neighbouring seeds, a seed on the turn above sits nearer
/// the target than the seed a half-span along the right turn; Newton then
/// converges faithfully on the wrong turn. So a scan keeps a handful of
/// candidates, one per basin (a cell that beats its eight neighbours),
/// and the refinement runs from each, nearest first, until one lands.
#[derive(Debug, Default, Clone, Copy)]
struct Starts {
    /// Nearest first.
    picks: [Option<Pick>; 4],
}

/// One basin of a seed scan.
#[derive(Debug, Clone, Copy)]
struct Pick {
    /// The grid cell the seed came from.
    cell: (usize, usize),
    /// The seed's parameters.
    at: (f64, f64),
    /// The seed's square distance to the target.
    d: f64,
}

impl Starts {
    /// Offer every cell of `row` (the grid's row `r`) that is no farther
    /// than any of its neighbours in the three rows around it.
    fn minima(&mut self, r: usize, before: Option<&Row>, row: &Row, after: Option<&Row>) {
        for (j, &(u, v, d)) in row.iter().enumerate() {
            if !d.is_finite() {
                continue;
            }
            let lo = j.saturating_sub(1);
            let hi = (j + 1).min(row.len() - 1);
            let beaten = |cells: &Row| cells[lo..=hi].iter().any(|c| c.2 < d);
            if beaten(row) || before.is_some_and(beaten) || after.is_some_and(beaten) {
                continue;
            }
            self.offer((r, j), (u, v), d);
        }
    }

    /// Consider a seed; it displaces an adjacent pick it beats, or the
    /// worst pick when it is nearer, and otherwise fills a free slot.
    fn offer(&mut self, cell: (usize, usize), at: (f64, f64), d: f64) {
        if !d.is_finite() {
            return;
        }
        let pick = Pick { cell, at, d };
        let adjacent = |a: (usize, usize)| a.0.abs_diff(cell.0) <= 1 && a.1.abs_diff(cell.1) <= 1;
        if let Some(slot) = self
            .picks
            .iter()
            .position(|p| p.is_some_and(|p| adjacent(p.cell)))
        {
            if self.picks[slot].is_some_and(|p| d < p.d) {
                self.picks[slot] = Some(pick);
                self.settle();
            }
            return;
        }
        if let Some(slot) = self.picks.iter().position(Option::is_none) {
            self.picks[slot] = Some(pick);
            self.settle();
        } else if self.picks[3].is_some_and(|p| d < p.d) {
            self.picks[3] = Some(pick);
            self.settle();
        }
    }

    /// Nearest first; a tie keeps the earlier pick.
    fn settle(&mut self) {
        self.picks.sort_by(|a, b| match (a, b) {
            (Some(a), Some(b)) => a.d.total_cmp(&b.d),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => std::cmp::Ordering::Equal,
        });
    }

    /// Newton from each pick, nearest first, keeping the closest foot; a
    /// foot within confusion ends the search, since nothing beats it.
    fn refine(
        self,
        surface: &SurfaceGeometry,
        target: Point,
        tol: Tolerances,
    ) -> OgeomResult<SurfaceProjection> {
        let ((ua, _), (va, _)) = surface.domain();
        let mut best: Option<SurfaceProjection> = None;
        for pick in self.picks.iter().flatten() {
            let found = refine_foot(surface, target, pick.at, tol)?;
            let better = best.as_ref().is_none_or(|b| found.distance < b.distance);
            if better {
                let done = found.distance <= tol.confusion();
                best = Some(found);
                if done {
                    break;
                }
            }
        }
        match best {
            Some(found) => Ok(found),
            None => refine_foot(surface, target, (ua, va), tol),
        }
    }
}

/// Where to seed a projection in each direction.
///
/// A fitted surface can carry hundreds of knot spans in one direction (a
/// thread flank swept hundreds of turns down a lead screw), and a
/// grid of sixteen or ninety-six seeds lands turns away from the nearest
/// point, where Newton converges faithfully onto the wrong flank. So a
/// patch is seeded *by its spans*, never by its domain: every span gets its
/// share of the caller's budget, one seed at the least, and a span a
/// thousand times wider than its neighbour gets no more for being wide.
/// Everything else has a domain that means what it says, and is seeded
/// evenly across it.
fn seed_lines(surface: &SurfaceGeometry, samples: usize) -> (Vec<f64>, Vec<f64>) {
    const CAP: usize = 4096;
    let base = samples.max(4);
    let ((ua, ub), (va, vb)) = surface.domain();
    let SurfaceGeometry::BSpline(spline) = surface else {
        return (spread(ua, ub, base), spread(va, vb, base));
    };
    // A patch's knots are where its shape is, and a file's knots are its
    // own business: an imported patch may put every knot but the first
    // inside the last small fraction of its domain, with its face in a
    // fraction of that. A grid spread evenly over the domain puts one seed
    // in the whole region the face lives in, and a projection seeded a
    // knot span away lands wherever Newton takes it, well off an edge that
    // sits on the surface. So the seeds follow the spans: each one gets
    // its share, however wide the file made it.
    (
        per_span(&breaks(spline.u_knots()), base, CAP),
        per_span(&breaks(spline.v_knots()), base, CAP),
    )
}

/// `count + 1` parameters evenly across `[from, to]`.
fn spread(from: f64, to: f64, count: usize) -> Vec<f64> {
    #[allow(clippy::cast_precision_loss)]
    (0..=count)
        .map(|i| from + (to - from) * (i as f64 / count as f64))
        .collect()
}

/// A knot vector's distinct values, without their multiplicities.
fn breaks(knots: &ogeom_math::KnotVector) -> Vec<f64> {
    knots.distinct().into_iter().map(|(at, _)| at).collect()
}

/// The budget shared out over the knot spans, ends included, once each.
fn per_span(knots: &[f64], budget: usize, cap: usize) -> Vec<f64> {
    let spans = knots.len().saturating_sub(1);
    if spans == 0 {
        return knots.to_vec();
    }
    let each = (budget / spans).min(cap / spans).max(1);
    let mut out = Vec::with_capacity(spans * each + 1);
    for pair in knots.windows(2) {
        out.extend(spread(pair[0], pair[1], each).into_iter().take(each));
    }
    out.push(knots[knots.len() - 1]);
    out
}

/// A surface's seeding grid, built once and asked many times.
///
/// [`project_on_surface`] evaluates the same grid of surface points for
/// every call: hundreds of evaluations per projection, identical each
/// time. A caller projecting *many* targets onto *one* surface builds the
/// grid once and each projection reduces to a nearest-seed scan plus the
/// Newton polish: the same seeds, the same refinement, the same answer to
/// the bit, at a fraction of the evaluations.
#[derive(Debug, Clone)]
pub struct SurfaceSeeds {
    /// `(parameters, point)` per cell, row by row; a gap where the surface
    /// would not evaluate.
    rows: Vec<Vec<(f64, f64, Option<Point>)>>,
}

impl SurfaceSeeds {
    /// Evaluate the grid `project_on_surface` would use, once.
    ///
    /// # Errors
    ///
    /// Never for a well-formed surface; evaluation failures leave gaps in
    /// the grid exactly as the per-call version tolerates them.
    pub fn over(surface: &SurfaceGeometry, samples: usize, tol: Tolerances) -> OgeomResult<Self> {
        let (us, vs) = seed_lines(surface, samples);
        let mut rows = Vec::with_capacity(us.len());
        for &u in &us {
            let mut row = Vec::with_capacity(vs.len());
            for &v in &vs {
                row.push((u, v, surface.point_at(u, v, tol).ok()));
            }
            rows.push(row);
        }
        Ok(Self { rows })
    }

    /// Project `target`, seeded from the stored grid, bit-identical to
    /// [`project_on_surface`] at the same sample count.
    ///
    /// # Errors
    ///
    /// As [`project_on_surface`].
    pub fn project(
        &self,
        surface: &SurfaceGeometry,
        target: Point,
        tol: Tolerances,
    ) -> OgeomResult<SurfaceProjection> {
        // The same closed-form foot [`project_on_surface`] takes first, so
        // the two answer alike to the bit.
        if let Some(seed) = closed_form_foot(surface, target, tol)
            && let Ok(found) = refine_foot(surface, target, seed, tol)
        {
            return Ok(found);
        }
        let mut scan = Scan::default();
        for row in &self.rows {
            scan.push_row(
                row.iter()
                    .map(|&(u, v, p)| {
                        (u, v, p.map_or(f64::INFINITY, |p| p.square_distance(target)))
                    })
                    .collect(),
            );
        }
        scan.finish().refine(surface, target, tol)
    }
}

/// The nearest point on a surface to `target`, starting from a guess.
///
/// [`project_on_surface`] brackets with a grid before refining. A caller
/// walking *along* something (the samples of a curve being projected into a
/// chart) already has a far better guess than any grid: where the previous
/// sample landed. Neighbouring samples of a curve are neighbouring points of
/// the surface, so the refinement starts inside the right basin and converges
/// in a few steps, and the grid's hundreds of evaluations per sample are not
/// spent at all.
///
/// The guess is load-bearing: this converges on the foot point nearest it, not
/// on the globally nearest one. A caller that cannot vouch for its guess
/// should check the reported distance and fall back to
/// [`project_on_surface`], which is what the exchange readers do.
///
/// # Errors
///
/// As [`project_on_surface`].
pub fn project_on_surface_from(
    surface: &SurfaceGeometry,
    target: Point,
    guess: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<SurfaceProjection> {
    refine_foot(surface, target, guess, tol)
}

/// Newton on the foot-point conditions from a starting parameter pair.
///
/// The iteration is box-constrained: the parameters never leave the
/// surface's domain. A periodic direction wraps; a bounded one clamps, and a
/// coordinate held against its bound by the step is pinned there while the
/// other one keeps solving on its own. Without that, every edge that runs
/// along a face's boundary (the outer helix of a thread flank sits exactly
/// on the flank's `u` bound) pushes the unconstrained foot a hair outside
/// the domain, and a solver that then rejects the whole answer hands back
/// its seed, turns away from the true foot.
///
/// The line search damps on the distance itself, which is the quantity a
/// projection minimises, so an iterate that stops improving is the answer
/// rather than a failure.
fn refine_foot(
    surface: &SurfaceGeometry,
    target: Point,
    start: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<SurfaceProjection> {
    let ((ua, ub), (va, vb)) = surface.domain();
    let periodic = (surface.is_periodic_u(), surface.is_periodic_v());
    let inside = |t: f64, a: f64, b: f64, wraps: bool| -> f64 {
        if wraps {
            a + (t - a).rem_euclid(b - a)
        } else {
            t.clamp(a, b)
        }
    };
    let square_distance = |u: f64, v: f64| -> f64 {
        surface
            .point_at(u, v, tol)
            .map_or(f64::INFINITY, |p| p.square_distance(target))
    };

    let mut x = (
        inside(start.0, ua, ub, periodic.0),
        inside(start.1, va, vb, periodic.1),
    );
    let mut best = (x.0, x.1, square_distance(x.0, x.1));

    for _ in 0..60 {
        // One evaluation, not three. Every value here comes from the same
        // jet, so they are consistent with each other, which is what a
        // Newton step needs. A surface with no second derivative (an
        // offset) steps by Gauss-Newton on its first: the curvature terms
        // drop out of the Jacobian, which still descends and converges to
        // the same foot, only more slowly.
        let jet = match surface.jet_at(x.0, x.1, tol) {
            Ok(jet) => jet,
            Err(_) => {
                let Ok((point, du, dv)) = surface.point_d1_at(x.0, x.1, tol) else {
                    break;
                };
                SurfaceJet {
                    point,
                    du,
                    dv,
                    d2u: ogeom_math::Vector::ZERO,
                    duv: ogeom_math::Vector::ZERO,
                    d2v: ogeom_math::Vector::ZERO,
                }
            }
        };
        let SurfaceJet {
            point: p,
            du,
            dv,
            d2u,
            duv,
            d2v,
        } = jet;
        let gap = p - target;
        // The foot point conditions: (S - target) . Su = 0 and (S - target) . Sv = 0.
        let r = [gap.dot(du), gap.dot(dv)];
        let j = [
            [du.dot(du) + gap.dot(d2u), du.dot(dv) + gap.dot(duv)],
            [du.dot(dv) + gap.dot(duv), dv.dot(dv) + gap.dot(d2v)],
        ];

        // A bounded coordinate sitting on its bound with the residual pushing
        // it further out is pinned: its condition cannot be met inside the
        // domain, and it drops out of the system.
        let pinned = |t: f64, a: f64, b: f64, wraps: bool, push: f64| -> bool {
            !wraps && ((t <= a && push < 0.0) || (t >= b && push > 0.0))
        };
        // The residual is the gradient of the half square distance, so the
        // descent pushes against it.
        let pin_u = pinned(x.0, ua, ub, periodic.0, -r[0]);
        let pin_v = pinned(x.1, va, vb, periodic.1, -r[1]);

        let free_norm = match (pin_u, pin_v) {
            (true, true) => 0.0,
            (true, false) => r[1].abs(),
            (false, true) => r[0].abs(),
            (false, false) => r[0].hypot(r[1]),
        };
        if free_norm <= tol.confusion() {
            break;
        }

        let delta = match (pin_u, pin_v) {
            (true, true) => break,
            (true, false) => {
                if j[1][1].abs() <= f64::EPSILON {
                    break;
                }
                [0.0, r[1] / j[1][1]]
            }
            (false, true) => {
                if j[0][0].abs() <= f64::EPSILON {
                    break;
                }
                [r[0] / j[0][0], 0.0]
            }
            (false, false) => {
                let Some(d) = solve_2x2(j, r) else {
                    break;
                };
                d
            }
        };
        if !delta[0].is_finite() || !delta[1].is_finite() {
            break;
        }

        // Damping: halve until the distance actually falls. Each candidate is
        // put back inside the domain first, so a step aimed past a bound
        // becomes a step to it.
        let mut scale = 1.0;
        let mut accepted = None;
        for _ in 0..30 {
            let candidate = (
                inside(delta[0].mul_add(-scale, x.0), ua, ub, periodic.0),
                inside(delta[1].mul_add(-scale, x.1), va, vb, periodic.1),
            );
            let d = square_distance(candidate.0, candidate.1);
            if d < best.2 {
                accepted = Some((candidate, d));
                break;
            }
            scale *= 0.5;
        }
        let Some((next, d)) = accepted else {
            break;
        };
        let step = (next.0 - x.0).hypot(next.1 - x.1);
        x = next;
        best = (x.0, x.1, d);
        if step <= tol.parametric() {
            break;
        }
    }

    let point = surface.point_at(best.0, best.1, tol)?;
    Ok(SurfaceProjection {
        parameters: (best.0, best.1),
        point,
        distance: point.distance(target),
    })
}

/// `j * d = r` for a two-by-two system, `None` when it is singular.
fn solve_2x2(j: [[f64; 2]; 2], r: [f64; 2]) -> Option<[f64; 2]> {
    let (row0, row1, rhs0, rhs1) = if j[0][0].abs() >= j[1][0].abs() {
        (j[0], j[1], r[0], r[1])
    } else {
        (j[1], j[0], r[1], r[0])
    };
    if row0[0].abs() <= f64::EPSILON * (row1[0].abs() + row0[1].abs()).max(1.0) {
        return None;
    }
    let factor = row1[0] / row0[0];
    let denom = factor.mul_add(-row0[1], row1[1]);
    if denom.abs() <= f64::EPSILON * row0[1].abs().max(1.0) {
        return None;
    }
    let d1 = factor.mul_add(-rhs0, rhs1) / denom;
    let d0 = d1.mul_add(-row0[1], rhs0) / row0[0];
    Some([d0, d1])
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::{CylinderSurface, OffsetSurface};
    use ogeom_math::{Cylinder, Frame};

    const T: Tolerances = Tolerances::millimetres();

    /// An offset surface has no second derivative, and its foot is still
    /// found: a point set off the offset drum along its normal projects
    /// back to the parameters it was set off from, at the distance it was
    /// set off by, wherever it lies between the seeds.
    #[test]
    fn a_point_off_an_offset_drum_projects_to_its_foot() {
        let basis: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 1.5, T).unwrap(), (-0.5, 5.5))
                .unwrap()
                .into();
        let offset = SurfaceGeometry::Offset(Box::new(OffsetSurface::new(basis, 0.5).unwrap()));
        for (u, v) in [(0.1, 5.0), (1.3, 0.37), (3.0, 2.2), (5.9, 4.91)] {
            let on = offset.point_at(u, v, T).unwrap();
            let target = Point::new(on.x * 1.15, on.y * 1.15, on.z);
            let found = project_on_surface(&offset, target, 24, T).unwrap();
            assert!(
                found.point.distance(on) < T.confusion(),
                "({u}, {v}) projected to {:?}",
                found.parameters
            );
            assert!((found.distance - 0.3).abs() < 1e-9, "{}", found.distance);
        }
    }
}
