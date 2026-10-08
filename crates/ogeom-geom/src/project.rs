//! The nearest point on a surface.
//!
//! A foot point is where the displacement from the surface to a target is
//! perpendicular to both tangents. It is found by Newton from a seed, and
//! the seed is what decides which foot: a closed form, a branch and bound
//! over a B-spline patch's spans, a grid over any other surface, a stored
//! copy of either for many targets, or a caller's own guess.

mod spans;

use crate::{Surface, SurfaceGeometry, SurfaceJet};
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_math::Point;
use spans::SpanBoxes;

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
/// A seed, then Newton on the two conditions that define a foot point: the
/// displacement from the surface to the target is perpendicular to both
/// tangents. An elementary or swept surface seeds in closed form or along
/// its profile. A B-spline patch seeds by branch and bound over its knot
/// spans, to within a hundredth of the confusion of the nearest point of
/// the whole patch, and ignores `samples`. Any other surface seeds from a grid of `samples`
/// per direction: the distance to a surface is generally multi-modal, and
/// the grid sets how fine a basin it can tell apart.
///
/// Where two feet of a B-spline patch are equally near to within a
/// hundredth of the confusion, either is a right answer, and the one
/// returned is fixed by the patch and the target: the one at the smaller
/// `u`, then `v`, among the points the search seeds Newton from.
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
    // An elementary surface has its nearest point in closed form, and a
    // swept one through its profile: where that foot lies inside the
    // surface's window, the grid below would find the same basin by a
    // thousand evaluations.
    if let Some(found) = closed_form(
        surface,
        target,
        Profile::of(surface, samples, tol).as_ref(),
        tol,
    ) {
        return Ok(found);
    }
    if let Some(boxes) = span_boxes(surface) {
        return boxes.project(surface, target, tol);
    }
    grid_projection(surface, target, samples, tol)
}

/// The span boxes of a B-spline patch whose hull bounds it; `None` for any
/// other surface.
fn span_boxes(surface: &SurfaceGeometry) -> Option<SpanBoxes> {
    match surface {
        SurfaceGeometry::BSpline(spline) => SpanBoxes::of(spline),
        _ => None,
    }
}

/// The nearest point from the seed grid alone.
fn grid_projection(
    surface: &SurfaceGeometry,
    target: Point,
    samples: usize,
    tol: Tolerances,
) -> OgeomResult<SurfaceProjection> {
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

/// The nearest point where a closed form or a profile scan finds it
/// inside the surface's window, polished by [`refine_foot`]; `None` for any
/// other surface, a foot outside the window, or a polish that fails.
///
/// `profile` is [`Profile::of`] the same surface.
fn closed_form(
    surface: &SurfaceGeometry,
    target: Point,
    profile: Option<&Profile>,
    tol: Tolerances,
) -> Option<SurfaceProjection> {
    if let Some(profile) = profile {
        return profile.project(surface, target, tol);
    }
    let seed = closed_form_foot(surface, target, surface.domain(), tol)?;
    refine_foot(surface, target, seed, tol).ok()
}

/// The parameters of an elementary surface's nearest point where they lie
/// in `window`: a plane's orthogonal foot, a cylinder's, cone's or torus's
/// for a point off the axis, a sphere's for a point off its centre. A
/// trimmed surface asks its basis, an angle shifted by whole turns into the
/// trim. `None` for any other surface, or a foot outside the window.
fn closed_form_foot(
    surface: &SurfaceGeometry,
    target: Point,
    window: ((f64, f64), (f64, f64)),
    tol: Tolerances,
) -> Option<(f64, f64)> {
    use ogeom_math::elementary;
    let (u, v) = match surface {
        SurfaceGeometry::Plane(p) => elementary::plane_parameters(&p.plane(), target),
        SurfaceGeometry::Cylinder(c) => {
            elementary::cylinder_parameters(&c.cylinder(), target, tol).ok()?
        }
        SurfaceGeometry::Cone(c) => elementary::cone_parameters(&c.cone(), target, tol).ok()?,
        SurfaceGeometry::Sphere(s) => {
            elementary::sphere_parameters(&s.sphere(), target, tol).ok()?
        }
        SurfaceGeometry::Torus(t) => elementary::torus_parameters(&t.torus(), target, tol).ok()?,
        SurfaceGeometry::Trimmed(t) => return closed_form_foot(t.basis(), target, window, tol),
        _ => return None,
    };
    let ((ua, ub), (va, vb)) = surface.domain();
    let u = into_window(u, window.0, surface.is_periodic_u().then_some(ub - ua))?;
    let v = into_window(v, window.1, surface.is_periodic_v().then_some(vb - va))?;
    Some((u, v))
}

/// `t` in `[from, to]`, shifted by whole periods where the direction has
/// one; `None` where no shift lands inside.
fn into_window(t: f64, (from, to): (f64, f64), period: Option<f64>) -> Option<f64> {
    if t >= from && t <= to {
        return Some(t);
    }
    let period = period.filter(|p| *p > 0.0)?;
    let shifted = ((from - t) / period).ceil().mul_add(period, t);
    (shifted >= from && shifted <= to).then_some(shifted)
}

/// A profile scan: per sample, its square distance across the sweep and
/// the seed it gives where the sweep coordinate is inside the window.
type ProfileScan = smallvec::SmallVec<[(f64, Option<(f64, f64)>); 64]>;

/// How a swept surface reduces to its profile.
#[derive(Debug, Clone, Copy)]
enum Sweep {
    /// An extrusion along a unit direction: `u` runs along the profile,
    /// `v` is the distance swept.
    Extrusion(ogeom_math::Vector),
    /// A revolution about an axis: `u` is the angle, `v` runs along the
    /// profile.
    Revolution(ogeom_math::Axis),
}

/// A swept surface's profile, sampled once for a one-dimensional foot scan.
///
/// The distance from a target to an extrusion's line through a profile
/// point is the distance across the sweep direction, and to a
/// revolution's circle through a profile point the distance in the
/// half-plane through the axis: each is the least over the sweep
/// coordinate, which then follows in closed form. So a scan along the
/// profile alone finds the basins the whole grid would, from one surface
/// line of samples instead of a grid of them.
#[derive(Debug, Clone)]
struct Profile {
    sweep: Sweep,
    /// The sweep coordinate of the sampled line: the window's first.
    at: f64,
    /// `(profile parameter, point on the line)` per sample, `None` where
    /// the surface would not evaluate.
    samples: Vec<(f64, Option<Point>)>,
}

impl Profile {
    /// The profile of an extrusion or a revolution, plain or trimmed,
    /// sampled as [`seed_lines`] seeds its profile direction; `None` for any
    /// other surface.
    fn of(surface: &SurfaceGeometry, samples: usize, tol: Tolerances) -> Option<Self> {
        let basis = match surface {
            SurfaceGeometry::Trimmed(t) => t.basis(),
            other => other,
        };
        let sweep = match basis {
            SurfaceGeometry::Extrusion(e) => Sweep::Extrusion(e.direction().vector()),
            SurfaceGeometry::Revolution(r) => Sweep::Revolution(r.axis()),
            _ => return None,
        };
        let (us, vs) = seed_lines(surface, samples);
        let ((ua, _), (va, _)) = surface.domain();
        let (at, along) = match sweep {
            Sweep::Extrusion(_) => (va, us),
            Sweep::Revolution(_) => (ua, vs),
        };
        let samples = along
            .into_iter()
            .map(|t| {
                let point = match sweep {
                    Sweep::Extrusion(_) => surface.point_at(t, at, tol),
                    Sweep::Revolution(_) => surface.point_at(at, t, tol),
                };
                (t, point.ok())
            })
            .collect();
        Some(Self { sweep, at, samples })
    }

    /// A sample's square distance to `target` across the sweep, and the
    /// surface parameters of the nearest point on its sweep line or circle;
    /// `None` for a target on a revolution's axis, where every angle is
    /// as near.
    fn reduce(
        &self,
        t: f64,
        point: Point,
        target: Point,
        tol: Tolerances,
    ) -> Option<(f64, (f64, f64))> {
        match self.sweep {
            Sweep::Extrusion(d) => {
                let gap = point - target;
                let along = gap.dot(d);
                let across = gap - d * along;
                Some((across.dot(across), (t, self.at - along)))
            }
            Sweep::Revolution(axis) => {
                let d = axis.direction.vector();
                let split = |p: Point| {
                    let off = p - axis.location;
                    let height = off.dot(d);
                    (off - d * height, height)
                };
                let (to, height_to) = split(target);
                let (from, height_from) = split(point);
                let (ring_to, ring_from) = (to.magnitude(), from.magnitude());
                if ring_to <= tol.confusion() {
                    return None;
                }
                let turn = from.cross(to).dot(d).atan2(from.dot(to));
                let gap = (ring_to - ring_from).hypot(height_to - height_from);
                Some((gap * gap, (self.at + turn, t)))
            }
        }
    }

    /// The nearest point, from every basin of the profile scan whose sweep
    /// coordinate lies in the window, nearest first; `None` where the
    /// nearest basin's lies outside it, or nothing evaluates.
    fn project(
        &self,
        surface: &SurfaceGeometry,
        target: Point,
        tol: Tolerances,
    ) -> Option<SurfaceProjection> {
        let ((ua, ub), (va, vb)) = surface.domain();
        let (sweep_window, period) = match self.sweep {
            Sweep::Extrusion(_) => ((va, vb), surface.is_periodic_v().then_some(vb - va)),
            Sweep::Revolution(_) => ((ua, ub), surface.is_periodic_u().then_some(ub - ua)),
        };
        let mut scan = ProfileScan::with_capacity(self.samples.len());
        for &(t, point) in &self.samples {
            let (d, seed) = match point {
                Some(point) => {
                    let (d, (u, v)) = self.reduce(t, point, target, tol)?;
                    // The sweep coordinate inside the window, or no seed.
                    let seed = match self.sweep {
                        Sweep::Extrusion(_) => into_window(v, sweep_window, period).map(|v| (u, v)),
                        Sweep::Revolution(_) => {
                            into_window(u, sweep_window, period).map(|u| (u, v))
                        }
                    };
                    (d, seed)
                }
                None => (f64::INFINITY, None),
            };
            scan.push((d, seed));
        }
        // The nearest sample decides: where its foot across the sweep is
        // outside the window, the nearest point is on the window's border,
        // and the grid finds it.
        let nearest = scan
            .iter()
            .filter(|s| s.0.is_finite())
            .min_by(|a, b| a.0.total_cmp(&b.0))?;
        nearest.1?;
        let mut starts = Starts::default();
        for (i, &(d, seed)) in scan.iter().enumerate() {
            let Some(at) = seed else { continue };
            let lo = i.saturating_sub(1);
            let hi = (i + 1).min(scan.len() - 1);
            if scan[lo..=hi].iter().any(|c| c.0 < d) {
                continue;
            }
            starts.offer((i, 0), at, d);
        }
        starts.refine(surface, target, tol).ok()
    }
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

/// A surface's seeding structure, built once and asked many times.
///
/// [`project_on_surface`] builds the same seeding structure for every
/// call: a B-spline patch's span boxes, or a grid of surface points
/// for a surface with no closed form. A caller projecting *many* targets
/// onto *one* surface builds it once, and each projection reduces to the
/// search plus the Newton polish: the same seeds, the same refinement, the
/// same answer to the bit.
#[derive(Debug, Clone)]
pub struct SurfaceSeeds {
    /// `(parameters, point)` per cell, row by row; a gap where the surface
    /// would not evaluate. Empty where span boxes answer.
    rows: Vec<Vec<(f64, f64, Option<Point>)>>,
    /// A B-spline patch's span boxes, which answer in place of the grid.
    spans: Option<SpanBoxes>,
    /// A swept surface's profile scan, which answers in place of the grid.
    profile: Option<Profile>,
}

impl SurfaceSeeds {
    /// Build the seeding structure `project_on_surface` would use, once.
    ///
    /// # Errors
    ///
    /// Never for a well-formed surface; evaluation failures leave gaps in
    /// the grid exactly as the per-call version tolerates them.
    pub fn over(surface: &SurfaceGeometry, samples: usize, tol: Tolerances) -> OgeomResult<Self> {
        let profile = Profile::of(surface, samples, tol);
        if let Some(spans) = span_boxes(surface) {
            return Ok(Self {
                rows: Vec::new(),
                spans: Some(spans),
                profile,
            });
        }
        let (us, vs) = seed_lines(surface, samples);
        let mut rows = Vec::with_capacity(us.len());
        for &u in &us {
            let mut row = Vec::with_capacity(vs.len());
            for &v in &vs {
                row.push((u, v, surface.point_at(u, v, tol).ok()));
            }
            rows.push(row);
        }
        Ok(Self {
            rows,
            spans: None,
            profile,
        })
    }

    /// Project `target`, seeded from the stored structure, bit-identical to
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
        // The same closed form [`project_on_surface`] takes first, so the
        // two answer alike to the bit.
        if let Some(found) = closed_form(surface, target, self.profile.as_ref(), tol) {
            return Ok(found);
        }
        if let Some(spans) = &self.spans {
            return spans.project(surface, target, tol);
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
        // The residual is the gap times a derivative, so it shrinks with the
        // chart's own speed: round a ring a sixth of a millimetre in radius,
        // a residual under the confusion still leaves the foot several times
        // the confusion along the surface. The gap's own tangential part, a
        // length, must be under it as well.
        let along = |r: f64, d: ogeom_math::Vector| -> f64 {
            let speed = d.magnitude();
            if speed > f64::MIN_POSITIVE {
                (r / speed).abs()
            } else {
                0.0
            }
        };
        let free_gap = match (pin_u, pin_v) {
            (true, true) => 0.0,
            (true, false) => along(r[1], dv),
            (false, true) => along(r[0], du),
            (false, false) => along(r[0], du).hypot(along(r[1], dv)),
        };
        if free_norm <= tol.confusion() && free_gap <= tol.confusion() {
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

    /// The closed forms and the profile scan answer at least as near as the
    /// grid they replace, on every surface they cover: a cone (both nappes
    /// in reach), a spindle torus, a trimmed drum whose window starts past a
    /// whole turn, an extrusion of a skew spline, a full and a partial
    /// revolution of a spline. Targets are spread over a box around each,
    /// off the axis and on it; each answer is a foot whatever the path, and
    /// a stored grid answers to the bit as the per-call projection does.
    #[test]
    fn closed_forms_are_no_farther_than_the_grid() {
        use crate::{
            ConeSurface, Curve, ExtrusionSurface, RevolutionSurface, TorusSurface, TrimmedSurface,
            curve::BSplineCurve,
        };
        use ogeom_math::{Axis, Cone, Direction, KnotVector, Torus};
        let tilted = Frame::new(
            Point::new(0.3, -0.2, 0.1),
            Direction::from_coords(0.2, 0.1, 1.0, T).unwrap(),
            Direction::from_coords(1.0, 0.0, 0.0, T).unwrap(),
            T,
        )
        .unwrap();
        let spline = Curve::BSpline(
            BSplineCurve::new(
                KnotVector::clamped_uniform(3, 6).unwrap(),
                vec![
                    Point::new(1.0, 0.0, -1.0),
                    Point::new(1.6, 0.4, -0.5),
                    Point::new(0.7, -0.3, 0.0),
                    Point::new(1.9, 0.2, 0.4),
                    Point::new(1.2, 0.5, 0.9),
                    Point::new(0.5, 0.1, 1.4),
                ],
                T,
            )
            .unwrap(),
        );
        let skew = Direction::from_coords(0.1, 0.3, 1.0, T).unwrap();
        let z = Axis::new(
            Point::ORIGIN,
            Direction::from_coords(0.0, 0.0, 1.0, T).unwrap(),
        );
        let drum: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(tilted, 1.0, T).unwrap(), (-1.0, 1.0))
                .unwrap()
                .into();
        let surfaces: Vec<SurfaceGeometry> = vec![
            ConeSurface::new(Cone::new(tilted, 0.4, 0.5, T).unwrap(), (-2.0, 1.5))
                .unwrap()
                .into(),
            TorusSurface::new(Torus::new(tilted, 1.0, 1.3, T).unwrap()).into(),
            SurfaceGeometry::Trimmed(Box::new(
                TrimmedSurface::new(drum, (7.0, 9.5), (-0.8, 0.6), T).unwrap(),
            )),
            SurfaceGeometry::Extrusion(Box::new(
                ExtrusionSurface::over(spline.clone(), skew, (-1.0, 1.5)).unwrap(),
            )),
            SurfaceGeometry::Revolution(Box::new(
                RevolutionSurface::new(spline.clone(), z, core::f64::consts::TAU).unwrap(),
            )),
            SurfaceGeometry::Revolution(Box::new(RevolutionSurface::new(spline, z, 2.0).unwrap())),
        ];
        for surface in &surfaces {
            let mut closed = 0;
            let seeds = SurfaceSeeds::over(surface, 24, T).unwrap();
            for i in 0..7 {
                for j in 0..7 {
                    for k in 0..5 {
                        let target = Point::new(
                            f64::from(i).mul_add(0.55, -1.6),
                            f64::from(j).mul_add(0.5, -1.4),
                            f64::from(k).mul_add(0.7, -1.3),
                        );
                        let found = project_on_surface(surface, target, 24, T).unwrap();
                        let profile = Profile::of(surface, 24, T);
                        if closed_form(surface, target, profile.as_ref(), T).is_some() {
                            closed += 1;
                        }
                        let grid = grid_projection(surface, target, 24, T).unwrap();
                        assert!(
                            found.distance <= grid.distance + 1e-9,
                            "{:?} at {target:?}: {} against the grid's {}",
                            surface.kind(),
                            found.distance,
                            grid.distance
                        );
                        assert_eq!(seeds.project(surface, target, T).unwrap(), found);
                    }
                }
            }
            assert!(closed > 0, "{:?} never took a closed form", surface.kind());
        }
    }

    /// A B-spline patch from a net of `nu` by `nv` points, each with a
    /// weight; uniform clamped knots unless `knots` gives its own.
    fn patch(
        degree: (usize, usize),
        knots: Option<(Vec<f64>, Vec<f64>)>,
        (nu, nv): (usize, usize),
        at: impl Fn(f64, f64) -> (Point, f64),
    ) -> SurfaceGeometry {
        use ogeom_math::{ControlGrid, KnotVector, Weighted};
        let (ku, kv) = match knots {
            Some((u, v)) => (
                KnotVector::new(u, degree.0).unwrap(),
                KnotVector::new(v, degree.1).unwrap(),
            ),
            None => (
                KnotVector::clamped_uniform(degree.0, nu).unwrap(),
                KnotVector::clamped_uniform(degree.1, nv).unwrap(),
            ),
        };
        #[allow(clippy::cast_precision_loss)]
        let net = (0..nu)
            .flat_map(|i| {
                (0..nv).map(move |j| (i as f64 / (nu - 1) as f64, j as f64 / (nv - 1) as f64))
            })
            .map(|(s, t)| {
                let (p, w) = at(s, t);
                Weighted::new(p, w, T).unwrap()
            })
            .collect();
        crate::BSplineSurface::rational(ku, kv, ControlGrid::new(net, nu, nv).unwrap())
            .unwrap()
            .into()
    }

    /// The nearest point by brute force: the nearest of a dense grid of
    /// surface points, `per` per span and direction, polished by Newton.
    fn brute_force(surface: &SurfaceGeometry, target: Point, per: usize) -> SurfaceProjection {
        let SurfaceGeometry::BSpline(spline) = surface else {
            unreachable!()
        };
        let lines = |k: &ogeom_math::KnotVector| -> Vec<f64> {
            let (a, b) = k.domain();
            let breaks: Vec<f64> = breaks(k)
                .into_iter()
                .filter(|x| *x >= a && *x <= b)
                .collect();
            let mut out: Vec<f64> = breaks
                .windows(2)
                .flat_map(|w| spread(w[0], w[1], per).into_iter().take(per))
                .collect();
            out.push(breaks[breaks.len() - 1]);
            out
        };
        let (us, vs) = (lines(spline.u_knots()), lines(spline.v_knots()));
        let mut best = (f64::INFINITY, (0.0, 0.0));
        for &u in &us {
            for &v in &vs {
                let d = surface.point_at(u, v, T).unwrap().distance(target);
                if d < best.0 {
                    best = (d, (u, v));
                }
            }
        }
        let polished = refine_foot(surface, target, best.1, T).unwrap();
        assert!(polished.distance <= best.0);
        polished
    }

    /// The branch and bound over a B-spline patch's spans finds the
    /// nearest point of the whole patch: never farther than a dense grid
    /// polished by Newton, on a finely knotted wavy patch, one with its
    /// knots crowded into a sliver of the domain, a rational patch with
    /// unclamped knots, a patch with an edge collapsed to a pole, a closed
    /// patch with targets beside its seam and beyond its open edges, and a
    /// patch whose two bumps stand equally near a target between them.
    /// Every foot is a point on the surface at its own parameters, and a
    /// stored set of span boxes answers as the per-call projection does, to the
    /// bit.
    #[test]
    fn a_spline_foot_is_the_nearest_point_of_the_whole_patch() {
        let wavy = patch((3, 3), None, (43, 43), |s, t| {
            let (x, y) = (20.0 * s, 20.0 * t);
            (
                Point::new(x, y, 0.8 * (0.7 * x).sin() * (0.45 * y).cos()),
                1.0,
            )
        });
        // Every knot but the ends inside the last thirtieth of `u`.
        let mut crowded_u = vec![0.0; 4];
        crowded_u.extend((1..20).map(|k| 0.97 + 0.03 * f64::from(k) / 20.0));
        crowded_u.extend([1.0; 4]);
        let crowded = patch(
            (3, 3),
            Some((crowded_u, vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0])),
            (23, 4),
            |s, t| {
                let (x, y) = (10.0 * s, 4.0 * t);
                (Point::new(x, y, (1.3 * x).sin() + 0.2 * y * y), 1.0)
            },
        );
        let unclamped = patch(
            (2, 3),
            Some((
                (0..12).map(f64::from).collect(),
                (0..13).map(|k| 0.5 * f64::from(k)).collect(),
            )),
            (9, 9),
            |s, t| {
                let (x, y) = (6.0 * s, 6.0 * t);
                let w = 0.6 + 1.2 * (s * t + 0.3 * (5.0 * s).sin().abs());
                (Point::new(x, y, (x - 3.0) * (y - 3.0) / 4.0), w)
            },
        );
        let pole = patch((3, 2), None, (8, 9), |s, t| {
            let (angle, r) = (std::f64::consts::PI * 1.5 * t, 3.0 * s);
            (
                Point::new(r * angle.cos(), r * angle.sin(), 2.0 - 2.0 * s),
                1.0,
            )
        });
        let seam = patch((3, 3), None, (13, 6), |s, t| {
            let angle = std::f64::consts::TAU * s;
            let r = 2.0 + 0.3 * (3.0 * angle).sin();
            (Point::new(r * angle.cos(), r * angle.sin(), 4.0 * t), 1.0)
        });
        let bumps = patch((3, 3), None, (21, 9), |s, t| {
            let bump = |c: f64| (-((s - c) * (s - c)) / 0.004).exp();
            (Point::new(s, 0.4 * t, 0.3 * (bump(0.3) + bump(0.7))), 1.0)
        });
        let targets = |lo: Point, hi: Point| -> Vec<Point> {
            (0..27)
                .map(|k| {
                    let f = |n: i32| f64::from((k / n) % 3) / 2.0;
                    Point::new(
                        lo.x + (hi.x - lo.x) * f(1),
                        lo.y + (hi.y - lo.y) * f(3),
                        lo.z + (hi.z - lo.z) * f(9),
                    )
                })
                .collect()
        };
        let cases: Vec<(&str, SurfaceGeometry, Vec<Point>)> = vec![
            (
                "wavy",
                wavy,
                targets(Point::new(-3.0, -2.0, -4.0), Point::new(23.0, 22.0, 4.0)),
            ),
            (
                "crowded",
                crowded,
                targets(Point::new(-1.0, -1.0, -2.0), Point::new(11.0, 5.0, 3.0)),
            ),
            (
                "unclamped",
                unclamped,
                targets(Point::new(-1.0, -1.0, -3.0), Point::new(7.0, 7.0, 3.0)),
            ),
            (
                "pole",
                pole,
                targets(Point::new(-3.5, -3.5, -1.0), Point::new(3.5, 3.5, 3.0)),
            ),
            (
                "seam",
                seam,
                [
                    targets(Point::new(1.0, -0.4, -1.0), Point::new(3.0, 0.4, 5.0)),
                    targets(Point::new(-3.0, -3.0, 1.0), Point::new(3.0, 3.0, 3.0)),
                ]
                .concat(),
            ),
            (
                "bumps",
                bumps,
                targets(Point::new(0.4, -0.2, 0.5), Point::new(0.6, 0.6, 1.5)),
            ),
        ];
        for (name, surface, targets) in &cases {
            let seeds = SurfaceSeeds::over(surface, 16, T).unwrap();
            for &target in targets {
                let found = project_on_surface(surface, target, 16, T).unwrap();
                let reference = brute_force(surface, target, 8);
                assert!(
                    found.distance <= reference.distance + T.confusion(),
                    "{name} at {target:?}: {} at {:?}, the brute force {} at {:?}",
                    found.distance,
                    found.parameters,
                    reference.distance,
                    reference.parameters
                );
                let (u, v) = found.parameters;
                let on = surface.point_at(u, v, T).unwrap();
                assert!(on.distance(found.point) < 1e-12, "{name} at {target:?}");
                assert_eq!(seeds.project(surface, target, T).unwrap(), found, "{name}");
            }
        }
    }

    /// Two feet equally near: the target stands over the middle of a patch
    /// symmetric about `u = 1/2`, between two equal bumps. Either foot is
    /// right; the one returned is the same on every call and from a stored
    /// set of span boxes, and it is on one bump or the other.
    #[test]
    fn a_tie_between_two_feet_resolves_the_same_way_every_time() {
        let bumps = patch((3, 3), None, (21, 9), |s, t| {
            let bump = |c: f64| (-((s - c) * (s - c)) / 0.004).exp();
            (Point::new(s, 0.4 * t, 0.3 * (bump(0.3) + bump(0.7))), 1.0)
        });
        let target = Point::new(0.5, 0.2, 0.9);
        let first = project_on_surface(&bumps, target, 16, T).unwrap();
        let reference = brute_force(&bumps, target, 8);
        assert!(first.distance <= reference.distance + T.confusion());
        let u = first.parameters.0;
        assert!((u - 0.5).abs() > 0.1, "landed between the bumps at {u}");
        let seeds = SurfaceSeeds::over(&bumps, 16, T).unwrap();
        for _ in 0..3 {
            assert_eq!(project_on_surface(&bumps, target, 16, T).unwrap(), first);
            assert_eq!(seeds.project(&bumps, target, T).unwrap(), first);
        }
    }

    /// A foot seeded beside a point on a thin drum lands on the point.
    ///
    /// Round a drum a sixth of a millimetre in radius, a radian of the chart
    /// is a sixth of a millimetre of the surface, and the residual (the gap
    /// times the chart's speed) is under the confusion while the foot still
    /// stands three confusions off along the surface. The seed is that
    /// close; the foot must not stay there.
    #[test]
    fn a_foot_on_a_thin_drum_converges_along_the_surface() {
        let drum: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 0.16, T).unwrap(), (-1.0, 1.0))
                .unwrap()
                .into();
        for (u, v) in [(0.5, 0.0), (3.0, 0.4), (6.0, -0.7)] {
            let target = drum.point_at(u, v, T).unwrap();
            let found = project_on_surface_from(&drum, target, (u + 2e-6, v), T).unwrap();
            assert!(
                found.distance <= T.confusion(),
                "({u}, {v}) stopped {:.2e} off at {:?}",
                found.distance,
                found.parameters
            );
        }
    }
}
