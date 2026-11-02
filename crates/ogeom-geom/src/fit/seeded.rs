//! Scattered points fitted over a seed surface's chart.
//!
//! [`fit_surface_scattered`](super::fit_surface_scattered) reads its
//! parameters off the cloud's principal plane, which folds anything that is
//! not a height field over that plane. Here they come from a seed surface
//! instead: each point's parameters are those of its nearest point on the
//! seed, so a band round a cylinder, a cap on a sphere or a sheet over a
//! previous fit unfolds the way the seed does.

use ogeom_core::{OgeomError, OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{KnotVector, Point};

use super::{fit_surface_scattered_at, precise};
use crate::project::{SurfaceSeeds, project_on_surface_from};
use crate::traits::Surface as _;
use crate::{BSplineSurface, SurfaceGeometry};

/// The seeding grid's samples a direction, for a seed with neither a closed
/// form nor span boxes: the count the projector's own callers use.
const SEED_SAMPLES: usize = 16;

/// The upper quartile of the standard normal distribution, `Φ⁻¹(3/4)`: the
/// median of the absolute value of a normal deviate, in deviations.
const NORMAL_QUARTILE: f64 = 0.674_489_750_196_081_7;

/// How far [`fit_surface_scattered_over`] refines its control net.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Target {
    /// Fit once, over the starting spans.
    Fixed,
    /// Split the spans holding a point farther from the fit than this
    /// distance, within the span bounds.
    Tolerance(f64),
}

/// What [`fit_surface_scattered_over`] fits with.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScatteredOptions {
    /// The degree in both directions.
    pub degree: usize,
    /// The knot spans to start from, `u` then `v`, evenly spaced.
    pub spans: (usize, usize),
    /// The most knot spans refinement may reach, `u` then `v`.
    pub most: (usize, usize),
    /// The thin-plate fairing weight, as in
    /// [`fit_surface_scattered_at`](super::fit_surface_scattered_at), over
    /// the unit square the fit's parameters run over.
    ///
    /// The energy pulls the fit towards a plane in those parameters, so it
    /// bends a curved shape as well as settling the controls no point
    /// reaches. A shape of extent `L` turning through `θ` radians has
    /// second derivatives of order `L θ` over the square; a weight well
    /// under `(tolerance / (L θ))²` leaves it alone.
    pub fairing: f64,
    /// Whether, and to what, the net is refined.
    pub target: Target,
}

/// Why [`fit_surface_scattered_over`] stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ended {
    /// [`Target::Fixed`]: there was nothing to refine to.
    Fixed,
    /// Every point is within the tolerance.
    Met,
    /// The spans that miss are at their bounds, or their cells hold too
    /// few points to settle a split.
    Bounds,
    /// A split would fit the noise: none takes more squared distance off
    /// than twice the noise's variance per control it adds.
    Noise,
    /// Only strays miss: every point beyond the tolerance stands farther
    /// off than the noise its cell shows reaches, and a split would bend
    /// the fit to it rather than follow a shape.
    Strays,
}

/// What [`fit_surface_scattered_over`] produced.
#[derive(Debug, Clone, PartialEq)]
pub struct ScatteredFit {
    /// The surface, over the unit square.
    pub surface: BSplineSurface,
    /// Each point's parameters: its nearest point on the surface, found
    /// from where the seed put it. The distances reported are to these.
    pub parameters: Vec<(f64, f64)>,
    /// The root mean square distance from the points to the surface.
    pub sigma: f64,
    /// The root mean square distance from each point to a fit, over the
    /// same knots and fairing, of the other half of the points (alternate
    /// points, by index). A fit that follows the shape predicts the points
    /// it did not see about as closely as those it did; a fit that follows
    /// the noise predicts them worse, and this stands well above
    /// [`sigma`](Self::sigma).
    pub cross_sigma: f64,
    /// The point farthest from the surface: its index and its distance.
    pub worst: (usize, f64),
    /// The knot spans, `u` then `v`.
    pub spans: (usize, usize),
    /// Why refinement stopped.
    pub ended: Ended,
}

impl ScatteredFit {
    /// Whether the tolerance was met.
    #[must_use]
    pub fn met(&self) -> bool {
        self.ended == Ended::Met
    }
}

/// A fit over given knots, its points at their feet on it.
struct Round {
    surface: BSplineSurface,
    u_knots: KnotVector,
    v_knots: KnotVector,
    /// Each point's foot on it, and the distance there.
    feet: Vec<(f64, f64)>,
    distances: Vec<f64>,
}

impl Round {
    fn square_sum(&self) -> f64 {
        self.distances.iter().map(|d| d * d).sum()
    }

    fn controls(&self) -> usize {
        self.u_knots.control_point_count() * self.v_knots.control_point_count()
    }
}

/// Fit a surface to scattered points, parameterised over a seed surface.
///
/// Each point takes the parameters of its nearest point on `seed`: a
/// closed form on a plane, cylinder, cone, sphere or torus (the whole
/// surface, whatever window it carries), the projector's search on anything
/// else, a previous fit included. A periodic direction is cut at the widest
/// gap the points leave round it, so a band across the seed's seam stays in
/// one piece; a cloud that runs all the way round is fitted open across
/// that gap. The parameters are then scaled onto the unit square. The seed
/// must unfold the points: a sphere seed's poles belong away from them,
/// where longitude still tells points apart.
///
/// The net is fitted at those parameters by
/// [`fit_surface_scattered_at`](super::fit_surface_scattered_at), faired as
/// it fairs. Every distance reported is from a point to its nearest point
/// on the fit, found from its parameters.
///
/// With [`Target::Tolerance`] the spans are then split where the fit
/// misses. A cell (a `u` span by a `v` span) misses where one of its points
/// stands farther from the fit than the tolerance and is not a stray, one
/// standing farther than the noise of the rest of its cell (their root mean
/// square distance) reaches over `m` points, `√(2 ln m)` times it
/// ([`Ended::Strays`]). A cell is split only where it holds at least as
/// many points as the `(degree + 1)²` controls that bear on it, fewer
/// leaving its controls to the fairing.
///
/// Each round splits the worst-missed cell's `u` span, its `v` span, or
/// both, at the middle, keeping the try with the least squared distances
/// plus twice the noise's variance per control (Mallows' criterion, the
/// variance estimated from the fit's median distance), so a shape that
/// turns along one direction only gains spans along that one. A cell
/// beside the worst can miss by the ringing of a fit that has not yet
/// caught it, which a split at the worst takes out. Where no try scores
/// below the fit so far, the split would fit noise rather than shape, and
/// refinement ends there ([`Ended::Noise`]), as it does for a tolerance
/// inside the noise's reach.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// there are fewer than eight points, a point is not finite, the degree is
/// nought, a span count is nought or its bound below it, the fairing is not
/// finite and non-negative, the tolerance is not a distance, a point has no
/// nearest point on the seed (on a cylinder's axis, at a sphere's centre),
/// or the points' parameters do not spread over both directions;
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if the whole
/// cloud cannot settle the net (raise the fairing).
pub fn fit_surface_scattered_over(
    points: &[Point],
    seed: &SurfaceGeometry,
    options: &ScatteredOptions,
    tol: Tolerances,
) -> OgeomResult<ScatteredFit> {
    let ScatteredOptions {
        degree,
        spans,
        most,
        fairing,
        target,
    } = *options;
    if points.len() < 8 {
        ogeom_bail!(
            Construction,
            "a seeded fit needs at least eight points, four to each half"
        );
    }
    if let Some(p) = points
        .iter()
        .find(|p| !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()))
    {
        ogeom_bail!(Construction, "the point {p:?} is not finite");
    }
    if degree == 0 {
        ogeom_bail!(Construction, "a fit needs a degree of at least one");
    }
    if spans.0 == 0 || spans.1 == 0 || most.0 < spans.0 || most.1 < spans.1 {
        ogeom_bail!(
            Construction,
            "{}x{} starting spans within {}x{} is not a span budget",
            spans.0,
            spans.1,
            most.0,
            most.1
        );
    }
    if !fairing.is_finite() || fairing < 0.0 {
        ogeom_bail!(
            Construction,
            "a fairing weight of {fairing} is not a weight"
        );
    }
    let tolerance = match target {
        Target::Fixed => None,
        Target::Tolerance(t) => {
            if !(t > 0.0 && t.is_finite()) {
                ogeom_bail!(Construction, "a tolerance of {t} is not a distance");
            }
            Some(t)
        }
    };

    let start = seeded_parameters(points, seed, tol)?;
    let uniform = |count: usize| KnotVector::clamped_uniform(degree, count + degree);
    let mut current = round(
        points,
        &start,
        (uniform(spans.0)?, uniform(spans.1)?),
        fairing,
        tol,
    )?;
    let Some(tolerance) = tolerance else {
        return finish(points, &start, current, Ended::Fixed, fairing, tol);
    };

    let judged = (degree + 1) * (degree + 1);
    loop {
        if current.distances.iter().all(|d| *d <= tolerance) {
            return finish(points, &start, current, Ended::Met, fairing, tol);
        }
        let (cells, strays) = cell_residuals(&current, tolerance);
        let mut missed: Vec<(f64, (usize, usize), usize)> = cells
            .iter()
            .filter_map(|c| c.missed.map(|w| (w, c.at, c.count)))
            .collect();
        if missed.is_empty() && strays {
            return finish(points, &start, current, Ended::Strays, fairing, tol);
        }
        missed.retain(|m| m.2 >= judged);
        missed.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        let room = (
            most.0 - span_count(&current.u_knots),
            most.1 - span_count(&current.v_knots),
        );
        let Some(&(_, (i, j), _)) = missed.first() else {
            return finish(points, &start, current, Ended::Bounds, fairing, tol);
        };
        let split_u: Vec<usize> = if room.0 > 0 { vec![i] } else { Vec::new() };
        let split_v: Vec<usize> = if room.1 > 0 { vec![j] } else { Vec::new() };
        if split_u.is_empty() && split_v.is_empty() {
            return finish(points, &start, current, Ended::Bounds, fairing, tol);
        }
        let mut tries: Vec<(&[usize], &[usize])> = Vec::with_capacity(3);
        if !split_u.is_empty() {
            tries.push((&split_u, &[]));
        }
        if !split_v.is_empty() {
            tries.push((&[], &split_v));
        }
        if !split_u.is_empty() && !split_v.is_empty() {
            tries.push((&split_u, &split_v));
        }
        // Mallows' criterion: a control fitted to noise alone takes the
        // noise's variance off the squared distances on average, so of two
        // fits the one whose squared distances plus twice that variance per
        // control is less predicts unseen points better. The variance is
        // read off the median distance, which the shape the fit has yet to
        // follow moves only where it moves half the points: for a normal
        // noise across the surface the median distance is `Φ⁻¹(3/4)` of
        // its deviation.
        let variance = (median(&current.distances) / NORMAL_QUARTILE).powi(2);
        let score = |r: &Round| 2.0f64.mul_add(variance * precise(r.controls()), r.square_sum());
        let mut best: Option<(f64, Round)> = None;
        for (along_u, along_v) in tries {
            let knots = (
                split(&current.u_knots, along_u)?,
                split(&current.v_knots, along_v)?,
            );
            let tried = round(points, &start, knots, fairing, tol)?;
            let scored = score(&tried);
            if best.as_ref().is_none_or(|(s, _)| scored < *s) {
                best = Some((scored, tried));
            }
        }
        match best {
            Some((scored, next)) if scored < score(&current) => current = next,
            _ => return finish(points, &start, current, Ended::Noise, fairing, tol),
        }
    }
}

/// The result, read off a round fitted at `at`.
fn finish(
    points: &[Point],
    at: &[(f64, f64)],
    round: Round,
    ended: Ended,
    fairing: f64,
    tol: Tolerances,
) -> OgeomResult<ScatteredFit> {
    let cross_sigma = cross_sigma(points, at, &round, fairing, tol)?;
    let spans = (span_count(&round.u_knots), span_count(&round.v_knots));
    let sigma = (round.square_sum() / precise(round.distances.len())).sqrt();
    let mut worst = (0, round.distances[0]);
    for (i, &d) in round.distances.iter().enumerate() {
        if d > worst.1 {
            worst = (i, d);
        }
    }
    Ok(ScatteredFit {
        surface: round.surface,
        parameters: round.feet,
        sigma,
        cross_sigma,
        worst,
        spans,
        ended,
    })
}

/// Each point's parameters on the seed, a periodic direction cut at the
/// points' widest gap round it, scaled onto the unit square.
fn seeded_parameters(
    points: &[Point],
    seed: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<Vec<(f64, f64)>> {
    use ogeom_math::elementary;
    let closed_form = |p: Point| -> Option<OgeomResult<(f64, f64)>> {
        Some(match seed {
            SurfaceGeometry::Plane(s) => Ok(elementary::plane_parameters(&s.plane(), p)),
            SurfaceGeometry::Cylinder(s) => elementary::cylinder_parameters(&s.cylinder(), p, tol),
            SurfaceGeometry::Cone(s) => elementary::cone_parameters(&s.cone(), p, tol),
            SurfaceGeometry::Sphere(s) => elementary::sphere_parameters(&s.sphere(), p, tol),
            SurfaceGeometry::Torus(s) => elementary::torus_parameters(&s.torus(), p, tol),
            _ => return None,
        })
    };
    let mut store: Option<SurfaceSeeds> = None;
    let mut us = Vec::with_capacity(points.len());
    let mut vs = Vec::with_capacity(points.len());
    for &p in points {
        let (u, v) = if let Some(found) = closed_form(p) {
            found?
        } else {
            let store = match &mut store {
                Some(store) => store,
                empty => empty.insert(SurfaceSeeds::over(seed, SEED_SAMPLES, tol)?),
            };
            store.project(seed, p, tol)?.parameters
        };
        us.push(u);
        vs.push(v);
    }
    let ((ua, ub), (va, vb)) = seed.domain();
    let us = unit_range(us, seed.is_periodic_u().then_some(ub - ua), "u")?;
    let vs = unit_range(vs, seed.is_periodic_v().then_some(vb - va), "v")?;
    Ok(us.into_iter().zip(vs).collect())
}

/// `values` onto `[0, 1]`: a periodic direction's first cut at the widest
/// gap between them round the period, so they run in one piece.
fn unit_range(mut values: Vec<f64>, period: Option<f64>, name: &str) -> OgeomResult<Vec<f64>> {
    if let Some(period) = period.filter(|p| *p > 0.0 && p.is_finite()) {
        for t in &mut values {
            *t = t.rem_euclid(period);
        }
        let mut sorted = values.clone();
        sorted.sort_by(f64::total_cmp);
        // The gap round the end of the period first, so an even spread
        // keeps the seed's own seam.
        let mut cut = None;
        let mut widest = sorted[0] + period - sorted[sorted.len() - 1];
        for w in sorted.windows(2) {
            if w[1] - w[0] > widest {
                widest = w[1] - w[0];
                cut = Some(w[1]);
            }
        }
        if let Some(cut) = cut {
            for t in &mut values {
                if *t < cut {
                    *t += period;
                }
            }
        }
    }
    let (lo, hi) = values
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), t| {
            (lo.min(*t), hi.max(*t))
        });
    if !(lo.is_finite() && hi.is_finite()) || hi <= lo {
        ogeom_bail!(
            Construction,
            "the points do not spread along the seed's {name}; the seed does not unfold them"
        );
    }
    let scaled: Vec<f64> = values
        .iter()
        .map(|t| ((t - lo) / (hi - lo)).clamp(0.0, 1.0))
        .collect();
    Ok(scaled)
}

/// The net fitted at `at` over the knots, each point's foot on it found
/// from there.
fn round(
    points: &[Point],
    at: &[(f64, f64)],
    (u_knots, v_knots): (KnotVector, KnotVector),
    fairing: f64,
    tol: Tolerances,
) -> OgeomResult<Round> {
    let surface = fit_surface_scattered_at(points, at, &u_knots, &v_knots, fairing, tol)?.curve;
    let geometry = SurfaceGeometry::BSpline(surface.clone());
    let mut feet = Vec::with_capacity(points.len());
    let mut distances = Vec::with_capacity(points.len());
    for (p, &guess) in points.iter().zip(at) {
        let foot = project_on_surface_from(&geometry, *p, guess, tol)?;
        feet.push(foot.parameters);
        distances.push(foot.distance);
    }
    Ok(Round {
        surface,
        u_knots,
        v_knots,
        feet,
        distances,
    })
}

/// The root mean square distance from each point to the fit of the other
/// half (alternate points by index) at `at`, over the round's knots. A half
/// the net cannot be settled from predicts nothing, and reads as infinitely
/// far.
fn cross_sigma(
    points: &[Point],
    at: &[(f64, f64)],
    round: &Round,
    fairing: f64,
    tol: Tolerances,
) -> OgeomResult<f64> {
    let mut sum = 0.0;
    for held in 0..2 {
        let (mut train, mut train_at) = (Vec::new(), Vec::new());
        for (i, (p, a)) in points.iter().zip(at).enumerate() {
            if i % 2 != held {
                train.push(*p);
                train_at.push(*a);
            }
        }
        let fitted = match fit_surface_scattered_at(
            &train,
            &train_at,
            &round.u_knots,
            &round.v_knots,
            fairing,
            tol,
        ) {
            Ok(f) => f,
            Err(OgeomError::Numeric(_)) => return Ok(f64::INFINITY),
            Err(e) => return Err(e),
        };
        let geometry = SurfaceGeometry::BSpline(fitted.curve);
        for (i, (p, a)) in points.iter().zip(at).enumerate() {
            if i % 2 == held {
                let d = project_on_surface_from(&geometry, *p, *a, tol)?.distance;
                sum += d * d;
            }
        }
    }
    Ok((sum / precise(points.len())).sqrt())
}

/// A cell's points: how many, and the farthest one beyond the tolerance
/// that is not a stray.
struct Cell {
    at: (usize, usize),
    count: usize,
    square_sum: f64,
    missed: Option<f64>,
}

/// The round's points binned by the knot cell their feet lie in, `u` span
/// then `v` span, the empty cells left out, and whether any point beyond
/// the tolerance is a stray.
///
/// A stray stands farther from the fit than noise as wide as the rest of
/// its cell (their root mean square distance) reaches over the whole
/// cloud: the largest of `m` normal deviates is under `√(2 ln m)` of their
/// deviation, and the rest of a cell, its shape misses and all, measures
/// that deviation from above.
fn cell_residuals(round: &Round, tolerance: f64) -> (Vec<Cell>, bool) {
    let (pu, pv) = (round.u_knots.degree(), round.v_knots.degree());
    let (nu, nv) = (span_count(&round.u_knots), span_count(&round.v_knots));
    let mut cells: Vec<Cell> = Vec::with_capacity(nu * nv);
    for i in 0..nu {
        for j in 0..nv {
            cells.push(Cell {
                at: (i, j),
                count: 0,
                square_sum: 0.0,
                missed: None,
            });
        }
    }
    let index: Vec<usize> = round
        .feet
        .iter()
        .map(|&(u, v)| {
            let i = round.u_knots.span_unchecked(u) - pu;
            let j = round.v_knots.span_unchecked(v) - pv;
            i * nv + j
        })
        .collect();
    for (&k, d) in index.iter().zip(&round.distances) {
        cells[k].count += 1;
        cells[k].square_sum += d * d;
    }
    let reach = (2.0 * precise(round.distances.len()).ln()).sqrt();
    let mut strays = false;
    for (&k, &d) in index.iter().zip(&round.distances) {
        if d <= tolerance {
            continue;
        }
        let cell = &mut cells[k];
        let rest = cell.count - 1;
        let stray =
            rest > 0 && d > reach * (((cell.square_sum - d * d).max(0.0)) / precise(rest)).sqrt();
        if stray {
            strays = true;
        } else if cell.missed.is_none_or(|w| d > w) {
            cell.missed = Some(d);
        }
    }
    cells.retain(|c| c.count > 0);
    (cells, strays)
}

/// The median of non-empty `values`: the mean of the middle two of an even
/// count.
fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let k = sorted.len() / 2;
    if sorted.len().is_multiple_of(2) {
        f64::midpoint(sorted[k - 1], sorted[k])
    } else {
        sorted[k]
    }
}

/// How many knot spans a vector's domain has.
fn span_count(knots: &KnotVector) -> usize {
    knots.distinct().len() - 1
}

/// The knot vector with a knot at the middle of each span listed.
fn split(knots: &KnotVector, spans: &[usize]) -> OgeomResult<KnotVector> {
    let distinct = knots.distinct();
    let mut out = knots.clone();
    for &s in spans {
        out = out.with_knot_inserted(f64::midpoint(distinct[s].0, distinct[s + 1].0), 1)?;
    }
    Ok(out)
}
