//! The stationary approaches between two geometries.
//!
//! *Elsewhere* this is `Extrema` and the `GeomAPI_Extrema*` family. The
//! consumer that drives it is minimum distance between shapes, which is
//! `BRepExtrema` there and lives in `ogeom-algo` here — this module answers for
//! the geometry, and the shape layer assembles the answer for topology.
//!
//! # What an extremum is, and what it is not
//!
//! An approach is *stationary* when the connecting vector is perpendicular to
//! every tangent it meets — the derivative of the squared distance is zero in
//! each parameter. Those are the only approaches this module reports. Two
//! kinds of candidate are deliberately not here:
//!
//! - **Domain-end candidates.** A pair of segments whose closest points are
//!   endpoint to endpoint has no interior stationary approach, and the answer
//!   comes back empty. The endpoints are *points*, and points against curves
//!   and surfaces are projections, which the caller owns — at the shape
//!   level, an edge's ends are vertices, and the vertex pairs cover exactly
//!   these candidates. Folding them in here would answer the shape question
//!   badly instead of the geometry question well.
//! - **A guessed point on a constant-distance locus.** Parallel lines,
//!   concentric circles, a sphere inside a sphere: the nearest distance is
//!   attained along a whole locus, and no isolated point is *the* answer.
//!   That is reported as [`Extrema::family`], the way the conventional
//!   kernel's `IsParallel` flag says the same thing.
//!
//! Every approach reported is verifiable on the spot: two parameter sets, the
//! two evaluated points, and the distance between them, which is the claim.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d, Surface, SurfaceGeometry};
use ogeom_math::{Point, Vector, solve};

/// One stationary approach between two geometries.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Approach<A, B> {
    /// The parameters on the first geometry.
    pub on_a: A,
    /// The parameters on the second.
    pub on_b: B,
    /// The evaluated point on the first geometry.
    pub point_a: Point,
    /// The evaluated point on the second.
    pub point_b: Point,
    /// The distance between them — the claim, checkable by evaluation.
    pub distance: f64,
}

/// Every stationary approach found, nearest first.
#[derive(Debug, Clone, PartialEq)]
pub struct Extrema<A, B> {
    /// Stationary approaches, sorted by distance. Nearest approaches first;
    /// stationary *farthest* points, where they exist in the interior, are at
    /// the back of the same list.
    pub approaches: Vec<Approach<A, B>>,
    /// Whether the nearest distance is attained along a locus rather than at
    /// an isolated point — parallel lines, concentric circles, coaxial
    /// quadrics. When set, the nearest approaches are representatives of the
    /// family, not the family.
    pub family: bool,
}

impl<A, B> Extrema<A, B> {
    /// The nearest stationary approach, if any is interior to both domains.
    #[must_use]
    pub fn nearest(&self) -> Option<&Approach<A, B>> {
        self.approaches.first()
    }
}

/// How hard the seeding looks.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExtremaOptions {
    /// How many segments a curve is sampled into.
    pub samples: usize,
    /// How finely a surface is sampled, per direction.
    pub grid: usize,
}

impl Default for ExtremaOptions {
    fn default() -> Self {
        Self {
            samples: 64,
            grid: 24,
        }
    }
}

/// A domain span beyond which sampling answers nothing.
///
/// An untrimmed line or plane spans ±1e9; sixty-four samples across that
/// resolve nothing a caller could want. The analytic line/line case answers
/// without sampling, and everything else is refused with instructions rather
/// than sampled into a wrong answer — the same stance `surface_bounds` takes
/// on an unbounded plane.
const WIDEST_DOMAIN: f64 = 1e8;

/// Two points closer than this many confusions are the same approach.
const DISTINCT: f64 = 1e2;

/// How many seeds a constant-distance configuration is allowed to spawn.
///
/// On a locus every sampled cell ties for local minimum. A few hundred
/// polished representatives are plenty to detect the family and report it;
/// thousands would only repeat them.
const MOST_SEEDS: usize = 256;

// --- curve / curve -----------------------------------------------------------

/// The stationary approaches between two curves.
///
/// Line/line is answered in closed form, parallel included; everything else
/// is seeded from a distance grid over both domains and polished by Newton on
/// the stationarity system.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the options
/// are unusable; [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if a curve's
/// domain is too wide to sample — trim it before asking.
pub fn extrema_curve_curve(
    a: &Curve,
    b: &Curve,
    options: ExtremaOptions,
    tol: Tolerances,
) -> OgeomResult<Extrema<f64, f64>> {
    if options.samples < 2 {
        ogeom_bail!(Construction, "seeding needs at least two samples");
    }
    if let (Curve::Line(la), Curve::Line(lb)) = (a, b) {
        return Ok(line_line(la, lb, tol));
    }
    for (name, curve) in [("first", a), ("second", b)] {
        let (lo, hi) = curve.domain();
        if hi - lo > WIDEST_DOMAIN {
            ogeom_bail!(
                Domain,
                "the {name} curve's domain spans {:.0e}; trim it before asking",
                hi - lo
            );
        }
    }

    let sa = sample_curve(a, options.samples, tol);
    let sb = sample_curve(b, options.samples, tol);
    if sa.len() < 2 || sb.len() < 2 {
        ogeom_bail!(Construction, "a curve failed to evaluate over its domain");
    }

    // Local extrema of the sampled distance field seed the polish. Ties count
    // — on a constant-distance locus everything ties, and the family is
    // exactly what the tied seeds go on to reveal.
    let mut seeds = Vec::new();
    for i in 0..sa.len() {
        for j in 0..sb.len() {
            let here = sa[i].1.square_distance(sb[j].1);
            let mut minimal = true;
            let mut maximal = true;
            let neighbours_a = sa
                .iter()
                .enumerate()
                .take((i + 2).min(sa.len()))
                .skip(i.saturating_sub(1));
            for (ni, near_a) in neighbours_a {
                let neighbours_b = sb
                    .iter()
                    .enumerate()
                    .take((j + 2).min(sb.len()))
                    .skip(j.saturating_sub(1));
                for (nj, near_b) in neighbours_b {
                    if ni == i && nj == j {
                        continue;
                    }
                    let there = near_a.1.square_distance(near_b.1);
                    if there < here {
                        minimal = false;
                    }
                    if there > here {
                        maximal = false;
                    }
                }
            }
            if minimal || maximal {
                seeds.push((sa[i].0, sb[j].0));
            }
        }
    }
    thin(&mut seeds);

    let mut approaches: Vec<Approach<f64, f64>> = Vec::new();
    for (seed_a, seed_b) in seeds {
        if let Some((t, s)) = stationary_curve_curve(a, b, seed_a, seed_b, tol) {
            let (Ok(pa), Ok(pb)) = (a.point_at(t, tol), b.point_at(s, tol)) else {
                continue;
            };
            keep(
                &mut approaches,
                Approach {
                    on_a: t,
                    on_b: s,
                    point_a: pa,
                    point_b: pb,
                    distance: pa.distance(pb),
                },
                tol,
            );
        }
    }
    Ok(finish(approaches, tol))
}

/// Newton on the two stationarity conditions of the squared distance.
pub(crate) fn stationary_curve_curve(
    a: &Curve,
    b: &Curve,
    seed_a: f64,
    seed_b: f64,
    tol: Tolerances,
) -> Option<(f64, f64)> {
    let system = |x: &[f64]| {
        let (t, s) = (fold_curve(a, x[0]), fold_curve(b, x[1]));
        let pa = a.point_at(t, tol).unwrap_or(Point::ORIGIN);
        let pb = b.point_at(s, tol).unwrap_or(Point::ORIGIN);
        let da = a.derivatives_at(t, 2, tol).unwrap_or_default();
        let db = b.derivatives_at(s, 2, tol).unwrap_or_default();
        let zero = Vector::ZERO;
        let (d1a, d2a) = (
            da.get(1).copied().unwrap_or(zero),
            da.get(2).copied().unwrap_or(zero),
        );
        let (d1b, d2b) = (
            db.get(1).copied().unwrap_or(zero),
            db.get(2).copied().unwrap_or(zero),
        );
        let gap = pa - pb;
        (
            vec![gap.dot(d1a), -gap.dot(d1b)],
            vec![
                vec![d1a.dot(d1a) + gap.dot(d2a), -d1a.dot(d1b)],
                vec![-d1a.dot(d1b), d1b.dot(d1b) - gap.dot(d2b)],
            ],
        )
    };
    let criteria = solve::Criteria {
        // The residual is a gradient of a squared distance, not a distance:
        // scale-squared, so the target is too.
        residual: tol.confusion() * tol.confusion(),
        step: tol.parametric(),
        max_iterations: 40,
    };
    let found = solve::newton_system(system, &[seed_a, seed_b], criteria).ok()?;
    Some((fold_curve(a, found.value[0]), fold_curve(b, found.value[1])))
}

/// Closed-form extrema between two lines.
fn line_line(
    a: &ogeom_geom::LineCurve,
    b: &ogeom_geom::LineCurve,
    tol: Tolerances,
) -> Extrema<f64, f64> {
    let (oa, da) = (a.axis().location, a.axis().direction.vector());
    let (ob, db) = (b.axis().location, b.axis().direction.vector());
    let cross = da.cross(db);
    let denominator = cross.dot(cross);

    if denominator <= tol.angular() * tol.angular() {
        // Parallel: constant distance, a family. The representative sits at
        // the middle of where the domains face each other; if they face each
        // other nowhere, the nearest is at ends the caller owns.
        let (a_lo, a_hi) = a.domain();
        let (b_lo, b_hi) = b.domain();
        let project = |p: Point| (p - oa).dot(da);
        let (s0, s1) = (project(ob + db * b_lo), project(ob + db * b_hi));
        let (lo, hi) = (s0.min(s1).max(a_lo), s0.max(s1).min(a_hi));
        if lo > hi {
            return Extrema {
                approaches: Vec::new(),
                family: true,
            };
        }
        let t = f64::midpoint(lo, hi);
        let pa = oa + da * t;
        let s = (pa - ob).dot(db);
        let pb = ob + db * s;
        return Extrema {
            approaches: vec![Approach {
                on_a: t,
                on_b: s,
                point_a: pa,
                point_b: pb,
                distance: pa.distance(pb),
            }],
            family: true,
        };
    }

    let between = ob - oa;
    let t = between.cross(db).dot(cross) / denominator;
    let s = between.cross(da).dot(cross) / denominator;
    let (a_lo, a_hi) = a.domain();
    let (b_lo, b_hi) = b.domain();
    if t < a_lo || t > a_hi || s < b_lo || s > b_hi {
        // The stationary approach exists on the unbounded lines but outside
        // these domains: within them the nearest is at an end.
        return Extrema {
            approaches: Vec::new(),
            family: false,
        };
    }
    let pa = oa + da * t;
    let pb = ob + db * s;
    Extrema {
        approaches: vec![Approach {
            on_a: t,
            on_b: s,
            point_a: pa,
            point_b: pb,
            distance: pa.distance(pb),
        }],
        family: false,
    }
}

// --- curve / surface ---------------------------------------------------------

/// The stationary approaches between a curve and a surface.
///
/// # Errors
///
/// As [`extrema_curve_curve`].
pub fn extrema_curve_surface(
    curve: &Curve,
    surface: &SurfaceGeometry,
    options: ExtremaOptions,
    tol: Tolerances,
) -> OgeomResult<Extrema<f64, (f64, f64)>> {
    if options.samples < 2 || options.grid < 2 {
        ogeom_bail!(Construction, "seeding needs at least two steps each way");
    }
    let (lo, hi) = curve.domain();
    if hi - lo > WIDEST_DOMAIN {
        ogeom_bail!(
            Domain,
            "the curve's domain spans {:.0e}; trim it before asking",
            hi - lo
        );
    }
    wide_surface_check(surface)?;

    let sc = sample_curve(curve, options.samples, tol);
    let ss = sample_surface(surface, options.grid, tol);
    if sc.len() < 2 || ss.is_empty() {
        ogeom_bail!(
            Construction,
            "a geometry failed to evaluate over its domain"
        );
    }

    // For each curve sample, its best and worst surface cells; local extrema
    // along the curve of those fields seed the polish.
    let mut best = Vec::with_capacity(sc.len());
    let mut worst = Vec::with_capacity(sc.len());
    for (_, p) in &sc {
        let mut near = (f64::INFINITY, (0.0, 0.0));
        let mut far = (f64::NEG_INFINITY, (0.0, 0.0));
        for (uv, q) in &ss {
            let d = p.square_distance(*q);
            if d < near.0 {
                near = (d, *uv);
            }
            if d > far.0 {
                far = (d, *uv);
            }
        }
        best.push(near);
        worst.push(far);
    }
    let mut seeds = Vec::new();
    for i in 0..sc.len() {
        let lower = i == 0 || best[i].0 <= best[i - 1].0;
        let upper = i + 1 == sc.len() || best[i].0 <= best[i + 1].0;
        if lower && upper {
            seeds.push((sc[i].0, best[i].1));
        }
        let lower = i == 0 || worst[i].0 >= worst[i - 1].0;
        let upper = i + 1 == sc.len() || worst[i].0 >= worst[i + 1].0;
        if lower && upper {
            seeds.push((sc[i].0, worst[i].1));
        }
    }
    thin(&mut seeds);

    let mut approaches: Vec<Approach<f64, (f64, f64)>> = Vec::new();
    for (seed_t, seed_uv) in seeds {
        if let Some((t, u, v)) = stationary_curve_surface(curve, surface, seed_t, seed_uv, tol) {
            let (Ok(pc), Ok(ps)) = (curve.point_at(t, tol), surface.point_at(u, v, tol)) else {
                continue;
            };
            keep(
                &mut approaches,
                Approach {
                    on_a: t,
                    on_b: (u, v),
                    point_a: pc,
                    point_b: ps,
                    distance: pc.distance(ps),
                },
                tol,
            );
        }
    }
    Ok(finish(approaches, tol))
}

/// Newton on the three stationarity conditions.
fn stationary_curve_surface(
    curve: &Curve,
    surface: &SurfaceGeometry,
    seed_t: f64,
    seed_uv: (f64, f64),
    tol: Tolerances,
) -> Option<(f64, f64, f64)> {
    let system = |x: &[f64]| {
        let t = fold_curve(curve, x[0]);
        let (u, v) = fold_surface(surface, x[1], x[2]);
        let pc = curve.point_at(t, tol).unwrap_or(Point::ORIGIN);
        let ps = surface.point_at(u, v, tol).unwrap_or(Point::ORIGIN);
        let dc = curve.derivatives_at(t, 2, tol).unwrap_or_default();
        let zero = Vector::ZERO;
        let (ct, ctt) = (
            dc.get(1).copied().unwrap_or(zero),
            dc.get(2).copied().unwrap_or(zero),
        );
        let (su, sv) = surface.d1_at(u, v, tol).unwrap_or((zero, zero));
        let (suu, suv, svv) = surface.d2_at(u, v, tol).unwrap_or((zero, zero, zero));
        let gap = pc - ps;
        (
            vec![gap.dot(ct), gap.dot(su), gap.dot(sv)],
            vec![
                vec![ct.dot(ct) + gap.dot(ctt), -su.dot(ct), -sv.dot(ct)],
                vec![
                    ct.dot(su),
                    -su.dot(su) + gap.dot(suu),
                    -sv.dot(su) + gap.dot(suv),
                ],
                vec![
                    ct.dot(sv),
                    -su.dot(sv) + gap.dot(suv),
                    -sv.dot(sv) + gap.dot(svv),
                ],
            ],
        )
    };
    let criteria = solve::Criteria {
        residual: tol.confusion() * tol.confusion(),
        step: tol.parametric(),
        max_iterations: 40,
    };
    let found = solve::newton_system(system, &[seed_t, seed_uv.0, seed_uv.1], criteria).ok()?;
    let t = fold_curve(curve, found.value[0]);
    let (u, v) = fold_surface(surface, found.value[1], found.value[2]);
    Some((t, u, v))
}

// --- surface / surface -------------------------------------------------------

/// The stationary approaches between two surfaces.
///
/// # Errors
///
/// As [`extrema_curve_curve`].
pub fn extrema_surface_surface(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: ExtremaOptions,
    tol: Tolerances,
) -> OgeomResult<Extrema<(f64, f64), (f64, f64)>> {
    if options.grid < 2 {
        ogeom_bail!(Construction, "seeding needs at least two steps each way");
    }
    wide_surface_check(a)?;
    wide_surface_check(b)?;

    let sa = sample_surface(a, options.grid, tol);
    let sb = sample_surface(b, options.grid, tol);
    if sa.is_empty() || sb.is_empty() {
        ogeom_bail!(Construction, "a surface failed to evaluate over its domain");
    }

    // Each side seeds from its own best-facing view of the other; the union
    // covers approaches either sampling would resolve.
    let mut seeds = Vec::new();
    for (uv_a, p) in &sa {
        let mut near = (f64::INFINITY, (0.0, 0.0));
        let mut far = (f64::NEG_INFINITY, (0.0, 0.0));
        for (uv_b, q) in &sb {
            let d = p.square_distance(*q);
            if d < near.0 {
                near = (d, *uv_b);
            }
            if d > far.0 {
                far = (d, *uv_b);
            }
        }
        seeds.push((*uv_a, near.1));
        seeds.push((*uv_a, far.1));
    }
    for (uv_b, q) in &sb {
        let mut near = (f64::INFINITY, (0.0, 0.0));
        for (uv_a, p) in &sa {
            let d = q.square_distance(*p);
            if d < near.0 {
                near = (d, *uv_a);
            }
        }
        seeds.push((near.1, *uv_b));
    }
    thin(&mut seeds);

    let mut approaches: Vec<Approach<(f64, f64), (f64, f64)>> = Vec::new();
    for (seed_a, seed_b) in seeds {
        if let Some((ua, va, ub, vb)) = stationary_surface_surface(a, b, seed_a, seed_b, tol) {
            let (Ok(pa), Ok(pb)) = (a.point_at(ua, va, tol), b.point_at(ub, vb, tol)) else {
                continue;
            };
            keep(
                &mut approaches,
                Approach {
                    on_a: (ua, va),
                    on_b: (ub, vb),
                    point_a: pa,
                    point_b: pb,
                    distance: pa.distance(pb),
                },
                tol,
            );
        }
    }
    Ok(finish(approaches, tol))
}

/// Newton on the four stationarity conditions.
fn stationary_surface_surface(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    seed_a: (f64, f64),
    seed_b: (f64, f64),
    tol: Tolerances,
) -> Option<(f64, f64, f64, f64)> {
    let system = |x: &[f64]| {
        let (ua, va) = fold_surface(a, x[0], x[1]);
        let (ub, vb) = fold_surface(b, x[2], x[3]);
        let zero = Vector::ZERO;
        let pa = a.point_at(ua, va, tol).unwrap_or(Point::ORIGIN);
        let pb = b.point_at(ub, vb, tol).unwrap_or(Point::ORIGIN);
        let (au, av) = a.d1_at(ua, va, tol).unwrap_or((zero, zero));
        let (auu, auv, avv) = a.d2_at(ua, va, tol).unwrap_or((zero, zero, zero));
        let (bu, bv) = b.d1_at(ub, vb, tol).unwrap_or((zero, zero));
        let (buu, buv, bvv) = b.d2_at(ub, vb, tol).unwrap_or((zero, zero, zero));
        let gap = pa - pb;
        (
            vec![gap.dot(au), gap.dot(av), gap.dot(bu), gap.dot(bv)],
            vec![
                vec![
                    au.dot(au) + gap.dot(auu),
                    au.dot(av) + gap.dot(auv),
                    -bu.dot(au),
                    -bv.dot(au),
                ],
                vec![
                    au.dot(av) + gap.dot(auv),
                    av.dot(av) + gap.dot(avv),
                    -bu.dot(av),
                    -bv.dot(av),
                ],
                vec![
                    au.dot(bu),
                    av.dot(bu),
                    -bu.dot(bu) + gap.dot(buu),
                    -bv.dot(bu) + gap.dot(buv),
                ],
                vec![
                    au.dot(bv),
                    av.dot(bv),
                    -bu.dot(bv) + gap.dot(buv),
                    -bv.dot(bv) + gap.dot(bvv),
                ],
            ],
        )
    };
    let criteria = solve::Criteria {
        residual: tol.confusion() * tol.confusion(),
        step: tol.parametric(),
        max_iterations: 40,
    };
    let found =
        solve::newton_system(system, &[seed_a.0, seed_a.1, seed_b.0, seed_b.1], criteria).ok()?;
    let (ua, va) = fold_surface(a, found.value[0], found.value[1]);
    let (ub, vb) = fold_surface(b, found.value[2], found.value[3]);
    Some((ua, va, ub, vb))
}

// --- shared machinery --------------------------------------------------------

fn wide_surface_check(surface: &SurfaceGeometry) -> OgeomResult<()> {
    let ((ua, ub), (va, vb)) = surface.domain();
    if ub - ua > WIDEST_DOMAIN || vb - va > WIDEST_DOMAIN {
        ogeom_bail!(
            Domain,
            "a surface domain spans more than {WIDEST_DOMAIN:.0e}; trim it before asking"
        );
    }
    Ok(())
}

fn sample_curve(curve: &Curve, samples: usize, tol: Tolerances) -> Vec<(f64, Point)> {
    let (lo, hi) = curve.domain();
    let mut out = Vec::with_capacity(samples + 1);
    for i in 0..=samples {
        #[allow(clippy::cast_precision_loss)]
        let t = lo + (hi - lo) * i as f64 / samples as f64;
        if let Ok(p) = curve.point_at(t, tol) {
            out.push((t, p));
        }
    }
    out
}

#[allow(clippy::type_complexity)]
fn sample_surface(
    surface: &SurfaceGeometry,
    grid: usize,
    tol: Tolerances,
) -> Vec<((f64, f64), Point)> {
    let ((ua, ub), (va, vb)) = surface.domain();
    let mut out = Vec::with_capacity((grid + 1) * (grid + 1));
    for i in 0..=grid {
        for j in 0..=grid {
            #[allow(clippy::cast_precision_loss)]
            let u = ua + (ub - ua) * i as f64 / grid as f64;
            #[allow(clippy::cast_precision_loss)]
            let v = va + (vb - va) * j as f64 / grid as f64;
            if let Ok(p) = surface.point_at(u, v, tol) {
                out.push(((u, v), p));
            }
        }
    }
    out
}

/// Cap the seed list, keeping an even spread.
fn thin<T>(seeds: &mut Vec<T>) {
    if seeds.len() <= MOST_SEEDS {
        return;
    }
    let step = seeds.len().div_ceil(MOST_SEEDS);
    let mut index = 0;
    seeds.retain(|_| {
        let kept = index % step == 0;
        index += 1;
        kept
    });
}

/// Add an approach unless one at the same pair of points is already known.
fn keep<A: Copy, B: Copy>(
    approaches: &mut Vec<Approach<A, B>>,
    candidate: Approach<A, B>,
    tol: Tolerances,
) {
    let reach = tol.confusion() * DISTINCT;
    if approaches.iter().any(|known| {
        known.point_a.distance(candidate.point_a) <= reach
            && known.point_b.distance(candidate.point_b) <= reach
    }) {
        return;
    }
    approaches.push(candidate);
}

/// Sort by distance and decide whether the nearest is a family.
fn finish<A: Copy, B: Copy>(mut approaches: Vec<Approach<A, B>>, tol: Tolerances) -> Extrema<A, B> {
    approaches.sort_by(|a, b| {
        a.distance
            .partial_cmp(&b.distance)
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    let family = match approaches.first() {
        None => false,
        Some(first) => {
            let near = tol.confusion().max(first.distance * 1e-9);
            let ties: Vec<&Approach<A, B>> = approaches
                .iter()
                .take_while(|a| a.distance - first.distance <= near)
                .collect();
            // Three or more equally-near approaches at genuinely different
            // places are not coincidence; they are a locus showing through
            // the sampling.
            ties.len() >= 3
                && ties
                    .iter()
                    .any(|a| a.point_a.distance(first.point_a) > tol.confusion() * DISTINCT * 10.0)
        }
    };
    Extrema { approaches, family }
}

fn fold_curve(curve: &Curve, t: f64) -> f64 {
    let (lo, hi) = curve.domain();
    if curve.is_periodic() {
        let span = hi - lo;
        if span > 0.0 {
            return lo + (t - lo).rem_euclid(span);
        }
    }
    t.clamp(lo, hi)
}

fn fold_surface(surface: &SurfaceGeometry, u: f64, v: f64) -> (f64, f64) {
    let ((ua, ub), (va, vb)) = surface.domain();
    let fold = |x: f64, lo: f64, hi: f64, periodic: bool| {
        if periodic {
            let span = hi - lo;
            if span > 0.0 {
                return lo + (x - lo).rem_euclid(span);
            }
        }
        x.clamp(lo, hi)
    };
    (
        fold(u, ua, ub, surface.is_periodic_u()),
        fold(v, va, vb, surface.is_periodic_v()),
    )
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use ogeom_geom::{CircleCurve, CylinderSurface, LineCurve, PlaneSurface, SphereSurface};
    use ogeom_math::{Circle, Cylinder, Direction, Frame, Plane, Sphere};

    const T: Tolerances = Tolerances::millimetres();

    fn segment(from: Point, to: Point) -> Curve {
        LineCurve::segment(from, to, T).unwrap().into()
    }

    fn circle_at(centre: Point, normal: Vector, radius: f64) -> Curve {
        CircleCurve::new(
            Circle::new(
                Frame::new(
                    centre,
                    Direction::new(normal, T).unwrap(),
                    Direction::from_cross(normal, Vector::new(0.3, 0.5, 0.9), T).unwrap(),
                    T,
                )
                .unwrap(),
                radius,
                T,
            )
            .unwrap(),
        )
        .into()
    }

    fn sphere_at(centre: Point, radius: f64) -> SurfaceGeometry {
        SphereSurface::new(Sphere::centred(centre, radius, T).unwrap()).into()
    }

    #[test]
    fn skew_segments_meet_the_closed_form() {
        // The x axis, and a line along y lifted by one: nearest distance one,
        // at the origin and at (0, 0, 1).
        let a = segment(Point::new(-5.0, 0.0, 0.0), Point::new(5.0, 0.0, 0.0));
        let b = segment(Point::new(0.0, -5.0, 1.0), Point::new(0.0, 5.0, 1.0));
        let found = extrema_curve_curve(&a, &b, ExtremaOptions::default(), T).unwrap();
        let nearest = found.nearest().unwrap();
        assert!((nearest.distance - 1.0).abs() < 1e-9);
        assert!(nearest.point_a.is_equal(Point::ORIGIN, T));
        assert!(nearest.point_b.is_equal(Point::new(0.0, 0.0, 1.0), T));
        assert!(!found.family);
    }

    #[test]
    fn endpoint_to_endpoint_nearness_is_the_callers_and_says_so() {
        // Collinear segments end to end: no interior stationary approach
        // exists, and pretending one did would misreport the geometry. The
        // endpoints are vertices at the shape level, and that is where this
        // answer lives.
        let a = segment(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0));
        let b = segment(Point::new(3.0, 0.0, 0.0), Point::new(5.0, 0.0, 0.0));
        let found = extrema_curve_curve(&a, &b, ExtremaOptions::default(), T).unwrap();
        assert!(found.approaches.is_empty());
    }

    #[test]
    fn parallel_lines_are_a_family_with_a_representative() {
        let a = segment(Point::new(-4.0, 0.0, 0.0), Point::new(4.0, 0.0, 0.0));
        let b = segment(Point::new(-2.0, 2.0, 0.0), Point::new(6.0, 2.0, 0.0));
        let found = extrema_curve_curve(&a, &b, ExtremaOptions::default(), T).unwrap();
        assert!(found.family);
        let nearest = found.nearest().unwrap();
        assert!((nearest.distance - 2.0).abs() < 1e-12);
        // The representative sits where the domains face each other.
        assert!(nearest.on_a >= -2.0 && nearest.on_a <= 8.0);
    }

    #[test]
    fn concentric_circles_are_a_family_found_by_sampling() {
        // No closed form handles this pair; the family shows through the
        // seeded path as many equally-near approaches at different places.
        let a = circle_at(Point::ORIGIN, Vector::Z, 3.0);
        let b = circle_at(Point::ORIGIN, Vector::Z, 1.0);
        let found = extrema_curve_curve(&a, &b, ExtremaOptions::default(), T).unwrap();
        assert!(found.family);
        assert!((found.nearest().unwrap().distance - 2.0).abs() < 1e-9);
    }

    #[test]
    fn a_tilted_circle_over_a_circle_has_isolated_extrema() {
        // Tilt one circle: the family collapses to isolated nearest and
        // farthest approaches.
        let a = circle_at(Point::new(0.0, 0.0, 2.0), Vector::new(0.3, 0.0, 1.0), 3.0);
        let b = circle_at(Point::ORIGIN, Vector::Z, 3.0);
        let found = extrema_curve_curve(&a, &b, ExtremaOptions::default(), T).unwrap();
        assert!(!found.family);
        let nearest = found.nearest().unwrap();
        // Verifiable on the spot: the claim is the distance between the
        // evaluated points.
        assert!((nearest.point_a.distance(nearest.point_b) - nearest.distance).abs() < 1e-12);
        assert!(nearest.distance < 2.0, "the tilt brings the rims closer");
    }

    #[test]
    fn a_segment_passing_a_sphere_finds_the_gap_to_it() {
        // A line at distance three from the centre of a unit sphere: nearest
        // approach two, on the line of shortest connection.
        let line = segment(Point::new(-5.0, 3.0, 0.0), Point::new(5.0, 3.0, 0.0));
        let ball = sphere_at(Point::ORIGIN, 1.0);
        let found = extrema_curve_surface(&line, &ball, ExtremaOptions::default(), T).unwrap();
        let nearest = found.nearest().unwrap();
        assert!((nearest.distance - 2.0).abs() < 1e-9);
        assert!(nearest.point_a.is_equal(Point::new(0.0, 3.0, 0.0), T));
        assert!(nearest.point_b.is_equal(Point::new(0.0, 1.0, 0.0), T));
    }

    #[test]
    fn a_circle_parallel_to_a_plane_is_a_family_above_it() {
        let ring = circle_at(Point::new(0.0, 0.0, 2.0), Vector::Z, 3.0);
        let ground: SurfaceGeometry = PlaneSurface::over(
            Plane::through(Point::ORIGIN, Direction::Z),
            (-8.0, 8.0),
            (-8.0, 8.0),
        )
        .unwrap()
        .into();
        let found = extrema_curve_surface(&ring, &ground, ExtremaOptions::default(), T).unwrap();
        assert!(found.family);
        assert!((found.nearest().unwrap().distance - 2.0).abs() < 1e-9);
    }

    #[test]
    fn two_spheres_apart_meet_along_the_line_of_centres() {
        let a = sphere_at(Point::ORIGIN, 1.0);
        let b = sphere_at(Point::new(5.0, 0.0, 0.0), 2.0);
        let found = extrema_surface_surface(&a, &b, ExtremaOptions::default(), T).unwrap();
        let nearest = found.nearest().unwrap();
        assert!((nearest.distance - 2.0).abs() < 1e-9);
        assert!(nearest.point_a.is_equal(Point::new(1.0, 0.0, 0.0), T));
        assert!(nearest.point_b.is_equal(Point::new(3.0, 0.0, 0.0), T));
        assert!(!found.family);
    }

    #[test]
    fn concentric_spheres_are_a_family() {
        let a = sphere_at(Point::ORIGIN, 1.0);
        let b = sphere_at(Point::ORIGIN, 3.0);
        let found = extrema_surface_surface(&a, &b, ExtremaOptions::default(), T).unwrap();
        assert!(found.family);
        assert!((found.nearest().unwrap().distance - 2.0).abs() < 1e-9);
    }

    #[test]
    fn a_cylinder_beside_a_plane_reports_the_ruling_gap_as_a_family() {
        // The nearest locus is a whole ruling of the cylinder.
        let drum: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 1.0, T).unwrap(), (-3.0, 3.0))
                .unwrap()
                .into();
        let wall: SurfaceGeometry = PlaneSurface::over(
            Plane::through(Point::new(4.0, 0.0, 0.0), Direction::X),
            (-8.0, 8.0),
            (-8.0, 8.0),
        )
        .unwrap()
        .into();
        let found = extrema_surface_surface(&drum, &wall, ExtremaOptions::default(), T).unwrap();
        assert!(found.family);
        assert!((found.nearest().unwrap().distance - 3.0).abs() < 1e-9);
    }

    #[test]
    fn an_untrimmed_line_is_refused_with_instructions() {
        let endless: Curve = LineCurve::new(ogeom_math::Axis {
            location: Point::ORIGIN,
            direction: Direction::X,
        })
        .into();
        let ring = circle_at(Point::ORIGIN, Vector::Z, 1.0);
        assert!(extrema_curve_curve(&endless, &ring, ExtremaOptions::default(), T).is_err());
        // Line against line has its closed form and needs no trimming.
        let other: Curve = LineCurve::new(ogeom_math::Axis {
            location: Point::new(0.0, 1.0, 0.0),
            direction: Direction::Y,
        })
        .into();
        assert!(extrema_curve_curve(&endless, &other, ExtremaOptions::default(), T).is_ok());
    }
}
