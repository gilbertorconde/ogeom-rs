//! Lofts that follow guide curves.
//!
//! A guided loft is a Gordon surface. The skin across the sections (each
//! section's control points interpolated across, as an unguided loft
//! does) is corrected, for each guide, by the guide's departure from that
//! skin's interpolation of the guide's own crossing points, the
//! corrections blended along the sections by the functions that
//! interpolate across the guides. For that, every section is paced so each
//! guide crosses it at one common parameter and every guide so it crosses
//! each section at one common parameter, both by a piecewise linear change
//! of parameter, which leaves every curve where it is.
//!
//! The sum is taken in homogeneous coordinates. Each guide is first given
//! the weight the sections carry where it crosses them, by multiplying it
//! through with the scalar spline that interpolates the ratio of the
//! weights (which leaves the guide where it is). At a crossing the guide's
//! term and the section's are then one homogeneous point, the corrections
//! vanish on every section, and the skin passes through every section
//! exactly and along every guide to within the distance by which the guide
//! misses the sections. Both are measured on the skin.

use ogeom_algo::{Built, History};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::BSplineCurve;
use ogeom_geom::Curve3d as _;
use ogeom_geom::Surface as _;
use ogeom_math::bspline::{self, Spline};
use ogeom_math::{Blend, KnotVector, Point, Weighted};
use ogeom_topo::{Model, Shape};

use super::{KNOT_SAME, Section, inverted, read_section, sheet, standard, surface_of};

/// Parameters closer than this to each other, or to an end of the unit
/// domain, are one parameter.
const ONE_PARAMETER: f64 = 1e-8;

/// Points per curve at which the skin is measured against its sections
/// and guides.
const MEASURED: u32 = 256;

/// A homogeneous spline: a rational curve's control points multiplied
/// through by their weights.
type Homogeneous = Spline<Weighted<Point>>;

/// The guided loft through `sections` following `guides`; see
/// [`make_loft_surface`](super::make_loft_surface).
#[allow(clippy::too_many_lines, reason = "one construction, spelled out")]
pub(super) fn guided_loft(
    model: &mut Model,
    sections: &[Shape],
    closed: bool,
    guides: &[Shape],
    ruled: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if ruled {
        ogeom_bail!(
            Construction,
            "a ruled loft surface is straight between its sections and follows no guide curves"
        );
    }
    if closed {
        ogeom_bail!(
            Construction,
            "a closed loft surface does not follow guide curves"
        );
    }
    if sections.len() < 2 {
        ogeom_bail!(
            Construction,
            "a loft surface needs at least 2 sections, given {}",
            sections.len()
        );
    }
    let mut read: Vec<Section> = Vec::with_capacity(sections.len());
    for (k, shape) in sections.iter().enumerate() {
        let section = read_section(model, shape, "loft section", tol)?;
        if section.edges.len() != 1 {
            ogeom_bail!(
                Construction,
                "a guided loft surface takes sections of one edge each; section {k} has {}",
                section.edges.len()
            );
        }
        read.push(section);
    }
    let closed_u = read[0].closed;
    for (k, s) in read.iter().enumerate() {
        if s.closed != closed_u {
            ogeom_bail!(
                Construction,
                "loft sections are all closed or all open; section {k} is {} and section 0 \
                 is not",
                if s.closed { "closed" } else { "open" }
            );
        }
    }
    let mut rails: Vec<BSplineCurve> = Vec::with_capacity(guides.len());
    for shape in guides {
        rails.push(joined(
            &read_section(model, shape, "loft guide", tol)?,
            tol,
        )?);
    }
    let near = tol.confusion() * 10.0;
    let (m, n) = (rails.len(), read.len());

    // Closed sections are cut open where the first guide crosses them.
    if closed_u {
        for s in &mut read {
            let (u, _, off) = crossing(&s.edges[0].curve, &rails[0], tol)?;
            if off <= near && u > ONE_PARAMETER && u < 1.0 - ONE_PARAMETER {
                s.edges[0].curve = standard(&s.edges[0].curve.reseamed_at(u, tol)?)?;
                s.edges[0].edge = None;
            }
        }
    }

    // Where each guide crosses each section: `(u, t)`, the section's
    // parameter and the guide's.
    let mut cross = vec![vec![(0.0, 0.0); n]; m];
    for (j, rail) in rails.iter().enumerate() {
        for (k, s) in read.iter().enumerate() {
            let (u, t, off) = crossing(&s.edges[0].curve, rail, tol)?;
            if off > near {
                ogeom_bail!(
                    Construction,
                    "guide {j} misses section {k} by {off:.3e}; a guide must cross every \
                     section within {near:.1e}"
                );
            }
            let u = if closed_u && u >= 1.0 - ONE_PARAMETER {
                0.0
            } else {
                at_ends(u)
            };
            cross[j][k] = (u, at_ends(t));
        }
    }
    for j in 0..m {
        if cross[j][n - 1].1 < cross[j][0].1 {
            rails[j] = reversed(&rails[j])?;
            for c in &mut cross[j] {
                c.1 = 1.0 - c.1;
            }
        }
        if cross[j]
            .windows(2)
            .any(|w| w[1].1 <= w[0].1 + ONE_PARAMETER)
        {
            ogeom_bail!(
                Construction,
                "guide {j} does not cross the sections in their order"
            );
        }
        let at_start = cross[j].iter().filter(|c| c.0 == 0.0).count();
        let at_end = cross[j].iter().filter(|c| c.0 == 1.0).count();
        if (at_start != 0 && at_start != n) || (at_end != 0 && at_end != n) {
            ogeom_bail!(
                Construction,
                "guide {j} crosses some sections at an end and others inside them"
            );
        }
    }
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&a, &b| cross[a][0].0.total_cmp(&cross[b][0].0));
    for w in order.windows(2) {
        let (before, after) = (&cross[w[0]], &cross[w[1]]);
        if let Some(k) = (0..n).find(|&k| after[k].0 <= before[k].0 + ONE_PARAMETER) {
            ogeom_bail!(
                Construction,
                "guides {} and {} cross section {k} in another order than section 0, or at one \
                 point",
                w[0],
                w[1]
            );
        }
    }

    // The common parameters: each guide's mean along the sections, each
    // section's mean along the guides (each guide's run between the first
    // section and the last taken as the unit).
    #[allow(clippy::cast_precision_loss)]
    let (sections_count, guides_count) = (n as f64, m as f64);
    let mut u_at: Vec<f64> = order
        .iter()
        .map(|&j| {
            let c0 = cross[j][0].0;
            if c0 == 0.0 || c0 == 1.0 {
                c0
            } else {
                cross[j].iter().map(|c| c.0).sum::<f64>() / sections_count
            }
        })
        .collect();
    let mut v_at: Vec<f64> = (0..n)
        .map(|k| {
            cross
                .iter()
                .map(|c| (c[k].1 - c[0].1) / (c[n - 1].1 - c[0].1))
                .sum::<f64>()
                / guides_count
        })
        .collect();
    v_at[0] = 0.0;
    v_at[n - 1] = 1.0;

    // Sections paced so guide `order[i]` crosses each at `u_at[i]`.
    for (k, s) in read.iter_mut().enumerate() {
        let mut pairs = vec![(0.0, 0.0)];
        for (i, &j) in order.iter().enumerate() {
            let u = cross[j][k].0;
            if u > 0.0 && u < 1.0 {
                pairs.push((u, u_at[i]));
            }
        }
        pairs.push((1.0, 1.0));
        let edge = &mut s.edges[0];
        if pairs.iter().any(|(a, b)| (a - b).abs() > KNOT_SAME) {
            edge.curve = repaced(&edge.curve, &pairs, tol)?;
            edge.paced = false;
        }
    }
    // Guides paced so each crosses section `k` at `v_at[k]`, trimmed to
    // the run between the first section and the last.
    let mut paced_rails: Vec<BSplineCurve> = Vec::with_capacity(m);
    for &j in &order {
        let pairs: Vec<(f64, f64)> = cross[j].iter().zip(&v_at).map(|(c, v)| (c.1, *v)).collect();
        paced_rails.push(repaced(&rails[j], &pairs, tol)?);
    }
    // Through closed sections the seam's guide stands at both ends.
    let mut names: Vec<usize> = order.clone();
    if closed_u {
        u_at.push(1.0);
        paced_rails.push(paced_rails[0].clone());
        names.push(order[0]);
    }

    let (u_knots, u_cardinals) = cardinals(&u_at, tol)?;
    let (v_knots, v_cardinals) = cardinals(&v_at, tol)?;
    let curves: Vec<Homogeneous> = read
        .iter()
        .map(|s| {
            let c = &s.edges[0].curve;
            (c.knots().clone(), c.control_points().to_vec())
        })
        .collect();

    // Each guide weighted as the sections are where it crosses them, and
    // its homogeneous points at the crossings.
    let mut weighted: Vec<Homogeneous> = Vec::with_capacity(paced_rails.len());
    let mut corners: Vec<Vec<Weighted<Point>>> = Vec::with_capacity(paced_rails.len());
    for (rail, &u) in paced_rails.iter().zip(&u_at) {
        let own: Homogeneous = (rail.knots().clone(), rail.control_points().to_vec());
        let mut ratios = Vec::with_capacity(n);
        let mut at = Vec::with_capacity(n);
        for (curve, &v) in curves.iter().zip(&v_at) {
            let on_section = bspline::evaluate(&curve.0, &curve.1, u, tol)?;
            let on_guide = bspline::evaluate(&own.0, &own.1, v, tol)?;
            let ratio = on_section.weight / on_guide.weight;
            ratios.push(ratio);
            at.push(on_guide.scale(ratio));
        }
        let first = ratios[0];
        let guide = if ratios
            .iter()
            .all(|r| (r - first).abs() <= 1e-12 * first.abs())
        {
            (own.0, own.1.iter().map(|p| p.scale(first)).collect())
        } else {
            let width = v_cardinals[0].len();
            let scale: Vec<f64> = (0..width)
                .map(|l| v_cardinals.iter().zip(&ratios).map(|(c, r)| c[l] * r).sum())
                .collect();
            product(&(v_knots.clone(), scale), &own, tol)?
        };
        weighted.push(guide);
        corners.push(at);
    }

    // Everything in one space along the sections and one across them.
    let along: Vec<&KnotVector> = curves
        .iter()
        .map(|c| &c.0)
        .chain(core::iter::once(&u_knots))
        .collect();
    let (u_degree, u_inner) = common_space(&along);
    let across: Vec<&KnotVector> = weighted
        .iter()
        .map(|g| &g.0)
        .chain(core::iter::once(&v_knots))
        .collect();
    let (v_degree, v_inner) = common_space(&across);
    let mut c = Vec::with_capacity(n);
    for curve in &curves {
        c.push(into_space(curve, u_degree, &u_inner, tol)?);
    }
    let mut mu = Vec::with_capacity(u_cardinals.len());
    for f in &u_cardinals {
        mu.push(into_space(&(u_knots.clone(), f.clone()), u_degree, &u_inner, tol)?.1);
    }
    let mut lv = Vec::with_capacity(n);
    for f in &v_cardinals {
        lv.push(into_space(&(v_knots.clone(), f.clone()), v_degree, &v_inner, tol)?.1);
    }
    let mut g = Vec::with_capacity(weighted.len());
    for guide in &weighted {
        g.push(into_space(guide, v_degree, &v_inner, tol)?.1);
    }
    let (net_u, net_v) = (c[0].1.len(), lv[0].len());
    // Each guide's departure from the skin's interpolation of its corners.
    let departures: Vec<Vec<Weighted<Point>>> = g
        .iter()
        .zip(&corners)
        .map(|(guide, at)| {
            (0..net_v)
                .map(|l| {
                    let across = at
                        .iter()
                        .zip(&lv)
                        .fold(Weighted::<Point>::zero(), |acc, (p, f)| {
                            acc.add(p.scale(f[l]))
                        });
                    guide[l].sub(across)
                })
                .collect()
        })
        .collect();
    let rows: Vec<Vec<Weighted<Point>>> = (0..net_v)
        .map(|l| {
            (0..net_u)
                .map(|i| {
                    let skin = c
                        .iter()
                        .zip(&lv)
                        .fold(Weighted::<Point>::zero(), |acc, (s, f)| {
                            acc.add(s.1[i].scale(f[l]))
                        });
                    departures
                        .iter()
                        .zip(&mu)
                        .fold(skin, |acc, (d, f)| acc.add(d[l].scale(f[i])))
                })
                .collect()
        })
        .collect();
    let v_space = KnotVector::new(space_knots(v_degree, &v_inner), v_degree)?;
    let surface = surface_of(&c[0].0, v_space, &rows, tol)?;

    // Measured: every section on its row, every guide on its column.
    let target = tol.confusion() * 10.0;
    for (k, (s, &v)) in read.iter().zip(&v_at).enumerate() {
        let worst = strays(|f| {
            Ok(surface
                .point_at(f, v, tol)?
                .distance(s.edges[0].curve.point_at(f, tol)?))
        })?;
        if worst > target {
            ogeom_bail!(
                NotDone,
                "the guided skin strays {worst:.3e} from section {k}"
            );
        }
    }
    for ((rail, &u), &j) in paced_rails.iter().zip(&u_at).zip(&names) {
        let worst = strays(|f| {
            Ok(surface
                .point_at(u, f, tol)?
                .distance(rail.point_at(f, tol)?))
        })?;
        if worst > target {
            ogeom_bail!(
                NotDone,
                "the guided skin strays {worst:.3e} from guide {j}, farther than {target:.1e}"
            );
        }
    }

    let shape = sheet(
        model,
        &read,
        &[vec![surface]],
        &[(0, n - 1)],
        false,
        false,
        tol,
    )?;
    let mut history = History::new();
    for input in sections.iter().chain(guides) {
        history.generate(input, shape.clone());
    }
    Ok(Built { shape, history })
}

/// A parameter on the unit domain, an end where it is within
/// [`ONE_PARAMETER`] of one.
fn at_ends(t: f64) -> f64 {
    if t <= ONE_PARAMETER {
        0.0
    } else if t >= 1.0 - ONE_PARAMETER {
        1.0
    } else {
        t
    }
}

/// The largest of `distance` at [`MEASURED`] steps over the unit domain.
fn strays(distance: impl Fn(f64) -> OgeomResult<f64>) -> OgeomResult<f64> {
    let mut worst: f64 = 0.0;
    for i in 0..=MEASURED {
        worst = worst.max(distance(f64::from(i) / f64::from(MEASURED))?);
    }
    Ok(worst)
}

/// A section's edges as one curve over `[0, 1]`, each edge raised to the
/// highest degree among them and its weights scaled to meet the one
/// before, joined end to end.
fn joined(section: &Section, tol: Tolerances) -> OgeomResult<BSplineCurve> {
    if section.edges.len() == 1 {
        return Ok(section.edges[0].curve.clone());
    }
    let degree = section
        .edges
        .iter()
        .map(|e| e.curve.degree())
        .max()
        .unwrap_or(1);
    let mut out: Option<Homogeneous> = None;
    for (i, e) in section.edges.iter().enumerate() {
        let mut c = e.curve.clone();
        while c.degree() < degree {
            c = c.elevated(tol)?;
        }
        #[allow(clippy::cast_precision_loss)]
        let knots = c.knots().reparameterized(i as f64, (i + 1) as f64)?;
        let mut control = c.control_points().to_vec();
        if let Some((_, before)) = &out {
            let k = before[before.len() - 1].weight / control[0].weight;
            control = control.iter().map(|p| p.scale(k)).collect();
        }
        let piece = (knots, control);
        out = Some(match out {
            None => piece,
            Some(before) => bspline::join(&before, &piece)?,
        });
    }
    let Some((knots, control)) = out else {
        ogeom_bail!(Construction, "a loft guide has no edges");
    };
    standard(&BSplineCurve::rational(knots, control)?)
}

/// The same curve run the other way over `[0, 1]`.
fn reversed(curve: &BSplineCurve) -> OgeomResult<BSplineCurve> {
    let (knots, control) = bspline::reverse(curve.knots(), curve.control_points());
    BSplineCurve::rational(knots, control)
}

/// Where a guide crosses a section, or comes nearest to it: the section's
/// parameter, the guide's, and the distance between the two points. Seeded
/// from the nearest pair of samples and polished by Gauss-Newton on the
/// difference of the two points, each parameter held to its domain.
fn crossing(
    section: &BSplineCurve,
    guide: &BSplineCurve,
    tol: Tolerances,
) -> OgeomResult<(f64, f64, f64)> {
    const ALONG_SECTION: u32 = 64;
    const ALONG_GUIDE: u32 = 256;
    let sample = |c: &BSplineCurve, count: u32| -> OgeomResult<Vec<(f64, Point)>> {
        (0..=count)
            .map(|i| {
                let f = f64::from(i) / f64::from(count);
                Ok((f, c.point_at(f, tol)?))
            })
            .collect()
    };
    let on_section = sample(section, ALONG_SECTION)?;
    let on_guide = sample(guide, ALONG_GUIDE)?;
    let mut seed = (0.0, 0.0, f64::INFINITY);
    for &(a, p) in &on_section {
        for &(b, q) in &on_guide {
            let d = p.distance(q);
            if d < seed.2 {
                seed = (a, b, d);
            }
        }
    }
    let (mut u, mut t) = (seed.0, seed.1);
    for _ in 0..64 {
        let r = section.point_at(u, tol)? - guide.point_at(t, tol)?;
        let (da, db) = (section.d1_at(u, tol)?, guide.d1_at(t, tol)?);
        let (aa, ab, bb) = (da.dot(da), da.dot(db), db.dot(db));
        let (ra, rb) = (r.dot(da), r.dot(db));
        let det = aa * bb - ab * ab;
        if det.is_nan() || det <= 1e-24 * aa * bb {
            break;
        }
        let du = (ab * rb - ra * bb) / det;
        let dt = (aa * rb - ab * ra) / det;
        let (nu, nt) = ((u + du).clamp(0.0, 1.0), (t + dt).clamp(0.0, 1.0));
        let moved = (nu - u).abs() + (nt - t).abs();
        (u, t) = (nu, nt);
        if moved <= 1e-16 {
            break;
        }
    }
    let off = section.point_at(u, tol)?.distance(guide.point_at(t, tol)?);
    Ok(if off <= seed.2 { (u, t, off) } else { seed })
}

/// The same curve with its parameter changed piecewise linearly: each
/// `pairs[i].0` (increasing over the curve's domain) carried onto
/// `pairs[i].1` (increasing), the curve outside the first and last
/// dropped. Each piece is the curve's own with its knots rescaled: the
/// same point set, run at another pace.
fn repaced(
    curve: &BSplineCurve,
    pairs: &[(f64, f64)],
    tol: Tolerances,
) -> OgeomResult<BSplineCurve> {
    let mut out: Option<Homogeneous> = None;
    for w in pairs.windows(2) {
        let ((a, c), (b, d)) = (w[0], w[1]);
        let piece = curve.segment((a, b), tol)?;
        let piece = (
            piece.knots().reparameterized(c, d)?,
            piece.control_points().to_vec(),
        );
        out = Some(match out {
            None => piece,
            Some(before) => bspline::join(&before, &piece)?,
        });
    }
    let Some((knots, control)) = out else {
        ogeom_bail!(Construction, "a curve paced over no pieces");
    };
    BSplineCurve::rational(knots, control)
}

/// The functions over `[0, 1]` that interpolate at `params`: function `i`
/// is one at `params[i]` and zero at every other. They share one clamped
/// space, of degree `params.len() - 1` up to three, its interior knots the
/// averages of the parameters; one parameter gives the constant one.
/// Returns the space's knots and each function's coefficients in it.
fn cardinals(params: &[f64], tol: Tolerances) -> OgeomResult<(KnotVector, Vec<Vec<f64>>)> {
    let n = params.len();
    if n == 1 {
        return Ok((
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1)?,
            vec![vec![1.0, 1.0]],
        ));
    }
    let degree = (n - 1).min(3);
    let mut knots = vec![0.0; degree + 1];
    #[allow(clippy::cast_precision_loss)]
    for j in 1..n - degree {
        knots.push(params[j..j + degree].iter().sum::<f64>() / degree as f64);
    }
    knots.extend(core::iter::repeat_n(1.0, degree + 1));
    let knots = KnotVector::new(knots, degree)?;
    let mut matrix = vec![vec![0.0; n]; n];
    for (row, &x) in matrix.iter_mut().zip(params) {
        let span = knots.span(x, tol)?;
        for (b, j) in knots.basis(span, x).iter().zip(span - degree..=span) {
            row[j] = *b;
        }
    }
    let inverse = inverted(matrix)?;
    let functions = (0..n)
        .map(|j| inverse.iter().map(|row| row[j]).collect())
        .collect();
    Ok((knots, functions))
}

/// Interior knots of a vector over `[0, 1]`, with their multiplicities.
fn interior(knots: &KnotVector) -> Vec<(f64, usize)> {
    knots
        .distinct()
        .into_iter()
        .filter(|(v, _)| *v > KNOT_SAME && *v < 1.0 - KNOT_SAME)
        .collect()
}

/// The least space over `[0, 1]` holding splines of every one of `spaces`:
/// the highest degree, and each interior knot at the highest multiplicity
/// it reaches once each space is raised to that degree.
fn common_space(spaces: &[&KnotVector]) -> (usize, Vec<(f64, usize)>) {
    let degree = spaces.iter().map(|k| k.degree()).max().unwrap_or(1);
    let mut union: Vec<(f64, usize)> = Vec::new();
    for space in spaces {
        let raise = degree - space.degree();
        for (value, mult) in interior(space) {
            match union
                .iter_mut()
                .find(|(v, _)| (*v - value).abs() <= KNOT_SAME)
            {
                Some(entry) => entry.1 = entry.1.max(mult + raise),
                None => union.push((value, mult + raise)),
            }
        }
    }
    union.sort_by(|a, b| a.0.total_cmp(&b.0));
    (degree, union)
}

/// The clamped knots over `[0, 1]` of a space of `degree` with `inner`
/// knots.
fn space_knots(degree: usize, inner: &[(f64, usize)]) -> Vec<f64> {
    let mut knots = vec![0.0; degree + 1];
    for &(value, mult) in inner {
        knots.extend(core::iter::repeat_n(value, mult));
    }
    knots.extend(core::iter::repeat_n(1.0, degree + 1));
    knots
}

/// A spline over `[0, 1]` restated in a larger space, unchanged: raised to
/// its degree and refined to its knots.
fn into_space<P: Blend>(
    spline: &Spline<P>,
    degree: usize,
    inner: &[(f64, usize)],
    tol: Tolerances,
) -> OgeomResult<Spline<P>> {
    let (mut knots, mut control) = spline.clone();
    while knots.degree() < degree {
        (knots, control) = bspline::elevate_degree(&knots, &control, tol)?;
    }
    for &(value, mult) in inner {
        let (at, have) = interior(&knots)
            .into_iter()
            .find(|(v, _)| (*v - value).abs() <= KNOT_SAME)
            .unwrap_or((value, 0));
        if have < mult {
            (knots, control) = bspline::insert_knot(&knots, &control, at, mult - have, tol)?;
        }
    }
    let target = space_knots(degree, inner);
    let same = target.len() == knots.knots().len()
        && target
            .iter()
            .zip(knots.knots())
            .all(|(a, b)| (a - b).abs() <= KNOT_SAME * 10.0);
    if !same {
        ogeom_bail!(
            Construction,
            "the guided loft's curves could not be brought to one knot vector"
        );
    }
    Ok((KnotVector::new(target, degree)?, control))
}

/// The product of a scalar spline and a homogeneous one over `[0, 1]`: of
/// the summed degree, each knot as smooth as the rougher factor is there,
/// found by interpolating the product at the space's Greville abscissae,
/// where it is exact.
fn product(scale: &Spline<f64>, curve: &Homogeneous, tol: Tolerances) -> OgeomResult<Homogeneous> {
    let (ds, dc) = (scale.0.degree(), curve.0.degree());
    let degree = ds + dc;
    let mut inner: Vec<(f64, usize)> = Vec::new();
    for (space, own) in [(&scale.0, ds), (&curve.0, dc)] {
        for (value, mult) in interior(space) {
            let smooth = own - mult;
            match inner
                .iter_mut()
                .find(|(v, _)| (*v - value).abs() <= KNOT_SAME)
            {
                Some(entry) => entry.1 = entry.1.max(degree - smooth),
                None => inner.push((value, degree - smooth)),
            }
        }
    }
    inner.sort_by(|a, b| a.0.total_cmp(&b.0));
    let knots = KnotVector::new(space_knots(degree, &inner), degree)?;
    let flat = knots.knots();
    let count = flat.len() - degree - 1;
    let mut matrix = vec![vec![0.0; count]; count];
    let mut values = Vec::with_capacity(count);
    for (i, row) in matrix.iter_mut().enumerate() {
        #[allow(clippy::cast_precision_loss)]
        let x = flat[i + 1..=i + degree].iter().sum::<f64>() / degree as f64;
        let span = knots.span(x, tol)?;
        for (b, j) in knots.basis(span, x).iter().zip(span - degree..=span) {
            row[j] = *b;
        }
        let s = bspline::evaluate(&scale.0, &scale.1, x, tol)?;
        values.push(bspline::evaluate(&curve.0, &curve.1, x, tol)?.scale(s));
    }
    let inverse = inverted(matrix)?;
    let control = inverse
        .iter()
        .map(|line| {
            line.iter()
                .zip(&values)
                .fold(Weighted::<Point>::zero(), |acc, (a, p)| {
                    acc.add(p.scale(*a))
                })
        })
        .collect();
    Ok((knots, control))
}
