//! Pcurves fitted by projection, for faces whose trims have no closed form.
//!
//! An exchange file's edge must end up with a curve in each bounding face's
//! parameters, or the face cannot be split or triangulated. Where the
//! curve/surface pair has a closed form the exact projection is used; where
//! it does not (a spline surface, mostly), the pcurve is *fitted at the
//! curve's own parameters*: sample the edge, project each sample into the
//! chart, fit the trace with the parameters held fixed, so the same-parameter
//! law holds by construction. An exact curve never carries a fitted pcurve
//! silently: the fit's error is returned, and the callers widen tolerances
//! and warn with it. That error
//! is reported as a *length*, in the model's own units: the fitted pcurve is
//! walked through the surface and compared against the trace it was fitted
//! to. A chart's units are whatever the file chose, and no single scale
//! converts them: a patch can span four microns across its `u` and ten
//! millimetres along its `v`.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::Surface as _;
use ogeom_geom::{Curve, PlanarCurve, SurfaceGeometry};
use ogeom_math::Point;

/// The fitted pcurve, its fit error as a length, whether the target was met,
/// the worst distance any sample sat from the surface, and the slop warning
/// to record when that distance is large enough to say out loud.
pub type FittedPcurve = OgeomResult<(PlanarCurve, f64, bool, f64, Option<String>)>;

pub(crate) fn chart_of(surface: &SurfaceGeometry, p: Point) -> Option<ogeom_math::Point2> {
    let tau = core::f64::consts::TAU;
    match surface {
        SurfaceGeometry::Plane(s) => {
            let l = s.plane().frame().to_local(p);
            Some(ogeom_math::Point2::new(l.x, l.y))
        }
        SurfaceGeometry::Cylinder(s) => {
            let l = s.cylinder().frame().to_local(p);
            Some(ogeom_math::Point2::new(l.y.atan2(l.x).rem_euclid(tau), l.z))
        }
        SurfaceGeometry::Cone(s) => {
            let l = s.cone().frame().to_local(p);
            Some(ogeom_math::Point2::new(l.y.atan2(l.x).rem_euclid(tau), l.z))
        }
        SurfaceGeometry::Sphere(s) => {
            let sphere = s.sphere();
            let l = sphere.frame().to_local(p);
            let lat = (l.z / sphere.radius()).clamp(-1.0, 1.0).asin();
            Some(ogeom_math::Point2::new(l.y.atan2(l.x).rem_euclid(tau), lat))
        }
        SurfaceGeometry::Torus(s) => {
            let torus = s.torus();
            let l = torus.frame().to_local(p);
            let u = l.y.atan2(l.x).rem_euclid(tau);
            let radial = l.x.hypot(l.y) - torus.major_radius();
            let v = l.z.atan2(radial).rem_euclid(tau);
            Some(ogeom_math::Point2::new(u, v))
        }
        _ => None,
    }
}

/// Where a point of the curve lands on the surface: its chart position and
/// how far off the surface it sat.
///
/// Analytic surfaces invert in closed form: grid seeding over a plane's
/// or cylinder's enormous stated extents lands microns off, and a fitted
/// pcurve inherits every micron. On a patch, where the previous sample
/// landed is a far better starting guess than any grid: consecutive
/// samples of a curve are neighbouring points of the surface. Trusted only
/// when it lands convincingly *on* the surface (the same bar the denser
/// reseed is judged against), so a guess that wandered into the wrong
/// basin, or a first sample with no predecessor, still pays for the grid.
fn land(
    surface: &SurfaceGeometry,
    p: Point,
    seed: Option<ogeom_math::Point2>,
    tol: Tolerances,
) -> OgeomResult<(ogeom_math::Point2, f64)> {
    if let Some(uv) = chart_of(surface, p) {
        let lifted = surface.point_at(uv.x, uv.y, tol)?;
        return Ok((uv, p.distance(lifted)));
    }
    let near = seed.and_then(|luv| {
        crate::measure::project_on_surface_from(surface, p, (luv.x, luv.y), tol).ok()
    });
    if let Some(close) = near.filter(|f| f.distance <= tol.confusion() * 1e5) {
        return Ok((
            ogeom_math::Point2::new(close.parameters.0, close.parameters.1),
            close.distance,
        ));
    }
    let mut projection = crate::measure::project_on_surface(surface, p, 24, tol)?;
    if projection.distance > tol.confusion() * 1e5 {
        // A miss this large on a spline surface is more often a projection
        // stuck in the wrong basin than real slop; seed denser before
        // believing it.
        let denser = crate::measure::project_on_surface(surface, p, 96, tol)?;
        if denser.distance < projection.distance {
            projection = denser;
        }
    }
    Ok((
        ogeom_math::Point2::new(projection.parameters.0, projection.parameters.1),
        projection.distance,
    ))
}

/// A sample of the curve landed on the surface: the parameter, the point,
/// where it landed in the chart, and how far off the surface it sat.
type Landed = (f64, Point, ogeom_math::Point2, f64);

/// Fit a pcurve by projection at the reader's own line: slop under a
/// millimetre is a file's error, honestly carried; anything past it is a
/// wrong pairing and refuses. The exchange readers call this; a healer
/// acting on instruction calls [`fit_projected_pcurve_capped`] with the
/// cap its caller chose.
///
/// # Errors
///
/// As [`fit_projected_pcurve_capped`], at the millimetre cap.
pub fn fit_projected_pcurve(
    curve: &Curve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> FittedPcurve {
    fit_projected_pcurve_capped(curve, range, surface, tol.confusion() * 1e7, tol)
}

/// As the reader's projected-pcurve fit, with the acceptance cap in the
/// caller's hands.
///
/// The reader draws its line at a millimetre (below it is a file's own
/// slop, above it a wrong pairing), but a *healer* acts on instruction, and
/// the instruction carries the cap. Returns the fitted pcurve, the fit's
/// reached error as a length (the widest the lifted fit misses where the
/// curve lands, between the samples as well as at them), whether it met its
/// target, the worst measured edge-to-surface offset, and the slop note
/// when that offset is worth saying out loud.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a
/// sample sits farther than `cap` from the surface, or the projection
/// cannot converge at all.
pub fn fit_projected_pcurve_capped(
    curve: &Curve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    cap: f64,
    tol: Tolerances,
) -> FittedPcurve {
    match fit_projected_pcurve_within(curve, range, surface, cap, tol)? {
        CappedFit::Fitted(fitted) => Ok(fitted),
        CappedFit::TooFar(worst_off) => ogeom_bail!(
            Construction,
            "the edge sits {worst_off:.2e} from the surface it should bound"
        ),
    }
}

/// What [`fit_projected_pcurve_within`] made of an edge.
#[derive(Debug, Clone)]
pub enum CappedFit {
    /// The fit, as [`fit_projected_pcurve_capped`] returns it.
    Fitted((PlanarCurve, f64, bool, f64, Option<String>)),
    /// Refused: some sample sits this far from the surface, past the cap.
    TooFar(f64),
}

/// [`fit_projected_pcurve_capped`], with a refusal at the cap returned as
/// the offset measured rather than as an error.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// the projection cannot converge at all.
pub fn fit_projected_pcurve_within(
    curve: &Curve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    cap: f64,
    tol: Tolerances,
) -> OgeomResult<CappedFit> {
    const SAMPLES: usize = 96;
    let mut parameters = Vec::with_capacity(SAMPLES + 1);
    let mut trace = Vec::with_capacity(SAMPLES + 1);
    let mut points = Vec::with_capacity(SAMPLES + 1);
    let mut offs = Vec::with_capacity(SAMPLES + 1);
    let mut previous: Option<(Point, ogeom_math::Point2)> = None;
    for i in 0..=SAMPLES {
        #[allow(clippy::cast_precision_loss)]
        let t = range.0 + (range.1 - range.0) * i as f64 / SAMPLES as f64;
        let p = curve.point_at(t, tol)?;
        let (uv, off) = land(surface, p, previous.map(|(_, luv)| luv), tol)?;
        previous = Some((p, uv));
        parameters.push(t);
        trace.push(uv);
        points.push(p);
        offs.push(off);
    }
    retry_stalled(surface, &points, &mut trace, &mut offs, tol);
    // The cap separates a file's own slop from an edge paired with the
    // wrong surface. Slop is routinely a micron or two and reaches a third
    // of a millimetre in exchange files, while a wrong pairing misses by
    // the distance between two different surfaces of the body, whole
    // millimetres. One millimetre stands between the slop a file carries and
    // the smallest wrong pairing plausible. Slop inside the
    // cap is accepted and *recorded*: the edge's tolerance is widened to
    // cover it, so the model says what it knows instead of refusing to
    // triangulate. Judged after the retry, so a projection that stalled
    // cannot refuse an edge that is actually on its surface.
    let worst_off = offs.iter().copied().fold(0.0_f64, f64::max);
    if worst_off > cap {
        return Ok(CappedFit::TooFar(worst_off));
    }
    let mut space_run = 0.0;
    let mut parameter_run = 0.0;
    for i in 1..trace.len() {
        space_run += points[i].distance(points[i - 1]);
        parameter_run += trace[i].distance(trace[i - 1]);
    }
    // A trace on a periodic chart may cross the seam mid-edge; unwrap it
    // pointwise so the fit sees a continuous curve. Closure, not
    // periodicity, is the right test: a skinned loft's wall is a clamped
    // B-spline that closes on itself without being periodic, and its
    // projections near the joining column land in either copy: both
    // answers are right pointwise, and only continuity chooses. It is
    // decided here for both exchange readers at once.
    let ((ua, ub), (va, vb)) = surface.domain();
    let spans = (
        if surface.is_periodic_u() || surface.is_closed_u(tol) {
            ub - ua
        } else {
            0.0
        },
        if surface.is_periodic_v() || surface.is_closed_v(tol) {
            vb - va
        } else {
            0.0
        },
    );
    for i in 1..trace.len() {
        if spans.0 > 0.0 {
            while trace[i].x - trace[i - 1].x > spans.0 * 0.5 {
                trace[i].x -= spans.0;
            }
            while trace[i].x - trace[i - 1].x < -spans.0 * 0.5 {
                trace[i].x += spans.0;
            }
        }
        if spans.1 > 0.0 {
            while trace[i].y - trace[i - 1].y > spans.1 * 0.5 {
                trace[i].y -= spans.1;
            }
            while trace[i].y - trace[i - 1].y < -spans.1 * 0.5 {
                trace[i].y += spans.1;
            }
        }
    }
    slide_into_chart(&mut trace, spans, ((ua, ub), (va, vb)));
    // Where the chart collapses (a sphere's pole, a cone's apex), the
    // u of a sample is atan2 of noise: the point determines no angle.
    // The *arc* does: a smooth curve through the pole approaches it at
    // a definite chart angle, which is the limit of its well-conditioned
    // neighbours. Samples whose u-direction has collapsed relative to
    // their v-direction are repaired by interpolating u between the
    // nearest sound samples, extrapolating at the ends.
    //
    // Weak is measured in millimetres, not against `dv`. A ratio calls a
    // direction weak whenever the *other* one is strong: on a patch whose
    // `v` is parameterised a thousand times more densely than its `u`,
    // every sample of every edge is called weak, its `u` held at one value,
    // and an edge across the patch fits as a single point. What makes a
    // direction degenerate is that crossing the whole of it moves the point
    // less than a micron; that question has an answer in length, and only
    // in length.
    let (u_span, _) = {
        let ((ua, ub), (va, vb)) = surface.domain();
        (ub - ua, vb - va)
    };
    let weak: Vec<bool> = trace
        .iter()
        .map(|uv| {
            surface
                .d1_at(uv.x, uv.y, tol)
                .is_ok_and(|(du, _)| du.magnitude() * u_span < tol.confusion() * 1e4)
        })
        .collect();
    if weak.iter().all(|w| *w) && !weak.is_empty() {
        // Not a row that collapses but a whole patch that does: a sliver
        // four microns wide and a tenth of a millimetre long, where `u` is
        // noise everywhere and the projector answers 0 at one sample and 1
        // at the next. A fit through that swings across the chart and its
        // controls are dragged back by hundreds of units. Any `u` describes
        // the same points to within the sliver's own width, so they all
        // take one: the middle of what the projections claimed, which is
        // the least arbitrary of the arbitrary answers.
        let mut claimed: Vec<f64> = trace.iter().map(|uv| uv.x).collect();
        claimed.sort_by(f64::total_cmp);
        let held = claimed[claimed.len() / 2];
        for uv in &mut trace {
            uv.x = held;
        }
    } else if weak.iter().any(|w| *w) && weak.iter().filter(|w| !**w).count() >= 2 {
        let strong: Vec<usize> = (0..trace.len()).filter(|&i| !weak[i]).collect();
        let u_span = if surface.is_periodic_u() {
            ua.max(ub) - ua.min(ub)
        } else {
            f64::INFINITY
        };
        for i in 0..trace.len() {
            if !weak[i] {
                continue;
            }
            let after = strong.iter().position(|&s| s > i);
            let (a, b) = match after {
                Some(0) => (strong[0], strong[1]),
                Some(k) => (strong[k - 1], strong[k]),
                None => (strong[strong.len() - 2], strong[strong.len() - 1]),
            };
            // A curve *through* the pole genuinely jumps its angle
            // there; only a run whose sound neighbours agree is noise
            // to smooth over.
            if a < i && i < b && (trace[b].x - trace[a].x).abs() > u_span * 0.25 {
                continue;
            }
            let (ta, tb) = (parameters[a], parameters[b]);
            let f = if (tb - ta).abs() <= f64::MIN_POSITIVE {
                0.0
            } else {
                (parameters[i] - ta) / (tb - ta)
            };
            trace[i].x = trace[a].x + (trace[b].x - trace[a].x) * f;
        }
    }

    // The tolerance carried into the chart through the trace's own
    // metric: the honest cheap version, refined by the fit's report.
    let scale = if space_run > tol.confusion() {
        parameter_run / space_run
    } else {
        1.0
    };
    let target = (tol.confusion() * 1e2 * scale).max(f64::MIN_POSITIVE);
    let fit_and_clamp = |parameters: &[f64],
                         trace: &[ogeom_math::Point2]|
     -> OgeomResult<ogeom_geom::fit::Fitted<ogeom_geom::BSpline2d>> {
        let fitted = ogeom_geom::fit::fit_points_2d_at(parameters, trace, 3, target, tol)?;
        // A least-squares fit wiggles past its samples at the ends, and a
        // surface with a *tight* stated window (an imported patch, not a
        // reader-built analytic with its enormous extents) refuses
        // evaluation a hair outside it. The control points clamp into the
        // window on the non-periodic axes: the curve lives in its controls'
        // hull, so the clamp is a guarantee. Whatever the clamp cost is not
        // hidden either: it is the clamped curve that is measured below.
        // Periodic axes stay free: an unwrapped trace crosses the seam on
        // purpose.
        let ((wa, wb), (va2, vb2)) = surface.domain();
        let clamp_u = !(surface.is_periodic_u() || surface.is_closed_u(tol));
        let clamp_v = !(surface.is_periodic_v() || surface.is_closed_v(tol));
        if !(clamp_u || clamp_v) {
            return Ok(fitted);
        }
        let mut moved = 0.0_f64;
        let knots = fitted.curve.knots().clone();
        let control: Vec<ogeom_math::Point2> = fitted
            .curve
            .control_points()
            .iter()
            .map(|w| {
                let p = w.point();
                let q = ogeom_math::Point2::new(
                    if clamp_u { p.x.clamp(wa, wb) } else { p.x },
                    if clamp_v { p.y.clamp(va2, vb2) } else { p.y },
                );
                moved = moved.max(p.distance(q));
                q
            })
            .collect();
        if moved > 0.0 {
            Ok(ogeom_geom::fit::Fitted {
                curve: ogeom_geom::BSpline2d::new(knots, control, tol)?,
                ..fitted
            })
        } else {
            Ok(fitted)
        }
    };
    // What the caller is told, as a length. The fitter reports its error
    // in *chart* units, and a chart's units are whatever the file chose: a
    // patch can span four microns across its `u` and ten millimetres along
    // its `v`, so no single scale converts the one number into the other,
    // and a control point dragged back into such a chart by the clamp reads
    // as hundreds of millimetres of mesh error on a face a tenth of a
    // millimetre across. So the fitted curve is walked instead, through the
    // surface, against the trace it was fitted to: that difference is the
    // fit's own and is measured where the mesh will be.
    //
    // Measured between the samples as well as at them. The samples are
    // spaced evenly in the curve's parameter, and a curve is free to run
    // many times faster at one end than the other: one that turns through
    // most of its bend inside the first interval leaves a cubic held only
    // at the interval's ends hooking well past it, off the face and across
    // its neighbouring ring, while every sample sits within a hundredth. A
    // midpoint the fit leaves is projected and joins the samples, and the
    // fit is asked again, a few rounds at most.
    let closed_form = points
        .first()
        .is_some_and(|p| chart_of(surface, *p).is_some());
    // A landing is believed only where it belongs: on the surface as
    // convincingly as its neighbours, and inside the chart interval
    // they span, widened by the interval itself. A projection that
    // settled in another basin, or on the far side of a seam, would
    // otherwise be fitted as if the curve went there. Beside its
    // neighbour on a periodic chart, as the trace was unwrapped; where
    // the chart collapses its angle is noise, and the neighbours' is
    // taken.
    let landing =
        |trace: &[ogeom_math::Point2], offs: &[f64], index: usize, tm: f64| -> Option<Landed> {
            let before = trace[index - 1];
            let after = trace[index];
            let p = curve.point_at(tm, tol).ok()?;
            let (mut uv, off) = land(surface, p, Some(before), tol).ok()?;
            if spans.0 > 0.0 {
                while uv.x - before.x > spans.0 * 0.5 {
                    uv.x -= spans.0;
                }
                while uv.x - before.x < -spans.0 * 0.5 {
                    uv.x += spans.0;
                }
            }
            if spans.1 > 0.0 {
                while uv.y - before.y > spans.1 * 0.5 {
                    uv.y -= spans.1;
                }
                while uv.y - before.y < -spans.1 * 0.5 {
                    uv.y += spans.1;
                }
            }
            if surface
                .d1_at(uv.x, uv.y, tol)
                .is_ok_and(|(du, _)| du.magnitude() * u_span < tol.confusion() * 1e4)
            {
                uv.x = 0.5 * (before.x + after.x);
            }
            let reach = before.distance(after).max(f64::EPSILON);
            let mid =
                ogeom_math::Point2::new(0.5 * (before.x + after.x), 0.5 * (before.y + after.y));
            let sound = offs[index - 1].max(offs[index]).max(tol.confusion() * 1e5);
            (off <= 2.0 * sound && mid.distance(uv) <= reach).then_some((tm, p, uv, off))
        };
    let deviation = |fitted: &ogeom_geom::fit::Fitted<ogeom_geom::BSpline2d>,
                     parameters: &[f64],
                     trace: &[ogeom_math::Point2],
                     offs: &[f64]|
     -> (f64, f64, Vec<Landed>) {
        let mut error = 0.0_f64;
        let mut between_all = 0.0_f64;
        let mut more = Vec::new();
        let mut landed_middles = Vec::new();
        let mut left = false;
        let landing = |index: usize, tm: f64| landing(trace, offs, index, tm);
        for (index, t) in parameters.iter().enumerate() {
            let Ok(at) = ogeom_geom::Curve2d::point_at(&fitted.curve, *t, tol) else {
                continue;
            };
            let (Ok(fitted_at), Ok(traced_at)) = (
                surface.point_at(at.x, at.y, tol),
                surface.point_at(trace[index].x, trace[index].y, tol),
            ) else {
                continue;
            };
            error = error.max(fitted_at.distance(traced_at));
            if index == 0 {
                continue;
            }
            // Probed at the quarters as well as the middle where the chart
            // inverts in closed form, which costs nothing: a hook sits
            // where the curve turns, wherever in the interval that is. On a
            // patch every landing is a projection, and the middle alone is
            // asked: a hook is the cubic's own excursion, broad across the
            // interval, and once the fit is asked again at twice the
            // samples the quarters of this round are the middles of the
            // next.
            let (ta, tb) = (parameters[index - 1], *t);
            let tm = 0.5 * (ta + tb);
            let quarters = [0.5 * (ta + tm), tm, 0.5 * (tm + tb)];
            let probes: &[f64] = if closed_form {
                &quarters
            } else {
                &quarters[1..2]
            };
            for &probe in probes {
                let Ok(at) = ogeom_geom::Curve2d::point_at(&fitted.curve, probe, tol) else {
                    continue;
                };
                let Ok(fitted_at) = surface.point_at(at.x, at.y, tol) else {
                    continue;
                };
                if !closed_form {
                    // Against the curve's own point first, which costs no
                    // projection: the curve sits about as far off the
                    // surface here as at the samples either side, so a
                    // fitted point within the bar and that slop of it is
                    // within the bar of where the curve lands. A hook worth
                    // the name is hundreds of times the slop; only a probe
                    // that misses by more than the slop is projected.
                    let Ok(p) = curve.point_at(probe, tol) else {
                        continue;
                    };
                    let slop = offs[index - 1].max(offs[index]);
                    let rough = fitted_at.distance(p);
                    if rough <= tol.confusion() * 1e4 + slop {
                        between_all = between_all.max((rough - slop).max(0.0));
                        continue;
                    }
                }
                let Some(landed) = landing(index, probe) else {
                    continue;
                };
                let Ok(traced_at) = surface.point_at(landed.2.x, landed.2.y, tol) else {
                    continue;
                };
                let between = fitted_at.distance(traced_at);
                between_all = between_all.max(between);
                if between > tol.confusion() * 1e4 {
                    left = true;
                }
                if probe == tm {
                    landed_middles.push(landed);
                }
            }
        }
        // The fit is asked again at twice the samples everywhere, the
        // middle of every interval joining them, projected here where the
        // rough test spared it.
        if left {
            more = landed_middles;
            if !closed_form {
                more = (1..parameters.len())
                    .filter_map(|index| {
                        landing(index, 0.5 * (parameters[index - 1] + parameters[index]))
                    })
                    .collect();
            }
        }
        (error, between_all, more)
    };
    // A micron between samples is the bar, a tenth of the finest chord a
    // mesh is asked for, and a hook worth the name is hundreds of times
    // that. Every fit pays one pass of probes; only the few that leave the
    // curve pay a refit, at twice the samples everywhere (a handful of new
    // samples in one interval draw the fitter's knots to themselves and the
    // curve wobbles on either side, while an even doubling keeps it
    // steady), and a curve the fit cannot follow, a corner inside the edge,
    // stops at a few hundred samples rather than doubling for ever.
    const DENSIFY: usize = 6;
    const MOST: usize = 512;
    let mut fitted = fit_and_clamp(&parameters, &trace)?;
    #[allow(
        clippy::type_complexity,
        reason = "a fit, its two misses and its sample count"
    )]
    let mut best: Option<(
        ogeom_geom::fit::Fitted<ogeom_geom::BSpline2d>,
        f64,
        f64,
        usize,
    )> = None;
    let mut round = 0;
    loop {
        let (at_samples, between, more) = deviation(&fitted, &parameters, &trace, &offs);
        // The best round stands, whichever it was: a refit is not obliged
        // to improve, and the fit handed on is the one measured closest.
        if best
            .as_ref()
            .is_none_or(|(_, a, b, _)| at_samples.max(between) < a.max(*b))
        {
            best = Some((fitted.clone(), at_samples, between, parameters.len()));
        }
        if more.is_empty() || round == DENSIFY || parameters.len() >= MOST {
            break;
        }
        for (tm, p, uv, off) in more {
            let k = parameters.partition_point(|&t| t < tm);
            parameters.insert(k, tm);
            trace.insert(k, uv);
            points.insert(k, p);
            offs.insert(k, off);
        }
        fitted = fit_and_clamp(&parameters, &trace)?;
        round += 1;
    }
    let (fitted, at_samples, between, own) =
        best.unwrap_or((fitted, f64::INFINITY, f64::INFINITY, 0));
    // The probes above decide where to refit; the fit handed on states the
    // widest miss its quarter points and their peaks show.
    let between = between.max(climbed_miss(
        &fitted.curve,
        (&parameters, own),
        &trace,
        &offs,
        closed_form,
        |t| curve.point_at(t, tol).ok(),
        surface,
        tol,
    ));
    let worst_off = offs.iter().copied().fold(worst_off, f64::max);
    // Each miss against its own bar: the samples hold the fit to a hair,
    // and the points between them to the micron that sent it back.
    let error = at_samples.max(between);
    let met = at_samples <= tol.confusion() * 1e2 && between <= tol.confusion() * 1e4;
    let slop = (worst_off > tol.confusion() * 1e3).then(|| {
        format!(
            "an edge sits up to {worst_off:.2e} from the surface it \
             bounds; the file's own slop, carried into the chart"
        )
    });
    Ok(CappedFit::Fitted((
        fitted.curve.into(),
        error,
        met,
        worst_off,
        slop,
    )))
}

/// The widest a fitted pcurve, lifted through `surface`, misses where the
/// curve lands between the samples of `parameters` (`trace` where they
/// landed, `offs` how far off the surface), and at those the fit was not
/// fitted through: the samples past the first `own`, which a refit added.
///
/// Every interval's quarter points are read, and the widest few readings,
/// those within a fifth of the widest, at most `PEAKS`, are climbed to the
/// tops of their peaks by golden sections between the quarters either
/// side. Where the chart inverts in closed form (`closed_form`) a reading
/// is the miss itself. On a patch every landing is a projection, and a
/// reading is the lifted fit's distance from the curve's own point
/// instead, which differs from the miss by no more than the curve's offset
/// from the surface there; a top is stated with the wider offset either
/// side of it added.
#[allow(
    clippy::too_many_arguments,
    reason = "the fit's context, each piece read"
)]
fn climbed_miss(
    pcurve: &ogeom_geom::BSpline2d,
    (parameters, own): (&[f64], usize),
    trace: &[ogeom_math::Point2],
    offs: &[f64],
    closed_form: bool,
    curve_at: impl Fn(f64) -> Option<Point>,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> f64 {
    const PEAKS: usize = 4;
    const NEAR_PEAK: f64 = 0.8;
    let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
    let lifted = |t: f64| -> Option<Point> {
        let at = ogeom_geom::Curve2d::point_at(pcurve, t, tol).ok()?;
        surface.point_at(at.x, at.y, tol).ok()
    };
    let reading = |t: f64| -> f64 {
        let on = if closed_form {
            curve_at(t)
                .and_then(|p| chart_of(surface, p))
                .and_then(|uv| surface.point_at(uv.x, uv.y, tol).ok())
        } else {
            curve_at(t)
        };
        match (lifted(t), on) {
            (Some(fit), Some(on)) => fit.distance(on),
            _ => 0.0,
        }
    };
    let mut widest = 0.0_f64;
    // The samples a refit added after this fit's round stand between its
    // own.
    if parameters.len() != own {
        for (&t, uv) in parameters.iter().zip(trace) {
            if let (Some(fit), Ok(traced)) = (lifted(t), surface.point_at(uv.x, uv.y, tol)) {
                widest = widest.max(fit.distance(traced));
            }
        }
    }
    // Each interval's widest quarter reading, where it stands, and the
    // quarter's length.
    let mut readings: Vec<(f64, usize, f64, f64)> = (1..parameters.len())
        .map(|index| {
            let (a, b) = (parameters[index - 1], parameters[index]);
            let quarter = (b - a) / 4.0;
            (1..=3)
                .map(|k| {
                    let t = quarter.mul_add(f64::from(k), a);
                    (reading(t), index, t, quarter)
                })
                .fold((f64::NEG_INFINITY, index, a, quarter), |held, r| {
                    if r.0 > held.0 { r } else { held }
                })
        })
        .collect();
    readings.sort_by(|a, b| b.0.total_cmp(&a.0));
    let tallest = readings.first().map_or(0.0, |r| r.0);
    for &(height, index, at, quarter) in readings
        .iter()
        .take_while(|r| r.0 > 0.0 && r.0 >= tallest * NEAR_PEAK)
        .take(PEAKS)
    {
        let (lo, hi) = (parameters[index - 1], parameters[index]);
        let (mut a, mut b) = ((at - quarter).max(lo), (at + quarter).min(hi));
        let (mut c, mut d) = (b - (b - a) * ratio, a + (b - a) * ratio);
        let (mut rc, mut rd) = (reading(c), reading(d));
        let mut top = height.max(rc).max(rd);
        for _ in 0..20 {
            if rc > rd {
                (b, d, rd) = (d, c, rc);
                c = b - (b - a) * ratio;
                rc = reading(c);
            } else {
                (a, c, rc) = (c, d, rd);
                d = a + (b - a) * ratio;
                rd = reading(d);
            }
            top = top.max(rc).max(rd);
        }
        let slop = if closed_form {
            0.0
        } else {
            offs[index - 1].max(offs[index])
        };
        widest = widest.max(top + slop);
    }
    widest
}

/// How far `pcurve`, lifted through `surface`, stands from `curve`, each
/// over its range and read at the same fraction along it.
///
/// This is the width an edge carrying both must state: the curve sitting
/// off the surface (an exact pcurve is the curve's projection, and a file's
/// line can lie microns off its plane) and a fit straying between the
/// samples it was fitted through both show here. The lifted point is
/// compared with the curve's point at the same fraction; unless `paced`
/// (an edge claiming same parameter is held to its claim), where that
/// misses by more than the confusion tolerance it is compared with the
/// nearest point of the curve within one sample's step, so a pcurve paced
/// a little differently is not charged for its pace. Samples the surface
/// cannot evaluate (past a pole) are passed over.
///
/// A fit strays most where its chart runs fast, as a curve passing close
/// by a sphere's pole sweeps across the chart in a few thousandths of its
/// length, and a stray that narrow falls between even samples. So the
/// widest few local peaks of the samples, and the two ends where a fit's
/// trace is seeded, are each searched for their summit between the samples
/// either side. Unless `paced`, the peaks are those of the lifted point's
/// distance square off the curve's tangent, which is the nearest-point gap
/// to first order: a pcurve paced a little differently strays along the
/// curve most where it strays across it least, and the peaks of the paced
/// distance stand away from those of the gap stated.
///
/// # Errors
///
/// Propagates a failure to evaluate the curve or the pcurve over its range.
pub fn lifted_gap(
    (curve, range): (&Curve, (f64, f64)),
    (pcurve, pcurve_range): (&PlanarCurve, (f64, f64)),
    surface: &SurfaceGeometry,
    paced: bool,
    tol: Tolerances,
) -> OgeomResult<f64> {
    use ogeom_geom::Curve2d as _;
    // Eight times the samples `check` takes, so its samples are among them.
    const SAMPLES: u32 = 256;
    // The peaks searched besides the ends: the widest stray and its
    // runners-up standing within `NEAR_PEAK` of it, at most `PEAKS`, since
    // the widest sample need not stand under the widest summit. A fit's
    // wiggle can be as narrow as these samples are apart, and a sample
    // then reads its summit a fifth low.
    const PEAKS: usize = 8;
    const NEAR_PEAK: f64 = 0.7;
    let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
    let n = f64::from(SAMPLES);
    let step = (range.1 - range.0) / n;
    let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
    // The lifted point at fraction `f` and its distance from the curve's
    // point there, or `None` where the surface cannot be evaluated.
    let paced_at = |f: f64| -> OgeomResult<Option<(Point, f64)>> {
        let uv = pcurve.point_at(
            (pcurve_range.1 - pcurve_range.0).mul_add(f, pcurve_range.0),
            tol,
        )?;
        let Ok(lifted) = surface.point_at(uv.x, uv.y, tol) else {
            return Ok(None);
        };
        let t = (range.1 - range.0).mul_add(f, range.0);
        Ok(Some((lifted, curve.point_at(t, tol)?.distance(lifted))))
    };
    // The part of the paced distance at fraction `f` square off the
    // curve's tangent there: the height the peaks are picked by.
    let across_at = |f: f64, (lifted, gap): (Point, f64)| -> OgeomResult<f64> {
        if paced {
            return Ok(gap);
        }
        let t = (range.1 - range.0).mul_add(f, range.0);
        let tangent = curve.d1_at(t, tol)?;
        let length = tangent.magnitude();
        if length <= f64::MIN_POSITIVE {
            return Ok(gap);
        }
        let along = (lifted - curve.point_at(t, tol)?).dot(tangent) / length;
        Ok(along.mul_add(-along, gap * gap).max(0.0).sqrt())
    };
    // The gap the edge must state at fraction `f`, given the lifted point
    // and its paced distance there: never more than the paced distance.
    let gap_at = |f: f64, (lifted, gap): (Point, f64)| -> OgeomResult<f64> {
        if paced || gap <= tol.confusion() {
            return Ok(gap);
        }
        let at = |s: f64| -> OgeomResult<f64> { Ok(curve.point_at(s, tol)?.distance(lifted)) };
        let t = (range.1 - range.0).mul_add(f, range.0);
        // A lifted point standing square off the curve's tangent is already
        // as near the curve's point here as any, to well under a part in a
        // thousand of the gap.
        let tangent = curve.d1_at(t, tol)?;
        let along = (lifted - curve.point_at(t, tol)?).dot(tangent).abs();
        if along <= gap * tangent.magnitude() * 1e-2 {
            return Ok(gap);
        }
        let (mut a, mut b) = ((t - step.abs()).max(lo), (t + step.abs()).min(hi));
        for _ in 0..24 {
            let (c, d) = (b - (b - a) * ratio, a + (b - a) * ratio);
            if at(c)? < at(d)? {
                b = d;
            } else {
                a = c;
            }
        }
        Ok(gap.min(at(f64::midpoint(a, b))?))
    };
    let samples = (0..=SAMPLES)
        .map(|i| paced_at(f64::from(i) / n))
        .collect::<OgeomResult<Vec<_>>>()?;
    let read = |i: usize| samples[i].map_or(0.0, |(_, gap)| gap);
    let across = (0..=SAMPLES)
        .zip(&samples)
        .map(|(i, sample)| sample.map_or(Ok(0.0), |sample| across_at(f64::from(i) / n, sample)))
        .collect::<OgeomResult<Vec<_>>>()?;

    // The widest sample, each measured fully only while its paced distance,
    // which bounds it, could still beat the widest found.
    let mut order: Vec<usize> = (0..samples.len()).collect();
    order.sort_by(|&a, &b| read(b).total_cmp(&read(a)));
    let mut widest = 0.0_f64;
    for &i in &order {
        let Some(sample) = samples[i] else { continue };
        if sample.1 <= widest {
            break;
        }
        #[allow(clippy::cast_precision_loss)]
        let f = i as f64 / n;
        widest = widest.max(gap_at(f, sample)?);
    }

    // Each summit, found by a golden-section search on the gap itself.
    let mut by_height: Vec<usize> = (0..samples.len()).collect();
    by_height.sort_by(|&a, &b| across[b].total_cmp(&across[a]));
    let tallest = across[by_height[0]];
    let mut peaks: Vec<usize> = by_height
        .into_iter()
        .filter(|&i| {
            across[i] > tol.confusion()
                && (i == 0 || across[i] >= across[i - 1])
                && (i + 1 == samples.len() || across[i] >= across[i + 1])
        })
        .take(PEAKS)
        .enumerate()
        .take_while(|&(k, i)| k < 3 || across[i] >= tallest * NEAR_PEAK)
        .map(|(_, i)| i)
        .collect();
    peaks.extend([0, samples.len() - 1]);
    for i in peaks {
        #[allow(clippy::cast_precision_loss)]
        let (mut a, mut b) = (
            i.saturating_sub(1) as f64 / n,
            (i + 1).min(samples.len() - 1) as f64 / n,
        );
        // Climbed on the height the peaks were picked by, which costs no
        // search along the curve, and the gap stated read at the top.
        let height = |f: f64| -> OgeomResult<f64> {
            paced_at(f)?.map_or(Ok(0.0), |sample| across_at(f, sample))
        };
        let (mut c, mut d) = (b - (b - a) * ratio, a + (b - a) * ratio);
        let (mut hc, mut hd) = (height(c)?, height(d)?);
        for _ in 0..20 {
            if hc > hd {
                (b, d, hd) = (d, c, hc);
                c = b - (b - a) * ratio;
                hc = height(c)?;
            } else {
                (a, c, hc) = (c, d, hd);
                d = a + (b - a) * ratio;
                hd = height(d)?;
            }
        }
        let top = if hc > hd { c } else { d };
        if let Some(sample) = paced_at(top)? {
            widest = widest.max(gap_at(top, sample)?);
        }
    }
    Ok(widest)
}

/// Re-project the samples a stalled projection left behind.
///
/// Where a chart collapses (a spline patch whose whole `v = 0` row is a
/// single point), the projector has no direction to move in, and it answers
/// with the pole's own parameters and the distance to it. A run of samples
/// pinned to such a row would be taken for the file's own boundary slop,
/// the edge widened to cover it, and a curve drawn through it. A sample
/// that landed badly is retried from a neighbour that landed well, the
/// same seeding the forward walk already trusts, run in both directions so
/// a run of them unwinds from whichever end is sound. The retry is kept
/// only when it lands closer, so it can never make an honest projection
/// worse: a file's real slop is left alone.
fn retry_stalled(
    surface: &SurfaceGeometry,
    points: &[Point],
    trace: &mut [ogeom_math::Point2],
    offs: &mut [f64],
    tol: Tolerances,
) {
    let sound = tol.confusion() * 1e5;
    if offs.iter().all(|off| *off <= sound) {
        return;
    }
    for backwards in [true, false] {
        let order: Vec<usize> = if backwards {
            (0..offs.len()).rev().collect()
        } else {
            (0..offs.len()).collect()
        };
        for i in order {
            if offs[i] <= sound {
                continue;
            }
            let Some(j) = (if backwards {
                i.checked_add(1)
            } else {
                i.checked_sub(1)
            }) else {
                continue;
            };
            if offs.get(j).is_none_or(|off| *off > sound) {
                continue;
            }
            let seed = (trace[j].x, trace[j].y);
            if let Ok(found) =
                crate::measure::project_on_surface_from(surface, points[i], seed, tol)
                && found.distance < offs[i]
            {
                trace[i] = ogeom_math::Point2::new(found.parameters.0, found.parameters.1);
                offs[i] = found.distance;
            }
        }
    }
}

/// Slide a trace back into its chart by whole turns.
///
/// Unwrapped for continuity, a trace can end up a whole turn outside the
/// chart it belongs to: a projection that starts near one edge of a closed
/// chart and walks off it keeps walking, and the surface then refuses to be
/// evaluated where its own trim lies, so the face draws as a hole.
///
/// A rigid shift keeps the trace exactly as continuous as the unwrap left
/// it and can only move it inward. One that genuinely spans more than a
/// turn has nowhere to go and is left alone; one that fits nowhere whole
/// takes the turn that centres it, which is the nearest thing to inside
/// there is.
fn slide_into_chart(
    trace: &mut [ogeom_math::Point2],
    spans: (f64, f64),
    domain: ((f64, f64), (f64, f64)),
) {
    let ((ua, ub), (va, vb)) = domain;
    for (across, span, lo, hi) in [(true, spans.0, ua, ub), (false, spans.1, va, vb)] {
        if span <= 0.0 {
            continue;
        }
        let read = |uv: &ogeom_math::Point2| if across { uv.x } else { uv.y };
        let (mut least, mut most) = (f64::INFINITY, f64::NEG_INFINITY);
        for uv in trace.iter() {
            least = least.min(read(uv));
            most = most.max(read(uv));
        }
        if !(least.is_finite() && most.is_finite()) || most - least > span {
            continue;
        }
        let turns = {
            let up = ((lo - least) / span).ceil();
            let down = ((hi - most) / span).floor();
            if up <= down {
                // Somewhere it fits whole; the nearest such turn.
                up.max(down.min(0.0))
            } else {
                (f64::midpoint(lo, hi) - f64::midpoint(least, most)) / span
            }
        };
        let turns = if turns.is_finite() {
            turns.round()
        } else {
            0.0
        };
        if turns == 0.0 {
            continue;
        }
        for uv in trace.iter_mut() {
            if across {
                uv.x += turns * span;
            } else {
                uv.y += turns * span;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use ogeom_math::Point2;

    /// A trace unwrapped clean off its chart is slid back by whole turns.
    ///
    /// A fitted `v` running from −2.5π to −π on a chart that stops at −π
    /// is continuous and outside. The surface refuses to be asked about
    /// it, and the face draws as a hole.
    #[test]
    fn a_trace_that_walked_off_its_chart_is_slid_back() {
        let pi = core::f64::consts::PI;
        let chart = ((0.0, 1.0), (-pi, pi));
        let turn = (0.0, 2.0 * pi);

        // A whole turn below: slid up, and as continuous as it was.
        let mut trace = vec![
            Point2::new(0.5, -2.5 * pi),
            Point2::new(0.5, -2.0 * pi),
            Point2::new(0.5, -1.5 * pi),
        ];
        super::slide_into_chart(&mut trace, turn, chart);
        assert!(
            trace.iter().all(|at| at.y >= -pi && at.y <= pi),
            "slid into the chart: {trace:?}"
        );
        for pair in trace.windows(2) {
            assert!(
                (pair[1].y - pair[0].y - 0.5 * pi).abs() < 1e-12,
                "and rigidly"
            );
        }

        // Already inside: untouched.
        let mut held = vec![Point2::new(0.5, -1.0), Point2::new(0.5, 1.0)];
        let was = held.clone();
        super::slide_into_chart(&mut held, turn, chart);
        assert_eq!(held, was);

        // Wider than a turn: nowhere to slide to, and left alone.
        let mut wide = vec![Point2::new(0.5, -4.0 * pi), Point2::new(0.5, 0.0)];
        let was = wide.clone();
        super::slide_into_chart(&mut wide, turn, chart);
        assert_eq!(wide, was);

        // An axis that does not close is not slid on at all.
        let mut across = vec![Point2::new(9.0, 0.0), Point2::new(9.5, 0.0)];
        let was = across.clone();
        super::slide_into_chart(&mut across, (0.0, 2.0 * pi), chart);
        assert_eq!(across, was);
    }

    use super::*;
    use ogeom_core::Tolerances;
    use ogeom_geom::{CircleCurve, CylinderSurface};
    use ogeom_math::{Circle, Cylinder, Frame};

    const T: Tolerances = Tolerances::millimetres();

    /// A sample the projector left at a pole is retried from its neighbour.
    ///
    /// The patch is `S(u, v) = v·C(u)`: its whole `v = 0` row is the origin,
    /// so a projection that reaches the pole has no direction left to move
    /// in and stops there, however far off it is. Seeding from the pole is
    /// shown stuck first (that is the trap the forward walk falls into,
    /// once per imported edge that starts on such a row), and the retry from a
    /// sound neighbour is shown to get out of it.
    #[test]
    fn a_sample_stalled_at_a_pole_is_retried_from_its_neighbour() {
        use ogeom_geom::{BSplineSurface, Surface as _};
        use ogeom_math::{ControlGrid, KnotVector};

        let mut control = Vec::new();
        for i in 0..4 {
            let across = -1.0 + 2.0 * f64::from(i) / 3.0;
            for j in 0..4 {
                let out = f64::from(j) / 3.0;
                control.push(Point::new(out, across * out, (0.5 + across * across) * out));
            }
        }
        let grid = ControlGrid::new(control, 4, 4).unwrap();
        let cone: SurfaceGeometry = BSplineSurface::new(
            KnotVector::clamped_uniform(3, 4).unwrap(),
            KnotVector::clamped_uniform(3, 4).unwrap(),
            &grid,
            T,
        )
        .unwrap()
        .into();
        assert!(
            cone.d1_at(0.5, 0.0, T).unwrap().0.magnitude() < 1e-12,
            "the v = 0 row is a pole"
        );

        let at = cone.point_at(1.0, 0.2, T).unwrap();
        let stuck = crate::measure::project_on_surface_from(&cone, at, (0.0, 0.0), T).unwrap();
        assert!(
            stuck.distance > 0.1,
            "a projection seeded at the pole stays there: {stuck:?}"
        );

        let mut trace = vec![
            ogeom_math::Point2::new(0.0, 0.0),
            ogeom_math::Point2::new(0.0, 0.0),
            ogeom_math::Point2::new(1.0, 0.5),
        ];
        let points = vec![
            cone.point_at(1.0, 0.1, T).unwrap(),
            at,
            cone.point_at(1.0, 0.5, T).unwrap(),
        ];
        let mut offs = vec![points[0].distance(Point::ORIGIN), stuck.distance, 0.0];
        retry_stalled(&cone, &points, &mut trace, &mut offs, T);
        for (index, off) in offs.iter().enumerate() {
            assert!(
                *off < T.confusion(),
                "sample {index} found its surface: {off:.3e}"
            );
        }
        assert!(
            (trace[1].x - 1.0).abs() < 1e-6 && (trace[1].y - 0.2).abs() < 1e-6,
            "and its own parameters: {:?}",
            trace[1]
        );

        // An honest miss is not a stall: nothing lands closer, so the
        // file's own slop survives the retry untouched.
        let adrift = Point::new(0.0, 0.0, -0.5);
        let mut honest = vec![trace[2], ogeom_math::Point2::new(1.0, 0.5)];
        let was = honest.clone();
        let mut misses = vec![0.0, 0.5];
        retry_stalled(&cone, &[points[2], adrift], &mut honest, &mut misses, T);
        assert_eq!(honest[1], was[1]);
        assert!((misses[1] - 0.5).abs() < 1e-12);
    }

    /// A direction is weak by what crossing it moves, not by its neighbour.
    ///
    /// The patch is a flat strip: `u` runs a millimetre across it and `v`
    /// runs twenty millimetres along it over a parameter span of a hundredth,
    /// two thousand times denser than `u`. Against `dv`, `du` looks weak at
    /// every sample, and a ratio test holds every `u` at one value: an edge
    /// a millimetre long across the strip fits as a single chart point.
    /// Crossing the whole of `u` moves the point a millimetre, which is the
    /// only thing "weak" can honestly mean, and it is not.
    #[test]
    fn a_direction_is_weak_by_what_crossing_it_moves() {
        use ogeom_geom::{BSplineSurface, LineCurve};
        use ogeom_math::{ControlGrid, KnotVector, Point};
        let mut control = Vec::new();
        for i in 0..2 {
            for j in 0..2 {
                control.push(Point::new(f64::from(i), 20.0 * f64::from(j), 0.0));
            }
        }
        let strip: SurfaceGeometry = BSplineSurface::new(
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 0.01, 0.01], 1).unwrap(),
            &ControlGrid::new(control, 2, 2).unwrap(),
            T,
        )
        .unwrap()
        .into();
        // Across the strip at v = 0.005 (the middle in space, 10 mm along).
        let across: Curve =
            LineCurve::segment(Point::new(0.0, 10.0, 0.0), Point::new(1.0, 10.0, 0.0), T)
                .unwrap()
                .into();
        let (pcurve, error, _, _, _) =
            fit_projected_pcurve(&across, (0.0, 1.0), &strip, T).unwrap();
        let a = ogeom_geom::Curve2d::point_at(&pcurve, 0.0, T).unwrap();
        let b = ogeom_geom::Curve2d::point_at(&pcurve, 1.0, T).unwrap();
        assert!(
            (b.x - a.x).abs() > 0.99,
            "the edge crosses the whole of u: {a:?} -> {b:?}"
        );
        assert!(error < 1e-6, "and fits: {error:.2e}");
    }

    /// A flat patch over the unit square, a spline surface whose chart
    /// inverts only by projection.
    fn flat_patch() -> SurfaceGeometry {
        use ogeom_geom::BSplineSurface;
        use ogeom_math::{ControlGrid, KnotVector, Point};
        let control = [(0.0, -1.0), (0.0, 1.0), (1.0, -1.0), (1.0, 1.0)]
            .map(|(x, y)| Point::new(x, y, 0.0))
            .to_vec();
        BSplineSurface::new(
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap(),
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap(),
            &ControlGrid::new(control, 2, 2).unwrap(),
            T,
        )
        .unwrap()
        .into()
    }

    /// A pcurve's error covers the points between its samples on a patch.
    ///
    /// The curve runs across the patch in a sine of a hundredth of a
    /// millimetre whose period is the samples' spacing: every sample and
    /// every middle between two stands on the sine's axis, and the fit
    /// through them is the axis. The quarter points stand on its crests,
    /// and the error stated is the crest's height.
    #[test]
    fn a_fit_on_a_patch_states_its_widest_miss_between_the_samples() {
        use ogeom_geom::{Curve2d as _, Curve3d as _, Surface as _};
        use ogeom_math::Point;
        let amplitude = 0.01;
        let sine = |t: f64| {
            Point::new(
                t,
                amplitude * (core::f64::consts::TAU * 96.0 * t).sin(),
                0.0,
            )
        };
        let ts: Vec<f64> = (0..=3840).map(|i| f64::from(i) / 3840.0).collect();
        let points: Vec<Point> = ts.iter().map(|t| sine(*t)).collect();
        let fitted = ogeom_geom::fit::fit_points_at(&ts, &points, 3, 1e-8, T).unwrap();
        let curve: Curve = fitted.curve.into();
        let patch = flat_patch();
        let (pcurve, error, met, _, _) =
            fit_projected_pcurve(&curve, (0.0, 1.0), &patch, T).unwrap();
        let mut worst = 0.0_f64;
        for i in 0..=20_000 {
            let t = f64::from(i) / 20_000.0;
            let uv = pcurve.point_at(t, T).unwrap();
            let lifted = patch.point_at(uv.x, uv.y, T).unwrap();
            worst = worst.max(lifted.distance(curve.point_at(t, T).unwrap()));
        }
        assert!(
            worst > 0.9 * amplitude,
            "the fit misses the crests: {worst:.3e}"
        );
        assert!(
            error >= worst,
            "the widest miss is stated: {error:.3e} for {worst:.3e}"
        );
        assert!(!met, "and the target is not met");
    }

    /// The gap stated where a pcurve is paced differently from its curve
    /// is found at its widest. The pcurve runs along a straight line of the
    /// patch a little ahead of it, and strays across it in bumps about as
    /// narrow as the measure's samples are apart, the tallest three
    /// quarters along. The paced distance is the lead's, which peaks
    /// halfway; the samples read the tallest bump a few percent low, and
    /// only a climb from a peak of the distance across reaches its top.
    #[test]
    fn a_pcurve_paced_differently_is_measured_at_its_widest_bump() {
        use ogeom_geom::{Curve2d as _, LineCurve};
        use ogeom_math::{Point, Point2};
        let edge: Curve =
            LineCurve::segment(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), T)
                .unwrap()
                .into();
        let patch = flat_patch();
        // u = t plus a slow lead; v (half of y on this patch) carries the
        // bumps, tallest three quarters along.
        let lead = |t: f64| 0.002 * (core::f64::consts::PI * t).sin();
        let across = |t: f64| {
            let bump = (core::f64::consts::TAU * 37.3 * t).sin().max(0.0).powi(8);
            1e-5 * (1.0 + 0.5 * (-((t - 0.75) / 0.05).powi(2)).exp()) * bump
        };
        let ts: Vec<f64> = (0..=4000).map(|i| f64::from(i) / 4000.0).collect();
        let chart: Vec<Point2> = ts
            .iter()
            .map(|&t| Point2::new(t + lead(t), 0.5 + 0.5 * across(t)))
            .collect();
        let pcurve: PlanarCurve = ogeom_geom::fit::fit_points_2d_at(&ts, &chart, 3, 1e-10, T)
            .unwrap()
            .curve
            .into();
        let mut truth = 0.0_f64;
        for i in 0..=40_000 {
            let t = f64::from(i) / 40_000.0;
            let uv = pcurve.point_at(t, T).unwrap();
            let lifted = patch.point_at(uv.x, uv.y, T).unwrap();
            // The edge is the line y = 0: the nearest-point gap is |y|.
            truth = truth.max(lifted.y.abs());
        }
        let gap = lifted_gap((&edge, (0.0, 1.0)), (&pcurve, (0.0, 1.0)), &patch, false, T).unwrap();
        assert!(
            gap >= 0.99 * truth && gap <= truth * 1.001 + 1e-12,
            "the widest bump, {truth:.4e}, is stated: {gap:.4e}"
        );
    }

    /// A boundary 0.3 mm off its surface fits, and says so.
    ///
    /// Exchange files carry boundary curves that far from the surfaces
    /// they trim, up to a third of a millimetre. The fit
    /// accepts anything under a millimetre and reports the offset, so the
    /// reader widens the edge's tolerance instead of leaving the face
    /// without a trim; a miss of whole millimetres (the signature of an
    /// edge paired with the wrong surface) still refuses.
    #[test]
    fn slop_under_a_millimetre_fits_and_is_reported() {
        let wall: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 10.0, T).unwrap(), (-50.0, 50.0))
                .unwrap()
                .into();
        let rim = |radius: f64| -> Curve {
            CircleCurve::new(Circle::new(Frame::WORLD, radius, T).unwrap()).into()
        };
        // 0.3 mm proud of the wall: every sample sits exactly that far off.
        let (_, _, _, worst_off, warning) =
            fit_projected_pcurve(&rim(10.3), (0.0, core::f64::consts::TAU), &wall, T).unwrap();
        assert!(
            (worst_off - 0.3).abs() < 1e-6,
            "the offset is measured: {worst_off}"
        );
        assert!(warning.is_some(), "slop this large is worth a warning");
        // 3 mm off is not slop; it is the wrong surface.
        assert!(
            fit_projected_pcurve(&rim(13.0), (0.0, core::f64::consts::TAU), &wall, T).is_err(),
            "a miss of millimetres still refuses"
        );
    }

    /// A fit is held between its samples, where a fast-running curve bends.
    ///
    /// The samples are spaced evenly in the curve's parameter, and this
    /// curve spends a twentieth of its parameter on a steep drop of a
    /// millimetre before running slowly round the drum for the rest: the
    /// whole of the drop, and the bend at its foot, fall inside the first
    /// sample interval. A cubic held only at the samples hooks past the
    /// curve there. Measured between the samples and refitted where it
    /// leaves them, the pcurve follows the curve everywhere.
    #[test]
    fn a_fit_is_held_between_its_samples_where_the_curve_runs_fast() {
        use ogeom_geom::{BSplineCurve, Curve, Curve2d, Curve3d, Surface, SurfaceGeometry};
        use ogeom_math::{KnotVector, Point};
        let radius = 10.0;
        let on = |angle: f64, z: f64| Point::new(radius * angle.cos(), radius * angle.sin(), z);
        let mut knots = vec![0.0, 0.0, 0.0, 0.0, 0.05];
        knots.extend((1..8).map(f64::from));
        knots.extend([8.0; 4]);
        let control = vec![
            on(0.12, 4.0),
            on(0.13, 3.7),
            on(0.14, 3.4),
            on(0.15, 3.1),
            on(0.10, 2.5),
            on(0.0, 1.5),
            on(-0.1, 0.5),
            on(-0.2, -0.5),
            on(-0.3, -1.5),
            on(-0.35, -2.5),
            on(-0.38, -3.3),
            on(-0.4, -4.0),
        ];
        let curve: Curve = BSplineCurve::new(KnotVector::new(knots, 3).unwrap(), control, T)
            .unwrap()
            .into();
        let drum: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(Frame::WORLD, radius, T).unwrap(), (-5.0, 5.0))
                .unwrap()
                .into();
        let (pcurve, error, _, _, _) = fit_projected_pcurve(&curve, (0.0, 8.0), &drum, T).unwrap();
        let mut worst = 0.0_f64;
        for i in 0..=2000 {
            let t = 8.0 * f64::from(i) / 2000.0;
            let p = curve.point_at(t, T).unwrap();
            let uv = pcurve.point_at(t, T).unwrap();
            let lifted = drum.point_at(uv.x, uv.y, T).unwrap();
            // The pcurve's own miss: the curve's point, pulled onto the drum.
            worst = worst.max(lifted.distance(on(p.y.atan2(p.x), p.z)));
        }
        assert!(
            worst <= T.confusion() * 3e4,
            "the pcurve leaves the curve by {worst:.3e} between the samples, \
             {error:.3e} reported"
        );
        assert!(
            error >= worst,
            "and the widest miss between the samples is reported: {error:.3e} for {worst:.3e}"
        );
    }
}
