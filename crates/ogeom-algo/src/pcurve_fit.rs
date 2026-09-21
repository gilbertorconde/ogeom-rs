//! Pcurves fitted by projection, for faces whose trims have no closed form.
//!
//! An exchange file's edge must end up with a curve in each bounding face's
//! parameters, or the face cannot be split or triangulated. Where the
//! curve/surface pair has a closed form the exact projection is used; where
//! it does not — a spline surface, mostly — the pcurve is *fitted at the
//! curve's own parameters*: sample the edge, project each sample into the
//! chart, fit the trace with the parameters held fixed, so the same-parameter
//! law holds by construction. This honours the standing decision that an
//! exact curve never carries a fitted pcurve silently: the fit's error is
//! returned, and the callers widen tolerances and warn with it. That error
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
/// The reader draws its line at a millimetre — below it is a file's own
/// slop, above it a wrong pairing — but a *healer* acts on instruction, and
/// the instruction carries the cap. Returns the fitted pcurve, the fit's
/// reached error as a length, whether it met its target, the worst measured
/// edge-to-surface offset, and the slop note when that offset is worth
/// saying out loud.
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
        let (uv, off) = match chart_of(surface, p) {
            // Analytic surfaces invert in closed form — grid seeding
            // over a plane's or cylinder's enormous stated extents lands
            // microns off, and a fitted pcurve inherits every micron.
            Some(uv) => {
                let lifted = surface.point_at(uv.x, uv.y, tol)?;
                (uv, p.distance(lifted))
            }
            None => {
                // Where the previous sample landed is a far better starting
                // guess than any grid: consecutive samples of a curve are
                // neighbouring points of the surface. Trusted only when it
                // lands convincingly *on* the surface — the same bar the
                // denser reseed below is judged against — so a guess that
                // wandered into the wrong basin, or the first sample, which
                // has no predecessor, still pays for the grid.
                let near = previous.and_then(|(_, luv)| {
                    crate::measure::project_on_surface_from(surface, p, (luv.x, luv.y), tol).ok()
                });
                if let Some(close) = near.filter(|f| f.distance <= tol.confusion() * 1e5) {
                    (
                        ogeom_math::Point2::new(close.parameters.0, close.parameters.1),
                        close.distance,
                    )
                } else {
                    let mut projection = crate::measure::project_on_surface(surface, p, 24, tol)?;
                    if projection.distance > tol.confusion() * 1e5 {
                        // A miss this large on a spline surface is more often
                        // a projection stuck in the wrong basin than real
                        // slop; seed denser before believing it.
                        let denser = crate::measure::project_on_surface(surface, p, 96, tol)?;
                        if denser.distance < projection.distance {
                            projection = denser;
                        }
                    }
                    (
                        ogeom_math::Point2::new(projection.parameters.0, projection.parameters.1),
                        projection.distance,
                    )
                }
            }
        };
        previous = Some((p, uv));
        parameters.push(t);
        trace.push(uv);
        points.push(p);
        offs.push(off);
    }
    retry_stalled(surface, &points, &mut trace, &mut offs, tol);
    // The cap separates a file's own slop from an edge paired with the
    // wrong surface. Slop is routinely a micron or two, but real community
    // exports carry as much as 0.34 mm (issue #15) — while a wrong pairing
    // misses by the distance between two different surfaces of the body,
    // whole millimetres. One millimetre stands between the worst slop
    // observed and the smallest wrong pairing plausible. Slop inside the
    // cap is accepted and *recorded*: the edge's tolerance is widened to
    // cover it, so the model says what it knows instead of refusing to
    // triangulate. Judged after the retry, so a projection that stalled
    // cannot refuse an edge that is actually on its surface.
    let worst_off = offs.iter().copied().fold(0.0_f64, f64::max);
    if worst_off > cap {
        ogeom_bail!(
            Construction,
            "the edge sits {worst_off:.2e} from the surface it should bound"
        );
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
    // projections near the joining column land in either copy — both
    // answers are right pointwise, and only continuity chooses. This is
    // the docs/PLAN.md F5 case, and it is decided here for both exchange
    // readers at once.
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
    // Where the chart collapses — a sphere's pole, a cone's apex — the
    // u of a sample is atan2 of noise: the point determines no angle.
    // The *arc* does: a smooth curve through the pole approaches it at
    // a definite chart angle, which is the limit of its well-conditioned
    // neighbours. Samples whose u-direction has collapsed relative to
    // their v-direction are repaired by interpolating u between the
    // nearest sound samples, extrapolating at the ends.
    //
    // Weak is measured in millimetres, not against `dv`. A ratio calls a
    // direction weak whenever the *other* one is strong, and a patch whose
    // `v` is parameterised a thousand times more densely than its `u` —
    // three millimetres over twelve thousandths of a unit, beside a unit of
    // `u` for a little over half a millimetre — had every sample of every
    // edge called weak, its `u` held at one value, and two edges half a
    // millimetre long fitted as a single point. What makes a direction
    // degenerate is that crossing the whole of it moves the point less than
    // a micron; that question has an answer in length, and only in length.
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
    // metric — the honest cheap version, refined by the fit's report.
    let scale = if space_run > tol.confusion() {
        parameter_run / space_run
    } else {
        1.0
    };
    let target = (tol.confusion() * 1e2 * scale).max(f64::MIN_POSITIVE);
    let fitted = ogeom_geom::fit::fit_points_2d_at(&parameters, &trace, 3, target, tol)?;
    // A least-squares fit wiggles past its samples at the ends, and a
    // surface with a *tight* stated window — an imported patch, not a
    // reader-built analytic with its enormous extents — refuses evaluation
    // a hair outside it. The control points clamp into the window on the
    // non-periodic axes: the curve lives in its controls' hull, so the
    // clamp is a guarantee. Whatever the clamp cost is not hidden either —
    // it is the clamped curve that is measured below. Periodic axes stay
    // free: an unwrapped trace crosses the seam on purpose.
    let fitted = {
        let ((wa, wb), (va2, vb2)) = surface.domain();
        let clamp_u = !(surface.is_periodic_u() || surface.is_closed_u(tol));
        let clamp_v = !(surface.is_periodic_v() || surface.is_closed_v(tol));
        if clamp_u || clamp_v {
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
                ogeom_geom::fit::Fitted {
                    curve: ogeom_geom::BSpline2d::new(knots, control, tol)?,
                    ..fitted
                }
            } else {
                fitted
            }
        } else {
            fitted
        }
    };
    // What the caller is told, as a length. The fitter reports its error
    // in *chart* units, and a chart's units are whatever the file chose: one
    // patch met in the wild spans four microns across its `u` and ten
    // millimetres along its `v`, so no single scale converts the one number
    // into the other — a control point dragged back into that chart by the
    // clamp read as seven hundred millimetres of mesh error, on a face a
    // tenth of a millimetre across. So the fitted curve is walked instead,
    // through the surface, against the trace it was fitted to: that
    // difference is the fit's own and is measured where the mesh will be.
    let mut error = 0.0_f64;
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
    }
    let met = error <= tol.confusion() * 1e2;
    let slop = (worst_off > tol.confusion() * 1e3).then(|| {
        format!(
            "an edge sits up to {worst_off:.2e} from the surface it \
             bounds; the file's own slop, carried into the chart"
        )
    });
    Ok((fitted.curve.into(), error, met, worst_off, slop))
}

/// Re-project the samples a stalled projection left behind.
///
/// Where a chart collapses — a spline patch whose whole `v = 0` row is a
/// single point — the projector has no direction to move in, and it answers
/// with the pole's own parameters and the distance to it. Four consecutive
/// samples of one imported edge came back pinned to such a row, the last of
/// them a tenth of a millimetre out; the reader repeated that as the file's
/// own boundary slop, widened the edge to cover it, and the fitter tried to
/// draw a curve through it. A sample that landed badly is retried from a
/// neighbour that landed well — the same seeding the forward walk already
/// trusts, run in both directions so a run of them unwinds from whichever
/// end is sound. The retry is kept only when it lands closer, so it can
/// never make an honest projection worse: a file's real slop is left alone.
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
/// evaluated where its own trim lies — an imported face whose
/// fitted v ran to −2.5π on a chart that stops at −π, and drew as a hole.
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
    /// A real assembly has a face whose fitted `v` ran from −2.5π to
    /// −π on a chart that stops at −π: continuous, outside, and the surface
    /// refuses to be asked about it, so the face drew as a hole.
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
    /// shown stuck first — that is the trap the forward walk falls into,
    /// once per imported edge that starts on such a row — and the retry from a
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
    /// runs twenty millimetres along it over a parameter span of a hundredth
    /// — two thousand times denser than `u`. Against `dv`, `du` looks weak at
    /// every sample, and a ratio test held every `u` at one value: an edge
    /// a millimetre long across the strip fitted as a single chart point.
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

    /// A boundary 0.3 mm off its surface fits, and says so.
    ///
    /// Community exports carry boundary curves that far from the surfaces
    /// they trim (issue #15: 0.12–0.34 mm on a real assembly). The fit
    /// accepts anything under a millimetre and reports the offset, so the
    /// reader widens the edge's tolerance instead of leaving the face
    /// without a trim; a miss of whole millimetres — the signature of an
    /// edge paired with the wrong surface — still refuses.
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
}
