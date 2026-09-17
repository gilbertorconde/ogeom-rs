//! Pcurves fitted by projection, for faces whose trims have no closed form.
//!
//! An exchange file's edge must end up with a curve in each bounding face's
//! parameters, or the face cannot be split or triangulated. Where the
//! curve/surface pair has a closed form the exact projection is used; where
//! it does not — a spline surface, mostly — the pcurve is *fitted at the
//! curve's own parameters*: sample the edge, project each sample into the
//! chart, fit the trace with the parameters held fixed, so the same-parameter
//! law holds by construction and the reported error is the true chart
//! deviation. This honours the standing decision that an exact curve never
//! carries a fitted pcurve silently: the fit's error is returned, and the
//! callers widen tolerances and warn with it.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::Surface as _;
use ogeom_geom::{Curve, PlanarCurve, SurfaceGeometry};
use ogeom_math::Point;

/// The fitted pcurve, its fit error, whether the target was met, the worst
/// distance any sample sat from the surface, and the slop warning to record
/// when that distance is large enough to say out loud.
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
/// reached error, whether it met its target, the worst measured
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
    let mut worst_off = 0.0_f64;
    let mut parameters = Vec::with_capacity(SAMPLES + 1);
    let mut trace = Vec::with_capacity(SAMPLES + 1);
    let mut space_run = 0.0;
    let mut parameter_run = 0.0;
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
        // The cap separates a file's own slop from an edge paired with
        // the wrong surface. Slop is routinely a micron or two, but real
        // community exports carry as much as 0.34 mm (issue #15) — while a
        // wrong pairing misses by the distance between two different
        // surfaces of the body, whole millimetres. One millimetre stands
        // between the worst slop observed and the smallest wrong pairing
        // plausible. Slop inside the cap is accepted and *recorded*: the
        // edge's tolerance is widened to cover it, so the model says what
        // it knows instead of refusing to triangulate.
        if off > cap {
            ogeom_bail!(
                Construction,
                "the edge sits {off:.2e} from the surface it should bound"
            );
        }
        worst_off = worst_off.max(off);
        if let Some((lp, luv)) = previous {
            space_run += p.distance(lp);
            parameter_run += uv.distance(luv);
        }
        previous = Some((p, uv));
        parameters.push(t);
        trace.push(uv);
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
    let weak: Vec<bool> = trace
        .iter()
        .map(|uv| {
            surface
                .d1_at(uv.x, uv.y, tol)
                .is_ok_and(|(du, dv)| du.magnitude() < dv.magnitude() * 1e-3)
        })
        .collect();
    if weak.iter().any(|w| *w) && weak.iter().filter(|w| !**w).count() >= 2 {
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
    // clamp is a guarantee, and the distance it moved joins the reported
    // error instead of being hidden. Periodic axes stay free — an unwrapped
    // trace crosses the seam on purpose.
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
                    error: fitted.error + moved,
                    met: fitted.met && fitted.error + moved <= target,
                }
            } else {
                fitted
            }
        } else {
            fitted
        }
    };
    let slop = (worst_off > tol.confusion() * 1e3).then(|| {
        format!(
            "an edge sits up to {worst_off:.2e} from the surface it \
             bounds; the file's own slop, carried into the chart"
        )
    });
    Ok((
        fitted.curve.into(),
        fitted.error,
        fitted.met,
        worst_off,
        slop,
    ))
}

/// Slide a trace back into its chart by whole turns.
///
/// Unwrapped for continuity, a trace can end up a whole turn outside the
/// chart it belongs to: a projection that starts near one edge of a closed
/// chart and walks off it keeps walking, and the surface then refuses to be
/// evaluated where its own trim lies — a face of the Voron assembly whose
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
    /// The Voron assembly has a face whose fitted `v` ran from −2.5π to
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
