//! Pcurves for imported faces, shared by the exchange readers.
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
pub(crate) type FittedPcurve = OgeomResult<(PlanarCurve, f64, bool, f64, Option<String>)>;

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

pub(crate) fn fit_projected_pcurve(
    curve: &Curve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
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
                    ogeom_algo::project_on_surface_from(surface, p, (luv.x, luv.y), tol).ok()
                });
                if let Some(close) = near.filter(|f| f.distance <= tol.confusion() * 1e5) {
                    (
                        ogeom_math::Point2::new(close.parameters.0, close.parameters.1),
                        close.distance,
                    )
                } else {
                    let mut projection = ogeom_algo::project_on_surface(surface, p, 24, tol)?;
                    if projection.distance > tol.confusion() * 1e5 {
                        // A miss this large on a spline surface is more often
                        // a projection stuck in the wrong basin than real
                        // slop; seed denser before believing it.
                        let denser = ogeom_algo::project_on_surface(surface, p, 96, tol)?;
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
        // The cap separates a file's own slop — routinely a micron or
        // two on these parts — from an edge paired with the wrong
        // surface, which misses by whole millimetres. Slop inside the
        // cap is accepted and *recorded*: the edge's tolerance is
        // widened to cover it, so the model says what it knows instead
        // of refusing to triangulate.
        if off > tol.confusion() * 1e6 {
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
