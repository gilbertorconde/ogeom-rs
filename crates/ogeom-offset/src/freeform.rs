//! The parallel of a free-form face of a solid, for the offset rebuild.
//!
//! A face on a B-spline, revolution, extrusion, trimmed or offset surface
//! has no parallel in its own family, so it moves to a fit: the surface's
//! points moved the distance along its normal, fitted at the chart's own
//! parameters as a sheet's offset is (see [`crate::thick`]), so the chart
//! and every pcurve's meaning carry over. Where the chart closes on itself
//! in a point (a pole), the normal there is the limit of the normals
//! approaching it, read from the second derivatives.
//!
//! The rebuild re-derives each edge where the moved faces meet, and a
//! moved neighbour meets this face past its old boundary wherever the
//! corner between them opens (a convex corner grown, a concave one
//! shrunk). The fit therefore covers the face's own window and runs on
//! past it: inside the surface's chart by the same reach the analytic
//! families are opened by, and past a side of the chart the face reaches
//! by continuing the surface there as itself (a B-spline's end spans
//! continued as their own polynomials, an extrusion along its direction
//! and along its profile's continuation). Only the part over the face's
//! own window is kept by the trim; the continuation says where the corner
//! falls.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, explore};

/// A chart window, `((u0, u1), (v0, v1))`.
type Window = ((f64, f64), (f64, f64));

/// How many samples a pcurve or a side of a window is read at.
const SAMPLES: u32 = 24;

/// The parallel of `face`'s free-form `surface` at `amount` along the
/// face's outward side (`sign` relating it to the surface's own normal),
/// opened by `reach` in space past the face so moved neighbours meet it.
/// Returns the fit and its measured deviation.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// offset reaches past the face's radius of curvature, the surface has no
/// normal somewhere the fit needs one, or the face covers no region of its
/// chart; [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the
/// finest fit misses its bound.
pub(crate) fn moved_free_form(
    model: &Model,
    face: &Shape,
    surface: &SurfaceGeometry,
    (sign, amount): (f64, f64),
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<(SurfaceGeometry, f64)> {
    let (base, carried) = match surface {
        SurfaceGeometry::Offset(o) => (o.basis(), o.distance()),
        other => (other, 0.0),
    };
    let base = match base {
        SurfaceGeometry::Trimmed(t) => t.basis(),
        other => other,
    };
    let total = carried + sign * amount;
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        ogeom_bail!(Construction, "face node holds no face data");
    };
    let own = face_window(model, face, surface, data.natural_restriction, tol)?;
    crate::thick::folds_within(base, total, amount, own, tol)?;

    let (opened, window) = opened(base, own, reach, tol)?;
    let moved_at = |u: f64, v: f64| -> OgeomResult<Option<Point>> {
        let p = opened.point_at(u, v, tol)?;
        Ok(limit_normal(&opened, u, v, total, tol)?.map(|n| p + n * total))
    };
    let target = tol.approximation();
    let bound = crate::thick::fit_bound(data.tolerance.get(), amount, tol);
    let (fitted, _, worst) = crate::thick::fitted_parallel(
        &opened,
        &opened,
        &moved_at,
        window,
        amount,
        (target, bound),
        tol,
    )?;
    Ok((fitted, worst))
}

/// The unit normal of `surface` at `(u, v)`, or where the chart is
/// degenerate there (a pole: one tangent vanishes along a whole side) the
/// limit of the normals approaching it from inside the chart. That limit
/// is read from the mixed derivative: beside a side `v = v0` where `du`
/// vanishes, `du ≈ (v - v0) duv`, so the normal tends to `duv × dv` turned
/// toward the chart's inside. It is kept only where it moves the pole the
/// same way, within the confusion distance over `along`, all across the
/// collapsed side; `None` where it does not, or no limit is read.
pub(crate) fn limit_normal(
    surface: &SurfaceGeometry,
    u: f64,
    v: f64,
    along: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Vector>> {
    if let Ok(n) = surface.normal_at(u, v, tol) {
        return Ok(Some(n.vector()));
    }
    let ((ua, ub), (va, vb)) = surface.domain();
    let read = |u: f64, v: f64| -> OgeomResult<Option<(Vector, bool)>> {
        let jet = surface.jet_at(u, v, tol)?;
        let (su, sv) = (jet.du.magnitude(), jet.dv.magnitude());
        let (raw, on_v_side) = if su <= sv {
            // `du` vanishes along the row `v`: the inside is toward the
            // domain's middle in `v`.
            let toward = if v <= f64::midpoint(va, vb) {
                1.0
            } else {
                -1.0
            };
            (jet.duv.cross(jet.dv) * toward, true)
        } else {
            let toward = if u <= f64::midpoint(ua, ub) {
                1.0
            } else {
                -1.0
            };
            (jet.du.cross(jet.duv) * toward, false)
        };
        let m = raw.magnitude();
        if m <= tol.confusion() {
            return Ok(None);
        }
        Ok(Some((raw / m, on_v_side)))
    };
    let Some((n, on_v_side)) = read(u, v)? else {
        return Ok(None);
    };
    let slack = tol.confusion() / along.abs().max(tol.confusion());
    for k in 0..=SAMPLES {
        let t = f64::from(k) / f64::from(SAMPLES);
        let (pu, pv) = if on_v_side {
            (ua + (ub - ua) * t, v)
        } else {
            (u, va + (vb - va) * t)
        };
        match read(pu, pv)? {
            Some((m, _)) if (m - n).magnitude() <= slack => {}
            _ => return Ok(None),
        }
    }
    Ok(Some(n))
}

/// The region of `surface`'s chart `face` covers: the box round its
/// pcurves on that surface, or round its edges projected where an edge
/// carries none, with a fiftieth of its size spare and capped to the
/// chart's domain where the chart does not repeat; the whole domain for a
/// face bounded by its surface's own extents.
fn face_window(
    model: &Model,
    face: &Shape,
    surface: &SurfaceGeometry,
    natural: bool,
    tol: Tolerances,
) -> OgeomResult<Window> {
    use ogeom_geom::{Curve2d as _, Curve3d as _};
    let domain = surface.domain();
    if natural {
        return Ok(domain);
    }
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        ogeom_bail!(Construction, "face node holds no face data");
    };
    let own = data.surface;
    let mut lo = (f64::INFINITY, f64::INFINITY);
    let mut hi = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut grow = |u: f64, v: f64| {
        lo = (lo.0.min(u), lo.1.min(v));
        hi = (hi.0.max(u), hi.1.max(v));
    };
    for edge in explore(model, face, Filter::OfType(ShapeType::Edge))? {
        let Some(edge_data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            ogeom_bail!(Construction, "edge node holds no edge data");
        };
        let mut charted = false;
        for repr in &edge_data.representations {
            let (curves, range) = match repr {
                EdgeRepr::PCurve {
                    curve,
                    surface,
                    range,
                    ..
                } if *surface == own => (vec![*curve], *range),
                EdgeRepr::Seam {
                    forward,
                    reversed,
                    surface,
                    range,
                    ..
                } if *surface == own => (vec![*forward, *reversed], *range),
                _ => continue,
            };
            for id in curves {
                let Some(pcurve) = model.geometry().pcurve(id) else {
                    ogeom_bail!(Dangling, "a pcurve is not in this model");
                };
                for k in 0..=SAMPLES {
                    let t = range.0 + (range.1 - range.0) * f64::from(k) / f64::from(SAMPLES);
                    let p = pcurve.point_at(t, tol)?;
                    grow(p.x, p.y);
                }
            }
            charted = true;
        }
        if charted || edge_data.degenerate {
            continue;
        }
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = edge_data.curve3d() else {
            continue;
        };
        let Some(curve) = model.geometry().curve(*curve) else {
            ogeom_bail!(Dangling, "curve is not in this model");
        };
        for k in 0..=SAMPLES {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / f64::from(SAMPLES);
            let found = ogeom_algo::project_on_surface(surface, curve.point_at(t, tol)?, 32, tol)?;
            grow(found.parameters.0, found.parameters.1);
        }
    }
    if ![lo.0, lo.1, hi.0, hi.1].iter().all(|x| x.is_finite()) {
        return Ok(domain);
    }
    let ((ua, ub), (va, vb)) = domain;
    let mu = (hi.0 - lo.0).max(tol.parametric()) * 0.02;
    let mv = (hi.1 - lo.1).max(tol.parametric()) * 0.02;
    let (mut u0, mut u1, mut v0, mut v1) = (lo.0 - mu, hi.0 + mu, lo.1 - mv, hi.1 + mv);
    if !surface.is_periodic_u() {
        u0 = u0.max(ua);
        u1 = u1.min(ub);
    }
    if !surface.is_periodic_v() {
        v0 = v0.max(va);
        v1 = v1.min(vb);
    }
    if u1 <= u0 || v1 <= v0 {
        ogeom_bail!(
            Construction,
            "a face covers no finite region of its surface's chart"
        );
    }
    Ok(((u0, u1), (v0, v1)))
}

/// `base` continued past each side of its chart that the face's window
/// `own` reaches, by `reach` in space, and the window opened by `reach`
/// on every side: inside the chart where the face stops short of its
/// side, onto the continuation where it reaches it. A side that is one
/// point (a pole), or that closes onto the opposite side, is not opened,
/// and a side the surface cannot be continued past stays at the chart's
/// side.
fn opened(
    base: &SurfaceGeometry,
    own: Window,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<(SurfaceGeometry, Window)> {
    let mut surface = base.clone();
    let mut window = own;
    for (along_u, at_end) in [(true, false), (true, true), (false, false), (false, true)] {
        let ((ua, ub), (va, vb)) = surface.domain();
        let periodic = if along_u {
            surface.is_periodic_u() || surface.is_closed_u(tol)
        } else {
            surface.is_periodic_v() || surface.is_closed_v(tol)
        };
        let ((u0, u1), (v0, v1)) = window;
        let (side, edge) = match (along_u, at_end) {
            (true, false) => (u0, ua),
            (true, true) => (u1, ub),
            (false, false) => (v0, va),
            (false, true) => (v1, vb),
        };
        if periodic {
            continue;
        }
        // The side's mean speed across it turns the reach into a parameter
        // span; a side whose points stand still along it is a pole.
        let (mut speed, mut sweep) = (0.0, 0.0);
        for k in 0..=SAMPLES {
            let t = f64::from(k) / f64::from(SAMPLES);
            let (u, v) = if along_u {
                (side, v0 + (v1 - v0) * t)
            } else {
                (u0 + (u1 - u0) * t, side)
            };
            let (du, dv) = surface.d1_at(u, v, tol)?;
            let (across, along) = if along_u { (du, dv) } else { (dv, du) };
            speed += across.magnitude();
            sweep += along.magnitude();
        }
        speed /= f64::from(SAMPLES + 1);
        sweep /= f64::from(SAMPLES + 1);
        if speed <= tol.confusion() || sweep <= tol.confusion() {
            continue;
        }
        let span = reach / speed;
        let reaches_edge = (side - edge).abs() <= tol.parametric();
        let opened_side = if reaches_edge {
            match continued(&surface, along_u, at_end, reach, tol)? {
                Some(longer) => {
                    surface = longer;
                    let ((na, nb), (nva, nvb)) = surface.domain();
                    match (along_u, at_end) {
                        (true, false) => na,
                        (true, true) => nb,
                        (false, false) => nva,
                        (false, true) => nvb,
                    }
                }
                None => edge,
            }
        } else if at_end {
            (side + span).min(edge)
        } else {
            (side - span).max(edge)
        };
        window = match (along_u, at_end) {
            (true, false) => ((opened_side, u1), (v0, v1)),
            (true, true) => ((u0, opened_side), (v0, v1)),
            (false, false) => ((u0, u1), (opened_side, v1)),
            (false, true) => ((u0, u1), (v0, opened_side)),
        };
    }
    Ok((surface, window))
}

/// `surface` continued past one side of its chart by `length` in space,
/// its chart kept where it was: a B-spline's end spans continued as their
/// own polynomials, an extrusion along its direction or along its
/// B-spline profile's continuation, a revolution along its profile's;
/// `None` for a side no such continuation exists for.
fn continued(
    surface: &SurfaceGeometry,
    along_u: bool,
    at_end: bool,
    length: f64,
    tol: Tolerances,
) -> OgeomResult<Option<SurfaceGeometry>> {
    let profile = |curve: &Curve| -> OgeomResult<Option<Curve>> {
        match curve {
            Curve::BSpline(b) if !b.is_periodic() => Ok(Some(Curve::BSpline(b.extended(
                at_end,
                length,
                b.degree(),
                tol,
            )?))),
            _ => Ok(None),
        }
    };
    Ok(match surface {
        SurfaceGeometry::BSpline(b) => {
            let degree = if along_u {
                b.u_knots().degree()
            } else {
                b.v_knots().degree()
            };
            Some(b.extended(along_u, at_end, length, degree, tol)?.into())
        }
        SurfaceGeometry::Extrusion(e) => {
            let (_, (e0, e1)) = surface.domain();
            if along_u {
                profile(e.curve())?
                    .map(|c| {
                        ogeom_geom::ExtrusionSurface::over(c, e.direction(), (e0, e1))
                            .map(|s| SurfaceGeometry::Extrusion(Box::new(s)))
                    })
                    .transpose()?
            } else {
                let extent = if at_end {
                    (e0, e1 + length)
                } else {
                    (e0 - length, e1)
                };
                Some(SurfaceGeometry::Extrusion(Box::new(
                    ogeom_geom::ExtrusionSurface::over(e.curve().clone(), e.direction(), extent)?,
                )))
            }
        }
        SurfaceGeometry::Revolution(r) if !along_u => {
            let ((a0, a1), _) = surface.domain();
            match profile(r.curve())? {
                Some(c) => {
                    let longer = ogeom_geom::RevolutionSurface::new(c, r.axis(), a1 - a0)?;
                    let same_turn = {
                        use ogeom_geom::Surface as _;
                        let ((b0, b1), _) = longer.domain();
                        (b0 - a0).abs() <= tol.parametric() && (b1 - a1).abs() <= tol.parametric()
                    };
                    same_turn.then(|| SurfaceGeometry::Revolution(Box::new(longer)))
                }
                None => None,
            }
        }
        _ => None,
    })
}
