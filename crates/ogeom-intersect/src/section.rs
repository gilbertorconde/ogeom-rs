//! Where two surfaces meet: the one call.
//!
//! Everything else in this crate is a stage — closed forms, seeding, tracing,
//! fitting. This is the function an application calls, and the one `ogeom-bool`
//! will build on: give it two surfaces, get back what they do to each other,
//! with the analytic path taken where it exists and the marched-and-fitted
//! path where it does not. The caller does not choose; the pair does.
//!
//! *Elsewhere* this is `GeomAPI_IntSS` over `IntPatch`/`GeomInt` — one entry
//! point hiding an analytic dispatch and a walking intersector.
//!
//! # What a section curve carries
//!
//! Three descriptions, because three consumers: the curve in space for the
//! edge, and a pcurve per surface for the faces — face splitting happens in
//! parameter space, and a curve a face cannot express is one it cannot be
//! split along. Analytic results carry exact pcurves where the projection has
//! a closed form and `None` where it does not; fitted results always carry
//! fitted pcurves, because the tracer recorded the parameters as it walked.
//!
//! A pcurve here is **same-parameter** with its 3D curve: evaluating either at
//! the same `t` lands on the same point of the intersection. That is the claim
//! `docs/DATA_MODEL.md` §6 makes edges carry, and it is arranged here by
//! construction — the 2D curves inherit the 3D curve's own parameterization —
//! rather than asserted and repaired later.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    Circle2d, Curve, Curve2d as _, Curve3d, Ellipse2d, Line2d, PlanarCurve, Surface,
    SurfaceGeometry,
};
use ogeom_math::{Circle2, Ellipse2, Frame2, Point, Point2};

use crate::approx::approximate_branch;
use crate::march::{Marching, branches, trace_tangential};
use crate::surface::{Meeting, surface_surface};

/// How to intersect, when the general path runs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntersectOptions {
    /// The tolerance the fitted curves are held to.
    pub tolerance: f64,
    /// The marching settings, for pairs with no closed form.
    pub marching: Marching,
}

impl Default for IntersectOptions {
    fn default() -> Self {
        Self {
            tolerance: 1e-6,
            marching: Marching::default(),
        }
    }
}

/// One curve of a section, with its parameter-space descriptions.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionCurve {
    /// The curve in space.
    pub curve: Curve,
    /// The curve in the first surface's parameter space, where it has one.
    ///
    /// Always present for a fitted curve. For an exact curve, present when the
    /// projection has a closed form — a line on a plane, a circle on the
    /// cylinder it wraps — and `None` where it does not, which is a statement
    /// about the projection rather than about the curve.
    pub on_a: Option<PlanarCurve>,
    /// The same, on the second surface.
    pub on_b: Option<PlanarCurve>,
    /// How far this curve may sit from the true intersection.
    ///
    /// Zero for an exact curve. For a fitted one, the trace's chord tolerance
    /// plus the fit's reported error — the sum of the stated parts.
    pub tolerance: f64,
    /// Whether the curve came from a closed form.
    pub exact: bool,
    /// Whether it is a closed loop.
    pub closed: bool,
    /// Whether the surfaces *touch* along this curve rather than crossing
    /// it.
    ///
    /// A tangential contact is a real curve — the two surfaces meet there,
    /// and a drawing has to show it — but it carries no boundary parity:
    /// neither surface passes through the other, so nothing is inside on
    /// one side and outside on the other. Consumers that classify by
    /// crossing must leave these out of that arithmetic; consumers that
    /// draw or measure contact want them.
    pub tangential: bool,
}

/// What two surfaces do to each other.
#[derive(Debug, Clone, PartialEq)]
pub enum SurfaceIntersection {
    /// They do not meet.
    ///
    /// From the general path this means *no crossing was found at the seeding
    /// resolution* — a branch thinner than the sampling grid is invisible to
    /// it, and [`coverage()`](crate::coverage()) is the instrument that checks.
    Apart,
    /// They touch at isolated points without crossing.
    Touching(Vec<Point>),
    /// They meet along these curves.
    Along(Vec<SectionCurve>),
    /// They are the same surface wherever they overlap.
    Same,
}

/// Where two surfaces meet.
///
/// The analytic path answers the pairs with closed forms — exactly, with
/// tolerance zero. Every other pair is seeded, traced and fitted to
/// `options.tolerance`. One call, and the pair decides the path.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the options
/// are unusable. A pair the marcher finds nothing for is [`Apart`], not an
/// error — see that variant for what it can and cannot claim.
///
/// [`Apart`]: SurfaceIntersection::Apart
pub fn intersect_surfaces(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: IntersectOptions,
    tol: Tolerances,
) -> OgeomResult<SurfaceIntersection> {
    if !options.tolerance.is_finite() || options.tolerance <= 0.0 {
        ogeom_bail!(
            Construction,
            "a tolerance of {} is not a distance",
            options.tolerance
        );
    }

    match surface_surface(a, b, tol) {
        Ok(Meeting::Apart) => Ok(SurfaceIntersection::Apart),
        Ok(Meeting::Same) => Ok(SurfaceIntersection::Same),
        Ok(Meeting::Touching(points)) => Ok(SurfaceIntersection::Touching(points)),
        Ok(Meeting::Along(curves)) => {
            let sections: Vec<SectionCurve> = curves
                .into_iter()
                .filter_map(|curve| exact_section(curve, a, b, tol))
                .collect();
            Ok(if sections.is_empty() {
                // Every curve fell outside the surfaces' stated extents: the
                // unbounded geometries meet, the surfaces as given do not.
                SurfaceIntersection::Apart
            } else {
                SurfaceIntersection::Along(sections)
            })
        }
        // No closed form for this pair: the statement that sends us marching.
        Err(_) => marched(a, b, options, tol),
    }
}

/// An exact curve dressed as a section, clipped to the surfaces it lies on.
///
/// The analytic layer works on the unbounded geometry — a plane and a cylinder
/// meet in unbounded lines — but the *surfaces* carry finite extents, and a
/// section running a billion units past both is not something an edge can be
/// built on. A line is clipped to the parameter interval where it is inside
/// both extents, through its exact pcurves; a curve wholly outside either
/// extent is dropped, or the boolean above would see a phantom edge on a
/// region the face does not have.
///
/// A *closed* curve partially outside an extent is kept whole: cutting it into
/// arcs is the restriction problem, and the restriction that matters is the
/// face's trim, which is §8's job — the extent here is only the surface's
/// parameterization window.
fn exact_section(
    curve: Curve,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<SectionCurve> {
    let closed = match &curve {
        Curve::Circle(_) | Curve::Ellipse(_) => true,
        _ => curve.is_closed(tol),
    };
    let range = curve.domain();
    let on_a = exact_pcurve(&curve, range, a, tol);
    let on_b = exact_pcurve(&curve, range, b, tol);

    if let Curve::Line(_) = &curve {
        // Clip through whichever pcurves exist; a missing pcurve leaves that
        // surface's extent unenforced, which errs long rather than wrong.
        let mut interval = curve.domain();
        if let Some(p) = &on_a {
            interval = intersect_intervals(interval, inside_box(p, a))?;
        }
        if let Some(p) = &on_b {
            interval = intersect_intervals(interval, inside_box(p, b))?;
        }
        let (lo, hi) = interval;
        let Curve::Line(line) = &curve else {
            unreachable!()
        };
        let clipped: Curve = ogeom_geom::LineCurve::over(line.axis(), lo, hi)
            .ok()?
            .into();
        let clip2 = |p: &PlanarCurve| -> Option<PlanarCurve> {
            let PlanarCurve::Line(l) = p else {
                return Some(p.clone());
            };
            Some(Line2d::over(l.axis(), lo, hi).ok()?.into())
        };
        let (ca, cb) = (on_a.as_ref().and_then(clip2), on_b.as_ref().and_then(clip2));
        let tangential = touching_along(&clipped, ca.as_ref(), cb.as_ref(), a, b, tol);
        return Some(SectionCurve {
            on_a: ca,
            on_b: cb,
            tolerance: 0.0,
            exact: true,
            closed: false,
            tangential,
            curve: clipped,
        });
    }

    // A closed curve: dropped only when wholly outside an extent it has a
    // pcurve to check against.
    for (pcurve, surface) in [(&on_a, a), (&on_b, b)] {
        if let Some(p) = pcurve
            && !touches_box(p, surface, tol)
        {
            return None;
        }
    }
    let tangential = touching_along(&curve, on_a.as_ref(), on_b.as_ref(), a, b, tol);
    Some(SectionCurve {
        on_a,
        on_b,
        tolerance: 0.0,
        exact: true,
        closed,
        tangential,
        curve,
    })
}

/// Whether the surfaces touch along an exact curve rather than crossing it:
/// their normals parallel at stations along its length.
///
/// Decided through the curve's own pcurves, which is where the normals can
/// be read without inverting anything. A curve missing a pcurve on either
/// surface is reported as a crossing — the honest default, since a section
/// nobody can place in a chart is one nothing can classify as contact
/// either.
fn touching_along(
    curve: &Curve,
    on_a: Option<&PlanarCurve>,
    on_b: Option<&PlanarCurve>,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> bool {
    let (Some(pa), Some(pb)) = (on_a, on_b) else {
        return false;
    };
    let (lo, hi) = curve.domain();
    for k in 0..5 {
        let t = (hi - lo).mul_add(f64::from(k) / 4.0, lo);
        let (Ok(ua), Ok(ub)) = (pa.point_at(t, tol), pb.point_at(t, tol)) else {
            return false;
        };
        let (Ok(na), Ok(nb)) = (a.normal_at(ua.x, ua.y, tol), b.normal_at(ub.x, ub.y, tol)) else {
            return false;
        };
        if na.vector().cross(nb.vector()).magnitude() > 1e-6 {
            return false;
        }
    }
    true
}

/// The parameter interval over which a 2D line stays inside a surface's
/// parameter box. `None` when it never enters.
fn inside_box(pcurve: &PlanarCurve, surface: &SurfaceGeometry) -> Option<(f64, f64)> {
    let PlanarCurve::Line(line) = pcurve else {
        return None;
    };
    let ((ua, ub), (va, vb)) = surface.domain();
    let axis = line.axis();
    let (o, d) = (axis.location, axis.direction.vector());

    // The slab test, one axis at a time.
    let mut lo = f64::NEG_INFINITY;
    let mut hi = f64::INFINITY;
    for (origin, direction, low, high) in [(o.x, d.x, ua, ub), (o.y, d.y, va, vb)] {
        if direction.abs() <= f64::MIN_POSITIVE {
            if origin < low || origin > high {
                return None;
            }
            continue;
        }
        let (a, b) = ((low - origin) / direction, (high - origin) / direction);
        let (near, far) = if a < b { (a, b) } else { (b, a) };
        lo = lo.max(near);
        hi = hi.min(far);
    }
    if lo >= hi {
        return None;
    }
    Some((lo, hi))
}

/// Whether any of a closed pcurve's samples lies inside the surface's box.
fn touches_box(pcurve: &PlanarCurve, surface: &SurfaceGeometry, tol: Tolerances) -> bool {
    use ogeom_geom::Curve2d;
    let ((ua, ub), (va, vb)) = surface.domain();
    let (lo, hi) = pcurve.domain();
    (0..=16).any(|i| {
        let t = lo + (hi - lo) * f64::from(i) / 16.0;
        pcurve.point_at(t, tol).is_ok_and(|p| {
            // Periodic directions always contain; only a bounded one excludes.
            let u_ok = surface.is_periodic_u() || (p.x >= ua && p.x <= ub);
            let v_ok = surface.is_periodic_v() || (p.y >= va && p.y <= vb);
            u_ok && v_ok
        })
    })
}

/// The overlap of two intervals. `None` when they miss.
fn intersect_intervals(a: (f64, f64), b: Option<(f64, f64)>) -> Option<(f64, f64)> {
    let b = b?;
    let (lo, hi) = (a.0.max(b.0), a.1.min(b.1));
    if lo >= hi {
        return None;
    }
    Some((lo, hi))
}

/// The general path: seed, trace, fit.
fn marched(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: IntersectOptions,
    tol: Tolerances,
) -> OgeomResult<SurfaceIntersection> {
    let traced = branches(a, b, options.marching, tol)?;
    if traced.is_empty() {
        return Ok(SurfaceIntersection::Apart);
    }
    let mut out = Vec::with_capacity(traced.len());
    let mut contacts: Vec<crate::march::Traced> = Vec::new();
    for branch in &traced {
        // A branch along which the two surfaces share their normal is a
        // tangency, not a crossing: the marcher's seeding cannot tell the
        // noise floor of a tangential valley from a genuine sign change, and
        // what it traces there is a stalled fragment of the valley, not a
        // section. The valley is still a curve, though, and the tangential
        // walker is the one that can follow it — so the fragment becomes a
        // seed rather than a discard, and what comes back is marked as
        // contact so nobody classifies by it.
        if branch_is_tangential(a, b, branch, tol)? {
            if let Some(contact) = walk_contact(a, b, branch, &contacts, options.marching, tol)? {
                contacts.push(contact);
            }
            continue;
        }
        let fitted = approximate_branch(a, b, branch, options.tolerance, tol)?;
        out.push(SectionCurve {
            curve: fitted.curve.into(),
            on_a: Some(fitted.on_a.into()),
            on_b: Some(fitted.on_b.into()),
            // The sum of the stated parts: the trace is within its chord of
            // the truth, the fit within its error of the trace.
            tolerance: options.marching.chord + fitted.fit_error,
            exact: false,
            closed: fitted.closed,
            tangential: false,
        });
    }
    for contact in &contacts {
        let fitted = approximate_branch(a, b, contact, options.tolerance, tol)?;
        out.push(SectionCurve {
            curve: fitted.curve.into(),
            on_a: Some(fitted.on_a.into()),
            on_b: Some(fitted.on_b.into()),
            tolerance: options.marching.chord + fitted.fit_error,
            exact: false,
            closed: fitted.closed,
            tangential: true,
        });
    }
    if out.is_empty() {
        return Ok(SurfaceIntersection::Apart);
    }
    Ok(SurfaceIntersection::Along(out))
}

/// Follow the contact a tangential fragment sits on, unless one already
/// traced covers it.
///
/// A tangential valley hands the crossing marcher several stalled fragments
/// — the seeds converge onto the contact from wherever they started and
/// wander there — so the fragments are candidates for *one* curve, not
/// several. A fragment whose middle already lies on a traced contact is one
/// of those repeats.
fn walk_contact(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    fragment: &crate::march::Traced,
    already: &[crate::march::Traced],
    marching: Marching,
    tol: Tolerances,
) -> OgeomResult<Option<crate::march::Traced>> {
    let middle = fragment.points.len() / 2;
    let Some(point) = fragment.points.get(middle).copied() else {
        return Ok(None);
    };
    for traced in already {
        // Traced points sit a step apart, so "on this curve" has to allow
        // half a step of gap to the nearest sample plus the chord budget.
        let spacing = traced
            .points
            .windows(2)
            .map(|w| w[0].distance(w[1]))
            .fold(0.0f64, f64::max);
        let near = traced
            .points
            .iter()
            .map(|p| p.distance(point))
            .fold(f64::INFINITY, f64::min);
        if near <= spacing.mul_add(0.5, marching.chord.max(tol.confusion())) {
            return Ok(None);
        }
    }
    let seed = crate::march::Contact {
        point,
        on_a: fragment.on_a[middle],
        on_b: fragment.on_b[middle],
    };
    // The walker refuses a seed that is not a contact; that refusal is an
    // answer, not a failure — the fragment simply had nothing to follow.
    // A walk that stalls where it started says the same thing in points:
    // too few to fit, so there is no contact curve to report here.
    Ok(trace_tangential(a, b, seed, marching, tol)
        .ok()
        .filter(|traced| traced.points.len() >= 4))
}

/// Whether a traced branch runs along a tangency of the two surfaces:
/// their normals parallel, sampled along its length.
fn branch_is_tangential(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    branch: &crate::march::Traced,
    tol: Tolerances,
) -> OgeomResult<bool> {
    use ogeom_geom::Surface as _;
    let count = branch.points.len();
    if count == 0 {
        return Ok(true);
    }
    for k in 0..5 {
        let i = (k * (count - 1)) / 4;
        let (ua, va) = branch.on_a[i.min(count - 1)];
        let (ub, vb) = branch.on_b[i.min(count - 1)];
        let (dau, dav) = a.d1_at(ua, va, tol)?;
        let (dbu, dbv) = b.d1_at(ub, vb, tol)?;
        let na = dau.cross(dav);
        let nb = dbu.cross(dbv);
        let (ma, mb) = (na.magnitude(), nb.magnitude());
        if ma <= tol.confusion() || mb <= tol.confusion() {
            continue;
        }
        if na.cross(nb).magnitude() / (ma * mb) > 1e-2 {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The exact pcurve of a curve lying on a surface, where the projection has
/// a closed form — `None` where it does not.
///
/// Public because the boolean's same-domain handling needs it: two faces on
/// one geometric surface may still carry different charts, and the other
/// face's boundary edges have to be spoken in this face's parameters before
/// they can split it.
#[must_use]
pub fn exact_pcurve_of(
    curve: &Curve,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    exact_pcurve(curve, curve.domain(), surface, tol)
}

/// As [`exact_pcurve_of`], with the parameter range the caller actually
/// uses.
///
/// A curve's chart image can depend on *which part* of the curve is meant: a
/// ruling on a cone crosses the apex, and its angle on the far nappe is half
/// a turn from its angle on the near one. The curve's own domain may span
/// both — an imported line's usually does — so a caller that knows its edge's
/// range must say so, or the exact projection may answer for the wrong side.
#[must_use]
pub fn exact_pcurve_over(
    curve: &Curve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    exact_pcurve(curve, range, surface, tol)
}

/// The exact pcurve of an analytic curve on an analytic surface, where the
/// projection has a closed form.
///
/// Same-parameter by construction: each 2D curve inherits the 3D curve's own
/// parameterization, so the two evaluate to the same point of the intersection
/// at the same `t`. The cases are the ones where that inheritance is exact —
/// anything else returns `None` rather than a fit, because an *exact* result
/// with a fitted pcurve would be a curve whose descriptions disagree by an
/// amount nothing on it records.
fn exact_pcurve(
    curve: &Curve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    match surface {
        SurfaceGeometry::Plane(p) => on_plane(curve, p.plane(), tol),
        SurfaceGeometry::Cylinder(c) => on_cylinder(curve, range, c.cylinder(), tol),
        SurfaceGeometry::Sphere(s) => on_sphere(curve, range, s.sphere(), tol),
        SurfaceGeometry::Torus(t) => on_torus(curve, t.torus(), tol),
        SurfaceGeometry::Cone(c) => on_cone(curve, range, c.cone(), tol),
        _ => None,
    }
}

/// The pcurve of a curve on a cone, for the two straight-line families.
///
/// A ruling — through the apex, on the surface — runs at constant `u`; a
/// circle perpendicular to the axis, centred on it, with the radius the cone
/// has at that height, runs at constant `v`. Both inherit the 3D curve's own
/// parameter, the circle with phase and winding exactly as the cylinder case.
/// The ruling's angle is measured over `range`, because the same line has
/// the opposite angle on the other side of the apex.
fn on_cone(
    curve: &Curve,
    range: (f64, f64),
    cone: ogeom_math::Cone,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    let frame = cone.frame();
    let axis_z = frame.z().vector();
    let tau = core::f64::consts::TAU;
    match curve {
        Curve::Circle(c) => {
            let circle = c.circle();
            if circle.frame().z().vector().cross(axis_z).magnitude() > tol.angular() {
                return None;
            }
            let local = frame.to_local(circle.centre());
            if local.x.hypot(local.y) > tol.confusion() {
                return None;
            }
            // The cone's radius at the circle's height must be the circle's.
            let expected = cone
                .half_angle()
                .tan()
                .mul_add(local.z, cone.reference_radius());
            if (expected - circle.radius()).abs() > tol.confusion() * 10.0 {
                return None;
            }
            let start = circle.centre() + circle.frame().x().vector() * circle.radius();
            let at = frame.to_local(start);
            let phase = at.y.atan2(at.x);
            let winding = circle.frame().z().vector().dot(axis_z).signum();
            let towards =
                ogeom_math::Direction2::new(ogeom_math::Vector2::new(winding, 0.0), tol).ok()?;
            Some(
                Line2d::over(
                    ogeom_math::Axis2::new(Point2::new(phase, local.z), towards),
                    0.0,
                    tau,
                )
                .ok()?
                .into(),
            )
        }
        Curve::Line(line) => {
            // A ruling: verified by sample, not assumed — three points on
            // the surface pin a line to it.
            let axis = line.axis();
            let on = |t: f64| {
                let p = axis.location + axis.direction.vector() * t;
                cone.distance_to(p) <= tol.confusion() * 10.0
            };
            if !on(0.0) || !on(1.0) || !on(-1.0) {
                return None;
            }
            // A ruling reaching the tip may be *stated* from the apex
            // itself — where the angle is atan2(0, 0), garbage — and its
            // own domain usually spans both nappes, where the angles differ
            // by half a turn. Measure the angle at whichever end of the
            // *used* range stands farthest from the axis: that is the side
            // the caller means.
            let (lo, hi) = if range.0.is_finite() && range.1.is_finite() && range.0 != range.1 {
                range
            } else {
                line.domain()
            };
            let mut local = frame.to_local(axis.location);
            for t in [lo, hi] {
                if !t.is_finite() {
                    continue;
                }
                let candidate = frame.to_local(axis.location + axis.direction.vector() * t);
                if candidate.x.hypot(candidate.y) > local.x.hypot(local.y) {
                    local = candidate;
                }
            }
            if local.x.hypot(local.y) <= tol.confusion() {
                return None;
            }
            let u = local.y.atan2(local.x).rem_euclid(tau);
            // Same-parameter exactly: a degree-one spline over the used
            // range maps t linearly onto the chart column, whatever rate
            // the slant climbs at.
            let v_at = |t: f64| {
                frame
                    .to_local(axis.location + axis.direction.vector() * t)
                    .z
            };
            let knots = ogeom_math::KnotVector::new(vec![lo, lo, hi, hi], 1).ok()?;
            Some(
                ogeom_geom::BSpline2d::new(
                    knots,
                    vec![Point2::new(u, v_at(lo)), Point2::new(u, v_at(hi))],
                    tol,
                )
                .ok()?
                .into(),
            )
        }
        _ => None,
    }
}

/// The pcurve of a circle on a torus, for the two families that are straight
/// lines in `(u, v)`.
///
/// A *parallel* — centred on the axis, in a plane perpendicular to it — runs
/// at constant `v`; a *tube circle* — minor radius, centred on the tube's
/// spine, in a plane through the axis — runs at constant `u`. Both inherit
/// the circle's own angle, phase and winding included, exactly as the
/// cylinder case does; the STEP reader is the consumer that forced the torus
/// into this list, fillet faces being tori more often than not.
fn on_torus(curve: &Curve, torus: ogeom_math::Torus, tol: Tolerances) -> Option<PlanarCurve> {
    let Curve::Circle(c) = curve else {
        return None;
    };
    let circle = c.circle();
    let frame = torus.frame();
    let axis_z = frame.z().vector();
    let normal = circle.frame().z().vector();
    let local = frame.to_local(circle.centre());
    let tau = core::f64::consts::TAU;

    // A parallel of the sweep.
    if normal.cross(axis_z).magnitude() <= tol.angular()
        && local.x.hypot(local.y) <= tol.confusion()
    {
        let sin_v = local.z / torus.minor_radius();
        let cos_v = (circle.radius() - torus.major_radius()) / torus.minor_radius();
        if (sin_v.hypot(cos_v) - 1.0).abs() > tol.confusion() {
            return None;
        }
        let v = sin_v.atan2(cos_v);
        let start = circle.centre() + circle.frame().x().vector() * circle.radius();
        let at = frame.to_local(start);
        let phase = at.y.atan2(at.x);
        let winding = normal.dot(axis_z).signum();
        let towards =
            ogeom_math::Direction2::new(ogeom_math::Vector2::new(winding, 0.0), tol).ok()?;
        return Some(
            Line2d::over(
                ogeom_math::Axis2::new(Point2::new(phase, v), towards),
                0.0,
                tau,
            )
            .ok()?
            .into(),
        );
    }

    // A circle of the tube.
    if (circle.radius() - torus.minor_radius()).abs() <= tol.confusion()
        && normal.dot(axis_z).abs() <= tol.angular()
        && (local.x.hypot(local.y) - torus.major_radius()).abs() <= tol.confusion()
        && local.z.abs() <= tol.confusion()
    {
        let u = local.y.atan2(local.x);
        let radial = frame.x().vector() * u.cos() + frame.y().vector() * u.sin();
        let xc = circle.frame().x().vector();
        let phase = xc.dot(axis_z).atan2(xc.dot(radial));
        let winding = normal.dot(radial.cross(axis_z)).signum();
        let towards =
            ogeom_math::Direction2::new(ogeom_math::Vector2::new(0.0, winding), tol).ok()?;
        return Some(
            Line2d::over(
                ogeom_math::Axis2::new(Point2::new(u, phase), towards),
                0.0,
                tau,
            )
            .ok()?
            .into(),
        );
    }
    None
}

/// Project a curve lying in a plane into the plane's own coordinates.
///
/// Exact for a line, a circle and an ellipse: the plane's frame is orthonormal,
/// so lengths and the curves' own parameterizations survive the projection
/// unchanged.
fn on_plane(curve: &Curve, plane: ogeom_math::Plane, tol: Tolerances) -> Option<PlanarCurve> {
    let frame = plane.frame();
    let flat = |p: Point| {
        let local = frame.to_local(p);
        Point2::new(local.x, local.y)
    };
    let flat_direction = |d: ogeom_math::Direction| {
        let tip = flat(frame.origin() + d.vector());
        ogeom_math::Direction2::new(tip - flat(frame.origin()), tol).ok()
    };
    match curve {
        Curve::Line(line) => {
            let axis = line.axis();
            let through = flat(axis.location);
            let direction = flat_direction(axis.direction)?;
            let (lo, hi) = line.domain();
            Some(
                Line2d::over(ogeom_math::Axis2::new(through, direction), lo, hi)
                    .ok()?
                    .into(),
            )
        }
        Curve::Circle(c) => {
            let circle = c.circle();
            let frame2 = Frame2::from_axes(
                flat(circle.centre()),
                flat_direction(circle.frame().x())?,
                flat_direction(circle.frame().y())?,
                tol,
            )
            .ok()?;
            Some(Circle2d::new(Circle2::new(frame2, circle.radius(), tol).ok()?).into())
        }
        Curve::Ellipse(e) => {
            let ellipse = e.ellipse();
            let frame2 = Frame2::from_axes(
                flat(ellipse.centre()),
                flat_direction(ellipse.frame().x())?,
                flat_direction(ellipse.frame().y())?,
                tol,
            )
            .ok()?;
            Some(
                Ellipse2d::new(
                    Ellipse2::new(frame2, ellipse.major_radius(), ellipse.minor_radius(), tol)
                        .ok()?,
                )
                .into(),
            )
        }
        Curve::BSpline(b) => {
            // Affine invariance: a (rational) B-spline in the plane projects
            // into the plane's own coordinates control point by control
            // point, knots and weights untouched — exact, and same-parameter
            // by construction.
            let control = b
                .control_points()
                .iter()
                .map(|w| ogeom_math::Weighted::new(flat((*w).point()), w.weight, tol))
                .collect::<Result<Vec<_>, _>>()
                .ok()?;
            Some(
                ogeom_geom::BSpline2d::rational(b.knots().clone(), control)
                    .ok()?
                    .into(),
            )
        }
        _ => None,
    }
}

/// The pcurve of a curve on a cylinder, where it is a straight line in
/// parameter space.
///
/// A line along the axis runs at constant `u`; a full circle around it runs at
/// constant `v`. Both are lines in `(u, v)`, exactly, and both inherit the 3D
/// curve's own parameter — height for the line, angle for the circle.
fn on_cylinder(
    curve: &Curve,
    range: (f64, f64),
    cylinder: ogeom_math::Cylinder,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    let axis = cylinder.axis();
    let frame = cylinder.frame();
    match curve {
        Curve::Line(line) => {
            // Parallel to the axis, on the surface.
            let direction = line.axis().direction;
            let along = direction.dot(axis.direction);
            if (along.abs() - 1.0).abs() > tol.angular() {
                return None;
            }
            let through = line.axis().location;
            if (axis.distance_to(through) - cylinder.radius()).abs() > tol.confusion() {
                return None;
            }
            let local = frame.to_local(through);
            let u = local.y.atan2(local.x).rem_euclid(core::f64::consts::TAU);
            // The 3D line's parameter is length from its origin; at constant u
            // the pcurve's `v` runs at the same rate, signed by whether the
            // line runs with the axis or against it.
            let (lo, hi) = line.domain();
            let start = Point2::new(u, local.z);
            let towards =
                ogeom_math::Direction2::new(ogeom_math::Vector2::new(0.0, along.signum()), tol)
                    .ok()?;
            Some(
                Line2d::over(ogeom_math::Axis2::new(start, towards), lo, hi)
                    .ok()?
                    .into(),
            )
        }
        Curve::Circle(c) => {
            let circle = c.circle();
            // Perpendicular to the axis, centred on it, of the same radius.
            if circle
                .frame()
                .z()
                .cross_with(axis.direction.vector())
                .magnitude()
                > tol.angular()
            {
                return None;
            }
            if axis.distance_to(circle.centre()) > tol.confusion() {
                return None;
            }
            if (circle.radius() - cylinder.radius()).abs() > tol.confusion() {
                return None;
            }
            let local = frame.to_local(circle.centre());
            // Where the circle's own angle zero sits in the cylinder's angle —
            // and which way its parameter runs around the axis. A section
            // circle inherits its winding from the pair that made it, and one
            // wound against the cylinder's `u` — a circle cut by a plane whose
            // normal opposes the axis — runs its pcurve in `-u`. Writing `+u`
            // unconditionally here was the bug the boolean's drill test found:
            // the pcurve evaluated half a turn away from the curve, and the
            // face's arrangement tore along a seam that was not there.
            let start = circle.centre() + circle.frame().x().vector() * circle.radius();
            let at = frame.to_local(start);
            let phase = at.y.atan2(at.x);
            let winding = circle.frame().z().dot(axis.direction).signum();
            let towards =
                ogeom_math::Direction2::new(ogeom_math::Vector2::new(winding, 0.0), tol).ok()?;
            Some(
                Line2d::over(
                    ogeom_math::Axis2::new(Point2::new(phase, local.z), towards),
                    0.0,
                    core::f64::consts::TAU,
                )
                .ok()?
                .into(),
            )
        }
        Curve::Ellipse(_) => {
            // An oblique plane's section: its plan projection is the
            // cylinder's own cross-section circle traced *uniformly*, so
            // the chart trace is u = s·t + φ, v = c₀ + a·cos t + b·sin t —
            // the trig-affine family. Derived from the curve's own
            // evaluations and verified by sample, never assumed.
            use ogeom_geom::Curve3d as _;
            let tau = core::f64::consts::TAU;
            let local = |t: f64| -> Option<ogeom_math::Point> {
                Some(frame.to_local(curve.point_at(t, tol).ok()?))
            };
            let l0 = local(0.0)?;
            let lq = local(tau / 4.0)?;
            let lh = local(tau / 2.0)?;
            // On the surface at all: plan radius must be the cylinder's.
            let r = cylinder.radius();
            for l in [&l0, &lq, &lh] {
                if (l.x.hypot(l.y) - r).abs() > tol.confusion() * 10.0 {
                    return None;
                }
            }
            let phase = l0.y.atan2(l0.x);
            // Winding from the quarter-turn sample: uniform tracing puts it
            // a quarter turn away, one side or the other.
            let uq = lq.y.atan2(lq.x);
            let step = (uq - phase).rem_euclid(tau);
            let winding = if (step - tau / 4.0).abs() < 1e-6 {
                1.0
            } else if (step - 3.0 * tau / 4.0).abs() < 1e-6 {
                -1.0
            } else {
                return None;
            };
            // Height coefficients from three samples.
            let c0 = f64::midpoint(l0.z, lh.z);
            let a = (l0.z - lh.z) / 2.0;
            let b = lq.z - c0;
            // The trig formula is global — cosine wraps, the linear angle
            // unwraps the chart — so the pcurve lives on whatever range the
            // edge actually spans, a loop crossing the period included.
            let candidate = ogeom_geom::Trig2d::new(
                Point2::new(phase, c0),
                ogeom_math::Vector2::new(winding, 0.0),
                ogeom_math::Vector2::new(0.0, a),
                ogeom_math::Vector2::new(0.0, b),
                range,
            )
            .ok()?;
            // The same-parameter law, verified at points the derivation
            // never touched, inside the range the edge will use.
            use ogeom_geom::Curve2d as _;
            for i in 0..7 {
                let t = range.0 + (range.1 - range.0) * (0.09 + 0.13 * f64::from(i)) / 0.91;
                let l = local(t)?;
                let chart = candidate.point_at(t, tol).ok()?;
                let du = (chart.x - l.y.atan2(l.x)).rem_euclid(tau);
                if du.min(tau - du) > 1e-9 {
                    return None;
                }
                if (chart.y - l.z).abs() > tol.confusion() * 10.0 {
                    return None;
                }
            }
            Some(PlanarCurve::Trig(candidate))
        }
        _ => None,
    }
}

/// The pcurve of half a meridian: a great circle through both poles,
/// restricted to one side of them.
///
/// The whole circle has no chart image a single curve can carry — its
/// longitude jumps by half a turn at each pole — but each *half* does, and it
/// is a straight line. Writing the circle's own parameter as `t` and the
/// sphere's axis as `Z = cos α·X + sin α·Y` in the circle's own frame, the
/// point's height above the equator is `r·cos(t − α)`, so the latitude is
/// `asin(cos(t − α))`, which on `t − α ∈ [0, π]` is exactly `π/2 − (t − α)` —
/// affine in `t`, with slope one. The longitude is constant on that half and
/// half a turn away on the other. So the pcurve is a vertical line in the
/// chart, sharing the circle's parameter exactly, and the caller's `range` is
/// what says which half is meant.
///
/// The half is not assumed: the returned line is lifted back through the
/// sphere at stations along the range and compared against the circle, so a
/// misread orientation is caught here rather than downstream.
fn on_meridian(
    curve: &ogeom_geom::CircleCurve,
    range: (f64, f64),
    sphere: ogeom_math::Sphere,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    let circle = curve.circle();
    // A reversed circle runs its own angle backwards, and the shifted angle
    // below is measured in the *curve's* parameter, so the sign travels with
    // it: the sweep flips and so do both the latitude's slope and which half
    // of the circle a range names.
    let sweep = if curve.is_reversed() { -1.0 } else { 1.0 };
    let frame = sphere.frame();
    let z = frame.z().vector();
    // A great circle: the sphere's own centre and radius, in a plane holding
    // the axis. Anything else is not a meridian.
    if circle.centre().distance(sphere.centre()) > tol.confusion() {
        return None;
    }
    if (circle.radius() - sphere.radius()).abs() > tol.confusion() {
        return None;
    }
    let (cx, cy) = (circle.frame().x().vector(), circle.frame().y().vector());
    let (xz, yz) = (cx.dot(z), cy.dot(z));
    // The axis must lie *in* the circle's plane, or the circle is neither a
    // parallel nor a meridian and has no closed-form chart image at all.
    if xz.hypot(yz) < 1.0 - tol.angular() {
        return None;
    }
    let raw_alpha = yz.atan2(xz);
    // `w` is the circle's own horizontal direction: the axis turned a quarter
    // turn within the circle's plane.
    let w = cx * -raw_alpha.sin() + cy * raw_alpha.cos();
    let local = frame.to_local(sphere.centre() + w);
    let longitude = local.y.atan2(local.x);

    let half = core::f64::consts::PI;
    let mid = f64::midpoint(range.0, range.1);
    // Where the range sits relative to the poles, in the shifted angle
    // `x = sweep·t − α` that measures the descent from the north pole.
    let x_mid = (sweep * mid - raw_alpha).rem_euclid(core::f64::consts::TAU);
    let x_mid = if x_mid > half {
        x_mid - core::f64::consts::TAU
    } else {
        x_mid
    };
    let span = sweep * (range.1 - range.0);
    let (mut x0, mut x1) = (x_mid - span / 2.0, x_mid + span / 2.0);
    if x0 > x1 {
        core::mem::swap(&mut x0, &mut x1);
    }
    // The turn count `α` was written with is what decides whether the
    // latitude comes out inside the chart or a whole turn away from it, so
    // the branch the range actually sits on is the one the line is built
    // from.
    let alpha = sweep.mul_add(mid, -x_mid);
    let slack = tol.parametric().max(1e-9);
    let (axis_point, towards) = if x0 >= -slack && x1 <= half + slack {
        // The descending half: latitude π/2 − (sweep·t − α), longitude
        // constant.
        (
            Point2::new(longitude, half.mul_add(0.5, alpha)),
            ogeom_math::Vector2::new(0.0, -sweep),
        )
    } else if x0 >= -half - slack && x1 <= slack {
        // The ascending half, half a turn round the chart.
        (
            Point2::new(longitude + half, half.mul_add(0.5, -alpha)),
            ogeom_math::Vector2::new(0.0, sweep),
        )
    } else {
        // The range straddles a pole: no one line covers it.
        return None;
    };
    let towards = ogeom_math::Direction2::new(towards, tol).ok()?;
    let margin = (range.1 - range.0) * 0.25;
    let line: PlanarCurve = Line2d::over(
        ogeom_math::Axis2::new(axis_point, towards),
        range.0 - margin,
        range.1 + margin,
    )
    .ok()?
    .into();

    // Measured, not assumed: the chart line lifted back through the sphere is
    // the circle it claims to be.
    for k in 0..=4 {
        let t = (range.1 - range.0).mul_add(f64::from(k) / 4.0, range.0);
        let uv = line.point_at(t, tol).ok()?;
        let lifted = ogeom_math::elementary::sphere_at(&sphere, uv.x, uv.y).point;
        let want = curve.point_at(t, tol).ok()?;
        if lifted.distance(want) > tol.confusion() {
            return None;
        }
    }
    Some(line)
}

/// The pcurve of a circle on a sphere: a parallel of latitude, or one half of
/// a meridian.
fn on_sphere(
    curve: &Curve,
    range: (f64, f64),
    sphere: ogeom_math::Sphere,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    let Curve::Circle(c) = curve else {
        return None;
    };
    let circle = c.circle();
    let frame = sphere.frame();
    // Perpendicular to the sphere's axis and centred on it: a parallel of
    // latitude, which is a horizontal line in (longitude, latitude).
    if circle
        .frame()
        .z()
        .cross_with(frame.z().vector())
        .magnitude()
        > tol.angular()
    {
        return on_meridian(c, range, sphere, tol);
    }
    let local = frame.to_local(circle.centre());
    if local.x.abs() > tol.confusion() || local.y.abs() > tol.confusion() {
        return None;
    }
    let latitude = (local.z / sphere.radius()).clamp(-1.0, 1.0).asin();
    // Sanity: the circle's radius must be the parallel's.
    if (circle.radius() - sphere.radius() * latitude.cos()).abs() > tol.confusion() {
        return None;
    }
    let start = circle.centre() + circle.frame().x().vector() * circle.radius();
    let at = frame.to_local(start);
    let phase = at.y.atan2(at.x);
    // Phase and winding exactly as the cylinder case: a parallel whose own
    // axis opposes the sphere's marches its angle *down* the longitude.
    let winding = circle.frame().z().vector().dot(frame.z().vector()).signum();
    let towards = ogeom_math::Direction2::new(ogeom_math::Vector2::new(winding, 0.0), tol).ok()?;
    Some(
        Line2d::over(
            ogeom_math::Axis2::new(Point2::new(phase, latitude), towards),
            0.0,
            core::f64::consts::TAU,
        )
        .ok()?
        .into(),
    )
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use ogeom_geom::{Curve2d, Curve3d, CylinderSurface, PlaneSurface, SphereSurface};
    use ogeom_math::{Cylinder, Direction, Frame, Plane, Sphere, Vector};

    const T: Tolerances = Tolerances::millimetres();

    fn sphere(centre: Point, radius: f64) -> SurfaceGeometry {
        SphereSurface::new(Sphere::centred(centre, radius, T).unwrap()).into()
    }

    fn cylinder(axis: Vector, radius: f64) -> SurfaceGeometry {
        let frame = Frame::new(
            Point::ORIGIN,
            Direction::new(axis, T).unwrap(),
            Direction::from_cross(axis, Vector::new(0.3, 0.5, 0.9), T).unwrap(),
            T,
        )
        .unwrap();
        CylinderSurface::new(Cylinder::new(frame, radius, T).unwrap(), (-4.0, 4.0))
            .unwrap()
            .into()
    }

    fn plane(origin: Point, normal: Vector) -> SurfaceGeometry {
        PlaneSurface::over(
            Plane::through(origin, Direction::new(normal, T).unwrap()),
            (-6.0, 6.0),
            (-6.0, 6.0),
        )
        .unwrap()
        .into()
    }

    /// Same-parameter: pcurve lifted through its surface equals the 3D curve,
    /// at the same parameter, everywhere sampled.
    fn assert_same_parameter(
        section: &SectionCurve,
        surface: &SurfaceGeometry,
        pcurve: &PlanarCurve,
        samples: usize,
    ) {
        let (lo, hi) = section.curve.domain();
        let (plo, phi) = pcurve.domain();
        assert!(
            (lo - plo).abs() < 1e-9 && (hi - phi).abs() < 1e-9,
            "domains disagree: [{lo}, {hi}] against [{plo}, {phi}]"
        );
        for i in 0..=samples {
            #[allow(clippy::cast_precision_loss)]
            let t = lo + (hi - lo) * i as f64 / samples as f64;
            let on_curve = section.curve.point_at(t, T).unwrap();
            let at = pcurve.point_at(t, T).unwrap();
            let lifted = surface.point_at(at.x, at.y, T).unwrap();
            assert!(
                on_curve.is_equal(lifted, T),
                "at t = {t}: curve {on_curve:?}, lifted {lifted:?}"
            );
        }
    }

    #[test]
    fn an_analytic_pair_comes_back_exact_with_matching_pcurves() {
        // A plane through a cylinder's axis: two lines, and every description
        // agrees at the same parameter — which is the claim edges carry and
        // booleans rely on.
        let drum = cylinder(Vector::Z, 2.0);
        let cut = plane(Point::ORIGIN, Vector::X);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&drum, &cut, IntersectOptions::default(), T).unwrap()
        else {
            panic!("a plane through a cylinder meets it along curves");
        };
        assert_eq!(curves.len(), 2);
        for section in &curves {
            assert!(section.exact);
            assert!((section.tolerance - 0.0).abs() < f64::EPSILON);
            let on_a = section.on_a.as_ref().expect("a line has a cylinder pcurve");
            let on_b = section.on_b.as_ref().expect("and a plane pcurve");
            assert_same_parameter(section, &drum, on_a, 50);
            assert_same_parameter(section, &cut, on_b, 50);
        }
    }

    #[test]
    fn an_oblique_cut_gives_the_ellipse_a_trig_pcurve_on_the_drum() {
        // The pcurve the deferred table owed: the oblique ellipse runs
        // linearly in the chart angle and sinusoidally in height — the
        // trig-affine family — exactly, same-parameter, both sides.
        let drum = cylinder(Vector::Z, 2.0);
        let angle: f64 = 0.5;
        let cut = plane(Point::ORIGIN, Vector::new(0.0, angle.sin(), angle.cos()));
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&drum, &cut, IntersectOptions::default(), T).unwrap()
        else {
            panic!("an oblique plane meets the cylinder along its ellipse");
        };
        assert_eq!(curves.len(), 1);
        let section = &curves[0];
        assert!(section.exact);
        assert!(matches!(section.curve, Curve::Ellipse(_)));
        let on_drum = section
            .on_a
            .as_ref()
            .expect("the oblique ellipse now carries its cylinder pcurve");
        assert!(
            matches!(on_drum, PlanarCurve::Trig(_)),
            "the chart trace is trig-affine: {on_drum:?}"
        );
        assert_same_parameter(section, &drum, on_drum, 60);
        let on_plane = section.on_b.as_ref().expect("and its plane pcurve");
        assert_same_parameter(section, &cut, on_plane, 60);
    }

    #[test]
    fn a_perpendicular_cut_gives_a_circle_with_a_straight_pcurve() {
        let drum = cylinder(Vector::Z, 2.0);
        let cut = plane(Point::new(0.0, 0.0, 1.0), Vector::Z);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&drum, &cut, IntersectOptions::default(), T).unwrap()
        else {
            panic!("expected curves");
        };
        assert_eq!(curves.len(), 1);
        let section = &curves[0];
        assert!(section.closed);
        assert!(matches!(section.curve, Curve::Circle(_)));
        // On the cylinder the circle is a horizontal line in (u, v).
        assert!(matches!(
            section.on_a.as_ref().unwrap(),
            PlanarCurve::Line(_)
        ));
        assert_same_parameter(section, &drum, section.on_a.as_ref().unwrap(), 60);
        assert_same_parameter(section, &cut, section.on_b.as_ref().unwrap(), 60);
    }

    #[test]
    fn coaxial_cylinder_and_sphere_give_circles_with_pcurves_on_both() {
        let drum = cylinder(Vector::Z, 1.5);
        let ball = sphere(Point::ORIGIN, 3.0);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&drum, &ball, IntersectOptions::default(), T).unwrap()
        else {
            panic!("expected curves");
        };
        assert_eq!(curves.len(), 2);
        for section in &curves {
            assert!(section.exact);
            assert_same_parameter(section, &drum, section.on_a.as_ref().unwrap(), 40);
            assert_same_parameter(section, &ball, section.on_b.as_ref().unwrap(), 40);
        }
    }

    fn torus(origin: Point, axis: Vector, major: f64, minor: f64) -> SurfaceGeometry {
        let frame = Frame::new(
            origin,
            Direction::new(axis, T).unwrap(),
            Direction::from_cross(axis, Vector::new(0.3, 0.5, 0.9), T).unwrap(),
            T,
        )
        .unwrap();
        ogeom_geom::TorusSurface::new(ogeom_math::Torus::new(frame, major, minor, T).unwrap())
            .into()
    }

    #[test]
    fn an_axis_normal_plane_meets_a_torus_in_two_parallels_with_pcurves() {
        let ring = torus(Point::ORIGIN, Vector::Z, 2.0, 0.5);
        let cut = plane(Point::new(0.0, 0.0, 0.3), Vector::Z);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&ring, &cut, IntersectOptions::default(), T).unwrap()
        else {
            panic!("an axis-normal plane through the tube meets it along curves");
        };
        assert_eq!(curves.len(), 2);
        let spread = 0.5_f64.mul_add(0.5, -(0.3 * 0.3)).sqrt();
        let mut radii: Vec<f64> = curves
            .iter()
            .map(|s| {
                let Curve::Circle(c) = &s.curve else {
                    panic!("a parallel is a circle");
                };
                c.circle().radius()
            })
            .collect();
        radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((radii[0] - (2.0 - spread)).abs() < 1e-12);
        assert!((radii[1] - (2.0 + spread)).abs() < 1e-12);
        for section in &curves {
            assert!(section.exact);
            assert_same_parameter(section, &ring, section.on_a.as_ref().unwrap(), 48);
            assert_same_parameter(section, &cut, section.on_b.as_ref().unwrap(), 48);
        }
    }

    #[test]
    fn the_plane_a_ball_rolls_on_touches_its_torus_along_the_circle_it_rolled() {
        // Tangency with length is reported as the curve it is — the way a
        // tangent plane reports its line on a cylinder — because the blend
        // machinery builds faces whose boundaries are exactly these circles,
        // and a Touching with no curve in it would read as a refusal upstream.
        let ring = torus(Point::ORIGIN, Vector::Z, 2.0, 0.5);
        let cut = plane(Point::new(0.0, 0.0, 0.5), Vector::Z);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&ring, &cut, IntersectOptions::default(), T).unwrap()
        else {
            panic!("the rolling plane touches along a circle, not at points");
        };
        assert_eq!(curves.len(), 1);
        let Curve::Circle(c) = &curves[0].curve else {
            panic!("the tangency is a circle");
        };
        assert!((c.circle().radius() - 2.0).abs() < 1e-12);
        assert_same_parameter(&curves[0], &ring, curves[0].on_a.as_ref().unwrap(), 48);
        assert_same_parameter(&curves[0], &cut, curves[0].on_b.as_ref().unwrap(), 48);
    }

    #[test]
    fn a_coaxial_cylinder_meets_a_torus_in_two_parallels_and_touches_in_one() {
        let ring = torus(Point::ORIGIN, Vector::Z, 2.0, 0.5);
        let drum = cylinder(Vector::Z, 2.2);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&drum, &ring, IntersectOptions::default(), T).unwrap()
        else {
            panic!("a coaxial cylinder through the tube meets it along curves");
        };
        assert_eq!(curves.len(), 2);
        for section in &curves {
            assert!(section.exact);
            let Curve::Circle(c) = &section.curve else {
                panic!("a parallel is a circle");
            };
            assert!((c.circle().radius() - 2.2).abs() < 1e-12);
            assert_same_parameter(section, &drum, section.on_a.as_ref().unwrap(), 48);
            assert_same_parameter(section, &ring, section.on_b.as_ref().unwrap(), 48);
        }

        // Tangent at the tube's outer equator: one circle, with both pcurves.
        let grazing = cylinder(Vector::Z, 2.5);
        let SurfaceIntersection::Along(touch) =
            intersect_surfaces(&grazing, &ring, IntersectOptions::default(), T).unwrap()
        else {
            panic!("the grazing cylinder touches along the equator");
        };
        assert_eq!(touch.len(), 1);
        assert_same_parameter(&touch[0], &grazing, touch[0].on_a.as_ref().unwrap(), 48);
        assert_same_parameter(&touch[0], &ring, touch[0].on_b.as_ref().unwrap(), 48);
    }

    #[test]
    fn coaxial_tori_are_the_same_or_meet_in_parallels() {
        let ring = torus(Point::ORIGIN, Vector::Z, 2.0, 0.5);
        assert!(matches!(
            intersect_surfaces(&ring, &ring.clone(), IntersectOptions::default(), T).unwrap(),
            SurfaceIntersection::Same
        ));

        // The same tube lifted half a radius: the profile circles cross
        // twice, and each crossing revolves into a parallel shared exactly.
        let lifted = torus(Point::new(0.0, 0.0, 0.5), Vector::Z, 2.0, 0.5);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&ring, &lifted, IntersectOptions::default(), T).unwrap()
        else {
            panic!("lifted coaxial tori meet along curves");
        };
        assert_eq!(curves.len(), 2);
        for section in &curves {
            assert!(section.exact);
            assert_same_parameter(section, &ring, section.on_a.as_ref().unwrap(), 48);
            assert_same_parameter(section, &lifted, section.on_b.as_ref().unwrap(), 48);
        }
    }

    #[test]
    fn a_pair_with_no_closed_form_comes_back_fitted_with_pcurves() {
        // Crossed cylinders: the marched path, end to end through one call.
        let a = cylinder(Vector::Z, 1.0);
        let b = cylinder(Vector::X, 1.6);
        let options = IntersectOptions {
            tolerance: 1e-5,
            marching: Marching {
                chord: 1e-5,
                ..Marching::default()
            },
        };
        let SurfaceIntersection::Along(curves) = intersect_surfaces(&a, &b, options, T).unwrap()
        else {
            panic!("crossed cylinders meet along curves");
        };
        assert_eq!(curves.len(), 2);
        for section in &curves {
            assert!(!section.exact);
            assert!(section.closed);
            assert!(
                section.tolerance <= 1e-5 + 1e-4,
                "got {}",
                section.tolerance
            );
            assert!(section.on_a.is_some() && section.on_b.is_some());

            // The fitted curve lies on both cylinders to its stated tolerance.
            let (lo, hi) = section.curve.domain();
            for i in 0..=200 {
                #[allow(clippy::cast_precision_loss)]
                let t = lo + (hi - lo) * f64::from(i) / 200.0;
                let p = section.curve.point_at(t, T).unwrap();
                let (SurfaceGeometry::Cylinder(x), SurfaceGeometry::Cylinder(y)) = (&a, &b) else {
                    unreachable!()
                };
                let off = x
                    .cylinder()
                    .distance_to(p)
                    .abs()
                    .max(y.cylinder().distance_to(p).abs());
                assert!(
                    off <= section.tolerance * 2.0,
                    "at t = {t} the fitted curve is {off:e} off, tolerance {}",
                    section.tolerance
                );
            }
        }
    }

    #[test]
    fn exact_lines_are_clipped_to_the_surfaces_extents() {
        // The analytic layer answers for the unbounded geometry; the surfaces
        // are finite. A section line a billion units long is not something an
        // edge can be built on, and one wholly outside the extents is a
        // phantom.
        let drum = cylinder(Vector::Z, 2.0);
        let cut = plane(Point::ORIGIN, Vector::X);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&drum, &cut, IntersectOptions::default(), T).unwrap()
        else {
            panic!("expected curves");
        };
        for section in &curves {
            let (lo, hi) = section.curve.domain();
            // Bounded by the cylinder's height, not by LINE_EXTENT.
            assert!(
                hi - lo <= 8.0 + 1e-9,
                "the line was not clipped: [{lo}, {hi}]"
            );
            let start = section.curve.point_at(lo, T).unwrap();
            let end = section.curve.point_at(hi, T).unwrap();
            assert!(start.z >= -4.0 - 1e-9 && end.z <= 4.0 + 1e-9);
        }

        // A circle at a height the bounded cylinder does not reach is not an
        // intersection of these surfaces, however truly the unbounded ones
        // meet there.
        let high = plane(Point::new(0.0, 0.0, 10.0), Vector::Z);
        assert_eq!(
            intersect_surfaces(&drum, &high, IntersectOptions::default(), T).unwrap(),
            SurfaceIntersection::Apart
        );
    }

    #[test]
    fn the_degenerate_answers_pass_through() {
        assert_eq!(
            intersect_surfaces(
                &sphere(Point::ORIGIN, 1.0),
                &sphere(Point::new(5.0, 0.0, 0.0), 1.0),
                IntersectOptions::default(),
                T
            )
            .unwrap(),
            SurfaceIntersection::Apart
        );
        assert_eq!(
            intersect_surfaces(
                &sphere(Point::ORIGIN, 1.0),
                &sphere(Point::ORIGIN, 1.0),
                IntersectOptions::default(),
                T
            )
            .unwrap(),
            SurfaceIntersection::Same
        );
        assert!(matches!(
            intersect_surfaces(
                &plane(Point::ORIGIN, Vector::Z),
                &sphere(Point::new(0.0, 0.0, 2.0), 2.0),
                IntersectOptions::default(),
                T
            )
            .unwrap(),
            SurfaceIntersection::Touching(ref p) if p.len() == 1
        ));
    }

    #[test]
    fn unusable_options_are_refused() {
        let a = sphere(Point::ORIGIN, 1.0);
        let b = plane(Point::ORIGIN, Vector::Z);
        for tolerance in [0.0, -1.0, f64::NAN] {
            let options = IntersectOptions {
                tolerance,
                ..IntersectOptions::default()
            };
            assert!(intersect_surfaces(&a, &b, options, T).is_err());
        }
    }

    #[test]
    fn a_circle_wound_against_the_axis_keeps_its_pcurve_same_parameter() {
        // The winding bug the boolean's drill test found: a plane whose
        // normal opposes the cylinder's axis cuts a circle wound against the
        // cylinder's `u`, and the pcurve must run in `-u` with it. Written
        // `+u` unconditionally, the pcurve evaluated half a turn away from
        // the curve and every face built on the section tore in parameter
        // space. Both windings are pinned by lifting the pcurve through the
        // surface and demanding the curve's own point back.
        let drum: SurfaceGeometry = CylinderSurface::new(
            Cylinder::new(
                Frame::new(Point::new(2.0, 2.0, -1.0), Direction::Z, Direction::X, T).unwrap(),
                0.5,
                T,
            )
            .unwrap(),
            (0.0, 3.0),
        )
        .unwrap()
        .into();
        for normal in [Direction::Z, -Direction::Z] {
            let frame = Frame::new(Point::ORIGIN, normal, Direction::X, T).unwrap();
            let ground: SurfaceGeometry =
                PlaneSurface::over(Plane::new(frame), (-4.0, 4.0), (-4.0, 4.0))
                    .unwrap()
                    .into();
            let met = intersect_surfaces(&ground, &drum, IntersectOptions::default(), T).unwrap();
            let SurfaceIntersection::Along(curves) = met else {
                panic!("a plane through a cylinder sections it");
            };
            for sc in &curves {
                let pcurve = sc
                    .on_b
                    .as_ref()
                    .expect("a circle on its cylinder has a pcurve");
                let (lo, hi) = sc.curve.domain();
                for i in 0..8 {
                    let t = lo + (hi - lo) * f64::from(i) / 8.0;
                    let p3 = sc.curve.point_at(t, T).unwrap();
                    let uv = pcurve.point_at(t, T).unwrap();
                    let lifted = drum
                        .point_at(uv.x.rem_euclid(core::f64::consts::TAU), uv.y, T)
                        .unwrap();
                    assert!(
                        p3.distance(lifted) < 1e-9,
                        "normal {normal:?}, t {t}: pcurve lifts {lifted:?} against {p3:?}"
                    );
                }
            }
        }
    }

    /// A plane through a ball's own axis cuts a meridian. The whole circle has
    /// no chart image — its longitude jumps half a turn at each pole — but
    /// each half is a straight line in the chart, exactly, at the circle's own
    /// parameter. Pinned by lifting the line back through the sphere and
    /// demanding the circle's point, on every half of every orientation.
    #[test]
    fn a_meridian_half_has_an_exact_line_for_a_pcurve() {
        use ogeom_geom::Surface as _;
        let half = core::f64::consts::PI;
        for (centre, radius) in [(Point::ORIGIN, 4.0), (Point::new(1.0, -2.0, 0.5), 1.25)] {
            let ball = sphere(centre, radius);
            let SurfaceGeometry::Sphere(s) = &ball else {
                panic!("a sphere surface");
            };
            // Three planes through the axis, at different azimuths, so the
            // constant longitude is not accidentally zero.
            for azimuth in [0.0_f64, 0.7, 2.4] {
                let normal = Vector::new(-azimuth.sin(), azimuth.cos(), 0.0);
                let cut = plane(centre, normal);
                let SurfaceIntersection::Along(curves) =
                    intersect_surfaces(&ball, &cut, IntersectOptions::default(), T).unwrap()
                else {
                    panic!("a plane through the centre meets the ball along a circle");
                };
                assert_eq!(curves.len(), 1, "one great circle");
                let circle = &curves[0].curve;
                assert!(curves[0].exact);
                // The whole circle has no chart image; each half does.
                assert!(
                    exact_pcurve_over(circle, circle.domain(), &ball, T).is_none(),
                    "the whole meridian has no single chart image"
                );
                for (lo, hi) in [(0.0, half), (half, 2.0 * half), (0.3, half - 0.1)] {
                    let pcurve = exact_pcurve_over(circle, (lo, hi), &ball, T)
                        .expect("half a meridian has an exact pcurve");
                    assert!(
                        matches!(pcurve, PlanarCurve::Line(_)),
                        "and it is a straight line in the chart"
                    );
                    for i in 0..=16 {
                        let t = (hi - lo).mul_add(f64::from(i) / 16.0, lo);
                        let want = circle.point_at(t, T).unwrap();
                        let uv = pcurve.point_at(t, T).unwrap();
                        assert!(
                            uv.y >= -half.mul_add(0.5, 1e-12) && uv.y <= half.mul_add(0.5, 1e-12),
                            "the latitude stays inside the chart: {}",
                            uv.y
                        );
                        let lifted = ball
                            .point_at(uv.x.rem_euclid(core::f64::consts::TAU), uv.y, T)
                            .unwrap();
                        assert!(
                            want.distance(lifted) < 1e-9,
                            "azimuth {azimuth}, t {t}: {lifted:?} against {want:?}"
                        );
                    }
                }
                // A range straddling a pole has none, and says so rather than
                // answering for one side.
                assert!(
                    exact_pcurve_over(circle, (half - 0.2, half + 0.2), &ball, T).is_none(),
                    "a range across a pole has no one line"
                );
                let _ = s;
            }
        }
    }
}
