//! Where two surfaces meet: the one call.
//!
//! Everything else in this crate is a stage — closed forms, seeding, tracing,
//! fitting. This is the function an application calls, and the one `og-bool`
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

use og_core::{OgResult, Tolerances, og_bail};
use og_geom::{Circle2d, Curve, Curve3d, Ellipse2d, Line2d, PlanarCurve, Surface, SurfaceGeometry};
use og_math::{Circle2, Ellipse2, Frame2, Point, Point2};

use crate::approx::approximate_branch;
use crate::march::{Marching, branches};
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
/// [`OgError::Construction`](og_core::OgError::Construction) if the options
/// are unusable. A pair the marcher finds nothing for is [`Apart`], not an
/// error — see that variant for what it can and cannot claim.
///
/// [`Apart`]: SurfaceIntersection::Apart
pub fn intersect_surfaces(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: IntersectOptions,
    tol: Tolerances,
) -> OgResult<SurfaceIntersection> {
    if !options.tolerance.is_finite() || options.tolerance <= 0.0 {
        og_bail!(
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
    let on_a = exact_pcurve(&curve, a, tol);
    let on_b = exact_pcurve(&curve, b, tol);

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
        let clipped: Curve = og_geom::LineCurve::over(line.axis(), lo, hi).ok()?.into();
        let clip2 = |p: &PlanarCurve| -> Option<PlanarCurve> {
            let PlanarCurve::Line(l) = p else {
                return Some(p.clone());
            };
            Some(Line2d::over(l.axis(), lo, hi).ok()?.into())
        };
        return Some(SectionCurve {
            on_a: on_a.as_ref().and_then(clip2),
            on_b: on_b.as_ref().and_then(clip2),
            tolerance: 0.0,
            exact: true,
            closed: false,
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
    Some(SectionCurve {
        on_a,
        on_b,
        tolerance: 0.0,
        exact: true,
        closed,
        curve,
    })
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
    use og_geom::Curve2d;
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
) -> OgResult<SurfaceIntersection> {
    let traced = branches(a, b, options.marching, tol)?;
    if traced.is_empty() {
        return Ok(SurfaceIntersection::Apart);
    }
    let mut out = Vec::with_capacity(traced.len());
    for branch in &traced {
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
        });
    }
    Ok(SurfaceIntersection::Along(out))
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
fn exact_pcurve(curve: &Curve, surface: &SurfaceGeometry, tol: Tolerances) -> Option<PlanarCurve> {
    match surface {
        SurfaceGeometry::Plane(p) => on_plane(curve, p.plane(), tol),
        SurfaceGeometry::Cylinder(c) => on_cylinder(curve, c.cylinder(), tol),
        SurfaceGeometry::Sphere(s) => on_sphere(curve, s.sphere(), tol),
        _ => None,
    }
}

/// Project a curve lying in a plane into the plane's own coordinates.
///
/// Exact for a line, a circle and an ellipse: the plane's frame is orthonormal,
/// so lengths and the curves' own parameterizations survive the projection
/// unchanged.
fn on_plane(curve: &Curve, plane: og_math::Plane, tol: Tolerances) -> Option<PlanarCurve> {
    let frame = plane.frame();
    let flat = |p: Point| {
        let local = frame.to_local(p);
        Point2::new(local.x, local.y)
    };
    let flat_direction = |d: og_math::Direction| {
        let tip = flat(frame.origin() + d.vector());
        og_math::Direction2::new(tip - flat(frame.origin()), tol).ok()
    };
    match curve {
        Curve::Line(line) => {
            let axis = line.axis();
            let through = flat(axis.location);
            let direction = flat_direction(axis.direction)?;
            let (lo, hi) = line.domain();
            Some(
                Line2d::over(og_math::Axis2::new(through, direction), lo, hi)
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
        _ => None,
    }
}

/// The pcurve of a curve on a cylinder, where it is a straight line in
/// parameter space.
///
/// A line along the axis runs at constant `u`; a full circle around it runs at
/// constant `v`. Both are lines in `(u, v)`, exactly, and both inherit the 3D
/// curve's own parameter — height for the line, angle for the circle.
fn on_cylinder(curve: &Curve, cylinder: og_math::Cylinder, tol: Tolerances) -> Option<PlanarCurve> {
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
                og_math::Direction2::new(og_math::Vector2::new(0.0, along.signum()), tol).ok()?;
            Some(
                Line2d::over(og_math::Axis2::new(start, towards), lo, hi)
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
            // Where the circle's own angle zero sits in the cylinder's angle.
            let start = circle.centre() + circle.frame().x().vector() * circle.radius();
            let at = frame.to_local(start);
            let phase = at.y.atan2(at.x);
            let towards = og_math::Direction2::new(og_math::Vector2::new(1.0, 0.0), tol).ok()?;
            Some(
                Line2d::over(
                    og_math::Axis2::new(Point2::new(phase, local.z), towards),
                    0.0,
                    core::f64::consts::TAU,
                )
                .ok()?
                .into(),
            )
        }
        _ => None,
    }
}

/// The pcurve of a circle of latitude on a sphere.
fn on_sphere(curve: &Curve, sphere: og_math::Sphere, tol: Tolerances) -> Option<PlanarCurve> {
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
        return None;
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
    let towards = og_math::Direction2::new(og_math::Vector2::new(1.0, 0.0), tol).ok()?;
    Some(
        Line2d::over(
            og_math::Axis2::new(Point2::new(phase, latitude), towards),
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
    use og_geom::{Curve2d, Curve3d, CylinderSurface, PlaneSurface, SphereSurface};
    use og_math::{Cylinder, Direction, Frame, Plane, Sphere, Vector};

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
}
