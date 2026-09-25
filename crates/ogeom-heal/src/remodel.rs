//! Remodelling: a shape's geometry restated wholesale, every trim
//! re-derived against what replaced its surface.
//!
//! [`restrict_degree`] brings every spline curve and surface under a degree
//! limit at a stated tolerance, for a format that caps it.
//! [`swept_to_elementary`] names the swept surfaces that are planes, drums,
//! cones, balls or tori as what they are: an extruded line or circle, a
//! revolved line or circle.

use ogeom_algo::{Built, project_on_surface, restate_geometry};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    ConeSurface, Curve, CylinderSurface, PlaneSurface, SphereSurface, Surface as _,
    SurfaceGeometry, TorusSurface,
};
use ogeom_math::{Axis, Cone, Cylinder, Direction, Frame, Plane, Point, Sphere, Torus, Vector};
use ogeom_topo::{Model, Shape};

/// `shape` with every spline curve and surface of degree above
/// `max_degree` refitted at that degree, each within `tolerance` of what
/// it replaces, and every trim re-derived. Analytic geometry has no degree
/// to restrict and stays as it is.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `max_degree` is zero or `tolerance` is not a distance;
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if a fit cannot
/// reach `tolerance`; as [`ogeom_algo::restate_geometry`].
pub fn restrict_degree(
    model: &mut Model,
    shape: &Shape,
    max_degree: usize,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if max_degree == 0 {
        ogeom_bail!(Construction, "a degree limit of zero leaves nothing");
    }
    if !(tolerance.is_finite() && tolerance > 0.0) {
        ogeom_bail!(Construction, "a tolerance of {tolerance} is not a distance");
    }
    let surface = |s: &SurfaceGeometry| -> OgeomResult<Option<(SurfaceGeometry, bool)>> {
        let SurfaceGeometry::BSpline(b) = s else {
            return Ok(None);
        };
        if b.u_knots().degree() <= max_degree && b.v_knots().degree() <= max_degree {
            return Ok(None);
        }
        let fitted = b.restricted_to_degree(max_degree, tolerance, tol)?;
        if !fitted.met {
            ogeom_bail!(
                NotDone,
                "a patch fitted at degree {max_degree} stays {} away",
                fitted.error
            );
        }
        let restated = SurfaceGeometry::BSpline(fitted.curve);
        let flipped = turned(s, &restated, tol)?;
        Ok(Some((restated, flipped)))
    };
    let curve = |c: &Curve, range: (f64, f64)| -> OgeomResult<Option<(Curve, (f64, f64))>> {
        let Curve::BSpline(b) = c else {
            return Ok(None);
        };
        if b.degree() <= max_degree {
            return Ok(None);
        }
        // The edge's piece, open: a ring's spline has no ends to fit
        // between, and its piece does.
        let piece = b.segment(range, tol)?;
        let fitted = piece.restricted_to_degree(max_degree, tolerance, tol)?;
        if !fitted.met {
            ogeom_bail!(
                NotDone,
                "a curve fitted at degree {max_degree} stays {} away",
                fitted.error
            );
        }
        let restated = Curve::BSpline(fitted.curve);
        let domain = ogeom_geom::Curve3d::domain(&restated);
        Ok(Some((restated, domain)))
    };
    restate_geometry(model, shape, &surface, &curve, tol)
}

/// `shape` with every surface of extrusion or revolution that is a plane,
/// drum, cone, ball or torus restated as it: a line or circle extruded, a
/// line or circle revolved in a plane through the axis. Every trim is
/// re-derived; the surfaces are the same point sets, so nothing moves.
///
/// # Errors
///
/// As [`ogeom_algo::restate_geometry`].
pub fn swept_to_elementary(
    model: &mut Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let surface = |s: &SurfaceGeometry| -> OgeomResult<Option<(SurfaceGeometry, bool)>> {
        let Some(candidate) = elementary(s, tol)? else {
            return Ok(None);
        };
        // Held to every sample: a wrong yes is a solid with the wrong
        // surface under every later operation.
        for p in interior_samples(s, 5, tol)? {
            if project_on_surface(&candidate, p, 16, tol)?.distance > tol.confusion() {
                return Ok(None);
            }
        }
        let flipped = turned(s, &candidate, tol)?;
        Ok(Some((candidate, flipped)))
    };
    let curve = |_: &Curve, _: (f64, f64)| -> OgeomResult<Option<(Curve, (f64, f64))>> { Ok(None) };
    restate_geometry(model, shape, &surface, &curve, tol)
}

/// The elementary surface a swept one proposes, unverified.
fn elementary(s: &SurfaceGeometry, tol: Tolerances) -> OgeomResult<Option<SurfaceGeometry>> {
    let samples = interior_samples(s, 5, tol)?;
    let unwrap = |c: &Curve| -> Curve {
        match c {
            Curve::Trimmed(t) => t.basis().clone(),
            other => other.clone(),
        }
    };
    let found = match s {
        SurfaceGeometry::Extrusion(e) => match unwrap(e.curve()) {
            Curve::Line(l) => {
                let axis = l.axis();
                let normal = axis.direction.vector().cross(e.direction().vector());
                let Ok(z) = Direction::new(normal, tol) else {
                    return Ok(None);
                };
                let frame = Frame::new(axis.location, z, axis.direction, tol)?;
                Some(plane_over(frame, &samples)?)
            }
            Curve::Circle(c) => {
                let circle = c.circle();
                let frame = circle.frame();
                if frame.z().vector().cross(e.direction().vector()).magnitude() > tol.angular() {
                    return Ok(None);
                }
                let axis = Frame::new(frame.origin(), e.direction(), frame.x(), tol)?;
                let heights = heights(&axis, &samples);
                Some(
                    CylinderSurface::new(Cylinder::new(axis, circle.radius(), tol)?, heights)?
                        .into(),
                )
            }
            _ => None,
        },
        SurfaceGeometry::Revolution(r) => {
            let axis = r.axis();
            match unwrap(r.curve()) {
                Curve::Line(l) => revolved_line(axis, l.axis(), &samples, tol)?,
                Curve::Circle(c) => revolved_circle(axis, c.circle(), tol)?,
                _ => None,
            }
        }
        _ => None,
    };
    Ok(found)
}

/// A line revolved in a plane through the axis: a drum where it runs with
/// the axis, a disc where it runs square to it, a cone otherwise.
fn revolved_line(
    axis: Axis,
    line: Axis,
    samples: &[Point],
    tol: Tolerances,
) -> OgeomResult<Option<SurfaceGeometry>> {
    let z = axis.direction.vector();
    let a = line.location;
    let b = a + line.direction.vector();
    let height = |p: Point| (p - axis.location).dot(z);
    let radius = |p: Point| p.distance(axis.project(p));
    // The profile's own start side names the chart's `u = 0`.
    let Some(off) = [a, b]
        .into_iter()
        .map(|p| p - axis.project(p))
        .find(|v| v.magnitude() > tol.confusion())
    else {
        return Ok(None);
    };
    let x = Direction::new(off, tol)?;
    let (dh, dr) = (height(b) - height(a), radius(b) - radius(a));
    if dh.abs() <= tol.confusion() {
        let frame = Frame::new(axis.project(a), axis.direction, x, tol)?;
        return Ok(Some(plane_over(frame, samples)?));
    }
    let frame = Frame::new(axis.project(a), axis.direction, x, tol)?;
    let heights = heights(&frame, samples);
    if dr.abs() <= tol.confusion() {
        return Ok(Some(
            CylinderSurface::new(Cylinder::new(frame, radius(a), tol)?, heights)?.into(),
        ));
    }
    let half_angle = (dr / dh).atan();
    let Ok(cone) = Cone::new(frame, radius(a), half_angle, tol) else {
        return Ok(None);
    };
    Ok(Some(ConeSurface::new(cone, heights)?.into()))
}

/// A circle revolved in a plane through the axis: a ball where its centre
/// is on the axis, a torus otherwise.
fn revolved_circle(
    axis: Axis,
    circle: ogeom_math::Circle,
    tol: Tolerances,
) -> OgeomResult<Option<SurfaceGeometry>> {
    let centre = circle.centre();
    let foot = axis.project(centre);
    let off = centre - foot;
    if off.magnitude() <= tol.confusion() {
        let x = Direction::new(circle.frame().x().vector(), tol)?;
        let x = if x.vector().cross(axis.direction.vector()).magnitude() > tol.angular() {
            x
        } else {
            Direction::new(circle.frame().y().vector(), tol)?
        };
        let side = x.vector() - axis.direction.vector() * x.vector().dot(axis.direction.vector());
        let frame = Frame::new(foot, axis.direction, Direction::new(side, tol)?, tol)?;
        return Ok(Some(
            SphereSurface::new(Sphere::new(frame, circle.radius(), tol)?).into(),
        ));
    }
    let frame = Frame::new(foot, axis.direction, Direction::new(off, tol)?, tol)?;
    let Ok(torus) = Torus::new(frame, off.magnitude(), circle.radius(), tol) else {
        return Ok(None);
    };
    Ok(Some(TorusSurface::new(torus).into()))
}

/// A plane in `frame`, bounded round the samples.
fn plane_over(frame: Frame, samples: &[Point]) -> OgeomResult<SurfaceGeometry> {
    let along = |p: Point, d: Vector| (p - frame.origin()).dot(d);
    let range = |d: Vector| {
        let (lo, hi) = samples
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), p| {
                (l.min(along(*p, d)), h.max(along(*p, d)))
            });
        let margin = (hi - lo) * 0.5 + 1.0;
        (lo - margin, hi + margin)
    };
    Ok(PlaneSurface::over(
        Plane::new(frame),
        range(frame.x().vector()),
        range(frame.y().vector()),
    )?
    .into())
}

/// The heights along a frame's `z` that the samples span, with room.
fn heights(frame: &Frame, samples: &[Point]) -> (f64, f64) {
    let z = frame.z().vector();
    let (lo, hi) = samples
        .iter()
        .map(|p| (*p - frame.origin()).dot(z))
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), v| {
            (l.min(v), h.max(v))
        });
    let margin = (hi - lo) * 0.5 + 1.0;
    (lo - margin, hi + margin)
}

/// Points over the inside of a surface's chart, corners included.
fn interior_samples(s: &SurfaceGeometry, n: usize, tol: Tolerances) -> OgeomResult<Vec<Point>> {
    let ((u0, u1), (v0, v1)) = s.domain();
    let mut out = Vec::with_capacity((n + 1) * (n + 1));
    for i in 0..=n {
        for j in 0..=n {
            #[allow(clippy::cast_precision_loss, reason = "a sample index")]
            let (fu, fv) = (i as f64 / n as f64, j as f64 / n as f64);
            out.push(s.point_at(u0 + (u1 - u0) * fu, v0 + (v1 - v0) * fv, tol)?);
        }
    }
    Ok(out)
}

/// Whether `new`'s normal points against `old`'s, where they meet at the
/// middle of `old`'s chart.
fn turned(old: &SurfaceGeometry, new: &SurfaceGeometry, tol: Tolerances) -> OgeomResult<bool> {
    let ((u0, u1), (v0, v1)) = old.domain();
    // Off the exact middle: a revolution's middle can sit on its axis.
    let (u, v) = (u0 + (u1 - u0) * 0.43, v0 + (v1 - v0) * 0.57);
    let p = old.point_at(u, v, tol)?;
    let n_old = old.normal_at(u, v, tol)?;
    let at = project_on_surface(new, p, 16, tol)?.parameters;
    let n_new = new.normal_at(at.0, at.1, tol)?;
    Ok(n_old.vector().dot(n_new.vector()) < 0.0)
}
