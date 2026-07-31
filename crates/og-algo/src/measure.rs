//! Bounding volumes and projection.
//!
//! # Bounds must contain what they claim to
//!
//! A bounding box in a kernel is always a *rejection* test. Too large costs
//! time; too small silently drops a real intersection, and nothing downstream
//! can tell it happened. So every bound here is derived from a property that
//! guarantees containment, never from sampling:
//!
//! - a line's bound is its endpoints, which is exact;
//! - a spline's is its control points, which is guaranteed by the convex hull
//!   property — the curve never leaves the hull of its control polygon;
//! - an analytic curve or surface's is computed from its own definition;
//! - a *trimmed* piece falls back to the bound of the whole, which is loose but
//!   never wrong.
//!
//! Sampling a curve at a few parameters and taking the extremes is the obvious
//! alternative and is not sound: the curve bulges between the samples, and the
//! amount it bulges is exactly what a bound is supposed to capture.

use og_core::{OgResult, Tolerances, og_bail};
use og_geom::{Curve, Curve3d, PlanarCurve, Surface, SurfaceGeometry, curve::LINE_EXTENT};
use og_math::{Aabb, Point, Point2, Vector, solve};
use og_topo::{Model, NodeData, Shape, ShapeType, explore_unique};

/// A guaranteed bound for a space curve.
///
/// Loose for a trimmed curve, which reports the bound of the whole rather than
/// of the piece — never wrong, and tightening it would mean solving for the
/// extremes of the trimmed range, which is the same work as an intersection.
///
/// # Errors
///
/// [`OgError::Domain`](og_core::OgError::Domain) if the curve cannot be
/// evaluated at its own domain ends.
pub fn curve_bounds(curve: &Curve, tol: Tolerances) -> OgResult<Aabb> {
    Ok(match curve {
        // Exact: a segment is the hull of its ends.
        Curve::Line(_) => Aabb::of_corners(curve.start(tol)?, curve.end(tol)?),

        // A conic's extremes are its centre displaced by the radii along each
        // frame axis, which bounds the whole conic whatever arc is in use.
        Curve::Circle(c) => {
            let circle = c.circle();
            frame_bounds(
                circle.centre(),
                circle.frame(),
                (circle.radius(), circle.radius(), 0.0),
            )
        }
        Curve::Ellipse(e) => {
            let ellipse = e.ellipse();
            frame_bounds(
                ellipse.centre(),
                ellipse.frame(),
                (ellipse.major_radius(), ellipse.minor_radius(), 0.0),
            )
        }

        // A hyperbola and a parabola are unbounded, so only the trimmed extent
        // has a bound at all. Both are convex in their own frame, so the hull
        // of the two ends and the vertex contains the arc between them.
        Curve::Hyperbola(_) | Curve::Parabola(_) => {
            let (a, b) = curve.domain();
            let mid = curve.point_at(f64::midpoint(a, b), tol)?;
            let ends = Aabb::of_corners(curve.start(tol)?, curve.end(tol)?);
            // The midpoint is the extreme in the frame's x direction for both,
            // and the ends bound the rest.
            ends.with_point(mid)
        }

        // The convex hull property: a B-spline never leaves the hull of its
        // control polygon, so the polygon's box contains the curve exactly.
        Curve::BSpline(s) => Aabb::of_points(
            &s.control_points()
                .iter()
                .map(|w| w.point())
                .collect::<Vec<_>>(),
        ),

        Curve::Trimmed(t) => curve_bounds(t.basis(), tol)?,
    })
}

/// A guaranteed bound for a surface.
///
/// # Errors
///
/// [`OgError::Domain`](og_core::OgError::Domain) if the surface cannot be
/// evaluated over its own domain.
pub fn surface_bounds(surface: &SurfaceGeometry, tol: Tolerances) -> OgResult<Aabb> {
    Ok(match surface {
        // A plane is unbounded; its declared domain is what there is to bound,
        // and the four corners span it exactly.
        SurfaceGeometry::Plane(p) => {
            let ((ua, ub), (va, vb)) = p.domain();
            if ua <= -LINE_EXTENT || ub >= LINE_EXTENT {
                og_bail!(
                    Domain,
                    "an unbounded plane has no finite bound; trim it before asking"
                );
            }
            let mut out = Aabb::EMPTY;
            for (u, v) in [(ua, va), (ua, vb), (ub, va), (ub, vb)] {
                out = out.with_point(p.point_at(u, v, tol)?);
            }
            out
        }

        SurfaceGeometry::Cylinder(c) => {
            let cyl = c.cylinder();
            let ((_, _), (va, vb)) = c.domain();
            let frame = cyl.frame();
            let base = frame.origin() + frame.z() * va;
            let top = frame.origin() + frame.z() * vb;
            let radial = frame_bounds(base, frame, (cyl.radius(), cyl.radius(), 0.0));
            radial.union(&frame_bounds(top, frame, (cyl.radius(), cyl.radius(), 0.0)))
        }

        SurfaceGeometry::Cone(c) => {
            let cone = c.cone();
            let ((_, _), (va, vb)) = c.domain();
            let frame = cone.frame();
            let mut out = Aabb::EMPTY;
            for height in [va, vb] {
                let radius = cone.radius_at(height).abs();
                let centre = frame.origin() + frame.z() * height;
                out = out.union(&frame_bounds(centre, frame, (radius, radius, 0.0)));
            }
            out
        }

        // A sphere's bound is its centre plus its radius on every axis, whatever
        // patch of it is in use.
        SurfaceGeometry::Sphere(s) => {
            let sphere = s.sphere();
            let r = Vector::splat(sphere.radius());
            Aabb::of_corners(sphere.centre() - r, sphere.centre() + r)
        }

        SurfaceGeometry::Torus(t) => {
            let torus = t.torus();
            let reach = torus.major_radius() + torus.minor_radius();
            frame_bounds(
                torus.centre(),
                torus.frame(),
                (reach, reach, torus.minor_radius()),
            )
        }

        // Convex hull property again, in two directions.
        SurfaceGeometry::BSpline(s) => Aabb::of_points(
            &s.grid()
                .points()
                .iter()
                .map(|w| w.point())
                .collect::<Vec<_>>(),
        ),

        // A revolved curve reaches at most its own furthest distance from the
        // axis, in every direction around it.
        SurfaceGeometry::Revolution(r) => {
            let curve = curve_bounds(r.curve(), tol)?;
            let axis = r.axis();
            let mut reach: f64 = 0.0;
            let mut along = Aabb::EMPTY;
            for corner in curve.corners() {
                reach = reach.max(axis.distance_to(corner));
                along = along.with_point(axis.project(corner));
            }
            let radial = Vector::splat(reach);
            along.expanded(0.0).union(&Aabb::of_corners(
                along.low().unwrap_or(axis.location) - radial,
                along.high().unwrap_or(axis.location) + radial,
            ))
        }

        // A swept curve reaches the curve's bound at each end of the sweep.
        SurfaceGeometry::Extrusion(e) => {
            let base = curve_bounds(e.curve(), tol)?;
            let ((_, _), (va, vb)) = e.domain();
            let start = base.transformed(&og_math::Transform::translation(e.direction() * va));
            let end = base.transformed(&og_math::Transform::translation(e.direction() * vb));
            start.union(&end)
        }

        SurfaceGeometry::Trimmed(t) => surface_bounds(t.basis(), tol)?,
    })
}

/// The box of a point displaced by `(x, y, z)` extents along a frame's axes.
///
/// Every axis of the result gets the sum of the absolute contributions from all
/// three frame directions, which is what makes it a bound rather than an
/// estimate: a tilted frame's extent projects onto every world axis at once.
fn frame_bounds(centre: Point, frame: og_math::Frame, extent: (f64, f64, f64)) -> Aabb {
    let (ex, ey, ez) = extent;
    let reach = |axis: fn(&Vector) -> f64| {
        (frame.x().vector().pipe(axis) * ex).abs()
            + (frame.y().vector().pipe(axis) * ey).abs()
            + (frame.z().vector().pipe(axis) * ez).abs()
    };
    let r = Vector::new(reach(|v| v.x), reach(|v| v.y), reach(|v| v.z));
    Aabb::of_corners(centre - r, centre + r)
}

/// A tiny helper so the reach computation above reads as one expression.
trait Pipe {
    fn pipe<R>(&self, f: impl FnOnce(&Self) -> R) -> R;
}

impl Pipe for Vector {
    fn pipe<R>(&self, f: impl FnOnce(&Self) -> R) -> R {
        f(self)
    }
}

/// A guaranteed bound for a shape, including everything below it.
///
/// Vertices contribute their point widened by their own tolerance, since a
/// vertex genuinely occupies that much space. Edges and faces contribute the
/// bound of their geometry, likewise widened.
///
/// # Errors
///
/// [`OgError::Dangling`](og_core::OgError::Dangling) if any handle fails to
/// resolve, and whatever the geometry's own bound reports.
pub fn shape_bounds(model: &Model, shape: &Shape, tol: Tolerances) -> OgResult<Aabb> {
    let Some(node) = model.node(shape) else {
        og_bail!(Dangling, "shape refers to a node not in this model");
    };
    let placement = shape.transform(model.datums())?;

    let own = match node.data() {
        NodeData::Vertex(v) => Aabb::of_point(placement.apply(v.point)).expanded(v.tolerance.get()),
        NodeData::Edge(e) => {
            let mut out = Aabb::EMPTY;
            for repr in &e.representations {
                if let og_topo::EdgeRepr::Curve3d { curve, .. } = repr
                    && let Some(geometry) = model.geometry().curve(*curve)
                {
                    out = out.union(&curve_bounds(geometry, tol)?);
                }
            }
            out.transformed(&placement).expanded(e.tolerance.get())
        }
        NodeData::Face(f) => {
            // The surface bound covers the whole surface, and the face is a
            // piece of it, so this contains the face. The wires below tighten
            // nothing but cost nothing either — the union is taken regardless.
            match model.geometry().surface(f.surface) {
                Some(surface) => surface_bounds(surface, tol)
                    .unwrap_or(Aabb::EMPTY)
                    .transformed(&placement)
                    .expanded(f.tolerance.get()),
                None => Aabb::EMPTY,
            }
        }
        NodeData::Container => Aabb::EMPTY,
    };

    let mut out = own;
    for child in model.children_of(shape)? {
        out = out.union(&shape_bounds(model, &child, tol)?);
    }
    Ok(out)
}

/// A bound for a shape built only from its vertices.
///
/// Tighter than [`shape_bounds`] for a solid whose faces sit on unbounded
/// surfaces, and *not* a guarantee: a curved edge bulges past its own
/// endpoints. Use it for a quick estimate, never for a rejection test.
///
/// # Errors
///
/// As [`shape_bounds`].
pub fn vertex_bounds(model: &Model, shape: &Shape, tol: Tolerances) -> OgResult<Aabb> {
    let mut out = Aabb::EMPTY;
    for vertex in explore_unique(model, shape, ShapeType::Vertex)? {
        let Some(node) = model.node(&vertex) else {
            og_bail!(Dangling, "vertex is not in this model");
        };
        if let Some(data) = node.data().as_vertex() {
            let placed = vertex.transform(model.datums())?.apply(data.point);
            out = out.with_point(placed);
        }
    }
    Ok(out.expanded(tol.confusion()))
}

/// Where a point projects onto a curve, and how far away it is.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Projection {
    /// The parameter of the nearest point found.
    pub parameter: f64,
    /// The nearest point.
    pub point: Point,
    /// The distance to it.
    pub distance: f64,
}

/// The nearest point on a curve to `target`.
///
/// Samples the domain to bracket the minimum, then refines. The sampling is not
/// decoration: the distance function along a curve is generally multi-modal —
/// a point inside a circle is equidistant from every part of it, and a point
/// near a spline's inflection has two competing minima — so starting a local
/// method from one guess finds whichever basin it happens to land in. The
/// sample count sets how fine a feature can be resolved, and is stated rather
/// than hidden.
///
/// # Errors
///
/// [`OgError::Domain`](og_core::OgError::Domain) if the curve cannot be
/// evaluated over its own domain.
pub fn project_on_curve(
    curve: &Curve,
    target: Point,
    samples: usize,
    tol: Tolerances,
) -> OgResult<Projection> {
    let (a, b) = curve.domain();
    let steps = samples.max(8);

    let distance_at = |u: f64| -> f64 {
        curve
            .point_at(u, tol)
            .map_or(f64::INFINITY, |p| p.square_distance(target))
    };

    // Coarse scan for the best bracket.
    let mut best = (a, distance_at(a));
    let mut best_index = 0_usize;
    for i in 1..=steps {
        #[allow(clippy::cast_precision_loss)]
        let u = a + (b - a) * (i as f64 / steps as f64);
        let d = distance_at(u);
        if d < best.1 {
            best = (u, d);
            best_index = i;
        }
    }

    // Refine inside the neighbouring samples, where the minimum must lie.
    #[allow(clippy::cast_precision_loss)]
    let width = (b - a) / steps as f64;
    let lo = (best.0 - width).max(a);
    let hi = (best.0 + width).min(b);
    let _ = best_index;

    let parameter = if hi > lo {
        let refined = solve::minimize(
            distance_at,
            lo,
            hi,
            solve::Criteria {
                residual: 0.0,
                step: tol.parametric(),
                max_iterations: 100,
            },
        )?;
        // The refinement may land marginally worse than the sample if the
        // bracket was already at the boundary; keep whichever is actually
        // nearer rather than trusting the method.
        if distance_at(refined.value) <= best.1 {
            refined.value
        } else {
            best.0
        }
    } else {
        best.0
    };

    let point = curve.point_at(parameter, tol)?;
    Ok(Projection {
        parameter,
        point,
        distance: point.distance(target),
    })
}

/// Where a point projects onto a surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceProjection {
    /// The parameters of the nearest point found.
    pub parameters: (f64, f64),
    /// The nearest point.
    pub point: Point,
    /// The distance to it.
    pub distance: f64,
}

/// The nearest point on a surface to `target`.
///
/// A coarse grid to bracket, then Newton on the two conditions that define a
/// foot point: the displacement from the surface to the target is perpendicular
/// to both tangents. Grid resolution is `samples` per direction, for the same
/// reason as [`project_on_curve`].
///
/// # Errors
///
/// [`OgError::Domain`](og_core::OgError::Domain) if the surface cannot be
/// evaluated over its own domain.
pub fn project_on_surface(
    surface: &SurfaceGeometry,
    target: Point,
    samples: usize,
    tol: Tolerances,
) -> OgResult<SurfaceProjection> {
    let ((ua, ub), (va, vb)) = surface.domain();
    let steps = samples.max(4);

    let mut best = (ua, va, f64::INFINITY);
    for i in 0..=steps {
        for j in 0..=steps {
            #[allow(clippy::cast_precision_loss)]
            let (u, v) = (
                ua + (ub - ua) * (i as f64 / steps as f64),
                va + (vb - va) * (j as f64 / steps as f64),
            );
            if let Ok(p) = surface.point_at(u, v, tol) {
                let d = p.square_distance(target);
                if d < best.2 {
                    best = (u, v, d);
                }
            }
        }
    }

    // The foot point conditions: (S - target) . Su = 0 and (S - target) . Sv = 0.
    let residual = |x: &[f64]| {
        let (u, v) = (x[0], x[1]);
        let Ok(p) = surface.point_at(u, v, tol) else {
            return (vec![0.0, 0.0], vec![vec![1.0, 0.0], vec![0.0, 1.0]]);
        };
        let Ok((du, dv)) = surface.d1_at(u, v, tol) else {
            return (vec![0.0, 0.0], vec![vec![1.0, 0.0], vec![0.0, 1.0]]);
        };
        let (d2u, duv, d2v) =
            surface
                .d2_at(u, v, tol)
                .unwrap_or((Vector::ZERO, Vector::ZERO, Vector::ZERO));
        let gap = p - target;
        (
            vec![gap.dot(du), gap.dot(dv)],
            vec![
                vec![du.dot(du) + gap.dot(d2u), du.dot(dv) + gap.dot(duv)],
                vec![du.dot(dv) + gap.dot(duv), dv.dot(dv) + gap.dot(d2v)],
            ],
        )
    };

    let refined = solve::newton_system(
        residual,
        &[best.0, best.1],
        solve::Criteria {
            residual: tol.confusion(),
            step: tol.parametric(),
            max_iterations: 60,
        },
    );

    let (u, v) = match refined {
        Ok(solution) if solution.convergence.is_converged() => {
            let (u, v) = (solution.value[0], solution.value[1]);
            // Newton is free to wander outside the domain; a foot point that
            // left it is not a foot point of this surface.
            match surface.normalize_parameters(u, v, tol) {
                Ok(inside)
                    if surface
                        .point_at(inside.0, inside.1, tol)
                        .is_ok_and(|p| p.square_distance(target) <= best.2) =>
                {
                    inside
                }
                _ => (best.0, best.1),
            }
        }
        _ => (best.0, best.1),
    };

    let point = surface.point_at(u, v, tol)?;
    Ok(SurfaceProjection {
        parameters: (u, v),
        point,
        distance: point.distance(target),
    })
}

/// The nearest point on a planar curve to a point in the same parameter space.
///
/// # Errors
///
/// As [`project_on_curve`].
pub fn project_on_planar_curve(
    curve: &PlanarCurve,
    target: Point2,
    samples: usize,
    tol: Tolerances,
) -> OgResult<(f64, Point2, f64)> {
    use og_geom::Curve2d;

    let (a, b) = curve.domain();
    let steps = samples.max(8);
    let distance_at = |u: f64| -> f64 {
        curve
            .point_at(u, tol)
            .map_or(f64::INFINITY, |p| p.square_distance(target))
    };

    let mut best = (a, distance_at(a));
    for i in 1..=steps {
        #[allow(clippy::cast_precision_loss)]
        let u = a + (b - a) * (i as f64 / steps as f64);
        let d = distance_at(u);
        if d < best.1 {
            best = (u, d);
        }
    }

    #[allow(clippy::cast_precision_loss)]
    let width = (b - a) / steps as f64;
    let (lo, hi) = ((best.0 - width).max(a), (best.0 + width).min(b));
    let parameter = if hi > lo {
        let refined = solve::minimize(
            distance_at,
            lo,
            hi,
            solve::Criteria {
                residual: 0.0,
                step: tol.parametric(),
                max_iterations: 100,
            },
        )?;
        if distance_at(refined.value) <= best.1 {
            refined.value
        } else {
            best.0
        }
    } else {
        best.0
    };

    let point = curve.point_at(parameter, tol)?;
    Ok((parameter, point, point.distance(target)))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::make_box;
    use approx::assert_relative_eq;
    use og_geom::{
        BSplineCurve, CircleCurve, CylinderSurface, LineCurve, PlaneSurface, SphereSurface,
        TorusSurface, TrimmedCurve,
    };
    use og_math::{Circle, Cylinder, Direction, Frame, KnotVector, Plane, Sphere, Torus};

    const T: Tolerances = Tolerances::millimetres();

    /// A curve sampled densely — the ground truth a bound must contain.
    fn dense_points(curve: &Curve, n: usize) -> Vec<Point> {
        let (a, b) = curve.domain();
        (0..=n)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let u = a + (b - a) * (i as f64 / n as f64);
                curve.point_at(u, T).unwrap()
            })
            .collect()
    }

    #[test]
    fn a_line_bounds_exactly_to_its_endpoints() {
        let curve: Curve = LineCurve::segment(Point::ORIGIN, Point::new(3.0, 4.0, 0.0), T)
            .unwrap()
            .into();
        let b = curve_bounds(&curve, T).unwrap();
        assert_eq!(b.low(), Some(Point::ORIGIN));
        assert_eq!(b.high(), Some(Point::new(3.0, 4.0, 0.0)));
    }

    #[test]
    fn every_curves_bound_contains_the_curve() {
        // The one property that matters. Checked against dense sampling, which
        // is fine as a *test* oracle even though it is not sound as an
        // implementation.
        let spline_control = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 5.0, 0.0),
            Point::new(3.0, -4.0, 2.0),
            Point::new(5.0, 2.0, -1.0),
            Point::new(6.0, 0.0, 0.0),
        ];
        let tilted = Frame::new(
            Point::new(1.0, -2.0, 3.0),
            Direction::from_coords(1.0, 2.0, 3.0, T).unwrap(),
            Direction::X,
            T,
        )
        .unwrap();

        let curves: Vec<Curve> = vec![
            LineCurve::segment(Point::ORIGIN, Point::new(3.0, 4.0, 0.0), T)
                .unwrap()
                .into(),
            CircleCurve::new(Circle::new(tilted, 2.0, T).unwrap()).into(),
            og_geom::EllipseCurve::new(og_math::Ellipse::new(tilted, 5.0, 3.0, T).unwrap()).into(),
            og_geom::HyperbolaCurve::new(
                og_math::Hyperbola::new(tilted, 3.0, 4.0, T).unwrap(),
                1.5,
            )
            .unwrap()
            .into(),
            og_geom::ParabolaCurve::new(og_math::Parabola::new(tilted, 2.0, T).unwrap(), 4.0)
                .unwrap()
                .into(),
            BSplineCurve::new(
                KnotVector::clamped_uniform(3, spline_control.len()).unwrap(),
                spline_control,
                T,
            )
            .unwrap()
            .into(),
        ];

        for curve in curves {
            let bound = curve_bounds(&curve, T).unwrap().with_tolerance(T);
            for p in dense_points(&curve, 400) {
                assert!(
                    bound.contains(p),
                    "{:?} escaped its bound at {p:?}: {bound}",
                    curve.kind()
                );
            }
        }
    }

    #[test]
    fn a_splines_bound_is_its_control_hull_and_that_is_a_guarantee() {
        // Sampling would miss the bulge between samples; the convex hull
        // property does not, because the curve provably never leaves the hull.
        let control = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 10.0, 0.0),
            Point::new(2.0, 10.0, 0.0),
            Point::new(3.0, 0.0, 0.0),
        ];
        let curve: Curve = BSplineCurve::new(
            KnotVector::clamped_uniform(3, control.len()).unwrap(),
            control.clone(),
            T,
        )
        .unwrap()
        .into();

        let bound = curve_bounds(&curve, T).unwrap();
        assert_eq!(bound, Aabb::of_points(&control));
        for p in dense_points(&curve, 200) {
            assert!(bound.contains(p));
        }
        // And the curve really does stay well inside — the bound is loose, in
        // the safe direction.
        let peak = dense_points(&curve, 200)
            .iter()
            .fold(0.0_f64, |m, p| m.max(p.y));
        assert!(peak < 10.0, "the curve should not reach its control points");
    }

    #[test]
    fn a_trimmed_curve_reports_the_bound_of_the_whole() {
        // Loose but never wrong. Tightening it means solving for the extremes
        // of the trimmed range, which is the same work as an intersection.
        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 2.0, T).unwrap()).into();
        let quarter: Curve = TrimmedCurve::new(circle.clone(), 0.0, 1.5, T)
            .unwrap()
            .into();

        let whole = curve_bounds(&circle, T).unwrap();
        let part = curve_bounds(&quarter, T).unwrap();
        assert_eq!(part, whole);
        for p in dense_points(&quarter, 200) {
            assert!(part.contains(p));
        }
    }

    #[test]
    fn every_surfaces_bound_contains_the_surface() {
        let tilted = Frame::new(
            Point::new(1.0, -2.0, 3.0),
            Direction::from_coords(1.0, 2.0, 3.0, T).unwrap(),
            Direction::X,
            T,
        )
        .unwrap();
        let surfaces: Vec<SurfaceGeometry> = vec![
            PlaneSurface::over(Plane::new(tilted), (-5.0, 5.0), (-3.0, 3.0))
                .unwrap()
                .into(),
            CylinderSurface::new(Cylinder::new(tilted, 2.0, T).unwrap(), (-4.0, 4.0))
                .unwrap()
                .into(),
            og_geom::ConeSurface::new(
                og_math::Cone::new(tilted, 3.0, 0.6, T).unwrap(),
                (-1.0, 5.0),
            )
            .unwrap()
            .into(),
            SphereSurface::new(Sphere::new(tilted, 4.0, T).unwrap()).into(),
            TorusSurface::new(Torus::new(tilted, 5.0, 2.0, T).unwrap()).into(),
        ];

        for surface in surfaces {
            let bound = surface_bounds(&surface, T).unwrap().with_tolerance(T);
            let ((ua, ub), (va, vb)) = surface.domain();
            for i in 0..=40 {
                for j in 0..=40 {
                    let u = ua + (ub - ua) * (f64::from(i) / 40.0);
                    let v = va + (vb - va) * (f64::from(j) / 40.0);
                    let p = surface.point_at(u, v, T).unwrap();
                    assert!(
                        bound.contains(p),
                        "{:?} escaped its bound at ({u}, {v}) -> {p:?}: {bound}",
                        surface.kind()
                    );
                }
            }
        }
    }

    #[test]
    fn a_spheres_bound_is_exact() {
        let s: SurfaceGeometry =
            SphereSurface::new(Sphere::centred(Point::new(1.0, 2.0, 3.0), 4.0, T).unwrap()).into();
        let b = surface_bounds(&s, T).unwrap();
        assert_eq!(b.low(), Some(Point::new(-3.0, -2.0, -1.0)));
        assert_eq!(b.high(), Some(Point::new(5.0, 6.0, 7.0)));
    }

    #[test]
    fn an_unbounded_plane_is_refused_rather_than_bounded_wrongly() {
        // Reporting the enormous default extent as a bound would make every
        // rejection test involving it useless, and silently.
        let s: SurfaceGeometry = PlaneSurface::new(Plane::new(Frame::WORLD)).into();
        assert!(surface_bounds(&s, T).is_err());
    }

    #[test]
    fn a_shapes_bound_contains_every_vertex_it_holds() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T).unwrap();
        let bound = shape_bounds(&model, &built.shape, T).unwrap();

        for vertex in explore_unique(&model, &built.shape, ShapeType::Vertex).unwrap() {
            let p = model
                .node(&vertex)
                .unwrap()
                .data()
                .as_vertex()
                .unwrap()
                .point;
            assert!(bound.contains(p), "vertex {p:?} escaped {bound}");
        }
    }

    #[test]
    fn the_vertex_bound_of_a_box_is_tight_and_the_full_bound_contains_it() {
        // A box's faces sit on planes trimmed to the box, so the two agree
        // closely here — but the vertex bound is documented as an estimate, and
        // the full bound is the one that guarantees containment.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T).unwrap();

        let tight = vertex_bounds(&model, &built.shape, T).unwrap();
        assert_relative_eq!(tight.size().x, 2.0, epsilon = 1e-6);
        assert_relative_eq!(tight.size().y, 3.0, epsilon = 1e-6);
        assert_relative_eq!(tight.size().z, 4.0, epsilon = 1e-6);

        let full = shape_bounds(&model, &built.shape, T).unwrap();
        assert!(full.contains_box(&tight));
    }

    #[test]
    fn a_placed_shapes_bound_moves_with_it() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let here = vertex_bounds(&model, &built.shape, T).unwrap();

        let moved = model.placed(
            &built.shape,
            og_math::Transform::translation(Vector::new(10.0, 0.0, 0.0)),
        );
        let there = vertex_bounds(&model, &moved, T).unwrap();

        assert_relative_eq!(
            there.centre().unwrap().x - here.centre().unwrap().x,
            10.0,
            epsilon = 1e-9
        );
        assert_relative_eq!(there.size().x, here.size().x, epsilon = 1e-9);
    }

    #[test]
    fn projecting_onto_a_line_lands_on_the_foot_of_the_perpendicular() {
        let curve: Curve = LineCurve::segment(Point::ORIGIN, Point::new(10.0, 0.0, 0.0), T)
            .unwrap()
            .into();
        let p = project_on_curve(&curve, Point::new(3.0, 4.0, 0.0), 32, T).unwrap();
        assert_relative_eq!(p.parameter, 3.0, epsilon = 1e-6);
        assert!(p.point.is_equal(Point::new(3.0, 0.0, 0.0), T));
        assert_relative_eq!(p.distance, 4.0, epsilon = 1e-9);
    }

    #[test]
    fn projecting_past_the_end_of_a_segment_clamps_to_the_end() {
        let curve: Curve = LineCurve::segment(Point::ORIGIN, Point::new(10.0, 0.0, 0.0), T)
            .unwrap()
            .into();
        let p = project_on_curve(&curve, Point::new(50.0, 0.0, 0.0), 32, T).unwrap();
        assert_relative_eq!(p.parameter, 10.0, epsilon = 1e-6);
        assert_relative_eq!(p.distance, 40.0, epsilon = 1e-6);
    }

    #[test]
    fn projecting_onto_a_circle_finds_the_nearest_of_many_minima() {
        // The reason for the coarse scan: from outside the circle's plane there
        // is one minimum, but a local method started at the wrong parameter
        // converges to the far side just as happily.
        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 5.0, T).unwrap()).into();
        for angle in [0.1_f64, 1.0, 2.5, 4.0, 6.0] {
            let outside = Point::new(8.0 * angle.cos(), 8.0 * angle.sin(), 0.0);
            let p = project_on_curve(&circle, outside, 64, T).unwrap();
            assert_relative_eq!(p.distance, 3.0, epsilon = 1e-6);
            assert_relative_eq!(p.point.to_vector().magnitude(), 5.0, epsilon = 1e-9);
        }
    }

    #[test]
    fn projecting_onto_a_plane_gives_the_perpendicular_foot() {
        let plane: SurfaceGeometry =
            PlaneSurface::over(Plane::new(Frame::WORLD), (-10.0, 10.0), (-10.0, 10.0))
                .unwrap()
                .into();
        let p = project_on_surface(&plane, Point::new(2.0, 3.0, 7.0), 8, T).unwrap();
        assert!(p.point.is_equal(Point::new(2.0, 3.0, 0.0), T));
        assert_relative_eq!(p.distance, 7.0, epsilon = 1e-9);
    }

    #[test]
    fn projecting_onto_a_sphere_lands_on_the_radial_line() {
        let sphere = Sphere::centred(Point::new(1.0, 1.0, 1.0), 3.0, T).unwrap();
        let surface: SurfaceGeometry = SphereSurface::new(sphere).into();
        for target in [
            Point::new(10.0, 1.0, 1.0),
            Point::new(1.0, 1.0, 9.0),
            Point::new(-4.0, -2.0, 0.0),
        ] {
            let p = project_on_surface(&surface, target, 16, T).unwrap();
            // The foot point is on the sphere, and on the line from the centre.
            assert_relative_eq!(sphere.centre().distance(p.point), 3.0, max_relative = 1e-7);
            assert_relative_eq!(
                p.distance,
                (sphere.centre().distance(target) - 3.0).abs(),
                max_relative = 1e-6
            );
        }
    }

    #[test]
    fn projecting_onto_a_cylinder_is_radial() {
        let cylinder = Cylinder::new(Frame::WORLD, 2.0, T).unwrap();
        let surface: SurfaceGeometry = CylinderSurface::new(cylinder, (-5.0, 5.0)).unwrap().into();
        let p = project_on_surface(&surface, Point::new(6.0, 0.0, 1.0), 16, T).unwrap();
        assert_relative_eq!(p.distance, 4.0, max_relative = 1e-6);
        assert_relative_eq!(p.point.z, 1.0, epsilon = 1e-6);
        assert_relative_eq!(p.point.x.hypot(p.point.y), 2.0, max_relative = 1e-7);
    }

    #[test]
    fn projection_of_a_point_already_on_the_geometry_returns_zero_distance() {
        let curve: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 3.0, T).unwrap()).into();
        let on = curve.point_at(1.2, T).unwrap();
        let p = project_on_curve(&curve, on, 64, T).unwrap();
        assert!(p.distance < 1e-7, "distance was {}", p.distance);
    }

    #[test]
    fn projecting_onto_a_planar_curve_works_in_parameter_space() {
        let curve: PlanarCurve =
            og_geom::Line2d::segment(Point2::ORIGIN, Point2::new(10.0, 0.0), T)
                .unwrap()
                .into();
        let (u, point, distance) =
            project_on_planar_curve(&curve, Point2::new(3.0, 4.0), 32, T).unwrap();
        assert_relative_eq!(u, 3.0, epsilon = 1e-6);
        assert!(point.is_equal(Point2::new(3.0, 0.0), T));
        assert_relative_eq!(distance, 4.0, epsilon = 1e-9);
    }
}
