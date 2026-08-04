//! Pipes and lofts: the sweeps whose surfaces are already in the vocabulary.
//!
//! A circular profile along a straight spine is a cylinder; around a full
//! circle it is a torus; along an arc it is a torus segment, built from two
//! half-tube patches and two meridian caps so that every edge is a circle
//! with a closed-form chart. A ruled loft between two parallel sections is
//! walls of planes and cones: segment to segment gives the planar quad,
//! coaxial circle to circle gives the frustum the cone primitive already
//! builds. The sweeps that need *new* surfaces — free-form spines, skew
//! ruled walls, smoothed skinning through many sections — are recorded in
//! the deferred table, not approximated here.

use og_algo::{
    Built, edge_vertices, make_cone, make_cylinder, make_edge_between, make_face_with_pcurves,
    make_solid, make_torus, make_vertex, sew,
};
use og_core::{OgResult, Tolerances, og_bail};
use og_geom::Curve3d as _;
use og_geom::{CircleCurve, Curve, Line2d, LineCurve, PlaneSurface, SurfaceGeometry, TorusSurface};
use og_math::{Circle, Direction, Frame, Plane, Point, Point2, Torus, Vector};
use og_topo::{EdgeRepr, Filter, Model, Shape, ShapeType, explore};

/// Sweep a circular profile of `radius` along a spine edge.
///
/// A straight spine gives a cylinder, a full circular spine a torus, an arc
/// a torus segment. The history generates the solid from the spine.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the spine is
/// not straight or circular, the profile radius is not usable, or the tube
/// would swallow its own spine.
pub fn make_pipe(
    model: &mut Model,
    spine: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgResult<Built> {
    if !radius.is_finite() || radius <= tol.confusion() {
        og_bail!(Construction, "a pipe of radius {radius} holds nothing");
    }
    let (curve, range) = {
        let Some(data) = model.node(spine).and_then(|n| n.data().as_edge()) else {
            og_bail!(Construction, "a pipe runs along an edge");
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            og_bail!(Construction, "the spine has no curve");
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            og_bail!(Dangling, "curve is not in this model");
        };
        (geometry.clone(), *range)
    };
    let mut built = match &curve {
        Curve::Line(line) => {
            let start = curve.point_at(range.0, tol)?;
            let length = range.1 - range.0;
            let frame = Frame::about(start, line.axis().direction);
            make_cylinder(model, frame, radius, length, tol)?
        }
        Curve::Circle(c) => {
            let circle = c.circle();
            if radius >= circle.radius() - tol.confusion() {
                og_bail!(
                    Construction,
                    "a tube of radius {radius} swallows its spine of radius {}",
                    circle.radius()
                );
            }
            let closed = curve
                .point_at(range.0, tol)?
                .distance(curve.point_at(range.1, tol)?)
                <= tol.confusion();
            if closed {
                make_torus(model, circle.frame(), circle.radius(), radius, tol)?
            } else {
                pipe_segment(model, circle, range, radius, tol)?
            }
        }
        _ => og_bail!(
            Construction,
            "a pipe along a free-form spine needs the sweep-surface \
             machinery — see the deferred table"
        ),
    };
    built.history.generate(spine, built.shape.clone());
    Ok(built)
}

/// The torus segment: two half-tube patches and two meridian caps.
///
/// The tube circles are framed so their own parameter *is* the torus tube
/// angle, which makes every tube pcurve a vertical line in the chart and the
/// two patches the clean rectangles `v ∈ [0, π]` and `[π, 2π]`. The outer
/// equator is then the seam between the halves across the period — one edge,
/// two chart rows — which is exactly what [`og_algo::attach_seam`] exists to
/// say.
fn pipe_segment(
    model: &mut Model,
    spine: Circle,
    range: (f64, f64),
    radius: f64,
    tol: Tolerances,
) -> OgResult<Built> {
    let frame = spine.frame();
    let (x, y, z) = (frame.x().vector(), frame.y().vector(), frame.z().vector());
    let major = spine.radius();
    let radial = |u: f64| x * u.cos() + y * u.sin();
    let tangent = |u: f64| x * -u.sin() + y * u.cos();
    let tube_point = |u: f64, v: f64| {
        frame.origin() + radial(u) * radius.mul_add(v.cos(), major) + z * (radius * v.sin())
    };
    let pi = core::f64::consts::PI;
    let tau = core::f64::consts::TAU;

    let torus: SurfaceGeometry = TorusSurface::new(Torus::new(frame, major, radius, tol)?).into();
    let surface_id = model.geometry_mut().add_surface(torus);

    // Vertices at the tube's v = 0 and v = π points of each end.
    let ends = [range.0, range.1];
    let mut verts: Vec<Vec<Shape>> = Vec::new();
    for &u in &ends {
        verts.push(vec![
            make_vertex(model, tube_point(u, 0.0)).shape,
            make_vertex(model, tube_point(u, pi)).shape,
        ]);
    }

    // The tube circles at each end, split at v = 0 and v = π, framed so that
    // the circle's parameter equals the torus tube angle: `z` against the
    // spine tangent makes the frame's `y` the torus's own axis.
    let mut tube_arcs: Vec<Vec<Shape>> = Vec::new();
    for (k, &u) in ends.iter().enumerate() {
        let centre = frame.origin() + radial(u) * major;
        let circle = Circle::new(
            Frame::new(
                centre,
                Direction::new(-tangent(u), tol)?,
                Direction::new(radial(u), tol)?,
                tol,
            )?,
            radius,
            tol,
        )?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        let arcs = vec![
            make_edge_between(
                model,
                curve.clone(),
                (0.0, pi),
                &verts[k][0],
                &verts[k][1],
                tol,
            )?
            .shape,
            make_edge_between(model, curve, (pi, tau), &verts[k][1], &verts[k][0], tol)?.shape,
        ];
        // In the chart both arcs run up the column at this end's angle.
        let column = Line2d::over(
            og_math::Axis2::new(Point2::new(u, 0.0), og_math::Direction2::Y),
            0.0,
            tau,
        )?;
        og_algo::attach_pcurve(
            model,
            &arcs[0],
            column.into(),
            surface_id,
            og_topo::Location::identity(),
            (0.0, pi),
        )?;
        og_algo::attach_pcurve(
            model,
            &arcs[1],
            column.into(),
            surface_id,
            og_topo::Location::identity(),
            (pi, tau),
        )?;
        tube_arcs.push(arcs);
    }

    // The long edges: the parallels at v = 0 and v = π, parameterized by the
    // spine's own angle.
    let parallel = |model: &mut Model, v: f64, from: &Shape, to: &Shape| -> OgResult<Shape> {
        let height = radius * v.sin();
        let ring = radius.mul_add(v.cos(), major);
        let circle = Circle::new(
            Frame::new(frame.origin() + z * height, frame.z(), frame.x(), tol)?,
            ring,
            tol,
        )?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        Ok(make_edge_between(model, curve, range, from, to, tol)?.shape)
    };
    let row = |v: f64| -> OgResult<Line2d> {
        Line2d::over(
            og_math::Axis2::new(Point2::new(0.0, v), og_math::Direction2::X),
            range.0 - 1.0,
            range.1 + 1.0,
        )
    };
    let inner = parallel(model, pi, &verts[0][1], &verts[1][1])?;
    og_algo::attach_pcurve(
        model,
        &inner,
        row(pi)?.into(),
        surface_id,
        og_topo::Location::identity(),
        range,
    )?;
    // The outer equator bounds both halves across the period: v = 2π for its
    // forward use under the upper patch, v = 0 for its reversed use under the
    // lower — a seam, said as one.
    let outer = parallel(model, 0.0, &verts[0][0], &verts[1][0])?;
    og_algo::attach_seam(
        model,
        &outer,
        row(tau)?.into(),
        row(0.0)?.into(),
        surface_id,
        og_topo::Location::identity(),
        range,
    )?;

    // The two half-tube patches, on the one registered surface.
    let lower = {
        let wire = og_algo::make_wire(
            model,
            &[
                tube_arcs[0][0].clone(),
                inner.clone(),
                tube_arcs[1][0].reversed(),
                outer.reversed(),
            ],
            tol,
        )?
        .shape;
        og_algo::make_face_on(model, surface_id, std::slice::from_ref(&wire), tol)?.shape
    };
    let upper = {
        let wire = og_algo::make_wire(
            model,
            &[
                tube_arcs[0][1].clone(),
                outer.clone(),
                tube_arcs[1][1].reversed(),
                inner.reversed(),
            ],
            tol,
        )?
        .shape;
        og_algo::make_face_on(model, surface_id, std::slice::from_ref(&wire), tol)?.shape
    };

    // The meridian caps, their outward normals along the spine and away from
    // the material between the ends.
    let mut caps: Vec<Shape> = Vec::new();
    for (k, &u) in ends.iter().enumerate() {
        let outward = if k == 0 { -tangent(u) } else { tangent(u) };
        let centre = frame.origin() + radial(u) * major;
        let plane = Plane::through(centre, Direction::new(outward, tol)?);
        let reach = (major + radius) * 2.0;
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        caps.push(
            make_face_with_pcurves(
                model,
                surface,
                &[vec![tube_arcs[k][0].clone(), tube_arcs[k][1].clone()]],
                tol,
            )?
            .shape,
        );
    }

    let faces = [lower, upper, caps[0].clone(), caps[1].clone()];
    let sewn = sew(model, &faces, tol)?;
    if sewn.shells.len() != 1 || !og_algo::is_shell_closed(model, &sewn.shells[0])? {
        og_bail!(Construction, "the pipe segment did not close");
    }
    make_solid(model, std::slice::from_ref(&sewn.shells[0]))
}

/// Loft two parallel closed sections into a solid, ruled.
///
/// Two coaxial circles give the cylinder or the cone frustum; two polygons
/// with the same corner count give planar walls. The sections pair edge by
/// edge in traversal order.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the sections
/// are not both circles or both polygons of the same count, are not on
/// parallel planes, or a ruled wall would be skew.
pub fn make_loft(
    model: &mut Model,
    bottom: &Shape,
    top: &Shape,
    tol: Tolerances,
) -> OgResult<Built> {
    for wire in [bottom, top] {
        if model.kind_of(wire)? != ShapeType::Wire {
            og_bail!(Construction, "a loft runs between wires");
        }
        if !og_algo::is_wire_closed(model, wire, tol)? {
            og_bail!(Construction, "a loft section must be closed");
        }
    }
    let circle_of = |model: &Model, wire: &Shape| -> OgResult<Option<Circle>> {
        let edges = explore(model, wire, Filter::OfType(ShapeType::Edge))?;
        if edges.len() != 1 {
            return Ok(None);
        }
        let Some(data) = model.node(&edges[0]).and_then(|n| n.data().as_edge()) else {
            return Ok(None);
        };
        let Some(EdgeRepr::Curve3d { curve, .. }) = data.curve3d() else {
            return Ok(None);
        };
        match model.geometry().curve(*curve) {
            Some(Curve::Circle(c)) => Ok(Some(c.circle())),
            _ => Ok(None),
        }
    };

    if let (Some(lower), Some(upper)) = (circle_of(model, bottom)?, circle_of(model, top)?) {
        // Coaxial circles: the revolved primitives already build these.
        let axis = lower.frame().z().vector();
        let rise = upper.centre() - lower.centre();
        let height = rise.dot(axis);
        if rise.cross(axis).magnitude() > tol.confusion() * 10.0 || height.abs() <= tol.confusion()
        {
            og_bail!(
                Construction,
                "lofted circles must be coaxial on parallel planes; the \
                 oblique loft needs the sweep machinery — see the deferred \
                 table"
            );
        }
        let frame = if height > 0.0 {
            Frame::new(lower.centre(), lower.frame().z(), lower.frame().x(), tol)?
        } else {
            Frame::new(
                lower.centre(),
                lower.frame().z().reversed(),
                lower.frame().x(),
                tol,
            )?
        };
        let mut built = if (lower.radius() - upper.radius()).abs() <= tol.confusion() {
            make_cylinder(model, frame, lower.radius(), height.abs(), tol)?
        } else {
            make_cone(
                model,
                frame,
                lower.radius(),
                upper.radius(),
                height.abs(),
                tol,
            )?
        };
        built.history.generate(bottom, built.shape.clone());
        built.history.generate(top, built.shape.clone());
        return Ok(built);
    }

    // Polygonal sections: matched corners, planar walls.
    let corners_of = |model: &Model, wire: &Shape| -> OgResult<Vec<Point>> {
        let mut out = Vec::new();
        for edge in explore(model, wire, Filter::OfType(ShapeType::Edge))? {
            let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                og_bail!(Construction, "a section edge holds no data");
            };
            let Some(EdgeRepr::Curve3d { curve, .. }) = data.curve3d() else {
                og_bail!(Construction, "a section edge has no curve");
            };
            let Some(Curve::Line(_)) = model.geometry().curve(*curve) else {
                og_bail!(
                    Construction,
                    "a mixed or curved section needs the skinning machinery — \
                     see the deferred table"
                );
            };
            let Some((sv, _)) = edge_vertices(model, &edge)? else {
                og_bail!(Construction, "a section edge has no vertices");
            };
            let Some(data) = model.node(&sv).and_then(|n| n.data().as_vertex()) else {
                og_bail!(Construction, "a section vertex holds no point");
            };
            out.push(sv.transform(model.datums())?.apply(data.point));
        }
        Ok(out)
    };
    let low = corners_of(model, bottom)?;
    let high = corners_of(model, top)?;
    if low.len() != high.len() {
        og_bail!(
            Construction,
            "lofted sections must have the same corner count, found {} and {}",
            low.len(),
            high.len()
        );
    }
    let n = low.len();
    let centroid = {
        let mut c = Vector::new(0.0, 0.0, 0.0);
        for p in low.iter().chain(high.iter()) {
            c += p.to_vector();
        }
        #[allow(clippy::cast_precision_loss)]
        let count = 2.0 * n as f64;
        Point::from_vector(c / count)
    };

    // Shared vertices and edges, then walls and caps referencing them.
    let vl: Vec<Shape> = low.iter().map(|p| make_vertex(model, *p).shape).collect();
    let vh: Vec<Shape> = high.iter().map(|p| make_vertex(model, *p).shape).collect();
    let seg = |model: &mut Model, a: (&Shape, Point), b: (&Shape, Point)| -> OgResult<Shape> {
        let line = LineCurve::segment(a.1, b.1, tol)?;
        let curve = Curve::Line(line);
        let domain = curve.domain();
        Ok(make_edge_between(model, curve, domain, a.0, b.0, tol)?.shape)
    };
    let mut low_edges = Vec::with_capacity(n);
    let mut high_edges = Vec::with_capacity(n);
    let mut rails = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        low_edges.push(seg(model, (&vl[i], low[i]), (&vl[j], low[j]))?);
        high_edges.push(seg(model, (&vh[i], high[i]), (&vh[j], high[j]))?);
        rails.push(seg(model, (&vl[i], low[i]), (&vh[i], high[i]))?);
    }

    let planar = |model: &mut Model, corners: &[Point], edges: Vec<Shape>| -> OgResult<Shape> {
        let normal = {
            let mut n = (corners[1] - corners[0]).cross(corners[2] - corners[0]);
            let m = n.magnitude();
            if m <= tol.confusion() {
                og_bail!(Construction, "a loft wall is degenerate");
            }
            n /= m;
            if n.dot(corners[0] - centroid) < 0.0 {
                -n
            } else {
                n
            }
        };
        for p in corners {
            if Plane::through(corners[0], Direction::new(normal, tol)?).distance_to(*p)
                > tol.confusion() * 10.0
            {
                og_bail!(
                    Construction,
                    "a skew ruled wall is not a plane; it needs the \
                     extrusion-surface loft — see the deferred table"
                );
            }
        }
        let plane = Plane::through(corners[0], Direction::new(normal, tol)?);
        let mut reach = 1.0_f64;
        for p in corners {
            reach = reach.max(p.distance(corners[0]) * 2.0);
        }
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        Ok(make_face_with_pcurves(model, surface, &[edges], tol)?.shape)
    };

    let mut faces: Vec<Shape> = Vec::with_capacity(n + 2);
    for i in 0..n {
        let j = (i + 1) % n;
        faces.push(planar(
            model,
            &[low[i], low[j], high[j], high[i]],
            vec![
                low_edges[i].clone(),
                rails[j].clone(),
                high_edges[i].reversed(),
                rails[i].reversed(),
            ],
        )?);
    }
    faces.push(planar(model, &low, low_edges.clone())?);
    faces.push(planar(model, &high, high_edges.clone())?);

    let sewn = sew(model, &faces, tol)?;
    if sewn.shells.len() != 1 || !og_algo::is_shell_closed(model, &sewn.shells[0])? {
        og_bail!(Construction, "the loft did not close");
    }
    let mut built = make_solid(model, std::slice::from_ref(&sewn.shells[0]))?;
    built.history.generate(bottom, built.shape.clone());
    built.history.generate(top, built.shape.clone());
    Ok(built)
}
