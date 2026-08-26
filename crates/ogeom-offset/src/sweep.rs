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
//! docs/PARITY.md (offset.sweeps), not approximated here.

use ogeom_algo::{
    Built, History, edge_vertices, make_cone, make_cylinder, make_edge, make_edge_between,
    make_face_with_pcurves, make_solid, make_torus, make_vertex, sew,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::Transformable as _;
use ogeom_geom::{
    CircleCurve, Curve, Line2d, LineCurve, PlaneSurface, SurfaceGeometry, TorusSurface,
};
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Point2, Torus, Transform, Vector};
use ogeom_topo::{EdgeData, EdgeRepr, Filter, Model, Shape, ShapeType, VertexData, explore};

/// Sweep a circular profile of `radius` along a spine edge.
///
/// A straight spine gives a cylinder, a full circular spine a torus, an arc
/// a torus segment. The history generates the solid from the spine.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the spine is
/// not straight or circular, the profile radius is not usable, or the tube
/// would swallow its own spine.
pub fn make_pipe(
    model: &mut Model,
    spine: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(Construction, "a pipe of radius {radius} holds nothing");
    }
    let (curve, range) = {
        let Some(data) = model.node(spine).and_then(|n| n.data().as_edge()) else {
            ogeom_bail!(Construction, "a pipe runs along an edge");
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            ogeom_bail!(Construction, "the spine has no curve");
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            ogeom_bail!(Dangling, "curve is not in this model");
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
                ogeom_bail!(
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
        _ => ogeom_bail!(
            Construction,
            "a pipe along a free-form spine needs the sweep-surface \
             machinery — docs/PARITY.md, offset.sweeps"
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
/// two chart rows — which is exactly what [`ogeom_algo::attach_seam`] exists to
/// say.
fn pipe_segment(
    model: &mut Model,
    spine: Circle,
    range: (f64, f64),
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
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
            ogeom_math::Axis2::new(Point2::new(u, 0.0), ogeom_math::Direction2::Y),
            0.0,
            tau,
        )?;
        ogeom_algo::attach_pcurve(
            model,
            &arcs[0],
            column.into(),
            surface_id,
            ogeom_topo::Location::identity(),
            (0.0, pi),
        )?;
        ogeom_algo::attach_pcurve(
            model,
            &arcs[1],
            column.into(),
            surface_id,
            ogeom_topo::Location::identity(),
            (pi, tau),
        )?;
        tube_arcs.push(arcs);
    }

    // The long edges: the parallels at v = 0 and v = π, parameterized by the
    // spine's own angle.
    let parallel = |model: &mut Model, v: f64, from: &Shape, to: &Shape| -> OgeomResult<Shape> {
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
    let row = |v: f64| -> OgeomResult<Line2d> {
        Line2d::over(
            ogeom_math::Axis2::new(Point2::new(0.0, v), ogeom_math::Direction2::X),
            range.0 - 1.0,
            range.1 + 1.0,
        )
    };
    let inner = parallel(model, pi, &verts[0][1], &verts[1][1])?;
    ogeom_algo::attach_pcurve(
        model,
        &inner,
        row(pi)?.into(),
        surface_id,
        ogeom_topo::Location::identity(),
        range,
    )?;
    // The outer equator bounds both halves across the period: v = 2π for its
    // forward use under the upper patch, v = 0 for its reversed use under the
    // lower — a seam, said as one.
    let outer = parallel(model, 0.0, &verts[0][0], &verts[1][0])?;
    ogeom_algo::attach_seam(
        model,
        &outer,
        row(tau)?.into(),
        row(0.0)?.into(),
        surface_id,
        ogeom_topo::Location::identity(),
        range,
    )?;

    // The two half-tube patches, on the one registered surface.
    let lower = {
        let wire = ogeom_algo::make_wire(
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
        ogeom_algo::make_face_on(model, surface_id, std::slice::from_ref(&wire), tol)?.shape
    };
    let upper = {
        let wire = ogeom_algo::make_wire(
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
        ogeom_algo::make_face_on(model, surface_id, std::slice::from_ref(&wire), tol)?.shape
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
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the pipe segment did not close");
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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the sections
/// are not both circles or both polygons of the same count, are not on
/// parallel planes, or a ruled wall would be skew.
pub fn make_loft(
    model: &mut Model,
    bottom: &Shape,
    top: &Shape,
    tol: Tolerances,
) -> OgeomResult<Built> {
    // A vertex for either section is the loft to a point: a cone over a
    // circle, an exact pyramid over a polygon.
    match (model.kind_of(bottom)?, model.kind_of(top)?) {
        (ShapeType::Wire, ShapeType::Vertex) => {
            return loft_to_point(model, bottom, top, tol);
        }
        (ShapeType::Vertex, ShapeType::Wire) => {
            let mut built = loft_to_point(model, top, bottom, tol)?;
            // The apex was named first: the same solid, the history the same.
            built.history.generate(bottom, built.shape.clone());
            return Ok(built);
        }
        _ => {}
    }
    for wire in [bottom, top] {
        if model.kind_of(wire)? != ShapeType::Wire {
            ogeom_bail!(Construction, "a loft runs between wires");
        }
        if !ogeom_algo::is_wire_closed(model, wire, tol)? {
            ogeom_bail!(Construction, "a loft section must be closed");
        }
    }
    let circle_of = |model: &Model, wire: &Shape| -> OgeomResult<Option<Circle>> {
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
            ogeom_bail!(
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
    let corners_of = |model: &Model, wire: &Shape| -> OgeomResult<Vec<Point>> {
        let mut out = Vec::new();
        for edge in explore(model, wire, Filter::OfType(ShapeType::Edge))? {
            let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                ogeom_bail!(Construction, "a section edge holds no data");
            };
            let Some(EdgeRepr::Curve3d { curve, .. }) = data.curve3d() else {
                ogeom_bail!(Construction, "a section edge has no curve");
            };
            let Some(Curve::Line(_)) = model.geometry().curve(*curve) else {
                ogeom_bail!(
                    Construction,
                    "a mixed or curved section needs the skinning machinery — \
                     docs/PARITY.md, offset.sweeps"
                );
            };
            let Some((sv, _)) = edge_vertices(model, &edge)? else {
                ogeom_bail!(Construction, "a section edge has no vertices");
            };
            let Some(data) = model.node(&sv).and_then(|n| n.data().as_vertex()) else {
                ogeom_bail!(Construction, "a section vertex holds no point");
            };
            out.push(sv.transform(model.datums())?.apply(data.point));
        }
        Ok(out)
    };
    let low = corners_of(model, bottom)?;
    let high = corners_of(model, top)?;
    if low.len() != high.len() {
        ogeom_bail!(
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
    let seg = |model: &mut Model, a: (&Shape, Point), b: (&Shape, Point)| -> OgeomResult<Shape> {
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

    let planar = |model: &mut Model, corners: &[Point], edges: Vec<Shape>| -> OgeomResult<Shape> {
        let normal = {
            let mut n = (corners[1] - corners[0]).cross(corners[2] - corners[0]);
            let m = n.magnitude();
            if m <= tol.confusion() {
                ogeom_bail!(Construction, "a loft wall is degenerate");
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
                ogeom_bail!(
                    Construction,
                    "a skew ruled wall is not a plane; loft skew sections \
                     through make_loft_skinned, aligned with hints — \
                     docs/PARITY.md, offset.loft"
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
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the loft did not close");
    }
    let mut built = make_solid(model, std::slice::from_ref(&sewn.shells[0]))?;
    built.history.generate(bottom, built.shape.clone());
    built.history.generate(top, built.shape.clone());
    Ok(built)
}

/// A skinned wall and the pieces a caller needs to close it: the rings at
/// both ends, their exact border curves off the control net, and the chart's
/// `u` window the ring pcurves span.
struct SkinnedWall {
    face: Shape,
    ring0: Shape,
    ring1: Shape,
    curve0: ogeom_geom::Curve,
    curve1: ogeom_geom::Curve,
    u_dom: (f64, f64),
}

/// The wall of a skin over a grid of section samples, closed the way round.
///
/// The wall is [`ogeom_geom::fit::fit_surface_grid`]'s surface with each row's
/// first sample repeated at its end: the row fits pin their ends, so the two
/// border control columns are *equal* and the seam closes exactly, not
/// within tolerance. The border iso-curves come straight off the control
/// net — the v-borders are the fitted sections, planar whenever the
/// sections are, which is what lets the caps be planes — and every pcurve
/// is an iso line in the fitted chart, same-parameter by construction.
fn skinned_wall(
    model: &mut Model,
    rows: &[Vec<Point>],
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<SkinnedWall> {
    use ogeom_geom::Surface as _;
    let mut closed_rows: Vec<Vec<Point>> = Vec::with_capacity(rows.len());
    for row in rows {
        let mut r = row.clone();
        r.push(row[0]);
        closed_rows.push(r);
    }
    let fitted = ogeom_geom::fit::fit_surface_grid(&closed_rows, 3, tolerance, tol)?;
    if !fitted.met {
        ogeom_bail!(
            NotDone,
            "the skin reached {} against a target of {tolerance}",
            fitted.error
        );
    }
    let surface = fitted.curve;
    let (u_knots, v_knots) = (surface.u_knots().clone(), surface.v_knots().clone());
    let (k, l, net) = {
        let grid = surface.grid();
        let net: Vec<Point> = grid.points().iter().map(|w| (*w).point()).collect();
        (grid.u_count(), grid.v_count(), net)
    };
    let point_at = |i: usize, j: usize| -> Point { net[i * l + j] };
    let (u_dom, v_dom) = surface.domain();

    // Border curves straight off the net: v-borders are the end sections,
    // the u-border is the seam.
    let border_v = |j: usize| -> OgeomResult<ogeom_geom::Curve> {
        let control: Vec<Point> = (0..k).map(|i| point_at(i, j)).collect();
        Ok(ogeom_geom::Curve::BSpline(ogeom_geom::BSplineCurve::new(
            u_knots.clone(),
            control,
            tol,
        )?))
    };
    let seam_curve = {
        let control: Vec<Point> = (0..l).map(|j| point_at(0, j)).collect();
        ogeom_geom::Curve::BSpline(ogeom_geom::BSplineCurve::new(
            v_knots.clone(),
            control,
            tol,
        )?)
    };

    let surface_geo: SurfaceGeometry = surface.into();
    let surface_id = model.geometry_mut().add_surface(surface_geo.clone());

    let ring0 = make_edge(model, border_v(0)?, u_dom, tol)?.shape;
    let ring1 = make_edge(model, border_v(l - 1)?, u_dom, tol)?.shape;
    let anchor0 = ogeom_algo::edge_vertices(model, &ring0)?
        .map(|(a, _)| a)
        .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "a skinned ring has no vertex"))?;
    let anchor1 = ogeom_algo::edge_vertices(model, &ring1)?
        .map(|(a, _)| a)
        .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "a skinned ring has no vertex"))?;
    let seam = make_edge_between(model, seam_curve, v_dom, &anchor0, &anchor1, tol)?.shape;

    // Pcurves: rows for the rings, both columns for the seam.
    let row_line = |v: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
        Ok(Line2d::over(
            ogeom_math::Axis2::new(Point2::new(0.0, v), ogeom_math::Direction2::X),
            u_dom.0 - 1.0,
            u_dom.1 + 1.0,
        )?
        .into())
    };
    let column_line = |u: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
        Ok(Line2d::over(
            ogeom_math::Axis2::new(Point2::new(u, 0.0), ogeom_math::Direction2::Y),
            v_dom.0 - 1.0,
            v_dom.1 + 1.0,
        )?
        .into())
    };
    ogeom_algo::attach_pcurve(
        model,
        &ring0,
        row_line(v_dom.0)?,
        surface_id,
        ogeom_topo::Location::identity(),
        u_dom,
    )?;
    ogeom_algo::attach_pcurve(
        model,
        &ring1,
        row_line(v_dom.1)?,
        surface_id,
        ogeom_topo::Location::identity(),
        u_dom,
    )?;
    ogeom_algo::attach_seam(
        model,
        &seam,
        column_line(u_dom.0)?,
        column_line(u_dom.1)?,
        surface_id,
        ogeom_topo::Location::identity(),
        v_dom,
    )?;

    let wall = {
        let wire = ogeom_algo::make_wire(
            model,
            &[
                ring0.clone(),
                seam.clone(),
                ring1.reversed(),
                seam.reversed(),
            ],
            tol,
        )?
        .shape;
        let face =
            ogeom_algo::make_face_on(model, surface_id, std::slice::from_ref(&wire), tol)?.shape;
        // Outward by measurement at the middle of the skin.
        let mid_u = f64::midpoint(u_dom.0, u_dom.1);
        let mid_v = f64::midpoint(v_dom.0, v_dom.1);
        let s_mid = surface_geo.point_at(mid_u, mid_v, tol)?;
        let (du, dv) = surface_geo.d1_at(mid_u, mid_v, tol)?;
        let centroid = {
            let mut c = Vector::new(0.0, 0.0, 0.0);
            let mut n = 0.0;
            for row in rows {
                for p in row {
                    c += p.to_vector();
                    n += 1.0;
                }
            }
            Point::from_vector(c / n)
        };
        if du.cross(dv).dot(s_mid - centroid) >= 0.0 {
            face
        } else {
            face.reversed()
        }
    };
    Ok(SkinnedWall {
        face: wall,
        ring0,
        ring1,
        curve0: border_v(0)?,
        curve1: border_v(l - 1)?,
        u_dom,
    })
}

/// A strip closed the *long* way: open across its own width, a smooth loop
/// along the sweep — one face of a faceted ring, [`skinned_wall`]'s
/// construction with the chart's roles swapped and the loop made C1 by
/// [`ogeom_geom::fit::fit_surface_grid_closed_v`]. The rails are the two
/// closed border loops; the seam is one station's column, used twice.
fn skinned_ring_strip(
    model: &mut Model,
    rows: &[Vec<Point>],
    outward_hint: Point,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    use ogeom_geom::Surface as _;
    // The loop: first row repeated at the end, as the closed fit demands.
    let mut looped: Vec<Vec<Point>> = rows.to_vec();
    looped.push(rows[0].clone());
    let fitted = ogeom_geom::fit::fit_surface_grid_closed_v(&looped, 3, tolerance, tol)?;
    if !fitted.met {
        ogeom_bail!(
            NotDone,
            "the ring strip reached {} against a target of {tolerance}",
            fitted.error
        );
    }
    let surface = fitted.curve;
    let (u_knots, v_knots) = (surface.u_knots().clone(), surface.v_knots().clone());
    let (k, l, net) = {
        let grid = surface.grid();
        let net: Vec<Point> = grid.points().iter().map(|w| (*w).point()).collect();
        (grid.u_count(), grid.v_count(), net)
    };
    let point_at = |i: usize, j: usize| -> Point { net[i * l + j] };
    let (u_dom, v_dom) = surface.domain();

    // The chart's roles, straight: `u` runs across the strip (open), `v`
    // around the loop (closed). The rails are v-curves — the closed border
    // loops at the two u-borders — and the seam is the u-row at the loop's
    // join, bounding the chart twice as every seam does.
    let rail_curve = |i: usize| -> OgeomResult<ogeom_geom::Curve> {
        let control: Vec<Point> = (0..l).map(|j| point_at(i, j)).collect();
        Ok(ogeom_geom::Curve::BSpline(ogeom_geom::BSplineCurve::new(
            v_knots.clone(),
            control,
            tol,
        )?))
    };
    let seam_curve = {
        let control: Vec<Point> = (0..k).map(|i| point_at(i, 0)).collect();
        ogeom_geom::Curve::BSpline(ogeom_geom::BSplineCurve::new(
            u_knots.clone(),
            control,
            tol,
        )?)
    };
    let surface_geo: SurfaceGeometry = surface.into();
    let surface_id = model.geometry_mut().add_surface(surface_geo.clone());

    let rail0 = make_edge(model, rail_curve(0)?, v_dom, tol)?.shape;
    let rail1 = make_edge(model, rail_curve(k - 1)?, v_dom, tol)?.shape;
    // Neighbouring strips fit the shared corner loop independently; the
    // rails agree only within the fits' own honesty, and the sew can only
    // join what the tolerances admit. Every border widens to the error the
    // fit reported — skinned_strip's own discipline.
    let slack = fitted.error + tol.confusion();
    for edge in [&rail0, &rail1] {
        model.widen(edge, ogeom_core::Tolerance::new(slack)?)?;
    }
    let anchor0 = ogeom_algo::edge_vertices(model, &rail0)?
        .map(|(a, _)| a)
        .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "a strip rail has no vertex"))?;
    let anchor1 = ogeom_algo::edge_vertices(model, &rail1)?
        .map(|(a, _)| a)
        .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "a strip rail has no vertex"))?;
    let seam = make_edge_between(model, seam_curve, u_dom, &anchor0, &anchor1, tol)?.shape;

    let row_line = |v: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
        Ok(Line2d::over(
            ogeom_math::Axis2::new(Point2::new(0.0, v), ogeom_math::Direction2::X),
            u_dom.0 - 1.0,
            u_dom.1 + 1.0,
        )?
        .into())
    };
    let column_line = |u: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
        Ok(Line2d::over(
            ogeom_math::Axis2::new(Point2::new(u, 0.0), ogeom_math::Direction2::Y),
            v_dom.0 - 1.0,
            v_dom.1 + 1.0,
        )?
        .into())
    };
    ogeom_algo::attach_pcurve(
        model,
        &rail0,
        column_line(u_dom.0)?,
        surface_id,
        ogeom_topo::Location::identity(),
        v_dom,
    )?;
    ogeom_algo::attach_pcurve(
        model,
        &rail1,
        column_line(u_dom.1)?,
        surface_id,
        ogeom_topo::Location::identity(),
        v_dom,
    )?;
    ogeom_algo::attach_seam(
        model,
        &seam,
        row_line(v_dom.0)?,
        row_line(v_dom.1)?,
        surface_id,
        ogeom_topo::Location::identity(),
        u_dom,
    )?;
    let wire = ogeom_algo::make_wire(
        model,
        &[
            rail0.clone(),
            seam.clone(),
            rail1.reversed(),
            seam.reversed(),
        ],
        tol,
    )?
    .shape;
    let face = ogeom_algo::make_face_on(model, surface_id, std::slice::from_ref(&wire), tol)?.shape;
    let mid_u = f64::midpoint(u_dom.0, u_dom.1);
    let mid_v = f64::midpoint(v_dom.0, v_dom.1);
    let s_mid = surface_geo.point_at(mid_u, mid_v, tol)?;
    let (du, dv) = surface_geo.d1_at(mid_u, mid_v, tol)?;
    Ok(if du.cross(dv).dot(s_mid - outward_hint) >= 0.0 {
        face
    } else {
        face.reversed()
    })
}

/// A solid skinned over a grid of section samples: [`skinned_wall`] with a
/// planar cap over each end ring.
fn skinned_solid(
    model: &mut Model,
    rows: &[Vec<Point>],
    cap_outward: (Vector, Vector),
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let wall = skinned_wall(model, rows, tolerance, tol)?;
    let u_dom = wall.u_dom;

    let cap = |model: &mut Model,
               ring: &Shape,
               curve: ogeom_geom::Curve,
               outward: Vector|
     -> OgeomResult<Shape> {
        let at = curve.point_at(u_dom.0, tol)?;
        let plane = Plane::through(at, Direction::new(outward, tol)?);
        let mut reach = 1.0_f64;
        for t in 0..8 {
            let p = curve.point_at(u_dom.0 + (u_dom.1 - u_dom.0) * f64::from(t) / 8.0, tol)?;
            reach = reach.max(p.distance(at) * 2.0);
        }
        let cap_surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let wire = ogeom_algo::make_wire(model, std::slice::from_ref(ring), tol)?.shape;
        let face =
            ogeom_algo::make_face(model, cap_surface.clone(), std::slice::from_ref(&wire), tol)?
                .shape;
        let id = {
            let Some(node) = model.node(&face) else {
                ogeom_bail!(Dangling, "the cap just built is not in this model");
            };
            let ogeom_topo::NodeData::Face(data) = node.data() else {
                ogeom_bail!(Construction, "the cap holds no face data");
            };
            data.surface
        };
        let Some(pcurve) = ogeom_intersect::exact_pcurve_of(&curve, &cap_surface, tol) else {
            ogeom_bail!(Construction, "a cap edge has no closed-form pcurve");
        };
        ogeom_algo::attach_pcurve(
            model,
            ring,
            pcurve,
            id,
            ogeom_topo::Location::identity(),
            u_dom,
        )?;
        Ok(face)
    };
    let cap0 = cap(model, &wall.ring0, wall.curve0.clone(), cap_outward.0)?;
    let cap1 = cap(model, &wall.ring1, wall.curve1.clone(), cap_outward.1)?;

    let faces = [wall.face, cap0, cap1];
    let sewn = sew(model, &faces, tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the skinned solid did not close");
    }
    make_solid(model, std::slice::from_ref(&sewn.shells[0]))
}

/// A solid skinned down to a point: [`skinned_wall`]'s construction with the
/// top ring replaced by the apex — a degenerate edge on one vertex, bounding
/// the chart's whole top row the way a cone's apex bounds a countersink.
/// One cap, at the open end; the apex end closes by construction.
fn skinned_solid_to_apex(
    model: &mut Model,
    rows: &[Vec<Point>],
    cap_outward: Vector,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    use ogeom_geom::Surface as _;
    let mut closed_rows: Vec<Vec<Point>> = Vec::with_capacity(rows.len());
    for row in rows {
        let mut r = row.clone();
        r.push(row[0]);
        closed_rows.push(r);
    }
    let fitted = ogeom_geom::fit::fit_surface_grid(&closed_rows, 3, tolerance, tol)?;
    if !fitted.met {
        ogeom_bail!(
            NotDone,
            "the skin reached {} against a target of {tolerance}",
            fitted.error
        );
    }
    let surface = fitted.curve;
    let (u_knots, v_knots) = (surface.u_knots().clone(), surface.v_knots().clone());
    let (k, l, net) = {
        let grid = surface.grid();
        let net: Vec<Point> = grid.points().iter().map(|w| (*w).point()).collect();
        (grid.u_count(), grid.v_count(), net)
    };
    let point_at = |i: usize, j: usize| -> Point { net[i * l + j] };
    let (u_dom, v_dom) = surface.domain();
    let apex = rows[rows.len() - 1][0];

    let ring_curve = {
        let control: Vec<Point> = (0..k).map(|i| point_at(i, 0)).collect();
        ogeom_geom::Curve::BSpline(ogeom_geom::BSplineCurve::new(
            u_knots.clone(),
            control,
            tol,
        )?)
    };
    let seam_curve = {
        let control: Vec<Point> = (0..l).map(|j| point_at(0, j)).collect();
        ogeom_geom::Curve::BSpline(ogeom_geom::BSplineCurve::new(
            v_knots.clone(),
            control,
            tol,
        )?)
    };
    let surface_geo: SurfaceGeometry = surface.into();
    let surface_id = model.geometry_mut().add_surface(surface_geo.clone());

    let ring0 = make_edge(model, ring_curve.clone(), u_dom, tol)?.shape;
    let anchor0 = ogeom_algo::edge_vertices(model, &ring0)?
        .map(|(a, _)| a)
        .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "a skinned ring has no vertex"))?;
    let apex_vertex = model.add_vertex(VertexData::new(apex));
    let apex_edge = {
        let mut data = EdgeData::new();
        data.degenerate = true;
        model.add_edge(data, &[apex_vertex.clone(), apex_vertex.clone()])?
    };
    let seam = make_edge_between(model, seam_curve, v_dom, &anchor0, &apex_vertex, tol)?.shape;

    let row_line = |v: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
        Ok(Line2d::over(
            ogeom_math::Axis2::new(Point2::new(0.0, v), ogeom_math::Direction2::X),
            u_dom.0 - 1.0,
            u_dom.1 + 1.0,
        )?
        .into())
    };
    let column_line = |u: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
        Ok(Line2d::over(
            ogeom_math::Axis2::new(Point2::new(u, 0.0), ogeom_math::Direction2::Y),
            v_dom.0 - 1.0,
            v_dom.1 + 1.0,
        )?
        .into())
    };
    ogeom_algo::attach_pcurve(
        model,
        &ring0,
        row_line(v_dom.0)?,
        surface_id,
        ogeom_topo::Location::identity(),
        u_dom,
    )?;
    // The apex bounds the chart's whole top row while covering no distance:
    // the degenerate edge carries the row's pcurve, exactly as a cone's apex
    // does after the reader synthesises it.
    ogeom_algo::attach_pcurve(
        model,
        &apex_edge,
        row_line(v_dom.1)?,
        surface_id,
        ogeom_topo::Location::identity(),
        u_dom,
    )?;
    ogeom_algo::attach_seam(
        model,
        &seam,
        column_line(u_dom.0)?,
        column_line(u_dom.1)?,
        surface_id,
        ogeom_topo::Location::identity(),
        v_dom,
    )?;

    let wall = {
        let wire = ogeom_algo::make_wire(
            model,
            &[
                ring0.clone(),
                seam.clone(),
                apex_edge.reversed(),
                seam.reversed(),
            ],
            tol,
        )?
        .shape;
        let face =
            ogeom_algo::make_face_on(model, surface_id, std::slice::from_ref(&wire), tol)?.shape;
        let mid_u = f64::midpoint(u_dom.0, u_dom.1);
        let mid_v = f64::midpoint(v_dom.0, v_dom.1);
        let s_mid = surface_geo.point_at(mid_u, mid_v, tol)?;
        let (du, dv) = surface_geo.d1_at(mid_u, mid_v, tol)?;
        let centroid = {
            let mut c = Vector::new(0.0, 0.0, 0.0);
            let mut n = 0.0;
            for row in rows {
                for p in row {
                    c += p.to_vector();
                    n += 1.0;
                }
            }
            Point::from_vector(c / n)
        };
        if du.cross(dv).dot(s_mid - centroid) >= 0.0 {
            face
        } else {
            face.reversed()
        }
    };

    // One cap, on the open end; the machinery is skinned_solid's, inlined
    // for the single ring.
    let cap = {
        let at = ring_curve.point_at(u_dom.0, tol)?;
        let plane = Plane::through(at, Direction::new(cap_outward, tol)?);
        let mut reach = 1.0_f64;
        for t in 0..8 {
            let p = ring_curve.point_at(u_dom.0 + (u_dom.1 - u_dom.0) * f64::from(t) / 8.0, tol)?;
            reach = reach.max(p.distance(at) * 2.0);
        }
        let cap_surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let wire = ogeom_algo::make_wire(model, std::slice::from_ref(&ring0), tol)?.shape;
        let face =
            ogeom_algo::make_face(model, cap_surface.clone(), std::slice::from_ref(&wire), tol)?
                .shape;
        let id = {
            let Some(node) = model.node(&face) else {
                ogeom_bail!(Dangling, "the cap just built is not in this model");
            };
            let ogeom_topo::NodeData::Face(data) = node.data() else {
                ogeom_bail!(Construction, "the cap holds no face data");
            };
            data.surface
        };
        let Some(pcurve) = ogeom_intersect::exact_pcurve_of(&ring_curve, &cap_surface, tol) else {
            ogeom_bail!(Construction, "a cap edge has no closed-form pcurve");
        };
        ogeom_algo::attach_pcurve(
            model,
            &ring0,
            pcurve,
            id,
            ogeom_topo::Location::identity(),
            u_dom,
        )?;
        face
    };

    let faces = [wall, cap];
    let sewn = sew(model, &faces, tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the skinned apex solid did not close");
    }
    make_solid(model, std::slice::from_ref(&sewn.shells[0]))
}

/// A solid skinned over a grid of sections that loops back on itself: the
/// wall is one face closed in both chart directions, no caps at all.
///
/// The `u` seam closes the way every skin's does — pinned row ends — and
/// the `v` loop closes through [`ogeom_geom::fit::fit_surface_grid_closed_v`],
/// C1 across the join. All four boundary traversals are two seam edges used
/// twice, anchored at one shared vertex, exactly as a torus bounds itself.
fn closed_skinned_solid(
    model: &mut Model,
    rows: &[Vec<Point>],
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let shell = closed_skinned_shell(model, rows, tolerance, tol)?;
    make_solid(model, std::slice::from_ref(&shell))
}

/// The closed skin as a shell, for callers assembling solids with voids —
/// a holed profile's ring is one outer shell and one per tunnel.
fn closed_skinned_shell(
    model: &mut Model,
    rows: &[Vec<Point>],
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    use ogeom_geom::Surface as _;
    let mut closed_rows: Vec<Vec<Point>> = Vec::with_capacity(rows.len() + 1);
    for row in rows {
        let mut r = row.clone();
        r.push(row[0]);
        closed_rows.push(r);
    }
    closed_rows.push(closed_rows[0].clone());
    let fitted = ogeom_geom::fit::fit_surface_grid_closed_v(&closed_rows, 3, tolerance, tol)?;
    if !fitted.met {
        ogeom_bail!(
            NotDone,
            "the closed skin reached {} against a target of {tolerance}",
            fitted.error
        );
    }
    let surface = fitted.curve;
    let (u_knots, v_knots) = (surface.u_knots().clone(), surface.v_knots().clone());
    let (k, l, net) = {
        let grid = surface.grid();
        let net: Vec<Point> = grid.points().iter().map(|w| (*w).point()).collect();
        (grid.u_count(), grid.v_count(), net)
    };
    let point_at = |i: usize, j: usize| -> Point { net[i * l + j] };
    let (u_dom, v_dom) = surface.domain();

    // Both seams straight off the net: the u-run at v's join, and the v-run
    // at u's.
    let along_u = {
        let control: Vec<Point> = (0..k).map(|i| point_at(i, 0)).collect();
        ogeom_geom::Curve::BSpline(ogeom_geom::BSplineCurve::new(u_knots, control, tol)?)
    };
    let along_v = {
        let control: Vec<Point> = (0..l).map(|j| point_at(0, j)).collect();
        ogeom_geom::Curve::BSpline(ogeom_geom::BSplineCurve::new(v_knots, control, tol)?)
    };
    let surface_geo: SurfaceGeometry = surface.into();
    let surface_id = model.geometry_mut().add_surface(surface_geo.clone());

    let u_edge = make_edge(model, along_u, u_dom, tol)?.shape;
    let Some((corner, _)) = ogeom_algo::edge_vertices(model, &u_edge)? else {
        ogeom_bail!(Construction, "the closed skin's seam has no vertex");
    };
    let v_edge = make_edge_between(model, along_v, v_dom, &corner, &corner, tol)?.shape;

    let row_line = |v: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
        Ok(Line2d::over(
            ogeom_math::Axis2::new(Point2::new(0.0, v), ogeom_math::Direction2::X),
            u_dom.0 - 1.0,
            u_dom.1 + 1.0,
        )?
        .into())
    };
    let column_line = |u: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
        Ok(Line2d::over(
            ogeom_math::Axis2::new(Point2::new(u, 0.0), ogeom_math::Direction2::Y),
            v_dom.0 - 1.0,
            v_dom.1 + 1.0,
        )?
        .into())
    };
    // The u-run is a seam in v — the same curve at both rows — and the
    // v-run a seam in u.
    ogeom_algo::attach_seam(
        model,
        &u_edge,
        row_line(v_dom.0)?,
        row_line(v_dom.1)?,
        surface_id,
        ogeom_topo::Location::identity(),
        u_dom,
    )?;
    ogeom_algo::attach_seam(
        model,
        &v_edge,
        column_line(u_dom.1)?,
        column_line(u_dom.0)?,
        surface_id,
        ogeom_topo::Location::identity(),
        v_dom,
    )?;

    let wire = ogeom_algo::make_wire(
        model,
        &[
            u_edge.clone(),
            v_edge.clone(),
            u_edge.reversed(),
            v_edge.reversed(),
        ],
        tol,
    )?
    .shape;
    let face = ogeom_algo::make_face_on(model, surface_id, std::slice::from_ref(&wire), tol)?.shape;
    let centroid = {
        let mut c = Vector::new(0.0, 0.0, 0.0);
        let mut n = 0.0;
        for row in rows {
            for p in row {
                c += p.to_vector();
                n += 1.0;
            }
        }
        Point::from_vector(c / n)
    };
    let mid_u = f64::midpoint(u_dom.0, u_dom.1);
    let mid_v = f64::midpoint(v_dom.0, v_dom.1);
    let s_mid = surface_geo.point_at(mid_u, mid_v, tol)?;
    let (du, dv) = surface_geo.d1_at(mid_u, mid_v, tol)?;
    let face = if du.cross(dv).dot(s_mid - centroid) >= 0.0 {
        face
    } else {
        face.reversed()
    };

    let sewn = sew(model, std::slice::from_ref(&face), tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the closed skin did not close");
    }
    Ok(sewn.shells[0].clone())
}

/// The loft to a point: a wire section closing onto a single apex vertex.
///
/// A circle takes the cone the revolved primitives already build, apex on
/// its axis or refused; a polygon takes exact planar triangle walls, sound
/// for *any* apex — a skew pyramid's walls are still triangles.
fn loft_to_point(
    model: &mut Model,
    section: &Shape,
    apex: &Shape,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !ogeom_algo::is_wire_closed(model, section, tol)? {
        ogeom_bail!(Construction, "a loft section must be closed");
    }
    let apex_point = {
        let Some(data) = model.node(apex).and_then(|n| n.data().as_vertex()) else {
            ogeom_bail!(Construction, "the apex vertex holds no data");
        };
        data.point
    };

    // The circular case: a cone, apex on the axis.
    let edges = explore(model, section, Filter::OfType(ShapeType::Edge))?;
    if edges.len() == 1
        && let Some(data) = model.node(&edges[0]).and_then(|n| n.data().as_edge())
        && let Some(EdgeRepr::Curve3d { curve, .. }) = data.curve3d()
        && let Some(Curve::Circle(c)) = model.geometry().curve(*curve)
    {
        let circle = c.circle();
        let axis = circle.frame().z().vector();
        let rise = apex_point - circle.centre();
        let height = rise.dot(axis);
        if rise.cross(axis).magnitude() > tol.confusion() * 10.0 {
            ogeom_bail!(
                Construction,
                "a circle lofts to a point on its own axis; the oblique cone \
                 needs the skinned machinery — docs/PARITY.md, offset.loft"
            );
        }
        if height.abs() <= tol.confusion() {
            ogeom_bail!(Construction, "the apex sits in the section's own plane");
        }
        let base = if height > 0.0 {
            circle.frame()
        } else {
            Frame::new(
                circle.centre(),
                -circle.frame().z(),
                circle.frame().x(),
                tol,
            )?
        };
        let mut built =
            ogeom_algo::make_cone(model, base, circle.radius(), 0.0, height.abs(), tol)?;
        built.history.generate(section, built.shape.clone());
        built.history.generate(apex, built.shape.clone());
        return Ok(built);
    }

    // The polygonal case: exact triangle walls to a shared apex.
    let mut corners: Vec<Point> = Vec::new();
    for edge in model.ordered_children_of(section)? {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            ogeom_bail!(Construction, "a section edge holds no data");
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            ogeom_bail!(Construction, "a section edge has no curve");
        };
        let Some(Curve::Line(line)) = model.geometry().curve(*curve).cloned() else {
            ogeom_bail!(
                Construction,
                "a mixed or curved section lofts to a point through the \
                 skinned machinery — docs/PARITY.md, offset.loft"
            );
        };
        let t = if edge.orientation() == ogeom_topo::Orientation::Reversed {
            range.1
        } else {
            range.0
        };
        corners.push(ogeom_geom::Curve::Line(line).point_at(t, tol)?);
    }
    if corners.len() < 3 {
        ogeom_bail!(Construction, "a pyramid needs at least three base corners");
    }
    let apex_vertex = ogeom_algo::make_vertex(model, apex_point).shape;
    let base_vertices: Vec<Shape> = corners
        .iter()
        .map(|p| ogeom_algo::make_vertex(model, *p).shape)
        .collect();
    let segment =
        |model: &mut Model, from: (&Shape, Point), to: (&Shape, Point)| -> OgeomResult<Shape> {
            let line = ogeom_geom::LineCurve::segment(from.1, to.1, tol)?;
            let curve: Curve = line.into();
            let domain = curve.domain();
            Ok(make_edge_between(model, curve, domain, from.0, to.0, tol)?.shape)
        };
    let count = corners.len();
    let mut base_edges = Vec::with_capacity(count);
    let mut rails = Vec::with_capacity(count);
    for i in 0..count {
        let next = (i + 1) % count;
        base_edges.push(segment(
            model,
            (&base_vertices[i], corners[i]),
            (&base_vertices[next], corners[next]),
        )?);
        rails.push(segment(
            model,
            (&base_vertices[i], corners[i]),
            (&apex_vertex, apex_point),
        )?);
    }
    let centroid = {
        let mut c = Vector::new(0.0, 0.0, 0.0);
        for p in &corners {
            c += p.to_vector();
        }
        #[allow(clippy::cast_precision_loss)]
        Point::from_vector(c / count as f64 / 4.0 * 3.0 + apex_point.to_vector() / 4.0)
    };
    let planar = |model: &mut Model, pts: [Point; 3], walk: Vec<Shape>| -> OgeomResult<Shape> {
        let n = (pts[1] - pts[0]).cross(pts[2] - pts[0]);
        let m = n.magnitude();
        if m <= tol.confusion() {
            ogeom_bail!(Construction, "a wall of the pyramid is degenerate");
        }
        let mut outward = n / m;
        if outward.dot(pts[0] - centroid) < 0.0 {
            outward = -outward;
        }
        let plane = ogeom_math::Plane::through(pts[0], Direction::new(outward, tol)?);
        let mut reach = 1.0_f64;
        for p in pts {
            reach = reach.max(p.distance(pts[0]) * 2.0);
        }
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let id = model.geometry_mut().add_surface(surface.clone());
        let signed = {
            let (du, dv) = {
                use ogeom_geom::Surface as _;
                surface.d1_at(0.0, 0.0, tol)?
            };
            du.cross(dv).dot(outward) >= 0.0
        };
        let mut wired = Vec::with_capacity(walk.len());
        for used in &walk {
            let (curve, range) = spine_curve_of(model, used)?;
            let Some(pcurve) = ogeom_intersect::exact_pcurve_of(&curve, &surface, tol) else {
                ogeom_bail!(Construction, "a wall edge has no closed-form pcurve");
            };
            ogeom_algo::attach_pcurve(
                model,
                used,
                pcurve,
                id,
                ogeom_topo::Location::identity(),
                range,
            )?;
            wired.push(used.clone());
        }
        let wire = ogeom_algo::make_wire(model, &wired, tol)?.shape;
        let face = ogeom_algo::make_face_on(model, id, std::slice::from_ref(&wire), tol)?.shape;
        Ok(if signed { face } else { face.reversed() })
    };
    let mut faces = Vec::with_capacity(count + 1);
    for i in 0..count {
        let next = (i + 1) % count;
        faces.push(planar(
            model,
            [corners[i], corners[next], apex_point],
            vec![
                base_edges[i].clone(),
                rails[next].clone(),
                rails[i].reversed(),
            ],
        )?);
    }
    // The base cap: all corners, wound against the walls.
    let base_walk: Vec<Shape> = (0..count).rev().map(|i| base_edges[i].reversed()).collect();
    faces.push({
        let n = (corners[1] - corners[0]).cross(corners[2] - corners[0]);
        let mut outward = n / n.magnitude();
        if outward.dot(corners[0] - centroid) < 0.0 {
            outward = -outward;
        }
        let plane = ogeom_math::Plane::through(corners[0], Direction::new(outward, tol)?);
        let mut reach = 1.0_f64;
        for p in &corners {
            reach = reach.max(p.distance(corners[0]) * 2.0);
        }
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let id = model.geometry_mut().add_surface(surface.clone());
        for used in &base_walk {
            let (curve, range) = spine_curve_of(model, used)?;
            let Some(pcurve) = ogeom_intersect::exact_pcurve_of(&curve, &surface, tol) else {
                ogeom_bail!(Construction, "a base edge has no closed-form pcurve");
            };
            ogeom_algo::attach_pcurve(
                model,
                used,
                pcurve,
                id,
                ogeom_topo::Location::identity(),
                range,
            )?;
        }
        let wire = ogeom_algo::make_wire(model, &base_walk, tol)?.shape;
        let face = ogeom_algo::make_face_on(model, id, std::slice::from_ref(&wire), tol)?.shape;
        let signed = {
            use ogeom_geom::Surface as _;
            let (du, dv) = surface.d1_at(0.0, 0.0, tol)?;
            du.cross(dv).dot(outward) >= 0.0
        };
        if signed { face } else { face.reversed() }
    });

    let sewn = sew(model, &faces, tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the pyramid did not close");
    }
    let mut built = make_solid(model, std::slice::from_ref(&sewn.shells[0]))?;
    built.history.generate(section, built.shape.clone());
    built.history.generate(apex, built.shape.clone());
    Ok(built)
}

/// Loft through sections with the start of each row named by the caller.
///
/// [`make_loft_skinned`] leaves alignment to each section's own traversal
/// start; this sibling takes one hint per section — a point near where its
/// row should begin — and rotates each sampling there, which is how a
/// caller untwists a loft whose wires happen to start in different places.
///
/// # Errors
///
/// As [`make_loft_skinned`], and additionally if the hints do not pair up
/// with the sections.
pub fn make_loft_skinned_aligned(
    model: &mut Model,
    sections: &[Shape],
    hints: &[Point],
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if hints.len() != sections.len() {
        ogeom_bail!(
            Construction,
            "{} hints against {} sections; each section names its own start",
            hints.len(),
            sections.len()
        );
    }
    if sections.len() < 2 {
        ogeom_bail!(Construction, "a loft needs at least two sections");
    }
    const AROUND: usize = 48;
    let mut rows: Vec<Vec<Point>> = Vec::with_capacity(sections.len());
    let mut planes: Vec<Plane> = Vec::with_capacity(sections.len());
    for (wire, hint) in sections.iter().zip(hints) {
        if model.kind_of(wire)? != ShapeType::Wire {
            ogeom_bail!(Construction, "a loft section is a closed wire");
        }
        if !ogeom_algo::is_wire_closed(model, wire, tol)? {
            ogeom_bail!(Construction, "a loft section must be closed");
        }
        let Some(plane) = ogeom_algo::find_plane(model, wire, tol)? else {
            ogeom_bail!(Construction, "a loft section must be planar");
        };
        planes.push(plane);
        rows.push(sample_wire_from(model, wire, AROUND, Some(*hint), tol)?);
    }
    let outward0 = {
        let towards = rows[1][0] - rows[0][0];
        let n = planes[0].normal().vector();
        if n.dot(towards) > 0.0 { -n } else { n }
    };
    let outward1 = {
        let towards = rows[rows.len() - 2][0] - rows[rows.len() - 1][0];
        let n = planes[planes.len() - 1].normal().vector();
        if n.dot(towards) > 0.0 { -n } else { n }
    };
    let mut built = skinned_solid(model, &rows, (outward0, outward1), tolerance, tol)?;
    for section in sections {
        built.history.generate(section, built.shape.clone());
    }
    Ok(built)
}

/// Loft a ring through closed planar sections that loop back to the first.
///
/// [`make_loft_skinned`]'s closed sibling: the sections are sampled the same
/// way, the skin runs through all of them and back to the start, C1 across
/// the loop, and there are no caps — the result bounds itself the way a
/// torus does. The sections are *not* repeated: the loop-back is the
/// construction's own.
///
/// The closed join costs freedom: a sparse loop fits only loosely, and the
/// refusal quotes the deviation it honestly reached. A loop that wants a
/// tight tolerance wants sections dense enough to bend around — in
/// practice, a dozen and up.
///
/// # Errors
///
/// As [`make_loft_skinned`], needing at least three sections;
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the closed
/// skin cannot reach the tolerance.
pub fn make_loft_skinned_closed(
    model: &mut Model,
    sections: &[Shape],
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if sections.len() < 3 {
        ogeom_bail!(Construction, "a closed loft needs at least three sections");
    }
    const AROUND: usize = 48;
    let mut rows: Vec<Vec<Point>> = Vec::with_capacity(sections.len());
    for wire in sections {
        if model.kind_of(wire)? != ShapeType::Wire {
            ogeom_bail!(Construction, "a loft section is a closed wire");
        }
        if !ogeom_algo::is_wire_closed(model, wire, tol)? {
            ogeom_bail!(Construction, "a loft section must be closed");
        }
        if ogeom_algo::find_plane(model, wire, tol)?.is_none() {
            ogeom_bail!(Construction, "a loft section must be planar");
        }
        rows.push(sample_wire(model, wire, AROUND, tol)?);
    }
    let mut built = closed_skinned_solid(model, &rows, tolerance, tol)?;
    for section in sections {
        built.history.generate(section, built.shape.clone());
    }
    Ok(built)
}

/// A skinned strip: one open patch of a sweep, with its border edges.
///
/// The wall of a *faceted* profile cannot be one closed skin — a fit cannot
/// speak a corner — so each profile edge sweeps its own strip, cornered at
/// the caller's shared vertices, and the strips weld along their rails by
/// the tolerance the fit honestly carries.
struct SkinnedStrip {
    face: Shape,
    /// The border along the first station, from `corners.0` to `corners.1`.
    bottom: Shape,
    /// The border along the last station, from `corners.2` to `corners.3`.
    top: Shape,
}

/// Skin an open grid of samples — stations by profile-edge samples — into
/// one strip. `corners` are the caller's vertices at (first station, edge
/// start), (first, end), (last, start), (last, end), shared with the
/// neighbouring strips so the wires chain.
fn skinned_strip(
    model: &mut Model,
    rows: &[Vec<Point>],
    corners: (&Shape, &Shape, &Shape, &Shape),
    outward_hint: Point,
    hole: bool,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<SkinnedStrip> {
    use ogeom_geom::Surface as _;
    let fitted = ogeom_geom::fit::fit_surface_grid(rows, 3, tolerance, tol)?;
    if !fitted.met {
        ogeom_bail!(
            NotDone,
            "the strip reached {} against a target of {tolerance}",
            fitted.error
        );
    }
    let error = fitted.error.max(tol.confusion());
    let surface = fitted.curve;
    let (u_knots, v_knots) = (surface.u_knots().clone(), surface.v_knots().clone());
    let (k, l, net) = {
        let grid = surface.grid();
        let net: Vec<Point> = grid.points().iter().map(|w| (*w).point()).collect();
        (grid.u_count(), grid.v_count(), net)
    };
    let point_at = |i: usize, j: usize| -> Point { net[i * l + j] };
    let (u_dom, v_dom) = surface.domain();

    let u_curve = |j: usize| -> OgeomResult<ogeom_geom::Curve> {
        let control: Vec<Point> = (0..k).map(|i| point_at(i, j)).collect();
        Ok(ogeom_geom::Curve::BSpline(ogeom_geom::BSplineCurve::new(
            u_knots.clone(),
            control,
            tol,
        )?))
    };
    let v_curve = |i: usize| -> OgeomResult<ogeom_geom::Curve> {
        let control: Vec<Point> = (0..l).map(|j| point_at(i, j)).collect();
        Ok(ogeom_geom::Curve::BSpline(ogeom_geom::BSplineCurve::new(
            v_knots.clone(),
            control,
            tol,
        )?))
    };
    let surface_geo: SurfaceGeometry = surface.into();
    let surface_id = model.geometry_mut().add_surface(surface_geo.clone());

    let (c00, c10, c01, c11) = corners;
    let bottom = make_edge_between(model, u_curve(0)?, u_dom, c00, c10, tol)?.shape;
    let top = make_edge_between(model, u_curve(l - 1)?, u_dom, c01, c11, tol)?.shape;
    let rail0 = make_edge_between(model, v_curve(0)?, v_dom, c00, c01, tol)?.shape;
    let rail1 = make_edge_between(model, v_curve(k - 1)?, v_dom, c10, c11, tol)?.shape;

    let row_line = |v: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
        Ok(Line2d::over(
            ogeom_math::Axis2::new(Point2::new(0.0, v), ogeom_math::Direction2::X),
            u_dom.0 - 1.0,
            u_dom.1 + 1.0,
        )?
        .into())
    };
    let column_line = |u: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
        Ok(Line2d::over(
            ogeom_math::Axis2::new(Point2::new(u, 0.0), ogeom_math::Direction2::Y),
            v_dom.0 - 1.0,
            v_dom.1 + 1.0,
        )?
        .into())
    };
    ogeom_algo::attach_pcurve(
        model,
        &bottom,
        row_line(v_dom.0)?,
        surface_id,
        ogeom_topo::Location::identity(),
        u_dom,
    )?;
    ogeom_algo::attach_pcurve(
        model,
        &top,
        row_line(v_dom.1)?,
        surface_id,
        ogeom_topo::Location::identity(),
        u_dom,
    )?;
    ogeom_algo::attach_pcurve(
        model,
        &rail0,
        column_line(u_dom.0)?,
        surface_id,
        ogeom_topo::Location::identity(),
        v_dom,
    )?;
    ogeom_algo::attach_pcurve(
        model,
        &rail1,
        column_line(u_dom.1)?,
        surface_id,
        ogeom_topo::Location::identity(),
        v_dom,
    )?;
    // The rails carry the fit's honest budget: the neighbouring strip fitted
    // the same transported corners independently, and the weld between them
    // is only as tight as both fits.
    for edge in [&bottom, &top, &rail0, &rail1] {
        model.widen(edge, ogeom_core::Tolerance::new(error)?)?;
    }

    let wire = ogeom_algo::make_wire(
        model,
        &[
            bottom.clone(),
            rail1.clone(),
            top.reversed(),
            rail0.reversed(),
        ],
        tol,
    )?
    .shape;
    let face = ogeom_algo::make_face_on(model, surface_id, std::slice::from_ref(&wire), tol)?.shape;
    let mid_u = f64::midpoint(u_dom.0, u_dom.1);
    let mid_v = f64::midpoint(v_dom.0, v_dom.1);
    let s_mid = surface_geo.point_at(mid_u, mid_v, tol)?;
    let (du, dv) = surface_geo.d1_at(mid_u, mid_v, tol)?;
    let natural_out = du.cross(dv).dot(s_mid - outward_hint) >= 0.0;
    let face = if natural_out == !hole {
        face
    } else {
        face.reversed()
    };
    Ok(SkinnedStrip { face, bottom, top })
}

/// Loft a solid through many closed planar sections, skinned smoothly.
///
/// The sections are sampled at matched arc-length fractions from their own
/// traversal starts — aligning those starts is the caller's authorship —
/// and the skin holds every section to `tolerance`. The caps are the first
/// and last sections' own planes.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if fewer than
/// two sections, a section is open or not planar;
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if
/// the skin cannot reach the tolerance.
pub fn make_loft_skinned(
    model: &mut Model,
    sections: &[Shape],
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if sections.len() < 2 {
        ogeom_bail!(Construction, "a loft needs at least two sections");
    }
    const AROUND: usize = 48;
    // A trailing vertex is the apex form: the skin narrows to a point and
    // the solid closes there without a cap.
    let apex = match model.kind_of(&sections[sections.len() - 1])? {
        ShapeType::Vertex => {
            if sections.len() < 2 {
                ogeom_bail!(
                    Construction,
                    "a loft to a point needs a section to start from"
                );
            }
            let Some(data) = model
                .node(&sections[sections.len() - 1])
                .and_then(|n| n.data().as_vertex())
            else {
                ogeom_bail!(Construction, "the apex vertex holds no point");
            };
            Some(data.point)
        }
        _ => None,
    };
    let wires = &sections[..sections.len() - usize::from(apex.is_some())];
    let mut rows: Vec<Vec<Point>> = Vec::with_capacity(sections.len());
    let mut cap_planes: Vec<Option<Plane>> = Vec::with_capacity(wires.len());
    for wire in wires {
        if model.kind_of(wire)? != ShapeType::Wire {
            ogeom_bail!(Construction, "a loft section is a wire");
        }
        if !ogeom_algo::is_wire_closed(model, wire, tol)? {
            ogeom_bail!(Construction, "a loft section must be closed");
        }
        // Planarity is a *cap's* requirement, not the fit's: only the
        // sections a cap will stand on must hold a plane. A wavy middle
        // section skins fine.
        cap_planes.push(ogeom_algo::find_plane(model, wire, tol)?);
        rows.push(sample_wire(model, wire, AROUND, tol)?);
    }
    let outward_at =
        |rows: &[Vec<Point>], planes: &[Option<Plane>], end: bool| -> OgeomResult<Vector> {
            let (i, j) = if end {
                (rows.len() - 1, rows.len() - 2)
            } else {
                (0, 1)
            };
            let Some(plane) = &planes[i] else {
                ogeom_bail!(
                    Construction,
                    "a loft's end section must be planar; a cap stands on it"
                );
            };
            let towards = rows[j][0] - rows[i][0];
            let n = plane.normal().vector();
            Ok(if n.dot(towards) > 0.0 { -n } else { n })
        };
    let mut built = if let Some(apex) = apex {
        if rows.len() < 2 {
            // One ring to a point is exact machinery's job when it can be;
            // the skin still needs two rows to shape the wall, so a middle
            // row is interpolated halfway toward the apex.
            let half: Vec<Point> = rows[0]
                .iter()
                .map(|p| Point::from_vector((p.to_vector() + apex.to_vector()) * 0.5))
                .collect();
            rows.push(half);
        }
        let outward0 = outward_at(&rows, &cap_planes, false)?;
        rows.push(vec![apex; AROUND]);
        skinned_solid_to_apex(model, &rows, outward0, tolerance, tol)?
    } else {
        let outward0 = outward_at(&rows, &cap_planes, false)?;
        let outward1 = outward_at(&rows, &cap_planes, true)?;
        skinned_solid(model, &rows, (outward0, outward1), tolerance, tol)?
    };
    for section in sections {
        built.history.generate(section, built.shape.clone());
    }
    Ok(built)
}

/// Sample a closed wire at `count` matched arc-length fractions.
fn sample_wire(
    model: &Model,
    wire: &Shape,
    count: usize,
    tol: Tolerances,
) -> OgeomResult<Vec<Point>> {
    sample_wire_from(model, wire, count, None, tol)
}

/// As [`sample_wire`], with the arc-length origin rotated to the dense
/// sample nearest `start_hint` — how a caller says which point of each
/// section rows up with which, instead of leaning on traversal starts.
fn sample_wire_from(
    model: &Model,
    wire: &Shape,
    count: usize,
    start_hint: Option<Point>,
    tol: Tolerances,
) -> OgeomResult<Vec<Point>> {
    // Dense polyline by traversal, then resample by cumulative length.
    let mut dense: Vec<Point> = Vec::new();
    for edge in explore(model, wire, Filter::OfType(ShapeType::Edge))? {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            ogeom_bail!(Construction, "a section edge holds no data");
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            ogeom_bail!(Construction, "a section edge has no curve");
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            ogeom_bail!(Dangling, "curve is not in this model");
        };
        let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
        for i in 0..64 {
            let f = f64::from(i) / 64.0;
            let t = if reversed {
                range.1 - (range.1 - range.0) * f
            } else {
                range.0 + (range.1 - range.0) * f
            };
            dense.push(geometry.point_at(t, tol)?);
        }
    }
    if let Some(hint) = start_hint {
        let mut best = 0usize;
        let mut held = f64::INFINITY;
        for (i, p) in dense.iter().enumerate() {
            let d = p.distance(hint);
            if d < held {
                held = d;
                best = i;
            }
        }
        dense.rotate_left(best);
    }
    let mut lengths = vec![0.0];
    for w in dense.windows(2) {
        let last = lengths[lengths.len() - 1];
        lengths.push(last + w[0].distance(w[1]));
    }
    let closing = dense[dense.len() - 1].distance(dense[0]);
    let total = lengths[lengths.len() - 1] + closing;
    let mut out = Vec::with_capacity(count);
    let mut cursor = 0usize;
    for s in 0..count {
        #[allow(clippy::cast_precision_loss)]
        let target = total * (s as f64) / (count as f64);
        while cursor + 1 < lengths.len() && lengths[cursor + 1] < target {
            cursor += 1;
        }
        let (a, b) = (dense[cursor], dense[(cursor + 1) % dense.len()]);
        let la = lengths[cursor];
        let lb = if cursor + 1 < lengths.len() {
            lengths[cursor + 1]
        } else {
            total
        };
        let f = if lb > la {
            (target - la) / (lb - la)
        } else {
            0.0
        };
        out.push(a + (b - a) * f.clamp(0.0, 1.0));
    }
    Ok(out)
}

/// Sweep a circular profile along a free-form spine, skinned.
///
/// Frames along the spine are rotation-minimizing (the double-reflection
/// construction), so the tube neither twists nor kinks where the spine
/// bends; the skin holds the sampled circles to `tolerance`, and the caps
/// sit perpendicular to the spine's ends.
///
/// # Errors
///
/// As [`make_pipe`], plus [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if
/// the skin cannot reach the
/// tolerance.
pub fn make_pipe_skinned(
    model: &mut Model,
    spine: &Shape,
    radius: f64,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(Construction, "a pipe of radius {radius} holds nothing");
    }
    let (curve, range) = {
        let Some(data) = model.node(spine).and_then(|n| n.data().as_edge()) else {
            ogeom_bail!(Construction, "a pipe runs along an edge");
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            ogeom_bail!(Construction, "the spine has no curve");
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            ogeom_bail!(Dangling, "curve is not in this model");
        };
        (geometry.clone(), *range)
    };
    const STATIONS: usize = 33;
    const AROUND: usize = 40;
    let mut stations: Vec<SpineStation> = Vec::with_capacity(STATIONS);
    for i in 0..STATIONS {
        #[allow(clippy::cast_precision_loss)]
        let t = range.0 + (range.1 - range.0) * (i as f64) / ((STATIONS - 1) as f64);
        let p = curve.point_at(t, tol)?;
        let d = curve.d1_at(t, tol)?;
        let m = d.magnitude();
        if m <= tol.confusion() {
            ogeom_bail!(Construction, "the spine is degenerate at {t}");
        }
        stations.push(SpineStation {
            at: p,
            tangent: d / m,
        });
    }
    let normals = rmf_normals(&stations);
    let mut rows: Vec<Vec<Point>> = Vec::with_capacity(STATIONS);
    for (i, station) in stations.iter().enumerate() {
        let x = normals[i];
        let y = station.tangent.cross(x);
        let mut row = Vec::with_capacity(AROUND);
        for a in 0..AROUND {
            #[allow(clippy::cast_precision_loss)]
            let ang = core::f64::consts::TAU * (a as f64) / (AROUND as f64);
            row.push(station.at + (x * ang.cos() + y * ang.sin()) * radius);
        }
        rows.push(row);
    }
    let mut built = skinned_solid(
        model,
        &rows,
        (-stations[0].tangent, stations[STATIONS - 1].tangent),
        tolerance,
        tol,
    )?;
    built.history.generate(spine, built.shape.clone());
    Ok(built)
}

/// One sampled spine station: where the spine is and which way it runs.
#[derive(Clone, Copy)]
struct SpineStation {
    at: Point,
    /// The unit tangent, in the direction of travel.
    tangent: Vector,
}

/// One profile wire's closed shell round the spine: smooth wires skin as a
/// single closed face, faceted ones as one ring strip per facet.
#[allow(clippy::too_many_arguments, reason = "one frame, spelled out")]
fn closed_loop_shell(
    model: &mut Model,
    profile_loop: &Shape,
    edges: &[Shape],
    smooth: bool,
    stations: &[SpineStation],
    normals: &[Vector],
    frame0: (Point, Vector),
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    const AROUND: usize = 40;
    let (origin, x0) = frame0;
    let t0 = stations[0].tangent;
    let y0 = t0.cross(x0);
    if !smooth {
        // A faceted profile: one ring strip per profile edge — a fit cannot
        // speak a corner, so each facet gets its own v-closed skin and the
        // strips sew along the corner loops they share within tolerance.
        const ALONG_EDGE: usize = 8;
        // Outward for a ring strip means away from the spine's own line,
        // not from the loop's centroid — a ring's inner side *faces* the
        // centroid. The hint is the station the strip's midpoint rides.
        let mid_station = stations[stations.len() / 2].at;
        let mut faces = Vec::with_capacity(edges.len());
        for edge in edges {
            let (curve, range) = spine_curve_of(model, edge)?;
            let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
            let mut flat_row: Vec<(f64, f64)> = Vec::with_capacity(ALONG_EDGE + 1);
            for kk in 0..=ALONG_EDGE {
                #[allow(clippy::cast_precision_loss)]
                let f = (kk as f64) / (ALONG_EDGE as f64);
                let t = if reversed {
                    range.1 - (range.1 - range.0) * f
                } else {
                    range.0 + (range.1 - range.0) * f
                };
                let p = curve.point_at(t, tol)?;
                flat_row.push(((p - origin).dot(x0), (p - origin).dot(y0)));
            }
            // rows[j = station][i = across the facet].
            let rows: Vec<Vec<Point>> = stations
                .iter()
                .enumerate()
                .map(|(i, station)| {
                    let x = normals[i];
                    let y = station.tangent.cross(x);
                    flat_row
                        .iter()
                        .map(|(a, b)| station.at + x * *a + y * *b)
                        .collect()
                })
                .collect();
            faces.push(skinned_ring_strip(
                model,
                &rows,
                mid_station,
                tolerance,
                tol,
            )?);
        }
        let sewn = sew(model, &faces, tol)?;
        if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
            ogeom_bail!(Construction, "the faceted ring did not close");
        }
        return Ok(sewn.shells[0].clone());
    }
    let samples = sample_wire(model, profile_loop, AROUND, tol)?;
    let flat: Vec<(f64, f64)> = samples
        .iter()
        .map(|p| ((*p - origin).dot(x0), (*p - origin).dot(y0)))
        .collect();
    let rows: Vec<Vec<Point>> = stations
        .iter()
        .enumerate()
        .map(|(i, station)| {
            let x = normals[i];
            let y = station.tangent.cross(x);
            flat.iter()
                .map(|(a, b)| station.at + x * *a + y * *b)
                .collect()
        })
        .collect();
    closed_skinned_shell(model, &rows, tolerance, tol)
}

/// Rotation-minimizing normals along the stations, by double reflection:
/// reflect in each chord's plane, then in the plane bisecting the tangents.
/// Self-contained — it needs only the station list — and shared by every
/// sweep that must not twist where its spine bends.
fn rmf_normals(stations: &[SpineStation]) -> Vec<Vector> {
    let mut normals: Vec<Vector> = Vec::with_capacity(stations.len());
    let t0 = stations[0].tangent;
    let seed = if t0.cross(ogeom_math::Vector::Z).magnitude() > 0.5 {
        ogeom_math::Vector::Z
    } else {
        ogeom_math::Vector::X
    };
    let n0 = {
        let v = seed - t0 * seed.dot(t0);
        v / v.magnitude()
    };
    normals.push(n0);
    for i in 1..stations.len() {
        let (p0, t0) = (stations[i - 1].at, stations[i - 1].tangent);
        let (p1, t1) = (stations[i].at, stations[i].tangent);
        let n = normals[i - 1];
        let v1 = p1 - p0;
        let c1 = v1.dot(v1);
        if c1 <= 1e-20 {
            // A corner's twin station: no travel to reflect through. The
            // normal carries straight across, which is exactly the parallel
            // transport a mitred joint's mirror symmetry needs.
            normals.push(n);
            continue;
        }
        let nl = n - v1 * (2.0 / c1 * v1.dot(n));
        let tl = t0 - v1 * (2.0 / c1 * v1.dot(t0));
        let v2 = t1 - tl;
        let c2 = v2.dot(v2);
        let next = if c2 > 1e-20 {
            nl - v2 * (2.0 / c2 * v2.dot(nl))
        } else {
            nl
        };
        normals.push(next / next.magnitude());
    }
    normals
}

/// Sweep a planar profile — a wire, or a face whose holes ride along —
/// down an arbitrary spine, one skinned wall per profile loop.
///
/// The spine may be a single edge or a wire of edges of any curve the
/// vocabulary evaluates — lines, arcs, splines, helices. Frames along it are
/// rotation-minimizing by default (the double-reflection construction), so
/// the profile neither twists nor kinks where the spine bends; `frenet`
/// asks for the Frenet frame instead, which turns with the spine's own
/// curvature — the law a thread wants. Stations are placed by each edge's
/// own turning, the skin holds every transported section to `tolerance`,
/// and the caps sit perpendicular to the spine's ends, holes and all.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// profile is not planar, leans along the spine, or does not sit at the
/// spine's start; if the spine is closed (the loop-back is the closed-skin
/// milestone — docs/PARITY.md, offset.sweeps); or if `frenet` is asked of a
/// spine that never bends.
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the skin
/// cannot reach the tolerance.
pub fn make_pipe_shell(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    frenet: bool,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    const AROUND: usize = 40;

    let stations = shell_stations(model, spine, tol)?;
    // Corners: twin stations standing on one point with different headings.
    let kinks: Vec<usize> = (0..stations.len() - 1)
        .filter(|&i| {
            stations[i].at.distance(stations[i + 1].at) <= tol.confusion()
                && (stations[i]
                    .tangent
                    .cross(stations[i + 1].tangent)
                    .magnitude()
                    > tol.angular()
                    || stations[i].tangent.dot(stations[i + 1].tangent) < 0.0)
        })
        .collect();
    if stations[0].at.distance(stations[stations.len() - 1].at) <= tol.confusion() * 10.0 {
        if !kinks.is_empty() {
            ogeom_bail!(
                Construction,
                "a closed spine's sharp corners are still owed their mitres \
                 — docs/PARITY.md, offset.sweeps"
            );
        }
        return closed_pipe_shell(model, profile, spine, stations, frenet, tolerance, tol);
    }
    if frenet && !kinks.is_empty() {
        ogeom_bail!(
            Construction,
            "a Frenet frame has no direction at a corner; sweep a cornered \
             spine with the rotation-minimizing frame"
        );
    }
    let normals = if frenet {
        frenet_normals(&stations, tol)?
    } else {
        rmf_normals(&stations)
    };
    // At each corner both twin sections are thrown onto the bisector plane
    // along their own tangents; the mirror symmetry of the rotation-
    // minimizing frame lands them on one ring, the mitre both runs share.
    let mitre: Vec<Option<(Point, Vector)>> = {
        let mut out: Vec<Option<(Point, Vector)>> = vec![None; stations.len()];
        for &k in &kinks {
            let n = stations[k].tangent + stations[k + 1].tangent;
            if n.magnitude() <= tol.angular() {
                ogeom_bail!(
                    Construction,
                    "the spine doubles straight back on itself; no mitre \
                     plane divides that corner"
                );
            }
            out[k] = Some((stations[k].at, n));
            out[k + 1] = Some((stations[k + 1].at, n));
        }
        out
    };
    let runs: Vec<(usize, usize)> = {
        let mut out = Vec::with_capacity(kinks.len() + 1);
        let mut start = 0;
        for &k in &kinks {
            out.push((start, k));
            start = k + 1;
        }
        out.push((start, stations.len() - 1));
        out
    };
    // A mitred end is a *shear*: the honest wall is the run's own surface
    // trimmed by the mitre plane, which for a straight leg is exactly the
    // ruled skin between its two end rings. A curved leg's trim is not a
    // loft of its rows, so a corner against one stays refused by name.
    let straight = |rs: usize, re: usize| -> bool {
        let t0 = stations[rs].tangent;
        (rs..=re).all(|i| stations[i].tangent.cross(t0).magnitude() <= tol.angular())
    };
    for &(rs, re) in &runs {
        if (mitre[rs].is_some() || mitre[re].is_some()) && !straight(rs, re) {
            ogeom_bail!(
                Construction,
                "a sharp corner against a curved leg is still owed its \
                 mitre; straight legs mitre exactly — docs/PARITY.md, \
                 offset.sweeps"
            );
        }
    }
    // The stations a run's skin actually interpolates: a straight run is
    // its two end rings, ruled — the trimmed prism itself.
    let run_stations = |rs: usize, re: usize| -> Vec<usize> {
        if straight(rs, re) {
            vec![rs, re]
        } else {
            (rs..=re).collect()
        }
    };

    // The profile's loops: a face contributes every wire, holes included;
    // a bare wire is one loop.
    let loops: Vec<Shape> = match model.kind_of(profile)? {
        ShapeType::Face => explore(model, profile, Filter::OfType(ShapeType::Wire))?,
        ShapeType::Wire => vec![profile.clone()],
        other => ogeom_bail!(
            Construction,
            "a pipe shell sweeps a planar wire or face, not a {other:?}"
        ),
    };
    if loops.is_empty() {
        ogeom_bail!(Construction, "the profile has no loop to sweep");
    }
    let Some(plane) = ogeom_algo::find_plane(model, profile, tol)? else {
        ogeom_bail!(Construction, "a pipe shell sweeps a planar profile");
    };
    let t0 = stations[0].tangent;
    if plane.normal().vector().cross(t0).magnitude() > 1e-9 {
        ogeom_bail!(
            Construction,
            "the profile leans along its spine; a pipe shell runs square to \
             the start"
        );
    }
    if plane.distance_to(stations[0].at) > tol.confusion() * 100.0 {
        ogeom_bail!(
            Construction,
            "the profile does not sit at the spine's start"
        );
    }

    // Transport: each loop expressed in the start frame's own 2D
    // coordinates, then re-expressed in every station's frame. A loop of one
    // smooth closed edge skins as one wall; a faceted loop skins one strip
    // per edge, cornered at shared vertices, because no single fit can
    // speak a corner.
    let x0 = normals[0];
    let y0 = t0.cross(x0);
    let origin = stations[0].at;
    let flat = |p: Point| -> (f64, f64) { ((p - origin).dot(x0), (p - origin).dot(y0)) };
    let place = |i: usize, (a, b): (f64, f64)| -> Point {
        let x = normals[i];
        let y = stations[i].tangent.cross(x);
        let p = stations[i].at + x * a + y * b;
        match mitre[i] {
            Some((corner, n)) => {
                let t = stations[i].tangent;
                p + t * ((corner - p).dot(n) / t.dot(n))
            }
            None => p,
        }
    };
    let last = stations.len() - 1;

    enum LoopWall {
        Ring {
            ring0: Shape,
            ring1: Shape,
        },
        Chain {
            bottoms: Vec<Shape>,
            tops: Vec<Shape>,
        },
    }
    let mut faces: Vec<Shape> = Vec::new();
    let mut ends: Vec<LoopWall> = Vec::with_capacity(loops.len());
    for (li, wire) in loops.iter().enumerate() {
        let hole = li != 0;
        let edges = model.ordered_children_of(wire)?;
        let single_smooth = edges.len() == 1 && {
            let (curve, _) = spine_curve_of(model, &edges[0])?;
            !matches!(curve, ogeom_geom::Curve::Line(_))
                && ogeom_algo::edge_vertices(model, &edges[0])?.is_some_and(|(a, b)| a.is_same(&b))
        };
        if single_smooth {
            let samples = sample_wire(model, wire, AROUND, tol)?;
            let flat_row: Vec<(f64, f64)> = samples.iter().map(|p| flat(*p)).collect();
            // One wall per smooth run: a fit across a corner speaks nothing,
            // and the twin stations put both runs' boundary rows on the one
            // mitred ring, where the sew joins them.
            let mut ring0: Option<Shape> = None;
            let mut ring1: Option<Shape> = None;
            for &(rs, re) in &runs {
                let rows: Vec<Vec<Point>> = run_stations(rs, re)
                    .into_iter()
                    .map(|i| flat_row.iter().map(|ab| place(i, *ab)).collect())
                    .collect();
                let wall = skinned_wall(model, &rows, tolerance, tol)?;
                faces.push(if hole {
                    wall.face.reversed()
                } else {
                    wall.face.clone()
                });
                if ring0.is_none() {
                    ring0 = Some(wall.ring0);
                }
                ring1 = Some(wall.ring1);
            }
            let (Some(ring0), Some(ring1)) = (ring0, ring1) else {
                ogeom_bail!(Construction, "the sweep produced no wall");
            };
            ends.push(LoopWall::Ring { ring0, ring1 });
        } else {
            // Shared corner vertices at both ends of every edge junction.
            let count = edges.len();
            let mut corner_flat: Vec<(f64, f64)> = Vec::with_capacity(count);
            for edge in &edges {
                let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
                    ogeom_bail!(Construction, "a profile edge has no vertices");
                };
                let start = if edge.orientation() == ogeom_topo::Orientation::Reversed {
                    b
                } else {
                    a
                };
                let Some(data) = model.node(&start).and_then(|n| n.data().as_vertex()) else {
                    ogeom_bail!(Construction, "a profile vertex holds no data");
                };
                corner_flat.push(flat(data.point));
            }
            let make_corners = |model: &mut Model, station: usize| -> Vec<Shape> {
                corner_flat
                    .iter()
                    .map(|ab| ogeom_algo::make_vertex(model, place(station, *ab)).shape)
                    .collect()
            };
            // Corner vertex sets at every run boundary; a kink's twin
            // stations land on the same mitred points, so both runs take
            // the same vertex objects.
            let mut corners_at: Vec<Option<Vec<Shape>>> = vec![None; stations.len()];
            corners_at[0] = Some(make_corners(model, 0));
            corners_at[last] = Some(make_corners(model, last));
            for &k in &kinks {
                let set = make_corners(model, k);
                corners_at[k] = Some(set.clone());
                corners_at[k + 1] = Some(set);
            }

            // The loop's own centroid line, for orienting each strip.
            let hint_flat = {
                let mut a = 0.0;
                let mut b = 0.0;
                for (fa, fb) in &corner_flat {
                    a += fa;
                    b += fb;
                }
                #[allow(clippy::cast_precision_loss)]
                let n = count as f64;
                (a / n, b / n)
            };

            let mut bottoms = Vec::with_capacity(count);
            let mut tops = Vec::with_capacity(count);
            for (ei, edge) in edges.iter().enumerate() {
                let (curve, range) = spine_curve_of(model, edge)?;
                let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
                const ALONG_EDGE: usize = 8;
                let mut flat_row: Vec<(f64, f64)> = Vec::with_capacity(ALONG_EDGE + 1);
                for k in 0..=ALONG_EDGE {
                    #[allow(clippy::cast_precision_loss)]
                    let f = (k as f64) / (ALONG_EDGE as f64);
                    let t = if reversed {
                        range.1 - (range.1 - range.0) * f
                    } else {
                        range.0 + (range.1 - range.0) * f
                    };
                    flat_row.push(flat(curve.point_at(t, tol)?));
                }
                let next = (ei + 1) % count;
                let mut bottom: Option<Shape> = None;
                let mut top: Option<Shape> = None;
                for &(rs, re) in &runs {
                    let rows: Vec<Vec<Point>> = run_stations(rs, re)
                        .into_iter()
                        .map(|i| flat_row.iter().map(|ab| place(i, *ab)).collect())
                        .collect();
                    let mid_i = usize::midpoint(rs, re);
                    let hint = {
                        let x = normals[mid_i];
                        let y = stations[mid_i].tangent.cross(x);
                        stations[mid_i].at + x * hint_flat.0 + y * hint_flat.1
                    };
                    let (Some(from), Some(to)) = (&corners_at[rs], &corners_at[re]) else {
                        ogeom_bail!(Construction, "a run boundary has no corners");
                    };
                    let strip = skinned_strip(
                        model,
                        &rows,
                        (&from[ei], &from[next], &to[ei], &to[next]),
                        hint,
                        hole,
                        tolerance,
                        tol,
                    )?;
                    faces.push(strip.face.clone());
                    if bottom.is_none() {
                        bottom = Some(strip.bottom);
                    }
                    top = Some(strip.top);
                }
                let (Some(bottom), Some(top)) = (bottom, top) else {
                    ogeom_bail!(Construction, "the sweep produced no strip");
                };
                bottoms.push(bottom);
                tops.push(top);
            }
            ends.push(LoopWall::Chain { bottoms, tops });
        }
    }

    // A cap per end: one plane, one wire per loop, each edge's pcurve the
    // exact projection of its control net into the plane's chart.
    for end in 0..2 {
        let (at, outward) = if end == 0 {
            (stations[0].at, -stations[0].tangent)
        } else {
            (stations[last].at, stations[last].tangent)
        };
        let cap_plane = Plane::through(at, Direction::new(outward, tol)?);
        let mut loop_edges: Vec<Vec<Shape>> = Vec::with_capacity(ends.len());
        for wall in &ends {
            loop_edges.push(match wall {
                LoopWall::Ring { ring0, ring1 } => {
                    vec![if end == 0 {
                        ring0.clone()
                    } else {
                        ring1.clone()
                    }]
                }
                LoopWall::Chain { bottoms, tops } => {
                    if end == 0 {
                        bottoms.clone()
                    } else {
                        tops.clone()
                    }
                }
            });
        }
        let mut reach = 1.0_f64;
        for edges in &loop_edges {
            for edge in edges {
                let (curve, range) = spine_curve_of(model, edge)?;
                for t in 0..8 {
                    let p =
                        curve.point_at(range.0 + (range.1 - range.0) * f64::from(t) / 8.0, tol)?;
                    reach = reach.max(p.distance(at) * 2.0);
                }
            }
        }
        let cap_surface: SurfaceGeometry =
            PlaneSurface::over(cap_plane, (-reach, reach), (-reach, reach))?.into();
        let mut wires: Vec<Shape> = Vec::with_capacity(loop_edges.len());
        for edges in &loop_edges {
            wires.push(ogeom_algo::make_wire(model, edges, tol)?.shape);
        }
        let face = ogeom_algo::make_face(model, cap_surface.clone(), &wires, tol)?.shape;
        let cap_id = {
            let Some(node) = model.node(&face) else {
                ogeom_bail!(Dangling, "the cap just built is not in this model");
            };
            let ogeom_topo::NodeData::Face(data) = node.data() else {
                ogeom_bail!(Construction, "the cap holds no face data");
            };
            data.surface
        };
        let frame = cap_plane.frame();
        for edges in &loop_edges {
            for edge in edges {
                let (curve, range) = spine_curve_of(model, edge)?;
                let ogeom_geom::Curve::BSpline(bs) = &curve else {
                    ogeom_bail!(Construction, "a swept ring is not a spline");
                };
                // A planar polynomial spline's chart image is the same-degree
                // spline of the projected control points — affine, so exact.
                let control2: Vec<Point2> = bs
                    .control_points()
                    .iter()
                    .map(|w| {
                        let local = frame.to_local(w.point());
                        Point2::new(local.x, local.y)
                    })
                    .collect();
                let pcurve: ogeom_geom::PlanarCurve =
                    ogeom_geom::BSpline2d::new(bs.knots().clone(), control2, tol)?.into();
                ogeom_algo::attach_pcurve(
                    model,
                    edge,
                    pcurve,
                    cap_id,
                    ogeom_topo::Location::identity(),
                    range,
                )?;
            }
        }
        faces.push(face);
    }

    let sewn = sew(model, &faces, tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the pipe shell did not close");
    }
    let mut built = make_solid(model, std::slice::from_ref(&sewn.shells[0]))?;
    built.history.generate(profile, built.shape.clone());
    built.history.generate(spine, built.shape.clone());
    Ok(built)
}

/// An edge's 3D curve and range, cloned out of the model.
fn spine_curve_of(model: &Model, edge: &Shape) -> OgeomResult<(ogeom_geom::Curve, (f64, f64))> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        ogeom_bail!(Construction, "an edge holds no data");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        ogeom_bail!(Construction, "an edge has no curve");
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    Ok((geometry.clone(), *range))
}

/// The pipe shell around a spine that loops back on itself: one wall,
/// closed both ways round, no caps at all.
///
/// The frames are rotation-minimizing with the loop's holonomy paid off:
/// transported round a closed spine, the frame comes home twisted by some
/// angle, and that twist is spread back along the arc so the last station's
/// frame *is* the first's — without it the closed fit fights a helical
/// grid. The profile must be one smooth closed loop; a faceted profile's
/// strips and a holed profile's nested shells are still owed, and the
/// Frenet law on a closed loop is not carried yet.
fn closed_pipe_shell(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    mut stations: Vec<SpineStation>,
    frenet: bool,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    // The walk visits the join twice; the loop owns it once.
    stations.pop();
    if stations.len() < 3 {
        ogeom_bail!(Construction, "a closed spine needs room to turn");
    }
    // A sharp corner turns the section through a finite angle over no arc at
    // all, which no skin can follow: the mitred ring is per-edge closed
    // strips, still owed.
    for i in 0..stations.len() {
        let next = &stations[(i + 1) % stations.len()];
        if stations[i].tangent.dot(next.tangent) < 0.9 {
            ogeom_bail!(
                Construction,
                "a closed spine with a sharp corner needs the mitred strips, \
                 which are still owed — docs/PARITY.md, offset.sweeps; round \
                 the corner and the ring sweeps"
            );
        }
    }

    let loops: Vec<Shape> = match model.kind_of(profile)? {
        ShapeType::Face => explore(model, profile, Filter::OfType(ShapeType::Wire))?,
        ShapeType::Wire => vec![profile.clone()],
        other => ogeom_bail!(
            Construction,
            "a pipe shell sweeps a planar wire or face, not a {other:?}"
        ),
    };
    // Every wire sweeps its own closed shell: the outer boundary first,
    // each hole a void tunnel inside it.
    let profile_loop = &loops[0];
    let edges = model.ordered_children_of(profile_loop)?;
    let smooth = edges.len() == 1
        && ogeom_algo::edge_vertices(model, &edges[0])?.is_some_and(|(a, b)| a.is_same(&b));
    let Some(plane) = ogeom_algo::find_plane(model, profile, tol)? else {
        ogeom_bail!(Construction, "a pipe shell sweeps a planar profile");
    };
    let t0 = stations[0].tangent;
    if plane.normal().vector().cross(t0).magnitude() > 1e-9 {
        ogeom_bail!(
            Construction,
            "the profile leans along its spine; a pipe shell runs square to \
             the start"
        );
    }
    if plane.distance_to(stations[0].at) > tol.confusion() * 100.0 {
        ogeom_bail!(
            Construction,
            "the profile does not sit at the spine's start"
        );
    }

    // Frames with the loop's mismatch paid off: carry once more back to the
    // start, read the twist between departure and return, and spread it
    // along the arc. Rotation-minimizing frames owe this for their
    // holonomy; the Frenet law owes it too, because straight stretches
    // carry the frame through by continuation and the continuation is
    // path-dependent. One reconciliation serves both.
    let mut normals: Vec<Vector> = {
        let mut extended = stations.clone();
        extended.push(stations[0]);
        if frenet {
            // Wired and measured: the Frenet frames reconcile at the join,
            // but the faceted strips built on them refuse to close on the
            // spines that would prove them. Refused until that is
            // understood rather than shipped hoping.
            ogeom_bail!(
                Construction,
                "the Frenet law on a closed spine is not carried yet; use \
                 the rotation-minimizing default — docs/PARITY.md, \
                 offset.sweeps"
            );
        }
        let carried = rmf_normals(&extended);
        let (n0, n_home) = (carried[0], carried[carried.len() - 1]);
        let twist = (n0.cross(n_home).dot(t0)).atan2(n0.dot(n_home));
        let mut lengths = vec![0.0_f64];
        for pair in extended.windows(2) {
            let last = lengths[lengths.len() - 1];
            lengths.push(last + pair[0].at.distance(pair[1].at));
        }
        let total = lengths[lengths.len() - 1];
        carried
            .iter()
            .take(stations.len())
            .enumerate()
            .map(|(i, n)| {
                let phi = -twist * lengths[i] / total;
                let t = extended[i].tangent;
                *n * phi.cos() + t.cross(*n) * phi.sin()
            })
            .collect()
    };
    for (n, station) in normals.iter_mut().zip(&stations) {
        // Re-square each corrected normal against its own tangent.
        let v = *n - station.tangent * n.dot(station.tangent);
        *n = v / v.magnitude();
    }

    let x0 = normals[0];
    let origin = stations[0].at;
    let mut shells: Vec<Shape> = Vec::with_capacity(loops.len());
    for (li, wire) in loops.iter().enumerate() {
        let wire_edges = model.ordered_children_of(wire)?;
        let wire_smooth = wire_edges.len() == 1
            && ogeom_algo::edge_vertices(model, &wire_edges[0])?
                .is_some_and(|(a, b)| a.is_same(&b));
        let shell = closed_loop_shell(
            model,
            wire,
            &wire_edges,
            wire_smooth,
            &stations,
            &normals,
            (origin, x0),
            tolerance,
            tol,
        )?;
        // A void's faces leave the material toward the tunnel: reversed
        // against the outward orientation every shell is built with.
        shells.push(if li == 0 { shell } else { shell.reversed() });
    }
    let mut built = make_solid(model, &shells)?;
    built.history.generate(profile, built.shape.clone());
    built.history.generate(spine, built.shape.clone());
    let _ = (smooth, profile_loop, edges);
    Ok(built)
}

/// Sample a spine — one edge or a wire of them — into stations, each edge
/// given a station count by its own turning.
fn shell_stations(model: &Model, spine: &Shape, tol: Tolerances) -> OgeomResult<Vec<SpineStation>> {
    let edges: Vec<Shape> = match model.kind_of(spine)? {
        ShapeType::Edge => vec![spine.clone()],
        ShapeType::Wire => model.ordered_children_of(spine)?,
        other => ogeom_bail!(
            Construction,
            "a pipe shell runs along an edge or a wire, not a {other:?}"
        ),
    };
    if edges.is_empty() {
        ogeom_bail!(Construction, "the spine has no edge to run along");
    }
    let mut stations: Vec<SpineStation> = Vec::new();
    for edge in &edges {
        let (curve, range) = {
            let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
                ogeom_bail!(Construction, "a spine edge holds no data");
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                ogeom_bail!(Construction, "a spine edge has no curve");
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                ogeom_bail!(Dangling, "curve is not in this model");
            };
            (geometry.clone(), *range)
        };
        let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
        // Stations by turning: sample tangents coarsely, sum the angles, and
        // give each edge enough stations that no step turns more than a few
        // degrees. A straight edge keeps a healthy minimum for the fit.
        let turning = {
            let mut sum = 0.0_f64;
            let mut last: Option<Vector> = None;
            for i in 0..=16 {
                let t = range.0 + (range.1 - range.0) * f64::from(i) / 16.0;
                let d = curve.d1_at(t, tol)?;
                let m = d.magnitude();
                if m <= tol.confusion() {
                    continue;
                }
                let u = d / m;
                if let Some(prev) = last {
                    sum += prev.dot(u).clamp(-1.0, 1.0).acos() * 16.0 / 16.0;
                }
                last = Some(u);
            }
            sum
        };
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let count = (turning / (core::f64::consts::TAU / 64.0)).ceil().max(8.0) as usize;
        for i in 0..=count {
            #[allow(clippy::cast_precision_loss)]
            let f = (i as f64) / (count as f64);
            let t = if reversed {
                range.1 - (range.1 - range.0) * f
            } else {
                range.0 + (range.1 - range.0) * f
            };
            let p = curve.point_at(t, tol)?;
            let d = curve.d1_at(t, tol)?;
            let m = d.magnitude();
            if m <= tol.confusion() {
                ogeom_bail!(Construction, "the spine is degenerate at {t}");
            }
            let tangent = if reversed { -(d / m) } else { d / m };
            if let Some(prev) = stations.last()
                && prev.at.distance(p) <= tol.confusion()
                && prev.tangent.cross(tangent).magnitude() <= tol.angular()
                && prev.tangent.dot(tangent) > 0.0
            {
                continue;
            }
            // A station coincident with the last but heading elsewhere is a
            // *corner*: both stations stay, a twin pair the sweep mitres.
            stations.push(SpineStation { at: p, tangent });
        }
    }
    if stations.len() < 2 {
        ogeom_bail!(Construction, "the spine collapses to a point");
    }
    Ok(stations)
}

/// Frenet normals: each station's frame turns with the spine's own
/// curvature, read from the tangents' finite differences. Straight runs
/// carry the last bending station's normal forward; a spine that never
/// bends has no Frenet frame at all and is refused by name.
fn frenet_normals(stations: &[SpineStation], tol: Tolerances) -> OgeomResult<Vec<Vector>> {
    let mut normals: Vec<Option<Vector>> = Vec::with_capacity(stations.len());
    for i in 0..stations.len() {
        let (before, after) = (
            &stations[i.saturating_sub(1)],
            &stations[(i + 1).min(stations.len() - 1)],
        );
        let dt = after.tangent - before.tangent;
        let t = stations[i].tangent;
        let bend = dt - t * dt.dot(t);
        let m = bend.magnitude();
        normals.push(if m > tol.angular().max(1e-9) {
            Some(bend / m)
        } else {
            None
        });
    }
    // Carry forward, then backward, so straight lead-ins take the first
    // bend's frame rather than none.
    let mut carried: Vec<Vector> = Vec::with_capacity(stations.len());
    let mut last: Option<Vector> = None;
    for n in &normals {
        if let Some(n) = n {
            last = Some(*n);
        }
        carried.push(last.unwrap_or(Vector::new(0.0, 0.0, 0.0)));
    }
    let mut ahead: Option<Vector> = None;
    for i in (0..stations.len()).rev() {
        if let Some(n) = normals[i] {
            ahead = Some(n);
        } else if carried[i].magnitude() < 0.5
            && let Some(n) = ahead
        {
            carried[i] = n;
        }
    }
    if carried.iter().any(|n| n.magnitude() < 0.5) {
        ogeom_bail!(
            Construction,
            "a straight spine has no Frenet frame; use the \
             rotation-minimizing default"
        );
    }
    Ok(carried)
}

// --- the evolved shape -------------------------------------------------------

/// One station of the spine: where it is, which way it runs, and how it gets
/// there.
struct Station {
    /// Where the traversal enters and leaves this edge.
    from: Point,
    to: Point,
    /// The unit tangent at each end, in the direction of travel.
    tangent_in: Vector,
    tangent_out: Vector,
    /// A straight run, or a turn about an axis through an angle.
    turn: Option<(ogeom_math::Axis, f64)>,
}

/// Sweep a profile along a spine, the way a moulding runs round a frame.
///
/// The spine is a **planar** wire, or a planar face whose outer wire is taken.
/// The profile is a wire standing in a plane that contains the spine's own
/// normal, positioned where the spine starts. What comes back is what the
/// profile sweeps out as it travels the spine, always square to it:
///
/// - a straight spine edge extrudes the profile — a prism;
/// - a circular one turns it about that arc's own axis — a revolution;
/// - and each corner between them turns it about the corner, through exactly
///   the angle the spine turns there, which is the join the 2D offset makes
///   for the same reason.
///
/// Every piece is exact: the surfaces are the ones a prism and a revolution
/// give for the profile's own curves, and nothing is fitted. The pieces are
/// then unioned, which is the assembly's real name — consecutive pieces meet
/// on the *same* placed profile, and a coincident face is what the boolean
/// identifies rather than probes across.
///
/// # Volume or shell
///
/// The result is always a volume, and which spine is given is what says
/// whether there is one to have. A **closed** profile bounds its own section
/// and sweeps a solid along either kind of spine. An **open** one does not,
/// and there is exactly one honest way to close it: against the plane a
/// **face** spine was drawn in, whose own plane the profile's two ends must
/// reach. An open profile along a wire spine is refused, and the refusal says
/// which spine would close it.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the spine is
/// not a planar wire or face, if it carries an edge that is neither straight
/// nor circular, if the profile is not planar, if the profile's plane does not
/// contain the spine's normal or does not cut across it — a profile that leans
/// or lies along is not square to the spine — or if an open profile has no
/// spine plane to close against.
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where a corner's turn
/// would sweep the profile across the corner itself, which no revolution can
/// express.
pub fn make_evolved(
    model: &mut Model,
    spine: &Shape,
    profile: &Shape,
    tol: Tolerances,
) -> OgeomResult<Built> {
    use ogeom_algo::{make_prism, make_revolution, transformed};

    let (wire, capped_by_spine_plane) = match model.kind_of(spine)? {
        ShapeType::Face => {
            let wires = explore(model, spine, Filter::OfType(ShapeType::Wire))?;
            let Some(outer) = wires.first().cloned() else {
                ogeom_bail!(Construction, "a face with no wire has no spine to run");
            };
            (outer, true)
        }
        ShapeType::Wire => (spine.clone(), false),
        other => ogeom_bail!(
            Construction,
            "a {other:?} is not a spine; sweep along a wire or a planar face"
        ),
    };

    let stations = spine_stations(model, &wire, tol)?;
    if stations.is_empty() {
        ogeom_bail!(Construction, "a spine with no edges goes nowhere");
    }
    let normal = spine_normal(&stations, tol)?;
    let (profile_origin, profile_normal) = profile_plane(model, profile, tol)?;
    if profile_normal.dot(normal.vector()).abs() > tol.angular() {
        ogeom_bail!(
            Construction,
            "the profile's plane must contain the spine's normal, or the \
             profile is not square to the spine it travels"
        );
    }
    let start = &stations[0];
    if profile_normal.cross(start.tangent_in).magnitude() > tol.angular() {
        ogeom_bail!(
            Construction,
            "the profile's plane must cut the spine across, not run along it: \
             the profile is not square to the spine it travels"
        );
    }
    let _ = profile_origin;
    let reference = (start.from, start.tangent_in);

    // The profile as a face, which is what makes each swept piece a *solid*
    // and the assembly a union rather than a hopeful sew. An open profile is
    // closed against the spine's own plane — which is exactly what a face
    // spine offers and a wire spine does not.
    let section = profile_face(
        model,
        profile,
        profile_normal,
        capped_by_spine_plane.then(|| Plane::through(start.from, normal)),
        tol,
    )?;

    let mut pieces: Vec<Shape> = Vec::new();
    for (index, station) in stations.iter().enumerate() {
        // The corner *before* this station, so the pieces come out in the
        // order the spine runs them.
        if index > 0 {
            let previous = &stations[index - 1];
            if let Some(piece) = corner_piece(
                model,
                &section,
                reference,
                previous.to,
                previous.tangent_out,
                station.tangent_in,
                normal,
                tol,
            )? {
                pieces.push(piece);
            }
        }
        let placed = transformed(
            model,
            &section,
            station_transform(reference, station.from, station.tangent_in, normal, tol)?,
        )?
        .shape;
        pieces.push(match station.turn {
            None => make_prism(model, &placed, station.to - station.from, tol)?.shape,
            Some((axis, angle)) => make_revolution(model, &placed, axis, angle, tol)?.shape,
        });
    }
    // A closed spine turns at the join between its last edge and its first
    // just as it does anywhere else.
    let last = &stations[stations.len() - 1];
    if last.to.distance(start.from) <= tol.confusion()
        && let Some(piece) = corner_piece(
            model,
            &section,
            reference,
            last.to,
            last.tangent_out,
            start.tangent_in,
            normal,
            tol,
        )?
    {
        pieces.push(piece);
    }

    // The union, in the order the spine runs: consecutive pieces meet on the
    // *same* placed profile, which is the coincident-face case the boolean
    // resolves by identifying it rather than by probing across it.
    let mut history = History::new();
    let mut shape = pieces[0].clone();
    for piece in &pieces[1..] {
        shape = ogeom_bool::fuse(model, &shape, piece, tol)?.shape;
    }
    history.generate(spine, shape.clone());
    history.generate(profile, shape.clone());
    Ok(Built::new(shape, history))
}

/// The profile as a face.
///
/// A closed profile bounds its own area. An open one does not, and there is
/// exactly one honest way to close it: against the plane the spine was given
/// as a face *in*, which is what a face spine says to do and a wire spine has
/// no answer for. The closing segment runs between the profile's two ends, and
/// both have to be on that plane or the profile does not reach it.
fn profile_face(
    model: &mut Model,
    profile: &Shape,
    profile_normal: Vector,
    against: Option<Plane>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    if model.kind_of(profile)? == ShapeType::Face {
        return Ok(profile.clone());
    }
    if model.kind_of(profile)? != ShapeType::Wire {
        ogeom_bail!(Construction, "a profile is a wire or a face");
    }
    let mut edges = ogeom_topo::explore(model, profile, Filter::OfType(ShapeType::Edge))?;
    let closed = ogeom_algo::is_wire_closed(model, profile, tol)?;
    if !closed {
        let Some(plane) = against else {
            ogeom_bail!(
                Construction,
                "an open profile sweeps a shell, not a volume; give the spine \
                 as a planar face for its plane to close the profile against, \
                 or close the profile itself"
            );
        };
        let [(from, v0), (to, v1)] = wire_ends(model, profile, tol)?;
        for end in [from, to] {
            if plane.signed_distance_to(end).abs() > tol.confusion() * 1e2 {
                ogeom_bail!(
                    Construction,
                    "an open profile is closed against the spine face's own \
                     plane, and this one does not reach it"
                );
            }
        }
        // Built on the profile's *own* end vertices, so the closed ring is a
        // wire rather than edges that merely touch.
        let line = LineCurve::new(ogeom_math::Axis {
            location: from,
            direction: Direction::new(to - from, tol)?,
        });
        edges.push(
            make_edge_between(model, line.into(), (0.0, from.distance(to)), &v0, &v1, tol)?.shape,
        );
    }
    let ordered = ogeom_algo::order_edges(model, &edges, tol)?;
    let mut bound = ogeom_math::Aabb::EMPTY;
    for edge in &ordered {
        bound = bound.union(&ogeom_algo::shape_bounds(model, edge, tol)?);
    }
    let Some(centre) = bound.centre() else {
        ogeom_bail!(Construction, "a profile with no extent sweeps nothing");
    };
    let reach = bound.diagonal().mul_add(2.0, 1.0);
    let plane = Plane::through(centre, Direction::new(profile_normal, tol)?);
    let surface = PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?;
    Ok(make_face_with_pcurves(model, surface.into(), &[ordered], tol)?.shape)
}

/// Where an open wire begins and ends: the point, and the vertex there.
fn wire_ends(model: &Model, wire: &Shape, tol: Tolerances) -> OgeomResult<[(Point, Shape); 2]> {
    let mut counts: Vec<(Point, Shape, usize)> = Vec::new();
    for edge in explore(model, wire, Filter::OfType(ShapeType::Edge))? {
        for v in explore(model, &edge, Filter::OfType(ShapeType::Vertex))? {
            let Some(data) = model.node(&v).and_then(|n| n.data().as_vertex()) else {
                continue;
            };
            let at = v.transform(model.datums())?.apply(data.point);
            match counts
                .iter_mut()
                .find(|(p, _, _)| p.distance(at) <= tol.confusion() * 10.0)
            {
                Some((_, _, n)) => *n += 1,
                None => counts.push((at, v.clone(), 1)),
            }
        }
    }
    let free: Vec<(Point, Shape)> = counts
        .into_iter()
        .filter(|(_, _, n)| *n == 1)
        .map(|(p, v, _)| (p, v))
        .collect();
    if free.len() != 2 {
        ogeom_bail!(
            Construction,
            "an open profile has exactly two ends; this one has {}",
            free.len()
        );
    }
    let mut ends = free.into_iter();
    let (Some(a), Some(b)) = (ends.next(), ends.next()) else {
        ogeom_bail!(Construction, "the profile lost an end between checks");
    };
    Ok([a, b])
}

/// The spine, edge by edge, in the order the wire runs it.
fn spine_stations(model: &Model, wire: &Shape, tol: Tolerances) -> OgeomResult<Vec<Station>> {
    let mut out = Vec::new();
    for edge in explore(model, wire, Filter::OfType(ShapeType::Edge))? {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            ogeom_bail!(Construction, "a spine edge is not in this model");
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            ogeom_bail!(Construction, "a spine edge with no curve runs nowhere");
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            ogeom_bail!(Dangling, "curve is not in this model");
        };
        let placed = geometry.transformed(&edge.transform(model.datums())?, tol)?;
        let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
        let (t0, t1) = if reversed {
            (range.1, range.0)
        } else {
            (range.0, range.1)
        };
        let sign = if reversed { -1.0 } else { 1.0 };
        let unit = |t: f64| -> OgeomResult<Vector> {
            let d = placed.d1_at(t, tol)? * sign;
            if d.magnitude() <= tol.confusion() {
                ogeom_bail!(Construction, "a spine edge has no direction at {t}");
            }
            Ok(d / d.magnitude())
        };
        let station = match &placed {
            Curve::Line(_) => Station {
                from: placed.point_at(t0, tol)?,
                to: placed.point_at(t1, tol)?,
                tangent_in: unit(t0)?,
                tangent_out: unit(t1)?,
                turn: None,
            },
            Curve::Circle(c) => {
                let circle = c.circle();
                let swept = (range.1 - range.0).abs();
                let axis = ogeom_math::Axis {
                    location: circle.centre(),
                    direction: if reversed {
                        -circle.frame().z()
                    } else {
                        circle.frame().z()
                    },
                };
                Station {
                    from: placed.point_at(t0, tol)?,
                    to: placed.point_at(t1, tol)?,
                    tangent_in: unit(t0)?,
                    tangent_out: unit(t1)?,
                    turn: Some((axis, swept)),
                }
            }
            other => ogeom_bail!(
                Construction,
                "a spine runs on straight and circular edges; a {:?} sweeps a \
                 surface this construction does not have",
                other.kind()
            ),
        };
        out.push(station);
    }
    Ok(out)
}

/// The spine's own normal, and the check that it has one.
///
/// Taken from the first turn the spine makes — a corner or an arc — because
/// that is exact, and then measured against every station: a spine that
/// leaves its own plane has no square profile to carry, and says so here
/// rather than by producing a shape nobody asked for.
fn spine_normal(stations: &[Station], tol: Tolerances) -> OgeomResult<Direction> {
    let mut best: Option<(f64, Vector)> = None;
    let mut consider = |a: Vector, b: Vector| {
        let cross = a.cross(b);
        let magnitude = cross.magnitude();
        if magnitude > best.map_or(tol.angular(), |(m, _)| m) {
            best = Some((magnitude, cross / magnitude));
        }
    };
    for (index, station) in stations.iter().enumerate() {
        consider(station.tangent_in, station.tangent_out);
        if index + 1 < stations.len() {
            consider(station.tangent_out, stations[index + 1].tangent_in);
        }
    }
    if stations.len() > 1 {
        consider(
            stations[stations.len() - 1].tangent_out,
            stations[0].tangent_in,
        );
    }
    let Some((_, normal)) = best else {
        ogeom_bail!(
            Construction,
            "a spine that never turns has no plane of its own; give the \
             profile's own orientation a spine with at least one corner or arc"
        );
    };
    for station in stations {
        for tangent in [station.tangent_in, station.tangent_out] {
            if tangent.dot(normal).abs() > tol.angular() {
                ogeom_bail!(
                    Construction,
                    "the spine leaves its own plane; an evolved sweep runs a \
                     planar spine"
                );
            }
        }
        if let Some((axis, _)) = &station.turn
            && axis.direction.vector().cross(normal).magnitude() > tol.angular()
        {
            ogeom_bail!(
                Construction,
                "a spine arc turns about an axis off the spine's own normal"
            );
        }
    }
    Direction::new(normal, tol)
}

/// The profile's plane: a point on it and its normal.
fn profile_plane(model: &Model, profile: &Shape, tol: Tolerances) -> OgeomResult<(Point, Vector)> {
    let mut points: Vec<Point> = Vec::new();
    for edge in explore(model, profile, Filter::OfType(ShapeType::Edge))? {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            ogeom_bail!(Dangling, "curve is not in this model");
        };
        let placed = geometry.transformed(&edge.transform(model.datums())?, tol)?;
        for i in 0..=8 {
            let t = range.0 + (range.1 - range.0) * f64::from(i) / 8.0;
            points.push(placed.point_at(t, tol)?);
        }
    }
    if points.len() < 3 {
        ogeom_bail!(Construction, "a profile needs an extent to sweep");
    }
    let origin = points[0];
    // The widest cross product among the sampled offsets: the plane's normal,
    // taken where it is best conditioned rather than from the first three
    // points that happen to be there.
    let mut best: Option<(f64, Vector)> = None;
    for (i, a) in points.iter().enumerate() {
        for b in points.iter().skip(i + 1) {
            let cross = (*a - origin).cross(*b - origin);
            let magnitude = cross.magnitude();
            if magnitude > best.map_or(tol.confusion(), |(m, _)| m) {
                best = Some((magnitude, cross / magnitude));
            }
        }
    }
    let Some((_, normal)) = best else {
        ogeom_bail!(Construction, "a profile with no area has no plane");
    };
    for p in &points {
        if (*p - origin).dot(normal).abs() > tol.confusion() * 1e2 {
            ogeom_bail!(Construction, "the profile is not planar");
        }
    }
    Ok((origin, normal))
}

/// The rigid motion that carries the profile from the spine's start to a
/// station: a turn about the spine's normal, then a translation.
fn station_transform(
    reference: (Point, Vector),
    at: Point,
    tangent: Vector,
    normal: Direction,
    tol: Tolerances,
) -> OgeomResult<Transform> {
    let (origin, from) = reference;
    let n = normal.vector();
    let angle = from.cross(tangent).dot(n).atan2(from.dot(tangent));
    let turn = if angle.abs() <= tol.angular() {
        Transform::IDENTITY
    } else {
        Transform::rotation(
            ogeom_math::Axis {
                location: origin,
                direction: normal,
            },
            angle,
        )
    };
    Ok(Transform::translation(at - origin) * turn)
}

/// The wedge a corner adds: the profile turned about the corner, through
/// exactly the angle the spine turns there.
///
/// `None` where the spine does not turn — two edges meeting smoothly leave no
/// wedge to fill.
#[allow(clippy::too_many_arguments)]
fn corner_piece(
    model: &mut Model,
    profile: &Shape,
    reference: (Point, Vector),
    corner: Point,
    incoming: Vector,
    outgoing: Vector,
    normal: Direction,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let n = normal.vector();
    let angle = incoming
        .cross(outgoing)
        .dot(n)
        .atan2(incoming.dot(outgoing));
    if angle.abs() <= tol.angular() {
        return Ok(None);
    }
    let placed = ogeom_algo::transformed(
        model,
        profile,
        station_transform(reference, corner, incoming, normal, tol)?,
    )?
    .shape;
    let axis = ogeom_math::Axis {
        location: corner,
        direction: if angle > 0.0 { normal } else { -normal },
    };
    let turned = ogeom_algo::make_revolution(model, &placed, axis, angle.abs(), tol);
    match turned {
        Ok(built) => Ok(Some(built.shape)),
        Err(_) => ogeom_bail!(
            NotDone,
            "the profile straddles the spine at a corner, so turning it about \
             that corner sweeps it through itself; there is no revolution for \
             that wedge"
        ),
    }
}
