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

use ogeom_algo::{
    Built, edge_vertices, make_cone, make_cylinder, make_edge, make_edge_between,
    make_face_with_pcurves, make_solid, make_torus, make_vertex, sew,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::{
    CircleCurve, Curve, Line2d, LineCurve, PlaneSurface, SurfaceGeometry, TorusSurface,
};
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Point2, Torus, Vector};
use ogeom_topo::{EdgeRepr, Filter, Model, Shape, ShapeType, explore};

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
                     see the deferred table"
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
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the loft did not close");
    }
    let mut built = make_solid(model, std::slice::from_ref(&sewn.shells[0]))?;
    built.history.generate(bottom, built.shape.clone());
    built.history.generate(top, built.shape.clone());
    Ok(built)
}

/// A solid skinned over a grid of section samples, closed the way round.
///
/// The wall is [`ogeom_geom::fit::fit_surface_grid`]'s surface with each row's
/// first sample repeated at its end: the row fits pin their ends, so the two
/// border control columns are *equal* and the seam closes exactly, not
/// within tolerance. The border iso-curves come straight off the control
/// net — the v-borders are the fitted sections, planar whenever the
/// sections are, which is what lets the caps be planes — and every pcurve
/// is an iso line in the fitted chart, same-parameter by construction.
fn skinned_solid(
    model: &mut Model,
    rows: &[Vec<Point>],
    cap_outward: (Vector, Vector),
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
    let cap0 = cap(model, &ring0, border_v(0)?, cap_outward.0)?;
    let cap1 = cap(model, &ring1, border_v(l - 1)?, cap_outward.1)?;

    let faces = [wall, cap0, cap1];
    let sewn = sew(model, &faces, tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the skinned solid did not close");
    }
    make_solid(model, std::slice::from_ref(&sewn.shells[0]))
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
    let mut rows: Vec<Vec<Point>> = Vec::with_capacity(sections.len());
    let mut planes: Vec<Plane> = Vec::with_capacity(sections.len());
    for wire in sections {
        if model.kind_of(wire)? != ShapeType::Wire {
            ogeom_bail!(Construction, "a loft section is a wire");
        }
        if !ogeom_algo::is_wire_closed(model, wire, tol)? {
            ogeom_bail!(Construction, "a loft section must be closed");
        }
        let Some(plane) = ogeom_algo::find_plane(model, wire, tol)? else {
            ogeom_bail!(Construction, "a loft section must be planar");
        };
        planes.push(plane);
        rows.push(sample_wire(model, wire, AROUND, tol)?);
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

/// Sample a closed wire at `count` matched arc-length fractions.
fn sample_wire(
    model: &Model,
    wire: &Shape,
    count: usize,
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
    // Rotation-minimizing frames by double reflection.
    let mut stations: Vec<(Point, Vector)> = Vec::with_capacity(STATIONS);
    for i in 0..STATIONS {
        #[allow(clippy::cast_precision_loss)]
        let t = range.0 + (range.1 - range.0) * (i as f64) / ((STATIONS - 1) as f64);
        let p = curve.point_at(t, tol)?;
        let d = curve.d1_at(t, tol)?;
        let m = d.magnitude();
        if m <= tol.confusion() {
            ogeom_bail!(Construction, "the spine is degenerate at {t}");
        }
        stations.push((p, d / m));
    }
    let mut normals: Vec<Vector> = Vec::with_capacity(STATIONS);
    {
        let t0 = stations[0].1;
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
        for i in 1..STATIONS {
            let (p0, t0) = stations[i - 1];
            let (p1, t1) = stations[i];
            let n = normals[i - 1];
            // Double reflection: reflect in the chord plane, then in the
            // plane bisecting the tangents.
            let v1 = p1 - p0;
            let c1 = v1.dot(v1);
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
    }
    let mut rows: Vec<Vec<Point>> = Vec::with_capacity(STATIONS);
    for (i, (p, t)) in stations.iter().enumerate() {
        let x = normals[i];
        let y = t.cross(x);
        let mut row = Vec::with_capacity(AROUND);
        for a in 0..AROUND {
            #[allow(clippy::cast_precision_loss)]
            let ang = core::f64::consts::TAU * (a as f64) / (AROUND as f64);
            row.push(*p + (x * ang.cos() + y * ang.sin()) * radius);
        }
        rows.push(row);
    }
    let mut built = skinned_solid(
        model,
        &rows,
        (-stations[0].1, stations[STATIONS - 1].1),
        tolerance,
        tol,
    )?;
    built.history.generate(spine, built.shape.clone());
    Ok(built)
}
