//! Constant-radius fillets, where the rolling ball's envelope is closed-form.
//!
//! A ball of constant radius rolling along a straight edge between two planar
//! faces stays tangent to both, and its envelope is exactly a cylinder; the
//! same ball rolling around the rim where a cylindrical wall meets its
//! perpendicular cap traces a torus. These are the blends whose *new* surface
//! costs nothing, and both are built the same way: the chamfer's wedge with
//! the bevel exchanged for the envelope — legs running along the faces to the
//! tangency lines, and the envelope between them, every pcurve exact. The
//! boolean then does what it did for the chamfer — the legs melt into the
//! solid's own faces by same-domain resolution, and the envelope stays as
//! the blend.

use crate::support::{
    Seat, apply_wedge, edge_curve, face_from_edges, planar_face, planar_seat, segment_between,
};
use ogeom_algo::{Built, make_edge, make_edge_between, make_revolution_band, make_vertex};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    CircleCurve, Curve, CylinderSurface, PlaneSurface, SurfaceGeometry, TorusSurface,
};
use ogeom_math::{Circle, Cylinder, Direction, Frame, Plane, Point, Torus, Vector};
use ogeom_topo::{Filter, Model, NodeData, Orientation, Shape, ShapeType, explore};

/// Round a convex edge of a solid with a constant-radius blend tangent to
/// both of its faces.
///
/// Two seats are recognised, and each has a closed-form envelope: a straight
/// edge between two planar faces, where the rolling ball traces a cylinder,
/// and a circular rim where a cylindrical wall meets a perpendicular planar
/// cap, where it traces a torus. In both, the result is the boolean
/// difference with a wedge whose curved face is that envelope, so the history
/// reads as a cut: the faces are modified into their trimmed pieces, and the
/// edge's neighbourhood gains the blend face.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the edge is
/// neither of the seats above, is concave, or `radius` is not a usable
/// length.
pub fn fillet_edge(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(Construction, "a fillet of radius {radius} rounds nothing");
    }
    let (curve, _) = edge_curve(model, edge)?;
    match curve {
        Curve::Line(_) => planar_fillet(model, solid, edge, radius, tol),
        Curve::Circle(c) => revolved_fillet(model, solid, edge, &c, radius, tol),
        _ => ogeom_bail!(
            Construction,
            "filleting an edge that is neither straight nor circular needs \
             the marching blend machinery"
        ),
    }
}

/// Round a straight convex edge with a blend whose radius runs linearly from
/// `start_radius` at the edge's start to `end_radius` at its end.
///
/// For a linear law on a straight edge between planes the rolling ball's
/// envelope is *exactly* a rational B-spline surface — degree one along the
/// edge, a rational quadratic arc across it, the control net affine in the
/// radius — so nothing here is fitted. The tangency lines are straight, the
/// legs stay planar, and the wedge subtracts through the boolean like its
/// constant-radius siblings.
///
/// # Errors
///
/// As [`fillet_edge`], for the straight planar seat only.
pub fn fillet_edge_variable(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    start_radius: f64,
    end_radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    for radius in [start_radius, end_radius] {
        if !radius.is_finite() || radius <= tol.confusion() {
            ogeom_bail!(Construction, "a fillet of radius {radius} rounds nothing");
        }
    }
    if (start_radius - end_radius).abs() <= tol.confusion() {
        return fillet_edge(model, solid, edge, start_radius, tol);
    }
    let seat = planar_seat(model, solid, edge, tol)?;
    let (first, second) = if seat
        .along
        .dot(seat.normals[0].cross(seat.normals[1]))
        .is_sign_positive()
    {
        (0, 1)
    } else {
        (1, 0)
    };
    let sign = if seat.convex { 1.0 } else { -1.0 };
    let n1 = seat.normals[first] * sign;
    let n2 = seat.normals[second] * sign;
    let a = seat.leg(first, tol)? * sign;
    let b = seat.leg(second, tol)? * sign;
    let bisector = {
        let u = a + b;
        let m = u.magnitude();
        if m <= tol.angular() {
            ogeom_bail!(Construction, "the faces meet too sharply to seat a fillet");
        }
        u / m
    };
    let depth = -n1.dot(bisector);
    if depth <= tol.angular() {
        ogeom_bail!(Construction, "the faces meet too sharply to seat a fillet");
    }
    let sweep = seat.along.dot(n1.cross(n2)).atan2(n1.dot(n2)) * sign;
    if sweep <= tol.angular() {
        ogeom_bail!(Construction, "the faces are parallel; there is no corner");
    }

    // The sections at the two ends: everything else is affine between them.
    // On a concave edge every direction above is already mirrored — the
    // arithmetic below cannot tell which case it serves.
    let radii = [start_radius, end_radius];
    let apex = [seat.start, seat.end];
    let centre = [
        apex[0] + bisector * (radii[0] / depth),
        apex[1] + bisector * (radii[1] / depth),
    ];
    let contact_a = [centre[0] + n1 * radii[0], centre[1] + n1 * radii[1]];
    let contact_b = [centre[0] + n2 * radii[0], centre[1] + n2 * radii[1]];

    // The unit arc from n1 to n2 as a rational quadratic: the middle control
    // point is where the end tangents meet, its weight the half-angle cosine.
    let half_cos = (sweep / 2.0).cos();
    let arc_mid = (n1 + n2) / (1.0 + n1.dot(n2));
    let weights = [1.0, half_cos, 1.0];

    // The blend: degree one by rational quadratic, the control net
    // c_i + r_i * V_j, weights carried across the rows unchanged. Every
    // section of the surface is that section's exact blend arc.
    let u_knots = ogeom_math::KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1)?;
    let v_knots = ogeom_math::KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2)?;
    let mut net: Vec<ogeom_math::Weighted<Point>> = Vec::with_capacity(6);
    for i in 0..2 {
        let arms = [n1, arc_mid, n2];
        for (j, arm) in arms.iter().enumerate() {
            net.push(ogeom_math::Weighted::new(
                centre[i] + *arm * radii[i],
                weights[j],
                tol,
            )?);
        }
    }
    let grid = ogeom_math::ControlGrid::new(net, 2, 3)?;
    let surface: SurfaceGeometry =
        ogeom_geom::BSplineSurface::rational(u_knots, v_knots, grid)?.into();
    let surface_id = model.geometry_mut().add_surface(surface.clone());

    // Shared vertices at the four tangency corners.
    let va = [
        make_vertex(model, contact_a[0]).shape,
        make_vertex(model, contact_a[1]).shape,
    ];
    let vb = [
        make_vertex(model, contact_b[0]).shape,
        make_vertex(model, contact_b[1]).shape,
    ];

    // The tangency lines, straight because the law is linear, with
    // degree-one pcurves mapping their arc length onto the chart's rows.
    let rail = |model: &mut Model,
                from: (&Shape, Point),
                to: (&Shape, Point),
                row: f64|
     -> OgeomResult<Shape> {
        let line = ogeom_geom::LineCurve::segment(from.1, to.1, tol)?;
        let curve = Curve::Line(line);
        let domain = ogeom_geom::Curve3d::domain(&curve);
        let built = make_edge_between(model, curve, domain, from.0, to.0, tol)?.shape;
        let knots = ogeom_math::KnotVector::new(vec![domain.0, domain.0, domain.1, domain.1], 1)?;
        let pcurve = ogeom_geom::BSpline2d::new(
            knots,
            vec![
                ogeom_math::Point2::new(0.0, row),
                ogeom_math::Point2::new(1.0, row),
            ],
            tol,
        )?;
        ogeom_algo::attach_pcurve(
            model,
            &built,
            pcurve.into(),
            surface_id,
            ogeom_topo::Location::identity(),
            domain,
        )?;
        Ok(built)
    };
    let rail_a = rail(model, (&va[0], contact_a[0]), (&va[1], contact_a[1]), 0.0)?;
    let rail_b = rail(model, (&vb[0], contact_b[0]), (&vb[1], contact_b[1]), 1.0)?;

    // The end arcs as the surface's own rational sections, so the cap edges
    // and the blend's chart speak the same parameter.
    let arc_edge = |model: &mut Model, i: usize| -> OgeomResult<Shape> {
        let knots = ogeom_math::KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2)?;
        let control = vec![
            ogeom_math::Weighted::new(contact_a[i], 1.0, tol)?,
            ogeom_math::Weighted::new(centre[i] + arc_mid * radii[i], half_cos, tol)?,
            ogeom_math::Weighted::new(contact_b[i], 1.0, tol)?,
        ];
        let curve = Curve::BSpline(ogeom_geom::BSplineCurve::rational(knots, control)?);
        let built = make_edge_between(model, curve, (0.0, 1.0), &va[i], &vb[i], tol)?.shape;
        #[allow(clippy::cast_precision_loss)]
        let column = ogeom_geom::Line2d::over(
            ogeom_math::Axis2::new(
                ogeom_math::Point2::new(i as f64, 0.0),
                ogeom_math::Direction2::Y,
            ),
            0.0,
            1.0,
        )?;
        ogeom_algo::attach_pcurve(
            model,
            &built,
            column.into(),
            surface_id,
            ogeom_topo::Location::identity(),
            (0.0, 1.0),
        )?;
        Ok(built)
    };
    let arc0 = arc_edge(model, 0)?;
    let arc1 = arc_edge(model, 1)?;

    // The blend face on the registered surface, oriented so its outward side
    // leaves the wedge — decided by measurement at the middle rather than by
    // convention.
    let blend = {
        let wire = ogeom_algo::make_wire(
            model,
            &[
                rail_a.clone(),
                arc1.clone(),
                rail_b.reversed(),
                arc0.reversed(),
            ],
            tol,
        )?
        .shape;
        let face =
            ogeom_algo::make_face_on(model, surface_id, std::slice::from_ref(&wire), tol)?.shape;
        use ogeom_geom::Surface as _;
        let s_mid = surface.point_at(0.5, 0.5, tol)?;
        let (du, dv) = surface.d1_at(0.5, 0.5, tol)?;
        let apex_mid = apex[0].midpoint(apex[1]);
        let outward = s_mid - apex_mid.midpoint(s_mid);
        if du.cross(dv).dot(outward) >= 0.0 {
            face
        } else {
            face.reversed()
        }
    };

    // Legs along the two faces, and caps closed by the end arcs. The caps go
    // through the generic builder: the rational arc projects into a plane
    // exactly, control point by control point.
    let leg_a = planar_face(
        model,
        &[apex[0], contact_a[0], contact_a[1], apex[1]],
        n1 * sign,
        tol,
    )?;
    let leg_b = planar_face(
        model,
        &[apex[0], contact_b[0], contact_b[1], apex[1]],
        n2 * sign,
        tol,
    )?;
    let cap = |model: &mut Model, i: usize, outward: Vector| -> OgeomResult<Built> {
        let plane = Plane::through(apex[i], Direction::new(outward, tol)?);
        let reach = (apex[i].distance(contact_a[i]) + apex[i].distance(contact_b[i]) + 1.0) * 2.0;
        let cap_surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let apex_v = make_vertex(model, apex[i]).shape;
        let edges = vec![
            crate::support::segment_between(
                model,
                (&apex_v, apex[i]),
                (&va[i], contact_a[i]),
                tol,
            )?,
            if i == 0 { arc0.clone() } else { arc1.clone() },
            crate::support::segment_between(
                model,
                (&vb[i], contact_b[i]),
                (&apex_v, apex[i]),
                tol,
            )?,
        ];
        ogeom_algo::make_face_with_pcurves(model, cap_surface, &[edges], tol)
    };
    let cap0 = cap(model, 0, -seat.along)?.shape;
    let cap1 = cap(model, 1, seat.along)?.shape;

    let faces = [leg_a, leg_b, cap0, cap1, blend];
    apply_wedge(model, solid, Some(edge), &faces, !seat.convex, tol)
}

/// The straight-edge fillet: the rolling ball's envelope is a cylinder.
fn planar_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let seat = planar_seat(model, solid, edge, tol)?;
    seated_fillet(model, solid, &seat, radius, Some(edge), tol)
}

/// The cylindrical blend of a seat, whatever found the seat.
///
/// A blend's construction cares about the seat — where the ball rolls and
/// which way the two planes face — and not at all about whether an edge of
/// the solid runs along it. An edge gives one; two faces that share nothing
/// give the same one through their planes' own intersection, and everything
/// from here down is common to both.
pub(crate) fn seated_fillet(
    model: &mut Model,
    solid: &Shape,
    seat: &Seat,
    radius: f64,
    edge: Option<&Shape>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    // Order the two faces so the blend arc sweeps positively about the edge
    // direction — the cylinder's parameterization and every arc below then
    // run from the first face's tangency line to the second's.
    let (first, second) = if seat
        .along
        .dot(seat.normals[0].cross(seat.normals[1]))
        .is_sign_positive()
    {
        (0, 1)
    } else {
        (1, 0)
    };
    let n1 = seat.normals[first];
    let n2 = seat.normals[second];
    // On a concave edge the legs mirror into the open dihedral, the ball
    // rolls there instead of in the material, and the contacts sit on the
    // faces' other side of the centre.
    let sign = if seat.convex { 1.0 } else { -1.0 };
    let a = seat.leg(first, tol)? * sign;
    let b = seat.leg(second, tol)? * sign;

    // The rolling ball's centre line: along the bisector, at the distance
    // where the ball touches both planes.
    let bisector = {
        let u = a + b;
        let m = u.magnitude();
        if m <= tol.angular() {
            ogeom_bail!(Construction, "the faces meet too sharply to seat a fillet");
        }
        u / m
    };
    let depth = (n1.dot(bisector)).abs();
    if depth <= tol.angular() {
        ogeom_bail!(Construction, "the faces meet too sharply to seat a fillet");
    }
    let centre = seat.start + bisector * (radius / depth);
    let sweep = seat.along.dot(n1.cross(n2)).atan2(n1.dot(n2));
    if sweep <= tol.angular() {
        ogeom_bail!(Construction, "the faces are parallel; there is no corner");
    }

    let travel = seat.end - seat.start;
    let length = travel.magnitude();
    let apex0 = seat.start;
    let apex1 = seat.end;
    let contact_a0 = centre + n1 * (radius * sign);
    let contact_b0 = centre + n2 * (radius * sign);
    let contact_a1 = contact_a0 + travel;
    let contact_b1 = contact_b0 + travel;

    let along_dir = Direction::new(seat.along, tol)?;
    let radial_dir = Direction::new(n1 * sign, tol)?;

    // An arc of the blend circle at height `h` along the edge, running from
    // the first face's contact vertex to the second's. Its frame is the
    // cylinder's frame translated up the axis, which is what makes its pcurve
    // on the cylinder a straight line in `(u, v)` with the arc's own
    // parameter.
    let arc = |model: &mut Model, h: f64, from: &Shape, to: &Shape| -> OgeomResult<Shape> {
        let frame = Frame::new(centre + seat.along * h, along_dir, radial_dir, tol)?;
        let circle = Circle::new(frame, radius, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        Ok(make_edge_between(model, curve, (0.0, sweep), from, to, tol)?.shape)
    };

    // The legs, coplanar with the solid's own faces: aligned with them on a
    // convex edge, where the melt keeps one copy, opposed on a concave one,
    // where the fuse cancels both.
    let leg_a = planar_face(
        model,
        &[apex0, contact_a0, contact_a1, apex1],
        n1 * sign,
        tol,
    )?;
    let leg_b = planar_face(
        model,
        &[apex0, contact_b0, contact_b1, apex1],
        n2 * sign,
        tol,
    )?;

    // The caps: two segments to the tangency points, closed by the arc.
    let cap = |model: &mut Model,
               apex: Point,
               ca: Point,
               cb: Point,
               h: f64,
               outward: Vector|
     -> OgeomResult<Shape> {
        let plane = Plane::through(apex, Direction::new(outward, tol)?);
        let reach = (apex.distance(ca).max(apex.distance(cb)) + radius).max(1.0) * 2.0;
        let surface = PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?;
        let apex_v = make_vertex(model, apex).shape;
        let ca_v = make_vertex(model, ca).shape;
        let cb_v = make_vertex(model, cb).shape;
        let edges = vec![
            segment_between(model, (&apex_v, apex), (&ca_v, ca), tol)?,
            arc(model, h, &ca_v, &cb_v)?,
            segment_between(model, (&cb_v, cb), (&apex_v, apex), tol)?,
        ];
        face_from_edges(model, surface.into(), &edges, tol)
    };
    let cap0 = cap(model, apex0, contact_a0, contact_b0, 0.0, -seat.along)?;
    let cap1 = cap(model, apex1, contact_a1, contact_b1, length, seat.along)?;

    // The blend face itself. A cylinder's natural normal points away from its
    // axis — into this wedge, whose material lies between the cylinder and
    // the apex — so the face enters the shell reversed.
    let blend = {
        let frame = Frame::new(centre, along_dir, radial_dir, tol)?;
        let surface = CylinderSurface::new(Cylinder::new(frame, radius, tol)?, (0.0, length))?;
        let a0_v = make_vertex(model, contact_a0).shape;
        let b0_v = make_vertex(model, contact_b0).shape;
        let a1_v = make_vertex(model, contact_a1).shape;
        let b1_v = make_vertex(model, contact_b1).shape;
        let edges = vec![
            arc(model, 0.0, &a0_v, &b0_v)?,
            segment_between(model, (&b0_v, contact_b0), (&b1_v, contact_b1), tol)?,
            arc(model, length, &a1_v, &b1_v)?.reversed(),
            segment_between(model, (&a0_v, contact_a0), (&a1_v, contact_a1), tol)?.reversed(),
        ];
        face_from_edges(model, surface.into(), &edges, tol)?
    };

    let faces = [leg_a, leg_b, cap0, cap1, blend.reversed()];
    apply_wedge(model, solid, edge, &faces, !seat.convex, tol)
}

/// The circular-rim fillet: a ball rolling around a rim where a cylindrical
/// wall meets a perpendicular planar cap traces a torus.
///
/// Four seats, one parameterization. With `sigma` the wall's outward radial
/// sign and `tau` telling whether the wall extends away from the cap's
/// outward side, the tube's centre circle sits at radius `R − sigma tau r`,
/// lifted `tau r` against the cap's normal — and `tau` alone decides whether
/// the wedge subtracts (the external rim and the hole's rim, both convex) or
/// fuses (the boss base and the blind hole's floor, both concave). The wedge
/// is always the same three revolved faces: a band of the wall, an annulus
/// of the cap, and the quarter-tube between them.
fn revolved_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    rim: &CircleCurve,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let rim_circle = rim.circle();
    let rim_centre = rim_circle.centre();
    let rim_radius = rim_circle.radius();
    let rim_axis = rim_circle.frame().z().vector();

    // The two faces at the rim: exactly one plane and one coaxial cylinder.
    let mut cap_normal: Option<Vector> = None;
    let mut wall: Option<(Shape, f64)> = None;
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let touches = explore(model, &face, Filter::OfType(ShapeType::Edge))?
            .iter()
            .any(|e| e.node() == edge.node());
        if !touches {
            continue;
        }
        let Some(node) = model.node(&face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let placement = face.transform(model.datums())?;
        let reversed = face.orientation() == Orientation::Reversed;
        match model.geometry().surface(data.surface) {
            Some(SurfaceGeometry::Plane(p)) => {
                let mut normal = placement.apply_vector(p.plane().normal().vector());
                if reversed {
                    normal = -normal;
                }
                if normal.cross(rim_axis).magnitude() > tol.angular()
                    || p.plane()
                        .distance_to(placement.inverse()?.apply(rim_centre))
                        > tol.confusion()
                {
                    ogeom_bail!(
                        Construction,
                        "the rim's cap is not the perpendicular plane through \
                         it; that seat needs the marching blend machinery"
                    );
                }
                cap_normal = Some(normal);
            }
            Some(SurfaceGeometry::Cylinder(c)) => {
                let cyl = c.cylinder();
                let axis_point = placement.apply(cyl.frame().origin());
                let axis_z = placement.apply_vector(cyl.frame().z().vector());
                let off_axis = {
                    let to_rim = rim_centre - axis_point;
                    (to_rim - axis_z * to_rim.dot(axis_z)).magnitude()
                };
                if axis_z.cross(rim_axis).magnitude() > tol.angular()
                    || off_axis > tol.confusion()
                    || (cyl.radius() - rim_radius).abs() > tol.confusion()
                {
                    ogeom_bail!(
                        Construction,
                        "the rim's wall is not the coaxial cylinder through \
                         it; that seat needs the marching blend machinery"
                    );
                }
                let sigma = if reversed { -1.0 } else { 1.0 };
                wall = Some((face.clone(), sigma));
            }
            Some(_) => ogeom_bail!(
                Construction,
                "the rim meets a face that is neither plane nor cylinder; \
                 that seat needs the marching blend machinery"
            ),
            None => ogeom_bail!(Dangling, "face refers to a surface not in this model"),
        }
    }
    let (Some(up_raw), Some((wall_face, sigma))) = (cap_normal, wall) else {
        ogeom_bail!(
            Construction,
            "a rim fillet needs the edge shared by one planar cap and one \
             cylindrical wall"
        );
    };
    let up = up_raw / up_raw.magnitude();

    // Which side of the cap the wall extends: read from the wall band's other
    // ring. `tau` positive means away from the cap's outward side — the rim
    // configurations — and negative means alongside it, the concave seats.
    let tau = {
        let mut side = None;
        for e in explore(model, &wall_face, Filter::OfType(ShapeType::Edge))? {
            if e.node() == edge.node() {
                continue;
            }
            // Any circular edge of the wall at another height says which
            // side the wall extends — whether the ring survives as one
            // closed edge or as the arcs a boolean rebuilt it into.
            let Ok((curve, _)) = edge_curve(model, &e) else {
                continue;
            };
            let Curve::Circle(c) = curve else {
                continue;
            };
            let lean = (c.circle().centre() - rim_centre).dot(up);
            if lean.abs() > tol.confusion() * 10.0 {
                side = Some(lean);
                break;
            }
        }
        let Some(lean) = side else {
            ogeom_bail!(
                Construction,
                "the rim's wall has no far ring to read its side from; the \
                 partial seat needs the marching blend machinery"
            );
        };
        if lean < 0.0 { 1.0 } else { -1.0 }
    };
    let additive = tau < 0.0;

    let tube_rho = sigma.mul_add(-(tau * radius), rim_radius);
    if tube_rho <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a fillet of radius {radius} swallows the axis of a rim of \
             radius {rim_radius}"
        );
    }
    let tube_level = rim_centre - up * (tau * radius);

    let up_dir = Direction::new(up, tol)?;
    let x_ref = rim_circle.frame().x();
    let frame_at = |origin: Point| Frame::new(origin, up_dir, x_ref, tol);

    let ring = |model: &mut Model, origin: Point, r: f64| -> OgeomResult<Shape> {
        let circle = Circle::new(frame_at(origin)?, r, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        let domain = ogeom_geom::Curve3d::domain(&curve);
        Ok(make_edge(model, curve, domain, tol)?.shape)
    };
    let apex_ring = ring(model, rim_centre, rim_radius)?;
    let wall_ring = ring(model, tube_level, rim_radius)?;
    let cap_ring = ring(model, rim_centre, tube_rho)?;

    // The band of the wall between the rim and the tangency parallel —
    // coincident with the solid's own wall, aligned when subtracting and
    // opposed when fusing, which is exactly what the melt needs either way.
    let wall_band = {
        let origin = rim_centre - up * (radius + 1.0);
        let surface: SurfaceGeometry = CylinderSurface::new(
            Cylinder::new(frame_at(origin)?, rim_radius, tol)?,
            (0.0, 2.0 * radius + 2.0),
        )?
        .into();
        let band = make_revolution_band(model, &surface, &wall_ring, &apex_ring, tol)?;
        if sigma * tau < 0.0 {
            band.reversed()
        } else {
            band
        }
    };

    // The blend: the quarter-tube between the tangencies, its natural normal
    // away from the tube's centre — into the wedge — so always reversed.
    let blend_band = {
        let surface: SurfaceGeometry =
            TorusSurface::new(Torus::new(frame_at(tube_level)?, tube_rho, radius, tol)?).into();
        make_revolution_band(model, &surface, &wall_ring, &cap_ring, tol)?.reversed()
    };

    // The annulus of the cap between the rim and the tangency circle.
    let annulus = {
        let plane = Plane::through(rim_centre, Direction::new(up * tau, tol)?);
        let reach = (rim_radius + radius) * 2.0;
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let outer = ogeom_algo::make_wire(model, std::slice::from_ref(&apex_ring), tol)?.shape;
        let inner = ogeom_algo::make_wire(model, std::slice::from_ref(&cap_ring), tol)?.shape;
        let face = ogeom_algo::make_face(model, surface.clone(), &[outer, inner], tol)?.shape;
        let surface_id = {
            let Some(node) = model.node(&face) else {
                ogeom_bail!(Dangling, "the face just built is not in this model");
            };
            let NodeData::Face(data) = node.data() else {
                ogeom_bail!(Construction, "the face holds no face data");
            };
            data.surface
        };
        for pedge in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
            let (curve, prange) = edge_curve(model, &pedge)?;
            let Some(pcurve) = ogeom_intersect::exact_pcurve_of(&curve, &surface, tol) else {
                ogeom_bail!(
                    Construction,
                    "an annulus edge has no closed-form pcurve on its plane"
                );
            };
            ogeom_algo::attach_pcurve(
                model,
                &pedge,
                pcurve,
                surface_id,
                ogeom_topo::Location::identity(),
                prange,
            )?;
        }
        face
    };

    let faces = [wall_band, annulus, blend_band];
    apply_wedge(model, solid, Some(edge), &faces, additive, tol)
}
