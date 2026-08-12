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
    Seat, apply_wedge, edge_curve, face_from_edges, planar_face, planar_seat, revolved_flanks,
    revolved_seat, segment_between,
};
use ogeom_algo::{Built, make_edge_between, make_revolution_band, make_vertex};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    CircleCurve, Curve, CylinderSurface, PlaneSurface, SurfaceGeometry, TorusSurface,
};
use ogeom_math::{Circle, Cylinder, Direction, Frame, Plane, Point, Torus, Vector};
use ogeom_topo::{Model, Shape};

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
    let (curve, crange) = edge_curve(model, edge, tol)?;
    let closed = ogeom_algo::edge_vertices(model, edge)?.is_some_and(|(a, b)| a.is_same(&b));
    match curve {
        Curve::Line(_) => planar_fillet(model, solid, edge, radius, tol),
        Curve::Circle(c) if closed => revolved_fillet(model, solid, edge, &c, radius, tol),
        Curve::Circle(c) => {
            // An open arc is two different seats: a piece a boolean split
            // off a full rim, whose blend must run the whole turn, or a rim
            // that genuinely stops — a stadium's rounded end. The wall
            // itself answers: probe it just below the arc's complement.
            let seat = revolved_seat(model, solid, edge, &c, tol)?;
            let frame = seat.frame_at(seat.centre, tol)?;
            let angle_of = |t: f64| -> OgeomResult<f64> {
                use ogeom_geom::Curve3d as _;
                let local = frame.to_local(c.point_at(t, tol)?);
                Ok(local.y.atan2(local.x))
            };
            let mid = f64::midpoint(
                angle_of(crange.0)?,
                angle_of(crange.0)?
                    + (angle_of(crange.1)? - angle_of(crange.0)?)
                        .rem_euclid(core::f64::consts::TAU),
            );
            let complement = mid + core::f64::consts::PI;
            let dir = frame.x().vector() * complement.cos() + frame.y().vector() * complement.sin();
            let eps = seat.radius * 1e-2;
            let probe = seat.centre + dir * seat.radius - seat.up * (seat.tau * eps);
            let deflection = ogeom_mesh::Deflection {
                chord: eps * 0.1,
                ..ogeom_mesh::Deflection::default()
            };
            let whole_turn =
                ogeom_algo::classify_on_face(model, &seat.wall_face, probe, deflection, tol)?
                    == ogeom_algo::Containment::In;
            if whole_turn {
                revolved_fillet(model, solid, edge, &c, radius, tol)
            } else {
                revolved_arc_fillet(model, solid, edge, &c, radius, tol)
            }
        }
        _ => crate::marched::marched_fillet(model, solid, edge, radius, tol),
    }
}

/// Round a chain of edges in one call, each blend applied in turn and every
/// later edge resolved through the history the earlier blends wrote.
///
/// At a tangent junction the neighbouring wedges' end caps stand in one
/// plane with one cross-section, and the boolean's same-domain resolution
/// melts them — the blends join without a seam face between them. Junctions
/// that are not tangent leave the wedges' caps standing, which is the honest
/// picture of a corner no single ball rolls around.
///
/// # Errors
///
/// As [`fillet_edge`] per edge, and additionally if an earlier blend
/// consumed or split a later edge — a chain whose members interfere is
/// refused rather than guessed at.
pub fn fillet_edges(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if edges.is_empty() {
        ogeom_bail!(Construction, "a chain of no edges rounds nothing");
    }
    let mut built: Option<Built> = None;
    for edge in edges {
        let (current, target) = match &built {
            None => (solid.clone(), edge.clone()),
            Some(b) => {
                let traced = b.history.trace(edge);
                match traced {
                    [] => ogeom_bail!(
                        Construction,
                        "an earlier blend in the chain consumed this edge; \
                         the chain's members interfere"
                    ),
                    [one] => {
                        let found = refind_edge(model, &b.shape, one, tol)?;
                        (b.shape.clone(), found)
                    }
                    _ => ogeom_bail!(
                        Construction,
                        "an earlier blend in the chain split this edge into \
                         {} pieces; the chain's members interfere",
                        traced.len()
                    ),
                }
            }
        };
        let mut step = fillet_edge(model, &current, &target, radius, tol)?;
        // The blend consumed the re-found stand-in; the caller's edge is
        // the same fact under its original name.
        if !target.is_same(edge) {
            step.history.delete(edge);
        }
        built = Some(match built {
            None => step,
            Some(prev) => Built {
                shape: step.shape.clone(),
                history: prev.history.then(&step.history),
            },
        });
    }
    Ok(built.unwrap_or_else(|| unreachable!()))
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
    let seat = revolved_seat(model, solid, edge, rim, tol)?;

    let tube_rho = seat.sigma.mul_add(-(seat.tau * radius), seat.radius);
    if tube_rho <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a fillet of radius {radius} swallows the axis of a rim of \
             radius {}",
            seat.radius
        );
    }
    let tube_level = seat.centre - seat.up * (seat.tau * radius);
    let flanks = revolved_flanks(model, &seat, radius, tube_rho, tol)?;

    // The blend: the quarter-tube between the tangencies, its natural normal
    // away from the tube's centre — into the wedge — so always reversed.
    let blend_band = {
        let surface: SurfaceGeometry = TorusSurface::new(Torus::new(
            seat.frame_at(tube_level, tol)?,
            tube_rho,
            radius,
            tol,
        )?)
        .into();
        make_revolution_band(model, &surface, &flanks.wall_ring, &flanks.cap_ring, tol)?.reversed()
    };

    let faces = [flanks.wall_band, flanks.annulus, blend_band];
    apply_wedge(model, solid, Some(edge), &faces, seat.additive(), tol)
}

/// The edge of `solid` standing where `edge` stood.
///
/// A boolean rebuilds every face it splits with fresh edges and its history
/// speaks of faces, not of them — so a chain's later edge is re-found by
/// geometry: the solid's edge whose own samples all lie on the sought
/// edge's curve. A trimmed survivor qualifies; that is exactly the part of
/// the edge still there to blend.
fn refind_edge(model: &Model, solid: &Shape, edge: &Shape, tol: Tolerances) -> OgeomResult<Shape> {
    use ogeom_geom::Curve3d as _;
    use ogeom_topo::ShapeType;
    let (curve, _) = edge_curve(model, edge, tol)?;
    let mut matches: Vec<Shape> = Vec::new();
    for candidate in ogeom_topo::explore_unique(model, solid, ShapeType::Edge)? {
        if candidate.is_same(edge) {
            return Ok(candidate);
        }
        let Ok((c_curve, c_range)) = edge_curve(model, &candidate, tol) else {
            continue;
        };
        let mut all_on = true;
        for i in 0..3 {
            let t = c_range.0 + (c_range.1 - c_range.0) * f64::from(i) / 2.0;
            let p = c_curve.point_at(t, tol)?;
            let projected = ogeom_algo::project_on_curve(&curve, p, 64, tol)?;
            if projected.distance > tol.confusion() * 100.0 {
                all_on = false;
                break;
            }
        }
        if all_on && !matches.iter().any(|m| m.is_same(&candidate)) {
            matches.push(candidate);
        }
    }
    match matches.as_slice() {
        [one] => Ok(one.clone()),
        [] => ogeom_bail!(
            Construction,
            "an earlier blend in the chain consumed this edge; the chain's \
             members interfere"
        ),
        _ => ogeom_bail!(
            Construction,
            "an earlier blend in the chain split this edge into {} pieces; \
             the chain's members interfere",
            matches.len()
        ),
    }
}

/// The rim fillet over an arc of the rim rather than its whole turn.
///
/// The same four-seat torus geometry as [`revolved_fillet`], restricted to
/// the arc's angular window and closed by two planar end caps in the
/// meridian half-planes at its ends. Every face is a chart rectangle or a
/// planar triangle, every edge an arc or a segment, every pcurve exact —
/// which is what lets a chain's neighbouring wedges melt at their shared
/// caps.
fn revolved_arc_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    rim: &CircleCurve,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    use ogeom_geom::Curve3d as _;
    let seat = revolved_seat(model, solid, edge, rim, tol)?;

    let tube_rho = seat.sigma.mul_add(-(seat.tau * radius), seat.radius);
    if tube_rho <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a fillet of radius {radius} swallows the axis of a rim of \
             radius {}",
            seat.radius
        );
    }
    let tube_level = seat.centre - seat.up * (seat.tau * radius);

    // The arc's window in the seat's own parameterization: the edge's stored
    // circle may be framed anyhow, so the ends are read back through the
    // seat frame and the window is the one that holds the arc's midpoint.
    let seat_frame = seat.frame_at(seat.centre, tol)?;
    let (curve, crange) = edge_curve(model, edge, tol)?;
    let angle_of = |t: f64| -> OgeomResult<f64> {
        let local = seat_frame.to_local(curve.point_at(t, tol)?);
        Ok(local.y.atan2(local.x))
    };
    let tau_turn = core::f64::consts::TAU;
    let raw0 = angle_of(crange.0)?;
    let raw1 = angle_of(crange.1)?;
    let mid = angle_of(f64::midpoint(crange.0, crange.1))?;
    let (theta0, theta1) = {
        let span = (raw1 - raw0).rem_euclid(tau_turn);
        let inside = (mid - raw0).rem_euclid(tau_turn) <= span;
        if inside {
            (raw0, raw0 + span)
        } else {
            (raw1, raw1 + (tau_turn - span))
        }
    };
    if theta1 - theta0 <= tol.angular() {
        ogeom_bail!(Construction, "the rim arc has no length to blend along");
    }
    let sweep = theta1 - theta0;

    // Corner points at each end: A on the rim, W down the wall, C in on the
    // cap; and the frames the arcs run on.
    let dir_at = |theta: f64| -> Vector {
        seat_frame.x().vector() * theta.cos() + seat_frame.y().vector() * theta.sin()
    };
    // Every wedge chart is anchored at the arc's start: the window then
    // lives in (0, sweep), clear of the chart's own seam whatever part of
    // the rim the arc covers.
    let wedge_frame_at = |origin: Point| -> OgeomResult<Frame> {
        Frame::new(
            origin,
            Direction::new(seat.up, tol)?,
            Direction::new(dir_at(theta0), tol)?,
            tol,
        )
    };
    let apex_at = |theta: f64| seat.centre + dir_at(theta) * seat.radius;
    let wall_at = |theta: f64| tube_level + dir_at(theta) * seat.radius;
    let cap_at = |theta: f64| seat.centre + dir_at(theta) * tube_rho;

    // Six shared corner vertices, nine shared edges: a wire chains through
    // vertex objects, and neighbouring faces sew through shared edges.
    let va0 = make_vertex(model, apex_at(theta0)).shape;
    let va1 = make_vertex(model, apex_at(theta1)).shape;
    let vw0 = make_vertex(model, wall_at(theta0)).shape;
    let vw1 = make_vertex(model, wall_at(theta1)).shape;
    let vc0 = make_vertex(model, cap_at(theta0)).shape;
    let vc1 = make_vertex(model, cap_at(theta1)).shape;

    let ring_arc = |model: &mut Model,
                    origin: Point,
                    r: f64,
                    from: &Shape,
                    to: &Shape|
     -> OgeomResult<Shape> {
        let circle = Circle::new(wedge_frame_at(origin)?, r, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        Ok(make_edge_between(model, curve, (0.0, sweep), from, to, tol)?.shape)
    };
    // The quarter meridian at `theta`, from the wall contact to the cap
    // contact — the tube's own circle, in the meridian plane.
    let quarter = |model: &mut Model, theta: f64, from: &Shape, to: &Shape| -> OgeomResult<Shape> {
        let dir = dir_at(theta);
        let centre_m = tube_level + dir * tube_rho;
        let x_dir = Direction::new(dir * (seat.sigma * seat.tau), tol)?;
        let z_dir = Direction::new(dir.cross(seat.up) * seat.sigma, tol)?;
        let frame = Frame::new(centre_m, z_dir, x_dir, tol)?;
        let circle = Circle::new(frame, radius, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        Ok(make_edge_between(
            model,
            curve,
            (0.0, core::f64::consts::FRAC_PI_2),
            from,
            to,
            tol,
        )?
        .shape)
    };

    let arc_apex = ring_arc(model, seat.centre, seat.radius, &va0, &va1)?;
    let arc_wall = ring_arc(model, tube_level, seat.radius, &vw0, &vw1)?;
    let arc_cap = ring_arc(model, seat.centre, tube_rho, &vc0, &vc1)?;
    let ruling0 = segment_between(model, (&va0, apex_at(theta0)), (&vw0, wall_at(theta0)), tol)?;
    let ruling1 = segment_between(model, (&vw1, wall_at(theta1)), (&va1, apex_at(theta1)), tol)?;
    let radial0 = segment_between(model, (&vc0, cap_at(theta0)), (&va0, apex_at(theta0)), tol)?;
    let radial1 = segment_between(model, (&va1, apex_at(theta1)), (&vc1, cap_at(theta1)), tol)?;
    let quarter0 = quarter(model, theta0, &vw0, &vc0)?;
    let quarter1 = quarter(model, theta1, &vw1, &vc1)?;

    // The wall patch: the solid wall's own chart between the contact ring
    // and the rim, over the window.
    let wall_patch = {
        let below = radius + 1.0;
        let origin = seat.centre - seat.up * (seat.tau * below);
        let surface: SurfaceGeometry = CylinderSurface::new(
            Cylinder::new(wedge_frame_at(origin)?, seat.radius, tol)?,
            (0.0, 2.0 * below),
        )?
        .into();
        let edges = vec![
            arc_wall.clone(),
            ruling1.clone(),
            arc_apex.reversed(),
            ruling0.clone(),
        ];
        let face = face_from_edges(model, surface, &edges, tol)?;
        if seat.sigma * seat.tau < 0.0 {
            face.reversed()
        } else {
            face
        }
    };

    // The cap sector: the cap plane between the rim and the contact ring.
    let cap_sector = {
        let plane = Plane::through(seat.centre, Direction::new(seat.up * seat.tau, tol)?);
        let reach = (seat.radius + radius + 1.0) * 2.0;
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let edges = vec![
            arc_apex.clone(),
            radial1.clone(),
            arc_cap.reversed(),
            radial0.clone(),
        ];
        face_from_edges(model, surface, &edges, tol)?
    };

    // The blend: the quarter-tube patch over the window, reversed as the
    // full rim's band is — its natural normal points into the wedge.
    let blend_patch = {
        let surface: SurfaceGeometry = TorusSurface::new(Torus::new(
            wedge_frame_at(tube_level)?,
            tube_rho,
            radius,
            tol,
        )?)
        .into();
        let edges = vec![
            arc_wall.clone(),
            quarter1.clone(),
            arc_cap.reversed(),
            quarter0.reversed(),
        ];
        face_from_edges(model, surface, &edges, tol)?.reversed()
    };

    // The end caps: planar triangles in the meridian half-planes, outward
    // along the rim's travel — behind it at the start, ahead at the end.
    let tangent_at = |theta: f64| seat.up.cross(dir_at(theta));
    let end_cap =
        |model: &mut Model, theta: f64, outward: Vector, edges: Vec<Shape>| -> OgeomResult<Shape> {
            let plane = Plane::through(apex_at(theta), Direction::new(outward, tol)?);
            let reach = (seat.radius + radius + 1.0) * 2.0;
            let surface: SurfaceGeometry =
                PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
            face_from_edges(model, surface, &edges, tol)
        };
    let cap_start = end_cap(
        model,
        theta0,
        -tangent_at(theta0),
        vec![ruling0.clone(), quarter0.clone(), radial0.clone()],
    )?;
    let cap_end = end_cap(
        model,
        theta1,
        tangent_at(theta1),
        vec![ruling1.clone(), radial1.clone(), quarter1.reversed()],
    )?;

    let faces = [wall_patch, cap_sector, blend_patch, cap_start, cap_end];
    apply_wedge(model, solid, Some(edge), &faces, seat.additive(), tol)
}
