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
    edge_curve, face_from_edges, planar_face, planar_seat, segment_between, subtract_wedge,
};
use og_algo::{
    Built, attach_pcurve, make_edge, make_edge_between, make_face, make_revolution_band,
    make_vertex, make_wire,
};
use og_core::{OgResult, Tolerances, og_bail};
use og_geom::{CircleCurve, Curve, CylinderSurface, PlaneSurface, SurfaceGeometry, TorusSurface};
use og_math::{Circle, Cylinder, Direction, Frame, Plane, Point, Torus, Vector};
use og_topo::{
    EdgeRepr, Filter, Location, Model, NodeData, Orientation, Shape, ShapeType, explore,
};

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
/// [`OgError::Construction`](og_core::OgError::Construction) if the edge is
/// neither of the seats above, is concave, or `radius` is not a usable
/// length.
pub fn fillet_edge(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgResult<Built> {
    if !radius.is_finite() || radius <= tol.confusion() {
        og_bail!(Construction, "a fillet of radius {radius} rounds nothing");
    }
    let (curve, _) = edge_curve(model, edge)?;
    match curve {
        Curve::Line(_) => planar_fillet(model, solid, edge, radius, tol),
        Curve::Circle(c) => revolved_fillet(model, solid, edge, &c, radius, tol),
        _ => og_bail!(
            Construction,
            "filleting an edge that is neither straight nor circular needs \
             the marching blend machinery"
        ),
    }
}

/// The straight-edge fillet: the rolling ball's envelope is a cylinder.
fn planar_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgResult<Built> {
    let seat = planar_seat(model, solid, edge, tol)?;

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
    let a = seat.leg(first, tol)?;
    let b = seat.leg(second, tol)?;

    // The rolling ball's centre line: along the bisector, at the distance
    // where the ball touches both planes.
    let bisector = {
        let u = a + b;
        let m = u.magnitude();
        if m <= tol.angular() {
            og_bail!(Construction, "the faces meet too sharply to seat a fillet");
        }
        u / m
    };
    let depth = -n1.dot(bisector);
    if depth <= tol.angular() {
        og_bail!(Construction, "the faces meet too sharply to seat a fillet");
    }
    let centre = seat.start + bisector * (radius / depth);
    let sweep = seat.along.dot(n1.cross(n2)).atan2(n1.dot(n2));
    if sweep <= tol.angular() {
        og_bail!(Construction, "the faces are parallel; there is no corner");
    }

    let travel = seat.end - seat.start;
    let length = travel.magnitude();
    let apex0 = seat.start;
    let apex1 = seat.end;
    let contact_a0 = centre + n1 * radius;
    let contact_b0 = centre + n2 * radius;
    let contact_a1 = contact_a0 + travel;
    let contact_b1 = contact_b0 + travel;

    let along_dir = Direction::new(seat.along, tol)?;
    let radial_dir = Direction::new(n1, tol)?;

    // An arc of the blend circle at height `h` along the edge, running from
    // the first face's contact vertex to the second's. Its frame is the
    // cylinder's frame translated up the axis, which is what makes its pcurve
    // on the cylinder a straight line in `(u, v)` with the arc's own
    // parameter.
    let arc = |model: &mut Model, h: f64, from: &Shape, to: &Shape| -> OgResult<Shape> {
        let frame = Frame::new(centre + seat.along * h, along_dir, radial_dir, tol)?;
        let circle = Circle::new(frame, radius, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        Ok(make_edge_between(model, curve, (0.0, sweep), from, to, tol)?.shape)
    };

    // The legs, coplanar with the solid's own faces and aligned with them —
    // the same-domain case the boolean resolves by melting them together.
    let leg_a = planar_face(model, &[apex0, contact_a0, contact_a1, apex1], n1, tol)?;
    let leg_b = planar_face(model, &[apex0, contact_b0, contact_b1, apex1], n2, tol)?;

    // The caps: two segments to the tangency points, closed by the arc.
    let cap = |model: &mut Model,
               apex: Point,
               ca: Point,
               cb: Point,
               h: f64,
               outward: Vector|
     -> OgResult<Shape> {
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
    subtract_wedge(model, solid, edge, &faces, tol)
}

/// The circular-rim fillet: a ball rolling around the rim where a cylindrical
/// wall meets a perpendicular planar cap traces a torus.
///
/// The wedge is the revolved cousin of the straight case's prism: a band of
/// the wall down to the tangency parallel, an annulus of the cap in to the
/// tangency circle, and the quarter-torus between them. Being a full
/// revolution it needs no end caps — and its bands carry seams, which is what
/// [`make_revolution_band`] exists to build correctly.
fn revolved_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    rim: &CircleCurve,
    radius: f64,
    tol: Tolerances,
) -> OgResult<Built> {
    let rim_circle = rim.circle();
    let rim_centre = rim_circle.centre();
    let rim_radius = rim_circle.radius();
    let rim_axis = rim_circle.frame().z().vector();
    if radius >= rim_radius - tol.confusion() {
        og_bail!(
            Construction,
            "a fillet of radius {radius} swallows the axis of a rim of \
             radius {rim_radius}"
        );
    }

    // The two faces at the rim: exactly one plane and one coaxial cylinder.
    let mut cap_normal: Option<Vector> = None;
    let mut wall_outward_radial: Option<bool> = None;
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let touches = explore(model, &face, Filter::OfType(ShapeType::Edge))?
            .iter()
            .any(|e| e.node() == edge.node());
        if !touches {
            continue;
        }
        let Some(node) = model.node(&face) else {
            og_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            og_bail!(Construction, "face node holds no face data");
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
                    og_bail!(
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
                    og_bail!(
                        Construction,
                        "the rim's wall is not the coaxial cylinder through \
                         it; that seat needs the marching blend machinery"
                    );
                }
                // A cylinder face's natural normal points away from its axis;
                // an external rim needs exactly that.
                wall_outward_radial = Some(!reversed);
            }
            Some(_) => og_bail!(
                Construction,
                "the rim meets a face that is neither plane nor cylinder; \
                 that seat needs the marching blend machinery"
            ),
            None => og_bail!(Dangling, "face refers to a surface not in this model"),
        }
    }
    let (Some(up_raw), Some(outward)) = (cap_normal, wall_outward_radial) else {
        og_bail!(
            Construction,
            "a rim fillet needs the edge shared by one planar cap and one \
             cylindrical wall"
        );
    };
    if !outward {
        og_bail!(
            Construction,
            "the rim's wall faces inward — a hole's rim — where the \
             subtractive wedge lies in empty space; the concave blend is \
             recorded in the deferred table"
        );
    }

    // Everything below is written with `up` the cap's outward normal: the
    // material sits below the cap and inside the wall, and the ball rolls in
    // the quarter between them.
    let up = up_raw / up_raw.magnitude();
    let up_dir = Direction::new(up, tol)?;
    let x_ref = if up.cross(Vector::X).magnitude() > 0.5 {
        Direction::new(up.cross(Vector::X), tol)?
    } else {
        Direction::new(up.cross(Vector::Y), tol)?
    };
    let frame_at = |origin: Point| Frame::new(origin, up_dir, x_ref, tol);
    let tube_centre = rim_centre - up * radius;

    let ring = |model: &mut Model, origin: Point, r: f64| -> OgResult<Shape> {
        let circle = Circle::new(frame_at(origin)?, r, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        let domain = og_geom::Curve3d::domain(&curve);
        Ok(make_edge(model, curve, domain, tol)?.shape)
    };
    let apex_ring = ring(model, rim_centre, rim_radius)?;
    let wall_ring = ring(model, tube_centre, rim_radius)?;
    let cap_ring = ring(model, rim_centre, rim_radius - radius)?;

    // The wall band, from the tangency parallel up to the rim — coincident
    // with the solid's own wall, which same-domain resolution melts.
    let wall_band = {
        let surface: SurfaceGeometry = CylinderSurface::new(
            Cylinder::new(frame_at(tube_centre)?, rim_radius, tol)?,
            (0.0, radius),
        )?
        .into();
        make_revolution_band(model, &surface, &wall_ring, &apex_ring, tol)?
    };

    // The blend itself: the quarter of the tube between the tangencies. Its
    // natural normal points away from the tube's centre circle — into the
    // wedge — so it enters the shell reversed.
    let blend_band = {
        let surface: SurfaceGeometry = TorusSurface::new(Torus::new(
            frame_at(tube_centre)?,
            rim_radius - radius,
            radius,
            tol,
        )?)
        .into();
        make_revolution_band(model, &surface, &wall_ring, &cap_ring, tol)?
    };

    // The cap annulus, rim to tangency circle, coplanar with the solid's cap.
    let annulus = {
        let plane = Plane::through(rim_centre, up_dir);
        let reach = rim_radius * 2.0;
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let outer = make_wire(model, std::slice::from_ref(&apex_ring), tol)?.shape;
        let inner = make_wire(model, std::slice::from_ref(&cap_ring), tol)?.shape;
        let face = make_face(model, surface.clone(), &[outer, inner], tol)?.shape;
        let surface_id = {
            let Some(node) = model.node(&face) else {
                og_bail!(Dangling, "the face just built is not in this model");
            };
            let NodeData::Face(data) = node.data() else {
                og_bail!(Construction, "the face holds no face data");
            };
            data.surface
        };
        for pedge in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
            let (curve, prange) = {
                let Some(node) = model.node(&pedge) else {
                    og_bail!(Dangling, "edge is not in this model");
                };
                let Some(data) = node.data().as_edge() else {
                    og_bail!(Construction, "edge node holds no edge data");
                };
                let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                    og_bail!(Construction, "an annulus edge has no 3D curve");
                };
                let Some(geometry) = model.geometry().curve(*curve) else {
                    og_bail!(Dangling, "curve is not in this model");
                };
                (geometry.clone(), *range)
            };
            let Some(pcurve) = og_intersect::exact_pcurve_of(&curve, &surface, tol) else {
                og_bail!(
                    Construction,
                    "an annulus edge has no closed-form pcurve on its plane"
                );
            };
            attach_pcurve(
                model,
                &pedge,
                pcurve,
                surface_id,
                Location::identity(),
                prange,
            )?;
        }
        face
    };

    let faces = [wall_band, annulus, blend_band.reversed()];
    subtract_wedge(model, solid, edge, &faces, tol)
}
