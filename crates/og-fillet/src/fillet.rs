//! Constant-radius fillets, where the rolling ball's envelope is closed-form.
//!
//! A ball of constant radius rolling along a straight edge between two planar
//! faces stays tangent to both, and its envelope is exactly a cylinder — the
//! first blend that needs a *new* surface, and the one whose new surface
//! costs nothing. The construction is the chamfer's wedge with the bevel
//! plane exchanged for that cylinder: two legs running along the faces to the
//! tangency lines, two caps closed by arcs, and the cylindrical face between
//! them, all with exact same-parameter pcurves. The boolean then does what it
//! did for the chamfer — the legs melt into the solid's own faces by
//! same-domain resolution, and the cylinder stays as the blend.

use crate::support::{face_from_edges, planar_face, planar_seat, segment_between, subtract_wedge};
use og_algo::{Built, make_edge_between, make_vertex};
use og_core::{OgResult, Tolerances, og_bail};
use og_geom::{CircleCurve, Curve, CylinderSurface, PlaneSurface};
use og_math::{Circle, Cylinder, Direction, Frame, Plane, Point, Vector};
use og_topo::{Model, Shape};

/// Round a straight convex edge of a solid with a constant-radius blend
/// tangent to both of its planar faces.
///
/// The result is the boolean difference with a wedge whose curved face is the
/// tangent cylinder, so the history reads as a cut: the two faces are
/// modified into their trimmed pieces, and the edge's neighbourhood gains the
/// blend face.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the edge is
/// not straight, not convex, not shared by exactly two planar faces of
/// `solid`, or `radius` is not a usable length.
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
