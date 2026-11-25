//! Each blended edge's bevel on its seat: the planes, rims and creases.

use super::*;

/// The two distances of `spec` on `edge`'s faces in `faces`' order, `None`
/// for an angle, whose second distance [`angled`] reads from the section.
pub(super) fn distances(curved: &Curved, faces: [usize; 2], spec: &Chamfer) -> Option<[f64; 2]> {
    match spec {
        Chamfer::Symmetric(d) => Some([*d, *d]),
        Chamfer::Distances {
            face,
            on_face,
            on_other,
        } => {
            let k = faces
                .iter()
                .position(|f| curved.faces[*f].shape.is_same(face))?;
            let mut d = [*on_other; 2];
            d[k] = *on_face;
            Some(d)
        }
        Chamfer::Angle { .. } => None,
    }
}

/// The contacts and bevel of `spec` on `edge`, `None` past the seats built
/// here. A straight edge between planes or a rim may be concave only where
/// `concave` says so; the rings build convex bevels alone.
#[allow(clippy::too_many_arguments, reason = "the seat and the request")]
pub(super) fn bevel_of(
    model: &mut Model,
    solid: &Shape,
    curved: &Curved,
    edge: Id,
    shape: &Shape,
    spec: &Chamfer,
    concave: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Strip>> {
    let Some(e) = curved.edges.get(&edge) else {
        return Ok(None);
    };
    let [f0, f1] = faces_of(curved, edge)[..] else {
        return Ok(None);
    };
    let faces = [f0, f1];
    match &e.curve {
        Curve::Line(_) => {
            // Which way each face runs from the edge, read on the face: a
            // face that is not convex can hold its other vertices on either
            // side of the edge's line.
            let t = f64::midpoint(e.range.0, e.range.1);
            let p = e.curve.point_at(t, tol)?;
            let into = [
                into_face(model, curved, (f0, edge), (p, t), tol)?,
                into_face(model, curved, (f1, edge), (p, t), tol)?,
            ];
            let d = match spec {
                Chamfer::Angle {
                    face,
                    distance,
                    angle,
                } => {
                    let Some(k) = faces
                        .iter()
                        .position(|f| curved.faces[*f].shape.is_same(face))
                    else {
                        return Ok(None);
                    };
                    let Some(other) =
                        angled(curved, edge, e, (faces, into), k, *distance, *angle, tol)
                    else {
                        return Ok(None);
                    };
                    let mut d = [other; 2];
                    d[k] = *distance;
                    d
                }
                _ => {
                    let Some(d) = distances(curved, faces, spec) else {
                        return Ok(None);
                    };
                    d
                }
            };
            if d.iter().any(|x| !(x.is_finite() && *x > tol.confusion())) {
                return Ok(None);
            }
            straight(curved, edge, e, (faces, into), d, concave, tol)
        }
        Curve::Circle(circle) => match legs(model, solid, shape, spec, tol) {
            Ok(Legs::Rim(seat, on_wall, on_cap)) => {
                rim(curved, edge, faces, &seat, on_wall, on_cap, concave, tol)
            }
            Ok(Legs::Sphere(seat, on_wall, on_cap)) => {
                sphere_rim(curved, edge, faces, &seat, on_wall, on_cap, tol)
            }
            _ => torus_rim(curved, edge, circle, faces, spec, tol),
        },
        _ => extruded_crease(curved, edge, e, faces, spec, tol),
    }
}

/// A distance-angle chamfer's second distance on a straight seat: in the
/// edge's section, the bevel leaves face `k`'s contact, `distance` in from
/// the edge, at `angle` from that face into the material, and runs to the
/// other face's trace, a line or the cylinder's circle; the distance there
/// is measured along it from the edge. `None` where the bevel does not
/// reach the other face, or a face is not a plane or a cylinder along the
/// edge.
#[allow(clippy::too_many_arguments, reason = "the seat's whole section")]
fn angled(
    curved: &Curved,
    edge: Id,
    e: &CurvedEdge,
    (faces, into): ([usize; 2], [Option<Vector>; 2]),
    k: usize,
    distance: f64,
    angle: f64,
    tol: Tolerances,
) -> Option<f64> {
    if !(angle.is_finite() && angle > tol.angular() && angle < std::f64::consts::PI)
        || !(distance.is_finite() && distance > tol.confusion())
    {
        return None;
    }
    let (a, b) = (curved.points[&e.ends.0], curved.points[&e.ends.1]);
    let u = (b - a) / a.distance(b);
    let at = a + (b - a) * 0.5;
    let normals = [
        outward(curved, faces[0], at)?,
        outward(curved, faces[1], at)?,
    ];
    // Each face's trace in the section: a line from the edge, or a circle.
    enum Trace {
        Line(Vector),
        Circle(Point, f64, Vector),
    }
    let trace = |j: usize| -> Option<Trace> {
        match &curved.faces[faces[j]].surface {
            SurfaceGeometry::Plane(_) => {
                let m = normals[j].cross(u);
                let towards = match into[j] {
                    Some(d) => d,
                    None => away(curved, faces[j], edge)? - at,
                };
                let m = if towards.dot(m) < 0.0 { -m } else { m };
                Some(Trace::Line(m / m.magnitude()))
            }
            SurfaceGeometry::Cylinder(c) => {
                let cyl = c.cylinder();
                let axis = cyl.axis();
                let (o, ax) = (axis.location, axis.direction.vector());
                if u.cross(ax).magnitude() > tol.angular() {
                    return None;
                }
                let centre = o + ax * (at - o).dot(ax);
                Some(Trace::Circle(centre, cyl.radius(), ax))
            }
            _ => None,
        }
    };
    let (named, other) = (trace(k)?, trace(1 - k)?);
    let into = |p: Point| -> Option<Vector> { Some(-outward(curved, faces[k], p)?) };
    // Where the bevel leaves the named face, and the way back along it.
    let (start, back) = match named {
        Trace::Line(m) => (at + m * distance, -m),
        Trace::Circle(centre, radius, ax) => {
            let r = at - centre;
            let turn = |sign: f64| {
                let w = ax.cross(r);
                centre + r * (distance / radius).cos() + w * (sign * (distance / radius).sin())
            };
            let n_other = normals[1 - k];
            let (plus, minus) = (turn(1.0), turn(-1.0));
            let p = if (plus - at).dot(n_other) < (minus - at).dot(n_other) {
                plus
            } else {
                minus
            };
            let w = ax.cross(p - centre);
            let w = w / w.magnitude();
            (p, if w.dot(at - p) > 0.0 { w } else { -w })
        }
    };
    let t = back * angle.cos() + into(start)? * angle.sin();
    match other {
        Trace::Line(m) => {
            let n = normals[1 - k];
            let s = n.dot(at - start) / n.dot(t);
            let hit = start + t * s;
            let lambda = (hit - at).dot(m);
            (s > 0.0 && lambda > tol.confusion()).then_some(lambda)
        }
        Trace::Circle(centre, radius, _) => {
            let w = start - centre;
            let (qa, qb, qc) = (t.dot(t), 2.0 * w.dot(t), w.dot(w) - radius * radius);
            let disc = qb * qb - 4.0 * qa * qc;
            if disc < 0.0 {
                return None;
            }
            let root = disc.sqrt();
            let s = [(-qb - root) / (2.0 * qa), (-qb + root) / (2.0 * qa)]
                .into_iter()
                .filter(|s| *s > tol.confusion())
                .fold(f64::INFINITY, f64::min);
            if !s.is_finite() {
                return None;
            }
            let hit = start + t * s;
            let (x, y) = (at - centre, hit - centre);
            let arc = radius
                * (x.dot(y) / (x.magnitude() * y.magnitude()))
                    .clamp(-1.0, 1.0)
                    .acos();
            (arc > tol.confusion()).then_some(arc)
        }
    }
}

/// A straight edge between two planes, or a plane and a cylinder it runs
/// along: each contact a segment the length of the edge and half again past
/// each end, so the corners at its ends lie on it.
#[allow(clippy::too_many_arguments, reason = "the seat's whole section")]
fn straight(
    curved: &Curved,
    edge: Id,
    e: &CurvedEdge,
    (faces, into): ([usize; 2], [Option<Vector>; 2]),
    d: [f64; 2],
    concave: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Strip>> {
    let (mut a, mut b) = (curved.points[&e.ends.0], curved.points[&e.ends.1]);
    // Between two planes the edge is their line, exactly: its ends taken
    // onto it, so a vertex standing off its faces (inexact input) does not
    // tilt the contacts off their faces.
    if let (SurfaceGeometry::Plane(p0), SurfaceGeometry::Plane(p1)) = (
        &curved.faces[faces[0]].surface,
        &curved.faces[faces[1]].surface,
    ) && let Some((on, dir)) = planes_line(&p0.plane(), &p1.plane(), a, tol)
    {
        let onto = |p: Point| on + dir * (p - on).dot(dir);
        (a, b) = (onto(a), onto(b));
    }
    let length = a.distance(b);
    if length <= tol.confusion() {
        return Ok(None);
    }
    let u = (b - a) / length;
    let run = u * (length / 2.0);
    let mid = a + (b - a) * 0.5;
    let (Some(n0), Some(n1)) = (
        outward(curved, faces[0], mid),
        outward(curved, faces[1], mid),
    ) else {
        return Ok(None);
    };
    let normals = [n0, n1];
    let mut lines: Vec<(Point, Point)> = Vec::with_capacity(2);
    for k in 0..2 {
        let (fi, other) = (faces[k], normals[1 - k]);
        let face = &curved.faces[fi];
        let shift: Vector = match &face.surface {
            SurfaceGeometry::Plane(_) => {
                let m = normals[k].cross(u);
                let towards = match into[k] {
                    Some(d) => d,
                    None => match away(curved, fi, edge) {
                        Some(inside) => inside - mid,
                        None => return Ok(None),
                    },
                };
                let m = if towards.dot(m) < 0.0 { -m } else { m };
                // A convex edge only, unless a concave one is asked for:
                // into the face goes under the other.
                if !concave && m.dot(other) >= 0.0 {
                    return Ok(None);
                }
                m / m.magnitude() * d[k]
            }
            SurfaceGeometry::Cylinder(c) => {
                let cyl = c.cylinder();
                let axis = cyl.axis();
                let (o, ax) = (axis.location, axis.direction.vector());
                if u.cross(ax).magnitude() > tol.angular() {
                    return Ok(None);
                }
                let radius = cyl.radius();
                let r = (mid - o) - ax * (mid - o).dot(ax);
                let angle = d[k] / radius;
                if angle >= std::f64::consts::FRAC_PI_2 {
                    return Ok(None);
                }
                // Round the axis, the way into the face: away from the other
                // face's outward side.
                let turn = |sign: f64| {
                    let w = ax.cross(r);
                    r * angle.cos() + w * (sign * angle.sin()) - r
                };
                let (plus, minus) = (turn(1.0), turn(-1.0));
                let shift = if plus.dot(other) < minus.dot(other) {
                    plus
                } else {
                    minus
                };
                if shift.dot(other) >= 0.0 {
                    return Ok(None);
                }
                shift
            }
            _ => return Ok(None),
        };
        lines.push((a + shift - run, b + shift + run));
    }
    let contacts = [
        Curve::from(LineCurve::segment(lines[0].0, lines[0].1, tol)?),
        Curve::from(LineCurve::segment(lines[1].0, lines[1].1, tol)?),
    ];
    // The bevel: the plane through both contacts.
    let across = lines[1].0 - lines[0].0;
    let n = u.cross(across);
    if n.magnitude() <= tol.confusion() {
        return Ok(None);
    }
    let frame = Frame::new(
        lines[0].0,
        Direction::new(n, tol)?,
        Direction::new(u, tol)?,
        tol,
    )?;
    Ok(Some(Strip {
        edge,
        faces,
        contacts,
        surface: PlaneSurface::new(Plane::new(frame)).into(),
        section: Section::Line,
    }))
}

/// The line two planes meet along: its point nearest `near` and its unit
/// direction; `None` for planes all but parallel.
fn planes_line(
    a: &ogeom_math::Plane,
    b: &ogeom_math::Plane,
    near: Point,
    tol: Tolerances,
) -> Option<(Point, Vector)> {
    let (na, nb) = (a.normal().vector(), b.normal().vector());
    let dir = na.cross(nb);
    let m = dir.magnitude();
    if m <= tol.angular() {
        return None;
    }
    let dir = dir / m;
    // The point of both planes nearest `near`: `near` moved onto each plane
    // in turn along the other's in-plane direction square to the line.
    let (da, db) = (
        na.dot(a.origin().to_vector()),
        nb.dot(b.origin().to_vector()),
    );
    let (ca, cb) = (da - na.dot(near.to_vector()), db - nb.dot(near.to_vector()));
    // Solve near + x na + y nb on both planes.
    let (k, aa, bb) = (na.dot(nb), na.dot(na), nb.dot(nb));
    let det = aa * bb - k * k;
    if det.abs() <= tol.angular() {
        return None;
    }
    let x = (ca * bb - cb * k) / det;
    let y = (cb * aa - ca * k) / det;
    Some((near + na * x + nb * y, dir))
}

/// A rim between a planar cap and a coaxial cylinder or cone, as the wedge
/// path reads it ([`legs`]): the ring along the cap and the ring down the
/// wall, and the cone through both.
#[allow(clippy::too_many_arguments, reason = "the seat's whole section")]
fn rim(
    curved: &Curved,
    edge: Id,
    faces: [usize; 2],
    seat: &RevolvedSeat,
    on_wall: f64,
    on_cap: f64,
    concave: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Strip>> {
    if seat.additive() && !concave {
        return Ok(None);
    }
    let cap = faces
        .iter()
        .position(|f| curved.faces[*f].shape.is_same(&seat.cap_face));
    let wall = faces
        .iter()
        .position(|f| curved.faces[*f].shape.is_same(&seat.wall_face));
    let (Some(cap), Some(wall)) = (cap, wall) else {
        return Ok(None);
    };
    let cap_rho = seat.sigma.mul_add(-(seat.tau * on_cap), seat.radius);
    let depth = wall_depth(seat, on_wall);
    let wall_rho = seat.radius_at(-seat.tau * depth);
    if cap_rho <= tol.confusion() || wall_rho <= tol.confusion() {
        return Ok(None);
    }
    let ring = |centre: Point, radius: f64| -> OgeomResult<Curve> {
        Ok(CircleCurve::new(Circle::new(seat.frame_at(centre, tol)?, radius, tol)?).into())
    };
    let mut contacts = [ring(seat.centre, cap_rho)?, ring(seat.centre, cap_rho)?];
    contacts[cap] = ring(seat.centre, cap_rho)?;
    contacts[wall] = ring(seat.centre + seat.up * (-seat.tau * depth), wall_rho)?;
    let surface = rim_bevel_surface(seat, depth, cap_rho, wall_rho, tol)?;
    Ok(Some(Strip {
        edge,
        faces,
        contacts,
        surface,
        section: Section::Line,
    }))
}

/// A rim between a planar cap and a sphere centred on its axis, as the
/// wedge path reads it ([`legs`]): the ring along the cap, the ring down
/// the sphere's great circle, and the cone through both. `None` where the
/// legs leave no such bevel.
fn sphere_rim(
    curved: &Curved,
    edge: Id,
    faces: [usize; 2],
    seat: &SphereSeat,
    on_wall: f64,
    on_cap: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Strip>> {
    let on = |face: &Shape| {
        faces
            .iter()
            .position(|f| curved.faces[*f].shape.is_same(face))
    };
    let (Some(cap), Some(wall)) = (on(&seat.cap_face), on(&seat.wall_face)) else {
        return Ok(None);
    };
    if cap == wall || on_wall >= seat.to_pole() {
        return Ok(None);
    }
    let Ok(bevel) = sphere_bevel_rings(seat, on_wall, on_cap, tol) else {
        return Ok(None);
    };
    let ring = |centre: Point, radius: f64| -> OgeomResult<Curve> {
        let frame = Frame::new(centre, Direction::new(seat.up, tol)?, seat.x_ref, tol)?;
        Ok(CircleCurve::new(Circle::new(frame, radius, tol)?).into())
    };
    let mut contacts = [
        ring(seat.centre, bevel.cap_rho)?,
        ring(seat.centre, bevel.cap_rho)?,
    ];
    contacts[wall] = ring(seat.centre + seat.up * bevel.wall_h, bevel.wall_rho)?;
    Ok(Some(Strip {
        edge,
        faces,
        contacts,
        surface: bevel.surface,
        section: Section::Line,
    }))
}

/// A rim between a planar cap square to a torus's axis and the torus, the
/// material inside its tube or outside it (a groove): in the meridian
/// half-plane through a point of the rim, `(rho, z)` from the torus's
/// centre along its axis, the leg on the cap runs straight along it and the
/// leg on the wall round the tube's circle, and the bevel is the cone
/// through the two contact rings. `None` where the faces are not such a
/// plane and torus, the rim is not coaxial, the edge is not convex, or the
/// legs leave no bevel in the material.
fn torus_rim(
    curved: &Curved,
    edge: Id,
    rim: &CircleCurve,
    faces: [usize; 2],
    spec: &Chamfer,
    tol: Tolerances,
) -> OgeomResult<Option<Strip>> {
    let (eps, ang) = (tol.confusion(), tol.angular());
    let surface = |k: usize| &curved.faces[faces[k]].surface;
    let (cap, wall) = match (surface(0), surface(1)) {
        (SurfaceGeometry::Plane(_), SurfaceGeometry::Torus(_)) => (0, 1),
        (SurfaceGeometry::Torus(_), SurfaceGeometry::Plane(_)) => (1, 0),
        _ => return Ok(None),
    };
    let (SurfaceGeometry::Plane(plane), SurfaceGeometry::Torus(torus)) =
        (surface(cap), surface(wall))
    else {
        return Ok(None);
    };
    // The material inside the tube, or outside it (a groove).
    let within: f64 = if curved.faces[faces[wall]].reversed {
        -1.0
    } else {
        1.0
    };
    let torus = torus.torus();
    let frame = torus.frame();
    let (o, a) = (frame.origin(), frame.z().vector());
    let (big, small) = (torus.major_radius(), torus.minor_radius());
    let circle = rim.circle();
    let centre = circle.centre();
    let off = centre - o;
    let up = plane.plane().normal().vector()
        * if curved.faces[faces[cap]].reversed {
            -1.0
        } else {
            1.0
        };
    if up.cross(a).magnitude() > ang
        || circle.frame().z().vector().cross(a).magnitude() > ang
        || (off - a * off.dot(a)).magnitude() > eps
    {
        return Ok(None);
    }
    let rim_rho = circle.radius();
    let rim_z = off.dot(a);
    if ((rim_rho - big).hypot(rim_z) - small).abs() > eps {
        return Ok(None);
    }
    // The meridian: the cap's outward normal, the wall's at the rim, and
    // which way along the cap from the rim the cap lies (toward the axis
    // where it has no other edge).
    let cap_up = up.dot(a).signum();
    let wall_out = (within * (rim_rho - big) / small, within * rim_z / small);
    let edge_ids: Vec<Id> = curved.faces[faces[cap]]
        .loops
        .iter()
        .flatten()
        .map(|(e, _)| *e)
        .filter(|e| *e != edge)
        .collect();
    let mut side = 0.0_f64;
    for e in &edge_ids {
        let c = &curved.edges[e];
        let p = c.curve.point_at(f64::midpoint(c.range.0, c.range.1), tol)?;
        let q = p - o;
        let gap = (q - a * q.dot(a)).magnitude() - rim_rho;
        if gap.abs() <= eps {
            // The rim's other arcs.
            continue;
        }
        let s = gap.signum();
        if side == 0.0 {
            side = s;
        } else if side != s {
            return Ok(None);
        }
    }
    if side == 0.0 {
        side = -1.0;
    }
    // Convex: the cap leaves the rim behind the wall's outward side.
    if side * wall_out.0 >= -ang {
        return Ok(None);
    }
    // Round the tube from the rim, away from the cap's outward side.
    let alpha = rim_z.atan2(rim_rho - big);
    let turn: f64 = if alpha.cos() * cap_up < 0.0 {
        1.0
    } else {
        -1.0
    };
    let on_tube = |t: f64| {
        let b = turn.mul_add(t / small, alpha);
        (small.mul_add(b.cos(), big), small * b.sin())
    };
    let named = |face: &Shape| {
        if curved.faces[faces[cap]].shape.is_same(face) {
            Some(true)
        } else if curved.faces[faces[wall]].shape.is_same(face) {
            Some(false)
        } else {
            None
        }
    };
    // The legs: along the cap, and round the tube.
    let (on_cap, on_wall) = match spec {
        Chamfer::Symmetric(d) => (*d, *d),
        Chamfer::Distances {
            face,
            on_face,
            on_other,
        } => match named(face) {
            Some(true) => (*on_face, *on_other),
            Some(false) => (*on_other, *on_face),
            None => return Ok(None),
        },
        Chamfer::Angle {
            face,
            distance,
            angle,
        } => {
            let (d, phi) = (*distance, *angle);
            match named(face) {
                Some(true) => {
                    // From the cap's contact back toward the rim and into
                    // the material, to where it meets the tube.
                    let c = (side.mul_add(d, rim_rho), rim_z);
                    let dir = (-side * phi.cos(), -cap_up * phi.sin());
                    let rel = (c.0 - big, c.1);
                    let b = rel.0.mul_add(dir.0, rel.1 * dir.1);
                    let k = rel.0.mul_add(rel.0, rel.1 * rel.1) - small * small;
                    let disc = b.mul_add(b, -k);
                    // Inside the tube the contact stands in it and the line
                    // leaves through the far root; outside, it first meets
                    // the tube at the near one.
                    if within * k >= 0.0 || disc < 0.0 {
                        return Ok(None);
                    }
                    let t = -b + within * disc.sqrt();
                    if t <= eps {
                        return Ok(None);
                    }
                    let w = (t.mul_add(dir.0, c.0), t.mul_add(dir.1, c.1));
                    let beta = w.1.atan2(w.0 - big);
                    let mut swept = turn * (beta - alpha);
                    if swept < 0.0 {
                        swept += core::f64::consts::TAU;
                    }
                    (d, swept * small)
                }
                Some(false) => {
                    // From the wall's contact back along the tube toward
                    // the rim and into the material, to the cap's line.
                    let w = on_tube(d);
                    let b = turn.mul_add(d / small, alpha);
                    let back = (turn * b.sin(), -turn * b.cos());
                    let into = (-within * b.cos(), -within * b.sin());
                    let dir = (
                        back.0.mul_add(phi.cos(), into.0 * phi.sin()),
                        back.1.mul_add(phi.cos(), into.1 * phi.sin()),
                    );
                    if dir.1.abs() <= ang {
                        return Ok(None);
                    }
                    let t = (rim_z - w.1) / dir.1;
                    if t <= eps {
                        return Ok(None);
                    }
                    let c_rho = t.mul_add(dir.0, w.0);
                    ((c_rho - rim_rho) * side, d)
                }
                None => return Ok(None),
            }
        }
    };
    if !(on_cap.is_finite() && on_cap > eps && on_wall.is_finite() && on_wall > eps)
        || on_wall >= core::f64::consts::PI * small
    {
        return Ok(None);
    }
    let c = (side.mul_add(on_cap, rim_rho), rim_z);
    let w = on_tube(on_wall);
    if c.0 <= eps || w.0 <= eps || (w.1 - c.1).abs() <= eps {
        return Ok(None);
    }
    // The chord from the wall's contact to the cap's meets the tube's
    // circle again `-2 (W - T).(C - W) / |C - W|^2` along it, which must
    // not fall between the two.
    let chord = (c.0 - w.0, c.1 - w.1);
    let spoke = (w.0 - big, w.1);
    let length2 = chord.0.mul_add(chord.0, chord.1 * chord.1);
    let again = -2.0 * spoke.0.mul_add(chord.0, spoke.1 * chord.1) / length2;
    if again * length2.sqrt() > eps && again < 1.0 {
        return Ok(None);
    }
    let x = circle.frame().x().vector();
    let x = Direction::new(x - a * x.dot(a), tol)?;
    let axis_dir = Direction::new(a, tol)?;
    let ring = |z: f64, radius: f64| -> OgeomResult<Curve> {
        let frame = Frame::new(o + a * z, axis_dir, x, tol)?;
        Ok(CircleCurve::new(Circle::new(frame, radius, tol)?).into())
    };
    let mut contacts = [ring(c.1, c.0)?, ring(c.1, c.0)?];
    contacts[wall] = ring(w.1, w.0)?;
    let surface = coaxial_cone(
        Frame::new(o + a * c.1, axis_dir, x, tol)?,
        c.0,
        (w.0, w.1 - c.1),
        tol,
    )?;
    Ok(Some(Strip {
        edge,
        faces,
        contacts,
        surface,
        section: Section::Line,
    }))
}

/// How many stations the cap's contact along an extruded crease is fitted
/// through.
const CREASE_STATIONS: u32 = 128;

/// How many stations a rational crease is moved down its wall through, as
/// a polynomial curve held within half the confusion distance of it.
const RATIONAL_STATIONS: u32 = 512;

/// How far an extruded crease's ruled bevel runs on past each contact, as
/// a share of its rulings' length.
const BEVEL_RUN_ON: f64 = 0.5;

/// A crease between a plane and the wall extruded from it square to it (a
/// prism's end cap and the side its profile's curve sweeps), the crease
/// any curve of the profile: in the section square to the crease the
/// corner is square, the leg on the wall runs straight down it and the leg
/// on the cap straight in along it. The wall's contact is the crease moved
/// down the wall, exactly on it; the cap's is the crease's offset in the
/// cap's plane, fitted through stations at the crease's own parameters, so
/// it lies in the plane; the bevel is the ruled surface between the two on
/// their shared knots, each ruling in its station's section, on the wall's
/// material side. `None` where the faces are not such a plane and wall, the
/// edge is not convex, or the offset's fit misses the bands' fit target.
fn extruded_crease(
    curved: &Curved,
    edge: Id,
    e: &CurvedEdge,
    faces: [usize; 2],
    spec: &Chamfer,
    tol: Tolerances,
) -> OgeomResult<Option<Strip>> {
    let (eps, ang) = (tol.confusion(), tol.angular());
    let surface = |k: usize| &curved.faces[faces[k]].surface;
    let (cap, wall) = match (surface(0), surface(1)) {
        (SurfaceGeometry::Plane(_), SurfaceGeometry::Extrusion(_)) => (0, 1),
        (SurfaceGeometry::Extrusion(_), SurfaceGeometry::Plane(_)) => (1, 0),
        _ => return Ok(None),
    };
    let (SurfaceGeometry::Plane(plane), SurfaceGeometry::Extrusion(extrusion)) =
        (surface(cap), surface(wall))
    else {
        return Ok(None);
    };
    let plane = plane.plane();
    let up = plane.normal().vector()
        * if curved.faces[faces[cap]].reversed {
            -1.0
        } else {
            1.0
        };
    if extrusion.direction().vector().cross(up).magnitude() > ang {
        return Ok(None);
    }
    let Curve::BSpline(crease) = &e.curve else {
        return Ok(None);
    };
    let (t0, t1) = e.range;
    let at = |k: u32| (t1 - t0).mul_add(f64::from(k) / f64::from(CREASE_STATIONS), t0);
    // A rational crease (a conic's arc) is moved down the wall as a
    // polynomial one through stations at its own parameters, held within
    // half the confusion distance of it between them, so the contact
    // lies on the wall as an exact one does.
    let polynomial;
    let crease = if crease.is_rational() {
        let fine = |k: u32| (t1 - t0).mul_add(f64::from(k) / f64::from(RATIONAL_STATIONS), t0);
        let mut parameters = Vec::with_capacity(RATIONAL_STATIONS as usize + 1);
        let mut points = Vec::with_capacity(RATIONAL_STATIONS as usize + 1);
        for k in 0..=RATIONAL_STATIONS {
            parameters.push(fine(k));
            points.push(e.curve.point_at(fine(k), tol)?);
        }
        let fitted = ogeom_geom::fit::fit_points_at(&parameters, &points, 3, 0.1 * eps, tol)?;
        for k in 0..RATIONAL_STATIONS {
            // Off the crease, not off its point at the same parameter:
            // the two need not run at one speed between stations.
            let at_middle = fitted
                .curve
                .point_at(f64::midpoint(fine(k), fine(k + 1)), tol)?;
            let near = param_on(&e.curve, (fine(k), fine(k + 1)), at_middle, tol)?;
            if e.curve.point_at(near, tol)?.distance(at_middle) > 0.5 * eps {
                return Ok(None);
            }
        }
        polynomial = fitted.curve;
        &polynomial
    } else if crease.knots().is_clamped() {
        crease
    } else {
        // A periodic crease's knots run on past its domain: the wall's
        // contact takes its knots, which must end where the cap's fit
        // does to share them.
        polynomial = clamped_over(crease, e.range, tol)?;
        &polynomial
    };
    // The crease in the cap's plane, and which way into the cap and down
    // the wall: toward the faces' other points.
    let mid = e.curve.point_at(at(CREASE_STATIONS / 2), tol)?;
    for k in 0..=CREASE_STATIONS {
        if plane
            .signed_distance_to(e.curve.point_at(at(k), tol)?)
            .abs()
            > eps
        {
            return Ok(None);
        }
    }
    let (Some(in_wall), Some(wall_out)) = (
        away(curved, faces[wall], edge),
        outward(curved, faces[wall], mid),
    ) else {
        return Ok(None);
    };
    // Convex: the wall runs down from the cap, under its outward side.
    if (in_wall - mid).dot(up) >= -eps {
        return Ok(None);
    }
    let across = |t: f64| -> OgeomResult<Vector> {
        let tangent = e.curve.d1_at(t, tol)?;
        let m = up.cross(tangent);
        Ok(m / m.magnitude())
    };
    // Into the cap from a convex crease: away from the wall's outward
    // side.
    let sign = if wall_out.dot(across(at(CREASE_STATIONS / 2))?) < 0.0 {
        1.0
    } else {
        -1.0
    };
    let named = |face: &Shape| {
        if curved.faces[faces[cap]].shape.is_same(face) {
            Some(true)
        } else if curved.faces[faces[wall]].shape.is_same(face) {
            Some(false)
        } else {
            None
        }
    };
    let (on_cap, on_wall) = match spec {
        Chamfer::Symmetric(d) => (*d, *d),
        Chamfer::Distances {
            face,
            on_face,
            on_other,
        } => match named(face) {
            Some(true) => (*on_face, *on_other),
            Some(false) => (*on_other, *on_face),
            None => return Ok(None),
        },
        // The corner is square: the other leg is the named one's along the
        // angle's tangent.
        Chamfer::Angle {
            face,
            distance,
            angle,
        } => match named(face) {
            Some(true) => (*distance, distance * angle.tan()),
            Some(false) => (distance * angle.tan(), *distance),
            None => return Ok(None),
        },
    };
    if !(on_cap.is_finite() && on_cap > eps && on_wall.is_finite() && on_wall > eps) {
        return Ok(None);
    }
    // The wall's contact: the crease moved down the wall.
    let mut on_side = ogeom_geom::BSplineCurve::new(
        crease.knots().clone(),
        crease
            .control_points()
            .iter()
            .map(|w| w.point() + up * -on_wall)
            .collect(),
        tol,
    )?;
    // The cap's: the offset through stations at the crease's parameters.
    let mut parameters = Vec::with_capacity(CREASE_STATIONS as usize + 1);
    let mut points = Vec::with_capacity(CREASE_STATIONS as usize + 1);
    for k in 0..=CREASE_STATIONS {
        let t = at(k);
        parameters.push(t);
        points.push(e.curve.point_at(t, tol)? + across(t)? * (sign * on_cap));
    }
    let fit = crate::marched::band_fit_target(tol);
    let fitted = ogeom_geom::fit::fit_points_at(&parameters, &points, 3, eps, tol)?;
    if !fitted.met && fitted.error > fit {
        return Ok(None);
    }
    let mut on_top = fitted.curve;
    // The fitted offset between its stations, measured at their middles.
    for k in 0..CREASE_STATIONS {
        let t = f64::midpoint(at(k), at(k + 1));
        let want = e.curve.point_at(t, tol)? + across(t)? * (sign * on_cap);
        if on_top.point_at(t, tol)?.distance(want) > fit {
            return Ok(None);
        }
    }
    if (on_top.domain().0 - on_side.domain().0).abs() > tol.parametric()
        || (on_top.domain().1 - on_side.domain().1).abs() > tol.parametric()
    {
        on_side = on_side.segment(on_top.domain(), tol)?;
    }
    let bevel = ruled_bevel(&mut on_top, &mut on_side, tol)?;
    let mut contacts = [Curve::from(on_top.clone()), Curve::from(on_top)];
    contacts[wall] = Curve::from(on_side);
    Ok(Some(Strip {
        edge,
        faces,
        contacts,
        surface: bevel,
        section: Section::Line,
    }))
}

/// The piece of `curve` over `range`, exactly, its knots clamped at both
/// ends: each end raised to the degree's multiplicity by knot insertion,
/// and the knots and control points past it dropped.
///
/// # Errors
///
/// As knot insertion, and building the piece.
fn clamped_over(
    curve: &ogeom_geom::BSplineCurve,
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<ogeom_geom::BSplineCurve> {
    let p = curve.degree();
    let mut raised = curve.clone();
    for t in [range.0, range.1] {
        let have = raised.knots().multiplicity_of(t);
        if have < p {
            raised = raised.with_knot_inserted(t, p - have, tol)?;
        }
    }
    let k = raised.knots().knots();
    // The last `p` copies of the start and the first `p` of the end; the
    // control point before the start's copies is the curve's point there.
    let first = k.partition_point(|x| *x <= range.0) - p;
    let last = k.partition_point(|x| *x < range.1) + p;
    let Some(start) = first.checked_sub(1) else {
        ogeom_core::ogeom_bail!(
            Construction,
            "a crease's start has no control point before it"
        );
    };
    let mut knots = Vec::with_capacity(last - first + 2);
    knots.push(range.0);
    knots.extend_from_slice(&k[first..last]);
    knots.push(range.1);
    let control: Vec<Point> = raised.control_points()[start..start + last - first + 1 - p]
        .iter()
        .map(|w| w.point())
        .collect();
    ogeom_geom::BSplineCurve::new(ogeom_math::KnotVector::new(knots, p)?, control, tol)
}

/// The ruled bevel between two contacts over one parameter range, each
/// ruling joining the two at one parameter: the two given shared knots and
/// degree, the rulings run on past both contacts by half their length, as
/// a swept bevel's band does, so a later blend reading the bevel's surface
/// past its edges finds it there. The contacts stay exact, the surface
/// linear across.
///
/// # Errors
///
/// As making the two compatible and building the surface.
fn ruled_bevel(
    a: &mut ogeom_geom::BSplineCurve,
    b: &mut ogeom_geom::BSplineCurve,
    tol: Tolerances,
) -> OgeomResult<SurfaceGeometry> {
    crate::bridge::compatible(a, b, tol)?;
    let knots = a.knots().clone();
    let rows = a.control_points().len();
    let mut grid = Vec::with_capacity(2 * rows);
    for (p, q) in a.control_points().iter().zip(b.control_points()) {
        let (p, q) = (p.point(), q.point());
        grid.push(p + (p - q) * BEVEL_RUN_ON);
        grid.push(q + (q - p) * BEVEL_RUN_ON);
    }
    Ok(ogeom_geom::BSplineSurface::new(
        knots,
        ogeom_math::KnotVector::clamped_uniform(1, 2)?,
        &ogeom_math::ControlGrid::new(grid, rows, 2)?,
        tol,
    )?
    .into())
}

/// The ruled band between two closed contacts as a marched band is laid
/// out: `u` across each ruling from the first contact to the second, `v`
/// round the loop on the contacts' shared knots. No run-on: a band round a
/// loop ends on its contacts, as the marched band a later round over it
/// reads does.
///
/// # Errors
///
/// As making the two compatible and building the surface.
fn across_band(
    a: &mut ogeom_geom::BSplineCurve,
    b: &mut ogeom_geom::BSplineCurve,
    tol: Tolerances,
) -> OgeomResult<SurfaceGeometry> {
    crate::bridge::compatible(a, b, tol)?;
    let knots = a.knots().clone();
    let cols = a.control_points().len();
    let mut grid = Vec::with_capacity(2 * cols);
    for p in a.control_points() {
        grid.push(p.point());
    }
    for q in b.control_points() {
        grid.push(q.point());
    }
    Ok(ogeom_geom::BSplineSurface::new(
        ogeom_math::KnotVector::clamped_uniform(1, 2)?,
        knots,
        &ogeom_math::ControlGrid::new(grid, 2, cols)?,
        tol,
    )?
    .into())
}

/// Whether every point of `face`'s mesh stands on `host`'s material side,
/// within `fit` past its surface at the foot of the point: a measure, not
/// a proof, for a host no certify region reads.
///
/// # Errors
///
/// As meshing the face and projecting onto the host.
pub(super) fn measured_within(
    model: &Model,
    face: &Shape,
    host: &super::curved::CurvedFace,
    fit: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let mesh = ogeom_mesh::triangulate_face(
        model,
        face,
        ogeom_mesh::Deflection {
            chord: fit / 4.0,
            ..ogeom_mesh::Deflection::default()
        },
        tol,
    )?;
    let sign = if host.reversed { -1.0 } else { 1.0 };
    for p in &mesh.positions {
        let foot = ogeom_geom::project::project_on_surface(&host.surface, *p, 16, tol)?;
        let (u, v) = foot.parameters;
        let (du, dv) = ogeom_geom::Surface::d1_at(&host.surface, u, v, tol)?;
        let n = du.cross(dv);
        let m = n.magnitude();
        if m <= 0.0 {
            return Ok(false);
        }
        if (*p - foot.point).dot(n / m) * sign > fit {
            return Ok(false);
        }
    }
    Ok(true)
}

/// How many stations a crease's contacts are read at on a marched seat.
const GENERAL_STATIONS: u32 = 64;

/// How far a general crease's open contact runs on past each end, as a
/// share of the chamfer's reach.
const GENERAL_RUN_ON: f64 = 0.05;

/// How far a chamfer reaches off its edge: its longest leg as given.
fn spec_reach(spec: &Chamfer) -> f64 {
    match spec {
        Chamfer::Symmetric(d) => *d,
        Chamfer::Distances {
            on_face, on_other, ..
        } => on_face.max(*on_other),
        Chamfer::Angle { distance, .. } => *distance,
    }
}

/// A crease of an open marched seat between any two faces the curved read
/// holds (a plane and a spline, two curved faces crossing): each contact
/// fitted through the marched bevel's stations at the crease's own
/// parameters, the farthest it stands off its face measured, at the
/// stations and between them, and the ruled bevel between the two. Along
/// with the strip, each contact's distance off its face, which its edges
/// are widened to. `None` where the seat is closed, a station does not
/// seat, or a contact stands off its face past the bands' fit target.
#[allow(clippy::too_many_arguments, reason = "the seat's whole reading")]
pub(super) fn general_crease(
    model: &mut Model,
    solid: &Shape,
    curved: &Curved,
    edge: Id,
    shape: &Shape,
    faces: [usize; 2],
    spec: &Chamfer,
    whole: Option<Point>,
    tol: Tolerances,
) -> OgeomResult<Option<(Strip, [f64; 2])>> {
    // Two planes meet in a straight edge, whose bevel is `bevel_of`'s; and
    // legs the wedge path's own reading refuses (an angle whose line never
    // meets the other face) are not marched.
    if faces
        .iter()
        .all(|f| matches!(curved.faces[*f].surface, SurfaceGeometry::Plane(_)))
        || legs(model, solid, shape, spec, tol).is_err()
    {
        return Ok(None);
    }
    // Each contact fitted through stations, four times as many again where
    // a fit strays past the bands' fit target between them.
    let fit = crate::marched::band_fit_target(tol);
    let mut fitted_contacts = None;
    for count in [
        GENERAL_STATIONS,
        4 * GENERAL_STATIONS,
        16 * GENERAL_STATIONS,
    ] {
        let whole = whole.or_else(|| {
            let c = &curved.edges[&edge];
            (c.ends.0 == c.ends.1).then(|| curved.points[&c.ends.0])
        });
        let Some(found) =
            crate::chamfer_marched::crease_contacts(model, solid, shape, spec, count, whole, tol)?
        else {
            return Ok(None);
        };
        // The seat's faces in the read's order.
        let mut order = [0usize; 2];
        for (k, fi) in faces.iter().enumerate() {
            let Some(j) = found
                .faces
                .iter()
                .position(|f| f.is_same(&curved.faces[*fi].shape))
            else {
                return Ok(None);
            };
            order[k] = j;
        }
        let mut curves = Vec::with_capacity(2);
        let mut widths = [0.0; 2];
        for k in 0..2 {
            let points = &found.points[order[k]];
            let fitted = if found.closed {
                ogeom_geom::fit::fit_points_at_closed(
                    &found.along,
                    points,
                    3,
                    tol.confusion(),
                    tol,
                )?
            } else {
                ogeom_geom::fit::fit_points_at(&found.along, points, 3, tol.confusion(), tol)?
            };
            let curve = fitted.curve;
            let surface = &curved.faces[faces[k]].surface;
            let mut off: f64 = fitted.error;
            for i in 0..found.along.len() - 1 {
                let t = f64::midpoint(found.along[i], found.along[i + 1]);
                let p = curve.point_at(t, tol)?;
                let foot = ogeom_geom::project::project_on_surface(surface, p, 16, tol)?;
                off = off.max(foot.distance);
            }
            if off > fit {
                break;
            }
            widths[k] = off + tol.confusion();
            curves.push(curve);
        }
        if curves.len() == 2 && !found.closed {
            // An open pair runs on a little past each end, so each crosses
            // the next contact where the crease turns at a vertex (two
            // pieces of a fitted seam meeting at a corner): made compatible,
            // each continued over one parameter span, about a share of the
            // reach on the faster; how far the run-on stands off its face
            // counts in its width.
            let (mut a, mut b) = (curves[0].clone(), curves[1].clone());
            crate::bridge::compatible(&mut a, &mut b, tol)?;
            for at_end in [false, true] {
                let (lo, hi) = a.domain();
                let t = if at_end { hi } else { lo };
                let (va, vb) = (a.d1_at(t, tol)?.magnitude(), b.d1_at(t, tol)?.magnitude());
                let span = GENERAL_RUN_ON * spec_reach(spec) / va.max(vb);
                a = a.extended(at_end, span * va, 2, tol)?;
                b = b.extended(at_end, span * vb, 2, tol)?;
            }
            // Each run-on's span is solved to its length, so the two agree
            // only to the solve: the second takes the first's knots, which
            // moves it by as little.
            let (ka, kb) = (a.knots().knots(), b.knots().knots());
            if ka.len() != kb.len()
                || ka.iter().zip(kb).any(|(x, y)| (x - y).abs() > 1e-6)
                || a.control_points().len() != b.control_points().len()
            {
                return Ok(None);
            }
            b = ogeom_geom::BSplineCurve::new(
                a.knots().clone(),
                b.control_points().iter().map(|w| w.point()).collect(),
                tol,
            )?;
            for (k, c) in [&a, &b].into_iter().enumerate() {
                let surface = &curved.faces[faces[k]].surface;
                let (lo, hi) = c.domain();
                for t in [lo, hi] {
                    let p = c.point_at(t, tol)?;
                    let foot = ogeom_geom::project::project_on_surface(surface, p, 16, tol)?;
                    widths[k] = widths[k].max(foot.distance + tol.confusion());
                }
            }
            if widths.iter().any(|w| *w > fit) {
                return Ok(None);
            }
            curves = vec![a, b];
        }
        if curves.len() == 2 {
            fitted_contacts = Some((curves, widths));
            break;
        }
    }
    let Some((curves, widths)) = fitted_contacts else {
        return Ok(None);
    };
    let (mut a, mut b) = (curves[0].clone(), curves[1].clone());
    let closed = a
        .point_at(a.domain().0, tol)?
        .distance(a.point_at(a.domain().1, tol)?)
        <= tol.confusion();
    let bevel = if closed {
        across_band(&mut a, &mut b, tol)?
    } else {
        ruled_bevel(&mut a, &mut b, tol)?
    };
    Ok(Some((
        Strip {
            edge,
            faces,
            contacts: [Curve::from(a), Curve::from(b)],
            surface: bevel,
            section: Section::Line,
        },
        widths,
    )))
}

/// The cone (a cylinder where the radii agree) about `frame`'s axis
/// through the ring of `radius` in its plane and the ring `(rho, dz)`
/// `dz` along its axis, its chart kept clear of its apex.
fn coaxial_cone(
    frame: Frame,
    radius: f64,
    (rho, dz): (f64, f64),
    tol: Tolerances,
) -> OgeomResult<SurfaceGeometry> {
    let slope = (rho - radius) / dz;
    let rows = (dz.min(0.0), dz.max(0.0));
    let depth = dz.abs();
    if slope.atan().abs() <= tol.angular() {
        let pad = 0.1 * depth;
        return Ok(ogeom_geom::CylinderSurface::new(
            ogeom_math::Cylinder::new(frame, radius, tol)?,
            (rows.0 - pad, rows.1 + pad),
        )?
        .into());
    }
    let apex = -radius / slope;
    let clear = (apex - rows.0).abs().min((apex - rows.1).abs());
    let pad = (0.1 * depth).min(0.5 * clear);
    let cone = ogeom_math::Cone::new(frame, radius, slope.atan(), tol)?;
    Ok(ogeom_geom::ConeSurface::new(cone, (rows.0 - pad, rows.1 + pad))?.into())
}
