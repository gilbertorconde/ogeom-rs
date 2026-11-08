//! A round of a rim that takes a narrow ring off whole.
//!
//! Where a ring beside a convex rim is narrower than the ball's setback (a
//! chamfered or small-rounded boss whose rim is rounded larger), no ball
//! touches both of the rim's faces inside the material. The ball that
//! rounds the rim then touches the cap and the wall either side of the
//! ring, and where the ring lies wholly in the corner that ball takes off,
//! the round is the blend of cap and wall, seated on the circle where the
//! cap's plane meets the wall, and the ring goes with the corner: the
//! revolved counterpart of [`crate::consumed`].

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, SurfaceGeometry};
use ogeom_math::{Aabb, Direction, Point, Vector};
use ogeom_topo::{Filter, Model, NodeData, Orientation, Shape, ShapeType, explore, explore_unique};

use crate::support::{RevolvedSeat, edge_curve, same_occurrence};

/// A point of the rim's half-plane section: its distance from the axis
/// and its height along the cap's outward normal from the cap's plane.
type Section = (f64, f64);

/// A face's surface, placed in the world, as the ring and the seat read it.
enum Placed {
    /// A plane through `on` with outward unit normal `normal`.
    Plane { on: Point, normal: Vector },
    /// A cylinder or cone about the axis through `origin` along `axis`, its
    /// radius `radius + growth z` at `z` along `axis` from `origin`; `sign`
    /// is `+1` with the material inside it.
    Wall {
        origin: Point,
        axis: Vector,
        radius: f64,
        growth: f64,
        sign: f64,
    },
    /// A torus about the axis through `origin` along `axis`.
    Torus {
        origin: Point,
        axis: Vector,
        major: f64,
        minor: f64,
    },
    /// Anything else.
    Other,
}

fn placed(model: &Model, face: &Shape) -> OgeomResult<Placed> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(Placed::Other);
    };
    let at = face.transform(model.datums())?;
    let reversed = face.orientation() == Orientation::Reversed;
    let sign = if reversed { -1.0 } else { 1.0 };
    Ok(match model.geometry().surface(data.surface) {
        Some(SurfaceGeometry::Plane(p)) => Placed::Plane {
            on: at.apply(p.plane().origin()),
            normal: at.apply_vector(p.plane().normal().vector()) * sign,
        },
        Some(SurfaceGeometry::Cylinder(c)) => {
            let cyl = c.cylinder();
            Placed::Wall {
                origin: at.apply(cyl.frame().origin()),
                axis: at.apply_vector(cyl.frame().z().vector()),
                radius: cyl.radius(),
                growth: 0.0,
                sign,
            }
        }
        Some(SurfaceGeometry::Cone(c)) => {
            let cone = c.cone();
            Placed::Wall {
                origin: at.apply(cone.frame().origin()),
                axis: at.apply_vector(cone.frame().z().vector()),
                radius: cone.reference_radius(),
                growth: cone.half_angle().tan(),
                sign,
            }
        }
        Some(SurfaceGeometry::Torus(t)) => {
            let torus = t.torus();
            Placed::Torus {
                origin: at.apply(torus.frame().origin()),
                axis: at.apply_vector(torus.frame().z().vector()),
                major: torus.major_radius(),
                minor: torus.minor_radius(),
            }
        }
        _ => Placed::Other,
    })
}

/// The faces of `solid` holding `edge`.
fn holders(model: &Model, solid: &Shape, edge: &Shape, tol: Tolerances) -> OgeomResult<Vec<Shape>> {
    let mut out = Vec::new();
    for face in explore_unique(model, solid, ShapeType::Face)? {
        if explore(model, &face, Filter::OfType(ShapeType::Edge))?
            .iter()
            .any(|e| same_occurrence(model, e, edge, tol))
        {
            out.push(face);
        }
    }
    Ok(out)
}

/// A face's edges other than its seams: a seam is the one edge used twice
/// in the face's boundary.
fn bounds(model: &Model, face: &Shape) -> OgeomResult<Vec<Shape>> {
    let all = explore(model, face, Filter::OfType(ShapeType::Edge))?;
    let mut out: Vec<Shape> = Vec::new();
    for e in &all {
        if all.iter().filter(|f| f.node() == e.node()).count() > 1 {
            continue;
        }
        if !out.iter().any(|o| o.node() == e.node()) {
            out.push(e.clone());
        }
    }
    Ok(out)
}

/// How far round its circle a circular edge turns.
fn sweep_of(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<f64> {
    let (_, range) = edge_curve(model, edge, tol)?;
    Ok((range.1 - range.0).abs())
}

/// The edges of `edges` on the circle `(centre, axis, radius)`, and the
/// rest.
fn on_circle(
    model: &Model,
    edges: &[Shape],
    (centre, axis, radius): (Point, Vector, f64),
    tol: Tolerances,
) -> OgeomResult<(Vec<Shape>, Vec<Shape>)> {
    let (mut on, mut off) = (Vec::new(), Vec::new());
    for e in edges {
        match circle_of(model, e, tol)? {
            Some((c, a, r))
                if c.distance(centre) <= tol.confusion()
                    && (r - radius).abs() <= tol.confusion()
                    && a.cross(axis).magnitude() <= tol.angular() =>
            {
                on.push(e.clone());
            }
            _ => off.push(e.clone()),
        }
    }
    Ok((on, off))
}

/// The one face of `solid` other than `strip` holding every edge of `run`.
fn beside(
    model: &Model,
    solid: &Shape,
    run: &[Shape],
    strip: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let mut found: Option<Shape> = None;
    for e in run {
        let others: Vec<Shape> = holders(model, solid, e, tol)?
            .into_iter()
            .filter(|f| !f.is_same(strip))
            .collect();
        let [face] = others.as_slice() else {
            return Ok(None);
        };
        match &found {
            Some(held) if !held.is_same(face) => return Ok(None),
            Some(_) => {}
            None => found = Some(face.clone()),
        }
    }
    Ok(found)
}

/// Whether the arcs of `run` turn their circle's whole turn.
fn whole_turn(model: &Model, run: &[Shape], tol: Tolerances) -> OgeomResult<bool> {
    // Each arc's sweep is read within the angular tolerance.
    let (mut turned, mut allowed) = (0.0, 0.0);
    for e in run {
        turned += sweep_of(model, e, tol)?;
        allowed += tol.angular();
    }
    Ok((turned - core::f64::consts::TAU).abs() <= allowed)
}

/// A circular edge's centre, unit axis and radius.
fn circle_of(
    model: &Model,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<(Point, Vector, f64)>> {
    let Ok((curve, _)) = edge_curve(model, edge, tol) else {
        return Ok(None);
    };
    let Curve::Circle(c) = curve else {
        return Ok(None);
    };
    let circle = c.circle();
    Ok(Some((
        circle.centre(),
        circle.frame().z().vector(),
        circle.radius(),
    )))
}

/// Whether the line through `origin` along `along` is the axis through
/// `centre` along `axis`.
fn coaxial(origin: Point, along: Vector, centre: Point, axis: Vector, tol: Tolerances) -> bool {
    let along = along / along.magnitude();
    let off = centre - origin;
    along.cross(axis).magnitude() <= tol.angular()
        && (off - along * off.dot(along)).magnitude() <= tol.confusion()
}

/// The seat of the cap and wall either side of a narrow ring beside `edge`,
/// where the ball of `radius` rolling between them takes the ring off
/// whole; `None` where it does not.
///
/// Every condition is read in closed form, in the rim's half-plane section.
/// The ring `F2`, a coaxial cone or torus, lies between the circle of `E`
/// and another coaxial circle, its arcs on each turning the whole turn, those
/// on `E`'s on the face `F1` and the others on `F3`, one
/// of `F1` and `F3` the cap (a plane square to the axis) and the other the
/// wall (a coaxial cylinder or cone) running away from the cap's outside.
/// With `V` the circle where the cap's plane meets the wall, `c` the ball's
/// centre a radius in from both and `T1`, `T3` its contacts:
///
/// - `E` stands on its face's leg from `V` to its contact, short of it,
///   and `E'` likewise on the other's;
/// - the ring's section, the segment from `E` to `E'` on a cone, or on a
///   torus the arc from `E` to `E'` facing `V`, the torus touching both
///   legs there, stays more than a radius from `c`.
///
/// The region the ring's section closes with the legs then lies in the
/// triangle `V E E'`, whose border stays off the ball, and the ball is not
/// within it, since it touches the legs past `E` and `E'`: so the region
/// lies in the corner the cut takes off. What else of the solid that corner
/// holds is bounded by the cap from `E` to its contact and the wall from
/// `E'` to its own: no edge of those faces but coaxial circles past the
/// contacts, and no other face's bound meets the corner's.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle
/// fails to resolve; a solid this does not read is `None`.
pub(crate) fn consumed_rim(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<RevolvedSeat>> {
    let Some((centre_e, axis, _)) = circle_of(model, edge, tol)? else {
        return Ok(None);
    };
    let faces = holders(model, solid, edge, tol)?;
    let [p, q] = faces.as_slice() else {
        return Ok(None);
    };
    for (strip, f1) in [(q, p), (p, q)] {
        if let Some(seat) = across(model, solid, edge, (centre_e, axis), strip, f1, radius, tol)? {
            return Ok(Some(seat));
        }
    }
    Ok(None)
}

/// The seat across `strip`, `edge` lying between it and `f1`.
#[allow(clippy::too_many_arguments, reason = "one question, all its data")]
#[allow(clippy::too_many_lines, reason = "one seat, read end to end")]
fn across(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    (centre_e, axis): (Point, Vector),
    strip: &Shape,
    f1: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<RevolvedSeat>> {
    let conf = tol.confusion();
    let ring = placed(model, strip)?;
    let strip_axis = match &ring {
        Placed::Wall {
            origin,
            axis: along,
            growth,
            ..
        } if *growth != 0.0 => (*origin, *along),
        Placed::Torus {
            origin,
            axis: along,
            ..
        } => (*origin, *along),
        _ => return Ok(None),
    };
    if !coaxial(strip_axis.0, strip_axis.1, centre_e, axis, tol) {
        return Ok(None);
    }
    // The ring: the arcs of `E`'s circle on one side, of one other
    // coaxial circle on the other, each side the whole turn, seams aside.
    let Some(circle_e) = circle_of(model, edge, tol)? else {
        return Ok(None);
    };
    let edges = bounds(model, strip)?;
    let (near, far) = on_circle(model, &edges, circle_e, tol)?;
    let Some(first_far) = far.first() else {
        return Ok(None);
    };
    let Some(circle_far) = circle_of(model, first_far, tol)? else {
        return Ok(None);
    };
    let (far, stray) = on_circle(model, &far, circle_far, tol)?;
    if !stray.is_empty()
        || !coaxial(centre_e, axis, circle_far.0, circle_far.1, tol)
        || !near.iter().any(|b| same_occurrence(model, b, edge, tol))
        || !whole_turn(model, &near, tol)?
        || !whole_turn(model, &far, tol)?
    {
        return Ok(None);
    }
    if !beside(model, solid, &near, strip, tol)?.is_some_and(|f| f.is_same(f1)) {
        return Ok(None);
    }
    let Some(f3) = beside(model, solid, &far, strip, tol)? else {
        return Ok(None);
    };
    let f3 = &f3;
    if f3.is_same(f1) {
        return Ok(None);
    }
    // One of the outer faces is the cap, the other the wall.
    let (cap_face, wall_face, cap_edge) = match (placed(model, f1)?, placed(model, f3)?) {
        (Placed::Plane { .. }, Placed::Wall { .. }) => (f1, f3, &near),
        (Placed::Wall { .. }, Placed::Plane { .. }) => (f3, f1, &far),
        _ => return Ok(None),
    };
    let wall_edge = if std::ptr::eq(cap_edge, &near) {
        &far
    } else {
        &near
    };
    let (cap_edges, wall_edges) = (cap_edge, wall_edge);
    let (Some(cap_edge), Some(wall_edge)) = (cap_edges.first(), wall_edges.first()) else {
        return Ok(None);
    };
    let Placed::Plane { on, normal: up } = placed(model, cap_face)? else {
        return Ok(None);
    };
    let Placed::Wall {
        origin,
        axis: wall_axis,
        radius: wall_radius,
        growth,
        sign: sigma,
    } = placed(model, wall_face)?
    else {
        return Ok(None);
    };
    let up = up / up.magnitude();
    if up.cross(axis).magnitude() > tol.angular()
        || !coaxial(origin, wall_axis, centre_e, axis, tol)
    {
        return Ok(None);
    }
    // The circle `V`: the axis's point on the cap's plane, and the wall's
    // radius there.
    let centre = centre_e + up * (on - centre_e).dot(up);
    let wall_up = wall_axis / wall_axis.magnitude();
    let rim_radius = growth.mul_add((centre - origin).dot(wall_up), wall_radius);
    if rim_radius <= conf {
        return Ok(None);
    }
    let slope = growth * wall_up.dot(up).signum();
    let section = |p: Point| -> Section {
        let off = p - centre;
        let h = off.dot(up);
        ((off - up * h).magnitude(), h)
    };
    // The rims' sections: each a point, the ring's centre on the axis.
    let point_of = |shape: &Shape| -> OgeomResult<Option<Section>> {
        use ogeom_geom::Curve3d as _;
        let (curve, range) = edge_curve(model, shape, tol)?;
        Ok(Some(section(curve.point_at(range.0, tol)?)))
    };
    let (Some(at_cap), Some(at_wall)) = (point_of(cap_edge)?, point_of(wall_edge)?) else {
        return Ok(None);
    };
    let Some(x_ref) = frame_x(model, edge, tol)? else {
        return Ok(None);
    };
    // The wall runs away from the cap's outside: the rim is convex.
    if at_wall.1 >= -conf {
        return Ok(None);
    }
    let seat = RevolvedSeat {
        centre,
        radius: rim_radius,
        slope,
        up,
        x_ref,
        sigma,
        tau: 1.0,
        cap_face: cap_face.clone(),
        wall_face: wall_face.clone(),
        closes: true,
    };
    let ball = seat.ball(radius);
    if ball.rho <= conf {
        return Ok(None);
    }
    let v: Section = (rim_radius, 0.0);
    let c: Section = (ball.rho, ball.height);
    let t_cap: Section = (ball.rho, 0.0);
    let t_wall: Section = (ball.wall_rho, ball.wall_height);
    let sub = |a: Section, b: Section| (a.0 - b.0, a.1 - b.1);
    let dot = |a: Section, b: Section| a.0.mul_add(b.0, a.1 * b.1);
    let norm = |a: Section| a.0.hypot(a.1);
    // How far along the leg from `V` to `t` a point stands, and how far off
    // the leg's line.
    let along_leg = |p: Section, t: Section| -> (f64, f64, f64) {
        let leg = sub(t, v);
        let reach = norm(leg);
        let unit = (leg.0 / reach, leg.1 / reach);
        let off = sub(p, v);
        let at = dot(off, unit);
        (at, (off.0 - unit.0 * at).hypot(off.1 - unit.1 * at), reach)
    };
    let on_leg = |p: Section, t: Section| {
        let (at, off, reach) = along_leg(p, t);
        off <= conf && at > conf && at < reach - conf
    };
    if !on_leg(at_cap, t_cap) || !on_leg(at_wall, t_wall) {
        return Ok(None);
    }
    // The ring's section, off the ball.
    let clear_of_ball = match ring {
        Placed::Wall { .. } => {
            let d = sub(at_wall, at_cap);
            let t = (dot(sub(c, at_cap), d) / dot(d, d)).clamp(0.0, 1.0);
            norm(sub(c, (at_cap.0 + d.0 * t, at_cap.1 + d.1 * t)))
        }
        Placed::Torus {
            origin,
            major,
            minor,
            ..
        } => {
            let o: Section = (major, (origin - centre).dot(up));
            // The torus touches both legs at the rims: its section's centre
            // a minor radius from each leg's line, the rims on its circle.
            let line_off = |t: Section, p: Section| along_leg(p, t).1;
            if (line_off(t_cap, o) - minor).abs() > conf
                || (line_off(t_wall, o) - minor).abs() > conf
                || (norm(sub(at_cap, o)) - minor).abs() > conf
                || (norm(sub(at_wall, o)) - minor).abs() > conf
            {
                return Ok(None);
            }
            // The arc facing `V`: the face holds its middle.
            let toward = sub(v, o);
            let mid = (
                o.0 + toward.0 * minor / norm(toward),
                o.1 + toward.1 * minor / norm(toward),
            );
            if !holds(
                model,
                strip,
                (centre, up, x_ref.vector()),
                mid,
                minor * 1e-3,
                tol,
            )? {
                return Ok(None);
            }
            // The arc's point nearest `c`, where the arc holds it; else an
            // end.
            let angle = |p: Section| (p.1 - o.1).atan2(p.0 - o.0);
            let (a, b, m) = (angle(at_cap), angle(at_wall), angle(mid));
            let to_c = sub(c, o);
            let near = angle((o.0 + to_c.0, o.1 + to_c.1));
            let within = |x: f64| {
                let turn = core::f64::consts::TAU;
                let span = (b - a).rem_euclid(turn);
                let forward = (m - a).rem_euclid(turn) <= span;
                let at = (x - a).rem_euclid(turn);
                if forward { at <= span } else { at >= span }
            };
            let ends = norm(sub(c, at_cap)).min(norm(sub(c, at_wall)));
            if within(near) {
                ends.min(norm(to_c) - minor)
            } else {
                ends
            }
        }
        _ => return Ok(None),
    };
    if clear_of_ball <= radius + conf {
        return Ok(None);
    }
    // The cap and the wall hold the band's contacts whole: every other edge
    // of each is a coaxial circle on its leg past the contact.
    for (face, own, t) in [
        (cap_face, cap_edges, t_cap),
        (wall_face, wall_edges, t_wall),
    ] {
        for e in bounds(model, face)? {
            if own.iter().any(|o| same_occurrence(model, &e, o, tol)) {
                continue;
            }
            let Some((c_e, a_e, _)) = circle_of(model, &e, tol)? else {
                return Ok(None);
            };
            if !coaxial(centre, up, c_e, a_e, tol) {
                return Ok(None);
            }
            let Some(p) = point_of(&e)? else {
                return Ok(None);
            };
            let (at, off, reach) = along_leg(p, t);
            if off > conf || at < reach - conf {
                return Ok(None);
            }
        }
    }
    // No other face of the solid comes near the corner the cut takes.
    let corner = {
        let x = seat.x_ref.vector();
        let y = up.cross(x);
        let mut points = Vec::new();
        for (rho, h) in [v, t_cap, t_wall] {
            for (dx, dy) in [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)] {
                points.push(centre + up * h + x * (rho * dx) + y * (rho * dy));
            }
        }
        Aabb::of_points(&points)
    };
    for face in explore_unique(model, solid, ShapeType::Face)? {
        if [cap_face, wall_face, strip]
            .iter()
            .any(|k| k.is_same(&face))
        {
            continue;
        }
        let overlap = ogeom_algo::shape_bounds(model, &face, tol)?.intersection(&corner);
        let size = overlap.size();
        if overlap.low().is_some() && size.x > conf && size.y > conf && size.z > conf {
            return Ok(None);
        }
    }
    Ok(Some(seat))
}

/// The `x` of a circular edge's frame, square to its axis.
fn frame_x(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<Option<Direction>> {
    let Ok((Curve::Circle(c), _)) = edge_curve(model, edge, tol) else {
        return Ok(None);
    };
    Ok(Some(c.circle().frame().x()))
}

/// Whether `face`, a ring round the axis through `centre` along `up`,
/// holds the point of section `p`, its mesh read to `chord`: on the
/// half-planes through `x`, a quarter turn on and a half turn on, none
/// outside and one inside, since a seam may run along any one of them.
fn holds(
    model: &Model,
    face: &Shape,
    (centre, up, x): (Point, Vector, Vector),
    p: Section,
    chord: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let deflection = ogeom_mesh::Deflection {
        chord,
        ..ogeom_mesh::Deflection::default()
    };
    let mut inside = false;
    for radial in [x, up.cross(x), -x] {
        let point = centre + radial * p.0 + up * p.1;
        match ogeom_algo::classify_on_face(model, face, point, deflection, tol)? {
            ogeom_algo::Containment::Out => return Ok(false),
            ogeom_algo::Containment::In => inside = true,
            ogeom_algo::Containment::On => {}
        }
    }
    Ok(inside)
}
