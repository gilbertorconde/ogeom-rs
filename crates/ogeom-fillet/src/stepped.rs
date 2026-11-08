//! A concave edge's round whose ball rests on what lies beyond the narrow
//! face beside it: a wall standing on a staircase, a flat in an inside
//! corner.
//!
//! The seam `E` between a face `F1` and a narrow strip `F2` is straight.
//! Beyond `F2` the solid runs on in strips `F2, F3, ..., Fk`, each bounded
//! by the far edge of the one before and its own far edge `E2, ..., Ek`,
//! every one a line parallel to `E` over the whole length of the crease.
//! The ball touching `F1` may clear every far edge but a convex `Ek`,
//! resting on it, or, past a concave `Ek`, touch the plane face beyond it
//! as well. Its centre then stands a radius off `F1`'s plane and a radius
//! off `Ek`'s line or the face's plane, which run parallel to `E`, so it
//! runs on a line parallel to `E` and every section square to the crease
//! is the same. The round is exactly the cylinder of the radius about that
//! line.
//!
//! The fill is the prism of one section: from the crease to the touch on
//! `F1`, the ball's arc to its other contact, and back to the crease along
//! the solid's own section through the far edges in between. Fused on, its
//! faces on the strips lie on them, so the strips the ball spans go into
//! the fill and the seam's round reaches over them.
//!
//! An edge is held only where it bounds no face of the seam (`k` of three
//! or more; a ball on `F2`'s own far edge is the rolled round of
//! [`crate::rolled`]), and a face only where the touch on it stands on it.
//! The walk stops at a face it has met, `F1` among them. The section must
//! be a simple polygon lying outside the solid: every side of it on a face
//! has that face's outward normal pointing into it, so the fill adds only
//! what lies off the faces it stands on (the walk does not reach round the
//! solid's far side). The crease ends on faces square to it, so each end
//! of the prism lies in one. Every ball is in free space, proved for each face of the solid:
//! behind `F1`'s plane,
//! behind the plane through the other contact square to the ball's radius
//! (the strips either side of a held edge and the end faces near it), at
//! a distance read in closed form where the face is square to the crease
//! (an end face: the centre's foot on its plane stands still), and
//! otherwise read along the run, an interval halved until the centre's
//! travel fits under its clearance.

use crate::consumed::{seat_of, segment};
use crate::rolled::{Flat, enters, face_distance, flat};
use crate::support::{Seat, edge_curve, same_occurrence};
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{CircleCurve, Curve, Curve3d as _, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Filter, Model, Shape, ShapeType, explore, explore_unique};

/// How many equal stretches of the run a face's clearance is first read
/// over before halving.
const STATIONS: usize = 64;

/// How many times an interval of the run is halved, at most, to show that
/// the ball stays clear of a face between two stations.
const HALVINGS: u32 = 24;

/// A stepped run: the ball's section at the crease's start and the run's
/// length along the crease.
pub(crate) struct Stepped {
    /// The crease's start.
    crease: Point,
    /// The touch on `F1` at the crease's start.
    touch: Point,
    /// The ball's centre at the crease's start.
    centre: Point,
    /// The section's corners on the solid from the contact on the held
    /// far edge back to the crease, exclusive of both ends' own points:
    /// the far edges in between, nearest the held edge first.
    between: Vec<Point>,
    /// The contact on the held far edge at the crease's start.
    contact: Point,
    /// The crease's unit direction and length.
    along: Vector,
    length: f64,
}

/// The run of a ball of `radius` rounding the concave straight `edge`
/// while it rests on the far edge of a strip beyond the narrow face beside
/// it; `None` where the edge is not of that kind or no such ball is in
/// free space all along.
///
/// # Errors
///
/// Whatever reading the solid reports.
#[allow(
    clippy::too_many_lines,
    reason = "one run, walked and proved end to end"
)]
pub(crate) fn stepped_run(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    seat: &Seat,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Stepped>> {
    let conf = tol.confusion();
    if seat.convex {
        return Ok(None);
    }
    let (curve, range) = edge_curve(model, edge, tol)?;
    if !matches!(curve, Curve::Line(_)) {
        return Ok(None);
    }
    let e0 = curve.point_at(range.0, tol)?;
    let e1 = curve.point_at(range.1, tol)?;
    let length = e0.distance(e1);
    if length <= conf {
        return Ok(None);
    }
    let t = (e1 - e0) / length;
    // Where along the crease a point stands, from its start.
    let at = |p: Point| (p - e0).dot(t);

    let faces = explore_unique(model, solid, ShapeType::Face)?;
    let mut flats = Vec::with_capacity(faces.len());
    for face in &faces {
        flats.push((face.clone(), flat(model, face, tol)?));
    }
    // Each end of the crease on a face square to it.
    for end in [0.0, length] {
        let square = flats.iter().any(|(_, found)| {
            found.as_ref().is_some_and(|f| {
                f.normal.cross(t).magnitude() <= tol.angular() && (at(f.origin) - end).abs() <= conf
            })
        });
        if !square {
            return Ok(None);
        }
    }

    for narrow in 0..2 {
        let f1 = &seat.faces[1 - narrow];
        let n1 = seat.normals[1 - narrow];
        let Some(f1_flat) = flat(model, f1, tol)? else {
            continue;
        };
        let into = -seat.leg(narrow, tol)?;
        let mut face = seat.faces[narrow].clone();
        let mut near = edge.clone();
        let mut corners: Vec<Point> = Vec::new();
        // Each strip's outward normal, from `F2` on.
        let mut strips: Vec<Vector> = Vec::new();
        let mut walked: Vec<Shape> = vec![f1.clone()];
        for step in 0..faces.len() {
            if walked.iter().any(|w| w.is_same(&face)) {
                break;
            }
            let Some(strip) = flat(model, &face, tol)? else {
                break;
            };
            strips.push(strip.normal);
            walked.push(face.clone());
            // The strip's far edge: its one other edge parallel to the
            // crease, off the near edge's line, over the crease's length.
            let mut far: Option<(Shape, Point)> = None;
            let mut count = 0;
            for g in explore_unique(model, &face, ShapeType::Edge)? {
                if same_occurrence(model, &g, &near, tol) {
                    continue;
                }
                let Some([p, q]) = segment(model, &g, tol)? else {
                    continue;
                };
                let d = q - p;
                if d.magnitude() <= conf || d.cross(t).magnitude() > tol.angular() * d.magnitude() {
                    continue;
                }
                count += 1;
                let (ap, aq) = (at(p), at(q));
                let spans = (ap.abs() <= conf && (aq - length).abs() <= conf)
                    || (aq.abs() <= conf && (ap - length).abs() <= conf);
                if spans {
                    far = Some((g, if ap.abs() <= conf { p } else { q }));
                }
            }
            let (Some((far_edge, foot)), 1) = (far, count) else {
                break;
            };
            let Some(far_seat) = seat_of(model, solid, &far_edge, tol)? else {
                break;
            };
            let Some(k) = far_seat.faces.iter().position(|f| !f.is_same(&face)) else {
                break;
            };
            if far_seat.faces.iter().all(|f| !f.is_same(&face)) {
                break;
            }
            if far_seat.convex
                && step >= 1
                && let Some(run) = held(
                    &flats,
                    (e0, t, length),
                    (f1, &f1_flat, n1, into),
                    (foot, &far_seat),
                    (&corners, &strips),
                    model,
                    edge,
                    radius,
                    tol,
                )?
            {
                return Ok(Some(run));
            }
            corners.push(foot);
            near = far_edge;
            face = far_seat.faces[k].clone();
            // Past a concave far edge, the ball may touch the face beyond
            // as well as `F1`: the face the strips run into.
            if !far_seat.convex
                && let Some(run) = face_held(
                    &flats,
                    (e0, t, length),
                    (f1, &f1_flat, n1, into),
                    (foot, &far_seat, k),
                    (&corners, &strips),
                    model,
                    edge,
                    radius,
                    tol,
                )?
            {
                return Ok(Some(run));
            }
        }
    }
    Ok(None)
}

/// The run where the ball touching `F1` rests on the convex far edge
/// through `foot`, if that ball is in free space all along: its centre a
/// radius off `F1`'s plane and off the edge's line on the narrow face's
/// side of the crease, its radius to the contact between the edge's two
/// faces' normals.
#[allow(clippy::too_many_arguments, reason = "one proof, all its data")]
fn held(
    flats: &[(Shape, Option<Flat>)],
    (e0, t, length): (Point, Vector, f64),
    (f1, f1_flat, n1, into): (&Shape, &Flat, Vector, Vector),
    (foot, far_seat): (Point, &Seat),
    (corners, strips): (&[Point], &[Vector]),
    model: &Model,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Stepped>> {
    let conf = tol.confusion();
    // In the section through the crease's start: the centre is
    // `e0 + n1 r + b x`, a radius from the foot.
    let b = t.cross(n1);
    let d = e0 + n1 * radius - foot;
    let half = d.dot(b);
    let disc = half * half - (d.dot(d) - radius * radius);
    if disc < 0.0 {
        return Ok(None);
    }
    let [na, nb] = far_seat.normals;
    let between = |u: Vector| {
        // `u = alpha na + beta nb`, both at least zero, in the section.
        let det = na.cross(nb).dot(t);
        if det.abs() <= tol.angular() {
            return false;
        }
        let alpha = u.cross(nb).dot(t) / det;
        let beta = na.cross(u).dot(t) / det;
        alpha >= -tol.angular() && beta >= -tol.angular()
    };
    for x in [-half - disc.sqrt(), -half + disc.sqrt()] {
        let centre = e0 + n1 * radius + b * x;
        if (centre - e0).dot(into) <= 0.0 {
            continue;
        }
        let reach = centre - foot;
        if !between(reach / radius) {
            continue;
        }
        let touch = centre - n1 * radius;
        // The touch stays on `F1`: no edge of it but the crease enters the
        // stretch from the crease to the line of contact.
        let leg = [e0, e0 + t * length, touch + t * length, touch];
        let mut off_face = false;
        for g in explore_unique(model, f1, ShapeType::Edge)? {
            if same_occurrence(model, &g, edge, tol) {
                continue;
            }
            let Some(side) = segment(model, &g, tol)? else {
                return Ok(None);
            };
            if enters(&leg, f1_flat.normal, side, conf) {
                off_face = true;
                break;
            }
        }
        if off_face {
            continue;
        }
        let mut normals = vec![Some(n1), None];
        normals.extend(strips.iter().rev().map(|n| Some(*n)));
        if !fills_outside(&ring(e0, touch, foot, corners), &normals, t, tol) {
            continue;
        }
        if !clear(
            model,
            flats,
            (centre, t, length),
            (e0, n1),
            (foot, reach),
            radius,
            tol,
        )? {
            continue;
        }
        return Ok(Some(Stepped {
            crease: e0,
            touch,
            centre,
            between: corners.iter().rev().copied().collect(),
            contact: foot,
            along: t,
            length,
        }));
    }
    Ok(None)
}

/// The run where the ball touching `F1` touches the plane face beyond the
/// concave far edge through `foot` as well (`far_seat.faces[k]`), if that
/// ball is in free space all along: its centre a radius off both planes on
/// the narrow face's side of the crease, each touch on its own face.
#[allow(clippy::too_many_arguments, reason = "one proof, all its data")]
fn face_held(
    flats: &[(Shape, Option<Flat>)],
    (e0, t, length): (Point, Vector, f64),
    (f1, f1_flat, n1, into): (&Shape, &Flat, Vector, Vector),
    (foot, far_seat, k): (Point, &Seat, usize),
    (corners, strips): (&[Point], &[Vector]),
    model: &Model,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Stepped>> {
    let conf = tol.confusion();
    let beyond = &far_seat.faces[k];
    let Some(beyond_flat) = flat(model, beyond, tol)? else {
        return Ok(None);
    };
    let n3 = far_seat.normals[k];
    // In the section through the crease's start: the centre is
    // `e0 + n1 r + b x`, a radius off the beyond face's plane.
    let b = t.cross(n1);
    let rate = n3.dot(b);
    if rate.abs() <= tol.angular() {
        return Ok(None);
    }
    let x = (radius - n3.dot(e0 + n1 * radius - foot)) / rate;
    let centre = e0 + n1 * radius + b * x;
    if (centre - e0).dot(into) <= 0.0 {
        return Ok(None);
    }
    let touch = centre - n1 * radius;
    let contact = centre - n3 * radius;
    // Each touch stands on its own face past the edge it leaves from: no
    // edge of the face but that one enters the stretch from the edge to the
    // line of contact.
    let into3 = -far_seat.leg(k, tol)?;
    if (contact - foot).dot(into3) <= conf {
        return Ok(None);
    }
    let far_edge = {
        let mut found = None;
        for g in explore_unique(model, beyond, ShapeType::Edge)? {
            if let Some([p, q]) = segment(model, &g, tol)?
                && (p - foot).cross(t).magnitude() <= conf
                && (q - foot).cross(t).magnitude() <= conf
            {
                found = Some(g);
                break;
            }
        }
        found
    };
    for (face, plane, from, to, skip) in [
        (f1, f1_flat, e0, touch, Some(edge.clone())),
        (beyond, &beyond_flat, foot, contact, far_edge),
    ] {
        let leg = [from, from + t * length, to + t * length, to];
        for g in explore_unique(model, face, ShapeType::Edge)? {
            if skip
                .as_ref()
                .is_some_and(|s| same_occurrence(model, &g, s, tol))
            {
                continue;
            }
            let Some(side) = segment(model, &g, tol)? else {
                return Ok(None);
            };
            if enters(&leg, plane.normal, side, conf) {
                return Ok(None);
            }
        }
    }
    let mut normals = vec![Some(n1), None, Some(n3)];
    normals.extend(strips.iter().rev().map(|n| Some(*n)));
    if !fills_outside(&ring(e0, touch, contact, corners), &normals, t, tol) {
        return Ok(None);
    }
    if !clear(
        model,
        flats,
        (centre, t, length),
        (e0, n1),
        (contact, n3 * radius),
        radius,
        tol,
    )? {
        return Ok(None);
    }
    Ok(Some(Stepped {
        crease: e0,
        touch,
        centre,
        between: corners.iter().rev().copied().collect(),
        contact,
        along: t,
        length,
    }))
}

/// The section's ring at the crease's start: the crease, the touch on
/// `F1`, the ball's other contact, and the corners back over the strips.
fn ring(crease: Point, touch: Point, contact: Point, corners: &[Point]) -> Vec<Point> {
    let mut ring = vec![crease, touch, contact];
    ring.extend(corners.iter().rev().copied());
    ring
}

/// Whether the section `ring`, its arc taken as the chord from its second
/// corner to its third, is a simple polygon square to `t` lying outside
/// the solid: each side on a face of the solid (`Some` of that face's
/// outward normal, from its corner to the next) has that normal pointing
/// into the polygon, so the fill adds only what lies off the solid's faces
/// there.
fn fills_outside(ring: &[Point], normals: &[Option<Vector>], t: Vector, tol: Tolerances) -> bool {
    let n = ring.len();
    if n < 3 || normals.len() != n {
        return false;
    }
    let conf = tol.confusion();
    let turn = (0..n).fold(Vector::ZERO, |sum, i| {
        let (p, q) = (ring[i], ring[(i + 1) % n]);
        sum + (p - ring[0]).cross(q - ring[0])
    });
    let sense = turn.dot(t);
    if sense.abs() <= conf * conf {
        return false;
    }
    let sense = sense.signum();
    for i in 0..n {
        let d = ring[(i + 1) % n] - ring[i];
        if d.magnitude() <= conf {
            return false;
        }
        if let Some(normal) = normals[i]
            && normal.dot(t.cross(d)) * sense <= tol.angular() * d.magnitude()
        {
            return false;
        }
    }
    // No two sides that do not share a corner meet.
    let (u, v) = {
        let d = ring[1] - ring[0];
        let u = d / d.magnitude();
        (u, t.cross(u))
    };
    let flat = |p: Point| ((p - ring[0]).dot(u), (p - ring[0]).dot(v));
    let cross = |a: (f64, f64), b: (f64, f64), c: (f64, f64)| {
        (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
    };
    for i in 0..n {
        for j in i + 1..n {
            if j == i + 1 || (i == 0 && j == n - 1) {
                continue;
            }
            let (a, b) = (flat(ring[i]), flat(ring[(i + 1) % n]));
            let (c, d) = (flat(ring[j]), flat(ring[(j + 1) % n]));
            let (d1, d2) = (cross(a, b, c), cross(a, b, d));
            let (d3, d4) = (cross(c, d, a), cross(c, d, b));
            if d1 * d2 <= 0.0 && d3 * d4 <= 0.0 {
                return false;
            }
        }
    }
    true
}

/// Whether every ball of the run, its centre on `centre + t s` for `s` in
/// `[0, length]`, lies in free space: no face of the solid comes within
/// its radius.
fn clear(
    model: &Model,
    flats: &[(Shape, Option<Flat>)],
    (centre, t, length): (Point, Vector, f64),
    (e0, n1): (Point, Vector),
    (foot, reach): (Point, Vector),
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let conf = tol.confusion();
    let at = |s: f64| centre + t * s;
    #[allow(clippy::cast_precision_loss)]
    let station = |i: usize| length * (i as f64) / (STATIONS as f64);
    for (face, found) in flats {
        let Some(found) = found else {
            // A curved face is cleared only by its box.
            let bounds = ogeom_algo::shape_bounds(model, face, tol)?;
            #[allow(clippy::cast_precision_loss)]
            let step = length / (STATIONS as f64);
            for i in 0..=STATIONS {
                if bounds.distance_to(at(station(i))) <= radius + step {
                    return Ok(false);
                }
            }
            continue;
        };
        let vertices: Vec<Point> = found.sides.iter().flatten().copied().collect();
        // Behind the touched face's plane, or behind the plane through the
        // contact square to the ball's radius, which the run keeps.
        if vertices.iter().all(|v| (*v - e0).dot(n1) <= conf)
            || vertices
                .iter()
                .all(|v| (*v - foot).dot(reach) <= conf * radius)
        {
            continue;
        }
        // Square to the run: the centre's foot on the face's plane stands
        // still, and the distance is least at the plane.
        if found.normal.cross(t).magnitude() <= tol.angular() {
            let height = (centre - found.origin).dot(found.normal);
            let on_plane = centre - found.normal * height;
            if face_distance(found, on_plane) < radius - conf {
                return Ok(false);
            }
            continue;
        }
        let clearance = |s: f64| face_distance(found, at(s)) - radius;
        let mut intervals: Vec<(f64, f64, f64, f64, u32)> = (0..STATIONS)
            .map(|i| {
                let (a, b) = (station(i), station(i + 1));
                (a, b, clearance(a), clearance(b), 0)
            })
            .collect();
        while let Some((a, b, ca, cb, depth)) = intervals.pop() {
            if ca < -conf || cb < -conf {
                return Ok(false);
            }
            if ca.min(cb) >= (b - a) * 0.5 {
                continue;
            }
            if depth >= HALVINGS {
                return Ok(false);
            }
            let m = (a + b) * 0.5;
            let cm = clearance(m);
            intervals.push((a, m, ca, cm, depth + 1));
            intervals.push((m, b, cm, cb, depth + 1));
        }
    }
    Ok(true)
}

/// Round `edge` with the stepped run's fill: the prism of its section,
/// fused on.
///
/// # Errors
///
/// As the boolean that fuses it on reports.
pub(crate) fn stepped_band(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    run: &Stepped,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let along = run.along;
    let mut ring: Vec<Point> = vec![run.crease, run.touch, run.contact];
    ring.extend(run.between.iter().copied());
    let vertices: Vec<Shape> = ring
        .iter()
        .map(|p| ogeom_algo::make_vertex(model, *p).shape)
        .collect();
    let line = |model: &mut Model, i: usize, j: usize| -> OgeomResult<Shape> {
        crate::support::segment_between(
            model,
            (&vertices[i], ring[i]),
            (&vertices[j], ring[j]),
            tol,
        )
    };
    // The ball's arc from the touch to the contact, the short way round,
    // which faces the crease.
    let x_dir = (run.touch - run.centre) / radius;
    let y_dir = along.cross(x_dir);
    let to = run.contact - run.centre;
    let angle = to.dot(y_dir).atan2(to.dot(x_dir));
    let normal = if angle >= 0.0 { along } else { -along };
    let frame = Frame::new(
        run.centre,
        Direction::new(normal, tol)?,
        Direction::new(x_dir, tol)?,
        tol,
    )?;
    let circle: Curve = CircleCurve::new(Circle::new(frame, radius, tol)?).into();
    let arc = ogeom_algo::make_edge_between(
        model,
        circle,
        (0.0, angle.abs()),
        &vertices[1],
        &vertices[2],
        tol,
    )?
    .shape;
    let mut edges = vec![line(model, 0, 1)?, arc];
    for i in 2..ring.len() {
        edges.push(line(model, i, (i + 1) % ring.len())?);
    }
    // Wound counter-clockwise about the section plane's normal.
    let turn = (0..ring.len()).fold(Vector::ZERO, |sum, i| {
        let (p, q) = (ring[i], ring[(i + 1) % ring.len()]);
        sum + (p - ring[0]).cross(q - ring[0])
    });
    let normal = if turn.dot(along) >= 0.0 {
        along
    } else {
        -along
    };
    let reach = ring
        .iter()
        .fold(radius, |m, p| m.max(p.distance(run.crease)))
        * 4.0;
    let plane: SurfaceGeometry = PlaneSurface::over(
        Plane::through(run.crease, Direction::new(normal, tol)?),
        (-reach, reach),
        (-reach, reach),
    )?
    .into();
    let loop_edges = if turn.dot(along) >= 0.0 {
        edges
    } else {
        edges.iter().rev().map(Shape::reversed).collect()
    };
    let section = ogeom_algo::make_face_with_pcurves(model, plane, &[loop_edges], tol)?.shape;
    let wedge = ogeom_algo::make_prism(model, &section, along * run.length, tol)?.shape;
    let wedge = explore(model, &wedge, Filter::OfType(ShapeType::Solid))?
        .into_iter()
        .next()
        .unwrap_or(wedge);
    crate::support::apply_wedge_solid(model, solid, Some(edge), &wedge, true, tol)
}
