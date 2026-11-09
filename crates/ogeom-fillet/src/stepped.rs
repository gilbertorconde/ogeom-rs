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
//! solid's far side).
//!
//! The crease ends on two planar faces, square to it or slanted, every
//! line of the section (the touch on `F1`, each far edge, the contact)
//! running from one end plane to the other. Every section square to the
//! crease is then the same, and the fill is the prism of the section over
//! the span the section's lines and arc reach between the end planes, cut
//! by them. Every ball over that span is in free space, proved for each
//! face of the solid: behind `F1`'s plane, behind the plane through the
//! other contact square to the ball's radius (the strips either side of a
//! held edge), and otherwise for a planar face by its least distance from
//! the centre's segment, read in closed form, so a face the balls touch
//! (an end face where the touch line meets it) is cleared exactly.

use crate::carried::trig_max;
use crate::consumed::{seat_of, segment};
use crate::rolled::{Flat, basis, enters, face_distance, flat, inside};
use crate::support::{Seat, edge_curve, same_occurrence};
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{CircleCurve, Curve, Curve3d as _, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Filter, Model, Shape, ShapeType, explore, explore_unique};

/// How many equal stretches of the run a curved face's box is read over.
const STATIONS: usize = 64;

/// The run's two end planes: each a point on it and its outward unit
/// normal, the crease running from the first to the second along `t`.
#[derive(Debug, Clone, Copy)]
struct Ends {
    e0: Point,
    t: Vector,
    length: f64,
    planes: [(Point, Vector); 2],
    /// Whether both planes stand square to the crease.
    square: bool,
}

impl Ends {
    /// Where the line through `p` along the crease meets end plane `i`: the
    /// `s` of `p + t s`. On square ends, read from the section through the
    /// crease's start, the start plane is at 0 and the end plane at the
    /// crease's length.
    fn meets(&self, i: usize, p: Point) -> f64 {
        if self.square {
            return if i == 0 { 0.0 } else { self.length };
        }
        let (o, m) = self.planes[i];
        (o - p).dot(m) / self.t.dot(m)
    }

    /// Whether `p` lies on end plane `i`.
    fn holds(&self, i: usize, p: Point, conf: f64) -> bool {
        let (o, m) = self.planes[i];
        (p - o).dot(m).abs() <= conf
    }
}

/// A stepped run: the ball's section in the plane square to the crease
/// through its start, the run's ends, and the span of the crease's
/// parameter the section's lines reach between the end planes.
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
    ends: Ends,
    /// The least and greatest `s` at which the section's lines and arc
    /// meet the end planes.
    span: (f64, f64),
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

    let faces = explore_unique(model, solid, ShapeType::Face)?;
    let mut flats = Vec::with_capacity(faces.len());
    for face in &faces {
        flats.push((face.clone(), flat(model, face, tol)?));
    }
    // Each end of the crease on one plane: the planar faces through it
    // that the crease leaves, not along it, all on that plane, facing back
    // along the crease at its start and on along it at its end.
    let mut planes = [(e0, -t), (e1, t)];
    for (i, end) in [e0, e1].into_iter().enumerate() {
        let mut found: Option<Vector> = None;
        for (face, plane) in &flats {
            let Some(plane) = plane else { continue };
            if seat.faces.iter().any(|f| f.is_same(face))
                || face_distance(plane, end) > conf
                || plane.normal.dot(t).abs() <= tol.angular()
            {
                continue;
            }
            match found {
                Some(m) if m.cross(plane.normal).magnitude() > tol.angular() => return Ok(None),
                _ => found = Some(plane.normal),
            }
        }
        let Some(m) = found else {
            return Ok(None);
        };
        let facing = if i == 0 { -m.dot(t) } else { m.dot(t) };
        if facing <= tol.angular() {
            return Ok(None);
        }
        planes[i] = (end, m);
    }
    let square = planes
        .iter()
        .all(|(_, m)| m.cross(t).magnitude() <= tol.angular());
    let ends = Ends {
        e0,
        t,
        length,
        planes,
        square,
    };

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
                // From end plane to end plane; its foot in the section
                // through the crease's start.
                let spans = (ends.holds(0, p, conf) && ends.holds(1, q, conf))
                    || (ends.holds(0, q, conf) && ends.holds(1, p, conf));
                if spans {
                    far = Some((g, p - t * (p - e0).dot(t)));
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
                    &ends,
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
                    &ends,
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
    ends: &Ends,
    (f1, f1_flat, n1, into): (&Shape, &Flat, Vector, Vector),
    (foot, far_seat): (Point, &Seat),
    (corners, strips): (&[Point], &[Vector]),
    model: &Model,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Stepped>> {
    let conf = tol.confusion();
    let (e0, t) = (ends.e0, ends.t);
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
        let leg = stretch(ends, e0, touch);
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
        let section = ring(e0, touch, foot, corners);
        if !fills_outside(&section, &normals, t, tol) {
            continue;
        }
        let span = run_span(ends, &section, centre, radius);
        if !clear(
            model,
            flats,
            (ends, centre, span),
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
            ends: *ends,
            span,
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
    ends: &Ends,
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
    let (e0, t) = (ends.e0, ends.t);
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
        let leg = stretch(ends, from, to);
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
    let section = ring(e0, touch, contact, corners);
    if !fills_outside(&section, &normals, t, tol) {
        return Ok(None);
    }
    let span = run_span(ends, &section, centre, radius);
    if !clear(
        model,
        flats,
        (ends, centre, span),
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
        ends: *ends,
        span,
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
pub(crate) fn fills_outside(
    ring: &[Point],
    normals: &[Option<Vector>],
    t: Vector,
    tol: Tolerances,
) -> bool {
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

/// Whether every ball of the run, its centre on `centre + t s` for `s` over
/// `span`, lies in free space: no face of the solid comes within its
/// radius. A planar face is read at its least distance from the centre's
/// segment, in closed form: to each side of the face, to each end of the
/// segment whose foot lands in the face, and none where the segment
/// pierces it. A curved face is cleared only by its box.
fn clear(
    model: &Model,
    flats: &[(Shape, Option<Flat>)],
    (ends, centre, span): (&Ends, Point, (f64, f64)),
    (e0, n1): (Point, Vector),
    (foot, reach): (Point, Vector),
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let conf = tol.confusion();
    let (from, to) = (centre + ends.t * span.0, centre + ends.t * span.1);
    for (face, found) in flats {
        let Some(found) = found else {
            let bounds = ogeom_algo::shape_bounds(model, face, tol)?;
            #[allow(clippy::cast_precision_loss)]
            let step = (span.1 - span.0) / (STATIONS as f64);
            for i in 0..=STATIONS {
                #[allow(clippy::cast_precision_loss)]
                let at = from + (to - from) * (i as f64 / STATIONS as f64);
                if bounds.distance_to(at) <= radius + step {
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
        if segment_face_distance(found, from, to) < radius - conf {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The least distance between the segment from `a` to `b` and a planar
/// face. The nearest point of the face to the segment lies on a side of the
/// face, or inside it as the foot of a point of the segment, whose height
/// over the plane is least at an end of the stretch of the segment over
/// the face: at an end of the segment, at a side (where the side is no
/// farther), or where the segment crosses the plane.
fn segment_face_distance(face: &Flat, a: Point, b: Point) -> f64 {
    let frame = basis(face.normal);
    let height = |p: Point| (p - face.origin).dot(face.normal);
    let over = |p: Point| inside(&face.sides, p - face.normal * height(p), frame);
    let mut least = face
        .sides
        .iter()
        .map(|&[p, q]| segment_segment_distance((a, b), (p, q)))
        .fold(f64::INFINITY, f64::min);
    for end in [a, b] {
        if over(end) {
            least = least.min(height(end).abs());
        }
    }
    let (ha, hb) = (height(a), height(b));
    if ha * hb < 0.0 && over(a + (b - a) * (ha / (ha - hb))) {
        least = 0.0;
    }
    least
}

/// The least distance between two segments: over the square of their
/// parameters, at the interior stationary point where the segments are
/// not parallel, else on one of its four sides, each a point's distance to
/// the other segment.
fn segment_segment_distance((a, b): (Point, Point), (c, d): (Point, Point)) -> f64 {
    let to_segment = |p: Point, (u, v): (Point, Point)| {
        let w = v - u;
        let l2 = w.dot(w);
        let k = if l2 > 0.0 {
            ((p - u).dot(w) / l2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        p.distance(u + w * k)
    };
    let mut least = to_segment(a, (c, d))
        .min(to_segment(b, (c, d)))
        .min(to_segment(c, (a, b)))
        .min(to_segment(d, (a, b)));
    let (e, f, r) = (b - a, d - c, a - c);
    let (aa, ee, bb) = (e.dot(e), f.dot(f), e.dot(f));
    let den = aa * ee - bb * bb;
    if den > 0.0 {
        let (cc, ff) = (e.dot(r), f.dot(r));
        let s = (bb * ff - cc * ee) / den;
        let t = (aa * ff - bb * cc) / den;
        if (0.0..=1.0).contains(&s) && (0.0..=1.0).contains(&t) {
            least = least.min((a + e * s).distance(c + f * t));
        }
    }
    least
}

/// The four corners of the stretch of a plane between the lines through
/// `from` and `to` along the crease, each run from end plane to end
/// plane.
fn stretch(ends: &Ends, from: Point, to: Point) -> [Point; 4] {
    let t = ends.t;
    [
        from + t * ends.meets(0, from),
        from + t * ends.meets(1, from),
        to + t * ends.meets(1, to),
        to + t * ends.meets(0, to),
    ]
}

/// The least and greatest `s` at which the lines along the crease through
/// the section's corners, and through the points of its arc from the touch
/// (`ring[1]`) to the contact (`ring[2]`) about `centre`, meet the end
/// planes. On each plane `s` is linear in the point, so over the arc it
/// is a sum of a cosine and a sine of the arc's angle, its extremes read
/// in closed form.
fn run_span(ends: &Ends, ring: &[Point], centre: Point, radius: f64) -> (f64, f64) {
    if ends.square {
        return (0.0, ends.length);
    }
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for p in ring {
        lo = lo.min(ends.meets(0, *p));
        hi = hi.max(ends.meets(1, *p));
    }
    let (sweep, x, y) = arc_frame(ends.t, centre, radius, ring[1], ring[2]);
    for (i, (o, m)) in ends.planes.iter().enumerate() {
        // s(theta) = base + p cos theta + q sin theta.
        let scale = 1.0 / ends.t.dot(*m);
        let base = (*o - centre).dot(*m) * scale;
        let (p, q) = (-radius * x.dot(*m) * scale, -radius * y.dot(*m) * scale);
        if i == 0 {
            lo = lo.min(base - trig_max(-p, -q, 0.0, sweep));
        } else {
            hi = hi.max(base + trig_max(p, q, 0.0, sweep));
        }
    }
    (lo, hi)
}

/// The ball's arc from `from` to `to` about `centre`, the short way round
/// in the section: its sweep and the frame its angle is read in, `x`
/// toward `from` and `y` a quarter turn on toward `to`.
fn arc_frame(
    t: Vector,
    centre: Point,
    radius: f64,
    from: Point,
    to: Point,
) -> (f64, Vector, Vector) {
    let x = (from - centre) / radius;
    let y = t.cross(x);
    let w = to - centre;
    let angle = w.dot(y).atan2(w.dot(x));
    if angle >= 0.0 {
        (angle, x, y)
    } else {
        (-angle, x, -y)
    }
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
    let along = run.ends.t;
    // The section where the run starts: the plane square to the crease at
    // the least `s` its lines reach, so the prism covers both end planes.
    let lift = along * run.span.0;
    let mut ring: Vec<Point> = vec![run.crease, run.touch, run.contact];
    ring.extend(run.between.iter().copied());
    let ring: Vec<Point> = ring.into_iter().map(|p| p + lift).collect();
    let (crease, touch, contact, centre) = (ring[0], ring[1], ring[2], run.centre + lift);
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
    let x_dir = (touch - centre) / radius;
    let y_dir = along.cross(x_dir);
    let to = contact - centre;
    let angle = to.dot(y_dir).atan2(to.dot(x_dir));
    let normal = if angle >= 0.0 { along } else { -along };
    let frame = Frame::new(
        centre,
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
    let reach = ring.iter().fold(radius, |m, p| m.max(p.distance(crease))) * 4.0;
    let plane: SurfaceGeometry = PlaneSurface::over(
        Plane::through(crease, Direction::new(normal, tol)?),
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
    let travel = along * (run.span.1 - run.span.0);
    let wedge = ogeom_algo::make_prism(model, &section, travel, tol)?.shape;
    let mut wedge = explore(model, &wedge, Filter::OfType(ShapeType::Solid))?
        .into_iter()
        .next()
        .unwrap_or(wedge);
    // On slanted ends the prism runs past them; what lies beyond each end
    // plane goes.
    if !run.ends.square {
        let span = (run.span.1 - run.span.0) + reach;
        for (origin, outward) in run.ends.planes {
            wedge = crate::support::trimmed_wedge(model, &wedge, (origin, -outward), span, tol)?;
        }
    }
    crate::support::apply_wedge_solid(model, solid, Some(edge), &wedge, true, tol)
}
