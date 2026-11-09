//! A round of a straight edge that takes a narrow face off whole.
//!
//! Where a planar face beside a convex straight edge is a strip narrower
//! than the ball's setback (a chamfer, a bevelled rib), no ball touches
//! both of the edge's faces inside the material. The ball that rounds the
//! edge then touches the two faces either side of the strip, and where the
//! strip lies wholly in the corner that ball takes off, the round is the
//! blend of those two faces, seated on the line their planes meet on, and
//! the strip goes with the corner. With the strip's ends square to the
//! line the round is that seat's; with them slanted it is the section's
//! prism cut back to the end planes ([`consumed_slanted`]).

use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_geom::Curve;
use ogeom_geom::Curve3d as _;
use ogeom_math::{Aabb, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use crate::support::{Seat, edge_curve, planar_seat, same_occurrence};

/// The seat between the faces either side of a narrow strip beside `edge`,
/// where the strip is taken off whole by the ball of `radius` rolling
/// between them; `None` where it is not.
///
/// Every condition is read in closed form. The strip `F2` lies between the
/// edge `E` and an edge `E'` parallel to it, both convex, `E` on the face
/// `F1` and `E'` on the face `F3`, all three planar, with the strip's two
/// ends on planes square to the line `F1` and `F3` meet on. The solid by
/// the strip is then a prism, and in its section, with `V` that line, `c`
/// the ball's centre a radius in from both planes and `T1`, `T3` its
/// contacts:
///
/// - `E` stands on the leg from `V` to `T1`, short of `T1`, and `E'` on the
///   leg from `V` to `T3`, short of `T3`; the strip's section, the segment
///   from `E` to `E'`, lies in the triangle `V T1 T3`;
/// - `c` stands more than a radius behind the strip's plane, so the segment
///   stays off the ball, on the corner's side of its arc.
///
/// So the strip lies in the corner between the two legs and the arc, which
/// the cut takes off. What else of the solid that corner holds is bounded
/// by `F1` from `E` to `T1`, `F3` from `E'` to `T3` and the end faces over
/// the section from the strip to the contacts: no edge of those faces may
/// enter those regions, so each is wholly face, and the band stands on its
/// faces and its caps flush with the ends. Nothing else of the solid comes
/// near: no other face's bound meets the corner's.
///
/// # Errors
///
/// [`OgeomError::Dangling`] if a handle fails to resolve, and whatever a
/// face's bound reports; a solid this does not read is `None`.
pub(crate) fn consumed_seat(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Seat>> {
    let Some(seat) = seat_of(model, solid, edge, tol)? else {
        return Ok(None);
    };
    if !seat.convex {
        return Ok(None);
    }
    for narrow in 0..2 {
        if let Some(found) = across(model, solid, edge, &seat, narrow, radius, tol)? {
            return Ok(Some(found));
        }
    }
    Ok(None)
}

/// The straight seat of `edge`, `None` where it has none.
pub(crate) fn seat_of(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Seat>> {
    match planar_seat(model, solid, edge, tol) {
        Ok(seat) => Ok(Some(seat)),
        Err(OgeomError::Construction(_)) => Ok(None),
        Err(e) => Err(e),
    }
}

/// A straight edge's two ends, placed in the world; `None` for a curved one.
pub(crate) fn segment(
    model: &Model,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<[Point; 2]>> {
    if model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .is_none_or(|d| d.curve3d().is_none())
    {
        return Ok(None);
    }
    let (curve, range) = edge_curve(model, edge, tol)?;
    if !matches!(curve, Curve::Line(_)) {
        return Ok(None);
    }
    Ok(Some([
        curve.point_at(range.0, tol)?,
        curve.point_at(range.1, tol)?,
    ]))
}

/// Whether two edges share a vertex.
fn share_vertex(model: &Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<bool> {
    let (Some((a0, a1)), Some((b0, b1))) = (
        ogeom_algo::edge_vertices(model, a)?,
        ogeom_algo::edge_vertices(model, b)?,
    ) else {
        return Ok(false);
    };
    Ok([&a0, &a1]
        .iter()
        .any(|v| same_occurrence(model, v, &b0, tol) || same_occurrence(model, v, &b1, tol)))
}

/// The face of `solid` other than `face` that `edge` bounds, its plane's
/// outward normal and the way it leaves the edge; `None` where there is not
/// exactly one or it is curved.
fn beyond(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vector, Vector)>> {
    let Some(seat) = seat_of(model, solid, edge, tol)? else {
        return Ok(None);
    };
    let Some(k) = seat.faces.iter().position(|f| f.is_same(face)) else {
        return Ok(None);
    };
    Ok(Some((
        seat.faces[1 - k].clone(),
        seat.normals[1 - k],
        seat.leg(1 - k, tol)?,
    )))
}

/// The seat across the strip on face `narrow` of `seat`, where the ball
/// takes the strip off.
fn across(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    seat: &Seat,
    narrow: usize,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Seat>> {
    let conf = tol.confusion();
    let strip = &seat.faces[narrow];
    let (f1, n1) = (&seat.faces[1 - narrow], seat.normals[1 - narrow]);
    let n2 = seat.normals[narrow];
    // The strip: four straight edges, `E`, the far edge `E'` parallel to
    // it, and one across at each end.
    let bounds = explore_unique(model, strip, ShapeType::Edge)?;
    if bounds.len() != 4 || !bounds.iter().any(|b| same_occurrence(model, b, edge, tol)) {
        return Ok(None);
    }
    let mut far: Option<Shape> = None;
    let mut ends: Vec<Shape> = Vec::new();
    for b in &bounds {
        if same_occurrence(model, b, edge, tol) {
            continue;
        }
        if segment(model, b, tol)?.is_none() {
            return Ok(None);
        }
        if share_vertex(model, b, edge, tol)? {
            ends.push(b.clone());
        } else if far.replace(b.clone()).is_some() {
            return Ok(None);
        }
    }
    let (Some(far), [end0, end1]) = (far, ends.as_slice()) else {
        return Ok(None);
    };
    let Some(far_seat) = seat_of(model, solid, &far, tol)? else {
        return Ok(None);
    };
    if !far_seat.convex || far_seat.along.cross(seat.along).magnitude() > tol.angular() {
        return Ok(None);
    }
    let Some(k) = far_seat.faces.iter().position(|f| f.is_same(strip)) else {
        return Ok(None);
    };
    let (f3, n3) = (&far_seat.faces[1 - k], far_seat.normals[1 - k]);
    if f3.is_same(f1) {
        return Ok(None);
    }
    // The line the outer planes meet on runs with the edge, both lying in
    // planes holding the edge's direction.
    let meet = n1.cross(n3);
    if meet.magnitude() <= tol.angular() {
        return Ok(None);
    }
    let along = seat.along;
    let length = seat.start.distance(seat.end);
    // The section through the edge's start, and the far edge's foot there.
    let e0 = seat.start;
    let e0_far = far_seat.start + far_seat.along * (e0 - far_seat.start).dot(far_seat.along);
    // The centre: a radius behind both planes, `n1.o = n3.o = -r` with `o`
    // square to the line, so `o = -r (n1 + n3) / (1 + n1.n3)`; and the
    // line's point in the section, `v`, on both planes.
    let k13 = n1.dot(n3);
    let to_line = {
        // `v = e0 + a n1 + b n3` with `n1.(v - e0) = 0`, `n3.(v - e0') = 0`.
        let h = n3.dot(e0_far - e0);
        let det = 1.0 - k13 * k13;
        (n1 * (-k13 * h) + n3 * h) / det
    };
    let v = e0 + to_line;
    let centre = v - (n1 + n3) * (radius / (1.0 + k13));
    let t1 = centre + n1 * radius;
    let t3 = centre + n3 * radius;
    // `E` on the leg from `v` to `t1`, short of `t1`, and `E'` likewise.
    let on_leg = |p: Point, t: Point| -> bool {
        let leg = t - v;
        let reach = leg.magnitude();
        let at = (p - v).dot(leg / reach);
        at > conf && at < reach - conf
    };
    if !on_leg(e0, t1) || !on_leg(e0_far, t3) {
        return Ok(None);
    }
    // The strip off the ball, on the corner's side of the arc.
    if n2.dot(centre - e0) >= -(radius + conf) {
        return Ok(None);
    }
    // Each end across on a plane square to the edge, the strip's two end
    // edges each at one end of `E` and of `E'` alike.
    let s_of = |p: Point| (p - e0).dot(along);
    let far_ends = [far_seat.start, far_seat.end].map(s_of);
    let (lo, hi) = (far_ends[0].min(far_ends[1]), far_ends[0].max(far_ends[1]));
    if lo.abs() > conf || (hi - length).abs() > conf {
        return Ok(None);
    }
    let mut caps: Vec<(Shape, f64)> = Vec::new();
    for end in [end0, end1] {
        let Some((cap, normal, leaves)) = beyond(model, solid, end, strip, tol)? else {
            return Ok(None);
        };
        if normal.cross(along).magnitude() > tol.angular() {
            return Ok(None);
        }
        let Some(points) = segment(model, end, tol)? else {
            return Ok(None);
        };
        // The end face lies on the ball's side of the strip.
        if (centre - e0).dot(leaves) <= 0.0 {
            return Ok(None);
        }
        caps.push((cap, s_of(points[0])));
    }
    // The outer faces hold the band's contacts and stop at their edges: no
    // edge of `F1` other than `E` enters the rectangle from the line to the
    // contact, nor of `F3` other than `E'` its own. Each face lies on the
    // contact's side of its edge, so it holds the part from the edge to the
    // contact whole and none of the part from the line to the edge.
    let rectangle = |foot: Point, contact: Point| {
        [
            foot,
            foot + along * length,
            contact + along * length,
            contact,
        ]
    };
    if (t1 - e0).dot(seat.leg(1 - narrow, tol)?) <= 0.0
        || (t3 - e0_far).dot(far_seat.leg(1 - k, tol)?) <= 0.0
    {
        return Ok(None);
    }
    if !clear(model, f1, &rectangle(v, t1), &[edge], tol)?
        || !clear(model, f3, &rectangle(v, t3), &[&far], tol)?
    {
        return Ok(None);
    }
    // The ends hold the band's caps likewise: no edge of an end face other
    // than the strip's end edge enters the section there, so the face holds
    // all of it from the strip to the contacts and none of it before.
    for ((cap, s), end) in caps.iter().zip([end0, end1]) {
        let shift = along * *s;
        let section = [v, t1, t3].map(|p| p + shift);
        if !clear(model, cap, &section, &[end], tol)? {
            return Ok(None);
        }
    }
    // No other face of the solid comes near the corner the cut takes: each
    // one's bound stays clear of the corner's.
    let corner = {
        let section = [v, t1, t3];
        let mut points = section.to_vec();
        points.extend(section.map(|p| p + along * length));
        Aabb::of_points(&points)
    };
    let known = [f1, strip, f3, &caps[0].0, &caps[1].0];
    for face in explore_unique(model, solid, ShapeType::Face)? {
        if known.iter().any(|k| k.is_same(&face)) {
            continue;
        }
        let overlap = ogeom_algo::shape_bounds(model, &face, tol)?.intersection(&corner);
        let size = overlap.size();
        if overlap.low().is_some() && size.x > conf && size.y > conf && size.z > conf {
            return Ok(None);
        }
    }
    let start = v;
    Ok(Some(Seat {
        start,
        end: start + along * length,
        along,
        normals: [n1, n3],
        faces: [f1.clone(), f3.clone()],
        convex: true,
    }))
}

/// Whether no edge of `face` but those in `skip` enters the open convex
/// polygon `corners`, which lies in the face's plane: every edge straight,
/// each clipped to the polygon shrunk by the confusion distance.
fn clear(
    model: &Model,
    face: &Shape,
    corners: &[Point],
    skip: &[&Shape],
    tol: Tolerances,
) -> OgeomResult<bool> {
    let conf = tol.confusion();
    let n = corners.len();
    // The polygon turns positively about this normal, so each side's
    // inward normal within the plane is the normal crossed with the side.
    let plane_normal = (corners[1] - corners[0]).cross(corners[2] - corners[0]);
    let mut sides: Vec<(Point, Vector)> = Vec::with_capacity(n);
    for i in 0..n {
        let (a, b) = (corners[i], corners[(i + 1) % n]);
        let inward = plane_normal.cross(b - a);
        let m = inward.magnitude();
        if m <= 0.0 {
            return Ok(false);
        }
        sides.push((a, inward / m));
    }
    for e in explore_unique(model, face, ShapeType::Edge)? {
        if skip.iter().any(|s| same_occurrence(model, &e, s, tol)) {
            continue;
        }
        let Some([p, q]) = segment(model, &e, tol)? else {
            return Ok(false);
        };
        // The part of `p + t (q - p)`, `t` in [0, 1], strictly inside.
        let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
        for (a, inward) in &sides {
            let start = (p - *a).dot(*inward) - conf;
            let rate = (q - p).dot(*inward);
            if rate == 0.0 {
                if start <= 0.0 {
                    hi = lo;
                }
                continue;
            }
            let t = -start / rate;
            if rate > 0.0 {
                lo = lo.max(t);
            } else {
                hi = hi.min(t);
            }
        }
        if lo < hi {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The round of the convex straight `edge` that takes the narrow strip
/// beside it off whole, as [`consumed_seat`] finds it, where the strip's
/// ends lie on planes slanted to the edge rather than square to it; `None`
/// where it is not of that kind (square ends are [`consumed_seat`]'s).
///
/// The section through the edge's start is as [`consumed_seat`] reads it:
/// the line `V` the outer planes meet on, the ball's centre a radius in
/// from both, its contacts `T1` and `T3`, `E` and `E'` on the legs short of
/// them, the strip off the ball. Each end of the strip lies on the plane of
/// the face across its end edge, and every line along the edge of that
/// section (`V`, `T1`, `T3`, `E`, `E'`) runs between the two planes; the
/// corner the cut takes off is the section's prism between them. Each
/// check of the square case is asked of that prism: `F1` and `F3` hold the
/// stretches from the line to the contacts between the planes, each end
/// face holds the corner's trace on its plane, and no other face's bound
/// meets the corner's. The cut is the section's prism, run past both
/// planes and cut back to them, taken off the solid.
///
/// # Errors
///
/// [`OgeomError::Dangling`] if a handle fails to resolve, and whatever a
/// face's bound or the boolean reports.
#[allow(clippy::too_many_lines, reason = "one corner, read and cut end to end")]
pub(crate) fn consumed_slanted(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<ogeom_algo::Built>> {
    let Some(seat) = seat_of(model, solid, edge, tol)? else {
        return Ok(None);
    };
    if !seat.convex {
        return Ok(None);
    }
    for narrow in 0..2 {
        if let Some(cut) = slanted_corner(model, solid, edge, &seat, narrow, radius, tol)? {
            return cut_corner(model, solid, edge, &cut, radius, tol).map(Some);
        }
    }
    Ok(None)
}

/// The corner a slanted-ended strip's round takes off: its section through
/// the edge's start, the edge's direction, and each end's plane as a point
/// on it and its outward normal.
struct SlantedCorner {
    v: Point,
    t1: Point,
    t3: Point,
    centre: Point,
    along: Vector,
    planes: [(Point, Vector); 2],
}

impl SlantedCorner {
    /// Where the line through `p` along the edge meets end plane `i`.
    fn meets(&self, i: usize, p: Point) -> Point {
        let (o, m) = self.planes[i];
        p + self.along * ((o - p).dot(m) / self.along.dot(m))
    }
}

/// The corner across the strip on face `narrow` of `seat`, where its ends
/// are slanted and the ball takes it off.
#[allow(
    clippy::too_many_lines,
    reason = "one corner, read and proved end to end"
)]
fn slanted_corner(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    seat: &Seat,
    narrow: usize,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<SlantedCorner>> {
    let conf = tol.confusion();
    let strip = &seat.faces[narrow];
    let (f1, n1) = (&seat.faces[1 - narrow], seat.normals[1 - narrow]);
    let n2 = seat.normals[narrow];
    let bounds = explore_unique(model, strip, ShapeType::Edge)?;
    if bounds.len() != 4 || !bounds.iter().any(|b| same_occurrence(model, b, edge, tol)) {
        return Ok(None);
    }
    let mut far: Option<Shape> = None;
    let mut ends: Vec<Shape> = Vec::new();
    for b in &bounds {
        if same_occurrence(model, b, edge, tol) {
            continue;
        }
        if segment(model, b, tol)?.is_none() {
            return Ok(None);
        }
        if share_vertex(model, b, edge, tol)? {
            ends.push(b.clone());
        } else if far.replace(b.clone()).is_some() {
            return Ok(None);
        }
    }
    let (Some(far), [end0, end1]) = (far, ends.as_slice()) else {
        return Ok(None);
    };
    let Some(far_seat) = seat_of(model, solid, &far, tol)? else {
        return Ok(None);
    };
    if !far_seat.convex || far_seat.along.cross(seat.along).magnitude() > tol.angular() {
        return Ok(None);
    }
    let Some(k) = far_seat.faces.iter().position(|f| f.is_same(strip)) else {
        return Ok(None);
    };
    let (f3, n3) = (&far_seat.faces[1 - k], far_seat.normals[1 - k]);
    if f3.is_same(f1) || n1.cross(n3).magnitude() <= tol.angular() {
        return Ok(None);
    }
    let along = seat.along;
    let e0 = seat.start;
    let e0_far = far_seat.start + far_seat.along * (e0 - far_seat.start).dot(far_seat.along);
    let k13 = n1.dot(n3);
    let to_line = {
        let h = n3.dot(e0_far - e0);
        let det = 1.0 - k13 * k13;
        (n1 * (-k13 * h) + n3 * h) / det
    };
    let v = e0 + to_line;
    let centre = v - (n1 + n3) * (radius / (1.0 + k13));
    let t1 = centre + n1 * radius;
    let t3 = centre + n3 * radius;
    let on_leg = |p: Point, t: Point| -> bool {
        let leg = t - v;
        let reach = leg.magnitude();
        let at = (p - v).dot(leg / reach);
        at > conf && at < reach - conf
    };
    if !on_leg(e0, t1) || !on_leg(e0_far, t3) {
        return Ok(None);
    }
    if n2.dot(centre - e0) >= -(radius + conf) {
        return Ok(None);
    }
    // Each end across on the plane of the face beyond its end edge, the
    // start's facing back along the edge and the end's on along it, at
    // least one slanted to it (square ends are the seat's own).
    let mut planes = [(e0, -along), (seat.end, along)];
    let mut caps: Vec<(Shape, usize)> = Vec::new();
    for end in [end0, end1] {
        let Some((cap, normal, leaves)) = beyond(model, solid, end, strip, tol)? else {
            return Ok(None);
        };
        let Some(points) = segment(model, end, tol)? else {
            return Ok(None);
        };
        if (centre - e0).dot(leaves) <= 0.0 {
            return Ok(None);
        }
        let i = usize::from((points[0] - e0).dot(along) > (seat.end - e0).dot(along) * 0.5);
        let facing = if i == 0 {
            -normal.dot(along)
        } else {
            normal.dot(along)
        };
        if facing <= tol.angular() {
            return Ok(None);
        }
        planes[i] = (points[0], normal);
        caps.push((cap, i));
    }
    if caps[0].1 == caps[1].1 {
        return Ok(None);
    }
    if planes
        .iter()
        .all(|(_, m)| m.cross(along).magnitude() <= tol.angular())
    {
        return Ok(None);
    }
    let corner = SlantedCorner {
        v,
        t1,
        t3,
        centre,
        along,
        planes,
    };
    // `E` and `E'` run from plane to plane.
    let on = |i: usize, p: Point| {
        let (o, m) = corner.planes[i];
        (p - o).dot(m).abs() <= conf
    };
    for line in [[seat.start, seat.end], [far_seat.start, far_seat.end]] {
        if !((on(0, line[0]) && on(1, line[1])) || (on(0, line[1]) && on(1, line[0]))) {
            return Ok(None);
        }
    }
    if (t1 - e0).dot(seat.leg(1 - narrow, tol)?) <= 0.0
        || (t3 - e0_far).dot(far_seat.leg(1 - k, tol)?) <= 0.0
    {
        return Ok(None);
    }
    // The outer faces hold the stretches from the line to their contacts
    // between the planes, and each end face the corner's trace on its plane.
    let stretch = |foot: Point, contact: Point| {
        [
            corner.meets(0, foot),
            corner.meets(1, foot),
            corner.meets(1, contact),
            corner.meets(0, contact),
        ]
    };
    if !clear(model, f1, &stretch(v, t1), &[edge], tol)?
        || !clear(model, f3, &stretch(v, t3), &[&far], tol)?
    {
        return Ok(None);
    }
    for ((cap, i), end) in caps.iter().zip([end0, end1]) {
        let trace = [v, t1, t3].map(|p| corner.meets(*i, p));
        if !clear(model, cap, &trace, &[end], tol)? {
            return Ok(None);
        }
    }
    // No other face of the solid comes near the corner the cut takes.
    let mut points = Vec::with_capacity(6);
    for i in 0..2 {
        points.extend([v, t1, t3].map(|p| corner.meets(i, p)));
    }
    let bound = Aabb::of_points(&points);
    let known = [f1, strip, f3, &caps[0].0, &caps[1].0];
    for face in explore_unique(model, solid, ShapeType::Face)? {
        if known.iter().any(|k| k.is_same(&face)) {
            continue;
        }
        let overlap = ogeom_algo::shape_bounds(model, &face, tol)?.intersection(&bound);
        let size = overlap.size();
        if overlap.low().is_some() && size.x > conf && size.y > conf && size.z > conf {
            return Ok(None);
        }
    }
    Ok(Some(corner))
}

/// Take the corner off: the prism of its section, from the line to each
/// contact and the ball's arc between them, run from the section's least
/// reach along the edge to its greatest and cut back to the end planes.
fn cut_corner(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    corner: &SlantedCorner,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<ogeom_algo::Built> {
    use ogeom_geom::{CircleCurve, PlaneSurface, SurfaceGeometry};
    use ogeom_math::{Circle, Direction, Frame, Plane};
    use ogeom_topo::{Filter, explore};
    let along = corner.along;
    // The arc lies in the triangle from the line to the contacts, so the
    // triangle's corners bound the reach.
    let reach = |i: usize, p: Point| (corner.meets(i, p) - p).dot(along);
    let section = [corner.v, corner.t1, corner.t3];
    let lo = section
        .iter()
        .map(|p| reach(0, *p))
        .fold(f64::INFINITY, f64::min);
    let hi = section
        .iter()
        .map(|p| reach(1, *p))
        .fold(f64::NEG_INFINITY, f64::max);
    let lift = along * lo;
    let [v, t1, t3] = section.map(|p| p + lift);
    let centre = corner.centre + lift;
    let vertices = [v, t1, t3].map(|p| ogeom_algo::make_vertex(model, p).shape);
    let leg1 = crate::support::segment_between(model, (&vertices[0], v), (&vertices[1], t1), tol)?;
    let leg3 = crate::support::segment_between(model, (&vertices[2], t3), (&vertices[0], v), tol)?;
    // The ball's arc from `T1` to `T3` the short way round, toward the line.
    let x = (t1 - centre) / radius;
    let y = along.cross(x);
    let to = t3 - centre;
    let angle = to.dot(y).atan2(to.dot(x));
    let normal = if angle >= 0.0 { along } else { -along };
    let frame = Frame::new(
        centre,
        Direction::new(normal, tol)?,
        Direction::new(x, tol)?,
        tol,
    )?;
    let circle = CircleCurve::new(Circle::new(frame, radius, tol)?).into();
    let arc = ogeom_algo::make_edge_between(
        model,
        circle,
        (0.0, angle.abs()),
        &vertices[1],
        &vertices[2],
        tol,
    )?
    .shape;
    let turn = (t1 - v).cross(t3 - v).dot(along);
    let face_normal = if turn >= 0.0 { along } else { -along };
    let size = (v.distance(t1) + v.distance(t3) + radius) * 4.0;
    let plane: SurfaceGeometry = PlaneSurface::over(
        Plane::through(v, Direction::new(face_normal, tol)?),
        (-size, size),
        (-size, size),
    )?
    .into();
    let ring = if turn >= 0.0 {
        vec![leg1, arc, leg3]
    } else {
        vec![leg3.reversed(), arc.reversed(), leg1.reversed()]
    };
    let face = ogeom_algo::make_face_with_pcurves(model, plane, &[ring], tol)?.shape;
    let prism = ogeom_algo::make_prism(model, &face, along * (hi - lo), tol)?.shape;
    let mut wedge = explore(model, &prism, Filter::OfType(ShapeType::Solid))?
        .into_iter()
        .next()
        .unwrap_or(prism);
    let span = (hi - lo) + size;
    for (origin, outward) in corner.planes {
        wedge = crate::support::trimmed_wedge(model, &wedge, (origin, -outward), span, tol)?;
    }
    crate::support::apply_wedge_solid(model, solid, Some(edge), &wedge, false, tol)
}
