//! The old slanted cut of a round that takes a narrow face off whole: the
//! section's prism cut back to the strip's slanted end planes, on the seat
//! [`crate::engine::consumed_seat`] reads.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_math::{Aabb, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use crate::engine::consumed_seat::{beyond, clear, seat_of, segment, share_vertex};
use crate::support::{Seat, same_occurrence};

/// The round of the convex straight `edge` that takes the narrow strip
/// beside it off whole, as [`consumed_seat`](crate::engine::consumed_seat::consumed_seat) finds it, where the strip's
/// ends lie on planes slanted to the edge rather than square to it; `None`
/// where it is not of that kind (square ends are [`consumed_seat`](crate::engine::consumed_seat::consumed_seat)'s).
///
/// The section through the edge's start is as [`consumed_seat`](crate::engine::consumed_seat::consumed_seat) reads it:
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
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle fails to resolve, and whatever a
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
        if slanted_corner(model, solid, edge, &seat, narrow, radius, tol)?.is_some() {
            return Err(crate::fillet::declined());
        }
    }
    Ok(None)
}

/// The corner a slanted-ended strip's round takes off: the edge's
/// direction, and each end's plane as a point on it and its outward normal.
struct SlantedCorner {
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
    let corner = SlantedCorner { along, planes };
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
