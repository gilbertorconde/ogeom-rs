//! A concave edge's round whose ball rolls on the far edge of the narrow
//! face beside it.
//!
//! Where two solids are fused at a slant, the seam between them is a
//! concave edge `E`, and the face of one of them beside it is often a
//! strip `F2` running between the seam and that solid's own convex edge
//! `E'`. Where the strip is narrower than the ball's setback, the ball
//! touching the seam's other face `F1` cannot touch the strip: it rests on
//! `E'` instead, on a line rather than a face. Its centre then stands a
//! radius off `F1`'s plane and a radius off `E'`'s line, on the plane
//! offset from `F1` and the cylinder about `E'`, which meet in an ellipse,
//! or in a line where `E'` runs parallel to `F1`. The ball touches `F1` at
//! its centre's foot on that plane and `E'` at its centre's foot on the
//! line, so the round runs from `F1` to `E'` itself, where it meets the
//! face `F3` beyond `E'` at a crease.
//!
//! Every station is closed-form. The ball's arc stands in the plane
//! square to the centre's path, which holds both contacts, and crosses the
//! crease where that plane does. The round is fitted through those arcs
//! with its rails set exactly on `F1`'s plane and on `E'`'s line, and held
//! to the band's fit target against the envelope of the balls; the fill
//! is the marched band's wedge, capped in the end arcs' planes.

use crate::consumed::{seat_of, segment};
use crate::march::{BlendStop, MarchedBlend, Sides};
use crate::marched::{ExactBand, Rail};
use crate::support::{Seat, edge_curve, same_occurrence};
use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_geom::Curve3d as _;
use ogeom_geom::SurfaceGeometry;
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, NodeData, Orientation, Shape, ShapeType, explore_unique};

/// How many equal steps of the path's parameter the run is solved at: the
/// stations the band is fitted through.
const STATIONS: usize = 64;

/// How many samples the path's parameter is searched at for the run's ends
/// over the arc of the ellipse on which the ball rests on the edge.
const SEARCH: usize = 4096;

/// How many times an interval of the run is halved, at most, to show that
/// the ball stays clear of a face between two stations.
const HALVINGS: u32 = 24;

/// How many chords, at most, the line of contact is read with: a run
/// needing more to hold it within the confusion is refused.
const MAX_CHORDS: f64 = 1e5;

/// The path of the ball's centre, a radius off the first face's plane and
/// a radius off the far edge's line.
#[derive(Debug, Clone, Copy)]
enum Path {
    /// The edge oblique to the plane: the centre at angle `s` round the
    /// edge is `base + along l(s) + radius (cos s u + sin s v)`, its foot
    /// on the edge `base + along l(s)`, with
    /// `l(s) = l0 - radius (cos s au + sin s av) / k`.
    Ellipse {
        base: Point,
        along: Vector,
        u: Vector,
        v: Vector,
        radius: f64,
        l0: f64,
        au: f64,
        av: f64,
        k: f64,
    },
    /// The edge parallel to the plane: the centre is `base + along s + w`.
    Line {
        base: Point,
        along: Vector,
        w: Vector,
    },
}

impl Path {
    /// Where along the far edge's line, from its base, the foot stands.
    fn foot_at(&self, s: f64) -> f64 {
        match *self {
            Self::Ellipse {
                radius,
                l0,
                au,
                av,
                k,
                ..
            } => l0 - radius * (s.cos() * au + s.sin() * av) / k,
            Self::Line { .. } => s,
        }
    }

    /// The ball's contact on the far edge.
    fn foot(&self, s: f64) -> Point {
        match *self {
            Self::Ellipse { base, along, .. } | Self::Line { base, along, .. } => {
                base + along * self.foot_at(s)
            }
        }
    }

    /// From the contact on the far edge to the centre.
    fn reach(&self, s: f64) -> Vector {
        match *self {
            Self::Ellipse { u, v, radius, .. } => (u * s.cos() + v * s.sin()) * radius,
            Self::Line { w, .. } => w,
        }
    }

    fn centre(&self, s: f64) -> Point {
        self.foot(s) + self.reach(s)
    }

    fn d1(&self, s: f64) -> Vector {
        match *self {
            Self::Ellipse {
                along,
                u,
                v,
                radius,
                au,
                av,
                k,
                ..
            } => {
                along * (radius * (s.sin() * au - s.cos() * av) / k)
                    + (v * s.cos() - u * s.sin()) * radius
            }
            Self::Line { along, .. } => along,
        }
    }

    fn d2(&self, s: f64) -> Vector {
        match *self {
            Self::Ellipse {
                along,
                u,
                v,
                radius,
                au,
                av,
                k,
                ..
            } => {
                along * (radius * (s.cos() * au + s.sin() * av) / k)
                    - (u * s.cos() + v * s.sin()) * radius
            }
            Self::Line { .. } => Vector::ZERO,
        }
    }

    /// A bound on the centre's speed and on its acceleration in the path's
    /// parameter, over the whole path.
    fn speed(&self) -> f64 {
        match *self {
            Self::Ellipse {
                radius, au, av, k, ..
            } => radius * (1.0 + au.hypot(av) / k.abs()),
            Self::Line { .. } => 1.0,
        }
    }

    /// The least and greatest of [`Path::foot_at`] over `[lo, hi]`: at the
    /// ends, or where the foot turns back, square to the edge's line.
    fn foot_span(&self, lo: f64, hi: f64) -> (f64, f64) {
        let (lo, hi) = (lo.min(hi), lo.max(hi));
        let mut values = vec![self.foot_at(lo), self.foot_at(hi)];
        if let Self::Ellipse { au, av, .. } = *self {
            let turn = av.atan2(au);
            let period = core::f64::consts::PI;
            let mut s = turn - period * ((turn - lo) / period).ceil();
            while s <= hi {
                if s >= lo {
                    values.push(self.foot_at(s));
                }
                s += period;
            }
        }
        values
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &x| {
                (a.min(x), b.max(x))
            })
    }
}

/// The run of the ball rolling on the far edge, from the edge's start to
/// its end, with what the band is built from.
pub(crate) struct Rolled {
    /// The ball's stations, ordered along the crease.
    blend: MarchedBlend,
    /// The crease's curve, placed.
    guide: ogeom_geom::Curve,
    /// The touched face's plane and the narrow face's, placed, each with
    /// its outward sign.
    hosts: [(SurfaceGeometry, f64); 2],
    /// The touched face's plane and the far edge's line.
    rails: [Rail; 2],
    /// The centre's path and the run's span of its parameter.
    path: Path,
    span: (f64, f64),
}

/// A planar face's plane as a point and its outward unit normal, and its
/// boundary as straight segments; `None` where it is curved or bounded by
/// anything else.
struct Flat {
    origin: Point,
    normal: Vector,
    sides: Vec<[Point; 2]>,
}

fn flat(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Option<Flat>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    let Some(SurfaceGeometry::Plane(plane)) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let placement = face.transform(model.datums())?;
    let mut normal = placement.apply_vector(plane.plane().normal().vector());
    if face.orientation() == Orientation::Reversed {
        normal = -normal;
    }
    let origin = placement.apply(plane.plane().origin());
    let mut sides = Vec::new();
    for e in explore_unique(model, face, ShapeType::Edge)? {
        let Some(points) = segment(model, &e, tol)? else {
            return Ok(None);
        };
        sides.push(points);
    }
    Ok(Some(Flat {
        origin,
        normal,
        sides,
    }))
}

/// A unit basis of the plane square to `normal`.
fn basis(normal: Vector) -> (Vector, Vector) {
    let helper = if normal.x.abs() < 0.9 {
        Vector::new(1.0, 0.0, 0.0)
    } else {
        Vector::new(0.0, 1.0, 0.0)
    };
    let e1 = normal.cross(helper);
    let e1 = e1 / e1.magnitude();
    (e1, normal.cross(e1))
}

/// Whether `p`, in the plane of `sides`, lies inside the region they
/// bound, read by the parity of the sides a ray along `e1` crosses.
fn inside(sides: &[[Point; 2]], p: Point, (e1, e2): (Vector, Vector)) -> bool {
    let at = |q: Point| ((q - p).dot(e1), (q - p).dot(e2));
    let mut odd = false;
    for [a, b] in sides {
        let ((ax, ay), (bx, by)) = (at(*a), at(*b));
        if (ay > 0.0) != (by > 0.0) {
            let x = ax + (bx - ax) * (0.0 - ay) / (by - ay);
            if x > 0.0 {
                odd = !odd;
            }
        }
    }
    odd
}

fn segment_distance(p: Point, [a, b]: [Point; 2]) -> f64 {
    let ab = b - a;
    let l2 = ab.dot(ab);
    let t = if l2 > 0.0 {
        ((p - a).dot(ab) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    p.distance(a + ab * t)
}

/// The distance from `p` to a planar face.
fn face_distance(face: &Flat, p: Point) -> f64 {
    let height = (p - face.origin).dot(face.normal);
    let foot = p - face.normal * height;
    if inside(&face.sides, foot, basis(face.normal)) {
        return height.abs();
    }
    face.sides
        .iter()
        .map(|s| segment_distance(p, *s))
        .fold(f64::INFINITY, f64::min)
}

/// Whether some stretch of the segment `side` lies inside the simple
/// polygon `ring` in the plane square to `normal`, more than `margin` from
/// its boundary.
fn enters(ring: &[Point], normal: Vector, [p, q]: [Point; 2], margin: f64) -> bool {
    let frame = basis(normal);
    let sides: Vec<[Point; 2]> = (0..ring.len())
        .map(|i| [ring[i], ring[(i + 1) % ring.len()]])
        .collect();
    let (e1, e2) = frame;
    let flat = |x: Point| ((x - p).dot(e1), (x - p).dot(e2));
    let (qx, qy) = flat(q);
    let mut cuts = vec![0.0, 1.0];
    for [a, b] in &sides {
        let ((ax, ay), (bx, by)) = (flat(*a), flat(*b));
        let (dx, dy) = (bx - ax, by - ay);
        let den = qx * dy - qy * dx;
        if den == 0.0 {
            continue;
        }
        let t = (ax * dy - ay * dx) / den;
        let u = (ax * qy - ay * qx) / den;
        if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
            cuts.push(t);
        }
    }
    cuts.sort_by(f64::total_cmp);
    for pair in cuts.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        if b - a <= 0.0 {
            continue;
        }
        for k in 1..8 {
            let x = p + (q - p) * (a + (b - a) * f64::from(k) / 8.0);
            let clear = sides.iter().all(|s| segment_distance(x, *s) > margin);
            if clear && inside(&sides, x, frame) {
                return true;
            }
        }
    }
    false
}

/// The run of a ball of `radius` rounding the concave straight `edge`
/// while it rests on the far edge of the narrow planar face beside it;
/// `None` where the edge is not of that kind.
///
/// # Errors
///
/// [`OgeomError::Construction`] naming why, where the ball does rest on
/// such a far edge but no sound round of this kind exists: the ball comes
/// down on a face either side of that edge or runs past its end along the
/// crease, touches its own face off that face, comes within its radius of
/// another face, or its sections cross; and whatever reading the solid
/// reports.
pub(crate) fn rolled_run(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Rolled>> {
    let Some(seat) = seat_of(model, solid, edge, tol)? else {
        return Ok(None);
    };
    if seat.convex {
        return Ok(None);
    }
    let mut named: Option<OgeomError> = None;
    for narrow in 0..2 {
        match across(model, solid, edge, &seat, narrow, radius, tol) {
            Ok(Some(run)) => return Ok(Some(run)),
            Ok(None) => {}
            Err(e @ OgeomError::Construction(_)) => {
                named.get_or_insert(e);
            }
            Err(e) => return Err(e),
        }
    }
    named.map_or(Ok(None), Err)
}

/// The far edge the ball meets crossing face `narrow` from the middle of
/// the crease, square to it, short of the setback; `None` where it meets
/// none, or the first edge it meets is not straight.
fn far_edge(
    model: &Model,
    edge: &Shape,
    seat: &Seat,
    narrow: usize,
    setback: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let face = &seat.faces[narrow];
    let into = -seat.leg(narrow, tol)?;
    let mid = seat.start.midpoint(seat.end);
    let normal = seat.normals[narrow];
    let mut best: Option<(f64, Shape)> = None;
    for g in explore_unique(model, face, ShapeType::Edge)? {
        if same_occurrence(model, &g, edge, tol) {
            continue;
        }
        let Some([a, b]) = segment(model, &g, tol)? else {
            return Ok(None);
        };
        let w = b - a;
        let den = into.cross(w).dot(normal);
        if den.abs() <= tol.angular() * w.magnitude() {
            continue;
        }
        let s = (a - mid).cross(w).dot(normal) / den;
        let u = (a - mid).cross(into).dot(normal) / den;
        if s > tol.confusion() && (0.0..=1.0).contains(&u) && best.as_ref().is_none_or(|b| s < b.0)
        {
            best = Some((s, g));
        }
    }
    Ok(best.and_then(|(s, g)| (s < setback).then_some(g)))
}

/// A face's placed surface, with its outward sign.
fn host(
    model: &Model,
    face: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<(SurfaceGeometry, f64)>> {
    use ogeom_geom::Transformable as _;
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let placed = surface
        .clone()
        .transformed(&face.transform(model.datums())?, tol)?;
    let sign = if face.orientation() == Orientation::Reversed {
        -1.0
    } else {
        1.0
    };
    Ok(Some((placed, sign)))
}

/// The run across face `narrow` of `seat`, where its far edge holds the
/// ball.
#[allow(clippy::too_many_lines, reason = "one run, read and proved end to end")]
fn across(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    seat: &Seat,
    narrow: usize,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Rolled>> {
    let conf = tol.confusion();
    let (n1, n2) = (seat.normals[1 - narrow], seat.normals[narrow]);
    let (f1, f2) = (&seat.faces[1 - narrow], &seat.faces[narrow]);
    let turn = n1.dot(n2).clamp(-1.0, 1.0).acos();
    let setback = radius * (turn * 0.5).tan();
    let Some(far) = far_edge(model, edge, seat, narrow, setback, tol)? else {
        return Ok(None);
    };
    let Some(far_seat) = seat_of(model, solid, &far, tol)? else {
        return Ok(None);
    };
    let Some(k3) = far_seat.faces.iter().position(|f| !f.is_same(f2)) else {
        return Ok(None);
    };
    if !far_seat.convex || far_seat.faces.iter().all(|f| !f.is_same(f2)) {
        return Ok(None);
    }
    let f3 = &far_seat.faces[k3];
    if f3.is_same(f1) {
        return Ok(None);
    }
    let n3 = far_seat.normals[k3];
    let (curve, (lo, hi)) = edge_curve(model, edge, tol)?;
    let e0 = curve.point_at(lo, tol)?;
    let e1 = curve.point_at(hi, tol)?;
    let length = e0.distance(e1);
    let t = (e1 - e0) / length;
    let mid = e0.midpoint(e1);
    let refuse = |why: &str| -> OgeomError {
        ogeom_core::ogeom_err!(
            Construction,
            "a blend of radius {radius} on the edge through {mid:?} rests on the far edge of \
             the narrow face beside it, and {why}; see docs/CAPABILITIES.md, fillet.edge-blends"
        )
    };

    // The far edge's line, and the frame round it: `u` the narrow face's
    // normal, `v` square to it toward the beyond face's, which stands at
    // `phi3` round from `u`. The ball rests on the edge alone where its
    // centre stands between the two normals.
    let base = far_seat.start;
    let along = far_seat.along;
    let reach = far_seat.start.distance(far_seat.end);
    let u = n2;
    let v = {
        let w = n3 - u * n3.dot(u);
        w / w.magnitude()
    };
    let phi3 = n3.dot(v).atan2(n3.dot(u));
    let k = n1.dot(along);
    let height = radius - n1.dot(base - e0);
    let (au, av) = (n1.dot(u), n1.dot(v));
    // The edge counts as parallel to the face where its height over the
    // face changes by less than the confusion along the crease.
    let path = if k.abs() * length > conf {
        Path::Ellipse {
            base,
            along,
            u,
            v,
            radius,
            l0: height / k,
            au,
            av,
            k,
        }
    } else {
        // Parallel: `radius (cos a au + sin a av) = height` at one angle
        // between the two normals, the height read at the crease's start.
        let level = height - k * (e0 - base).dot(along);
        let m = au.hypot(av);
        let c = level / (radius * m);
        if !(-1.0..=1.0).contains(&c) {
            return Err(refuse(
                "no ball rests on that edge and touches its own face",
            ));
        }
        let mid_angle = av.atan2(au);
        let angles = [mid_angle - c.acos(), mid_angle + c.acos()];
        let held: Vec<f64> = angles
            .iter()
            .map(|a| a.rem_euclid(core::f64::consts::TAU))
            .filter(|a| (0.0..=phi3).contains(a))
            .collect();
        let [angle] = held.as_slice() else {
            return Err(refuse("no one ball rests on that edge alone"));
        };
        Path::Line {
            base,
            along,
            w: (u * angle.cos() + v * angle.sin()) * radius,
        }
    };
    // Where the arc's plane, square to the path, crosses the crease.
    let crossing = |s: f64| -> Option<f64> {
        let d = path.d1(s);
        let den = d.dot(t);
        (den != 0.0).then(|| d.dot(path.centre(s) - e0) / den)
    };
    let span = match path {
        Path::Line { .. } => {
            // The crease runs with the path: the run is the crease's span
            // along it.
            let at = |x: f64| (e0 + t * x - path.centre(0.0)).dot(along);
            (at(0.0), at(length))
        }
        Path::Ellipse { .. } => {
            let Some(span) = ellipse_run(&crossing, phi3, length) else {
                let at_start = crossing(0.0).is_some_and(|x| x > 0.0 && x < length);
                let at_end = crossing(phi3).is_some_and(|x| x > 0.0 && x < length);
                return Err(refuse(if at_start {
                    "the ball comes down on that face along the crease"
                } else if at_end {
                    "the ball comes onto the face beyond that edge along the crease"
                } else {
                    "the ball does not rest on that edge over the whole crease"
                }));
            };
            span
        }
    };
    let (sa, sb) = span;
    #[allow(clippy::cast_precision_loss)]
    let station = |i: usize| sa + (sb - sa) * (i as f64) / (STATIONS as f64);

    // The contact on the far edge stays on it, and moves one way along it.
    let (foot_lo, foot_hi) = path.foot_span(sa, sb);
    if foot_lo < -conf || foot_hi > reach + conf {
        return Err(refuse("the ball runs past that edge's end"));
    }
    let ends = [path.foot_at(sa), path.foot_at(sb)];
    if (foot_lo - ends[0].min(ends[1])).abs() > conf
        || (foot_hi - ends[0].max(ends[1])).abs() > conf
    {
        return Err(refuse("the ball's contact turns back along that edge"));
    }

    // The sections do not cross: the crease's crossing moves one way along
    // it, and each section's corners lie inside the curvature of the path,
    // so no point of one section lies on another's plane.
    let fine = STATIONS * 16;
    let mut last: Option<f64> = None;
    for i in 0..=fine {
        #[allow(clippy::cast_precision_loss)]
        let s = sa + (sb - sa) * (i as f64) / (fine as f64);
        let Some(x) = crossing(s) else {
            return Err(refuse("the ball's sections cross"));
        };
        if last.is_some_and(|l| x <= l) {
            return Err(refuse("the ball's sections cross"));
        }
        last = Some(x);
        let c = path.centre(s);
        let corners = [e0 + t * x, c - n1 * radius, path.foot(s)];
        let bend = corners
            .iter()
            .map(|p| path.d2(s).dot(*p - c))
            .fold(f64::NEG_INFINITY, f64::max);
        let speed = path.d1(s);
        if bend >= speed.dot(speed) {
            return Err(refuse("the ball's sections cross"));
        }
    }

    // The legs stand on their faces: the touched face holds the stretch of
    // its plane from the crease to the line of contact, and the narrow face
    // the stretch from the crease to the far edge. No edge of either face
    // but the crease and the far edge enters its stretch, and each face
    // lies on its stretch's side of the crease. The line of contact is held
    // by chords close enough that none strays from it by the confusion.
    let x_of = |s: f64| crossing(s).map(|x| e0 + t * x);
    let (Some(xa), Some(xb)) = (x_of(sa), x_of(sb)) else {
        return Err(refuse("the ball's sections cross"));
    };
    let into1 = -seat.leg(1 - narrow, tol)?;
    let into2 = -seat.leg(narrow, tol)?;
    let middle = (sa + sb) * 0.5;
    let Some(x_mid) = x_of(middle) else {
        return Err(refuse("the ball's sections cross"));
    };
    if (path.centre(middle) - n1 * radius - x_mid).dot(into1) <= 0.0
        || (path.foot(middle) - x_mid).dot(into2) <= 0.0
    {
        return Err(refuse("the round stands off its faces"));
    }
    let chords = {
        let bend = path.speed() * (sb - sa).abs();
        // A chord of `h` in the parameter strays at most `bend h^2 / 8`.
        let steps = ((bend * (sb - sa).abs() / (8.0 * conf)).sqrt()).ceil();
        if !steps.is_finite() || steps > MAX_CHORDS {
            return Err(refuse("the line of contact is too long to read"));
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let steps = (steps as usize).max(STATIONS);
        steps
    };
    let mut leg1 = vec![xa, xb];
    for i in (0..=chords).rev() {
        #[allow(clippy::cast_precision_loss)]
        let s = sa + (sb - sa) * (i as f64) / (chords as f64);
        leg1.push(path.centre(s) - n1 * radius);
    }
    let leg2 = vec![xa, xb, path.foot(sb), path.foot(sa)];
    for (face, ring, normal, skip) in [
        (f1, &leg1, n1, vec![edge.clone()]),
        (f2, &leg2, n2, vec![edge.clone(), far.clone()]),
    ] {
        for e in explore_unique(model, face, ShapeType::Edge)? {
            if skip.iter().any(|s| same_occurrence(model, &e, s, tol)) {
                continue;
            }
            let Some(side) = segment(model, &e, tol)? else {
                return Err(refuse("a face the round stands on has a curved edge"));
            };
            if enters(ring, normal, side, conf) {
                return Err(if face.is_same(f1) {
                    refuse("the ball's contact runs off the face it touches")
                } else {
                    refuse("the round's foot runs off the narrow face")
                });
            }
        }
    }

    // Every ball lies in free space: no face of the solid comes inside it.
    // A face on or behind the touched face's plane stands a radius off
    // every centre.
    // The narrow face lies on the crease's side of the far edge and the
    // beyond face on its own side, and the centre stands between their
    // normals, so each stands a radius off it. Every other face is read
    // at the stations, its distance falling by at most the centre's travel
    // between them, and an interval is halved until that travel fits under
    // its clearance.
    let speed = path.speed();
    let one_side = |flat: &Flat, inward: Vector| {
        flat.sides
            .iter()
            .flatten()
            .all(|p| (*p - base).dot(inward) >= -conf)
    };
    let to_narrow = {
        let w = n2.cross(along);
        let w = w / w.magnitude();
        if (e0 - base).dot(w) >= 0.0 { w } else { -w }
    };
    let to_beyond = far_seat.leg(k3, tol)?;
    for face in explore_unique(model, solid, ShapeType::Face)? {
        let Some(flat) = flat(model, &face, tol)? else {
            let bounds = ogeom_algo::shape_bounds(model, &face, tol)?;
            for i in 0..=STATIONS {
                let c = path.centre(station(i));
                let step = speed * (station(1) - station(0)).abs();
                if bounds.distance_to(c) <= radius + step {
                    return Err(refuse("a curved face stands within reach of the ball"));
                }
            }
            continue;
        };
        let behind_first = flat
            .sides
            .iter()
            .flatten()
            .all(|p| (*p - e0).dot(n1) <= conf);
        if behind_first
            || face.is_same(f2) && one_side(&flat, to_narrow)
            || face.is_same(f3) && one_side(&flat, to_beyond)
        {
            continue;
        }
        let clearance = |s: f64| face_distance(&flat, path.centre(s)) - radius;
        let mut intervals: Vec<(f64, f64, f64, f64, u32)> = Vec::new();
        for i in 0..STATIONS {
            let (a, b) = (station(i), station(i + 1));
            intervals.push((a, b, clearance(a), clearance(b), 0));
        }
        while let Some((a, b, ca, cb, depth)) = intervals.pop() {
            if ca < -conf || cb < -conf {
                return Err(refuse("the ball comes within its radius of another face"));
            }
            if ca.min(cb) >= speed * (b - a).abs() * 0.5 {
                continue;
            }
            if depth >= HALVINGS {
                return Err(refuse("the ball comes within its radius of another face"));
            }
            let m = (a + b) * 0.5;
            let cm = clearance(m);
            intervals.push((a, m, ca, cm, depth + 1));
            intervals.push((m, b, cm, cb, depth + 1));
        }
    }

    // The stations, ordered along the crease.
    let Some((first_host, first_sign)) = host(model, f1, tol)? else {
        return Ok(None);
    };
    let Some((second_host, second_sign)) = host(model, f2, tol)? else {
        return Ok(None);
    };
    let mut blend = MarchedBlend {
        spine: Vec::with_capacity(STATIONS + 1),
        on_first: Vec::with_capacity(STATIONS + 1),
        on_second: Vec::with_capacity(STATIONS + 1),
        touch_first: Vec::with_capacity(STATIONS + 1),
        touch_second: Vec::with_capacity(STATIONS + 1),
        along: Vec::with_capacity(STATIONS + 1),
        sides: Sides {
            first: 1,
            second: 1,
        },
        stopped: BlendStop::RanPastTheGuide,
    };
    for i in 0..=STATIONS {
        let s = station(i);
        let Some(x) = crossing(s) else {
            return Err(refuse("the ball's sections cross"));
        };
        let c = path.centre(s);
        let touch = c - n1 * radius;
        let foot = path.foot(s);
        blend.along.push(lo + (hi - lo) * x / length);
        blend.spine.push(c);
        blend.touch_first.push(touch);
        blend.touch_second.push(foot);
        blend
            .on_first
            .push(ogeom_algo::project_on_surface(&first_host, touch, 16, tol)?.parameters);
        blend
            .on_second
            .push(ogeom_algo::project_on_surface(&second_host, foot, 16, tol)?.parameters);
    }
    Ok(Some(Rolled {
        blend,
        guide: curve,
        hosts: [(first_host, first_sign), (second_host, second_sign)],
        rails: [Rail::Plane(e0, n1), Rail::Line(base, along)],
        path,
        span,
    }))
}

/// The span of the ellipse's angle, within `[0, phi3]` where the ball
/// rests on the edge alone, over which the arcs' planes cross the crease
/// from its start (`crossing` 0) to its end (`length`), one way along it
/// and without the crossing running off to infinity between; `None` where
/// no one such span exists.
fn ellipse_run(
    crossing: &dyn Fn(f64) -> Option<f64>,
    phi3: f64,
    length: f64,
) -> Option<(f64, f64)> {
    #[allow(clippy::cast_precision_loss)]
    let angle = |i: usize| phi3 * (i as f64) / (SEARCH as f64);
    let samples: Vec<Option<f64>> = (0..=SEARCH).map(|i| crossing(angle(i))).collect();
    // A root of `crossing - target` between two samples, by bisection.
    let root = |i: usize, target: f64| -> f64 {
        let (mut a, mut b) = (angle(i), angle(i + 1));
        let fa = samples[i].unwrap_or(0.0) - target;
        for _ in 0..80 {
            let m = (a + b) * 0.5;
            let fm = crossing(m).unwrap_or(0.0) - target;
            if (fm > 0.0) == (fa > 0.0) {
                a = m;
            } else {
                b = m;
            }
        }
        (a + b) * 0.5
    };
    // Each stretch of samples over which the crossing is finite and moves
    // one way, and the targets it passes.
    let mut found: Vec<(f64, f64)> = Vec::new();
    let mut i = 0;
    while i < SEARCH {
        let (Some(a), Some(b)) = (samples[i], samples[i + 1]) else {
            i += 1;
            continue;
        };
        let rising = b > a;
        let mut j = i;
        let mut ends: [Option<f64>; 2] = [None, None];
        while j < SEARCH {
            let (Some(p), Some(q)) = (samples[j], samples[j + 1]) else {
                break;
            };
            if (q > p) != rising || (q - p).abs() > length {
                break;
            }
            for (slot, target) in [(0, 0.0), (1, length)] {
                if (p - target) * (q - target) <= 0.0 && ends[slot].is_none() {
                    ends[slot] = Some(root(j, target));
                }
            }
            j += 1;
        }
        if let [Some(s0), Some(s1)] = ends {
            found.push((s0, s1));
        }
        i = j.max(i + 1);
    }
    match found.as_slice() {
        [one] => Some(*one),
        _ => None,
    }
}

/// Round `edge` with the ball resting on the far edge of the narrow face
/// beside it, where [`rolled_run`] finds that run; `None` where the edge
/// is not of that kind.
///
/// # Errors
///
/// As [`rolled_run`], and [`OgeomError::NotDone`] where the band's fit
/// stands further off the balls' envelope than its target.
pub(crate) fn rolled_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some(run) = rolled_run(model, solid, edge, radius, tol)? else {
        return Ok(None);
    };
    rolled_band(model, solid, edge, &run, radius, tol).map(Some)
}

/// The fill of a run found by [`rolled_run`], applied to the solid.
pub(crate) fn rolled_band(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    run: &Rolled,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let path = run.path;
    let (sa, sb) = run.span;
    let n = run.blend.len();
    // The distance from a point to the envelope: from its nearest centre on
    // the run, found by golden section about the nearest station.
    let off = |p: Point| -> f64 {
        let nearest = (0..n)
            .min_by(|&i, &j| {
                p.distance(run.blend.spine[i])
                    .total_cmp(&p.distance(run.blend.spine[j]))
            })
            .unwrap_or(0);
        #[allow(clippy::cast_precision_loss)]
        let at = |i: f64| sa + (sb - sa) * i / ((n - 1) as f64);
        #[allow(clippy::cast_precision_loss)]
        let (mut a, mut b) = (at(nearest as f64 - 1.0), at(nearest as f64 + 1.0));
        let f = |s: f64| p.distance(path.centre(s));
        let g = 0.5 * (5.0_f64.sqrt() - 1.0);
        let (mut x1, mut x2) = (b - g * (b - a), a + g * (b - a));
        let (mut f1, mut f2) = (f(x1), f(x2));
        for _ in 0..80 {
            if f1 < f2 {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = b - g * (b - a);
                f1 = f(x1);
            } else {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + g * (b - a);
                f2 = f(x2);
            }
        }
        (f1.min(f2) - radius).abs()
    };
    let exact = ExactBand {
        rails: run.rails,
        off: &off,
    };
    let [(first, first_sign), (second, second_sign)] = &run.hosts;
    crate::marched::build_held_band(
        model,
        solid,
        Some(edge),
        &run.blend,
        &run.guide,
        [(first, *first_sign), (second, *second_sign)],
        radius,
        false,
        [false, false],
        Some(&exact),
        tol,
    )
}
