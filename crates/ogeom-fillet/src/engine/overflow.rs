//! Blend into round (ADR 0001, amendment "overflow handling"): a round
//! whose ball, rolling on its floor, reaches an existing round of the
//! body rolls on over that round instead of trimming it; and a bevel whose
//! leg on the floor reaches such a round runs on over it.
//!
//! Built here for a boss's foot: the closed rim between a planar floor and
//! a wall of revolution about an axis square to it (a cylinder, or a cone
//! narrowing away from the floor), where the foot round's contact circle
//! crosses straight edges of the floor along which cylindrical rounds meet
//! the floor tangentially, convex or concave, over disjoint arcs. Where the
//! contact lies on the floor the round is the torus its ball sweeps; past a
//! tangent line the ball rests on the wall and that round, and its centre
//! follows the intersection of their offsets (the wall's offset about its
//! axis and the cylinder about the round's), so the pieces meet G1 along
//! the section where the ball touches the tangent line. The fill between
//! the round, the wall, the floor and the existing rounds is one closed
//! tool, fused into the solid.
//!
//! A bevel of the same foot (a cylindrical wall, one convex round) is the
//! cone between its two legs over the floor. Past the tangent line its leg
//! on the floor is where the wall's offset meets the existing round, and
//! its leg on the wall where the floor's offset meets the wall, so each
//! section is the straight ruling between them in the plane through the
//! wall's axis; the two pieces share that ruling where the leg crosses the
//! tangent line, and meet G1 there since the existing round leaves the
//! floor tangentially.

use core::f64::consts::{PI, TAU};

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{
    CircleCurve, ConeSurface, Curve, Curve3d as _, LineCurve, Surface as _, SurfaceGeometry,
    TorusSurface,
};
use ogeom_math::{Aabb, Circle, Cone, Direction, Frame, Point, Torus, Vector};
use ogeom_topo::{Model, Shape};

use super::curved::{Curved, SheetCurve, SheetRings, read};
use super::trimmed::{faces_of, outward};
use super::{Id, id};

/// Stations along each piece over an existing round.
const STATIONS: usize = 128;
/// Points across each of its sections.
const ACROSS: usize = 129;
/// Stations per arc where the fill's room and the balls are proved clear.
const PROOF: usize = 96;
/// The least angle between two crossings' arcs, so the torus between them
/// is a face.
const APART: f64 = 1e-2;
/// The wall's steepest lean from square to the floor, as the change of its
/// radius per unit of height.
const LEAN: f64 = 1.0;

/// What the blend sweeps along the foot.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Size {
    /// A round of this radius: the ball's arc between its contacts.
    Round(f64),
    /// A bevel, `out` along the floor and `rise` up the wall: the straight
    /// ruling between its legs.
    Bevel {
        /// From the wall along the floor.
        out: f64,
        /// Up the wall from the floor.
        rise: f64,
    },
}

impl Size {
    /// How far the blend reaches along the floor from the wall.
    const fn out(self) -> f64 {
        match self {
            Self::Round(r) => r,
            Self::Bevel { out, .. } => out,
        }
    }

    /// How far it reaches up the wall over the floor.
    const fn rise(self) -> f64 {
        match self {
            Self::Round(r) => r,
            Self::Bevel { rise, .. } => rise,
        }
    }
}

/// One straight edge of the floor whose round the contact circle crosses.
#[derive(Debug, Clone)]
struct Crossing {
    line: Id,
    round: usize,
    /// Square to the tangent line, from the boss toward it, in the floor.
    across: Vector,
    /// The angle of `across` about the wall's axis.
    angle: f64,
    /// From the foot circle's centre to the tangent line.
    gap: f64,
    /// The existing round's radius.
    round_radius: f64,
    /// Whether the round is concave: its axis stands above the floor.
    concave: bool,
}

/// A foot blend reaching straight rounds tangent to its floor.
#[derive(Debug, Clone)]
pub(crate) struct Reach {
    curved: Curved,
    edges: Vec<Id>,
    wall: usize,
    floor: usize,
    /// The foot circle's centre, on the floor.
    base: Point,
    /// The floor's outward normal.
    up: Vector,
    /// The first crossing's `across`: angle zero about the axis.
    ex: Vector,
    /// `up x ex`.
    ey: Vector,
    /// The foot circle's radius.
    foot_radius: f64,
    /// The wall's radius change per unit of height: zero for a cylinder,
    /// negative for a cone narrowing away from the floor.
    slope: f64,
    /// By angle, over disjoint arcs.
    crossings: Vec<Crossing>,
    /// The foot blend's size.
    size: Size,
    /// How far it reaches along the floor: the round's radius, or the
    /// bevel's leg there.
    radius: f64,
}

/// `a` brought into `(-pi, pi]`.
fn wrapped(a: f64) -> f64 {
    let w = a.rem_euclid(TAU);
    if w > PI { w - TAU } else { w }
}

impl Reach {
    fn radial(&self, theta: f64) -> Vector {
        self.ex * theta.cos() + self.ey * theta.sin()
    }

    /// The ball centre's distance from the axis at height `z`: on the
    /// wall's offset. A bevel's leg on the floor stands `out` from the wall.
    fn spine_radius(&self, z: f64) -> f64 {
        match self.size {
            Size::Round(r) => self
                .slope
                .mul_add(z, r.mul_add(self.slope.hypot(1.0), self.foot_radius)),
            Size::Bevel { out, .. } => self.foot_radius + out,
        }
    }

    /// The contact circle's radius on the floor.
    fn reach_radius(&self) -> f64 {
        self.spine_radius(self.radius)
    }

    /// The angle either side of a crossing's `across` within which the
    /// ball rests on its round.
    fn half_angle(&self, c: &Crossing) -> f64 {
        (c.gap / self.reach_radius()).acos()
    }

    /// The arc of `c` about the axis: from its first seat change to its
    /// second.
    fn span(&self, c: &Crossing) -> (f64, f64) {
        let h = self.half_angle(c);
        (c.angle - h, c.angle + h)
    }

    /// The crossing whose arc holds `theta`, and the angle from its middle.
    fn crossing_at(&self, theta: f64) -> Option<(&Crossing, f64)> {
        self.crossings.iter().find_map(|c| {
            let d = wrapped(theta - c.angle);
            (d.abs() < self.half_angle(c)).then_some((c, d))
        })
    }

    /// The ball centre's height over a round `s` past its tangent line.
    fn height_on(&self, c: &Crossing, s: f64) -> f64 {
        let r = self.radius;
        if s <= 0.0 {
            return r;
        }
        let ra = c.round_radius;
        if c.concave {
            ra - (ra - r).mul_add(ra - r, -s * s).max(0.0).sqrt()
        } else {
            (ra + r).mul_add(ra + r, -s * s).max(0.0).sqrt() - ra
        }
    }

    /// The ball centre's height above the floor at `theta`, where the
    /// wall's offset and the floor's meet, or the height of the bevel's leg
    /// on the wall; with the residual of the solve.
    fn solved_height(&self, theta: f64) -> (f64, f64) {
        if let Size::Bevel { rise, .. } = self.size {
            return (rise, 0.0);
        }
        let Some((c, d)) = self.crossing_at(theta) else {
            return (self.radius, 0.0);
        };
        // The offset of the wall leans by `slope`, so the height is a fixed
        // point of the round's offset over the spine radius it gives.
        let at = |z: f64| self.height_on(c, self.spine_radius(z).mul_add(d.cos(), -c.gap));
        let mut z = self.radius;
        for _ in 0..200 {
            let next = at(z);
            let step = (next - z).abs();
            z = next;
            if step <= f64::EPSILON * (1.0 + z.abs()) {
                break;
            }
        }
        (z, (at(z) - z).abs())
    }

    fn height(&self, theta: f64) -> f64 {
        self.solved_height(theta).0
    }

    fn centre(&self, theta: f64) -> Point {
        let z = self.height(theta);
        self.base + self.radial(theta) * self.spine_radius(z) + self.up * z
    }

    /// The wall's normal into the fill at `theta`.
    fn wall_normal(&self, theta: f64) -> Vector {
        (self.radial(theta) - self.up * self.slope) / self.slope.hypot(1.0)
    }

    fn wall_touch(&self, theta: f64) -> Point {
        match self.size {
            Size::Round(r) => self.centre(theta) - self.wall_normal(theta) * r,
            Size::Bevel { rise, .. } => {
                self.base + self.radial(theta) * self.foot_radius + self.up * rise
            }
        }
    }

    /// The point of a crossing's tangent line under the ball's centre
    /// `centre`, and the round's axis below or above it.
    fn round_axis(&self, c: &Crossing, centre: Point) -> (Point, Point) {
        let line = self.base + c.across * c.gap;
        let along = self.up.cross(c.across);
        let under = line + along * (centre - line).dot(along);
        let side = if c.concave { 1.0 } else { -1.0 };
        (under, under + self.up * (side * c.round_radius))
    }

    /// Where the ball touches the floor or the existing round at `theta`,
    /// or where the bevel's leg there lies: on the wall's offset, straight
    /// below `on_line` on the floor or on the existing round.
    fn floor_touch(&self, theta: f64) -> Point {
        let found = self.crossing_at(theta);
        if let Size::Bevel { .. } = self.size {
            let Some((c, d)) = found else {
                return self.on_line(theta);
            };
            let s = self.reach_radius().mul_add(d.cos(), -c.gap);
            let r = c.round_radius;
            return self.on_line(theta) + self.up * ((r * r - s * s).sqrt() - r);
        }
        let centre = self.centre(theta);
        let Some((c, _)) = found else {
            return centre - self.up * self.radius;
        };
        let (_, axis) = self.round_axis(c, centre);
        let w = centre - axis;
        axis + w * (c.round_radius / w.magnitude())
    }

    fn foot(&self, theta: f64) -> Point {
        self.base + self.radial(theta) * self.foot_radius
    }

    /// The point where the ball's contact meets a tangent line at a seat
    /// change: on the contact circle.
    fn on_line(&self, theta: f64) -> Point {
        self.base + self.radial(theta) * self.reach_radius()
    }

    /// The section at `theta`, from the wall's touch to the floor's: the
    /// ball's arc or the bevel's ruling, at `n` points.
    fn section(&self, theta: f64, n: usize, tol: Tolerances) -> Option<Vec<Point>> {
        let (w, f) = (self.wall_touch(theta), self.floor_touch(theta));
        #[expect(clippy::cast_precision_loss, reason = "a few points")]
        let at = |k: usize| k as f64 / (n - 1) as f64;
        if let Size::Bevel { .. } = self.size {
            return Some((0..n).map(|k| w + (f - w) * at(k)).collect());
        }
        let c = self.centre(theta);
        let a = (w - c) / self.radius;
        let b = (f - c) / self.radius;
        let cross = a.cross(b);
        let m = cross.magnitude();
        if m <= tol.angular() {
            return None;
        }
        let axis = cross / m;
        let sweep = a.dot(b).clamp(-1.0, 1.0).acos();
        Some(
            (0..n)
                .map(|k| {
                    let t = sweep * at(k);
                    c + (a * t.cos() + axis.cross(a) * t.sin()) * self.radius
                })
                .collect(),
        )
    }

    /// A point inside the blend's face at `theta`, and the fill's outward
    /// direction there: away from the corner the fill stands in.
    fn blend_mid(&self, theta: f64) -> (Point, Vector) {
        let (w, f) = (self.wall_touch(theta), self.floor_touch(theta));
        if let Size::Bevel { .. } = self.size {
            let d = f - w;
            let radial = self.radial(theta);
            let n = radial * -d.dot(self.up) + self.up * d.dot(radial);
            let n = n / n.magnitude();
            let mid = w + d * 0.5;
            let away = mid - self.foot(theta);
            return (mid, if n.dot(away) < 0.0 { -n } else { n });
        }
        let c = self.centre(theta);
        let m = (w - c) + (f - c);
        let m = m / m.magnitude();
        (c + m * self.radius, -m)
    }

    /// The angles at which the fill's room and the balls are proved: each
    /// crossing's arc and each arc between, `PROOF` steps apiece, in
    /// order round the axis.
    fn proof_angles(&self) -> Vec<f64> {
        let mut out = Vec::new();
        let n = self.crossings.len();
        for k in 0..n {
            let (a, b) = self.span(&self.crossings[k]);
            let (next, _) = self.span(&self.crossings[(k + 1) % n]);
            let around = (next - b).rem_euclid(TAU);
            for (from, sweep) in [(a, b - a), (b, around)] {
                for i in 0..PROOF {
                    #[expect(clippy::cast_precision_loss, reason = "a hundred stations")]
                    out.push(sweep.mul_add(i as f64 / PROOF as f64, from));
                }
            }
        }
        out
    }
}

/// The plane face and the wall's face (a cylinder, or for a round a cone
/// too) `faces` holds, in that order.
fn plane_and_wall(curved: &Curved, faces: &[usize], cone: bool) -> Option<(usize, usize)> {
    let [a, b] = faces[..] else {
        return None;
    };
    let wall = |s: &SurfaceGeometry| {
        matches!(s, SurfaceGeometry::Cylinder(_)) || (cone && matches!(s, SurfaceGeometry::Cone(_)))
    };
    match (&curved.faces[a].surface, &curved.faces[b].surface) {
        (SurfaceGeometry::Plane(_), s) if wall(s) => Some((a, b)),
        (s, SurfaceGeometry::Plane(_)) if wall(s) => Some((b, a)),
        _ => None,
    }
}

/// `edges` read as a boss's foot rim blended at `size` whose blend
/// reaches straight rounds tangent to the floor; `None` for any other
/// request, or where anything else of the solid stands within reach of the
/// fill or of the balls.
///
/// # Errors
///
/// As reading the solid.
#[expect(clippy::too_many_lines, reason = "one reading of the foot, in order")]
pub(crate) fn reach(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    size: Size,
    tol: Tolerances,
) -> OgeomResult<Option<Reach>> {
    let eps = tol.confusion();
    let ang = tol.angular();
    let (radius, rise) = (size.out(), size.rise());
    if edges.is_empty() || [radius, rise].iter().any(|d| !d.is_finite() || *d <= eps) {
        return Ok(None);
    }
    let round = matches!(size, Size::Round(_));
    // A plane and a wall either side of the first edge, before the solid is
    // read.
    let kinds: Vec<Option<SurfaceGeometry>> =
        ogeom_topo::ancestors_of(model, solid, &edges[0], ogeom_topo::ShapeType::Face)?
            .iter()
            .map(|f| {
                model
                    .node(f)
                    .and_then(|n| n.data().as_face())
                    .and_then(|d| model.geometry().surface(d.surface))
                    .cloned()
            })
            .collect();
    let count = |want: fn(&SurfaceGeometry) -> bool| {
        kinds
            .iter()
            .filter(|k| k.as_ref().is_some_and(want))
            .count()
    };
    let plane = count(|s| matches!(s, SurfaceGeometry::Plane(_)));
    let walls = if round {
        count(|s| matches!(s, SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_)))
    } else {
        count(|s| matches!(s, SurfaceGeometry::Cylinder(_)))
    };
    if kinds.len() != 2 || plane != 1 || walls != 1 {
        return Ok(None);
    }
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    let mut asked: Vec<Id> = edges.iter().map(id).collect();
    asked.sort_unstable();
    asked.dedup();
    if asked.iter().any(|e| !curved.edges.contains_key(e)) {
        return Ok(None);
    }
    let Some((floor, wall)) = plane_and_wall(&curved, &faces_of(&curved, asked[0]), round) else {
        return Ok(None);
    };
    // Every edge the floor and the wall share is asked, and nothing else.
    let mut shared: Vec<Id> = curved
        .edges
        .keys()
        .copied()
        .filter(|e| {
            let mut f = faces_of(&curved, *e);
            f.sort_unstable();
            let mut pair = vec![floor, wall];
            pair.sort_unstable();
            f == pair
        })
        .collect();
    shared.sort_unstable();
    if shared != asked {
        return Ok(None);
    }
    let SurfaceGeometry::Plane(plane) = &curved.faces[floor].surface else {
        return Ok(None);
    };
    let plane = plane.plane();
    let up = plane.normal().vector()
        * if curved.faces[floor].reversed {
            -1.0
        } else {
            1.0
        };
    // The wall's axis, and its radius as a line in the height over the floor.
    let (axis, at_floor, slope) = match &curved.faces[wall].surface {
        SurfaceGeometry::Cylinder(c) => (c.cylinder().axis(), None, 0.0),
        SurfaceGeometry::Cone(c) => {
            let cone = c.cone();
            let z = cone.frame().z().vector();
            let lean = cone.half_angle().tan();
            let base = cone.axis().location;
            let base = base - up * (base - plane.origin()).dot(up);
            let along = (base - cone.frame().origin()).dot(z);
            (
                cone.axis(),
                Some(lean.mul_add(along, cone.reference_radius())),
                lean * z.dot(up),
            )
        }
        _ => return Ok(None),
    };
    if axis.direction.vector().cross(up).magnitude() > ang || !(-LEAN..=0.0).contains(&slope) {
        return Ok(None);
    }
    let foot_radius = match (&curved.faces[wall].surface, at_floor) {
        (SurfaceGeometry::Cylinder(c), _) => c.cylinder().radius(),
        (_, Some(r)) => r,
        _ => return Ok(None),
    };
    if foot_radius <= eps * 10.0 {
        return Ok(None);
    }
    let base = axis.location - up * (axis.location - plane.origin()).dot(up);
    // The foot circle closes round the axis on the floor: every edge lies
    // on it, and together they turn once round. An open arc of a foot is
    // not taken.
    let mut turned = 0.0;
    for e in &asked {
        let edge = &curved.edges[e];
        let mut last: Option<Vector> = None;
        for k in 0..=16 {
            let t = edge.range.0 + (edge.range.1 - edge.range.0) * f64::from(k) / 16.0;
            let p = edge.curve.point_at(t, tol)?;
            let q = p - base;
            if q.dot(up).abs() > eps * 10.0 || (q.magnitude() - foot_radius).abs() > eps * 10.0 {
                return Ok(None);
            }
            if let Some(l) = last {
                turned += l.cross(q).magnitude().atan2(l.dot(q));
            }
            last = Some(q);
        }
    }
    if (turned - TAU).abs() > 1e-6 {
        return Ok(None);
    }
    // A boss: the wall's material inside its axis' reach, standing above
    // the floor.
    let some_foot = curved.edges[&asked[0]]
        .curve
        .point_at(curved.edges[&asked[0]].range.0, tol)?;
    let out = outward(&curved, wall, some_foot).unwrap_or(Vector::ZERO);
    if out.dot(some_foot - base) <= 0.0 {
        return Ok(None);
    }
    let offset = match size {
        Size::Round(r) => r * slope.hypot(1.0) + slope * r,
        Size::Bevel { out, .. } => out,
    };
    let reach_radius = foot_radius + offset;
    // The floor's straight edges along which a round meets it tangentially,
    // those the contact circle crosses.
    let mut crossings: Vec<Crossing> = Vec::new();
    for (e, _) in curved.faces[floor].loops.iter().flatten() {
        if asked.contains(e) {
            continue;
        }
        let edge = &curved.edges[e];
        let Curve::Line(line) = &edge.curve else {
            continue;
        };
        let others: Vec<usize> = faces_of(&curved, *e)
            .into_iter()
            .filter(|f| *f != floor)
            .collect();
        let [round_face] = others[..] else {
            continue;
        };
        let SurfaceGeometry::Cylinder(c) = &curved.faces[round_face].surface else {
            continue;
        };
        let along = line.axis().direction.vector();
        let round_axis = c.cylinder().axis();
        if round_axis.direction.vector().cross(along).magnitude() > ang {
            continue;
        }
        let round_radius = c.cylinder().radius();
        let p = line.axis().location;
        let off = p - round_axis.location;
        let off = off - along * off.dot(along);
        let concave = if (off - up * round_radius).magnitude() <= eps * 10.0 {
            false
        } else if round && (off + up * round_radius).magnitude() <= eps * 10.0 {
            true
        } else {
            continue;
        };
        if outward(&curved, round_face, p).is_none_or(|n| n.dot(up) < 1.0 - ang) {
            continue;
        }
        let to = p - base;
        let to = to - up * to.dot(up);
        let to = to - along * to.dot(along);
        let gap = to.magnitude();
        if gap <= foot_radius + eps * 10.0 || gap >= reach_radius - eps * 10.0 {
            continue;
        }
        let past = reach_radius - gap;
        // A bevel's leg stays on the existing round, short of where it
        // turns down from the floor; a ball resting in a concave round
        // stays inside its cylinder, clear of where it ends.
        let seated = match size {
            Size::Round(r) => {
                !concave || (r < round_radius && past < round_radius - r - eps * 10.0)
            }
            Size::Bevel { .. } => past < round_radius - eps * 10.0,
        };
        if !seated {
            return Ok(None);
        }
        let across = to / gap;
        crossings.push(Crossing {
            line: *e,
            round: round_face,
            across,
            angle: 0.0,
            gap,
            round_radius,
            concave,
        });
    }
    if crossings.is_empty() || (!round && (crossings.len() > 1 || crossings[0].concave)) {
        return Ok(None);
    }
    let ex = crossings[0].across;
    let ey = up.cross(ex);
    for c in &mut crossings {
        c.angle = c.across.dot(ey).atan2(c.across.dot(ex));
    }
    crossings.sort_by(|a, b| a.angle.total_cmp(&b.angle));
    let reach = Reach {
        curved,
        edges: asked,
        wall,
        floor,
        base,
        up,
        ex,
        ey,
        foot_radius,
        slope,
        crossings,
        size,
        radius,
    };
    // The crossings' arcs stand apart, so a torus lies between each two.
    let n = reach.crossings.len();
    if n > 1 {
        for k in 0..n {
            let (_, b) = reach.span(&reach.crossings[k]);
            let (a, _) = reach.span(&reach.crossings[(k + 1) % n]);
            if (a - b).rem_euclid(TAU) < APART || (a - b).rem_euclid(TAU) > TAU - APART {
                return Ok(None);
            }
        }
    }
    // The centres solve where the offsets meet, and the wall stands above
    // every touch on it.
    let mut highest = rise;
    for theta in reach.proof_angles() {
        let (_, residual) = reach.solved_height(theta);
        if residual > eps * 1e-3 {
            return Ok(None);
        }
        highest = highest.max((reach.wall_touch(theta) - base).dot(up));
    }
    let c = &reach.curved;
    let above = c.faces[wall].loops.iter().flatten().any(|(e, _)| {
        let edge = &c.edges[e];
        [edge.ends.0, edge.ends.1]
            .iter()
            .any(|v| (c.points[v] - base).dot(up) > highest + eps * 10.0)
    });
    if !above || !clear(model, &reach, tol)? {
        return Ok(None);
    }
    Ok(Some(reach))
}

/// Whether nothing of the solid but the wall, the floor and the rounds
/// crossed stands within the room the fill can take up, and no other face
/// comes within the radius of a ball's centre.
///
/// The room is the union of boxes, one between each two stations round
/// the axis, each holding both stations' sections, their feet on the wall
/// and the floor under them, padded by the chord's sag. A face whose box
/// meets a ball's is measured from each centre near it.
fn clear(model: &Model, reach: &Reach, tol: Tolerances) -> OgeomResult<bool> {
    let eps = tol.confusion();
    let up = reach.up;
    let margin = reach.radius * 1e-3 + eps * 10.0;
    let hi = reach.size.rise().max(reach.radius);
    let angles = reach.proof_angles();
    let n = angles.len();
    // Each station's points and, for a round, its ball's centre.
    let mut stations: Vec<(Vec<Point>, Option<Point>)> = Vec::with_capacity(n);
    for &theta in &angles {
        let (foot, wall) = (reach.foot(theta), reach.wall_touch(theta));
        let floor = reach.floor_touch(theta);
        let top = hi.max((wall - reach.base).dot(up));
        let mut points = vec![foot, foot + up * top, wall, floor, reach.on_line(theta)];
        points.push(reach.on_line(theta) + up * top);
        if let Some((c, d)) = reach.crossing_at(theta) {
            // Where the radial crosses the tangent line.
            points.push(reach.base + reach.radial(theta) * (c.gap / d.cos()));
        }
        if let Some(arc) = reach.section(theta, 9, tol) {
            points.extend(arc);
        }
        let centre = matches!(reach.size, Size::Round(_)).then(|| reach.centre(theta));
        stations.push((points, centre));
    }
    let mut room: Vec<Aabb> = Vec::with_capacity(n);
    let mut balls: Vec<(Aabb, Point, f64)> = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let step = (angles[j] - angles[i]).rem_euclid(TAU);
        let far = reach.reach_radius() + reach.radius;
        let sag = far * (1.0 - (0.5 * step).cos()) + margin;
        let mut points = stations[i].0.clone();
        points.extend(stations[j].0.iter().copied());
        let b = Aabb::of_points(&points);
        if let (Some(lo), Some(hi)) = (b.low(), b.high()) {
            let pad = Vector::new(sag, sag, sag);
            room.push(Aabb::of_points(&[lo - pad, hi + pad]));
        }
        if let (Some(a), Some(b)) = (stations[i].1, stations[j].1) {
            // Every ball between the two lies within their spacing of the
            // first.
            let half = a.distance(b) + margin;
            let reachable = reach.radius + half;
            let pad = Vector::new(reachable, reachable, reachable);
            balls.push((
                Aabb::of_points(&[a - pad, a + pad, b - pad, b + pad]),
                a,
                half,
            ));
        }
    }
    let meets = |bounds: &Aabb| room.iter().any(|r| r.intersects(bounds));
    let c = &reach.curved;
    let seated: Vec<usize> = [reach.wall, reach.floor]
        .into_iter()
        .chain(reach.crossings.iter().map(|x| x.round))
        .collect();
    let lines: Vec<Id> = reach.crossings.iter().map(|x| x.line).collect();
    for (fi, face) in c.faces.iter().enumerate() {
        if seated.contains(&fi) {
            continue;
        }
        let bounds = ogeom_algo::shape_bounds(model, &face.shape, tol)?;
        if meets(&bounds) {
            return Ok(false);
        }
        let near: Vec<&(Aabb, Point, f64)> = balls
            .iter()
            .filter(|(b, _, _)| b.intersects(&bounds))
            .collect();
        if near.is_empty() {
            continue;
        }
        let measure = ogeom_algo::PointDistance::of(model, &face.shape, tol)?;
        for (_, centre, half) in near {
            let (d, _) = measure.nearest(*centre, tol)?;
            if d < reach.radius + half {
                return Ok(false);
            }
        }
    }
    for fi in seated {
        for (e, _) in c.faces[fi].loops.iter().flatten() {
            if reach.edges.contains(e) || lines.contains(e) {
                continue;
            }
            // A seam of the wall is the wall's own.
            if faces_of(c, *e) == [reach.wall] {
                continue;
            }
            if meets(&ogeom_algo::shape_bounds(model, &c.edges[e].shape, tol)?) {
                return Ok(false);
            }
        }
    }
    for line in &lines {
        let edge = &c.edges[line];
        for v in [edge.ends.0, edge.ends.1] {
            if room.iter().any(|r| r.contains(c.points[&v])) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}
/// The arc of `radius` about `centre` in the plane square to `z`, from the
/// direction `x` through `sweep`.
pub(super) fn arc(
    centre: Point,
    z: Vector,
    x: Vector,
    radius: f64,
    sweep: f64,
    tol: Tolerances,
) -> OgeomResult<(Curve, (f64, f64))> {
    let frame = Frame::new(
        centre,
        Direction::new(z, tol)?,
        Direction::new(x, tol)?,
        tol,
    )?;
    Ok((
        CircleCurve::new(Circle::new(frame, radius, tol)?).into(),
        (0.0, sweep),
    ))
}

pub(super) fn segment(a: Point, b: Point, tol: Tolerances) -> OgeomResult<(Curve, (f64, f64))> {
    Ok((
        Curve::Line(LineCurve::segment(a, b, tol)?),
        (0.0, a.distance(b)),
    ))
}

/// The piece of the blend over a crossing's round: its section at each
/// station from its wall touch to its touch on the round, fitted.
fn strip(
    reach: &Reach,
    crossing: &Crossing,
    tol: Tolerances,
) -> OgeomResult<Option<ogeom_geom::BSplineSurface>> {
    let (a, b) = reach.span(crossing);
    let mut stations: Vec<Vec<Point>> = Vec::with_capacity(STATIONS + 1);
    for i in 0..=STATIONS {
        #[expect(clippy::cast_precision_loss, reason = "a hundred stations")]
        let theta = (b - a).mul_add(i as f64 / STATIONS as f64, a);
        let Some(row) = reach.section(theta, ACROSS, tol) else {
            return Ok(None);
        };
        stations.push(row);
    }
    let target = tol.confusion() * 0.05;
    let Ok(fitted) = ogeom_geom::fit::fit_surface_grid_chordal(&stations, 3, target, tol) else {
        return Ok(None);
    };
    Ok((fitted.error <= tol.confusion()).then_some(fitted.curve))
}
/// The boundary iso-curve of `surface` between `points[ends.0]` and
/// `points[ends.1]`, as a sheet curve over those ends.
pub(super) fn boundary_from(
    surface: &ogeom_geom::BSplineSurface,
    points: &[Point],
    ends: (usize, usize),
    tol: Tolerances,
) -> OgeomResult<Option<SheetCurve>> {
    let (from, to) = (points[ends.0], points[ends.1]);
    let ((u0, u1), (v0, v1)) = surface.domain();
    let mut candidates = Vec::with_capacity(4);
    for u in [u0, u1] {
        candidates.push(surface.iso_u_curve(u, tol)?);
    }
    for v in [v0, v1] {
        candidates.push(surface.iso_v_curve(v, tol)?);
    }
    let reach = tol.confusion() * 10.0;
    for c in candidates {
        let (a, b) = c.domain();
        let (pa, pb) = (c.point_at(a, tol)?, c.point_at(b, tol)?);
        if pa.distance(from) <= reach && pb.distance(to) <= reach {
            return Ok(Some(SheetCurve {
                ends,
                curve: c.into(),
                range: (a, b),
            }));
        }
        if pb.distance(from) <= reach && pa.distance(to) <= reach {
            return Ok(Some(SheetCurve {
                ends: (ends.1, ends.0),
                curve: c.into(),
                range: (a, b),
            }));
        }
    }
    Ok(None)
}

/// The parameters' periods of `surface`, where they close.
fn periods(surface: &SurfaceGeometry) -> (Option<f64>, Option<f64>) {
    let turn = core::f64::consts::TAU;
    match surface {
        SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_) => (Some(turn), None),
        SurfaceGeometry::Torus(_) => (Some(turn), Some(turn)),
        SurfaceGeometry::Sphere(_) => (Some(turn), None),
        _ => (None, None),
    }
}

/// One face of a tool: its surface, its loop by curve, and a point of it
/// with the tool's outward direction there.
pub(super) struct ToolFace {
    pub(super) surface: SurfaceGeometry,
    pub(super) curves: Vec<usize>,
    pub(super) at: Point,
    pub(super) out: Vector,
}

/// The closed solid `faces` bound over `points` and `curves`, each loop
/// walked so its face lies on its left about the surface's normal (by the
/// loop's signed area in the surface's parameters), each face turned so
/// it faces `out`. `None` where a loop will not close, the shell is open,
/// or the solid holds no volume.
pub(super) fn closed(
    model: &mut Model,
    points: &[Point],
    curves: &[SheetCurve],
    faces: &[ToolFace],
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
    let mut sheet: Vec<SheetRings> = Vec::with_capacity(faces.len());
    for face in faces {
        let first = face.curves[0];
        let mut walk = vec![(first, true)];
        let mut at = curves[first].ends.1;
        let mut left: Vec<usize> = face.curves[1..].to_vec();
        while !left.is_empty() {
            let Some(k) = left
                .iter()
                .position(|c| curves[*c].ends.0 == at || curves[*c].ends.1 == at)
            else {
                return Ok(None);
            };
            let c = left.remove(k);
            let forward = curves[c].ends.0 == at;
            at = if forward {
                curves[c].ends.1
            } else {
                curves[c].ends.0
            };
            walk.push((c, forward));
        }
        if at != curves[first].ends.0 {
            return Ok(None);
        }
        // The loop's signed area in the surface's parameters, each period
        // unwrapped along the walk.
        let (pu, pv) = periods(&face.surface);
        let mut uv: Vec<(f64, f64)> = Vec::new();
        for (c, forward) in &walk {
            let (lo, hi) = curves[*c].range;
            for k in 0..24 {
                let s = f64::from(k) / 24.0;
                let t = if *forward {
                    lo + (hi - lo) * s
                } else {
                    hi - (hi - lo) * s
                };
                let p = curves[*c].curve.point_at(t, tol)?;
                let foot = ogeom_geom::project::project_on_surface(&face.surface, p, 16, tol)?;
                let (mut u, mut v) = foot.parameters;
                if let Some(&(lu, lv)) = uv.last() {
                    if let Some(period) = pu {
                        u += period * ((lu - u) / period).round();
                    }
                    if let Some(period) = pv {
                        v += period * ((lv - v) / period).round();
                    }
                }
                uv.push((u, v));
            }
        }
        let n = uv.len();
        let area: f64 = (0..n)
            .map(|i| {
                let (a, b) = (uv[i], uv[(i + 1) % n]);
                a.0 * b.1 - b.0 * a.1
            })
            .sum();
        if area.abs() <= f64::EPSILON {
            return Ok(None);
        }
        if area < 0.0 {
            walk = walk.iter().rev().map(|(c, f)| (*c, !f)).collect();
        }
        let foot = ogeom_geom::project::project_on_surface(&face.surface, face.at, 16, tol)?;
        let (u, v) = foot.parameters;
        let (du, dv) = face.surface.d1_at(u, v, tol)?;
        sheet.push(SheetRings {
            surface: face.surface.clone(),
            rings: vec![walk],
            reversed: du.cross(dv).dot(face.out) < 0.0,
        });
    }
    let Ok(built) = super::curved::build(model, points, curves, &sheet, tol) else {
        return Ok(None);
    };
    let shell = ogeom_algo::make_shell(model, &built.faces)?.shape;
    if !ogeom_algo::is_shell_closed(model, &shell)? {
        return Ok(None);
    }
    let solid = ogeom_algo::make_solid(model, &[shell])?.shape;
    let volume =
        ogeom_algo::volume_properties(model, &solid, ogeom_mesh::Deflection::default(), tol)?.mass;
    if volume <= 0.0 {
        return Ok(None);
    }
    Ok(Some((solid, built.faces)))
}
/// The blend's face over the floor: the torus the ball sweeps, or the cone
/// between the bevel's legs.
fn over_floor(reach: &Reach, tol: Tolerances) -> OgeomResult<SurfaceGeometry> {
    let (r, rr, up) = (reach.radius, reach.foot_radius, reach.up);
    let Size::Bevel { rise, .. } = reach.size else {
        return Ok(TorusSurface::new(Torus::new(
            Frame::new(
                reach.base + up * r,
                Direction::new(up, tol)?,
                Direction::new(reach.ex, tol)?,
                tol,
            )?,
            reach.reach_radius(),
            r,
            tol,
        )?)
        .into());
    };
    // The cone's radius falls from `rr + r` on the floor to `rr` at the
    // leg on the wall; its apex stands past that leg.
    let frame = Frame::new(
        reach.base,
        Direction::new(up, tol)?,
        Direction::new(reach.ex, tol)?,
        tol,
    )?;
    let apex = rise * (rr + r) / r;
    let pad = (0.1 * rise).min(0.5 * (apex - rise));
    let cone = Cone::new(frame, rr + r, (-r / rise).atan(), tol)?;
    Ok(ConeSurface::new(cone, (-pad, rise + pad))?.into())
}

/// The curves at one seat change: up the wall, out along the floor, and
/// the section between them; with the points they join.
#[derive(Clone, Copy)]
struct Seat {
    foot: usize,
    wall: usize,
    line: usize,
    up: usize,
    ray: usize,
    section: usize,
}

/// The fill of `reach` as one closed tool, with its blend faces.
#[expect(clippy::too_many_lines, reason = "one tool, face by face")]
fn tool(
    model: &mut Model,
    reach: &Reach,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
    let r = reach.radius;
    let rr = reach.foot_radius;
    let big = reach.reach_radius();
    let up = reach.up;
    let n = reach.crossings.len();
    let mut points: Vec<Point> = Vec::new();
    let mut curves: Vec<SheetCurve> = Vec::new();
    let push = |curves: &mut Vec<SheetCurve>,
                ends: (usize, usize),
                (curve, range): (Curve, (f64, f64))| {
        curves.push(SheetCurve { ends, curve, range });
        curves.len() - 1
    };
    // The seats at each crossing's two seat changes.
    let mut seats: Vec<[Seat; 2]> = Vec::with_capacity(n);
    for c in &reach.crossings {
        let (a, b) = reach.span(c);
        let mut pair = [Seat {
            foot: 0,
            wall: 0,
            line: 0,
            up: 0,
            ray: 0,
            section: 0,
        }; 2];
        for (k, theta) in [a, b].into_iter().enumerate() {
            let at = points.len();
            points.extend([
                reach.foot(theta),
                reach.wall_touch(theta),
                reach.on_line(theta),
            ]);
            let (f, w, p) = (at, at + 1, at + 2);
            let radial = reach.radial(theta);
            let up_curve = push(&mut curves, (f, w), segment(points[f], points[w], tol)?);
            let ray = push(&mut curves, (f, p), segment(points[f], points[p], tol)?);
            let section = match reach.size {
                Size::Round(_) => {
                    let centre = reach.centre(theta);
                    let x = (points[w] - centre) / r;
                    let sweep = x.dot(-up).clamp(-1.0, 1.0).acos();
                    arc(centre, radial.cross(up), x, r, sweep, tol)?
                }
                Size::Bevel { .. } => segment(points[w], points[p], tol)?,
            };
            let section = push(&mut curves, (w, p), section);
            pair[k] = Seat {
                foot: f,
                wall: w,
                line: p,
                up: up_curve,
                ray,
                section,
            };
        }
        seats.push(pair);
    }
    let curved = &reach.curved;
    let wall = curved.faces[reach.wall].surface.clone();
    let floor = curved.faces[reach.floor].surface.clone();
    let blend = over_floor(reach, tol)?;
    let mut faces: Vec<ToolFace> = Vec::new();
    let mut blends: Vec<usize> = Vec::new();
    // Round the floor, from each crossing to the next: the foot, the
    // wall's touch and the floor's touch.
    for k in 0..n {
        let (_, from) = reach.span(&reach.crossings[k]);
        let (to, _) = reach.span(&reach.crossings[(k + 1) % n]);
        let around = (to - from).rem_euclid(TAU);
        let (s, e) = (seats[k][1], seats[(k + 1) % n][0]);
        let x = reach.radial(from);
        let touch = points[s.wall] - reach.base;
        let (height, inner) = (touch.dot(up), (touch - up * touch.dot(up)).magnitude());
        let foot = push(
            &mut curves,
            (s.foot, e.foot),
            arc(reach.base, up, x, rr, around, tol)?,
        );
        let top = push(
            &mut curves,
            (s.wall, e.wall),
            arc(reach.base + up * height, up, x, inner, around, tol)?,
        );
        let out = push(
            &mut curves,
            (s.line, e.line),
            arc(reach.base, up, x, big, around, tol)?,
        );
        let mid = 0.5f64.mul_add(around, from);
        faces.push(ToolFace {
            surface: wall.clone(),
            curves: vec![foot, s.up, top, e.up],
            at: reach.foot(mid) + up * (height * 0.5),
            out: -reach.radial(mid),
        });
        faces.push(ToolFace {
            surface: floor.clone(),
            curves: vec![foot, s.ray, out, e.ray],
            at: reach.base + reach.radial(mid) * (0.5 * (rr + big)),
            out: -up,
        });
        let (at, out_dir) = reach.blend_mid(mid);
        blends.push(faces.len());
        faces.push(ToolFace {
            surface: blend.clone(),
            curves: vec![top, s.section, out, e.section],
            at,
            out: out_dir,
        });
    }
    // Over each existing round.
    for (k, c) in reach.crossings.iter().enumerate() {
        let Some(band) = strip(reach, c, tol)? else {
            return Ok(None);
        };
        let (s, e) = (seats[k][0], seats[k][1]);
        let (a, b) = reach.span(c);
        let foot = push(
            &mut curves,
            (s.foot, e.foot),
            arc(reach.base, up, reach.radial(a), rr, b - a, tol)?,
        );
        let Some(wall_side) = boundary_from(&band, &points, (s.wall, e.wall), tol)? else {
            return Ok(None);
        };
        let wall_b = push(
            &mut curves,
            wall_side.ends,
            (wall_side.curve, wall_side.range),
        );
        let Some(round_side) = boundary_from(&band, &points, (s.line, e.line), tol)? else {
            return Ok(None);
        };
        let round_b = push(
            &mut curves,
            round_side.ends,
            (round_side.curve, round_side.range),
        );
        let line_b = push(
            &mut curves,
            (s.line, e.line),
            segment(points[s.line], points[e.line], tol)?,
        );
        let line_at = reach.base + c.across * c.gap;
        let (_, axis) = reach.round_axis(c, line_at);
        let on_round = {
            let t = reach.floor_touch(c.angle) - axis;
            let l = line_at - axis;
            let m = t / t.magnitude() + l / l.magnitude();
            m / m.magnitude()
        };
        let (band_at, band_out) = reach.blend_mid(c.angle);
        faces.push(ToolFace {
            surface: wall.clone(),
            curves: vec![foot, e.up, wall_b, s.up],
            at: reach.foot(c.angle) + up * (0.5 * (reach.wall_touch(c.angle) - reach.base).dot(up)),
            out: -c.across,
        });
        faces.push(ToolFace {
            surface: floor.clone(),
            curves: vec![foot, e.ray, line_b, s.ray],
            at: reach.base + c.across * (0.5 * (rr + c.gap)),
            out: -up,
        });
        faces.push(ToolFace {
            surface: curved.faces[c.round].surface.clone(),
            curves: vec![round_b, line_b],
            at: axis + on_round * c.round_radius,
            // Into the solid: toward a convex round's axis, away from a
            // concave one's.
            out: if c.concave { on_round } else { -on_round },
        });
        blends.push(faces.len());
        faces.push(ToolFace {
            surface: band.into(),
            curves: vec![wall_b, e.section, round_b, s.section],
            at: band_at,
            out: band_out,
        });
    }
    let Some((solid, made)) = closed(model, &points, &curves, &faces, tol)? else {
        return Ok(None);
    };
    Ok(Some((
        solid,
        blends.iter().map(|i| made[*i].clone()).collect(),
    )))
}

/// The foot blend of `edges` at `size` blended into the rounds it
/// reaches, fused into `solid` as one tool; `None` where `reach` does not
/// hold the request or the tool will not build, apply, or prove.
///
/// # Errors
///
/// As reading the solid and the booleans, other than their refusals.
pub(crate) fn blended(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    size: Size,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some(reach) = reach(model, solid, edges, size, tol)? else {
        return Ok(None);
    };
    let Some((fill, made)) = tool(model, &reach, tol)? else {
        return Ok(None);
    };
    let deflection = ogeom_mesh::Deflection::default();
    let before = ogeom_algo::volume_properties(model, solid, deflection, tol)?.mass;
    let added = ogeom_algo::volume_properties(model, &fill, deflection, tol)?.mass;
    let mut named = Vec::with_capacity(made.len() * edges.len());
    for e in edges {
        for face in &made {
            named.push((face.clone(), e.clone()));
        }
    }
    let tools = [super::fallback::Tool {
        solid: fill,
        sense: super::fallback::Sense::Fuse,
        faces: named,
    }];
    let Some(built) = super::fallback::applied(model, solid, edges, &tools, tol)? else {
        return Ok(None);
    };
    let built = super::merged(model, solid, built, tol)?;
    // The fill stands outside the solid, touching it only on the faces it
    // was closed on: the fuse adds its whole volume.
    let after = ogeom_algo::volume_properties(model, &built.shape, deflection, tol)?.mass;
    if (after - before - added).abs() > 1e-6 * (1.0 + before.abs()) {
        return Ok(None);
    }
    if !ogeom_algo::check(model, &built.shape, tol)?.is_valid() {
        return Ok(None);
    }
    Ok(Some(built))
}

/// Whether `edges` at `size` are a foot blend that blend into round
/// takes ([`reach`]).
///
/// # Errors
///
/// As reading the solid.
pub(crate) fn reach_found(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    size: Size,
    tol: Tolerances,
) -> OgeomResult<bool> {
    Ok(reach(model, solid, edges, size, tol)?.is_some())
}

/// The edges of `specs` and the one bevel they ask, where every spec asks
/// the same legs on a rim between a cylinder's wall and a plane square to
/// it; `None` for any other request.
///
/// # Errors
///
/// As reading the seats, other than a spec's own refusal, which the
/// engine's construction reports.
pub(crate) fn bevel_size(
    model: &Model,
    solid: &Shape,
    specs: &[(Shape, crate::Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<(Vec<Shape>, Size)>> {
    let eps = tol.confusion();
    let mut legs: Option<(f64, f64)> = None;
    for (edge, spec) in specs {
        let Ok(crate::chamfer::Legs::Rim(seat, on_wall, on_cap)) =
            crate::chamfer::legs(model, solid, edge, spec, tol)
        else {
            return Ok(None);
        };
        if seat.slope.abs() > tol.angular() {
            return Ok(None);
        }
        match legs {
            None => legs = Some((on_wall, on_cap)),
            Some((w, c)) if (w - on_wall).abs() <= eps && (c - on_cap).abs() <= eps => {}
            Some(_) => return Ok(None),
        }
    }
    let Some((rise, out)) = legs else {
        return Ok(None);
    };
    let edges = specs.iter().map(|(e, _)| e.clone()).collect();
    Ok(Some((edges, Size::Bevel { out, rise })))
}

/// Whether `built`, the rounds of `edges` on `solid`, changed a face of
/// `solid` other than the edges' own and those at their ends: a face the
/// rounds reached and trimmed through. A face is unchanged where the result
/// holds a face of its area and its bounds on whose surface a point of it
/// lies.
///
/// # Errors
///
/// As exploring the shapes and measuring the faces.
pub(crate) fn changes_others(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    built: &Built,
    tol: Tolerances,
) -> OgeomResult<bool> {
    use ogeom_topo::{ShapeType, ancestors_of, explore_unique};
    let eps = tol.confusion() * 10.0;
    let deflection = ogeom_mesh::Deflection::default();
    let mut own: Vec<Shape> = Vec::new();
    for edge in edges {
        own.extend(ancestors_of(model, solid, edge, ShapeType::Face)?);
        for vertex in explore_unique(model, edge, ShapeType::Vertex)? {
            own.extend(ancestors_of(model, solid, &vertex, ShapeType::Face)?);
        }
    }
    let read = |face: &Shape| -> OgeomResult<(f64, Aabb)> {
        Ok((
            ogeom_algo::surface_properties(model, face, deflection, tol)?.mass,
            ogeom_algo::shape_bounds(model, face, tol)?,
        ))
    };
    let same_box = |a: &Aabb, b: &Aabb| match (a.low(), a.high(), b.low(), b.high()) {
        (Some(al), Some(ah), Some(bl), Some(bh)) => al.distance(bl).max(ah.distance(bh)) <= eps,
        _ => false,
    };
    let made: Vec<(Shape, f64, Aabb)> = explore_unique(model, &built.shape, ShapeType::Face)?
        .into_iter()
        .map(|f| read(&f).map(|(area, bounds)| (f, area, bounds)))
        .collect::<OgeomResult<_>>()?;
    for face in explore_unique(model, solid, ShapeType::Face)? {
        if own.iter().any(|o| o.is_same(&face)) {
            continue;
        }
        let (area, bounds) = read(&face)?;
        let (p, _) = ogeom_algo::face_normal(model, &face, tol)?;
        let mut kept = false;
        for (g, now, around) in &made {
            if (area - now).abs() > 1e-9 * (1.0 + area) || !same_box(&bounds, around) {
                continue;
            }
            let Some(surface) = model
                .node(g)
                .and_then(|x| x.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
            else {
                continue;
            };
            let on = ogeom_geom::project::project_on_surface(surface, p, 16, tol)?.distance;
            if on <= eps {
                kept = true;
                break;
            }
        }
        if !kept {
            return Ok(true);
        }
    }
    Ok(false)
}
