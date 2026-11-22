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
//! A round's rim may also be an open arc of the foot circle whose two ends
//! stand on planes through the wall's axis, the rounds it reaches inside
//! the arc: each end's plane cuts the fill on the section there, which
//! lies in it, so the tool is capped flush in that plane (amendment "a rim
//! arc ending obliquely on a flat").
//!
//! On a cylinder's closed rim the arcs may meet: where convex rounds of
//! one radius run on into each other round a corner of the floor (a torus
//! about an axis square to the floor, or a sphere), the top past the floor
//! lies at that radius from the floor's region sunk by it, so the ball's
//! centre stands over the nearest point of the floor's outline and rests
//! on the round there, one fitted piece per round, passing from one to the
//! next over the edge between them.
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
/// The angle under which a bevel's leg on the floor grazes a tangent line
/// it touches or crosses: the floor and the round are left wedges that
/// sharp at the crossings, which the booleans do not resolve reliably.
const GRAZE: f64 = 1e-2;

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

/// An open arc of the foot: where it starts about the axis, how far it
/// turns, and the planes through the axis its two ends stand on.
#[derive(Debug, Clone, Copy)]
struct Open {
    /// The start's angle about the axis.
    from: f64,
    /// Counter-clockwise about the floor's normal, under a whole turn.
    sweep: f64,
    /// The end faces at `from` and at `from + sweep`.
    faces: [usize; 2],
}

/// One piece of the floor's outline, seen from above, along which a round
/// meets the floor tangentially.
#[derive(Debug, Clone, Copy)]
enum Outline {
    /// A straight edge from `from` to `to`; `out` lies in the floor, square
    /// to it, away from the floor.
    Line { from: Point, to: Point, out: Vector },
    /// An arc about `centre` on the floor, the floor inside it, turning
    /// `sweep` counter-clockwise about the floor's normal from the direction
    /// `start`.
    Arc {
        centre: Point,
        radius: f64,
        start: Vector,
        sweep: f64,
    },
    /// A corner of the floor rounded by a sphere about the point under it.
    Corner(Point),
}

/// A piece of the outline, the round's face past it, and the floor's
/// vertices at its ends in the loop's order.
#[derive(Debug, Clone, Copy)]
struct Element {
    outline: Outline,
    face: usize,
    /// The floor's edge along it; none at a corner.
    edge: Option<Id>,
    ends: (Id, Id),
}

impl Element {
    /// For a point `p` of the floor's plane: its foot on the outline, the
    /// direction there away from the floor, and how far `p` stands past it
    /// (negative on the floor's side). Where `strict`, `None` where the foot
    /// falls off this piece.
    fn seat(&self, p: Point, up: Vector, strict: bool) -> Option<(Point, Vector, f64)> {
        match self.outline {
            Outline::Line { from, to, out } => {
                let along = to - from;
                let t = (p - from).dot(along) / along.dot(along);
                if strict && !(0.0..=1.0).contains(&t) {
                    return None;
                }
                let q = from + along * t;
                Some((q, out, (p - q).dot(out)))
            }
            Outline::Arc {
                centre,
                radius,
                start,
                sweep,
            } => {
                let v = p - centre;
                let v = v - up * v.dot(up);
                let m = v.magnitude();
                if m <= f64::EPSILON * radius {
                    return None;
                }
                let w = v / m;
                let turn = up.dot(start.cross(w)).atan2(start.dot(w)).rem_euclid(TAU);
                if strict && turn > sweep {
                    return None;
                }
                Some((centre + w * radius, w, m - radius))
            }
            Outline::Corner(at) => {
                let v = p - at;
                let m = v.magnitude();
                (m > 0.0).then(|| (at, v / m, m))
            }
        }
    }
}

/// The piece of `chain` whose foot stands nearest `p`, with that foot; a
/// corner gives way to a piece beside it at the same distance.
fn nearest(chain: &[Element], p: Point, up: Vector) -> Option<(usize, (Point, Vector, f64))> {
    let mut best: Option<(usize, (Point, Vector, f64))> = None;
    for (k, e) in chain.iter().enumerate() {
        let Some(seat) = e.seat(p, up, true) else {
            continue;
        };
        let corner = matches!(e.outline, Outline::Corner(_));
        let better = best.is_none_or(|(_, (_, _, held))| {
            let (d, held) = (seat.2.abs(), held.abs());
            let slack = 1e-12 * (1.0 + held);
            if corner {
                d < held - slack
            } else {
                d <= held + slack
            }
        });
        if better {
            best = Some((k, seat));
        }
    }
    best
}

/// A stretch of the floor's outline the contact circle crosses, along
/// which rounds of one radius meet the floor tangentially and one another
/// G1: the ball leaves the floor at the first seat, rests on each piece's
/// round in turn, changing at the seats between, and comes back to the
/// floor at the last.
#[derive(Debug, Clone)]
struct Run {
    round_radius: f64,
    /// The pieces the ball rests on, in order about the axis.
    elements: Vec<Element>,
    /// About the axis, increasing, one more than the pieces.
    seats: Vec<f64>,
    /// At each seat between two pieces: the floor's vertex there, and the
    /// solid's edge between the two rounds' faces running down from it.
    junctions: Vec<(Id, Id)>,
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
    /// By angle, over disjoint arcs; along an open arc, in its order and
    /// inside it.
    crossings: Vec<Crossing>,
    /// Stretches of the floor's outline the contact circle crosses where
    /// rounds meet one another; where this holds any, `crossings` holds
    /// none.
    runs: Vec<Run>,
    /// The rim's ends, where it is an open arc.
    open: Option<Open>,
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

    /// `theta` turned to lie at or past an open arc's start, within a turn
    /// of it; as it is round a closed rim.
    fn unwrapped(&self, theta: f64) -> f64 {
        self.open
            .map_or(theta, |o| o.from + (theta - o.from).rem_euclid(TAU))
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

    /// The run whose stretch holds `theta` strictly, the piece there, and
    /// the foot on it of the ball centre's point on the floor.
    fn run_at(&self, theta: f64) -> Option<(&Run, (Point, Vector, f64))> {
        let p = self.on_line(theta);
        self.runs.iter().find_map(|run| {
            let from = run.seats[0];
            let t = from + (theta - from).rem_euclid(TAU);
            let last = run.seats[run.seats.len() - 1];
            if t <= from || t >= last {
                return None;
            }
            let k = run.seats[1..].iter().position(|s| t <= *s)?;
            Some((run, run.elements[k].seat(p, self.up, false)?))
        })
    }

    /// The ball centre's height over the floor where it stands `d` past a
    /// run's outline.
    fn run_height(&self, run: &Run, d: f64) -> f64 {
        let a = run.round_radius;
        (a + self.radius)
            .mul_add(a + self.radius, -d * d)
            .max(0.0)
            .sqrt()
            - a
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
        if let Some((run, (_, _, d))) = self.run_at(theta) {
            return (self.run_height(run, d), 0.0);
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
        if let Some((run, (q, _, _))) = self.run_at(theta) {
            let axis = q - self.up * run.round_radius;
            let w = centre - axis;
            return axis + w * (run.round_radius / w.magnitude());
        }
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
        if !self.runs.is_empty() {
            // The runs stand apart in order round the closed rim.
            let cuts: Vec<f64> = self
                .runs
                .iter()
                .flat_map(|r| r.seats.iter().copied())
                .collect();
            for k in 0..cuts.len() {
                let (a, b) = (cuts[k], cuts.get(k + 1).copied().unwrap_or(cuts[0] + TAU));
                for i in 0..PROOF {
                    #[expect(clippy::cast_precision_loss, reason = "a hundred stations")]
                    out.push((b - a).mul_add(i as f64 / PROOF as f64, a));
                }
            }
            return out;
        }
        if let Some(open) = self.open {
            let mut cuts = vec![open.from];
            for c in &self.crossings {
                let (a, b) = self.span(c);
                let a = self.unwrapped(a);
                cuts.extend([a, a + (b - a)]);
            }
            cuts.push(open.from + open.sweep);
            for w in cuts.windows(2) {
                for i in 0..PROOF {
                    #[expect(clippy::cast_precision_loss, reason = "a hundred stations")]
                    out.push((w[1] - w[0]).mul_add(i as f64 / PROOF as f64, w[0]));
                }
            }
            out.push(open.from + open.sweep);
            return out;
        }
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
pub(crate) fn reach(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    size: Size,
    tol: Tolerances,
) -> OgeomResult<Option<Reach>> {
    read_reach(model, solid, edges, size, &mut None, tol)
}

/// [`reach`], setting `graze` to the angle at which a bevel's leg on the
/// floor meets a tangent line it touches or crosses under [`GRAZE`], nought
/// where it touches, whether or not the reading takes the request.
#[expect(clippy::too_many_lines, reason = "one reading of the foot, in order")]
fn read_reach(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    size: Size,
    graze: &mut Option<f64>,
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
    // Every edge lies on the foot circle on the floor. Together they turn
    // once round the axis, or, for a round, they are an open arc whose
    // ends stand on planes through the axis.
    let mut turned = 0.0;
    let mut samples: Vec<Vector> = Vec::new();
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
            samples.push(q);
        }
    }
    let mut once: Vec<Id> = Vec::new();
    for e in &asked {
        let (a, b) = curved.edges[e].ends;
        for v in [a, b] {
            if let Some(k) = once.iter().position(|x| *x == v) {
                once.remove(k);
            } else {
                once.push(v);
            }
        }
    }
    let arc_ends = match once[..] {
        [] if (turned - TAU).abs() <= 1e-6 => None,
        [a, b] if round && turned < TAU - APART => {
            let ends = (curved.points[&a] - base, curved.points[&b] - base);
            let Some(found) = open_ends(
                &curved,
                (floor, wall),
                (a, b),
                ends,
                (&samples, turned),
                (axis, up),
                tol,
            ) else {
                return Ok(None);
            };
            Some(found)
        }
        _ => return Ok(None),
    };
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
    let foot = Foot {
        curved,
        asked,
        wall,
        floor,
        base,
        up,
        foot_radius,
        slope,
        arc_ends,
        reach_radius: foot_radius + offset,
    };
    if let Some(found) = crossing_reach(model, foot.clone(), size, graze, tol)? {
        return Ok(Some(found));
    }
    run_reach(model, foot, size, tol)
}

/// A boss's foot rim read off the solid, before the rounds its blend
/// reaches.
#[derive(Debug, Clone)]
struct Foot {
    curved: Curved,
    asked: Vec<Id>,
    wall: usize,
    floor: usize,
    base: Point,
    up: Vector,
    foot_radius: f64,
    slope: f64,
    /// An open arc's start, sweep and end faces.
    arc_ends: Option<(Vector, f64, [usize; 2])>,
    /// The contact circle's radius on the floor.
    reach_radius: f64,
}

/// [`read_reach`] where the contact circle crosses straight rounds over
/// arcs standing apart.
#[expect(clippy::too_many_lines, reason = "one reading of the rounds, in order")]
fn crossing_reach(
    model: &Model,
    foot: Foot,
    size: Size,
    graze: &mut Option<f64>,
    tol: Tolerances,
) -> OgeomResult<Option<Reach>> {
    let Foot {
        curved,
        asked,
        wall,
        floor,
        base,
        up,
        foot_radius,
        slope,
        arc_ends,
        reach_radius,
    } = foot;
    let eps = tol.confusion();
    let ang = tol.angular();
    let (radius, rise) = (size.out(), size.rise());
    let round = matches!(size, Size::Round(_));
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
        if gap <= foot_radius + eps * 10.0 {
            continue;
        }
        if !round && !concave && gap < reach_radius + eps * 10.0 {
            let angle = (gap / reach_radius).min(1.0).acos();
            if angle < GRAZE {
                *graze = Some(graze.map_or(angle, |g| g.min(angle)));
            }
        }
        if gap >= reach_radius - eps * 10.0 {
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
    let open = arc_ends.map(|(start, sweep, faces)| Open {
        from: start.dot(ey).atan2(start.dot(ex)),
        sweep,
        faces,
    });
    if let Some(o) = open {
        // Each crossing's arc stands inside the open arc, clear of its ends,
        // or wholly off it; one reaching an end would need the end's plane
        // to cut the fitted piece.
        let mut inside = Vec::with_capacity(crossings.len());
        for c in crossings {
            let half = (c.gap / reach_radius).acos();
            let mid = (c.angle - o.from).rem_euclid(TAU);
            let (lo, hi) = (mid - half, mid + half);
            if lo >= APART && hi <= o.sweep - APART {
                inside.push(c);
            } else if lo <= o.sweep + APART || hi >= TAU - APART {
                return Ok(None);
            }
        }
        inside.sort_by(|a, b| {
            let at = |c: &Crossing| (c.angle - o.from).rem_euclid(TAU);
            at(a).total_cmp(&at(b))
        });
        crossings = inside;
        if crossings.is_empty() {
            return Ok(None);
        }
    }
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
        runs: Vec::new(),
        open,
        size,
        radius,
    };
    // The crossings' arcs stand apart, so a torus lies between each two.
    let n = reach.crossings.len();
    if n > 1 {
        let pairs = if reach.open.is_some() { n - 1 } else { n };
        for k in 0..pairs {
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

/// The floor's edge `edge`, walked from `start` to `end`, read as a piece
/// of a round's outline: a line along a convex cylinder tangent to the
/// floor, or an arc along the outside of a torus about an axis square to
/// the floor, tangent to it; with the round's radius.
fn element(
    curved: &Curved,
    floor: usize,
    edge: Id,
    (start, end): (Id, Id),
    up: Vector,
    tol: Tolerances,
) -> OgeomResult<Option<(Element, f64)>> {
    let eps = tol.confusion() * 10.0;
    let ang = tol.angular();
    let others: Vec<usize> = faces_of(curved, edge)
        .into_iter()
        .filter(|f| *f != floor)
        .collect();
    let [face] = others[..] else {
        return Ok(None);
    };
    let e = &curved.edges[&edge];
    let (from, to) = (curved.points[&start], curved.points[&end]);
    // How far, square to the outline in the floor, the round's face reaches
    // past it: its edges sampled.
    let mut samples: Vec<Point> = Vec::new();
    for (x, _) in curved.faces[face].loops.iter().flatten() {
        let other = &curved.edges[x];
        for k in 0..=8 {
            let t = other.range.0 + (other.range.1 - other.range.0) * f64::from(k) / 8.0;
            samples.push(other.curve.point_at(t, tol)?);
        }
    }
    let flat = |v: Vector| v - up * v.dot(up);
    let found = match (&e.curve, &curved.faces[face].surface) {
        (Curve::Line(line), SurfaceGeometry::Cylinder(c)) => {
            let along = line.axis().direction.vector();
            let axis = c.cylinder().axis();
            let a = c.cylinder().radius();
            let off = from - axis.location;
            let off = off - along * off.dot(along);
            if axis.direction.vector().cross(along).magnitude() > ang
                || (off - up * a).magnitude() > eps
                || outward(curved, face, from).is_none_or(|n| n.dot(up) < 1.0 - ang)
            {
                return Ok(None);
            }
            let side = up.cross(along);
            let Some(far) = samples
                .iter()
                .map(|q| flat(*q - from).dot(side))
                .max_by(|x, y| x.abs().total_cmp(&y.abs()))
            else {
                return Ok(None);
            };
            if far.abs() <= eps {
                return Ok(None);
            }
            let out = side * far.signum();
            (Outline::Line { from, to, out }, a)
        }
        (Curve::Circle(circle), SurfaceGeometry::Torus(t)) => {
            let circle = circle.circle();
            let torus = t.torus();
            let (centre, rho) = (circle.centre(), circle.radius());
            let a = torus.minor_radius();
            let z = torus.frame().z().vector();
            let core = torus.frame().origin();
            if z.cross(up).magnitude() > ang
                || circle.frame().z().vector().cross(up).magnitude() > ang
                || core.distance(centre - up * a) > eps
                || (torus.major_radius() - rho).abs() > eps
                || rho <= eps
            {
                return Ok(None);
            }
            // The round runs down the outside of the torus: the floor lies
            // inside the arc.
            let reach = samples
                .iter()
                .map(|q| flat(*q - centre).magnitude() - rho)
                .fold(0.0_f64, |m, d| if d.abs() > m.abs() { d } else { m });
            if reach <= eps {
                return Ok(None);
            }
            let mid = e.curve.point_at(0.5 * (e.range.0 + e.range.1), tol)?;
            let dir = |p: Point| {
                let v = flat(p - centre);
                v / v.magnitude()
            };
            let ccw = |a: Vector, b: Vector| up.dot(a.cross(b)).atan2(a.dot(b)).rem_euclid(TAU);
            let (s, m, f) = (dir(from), dir(mid), dir(to));
            let (start, sweep) = if ccw(s, m) < ccw(s, f) {
                (s, ccw(s, f))
            } else {
                (f, ccw(f, s))
            };
            (
                Outline::Arc {
                    centre,
                    radius: rho,
                    start,
                    sweep,
                },
                a,
            )
        }
        _ => return Ok(None),
    };
    Ok(Some((
        Element {
            outline: found.0,
            face,
            edge: Some(edge),
            ends: (start, end),
        },
        found.1,
    )))
}

/// The direction away from the floor at the vertex `at` of the piece `e`.
fn out_at(e: &Element, at: Point, up: Vector) -> Vector {
    e.seat(at, up, false).map_or(Vector::ZERO, |(_, n, _)| n)
}

/// The floor's outline read as chains of round pieces of one radius, each
/// joined to the next G1 or across a convex corner a sphere of that radius
/// rounds; with the radius and whether the chain closes on itself.
fn chains(
    curved: &Curved,
    floor: usize,
    rim: &[Id],
    up: Vector,
    tol: Tolerances,
) -> OgeomResult<Vec<(Vec<Element>, f64, bool)>> {
    let eps = tol.confusion() * 10.0;
    let ang = tol.angular();
    let mut out = Vec::new();
    for ring in &curved.faces[floor].loops {
        if ring.iter().any(|(e, _)| rim.contains(e)) {
            continue;
        }
        // The loop's pieces in order; `None` where it breaks.
        let mut items: Vec<Option<(Element, f64)>> = Vec::new();
        let mut walked: Vec<(Id, Id)> = Vec::new();
        for (e, forward) in ring {
            let (a, b) = curved.edges[e].ends;
            walked.push(if *forward { (a, b) } else { (b, a) });
        }
        let n = walked.len();
        if (0..n).any(|k| walked[k].1 != walked[(k + 1) % n].0) {
            continue;
        }
        let pieces: Vec<Option<(Element, f64)>> = ring
            .iter()
            .zip(&walked)
            .map(|((e, _), ends)| element(curved, floor, *e, *ends, up, tol))
            .collect::<OgeomResult<_>>()?;
        for k in 0..n {
            items.push(pieces[k]);
            let (Some((x, a)), Some((y, b))) = (pieces[k], pieces[(k + 1) % n]) else {
                items.push(None);
                continue;
            };
            let v = walked[k].1;
            let at = curved.points[&v];
            if (a - b).abs() > eps {
                items.push(None);
                continue;
            }
            let (nx, ny) = (out_at(&x, at, up), out_at(&y, at, up));
            if nx.cross(ny).magnitude() <= ang && nx.dot(ny) > 0.0 {
                continue;
            }
            // A convex corner between two lines, rounded by a sphere about
            // the point under it.
            let convex = match (x.outline, y.outline) {
                (Outline::Line { from, .. }, Outline::Line { to, .. }) => {
                    (to - at).dot(nx) < -eps && (from - at).dot(ny) < -eps
                }
                _ => false,
            };
            let sphere = (0..curved.faces.len()).find(|f| {
                ![floor, x.face, y.face].contains(f)
                    && curved.faces[*f].loops.iter().flatten().any(|(e, _)| {
                        let ends = curved.edges[e].ends;
                        ends.0 == v || ends.1 == v
                    })
                    && matches!(&curved.faces[*f].surface, SurfaceGeometry::Sphere(s)
                        if s.sphere().centre().distance(at - up * a) <= eps
                            && (s.sphere().radius() - a).abs() <= eps)
            });
            match sphere {
                Some(face) if convex => items.push(Some((
                    Element {
                        outline: Outline::Corner(at),
                        face,
                        edge: None,
                        ends: (v, v),
                    },
                    a,
                ))),
                _ => items.push(None),
            }
        }
        // Split at the breaks, the chain before the first break joined to
        // the one after the last.
        let Some(first) = items.iter().position(Option::is_none) else {
            let a = items[0].map_or(0.0, |(_, a)| a);
            out.push((items.iter().flatten().map(|(e, _)| *e).collect(), a, true));
            continue;
        };
        let m = items.len();
        let mut chain: Vec<(Element, f64)> = Vec::new();
        for k in 1..=m {
            match items[(first + k) % m] {
                Some(x) => chain.push(x),
                None if !chain.is_empty() => {
                    let a = chain[0].1;
                    out.push((chain.iter().map(|(e, _)| *e).collect(), a, false));
                    chain.clear();
                }
                None => {}
            }
        }
    }
    Ok(out)
}

/// Where along a turn of the contact circle a ball over `chain` stands
/// past its outline: each stretch from the seat where it leaves the floor,
/// through the seats where its foot passes from one piece to the next, to
/// where it comes back, with the pieces in order. `None` where the
/// circle stands past it all the way round, or its feet do not pass the
/// pieces in turn.
fn stretches(
    chain: &[Element],
    closes: bool,
    at: impl Fn(f64) -> Point,
    up: Vector,
) -> Option<Vec<(Vec<f64>, Vec<usize>)>> {
    const SAMPLES: usize = 2048;
    let state = |theta: f64| -> Option<usize> {
        nearest(chain, at(theta), up).and_then(|(k, (_, _, d))| (d > 0.0).then_some(k))
    };
    // Between `a` in one state and `b` in another, where the state changes.
    let change = |mut a: f64, mut b: f64| {
        let held = state(a);
        for _ in 0..80 {
            let mid = 0.5 * (a + b);
            if state(mid) == held {
                a = mid;
            } else {
                b = mid;
            }
        }
        0.5 * (a + b)
    };
    #[expect(clippy::cast_precision_loss, reason = "two thousand samples")]
    let angle = |k: usize| TAU * k as f64 / SAMPLES as f64;
    let states: Vec<Option<usize>> = (0..SAMPLES).map(|k| state(angle(k))).collect();
    let start = states.iter().position(Option::is_none)?;
    let mut out: Vec<(Vec<f64>, Vec<usize>)> = Vec::new();
    let mut open: Option<(Vec<f64>, Vec<usize>)> = None;
    for i in 1..=SAMPLES {
        let (k0, k1) = (start + i - 1, start + i);
        let (was, now) = (states[k0 % SAMPLES], states[k1 % SAMPLES]);
        if was == now {
            continue;
        }
        let theta = change(angle(k0), angle(k1));
        // Leaving the floor or coming back, the circle crosses the outline
        // there, rather than running off the chain's end.
        if was.is_none() != now.is_none() {
            let past = if was.is_none() { 1e-9 } else { -1e-9 };
            let crosses =
                nearest(chain, at(theta + past), up).is_some_and(|(_, (_, _, d))| d.abs() <= 1e-6);
            if !crosses {
                return None;
            }
        }
        match (was, now, open.as_mut()) {
            (None, Some(k), None) => open = Some((vec![theta], vec![k])),
            (Some(_), Some(k), Some((seats, pieces))) => {
                seats.push(theta);
                pieces.push(k);
            }
            (Some(_), None, Some(_)) => {
                let (mut seats, pieces) = open.take()?;
                seats.push(theta);
                out.push((seats, pieces));
            }
            _ => return None,
        }
    }
    // The pieces follow one another along the chain, one way.
    let n = chain.len();
    for (_, pieces) in &out {
        let step = |a: usize, b: usize| {
            if b == (a + 1) % n && (closes || b == a + 1) {
                Some(1)
            } else if a == (b + 1) % n && (closes || a == b + 1) {
                Some(-1)
            } else {
                None
            }
        };
        let steps: Option<Vec<i32>> = pieces.windows(2).map(|w| step(w[0], w[1])).collect();
        let steps = steps?;
        if steps.windows(2).any(|w| w[0] != w[1]) {
            return None;
        }
    }
    Some(out)
}

/// [`read_reach`] where the contact circle crosses a stretch of the floor's
/// outline along which rounds of one radius meet one another, at a corner
/// a torus or a sphere rounds: a round of a closed rim on a cylinder's
/// foot, the rounds convex.
#[expect(clippy::too_many_lines, reason = "one reading of the runs, in order")]
fn run_reach(model: &Model, foot: Foot, size: Size, tol: Tolerances) -> OgeomResult<Option<Reach>> {
    let eps = tol.confusion();
    let Size::Round(radius) = size else {
        return Ok(None);
    };
    if foot.slope != 0.0 || foot.arc_ends.is_some() {
        return Ok(None);
    }
    let Foot {
        curved,
        asked,
        wall,
        floor,
        base,
        up,
        foot_radius,
        reach_radius,
        ..
    } = foot;
    let first = curved.edges[&asked[0]]
        .curve
        .point_at(curved.edges[&asked[0]].range.0, tol)?;
    let ex = first - base;
    let ex = (ex - up * ex.dot(up)) / (ex - up * ex.dot(up)).magnitude();
    let ey = up.cross(ex);
    let at = |theta: f64| base + (ex * theta.cos() + ey * theta.sin()) * reach_radius;
    let mut runs: Vec<Run> = Vec::new();
    let mut corner = false;
    for (chain, round_radius, closes) in chains(&curved, floor, &asked, up, tol)? {
        let Some(found) = stretches(&chain, closes, at, up) else {
            return Ok(None);
        };
        for (seats, pieces) in found {
            // Increasing along the stretch, from a start within the first
            // turn.
            let shift = seats[0] - seats[0].rem_euclid(TAU);
            let seats: Vec<f64> = seats.iter().map(|s| s - shift).collect();
            let elements: Vec<Element> = pieces.iter().map(|k| chain[*k]).collect();
            if seats.windows(2).any(|w| w[1] - w[0] < APART) {
                return Ok(None);
            }
            let mut junctions = Vec::with_capacity(elements.len() - 1);
            for w in elements.windows(2) {
                let v = if w[0].ends.1 == w[1].ends.0 {
                    w[0].ends.1
                } else {
                    w[0].ends.0
                };
                let shared = curved.edges.iter().find_map(|(id, e)| {
                    let mut f = faces_of(&curved, *id);
                    f.sort_unstable();
                    let mut pair = vec![w[0].face, w[1].face];
                    pair.sort_unstable();
                    (f == pair && (e.ends.0 == v || e.ends.1 == v)).then_some(*id)
                });
                let Some(shared) = shared else {
                    return Ok(None);
                };
                junctions.push((v, shared));
            }
            corner |= elements.len() > 1 || !matches!(elements[0].outline, Outline::Line { .. });
            runs.push(Run {
                round_radius,
                elements,
                seats,
                junctions,
            });
        }
    }
    // Only a stretch through a corner is read here; one along a single
    // round is the crossings' own.
    if !corner {
        return Ok(None);
    }
    runs.sort_by(|a, b| a.seats[0].total_cmp(&b.seats[0]));
    let n = runs.len();
    for k in 0..n {
        let end = runs[k].seats[runs[k].seats.len() - 1];
        let next = runs[(k + 1) % n].seats[0] + if k + 1 == n { TAU } else { 0.0 };
        if next - end < APART {
            return Ok(None);
        }
    }
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
        slope: 0.0,
        crossings: Vec::new(),
        runs,
        open: None,
        size,
        radius,
    };
    // Every ball stays clear of the bottom of the run's rounds, and its
    // contact at each seat between two of them lies on the edge they share.
    let mut highest = radius;
    for theta in reach.proof_angles() {
        if let Some((run, (_, _, d))) = reach.run_at(theta)
            && d >= (run.round_radius + radius) * (1.0 - 1e-6)
        {
            return Ok(None);
        }
        highest = highest.max((reach.wall_touch(theta) - base).dot(up));
    }
    for run in &reach.runs {
        for (k, (_, edge)) in run.junctions.iter().enumerate() {
            let touch = reach.floor_touch(run.seats[k + 1]);
            let measure =
                ogeom_algo::PointDistance::of(model, &reach.curved.edges[edge].shape, tol)?;
            if measure.nearest(touch, tol)?.0 > eps * 10.0 {
                return Ok(None);
            }
        }
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

/// The open arc of a foot rim between the vertices `vertices`, standing at
/// `ends` from the foot circle's centre: its start, its sweep
/// counter-clockwise about the floor's normal, and the end faces at its
/// start and at its end. Each end's one face besides the floor and the
/// wall must be a plane holding the wall's axis; `None` otherwise, or where
/// the rim's `samples`, turning `turned` in all, do not run one way round
/// from one end to the other.
fn open_ends(
    curved: &Curved,
    (floor, wall): (usize, usize),
    vertices: (Id, Id),
    ends: (Vector, Vector),
    (samples, turned): (&[Vector], f64),
    (axis, up): (ogeom_math::Axis, Vector),
    tol: Tolerances,
) -> Option<(Vector, f64, [usize; 2])> {
    let eps = tol.confusion() * 10.0;
    let end_face = |v: Id| -> Option<usize> {
        let touching: Vec<usize> = (0..curved.faces.len())
            .filter(|f| *f != floor && *f != wall)
            .filter(|f| {
                curved.faces[*f].loops.iter().flatten().any(|(e, _)| {
                    curved
                        .edges
                        .get(e)
                        .is_some_and(|x| x.ends.0 == v || x.ends.1 == v)
                })
            })
            .collect();
        let [f] = touching[..] else {
            return None;
        };
        let SurfaceGeometry::Plane(plane) = &curved.faces[f].surface else {
            return None;
        };
        let plane = plane.plane();
        let n = plane.normal().vector();
        let holds = (axis.location - plane.origin()).dot(n).abs() <= eps
            && n.dot(axis.direction.vector()).abs() <= tol.angular();
        holds.then_some(f)
    };
    let faces = [end_face(vertices.0)?, end_face(vertices.1)?];
    let ccw = |from: Vector, to: Vector| up.dot(from.cross(to)).atan2(from.dot(to)).rem_euclid(TAU);
    let runs = |start: Vector, end: Vector| {
        (ccw(start, end) - turned).abs() <= 1e-6
            && samples.iter().all(|q| {
                let a = ccw(start, *q);
                a <= turned + 1e-6 || a >= TAU - 1e-6
            })
    };
    if runs(ends.0, ends.1) {
        Some((ends.0, turned, faces))
    } else if runs(ends.1, ends.0) {
        Some((ends.1, turned, [faces[1], faces[0]]))
    } else {
        None
    }
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
        if let Some((_, (q, _, _))) = reach.run_at(theta) {
            points.push(q);
        }
        if let Some(arc) = reach.section(theta, 9, tol) {
            points.extend(arc);
        }
        let centre = matches!(reach.size, Size::Round(_)).then(|| reach.centre(theta));
        stations.push((points, centre));
    }
    let mut room: Vec<Aabb> = Vec::with_capacity(n);
    let mut balls: Vec<(Aabb, Point, f64)> = Vec::with_capacity(n);
    let pairs = if reach.open.is_some() { n - 1 } else { n };
    for i in 0..pairs {
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
    let ends: Vec<(usize, f64)> = reach.open.map_or_else(Vec::new, |o| {
        vec![(o.faces[0], o.from), (o.faces[1], o.from + o.sweep)]
    });
    let mut seated: Vec<usize> = [reach.wall, reach.floor]
        .into_iter()
        .chain(reach.crossings.iter().map(|x| x.round))
        .chain(
            reach
                .runs
                .iter()
                .flat_map(|r| r.elements.iter().map(|e| e.face)),
        )
        .collect();
    seated.sort_unstable();
    seated.dedup();
    for (face, _) in &ends {
        if !seated.contains(face) {
            seated.push(*face);
        }
    }
    // At an end, the cap's two straight sides: up the wall from the foot,
    // and out along the floor. An edge of a seated face on either line
    // bounds the cap's plane there, not the fill's room.
    let on_cap_side = |e: &Id| -> bool {
        let edge = &c.edges[e];
        if !matches!(edge.curve, Curve::Line(_)) {
            return false;
        }
        let (a, b) = (c.points[&edge.ends.0], c.points[&edge.ends.1]);
        for (_, theta) in &ends {
            let foot = reach.foot(*theta);
            for along in [reach.wall_touch(*theta) - foot, reach.radial(*theta)] {
                let along = along / along.magnitude();
                let off = |p: Point| {
                    let d = p - foot;
                    (d - along * d.dot(along)).magnitude()
                };
                if off(a) <= eps * 10.0 && off(b) <= eps * 10.0 {
                    return true;
                }
            }
        }
        false
    };
    let lines: Vec<Id> = reach.crossings.iter().map(|x| x.line).collect();
    // A run's outline and the edges between its rounds bound the fill.
    let along: Vec<Id> = reach
        .runs
        .iter()
        .flat_map(|r| {
            r.elements
                .iter()
                .filter_map(|e| e.edge)
                .chain(r.junctions.iter().map(|(_, e)| *e))
        })
        .collect();
    let corners: Vec<Id> = reach
        .runs
        .iter()
        .flat_map(|r| r.junctions.iter().map(|(v, _)| *v))
        .collect();
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
            if reach.edges.contains(e) || lines.contains(e) || along.contains(e) {
                continue;
            }
            // A seam of the wall is the wall's own.
            if faces_of(c, *e) == [reach.wall] || on_cap_side(e) {
                continue;
            }
            if meets(&ogeom_algo::shape_bounds(model, &c.edges[e].shape, tol)?) {
                return Ok(false);
            }
        }
    }
    for line in lines.iter().chain(&along) {
        let edge = &c.edges[line];
        for v in [edge.ends.0, edge.ends.1] {
            if !corners.contains(&v) && room.iter().any(|r| r.contains(c.points[&v])) {
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
    strip_over(reach, (a, b), tol)
}

/// The piece of the blend between the angles `a` and `b` about the axis,
/// over one round: its section at each station, fitted.
fn strip_over(
    reach: &Reach,
    (a, b): (f64, f64),
    tol: Tolerances,
) -> OgeomResult<Option<ogeom_geom::BSplineSurface>> {
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
    closed_within(model, points, curves, faces, &[], tol)
}

/// [`closed`], each curve widened to its width in `widths` (by the curves'
/// order, missing or nought for none), as [`super::curved::build_with_widths`]
/// widens them: curves read off the solid's own edges, which meet within
/// those edges' tolerances.
pub(super) fn closed_within(
    model: &mut Model,
    points: &[Point],
    curves: &[SheetCurve],
    faces: &[ToolFace],
    widths: &[f64],
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
    let Ok(built) = super::curved::build_with_widths(model, points, curves, &sheet, widths, tol)
    else {
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
    // The seats in order along the rim: an open arc's start, each
    // crossing's two seat changes, and an open arc's end.
    let mut angles: Vec<f64> = Vec::with_capacity(2 * n + 2);
    if let Some(o) = reach.open {
        angles.push(o.from);
    }
    for c in &reach.crossings {
        let (a, b) = reach.span(c);
        let a = reach.unwrapped(a);
        angles.extend([a, a + (b - a)]);
    }
    if let Some(o) = reach.open {
        angles.push(o.from + o.sweep);
    }
    // A run's seats, the ball between two of its rounds at the inner ones.
    let mut inner: Vec<bool> = vec![false; angles.len()];
    let mut run_seats: Vec<(usize, usize)> = Vec::with_capacity(reach.runs.len());
    for run in &reach.runs {
        let k = run.seats.len();
        run_seats.push((angles.len(), angles.len() + k - 1));
        angles.extend(&run.seats);
        inner.extend((0..k).map(|i| i > 0 && i + 1 < k));
    }
    let mut ring: Vec<Seat> = Vec::with_capacity(angles.len());
    for (&theta, &inner) in angles.iter().zip(&inner) {
        let at = points.len();
        points.extend([
            reach.foot(theta),
            reach.wall_touch(theta),
            if inner {
                reach.floor_touch(theta)
            } else {
                reach.on_line(theta)
            },
        ]);
        let (f, w, p) = (at, at + 1, at + 2);
        let radial = reach.radial(theta);
        let up_curve = push(&mut curves, (f, w), segment(points[f], points[w], tol)?);
        let ray = if inner {
            usize::MAX
        } else {
            push(&mut curves, (f, p), segment(points[f], points[p], tol)?)
        };
        let section = match reach.size {
            Size::Round(_) if inner => {
                // Between two rounds the contact leaves the plane through
                // the axis: the arc lies in the plane of the centre and
                // both touches.
                let centre = reach.centre(theta);
                let x = (points[w] - centre) / r;
                let y = (points[p] - centre) / r;
                let sweep = x.dot(y).clamp(-1.0, 1.0).acos();
                arc(centre, x.cross(y), x, r, sweep, tol)?
            }
            Size::Round(_) => {
                let centre = reach.centre(theta);
                let x = (points[w] - centre) / r;
                let sweep = x.dot(-up).clamp(-1.0, 1.0).acos();
                arc(centre, radial.cross(up), x, r, sweep, tol)?
            }
            Size::Bevel { .. } => segment(points[w], points[p], tol)?,
        };
        let section = push(&mut curves, (w, p), section);
        ring.push(Seat {
            foot: f,
            wall: w,
            line: p,
            up: up_curve,
            ray,
            section,
        });
    }
    // The torus pieces run from seat `2k` to `2k + 1` along an open arc,
    // and from seat `2k + 1` to the next round a closed rim; the fitted
    // pieces over the rounds between the others.
    let m = ring.len();
    let (over_floor_pieces, over_rounds): (Vec<_>, Vec<_>) = if !reach.runs.is_empty() {
        let k = run_seats.len();
        (
            (0..k)
                .map(|i| (run_seats[i].1, run_seats[(i + 1) % k].0))
                .collect(),
            Vec::new(),
        )
    } else if reach.open.is_some() {
        (
            (0..=n).map(|k| (2 * k, 2 * k + 1)).collect(),
            (0..n).map(|k| (2 * k + 1, 2 * k + 2)).collect(),
        )
    } else {
        (
            (0..n).map(|k| (2 * k + 1, (2 * k + 2) % m)).collect(),
            (0..n).map(|k| (2 * k, 2 * k + 1)).collect(),
        )
    };
    let curved = &reach.curved;
    let wall = curved.faces[reach.wall].surface.clone();
    let floor = curved.faces[reach.floor].surface.clone();
    let blend = over_floor(reach, tol)?;
    let mut faces: Vec<ToolFace> = Vec::new();
    let mut blends: Vec<usize> = Vec::new();
    // Round the floor, from each crossing to the next: the foot, the
    // wall's touch and the floor's touch.
    for &(i, j) in &over_floor_pieces {
        let (from, to) = (angles[i], angles[j]);
        let around = (to - from).rem_euclid(TAU);
        let (s, e) = (ring[i], ring[j]);
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
    for (c, &(i, j)) in reach.crossings.iter().zip(&over_rounds) {
        let Some(band) = strip(reach, c, tol)? else {
            return Ok(None);
        };
        let (s, e) = (ring[i], ring[j]);
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
    for (run, &(first, last)) in reach.runs.iter().zip(&run_seats) {
        if !run_pieces(
            reach,
            run,
            (first, last),
            &ring,
            (&mut points, &mut curves),
            (&mut faces, &mut blends),
            tol,
        )? {
            return Ok(None);
        }
    }
    // An open arc's caps, in its end faces' planes: the corner between the
    // wall and the floor less the ball's section there, facing away from
    // the arc.
    if let Some(o) = reach.open {
        for (k, (seat, theta, side)) in [
            (ring[0], angles[0], -1.0),
            (ring[m - 1], angles[m - 1], 1.0),
        ]
        .into_iter()
        .enumerate()
        {
            let (f, w, p) = (points[seat.foot], points[seat.wall], points[seat.line]);
            faces.push(ToolFace {
                surface: curved.faces[o.faces[k]].surface.clone(),
                curves: vec![seat.up, seat.section, seat.ray],
                at: f + (w - f) * 0.15 + (p - f) * 0.15,
                out: up.cross(reach.radial(theta)) * side,
            });
        }
    }
    let Some((solid, made)) = closed(model, &points, &curves, &faces, tol)? else {
        return Ok(None);
    };
    Ok(Some((
        solid,
        blends.iter().map(|i| made[*i].clone()).collect(),
    )))
}

/// The tool's faces over the run `run`, between the seats `first` and
/// `last` of `ring`: the wall and the floor under the whole stretch, and
/// over each of its rounds the fitted piece of the blend and the round's
/// face under it, bounded where the ball passes from one round to the next
/// by the edge between them. `false` where a piece will not fit.
fn run_pieces(
    reach: &Reach,
    run: &Run,
    (first, last): (usize, usize),
    ring: &[Seat],
    (points, curves): (&mut Vec<Point>, &mut Vec<SheetCurve>),
    (faces, blends): (&mut Vec<ToolFace>, &mut Vec<usize>),
    tol: Tolerances,
) -> OgeomResult<bool> {
    let up = reach.up;
    let a = run.round_radius;
    let curved = &reach.curved;
    let push = |curves: &mut Vec<SheetCurve>,
                ends: (usize, usize),
                (curve, range): (Curve, (f64, f64))| {
        curves.push(SheetCurve { ends, curve, range });
        curves.len() - 1
    };
    let (s, e) = (ring[first], ring[last]);
    let (from, to) = (run.seats[0], run.seats[run.seats.len() - 1]);
    let foot = push(
        curves,
        (s.foot, e.foot),
        arc(
            reach.base,
            up,
            reach.radial(from),
            reach.foot_radius,
            to - from,
            tol,
        )?,
    );
    // The floor's vertex at each inner seat: one point where both seats
    // either side of a sphere's corner stand on it.
    let mut corners: Vec<usize> = Vec::with_capacity(run.junctions.len());
    for (k, (v, _)) in run.junctions.iter().enumerate() {
        if k > 0 && run.junctions[k - 1].0 == *v {
            corners.push(corners[k - 1]);
        } else {
            points.push(curved.points[v]);
            corners.push(points.len() - 1);
        }
    }
    // Down the edge between two rounds, from the floor's vertex to the
    // ball's contact: an arc of the rounds' radius about the point under
    // the vertex, in the plane square to the outline there.
    let mut downs: Vec<usize> = Vec::with_capacity(corners.len());
    for (k, &corner) in corners.iter().enumerate() {
        let touch = ring[first + k + 1].line;
        let j = points[corner];
        let beside = if matches!(run.elements[k].outline, Outline::Corner(_)) {
            run.elements[k + 1]
        } else {
            run.elements[k]
        };
        let n = out_at(&beside, j, up);
        let centre = j - up * a;
        let sweep = ((points[touch] - centre) / a)
            .dot(up)
            .clamp(-1.0, 1.0)
            .acos();
        downs.push(push(
            curves,
            (corner, touch),
            arc(centre, up.cross(n), up, a, sweep, tol)?,
        ));
    }
    let count = run.elements.len();
    let mut wall_loop = vec![foot, e.up, s.up];
    let mut floor_loop = vec![foot, e.ray, s.ray];
    for (i, el) in run.elements.iter().enumerate() {
        let (sa, sb) = (ring[first + i], ring[first + i + 1]);
        let (ta, tb) = (run.seats[i], run.seats[i + 1]);
        let Some(band) = strip_over(reach, (ta, tb), tol)? else {
            return Ok(false);
        };
        let Some(wall_side) = boundary_from(&band, points, (sa.wall, sb.wall), tol)? else {
            return Ok(false);
        };
        let wall_b = push(curves, wall_side.ends, (wall_side.curve, wall_side.range));
        let Some(round_side) = boundary_from(&band, points, (sa.line, sb.line), tol)? else {
            return Ok(false);
        };
        let round_b = push(
            curves,
            round_side.ends,
            (round_side.curve, round_side.range),
        );
        wall_loop.push(wall_b);
        let pa = if i == 0 { sa.line } else { corners[i - 1] };
        let pb = if i + 1 == count { sb.line } else { corners[i] };
        let mut patch = vec![round_b];
        if i > 0 {
            patch.push(downs[i - 1]);
        }
        if i + 1 < count {
            patch.push(downs[i]);
        }
        let along = match el.outline {
            Outline::Corner(_) => None,
            Outline::Line { .. } => Some(segment(points[pa], points[pb], tol)?),
            Outline::Arc { centre, radius, .. } => {
                let flat = |p: Point| {
                    let v = p - centre;
                    let v = v - up * v.dot(up);
                    v / v.magnitude()
                };
                let (u, w) = (flat(points[pa]), flat(points[pb]));
                let turn = up.dot(u.cross(w)).atan2(u.dot(w));
                let z = if turn >= 0.0 { up } else { -up };
                Some(arc(centre, z, u, radius, turn.abs(), tol)?)
            }
        };
        if let Some(along) = along {
            let c = push(curves, (pa, pb), along);
            patch.push(c);
            floor_loop.push(c);
        }
        let mid = 0.5 * (ta + tb);
        let Some((q, _, _)) = el.seat(reach.on_line(mid), up, false) else {
            return Ok(false);
        };
        let axis = q - up * a;
        let m = (reach.floor_touch(mid) - axis) / a + up;
        let m = m / m.magnitude();
        faces.push(ToolFace {
            surface: curved.faces[el.face].surface.clone(),
            curves: patch,
            at: axis + m * a,
            // Into the solid, toward the round's core.
            out: -m,
        });
        let (band_at, band_out) = reach.blend_mid(mid);
        blends.push(faces.len());
        faces.push(ToolFace {
            surface: band.into(),
            curves: vec![wall_b, sb.section, round_b, sa.section],
            at: band_at,
            out: band_out,
        });
    }
    let mid = 0.5 * (from + to);
    let radial = reach.radial(mid);
    // Where the outline crosses the radial at `mid`, between the foot on
    // the floor and the contact circle past the outline.
    let (mut inside, mut outside) = (reach.foot_radius, reach.reach_radius());
    for _ in 0..60 {
        let t = 0.5 * (inside + outside);
        let past = nearest(&run.elements, reach.base + radial * t, up)
            .is_some_and(|(_, (_, _, d))| d > 0.0);
        if past {
            outside = t;
        } else {
            inside = t;
        }
    }
    faces.push(ToolFace {
        surface: curved.faces[reach.wall].surface.clone(),
        curves: wall_loop,
        at: reach.foot(mid) + up * (0.5 * (reach.wall_touch(mid) - reach.base).dot(up)),
        out: -radial,
    });
    faces.push(ToolFace {
        surface: curved.faces[reach.floor].surface.clone(),
        curves: floor_loop,
        at: reach.base + radial * (0.5 * (reach.foot_radius + inside)),
        out: -up,
    });
    Ok(true)
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

/// The angle at which the bevel `size` of the foot rim `edges` grazes the
/// tangent line of a straight convex round its leg on the floor touches or
/// crosses ([`GRAZE`]); `None` where the foot is not one blend into round
/// reads, or its leg meets no tangent line that sharply.
///
/// # Errors
///
/// As reading the solid.
pub(crate) fn grazes(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    size: Size,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let mut graze = None;
    read_reach(model, solid, edges, size, &mut graze, tol)?;
    Ok(graze)
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
