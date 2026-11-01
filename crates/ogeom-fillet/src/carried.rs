//! A concave edge's round whose ball changes what it rests on along the
//! edge.
//!
//! [`crate::rolled`] rounds a concave seam `E` whose ball touches one face
//! and rests on the far edge of the narrow face beside it all along. Where
//! the narrow face widens past the setback, narrows the other way, or ends
//! in a corner, the ball's contact on a side of `E` changes kind partway:
//! it holds a face, an edge of that face, or a corner of it. The round is
//! then a chain of pieces, each the envelope of a ball holding one thing on
//! each side, its centre on the curve those two holds leave it:
//!
//! - a face on each side: the line a radius off both planes, the ordinary
//!   cylinder between them;
//! - a face and an edge of the other face: the ellipse where the plane a
//!   radius off the face meets the cylinder a radius about the edge, the
//!   rolled round of [`crate::rolled`];
//! - a face and a corner of the other face: the circle where that plane
//!   meets the sphere a radius about the corner, a spindle torus's patch
//!   whose arcs all end on the corner.
//!
//! The chain changes piece where the ball touches three things at once, at
//! a point read in closed form: where the line meets the plane a radius
//! over a face's edge (the contact on that face reaching the edge), where
//! the ellipse's angle reaches the face's normal (the ball coming down on
//! the face), where the contact on an edge reaches its end, and where the
//! centre on the circle comes square to the next edge out of the corner.
//! At each such point the two curves are tangent: the plane, cylinder or
//! sphere the centre leaves is tangent there to the one it joins, so both
//! curves, cut from the same plane, touch. The sections of both pieces
//! there are one plane, holding one arc, and the pieces meet on it with a
//! common centre and tangent.
//!
//! The run is found from the one ball whose section crosses the crease's
//! middle with its contacts on their faces, edges and corners, in free
//! space, an edge it rests on being the far edge met crossing its face
//! square to the crease there; from that ball the run is walked to each
//! end, piece by piece, each ending at its first switch the way its
//! section moves. Every piece is then proved as the rolled run is: its
//! sections do not cross, its contacts stay on their faces and edges, and
//! every ball lies in free space.
//!
//! Where the ball comes onto the face beyond an edge it rests on, it no
//! longer holds anything of the edge's own faces on that side: it touches
//! the face it held and the face beyond, and rounds the corner between
//! them, which is the round of the seam between them, not of this edge.
//! That is refused by name, as is a ball holding an edge or a corner on
//! both sides at once.

use crate::consumed::{seat_of, segment};
use crate::march::{BlendStop, MarchedBlend, Sides};
use crate::marched::{ExactBand, Rail};
use crate::rolled::{Flat, basis, enters, face_distance, flat, host, inside};
use crate::support::{Seat, edge_curve, same_occurrence};
use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_geom::Curve3d as _;
use ogeom_geom::SurfaceGeometry;
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

/// Why a ball that comes onto the face beyond an edge it rests on is
/// refused: it then rounds the corner between that face and the one it
/// touches, which is another seam's round.
pub(crate) const BEYOND: &str = "the ball comes onto the face beyond that edge along the \
                                 crease, rounding the corner between that face and the one \
                                 it touches rather than this edge";

/// How many equal steps of each piece's parameter its band is fitted
/// through.
const STATIONS: usize = 64;

/// How many samples a piece's parameter is read at for its sections'
/// crossings of the crease, one way along it, and for a start.
const SEARCH: usize = 4096;

/// How many times an interval of a piece is halved, at most, to show that
/// the ball stays clear of a face between two stations.
const HALVINGS: u32 = 24;

/// How many chords, at most, a curved line of contact is read with.
const MAX_CHORDS: f64 = 1e5;

/// What the ball holds on one side of the crease.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Hold {
    /// The side's face, inside it.
    Face,
    /// An edge of the side's face, by its index among the side's edges.
    Edge(usize),
    /// A corner of the side's face.
    Corner(Point),
}

/// The curve the ball's centre runs on over a piece.
#[derive(Debug, Clone, Copy)]
enum Track {
    /// `base + along s`, `along` a unit vector.
    Line { base: Point, along: Vector },
    /// A radius off a plane and a radius off a line: the centre at angle
    /// `s` round the line is `base + along l(s) + radius (cos s u + sin s
    /// v)`, with `l(s) = l0 - radius (cos s au + sin s av) / k`.
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
    /// `centre + rho (cos s e1 + sin s e2)`.
    Circle {
        centre: Point,
        e1: Vector,
        e2: Vector,
        rho: f64,
    },
}

impl Track {
    fn centre(&self, s: f64) -> Point {
        match *self {
            Self::Line { base, along } => base + along * s,
            Self::Ellipse {
                base,
                along,
                u,
                v,
                radius,
                l0,
                au,
                av,
                k,
            } => {
                let (sn, cs) = s.sin_cos();
                base + along * (l0 - radius * (cs * au + sn * av) / k) + (u * cs + v * sn) * radius
            }
            Self::Circle {
                centre,
                e1,
                e2,
                rho,
            } => {
                let (sn, cs) = s.sin_cos();
                centre + (e1 * cs + e2 * sn) * rho
            }
        }
    }

    fn d1(&self, s: f64) -> Vector {
        let (sn, cs) = s.sin_cos();
        match *self {
            Self::Line { along, .. } => along,
            Self::Ellipse {
                along,
                u,
                v,
                radius,
                au,
                av,
                k,
                ..
            } => along * (radius * (sn * au - cs * av) / k) + (v * cs - u * sn) * radius,
            Self::Circle { e1, e2, rho, .. } => (e2 * cs - e1 * sn) * rho,
        }
    }

    fn d2(&self, s: f64) -> Vector {
        let (sn, cs) = s.sin_cos();
        match *self {
            Self::Line { .. } => Vector::ZERO,
            Self::Ellipse {
                along,
                u,
                v,
                radius,
                au,
                av,
                k,
                ..
            } => along * (radius * (cs * au + sn * av) / k) - (u * cs + v * sn) * radius,
            Self::Circle { e1, e2, rho, .. } => -(e1 * cs + e2 * sn) * rho,
        }
    }

    /// A bound on the centre's speed and on its acceleration in the
    /// parameter, over the whole track.
    fn speed(&self) -> f64 {
        match *self {
            Self::Line { .. } => 1.0,
            Self::Ellipse {
                radius, au, av, k, ..
            } => radius * (1.0 + au.hypot(av) / k.abs()),
            Self::Circle { rho, .. } => rho,
        }
    }

    /// The parameter at which the track passes through `c`, which it holds.
    fn at(&self, c: Point) -> f64 {
        match *self {
            Self::Line { base, along } => (c - base).dot(along),
            Self::Ellipse {
                base, along, u, v, ..
            } => {
                let w = c - (base + along * (c - base).dot(along));
                w.dot(v).atan2(w.dot(u))
            }
            Self::Circle { centre, e1, e2, .. } => {
                let w = c - centre;
                w.dot(e2).atan2(w.dot(e1))
            }
        }
    }
}

/// One side of the crease: its face, and the edges of that face other
/// than the crease, each with the outward normal of the face beyond it
/// where it is a convex planar seam.
struct Side {
    face: Shape,
    normal: Vector,
    /// Into the face from the crease, in its plane.
    into: Vector,
    edges: Vec<SideEdge>,
}

struct SideEdge {
    shape: Shape,
    ends: [Point; 2],
    beyond: Option<Vector>,
}

/// One piece of the chain: what the ball holds on each side, the track
/// its centre runs on, and the span of the track's parameter, from the
/// crease's start side to its end side.
#[derive(Debug, Clone, Copy)]
struct Piece {
    holds: [Hold; 2],
    track: Track,
    span: (f64, f64),
}

/// The run of a ball carried across changes of what it holds, from the
/// crease's start to its end, with what its bands are built from.
pub(crate) struct Carried {
    pieces: Vec<Piece>,
    sides: [Side; 2],
    guide: ogeom_geom::Curve,
    range: (f64, f64),
    e0: Point,
    t: Vector,
    length: f64,
    hosts: [(SurfaceGeometry, f64); 2],
}

/// The pieces' crease geometry and what they read, shared by the walk and
/// its proofs.
struct Crease<'a> {
    sides: &'a [Side; 2],
    e0: Point,
    t: Vector,
    length: f64,
    radius: f64,
    tol: Tolerances,
}

/// What follows a switch: the holds of the next piece, or why the ball
/// cannot go on.
enum Next {
    Holds([Hold; 2]),
    Refuse(&'static str),
}

/// The roots of `a cos s + b sin s = c`, unwrapped to stand within a turn
/// either side of `near`.
fn trig_roots(a: f64, b: f64, c: f64, near: f64) -> Vec<f64> {
    let m = a.hypot(b);
    if m == 0.0 || c.abs() > m {
        return Vec::new();
    }
    let phase = b.atan2(a);
    let half = (c / m).clamp(-1.0, 1.0).acos();
    let tau = core::f64::consts::TAU;
    let mut out = Vec::new();
    for root in [phase - half, phase + half] {
        let base = root - tau * ((root - near) / tau).round();
        for k in [-1.0, 0.0, 1.0] {
            out.push(base + tau * k);
        }
    }
    out
}

/// The greatest of `a cos s + b sin s` over `[lo, hi]`.
fn trig_max(a: f64, b: f64, lo: f64, hi: f64) -> f64 {
    let (lo, hi) = (lo.min(hi), lo.max(hi));
    let m = a.hypot(b);
    let phase = b.atan2(a);
    let tau = core::f64::consts::TAU;
    let peak = phase + tau * ((lo - phase) / tau).ceil();
    if peak <= hi {
        return m;
    }
    let f = |s: f64| a * s.cos() + b * s.sin();
    f(lo).max(f(hi))
}

impl Crease<'_> {
    fn conf(&self) -> f64 {
        self.tol.confusion()
    }

    /// The ball's contact on side `i` holding `hold` with its centre at `c`.
    fn contact(&self, i: usize, hold: Hold, c: Point) -> Point {
        match hold {
            Hold::Face => c - self.sides[i].normal * self.radius,
            Hold::Edge(j) => {
                let [p, q] = self.sides[i].edges[j].ends;
                let d = (q - p) / p.distance(q);
                p + d * (c - p).dot(d)
            }
            Hold::Corner(w) => w,
        }
    }

    /// Where along the crease the plane square to the track at `s` crosses
    /// it, and how fast that moves with `s`.
    fn crossing(&self, track: &Track, s: f64) -> Option<(f64, f64)> {
        let (c, d, dd) = (track.centre(s), track.d1(s), track.d2(s));
        let den = d.dot(self.t);
        if den == 0.0 {
            return None;
        }
        let num = d.dot(c - self.e0);
        let rate = ((dd.dot(c - self.e0) + d.dot(d)) * den - num * dd.dot(self.t)) / (den * den);
        Some((num / den, rate))
    }

    /// The side holding something other than its face, where one does.
    fn odd_side(holds: [Hold; 2]) -> Option<usize> {
        match holds {
            [Hold::Face, Hold::Face] => None,
            [Hold::Face, _] => Some(1),
            _ => Some(0),
        }
    }

    /// The frame round edge `j` of side `i` the ball rests on: the edge's
    /// start and unit direction, its length, `u` side `i`'s normal, `v`
    /// square to it toward the beyond face's normal, which stands at
    /// `phi3` round from `u`.
    fn edge_frame(
        &self,
        i: usize,
        j: usize,
    ) -> Result<(Point, Vector, f64, Vector, Vector, f64), &'static str> {
        let edge = &self.sides[i].edges[j];
        let Some(n3) = edge.beyond else {
            return Err("the ball comes to an edge of that face that cannot hold it");
        };
        let [p, q] = edge.ends;
        let reach = p.distance(q);
        let along = (q - p) / reach;
        let u = self.sides[i].normal;
        let w = n3 - u * n3.dot(u);
        let m = w.magnitude();
        if m <= self.tol.angular() {
            return Err("the ball comes to an edge of that face that cannot hold it");
        }
        let v = w / m;
        Ok((p, along, reach, u, v, n3.dot(v).atan2(n3.dot(u))))
    }

    /// The track the centre runs on holding `holds`.
    fn track(&self, holds: [Hold; 2]) -> Result<Track, &'static str> {
        let r = self.radius;
        match holds {
            [Hold::Face, Hold::Face] => {
                // A radius off both planes, square to the crease.
                let (n0, n1) = (self.sides[0].normal, self.sides[1].normal);
                let w = (n0 + n1) * (r / (1.0 + n0.dot(n1)));
                Ok(Track::Line {
                    base: self.e0 + w,
                    along: self.t,
                })
            }
            [Hold::Face, Hold::Edge(j)] | [Hold::Edge(j), Hold::Face] => {
                let i = usize::from(holds[1] != Hold::Face);
                let n1 = self.sides[1 - i].normal;
                let (base, along, _, u, v, _) = self.edge_frame(i, j)?;
                let k = n1.dot(along);
                if k.abs() <= self.tol.angular() {
                    return Err(
                        "the ball rests on an edge running parallel to the face it touches",
                    );
                }
                Ok(Track::Ellipse {
                    base,
                    along,
                    u,
                    v,
                    radius: r,
                    l0: (r - n1.dot(base - self.e0)) / k,
                    au: n1.dot(u),
                    av: n1.dot(v),
                    k,
                })
            }
            [Hold::Face, Hold::Corner(w)] | [Hold::Corner(w), Hold::Face] => {
                let i = usize::from(holds[1] != Hold::Face);
                let n1 = self.sides[1 - i].normal;
                // A radius over the touched face and a radius from the
                // corner, which stands `h` over that face.
                let h = (w - self.e0).dot(n1);
                if h <= 0.0 || h >= 2.0 * r {
                    return Err("the ball cannot touch its face and that corner at once");
                }
                let (e1, e2) = basis(n1);
                Ok(Track::Circle {
                    centre: w + n1 * (r - h),
                    e1,
                    e2,
                    rho: (h * (2.0 * r - h)).sqrt(),
                })
            }
            _ => Err("the ball rests on the far edges of both faces at once"),
        }
    }

    /// Where the contact on side `o`'s face, `c - radius n_o`, crosses an
    /// edge of that face ahead of `from`, given that contact's projection
    /// on the edge's in-plane normal `m` as `a cos s + b sin s + offset`
    /// from the edge's line: what holds the ball past it, an edge of each
    /// face or a face of neither, either refused.
    fn touched_exits(
        &self,
        o: usize,
        track: &Track,
        project: &dyn Fn(Vector, Point) -> (f64, f64, f64),
        ahead: &dyn Fn(f64) -> bool,
        from: f64,
        found: &mut Vec<(f64, Next)>,
    ) {
        let side = &self.sides[o];
        for edge in &side.edges {
            let [p, q] = edge.ends;
            let m = side.normal.cross(q - p);
            let m = m / m.magnitude();
            let (a, b, offset) = project(m, p);
            for s in trig_roots(a, b, -offset, from) {
                if !ahead(s) {
                    continue;
                }
                let at = self.contact(o, Hold::Face, track.centre(s));
                let w = (at - p).dot(q - p) / (q - p).dot(q - p);
                if (0.0..=1.0).contains(&w) {
                    found.push((
                        s,
                        Next::Refuse(if edge.beyond.is_some() {
                            "the ball rests on the far edges of both faces at once"
                        } else {
                            "the ball's contact runs off the face it touches"
                        }),
                    ));
                }
            }
        }
    }

    /// Every switch of `piece`'s track ahead of `from` in direction `sense`
    /// of its parameter, with what follows it.
    #[allow(clippy::too_many_lines, reason = "one case per kind of piece")]
    fn switches(
        &self,
        holds: [Hold; 2],
        track: &Track,
        from: f64,
        sense: f64,
    ) -> Result<Vec<(f64, Next)>, &'static str> {
        let conf = self.conf();
        let ahead = |s: f64| (s - from) * sense > conf / track.speed();
        let mut found: Vec<(f64, Next)> = Vec::new();
        match (Self::odd_side(holds), *track) {
            (None, Track::Line { base, .. }) => {
                // Each contact runs along the crease at the setback, and
                // leaves its face where it crosses one of the face's edges.
                for (i, side) in self.sides.iter().enumerate() {
                    let start = base - side.normal * self.radius;
                    let across = (start - self.e0).dot(side.into);
                    for (j, edge) in side.edges.iter().enumerate() {
                        let [p, q] = edge.ends;
                        let (py, qy) = ((p - self.e0).dot(side.into), (q - self.e0).dot(side.into));
                        if py == qy {
                            continue;
                        }
                        let w = (across - py) / (qy - py);
                        if !(0.0..=1.0).contains(&w) {
                            continue;
                        }
                        let at = (p + (q - p) * w - start).dot(self.t);
                        if !ahead(at) {
                            continue;
                        }
                        let off_end = p.distance(q) * w.min(1.0 - w);
                        let next = if off_end <= conf {
                            Next::Refuse("the ball's contact runs onto a corner of a face")
                        } else {
                            let mut next = [Hold::Face; 2];
                            next[i] = Hold::Edge(j);
                            Next::Holds(next)
                        };
                        found.push((at, next));
                    }
                }
            }
            (
                Some(i),
                Track::Ellipse {
                    radius,
                    l0,
                    au,
                    av,
                    k,
                    ..
                },
            ) => {
                let Hold::Edge(j) = holds[i] else {
                    return Err("the ball's holds do not match its track");
                };
                let (p, _, reach, _, _, phi3) = self.edge_frame(i, j)?;
                let q = self.sides[i].edges[j].ends[1];
                if ahead(0.0) {
                    let mut next = holds;
                    next[i] = Hold::Face;
                    found.push((0.0, Next::Holds(next)));
                }
                if ahead(phi3) {
                    found.push((phi3, Next::Refuse(BEYOND)));
                }
                // The contact on the touched face reaches one of its edges.
                if let Track::Ellipse {
                    base, along, u, v, ..
                } = *track
                {
                    let o = 1 - i;
                    let n = self.sides[o].normal;
                    let project = |m: Vector, p: Point| {
                        let am = along.dot(m);
                        (
                            radius * (u.dot(m) - au * am / k),
                            radius * (v.dot(m) - av * am / k),
                            (base - n * radius - p).dot(m) + l0 * am,
                        )
                    };
                    self.touched_exits(o, track, &project, &ahead, from, &mut found);
                }
                // The foot reaches an end of the edge: `l(s) = 0` or the
                // edge's length.
                for (at_end, corner) in [(0.0, p), (reach, q)] {
                    for s in trig_roots(au, av, (l0 - at_end) * k / radius, from) {
                        if ahead(s) && (0.0..=phi3).contains(&s) {
                            let mut next = holds;
                            next[i] = Hold::Corner(corner);
                            found.push((s, Next::Holds(next)));
                        }
                    }
                }
            }
            (
                Some(i),
                Track::Circle {
                    centre,
                    e1,
                    e2,
                    rho,
                },
            ) => {
                let Hold::Corner(w) = holds[i] else {
                    return Err("the ball's holds do not match its track");
                };
                // The contact on the touched face reaches one of its edges.
                let o = 1 - i;
                let n = self.sides[o].normal;
                let project = |m: Vector, p: Point| {
                    (
                        rho * e1.dot(m),
                        rho * e2.dot(m),
                        (centre - n * self.radius - p).dot(m),
                    )
                };
                self.touched_exits(o, track, &project, &ahead, from, &mut found);
                // Out of the corner, the ball comes square to an edge of the
                // face leaving it, and goes on resting on that edge.
                for (j, edge) in self.sides[i].edges.iter().enumerate() {
                    let [p, q] = edge.ends;
                    let d = if p.distance(w) <= conf {
                        (q - p) / p.distance(q)
                    } else if q.distance(w) <= conf {
                        (p - q) / p.distance(q)
                    } else {
                        continue;
                    };
                    let (a, b, c) = (rho * e1.dot(d), rho * e2.dot(d), -(centre - w).dot(d));
                    for s in trig_roots(a, b, c, from) {
                        // Entering the edge: the foot moving onto it.
                        let rate = rho * (e2.dot(d) * s.cos() - e1.dot(d) * s.sin());
                        if ahead(s) && rate * sense > 0.0 {
                            let mut next = holds;
                            next[i] = Hold::Edge(j);
                            found.push((s, Next::Holds(next)));
                        }
                    }
                }
            }
            _ => return Err("the ball's holds do not match its track"),
        }
        Ok(found)
    }
}

/// The run of a ball of `radius` rounding the concave straight `edge`
/// while it changes what it holds on a side, where it rests on the far
/// edge of a narrow face beside it somewhere along it; `None` where the
/// edge is not of that kind, no one ball holds the crease's middle, or one
/// piece holds the whole crease (the rolled run's own, refused there).
///
/// # Errors
///
/// [`OgeomError::Construction`] naming why, where the run starts but no
/// sound chain carries it to the crease's end: the ball comes onto the
/// face beyond an edge it rests on, holds an edge or a corner on both
/// sides at once, comes within its radius of another face, its contacts
/// leave their faces or edges, or its sections cross; and whatever reading
/// the solid reports.
#[allow(
    clippy::too_many_lines,
    reason = "one run, walked and proved end to end"
)]
pub(crate) fn carried_run(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    seat: &Seat,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Carried>> {
    let conf = tol.confusion();
    let mut sides = Vec::with_capacity(2);
    for i in 0..2 {
        let face = seat.faces[i].clone();
        if flat(model, &face, tol)?.is_none() {
            return Ok(None);
        }
        let mut edges = Vec::new();
        for e in explore_unique(model, &face, ShapeType::Edge)? {
            if same_occurrence(model, &e, edge, tol) {
                continue;
            }
            let Some(ends) = segment(model, &e, tol)? else {
                return Ok(None);
            };
            if ends[0].distance(ends[1]) <= conf {
                continue;
            }
            let beyond = match seat_of(model, solid, &e, tol)? {
                Some(s) if s.convex && s.faces.iter().any(|f| f.is_same(&face)) => s
                    .faces
                    .iter()
                    .position(|f| !f.is_same(&face))
                    .map(|k| s.normals[k]),
                _ => None,
            };
            edges.push(SideEdge {
                shape: e,
                ends,
                beyond,
            });
        }
        sides.push(Side {
            face,
            normal: seat.normals[i],
            into: -seat.leg(i, tol)?,
            edges,
        });
    }
    let sides: [Side; 2] = [sides.remove(0), sides.remove(0)];
    let (curve, range) = edge_curve(model, edge, tol)?;
    let e0 = curve.point_at(range.0, tol)?;
    let e1 = curve.point_at(range.1, tol)?;
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
    let crease = Crease {
        sides: &sides,
        e0,
        t,
        length,
        radius,
        tol,
    };

    // The ball holding the crease's middle: each kind of piece whose
    // section crosses the crease there with its contacts on their faces,
    // edges and corners, in free space. A run no one ball holds there is
    // not of this kind.
    let mut reach = Vec::new();
    for face in explore_unique(model, solid, ShapeType::Face)? {
        reach.push(match flat(model, &face, tol)? {
            Some(flat) => Reach::Flat(flat),
            None => Reach::Bounds(ogeom_algo::shape_bounds(model, &face, tol)?),
        });
    }
    let Some((holds, from)) = start(&crease, &reach, length * 0.5) else {
        return Ok(None);
    };

    // Walk the crease from there to each end, piece by piece, and join the
    // two halves of the piece they start on.
    let back = walk(&crease, holds, from, 0.0).map_err(refuse)?;
    let ahead = walk(&crease, holds, from, length).map_err(refuse)?;
    let mut pieces: Vec<Piece> = back
        .iter()
        .rev()
        .map(|p| Piece {
            span: (p.span.1, p.span.0),
            ..*p
        })
        .collect();
    let joined = pieces.pop().map(|p| p.span.0);
    for (k, piece) in ahead.into_iter().enumerate() {
        match (k, joined) {
            (0, Some(left)) => pieces.push(Piece {
                span: (left, piece.span.1),
                ..piece
            }),
            _ => pieces.push(piece),
        }
    }
    for pair in pieces.windows(2) {
        if pair[0].holds == pair[1].holds
            || pieces.iter().filter(|p| p.holds == pair[1].holds).count() > 1
        {
            return Err(refuse("the ball comes back to what it held before"));
        }
    }
    if pieces.len() < 2 {
        // One piece is the rolled round's own run, refused there.
        return Ok(None);
    }
    prove(model, solid, edge, &crease, &pieces, &refuse)?;
    let Some(first) = host(model, &sides[0].face, tol)? else {
        return Ok(None);
    };
    let Some(second) = host(model, &sides[1].face, tol)? else {
        return Ok(None);
    };
    Ok(Some(Carried {
        pieces,
        sides,
        guide: curve,
        range,
        e0,
        t,
        length,
        hosts: [first, second],
    }))
}

/// The pieces from the ball holding `holds` at its track's parameter
/// `from` to where its section reaches `target` along the crease, each
/// ending at its first switch toward `target`, with spans from `from`'s
/// side to `target`'s. No piece holds what one before it held, so the walk
/// ends.
fn walk(
    crease: &Crease<'_>,
    mut holds: [Hold; 2],
    mut from: f64,
    target: f64,
) -> Result<Vec<Piece>, &'static str> {
    let conf = crease.conf();
    let mut pieces: Vec<Piece> = Vec::new();
    loop {
        if pieces.iter().any(|p| p.holds == holds) {
            return Err("the ball comes back to what it held before");
        }
        let track = crease.track(holds)?;
        let Some((x_from, rate)) = crease.crossing(&track, from) else {
            return Err("the ball's sections cross");
        };
        // The way along the crease toward the target, and the way of the
        // track's parameter that moves the section so.
        let way = (target - x_from).signum();
        let sense = rate.signum() * way;
        if rate == 0.0 || way == 0.0 {
            return Err("the ball's sections cross");
        }
        let mut switches = crease.switches(holds, &track, from, sense)?;
        switches.sort_by(|a, b| ((a.0 - from) * sense).total_cmp(&((b.0 - from) * sense)));
        // An edge holds the ball only between the side's normal and the
        // beyond face's: a piece leaving that arc at either end turns its
        // sections back.
        if let (Some(i), Track::Ellipse { .. }) = (Crease::odd_side(holds), track) {
            let Hold::Edge(j) = holds[i] else {
                return Err("the ball's holds do not match its track");
            };
            let (.., phi3) = crease.edge_frame(i, j)?;
            let step = from + sense * conf / track.speed();
            if switches.is_empty() || !(0.0..=phi3).contains(&step) {
                return Err("the ball's sections cross");
            }
        }
        // How far the piece may run: to its first switch, or half a turn
        // where its parameter is an angle and nothing switches.
        let limit = match (switches.first(), track) {
            (Some((s, _)), _) => *s,
            (None, Track::Line { .. }) => {
                from + sense * (crease.length * 2.0 + crease.radius * 4.0)
            }
            (None, _) => from + sense * core::f64::consts::PI,
        };
        // The crossing moves one way along the crease from `from` to the
        // limit, and the piece ends where it reaches the target first.
        let at = |i: usize| {
            #[allow(clippy::cast_precision_loss)]
            let f = (i as f64) / (SEARCH as f64);
            from + (limit - from) * f
        };
        let mut last = x_from;
        let mut end: Option<f64> = None;
        for i in 1..=SEARCH {
            let Some((x, _)) = crease.crossing(&track, at(i)) else {
                return Err("the ball's sections cross");
            };
            if (x - last) * way <= 0.0 {
                return Err("the ball's sections cross");
            }
            if (x - target) * way >= 0.0 {
                let (mut a, mut b) = (at(i - 1), at(i));
                for _ in 0..80 {
                    let m = (a + b) * 0.5;
                    let xm = crease.crossing(&track, m).map_or(f64::NAN, |c| c.0);
                    if (xm - target) * way < 0.0 {
                        a = m;
                    } else {
                        b = m;
                    }
                }
                end = Some((a + b) * 0.5);
                break;
            }
            last = x;
        }
        if let Some(end) = end {
            pieces.push(Piece {
                holds,
                track,
                span: (from, end),
            });
            return Ok(pieces);
        }
        let Some((switch, next)) = switches.into_iter().next() else {
            return Err("the ball does not reach the crease's end");
        };
        pieces.push(Piece {
            holds,
            track,
            span: (from, switch),
        });
        let next = match next {
            Next::Holds(next) => next,
            Next::Refuse(why) => return Err(why),
        };
        // The next piece starts where this one ends, on its own track,
        // which passes through the same centre tangent to this one.
        let c = track.centre(switch);
        let next_track = crease.track(next)?;
        let s = next_track.at(c);
        if next_track.centre(s).distance(c) > conf {
            return Err("its pieces do not meet where it changes hold");
        }
        holds = next;
        from = s;
    }
}

/// A face of the solid as the start's check of free space reads it: a
/// planar face exactly, a curved one by its box.
enum Reach {
    Flat(Flat),
    Bounds(ogeom_math::Aabb),
}

/// The one ball whose section crosses the crease at `x` along it, and the
/// piece it is on: its holds and its track's parameter there. `None` where
/// no ball, or more than one, does.
fn start(crease: &Crease<'_>, faces: &[Reach], x: f64) -> Option<([Hold; 2], f64)> {
    let (conf, radius) = (crease.conf(), crease.radius);
    let mut found: Vec<([Hold; 2], f64)> = Vec::new();
    // A face's contact inside it, clear of its boundary.
    let holds_face = |i: usize, p: Point| -> bool {
        let side = &crease.sides[i];
        let sides: Vec<[Point; 2]> = side
            .edges
            .iter()
            .map(|e| e.ends)
            .chain(std::iter::once([
                crease.e0,
                crease.e0 + crease.t * crease.length,
            ]))
            .collect();
        (p - crease.e0).dot(side.into) > conf
            && inside(&sides, p, basis(side.normal))
            && sides
                .iter()
                .all(|[a, b]| crate::rolled::segment_distance(p, [*a, *b]) > conf)
    };
    // A ball clear of every face of the solid.
    let free = |c: Point| -> bool {
        faces.iter().all(|face| match face {
            Reach::Flat(flat) => face_distance(flat, c) >= radius - conf,
            Reach::Bounds(bounds) => bounds.distance_to(c) > radius,
        })
    };
    // The edge of each face met first crossing it from the crease at `x`,
    // square to the crease: the far edge there, the only one the ball
    // rounding the crease rests on, or ends on at a corner.
    let far: Vec<Option<usize>> = (0..2)
        .map(|i| {
            let side = &crease.sides[i];
            let from = crease.e0 + crease.t * x;
            let mut best: Option<(f64, usize)> = None;
            for (j, edge) in side.edges.iter().enumerate() {
                let [p, q] = edge.ends;
                let w = q - p;
                let den = side.into.cross(w).dot(side.normal);
                if den.abs() <= crease.tol.angular() * w.magnitude() {
                    continue;
                }
                let along = (p - from).cross(w).dot(side.normal) / den;
                let at = (p - from).cross(side.into).dot(side.normal) / den;
                if along > conf && (0.0..=1.0).contains(&at) && best.is_none_or(|(b, _)| along < b)
                {
                    best = Some((along, j));
                }
            }
            best.map(|(_, j)| j)
        })
        .collect();
    let x1 = crease.e0 + crease.t * crease.length;
    let mut candidates: Vec<[Hold; 2]> = vec![[Hold::Face, Hold::Face]];
    for i in 0..2 {
        let Some(j) = far[i] else {
            continue;
        };
        let edge = &crease.sides[i].edges[j];
        let mut h = [Hold::Face; 2];
        h[i] = Hold::Edge(j);
        candidates.push(h);
        let mut corners: Vec<Point> = Vec::new();
        {
            for w in edge.ends {
                if w.distance(crease.e0) > conf
                    && w.distance(x1) > conf
                    && corners.iter().all(|c| c.distance(w) > conf)
                {
                    corners.push(w);
                }
            }
        }
        for w in corners {
            let mut h = [Hold::Face; 2];
            h[i] = Hold::Corner(w);
            candidates.push(h);
        }
    }
    for holds in candidates {
        let Ok(track) = crease.track(holds) else {
            continue;
        };
        // The parameter range a piece of this kind holds on: the line's
        // whole reach, the ellipse's arc between the side's normal and the
        // beyond face's, the circle's whole turn.
        let pi = core::f64::consts::PI;
        let (lo, hi) = match (Crease::odd_side(holds), &track) {
            (None, _) => (-crease.length * 2.0, crease.length * 2.0),
            (Some(i), Track::Ellipse { .. }) => {
                let Hold::Edge(j) = holds[i] else { continue };
                let Ok((.., phi3)) = crease.edge_frame(i, j) else {
                    continue;
                };
                (0.0, phi3)
            }
            (Some(_), Track::Circle { .. }) => (-pi, pi),
            _ => continue,
        };
        #[allow(clippy::cast_precision_loss)]
        let at = |k: usize| lo + (hi - lo) * (k as f64) / (SEARCH as f64);
        let mut previous: Option<f64> = None;
        for k in 0..=SEARCH {
            let along = crease.crossing(&track, at(k)).map(|c| c.0 - x);
            if let (Some(p), Some(q)) = (previous, along)
                && (p < 0.0) != (q < 0.0)
                && (p - q).abs() < crease.length
            {
                let (mut a, mut b) = (at(k - 1), at(k));
                for _ in 0..80 {
                    let m = (a + b) * 0.5;
                    let xm = crease.crossing(&track, m).map_or(f64::NAN, |c| c.0 - x);
                    if (xm > 0.0) == (p > 0.0) {
                        a = m;
                    } else {
                        b = m;
                    }
                }
                let s = (a + b) * 0.5;
                let c = track.centre(s);
                let held = (0..2).all(|i| match holds[i] {
                    Hold::Face => holds_face(i, crease.contact(i, Hold::Face, c)),
                    Hold::Edge(j) => {
                        let [p, q] = crease.sides[i].edges[j].ends;
                        let foot = crease.contact(i, holds[i], c);
                        let w = (foot - p).dot(q - p) / (q - p).dot(q - p);
                        w * p.distance(q) > conf && (1.0 - w) * p.distance(q) > conf
                    }
                    // The corner alone: the ball's radius to it turned away
                    // from both of the face's edges leaving it.
                    Hold::Corner(w) => crease.sides[i].edges.iter().all(|e| {
                        let [p, q] = e.ends;
                        let d = if p.distance(w) <= conf {
                            q - p
                        } else if q.distance(w) <= conf {
                            p - q
                        } else {
                            return true;
                        };
                        (c - w).dot(d) / d.magnitude() < -conf
                    }),
                }) && free(c);
                if held {
                    found.push((holds, s));
                }
            }
            previous = along;
        }
    }
    match found.as_slice() {
        [one] => Some(*one),
        _ => None,
    }
}

/// The proofs the chain stands on, each refused by name: the sections do
/// not cross, the contacts stay on their faces and edges, and every ball
/// lies in free space.
#[allow(
    clippy::too_many_lines,
    reason = "one proof per condition, read end to end"
)]
fn prove(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    crease: &Crease<'_>,
    pieces: &[Piece],
    refuse: &dyn Fn(&str) -> OgeomError,
) -> OgeomResult<()> {
    let (conf, radius) = (crease.conf(), crease.radius);
    let station = |piece: &Piece, i: usize, of: usize| {
        #[allow(clippy::cast_precision_loss)]
        let f = (i as f64) / (of as f64);
        piece.span.0 + (piece.span.1 - piece.span.0) * f
    };
    // The sections do not cross: within each piece the crease's crossing
    // moves one way, and each section's corners lie inside the curvature
    // of the path, so no point of one section lies on another's plane.
    let mut last: Option<f64> = None;
    for piece in pieces {
        let fine = STATIONS * 16;
        for i in 0..=fine {
            let s = station(piece, i, fine);
            let Some((x, _)) = crease.crossing(&piece.track, s) else {
                return Err(refuse("the ball's sections cross"));
            };
            if last.is_some_and(|l| x < l - conf) {
                return Err(refuse("the ball's sections cross"));
            }
            last = Some(x);
            let c = piece.track.centre(s);
            let corners = [
                crease.e0 + crease.t * x,
                crease.contact(0, piece.holds[0], c),
                crease.contact(1, piece.holds[1], c),
            ];
            // The plane's line of rest stands where `c''.(p - c) = |c'|^2`;
            // a corner on it within the confusion is touched, not crossed.
            let dd = piece.track.d2(s);
            let bend = corners
                .iter()
                .map(|p| dd.dot(*p - c))
                .fold(f64::NEG_INFINITY, f64::max);
            let speed = piece.track.d1(s);
            let across = (dd - speed * (dd.dot(speed) / speed.dot(speed))).magnitude();
            if bend - speed.dot(speed) > conf * across {
                return Err(refuse("the ball's sections cross"));
            }
        }
    }

    // The contacts on an edge stay on it and move one way along it, the
    // foot's reach over a piece read at its ends and where it turns.
    for piece in pieces {
        for i in 0..2 {
            let Hold::Edge(j) = piece.holds[i] else {
                continue;
            };
            let [p, q] = crease.sides[i].edges[j].ends;
            let reach = p.distance(q);
            let d = (q - p) / reach;
            let foot =
                |s: f64| (crease.contact(i, piece.holds[i], piece.track.centre(s)) - p).dot(d);
            let (a, b) = piece.span;
            let mut values = vec![foot(a), foot(b)];
            if let Track::Ellipse { au, av, .. } = piece.track {
                // The foot turns where `l'(s) = 0`: `sin s au = cos s av`.
                let turn = av.atan2(au);
                let tau = core::f64::consts::PI;
                let (lo, hi) = (a.min(b), a.max(b));
                let mut s = turn - tau * ((turn - lo) / tau).ceil();
                while s <= hi {
                    if s >= lo {
                        values.push(foot(s));
                    }
                    s += tau;
                }
            }
            let lo = values.iter().copied().fold(f64::INFINITY, f64::min);
            let hi = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            if lo < -conf || hi > reach + conf {
                return Err(refuse("the ball runs past that edge's end"));
            }
            let (fa, fb) = (values[0], values[1]);
            if (lo - fa.min(fb)).abs() > conf || (hi - fa.max(fb)).abs() > conf {
                return Err(refuse("the ball's contact turns back along that edge"));
            }
        }
    }

    // The legs stand on their faces: each face holds the stretch of its
    // plane from the crease to its line of contact. No edge of the face
    // but the crease and the edges the ball rests on enters that stretch,
    // and each face lies on its stretch's side of the crease. A curved line
    // of contact is read with chords close enough that none strays from it
    // by the confusion.
    let x0 = crease.e0;
    let x1 = crease.e0 + crease.t * crease.length;
    for (i, side) in crease.sides.iter().enumerate() {
        let mut ring = vec![x0, x1];
        for piece in pieces.iter().rev() {
            let (a, b) = piece.span;
            let curved = piece.holds[i] == Hold::Face && !matches!(piece.track, Track::Line { .. });
            let chords = if curved {
                let bend = piece.track.speed() * (b - a).abs();
                let steps = ((bend * (b - a).abs() / (8.0 * conf)).sqrt()).ceil();
                if !steps.is_finite() || steps > MAX_CHORDS {
                    return Err(refuse("the line of contact is too long to read"));
                }
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let steps = (steps as usize).max(STATIONS);
                steps
            } else {
                1
            };
            for k in (0..=chords).rev() {
                #[allow(clippy::cast_precision_loss)]
                let s = a + (b - a) * (k as f64) / (chords as f64);
                let p = crease.contact(i, piece.holds[i], piece.track.centre(s));
                if ring.last().is_none_or(|l: &Point| l.distance(p) > conf) {
                    ring.push(p);
                }
            }
        }
        let held: Vec<&Shape> = pieces
            .iter()
            .filter_map(|p| match p.holds[i] {
                Hold::Edge(j) => Some(&side.edges[j].shape),
                _ => None,
            })
            .collect();
        for e in explore_unique(model, &side.face, ShapeType::Edge)? {
            if same_occurrence(model, &e, edge, crease.tol)
                || held
                    .iter()
                    .any(|h| same_occurrence(model, &e, h, crease.tol))
            {
                continue;
            }
            let Some(segment) = segment(model, &e, crease.tol)? else {
                return Err(refuse("a face the round stands on has a curved edge"));
            };
            if enters(&ring, side.normal, segment, conf) {
                return Err(refuse(if pieces.iter().any(|p| p.holds[i] == Hold::Face) {
                    "the ball's contact runs off the face it touches"
                } else {
                    "the round's foot runs off the narrow face"
                }));
            }
        }
        // Each contact stands on its face's side of the crease.
        for piece in pieces {
            let s = (piece.span.0 + piece.span.1) * 0.5;
            let Some((x, _)) = crease.crossing(&piece.track, s) else {
                return Err(refuse("the ball's sections cross"));
            };
            let p = crease.contact(i, piece.holds[i], piece.track.centre(s));
            if (p - (crease.e0 + crease.t * x)).dot(side.into) <= 0.0 {
                return Err(refuse("the round stands off its faces"));
            }
        }
    }

    // Every ball lies in free space: no face of the solid comes inside it.
    // A face is a radius off every centre of a piece where it lies behind
    // the plane through one of the ball's contacts square to the ball's
    // radius there: behind the touched face's plane, or behind the plane
    // through the contact on an edge or a corner square to the ball's
    // radius to it, read over the piece in closed form. Every other face is
    // read at the stations, its distance falling by at most the centre's
    // travel between them, and an interval is halved until that travel
    // fits under its clearance.
    let mut flats: Vec<(Shape, Option<Flat>)> = Vec::new();
    for face in explore_unique(model, solid, ShapeType::Face)? {
        let flat = flat(model, &face, crease.tol)?;
        flats.push((face, flat));
    }
    for piece in pieces {
        let speed = piece.track.speed();
        let (sa, sb) = piece.span;
        for (face, found) in &flats {
            let Some(found) = found else {
                let bounds = ogeom_algo::shape_bounds(model, face, crease.tol)?;
                #[allow(clippy::cast_precision_loss)]
                let step = speed * (sb - sa).abs() / (STATIONS as f64);
                for i in 0..=STATIONS {
                    let c = piece.track.centre(station(piece, i, STATIONS));
                    if bounds.distance_to(c) <= radius + step {
                        return Err(refuse("a curved face stands within reach of the ball"));
                    }
                }
                continue;
            };
            let vertices: Vec<Point> = found.sides.iter().flatten().copied().collect();
            if behind(crease, piece, &vertices) {
                continue;
            }
            let clearance = |s: f64| face_distance(found, piece.track.centre(s)) - radius;
            let mut intervals: Vec<(f64, f64, f64, f64, u32)> = Vec::new();
            for i in 0..STATIONS {
                let (a, b) = (station(piece, i, STATIONS), station(piece, i + 1, STATIONS));
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
    }
    Ok(())
}

/// Whether every vertex of a planar face lies behind the plane through one
/// of `piece`'s contacts square to the ball's radius there, over the whole
/// piece: then every point of the face does, and stands at least a radius
/// from every centre.
fn behind(crease: &Crease<'_>, piece: &Piece, vertices: &[Point]) -> bool {
    let conf = crease.conf();
    let (a, b) = piece.span;
    (0..2).any(|i| match piece.holds[i] {
        Hold::Face => {
            let n = crease.sides[i].normal;
            vertices.iter().all(|v| (*v - crease.e0).dot(n) <= conf)
        }
        Hold::Edge(j) => {
            // From the foot, the ball's radius is `cos s u + sin s v`, square
            // to the edge.
            let Ok((p, _, _, u, v, _)) = crease.edge_frame(i, j) else {
                return false;
            };
            vertices
                .iter()
                .all(|x| trig_max(u.dot(*x - p), v.dot(*x - p), a, b) <= conf)
        }
        Hold::Corner(w) => {
            let Track::Circle {
                centre,
                e1,
                e2,
                rho,
            } = piece.track
            else {
                return false;
            };
            vertices.iter().all(|x| {
                let base = (centre - w).dot(*x - w);
                base + trig_max(rho * e1.dot(*x - w), rho * e2.dot(*x - w), a, b)
                    <= conf * crease.radius
            })
        }
    })
}

/// Round `edge` with the ball carried across its changes of hold, where
/// [`carried_run`] finds that run.
///
/// # Errors
///
/// As [`carried_run`], and [`OgeomError::NotDone`] where a band's fit
/// stands further off its balls' envelope than its target.
pub(crate) fn carried_band(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    run: &Carried,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let crease = Crease {
        sides: &run.sides,
        e0: run.e0,
        t: run.t,
        length: run.length,
        radius,
        tol,
    };
    let (lo, hi) = run.range;
    let mut blends = Vec::with_capacity(run.pieces.len());
    for piece in &run.pieces {
        let (sa, sb) = piece.span;
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
            #[allow(clippy::cast_precision_loss)]
            let s = sa + (sb - sa) * (i as f64) / (STATIONS as f64);
            let Some((x, _)) = crease.crossing(&piece.track, s) else {
                ogeom_core::ogeom_bail!(Construction, "a carried round's sections cross");
            };
            let c = piece.track.centre(s);
            let first = crease.contact(0, piece.holds[0], c);
            let second = crease.contact(1, piece.holds[1], c);
            blend.along.push(lo + (hi - lo) * x / run.length);
            blend.spine.push(c);
            blend.touch_first.push(first);
            blend.touch_second.push(second);
            blend
                .on_first
                .push(ogeom_algo::project_on_surface(&run.hosts[0].0, first, 16, tol)?.parameters);
            blend
                .on_second
                .push(ogeom_algo::project_on_surface(&run.hosts[1].0, second, 16, tol)?.parameters);
        }
        blends.push(blend);
    }
    let rails: Vec<[Rail; 2]> = run
        .pieces
        .iter()
        .map(|piece| {
            let rail = |i: usize| match piece.holds[i] {
                Hold::Face => Rail::Plane(run.e0, run.sides[i].normal),
                Hold::Edge(j) => {
                    let [p, q] = run.sides[i].edges[j].ends;
                    Rail::Line(p, (q - p) / p.distance(q))
                }
                Hold::Corner(w) => Rail::Point(w),
            };
            [rail(0), rail(1)]
        })
        .collect();
    // The distance from a point to a piece's envelope: from its nearest
    // centre on the piece, found by golden section about the nearest
    // station.
    let offs: Vec<Box<dyn Fn(Point) -> f64 + '_>> = run
        .pieces
        .iter()
        .zip(&blends)
        .map(|(piece, blend)| {
            let piece = *piece;
            let spine = blend.spine.clone();
            Box::new(move |p: Point| -> f64 {
                let n = spine.len();
                let nearest = (0..n)
                    .min_by(|&i, &j| p.distance(spine[i]).total_cmp(&p.distance(spine[j])))
                    .unwrap_or(0);
                let (sa, sb) = piece.span;
                #[allow(clippy::cast_precision_loss)]
                let at = |i: f64| sa + (sb - sa) * i / ((n - 1) as f64);
                #[allow(clippy::cast_precision_loss)]
                let (mut a, mut b) = (at(nearest as f64 - 1.0), at(nearest as f64 + 1.0));
                let f = |s: f64| p.distance(piece.track.centre(s));
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
            }) as Box<dyn Fn(Point) -> f64 + '_>
        })
        .collect();
    let exact: Vec<ExactBand<'_>> = rails
        .iter()
        .zip(&offs)
        .map(|(rails, off)| ExactBand {
            rails: *rails,
            off: off.as_ref(),
        })
        .collect();
    let [(first, first_sign), (second, second_sign)] = &run.hosts;
    crate::marched::build_chain_band(
        model,
        solid,
        Some(edge),
        &blends,
        &exact,
        &run.guide,
        [(first, *first_sign), (second, *second_sign)],
        radius,
        tol,
    )
}
