//! Hold edge (ADR 0001, amendment "overflow handling"): a round whose
//! contact on one of its faces, `a`, would run past a sharp edge of that
//! face keeps the edge. The ball touches the other face, `b`, and rests on
//! the nearest point of `a`'s region: `a` itself where its contact lies on
//! it (the band), else a sharp edge of `a` (the ball's foot inside the
//! edge) or a vertex where two such edges meet (the foot at the vertex).
//! Its centre stays on `b`'s offset at the radius from that seat, and the
//! round is the canal surface those balls sweep, one piece per seat.
//!
//! Built for a straight edge, or a chain of edges along one line, between
//! two planes, the held edges lines or circles lying in `a`'s plane, the
//! round ending on end faces within 60 degrees of square to its edge. Each
//! ball's arc runs from its contact on `b` to its foot, on its
//! characteristic circle. Where the seat changes the ball touches both
//! seats at once and the two offsets are tangent, so neighbouring pieces
//! share that ball's characteristic circle and meet G1 on it. From the band
//! to an edge the spine runs along the edge there, so that circle is the
//! band's own section; from an edge to a vertex it is the ball's arc in its
//! plane square to the spine. A piece is fitted through its balls' arcs,
//! or, where it reaches an end face, through its sections by planes turning
//! to that face's plane. Past an end face the centres run on as they are
//! and the end face's plane cuts the round (amendment "M5f"). A ball
//! resting on two edges at once, its centre turning a corner, is not held.
//!
//! The region between the edge and the round is one closed tool, cut from
//! the solid; it is built only where nothing else of the solid stands
//! within its reach, the cut takes the tool's whole volume, the result is
//! valid and the held edges are still there under every ball's foot.

use core::f64::consts::TAU;

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve3d as _, CylinderSurface, SurfaceGeometry};
use ogeom_math::{Aabb, Cylinder, Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape};

use super::curved::{Curved, SheetCurve, read};
use super::overflow::{ToolFace, arc, boundary_from, closed, segment};
use super::trimmed::{faces_of, outward};
use super::{Id, id};

/// Planar sections through each fitted piece.
const STATIONS: usize = 160;
/// Points across each section.
const ACROSS: usize = 65;
/// Points along an edge, checked against the feet.
const SAMPLES: usize = 64;
/// Stations along the edge where the ball's seat is read.
const READINGS: usize = 1024;

/// The curve of a sharp edge of `a`: a line through `at` along `along`, or
/// the circle of `radius` about `centre` in the plane square to `normal`,
/// its angle measured from `x` toward `y`.
#[derive(Debug, Clone, Copy)]
enum Cliff {
    Line {
        at: Point,
        along: Vector,
    },
    Circle {
        centre: Point,
        normal: Vector,
        radius: f64,
        x: Vector,
        y: Vector,
    },
}

impl Cliff {
    /// The point of the curve nearest `p`, `None` on a circle's axis.
    fn foot(self, p: Point) -> Option<Point> {
        match self {
            Self::Line { at, along } => Some(at + along * (p - at).dot(along)),
            Self::Circle {
                centre,
                normal,
                radius,
                ..
            } => {
                let q = p - centre;
                let radial = q - normal * q.dot(normal);
                let m = radial.magnitude();
                (m > 0.0).then(|| centre + radial * (radius / m))
            }
        }
    }

    /// Where `p` lies along the curve: along a line, or its angle about a
    /// circle.
    fn param(self, p: Point) -> f64 {
        match self {
            Self::Line { at, along } => (p - at).dot(along),
            Self::Circle { centre, x, y, .. } => {
                let q = p - centre;
                q.dot(y).atan2(q.dot(x))
            }
        }
    }
}

/// Where the ball rests besides `b`: on `a`, on a sharp edge of it (by
/// index into the hold's edges), or on a vertex where two of them meet (by
/// index into its corners).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Seat {
    Face,
    Rim(usize),
    Corner(usize),
}

/// A sharp edge of `a` the ball may rest on: its curve and the stretch of
/// it the edge covers, from `range.0` at its first end to `range.1` (an
/// angle growing round a circle), and at each end the corner where it
/// stops, or `None` where it runs on to an end face, past which the round
/// runs on as the edge's curve does.
#[derive(Debug, Clone)]
struct Rim {
    id: Id,
    cliff: Cliff,
    range: (f64, f64),
    whole: bool,
    stops: [Option<usize>; 2],
}

/// Where along a rim a parameter lies.
enum Along {
    Inside,
    Stopped(usize),
}

impl Rim {
    /// `t` placed on the rim: within its stretch, or carried on past an
    /// end it runs on at, else stopped at the nearer corner.
    fn place(&self, t: f64) -> Along {
        match self.cliff {
            Cliff::Line { .. } => {
                if t < self.range.0
                    && let Some(c) = self.stops[0]
                {
                    return Along::Stopped(c);
                }
                if t > self.range.1
                    && let Some(c) = self.stops[1]
                {
                    return Along::Stopped(c);
                }
                Along::Inside
            }
            Cliff::Circle { .. } => {
                if self.whole {
                    return Along::Inside;
                }
                let rel = (t - self.range.0).rem_euclid(TAU);
                let sweep = self.range.1 - self.range.0;
                if rel <= sweep {
                    return Along::Inside;
                }
                let (over, under) = (rel - sweep, TAU - rel);
                let end = usize::from(over < under);
                self.stops[end].map_or(Along::Inside, Along::Stopped)
            }
        }
    }

    /// Whether `t` lies within the edge itself, not past an end.
    fn covers(&self, t: f64, slack: f64) -> bool {
        match self.cliff {
            Cliff::Line { .. } => t >= self.range.0 - slack && t <= self.range.1 + slack,
            Cliff::Circle { .. } => {
                self.whole
                    || (t - self.range.0 + slack).rem_euclid(TAU)
                        <= self.range.1 - self.range.0 + 2.0 * slack
            }
        }
    }
}

/// A vertex where two sharp edges of `a` meet.
#[derive(Debug, Clone, Copy)]
struct Corner {
    id: Id,
    at: Point,
}

/// A stretch of the edge over which the ball keeps one seat.
#[derive(Debug, Clone, Copy)]
struct Piece {
    seat: Seat,
    from: f64,
    to: f64,
}

/// A straight round whose contact on face `a` runs past sharp edges of it,
/// read piece by piece.
#[derive(Debug, Clone)]
pub(crate) struct Hold {
    curved: Curved,
    /// The round's edges, end to end along one line.
    edges: Vec<Id>,
    /// A face of each plane, and every face of each along the round.
    a: usize,
    b: usize,
    a_faces: Vec<usize>,
    b_faces: Vec<usize>,
    ends: [usize; 2],
    /// The end faces' outward normals.
    end_normals: [Vector; 2],
    /// The edge's start, its unit direction and its length.
    origin: Point,
    along: Vector,
    len: f64,
    /// The outward normals of `a` and `b`.
    na: Vector,
    nb: Vector,
    /// In `a`, square to the edge, pointing from the edge into `a`.
    across: Vector,
    radius: f64,
    rims: Vec<Rim>,
    corners: Vec<Corner>,
    pieces: Vec<Piece>,
    /// Whether a stretch held reaches an end of the band.
    pub(crate) at_end: bool,
}

/// A ball of the round: its centre, its contact on `b`, and its touch on
/// `a` or its foot on a held edge or vertex.
#[derive(Debug, Clone, Copy)]
struct Ball {
    centre: Point,
    on_b: Point,
    on_a: Point,
}

impl Hold {
    fn on_edge(&self, x: f64) -> Point {
        self.origin + self.along * x
    }

    fn band(&self, x: f64) -> Point {
        let k = -self.radius / (1.0 + self.na.dot(self.nb));
        self.on_edge(x) + (self.na + self.nb) * k
    }

    /// How far into `a` the band's contact lies from the edge.
    fn contact_depth(&self) -> f64 {
        (self.band(0.0) + self.na * self.radius - self.origin).dot(self.across)
    }

    /// Whether the point `x` along the edge and `y` into `a` lies on `a`:
    /// whether the walk to it from the edge crosses the sharp edges of `a`
    /// an even number of times.
    fn inside_at(&self, x: f64, y: f64) -> bool {
        let start = self.contact_depth() * 1e-9;
        if y <= start {
            return true;
        }
        let base = self.on_edge(x);
        let mut crossings = 0;
        for rim in &self.rims {
            match rim.cliff {
                Cliff::Line { at, along } => {
                    let far = 10.0 * (self.len + self.radius);
                    let lo = if rim.stops[0].is_none() {
                        rim.range.0 - far
                    } else {
                        rim.range.0
                    };
                    let hi = if rim.stops[1].is_none() {
                        rim.range.1 + far
                    } else {
                        rim.range.1
                    };
                    let (p, q) = (at + along * lo, at + along * hi);
                    let (px, qx) = (
                        (p - self.origin).dot(self.along) - x,
                        (q - self.origin).dot(self.along) - x,
                    );
                    if px * qx >= 0.0 {
                        continue;
                    }
                    let (py, qy) = (
                        (p - self.origin).dot(self.across),
                        (q - self.origin).dot(self.across),
                    );
                    let cy = py + (qy - py) * px / (px - qx);
                    if cy > start && cy < y {
                        crossings += 1;
                    }
                }
                Cliff::Circle { centre, radius, .. } => {
                    let d = base - centre;
                    let half = d.dot(self.across);
                    let disc = half * half - (d.dot(d) - radius * radius);
                    if disc <= 0.0 {
                        continue;
                    }
                    for root in [-half - disc.sqrt(), -half + disc.sqrt()] {
                        if root <= start || root >= y {
                            continue;
                        }
                        let p = base + self.across * root;
                        if matches!(rim.place(rim.cliff.param(p)), Along::Inside) {
                            crossings += 1;
                        }
                    }
                }
            }
        }
        crossings % 2 == 0
    }

    /// The distance from `p` to its seat, the seat's nearest point and the
    /// seat: the one given, or the nearest of `a`'s region.
    fn support(&self, p: Point, seat: Option<Seat>) -> Option<(f64, Point, Seat)> {
        let on_plane = |p: Point| p - self.na * (p - self.origin).dot(self.na);
        match seat {
            Some(Seat::Face) => {
                let f = on_plane(p);
                Some((p.distance(f), f, Seat::Face))
            }
            Some(Seat::Rim(i)) => {
                let f = self.rims[i].cliff.foot(p)?;
                Some((p.distance(f), f, Seat::Rim(i)))
            }
            Some(Seat::Corner(k)) => {
                let f = self.corners[k].at;
                Some((p.distance(f), f, Seat::Corner(k)))
            }
            None => {
                let q = p - self.origin;
                if self.inside_at(q.dot(self.along), q.dot(self.across)) {
                    let f = on_plane(p);
                    return Some((p.distance(f), f, Seat::Face));
                }
                let mut best: Option<(f64, Point, Seat)> = None;
                for (i, rim) in self.rims.iter().enumerate() {
                    let Some(f) = rim.cliff.foot(p) else {
                        continue;
                    };
                    let (f, seat) = match rim.place(rim.cliff.param(f)) {
                        Along::Inside => (f, Seat::Rim(i)),
                        Along::Stopped(k) => (self.corners[k].at, Seat::Corner(k)),
                    };
                    let d = p.distance(f);
                    if best.is_none_or(|b| d < b.0) {
                        best = Some((d, f, seat));
                    }
                }
                best
            }
        }
    }

    /// The ball at `x` along the spine on `seat`, or on whatever it rests
    /// on there: on `b`'s offset at the radius from its seat, found from
    /// the band's centre toward `a`'s plane.
    fn ball(&self, x: f64, seat: Option<Seat>) -> Option<(Ball, Seat)> {
        let r = self.radius;
        let c = self.band(x);
        let band = Ball {
            centre: c,
            on_b: c + self.nb * r,
            on_a: c + self.na * r,
        };
        if seat == Some(Seat::Face) || (seat.is_none() && self.inside_at(x, self.contact_depth())) {
            return Some((band, Seat::Face));
        }
        let g = self.along.cross(self.nb);
        let g = if g.dot(self.na) < 0.0 { -g } else { g };
        let reach = r / g.dot(self.na);
        let f = |s: f64| self.support(c + g * s, seat).map(|(d, _, _)| d - r);
        let s = root(f, 0.0, reach)?;
        let centre = c + g * s;
        let (_, foot, found) = self.support(centre, seat)?;
        Some((
            Ball {
                centre,
                on_b: centre + self.nb * r,
                on_a: foot,
            },
            found,
        ))
    }

    /// The point a fraction `u` along the arc of the ball at `t` on
    /// `seat`, from its contact on `b` to its touch on `a` or its foot.
    fn point(&self, t: f64, u: f64, seat: Seat) -> Option<Point> {
        let (ball, _) = self.ball(t, Some(seat))?;
        let a = (ball.on_b - ball.centre) / self.radius;
        let e = (ball.on_a - ball.centre) / self.radius;
        let cross = a.cross(e);
        let m = cross.magnitude();
        if m <= 1e-12 {
            return None;
        }
        let sweep = a.dot(e).clamp(-1.0, 1.0).acos() * u;
        let w = (cross / m).cross(a);
        Some(ball.centre + (a * sweep.cos() + w * sweep.sin()) * self.radius)
    }

    /// The section of the piece on `seat` by the plane through `at` square
    /// to `normal`, near `x`, at `n` points from its contact on `b` to its
    /// foot: each point the one of its fraction along a ball's arc that
    /// lies in that plane, the ball found along the spine. The piece's
    /// balls run on past its ends, so a section at a seat change is the
    /// ball's there, which the two pieces share.
    fn section(
        &self,
        seat: Seat,
        at: Point,
        normal: Vector,
        x: f64,
        n: usize,
    ) -> Option<Vec<Point>> {
        let mut row = Vec::with_capacity(n);
        for k in 0..n {
            #[expect(clippy::cast_precision_loss, reason = "a few dozen points")]
            let u = k as f64 / (n - 1) as f64;
            let gap = |t: f64| self.point(t, u, seat).map(|p| (p - at).dot(normal));
            // The nearest bracket that holds the plane, walking out from
            // `x` both ways: a piece's balls reach only so far past its
            // ends, and a foot may trail its ball's station.
            let step = self.radius / 32.0;
            let mut bracket = None;
            let mut sides = [(x, gap(x)), (x, gap(x))];
            'walk: for k in 1..=96 {
                for (side, sign) in sides.iter_mut().zip([-1.0, 1.0]) {
                    let t = x + sign * step * f64::from(k);
                    let g = gap(t);
                    if let (Some(a), Some(b)) = (side.1, g)
                        && (a <= 0.0) != (b <= 0.0)
                    {
                        bracket = Some((side.0.min(t), side.0.max(t)));
                        break 'walk;
                    }
                    *side = (t, g);
                }
            }
            let (lo, hi) = bracket?;
            row.push(self.point(root(gap, lo, hi)?, u, seat)?);
        }
        Some(row)
    }

    /// The seat at `x`.
    fn seat_at(&self, x: f64) -> Option<Seat> {
        self.ball(x, None).map(|(_, s)| s)
    }

    /// The edge read into pieces of one seat each, every seat change one
    /// the round meets G1 across: the band onto an edge it crosses, or an
    /// edge onto a vertex where it stops. `None` where the seat cannot be
    /// read, or changes any other way (the ball resting on two edges at
    /// once, its centre turning a corner).
    fn read_pieces(&self) -> Option<Vec<Piece>> {
        #[expect(clippy::cast_precision_loss, reason = "a thousand readings")]
        let at = |k: usize| self.len * k as f64 / READINGS as f64;
        let mut groups: Vec<(Seat, usize, usize)> = Vec::new();
        for k in 0..=READINGS {
            let seat = self.seat_at(at(k))?;
            match groups.last_mut() {
                Some(g) if g.0 == seat => g.2 = k,
                _ => groups.push((seat, k, k)),
            }
        }
        let bisect = |inside: f64, outside: f64, seat: Seat| -> Option<f64> {
            let (mut i, mut o) = (inside, outside);
            for _ in 0..60 {
                let m = 0.5 * (i + o);
                if self.seat_at(m)? == seat {
                    i = m;
                } else {
                    o = m;
                }
            }
            Some(0.5 * (i + o))
        };
        let mut pieces: Vec<Piece> = Vec::new();
        let mut from = 0.0;
        for pair in groups.windows(2) {
            let (left, right) = (pair[0], pair[1]);
            let (x0, x1) = (at(left.2), at(right.1));
            let end = bisect(x0, x1, left.0)?;
            let start = bisect(x1, x0, right.0)?;
            pieces.push(Piece {
                seat: left.0,
                from,
                to: end,
            });
            if start - end > self.len * 1e-9 {
                // A seat between the readings: it holds over the gap.
                let seat = self.seat_at(0.5 * (end + start))?;
                for k in 1..8 {
                    let x = end + (start - end) * f64::from(k) / 8.0;
                    if self.seat_at(x)? != seat {
                        return None;
                    }
                }
                pieces.push(Piece {
                    seat,
                    from: end,
                    to: start,
                });
                from = start;
            } else {
                let mid = 0.5 * (end + start);
                if let Some(last) = pieces.last_mut() {
                    last.to = mid;
                }
                from = mid;
            }
        }
        let last = groups.last()?;
        pieces.push(Piece {
            seat: last.0,
            from,
            to: self.len,
        });
        let fits = pieces.windows(2).all(|w| match (w[0].seat, w[1].seat) {
            (Seat::Face, Seat::Rim(_)) | (Seat::Rim(_), Seat::Face) => true,
            (Seat::Rim(i), Seat::Corner(k)) | (Seat::Corner(k), Seat::Rim(i)) => {
                self.rims[i].stops.contains(&Some(k))
            }
            _ => false,
        });
        let short = pieces.iter().any(|p| p.to - p.from <= self.radius * 1e-4);
        (fits && !short).then_some(pieces)
    }

    /// The plane a piece ends in at break `k`, at `x`: the end face's at
    /// an end; square to the edge where the band meets an edge, the spine
    /// running along it there; else the plane of the ball's arc, square to
    /// the spine. Each normal points along the edge.
    fn plane_at(&self, x: f64, k: usize) -> Option<(Point, Vector)> {
        let last = self.pieces.len();
        if k == 0 {
            return Some((self.on_edge(x), -self.end_normals[0]));
        }
        if k == last {
            return Some((self.on_edge(x), self.end_normals[1]));
        }
        let square = self.pieces[k - 1].seat == Seat::Face || self.pieces[k].seat == Seat::Face;
        let seat = if k == 0 {
            self.pieces[0].seat
        } else {
            self.pieces[k - 1].seat
        };
        let (ball, _) = self.ball(x, Some(seat))?;
        if square {
            return Some((ball.centre, self.along));
        }
        let h = self.radius * 1e-5;
        let (ahead, _) = self.ball(x + h, Some(seat))?;
        let (behind, _) = self.ball(x - h, Some(seat))?;
        let t = ahead.centre - behind.centre;
        let m = t.magnitude();
        (m > 0.0).then(|| (ball.centre, t / m))
    }

    /// Where the end face at break `k` (the first or the last) cuts the
    /// round: its contact on `b` and its touch on `a` or its foot.
    fn end_cut(&self, k: usize) -> Option<(Point, Point)> {
        let last = self.pieces.len();
        let (x, piece) = if k == 0 {
            (0.0, self.pieces[0])
        } else {
            (self.len, self.pieces[last - 1])
        };
        let (at, normal) = self.plane_at(x, k)?;
        if piece.seat == Seat::Face {
            // The band's contact lines, straight along the edge, met by
            // the end face's plane.
            let (ball, _) = self.ball(0.0, Some(Seat::Face))?;
            let meet = |p: Point| p + self.along * ((at - p).dot(normal) / self.along.dot(normal));
            return Some((meet(ball.on_b), meet(ball.on_a)));
        }
        let row = self.section(piece.seat, at, normal, x, 2)?;
        Some((row[0], row[1]))
    }

    /// Whether every ball's foot on a held edge lies on one of `runs`,
    /// each the points along an edge of the result on that edge's curve.
    fn feet_on(&self, runs: &[(usize, Vec<Point>)]) -> bool {
        let slack = self.radius * 1e-6;
        self.pieces.iter().all(|piece| {
            let Seat::Rim(i) = piece.seat else {
                return true;
            };
            let rim = &self.rims[i];
            let Some((middle, _)) = self.ball(0.5 * (piece.from + piece.to), Some(piece.seat))
            else {
                return false;
            };
            // Angles about a circle measured from the middle foot's, so a
            // run does not straddle the cut at a half turn.
            let reference = rim.cliff.param(middle.on_a);
            let rel = |t: f64| match rim.cliff {
                Cliff::Line { .. } => t,
                Cliff::Circle { .. } => (t - reference + TAU / 2.0).rem_euclid(TAU) - TAU / 2.0,
            };
            let coords: Vec<Vec<f64>> = runs
                .iter()
                .filter(|(r, _)| *r == i)
                .map(|(_, ps)| ps.iter().map(|p| rel(rim.cliff.param(*p))).collect())
                .collect();
            (0..=32).all(|k| {
                let x = piece.from + (piece.to - piece.from) * f64::from(k) / 32.0;
                let Some((ball, _)) = self.ball(x, Some(piece.seat)) else {
                    return false;
                };
                let t = rim.cliff.param(ball.on_a);
                // A foot past the end face, on the edge's curve carried on,
                // is cut off with the round.
                if !rim.covers(t, slack) {
                    return true;
                }
                let c = rel(t);
                coords.iter().any(|run| {
                    run.windows(2).any(|w| {
                        (w[1] - w[0]).abs() < core::f64::consts::FRAC_PI_2
                            && c >= w[0].min(w[1]) - slack
                            && c <= w[0].max(w[1]) + slack
                    })
                })
            })
        })
    }
}

/// The root of `f` between `lo` and `hi`, where its signs differ, by
/// regula falsi with the Illinois step; `None` where they do not, or `f`
/// has no value on the way.
fn root(f: impl Fn(f64) -> Option<f64>, lo: f64, hi: f64) -> Option<f64> {
    let (mut a, mut b) = (lo, hi);
    let (mut fa, mut fb) = (f(a)?, f(b)?);
    if fa == 0.0 {
        return Some(a);
    }
    if fb == 0.0 {
        return Some(b);
    }
    if (fa > 0.0) == (fb > 0.0) {
        return None;
    }
    let mut side = 0;
    for _ in 0..100 {
        let x = (a * fb - b * fa) / (fb - fa);
        let fx = f(x)?;
        if fx == 0.0 || (b - a).abs() <= 1e-14 * (1.0 + x.abs()) {
            return Some(x);
        }
        if (fx > 0.0) == (fb > 0.0) {
            b = x;
            fb = fx;
            if side == -1 {
                fa *= 0.5;
            }
            side = -1;
        } else {
            a = x;
            fa = fx;
            if side == 1 {
                fb *= 0.5;
            }
            side = 1;
        }
        if fx.abs() <= 1e-15 {
            return Some(x);
        }
    }
    Some((a * fb - b * fa) / (fb - fa))
}

/// The curve of the edge `held` of the plane with normal `plane_normal`,
/// where it is a line or a circle lying in that plane.
fn cliff_of(curved: &Curved, held: Id, plane_normal: Vector, eps: f64) -> Option<Cliff> {
    let edge = &curved.edges[&held];
    match &edge.curve {
        Curve::Line(line) => {
            let axis = line.axis();
            let along = axis.direction.vector();
            (along.dot(plane_normal).abs() <= eps).then_some(Cliff::Line {
                at: axis.location,
                along,
            })
        }
        Curve::Circle(c) => {
            let circle = c.circle();
            let frame = circle.frame();
            let normal = frame.z().vector();
            (normal.cross(plane_normal).magnitude() <= eps).then_some(Cliff::Circle {
                centre: circle.centre(),
                normal,
                radius: circle.radius(),
                x: frame.x().vector(),
                y: frame.y().vector(),
            })
        }
        _ => None,
    }
}

/// The faces other than `ours` at vertex `v`.
fn others_at(curved: &Curved, v: Id, ours: &[usize]) -> Vec<usize> {
    (0..curved.faces.len())
        .filter(|fi| !ours.contains(fi))
        .filter(|fi| touches(curved, *fi, v))
        .collect()
}

/// Whether face `fi` has `v` among its vertices.
fn touches(curved: &Curved, fi: usize, v: Id) -> bool {
    curved.faces[fi].loops.iter().flatten().any(|(e, _)| {
        let c = &curved.edges[e];
        c.ends.0 == v || c.ends.1 == v
    })
}

/// The sharp edges of the faces of `a`'s plane near the band, as rims,
/// and the vertices where two of them meet, as corners; `None` where an
/// edge of theirs within reach of the band is neither a line nor a circle
/// in that plane. The round's own edges, the seams between those faces and
/// their edges on an end face are not rims.
#[allow(clippy::too_many_arguments, reason = "the band's reading, passed on")]
fn rims_of(
    model: &Model,
    curved: &Curved,
    skip: &[Id],
    a: &[usize],
    ends: [usize; 2],
    na: Vector,
    room: &Aabb,
    tol: Tolerances,
) -> OgeomResult<Option<(Vec<Rim>, Vec<Corner>)>> {
    let mut rims = Vec::new();
    let mut corners: Vec<Corner> = Vec::new();
    let mut seen: Vec<Id> = Vec::new();
    for (e, _) in a
        .iter()
        .flat_map(|fi| curved.faces[*fi].loops.iter().flatten())
    {
        if skip.contains(e) || seen.contains(e) {
            continue;
        }
        seen.push(*e);
        let between = faces_of(curved, *e);
        if between.iter().any(|f| ends.contains(f)) || between.iter().all(|f| a.contains(f)) {
            continue;
        }
        let Some(cliff) = cliff_of(curved, *e, na, tol.angular()) else {
            let edge = &curved.edges[e];
            if ogeom_algo::shape_bounds(model, &edge.shape, tol)?.intersects(room) {
                return Ok(None);
            }
            continue;
        };
        let edge = &curved.edges[e];
        let (p0, p1) = (curved.points[&edge.ends.0], curved.points[&edge.ends.1]);
        let (cliff, range, whole) = match cliff {
            Cliff::Line { .. } => {
                let len = p0.distance(p1);
                if len <= tol.confusion() {
                    continue;
                }
                (
                    Cliff::Line {
                        at: p0,
                        along: (p1 - p0) / len,
                    },
                    (0.0, len),
                    false,
                )
            }
            Cliff::Circle { .. } => {
                let start = cliff.param(p0);
                let sweep = edge.range.1 - edge.range.0;
                (cliff, (start, start + sweep), edge.ends.0 == edge.ends.1)
            }
        };
        let mut stops = [None, None];
        for (k, v) in [edge.ends.0, edge.ends.1].into_iter().enumerate() {
            if whole || ends.iter().any(|f| touches(curved, *f, v)) {
                continue;
            }
            let at = corners.iter().position(|c| c.id == v).unwrap_or_else(|| {
                corners.push(Corner {
                    id: v,
                    at: curved.points[&v],
                });
                corners.len() - 1
            });
            stops[k] = Some(at);
        }
        rims.push(Rim {
            id: *e,
            cliff,
            range,
            whole,
            stops,
        });
    }
    Ok(Some((rims, corners)))
}

/// `edges` read as one straight chain: lines end to end along one line,
/// its first and last vertices and the points there; `None` for any other.
fn chain_of(curved: &Curved, ids: &[Id], eps: f64) -> Option<(Id, Id)> {
    let mut ends: Vec<(Id, usize)> = Vec::new();
    for e in ids {
        let edge = curved.edges.get(e)?;
        let Curve::Line(_) = &edge.curve else {
            return None;
        };
        if edge.ends.0 == edge.ends.1 {
            return None;
        }
        for v in [edge.ends.0, edge.ends.1] {
            match ends.iter_mut().find(|(w, _)| *w == v) {
                Some(slot) => slot.1 += 1,
                None => ends.push((v, 1)),
            }
        }
    }
    if ends.iter().any(|(_, n)| *n > 2) {
        return None;
    }
    let [first, last] = ends
        .iter()
        .filter(|(_, n)| *n == 1)
        .map(|(v, _)| *v)
        .collect::<Vec<_>>()[..]
    else {
        return None;
    };
    let (p, q) = (curved.points[&first], curved.points[&last]);
    let len = p.distance(q);
    if len <= eps {
        return None;
    }
    let along = (q - p) / len;
    let straight = ends.iter().all(|(v, _)| {
        let w = curved.points[v] - p;
        (w - along * w.dot(along)).magnitude() <= eps
    });
    // One chain, not two pieces: its length is its edges' together.
    let total: f64 = ids
        .iter()
        .map(|e| {
            let edge = &curved.edges[e];
            curved.points[&edge.ends.0].distance(curved.points[&edge.ends.1])
        })
        .sum();
    (straight && (total - len).abs() <= eps * 10.0).then_some((first, last))
}

/// `edges` read as one straight round of `radius` between two planes
/// whose contact on one runs past sharp edges of it, held as hold edge
/// holds them; `None` for any other request, or where anything else of the
/// solid stands within its reach. The edges may be a chain along one line,
/// each between faces of the same two planes.
///
/// # Errors
///
/// As reading the solid.
#[allow(clippy::too_many_lines, reason = "the request read in turn")]
pub(crate) fn reach(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Hold>> {
    let eps = tol.confusion();
    let ang = tol.angular();
    if edges.is_empty() || !radius.is_finite() || radius <= eps {
        return Ok(None);
    }
    for edge in edges {
        let kinds = ogeom_topo::ancestors_of(model, solid, edge, ogeom_topo::ShapeType::Face)?;
        let planes = kinds.iter().all(|f| {
            model
                .node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(|s| matches!(s, SurfaceGeometry::Plane(_)))
        });
        if kinds.len() != 2 || !planes {
            return Ok(None);
        }
        if crate::setback::edge_is_convex(model, solid, edge, radius, tol)? != Some(true) {
            return Ok(None);
        }
    }
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    let ids: Vec<Id> = edges.iter().map(id).collect();
    if ids.iter().any(|e| !curved.edges.contains_key(e)) {
        return Ok(None);
    }
    let Some((v0, v1)) = chain_of(&curved, &ids, eps) else {
        return Ok(None);
    };
    let (p0, p1) = (curved.points[&v0], curved.points[&v1]);
    let len = p0.distance(p1);
    let along = (p1 - p0) / len;
    let [f0, f1] = faces_of(&curved, ids[0])[..] else {
        return Ok(None);
    };
    let normal = |fi: usize| outward(&curved, fi, p0);
    let (Some(n0), Some(n1)) = (normal(f0), normal(f1)) else {
        return Ok(None);
    };
    // Every edge of the chain lies between a face of each plane.
    let mut sides: [Vec<usize>; 2] = [Vec::new(), Vec::new()];
    for e in &ids {
        let [g, h] = faces_of(&curved, *e)[..] else {
            return Ok(None);
        };
        for fi in [g, h] {
            let Some(n) = normal(fi) else {
                return Ok(None);
            };
            if n.dot(n0) >= 1.0 - ang {
                sides[0].push(fi);
            } else if n.dot(n1) >= 1.0 - ang {
                sides[1].push(fi);
            } else {
                return Ok(None);
            }
        }
    }
    for side in &mut sides {
        side.sort_unstable();
        side.dedup();
    }
    let ours: Vec<usize> = sides.concat();
    if sides[0].iter().any(|f| sides[1].contains(f)) {
        return Ok(None);
    }
    // Inside the chain nothing else meets it.
    for e in &ids {
        let edge = &curved.edges[e];
        for v in [edge.ends.0, edge.ends.1] {
            if v != v0 && v != v1 && !others_at(&curved, v, &ours).is_empty() {
                return Ok(None);
            }
        }
    }
    // The end faces: one plane at each end, its normal within 60 degrees of
    // the edge, the way out of the band.
    let mut ends = [0; 2];
    let mut end_normals = [along; 2];
    for (k, v) in [v0, v1].into_iter().enumerate() {
        let [end] = others_at(&curved, v, &ours)[..] else {
            return Ok(None);
        };
        let SurfaceGeometry::Plane(_) = &curved.faces[end].surface else {
            return Ok(None);
        };
        let Some(n) = outward(&curved, end, curved.points[&v]) else {
            return Ok(None);
        };
        let out = if k == 0 { -along } else { along };
        if n.dot(out) < 0.5 {
            return Ok(None);
        }
        ends[k] = end;
        end_normals[k] = n;
    }
    let mut found: Vec<Hold> = Vec::new();
    for (sa, na, nb) in [(0, n0, n1), (1, n1, n0)] {
        let (a_faces, b_faces) = (sides[sa].clone(), sides[1 - sa].clone());
        let across = {
            let v = along.cross(na);
            let k = -radius / (1.0 + na.dot(nb));
            let contact = (na + nb) * k + na * radius;
            if contact.dot(v) < 0.0 { -v } else { v }
        };
        // The room the band's balls take up, with a radius to spare.
        let spare = Vector::new(radius, radius, radius) * 2.0;
        let room = Aabb::of_points(&[p0 - spare, p0 + spare, p1 - spare, p1 + spare]);
        let Some((rims, corners)) = rims_of(model, &curved, &ids, &a_faces, ends, na, &room, tol)?
        else {
            continue;
        };
        let mut hold = Hold {
            curved: curved.clone(),
            edges: ids.clone(),
            a: a_faces[0],
            b: b_faces[0],
            a_faces,
            b_faces,
            ends,
            end_normals,
            origin: p0,
            along,
            len,
            na,
            nb,
            across,
            radius,
            rims,
            corners,
            pieces: Vec::new(),
            at_end: false,
        };
        let Some(pieces) = hold.read_pieces() else {
            continue;
        };
        let (Some(first), Some(last)) = (pieces.first(), pieces.last()) else {
            continue;
        };
        // No seat but the band: nothing held. A run over the whole band
        // would take face `a` whole: the round is then too large for its
        // seat, as the largest-radius readings say, not held.
        if pieces.iter().all(|p| p.seat == Seat::Face)
            || pieces.iter().all(|p| p.seat != Seat::Face)
        {
            continue;
        }
        hold.at_end = first.seat != Seat::Face || last.seat != Seat::Face;
        hold.pieces = pieces;
        found.push(hold);
    }
    let [hold] = &found[..] else {
        return Ok(None);
    };
    // Every edge held is sharp and convex, and so is every edge at a held
    // vertex but those of `a`: the faces past them fall away.
    let mut sharp: Vec<Id> = Vec::new();
    for piece in &hold.pieces {
        match piece.seat {
            Seat::Face => {}
            Seat::Rim(i) => sharp.push(hold.rims[i].id),
            Seat::Corner(k) => {
                let v = hold.corners[k].id;
                for (e, c) in &hold.curved.edges {
                    let in_a = faces_of(&hold.curved, *e)
                        .iter()
                        .any(|f| hold.a_faces.contains(f));
                    if !in_a && (c.ends.0 == v || c.ends.1 == v) {
                        sharp.push(*e);
                    }
                }
            }
        }
    }
    sharp.sort_unstable();
    sharp.dedup();
    for e in sharp {
        let Some(shape) = hold.curved.shapes.get(&e).cloned() else {
            return Ok(None);
        };
        if crate::setback::edge_is_convex(model, solid, &shape, radius, tol)? != Some(true) {
            return Ok(None);
        }
    }
    if !clear(model, hold, tol)? {
        return Ok(None);
    }
    Ok(Some(hold.clone()))
}

/// Whether nothing of the solid but the round's faces, its end faces and
/// the faces past its held edges and vertices stands within the box the
/// tool can take up, and no edge of theirs but the round's, the held
/// edges and those between them.
fn clear(model: &Model, hold: &Hold, tol: Tolerances) -> OgeomResult<bool> {
    let eps = tol.confusion();
    let mut points = vec![hold.on_edge(0.0), hold.on_edge(hold.len)];
    let mut at: Vec<f64> = (0..=64).map(|k| hold.len * f64::from(k) / 64.0).collect();
    at.extend(hold.pieces.iter().map(|p| p.to));
    for x in at {
        let Some((ball, _)) = hold.ball(x, None) else {
            return Ok(false);
        };
        points.extend([ball.on_a, ball.on_b]);
    }
    // An end face leaning over the band carries the round past the edge's
    // end.
    for k in [0, hold.pieces.len()] {
        let Some((pb, pa)) = hold.end_cut(k) else {
            return Ok(false);
        };
        points.extend([pb, pa]);
    }
    let room = Aabb::of_points(&points);
    let margin = hold.radius * 1e-3 + eps * 10.0;
    let (Some(lo), Some(hi)) = (room.low(), room.high()) else {
        return Ok(false);
    };
    let pad = Vector::new(margin, margin, margin);
    let room = Aabb::of_points(&[lo - pad, hi + pad]);
    let c = &hold.curved;
    let mut ours = [
        hold.a_faces.clone(),
        hold.b_faces.clone(),
        hold.ends.to_vec(),
    ]
    .concat();
    let mut held = hold.edges.clone();
    for piece in &hold.pieces {
        match piece.seat {
            Seat::Face => {}
            Seat::Rim(i) => {
                held.push(hold.rims[i].id);
                ours.extend(faces_of(c, hold.rims[i].id));
            }
            Seat::Corner(k) => {
                let v = hold.corners[k].id;
                ours.extend((0..c.faces.len()).filter(|fi| touches(c, *fi, v)));
            }
        }
    }
    ours.sort_unstable();
    ours.dedup();
    for (fi, face) in c.faces.iter().enumerate() {
        if ours.contains(&fi) {
            continue;
        }
        if ogeom_algo::shape_bounds(model, &face.shape, tol)?.intersects(&room) {
            return Ok(false);
        }
    }
    for &fi in &ours {
        for (e, _) in c.faces[fi].loops.iter().flatten() {
            if held.contains(e) {
                continue;
            }
            let between = faces_of(c, *e);
            if between.iter().all(|f| ours.contains(f)) && between.len() == 2 {
                continue;
            }
            if between == [fi] {
                continue;
            }
            if ogeom_algo::shape_bounds(model, &c.edges[e].shape, tol)?.intersects(&room) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// The piece on `seat` over `from..to` fitted through its sections,
/// closer together toward the ends. A piece reaching an end face (`planes`
/// given) is cut by planes turning from the plane at its start to the one
/// at its end, so its end lies in the end face's plane; any other is
/// fitted through its balls' arcs, which at each end are the arcs it
/// shares with its neighbours.
fn fitted(
    hold: &Hold,
    seat: Seat,
    (from, to): (f64, f64),
    planes: Option<((Point, Vector), (Point, Vector))>,
    tol: Tolerances,
) -> OgeomResult<Option<ogeom_geom::BSplineSurface>> {
    let centre = |x: f64| hold.ball(x, Some(seat)).map(|(b, _)| b.centre);
    let (Some(c0), Some(c1)) = (centre(from), centre(to)) else {
        return Ok(None);
    };
    let mut rows = Vec::with_capacity(STATIONS + 1);
    for i in 0..=STATIONS {
        #[expect(clippy::cast_precision_loss, reason = "a hundred stations")]
        let s = 0.5 * (1.0 - (core::f64::consts::PI * i as f64 / STATIONS as f64).cos());
        let x = from + (to - from) * s;
        let row = if let Some(((p0, n0), (p1, n1))) = planes {
            let n = n0 * (1.0 - s) + n1 * s;
            let m = n.magnitude();
            let Some(c) = centre(x) else {
                return Ok(None);
            };
            if m <= 0.0 {
                return Ok(None);
            }
            // The plane through the station's centre, carried to each end's
            // plane as the station reaches it.
            let at = c + (p0 - c0) * (1.0 - s) + (p1 - c1) * s;
            hold.section(seat, at, n / m, x, ACROSS)
        } else {
            (0..ACROSS)
                .map(|k| {
                    #[expect(clippy::cast_precision_loss, reason = "a few dozen points")]
                    let u = k as f64 / (ACROSS - 1) as f64;
                    hold.point(x, u, seat)
                })
                .collect()
        };
        let Some(row) = row else {
            return Ok(None);
        };
        rows.push(row);
    }
    let target = tol.confusion() * 0.5;
    let Ok(fitted) = ogeom_geom::fit::fit_surface_grid_chordal(&rows, 3, target, tol) else {
        return Ok(None);
    };
    Ok((fitted.error <= tol.confusion()).then_some(fitted.curve))
}

/// A curve and the range of it a sheet edge covers.
type Arc = (Curve, (f64, f64));

/// The arc from `ends.0` to `ends.1` of the ellipse the plane through `at`
/// square to `normal` cuts from the band's cylinder, the quarter facing the
/// edge, as an exact rational spline, and whether it runs from `ends.1`
/// to `ends.0` instead.
fn band_end(
    hold: &Hold,
    at: Point,
    normal: Vector,
    ends: (Point, Point),
    tol: Tolerances,
) -> OgeomResult<Option<(bool, Arc)>> {
    let r = hold.radius;
    let axis = hold.band(0.0);
    let slope = hold.along.dot(normal);
    if slope.abs() <= tol.angular() {
        return Ok(None);
    }
    let centre = axis + hold.along * ((at - axis).dot(normal) / slope);
    let minor = hold.along.cross(normal);
    let minor = minor / minor.magnitude();
    let major = normal.cross(minor);
    let ellipse = ogeom_math::Ellipse::new(
        Frame::new(
            centre,
            Direction::new(normal, tol)?,
            Direction::new(major, tol)?,
            tol,
        )?,
        r / slope.abs(),
        r,
        tol,
    )?;
    let frame = ellipse.frame();
    let (x, y) = (frame.x().vector(), frame.y().vector());
    let angle = |p: Point| {
        let q = p - centre;
        (q.dot(y) / r).atan2(q.dot(x) * slope.abs() / r)
    };
    let (from, to) = (angle(ends.0), angle(ends.1));
    // The quarter between the contacts, less than a half turn either way.
    let mut sweep = (to - from).rem_euclid(TAU);
    let mut start = from;
    let reversed = sweep > core::f64::consts::PI;
    if reversed {
        sweep = TAU - sweep;
        start = to;
    }
    let curve = Curve::Ellipse(ogeom_geom::EllipseCurve::new(ellipse));
    let spline = curve.to_bspline_over((start, start + sweep), tol)?;
    let range = ogeom_geom::Curve3d::domain(&spline);
    let curve = Curve::BSpline(spline);
    let (p, q) = (curve.point_at(range.0, tol)?, curve.point_at(range.1, tol)?);
    let (p0, p1) = if reversed { (ends.1, ends.0) } else { ends };
    let near = tol.confusion() * 10.0;
    Ok((p.distance(p0) <= near && q.distance(p1) <= near).then_some((reversed, (curve, range))))
}

/// The tool: the region between the edge and the round, closed by its
/// two faces' planes and its end faces' planes, with its round faces.
#[allow(
    clippy::too_many_lines,
    reason = "one closed tool, its pieces walked in turn"
)]
fn tool(
    model: &mut Model,
    hold: &Hold,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
    let r = hold.radius;
    let pieces = &hold.pieces;
    let last = pieces.len();
    let mut breaks = vec![0.0];
    breaks.extend(pieces.iter().map(|p| p.to));
    let mut planes = Vec::with_capacity(breaks.len());
    for (k, &x) in breaks.iter().enumerate() {
        let Some(plane) = hold.plane_at(x, k) else {
            return Ok(None);
        };
        planes.push(plane);
    }
    // Each piece off the band, fitted once.
    let mut surfaces = Vec::with_capacity(last);
    for (k, piece) in pieces.iter().enumerate() {
        surfaces.push(if piece.seat == Seat::Face {
            None
        } else {
            // A vertex piece's balls all pass through the vertex: it is
            // fitted through their arcs, and not built where an end face
            // would cut it.
            let reaches_end = k == 0 || k + 1 == last;
            if reaches_end && matches!(piece.seat, Seat::Corner(_)) {
                return Ok(None);
            }
            let planes = reaches_end.then(|| (planes[k], planes[k + 1]));
            let Some(s) = fitted(hold, piece.seat, (piece.from, piece.to), planes, tol)? else {
                return Ok(None);
            };
            Some(s)
        });
    }
    let mut points = vec![hold.on_edge(0.0), hold.on_edge(hold.len)];
    let mut corner_points: Vec<(usize, usize)> = Vec::new();
    // At each break, its contact on `b` and on `a` or the held edge.
    let mut at_b = Vec::with_capacity(breaks.len());
    let mut at_a = Vec::with_capacity(breaks.len());
    for (k, &x) in breaks.iter().enumerate() {
        let piece = if k == last { pieces[k - 1] } else { pieces[k] };
        let end = k == 0 || k == last;
        let (pb, pa) = if end {
            let Some(cut) = hold.end_cut(k) else {
                return Ok(None);
            };
            cut
        } else {
            let Some((ball, _)) = hold.ball(x, Some(piece.seat)) else {
                return Ok(None);
            };
            (ball.on_b, ball.on_a)
        };
        at_b.push(points.len());
        points.push(pb);
        // A held vertex is one point for every piece reaching it.
        let corner = [
            k.checked_sub(1).map(|j| pieces[j].seat),
            pieces.get(k).map(|p| p.seat),
        ]
        .into_iter()
        .flatten()
        .find_map(|s| match s {
            Seat::Corner(c) => Some(c),
            _ => None,
        });
        if let Some(c) = corner {
            if let Some(&(_, i)) = corner_points.iter().find(|(cc, _)| *cc == c) {
                at_a.push(i);
            } else {
                corner_points.push((c, points.len()));
                at_a.push(points.len());
                points.push(hold.corners[c].at);
            }
        } else {
            at_a.push(points.len());
            points.push(pa);
        }
    }
    let (e0, e1) = (0, 1);
    let mut curves: Vec<SheetCurve> = Vec::new();
    let mut push = |ends: (usize, usize), (curve, range): (Curve, (f64, f64))| {
        curves.push(SheetCurve { ends, curve, range });
        curves.len() - 1
    };
    let edge_line = push((e0, e1), segment(points[e0], points[e1], tol)?);
    let b_start = push((e0, at_b[0]), segment(points[e0], points[at_b[0]], tol)?);
    let a_start = push((e0, at_a[0]), segment(points[e0], points[at_a[0]], tol)?);
    let b_end = push(
        (e1, at_b[last]),
        segment(points[e1], points[at_b[last]], tol)?,
    );
    let a_end = push(
        (e1, at_a[last]),
        segment(points[e1], points[at_a[last]], tol)?,
    );
    // The section at each break: the band's arc, the arc of the ball two
    // pieces share, or a fitted piece's end in an end face's plane. Where a
    // piece on a vertex meets the next it is that piece's own end, which
    // its surface holds exactly up to the vertex it closes to.
    let mut sections = Vec::with_capacity(breaks.len());
    for (k, &x) in breaks.iter().enumerate() {
        let ends = (at_b[k], at_a[k]);
        let piece = if k == last { k - 1 } else { k };
        let on_corner = [k.checked_sub(1), (k < last).then_some(k)]
            .into_iter()
            .flatten()
            .find(|j| matches!(pieces[*j].seat, Seat::Corner(_)));
        let own = if k == 0 || k == last {
            Some(piece)
        } else {
            on_corner
        };
        if let Some(j) = own
            && let Some(surface) = &surfaces[j]
        {
            let Some(side) = boundary_from(surface, &points, ends, tol)? else {
                return Ok(None);
            };
            sections.push(push(side.ends, (side.curve, side.range)));
        } else {
            let seat = if k == 0 {
                pieces[0].seat
            } else {
                pieces[k - 1].seat
            };
            let (at, normal) = planes[k];
            let oblique =
                (k == 0 || k == last) && normal.cross(hold.along).magnitude() > tol.angular();
            if oblique {
                // The band cut by an end face leaning over it: an arc of the
                // ellipse its plane cuts from the band's cylinder.
                let Some((flip, cut)) =
                    band_end(hold, at, normal, (points[at_b[k]], points[at_a[k]]), tol)?
                else {
                    return Ok(None);
                };
                sections.push(push(if flip { (ends.1, ends.0) } else { ends }, cut));
                continue;
            }
            let Some((ball, _)) = hold.ball(x, Some(seat)) else {
                return Ok(None);
            };
            let (u, w) = (ball.on_b - ball.centre, points[at_a[k]] - ball.centre);
            let z = u.cross(w);
            let sweep = (u.dot(w) / (r * r)).clamp(-1.0, 1.0).acos();
            sections.push(push(ends, arc(ball.centre, z, u, r, sweep, tol)?));
        }
    }
    // Along each piece, its contact on `b` and on `a` or the held edge; a
    // piece on a vertex has its one foot there.
    let mut along_b = Vec::with_capacity(last);
    let mut along_a = Vec::with_capacity(last);
    for (k, surface) in surfaces.iter().enumerate() {
        let (b0, b1, a0, a1) = (at_b[k], at_b[k + 1], at_a[k], at_a[k + 1]);
        if let Some(surface) = surface {
            let Some(sb) = boundary_from(surface, &points, (b0, b1), tol)? else {
                return Ok(None);
            };
            along_b.push(push(sb.ends, (sb.curve, sb.range)));
            if a0 == a1 {
                along_a.push(None);
            } else {
                let Some(sa) = boundary_from(surface, &points, (a0, a1), tol)? else {
                    return Ok(None);
                };
                along_a.push(Some(push(sa.ends, (sa.curve, sa.range))));
            }
        } else {
            along_b.push(push((b0, b1), segment(points[b0], points[b1], tol)?));
            along_a.push(Some(push((a0, a1), segment(points[a0], points[a1], tol)?)));
        }
    }
    let first = pieces[0];
    let mid = 0.5 * (first.from + first.to);
    let Some((ball, _)) = hold.ball(mid, Some(first.seat)) else {
        return Ok(None);
    };
    let (mb, ma) = (ball.on_b, ball.on_a);
    let on_mid = hold.on_edge(mid);
    let curved = &hold.curved;
    let mut faces = vec![
        ToolFace {
            surface: curved.faces[hold.b].surface.clone(),
            curves: [vec![edge_line, b_end], along_b.clone(), vec![b_start]].concat(),
            at: on_mid + (mb - on_mid) * 0.5,
            out: hold.nb,
        },
        ToolFace {
            surface: curved.faces[hold.a].surface.clone(),
            curves: [
                vec![edge_line, a_end],
                along_a.iter().flatten().copied().collect(),
                vec![a_start],
            ]
            .concat(),
            at: on_mid + (ma - on_mid) * 0.5,
            out: hold.na,
        },
    ];
    for (k, end) in [(0, hold.ends[0]), (last, hold.ends[1])] {
        let (pe, out) = if k == 0 {
            (points[e0], hold.end_normals[0])
        } else {
            (points[e1], hold.end_normals[1])
        };
        let (pb, pa) = (points[at_b[k]], points[at_a[k]]);
        faces.push(ToolFace {
            surface: curved.faces[end].surface.clone(),
            curves: if k == 0 {
                vec![b_start, sections[0], a_start]
            } else {
                vec![b_end, sections[last], a_end]
            },
            at: pe + (pb - pe) * 0.2 + (pa - pe) * 0.2,
            out,
        });
    }
    let round_from = faces.len();
    for (k, piece) in pieces.iter().enumerate() {
        let x = 0.5 * (piece.from + piece.to);
        let Some((ball, _)) = hold.ball(x, Some(piece.seat)) else {
            return Ok(None);
        };
        let Some(p) = hold.point(x, 0.5, piece.seat) else {
            return Ok(None);
        };
        let surface = if let Some(s) = &surfaces[k] {
            s.clone().into()
        } else {
            let frame = Frame::new(
                hold.band(0.0),
                Direction::new(hold.along, tol)?,
                Direction::new(hold.nb, tol)?,
                tol,
            )?;
            CylinderSurface::new(
                Cylinder::new(frame, r, tol)?,
                (-4.0 * r, hold.len + 4.0 * r),
            )?
            .into()
        };
        let mut loop_ = vec![along_b[k], sections[k + 1]];
        loop_.extend(along_a[k]);
        loop_.push(sections[k]);
        faces.push(ToolFace {
            surface,
            curves: loop_,
            at: p,
            out: (ball.centre - p) / r,
        });
    }
    let Some((solid, made)) = closed(model, &points, &curves, &faces, tol)? else {
        return Ok(None);
    };
    Ok(Some((solid, made[round_from..].to_vec())))
}

/// The round of `hold` cut from `solid` as one tool; `None` where the tool
/// will not build, apply or prove.
///
/// # Errors
///
/// As the booleans, other than their refusals.
pub(crate) fn held(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    hold: &Hold,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some((cut, rounds)) = tool(model, hold, tol)? else {
        return Ok(None);
    };
    let deflection = ogeom_mesh::Deflection::default();
    let before = ogeom_algo::volume_properties(model, solid, deflection, tol)?.mass;
    let taken = ogeom_algo::volume_properties(model, &cut, deflection, tol)?.mass;
    let named = rounds
        .iter()
        .flat_map(|f| edges.iter().map(move |e| (f.clone(), e.clone())))
        .collect();
    let tools = [super::fallback::Tool {
        solid: cut,
        sense: super::fallback::Sense::Cut,
        faces: named,
    }];
    let Some(built) = super::fallback::applied(model, solid, edges, &tools, tol)? else {
        return Ok(None);
    };
    let built = super::merged(model, solid, built, tol)?;
    // The tool stands inside the solid: the cut takes its whole volume.
    let after = ogeom_algo::volume_properties(model, &built.shape, deflection, tol)?.mass;
    if (before - after - taken).abs() > 1e-6 * (1.0 + before.abs()) {
        return Ok(None);
    }
    if !ogeom_algo::check(model, &built.shape, tol)?.is_valid() {
        return Ok(None);
    }
    if !kept(model, &built.shape, hold, tol)? {
        return Ok(None);
    }
    Ok(Some(built))
}

/// Whether the result's edges lying on the held edges' curves still hold
/// the balls' feet over every piece resting on one.
fn kept(model: &Model, shape: &Shape, hold: &Hold, tol: Tolerances) -> OgeomResult<bool> {
    use ogeom_topo::{ShapeType, explore_unique};
    let near = tol.confusion() * 100.0;
    let mut runs: Vec<(usize, Vec<Point>)> = Vec::new();
    for e in explore_unique(model, shape, ShapeType::Edge)? {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(model, &e, tol)? else {
            continue;
        };
        let mut samples = Vec::with_capacity(SAMPLES + 1);
        for k in 0..=SAMPLES {
            #[expect(clippy::cast_precision_loss, reason = "a few dozen samples")]
            let p = c.point_at(a + (b - a) * k as f64 / SAMPLES as f64, tol)?;
            samples.push(p);
        }
        for (i, rim) in hold.rims.iter().enumerate() {
            let on = samples.iter().all(|p| {
                (*p - hold.origin).dot(hold.na).abs() <= near
                    && rim.cliff.foot(*p).is_some_and(|f| f.distance(*p) <= near)
            });
            if on {
                runs.push((i, samples.clone()));
            }
        }
    }
    Ok(hold.feet_on(&runs))
}
