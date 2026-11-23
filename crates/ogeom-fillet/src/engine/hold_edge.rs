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
//! a plane `a` and a plane or a cylinder `b` whose axis runs along the
//! edge, the held edges any curves lying in `a`'s plane, the round ending
//! on end faces within 60 degrees of square to its edge. On a cylinder the
//! band's centres still lie on a line, and
//! off the band a ball's centre rises round the axis on `b`'s offset. Each
//! ball's arc runs from its contact on `b` to its foot, on its
//! characteristic circle. Where the seat changes the ball touches both
//! seats at once and the two offsets are tangent, so neighbouring pieces
//! share that ball's characteristic circle and meet G1 on it. From the band
//! to an edge the spine runs along the edge there, so that circle is the
//! band's own section; from an edge to a vertex, or to the next edge where
//! the two meet tangentially and the foot passes through their vertex, it
//! is the ball's arc in its plane square to the spine. A piece is fitted through its balls' arcs,
//! or, where it reaches an end face, through its sections by planes turning
//! to that face's plane about the line the two planes share. Past an end face the centres run on as they are
//! and the end face's plane cuts the round (amendment "M5f"). Where that
//! is not one piece's section (the balls past the face that still reach
//! inside it change seat, the piece there rests on a vertex, or the face
//! leans on a cylinder), the round runs on to where it stands wholly past
//! the face and the cut trims it there. Where the
//! nearest point jumps from one edge to another, the ball resting on both
//! at once, its centre's path turns a corner: the round there is that
//! ball's own surface between its two arcs, up to `a`'s plane, which meets
//! both pieces G1, and the faces past the two edges trim the round where
//! it passes over them.
//!
//! The region between the edge and the round is one closed tool, cut from
//! the solid; it is built only where nothing else of the solid stands
//! within its reach, the cut takes the tool's whole volume (where the ball
//! jumps between edges, its volume within the solid), the result is valid
//! and the held edges are still there under every ball's foot.

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

/// The curve of a sharp edge of `a`: a line through `at` along `along`,
/// the circle of `radius` about `centre` in the plane square to `normal`,
/// its angle measured from `x` toward `y`, or any other curve in `a`'s
/// plane over the range the edge covers, read at `samples` for the nearest
/// point's start, `closed` where the edge is the whole curve. An open
/// curve of the third kind is carried on past each end along its tangent
/// there, its parameter growing at the rate it has at that end.
#[derive(Debug, Clone)]
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
    General {
        curve: Curve,
        range: (f64, f64),
        samples: Vec<(f64, Point)>,
        closed: bool,
        tol: Tolerances,
    },
}

/// Points read along a curve of the third kind for its nearest point.
const READ_ALONG: usize = 256;

impl Cliff {
    /// The point of the curve nearest `p`, `None` on a circle's axis.
    fn foot(&self, p: Point) -> Option<Point> {
        match *self {
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
            Self::General { .. } => self.nearest(p).map(|(_, f)| f),
        }
    }

    /// Where `p` lies along the curve: along a line, its angle about a
    /// circle, or the parameter of its nearest point.
    fn param(&self, p: Point) -> f64 {
        match *self {
            Self::Line { at, along } => (p - at).dot(along),
            Self::Circle { centre, x, y, .. } => {
                let q = p - centre;
                q.dot(y).atan2(q.dot(x))
            }
            Self::General { range, .. } => self.nearest(p).map_or(range.0, |(t, _)| t),
        }
    }

    /// The point at `t` of a curve of the third kind, its first and second
    /// derivatives.
    fn general_at(&self, t: f64) -> Option<(Point, Vector, Vector)> {
        let Self::General { curve, tol, .. } = self else {
            return None;
        };
        let d = curve.derivatives_at(t, 2, *tol).ok()?;
        let [p, d1, d2] = d[..] else {
            return None;
        };
        Some((Point::new(p.x, p.y, p.z), d1, d2))
    }

    /// The end `k` of an open curve of the third kind: its point, its unit
    /// tangent toward growing parameter and its parameter's rate there.
    fn general_end(&self, k: usize) -> Option<(Point, Vector, f64)> {
        let Self::General { range, .. } = self else {
            return None;
        };
        let t = if k == 0 { range.0 } else { range.1 };
        let (p, d1, _) = self.general_at(t)?;
        let m = d1.magnitude();
        (m > 0.0).then(|| (p, d1 / m, m))
    }

    /// The parameter and point of a curve of the third kind nearest `p`:
    /// from the nearest of its samples, by Newton's step kept within the
    /// samples either side; past an open end, on the tangent there.
    fn nearest(&self, p: Point) -> Option<(f64, Point)> {
        let Self::General {
            range,
            samples,
            closed,
            ..
        } = self
        else {
            return None;
        };
        let k = (0..samples.len()).min_by(|i, j| {
            p.distance(samples[*i].1)
                .total_cmp(&p.distance(samples[*j].1))
        })?;
        let n = samples.len() - 1;
        let (mut lo, mut hi) = (samples[k.saturating_sub(1)].0, samples[(k + 1).min(n)].0);
        let g = |t: f64| {
            self.general_at(t)
                .map(|(c, d1, d2)| ((c - p).dot(d1), d1.dot(d1) + (c - p).dot(d2)))
        };
        let (glo, ghi) = (g(lo)?.0, g(hi)?.0);
        let mut t = if glo >= 0.0 {
            lo
        } else if ghi <= 0.0 {
            hi
        } else {
            let span = range.1 - range.0;
            let mut t = samples[k].0;
            for _ in 0..60 {
                let (v, slope) = g(t)?;
                if v == 0.0 {
                    break;
                }
                if v > 0.0 {
                    hi = t;
                } else {
                    lo = t;
                }
                let next = if slope > 0.0 { t - v / slope } else { f64::NAN };
                let next = if next > lo && next < hi {
                    next
                } else {
                    0.5 * (lo + hi)
                };
                let done = (next - t).abs() <= 1e-15 * span || hi - lo <= 1e-15 * span;
                t = next;
                if done {
                    break;
                }
            }
            t
        };
        if !closed {
            for (end, sign) in [(0, -1.0), (1, 1.0)] {
                let at_end = if end == 0 { range.0 } else { range.1 };
                if t != at_end {
                    continue;
                }
                let (q, tangent, rate) = self.general_end(end)?;
                let s = (p - q).dot(tangent);
                if s * sign > 0.0 {
                    return Some((at_end + s / rate, q + tangent * s));
                }
            }
        }
        t = t.clamp(range.0, range.1);
        Some((t, self.general_at(t)?.0))
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
            Cliff::General { closed: true, .. } => Along::Inside,
            Cliff::Line { .. } | Cliff::General { .. } => {
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

    /// The unit tangent of the rim's curve at `p` on it, `None` on a
    /// circle's axis.
    fn tangent(&self, p: Point) -> Option<Vector> {
        match self.cliff {
            Cliff::Line { along, .. } => Some(along),
            Cliff::Circle { centre, normal, .. } => {
                let t = normal.cross(p - centre);
                let m = t.magnitude();
                (m > 0.0).then(|| t / m)
            }
            Cliff::General { range, .. } => {
                let t = self.cliff.param(p);
                if t <= range.0 {
                    return self.cliff.general_end(0).map(|(_, d, _)| d);
                }
                if t >= range.1 {
                    return self.cliff.general_end(1).map(|(_, d, _)| d);
                }
                let (_, d1, _) = self.cliff.general_at(t)?;
                let m = d1.magnitude();
                (m > 0.0).then(|| d1 / m)
            }
        }
    }

    /// Whether `t` lies within the edge itself, not past an end.
    fn covers(&self, t: f64, slack: f64) -> bool {
        match &self.cliff {
            Cliff::Line { .. } => t >= self.range.0 - slack && t <= self.range.1 + slack,
            Cliff::General {
                closed, samples, ..
            } => {
                // The slack, a length, as a parameter span.
                let length: f64 = samples.windows(2).map(|w| w[0].1.distance(w[1].1)).sum();
                let slack = slack * (self.range.1 - self.range.0) / length.max(f64::MIN_POSITIVE);
                *closed || (t >= self.range.0 - slack && t <= self.range.1 + slack)
            }
            Cliff::Circle { .. } => {
                self.whole
                    || (t - self.range.0 + slack).rem_euclid(TAU)
                        <= self.range.1 - self.range.0 + 2.0 * slack
            }
        }
    }
}

/// A vertex where two sharp edges of `a` meet, and whether they meet
/// tangentially there, so the ball's foot passes from one to the other
/// without resting on the vertex.
#[derive(Debug, Clone, Copy)]
struct Corner {
    id: Id,
    at: Point,
    smooth: bool,
}

/// A stretch of the edge over which the ball keeps one seat, and whether
/// the foot jumps onto it from the previous piece's edge.
#[derive(Debug, Clone, Copy)]
struct Piece {
    seat: Seat,
    from: f64,
    to: f64,
    jump: bool,
}

/// The ball where the foot jumps from one edge to another: its centre, its
/// contact on `b`, its feet on the two edges and the spine's unit tangent
/// on each side, and the arc in `a`'s plane between the feet where the
/// ball's surface meets that plane, about `hub` square to `normal`, from
/// the direction `from` through `sweep` at radius `reach`.
#[derive(Debug, Clone, Copy)]
struct Jump {
    centre: Point,
    on_b: Point,
    feet: [Point; 2],
    spines: [Vector; 2],
    hub: Point,
    normal: Vector,
    from: Vector,
    reach: f64,
    sweep: f64,
}

/// The face `b` the ball touches throughout: a plane, or a cylinder whose
/// axis runs along the edge, through `axis`, of `radius`, `sense` 1 where
/// the solid lies inside it and -1 where it lies outside. Either way the
/// band's centres lie on a line along the edge.
#[derive(Debug, Clone, Copy)]
enum Front {
    Plane,
    Cylinder {
        axis: Point,
        radius: f64,
        sense: f64,
    },
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
    /// The outward normals of `a` and of `b` at the edge, and `b` itself.
    na: Vector,
    nb: Vector,
    front: Front,
    /// The band's centre at the edge's start.
    centre: Point,
    /// In `a`, square to the edge, pointing from the edge into `a`.
    across: Vector,
    radius: f64,
    rims: Vec<Rim>,
    corners: Vec<Corner>,
    pieces: Vec<Piece>,
    /// Whether a stretch held reaches an end of the band.
    pub(crate) at_end: bool,
    /// At each end run on past its end face, that end face's point on the
    /// edge and outward normal: the round there ends wholly past the end
    /// face, which the cut trims it by.
    ends_at: [Option<(Point, Vector)>; 2],
    /// The stretch of the spine along the edge itself.
    real: (f64, f64),
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
        self.centre + self.along * x
    }

    /// The outward normal of `b` at the point of it nearest `p`.
    fn nb_at(&self, p: Point) -> Vector {
        match self.front {
            Front::Plane => self.nb,
            Front::Cylinder { axis, sense, .. } => {
                let q = p - axis;
                let radial = q - self.along * q.dot(self.along);
                radial * (sense / radial.magnitude())
            }
        }
    }

    /// The curve on `b` from `from` to `to`, both on it in one plane square
    /// to the edge: a segment on a plane, an arc about a cylinder's axis.
    fn across_b(&self, from: Point, to: Point, tol: Tolerances) -> OgeomResult<Arc> {
        match self.front {
            Front::Plane => segment(from, to, tol),
            Front::Cylinder { axis, radius, .. } => {
                let hub = axis + self.along * (from - axis).dot(self.along);
                let (u, w) = (from - hub, to - hub);
                let z = u.cross(w);
                let sweep = (u.dot(w) / (u.magnitude() * w.magnitude()))
                    .clamp(-1.0, 1.0)
                    .acos();
                arc(hub, z, u, radius, sweep, tol)
            }
        }
    }

    /// Points along the curve on `b` from `from` to `to`, as
    /// [`Self::across_b`] runs.
    fn across_b_points(&self, from: Point, to: Point) -> Vec<Point> {
        match self.front {
            Front::Plane => vec![from, to],
            Front::Cylinder { axis, radius, .. } => {
                let hub = axis + self.along * (from - axis).dot(self.along);
                (0..=8)
                    .map(|k| {
                        let s = f64::from(k) / 8.0;
                        let v = (from - hub) * (1.0 - s) + (to - hub) * s;
                        hub + v * (radius / v.magnitude())
                    })
                    .collect()
            }
        }
    }

    /// The point `s` along the path of the centres touching `b` that rises
    /// from the band's centre `c` toward `a`'s plane, and how far that path
    /// runs before it reaches the plane or turns back from it: a line on a
    /// plane's offset, an arc about a cylinder's axis.
    fn rise(&self, c: Point) -> Option<(impl Fn(f64) -> Point, f64)> {
        let (along, na) = (self.along, self.na);
        let height = move |p: Point| (p - self.origin).dot(na);
        let (step, foot, v0, w0) = match self.front {
            Front::Plane => {
                let g = along.cross(self.nb);
                let g = if g.dot(na) < 0.0 { -g } else { g };
                (Some(g), c, Vector::new(0.0, 0.0, 0.0), g)
            }
            Front::Cylinder { axis, .. } => {
                let foot = axis + along * (c - axis).dot(along);
                let v0 = c - foot;
                let w0 = along.cross(v0);
                let w0 = if w0.dot(na) < 0.0 { -w0 } else { w0 };
                (None, foot, v0, w0)
            }
        };
        let at = move |s: f64| match step {
            Some(g) => foot + g * s,
            None => foot + v0 * s.cos() + w0 * s.sin(),
        };
        let reach = match step {
            Some(g) => -height(c) / g.dot(na),
            None => {
                // The arc's highest point, or where it reaches the plane
                // before that.
                let top = w0.dot(na).atan2(v0.dot(na));
                if top <= 0.0 {
                    return None;
                }
                if height(at(top)) < 0.0 {
                    top
                } else {
                    root(|s| Some(height(at(s))), 0.0, top)?
                }
            }
        };
        (reach > 0.0).then_some((at, reach))
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
        let far = 10.0 * (self.len + self.radius);
        // Whether the walk crosses the segment from `p` to `q`.
        let crosses = |p: Point, q: Point| {
            let (px, qx) = (
                (p - self.origin).dot(self.along) - x,
                (q - self.origin).dot(self.along) - x,
            );
            if px * qx >= 0.0 {
                return false;
            }
            let (py, qy) = (
                (p - self.origin).dot(self.across),
                (q - self.origin).dot(self.across),
            );
            let cy = py + (qy - py) * px / (px - qx);
            cy > start && cy < y
        };
        let mut crossings = 0;
        for rim in &self.rims {
            match &rim.cliff {
                &Cliff::Line { at, along } => {
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
                    if crosses(at + along * lo, at + along * hi) {
                        crossings += 1;
                    }
                }
                Cliff::General {
                    samples, closed, ..
                } => {
                    // Where the curve meets the plane of the walk, between
                    // samples either side of it, and on the tangents an
                    // open end is carried on along.
                    let gap = |t: f64| {
                        rim.cliff
                            .general_at(t)
                            .map(|(p, _, _)| (p - base).dot(self.along))
                    };
                    for w in samples.windows(2) {
                        let (a, b) = (
                            (w[0].1 - base).dot(self.along),
                            (w[1].1 - base).dot(self.along),
                        );
                        if (a <= 0.0) == (b <= 0.0) {
                            continue;
                        }
                        let Some(t) = root(gap, w[0].0, w[1].0) else {
                            continue;
                        };
                        let Some((p, _, _)) = rim.cliff.general_at(t) else {
                            continue;
                        };
                        let cy = (p - self.origin).dot(self.across);
                        if cy > start && cy < y {
                            crossings += 1;
                        }
                    }
                    if !closed {
                        for (k, sign) in [(0, -1.0), (1, 1.0)] {
                            if rim.stops[k].is_some() {
                                continue;
                            }
                            let Some((q, tangent, _)) = rim.cliff.general_end(k) else {
                                continue;
                            };
                            if crosses(q, q + tangent * (sign * far)) {
                                crossings += 1;
                            }
                        }
                    }
                }
                &Cliff::Circle { centre, radius, .. } => {
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
                        // Past a tangent vertex the next edge's own foot is
                        // as near.
                        Along::Stopped(k) if self.corners[k].smooth => continue,
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
            on_b: c + self.nb_at(c) * r,
            on_a: c + self.na * r,
        };
        if seat == Some(Seat::Face) || (seat.is_none() && self.inside_at(x, self.contact_depth())) {
            return Some((band, Seat::Face));
        }
        let (path, reach) = self.rise(c)?;
        let f = |s: f64| self.support(path(s), seat).map(|(d, _, _)| d - r);
        let s = root(f, 0.0, reach)?;
        let centre = path(s);
        let (_, foot, found) = self.support(centre, seat)?;
        Some((
            Ball {
                centre,
                on_b: centre + self.nb_at(centre) * r,
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
    /// the round meets G1 across: the band onto an edge it crosses, an edge
    /// onto a vertex where it stops, or an edge onto the next where the two
    /// meet tangentially, the foot passing through their vertex, or an
    /// edge onto another the foot jumps to, the ball resting on both, its
    /// own sphere meeting the pieces either side G1. `None` where the seat
    /// cannot be read, or changes any other way.
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
                jump: false,
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
                    jump: false,
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
            jump: false,
        });
        let mut fits = true;
        for k in 1..pieces.len() {
            let (before, x) = (pieces[k - 1].seat, pieces[k - 1].to);
            fits &= match (before, pieces[k].seat) {
                // The band's contact crosses the edge there.
                (Seat::Face, Seat::Rim(i)) | (Seat::Rim(i), Seat::Face) => {
                    let contact = self.band(x) + self.na * self.radius;
                    self.rims[i]
                        .cliff
                        .foot(contact)
                        .is_some_and(|f| f.distance(contact) <= self.radius * 1e-6)
                }
                (Seat::Rim(i), Seat::Corner(c)) | (Seat::Corner(c), Seat::Rim(i)) => {
                    self.rims[i].stops.contains(&Some(c))
                }
                (Seat::Rim(i), Seat::Rim(j)) => {
                    let jump = !self.passes(i, j, x) && self.jump(x, i, j).is_some();
                    pieces[k].jump = jump;
                    jump || self.passes(i, j, x)
                }
                _ => false,
            };
        }
        let short = pieces.iter().any(|p| p.to - p.from <= self.radius * 1e-4);
        (fits && !short).then_some(pieces)
    }

    /// Whether the foot of the ball at `x` passes from rim `i` onto rim `j`
    /// through a tangent vertex where both stop: the foot on each lies at
    /// that vertex. Elsewhere the nearest point jumps between them.
    fn passes(&self, i: usize, j: usize, x: f64) -> bool {
        let near = self.radius * 1e-6;
        let [ri, rj] = [&self.rims[i], &self.rims[j]];
        ri.stops.iter().flatten().any(|&k| {
            let corner = self.corners[k];
            corner.smooth
                && rj.stops.contains(&Some(k))
                && [Seat::Rim(i), Seat::Rim(j)].into_iter().all(|seat| {
                    self.ball(x, Some(seat))
                        .is_some_and(|(b, _)| b.on_a.distance(corner.at) <= near)
                })
        })
    }

    /// The unit tangent of the path of the centres on `seat` at `x`.
    fn spine(&self, x: f64, seat: Seat) -> Option<Vector> {
        let h = self.radius * 1e-5;
        let (ahead, _) = self.ball(x + h, Some(seat))?;
        let (behind, _) = self.ball(x - h, Some(seat))?;
        let t = ahead.centre - behind.centre;
        let m = t.magnitude();
        (m > 0.0).then(|| t / m)
    }

    /// The ball at `x` where the foot jumps from rim `i` to rim `j`: one
    /// centre the radius from both, two feet apart, the path of the centres
    /// turning there, and each foot clear of the other side's next balls,
    /// so the ball's surface between its two arcs (the lune facing the way
    /// the path turns from) is the round there. `None` for any other change.
    fn jump(&self, x: f64, i: usize, j: usize) -> Option<Jump> {
        let r = self.radius;
        let near = r * 1e-6;
        let (bi, _) = self.ball(x, Some(Seat::Rim(i)))?;
        let (bj, _) = self.ball(x, Some(Seat::Rim(j)))?;
        let (ti, tj) = (self.spine(x, Seat::Rim(i))?, self.spine(x, Seat::Rim(j))?);
        let c = bi.centre;
        let (ui, uj) = ((bi.on_a - c) / r, (bj.on_a - c) / r);
        let apart = c.distance(bj.centre) <= near
            && bi.on_a.distance(bj.on_a) > near
            && ti.cross(tj).magnitude() > 1e-6;
        if !apart || ui.dot(tj) > 1e-7 || uj.dot(ti) < -1e-7 {
            return None;
        }
        let hub = c - self.na * (c - self.origin).dot(self.na);
        let (fi, fj) = (bi.on_a - hub, bj.on_a - hub);
        let reach = fi.magnitude();
        if reach <= near || (fj.magnitude() - reach).abs() > near {
            return None;
        }
        let from = fi / reach;
        let turn = fj
            .dot(self.na.cross(from))
            .atan2(fj.dot(from))
            .rem_euclid(TAU);
        // The way round from one foot to the other that stays on the lune.
        [(self.na, turn), (-self.na, TAU - turn)]
            .into_iter()
            .find(|(normal, sweep)| {
                let half = 0.5 * sweep;
                let mid = hub + (from * half.cos() + normal.cross(from) * half.sin()) * reach;
                (mid - c).dot(ti) >= -near && (mid - c).dot(tj) <= near
            })
            .map(|(normal, sweep)| Jump {
                centre: c,
                on_b: bi.on_b,
                feet: [bi.on_a, bj.on_a],
                spines: [ti, tj],
                hub,
                normal,
                from,
                reach,
                sweep,
            })
    }

    /// The jump at the start of piece `k`, where it has one.
    fn jump_into(&self, k: usize) -> Option<Option<Jump>> {
        let piece = self.pieces[k];
        if !piece.jump {
            return Some(None);
        }
        let (Seat::Rim(i), Seat::Rim(j)) = (self.pieces[k - 1].seat, piece.seat) else {
            return None;
        };
        self.jump(piece.from, i, j).map(Some)
    }

    /// The one point every ball of piece `k` sets its foot on, where there
    /// is one: a vertex seat, or an edge square to the spine.
    fn pin(&self, k: usize) -> Option<Point> {
        let piece = self.pieces[k];
        match piece.seat {
            Seat::Face => None,
            Seat::Corner(c) => Some(self.corners[c].at),
            Seat::Rim(_) => {
                let foot = |x: f64| self.ball(x, Some(piece.seat)).map(|(b, _)| b.on_a);
                let mid = 0.5 * (piece.from + piece.to);
                let (p, q, m) = (foot(piece.from)?, foot(piece.to)?, foot(mid)?);
                let near = self.radius * 1e-9;
                (p.distance(q) <= near && p.distance(m) <= near).then_some(m)
            }
        }
    }

    /// The rims from `i` to `j` through the vertices where they stop, and
    /// those vertices, by the fewest rims; `None` where no walk joins them.
    fn walk(&self, i: usize, j: usize) -> Option<(Vec<usize>, Vec<usize>)> {
        let mut came: Vec<Option<(usize, usize)>> = vec![None; self.rims.len()];
        let mut queue = std::collections::VecDeque::from([i]);
        let mut seen = vec![false; self.rims.len()];
        seen[i] = true;
        while let Some(a) = queue.pop_front() {
            if a == j {
                let (mut rims, mut corners) = (vec![j], Vec::new());
                let mut at = j;
                while let Some((prev, corner)) = came[at] {
                    rims.push(prev);
                    corners.push(corner);
                    at = prev;
                }
                return Some((rims, corners));
            }
            for &c in self.rims[a].stops.iter().flatten() {
                for (b, rim) in self.rims.iter().enumerate() {
                    if !seen[b] && rim.stops.contains(&Some(c)) {
                        seen[b] = true;
                        came[b] = Some((a, c));
                        queue.push_back(b);
                    }
                }
            }
        }
        None
    }

    /// The plane a piece on `seat` ends in at break `k`, at `x`: the end
    /// face's at an end; square to the edge where the band meets an edge,
    /// the spine running along it there; else the plane of the ball's arc,
    /// square to the spine. Each normal points along the edge.
    fn plane_at(&self, x: f64, k: usize, seat: Seat) -> Option<(Point, Vector)> {
        let last = self.pieces.len();
        if k == 0 {
            return Some((self.on_edge(x), -self.end_normals[0]));
        }
        if k == last {
            return Some((self.on_edge(x), self.end_normals[1]));
        }
        let square = self.pieces[k - 1].seat == Seat::Face || self.pieces[k].seat == Seat::Face;
        let (ball, _) = self.ball(x, Some(seat))?;
        if square {
            return Some((ball.centre, self.along));
        }
        Some((ball.centre, self.spine(x, seat)?))
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
        let (at, normal) = self.plane_at(x, k, piece.seat)?;
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

    /// The balls past the end `k` of the edge that still reach inside its
    /// end face, out to where the round stands wholly past that face: the
    /// station from which every ball's arc, and the region between it and
    /// the edge, lie past the end face's plane on to the end of the scan,
    /// and the seat of each station before it that reaches inside. `None`
    /// where a ball there cannot be found, or none stands past.
    fn past_end(&self, k: usize) -> Option<(f64, Vec<Seat>)> {
        let r = self.radius;
        let (from, out) = if k == 0 {
            (self.real.0, -1.0)
        } else {
            (self.real.1, 1.0)
        };
        let end = self.on_edge(from);
        let normal = self.end_normals[k];
        let lean = normal.dot(self.along * out);
        let margin = r * 1e-3;
        // Far enough for a ball's centre, standing off the edge by at most
        // the band's offset and its rise, to stand a radius past the face.
        let offset = (self.band(0.0) - self.origin).magnitude() + r;
        let span = (r + offset) / lean + r;
        let step = r / 64.0;
        let mut seats = Vec::new();
        let mut clear: Option<f64> = None;
        let mut j = 1.0;
        while step * j <= span {
            let t = from + out * step * j;
            let (ball, seat) = self.ball(t, None)?;
            let on_edge = self.on_edge(t);
            let mut points = self.across_b_points(on_edge, ball.on_b);
            points.push(ball.on_a);
            for u in 1..8 {
                points.push(self.point(t, f64::from(u) / 8.0, seat)?);
            }
            if points.iter().all(|p| (*p - end).dot(normal) >= margin) {
                clear.get_or_insert(t);
            } else {
                clear = None;
                seats.push(seat);
            }
            j += 1.0;
        }
        Some((clear?, seats))
    }

    /// The hold with each end that needs it run on past its end face: an
    /// end whose balls past the face that still reach inside it change
    /// seat, or whose piece there rests on one point, or one `forced`. The
    /// round then runs on to where it stands wholly past the end face and
    /// ends there ([`Self::end_past`]), and the cut trims it by the end
    /// face. `None` where the balls
    /// past an end cannot be read, or the round run on does not read into
    /// pieces.
    fn run_on(mut self, forced: [bool; 2]) -> Option<Self> {
        let mut ends: [Option<(f64, Vector)>; 2] = [None; 2];
        let last = self.pieces.len() - 1;
        for (k, index) in [(0, 0), (1, last)] {
            let (at, seats) = self.past_end(k)?;
            let seat = self.pieces[index].seat;
            let simple =
                !forced[k] && self.pin(index).is_none() && seats.iter().all(|s| *s == seat);
            if !simple {
                ends[k] = Some(self.end_past(k, at)?);
            }
        }
        if ends.iter().all(Option::is_none) {
            return Some(self);
        }
        let by = ends.map(|e| e.map_or(0.0, |(by, _)| by));
        for (k, end) in ends.iter().enumerate() {
            if let Some((_, normal)) = end {
                let x = if k == 0 { self.real.0 } else { self.real.1 };
                self.ends_at[k] = Some((self.on_edge(x), self.end_normals[k]));
                self.end_normals[k] = *normal;
            }
        }
        self.origin -= self.along * by[0];
        self.centre -= self.along * by[0];
        self.real = (by[0], by[0] + self.len);
        self.len += by[0] + by[1];
        self.pieces = self.read_pieces()?;
        let last = self.pieces.len() - 1;
        for (k, index) in [(0, 0), (1, last)] {
            if ends[k].is_some() && self.pin(index).is_some() {
                return None;
            }
        }
        let (a, b) = self.real;
        let near = self.len * 1e-9;
        self.at_end = self
            .pieces
            .iter()
            .any(|p| p.seat != Seat::Face && (p.from < a + near || p.to > b - near));
        Some(self)
    }

    /// Where the round run on past end `k` ends, from the station `at`
    /// where it stands wholly past the end face, up to a quarter radius
    /// further as far as balls stand there: how far past the edge's end
    /// the edge's line meets the end's plane, and that plane's outward
    /// normal. On the band the plane stands square to the edge; elsewhere
    /// it is the plane of the last ball's arc, square to the path of the
    /// centres, so the round's end is that arc and the region between it
    /// and the edge one planar face. `None` on a cylinder off the band,
    /// where that plane leans more than 60 degrees from square, where the
    /// seat changes between the ball and where the plane meets the edge,
    /// or where the end so placed is not wholly past the end face.
    fn end_past(&self, k: usize, at: f64) -> Option<(f64, Vector)> {
        let r = self.radius;
        let (from, out) = if k == 0 {
            (self.real.0, -1.0)
        } else {
            (self.real.1, 1.0)
        };
        let step = r / 64.0;
        let mut x = at;
        while (x - at) * out < 0.25 * r && self.ball(x + out * step, None).is_some() {
            x += out * step;
        }
        let (ball, seat) = self.ball(x, None)?;
        let (meet, normal) = if seat == Seat::Face {
            (x, self.along * out)
        } else {
            if matches!(self.front, Front::Cylinder { .. }) {
                return None;
            }
            let spine = self.spine(x, seat)? * out;
            if spine.dot(self.along * out) < 0.5 {
                return None;
            }
            let meet = x + (ball.centre - self.on_edge(x)).dot(spine) / self.along.dot(spine);
            (meet, spine)
        };
        let (lo, hi) = (x.min(meet), x.max(meet));
        let mut t = lo;
        while t < hi {
            if self.ball(t, None)?.1 != seat {
                return None;
            }
            t += step;
        }
        let end = self.on_edge(from);
        let real = self.end_normals[k];
        let margin = r * 1e-3;
        let mut points = self.across_b_points(self.on_edge(meet), ball.on_b);
        points.push(ball.on_a);
        for u in 1..8 {
            points.push(self.point(x, f64::from(u) / 8.0, seat)?);
        }
        if points.iter().any(|p| (*p - end).dot(real) < margin) {
            return None;
        }
        Some(((meet - from) * out, normal))
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
                Cliff::Line { .. } | Cliff::General { .. } => t,
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
/// where it lies in that plane, over the range the edge covers.
fn cliff_of(curved: &Curved, held: Id, plane_normal: Vector, tol: Tolerances) -> Option<Cliff> {
    let eps = tol.angular();
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
        curve => {
            let (lo, hi) = (
                edge.range.0.min(edge.range.1),
                edge.range.0.max(edge.range.1),
            );
            if hi - lo <= 0.0 {
                return None;
            }
            let mut samples = Vec::with_capacity(READ_ALONG + 1);
            for k in 0..=READ_ALONG {
                #[expect(clippy::cast_precision_loss, reason = "a few hundred samples")]
                let t = lo + (hi - lo) * k as f64 / READ_ALONG as f64;
                samples.push((t, curve.point_at(t, tol).ok()?));
            }
            let first = samples[0].1;
            samples
                .iter()
                .all(|(_, p)| (*p - first).dot(plane_normal).abs() <= tol.confusion())
                .then(|| Cliff::General {
                    curve: curve.clone(),
                    range: (lo, hi),
                    samples,
                    closed: edge.ends.0 == edge.ends.1,
                    tol,
                })
        }
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
/// edge of theirs within reach of the band does not lie in that plane. The round's own edges, the seams between those faces and
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
        let Some(cliff) = cliff_of(curved, *e, na, tol) else {
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
            Cliff::General { range, closed, .. } => (cliff, range, closed),
        };
        // The vertex at each end of the range: a curve of the third kind
        // may run its range from the edge's last vertex to its first.
        let mut ends_at = [edge.ends.0, edge.ends.1];
        if matches!(cliff, Cliff::General { .. }) && edge.range.0 > edge.range.1 {
            ends_at.swap(0, 1);
        }
        let mut stops = [None, None];
        for (k, v) in ends_at.into_iter().enumerate() {
            if whole || ends.iter().any(|f| touches(curved, *f, v)) {
                continue;
            }
            let at = corners.iter().position(|c| c.id == v).unwrap_or_else(|| {
                corners.push(Corner {
                    id: v,
                    at: curved.points[&v],
                    smooth: false,
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
    for (k, corner) in corners.iter_mut().enumerate() {
        let tangents: Vec<Option<Vector>> = rims
            .iter()
            .filter(|r| r.stops.contains(&Some(k)))
            .map(|r| r.tangent(corner.at))
            .collect();
        corner.smooth = matches!(
            tangents[..],
            [Some(s), Some(t)] if s.cross(t).magnitude() <= tol.angular()
        );
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

/// The band's centre at `p0` on the edge: the radius inside `a`'s plane
/// (outward normal `na`) and inside `b`, on `b`'s offset; on a cylinder, of
/// the two points where the offsets meet, the one nearer the edge. `None`
/// where they do not meet, or the cylinder's offset closes up.
fn band_centre(
    p0: Point,
    along: Vector,
    na: Vector,
    nb: Vector,
    front: Front,
    r: f64,
    eps: f64,
) -> Option<Point> {
    match front {
        Front::Plane => {
            let k = -r / (1.0 + na.dot(nb));
            Some(p0 + (na + nb) * k)
        }
        Front::Cylinder {
            axis,
            radius,
            sense,
        } => {
            let rho = radius - sense * r;
            if rho <= eps {
                return None;
            }
            let o = p0 - na * r;
            let w = along.cross(na);
            let q = o - axis;
            let d = q - along * q.dot(along);
            let half = d.dot(w);
            let disc = half * half - (d.dot(d) - rho * rho);
            if disc < 0.0 {
                return None;
            }
            [-half - disc.sqrt(), -half + disc.sqrt()]
                .into_iter()
                .map(|s| o + w * s)
                .min_by(|a, b| a.distance(p0).total_cmp(&b.distance(p0)))
        }
    }
}

/// `edges` read as one straight round of `radius` between a plane and a
/// plane or a cylinder along the edge, whose contact on the plane runs
/// past sharp edges of it, held as hold edge holds them; `None` for any
/// other request, or where anything else of the solid stands within its
/// reach. The edges may be a chain along one line, each between faces of
/// the same two surfaces.
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
                .is_some_and(|s| {
                    matches!(s, SurfaceGeometry::Plane(_) | SurfaceGeometry::Cylinder(_))
                })
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
    // Each side's faces: planes, or faces of one cylinder whose axis runs
    // along the edge.
    let front_of = |side: &[usize]| -> Option<Front> {
        let mut found: Option<(Point, f64, f64)> = None;
        let mut plane = false;
        for fi in side {
            match &curved.faces[*fi].surface {
                SurfaceGeometry::Plane(_) => plane = true,
                SurfaceGeometry::Cylinder(c) => {
                    let axis = c.cylinder().axis();
                    if axis.direction.vector().cross(along).magnitude() > ang {
                        return None;
                    }
                    let sense = if curved.faces[*fi].reversed {
                        -1.0
                    } else {
                        1.0
                    };
                    let this = (axis.location, c.cylinder().radius(), sense);
                    if let Some((at, r, s)) = found {
                        let q = this.0 - at;
                        let off = (q - along * q.dot(along)).magnitude();
                        if off > eps || (this.1 - r).abs() > eps || s != this.2 {
                            return None;
                        }
                    }
                    found = Some(this);
                }
                _ => return None,
            }
        }
        match found {
            None => Some(Front::Plane),
            Some(_) if plane => None,
            Some((axis, radius, sense)) => Some(Front::Cylinder {
                axis,
                radius,
                sense,
            }),
        }
    };
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
        let (Some(Front::Plane), Some(front)) = (front_of(&a_faces), front_of(&b_faces)) else {
            continue;
        };
        let Some(centre) = band_centre(p0, along, na, nb, front, radius, eps) else {
            continue;
        };
        // On a cylinder an end face leaning over the band cuts `b` in an
        // ellipse: the round runs on past it, cut by the boolean.
        let leaning = [
            end_normals[0].dot(-along) < 1.0 - ang,
            end_normals[1].dot(along) < 1.0 - ang,
        ];
        let forced = leaning.map(|l| l && matches!(front, Front::Cylinder { .. }));
        let across = {
            let v = along.cross(na);
            let contact = centre + na * radius - p0;
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
            front,
            centre,
            across,
            radius,
            rims,
            corners,
            pieces: Vec::new(),
            at_end: false,
            ends_at: [None; 2],
            real: (0.0, len),
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
        let Some(hold) = hold.run_on(forced) else {
            continue;
        };
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
        points.push(ball.on_a);
        // A cylinder bulges between the edge and the contact on it.
        points.extend(hold.across_b_points(hold.on_edge(x), ball.on_b));
    }
    // An end face leaning over the band carries the round past the edge's
    // end.
    for k in [0, hold.pieces.len()] {
        let Some((pb, pa)) = hold.end_cut(k) else {
            return Ok(false);
        };
        points.extend([pb, pa]);
    }
    // Where the foot jumps, the ball's arc in `a`'s plane bulges toward the
    // corner between the two edges.
    for k in 1..hold.pieces.len() {
        let Some(jump) = hold.jump_into(k) else {
            return Ok(false);
        };
        if let Some(j) = jump {
            let w = j.normal.cross(j.from);
            points.extend((0..=16).map(|s| {
                let t = j.sweep * f64::from(s) / 16.0;
                j.hub + (j.from * t.cos() + w * t.sin()) * j.reach
            }));
        }
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
    // Where the foot jumps, the round passes over the edges between the two
    // it rests on, and the faces past them trim it.
    for (k, pair) in hold.pieces.windows(2).enumerate() {
        let (Seat::Rim(i), Seat::Rim(j)) = (pair[0].seat, pair[1].seat) else {
            continue;
        };
        if !hold.pieces[k + 1].jump {
            continue;
        }
        let Some((rims, corners)) = hold.walk(i, j) else {
            return Ok(false);
        };
        for r in rims {
            held.push(hold.rims[r].id);
            ours.extend(faces_of(c, hold.rims[r].id));
        }
        for k in corners {
            let v = hold.corners[k].id;
            ours.extend((0..c.faces.len()).filter(|fi| touches(c, *fi, v)));
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
/// given) is cut by the planes of [`Pencil`] from the plane at its start to
/// the one at its end, so its end lies in the end face's plane; any other
/// is fitted through its balls' arcs, which at each end are the arcs it
/// shares with its neighbours. `None` where two neighbouring sections
/// cross.
fn fitted(
    hold: &Hold,
    seat: Seat,
    (from, to): (f64, f64),
    planes: Option<((Point, Vector), (Point, Vector))>,
    tol: Tolerances,
) -> OgeomResult<Option<ogeom_geom::BSplineSurface>> {
    let pencil = match planes {
        Some((p0, p1)) => {
            let Some(pencil) = Pencil::new(p0, p1) else {
                return Ok(None);
            };
            Some(pencil)
        }
        None => None,
    };
    let mut rows: Vec<Vec<Point>> = Vec::with_capacity(STATIONS + 1);
    let mut cut_by: Option<(Point, Vector)> = None;
    for i in 0..=STATIONS {
        #[expect(clippy::cast_precision_loss, reason = "a hundred stations")]
        let s = 0.5 * (1.0 - (core::f64::consts::PI * i as f64 / STATIONS as f64).cos());
        let row = if let Some(pencil) = &pencil {
            let (at, normal) = pencil.at(s);
            // The ball whose arc the plane meets lies near where it meets
            // the edge.
            let along = normal.dot(hold.along);
            if along <= 0.0 {
                return Ok(None);
            }
            let x = (at - hold.origin).dot(normal) / along;
            // Each section lies wholly ahead of the one before: the
            // sections do not cross, and the rails run one way.
            let row = hold.section(seat, at, normal, x, ACROSS);
            if let (Some(row), Some((at, normal))) = (&row, cut_by)
                && row.iter().any(|p| (*p - at).dot(normal) <= 0.0)
            {
                return Ok(None);
            }
            cut_by = Some((at, normal));
            row
        } else {
            let x = from + (to - from) * s;
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

/// The planes turning from one plane to another about the line the two
/// share: any two of them meet only on that line, so the sections they cut
/// from a piece the line stays clear of never cross. Each plane is a point
/// on it and its unit normal; between parallel planes, the plane moved
/// from the first to the second.
struct Pencil {
    /// A point on the line the planes share, or on the first plane.
    at: Point,
    /// The first plane's unit normal.
    normal: Vector,
    /// The unit normal a quarter turn on from `normal` about the line.
    turn: Vector,
    /// The angle between the two normals, nought where they are parallel.
    angle: f64,
    /// Between parallel planes, the move from the first to the second.
    shift: Vector,
}

impl Pencil {
    /// The planes from `(p0, n0)` to `(p1, n1)`, the normals pointing the
    /// same way along the piece; `None` where a normal is null or the two
    /// face each other.
    fn new((p0, n0): (Point, Vector), (p1, n1): (Point, Vector)) -> Option<Self> {
        let (m0, m1) = (n0.magnitude(), n1.magnitude());
        if m0 <= 0.0 || m1 <= 0.0 {
            return None;
        }
        let (n0, n1) = (n0 / m0, n1 / m1);
        let cos = n0.dot(n1);
        if cos <= 0.0 {
            return None;
        }
        let axis = n0.cross(n1);
        let sin = axis.magnitude();
        if sin <= 1e-12 {
            return Some(Self {
                at: p0,
                normal: n0,
                turn: Vector::ZERO,
                angle: 0.0,
                shift: n0 * (p1 - p0).dot(n0),
            });
        }
        // The point of the shared line nearest the origin: the mix of `n0`
        // and `n1` lying in both planes.
        let (h0, h1) = (n0.dot(p0 - Point::ORIGIN), n1.dot(p1 - Point::ORIGIN));
        let s2 = sin * sin;
        let at = Point::ORIGIN + n0 * ((h0 - h1 * cos) / s2) + n1 * ((h1 - h0 * cos) / s2);
        Some(Self {
            at,
            normal: n0,
            turn: (axis / sin).cross(n0),
            angle: sin.atan2(cos),
            shift: Vector::ZERO,
        })
    }

    /// The plane a fraction `s` of the way from the first to the second.
    fn at(&self, s: f64) -> (Point, Vector) {
        if self.angle == 0.0 {
            return (self.at + self.shift * s, self.normal);
        }
        let (sin, cos) = (self.angle * s).sin_cos();
        (self.at, self.normal * cos + self.turn * sin)
    }
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

/// The surface of `b` over the tool's reach: a cylinder's height range
/// grown to cover the round run on past its end faces.
fn front_surface(hold: &Hold) -> OgeomResult<SurfaceGeometry> {
    let surface = &hold.curved.faces[hold.b].surface;
    let SurfaceGeometry::Cylinder(c) = surface else {
        return Ok(surface.clone());
    };
    if hold.ends_at.iter().all(Option::is_none) {
        return Ok(surface.clone());
    }
    let cylinder = c.cylinder();
    let axis = cylinder.axis();
    let height = |p: Point| (p - axis.location).dot(axis.direction.vector());
    let (h0, h1) = (height(hold.on_edge(0.0)), height(hold.on_edge(hold.len)));
    let pad = 4.0 * hold.radius;
    Ok(CylinderSurface::new(cylinder, (h0.min(h1) - pad, h0.max(h1) + pad))?.into())
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
    // The pieces either side of each break.
    let sides = |k: usize| (k.saturating_sub(1), k.min(last - 1));
    let mut jumps = Vec::with_capacity(breaks.len());
    for k in 0..breaks.len() {
        let jump = if k == 0 || k == last {
            None
        } else {
            let Some(jump) = hold.jump_into(k) else {
                return Ok(None);
            };
            jump
        };
        jumps.push(jump);
    }
    let pins: Vec<Option<Point>> = (0..last).map(|k| hold.pin(k)).collect();
    // An edge square to the spine holds every foot on one point. The piece
    // fitted through those balls' arcs meets the band and the jump's
    // sphere further from tangent than a fitted piece may. Not built.
    if pieces
        .iter()
        .zip(&pins)
        .any(|(p, pin)| matches!(p.seat, Seat::Rim(_)) && pin.is_some())
    {
        return Ok(None);
    }
    // Each piece off the band, fitted once.
    let mut surfaces = Vec::with_capacity(last);
    for (k, piece) in pieces.iter().enumerate() {
        surfaces.push(if piece.seat == Seat::Face {
            None
        } else {
            // A piece whose balls all pass through one point is fitted
            // through their arcs. At an end it is built only run on past
            // the end face, which the cut trims it by; ended in the end
            // piece's section planes it is not.
            let reaches_end = k == 0 || k + 1 == last;
            let ended = (k == 0 && hold.ends_at[0].is_none())
                || (k + 1 == last && hold.ends_at[1].is_none());
            if pins[k].is_some() && ended {
                return Ok(None);
            }
            let planes = if reaches_end {
                let (Some(p0), Some(p1)) = (
                    hold.plane_at(piece.from, k, piece.seat),
                    hold.plane_at(piece.to, k + 1, piece.seat),
                ) else {
                    return Ok(None);
                };
                Some((p0, p1))
            } else {
                None
            };
            let Some(s) = fitted(hold, piece.seat, (piece.from, piece.to), planes, tol)? else {
                return Ok(None);
            };
            Some(s)
        });
    }
    let mut points = vec![hold.on_edge(0.0), hold.on_edge(hold.len)];
    // A point every ball of a piece sets its foot on is one point for
    // every break that piece reaches.
    let mut pinned: Vec<Option<usize>> = vec![None; last];
    let foot = |j: usize, at: Point, points: &mut Vec<Point>, pinned: &mut Vec<Option<usize>>| {
        let Some(pin) = pins[j] else {
            points.push(at);
            return points.len() - 1;
        };
        *pinned[j].get_or_insert_with(|| {
            points.push(pin);
            points.len() - 1
        })
    };
    // At each break, its contact on `b`, and its feet on `a` or the held
    // edges, the piece before's and the piece after's, two apart where the
    // foot jumps.
    let mut at_b = Vec::with_capacity(breaks.len());
    let mut at_a: Vec<(usize, usize)> = Vec::with_capacity(breaks.len());
    for (k, &x) in breaks.iter().enumerate() {
        let (before, after) = sides(k);
        if let Some(jump) = jumps[k] {
            at_b.push(points.len());
            points.push(jump.on_b);
            let left = foot(before, jump.feet[0], &mut points, &mut pinned);
            let right = foot(after, jump.feet[1], &mut points, &mut pinned);
            at_a.push((left, right));
            continue;
        }
        let (pb, pa) = if k == 0 || k == last {
            let Some(cut) = hold.end_cut(k) else {
                return Ok(None);
            };
            cut
        } else {
            let Some((ball, _)) = hold.ball(x, Some(pieces[after].seat)) else {
                return Ok(None);
            };
            (ball.on_b, ball.on_a)
        };
        at_b.push(points.len());
        points.push(pb);
        let one = match (
            pins[before].is_some() && k > 0,
            pins[after].is_some() && k < last,
        ) {
            (true, true) => {
                let i = foot(before, pa, &mut points, &mut pinned);
                pinned[after] = Some(i);
                i
            }
            (true, false) => foot(before, pa, &mut points, &mut pinned),
            (false, true) => foot(after, pa, &mut points, &mut pinned),
            (false, false) => {
                points.push(pa);
                points.len() - 1
            }
        };
        at_a.push((one, one));
    }
    let (e0, e1) = (0, 1);
    let mut curves: Vec<SheetCurve> = Vec::new();
    let mut push = |ends: (usize, usize), (curve, range): (Curve, (f64, f64))| {
        curves.push(SheetCurve { ends, curve, range });
        curves.len() - 1
    };
    let edge_line = push((e0, e1), segment(points[e0], points[e1], tol)?);
    let b_start = push(
        (e0, at_b[0]),
        hold.across_b(points[e0], points[at_b[0]], tol)?,
    );
    let a_start = push(
        (e0, at_a[0].1),
        segment(points[e0], points[at_a[0].1], tol)?,
    );
    let b_end = push(
        (e1, at_b[last]),
        hold.across_b(points[e1], points[at_b[last]], tol)?,
    );
    let a_end = push(
        (e1, at_a[last].0),
        segment(points[e1], points[at_a[last].0], tol)?,
    );
    // The section at each break, the piece before's and the piece after's:
    // the band's arc, the arc of the ball two pieces share, or a fitted
    // piece's end in an end face's plane. A piece whose feet stay on one
    // point ends in its own boundary, which its surface holds exactly up
    // to the point it closes to. Where the foot jumps, the two are the
    // ball's arcs to its two feet.
    let mut sections: Vec<(usize, usize)> = Vec::with_capacity(breaks.len());
    for (k, &x) in breaks.iter().enumerate() {
        let (before, after) = sides(k);
        let end = k == 0 || k == last;
        let mut side = |j: usize, a: usize| -> OgeomResult<Option<usize>> {
            let ends = (at_b[k], a);
            let own = end || pins[j].is_some() || jumps[k].is_some();
            if own && let Some(surface) = &surfaces[j] {
                let Some(side) = boundary_from(surface, &points, ends, tol)? else {
                    return Ok(None);
                };
                return Ok(Some(push(side.ends, (side.curve, side.range))));
            }
            let (at, normal) = if k == 0 {
                (hold.on_edge(x), -hold.end_normals[0])
            } else if k == last {
                (hold.on_edge(x), hold.end_normals[1])
            } else {
                (points[at_b[k]], hold.along)
            };
            if end && normal.cross(hold.along).magnitude() > tol.angular() {
                // The band cut by an end face leaning over it: an arc of the
                // ellipse its plane cuts from the band's cylinder.
                let Some((flip, cut)) =
                    band_end(hold, at, normal, (points[ends.0], points[ends.1]), tol)?
                else {
                    return Ok(None);
                };
                return Ok(Some(push(if flip { (ends.1, ends.0) } else { ends }, cut)));
            }
            let Some((ball, _)) = hold.ball(x, Some(pieces[j].seat)) else {
                return Ok(None);
            };
            let (u, w) = (ball.on_b - ball.centre, points[a] - ball.centre);
            let z = u.cross(w);
            let sweep = (u.dot(w) / (r * r)).clamp(-1.0, 1.0).acos();
            Ok(Some(push(ends, arc(ball.centre, z, u, r, sweep, tol)?)))
        };
        let pair = if jumps[k].is_some() {
            let (Some(left), Some(right)) = (side(before, at_a[k].0)?, side(after, at_a[k].1)?)
            else {
                return Ok(None);
            };
            (left, right)
        } else {
            // Away from the ends, a piece on one point holds the section.
            let j = if end {
                if k == 0 { after } else { before }
            } else if pins[before].is_some() {
                before
            } else {
                after
            };
            let Some(one) = side(j, at_a[k].0)? else {
                return Ok(None);
            };
            (one, one)
        };
        sections.push(pair);
    }
    // Along each piece, its contact on `b` and on `a` or the held edge; a
    // piece on one point has its one foot there.
    let mut along_b = Vec::with_capacity(last);
    let mut along_a = Vec::with_capacity(last);
    for (k, surface) in surfaces.iter().enumerate() {
        let (b0, b1, a0, a1) = (at_b[k], at_b[k + 1], at_a[k].1, at_a[k + 1].0);
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
    // Where the foot jumps, the arc in `a`'s plane between the two feet.
    let mut tops = Vec::new();
    for (k, jump) in jumps.iter().enumerate() {
        if let Some(j) = jump {
            let curve = arc(j.hub, j.normal, j.from, j.reach, j.sweep, tol)?;
            tops.push((k, push(at_a[k], curve)));
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
            surface: front_surface(hold)?,
            curves: [vec![edge_line, b_end], along_b.clone(), vec![b_start]].concat(),
            at: on_mid + (mb - on_mid) * 0.5,
            out: hold.nb_at(on_mid + (mb - on_mid) * 0.5),
        },
        ToolFace {
            surface: curved.faces[hold.a].surface.clone(),
            curves: [
                vec![edge_line, a_end],
                along_a.iter().flatten().copied().collect(),
                tops.iter().map(|(_, c)| *c).collect(),
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
        let (pb, pa) = (points[at_b[k]], points[at_a[k].0]);
        let side = usize::from(k != 0);
        let surface = if hold.ends_at[side].is_some() {
            // Run on past the end face: the plane the round ends in.
            let x = hold.na - out * hold.na.dot(out);
            let frame = Frame::new(pe, Direction::new(out, tol)?, Direction::new(x, tol)?, tol)?;
            SurfaceGeometry::Plane(ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(frame)))
        } else {
            curved.faces[end].surface.clone()
        };
        faces.push(ToolFace {
            surface,
            curves: if k == 0 {
                vec![b_start, sections[0].1, a_start]
            } else {
                vec![b_end, sections[last].0, a_end]
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
        let mut loop_ = vec![along_b[k], sections[k + 1].0];
        loop_.extend(along_a[k]);
        loop_.push(sections[k].1);
        faces.push(ToolFace {
            surface,
            curves: loop_,
            at: p,
            out: (ball.centre - p) / r,
        });
    }
    // Where the foot jumps, the ball's own sphere between its two arcs, its
    // poles along the spine's turn and its seam behind the lune.
    for (k, top) in tops {
        let Some(j) = jumps[k] else {
            return Ok(None);
        };
        let w = j.normal.cross(j.from);
        let half = 0.5 * j.sweep;
        let on_top = j.hub + (j.from * half.cos() + w * half.sin()) * j.reach;
        let lune = (on_top - j.centre) + (j.on_b - j.centre);
        let lune = lune / lune.magnitude();
        let pole = j.spines[0] + j.spines[1];
        let pole = pole / pole.magnitude();
        let seam = -(lune - pole * lune.dot(pole));
        let frame = Frame::new(
            j.centre,
            Direction::new(pole, tol)?,
            Direction::new(seam, tol)?,
            tol,
        )?;
        let sphere = ogeom_math::Sphere::new(frame, r, tol)?;
        faces.push(ToolFace {
            surface: SurfaceGeometry::Sphere(ogeom_geom::SphereSurface::new(sphere)),
            curves: vec![sections[k].0, top, sections[k].1],
            at: j.centre + lune * r,
            out: -lune,
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
    // Where the foot jumps the tool passes over the faces past the edges,
    // and where it runs on past an end face over the face; the cut takes
    // its volume within the solid.
    let inside = if hold.pieces.iter().any(|p| p.jump) || hold.ends_at.iter().any(Option::is_some) {
        let Ok(within) = ogeom_bool::common(model, &cut, solid, tol) else {
            return Ok(None);
        };
        within.shape
    } else {
        cut.clone()
    };
    let taken = ogeom_algo::volume_properties(model, &inside, deflection, tol)?.mass;
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
