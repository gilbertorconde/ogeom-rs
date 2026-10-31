//! Face blends: `blend_faces` and `fillet_faces` on generated solids,
//! judged against the rolling ball worked out from the inputs alone.
//!
//! Every case knows its solid analytically: the two faces' supports (planes,
//! cylinders, cones), every other face, and the strip a chamfer, an old
//! round or a wedge left between the two. From that alone it finds the
//! balls of the blend's radius touching both supports on the side the
//! operation rolls on, and asks three things of them.
//!
//! - Does a ball roll there at all? Its contacts must lie on the two faces,
//!   read to the chord the kernel's blends read faces at ([`reading`]); and
//!   no other face of the input may come inside a ball's lune, the part of
//!   it between its centre and its arc, where the round's material (or its
//!   void) must be. A face poking into the back of a ball, away from the
//!   corner, changes nothing the round touches, as an edge round's ball may
//!   poke through a thin plate. The lunes are searched by branch and bound
//!   over the face's parameters and the run's, each box's gap bounded below
//!   by Lipschitz constants of the face, of the centre's path and of the
//!   contacts' turn; a face found deeper than [`TOUCH`] inside a lune
//!   blocks the round, and one held clear of every lune by more than
//!   [`TOUCH`] lets it be. Faces the run ends on (a prism's caps) are not
//!   asked.
//! - Where the kernel built a round, does every new face lie on the balls'
//!   envelope? Each point sampled on a face that lies on no input surface
//!   (the mesh's vertices, and the surface at its triangles' centres) must
//!   stand the radius from the nearest centre within the face's stated
//!   tolerance, and inside the arc between that ball's contacts.
//! - Is the volume right? A round takes off or adds the region between the
//!   ball's arc and the corner of the two supports, the strip with it. For
//!   a section swept straight (planes, and cylinders sharing a direction)
//!   that is the section's area times the run, in closed form for two
//!   planes and integrated on the section's lines and arcs otherwise;
//!   turned about an axis it is Pappus' ring. A corner with no shared
//!   direction or axis is integrated over the balls' normal sections along
//!   the centre's closed path, each a fan of rays from the centre between
//!   the contacts, out to the first support.
//!
//! A case is `ok` where the kernel refused a corner no ball rolls in, or
//! built the right shape where one does. A shape built where no ball rolls
//! is `wrong`, as is one off its envelope or its volume. Refusing a corner a
//! ball does roll in is an honest `refused`, and its note says so. A
//! marched run whose contact leaves a face over part of it may end in caps,
//! which nothing here describes: there only `check` judges the result.
//! `fillet_faces` returns no solid, so its rounds are judged on the
//! envelope alone.
//!
//! Strips built by the kernel on a marched corner (an old round from
//! `fillet_edges`) have no closed form here; a smaller ball's round lies
//! wholly inside a larger one's corner (the morphological closing and
//! opening are monotone in the ball), so such a strip blocks a round of
//! radius `r` exactly where its own radius exceeds `r`.

use core::f64::consts::{PI, TAU};
use std::sync::OnceLock;

use ogeom::math::{Direction, Frame, Point, Point2, Vector, Vector2};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

use super::{Outcome, Rng, T, Verdict, across, allowance, judged, measure};

type P2 = Point2;
type V2 = Vector2;

/// How deep a point of another face may stand inside a ball and still count
/// as touching it: ten times the confusion distance the kernel refuses at.
/// The lens a ball cuts that deep holds under 1e-9 of a millimetre-scale
/// run, below every volume slack here.
const TOUCH: f64 = 1e-6;

/// How far past its face's edge a contact may stand and still be read as
/// on it, either way: the chord the kernel's blends read faces at for
/// whether a ball touches them, a thousandth of the radius.
fn reading(r: f64) -> f64 {
    (r * 1e-3).max(TOUCH)
}

/// The families, each its own generator; a case's key names one.
pub(super) const FAMILIES: [&str; 13] = [
    "planes", "notch", "dshaft", "rod", "drums", "rim", "bore", "boss", "frustum", "conefoot",
    "leaning", "tee", "slanted",
];

// --- numerics --------------------------------------------------------------

/// Gauss-Legendre nodes and weights on `[-1, 1]`, 32 of them: exact for
/// polynomials to degree 63, and to rounding on the smooth integrands here.
fn gauss() -> &'static [(f64, f64)] {
    static NODES: OnceLock<Vec<(f64, f64)>> = OnceLock::new();
    NODES.get_or_init(|| {
        let n = 32;
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let mut x = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
            let mut slope = 1.0;
            for _ in 0..100 {
                let (mut p0, mut p1) = (1.0, x);
                for k in 2..=n {
                    let k = k as f64;
                    let p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k;
                    p0 = p1;
                    p1 = p2;
                }
                slope = n as f64 * (x * p1 - p0) / (x * x - 1.0);
                let step = p1 / slope;
                x -= step;
                if step.abs() < 1e-16 {
                    break;
                }
            }
            out.push((x, 2.0 / ((1.0 - x * x) * slope * slope)));
        }
        out
    })
}

fn integrate(a: f64, b: f64, mut f: impl FnMut(f64) -> f64) -> f64 {
    let (mid, half) = (f64::midpoint(a, b), (b - a) / 2.0);
    gauss()
        .iter()
        .map(|(x, w)| w * f(half.mul_add(*x, mid)))
        .sum::<f64>()
        * half
}

/// The angle in `(-pi, pi]` equal to `a` modulo a turn.
fn wrap(a: f64) -> f64 {
    let mut a = a % TAU;
    if a > PI {
        a -= TAU;
    } else if a <= -PI {
        a += TAU;
    }
    a
}

fn angle_between(a: Vector, b: Vector) -> f64 {
    a.cross(b).magnitude().atan2(a.dot(b))
}

/// Whether a ball rolls in the corner, as the inputs decide it.
#[derive(Clone, Debug)]
enum Rolls {
    Yes,
    No(String),
    /// Not decided within the bounds' budget, or at a margin inside
    /// [`TOUCH`]: either answer is taken.
    Unsure(String),
    /// A contact leaves a face over part of a marched run: the round may
    /// stop there in a cap, which no closed form here describes, so only
    /// what `check` says of it is judged.
    Capped(String),
}

impl Rolls {
    /// The more decisive of two findings: any block blocks.
    fn and(self, other: Self) -> Self {
        match (self, other) {
            (Self::No(a), _) | (_, Self::No(a)) => Self::No(a),
            (Self::Capped(a), _) | (_, Self::Capped(a)) => Self::Capped(a),
            (Self::Unsure(a), _) | (_, Self::Unsure(a)) => Self::Unsure(a),
            _ => Self::Yes,
        }
    }

    /// A margin read as a finding: room past `band`, a block short of
    /// `-band`, and either answer taken between.
    fn of_margin(margin: f64, band: f64, why: impl FnOnce() -> String) -> Self {
        if margin > band {
            Self::Yes
        } else if margin >= -band {
            Self::Unsure(why())
        } else {
            Self::No(why())
        }
    }
}

// --- the section: lines and circles in a plane -----------------------------

/// A support's trace in the section, signed distance positive away from the
/// material.
#[derive(Clone, Copy, Debug)]
enum Curve2 {
    Line {
        at: P2,
        out: V2,
    },
    Circle {
        centre: P2,
        radius: f64,
        inside: bool,
    },
}

impl Curve2 {
    /// The line through `p` and `q`, the material on the side of `m`.
    fn line(p: P2, q: P2, m: P2) -> Self {
        let d = q - p;
        let d = d / d.magnitude();
        let mut out = V2::new(-d.y, d.x);
        if (m - p).dot(out) > 0.0 {
            out = -out;
        }
        Self::Line { at: p, out }
    }

    fn signed(self, q: P2) -> f64 {
        match self {
            Self::Line { at, out } => out.dot(q - at),
            Self::Circle {
                centre,
                radius,
                inside,
            } => {
                let d = q.distance(centre);
                if inside { d - radius } else { radius - d }
            }
        }
    }

    fn foot(self, q: P2) -> P2 {
        match self {
            Self::Line { out, .. } => q - out * self.signed(q),
            Self::Circle { centre, radius, .. } => {
                let d = q - centre;
                centre + d * (radius / d.magnitude())
            }
        }
    }

    /// The locus where the signed distance is `d`.
    fn offset(self, d: f64) -> Self {
        match self {
            Self::Line { at, out } => Self::Line {
                at: at + out * d,
                out,
            },
            Self::Circle {
                centre,
                radius,
                inside,
            } => Self::Circle {
                centre,
                radius: if inside { radius + d } else { radius - d },
                inside,
            },
        }
    }

    /// The piece of this trace from `p` to `q`, the short way.
    fn piece(self, p: P2, q: P2) -> Piece {
        match self {
            Self::Line { .. } => Piece::Seg(p, q),
            Self::Circle { centre, radius, .. } => Piece::arc(centre, radius, p, q),
        }
    }
}

/// The meeting of two traces nearest `near`.
fn meet(a: Curve2, b: Curve2, near: P2) -> Option<P2> {
    let line_circle = |at: P2, out: V2, centre: P2, radius: f64| -> Vec<P2> {
        let along = V2::new(-out.y, out.x);
        let w = at - centre;
        let (b, c) = (along.dot(w), w.dot(w) - radius * radius);
        let disc = b * b - c;
        if disc < 0.0 || radius <= 0.0 {
            return Vec::new();
        }
        let s = disc.sqrt();
        vec![at + along * (-b - s), at + along * (-b + s)]
    };
    let found: Vec<P2> = match (a, b) {
        (Curve2::Line { at: pa, out: na }, Curve2::Line { at: pb, out: nb }) => {
            let det = na.cross(nb);
            if det.abs() < 1e-12 {
                return None;
            }
            let (ca, cb) = (na.dot(pa.to_vector()), nb.dot(pb.to_vector()));
            vec![P2::new(
                (ca * nb.y - cb * na.y) / det,
                (na.x * cb - nb.x * ca) / det,
            )]
        }
        (Curve2::Line { at, out }, Curve2::Circle { centre, radius, .. })
        | (Curve2::Circle { centre, radius, .. }, Curve2::Line { at, out }) => {
            line_circle(at, out, centre, radius)
        }
        (
            Curve2::Circle {
                centre: c0,
                radius: r0,
                ..
            },
            Curve2::Circle {
                centre: c1,
                radius: r1,
                ..
            },
        ) => {
            let d = c0.distance(c1);
            if r0 <= 0.0 || r1 <= 0.0 || d > r0 + r1 || d < (r0 - r1).abs() || d < 1e-12 {
                return None;
            }
            let along = (c1 - c0) / d;
            let x = (d * d + r0 * r0 - r1 * r1) / (2.0 * d);
            let y = (r0 * r0 - x * x).max(0.0).sqrt();
            let side = V2::new(-along.y, along.x);
            vec![c0 + along * x + side * y, c0 + along * x - side * y]
        }
    };
    found
        .into_iter()
        .min_by(|p, q| p.distance(near).total_cmp(&q.distance(near)))
}

/// A piece of a section boundary.
#[derive(Clone, Copy, Debug)]
enum Piece {
    Seg(P2, P2),
    Arc {
        centre: P2,
        radius: f64,
        from: f64,
        to: f64,
    },
}

impl Piece {
    /// The arc about `centre` from `p` to `q`, the short way round.
    fn arc(centre: P2, radius: f64, p: P2, q: P2) -> Self {
        let from = (p.y - centre.y).atan2(p.x - centre.x);
        let to = from + wrap((q.y - centre.y).atan2(q.x - centre.x) - from);
        Self::Arc {
            centre,
            radius,
            from,
            to,
        }
    }

    fn at(self, t: f64) -> P2 {
        match self {
            Self::Seg(a, b) => a + (b - a) * t,
            Self::Arc {
                centre,
                radius,
                from,
                to,
            } => {
                let a = (to - from).mul_add(t, from);
                centre + V2::new(a.cos(), a.sin()) * radius
            }
        }
    }

    fn velocity(self, t: f64) -> V2 {
        match self {
            Self::Seg(a, b) => b - a,
            Self::Arc {
                radius, from, to, ..
            } => {
                let a = (to - from).mul_add(t, from);
                V2::new(-a.sin(), a.cos()) * (radius * (to - from))
            }
        }
    }

    fn length(self) -> f64 {
        match self {
            Self::Seg(a, b) => a.distance(b),
            Self::Arc {
                radius, from, to, ..
            } => radius * (to - from).abs(),
        }
    }

    /// The distance from `q` to this piece's line or circle.
    fn off_carrier(self, q: P2) -> f64 {
        match self {
            Self::Seg(a, b) => {
                let d = b - a;
                (d.cross(q - a) / d.magnitude()).abs()
            }
            Self::Arc { centre, radius, .. } => (q.distance(centre) - radius).abs(),
        }
    }
}

/// The signed area of a closed boundary and its first moment about the
/// section's `y` axis (`x` times the area), both by Green's theorem.
fn region(pieces: &[Piece]) -> (f64, f64) {
    let (mut area, mut moment) = (0.0, 0.0);
    for p in pieces {
        area += integrate(0.0, 1.0, |t| {
            let (q, v) = (p.at(t), p.velocity(t));
            0.5 * (q.x * v.y - q.y * v.x)
        });
        moment += integrate(0.0, 1.0, |t| {
            let (q, v) = (p.at(t), p.velocity(t));
            0.5 * q.x * q.x * v.y
        });
    }
    (area, moment)
}

/// The planes through a ball's centre and each contact, as the unit
/// normals pointing toward the other contact: the ball's lune between them
/// is the part of it between the centre and its arc.
fn lune_normals<V>(ua: V, ub: V) -> [V; 2]
where
    V: Copy
        + core::ops::Sub<Output = V>
        + core::ops::Mul<f64, Output = V>
        + core::ops::Div<f64, Output = V>,
    V: Dot,
{
    let k = ua.dot(ub);
    let m0 = ub - ua * k;
    let m1 = ua - ub * k;
    [m0 / m0.dot(m0).sqrt(), m1 / m1.dot(m1).sqrt()]
}

trait Dot {
    fn dot(self, other: Self) -> f64;
}

impl Dot for V2 {
    fn dot(self, other: Self) -> f64 {
        Vector2::dot(self, other)
    }
}

impl Dot for Vector {
    fn dot(self, other: Self) -> f64 {
        Vector::dot(self, other)
    }
}

/// How far `w` (a point less the ball's centre) stands outside the ball's
/// lune: the largest of its distance past the sphere and past either plane,
/// the latter held to at most `r`. One-Lipschitz in the point, negative
/// strictly inside the lune, and below `-d` only where the point is more
/// than `d` inside the ball and both planes.
fn lune_gap<V: Dot + Copy>(w: V, r: f64, m: [V; 2]) -> f64 {
    let ball = w.dot(w).sqrt() - r;
    ball.max((-w.dot(m[0])).min(r)).max((-w.dot(m[1])).min(r))
}

/// Whether a point of `piece` lies inside the lune of the disc of radius
/// `r` about `centre` between its contacts, by halving the piece under its
/// length as Lipschitz constant.
fn piece_inside(piece: Piece, centre: P2, contacts: [P2; 2], r: f64) -> Rolls {
    let unit = |p: P2| {
        let d = p - centre;
        d / d.magnitude()
    };
    let m = lune_normals(unit(contacts[0]), unit(contacts[1]));
    let length = piece.length();
    let mut stack: Vec<(f64, f64)> = (0..16)
        .map(|i| (f64::from(i) / 16.0, f64::from(i + 1) / 16.0))
        .collect();
    let mut budget = 200_000;
    while let Some((a, b)) = stack.pop() {
        budget -= 1;
        if budget == 0 {
            return Rolls::Unsure("a face's distance to the ball was not bounded".into());
        }
        let mid = f64::midpoint(a, b);
        let q = piece.at(mid);
        let gap = lune_gap(q - centre, r, m);
        if gap < -TOUCH {
            return Rolls::No(format!(
                "another face stands {:.3e} inside the ball at ({:.6}, {:.6}) of the section",
                -gap, q.x, q.y
            ));
        }
        if gap - length * (b - a) / 2.0 >= -TOUCH {
            continue;
        }
        stack.push((a, mid));
        stack.push((mid, b));
    }
    Rolls::Yes
}

// --- the section swept -----------------------------------------------------

/// How a section makes the solid: swept straight along a frame's `z`, or
/// turned about it (the section's `x` the radius, `y` the height).
#[derive(Clone, Copy)]
enum Sweep {
    Prism { frame: Frame, length: f64 },
    Revolve { frame: Frame },
}

impl Sweep {
    fn frame(self) -> Frame {
        match self {
            Self::Prism { frame, .. } | Self::Revolve { frame } => frame,
        }
    }

    fn section(self, p: Point) -> P2 {
        let f = self.frame();
        let d = p - f.origin();
        match self {
            Self::Prism { .. } => P2::new(d.dot(f.x().vector()), d.dot(f.y().vector())),
            Self::Revolve { .. } => {
                let h = d.dot(f.z().vector());
                P2::new((d - f.z().vector() * h).magnitude(), h)
            }
        }
    }

    /// The point of the section at `q`, `w` along the prism or `w` radians
    /// round the axis from the frame's `x`.
    fn place(self, q: P2, w: f64) -> Point {
        let f = self.frame();
        let (x, y, z) = (f.x().vector(), f.y().vector(), f.z().vector());
        match self {
            Self::Prism { .. } => f.origin() + x * q.x + y * q.y + z * w,
            Self::Revolve { .. } => f.origin() + z * q.y + (x * w.cos() + y * w.sin()) * q.x,
        }
    }
}

/// A face's reach along its support: how far a contact `q` stands inside
/// the face (negative outside it).
type Reach2 = Box<dyn Fn(P2) -> f64>;

/// A line support's reach: the contact's distance from `v` along `dir`
/// between `lo` and `hi`.
fn along(v: P2, dir: V2, lo: f64, hi: f64) -> Reach2 {
    Box::new(move |q| {
        let s = (q - v).dot(dir);
        (s - lo).min(hi - s)
    })
}

/// A circle support's reach: the arc length from `v` round `centre`, the
/// way `sense` gives, between `lo` and `hi`.
fn round_from(centre: P2, radius: f64, v: P2, sense: f64, lo: f64, hi: f64) -> Reach2 {
    let start = (v.y - centre.y).atan2(v.x - centre.x);
    Box::new(move |q| {
        let s = sense * wrap((q.y - centre.y).atan2(q.x - centre.x) - start) * radius;
        (s - lo).min(hi - s)
    })
}

/// The corner as a section: two supports meeting at `crease`, the faces'
/// reach on each, the section's other boundary pieces, and what the solid
/// is made of.
struct Section {
    sweep: Sweep,
    supports: [Curve2; 2],
    crease: P2,
    reach: [Reach2; 2],
    others: Vec<Piece>,
    /// The second crease the two faces share, as the reflection taking the
    /// first to it.
    mirror: Option<fn(P2) -> P2>,
    /// Section carriers of every input face, for telling new faces from
    /// old: the supports, the other pieces' lines and circles.
    carriers: Vec<Piece>,
    /// Planes of input faces outside the section (a prism's caps, a
    /// block's sides about a revolved boss).
    planes: Vec<(Point, Vector)>,
    /// Where the corner's two supports are planes, the angle between the
    /// faces across the ball's side: its round has a closed form.
    angle: Option<f64>,
    /// The ball's distance from the axis must stay under this, a revolved
    /// part's non-revolved faces beyond it.
    within: Option<f64>,
}

/// A ball of the section: its centre and contacts.
struct Ball2 {
    centre: P2,
    contacts: [P2; 2],
}

impl Section {
    fn ball(&self, sides: [f64; 2], r: f64) -> Option<Ball2> {
        let centre = meet(
            self.supports[0].offset(sides[0] * r),
            self.supports[1].offset(sides[1] * r),
            self.crease,
        )?;
        let contacts = [self.supports[0].foot(centre), self.supports[1].foot(centre)];
        // The other root of a line and a circle can lie far off; the ball in
        // this corner touches both supports near the crease.
        let scale = 1.0 + self.crease.to_vector().magnitude();
        for (i, s) in self.supports.iter().enumerate() {
            if (s.signed(centre) - sides[i] * r).abs() > 1e-9 * scale {
                return None;
            }
        }
        Some(Ball2 { centre, contacts })
    }

    fn rolls(&self, sides: [f64; 2], r: f64, solid: bool) -> Rolls {
        let Some(ball) = self.ball(sides, r) else {
            return Rolls::No("no ball touches both supports beside the crease".into());
        };
        let mut found = Rolls::Yes;
        for i in 0..2 {
            let q = ball.contacts[i];
            let margin = (self.reach[i])(q);
            found = found.and(Rolls::of_margin(margin, reading(r), || {
                format!(
                    "the ball's contact ({:.6}, {:.6}) is {:.3e} off face {}",
                    q.x,
                    q.y,
                    -margin,
                    "ab".as_bytes()[i] as char
                )
            }));
        }
        if solid {
            for p in &self.others {
                found = found.and(piece_inside(*p, ball.centre, ball.contacts, r));
            }
        }
        found
    }

    /// The section of the corner between the ball and the crease.
    fn corner(&self, ball: &Ball2) -> [Piece; 3] {
        let [a, b] = self.supports;
        [
            a.piece(self.crease, ball.contacts[0]),
            Piece::arc(
                ball.centre,
                (ball.centre - ball.contacts[0]).magnitude(),
                ball.contacts[0],
                ball.contacts[1],
            ),
            b.piece(ball.contacts[1], self.crease),
        ]
    }

    /// The volume between the ball's arc and the crease over the whole run,
    /// every crease the faces share.
    fn corner_volume(&self, ball: &Ball2, r: f64) -> f64 {
        let creases = if self.mirror.is_some() { 2.0 } else { 1.0 };
        let (area, moment) = match self.angle {
            Some(angle) => (
                r * r * (1.0 / (angle / 2.0).tan() - (PI - angle) / 2.0),
                0.0,
            ),
            None => region(&self.corner(ball)),
        };
        creases
            * match self.sweep {
                Sweep::Prism { length, .. } => area.abs() * length,
                Sweep::Revolve { .. } => TAU * moment.abs(),
            }
    }

    /// How far `p` stands off the envelope of the balls and outside the arc
    /// between the nearest ball's contacts.
    fn off_round(&self, ball: &Ball2, p: Point) -> (f64, f64) {
        let q = self.sweep.section(p);
        let mut balls = vec![(ball.centre, ball.contacts)];
        if let Some(m) = self.mirror {
            balls.push((m(ball.centre), ball.contacts.map(m)));
        }
        let (centre, contacts) = balls
            .into_iter()
            .min_by(|a, b| a.0.distance(q).total_cmp(&b.0.distance(q)))
            .unwrap();
        let d = q.distance(centre);
        let r = centre.distance(contacts[0]);
        let to = |p: P2| (p - centre).to_3d();
        let u = to(q);
        let outside = angle_between(u, to(contacts[0])) + angle_between(u, to(contacts[1]))
            - angle_between(to(contacts[0]), to(contacts[1]));
        (d - r, outside)
    }

    fn on_input(&self, p: Point) -> f64 {
        let q = self.sweep.section(p);
        let mut near = f64::INFINITY;
        for c in &self.carriers {
            near = near.min(c.off_carrier(q));
        }
        for (at, n) in &self.planes {
            near = near.min((p - *at).dot(*n).abs());
        }
        near
    }
}

// --- a corner with no shared direction or axis ----------------------------

/// A support in space, signed distance positive away from the material.
#[derive(Clone, Copy)]
enum Surf {
    Plane {
        at: Point,
        out: Vector,
    },
    Cylinder {
        at: Point,
        axis: Vector,
        radius: f64,
        inside: bool,
    },
}

impl Surf {
    fn radial(at: Point, axis: Vector, p: Point) -> Vector {
        let w = p - at;
        w - axis * w.dot(axis)
    }

    fn signed(self, p: Point) -> f64 {
        match self {
            Self::Plane { at, out } => out.dot(p - at),
            Self::Cylinder {
                at,
                axis,
                radius,
                inside,
            } => {
                let d = Self::radial(at, axis, p).magnitude();
                if inside { d - radius } else { radius - d }
            }
        }
    }

    fn gradient(self, p: Point) -> Vector {
        match self {
            Self::Plane { out, .. } => out,
            Self::Cylinder {
                at, axis, inside, ..
            } => {
                let w = Self::radial(at, axis, p);
                let w = w / w.magnitude();
                if inside { w } else { -w }
            }
        }
    }

    fn foot(self, p: Point) -> Point {
        p - self.gradient(p) * self.signed(p)
    }

    /// The first distance past `from` along the ray from `c` in direction
    /// `u` where the ray meets this surface.
    fn hit(self, c: Point, u: Vector, from: f64) -> Option<f64> {
        let from = from * (1.0 - 1e-9);
        match self {
            Self::Plane { out, .. } => {
                let d = out.dot(u);
                if d.abs() < 1e-15 {
                    return None;
                }
                let s = -self.signed(c) / d;
                (s >= from).then_some(s)
            }
            Self::Cylinder {
                at, axis, radius, ..
            } => {
                let w = Self::radial(at, axis, c);
                let up = u - axis * u.dot(axis);
                let a = up.dot(up);
                if a < 1e-18 {
                    return None;
                }
                let b = w.dot(up);
                let disc = b.mul_add(b, -a * (w.dot(w) - radius * radius));
                let disc = if disc < 0.0 && disc > -1e-12 * b * b {
                    0.0
                } else {
                    disc
                };
                if disc < 0.0 {
                    return None;
                }
                let s = disc.sqrt();
                [(-b - s) / a, (-b + s) / a]
                    .into_iter()
                    .filter(|t| *t >= from)
                    .reduce(f64::min)
            }
        }
    }
}

/// A patch of another face, over `[0, 1]` squared.
#[derive(Clone, Copy)]
enum Patch {
    Quad {
        at: Point,
        u: Vector,
        v: Vector,
    },
    Tri {
        a: Point,
        b: Point,
        c: Point,
    },
    Disk {
        centre: Point,
        e1: Vector,
        e2: Vector,
        radius: f64,
    },
}

impl Patch {
    fn at(self, u: f64, v: f64) -> Point {
        match self {
            Self::Quad { at, u: du, v: dv } => at + du * u + dv * v,
            Self::Tri { a, b, c } => a + (b - a) * u + (c - b) * (u * v),
            Self::Disk {
                centre,
                e1,
                e2,
                radius,
            } => {
                let t = TAU * v;
                centre + (e1 * t.cos() + e2 * t.sin()) * (radius * u)
            }
        }
    }

    /// Bounds on the patch's speed in `u` and in `v`.
    fn lipschitz(self) -> (f64, f64) {
        match self {
            Self::Quad { u, v, .. } => (u.magnitude(), v.magnitude()),
            Self::Tri { a, b, c } => (
                (b - a).magnitude() + (c - b).magnitude(),
                (c - b).magnitude(),
            ),
            Self::Disk { radius, .. } => (radius, TAU * radius),
        }
    }

    fn plane(self) -> (Point, Vector) {
        let (o, n) = match self {
            Self::Quad { at, u, v } => (at, u.cross(v)),
            Self::Tri { a, b, c } => (a, (b - a).cross(c - a)),
            Self::Disk { centre, e1, e2, .. } => (centre, e1.cross(e2)),
        };
        (o, n / n.magnitude())
    }
}

/// The centres of a closed run, `t` once round `[0, 2 pi)`.
struct Spine {
    at: Box<dyn Fn(f64) -> Point>,
    /// A bound on the centre's speed in `t`.
    speed: f64,
}

/// How the lunes of a run turn: the directions to the contacts at `t`,
/// and a bound on the turn of the lune's planes per unit `t`.
struct Lunes<'a> {
    spine: &'a Spine,
    supports: [Surf; 2],
    sides: [f64; 2],
    turn: f64,
}

impl<'a> Lunes<'a> {
    fn new(spine: &'a Spine, supports: [Surf; 2], sides: [f64; 2], r: f64) -> Self {
        // A contact direction is a support's unit normal at the centre: a
        // plane's does not turn, and a cylinder's turns at the centre's
        // speed over its distance from the axis, which on the offset is the
        // cylinder's radius moved by the radius.
        let mut omega = 0.0;
        for (s, side) in supports.iter().zip(sides) {
            if let Surf::Cylinder { radius, inside, .. } = s {
                let reach = if *inside {
                    radius + side * r
                } else {
                    radius - side * r
                };
                omega += spine.speed / reach;
            }
        }
        // The angle between the contacts turns no faster than both together;
        // its least over the run is the least at the stations less that turn
        // over half a station.
        let mut least = f64::INFINITY;
        for k in 0..STATIONS {
            let [a, b] = Self::contacts(spine, supports, sides, Marched::station(k));
            least = least.min(angle_between(a, b));
        }
        let least = least - omega * TAU / STATIONS as f64 / 2.0;
        // m = w / |w| with w = b - a (a.b): |w'| <= 2 (|a'| + |b'|), |w| the
        // sine of the angle, so |m'| <= 4 omega / sin.
        let turn = if least > 0.0 {
            4.0 * omega / least.sin()
        } else {
            f64::INFINITY
        };
        Self {
            spine,
            supports,
            sides,
            turn,
        }
    }

    fn contacts(spine: &Spine, supports: [Surf; 2], sides: [f64; 2], t: f64) -> [Vector; 2] {
        let c = (spine.at)(t);
        [
            supports[0].gradient(c) * -sides[0],
            supports[1].gradient(c) * -sides[1],
        ]
    }

    fn gap(&self, q: Point, t: f64, r: f64) -> f64 {
        let [a, b] = Self::contacts(self.spine, self.supports, self.sides, t);
        lune_gap(q - (self.spine.at)(t), r, lune_normals(a, b))
    }
}

/// Whether a point of `patch` comes inside the lune of a ball of radius
/// `r` centred on the spine, by halving boxes of the patch's two parameters
/// and the run's under their Lipschitz bounds. Where the lune's gap is held
/// by its planes the point is within twice the radius of the centre, so the
/// gap moves with `t` at most the centre's speed and twice the radius times
/// the planes' turn.
fn patch_inside(patch: Patch, lunes: &Lunes, r: f64) -> Rolls {
    let (lu, lv) = patch.lipschitz();
    let lt = 2.0f64.mul_add(r * lunes.turn, lunes.spine.speed);
    if !lt.is_finite() {
        return Rolls::Unsure("the contacts' directions meet along the run".into());
    }
    let mut stack: Vec<[f64; 6]> = Vec::new();
    let (nu, nt) = (4, 96);
    for i in 0..nu {
        for j in 0..nu {
            for k in 0..nt {
                let (i, j, k) = (f64::from(i), f64::from(j), f64::from(k));
                let n = f64::from(nu);
                let m = f64::from(nt);
                stack.push([
                    i / n,
                    (i + 1.0) / n,
                    j / n,
                    (j + 1.0) / n,
                    TAU * k / m,
                    TAU * (k + 1.0) / m,
                ]);
            }
        }
    }
    let mut budget = 3_000_000_u32;
    while let Some([u0, u1, v0, v1, t0, t1]) = stack.pop() {
        budget -= 1;
        if budget == 0 {
            return Rolls::Unsure("a face's distance to the balls was not bounded".into());
        }
        let (u, v, t) = (
            f64::midpoint(u0, u1),
            f64::midpoint(v0, v1),
            f64::midpoint(t0, t1),
        );
        let q = patch.at(u, v);
        let gap = lunes.gap(q, t, r);
        if gap < -TOUCH {
            return Rolls::No(format!(
                "another face stands {:.3e} inside a ball at ({:.6}, {:.6}, {:.6})",
                -gap, q.x, q.y, q.z
            ));
        }
        let (su, sv, st) = (lu * (u1 - u0), lv * (v1 - v0), lt * (t1 - t0));
        if gap - (su + sv + st) / 2.0 >= -TOUCH {
            continue;
        }
        if su >= sv && su >= st {
            stack.push([u0, u, v0, v1, t0, t1]);
            stack.push([u, u1, v0, v1, t0, t1]);
        } else if sv >= st {
            stack.push([u0, u1, v0, v, t0, t1]);
            stack.push([u0, u1, v, v1, t0, t1]);
        } else {
            stack.push([u0, u1, v0, v1, t0, t]);
            stack.push([u0, u1, v0, v1, t, t1]);
        }
    }
    Rolls::Yes
}

/// A face's reach in space: how far a contact stands inside it.
type Reach3 = Box<dyn Fn(Point) -> f64>;

/// The spine of a corner for the given sides (the centre's signed distance
/// from each support over the radius) and radius.
type SpineOf = Box<dyn Fn([f64; 2], f64) -> Spine>;

/// A corner in space: two supports, the centres' path for any side and
/// radius, the faces' reach, the other faces.
struct Marched {
    supports: [Surf; 2],
    spine: SpineOf,
    reach: [Reach3; 2],
    others: Vec<Patch>,
    /// Planes of input faces other than the patches'.
    planes: Vec<(Point, Vector)>,
    /// An old round the kernel made between the two faces, by its radius.
    strip: Option<f64>,
    /// How far a point stands inside a strip fused on the faces (negative
    /// outside it): a contact there is off its face.
    cover: Option<Reach3>,
}

/// Stations the run is sampled at for the contacts' reach and the nearest
/// centre's seed.
const STATIONS: usize = 720;

impl Marched {
    fn station(t: usize) -> f64 {
        TAU * t as f64 / STATIONS as f64
    }

    /// How far a contact `q` stands inside face `i`: its reach, less how
    /// far it stands inside a strip lying on the face.
    fn margin(&self, i: usize, q: Point) -> f64 {
        let reach = (self.reach[i])(q);
        self.cover
            .as_ref()
            .map_or(reach, |cover| reach.min(-cover(q)))
    }

    fn rolls(&self, spine: &Spine, sides: [f64; 2], r: f64, solid: bool) -> Rolls {
        // Stations where a contact is off its face, and the last of them.
        let mut off = [0_usize; 2];
        let mut unsure = None;
        let mut seen = None;
        for k in 0..STATIONS {
            let c = (spine.at)(Self::station(k));
            for i in 0..2 {
                let s = self.supports[i];
                if !c.is_finite()
                    || (s.signed(c) - sides[i] * r).abs() > 1e-9 * (1.0 + c.to_vector().magnitude())
                {
                    return Rolls::Unsure("the centres' path left the offsets".into());
                }
                let q = s.foot(c);
                let margin = self.margin(i, q);
                if margin < -reading(r) {
                    off[i] += 1;
                    seen = Some((i, q));
                } else if margin <= reading(r) {
                    unsure = Some(format!(
                        "a contact touches the edge of face {}",
                        "ab".as_bytes()[i] as char
                    ));
                }
            }
        }
        if let Some((i, q)) = seen {
            let at = format!(
                "the ball's contact ({:.6}, {:.6}, {:.6}) is off face {}",
                q.x,
                q.y,
                q.z,
                "ab".as_bytes()[i] as char
            );
            return if off.contains(&STATIONS) {
                Rolls::No(format!("{at}, all round the run"))
            } else {
                Rolls::Capped(at)
            };
        }
        let mut found = unsure.map_or(Rolls::Yes, Rolls::Unsure);
        if !solid {
            return found;
        }
        if let Some(rho) = self.strip {
            found = found.and(Rolls::of_margin(r - rho, TOUCH, || {
                format!("an old round of radius {rho} stands inside the ball")
            }));
        }
        let lunes = Lunes::new(spine, self.supports, sides, r);
        for p in &self.others {
            found = found.and(patch_inside(*p, &lunes, r));
            if matches!(found, Rolls::No(_)) {
                return found;
            }
        }
        found
    }

    /// The ball's frame at `t`: centre, unit tangent, the tangent's turn
    /// per unit `t`, speed, and the unit directions to the two contacts.
    fn frame(
        &self,
        spine: &Spine,
        sides: [f64; 2],
        t: f64,
    ) -> (Point, Vector, Vector, f64, [Vector; 2]) {
        // Fourth-order central differences: their error, h^4 times the
        // path's fifth derivative, and their rounding, 1e-16 over h, both
        // stay near 1e-12 at h = 1e-3.
        let h = 1e-3;
        let slope = |f: &dyn Fn(f64) -> Vector, t: f64| {
            (f(t - 2.0 * h) - f(t + 2.0 * h) + (f(t + h) - f(t - h)) * 8.0) / (12.0 * h)
        };
        let at = |t: f64| (spine.at)(t).to_vector();
        let velocity = |t: f64| slope(&at, t);
        let tangent = |t: f64| {
            let v = velocity(t);
            v / v.magnitude()
        };
        let c = (spine.at)(t);
        let v = velocity(t);
        let tan = v / v.magnitude();
        let turn = slope(&tangent, t);
        let to = |i: usize| {
            let g = self.supports[i].gradient(c) * -sides[i];
            let g = g - tan * g.dot(tan);
            g / g.magnitude()
        };
        (c, tan, turn, v.magnitude(), [to(0), to(1)])
    }

    /// The volume between the balls' arcs and the two supports, over the
    /// run: each normal section a fan of rays from the centre between the
    /// contacts, from the ball out to the first support, weighted by the
    /// tube's Jacobian `s (|c'| - s <u, T'>)`.
    fn corner_volume(&self, spine: &Spine, sides: [f64; 2], r: f64) -> Option<f64> {
        let n = 384;
        let mut total = 0.0;
        for k in 0..n {
            let t = TAU * f64::from(k) / f64::from(n);
            let (c, _, turn, speed, [ua, ub]) = self.frame(spine, sides, t);
            let sweep = angle_between(ua, ub);
            let e2 = ub - ua * ub.dot(ua);
            let e2 = e2 / e2.magnitude();
            let ray = |phi: f64| ua * phi.cos() + e2 * phi.sin();
            let reach = |phi: f64, i: usize| {
                self.supports[i]
                    .hit(c, ray(phi), r)
                    .unwrap_or(f64::INFINITY)
            };
            // Where the first support met changes from a to b.
            let (mut lo, mut hi) = (0.0, sweep);
            for _ in 0..100 {
                let mid = f64::midpoint(lo, hi);
                if reach(mid, 0) <= reach(mid, 1) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let kink = f64::midpoint(lo, hi);
            let mut bad = false;
            let mut slice = |a: f64, b: f64, i: usize| {
                integrate(a, b, |phi| {
                    let s = reach(phi, i);
                    let bend = ray(phi).dot(turn);
                    if !s.is_finite() || s * bend >= speed {
                        bad = true;
                        return 0.0;
                    }
                    speed * (s * s - r * r) / 2.0 - bend * (s * s * s - r * r * r) / 3.0
                })
            };
            let section = slice(0.0, kink, 0) + slice(kink, sweep, 1);
            if bad {
                return None;
            }
            total += section;
        }
        Some(total * TAU / f64::from(n))
    }

    /// The nearest centre to `p` and its parameter.
    fn nearest(spine: &Spine, p: Point) -> f64 {
        let mut best = (f64::INFINITY, 0.0);
        for k in 0..STATIONS {
            let t = Self::station(k);
            let d = (spine.at)(t).distance(p);
            if d < best.0 {
                best = (d, t);
            }
        }
        let step = TAU / STATIONS as f64;
        let (mut a, mut b) = (best.1 - step, best.1 + step);
        let golden = (5.0_f64.sqrt() - 1.0) / 2.0;
        for _ in 0..80 {
            let (x, y) = (b - golden * (b - a), a + golden * (b - a));
            if (spine.at)(x).distance(p) < (spine.at)(y).distance(p) {
                b = y;
            } else {
                a = x;
            }
        }
        f64::midpoint(a, b)
    }

    fn off_round(&self, spine: &Spine, sides: [f64; 2], r: f64, p: Point) -> (f64, f64) {
        let t = Self::nearest(spine, p);
        let (c, _, _, _, [ua, ub]) = self.frame(spine, sides, t);
        let u = p - c;
        let outside = angle_between(u, ua) + angle_between(u, ub) - angle_between(ua, ub);
        (u.magnitude() - r, outside)
    }

    fn on_input(&self, p: Point) -> f64 {
        let mut near = self.supports[0]
            .signed(p)
            .abs()
            .min(self.supports[1].signed(p).abs());
        for patch in &self.others {
            let (o, n) = patch.plane();
            near = near.min((p - o).dot(n).abs());
        }
        for (o, n) in &self.planes {
            near = near.min((p - *o).dot(*n).abs());
        }
        near
    }
}

// --- a case ----------------------------------------------------------------

enum Geometry {
    Section(Section),
    Marched(Marched),
}

/// A generated solid, its two faces, and what the case knows of it.
struct Scene {
    model: Model,
    solid: Shape,
    faces: [Shape; 2],
    what: String,
    geometry: Geometry,
    /// The ball rolls inside the material (a convex corner, taken off) or
    /// outside it (a concave one, filled).
    convex: bool,
    /// The solid's volume without its strip.
    sharp: f64,
    /// The solid's own volume in closed form, where it has one.
    input: Option<f64>,
    /// Where the radius is drawn.
    radii: (f64, f64),
}

/// The face of `solid` holding `p`.
fn face_at(model: &Model, solid: &Shape, p: Point) -> Option<Shape> {
    let deflection = Deflection::with_chord(1e-3).ok()?;
    explore_unique(model, solid, ShapeType::Face)
        .ok()?
        .into_iter()
        .find(|f| {
            matches!(
                ogeom::algo::classify_on_face(model, f, p, deflection, T),
                Ok(ogeom::algo::Containment::In)
            )
        })
}

fn shared_edges(model: &Model, a: &Shape, b: &Shape) -> Vec<Shape> {
    let theirs = explore_unique(model, b, ShapeType::Edge).unwrap_or_default();
    explore_unique(model, a, ShapeType::Edge)
        .unwrap_or_default()
        .into_iter()
        .filter(|e| theirs.iter().any(|t| t.node() == e.node()))
        .collect()
}

/// A face's loosest stated tolerance: its own, its edges', its vertices'.
fn loosest(model: &Model, face: &Shape) -> f64 {
    let mut worst = model
        .node(face)
        .and_then(|n| n.data().as_face().map(|d| d.tolerance.get()))
        .unwrap_or(0.0);
    for e in explore_unique(model, face, ShapeType::Edge).unwrap_or_default() {
        if let Some(d) = model.node(&e).and_then(|n| n.data().as_edge()) {
            worst = worst.max(d.tolerance.get());
        }
    }
    for v in explore_unique(model, face, ShapeType::Vertex).unwrap_or_default() {
        if let Some(d) = model.node(&v).and_then(|n| n.data().as_vertex()) {
            worst = worst.max(d.tolerance.get());
        }
    }
    worst
}

fn input_error(what: &str, e: impl std::fmt::Display) -> Verdict {
    Verdict::of(Outcome::Refused, format!("input: {what}: {e}"))
}

/// The balls a case judges against, for one side and radius.
enum Balls<'a> {
    Section(&'a Section, Ball2),
    Marched(&'a Marched, Spine),
}

impl Balls<'_> {
    fn off_round(&self, sides: [f64; 2], r: f64, p: Point) -> (f64, f64) {
        match self {
            Self::Section(s, ball) => s.off_round(ball, p),
            Self::Marched(m, spine) => m.off_round(spine, sides, r, p),
        }
    }

    fn on_input(&self, p: Point) -> f64 {
        match self {
            Self::Section(s, _) => s.on_input(p),
            Self::Marched(m, _) => m.on_input(p),
        }
    }
}

/// Every face of `result` that lies on no input surface must lie on the
/// balls' envelope, inside the arcs between their contacts, within its
/// stated tolerance, and a `sheet` (which has no volume to say it) must
/// hold one. The faces' tolerance times area, summed, comes back: what the
/// new faces may move the volume by.
fn judge_faces(
    model: &Model,
    result: &Shape,
    balls: &Balls,
    (sides, r): ([f64; 2], f64),
    sheet: bool,
    what: &str,
) -> Result<f64, Verdict> {
    let deflection = Deflection::with_chord(r * 0.02).map_err(|e| input_error(what, e))?;
    let mut doubt = 0.0;
    let mut rounds = 0;
    for face in explore_unique(model, result, ShapeType::Face).unwrap_or_default() {
        // A face that will not mesh is left to the volume.
        let Ok(mesh) = ogeom::mesh::triangulate_face(model, &face, deflection, T) else {
            continue;
        };
        // The mesh's vertices, and the surface at each triangle's centre in
        // its chart: a band can mesh with every vertex on its boundary.
        let mut samples = mesh.positions.clone();
        if let (Some(data), Ok(placement)) = (
            model.node(&face).and_then(|n| n.data().as_face()),
            face.transform(model.datums()),
        ) && let Some(surface) = model.geometry().surface(data.surface)
            && mesh.parameters.len() == mesh.positions.len()
        {
            use ogeom::geom::Surface as _;
            for t in &mesh.triangles {
                let [a, b, c] = t.map(|i| mesh.parameters[i as usize]);
                let (u, v) = ((a.0 + b.0 + c.0) / 3.0, (a.1 + b.1 + c.1) / 3.0);
                if let Ok(p) = surface.point_at(u, v, T) {
                    samples.push(placement.apply(p));
                }
            }
        }
        let tol = loosest(model, &face);
        let scale = 1.0
            + mesh
                .positions
                .iter()
                .map(|p| p.to_vector().magnitude())
                .fold(0.0, f64::max);
        let slack = tol + 1e-9 * scale;
        if samples.iter().all(|p| balls.on_input(*p) <= slack) {
            continue;
        }
        rounds += 1;
        let (mut worst, mut wide, mut at) = (0.0_f64, 0.0_f64, Point::ORIGIN);
        for p in &samples {
            let (off, outside) = balls.off_round(sides, r, *p);
            if off.abs() > worst {
                worst = off.abs();
                at = *p;
            }
            wide = wide.max(outside);
        }
        if worst > slack {
            return Err(Verdict::of(
                Outcome::Wrong,
                format!(
                    "{what}: a new face stands {worst:.3e} off the balls' envelope at \
                     ({:.6}, {:.6}, {:.6}), against its tolerance {tol:.1e}",
                    at.x, at.y, at.z
                ),
            ));
        }
        let angle = 2.0 * slack / r + 1e-9;
        if wide > angle {
            return Err(Verdict::of(
                Outcome::Wrong,
                format!("{what}: a new face runs {wide:.3e} rad past the ball's contacts"),
            ));
        }
        let area = ogeom::algo::surface_properties(
            model,
            &face,
            Deflection::with_chord(super::CHORD).unwrap(),
            T,
        )
        .map_or(0.0, |p| p.mass);
        doubt += tol * area;
    }
    if rounds == 0 && sheet {
        return Err(Verdict::of(
            Outcome::Wrong,
            format!("{what}: no new face; the result is the input's faces alone"),
        ));
    }
    Ok(doubt)
}

/// One case of a family. Every note but an `ok`'s leads with what the
/// operation said or what the oracle found, and ends with the case.
pub(super) fn faceblend_case(family: &str, rng: &mut Rng) -> Verdict {
    let mut described = String::new();
    let mut verdict = run_case(family, rng, &mut described);
    if verdict.outcome != Outcome::Ok && !described.is_empty() {
        verdict.note = format!("{}; case: {described}", verdict.note);
    }
    verdict
}

fn run_case(family: &str, rng: &mut Rng, described: &mut String) -> Verdict {
    let scene = match family {
        "planes" => planes(rng, false),
        "notch" => planes(rng, true),
        "dshaft" => dshaft(rng),
        "rod" => rod(rng),
        "drums" => drums(rng),
        "rim" => rim(rng),
        "bore" => bore(rng),
        "boss" => boss(rng),
        "frustum" => frustum(rng),
        "conefoot" => cone_foot(rng),
        "leaning" => leaning(rng),
        "tee" => tee(rng),
        "slanted" => slanted(rng),
        _ => return Verdict::of(Outcome::Refused, format!("no family {family}")),
    };
    let mut scene = match scene {
        Ok(s) => s,
        Err(v) => return v,
    };
    described.clone_from(&scene.what);
    // The input itself, against its closed form.
    let start = match measure(&scene.model, &scene.solid) {
        Ok(m) => m,
        Err(e) => return input_error("volume", e),
    };
    if let Some(want) = scene.input {
        let slack = allowance(want, &[&start]);
        if (start.volume - want).abs() > slack {
            return Verdict::of(
                Outcome::Wrong,
                format!(
                    "input volume {:.9} against {want:.9}, off {:.3e}",
                    start.volume,
                    start.volume - want
                ),
            );
        }
    }

    let sheet = rng.chance(0.25);
    // Trimming a cone standing on a plane returns a cone face whose
    // pcurves miss their curves; those sheets are asked untrimmed.
    let trim = rng.chance(0.5) && !matches!(family, "frustum" | "conefoot");
    let blend_sides = if scene.convex {
        [-1.0, -1.0]
    } else {
        [1.0, 1.0]
    };
    let sides = if sheet { [1.0, 1.0] } else { blend_sides };
    let r = draw_radius(&scene, sides, rng);
    let what = format!(
        "{} r {r}",
        if sheet {
            if trim {
                "fillet_faces trimmed"
            } else {
                "fillet_faces"
            }
        } else {
            "blend_faces"
        }
    );

    let (rolls, balls) = match &scene.geometry {
        Geometry::Section(s) => {
            let rolls = s.rolls(sides, r, !sheet);
            let Some(ball) = s.ball(sides, r) else {
                return judge_refusal(&mut scene, sheet, trim, r, &what, &rolls);
            };
            (rolls, Balls::Section(s, ball))
        }
        Geometry::Marched(m) => {
            let spine = (m.spine)(sides, r);
            (m.rolls(&spine, sides, r, !sheet), Balls::Marched(m, spine))
        }
    };
    if let Rolls::No(_) = rolls {
        return judge_refusal(&mut scene, sheet, trim, r, &what, &rolls);
    }

    let [a, b] = scene.faces.clone();
    let result = if sheet {
        ogeom::fillet::fillet_faces(&mut scene.model, &a, &b, r, trim, T)
    } else {
        ogeom::fillet::blend_faces(&mut scene.model, &scene.solid, &a, &b, r, T)
    };
    let note = match &rolls {
        Rolls::Unsure(why) => format!(" (undecided whether a ball rolls: {why})"),
        Rolls::Capped(why) => format!(" (a run the round may cap: {why})"),
        _ => " (a ball rolls there)".to_string(),
    };
    if sheet {
        let built = match result {
            Ok(b) => b,
            Err(e) => return Verdict::of(Outcome::Refused, format!("{what}: {e}{note}")),
        };
        let diagnosis = match ogeom::algo::check(&scene.model, &built.shape, T) {
            Ok(d) => d,
            Err(e) => {
                return Verdict::of(
                    Outcome::Invalid,
                    format!("{what}: check cannot read it: {e}"),
                );
            }
        };
        if let Some(first) = sheet_problem(&diagnosis) {
            return Verdict::of(Outcome::Invalid, format!("{what}: invalid: {first}"));
        }
        if let Rolls::Capped(_) = rolls {
            return Verdict::ok();
        }
        return match judge_faces(&scene.model, &built.shape, &balls, (sides, r), true, &what) {
            Ok(_) => Verdict::ok(),
            Err(v) => v,
        };
    }
    let (shape, m) = match judged(&scene.model, &what, result) {
        Ok(r) => r,
        Err(mut v) => {
            if v.outcome == Outcome::Refused {
                v.note.push_str(&note);
            }
            return v;
        }
    };
    if let Rolls::Capped(_) = rolls {
        return Verdict::ok().read_from(&[("result", &m)]);
    }
    let doubt = match judge_faces(&scene.model, &shape, &balls, (sides, r), false, &what) {
        Ok(d) => d,
        Err(v) => return v,
    };
    let corner = match &balls {
        Balls::Section(s, ball) => Some(s.corner_volume(ball, r)),
        Balls::Marched(mm, spine) => mm.corner_volume(spine, sides, r),
    };
    let Some(corner) = corner else {
        return Verdict::ok().read_from(&[("result", &m)]);
    };
    let want = if scene.convex {
        scene.sharp - corner
    } else {
        scene.sharp + corner
    };
    let slack = allowance(want, &[&m, &start]) + doubt + corner * 1e-9;
    if (m.volume - want).abs() > slack {
        return Verdict::of(
            Outcome::Wrong,
            format!(
                "{what}: volume {:.9} against {want:.9} (the sharp {:.9} {} the corner {corner:.9}), \
                 off {:.3e} against {slack:.3e}",
                m.volume,
                scene.sharp,
                if scene.convex { "less" } else { "with" },
                m.volume - want
            ),
        );
    }
    Verdict::ok().read_from(&[("input", &start), ("result", &m)])
}

/// A case where no ball rolls: a refusal is right, and anything built is
/// wrong unless it is invalid (then it is that).
fn judge_refusal(
    scene: &mut Scene,
    sheet: bool,
    trim: bool,
    r: f64,
    what: &str,
    rolls: &Rolls,
) -> Verdict {
    let why = match rolls {
        Rolls::No(why) | Rolls::Unsure(why) | Rolls::Capped(why) => why.clone(),
        Rolls::Yes => "no ball".to_string(),
    };
    let [a, b] = scene.faces.clone();
    let result = if sheet {
        ogeom::fillet::fillet_faces(&mut scene.model, &a, &b, r, trim, T)
    } else {
        ogeom::fillet::blend_faces(&mut scene.model, &scene.solid, &a, &b, r, T)
    };
    let Ok(built) = result else {
        return Verdict::ok();
    };
    let problem = ogeom::algo::check(&scene.model, &built.shape, T)
        .ok()
        .and_then(|d| {
            if sheet {
                sheet_problem(&d)
            } else {
                d.problems.first().map(|p| p.what.clone())
            }
        });
    match problem {
        Some(first) => Verdict::of(
            Outcome::Invalid,
            format!("{what}: invalid: {first} (and no ball rolls: {why})"),
        ),
        _ => Verdict::of(
            Outcome::Wrong,
            format!("{what}: built a valid shape where no ball rolls: {why}"),
        ),
    }
}

/// The first problem `check` finds with a sheet: an open shell is what
/// `fillet_faces` returns, so the note that it is open is no problem.
fn sheet_problem(d: &ogeom::algo::Diagnosis) -> Option<String> {
    d.problems
        .iter()
        .find(|p| !p.what.contains("so the shell is open along them"))
        .map(|p| p.what.clone())
}

/// A radius from the scene's range; a third of the time one just either
/// side of where a ball starts or stops rolling, if the range holds one.
fn draw_radius(scene: &Scene, sides: [f64; 2], rng: &mut Rng) -> f64 {
    let (lo, hi) = scene.radii;
    let feasible = |r: f64| -> bool {
        match &scene.geometry {
            Geometry::Section(s) => matches!(s.rolls(sides, r, true), Rolls::Yes),
            Geometry::Marched(m) => {
                let spine = (m.spine)(sides, r);
                // The contacts and the old round alone: the cheap part.
                let mut ok = true;
                for k in (0..STATIONS).step_by(4) {
                    let c = (spine.at)(Marched::station(k));
                    for i in 0..2 {
                        ok &= m.margin(i, m.supports[i].foot(c)) > 0.0;
                    }
                }
                ok && m.strip.is_none_or(|rho| rho < r)
            }
        }
    };
    let fits = |r: f64| match &scene.geometry {
        Geometry::Section(s) => s
            .within
            .is_none_or(|w| s.ball(sides, r).is_some_and(|b| b.centre.x + r < w)),
        Geometry::Marched(_) => true,
    };
    let mut r = rng.range(lo, hi);
    if rng.chance(0.35) {
        let n = 48;
        let grid: Vec<f64> = (0..=n)
            .map(|i| lo + (hi - lo) * f64::from(i) / f64::from(n))
            .collect();
        let marks: Vec<bool> = grid.iter().map(|r| feasible(*r)).collect();
        let changes: Vec<usize> = (0..n as usize)
            .filter(|i| marks[*i] != marks[i + 1])
            .collect();
        if !changes.is_empty() {
            let i = *rng.pick(&changes);
            let (mut a, mut b) = (grid[i], grid[i + 1]);
            for _ in 0..60 {
                let mid = f64::midpoint(a, b);
                if feasible(mid) == marks[i] {
                    a = mid;
                } else {
                    b = mid;
                }
            }
            let edge = f64::midpoint(a, b);
            let nudge = 10f64.powf(rng.range(-5.0, -2.0));
            r = edge
                * if rng.chance(0.5) {
                    1.0 + nudge
                } else {
                    1.0 - nudge
                };
        }
    }
    for _ in 0..40 {
        if fits(r) {
            break;
        }
        r *= 0.85;
    }
    r
}

// --- generators ------------------------------------------------------------

/// A frame in a line: its origin, its axis and its `x`.
fn brief(f: &Frame) -> String {
    let (o, z, x) = (f.origin(), f.z().vector(), f.x().vector());
    format!(
        "at ({:.6}, {:.6}, {:.6}) along ({:.6}, {:.6}, {:.6}) x ({:.6}, {:.6}, {:.6})",
        o.x, o.y, o.z, z.x, z.y, z.z, x.x, x.y, x.z
    )
}

fn random_frame(rng: &mut Rng) -> Result<Frame, Verdict> {
    let at = Point::new(
        rng.range(-3.0, 3.0),
        rng.range(-3.0, 3.0),
        rng.range(-3.0, 3.0),
    );
    let z = Vector::new(
        rng.range(-1.0, 1.0),
        rng.range(-1.0, 1.0),
        rng.range(-1.0, 1.0),
    );
    let z = Direction::new(z, T).map_err(|_| Verdict::of(Outcome::Refused, "degenerate axis"))?;
    let x = Direction::new(across(z.vector(), rng), T)
        .map_err(|_| Verdict::of(Outcome::Refused, "degenerate axis"))?;
    Frame::new(at, z, x, T).map_err(|e| input_error("frame", e))
}

/// The frame at `at` with axis `z` and `x` toward `x`.
fn frame_in(at: Point, z: Vector, x: Vector) -> Result<Frame, Verdict> {
    let z = Direction::new(z, T).map_err(|e| input_error("frame", e))?;
    let x = Direction::new(x, T).map_err(|e| input_error("frame", e))?;
    Frame::new(at, z, x, T).map_err(|e| input_error("frame", e))
}

/// A world vector from `frame`'s own coordinates.
fn dir(frame: &Frame, x: f64, y: f64, z: f64) -> Vector {
    frame.x().vector() * x + frame.y().vector() * y + frame.z().vector() * z
}

fn place(frame: &Frame, x: f64, y: f64, z: f64) -> Point {
    frame.origin() + dir(frame, x, y, z)
}

/// The prism of a section polygon in `frame`'s `xy`, from `z0` along `z`
/// by `length`.
fn prism(
    model: &mut Model,
    frame: &Frame,
    section: &[P2],
    z0: f64,
    length: f64,
) -> ogeom::core::OgeomResult<Shape> {
    let mut section = section.to_vec();
    let area: f64 = (0..section.len())
        .map(|i| {
            let (p, q) = (section[i], section[(i + 1) % section.len()]);
            p.x * q.y - q.x * p.y
        })
        .sum();
    if area < 0.0 {
        section.reverse();
    }
    let points: Vec<Point> = section.iter().map(|q| place(frame, q.x, q.y, z0)).collect();
    let wire = ogeom::algo::make_polygon(model, &points, true, T)?.shape;
    let base = Frame::new(place(frame, 0.0, 0.0, z0), frame.z(), frame.x(), T)?;
    let plane = ogeom::geom::PlaneSurface::new(ogeom::math::Plane::new(base));
    let face = ogeom::algo::make_face(model, plane.into(), &[wire], T)?.shape;
    Ok(ogeom::algo::make_prism(model, &face, frame.z().vector() * length, T)?.shape)
}

fn polygon_area(section: &[P2]) -> f64 {
    (0..section.len())
        .map(|i| {
            let (p, q) = (section[i], section[(i + 1) % section.len()]);
            p.x * q.y - q.x * p.y
        })
        .sum::<f64>()
        .abs()
        / 2.0
}

/// The shape a boolean or a blend built, or why not as a refusal of the
/// input.
fn built(what: &str, r: ogeom::core::OgeomResult<ogeom::algo::Built>) -> Result<Shape, Verdict> {
    r.map(|b| b.shape).map_err(|e| input_error(what, e))
}

fn find_face(model: &Model, solid: &Shape, p: Point, what: &str) -> Result<Shape, Verdict> {
    face_at(model, solid, p).ok_or_else(|| {
        Verdict::of(
            Outcome::Refused,
            format!(
                "input: {what}: no face at ({:.6}, {:.6}, {:.6})",
                p.x, p.y, p.z
            ),
        )
    })
}

/// Two planes meeting at the crease `V` of a prism, at a random angle:
/// a parallelogram's corner (convex) or a notch's (concave), with a strip
/// between the faces or not: a chamfer or an old round by the kernel, or a
/// wedge cut off or filled in by a prism tool, plane or B-spline, over the
/// whole run or part of it.
fn planes(rng: &mut Rng, concave: bool) -> Result<Scene, Verdict> {
    let frame = random_frame(rng)?;
    let mut model = Model::new();
    let angle = rng.range(25.0, 155.0).to_radians();
    let (la, lb, length) = (
        rng.range(2.0, 6.0),
        rng.range(2.0, 6.0),
        rng.range(2.0, 6.0),
    );
    let v = P2::new(0.0, 0.0);
    let (da, db) = (V2::new(1.0, 0.0), V2::new(angle.cos(), angle.sin()));
    let (a_end, b_end) = (v + da * la, v + db * lb);
    let h = rng.range(1.0, 3.0);
    let (section, inner, mut others) = if concave {
        // Free space in the wedge between da and db; the material a band
        // of depth h behind both faces.
        let (na, nb) = (V2::new(0.0, 1.0), V2::new(angle.sin(), -angle.cos()));
        let k = v - (na + nb) * (h / (1.0 - angle.cos()));
        let (a2, b2) = (a_end - na * h, b_end - nb * h);
        let poly = vec![v, a_end, a2, k, b2, b_end];
        let others = vec![
            Piece::Seg(a_end, a2),
            Piece::Seg(a2, k),
            Piece::Seg(k, b2),
            Piece::Seg(b2, b_end),
        ];
        (poly, v + (k - v) * 0.5, others)
    } else {
        let far = a_end + (b_end - v);
        let poly = vec![v, a_end, far, b_end];
        let others = vec![Piece::Seg(a_end, far), Piece::Seg(far, b_end)];
        (poly, v + (far - v) * 0.5, others)
    };
    let supports = [Curve2::line(v, a_end, inner), Curve2::line(v, b_end, inner)];
    let sharp_area = polygon_area(&section);
    let solid =
        prism(&mut model, &frame, &section, 0.0, length).map_err(|e| input_error("prism", e))?;

    let kinds: &[&str] = if concave {
        &[
            "none",
            "round",
            "fill",
            "fill-part",
            "spline",
            "spline-part",
        ]
    } else {
        &[
            "none",
            "chamfer",
            "round",
            "cut",
            "cut-part",
            "spline",
            "spline-part",
        ]
    };
    let strip = *rng.pick(kinds);
    let short = la.min(lb);
    let (d1, d2) = (short * rng.range(0.05, 0.4), short * rng.range(0.05, 0.4));
    let (d1, d2) = if strip == "chamfer" && rng.chance(0.4) {
        (d1, d1)
    } else {
        (d1, d2)
    };
    let setback = short * rng.range(0.05, 0.4);
    let rho = setback * (angle / 2.0).tan();
    let (c1, c2) = (v + da * d1, v + db * d2);
    let span = {
        let z0 = length * rng.range(0.1, 0.45);
        (z0, z0 + length * rng.range(0.05, 0.4))
    };
    let partial = strip.ends_with("part");
    let mid = |model: &Model, solid: &Shape, lo: f64, lo_b: f64| -> Result<[Shape; 2], Verdict> {
        let z = if partial {
            f64::midpoint(span.0, span.1)
        } else {
            length / 2.0
        };
        let pa = v + da * f64::midpoint(lo, la);
        let pb = v + db * f64::midpoint(lo_b, lb);
        Ok([
            find_face(model, solid, place(&frame, pa.x, pa.y, z), "face a")?,
            find_face(model, solid, place(&frame, pb.x, pb.y, z), "face b")?,
        ])
    };
    let mut planes = vec![
        (frame.origin(), frame.z().vector()),
        (place(&frame, 0.0, 0.0, length), frame.z().vector()),
    ];
    let (solid, lo, input, what) = match strip {
        "none" => (solid, (0.0, 0.0), sharp_area * length, String::new()),
        "chamfer" => {
            let [fa, fb] = mid(&model, &solid, 0.0, 0.0)?;
            let edge = shared_edges(&model, &fa, &fb);
            if edge.len() != 1 {
                return Err(Verdict::of(
                    Outcome::Refused,
                    "input: the corner is not one edge",
                ));
            }
            let s = built(
                "chamfer",
                ogeom::fillet::chamfer_edge_distances(&mut model, &solid, &edge[0], &fa, d1, d2, T),
            )?;
            others.push(Piece::Seg(c1, c2));
            let cut = d1 * d2 * angle.sin() / 2.0;
            (
                s,
                (d1, d2),
                (sharp_area - cut) * length,
                format!(", chamfered {d1} by {d2}"),
            )
        }
        "round" => {
            let [fa, fb] = mid(&model, &solid, 0.0, 0.0)?;
            let edge = shared_edges(&model, &fa, &fb);
            if edge.len() != 1 {
                return Err(Verdict::of(
                    Outcome::Refused,
                    "input: the corner is not one edge",
                ));
            }
            let s = built(
                "round",
                ogeom::fillet::fillet_edge(&mut model, &solid, &edge[0], rho, T),
            )?;
            let side = if concave { 1.0 } else { -1.0 };
            let centre = meet(
                supports[0].offset(side * rho),
                supports[1].offset(side * rho),
                v,
            )
            .unwrap();
            let (ta, tb) = (supports[0].foot(centre), supports[1].foot(centre));
            others.push(Piece::arc(centre, rho, ta, tb));
            let corner = rho * rho * (1.0 / (angle / 2.0).tan() - (PI - angle) / 2.0);
            let input = if concave {
                sharp_area + corner
            } else {
                sharp_area - corner
            };
            (
                s,
                (setback, setback),
                input * length,
                format!(", rounded {rho}"),
            )
        }
        _ => {
            // A wedge tool past the chamfer line c1 c2, out beyond the
            // corner (convex: cut) or deep behind it (concave: fused).
            let k = 0.25;
            let (p1, p2) = (c1 + (c1 - c2) * k, c2 + (c2 - c1) * k);
            let outward = match (supports[0], supports[1]) {
                (Curve2::Line { out: na, .. }, Curve2::Line { out: nb, .. }) => na + nb,
                _ => unreachable!(),
            };
            let p3 = if concave {
                // Seven tenths of the way to the notch's inner corner.
                v - outward * (0.7 * h / (1.0 - angle.cos()))
            } else {
                v + outward * (3.0 * la.max(lb))
            };
            let inside_tri = |q: P2| {
                let s = |a: P2, b: P2| (b - a).cross(q - a);
                let (x, y, z) = (s(p1, p2), s(p2, p3), s(p3, p1));
                (x > 0.0 && y > 0.0 && z > 0.0) || (x < 0.0 && y < 0.0 && z < 0.0)
            };
            if !inside_tri(v) {
                return Err(Verdict::of(
                    Outcome::Refused,
                    "input: the wedge tool misses the corner",
                ));
            }
            // A cut tool runs past the caps; a fused one stops at them, or it
            // would add material beyond the part.
            let (z0, extent) = if partial {
                (span.0, span.1 - span.0)
            } else if concave {
                (0.0, length)
            } else {
                (-1.0, length + 2.0)
            };
            let mut tool = prism(&mut model, &frame, &[p1, p2, p3], z0, extent)
                .map_err(|e| input_error("wedge", e))?;
            if strip.starts_with("spline") {
                tool = built("to_nurbs", ogeom::algo::to_nurbs(&mut model, &tool, T))?;
            }
            let s = if concave {
                built("fuse", ogeom::boolean::fuse(&mut model, &solid, &tool, T))?
            } else {
                built("cut", ogeom::boolean::cut(&mut model, &solid, &tool, T))?
            };
            others.push(Piece::Seg(c1, c2));
            let wedge = d1 * d2 * angle.sin() / 2.0;
            let run = if partial { span.1 - span.0 } else { length };
            if partial {
                planes.push((place(&frame, 0.0, 0.0, span.0), frame.z().vector()));
                planes.push((place(&frame, 0.0, 0.0, span.1), frame.z().vector()));
            }
            let input = if concave {
                sharp_area * length + wedge * run
            } else {
                sharp_area * length - wedge * run
            };
            (
                s,
                (d1, d2),
                input,
                format!(
                    ", a {} wedge {d1} by {d2}{}",
                    if strip.starts_with("spline") {
                        "B-spline"
                    } else {
                        "plane"
                    },
                    if partial {
                        format!(" over z {:.4}..{:.4}", span.0, span.1)
                    } else {
                        String::new()
                    }
                ),
            )
        }
    };
    let faces = mid(&model, &solid, lo.0, lo.1)?;
    // A strip over part of the run takes the faces back to its edges over
    // that part, so a contact must stand past them all along.
    let reach_lo = lo;
    let mut carriers: Vec<Piece> = others.clone();
    carriers.push(Piece::Seg(v, a_end));
    carriers.push(Piece::Seg(v, b_end));
    let what = format!(
        "{} prism {la:.4} by {lb:.4} at {:.2} deg, {length:.4} long, frame {}{what}",
        if concave { "notch" } else { "corner" },
        angle.to_degrees(),
        brief(&frame),
    );
    let top = short * (angle / 2.0).tan();
    Ok(Scene {
        model,
        solid,
        faces,
        what,
        geometry: Geometry::Section(Section {
            sweep: Sweep::Prism { frame, length },
            supports,
            crease: v,
            reach: [along(v, da, reach_lo.0, la), along(v, db, reach_lo.1, lb)],
            others,
            mirror: None,
            carriers,
            planes,
            angle: Some(angle),
            within: None,
        }),
        convex: !concave,
        sharp: sharp_area * length,
        input: Some(input),
        radii: (top * 0.03, top * 1.1),
    })
}

/// The section's ball of radius `rho` in its corner, for an old round: its
/// arc, and each tangent point's distance from the crease along its face.
fn old_round(supports: [Curve2; 2], crease: P2, side: f64, rho: f64) -> Option<(Piece, [P2; 2])> {
    let centre = meet(
        supports[0].offset(side * rho),
        supports[1].offset(side * rho),
        crease,
    )?;
    let t = [supports[0].foot(centre), supports[1].foot(centre)];
    Some((Piece::arc(centre, rho, t[0], t[1]), t))
}

/// A drum along `z` cut flat along its length, a plane parallel to its
/// axis: two convex creases between the flat and the wall, both blended.
fn dshaft(rng: &mut Rng) -> Result<Scene, Verdict> {
    let frame = random_frame(rng)?;
    let mut model = Model::new();
    let (radius, length) = (rng.range(1.5, 4.0), rng.range(2.0, 6.0));
    let e = radius * rng.range(-0.6, 0.6);
    let xv = (radius * radius - e * e).sqrt();
    let v = P2::new(xv, e);
    let drum = built(
        "drum",
        ogeom::algo::make_cylinder(
            &mut model,
            frame_in(frame.origin(), frame.z().vector(), -frame.y().vector())?,
            radius,
            length,
            T,
        ),
    )?;
    let lid = frame_in(
        place(&frame, -radius - 1.0, e, -1.0),
        frame.z().vector(),
        frame.x().vector(),
    )?;
    let lid = built(
        "box",
        ogeom::algo::make_box(
            &mut model,
            lid,
            (2.0 * radius + 2.0, radius + 2.0, length + 2.0),
            T,
        ),
    )?;
    let mut solid = built("cut", ogeom::boolean::cut(&mut model, &drum, &lid, T))?;
    let inner = P2::new(0.0, (e - radius) / 2.0);
    let supports = [
        Curve2::line(P2::new(-xv, e), v, inner),
        Curve2::Circle {
            centre: P2::new(0.0, 0.0),
            radius,
            inside: true,
        },
    ];
    let segment = radius * radius * (e / radius).acos() - e * xv;
    let sharp = (PI * radius * radius - segment) * length;
    let theta_v = e.atan2(xv);
    let mut lo = [0.0, 0.0];
    let mut others = Vec::new();
    let mut input = sharp;
    let mut what = String::new();
    let [fa, fb] = [
        find_face(&model, &solid, place(&frame, 0.0, e, length / 2.0), "flat")?,
        find_face(
            &model,
            &solid,
            place(&frame, 0.0, -radius, length / 2.0),
            "wall",
        )?,
    ];
    if rng.chance(0.5) {
        let rho = radius * rng.range(0.03, 0.4);
        let edges = shared_edges(&model, &fa, &fb);
        solid = built(
            "round",
            ogeom::fillet::fillet_edges(&mut model, &solid, &edges, rho, T),
        )?;
        if let Some((arc, t)) = old_round(supports, v, -1.0, rho) {
            let (area, _) = region(&[supports[0].piece(v, t[0]), arc, supports[1].piece(t[1], v)]);
            input -= 2.0 * area.abs() * length;
            lo = [xv - t[0].x, (theta_v - t[1].y.atan2(t[1].x)) * radius];
            others.push(arc);
        }
        what = format!(", creases rounded {rho}");
    }
    let faces = [
        find_face(&model, &solid, place(&frame, 0.0, e, length / 2.0), "flat")?,
        find_face(
            &model,
            &solid,
            place(&frame, 0.0, -radius, length / 2.0),
            "wall",
        )?,
    ];
    let mut carriers = others.clone();
    carriers.push(Piece::Seg(P2::new(-xv, e), v));
    carriers.push(Piece::Arc {
        centre: P2::new(0.0, 0.0),
        radius,
        from: 0.0,
        to: TAU,
    });
    Ok(Scene {
        model,
        solid,
        faces,
        what: format!(
            "D-shaft r {radius:.4} flat at {e:.4}, {length:.4} long, frame {}{what}",
            brief(&frame)
        ),
        geometry: Geometry::Section(Section {
            sweep: Sweep::Prism { frame, length },
            supports,
            crease: v,
            reach: [
                along(v, V2::new(-1.0, 0.0), lo[0], xv),
                round_from(
                    P2::new(0.0, 0.0),
                    radius,
                    v,
                    -1.0,
                    lo[1],
                    (theta_v + PI / 2.0) * radius,
                ),
            ],
            others,
            mirror: Some(|p| P2::new(-p.x, p.y)),
            carriers,
            planes: vec![
                (frame.origin(), frame.z().vector()),
                (place(&frame, 0.0, 0.0, length), frame.z().vector()),
            ],
            angle: None,
            within: None,
        }),
        convex: true,
        sharp,
        input: Some(input),
        radii: (radius * 0.02, radius * 0.7),
    })
}

/// A rod along `z` half sunk in a block: the block's top meets its wall
/// along two concave creases; the wall and the top beside one crease are
/// blended.
fn rod(rng: &mut Rng) -> Result<Scene, Verdict> {
    let frame = random_frame(rng)?;
    let mut model = Model::new();
    let (radius, length) = (rng.range(1.0, 3.0), rng.range(2.0, 6.0));
    let yc = radius * rng.range(-0.6, 0.6);
    let width = 2.0 * radius + rng.range(3.0, 8.0);
    let depth = (radius - yc) + rng.range(0.5, 2.0);
    let block = frame_in(
        place(&frame, -width / 2.0, -depth, 0.0),
        frame.z().vector(),
        frame.x().vector(),
    )?;
    let block = built(
        "block",
        ogeom::algo::make_box(&mut model, block, (width, depth, length), T),
    )?;
    let axis = frame_in(
        place(&frame, 0.0, yc, 0.0),
        frame.z().vector(),
        -frame.y().vector(),
    )?;
    let bar = built(
        "rod",
        ogeom::algo::make_cylinder(&mut model, axis, radius, length, T),
    )?;
    let mut solid = built("fuse", ogeom::boolean::fuse(&mut model, &block, &bar, T))?;
    let xv = (radius * radius - yc * yc).sqrt();
    let v = P2::new(xv, 0.0);
    let centre = P2::new(0.0, yc);
    let supports = [
        Curve2::Circle {
            centre,
            radius,
            inside: true,
        },
        Curve2::line(
            P2::new(0.0, 0.0),
            P2::new(1.0, 0.0),
            P2::new(0.0, -depth / 2.0),
        ),
    ];
    let d = -yc;
    let cap = radius * radius * (d / radius).acos() - d * (radius * radius - d * d).sqrt();
    let sharp = (width * depth + cap) * length;
    let theta_v = (-yc).atan2(xv);
    let wall_at = |model: &Model, solid: &Shape| {
        find_face(
            model,
            solid,
            place(&frame, 0.0, yc + radius, length / 2.0),
            "rod wall",
        )
    };
    let top_at = |model: &Model, solid: &Shape| {
        find_face(
            model,
            solid,
            place(&frame, (xv + width / 2.0) / 2.0, 0.0, length / 2.0),
            "top",
        )
    };
    let mut lo = [0.0, 0.0];
    let mut others = vec![Piece::Seg(
        P2::new(width / 2.0, 0.0),
        P2::new(width / 2.0, -depth),
    )];
    let mut input = sharp;
    let mut what = String::new();
    if rng.chance(0.5) {
        let rho = radius * rng.range(0.03, 0.5);
        let edges = shared_edges(&model, &wall_at(&model, &solid)?, &top_at(&model, &solid)?);
        solid = built(
            "round",
            ogeom::fillet::fillet_edges(&mut model, &solid, &edges, rho, T),
        )?;
        if let Some((arc, t)) = old_round(supports, v, 1.0, rho) {
            let (area, _) = region(&[supports[0].piece(v, t[0]), arc, supports[1].piece(t[1], v)]);
            input += area.abs() * length;
            lo = [
                (wrap((t[0].y - yc).atan2(t[0].x) - theta_v)) * radius,
                t[1].x - xv,
            ];
            others.push(arc);
        }
        what = format!(", crease rounded {rho}");
    }
    let faces = [wall_at(&model, &solid)?, top_at(&model, &solid)?];
    let mut carriers = others.clone();
    carriers.push(Piece::Seg(P2::new(-width / 2.0, 0.0), v));
    carriers.push(Piece::Seg(
        P2::new(-width / 2.0, 0.0),
        P2::new(-width / 2.0, -depth),
    ));
    carriers.push(Piece::Seg(
        P2::new(-width / 2.0, -depth),
        P2::new(width / 2.0, -depth),
    ));
    carriers.push(Piece::Arc {
        centre,
        radius,
        from: 0.0,
        to: TAU,
    });
    Ok(Scene {
        model,
        solid,
        faces,
        what: format!(
            "rod r {radius:.4} at {yc:.4} in a block {width:.4} by {depth:.4}, {length:.4} long, frame {}{what}",
            brief(&frame)
        ),
        geometry: Geometry::Section(Section {
            sweep: Sweep::Prism { frame, length },
            supports,
            crease: v,
            reach: [
                round_from(centre, radius, v, 1.0, lo[0], (PI - 2.0 * theta_v) * radius),
                along(v, V2::new(1.0, 0.0), lo[1], width / 2.0 - xv),
            ],
            others,
            mirror: None,
            carriers,
            planes: vec![
                (frame.origin(), frame.z().vector()),
                (place(&frame, 0.0, 0.0, length), frame.z().vector()),
            ],
            angle: None,
            within: None,
        }),
        convex: false,
        sharp,
        input: Some(input),
        radii: (radius * 0.03, radius * 1.2),
    })
}

/// Two parallel drums fused side by side: their walls meet along two
/// concave creases, both blended.
fn drums(rng: &mut Rng) -> Result<Scene, Verdict> {
    let frame = random_frame(rng)?;
    let mut model = Model::new();
    let r1 = rng.range(1.5, 4.0);
    let r2 = r1 * rng.range(0.4, 1.0);
    let gap = rng.range((r1 - r2) + 0.25 * r2, r1 + r2 - 0.25 * r2);
    let length = rng.range(2.0, 6.0);
    let one = frame_in(frame.origin(), frame.z().vector(), -frame.x().vector())?;
    let two = frame_in(
        place(&frame, gap, 0.0, 0.0),
        frame.z().vector(),
        frame.x().vector(),
    )?;
    let one = built(
        "drum",
        ogeom::algo::make_cylinder(&mut model, one, r1, length, T),
    )?;
    let two = built(
        "drum",
        ogeom::algo::make_cylinder(&mut model, two, r2, length, T),
    )?;
    let mut solid = built("fuse", ogeom::boolean::fuse(&mut model, &one, &two, T))?;
    let (c1, c2) = (P2::new(0.0, 0.0), P2::new(gap, 0.0));
    let x = (gap * gap + r1 * r1 - r2 * r2) / (2.0 * gap);
    let v = P2::new(x, (r1 * r1 - x * x).sqrt());
    let supports = [
        Curve2::Circle {
            centre: c1,
            radius: r1,
            inside: true,
        },
        Curve2::Circle {
            centre: c2,
            radius: r2,
            inside: true,
        },
    ];
    let heron = ((-gap + r1 + r2) * (gap + r1 - r2) * (gap - r1 + r2) * (gap + r1 + r2)).sqrt();
    let lens = r1 * r1 * ((gap * gap + r1 * r1 - r2 * r2) / (2.0 * gap * r1)).acos()
        + r2 * r2 * ((gap * gap + r2 * r2 - r1 * r1) / (2.0 * gap * r2)).acos()
        - heron / 2.0;
    let sharp = (PI * (r1 * r1 + r2 * r2) - lens) * length;
    let (t1, t2) = (v.y.atan2(v.x), v.y.atan2(v.x - gap));
    let walls = |model: &Model, solid: &Shape| -> Result<[Shape; 2], Verdict> {
        Ok([
            find_face(
                model,
                solid,
                place(&frame, -r1, 0.0, length / 2.0),
                "wall one",
            )?,
            find_face(
                model,
                solid,
                place(&frame, gap + r2, 0.0, length / 2.0),
                "wall two",
            )?,
        ])
    };
    let mut lo = [0.0, 0.0];
    let mut others = Vec::new();
    let mut input = sharp;
    let mut what = String::new();
    if rng.chance(0.5) {
        let rho = r2 * rng.range(0.03, 0.5);
        let [fa, fb] = walls(&model, &solid)?;
        let edges = shared_edges(&model, &fa, &fb);
        solid = built(
            "round",
            ogeom::fillet::fillet_edges(&mut model, &solid, &edges, rho, T),
        )?;
        if let Some((arc, t)) = old_round(supports, v, 1.0, rho) {
            let (area, _) = region(&[supports[0].piece(v, t[0]), arc, supports[1].piece(t[1], v)]);
            input += 2.0 * area.abs() * length;
            lo = [
                wrap(t[0].y.atan2(t[0].x) - t1) * r1,
                -wrap((t[1].y).atan2(t[1].x - gap) - t2) * r2,
            ];
            others.push(arc);
        }
        what = format!(", creases rounded {rho}");
    }
    let faces = walls(&model, &solid)?;
    let mut carriers = others.clone();
    carriers.push(Piece::Arc {
        centre: c1,
        radius: r1,
        from: 0.0,
        to: TAU,
    });
    carriers.push(Piece::Arc {
        centre: c2,
        radius: r2,
        from: 0.0,
        to: TAU,
    });
    Ok(Scene {
        model,
        solid,
        faces,
        what: format!(
            "drums r {r1:.4} and {r2:.4} {gap:.4} apart, {length:.4} long, frame {}{what}",
            brief(&frame)
        ),
        geometry: Geometry::Section(Section {
            sweep: Sweep::Prism { frame, length },
            supports,
            crease: v,
            reach: [
                round_from(c1, r1, v, 1.0, lo[0], (PI - t1) * r1),
                round_from(c2, r2, v, -1.0, lo[1], t2 * r2),
            ],
            others,
            mirror: Some(|p| P2::new(p.x, -p.y)),
            carriers,
            planes: vec![
                (frame.origin(), frame.z().vector()),
                (place(&frame, 0.0, 0.0, length), frame.z().vector()),
            ],
            angle: None,
            within: None,
        }),
        convex: false,
        sharp,
        input: Some(input),
        radii: (r2 * 0.03, r2 * 1.0),
    })
}

/// The rim edge of a revolved part: the edges `a` and `b` share.
fn rim_edges(model: &Model, a: &Shape, b: &Shape) -> Result<Vec<Shape>, Verdict> {
    let edges = shared_edges(model, a, b);
    if edges.is_empty() {
        return Err(Verdict::of(
            Outcome::Refused,
            "input: the faces share no rim",
        ));
    }
    Ok(edges)
}

/// A section turned about its frame's axis, with the closed form of what a
/// strip takes or adds: by Pappus on the strip's section.
fn ring(pieces: &[Piece]) -> f64 {
    TAU * region(pieces).1.abs()
}

/// A drum's top rim (convex, plane square to a cylinder's axis), bare,
/// chamfered unevenly, or rounded.
fn rim(rng: &mut Rng) -> Result<Scene, Verdict> {
    let frame = random_frame(rng)?;
    let mut model = Model::new();
    let (radius, height) = (rng.range(1.5, 5.0), rng.range(1.5, 5.0));
    let mut solid = built(
        "drum",
        ogeom::algo::make_cylinder(&mut model, frame, radius, height, T),
    )?;
    let v = P2::new(radius, height);
    let inner = P2::new(radius / 2.0, height / 2.0);
    let supports = [
        Curve2::line(P2::new(radius, 0.0), v, inner),
        Curve2::line(P2::new(0.0, height), v, inner),
    ];
    let sharp = PI * radius * radius * height;
    let sweep = Sweep::Revolve { frame };
    let faces_of = |model: &Model, solid: &Shape, lo: [f64; 2]| -> Result<[Shape; 2], Verdict> {
        Ok([
            find_face(
                model,
                solid,
                sweep.place(P2::new(radius, (height - lo[0]) / 2.0), 2.0),
                "wall",
            )?,
            find_face(
                model,
                solid,
                sweep.place(P2::new((radius - lo[1]) / 2.0, height), 2.0),
                "top",
            )?,
        ])
    };
    let mut lo = [0.0, 0.0];
    let mut others = vec![Piece::Seg(P2::new(0.0, 0.0), P2::new(radius, 0.0))];
    let mut input = sharp;
    let mut what = String::new();
    let small = radius.min(height);
    let pick = rng.below(3);
    if pick > 0 {
        let [wall, top] = faces_of(&model, &solid, lo)?;
        let edges = rim_edges(&model, &wall, &top)?;
        if pick == 1 {
            let (dw, dt) = (small * rng.range(0.03, 0.35), small * rng.range(0.03, 0.35));
            solid = built(
                "chamfer",
                ogeom::fillet::chamfer_edge_distances(
                    &mut model, &solid, &edges[0], &wall, dw, dt, T,
                ),
            )?;
            let (c1, c2) = (P2::new(radius, height - dw), P2::new(radius - dt, height));
            input -= ring(&[Piece::Seg(v, c1), Piece::Seg(c1, c2), Piece::Seg(c2, v)]);
            others.push(Piece::Seg(c1, c2));
            lo = [dw, dt];
            what = format!(", chamfered {dw} down the wall by {dt} across the top");
        } else {
            let rho = small * rng.range(0.03, 0.35);
            solid = built(
                "round",
                ogeom::fillet::fillet_edge(&mut model, &solid, &edges[0], rho, T),
            )?;
            let (arc, t) = old_round(supports, v, -1.0, rho).unwrap();
            input -= ring(&[supports[0].piece(v, t[0]), arc, supports[1].piece(t[1], v)]);
            others.push(arc);
            lo = [rho, rho];
            what = format!(", rounded {rho}");
        }
    }
    let faces = faces_of(&model, &solid, lo)?;
    let mut carriers = others.clone();
    carriers.push(Piece::Seg(P2::new(radius, 0.0), v));
    carriers.push(Piece::Seg(P2::new(0.0, height), v));
    Ok(Scene {
        model,
        solid,
        faces,
        what: format!(
            "drum r {radius:.4} h {height:.4}, frame {}{what}",
            brief(&frame)
        ),
        geometry: Geometry::Section(Section {
            sweep,
            supports,
            crease: v,
            reach: [
                along(v, V2::new(0.0, -1.0), lo[0], height),
                along(v, V2::new(-1.0, 0.0), lo[1], radius),
            ],
            others,
            mirror: None,
            carriers,
            planes: Vec::new(),
            angle: None,
            within: None,
        }),
        convex: true,
        sharp,
        input: Some(input),
        radii: (small * 0.02, small * 0.6),
    })
}

/// A block's sides and bottom about a revolved feature at its frame's axis:
/// the block is `w` by `d` centred on the axis, its bottom at height 0.
fn block_planes(frame: &Frame, w: f64, d: f64) -> Vec<(Point, Vector)> {
    vec![
        (place(frame, w / 2.0, 0.0, 0.0), frame.x().vector()),
        (place(frame, -w / 2.0, 0.0, 0.0), frame.x().vector()),
        (place(frame, 0.0, d / 2.0, 0.0), frame.y().vector()),
        (place(frame, 0.0, -d / 2.0, 0.0), frame.y().vector()),
        (frame.origin(), frame.z().vector()),
    ]
}

fn centred_block(
    model: &mut Model,
    frame: &Frame,
    w: f64,
    d: f64,
    h: f64,
) -> Result<Shape, Verdict> {
    let at = frame_in(
        place(frame, -w / 2.0, -d / 2.0, 0.0),
        frame.z().vector(),
        frame.x().vector(),
    )?;
    built("block", ogeom::algo::make_box(model, at, (w, d, h), T))
}

/// A block bored through, the bore's top rim (convex: a cylinder's inside
/// and a plane square to it), bare, chamfered or rounded.
fn bore(rng: &mut Rng) -> Result<Scene, Verdict> {
    let frame = random_frame(rng)?;
    let mut model = Model::new();
    let radius = rng.range(0.8, 2.5);
    let (w, d, height) = (
        2.0 * (radius + rng.range(1.5, 4.0)),
        2.0 * (radius + rng.range(1.5, 4.0)),
        rng.range(1.5, 5.0),
    );
    let half = w.min(d) / 2.0;
    let block = centred_block(&mut model, &frame, w, d, height)?;
    let drill = frame_in(
        place(&frame, 0.0, 0.0, -1.0),
        frame.z().vector(),
        frame.x().vector(),
    )?;
    let drill = built(
        "drill",
        ogeom::algo::make_cylinder(&mut model, drill, radius, height + 2.0, T),
    )?;
    let mut solid = built("bore", ogeom::boolean::cut(&mut model, &block, &drill, T))?;
    let v = P2::new(radius, height);
    let inner = P2::new(radius + 0.5, height / 2.0);
    let supports = [
        Curve2::line(P2::new(radius, 0.0), v, inner),
        Curve2::line(v, P2::new(radius + 1.0, height), inner),
    ];
    let sharp = (w * d - PI * radius * radius) * height;
    let sweep = Sweep::Revolve { frame };
    let faces_of = |model: &Model, solid: &Shape, lo: [f64; 2]| -> Result<[Shape; 2], Verdict> {
        Ok([
            find_face(
                model,
                solid,
                sweep.place(P2::new(radius, (height - lo[0]) / 2.0), 2.0),
                "bore wall",
            )?,
            find_face(
                model,
                solid,
                sweep.place(P2::new(f64::midpoint(radius + lo[1], half), height), 2.0),
                "top",
            )?,
        ])
    };
    let mut lo = [0.0, 0.0];
    let mut others = vec![Piece::Seg(P2::new(radius, 0.0), P2::new(half, 0.0))];
    let mut input = sharp;
    let mut what = String::new();
    let small = height.min(half - radius);
    let pick = rng.below(3);
    if pick > 0 {
        let [wall, top] = faces_of(&model, &solid, lo)?;
        let edges = rim_edges(&model, &wall, &top)?;
        if pick == 1 {
            let c = small * rng.range(0.03, 0.35);
            solid = built(
                "chamfer",
                ogeom::fillet::chamfer_edges(&mut model, &solid, &edges, c, T),
            )?;
            let (c1, c2) = (P2::new(radius, height - c), P2::new(radius + c, height));
            input -= ring(&[Piece::Seg(v, c1), Piece::Seg(c1, c2), Piece::Seg(c2, v)]);
            others.push(Piece::Seg(c1, c2));
            lo = [c, c];
            what = format!(", chamfered {c}");
        } else {
            let rho = small * rng.range(0.03, 0.35);
            solid = built(
                "round",
                ogeom::fillet::fillet_edges(&mut model, &solid, &edges, rho, T),
            )?;
            let (arc, t) = old_round(supports, v, -1.0, rho).unwrap();
            input -= ring(&[supports[0].piece(v, t[0]), arc, supports[1].piece(t[1], v)]);
            others.push(arc);
            lo = [rho, rho];
            what = format!(", rounded {rho}");
        }
    }
    let faces = faces_of(&model, &solid, lo)?;
    let mut carriers = others.clone();
    carriers.push(Piece::Seg(P2::new(radius, 0.0), v));
    carriers.push(Piece::Seg(v, P2::new(half, height)));
    Ok(Scene {
        model,
        solid,
        faces,
        what: format!(
            "block {w:.4} by {d:.4} by {height:.4} bored r {radius:.4}, frame {}{what}",
            brief(&frame)
        ),
        geometry: Geometry::Section(Section {
            sweep,
            supports,
            crease: v,
            reach: [
                along(v, V2::new(0.0, -1.0), lo[0], height),
                along(v, V2::new(1.0, 0.0), lo[1], half - radius),
            ],
            others,
            mirror: None,
            carriers,
            planes: block_planes(&frame, w, d),
            angle: None,
            within: Some(half),
        }),
        convex: true,
        sharp,
        input: Some(input),
        radii: (small * 0.02, small * 0.6),
    })
}

/// A drum fused standing on a block: its foot (concave, a cylinder and a
/// plane square to its axis), bare, rounded, or filled by a cone fused on
/// (plane or B-spline).
fn boss(rng: &mut Rng) -> Result<Scene, Verdict> {
    let frame = random_frame(rng)?;
    let mut model = Model::new();
    let radius = rng.range(1.0, 3.0);
    let (w, d) = (
        2.0 * (radius + rng.range(2.0, 5.0)),
        2.0 * (radius + rng.range(2.0, 5.0)),
    );
    let (base, tall) = (rng.range(1.0, 3.0), rng.range(1.5, 4.0));
    let half = w.min(d) / 2.0;
    let block = centred_block(&mut model, &frame, w, d, base)?;
    let post = frame_in(
        place(&frame, 0.0, 0.0, base - 0.5),
        frame.z().vector(),
        frame.x().vector(),
    )?;
    let post = built(
        "post",
        ogeom::algo::make_cylinder(&mut model, post, radius, tall + 0.5, T),
    )?;
    let mut solid = built("fuse", ogeom::boolean::fuse(&mut model, &block, &post, T))?;
    let v = P2::new(radius, base);
    let supports = [
        Curve2::line(
            v,
            P2::new(radius, base + tall),
            P2::new(radius / 2.0, base + tall / 2.0),
        ),
        Curve2::line(
            v,
            P2::new(radius + 1.0, base),
            P2::new(radius + 1.0, base / 2.0),
        ),
    ];
    let sharp = w * d * base + PI * radius * radius * tall;
    let sweep = Sweep::Revolve { frame };
    let faces_of = |model: &Model, solid: &Shape, lo: [f64; 2]| -> Result<[Shape; 2], Verdict> {
        Ok([
            find_face(
                model,
                solid,
                sweep.place(P2::new(radius, base + f64::midpoint(lo[0], tall)), 2.0),
                "post wall",
            )?,
            find_face(
                model,
                solid,
                sweep.place(P2::new(f64::midpoint(radius + lo[1], half), base), 2.0),
                "top",
            )?,
        ])
    };
    let mut lo = [0.0, 0.0];
    let mut others = vec![Piece::Seg(
        P2::new(0.0, base + tall),
        P2::new(radius, base + tall),
    )];
    let mut input = sharp;
    let mut what = String::new();
    let small = tall.min(half - radius);
    let pick = rng.below(4);
    if pick == 1 {
        let [wall, top] = faces_of(&model, &solid, lo)?;
        let edges = rim_edges(&model, &wall, &top)?;
        let rho = small * rng.range(0.03, 0.35);
        solid = built(
            "round",
            ogeom::fillet::fillet_edges(&mut model, &solid, &edges, rho, T),
        )?;
        let (arc, t) = old_round(supports, v, 1.0, rho).unwrap();
        input += ring(&[supports[0].piece(v, t[0]), arc, supports[1].piece(t[1], v)]);
        others.push(arc);
        lo = [rho, rho];
        what = format!(", foot rounded {rho}");
    } else if pick > 1 {
        // A cone through (radius + dt, base) and (radius, base + dw), run on
        // past both into the block and the post.
        let (dw, dt) = (small * rng.range(0.03, 0.35), small * rng.range(0.03, 0.35));
        let slope = dt / dw;
        // Run on into the block no further than keeps the cone's foot
        // inside the block's sides.
        let below = (0.3 * base).min((half - radius - dt) * 0.5 / slope);
        let above = 0.3 * dw;
        let foot = frame_in(
            place(&frame, 0.0, 0.0, base - below),
            frame.z().vector(),
            frame.x().vector(),
        )?;
        let mut cone = built(
            "cone",
            ogeom::algo::make_cone(
                &mut model,
                foot,
                radius + dt + slope * below,
                radius - slope * above,
                below + dw + above,
                T,
            ),
        )?;
        if pick == 3 {
            cone = built("to_nurbs", ogeom::algo::to_nurbs(&mut model, &cone, T))?;
        }
        solid = built("fuse", ogeom::boolean::fuse(&mut model, &solid, &cone, T))?;
        let (c1, c2) = (P2::new(radius, base + dw), P2::new(radius + dt, base));
        input += ring(&[Piece::Seg(v, c1), Piece::Seg(c1, c2), Piece::Seg(c2, v)]);
        others.push(Piece::Seg(c1, c2));
        lo = [dw, dt];
        what = format!(
            ", foot filled by a {} cone {dw} up the post by {dt} across the top",
            if pick == 3 { "B-spline" } else { "plane" }
        );
    }
    let faces = faces_of(&model, &solid, lo)?;
    let mut carriers = others.clone();
    carriers.push(Piece::Seg(v, P2::new(radius, base + tall)));
    carriers.push(Piece::Seg(v, P2::new(half, base)));
    carriers.push(Piece::Seg(
        P2::new(0.0, base + tall),
        P2::new(radius, base + tall),
    ));
    Ok(Scene {
        model,
        solid,
        faces,
        what: format!(
            "post r {radius:.4} h {tall:.4} on a block {w:.4} by {d:.4} by {base:.4}, frame {}{what}",
            brief(&frame)
        ),
        geometry: Geometry::Section(Section {
            sweep,
            supports,
            crease: v,
            reach: [
                along(v, V2::new(0.0, 1.0), lo[0], tall),
                along(v, V2::new(1.0, 0.0), lo[1], half - radius),
            ],
            others,
            mirror: None,
            carriers,
            planes: block_planes(&frame, w, d),
            angle: None,
            within: Some(half),
        }),
        convex: false,
        sharp,
        input: Some(input),
        radii: (small * 0.02, small * 0.6),
    })
}

/// A cone frustum's top rim (convex, a plane and a cone about one axis),
/// bare or rounded.
fn frustum(rng: &mut Rng) -> Result<Scene, Verdict> {
    let frame = random_frame(rng)?;
    let mut model = Model::new();
    let r0 = rng.range(1.5, 5.0);
    let r1 = r0
        * if rng.chance(0.5) {
            rng.range(0.3, 0.85)
        } else {
            rng.range(1.15, 1.6)
        };
    let height = rng.range(1.5, 5.0);
    let mut solid = built(
        "cone",
        ogeom::algo::make_cone(&mut model, frame, r0, r1, height, T),
    )?;
    let v = P2::new(r1, height);
    let inner = P2::new(r0.min(r1) / 2.0, height / 2.0);
    let low = P2::new(r0, 0.0);
    let supports = [
        Curve2::line(low, v, inner),
        Curve2::line(P2::new(0.0, height), v, inner),
    ];
    let slant = low.distance(v);
    let down = (low - v) / slant;
    let sharp = PI * height * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0;
    let sweep = Sweep::Revolve { frame };
    let faces_of = |model: &Model, solid: &Shape, lo: [f64; 2]| -> Result<[Shape; 2], Verdict> {
        let wall = v + down * f64::midpoint(lo[0], slant);
        Ok([
            find_face(model, solid, sweep.place(wall, 2.0), "wall")?,
            find_face(
                model,
                solid,
                sweep.place(P2::new((r1 - lo[1]) / 2.0, height), 2.0),
                "top",
            )?,
        ])
    };
    let mut lo = [0.0, 0.0];
    let mut others = vec![Piece::Seg(P2::new(0.0, 0.0), low)];
    let mut input = sharp;
    let mut what = String::new();
    let small = r1.min(height);
    if rng.chance(0.5) {
        let [wall, top] = faces_of(&model, &solid, lo)?;
        let edges = rim_edges(&model, &wall, &top)?;
        let rho = small * rng.range(0.03, 0.3);
        solid = built(
            "round",
            ogeom::fillet::fillet_edges(&mut model, &solid, &edges, rho, T),
        )?;
        let (arc, t) = old_round(supports, v, -1.0, rho).unwrap();
        input -= ring(&[supports[0].piece(v, t[0]), arc, supports[1].piece(t[1], v)]);
        others.push(arc);
        lo = [(t[0] - v).dot(down), r1 - t[1].x];
        what = format!(", rounded {rho}");
    }
    let faces = faces_of(&model, &solid, lo)?;
    let mut carriers = others.clone();
    carriers.push(Piece::Seg(low, v));
    carriers.push(Piece::Seg(P2::new(0.0, height), v));
    Ok(Scene {
        model,
        solid,
        faces,
        what: format!(
            "frustum r {r0:.4} to {r1:.4} h {height:.4}, frame {}{what}",
            brief(&frame)
        ),
        geometry: Geometry::Section(Section {
            sweep,
            supports,
            crease: v,
            reach: [
                along(v, down, lo[0], slant),
                along(v, V2::new(-1.0, 0.0), lo[1], r1),
            ],
            others,
            mirror: None,
            carriers,
            planes: Vec::new(),
            angle: None,
            within: None,
        }),
        convex: true,
        sharp,
        input: Some(input),
        radii: (small * 0.02, small * 0.6),
    })
}

/// A cone frustum fused standing on a block: its foot (concave, a cone and
/// a plane square to its axis), bare or rounded.
fn cone_foot(rng: &mut Rng) -> Result<Scene, Verdict> {
    let frame = random_frame(rng)?;
    let mut model = Model::new();
    let r0 = rng.range(1.0, 3.0);
    let r1 = r0 * rng.range(0.3, 1.3);
    let tall = rng.range(1.5, 4.0);
    let (w, d, base) = (
        2.0 * (r0.max(r1) + rng.range(2.0, 5.0)),
        2.0 * (r0.max(r1) + rng.range(2.0, 5.0)),
        rng.range(1.0, 3.0),
    );
    let half = w.min(d) / 2.0;
    let block = centred_block(&mut model, &frame, w, d, base)?;
    let below = 0.3 * base;
    let foot = frame_in(
        place(&frame, 0.0, 0.0, base - below),
        frame.z().vector(),
        frame.x().vector(),
    )?;
    let flare = (r0 - r1) / tall;
    let cone = built(
        "cone",
        ogeom::algo::make_cone(&mut model, foot, r0 + flare * below, r1, tall + below, T),
    )?;
    let mut solid = built("fuse", ogeom::boolean::fuse(&mut model, &block, &cone, T))?;
    let v = P2::new(r0, base);
    let top_end = P2::new(r1, base + tall);
    let supports = [
        Curve2::line(v, top_end, P2::new(r0.min(r1) / 2.0, base + tall / 2.0)),
        Curve2::line(v, P2::new(r0 + 1.0, base), P2::new(r0 + 1.0, base / 2.0)),
    ];
    let slant = v.distance(top_end);
    let up = (top_end - v) / slant;
    let sharp = w * d * base + PI * tall * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0;
    let sweep = Sweep::Revolve { frame };
    let faces_of = |model: &Model, solid: &Shape, lo: [f64; 2]| -> Result<[Shape; 2], Verdict> {
        Ok([
            find_face(
                model,
                solid,
                sweep.place(v + up * f64::midpoint(lo[0], slant), 2.0),
                "cone wall",
            )?,
            find_face(
                model,
                solid,
                sweep.place(P2::new(f64::midpoint(r0 + lo[1], half), base), 2.0),
                "top",
            )?,
        ])
    };
    let mut lo = [0.0, 0.0];
    let mut others = vec![Piece::Seg(P2::new(0.0, base + tall), top_end)];
    let mut input = sharp;
    let mut what = String::new();
    let small = slant.min(half - r0);
    if rng.chance(0.5) {
        let [wall, top] = faces_of(&model, &solid, lo)?;
        let edges = rim_edges(&model, &wall, &top)?;
        let rho = small * rng.range(0.03, 0.3);
        solid = built(
            "round",
            ogeom::fillet::fillet_edges(&mut model, &solid, &edges, rho, T),
        )?;
        let (arc, t) = old_round(supports, v, 1.0, rho).unwrap();
        input += ring(&[supports[0].piece(v, t[0]), arc, supports[1].piece(t[1], v)]);
        others.push(arc);
        lo = [(t[0] - v).dot(up), t[1].x - r0];
        what = format!(", foot rounded {rho}");
    }
    let faces = faces_of(&model, &solid, lo)?;
    let mut carriers = others.clone();
    carriers.push(Piece::Seg(v, top_end));
    carriers.push(Piece::Seg(v, P2::new(half, base)));
    Ok(Scene {
        model,
        solid,
        faces,
        what: format!(
            "cone r {r0:.4} to {r1:.4} h {tall:.4} on a block {w:.4} by {d:.4} by {base:.4}, frame {}{what}",
            brief(&frame)
        ),
        geometry: Geometry::Section(Section {
            sweep,
            supports,
            crease: v,
            reach: [
                along(v, up, lo[0], slant),
                along(v, V2::new(1.0, 0.0), lo[1], half - r0),
            ],
            others,
            mirror: None,
            carriers,
            planes: block_planes(&frame, w, d),
            angle: None,
            within: Some(half),
        }),
        convex: false,
        sharp,
        input: Some(input),
        radii: (small * 0.02, small * 0.6),
    })
}

/// `face` of `shape` split along where the plane through `on` square to
/// `normal` crosses it.
fn split_by_plane(
    model: &mut Model,
    shape: &Shape,
    face: &Shape,
    on: Point,
    normal: Vector,
) -> Result<Shape, Verdict> {
    use ogeom::geom::Curve3d as _;
    let data = model
        .node(face)
        .and_then(|n| n.data().as_face())
        .ok_or_else(|| Verdict::of(Outcome::Refused, "input: split: not a face"))?;
    let surface = model
        .geometry()
        .surface(data.surface)
        .ok_or_else(|| Verdict::of(Outcome::Refused, "input: split: no surface"))?
        .clone();
    let normal = Direction::new(normal, T).map_err(|e| input_error("split", e))?;
    let plane: ogeom::geom::SurfaceGeometry = ogeom::geom::PlaneSurface::over(
        ogeom::math::Plane::through(on, normal),
        (-200.0, 200.0),
        (-200.0, 200.0),
    )
    .map_err(|e| input_error("split", e))?
    .into();
    let Ok(ogeom::intersect::SurfaceIntersection::Along(curves)) =
        ogeom::intersect::intersect_surfaces(
            &plane,
            &surface,
            ogeom::intersect::IntersectOptions::default(),
            T,
        )
    else {
        return Err(Verdict::of(
            Outcome::Refused,
            "input: split: the plane misses the face",
        ));
    };
    let mut edges = Vec::new();
    for c in curves {
        let range = c.curve.domain();
        edges.push(built(
            "split edge",
            ogeom::algo::make_edge(model, c.curve, range, T),
        )?);
    }
    built(
        "split",
        ogeom::heal::split_face(
            model,
            shape,
            face,
            &edges,
            ogeom::heal::Projection::OnFace,
            T,
        ),
    )
}

/// A drum leaning on a block, its foot a closed ellipse on the top: bare
/// (the wall and the top share the foot's edges), its wall split a little
/// above the foot, its foot rounded small, or a wedge fused into the foot
/// over part of its run.
fn leaning(rng: &mut Rng) -> Result<Scene, Verdict> {
    let g = random_frame(rng)?;
    let mut model = Model::new();
    let radius = rng.range(2.0, 4.0);
    let lean = rng.range(3.0, 25.0).to_radians();
    let toward = rng.range(0.0, TAU);
    let side = 6.0 * radius + 6.0;
    // The drum's base sunk in the block far enough that its tilted rim
    // stays under the top, and the block deep enough to hold it.
    let tip = radius * lean.sin();
    let depth = tip + 0.3;
    let base = depth + tip + rng.range(0.5, 2.0);
    let rise = rng.range(3.0, 6.0);
    let foot_at = (
        side / 2.0 + rng.range(-1.0, 1.0),
        side / 2.0 + rng.range(-1.0, 1.0),
    );
    let foot = place(&g, foot_at.0, foot_at.1, base);
    let axis = dir(
        &g,
        lean.sin() * toward.cos(),
        lean.sin() * toward.sin(),
        lean.cos(),
    );
    let up = g.z().vector();
    let e = across(axis, rng);
    let e = (e - axis * e.dot(axis)) / (e - axis * e.dot(axis)).magnitude();
    let f = axis.cross(e);
    let block = frame_in(g.origin(), up, g.x().vector())?;
    let block = built(
        "block",
        ogeom::algo::make_box(&mut model, block, (side, side, base), T),
    )?;
    let sunk = depth / axis.dot(up);
    let drum = frame_in(foot - axis * sunk, axis, e)?;
    let drum = built(
        "drum",
        ogeom::algo::make_cylinder(&mut model, drum, radius, sunk + rise, T),
    )?;
    let sharp_solid = built("fuse", ogeom::boolean::fuse(&mut model, &block, &drum, T))?;
    let sharp = measure(&model, &sharp_solid)
        .map_err(|e| input_error("volume", e))?
        .volume;
    let wall_point = foot + axis * (rise * 0.8) + (e * 2.0_f64.cos() + f * 2.0_f64.sin()) * radius;
    let top_point = place(&g, side * 0.05, side * 0.05, base);
    let mut solid = sharp_solid;
    let mut split = None;
    let mut strip = None;
    let mut cover: Option<Reach3> = None;
    let mut others = vec![
        Patch::Quad {
            at: g.origin(),
            u: dir(&g, side, 0.0, 0.0),
            v: dir(&g, 0.0, 0.0, base),
        },
        Patch::Quad {
            at: g.origin(),
            u: dir(&g, 0.0, side, 0.0),
            v: dir(&g, 0.0, 0.0, base),
        },
        Patch::Quad {
            at: place(&g, side, 0.0, 0.0),
            u: dir(&g, 0.0, side, 0.0),
            v: dir(&g, 0.0, 0.0, base),
        },
        Patch::Quad {
            at: place(&g, 0.0, side, 0.0),
            u: dir(&g, side, 0.0, 0.0),
            v: dir(&g, 0.0, 0.0, base),
        },
        Patch::Disk {
            centre: foot + axis * rise,
            e1: e,
            e2: f,
            radius,
        },
    ];
    let mut what = String::new();
    match rng.below(4) {
        1 => {
            let s = rng.range(0.05, 0.6);
            let wall = find_face(&model, &solid, wall_point, "wall")?;
            solid = split_by_plane(&mut model, &solid, &wall, foot + up * s, up)?;
            split = Some(s);
            what = format!(", wall split {s} above the top");
        }
        2 => {
            let rho = radius * rng.range(0.03, 0.3);
            let wall = find_face(&model, &solid, wall_point, "wall")?;
            let top = find_face(&model, &solid, top_point, "top")?;
            let edges = shared_edges(&model, &wall, &top);
            solid = built(
                "round",
                ogeom::fillet::fillet_edges(&mut model, &solid, &edges, rho, T),
            )?;
            strip = Some(rho);
            what = format!(", foot rounded {rho}");
        }
        3 => {
            // A wedge in the vertical plane through the axis's foot at
            // azimuth `at`, its slope from `leg_t` out across the top to
            // `leg_w` up the wall, `width` wide.
            let at = rng.range(0.0, TAU);
            let out = dir(&g, at.cos(), at.sin(), 0.0);
            let across_ = up.cross(out);
            // The crease point at this azimuth: where the ray from the foot
            // along `out` meets the wall.
            let p0 = foot + out * (radius / (out - axis * out.dot(axis)).magnitude());
            let (leg_t, leg_w) = (radius * rng.range(0.05, 0.4), radius * rng.range(0.05, 0.4));
            let width = radius * rng.range(0.2, 1.0);
            let delta = 0.2;
            let back = delta + leg_w * lean.tan() + 0.1;
            let q = |a: f64, z: f64| p0 + out * a + up * z;
            let (q0, q1, q2) = (
                q(-back, -delta),
                q(leg_t * (1.0 + delta / leg_w), -delta),
                q(-back, leg_w * (1.0 + back / leg_t)),
            );
            let section = [q0, q1, q2].map(|p| p - across_ * (width / 2.0));
            let wire = built(
                "wedge wire",
                ogeom::algo::make_polygon(&mut model, &section, true, T),
            )?;
            let normal = Direction::new(across_, T).map_err(|e| input_error("wedge", e))?;
            let plane =
                ogeom::geom::PlaneSurface::new(ogeom::math::Plane::through(section[0], normal));
            let face = built(
                "wedge face",
                ogeom::algo::make_face(&mut model, plane.into(), &[wire], T),
            )?;
            let wedge = built(
                "wedge",
                ogeom::algo::make_prism(&mut model, &face, across_ * width, T),
            )?;
            solid = built("fuse", ogeom::boolean::fuse(&mut model, &solid, &wedge, T))?;
            // A contact inside the wedge is off its face: how far inside
            // the wedge's section and between its ends a point stands.
            let tri = [
                (-back, -delta),
                (leg_t * (1.0 + delta / leg_w), -delta),
                (-back, leg_w * (1.0 + back / leg_t)),
            ];
            cover = Some(Box::new(move |p: Point| {
                let d = p - p0;
                let (a, z, c) = (d.dot(out), d.dot(up), d.dot(across_));
                let centre = (
                    (tri[0].0 + tri[1].0 + tri[2].0) / 3.0,
                    (tri[0].1 + tri[1].1 + tri[2].1) / 3.0,
                );
                let mut inside = width / 2.0 - c.abs();
                for k in 0..3 {
                    let (p, q) = (tri[k], tri[(k + 1) % 3]);
                    let (ex, ez) = (q.0 - p.0, q.1 - p.1);
                    let len = ex.hypot(ez);
                    let mut n = (-ez / len, ex / len);
                    if n.0 * (centre.0 - p.0) + n.1 * (centre.1 - p.1) < 0.0 {
                        n = (-n.0, -n.1);
                    }
                    inside = inside.min(n.0 * (a - p.0) + n.1 * (z - p.1));
                }
                inside
            }));
            let shift = across_ * width;
            others.push(Patch::Quad {
                at: section[1],
                u: section[2] - section[1],
                v: shift,
            });
            others.push(Patch::Tri {
                a: section[0],
                b: section[1],
                c: section[2],
            });
            others.push(Patch::Tri {
                a: section[0] + shift,
                b: section[1] + shift,
                c: section[2] + shift,
            });
            others.push(Patch::Quad {
                at: section[0],
                u: section[1] - section[0],
                v: shift,
            });
            others.push(Patch::Quad {
                at: section[0],
                u: section[2] - section[0],
                v: shift,
            });
            what = format!(
                ", a wedge {leg_t} across the top by {leg_w} up the wall, {width} wide at azimuth {at}"
            );
        }
        _ => {}
    }
    let faces = [
        find_face(&model, &solid, wall_point, "wall")?,
        find_face(&model, &solid, top_point, "top")?,
    ];
    let supports = [
        Surf::Cylinder {
            at: foot,
            axis,
            radius,
            inside: true,
        },
        Surf::Plane { at: foot, out: up },
    ];
    let tan = lean.tan();
    let spine: SpineOf = Box::new(move |sides, r| {
        let big = radius + sides[0] * r;
        let lift = sides[1] * r;
        let (ez, fz, az) = (e.dot(up), f.dot(up), axis.dot(up));
        Spine {
            at: Box::new(move |t: f64| {
                let w = (lift - big * (t.cos() * ez + t.sin() * fz)) / az;
                foot + axis * w + (e * t.cos() + f * t.sin()) * big
            }),
            speed: big * (1.0 + tan),
        }
    });
    let (g2, side2) = (g, side);
    let reach: [Reach3; 2] = [
        Box::new(move |q: Point| {
            let high = rise - (q - foot).dot(axis);
            match split {
                Some(s) => high.min((q - foot).dot(up) - s),
                None => high,
            }
        }),
        Box::new(move |q: Point| {
            let d = q - g2.origin();
            let (x, y) = (d.dot(g2.x().vector()), d.dot(g2.y().vector()));
            x.min(y).min(side2 - x).min(side2 - y)
        }),
    ];
    let planes = vec![(g.origin(), up), (foot + axis * rise, axis)];
    Ok(Scene {
        model,
        solid,
        faces,
        what: format!(
            "drum r {radius:.4} leaning {:.3} deg toward {toward:.4} on a block {side:.4} by {base:.4}, foot at {foot_at:?}, frame {}{what}",
            lean.to_degrees(),
            brief(&g)
        ),
        geometry: Geometry::Marched(Marched {
            supports,
            spine,
            reach,
            others,
            planes,
            strip,
            cover,
        }),
        convex: false,
        sharp,
        input: None,
        radii: (radius * 0.04, radius * 0.6),
    })
}

/// A tee: a pipe fused into a bigger one, its axis crossing the bigger's
/// at a slant and off to one side: the walls meet along a closed concave
/// crease. Bare, or the crease rounded small.
fn tee(rng: &mut Rng) -> Result<Scene, Verdict> {
    let g = random_frame(rng)?;
    let mut model = Model::new();
    let r1 = rng.range(2.0, 4.0);
    let r2 = r1 * rng.range(0.3, 0.6);
    let slant = rng.range(0.0, 25.0).to_radians();
    let off = (r1 - r2) * rng.range(-0.5, 0.5);
    let reach_x = r1 + r2 + 4.0;
    let tall = r1 + 2.0 * r2 + rng.range(2.0, 4.0);
    let big_frame = frame_in(
        place(&g, -reach_x, 0.0, 0.0),
        g.x().vector(),
        -g.z().vector(),
    )?;
    let big = built(
        "pipe",
        ogeom::algo::make_cylinder(&mut model, big_frame, r1, 2.0 * reach_x, T),
    )?;
    let a2 = dir(&g, slant.sin(), 0.0, slant.cos());
    let e2 = dir(&g, slant.cos(), 0.0, -slant.sin());
    let f2 = g.y().vector();
    let start = place(&g, 0.0, off, 0.0);
    let small_frame = frame_in(start, a2, -e2)?;
    let small = built(
        "branch",
        ogeom::algo::make_cylinder(&mut model, small_frame, r2, tall, T),
    )?;
    let sharp_solid = built("fuse", ogeom::boolean::fuse(&mut model, &big, &small, T))?;
    let sharp = measure(&model, &sharp_solid)
        .map_err(|e| input_error("volume", e))?
        .volume;
    let big_point = place(&g, reach_x * 0.8, r1, 0.0);
    let small_point = start + a2 * (tall * 0.85) + f2 * r2;
    let mut solid = sharp_solid;
    let mut strip = None;
    let mut what = String::new();
    if rng.chance(0.5) {
        let rho = r2 * rng.range(0.03, 0.3);
        let a = find_face(&model, &solid, big_point, "pipe wall")?;
        let b = find_face(&model, &solid, small_point, "branch wall")?;
        let edges = shared_edges(&model, &a, &b);
        solid = built(
            "round",
            ogeom::fillet::fillet_edges(&mut model, &solid, &edges, rho, T),
        )?;
        strip = Some(rho);
        what = format!(", crease rounded {rho}");
    }
    let faces = [
        find_face(&model, &solid, big_point, "pipe wall")?,
        find_face(&model, &solid, small_point, "branch wall")?,
    ];
    let supports = [
        Surf::Cylinder {
            at: g.origin(),
            axis: g.x().vector(),
            radius: r1,
            inside: true,
        },
        Surf::Cylinder {
            at: start,
            axis: a2,
            radius: r2,
            inside: true,
        },
    ];
    let (gx, gy, gz, o) = (g.x().vector(), g.y().vector(), g.z().vector(), g.origin());
    let spine: SpineOf = Box::new(move |sides, r| {
        let (b1, b2) = (r1 + sides[0] * r, r2 + sides[1] * r);
        let (sb, cb) = slant.sin_cos();
        let zmin = (b1 * b1 - (off.abs() + b2).powi(2)).max(1e-12).sqrt();
        let dz = (off.abs() + b2) * b2 / zmin;
        let dw = (dz + b2 * sb) / cb;
        Spine {
            at: Box::new(move |t: f64| {
                let y = off + b2 * t.sin();
                let z = (b1 * b1 - y * y).max(0.0).sqrt();
                let w = (z + b2 * t.cos() * sb) / cb;
                let x = w * sb + b2 * t.cos() * cb;
                o + gx * x + gy * y + gz * z
            }),
            speed: dw + b2,
        }
    });
    let reach: [Reach3; 2] = [
        Box::new(move |q: Point| reach_x - (q - o).dot(gx).abs()),
        Box::new(move |q: Point| {
            let w = (q - start).dot(a2);
            w.min(tall - w)
        }),
    ];
    let others = vec![
        Patch::Disk {
            centre: start + a2 * tall,
            e1: e2,
            e2: f2,
            radius: r2,
        },
        Patch::Disk {
            centre: place(&g, reach_x, 0.0, 0.0),
            e1: gy,
            e2: gz,
            radius: r1,
        },
        Patch::Disk {
            centre: place(&g, -reach_x, 0.0, 0.0),
            e1: gy,
            e2: gz,
            radius: r1,
        },
    ];
    Ok(Scene {
        model,
        solid,
        faces,
        what: format!(
            "pipe r {r1:.4} with a branch r {r2:.4} at {:.3} deg, {off:.4} off, frame {}{what}",
            slant.to_degrees(),
            brief(&g)
        ),
        geometry: Geometry::Marched(Marched {
            supports,
            spine,
            reach,
            others,
            planes: Vec::new(),
            strip,
            cover: None,
        }),
        convex: false,
        sharp,
        input: None,
        radii: (r2 * 0.04, r2 * 0.6),
    })
}

/// A drum cut off by a slanted plane: its top rim (convex, an ellipse
/// where a plane crosses a cylinder at a slant), bare, its wall split a
/// little below the rim, or the rim rounded small.
fn slanted(rng: &mut Rng) -> Result<Scene, Verdict> {
    let g = random_frame(rng)?;
    let mut model = Model::new();
    let radius = rng.range(1.5, 4.0);
    let tilt = rng.range(5.0, 35.0).to_radians();
    let toward = rng.range(0.0, TAU);
    let cap_at = rng.range(3.0, 6.0) + radius * tilt.tan();
    let height = cap_at + radius * tilt.tan() + 2.0;
    let n = dir(
        &g,
        tilt.sin() * toward.cos(),
        tilt.sin() * toward.sin(),
        tilt.cos(),
    );
    let p = place(&g, 0.0, 0.0, cap_at);
    let drum = built(
        "drum",
        ogeom::algo::make_cylinder(&mut model, g, radius, height, T),
    )?;
    let size = 4.0 * radius + height + 4.0;
    let ux = across(n, rng);
    let ux = (ux - n * ux.dot(n)) / (ux - n * ux.dot(n)).magnitude();
    let uy = n.cross(ux);
    let lid = frame_in(p - ux * (size / 2.0) - uy * (size / 2.0), n, ux)?;
    let lid = built(
        "box",
        ogeom::algo::make_box(&mut model, lid, (size, size, size), T),
    )?;
    let sharp_solid = built("cut", ogeom::boolean::cut(&mut model, &drum, &lid, T))?;
    let sharp = measure(&model, &sharp_solid)
        .map_err(|e| input_error("volume", e))?
        .volume;
    let wall_point = place(&g, radius * 2.0_f64.cos(), radius * 2.0_f64.sin(), 0.5);
    let mut solid = sharp_solid;
    let mut split = None;
    let mut strip = None;
    let mut what = String::new();
    match rng.below(3) {
        1 => {
            let s = rng.range(0.05, 0.6);
            let wall = find_face(&model, &solid, wall_point, "wall")?;
            solid = split_by_plane(&mut model, &solid, &wall, p - n * s, n)?;
            split = Some(s);
            what = format!(", wall split {s} below the rim");
        }
        2 => {
            let rho = radius * rng.range(0.03, 0.3);
            let wall = find_face(&model, &solid, wall_point, "wall")?;
            let cap = find_face(&model, &solid, p, "cap")?;
            let edges = shared_edges(&model, &wall, &cap);
            solid = built(
                "round",
                ogeom::fillet::fillet_edges(&mut model, &solid, &edges, rho, T),
            )?;
            strip = Some(rho);
            what = format!(", rim rounded {rho}");
        }
        _ => {}
    }
    let faces = [
        find_face(&model, &solid, wall_point, "wall")?,
        find_face(&model, &solid, p, "cap")?,
    ];
    let (gz, o) = (g.z().vector(), g.origin());
    let (gx, gy) = (g.x().vector(), g.y().vector());
    let supports = [
        Surf::Cylinder {
            at: o,
            axis: gz,
            radius,
            inside: true,
        },
        Surf::Plane { at: p, out: n },
    ];
    let tan = tilt.tan();
    let spine: SpineOf = Box::new(move |sides, r| {
        let big = radius + sides[0] * r;
        let lift = sides[1] * r;
        let (nx, ny, nz) = (n.dot(gx), n.dot(gy), n.dot(gz));
        Spine {
            at: Box::new(move |t: f64| {
                let z = cap_at + (lift - big * (nx * t.cos() + ny * t.sin())) / nz;
                o + gx * (big * t.cos()) + gy * (big * t.sin()) + gz * z
            }),
            speed: big * (1.0 + tan),
        }
    });
    let reach: [Reach3; 2] = [
        Box::new(move |q: Point| {
            let below = -(q - p).dot(n);
            let low = (q - o).dot(gz);
            match split {
                Some(s) => (below - s).min(low),
                None => below.min(low),
            }
        }),
        Box::new(move |q: Point| {
            let w = q - o;
            radius - (w - gz * w.dot(gz)).magnitude()
        }),
    ];
    let others = vec![Patch::Disk {
        centre: o,
        e1: gx,
        e2: gy,
        radius,
    }];
    Ok(Scene {
        model,
        solid,
        faces,
        what: format!(
            "drum r {radius:.4} cut at {cap_at:.4} tilted {:.3} deg toward {toward:.4}, frame {}{what}",
            tilt.to_degrees(),
            brief(&g)
        ),
        geometry: Geometry::Marched(Marched {
            supports,
            spine,
            reach,
            others,
            planes: Vec::new(),
            strip,
            cover: None,
        }),
        convex: true,
        sharp,
        input: None,
        radii: (radius * 0.03, radius * 0.5),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The section integral of the corner between two planes is the closed
    /// form, and Pappus on a square's corner is the known ring.
    #[test]
    fn the_section_integrals_match_their_closed_forms() {
        for angle in [0.5_f64, 1.2, 2.5] {
            let v = P2::new(0.0, 0.0);
            let (a, b) = (
                P2::new(3.0, 0.0),
                P2::new(3.0 * angle.cos(), 3.0 * angle.sin()),
            );
            let m = P2::new((a.x + b.x) / 3.0, (a.y + b.y) / 3.0);
            let supports = [Curve2::line(v, a, m), Curve2::line(v, b, m)];
            let r = 0.4;
            let c = meet(supports[0].offset(-r), supports[1].offset(-r), v).unwrap();
            let (ta, tb) = (supports[0].foot(c), supports[1].foot(c));
            let (area, _) = region(&[
                Piece::Seg(v, ta),
                Piece::arc(c, r, ta, tb),
                Piece::Seg(tb, v),
            ]);
            let want = r * r * (1.0 / (angle / 2.0).tan() - (PI - angle) / 2.0);
            assert!((area.abs() - want).abs() < 1e-13, "{area} against {want}");
        }
        // A drum's rim of radius 10 rounded 1: the ring 2 pi (R - k r) (1 - pi/4) r^2.
        let (radius, r) = (10.0, 1.0);
        let v = P2::new(radius, 5.0);
        let m = P2::new(5.0, 2.0);
        let s = [
            Curve2::line(P2::new(radius, 0.0), v, m),
            Curve2::line(P2::new(0.0, 5.0), v, m),
        ];
        let c = meet(s[0].offset(-r), s[1].offset(-r), v).unwrap();
        let (ta, tb) = (s[0].foot(c), s[1].foot(c));
        let got = ring(&[
            Piece::Seg(v, ta),
            Piece::arc(c, r, ta, tb),
            Piece::Seg(tb, v),
        ]);
        let k = (10.0 - 3.0 * PI) / (3.0 * (4.0 - PI));
        let want = TAU * (radius - k * r) * (1.0 - PI / 4.0) * r * r;
        assert!((got - want).abs() < 1e-12, "{got} against {want}");
    }

    /// The marched integral over an upright drum's foot is the torus ring
    /// by Pappus: the slice integral needs no shared axis, so the upright
    /// case is a check on it.
    #[test]
    fn the_slice_integral_is_pappus_on_an_upright_foot() {
        let (radius, r) = (6.0, 2.0);
        let foot = Point::new(0.0, 0.0, 10.0);
        let axis = Vector::new(0.0, 0.0, 1.0);
        let (e, f) = (Vector::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0));
        let m = Marched {
            supports: [
                Surf::Cylinder {
                    at: foot,
                    axis,
                    radius,
                    inside: true,
                },
                Surf::Plane {
                    at: foot,
                    out: axis,
                },
            ],
            spine: Box::new(move |sides, r| {
                let big = radius + sides[0] * r;
                let lift = sides[1] * r;
                Spine {
                    at: Box::new(move |t: f64| {
                        foot + axis * lift + (e * t.cos() + f * t.sin()) * big
                    }),
                    speed: big,
                }
            }),
            reach: [Box::new(|_| 1.0), Box::new(|_| 1.0)],
            others: Vec::new(),
            planes: Vec::new(),
            strip: None,
            cover: None,
        };
        let spine = (m.spine)([1.0, 1.0], r);
        let got = m.corner_volume(&spine, [1.0, 1.0], r).unwrap();
        let k = (10.0 - 3.0 * PI) / (3.0 * (4.0 - PI));
        let want = TAU * (radius + k * r) * (1.0 - PI / 4.0) * r * r;
        assert!((got - want).abs() < want * 1e-9, "{got} against {want}");
        // And convex, inside the drum's rim.
        let spine = (m.spine)([-1.0, -1.0], r);
        let got = m.corner_volume(&spine, [-1.0, -1.0], r).unwrap();
        let want = TAU * (radius - k * r) * (1.0 - PI / 4.0) * r * r;
        assert!((got - want).abs() < want * 1e-9, "{got} against {want}");
    }
}
