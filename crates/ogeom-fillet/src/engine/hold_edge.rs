//! Hold edge (ADR 0001, amendment "overflow handling"): a round whose
//! contact on one of its faces would run past a sharp edge of that face
//! keeps the edge. Where the contact lies past it, the ball touches the
//! other face and rests on the edge: its centre stays on the other face's
//! offset at the radius from the edge's curve, and the round there is the
//! canal surface those balls sweep, through the held edge.
//!
//! Built for a straight edge between two planes, the held edge a line or a
//! closed circle lying in the plane it bounds, the round ending flush on
//! end faces square to its edge. Each ball's arc runs from its contact on
//! the other face to its foot on the held edge, in the plane through its
//! centre square to the spine (both lie on its characteristic circle), so
//! the section where the contact crosses the held edge is the band's own
//! and the canal meets the band G1 there. The canal is fitted through its
//! sections by the planes square to the edge, so its ends lie in an end
//! face's plane or on a band's section. Past an end face the centres run
//! on as they are and the end face cuts the round (amendment "M5f").
//!
//! The region between the edge and the round is one closed tool, cut from
//! the solid; it is built only where nothing else of the solid stands
//! within its reach, the cut takes the tool's whole volume, the result is
//! valid and the held edge is still there over the run.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve3d as _, CylinderSurface, SurfaceGeometry};
use ogeom_math::{Aabb, Cylinder, Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape};

use super::curved::{Curved, SheetCurve, read};
use super::overflow::{ToolFace, arc, boundary_from, closed, segment};
use super::trimmed::{faces_of, outward};
use super::{Id, id};

/// Planar sections through the canal.
const STATIONS: usize = 160;
/// Points across each section.
const ACROSS: usize = 65;

/// The sharp edge held: a line through `at` along `along`, or the circle
/// of `radius` about `centre` in the plane square to `normal`.
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
            } => {
                let q = p - centre;
                let radial = q - normal * q.dot(normal);
                let m = radial.magnitude();
                (m > 0.0).then(|| centre + radial * (radius / m))
            }
        }
    }

    fn distance(self, p: Point) -> Option<f64> {
        self.foot(p).map(|f| f.distance(p))
    }

    /// Which side of the curve the point `p` of its plane lies on: the
    /// sign of its offset across a line, inside or outside a circle.
    fn side(self, p: Point, plane_normal: Vector) -> bool {
        match self {
            Self::Line { at, along } => along.cross(plane_normal).dot(p - at) > 0.0,
            Self::Circle {
                centre,
                normal,
                radius,
            } => {
                let q = p - centre;
                (q - normal * q.dot(normal)).magnitude() < radius
            }
        }
    }
}

/// A straight round whose contact on face `a` runs past the held edge
/// over `run`.
#[derive(Debug, Clone)]
pub(crate) struct Hold {
    curved: Curved,
    edge: Id,
    a: usize,
    b: usize,
    held: Id,
    ends: [usize; 2],
    /// The edge's start, its unit direction and its length.
    origin: Point,
    along: Vector,
    len: f64,
    /// The outward normals of `a` and `b`.
    na: Vector,
    nb: Vector,
    radius: f64,
    cliff: Cliff,
    /// Where along the edge the contact lies past the held edge.
    run: (f64, f64),
    /// Whether the run reaches one end of the band and not the other.
    pub(crate) at_end: bool,
}

/// A ball of the round: its centre, its contact on `b`, and its touch on
/// `a` or its foot on the held edge.
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

    /// Whether the band's contact on `a` at `x` lies past the held edge.
    fn past(&self, x: f64) -> bool {
        let touch = self.band(x) + self.na * self.radius;
        self.cliff.side(touch, self.na) != self.cliff.side(self.on_edge(x), self.na)
    }

    fn in_run(&self, x: f64) -> bool {
        x > self.run.0 && x < self.run.1
    }

    /// The ball at `x` along the spine: on both faces off the run; on it,
    /// on `b`'s offset at the radius from the held edge, found from the
    /// band's centre toward `a`'s plane.
    fn ball(&self, x: f64, hold: bool) -> Option<Ball> {
        let r = self.radius;
        let c = self.band(x);
        if !hold {
            return Some(Ball {
                centre: c,
                on_b: c + self.nb * r,
                on_a: c + self.na * r,
            });
        }
        let g = self.along.cross(self.nb);
        let g = if g.dot(self.na) < 0.0 { -g } else { g };
        let reach = r / g.dot(self.na);
        let f = |s: f64| self.cliff.distance(c + g * s).map(|d| d - r);
        let s = root(f, 0.0, reach)?;
        let centre = c + g * s;
        Some(Ball {
            centre,
            on_b: centre + self.nb * r,
            on_a: self.cliff.foot(centre)?,
        })
    }

    /// Where `p` lies along the held edge's curve: along a line, or the
    /// angle about a circle from `reference`.
    fn coordinate(&self, p: Point, reference: Vector) -> f64 {
        match self.cliff {
            Cliff::Line { at, along } => (p - at).dot(along),
            Cliff::Circle { centre, normal, .. } => {
                let q = p - centre;
                let side = normal.cross(reference);
                q.dot(side).atan2(q.dot(reference))
            }
        }
    }

    /// Whether every ball's foot over the run lies within one of `runs`,
    /// each the points along an edge on the held edge's curve.
    fn feet_on(&self, runs: &[Vec<Point>]) -> bool {
        let (lo, hi) = self.run;
        let Some(middle) = self.ball(0.5 * (lo + hi), true) else {
            return false;
        };
        let reference = match self.cliff {
            Cliff::Line { along, .. } => along,
            Cliff::Circle { centre, normal, .. } => {
                let q = middle.on_a - centre;
                let q = q - normal * q.dot(normal);
                q / q.magnitude().max(f64::MIN_POSITIVE)
            }
        };
        let coords: Vec<Vec<f64>> = runs
            .iter()
            .map(|r| r.iter().map(|p| self.coordinate(*p, reference)).collect())
            .collect();
        let slack = self.radius * 1e-6;
        (0..=32).all(|k| {
            let x = lo + (hi - lo) * f64::from(k) / 32.0;
            let Some(ball) = self.ball(x, true) else {
                return false;
            };
            let c = self.coordinate(ball.on_a, reference);
            coords.iter().any(|run| {
                run.windows(2).any(|w| {
                    (w[1] - w[0]).abs() < core::f64::consts::FRAC_PI_2
                        && c >= w[0].min(w[1]) - slack
                        && c <= w[0].max(w[1]) + slack
                })
            })
        })
    }

    fn ball_at(&self, x: f64) -> Option<Ball> {
        self.ball(x, self.in_run(x))
    }

    /// The point a fraction `u` along the arc of the ball at `t`, from its
    /// contact on `b` to its touch on `a` or the held edge.
    fn point(&self, t: f64, u: f64, hold: bool) -> Option<Point> {
        let ball = self.ball(t, hold)?;
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

    /// The canal's section by the plane square to the edge at `x`, at
    /// `n` points from its contact on `b` to its foot on the held edge:
    /// each point the one of its fraction along a ball's arc that lies in
    /// that plane, the ball found by bisection along the spine. The
    /// canal's balls run on past the run's ends, so a section at a seat
    /// change is the band's.
    fn section(&self, x: f64, n: usize) -> Option<Vec<Point>> {
        let mut row = Vec::with_capacity(n);
        for k in 0..n {
            #[expect(clippy::cast_precision_loss, reason = "a few dozen points")]
            let u = k as f64 / (n - 1) as f64;
            let gap = |t: f64| {
                self.point(t, u, true)
                    .map(|p| (p - self.origin).dot(self.along) - x)
            };
            // The nearest bracket that holds the plane: the canal's balls
            // reach only so far past the run.
            let mut bracket = None;
            for k in [0.02, 0.1, 0.3, 1.0, 3.0] {
                let span = k * self.radius;
                if let (Some(a), Some(b)) = (gap(x - span), gap(x + span))
                    && a <= 0.0
                    && b >= 0.0
                {
                    bracket = Some((x - span, x + span));
                    break;
                }
            }
            let (lo, hi) = bracket?;
            row.push(self.point(root(gap, lo, hi)?, u, true)?);
        }
        Some(row)
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

/// `faces` read for the held edge: the plane at index 0 whose edge `held`
/// the contact runs past, and the other.
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
            let normal = circle.frame().z().vector();
            (normal.cross(plane_normal).magnitude() <= eps).then_some(Cliff::Circle {
                centre: circle.centre(),
                normal,
                radius: circle.radius(),
            })
        }
        _ => None,
    }
}

/// The faces other than `a` and `b` at vertex `v`.
fn others_at(curved: &Curved, v: Id, a: usize, b: usize) -> Vec<usize> {
    let mut found: Vec<usize> = (0..curved.faces.len())
        .filter(|fi| *fi != a && *fi != b)
        .filter(|fi| {
            curved.faces[*fi].loops.iter().flatten().any(|(e, _)| {
                let c = &curved.edges[e];
                c.ends.0 == v || c.ends.1 == v
            })
        })
        .collect();
    found.sort_unstable();
    found.dedup();
    found
}

/// `edges` read as one straight round of `radius` between two planes
/// whose contact on one runs past a sharp edge of it, held as hold edge
/// holds it; `None` for any other request, or where anything else of the
/// solid stands within its reach.
///
/// # Errors
///
/// As reading the solid.
pub(crate) fn reach(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Hold>> {
    let eps = tol.confusion();
    let ang = tol.angular();
    let [edge] = edges else {
        return Ok(None);
    };
    if !radius.is_finite() || radius <= eps {
        return Ok(None);
    }
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
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    let eid = id(edge);
    let Some(e) = curved.edges.get(&eid) else {
        return Ok(None);
    };
    let Curve::Line(_) = &e.curve else {
        return Ok(None);
    };
    let (p0, p1) = (curved.points[&e.ends.0], curved.points[&e.ends.1]);
    let len = p0.distance(p1);
    if len <= eps {
        return Ok(None);
    }
    let along = (p1 - p0) / len;
    let [f0, f1] = faces_of(&curved, eid)[..] else {
        return Ok(None);
    };
    let normal = |fi: usize| outward(&curved, fi, p0);
    let (Some(n0), Some(n1)) = (normal(f0), normal(f1)) else {
        return Ok(None);
    };
    // The end faces: one plane square to the edge at each end.
    let mut ends = [0; 2];
    for (k, v) in [e.ends.0, e.ends.1].into_iter().enumerate() {
        let [end] = others_at(&curved, v, f0, f1)[..] else {
            return Ok(None);
        };
        let SurfaceGeometry::Plane(_) = &curved.faces[end].surface else {
            return Ok(None);
        };
        let Some(n) = outward(&curved, end, curved.points[&v]) else {
            return Ok(None);
        };
        if n.cross(along).magnitude() > ang {
            return Ok(None);
        }
        ends[k] = end;
    }
    let mut found: Vec<Hold> = Vec::new();
    for (a, b, na, nb) in [(f0, f1, n0, n1), (f1, f0, n1, n0)] {
        for (held, _) in curved.faces[a].loops.iter().flatten() {
            if *held == eid {
                continue;
            }
            if faces_of(&curved, *held)
                .iter()
                .any(|f| ends.contains(f) || *f == b)
            {
                continue;
            }
            let Some(cliff) = cliff_of(&curved, *held, na, ang) else {
                continue;
            };
            let hold = Hold {
                curved: curved.clone(),
                edge: eid,
                a,
                b,
                held: *held,
                ends,
                origin: p0,
                along,
                len,
                na,
                nb,
                radius,
                cliff,
                run: (0.0, 0.0),
                at_end: false,
            };
            // A run over the whole band would take face `a` whole: the
            // round is then too large for its seat, as the largest-radius
            // readings say, not held.
            let Some(run) = run_of(&hold).filter(|r| r.0 > eps || r.1 < len - eps) else {
                continue;
            };
            let hold = Hold {
                run,
                at_end: (run.0 <= eps) != (run.1 >= len - eps),
                ..hold
            };
            if feet_held(&hold, tol)? {
                found.push(hold);
            }
        }
    }
    let [hold] = &found[..] else {
        return Ok(None);
    };
    let held_edge = curved.shapes.get(&hold.held).cloned();
    let Some(held_edge) = held_edge else {
        return Ok(None);
    };
    // The held edge is sharp and convex: the face past it falls away.
    if crate::setback::edge_is_convex(model, solid, &held_edge, radius, tol)? != Some(true) {
        return Ok(None);
    }
    if !clear(model, hold, tol)? {
        return Ok(None);
    }
    Ok(Some(hold.clone()))
}

/// Whether the balls' feet over the run stay on the held edge itself, but
/// within the radius of an end face the run reaches, where the end face
/// cuts the canal.
fn feet_held(hold: &Hold, tol: Tolerances) -> OgeomResult<bool> {
    let eps = tol.confusion();
    let edge = &hold.curved.edges[&hold.held];
    let mut along = Vec::with_capacity(SAMPLES + 1);
    for k in 0..=SAMPLES {
        #[expect(clippy::cast_precision_loss, reason = "a few dozen samples")]
        let t = edge.range.0 + (edge.range.1 - edge.range.0) * k as f64 / SAMPLES as f64;
        along.push(edge.curve.point_at(t, tol)?);
    }
    let (lo, hi) = hold.run;
    let inner = Hold {
        run: (
            if lo <= eps { hold.radius } else { lo },
            if hi >= hold.len - eps {
                hold.len - hold.radius
            } else {
                hi
            },
        ),
        ..hold.clone()
    };
    Ok(inner.run.1 <= inner.run.0 || inner.feet_on(&[along]))
}

/// The one stretch of the edge where the band's contact lies past the held
/// edge, its ends found by bisection; `None` where there is none, or more
/// than one.
fn run_of(hold: &Hold) -> Option<(f64, f64)> {
    const N: usize = 512;
    #[expect(clippy::cast_precision_loss, reason = "a few hundred samples")]
    let at = |k: usize| hold.len * k as f64 / N as f64;
    let flags: Vec<bool> = (0..=N).map(|k| hold.past(at(k))).collect();
    let first = flags.iter().position(|f| *f)?;
    let last = flags.iter().rposition(|f| *f)?;
    if flags[first..=last].iter().any(|f| !f) {
        return None;
    }
    let edge = |inside: usize, outside: usize| {
        let (mut i, mut o) = (at(inside), at(outside));
        for _ in 0..80 {
            let m = 0.5 * (i + o);
            if hold.past(m) {
                i = m;
            } else {
                o = m;
            }
        }
        0.5 * (i + o)
    };
    let lo = if first == 0 {
        0.0
    } else {
        edge(first, first - 1)
    };
    let hi = if last == N {
        hold.len
    } else {
        edge(last, last + 1)
    };
    (hi - lo > hold.radius * 1e-3).then_some((lo, hi))
}

/// Whether nothing of the solid but the round's faces, its end faces and
/// the face past the held edge stands within the box the tool can take up,
/// and no edge of theirs but the round's, the held edge and those between
/// them.
fn clear(model: &Model, hold: &Hold, tol: Tolerances) -> OgeomResult<bool> {
    let eps = tol.confusion();
    let mut points = vec![hold.on_edge(0.0), hold.on_edge(hold.len)];
    for k in 0..=64 {
        let x = hold.len * f64::from(k) / 64.0;
        let Some(ball) = hold.ball_at(x) else {
            return Ok(false);
        };
        points.extend([ball.on_a, ball.on_b]);
    }
    let room = Aabb::of_points(&points);
    let margin = hold.radius * 1e-3 + eps * 10.0;
    let (Some(lo), Some(hi)) = (room.low(), room.high()) else {
        return Ok(false);
    };
    let pad = Vector::new(margin, margin, margin);
    let room = Aabb::of_points(&[lo - pad, hi + pad]);
    let c = &hold.curved;
    let [past] = faces_of(c, hold.held)
        .into_iter()
        .filter(|f| *f != hold.a)
        .collect::<Vec<_>>()[..]
    else {
        return Ok(false);
    };
    let ours = [hold.a, hold.b, hold.ends[0], hold.ends[1], past];
    for (fi, face) in c.faces.iter().enumerate() {
        if ours.contains(&fi) {
            continue;
        }
        if ogeom_algo::shape_bounds(model, &face.shape, tol)?.intersects(&room) {
            return Ok(false);
        }
    }
    for fi in ours {
        for (e, _) in c.faces[fi].loops.iter().flatten() {
            if *e == hold.edge || *e == hold.held {
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

/// The canal over `from..to` fitted through its planar sections, closer
/// together toward the ends.
fn canal(
    hold: &Hold,
    from: f64,
    to: f64,
    tol: Tolerances,
) -> OgeomResult<Option<ogeom_geom::BSplineSurface>> {
    let mut rows = Vec::with_capacity(STATIONS + 1);
    for i in 0..=STATIONS {
        #[expect(clippy::cast_precision_loss, reason = "a hundred stations")]
        let s = i as f64 / STATIONS as f64;
        let x = from + (to - from) * 0.5 * (1.0 - (core::f64::consts::PI * s).cos());
        let Some(row) = hold.section(x, ACROSS) else {
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
    let eps = tol.confusion();
    let r = hold.radius;
    let (lo, hi) = hold.run;
    // The breaks along the edge, and whether the piece after each is the
    // canal.
    let mut breaks = vec![0.0];
    let mut canal_after = Vec::new();
    if lo > eps {
        breaks.push(lo);
        canal_after.push(false);
    }
    canal_after.push(true);
    if hi < hold.len - eps {
        breaks.push(hi);
        canal_after.push(false);
    }
    breaks.push(hold.len);
    // The canal, once, over the run.
    let Some(surface) = canal(hold, lo, hi, tol)? else {
        return Ok(None);
    };
    let mut points = vec![hold.on_edge(0.0), hold.on_edge(hold.len)];
    // At each break, its contact on `b` and on `a` or the held edge.
    let mut at_b = Vec::with_capacity(breaks.len());
    let mut at_a = Vec::with_capacity(breaks.len());
    for &x in &breaks {
        let canal_end = (x <= eps && lo <= eps) || (x >= hold.len - eps && hi >= hold.len - eps);
        let (pb, pa) = if canal_end {
            let Some(row) = hold.section(x, 2) else {
                return Ok(None);
            };
            (row[0], row[1])
        } else {
            let Some(ball) = hold.ball(x, false) else {
                return Ok(None);
            };
            (ball.on_b, ball.on_a)
        };
        at_b.push(points.len());
        points.push(pb);
        at_a.push(points.len());
        points.push(pa);
    }
    let (e0, e1) = (0, 1);
    let last = breaks.len() - 1;
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
    // The section at each break: the band's arc, or the canal's end.
    let mut sections = Vec::with_capacity(breaks.len());
    for (k, &x) in breaks.iter().enumerate() {
        let ends = (at_b[k], at_a[k]);
        let canal_end = (k == 0 && lo <= eps) || (k == last && hi >= hold.len - eps);
        if canal_end {
            let Some(side) = boundary_from(&surface, &points, ends, tol)? else {
                return Ok(None);
            };
            sections.push(push(side.ends, (side.curve, side.range)));
        } else {
            let Some(ball) = hold.ball(x, false) else {
                return Ok(None);
            };
            let (u, w) = (ball.on_b - ball.centre, ball.on_a - ball.centre);
            let z = u.cross(w);
            let sweep = (u.dot(w) / (r * r)).clamp(-1.0, 1.0).acos();
            sections.push(push(ends, arc(ball.centre, z, u, r, sweep, tol)?));
        }
    }
    // Along each piece, its contact on `b` and on `a` or the held edge.
    let mut along_b = Vec::new();
    let mut along_a = Vec::new();
    for (k, canal_piece) in canal_after.iter().enumerate() {
        let (b0, b1, a0, a1) = (at_b[k], at_b[k + 1], at_a[k], at_a[k + 1]);
        if *canal_piece {
            let (Some(sb), Some(sa)) = (
                boundary_from(&surface, &points, (b0, b1), tol)?,
                boundary_from(&surface, &points, (a0, a1), tol)?,
            ) else {
                return Ok(None);
            };
            along_b.push(push(sb.ends, (sb.curve, sb.range)));
            along_a.push(push(sa.ends, (sa.curve, sa.range)));
        } else {
            along_b.push(push((b0, b1), segment(points[b0], points[b1], tol)?));
            along_a.push(push((a0, a1), segment(points[a0], points[a1], tol)?));
        }
    }
    let curved = &hold.curved;
    let mid = 0.5 * (breaks[0] + breaks[1]);
    let first_canal = canal_after[0];
    let (mb, ma) = if first_canal {
        let Some(row) = hold.section(mid, 2) else {
            return Ok(None);
        };
        (row[0], row[1])
    } else {
        let Some(ball) = hold.ball(mid, false) else {
            return Ok(None);
        };
        (ball.on_b, ball.on_a)
    };
    let on_mid = hold.on_edge(mid);
    let mut faces = vec![
        ToolFace {
            surface: curved.faces[hold.b].surface.clone(),
            curves: [vec![edge_line, b_end], along_b.clone(), vec![b_start]].concat(),
            at: on_mid + (mb - on_mid) * 0.5,
            out: hold.nb,
        },
        ToolFace {
            surface: curved.faces[hold.a].surface.clone(),
            curves: [vec![edge_line, a_end], along_a.clone(), vec![a_start]].concat(),
            at: on_mid + (ma - on_mid) * 0.5,
            out: hold.na,
        },
    ];
    for (k, end) in [(0, hold.ends[0]), (last, hold.ends[1])] {
        let pe = if k == 0 { points[e0] } else { points[e1] };
        let (pb, pa) = (points[at_b[k]], points[at_a[k]]);
        faces.push(ToolFace {
            surface: curved.faces[end].surface.clone(),
            curves: if k == 0 {
                vec![b_start, sections[0], a_start]
            } else {
                vec![b_end, sections[last], a_end]
            },
            at: pe + (pb - pe) * 0.2 + (pa - pe) * 0.2,
            out: if k == 0 { -hold.along } else { hold.along },
        });
    }
    let round_from = faces.len();
    for (k, canal_piece) in canal_after.iter().enumerate() {
        let x = 0.5 * (breaks[k] + breaks[k + 1]);
        let Some(ball) = hold.ball(x, *canal_piece) else {
            return Ok(None);
        };
        let Some(p) = hold.point(x, 0.5, *canal_piece) else {
            return Ok(None);
        };
        let surface = if *canal_piece {
            surface.clone().into()
        } else {
            let frame = Frame::new(
                hold.band(0.0),
                Direction::new(hold.along, tol)?,
                Direction::new(hold.nb, tol)?,
                tol,
            )?;
            CylinderSurface::new(Cylinder::new(frame, r, tol)?, (-r, hold.len + r))?.into()
        };
        faces.push(ToolFace {
            surface,
            curves: vec![along_b[k], sections[k + 1], along_a[k], sections[k]],
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

/// Whether the result's edges lying on the held edge's curve still hold
/// the balls' feet over the run.
fn kept(model: &Model, shape: &Shape, hold: &Hold, tol: Tolerances) -> OgeomResult<bool> {
    use ogeom_topo::{ShapeType, explore_unique};
    let near = tol.confusion() * 100.0;
    let mut runs: Vec<Vec<Point>> = Vec::new();
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
        let on = samples.iter().all(|p| {
            (*p - hold.origin).dot(hold.na).abs() <= near
                && hold.cliff.distance(*p).is_some_and(|d| d <= near)
        });
        if on {
            runs.push(samples);
        }
    }
    Ok(hold.feet_on(&runs))
}

/// Points along an edge, checked against the feet.
const SAMPLES: usize = 64;
