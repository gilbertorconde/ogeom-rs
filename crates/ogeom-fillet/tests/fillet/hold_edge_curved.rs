//! Hold edge where the ball runs onto a cylinder beside the edge, and
//! where the held edge is neither a line nor a circle (an ellipse, a
//! B-spline). The reference is the definition: the ball touches the front
//! and rests on the nearest point of the top's region.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use ogeom_core::Tolerances;
use ogeom_fillet::{Overflow, fillet_edges, fillet_edges_overflow};
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The round's radius on every part.
const R: f64 = 4.0;

/// The drum's axis `(y, z)`, along `x`, and its radius: through the edge
/// `y = z = 0`, the solid inside it.
const DRUM: (f64, f64, f64) = (8.0, -6.0, 10.0);

/// The pocket in the drum part's top: centre `(x, y)` and radius.
const POCKET: (f64, f64, f64) = (15.0, 6.5, 5.0);

/// The ellipse whose lower half is the elliptic part's cliff edge: centre
/// `(x, y)`, semi-axes along `x` and `y`.
const ELLIPSE: (f64, f64, f64, f64) = (10.0, 9.0, 14.0, 6.0);

/// The cubic Bezier in plan that is the spline part's cliff edge, carried
/// past both ends of the part.
const BEZIER: [(f64, f64); 4] = [(-1.0, 5.6), (6.0, 0.6), (14.0, 1.6), (21.0, 5.2)];

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn edges_where(m: &Model, s: &Shape, on: impl Fn(Point) -> bool) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, e, T).unwrap() else {
                return false;
            };
            (0..=10).all(|k| on(c.point_at(a + (b - a) * f64::from(k) / 10.0, T).unwrap()))
        })
        .collect()
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

fn surface_of(m: &Model, face: &Shape) -> SurfaceGeometry {
    let data = m.node(face).unwrap().data().as_face().unwrap();
    m.geometry().surface(data.surface).unwrap().clone()
}

/// The front top edge `y = z = 0` of `s`, the one edge asked.
fn front_edge(m: &Model, s: &Shape) -> Vec<Shape> {
    let edge = edges_where(m, s, |p| p.y.abs() < 1e-9 && p.z.abs() < 1e-9);
    assert_eq!(edge.len(), 1, "the front top edge");
    edge
}

/// A part 30 long whose top `z = 0` runs from the front top edge to a
/// back face at `y = 14`, its front the drum [`DRUM`] bulging out below
/// the edge, the bottom at `z = -12`; with `pocket`, a pocket [`POCKET`]
/// sunk 8 into its top.
fn drum(m: &mut Model, pocket: bool) -> Shape {
    let block = ogeom_algo::make_box(m, at(0.0, -5.0, -12.0), (30.0, 19.0, 12.0), T)
        .unwrap()
        .shape;
    let (ay, az, ar) = DRUM;
    let axis = Frame::new(Point::new(-1.0, ay, az), Direction::X, Direction::Y, T).unwrap();
    let drum = ogeom_algo::make_cylinder(m, axis, ar, 32.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::common(m, &block, &drum, T).unwrap().shape;
    if !pocket {
        return s;
    }
    let (px, py, pr) = POCKET;
    let cutter = ogeom_algo::make_cylinder(m, at(px, py, -8.0), pr, 9.0, T)
        .unwrap()
        .shape;
    ogeom_bool::cut(m, &s, &cutter, T).unwrap().shape
}

/// A block 20 by 12 by 15 under `z = 0`, its top lying in front of the
/// curve `edge` in plan (from `x = -1` to `x = 21`), a cliff 6 high falling
/// from that edge to a ledge. The cliff is cut by a prism standing past the
/// block's ends and back.
fn cliffed(m: &mut Model, edge: ogeom_geom::Curve, range: (f64, f64)) -> Shape {
    use ogeom_geom::{Curve, LineCurve, PlaneSurface};
    let block = ogeom_algo::make_box(m, at(0.0, 0.0, -15.0), (20.0, 12.0, 15.0), T)
        .unwrap()
        .shape;
    let floor = -6.0;
    let (p0, p1) = (
        edge.point_at(range.0, T).unwrap(),
        edge.point_at(range.1, T).unwrap(),
    );
    let p = |(x, y): (f64, f64)| Point::new(x, y, floor);
    let plan = [p0, p1, p((21.0, 13.0)), p((-1.0, 13.0))];
    let corners: Vec<Shape> = plan
        .iter()
        .map(|q| ogeom_algo::make_vertex(m, *q).shape)
        .collect();
    let mut edges = Vec::new();
    for k in 0..plan.len() {
        let (a, b) = (plan[k], plan[(k + 1) % plan.len()]);
        let (va, vb) = (&corners[k], &corners[(k + 1) % plan.len()]);
        let built = if k == 0 {
            ogeom_algo::make_edge_between(m, edge.clone(), range, va, vb, T)
        } else {
            let line = Curve::Line(LineCurve::segment(a, b, T).unwrap());
            ogeom_algo::make_edge_between(m, line, (0.0, a.distance(b)), va, vb, T)
        };
        edges.push(built.unwrap().shape);
    }
    let wire = ogeom_algo::make_wire(m, &edges, T).unwrap().shape;
    let plane = PlaneSurface::new(ogeom_math::Plane::new(at(0.0, 0.0, floor)));
    let face = ogeom_algo::make_face(m, SurfaceGeometry::Plane(plane), &[wire], T)
        .unwrap()
        .shape;
    let cutter = ogeom_algo::make_prism(m, &face, Vector::new(0.0, 0.0, 7.0), T)
        .unwrap()
        .shape;
    ogeom_bool::cut(m, &block, &cutter, T).unwrap().shape
}

/// [`cliffed`] behind the lower half of [`ELLIPSE`].
fn elliptic(m: &mut Model) -> Shape {
    let (cx, cy, ax, by) = ELLIPSE;
    let alpha = ((cx + 1.0) / ax).acos();
    let ellipse = ogeom_math::Ellipse::new(at(cx, cy, -6.0), ax, by, T).unwrap();
    cliffed(
        m,
        ogeom_geom::EllipseCurve::new(ellipse).into(),
        (PI + alpha, TAU - alpha),
    )
}

/// [`cliffed`] behind [`BEZIER`], a polynomial B-spline.
fn splined(m: &mut Model) -> Shape {
    let control = BEZIER
        .iter()
        .map(|(x, y)| Point::new(*x, *y, -6.0))
        .collect();
    let knots = ogeom_math::KnotVector::clamped_uniform(3, 4).unwrap();
    let curve = ogeom_geom::BSplineCurve::new(knots, control, T).unwrap();
    cliffed(m, ogeom_geom::Curve::BSpline(curve), (0.0, 1.0))
}

/// The top's region, in its plane `z = 0`: where it gives way past the
/// front top edge.
#[derive(Clone, Copy)]
enum Top {
    /// The top holds [`POCKET`].
    Pocket,
    /// The top lies in front of the lower half of [`ELLIPSE`].
    Ellipse,
    /// The top lies in front of [`BEZIER`].
    Bezier,
}

/// The cubic Bezier [`BEZIER`] at `u` and its first two derivatives.
fn bezier(u: f64) -> [(f64, f64); 3] {
    let [p0, p1, p2, p3] = BEZIER;
    let v = 1.0 - u;
    let mix = |f: &dyn Fn((f64, f64)) -> f64| {
        let (a, b, c, d) = (f(p0), f(p1), f(p2), f(p3));
        [
            v * v * v * a + 3.0 * v * v * u * b + 3.0 * v * u * u * c + u * u * u * d,
            3.0 * (v * v * (b - a) + 2.0 * v * u * (c - b) + u * u * (d - c)),
            6.0 * (v * (c - 2.0 * b + a) + u * (d - 2.0 * c + b)),
        ]
    };
    let (x, y) = (mix(&|p| p.0), mix(&|p| p.1));
    [(x[0], y[0]), (x[1], y[1]), (x[2], y[2])]
}

impl Top {
    /// The curve the top gives way along at `t`, and its first two
    /// derivatives: the ellipse's lower half by its angle, the Bezier by
    /// its parameter.
    fn curve(self, t: f64) -> [(f64, f64); 3] {
        match self {
            Self::Pocket => unreachable!("the pocket's distances are in closed form"),
            Self::Ellipse => {
                let (cx, cy, ax, by) = ELLIPSE;
                [
                    (cx + ax * t.cos(), cy + by * t.sin()),
                    (-ax * t.sin(), by * t.cos()),
                    (-ax * t.cos(), -by * t.sin()),
                ]
            }
            Self::Bezier => bezier(t),
        }
    }

    /// The parameter where the curve stands at `x`, by bisection.
    fn at_x(self, x: f64) -> f64 {
        let (mut lo, mut hi) = match self {
            Self::Ellipse => (PI, TAU),
            _ => (0.0, 1.0),
        };
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if self.curve(mid)[0].0 < x {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// How far behind the front the top runs at `x` before it gives way,
    /// up to `R` past the front's edge: the faces below its edges stand
    /// square to it, deeper than any ball reaches.
    fn wall(self, x: f64) -> f64 {
        match self {
            Self::Pocket => {
                let (px, py, pr) = POCKET;
                let d = x - px;
                if d.abs() >= pr {
                    f64::INFINITY
                } else {
                    py - (pr * pr - d * d).sqrt()
                }
            }
            _ => self.curve(self.at_x(x))[0].1,
        }
    }

    /// The distance in the top's plane from `(x, y)` to the top, and the
    /// nearest point of it: `(x, y)` itself on the top.
    fn nearest(self, x: f64, y: f64) -> (f64, (f64, f64)) {
        match self {
            Self::Pocket => {
                let (px, py, pr) = POCKET;
                let (dx, dy) = (x - px, y - py);
                let m = dx.hypot(dy);
                if m >= pr {
                    return (0.0, (x, y));
                }
                (pr - m, (px + dx * pr / m, py + dy * pr / m))
            }
            _ => {
                if y <= self.wall(x) {
                    return (0.0, (x, y));
                }
                // Newton on the curve from the point at the same `x`.
                let mut t = self.at_x(x);
                for _ in 0..50 {
                    let [(cx, cy), (dx, dy), (ex, ey)] = self.curve(t);
                    let g = (cx - x) * dx + (cy - y) * dy;
                    let slope = dx * dx + dy * dy + (cx - x) * ex + (cy - y) * ey;
                    let step = g / slope;
                    t -= step;
                    if step.abs() < 1e-15 {
                        break;
                    }
                }
                let (cx, cy) = self.curve(t)[0];
                ((cx - x).hypot(cy - y), (cx, cy))
            }
        }
    }
}

/// The face the ball touches throughout: the plane `y = 0` (material
/// behind) or the drum [`DRUM`].
#[derive(Clone, Copy)]
enum Front {
    Plane,
    Drum,
}

impl Front {
    fn inside(self, p: Point) -> bool {
        match self {
            Self::Plane => p.y >= 0.0,
            Self::Drum => {
                let (ay, az, ar) = DRUM;
                (p.y - ay).hypot(p.z - az) <= ar
            }
        }
    }

    /// Where the front stands at height `z`.
    fn at(self, z: f64) -> f64 {
        match self {
            Self::Plane => 0.0,
            Self::Drum => {
                let (ay, az, ar) = DRUM;
                ay - (ar * ar - (z - az) * (z - az)).sqrt()
            }
        }
    }
}

/// A ball of the reference: its centre, its contact on the front and its
/// touch on the top or its foot.
#[derive(Clone, Copy)]
struct Ball {
    centre: Point,
    on_b: Point,
    on_a: Point,
}

/// The independent reference for a round of radius `R` on the edge along
/// `x` between the top `z = 0` (material below) and the front. The ball
/// touches the front, and stands as near the top's plane as it may while
/// it stays `R` from the top's region: on the plane front at depth
/// `sqrt(R^2 - h^2)`, `h` the distance from its contact `(x, R)` to the
/// region in the top's plane; on the drum, its centre on the circle `R`
/// inside the drum, risen from the band's centre until its distance to the
/// region, `sqrt(h^2 + z^2)`, is `R`. A point of the solid is cut away
/// where it lies farther than `R` from every centre, between the arcs of
/// its nearest ball to its contact on the front and its touch on the top
/// or its foot.
struct Reference {
    top: Top,
    front: Front,
}

impl Reference {
    fn ball(&self, t: f64) -> Ball {
        match self.front {
            Front::Plane => {
                let (h, (fx, fy)) = self.top.nearest(t, R);
                let h = h.min(R);
                let centre = Point::new(t, R, -(R * R - h * h).sqrt());
                Ball {
                    centre,
                    on_b: Point::new(t, 0.0, centre.z),
                    on_a: Point::new(fx, fy, 0.0),
                }
            }
            Front::Drum => {
                let (ay, az, ar) = DRUM;
                let rho = ar - R;
                let place = |phi: f64| (ay + rho * phi.cos(), az + rho * phi.sin());
                let band = (2.0_f64).atan2(-(rho * rho - 4.0).sqrt());
                let gap = |phi: f64| {
                    let (y, z) = place(phi);
                    let (h, _) = self.top.nearest(t, y);
                    h * h + z * z - R * R
                };
                let phi = if gap(band) <= 0.0 {
                    band
                } else {
                    let (mut lo, mut hi) = (FRAC_PI_2, band);
                    for _ in 0..80 {
                        let mid = 0.5 * (lo + hi);
                        if gap(mid) > 0.0 {
                            hi = mid;
                        } else {
                            lo = mid;
                        }
                    }
                    0.5 * (lo + hi)
                };
                let (y, z) = place(phi);
                let (_, (fx, fy)) = self.top.nearest(t, y);
                Ball {
                    centre: Point::new(t, y, z),
                    on_b: Point::new(t, ay + ar * phi.cos(), az + ar * phi.sin()),
                    on_a: Point::new(fx, fy, 0.0),
                }
            }
        }
    }

    fn cut(&self, p: Point) -> bool {
        if p.z > 0.0 || !self.front.inside(p) {
            return false;
        }
        let d = |t: f64| p.distance(self.ball(t).centre);
        let (mut lo, mut hi) = (p.x - 2.0 * R, p.x + 2.0 * R);
        let mut best = p.x;
        let mut least = f64::INFINITY;
        for k in 0..=128 {
            let t = lo + (hi - lo) * f64::from(k) / 128.0;
            if d(t) < least {
                least = d(t);
                best = t;
            }
        }
        let step = (hi - lo) / 128.0;
        lo = best - step;
        hi = best + step;
        let g = 0.5 * (5.0_f64.sqrt() - 1.0);
        for _ in 0..60 {
            let (a, b) = (hi - g * (hi - lo), lo + g * (hi - lo));
            if d(a) < d(b) {
                hi = b;
            } else {
                lo = a;
            }
        }
        let ball = self.ball(0.5 * (lo + hi));
        let w = p - ball.centre;
        let (eb, ea) = (ball.on_b - ball.centre, ball.on_a - ball.centre);
        let n = eb.cross(ea);
        w.magnitude() > R && eb.cross(w).dot(n) >= 0.0 && w.cross(ea).dot(n) >= 0.0
    }

    /// How far behind the front the cut reaches at `(x, z)`, by bisection.
    fn depth(&self, x: f64, z: f64) -> f64 {
        let y0 = self.front.at(z);
        let (mut lo, mut hi) = (0.0, R);
        if !self.cut(Point::new(x, y0 + 1e-12, z)) {
            return 0.0;
        }
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            if self.cut(Point::new(x, y0 + mid, z)) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// The volume cut away: in `x` over `ranges`, each Gauss-Legendre of
    /// `n` points; down the front from the top to the ball's contact there,
    /// in `s` with `z = -s^2`, each column's depth as far as the solid runs
    /// (to the top's [`Top::wall`]). Where the cut reaches past that, the
    /// column is split at the `s` where it stops doing so.
    fn volume(&self, ranges: &[(f64, f64)], n: usize) -> f64 {
        let g = gauss(n);
        let mut sum = 0.0;
        for &(x0, x1) in ranges {
            let hx = 0.5 * (x1 - x0);
            for (gx, wx) in &g {
                let x = x0 + hx * (gx + 1.0);
                let s1 = (-self.ball(x).on_b.z).sqrt();
                let wall = self.top.wall(x);
                let room = |s: f64| wall - self.front.at(-s * s);
                let past = |s: f64| self.depth(x, -s * s) > room(s);
                let mut splits = vec![0.0];
                if past(0.0) && !past(s1) {
                    let (mut lo, mut hi) = (0.0, s1);
                    for _ in 0..50 {
                        let mid = 0.5 * (lo + hi);
                        if past(mid) {
                            lo = mid;
                        } else {
                            hi = mid;
                        }
                    }
                    splits.push(0.5 * (lo + hi));
                }
                splits.push(s1);
                for w in splits.windows(2) {
                    let hs = 0.5 * (w[1] - w[0]);
                    for (gs, ws) in &g {
                        let s = w[0] + hs * (gs + 1.0);
                        let d = self.depth(x, -s * s).min(room(s));
                        sum += wx * ws * hx * hs * d * 2.0 * s;
                    }
                }
            }
        }
        sum
    }

    /// Where the band's contact on the top, `(x, y)`, crosses the region's
    /// edge, by bisection between `lo` (on the top) and `hi` (past it).
    fn change(&self, lo: f64, hi: f64) -> f64 {
        let y = self.ball(lo).on_a.y;
        let on = |x: f64| self.top.nearest(x, y).0 == 0.0;
        let (mut a, mut b) = (lo, hi);
        assert!(on(a) && !on(b), "a seat change between {lo} and {hi}");
        for _ in 0..80 {
            let mid = 0.5 * (a + b);
            if on(mid) {
                a = mid;
            } else {
                b = mid;
            }
        }
        0.5 * (a + b)
    }
}

/// Gauss-Legendre nodes and weights on [-1, 1], by Newton on the
/// polynomial's roots.
fn gauss(n: usize) -> Vec<(f64, f64)> {
    #[expect(clippy::cast_precision_loss, reason = "a few dozen nodes")]
    let nf = n as f64;
    (0..n)
        .map(|i| {
            #[expect(clippy::cast_precision_loss, reason = "a few dozen nodes")]
            let mut x = (PI * (i as f64 + 0.75) / (nf + 0.5)).cos();
            let mut d = 0.0;
            for _ in 0..100 {
                let (mut p0, mut p1) = (1.0, x);
                for k in 2..=n {
                    #[expect(clippy::cast_precision_loss, reason = "a few dozen nodes")]
                    let kf = k as f64;
                    let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
                    p0 = p1;
                    p1 = p2;
                }
                d = nf * (x * p1 - p0) / (x * x - 1.0);
                let step = p1 / d;
                x -= step;
                if step.abs() < 1e-16 {
                    break;
                }
            }
            (x, 2.0 / ((1.0 - x * x) * d * d))
        })
        .collect()
}

/// The largest angle between the normals either side of each edge a fitted
/// face shares with another fitted face, a cylinder along `x` (the band or
/// the drum) or the front `y = 0`, sampled along the edge, and how many
/// such edges.
fn worst_join(m: &Model, s: &Shape) -> (f64, usize) {
    let mut worst = 0.0_f64;
    let mut joins = 0;
    for e in explore_unique(m, s, ShapeType::Edge).unwrap() {
        let faces = ancestors_of(m, s, &e, ShapeType::Face).unwrap();
        let [a, b] = &faces[..] else {
            continue;
        };
        let (sa, sb) = (surface_of(m, a), surface_of(m, b));
        let other = match (&sa, &sb) {
            (SurfaceGeometry::BSpline(_), o) | (o, SurfaceGeometry::BSpline(_)) => o,
            _ => continue,
        };
        let Some((c, (lo, hi))) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap() else {
            continue;
        };
        let mid = c.point_at(0.5 * (lo + hi), T).unwrap();
        let tangent_kind = match other {
            SurfaceGeometry::BSpline(_) => true,
            SurfaceGeometry::Cylinder(cy) => {
                cy.cylinder()
                    .axis()
                    .direction
                    .vector()
                    .cross(Vector::X)
                    .magnitude()
                    < 1e-9
            }
            SurfaceGeometry::Plane(_) => mid.y.abs() < 1e-6,
            _ => false,
        };
        if !tangent_kind {
            continue;
        }
        joins += 1;
        for k in 1..20 {
            let p = c.point_at(lo + (hi - lo) * f64::from(k) / 20.0, T).unwrap();
            let normal = |g: &SurfaceGeometry| -> Vector {
                let foot = ogeom_geom::project::project_on_surface(g, p, 16, T).unwrap();
                let (u, v) = foot.parameters;
                g.normal_at(u, v, T).unwrap().vector()
            };
            let angle = normal(&sa).dot(normal(&sb)).abs().clamp(0.0, 1.0).acos();
            worst = worst.max(angle);
        }
    }
    (worst, joins)
}

/// The result's edges on the top's plane within `near` of the curve the
/// top gives way along: how far, at most, any of their points strays from
/// that curve (a tenth of `near` at most for a kept edge), and the
/// x-range they cover behind the front by less than `R`.
fn kept(m: &Model, s: &Shape, top: Top, near: f64) -> (f64, (f64, f64)) {
    let off = |p: Point| {
        let (h, (fx, fy)) = top.nearest(p.x, p.y);
        // On the pocket's rim the distance from inside is `h`; on a cliff
        // the nearest point is found from either side.
        let planar = if h > 0.0 {
            h
        } else {
            (fx - p.x).hypot(fy - p.y).max(match top {
                Top::Pocket => {
                    let (px, py, pr) = POCKET;
                    (p.x - px).hypot(p.y - py) - pr
                }
                _ => p.y - top.wall(p.x),
            })
        };
        planar.abs().max(p.z.abs())
    };
    let mut worst = 0.0_f64;
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for e in edges_where(m, s, |p| off(p) < near) {
        let (c, (ta, tb)) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap().unwrap();
        for k in 0..=50 {
            let p = c.point_at(ta + (tb - ta) * f64::from(k) / 50.0, T).unwrap();
            worst = worst.max(off(p));
            if p.y < R + 1e-6 {
                lo = lo.min(p.x);
                hi = hi.max(p.x);
            }
        }
    }
    (worst, (lo, hi))
}

/// The checks every part shares: valid and closed, every join of a fitted
/// piece tangent within 1e-5 rad, and the curve the top gives way along
/// kept within 1e-9 from `from` to `to`; with `heavy`, the volume cut away
/// matching the reference split at `ranges` within its own bound (the
/// change from 24 to 48 points per range) plus 1e-6.
#[allow(clippy::too_many_arguments, reason = "one part's checks")]
fn holds(
    heavy: bool,
    m: &Model,
    before: f64,
    built: &ogeom_algo::Built,
    edge: &Shape,
    reference: &Reference,
    ranges: &[(f64, f64)],
    (from, to): (f64, f64),
) -> f64 {
    let found = ogeom_algo::check(m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    for shell in explore_unique(m, &built.shape, ShapeType::Shell).unwrap() {
        assert!(ogeom_algo::is_shell_closed(m, &shell).unwrap(), "closed");
    }
    assert!(!built.history.generated(edge).is_empty(), "the round");

    let (crease, joins) = worst_join(m, &built.shape);
    eprintln!("worst join {crease:.1e} rad over {joins} edges");
    assert!(joins >= 2, "the fitted pieces' joins");
    assert!(crease < 1e-5, "the fitted pieces' joins crease by {crease}");

    let (worst, (lo, hi)) = kept(m, &built.shape, reference.top, 1e-6);
    eprintln!("held edge kept over x {lo:.6} .. {hi:.6}, off its curve by {worst:.1e}");
    assert!(worst < 1e-9, "the held edge strays {worst} off its curve");
    assert!(
        lo <= from + 1e-6 && hi >= to - 1e-6,
        "the held edge kept over {lo}..{hi}, the feet over {from}..{to}"
    );
    let got = before - volume(m, &built.shape);
    if !heavy {
        return got;
    }
    let (coarse, fine) = (reference.volume(ranges, 24), reference.volume(ranges, 48));
    let bound = (fine - coarse).abs();
    eprintln!("cut away {got:.9}, reference {fine:.9}, bound {bound:.1e}");
    assert!(
        (got - fine).abs() < 1e-6 + bound,
        "cut away {got} against the reference {fine} (within {bound})"
    );
    got
}

/// `[a, b]` and the seat changes `c`, `d` inside it as ranges, each piece
/// halved.
fn split(a: f64, c: f64, d: f64, b: f64) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for (p, q) in [(a, c), (c, d), (d, b)] {
        let m = 0.5 * (p + q);
        out.extend([(p, m), (m, q)]);
    }
    out
}

/// The x-range where a ball rests on the held edge, and the band's piece
/// ends either side, by the reference.
fn run(reference: &Reference, len: f64) -> (f64, f64) {
    let mid = 0.5 * len;
    (
        reference.change(0.0, mid),
        // From the far end, back toward the middle.
        {
            let y = reference.ball(len).on_a.y;
            let on = |x: f64| reference.top.nearest(x, y).0 == 0.0;
            let (mut a, mut b) = (len, mid);
            assert!(on(a) && !on(b), "a seat change before the end");
            for _ in 0..80 {
                let m = 0.5 * (a + b);
                if on(m) {
                    a = m;
                } else {
                    b = m;
                }
            }
            0.5 * (a + b)
        },
    )
}

/// The ball runs onto a cylinder: the front top edge of the drum part
/// rounded at 4, its front the drum bulging below it. The band's centre
/// stands 4 under the top and 6 from the drum's axis, its contact on the
/// top 2.34 behind the edge; the pocket's rim comes within 1.5 of the
/// edge, so between `x = 12.22` and `x = 17.78` the ball rests on the rim,
/// its centre risen round the drum's axis. The rim is kept on its circle,
/// and the fitted piece meets the band and the drum G1.
fn the_ball_runs_onto_a_drum_case(heavy: bool) -> f64 {
    let mut m = Model::new();
    let s = drum(&mut m, true);
    let edge = front_edge(&m, &s);
    let before = volume(&m, &s);
    let built = fillet_edges(&mut m, &s, &edge, R, T).unwrap();
    let reference = Reference {
        top: Top::Pocket,
        front: Front::Drum,
    };
    let (c, d) = run(&reference, 30.0);
    eprintln!("held over {c} .. {d}");
    holds(
        heavy,
        &m,
        before,
        &built,
        &edge[0],
        &reference,
        &split(0.0, c, d, 30.0),
        (c, d),
    )
}

#[test]
fn hold_edge_holds_a_rim_beside_a_drum() {
    the_ball_runs_onto_a_drum_case(false);
}

#[test]
#[ignore = "heavy"]
fn a_drum_beside_the_edge_cuts_away_the_reference() {
    the_ball_runs_onto_a_drum_case(true);
}

/// The reference reproduces the drum's band in closed form: far from the
/// pocket the cut's section is the triangle from the edge to the ball's
/// two contacts, less the ball's segment over the chord between them, plus
/// the drum's segment over the chord from the edge to its contact.
#[test]
fn the_drum_reference_reads_the_band() {
    let reference = Reference {
        top: Top::Pocket,
        front: Front::Drum,
    };
    let ball = reference.ball(-20.0);
    let (ar, a, b) = (DRUM.2, ball.on_a, ball.on_b);
    let triangle = 0.5 * (a.y * b.z - b.y * a.z).abs();
    let segment = |radius: f64, chord: f64| {
        let theta = 2.0 * (0.5 * chord / radius).asin();
        0.5 * radius * radius * (theta - theta.sin())
    };
    let section =
        triangle - segment(R, a.distance(b)) + segment(ar, b.distance(Point::new(-20.0, 0.0, 0.0)));
    let sampled = reference.volume(&[(-21.0, -19.0)], 24) / 2.0;
    eprintln!("band section {section:.12}, sampled {sampled:.12}");
    assert!(
        (section - sampled).abs() < 1e-9,
        "{section} against {sampled}"
    );
}

/// A held edge that is an ellipse: the front top edge of the elliptic part
/// rounded at 4. Between `x = 2.26` and `x = 17.74` the band's contact
/// `(x, 4)` lies past the cliff edge, and the ball rests on it. The
/// ellipse is kept on its curve, and the fitted piece meets the band and
/// the front G1.
fn an_elliptic_held_edge_case(heavy: bool) -> f64 {
    let mut m = Model::new();
    let s = elliptic(&mut m);
    general_case(heavy, &mut m, &s, Top::Ellipse)
}

/// A held edge that is a polynomial B-spline: [`an_elliptic_held_edge_case`]
/// behind [`BEZIER`].
fn a_spline_held_edge_case(heavy: bool) -> f64 {
    let mut m = Model::new();
    let s = splined(&mut m);
    general_case(heavy, &mut m, &s, Top::Bezier)
}

fn general_case(heavy: bool, m: &mut Model, s: &Shape, top: Top) -> f64 {
    let edge = front_edge(m, s);
    let before = volume(m, s);
    let built = fillet_edges(m, s, &edge, R, T).unwrap();
    let reference = Reference {
        top,
        front: Front::Plane,
    };
    let (c, d) = run(&reference, 20.0);
    eprintln!("held over {c} .. {d}");
    holds(
        heavy,
        m,
        before,
        &built,
        &edge[0],
        &reference,
        &split(0.0, c, d, 20.0),
        (c, d),
    )
}

#[test]
fn hold_edge_holds_an_ellipse() {
    an_elliptic_held_edge_case(false);
}

#[test]
#[ignore = "heavy"]
fn an_elliptic_held_edge_cuts_away_the_reference() {
    an_elliptic_held_edge_case(true);
}

#[test]
fn hold_edge_holds_a_spline() {
    a_spline_held_edge_case(false);
}

#[test]
#[ignore = "heavy"]
fn a_spline_held_edge_cuts_away_the_reference() {
    a_spline_held_edge_case(true);
}

/// With hold edge off the parts go on to trim through as before: the drum
/// part's round notches the rim (cutting more than the held round), the
/// elliptic part's is refused by trim through by name; with trim through
/// off too, either is refused naming hold edge.
#[test]
fn hold_edge_off_trims_through_as_before() {
    let mut off = Overflow::default();
    off.hold_edge = false;
    let mut neither = off;
    neither.trim_through = false;

    let mut m = Model::new();
    let s = drum(&mut m, true);
    let edge = front_edge(&m, &s);
    let before = volume(&m, &s);
    let held = fillet_edges(&mut m, &s, &edge, R, T).unwrap().shape;
    let held = before - volume(&m, &held);
    let trimmed = fillet_edges_overflow(&mut m, &s, &edge, R, off, T).unwrap();
    assert!(ogeom_algo::check(&m, &trimmed.shape, T).unwrap().is_valid());
    let notch = before - volume(&m, &trimmed.shape);
    eprintln!("held {held:.6}, trimmed through {notch:.6}");
    assert!(notch > held + 1e-2, "{notch} against {held}");
    // The notch takes the rim where the ball would rest on it.
    let (px, py, pr) = POCKET;
    let on_rim = |p: Point| p.z.abs() < 1e-6 && ((p.x - px).hypot(p.y - py) - pr).abs() < 1e-6;
    let left = edges_where(&m, &trimmed.shape, on_rim).iter().any(|e| {
        let (c, (a, b)) = ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap();
        (0..=50).any(|k| {
            let p = c.point_at(a + (b - a) * f64::from(k) / 50.0, T).unwrap();
            p.y < py && (p.x - px).abs() < 1.0
        })
    });
    assert!(!left, "trim through keeps the rim in front of the pocket");
    let err = fillet_edges_overflow(&mut m, &s, &edge, R, neither, T).unwrap_err();
    assert!(err.to_string().contains("hold edge is off"), "{err}");

    for make in [elliptic, splined] {
        let mut m = Model::new();
        let s = make(&mut m);
        let edge = front_edge(&m, &s);
        let err = fillet_edges_overflow(&mut m, &s, &edge, R, off, T).unwrap_err();
        assert!(err.to_string().contains("sets back"), "{err}");
        let err = fillet_edges_overflow(&mut m, &s, &edge, R, neither, T).unwrap_err();
        assert!(err.to_string().contains("hold edge is off"), "{err}");
    }
}

/// The same request with hold edge on and off: the same solid's volume, or
/// the same refusal.
fn as_with_hold_edge_off(m: &mut Model, s: &Shape, edge: &[Shape]) {
    let mut off = Overflow::default();
    off.hold_edge = false;
    let on = fillet_edges(m, s, edge, R, T).map(|b| volume(m, &b.shape));
    let without = fillet_edges_overflow(m, s, edge, R, off, T).map(|b| volume(m, &b.shape));
    eprintln!("on {on:?}, off {without:?}");
    match (on, without) {
        (Ok(a), Ok(b)) => assert!((a - b).abs() < 1e-9, "{a} against {b}"),
        (Err(a), Err(b)) => assert_eq!(a.to_string(), b.to_string()),
        (a, b) => panic!("hold edge on gives {a:?}, off {b:?}"),
    }
}

/// A cone beside the edge (the edge a ruling of it) is not held: its offset
/// and the top's meet in a conic, not a line along the edge, so the band
/// is not straight. The request goes on to trim through, as with hold edge
/// off.
#[test]
fn a_cone_beside_the_edge_goes_on_to_trim_through() {
    let mut m = Model::new();
    let block = ogeom_algo::make_box(&mut m, at(0.0, -5.0, -12.0), (30.0, 19.0, 12.0), T)
        .unwrap()
        .shape;
    // The edge is the ruling through the apex at `x = -40`; the axis leans
    // a quarter's slope from it, toward the drum's axis.
    let slope: f64 = 0.25;
    let beta = slope.atan();
    let (ay, az, ar) = DRUM;
    let lean = Vector::new(0.0, ay, az) / ar;
    let axis = Vector::X * beta.cos() + lean * beta.sin();
    let apex = Point::new(-40.0, 0.0, 0.0);
    let (s0, s1) = (36.0, 75.0);
    let frame = Frame::new(
        apex + axis * s0,
        Direction::new(axis, T).unwrap(),
        Direction::new(axis.cross(Vector::Y), T).unwrap(),
        T,
    )
    .unwrap();
    let cone = ogeom_algo::make_cone(&mut m, frame, s0 * slope, s1 * slope, s1 - s0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::common(&mut m, &block, &cone, T).unwrap().shape;
    let (px, py, pr) = POCKET;
    let cutter = ogeom_algo::make_cylinder(&mut m, at(px, py, -8.0), pr, 9.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(&mut m, &s, &cutter, T).unwrap().shape;
    let edge = front_edge(&m, &s);
    as_with_hold_edge_off(&mut m, &s, &edge);
}

/// A held edge on a curved face is not held: a bore into the drum's front
/// whose rim the band's contact on the drum runs past. The rim is not a
/// curve in a plane, so no ball rests on it; the request goes on to trim
/// through, as with hold edge off.
#[test]
fn a_rim_on_the_drum_goes_on_to_trim_through() {
    let mut m = Model::new();
    let s = drum(&mut m, false);
    let bore = Frame::new(Point::new(15.0, -5.0, -4.5), Direction::Y, Direction::X, T).unwrap();
    let cutter = ogeom_algo::make_cylinder(&mut m, bore, 2.0, 10.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(&mut m, &s, &cutter, T).unwrap().shape;
    let edge = front_edge(&m, &s);
    as_with_hold_edge_off(&mut m, &s, &edge);
}
