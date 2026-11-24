//! Blend into round beyond a cylinder's foot on one convex round: a foot
//! round reaching two rounds at a corner, a conical boss's foot, and a
//! foot reaching a concave round. Each is measured against the
//! morphological closing, sampled without the construction.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use ogeom_core::Tolerances;
use ogeom_fillet::{Overflow, fillet_edges, fillet_edges_overflow};
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The foot round's radius on every part.
const FOOT: f64 = 3.0;

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn edges_where(m: &Model, s: &Shape, on: impl Fn(Point) -> bool) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            [0.0, 0.3, 0.7, 1.0]
                .iter()
                .all(|t| on(c.point_at(a + (b - a) * t, T).unwrap()))
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

/// The rim where the floor at height `z` meets a boss whose foot circle
/// has radius `r` about `(x, y)`.
fn foot_rim(m: &Model, s: &Shape, (x, y, z): (f64, f64, f64), r: f64) -> Vec<Shape> {
    let rim = edges_where(m, s, |p| {
        (p.z - z).abs() < 1e-9 && ((p.x - x).hypot(p.y - y) - r).abs() < 1e-9
    });
    assert!(!rim.is_empty(), "the foot rim");
    rim
}

/// A block 40 by 30 by 15 with the top edges along `y = 0`, and along
/// `x = 0` too where `corner`, rounded at 5.
fn rounded_block(m: &mut Model, corner: bool) -> Shape {
    let block = ogeom_algo::make_box(m, at(0.0, 0.0, 0.0), (40.0, 30.0, 15.0), T)
        .unwrap()
        .shape;
    let top = edges_where(m, &block, |p| {
        (p.y.abs() < 1e-9 || (corner && p.x.abs() < 1e-9)) && (p.z - 15.0).abs() < 1e-9
    });
    fillet_edges(m, &block, &top, 5.0, T).unwrap().shape
}

/// A straight round of the floor, seen from above: a point of its tangent
/// line, the direction square to it past the line (away from the boss),
/// its radius, and whether it is concave.
#[derive(Clone, Copy)]
struct Round {
    p: (f64, f64),
    n: (f64, f64),
    radius: f64,
    concave: bool,
}

impl Round {
    fn past(&self, x: f64, y: f64) -> f64 {
        (x - self.p.0).mul_add(self.n.0, (y - self.p.1) * self.n.1)
    }
}

/// The independent reference: the fill a ball of radius `FOOT` leaves in a
/// boss's foot, as the morphological closing defines it, in the frame of
/// the floor (height `z` over it). A point outside the solid is in the fill
/// where its distance from the set `Free` of centres clear of the solid is
/// more than the radius. `Free` is bounded near the foot by the wall's
/// offset (a cone, or a cylinder, about the boss's axis) and the floor's
/// offset (the plane at the radius, the cylinders about the rounds' axes
/// past their tangent lines), which meet along the spine. The distance is
/// the least of the projections onto either boundary that land on the
/// other's side, and the distance from the spine, by sampling and
/// golden-section search.
struct Closing {
    /// The boss's axis on the floor.
    axis: (f64, f64),
    /// The foot circle's radius.
    foot: f64,
    /// The wall's radius change per unit of height.
    slope: f64,
    rounds: Vec<Round>,
    /// Where `rounds` holds two convex rounds of one radius, square to each
    /// other, that meet round the top's corner: the corner's radius seen
    /// from above (nought where a sphere rounds it). The top is then the
    /// set of points within the rounds' radius of the prism under the
    /// floor's region (the half planes inside both tangent lines, its
    /// corner rounded at this radius) sunk by that radius, so a point's
    /// distance from it comes from its distance from that region.
    corner: Option<f64>,
}

impl Closing {
    /// Over a rounded corner's top: where `(x, y)` stands outside the
    /// floor's region, how far, and its foot on the region's outline.
    fn outside(&self, x: f64, y: f64) -> Option<(f64, (f64, f64))> {
        let rho = self.corner?;
        let [a, b] = [self.rounds[0], self.rounds[1]];
        // The region shrunk by `rho` is the quadrant inside both lines
        // moved in by `rho`; its foot pulls the point back square to each.
        let (sa, sb) = ((a.past(x, y) + rho).max(0.0), (b.past(x, y) + rho).max(0.0));
        let (fx, fy) = (
            sa.mul_add(-a.n.0, sb.mul_add(-b.n.0, x)),
            sa.mul_add(-a.n.1, sb.mul_add(-b.n.1, y)),
        );
        let to = sa.hypot(sb);
        if to <= rho {
            return None;
        }
        let k = rho / to;
        Some((to - rho, ((x - fx).mul_add(k, fx), (y - fy).mul_add(k, fy))))
    }

    /// The rounds' radius at a rounded corner.
    fn corner_radius(&self) -> f64 {
        self.rounds[0].radius
    }

    fn rho(&self, x: f64, y: f64) -> f64 {
        (x - self.axis.0).hypot(y - self.axis.1)
    }

    /// The round past whose tangent line `(x, y)` lies, and how far past.
    fn past(&self, x: f64, y: f64) -> Option<(Round, f64)> {
        self.rounds
            .iter()
            .map(|r| (*r, r.past(x, y)))
            .filter(|(_, s)| *s > 0.0)
            .max_by(|a, b| a.1.total_cmp(&b.1))
    }

    /// The floor's height at `(x, y)`.
    fn floor(&self, x: f64, y: f64) -> f64 {
        if self.corner.is_some() {
            let a = self.corner_radius();
            return self.outside(x, y).map_or(0.0, |(d, _)| {
                let d = d.min(a);
                a.mul_add(a, -d * d).sqrt() - a
            });
        }
        match self.past(x, y) {
            None => 0.0,
            Some((r, s)) => {
                let s = s.min(r.radius);
                let rise = r.radius - r.radius.mul_add(r.radius, -s * s).sqrt();
                if r.concave { rise } else { -rise }
            }
        }
    }

    /// A centre's clearance from the floor less the radius.
    fn above(&self, c: Point) -> f64 {
        if self.corner.is_some() {
            let a = self.corner_radius();
            return self
                .outside(c.x, c.y)
                .map_or(c.z - FOOT, |(d, _)| d.hypot(c.z + a) - a - FOOT);
        }
        match self.past(c.x, c.y) {
            None => c.z - FOOT,
            Some((r, s)) if r.concave => r.radius - s.hypot(c.z - r.radius) - FOOT,
            Some((r, s)) => s.hypot(c.z + r.radius) - r.radius - FOOT,
        }
    }

    /// A centre's clearance from the wall less the radius.
    fn beside(&self, c: Point) -> f64 {
        self.slope.mul_add(-c.z, self.rho(c.x, c.y) - self.foot) / self.slope.hypot(1.0) - FOOT
    }

    fn free(&self, c: Point) -> bool {
        self.beside(c) >= -1e-12 && self.above(c) >= -1e-12
    }

    /// The wall's offset: its distance from the axis at height `z`.
    fn offset_radius(&self, z: f64) -> f64 {
        self.slope
            .mul_add(z, FOOT.mul_add(self.slope.hypot(1.0), self.foot))
    }

    /// The floor's offset: a centre's height at `(x, y)`.
    fn offset_height(&self, x: f64, y: f64) -> f64 {
        if self.corner.is_some() {
            let a = self.corner_radius();
            return self
                .outside(x, y)
                .map_or(FOOT, |(d, _)| (a + FOOT).powi(2).sub_sq(d).sqrt() - a);
        }
        match self.past(x, y) {
            None => FOOT,
            Some((r, s)) if r.concave => r.radius - (r.radius - FOOT).powi(2).sub_sq(s).sqrt(),
            Some((r, s)) => (r.radius + FOOT).powi(2).sub_sq(s).sqrt() - r.radius,
        }
    }

    /// The spine at angle `t` about the axis, where the two offsets meet.
    fn spine(&self, t: f64) -> Point {
        let mut z = FOOT;
        let mut q = Point::new(0.0, 0.0, 0.0);
        for _ in 0..100 {
            let rho = self.offset_radius(z);
            let (x, y) = (
                rho.mul_add(t.cos(), self.axis.0),
                rho.mul_add(t.sin(), self.axis.1),
            );
            let next = self.offset_height(x, y);
            q = Point::new(x, y, next);
            if (next - z).abs() <= 1e-15 {
                break;
            }
            z = next;
        }
        q
    }

    fn distance_to_free(&self, p: Point) -> f64 {
        if self.free(p) {
            return 0.0;
        }
        let mut best = f64::INFINITY;
        let rho = self.rho(p.x, p.y);
        let (ux, uy) = ((p.x - self.axis.0) / rho, (p.y - self.axis.1) / rho);
        // Onto the wall's offset, in the plane through the axis.
        let k = self.slope;
        let r0 = self.offset_radius(0.0);
        let t = (rho - r0).mul_add(k, p.z) / k.mul_add(k, 1.0);
        let along = k.mul_add(t, r0);
        let q = Point::new(
            along.mul_add(ux, self.axis.0),
            along.mul_add(uy, self.axis.1),
            t,
        );
        if self.above(q) >= -1e-12 {
            best = best.min(p.distance(q));
        }
        // Onto the floor's offset over a rounded corner: from the nearest
        // point of the sunk prism, out to the rounds' radius and the ball's.
        if self.corner.is_some() {
            let a = self.corner_radius();
            let q = match self.outside(p.x, p.y) {
                None => Point::new(p.x, p.y, FOOT),
                Some((_, (fx, fy))) => {
                    let m = Point::new(fx, fy, -a);
                    m + (p - m) * ((a + FOOT) / p.distance(m))
                }
            };
            if self.beside(q) >= -1e-12 {
                best = best.min(p.distance(q));
            }
        }
        // Onto the floor's offset: the plane, and each round's cylinder.
        let flat = Point::new(p.x, p.y, FOOT);
        if self.corner.is_none() && self.past(p.x, p.y).is_none() && self.beside(flat) >= -1e-12 {
            best = best.min(p.distance(flat));
        }
        for r in self.rounds.iter().filter(|_| self.corner.is_none()) {
            let s = r.past(p.x, p.y);
            let (centre, size) = if r.concave {
                (r.radius, r.radius - FOOT)
            } else {
                (-r.radius, r.radius + FOOT)
            };
            let (ds, dz) = (s, p.z - centre);
            let m = ds.hypot(dz);
            if m <= 0.0 {
                continue;
            }
            let (qs, qz) = (ds * size / m, centre + dz * size / m);
            let q = Point::new(
                (qs - s).mul_add(r.n.0, p.x),
                (qs - s).mul_add(r.n.1, p.y),
                qz,
            );
            let own = self.past(q.x, q.y).is_some_and(|(o, _)| o.p == r.p);
            if qs >= 0.0 && own && self.beside(q) >= -1e-12 {
                best = best.min(p.distance(q));
            }
        }
        // The spine, near the point's own angle.
        let t0 = (p.y - self.axis.1).atan2(p.x - self.axis.0);
        let d = |t: f64| p.distance(self.spine(t));
        let (mut lo, mut hi) = (t0 - 1.2, t0 + 1.2);
        let mut at = t0;
        let mut least = f64::INFINITY;
        for k in 0..=96 {
            let t = lo + (hi - lo) * f64::from(k) / 96.0;
            if d(t) < least {
                least = d(t);
                at = t;
            }
        }
        lo = at - 2.4 / 96.0;
        hi = at + 2.4 / 96.0;
        let g = 0.5 * (5.0_f64.sqrt() - 1.0);
        for _ in 0..60 {
            let (a, b) = (hi - g * (hi - lo), lo + g * (hi - lo));
            if d(a) < d(b) {
                hi = b;
            } else {
                lo = a;
            }
        }
        best.min(d(0.5 * (lo + hi)))
    }

    /// The point at `rho` from the axis at angle `t` on the solid's
    /// surface: on the wall below the wall's foot circle, on the floor
    /// beyond it.
    fn ground(&self, rho: f64, t: f64) -> Point {
        let (x, y) = (
            rho.mul_add(t.cos(), self.axis.0),
            rho.mul_add(t.sin(), self.axis.1),
        );
        let z = if rho < self.foot {
            (rho - self.foot) / self.slope
        } else {
            self.floor(x, y)
        };
        Point::new(x, y, z)
    }

    fn in_fill(&self, p: Point) -> bool {
        self.distance_to_free(p) > FOOT
    }

    /// The fill's depth over the solid at `rho`, `t`, by bisection on the
    /// height where the distance from `Free` comes down to the radius.
    fn column(&self, rho: f64, t: f64) -> f64 {
        let g = self.ground(rho, t);
        if !self.in_fill(g) {
            return 0.0;
        }
        let (mut lo, mut hi) = (0.0, 3.0 * FOOT);
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            if self.in_fill(g + Vector::new(0.0, 0.0, mid)) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Where along the radius at `t` the fill ends on the solid's surface,
    /// between `inside` (in the fill) and `outside`. Past that point the
    /// surface lies at exactly the radius from `Free`, so the test asks for
    /// a margin: where the fill ends tangentially, its depth a step `d`
    /// short of the end is about `d^2 / 2r`, so the margin of 1e-10 moves the
    /// end by 2.5e-5 and the volume by far less.
    fn edge(&self, t: f64, mut inside: f64, mut outside: f64) -> f64 {
        for _ in 0..50 {
            let mid = 0.5 * (inside + outside);
            if self.distance_to_free(self.ground(mid, t)) > FOOT + 1e-10 {
                inside = mid;
            } else {
                outside = mid;
            }
        }
        0.5 * (inside + outside)
    }

    /// The fill's volume: at each angle about the axis, the columns
    /// integrated out from where the fill starts on the wall to where it
    /// ends on the floor, each angular range split where the ball's seat
    /// changes; Gauss-Legendre of `n` points per range. Over a cylinder's
    /// foot the columns run `rho = out - (out - foot) cos s`, smooth where
    /// the round meets the wall; up a cone's wall and out along the floor
    /// they are split at the foot circle.
    /// Over a rounded corner: the angles about the axis at which the ball's
    /// centre, on the circle of radius `reach`, changes seat (between the
    /// floor, either round and the corner), found by sampling and
    /// bisection.
    fn corner_cuts(&self, reach: f64) -> Vec<f64> {
        let seat = |t: f64| -> u8 {
            let (x, y) = (
                reach.mul_add(t.cos(), self.axis.0),
                reach.mul_add(t.sin(), self.axis.1),
            );
            let rho = self.corner.unwrap_or(0.0);
            if self.outside(x, y).is_none() {
                return 0;
            }
            match (
                self.rounds[0].past(x, y) + rho > 0.0,
                self.rounds[1].past(x, y) + rho > 0.0,
            ) {
                (true, true) => 3,
                (true, false) => 1,
                _ => 2,
            }
        };
        let samples = 4096;
        let mut cuts = Vec::new();
        for k in 0..samples {
            let (mut a, mut b) = (
                TAU * f64::from(k) / f64::from(samples),
                TAU * f64::from(k + 1) / f64::from(samples),
            );
            let held = seat(a);
            if seat(b) == held {
                continue;
            }
            for _ in 0..80 {
                let mid = 0.5 * (a + b);
                if seat(mid) == held {
                    a = mid;
                } else {
                    b = mid;
                }
            }
            cuts.push(0.5 * (a + b));
        }
        cuts
    }

    fn volume(&self, n: usize) -> f64 {
        self.volume_within(n, None)
    }

    /// [`Closing::volume`] within the wedge of angles from `from` to `to`
    /// about the axis, where `within` gives one: the fill of an open arc
    /// whose ends' planes through the axis cut the round.
    fn volume_within(&self, n: usize, within: Option<(f64, f64)>) -> f64 {
        let g = gauss(n);
        let reach = self.offset_radius(FOOT);
        let mut cuts = if self.corner.is_some() {
            self.corner_cuts(reach)
        } else {
            Vec::new()
        };
        for r in self.rounds.iter().filter(|_| self.corner.is_none()) {
            let gap = -r.past(self.axis.0, self.axis.1);
            if gap < reach {
                let mid = r.n.1.atan2(r.n.0);
                let half = (gap / reach).acos();
                cuts.extend([mid - half, mid + half]);
            }
        }
        if cuts.is_empty() {
            cuts.push(0.0);
        }
        let ranges: Vec<(f64, f64)> = if let Some((from, to)) = within {
            let mut c: Vec<f64> = cuts
                .iter()
                .map(|a| from + (a - from).rem_euclid(TAU))
                .filter(|a| *a < to)
                .collect();
            c.extend([from, to]);
            c.sort_by(f64::total_cmp);
            c.windows(2).map(|w| (w[0], w[1])).collect()
        } else {
            let mut c: Vec<f64> = cuts.iter().map(|a| a.rem_euclid(TAU)).collect();
            c.sort_by(f64::total_cmp);
            (0..c.len())
                .map(|i| (c[i], c.get(i + 1).copied().unwrap_or(c[0] + TAU)))
                .collect()
        };
        let mut sum = 0.0;
        for (a0, a1) in ranges {
            let ha = 0.5 * (a1 - a0);
            for (ga, wa) in &g {
                let t = a0 + ha * (ga + 1.0);
                let out = self.edge(t, self.foot, reach + FOOT);
                let mut ring = 0.0;
                if self.slope == 0.0 {
                    for (gs, ws) in &g {
                        let s = FRAC_PI_2 * 0.5 * (gs + 1.0);
                        let rho = (out - self.foot).mul_add(-s.cos(), out);
                        let drho = (out - self.foot) * s.sin() * FRAC_PI_2 * 0.5;
                        ring += ws * self.column(rho, t) * rho * drho;
                    }
                } else {
                    let low = self.foot + self.slope * 3.0 * FOOT;
                    let start = self.edge(t, self.foot, low);
                    for (from, to) in [(start, self.foot), (self.foot, out)] {
                        let h = 0.5 * (to - from);
                        for (gs, ws) in &g {
                            let rho = from + h * (gs + 1.0);
                            ring += ws * self.column(rho, t) * rho * h;
                        }
                    }
                }
                sum += wa * ring * ha;
            }
        }
        sum
    }
}

trait SubSq {
    fn sub_sq(self, s: f64) -> f64;
}

impl SubSq for f64 {
    fn sub_sq(self, s: f64) -> f64 {
        s.mul_add(-s, self).max(0.0)
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

/// The largest angle between the normals of the two faces either side of
/// each edge bounding a face on a fitted surface, sampled along the edge's
/// own curve and projected onto each face's surface; with how many fitted
/// faces there are.
fn worst_crease(m: &Model, s: &Shape) -> (f64, usize) {
    let fitted = |g: &SurfaceGeometry| matches!(g, SurfaceGeometry::BSpline(_));
    let pieces = explore_unique(m, s, ShapeType::Face)
        .unwrap()
        .iter()
        .filter(|f| fitted(&surface_of(m, f)))
        .count();
    let mut worst = 0.0_f64;
    for e in explore_unique(m, s, ShapeType::Edge).unwrap() {
        let faces = ancestors_of(m, s, &e, ShapeType::Face).unwrap();
        let [a, b] = &faces[..] else {
            continue;
        };
        let (sa, sb) = (surface_of(m, a), surface_of(m, b));
        if !fitted(&sa) && !fitted(&sb) {
            continue;
        }
        let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap().unwrap();
        for k in 1..20 {
            let p = c.point_at(lo + (hi - lo) * f64::from(k) / 20.0, T).unwrap();
            let normal = |g: &SurfaceGeometry| -> Vector {
                let foot = ogeom_geom::project::project_on_surface(g, p, 16, T).unwrap();
                let (u, v) = foot.parameters;
                g.normal_at(u, v, T).unwrap().vector()
            };
            let (na, nb) = (normal(&sa), normal(&sb));
            worst = worst.max(na.dot(nb).abs().clamp(0.0, 1.0).acos());
        }
    }
    (worst, pieces)
}

/// A part whose boss's foot round blends into the rounds it reaches: the
/// solid, the foot rim, the reference for its fill, and how many rounds the
/// foot reaches.
struct Case {
    m: Model,
    s: Shape,
    foot: Vec<Shape>,
    reference: Closing,
    rounds: usize,
}

fn convex_edge_round() -> Round {
    Round {
        p: (0.0, 5.0),
        n: (0.0, -1.0),
        radius: 5.0,
        concave: false,
    }
}

/// A boss r 6 at `(12, 12)` by the corner where the block's two rounded
/// top edges meet: its foot round's contact circle (radius 9) crosses both
/// tangent lines, `x = 5` and `y = 5`, by 2, over two arcs apart (the
/// corner of the top stands 9.9 from the axis).
fn corner() -> Case {
    let mut m = Model::new();
    let rounded = rounded_block(&mut m, true);
    let boss = ogeom_algo::make_cylinder(&mut m, at(12.0, 12.0, 15.0), 6.0, 15.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(&mut m, &rounded, &boss, T).unwrap().shape;
    let foot = foot_rim(&m, &s, (12.0, 12.0, 15.0), 6.0);
    let side = Round {
        p: (5.0, 0.0),
        n: (-1.0, 0.0),
        ..convex_edge_round()
    };
    Case {
        m,
        s,
        foot,
        reference: Closing {
            axis: (12.0, 12.0),
            foot: 6.0,
            slope: 0.0,
            rounds: vec![convex_edge_round(), side],
            corner: None,
        },
        rounds: 2,
    }
}

/// The wall's lean of the conical boss: its radius falls from 6 to 4 over
/// its height of 15.
const LEAN: f64 = -2.0 / 15.0;

/// A conical boss at `(20, 12)` by the block's rounded top edge: its foot
/// round's contact circle (radius `6 + 3 sqrt(1 + k^2) + 3 k`) crosses the
/// tangent line `y = 5`.
fn conical() -> Case {
    let mut m = Model::new();
    let rounded = rounded_block(&mut m, false);
    let boss = ogeom_algo::make_cone(&mut m, at(20.0, 12.0, 15.0), 6.0, 4.0, 15.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(&mut m, &rounded, &boss, T).unwrap().shape;
    let foot = foot_rim(&m, &s, (20.0, 12.0, 15.0), 6.0);
    Case {
        m,
        s,
        foot,
        reference: Closing {
            axis: (20.0, 12.0),
            foot: 6.0,
            slope: LEAN,
            rounds: vec![convex_edge_round()],
            corner: None,
        },
        rounds: 1,
    }
}

/// A boss r 6 on the lower floor of a step, at `(20, 21)`; the step's
/// foot (the concave edge at `y = 8`) is rounded at 6, so its round meets
/// the floor along `y = 14`. The foot round's contact circle crosses that
/// line by 2, and the ball climbs the concave round, staying inside its
/// cylinder.
fn concave() -> Case {
    let mut m = Model::new();
    let base = ogeom_algo::make_box(&mut m, at(0.0, 0.0, 0.0), (40.0, 40.0, 10.0), T)
        .unwrap()
        .shape;
    let step = ogeom_algo::make_box(&mut m, at(0.0, 0.0, 10.0), (40.0, 8.0, 10.0), T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(&mut m, &base, &step, T).unwrap().shape;
    let valley = edges_where(&m, &s, |p| {
        (p.y - 8.0).abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9
    });
    let s = fillet_edges(&mut m, &s, &valley, 6.0, T).unwrap().shape;
    let boss = ogeom_algo::make_cylinder(&mut m, at(20.0, 21.0, 10.0), 6.0, 15.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(&mut m, &s, &boss, T).unwrap().shape;
    let foot = foot_rim(&m, &s, (20.0, 21.0, 10.0), 6.0);
    Case {
        m,
        s,
        foot,
        reference: Closing {
            axis: (20.0, 21.0),
            foot: 6.0,
            slope: 0.0,
            rounds: vec![Round {
                p: (0.0, 14.0),
                n: (0.0, -1.0),
                radius: 6.0,
                concave: true,
            }],
            corner: None,
        },
        rounds: 1,
    }
}

/// The case's foot rounded at `FOOT` by default: valid and closed, one
/// fitted face per round reached, meeting its neighbours G1. Blend into
/// round off falls to trim through, which refuses the foot by name as
/// before, and every behaviour off refuses it naming blend into round.
fn blends_smoothly(mut case: Case) {
    let (m, s, foot) = (&mut case.m, &case.s, &case.foot);
    let built = fillet_edges(m, s, foot, FOOT, T).unwrap();
    let found = ogeom_algo::check(m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    for shell in explore_unique(m, &built.shape, ShapeType::Shell).unwrap() {
        assert!(ogeom_algo::is_shell_closed(m, &shell).unwrap(), "closed");
    }
    for e in foot {
        assert!(!built.history.generated(e).is_empty(), "the rim's round");
    }
    let (crease, fitted) = worst_crease(m, &built.shape);
    println!("worst crease {crease:.2e} over {fitted} fitted faces");
    assert_eq!(fitted, case.rounds, "one fitted piece per round reached");
    assert!(crease < 1e-5, "the round's joins crease by {crease}");

    let mut off = Overflow::default();
    off.blend_into_round = false;
    let err = fillet_edges_overflow(m, s, foot, FOOT, off, T).unwrap_err();
    assert!(
        err.to_string().contains("does not reach so far"),
        "trim through's refusal: {err}"
    );
    let err = fillet_edges_overflow(m, s, foot, FOOT, Overflow::none(), T).unwrap_err();
    assert!(err.to_string().contains("blend into round"), "{err}");
}

/// The fill the case's foot round adds matches the reference within the
/// reference's bound plus 1e-6; the bound is the change from 32 to 64
/// points per range.
fn fills_as_the_closing(mut case: Case) {
    let (m, s, foot) = (&mut case.m, &case.s, &case.foot);
    let before = volume(m, s);
    let built = fillet_edges(m, s, foot, FOOT, T).unwrap();
    let (coarse, fine) = (case.reference.volume(32), case.reference.volume(64));
    let bound = (fine - coarse).abs().max(1e-9);
    let got = volume(m, &built.shape) - before;
    println!(
        "fill {got:.9}, reference {fine:.9}, error {:.2e}, bound {bound:.2e}",
        (got - fine).abs()
    );
    assert!(
        (got - fine).abs() < 1e-6 + bound,
        "the fill {got} against the reference {fine} (within {bound})"
    );
}

/// At the corner the round blends into both rounds; between them it is
/// the torus.
#[test]
fn a_foot_round_blends_into_two_rounds_at_a_corner() {
    blends_smoothly(corner());
}

#[test]
#[ignore = "heavy: the sampled reference takes about 10 s"]
fn a_foot_round_into_two_rounds_fills_as_the_closing() {
    fills_as_the_closing(corner());
}

#[test]
fn a_conical_boss_s_foot_round_blends_into_the_edge_round() {
    blends_smoothly(conical());
}

/// The reference is first checked against Pappus on a conical boss clear
/// of the round: the section between the wall, the floor and the ball's
/// arc is the quadrilateral of its foot, its touches and its centre, less
/// the ball's sector.
#[test]
#[ignore = "heavy: the sampled reference takes about 10 s"]
fn a_conical_boss_s_foot_round_fills_as_the_closing() {
    let k = LEAN;
    let clear = Closing {
        axis: (20.0, 40.0),
        foot: 6.0,
        slope: k,
        rounds: Vec::new(),
        corner: None,
    };
    let pappus = {
        let q = k.hypot(1.0);
        let big = FOOT.mul_add(k, FOOT.mul_add(q, 6.0));
        let wall = (big - FOOT / q, FOOT + FOOT * k / q);
        let quad = [(6.0, 0.0), (big, 0.0), (big, FOOT), wall];
        let (mut area, mut moment) = (0.0, 0.0);
        for i in 0..4 {
            let (a, b) = (quad[i], quad[(i + 1) % 4]);
            let cross = a.0 * b.1 - b.0 * a.1;
            area += 0.5 * cross;
            moment += (a.0 + b.0) * cross / 6.0;
        }
        if area < 0.0 {
            moment = -moment;
        }
        let (from, to) = (k.atan2(-1.0), -FRAC_PI_2);
        let sector =
            big * FOOT * FOOT * 0.5 * (to - from) + FOOT.powi(3) / 3.0 * (to.sin() - from.sin());
        TAU * (moment - sector)
    };
    let self_check = clear.volume(32);
    assert!(
        (self_check - pappus).abs() < 1e-6,
        "the reference {self_check} against Pappus {pappus}"
    );
    fills_as_the_closing(conical());
}

#[test]
fn a_foot_round_blends_into_a_concave_round() {
    blends_smoothly(concave());
}

#[test]
#[ignore = "heavy: the sampled reference takes about 10 s"]
fn a_foot_round_into_a_concave_round_fills_as_the_closing() {
    fills_as_the_closing(concave());
}

/// The part of the cone's place with its boss a cylinder r 6, cut away
/// beyond the plane `y = cut`: where `cut` is 12 the cut
/// holds the axis, the foot is a half circle and a straight edge, and the
/// half circle's ends stand on the flat. Its round reaches the edge round
/// over the middle of the arc.
fn cut_boss(cut: f64) -> (Model, Shape, Vec<Shape>) {
    let mut m = Model::new();
    let rounded = rounded_block(&mut m, false);
    let away = ogeom_algo::make_box(&mut m, at(0.0, cut, 15.0), (40.0, 30.0 - cut, 20.0), T)
        .unwrap()
        .shape;
    cut_boss_on(m, &rounded, &away)
}

/// The block `rounded` with a cylinder r 6, 15 tall, about `(20, 12)` on
/// its top, `away` cut from the cylinder before it is fused on, so the top
/// stays one face; with the foot's arc.
fn cut_boss_on(mut m: Model, rounded: &Shape, away: &Shape) -> (Model, Shape, Vec<Shape>) {
    let boss = ogeom_algo::make_cylinder(&mut m, at(20.0, 12.0, 15.0), 6.0, 15.0, T)
        .unwrap()
        .shape;
    let boss = ogeom_bool::cut(&mut m, &boss, away, T).unwrap().shape;
    let s = ogeom_bool::fuse(&mut m, rounded, &boss, T).unwrap().shape;
    let arc = foot_rim(&m, &s, (20.0, 12.0, 15.0), 6.0);
    (m, s, arc)
}

/// A half boss r 6 against a wall: the block's top behind `y = 12` rises
/// 15 as a wall, and a cylinder r 6, 10 tall, about `(20, 12)` stands
/// half out of it, so its foot is a half circle ending on the wall's
/// plane, which holds the axis.
fn boss_on_a_wall() -> (Model, Shape, Vec<Shape>) {
    let mut m = Model::new();
    let rounded = rounded_block(&mut m, false);
    let wall = ogeom_algo::make_box(&mut m, at(0.0, 12.0, 15.0), (40.0, 18.0, 15.0), T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(&mut m, &rounded, &wall, T).unwrap().shape;
    let boss = ogeom_algo::make_cylinder(&mut m, at(20.0, 12.0, 15.0), 6.0, 10.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(&mut m, &s, &boss, T).unwrap().shape;
    let arc = foot_rim(&m, &s, (20.0, 12.0, 15.0), 6.0);
    (m, s, arc)
}

/// The open arc of either part: the half circle in front of `y = 12`,
/// from angle pi to 2 pi about the axis, its round reaching the edge
/// round along `y = 5`.
fn open_arc(part: (Model, Shape, Vec<Shape>)) -> Case {
    let (m, s, foot) = part;
    Case {
        m,
        s,
        foot,
        reference: Closing {
            axis: (20.0, 12.0),
            foot: 6.0,
            slope: 0.0,
            rounds: vec![convex_edge_round()],
            corner: None,
        },
        rounds: 1,
    }
}

/// The open arc's round: valid and closed, one fitted face meeting its
/// neighbours G1. Blend into round off falls to trim through's result,
/// which differs, or its refusal by name; every behaviour off refuses it,
/// naming blend into round.
fn open_arc_blends(mut case: Case) {
    let (m, s, foot) = (&mut case.m, &case.s, &case.foot);
    let built = fillet_edges(m, s, foot, FOOT, T).unwrap();
    let found = ogeom_algo::check(m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    for shell in explore_unique(m, &built.shape, ShapeType::Shell).unwrap() {
        assert!(ogeom_algo::is_shell_closed(m, &shell).unwrap(), "closed");
    }
    for e in foot {
        assert!(!built.history.generated(e).is_empty(), "the rim's round");
    }
    let (crease, fitted) = worst_crease(m, &built.shape);
    println!("worst crease {crease:.2e} over {fitted} fitted faces");
    assert_eq!(fitted, case.rounds, "one fitted piece per round reached");
    assert!(crease < 1e-5, "the round's joins crease by {crease}");

    let mut off = Overflow::default();
    off.blend_into_round = false;
    match fillet_edges_overflow(m, s, foot, FOOT, off, T) {
        Ok(next) => assert!(
            (volume(m, &next.shape) - volume(m, &built.shape)).abs() > 1e-6,
            "trim through's own result"
        ),
        Err(err) => assert!(err.to_string().contains("reach"), "by name: {err}"),
    }
    let err = fillet_edges_overflow(m, s, foot, FOOT, Overflow::none(), T).unwrap_err();
    assert!(err.to_string().contains("blend into round"), "{err}");
}

/// The cut boss's open arc blends into the edge round, its round ending
/// flush with the flat.
#[test]
fn an_open_arc_of_a_cut_boss_blends_into_the_edge_round() {
    open_arc_blends(open_arc(cut_boss(12.0)));
}

/// The half boss's open arc against the wall blends into the edge round,
/// its round ending on the wall.
#[test]
fn an_open_arc_against_a_wall_blends_into_the_edge_round() {
    open_arc_blends(open_arc(boss_on_a_wall()));
}

/// The reference for an open arc is the closing of the whole boss within
/// the arc's wedge, pi to 2 pi: the end planes hold the axis and cut the
/// round. It is first checked against Pappus on a boss clear of the round:
/// half the turn of the section, the square of the radius less the
/// quarter disc, about the axis.
#[test]
#[ignore = "heavy: the sampled reference takes about 10 s"]
fn open_arcs_fill_as_the_closing_within_their_wedge() {
    let clear = Closing {
        axis: (20.0, 40.0),
        foot: 6.0,
        slope: 0.0,
        rounds: Vec::new(),
        corner: None,
    };
    let r = FOOT;
    let area = r * r * (1.0 - PI / 4.0);
    // The section's centroid from the wall: `r (10 - 3 pi) / (12 - 3 pi)`.
    let from_wall = r * 3.0f64.mul_add(-PI, 10.0) / 3.0f64.mul_add(-PI, 12.0);
    let pappus = PI * area * (6.0 + from_wall);
    let self_check = clear.volume_within(32, Some((PI, TAU)));
    assert!(
        (self_check - pappus).abs() < 1e-6,
        "the reference {self_check} against Pappus {pappus}"
    );
    for part in [cut_boss(12.0), boss_on_a_wall()] {
        let mut case = open_arc(part);
        let (m, s, foot) = (&mut case.m, &case.s, &case.foot);
        let before = volume(m, s);
        let built = fillet_edges(m, s, foot, FOOT, T).unwrap();
        let wedge = Some((PI, TAU));
        let (coarse, fine) = (
            case.reference.volume_within(32, wedge),
            case.reference.volume_within(64, wedge),
        );
        let bound = (fine - coarse).abs().max(1e-9);
        let got = volume(m, &built.shape) - before;
        println!(
            "fill {got:.9}, reference {fine:.9}, error {:.2e}, bound {bound:.2e}",
            (got - fine).abs()
        );
        assert!(
            (got - fine).abs() < 1e-6 + bound,
            "the fill {got} against the reference {fine} (within {bound})"
        );
    }
}

/// What blend into round does not take of an open arc gives what blend
/// into round off gives.
fn not_blended(mut m: Model, s: &Shape, arc: &[Shape]) {
    let mut off = Overflow::default();
    off.blend_into_round = false;
    let default = fillet_edges(&mut m, s, arc, FOOT, T);
    let trimmed = fillet_edges_overflow(&mut m, s, arc, FOOT, off, T);
    match (default, trimmed) {
        (Err(a), Err(b)) => assert_eq!(a.to_string(), b.to_string()),
        (Ok(a), Ok(b)) => {
            assert!((volume(&m, &a.shape) - volume(&m, &b.shape)).abs() < 1e-9);
        }
        (a, b) => panic!("the default {a:?} against trim through {b:?}"),
    }
}

/// A boss cut by a plane beside its axis, `y = 13`: the arc's ends stand
/// on a plane the round's sections do not lie in, so its end would need
/// the plane to trim the round obliquely. Not blended.
#[test]
fn an_open_arc_ending_off_the_axis_is_not_blended() {
    let (m, s, arc) = cut_boss(13.0);
    not_blended(m, &s, &arc);
}

/// A boss cut by `x = 20` from the floor up: the arc runs from angle
/// pi / 2 to 3 pi / 2, and the round's seat on the edge round straddles
/// its end at 3 pi / 2. Not blended.
#[test]
fn an_open_arc_ending_on_its_round_s_seat_is_not_blended() {
    let mut m = Model::new();
    let rounded = rounded_block(&mut m, false);
    let away = ogeom_algo::make_box(&mut m, at(20.0, 0.0, 15.0), (20.0, 30.0, 20.0), T)
        .unwrap()
        .shape;
    let (m, s, arc) = cut_boss_on(m, &rounded, &away);
    not_blended(m, &s, &arc);
}

/// The block of [`rounded_block`] whose top edges along `y = 0` and
/// `x = 0` meet round a rounded corner: with `vertical` its upright edge
/// at the origin is rounded at that radius first, so the top's two rounds
/// r 5 run on round the corner as one chain, a torus between them;
/// without, the three edges at the corner are rounded together at 5 and
/// a sphere rounds the corner.
fn corner_block(m: &mut Model, vertical: Option<f64>) -> Shape {
    let block = ogeom_algo::make_box(m, at(0.0, 0.0, 0.0), (40.0, 30.0, 15.0), T)
        .unwrap()
        .shape;
    let upright = |p: Point| p.x.abs() < 1e-9 && p.y.abs() < 1e-9;
    let Some(rho) = vertical else {
        let edges = edges_where(m, &block, |p| {
            ((p.y.abs() < 1e-9 || p.x.abs() < 1e-9) && (p.z - 15.0).abs() < 1e-9) || upright(p)
        });
        return fillet_edges(m, &block, &edges, 5.0, T).unwrap().shape;
    };
    let edge = edges_where(m, &block, upright);
    let s = fillet_edges(m, &block, &edge, rho, T).unwrap().shape;
    let top = edges_where(m, &s, |p| {
        (p.z - 15.0).abs() < 1e-9
            && (p.x.abs() < 1e-9
                || p.y.abs() < 1e-9
                || ((p.x - rho).hypot(p.y - rho) - rho).abs() < 1e-9)
    });
    fillet_edges(m, &s, &top, 5.0, T).unwrap().shape
}

/// A boss r `r` at `(x, y)` on [`corner_block`], its foot round's contact
/// circle crossing both tangent lines and the corner between them; the
/// corner is rounded at `rho` seen from above (nought for the sphere).
fn crossing(vertical: Option<f64>, (x, y): (f64, f64), r: f64) -> Case {
    let mut m = Model::new();
    let block = corner_block(&mut m, vertical);
    let boss = ogeom_algo::make_cylinder(&mut m, at(x, y, 15.0), r, 15.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(&mut m, &block, &boss, T).unwrap().shape;
    let foot = foot_rim(&m, &s, (x, y, 15.0), r);
    let side = Round {
        p: (5.0, 0.0),
        n: (-1.0, 0.0),
        ..convex_edge_round()
    };
    Case {
        m,
        s,
        foot,
        reference: Closing {
            axis: (x, y),
            foot: r,
            slope: 0.0,
            rounds: vec![convex_edge_round(), side],
            corner: Some(vertical.map_or(0.0, |v| v - 5.0)),
        },
        rounds: 3,
    }
}

/// The upright edge rounded at 12, so the top's outline turns round an arc
/// r 7 about `(12, 12)` and a torus joins the two rounds. A boss r 5 at
/// `(11, 11)`: its contact circle (radius 8) crosses `y = 5` at
/// `x = 16.3`, passes 2.4 beyond the arc and comes back across `x = 5`, so
/// the ball rests on the first round, the torus and the second in turn.
fn torus_corner() -> Case {
    crossing(Some(12.0), (11.0, 11.0), 5.0)
}

/// The three edges at the corner rounded together at 5, a sphere about
/// `(5, 5, 10)` rounding it. A boss r 4 at `(9.5, 10)`: its contact circle
/// (radius 7) crosses both tangent lines and passes 0.27 beyond the top's
/// corner at `(5, 5)`, where the ball rests on the sphere.
fn sphere_corner() -> Case {
    crossing(None, (9.5, 10.0), 4.0)
}

/// Where the two rounds meet round a torus, the ball passes from the
/// first round to the torus and on to the second: three fitted pieces
/// meeting each other and the torus pieces G1.
#[test]
fn a_foot_round_blends_across_a_torus_corner() {
    blends_smoothly(torus_corner());
}

/// The same across a sphere rounding the corner.
#[test]
fn a_foot_round_blends_across_a_sphere_corner() {
    blends_smoothly(sphere_corner());
}

#[test]
#[ignore = "heavy: the sampled reference takes about 20 s"]
fn a_foot_round_across_a_torus_corner_fills_as_the_closing() {
    fills_as_the_closing(torus_corner());
}

#[test]
#[ignore = "heavy: the sampled reference takes about 20 s"]
fn a_foot_round_across_a_sphere_corner_fills_as_the_closing() {
    fills_as_the_closing(sphere_corner());
}

/// The two top edges rounded alone meet in a crease, the upright edge left
/// sharp: there the ball would rest on the crease between the rounds, an
/// edge seat blend into round does not build. The boss of
/// [`sphere_corner`] there gets what blend into round off gives.
#[test]
fn a_foot_round_reaching_a_creased_corner_is_not_blended() {
    let mut m = Model::new();
    let block = rounded_block(&mut m, true);
    let boss = ogeom_algo::make_cylinder(&mut m, at(9.5, 10.0, 15.0), 4.0, 15.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(&mut m, &block, &boss, T).unwrap().shape;
    let foot = foot_rim(&m, &s, (9.5, 10.0, 15.0), 4.0);
    not_blended(m, &s, &foot);
}
