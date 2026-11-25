//! Hold edge (`Overflow::hold_edge`, `Overflow::hold_edge_at_end`): a
//! straight round whose ball, rolled on its two faces, would roll over a
//! sharp edge of one of them. The edge is kept: the ball there touches the
//! other face and rests on the edge, sweeping a canal surface that meets
//! the band G1 where the contact crosses the edge.
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom_core::Tolerances;
use ogeom_fillet::{Overflow, fillet_edges, fillet_edges_overflow};
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The round's radius on both parts.
const R: f64 = 4.0;

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn edges_where(m: &Model, s: &Shape, on: impl Fn(Point) -> bool) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
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

/// The solid lofted between two polygons in the planes `x = x0` and
/// `x = x1`, corner `k` of one joined to corner `k` of the other, every
/// side planar.
fn loft(m: &mut Model, x0: f64, a: &[(f64, f64)], x1: f64, b: &[(f64, f64)]) -> Shape {
    let n = a.len();
    let mut points: Vec<Point> = a.iter().map(|(y, z)| Point::new(x0, *y, *z)).collect();
    points.extend(b.iter().map(|(y, z)| Point::new(x1, *y, *z)));
    for flip in [false, true] {
        let mut ends = [(0..n).collect::<Vec<_>>(), (n..2 * n).rev().collect()];
        let mut sides: Vec<Vec<usize>> = (0..n)
            .map(|i| vec![i, n + i, n + (i + 1) % n, (i + 1) % n])
            .collect();
        if flip {
            ends.iter_mut().for_each(|r| r.reverse());
            sides.iter_mut().for_each(|r| r.reverse());
        }
        let rings: Vec<Vec<usize>> = ends.into_iter().chain(sides).collect();
        if let Ok((built, _)) = ogeom_algo::make_faceted_solid(m, &points, &rings, &[], T) {
            return built.shape;
        }
    }
    panic!("the loft's rings enclose no volume either way round");
}

/// The sharp edge a part's top face ends at, beside the rounded edge: a
/// straight cliff edge through `(x, y)` along `(1, slope)`, or a pocket's
/// rim, the circle of `radius` about `(x, y)`.
#[derive(Clone, Copy)]
enum Cliff {
    Line { y0: f64, slope: f64 },
    Rim { x: f64, y: f64, radius: f64 },
}

impl Cliff {
    /// Whether the point `(x, y)` of the top's plane lies past the edge,
    /// off the top.
    fn past(self, x: f64, y: f64) -> bool {
        match self {
            Self::Line { y0, slope } => y > slope.mul_add(x, y0),
            Self::Rim {
                x: cx,
                y: cy,
                radius,
            } => (x - cx).hypot(y - cy) < radius,
        }
    }

    /// The distance in the top's plane from `(x, y)` to the edge's curve.
    fn distance(self, x: f64, y: f64) -> f64 {
        match self {
            Self::Line { y0, slope } => (y - slope.mul_add(x, y0)).abs() / slope.hypot(1.0),
            Self::Rim {
                x: cx,
                y: cy,
                radius,
            } => ((x - cx).hypot(y - cy) - radius).abs(),
        }
    }

    /// Whether `p` lies on the edge's curve, in the plane `z = 0`.
    fn holds(self, p: Point) -> bool {
        p.z.abs() < 1e-6 && self.distance(p.x, p.y) < 1e-6
    }
}

/// The independent reference for a round of radius `R` on the edge along
/// `x` from 0 to `len` between the top `z = 0` (material below) and the
/// front `y = 0` (material behind), the top ending at `cliff`. It is the
/// definition, sampled: the ball's centre stands `R` behind the front, at
/// `R` under the top where its contact `(x, R, 0)` lies on the top, else
/// where it is `R` from the cliff's curve (closed form, the curve lying in
/// the top's plane). A point of the corner is cut away where it lies
/// farther than `R` from every centre, on the corner's side of its nearest
/// one (behind and above it). Past the part's ends the centres run on as
/// they are, the end faces cutting the round.
struct Reference {
    cliff: Cliff,
    len: f64,
}

impl Reference {
    fn centre(&self, t: f64) -> Point {
        if !self.cliff.past(t, R) {
            return Point::new(t, R, -R);
        }
        let h = self.cliff.distance(t, R);
        Point::new(t, R, -(R * R - h * h).max(0.0).sqrt())
    }

    fn cut(&self, p: Point) -> bool {
        if p.z > 0.0 || p.y < 0.0 {
            return false;
        }
        let d = |t: f64| p.distance(self.centre(t));
        let (mut lo, mut hi) = (p.x - 2.0 * R, p.x + 2.0 * R);
        let mut best = p.x;
        let mut least = f64::INFINITY;
        for k in 0..=64 {
            let t = lo + (hi - lo) * f64::from(k) / 64.0;
            if d(t) < least {
                least = d(t);
                best = t;
            }
        }
        let step = (hi - lo) / 64.0;
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
        let c = self.centre(0.5 * (lo + hi));
        let w = p - c;
        w.magnitude() > R && w.y <= 0.0 && w.z >= 0.0
    }

    /// How far behind the front the cut reaches at `(x, z)`, by bisection.
    fn depth(&self, x: f64, z: f64) -> f64 {
        let (mut lo, mut hi) = (0.0, R);
        if !self.cut(Point::new(x, 1e-12, z)) {
            return 0.0;
        }
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            if self.cut(Point::new(x, mid, z)) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// The volume cut away: in `x` over `ranges`, each Gauss-Legendre of
    /// `n` points; down the front from the top to the ball's contact there,
    /// in `s` with `z = -s^2` (smooth where the band meets the top
    /// tangentially), each column's depth.
    fn volume(&self, ranges: &[(f64, f64)], n: usize) -> f64 {
        let g = gauss(n);
        let mut sum = 0.0;
        for &(x0, x1) in ranges {
            let hx = 0.5 * (x1 - x0);
            for (gx, wx) in &g {
                let x = x0 + hx * (gx + 1.0);
                let s1 = (-self.centre(x).z).sqrt();
                let hs = 0.5 * s1;
                for (gs, ws) in &g {
                    let s = hs * (gs + 1.0);
                    sum += wx * ws * hx * hs * self.depth(x, -s * s) * 2.0 * s;
                }
            }
        }
        debug_assert!(self.len > 0.0);
        sum
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

/// The largest angle between the normals either side of each edge where a
/// fitted face meets a cylinder or a plane tangent to it (the band, the
/// front), sampled along the edge.
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
        let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap().unwrap();
        let mid = c.point_at(0.5 * (lo + hi), T).unwrap();
        // The front `y = 0` or a band about a line along `x`; not the top,
        // the pocket's wall or the cliff, where the held edge is sharp.
        let tangent_kind = match other {
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

/// The x-range the result's edges lying on the cliff's curve cover, and
/// how many there are.
fn held_span(m: &Model, s: &Shape, cliff: Cliff) -> (f64, f64, usize) {
    let held = edges_where(m, s, |p| cliff.holds(p));
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for e in &held {
        let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
        for k in 0..=50 {
            let p = c.point_at(a + (b - a) * f64::from(k) / 50.0, T).unwrap();
            if p.y < R + 1e-6 {
                lo = lo.min(p.x);
                hi = hi.max(p.x);
            }
        }
    }
    (lo, hi, held.len())
}

/// A block 30 by 20 by 15 under the plane `z = 0`, a pocket of radius 5 and
/// depth 8 sunk into its top about `(15, 7)`, its front top edge along
/// `y = 0` asked.
fn pocketed(m: &mut Model) -> (Shape, Vec<Shape>) {
    let block = ogeom_algo::make_box(m, at(0.0, 0.0, -15.0), (30.0, 20.0, 15.0), T)
        .unwrap()
        .shape;
    let pocket = ogeom_algo::make_cylinder(m, at(15.0, 7.0, -8.0), 5.0, 9.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(m, &block, &pocket, T).unwrap().shape;
    let edge = edges_where(m, &s, |p| p.y.abs() < 1e-9 && p.z.abs() < 1e-9);
    assert_eq!(edge.len(), 1, "the front top edge");
    (s, edge)
}

const RIM: Cliff = Cliff::Rim {
    x: 15.0,
    y: 7.0,
    radius: 5.0,
};

/// A block 20 long whose top `z = 0` narrows from 6 deep at `x = 0` to 2 at
/// `x = 20`, a cliff 6 high dropping from its back edge to a ledge, the
/// block 12 deep and 15 tall below; its front top edge along `y = 0`
/// asked.
fn cliffed(m: &mut Model) -> (Shape, Vec<Shape>) {
    let end0 = [
        (0.0, 0.0),
        (0.0, -15.0),
        (12.0, -15.0),
        (12.0, -6.0),
        (6.0, -6.0),
        (6.0, 0.0),
    ];
    let end1 = [
        (0.0, 0.0),
        (0.0, -15.0),
        (12.0, -15.0),
        (12.0, -6.0),
        (2.0, -6.0),
        (2.0, 0.0),
    ];
    let s = loft(m, 0.0, &end0, 20.0, &end1);
    let edge = edges_where(m, &s, |p| p.y.abs() < 1e-9 && p.z.abs() < 1e-9);
    assert_eq!(edge.len(), 1, "the front top edge");
    (s, edge)
}

const SLOPE: Cliff = Cliff::Line {
    y0: 6.0,
    slope: -0.2,
};

/// The checks both parts share: valid and closed, the cliff's edge kept on
/// its curve over `held`, the volume cut away matching the reference split
/// at `ranges` within its own bound (the change from 24 to 48 points per
/// range) plus 1e-6, and the canal meeting the band and the front G1.
fn holds(
    m: &Model,
    before: f64,
    built: &ogeom_algo::Built,
    edge: &Shape,
    reference: &Reference,
    held: (f64, f64),
    ranges: &[(f64, f64)],
) {
    let (cliff, len) = (reference.cliff, reference.len);
    let found = ogeom_algo::check(m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    for shell in explore_unique(m, &built.shape, ShapeType::Shell).unwrap() {
        assert!(ogeom_algo::is_shell_closed(m, &shell).unwrap(), "closed");
    }
    assert!(!built.history.generated(edge).is_empty(), "the round");

    let (lo, hi, count) = held_span(m, &built.shape, cliff);
    assert!(count > 0, "the held edge is kept");
    assert!(
        lo <= held.0 + 1e-6 && hi >= held.1 - 1e-6,
        "the held edge spans {lo}..{hi}, not {held:?}"
    );

    let (coarse, fine) = (reference.volume(ranges, 24), reference.volume(ranges, 48));
    let bound = (fine - coarse).abs();
    let got = before - volume(m, &built.shape);
    assert!(
        (got - fine).abs() < 1e-6 + bound,
        "cut away {got} against the reference {fine} (within {bound})"
    );
    // Not the band trimmed through: the ball rests on the held edge.
    let band = len * R * R * (1.0 - PI / 4.0);
    assert!((got - band).abs() > 1e-2, "{got} against the band's {band}");

    let (crease, joins) = worst_join(m, &built.shape);
    assert!(joins >= 2, "the canal's joins to the band and the front");
    assert!(crease < 1e-5, "the canal's joins crease by {crease}");
}

/// The reference reproduces the band's closed form, `len R^2 (1 - pi/4)`,
/// where the top reaches past the contact everywhere.
#[test]
fn the_reference_reproduces_the_band() {
    let clear = Reference {
        cliff: Cliff::Line {
            y0: 9.0,
            slope: 0.0,
        },
        len: 10.0,
    };
    let want = 10.0 * R * R * (1.0 - PI / 4.0);
    let got = clear.volume(&[(0.0, 10.0)], 24);
    assert!((got - want).abs() < 1e-9, "{got} against {want}");
}

/// Hold edge: the front top edge of a block rounded at 4 beside a pocket
/// whose rim comes within 2 of it. Between `x = 11` and `x = 19` the
/// band's contact `y = 4` falls inside the rim: the ball rests on the rim
/// there and the rim is kept.
#[test]
fn hold_edge_keeps_a_pocket_rim_the_ball_would_roll_over() {
    let mut m = Model::new();
    let (s, edge) = pocketed(&mut m);
    let before = volume(&m, &s);
    let built = fillet_edges(&mut m, &s, &edge, R, T).unwrap();
    holds(
        &m,
        before,
        &built,
        &edge[0],
        &Reference {
            cliff: RIM,
            len: 30.0,
        },
        (11.0, 19.0),
        &[(0.0, 11.0), (11.0, 19.0), (19.0, 30.0)],
    );
}

/// Hold edge at the end: the front top edge of a block whose top narrows
/// below the radius past `x = 10`. From there to the end the ball rests on
/// the cliff edge, which is kept, and the end face cuts the canal.
#[test]
fn hold_edge_at_the_end_keeps_a_cliff_edge_the_band_s_end_runs_over() {
    let mut m = Model::new();
    let (s, edge) = cliffed(&mut m);
    let before = volume(&m, &s);
    let built = fillet_edges(&mut m, &s, &edge, R, T).unwrap();
    holds(
        &m,
        before,
        &built,
        &edge[0],
        &Reference {
            cliff: SLOPE,
            len: 20.0,
        },
        (10.0, 20.0),
        &[(0.0, 10.0), (10.0, 20.0)],
    );
}

fn with(f: impl FnOnce(&mut Overflow)) -> Overflow {
    let mut o = Overflow::default();
    f(&mut o);
    o
}

/// Whether a result edge on the cliff's curve passes nearer the rounded
/// edge than `y`: the held edge kept where the ball rests on it.
fn holds_near(m: &Model, s: &Shape, cliff: Cliff, y: f64) -> bool {
    edges_where(m, s, |p| cliff.holds(p)).iter().any(|e| {
        let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
        (0..=20).any(|k| c.point_at(a + (b - a) * f64::from(k) / 20.0, T).unwrap().y < y)
    })
}

/// Each option switched off: hold edge off gives trim through's notch (the
/// band runs on through the pocket's wall and cuts the rim away), and with
/// trim through off too the request is refused naming hold edge. Hold edge
/// at the end does not bear on a run inside the band, and trim through off
/// alone leaves hold edge to build.
#[test]
fn hold_edge_off_falls_to_trim_through_or_is_refused() {
    let mut m = Model::new();
    let (s, edge) = pocketed(&mut m);
    let off =
        fillet_edges_overflow(&mut m, &s, &edge, R, with(|o| o.hold_edge = false), T).unwrap();
    assert!(ogeom_algo::check(&m, &off.shape, T).unwrap().is_valid());
    assert!(
        !holds_near(&m, &off.shape, RIM, 3.0),
        "trim through cuts the rim away where the band crosses it"
    );
    let both = with(|o| {
        o.hold_edge_at_end = false;
        o.trim_through = false;
    });
    let kept = fillet_edges_overflow(&mut m, &s, &edge, R, both, T).unwrap();
    assert!(
        holds_near(&m, &kept.shape, RIM, 3.0),
        "hold edge keeps the rim"
    );

    let both = with(|o| {
        o.hold_edge = false;
        o.trim_through = false;
    });
    let err = fillet_edges_overflow(&mut m, &s, &edge, R, both, T).unwrap_err();
    assert!(err.to_string().contains("hold edge is off"), "{err}");
    let err = fillet_edges_overflow(&mut m, &s, &edge, R, Overflow::none(), T).unwrap_err();
    assert!(err.to_string().contains("hold edge is off"), "{err}");
}

/// At the end: hold edge at the end off, or hold edge off, gives trim
/// through's result, and with trim through off too it is refused naming
/// hold edge at the end; trim through off alone leaves hold edge at the end
/// to build.
#[test]
fn hold_edge_at_end_off_falls_to_trim_through_or_is_refused() {
    let mut m = Model::new();
    let (s, edge) = cliffed(&mut m);
    let end_off = fillet_edges_overflow(
        &mut m,
        &s,
        &edge,
        R,
        with(|o| o.hold_edge_at_end = false),
        T,
    )
    .unwrap();
    let off =
        fillet_edges_overflow(&mut m, &s, &edge, R, with(|o| o.hold_edge = false), T).unwrap();
    for trimmed in [&end_off, &off] {
        assert!(ogeom_algo::check(&m, &trimmed.shape, T).unwrap().is_valid());
        assert!(
            !holds_near(&m, &trimmed.shape, SLOPE, 3.0),
            "trim through cuts the cliff edge away past the band's contact"
        );
    }
    assert!((volume(&m, &end_off.shape) - volume(&m, &off.shape)).abs() < 1e-9);

    let both = with(|o| {
        o.hold_edge_at_end = false;
        o.trim_through = false;
    });
    let err = fillet_edges_overflow(&mut m, &s, &edge, R, both, T).unwrap_err();
    assert!(
        err.to_string().contains("hold edge at the end is off"),
        "{err}"
    );
    let trim_off =
        fillet_edges_overflow(&mut m, &s, &edge, R, with(|o| o.trim_through = false), T).unwrap();
    assert!(
        holds_near(&m, &trim_off.shape, SLOPE, 3.0),
        "hold edge at the end keeps the cliff edge"
    );
}
