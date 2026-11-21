//! Hold edge where the ball's seat changes more than once: it rests on a
//! held edge, on the vertex where that edge ends, on the next held edge,
//! or on two pockets' rims in turn. The reference is the definition: the
//! ball touches the front and rests on the nearest point of the top's
//! region.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom_core::Tolerances;
use ogeom_fillet::{Overflow, fillet_edges, fillet_edges_overflow};
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The round's radius on every part.
const R: f64 = 4.0;

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn edges_where(m: &Model, s: &Shape, on: impl Fn(Point) -> bool) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            // A degenerate edge (where a piece closes to a vertex) has no
            // curve.
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

/// The faceted solid over `points` and `rings`, the rings taken as given
/// or every one reversed, whichever encloses a volume.
fn faceted(m: &mut Model, points: &[Point], rings: &[Vec<usize>]) -> Shape {
    for flip in [false, true] {
        let mut rings = rings.to_vec();
        if flip {
            rings.iter_mut().for_each(|r| r.reverse());
        }
        if let Ok((built, _)) = ogeom_algo::make_faceted_solid(m, points, &rings, &[], T) {
            return built.shape;
        }
    }
    panic!("the rings enclose no volume either way round");
}

/// The top `z = 0` of a part, the front `y = 0` below it: where it ends,
/// and the distance in its plane from a point to it.
#[derive(Clone, Copy)]
enum Top {
    /// The top lies in front of the polyline through these points, its
    /// first and last pieces carried on past their ends.
    Below(&'static [(f64, f64)]),
    /// The top holds these pockets, each `(x, y, radius)`.
    Pockets(&'static [(f64, f64, f64)]),
}

/// The distance from `(x, y)` to the segment `a b`, or to the line through
/// it carried on past `a` or `b` where `open` says so.
fn to_piece(x: f64, y: f64, a: (f64, f64), b: (f64, f64), open: (bool, bool)) -> f64 {
    let (dx, dy) = (b.0 - a.0, b.1 - a.1);
    let mut t = ((x - a.0) * dx + (y - a.1) * dy) / dx.mul_add(dx, dy * dy);
    if !open.0 {
        t = t.max(0.0);
    }
    if !open.1 {
        t = t.min(1.0);
    }
    (x - t.mul_add(dx, a.0)).hypot(y - t.mul_add(dy, a.1))
}

impl Top {
    /// The distance in the top's plane from `(x, y)` to the top, nought on
    /// it.
    fn distance(self, x: f64, y: f64) -> f64 {
        match self {
            Self::Below(line) => {
                let n = line.len();
                let k = (1..n - 1).find(|k| x < line[*k].0).unwrap_or(n - 1);
                let (a, b) = (line[k - 1], line[k]);
                let over = (b.1 - a.1).mul_add((x - a.0) / (b.0 - a.0), a.1);
                if y <= over {
                    return 0.0;
                }
                (1..n)
                    .map(|k| to_piece(x, y, line[k - 1], line[k], (k == 1, k == n - 1)))
                    .fold(f64::INFINITY, f64::min)
            }
            Self::Pockets(pockets) => pockets
                .iter()
                .map(|(cx, cy, r)| (r - (x - cx).hypot(y - cy)).max(0.0))
                .fold(0.0, f64::max),
        }
    }
}

/// The independent reference for a round of radius `R` on the edge along
/// `x` from 0 to `len` between the top `z = 0` (material below) and the
/// front `y = 0` (material behind). The ball's centre stands `R` behind the
/// front, as high as it may while it stays `R` from the top's region: at
/// depth `sqrt(R^2 - h^2)`, `h` the distance from its contact `(x, R)` to
/// the top in the top's plane. A point of the corner is cut away where it
/// lies farther than `R` from every centre, on the corner's side of its
/// nearest one. Past the part's ends the centres run on as they are.
struct Reference {
    top: Top,
}

impl Reference {
    fn centre(&self, t: f64) -> Point {
        let h = self.top.distance(t, R).min(R);
        Point::new(t, R, -(R * R - h * h).sqrt())
    }

    fn cut(&self, p: Point) -> bool {
        if p.z > 0.0 || p.y < 0.0 {
            return false;
        }
        let d = |t: f64| p.distance(self.centre(t));
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
    /// in `s` with `z = -s^2`, each column's depth.
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

/// The largest angle between the normals either side of each edge a fitted
/// face shares with another fitted face, a cylinder along `x` (the band) or
/// the front `y = 0`, sampled along the edge, and how many such edges.
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

/// The x-range covered by the result's edges lying on the top's plane
/// along the line through `a` and `b`, behind the front by less than `R`.
fn kept_along(m: &Model, s: &Shape, a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    let on = |p: Point| p.z.abs() < 1e-6 && to_piece(p.x, p.y, a, b, (true, true)) < 1e-6;
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for e in edges_where(m, s, on) {
        let (c, (ta, tb)) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap().unwrap();
        for k in 0..=50 {
            let p = c.point_at(ta + (tb - ta) * f64::from(k) / 50.0, T).unwrap();
            if p.y < R + 1e-6 {
                lo = lo.min(p.x);
                hi = hi.max(p.x);
            }
        }
    }
    (lo, hi)
}

/// The checks every part shares: valid and closed, and every join of a
/// fitted piece tangent within 1e-5 rad; with `heavy`, the volume cut away
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
    len: f64,
) {
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
    if !heavy {
        return;
    }
    let (coarse, fine) = (reference.volume(ranges, 24), reference.volume(ranges, 48));
    let bound = (fine - coarse).abs();
    let got = before - volume(m, &built.shape);
    eprintln!("cut away {got:.9}, reference {fine:.9}, bound {bound:.1e}");
    assert!(
        (got - fine).abs() < 1e-6 + bound,
        "cut away {got} against the reference {fine} (within {bound})"
    );
    let band = len * R * R * (1.0 - PI / 4.0);
    assert!((got - band).abs() > 1e-2, "{got} against the band's {band}");
}

/// The top's back edge of the kinked part: down from 4.6 deep at `x = 0`
/// to 2.5 at `x = 14`, then more steeply to 0.4 at `x = 20`.
const KINK: [(f64, f64); 3] = [(0.0, 4.6), (14.0, 2.5), (20.0, 0.4)];

/// A part 20 long whose top `z = 0` lies in front of [`KINK`], a cliff 6
/// high falling from that edge to a ledge, the block 12 deep and 15 tall
/// below; its front top edge along `y = 0` asked.
fn kinked(m: &mut Model) -> (Shape, Vec<Shape>) {
    let p = |x: f64, y: f64, z: f64| Point::new(x, y, z);
    let points = [
        p(0.0, 0.0, 0.0),    // 0
        p(20.0, 0.0, 0.0),   // 1
        p(20.0, 0.4, 0.0),   // 2
        p(14.0, 2.5, 0.0),   // 3
        p(0.0, 4.6, 0.0),    // 4
        p(20.0, 0.4, -6.0),  // 5
        p(14.0, 2.5, -6.0),  // 6
        p(0.0, 4.6, -6.0),   // 7
        p(20.0, 12.0, -6.0), // 8
        p(0.0, 12.0, -6.0),  // 9
        p(0.0, 0.0, -15.0),  // 10
        p(20.0, 0.0, -15.0), // 11
        p(20.0, 12.0, -15.0),
        p(0.0, 12.0, -15.0),
    ];
    let rings = vec![
        vec![0, 1, 2, 3, 4],
        vec![10, 13, 12, 11],
        vec![0, 10, 11, 1],
        vec![1, 11, 12, 8, 5, 2],
        vec![0, 4, 7, 9, 13, 10],
        vec![3, 6, 7, 4],
        vec![2, 5, 6, 3],
        vec![5, 8, 9, 7, 6],
        vec![8, 12, 13, 9],
    ];
    let s = faceted(m, &points, &rings);
    let edge = edges_where(m, &s, |p| p.y.abs() < 1e-9 && p.z.abs() < 1e-9);
    assert_eq!(edge.len(), 1, "the front top edge");
    (s, edge)
}

/// The ball rests on a vertex: the front top edge of [`kinked`] rounded at
/// 4. Past `x = 4` its contact `y = 4` lies behind the top, and the ball
/// rests on the first cliff edge; between `x = 14.225` and `14.525` its
/// nearest point of the top is the kink's vertex `(14, 2.5)`, so its centre
/// stays 4 from that vertex; past it, the ball rests on the second cliff
/// edge to the end. Both cliff edges are kept, and the three fitted pieces
/// meet each other and the band G1.
fn hold_edge_rolls_the_ball_round_the_vertex_where_a_held_edge_ends_case(heavy: bool) {
    let mut m = Model::new();
    let (s, edge) = kinked(&mut m);
    let before = volume(&m, &s);
    let built = fillet_edges(&mut m, &s, &edge, R, T).unwrap();
    holds(
        heavy,
        &m,
        before,
        &built,
        &edge[0],
        &Reference {
            top: Top::Below(&KINK),
        },
        &[
            (0.0, 4.0),
            (4.0, 12.0),
            (12.0, 13.5),
            (13.5, 14.225),
            (14.225, 14.525),
            (14.525, 15.25),
            (15.25, 16.5),
            (16.5, 20.0),
        ],
        20.0,
    );
    let (lo, hi) = kept_along(&m, &built.shape, KINK[0], KINK[1]);
    assert!(
        lo <= 4.0 + 1e-6 && hi >= 14.0 - 1e-6,
        "first edge {lo}..{hi}"
    );
    let (lo, hi) = kept_along(&m, &built.shape, KINK[1], KINK[2]);
    assert!(
        lo <= 14.0 + 1e-6 && hi >= 20.0 - 1e-6,
        "second edge {lo}..{hi}"
    );
}

#[test]
fn hold_edge_rolls_the_ball_round_the_vertex_where_a_held_edge_ends() {
    hold_edge_rolls_the_ball_round_the_vertex_where_a_held_edge_ends_case(false);
}

/// The vertex seat's volume against the reference, split at its seat
/// changes and finer round them.
#[test]
#[ignore = "heavy"]
fn a_vertex_seat_cuts_away_the_reference() {
    hold_edge_rolls_the_ball_round_the_vertex_where_a_held_edge_ends_case(true);
}

/// The run reaches the part's end: with hold edge at the end off, or hold
/// edge off, the request goes on to trim through, which refuses this part
/// by name (its band would set back across a top too shallow for it), and
/// with trim through off too it is refused naming hold edge at the end.
#[test]
fn a_vertex_seat_at_the_end_follows_hold_edge_at_end() {
    let mut m = Model::new();
    let (s, edge) = kinked(&mut m);
    let mut end_off = Overflow::default();
    end_off.hold_edge_at_end = false;
    let mut off = Overflow::default();
    off.hold_edge = false;
    let a = fillet_edges_overflow(&mut m, &s, &edge, R, end_off, T).unwrap_err();
    let b = fillet_edges_overflow(&mut m, &s, &edge, R, off, T).unwrap_err();
    assert_eq!(a.to_string(), b.to_string());
    assert!(a.to_string().contains("sets back"), "{a}");
    end_off.trim_through = false;
    let err = fillet_edges_overflow(&mut m, &s, &edge, R, end_off, T).unwrap_err();
    assert!(
        err.to_string().contains("hold edge at the end is off"),
        "{err}"
    );
}

/// Two pockets of radius 5, sunk 8, about `(10, 7)` and `(28, 7)`.
const POCKETS: [(f64, f64, f64); 2] = [(10.0, 7.0, 5.0), (28.0, 7.0, 5.0)];

/// A block 40 by 20 by 15 under the plane `z = 0`, [`POCKETS`] sunk into
/// its top, its front top edge along `y = 0` asked.
fn two_pockets(m: &mut Model) -> (Shape, Vec<Shape>) {
    let mut s = ogeom_algo::make_box(m, at(0.0, 0.0, -15.0), (40.0, 20.0, 15.0), T)
        .unwrap()
        .shape;
    for (x, y, r) in POCKETS {
        let pocket = ogeom_algo::make_cylinder(m, at(x, y, -8.0), r, 9.0, T)
            .unwrap()
            .shape;
        s = ogeom_bool::cut(m, &s, &pocket, T).unwrap().shape;
    }
    let edge = edges_where(m, &s, |p| p.y.abs() < 1e-9 && p.z.abs() < 1e-9);
    assert_eq!(edge.len(), 1, "the front top edge");
    (s, edge)
}

/// Two held runs on one band: the front top edge of [`two_pockets`]
/// rounded at 4. Between `x = 6` and `14`, and between `24` and `32`, the
/// contact `y = 4` falls inside a rim and the ball rests on it; both rims
/// are kept there.
fn hold_edge_holds_two_runs_on_one_band_case(heavy: bool) {
    let mut m = Model::new();
    let (s, edge) = two_pockets(&mut m);
    let before = volume(&m, &s);
    let built = fillet_edges(&mut m, &s, &edge, R, T).unwrap();
    holds(
        heavy,
        &m,
        before,
        &built,
        &edge[0],
        &Reference {
            top: Top::Pockets(&POCKETS),
        },
        &[
            (0.0, 6.0),
            (6.0, 14.0),
            (14.0, 24.0),
            (24.0, 32.0),
            (32.0, 40.0),
        ],
        40.0,
    );
    for (cx, cy, r) in POCKETS {
        let on = |p: Point| p.z.abs() < 1e-6 && ((p.x - cx).hypot(p.y - cy) - r).abs() < 1e-6;
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for e in edges_where(&m, &built.shape, on) {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(&m, &e, T).unwrap().unwrap();
            for k in 0..=50 {
                let p = c.point_at(a + (b - a) * f64::from(k) / 50.0, T).unwrap();
                if p.y < R + 1e-6 {
                    lo = lo.min(p.x);
                    hi = hi.max(p.x);
                }
            }
        }
        // The feet run from `(cx - 4, 4)` round the rim's front to
        // `(cx + 4, 4)`.
        assert!(
            lo <= cx - 4.0 + 1e-3 && hi >= cx + 4.0 - 1e-3,
            "the rim about x = {cx} is kept over {lo}..{hi}"
        );
    }
}

#[test]
fn hold_edge_holds_two_runs_on_one_band() {
    hold_edge_holds_two_runs_on_one_band_case(false);
}

/// The two runs' volume against the reference.
#[test]
#[ignore = "heavy"]
fn two_runs_cut_away_the_reference() {
    hold_edge_holds_two_runs_on_one_band_case(true);
}

/// A square pocket beside the edge: its contact crosses the pocket's side
/// rims square to the edge, and between them the ball would rest on two
/// rims at once, its centre turning a corner. That is not held: the request
/// gets trim through's solid, the same as with hold edge off.
#[test]
fn a_ball_on_two_rims_at_once_is_not_held() {
    let mut m = Model::new();
    let block = ogeom_algo::make_box(&mut m, at(0.0, 0.0, -15.0), (30.0, 20.0, 15.0), T)
        .unwrap()
        .shape;
    let pocket = ogeom_algo::make_box(&mut m, at(10.0, 2.0, -8.0), (10.0, 8.0, 9.0), T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(&mut m, &block, &pocket, T).unwrap().shape;
    let edge = edges_where(&m, &s, |p| p.y.abs() < 1e-9 && p.z.abs() < 1e-9);
    assert_eq!(edge.len(), 1, "the front top edge");
    let built = fillet_edges(&mut m, &s, &edge, R, T).unwrap();
    let mut off = Overflow::default();
    off.hold_edge = false;
    let trimmed = fillet_edges_overflow(&mut m, &s, &edge, R, off, T).unwrap();
    assert!(ogeom_algo::check(&m, &built.shape, T).unwrap().is_valid());
    assert!((volume(&m, &built.shape) - volume(&m, &trimmed.shape)).abs() < 1e-9);
}

/// One pocket of radius 5, sunk 8, about `(15, 7)`.
const POCKET: [(f64, f64, f64); 1] = [(15.0, 7.0, 5.0)];

/// The pocket's reference, the run over `x` from 11 to 19 on a band 30
/// long: the band's ends square to it.
fn pocket_reference() -> (f64, f64) {
    let reference = Reference {
        top: Top::Pockets(&POCKET),
    };
    let ranges = [(0.0, 11.0), (11.0, 19.0), (19.0, 30.0)];
    let (coarse, fine) = (reference.volume(&ranges, 24), reference.volume(&ranges, 48));
    (fine, (fine - coarse).abs())
}

/// A chain: a block 30 by 20 by 15 under `z = 0` whose top is two faces
/// split at `x = 5` and whose front is two split at `x = 25`, so its front
/// top edge is three edges end to end; the pocket of [`POCKET`] sunk into
/// the top. All three asked together round as the one edge would: the
/// ball rests on the rim between `x = 11` and `19`, the rim is kept, and
/// the volume cut away is the one edge's reference.
fn hold_edge_holds_a_chain_of_edges_along_one_line_case(heavy: bool) {
    let mut m = Model::new();
    let p = |x: f64, y: f64, z: f64| Point::new(x, y, z);
    let points = [
        p(0.0, 0.0, 0.0),
        p(5.0, 0.0, 0.0),
        p(25.0, 0.0, 0.0),
        p(30.0, 0.0, 0.0),
        p(0.0, 20.0, 0.0),
        p(5.0, 20.0, 0.0),
        p(30.0, 20.0, 0.0),
        p(0.0, 0.0, -15.0),
        p(25.0, 0.0, -15.0),
        p(30.0, 0.0, -15.0),
        p(0.0, 20.0, -15.0),
        p(30.0, 20.0, -15.0),
    ];
    let rings = vec![
        vec![0, 1, 5, 4],
        vec![1, 2, 3, 6, 5],
        vec![7, 10, 11, 9, 8],
        vec![0, 7, 8, 2, 1],
        vec![2, 8, 9, 3],
        vec![6, 11, 10, 4, 5],
        vec![0, 4, 10, 7],
        vec![3, 9, 11, 6],
    ];
    let block = faceted(&mut m, &points, &rings);
    let pocket = ogeom_algo::make_cylinder(&mut m, at(15.0, 7.0, -8.0), 5.0, 9.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(&mut m, &block, &pocket, T).unwrap().shape;
    let chain = edges_where(&m, &s, |p| p.y.abs() < 1e-9 && p.z.abs() < 1e-9);
    assert_eq!(chain.len(), 3, "the front top edge in three");
    let before = volume(&m, &s);
    let built = fillet_edges(&mut m, &s, &chain, R, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    for e in &chain {
        assert!(!built.history.generated(e).is_empty(), "the round");
    }
    let (crease, joins) = worst_join(&m, &built.shape);
    eprintln!("chain: worst join {crease:.1e} rad over {joins} edges");
    assert!(joins >= 2 && crease < 1e-5, "{crease} over {joins}");
    let on = |p: Point| p.z.abs() < 1e-6 && ((p.x - 15.0).hypot(p.y - 7.0) - 5.0).abs() < 1e-6;
    assert!(
        !edges_where(&m, &built.shape, on).is_empty(),
        "the rim is kept"
    );
    if heavy {
        let (want, bound) = pocket_reference();
        let got = before - volume(&m, &built.shape);
        eprintln!("chain: cut away {got:.9}, reference {want:.9}, bound {bound:.1e}");
        assert!((got - want).abs() < 1e-6 + bound, "{got} against {want}");
    }
}

#[test]
fn hold_edge_holds_a_chain_of_edges_along_one_line() {
    hold_edge_holds_a_chain_of_edges_along_one_line_case(false);
}

/// The chain's volume against the one edge's reference.
#[test]
#[ignore = "heavy"]
fn a_chain_cuts_away_the_one_edge_s_reference() {
    hold_edge_holds_a_chain_of_edges_along_one_line_case(true);
}

/// The spandrel's first moment about the edge, `R^3 (5/6 - pi/4)`: about
/// either face's plane, the band's cross-section being symmetric.
fn spandrel_moment() -> f64 {
    R * R * R * (5.0 / 6.0 - PI / 4.0)
}

/// Oblique end faces: the block and pocket of the chain's case, one end
/// face the plane `x = y / 4` and the other `x = 30 + 3 z / 10`. The band
/// runs on past the edge's ends and each end face's plane cuts it in an
/// ellipse. Near the ends the cut is the band's spandrel over the length
/// between the end planes, `30 + 0.3 z - 0.25 y`, so the volume is the
/// square-ended reference less `0.55` times the spandrel's moment.
fn hold_edge_holds_a_band_ending_on_oblique_end_faces_case(heavy: bool) {
    let mut m = Model::new();
    let p = |x: f64, y: f64, z: f64| Point::new(x, y, z);
    let points = [
        p(0.0, 0.0, 0.0),
        p(30.0, 0.0, 0.0),
        p(30.0, 20.0, 0.0),
        p(5.0, 20.0, 0.0),
        p(0.0, 0.0, -15.0),
        p(25.5, 0.0, -15.0),
        p(25.5, 20.0, -15.0),
        p(5.0, 20.0, -15.0),
    ];
    let rings = vec![
        vec![0, 1, 2, 3],
        vec![4, 7, 6, 5],
        vec![0, 4, 5, 1],
        vec![2, 6, 7, 3],
        vec![0, 3, 7, 4],
        vec![1, 5, 6, 2],
    ];
    let block = faceted(&mut m, &points, &rings);
    let pocket = ogeom_algo::make_cylinder(&mut m, at(15.0, 7.0, -8.0), 5.0, 9.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(&mut m, &block, &pocket, T).unwrap().shape;
    let edge = edges_where(&m, &s, |p| p.y.abs() < 1e-9 && p.z.abs() < 1e-9);
    assert_eq!(edge.len(), 1, "the front top edge");
    let before = volume(&m, &s);
    let built = fillet_edges(&mut m, &s, &edge, R, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    for shell in explore_unique(&m, &built.shape, ShapeType::Shell).unwrap() {
        assert!(ogeom_algo::is_shell_closed(&m, &shell).unwrap(), "closed");
    }
    let (crease, joins) = worst_join(&m, &built.shape);
    eprintln!("oblique band ends: worst join {crease:.1e} rad over {joins} edges");
    assert!(joins >= 2 && crease < 1e-5, "{crease} over {joins}");
    if heavy {
        let (square, bound) = pocket_reference();
        let want = square - 0.55 * spandrel_moment();
        let got = before - volume(&m, &built.shape);
        eprintln!("oblique band ends: cut away {got:.9}, reference {want:.9}, bound {bound:.1e}");
        assert!((got - want).abs() < 1e-6 + bound, "{got} against {want}");
    }
}

#[test]
fn hold_edge_holds_a_band_ending_on_oblique_end_faces() {
    hold_edge_holds_a_band_ending_on_oblique_end_faces_case(false);
}

/// The oblique band ends' volume against the reference.
#[test]
#[ignore = "heavy"]
fn oblique_band_ends_cut_away_the_reference() {
    hold_edge_holds_a_band_ending_on_oblique_end_faces_case(true);
}

/// The solid between the polygons `a` and `b`, corner `k` of one joined to
/// corner `k` of the other, every side planar.
fn prism(m: &mut Model, a: &[Point], b: &[Point]) -> Shape {
    let n = a.len();
    let mut points = a.to_vec();
    points.extend_from_slice(b);
    let mut rings = vec![(0..n).collect::<Vec<_>>(), (n..2 * n).rev().collect()];
    rings.extend((0..n).map(|i| vec![i, n + i, n + (i + 1) % n, (i + 1) % n]));
    faceted(m, &points, &rings)
}

/// The cliff's slope on the part whose far end leans.
const SLOPE: [(f64, f64); 2] = [(0.0, 6.0), (20.0, 2.0)];

/// Hold edge at an oblique end: a block whose top narrows from 6 deep at
/// `x = 0` to 2 at `x = 20`, a cliff 6 high falling from its back edge,
/// its far end face the plane `x = 20 + z / 4`. Past `x = 10` the ball
/// rests on the cliff edge to the end, and the leaning end face cuts that
/// canal. The reference integrates the definition up to the end plane:
/// near the end, over `x = 20 - q^2 / 4`, the column's height is `q^2`.
fn hold_edge_holds_a_canal_ending_on_an_oblique_end_face_case(heavy: bool) {
    const K: f64 = 0.25;
    let mut m = Model::new();
    let at_x0 = |y: f64, z: f64| Point::new(0.0, y, z);
    let at_x1 = |y: f64, z: f64| Point::new(K.mul_add(z, 20.0), y, z);
    let near = [
        at_x0(0.0, 0.0),
        at_x0(0.0, -15.0),
        at_x0(12.0, -15.0),
        at_x0(12.0, -6.0),
        at_x0(6.0, -6.0),
        at_x0(6.0, 0.0),
    ];
    // The cliff's foot meets the leaning end at `x = 18.5`, where the
    // cliff stands 2.3 behind the front.
    let far = [
        at_x1(0.0, 0.0),
        at_x1(0.0, -15.0),
        at_x1(12.0, -15.0),
        at_x1(12.0, -6.0),
        at_x1(2.3, -6.0),
        at_x1(2.0, 0.0),
    ];
    let s = prism(&mut m, &near, &far);
    let edge = edges_where(&m, &s, |p| p.y.abs() < 1e-9 && p.z.abs() < 1e-9);
    assert_eq!(edge.len(), 1, "the front top edge");
    let before = volume(&m, &s);
    let built = fillet_edges(&mut m, &s, &edge, R, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    for shell in explore_unique(&m, &built.shape, ShapeType::Shell).unwrap() {
        assert!(ogeom_algo::is_shell_closed(&m, &shell).unwrap(), "closed");
    }
    let (crease, joins) = worst_join(&m, &built.shape);
    eprintln!("oblique canal end: worst join {crease:.1e} rad over {joins} edges");
    assert!(joins >= 2 && crease < 1e-5, "{crease} over {joins}");
    let (lo, hi) = kept_along(&m, &built.shape, SLOPE[0], SLOPE[1]);
    assert!(
        lo <= 10.0 + 1e-6 && hi >= 19.9,
        "the cliff edge kept over {lo}..{hi}"
    );
    if !heavy {
        return;
    }
    let reference = Reference {
        top: Top::Below(&SLOPE),
    };
    // Where the end plane first clips a column: its height `(20 - x) / K`
    // meets the depth of the ball's contact on the front.
    let (mut lo, mut hi) = (10.0, 20.0);
    for _ in 0..80 {
        let x = 0.5 * (lo + hi);
        if -reference.centre(x).z < (20.0 - x) / K {
            lo = x;
        } else {
            hi = x;
        }
    }
    let x_star = 0.5 * (lo + hi);
    let volume_at = |n: usize| {
        let open = reference.volume(&[(0.0, 10.0), (10.0, x_star)], n);
        let g = gauss(n);
        let q1 = ((20.0 - x_star) / K).sqrt();
        let mut tail = 0.0;
        for (gq, wq) in &g {
            let q = 0.5 * q1 * (gq + 1.0);
            let x = (K * q).mul_add(-q, 20.0);
            for (gs, ws) in &g {
                let s = 0.5 * q * (gs + 1.0);
                tail += wq
                    * ws
                    * 0.5
                    * q1
                    * 0.5
                    * q
                    * reference.depth(x, -s * s)
                    * 2.0
                    * s
                    * 2.0
                    * K
                    * q;
            }
        }
        open + tail
    };
    let (coarse, fine) = (volume_at(24), volume_at(48));
    let bound = (fine - coarse).abs();
    let got = before - volume(&m, &built.shape);
    eprintln!("oblique canal end: cut away {got:.9}, reference {fine:.9}, bound {bound:.1e}");
    assert!((got - fine).abs() < 1e-6 + bound, "{got} against {fine}");
}

#[test]
fn hold_edge_holds_a_canal_ending_on_an_oblique_end_face() {
    hold_edge_holds_a_canal_ending_on_an_oblique_end_face_case(false);
}

/// The oblique canal end's volume against the reference.
#[test]
#[ignore = "heavy"]
fn an_oblique_canal_end_cuts_away_the_reference() {
    hold_edge_holds_a_canal_ending_on_an_oblique_end_face_case(true);
}
