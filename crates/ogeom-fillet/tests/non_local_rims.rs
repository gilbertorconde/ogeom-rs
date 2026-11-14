//! A rim's round reaching a face it does not end at, through the blend
//! engine's one boolean (ADR 0001, M5a and M5c): each arc's strip closed
//! into its tool (the strip, its flap on each face, a cap where the arcs
//! meet) and cut or fused.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{FRAC_PI_4, PI};

use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

/// The rim of radius `r` at height `z`: every edge on that circle.
fn rim(m: &Model, s: &Shape, z: f64, r: f64) -> Vec<Shape> {
    let on = |p: Point| (p.z - z).abs() < 1e-6 && (p.x.hypot(p.y) - r).abs() < 1e-6;
    let found: Vec<Shape> = explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            [0.0, 0.3, 0.7, 1.0]
                .iter()
                .all(|t| on(ogeom_geom::Curve3d::point_at(&c, a + (b - a) * t, T).unwrap()))
        })
        .collect();
    assert!(!found.is_empty(), "the rim");
    found
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

/// The integral of `f` over the rectangle `x0..x1` by `y0..y1`.
fn over(f: impl Fn(f64, f64) -> f64, (x0, x1): (f64, f64), (y0, y1): (f64, f64)) -> f64 {
    let g = gauss(48);
    let (hx, hy) = (0.5 * (x1 - x0), 0.5 * (y1 - y0));
    let mut sum = 0.0;
    for (a, wa) in &g {
        for (b, wb) in &g {
            sum += wa * wb * f(x0 + hx * (a + 1.0), y0 + hy * (b + 1.0));
        }
    }
    sum * hx * hy
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// A plate 40 by 40 by 4 with a boss of radius 5, 5 tall, standing on its
/// middle, less a pocket `0.5` square (`x` from `x0`, `y` within a quarter
/// of nought) from `z0` to `z1`.
fn bossed(m: &mut Model, x0: f64, (z0, z1): (f64, f64)) -> Shape {
    let plate = ogeom_algo::make_box(m, at(-20.0, -20.0, 0.0), (40.0, 40.0, 4.0), T)
        .unwrap()
        .shape;
    let boss = ogeom_algo::make_cylinder(m, at(0.0, 0.0, 4.0), 5.0, 5.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(m, &plate, &boss, T).unwrap().shape;
    let pocket = ogeom_algo::make_box(m, at(x0, -0.25, z0), (0.5, 0.5, z1 - z0), T)
        .unwrap()
        .shape;
    ogeom_bool::cut(m, &s, &pocket, T).unwrap().shape
}

/// The rounded part, valid, the rim generating its strip, its volume `want`.
fn rounds_to(m: &mut Model, s: &Shape, edges: &[Shape], want: f64) {
    let built = ogeom_fillet::fillet_edges(m, s, edges, 1.0, T).unwrap();
    let found = ogeom_algo::check(m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    let got = volume(m, &built.shape);
    assert!((got - want).abs() < 1e-6, "{got} against {want}");
    for edge in edges {
        assert!(!built.history.generated(edge).is_empty(), "the rim's strip");
    }
}

/// The boss's top rim rounded at 1 where a pocket is sunk 3 into the boss
/// beside it (`x` 4.3 to 4.8): the round's spandrel reaches into the
/// pocket, whose walls are no face of the rim's. The spandrel is cut by the
/// one boolean: the part less the pocket and the spandrel's ring (its
/// section `(1 - pi/4)` about a centroid `(10 - 3 pi) / (12 - 3 pi)` in from
/// the rim, by Pappus), less what the two share (the spandrel's height
/// `1 - sqrt(1 - (rho - 4)^2)` over the pocket's floor plan, by quadrature).
#[test]
fn a_boss_s_top_rim_rounds_past_a_pocket_beside_it() {
    let mut m = Model::new();
    let s = bossed(&mut m, 4.3, (6.0, 10.0));
    let edges = rim(&m, &s, 9.0, 5.0);
    let section = 1.0 - FRAC_PI_4;
    let centroid = (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI);
    let ring = 2.0 * PI * (5.0 - centroid) * section;
    let shared = over(
        |x, y| {
            let rho = x.hypot(y);
            1.0 - (1.0 - (rho - 4.0).powi(2)).sqrt()
        },
        (4.3, 4.8),
        (-0.25, 0.25),
    );
    let want = 6400.0 + 125.0 * PI - 0.75 - ring + shared;
    rounds_to(&mut m, &s, &edges, want);
}

/// The boss's foot rounded at 1 where a pocket goes through the plate
/// beside it (`x` 5.3 to 5.8): the round's fill reaches over the pocket's
/// mouth, whose walls are no face of the rim's. The fill is fused by the
/// one boolean, keeping its surface: it stands above the plate, so it roofs
/// the mouth where it reaches over it and the pocket below is kept. The
/// part with the fill's whole ring (its section `(1 - pi/4)` about a
/// centroid `(10 - 3 pi) / (12 - 3 pi)` out from the foot, by Pappus), and
/// over the pocket the result stands as high as the fill (`1 - sqrt(1 -
/// (rho - 6)^2)` over the plate).
#[test]
fn a_boss_s_foot_rounds_over_a_pocket_beside_it() {
    let mut m = Model::new();
    let s = bossed(&mut m, 5.3, (-1.0, 5.0));
    let edges = rim(&m, &s, 4.0, 5.0);
    let section = 1.0 - FRAC_PI_4;
    let centroid = (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI);
    let ring = 2.0 * PI * (5.0 + centroid) * section;
    let want = 6400.0 + 125.0 * PI - 1.0 + ring;
    rounds_to(&mut m, &s, &edges, want);
}

/// A drum of radius 5, 10 tall, alone, its top rim rounded at 1 where a
/// pocket is sunk 3 into its top beside the rim (`x` 4.3 to 4.8): each of
/// the rim's arcs closed into its tool, capped where the arcs meet. The
/// drum less the pocket and the spandrel's ring, less what the two share,
/// as for the boss.
#[test]
fn a_drum_s_top_rim_rounds_past_a_pocket_beside_it() {
    let mut m = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut m, at(0.0, 0.0, 0.0), 5.0, 10.0, T)
        .unwrap()
        .shape;
    let pocket = ogeom_algo::make_box(&mut m, at(4.3, -0.25, 7.0), (0.5, 0.5, 4.0), T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(&mut m, &drum, &pocket, T).unwrap().shape;
    let edges = rim(&m, &s, 10.0, 5.0);
    let section = 1.0 - FRAC_PI_4;
    let centroid = (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI);
    let ring = 2.0 * PI * (5.0 - centroid) * section;
    let shared = over(
        |x, y| {
            let rho = x.hypot(y);
            1.0 - (1.0 - (rho - 4.0).powi(2)).sqrt()
        },
        (4.3, 4.8),
        (-0.25, 0.25),
    );
    let want = 250.0 * PI - 0.75 - ring + shared;
    rounds_to(&mut m, &s, &edges, want);
}
