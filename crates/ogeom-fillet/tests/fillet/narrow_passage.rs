//! Narrow passages (`Overflow::narrow_blend_made`, `narrow_whole_face`,
//! `narrow_partial_face`): concave rounds of a floor's rim where the zone
//! between them is tighter than the radius. Each round keeps its constant
//! section and runs on through the zone, and the rounds meet each other
//! where they cross: the fill is the union of the bands' fills.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{FRAC_1_SQRT_2, PI};

use ogeom_core::Tolerances;
use ogeom_fillet::{Overflow, fillet_edges, fillet_edges_overflow};
use ogeom_geom::{
    CircleCurve, Curve, Curve3d as _, LineCurve, PlaneSurface, Surface as _, SurfaceGeometry,
};
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The round's radius in every case.
const R: f64 = 4.0;
/// The floor's height in every part.
const FLOOR: f64 = 7.0;

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

fn valid(m: &Model, s: &Shape) -> bool {
    ogeom_algo::check(m, s, T).unwrap().is_valid()
}

fn surface_of(m: &Model, face: &Shape) -> SurfaceGeometry {
    let data = m.node(face).unwrap().data().as_face().unwrap();
    m.geometry().surface(data.surface).unwrap().clone()
}

/// Every edge lying in the floor's plane.
fn in_floor(m: &Model, s: &Shape) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            (0..=10).all(|k| {
                let p = c.point_at(a + (b - a) * f64::from(k) / 10.0, T).unwrap();
                (p.z - FLOOR).abs() < 1e-9
            })
        })
        .collect()
}

/// The floor's rim edges along its walls: a slot's floor ends on end
/// faces, whose edges with the floor are left sharp.
fn rim(m: &Model, s: &Shape) -> Vec<Shape> {
    in_floor(m, s)
        .into_iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            let (p, q) = (c.point_at(a, T).unwrap(), c.point_at(b, T).unwrap());
            let across = (p.x - q.x).abs() < 1e-9;
            !(across && (p.x.abs() < 1e-9 || (p.x - 20.0).abs() < 1e-9))
        })
        .collect()
}

/// A piece of a plan outline at the floor's height: a segment, or the arc
/// about `c` from `from` to `to` counter-clockwise.
enum Piece {
    Line(f64, f64, f64, f64),
    Arc {
        c: (f64, f64),
        from: (f64, f64),
        to: (f64, f64),
    },
}

/// The solid swept up by `height` from the plane outline `pieces` at the
/// floor's height, walked counter-clockwise end to end.
fn outline(m: &mut Model, pieces: &[Piece], height: f64) -> Shape {
    let start = |p: &Piece| match p {
        Piece::Line(x, y, _, _) => (*x, *y),
        Piece::Arc { from, .. } => *from,
    };
    let corners: Vec<Shape> = pieces
        .iter()
        .map(|p| {
            let (x, y) = start(p);
            ogeom_algo::make_vertex(m, Point::new(x, y, FLOOR)).shape
        })
        .collect();
    let n = pieces.len();
    let mut edges = Vec::with_capacity(n);
    for (k, piece) in pieces.iter().enumerate() {
        let (va, vb) = (&corners[k], &corners[(k + 1) % n]);
        let edge = match piece {
            Piece::Line(x0, y0, x1, y1) => {
                let (a, b) = (Point::new(*x0, *y0, FLOOR), Point::new(*x1, *y1, FLOOR));
                let line = Curve::Line(LineCurve::segment(a, b, T).unwrap());
                ogeom_algo::make_edge_between(m, line, (0.0, a.distance(b)), va, vb, T)
            }
            Piece::Arc { c, from, to } => {
                let radius = (from.0 - c.0).hypot(from.1 - c.1);
                let a0 = (from.1 - c.1).atan2(from.0 - c.0);
                let mut a1 = (to.1 - c.1).atan2(to.0 - c.0);
                while a1 <= a0 {
                    a1 += 2.0 * PI;
                }
                let circle = Circle::new(at(c.0, c.1, FLOOR), radius, T).unwrap();
                ogeom_algo::make_edge_between(
                    m,
                    CircleCurve::new(circle).into(),
                    (a0, a1),
                    va,
                    vb,
                    T,
                )
            }
        };
        edges.push(edge.unwrap().shape);
    }
    let wire = ogeom_algo::make_wire(m, &edges, T).unwrap().shape;
    let plane = PlaneSurface::new(Plane::new(at(0.0, 0.0, FLOOR)));
    let face = ogeom_algo::make_face(m, SurfaceGeometry::Plane(plane), &[wire], T)
        .unwrap()
        .shape;
    ogeom_algo::make_prism(m, &face, Vector::new(0.0, 0.0, height), T)
        .unwrap()
        .shape
}

/// A block 40 x 30 x 15 with the pocket `pieces` sunk to the floor.
fn pocket(m: &mut Model, pieces: &[Piece]) -> Shape {
    let block = ogeom_algo::make_box(m, at(0.0, 0.0, 0.0), (40.0, 30.0, 15.0), T)
        .unwrap()
        .shape;
    let tool = outline(m, pieces, 10.0);
    ogeom_bool::cut(m, &block, &tool, T).unwrap().shape
}

/// A block 20 long with a through slot `w` wide from `y = 10`, its floor
/// at the floor's height: the U section swept along `x`.
fn slot(m: &mut Model, w: f64) -> Shape {
    slot_deep(m, w, 15.0 - FLOOR)
}

/// [`slot`], the slot `depth` deep.
fn slot_deep(m: &mut Model, w: f64, depth: f64) -> Shape {
    let (y0, y1) = (10.0, 10.0 + w);
    let top = FLOOR + depth;
    let section = [
        (0.0, 0.0),
        (30.0, 0.0),
        (30.0, top),
        (y1, top),
        (y1, FLOOR),
        (y0, FLOOR),
        (y0, top),
        (0.0, top),
    ];
    let n = section.len();
    let mut points: Vec<Point> = section
        .iter()
        .map(|(y, z)| Point::new(0.0, *y, *z))
        .collect();
    points.extend(section.iter().map(|(y, z)| Point::new(20.0, *y, *z)));
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
    panic!("the slot's rings enclose no volume either way round");
}

/// The stadium pocket: `x` from 10 to 30, `y` from 12 to 18, its ends
/// half cylinders of radius 3.
fn stadium(m: &mut Model) -> Shape {
    pocket(
        m,
        &[
            Piece::Line(13.0, 12.0, 27.0, 12.0),
            Piece::Arc {
                c: (27.0, 15.0),
                from: (27.0, 12.0),
                to: (27.0, 18.0),
            },
            Piece::Line(27.0, 18.0, 13.0, 18.0),
            Piece::Arc {
                c: (13.0, 15.0),
                from: (13.0, 18.0),
                to: (13.0, 12.0),
            },
        ],
    )
}

/// The rectangular pocket `x` 10 to 30, `y` 8 to 22, each corner an arc
/// about a centre `m` in from both walls, meeting them `k` in
/// from the corner: tangent where `m = k = rc`, else at kinks.
fn cornered(model: &mut Model, m: f64, k: f64) -> Shape {
    let (x0, x1, y0, y1) = (10.0, 30.0, 8.0, 22.0);
    pocket(
        model,
        &[
            Piece::Line(x0 + k, y0, x1 - k, y0),
            Piece::Arc {
                c: (x1 - m, y0 + m),
                from: (x1 - k, y0),
                to: (x1, y0 + k),
            },
            Piece::Line(x1, y0 + k, x1, y1 - k),
            Piece::Arc {
                c: (x1 - m, y1 - m),
                from: (x1, y1 - k),
                to: (x1 - k, y1),
            },
            Piece::Line(x1 - k, y1, x0 + k, y1),
            Piece::Arc {
                c: (x0 + m, y1 - m),
                from: (x0 + k, y1),
                to: (x0, y1 - k),
            },
            Piece::Line(x0, y1 - k, x0, y0 + k),
            Piece::Arc {
                c: (x0 + m, y0 + m),
                from: (x0, y0 + k),
                to: (x0 + k, y0),
            },
        ],
    )
}

// The independent reference. Every case's fill is the union of the bands'
// fills, cut to the passage's air: over a point of the floor at distance
// `s` from the nearest wall plane, the fill stands `H(s)` tall, the height
// of the round's quarter circle there. The volume is integrated over the
// air's plan by Gauss-Legendre quadrature in one variable, with `s = t^2`
// where a square root starts at `s = 0`.

/// The round's height over the floor at distance `s` from its wall.
fn h(s: f64) -> f64 {
    if s >= R {
        return 0.0;
    }
    R - (R * R - (R - s).powi(2)).max(0.0).sqrt()
}

/// Gauss-Legendre nodes and weights on `[-1, 1]`.
fn gauss(n: usize) -> Vec<(f64, f64)> {
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        #[expect(clippy::cast_precision_loss, reason = "a few dozen nodes")]
        let mut x = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 1.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                #[expect(clippy::cast_precision_loss, reason = "a few dozen nodes")]
                let kf = k as f64;
                let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
                p0 = p1;
                p1 = p2;
            }
            #[expect(clippy::cast_precision_loss, reason = "a few dozen nodes")]
            let nf = n as f64;
            dp = nf * (x * p1 - p0) / (x * x - 1.0);
            let step = p1 / dp;
            x -= step;
            if step.abs() < 1e-16 {
                break;
            }
        }
        out.push((x, 2.0 / ((1.0 - x * x) * dp * dp)));
    }
    out
}

/// `f(s)` integrated over `s` in `[0, top]` with `s = t^2`, on `n` points.
fn integral(f: impl Fn(f64) -> f64, top: f64, n: usize) -> f64 {
    let half = top.sqrt() / 2.0;
    gauss(n)
        .into_iter()
        .map(|(x, w)| {
            let t = half * (x + 1.0);
            w * half * f(t * t) * 2.0 * t
        })
        .sum()
}

/// A reference and its own error bound, the change from `n` to `2n` points.
fn reference(f: impl Fn(usize) -> f64) -> (f64, f64) {
    let (coarse, fine) = (f(24), f(48));
    (fine, (fine - coarse).abs())
}

/// The fill of a rectangle `wide` by `long` whose walls are all rounded:
/// the strip `s` in from its walls is `2 (wide + long) - 8 s` long.
fn rectangle_fill(wide: f64, long: f64, n: usize) -> f64 {
    integral(|s| h(s) * 8.0f64.mul_add(-s, 2.0 * (wide + long)), R, n)
}

/// The fill over one corner's material, where an arc of radius `rc` about
/// the point `m` in from both walls bounds it: by symmetry twice the part
/// nearer the first wall, its height `H` of the distance `a` to it.
fn corner_fill(rc: f64, m: f64, n: usize) -> f64 {
    let diagonal = rc.mul_add(-FRAC_1_SQRT_2, m);
    let below = |a: f64| m - (rc * rc - (a - m).powi(2)).max(0.0).sqrt();
    2.0 * integral(|a| h(a) * (below(a) - a).max(0.0), diagonal, n)
}

/// The largest angle between the normals either side of each edge where a
/// round (a cylinder of radius `R`) meets a wall (a plane square to the
/// floor and holding the round's axis direction), sampled along the edge; and the joins' count. These are G1.
fn worst_wall_join(m: &Model, s: &Shape) -> (f64, usize) {
    let mut worst = 0.0_f64;
    let mut joins = 0;
    for e in explore_unique(m, s, ShapeType::Edge).unwrap() {
        let faces = ancestors_of(m, s, &e, ShapeType::Face).unwrap();
        let [a, b] = &faces[..] else {
            continue;
        };
        let (sa, sb) = (surface_of(m, a), surface_of(m, b));
        let wall = match (&sa, &sb) {
            (SurfaceGeometry::Cylinder(c), SurfaceGeometry::Plane(p))
            | (SurfaceGeometry::Plane(p), SurfaceGeometry::Cylinder(c))
                if (c.cylinder().radius() - R).abs() < 1e-9 =>
            {
                let n = p.plane().normal().vector();
                n.dot(Vector::Z).abs() < 1e-9
                    && n.dot(c.cylinder().axis().direction.vector()).abs() < 1e-9
            }
            _ => false,
        };
        if !wall {
            continue;
        }
        joins += 1;
        let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap().unwrap();
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

/// Whether any face of `s` lies in the floor's plane.
fn floor_left(m: &Model, s: &Shape) -> bool {
    explore_unique(m, s, ShapeType::Face)
        .unwrap()
        .iter()
        .any(|f| {
            matches!(surface_of(m, f), SurfaceGeometry::Plane(p)
            if p.plane().normal().vector().cross(Vector::Z).magnitude() < 1e-9
                && (p.plane().frame().origin().z - FLOOR).abs() < 1e-9)
        })
}

/// Builds the case, checks it and its oracle, and reports. Returns the
/// fill the round added.
fn check_case(name: &str, m: &mut Model, s: &Shape, overflow: Overflow, oracle: (f64, f64)) -> f64 {
    let edges = rim(m, s);
    let before = volume(m, s);
    let built = fillet_edges_overflow(m, s, &edges, R, overflow, T).unwrap();
    assert!(valid(m, &built.shape), "{name}: the result checks invalid");
    assert!(
        ogeom_algo::is_shell_closed(
            m,
            &explore_unique(m, &built.shape, ShapeType::Shell).unwrap()[0]
        )
        .unwrap()
    );
    let added = volume(m, &built.shape) - before;
    let (expected, bound) = oracle;
    let (worst, joins) = worst_wall_join(m, &built.shape);
    println!(
        "{name}: fill {added:.9}, reference {expected:.9}, error {:.1e}, bound {bound:.1e}, \
         worst G1 {worst:.1e} rad over {joins} round-to-wall edges",
        (added - expected).abs()
    );
    assert!(
        (added - expected).abs() <= 1e-6 + bound,
        "{name}: fill {added} against {expected}"
    );
    assert!(
        joins > 0 && worst < 1e-5,
        "{name}: rounds meet their walls at {worst} rad"
    );
    added
}

#[test]
fn a_through_slot_narrower_than_twice_the_radius() {
    let w = 6.0;
    let mut m = Model::new();
    let s = slot(&mut m, w);
    // Closed form: the two quarter circles' spandrels meet mid-slot.
    let g = |u: f64| (u * (R * R - u * u).max(0.0).sqrt() + R * R * (u / R).asin()) / 2.0;
    let area = 2.0 * R.mul_add(w / 2.0, -(g(R) - g(R - w / 2.0)));
    let added = check_case(
        "through slot w 6",
        &mut m,
        &s,
        Overflow::default(),
        (area * 20.0, 0.0),
    );
    assert!(added > 0.0);
}

#[test]
fn a_through_slot_keeps_its_wider_round() {
    // Twice the radius or wider is no narrow passage: every strategy off
    // builds the same two bands as the default.
    for overflow in [Overflow::default(), narrow_off()] {
        let mut m = Model::new();
        let s = slot(&mut m, 10.0);
        let edges = rim(&m, &s);
        let before = volume(&m, &s);
        let built = fillet_edges_overflow(&mut m, &s, &edges, R, overflow, T).unwrap();
        let added = volume(&m, &built.shape) - before;
        let spandrel = R * R - PI * R * R / 4.0;
        assert!((added - 2.0 * spandrel * 20.0).abs() < 1e-6);
    }
}

#[test]
fn a_stadium_pocket_narrower_than_twice_the_radius() {
    let mut m = Model::new();
    let s = stadium(&mut m);
    assert_eq!(rim(&m, &s).len(), 4);
    let oracle = reference(|n| {
        let (rc, straight) = (3.0, 14.0);
        2.0 * integral(
            |s| h(s) * 2.0f64.mul_add((rc * rc - (rc - s).powi(2)).max(0.0).sqrt(), straight),
            rc,
            n,
        )
    });
    check_case("stadium w 6", &mut m, &s, Overflow::default(), oracle);
    let mut m = Model::new();
    let s = stadium(&mut m);
    let edges = rim(&m, &s);
    let built = fillet_edges(&mut m, &s, &edges, R, T).unwrap();
    assert!(!floor_left(&m, &built.shape), "the floor is taken whole");
}

#[test]
fn a_pocket_whose_corner_rounds_are_tighter_than_the_radius() {
    let mut m = Model::new();
    let s = cornered(&mut m, 2.0, 2.0);
    assert_eq!(rim(&m, &s).len(), 8);
    let oracle = reference(|n| rectangle_fill(20.0, 14.0, n) - 4.0 * corner_fill(2.0, 2.0, n));
    check_case(
        "pocket corners rc 2",
        &mut m,
        &s,
        Overflow::default(),
        oracle,
    );
}

/// The kinked corner: radius 3, meeting each wall 2 from the corner.
fn kinked_centre() -> f64 {
    let (rc, k) = (3.0_f64, 2.0_f64);
    (k + 2.0f64.mul_add(rc * rc, -(k * k)).sqrt()) / 2.0
}

#[test]
fn a_pocket_whose_corners_are_kinked_arcs() {
    let mut m = Model::new();
    let s = cornered(&mut m, kinked_centre(), 2.0);
    assert_eq!(rim(&m, &s).len(), 8);
    let oracle =
        reference(|n| rectangle_fill(20.0, 14.0, n) - 4.0 * corner_fill(3.0, kinked_centre(), n));
    check_case(
        "pocket kinked corners rc 3",
        &mut m,
        &s,
        Overflow::default(),
        oracle,
    );
}

/// Every narrow-passage strategy off, the rest as by default.
fn narrow_off() -> Overflow {
    let mut overflow = Overflow::default();
    overflow.narrow_blend_made = false;
    overflow.narrow_whole_face = false;
    overflow.narrow_partial_face = false;
    overflow
}

fn refusal(m: &mut Model, s: &Shape, overflow: Overflow) -> String {
    let edges = rim(m, s);
    fillet_edges_overflow(m, s, &edges, R, overflow, T)
        .expect_err("refused by name")
        .to_string()
}

#[test]
fn each_strategy_off_passes_a_blend_made_passage_on() {
    let reference_fill = {
        let mut m = Model::new();
        let s = slot(&mut m, 6.0);
        let before = volume(&m, &s);
        let edges = rim(&m, &s);
        let built = fillet_edges(&mut m, &s, &edges, R, T).unwrap();
        volume(&m, &built.shape) - before
    };
    let mut blend_off = Overflow::default();
    blend_off.narrow_blend_made = false;
    let mut partial_only = blend_off;
    partial_only.narrow_whole_face = false;
    for overflow in [blend_off, partial_only] {
        let mut m = Model::new();
        let s = slot(&mut m, 6.0);
        let before = volume(&m, &s);
        let edges = rim(&m, &s);
        let built = fillet_edges_overflow(&mut m, &s, &edges, R, overflow, T).unwrap();
        assert!((volume(&m, &built.shape) - before - reference_fill).abs() < 1e-9);
    }
    let mut m = Model::new();
    let s = slot(&mut m, 6.0);
    let why = refusal(&mut m, &s, narrow_off());
    assert!(
        why.contains("every narrow-passage strategy is off"),
        "{why}"
    );
    let mut m = Model::new();
    let s = slot(&mut m, 6.0);
    assert!(refusal(&mut m, &s, Overflow::none()).contains("narrow-passage"));
}

#[test]
fn a_kinked_corner_needs_whole_face_or_partial_face() {
    let mut blend_only = narrow_off();
    blend_only.narrow_blend_made = true;
    let mut m = Model::new();
    let s = cornered(&mut m, kinked_centre(), 2.0);
    let why = refusal(&mut m, &s, blend_only);
    assert!(why.contains("whole face and partial face are off"), "{why}");
    let mut partial_only = narrow_off();
    partial_only.narrow_partial_face = true;
    let mut m = Model::new();
    let s = cornered(&mut m, kinked_centre(), 2.0);
    let before = volume(&m, &s);
    let edges = rim(&m, &s);
    let built = fillet_edges_overflow(&mut m, &s, &edges, R, partial_only, T).unwrap();
    let expected = rectangle_fill(20.0, 14.0, 48) - 4.0 * corner_fill(3.0, kinked_centre(), 48);
    assert!((volume(&m, &built.shape) - before - expected).abs() < 1e-5);
    let mut m = Model::new();
    let s = cornered(&mut m, kinked_centre(), 2.0);
    let why = refusal(&mut m, &s, narrow_off());
    assert!(why.contains("whole face and partial face are off"), "{why}");
}

/// `f` integrated over `[lo, hi]` on `n` Gauss-Legendre points.
fn plain(f: impl Fn(f64) -> f64, lo: f64, hi: f64, n: usize) -> f64 {
    let half = (hi - lo) / 2.0;
    gauss(n)
        .into_iter()
        .map(|(x, w)| w * half * f(half.mul_add(x + 1.0, lo)))
        .sum()
}

#[test]
fn the_reference_reproduces_the_closed_forms() {
    // Walls with no corner faces: the strip integral of the walls'
    // perimeter against its spandrels' closed form, before the corners.
    let spandrel = R * R - PI * R * R / 4.0;
    let along = integral(|s| h(s) * 2.0 * (20.0 + 14.0), R, 48);
    assert!((along - 2.0 * (20.0 + 14.0) * spandrel).abs() < 1e-9);
    // A corner's material: the symmetric reading against both halves
    // integrated directly, `H(min(a, c))` over `c` under the arc at each
    // `a`, split where the arc crosses the diagonal.
    for (rc, m, k) in [(2.0, 2.0, 2.0), (3.0, kinked_centre(), 2.0)] {
        let under = |a: f64| m - (rc * rc - (a - m).powi(2)).max(0.0).sqrt();
        let column = |a: f64| {
            let top = under(a);
            let low = integral(h, a.min(top), 24);
            low + h(a) * (top - a).max(0.0)
        };
        let diagonal = rc.mul_add(-FRAC_1_SQRT_2, m);
        let direct = integral(column, diagonal, 48) + plain(column, diagonal, k, 48);
        let symmetric = corner_fill(rc, m, 48);
        assert!((symmetric - direct).abs() < 1e-9, "{symmetric} {direct}");
    }
}

#[test]
fn passages_the_union_does_not_round_are_refused_by_name() {
    // A slot no wider than the radius: a band would cross the wall it
    // faces. A slot shallower than the radius: the walls end inside the
    // bands' room. Both go on to the engine, which refuses them by name.
    for (w, depth) in [(3.0, 8.0), (6.0, 3.0)] {
        let mut m = Model::new();
        let s = slot_deep(&mut m, w, depth);
        let edges = rim(&m, &s);
        assert_eq!(edges.len(), 2);
        let why = fillet_edges(&mut m, &s, &edges, R, T)
            .expect_err("refused by name")
            .to_string();
        assert!(why.contains("radius 4") || why.contains("blend"), "{why}");
    }
}

/// `fillet_edges_overflow` on the slot `w` wide, run on its own thread
/// and given `limit` to answer: the fill it added, or its refusal.
fn bounded(w: f64, overflow: Overflow, limit: std::time::Duration) -> Result<f64, String> {
    bounded_on(w, overflow, false, limit)
}

/// [`bounded`], on the slot's first floor edge alone where `one`.
fn bounded_on(
    w: f64,
    overflow: Overflow,
    one: bool,
    limit: std::time::Duration,
) -> Result<f64, String> {
    let (send, receive) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let mut m = Model::new();
        let s = slot(&mut m, w);
        let mut edges = rim(&m, &s);
        if one {
            edges.truncate(1);
        }
        let before = volume(&m, &s);
        let answer = fillet_edges_overflow(&mut m, &s, &edges, R, overflow, T).map(|built| {
            assert!(valid(&m, &built.shape), "the result checks invalid");
            let (worst, joins) = worst_wall_join(&m, &built.shape);
            println!("  worst G1 {worst:.1e} rad over {joins} round-to-wall edges");
            assert!(
                joins > 0 && worst < 1e-5,
                "rounds meet their walls at {worst} rad"
            );
            volume(&m, &built.shape) - before
        });
        let _ = send.send(answer.map_err(|e| e.to_string()));
    });
    receive
        .recv_timeout(limit)
        .unwrap_or_else(|_| panic!("the slot {w} wide gave no answer within {limit:?}"))
}

#[test]
fn a_slot_exactly_twice_the_radius_wide_answers_in_time() {
    // The two contacts meet on the floor's centre line: the floor is taken
    // whole and the two quarter circles make one half circle, the trough's
    // section, its area twice the spandrel's.
    let limit = std::time::Duration::from_secs(60);
    let trough = 2.0 * (R * R - PI * R * R / 4.0) * 20.0;
    let added = bounded(2.0 * R, Overflow::default(), limit).expect("the trough builds");
    println!("2r slot: fill {added:.9} against {trough:.9}");
    assert!(
        (added - trough).abs() < 1e-6,
        "fill {added} against {trough}"
    );
    // Every narrow-passage strategy off: an answer in time all the same, a
    // refusal by name or the same trough.
    let off = bounded(2.0 * R, narrow_off(), limit);
    println!("2r slot, narrow strategies off: {off:?}");
    match off {
        Ok(added) => assert!(
            (added - trough).abs() < 1e-6,
            "fill {added} against {trough}"
        ),
        Err(why) => assert!(!why.is_empty()),
    }
}

#[test]
fn one_round_of_a_slot_twice_the_radius_wide_answers_in_time() {
    // The ball stands tangent to the far wall all along: the proof that
    // it clears that wall has no margin anywhere and must still end. The
    // band alone fits (its contact lies on the floor's centre line), so
    // it builds the one spandrel, or is refused by name.
    let limit = std::time::Duration::from_secs(60);
    let spandrel = (R * R - PI * R * R / 4.0) * 20.0;
    let answer = bounded_on(2.0 * R, Overflow::default(), true, limit);
    println!("one round of the 2r slot: {answer:?}");
    if let Ok(added) = answer {
        assert!(
            (added - spandrel).abs() < 1e-6,
            "fill {added} against {spandrel}"
        );
    }
}
