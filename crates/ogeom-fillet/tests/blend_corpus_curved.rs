//! The acceptance corpus for blends on curved faces (ADR 0001, M3, §12).
//!
//! Each case asks the public calls and holds the answer to an exact oracle:
//! Pappus for a rim's torus (the corner's area swept round its centroid),
//! a prism's section for a straight edge running along a curved face, and
//! closed forms for chamfers; where the result holds a general vertex patch
//! (the cusp S-blend), tangency within the bands' fit target instead. Every
//! case is also held to:
//!
//! - every order of the request giving the bit-identical solid;
//! - validity;
//! - every face proved within the original's box, a concave round's fill
//!   only where the case says.
//!
//! The cases run only with `OGEOM_M3_CORPUS` set, until rounds on curved
//! faces route through the engine (ADR 0001, §11): today's path fails each,
//! refusing it or answering differently by the request's order, as each
//! case's documentation records. Without it each returns at once.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom_algo::HalfSpaces;
use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// How far a face may stand past the region it is held to, and how far a
/// blend may meet its faces off tangent: the bands' fit target.
const FIT: f64 = 2e-4;

/// Whether the corpus is asked for: until rounds on curved faces route
/// through the engine, only where `OGEOM_M3_CORPUS` is set.
fn asked() -> bool {
    std::env::var_os("OGEOM_M3_CORPUS").is_some()
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(
        m,
        s,
        ogeom_mesh::Deflection {
            chord: 1e-5,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

fn point_of(m: &Model, v: &Shape) -> Point {
    m.node(v)
        .and_then(|n| n.data().as_vertex().map(|d| d.point))
        .unwrap()
}

/// A frame at `(x, y, z)`, its `z` up.
fn up_at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(
        Point::new(x, y, z),
        Direction::new(Vector::new(0.0, 0.0, 1.0), T).unwrap(),
        Direction::new(Vector::new(1.0, 0.0, 0.0), T).unwrap(),
        T,
    )
    .unwrap()
}

/// The circle edges of `s` whose centres stand at height `z`.
fn circles_at(m: &Model, s: &Shape, z: f64) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            matches!(
                ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap().0,
                ogeom_geom::Curve::Circle(c) if (c.circle().centre().z - z).abs() < 1e-9
            )
        })
        .collect()
}

/// The edges of `s` both of whose ends `keep`.
fn edges_where(m: &Model, s: &Shape, keep: impl Fn(Point, Point) -> bool) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
            keep(point_of(m, &a), point_of(m, &b))
        })
        .collect()
}

/// Every vertex of `s` as bits, sorted, and its volume's bits.
fn bits(m: &Model, s: &Shape) -> (u64, Vec<[u64; 3]>) {
    let mut out: Vec<[u64; 3]> = explore_unique(m, s, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| {
            let p = point_of(m, v);
            [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
        })
        .collect();
    out.sort_unstable();
    (volume(m, s).to_bits(), out)
}

/// The guards every case is held to: validity, and every face proved
/// within `original`'s box grown by `fill` (a concave round's reach past
/// the original's faces, nought where the request fills nothing).
fn guards(m: &Model, original: &Shape, result: &Shape, fill: f64) {
    let found = ogeom_algo::check(m, result, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    let corners = ogeom_algo::shape_bounds(m, original, T).unwrap().corners();
    let (mut lo, mut hi) = (corners[0], corners[0]);
    for c in &corners {
        lo = Point::new(lo.x.min(c.x), lo.y.min(c.y), lo.z.min(c.z));
        hi = Point::new(hi.x.max(c.x), hi.y.max(c.y), hi.z.max(c.z));
    }
    let mut sides = Vec::new();
    for n in [
        Vector::new(1.0, 0.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
        Vector::new(0.0, 0.0, 1.0),
    ] {
        sides.push(HalfSpaces::Below {
            normal: n,
            at: hi + n * fill,
        });
        sides.push(HalfSpaces::Below {
            normal: -n,
            at: lo - n * fill,
        });
    }
    let bounds = HalfSpaces::Every(sides);
    for face in explore_unique(m, result, ShapeType::Face).unwrap() {
        let proved = ogeom_algo::certify_face_within(m, &face, &bounds, FIT, T).unwrap();
        assert!(proved.is_ok(), "a face stands outside the box: {proved:?}");
    }
}

/// `round` asked in the order given and reversed: the same solid to the
/// bit, held to the guards; the first.
fn every_order(
    m: &mut Model,
    original: &Shape,
    asked: &[(Shape, f64)],
    fill: f64,
) -> ogeom_algo::Built {
    let built = ogeom_fillet::fillet_edges_with(m, original, asked, T).unwrap();
    let reversed: Vec<(Shape, f64)> = asked.iter().rev().cloned().collect();
    let back = ogeom_fillet::fillet_edges_with(m, original, &reversed, T).unwrap();
    assert_eq!(
        bits(m, &built.shape),
        bits(m, &back.shape),
        "the request's order changes the solid"
    );
    guards(m, original, &built.shape, fill);
    built
}

/// [`every_order`] for chamfers of one distance.
fn every_order_chamfered(m: &mut Model, original: &Shape, edges: &[Shape], d: f64) -> Shape {
    let built = ogeom_fillet::chamfer_edges(m, original, edges, d, T).unwrap();
    let reversed: Vec<Shape> = edges.iter().rev().cloned().collect();
    let back = ogeom_fillet::chamfer_edges(m, original, &reversed, d, T).unwrap();
    assert_eq!(
        bits(m, &built.shape),
        bits(m, &back.shape),
        "the request's order changes the solid"
    );
    guards(m, original, &built.shape, 0.0);
    built.shape
}

/// Every blend face of `built` (a face the history makes of an asked edge)
/// meets the faces its edge was between tangent within [`FIT`] (those of
/// `original` beside the edge, as the history carries them), and every
/// blend face it touches. A round's end on a face across its edge (a flush
/// end on a cap) is a crease by design.
fn tangent(m: &Model, original: &Shape, built: &ogeom_algo::Built, asked: &[(Shape, f64)]) {
    let blends: Vec<Shape> = asked
        .iter()
        .flat_map(|(e, _)| built.history.generated(e).to_vec())
        .collect();
    for (edge, _) in asked {
        let mut beside: Vec<Shape> = ogeom_topo::ancestors_of(m, original, edge, ShapeType::Face)
            .unwrap()
            .iter()
            .flat_map(|f| {
                let mut images = built.history.modified(f).to_vec();
                images.push(f.clone());
                images
            })
            .collect();
        beside.extend(blends.iter().cloned());
        for face in built.history.generated(edge) {
            let contacts = ogeom_fillet::analyse_blend(m, &built.shape, face, 17, T).unwrap();
            for c in contacts
                .iter()
                .filter(|c| beside.iter().any(|f| f.is_same(&c.neighbour)))
            {
                assert!(c.tangency_error <= FIT, "a crease across a blend: {c:?}");
            }
        }
    }
}

/// The region a round of `r` takes off, or fills into, a corner whose two
/// faces leave it along the unit directions `a` and `b` of a meridian
/// (`(radial, axial)`), the region between them: its area, and its
/// centroid's offset from the corner, as `(radial, axial)`. The region is
/// the kite from the corner to the two contacts and the ball's centre, less
/// the ball's sector between the contacts.
fn corner_region(a: (f64, f64), b: (f64, f64), r: f64) -> (f64, (f64, f64)) {
    let opening = (a.0 * b.0 + a.1 * b.1).clamp(-1.0, 1.0).acos();
    let half = opening / 2.0;
    let leg = r / half.tan();
    let to_centre = r / half.sin();
    let bisector = {
        let (x, y) = (a.0 + b.0, a.1 + b.1);
        let n = x.hypot(y);
        (x / n, y / n)
    };
    let kite = leg * r;
    let kite_at = (leg * half.cos() + to_centre) / 3.0;
    let sweep = PI - opening;
    let sector = sweep * r * r / 2.0;
    let sector_at = to_centre - 4.0 * r * (sweep / 2.0).sin() / (3.0 * sweep);
    let area = kite - sector;
    let at = (kite * kite_at - sector * sector_at) / area;
    (area, (bisector.0 * at, bisector.1 * at))
}

/// The volume a rim's round sweeps at a corner `radius` off the axis,
/// [`corner_region`]'s area round its centroid: Pappus.
fn rim_swept(radius: f64, a: (f64, f64), b: (f64, f64), r: f64) -> f64 {
    let (area, (dr, _)) = corner_region(a, b, r);
    area * 2.0 * PI * (radius + dr)
}

/// The volume a chamfer of `d` along both faces sweeps at a rim `radius` off
/// the axis, its faces leaving along `a` and `b`: the triangle from the
/// corner to the two contacts, round its centroid.
fn rim_chamfered(radius: f64, a: (f64, f64), b: (f64, f64), d: f64) -> f64 {
    let cross = (a.0 * b.1 - a.1 * b.0).abs();
    let area = d * d * cross / 2.0;
    let dr = d * (a.0 + b.0) / 3.0;
    area * 2.0 * PI * (radius + dr)
}

/// Twice the signed area a closed boundary of segments and arcs encloses,
/// by Green's theorem: `x dy - y dx` summed over each piece.
enum Piece {
    Segment((f64, f64), (f64, f64)),
    /// An arc about a centre, of a radius, between two angles.
    Arc((f64, f64), f64, f64, f64),
}

fn twice_area(boundary: &[Piece]) -> f64 {
    boundary
        .iter()
        .map(|piece| match *piece {
            Piece::Segment(p, q) => p.0 * q.1 - q.0 * p.1,
            Piece::Arc(c, rho, t0, t1) => {
                rho * (c.0 * (t1.sin() - t0.sin()) - c.1 * (t1.cos() - t0.cos()))
                    + rho * rho * (t1 - t0)
            }
        })
        .sum()
}

/// The area a round of `r` takes off the corner where the line `x = 0`
/// meets the circle of radius `big` about the origin at `(0, big)`, the
/// material `x >= 0` inside the circle: the corner, down the line to the
/// ball's contact, back over the ball's top (its angle falling from the
/// line's side to the circle's) to its contact on the circle, and along
/// the circle to the corner.
fn line_circle_corner(big: f64, r: f64) -> f64 {
    let centre = (r, ((big - r) * (big - r) - r * r).sqrt());
    let on_line = (0.0, centre.1);
    let on_circle = (centre.0 * big / (big - r), centre.1 * big / (big - r));
    let angle = |c: (f64, f64), p: (f64, f64)| (p.1 - c.1).atan2(p.0 - c.0);
    let boundary = [
        Piece::Segment((0.0, big), on_line),
        Piece::Arc(centre, r, angle(centre, on_line), angle(centre, on_circle)),
        Piece::Arc((0.0, 0.0), big, angle((0.0, 0.0), on_circle), PI / 2.0),
    ];
    let twice = twice_area(&boundary);
    (twice / 2.0).abs()
}

/// The area a chamfer `d` along both faces takes off [`line_circle_corner`]'s
/// corner, the distance on the circle measured along it.
fn line_circle_chamfer(big: f64, d: f64) -> f64 {
    let on_line = (0.0, big - d);
    let turn = d / big;
    let on_circle = (big * turn.sin(), big * turn.cos());
    let boundary = [
        Piece::Segment((0.0, big), on_line),
        Piece::Segment(on_line, on_circle),
        Piece::Arc((0.0, 0.0), big, PI / 2.0 - turn, PI / 2.0),
    ];
    (twice_area(&boundary) / 2.0).abs()
}

// Directions leaving a rim's corner in the meridian, `(radial, axial)`.
const IN: (f64, f64) = (-1.0, 0.0);
const OUT: (f64, f64) = (1.0, 0.0);
const DOWN: (f64, f64) = (0.0, -1.0);
const UP: (f64, f64) = (0.0, 1.0);

/// A drum of radius 2 only 1 tall, both rims rounded at 0.49, the walls'
/// two rings 0.02 apart: each rim takes its Pappus ring. On main: built,
/// the volume right, but the two orders give different solids.
#[test]
#[ignore = "heavy"]
fn a_short_drum_s_rims_take_their_rings_whatever_the_order() {
    if !asked() {
        return;
    }
    let mut m = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 2.0, 1.0, T)
        .unwrap()
        .shape;
    let mut asked: Vec<(Shape, f64)> = circles_at(&m, &drum, 1.0)
        .into_iter()
        .map(|e| (e, 0.49))
        .collect();
    asked.extend(circles_at(&m, &drum, 0.0).into_iter().map(|e| (e, 0.49)));
    let built = every_order(&mut m, &drum, &asked, 0.0);
    let want = PI * 4.0 - rim_swept(2.0, IN, DOWN, 0.49) - rim_swept(2.0, IN, UP, 0.49);
    let v = volume(&m, &built.shape);
    assert!((v - want).abs() < 1e-6, "{v} against {want}");
    tangent(&m, &drum, &built, &asked);
}

/// A tube of radius 3 bored to 2, 4 tall, all four rims at 0.3: the outer
/// two inside the wall, the bore's two outside it. On main: built, the
/// volume right, the two orders different.
#[test]
#[ignore = "heavy"]
fn a_tube_s_four_rims_take_their_rings_whatever_the_order() {
    if !asked() {
        return;
    }
    let mut m = Model::new();
    let outer = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 4.0, T)
        .unwrap()
        .shape;
    let inner = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, -1.0), 2.0, 6.0, T)
        .unwrap()
        .shape;
    let tube = ogeom_bool::cut(&mut m, &outer, &inner, T).unwrap().shape;
    let mut asked: Vec<(Shape, f64)> = circles_at(&m, &tube, 4.0)
        .into_iter()
        .map(|e| (e, 0.3))
        .collect();
    asked.extend(circles_at(&m, &tube, 0.0).into_iter().map(|e| (e, 0.3)));
    let built = every_order(&mut m, &tube, &asked, 0.0);
    let want = PI * 5.0 * 4.0
        - rim_swept(3.0, IN, DOWN, 0.3)
        - rim_swept(3.0, IN, UP, 0.3)
        - rim_swept(2.0, OUT, DOWN, 0.3)
        - rim_swept(2.0, OUT, UP, 0.3);
    let v = volume(&m, &built.shape);
    assert!((v - want).abs() < 1e-6, "{v} against {want}");
    tangent(&m, &tube, &built, &asked);
}

/// A conical tube, a frustum narrowing from radius 3 to 2 over 4 bored to
/// 1, all four rims at 0.3: the cone's rims meet their caps off square.
/// On main: built, the two orders different.
#[test]
#[ignore = "heavy"]
fn a_conical_tube_s_rims_meet_their_caps_off_square_whatever_the_order() {
    if !asked() {
        return;
    }
    let mut m = Model::new();
    let outer = ogeom_algo::make_cone(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 2.0, 4.0, T)
        .unwrap()
        .shape;
    let inner = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, -1.0), 1.0, 6.0, T)
        .unwrap()
        .shape;
    let cone = ogeom_bool::cut(&mut m, &outer, &inner, T).unwrap().shape;
    let mut asked: Vec<(Shape, f64)> = circles_at(&m, &cone, 4.0)
        .into_iter()
        .map(|e| (e, 0.3))
        .collect();
    asked.extend(circles_at(&m, &cone, 0.0).into_iter().map(|e| (e, 0.3)));
    let built = every_order(&mut m, &cone, &asked, 0.0);
    // The wall leaves the top rim down and out, the bottom rim up and in.
    let slope = {
        let n = 1.0f64.hypot(4.0);
        (1.0 / n, 4.0 / n)
    };
    let frustum = PI * 4.0 / 3.0 * (9.0 + 6.0 + 4.0);
    let want = frustum
        - PI * 4.0
        - rim_swept(2.0, IN, (slope.0, -slope.1), 0.3)
        - rim_swept(3.0, IN, (-slope.0, slope.1), 0.3)
        - rim_swept(1.0, OUT, DOWN, 0.3)
        - rim_swept(1.0, OUT, UP, 0.3);
    let v = volume(&m, &built.shape);
    assert!((v - want).abs() < 1e-6, "{v} against {want}");
    tangent(&m, &cone, &built, &asked);
}

/// A block with a countersunk bore: radius 1.5 through it, a cone opening
/// to 3.5 at the top over 2. Its three rims at 0.3: the top one on the
/// cone, off square; the one where the cone meets the bore, between two
/// curved walls; the bore's bottom one. On main: refused in one order ("a
/// face boundary strand dangles"), built in the other.
#[test]
#[ignore = "heavy"]
fn a_countersink_s_rims_round_whatever_the_order() {
    if !asked() {
        return;
    }
    let mut m = Model::new();
    let block = ogeom_algo::make_box(&mut m, Frame::WORLD, (10.0, 8.0, 4.0), T)
        .unwrap()
        .shape;
    let bore = ogeom_algo::make_cylinder(&mut m, up_at(5.0, 4.0, -1.0), 1.5, 6.0, T)
        .unwrap()
        .shape;
    let sink = ogeom_algo::make_cone(&mut m, up_at(5.0, 4.0, 2.0), 1.5, 4.5, 3.0, T)
        .unwrap()
        .shape;
    let holed = ogeom_bool::cut(&mut m, &block, &bore, T).unwrap().shape;
    let sunk = ogeom_bool::cut(&mut m, &holed, &sink, T).unwrap().shape;
    let before = volume(&m, &sunk);
    let mut asked: Vec<(Shape, f64)> = circles_at(&m, &sunk, 4.0)
        .into_iter()
        .map(|e| (e, 0.3))
        .collect();
    asked.extend(circles_at(&m, &sunk, 2.0).into_iter().map(|e| (e, 0.3)));
    asked.extend(circles_at(&m, &sunk, 0.0).into_iter().map(|e| (e, 0.3)));
    let built = every_order(&mut m, &sunk, &asked, 0.0);
    // The cone leaves the top rim down and in at 45 degrees, and the rim
    // where it meets the bore up and out, the bore down: two walls, 135
    // degrees apart through the material.
    let s = core::f64::consts::FRAC_1_SQRT_2;
    let want = before
        - rim_swept(3.5, OUT, (-s, -s), 0.3)
        - rim_swept(1.5, (s, s), DOWN, 0.3)
        - rim_swept(1.5, OUT, UP, 0.3);
    let v = volume(&m, &built.shape);
    assert!((v - want).abs() < 1e-6, "{v} against {want}");
    tangent(&m, &sunk, &built, &asked);
}

/// Half a drum (a D-shaped block, radius 3, 10 long), its two straight
/// edges between the flat and the curved face rounded at 0.5: the half
/// disc less the line-and-circle corner twice, along the 10. On main:
/// built, the two orders different.
#[test]
#[ignore = "heavy"]
fn a_d_block_s_straight_edges_round_along_its_curved_face() {
    if !asked() {
        return;
    }
    let mut m = Model::new();
    let full = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 10.0, T)
        .unwrap()
        .shape;
    let half = ogeom_algo::make_box(&mut m, up_at(-4.0, -4.0, -1.0), (4.0, 8.0, 12.0), T)
        .unwrap()
        .shape;
    let d = ogeom_bool::cut(&mut m, &full, &half, T).unwrap().shape;
    let asked: Vec<(Shape, f64)> = edges_where(&m, &d, |p, q| {
        p.x.abs() < 1e-9 && q.x.abs() < 1e-9 && (p.z - q.z).abs() > 1.0
    })
    .into_iter()
    .map(|e| (e, 0.5))
    .collect();
    assert_eq!(asked.len(), 2);
    let built = every_order(&mut m, &d, &asked, 0.0);
    let want = (PI * 9.0 / 2.0 - 2.0 * line_circle_corner(3.0, 0.5)) * 10.0;
    let v = volume(&m, &built.shape);
    assert!((v - want).abs() < 1e-6, "{v} against {want}");
    tangent(&m, &d, &built, &asked);
}

/// A stadium prism (a 10 by 4 box between two drums of radius 2, 3 tall),
/// its top loop of two straight edges and two half rims rounded at 0.5 as
/// one tangent chain: the corner's area along the straight edges and round
/// each half rim by Pappus. On main: refused ("the edge's faces are
/// tangent").
#[test]
#[ignore = "heavy"]
fn a_stadium_s_top_loop_rounds_as_one_chain() {
    if !asked() {
        return;
    }
    let mut m = Model::new();
    let st = stadium(&mut m);
    let asked: Vec<(Shape, f64)> = edges_where(&m, &st, |p, q| {
        (p.z - 3.0).abs() < 1e-9 && (q.z - 3.0).abs() < 1e-9
    })
    .into_iter()
    .map(|e| (e, 0.5))
    .collect();
    let built = every_order(&mut m, &st, &asked, 0.0);
    let (area, _) = corner_region(IN, DOWN, 0.5);
    let body = 10.0 * 4.0 * 3.0 + PI * 4.0 * 3.0;
    let want = body - 2.0 * area * 10.0 - 2.0 * rim_swept(2.0, IN, DOWN, 0.5) / 2.0;
    let v = volume(&m, &built.shape);
    assert!((v - want).abs() < 1e-6, "{v} against {want}");
    tangent(&m, &st, &built, &asked);
}

/// The stadium prism of [`a_stadium_s_top_loop_rounds_as_one_chain`].
fn stadium(m: &mut Model) -> Shape {
    let core = ogeom_algo::make_box(m, up_at(-5.0, -2.0, 0.0), (10.0, 4.0, 3.0), T)
        .unwrap()
        .shape;
    let left = ogeom_algo::make_cylinder(m, up_at(-5.0, 0.0, 0.0), 2.0, 3.0, T)
        .unwrap()
        .shape;
    let right = ogeom_algo::make_cylinder(m, up_at(5.0, 0.0, 0.0), 2.0, 3.0, T)
        .unwrap()
        .shape;
    let st = ogeom_bool::fuse(m, &core, &left, T).unwrap().shape;
    ogeom_bool::fuse(m, &st, &right, T).unwrap().shape
}

/// A 40 by 40 by 5 plate with a boss of radius 6 and 5 tall in its middle,
/// every edge rounded at 0.5: the plate as the rounded box, the boss's top
/// rim taking its ring and its foot filling outside it. On main: refused
/// ("a marched fillet needs an edge shared by exactly two faces").
#[test]
#[ignore = "heavy"]
fn a_bossed_plate_rounds_on_every_edge() {
    if !asked() {
        return;
    }
    let mut m = Model::new();
    let plate = ogeom_algo::make_box(&mut m, up_at(-20.0, -20.0, 0.0), (40.0, 40.0, 5.0), T)
        .unwrap()
        .shape;
    let boss = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 5.0), 6.0, 5.0, T)
        .unwrap()
        .shape;
    let bossed = ogeom_bool::fuse(&mut m, &plate, &boss, T).unwrap().shape;
    let r = 0.5;
    let asked: Vec<(Shape, f64)> = explore_unique(&m, &bossed, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .map(|e| (e, r))
        .collect();
    let built = every_order(&mut m, &bossed, &asked, r);
    let (a, b, c) = (40.0 - 2.0 * r, 40.0 - 2.0 * r, 5.0 - 2.0 * r);
    let rounded_plate = a * b * c
        + 2.0 * r * (a * b + b * c + c * a)
        + PI * r * r * (a + b + c)
        + 4.0 / 3.0 * PI * r * r * r;
    let want =
        rounded_plate + PI * 36.0 * 5.0 - rim_swept(6.0, IN, DOWN, r) + rim_swept(6.0, OUT, UP, r);
    let v = volume(&m, &built.shape);
    assert!((v - want).abs() < 1e-6, "{v} against {want}");
    tangent(&m, &bossed, &built, &asked);
}

/// The bossed plate of [`a_bossed_plate_rounds_on_every_edge`], its four
/// upright edges rounded at 1 and every other edge at 0.5: no one ball fits
/// a plate corner, so each closes on a vertex patch. The plate's part of
/// the solid is the plain box rounded alike through the planar engine (the
/// boss stands far from its corners), the boss's by Pappus; every patch
/// meets its rounds tangent within the fit target and stays within the
/// faces round its vertex. On main: refused ("a marched fillet needs an
/// edge shared by exactly two faces").
#[test]
#[ignore = "heavy"]
fn a_bossed_plate_s_corners_close_where_no_ball_fits() {
    if !asked() {
        return;
    }
    let mut m = Model::new();
    let plate = ogeom_algo::make_box(&mut m, up_at(-20.0, -20.0, 0.0), (40.0, 40.0, 5.0), T)
        .unwrap()
        .shape;
    let boss = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 5.0), 6.0, 5.0, T)
        .unwrap()
        .shape;
    let bossed = ogeom_bool::fuse(&mut m, &plate, &boss, T).unwrap().shape;
    let (r, upright) = (0.5, 1.0);
    let radius_of = |m: &Model, e: &Shape| {
        let ends = ogeom_topo::explore_unique(m, e, ShapeType::Vertex).unwrap();
        let (p, q) = (point_of(m, &ends[0]), point_of(m, &ends[ends.len() - 1]));
        if ends.len() == 2 && (p.z - q.z).abs() > 1.0 {
            upright
        } else {
            r
        }
    };
    let asked: Vec<(Shape, f64)> = explore_unique(&m, &bossed, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .map(|e| {
            let k = radius_of(&m, &e);
            (e, k)
        })
        .collect();
    let built = every_order(&mut m, &bossed, &asked, r);
    let plain: Vec<(Shape, f64)> = explore_unique(&m, &plate, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .map(|e| {
            let k = radius_of(&m, &e);
            (e, k)
        })
        .collect();
    let box_rounded = ogeom_fillet::fillet_edges_with(&mut m, &plate, &plain, T).unwrap();
    let want = volume(&m, &box_rounded.shape) + PI * 36.0 * 5.0 - rim_swept(6.0, IN, DOWN, r)
        + rim_swept(6.0, OUT, UP, r);
    let v = volume(&m, &built.shape);
    assert!((v - want).abs() < 1e-6, "{v} against {want}");
    tangent(&m, &bossed, &built, &asked);
    let mut patches = 0;
    for corner in explore_unique(&m, &bossed, ShapeType::Vertex).unwrap() {
        for face in built.history.generated(&corner) {
            patches += 1;
            let contacts = ogeom_fillet::analyse_blend(&m, &built.shape, face, 17, T).unwrap();
            for c in &contacts {
                assert!(c.tangency_error <= FIT, "a crease round a patch: {c:?}");
            }
        }
    }
    assert!(patches >= 8, "{patches} patches");
}

/// A plate of 40 by 40 by 5, a boss of radius 6 at (14, 0)
/// tangent to the side x = 20, the foot and the top edge along x = 20
/// rounded as one chain, the S-blend at the cusp where they meet. No closed
/// form: every blend tangent within the fit target, in either order, at
/// 0.5 and at 1.25. On main: refused at both ("the filling misses its
/// tolerance", 2.1e-4 and 3.0e-4 rad).
#[test]
#[ignore = "heavy"]
fn the_cusp_s_blend_closes_tangent_at_radii_other_than_one() {
    if !asked() {
        return;
    }
    for r in [0.5, 1.25] {
        let mut m = Model::new();
        let plate = ogeom_algo::make_box(&mut m, up_at(-20.0, -20.0, 0.0), (40.0, 40.0, 5.0), T)
            .unwrap()
            .shape;
        let boss = ogeom_algo::make_cylinder(&mut m, up_at(14.0, 0.0, 5.0), 6.0, 5.0, T)
            .unwrap()
            .shape;
        let part = ogeom_bool::fuse(&mut m, &plate, &boss, T).unwrap().shape;
        let mut asked: Vec<(Shape, f64)> = circles_at(&m, &part, 5.0)
            .into_iter()
            .map(|e| (e, r))
            .collect();
        asked.extend(
            edges_where(&m, &part, |p, q| {
                (p.x - 20.0).abs() < 1e-9
                    && (q.x - 20.0).abs() < 1e-9
                    && (p.z - 5.0).abs() < 1e-9
                    && (q.z - 5.0).abs() < 1e-9
            })
            .into_iter()
            .map(|e| (e, r)),
        );
        let built = every_order(&mut m, &part, &asked, r);
        tangent(&m, &part, &built, &asked);
    }
}

/// The tube of [`a_tube_s_four_rims_take_their_rings_whatever_the_order`],
/// its four rims chamfered at 0.3: each rim's triangle round its centroid.
/// On main: built, the two orders different.
#[test]
#[ignore = "heavy"]
fn a_tube_s_rims_chamfer_whatever_the_order() {
    if !asked() {
        return;
    }
    let mut m = Model::new();
    let outer = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 4.0, T)
        .unwrap()
        .shape;
    let inner = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, -1.0), 2.0, 6.0, T)
        .unwrap()
        .shape;
    let tube = ogeom_bool::cut(&mut m, &outer, &inner, T).unwrap().shape;
    let mut edges = circles_at(&m, &tube, 4.0);
    edges.extend(circles_at(&m, &tube, 0.0));
    let built = every_order_chamfered(&mut m, &tube, &edges, 0.3);
    let want = PI * 5.0 * 4.0
        - rim_chamfered(3.0, IN, DOWN, 0.3)
        - rim_chamfered(3.0, IN, UP, 0.3)
        - rim_chamfered(2.0, OUT, DOWN, 0.3)
        - rim_chamfered(2.0, OUT, UP, 0.3);
    let v = volume(&m, &built);
    assert!((v - want).abs() < 1e-6, "{v} against {want}");
}

/// The stadium's top loop chamfered at 0.5: the triangle along the straight
/// edges and round each half rim. On main: refused.
#[test]
#[ignore = "heavy"]
fn a_stadium_s_top_loop_chamfers_as_one_chain() {
    if !asked() {
        return;
    }
    let mut m = Model::new();
    let st = stadium(&mut m);
    let edges = edges_where(&m, &st, |p, q| {
        (p.z - 3.0).abs() < 1e-9 && (q.z - 3.0).abs() < 1e-9
    });
    let built = every_order_chamfered(&mut m, &st, &edges, 0.5);
    let body = 10.0 * 4.0 * 3.0 + PI * 4.0 * 3.0;
    let want = body - 2.0 * 0.125 * 10.0 - 2.0 * rim_chamfered(2.0, IN, DOWN, 0.5) / 2.0;
    let v = volume(&m, &built);
    assert!((v - want).abs() < 1e-6, "{v} against {want}");
}

/// The D-shaped block's straight edges chamfered at 0.5, the distance on
/// the curved face measured along it: the half disc less the corner
/// between the line, the chord and the arc, twice, along the 10. On main:
/// refused ("blending an edge of a curved face needs the marching blend
/// machinery").
#[test]
#[ignore = "heavy"]
fn a_d_block_s_straight_edges_chamfer_along_its_curved_face() {
    if !asked() {
        return;
    }
    let mut m = Model::new();
    let full = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 10.0, T)
        .unwrap()
        .shape;
    let half = ogeom_algo::make_box(&mut m, up_at(-4.0, -4.0, -1.0), (4.0, 8.0, 12.0), T)
        .unwrap()
        .shape;
    let d = ogeom_bool::cut(&mut m, &full, &half, T).unwrap().shape;
    let edges = edges_where(&m, &d, |p, q| {
        p.x.abs() < 1e-9 && q.x.abs() < 1e-9 && (p.z - q.z).abs() > 1.0
    });
    let built = every_order_chamfered(&mut m, &d, &edges, 0.5);
    let want = (PI * 9.0 / 2.0 - 2.0 * line_circle_chamfer(3.0, 0.5)) * 10.0;
    let v = volume(&m, &built);
    assert!((v - want).abs() < 1e-6, "{v} against {want}");
}

/// A drilled flange: a disc of radius 20, 4 thick, four bolt holes of
/// radius 2 on a circle of 14, every edge rounded at 0.5. Each rim takes
/// its ring about its own axis. On main: refused ("a marched fillet needs
/// an edge shared by exactly two faces").
#[test]
#[ignore = "heavy"]
fn a_drilled_flange_rounds_on_every_edge() {
    if !asked() {
        return;
    }
    let mut m = Model::new();
    let mut flange = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 20.0, 4.0, T)
        .unwrap()
        .shape;
    for k in 0..4 {
        let at = core::f64::consts::FRAC_PI_2 * f64::from(k);
        let hole = ogeom_algo::make_cylinder(
            &mut m,
            up_at(14.0 * at.cos(), 14.0 * at.sin(), -1.0),
            2.0,
            6.0,
            T,
        )
        .unwrap()
        .shape;
        flange = ogeom_bool::cut(&mut m, &flange, &hole, T).unwrap().shape;
    }
    let r = 0.5;
    let asked: Vec<(Shape, f64)> = explore_unique(&m, &flange, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            matches!(
                ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap().0,
                ogeom_geom::Curve::Circle(_)
            )
        })
        .map(|e| (e, r))
        .collect();
    let built = every_order(&mut m, &flange, &asked, 0.0);
    let body = PI * 400.0 * 4.0 - 4.0 * PI * 4.0 * 4.0;
    let want = body
        - rim_swept(20.0, IN, DOWN, r)
        - rim_swept(20.0, IN, UP, r)
        - 4.0 * (rim_swept(2.0, OUT, DOWN, r) + rim_swept(2.0, OUT, UP, r));
    let v = volume(&m, &built.shape);
    assert!((v - want).abs() < 1e-6, "{v} against {want}");
    tangent(&m, &flange, &built, &asked);
}

/// The oracles' own closed forms against cases worked by hand: a square
/// corner's round is `(1 - pi/4) r^2` with its centroid `r (10 - 3 pi) /
/// (12 - 3 pi)` along each leg; a chamfer's triangle is `d^2 / 2`; a round
/// of a line and a circle matches the region counted on a fine grid; and
/// the chamfer of a line and a circle is its triangle and the circle's
/// segment beyond its chord.
#[test]
fn the_oracles_agree_with_cases_worked_by_hand() {
    let r = 0.7;
    let (area, (dr, dz)) = corner_region(IN, DOWN, r);
    assert!((area - (1.0 - PI / 4.0) * r * r).abs() < 1e-14);
    let x = r * (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI);
    assert!(
        (dr + x).abs() < 1e-14 && (dz + x).abs() < 1e-14,
        "{dr} {dz} against {x}"
    );
    let tri = rim_chamfered(1.0, IN, DOWN, 0.3) / (2.0 * PI * (1.0 - 0.1));
    assert!((tri - 0.045).abs() < 1e-14);
    // The line-and-circle corner against the region counted on a grid:
    // the points of the half disc beyond the ball, in the wedge between the
    // ball's two contacts.
    let (big, r): (f64, f64) = (3.0, 0.5);
    let centre = (r, ((big - r) * (big - r) - r * r).sqrt());
    let to_circle = centre.1.atan2(centre.0);
    let n = 4000;
    let (w, h) = (centre.0 * big / (big - r), big - centre.1);
    let cell = (w / f64::from(n)) * (h / f64::from(n));
    let mut counted = 0.0;
    for i in 0..n {
        for j in 0..n {
            let x = (f64::from(i) + 0.5) * w / f64::from(n);
            let y = centre.1 + (f64::from(j) + 0.5) * h / f64::from(n);
            let (dx, dy) = (x - centre.0, y - centre.1);
            let angle = dy.atan2(dx);
            if x * x + y * y <= big * big && dx.hypot(dy) > r && angle > to_circle && angle <= PI {
                counted += cell;
            }
        }
    }
    let exact = line_circle_corner(big, r);
    assert!((exact - counted).abs() < 1e-5, "{exact} against {counted}");
    // The chamfer's: the triangle and the circle's segment beyond its chord.
    let turn = 0.5 / big;
    let segment = big * big * (turn - turn.sin()) / 2.0;
    let triangle = {
        let (p, q) = ((0.0, big - 0.5), (big * turn.sin(), big * turn.cos()));
        ((0.0 - q.0) * (p.1 - q.1) - (p.0 - q.0) * (big - q.1)).abs() / 2.0
    };
    let flat = line_circle_chamfer(big, 0.5);
    assert!((flat - (triangle + segment)).abs() < 1e-14, "{flat}");
}
