//! Non-local rounds the strips' own tools do not close, through the general
//! route of the blend engine's one boolean (ADR 0001, M5c): the far feature
//! taken off a copy of the part, the rounds built on the copy as the local
//! request they are there, and the copy less its rounding cut, the rounding
//! less the copy fused, on the part as given.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{FRAC_PI_4, PI};

use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn block(m: &mut Model, o: (f64, f64, f64), size: (f64, f64, f64)) -> Shape {
    ogeom_algo::make_box(m, at(o.0, o.1, o.2), size, T)
        .unwrap()
        .shape
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// The edges of `s` with an end at `p`.
fn edges_at(m: &Model, s: &Shape, p: Point) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
            let q = |v: &Shape| m.node(v).unwrap().data().as_vertex().unwrap().point;
            q(&a).distance(p) < 1e-6 || q(&b).distance(p) < 1e-6
        })
        .collect()
}

/// The area under `sqrt(r^2 - u^2)` from nought to `u`.
fn under(r: f64, u: f64) -> f64 {
    0.5 * (u * (r * r - u * u).sqrt() + r * r * (u / r).asin())
}

/// Rounded through the engine's one boolean: valid, every edge generating
/// its strip; the volume.
fn rounded(m: &mut Model, s: &Shape, edges: &[Shape], r: f64) -> f64 {
    let (built, work) =
        ogeom_bool::work::counted(|| ogeom_fillet::fillet_edges(m, s, edges, r, T).unwrap());
    assert!(work.face_pairs > 0, "through the booleans");
    let found = ogeom_algo::check(m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    for e in edges {
        assert!(!built.history.generated(e).is_empty(), "each strip");
    }
    volume(m, &built.shape)
}

/// A 20 by 20 by 10 block with a pocket 16 by 16 by 8 sunk 7 into its top,
/// a pillar 4 by 0.4 by 3 standing on the pocket's floor 0.3 in front of
/// its wall `y = 2`; the three concave edges at the pocket's corner (2, 2,
/// 3) rounded at 1, closing on the corner's sphere, a vertex not convex all
/// round. The fill along the wall's foot reaches into the pillar. The block
/// rounded alike with no pillar (the engine's own corner), with the
/// pillar, less what the fill and the pillar share (its height `1 - sqrt(1
/// - u^2)` at `u = 3 - y` behind the ball's centre, across the pillar).
#[test]
fn a_concave_corner_s_rounds_reaching_a_pillar_go_through_the_general_route() {
    let pocketed = |m: &mut Model| {
        let solid = block(m, (0.0, 0.0, 0.0), (20.0, 20.0, 10.0));
        let pocket = block(m, (2.0, 2.0, 3.0), (16.0, 16.0, 8.0));
        ogeom_bool::cut(m, &solid, &pocket, T).unwrap().shape
    };
    let corner = Point::new(2.0, 2.0, 3.0);
    let mut plain = Model::new();
    let p = pocketed(&mut plain);
    let edges = edges_at(&plain, &p, corner);
    assert_eq!(edges.len(), 3);
    let (built, work) = ogeom_bool::work::counted(|| {
        ogeom_fillet::fillet_edges(&mut plain, &p, &edges, 1.0, T).unwrap()
    });
    assert_eq!(work.face_pairs, 0, "the corner alone, the engine's own");
    let alone = volume(&plain, &built.shape);

    let mut m = Model::new();
    let p = pocketed(&mut m);
    let pillar = block(&mut m, (8.0, 2.3, 3.0), (4.0, 0.4, 3.0));
    let s = ogeom_bool::fuse(&mut m, &p, &pillar, T).unwrap().shape;
    let edges = edges_at(&m, &s, corner);
    assert_eq!(edges.len(), 3);
    let shared = 4.0 * (0.4 - (under(1.0, 0.7) - under(1.0, 0.3)));
    let want = alone + 4.0 * 0.4 * 3.0 - shared;
    let got = rounded(&mut m, &s, &edges, 1.0);
    assert!((got - want).abs() < 1e-6, "{got} against {want}");
}

/// A 40 by 40 by 5 plate with a boss of radius 6 on its middle (a curved
/// solid), the three edges at its corner (20, -20, 5) rounded at 1 on the
/// corner's sphere, a hole 0.5 by 0.07 drilled 0.3 into the top within the
/// front edge's spandrel. The strips meet at an angle, which no strip's
/// tool closes. The hole goes with the spandrel: the plate less each edge's
/// section `1 - pi/4` short of the corner's cube and the cube less the
/// sphere's octant, the boss whole.
#[test]
fn a_curved_part_s_corner_rounds_past_a_hole_through_the_general_route() {
    let mut m = Model::new();
    let plate = block(&mut m, (-20.0, -20.0, 0.0), (40.0, 40.0, 5.0));
    let boss = ogeom_algo::make_cylinder(&mut m, at(0.0, 0.0, 5.0), 6.0, 5.0, T)
        .unwrap()
        .shape;
    let bossed = ogeom_bool::fuse(&mut m, &plate, &boss, T).unwrap().shape;
    let hole = block(&mut m, (0.0, -19.95, 4.7), (0.5, 0.07, 1.0));
    let s = ogeom_bool::cut(&mut m, &bossed, &hole, T).unwrap().shape;
    let edges = edges_at(&m, &s, Point::new(20.0, -20.0, 5.0));
    assert_eq!(edges.len(), 3);
    let section = 1.0 - FRAC_PI_4;
    let taken = section * (40.0 + 40.0 + 5.0 - 3.0) + (1.0 - PI / 6.0);
    let want = 8000.0 + PI * 36.0 * 5.0 - taken;
    let got = rounded(&mut m, &s, &edges, 1.0);
    assert!((got - want).abs() < 1e-6, "{got} against {want}");
}

/// The half frustum of `cone_ruling.rs` (radius 6 down to 3 over 6, cut
/// through its axis), one ruling rounded at 1/2 where a hole 0.04 square is
/// drilled 1 into the cut face within its spandrel. A strip along a cone,
/// whose end faces stand askew to it, no strip's tool closes. The hole goes
/// with the spandrel: the half frustum rounded alike with no hole.
#[test]
fn a_cone_s_ruling_rounds_past_a_hole_through_the_general_route() {
    let half = |m: &mut Model| {
        let cone = ogeom_algo::make_cone(m, Frame::WORLD, 6.0, 3.0, 6.0, T)
            .unwrap()
            .shape;
        let tool = block(m, (-10.0, 0.0, -1.0), (20.0, 10.0, 8.0));
        ogeom_bool::cut(m, &cone, &tool, T).unwrap().shape
    };
    let ruling = |m: &Model, s: &Shape| {
        explore_unique(m, s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let (a, b) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
                let q = |v: &Shape| m.node(v).unwrap().data().as_vertex().unwrap().point;
                let (a, b) = (q(&a), q(&b));
                a.y.abs() < 1e-9 && b.y.abs() < 1e-9 && a.x > 1.0 && b.x > 1.0
            })
            .expect("the ruling at x > 0")
    };
    let mut plain = Model::new();
    let s = half(&mut plain);
    let e = ruling(&plain, &s);
    let built = ogeom_fillet::fillet_edges(&mut plain, &s, &[e], 0.5, T).unwrap();
    let want = volume(&plain, &built.shape);

    let mut m = Model::new();
    let s = half(&mut m);
    // Within the spandrel at height 3, where the ruling stands at x = 4.5:
    // just in from the edge on the cut face, under 0.1 deep.
    let hole = block(&mut m, (4.38, -0.06, 2.9), (0.04, 0.07, 0.2));
    let s = ogeom_bool::cut(&mut m, &s, &hole, T).unwrap().shape;
    let e = ruling(&m, &s);
    let got = rounded(&mut m, &s, &[e], 0.5);
    assert!((got - want).abs() < 1e-6, "{got} against {want}");
}
