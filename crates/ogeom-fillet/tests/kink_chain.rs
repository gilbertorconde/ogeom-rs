//! Tangent chains (ADR 0001, amendment "tangent chains stop at a kink"):
//! which edges a pick runs on to, and what a round ending at a kink
//! takes.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom_core::Tolerances;
use ogeom_fillet::{TANGENT_CHAIN_ANGLE, fillet_edges, tangent_chain};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// The prism of `profile` in the plane z = 0, `height` up.
fn prism(m: &mut Model, profile: &[(f64, f64)], height: f64) -> Shape {
    let points: Vec<Point> = profile
        .iter()
        .map(|&(x, y)| Point::new(x, y, 0.0))
        .collect();
    let wire = ogeom_algo::make_polygon(m, &points, true, T).unwrap().shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-1000.0, 1000.0),
        (-1000.0, 1000.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(m, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let swept = ogeom_algo::make_prism(m, &face, Vector::new(0.0, 0.0, height), T)
        .unwrap()
        .shape;
    ogeom_algo::baked_shape(m, &swept, T).unwrap().shape
}

/// The edge of `s` whose middle stands nearest `p`.
fn edge_near(m: &Model, s: &Shape, p: Point) -> Shape {
    let middle = |e: &Shape| {
        let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
        c.point_at(0.5 * (a + b), T).unwrap()
    };
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| middle(a).distance(p).total_cmp(&middle(b).distance(p)))
        .unwrap()
}

/// A block 10 high whose front y = 0 runs from x = 0 to 20 and then turns
/// `kink` radians into the material, on to x = `far`; its back is y = 10.
/// With the front's two top edges.
fn kinked(m: &mut Model, kink: f64, far: f64) -> (Shape, Shape, Shape) {
    let s = prism(
        m,
        &[
            (0.0, 0.0),
            (20.0, 0.0),
            (far, (far - 20.0) * kink.tan()),
            (far, 10.0),
            (0.0, 10.0),
        ],
        10.0,
    );
    let first = edge_near(m, &s, Point::new(10.0, 0.0, 10.0));
    let second = edge_near(
        m,
        &s,
        Point::new(0.5 * (20.0 + far), 0.5 * (far - 20.0) * kink.tan(), 10.0),
    );
    (s, first, second)
}

fn same(a: &[Shape], b: &[Shape]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.is_same(y))
}

/// A 4 degree kink stops the chain at the default angle, as it stops it in
/// every mainstream system's documented tangency tolerance (one degree at
/// most); a 5 degree angle carries it on, and a half-degree kink is within
/// the default. The front's top edges meet the uprights square, so nothing
/// else joins.
#[test]
fn a_kink_wider_than_the_angle_stops_the_chain() {
    let mut m = Model::new();
    let (s, first, second) = kinked(&mut m, 4.0_f64.to_radians(), 40.0);
    let chain =
        tangent_chain(&m, &s, std::slice::from_ref(&first), TANGENT_CHAIN_ANGLE, T).unwrap();
    assert!(
        same(&chain, std::slice::from_ref(&first)),
        "{} edges",
        chain.len()
    );
    let wide = tangent_chain(
        &m,
        &s,
        std::slice::from_ref(&first),
        5.0_f64.to_radians(),
        T,
    )
    .unwrap();
    assert!(same(&wide, &[first, second]), "{} edges", wide.len());

    let mut m = Model::new();
    let (s, first, second) = kinked(&mut m, 0.5_f64.to_radians(), 40.0);
    let chain =
        tangent_chain(&m, &s, std::slice::from_ref(&first), TANGENT_CHAIN_ANGLE, T).unwrap();
    assert!(same(&chain, &[first, second]), "{} edges", chain.len());
}

/// A block 30 × 20 × 10 with its four uprights rounded 4: one straight
/// edge of its top brings the whole top loop, four sides and four arcs,
/// and nothing of the bottom or the uprights' seams.
#[test]
fn a_rounded_outline_s_side_brings_its_whole_loop() {
    let mut m = Model::new();
    let block = ogeom_algo::make_box(
        &mut m,
        Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap(),
        (30.0, 20.0, 10.0),
        T,
    )
    .unwrap()
    .shape;
    let uprights: Vec<Shape> = [(0.0, 0.0), (30.0, 0.0), (30.0, 20.0), (0.0, 20.0)]
        .iter()
        .map(|&(x, y)| edge_near(&m, &block, Point::new(x, y, 5.0)))
        .collect();
    let rounded = fillet_edges(&mut m, &block, &uprights, 4.0, T)
        .unwrap()
        .shape;
    let side = edge_near(&m, &rounded, Point::new(15.0, 0.0, 10.0));
    let chain = tangent_chain(&m, &rounded, &[side], TANGENT_CHAIN_ANGLE, T).unwrap();
    assert_eq!(chain.len(), 8, "the top loop");
    for e in &chain {
        let bb = ogeom_algo::shape_bounds(&m, e, T).unwrap();
        let (l, h) = (bb.low().unwrap(), bb.high().unwrap());
        assert!(
            (l.z - 10.0).abs() < 1e-6 && (h.z - 10.0).abs() < 1e-6,
            "{l:?} {h:?}"
        );
    }
}

/// `∫₀^m A(s) ds` for the unit round's spandrel section `A(s)`: what of
/// it lies at `s` or more from the face its band leaves (the area between
/// the quarter circle and its corner square, over `y ≥ s`).
fn spandrel_integral(m: f64) -> f64 {
    let a = 1.0 - m;
    let disc = 0.5 * (a * (1.0 - a * a).sqrt() + a.asin());
    0.5 * (1.0 - a * a) - PI / 4.0 + a * disc + (1.0 - a * a).powf(1.5) / 3.0
}

/// The first top edge alone rounded 1, ending at the kink (ADR 0001,
/// amendment "corners where a round's rail misses"): the band runs on past
/// the vertex while the turned front trims it, until the trim reaches its
/// contact on the top, `cot(kink)` past the vertex, or the block's far end
/// stops it first. At `x` past the vertex the band takes what of its
/// spandrel lies at `x tan(kink)` or more from the front, so it takes
/// `20 (1 − π/4) + cot(kink) ∫₀^m A(s) ds`, `m` the lesser of 1 and
/// `(far − 20) tan(kink)`. The same whether the far end stops it or not.
#[test]
fn a_round_ending_at_a_kink_runs_on_until_the_turned_front_trims_it() {
    for (degrees, far) in [(4.0_f64, 40.0), (4.0, 30.0), (2.0, 40.0), (0.9, 40.0)] {
        let kink = degrees.to_radians();
        let mut m = Model::new();
        let (s, first, _) = kinked(&mut m, kink, far);
        let before = volume(&m, &s);
        let built = fillet_edges(&mut m, &s, &[first], 1.0, T).unwrap();
        let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{degrees}°, {far}: {:?}", found.problems);
        let reach = ((far - 20.0) * kink.tan()).min(1.0);
        let want = 20.0 * (1.0 - PI / 4.0) + spandrel_integral(reach) / kink.tan();
        let got = before - volume(&m, &built.shape);
        assert!(
            (got - want).abs() <= 1e-6 * want,
            "{degrees}°, {far}: took {got} against {want}"
        );
    }
}

/// At a 4 degree kink with the far end at 40, the first top edge picked
/// with tangent propagation brings no other edge, and the
/// round of that chain takes `20 (1 − π/4) + cot(4°)(5/6 − π/4)`,
/// 4.977542, not the 4.292037 of a band capped square at the kink.
#[test]
fn a_tangent_pick_at_a_four_degree_kink_rounds_its_first_edge_alone() {
    let kink = 4.0_f64.to_radians();
    let mut m = Model::new();
    let (s, first, _) = kinked(&mut m, kink, 40.0);
    let chain = tangent_chain(&m, &s, &[first], TANGENT_CHAIN_ANGLE, T).unwrap();
    assert_eq!(chain.len(), 1);
    let before = volume(&m, &s);
    let built = fillet_edges(&mut m, &s, &chain, 1.0, T).unwrap();
    let want = 20.0 * (1.0 - PI / 4.0) + (5.0 / 6.0 - PI / 4.0) / kink.tan();
    assert!((want - 4.977542).abs() < 5e-7, "{want}");
    let got = before - volume(&m, &built.shape);
    assert!((got - want).abs() <= 1e-6 * want, "took {got}");
}
