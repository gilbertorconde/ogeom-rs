//! Tangent chains (ADR 0001, amendments "tangent chains stop at a kink"
//! and "tangent chains through fading creases"): which edges a pick runs
//! on to, and what a round ending at a kink takes.
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom_core::Tolerances;
use ogeom_fillet::{
    ChainRules, TANGENT_CHAIN_ANGLE, fillet_edges, tangent_chain, tangent_chain_with,
};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Vector};
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

/// `a` and `b` hold the same edges, in any order.
fn same_set(a: &[Shape], b: &[Shape]) -> bool {
    a.len() == b.len() && a.iter().all(|x| b.iter().any(|y| x.is_same(y)))
}

/// 33 points along `e`, ends included.
fn samples(m: &Model, e: &Shape) -> Vec<Point> {
    let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
    (0..=32)
        .map(|i| c.point_at(a + (b - a) * f64::from(i) / 32.0, T).unwrap())
        .collect()
}

/// The unit direction `e` leaves its end at `p` in.
fn leaving(m: &Model, e: &Shape, p: Point) -> Vector {
    let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
    let (u, sign) = if c.point_at(a, T).unwrap().distance(p) < 1e-6 {
        (a, 1.0)
    } else {
        assert!(c.point_at(b, T).unwrap().distance(p) < 1e-6, "{p:?}");
        (b, -1.0)
    };
    let d = c.d1_at(u, T).unwrap() * sign;
    d / d.magnitude()
}

fn angle(a: Vector, b: Vector) -> f64 {
    a.dot(b).clamp(-1.0, 1.0).acos()
}

/// Two drums of radius 5 crossing at right angles, about the x and y axes,
/// 30 long, fused: their seam is the two ellipses where the drums meet, in
/// the planes `y = x` and `y = -x`, crossing at `(0, 0, ±5)` where the
/// drums touch tangentially.
fn crossed_drums(m: &mut Model) -> Shape {
    let x = Frame::new(Point::new(-15.0, 0.0, 0.0), Direction::X, Direction::Y, T).unwrap();
    let y = Frame::new(Point::new(0.0, -15.0, 0.0), Direction::Y, Direction::Z, T).unwrap();
    let a = ogeom_algo::make_cylinder(m, x, 5.0, 30.0, T).unwrap().shape;
    let b = ogeom_algo::make_cylinder(m, y, 5.0, 30.0, T).unwrap().shape;
    ogeom_bool::fuse(m, &a, &b, T).unwrap().shape
}

/// Picking one arc of the drums' seam takes its ellipse by default, and
/// the whole seam through fading creases. The oracle is the drums' own
/// equations: the seam is every edge lying on both `y² + z² = 25` and
/// `x² + z² = 25`, an ellipse the part of it in its plane. At the two
/// points the ellipses cross, the arcs of the two ellipses meet at a right
/// angle (measured from the edges' tangents), and each drum's outward
/// normal, `(0, y, z)/5` and `(x, 0, z)/5`, agrees with the other's within
/// the bound (measured at the vertex the edges share).
#[test]
fn a_seam_through_fading_creases_is_one_chain_with_through_fading() {
    let mut m = Model::new();
    let s = crossed_drums(&mut m);
    let edges = explore_unique(&m, &s, ShapeType::Edge).unwrap();
    let on = |e: &Shape, f: &dyn Fn(Point) -> bool| samples(&m, e).into_iter().all(f);
    let drums =
        |p: Point| (p.y.hypot(p.z) - 5.0).abs() < 1e-6 && (p.x.hypot(p.z) - 5.0).abs() < 1e-6;
    let seam: Vec<Shape> = edges.iter().filter(|e| on(e, &drums)).cloned().collect();
    let rising: Vec<Shape> = seam
        .iter()
        .filter(|e| on(e, &|p: Point| (p.x - p.y).abs() < 1e-6))
        .cloned()
        .collect();
    assert!(
        seam.len() > rising.len() && !rising.is_empty(),
        "{} {}",
        seam.len(),
        rising.len()
    );

    let picked = edge_near(&m, &s, Point::new(3.535_533_9, 3.535_533_9, 3.535_533_9));
    assert!(rising.iter().any(|e| e.is_same(&picked)));
    let default = tangent_chain(
        &m,
        &s,
        std::slice::from_ref(&picked),
        TANGENT_CHAIN_ANGLE,
        T,
    )
    .unwrap();
    assert!(same_set(&default, &rising), "{} edges", default.len());
    let rules = ChainRules {
        through_fading: true,
        ..ChainRules::default()
    };
    let all = tangent_chain_with(&m, &s, std::slice::from_ref(&picked), rules, T).unwrap();
    assert!(
        same_set(&all, &seam),
        "{} of {} edges",
        all.len(),
        seam.len()
    );

    for z in [5.0, -5.0] {
        let p = Point::new(0.0, 0.0, z);
        let normals = [
            Vector::new(0.0, p.y, p.z) / 5.0,
            Vector::new(p.x, 0.0, p.z) / 5.0,
        ];
        let fade = angle(normals[0], normals[1]);
        eprintln!("z = {z}: drums' normals {fade} rad apart, bound {TANGENT_CHAIN_ANGLE}");
        assert!(fade <= TANGENT_CHAIN_ANGLE);
        let at: Vec<&Shape> = seam
            .iter()
            .filter(|e| samples(&m, e).iter().any(|q| q.distance(p) < 1e-6))
            .collect();
        assert_eq!(at.len(), 4, "four arcs meet at {p:?}");
        let across: Vec<f64> = at
            .iter()
            .filter(|e| !rising.iter().any(|r| r.is_same(e)))
            .map(|e| {
                angle(leaving(&m, e, p), leaving(&m, &rising[0], p))
                    .min(angle(leaving(&m, e, p), -leaving(&m, &rising[0], p)))
            })
            .collect();
        // The other ellipse's arcs leave square to this one's: a turn of
        // π/2, far past the bound, so they join by fading alone.
        for turn in across {
            eprintln!("z = {z}: arcs across turn {turn} rad");
            assert!((turn - PI / 2.0).abs() < 1e-6, "{turn}");
        }
    }
}

/// A box's top edge meets its neighbours square, creased all along: no
/// rule takes another edge.
#[test]
fn a_box_edge_stops_at_its_corners_under_every_rule() {
    let mut m = Model::new();
    let block = ogeom_algo::make_box(
        &mut m,
        Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap(),
        (30.0, 20.0, 10.0),
        T,
    )
    .unwrap()
    .shape;
    let edge = edge_near(&m, &block, Point::new(15.0, 0.0, 10.0));
    let rules = ChainRules {
        through_fading: true,
        through_cusps: true,
        ..ChainRules::default()
    };
    let chain = tangent_chain_with(&m, &block, std::slice::from_ref(&edge), rules, T).unwrap();
    assert!(same(&chain, &[edge]), "{} edges", chain.len());
}

/// A prism 5 high on a profile with a cusp at the origin: the line from it
/// to `(10, 0)`, a line on to `(6, 2)`, and back to the origin the arc of
/// the circle about `(0, 10)` of radius 10, which leaves the origin along
/// `+x` as the first line does. With the prism's top edges on the first
/// line and the arc.
fn cusp_prism(m: &mut Model) -> (Shape, Shape, Shape) {
    let tip = ogeom_algo::make_vertex(m, Point::new(0.0, 0.0, 0.0)).shape;
    let foot = ogeom_algo::make_vertex(m, Point::new(10.0, 0.0, 0.0)).shape;
    let top = ogeom_algo::make_vertex(m, Point::new(6.0, 2.0, 0.0)).shape;
    let seg = |m: &mut Model, from: (&Shape, Point), to: (&Shape, Point)| {
        let curve =
            ogeom_geom::Curve::Line(ogeom_geom::LineCurve::segment(from.1, to.1, T).unwrap());
        let domain = ogeom_geom::Curve3d::domain(&curve);
        ogeom_algo::make_edge_between(m, curve, domain, from.0, to.0, T)
            .unwrap()
            .shape
    };
    let base = seg(
        m,
        (&tip, Point::ORIGIN),
        (&foot, Point::new(10.0, 0.0, 0.0)),
    );
    let slope = seg(
        m,
        (&foot, Point::new(10.0, 0.0, 0.0)),
        (&top, Point::new(6.0, 2.0, 0.0)),
    );
    let frame = Frame::new(Point::new(0.0, 10.0, 0.0), Direction::Z, -Direction::Y, T).unwrap();
    let circle = ogeom_geom::Curve::Circle(ogeom_geom::CircleCurve::new(
        Circle::new(frame, 10.0, T).unwrap(),
    ));
    let arc = ogeom_algo::make_edge_between(m, circle, (0.0, 0.6_f64.asin()), &tip, &top, T)
        .unwrap()
        .shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face_with_pcurves(
        m,
        plane.into(),
        &[vec![base, slope, arc.reversed()]],
        T,
    )
    .unwrap()
    .shape;
    let swept = ogeom_algo::make_prism(m, &face, Vector::new(0.0, 0.0, 5.0), T)
        .unwrap()
        .shape;
    let s = ogeom_algo::baked_shape(m, &swept, T).unwrap().shape;
    let line = edge_near(m, &s, Point::new(5.0, 0.0, 5.0));
    let arc = edge_near(m, &s, Point::new(3.0, 10.0 - 91.0_f64.sqrt(), 5.0));
    (s, line, arc)
}

/// The cusp pair joins only with `through_cusps`: at the cusp the two
/// edges leave the same way (measured from their tangents, within the
/// bound), and each crease stays square there, so fading takes neither.
#[test]
fn a_cusp_pair_is_one_chain_only_with_through_cusps() {
    let mut m = Model::new();
    let (s, line, arc) = cusp_prism(&mut m);
    let tip = Point::new(0.0, 0.0, 5.0);
    let turn = angle(leaving(&m, &line, tip), leaving(&m, &arc, tip));
    eprintln!("cusp: the edges leave {turn} rad apart, bound {TANGENT_CHAIN_ANGLE}");
    assert!(turn <= TANGENT_CHAIN_ANGLE, "{turn}");

    let pick = std::slice::from_ref(&line);
    let default = tangent_chain(&m, &s, pick, TANGENT_CHAIN_ANGLE, T).unwrap();
    assert!(same(&default, pick), "{} edges", default.len());
    let fading = ChainRules {
        through_fading: true,
        ..ChainRules::default()
    };
    let faded = tangent_chain_with(&m, &s, pick, fading, T).unwrap();
    assert!(same(&faded, pick), "{} edges", faded.len());
    let cusps = ChainRules {
        through_cusps: true,
        ..ChainRules::default()
    };
    let both = tangent_chain_with(&m, &s, pick, cusps, T).unwrap();
    assert!(same(&both, &[line, arc]), "{} edges", both.len());
}
