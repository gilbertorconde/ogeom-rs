//! Rolling corners where curved faces take part (ADR 0001, amendment
//! "rolling corners on curved faces"): every request `Setback` builds on
//! pocket outlines with arcs, tangent chains, a pyramid's apex, rounds on
//! rounds and a kidney pocket builds with `Rolling` too, sound, to the
//! closed form where there is one, and alike in any edge order.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom_core::Tolerances;
use ogeom_fillet::{FilletCorner, fillet_edges, fillet_edges_corner};
use ogeom_geom::{CircleCurve, Curve, Curve3d as _, LineCurve};
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// What a round of radius 1 takes off or fills in per unit length of a
/// right-angled edge: the square less the quarter disc.
const Q: f64 = 1.0 - PI / 4.0;

/// How far that section's centroid stands from each of the edge's faces.
fn qx() -> f64 {
    (5.0 / 6.0 - PI / 4.0) / Q
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

/// A piece of a profile in the plane z = 0, in order round it.
#[derive(Clone, Copy)]
enum Seg {
    Line([f64; 2], [f64; 2]),
    /// From, through, to.
    Arc([f64; 2], [f64; 2], [f64; 2]),
}

impl Seg {
    fn start(self) -> [f64; 2] {
        match self {
            Seg::Line(a, _) | Seg::Arc(a, _, _) => a,
        }
    }
}

/// The circle through three points: its centre and radius.
fn circle_through(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> ([f64; 2], f64) {
    let d = 2.0 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
    let n = |p: [f64; 2]| p[0] * p[0] + p[1] * p[1];
    let x = (n(a) * (b[1] - c[1]) + n(b) * (c[1] - a[1]) + n(c) * (a[1] - b[1])) / d;
    let y = (n(a) * (c[0] - b[0]) + n(b) * (a[0] - c[0]) + n(c) * (b[0] - a[0])) / d;
    ([x, y], ((a[0] - x).powi(2) + (a[1] - y).powi(2)).sqrt())
}

/// The prism of the closed profile `segs` (anticlockwise) from z = `z0`,
/// `height` up.
fn prism(m: &mut Model, segs: &[Seg], z0: f64, height: f64) -> Shape {
    let p = |q: [f64; 2]| Point::new(q[0], q[1], z0);
    let vertices: Vec<Shape> = segs
        .iter()
        .map(|s| ogeom_algo::make_vertex(m, p(s.start())).shape)
        .collect();
    let n = segs.len();
    let mut edges = Vec::with_capacity(n);
    for (i, seg) in segs.iter().enumerate() {
        let (from, to) = (&vertices[i], &vertices[(i + 1) % n]);
        let edge = match *seg {
            Seg::Line(a, b) => {
                let curve = Curve::Line(LineCurve::segment(p(a), p(b), T).unwrap());
                let range = curve.domain();
                ogeom_algo::make_edge_between(m, curve, range, from, to, T).unwrap()
            }
            Seg::Arc(a, b, c) => {
                let (o, r) = circle_through(a, b, c);
                let turn = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]);
                let normal = if turn > 0.0 {
                    Direction::Z
                } else {
                    Direction::new(Vector::new(0.0, 0.0, -1.0), T).unwrap()
                };
                let x = Direction::new(Vector::new(a[0] - o[0], a[1] - o[1], 0.0), T).unwrap();
                let frame = Frame::new(p(o), normal, x, T).unwrap();
                let circle = Circle::new(frame, r, T).unwrap();
                let curve = Curve::Circle(CircleCurve::new(circle));
                // The sweep from `a` to `c` about the normal.
                let to_end = Vector::new(c[0] - o[0], c[1] - o[1], 0.0);
                let xv = x.vector();
                let yv = normal.vector().cross(xv);
                let mut sweep = to_end.dot(yv).atan2(to_end.dot(xv));
                if sweep <= 1e-12 {
                    sweep += 2.0 * PI;
                }
                ogeom_algo::make_edge_between(m, curve, (0.0, sweep), from, to, T).unwrap()
            }
        };
        edges.push(edge.shape);
    }
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::new(0.0, 0.0, z0), Direction::Z),
        (-1000.0, 1000.0),
        (-1000.0, 1000.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face_with_pcurves(m, plane.into(), &[edges], T)
        .unwrap()
        .shape;
    let swept = ogeom_algo::make_prism(m, &face, Vector::new(0.0, 0.0, height), T)
        .unwrap()
        .shape;
    ogeom_algo::baked_shape(m, &swept, T).unwrap().shape
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Seg> {
    vec![
        Seg::Line([x0, y0], [x1, y0]),
        Seg::Line([x1, y0], [x1, y1]),
        Seg::Line([x1, y1], [x0, y1]),
        Seg::Line([x0, y1], [x0, y0]),
    ]
}

/// The only solid of `s`.
fn solid_of(m: &Model, s: &Shape) -> Shape {
    let solids = explore_unique(m, s, ShapeType::Solid).unwrap();
    assert_eq!(solids.len(), 1);
    solids[0].clone()
}

/// A block 30 × 30 × 10 with the profile `pocket` cut 5 deep into its top.
fn pocketed(m: &mut Model, pocket: &[Seg]) -> Shape {
    let block = prism(m, &rect(0.0, 0.0, 30.0, 30.0), 0.0, 10.0);
    let tool = prism(m, pocket, 5.0, 6.0);
    let cut = ogeom_bool::cut(m, &block, &tool, T).unwrap().shape;
    solid_of(m, &cut)
}

/// The edges of the face of `s` holding `p`.
fn edges_of_face_at(m: &Model, s: &Shape, p: Point) -> Vec<Shape> {
    let face = explore_unique(m, s, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            ogeom_algo::classify_on_face(m, f, p, ogeom_mesh::Deflection::default(), T)
                .is_ok_and(|c| c == ogeom_algo::Containment::In)
        })
        .unwrap();
    explore_unique(m, &face, ShapeType::Edge).unwrap()
}

/// The edge of `s` passing nearest `q`, 32 samples along each.
fn edge_near(m: &Model, s: &Shape, q: Point) -> Shape {
    let d = |e: &Shape| {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, e, T).unwrap() else {
            return f64::INFINITY;
        };
        (0..=32)
            .map(|i| {
                c.point_at(a + (b - a) * f64::from(i) / 32.0, T)
                    .unwrap()
                    .distance(q)
            })
            .fold(f64::INFINITY, f64::min)
    };
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| d(a).total_cmp(&d(b)))
        .unwrap()
}

/// `edges` in a shuffled order, by a fixed `seed`.
fn shuffled(edges: &[Shape], seed: u64) -> Vec<Shape> {
    let mut out = edges.to_vec();
    let mut x = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
    for i in (1..out.len()).rev() {
        x = x
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let j = usize::try_from(x >> 33).unwrap() % (i + 1);
        out.swap(i, j);
    }
    out
}

/// What a case should come to with rolling corners.
enum Want {
    /// This volume, within a millionth of it.
    Volume(f64),
    /// The volume `Setback` gives, within a millionth of it.
    AsSetback,
}

/// `edges` of `s` rounded `r` with `Setback`, then with `Rolling` in three
/// orders: each builds, sound, to `want`, and the orders agree.
fn rolls(name: &str, m: &Model, s: &Shape, edges: &[Shape], r: f64, want: Want) {
    let mut setback = m.clone();
    let built = fillet_edges(&mut setback, s, edges, r, T)
        .unwrap_or_else(|e| panic!("{name}: Setback refused: {e}"));
    let by_setback = volume(&setback, &built.shape);
    let want = match want {
        Want::Volume(v) => v,
        Want::AsSetback => by_setback,
    };
    let mut first: Option<f64> = None;
    for seed in 0..3 {
        let mut rolled = m.clone();
        let order = shuffled(edges, seed);
        let built = fillet_edges_corner(&mut rolled, s, &order, r, FilletCorner::Rolling, T)
            .unwrap_or_else(|e| panic!("{name}, order {seed}: Rolling refused: {e}"));
        let found = ogeom_algo::check(&rolled, &built.shape, T).unwrap();
        assert!(
            found.is_valid(),
            "{name}, order {seed}: {:?}",
            found.problems
        );
        let got = volume(&rolled, &built.shape);
        assert!(
            (got - want).abs() <= 1e-6 * want.abs(),
            "{name}, order {seed}: volume {got} against {want}"
        );
        match first {
            None => first = Some(got),
            Some(v) => assert!(
                (got - v).abs() <= 1e-9 * v.abs(),
                "{name}: order {seed} gives {got} against {v}"
            ),
        }
        // Every new face is the ball's: a plane, a cylinder, a cone, a
        // sphere or a torus, never a patch.
        for face in explore_unique(&rolled, &built.shape, ShapeType::Face).unwrap() {
            let surface = rolled
                .node(&face)
                .and_then(|n| n.data().as_face())
                .and_then(|d| rolled.geometry().surface(d.surface))
                .unwrap();
            assert!(
                !matches!(surface, ogeom_geom::SurfaceGeometry::BSpline(_)),
                "{name}, order {seed}: a corner patch"
            );
        }
    }
}

/// The obround pocket's floor loop: two lines and two half circles of
/// radius 3, tangent all round, rounded 1 (concave). The bands run on
/// through every junction: the pocket plus the section swept round the
/// loop, its centroid `qx` in from the wall.
#[test]
fn an_obround_pocket_s_floor_rolls_round_its_tangent_loop() {
    let mut m = Model::new();
    let s = pocketed(
        &mut m,
        &[
            Seg::Line([8.0, 12.0], [22.0, 12.0]),
            Seg::Arc([22.0, 12.0], [25.0, 15.0], [22.0, 18.0]),
            Seg::Line([22.0, 18.0], [8.0, 18.0]),
            Seg::Arc([8.0, 18.0], [5.0, 15.0], [8.0, 12.0]),
        ],
    );
    let floor = edges_of_face_at(&m, &s, Point::new(15.0, 15.0, 5.0));
    assert_eq!(floor.len(), 4);
    let base = 9000.0 - (14.0 * 6.0 + PI * 9.0) * 5.0;
    let want = base + Q * (28.0 + 2.0 * PI * (3.0 - qx()));
    rolls("obround floor", &m, &s, &floor, 1.0, Want::Volume(want));
}

/// The kidney pocket's floor loop: four tangent arcs (radii 12 and 6 about
/// (15, 8) over 120 degrees, ends of radius 3), rounded 1.
#[test]
fn a_kidney_pocket_s_floor_rolls_round_its_arcs() {
    let (c, a0, a1) = ([15.0, 8.0], PI / 6.0, 5.0 * PI / 6.0);
    let p = |r: f64, a: f64| [c[0] + r * a.cos(), c[1] + r * a.sin()];
    let am = 0.5 * (a0 + a1);
    let end = |a: f64, out: f64| {
        let q = p(9.0, a);
        [q[0] - 3.0 * out * a.sin(), q[1] + 3.0 * out * a.cos()]
    };
    let mut m = Model::new();
    let s = pocketed(
        &mut m,
        &[
            Seg::Arc(p(12.0, a0), p(12.0, am), p(12.0, a1)),
            Seg::Arc(p(12.0, a1), end(a1, 1.0), p(6.0, a1)),
            Seg::Arc(p(6.0, a1), p(6.0, am), p(6.0, a0)),
            Seg::Arc(p(6.0, a0), end(a0, -1.0), p(12.0, a0)),
        ],
    );
    let floor = edges_of_face_at(&m, &s, Point::new(15.0, 17.0, 5.0));
    assert_eq!(floor.len(), 4);
    let base = 9000.0 - (PI / 3.0 * (144.0 - 36.0) + PI * 9.0) * 5.0;
    let want = base + Q * (12.0 * PI + 2.0 * PI * (3.0 - qx()));
    rolls("kidney floor", &m, &s, &floor, 1.0, Want::Volume(want));
}

/// A square pocket 10 × 10 × 5: its floor loop rounded 1. At each corner
/// two concave rounds meet beside the upright left sharp and cross: the
/// pocket plus the section along 40, less what the crossings share.
#[test]
fn a_square_pocket_s_floor_rolls_as_its_bands_cross() {
    let mut m = Model::new();
    let s = pocketed(&mut m, &rect(10.0, 10.0, 20.0, 20.0));
    let floor = edges_of_face_at(&m, &s, Point::new(15.0, 15.0, 5.0));
    assert_eq!(floor.len(), 4);
    let want = 8500.0 + Q * (40.0 - 8.0 * qx());
    rolls(
        "square pocket floor",
        &m,
        &s,
        &floor,
        1.0,
        Want::Volume(want),
    );
}

/// The obround pocket's rim loop: the same tangent loop on the block's
/// top, rounded 1 (convex): the section swept round the loop, its
/// centroid `qx` out from the wall.
#[test]
fn an_obround_pocket_s_rim_rolls_round_its_tangent_loop() {
    let mut m = Model::new();
    let s = pocketed(
        &mut m,
        &[
            Seg::Line([8.0, 12.0], [22.0, 12.0]),
            Seg::Arc([22.0, 12.0], [25.0, 15.0], [22.0, 18.0]),
            Seg::Line([22.0, 18.0], [8.0, 18.0]),
            Seg::Arc([8.0, 18.0], [5.0, 15.0], [8.0, 12.0]),
        ],
    );
    let rim: Vec<Shape> = [(15.0, 12.0), (25.0, 15.0), (15.0, 18.0), (5.0, 15.0)]
        .iter()
        .map(|&(x, y)| edge_near(&m, &s, Point::new(x, y, 10.0)))
        .collect();
    let base = 9000.0 - (14.0 * 6.0 + PI * 9.0) * 5.0;
    let want = base - Q * (28.0 + 2.0 * PI * (3.0 + qx()));
    rolls("obround rim", &m, &s, &rim, 1.0, Want::Volume(want));
}

/// A block 10 × 10 × 10 notched at an upright by a cylinder of radius 4:
/// the three convex edges at the corner where the top, the back and the
/// notch meet, rounded 1, close on the ball's sphere touching the two
/// planes and the cylinder.
#[test]
fn a_notched_corner_s_three_rounds_roll_on_the_ball() {
    let mut m = Model::new();
    let body = prism(&mut m, &rect(0.0, 0.0, 10.0, 10.0), 0.0, 10.0);
    let drill = ogeom_algo::make_cylinder(
        &mut m,
        Frame::new(Point::new(10.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap(),
        4.0,
        12.0,
        T,
    )
    .unwrap()
    .shape;
    let notched = ogeom_bool::cut(&mut m, &body, &drill, T).unwrap().shape;
    let s = solid_of(&m, &notched);
    let corner = vec![
        edge_near(&m, &s, Point::new(3.0, 10.0, 10.0)),
        edge_near(&m, &s, Point::new(6.0, 10.0, 5.0)),
        edge_near(&m, &s, Point::new(10.0 - 4.0 * 0.8, 10.0 - 4.0 * 0.6, 10.0)),
    ];
    rolls("notched corner", &m, &s, &corner, 1.0, Want::AsSetback);
}

/// A block 20 × 20 × 10 with its top loop rounded 2, then its uprights
/// rounded 1: each upright's band runs up into the top's bands, meeting
/// them on their cylinders.
#[test]
#[ignore = "heavy"]
fn uprights_rounded_after_the_top_roll_into_its_rounds() {
    let mut m = Model::new();
    let block = prism(&mut m, &rect(0.0, 0.0, 20.0, 20.0), 0.0, 10.0);
    let top = edges_of_face_at(&m, &block, Point::new(10.0, 10.0, 10.0));
    let rounded = fillet_edges(&mut m, &block, &top, 2.0, T).unwrap().shape;
    let uprights: Vec<Shape> = [(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0)]
        .iter()
        .map(|&(x, y)| edge_near(&m, &rounded, Point::new(x, y, 4.0)))
        .collect();
    rolls(
        "uprights after the top",
        &m,
        &rounded,
        &uprights,
        1.0,
        Want::AsSetback,
    );
}

/// A plate 20 × 20 × 4 under a D-shaped boss (a half disc of radius 4 on
/// its chord, 6 high): its foot loop and the two uprights at the chord's
/// ends rounded 1. At each foot corner two concave rounds and a convex one
/// meet with the boss's cylinder there; the ball's sweep round the
/// upright's end is no torus on a curved face, and `Setback` refuses the
/// corner too. Rolling refuses it by name, naming the vertex.
#[test]
fn a_mixed_corner_on_a_curved_face_is_refused_by_name() {
    let mut m = Model::new();
    let plate = prism(&mut m, &rect(0.0, 0.0, 20.0, 20.0), 0.0, 4.0);
    let boss = prism(
        &mut m,
        &[
            Seg::Line([6.0, 10.0], [14.0, 10.0]),
            Seg::Arc([14.0, 10.0], [10.0, 14.0], [6.0, 10.0]),
        ],
        4.0,
        6.0,
    );
    let s = ogeom_bool::fuse(&mut m, &plate, &boss, T).unwrap().shape;
    let s = solid_of(&m, &s);
    let edges: Vec<Shape> = [
        Point::new(10.0, 10.0, 4.0),
        Point::new(10.0, 14.0, 4.0),
        Point::new(6.0, 10.0, 7.0),
        Point::new(14.0, 10.0, 7.0),
    ]
    .iter()
    .map(|&p| edge_near(&m, &s, p))
    .collect();
    assert!(fillet_edges(&mut m.clone(), &s, &edges, 1.0, T).is_err());
    let why = fillet_edges_corner(&mut m, &s, &edges, 1.0, FilletCorner::Rolling, T)
        .unwrap_err()
        .to_string();
    assert!(
        why.contains("does not close the vertex at") && why.contains("curved"),
        "{why}"
    );
}

/// A rounded rectangle 30 × 20 (corners of radius 4), 10 high: its top
/// loop rounded 2, a tangent chain of lines and arcs. The section swept
/// round the loop, its centroid `2 qx` in from the edge.
#[test]
fn a_rounded_rectangle_s_top_loop_rolls_on_through_its_arcs() {
    let s2 = 4.0 * core::f64::consts::FRAC_1_SQRT_2;
    let arc = |a: [f64; 2], c: [f64; 2], b: [f64; 2]| {
        Seg::Arc(
            a,
            [
                c[0] + s2 * (a[0] + b[0] - 2.0 * c[0]).signum(),
                c[1] + s2 * (a[1] + b[1] - 2.0 * c[1]).signum(),
            ],
            b,
        )
    };
    let mut m = Model::new();
    let s = prism(
        &mut m,
        &[
            Seg::Line([4.0, 0.0], [26.0, 0.0]),
            arc([26.0, 0.0], [26.0, 4.0], [30.0, 4.0]),
            Seg::Line([30.0, 4.0], [30.0, 16.0]),
            arc([30.0, 16.0], [26.0, 16.0], [26.0, 20.0]),
            Seg::Line([26.0, 20.0], [4.0, 20.0]),
            arc([4.0, 20.0], [4.0, 16.0], [0.0, 16.0]),
            Seg::Line([0.0, 16.0], [0.0, 4.0]),
            arc([0.0, 4.0], [4.0, 4.0], [4.0, 0.0]),
        ],
        0.0,
        10.0,
    );
    let top = edges_of_face_at(&m, &s, Point::new(15.0, 10.0, 10.0));
    assert_eq!(top.len(), 8);
    let base = (600.0 - (4.0 - PI) * 16.0) * 10.0;
    let want = base - 4.0 * Q * (68.0 + 2.0 * PI * (4.0 - 2.0 * qx()));
    rolls(
        "rounded rectangle top",
        &m,
        &s,
        &top,
        2.0,
        Want::Volume(want),
    );
}

/// A block 20 × 20 × 10 with its uprights rounded 3, then its top loop
/// rounded 1: four lines and the four arcs the rounds left, tangent all
/// round, where the rounds meet on the uprights' cylinders.
#[test]
fn a_round_on_rounds_rolls_round_the_rounded_block() {
    let mut m = Model::new();
    let block = prism(&mut m, &rect(0.0, 0.0, 20.0, 20.0), 0.0, 10.0);
    let uprights: Vec<Shape> = [(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0)]
        .iter()
        .map(|&(x, y)| edge_near(&m, &block, Point::new(x, y, 5.0)))
        .collect();
    let rounded = fillet_edges(&mut m, &block, &uprights, 3.0, T)
        .unwrap()
        .shape;
    let top = edges_of_face_at(&m, &rounded, Point::new(10.0, 10.0, 10.0));
    assert_eq!(top.len(), 8);
    let base = (400.0 - (4.0 - PI) * 9.0) * 10.0;
    let want = base - Q * (56.0 + 2.0 * PI * (3.0 - qx()));
    rolls(
        "round on rounds",
        &m,
        &rounded,
        &top,
        1.0,
        Want::Volume(want),
    );
}

/// A square pyramid, base 20 × 20 and apex 10 over its middle: its four
/// edges at the apex rounded 1, the apex the ball's own sphere where it
/// touches all four faces.
#[test]
fn a_pyramid_s_apex_of_four_rounds_rolls_on_the_ball() {
    let mut m = Model::new();
    let points = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(20.0, 0.0, 0.0),
        Point::new(20.0, 20.0, 0.0),
        Point::new(0.0, 20.0, 0.0),
        Point::new(10.0, 10.0, 10.0),
    ];
    let rings = vec![
        vec![0, 3, 2, 1],
        vec![0, 1, 4],
        vec![1, 2, 4],
        vec![2, 3, 4],
        vec![3, 0, 4],
    ];
    let s = ogeom_algo::make_polyhedron(&mut m, &points, &rings, T)
        .unwrap()
        .shape;
    let apex: Vec<Shape> = [(5.0, 5.0), (15.0, 5.0), (15.0, 15.0), (5.0, 15.0)]
        .iter()
        .map(|&(x, y)| edge_near(&m, &s, Point::new(x, y, 5.0)))
        .collect();
    rolls("pyramid apex", &m, &s, &apex, 1.0, Want::AsSetback);
}
