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

/// `edges` of `s` rounded `r` with `Rolling` in three orders, where
/// `Setback` refuses: each builds, sound, to `want` within a millionth,
/// the orders agree, every new face is the ball's (no patch), and every
/// edge the rounds leave between two of their own faces is smooth.
fn rolls_alone(name: &str, m: &Model, s: &Shape, edges: &[Shape], r: f64, want: f64) {
    rolls_to(name, m, s, edges, r, want, true);
}

/// [`rolls_alone`], each edge left between two curved faces checked
/// smooth only where `smooth`: a band ending against a curved wall at a
/// sharp corner meets it along a crease.
fn rolls_to(name: &str, m: &Model, s: &Shape, edges: &[Shape], r: f64, want: f64, smooth: bool) {
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
        if smooth {
            smooth_between_rounds(name, &rolled, &built.shape);
        }
    }
}

/// Every edge of `s` between two faces neither of them a plane is G1:
/// the faces' normals agree at its middle within a millionth of a radian.
/// A plane meets a round only along the round's contact, G1 by its
/// construction; two curved faces meet where the rounds join.
fn smooth_between_rounds(name: &str, m: &Model, s: &Shape) {
    let faces = explore_unique(m, s, ShapeType::Face).unwrap();
    let curved = |f: &Shape| {
        let surface = m
            .node(f)
            .and_then(|n| n.data().as_face())
            .and_then(|d| m.geometry().surface(d.surface))
            .unwrap();
        !matches!(surface, ogeom_geom::SurfaceGeometry::Plane(_))
    };
    for edge in explore_unique(m, s, ShapeType::Edge).unwrap() {
        let on: Vec<&Shape> = faces
            .iter()
            .filter(|f| {
                explore_unique(m, f, ShapeType::Edge)
                    .unwrap()
                    .iter()
                    .any(|e| e.is_same(&edge))
            })
            .collect();
        let [a, b] = on.as_slice() else { continue };
        if !(curved(a) && curved(b)) {
            continue;
        }
        let Some((c, (t0, t1))) = ogeom_algo::placed_edge_curve(m, &edge, T).unwrap() else {
            continue;
        };
        let mid = c.point_at(0.5 * (t0 + t1), T).unwrap();
        let normal = |f: &Shape| {
            use ogeom_geom::{Surface as _, Transformable as _};
            let surface = m
                .node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| m.geometry().surface(d.surface))
                .unwrap()
                .clone()
                .transformed(&f.transform(m.datums()).unwrap(), T)
                .unwrap();
            let (u, v) = ogeom_algo::project_on_surface(&surface, mid, 16, T)
                .unwrap()
                .parameters;
            surface.normal_at(u, v, T).unwrap().vector()
        };
        // Up to the faces' senses: a crease is a turn, not a flip.
        let turn = normal(a).dot(normal(b)).abs().clamp(0.0, 1.0).acos();
        assert!(turn <= 1e-6, "{name}: a crease of {turn} rad at {mid:?}");
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

/// The D of radius `big` on its chord with its two corners rounded `r`
/// (each round tangent to the chord and to the arc): its area and its
/// perimeter.
fn rounded_d(big: f64, r: f64) -> (f64, f64) {
    let x = ((big - r).powi(2) - r * r).sqrt();
    let phi = (r / x).atan();
    let corner = big * big * phi / 2.0 - x * r / 2.0 - r * r * (PI / 2.0 + phi) / 2.0;
    let area = PI * big * big / 2.0 - 2.0 * corner;
    let perimeter = 2.0 * x + 2.0 * r * (PI / 2.0 + phi) + big * (PI - 2.0 * phi);
    (area, perimeter)
}

/// What a round of radius `r` adds round (or takes within) a convex
/// outline of perimeter `p`: at `t` from the round's flat face its section
/// is the outline grown (or shrunk) by `w = r - sqrt(r² - (r - t)²)`, so
/// over the round's height the area changes by `p ∫w ± π ∫w²`, with
/// `∫w = r²(1 - π/4)` and `∫w² = r³(5/3 - π/2)`.
fn grown(p: f64, r: f64, outward: bool) -> f64 {
    let (w, w2) = (r * r * (1.0 - PI / 4.0), r.powi(3) * (5.0 / 3.0 - PI / 2.0));
    if outward {
        p * w + PI * w2
    } else {
        p * w - PI * w2
    }
}

/// A plate 20 × 20 × 4 under a D-shaped boss: a half disc of radius 4 on
/// its chord, 6 high.
fn d_boss(m: &mut Model) -> Shape {
    let plate = prism(m, &rect(0.0, 0.0, 20.0, 20.0), 0.0, 4.0);
    let boss = prism(
        m,
        &[
            Seg::Line([6.0, 10.0], [14.0, 10.0]),
            Seg::Arc([14.0, 10.0], [10.0, 14.0], [6.0, 10.0]),
        ],
        4.0,
        6.0,
    );
    let s = ogeom_bool::fuse(m, &plate, &boss, T).unwrap().shape;
    solid_of(m, &s)
}

/// The D boss's foot loop and its two uprights.
fn d_boss_foot_and_uprights(m: &Model, s: &Shape) -> Vec<Shape> {
    [
        Point::new(10.0, 10.0, 4.0),
        Point::new(10.0, 14.0, 4.0),
        Point::new(6.0, 10.0, 7.0),
        Point::new(14.0, 10.0, 7.0),
    ]
    .iter()
    .map(|&p| edge_near(m, s, p))
    .collect()
}

/// The D boss's foot loop and its two uprights rounded 1. At each foot
/// corner the two concave rounds and the convex upright meet with the
/// boss's cylinder there; the ball rolls round the upright band's end on
/// the torus about its axis, major radius 2, so every section of the boss
/// at `t` above the plate within the round is the rounded D grown by the
/// round's offset there (ADR 0001, amendment "mixed rolling corners on
/// curved faces"). `Setback` refuses the corner.
#[test]
fn a_d_boss_s_mixed_foot_corners_roll_round_the_uprights() {
    let mut m = Model::new();
    let s = d_boss(&mut m);
    let edges = d_boss_foot_and_uprights(&m, &s);
    assert!(fillet_edges(&mut m.clone(), &s, &edges, 1.0, T).is_err());
    let (area, perimeter) = rounded_d(4.0, 1.0);
    let want = 1600.0 + area * 6.0 + grown(perimeter, 1.0, true);
    rolls_alone("D boss foot", &m, &s, &edges, 1.0, want);
}

/// The same with the boss's top loop too: the top corners are of one
/// sense and close on the ball's sphere, the sections there the rounded D
/// shrunk by the top round's offset.
#[test]
fn a_d_boss_rounded_all_round_rolls_at_its_foot_and_top() {
    let mut m = Model::new();
    let s = d_boss(&mut m);
    let mut edges = d_boss_foot_and_uprights(&m, &s);
    edges.extend(edges_of_face_at(&m, &s, Point::new(10.0, 12.0, 10.0)));
    assert_eq!(edges.len(), 6);
    let (area, perimeter) = rounded_d(4.0, 1.0);
    let want = 1600.0 + area * 6.0 + grown(perimeter, 1.0, true) - grown(perimeter, 1.0, false);
    rolls_alone("D boss all round", &m, &s, &edges, 1.0, want);
}

/// A D pocket (a half disc of radius 5 on its chord, 5 deep): its rim
/// and the two uprights at the chord's ends rounded 1. At each rim corner
/// two convex rounds and the concave upright meet with the pocket's
/// cylinder: the void's section at `t` below the top within the round is
/// the rounded D grown by the round's offset there.
#[test]
fn a_d_pocket_s_mixed_rim_corners_roll_round_the_uprights() {
    let mut m = Model::new();
    let s = pocketed(
        &mut m,
        &[
            Seg::Line([10.0, 12.0], [20.0, 12.0]),
            Seg::Arc([20.0, 12.0], [15.0, 17.0], [10.0, 12.0]),
        ],
    );
    // The boolean leaves the rim's arc in two pieces.
    let h = 5.0 * core::f64::consts::FRAC_1_SQRT_2;
    let edges: Vec<Shape> = [
        Point::new(15.0, 12.0, 10.0),
        Point::new(15.0 + h, 12.0 + h, 10.0),
        Point::new(15.0 - h, 12.0 + h, 10.0),
        Point::new(10.0, 12.0, 7.0),
        Point::new(20.0, 12.0, 7.0),
    ]
    .iter()
    .map(|&p| edge_near(&m, &s, p))
    .collect();
    let (area, perimeter) = rounded_d(5.0, 1.0);
    let want = 9000.0 - area * 5.0 - grown(perimeter, 1.0, true);
    rolls_alone("D pocket rim", &m, &s, &edges, 1.0, want);
}

/// The D boss's foot line and one upright rounded without the foot arc:
/// a convex and a concave round meet beside an edge left sharp, which the
/// ball's sweep does not close. Refused by name, naming the vertex.
#[test]
fn a_mixed_corner_beside_an_edge_left_sharp_on_a_curved_face_is_refused_by_name() {
    let mut m = Model::new();
    let s = d_boss(&mut m);
    let edges = vec![
        edge_near(&m, &s, Point::new(10.0, 10.0, 4.0)),
        edge_near(&m, &s, Point::new(14.0, 10.0, 7.0)),
    ];
    let why = fillet_edges_corner(&mut m, &s, &edges, 1.0, FilletCorner::Rolling, T)
        .unwrap_err()
        .to_string();
    assert!(
        why.contains("does not close the vertex at (14.000000, 10.000000, 4.000000)")
            && why.contains("left sharp"),
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

/// A block 30 × 30 × 10 with a D pocket cut 5 deep: a half disc of radius
/// 5 on its chord from (10, 12) to (20, 12).
fn d_pocket(m: &mut Model) -> Shape {
    pocketed(
        m,
        &[
            Seg::Line([10.0, 12.0], [20.0, 12.0]),
            Seg::Arc([20.0, 12.0], [15.0, 17.0], [10.0, 12.0]),
        ],
    )
}

/// `∫₀¹ f(s) ds` by Simpson's rule after `s = u²`, which takes the
/// `√s` a round's section has at its contact out of the integrand.
fn integral_01(f: impl Fn(f64) -> f64) -> f64 {
    let n = 20_000;
    let h = 1.0 / f64::from(n);
    let g = |u: f64| f(u * u) * 2.0 * u;
    let mut sum = g(0.0) + g(1.0);
    for i in 1..n {
        let w = if i % 2 == 1 { 4.0 } else { 2.0 };
        sum += w * g(f64::from(i) * h);
    }
    sum * h / 3.0
}

/// What the D pocket's floor line rounded 1 fills: at `s` in from the
/// chord its section is `h(s) = 1 - √(2s - s²)` high and runs across the
/// D, `2√(25 - s²)` long.
fn d_floor_line_fill() -> f64 {
    integral_01(|s| 2.0 * (25.0 - s * s).sqrt() * (1.0 - (2.0 * s - s * s).sqrt()))
}

/// The D pocket's floor line rounded 1 alone (concave), its band ending at
/// each end against the pocket's cylinder at a sharp corner. The band runs
/// on and the wall trims it.
#[test]
fn a_d_pocket_s_floor_line_alone_runs_on_until_the_wall_trims_it() {
    let mut m = Model::new();
    let s = d_pocket(&mut m);
    let edges = vec![edge_near(&m, &s, Point::new(15.0, 12.0, 5.0))];
    let want = 9000.0 - 12.5 * PI * 5.0 + d_floor_line_fill();
    rolls_to("D pocket floor line", &m, &s, &edges, 1.0, want, false);
    let mut setback = m.clone();
    let built = fillet_edges(&mut setback, &s, &edges, 1.0, T).unwrap();
    assert!(
        ogeom_algo::check(&setback, &built.shape, T)
            .unwrap()
            .is_valid()
    );
    let got = volume(&setback, &built.shape);
    assert!(
        (got - want).abs() <= 1e-6 * want,
        "Setback: {got} against {want}"
    );
}

/// The D pocket's whole floor loop rounded 1: the line's cylinder band
/// crosses the arc's torus band at each corner, the upright there left
/// sharp, and the fills take their union. The arc's fill is half the
/// section turned about the axis (Pappus, its centroid `qx` in from the
/// wall). Near each corner both fill the points `s` in from the chord and
/// `t` in from the cylinder, up to `h(max(s, t))`; in `(s, t)` the plan's
/// area element is `(5 - t) / √((5 - t)² - s²)`.
#[test]
fn a_d_pocket_s_floor_loop_rolls_as_its_line_and_arc_bands_cross() {
    let mut m = Model::new();
    let s = d_pocket(&mut m);
    let edges = edges_of_face_at(&m, &s, Point::new(15.0, 14.0, 5.0));
    assert_eq!(edges.len(), 2);
    let h = |m: f64| 1.0 - (2.0 * m - m * m).sqrt();
    let area = |s: f64, t: f64| (5.0 - t) / ((5.0 - t).powi(2) - s * s).sqrt();
    let simpson = |f: &dyn Fn(f64) -> f64, b: f64| {
        let n = 200;
        let step = b / f64::from(n);
        let mut sum = f(0.0) + f(b);
        for i in 1..n {
            sum += if i % 2 == 1 { 4.0 } else { 2.0 } * f(f64::from(i) * step);
        }
        sum * step / 3.0
    };
    let shared =
        integral_01(|m| h(m) * (simpson(&|t| area(m, t), m) + simpson(&|s| area(s, m), m)));
    let arc = PI * (5.0 - qx()) * Q;
    let want = 9000.0 - 12.5 * PI * 5.0 + d_floor_line_fill() + arc - 2.0 * shared;
    rolls(
        "D pocket floor loop",
        &m,
        &s,
        &edges,
        1.0,
        Want::Volume(want),
    );
}

/// A piece of a plan outline for the oracle: a line, or an arc about a
/// centre with its radius, anticlockwise where `ccw`.
#[derive(Clone, Copy)]
enum Piece {
    Line([f64; 2], [f64; 2]),
    Arc([f64; 2], f64, bool),
}

/// A piece's carrier offset in plan: a line (a point, a unit direction) or
/// a circle (a centre, a radius).
type Carrier = Result<([f64; 2], [f64; 2]), ([f64; 2], f64)>;

/// The pieces of `segs`, an anticlockwise outline.
fn pieces(segs: &[Seg]) -> Vec<Piece> {
    segs.iter()
        .map(|s| match *s {
            Seg::Line(a, b) => Piece::Line(a, b),
            Seg::Arc(a, b, c) => {
                let (o, r) = circle_through(a, b, c);
                let turn = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]);
                Piece::Arc(o, r, turn > 0.0)
            }
        })
        .collect()
}

/// The area and perimeter of the anticlockwise outline `segs` with each
/// corner rounded `|r|` in plan: each corner's circle centred where the
/// two pieces' carriers offset `r` inwards meet nearest the corner (a
/// negative `r` offsets them outwards, rounding reflex corners), its arc
/// between its feet on them. By Green's theorem piece by piece.
fn rounded_outline(segs: &[Seg], r: f64) -> (f64, f64) {
    let ps = pieces(segs);
    let n = ps.len();
    let sub = |a: [f64; 2], b: [f64; 2]| [a[0] - b[0], a[1] - b[1]];
    let len = |a: [f64; 2]| a[0].hypot(a[1]);
    // The offset carrier of a piece: a line (point, unit direction) or a
    // circle (centre, radius).
    let offset = |p: Piece| -> Carrier {
        match p {
            Piece::Line(a, b) => {
                let d = sub(b, a);
                let l = len(d);
                let u = [d[0] / l, d[1] / l];
                Ok(([a[0] - u[1] * r, a[1] + u[0] * r], u))
            }
            Piece::Arc(o, rad, ccw) => Err((o, if ccw { rad - r } else { rad + r })),
        }
    };
    let foot = |p: Piece, f: [f64; 2]| -> [f64; 2] {
        match p {
            Piece::Line(a, b) => {
                let d = sub(b, a);
                let t = (sub(f, a)[0] * d[0] + sub(f, a)[1] * d[1]) / (d[0] * d[0] + d[1] * d[1]);
                [a[0] + t * d[0], a[1] + t * d[1]]
            }
            Piece::Arc(o, rad, _) => {
                let v = sub(f, o);
                let l = len(v);
                [o[0] + v[0] / l * rad, o[1] + v[1] / l * rad]
            }
        }
    };
    let meet = |a: Carrier, b: Carrier, near: [f64; 2]| -> [f64; 2] {
        let line_circle = |(p, u): ([f64; 2], [f64; 2]), (o, rho): ([f64; 2], f64)| {
            let w = sub(p, o);
            let bq = w[0] * u[0] + w[1] * u[1];
            let c = w[0] * w[0] + w[1] * w[1] - rho * rho;
            let s = (bq * bq - c).sqrt();
            [-bq - s, -bq + s].map(|t| [p[0] + t * u[0], p[1] + t * u[1]])
        };
        let cands: Vec<[f64; 2]> = match (a, b) {
            (Ok((p, u)), Ok((q, v))) => {
                let den = u[0] * v[1] - u[1] * v[0];
                let w = sub(q, p);
                let t = (w[0] * v[1] - w[1] * v[0]) / den;
                vec![[p[0] + t * u[0], p[1] + t * u[1]]]
            }
            (Ok(l), Err(c)) | (Err(c), Ok(l)) => line_circle(l, c).to_vec(),
            (Err((o1, r1)), Err((o2, r2))) => {
                let d = sub(o2, o1);
                let dl = len(d);
                let a = (r1 * r1 - r2 * r2 + dl * dl) / (2.0 * dl);
                let h = (r1 * r1 - a * a).sqrt();
                let m = [o1[0] + a * d[0] / dl, o1[1] + a * d[1] / dl];
                vec![
                    [m[0] - h * d[1] / dl, m[1] + h * d[0] / dl],
                    [m[0] + h * d[1] / dl, m[1] - h * d[0] / dl],
                ]
            }
        };
        cands
            .into_iter()
            .min_by(|x, y| len(sub(*x, near)).total_cmp(&len(sub(*y, near))))
            .unwrap()
    };
    // Each corner k sits at the start of piece k: its centre and its feet
    // on piece k - 1 and piece k.
    let corners: Vec<([f64; 2], [f64; 2], [f64; 2])> = (0..n)
        .map(|k| {
            let (prev, next) = (ps[(k + n - 1) % n], ps[k]);
            let f = meet(offset(prev), offset(next), segs[k].start());
            (f, foot(prev, f), foot(next, f))
        })
        .collect();
    let ang = |o: [f64; 2], p: [f64; 2]| (p[1] - o[1]).atan2(p[0] - o[0]);
    let wrap = |mut d: f64, ccw: bool| {
        if ccw {
            while d < 0.0 {
                d += 2.0 * PI;
            }
        } else {
            while d > 0.0 {
                d -= 2.0 * PI;
            }
        }
        d
    };
    // ∫(x dy - y dx)/2 and the length of an arc about `o`, radius `rho`,
    // from `a` to `b` turning the signed `sweep`.
    let arc = |o: [f64; 2], rho: f64, a: [f64; 2], sweep: f64| {
        let t0 = ang(o, a);
        let t1 = t0 + sweep;
        let twice = rho * o[0] * (t1.sin() - t0.sin()) - rho * o[1] * (t1.cos() - t0.cos())
            + rho * rho * sweep;
        (twice / 2.0, rho * sweep.abs())
    };
    let (mut area, mut perimeter) = (0.0, 0.0);
    for k in 0..n {
        let (_, _, from) = corners[k];
        let (_, to, _) = corners[(k + 1) % n];
        match ps[k] {
            Piece::Line(..) => {
                area += (from[0] * to[1] - to[0] * from[1]) / 2.0;
                perimeter += len(sub(to, from));
            }
            Piece::Arc(o, rho, ccw) => {
                let (a, l) = arc(o, rho, from, wrap(ang(o, to) - ang(o, from), ccw));
                area += a;
                perimeter += l;
            }
        }
        // The corner round at the start of piece k + 1 turns anticlockwise.
        let (gc, gfrom, gto) = corners[(k + 1) % n];
        let (a, l) = arc(
            gc,
            r.abs(),
            gfrom,
            wrap(ang(gc, gto) - ang(gc, gfrom), r > 0.0),
        );
        area += a;
        perimeter += l;
    }
    (area, perimeter)
}

/// The rounded D's area and perimeter agree with the closed form.
#[test]
fn the_plan_oracle_rounds_the_d_as_the_closed_form() {
    let d = [
        Seg::Line([-4.0, 0.0], [4.0, 0.0]),
        Seg::Arc([4.0, 0.0], [0.0, 4.0], [-4.0, 0.0]),
    ];
    let (a, p) = rounded_outline(&d, 1.0);
    let (wa, wp) = rounded_d(4.0, 1.0);
    assert!(
        (a - wa).abs() < 1e-12 && (p - wp).abs() < 1e-12,
        "{a} {p} against {wa} {wp}"
    );
}

/// A plate 20 × 20 × 4 under a boss 6 high of the anticlockwise outline
/// `segs`.
fn boss_on_plate(m: &mut Model, segs: &[Seg]) -> Shape {
    let plate = prism(m, &rect(0.0, 0.0, 20.0, 20.0), 0.0, 4.0);
    let boss = prism(m, segs, 4.0, 6.0);
    let s = ogeom_bool::fuse(m, &plate, &boss, T).unwrap().shape;
    solid_of(m, &s)
}

/// A point inside each piece of `segs` at height `z`.
fn piece_points(segs: &[Seg], z: f64) -> Vec<Point> {
    segs.iter()
        .map(|s| match *s {
            Seg::Line(a, b) => Point::new(0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), z),
            Seg::Arc(_, b, _) => Point::new(b[0], b[1], z),
        })
        .collect()
}

/// The boss of `segs` on the plate, its foot loop and uprights rounded 1
/// with `Rolling`: the section within the foot round is the rounded
/// outline grown by the round's offset.
fn boss_foot_rolls(name: &str, segs: &[Seg]) {
    let mut m = Model::new();
    let s = boss_on_plate(&mut m, segs);
    let mut edges: Vec<Shape> = piece_points(segs, 4.0)
        .into_iter()
        .map(|p| edge_near(&m, &s, p))
        .collect();
    for seg in segs {
        let a = seg.start();
        edges.push(edge_near(&m, &s, Point::new(a[0], a[1], 7.0)));
    }
    let (area, perimeter) = rounded_outline(segs, 1.0);
    let want = 1600.0 + area * 6.0 + grown(perimeter, 1.0, true);
    rolls_alone(name, &m, &s, &edges, 1.0, want);
}

/// A lens boss: two arcs of radius 5 meeting at corners on y = 10.
#[test]
fn a_lens_boss_s_mixed_foot_corners_roll_round_the_uprights() {
    boss_foot_rolls(
        "lens boss foot",
        &[
            Seg::Arc([14.0, 10.0], [10.0, 12.0], [6.0, 10.0]),
            Seg::Arc([6.0, 10.0], [10.0, 8.0], [14.0, 10.0]),
        ],
    );
}

/// A square boss whose top side bows inwards: an arc of radius 5 into the
/// boss, meeting the uprights' lines at corners.
#[test]
fn a_boss_with_an_inward_arc_rolls_at_its_foot() {
    let sag = 5.0 - 21.0_f64.sqrt();
    boss_foot_rolls(
        "inward arc boss foot",
        &[
            Seg::Line([6.0, 6.0], [14.0, 6.0]),
            Seg::Line([14.0, 6.0], [14.0, 14.0]),
            Seg::Arc([14.0, 14.0], [10.0, 14.0 - 2.0 * sag], [6.0, 14.0]),
            Seg::Line([6.0, 14.0], [6.0, 6.0]),
        ],
    );
}

/// The block's pocket of `segs`, its rim and uprights rounded 1 with
/// `Rolling`: the void's section within the rim round is the rounded
/// outline grown by the round's offset.
fn pocket_rim_rolls(name: &str, segs: &[Seg]) {
    let mut m = Model::new();
    let s = pocketed(&mut m, segs);
    let mut edges: Vec<Shape> = Vec::new();
    for e in explore_unique(&m, &s, ShapeType::Edge).unwrap() {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(&m, &e, T).unwrap() else {
            continue;
        };
        let (p, q) = (c.point_at(a, T).unwrap(), c.point_at(b, T).unwrap());
        let mid = c.point_at(0.5 * (a + b), T).unwrap();
        let rim =
            (p.z - 10.0).abs() < 1e-9 && (q.z - 10.0).abs() < 1e-9 && (mid.z - 10.0).abs() < 1e-9;
        let inner = mid.x > 0.5 && mid.x < 29.5 && mid.y > 0.5 && mid.y < 29.5;
        let upright = (p.x - q.x).abs() < 1e-9 && (p.y - q.y).abs() < 1e-9 && mid.z > 5.0;
        if inner && (rim || upright) {
            edges.push(e);
        }
    }
    let (area, perimeter) = rounded_outline(segs, 1.0);
    let want = 9000.0 - area * 5.0 - grown(perimeter, 1.0, true);
    rolls_alone(name, &m, &s, &edges, 1.0, want);
}

/// A lens pocket: two arcs of radius 6 meeting at corners on y = 15.
#[test]
fn a_lens_pocket_s_mixed_rim_corners_roll_round_the_uprights() {
    let h = 6.0 - 20.0_f64.sqrt();
    pocket_rim_rolls(
        "lens pocket rim",
        &[
            Seg::Arc([19.0, 15.0], [15.0, 15.0 + 2.0 * h], [11.0, 15.0]),
            Seg::Arc([11.0, 15.0], [15.0, 15.0 - 2.0 * h], [19.0, 15.0]),
        ],
    );
}

/// A square pocket whose far side bows into the void: an arc of radius 6.
#[test]
fn a_pocket_with_an_inward_arc_rolls_at_its_rim() {
    let sag = 6.0 - 32.0_f64.sqrt();
    pocket_rim_rolls(
        "inward arc pocket rim",
        &[
            Seg::Line([11.0, 11.0], [19.0, 11.0]),
            Seg::Line([19.0, 11.0], [19.0, 19.0]),
            Seg::Arc([19.0, 19.0], [15.0, 19.0 - sag], [11.0, 19.0]),
            Seg::Line([11.0, 19.0], [11.0, 11.0]),
        ],
    );
}

/// What a round of radius `r` takes within a convex-cornered outline of
/// perimeter `p` where the section shrinks: `p ∫w - π ∫w²`, the outline
/// smooth and every concave arc of it wider than `r`.
fn shrunk(p: f64, r: f64) -> f64 {
    grown(p, r, false)
}

/// Two discs of radius 5 about (7, 10) and (13, 10), fused: a figure
/// eight whose waist corners are reflex, anticlockwise from the lower
/// waist corner.
fn figure_eight(dx: f64, dy: f64) -> Vec<Seg> {
    let w = 4.0;
    let p = |x: f64, y: f64| [x + dx, y + dy];
    vec![
        Seg::Arc(p(10.0, 10.0 - w), p(18.0, 10.0), p(10.0, 10.0 + w)),
        Seg::Arc(p(10.0, 10.0 + w), p(2.0, 10.0), p(10.0, 10.0 - w)),
    ]
}

/// The figure-eight boss's top loop and its two waist uprights rounded 1:
/// at each top corner two convex rounds meet the concave upright on the
/// boss's cylinders; the section within the top round is the rounded
/// outline shrunk by the round's offset.
#[test]
fn a_figure_eight_boss_s_top_corners_roll_round_its_waist() {
    let segs = figure_eight(0.0, 0.0);
    let mut m = Model::new();
    let s = boss_on_plate(&mut m, &segs);
    let mut edges: Vec<Shape> = piece_points(&segs, 10.0)
        .into_iter()
        .map(|p| edge_near(&m, &s, p))
        .collect();
    edges.push(edge_near(&m, &s, Point::new(10.0, 6.0, 7.0)));
    edges.push(edge_near(&m, &s, Point::new(10.0, 14.0, 7.0)));
    let (area, perimeter) = rounded_outline(&segs, -1.0);
    let want = 1600.0 + area * 6.0 - shrunk(perimeter, 1.0);
    rolls_alone("figure-eight boss top", &m, &s, &edges, 1.0, want);
}

/// The figure-eight pocket's floor loop and its two waist uprights
/// rounded 1: the void's section within the floor round is the rounded
/// outline shrunk by the round's offset.
#[test]
fn a_figure_eight_pocket_s_floor_corners_roll_round_its_waist() {
    let segs = figure_eight(5.0, 5.0);
    let mut m = Model::new();
    let s = pocketed(&mut m, &segs);
    let mut edges: Vec<Shape> = piece_points(&segs, 5.0)
        .into_iter()
        .map(|p| edge_near(&m, &s, p))
        .collect();
    edges.push(edge_near(&m, &s, Point::new(15.0, 11.0, 7.0)));
    edges.push(edge_near(&m, &s, Point::new(15.0, 19.0, 7.0)));
    let (area, perimeter) = rounded_outline(&segs, -1.0);
    let want = 9000.0 - area * 5.0 + shrunk(perimeter, 1.0);
    rolls_alone("figure-eight pocket floor", &m, &s, &edges, 1.0, want);
}

/// The lens boss rounded all round: foot, uprights and top.
#[test]
fn a_lens_boss_rounded_all_round_rolls_at_its_foot_and_top() {
    let segs = [
        Seg::Arc([14.0, 10.0], [10.0, 12.0], [6.0, 10.0]),
        Seg::Arc([6.0, 10.0], [10.0, 8.0], [14.0, 10.0]),
    ];
    let mut m = Model::new();
    let s = boss_on_plate(&mut m, &segs);
    let mut edges: Vec<Shape> = Vec::new();
    for z in [4.0, 10.0] {
        edges.extend(
            piece_points(&segs, z)
                .into_iter()
                .map(|p| edge_near(&m, &s, p)),
        );
    }
    for seg in &segs {
        let a = seg.start();
        edges.push(edge_near(&m, &s, Point::new(a[0], a[1], 7.0)));
    }
    let (area, perimeter) = rounded_outline(&segs, 1.0);
    let want = 1600.0 + area * 6.0 + grown(perimeter, 1.0, true) - shrunk(perimeter, 1.0);
    rolls_alone("lens boss all round", &m, &s, &edges, 1.0, want);
}

/// The D pocket rounded all over: floor loop, rim and uprights. The floor
/// corners are of one sense and close on the ball's sphere; the void's
/// section shrinks within the floor round and grows within the rim round.
#[test]
fn a_d_pocket_rounded_all_over_rolls_at_its_floor_and_rim() {
    let segs = [
        Seg::Line([10.0, 12.0], [20.0, 12.0]),
        Seg::Arc([20.0, 12.0], [15.0, 17.0], [10.0, 12.0]),
    ];
    let mut m = Model::new();
    let s = pocketed(&mut m, &segs);
    let mut edges = edges_of_face_at(&m, &s, Point::new(15.0, 14.0, 5.0));
    let h = 5.0 * core::f64::consts::FRAC_1_SQRT_2;
    edges.extend(
        [
            Point::new(15.0, 12.0, 10.0),
            Point::new(15.0 + h, 12.0 + h, 10.0),
            Point::new(15.0 - h, 12.0 + h, 10.0),
            Point::new(10.0, 12.0, 7.0),
            Point::new(20.0, 12.0, 7.0),
        ]
        .iter()
        .map(|&p| edge_near(&m, &s, p)),
    );
    let (area, perimeter) = rounded_outline(&segs, 1.0);
    let want = 9000.0 - area * 5.0 - grown(perimeter, 1.0, true) + shrunk(perimeter, 1.0);
    rolls_alone("D pocket all over", &m, &s, &edges, 1.0, want);
}

/// The D pocket's floor loop and uprights rounded 1: at each floor corner
/// three concave rounds meet with the pocket's cylinder and close on the
/// ball's sphere, built from the void's side (ADR 0001, amendment
/// "rolling corners on arcs, and concave corners from the void's side").
/// `Setback` and `Rolling` agree; the void's section within the floor
/// round is the rounded D shrunk by the round's offset.
#[test]
fn a_d_pocket_s_floor_and_uprights_close_on_the_ball_at_the_floor() {
    let mut m = Model::new();
    let s = d_pocket(&mut m);
    let mut edges = edges_of_face_at(&m, &s, Point::new(15.0, 14.0, 5.0));
    edges.push(edge_near(&m, &s, Point::new(10.0, 12.0, 7.0)));
    edges.push(edge_near(&m, &s, Point::new(20.0, 12.0, 7.0)));
    let segs = [
        Seg::Line([10.0, 12.0], [20.0, 12.0]),
        Seg::Arc([20.0, 12.0], [15.0, 17.0], [10.0, 12.0]),
    ];
    let (area, perimeter) = rounded_outline(&segs, 1.0);
    let want = 9000.0 - area * 5.0 + shrunk(perimeter, 1.0);
    rolls(
        "D pocket floor and uprights",
        &m,
        &s,
        &edges,
        1.0,
        Want::Volume(want),
    );
    rolls_alone("D pocket floor and uprights", &m, &s, &edges, 1.0, want);
}

/// The lens pocket rounded all over: floor, rim and uprights.
#[test]
fn a_lens_pocket_rounded_all_over_rolls_at_its_floor_and_rim() {
    let h = 6.0 - 20.0_f64.sqrt();
    let segs = [
        Seg::Arc([19.0, 15.0], [15.0, 15.0 + 2.0 * h], [11.0, 15.0]),
        Seg::Arc([11.0, 15.0], [15.0, 15.0 - 2.0 * h], [19.0, 15.0]),
    ];
    let mut m = Model::new();
    let s = pocketed(&mut m, &segs);
    let mut edges = edges_of_face_at(&m, &s, Point::new(15.0, 15.0, 5.0));
    for e in explore_unique(&m, &s, ShapeType::Edge).unwrap() {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(&m, &e, T).unwrap() else {
            continue;
        };
        let mid = c.point_at(0.5 * (a + b), T).unwrap();
        let (p, q) = (c.point_at(a, T).unwrap(), c.point_at(b, T).unwrap());
        let inner = mid.x > 0.5 && mid.x < 29.5 && mid.y > 0.5 && mid.y < 29.5;
        let rim = (mid.z - 10.0).abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9;
        let upright = (p.x - q.x).abs() < 1e-9 && (p.y - q.y).abs() < 1e-9 && mid.z > 5.0;
        if inner && (rim || upright) {
            edges.push(e);
        }
    }
    let (area, perimeter) = rounded_outline(&segs, 1.0);
    let want = 9000.0 - area * 5.0 - grown(perimeter, 1.0, true) + shrunk(perimeter, 1.0);
    rolls_alone("lens pocket all over", &m, &s, &edges, 1.0, want);
}

/// A frame at `p` about `z`.
fn frame_at(p: Point, z: Direction) -> Frame {
    let x = if z == Direction::Z {
        Direction::X
    } else {
        Direction::Z
    };
    Frame::new(p, z, x, T).unwrap()
}

/// `edges` of `s` rounded `r` with `Rolling` and `Setback`: both refuse,
/// `Rolling` naming the vertex at `at` and saying `why`.
fn refused_at(m: &Model, s: &Shape, edges: &[Shape], r: f64, at: &str, why: &str) {
    let got = fillet_edges_corner(&mut m.clone(), s, edges, r, FilletCorner::Rolling, T)
        .unwrap_err()
        .to_string();
    assert!(
        got.contains(&format!("does not close the vertex at ({at})")) && got.contains(why),
        "{got}"
    );
    assert!(fillet_edges(&mut m.clone(), s, edges, r, T).is_err());
}

/// A half cylinder of radius 4 lying on a plate, its axis on the plate's
/// top: its two end arcs (convex) and its foot (concave) rounded 1. At
/// each foot corner the lone round is the arc, whose band is a torus; the
/// ball rolling round its end has its centre on a section of the torus
/// offset by the radius, not a circle, so the sweep is no torus and the
/// corner is refused by name.
#[test]
fn a_lone_round_on_an_arc_is_refused_by_name() {
    let mut m = Model::new();
    let plate = ogeom_algo::make_box(
        &mut m,
        frame_at(Point::new(0.0, 0.0, 0.0), Direction::Z),
        (20.0, 20.0, 6.0),
        T,
    )
    .unwrap()
    .shape;
    let drum = ogeom_algo::make_cylinder(
        &mut m,
        frame_at(Point::new(5.0, 10.0, 6.0), Direction::X),
        4.0,
        10.0,
        T,
    )
    .unwrap()
    .shape;
    let fused = ogeom_bool::fuse(&mut m, &plate, &drum, T).unwrap().shape;
    let s = solid_of(&m, &fused);
    let edges: Vec<Shape> = [
        Point::new(5.0, 10.0, 10.0),
        Point::new(15.0, 10.0, 10.0),
        Point::new(5.0, 10.0, 6.0),
        Point::new(15.0, 10.0, 6.0),
        Point::new(10.0, 6.0, 6.0),
        Point::new(10.0, 14.0, 6.0),
    ]
    .iter()
    .map(|&p| edge_near(&m, &s, p))
    .collect();
    refused_at(
        &m,
        &s,
        &edges,
        1.0,
        "5.000000, 6.000000, 6.000000",
        "the lone round's edge is curved",
    );
}

/// A rib fused against a drum, lower than the drum: the rib's top edge
/// (convex) and its two edges on the drum (concave) rounded 0.5. The face
/// across from the lone round is the drum: refused by name.
#[test]
fn a_mixed_corner_on_a_curved_cap_is_refused_by_name() {
    let mut m = Model::new();
    let plate = ogeom_algo::make_box(
        &mut m,
        frame_at(Point::new(0.0, 0.0, 0.0), Direction::Z),
        (30.0, 20.0, 4.0),
        T,
    )
    .unwrap()
    .shape;
    let drum = ogeom_algo::make_cylinder(
        &mut m,
        frame_at(Point::new(10.0, 10.0, 4.0), Direction::Z),
        5.0,
        8.0,
        T,
    )
    .unwrap()
    .shape;
    let rib = ogeom_algo::make_box(
        &mut m,
        frame_at(Point::new(12.0, 9.0, 4.0), Direction::Z),
        (12.0, 2.0, 4.0),
        T,
    )
    .unwrap()
    .shape;
    let a = ogeom_bool::fuse(&mut m, &plate, &drum, T).unwrap().shape;
    let a = solid_of(&m, &a);
    let fused = ogeom_bool::fuse(&mut m, &a, &rib, T).unwrap().shape;
    let s = solid_of(&m, &fused);
    let x = 24.0_f64.sqrt() + 10.0;
    let edges: Vec<Shape> = [
        Point::new(18.0, 11.0, 8.0),
        Point::new(x, 11.0, 6.0),
        Point::new(x - 0.05, 10.5, 8.0),
    ]
    .iter()
    .map(|&p| edge_near(&m, &s, p))
    .collect();
    refused_at(
        &m,
        &s,
        &edges,
        0.5,
        "14.898979, 11.000000, 8.000000",
        "the face across from the lone round is curved",
    );
}
