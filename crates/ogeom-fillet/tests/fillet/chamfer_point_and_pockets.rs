//! `ChamferCorner::Point` where no three bevels meet, and concave bevels
//! (ADR 0001, amendment "chamfers the engine refused on 0.19.0"): requests
//! with no vertex of three bevels asked with `ChamferCorner::Point` build
//! exactly what the plain call builds, and concave bevels round a pocket
//! floor's tangent loop and on straight edges between planes build, to the
//! section swept along each edge (Pappus round each arc).
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom_core::Tolerances;
use ogeom_fillet::{Chamfer, ChamferCorner, chamfer_edges_checked, chamfer_edges_cornered};
use ogeom_geom::{CircleCurve, Curve, Curve3d as _, LineCurve};
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

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

fn frame_at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

/// The only solid of `s`.
fn solid_of(m: &Model, s: &Shape) -> Shape {
    let solids = explore_unique(m, s, ShapeType::Solid).unwrap();
    assert_eq!(solids.len(), 1);
    solids[0].clone()
}

/// The edge of `s` passing nearest `q`, 64 samples along each.
fn edge_near(m: &Model, s: &Shape, q: Point) -> Shape {
    let d = |e: &Shape| {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, e, T).unwrap() else {
            return f64::INFINITY;
        };
        (0..=64)
            .map(|i| {
                c.point_at(a + (b - a) * f64::from(i) / 64.0, T)
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

fn symmetric(edges: &[Shape], d: f64) -> Vec<(Shape, Chamfer)> {
    edges
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(d)))
        .collect()
}

/// `specs` on `s` by the plain call: it builds, sound; its volume.
fn plain(name: &str, m: &Model, s: &Shape, specs: &[(Shape, Chamfer)]) -> f64 {
    let mut m = m.clone();
    let built = chamfer_edges_checked(&mut m, s, specs, T)
        .unwrap_or_else(|e| panic!("{name}: the plain call refused: {}", e.error));
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{name}, plain: {:?}", found.problems);
    volume(&m, &built.shape)
}

/// `specs` on `s` with `Point` corners: it builds, sound; its volume.
fn pointed(name: &str, m: &Model, s: &Shape, specs: &[(Shape, Chamfer)]) -> f64 {
    let mut m = m.clone();
    let built = chamfer_edges_cornered(&mut m, s, specs, ChamferCorner::Point, T)
        .unwrap_or_else(|e| panic!("{name}: Point refused: {e}"));
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{name}, Point: {:?}", found.problems);
    volume(&m, &built.shape)
}

/// Both calls build `specs`, sound, to `want` within a millionth, and
/// `Point` to exactly the plain call's volume.
fn both_build(name: &str, m: &Model, s: &Shape, specs: &[(Shape, Chamfer)], want: f64) {
    let by_plain = plain(name, m, s, specs);
    assert!(
        (by_plain - want).abs() <= 1e-6 * want,
        "{name}, plain: volume {by_plain} against {want}"
    );
    let by_point = pointed(name, m, s, specs);
    assert!(
        by_point.to_bits() == by_plain.to_bits(),
        "{name}: Point gives {by_point}, the plain call {by_plain}"
    );
}

/// The volume the triangle `a b c` of a section in the (radius, height)
/// half plane sweeps about the axis, by Pappus.
fn pappus(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    let area = 0.5 * ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs();
    area * 2.0 * PI * (a.0 + b.0 + c.0) / 3.0
}

fn plate(m: &mut Model) -> Shape {
    ogeom_algo::make_box(m, frame_at(-20.0, -20.0, 0.0), (40.0, 40.0, 5.0), T)
        .unwrap()
        .shape
}

/// A 40 × 40 × 5 plate under a cone of radius 6 to 2, 8 high, its axis
/// through the origin.
fn plate_and_cone(m: &mut Model) -> Shape {
    let base = plate(m);
    let cone = ogeom_algo::make_cone(m, frame_at(0.0, 0.0, 5.0), 6.0, 2.0, 8.0, T)
        .unwrap()
        .shape;
    let fused = ogeom_bool::fuse(m, &base, &cone, T).unwrap().shape;
    solid_of(m, &fused)
}

fn plate_and_cone_volume() -> f64 {
    8000.0 + PI * 8.0 / 3.0 * (36.0 + 12.0 + 4.0)
}

/// A point `s` along the cone's wall from `(r, z)`, up (`+1`) or down.
fn on_cone_wall((r, z): (f64, f64), s: f64) -> (f64, f64) {
    let k = 80.0_f64.sqrt();
    (r - 4.0 * s / k, z + 8.0 * s / k)
}

#[test]
fn a_cone_foot_rim_chamfers_with_point_corners_as_the_plain_call() {
    let mut m = Model::new();
    let s = plate_and_cone(&mut m);
    let rim = edge_near(&m, &s, Point::new(6.0, 0.0, 5.0));
    let d = 0.5;
    let foot = (6.0, 5.0);
    let filled = pappus(foot, (6.0 + d, 5.0), on_cone_wall(foot, d));
    let want = plate_and_cone_volume() + filled;
    both_build("cone foot rim", &m, &s, &symmetric(&[rim], d), want);
}

#[test]
fn a_cone_top_rim_chamfers_with_point_corners_as_the_plain_call() {
    let mut m = Model::new();
    let s = plate_and_cone(&mut m);
    let rim = edge_near(&m, &s, Point::new(2.0, 0.0, 13.0));
    let d = 0.5;
    let top = (2.0, 13.0);
    let taken = pappus(top, (2.0 - d, 13.0), on_cone_wall(top, -d));
    let want = plate_and_cone_volume() - taken;
    both_build("cone top rim", &m, &s, &symmetric(&[rim], d), want);
}

#[test]
fn a_boss_foot_rim_chamfers_with_point_corners_as_the_plain_call() {
    let mut m = Model::new();
    let base = plate(&mut m);
    let boss = ogeom_algo::make_cylinder(&mut m, frame_at(0.0, 0.0, 5.0), 6.0, 8.0, T)
        .unwrap()
        .shape;
    let fused = ogeom_bool::fuse(&mut m, &base, &boss, T).unwrap().shape;
    let s = solid_of(&m, &fused);
    let rim = edge_near(&m, &s, Point::new(6.0, 0.0, 5.0));
    let d = 1.5;
    let want = 8000.0 + PI * 36.0 * 8.0 + pappus((6.0, 5.0), (6.0 + d, 5.0), (6.0, 5.0 + d));
    both_build("boss foot rim", &m, &s, &symmetric(&[rim], d), want);
}

/// A 20 × 20 × 10 block with a 20 × 5 × 1 step fused on its top along
/// y 0..5: the fuse leaves seams where the step's sides run on from the
/// block's.
fn seamed_block(m: &mut Model) -> Shape {
    let block = ogeom_algo::make_box(m, frame_at(0.0, 0.0, 0.0), (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let step = ogeom_algo::make_box(m, frame_at(0.0, 0.0, 10.0), (20.0, 5.0, 1.0), T)
        .unwrap()
        .shape;
    let fused = ogeom_bool::fuse(m, &block, &step, T).unwrap().shape;
    solid_of(m, &fused)
}

/// Three convex edges of the seamed block, each alone and together, bevel
/// by `Point` as by the plain call, to the closed form.
#[test]
fn edges_of_a_seamed_block_chamfer_with_point_corners_as_the_plain_call() {
    for d in [0.5, 2.0] {
        let mut m = Model::new();
        let s = seamed_block(&mut m);
        let top_back = edge_near(&m, &s, Point::new(10.0, 20.0, 10.0));
        let step_back = edge_near(&m, &s, Point::new(10.0, 5.0, 11.0));
        let upright = edge_near(&m, &s, Point::new(20.0, 20.0, 5.0));
        let section = 0.5 * d * d;
        // The step's back face is 1 high: a bevel of 2 on its top edge
        // would cut below the floor beyond it.
        let step_fits = d <= 1.0;
        let mut cases: Vec<(&str, Vec<Shape>, f64)> = vec![
            ("top back", vec![top_back.clone()], section * 20.0),
            ("upright", vec![upright.clone()], section * 10.0),
            (
                "top back and upright",
                vec![top_back.clone(), upright.clone()],
                // The two meet at (20, 20, 10), sharing d³/3.
                section * 30.0 - d * d * d / 3.0,
            ),
        ];
        if !step_fits {
            // Both calls refuse it alike (ADR 0001, amendment "a bevel stops
            // at a concave floor").
            let specs = symmetric(std::slice::from_ref(&step_back), d);
            let plain = chamfer_edges_checked(&mut m.clone(), &s, &specs, T)
                .expect_err("the plain call refuses the step at 2")
                .error;
            let point = chamfer_edges_cornered(&mut m.clone(), &s, &specs, ChamferCorner::Point, T)
                .expect_err("Point refuses the step at 2");
            assert_eq!(plain.to_string(), point.to_string());
        } else {
            cases.push(("step back", vec![step_back.clone()], section * 20.0));
            cases.push((
                "all three",
                vec![top_back, step_back, upright],
                section * 50.0 - d * d * d / 3.0,
            ));
        }
        for (name, edges, taken) in cases {
            let name = format!("seamed block, {name}, {d}");
            both_build(&name, &m, &s, &symmetric(&edges, d), 4100.0 - taken);
        }
    }
}

/// The seamed block's corner (20, 20, 10), its three edges bevelled 0.5:
/// `Point` mitres it on the block read with its seams dissolved, the three
/// planes cutting the corner (pairwise overlaps `d³/3`, all three `d³/4`),
/// and the plain call's facet takes `d³/12` more.
#[test]
fn a_seamed_block_s_corner_of_three_mitres_with_point_corners() {
    let mut m = Model::new();
    let s = seamed_block(&mut m);
    let edges: Vec<Shape> = [(10.0, 20.0, 10.0), (20.0, 20.0, 5.0), (20.0, 12.0, 10.0)]
        .iter()
        .map(|&(x, y, z)| edge_near(&m, &s, Point::new(x, y, z)))
        .collect();
    let d = 0.5;
    let specs = symmetric(&edges, d);
    let mitred = 4100.0 - (0.5 * d * d * (20.0 + 10.0 + 15.0) - d * d * d + d * d * d / 4.0);
    let by_point = pointed("seamed corner", &m, &s, &specs);
    assert!(
        (by_point - mitred).abs() <= 1e-6 * mitred,
        "seamed corner, Point: volume {by_point} against {mitred}"
    );
    let by_plain = plain("seamed corner", &m, &s, &specs);
    let faceted = mitred - d * d * d / 12.0;
    assert!(
        (by_plain - faceted).abs() <= 1e-6 * faceted,
        "seamed corner, plain: volume {by_plain} against {faceted}"
    );
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

/// A block 30 × 30 × 10 with the profile `pocket` cut 5 deep into its top.
fn pocketed(m: &mut Model, pocket: &[Seg]) -> Shape {
    let block = ogeom_algo::make_box(m, frame_at(0.0, 0.0, 0.0), (30.0, 30.0, 10.0), T)
        .unwrap()
        .shape;
    let tool = prism(m, pocket, 5.0, 6.0);
    let cut = ogeom_bool::cut(m, &block, &tool, T).unwrap().shape;
    solid_of(m, &cut)
}

/// The obround pocket's floor loop, two lines and two half circles of
/// radius 3 tangent all round, bevelled 1 (concave): the pocket plus the
/// section `d²/2` along the lines and, by Pappus, round the arcs with its
/// centroid `d/3` in from the wall.
#[test]
fn an_obround_pocket_s_floor_loop_chamfers() {
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
    let d = 1.0;
    let base = 9000.0 - (14.0 * 6.0 + PI * 9.0) * 5.0;
    let want = base + 0.5 * d * d * (28.0 + 2.0 * PI * (3.0 - d / 3.0));
    both_build("obround floor", &m, &s, &symmetric(&floor, d), want);
}

/// The kidney pocket's floor loop, four tangent arcs (radii 12 and 6
/// about (15, 8) over 120 degrees, ends of radius 3), bevelled 1: the
/// centroid `d/3` inside the outer arc and the ends, outside the inner.
#[test]
fn a_kidney_pocket_s_floor_loop_chamfers() {
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
    let d = 1.0;
    let base = 9000.0 - (PI / 3.0 * (144.0 - 36.0) + PI * 9.0) * 5.0;
    let third = 2.0 * PI / 3.0;
    let filled = 0.5
        * d
        * d
        * (third * (12.0 - d / 3.0) + third * (6.0 + d / 3.0) + 2.0 * PI * (3.0 - d / 3.0));
    both_build("kidney floor", &m, &s, &symmetric(&floor, d), base + filled);
}

/// Two flanges 4 thick, 30 high and 10 deep, 36 apart, and a bar 6 high
/// between them rising 13 across: the four concave straight edges where
/// the bar meets the flanges, bevelled 1, each filling `sin θ / 2` of
/// section along its 10, `θ` the angle its faces make.
#[test]
fn a_slanted_bar_s_concave_edges_chamfer() {
    let mut m = Model::new();
    let flange = |m: &mut Model, x: f64| {
        ogeom_algo::make_box(m, frame_at(x, 0.0, 0.0), (4.0, 10.0, 30.0), T)
            .unwrap()
            .shape
    };
    let left = flange(&mut m, 0.0);
    let right = flange(&mut m, 40.0);
    // The bar's profile in the plane y = 0, run 2 into each flange, swept
    // 10 along y.
    let z = |x: f64| 5.0 + 13.0 * (x - 4.0) / 36.0;
    let outline = [
        Point::new(2.0, 0.0, z(2.0)),
        Point::new(42.0, 0.0, z(42.0)),
        Point::new(42.0, 0.0, z(42.0) + 6.0),
        Point::new(2.0, 0.0, z(2.0) + 6.0),
    ];
    let wire = ogeom_algo::make_polygon(&mut m, &outline, true, T)
        .unwrap()
        .shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Y),
        (-1000.0, 1000.0),
        (-1000.0, 1000.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(&mut m, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let bar = ogeom_algo::make_prism(&mut m, &face, Vector::new(0.0, 10.0, 0.0), T)
        .unwrap()
        .shape;
    let bar = ogeom_algo::baked_shape(&mut m, &bar, T).unwrap().shape;
    let one = ogeom_bool::fuse(&mut m, &left, &bar, T).unwrap().shape;
    let both = ogeom_bool::fuse(&mut m, &one, &right, T).unwrap().shape;
    let s = solid_of(&m, &both);
    let base = 2400.0 + 36.0 * 6.0 * 10.0;
    let edges: Vec<Shape> = [
        (4.0, z(4.0)),
        (4.0, z(4.0) + 6.0),
        (40.0, z(40.0)),
        (40.0, z(40.0) + 6.0),
    ]
    .iter()
    .map(|&(x, h)| edge_near(&m, &s, Point::new(x, 5.0, h)))
    .collect();
    let a = (13.0_f64 / 36.0).atan();
    let d = 1.0;
    let want = base + 2.0 * 10.0 * 0.5 * d * d * ((PI / 2.0 + a).sin() + (PI / 2.0 - a).sin());
    both_build("slanted bar", &m, &s, &symmetric(&edges, d), want);
}
