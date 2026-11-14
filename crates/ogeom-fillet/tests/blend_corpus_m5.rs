//! The acceptance corpus for roll-over and non-local blends (ADR 0001, M5,
//! §12; §13 question 8): rounds whose contact leaves its face roll onto the
//! edge it crosses, the ball resting on the edge until the contact comes
//! back, its centre a radius off the other face and a radius from the edge.
//! Each case is built through the engine with no boolean, valid, the same
//! whatever the request's order, and measured against the exact rolling
//! envelope, sampled from the closed form of each ball, never from the
//! builder's curves.
//!
//! The cases run only with `OGEOM_M5_CORPUS` set, as M3's corpus does until
//! its classes route through the engine (ADR 0001, §11); without it each
//! returns at once. The round side's cases are #973's and #974's; the
//! chamfer side adds its own.
#![allow(clippy::unwrap_used, reason = "test constructions")]

use ogeom_core::Tolerances;
use ogeom_geom::{CircleCurve, Curve, LineCurve, PlaneSurface};
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// Whether the corpus is asked for: only where `OGEOM_M5_CORPUS` is set.
fn asked() -> bool {
    std::env::var_os("OGEOM_M5_CORPUS").is_some()
}

fn up_at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn point_of(model: &Model, v: &Shape) -> Point {
    v.transform(model.datums())
        .unwrap()
        .apply(model.node(v).unwrap().data().as_vertex().unwrap().point)
}

/// The edges of `solid` both of whose ends `keep`.
fn edges_where(model: &Model, solid: &Shape, keep: impl Fn(Point) -> bool) -> Vec<Shape> {
    explore_unique(model, solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            ogeom_algo::edge_vertices(model, e)
                .unwrap()
                .is_some_and(|(a, b)| keep(point_of(model, &a)) && keep(point_of(model, &b)))
        })
        .collect()
}

/// The arc of the ball about `centre` of `radius` from `a` to `b`, at
/// `n + 1` points.
fn arc(centre: Point, a: Point, b: Point, radius: f64, n: u32) -> Vec<Point> {
    let (p, q) = ((a - centre) / radius, (b - centre) / radius);
    let axis = p.cross(q) / p.cross(q).magnitude();
    let sweep = p.dot(q).clamp(-1.0, 1.0).acos();
    (0..=n)
        .map(|k| {
            let t = sweep * f64::from(k) / f64::from(n);
            centre + (p * t.cos() + axis.cross(p) * t.sin()) * radius
        })
        .collect()
}

/// How far `p` stands from the nearest face of `shape` whose trim holds
/// its foot: each face asked directly, as a foot within a parameter's
/// epsilon of a trim's edge (both sides of a change) is what a scan of the
/// whole shape can miss.
fn off(model: &Model, shape: &Shape, p: Point) -> f64 {
    use ogeom_geom::Transformable as _;
    let mut best = f64::INFINITY;
    for f in explore_unique(model, shape, ShapeType::Face).unwrap() {
        let data = model.node(&f).unwrap().data().as_face().unwrap();
        let surface = model
            .geometry()
            .surface(data.surface)
            .unwrap()
            .clone()
            .transformed(&f.transform(model.datums()).unwrap(), T)
            .unwrap();
        let foot = ogeom_algo::project_on_surface(&surface, p, 32, T).unwrap();
        if foot.distance >= best {
            continue;
        }
        let fine = ogeom_mesh::Deflection {
            chord: T.confusion(),
            ..Default::default()
        };
        let trim = ogeom_algo::classify::FaceTrim::of(model, &f, fine, T).unwrap();
        if trim
            .holds(model, foot.parameters, T.confusion(), T)
            .unwrap()
            .is_inside_or_on()
        {
            best = foot.distance;
        }
    }
    best
}

/// Rounds `edges` of `solid` at `radius`, asserting the engine built it with
/// no boolean and the solid is valid.
fn rounded_by_engine(model: &mut Model, solid: &Shape, edges: &[Shape], radius: f64) -> Shape {
    let (built, work) = ogeom_bool::work::counted(|| {
        ogeom_fillet::fillet_edges(model, solid, edges, radius, T).unwrap()
    });
    assert_eq!(work.face_pairs, 0, "built with no boolean");
    assert_eq!(work.seated_cuts, 0, "nothing trimmed");
    let check = ogeom_algo::check(model, &built.shape, T).unwrap();
    assert!(check.is_valid(), "{:?}", check.problems);
    built.shape
}

/// #974's shelf: a cylinder of radius 10, 26 tall, fused with a block 12.6
/// by 10 by 21 beside it, the block's top front edge rounded at 3. The top
/// is narrower than the round where the cylinder's rim cuts it (|y| < 2.8):
/// there the ball rests on the rim, its centre (9.6, y, 21 - h), h the
/// height that holds it 3 from the rim's point (96 / rho, 10 y / rho, 21),
/// rho = |(9.6, y)|. Either side it is the plain round, centre (9.6, y,
/// 18). Every ball's arc from the front to its other touch lies on the
/// result within 2e-6, both sides of each change included.
#[test]
fn a_shelf_s_top_edge_rolls_on_the_cylinder_s_rim() {
    if !asked() {
        return;
    }
    let mut model = Model::new();
    let cylinder = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 10.0, 26.0, T)
        .unwrap()
        .shape;
    let block = ogeom_algo::make_box(&mut model, up_at(0.0, -5.0, 0.0), (12.6, 10.0, 21.0), T)
        .unwrap()
        .shape;
    let solid = ogeom_bool::fuse(&mut model, &cylinder, &block, T)
        .unwrap()
        .shape;
    let edges = edges_where(&model, &solid, |p| {
        (p.x - 12.6).abs() < 1e-7 && (p.z - 21.0).abs() < 1e-7
    });
    assert_eq!(edges.len(), 1);
    let shape = rounded_by_engine(&mut model, &solid, &edges, 3.0);
    for y in [
        -5.0_f64, -2.80001, -2.8, -2.79999, -1.3, 0.0, 1.3, 2.79999, 2.8, 2.80001, 5.0,
    ] {
        let rho = 9.6_f64.hypot(y);
        let over = (10.0 - rho).max(0.0);
        let centre = Point::new(9.6, y, 21.0 - (9.0 - over * over).sqrt());
        let top = if over > 0.0 {
            Point::new(96.0 / rho, 10.0 * y / rho, 21.0)
        } else {
            Point::new(9.6, y, 21.0)
        };
        let front = Point::new(12.6, y, centre.z);
        for p in arc(centre, front, top, 3.0, 16) {
            let d = off(&model, &shape, p);
            assert!(d < 2e-6, "y {y}: {p:?} is {d} off the result");
        }
    }
}

/// #973's boss foot: a stadium boss (two half circles of radius 0.6 at
/// y = ±5 joined by straight sides at x = ±0.6) fused on a plate 20 by 12,
/// its foot rounded at 1. Round each half circle the ball's centre stands
/// 1.6 off its axis; past the plate's edge (|y| = 6) it rests on that edge,
/// at the height that holds it 1 from it. Every ball's arc from the boss to
/// its other touch lies on the result within 1e-5, at the changes too.
#[test]
fn a_stadium_boss_s_foot_rolls_on_the_plate_s_edges() {
    if !asked() {
        return;
    }
    let mut model = Model::new();
    let corners = [
        Point::new(0.6, -5.0, 0.0),
        Point::new(0.6, 5.0, 0.0),
        Point::new(-0.6, 5.0, 0.0),
        Point::new(-0.6, -5.0, 0.0),
    ];
    let vertices: Vec<Shape> = corners
        .iter()
        .map(|&p| ogeom_algo::make_vertex(&mut model, p).shape)
        .collect();
    let mut profile = Vec::new();
    for i in 0..4 {
        let (curve, range) = if i % 2 == 0 {
            let line = LineCurve::segment(corners[i], corners[(i + 1) % 4], T).unwrap();
            (Curve::Line(line), (0.0, 10.0))
        } else {
            let y = if i == 1 { 5.0 } else { -5.0 };
            let lo = if i == 1 { 0.0 } else { core::f64::consts::PI };
            let circle = Circle::new(up_at(0.0, y, 0.0), 0.6, T).unwrap();
            (
                Curve::Circle(CircleCurve::new(circle)),
                (lo, lo + core::f64::consts::PI),
            )
        };
        profile.push(
            ogeom_algo::make_edge_between(
                &mut model,
                curve,
                range,
                &vertices[i],
                &vertices[(i + 1) % 4],
                T,
            )
            .unwrap()
            .shape,
        );
    }
    let wire = ogeom_algo::make_wire(&mut model, &profile, T)
        .unwrap()
        .shape;
    let face = ogeom_algo::make_face(
        &mut model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let boss = ogeom_algo::make_prism(&mut model, &face, Vector::new(0.0, 0.0, 10.0), T)
        .unwrap()
        .shape;
    let plate = ogeom_algo::make_box(&mut model, up_at(-10.0, -6.0, -2.0), (20.0, 12.0, 2.0), T)
        .unwrap()
        .shape;
    let solid = ogeom_bool::fuse(&mut model, &plate, &boss, T)
        .unwrap()
        .shape;
    let edges = edges_where(&model, &solid, |p| {
        p.z.abs() < T.confusion() && p.x.abs() < 1.0 && p.y.abs() < 5.7
    });
    assert_eq!(edges.len(), 4);
    // The request the other way round gives the same solid.
    let mut other = model.clone();
    let mut reversed = edges.clone();
    reversed.reverse();
    let turned = rounded_by_engine(&mut other, &solid, &reversed, 1.0);
    let shape = rounded_by_engine(&mut model, &solid, &edges, 1.0);
    let volume = |m: &Model, s: &Shape| {
        ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
            .unwrap()
            .mass
    };
    let faces = |m: &Model, s: &Shape| explore_unique(m, s, ShapeType::Face).unwrap().len();
    assert_eq!(
        faces(&model, &shape),
        faces(&other, &turned),
        "faces either order"
    );
    assert!(
        (volume(&model, &shape) - volume(&other, &turned)).abs() <= 1e-9,
        "volume either order"
    );
    let near = ogeom_algo::PointDistance::of(&model, &shape, T).unwrap();
    let switch = (1.0_f64 / 1.6).asin();
    for theta in [
        0.0,
        switch - 1e-6,
        switch,
        switch + 1e-6,
        core::f64::consts::FRAC_PI_2,
        core::f64::consts::PI - switch,
        core::f64::consts::PI,
    ] {
        let (sin, cos) = theta.sin_cos();
        // Past the edge by `over` along y, the ball rests on it.
        let over = (1.6 * sin - 1.0).max(0.0);
        let height = (1.0 - over * over).sqrt();
        for end in [1.0_f64, -1.0] {
            let centre = Point::new(1.6 * cos, end * (5.0 + 1.6 * sin), height);
            let wall = Point::new(0.6 * cos, end * (5.0 + 0.6 * sin), height);
            let base = if over > 0.0 {
                Point::new(1.6 * cos, end * 6.0, 0.0)
            } else {
                Point::new(1.6 * cos, end * (5.0 + 1.6 * sin), 0.0)
            };
            for p in arc(centre, wall, base, 1.0, 16) {
                let d = near.nearest(p, T).unwrap().0;
                assert!(
                    d < 1e-5,
                    "theta {theta} end {end}: {p:?} is {d} off the result"
                );
            }
        }
    }
}
