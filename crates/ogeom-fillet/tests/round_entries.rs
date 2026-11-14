//! What a round request may be asked on, through the blend engine: a solid
//! with voids (each shell rounded as a solid of its own, a round reaching a
//! void refused by name), a compound of solids (each rounded alone), a
//! solid with placed shapes (rounded baked, its history carried back); and
//! #974's two negatives refused by the engine itself, by name.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::FRAC_PI_4;

use ogeom_core::Tolerances;
use ogeom_geom::{LineCurve, PlaneSurface};
use ogeom_math::{Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn up_at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn point_of(model: &Model, v: &Shape) -> Point {
    v.transform(model.datums())
        .unwrap()
        .apply(model.node(v).unwrap().data().as_vertex().unwrap().point)
}

/// The edges of `shape` both of whose ends `keep`.
fn edges_where(model: &Model, shape: &Shape, keep: impl Fn(Point) -> bool) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            ogeom_algo::edge_vertices(model, e)
                .unwrap()
                .is_some_and(|(a, b)| keep(point_of(model, &a)) && keep(point_of(model, &b)))
        })
        .collect()
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom_algo::volume_properties(model, shape, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// #974's shelf (a cylinder of radius 10, 26 tall, fused with a block 12.6
/// by 10 by 21 beside it), its block's top front edge at `x = 12.6`, `z =
/// 21`.
fn shelf(model: &mut Model) -> Shape {
    let cylinder = ogeom_algo::make_cylinder(model, Frame::WORLD, 10.0, 26.0, T)
        .unwrap()
        .shape;
    let block = ogeom_algo::make_box(model, up_at(0.0, -5.0, 0.0), (12.6, 10.0, 21.0), T)
        .unwrap()
        .shape;
    ogeom_bool::fuse(model, &cylinder, &block, T).unwrap().shape
}

fn shelf_edge(model: &Model, solid: &Shape) -> Vec<Shape> {
    edges_where(model, solid, |p| {
        (p.x - 12.6).abs() < 1e-7 && (p.z - 21.0).abs() < 1e-7
    })
}

/// #974's shelf with a small spherical cavity just under its rounded edge,
/// where the round's ball sweeps: the round reaches into the void, refused
/// by the engine by name (not handed to another path).
#[test]
fn a_round_reaching_into_a_void_is_refused_by_name() {
    let mut model = Model::new();
    let solid = shelf(&mut model);
    let hole = ogeom_algo::make_sphere(&mut model, up_at(11.7, 0.2, 20.0), 0.3, T)
        .unwrap()
        .shape;
    let drilled = ogeom_bool::cut(&mut model, &solid, &hole, T).unwrap();
    let edges: Vec<Shape> = shelf_edge(&model, &solid)
        .iter()
        .map(|e| drilled.history.trace(e)[0].clone())
        .collect();
    let refusal = ogeom_fillet::fillet_edges(&mut model, &drilled.shape, &edges, 3.0, T)
        .expect_err("refused");
    assert!(
        refusal.to_string().contains("a round reaches into a void"),
        "{refusal}"
    );
}

/// #974's other negative: the cylinder above the shelf cut away past `y =
/// 1`, so its rim ends there while the shelf runs on. The ball resting on
/// the rim would roll past the rim's end, a vertex seat the engine refuses
/// by name.
#[test]
fn a_round_whose_ball_would_roll_past_its_rim_s_end_is_refused_by_name() {
    let mut model = Model::new();
    let cylinder = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 10.0, 26.0, T)
        .unwrap()
        .shape;
    let tool = ogeom_algo::make_box(&mut model, up_at(-11.0, 1.0, 21.0), (22.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let cylinder = ogeom_bool::cut(&mut model, &cylinder, &tool, T)
        .unwrap()
        .shape;
    let block = ogeom_algo::make_box(&mut model, up_at(0.0, -5.0, 0.0), (12.6, 10.0, 21.0), T)
        .unwrap()
        .shape;
    let solid = ogeom_bool::fuse(&mut model, &cylinder, &block, T)
        .unwrap()
        .shape;
    let (unified, _) = ogeom_heal::unify_same_domain(&mut model, &solid, T).unwrap();
    let solid = unified.shape;
    let edges = shelf_edge(&model, &solid);
    assert_eq!(edges.len(), 1);
    let refusal =
        ogeom_fillet::fillet_edges(&mut model, &solid, &edges, 3.0, T).expect_err("refused");
    assert!(
        refusal.to_string().contains("roll past a vertex"),
        "{refusal}"
    );
}

/// A 10 by 8 by 6 box with a spherical void of radius 1 at its middle, a
/// top edge rounded at 1: the outer shell rounds as a solid alone and the
/// void comes back as it was. The box less the void less the spandrel
/// `(1 - pi/4) r^2` along the 10 long edge, with no boolean.
#[test]
fn a_solid_with_a_void_the_round_does_not_reach_keeps_it() {
    let mut model = Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 8.0, 6.0), T)
        .unwrap()
        .shape;
    let ball = ogeom_algo::make_sphere(&mut model, up_at(5.0, 4.0, 3.0), 1.0, T)
        .unwrap()
        .shape;
    let voided = ogeom_bool::cut(&mut model, &block, &ball, T).unwrap();
    let shape = voided.shape;
    assert_eq!(model.children_of(&shape).unwrap().len(), 2, "a void");
    let edges = edges_where(&model, &shape, |p| {
        p.y.abs() < 1e-9 && (p.z - 6.0).abs() < 1e-9
    });
    assert_eq!(edges.len(), 1);
    let before = volume(&model, &shape);
    let (built, work) = ogeom_bool::work::counted(|| {
        ogeom_fillet::fillet_edges(&mut model, &shape, &edges, 1.0, T).unwrap()
    });
    assert_eq!(work.face_pairs, 0, "built with no boolean");
    let check = ogeom_algo::check(&model, &built.shape, T).unwrap();
    assert!(check.is_valid(), "{:?}", check.problems);
    assert_eq!(
        model.children_of(&built.shape).unwrap().len(),
        2,
        "the void kept"
    );
    let want = before - (1.0 - FRAC_PI_4) * 10.0;
    let got = volume(&model, &built.shape);
    assert!((got - want).abs() < 1e-6, "{got} against {want}");
}

/// Two boxes side by side in one compound, an edge of each rounded at 1:
/// each box rounded alone, the compound gathered again, each the box less
/// its edge's spandrel, with no boolean.
#[test]
fn a_compound_s_solids_round_each_alone() {
    let mut model = Model::new();
    let a = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 8.0, 6.0), T)
        .unwrap()
        .shape;
    let b = ogeom_algo::make_box(&mut model, up_at(20.0, 0.0, 0.0), (4.0, 5.0, 6.0), T)
        .unwrap()
        .shape;
    let both = ogeom_algo::make_compound(&mut model, &[a.clone(), b.clone()])
        .unwrap()
        .shape;
    let mut edges = edges_where(&model, &a, |p| p.y.abs() < 1e-9 && (p.z - 6.0).abs() < 1e-9);
    edges.extend(edges_where(&model, &b, |p| {
        (p.x - 20.0).abs() < 1e-9 && (p.z - 6.0).abs() < 1e-9
    }));
    assert_eq!(edges.len(), 2);
    let (built, work) = ogeom_bool::work::counted(|| {
        ogeom_fillet::fillet_edges(&mut model, &both, &edges, 1.0, T).unwrap()
    });
    assert_eq!(work.face_pairs, 0, "built with no boolean");
    let parts = model.children_of(&built.shape).unwrap();
    assert_eq!(parts.len(), 2);
    let spandrel = 1.0 - FRAC_PI_4;
    let mut volumes: Vec<f64> = parts.iter().map(|p| volume(&model, p)).collect();
    volumes.sort_by(f64::total_cmp);
    let want = [120.0 - spandrel * 5.0, 480.0 - spandrel * 10.0];
    for (got, want) in volumes.iter().zip(want) {
        assert!((got - want).abs() < 1e-6, "{got} against {want}");
    }
    for p in &parts {
        assert!(ogeom_algo::check(&model, p, T).unwrap().is_valid());
    }
}

/// A prism whose far cap is its profile placed (the same nodes standing
/// twice): a 6 by 4 rectangle swept 5 up, an edge of its far cap rounded
/// at 1. Rounded with its placements baked in, with no boolean: the prism
/// less the spandrel along the 6 long edge, and the far cap's history
/// carried back to the faces asked.
#[test]
fn a_prism_s_placed_far_cap_rounds_baked() {
    let mut model = Model::new();
    let corners = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(6.0, 0.0, 0.0),
        Point::new(6.0, 4.0, 0.0),
        Point::new(0.0, 4.0, 0.0),
    ];
    let vertices: Vec<Shape> = corners
        .iter()
        .map(|&p| ogeom_algo::make_vertex(&mut model, p).shape)
        .collect();
    let mut sides = Vec::new();
    for i in 0..4 {
        let (p, q) = (corners[i], corners[(i + 1) % 4]);
        let line = LineCurve::segment(p, q, T).unwrap();
        let range = ogeom_geom::Curve3d::domain(&line);
        sides.push(
            ogeom_algo::make_edge_between(
                &mut model,
                line.into(),
                range,
                &vertices[i],
                &vertices[(i + 1) % 4],
                T,
            )
            .unwrap()
            .shape,
        );
    }
    let wire = ogeom_algo::make_wire(&mut model, &sides, T).unwrap().shape;
    let face = ogeom_algo::make_face(
        &mut model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let prism = ogeom_algo::make_prism(&mut model, &face, Vector::new(0.0, 0.0, 5.0), T)
        .unwrap()
        .shape;
    let edges = edges_where(&model, &prism, |p| {
        p.y.abs() < 1e-9 && (p.z - 5.0).abs() < 1e-9
    });
    assert_eq!(edges.len(), 1);
    let (built, work) = ogeom_bool::work::counted(|| {
        ogeom_fillet::fillet_edges(&mut model, &prism, &edges, 1.0, T).unwrap()
    });
    assert_eq!(work.face_pairs, 0, "built with no boolean");
    let check = ogeom_algo::check(&model, &built.shape, T).unwrap();
    assert!(check.is_valid(), "{:?}", check.problems);
    let want = 120.0 - (1.0 - FRAC_PI_4) * 6.0;
    let got = volume(&model, &built.shape);
    assert!((got - want).abs() < 1e-6, "{got} against {want}");
    // The rounded edge generated the strip, through the baking.
    assert!(
        !built.history.generated(&edges[0]).is_empty(),
        "the strip from its edge"
    );
}
