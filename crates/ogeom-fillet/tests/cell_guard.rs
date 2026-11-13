//! A setback corner cuts a cell out round its vertex and sews back only
//! the corner's own patches and trims: anything else of the solid in that
//! cell (here a 1 × 1 × 1 void) would be lost from a solid that still
//! checks valid. Each corner that cuts a cell refuses it by name instead.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{Blend, Chamfer, ChamferCorner, Setbacks};
use ogeom_geom::{Continuity, Curve3d as _};
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn boxed(m: &mut Model, at: [f64; 3], size: [f64; 3]) -> Shape {
    let frame = Frame::new(
        Point::new(at[0], at[1], at[2]),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    ogeom_algo::make_box(m, frame, (size[0], size[1], size[2]), T)
        .unwrap()
        .shape
}

/// `solid` less a 1 × 1 × 1 void at `at`, as a second shell.
fn voided(m: &mut Model, solid: &Shape, at: [f64; 3]) -> Shape {
    voided_sized(m, solid, at, 1.0)
}

/// `solid` less a cube void of side `side` at `at`, as a second shell.
fn voided_sized(m: &mut Model, solid: &Shape, at: [f64; 3], side: f64) -> Shape {
    let void = boxed(m, at, [side, side, side]);
    let holed = ogeom_bool::cut(m, solid, &void, T).unwrap().shape;
    assert_eq!(
        explore_unique(m, &holed, ShapeType::Shell).unwrap().len(),
        2,
        "the void is a shell of its own"
    );
    holed
}

fn vertex_at(m: &Model, s: &Shape, p: [f64; 3]) -> Shape {
    let q = Point::new(p[0], p[1], p[2]);
    explore_unique(m, s, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .find(|v| {
            m.node(v)
                .and_then(|n| n.data().as_vertex().map(|d| d.point))
                .is_some_and(|x| x.distance(q) < 1e-9)
        })
        .expect("a vertex there")
}

/// The edge of `s` passing nearest `p`.
fn edge_near(m: &Model, s: &Shape, p: [f64; 3]) -> Shape {
    let q = Point::new(p[0], p[1], p[2]);
    let d = |e: &Shape| {
        let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
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

fn even(setback: f64) -> Setbacks {
    Setbacks {
        per_edge: Vec::new(),
        default: setback,
        continuity: Continuity::G1,
    }
}

fn refused(what: &str, got: ogeom_core::OgeomResult<ogeom_algo::Built>, m: &Model) {
    match got {
        Ok(b) => panic!(
            "{what}: built, {} shells, the void lost",
            explore_unique(m, &b.shape, ShapeType::Shell).unwrap().len()
        ),
        Err(e) => assert!(
            e.to_string()
                .contains("another feature of the solid lies inside the setback corner's cell"),
            "{what}: {e}"
        ),
    }
}

/// The 20 × 20 × 10 block with the void at (2, 2, 6), inside the cell of
/// its top corner at (0, 0, 10) set back 5.
fn block_with_void(m: &mut Model) -> Shape {
    let block = boxed(m, [0.0, 0.0, 0.0], [20.0, 20.0, 10.0]);
    voided(m, &block, [2.0, 2.0, 6.0])
}

/// The corner's three edges at (0, 0, 10): along x, along y, down z.
fn top_corner_edges(m: &Model, s: &Shape) -> [Shape; 3] {
    [
        edge_near(m, s, [10.0, 0.0, 10.0]),
        edge_near(m, s, [0.0, 10.0, 10.0]),
        edge_near(m, s, [0.0, 0.0, 5.0]),
    ]
}

#[test]
fn a_setback_corner_with_a_void_in_its_cell_is_refused_by_name() {
    let mut m = Model::new();
    let s = block_with_void(&mut m);
    let v = vertex_at(&m, &s, [0.0, 0.0, 10.0]);
    let got = ogeom_fillet::round_vertex_setback(&mut m, &s, &v, 2.0, &even(5.0), T);
    refused("round_vertex_setback", got, &m);
}

#[test]
fn setback_corners_through_fillet_edges_with_a_void_in_a_cell_are_refused_by_name() {
    let mut m = Model::new();
    let s = block_with_void(&mut m);
    let v = vertex_at(&m, &s, [0.0, 0.0, 10.0]);
    let block_edges: Vec<Shape> = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let b = ogeom_algo::shape_bounds(&m, e, T).unwrap();
            let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
            // The block's own edges, not the void's.
            [lo.x, lo.y, lo.z, hi.x, hi.y, hi.z]
                .iter()
                .all(|c| c.abs() < 1e-6 || (c - 20.0).abs() < 1e-6 || (c - 10.0).abs() < 1e-6)
        })
        .collect();
    assert_eq!(block_edges.len(), 12);
    let got =
        ogeom_fillet::fillet_edges_cornered(&mut m, &s, &block_edges, 2.0, &[(v, even(5.0))], T);
    refused("fillet_edges_cornered", got, &m);
}

/// The void of the three-radii corner: a cube of 0.4 at (1.7, 1.7, 8.0),
/// inside its smaller cell, clear of the bands and of the balls' centres.
const SPOT: [f64; 4] = [1.7, 1.7, 8.0, 0.4];

#[test]
fn a_corner_of_three_radii_with_a_void_in_its_cell_is_refused_by_name() {
    let mut m = Model::new();
    let block = boxed(&mut m, [0.0, 0.0, 0.0], [20.0, 20.0, 10.0]);
    let s = voided_sized(&mut m, &block, [SPOT[0], SPOT[1], SPOT[2]], SPOT[3]);
    let [x, y, z] = top_corner_edges(&m, &s);
    let got = ogeom_fillet::fillet_edges_with(&mut m, &s, &[(x, 1.0), (y, 1.5), (z, 1.25)], T);
    refused("fillet_edges_with", got, &m);
}

/// The box corner at the origin of `blend_edges`' own tests, with the void
/// at (2.5, 2.5, 2.5) inside the corner's cell, clear of every band and of
/// the ball's centre.
fn origin_corner_with_void(m: &mut Model) -> (Shape, [Shape; 3]) {
    let block = boxed(m, [0.0, 0.0, 0.0], [20.0, 20.0, 10.0]);
    let s = voided(m, &block, [2.5, 2.5, 2.5]);
    let e = [
        edge_near(m, &s, [10.0, 0.0, 0.0]),
        edge_near(m, &s, [0.0, 10.0, 0.0]),
        edge_near(m, &s, [0.0, 0.0, 5.0]),
    ];
    (s, e)
}

/// `request`'s bevels and rounds on the voided block's corner through
/// `blend_edges`: built, valid, the void kept as a shell of its own, and
/// the same solid as the bevels' own call and the rounds' own call on what
/// it leaves (ADR 0001, section 13, question 3: two features, so no corner
/// cell reaches the void, which stands clear of every band).
fn the_bevels_then_the_rounds_keep_the_void(what: &str, kinds: &[(usize, bool)]) {
    let blend_of = |round: bool| {
        if round {
            Blend::Round(2.0)
        } else {
            Blend::Chamfer(Chamfer::Symmetric(2.0))
        }
    };
    let mut m = Model::new();
    let (s, e) = origin_corner_with_void(&mut m);
    let request: Vec<(Shape, Blend)> = kinds
        .iter()
        .map(|(k, round)| (e[*k].clone(), blend_of(*round)))
        .collect();
    let got = ogeom_fillet::blend_edges(&mut m, &s, &request, ChamferCorner::Point, T)
        .unwrap_or_else(|why| panic!("{what}: {why}"));
    assert!(
        ogeom_algo::check(&m, &got.shape, T).unwrap().is_valid(),
        "{what}"
    );
    assert_eq!(
        explore_unique(&m, &got.shape, ShapeType::Shell)
            .unwrap()
            .len(),
        2,
        "{what}: the void kept"
    );
    let mut n = Model::new();
    let (s2, e2) = origin_corner_with_void(&mut n);
    let bevels: Vec<(Shape, Chamfer)> = kinds
        .iter()
        .filter(|(_, round)| !*round)
        .map(|(k, _)| (e2[*k].clone(), Chamfer::Symmetric(2.0)))
        .collect();
    let beveled =
        ogeom_fillet::chamfer_edges_cornered(&mut n, &s2, &bevels, ChamferCorner::Point, T)
            .unwrap();
    let rounds: Vec<(Shape, f64)> = kinds
        .iter()
        .filter(|(_, round)| *round)
        .map(|(k, _)| (beveled.history.trace(&e2[*k])[0].clone(), 2.0))
        .collect();
    let separate = ogeom_fillet::fillet_edges_with(&mut n, &beveled.shape, &rounds, T).unwrap();
    let volume = |m: &Model, s: &Shape| {
        ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
            .unwrap()
            .mass
    };
    assert_eq!(
        volume(&m, &got.shape).to_bits(),
        volume(&n, &separate.shape).to_bits(),
        "{what}: the two features in turn"
    );
}

#[test]
fn a_corner_of_rounds_and_bevels_with_a_void_near_it_is_the_bevels_then_the_rounds() {
    the_bevels_then_the_rounds_keep_the_void(
        "a bevel and two rounds",
        &[(0, false), (1, true), (2, true)],
    );
}

#[test]
fn a_round_and_a_bevel_beside_a_sharp_edge_with_a_void_near_them_are_the_bevel_then_the_round() {
    the_bevels_then_the_rounds_keep_the_void(
        "a round and a bevel beside a sharp edge",
        &[(0, false), (1, true)],
    );
}

#[test]
fn a_setback_corner_with_an_edge_left_sharp_and_a_void_in_its_cell_is_refused_by_name() {
    let mut m = Model::new();
    let s = block_with_void(&mut m);
    let v = vertex_at(&m, &s, [0.0, 0.0, 10.0]);
    let [_, _, down] = top_corner_edges(&m, &s);
    let got = ogeom_fillet::round_vertex_setback_sharp(&mut m, &s, &v, &down, 2.0, &even(5.0), T);
    refused("round_vertex_setback_sharp", got, &m);
}

#[test]
fn a_convex_and_a_concave_round_with_a_void_in_their_cell_are_refused_by_name() {
    let mut m = Model::new();
    let base = boxed(&mut m, [0.0, 0.0, 0.0], [20.0, 20.0, 5.0]);
    let upright = boxed(&mut m, [0.0, 0.0, 0.0], [5.0, 20.0, 20.0]);
    let bracket = ogeom_bool::fuse(&mut m, &base, &upright, T).unwrap().shape;
    let s = voided(&mut m, &bracket, [2.5, 1.0, 2.5]);
    let edges = vec![
        edge_near(&m, &s, [12.5, 0.0, 5.0]),
        edge_near(&m, &s, [5.0, 10.0, 5.0]),
        edge_near(&m, &s, [5.0, 0.0, 12.5]),
    ];
    // All three edges rounded, and the upright's front edge left sharp,
    // in either order.
    for asked in [
        edges.clone(),
        edges[..2].to_vec(),
        vec![edges[1].clone(), edges[0].clone()],
    ] {
        let mut own = m.clone();
        let got = ogeom_fillet::fillet_edges(&mut own, &s, &asked, 2.0, T);
        refused(
            &format!("fillet_edges at a signed corner, {} edges", asked.len()),
            got,
            &own,
        );
    }
}
