//! Rounds and bevels in one request, closed as one corner where they meet.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{Blend, Chamfer, ChamferCorner, blend_edges};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn vol(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// A 20 × 20 × 10 box and its edges at the origin along x, y and z.
fn corner(m: &mut Model) -> (Shape, [Shape; 3]) {
    let frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
    let b = ogeom_algo::make_box(m, frame, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let along = |axis: usize| {
        explore_unique(m, &b, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let (c, (s, t)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
                let (p, q) = (c.point_at(s, T).unwrap(), c.point_at(t, T).unwrap());
                let at = |p: Point| p.distance(Point::ORIGIN) < 1e-9;
                let other = if at(p) {
                    q
                } else if at(q) {
                    p
                } else {
                    return false;
                };
                let v = [other.x, other.y, other.z];
                (0..3).all(|i| (i == axis) == (v[i].abs() > 1e-9))
            })
            .unwrap()
    };
    let e = [along(0), along(1), along(2)];
    (b, e)
}

/// The request on the box corner, `kinds[i]` for the edge along axis `i`:
/// a bevel `Symmetric(2)`, a round of 2, or none; in `order`.
fn request(e: &[Shape; 3], kinds: [usize; 3], order: [usize; 3]) -> Vec<(Shape, Blend)> {
    order
        .iter()
        .filter(|&&i| kinds[i] != 2)
        .map(|&i| {
            let blend = if kinds[i] == 0 {
                Blend::Chamfer(Chamfer::Symmetric(2.0))
            } else {
                Blend::Round(2.0)
            };
            (e[i].clone(), blend)
        })
        .collect()
}

/// Every vertex's point, sorted, as bits.
fn points(m: &Model, s: &Shape) -> Vec<[u64; 3]> {
    let mut out: Vec<[u64; 3]> = explore_unique(m, s, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| {
            let p = m.node(v).unwrap().data().as_vertex().unwrap().point;
            [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
        })
        .collect();
    out.sort_unstable();
    out
}

/// The corner `kinds` asks, in every order: valid, and the same solid to
/// the bit as the bevels' own call and then the rounds' own call on what it
/// leaves (ADR 0001, section 13, question 3: two features, the bevels
/// first).
fn is_the_bevels_then_the_rounds(kinds: [usize; 3]) {
    let mut s = Model::new();
    let (sb, se) = corner(&mut s);
    let bevels: Vec<(Shape, Chamfer)> = (0..3)
        .filter(|&i| kinds[i] == 0)
        .map(|i| (se[i].clone(), Chamfer::Symmetric(2.0)))
        .collect();
    let beveled =
        ogeom_fillet::chamfer_edges_cornered(&mut s, &sb, &bevels, ChamferCorner::Point, T)
            .unwrap();
    let rounds: Vec<(Shape, f64)> = (0..3)
        .filter(|&i| kinds[i] == 1)
        .map(|i| (beveled.history.trace(&se[i])[0].clone(), 2.0))
        .collect();
    let separate = ogeom_fillet::fillet_edges_with(&mut s, &beveled.shape, &rounds, T)
        .unwrap_or_else(|e| panic!("{kinds:?}: the rounds on the bevels: {e}"));
    let want = (
        vol(&s, &separate.shape).to_bits(),
        points(&s, &separate.shape),
    );
    for order in [
        [0usize, 1, 2],
        [2, 1, 0],
        [1, 0, 2],
        [0, 2, 1],
        [1, 2, 0],
        [2, 0, 1],
    ] {
        let mut m = Model::new();
        let (b, e) = corner(&mut m);
        let r = blend_edges(
            &mut m,
            &b,
            &request(&e, kinds, order),
            ChamferCorner::Point,
            T,
        )
        .unwrap_or_else(|why| panic!("{kinds:?} {order:?}: {why}"));
        let found = ogeom_algo::check(&m, &r.shape, T).unwrap();
        assert!(
            found.is_valid(),
            "{kinds:?} {order:?}: {:?}",
            found.problems
        );
        let key = (vol(&m, &r.shape).to_bits(), points(&m, &r.shape));
        assert_eq!(key, want, "{kinds:?}: order {order:?}");
    }
}

/// A bevel `Symmetric(2)` and two rounds of 2 meeting at a box's corner.
#[test]
fn a_bevel_and_two_rounds_meeting_are_the_bevel_then_the_rounds_whatever_the_order() {
    is_the_bevels_then_the_rounds([0, 1, 1]);
}

/// Two bevels and a round meeting at a box's corner.
#[test]
fn two_bevels_and_a_round_meeting_are_the_bevels_then_the_round_whatever_the_order() {
    is_the_bevels_then_the_rounds([0, 0, 1]);
}

/// A round meeting a bevel at a box's corner, the third edge left sharp,
/// each edge taking each part in turn.
#[test]
fn a_round_meeting_a_bevel_beside_a_sharp_edge_is_the_bevel_then_the_round() {
    for kinds in [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        is_the_bevels_then_the_rounds(kinds);
    }
}

#[test]
fn a_request_of_one_kind_is_that_kind_s_own_call() {
    let mut m = Model::new();
    let (b, e) = corner(&mut m);
    let mut plain = m.clone();
    let rounds: Vec<(Shape, f64)> = e.iter().map(|x| (x.clone(), 2.0)).collect();
    let own = ogeom_fillet::fillet_edges_with(&mut plain, &b, &rounds, T).unwrap();
    let blended = blend_edges(
        &mut m,
        &b,
        &request(&e, [1, 1, 1], [0, 1, 2]),
        ChamferCorner::Point,
        T,
    )
    .unwrap();
    assert_eq!(
        vol(&plain, &own.shape).to_bits(),
        vol(&m, &blended.shape).to_bits()
    );
    assert_eq!(points(&plain, &own.shape), points(&m, &blended.shape));

    let mut m = Model::new();
    let (b, e) = corner(&mut m);
    let mut plain = m.clone();
    let bevels: Vec<(Shape, Chamfer)> = e
        .iter()
        .map(|x| (x.clone(), Chamfer::Symmetric(2.0)))
        .collect();
    let own =
        ogeom_fillet::chamfer_edges_cornered(&mut plain, &b, &bevels, ChamferCorner::Facet, T)
            .unwrap();
    let blended = blend_edges(
        &mut m,
        &b,
        &request(&e, [0, 0, 0], [0, 1, 2]),
        ChamferCorner::Facet,
        T,
    )
    .unwrap();
    assert_eq!(
        vol(&plain, &own.shape).to_bits(),
        vol(&m, &blended.shape).to_bits()
    );
    assert_eq!(points(&plain, &own.shape), points(&m, &blended.shape));
}

#[test]
fn a_round_and_a_bevel_that_do_not_meet_are_their_own_calls_in_turn() {
    // The bevel along x at the origin and a round on the top edge along x
    // at the far side: no vertex shared, so the bevel and then the round,
    // each as its own call builds it.
    let edge_at = |m: &Model, b: &Shape, p: Point| {
        explore_unique(m, b, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let (c, (s, t)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
                c.point_at(f64::midpoint(s, t), T).unwrap().distance(p) < 1e-9
            })
            .unwrap()
    };
    let mut m = Model::new();
    let (b, e) = corner(&mut m);
    let far = edge_at(&m, &b, Point::new(10.0, 20.0, 10.0));
    let mut plain = m.clone();
    let beveled = ogeom_fillet::chamfer_edges_with(
        &mut plain,
        &b,
        &[(e[0].clone(), Chamfer::Symmetric(2.0))],
        T,
    )
    .unwrap();
    let traced = beveled.history.trace(&far)[0].clone();
    let own =
        ogeom_fillet::fillet_edges_with(&mut plain, &beveled.shape, &[(traced, 2.0)], T).unwrap();
    let blended = blend_edges(
        &mut m,
        &b,
        &[
            (far.clone(), Blend::Round(2.0)),
            (e[0].clone(), Blend::Chamfer(Chamfer::Symmetric(2.0))),
        ],
        ChamferCorner::Point,
        T,
    )
    .unwrap();
    assert_eq!(
        vol(&plain, &own.shape).to_bits(),
        vol(&m, &blended.shape).to_bits()
    );
    assert_eq!(points(&plain, &own.shape), points(&m, &blended.shape));
}
