//! A convex and a concave round meeting at one vertex: the L-bracket's
//! front edge at its re-entrant edge, (5, 0, 5).
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn edge_near(m: &Model, s: &Shape, q: Point) -> Shape {
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

/// The L-bracket: a base 20 x 20 x 5 and an upright 5 x 20 x 20.
fn bracket(m: &mut Model) -> Shape {
    let bx = |m: &mut Model, s: (f64, f64, f64)| {
        ogeom_algo::make_box(
            m,
            Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap(),
            s,
            T,
        )
        .unwrap()
        .shape
    };
    let a = bx(m, (20.0, 20.0, 5.0));
    let b = bx(m, (5.0, 20.0, 20.0));
    ogeom_bool::fuse(m, &a, &b, T).unwrap().shape
}

/// The bracket's three edges at (5, 0, 5): the front edge along the base's
/// top (convex), the re-entrant edge (concave) and the upright's front
/// edge (convex).
fn corner_edges(m: &Model, l: &Shape) -> [Shape; 3] {
    [
        edge_near(m, l, Point::new(12.5, 0.0, 5.0)),
        edge_near(m, l, Point::new(5.0, 10.0, 5.0)),
        edge_near(m, l, Point::new(5.0, 0.0, 12.5)),
    ]
}

/// The issue's request: the front edge and the re-entrant edge rounded,
/// the upright's front edge left sharp. The patch that closes the opening
/// through the sharp edge's point, tangent to every band end and face,
/// stands 0.10 to 0.21 into the concave round's own ball (about a tenth of
/// the radius), where neither the bracket nor any round of the request
/// puts material: refused by name, in either order, at every radius.
#[test]
fn a_convex_and_a_concave_round_beside_an_edge_left_sharp_are_refused_where_the_patch_leaves_the_solid()
 {
    for radius in [2.0, 1.0, 1.5] {
        for flip in [false, true] {
            let mut m = Model::new();
            let l = bracket(&mut m);
            let [front, crease, _] = corner_edges(&m, &l);
            let edges = if flip {
                vec![crease, front]
            } else {
                vec![front, crease]
            };
            let why = ogeom_fillet::fillet_edges(&mut m, &l, &edges, radius, T)
                .err()
                .map(|e| e.to_string())
                .unwrap_or_default();
            assert!(
                why.contains("outside what the solid and its concave rounds hold"),
                "radius {radius} flip {flip}: {why}"
            );
        }
    }
}

/// All three edges at the vertex rounded: the one patch over the opening,
/// tangent to every band end and face, stands 0.19 to 0.46 into the
/// concave round's own ball, where neither the bracket nor any round of
/// the request puts material: refused by name, whatever the order.
#[test]
fn three_rounds_of_both_kinds_at_a_vertex_are_refused_where_the_patch_leaves_the_solid() {
    let orders = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    for order in orders {
        let mut m = Model::new();
        let l = bracket(&mut m);
        let edges = corner_edges(&m, &l);
        let asked: Vec<Shape> = order.iter().map(|&k| edges[k].clone()).collect();
        let why = ogeom_fillet::fillet_edges(&mut m, &l, &asked, 2.0, T)
            .err()
            .map(|e| e.to_string())
            .unwrap_or_default();
        assert!(
            why.contains("outside what the solid and its concave rounds hold"),
            "{order:?}: {why}"
        );
    }
}

/// Bands as wide as the upright is thick: the cell the corner cuts reaches
/// twice the widest band behind the vertex, through the upright's far
/// face, which nothing would put back. Refused by name, with all three
/// edges rounded or the upright's front edge left sharp, as the cell
/// guard refuses any other feature inside the cell.
#[test]
fn a_signed_corner_whose_cell_reaches_through_the_wall_is_refused_by_name() {
    for three in [false, true] {
        let mut m = Model::new();
        let l = bracket(&mut m);
        let [front, crease, upright] = corner_edges(&m, &l);
        let mut edges = vec![front, crease];
        if three {
            edges.push(upright);
        }
        let err = ogeom_fillet::fillet_edges(&mut m, &l, &edges, 3.0, T)
            .map(|_| ())
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("another feature of the solid lies inside the setback corner's cell"),
            "three {three}: {err}"
        );
    }
}
