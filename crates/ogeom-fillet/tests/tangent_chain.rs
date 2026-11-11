//! Rounds where a flat wall runs tangent into a curved one, the face
//! between them narrowing to a cusp at the tangent point.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
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

/// A 40 x 26 x 5 plate, x from -20 to 20 and y from -6 to 20, fused with a
/// boss of radius 6 about the z axis, 10 high: the plate's side at y = -6
/// runs tangent into the boss's wall at (0, -6).
fn tangent_boss(m: &mut Model) -> Shape {
    let plate = ogeom_algo::make_box(m, at(-20.0, -6.0, 0.0), (40.0, 26.0, 5.0), T)
        .unwrap()
        .shape;
    let boss = ogeom_algo::make_cylinder(m, at(0.0, 0.0, 0.0), 6.0, 10.0, T)
        .unwrap()
        .shape;
    ogeom_bool::fuse(m, &plate, &boss, T).unwrap().shape
}

/// How far past the plate's side `y = -6` the result stands above the
/// plate's top, outside the boss: read on every face's mesh, each point
/// classified against the side's plane. Nought where nothing hangs over.
fn overhang(m: &Model, s: &Shape) -> f64 {
    let mut worst = 0.0_f64;
    for face in explore_unique(m, s, ShapeType::Face).unwrap() {
        let mesh = ogeom_mesh::triangulate_face(
            m,
            &face,
            ogeom_mesh::Deflection {
                chord: 1e-3,
                ..ogeom_mesh::Deflection::default()
            },
            T,
        )
        .unwrap();
        for p in &mesh.positions {
            if p.z >= 5.0 - 1e-9 && p.x.hypot(p.y) > 6.0 + 1e-6 {
                worst = worst.max(-6.0 - p.y);
            }
        }
    }
    worst
}

/// `fillet_edges` on `asked` (points on the edges) of the tangent boss at
/// radius 1, in both orders: refused by name as a concave band laying its
/// fill past a face, never built hanging past the plate's side.
fn refused_past_the_side(asked: &[Point]) {
    for flip in [false, true] {
        let mut m = Model::new();
        let s = tangent_boss(&mut m);
        let mut edges: Vec<Shape> = asked.iter().map(|p| edge_near(&m, &s, *p)).collect();
        if flip {
            edges.reverse();
        }
        match ogeom_fillet::fillet_edges(&mut m, &s, &edges, 1.0, T) {
            Ok(b) => panic!(
                "flip {flip}: built, standing {} past the plate's side",
                overhang(&m, &b.shape)
            ),
            Err(why) => assert!(
                why.to_string()
                    .contains("would lay its fill past that face's side over nothing"),
                "flip {flip}: {why}"
            ),
        }
    }
}

/// The boss's foot alone, rounded at 1: within about 3.5 of (0, -6, 5) the
/// plate's top between the foot and its side narrows below the band's
/// width, so the foot's concave band would hang past the side over nothing
/// (before, built so, a full radius out at (0, -7, 5)). Refused by name,
/// the whole foot or its half beside the cusp.
#[test]
fn a_foot_round_on_a_face_narrowing_to_a_cusp_is_refused_not_hung_past_it() {
    let left = Point::new(-4.242_640_687, 4.242_640_687, 5.0);
    let right = Point::new(4.242_640_687, 4.242_640_687, 5.0);
    refused_past_the_side(&[left]);
    refused_past_the_side(&[left, right]);
}

/// The foot on the side x < 0 (concave) and the plate's top edge on the
/// side x > 0 (convex), a tangent chain running on through the cusp at
/// (0, -6, 5): the foot's band hangs past the side as it does alone, and
/// the request is refused by that name in either order.
#[test]
fn a_tangent_chain_through_the_cusp_is_refused_by_the_foot_s_overhang() {
    refused_past_the_side(&[
        Point::new(-4.242_640_687, 4.242_640_687, 5.0),
        Point::new(10.0, -6.0, 5.0),
    ]);
}

/// printCAD's corpus case: a 40 x 40 x 5 plate under a boss of radius 6 at
/// (14, 0), 8 high, its wall tangent to the plate's side at x = 20, the
/// foot and the plate's top edge at x = 20 rounded at 1 as one tangent
/// chain. Either side of (20, 0, 5) the plate's top narrows to a cusp
/// between the edge and the foot, and a sound round there must turn the
/// edge's round into the foot's through it. Before 0.18.0 the foot's band
/// ran on round the boss past the plate's side, nothing under it (material
/// at (20.95, 0, 5.001), outside the solid); the request is refused by
/// name in either order until that turn is built.
#[test]
fn a_foot_and_an_edge_through_the_cusp_where_a_boss_runs_tangent_into_a_plate_are_refused() {
    for flip in [false, true] {
        let mut m = Model::new();
        let plate = ogeom_algo::make_box(&mut m, at(-20.0, -20.0, 0.0), (40.0, 40.0, 5.0), T)
            .unwrap()
            .shape;
        let boss = ogeom_algo::make_cylinder(&mut m, at(14.0, 0.0, 5.0), 6.0, 8.0, T)
            .unwrap()
            .shape;
        let s = ogeom_bool::fuse(&mut m, &plate, &boss, T).unwrap().shape;
        let mut chain = Vec::new();
        for e in explore_unique(&m, &s, ShapeType::Edge).unwrap() {
            let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(&m, &e, T).unwrap() else {
                continue;
            };
            let (p, q) = (c.point_at(a, T).unwrap(), c.point_at(b, T).unwrap());
            let mid = c.point_at(0.5 * (a + b), T).unwrap();
            let on_top = [p, q, mid].iter().all(|x| (x.z - 5.0).abs() < 1e-6);
            let foot = ((mid.x - 14.0).hypot(mid.y) - 6.0).abs() < 1e-6;
            let side = [p, q, mid].iter().all(|x| (x.x - 20.0).abs() < 1e-6);
            if on_top && (foot || side) {
                chain.push(e);
            }
        }
        assert_eq!(
            chain.len(),
            4,
            "two arcs of the foot, two halves of the edge"
        );
        if flip {
            chain.reverse();
        }
        let refused = ogeom_fillet::fillet_edges(&mut m, &s, &chain, 1.0, T);
        let Err(why) = refused else {
            panic!("flip {flip}: built");
        };
        let why = why.to_string();
        assert!(
            why.contains("close as one corner only where")
                || why.contains("would lay its fill past that face's side over nothing"),
            "flip {flip}: {why}"
        );
    }
}

/// The issue's arm and hub: an arm x 0..40, y -6..6, z 0..5 fused with a
/// hub of radius 6 about z, 10 high, the arm's sides tangent to the hub's
/// wall. Its top side edge (convex) and the arc where its top meets the hub
/// (concave) leave (0, 6, 5) the same way, the arm's top between them a
/// cusp, so the two bands would share it: refused by name in either order.
#[test]
fn the_arm_s_side_edge_and_its_arc_into_the_hub_are_refused_in_either_order() {
    for flip in [false, true] {
        let mut m = Model::new();
        let arm = ogeom_algo::make_box(&mut m, at(0.0, -6.0, 0.0), (40.0, 12.0, 5.0), T)
            .unwrap()
            .shape;
        let hub = ogeom_algo::make_cylinder(&mut m, at(0.0, 0.0, 0.0), 6.0, 10.0, T)
            .unwrap()
            .shape;
        let s = ogeom_bool::fuse(&mut m, &arm, &hub, T).unwrap().shape;
        let mut edges = vec![
            edge_near(&m, &s, Point::new(20.0, 6.0, 5.0)),
            edge_near(&m, &s, Point::new(4.242_640_687, 4.242_640_687, 5.0)),
        ];
        if flip {
            edges.reverse();
        }
        let Err(why) = ogeom_fillet::fillet_edges(&mut m, &s, &edges, 1.0, T) else {
            panic!("flip {flip}: built");
        };
        let why = why.to_string();
        assert!(
            why.contains("close as one corner only where")
                || why.contains("would lay its fill past that face's side over nothing"),
            "flip {flip}: {why}"
        );
    }
}
