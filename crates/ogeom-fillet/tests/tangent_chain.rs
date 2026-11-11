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

fn volume(m: &Model, s: &Shape) -> f64 {
    let deflection = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    ogeom_algo::volume_properties(m, s, deflection, T)
        .unwrap()
        .mass
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

/// Every edge of the boss's foot, the circle of radius 6 at `z = 5`.
fn foot_arcs(m: &Model, s: &Shape) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, e, T).unwrap() else {
                return false;
            };
            [a, 0.5 * (a + b), b].iter().all(|t| {
                let p = c.point_at(*t, T).unwrap();
                (p.z - 5.0).abs() < 1e-9 && (p.x.hypot(p.y) - 6.0).abs() < 1e-9
            })
        })
        .collect()
}

/// The volume of the foot's concave band of radius 1 about the boss's axis
/// standing past the plane `y = -6`, over the whole turn: the band's
/// section (the corner square less the quarter disc, of height
/// `1 - sqrt(1 - (7 - rho)^2)` at `rho`) swept over the angles where
/// `rho sin(phi) < -6`, `pi - 2 asin(6 / rho)` of them, at radius `rho`.
fn overhang_volume() -> f64 {
    let n = 2_000_000_u32;
    let mut sum = 0.0;
    for k in 0..n {
        let rho = 6.0 + (f64::from(k) + 0.5) / f64::from(n);
        let height = 1.0 - (1.0 - (7.0 - rho) * (7.0 - rho)).sqrt();
        sum += rho * height * (std::f64::consts::PI - 2.0 * (6.0 / rho).asin());
    }
    sum / f64::from(n)
}

/// The boss's foot alone, rounded at 1. Within about 3.5 of (0, -6, 5) the
/// plate's top between the foot and its side narrows below the band's
/// width, so the band runs past the side; it is cut off at the side's
/// plane, its fill ending flush with the side. In either order, the half
/// foot beside the cusp and the whole foot: valid, no point past the side,
/// and the volume the plate and boss plus the band (Pappus) less its part
/// past the side.
#[test]
fn a_foot_round_on_a_face_narrowing_to_a_cusp_is_cut_off_at_the_face_s_side() {
    use std::f64::consts::PI;
    let left = Point::new(-4.242_640_687, 4.242_640_687, 5.0);
    let area = 1.0 - PI / 4.0;
    let offset = (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI);
    let base = 40.0 * 26.0 * 5.0 + PI * 36.0 * 5.0;
    let past = overhang_volume();
    for (whole, turn, share) in [(false, PI, 0.5), (true, 2.0 * PI, 1.0)] {
        let want = base + area * (6.0 + offset) * turn - past * share;
        let mut first: Option<f64> = None;
        for flip in [false, true] {
            let mut m = Model::new();
            let s = tangent_boss(&mut m);
            let mut edges = if whole {
                foot_arcs(&m, &s)
            } else {
                vec![edge_near(&m, &s, left)]
            };
            if flip {
                edges.reverse();
            }
            let rounded = ogeom_fillet::fillet_edges(&mut m, &s, &edges, 1.0, T)
                .unwrap_or_else(|e| panic!("whole {whole}, flip {flip}: {e}"))
                .shape;
            let found = ogeom_algo::check(&m, &rounded, T).unwrap();
            assert!(found.is_valid(), "{:?}", found.problems);
            let out = overhang(&m, &rounded);
            let got = volume(&m, &rounded);
            eprintln!(
                "whole {whole}, flip {flip}: {out} past the side, volume {got} against {want}"
            );
            assert!(out <= 1e-6, "{out} past the plate's side");
            assert!((got - want).abs() < 2e-4, "{got} against {want}");
            // The foot's arcs are rounded in the order asked, as before the
            // cut; the two orders agree on the solid's volume.
            match first {
                None => first = Some(got),
                Some(held) => assert!(
                    (held - got).abs() < 1e-9 * got,
                    "flip {flip}: {held} against {got}"
                ),
            }
        }
    }
}

/// The foot on the side x < 0 (concave) and the plate's top edge on the
/// side x > 0 (convex), a tangent chain running on through the cusp at
/// (0, -6, 5). The two bands lie either side of the vertex; the foot's is
/// cut off at the plate's side as it is alone. Both orders give the same
/// solid, valid, every point of it inside the plate and boss or in the
/// foot's fill over the plate's top, and the volume the plate and boss
/// plus the foot's band less its part past the side, less the edge's band
/// over its 20.
#[test]
fn a_tangent_chain_through_the_cusp_rounds_either_side_of_it_in_either_order() {
    use std::f64::consts::PI;
    let r = 1.0_f64;
    let area = r * r * (1.0 - PI / 4.0);
    let offset = r * (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI);
    let want = 40.0 * 26.0 * 5.0 + PI * 36.0 * 5.0 + area * PI * (6.0 + offset)
        - 0.5 * overhang_volume()
        - area * 20.0;
    let mut first: Option<(u64, Vec<[u64; 3]>)> = None;
    for flip in [false, true] {
        let mut m = Model::new();
        let s = tangent_boss(&mut m);
        let foot = edge_near(&m, &s, Point::new(-4.242_640_687, 4.242_640_687, 5.0));
        let side = edge_near(&m, &s, Point::new(10.0, -6.0, 5.0));
        let asked = if flip {
            vec![side, foot]
        } else {
            vec![foot, side]
        };
        let rounded = ogeom_fillet::fillet_edges(&mut m, &s, &asked, r, T)
            .unwrap_or_else(|e| panic!("flip {flip}: {e}"))
            .shape;
        let found = ogeom_algo::check(&m, &rounded, T).unwrap();
        assert!(found.is_valid(), "flip {flip}: {:?}", found.problems);
        let got = volume(&m, &rounded);
        eprintln!("flip {flip}: volume {got} against {want}");
        assert!(
            (got - want).abs() < 2e-4,
            "flip {flip}: {got} against {want}"
        );
        let inside = |p: Point| -> bool {
            let e = 1e-6;
            let plate = p.x.abs() <= 20.0 + e
                && p.y >= -6.0 - e
                && p.y <= 20.0 + e
                && p.z >= -e
                && p.z <= 5.0 + e;
            let rho = p.x.hypot(p.y);
            let boss = rho <= 6.0 + e && p.z >= -e && p.z <= 10.0 + e;
            let fill = p.x <= e
                && p.y >= -6.0 - e
                && p.z >= 5.0 - e
                && p.z <= 5.0 + r + e
                && rho >= 6.0 - e
                && rho <= 6.0 + r + e
                && (rho - 6.0 - r).hypot(p.z - 5.0 - r) >= r - 1e-4;
            plate || boss || fill
        };
        for face in explore_unique(&m, &rounded, ShapeType::Face).unwrap() {
            let mesh = ogeom_mesh::triangulate_face(
                &m,
                &face,
                ogeom_mesh::Deflection {
                    chord: 1e-3,
                    ..ogeom_mesh::Deflection::default()
                },
                T,
            )
            .unwrap();
            for p in &mesh.positions {
                assert!(
                    inside(*p),
                    "flip {flip}: outside the solid and the fill: {p:?}"
                );
            }
        }
        let key = (got.to_bits(), points(&m, &rounded));
        match &first {
            None => first = Some(key),
            Some(held) => assert_eq!(held, &key, "flip {flip}: another solid"),
        }
    }
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
