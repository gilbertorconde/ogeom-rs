//! Setback corners: the bands at a vertex stopped a chosen distance along
//! their edges, and one patch over the opening between them.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{Setbacks, analyse_blend, fillet_edges, round_vertex_setback};
use ogeom_geom::{Continuity, Curve3d as _};
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

/// The edge of `s` passing nearest `p`, 64 samples along each.
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

/// The vertex of `s` standing at `p`.
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

fn volume(m: &Model, s: &Shape) -> f64 {
    let deflection = ogeom_mesh::Deflection {
        chord: 1e-3,
        ..ogeom_mesh::Deflection::default()
    };
    ogeom_algo::volume_properties(m, s, deflection, T)
        .unwrap()
        .mass
}

fn valid(m: &Model, s: &Shape) {
    let diagnosis = ogeom_algo::check(m, s, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
}

/// A 20 x 20 x 10 block and its corner at the origin's top, (0, 0, 10),
/// with the three edges leaving it.
fn block() -> (Model, Shape, Shape, Vec<Shape>) {
    let mut m = Model::new();
    let b = ogeom_algo::make_box(&mut m, at(0.0, 0.0, 0.0), (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edges: Vec<Shape> = [[10.0, 0.0, 10.0], [0.0, 10.0, 10.0], [0.0, 0.0, 5.0]]
        .iter()
        .map(|p| edge_near(&m, &b, *p))
        .collect();
    let v = vertex_at(&m, &b, [0.0, 0.0, 10.0]);
    (m, b, v, edges)
}

fn even(setback: f64, continuity: Continuity) -> Setbacks {
    Setbacks {
        per_edge: Vec::new(),
        default: setback,
        continuity,
    }
}

#[test]
fn a_setback_of_the_radius_is_the_ball_s_corner() {
    let (mut m, b, v, edges) = block();
    let setback = round_vertex_setback(&mut m, &b, &v, 2.0, &even(2.0, Continuity::G1), T)
        .unwrap()
        .shape;
    let ball = fillet_edges(&mut m, &b, &edges, 2.0, T).unwrap().shape;
    valid(&m, &setback);
    let (a, c) = (volume(&m, &setback), volume(&m, &ball));
    assert!((a - c).abs() < 1e-12, "{a} against {c}");
}

/// The patch a setback corner of `continuity` builds at the block's
/// corner, its bands stopped 5 from the vertex, and the patch's worst
/// tangency against the faces round it.
fn five_back(continuity: Continuity) -> (Model, Shape, Shape, f64) {
    let (mut m, b, v, _) = block();
    let rounded = round_vertex_setback(&mut m, &b, &v, 2.0, &even(5.0, continuity), T).unwrap();
    valid(&m, &rounded.shape);
    let patch = rounded.history.generated(&v)[0].clone();
    let contacts = analyse_blend(&m, &rounded.shape, &patch, 17, T).unwrap();
    assert_eq!(
        contacts.len(),
        6,
        "three band ends and three faces: {contacts:?}"
    );
    let tangency = contacts
        .iter()
        .map(|c| c.tangency_error)
        .fold(0.0, f64::max);
    (m, rounded.shape, patch, tangency)
}

#[test]
fn a_setback_corner_stops_each_band_at_its_setback_tangent_to_the_patch() {
    let (m, rounded, patch, tangency) = five_back(Continuity::G1);
    assert!(
        tangency < 1e-3,
        "the patch meets its neighbours within {tangency} rad"
    );
    // Every point of the patch's boundary on a band's end lies 5 along
    // that band's edge from the vertex: on x = 5, y = 5 or z = 10 - 5.
    for e in explore_unique(&m, &patch, ShapeType::Edge).unwrap() {
        let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(&m, &e, T).unwrap().unwrap();
        let p = c.point_at(f64::midpoint(lo, hi), T).unwrap();
        let on_face = p.x.abs() < 1e-6 || p.y.abs() < 1e-6 || (p.z - 10.0).abs() < 1e-6;
        let on_end =
            (p.x - 5.0).abs() < 1e-6 || (p.y - 5.0).abs() < 1e-6 || (p.z - 5.0).abs() < 1e-6;
        assert!(on_face || on_end, "a patch edge at {p:?}");
    }
    // The setback corner takes off more than the ball's corner does.
    let (mut m2, b, v, edges) = block();
    let ball = fillet_edges(&mut m2, &b, &edges, 2.0, T).unwrap().shape;
    let _ = v;
    assert!(volume(&m, &rounded) < volume(&m2, &ball));
}

#[test]
fn a_setback_corner_asked_g2_is_refused_by_name_while_bands_are_circular() {
    // At a band's arc end on a face the arc bends the patch at 1/r along
    // itself; the face asks no bend across its trim, which runs along
    // the arc's direction there.
    let (mut m, b, v, _) = block();
    let refused = round_vertex_setback(&mut m, &b, &v, 2.0, &even(5.0, Continuity::G2), T);
    let Err(e) = refused else {
        panic!("a G2 setback corner built against circular bands");
    };
    assert!(
        e.to_string()
            .contains("curvature continuous with their faces"),
        "{e}"
    );
}

#[test]
fn a_setback_corner_rounds_a_square_pyramid_s_apex() {
    let mut m = Model::new();
    let base = ogeom_algo::make_box(&mut m, at(-10.0, -10.0, -10.0), (20.0, 20.0, 20.0), T)
        .unwrap()
        .shape;
    let mut solid = base;
    // Four planes leaning in from the base's top edges to the apex at
    // (0, 0, 10), each cutting a half-space off.
    for (nx, ny) in [(1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)] {
        let normal = ogeom_math::Vector::new(nx, ny, 1.0);
        let across = ogeom_math::Vector::new(-ny, nx, 0.0);
        let side = normal.cross(across) / normal.cross(across).magnitude();
        let corner = Point::new(0.0, 0.0, 10.0) - across * 40.0 - side * 40.0;
        let frame = Frame::new(
            corner,
            Direction::new(normal, T).unwrap(),
            Direction::new(across, T).unwrap(),
            T,
        )
        .unwrap();
        let away = ogeom_algo::make_box(&mut m, frame, (80.0, 80.0, 40.0), T)
            .unwrap()
            .shape;
        solid = ogeom_bool::cut(&mut m, &solid, &away, T).unwrap().shape;
    }
    let apex = vertex_at(&m, &solid, [0.0, 0.0, 10.0]);
    let rounded =
        round_vertex_setback(&mut m, &solid, &apex, 1.0, &even(4.0, Continuity::G1), T).unwrap();
    valid(&m, &rounded.shape);
}

#[test]
fn a_setback_shorter_than_the_ball_s_is_refused() {
    let (mut m, b, v, _) = block();
    let refused = round_vertex_setback(&mut m, &b, &v, 2.0, &even(1.5, Continuity::G1), T);
    let Err(e) = refused else {
        panic!("a setback inside the ball's corner rounded");
    };
    assert!(
        e.to_string().contains("inside the ball's own corner"),
        "{e}"
    );
}

fn probe_pyramid(m: &mut Model) -> Shape {
    let base = ogeom_algo::make_box(m, at(-10.0, -10.0, -10.0), (20.0, 20.0, 20.0), T)
        .unwrap()
        .shape;
    let mut solid = base;
    for (nx, ny) in [(1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)] {
        let normal = ogeom_math::Vector::new(nx, ny, 1.0);
        let across = ogeom_math::Vector::new(-ny, nx, 0.0);
        let side = normal.cross(across) / normal.cross(across).magnitude();
        let corner = Point::new(0.0, 0.0, 10.0) - across * 40.0 - side * 40.0;
        let frame = Frame::new(
            corner,
            Direction::new(normal, T).unwrap(),
            Direction::new(across, T).unwrap(),
            T,
        )
        .unwrap();
        let away = ogeom_algo::make_box(m, frame, (80.0, 80.0, 40.0), T)
            .unwrap()
            .shape;
        solid = ogeom_bool::cut(m, &solid, &away, T).unwrap().shape;
    }
    solid
}

#[test]
fn probe_variants() {
    for (r, s, uneven) in [
        (1.0, 4.0, None),
        (1.0, 6.0, None),
        (1.5, 4.0, None),
        (1.0, 3.0, None),
        (1.0, 4.0, Some(6.0)),
        (1.0, 3.0, Some(5.0)),
        (0.5, 4.0, None),
        (1.0, 8.0, None),
    ] {
        let mut m = Model::new();
        let solid = probe_pyramid(&mut m);
        let apex = vertex_at(&m, &solid, [0.0, 0.0, 10.0]);
        let mut sb = even(s, Continuity::G1);
        if let Some(other) = uneven {
            let edges: Vec<Shape> = explore_unique(&m, &solid, ShapeType::Edge)
                .unwrap()
                .into_iter()
                .filter(|e| {
                    let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
                    a.is_same(&apex) || b.is_same(&apex)
                })
                .collect();
            sb.per_edge = vec![(edges[0].clone(), other), (edges[2].clone(), other)];
        }
        match round_vertex_setback(&mut m, &solid, &apex, r, &sb, T) {
            Ok(b) => eprintln!(
                "PROBE r {r} s {s} uneven {uneven:?}: ok faces {} valid {}",
                b.history.generated(&apex).len(),
                ogeom_algo::check(&m, &b.shape, T).unwrap().is_valid()
            ),
            Err(e) => eprintln!("PROBE r {r} s {s} uneven {uneven:?}: {e}"),
        }
    }
}

#[test]
fn probe_mixed() {
    let quads = |m: &Model, s: &Shape| {
        explore_unique(m, s, ShapeType::Face)
            .unwrap()
            .into_iter()
            .filter(|f| {
                m.node(f)
                    .and_then(|n| n.data().as_face())
                    .and_then(|d| m.geometry().surface(d.surface))
                    .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::BSpline(_)))
            })
            .count()
    };
    for radii in [
        [1.0, 3.0, 2.0],
        [1.0, 2.5, 2.0],
        [1.0, 2.2, 2.0],
        [1.5, 3.0, 2.0],
        [1.0, 3.0, 3.0],
        [1.0, 1.0, 3.0],
    ] {
        let (mut m, b, _, edges) = block();
        let asked: Vec<(Shape, f64)> = edges.iter().cloned().zip(radii).collect();
        match ogeom_fillet::fillet_edges_with(&mut m, &b, &asked, T) {
            Ok(r) => eprintln!(
                "PROBE box {radii:?}: ok bspline faces {} valid {}",
                quads(&m, &r.shape),
                ogeom_algo::check(&m, &r.shape, T).unwrap().is_valid()
            ),
            Err(e) => eprintln!(
                "PROBE box {radii:?}: {}",
                e.to_string().chars().take(120).collect::<String>()
            ),
        }
    }
    for radii in [
        [1.0, 1.5, 1.0, 1.5],
        [1.0, 2.0, 1.0, 2.0],
        [1.0, 1.2, 1.4, 1.6],
        [0.5, 1.5, 0.5, 1.5],
    ] {
        let mut m = Model::new();
        let solid = probe_pyramid(&mut m);
        let apex = vertex_at(&m, &solid, [0.0, 0.0, 10.0]);
        let edges: Vec<Shape> = explore_unique(&m, &solid, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
                a.is_same(&apex) || b.is_same(&apex)
            })
            .collect();
        let asked: Vec<(Shape, f64)> = edges.iter().cloned().zip(radii).collect();
        match ogeom_fillet::fillet_edges_with(&mut m, &solid, &asked, T) {
            Ok(r) => eprintln!(
                "PROBE pyr {radii:?}: ok bspline faces {} valid {}",
                quads(&m, &r.shape),
                ogeom_algo::check(&m, &r.shape, T).unwrap().is_valid()
            ),
            Err(e) => eprintln!(
                "PROBE pyr {radii:?}: {}",
                e.to_string().chars().take(120).collect::<String>()
            ),
        }
    }
}
