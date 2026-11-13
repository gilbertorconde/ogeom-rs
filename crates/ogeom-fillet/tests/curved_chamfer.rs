//! Chamfers of solids with curved faces through the blend engine: what the
//! history tells of each element, beside the M3 corpus's volumes.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn up_at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(
        Point::new(x, y, z),
        Direction::new(Vector::new(0.0, 0.0, 1.0), T).unwrap(),
        Direction::new(Vector::new(1.0, 0.0, 0.0), T).unwrap(),
        T,
    )
    .unwrap()
}

fn point_of(m: &Model, v: &Shape) -> Point {
    m.node(v)
        .and_then(|n| n.data().as_vertex().map(|d| d.point))
        .unwrap()
}

fn surface_kind(m: &Model, face: &Shape) -> ogeom_geom::SurfaceKind {
    use ogeom_geom::Surface as _;
    let data = m.node(face).unwrap().data().as_face().unwrap();
    m.geometry().surface(data.surface).unwrap().kind()
}

/// Every edge chamfered is deleted and generates one face on `bevel`, every
/// face of the solid is modified into one face, and the result is valid.
fn history_holds(
    m: &Model,
    solid: &Shape,
    built: &ogeom_algo::Built,
    edges: &[Shape],
    bevel: ogeom_geom::SurfaceKind,
) {
    assert!(ogeom_algo::check(m, &built.shape, T).unwrap().is_valid());
    for e in edges {
        assert!(built.history.is_deleted(e), "a chamfered edge left");
        let made = built.history.generated(e);
        assert_eq!(made.len(), 1, "a chamfered edge without its one bevel");
        assert_eq!(
            surface_kind(m, &made[0]),
            bevel,
            "a bevel on another surface"
        );
    }
    for f in explore_unique(m, solid, ShapeType::Face).unwrap() {
        assert_eq!(
            built.history.modified(&f).len(),
            1,
            "a face not carried to one face"
        );
    }
}

/// A tube of radius 3 bored to 2, 4 tall, its four rims chamfered 0.3:
/// each rim is deleted into one cone, each face carried to one face.
#[test]
fn a_tube_s_rims_chamfer_into_cones_their_faces_carried() {
    let mut m = Model::new();
    let outer = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 4.0, T)
        .unwrap()
        .shape;
    let inner = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, -1.0), 2.0, 6.0, T)
        .unwrap()
        .shape;
    let tube = ogeom_bool::cut(&mut m, &outer, &inner, T).unwrap().shape;
    let rims: Vec<Shape> = explore_unique(&m, &tube, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            matches!(
                ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap().0,
                ogeom_geom::Curve::Circle(_)
            )
        })
        .collect();
    assert!(rims.len() >= 4);
    let built = ogeom_fillet::chamfer_edges(&mut m, &tube, &rims, 0.3, T).unwrap();
    history_holds(&m, &tube, &built, &rims, ogeom_geom::SurfaceKind::Cone);
}

/// The D-shaped block (a drum of radius 3 cut in half along its axis, 10
/// long), its two straight edges chamfered 0.5: each edge deleted into one
/// plane, each face carried, and each cap's corners where the bevel ends
/// gone, the cap running round the bevel's end there.
#[test]
fn a_d_block_s_straight_edges_chamfer_into_planes_ending_on_the_caps() {
    let mut m = Model::new();
    let full = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 10.0, T)
        .unwrap()
        .shape;
    let half = ogeom_algo::make_box(&mut m, up_at(-4.0, -4.0, -1.0), (4.0, 8.0, 12.0), T)
        .unwrap()
        .shape;
    let d = ogeom_bool::cut(&mut m, &full, &half, T).unwrap().shape;
    let straight: Vec<Shape> = explore_unique(&m, &d, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
            let (p, q) = (point_of(&m, &a), point_of(&m, &b));
            p.x.abs() < 1e-9 && q.x.abs() < 1e-9 && (p.z - q.z).abs() > 1.0
        })
        .collect();
    assert_eq!(straight.len(), 2);
    let built = ogeom_fillet::chamfer_edges(&mut m, &d, &straight, 0.5, T).unwrap();
    history_holds(&m, &d, &built, &straight, ogeom_geom::SurfaceKind::Plane);
    for e in &straight {
        let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
        for v in [a, b] {
            assert!(
                built.history.is_deleted(&v),
                "a corner the bevel takes left"
            );
        }
    }
}
