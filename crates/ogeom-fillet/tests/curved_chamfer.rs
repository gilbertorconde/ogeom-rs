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

/// The D-shaped block's straight edges chamfered 0.5 along the flat at 45
/// degrees from it: in the section the bevel leaves `(0, 2.5)` towards
/// the edge and into the material, meeting the wall of radius 3 where
/// `2w^2 + 5w - 2.75 = 0` for `w` its run along each axis. The half disc
/// less, twice along the 10, the triangle from the corner to the two
/// contacts and the wall's segment beyond the chord.
#[test]
fn a_d_block_s_edges_chamfer_by_distance_and_angle_on_the_flat() {
    let mut m = Model::new();
    let full = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 10.0, T)
        .unwrap()
        .shape;
    let half = ogeom_algo::make_box(&mut m, up_at(-4.0, -4.0, -1.0), (4.0, 8.0, 12.0), T)
        .unwrap()
        .shape;
    let d = ogeom_bool::cut(&mut m, &full, &half, T).unwrap().shape;
    let flat = explore_unique(&m, &d, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let (_, n) = ogeom_algo::face_normal(&m, f, T).unwrap();
            n.x < -0.5
        })
        .unwrap();
    let specs: Vec<(Shape, ogeom_fillet::Chamfer)> = explore_unique(&m, &d, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
            let (p, q) = (point_of(&m, &a), point_of(&m, &b));
            p.x.abs() < 1e-9 && q.x.abs() < 1e-9 && (p.z - q.z).abs() > 1.0
        })
        .map(|e| {
            (
                e,
                ogeom_fillet::Chamfer::Angle {
                    face: flat.clone(),
                    distance: 0.5,
                    angle: core::f64::consts::FRAC_PI_4,
                },
            )
        })
        .collect();
    assert_eq!(specs.len(), 2);
    let built = ogeom_fillet::chamfer_edges_with(&mut m, &d, &specs, T).unwrap();
    assert!(ogeom_algo::check(&m, &built.shape, T).unwrap().is_valid());
    for (e, _) in &specs {
        assert_eq!(
            built.history.generated(e).len(),
            1,
            "an edge without its bevel"
        );
    }
    let w = (-5.0 + (25.0_f64 + 22.0).sqrt()) / 4.0;
    let (corner, on_flat, on_wall) = ((0.0, 3.0), (0.0, 2.5), (w, 2.5 + w));
    let triangle = ((on_flat.0 - corner.0) * (on_wall.1 - corner.1)
        - (on_flat.1 - corner.1) * (on_wall.0 - corner.0))
        .abs()
        / 2.0;
    let phi = (w / 3.0).asin();
    let segment = 9.0 / 2.0 * (phi - phi.sin());
    let want = (core::f64::consts::PI * 9.0 / 2.0 - 2.0 * (triangle + segment)) * 10.0;
    let got = volume(&m, &built.shape);
    assert!((got - want).abs() < 1e-6, "volume {got} against {want}");
}

/// A drum's top rim chamfered 0.3 along the cap at 45 degrees from it:
/// on its square rim the same bevel as 0.3 along both faces, the same
/// solid's volume.
#[test]
fn a_rim_chamfers_by_distance_and_angle_as_its_square_corner_says() {
    let make = |m: &mut Model| {
        ogeom_algo::make_cylinder(m, up_at(0.0, 0.0, 0.0), 3.0, 4.0, T)
            .unwrap()
            .shape
    };
    let rim_and_cap = |m: &Model, s: &Shape| {
        let rim = explore_unique(m, s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let (a, _) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
                matches!(
                    ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap().0,
                    ogeom_geom::Curve::Circle(_)
                ) && point_of(m, &a).z > 3.0
            })
            .unwrap();
        let cap = explore_unique(m, s, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| ogeom_algo::face_normal(m, f, T).unwrap().1.z > 0.5)
            .unwrap();
        (rim, cap)
    };
    let mut m = Model::new();
    let drum = make(&mut m);
    let (rim, cap) = rim_and_cap(&m, &drum);
    let angled = ogeom_fillet::chamfer_edges_with(
        &mut m,
        &drum,
        &[(
            rim.clone(),
            ogeom_fillet::Chamfer::Angle {
                face: cap,
                distance: 0.3,
                angle: core::f64::consts::FRAC_PI_4,
            },
        )],
        T,
    )
    .unwrap();
    let square = ogeom_fillet::chamfer_edges(&mut m, &drum, &[rim], 0.3, T).unwrap();
    assert!(ogeom_algo::check(&m, &angled.shape, T).unwrap().is_valid());
    let (a, b) = (volume(&m, &angled.shape), volume(&m, &square.shape));
    assert!((a - b).abs() < 1e-9 * b, "{a} against {b}");
}
