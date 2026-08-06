//! The screen reaches into the model: a ray down a bore, a marquee over a
//! pocket, an aperture on a corner.
//!
//! Every claim here is settled against the model's own geometry — the struck
//! face *is* the drill's cylinder, at the drill's radius — rather than
//! against some other module's opinion of what the part contains. A pick that
//! agreed with a feature recognizer and disagreed with the surface it hit
//! would be wrong in the way that matters.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::SurfaceGeometry;
use ogeom_math::{Direction, Frame, Point, Point2, Vector};
use ogeom_mesh::Deflection;
use ogeom_select::{Marquee, PickKind, Pickable, Ray};
use ogeom_topo::{Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A block with a through bore and a milled pocket.
fn worked_block(model: &mut Model) -> Shape {
    let block = ogeom_algo::make_box(model, Frame::WORLD, (40.0, 30.0, 12.0), T)
        .unwrap()
        .shape;
    let drill_frame =
        Frame::new(Point::new(30.0, 15.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom_algo::make_cylinder(model, drill_frame, 4.0, 14.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom_bool::cut(model, &block, &drill, T).unwrap().shape;
    let mill_frame = Frame::new(Point::new(5.0, 5.0, 7.0), Direction::Z, Direction::X, T).unwrap();
    let mill = ogeom_algo::make_box(model, mill_frame, (12.0, 12.0, 6.0), T)
        .unwrap()
        .shape;
    ogeom_bool::cut(model, &drilled, &mill, T).unwrap().shape
}

/// The planar face of `shape` whose plane passes through `on`.
fn planar_face_at(model: &Model, shape: &Shape, on: Point) -> Shape {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let NodeData::Face(data) = model.node(f).unwrap().data() else {
                return false;
            };
            let Some(SurfaceGeometry::Plane(p)) = model.geometry().surface(data.surface) else {
                return false;
            };
            p.plane().distance_to(on).abs() < 1e-9
        })
        .expect("a planar face there")
}

/// A ray sent down the bore strikes the bore wall, and the wall is the
/// cylinder the drill cut, at the radius it was cut with.
#[test]
fn a_ray_down_the_bore_strikes_the_drills_own_cylinder() {
    let mut model = Model::new();
    let part = worked_block(&mut model);
    let scene = Pickable::build(&model, &part, Deflection::default(), T).unwrap();

    // Into the bore at a slant: straight down the void would graze along the
    // wall forever and strike nothing.
    let hit = scene
        .pick_first(
            Ray {
                origin: Point::new(30.0, 15.0, 20.0),
                direction: Vector::new(0.0, 0.3, -1.0),
            },
            0.0,
        )
        .expect("the ray strikes the bore");
    let struck = scene.triangle_face(hit.triangle).unwrap();

    let NodeData::Face(data) = model.node(struck).unwrap().data() else {
        panic!("a pick names a face");
    };
    let Some(SurfaceGeometry::Cylinder(c)) = model.geometry().surface(data.surface) else {
        panic!("the bore wall is a cylinder, not {:?}", data.surface);
    };
    assert!(
        (c.cylinder().radius() - 4.0).abs() < 1e-9,
        "radius {}",
        c.cylinder().radius()
    );
    // And the axis it was drilled about.
    let axis = c.cylinder().frame().origin();
    assert!((axis.x - 30.0).abs() < 1e-9 && (axis.y - 15.0).abs() < 1e-9);
}

/// A marquee over the pocket corner reaches its floor — the face whose plane
/// is the depth the pocket was milled to.
#[test]
fn a_marquee_over_the_pocket_reaches_its_floor() {
    let mut model = Model::new();
    let part = worked_block(&mut model);
    let scene = Pickable::build(&model, &part, Deflection::default(), T).unwrap();
    let floor = planar_face_at(&model, &part, Point::new(11.0, 11.0, 7.0));

    let picked = scene.select_rectangle(
        &Frame::WORLD,
        Point2::new(4.0, 4.0),
        Point2::new(18.5, 18.5),
        Marquee::Crossing,
    );
    assert!(
        picked.iter().any(|s| s.is_same(&floor)),
        "the marquee reaches the pocket floor"
    );
}

/// Sub-shape granularity: an aperture near a pocket corner resolves to a
/// vertex rather than to the face the corner sits on.
#[test]
fn an_aperture_on_a_corner_resolves_to_a_vertex() {
    let mut model = Model::new();
    let part = worked_block(&mut model);
    let scene = Pickable::build(&model, &part, Deflection::default(), T).unwrap();

    let corner = scene
        .pick_first(
            Ray {
                origin: Point::new(5.01, 5.01, 20.0),
                direction: Vector::new(0.0, 0.0, -1.0),
            },
            0.2,
        )
        .expect("the corner is pickable");
    assert_eq!(corner.kind, PickKind::Vertex);
}
