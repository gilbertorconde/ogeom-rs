//! Sectioning a face or an open shell by a plane: the curves where its
//! faces cross the plane, trimmed to the faces, with no solid behind them.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::{Curve3d as _, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Plane, Point};
use ogeom::topo::{EdgeRepr, Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn edges(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Edge).unwrap()
}

/// An edge's two ends.
fn ends(model: &Model, edge: &Shape) -> (Point, Point) {
    let data = model.node(edge).unwrap().data().as_edge().unwrap().clone();
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        panic!("an edge with no curve")
    };
    let c = model.geometry().curve(*curve).unwrap();
    (
        c.point_at(range.0, T).unwrap(),
        c.point_at(range.1, T).unwrap(),
    )
}

/// The cube's face whose points all have the given x.
fn face_at_x(model: &Model, solid: &Shape, x: f64) -> Shape {
    explore_unique(model, solid, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let b = ogeom::algo::shape_bounds(model, f, T).unwrap();
            (b.low().unwrap().x - x).abs() < 1e-6 && (b.high().unwrap().x - x).abs() < 1e-6
        })
        .unwrap()
}

fn half_space_below(model: &mut Model, height: f64) -> Shape {
    let plane = Plane::through(Point::new(0.0, 0.0, height), Direction::Z);
    let surface = PlaneSurface::over(plane, (-1e6, 1e6), (-1e6, 1e6)).unwrap();
    let face = ogeom::algo::make_natural_face(model, SurfaceGeometry::Plane(surface))
        .unwrap()
        .shape;
    ogeom::algo::make_half_space(model, &face, Point::new(0.0, 0.0, height - 1.0), T)
        .unwrap()
        .shape
}

#[test]
fn a_face_on_its_own_is_sectioned_by_a_half_space() {
    let mut model = Model::new();
    let cube = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let side = face_at_x(&model, &cube, 0.0);
    let half = half_space_below(&mut model, 4.0);
    let section = ogeom::boolean::section(&mut model, &side, &half, T)
        .unwrap()
        .shape;
    let found = edges(&model, &section);
    assert_eq!(found.len(), 1);
    let (a, b) = ends(&model, &found[0]);
    let (a, b) = if a.y <= b.y { (a, b) } else { (b, a) };
    assert!(a.distance(Point::new(0.0, 0.0, 4.0)) < 1e-9, "{a:?}");
    assert!(b.distance(Point::new(0.0, 10.0, 4.0)) < 1e-9, "{b:?}");
}

#[test]
fn a_plane_that_misses_a_face_sections_nothing_and_one_holding_it_its_boundary() {
    let mut model = Model::new();
    let cube = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let side = face_at_x(&model, &cube, 0.0);
    let clear = Plane::through(Point::new(0.0, 0.0, 12.0), Direction::Z);
    let none = ogeom::boolean::section_face(&mut model, &side, &clear, T)
        .unwrap()
        .shape;
    assert!(edges(&model, &none).is_empty());
    let holding = Plane::through(Point::ORIGIN, Direction::X);
    let boundary = ogeom::boolean::section_face(&mut model, &side, &holding, T)
        .unwrap()
        .shape;
    assert_eq!(edges(&model, &boundary).len(), 4);
}

#[test]
fn a_drum_s_side_sections_to_its_circle() {
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 3.0, 10.0, T)
        .unwrap()
        .shape;
    let wall = explore_unique(&model, &drum, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            matches!(
                model.geometry().surface(data.surface),
                Some(SurfaceGeometry::Cylinder(_))
            )
        })
        .unwrap();
    let plane = Plane::through(Point::new(0.0, 0.0, 6.0), Direction::Z);
    let section = ogeom::boolean::section_face(&mut model, &wall, &plane, T)
        .unwrap()
        .shape;
    let found = edges(&model, &section);
    assert!(!found.is_empty());
    let length: f64 = found
        .iter()
        .map(|e| {
            let data = model.node(e).unwrap().data().as_edge().unwrap().clone();
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                panic!()
            };
            let c = model.geometry().curve(*curve).unwrap();
            (0..400)
                .map(|k| {
                    let t = |k: i32| range.0 + (range.1 - range.0) * f64::from(k) / 400.0;
                    c.point_at(t(k), T)
                        .unwrap()
                        .distance(c.point_at(t(k + 1), T).unwrap())
                })
                .sum::<f64>()
        })
        .sum();
    let want = core::f64::consts::TAU * 3.0;
    assert!(
        (length - want).abs() < want * 1e-4,
        "{length} against {want}"
    );
}

#[test]
fn an_open_shell_s_sections_meet_at_the_edge_its_faces_share() {
    let mut model = Model::new();
    let cube = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let two = [
        face_at_x(&model, &cube, 0.0),
        face_at_x(&model, &cube, 10.0),
    ];
    let front = explore_unique(&model, &cube, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let b = ogeom::algo::shape_bounds(&model, f, T).unwrap();
            b.high().unwrap().y.abs() < 1e-6
        })
        .unwrap();
    let shell = ogeom::algo::make_shell(&mut model, &[two[0].clone(), front, two[1].clone()])
        .unwrap()
        .shape;
    let plane = Plane::through(Point::new(0.0, 0.0, 4.0), Direction::Z);
    let section = ogeom::boolean::section_face(&mut model, &shell, &plane, T)
        .unwrap()
        .shape;
    assert_eq!(edges(&model, &section).len(), 3);
    // Three edges end to end: four distinct vertices.
    assert_eq!(
        explore_unique(&model, &section, ShapeType::Vertex)
            .unwrap()
            .len(),
        4
    );
}
