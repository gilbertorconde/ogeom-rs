//! Where a curve first and last meets a shape, against closed forms: a
//! line and a helix on a placed box and on a cylinder, and a line that
//! crosses a face's surface only outside its trim.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_algo::{CurveHit, Hit, curve_hit, make_box, make_cylinder};
use ogeom_core::Tolerances;
use ogeom_geom::{Curve, HelixCurve, LineCurve};
use ogeom_math::{Axis, Direction, Frame, Point};
use ogeom_topo::{Model, NodeData, Shape};

const T: Tolerances = Tolerances::millimetres();

/// The bound a closed-form hit is held to.
const EXACT: f64 = 1e-9;

fn line(from: Point, direction: Direction, length: f64) -> Curve {
    LineCurve::over(
        Axis {
            location: from,
            direction,
        },
        0.0,
        length,
    )
    .unwrap()
    .into()
}

fn hit(model: &Model, shape: &Shape, curve: &Curve, range: (f64, f64), which: Hit) -> CurveHit {
    curve_hit(model, shape, curve, range, which, T)
        .unwrap()
        .expect("a hit")
}

/// The constant coordinate `axis` (0, 1, 2) of a planar face's vertices,
/// or `None` where they do not share one.
fn level(model: &Model, face: &Shape, axis: usize) -> Option<f64> {
    let vertices = ogeom_topo::explore_unique(model, face, ogeom_topo::ShapeType::Vertex).unwrap();
    let values: Vec<f64> = vertices
        .iter()
        .map(|v| {
            let Some(NodeData::Vertex(data)) = model.node(v).map(|n| n.data()) else {
                unreachable!()
            };
            let p = v.transform(model.datums()).unwrap().apply(data.point);
            [p.x, p.y, p.z][axis]
        })
        .collect();
    let first = *values.first()?;
    values
        .iter()
        .all(|x| (x - first).abs() < 1e-9)
        .then_some(first)
}

/// A line along X through a box placed at (3, 1, 0), 10 x 6 x 4: it enters
/// at x = 3 and leaves at x = 13, and names the faces there.
#[test]
fn a_line_through_a_placed_box_enters_and_leaves_by_its_end_faces() {
    let mut model = Model::new();
    let frame = Frame::new(Point::new(3.0, 1.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let block = make_box(&mut model, frame, (10.0, 6.0, 4.0), T)
        .unwrap()
        .shape;
    let ray = line(Point::new(-5.0, 2.5, 1.5), Direction::X, 40.0);
    let first = hit(&model, &block, &ray, (0.0, 40.0), Hit::First);
    assert!((first.parameter - 8.0).abs() < EXACT, "{}", first.parameter);
    assert!(
        first.point.distance(Point::new(3.0, 2.5, 1.5)) < EXACT,
        "{:?}",
        first.point
    );
    assert_eq!(level(&model, &first.face, 0), Some(3.0));
    let last = hit(&model, &block, &ray, (0.0, 40.0), Hit::Last);
    assert!((last.parameter - 18.0).abs() < EXACT, "{}", last.parameter);
    assert!(
        last.point.distance(Point::new(13.0, 2.5, 1.5)) < EXACT,
        "{:?}",
        last.point
    );
    assert_eq!(level(&model, &last.face, 0), Some(13.0));
    // A range that starts inside the box meets only the far face.
    let inside = hit(&model, &block, &ray, (10.0, 40.0), Hit::First);
    assert!((inside.parameter - 18.0).abs() < EXACT);
}

/// A line that crosses the planes of the box's faces only outside their
/// trims meets nothing, though `intersect_curve_surface` finds piercings
/// of those planes.
#[test]
fn a_line_crossing_only_a_faces_carrier_meets_nothing() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (10.0, 6.0, 4.0), T)
        .unwrap()
        .shape;
    let ray = line(Point::new(-5.0, 20.0, 1.0), Direction::X, 40.0);
    assert!(
        curve_hit(&model, &block, &ray, (0.0, 40.0), Hit::First, T)
            .unwrap()
            .is_none()
    );
    // Through the box's top plane beside the box.
    let diagonal = Direction::new(ogeom_math::Vector::new(1.0, 0.0, 1.0), T).unwrap();
    let ray = line(Point::new(-8.0, 3.0, 0.0), diagonal, 20.0);
    assert!(
        curve_hit(&model, &block, &ray, (0.0, 20.0), Hit::First, T)
            .unwrap()
            .is_none()
    );
}

/// A line along X at y = 2, z = 1 through the cylinder of radius 4 about Z:
/// it pierces the side at x = -sqrt(12) and x = sqrt(12).
#[test]
fn a_line_through_a_cylinder_pierces_its_side_where_the_quadratic_says() {
    let mut model = Model::new();
    let rod = make_cylinder(&mut model, Frame::WORLD, 4.0, 10.0, T)
        .unwrap()
        .shape;
    let ray = line(Point::new(-20.0, 2.0, 1.0), Direction::X, 40.0);
    let root = 12.0_f64.sqrt();
    let first = hit(&model, &rod, &ray, (0.0, 40.0), Hit::First);
    assert!(
        (first.parameter - (20.0 - root)).abs() < EXACT,
        "{}",
        first.parameter
    );
    assert!(first.point.distance(Point::new(-root, 2.0, 1.0)) < EXACT);
    let last = hit(&model, &rod, &ray, (0.0, 40.0), Hit::Last);
    assert!(
        (last.parameter - (20.0 + root)).abs() < EXACT,
        "{}",
        last.parameter
    );
    assert!(last.point.distance(Point::new(root, 2.0, 1.0)) < EXACT);
}

/// A helix of radius 5 and pitch 2 about Z starts at (5, 0, 0) beside a
/// box filling x <= 3: it first reaches x = 3 at the angle acos(3 / 5),
/// at (3, 4, acos(0.6) / pi), and last leaves it on its third turn,
/// before its end, at the angle 6 pi - acos(0.6), at (3, -4, ...).
#[test]
fn a_helix_meets_a_box_where_its_angle_reaches_the_wall() {
    let mut model = Model::new();
    let frame = Frame::new(
        Point::new(-20.0, -10.0, -1.0),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    let block = make_box(&mut model, frame, (23.0, 20.0, 20.0), T)
        .unwrap()
        .shape;
    let helix: Curve = HelixCurve::new(Frame::WORLD, 5.0, 2.0, 3.0).unwrap().into();
    let range = (0.0, 6.0 * core::f64::consts::PI);
    let angle = 0.6_f64.acos();
    let first = hit(&model, &block, &helix, range, Hit::First);
    assert!(
        (first.parameter - angle).abs() < EXACT,
        "{}",
        first.parameter
    );
    let want = Point::new(3.0, 4.0, angle / core::f64::consts::PI);
    assert!(first.point.distance(want) < EXACT, "{:?}", first.point);
    assert_eq!(level(&model, &first.face, 0), Some(3.0));
    let last = hit(&model, &block, &helix, range, Hit::Last);
    let back = range.1 - angle;
    assert!((last.parameter - back).abs() < EXACT, "{}", last.parameter);
    let want = Point::new(3.0, -4.0, back / core::f64::consts::PI);
    assert!(last.point.distance(want) < EXACT, "{:?}", last.point);
}

/// A helix of radius 5 and pitch 2 rising from z = -2 through a cylinder of
/// radius 8 standing on z = 0, 3 tall: it enters the bottom cap after one
/// turn, at (5, 0, 0), and leaves by the top cap at the angle 5 pi, at
/// (-5, 0, 3).
#[test]
fn a_helix_meets_a_cylinder_by_its_caps() {
    let mut model = Model::new();
    let rod = make_cylinder(&mut model, Frame::WORLD, 8.0, 3.0, T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(0.0, 0.0, -2.0), Direction::Z, Direction::X, T).unwrap();
    let helix: Curve = HelixCurve::new(frame, 5.0, 2.0, 4.0).unwrap().into();
    let range = (0.0, 8.0 * core::f64::consts::PI);
    let pi = core::f64::consts::PI;
    let first = hit(&model, &rod, &helix, range, Hit::First);
    assert!(
        (first.parameter - 2.0 * pi).abs() < EXACT,
        "{}",
        first.parameter
    );
    assert!(
        first.point.distance(Point::new(5.0, 0.0, 0.0)) < EXACT,
        "{:?}",
        first.point
    );
    assert_eq!(level(&model, &first.face, 2), Some(0.0));
    let last = hit(&model, &rod, &helix, range, Hit::Last);
    assert!(
        (last.parameter - 5.0 * pi).abs() < EXACT,
        "{}",
        last.parameter
    );
    assert!(
        last.point.distance(Point::new(-5.0, 0.0, 3.0)) < EXACT,
        "{:?}",
        last.point
    );
    assert_eq!(level(&model, &last.face, 2), Some(3.0));
}

/// A range outside the curve's domain is refused, not clamped.
#[test]
fn a_range_off_the_curve_is_refused() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
        .unwrap()
        .shape;
    let ray = line(Point::new(-5.0, 0.5, 0.5), Direction::X, 10.0);
    assert!(curve_hit(&model, &block, &ray, (0.0, 20.0), Hit::First, T).is_err());
    assert!(curve_hit(&model, &block, &ray, (3.0, 3.0), Hit::First, T).is_err());
}
