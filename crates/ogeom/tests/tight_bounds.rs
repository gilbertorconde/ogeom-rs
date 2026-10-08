//! How big a body is: the smallest box holding it, not the carriers'.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    face_bounds, make_box, make_cylinder, make_face, make_polygon, make_revolution, make_sphere,
    tight_bounds,
};
use ogeom::core::Tolerances;
use ogeom::geom::PlaneSurface;
use ogeom::math::{Aabb, Axis, Direction, Frame, Plane, Point, Transform, Vector};
use ogeom::mesh::Deflection;
use ogeom::offset::move_faces;
use ogeom::topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn near(a: Point, b: Point) {
    assert!(a.distance(b) < 1e-6, "{a:?} against {b:?}");
}

#[test]
fn a_revolved_tube_reaches_its_outer_radius() {
    let mut model = Model::new();
    let pts =
        [(4.0, 0.0), (6.0, 0.0), (6.0, 10.0), (4.0, 10.0)].map(|(x, y)| Point::new(x, y, 0.0));
    let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
    let face = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Y,
    };
    let tube = make_revolution(&mut model, &face, axis, core::f64::consts::TAU, T)
        .unwrap()
        .shape;
    let b = tight_bounds(&model, &tube, T).unwrap();
    near(b.low().unwrap(), Point::new(-6.0, 0.0, -6.0));
    near(b.high().unwrap(), Point::new(6.0, 10.0, 6.0));
}

#[test]
fn a_ball_and_a_box_are_as_big_as_they_are() {
    let mut model = Model::new();
    let ball = make_sphere(&mut model, Frame::WORLD, 5.0, T).unwrap().shape;
    let b = tight_bounds(&model, &ball, T).unwrap();
    near(b.low().unwrap(), Point::new(-5.0, -5.0, -5.0));
    near(b.high().unwrap(), Point::new(5.0, 5.0, 5.0));
    let block = make_box(&mut model, Frame::WORLD, (3.0, 4.0, 5.0), T)
        .unwrap()
        .shape;
    let b = tight_bounds(&model, &block, T).unwrap();
    near(b.low().unwrap(), Point::ORIGIN);
    near(b.high().unwrap(), Point::new(3.0, 4.0, 5.0));
}

/// A square plate drilled through with a grid of holes of radius 1.5,
/// `pitch` apart and from the edges.
fn drilled_plate(model: &mut Model, rows: u32, pitch: f64, height: f64) -> ogeom::topo::Shape {
    let size = pitch * f64::from(rows) + pitch;
    let block = make_box(model, Frame::WORLD, (size, size, height), T)
        .unwrap()
        .shape;
    let mut pins = Vec::new();
    for i in 0..rows {
        for j in 0..rows {
            let at = Point::new(
                pitch + pitch * f64::from(i),
                pitch + pitch * f64::from(j),
                -height,
            );
            let frame = Frame::new(at, Direction::Z, Direction::X, T).unwrap();
            pins.push(
                make_cylinder(model, frame, 1.5, 3.0 * height, T)
                    .unwrap()
                    .shape,
            );
        }
    }
    let pins = model.add_compound(&pins).unwrap();
    ogeom::boolean::cut(model, &block, &pins, T).unwrap().shape
}

/// Every face's kept box holds a fine mesh of the face and stands off its
/// extremes by no more than the mesh's chord and the face's tolerance.
fn holds_its_faces(model: &Model, shape: &ogeom::topo::Shape) {
    let chord = 1e-3;
    for face in explore_unique(model, shape, ShapeType::Face).unwrap() {
        let kept = face_bounds(model, &face).unwrap();
        assert_eq!(model.face_bounds(&face), Some(kept));
        let tolerance = model.tolerance_of(&face).unwrap().unwrap().get();
        let mesh =
            ogeom::mesh::triangulate_face(model, &face, Deflection::with_chord(chord).unwrap(), T)
                .unwrap();
        let reached = Aabb::of_points(&mesh.positions);
        let (low, high) = (kept.low().unwrap(), kept.high().unwrap());
        let (near_low, near_high) = (reached.low().unwrap(), reached.high().unwrap());
        let slack = chord + tolerance + 1e-6;
        for axis in 0..3 {
            let pick = |p: Point| [p.x, p.y, p.z][axis];
            assert!(
                pick(low) <= pick(near_low) + 1e-9,
                "{kept} against {reached}"
            );
            assert!(
                pick(high) >= pick(near_high) - 1e-9,
                "{kept} against {reached}"
            );
            assert!(
                pick(near_low) - pick(low) <= slack,
                "{kept} against {reached}"
            );
            assert!(
                pick(high) - pick(near_high) <= slack,
                "{kept} against {reached}"
            );
        }
    }
}

#[test]
fn each_face_keeps_a_box_that_holds_it_and_no_more() {
    let mut model = Model::new();
    let plate = drilled_plate(&mut model, 2, 10.0, 5.0);
    // A boolean's faces on planes and drums come out with their boxes.
    for face in explore_unique(&model, &plate, ShapeType::Face).unwrap() {
        assert!(model.face_bounds(&face).is_some());
    }
    let b = tight_bounds(&model, &plate, T).unwrap();
    near(b.low().unwrap(), Point::ORIGIN);
    near(b.high().unwrap(), Point::new(30.0, 30.0, 5.0));
    holds_its_faces(&model, &plate);

    let ball = make_sphere(&mut model, Frame::WORLD, 5.0, T).unwrap().shape;
    let any = explore_unique(&model, &ball, ShapeType::Face).unwrap()[0].clone();
    assert_eq!(model.face_bounds(&any), None);
    holds_its_faces(&model, &ball);
}

#[test]
fn a_turned_occurrence_is_bounded_where_it_stands() {
    let mut model = Model::new();
    let pts =
        [(4.0, 0.0), (6.0, 0.0), (6.0, 10.0), (4.0, 10.0)].map(|(x, y)| Point::new(x, y, 0.0));
    let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
    let face = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Y,
    };
    let tube = make_revolution(&mut model, &face, axis, core::f64::consts::TAU, T)
        .unwrap()
        .shape;
    tight_bounds(&model, &tube, T).unwrap();
    // Turned about its own axis the tube fills the same box, which the box
    // of its turned box would overshoot by the diagonal.
    let turned = model.placed(
        &tube,
        Transform::rotation(axis, core::f64::consts::FRAC_PI_4),
    );
    let b = tight_bounds(&model, &turned, T).unwrap();
    near(b.low().unwrap(), Point::new(-6.0, 0.0, -6.0));
    near(b.high().unwrap(), Point::new(6.0, 10.0, 6.0));
}

#[test]
fn a_face_moved_by_move_faces_reports_its_new_box() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    tight_bounds(&model, &block, T).unwrap();
    let top = explore_unique(&model, &block, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| face_bounds(&model, f).unwrap().low().unwrap().z > 9.0)
        .unwrap();
    let lift = Transform::translation(Vector::new(0.0, 0.0, 2.0));
    let raised = move_faces(&mut model, &block, &[top], &lift, T)
        .unwrap()
        .shape;
    let b = tight_bounds(&model, &raised, T).unwrap();
    near(b.high().unwrap(), Point::new(10.0, 10.0, 12.0));
    let tops: Vec<_> = explore_unique(&model, &raised, ShapeType::Face)
        .unwrap()
        .into_iter()
        .map(|f| face_bounds(&model, &f).unwrap())
        .filter(|b| b.low().unwrap().z > 11.0)
        .collect();
    assert_eq!(tops.len(), 1);
    near(tops[0].high().unwrap(), Point::new(10.0, 10.0, 12.0));
    let b = tight_bounds(&model, &block, T).unwrap();
    near(b.high().unwrap(), Point::new(10.0, 10.0, 10.0));
}

/// The plate of 582 faces: a box with a 24 by 24 grid of holes.
#[test]
#[ignore = "heavy"]
fn a_second_tight_bounds_of_a_drilled_plate_costs_a_tenth_of_the_first() {
    let mut model = Model::new();
    let plate = drilled_plate(&mut model, 24, 10.0, 10.0);
    assert_eq!(
        explore_unique(&model, &plate, ShapeType::Face)
            .unwrap()
            .len(),
        582
    );
    let time = |model: &Model| {
        let started = ogeom::core::clock::Instant::now();
        let b = tight_bounds(model, &plate, T).unwrap();
        near(b.high().unwrap(), Point::new(250.0, 250.0, 10.0));
        started.elapsed().as_secs_f64()
    };
    // The least of a few runs each way, every first one on a copy whose
    // faces have forgotten their boxes: the boolean keeps them, and a face
    // handed out for editing forgets its own.
    let forgotten = || {
        let mut copy = model.clone();
        for face in explore_unique(&model, &plate, ShapeType::Face).unwrap() {
            let _ = copy.node_mut(&face);
        }
        copy
    };
    let first = (0..5)
        .map(|_| time(&forgotten()))
        .fold(f64::INFINITY, f64::min);
    time(&model);
    let second = (0..5).map(|_| time(&model)).fold(f64::INFINITY, f64::min);
    assert!(
        second * 10.0 < first,
        "first {first:.6} s, second {second:.6} s"
    );
    holds_its_faces(&model, &plate);
}
