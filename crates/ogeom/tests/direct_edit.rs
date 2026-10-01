//! Editing a solid's faces directly: some faces offset or moved, the faces
//! around them following, the solid staying closed.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    Built, check, face_normal, make_box, make_cylinder, make_sphere, volume_properties,
};
use ogeom::core::{OgeomResult, Tolerances};
use ogeom::geom::{Surface as _, SurfaceKind};
use ogeom::math::{Axis, Direction, Frame, Point, Transform, Vector};
use ogeom::mesh::Deflection;
use ogeom::offset::{move_faces, offset_faces};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// The faces whose middle and outward normal pass `pick`.
fn faces_where(model: &Model, shape: &Shape, pick: impl Fn(Point, Vector) -> bool) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            let (p, n) = face_normal(model, f, T).unwrap();
            pick(p, n)
        })
        .collect()
}

/// The faces on surfaces of `kind`.
fn faces_of_kind(model: &Model, shape: &Shape, kind: SurfaceKind) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            model.geometry().surface(data.surface).unwrap().kind() == kind
        })
        .collect()
}

/// The edit's result, valid and of volume `want`, with `faces` faces.
fn sound(model: &Model, built: OgeomResult<Built>, want: f64, faces: usize, at: &str) -> Built {
    let built = built.unwrap_or_else(|e| panic!("{at}: {e}"));
    let diagnosis = check(model, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{at}: {diagnosis}");
    let v = volume(model, &built.shape);
    assert!((v - want).abs() < want * 1e-9, "{at}: {v} against {want}");
    let count = explore_unique(model, &built.shape, ShapeType::Face)
        .unwrap()
        .len();
    assert_eq!(count, faces, "{at}");
    built
}

fn block(model: &mut Model) -> Shape {
    make_box(model, Frame::WORLD, (20.0, 10.0, 5.0), T)
        .unwrap()
        .shape
}

/// A box's top offset out and in: the sides follow, and every face maps to
/// the face it became.
#[test]
fn a_box_top_offset_moves_its_sides_with_it() {
    let mut model = Model::new();
    let block = block(&mut model);
    let top = faces_where(&model, &block, |p, n| (p.z - 5.0).abs() < 1e-9 && n.z > 0.5);
    for (distance, want) in [(2.0, 1400.0), (-2.0, 600.0)] {
        let edited = offset_faces(&mut model, &block, &top, distance, T);
        let built = sound(&model, edited, want, 6, &format!("top by {distance}"));
        for face in explore_unique(&model, &block, ShapeType::Face).unwrap() {
            assert_eq!(built.history.trace(&face).len(), 1);
        }
    }
}

/// A side moved along its normal is the side offset by as much.
#[test]
fn moving_a_face_along_its_normal_is_offsetting_it() {
    let mut model = Model::new();
    let block = block(&mut model);
    let side = faces_where(&model, &block, |p, n| {
        (p.x - 20.0).abs() < 1e-9 && n.x > 0.5
    });
    let offset = offset_faces(&mut model, &block, &side, 3.0, T);
    sound(&model, offset, 1150.0, 6, "offset");
    let shift = Transform::translation(Vector::new(3.0, 0.0, 0.0));
    let moved = move_faces(&mut model, &block, &side, &shift, T);
    sound(&model, moved, 1150.0, 6, "moved");
}

/// A through bore offset into the material and out of it: its radius
/// changes, and the planes it pierces follow.
#[test]
fn a_bore_offset_changes_its_radius() {
    let mut model = Model::new();
    let block = block(&mut model);
    let seat = Frame::new(Point::new(10.0, 5.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let pin = make_cylinder(&mut model, seat, 2.0, 7.0, T).unwrap().shape;
    let bored = ogeom::boolean::cut(&mut model, &block, &pin, T)
        .unwrap()
        .shape;
    let bore = faces_of_kind(&model, &bored, SurfaceKind::Cylinder);
    // Its outward normal points at the axis: into the material is wider.
    for (distance, radius) in [(-0.5, 2.5), (0.5, 1.5)] {
        let edited = offset_faces(&mut model, &bored, &bore, distance, T);
        let want = 1000.0 - PI * radius * radius * 5.0;
        sound(&model, edited, want, 7, &format!("bore by {distance}"));
    }
}

/// A boss on a block: its top raised, its wall thickened, and the whole
/// boss shifted along the block, the block's top following each time.
#[test]
fn a_boss_is_raised_thickened_and_shifted() {
    let mut model = Model::new();
    let block = block(&mut model);
    let seat = Frame::new(Point::new(10.0, 5.0, 5.0), Direction::Z, Direction::X, T).unwrap();
    let post = make_cylinder(&mut model, seat, 2.0, 3.0, T).unwrap().shape;
    let bossed = ogeom::boolean::fuse(&mut model, &block, &post, T)
        .unwrap()
        .shape;
    let top = faces_where(&model, &bossed, |p, n| {
        (p.z - 8.0).abs() < 1e-9 && n.z > 0.5
    });
    let wall = faces_of_kind(&model, &bossed, SurfaceKind::Cylinder);
    let raised = offset_faces(&mut model, &bossed, &top, 1.0, T);
    sound(&model, raised, 1000.0 + PI * 4.0 * 4.0, 8, "raised");
    let thicker = offset_faces(&mut model, &bossed, &wall, 0.5, T);
    sound(&model, thicker, 1000.0 + PI * 6.25 * 3.0, 8, "thickened");
    let boss: Vec<Shape> = top.iter().chain(&wall).cloned().collect();
    let shift = Transform::translation(Vector::new(2.0, 0.0, 0.0));
    let shifted = move_faces(&mut model, &bossed, &boss, &shift, T);
    sound(&model, shifted, 1000.0 + PI * 4.0 * 3.0, 8, "shifted");
}

/// A top turned about a line through its middle: the volume it gains on
/// one side it loses on the other.
#[test]
fn a_face_turned_about_its_middle_keeps_the_volume() {
    let mut model = Model::new();
    let block = block(&mut model);
    let top = faces_where(&model, &block, |p, n| (p.z - 5.0).abs() < 1e-9 && n.z > 0.5);
    let hinge = Axis::new(Point::new(0.0, 5.0, 5.0), Direction::X);
    let turn = Transform::rotation(hinge, 5f64.to_radians());
    let turned = move_faces(&mut model, &block, &top, &turn, T);
    sound(&model, turned, 1000.0, 6, "turned");
}

/// A ball, its one face offset: a larger ball.
#[test]
fn a_ball_offset_is_a_larger_ball() {
    let mut model = Model::new();
    let ball = make_sphere(&mut model, Frame::WORLD, 3.0, T).unwrap().shape;
    let face = explore_unique(&model, &ball, ShapeType::Face).unwrap();
    let grown = offset_faces(&mut model, &ball, &face, 1.0, T);
    sound(&model, grown, 4.0 / 3.0 * PI * 64.0, 1, "grown");
}

/// A face driven past the one across from it, or into the rest of the
/// solid, is refused by name, and so is a motion that is not rigid.
#[test]
fn an_edit_that_breaks_the_solid_is_refused() {
    let mut model = Model::new();
    let block = block(&mut model);
    let top = faces_where(&model, &block, |p, n| (p.z - 5.0).abs() < 1e-9 && n.z > 0.5);
    let through = offset_faces(&mut model, &block, &top, -6.0, T).unwrap_err();
    assert!(through.to_string().contains("inside out"), "{through}");

    let seat = Frame::new(Point::new(10.0, 5.0, 5.0), Direction::Z, Direction::X, T).unwrap();
    let post = make_cylinder(&mut model, seat, 2.0, 3.0, T).unwrap().shape;
    let bossed = ogeom::boolean::fuse(&mut model, &block, &post, T)
        .unwrap()
        .shape;
    let boss: Vec<Shape> = faces_where(&model, &bossed, |p, n| {
        (p.z - 8.0).abs() < 1e-9 && n.z > 0.5
    })
    .into_iter()
    .chain(faces_of_kind(&model, &bossed, SurfaceKind::Cylinder))
    .collect();
    let off_the_edge = Transform::translation(Vector::new(9.0, 0.0, 0.0));
    let crossing = move_faces(&mut model, &bossed, &boss, &off_the_edge, T).unwrap_err();
    assert!(crossing.to_string().contains("cross itself"), "{crossing}");

    let scale = Transform::scaling(Point::ORIGIN, 2.0, T).unwrap();
    let reshaped = move_faces(&mut model, &block, &top, &scale, T).unwrap_err();
    assert!(
        reshaped.to_string().contains("translation or a rotation"),
        "{reshaped}"
    );
}
