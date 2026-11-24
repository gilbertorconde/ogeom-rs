//! Faces the blend engine builds round a closed surface, blended again: a
//! bore's wall and top after its rim is chamfered, a post's wall and the
//! block's top after its foot is rounded. Each face's rings must close in
//! its surface's chart for the face blend's cuts to find them.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::SurfaceKind;
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A frame at `(x, y, z)` with its axis along `axis`.
fn frame(x: f64, y: f64, z: f64, axis: Vector) -> Frame {
    let axis = Direction::new(axis, T).unwrap();
    let across = Direction::new(axis.vector().cross(Vector::new(0.3, 0.5, 0.7)), T).unwrap();
    Frame::new(Point::new(x, y, z), axis, across, T).unwrap()
}

/// A leaning axis, so no chart is the one a world-aligned solid gets.
fn leaning() -> Vector {
    Vector::new(0.300_882, -0.814_548, -0.495_965)
}

fn kind(m: &Model, face: &Shape) -> SurfaceKind {
    use ogeom_geom::Surface as _;
    let data = m.node(face).unwrap().data().as_face().unwrap();
    m.geometry().surface(data.surface).unwrap().kind()
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// The face of `want` holding the point `at`, a point of its surface.
fn face_at(m: &Model, solid: &Shape, want: SurfaceKind, at: Point) -> Shape {
    let found: Vec<Shape> = explore_unique(m, solid, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| kind(m, f) == want)
        .filter(|f| {
            ogeom_algo::classify_on_face(m, f, at, ogeom_mesh::Deflection::default(), T).unwrap()
                == ogeom_algo::Containment::In
        })
        .collect();
    let [face] = found.as_slice() else {
        panic!("{} faces of {want:?} at {at:?}", found.len());
    };
    face.clone()
}

/// The circular edges between `a` and `b`.
fn rim(m: &Model, a: &Shape, b: &Shape) -> Vec<Shape> {
    let of_b = explore_unique(m, b, ShapeType::Edge).unwrap();
    explore_unique(m, a, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| of_b.iter().any(|x| x.is_same(e)))
        .collect()
}

/// A block 10 by 10 by 3 on a leaning axis, bored through with radius 1.5,
/// its top rim chamfered 0.1; then the bore's wall and the top blended at
/// 0.6. The ball rolls in the corner the chamfer leaves between them, and
/// the round takes off more than the chamfer did.
#[test]
fn a_chamfered_bore_s_wall_and_top_blend_again() {
    let mut m = Model::new();
    let base = frame(0.7, 1.4, -2.6, leaning());
    let axis = base.z().vector();
    let o = base.origin();
    let (x, y) = (base.x().vector(), base.y().vector());
    let corner = o - x * 5.0 - y * 5.0;
    let block_frame = Frame::new(corner, base.z(), base.x(), T).unwrap();
    let block = ogeom_algo::make_box(&mut m, block_frame, (10.0, 10.0, 3.0), T)
        .unwrap()
        .shape;
    let drill = ogeom_algo::make_cylinder(
        &mut m,
        Frame::new(o - axis * 1.0, base.z(), base.x(), T).unwrap(),
        1.5,
        5.0,
        T,
    )
    .unwrap()
    .shape;
    let bored = ogeom_bool::cut(&mut m, &block, &drill, T).unwrap().shape;
    let top_at = o + axis * 3.0;
    let wall = face_at(
        &m,
        &bored,
        SurfaceKind::Cylinder,
        o + axis * 1.5 + (x * 0.6 + y * 0.8) * 1.5,
    );
    let top = face_at(&m, &bored, SurfaceKind::Plane, top_at + x * 3.0);
    let edges = rim(&m, &wall, &top);
    assert!(!edges.is_empty());
    let chamfered = ogeom_fillet::chamfer_edges(&mut m, &bored, &edges, 0.1, T)
        .unwrap()
        .shape;
    let before = volume(&m, &chamfered);
    let wall = face_at(
        &m,
        &chamfered,
        SurfaceKind::Cylinder,
        o + axis * 1.5 + (x * 0.6 + y * 0.8) * 1.5,
    );
    let top = face_at(&m, &chamfered, SurfaceKind::Plane, top_at + x * 3.0);
    let blended = ogeom_fillet::blend_faces(&mut m, &chamfered, &wall, &top, 0.6, T).unwrap();
    assert!(ogeom_algo::check(&m, &blended.shape, T).unwrap().is_valid());
    let after = volume(&m, &blended.shape);
    // The round of 0.6 over the sharp corner takes off (1 - pi/4) 0.36 of
    // section, the chamfer half of 0.01; turned about the axis at a radius
    // past the bore's.
    assert!(after < before, "{after} against {before}");
    let section = (1.0 - core::f64::consts::FRAC_PI_4) * 0.36 - 0.005;
    let most = section * core::f64::consts::TAU * (1.5 + 0.6);
    assert!(before - after < most, "{} past {most}", before - after);
}

/// A post of radius 2, 3.5 tall, fused on a block 12 by 8 by 1 on a
/// leaning axis, its foot rounded 0.1; then the post's wall and the block's
/// top blended at 0.25. The ball rolls in the corner the round leaves, and
/// the blend fills more than the round did.
#[test]
fn a_rounded_post_s_foot_blends_again() {
    let mut m = Model::new();
    let base = frame(-2.4, 1.9, 2.5, leaning());
    let axis = base.z().vector();
    let o = base.origin();
    let (x, y) = (base.x().vector(), base.y().vector());
    let corner = o - x * 6.0 - y * 4.0;
    let block_frame = Frame::new(corner, base.z(), base.x(), T).unwrap();
    let block = ogeom_algo::make_box(&mut m, block_frame, (12.0, 8.0, 1.0), T)
        .unwrap()
        .shape;
    let post = ogeom_algo::make_cylinder(
        &mut m,
        Frame::new(o + axis * 0.5, base.z(), base.x(), T).unwrap(),
        2.0,
        4.0,
        T,
    )
    .unwrap()
    .shape;
    let solid = ogeom_bool::fuse(&mut m, &block, &post, T).unwrap().shape;
    let foot_at = o + axis * 1.0;
    let wall = face_at(
        &m,
        &solid,
        SurfaceKind::Cylinder,
        foot_at + axis * 1.0 + (x * 0.6 + y * 0.8) * 2.0,
    );
    let top = face_at(&m, &solid, SurfaceKind::Plane, foot_at + x * 4.0);
    let edges = rim(&m, &wall, &top);
    assert!(!edges.is_empty());
    let rounded = ogeom_fillet::fillet_edges(&mut m, &solid, &edges, 0.1, T)
        .unwrap()
        .shape;
    let before = volume(&m, &rounded);
    let wall = face_at(
        &m,
        &rounded,
        SurfaceKind::Cylinder,
        foot_at + axis * 2.0 + (x * 0.6 + y * 0.8) * 2.0,
    );
    let top = face_at(&m, &rounded, SurfaceKind::Plane, foot_at + x * 4.0);
    let blended = ogeom_fillet::blend_faces(&mut m, &rounded, &wall, &top, 0.25, T).unwrap();
    assert!(ogeom_algo::check(&m, &blended.shape, T).unwrap().is_valid());
    let after = volume(&m, &blended.shape);
    // The fill of 0.25 in the sharp corner adds (1 - pi/4) 0.0625 of
    // section, the round of 0.1 had (1 - pi/4) 0.01; turned about the axis
    // at a radius past the post's.
    assert!(after > before, "{after} against {before}");
    let section = (1.0 - core::f64::consts::FRAC_PI_4) * (0.0625 - 0.01);
    let most = section * core::f64::consts::TAU * (2.0 + 0.25);
    assert!(after - before < most, "{} past {most}", after - before);
}
