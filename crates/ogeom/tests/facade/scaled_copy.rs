//! A uniformly scaled copy is the original grown: a cube 1.5 times the
//! size, whole whether it stands apart from the original or overlaps it.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, general_transformed_shape, make_face, make_polygon, make_prism, shape_bounds,
    volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::PlaneSurface;
use ogeom::math::{Frame, GeneralTransform, Matrix3, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    let d = check(model, shape, T).unwrap();
    assert!(d.is_valid(), "{:?}", d.problems);
    volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

#[test]
fn a_scaled_copy_standing_apart_is_whole() {
    let mut model = Model::new();
    // A pad: the square swept up, its far cap its near cap moved.
    let pts = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)].map(|(x, y)| Point::new(x, y, 0.0));
    let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
    let face = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let cube = make_prism(&mut model, &face, Vector::new(0.0, 0.0, 4.0), T)
        .unwrap()
        .shape;
    let grow = GeneralTransform {
        linear: Matrix3::scaling(1.5),
        translation: Vector::new(10.0, 0.0, 0.0),
    };
    let copy = general_transformed_shape(&mut model, &cube, &grow, T)
        .unwrap()
        .shape;
    let v = volume(&model, &copy);
    assert!((v - 216.0).abs() < 216.0 * 1e-9, "the copy alone: {v}");
    let both = ogeom::boolean::fuse(&mut model, &cube, &copy, T)
        .unwrap()
        .shape;
    let v = volume(&model, &both);
    assert!((v - 280.0).abs() < 280.0 * 1e-4, "{v}");
    let b = shape_bounds(&model, &both, T).unwrap();
    let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
    assert!(
        lo.x.abs() < 1e-3 && (hi.x - 16.0).abs() < 1e-3 && (hi.y - 6.0).abs() < 1e-3,
        "{lo:?} {hi:?}"
    );
}

#[test]
fn a_pad_placed_at_one_and_a_half_times_measures_so() {
    let mut model = Model::new();
    let pts = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)].map(|(x, y)| Point::new(x, y, 0.0));
    let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
    let face = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let pad = make_prism(&mut model, &face, Vector::new(0.0, 0.0, 4.0), T)
        .unwrap()
        .shape;
    let grow = ogeom::math::Transform::from_parts(
        Matrix3::IDENTITY,
        1.5,
        Vector::new(10.0, 0.0, 0.0),
        1e-12,
    )
    .unwrap();
    let placed = ogeom::algo::transformed(&mut model, &pad, grow)
        .unwrap()
        .shape;
    let v = volume_properties(&model, &placed, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((v - 216.0).abs() < 216.0 * 1e-9, "{v}");
}
