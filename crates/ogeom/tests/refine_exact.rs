//! Refining a solid never turns an exact measurement into an approximate
//! one: the merged face of coplanar planes is a plane face, measured in
//! closed form like the faces it replaced.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{make_box, make_cylinder, volume_properties};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn exact_volume(model: &Model, shape: &Shape) -> f64 {
    let p = volume_properties(model, shape, Deflection::default(), T).unwrap();
    assert_eq!(p.deflection, 0.0, "measured from a mesh: {}", p.mass);
    p.mass
}

/// A pad: the rectangle at height `z` swept up by `h`.
fn pad(model: &mut Model, x: f64, y: f64, z: f64, h: f64) -> Shape {
    use ogeom::algo::{make_face, make_polygon, make_prism};
    let pts = [(0.0, 0.0), (x, 0.0), (x, y), (0.0, y)].map(|(a, b)| Point::new(a, b, z));
    let wire = make_polygon(model, &pts, true, T).unwrap().shape;
    let frame = Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X, T).unwrap();
    let plane = ogeom::geom::PlaneSurface::new(ogeom::math::Plane::new(frame));
    let face = make_face(model, plane.into(), &[wire], T).unwrap().shape;
    make_prism(model, &face, ogeom::math::Vector::new(0.0, 0.0, h), T)
        .unwrap()
        .shape
}

#[test]
fn a_refined_stack_of_pads_measures_exactly() {
    let mut model = Model::new();
    let base = pad(&mut model, 20.0, 20.0, 0.0, 10.0);
    let step = pad(&mut model, 5.0, 5.0, 10.0, 5.0);
    let fused = ogeom::boolean::fuse(&mut model, &base, &step, T)
        .unwrap()
        .shape;
    let refined = ogeom::heal::unify_same_domain(&mut model, &fused, T)
        .unwrap()
        .0
        .shape;
    let v = exact_volume(&model, &refined);
    assert!((v - 4125.0).abs() < 1e-9, "{v}");
}

#[test]
fn a_refined_solid_measures_exactly() {
    let mut model = Model::new();
    let base = make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(0.0, 0.0, 10.0), Direction::Z, Direction::X, T).unwrap();
    let step = make_box(&mut model, at, (5.0, 5.0, 5.0), T).unwrap().shape;
    let fused = ogeom::boolean::fuse(&mut model, &base, &step, T)
        .unwrap()
        .shape;
    let refined = ogeom::heal::unify_same_domain(&mut model, &fused, T)
        .unwrap()
        .0
        .shape;
    let v = exact_volume(&model, &refined);
    assert!((v - 4125.0).abs() < 1e-9, "{v}");

    let side = Frame::new(Point::new(20.0, 10.0, 5.0), Direction::X, Direction::Y, T).unwrap();
    let post = make_cylinder(&mut model, side, 2.0, 3.0, T).unwrap().shape;
    let more = ogeom::boolean::fuse(&mut model, &fused, &post, T)
        .unwrap()
        .shape;
    let refined = ogeom::heal::unify_same_domain(&mut model, &more, T)
        .unwrap()
        .0
        .shape;
    let want = 4125.0 + core::f64::consts::PI * 4.0 * 3.0;
    let v = exact_volume(&model, &refined);
    assert!((v - want).abs() < want * 1e-9, "{v} against {want}");
}
