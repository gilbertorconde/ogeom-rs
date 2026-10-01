//! Spline faces with any trim measure exactly on their surfaces: a box
//! converted to splines and drilled keeps its volume to rounding, without
//! a mesh.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{make_box, make_cylinder, to_nurbs, volume_properties};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::Model;

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_drilled_spline_box_measures_exactly() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let block = to_nurbs(&mut model, &block, T).unwrap().shape;
    let at = Frame::new(Point::new(5.0, 5.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = make_cylinder(&mut model, at, 2.0, 12.0, T).unwrap().shape;
    let holed = ogeom::boolean::cut(&mut model, &block, &drill, T)
        .unwrap()
        .shape;
    let p = volume_properties(&model, &holed, Deflection::default(), T).unwrap();
    let want = 1000.0 - core::f64::consts::PI * 4.0 * 10.0;
    assert_eq!(p.deflection, 0.0, "measured from a mesh: {}", p.mass);
    assert!(
        (p.mass - want).abs() < want * 1e-6,
        "{} against {want}",
        p.mass
    );
}

/// A pointed cone converted to splines, or stretched unevenly, keeps its
/// apex: the conversion's chart stops there rather than running on into the
/// cone's other nappe.
#[test]
fn a_pointed_cone_converts_up_to_its_apex() {
    let pointed = std::f64::consts::PI * 36.0 * 9.0 / 3.0;
    let mut model = Model::new();
    let cone = ogeom::algo::make_cone(&mut model, Frame::WORLD, 6.0, 0.0, 9.0, T)
        .unwrap()
        .shape;
    let splines = to_nurbs(&mut model, &cone, T).unwrap().shape;
    let fine = Deflection::with_chord(1e-3).unwrap();
    let measured = volume_properties(&model, &splines, fine, T).unwrap().mass;
    assert!(
        (measured - pointed).abs() < pointed * 1e-6,
        "{measured} against {pointed}"
    );
    let stretch = ogeom::math::GeneralTransform::scaling_xyz(2.0, 1.0, 1.0);
    let stretched = ogeom::algo::general_transformed_shape(&mut model, &cone, &stretch, T)
        .unwrap()
        .shape;
    let measured = volume_properties(&model, &stretched, fine, T).unwrap().mass;
    assert!(
        (measured - 2.0 * pointed).abs() < pointed * 1e-6,
        "{measured} against {}",
        2.0 * pointed
    );
    let top = ogeom::algo::tight_bounds(&model, &stretched, T)
        .unwrap()
        .high()
        .unwrap()
        .z;
    assert!((top - 9.0).abs() < 1e-6, "the apex stands at {top}");
}
