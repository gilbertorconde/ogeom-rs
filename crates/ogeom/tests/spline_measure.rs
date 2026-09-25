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
