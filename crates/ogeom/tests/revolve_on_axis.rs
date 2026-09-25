//! A profile touching its axis revolves through any angle: the edge on the
//! axis is its own image, shared by the start and end faces.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_face, make_polygon, make_revolution, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::PlaneSurface;
use ogeom::math::{Axis, Direction, Frame, Plane, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::Model;

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_rectangle_on_its_axis_revolves_three_quarters() {
    for degrees in [90.0_f64, 270.0, 360.0] {
        let mut model = Model::new();
        let pts =
            [(0.0, 0.0), (5.0, 0.0), (5.0, 10.0), (0.0, 10.0)].map(|(x, y)| Point::new(x, y, 0.0));
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
        let solid = make_revolution(&mut model, &face, axis, degrees.to_radians(), T)
            .unwrap_or_else(|e| panic!("{degrees}: {e}"))
            .shape;
        let diagnosis = check(&model, &solid, T).unwrap();
        assert!(diagnosis.is_valid(), "{degrees}: {:?}", diagnosis.problems);
        let v = volume_properties(&model, &solid, Deflection::default(), T)
            .unwrap()
            .mass;
        let want = degrees / 360.0 * core::f64::consts::PI * 25.0 * 10.0;
        assert!(
            (v - want).abs() < want * 1e-4,
            "{degrees}: {v} against {want}"
        );
    }
}

#[test]
fn a_trapezoid_on_its_axis_revolves_to_part_of_a_cone() {
    let mut model = Model::new();
    let pts =
        [(0.0, 0.0), (5.0, 0.0), (3.0, 10.0), (0.0, 10.0)].map(|(x, y)| Point::new(x, y, 0.0));
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
    let solid = make_revolution(&mut model, &face, axis, 200f64.to_radians(), T)
        .unwrap()
        .shape;
    let diagnosis = check(&model, &solid, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let v = volume_properties(&model, &solid, Deflection::default(), T)
        .unwrap()
        .mass;
    let frustum = core::f64::consts::PI * 10.0 / 3.0 * (25.0 + 15.0 + 9.0);
    let want = frustum * 200.0 / 360.0;
    assert!((v - want).abs() < want * 1e-4, "{v} against {want}");
}
