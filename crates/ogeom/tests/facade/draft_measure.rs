//! A wall drafted by a small angle measures like one drafted by a large.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{make_face, make_polygon, make_prism, surface_properties, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::PlaneSurface;
use ogeom::math::{Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_box_drafted_a_degree_and_a_half_measures() {
    for degrees in [1.5_f64, 10.0] {
        let mut model = Model::new();
        let pts = [(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0)]
            .map(|(x, y)| Point::new(x, y, 0.0));
        let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
        let face = make_face(
            &mut model,
            PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
            &[wire],
            T,
        )
        .unwrap()
        .shape;
        let block = make_prism(&mut model, &face, Vector::new(0.0, 0.0, 10.0), T)
            .unwrap()
            .shape;
        let wall = explore_unique(&model, &block, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                let b = ogeom::algo::shape_bounds(&model, f, T).unwrap();
                b.high().unwrap().x < 1e-3
            })
            .expect("the wall at x = 0");
        let neutral = Plane::through(Point::new(10.0, 10.0, 0.0), -Direction::Z);
        let drafted = ogeom::offset::apply_draft(
            &mut model,
            &block,
            &[wall],
            neutral,
            Direction::Z,
            degrees.to_radians(),
            T,
        )
        .unwrap()
        .shape;
        surface_properties(&model, &drafted, Deflection::default(), T).unwrap();
        let v = volume_properties(&model, &drafted, Deflection::default(), T)
            .unwrap()
            .mass;
        let wedge = 0.5 * 10.0 * 10.0 * degrees.to_radians().tan() * 20.0;
        assert!(
            (v - (4000.0 - wedge)).abs() < 1e-3 || (v - (4000.0 + wedge)).abs() < 1e-3,
            "{degrees}: {v}"
        );
    }
}
