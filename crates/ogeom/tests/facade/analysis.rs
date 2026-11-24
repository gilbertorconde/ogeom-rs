//! Self-intersection: a shape asked whether it crosses itself, and answering
//! with the crossings rather than with a flag.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::{LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Plane, Point};
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn crossing_sheets_confess_and_a_box_does_not() {
    let mut model = Model::new();
    let box_shape = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    assert!(
        ogeom::algo::check_self_intersection(&model, &box_shape, T)
            .unwrap()
            .is_empty(),
        "a valid solid interferes with nothing"
    );

    // Two square sheets crossing like an X, sharing no topology.
    let sheet = |model: &mut Model, frame: Frame| -> Shape {
        let corners = [
            frame.origin() + frame.x().vector() * -5.0 + frame.y().vector() * -5.0,
            frame.origin() + frame.x().vector() * 5.0 + frame.y().vector() * -5.0,
            frame.origin() + frame.x().vector() * 5.0 + frame.y().vector() * 5.0,
            frame.origin() + frame.x().vector() * -5.0 + frame.y().vector() * 5.0,
        ];
        let vertices: Vec<Shape> = corners
            .iter()
            .map(|c| ogeom::algo::make_vertex(model, *c).shape)
            .collect();
        let edges: Vec<Shape> = (0..4)
            .map(|i| {
                let (a, b) = (corners[i], corners[(i + 1) % 4]);
                ogeom::algo::make_edge_between(
                    model,
                    LineCurve::segment(a, b, T).unwrap().into(),
                    (0.0, a.distance(b)),
                    &vertices[i],
                    &vertices[(i + 1) % 4],
                    T,
                )
                .unwrap()
                .shape
            })
            .collect();
        ogeom::algo::make_face_with_pcurves(
            model,
            SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(frame))),
            &[edges],
            T,
        )
        .unwrap()
        .shape
    };
    let flat = sheet(&mut model, Frame::WORLD);
    let tilted_frame = Frame::new(
        Point::new(0.0, 0.0, -2.0),
        Direction::from_coords(0.0, 1.0, 0.2, T).unwrap(),
        Direction::X,
        T,
    )
    .unwrap();
    let tilted = sheet(&mut model, tilted_frame);
    let pair = ogeom::algo::make_compound(&mut model, &[flat, tilted])
        .unwrap()
        .shape;

    let crossings = ogeom::algo::check_self_intersection(&model, &pair, T).unwrap();
    assert_eq!(crossings.len(), 1, "the sheets cross once: {crossings:?}");
}
