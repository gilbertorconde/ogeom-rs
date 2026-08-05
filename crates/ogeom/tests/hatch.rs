//! Hatching a trimmed face: the lines stop at the boundary and skip the
//! hole, and their total length is the area divided by the spacing —
//! Cavalieri doing the bookkeeping.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::{LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Plane, Point};
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn square_edges(model: &mut Model, corners: [Point; 4]) -> Vec<Shape> {
    let vertices: Vec<Shape> = corners
        .iter()
        .map(|c| ogeom::algo::make_vertex(model, *c).shape)
        .collect();
    (0..4)
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
        .collect()
}

#[test]
fn hatch_lines_stop_at_the_boundary_and_skip_the_hole() {
    let mut model = Model::new();
    let outer = square_edges(
        &mut model,
        [
            Point::new(0.0, 0.0, 0.0),
            Point::new(10.0, 0.0, 0.0),
            Point::new(10.0, 10.0, 0.0),
            Point::new(0.0, 10.0, 0.0),
        ],
    );
    let hole = square_edges(
        &mut model,
        [
            Point::new(4.0, 4.0, 0.0),
            Point::new(6.0, 4.0, 0.0),
            Point::new(6.0, 6.0, 0.0),
            Point::new(4.0, 6.0, 0.0),
        ],
    );
    let face = ogeom::algo::make_face_with_pcurves(
        &mut model,
        SurfaceGeometry::Plane(PlaneSurface::new(Plane::XY)),
        &[outer, hole],
        T,
    )
    .unwrap()
    .shape;

    let segments = ogeom::mesh::hatch_face(
        &model,
        &face,
        0.0,
        1.0,
        ogeom::mesh::Deflection::default(),
        T,
    )
    .unwrap();

    // Ten scanlines at half-offsets; the two through the hole band split in
    // two, so twelve segments whose lengths sum to the face's area.
    assert_eq!(segments.len(), 12, "{segments:?}");
    let total: f64 = segments.iter().map(|[a, b]| a.distance(*b)).sum();
    assert!(
        (total - 96.0).abs() < 1e-9,
        "hatch length {total} against area 96"
    );

    // At forty-five degrees the count changes but Cavalieri holds within
    // the sampling of the corners.
    let slanted = ogeom::mesh::hatch_face(
        &model,
        &face,
        core::f64::consts::FRAC_PI_4,
        0.25,
        ogeom::mesh::Deflection::default(),
        T,
    )
    .unwrap();
    let slant_total: f64 = slanted.iter().map(|[a, b]| a.distance(*b)).sum();
    assert!(
        (slant_total * 0.25 - 96.0).abs() < 3.0,
        "slanted hatch estimates the area: {}",
        slant_total * 0.25
    );
}
