//! §J of `docs/PLAN.md`: the medial axis, held to closed forms — the
//! rectangle's roof line and the triangle's incenter.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Frame, Point};
use ogeom_topo::Model;

const T: Tolerances = Tolerances::millimetres();

fn face_over(model: &mut Model, corners: &[Point]) -> ogeom_topo::Shape {
    ogeom_algo::make_polygon(model, corners, true, T)
        .unwrap()
        .shape
}

fn planar_face(model: &mut Model, corners: &[Point]) -> ogeom_topo::Shape {
    use ogeom_geom::{PlaneSurface, SurfaceGeometry};
    use ogeom_math::Plane;
    let wire = face_over(model, corners);
    let edges = ogeom_topo::explore(
        model,
        &wire,
        ogeom_topo::Filter::OfType(ogeom_topo::ShapeType::Edge),
    )
    .unwrap();
    ogeom_algo::make_face_with_pcurves(
        model,
        SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(Frame::WORLD))),
        &[edges],
        T,
    )
    .unwrap()
    .shape
}

/// A 20×10 rectangle's axis: four corner diagonals at forty-five degrees
/// meeting the roof line, whose clearance is the half height everywhere.
#[test]
fn a_rectangles_axis_is_its_roof_line() {
    let mut model = Model::new();
    let face = planar_face(
        &mut model,
        &[
            Point::new(0.0, 0.0, 0.0),
            Point::new(20.0, 0.0, 0.0),
            Point::new(20.0, 10.0, 0.0),
            Point::new(0.0, 10.0, 0.0),
        ],
    );
    let axis = ogeom_algo::medial_axis(&model, &face, T).unwrap();
    // Four diagonals plus the ridge.
    assert_eq!(axis.segments.len(), 5, "{:?}", axis.segments);
    let total: f64 = axis.segments.iter().map(|(a, b)| a.distance(*b)).sum();
    let want = 4.0 * 5.0 * core::f64::consts::SQRT_2 + 10.0;
    assert!((total - want).abs() < 1e-9, "{total} against {want}");
    // The ridge runs between the two meets at half height.
    let meets = [Point::new(5.0, 5.0, 0.0), Point::new(15.0, 5.0, 0.0)];
    for meet in meets {
        assert!(
            axis.segments
                .iter()
                .any(|(a, b)| a.distance(meet) < 1e-9 || b.distance(meet) < 1e-9),
            "a meet at {meet:?}"
        );
    }
    let deepest = axis.clearance.iter().copied().fold(0.0, f64::max);
    assert!((deepest - 5.0).abs() < 1e-9, "half the height: {deepest}");
}

/// A 3–4–5 right triangle's axis meets at the incenter, and the clearance
/// there is the inradius, (a + b − c) / 2 = 1 — a closed form older than
/// the field.
#[test]
fn a_triangles_axis_meets_at_the_incenter() {
    let mut model = Model::new();
    let face = planar_face(
        &mut model,
        &[
            Point::new(0.0, 0.0, 0.0),
            Point::new(4.0, 0.0, 0.0),
            Point::new(0.0, 3.0, 0.0),
        ],
    );
    let axis = ogeom_algo::medial_axis(&model, &face, T).unwrap();
    let incenter = Point::new(1.0, 1.0, 0.0);
    for (a, b) in &axis.segments {
        assert!(
            a.distance(incenter) < 1e-9 || b.distance(incenter) < 1e-9,
            "every branch reaches the incenter: {a:?} -> {b:?}"
        );
    }
    let deepest = axis.clearance.iter().copied().fold(0.0, f64::max);
    assert!((deepest - 1.0).abs() < 1e-9, "the inradius: {deepest}");
}

/// The refusals, by name: holes, reflex corners and arcs each change the
/// mathematics, and a wrong axis gouges quietly.
#[test]
fn what_the_axis_cannot_yet_hold_is_refused_by_name() {
    let mut model = Model::new();
    let face = planar_face(
        &mut model,
        &[
            Point::new(0.0, 0.0, 0.0),
            Point::new(10.0, 0.0, 0.0),
            Point::new(10.0, 10.0, 0.0),
            Point::new(5.0, 3.0, 0.0), // reflex
            Point::new(0.0, 10.0, 0.0),
        ],
    );
    let err = ogeom_algo::medial_axis(&model, &face, T).unwrap_err();
    assert!(err.to_string().contains("reflex"), "{err}");
}
