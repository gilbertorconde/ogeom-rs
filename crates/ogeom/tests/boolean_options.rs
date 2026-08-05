//! The boolean's API completions: scaled placements bake automatically and
//! run the analytic pipeline.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Transform};
use ogeom::mesh::Deflection;
use ogeom::topo::Model;

const T: Tolerances = Tolerances::millimetres();

/// A placement with scale no longer refuses: the boolean bakes it — every
/// surface restated exactly in its own analytic vocabulary, pcurves
/// re-derived — and the drilled result measures.
#[test]
fn a_scaled_placement_bakes_and_the_boolean_runs() {
    let mut model = Model::new();
    let unit = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let scaled = ogeom::algo::transformed(
        &mut model,
        &unit,
        Transform::scaling(Point::ORIGIN, 2.0, T).unwrap(),
    )
    .unwrap()
    .shape;

    // The bake alone is exact: the doubled box measures eight thousand.
    let baked = ogeom::algo::baked_shape(&mut model, &scaled, T)
        .unwrap()
        .shape;
    let v = ogeom::algo::volume_properties(&model, &baked, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((v - 8000.0).abs() < 1e-6, "the bake is exact: {v}");

    // And the boolean takes the scaled shape directly.
    let f = Frame::new(Point::new(10.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, f, 3.0, 22.0, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(&mut model, &scaled, &drill, T)
        .unwrap()
        .shape;
    let v = ogeom::algo::volume_properties(&model, &cut, Deflection::with_chord(1e-2).unwrap(), T)
        .unwrap()
        .mass;
    let expected = 8000.0 - core::f64::consts::PI * 9.0 * 20.0;
    assert!(
        (v - expected).abs() / expected < 1e-3,
        "the scaled cut measures: {v} vs {expected}"
    );
}

/// A half space finally has its consumer: cutting a box with the solid on
/// one side of a plane slices it exactly at the plane, and the union with
/// a half space refuses as the unbounded thing it would be.
#[test]
fn a_half_space_cuts_and_refuses_to_fuse() {
    use ogeom::geom::{LineCurve, PlaneSurface, SurfaceGeometry};
    use ogeom::math::Plane;
    use ogeom::topo::Shape;

    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;

    // A face on the plane z = 4, and the half space above it.
    let corners = [
        Point::new(-5.0, -5.0, 4.0),
        Point::new(15.0, -5.0, 4.0),
        Point::new(15.0, 15.0, 4.0),
        Point::new(-5.0, 15.0, 4.0),
    ];
    let vertices: Vec<Shape> = corners
        .iter()
        .map(|c| ogeom::algo::make_vertex(&mut model, *c).shape)
        .collect();
    let edges: Vec<Shape> = (0..4)
        .map(|i| {
            let (a, b) = (corners[i], corners[(i + 1) % 4]);
            ogeom::algo::make_edge_between(
                &mut model,
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
    let plane_frame = Frame::new(Point::new(0.0, 0.0, 4.0), Direction::Z, Direction::X, T).unwrap();
    let face = ogeom::algo::make_face_with_pcurves(
        &mut model,
        SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(plane_frame))),
        &[edges],
        T,
    )
    .unwrap()
    .shape;
    let upper = ogeom::algo::make_half_space(&mut model, &face, Point::new(0.0, 0.0, 100.0), T)
        .unwrap()
        .shape;

    // Cutting away everything above z = 4 leaves the 10x10x4 slab.
    let slab = ogeom::boolean::cut(&mut model, &block, &upper, T)
        .unwrap()
        .shape;
    let v = ogeom::algo::volume_properties(&model, &slab, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((v - 400.0).abs() < 1e-6, "the slab below the plane: {v}");

    // Common keeps the upper part.
    let cap = ogeom::boolean::common(&mut model, &block, &upper, T)
        .unwrap()
        .shape;
    let v = ogeom::algo::volume_properties(&model, &cap, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((v - 600.0).abs() < 1e-6, "the cap above the plane: {v}");

    // And the union is refused as unbounded, by name.
    let err = ogeom::boolean::fuse(&mut model, &block, &upper, T).unwrap_err();
    assert!(err.to_string().contains("unbounded"), "{err}");
}
