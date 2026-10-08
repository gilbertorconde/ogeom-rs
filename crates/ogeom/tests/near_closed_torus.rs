//! A torus whose hole nearly closes costs what any torus costs to mesh, to
//! check and to build by revolution: the parallel round the inner equator
//! shrinks with the hole, so a fixed chord and angle ask about as many
//! points of it as of an open torus.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{PI, TAU};
use ogeom::core::clock::Instant;

use ogeom::algo::{
    check, make_edge, make_face, make_revolution, make_torus, make_wire, volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve, PlaneSurface};
use ogeom::math::{Axis, Circle, Direction, Frame, Plane, Point};
use ogeom::mesh::{Deflection, triangulate_face};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const MINOR: f64 = 5.0;

/// The disc a torus of major radius `major` is the revolution of: radius
/// [`MINOR`], centred at `(major, 0, 0)` in the XZ plane.
fn disc(model: &mut Model, major: f64) -> Shape {
    let frame = Frame::new(Point::new(major, 0.0, 0.0), Direction::Y, Direction::X, T).unwrap();
    let circle = Circle::new(frame, MINOR, T).unwrap();
    let edge = make_edge(
        model,
        Curve::Circle(CircleCurve::new(circle)),
        (0.0, TAU),
        T,
    )
    .unwrap()
    .shape;
    let wire = make_wire(model, &[edge], T).unwrap().shape;
    make_face(
        model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape
}

/// The torus's own volume, `2 pi^2 R r^2`.
fn torus_volume(major: f64) -> f64 {
    2.0 * PI * PI * major * MINOR * MINOR
}

fn assert_valid_with_volume(model: &Model, solid: &Shape, major: f64, what: &str) {
    let report = check(model, solid, T).unwrap();
    assert!(report.is_valid(), "{what}: {report:?}");
    let volume = volume_properties(model, solid, Deflection::default(), T)
        .unwrap()
        .mass;
    let want = torus_volume(major);
    assert!(
        (volume - want).abs() <= want * 0.01,
        "{what}: volume {volume}, want {want}"
    );
}

#[test]
fn a_torus_whose_hole_nearly_closes_meshes_like_an_open_one() {
    let axes = Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
    for gap in [0.1, 0.01, 0.001] {
        let major = MINOR + gap;
        let mut model = Model::new();
        let torus = make_torus(&mut model, axes, major, MINOR, T).unwrap().shape;
        let face = explore_unique(&model, &torus, ShapeType::Face)
            .unwrap()
            .remove(0);
        let mesh = triangulate_face(&model, &face, Deflection::default(), T).unwrap();
        assert!(
            mesh.triangles.len() < 20_000,
            "gap {gap}: {} triangles",
            mesh.triangles.len()
        );
        // A polyhedron inscribed at the default deflection falls short of
        // the true area by about a percent, as an open torus's does.
        let area = 4.0 * PI * PI * major * MINOR;
        assert!(
            (mesh.area() - area).abs() <= area * 0.02,
            "gap {gap}: area {}, want {area}",
            mesh.area()
        );
    }
}

#[test]
fn a_nearly_closed_torus_checks_and_revolves_in_good_time() {
    let axes = Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
    for gap in [0.1, 0.01, 0.001] {
        let major = MINOR + gap;
        let started = Instant::now();
        let mut model = Model::new();
        let torus = make_torus(&mut model, axes, major, MINOR, T).unwrap().shape;
        assert_valid_with_volume(&model, &torus, major, &format!("torus, gap {gap}"));

        let mut model = Model::new();
        let face = disc(&mut model, major);
        let axis = Axis::new(Point::ORIGIN, Direction::Z);
        let solid = make_revolution(&mut model, &face, axis, TAU, T)
            .unwrap()
            .shape;
        assert_valid_with_volume(&model, &solid, major, &format!("revolution, gap {gap}"));
        let took = started.elapsed();
        assert!(took.as_secs_f64() < 1.0, "gap {gap}: took {took:?}");
    }
}
