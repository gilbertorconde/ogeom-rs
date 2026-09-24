//! Mass properties of faces the closed-form rectangles and discs do not
//! cover: an elliptic wall, a plane bounded by an ellipse or by half of
//! one. Each is integrated round its own chart boundary, so the answer
//! does not depend on the deflection passed in.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    make_edge, make_edge_between, make_face_with_pcurves, make_prism, surface_properties,
    volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{Curve, Curve3d, EllipseCurve, LineCurve, PlaneSurface};
use ogeom::math::{Ellipse, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, VertexData};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn ellipse(a: f64, b: f64) -> Curve {
    Curve::Ellipse(EllipseCurve::new(
        Ellipse::new(Frame::WORLD, a, b, T).unwrap(),
    ))
}

fn prism_of(model: &mut Model, boundary: Vec<Shape>, height: f64) -> Shape {
    let plane = PlaneSurface::over(Plane::XY, (-20.0, 20.0), (-20.0, 20.0)).unwrap();
    let face = make_face_with_pcurves(model, plane.into(), &[boundary], T)
        .unwrap()
        .shape;
    make_prism(model, &face, Vector::new(0.0, 0.0, height), T)
        .unwrap()
        .shape
}

/// The perimeter of an ellipse, by adaptive quadrature of its speed.
fn perimeter(a: f64, b: f64) -> f64 {
    ogeom::math::integrate(
        |t: f64| (a * t.sin()).hypot(b * t.cos()),
        0.0,
        2.0 * PI,
        1e-13,
    )
    .unwrap()
}

/// A 20 by 10 ellipse padded 4 high measures to rounding at the default
/// deflection, and says it was not meshed.
#[test]
fn an_elliptic_prism_measures_within_a_part_in_a_thousand() {
    let mut model = Model::new();
    let curve = ellipse(10.0, 5.0);
    let range = curve.domain();
    let edge = make_edge(&mut model, curve, range, T).unwrap().shape;
    let prism = prism_of(&mut model, vec![edge], 4.0);

    let v = volume_properties(&model, &prism, Deflection::default(), T).unwrap();
    let exact = PI * 10.0 * 5.0 * 4.0;
    assert!((v.mass - exact).abs() < 1e-9 * exact, "{}", v.mass);
    assert_eq!(v.deflection, 0.0, "integrated, not meshed");
    let centre = v.centre;
    assert!(
        centre.distance(Point::new(0.0, 0.0, 2.0)) < 1e-9,
        "{centre:?}"
    );

    let s = surface_properties(&model, &prism, Deflection::default(), T).unwrap();
    let area = 2.0 * PI * 10.0 * 5.0 + perimeter(10.0, 5.0) * 4.0;
    assert!(
        (s.mass - area).abs() < 1e-9 * area,
        "{} against {area}",
        s.mass
    );
}

/// Half of that ellipse, closed along its major axis: an arc and a line
/// on the plane, an elliptic wall and a flat one round the side.
#[test]
fn half_an_elliptic_prism_measures_half() {
    let mut model = Model::new();
    let (east, west) = (Point::new(10.0, 0.0, 0.0), Point::new(-10.0, 0.0, 0.0));
    let (e, w) = (
        model.add_vertex(VertexData::new(east)),
        model.add_vertex(VertexData::new(west)),
    );
    let arc = make_edge_between(&mut model, ellipse(10.0, 5.0), (0.0, PI), &e, &w, T)
        .unwrap()
        .shape;
    let chord = LineCurve::segment(west, east, T).unwrap();
    let range = chord.domain();
    let base = make_edge_between(&mut model, Curve::Line(chord), range, &w, &e, T)
        .unwrap()
        .shape;
    let half = prism_of(&mut model, vec![arc, base], 4.0);

    let v = volume_properties(&model, &half, Deflection::default(), T).unwrap();
    let exact = PI * 10.0 * 5.0 * 4.0 / 2.0;
    assert!((v.mass - exact).abs() < 1e-9 * exact, "{}", v.mass);
    assert_eq!(v.deflection, 0.0);
    // The centroid of a half ellipse stands 4b/(3 pi) off its major axis.
    let centre = v.centre;
    let y = 4.0 * 5.0 / (3.0 * PI);
    assert!(
        centre.distance(Point::new(0.0, y, 2.0)) < 1e-9,
        "{centre:?}"
    );
}

/// A drum with a quarter cut away, along its seam and across its axis:
/// its caps are concave, three quarters of a disc, and its wall keeps the
/// seam down one column only. It measures exactly, which needs every face
/// read the right way out.
#[test]
fn a_three_quarter_drum_measures_exactly() {
    use ogeom::math::Direction;
    let mut model = Model::new();
    let below = Frame::new(Point::new(0.0, 0.0, -3.0), Direction::Z, Direction::X, T).unwrap();
    let drum = ogeom::algo::make_cylinder(&mut model, below, 5.0, 6.0, T)
        .unwrap()
        .shape;
    let corner = Frame::new(Point::new(0.0, 0.0, -20.0), Direction::Z, Direction::X, T).unwrap();
    let quarter = ogeom::algo::make_box(&mut model, corner, (20.0, 20.0, 40.0), T)
        .unwrap()
        .shape;
    let rest = ogeom::boolean::cut(&mut model, &drum, &quarter, T)
        .unwrap()
        .shape;
    let v = volume_properties(&model, &rest, Deflection::default(), T).unwrap();
    let exact = 0.75 * PI * 25.0 * 6.0;
    assert_eq!(v.deflection, 0.0, "integrated, not meshed");
    assert!((v.mass - exact).abs() < 1e-9 * exact, "{}", v.mass);
}
