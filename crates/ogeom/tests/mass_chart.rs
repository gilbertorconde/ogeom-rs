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

/// The area of the cone wall of a frustum of radius 6 at its base and 3 at
/// its top, 10 high, that a drill of radius `r` lying along y on its base
/// through its rim takes away: the part of the wall within `r` of the line
/// `x = 6, z = r`. At height `z` the wall's circle of radius `rho` is inside
/// the drill where `rho cos(theta)` exceeds `6 - s`, `s` the drill's half
/// width there, up to the height `2r / 1.09` where the drill leaves the
/// wall. Its square-root ends are smoothed by a substitution at each.
fn drilled_cone_wall(r: f64) -> f64 {
    let slant = 1.09_f64.sqrt();
    let width = |z: f64| {
        let rho = 6.0 - 0.3 * z;
        let s = (2.0 * r * z - z * z).max(0.0).sqrt();
        2.0 * ((6.0 - s) / rho).clamp(-1.0, 1.0).acos() * rho * slant
    };
    let top = 2.0 * r / 1.09;
    let half = 0.5 * top;
    let low = ogeom::math::integrate(
        |t| width(half * t.powi(4)) * 4.0 * half * t.powi(3),
        0.0,
        1.0,
        1e-13,
    )
    .unwrap();
    let high = ogeom::math::integrate(
        |t| width(top - half * t * t) * 2.0 * half * t,
        0.0,
        1.0,
        1e-13,
    )
    .unwrap();
    low + high
}

/// A frustum drilled by a drill lying on its base, the drill's wall
/// touching the base's circle at the cone wall's seam. The section on the
/// wall is fitted, its pieces' ends a few hundred-thousandths off the
/// vertices they meet, and the cut's wall wraps the whole chart round the
/// hole. The cut's wall measures the lateral area less the hole's, the
/// common's wall the hole's, and cut and common add up to the frustum to
/// a part in a billion.
#[test]
fn a_frustum_drilled_through_its_rim_measures_its_hole() {
    use ogeom::algo::{make_cone, make_cylinder};
    use ogeom::math::Direction;
    use ogeom::topo::{NodeData, ShapeType, explore_unique};
    let fine = Deflection::with_chord(1e-3).unwrap();
    let cone_walls = |model: &Model, shape: &Shape| -> f64 {
        explore_unique(model, shape, ShapeType::Face)
            .unwrap()
            .iter()
            .filter(|face| {
                let NodeData::Face(data) = model.node(face).unwrap().data() else {
                    return false;
                };
                matches!(
                    model.geometry().surface(data.surface),
                    Some(ogeom::geom::SurfaceGeometry::Cone(_))
                )
            })
            .map(|face| surface_properties(model, face, fine, T).unwrap().mass)
            .sum()
    };
    let lateral = PI * 9.0 * 109.0_f64.sqrt();
    let frustum = PI * 10.0 / 3.0 * (36.0 + 18.0 + 9.0);
    for r in [0.3, 0.5] {
        let mut model = Model::new();
        let part = make_cone(&mut model, Frame::WORLD, 6.0, 3.0, 10.0, T)
            .unwrap()
            .shape;
        let frame = Frame::new(
            Point::new(6.0, -8.0, r),
            Direction::Y,
            Direction::new(Vector::new(0.0, 0.0, -1.0), T).unwrap(),
            T,
        )
        .unwrap();
        let drill = make_cylinder(&mut model, frame, r, 16.0, T).unwrap().shape;
        let cut = ogeom::boolean::cut(&mut model, &part, &drill, T)
            .unwrap()
            .shape;
        let common = ogeom::boolean::common(&mut model, &part, &drill, T)
            .unwrap()
            .shape;

        let hole = drilled_cone_wall(r);
        let (kept, taken) = (cone_walls(&model, &cut), cone_walls(&model, &common));
        // The two walls share the section's pcurves, so between them they
        // cover the lateral area to rounding.
        assert!(
            (kept + taken - lateral).abs() < 1e-9 * lateral,
            "r {r}: the walls {kept} + {taken} against {lateral}"
        );
        // Each against the hole: the section is fitted to within 2e-5 of
        // the true one along its 2 mm or so, which moves the hole's area by
        // a few millionths.
        assert!(
            (kept - (lateral - hole)).abs() < 5e-6,
            "r {r}: the cut's wall {kept} against {}",
            lateral - hole
        );
        assert!(
            (taken - hole).abs() < 5e-6,
            "r {r}: the common's wall {taken} against {hole}"
        );

        let (v_cut, v_common) = (
            volume_properties(&model, &cut, fine, T).unwrap(),
            volume_properties(&model, &common, fine, T).unwrap(),
        );
        assert_eq!(v_cut.deflection, 0.0, "integrated, not meshed");
        assert_eq!(v_common.deflection, 0.0, "integrated, not meshed");
        let sum = v_cut.mass + v_common.mass;
        assert!(
            (sum - frustum).abs() < 1e-9 * frustum,
            "r {r}: cut {} + common {} against {frustum}",
            v_cut.mass,
            v_common.mass
        );
    }
}
