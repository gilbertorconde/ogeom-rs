//! Circular pipes along wires of several edges, and profiles brought onto
//! their spine before a pipe shell sweeps them, measured against closed
//! forms.
//!
//! A disc of radius `r` square to a spine, its centre on it, holds by
//! Pappus `π r²` times the spine's length, mitred corners included: each
//! leg is a cylinder cut on planes through its axis's ends.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, make_edge, make_edge_between, make_face, make_polygon, make_vertex, make_wire,
    volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve3d as _, LineCurve, PlaneSurface};
use ogeom::math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::offset::{PipePlacement, make_pipe, make_pipe_shell, place_on_spine};
use ogeom::topo::{Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// The volume and centroid, and whether they were integrated exactly.
fn measure(model: &Model, shape: &Shape) -> (f64, Point, bool) {
    let p = volume_properties(model, shape, Deflection::default(), T).unwrap();
    (p.mass, p.centre, p.deflection == 0.0)
}

/// Up `z` to `(0, 0, 5)`, a quarter arc of radius 3 over to `(3, 0, 8)`
/// heading `+x`, and on along `x` to `(8, 0, 8)`, each edge tangent to the
/// next; then, where `turn`, a sharp quarter turn along `y` to
/// `(8, 6, 8)`.
fn line_arc_line(model: &mut Model, turn: bool) -> Shape {
    let points = [
        Point::ORIGIN,
        Point::new(0.0, 0.0, 5.0),
        Point::new(3.0, 0.0, 8.0),
        Point::new(8.0, 0.0, 8.0),
        Point::new(8.0, 6.0, 8.0),
    ];
    let v: Vec<Shape> = points
        .iter()
        .map(|p| make_vertex(model, *p).shape)
        .collect();
    let segment = |model: &mut Model, i: usize| {
        let s = LineCurve::segment(points[i], points[i + 1], T).unwrap();
        let range = s.domain();
        make_edge_between(model, s.into(), range, &v[i], &v[i + 1], T)
            .unwrap()
            .shape
    };
    let frame = Frame::new(Point::new(3.0, 0.0, 5.0), Direction::Y, -Direction::X, T).unwrap();
    let arc = CircleCurve::new(Circle::new(frame, 3.0, T).unwrap());
    let mut edges = vec![segment(model, 0)];
    edges.push(
        make_edge_between(model, arc.into(), (0.0, PI / 2.0), &v[1], &v[2], T)
            .unwrap()
            .shape,
    );
    edges.push(segment(model, 2));
    if turn {
        edges.push(segment(model, 3));
    }
    make_wire(model, &edges, T).unwrap().shape
}

fn analytic_only(model: &Model, shape: &Shape) -> bool {
    use ogeom::geom::SurfaceGeometry as S;
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .iter()
        .all(|f| {
            let Some(NodeData::Face(d)) = model.node(f).map(|n| n.data()) else {
                return false;
            };
            matches!(
                model.geometry().surface(d.surface),
                Some(S::Plane(_) | S::Cylinder(_) | S::Torus(_))
            )
        })
}

/// A pipe of radius `r` along three straight legs, 10, 10 and 7 long,
/// mitred at both corners: `π r² · 27`, exactly. Along a line, a quarter
/// arc of radius 3 and a line, meeting tangent: `π r² (10 + 3π/2)`,
/// exactly, by Pappus.
#[test]
fn a_pipe_along_a_wire_holds_pappus_s_volume_exactly() {
    let r = 0.5;
    let mut model = Model::new();
    let legs = make_polygon(
        &mut model,
        &[
            Point::ORIGIN,
            Point::new(0.0, 0.0, 10.0),
            Point::new(10.0, 0.0, 10.0),
            Point::new(10.0, 7.0, 10.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    let smooth = line_arc_line(&mut model, false);
    for (spine, length) in [(legs, 27.0), (smooth, 10.0 + 1.5 * PI)] {
        let pipe = make_pipe(&mut model, &spine, r, T).unwrap().shape;
        let found = check(&model, &pipe, T).unwrap();
        assert!(found.is_valid(), "{found}");
        assert!(analytic_only(&model, &pipe));
        let (v, _, exact) = measure(&model, &pipe);
        assert!(exact);
        let want = PI * r * r * length;
        assert!((v - want).abs() < want * 1e-12, "{v} against {want}");
    }
}

#[test]
fn a_pipe_along_a_wire_refuses_what_it_does_not_sweep_exactly_by_name() {
    let mut model = Model::new();
    let turned = line_arc_line(&mut model, true);
    let err = make_pipe(&mut model, &turned, 0.5, T).unwrap_err();
    assert!(err.to_string().contains("turns no sharp corner"), "{err}");
    let tight = line_arc_line(&mut model, false);
    let err = make_pipe(&mut model, &tight, 3.0, T).unwrap_err();
    assert!(err.to_string().contains("swallows its spine"), "{err}");
}

/// A `2 × 2` square drawn about `(5, 3, 1)` on the plane facing
/// `(1, 1, 1)`, brought onto a spine up `z` from the origin by contact and
/// correction: the prism `4 · 4`, its centroid half way up the spine, at
/// `(0, 0, 2)`, exactly.
#[test]
fn a_leaning_square_off_the_spine_is_brought_onto_it_and_swept() {
    let mut model = Model::new();
    let at = Point::new(5.0, 3.0, 1.0);
    let normal = Direction::new(Vector::new(1.0, 1.0, 1.0), T).unwrap();
    let frame = Frame::about(at, normal);
    let (x, y) = (frame.x().vector(), frame.y().vector());
    let corners =
        [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)].map(|(a, b)| at + x * a + y * b);
    let ring = make_polygon(&mut model, &corners, true, T).unwrap().shape;
    let square = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &[ring],
        T,
    )
    .unwrap()
    .shape;
    let line = LineCurve::segment(Point::ORIGIN, Point::new(0.0, 0.0, 4.0), T).unwrap();
    let range = line.domain();
    let spine = make_edge(&mut model, line.into(), range, T).unwrap().shape;

    let leaning = make_pipe_shell(&mut model, &square, &spine, false, 1e-3, T);
    assert!(leaning.is_err(), "a leaning square off the spine swept");
    let placed = place_on_spine(
        &mut model,
        &square,
        &spine,
        PipePlacement {
            contact: true,
            correction: true,
        },
        T,
    )
    .unwrap()
    .shape;
    let pipe = make_pipe_shell(&mut model, &placed, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    let found = check(&model, &pipe, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let (v, centre, exact) = measure(&model, &pipe);
    assert!(exact);
    assert!((v - 16.0).abs() < 16.0 * 1e-12, "{v}");
    assert!(
        centre.distance(Point::new(0.0, 0.0, 2.0)) < 10.0 * T.confusion(),
        "{centre:?}"
    );
}

/// Contact alone moves a wire unturned: a `2 × 4` rectangle about
/// `(3, 4, 0)` in the `xy` plane lands with its corners at `(±1, ±2, 0)`.
/// Correction alone turns a square leaning a sixth of a turn about `x`
/// through the spine's start, its centroid there already, square to the
/// spine up `z`: the prism `4 · 3`, exactly.
#[test]
fn contact_and_correction_each_move_the_profile_as_named() {
    let mut model = Model::new();
    let line = LineCurve::segment(Point::ORIGIN, Point::new(0.0, 0.0, 3.0), T).unwrap();
    let range = line.domain();
    let spine = make_edge(&mut model, line.into(), range, T).unwrap().shape;

    let rectangle = make_polygon(
        &mut model,
        &[(2.0, 2.0), (4.0, 2.0), (4.0, 6.0), (2.0, 6.0)].map(|(x, y)| Point::new(x, y, 0.0)),
        true,
        T,
    )
    .unwrap()
    .shape;
    let moved = place_on_spine(
        &mut model,
        &rectangle,
        &spine,
        PipePlacement {
            contact: true,
            correction: false,
        },
        T,
    )
    .unwrap()
    .shape;
    assert_eq!(model.kind_of(&moved).unwrap(), ShapeType::Wire);
    for v in explore_unique(&model, &moved, ShapeType::Vertex).unwrap() {
        let Some(NodeData::Vertex(d)) = model.node(&v).map(|n| n.data()) else {
            panic!("a vertex holds no point");
        };
        let p = d.point;
        assert!(
            (p.x.abs() - 1.0).abs() < 1e-12 && (p.y.abs() - 2.0).abs() < 1e-12 && p.z.abs() < 1e-12,
            "{p:?}"
        );
    }

    let tilt = PI / 6.0;
    let normal = Direction::new(Vector::new(0.0, -tilt.sin(), tilt.cos()), T).unwrap();
    let frame = Frame::new(Point::ORIGIN, normal, Direction::X, T).unwrap();
    let (x, y) = (frame.x().vector(), frame.y().vector());
    let corners = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
        .map(|(a, b)| Point::ORIGIN + x * a + y * b);
    let ring = make_polygon(&mut model, &corners, true, T).unwrap().shape;
    let square = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &[ring],
        T,
    )
    .unwrap()
    .shape;
    let turned = place_on_spine(
        &mut model,
        &square,
        &spine,
        PipePlacement {
            contact: false,
            correction: true,
        },
        T,
    )
    .unwrap()
    .shape;
    let pipe = make_pipe_shell(&mut model, &turned, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    let (v, _, exact) = measure(&model, &pipe);
    assert!(exact);
    assert!((v - 12.0).abs() < 12.0 * 1e-12, "{v}");
}
