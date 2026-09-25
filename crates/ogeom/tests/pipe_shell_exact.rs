//! Pipe shells that have closed forms: a flat profile edge down a straight
//! leg is a plane, a circle down a line a drum; a ring walked by reversed
//! edges sweeps as the same ring walked forward; a profile square to a
//! curved spine's start sweeps along it.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{make_edge, make_face, make_polygon, make_wire, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve, Curve3d as _, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Circle, Direction, Frame, Plane, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn line(model: &mut Model, a: Point, b: Point) -> Shape {
    let curve = LineCurve::segment(a, b, T).unwrap();
    let range = curve.domain();
    make_edge(model, Curve::Line(curve), range, T)
        .unwrap()
        .shape
}

fn xy_face(model: &mut Model, wire: Shape) -> Shape {
    let plane = PlaneSurface::new(Plane::new(Frame::WORLD));
    make_face(model, plane.into(), &[wire], T).unwrap().shape
}

fn square(model: &mut Model, half: f64) -> Shape {
    let pts = [(-half, -half), (half, -half), (half, half), (-half, half)]
        .map(|(x, y)| Point::new(x, y, 0.0));
    let wire = make_polygon(model, &pts, true, T).unwrap().shape;
    xy_face(model, wire)
}

/// Volume, and whether it was integrated exactly.
fn measure(model: &Model, shape: &Shape) -> (f64, bool) {
    let p = volume_properties(model, shape, Deflection::default(), T).unwrap();
    (p.mass, p.deflection == 0.0)
}

fn planes_only(model: &Model, shape: &Shape) -> bool {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .iter()
        .all(|f| {
            let Some(NodeData::Face(d)) = model.node(f).map(|n| n.data()) else {
                return false;
            };
            matches!(
                model.geometry().surface(d.surface),
                Some(SurfaceGeometry::Plane(_))
            )
        })
}

#[test]
fn a_square_down_an_l_is_planes_and_measures_exactly() {
    let mut model = Model::new();
    let profile = square(&mut model, 2.0);
    let spine = make_polygon(
        &mut model,
        &[
            Point::ORIGIN,
            Point::new(0.0, 0.0, 20.0),
            Point::new(20.0, 0.0, 20.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    let pipe = ogeom::offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    assert!(planes_only(&model, &pipe));
    let (v, exact) = measure(&model, &pipe);
    assert!(exact, "measured from a mesh");
    assert!((v - 640.0).abs() < 640.0 * 1e-9, "{v}");
}

#[test]
fn a_ring_of_reversed_edges_sweeps_as_the_forward_ring() {
    let mut model = Model::new();
    let pts =
        [(-2.0, -2.0), (2.0, -2.0), (2.0, 2.0), (-2.0, 2.0)].map(|(x, y)| Point::new(x, y, 0.0));
    let forward = make_polygon(&mut model, &pts, true, T).unwrap().shape;
    let edges = model.ordered_children_of(&forward).unwrap();
    let reversed: Vec<Shape> = edges.iter().rev().map(Shape::reversed).collect();
    let wire = make_wire(&mut model, &reversed, T).unwrap().shape;
    let profile = xy_face(&mut model, wire);
    let e = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, 20.0));
    let spine = make_wire(&mut model, &[e], T).unwrap().shape;
    let pipe = ogeom::offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    let (v, _) = measure(&model, &pipe);
    assert!((v - 320.0).abs() < 320.0 * 1e-9, "{v}");
}

#[test]
fn a_circle_down_a_line_is_a_drum() {
    let mut model = Model::new();
    let circle = Circle::new(Frame::WORLD, 2.0, T).unwrap();
    let e = make_edge(
        &mut model,
        CircleCurve::new(circle).into(),
        (0.0, core::f64::consts::TAU),
        T,
    )
    .unwrap()
    .shape;
    let wire = make_wire(&mut model, &[e], T).unwrap().shape;
    let profile = xy_face(&mut model, wire);
    let s = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, 20.0));
    let spine = make_wire(&mut model, &[s], T).unwrap().shape;
    let pipe = ogeom::offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    let (v, exact) = measure(&model, &pipe);
    let want = core::f64::consts::PI * 4.0 * 20.0;
    assert!(exact, "measured from a mesh");
    assert!((v - want).abs() < want * 1e-9, "{v} against {want}");
}

#[test]
fn a_square_square_to_an_arcs_start_sweeps_along_it() {
    let quarter = core::f64::consts::FRAC_PI_2;
    // A quarter circle from the origin to (10, 0, 10), centred on
    // (10, 0, 0) in the XZ plane: its start tangent is +Z. Walked forward,
    // and as the same arc stored the other way round and used reversed.
    for backward in [false, true] {
        let mut model = Model::new();
        let profile = square(&mut model, 2.0);
        let (x, range) = if backward {
            (Direction::Z, (0.0, quarter))
        } else {
            (-Direction::X, (0.0, quarter))
        };
        let frame = Frame::new(
            Point::new(10.0, 0.0, 0.0),
            if backward {
                -Direction::Y
            } else {
                Direction::Y
            },
            x,
            T,
        )
        .unwrap();
        let circle = Circle::new(frame, 10.0, T).unwrap();
        let edge = make_edge(&mut model, CircleCurve::new(circle).into(), range, T)
            .unwrap()
            .shape;
        let edge = if backward { edge.reversed() } else { edge };
        let spine = make_wire(&mut model, &[edge], T).unwrap().shape;
        let pipe = ogeom::offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T)
            .unwrap_or_else(|e| panic!("backward {backward}: {e}"))
            .shape;
        let (v, exact) = measure(&model, &pipe);
        let want = 16.0 * quarter * 10.0;
        assert!(exact, "measured from a mesh");
        assert!((v - want).abs() < want * 1e-9, "{v} against {want}");
    }
}

#[test]
fn a_profile_square_to_a_conical_helix_sweeps() {
    let mut model = Model::new();
    let frame = Frame::WORLD;
    let taper = 5.0 * 10f64.to_radians().tan();
    let helix = ogeom::geom::HelixCurve::conical(
        frame,
        10.0,
        5.0,
        taper,
        0.0,
        2.0 * core::f64::consts::TAU,
    )
    .unwrap();
    let range = helix.domain();
    let start = helix.point_at(range.0, T).unwrap();
    let tangent = helix.d1_at(range.0, T).unwrap();
    let edge = make_edge(&mut model, helix.into(), range, T).unwrap().shape;
    let spine = make_wire(&mut model, &[edge], T).unwrap().shape;
    // A small square square to the start tangent, centred on the start.
    let normal = Direction::new(tangent, T).unwrap();
    let plane = Plane::through(start, normal);
    let (u, w) = (plane.frame().x().vector(), plane.frame().y().vector());
    let pts: Vec<Point> = [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]
        .iter()
        .map(|(a, b)| start + u * *a + w * *b)
        .collect();
    let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
    let face = make_face(&mut model, PlaneSurface::new(plane).into(), &[wire], T)
        .unwrap()
        .shape;
    ogeom::offset::make_pipe_shell(&mut model, &face, &spine, false, 1e-3, T).unwrap();
}

#[test]
fn a_placed_profile_square_to_an_arc_sweeps_along_it() {
    let mut model = Model::new();
    let profile = square(&mut model, 2.0);
    let profile = ogeom::algo::transformed(&mut model, &profile, ogeom::math::Transform::IDENTITY)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(10.0, 0.0, 0.0), Direction::Y, -Direction::X, T).unwrap();
    let circle = Circle::new(frame, 10.0, T).unwrap();
    let quarter = core::f64::consts::FRAC_PI_2;
    let edge = make_edge(
        &mut model,
        CircleCurve::new(circle).into(),
        (0.0, quarter),
        T,
    )
    .unwrap()
    .shape;
    let spine = make_wire(&mut model, &[edge], T).unwrap().shape;
    let pipe = ogeom::offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    let (v, exact) = measure(&model, &pipe);
    let want = 16.0 * quarter * 10.0;
    assert!(exact);
    assert!((v - want).abs() < want * 1e-9, "{v} against {want}");
}

#[test]
/// A square, holed or not, round a closed ellipse under either frame law:
/// centred on the spine, it sweeps its area times the ellipse's length, the
/// hole a void tunnel whichever way its ring is wound.
fn a_holed_profile_rounds_an_elliptic_ring() {
    let length = {
        let n = 20_000;
        (0..n)
            .map(|i| {
                let t = core::f64::consts::TAU * (f64::from(i) + 0.5) / f64::from(n);
                (20.0 * t.sin()).hypot(15.0 * t.cos())
            })
            .sum::<f64>()
            * core::f64::consts::TAU
            / f64::from(n)
    };
    for frenet in [false, true] {
        for holed in [false, true] {
            let mut model = Model::new();
            // An ellipse 20 by 15 in the XY plane, closed.
            let ellipse = ogeom::math::Ellipse::new(Frame::WORLD, 20.0, 15.0, T).unwrap();
            let e = make_edge(
                &mut model,
                ogeom::geom::EllipseCurve::new(ellipse).into(),
                (0.0, core::f64::consts::TAU),
                T,
            )
            .unwrap()
            .shape;
            let spine = make_wire(&mut model, &[e], T).unwrap().shape;
            // A 4 x 4 square in the XZ plane at (20, 0, 0), square to +Y.
            let pts = [(-2.0, -2.0), (2.0, -2.0), (2.0, 2.0), (-2.0, 2.0)]
                .map(|(a, b)| Point::new(20.0 + a, 0.0, b));
            let outer = make_polygon(&mut model, &pts, true, T).unwrap().shape;
            let frame =
                Frame::new(Point::new(20.0, 0.0, 0.0), Direction::Y, Direction::X, T).unwrap();
            let mut wires = vec![outer];
            if holed {
                let c = Circle::new(frame, 1.0, T).unwrap();
                let he = make_edge(
                    &mut model,
                    CircleCurve::new(c).into(),
                    (0.0, core::f64::consts::TAU),
                    T,
                )
                .unwrap()
                .shape;
                wires.push(make_wire(&mut model, &[he], T).unwrap().shape);
            }
            let face = make_face(
                &mut model,
                PlaneSurface::new(Plane::new(frame)).into(),
                &wires,
                T,
            )
            .unwrap()
            .shape;
            let r = ogeom::offset::make_pipe_shell(&mut model, &face, &spine, frenet, 1e-3, T);
            let pipe = r
                .unwrap_or_else(|e| panic!("frenet {frenet} holed {holed}: {e}"))
                .shape;
            let (v, _) = measure(&model, &pipe);
            let area = if holed {
                16.0 - core::f64::consts::PI
            } else {
                16.0
            };
            let want = area * length;
            assert!(
                (v - want).abs() < want * 2e-3,
                "frenet {frenet} holed {holed}: {v} against {want}"
            );
        }
    }
}
