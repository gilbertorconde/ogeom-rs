//! A prism of a planar face standing far from the origin on a tilted plane
//! is the same solid as the one built at the origin and moved there: valid,
//! facing out, with the volume of its profile times its height.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom::algo::{
    check, make_edge, make_edge_between, make_face, make_polygon, make_prism, make_wire,
    volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve, Curve3d, LineCurve, PlaneSurface};
use ogeom::math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, VertexData};

const T: Tolerances = Tolerances::millimetres();

/// The profile: a 40 by 20 rectangle in the plane through `(o, -o, o)`
/// square to `(1, 1, 1)`, its top side optionally a half circle of radius
/// 20, optionally holed by a circle of radius 5. Returns the face and its
/// area.
fn profile(model: &mut Model, o: f64, arc: bool, hole: bool) -> (Shape, f64) {
    let n = Direction::new(Vector::new(1.0, 1.0, 1.0), T).unwrap();
    let fr = Frame::new(Point::new(o, -o, o), n, Direction::X, T).unwrap();
    let at = |x: f64, y: f64| fr.origin() + fr.x().vector() * x + fr.y().vector() * y;
    let mut area = 800.0;
    let outer = if arc {
        area += 200.0 * PI;
        let corners = [(0.0, 20.0), (0.0, 0.0), (40.0, 0.0), (40.0, 20.0)].map(|(x, y)| at(x, y));
        let vertices: Vec<Shape> = corners
            .iter()
            .map(|p| model.add_vertex(VertexData::new(*p)))
            .collect();
        let mut edges = Vec::new();
        for i in 0..3 {
            let line = LineCurve::segment(corners[i], corners[i + 1], T).unwrap();
            let range = line.domain();
            edges.push(
                make_edge_between(
                    model,
                    Curve::Line(line),
                    range,
                    &vertices[i],
                    &vertices[i + 1],
                    T,
                )
                .unwrap()
                .shape,
            );
        }
        let top = Circle::new(
            Frame::new(at(20.0, 20.0), fr.z(), fr.x(), T).unwrap(),
            20.0,
            T,
        )
        .unwrap();
        edges.push(
            make_edge_between(
                model,
                Curve::Circle(CircleCurve::new(top)),
                (0.0, PI),
                &vertices[3],
                &vertices[0],
                T,
            )
            .unwrap()
            .shape,
        );
        make_wire(model, &edges, T).unwrap().shape
    } else {
        let corners = [(0.0, 0.0), (40.0, 0.0), (40.0, 20.0), (0.0, 20.0)].map(|(x, y)| at(x, y));
        make_polygon(model, &corners, true, T).unwrap().shape
    };
    let mut wires = vec![outer];
    if hole {
        area -= 25.0 * PI;
        let circle = Circle::new(
            Frame::new(at(20.0, 15.0), fr.z(), fr.x(), T).unwrap(),
            5.0,
            T,
        )
        .unwrap();
        let edge = make_edge(
            model,
            Curve::Circle(CircleCurve::new(circle)),
            (0.0, 2.0 * PI),
            T,
        )
        .unwrap()
        .shape;
        wires.push(make_wire(model, &[edge], T).unwrap().shape.reversed());
    }
    let face = make_face(model, PlaneSurface::new(Plane::new(fr)).into(), &wires, T)
        .unwrap()
        .shape;
    (face, area)
}

#[test]
fn a_far_tilted_prism_is_valid_and_measures_its_profile_times_its_height() {
    let n = Vector::new(1.0, 1.0, 1.0) / 3.0_f64.sqrt();
    let mut failures = Vec::new();
    for o in [0.0, 1e6, 5e6, 7e6, 1e7, 2e7] {
        for (arc, hole) in [(false, false), (true, false), (false, true), (true, true)] {
            let mut model = Model::new();
            let (face, area) = profile(&mut model, o, arc, hole);
            let tag = format!("o {o:e}, arc {arc}, hole {hole}");
            let prism = match make_prism(&mut model, &face, n * 30.0, T) {
                Ok(built) => built.shape,
                Err(e) => {
                    failures.push(format!("{tag}: make_prism: {e}"));
                    continue;
                }
            };
            let report = check(&model, &prism, T).unwrap();
            if !report.is_valid() {
                failures.push(format!("{tag}: check: {report:?}"));
                continue;
            }
            match volume_properties(&model, &prism, Deflection::default(), T) {
                Ok(props) => {
                    let want = 30.0 * area;
                    let error = (props.mass - want).abs() / want;
                    if error > 1e-6 {
                        failures.push(format!("{tag}: volume {} not {want}", props.mass));
                    }
                }
                Err(e) => failures.push(format!("{tag}: volume: {e}")),
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}
