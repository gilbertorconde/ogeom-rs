//! A swept solid faces outward whichever way its profile was walked.
//!
//! A closed wire on a plane bounds one region however it is walked, but a
//! prism reads each wall's side off its edge's direction. Walls read off
//! the edge direction alone face into the material for a square walked
//! clockwise about the travel while its caps face out: a closed shell the
//! checker accepts, whose volume is refused as wound inward. So every
//! ring's walls follow the ring's winding about the travel: the outer ring
//! turning positively, a hole the other way.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

#[path = "support/walks.rs"]
mod walks;

use ogeom::algo::{
    face_normal, make_edge_between, make_face_with_pcurves, make_prism, volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{Curve, Curve3d, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Filter, Model, Shape, ShapeType, VertexData, explore};

const T: Tolerances = Tolerances::millimetres();

/// A closed ring of segments through `corners` on the XY plane, in order.
fn ring(model: &mut Model, corners: &[(f64, f64)]) -> Vec<Shape> {
    let points: Vec<Point> = corners
        .iter()
        .map(|(x, y)| Point::new(*x, *y, 0.0))
        .collect();
    let vertices: Vec<Shape> = points
        .iter()
        .map(|p| model.add_vertex(VertexData::new(*p)))
        .collect();
    (0..points.len())
        .map(|i| {
            let j = (i + 1) % points.len();
            let curve = LineCurve::segment(points[i], points[j], T).unwrap();
            let range = curve.domain();
            make_edge_between(
                model,
                Curve::Line(curve),
                range,
                &vertices[i],
                &vertices[j],
                T,
            )
            .unwrap()
            .shape
        })
        .collect()
}

/// A planar face on XY from rings of corners, the first the outer one.
fn profile(model: &mut Model, rings: &[Vec<(f64, f64)>]) -> Shape {
    let wires: Vec<Vec<Shape>> = rings.iter().map(|r| ring(model, r)).collect();
    let plane = Plane::new(Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap());
    let surface = PlaneSurface::over(plane, (-1.0, 21.0), (-1.0, 21.0)).unwrap();
    make_face_with_pcurves(model, SurfaceGeometry::Plane(surface), &wires, T)
        .unwrap()
        .shape
}

fn square(clockwise: bool) -> Vec<(f64, f64)> {
    let mut corners = vec![(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0)];
    if clockwise {
        corners.reverse();
    }
    corners
}

/// The same box, whichever way its profile was walked.
#[test]
fn a_prism_from_a_clockwise_square_faces_out_everywhere() {
    for clockwise in [false, true] {
        let mut model = Model::new();
        let face = profile(&mut model, &[square(clockwise)]);
        let solid = make_prism(&mut model, &face, Vector::new(0.0, 0.0, 10.0), T)
            .unwrap()
            .shape;
        let faces = explore(&model, &solid, Filter::OfType(ShapeType::Face)).unwrap();
        assert_eq!(faces.len(), 6);
        let centre = Point::new(10.0, 10.0, 5.0);
        let outward = faces
            .iter()
            .filter(|f| {
                let (at, normal) = face_normal(&model, f, T).unwrap();
                normal.dot(at - centre) > 0.0
            })
            .count();
        assert_eq!(outward, 6, "every face presents away from the centre");
        let volume = volume_properties(&model, &solid, Deflection::default(), T)
            .unwrap()
            .mass;
        assert!((volume - 4000.0).abs() < 1e-3, "volume {volume}");
    }
}

/// A holed profile: either ring walked either way, swept either way along
/// the normal, is the same square tube.
#[test]
fn a_holed_prism_faces_out_whichever_way_each_ring_was_walked() {
    let hole = |clockwise: bool| {
        let mut corners = vec![(5.0, 5.0), (15.0, 5.0), (15.0, 15.0), (5.0, 15.0)];
        if clockwise {
            corners.reverse();
        }
        corners
    };
    for outer_cw in [false, true] {
        for hole_cw in [false, true] {
            for dz in [10.0, -10.0] {
                let mut model = Model::new();
                let face = profile(&mut model, &[square(outer_cw), hole(hole_cw)]);
                let solid = make_prism(&mut model, &face, Vector::new(0.0, 0.0, dz), T)
                    .unwrap()
                    .shape;
                assert!(ogeom::algo::check(&model, &solid, T).unwrap().is_valid());
                let volume = volume_properties(&model, &solid, Deflection::default(), T)
                    .unwrap()
                    .mass;
                assert!(
                    (volume - 3000.0).abs() < 1e-3,
                    "outer cw {outer_cw}, hole cw {hole_cw}, travel {dz}: volume {volume}"
                );
            }
        }
    }
}

/// A pipe along part of a circle is a torus segment whose two half-tube
/// patches and two end discs walk every edge once each way: valid, and
/// the volume Pappus names, wherever the arc starts.
#[test]
fn a_torus_segment_pipe_walks_each_edge_once_each_way() {
    use ogeom::geom::CircleCurve;
    use ogeom::math::Circle;
    let (major, minor) = (5.0, 1.0);
    for (from, to) in [(0.0, core::f64::consts::FRAC_PI_2), (1.0, 3.5)] {
        let mut model = Model::new();
        let frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
        let spine = ogeom::algo::make_edge(
            &mut model,
            CircleCurve::new(Circle::new(frame, major, T).unwrap()).into(),
            (from, to),
            T,
        )
        .unwrap()
        .shape;
        let solid = ogeom::offset::make_pipe(&mut model, &spine, minor, T)
            .unwrap()
            .shape;
        assert_eq!(
            walks::edges_walked_one_way(&model, &solid),
            0,
            "from {from}"
        );
        let found = ogeom::algo::check(&model, &solid, T).unwrap();
        assert!(found.is_valid(), "from {from}: {found}");
        let want = core::f64::consts::PI * minor * minor * major * (to - from);
        let fine = Deflection {
            chord: 1e-3,
            ..Deflection::default()
        };
        let v = volume_properties(&model, &solid, fine, T).unwrap().mass;
        assert!(
            (v - want).abs() < want * 1e-4,
            "from {from}: {v} against {want}"
        );
    }
}
