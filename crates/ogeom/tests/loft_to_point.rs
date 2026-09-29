//! Lofts and pipes closing to a point through several sections.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{make_edge, make_polygon, make_wire, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::CircleCurve;
use ogeom::math::{Circle, Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, VertexData, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

fn rectangle(model: &mut Model, (x0, y0): (f64, f64), (x1, y1): (f64, f64), z: f64) -> Shape {
    let corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)].map(|(x, y)| Point::new(x, y, z));
    make_polygon(model, &corners, true, T).unwrap().shape
}

fn circle(model: &mut Model, radius: f64, z: f64) -> Shape {
    let frame = Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X, T).unwrap();
    let curve = CircleCurve::new(Circle::new(frame, radius, T).unwrap());
    let edge = make_edge(model, curve.into(), (0.0, 2.0 * PI), T)
        .unwrap()
        .shape;
    make_wire(model, &[edge], T).unwrap().shape
}

fn has_vertex_at(model: &Model, shape: &Shape, at: Point) -> bool {
    explore_unique(model, shape, ShapeType::Vertex)
        .unwrap()
        .iter()
        .any(|v| {
            let p = model.node(v).unwrap().data().as_vertex().unwrap().point;
            p.distance(at) <= 1e-9
        })
}

/// Two rectangles and a point: the skin passes through the middle section
/// and closes on the point. The middle section is one of the skin's own
/// rows, so it is probed by classifying points either side of it.
#[test]
fn a_skinned_loft_of_rectangles_closes_to_a_point() {
    let mut model = Model::new();
    let base = rectangle(&mut model, (0.0, 0.0), (10.0, 20.0), 0.0);
    let middle = rectangle(&mut model, (1.0, 2.0), (9.0, 18.0), 5.0);
    let apex = model.add_vertex(VertexData::new(Point::new(5.0, 10.0, 10.0)));
    let loft = ogeom::offset::make_loft_skinned(&mut model, &[base, middle, apex], 1e-3, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &loft, T).unwrap().is_valid());
    assert!(has_vertex_at(&model, &loft, Point::new(5.0, 10.0, 10.0)));
    // The skin passes through the middle rectangle: just inside its corners
    // and sides is material, just outside is not.
    let fine = Deflection::with_chord(1e-3).unwrap();
    let side = |at: Point| ogeom::algo::classify_in_solid(&model, &loft, at, fine, T).unwrap();
    for (x, y) in [(1.0_f64, 2.0_f64), (9.0, 18.0), (5.0, 2.0), (9.0, 10.0)] {
        let inward = Point::new(
            x + (5.0 - x).signum() * 0.02,
            y + (10.0 - y).signum() * 0.02,
            5.0,
        );
        let outward = Point::new(
            x - (5.0 - x).signum() * 0.02,
            y - (10.0 - y).signum() * 0.02,
            5.0,
        );
        assert_eq!(side(inward), ogeom::algo::Containment::In, "{inward:?}");
        assert_eq!(side(outward), ogeom::algo::Containment::Out, "{outward:?}");
    }
}

/// Two coaxial circles and a point on their axis: a solid of revolution
/// about the axis, closing exactly on the point. Its volume is the one its
/// own meridian (read off the seam, which runs along y = 0) sweeps.
#[test]
fn a_skinned_loft_of_circles_closes_to_a_point_on_the_axis() {
    let mut model = Model::new();
    let base = circle(&mut model, 5.0, 0.0);
    let middle = circle(&mut model, 3.0, 5.0);
    let apex = model.add_vertex(VertexData::new(Point::new(0.0, 0.0, 10.0)));
    let loft = ogeom::offset::make_loft_skinned(&mut model, &[base, middle, apex], 1e-3, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &loft, T).unwrap().is_valid());
    assert!(has_vertex_at(&model, &loft, Point::new(0.0, 0.0, 10.0)));
    let mesh =
        ogeom::mesh::triangulate(&model, &loft, Deflection::with_chord(1e-3).unwrap(), T).unwrap();
    let mut meridian: Vec<(f64, f64)> = mesh
        .positions
        .iter()
        .filter(|p| p.y.abs() < 1e-9 && p.x > 0.0 && p.z > 1e-9)
        .map(|p| (p.z, p.x))
        .collect();
    meridian.sort_by(|a, b| a.0.total_cmp(&b.0));
    meridian.insert(0, (0.0, 5.0));
    meridian.push((10.0, 0.0));
    let swept: f64 = meridian
        .windows(2)
        .map(|w| PI * (w[1].0 - w[0].0) * (w[0].1 * w[0].1 + w[1].1 * w[1].1) / 2.0)
        .sum();
    let v = volume(&model, &loft);
    assert!((v - swept).abs() < swept * 1e-3, "{v} against {swept}");
}
