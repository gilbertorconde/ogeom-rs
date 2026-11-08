//! A corner of three bevels closed as a point or a facet.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{Chamfer, ChamferCorner, chamfer_edges_cornered, chamfer_edges_with};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn vol(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

fn faces(m: &Model, s: &Shape) -> usize {
    explore_unique(m, s, ShapeType::Face).unwrap().len()
}

/// A 20 × 20 × 10 box and its three edges at the origin, each with `spec`.
fn corner(
    m: &mut Model,
    spec: impl Fn(usize, &Model, &Shape) -> Chamfer,
) -> (Shape, Vec<(Shape, Chamfer)>) {
    let frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
    let b = ogeom_algo::make_box(m, frame, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let at_origin: Vec<Shape> = explore_unique(m, &b, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (s, t)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            [s, t]
                .iter()
                .any(|&u| c.point_at(u, T).unwrap().distance(Point::ORIGIN) < 1e-9)
        })
        .collect();
    assert_eq!(at_origin.len(), 3);
    let specs = at_origin
        .iter()
        .enumerate()
        .map(|(n, e)| (e.clone(), spec(n, m, &b)))
        .collect();
    (b, specs)
}

#[test]
fn a_point_corner_is_the_mitre() {
    let mut m = Model::new();
    let (b, specs) = corner(&mut m, |_, _, _| Chamfer::Symmetric(2.0));
    let mut plain = m.clone();
    let a = chamfer_edges_with(&mut plain, &b, &specs, T).unwrap();
    let p = chamfer_edges_cornered(&mut m, &b, &specs, ChamferCorner::Point, T).unwrap();
    assert_eq!(faces(&plain, &a.shape), faces(&m, &p.shape));
    assert_eq!(vol(&plain, &a.shape).to_bits(), vol(&m, &p.shape).to_bits());
}

/// At the box's corner the bevels `x + y = 2`, `y + z = 2` and `x + z = 2`
/// meet at `(1, 1, 1)`, and their contact lines cross on the faces at
/// `(2, 2, 0)`, `(0, 2, 2)` and `(2, 0, 2)`: the facet `x + y + z = 4`
/// takes the tetrahedron between them off the mitre.
#[test]
fn a_facet_corner_takes_the_mitre_s_tetrahedron_off() {
    let mut m = Model::new();
    let (b, specs) = corner(&mut m, |_, _, _| Chamfer::Symmetric(2.0));
    let mut plain = m.clone();
    let point = chamfer_edges_with(&mut plain, &b, &specs, T).unwrap();
    let facet = chamfer_edges_cornered(&mut m, &b, &specs, ChamferCorner::Facet, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &facet.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert_eq!(faces(&m, &facet.shape), faces(&plain, &point.shape) + 1);
    let apex = Point::new(1.0, 1.0, 1.0);
    let [a, b, c] = [
        Point::new(2.0, 2.0, 0.0),
        Point::new(0.0, 2.0, 2.0),
        Point::new(2.0, 0.0, 2.0),
    ]
    .map(|q| q - apex);
    let tetrahedron = a.dot(b.cross(c)).abs() / 6.0;
    let expected = vol(&plain, &point.shape) - tetrahedron;
    let v = vol(&m, &facet.shape);
    assert!(
        (v - expected).abs() < 1e-9 * expected,
        "{v} against {expected}"
    );
    // The facet: the one face within the corner's 2 mm cube.
    let on_facet = explore_unique(&m, &facet.shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            let bounds = ogeom_algo::shape_bounds(&m, f, T).unwrap();
            let hi = bounds.high().unwrap();
            hi.x <= 2.0 + 1e-3 && hi.y <= 2.0 + 1e-3 && hi.z <= 2.0 + 1e-3
        })
        .count();
    assert_eq!(on_facet, 1);
}

#[test]
fn a_facet_closes_unequal_bevels() {
    // Three bevels of 1, 1.5 and 2: the facet still stands on the three
    // crossings, each bevel ending on it.
    let mut m = Model::new();
    let (b, specs) = corner(&mut m, |n, _, _| Chamfer::Symmetric([1.0, 1.5, 2.0][n]));
    let mut plain = m.clone();
    let point = chamfer_edges_with(&mut plain, &b, &specs, T).unwrap();
    let facet = chamfer_edges_cornered(&mut m, &b, &specs, ChamferCorner::Facet, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &facet.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert_eq!(faces(&m, &facet.shape), faces(&plain, &point.shape) + 1);
    assert!(vol(&m, &facet.shape) < vol(&plain, &point.shape));
}
