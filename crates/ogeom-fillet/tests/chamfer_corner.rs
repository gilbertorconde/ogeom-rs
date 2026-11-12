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

/// The default corner is the facet (ADR 0001, §13): `chamfer_edges_with`
/// builds what [`ChamferCorner::Facet`] builds, to the bit. The mitre is
/// [`ChamferCorner::Point`], one face fewer: no facet, the three bevels
/// meeting at the point their planes share.
#[test]
fn the_default_corner_is_the_facet_and_the_point_corner_the_mitre() {
    let mut m = Model::new();
    let (b, specs) = corner(&mut m, |_, _, _| Chamfer::Symmetric(2.0));
    let mut plain = m.clone();
    let a = chamfer_edges_with(&mut plain, &b, &specs, T).unwrap();
    let mut faceted = m.clone();
    let f = chamfer_edges_cornered(&mut faceted, &b, &specs, ChamferCorner::Facet, T).unwrap();
    let p = chamfer_edges_cornered(&mut m, &b, &specs, ChamferCorner::Point, T).unwrap();
    assert_eq!(faces(&plain, &a.shape), faces(&faceted, &f.shape));
    assert_eq!(
        vol(&plain, &a.shape).to_bits(),
        vol(&faceted, &f.shape).to_bits()
    );
    assert_eq!(faces(&m, &p.shape) + 1, faces(&plain, &a.shape));
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
    let point = chamfer_edges_cornered(&mut plain, &b, &specs, ChamferCorner::Point, T).unwrap();
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
    let point = chamfer_edges_cornered(&mut plain, &b, &specs, ChamferCorner::Point, T).unwrap();
    let facet = chamfer_edges_cornered(&mut m, &b, &specs, ChamferCorner::Facet, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &facet.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert_eq!(faces(&m, &facet.shape), faces(&plain, &point.shape) + 1);
    assert!(vol(&m, &facet.shape) < vol(&plain, &point.shape));
}

/// A square pyramid, base 20 × 20 and apex 12 over its middle, and its four
/// edges at the apex, each with `Symmetric(d)`.
fn pyramid(m: &mut Model, d: f64) -> (Shape, Vec<(Shape, Chamfer)>) {
    let points = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(20.0, 0.0, 0.0),
        Point::new(20.0, 20.0, 0.0),
        Point::new(0.0, 20.0, 0.0),
        Point::new(10.0, 10.0, 12.0),
    ];
    let rings = vec![
        vec![0, 3, 2, 1],
        vec![0, 1, 4],
        vec![1, 2, 4],
        vec![2, 3, 4],
        vec![3, 0, 4],
    ];
    let p = ogeom_algo::make_polyhedron(m, &points, &rings, T)
        .unwrap()
        .shape;
    let apex = points[4];
    let specs = explore_unique(m, &p, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (s, t)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            [s, t]
                .iter()
                .any(|&u| c.point_at(u, T).unwrap().distance(apex) < 1e-9)
        })
        .map(|e| (e, Chamfer::Symmetric(d)))
        .collect();
    (p, specs)
}

#[test]
fn a_facet_closes_a_pyramid_s_apex_of_four_bevels_in_a_square() {
    // On each side face, between apex edges `u` and `w` at angle `phi`, the
    // two contact lines `d` from those edges cross on the bisector, `d /
    // sin(phi / 2)` from the apex: the facet's four corners, a square. The
    // mitred bevels meet on the axis, where bevel `i`'s plane, through its
    // two contact lines, crosses it. The facet takes off the square pyramid
    // on the facet with that point its apex.
    use ogeom_math::Vector;
    let d = 2.0;
    let mut m = Model::new();
    let (p, specs) = pyramid(&mut m, d);
    assert_eq!(specs.len(), 4);
    let mut plain = m.clone();
    let mitre = chamfer_edges_cornered(&mut plain, &p, &specs, ChamferCorner::Point, T).unwrap();
    let facet = chamfer_edges_cornered(&mut m, &p, &specs, ChamferCorner::Facet, T).unwrap();
    assert!(ogeom_algo::check(&m, &facet.shape, T).unwrap().is_usable());
    assert_eq!(faces(&m, &facet.shape), faces(&plain, &mitre.shape) + 1);

    let apex = Point::new(10.0, 10.0, 12.0);
    let base = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(20.0, 0.0, 0.0),
        Point::new(20.0, 20.0, 0.0),
        Point::new(0.0, 20.0, 0.0),
    ];
    let unit = |v: Vector| v / v.magnitude();
    let edges: Vec<Vector> = base.iter().map(|b| unit(*b - apex)).collect();
    let feet: Vec<Point> = (0..4)
        .map(|i| {
            let (u, w) = (edges[i], edges[(i + 1) % 4]);
            let phi = u.dot(w).acos();
            apex + unit(u + w) * (d / (phi / 2.0).sin())
        })
        .collect();
    // Bevel 0's plane: through the lines `d` off edge 0 in its two faces.
    let into = |u: Vector, w: Vector| unit(w - u * u.dot(w));
    let (u, w1, w2) = (edges[0], into(edges[0], edges[1]), into(edges[0], edges[3]));
    let normal = u.cross(w2 - w1);
    let down = Vector::new(0.0, 0.0, -1.0);
    let t = normal.dot(w1 * d) / normal.dot(down);
    let point = apex + down * t;
    let side = (feet[1] - feet[0]).magnitude();
    let height = point.z - feet[0].z;
    let piece = side * side * height / 3.0;
    let expected = vol(&plain, &mitre.shape) - piece;
    let got = vol(&m, &facet.shape);
    assert!(
        (got - expected).abs() <= expected * 1e-9,
        "{got} against {expected}"
    );
}

#[test]
fn a_facet_through_four_points_off_one_plane_is_refused() {
    // The apex moved off the middle: the four crossings no longer share a
    // plane, and no facet passes through them all: refused by name.
    let mut m = Model::new();
    let points = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(20.0, 0.0, 0.0),
        Point::new(20.0, 20.0, 0.0),
        Point::new(0.0, 20.0, 0.0),
        Point::new(5.0, 8.0, 12.0),
    ];
    let rings = vec![
        vec![0, 3, 2, 1],
        vec![0, 1, 4],
        vec![1, 2, 4],
        vec![2, 3, 4],
        vec![3, 0, 4],
    ];
    let p = ogeom_algo::make_polyhedron(&mut m, &points, &rings, T)
        .unwrap()
        .shape;
    let specs: Vec<(Shape, Chamfer)> = explore_unique(&m, &p, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (s, t)) = ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap();
            [s, t]
                .iter()
                .any(|&u| c.point_at(u, T).unwrap().distance(points[4]) < 1e-9)
        })
        .map(|e| (e, Chamfer::Symmetric(2.0)))
        .collect();
    let refused = chamfer_edges_cornered(&mut m, &p, &specs, ChamferCorner::Facet, T)
        .expect_err("the facet should be refused")
        .to_string();
    assert!(refused.contains("off one plane"), "{refused}");
}
