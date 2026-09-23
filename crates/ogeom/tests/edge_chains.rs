//! Edges treated as one operation: chamfers mitred along a chain, and three
//! or more fillets closing their shared corner with the rolling ball's
//! patch.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{make_edge_between, make_face_with_pcurves, make_prism, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::{Curve, Curve3d, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Filter, Model, Shape, ShapeType, VertexData, explore, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A 20 × 20 × 10 box swept from a square on XY: its top vertices and
/// edges are its bottom ones placed by the sweep's travel, which is what
/// a caller's prism hands over.
fn prism_box(model: &mut Model) -> Shape {
    let corners = [(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0)];
    let points: Vec<Point> = corners
        .iter()
        .map(|(x, y)| Point::new(*x, *y, 0.0))
        .collect();
    let vertices: Vec<Shape> = points
        .iter()
        .map(|p| model.add_vertex(VertexData::new(*p)))
        .collect();
    let edges: Vec<Shape> = (0..4)
        .map(|i| {
            let j = (i + 1) % 4;
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
        .collect();
    let plane = Plane::new(Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap());
    let surface = PlaneSurface::over(plane, (-1.0, 21.0), (-1.0, 21.0)).unwrap();
    let face = make_face_with_pcurves(model, SurfaceGeometry::Plane(surface), &[edges], T)
        .unwrap()
        .shape;
    make_prism(model, &face, Vector::new(0.0, 0.0, 10.0), T)
        .unwrap()
        .shape
}

/// The edge whose placed midpoint is nearest `at`.
fn edge_near(model: &Model, solid: &Shape, at: Point) -> Shape {
    let mid = |e: &Shape| -> Point {
        let data = model.node(e).unwrap().data().as_edge().unwrap().clone();
        let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            return Point::new(f64::INFINITY, 0.0, 0.0);
        };
        let local = model
            .geometry()
            .curve(*curve)
            .unwrap()
            .point_at(f64::midpoint(range.0, range.1), T)
            .unwrap();
        e.transform(model.datums()).unwrap().apply(local)
    };
    explore_unique(model, solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| mid(a).distance(at).total_cmp(&mid(b).distance(at)))
        .unwrap()
}

fn faces(model: &Model, shape: &Shape) -> usize {
    explore(model, shape, Filter::OfType(ShapeType::Face))
        .unwrap()
        .len()
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// The acceptance of the chamfer chain: the box's four top edges bevelled
/// by one millimetre as one operation mitre at every corner: ten faces,
/// the four prisms less their four corner overlaps.
#[test]
fn a_chamfer_chain_mitres_at_shared_vertices() {
    let mut model = Model::new();
    let solid = prism_box(&mut model);
    let top: Vec<Shape> = [
        Point::new(10.0, 0.0, 10.0),
        Point::new(20.0, 10.0, 10.0),
        Point::new(10.0, 20.0, 10.0),
        Point::new(0.0, 10.0, 10.0),
    ]
    .iter()
    .map(|at| edge_near(&model, &solid, *at))
    .collect();
    let bevelled = ogeom::fillet::chamfer_edges(&mut model, &solid, &top, 1.0, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &bevelled, T).unwrap().is_valid());
    assert_eq!(faces(&model, &bevelled), 10);
    let expected = 4000.0 - 4.0 * (0.5 * 20.0) + 4.0 / 3.0;
    let measured = volume(&model, &bevelled);
    assert!(
        (measured - expected).abs() < 1e-2,
        "volume {measured} against {expected}"
    );
}

/// Three bevels at a convex corner meet at one point, and every edge of
/// the box bevelled at once is the same inclusion–exclusion: the prisms,
/// less their pairwise overlaps at each corner, plus the triple overlap,
/// a quarter of the cube of the distance.
#[test]
fn chamfer_chains_meet_three_at_a_corner() {
    let mut model = Model::new();
    let solid = prism_box(&mut model);
    let corner: Vec<Shape> = [
        Point::new(10.0, 20.0, 10.0),
        Point::new(20.0, 10.0, 10.0),
        Point::new(20.0, 20.0, 5.0),
    ]
    .iter()
    .map(|at| edge_near(&model, &solid, *at))
    .collect();
    let bevelled = ogeom::fillet::chamfer_edges(&mut model, &solid, &corner, 1.0, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &bevelled, T).unwrap().is_valid());
    let expected = 4000.0 - (0.5 * (20.0 + 20.0 + 10.0) - 3.0 / 3.0 + 0.25);
    let measured = volume(&model, &bevelled);
    assert!(
        (measured - expected).abs() < 1e-2,
        "volume {measured} against {expected}"
    );

    let every = explore_unique(&model, &solid, ShapeType::Edge).unwrap();
    assert_eq!(every.len(), 12);
    let bevelled = ogeom::fillet::chamfer_edges(&mut model, &solid, &every, 1.0, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &bevelled, T).unwrap().is_valid());
    assert_eq!(faces(&model, &bevelled), 18, "six walls and twelve bevels");
    let expected = 4000.0 - (0.5 * (8.0 * 20.0 + 4.0 * 10.0) - 8.0 + 8.0 * 0.25);
    let measured = volume(&model, &bevelled);
    assert!(
        (measured - expected).abs() < 1e-2,
        "volume {measured} against {expected}"
    );
}

/// The acceptance of the fillet corner: three fillets at a box's top corner
/// asked together close it with the octant of the sphere a radius in from
/// all three faces: ten faces, and inside the corner cube nothing but that
/// sphere. The corner is the prism's far end, placed by the sweep, so the
/// corner tool must read where the vertex stands and not where its node
/// was built.
#[test]
fn three_fillets_at_a_corner_close_it_with_a_sphere() {
    let mut model = Model::new();
    let solid = prism_box(&mut model);
    let corner: Vec<Shape> = [
        Point::new(10.0, 20.0, 10.0),
        Point::new(20.0, 10.0, 10.0),
        Point::new(20.0, 20.0, 5.0),
    ]
    .iter()
    .map(|at| edge_near(&model, &solid, *at))
    .collect();
    let rounded = ogeom::fillet::fillet_edges(&mut model, &solid, &corner, 2.0, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &rounded, T).unwrap().is_valid());
    assert_eq!(faces(&model, &rounded), 10);

    let chord = 1e-2;
    let mesh = ogeom::mesh::triangulate(
        &model,
        &rounded,
        Deflection {
            chord,
            ..Deflection::default()
        },
        T,
    )
    .unwrap();
    let centre = Point::new(18.0, 18.0, 8.0);
    let inside =
        |p: Point| p.x > 18.0 && p.x < 20.0 && p.y > 18.0 && p.y < 20.0 && p.z > 8.0 && p.z < 10.0;
    let mut seen = 0;
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        let g = Point::new(
            (a.x + b.x + c.x) / 3.0,
            (a.y + b.y + c.y) / 3.0,
            (a.z + b.z + c.z) / 3.0,
        );
        if !inside(g) {
            continue;
        }
        seen += 1;
        for p in [a, b, c] {
            assert!(
                (p.distance(centre) - 2.0).abs() < 1e-6,
                "a vertex inside the corner off the sphere: {p:?}"
            );
        }
        assert!(
            (g.distance(centre) - 2.0).abs() < 2.0 * chord,
            "a triangle inside the corner off the sphere: {g:?}"
        );
    }
    assert!(seen > 0, "the corner is drawn");
}

/// At a vertex of more edges the corner still closes: a square pyramid's
/// four slopes filleted together round the apex with one sphere, and a
/// rectangular one's (no single ball touches its four slopes) with two
/// spheres and a cylinder between them.
#[test]
fn fillet_chains_close_an_apex() {
    for (half_y, spheres, cylinders) in [(10.0, 1, 4), (5.0, 2, 5)] {
        let mut model = Model::new();
        let base = [
            Point::new(-10.0, -half_y, 0.0),
            Point::new(10.0, -half_y, 0.0),
            Point::new(10.0, half_y, 0.0),
            Point::new(-10.0, half_y, 0.0),
        ];
        let apex = Point::new(0.0, 0.0, 15.0);
        let polygon = ogeom::algo::make_polygon(&mut model, &base, true, T)
            .unwrap()
            .shape;
        let tip = ogeom::algo::make_vertex(&mut model, apex).shape;
        let pyramid = ogeom::offset::make_loft(&mut model, &polygon, &tip, T)
            .unwrap()
            .shape;
        let slopes: Vec<Shape> = base
            .iter()
            .map(|b| edge_near(&model, &pyramid, *b + (apex - *b) * 0.5))
            .collect();
        let rounded = ogeom::fillet::fillet_edges(&mut model, &pyramid, &slopes, 1.5, T)
            .unwrap()
            .shape;
        assert!(ogeom::algo::check(&model, &rounded, T).unwrap().is_valid());
        let mut counts = (0, 0);
        for face in explore_unique(&model, &rounded, ShapeType::Face).unwrap() {
            let data = model.node(&face).unwrap().data().as_face().unwrap();
            match model.geometry().surface(data.surface) {
                Some(SurfaceGeometry::Sphere(_)) => counts.0 += 1,
                Some(SurfaceGeometry::Cylinder(_)) => counts.1 += 1,
                _ => {}
            }
        }
        assert_eq!(counts, (spheres, cylinders), "half width {half_y}");
    }
}
