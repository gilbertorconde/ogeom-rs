//! `fix_shape`: one pass over a shape nobody promised was well-formed.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, make_box, make_compound, make_face, make_polygon, make_vertex, surface_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{PlaneSurface, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Plane, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{FaceData, Location, Model, NodeData, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A wire whose edges were handed over in the wrong order is put in the
/// order that walks them.
#[test]
fn a_wire_out_of_order_is_reordered() {
    let mut model = Model::new();
    let square = make_polygon(
        &mut model,
        &[
            Point::new(0.0, 0.0, 0.0),
            Point::new(4.0, 0.0, 0.0),
            Point::new(4.0, 4.0, 0.0),
            Point::new(0.0, 4.0, 0.0),
        ],
        true,
        T,
    )
    .unwrap()
    .shape;
    let mut edges = model.ordered_children_of(&square).unwrap();
    edges.swap(1, 3);
    // The model takes a bag; only the checker knows it is no path. The
    // face is built straight on it, since the builder would refuse.
    let shuffled = model.add_wire(&edges).unwrap();
    let plane = SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(Frame::WORLD)));
    let surface = model.geometry_mut().add_surface(plane);
    let face = model
        .add_face(FaceData::new(surface, Location::default()), &[shuffled])
        .unwrap();
    let before = check(&model, &face, T).unwrap();
    assert!(!before.is_usable(), "the checker sees the gap: {before}");

    let fixed = ogeom::heal::fix_shape(&mut model, &face, T).unwrap();
    assert_eq!(fixed.report.wires_reordered, 1);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let area = surface_properties(&model, &fixed.shape, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((area - 16.0).abs() < 1e-9, "the same square: {area}");
}

/// An edge shorter than its vertices' tolerances collapses to one vertex.
#[test]
fn a_tiny_edge_is_collapsed() {
    let mut model = Model::new();
    // A square with a fifth corner a micron from the fourth.
    let ring = make_polygon(
        &mut model,
        &[
            Point::new(0.0, 0.0, 0.0),
            Point::new(4.0, 0.0, 0.0),
            Point::new(4.0, 4.0, 0.0),
            Point::new(0.0, 4.0, 0.0),
            Point::new(0.0, 4.0 - 1e-6, 0.0),
        ],
        true,
        T,
    )
    .unwrap()
    .shape;
    // The model admits the two corners are one place: a vertex tolerance
    // of ten microns, as an imported vertex might carry.
    let vertices = explore_unique(&model, &ring, ShapeType::Vertex).unwrap();
    let near = vertices
        .iter()
        .find(|v| {
            model
                .node(v)
                .and_then(|n| n.data().as_vertex())
                .is_some_and(|d| d.point.distance(Point::new(0.0, 4.0, 0.0)) < 1e-9)
        })
        .expect("the fourth corner");
    if let Some(node) = model.node_mut(near)
        && let NodeData::Vertex(data) = node.data_mut()
    {
        data.tolerance = data.tolerance.widen_to(1e-5);
    }
    let plane = SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(Frame::WORLD)));
    let face = make_face(&mut model, plane, &[ring], T).unwrap().shape;
    assert_eq!(
        explore_unique(&model, &face, ShapeType::Edge)
            .unwrap()
            .len(),
        5
    );

    let fixed = ogeom::heal::fix_shape(&mut model, &face, T).unwrap();
    assert_eq!(fixed.report.edges_collapsed, 1);
    assert_eq!(
        explore_unique(&model, &fixed.shape, ShapeType::Edge)
            .unwrap()
            .len(),
        4,
        "a square again"
    );
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let area = surface_properties(&model, &fixed.shape, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((area - 16.0).abs() < 1e-4, "the same square: {area}");
}

/// A compound of loose faces is sewn into the shell they make.
#[test]
fn loose_faces_are_sewn() {
    let mut model = Model::new();
    // Six faces of a unit cube, each with vertices and edges of its own:
    // nothing shared but position.
    let p = |x: f64, y: f64, z: f64| Point::new(x, y, z);
    let quads: [[Point; 4]; 6] = [
        [p(0., 0., 0.), p(0., 1., 0.), p(1., 1., 0.), p(1., 0., 0.)],
        [p(0., 0., 1.), p(1., 0., 1.), p(1., 1., 1.), p(0., 1., 1.)],
        [p(0., 0., 0.), p(1., 0., 0.), p(1., 0., 1.), p(0., 0., 1.)],
        [p(0., 1., 0.), p(0., 1., 1.), p(1., 1., 1.), p(1., 1., 0.)],
        [p(0., 0., 0.), p(0., 0., 1.), p(0., 1., 1.), p(0., 1., 0.)],
        [p(1., 0., 0.), p(1., 1., 0.), p(1., 1., 1.), p(1., 0., 1.)],
    ];
    let mut faces = Vec::new();
    for quad in &quads {
        let ring = make_polygon(&mut model, quad, true, T).unwrap().shape;
        let normal = Direction::new((quad[1] - quad[0]).cross(quad[2] - quad[0]), T).unwrap();
        let reference = Direction::new(quad[1] - quad[0], T).unwrap();
        let frame = Frame::new(quad[0], normal, reference, T).unwrap();
        let plane = SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(frame)));
        faces.push(make_face(&mut model, plane, &[ring], T).unwrap().shape);
    }
    let bag = make_compound(&mut model, &faces).unwrap().shape;

    let fixed = ogeom::heal::fix_shape(&mut model, &bag, T).unwrap();
    let (joined, free) = fixed
        .report
        .sewn
        .expect("sewing ran over a compound of faces");
    assert_eq!(free, 0, "the box closes");
    assert!(
        joined >= 12,
        "every edge of the box was found shared: {joined}"
    );
    let shells = explore_unique(&model, &fixed.shape, ShapeType::Shell).unwrap();
    assert_eq!(shells.len(), 1);
    assert!(fixed.report.after.is_usable(), "{}", fixed.report.after);
}

/// A shape with nothing wrong is handed back as it is.
#[test]
fn a_sound_shape_is_left_alone() {
    let mut model = Model::new();
    let solid = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T)
        .unwrap()
        .shape;
    let fixed = ogeom::heal::fix_shape(&mut model, &solid, T).unwrap();
    assert!(fixed.shape.is_same(&solid), "the same node");
    assert_eq!(fixed.report.wires_reordered, 0);
    assert_eq!(fixed.report.edges_collapsed, 0);
    assert_eq!(fixed.report.edges_trimmed, 0);
    assert!(fixed.report.sewn.is_none(), "a solid is not sewn");
    assert!(fixed.report.after.is_valid());
    let _ = make_vertex(&mut model, Point::ORIGIN);
}
