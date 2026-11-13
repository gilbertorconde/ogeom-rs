//! Corner history identifies the new patch by the source vertex.
#![allow(clippy::unwrap_used, reason = "test constructions")]
use ogeom_core::Tolerances;
use ogeom_geom::SurfaceGeometry;
use ogeom_math::{Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};
const T: Tolerances = Tolerances::millimetres();

fn assert_patch(model: &Model, built: &ogeom_algo::Built, vertex: &Shape, centre: Point) {
    let faces = explore_unique(model, &built.shape, ShapeType::Face).unwrap();
    let patches: Vec<_> = built
        .history
        .generated(vertex)
        .iter()
        .filter(|image| faces.iter().any(|face| face.is_same(image)))
        .collect();
    let mut spheres = 0;
    for patch in patches {
        let surface = model.node(patch).unwrap().data().as_face().unwrap().surface;
        if let Some(SurfaceGeometry::Sphere(sphere)) = model.geometry().surface(surface) {
            spheres += 1;
            // A box's unit corner ball is exactly one unit from each wall.
            assert!(sphere.sphere().centre().distance(centre) < T.confusion());
            assert!((sphere.sphere().radius() - 1.0).abs() < T.confusion());
        }
    }
    assert_eq!(
        spheres, 1,
        "the source vertex must generate its own corner sphere"
    );
}

#[test]
fn a_vertex_generates_its_corner_faces_through_composed_history() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 10.0, 5.0), T)
        .unwrap()
        .shape;
    let vertices = explore_unique(&model, &solid, ShapeType::Vertex).unwrap();
    for vertex in &vertices {
        let point = model
            .node(vertex)
            .unwrap()
            .data()
            .as_vertex()
            .unwrap()
            .point;
        let centre = Point::new(
            if point.x < 10.0 { 1.0 } else { 19.0 },
            if point.y < 5.0 { 1.0 } else { 9.0 },
            if point.z < 2.5 { 1.0 } else { 4.0 },
        );
        let corner = ogeom_fillet::round_vertex(&mut model, &solid, vertex, 1.0, T).unwrap();
        assert_patch(&model, &corner, vertex, centre);
        let edges: Vec<_> = explore_unique(&model, &solid, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|edge| {
                ogeom_algo::edge_vertices(&model, edge)
                    .unwrap()
                    .is_some_and(|(a, b)| a.is_same(vertex) || b.is_same(vertex))
            })
            .collect();
        let rounded = ogeom_fillet::fillet_edges(&mut model, &solid, &edges, 1.0, T).unwrap();
        assert_patch(&model, &rounded, vertex, centre);
    }
    let edges = explore_unique(&model, &solid, ShapeType::Edge).unwrap();
    let rounded = ogeom_fillet::fillet_edges(&mut model, &solid, &edges, 1.0, T).unwrap();
    for vertex in &vertices {
        let point = model
            .node(vertex)
            .unwrap()
            .data()
            .as_vertex()
            .unwrap()
            .point;
        let centre = Point::new(
            if point.x < 10.0 { 1.0 } else { 19.0 },
            if point.y < 5.0 { 1.0 } else { 9.0 },
            if point.z < 2.5 { 1.0 } else { 4.0 },
        );
        assert_patch(&model, &rounded, vertex, centre);
    }
}
