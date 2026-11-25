//! A closed wire reached reversed, such as the boundary of a face used
//! against its stored sense, is still closed, and a loft takes it as a
//! section whichever way it is reached.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Frame, Point};
use ogeom_topo::{Filter, Model, Shape, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

fn point_of(model: &Model, vertex: &Shape) -> Point {
    let data = model
        .node(vertex)
        .and_then(|n| n.data().as_vertex())
        .unwrap();
    vertex.transform(model.datums()).unwrap().apply(data.point)
}

fn z_of(model: &Model, vertex: &Shape) -> f64 {
    point_of(model, vertex).z
}

/// The section's corners in its walk, from its start, drawn half size
/// about (5, 10) and 10 below: a loft takes each section in its own sense
/// from its own start, so the far section follows the near one's walk.
fn shrunk_below(model: &mut Model, section: &Shape) -> Shape {
    let corners: Vec<Point> = model
        .ordered_children_of(section)
        .unwrap()
        .iter()
        .map(|edge| {
            let (start, _) = ogeom_algo::edge_vertices(model, edge).unwrap().unwrap();
            let p = point_of(model, &start);
            Point::new(5.0 + (p.x - 5.0) / 2.0, 10.0 + (p.y - 10.0) / 2.0, -10.0)
        })
        .collect();
    ogeom_algo::make_polygon(model, &corners, true, T)
        .unwrap()
        .shape
}

/// The wire of a 10 x 20 x 5 box's face at z = 0, as the shell reaches it.
/// Its edges are stored in one order; reversed, the walk runs the other way.
fn bottom_wire(model: &mut Model) -> Shape {
    let block = ogeom_algo::make_box(model, Frame::WORLD, (10.0, 20.0, 5.0), T)
        .unwrap()
        .shape;
    let bottom = explore(model, &block, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find(|face| {
            explore(model, face, Filter::OfType(ShapeType::Vertex))
                .unwrap()
                .iter()
                .all(|v| z_of(model, v).abs() < 1e-9)
        })
        .unwrap();
    model
        .children_of(&bottom)
        .unwrap()
        .into_iter()
        .find(|c| model.kind_of(c) == Ok(ShapeType::Wire))
        .unwrap()
}

#[test]
fn a_face_s_wire_is_closed_either_way() {
    let mut model = Model::new();
    let wire = bottom_wire(&mut model);
    for wire in [wire.clone(), wire.reversed()] {
        assert!(
            ogeom_algo::is_wire_closed(&model, &wire, T).unwrap(),
            "{:?}",
            wire.orientation()
        );
    }
}

#[test]
fn a_loft_from_a_face_s_wire_either_way_is_the_frustum() {
    let mut model = Model::new();
    let wire = bottom_wire(&mut model);
    for (section, ruled) in [
        (wire.clone(), true),
        (wire.clone(), false),
        (wire.reversed(), true),
        (wire.reversed(), false),
    ] {
        let sense = section.orientation();
        let rectangle = shrunk_below(&mut model, &section);
        let result =
            ogeom_offset::make_loft_solid(&mut model, &[section, rectangle], &[], ruled, T)
                .unwrap_or_else(|e| panic!("{sense:?} ruled {ruled}: {e}"));
        let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "{sense:?} ruled {ruled}: {:?}",
            diagnosis.problems
        );
        // The pyramidal frustum: h/3 (A1 + A2 + sqrt(A1 A2)).
        let expected = 10.0 / 3.0 * (200.0 + 50.0 + 100.0);
        let measured = ogeom_algo::volume_properties(
            &model,
            &result.shape,
            ogeom_mesh::Deflection::default(),
            T,
        )
        .unwrap()
        .mass;
        assert!(
            (measured - expected).abs() < 1e-6 * expected,
            "{sense:?} ruled {ruled}: loft volume {measured} against {expected}"
        );
    }
}
