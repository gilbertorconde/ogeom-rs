//! An operation run on a boolean's result leaves the boolean's operands as
//! they were, though the result holds the faces the boolean set aside as
//! the operand's own nodes.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use std::collections::HashSet;

use ogeom::algo::{check, make_box, make_cylinder};
use ogeom::core::{Tolerance, Tolerances};
use ogeom::math::{Direction, Frame, Point};
use ogeom::topo::{EdgeRepr, Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A block 20 by 20 by 10 with its far upright edge (at x = y = 20) and
/// that edge's vertices stating 1e-3, and the block with a pin drilled
/// through it near the opposite corner: the far faces are set aside.
fn drilled() -> (Model, Shape, Shape) {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let far = edge_where(&model, &block, |p| p.x == 20.0 && p.y == 20.0);
    model.widen(&far, Tolerance::new(1e-3).unwrap()).unwrap();
    let frame = Frame::new(Point::new(3.0, 3.0, -5.0), Direction::Z, Direction::X, T).unwrap();
    let pin = make_cylinder(&mut model, frame, 1.5, 20.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom::boolean::cut(&mut model, &block, &pin, T)
        .unwrap()
        .shape;
    let block_faces = nodes_of(&model, &block, ShapeType::Face);
    let shared = nodes_of(&model, &drilled, ShapeType::Face)
        .intersection(&block_faces)
        .count();
    assert!(shared >= 2, "the far faces come through as the block's own");
    (model, block, drilled)
}

fn nodes_of(model: &Model, shape: &Shape, kind: ShapeType) -> HashSet<ogeom::topo::TShapeId> {
    explore_unique(model, shape, kind)
        .unwrap()
        .iter()
        .map(Shape::node)
        .collect()
}

fn points_of(model: &Model, shape: &Shape) -> Vec<Point> {
    explore_unique(model, shape, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| {
            let at = model.node(v).unwrap().data().as_vertex().unwrap().point;
            v.transform(model.datums()).unwrap().apply(at)
        })
        .collect()
}

/// The one edge of `shape` whose ends both pass `keep`.
fn edge_where(model: &Model, shape: &Shape, keep: impl Fn(Point) -> bool) -> Shape {
    let found: Vec<Shape> = explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| points_of(model, e).into_iter().all(&keep))
        .collect();
    assert_eq!(found.len(), 1);
    found[0].clone()
}

/// The one face of `shape` whose vertices all pass `keep`.
fn face_where(model: &Model, shape: &Shape, keep: impl Fn(Point) -> bool) -> Shape {
    let found: Vec<Shape> = explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| points_of(model, f).into_iter().all(&keep))
        .collect();
    assert_eq!(found.len(), 1);
    found[0].clone()
}

/// Every node below `root` as it stands: its data and children, and the
/// curves, pcurves and surfaces it names.
fn snapshot(model: &Model, root: &Shape) -> Vec<String> {
    let mut out = Vec::new();
    let mut seen = HashSet::new();
    let mut stack = vec![root.node()];
    while let Some(id) = stack.pop() {
        if !seen.insert(id) {
            continue;
        }
        let node = model.node_by_id(id).unwrap();
        out.push(format!("{id:?} {node:?}"));
        let geometry = model.geometry();
        match node.data() {
            NodeData::Edge(edge) => {
                for repr in &edge.representations {
                    match repr {
                        EdgeRepr::Curve3d { curve, .. } => {
                            out.push(format!("{id:?} {:?}", geometry.curve(*curve)));
                        }
                        EdgeRepr::PCurve { curve, .. } => {
                            out.push(format!("{id:?} {:?}", geometry.pcurve(*curve)));
                        }
                        EdgeRepr::Seam {
                            forward, reversed, ..
                        } => {
                            out.push(format!(
                                "{id:?} {:?} {:?}",
                                geometry.pcurve(*forward),
                                geometry.pcurve(*reversed)
                            ));
                        }
                        _ => {}
                    }
                }
            }
            NodeData::Face(face) => {
                out.push(format!("{id:?} {:?}", geometry.surface(face.surface)));
            }
            _ => {}
        }
        stack.extend(node.children().iter().map(Shape::node));
    }
    out.sort();
    out
}

/// `edit` on the drilled block leaves the block, and whatever else is
/// named, exactly as it was and valid.
fn leaves_the_block(edit: impl FnOnce(&mut Model, &Shape, &Shape) -> Vec<Shape>) {
    let (mut model, block, drilled) = drilled();
    let before = snapshot(&model, &block);
    let named = edit(&mut model, &block, &drilled);
    assert_eq!(snapshot(&model, &block), before, "the block changed");
    assert!(check(&model, &block, T).unwrap().is_valid());
    for shape in named {
        assert!(check(&model, &shape, T).unwrap().is_usable());
    }
}

fn all_same_parameter(model: &Model, shape: &Shape) -> bool {
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .iter()
        .all(|e| {
            model
                .node(e)
                .unwrap()
                .data()
                .as_edge()
                .unwrap()
                .same_parameter()
        })
}

#[test]
fn repairing_same_parameter_on_the_result_leaves_the_operand() {
    leaves_the_block(|model, block, drilled| {
        assert!(!all_same_parameter(model, block));
        ogeom::heal::repair_same_parameter(model, drilled, T).unwrap();
        assert!(all_same_parameter(model, drilled));
        vec![drilled.clone()]
    });
}

#[test]
fn reducing_tolerances_on_the_result_leaves_the_operand() {
    leaves_the_block(|model, _, drilled| {
        let shrunk = ogeom::heal::reduce_tolerances(model, drilled, T).unwrap();
        assert!(shrunk >= 1, "the far edge's 1e-3 comes down");
        let far = edge_where(model, drilled, |p| p.x == 20.0 && p.y == 20.0);
        let edge = model.node(&far).unwrap().data().as_edge().unwrap();
        assert!(edge.tolerance.get() < 1e-3);
        vec![drilled.clone()]
    });
}

#[test]
fn fixing_the_result_leaves_the_operand() {
    leaves_the_block(|model, _, drilled| {
        let fixed = ogeom::heal::fix_shape(model, drilled, T).unwrap();
        assert!(fixed.report.tolerances_reduced >= 1);
        vec![fixed.shape]
    });
}

#[test]
fn sewing_a_face_of_the_result_leaves_the_operand() {
    leaves_the_block(|model, _, drilled| {
        // The result's far side at x = 20 and the top of a block standing
        // against it, which meet along the side's top edge.
        let side = face_where(model, drilled, |p| p.x == 20.0);
        let frame = Frame::new(Point::new(20.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
        let beside = make_box(model, frame, (10.0, 20.0, 10.0), T).unwrap().shape;
        let top = face_where(model, &beside, |p| p.z == 10.0);
        let sewn = ogeom::algo::sew(model, &[side, top], T).unwrap();
        assert_eq!(sewn.joined, 1);
        sewn.shells
    });
}

#[test]
fn filleting_the_result_leaves_the_operand_and_the_result() {
    leaves_the_block(|model, _, drilled| {
        let far = edge_where(model, drilled, |p| p.x == 20.0 && p.y == 20.0);
        let before = snapshot(model, drilled);
        let rounded = ogeom::fillet::fillet_edges(model, drilled, &[far], 2.0, T).unwrap();
        assert_eq!(
            snapshot(model, drilled),
            before,
            "the drilled block changed"
        );
        vec![rounded.shape]
    });
}

#[test]
fn tessellating_the_result_leaves_the_operand() {
    leaves_the_block(|model, _, drilled| {
        let meshed =
            ogeom::mesh::tessellate(model, drilled, ogeom::mesh::Deflection::default(), T).unwrap();
        assert!(meshed.faces > 0);
        for face in explore_unique(model, drilled, ShapeType::Face).unwrap() {
            assert!(ogeom::mesh::triangulation_of(model, &face).is_some());
        }
        vec![drilled.clone()]
    });
}

#[test]
fn filleting_an_edge_beside_a_rebuilt_face_leaves_the_operand_and_the_result() {
    leaves_the_block(|model, _, drilled| {
        // Between the far side, set aside, and the top, which the pin
        // drills through.
        let far = edge_where(model, drilled, |p| p.x == 20.0 && p.z == 10.0);
        let before = snapshot(model, drilled);
        let rounded = ogeom::fillet::fillet_edges(model, drilled, &[far], 2.0, T).unwrap();
        assert_eq!(
            snapshot(model, drilled),
            before,
            "the drilled block changed"
        );
        vec![rounded.shape]
    });
}

#[test]
fn chamfering_and_rounding_a_corner_of_the_result_leave_the_operand_and_the_result() {
    leaves_the_block(|model, _, drilled| {
        let before = snapshot(model, drilled);
        let far = edge_where(model, drilled, |p| p.x == 20.0 && p.y == 20.0);
        let bevelled = ogeom::fillet::chamfer_edge(model, drilled, &far, 1.0, T).unwrap();
        let corner: Vec<Shape> = [
            |p: Point| p.x == 20.0 && p.y == 20.0,
            |p: Point| p.x == 20.0 && p.z == 10.0,
            |p: Point| p.y == 20.0 && p.z == 10.0,
        ]
        .into_iter()
        .map(|keep| edge_where(model, drilled, keep))
        .collect();
        let rounded = ogeom::fillet::fillet_edges(model, drilled, &corner, 1.0, T).unwrap();
        assert_eq!(
            snapshot(model, drilled),
            before,
            "the drilled block changed"
        );
        vec![bevelled.shape, rounded.shape]
    });
}
