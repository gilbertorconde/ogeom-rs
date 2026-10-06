//! Which way a shape's faces walk the edges they share: in a closed shell
//! each face keeps its material on the left of its rings, so every edge
//! between two faces is walked once each way.

use std::collections::HashMap;

use ogeom::topo::{Filter, Location, Model, Orientation, Shape, ShapeType, TShapeId, explore};

/// How many edges between two faces of `shape` both faces walk the same
/// way: none where every face keeps its material on the left of its rings.
/// A degenerate edge (a pole, an apex) bounds nothing and is not counted.
pub fn edges_walked_one_way(model: &Model, shape: &Shape) -> usize {
    let mut walks: HashMap<(TShapeId, Location), (usize, usize)> = HashMap::new();
    for face in explore(model, shape, Filter::OfType(ShapeType::Face)).unwrap() {
        for edge in explore(model, &face, Filter::OfType(ShapeType::Edge)).unwrap() {
            if model
                .node(&edge)
                .and_then(|n| n.data().as_edge())
                .is_some_and(|d| d.degenerate)
            {
                continue;
            }
            let walk = walks
                .entry((edge.node(), edge.location().clone()))
                .or_default();
            walk.0 += 1;
            if edge.orientation() == Orientation::Forward {
                walk.1 += 1;
            }
        }
    }
    walks
        .values()
        .filter(|(uses, forward)| *uses == 2 && *forward != 1)
        .count()
}
