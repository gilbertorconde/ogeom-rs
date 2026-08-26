//! The reshape framework: record substitutions, apply them over a shape in
//! one pass.
//!
//! Healing operations want to say "this edge becomes that one, this face
//! goes away" and have the change ripple upward — every wire holding the
//! edge rebuilt, every face holding the wire, up to the solid — without
//! each fix reimplementing the traversal. A [`Reshape`] collects the
//! requests; [`Reshape::apply`] rebuilds bottom-up, sharing rebuilt nodes
//! so a substituted edge is one new node however many faces reach it, and
//! reports what became of every input through the ordinary history.

use std::collections::HashMap;

use ogeom_algo::{Built, History};
use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_topo::{Model, NodeData, Orientation, Shape, ShapeType, TShapeId};

/// A batch of substitutions, applied in one rebuild.
#[derive(Debug, Default)]
pub struct Reshape {
    /// `None` removes the node; `Some` replaces it.
    requests: HashMap<TShapeId, Option<Shape>>,
}

impl Reshape {
    /// An empty batch.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Replace every occurrence of `old` with `new`.
    pub fn replace(&mut self, old: &Shape, new: Shape) {
        self.requests.insert(old.node(), Some(new));
    }

    /// Remove every occurrence of `old`.
    pub fn remove(&mut self, old: &Shape) {
        self.requests.insert(old.node(), None);
    }

    /// How many replacements are staged.
    #[must_use]
    pub fn len(&self) -> usize {
        self.requests.len()
    }

    /// Whether anything is requested.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.requests.is_empty()
    }

    /// Apply the batch over `shape`, rebuilding what the substitutions
    /// touch and sharing everything they do not.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
    /// a rebuilt container ends up empty where the model forbids it, or a
    /// substitution's kind does not fit its slot.
    pub fn apply(&self, model: &mut Model, shape: &Shape) -> OgeomResult<Built> {
        let mut memo: HashMap<TShapeId, Option<Shape>> = HashMap::new();
        let mut history = History::new();
        let Some(result) = self.rebuilt(model, shape, &mut memo, &mut history)? else {
            ogeom_bail!(Construction, "the reshape removed the shape itself");
        };
        for (old, request) in &self.requests {
            let stand_in = Shape::of(*old);
            match request {
                Some(new) => history.modify(&stand_in, new.clone()),
                None => history.delete(&stand_in),
            }
        }
        Ok(Built::new(result, history))
    }

    /// The rebuilt occurrence of `shape`, `None` if it is removed.
    fn rebuilt(
        &self,
        model: &mut Model,
        shape: &Shape,
        memo: &mut HashMap<TShapeId, Option<Shape>>,
        history: &mut History,
    ) -> OgeomResult<Option<Shape>> {
        // A direct request wins, orientation carried from the occurrence.
        if let Some(request) = self.requests.get(&shape.node()) {
            return Ok(request.as_ref().map(|new| {
                if shape.orientation() == Orientation::Reversed {
                    new.reversed()
                } else {
                    new.clone()
                }
            }));
        }
        if let Some(held) = memo.get(&shape.node()) {
            return Ok(held.as_ref().map(|new| {
                if shape.orientation() == Orientation::Reversed {
                    new.reversed()
                } else {
                    new.clone()
                }
            }));
        }

        // Rebuild children; if none changed, the node itself is shared.
        let children = model.children_of(shape)?;
        let mut rebuilt_children = Vec::with_capacity(children.len());
        let mut changed = false;
        for child in &children {
            match self.rebuilt(model, child, memo, history)? {
                Some(new) => {
                    if new.node() != child.node() {
                        changed = true;
                    }
                    rebuilt_children.push(new);
                }
                None => changed = true,
            }
        }
        if !changed {
            memo.insert(shape.node(), Some(shape.clone()));
            return Ok(Some(shape.clone()));
        }

        let kind = model.kind_of(shape)?;
        let data = {
            let Some(node) = model.node(shape) else {
                ogeom_bail!(Construction, "shape is not in this model");
            };
            node.data().clone()
        };
        let fresh = match (kind, data) {
            (_, NodeData::Face(face)) => {
                if rebuilt_children.is_empty() {
                    memo.insert(shape.node(), None);
                    return Ok(None);
                }
                model.add_face(*face, &rebuilt_children)?
            }
            (_, NodeData::Edge(edge)) => {
                if rebuilt_children.is_empty() {
                    memo.insert(shape.node(), None);
                    return Ok(None);
                }
                model.add_edge(*edge, &rebuilt_children)?
            }
            (ShapeType::Wire, NodeData::Container) => {
                if rebuilt_children.is_empty() {
                    memo.insert(shape.node(), None);
                    return Ok(None);
                }
                model.add_wire(&rebuilt_children)?
            }
            (ShapeType::Shell, NodeData::Container) => {
                if rebuilt_children.is_empty() {
                    memo.insert(shape.node(), None);
                    return Ok(None);
                }
                model.add_shell(&rebuilt_children)?
            }
            (ShapeType::Solid, NodeData::Container) => {
                if rebuilt_children.is_empty() {
                    memo.insert(shape.node(), None);
                    return Ok(None);
                }
                model.add_solid(&rebuilt_children)?
            }
            (ShapeType::Compound, NodeData::Container) => model.add_compound(&rebuilt_children)?,
            (other, _) => {
                ogeom_bail!(
                    Construction,
                    "a {other:?} cannot be rebuilt around substituted children"
                );
            }
        };
        let oriented = if shape.orientation() == Orientation::Reversed {
            fresh.reversed()
        } else {
            fresh.clone()
        };
        history.modify(shape, fresh.clone());
        memo.insert(shape.node(), Some(fresh));
        Ok(Some(oriented))
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use ogeom_core::Tolerances;
    use ogeom_math::{Frame, Point};
    use ogeom_topo::explore_unique;

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn a_substituted_vertex_ripples_to_the_top_and_shares_the_rest() {
        let mut model = Model::new();
        let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T)
            .unwrap()
            .shape;
        let vertices = explore_unique(&model, &solid, ShapeType::Vertex).unwrap();
        let corner = vertices
            .iter()
            .find(|v| {
                model
                    .node(v)
                    .and_then(|n| n.data().as_vertex())
                    .is_some_and(|data| data.point.distance(Point::ORIGIN) < 1e-9)
            })
            .expect("the origin corner")
            .clone();
        let moved = ogeom_algo::make_vertex(&mut model, Point::new(0.1, 0.0, 0.0)).shape;

        let mut reshape = Reshape::new();
        reshape.replace(&corner, moved.clone());
        let rebuilt = reshape.apply(&mut model, &solid).unwrap();

        // The new solid holds the new vertex and none of the old one.
        let after = explore_unique(&model, &rebuilt.shape, ShapeType::Vertex).unwrap();
        assert!(after.iter().any(|v| v.node() == moved.node()));
        assert!(after.iter().all(|v| v.node() != corner.node()));
        // Only the three faces at that corner rebuilt; the rest shared.
        let before_faces = explore_unique(&model, &solid, ShapeType::Face).unwrap();
        let after_faces = explore_unique(&model, &rebuilt.shape, ShapeType::Face).unwrap();
        let shared = after_faces
            .iter()
            .filter(|f| before_faces.iter().any(|b| b.node() == f.node()))
            .count();
        assert_eq!(shared, 3, "the untouched half of the box is the same nodes");
        // History names the substitution and the ripple.
        assert!(!rebuilt.history.modified(&solid).is_empty());
    }
}
