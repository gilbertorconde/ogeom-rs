//! The reshape framework: record substitutions, apply them over a shape in
//! one pass.
//!
//! Healing operations want to say "this edge becomes that one, this face
//! goes away" and have the change ripple upward (every wire holding the
//! edge rebuilt, every face holding the wire, up to the solid) without
//! each fix reimplementing the traversal. A [`Reshape`] collects the
//! requests; [`Reshape::apply`] rebuilds bottom-up, sharing rebuilt nodes
//! so a substituted edge is one new node however many faces reach it, and
//! reports what became of every input through the ordinary history.
//!
//! A node can stand in a shape at more than one placement (an assembly
//! placing one part twice). A request on a node ([`Reshape::replace`]) is
//! stated in the frame of the occurrence it names, and every occurrence of
//! the node takes it, carried by the placement between the two: each copy
//! changes alike. A request on one occurrence ([`Reshape::replace_at`])
//! changes that copy alone; every other keeps its nodes. A node rebuilt
//! around changed children is made in its own frame and placed where its
//! occurrence stands, and occurrences rebuilt alike share the new node.

use ogeom_core::FastMap;

use ogeom_algo::{Built, History};
use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_topo::{Location, Model, NodeData, Orientation, SameKey, Shape, ShapeType, TShapeId};

/// A batch of substitutions, applied in one rebuild.
#[derive(Debug, Default)]
pub struct Reshape {
    /// Per node, the occurrence named and what it becomes: `None` removes
    /// it, `Some` replaces it, stated in the named occurrence's frame.
    requests: FastMap<TShapeId, (Location, Option<Shape>)>,
    /// Per occurrence alone (node and placement), what it becomes, stated
    /// in the frame of the shape the batch is applied to.
    at: FastMap<SameKey, Option<Shape>>,
    /// Nodes cut into pieces, in their forward order, with the placement
    /// of the occurrence named.
    splits: FastMap<TShapeId, (Location, Vec<Shape>)>,
}

impl Reshape {
    /// An empty batch.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Replace every occurrence of `old`'s node with `new`, stated in the
    /// frame `old` is given in: an occurrence placed elsewhere takes `new`
    /// carried by the placement taking `old` onto it.
    pub fn replace(&mut self, old: &Shape, new: Shape) {
        self.requests
            .insert(old.node(), (old.location().clone(), Some(new)));
    }

    /// Replace every occurrence of `old`'s node with `pieces`, in order: an
    /// edge cut in a wire, a face cut in a shell. A reversed occurrence
    /// takes the pieces reversed, in reverse order. The pieces are stated
    /// in `old`'s frame, as in [`Reshape::replace`].
    pub fn split(&mut self, old: &Shape, pieces: Vec<Shape>) {
        self.splits
            .insert(old.node(), (old.location().clone(), pieces));
    }

    /// Remove every occurrence of `old`'s node.
    pub fn remove(&mut self, old: &Shape) {
        self.requests
            .insert(old.node(), (old.location().clone(), None));
    }

    /// Replace `occurrence` alone, at the placement it has in the shape
    /// the batch is applied to, with `new`, stated in that shape's frame.
    /// Other occurrences of the node keep it. Takes precedence over a
    /// request on the node.
    pub fn replace_at(&mut self, occurrence: &Shape, new: Shape) {
        self.at.insert(SameKey(occurrence.clone()), Some(new));
    }

    /// Remove `occurrence` alone, as [`Reshape::replace_at`] names it.
    pub fn remove_at(&mut self, occurrence: &Shape) {
        self.at.insert(SameKey(occurrence.clone()), None);
    }

    /// How many replacements are staged.
    #[must_use]
    pub fn len(&self) -> usize {
        self.requests.len() + self.splits.len() + self.at.len()
    }

    /// Whether anything is requested.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.requests.is_empty() && self.splits.is_empty() && self.at.is_empty()
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
        let mut memo = Memo::default();
        let mut history = History::new();
        let Some(result) = self.rebuilt(model, shape, &mut memo, &mut history)? else {
            ogeom_bail!(Construction, "the reshape removed the shape itself");
        };
        for (old, (from, request)) in &self.requests {
            let stand_in = Shape::new(*old, from.clone(), Orientation::Forward);
            match request {
                Some(new) => history.modify(&stand_in, new.clone()),
                None => history.delete(&stand_in),
            }
        }
        for (occurrence, request) in &self.at {
            match request {
                Some(new) => history.modify(&occurrence.0, new.clone()),
                None => history.delete(&occurrence.0),
            }
        }
        for (old, (from, pieces)) in &self.splits {
            let stand_in = Shape::new(*old, from.clone(), Orientation::Forward);
            for piece in pieces {
                history.modify(&stand_in, piece.clone());
            }
        }
        Ok(Built::new(result, history))
    }

    /// The rebuilt occurrence of `shape`, `None` if it is removed.
    fn rebuilt(
        &self,
        model: &mut Model,
        shape: &Shape,
        memo: &mut Memo,
        history: &mut History,
    ) -> OgeomResult<Option<Shape>> {
        let oriented = |new: &Shape| {
            if shape.orientation() == Orientation::Reversed {
                new.reversed()
            } else {
                new.clone()
            }
        };
        // A direct request wins, orientation carried from the occurrence.
        let key = SameKey(shape.clone());
        if let Some(request) = self.at.get(&key) {
            return Ok(request.as_ref().map(oriented));
        }
        if let Some((from, request)) = self.requests.get(&shape.node()) {
            return Ok(request
                .as_ref()
                .map(|new| oriented(&carried(new, from, shape.location()))));
        }
        if let Some(held) = memo.done.get(&key) {
            return Ok(held.as_ref().map(oriented));
        }

        // Rebuild children; if none changed, the node itself is shared.
        // Read through the forward occurrence: `children_of` composes the
        // occurrence's orientation into every child, and the rebuilt node is
        // oriented as the occurrence once, below. Read through a reversed
        // occurrence, a reversed edge rebuilt around a substituted vertex
        // comes out with its start and end swapped and is then reversed
        // again, so the wire it sits in walks it against its neighbours.
        let forward = if shape.orientation() == Orientation::Reversed {
            shape.reversed()
        } else {
            shape.clone()
        };
        let children = model.children_of(&forward)?;
        let mut rebuilt_children = Vec::with_capacity(children.len());
        let mut changed = false;
        for child in &children {
            if let Some((from, pieces)) = self.splits.get(&child.node()) {
                changed = true;
                let pieces = pieces
                    .iter()
                    .map(|piece| carried(piece, from, child.location()));
                if child.orientation() == Orientation::Reversed {
                    let pieces: Vec<Shape> = pieces.collect();
                    rebuilt_children.extend(pieces.iter().rev().map(Shape::reversed));
                } else {
                    rebuilt_children.extend(pieces);
                }
                continue;
            }
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
            memo.done.insert(key, Some(forward));
            return Ok(Some(shape.clone()));
        }

        // The node is made in its own frame, its children read relative to
        // the occurrence, and placed where the occurrence stands.
        let own = forward.location().clone();
        if !own.is_identity() {
            let back = own.inverted();
            for child in &mut rebuilt_children {
                *child = Shape::new(
                    child.node(),
                    back.then(child.location()),
                    child.orientation(),
                );
            }
        }
        let alike = (shape.node(), rebuilt_children);
        if let Some(fresh) = memo.made.get(&alike) {
            let placed = fresh.located(own);
            history.modify(shape, placed.clone());
            memo.done.insert(key, Some(placed.clone()));
            return Ok(Some(oriented(&placed)));
        }
        let rebuilt_children = &alike.1;

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
                    memo.done.insert(key, None);
                    return Ok(None);
                }
                model.add_face(*face, rebuilt_children)?
            }
            (_, NodeData::Edge(edge)) => {
                if rebuilt_children.is_empty() {
                    memo.done.insert(key, None);
                    return Ok(None);
                }
                model.add_edge(*edge, rebuilt_children)?
            }
            (ShapeType::Wire, NodeData::Container) => {
                if rebuilt_children.is_empty() {
                    memo.done.insert(key, None);
                    return Ok(None);
                }
                model.add_wire(rebuilt_children)?
            }
            (ShapeType::Shell, NodeData::Container) => {
                if rebuilt_children.is_empty() {
                    memo.done.insert(key, None);
                    return Ok(None);
                }
                model.add_shell(rebuilt_children)?
            }
            (ShapeType::Solid, NodeData::Container) => {
                if rebuilt_children.is_empty() {
                    memo.done.insert(key, None);
                    return Ok(None);
                }
                model.add_solid(rebuilt_children)?
            }
            (ShapeType::Compound, NodeData::Container) => model.add_compound(rebuilt_children)?,
            (other, _) => {
                ogeom_bail!(
                    Construction,
                    "a {other:?} cannot be rebuilt around substituted children"
                );
            }
        };
        memo.made.insert(alike, fresh.clone());
        let placed = if own.is_identity() {
            fresh
        } else {
            fresh.located(own)
        };
        history.modify(shape, placed.clone());
        memo.done.insert(key, Some(placed.clone()));
        Ok(Some(oriented(&placed)))
    }
}

/// What [`Reshape::rebuilt`] has answered: per occurrence, its rebuilt
/// form; per node and its children as rebuilt (relative to it), the node
/// made, so occurrences rebuilt alike share it.
#[derive(Default)]
struct Memo {
    done: FastMap<SameKey, Option<Shape>>,
    made: FastMap<(TShapeId, Vec<Shape>), Shape>,
}

/// `new`, stated in the frame of an occurrence at `from`, as it stands for
/// an occurrence at `to`.
fn carried(new: &Shape, from: &Location, to: &Location) -> Shape {
    let by = to.then(&from.inverted());
    if by.is_identity() {
        return new.clone();
    }
    Shape::new(new.node(), by.then(new.location()), new.orientation())
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use ogeom_core::Tolerances;
    use ogeom_geom::{PlaneSurface, SurfaceGeometry};
    use ogeom_math::{Frame, Plane, Point, Transform, Vector};
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

    /// A 2 mm square face and a copy of it placed 10 above, in one
    /// compound; the copy, and the square's corner at the origin.
    fn square_and_copy(model: &mut Model) -> (Shape, Shape, Shape) {
        let corners =
            [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)].map(|(x, y)| Point::new(x, y, 0.0));
        let ring = ogeom_algo::make_polygon(model, &corners, true, T)
            .unwrap()
            .shape;
        let plane = SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(Frame::WORLD)));
        let face = ogeom_algo::make_face(model, plane, &[ring], T)
            .unwrap()
            .shape;
        let up = model.placed(&face, Transform::translation(Vector::new(0.0, 0.0, 10.0)));
        let corner = explore_unique(model, &face, ShapeType::Vertex)
            .unwrap()
            .into_iter()
            .find(|v| placed_point(model, v).distance(Point::ORIGIN) < 1e-12)
            .unwrap();
        let root = model.add_compound(&[face, up.clone()]).unwrap();
        (root, up, corner)
    }

    fn placed_point(model: &Model, vertex: &Shape) -> Point {
        let point = model
            .node(vertex)
            .unwrap()
            .data()
            .as_vertex()
            .unwrap()
            .point;
        vertex.transform(model.datums()).unwrap().apply(point)
    }

    /// Every vertex of `shape` where it stands, by height then x.
    fn points(model: &Model, shape: &Shape) -> Vec<(f64, f64)> {
        let mut out: Vec<(f64, f64)> = explore_unique(model, shape, ShapeType::Vertex)
            .unwrap()
            .iter()
            .map(|v| placed_point(model, v))
            .map(|p| (p.z, p.x + 10.0 * p.y))
            .collect();
        out.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
        out
    }

    /// A vertex replaced by node reaches both copies, each where it stands,
    /// and the copies share the face rebuilt around it.
    #[test]
    fn a_node_replaced_changes_every_placed_copy_where_it_stands() {
        let mut model = Model::new();
        let (root, _, corner) = square_and_copy(&mut model);
        let moved = ogeom_algo::make_vertex(&mut model, Point::new(0.5, 0.0, 0.0)).shape;
        let mut reshape = Reshape::new();
        reshape.replace(&corner, moved);
        let built = reshape.apply(&mut model, &root).unwrap().shape;
        let want = [0.5, 2.0, 20.0, 22.0];
        let got = points(&model, &built);
        assert_eq!(got.len(), 8);
        for (k, (z, at)) in got.iter().enumerate() {
            let (wz, wat) = (if k < 4 { 0.0 } else { 10.0 }, want[k % 4]);
            assert!(
                (z - wz).abs() < 1e-12 && (at - wat).abs() < 1e-12,
                "{got:?}"
            );
        }
        let faces = explore_unique(&model, &built, ShapeType::Face).unwrap();
        assert_eq!(faces.len(), 2);
        assert_eq!(faces[0].node(), faces[1].node());
    }

    /// A vertex replaced at one copy alone changes that copy, where it
    /// stands; the other keeps its nodes and its corner.
    #[test]
    fn an_occurrence_replaced_changes_that_copy_alone() {
        let mut model = Model::new();
        let (root, up, corner) = square_and_copy(&mut model);
        let at_copy = corner.moved(up.location());
        let moved = ogeom_algo::make_vertex(&mut model, Point::new(0.5, 0.0, 0.0)).shape;
        let mut reshape = Reshape::new();
        reshape.replace_at(&at_copy, moved.located(up.location().clone()));
        let built = reshape.apply(&mut model, &root).unwrap().shape;
        let got = points(&model, &built);
        let want = [(0.0, 0.0), (0.0, 2.0), (0.0, 20.0), (0.0, 22.0)]
            .into_iter()
            .chain([0.5, 2.0, 20.0, 22.0].map(|x| (10.0, x)));
        for ((z, at), (wz, wat)) in got.iter().zip(want) {
            assert!(
                (z - wz).abs() < 1e-12 && (at - wat).abs() < 1e-12,
                "{got:?}"
            );
        }
        let faces = explore_unique(&model, &built, ShapeType::Face).unwrap();
        let below = faces.iter().find(|f| f.location().is_identity()).unwrap();
        let original = &model.children_of(&root).unwrap()[0];
        assert_eq!(
            below.node(),
            original.node(),
            "the copy below keeps its nodes"
        );
        assert_ne!(faces[0].node(), faces[1].node());
    }
}
