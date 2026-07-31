//! Placing and duplicating shapes.
//!
//! # A rigid move copies nothing
//!
//! [`transformed`] returns the *same* topology at a different
//! [`Location`]. No node is created, no curve is
//! re-evaluated, and the result compares equal to the original under
//! [`Shape::is_partner`] — which is how "these thousand bolts are the same
//! bolt" stays a fact the model knows rather than one an application has to
//! remember (`docs/DATA_MODEL.md` §2, §3).
//!
//! That is only sound because a location is a rigid motion with a uniform
//! scale. Such a motion carries a line to a line and a circle to a circle, so
//! the geometry underneath still describes the moved shape. An affine transform
//! that shears or scales unevenly does not: it carries a circle to an ellipse,
//! and no amount of placement makes a circle record that. The type system says
//! so — [`transformed`] takes a [`Transform`], which is a similarity by
//! construction, and a general affine transform is a different type it will not
//! accept. Applying one means rebuilding the geometry, which is not written
//! yet; see the deferred list in `docs/SCOPE.md`.
//!
//! # A copy is for editing, not for moving
//!
//! [`copied`] duplicates the topology so the two can diverge. It shares the
//! *geometry* — curves and surfaces are immutable values in an arena, so two
//! shapes naming one circle can never disagree about it, and copying it would
//! only make the model larger. What a copy buys is independent topology:
//! tolerances, representations and children that one shape can change without
//! the other seeing it.

use std::collections::HashMap;

use og_core::{OgResult, og_bail};
use og_math::Transform;
use og_topo::{Location, Model, NodeData, Shape, ShapeType, TShapeId};

use crate::history::{Built, History};

/// Roles this module assigns.
pub mod roles {
    use og_core::Role;

    /// An entity that is a copy of another.
    pub const COPY: Role = Role::op_defined(30);
}

/// Move a shape by a rigid motion, sharing everything.
///
/// Cheap and exact: the result names the same topology nodes and the same
/// geometry, at a new placement.
///
/// There is no way to pass something that is *not* a placement. [`Transform`]
/// is a similarity by construction — rigid motion with a uniform scale — and a
/// shear or a non-uniform scale is a
/// [`GeneralTransform`](og_math::GeneralTransform), which this does not accept.
/// That is deliberate: such a transform carries a circle to an ellipse, and
/// recording it as a placement would leave every circle in the shape claiming
/// to be a circle while sitting on an ellipse. The type refuses it, so no
/// runtime check has to.
///
/// # Errors
///
/// [`OgError::Dangling`](og_core::OgError::Dangling) if the shape does not
/// resolve in this model.
pub fn transformed(model: &mut Model, shape: &Shape, transform: Transform) -> OgResult<Built> {
    if model.node(shape).is_none() {
        og_bail!(Dangling, "shape refers to a node not in this model");
    }
    model.begin_operation();
    let datum = model.add_datum(transform);
    let moved = shape.moved(&Location::of(datum));

    let mut history = History::new();
    history.modify(shape, moved.clone());
    Ok(Built::new(moved, history))
}

/// Duplicate a shape's topology so the two can be edited apart.
///
/// Geometry is shared, not duplicated — see the module docs.
///
/// # Errors
///
/// [`OgError::Dangling`](og_core::OgError::Dangling) if the shape does not
/// resolve in this model.
pub fn copied(model: &mut Model, shape: &Shape) -> OgResult<Built> {
    if model.node(shape).is_none() {
        og_bail!(Dangling, "shape refers to a node not in this model");
    }
    model.begin_operation();

    let mut done: HashMap<TShapeId, Shape> = HashMap::new();
    let mut history = History::new();
    let root = duplicate(model, shape, &mut done, &mut history)?;
    Ok(Built::new(root, history))
}

/// Copy one node and everything below it, memoized.
///
/// The memo is not an optimization. A shared edge appears under two faces, and
/// copying it twice would give the copy two edges where the original had one —
/// the shell would then be open along every shared boundary, and nothing about
/// the geometry would say why.
fn duplicate(
    model: &mut Model,
    shape: &Shape,
    done: &mut HashMap<TShapeId, Shape>,
    history: &mut History,
) -> OgResult<Shape> {
    if let Some(existing) = done.get(&shape.node()) {
        return Ok(existing
            .moved(shape.location())
            .composed(shape.orientation()));
    }

    let Some(node) = model.node(shape) else {
        og_bail!(Dangling, "shape refers to a node not in this model");
    };
    let kind = node.kind();
    let data = node.data().clone();

    // Children first: a node is built from the shapes below it, so they have to
    // exist before it does.
    let mut children = Vec::new();
    for child in model.children_of(shape)? {
        children.push(duplicate(model, &child, done, history)?);
    }

    let fresh = match (kind, data) {
        (ShapeType::Vertex, NodeData::Vertex(v)) => model.add_vertex(v),
        (ShapeType::Edge, NodeData::Edge(e)) => model.add_edge(*e, &children)?,
        (ShapeType::Wire, _) => model.add_wire(&children)?,
        (ShapeType::Face, NodeData::Face(f)) => model.add_face(*f, &children)?,
        (ShapeType::Shell, _) => model.add_shell(&children)?,
        (ShapeType::Solid, _) => model.add_solid(&children)?,
        (ShapeType::CompSolid, _) => model.add_compsolid(&children)?,
        (ShapeType::Compound, _) => model.add_compound(&children)?,
        (other, _) => og_bail!(
            Construction,
            "a {other:?} node does not hold the data its kind requires, so it \
             cannot be copied"
        ),
    };

    // The copy names what it came from, so a reference into the original still
    // resolves after the copy is edited.
    let bare = Shape::of(fresh.node());
    model.set_derived(&bare, std::slice::from_ref(shape), roles::COPY)?;
    history.modify(shape, bare.clone());
    done.insert(shape.node(), bare.clone());

    Ok(bare.moved(shape.location()).composed(shape.orientation()))
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::check::check;
    use crate::mass::volume_properties;
    use crate::{make_box, make_cylinder};
    use approx::assert_relative_eq;
    use og_core::Tolerances;
    use og_math::{Axis, Direction, Frame, Point, Vector};
    use og_mesh::Deflection;
    use og_topo::explore_unique;

    const T: Tolerances = Tolerances::millimetres();

    fn deflection() -> Deflection {
        Deflection {
            chord: 0.01,
            ..Deflection::default()
        }
    }

    #[test]
    fn a_rigid_move_creates_no_topology_at_all() {
        // The whole point. A thousand identical bolts should cost one bolt and
        // a thousand placements, not a thousand bolts.
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T)
            .unwrap()
            .shape;
        let before = model.node_count();

        let moved = transformed(&mut model, &solid, Transform::translation(Vector::X * 10.0))
            .unwrap()
            .shape;

        assert_eq!(model.node_count(), before, "a placement copied something");
        assert!(moved.is_partner(&solid), "the same topology, elsewhere");
        assert!(!moved.is_same(&solid), "but at a different placement");
    }

    #[test]
    fn a_moved_shape_measures_the_same_and_sits_elsewhere() {
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T)
            .unwrap()
            .shape;
        let offset = Vector::new(10.0, -20.0, 30.0);
        let moved = transformed(&mut model, &solid, Transform::translation(offset))
            .unwrap()
            .shape;

        let here = volume_properties(&model, &solid, deflection(), T).unwrap();
        let there = volume_properties(&model, &moved, deflection(), T).unwrap();
        assert_relative_eq!(here.mass, there.mass, epsilon = 1e-9);
        assert!(there.centre.distance(here.centre + offset) < 1e-9);

        assert!(check(&model, &moved, T).unwrap().is_valid());
    }

    #[test]
    fn a_rotation_is_a_placement_and_a_reflection_still_is() {
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T)
            .unwrap()
            .shape;
        for transform in [
            Transform::rotation(Axis::Z, 0.7),
            Transform::scaling(Point::ORIGIN, 2.0, T).unwrap(),
            Transform::plane_mirror(Point::ORIGIN, Direction::Z),
        ] {
            assert!(
                transformed(&mut model, &solid, transform).is_ok(),
                "a rigid or uniformly scaled motion should be a placement"
            );
        }
    }

    #[test]
    fn a_uniform_scale_scales_the_volume_by_its_cube() {
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let bigger = transformed(
            &mut model,
            &solid,
            Transform::scaling(Point::ORIGIN, 3.0, T).unwrap(),
        )
        .unwrap()
        .shape;
        let props = volume_properties(&model, &bigger, deflection(), T).unwrap();
        assert_relative_eq!(props.mass, 27.0, epsilon = 1e-9);
    }

    #[test]
    fn a_copy_has_its_own_topology_and_the_same_geometry() {
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T)
            .unwrap()
            .shape;
        let (curves, pcurves, surfaces) = model.geometry().counts();

        let copy = copied(&mut model, &solid).unwrap().shape;
        assert!(!copy.is_partner(&solid), "a copy is not the same topology");
        assert_eq!(
            model.geometry().counts(),
            (curves, pcurves, surfaces),
            "geometry is immutable and shared; copying it buys nothing"
        );

        // And it is a whole, valid solid, not a shell of loose faces.
        assert!(check(&model, &copy, T).unwrap().is_valid());
        let props = volume_properties(&model, &copy, deflection(), T).unwrap();
        assert_relative_eq!(props.mass, 8.0, epsilon = 1e-9);
    }

    #[test]
    fn a_copy_keeps_shared_edges_shared() {
        // The memo is what makes this true. Copying each edge once per face
        // would give the copy twice the edges, and the shell would be open
        // along every one of them with nothing in the geometry to say so.
        let mut model = Model::new();
        let solid = make_cylinder(&mut model, Frame::WORLD, 2.0, 3.0, T)
            .unwrap()
            .shape;
        let copy = copied(&mut model, &solid).unwrap().shape;

        for kind in [
            ShapeType::Face,
            ShapeType::Edge,
            ShapeType::Vertex,
            ShapeType::Wire,
        ] {
            assert_eq!(
                explore_unique(&model, &copy, kind).unwrap().len(),
                explore_unique(&model, &solid, kind).unwrap().len(),
                "the copy has a different number of {kind:?}"
            );
        }
        assert!(check(&model, &copy, T).unwrap().is_valid());
    }

    #[test]
    fn a_copy_names_what_it_came_from() {
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let built = copied(&mut model, &solid).unwrap();

        assert_eq!(
            model
                .provenance_of(&built.shape)
                .and_then(og_core::Provenance::role),
            Some(roles::COPY)
        );
        // And the history says the original became it, so a reference into the
        // original still resolves.
        assert_eq!(
            built.history.modified(&solid),
            std::slice::from_ref(&built.shape)
        );
        assert!(!built.history.is_deleted(&solid));
    }

    #[test]
    fn placing_or_copying_a_stranger_is_an_error() {
        let mut other = Model::new();
        for _ in 0..4 {
            other.add_vertex(og_topo::VertexData::new(Point::ORIGIN));
        }
        let beyond = other.add_vertex(og_topo::VertexData::new(Point::ORIGIN));

        let mut empty = Model::new();
        assert!(transformed(&mut empty, &beyond, Transform::IDENTITY).is_err());
        assert!(copied(&mut empty, &beyond).is_err());
    }
}
