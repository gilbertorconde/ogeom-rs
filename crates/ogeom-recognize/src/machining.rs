//! From recognized features to the operations that would cut them.
//!
//! Recognition says what a shape *has*; this says what a machine would have
//! to do about it. The mapping is not clever and should not be: a hole is
//! drilled, a pocket is milled, a chamfer is cut with a chamfer tool, and a
//! round is either a tool radius or a shape the cutter leaves behind. What
//! is worth doing carefully is the *bookkeeping* — the direction each
//! operation approaches from, the largest tool that fits, and the order
//! imposed by the feature tree, since a hole through a pocket floor cannot
//! be drilled before the pocket exists.
//!
//! Every number here comes from the recognized geometry. Where a number
//! cannot be had from it — the depth of a pocket whose floor the recognizer
//! found but whose opening it did not — the field says so by being absent
//! rather than by holding a guess.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::SurfaceGeometry;
use ogeom_math::{Aabb, Direction, Vector};
use ogeom_topo::{Model, NodeData, Shape, ShapeType, explore_unique};

use crate::recognize::{Feature, FeatureNode, HoleKind};

/// What a machine would do to make one feature.
#[derive(Debug, Clone, PartialEq)]
pub enum Operation {
    /// A drilled hole.
    Drill {
        /// The bore diameter — the smallest in the chain, which is what a
        /// drill has to fit.
        diameter: f64,
        /// The direction the drill enters along.
        approach: Direction,
        /// How deep, where the geometry says. A through hole says how far
        /// it runs; a blind one, where its floor is.
        depth: f64,
        /// Whether it breaks through.
        through: bool,
        /// Whether the chain steps — a counterbore — or opens on a cone,
        /// which are second operations on the same axis.
        counterbored: bool,
        /// Whether a cone opens the entry.
        countersunk: bool,
    },
    /// A milled depression.
    Mill {
        /// The direction the cutter comes down.
        approach: Direction,
        /// How deep the floor sits below the material it was cut from,
        /// where that can be measured.
        depth: Option<f64>,
        /// The largest end mill that reaches every corner: twice the
        /// smallest inside radius the floor's boundary turns through, or
        /// `None` where the floor's corners are sharp and no round tool
        /// leaves that shape.
        largest_tool: Option<f64>,
        /// Whether the floor is an obround — a slot rather than a pocket.
        slot: bool,
    },
    /// A bevel.
    Chamfer,
    /// A blend the cutter leaves: a round of this radius is what an end
    /// mill of twice it produces at an inside corner, and what a ball or
    /// bull-nose tool produces at an outside one.
    Blend {
        /// The rolling radius.
        radius: f64,
        /// Inside corner or outside.
        concave: bool,
    },
    /// Material left standing: what is around a boss has to come off.
    ClearAround {
        /// The direction the cutter comes down.
        approach: Direction,
        /// How far it stands proud, where that can be measured.
        height: Option<f64>,
    },
}

/// One step of a plan: an operation and the feature it comes from.
#[derive(Debug, Clone)]
pub struct Step {
    /// What to do.
    pub operation: Operation,
    /// The feature it makes.
    pub feature: Feature,
    /// How many features it sits under: a hole through a pocket's floor is
    /// at depth one, and everything at depth zero can be done first.
    pub after: usize,
}

/// A machining plan: the recognized tree read as operations, parents first.
///
/// The order is the tree's own: a feature cut into another's face cannot be
/// cut before it exists. Within one level the order is the recognizer's,
/// which is deterministic.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a
/// feature's faces cannot be read out of the model.
pub fn manufacturing_plan(
    model: &Model,
    tree: &[FeatureNode],
    tol: Tolerances,
) -> OgeomResult<Vec<Step>> {
    let mut out = Vec::new();
    for node in tree {
        walk(model, node, 0, &mut out, tol)?;
    }
    Ok(out)
}

fn walk(
    model: &Model,
    node: &FeatureNode,
    after: usize,
    out: &mut Vec<Step>,
    tol: Tolerances,
) -> OgeomResult<()> {
    out.push(Step {
        operation: operation_for(model, &node.feature, tol)?,
        feature: node.feature.clone(),
        after,
    });
    for child in &node.children {
        walk(model, child, after + 1, out, tol)?;
    }
    Ok(())
}

/// The operation one feature implies.
fn operation_for(model: &Model, feature: &Feature, tol: Tolerances) -> OgeomResult<Operation> {
    Ok(match feature {
        Feature::Hole(hole) => {
            let bound = extent(model, &hole.faces)?;
            // How far the bore runs along its own axis: the extent of the
            // faces it claims, measured along it.
            let depth = along(&bound, hole.axis.direction.vector());
            Operation::Drill {
                diameter: hole.radius * 2.0,
                // A drill enters against the axis's own sense only if the
                // material is that way; the axis as recognized points along
                // the bore, and entering along it is the convention here.
                approach: hole.axis.direction,
                depth,
                through: matches!(hole.kind, HoleKind::Through),
                counterbored: hole.counterbored,
                countersunk: hole.countersunk,
            }
        }
        Feature::Pocket(pocket) => {
            let normal = face_normal(model, &pocket.floor, tol)?;
            let walls = extent(model, &pocket.walls)?;
            let floor = extent(model, std::slice::from_ref(&pocket.floor))?;
            let depth = (!pocket.walls.is_empty())
                .then(|| along(&walls, normal.vector()) - along(&floor, normal.vector()));
            Operation::Mill {
                approach: -normal,
                depth,
                largest_tool: smallest_inside_radius(model, &pocket.floor, tol)?.map(|r| r * 2.0),
                slot: pocket.slot,
            }
        }
        Feature::Boss(boss) => {
            let normal = face_normal(model, &boss.top, tol)?;
            let walls = extent(model, &boss.walls)?;
            let top = extent(model, std::slice::from_ref(&boss.top))?;
            let height = (!boss.walls.is_empty())
                .then(|| along(&top, normal.vector()) - along(&walls, normal.vector()));
            Operation::ClearAround {
                approach: -normal,
                height,
            }
        }
        Feature::Chamfer(_) => Operation::Chamfer,
        Feature::Fillet(f) => Operation::Blend {
            radius: f.radius,
            concave: f.concave,
        },
        Feature::PartialRound(r) => Operation::Blend {
            radius: r.radius,
            concave: false,
        },
    })
}

/// The world bound of a set of faces.
fn extent(model: &Model, faces: &[Shape]) -> OgeomResult<Aabb> {
    let mut bound = Aabb::EMPTY;
    for face in faces {
        for vertex in explore_unique(model, face, ShapeType::Vertex)? {
            if let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) {
                bound = bound.with_point(vertex.transform(model.datums())?.apply(data.point));
            }
        }
    }
    Ok(bound)
}

/// How far a bound reaches along a direction.
fn along(bound: &Aabb, direction: Vector) -> f64 {
    let (Some(low), Some(high)) = (bound.low(), bound.high()) else {
        return 0.0;
    };
    let mut span = 0.0f64;
    for corner in Aabb::of_corners(low, high).corners() {
        span = span.max((corner - low).dot(direction));
    }
    span
}

/// A planar face's outward normal.
fn face_normal(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Direction> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data().clone()) else {
        ogeom_core::ogeom_bail!(Construction, "expected a face");
    };
    let Some(SurfaceGeometry::Plane(plane)) = model.geometry().surface(data.surface) else {
        ogeom_core::ogeom_bail!(
            Construction,
            "a pocket floor or a boss top is planar; this one is not"
        );
    };
    let placement = face.transform(model.datums())?;
    let mut normal = placement.apply_vector(plane.plane().normal().vector());
    if face.orientation() == ogeom_topo::Orientation::Reversed {
        normal = -normal;
    }
    Direction::new(normal, tol)
}

/// The smallest radius the floor's own boundary turns through, which is
/// what limits the tool. `None` when every corner is sharp — no round tool
/// leaves that shape, and saying so is more use than a number.
fn smallest_inside_radius(
    model: &Model,
    floor: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let mut least: Option<f64> = None;
    for edge in explore_unique(model, floor, ShapeType::Edge)? {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(ogeom_topo::EdgeRepr::Curve3d { curve, .. }) = data.curve3d() else {
            continue;
        };
        let Some(ogeom_geom::Curve::Circle(circle)) = model.geometry().curve(*curve) else {
            continue;
        };
        let radius = circle.circle().radius();
        if radius > tol.confusion() && least.is_none_or(|held| radius < held) {
            least = Some(radius);
        }
    }
    Ok(least)
}
