//! Feature recognition: reading design intent back out of raw topology.
//!
//! An imported solid is dumb — faces and edges with no memory of the
//! operations that made them. But the operations leave fingerprints: a
//! drilled hole is a concave full cylinder, a fillet is a blend tangent to
//! both its neighbours, a chamfer is a bevel that is tangent to neither, a
//! pocket is a floor whose every boundary edge folds inward. Recognition
//! walks the topology, measures those fingerprints, and returns a feature
//! list — the beginning of a feature tree for a model that arrived without
//! one.
//!
//! Convexity is the classical reading: at a shared edge, the sign of
//! `(n1 x n2) . t` — the two outward normals against the edge as the first
//! face traverses it — says whether the surface turns outward or folds in,
//! and parallel normals are a smooth join with no fold to name.

mod defeature;
mod machining;
mod recognize;

pub use defeature::remove_feature;
pub use machining::{Operation, Step, manufacturing_plan};
pub use recognize::{
    Boss, Chamfer, Feature, FeatureNode, Fillet, Hole, HoleKind, PartialRound, Pocket,
    feature_tree, recognize,
};
