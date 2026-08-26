//! Shape healing — analysis, fixing and upgrading.
//!
//! *Elsewhere:* `ShapeFix`, `ShapeAnalysis`, `ShapeUpgrade`, `ShapeBuild`,
//! `ShapeConstruct` and `ShapeExtend`.
//!
//! Not optional and not deferrable: essentially every real-world STEP or IGES file
//! needs healing before it can be modeled with. A kernel that cannot survive
//! imperfect imported geometry is unusable regardless of how good its booleans
//! are.

pub mod fix;
pub mod reanchor;
pub mod same_parameter;

pub mod canonical;
pub mod reshape;
pub mod upgrade;

pub use canonical::{CanonicalReport, Simplified, canonical_simplify, recognize_surface};
pub use fix::{FixedTrims, ReanchoredBoundaries, fix_face_pcurves, reanchor_boundaries};
pub use reanchor::reanchor_periodic_rings;
pub use reshape::Reshape;
pub use same_parameter::{SameParameterReport, repair_same_parameter};
pub use upgrade::{merge_edges, reduce_tolerances, unify_same_domain};
