//! Shape healing — analysis, fixing and upgrading.
//!
//! *Elsewhere:* `ShapeFix`, `ShapeAnalysis`, `ShapeUpgrade`, `ShapeBuild`,
//! `ShapeConstruct` and `ShapeExtend`.
//!
//! Not optional and not deferrable: essentially every real-world STEP or IGES file
//! needs healing before it can be modeled with. A kernel that cannot survive
//! imperfect imported geometry is unusable regardless of how good its booleans
//! are.

pub mod reanchor;
pub mod same_parameter;

pub mod canonical;
pub mod reshape;

pub use canonical::{Canonical, Recognized, mesh_to_brep, recognize_points, recognize_surface};
pub use reanchor::reanchor_periodic_rings;
pub use reshape::Reshape;
pub use same_parameter::{SameParameterReport, repair_same_parameter};
