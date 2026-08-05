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

pub use reanchor::reanchor_periodic_rings;
