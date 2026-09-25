//! Shape healing: analysis, fixing and upgrading.
//!
//! *Elsewhere:* `ShapeFix`, `ShapeAnalysis`, `ShapeUpgrade`, `ShapeBuild`,
//! `ShapeConstruct` and `ShapeExtend`.
//!
//! Not optional and not deferrable: essentially every real-world STEP or IGES file
//! needs healing before it can be modeled with. A kernel that cannot survive
//! imperfect imported geometry is unusable regardless of how good its booleans
//! are.

pub mod divide;
pub mod fix;
pub mod fix_shape;
pub mod reanchor;
pub mod same_parameter;
pub mod small;

pub mod canonical;
pub mod reshape;
pub mod upgrade;

pub use canonical::{CanonicalReport, Simplified, canonical_simplify, recognize_surface};
pub use divide::{
    IsoLine, divide_by_angle, divide_by_area, divide_by_continuity, divide_face, to_bezier,
};
pub use fix::{FixedTrims, ReanchoredBoundaries, fix_face_pcurves, reanchor_boundaries};
pub use fix_shape::{FixReport, Fixed, fix_shape};
pub use reanchor::reanchor_periodic_rings;
pub use reshape::Reshape;
pub use same_parameter::{SameParameterReport, repair_same_parameter};
pub use small::{SmallFaces, fix_small_faces, remove_small_solids};
pub use upgrade::{merge_edges, reduce_tolerances, unify_same_domain};
