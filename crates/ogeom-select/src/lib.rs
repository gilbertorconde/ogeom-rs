//! Spatial acceleration and selection.
//!
//! BVH construction over tessellated and analytic geometry; ray picking with
//! sub-shape granularity; rectangle and polygon selection; depth sorting; level of
//! detail; and a stable mapping from a triangle back to the topological entity that
//! produced it.
//!
//! Not a renderer. `tools/ogeom-view` consumes this crate; it is not part of it.

mod pick;

pub use pick::{FaceDraft, FaceThickness, Hit, Marquee, PickHierarchy, PickKind, Pickable, Ray};
