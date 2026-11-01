//! Picking and selection.
//!
//! Ray picking over the kernel's triangle hierarchy (`ogeom_mesh::Bvh`) with
//! sub-shape granularity; rectangle and polygon selection; depth sorting; level of
//! detail; and a stable mapping from a triangle back to the topological entity that
//! produced it.
//!
//! Not a renderer. `tools/ogeom-view` consumes this crate; it is not part of it.

mod pick;

pub use pick::{FaceDraft, FaceThickness, Hit, Marquee, PickHierarchy, PickKind, Pickable, Ray};
