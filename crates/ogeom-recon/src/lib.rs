//! A B-rep solid from a triangle mesh: the mesh's regions recognized as
//! planes, canonical surfaces, extrusions, revolutions or free-form
//! patches, and built into the faces of a solid ([`solid_from_mesh`]).
//!
//! Also the cut of a shape by a region of a mesh
//! ([`cut_by_mesh_region`]), which builds its slab through
//! [`solid_from_mesh`].

// Debug output goes through `ogeom_core::diag`, never to stdout or stderr.
#![cfg_attr(
    not(test),
    warn(clippy::print_stderr, clippy::print_stdout, clippy::dbg_macro)
)]

mod mesh_region;
pub mod mesh_solid;
mod recognize_patch;
mod recognize_swept;

pub use mesh_region::cut_by_mesh_region;
pub use mesh_solid::{
    FallbackReason, FitConstraints, MeshRegion, MeshRegions, MeshSolid, MeshSolidOptions,
    MeshSolidReport, RegionFallback, RegionId, RegionRefusal, SurfaceKind, refine_solid,
    solid_from_mesh,
};
