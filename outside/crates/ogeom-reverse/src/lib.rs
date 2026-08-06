//! Reverse engineering: mesh → B-rep, and the canonical recognition that
//! decides what the curved parts of a mesh *are*.
//!
//! This is the inverse of tessellation, and it is not the kernel's job. A
//! kernel takes exact geometry and produces triangles; recovering exact
//! geometry from triangles is a separate discipline, with a separate failure
//! mode — every answer it gives is a *decision* about information the mesh no
//! longer carries, and the reference kernel's modelling modules do not attempt
//! one. See `docs/SCOPE.md`.
//!
//! Two halves. [`reconstruct`] does the topology: it groups triangles into
//! regions, decides which chains of mesh edges were one curve, and builds
//! wires, faces, a shell and a solid from them. [`canonical`] does the
//! geometry: given samples with normals, it decides whether a patch *is* a
//! plane, a cylinder, a cone, a sphere or a torus — and answers nothing rather
//! than fitting one that merely resembles it.
//!
//! Both halves refuse by name rather than guessing. A mesh no segmentation of
//! which recognizes comes back as an error saying so, because a solid built on
//! a guessed surface looks right, measures nearly right, and is wrong
//! underneath every operation that follows.

pub mod canonical;
pub mod reconstruct;

pub use canonical::{Canonical, Recognized, mesh_to_brep, recognize_points, recognize_surface};
pub use reconstruct::{CREASE, Region, SurfaceRecognizer, planar_regions, to_brep, to_brep_with};
