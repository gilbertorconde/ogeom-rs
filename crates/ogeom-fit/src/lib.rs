//! Fitting and reading of measured points and meshes: primitive fits, zones
//! and region growing, surfaces fitted to curve networks, registration, and
//! the outlines, segments and resting faces read off a mesh.
//!
//! Curve fitting through points lives in `ogeom_geom::fit`; `fit` re-exports
//! it.

// Debug output goes through `ogeom_core::diag`, never to stdout or stderr.
#![cfg_attr(
    not(test),
    warn(clippy::print_stderr, clippy::print_stdout, clippy::dbg_macro)
)]

mod cutoff;
pub mod fit;
pub mod held;
mod hull3;
pub mod network;
pub mod orient;
pub mod outline;
pub mod recognize;
pub mod register;
pub mod rest;
pub mod segments_2d;

pub use cutoff::Cutoff;
pub use held::{NetOutside, net_held_within};
pub use network::{
    Across, ChainedPatch, CurveDef, IrregularCorner, LoopSide, Network, NetworkOptions, PatchDef,
    PatchFit, PatchSide, fit_loop_quads, fit_network, fit_network_chained,
};
pub use outline::{Outline, Piece2, Wire2d, fit_outline};
pub use recognize::{Canonical, Recognized, recognize_points};
pub use rest::{Rest, resting_faces, resting_faces_of};
pub use segments_2d::{Segment2, Segments, fit_segments_2d};
