//! ogeom: a boundary-representation CAD kernel.
//!
//! This umbrella crate is the public API. Depend on this, not on the `og-*` crates
//! individually — their boundaries are an implementation detail and will move.
//!
//! A from-scratch Rust implementation: no binding, no wrapper, no dependency on
//! any existing CAD kernel. It keeps the *data model* the field converged on —
//! shared topology with location chains, per-entity tolerances,
//! multi-representation edges, operation history — because that model is what
//! makes the algorithms possible, and discards the 1990s C++ around it. See
//! `docs/DATA_MODEL.md`.

pub use ogeom_algo as algo;
pub use ogeom_bool as boolean;
pub use ogeom_core as core;
pub use ogeom_doc as doc;
pub use ogeom_fillet as fillet;
pub use ogeom_geom as geom;
pub use ogeom_heal as heal;
pub use ogeom_hlr as hlr;
pub use ogeom_intersect as intersect;
pub use ogeom_io as io;
pub use ogeom_math as math;
pub use ogeom_mesh as mesh;
pub use ogeom_offset as offset;
pub use ogeom_topo as topo;

/// The kernel version, as reported by front ends.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
