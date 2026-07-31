//! openGeometry: a boundary-representation CAD kernel.
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

pub use og_algo as algo;
pub use og_bool as boolean;
pub use og_core as core;
pub use og_doc as doc;
pub use og_fillet as fillet;
pub use og_geom as geom;
pub use og_heal as heal;
pub use og_hlr as hlr;
pub use og_intersect as intersect;
pub use og_io as io;
pub use og_math as math;
pub use og_mesh as mesh;
pub use og_offset as offset;
pub use og_select as select;
pub use og_sketch as sketch;
pub use og_topo as topo;

/// The kernel version, as reported by front ends.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
