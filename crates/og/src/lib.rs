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

/// The kernel version, as reported by front ends.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
