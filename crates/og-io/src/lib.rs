//! Data exchange.
//!
//! *Elsewhere:* the shape serialization, STEP, IGES, glTF and transfer-framework
//! packages.
//!
//! The native format comes first, in Phase 1, because it is what the test suite
//! round-trips through — not because anyone wants it.
//!
//! STEP is the one with real external value: there is no production-grade Rust STEP
//! importer today. Parsing is the easy 20% (`ruststep`/`espr` cover it); the other
//! 80% is semantic mapping onto our topology, unit and assembly-transform handling,
//! and surviving spec-violating output from commercial CAD systems.

pub mod stl;

pub use stl::{Encoding, read, write};
