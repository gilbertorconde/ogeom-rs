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
//!
//! # Two formats, two things
//!
//! [`native`] carries a whole document — topology, geometry, placements,
//! tolerances, provenance — and reads back as the same model, handles and all.
//! [`stl`] carries a triangle soup and loses everything else, which is the
//! format's nature rather than a shortcoming of the writer.
//!
//! The crate root re-exports STL's [`read`] and [`write()`] for the common case;
//! the native pair are reached as [`native::read`] and
//! [`native::write()`], because "write this shape" means something different for
//! each and a name that did not say which would be the wrong kind of convenience.

pub mod native;
pub mod step;
pub mod stl;

pub use step::{StepImport, StepReport, read_step};
pub use stl::{Encoding, read, write};
