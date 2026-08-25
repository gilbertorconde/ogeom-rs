//! Foundation layer: everything else in the workspace is shaped by this crate.
//!
//! *Elsewhere:* the `Standard` package, `Precision`, and the parts of a
//! collections library that carry semantics rather than just storage.
//!
//! Four responsibilities:
//!
//! - **Arenas.** [`Arena`] and [`Key`] — typed generational index arenas replace
//!   the intrusive reference counting conventional kernels use. No cycle hazard,
//!   cache-friendly, and — critically — arena keys are what make stable entity
//!   identity possible at all.
//! - **Identity.** [`EntityId`] plus [`Provenance`]. Conventionally a topology
//!   node's identity is its address, so every modeling operation produces new ones
//!   and downstream references break; that *is* the topological naming problem. We
//!   record where an entity came from at creation time instead of reconstructing it
//!   from history maps afterwards.
//! - **Errors.** [`OgeomError`] covers the failure vocabulary a kernel needs, lining
//!   up with the categories applications already handle, but as a plain [`Result`]
//!   — no exceptions, and no conversion of hardware signals into throwable
//!   objects.
//! - **Numerics.** [`Tolerances`], [`Tolerance`] and the [`Predicates`] trait, so
//!   the robustness strategy can change without rewriting a single algorithm.
//!
//! See `docs/DATA_MODEL.md` for the normative invariants. Section numbers in this
//! crate's documentation refer to it.

pub mod arena;
pub mod error;
pub mod id;
pub mod parallel;
pub mod predicates;
pub mod progress;
pub mod tolerance;

pub use arena::{Arena, Key, UNSCOPED};
pub use error::{Cause, OgeomError, OgeomResult};
pub use id::{EntityId, OpId, Provenance, ProvenanceTable, Role, SourceId};
pub use predicates::{Exact, Fast, P2, P3, Predicates, Sign};
pub use progress::{Canceller, Stage, Watch};
pub use tolerance::{Tolerance, Tolerances, check_containment};

/// The predicate implementation algorithms use unless told otherwise.
///
/// Exact by default: a wrong combinatorial decision in a triangulation or a
/// boolean is not a slightly-off number, it is an inconsistent structure that
/// fails somewhere else entirely. Opt into [`Fast`] with a measurement in hand.
pub type DefaultPredicates = Exact;
