//! The corpus suites of the `ogeom` facade: the tests that read `tests/corpus/`
//! at the repository root, which the published package leaves out.

#[path = "../support/walks.rs"]
#[allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
mod walks;

mod assembly_corpus;
mod exchange;
mod exchange_tolerances;
mod foundations;
mod imported_corpus;
mod m5_close;
mod oblique_drills;
mod orient_corpus;
mod pmi_corpus;
mod seamless_faces;
mod wound_rims;
