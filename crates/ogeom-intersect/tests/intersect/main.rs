//! The integration tests of `ogeom-intersect`, in one test program so the
//! crate links once.

#[path = "../support/mod.rs"]
#[allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
mod support;

mod branch_points;
mod bspline_boundary_pcurve;
mod hard;
mod instruments;
mod plane_cone;
mod ssi_corpus;
mod stress;
