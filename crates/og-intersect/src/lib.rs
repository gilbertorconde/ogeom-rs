//! Intersection. Analytic special cases first (plane/plane, plane/cylinder,
//! plane/sphere, cylinder/cylinder, cone/*, torus/*), then a general marching
//! intersector with an approximation stage that fits the resulting polyline to a
//! B-spline.
//!
//! *Elsewhere:* `IntPatch`, `GeomInt`, `IntAna`, `IntCurve`, `IntCurveSurface`,
//! `IntWalk`, `IntSurf`, `ApproxInt` and `Extrema`.
//!
//! **This crate is the project's single largest risk.** No Rust equivalent exists,
//! and it is where prior open-source B-rep efforts have failed. Its quality is
//! gated by a standalone benchmark against analytic ground truth and published
//! datasets before the boolean pipeline in `og-bool` is started.

pub mod benchmark;
pub mod march;
pub mod surface;

pub use benchmark::{Measured, Report, measure, measure_all};
pub use march::{Contact, Marching, Stopped, Traced, branches, seeds, trace};
pub use surface::{Meeting, surface_surface};
