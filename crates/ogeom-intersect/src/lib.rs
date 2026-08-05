//! Intersection: where curves and surfaces meet.
//!
//! *Elsewhere:* `IntPatch`, `GeomInt`, `IntAna`, `IntCurve`, `IntCurveSurface`,
//! `IntWalk`, `IntSurf`, `ApproxInt` and `Extrema`.
//!
//! The entry points, by what is being intersected:
//!
//! - [`intersect_surfaces`] — surface/surface, the one `ogeom-bool` builds on.
//!   Analytic closed forms where they exist, with exact same-parameter
//!   pcurves; the marching-and-fitting pipeline where they do not, with the
//!   tolerance reported as a sum of stated parts.
//! - [`intersect_curves_2d`] / [`intersect_curves`] — curve/curve in the
//!   plane and in space. The planar one is what boolean face splitting runs
//!   on; the spatial one reports the gap each crossing achieved, because
//!   space curves generically miss.
//! - [`intersect_curve_surface`] — curve/surface, the well-posed system, and
//!   what the exact point-in-solid classifier will cast its rays with.
//!
//! The stages underneath are public because each is separately measurable:
//! [`surface_surface`] (the closed forms), [`seeds`]/[`trace`]/[`branches`]
//! (the marcher), [`approximate_branch`] (polyline to curves).
//!
//! # The instruments
//!
//! **This crate was the project's single largest risk** — no Rust equivalent
//! existed, and it is where prior open-source B-rep efforts failed — so it
//! carries its own instruments rather than trusting itself. [`measure_all`]
//! scores accuracy: every point of every reported curve against both surfaces,
//! with ground truth being the surfaces themselves. [`coverage()`] scores
//! completeness: the surfaces are asked, by signed distance and the
//! intermediate value theorem, where the intersection *must* be, and every
//! such cell had better be reached by some branch — because an intersector
//! that finds one circle of two scores perfectly on accuracy alone. Both
//! instruments have negative controls in the test suite: they demonstrably
//! fail when something is genuinely missing.
//!
//! What the gate still lacks is an input no one here generated: a published
//! corpus. Until that lands, the numbers say the machinery is sound on what it
//! has seen, and the deferred table says so in as many words.

pub mod approx;
pub mod benchmark;
pub mod coverage;
pub mod curve_surface;
pub mod curves;
pub mod extrema;
pub mod march;
pub mod section;
pub mod surface;

pub use approx::{IntersectionCurve, approximate_branch};
pub use benchmark::{Measured, Report, measure, measure_all};
pub use coverage::{Coverage, coverage};
pub use curve_surface::{
    CurveSurfaceIntersection, CurveSurfaceOptions, Piercing, intersect_curve_surface,
};
pub use curves::{
    Crossing, CurveCurveOptions, CurveIntersection, Overlap, intersect_curves, intersect_curves_2d,
};
pub use extrema::{
    Approach, Extrema, ExtremaOptions, extrema_curve_curve, extrema_curve_surface,
    extrema_surface_surface,
};
pub use march::{Contact, Marching, Stopped, Traced, branches, seeds, trace};
pub use section::{
    IntersectOptions, SectionCurve, SurfaceIntersection, exact_pcurve_of, exact_pcurve_over,
    intersect_surfaces,
};
pub use surface::{Meeting, surface_surface};
