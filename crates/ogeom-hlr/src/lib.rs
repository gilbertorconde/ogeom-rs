//! Hidden line removal: projecting 3D shapes to annotated 2D edge sets for
//! technical drawings.

// Debug output goes through `ogeom_core::diag`, never to stdout or stderr.
#![cfg_attr(
    not(test),
    warn(clippy::print_stderr, clippy::print_stdout, clippy::dbg_macro)
)]

mod crossings;
pub mod exact;
pub mod project;
pub mod section;

pub use project::{Drawing, DrawnCurve, Source, View, Visibility, project};
pub use section::{SectionView, broken_section, half_section, hatch, section};
