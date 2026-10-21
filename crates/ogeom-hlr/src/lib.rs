//! Hidden line removal: projecting 3D shapes to annotated 2D edge sets for
//! technical drawings.

mod crossings;
pub mod exact;
pub mod project;
pub mod section;

pub use project::{Drawing, DrawnCurve, Source, View, Visibility, project};
pub use section::{SectionView, broken_section, half_section, hatch, section};
