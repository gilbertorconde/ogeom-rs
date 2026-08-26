//! Hidden line removal — projecting 3D shapes to annotated 2D edge sets for
//! technical drawings.
//!
//! *Elsewhere:* `HLRAlgo`, `HLRBRep`, `HLRTopoBRep` and `HLRAppli`.
//!
//! Nothing in the modeling path depends on it, but producing a 2D drawing from a
//! 3D model is a core CAD capability and not optional for a complete kernel.

pub mod exact;
pub mod project;
pub mod section;

pub use project::{Drawing, DrawnCurve, Source, View, Visibility, project};
pub use section::{SectionView, broken_section, half_section, hatch, section};
