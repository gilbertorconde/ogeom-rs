//! Blending: constant- and variable-radius edge fillets, vertex blends, chamfers,
//! and their 2D counterparts.
//!
//! *Elsewhere:* `ChFi2d`, `ChFi3d`, `ChFiDS`, `Blend`, `BlendFunc`, `BRepBlend`
//! and `BRepFilletAPI`.
//!
//! Along with offsetting, one of the two areas where every kernel is fragile.
//! Robustness here is a differentiator, not a checkbox.

pub mod analyse;
pub mod chamfer;
pub mod corner;
pub mod corner2d;
pub mod facepair;
pub mod fillet;
pub mod march;
mod marched;
mod support;

pub use analyse::{BlendContact, analyse_blend};
pub use chamfer::{chamfer_edge, chamfer_edge_angle, chamfer_edge_distances};
pub use corner::round_vertex;
pub use corner2d::{chamfer_corner_2d, fillet_corner_2d};
pub use facepair::blend_faces;
pub use fillet::{fillet_edge, fillet_edge_variable, fillet_edges};
pub use march::{
    BlendStop, MarchedBlend, Sides, march_blend, march_blend_seeded, march_blend_sided,
};
