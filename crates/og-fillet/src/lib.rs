//! Blending: constant- and variable-radius edge fillets, vertex blends, chamfers,
//! and their 2D counterparts.
//!
//! *Elsewhere:* `ChFi2d`, `ChFi3d`, `ChFiDS`, `Blend`, `BlendFunc`, `BRepBlend`
//! and `BRepFilletAPI`.
//!
//! Along with offsetting, one of the two areas where every kernel is fragile.
//! Robustness here is a differentiator, not a checkbox.

pub mod chamfer;
pub mod fillet;
mod support;

pub use chamfer::chamfer_edge;
pub use fillet::fillet_edge;
