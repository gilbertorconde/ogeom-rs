//! Offsetting and sweeping: offset shape, thicken/shell, pipe, pipe-shell,
//! loft/thru-sections, draft angle, evolved, filling, normal projection.
//!
//! *Elsewhere:* `BRepOffset`, `BRepOffsetAPI`, `BRepFill`, `BRepSweep`,
//! `GeomFill`, `GeomPlate` and `BRepFeat`.
//!
//! The other perennially fragile area. Robustness here is a differentiator.

pub mod draft;
pub mod feature;
pub mod fill;
pub mod project;
pub mod shape;
pub mod sweep;
pub mod wire2d;

pub use draft::apply_draft;
pub use feature::{Feature, feature_prism, feature_revol, feature_rib, feature_slot};
pub use fill::make_filling;
pub use project::{Projected, normal_projection};
pub use shape::{make_thick_solid, offset_shape};
pub use sweep::{
    make_evolved, make_loft, make_loft_skinned, make_loft_skinned_closed, make_pipe,
    make_pipe_shell, make_pipe_skinned,
};
pub use wire2d::{Join, offset_wire};
