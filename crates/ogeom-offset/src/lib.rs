//! Offsetting and sweeping: offset shape, thicken/shell, pipe, pipe-shell,
//! loft/thru-sections, draft angle, evolved, filling, normal projection.
//!
//! *Elsewhere:* `BRepOffset`, `BRepOffsetAPI`, `BRepFill`, `BRepSweep`,
//! `GeomFill`, `GeomPlate` and `BRepFeat`.

pub mod draft;
pub mod feature;
pub mod fill;
pub(crate) mod fixed;
pub mod middle;
pub mod project;
pub mod shape;
pub mod sweep;
pub mod wire2d;

pub use draft::apply_draft;
pub use feature::{Feature, feature_prism, feature_revol, feature_rib, feature_slot};
pub use fill::make_filling;
pub use middle::{MiddlePath, middle_path};
pub use project::{Projected, normal_projection};
pub use shape::{make_thick_solid, make_thick_solid_with, move_faces, offset_faces, offset_shape};
pub use sweep::{
    PipeCorners, PipeLaw, make_evolved, make_helical_sweep, make_loft, make_loft_skinned,
    make_loft_skinned_aligned, make_loft_skinned_closed, make_pipe, make_pipe_sections,
    make_pipe_shell, make_pipe_shell_law, make_pipe_shell_with, make_pipe_skinned,
    make_revolution_until,
};
pub use wire2d::{Join, offset_wire};
