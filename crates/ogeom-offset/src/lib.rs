//! Offsetting and sweeping: offset shape, thicken/shell, pipe, pipe-shell,
//! loft/thru-sections, draft angle, evolved, filling, normal projection.

// Debug output goes through `ogeom_core::diag`, never to stdout or stderr.
#![cfg_attr(
    not(test),
    warn(clippy::print_stderr, clippy::print_stdout, clippy::dbg_macro)
)]

pub mod draft;
pub mod feature;
pub mod fill;
pub(crate) mod fixed;
mod freeform;
pub mod loft;
pub mod middle;
pub mod open_chain;
mod prism_until;
pub mod project;
pub mod shape;
pub mod sheet;
pub mod sweep;
pub mod thick;
pub mod wire2d;

pub use draft::apply_draft;
pub use feature::{Feature, feature_prism, feature_revol, feature_rib, feature_slot};
pub use fill::{FillBoundary, FillSide, Filled, make_filling, make_filling_n};
pub use loft::make_loft_solid;
pub use middle::{MiddlePath, middle_path};
pub use open_chain::{Cap, make_rib, thicken_open_wire};
pub use prism_until::make_prism_until;
pub use project::{Projected, normal_projection};
pub use shape::{
    Walls, make_thick_solid, make_thick_solid_sided, make_thick_solid_with, move_faces,
    offset_faces, offset_shape,
};
pub use sheet::{make_loft_surface, make_ruled, make_sweep_surface, make_sweep_two_rails};
pub use sweep::{
    PipeContact, PipeCorners, PipeEvolution, PipeLaw, PipeLaws, PipePlacement, SpineFunction,
    SpineLaw, make_evolved, make_helical_sweep, make_loft, make_loft_skinned,
    make_loft_skinned_aligned, make_loft_skinned_closed, make_pipe, make_pipe_sections,
    make_pipe_shell, make_pipe_shell_evolving, make_pipe_shell_law, make_pipe_shell_laws,
    make_pipe_shell_with, make_pipe_skinned, make_revolution_until, place_on_spine,
};
pub use thick::{make_thick_sheet, make_thick_sheet_sided, offset_sheet};
pub use wire2d::{Join, offset_wire};
