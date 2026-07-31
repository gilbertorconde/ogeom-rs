//! Shape construction and query — everything that does not require surface/surface
//! intersection.
//!
//! *Elsewhere:* the `BRepBuilderAPI`, `BRepPrimAPI`, `BRepAdaptor`, `BRepTools`,
//! `BRepLib`, `BRepGProp`, `BRepCheck`, `BRepClass3d`, `BRepExtrema`, `GeomAPI`,
//! `GCPnts` and `GProp` families, plus the classical curve constructors.
//!
//! Every operation in this crate emits history (`generated` / `modified` /
//! `is_deleted`) from the start. That is not optional and cannot be retrofitted:
//! downstream stable naming is built directly on it.

pub mod build;
pub mod history;
pub mod primitive;

pub use build::{
    attach_pcurve, edge_vertices, is_shell_closed, is_wire_closed, make_edge, make_edge_between,
    make_face, make_natural_face, make_shell, make_solid, make_vertex, make_wire,
};
pub use history::{Built, History};
pub use primitive::make_box;
