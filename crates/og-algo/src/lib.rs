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
pub mod check;
pub mod classify;
pub mod fit;
pub mod history;
pub mod length;
pub mod mass;
pub mod measure;
pub mod place;
pub mod primitive;
pub mod proximity;
pub mod reconstruct;
pub mod sew;
pub mod sweep;

pub use build::{
    attach_pcurve, attach_seam, edge_vertices, find_plane, is_shell_closed, is_wire_closed,
    make_edge, make_edge_between, make_face, make_face_on, make_face_with_pcurves,
    make_natural_face, make_polygon, make_revolution_band, make_shell, make_solid, make_vertex,
    make_wire, surface_iso_u_curve,
};
pub use check::{Diagnosis, Problem, Severity, check, check_tessellation};
pub use classify::{Containment, classify_in_solid, classify_in_solid_exact, classify_on_face};
pub use fit::{Spacing, approximate, approximate_within, interpolate};
pub use history::{Built, History};
pub use length::{curve_length, parameter_at_length, points_by_count, points_by_spacing};
pub use mass::{MassProperties, linear_properties, surface_properties, volume_properties};
pub use measure::{
    Obb, Projection, SurfaceProjection, curve_bounds, oriented_bounds, project_on_curve,
    project_on_planar_curve, project_on_surface, relative_deflection, shape_bounds, surface_bounds,
    vertex_bounds,
};
pub use place::{copied, transformed};
pub use primitive::{
    make_box, make_cone, make_cylinder, make_half_space, make_sphere, make_torus, make_wedge,
};
pub use proximity::{ClosestPair, ShapeDistance, distance_between_shapes};
pub use reconstruct::{CREASE, Region, planar_regions, to_brep};
pub use sew::{Sewn, make_wire_unordered, order_edges, sew};
pub use sweep::{make_prism, make_revolution};
