//! Shape construction and query: everything that does not require surface/surface
//! intersection.
//!
//! *Elsewhere:* the `BRepBuilderAPI`, `BRepPrimAPI`, `BRepTools`,
//! `BRepLib`, `BRepGProp`, `BRepCheck`, `BRepClass3d`, `BRepExtrema`, `GeomAPI`,
//! `GCPnts` and `GProp` families, plus the classical curve constructors.
//! The adaptor layer is `ogeom-geom`'s: evaluating topology as geometry
//! belongs with the geometry traits.
//!
//! Every operation in this crate emits history (`generated` / `modified` /
//! `is_deleted`): downstream stable naming is built directly on it.

mod bins;
pub mod build;
pub mod check;
pub mod classify;
pub mod convert;
pub mod extend;
pub mod fit;
pub mod history;
pub mod length;
pub mod mass;
mod mass_chart;
pub mod measure;
pub mod medial;
pub mod medial_graph;
pub mod mesh_solid;
pub mod pcurve_fit;
pub mod place;
pub mod primitive;
pub mod project_plane;
pub mod proximity;
pub mod recognize;
mod recognize_swept;
pub mod sew;
pub mod sweep;
pub mod text;
pub mod tight;

pub use build::{
    attach_pcurve, attach_seam, chain_wire_branches, closed_at_poles, edge_vertices, find_plane,
    is_shell_closed, is_wire_closed, make_apex_band, make_band_between, make_band_of_rings,
    make_compound, make_compsolid, make_edge, make_edge_between, make_face, make_face_on,
    make_face_with_pcurves, make_natural_face, make_polygon, make_revolution_band, make_shell,
    make_solid, make_vertex, make_wire, rings_are_parallels, surface_iso_u_curve,
};
pub use check::{
    Diagnosis, Problem, Severity, check, check_self_intersection, check_self_intersection_near,
    check_tessellation, inside_out_faces, restore_containment,
};
pub use classify::{
    Containment, SolidBoundary, SolidMesh, classify_in_solid, classify_in_solid_exact,
    classify_in_solid_exact_banded, classify_on_face,
};
pub use convert::{
    CurveRestatement, SurfaceRestatement, baked_shape, general_transformed_shape, normals_oppose,
    restate_geometry, to_nurbs, to_nurbs_within,
};
pub use extend::{Extension, extend_face};
pub use fit::{Spacing, approximate, approximate_within, interpolate};
pub use history::{Built, History};
pub use length::{curve_length, parameter_at_length, points_by_count, points_by_spacing};
pub use mass::{MassProperties, linear_properties, surface_properties, volume_properties};
pub use measure::{
    Obb, Projection, SurfaceProjection, SurfaceSeeds, curve_bounds, face_normal, oriented_bounds,
    project_on_curve, project_on_planar_curve, project_on_surface, project_on_surface_from,
    relative_deflection, shape_bounds, surface_bounds, vertex_bounds, widened_to_hold,
};
pub use medial::{MedialAxis, medial_axis};
pub use medial_graph::{MedialBranch, MedialGraph, MedialSite, MedialVertex, medial_graph};
pub use mesh_solid::{
    MeshSolid, MeshSolidOptions, MeshSolidReport, refine_solid, single_precision_quantum,
    solid_from_mesh,
};
pub use place::{copied, transformed};
pub use primitive::{
    make_box, make_cone, make_cylinder, make_half_space, make_hexahedron, make_parallelepiped,
    make_polyhedron, make_sphere, make_torus, make_wedge,
};
pub use project_plane::{ProjectedCurve, project_edge_onto_plane};
pub use proximity::{ClosestPair, ShapeDistance, distance_between_shapes};
pub use recognize::{Canonical, Recognized, recognize_points};
pub use sew::{Sewn, make_wire_unordered, order_edges, sew, sew_within};
pub use sweep::{make_prism, make_prism_tapered, make_revolution};
pub use text::make_text;
pub use tight::tight_bounds;
