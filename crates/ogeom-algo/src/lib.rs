//! Shape construction and query: everything that does not require surface/surface
//! intersection.
//!
//! The adaptor layer is `ogeom-geom`'s: evaluating topology as geometry
//! belongs with the geometry traits.
//!
//! Every operation in this crate emits history (`generated` / `modified` /
//! `is_deleted`): downstream stable naming is built directly on it.

// Debug output goes through `ogeom_core::diag`, never to stdout or stderr.
#![cfg_attr(
    not(test),
    warn(clippy::print_stderr, clippy::print_stdout, clippy::dbg_macro)
)]

mod bins;
pub mod build;
pub mod certify;
mod chart_cross;
pub mod check;
pub mod classify;
pub mod convert;
pub mod derived;
pub mod distance;
pub mod extend;
mod fill_n;
mod fill_patch;
pub mod frames;
pub mod history;
pub mod hit;
pub mod length;
pub mod mass;
mod mass_chart;
pub mod measure;
pub mod medial;
pub mod medial_graph;
pub mod mesh_solid;
pub mod out_and_back;
pub mod pcurve_fit;
pub mod pcurve_gap;
pub mod place;
mod planar_faces;
pub mod pole;
pub mod primitive;
pub mod project_plane;
pub mod proximity;
pub mod recognize;
mod recognize_patch;
mod recognize_swept;
mod reshape;
mod seam_columns;
mod self_cross;
pub mod sew;
pub mod sheet;
mod sheet_flat;
mod sphere_chart;
pub mod sweep;
pub mod text;
pub mod thickness;
pub mod tight;
pub mod traced;
mod upgrade;
pub mod volume_difference;
mod winding;

pub use build::{
    attach_pcurve, attach_seam, belted_wire, chain_wire_branches, closed_at_poles,
    covers_whole_surface, edge_vertices, edges_turn, find_plane, is_shell_closed, is_wire_closed,
    make_apex_band, make_band_between, make_band_of_rings, make_compound, make_compsolid,
    make_edge, make_edge_between, make_face, make_face_on, make_face_with_pcurves,
    make_natural_face, make_polygon, make_profile_revolution_band, make_revolution_band,
    make_shell, make_solid, make_vertex, make_wire, outer_solid, placed_edge_curve,
    rings_are_parallels, split_ring_on_column, surface_iso_u_curve, trimmed_or_fitted_where_bare,
};
pub use certify::{HalfSpaces, Uncertified, certify_face_within};
pub use chart_cross::{ChartCrossing, chart_crossings};
pub use check::{
    Diagnosis, Problem, Severity, chart_breaks, check, check_kept, check_self_intersection,
    check_self_intersection_near, check_tessellation, inside_out_faces, inside_out_faces_kept,
    restore_containment,
};
pub use classify::{
    Containment, ProbeCache, SolidBoundary, SolidMesh, classify_in_solid, classify_in_solid_exact,
    classify_in_solid_exact_banded, classify_on_face,
};
pub use convert::{
    CurveRestatement, SurfaceRestatement, baked_shape, baked_where_placed,
    general_transformed_shape, normals_oppose, restate_geometry, to_nurbs, to_nurbs_within,
};
pub use extend::{Extension, extend_face};
pub use fill_n::{FillBoundary, FillSide, Filled, make_filling_n};
// The blend engine's side of the filling: held inside a region, or fitted
// once on a given net.
#[doc(hidden)]
pub use fill_n::{
    Region as FillRegion, make_filling_n_fixed, make_filling_n_region, make_filling_n_within,
};
pub use frames::{CurveFrame, CurveFrames, rotation_minimizing_frames};
pub use history::{Built, History, on_own_nodes, on_own_nodes_with, unshare};
pub use hit::{CurveHit, Hit, curve_hit};
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
    FallbackReason, FitConstraints, MeshRegion, MeshRegions, MeshSolid, MeshSolidOptions,
    MeshSolidReport, RegionFallback, RegionId, RegionRefusal, SurfaceKind, refine_solid,
    solid_from_mesh,
};
/// Moved to `ogeom_geom::fit`. These paths are kept until 0.23.0, which
/// removes them.
pub use ogeom_geom::fit::{Spacing, approximate, approximate_within, interpolate};
/// Moved to `ogeom_mesh`. This path is kept until 0.23.0, which removes it.
pub use ogeom_mesh::single_precision_quantum;
pub use out_and_back::{
    Notch, OutAndBack, edge_uses, is_tail, notch_of, out_and_back, placed_edge_uses,
};
pub use pcurve_gap::{
    edge_pcurve_gap, edges_apart, pcurve_gap, state_pcurve_gaps, state_pcurve_gaps_of,
};
pub use place::{copied, transformed};
pub use planar_faces::{PlanarFaces, faces_from_planar_wires};
pub use pole::pole_crossings;
pub use primitive::{
    make_box, make_cone, make_cylinder, make_faceted_solid, make_faceted_solid_with_holes,
    make_half_space, make_hexahedron, make_parallelepiped, make_polyhedron, make_sphere,
    make_torus, make_wedge, torus_face_seamed_at,
};
pub use project_plane::{ProjectedCurve, project_edge_onto_plane};
pub use proximity::{ClosestPair, PointDistance, ShapeDistance, distance_between_shapes};
pub use recognize::{Canonical, Recognized, recognize_points};
pub use reshape::Reshape;
pub use seam_columns::join_seam_columns;
pub use self_cross::{SelfCrossing, self_crossings};
// Reached from the crates above that convert and fit meshes: a wire's and an
// edge's checks, the volume as the faces are flagged and what it keeps, the
// face crossings on kept meshes, distances to a mesh whose normals are
// gathered once, and points welded on a grid.
#[doc(hidden)]
pub use check::{check_wire, edge_spanned_past};
#[doc(hidden)]
pub use distance::{Pseudonormals, to_mesh_with};
#[doc(hidden)]
pub use mass::{VolumeKept, volume_as_flagged};
#[doc(hidden)]
pub use mesh_solid::weld::weld_points;
#[doc(hidden)]
pub use self_cross::crossings as self_crossings_kept;
pub use sew::{
    Agreement, Sewn, agreeing_orientation, make_wire_unordered, order_edges, sew, sew_around,
    sew_around_within, sew_within,
};
pub use sheet::{SheetDistortion, SheetMethod, SheetPlacement, SurfaceSheet};
pub use sphere_chart::{POLE_MARGIN, Recharted, rechart_clear_of_poles};
// The booleans' side of the pole re-chart: the ball operands turned clear
// of the other operand, and how near a pole may come before they are.
#[doc(hidden)]
pub use sphere_chart::{TOUCHING as POLE_TOUCHING, turned as operands_clear_of_poles};
pub use sweep::{make_prism, make_prism_tapered, make_revolution};
pub use text::make_text;
pub use thickness::{Thickness, thickness, thickness_of_shape};
pub use tight::{extent_along, face_bounds, tight_bounds};
pub use upgrade::{
    UnifyOptions, merge_edges, reduce_tolerances, same_domain, unify_same_domain,
    unify_same_domain_around, unify_same_domain_with,
};
// The healing's and the many-operand union's side of the upgrade: edges
// joined at named vertices only, and an edge's measured deviation.
#[doc(hidden)]
pub use upgrade::{measure_edge, merge_edges_at};
pub use volume_difference::{Form, Material, Region, volume_difference};
