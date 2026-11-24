//! ogeom: a boundary-representation CAD kernel.
//!
//! This umbrella crate is the public API. Depend on this, not on the `ogeom-*`
//! crates individually. Their boundaries are an implementation detail and will
//! move.
//!
//! The modules below, and the items each one names, are that API. Every module
//! lists its items one by one, so an `ogeom::<module>::<item>` path stays where
//! it is when the item moves to another crate underneath. Submodules such as
//! [`algo::check`](mod@algo::check) or [`io::step`](mod@io::step) are re-exported whole.
//!
//! By task:
//!
//! | Task | Modules |
//! |---|---|
//! | Modelling | [`topo`] (the data model), [`geom`], [`math`], [`algo`] (primitives, sweeps, sewing), [`offset`] |
//! | Booleans | [`boolean`], [`intersect`] |
//! | Blends | [`fillet`] |
//! | Exchange | [`io`], [`doc`] |
//! | Mesh | [`mesh`] (tessellation), [`algo::mesh_solid`] (mesh to solid) |
//! | Analysis | [`algo`] (measurement, mass properties, classification, checks), [`hlr`] |
//! | Repair | [`heal`] |
//! | Foundation | [`core`] (errors, tolerances) |
//!
//! A from-scratch Rust implementation: no binding, no wrapper, no dependency on
//! any existing CAD kernel. It keeps the *data model* the field converged on (
//! shared topology with location chains, per-entity tolerances,
//! multi-representation edges, operation history) because that model is what
//! makes the algorithms possible. See `docs/DATA_MODEL.md`.

// Debug output goes through `ogeom_core::diag`, never to stdout or stderr.
#![cfg_attr(
    not(test),
    warn(clippy::print_stderr, clippy::print_stdout, clippy::dbg_macro)
)]

/// Construction, measurement, classification, sewing and mesh to solid.
pub mod algo {
    pub use ogeom_algo::{
        Across, Agreement, Built, Canonical, ChainedPatch, ChartCrossing, ClosestPair, Containment,
        CurveDef, CurveFrame, CurveFrames, CurveHit, CurveRestatement, Cutoff, Diagnosis,
        Extension, FallbackReason, FillBoundary, FillSide, Filled, FitConstraints, Form,
        HalfSpaces, History, Hit, IrregularCorner, LoopSide, MassProperties, Material, MedialAxis,
        MedialBranch, MedialGraph, MedialSite, MedialVertex, MeshRegion, MeshRegions, MeshSolid,
        MeshSolidOptions, MeshSolidReport, NetOutside, Network, NetworkOptions, Notch, Obb,
        OutAndBack, Outline, POLE_MARGIN, PatchDef, PatchFit, PatchSide, Piece2, PlanarFaces,
        PointDistance, ProbeCache, Problem, ProjectedCurve, Projection, Recharted, Recognized,
        Region, RegionFallback, RegionId, RegionRefusal, Reshape, Rest, Segment2, Segments,
        SelfCrossing, Severity, Sewn, ShapeDistance, SheetDistortion, SheetMethod, SheetPlacement,
        SolidBoundary, SolidMesh, SurfaceKind, SurfaceProjection, SurfaceRestatement, SurfaceSeeds,
        SurfaceSheet, Thickness, Uncertified, UnifyOptions, Wire2d, agreeing_orientation,
        attach_pcurve, attach_seam, baked_shape, baked_where_placed, belted_wire, build, certify,
        certify_face_within, chain_wire_branches, chart_breaks, chart_crossings, check, check_kept,
        check_self_intersection, check_self_intersection_near, check_tessellation, classify,
        classify_in_solid, classify_in_solid_exact, classify_in_solid_exact_banded,
        classify_on_face, closed_at_poles, convert, copied, covers_whole_surface, curve_bounds,
        curve_hit, curve_length, derived, distance, distance_between_shapes, edge_pcurve_gap,
        edge_uses, edge_vertices, edges_apart, edges_turn, extend, extend_face, extent_along,
        face_bounds, face_normal, faces_from_planar_wires, find_plane, fit, fit_loop_quads,
        fit_network, fit_network_chained, fit_outline, fit_segments_2d, frames,
        general_transformed_shape, held, history, hit, inside_out_faces, inside_out_faces_kept,
        is_shell_closed, is_tail, is_wire_closed, join_seam_columns, length, linear_properties,
        make_apex_band, make_band_between, make_band_of_rings, make_box, make_compound,
        make_compsolid, make_cone, make_cylinder, make_edge, make_edge_between, make_face,
        make_face_on, make_face_with_pcurves, make_faceted_solid, make_faceted_solid_with_holes,
        make_filling_n, make_half_space, make_hexahedron, make_natural_face, make_parallelepiped,
        make_polygon, make_polyhedron, make_prism, make_prism_tapered,
        make_profile_revolution_band, make_revolution, make_revolution_band, make_shell,
        make_solid, make_sphere, make_text, make_torus, make_vertex, make_wedge, make_wire,
        make_wire_unordered, mass, measure, medial, medial_axis, medial_graph, merge_edges,
        mesh_solid, net_held_within, network, normals_oppose, notch_of, on_own_nodes,
        on_own_nodes_with, order_edges, orient, oriented_bounds, out_and_back, outer_solid,
        outline, parameter_at_length, pcurve_fit, pcurve_gap, place, placed_edge_curve,
        placed_edge_uses, points_by_count, points_by_spacing, pole, pole_crossings, primitive,
        project_edge_onto_plane, project_on_curve, project_on_planar_curve, project_on_surface,
        project_on_surface_from, project_plane, proximity, rechart_clear_of_poles, recognize,
        recognize_points, reduce_tolerances, refine_solid, register, relative_deflection, rest,
        restate_geometry, resting_faces, resting_faces_of, restore_containment,
        rings_are_parallels, rotation_minimizing_frames, same_domain, segments_2d, self_crossings,
        sew, sew_around, sew_around_within, sew_within, shape_bounds, sheet, solid_from_mesh,
        split_ring_on_column, state_pcurve_gaps, state_pcurve_gaps_of, surface_bounds,
        surface_iso_u_curve, surface_properties, sweep, text, thickness, thickness_of_shape, tight,
        tight_bounds, to_nurbs, to_nurbs_within, torus_face_seamed_at, traced, transformed,
        trimmed_or_fitted_where_bare, unify_same_domain, unify_same_domain_around,
        unify_same_domain_with, unshare, vertex_bounds, volume_difference, volume_properties,
        widened_to_hold,
    };
    pub use ogeom_geom::fit::{Spacing, approximate, approximate_within, interpolate};
    pub use ogeom_mesh::single_precision_quantum;
}

/// Booleans, sections, splits and face removal.
pub mod boolean {
    pub use ogeom_bool::{
        Cells, POLE_MARGIN, POLE_TOUCHING, Recharted, cells, common, common_fuzzy, cut,
        cut_by_mesh_region, cut_fuzzy, cut_seated, fuse, fuse_all, fuse_fuzzy, fuse_seated,
        make_periodic, make_volume, note_meeting_along_edge, rechart_clear_of_poles, remove_faces,
        section, section_face, self_union, split_sheet,
    };
}

/// Arenas, identity, errors, tolerances, predicates and diagnostics.
pub mod core {
    pub use ogeom_core::{
        Arena, Canceller, Cause, DefaultPredicates, EntityId, Exact, Fast, FastHasher, FastMap,
        FastSet, Key, OgeomError, OgeomResult, OpId, P2, P3, Predicates, Provenance,
        ProvenanceTable, Role, Sign, SourceId, Stage, Tolerance, Tolerances, UNSCOPED, Watch,
        arena, check_containment, clock, collections, diag, error, id, ogeom_bail, ogeom_err,
        parallel, predicates, progress, tolerance,
    };
}

/// Documents: assemblies, instances, colours, names and PMI.
pub mod doc {
    pub use ogeom_doc::{
        Annotated, Callout, Colour, Construction, ConstructionGeometry, Datum, DatumModifier,
        DatumModifierKind, DatumReference, DatumTarget, DatumTargetKind, Dimension, Document,
        GeometricTolerance, Instance, LEAST_MATERIAL, Layer, LayerId, MAXIMUM_MATERIAL, Material,
        MaterialId, MeasureKind, Note, Occurrence, Pmi, Product, ProductId, ProductKind, Property,
        PropertyValue, Texture, TextureId, TextureMapping, ValidationProperties, View, attributes,
        construction, pmi, structure, view,
    };
}

/// Fillets, chamfers and vertex blends.
pub mod fillet {
    #[allow(
        deprecated,
        reason = "re-exported for the release they are deprecated in"
    )]
    pub use ogeom_fillet::{
        Blend, BlendContact, BlendLimit, BlendRefusal, BlendRequest, BlendSection, BlendSizing,
        BlendStop, ChainRules, Chamfer, ChamferCorner, DistanceLawSize, End, FaceCurvature,
        FillBoundary, FillSide, Filled, FilletCorner, LIMIT_PROBES, LIMIT_SAMPLES, LIMIT_TOLERANCE,
        LawSize, LawStation, MarchedBlend, Overflow, RadiusLaw, RoundProfile, Setbacks, Sides,
        TANGENT_CHAIN_ANGLE, Vary, analyse, analyse_blend, blend_edges, blend_edges_overflow,
        blend_faces, blend_limit, blend_seat, bridge, chamfer, chamfer_corner_2d, chamfer_edge,
        chamfer_edge_angle, chamfer_edge_distances, chamfer_edges, chamfer_edges_checked,
        chamfer_edges_cornered, chamfer_edges_overflow, chamfer_edges_overflow_checked,
        chamfer_edges_partial, chamfer_edges_with, corner2d, face_curvature_at,
        face_curvature_samples, facepair, fillet, fillet_corner_2d, fillet_edge,
        fillet_edge_distances, fillet_edge_distances_overflow, fillet_edge_distances_variable,
        fillet_edge_variable, fillet_edges, fillet_edges_checked, fillet_edges_corner,
        fillet_edges_corner_checked, fillet_edges_corner_overflow,
        fillet_edges_corner_overflow_checked, fillet_edges_cornered,
        fillet_edges_distances_variable, fillet_edges_overflow, fillet_edges_partial,
        fillet_edges_variable, fillet_edges_variable_overflow, fillet_edges_with, fillet_faces,
        fillet_sheet_edges, largest_chamfer_distance, largest_chamfer_distance_all,
        largest_chamfer_distance_all_overflow, largest_chamfer_distance_overflow,
        largest_fillet_distances, largest_fillet_radius, largest_fillet_radius_all,
        largest_fillet_radius_overflow, largest_law_distances, largest_law_radius,
        largest_law_scale, make_blend_curve, make_blend_surface, make_filling_n, march,
        march_blend, march_blend_seeded, march_blend_sided, round_vertex, round_vertex_setback,
        round_vertex_setback_sharp, sheet, tangent_chain, tangent_chain_with,
    };
}

/// Curves and surfaces.
pub mod geom {
    pub use ogeom_geom::{
        BSpline2d, BSplineCurve, BSplineSurface, Circle2d, CircleCurve, ConeSurface, Continuity,
        Curve, Curve2d, Curve3d, CurveKind, CurveOnSurface, CylinderSurface, Ellipse2d,
        EllipseCurve, ExtrusionSurface, HelixCurve, Hyperbola2d, HyperbolaCurve, Line2d, LineCurve,
        Offset2d, OffsetCurve, OffsetSurface, ParabolaCurve, PlanarCurve, PlaneSurface, Reversible,
        RevolutionSurface, SphereSurface, Surface, SurfaceCurvature, SurfaceGeometry, SurfaceJet,
        SurfaceKind, SurfaceProjection, SurfaceSeeds, TorusSurface, Transformable, Trig2d,
        Trimmed2d, TrimmedCurve, TrimmedSurface, carried, convert, curve, curve2d, fit, project,
        project_on_surface, project_on_surface_from, surface, tangent_angle, traits,
    };
}

/// Healing, simplification and upgrading of imported shapes.
pub mod heal {
    pub use ogeom_heal::{
        CanonicalReport, FixReport, Fixed, FixedTrims, FreeBounds, IsoLine, PieceOrigin,
        Projection, ReanchoredBoundaries, Reshape, SameParameterReport, Simplified, SmallFaces,
        Split, SplitOptions, UnifyOptions, canonical, canonical_simplify, divide, divide_by_angle,
        divide_by_area, divide_by_continuity, divide_face, fix, fix_face_pcurves, fix_shape,
        fix_small_faces, free_bounds, join_seam_columns, merge_edges, reanchor,
        reanchor_boundaries, reanchor_periodic_rings, recognize_surface, reduce_tolerances,
        remodel, remove_small_solids, repair_same_parameter, reshape, restrict_degree, same_domain,
        same_parameter, seam, seam_periodic_faces, small, solid, solids_from_shells, split,
        split_face, split_face_with, swept_to_elementary, to_bezier, unify_same_domain,
        unify_same_domain_around, unify_same_domain_with, upgrade,
    };
}

/// Hidden line removal and drawings.
pub mod hlr {
    pub use ogeom_hlr::{
        Drawing, DrawnCurve, SectionView, Source, View, Visibility, broken_section, exact,
        half_section, hatch, project, section,
    };
}

/// Curve and surface intersection.
pub mod intersect {
    pub use ogeom_intersect::{
        Approach, Condition, Contact, Crossing, CurveCurveOptions, CurveIntersection,
        CurveSurfaceIntersection, CurveSurfaceOptions, Extrema, ExtremaOptions, IntersectOptions,
        IntersectionCurve, Marching, Meeting, Overlap, ParameterWindow, Piercing, PreparedSurface,
        SectionCurve, Stopped, SurfaceIntersection, Traced, Walked, approx, approximate_branch,
        branches, branches_cut_at_poles, contact, crosses_extrusion_rim, curve_surface, curves,
        cut_at_poles, drum_onto_plane, exact_pcurve_of, exact_pcurve_over, extrema,
        extrema_curve_curve, extrema_curve_surface, extrema_surface_surface, follow,
        intersect_curve_surface, intersect_curves, intersect_curves_2d, intersect_surfaces,
        intersect_surfaces_within, march, passes_a_pole, plane_extrusion_at_rims, section, seeds,
        surface, surface_surface, touch_curve_surface, trace, trace_tangential, walk, walk_one_way,
    };
}

/// File formats: STEP, IGES, mesh formats, DXF and the native format.
pub mod io {
    pub use ogeom_io::{
        Encoding, ExportMesh, IgesExport, IgesImport, IgesReport, ImportedMesh, MeshColours,
        ObjImport, ObjectRole, ObjectType, PlyImport, PointCloud, Scalar, StepImport, StepReport,
        StlMesh, ThreeMfImport, ThreeMfObject, ThreeMfVolume, UntrimmedFace, WarningSummary,
        XyzColumn, brep, cloud, dxf, iges, json, mesh_formats, native, read, read_3mf, read_glb,
        read_gltf, read_iges, read_iges_file, read_obj, read_obj_file, read_obj_with_materials,
        read_ply, read_ply_bytes, read_pts, read_step, read_vrml, read_with_quantum, read_xyz,
        step, stl, threemf, vrml, write, write_3mf, write_glb, write_iges, write_iges_reported,
        write_obj, write_ply, write_ply_binary, write_step, write_xyz,
    };
}

/// Points, vectors, transforms, B-spline basis and solvers.
pub mod math {
    pub use ogeom_math::{
        Aabb, Axis, Axis2, BasisValues, Bisector2, Blend, Circle, Circle2, Cone, ControlGrid,
        Convergence, Criteria, CurvePoint, Cylinder, Direction, Direction2, Ellipse, Ellipse2,
        Frame, Frame2, GeneralTransform, Handedness, Hyperbola, Hyperbola2, Interval, KnotVector,
        Matrix2, Matrix3, Minimum, Parabola, Parabola2, Placement, Plane, Point, Point2,
        Quaternion, Solution, SparseMatrix, Sphere, SurfacePoint, SystemSolution, TangentCircle,
        Target2, Torus, TorusKind, Transform, Transform2, TransformKind, Vector, Vector2, Weighted,
        bisector, bounds, bspline, circles_of_radius_tangent_to_two, circles_tangent_to_three,
        conic, construct2d, direction, elementary, frame, gauss_kronrod, gauss_kronrod_21_rule,
        gauss_legendre, gauss_legendre_rule, global_minimum, integrate, integrate_2d, interval,
        knots, least_squares_cgnr, lines_tangent_to_two_circles, matrix, minimize_local, optimize,
        point, quadric, quaternion, solve, sparse, swarm_minimum, transform, vector,
    };
}

/// Tessellation and triangle-mesh operations.
pub mod mesh {
    pub use ogeom_mesh::{
        BoolKind, Bvh, Cap, Chain, Component, Curvature, CutReport, Deflection, Deviation, Edge2,
        EdgeChords, Editable, FaceMesh, FaceMeshCache, FaceMeshCacheStats, Fairing, Fill,
        FillReport, Grid, Hit, Hollowed, Lies, Line, Loop2, MeshBoolean, NEIGHBOURS, Neighbour,
        OFFSET_EDGE_ANGLE, OFFSET_FARTHEST, OPTIMAL_UP_TO, OffsetReport, PASSES, PointIndex,
        Polyline, Quantum, RIM_CORNER, Remesh, RemeshReport, Resolved, Rim, RimId, Section, Side,
        Simplified, SimplifyReport, Smooth, SmoothReport, SubdivideReport, Subdivision, Target,
        Tessellated, Thicken, Thickened, WeldedMesh, attach, boolean, boolean_with_shape, bvh,
        curvature, discretize, discretize_on_surface, discretize_planar, edge_chords_for, editable,
        face_boundary, face_boundary_in_chart, face_boundary_in_chart_met, face_meshes_for, hatch,
        hatch_face, hollow, inside_boundary, inside_boundary_with, is_straight, is_straight_planar,
        open_chart_ring, parity, planar, point_index, polyline_of, polyline_of_edge, region,
        resolve_self_intersections, resolve_self_intersections_with, section, sections,
        self_intersections, shadow, simplify, single_precision_quantum, tessellate, thicken,
        triangulate, triangulate_face, triangulate_face_kept, triangulate_face_with,
        triangulate_with_chords, triangulate_with_face_meshes, triangulate_with_faces,
        triangulation_of,
    };
}

/// Offsets, shells, sweeps, lofts and drafts.
pub mod offset {
    pub use ogeom_offset::{
        Cap, Feature, FillBoundary, FillSide, Filled, Join, MiddlePath, PipeContact, PipeCorners,
        PipeEvolution, PipeLaw, PipeLaws, PipePlacement, Projected, SpineFunction, SpineLaw,
        ThreadForm, Walls, apply_draft, draft, feature, feature_prism, feature_revol, feature_rib,
        feature_slot, fill, loft, make_evolved, make_filling, make_filling_n, make_helical_sweep,
        make_loft, make_loft_skinned, make_loft_skinned_aligned, make_loft_skinned_closed,
        make_loft_solid, make_loft_surface, make_pipe, make_pipe_sections, make_pipe_shell,
        make_pipe_shell_evolving, make_pipe_shell_law, make_pipe_shell_laws, make_pipe_shell_with,
        make_pipe_skinned, make_prism_until, make_revolution_until, make_rib, make_ruled,
        make_sweep_surface, make_sweep_two_rails, make_thick_sheet, make_thick_sheet_sided,
        make_thick_solid, make_thick_solid_sided, make_thick_solid_with, make_thread, middle,
        middle_path, move_faces, normal_projection, offset_faces, offset_shape, offset_sheet,
        offset_wire, open_chain, place_on_spine, project, shape, sheet, sweep, thick,
        thicken_open_wire, thread, wire2d,
    };
}

/// The B-rep data model: shapes, topology, locations and the model that owns them.
pub mod topo {
    pub use ogeom_topo::{
        Absorbed, CurveId, Datum, DatumId, DatumStore, EdgeData, EdgeRepr, FaceData, Filter,
        GeometryStore, Location, MeshFinding, Model, ModelParts, NodeData, Orientation, PCurveId,
        PartnerKey, SameKey, Shape, ShapeType, SurfaceId, TShape, TShapeId, Triangulation,
        TriangulationId, Unshared, VertexData, WidenMark, ancestors_of, check_containment,
        enforce_containment, entity, explore, explore_unique, face_content, location, model, shape,
        tessellation,
    };
}

/// The kernel version, as reported by front ends.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
