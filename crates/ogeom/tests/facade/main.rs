//! The integration tests of the `ogeom` facade, in one test program so the
//! crate links once.

#[path = "../support/cpu_time.rs"]
#[allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
mod cpu_time;
#[path = "../support/json.rs"]
#[allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
mod json;
#[path = "../support/pcurves.rs"]
#[allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
mod pcurves;
#[path = "../support/truth.rs"]
#[allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
mod truth;
#[path = "../support/walks.rs"]
#[allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
mod walks;

mod absorb;
mod analysis;
mod blend_after_refusal;
mod blends;
mod book;
mod boolean_chart_offsets;
mod boolean_contact;
mod boolean_contact_column;
mod boolean_far_nappe;
mod boolean_folded_rings;
mod boolean_grazing_drill;
mod boolean_history;
mod boolean_junction_trim;
mod boolean_near_coplanar;
mod boolean_nesting;
mod boolean_pole_drill;
mod boolean_rounded_drill;
mod boolean_scale;
mod boolean_tangent_bands;
mod bore_rim;
mod brep_format;
mod brep_saved;
mod built_by_hand;
mod canonical;
mod check_orientation;
mod classify_curved_seam;
mod classify_loose_rims;
mod cut_by_mesh_region;
mod defeature;
mod direct_edit;
mod divide;
mod document;
mod downstream_paths;
mod draft_measure;
mod earlier_solids;
mod edge_chains;
mod evolved;
mod exact_hlr;
mod face_meshes;
mod far_prism;
mod features;
mod filling;
mod fix_shape;
mod flush_ends;
mod fuse_all;
mod gltf_reading;
mod groove;
mod half_dome_sew;
mod half_space_curved;
mod hatch;
mod heal_trims;
mod helical_sweep;
mod helix;
mod iges_roundtrip;
mod iges_surface_decks;
mod instanced;
mod local_boolean;
mod loft;
mod loft_exact;
mod loft_solid;
mod loft_to_point;
mod mangled_files;
mod mass_chart;
mod measurement;
mod mechanical_part;
mod mesh_chart_bow;
mod mesh_chord_rims;
mod mesh_corners_off_vertices;
mod mesh_corpus;
mod mesh_editable;
mod mesh_offset;
mod mesh_section;
mod mesh_shadow;
mod mesh_smooth;
mod mesh_solid;
mod mesh_solid_checked;
mod mesh_steps;
mod mesh_touching_rims;
mod mirror_exchange;
mod mirrored;
mod native_subset;
mod near_closed_torus;
mod network;
mod nonmanifold;
mod offset_solid;
mod pipe_corners;
mod pipe_evolving;
mod pipe_laws;
mod pipe_placed;
mod pipe_rings;
mod pipe_shell_exact;
mod project_plane;
mod refine_exact;
mod remodel;
mod retain_reachable;
mod reversed_holed_face;
mod revolution_until;
mod revolve_on_axis;
mod saved_views;
mod scaled_copy;
mod section_face;
mod sheet_fillets;
mod sheet_inverse;
mod sheet_offset;
mod sheet_wrap;
mod sphere_exchange;
mod sphere_pole_trim;
mod spline_measure;
mod step_callouts;
mod step_conformance;
mod step_extrusions;
mod step_geometry_kinds;
mod step_surface_curves;
mod step_unseamed_faces;
mod surface_lofts;
mod sweep_orientation;
mod sweep_profile_pcurves;
mod thickness_own_faces;
mod thread_groove;
mod tight_bounds;
mod trimmed_hole_mass;
mod uneven_scale;
mod unit_scales;
mod upgrade;
mod void_exchange;
mod vrml_read;
