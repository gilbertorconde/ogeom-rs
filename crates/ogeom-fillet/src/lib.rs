//! Blending: constant- and variable-radius edge fillets, vertex blends, chamfers,
//! and their 2D counterparts.

pub mod analyse;
mod asym;
mod asym_law;
mod between;
mod blame;
mod blend;
pub mod bridge;
pub mod chamfer;
mod chamfer_corner;
mod chamfer_marched;
mod chamfer_overflow;
mod chamfer_sphere;
mod compound;
pub mod corner2d;
mod cusp;
mod engine;
pub mod facepair;
pub mod fillet;
mod fits;
mod law;
mod limit;
pub mod march;
mod marched;
mod pair_marched;
mod pinched;
mod rim_arc;
#[allow(
    dead_code,
    reason = "the blend engine's sections and vertex blends, which its phases call"
)]
mod round;
mod setback;
mod sharp;
pub mod sheet;
mod sheet_curved;
mod support;
mod tangent_chain;
mod turned_pole;

pub use analyse::{BlendContact, analyse_blend, face_curvature_samples};
pub use asym::{
    RoundProfile, fillet_edge_distances, fillet_edge_distances_overflow, largest_fillet_distances,
};
pub use asym_law::{fillet_edge_distances_variable, fillet_edges_distances_variable};
pub use blame::{
    BlendRefusal, chamfer_edges_checked, chamfer_edges_partial, fillet_edges_checked,
    fillet_edges_corner_checked, fillet_edges_partial,
};
pub use blend::{Blend, blend_edges, blend_edges_overflow};
pub use bridge::{End, make_blend_curve, make_blend_surface};
pub use chamfer::{
    Chamfer, chamfer_edge, chamfer_edge_angle, chamfer_edge_distances, chamfer_edges,
    chamfer_edges_with, largest_chamfer_distance, largest_chamfer_distance_all,
    largest_chamfer_distance_all_overflow, largest_chamfer_distance_overflow,
};
pub use chamfer_corner::{ChamferCorner, chamfer_edges_cornered};
pub use chamfer_overflow::chamfer_edges_overflow;
pub use corner2d::{chamfer_corner_2d, fillet_corner_2d};
pub use engine::vertex::round_vertex;
pub use facepair::blend_faces;
pub use fillet::{
    FilletCorner, Overflow, fillet_edge, fillet_edge_variable, fillet_edges, fillet_edges_corner,
    fillet_edges_corner_overflow, fillet_edges_overflow, fillet_edges_variable,
    fillet_edges_variable_overflow, fillet_edges_with,
};
pub use fits::{largest_fillet_radius, largest_fillet_radius_all, largest_fillet_radius_overflow};
pub use law::RadiusLaw;
pub use limit::{
    BlendLimit, BlendSizing, DistanceLawSize, LIMIT_PROBES, LIMIT_SAMPLES, LIMIT_TOLERANCE,
    LawSize, LawStation, blend_limit, blend_seat, largest_law_distances, largest_law_radius,
    largest_law_scale,
};
pub use march::{
    BlendStop, MarchedBlend, Sides, march_blend, march_blend_seeded, march_blend_sided,
};
pub use ogeom_algo::{FillBoundary, FillSide, Filled, make_filling_n};
pub use setback::{Setbacks, fillet_edges_cornered, round_vertex_setback};
pub use sharp::round_vertex_setback_sharp;
pub use sheet::{fillet_faces, fillet_sheet_edges};
pub use tangent_chain::{ChainRules, TANGENT_CHAIN_ANGLE, tangent_chain, tangent_chain_with};
