//! Blending: constant- and variable-radius edge fillets, vertex blends, chamfers,
//! and their 2D counterparts.

pub mod analyse;
mod between;
pub mod bridge;
pub mod chamfer;
mod consumed;
pub mod corner;
pub mod corner2d;
mod corner_curved;
pub mod facepair;
pub mod fillet;
mod junction;
pub mod march;
mod marched;
mod pair_marched;
mod pinched;
mod rolled;
mod ruled;
pub mod sheet;
mod sheet_curved;
mod support;
mod turned_pole;

pub use analyse::{BlendContact, analyse_blend, face_curvature_samples};
pub use bridge::{End, make_blend_curve, make_blend_surface};
pub use chamfer::{
    Chamfer, chamfer_edge, chamfer_edge_angle, chamfer_edge_distances, chamfer_edges,
    chamfer_edges_with,
};
pub use corner::round_vertex;
pub use corner2d::{chamfer_corner_2d, fillet_corner_2d};
pub use facepair::blend_faces;
pub use fillet::{fillet_edge, fillet_edge_variable, fillet_edges};
pub use march::{
    BlendStop, MarchedBlend, Sides, march_blend, march_blend_seeded, march_blend_sided,
};
pub use sheet::{fillet_faces, fillet_sheet_edges};
