//! Blending: constant- and variable-radius edge fillets, vertex blends, chamfers,
//! and their 2D counterparts.
//!
//! *Elsewhere:* `ChFi2d`, `ChFi3d`, `ChFiDS`, `Blend`, `BlendFunc`, `BRepBlend`
//! and `BRepFilletAPI`.

pub mod analyse;
pub mod bridge;
pub mod chamfer;
pub mod corner;
pub mod corner2d;
mod corner_curved;
pub mod facepair;
pub mod fillet;
pub mod march;
mod marched;
mod pinched;
mod ruled;
mod support;

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
