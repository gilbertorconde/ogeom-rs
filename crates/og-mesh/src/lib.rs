//! Tessellation. Discretize edges to a deflection tolerance, then build a
//! constrained Delaunay triangulation in each face's (u,v) parametric domain,
//! refine it, and attach the result to the face.
//!
//! *Elsewhere:* `BRepMesh` and the `Poly` triangulation types.
//!
//! Uses `spade` for CDT and refinement rather than ear-clipping: ear-clipping offers
//! no quality guarantees and cannot insert Steiner points, which curved surfaces
//! need.

pub mod attach;
pub mod discretize;
pub mod triangulate;

pub use attach::{Tessellated, polyline_of, tessellate, triangulation_of};
pub use discretize::{
    Deflection, Polyline, discretize, discretize_planar, is_straight, is_straight_planar,
};
pub use triangulate::{
    face_boundary, inside_boundary, polyline_of_edge, triangulate, triangulate_face,
};
