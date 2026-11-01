//! Tessellation. Discretize edges to a deflection tolerance, then build a
//! constrained Delaunay triangulation in each face's (u,v) parametric domain,
//! refine it, and attach the result to the face.
//!
//! Uses `spade` for CDT and refinement rather than ear-clipping: ear-clipping offers
//! no quality guarantees and cannot insert Steiner points, which curved surfaces
//! need.
//!
//! Also the spatial indexes over meshes and point clouds that inspection
//! asks of them: [`Bvh`] over triangles and [`PointIndex`] over points.

pub mod attach;
pub mod bvh;
pub mod discretize;
pub mod hatch;
pub mod point_index;
pub mod simplify;
mod tails;
pub mod triangulate;

pub use attach::{Tessellated, polyline_of, tessellate, triangulation_of};
pub use bvh::{Bvh, Hit};
pub use discretize::{
    Deflection, Polyline, discretize, discretize_on_surface, discretize_planar, is_straight,
    is_straight_planar,
};
pub use hatch::hatch_face;
pub use point_index::{Neighbour, PointIndex};
pub use simplify::{Simplified, Target, simplify};
pub use triangulate::{
    EdgeChords, FaceMesh, FaceMeshCache, edge_chords_for, face_boundary, face_boundary_in_chart,
    face_boundary_in_chart_met, face_meshes_for, inside_boundary, inside_boundary_with,
    open_chart_ring, polyline_of_edge, triangulate, triangulate_face, triangulate_face_kept,
    triangulate_face_with, triangulate_with_chords, triangulate_with_face_meshes,
};
