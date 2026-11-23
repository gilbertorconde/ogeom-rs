//! What the mesh readers share: the colours a read carries, and the
//! triangulation every reader builds from a file's vertices, normals and
//! triangles.

use ogeom_math::{Point, Vector};
use ogeom_topo::Triangulation;

/// The colours a mesh read carries, in the most compact form that says
/// exactly what the file said.
///
/// A file colours a mesh at the corners of its triangles, or at its
/// vertices, or per face. Where every corner meeting at a vertex has the
/// same colour, they are given per vertex; else, where each triangle's
/// three corners agree, per triangle; else per corner.
#[derive(Debug, Clone, PartialEq)]
pub enum MeshColours {
    /// One RGBA colour per vertex of the mesh.
    PerVertex(Vec<[f64; 4]>),
    /// One RGBA colour per triangle.
    PerTriangle(Vec<[f64; 4]>),
    /// Three RGBA colours per triangle, at its corners in order.
    PerCorner(Vec<[[f64; 4]; 3]>),
}

impl MeshColours {
    /// The colour at corner `corner` (0, 1 or 2) of triangle `triangle` of
    /// `mesh`, whichever form the colours take; `None` past the end.
    #[must_use]
    pub fn at(&self, mesh: &Triangulation, triangle: usize, corner: usize) -> Option<[f64; 4]> {
        match self {
            Self::PerVertex(colours) => {
                let vertex = *mesh.triangles.get(triangle)?.get(corner)?;
                colours.get(vertex as usize).copied()
            }
            Self::PerTriangle(colours) => (corner < 3).then(|| colours.get(triangle).copied())?,
            Self::PerCorner(colours) => colours.get(triangle)?.get(corner).copied(),
        }
    }

    /// The compact form of colours given at every corner of `triangles`,
    /// over a mesh of `vertex_count` vertices.
    pub(crate) fn from_corners(
        triangles: &[[u32; 3]],
        vertex_count: usize,
        corners: Vec<[[f64; 4]; 3]>,
    ) -> Self {
        let mut per_vertex: Vec<Option<[f64; 4]>> = vec![None; vertex_count];
        let mut agree = true;
        'corners: for (triangle, colours) in triangles.iter().zip(&corners) {
            for (&v, colour) in triangle.iter().zip(colours) {
                match per_vertex.get_mut(v as usize) {
                    Some(slot @ None) => *slot = Some(*colour),
                    Some(Some(held)) if held == colour => {}
                    _ => {
                        agree = false;
                        break 'corners;
                    }
                }
            }
        }
        if agree && let Some(vertex) = per_vertex.into_iter().collect::<Option<Vec<_>>>() {
            return Self::PerVertex(vertex);
        }
        if corners.iter().all(|[a, b, c]| a == b && b == c) {
            return Self::PerTriangle(corners.into_iter().map(|[a, _, _]| a).collect());
        }
        Self::PerCorner(corners)
    }
}

/// Whether a normal slot is still waiting to be accumulated into.
///
/// A normal the file gave is a unit vector; a slot nobody has filled starts
/// at zero and grows by the faces around it. Telling them apart by length is
/// enough, and it means a file that gives *some* normals keeps them while
/// the rest are worked out.
pub(crate) fn held_is_unset(slot: &ogeom_math::Vector) -> bool {
    !(0.999..=1.001).contains(&slot.magnitude())
}

/// A mesh read from a file: its vertices, their normals and its triangles,
/// with the zero surface parameters a mesh from no surface holds.
pub(crate) fn imported(
    positions: Vec<Point>,
    normals: Vec<Vector>,
    triangles: Vec<[u32; 3]>,
    deflection_met: bool,
) -> Triangulation {
    let zeros = vec![(0.0, 0.0); positions.len()];
    let mut mesh = Triangulation {
        positions,
        normals,
        triangles,
        deflection_met,
        ..Triangulation::default()
    };
    mesh.set_vertex_parameters(zeros);
    mesh
}

/// Area-weighted vertex normals, for a file that gave none.
pub(crate) fn normals_from_triangles(positions: &[Point], triangles: &[[u32; 3]]) -> Vec<Vector> {
    let mut out = vec![Vector::new(0.0, 0.0, 0.0); positions.len()];
    for [a, b, c] in triangles {
        let (pa, pb, pc) = (
            positions[*a as usize],
            positions[*b as usize],
            positions[*c as usize],
        );
        // Not normalized: the cross product's length is twice the triangle's
        // area, which is exactly the weight a vertex normal wants.
        let n = (pb - pa).cross(pc - pa);
        for &k in &[*a, *b, *c] {
            out[k as usize] += n;
        }
    }
    for n in &mut out {
        let magnitude = n.magnitude();
        *n = if magnitude > 0.0 {
            *n / magnitude
        } else {
            Vector::new(0.0, 0.0, 1.0)
        };
    }
    out
}
