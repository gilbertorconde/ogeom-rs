//! Distances to a triangle mesh, through its hierarchy.

use ogeom_core::{FastMap, OgeomResult, parallel};
use ogeom_math::{Point, Vector};
use ogeom_mesh::{Bvh, Hit};

use super::{Closest, DistanceOptions, MeshFeature, signed, unit, unit_or_zero};

/// The signed distance from each point to the mesh `index` was built over,
/// with the closest point, the feature it lies on and the normal there; in
/// the points' order, `None` where nothing is within
/// [`DistanceOptions::max`].
///
/// The closest point is [`Bvh::nearest`]'s, and with a facing limit
/// [`Bvh::nearest_where`]'s over the triangles facing within the limit of
/// the point's normal: exactly what brute force over the same triangles
/// finds, ties to the lowest triangle. A point that is not finite, or with
/// a facing limit a normal that has no direction, answers `None`. See the
/// [module](super) for how the sign is read.
///
/// # Errors
///
/// As [`DistanceOptions::check`];
/// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) if there are
/// normals but not one per point, and
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if there is a
/// facing limit and no normals.
pub fn to_mesh(
    index: &Bvh,
    points: &[Point],
    normals: Option<&[Vector]>,
    options: &DistanceOptions,
) -> OgeomResult<Vec<Option<Closest<MeshFeature>>>> {
    to_mesh_with(index, &Pseudonormals::of(index), points, normals, options)
}

/// As [`to_mesh`], with the mesh's normals already gathered from `index`,
/// for a caller that measures to the same mesh many times.
pub(crate) fn to_mesh_with(
    index: &Bvh,
    pseudo: &Pseudonormals,
    points: &[Point],
    normals: Option<&[Vector]>,
    options: &DistanceOptions,
) -> OgeomResult<Vec<Option<Closest<MeshFeature>>>> {
    let facing = options.checked(points.len(), normals)?;
    let hits = match (facing, normals) {
        (Some(cos), Some(normals)) => {
            let wanted: Vec<Option<Vector>> = normals.iter().map(|&n| unit(n)).collect();
            index.nearest_many_where(points, options.max, |i, t| {
                wanted[i].is_some_and(|m| pseudo.faces(t, m, cos))
            })
        }
        _ => index.nearest_many(points, options.max),
    };
    Ok(parallel::map_ordered(&hits, |i, hit| {
        hit.map(|hit| pseudo.closest(points[i], &hit, options.signed))
    }))
}

/// A mesh's normals at its faces, edges and vertices, the corners welded
/// where their positions are equal.
pub(crate) struct Pseudonormals {
    /// Each triangle's unit normal, zero for a degenerate one.
    face: Vec<Vector>,
    /// Each triangle's edges' pseudonormals, edge `k` from corner `k`: the
    /// sum of the unit normals of every triangle on the welded edge.
    edge: Vec<[Vector; 3]>,
    /// Each triangle's corners' welded vertices.
    corner: Vec<[u32; 3]>,
    /// Each welded vertex's pseudonormal: its triangles' unit normals
    /// weighted by their angles at it.
    vertex: Vec<Vector>,
}

/// A triangle's corners, or its edges by the corner each starts at.
const CORNERS: [u8; 3] = [0, 1, 2];

/// A position's key: its coordinates' bits, with negative zero as zero so
/// that equal positions weld.
fn key(p: Point) -> [u64; 3] {
    [p.x, p.y, p.z].map(|c| (c + 0.0).to_bits())
}

impl Pseudonormals {
    pub(crate) fn of(index: &Bvh) -> Self {
        let count = index.triangle_count();
        let triangles: Vec<[Point; 3]> = (0..count)
            .filter_map(|t| index.triangle(u32::try_from(t).ok()?))
            .collect();
        let mut welded: FastMap<[u64; 3], u32> = FastMap::default();
        let mut vertex: Vec<Vector> = Vec::new();
        let mut face = Vec::with_capacity(count);
        let mut corner = Vec::with_capacity(count);
        for t in &triangles {
            let ids = t.map(|p| {
                let next = u32::try_from(vertex.len()).unwrap_or(u32::MAX);
                let id = *welded.entry(key(p)).or_insert(next);
                if id == next {
                    vertex.push(Vector::ZERO);
                }
                id
            });
            let normal = unit_or_zero((t[1] - t[0]).cross(t[2] - t[0]));
            for k in 0..3 {
                let (a, b, c) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
                let (u, v) = (b - a, c - a);
                let angle = u.cross(v).magnitude().atan2(u.dot(v));
                vertex[ids[k] as usize] += normal * angle;
            }
            face.push(normal);
            corner.push(ids);
        }
        let mut edges: FastMap<(u32, u32), Vector> = FastMap::default();
        let edge_of = |ids: &[u32; 3], k: usize| {
            let (a, b) = (ids[k], ids[(k + 1) % 3]);
            (a.min(b), a.max(b))
        };
        for (ids, &normal) in corner.iter().zip(&face) {
            for k in 0..3 {
                let sum = edges.entry(edge_of(ids, k)).or_insert(Vector::ZERO);
                *sum += normal;
            }
        }
        let edge = corner
            .iter()
            .map(|ids| [0, 1, 2].map(|k| edges[&edge_of(ids, k)]))
            .collect();
        Self {
            face,
            edge,
            corner,
            vertex,
        }
    }

    /// Whether triangle `t` faces within the limit whose cosine is `cos` of
    /// the unit normal `m`. A degenerate triangle faces nowhere.
    fn faces(&self, t: u32, m: Vector, cos: f64) -> bool {
        let n = self.face[t as usize];
        n != Vector::ZERO && n.dot(m) >= cos
    }

    /// The hit as the feature it lies on, by its barycentric weights: a
    /// corner where two weights are zero, the edge between the other two
    /// where one is, the face where none is.
    fn closest(&self, p: Point, hit: &Hit, sign: bool) -> Closest<MeshFeature> {
        let t = hit.triangle;
        let at = t as usize;
        let zero = hit.barycentric.map(|w| w == 0.0);
        let (feature, normal) = match zero {
            [false, false, false] => (MeshFeature::Face(t), self.face[at]),
            [true, true, _] | [true, _, true] | [_, true, true] => {
                let k = zero.iter().position(|z| !z).unwrap_or(0);
                (
                    MeshFeature::Vertex(t, CORNERS[k]),
                    self.vertex[self.corner[at][k] as usize],
                )
            }
            _ => {
                let j = zero.iter().position(|&z| z).unwrap_or(0);
                let k = (j + 1) % 3;
                (MeshFeature::Edge(t, CORNERS[k]), self.edge[at][k])
            }
        };
        Closest {
            distance: signed(hit.distance, p, hit.point, normal, sign),
            point: hit.point,
            feature,
            normal: unit_or_zero(normal),
        }
    }
}
