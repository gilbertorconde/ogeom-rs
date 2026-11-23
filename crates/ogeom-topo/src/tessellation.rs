//! Cached tessellation: the polyline and triangle forms of exact geometry.
//!
//! A triangulation is a *representation* of a face, not a replacement for it
//! (`docs/DATA_MODEL.md` §6). It lives here, beside the entity data, because
//! that is what it belongs to: a face holds one the way an edge holds a
//! pcurve, and the algorithms that build it live a layer up in `ogeom-mesh`.
//!
//! Everything here is plain data with the queries that read it. Nothing here
//! decides how finely to sample anything.

use ogeom_core::Tolerances;
use ogeom_math::{Aabb, Point, Vector};

/// A triangulated surface.
///
/// Triangle corners carry their surface parameters as well as their
/// vertices' positions, so a caller can ask the exact surface about a
/// triangulated point rather than only the approximation.
///
/// Parameters and normals are kept per corner, not per vertex, because a
/// point of the surface can have more than one: a cone's apex or a
/// sphere's pole is one vertex, the fan round it a ring of real triangles,
/// and the `u` it has depends on the triangle it is read from. A vertex on
/// a seam has a `u` at either end of the period in the same way.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Triangulation {
    /// Vertex positions.
    pub positions: Vec<Point>,
    /// Outward unit normals, one per vertex, or none for a mesh read
    /// without them. Where a vertex has no single normal (a cone's tip),
    /// this is one of its corners' and [`Self::corner_normals`] has each.
    pub normals: Vec<Vector>,
    /// The surface parameters each vertex came from, or none for a mesh
    /// that did not come from a surface. A vertex with more than one (a
    /// cone's apex, a sphere's pole) holds one of them here.
    #[deprecated(
        since = "0.20.0",
        note = "a vertex at an apex or a pole has no single (u, v); read                 `corner_parameters` (or `Triangulation::triangle_parameters`)                 instead. Still filled for one release."
    )]
    pub parameters: Vec<(f64, f64)>,
    /// Triangles, as indices into the vertex arrays, wound counter-clockwise
    /// about the outward normal.
    pub triangles: Vec<[u32; 3]>,
    /// The surface parameters at each triangle's corners, in the order of
    /// its indices, one entry per triangle; or none for a mesh that did not
    /// come from a surface.
    ///
    /// The mean of a triangle's three is a point inside the triangle's
    /// patch of the chart, so the surface there is the point the triangle
    /// stands in for. At an apex or a pole, where any `u` names the same
    /// point, the corner holds the apex's own `v` and the middle of the
    /// other two corners' `u`.
    pub corner_parameters: Vec<[(f64, f64); 3]>,
    /// Outward unit normals at each triangle's corners, in the order of its
    /// indices, one entry per triangle; or none for a mesh without normals.
    /// At a cone's tip each corner has the normal of the ruling its
    /// triangle meets the tip along.
    pub corner_normals: Vec<[Vector; 3]>,
    /// Whether every face met its requested deflection.
    pub deflection_met: bool,
    /// What the mesher found about its faces that it could not mesh as
    /// asked, each by name; empty for a mesh that did not come from one.
    pub findings: Vec<MeshFinding>,
}

/// Something a mesher found about a face that kept it from meshing the face
/// as asked, named, with where.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum MeshFinding {
    /// The face's surface folds at a point inside the face: its two tangents
    /// run parallel there, neither of them vanishing, and its normal turns
    /// over round the point. No flat triangle across the point faces the
    /// way the surface does on every side of it.
    SurfaceFolds {
        /// The surface parameters of the point.
        at: (f64, f64),
        /// The point, placed where the face is.
        point: Point,
    },
    /// The chord repair stopped at its cap on rounds while it was still
    /// refining: somewhere the mesh may stand off the face by more than the
    /// chord.
    RepairCapped {
        /// The surface parameters of a point it would have refined at next.
        at: (f64, f64),
        /// That point on the surface, placed where the face is.
        point: Point,
    },
}

impl Triangulation {
    /// An empty mesh.
    #[must_use]
    pub fn new() -> Self {
        Self {
            deflection_met: true,
            ..Self::default()
        }
    }

    /// Number of vertices.
    #[must_use]
    pub fn vertex_count(&self) -> usize {
        self.positions.len()
    }

    /// Number of triangles.
    #[must_use]
    pub fn triangle_count(&self) -> usize {
        self.triangles.len()
    }

    /// Whether the mesh holds no triangles.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.triangles.is_empty()
    }

    /// The surface parameters at triangle `index`'s corners, in the order
    /// of its indices.
    ///
    /// Read from [`Self::corner_parameters`] where it holds one entry per
    /// triangle, otherwise from the vertices' own, for a mesh built with
    /// those alone. `None` for a mesh with neither, or past the triangles.
    #[must_use]
    pub fn triangle_parameters(&self, index: usize) -> Option<[(f64, f64); 3]> {
        if self.corner_parameters.len() == self.triangles.len() {
            return self.corner_parameters.get(index).copied();
        }
        let t = self.triangles.get(index)?;
        #[allow(deprecated, reason = "the fallback for a mesh built before corners")]
        let per_vertex = &self.parameters;
        Some([
            *per_vertex.get(t[0] as usize)?,
            *per_vertex.get(t[1] as usize)?,
            *per_vertex.get(t[2] as usize)?,
        ])
    }

    /// The normals at triangle `index`'s corners, as
    /// [`Self::triangle_parameters`] reads the parameters.
    #[must_use]
    pub fn triangle_normals(&self, index: usize) -> Option<[Vector; 3]> {
        if self.corner_normals.len() == self.triangles.len() {
            return self.corner_normals.get(index).copied();
        }
        let t = self.triangles.get(index)?;
        Some([
            *self.normals.get(t[0] as usize)?,
            *self.normals.get(t[1] as usize)?,
            *self.normals.get(t[2] as usize)?,
        ])
    }

    /// Whether every triangle has its corners' surface parameters.
    #[must_use]
    pub fn has_parameters(&self) -> bool {
        #[allow(deprecated, reason = "the fallback for a mesh built before corners")]
        let per_vertex = self.parameters.len() == self.positions.len();
        !self.triangles.is_empty()
            && (self.corner_parameters.len() == self.triangles.len() || per_vertex)
    }

    /// The surface parameters a corner of each vertex has, `None` for a
    /// vertex no triangle uses or a mesh without parameters.
    ///
    /// A vertex at an apex or a pole has more than one; this is the first
    /// triangle's that uses it, which names the same point of the surface.
    /// For a value tied to a triangle, read [`Self::triangle_parameters`].
    #[must_use]
    pub fn vertex_parameters(&self) -> Vec<Option<(f64, f64)>> {
        let mut out = vec![None; self.positions.len()];
        for index in 0..self.triangles.len() {
            let Some(corners) = self.triangle_parameters(index) else {
                break;
            };
            for (&v, uv) in self.triangles[index].iter().zip(corners) {
                if let Some(slot) = out.get_mut(v as usize)
                    && slot.is_none()
                {
                    *slot = Some(uv);
                }
            }
        }
        out
    }

    /// Set every vertex's surface parameters, and each corner's from its
    /// vertex: for a mesh whose vertices each have one `(u, v)`. Call it
    /// once the triangles are final.
    pub fn set_vertex_parameters(&mut self, parameters: Vec<(f64, f64)>) {
        self.corner_parameters = self
            .triangles
            .iter()
            .map(|t| t.map(|i| parameters.get(i as usize).copied().unwrap_or((0.0, 0.0))))
            .collect();
        if parameters.len() != self.positions.len() {
            self.corner_parameters.clear();
        }
        #[allow(deprecated, reason = "kept filled for the release it is deprecated in")]
        {
            self.parameters = parameters;
        }
    }

    /// Set each corner's normal from its vertex's, for a mesh whose vertices
    /// each have one. Call it once the triangles are final.
    pub fn set_corner_normals_from_vertices(&mut self) {
        self.corner_normals = if self.normals.len() == self.positions.len() {
            self.triangles
                .iter()
                .map(|t| t.map(|i| self.normals[i as usize]))
                .collect()
        } else {
            Vec::new()
        };
    }

    /// The corner arrays for triangles taken from this mesh's `sources`,
    /// each array empty unless this mesh has one entry per triangle.
    fn corners_of(&self, sources: &[usize]) -> Corners {
        (
            picked(&self.corner_parameters, self.triangles.len(), sources),
            picked(&self.corner_normals, self.triangles.len(), sources),
        )
    }

    /// Every triangle's corner parameters, read as
    /// [`Self::triangle_parameters`] reads them; `None` where some triangle
    /// has none.
    fn all_corner_parameters(&self) -> Option<Vec<[(f64, f64); 3]>> {
        if self.corner_parameters.len() == self.triangles.len() {
            return Some(self.corner_parameters.clone());
        }
        (0..self.triangles.len())
            .map(|i| self.triangle_parameters(i))
            .collect()
    }

    /// Every triangle's corner normals, as [`Self::all_corner_parameters`].
    fn all_corner_normals(&self) -> Option<Vec<[Vector; 3]>> {
        if self.corner_normals.len() == self.triangles.len() {
            return Some(self.corner_normals.clone());
        }
        (0..self.triangles.len())
            .map(|i| self.triangle_normals(i))
            .collect()
    }

    /// The bounding box of the vertices.
    #[must_use]
    pub fn bounds(&self) -> Aabb {
        Aabb::of_points(&self.positions)
    }

    /// The total area of the triangles.
    ///
    /// An *under*estimate of the surface's own area for a convex patch, since a
    /// triangle chord-cuts the surface it spans. It converges from below as the
    /// deflection tightens.
    #[must_use]
    pub fn area(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| self.positions[i as usize]);
                (b - a).cross(c - a).magnitude() * 0.5
            })
            .sum()
    }

    /// The signed volume enclosed, by the divergence theorem.
    ///
    /// Meaningful only for a mesh that is closed and consistently wound
    /// outward: each triangle contributes the signed volume of the tetrahedron
    /// it forms with the mesh's first position, and the contributions cancel
    /// except over the enclosed region. An open mesh gives a number with no
    /// meaning, and a mesh wound inward gives the negative, which is why
    /// [`Triangulation::is_closed`] exists to be asked first.
    ///
    /// The apex is a point of the mesh rather than the origin so a mesh far
    /// from the origin sums tetrahedra of its own size, not slivers of huge
    /// ones whose volumes cancel to rounding.
    #[must_use]
    pub fn volume(&self) -> f64 {
        let Some(&apex) = self.positions.first() else {
            return 0.0;
        };
        self.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| self.positions[i as usize] - apex);
                a.dot(b.cross(c)) / 6.0
            })
            .sum()
    }

    /// Whether every triangle edge is crossed as often one way as the other.
    ///
    /// The mesh equivalent of a closed shell, and the precondition for
    /// [`Triangulation::volume`] meaning anything. It is exactly what the
    /// divergence theorem needs: the surface has no boundary, and it is
    /// wound consistently, so each triangle's contribution cancels against
    /// its neighbours' except over the region enclosed.
    ///
    /// Counting *directed* edges rather than undirected ones is what makes
    /// this the right question, and it is stricter and looser than the
    /// obvious test in the two different ways that matter.
    ///
    /// Stricter: two triangles sharing an edge and winding the *same* way
    /// round it traverse it twice in the same direction. The edge is used
    /// twice, so a count of uses calls it closed, and the volume that comes
    /// out is wrong because one of the two faces is inside out.
    ///
    /// Looser: an edge may legitimately carry four triangles. Where two
    /// faces meet along a short edge that discretizes into several segments,
    /// each can fill the sliver between the polyline and its own chord, and
    /// the chord then belongs to both: four triangles round one edge, two
    /// crossing each way. There is no hole there and the volume is right.
    /// Demanding exactly two refuses a mesh for being non-manifold when
    /// nothing was asked about manifoldness.
    ///
    /// This agrees with the topology side:
    /// [`is_shell_closed`](../../ogeom_algo/fn.is_shell_closed.html) counts an
    /// edge's uses and accepts any even number, and the two halves of the
    /// kernel should not mean different things by the same word.
    #[must_use]
    pub fn is_closed(&self) -> bool {
        use ogeom_core::FastMap;
        let mut balance: FastMap<(u32, u32), i64> = FastMap::default();
        for t in &self.triangles {
            for i in 0..3 {
                let (a, b) = (t[i], t[(i + 1) % 3]);
                // One key per undirected edge; the direction decides the sign.
                let (key, step) = if a <= b { ((a, b), 1) } else { ((b, a), -1) };
                *balance.entry(key).or_default() += step;
            }
        }
        !balance.is_empty() && balance.values().all(|&n| n == 0)
    }

    /// Weld only the mesh's *border* vertices, within `reach`.
    ///
    /// The second pass after [`Triangulation::welded`]: interior edges are
    /// already manifold, and touching them at a widened tolerance would eat
    /// real features. Borders are where imported slop lives (an edge's curve
    /// and its neighbour's disagree by the file's own tolerance, which the
    /// model records on the edge), so only vertices on unmatched triangle
    /// edges are candidates, merged to their nearest counterpart within
    /// `reach`.
    #[must_use]
    pub fn border_welded(&self, reach: f64) -> Self {
        self.border_welded_with_sources(reach).0
    }

    /// As [`Self::border_welded`], with the input triangle each output
    /// triangle is.
    #[must_use]
    pub fn border_welded_with_sources(&self, reach: f64) -> (Self, Vec<usize>) {
        use ogeom_core::FastMap;
        if !reach.is_finite() || reach <= 0.0 {
            return (self.clone(), (0..self.triangles.len()).collect());
        }
        // Border vertices: endpoints of triangle edges used an odd number of
        // times.
        let uses = edge_uses(&self.triangles);
        let mut border: Vec<u32> = uses
            .iter()
            .filter(|&(_, &n)| n % 2 == 1)
            .flat_map(|(&(a, b), _)| [a, b])
            .collect();
        border.sort_unstable();
        border.dedup();
        if border.is_empty() {
            return (self.clone(), (0..self.triangles.len()).collect());
        }

        // Cluster border vertices within reach, first-seen wins, checked
        // against the cluster representative so chains cannot creep.
        let cell = reach.max(f64::MIN_POSITIVE);
        let key = |p: Point| {
            #[allow(clippy::cast_possible_truncation)]
            (
                (p.x / cell).round() as i64,
                (p.y / cell).round() as i64,
                (p.z / cell).round() as i64,
            )
        };
        let mut buckets: FastMap<(i64, i64, i64), Vec<u32>> = FastMap::default();
        #[allow(clippy::cast_possible_truncation)]
        let mut remap: Vec<u32> = (0..self.positions.len() as u32).collect();
        for &v in &border {
            let p = self.positions[v as usize];
            let (kx, ky, kz) = key(p);
            let mut found = None;
            'search: for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        for &candidate in buckets
                            .get(&(
                                kx.saturating_add(dx),
                                ky.saturating_add(dy),
                                kz.saturating_add(dz),
                            ))
                            .map_or(&[][..], Vec::as_slice)
                        {
                            if self.positions[candidate as usize].distance(p) <= reach {
                                found = Some(candidate);
                                break 'search;
                            }
                        }
                    }
                }
            }
            match found {
                Some(rep) => remap[v as usize] = rep,
                None => buckets.entry((kx, ky, kz)).or_default().push(v),
            }
        }

        let mut out = Self::new();
        out.deflection_met = self.deflection_met;
        out.findings.clone_from(&self.findings);
        // Every vertex is kept, merged ones included, and the triangles are
        // remapped onto the representatives.
        out.positions = self.positions.clone();
        out.normals = self.normals.clone();
        #[allow(deprecated, reason = "kept filled for the release it is deprecated in")]
        {
            out.parameters = self.parameters.clone();
        }
        let mut sources = Vec::new();
        for (source, t) in self.triangles.iter().enumerate() {
            let mapped = t.map(|i| remap[i as usize]);
            if mapped[0] != mapped[1] && mapped[1] != mapped[2] && mapped[2] != mapped[0] {
                out.triangles.push(mapped);
                sources.push(source);
            }
        }
        (out.corner_parameters, out.corner_normals) = self.corners_of(&sources);
        (out, sources)
    }

    /// Split border segments at border vertices that lie on them.
    ///
    /// The T-junction repair that follows [`Triangulation::border_welded`]:
    /// after welding, two faces' border chains share their vertices but may
    /// subdivide the same stretch differently: one face's segment spans two
    /// of its neighbour's. Splitting the long segment *at the neighbour's own
    /// vertex index* makes the chains segment-for-segment identical, which is
    /// what closure counts. No positions move and none are added.
    #[must_use]
    pub fn border_stitched(&self, reach: f64) -> Self {
        self.border_stitched_with_sources(reach).0
    }

    /// As [`Self::border_stitched`], retaining the source triangle of every output triangle.
    #[must_use]
    pub fn border_stitched_with_sources(&self, reach: f64) -> (Self, Vec<usize>) {
        use ogeom_core::FastMap;
        if !reach.is_finite() || reach <= 0.0 {
            return (self.clone(), (0..self.triangles.len()).collect());
        }
        let uses = edge_uses(&self.triangles);
        let border_edges: Vec<(u32, u32)> = uses
            .iter()
            .filter(|&(_, &n)| n % 2 == 1)
            .map(|(&e, _)| e)
            .collect();
        if border_edges.is_empty() {
            return (self.clone(), (0..self.triangles.len()).collect());
        }
        let mut border_vertices: Vec<u32> =
            border_edges.iter().flat_map(|&(a, b)| [a, b]).collect();
        border_vertices.sort_unstable();
        border_vertices.dedup();
        // Cells about a border segment long, so a segment's box spans few.
        #[allow(clippy::cast_precision_loss, reason = "a mean over a count")]
        let spacing = border_edges
            .iter()
            .map(|&(a, b)| self.positions[a as usize].distance(self.positions[b as usize]))
            .sum::<f64>()
            / border_edges.len() as f64;
        let grid = PointGrid::over(&self.positions, &border_vertices, spacing, reach);

        // For every border segment, the border vertices sitting on its
        // interior, ordered along it (by vertex where two stand level).
        let mut splits: FastMap<(u32, u32), Vec<u32>> = FastMap::default();
        let mut near: Vec<u32> = Vec::new();
        for &(a, b) in &border_edges {
            let (pa, pb) = (self.positions[a as usize], self.positions[b as usize]);
            let d = pb - pa;
            let l2 = d.dot(d);
            if l2 <= 0.0 {
                continue;
            }
            near.clear();
            match &grid {
                Some(grid) => grid.within(pa, pb, reach, &mut near),
                None => near.extend_from_slice(&border_vertices),
            }
            let mut on: Vec<(f64, u32)> = near
                .iter()
                .filter(|&&v| v != a && v != b)
                .filter_map(|&v| {
                    let p = self.positions[v as usize];
                    let t = (p - pa).dot(d) / l2;
                    if !(0.001..=0.999).contains(&t) {
                        return None;
                    }
                    ((pa + d * t).distance(p) <= reach).then_some((t, v))
                })
                .collect();
            if on.is_empty() {
                continue;
            }
            on.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
            splits.insert((a, b), on.into_iter().map(|(_, v)| v).collect());
        }
        if splits.is_empty() {
            return (self.clone(), (0..self.triangles.len()).collect());
        }

        let mut out = Self::new();
        out.deflection_met = self.deflection_met;
        out.findings.clone_from(&self.findings);
        out.positions = self.positions.clone();
        out.normals = self.normals.clone();
        #[allow(deprecated, reason = "kept filled for the release it is deprecated in")]
        {
            out.parameters = self.parameters.clone();
        }
        let mut sources = Vec::new();
        for (source, t) in self.triangles.iter().enumerate() {
            // Split each original side from the opposite corner. A fan from
            // the ring's first vertex would leave its incident sides unsplit
            // and emit collinear triangles along them.
            let sides = std::array::from_fn(|i| {
                let (a, b) = (t[i], t[(i + 1) % 3]);
                splits
                    .get(&(a.min(b), a.max(b)))
                    .map_or(&[][..], Vec::as_slice)
            });
            let nearest = nearest_split_sides(t, sides, |i, v| {
                let (a, b) = (t[i], t[(i + 1) % 3]);
                let (pa, pb, p) = (
                    self.positions[a as usize],
                    self.positions[b as usize],
                    self.positions[v as usize],
                );
                let d = pb - pa;
                let along = (p - pa).dot(d) / d.dot(d);
                (pa + d * along).distance(p)
            });
            let mut pieces = vec![*t];
            for i in 0..3 {
                let (a, b) = (t[i], t[(i + 1) % 3]);
                if let Some(vs) = splits.get(&(a.min(b), a.max(b))) {
                    // At a narrow corner a vertex can be within reach of two
                    // sides. Insert it only on the nearest: inserting it on
                    // both turns the polygon back through the same vertex.
                    let vs: Vec<_> = vs
                        .iter()
                        .copied()
                        .filter(|v| nearest.get(v) == Some(&i))
                        .collect();
                    if vs.is_empty() {
                        continue;
                    }
                    let mut chain = vec![a];
                    if a < b {
                        chain.extend(vs.iter().copied());
                    } else {
                        chain.extend(vs.iter().rev().copied());
                    }
                    chain.push(b);
                    let Some(at) = pieces
                        .iter()
                        .position(|piece| (0..3).any(|k| piece[k] == a && piece[(k + 1) % 3] == b))
                    else {
                        panic!(
                            "source triangle {source}: original side {a}->{b} is missing from current pieces {pieces:?}"
                        );
                    };
                    let piece = pieces.remove(at);
                    let Some(&opposite) = piece.iter().find(|&&v| v != a && v != b) else {
                        panic!(
                            "source triangle {source}: piece {piece:?} has no corner opposite side {a}->{b}"
                        );
                    };
                    pieces.extend(chain.windows(2).map(|side| [side[0], side[1], opposite]));
                }
            }
            for tri in pieces {
                out.triangles.push(tri);
                sources.push(source);
            }
        }
        self.split_corners(&mut out, &sources);
        (out, sources)
    }

    /// Give `out`'s triangles, each a piece of this mesh's triangle
    /// `sources[i]`, their corners' parameters and normals: a corner of the
    /// source keeps its own, a vertex put on a side of it is interpolated
    /// along that side.
    fn split_corners(&self, out: &mut Self, sources: &[usize]) {
        let full_p = !self.corner_parameters.is_empty()
            && self.corner_parameters.len() == self.triangles.len();
        let full_n =
            !self.corner_normals.is_empty() && self.corner_normals.len() == self.triangles.len();
        if !full_p && !full_n {
            return;
        }
        for (piece, &source) in out.triangles.iter().zip(sources) {
            let t = self.triangles[source];
            let weights = piece.map(|v| self.weights_in(t, v));
            if full_p {
                let uv = self.corner_parameters[source];
                out.corner_parameters.push(weights.map(|w| {
                    (
                        w[0] * uv[0].0 + w[1] * uv[1].0 + w[2] * uv[2].0,
                        w[0] * uv[0].1 + w[1] * uv[1].1 + w[2] * uv[2].1,
                    )
                }));
            }
            if full_n {
                let n = self.corner_normals[source];
                out.corner_normals.push(weights.map(|w| {
                    let blend = n[0] * w[0] + n[1] * w[1] + n[2] * w[2];
                    let length = blend.magnitude();
                    if length > 0.0 { blend / length } else { n[0] }
                }));
            }
        }
    }

    /// Vertex `v`'s weights on triangle `t`'s corners: one on the corner it
    /// is, else its nearest point's on the side it lies nearest.
    fn weights_in(&self, t: [u32; 3], v: u32) -> [f64; 3] {
        if let Some(k) = t.iter().position(|&c| c == v) {
            let mut w = [0.0; 3];
            w[k] = 1.0;
            return w;
        }
        let p = self.positions[v as usize];
        let mut best = (f64::INFINITY, [1.0, 0.0, 0.0]);
        for k in 0..3 {
            let (a, b) = (
                self.positions[t[k] as usize],
                self.positions[t[(k + 1) % 3] as usize],
            );
            let d = b - a;
            let l2 = d.dot(d);
            let along = if l2 > 0.0 {
                ((p - a).dot(d) / l2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let off = (a + d * along).distance(p);
            if off < best.0 {
                let mut w = [0.0; 3];
                w[k] = 1.0 - along;
                w[(k + 1) % 3] = along;
                best = (off, w);
            }
        }
        best.1
    }

    /// This mesh with its folds cancelled and its cracks sealed: the last
    /// pass over a mesh welded from faces that each met their edges.
    ///
    /// A fold is two triangles on the same three vertices facing opposite
    /// ways, left where a sliver face collapses in the weld; they enclose
    /// nothing, and both go. A crack is a loop of border edges whose mean
    /// width (twice its area over its perimeter) is within `width`: two
    /// faces sampling a shared corner differently leave one, narrower than
    /// the chord they were drawn to. It is fanned shut from one of its
    /// corners, each new triangle crossing a border edge the other way from
    /// the triangle already on it. A loop wider than that, or a border
    /// vertex with more than one way on, is a real opening and stays.
    #[must_use]
    pub fn sealed(&self, width: f64) -> Self {
        self.sealed_with_sources(width).0
    }

    /// As [`Self::sealed`], with an input triangle for each output triangle:
    /// the one it is, or for a triangle fanned across a crack, the one on
    /// the other side of the crack edge it closes, so a crack's triangles
    /// belong to the triangles round it.
    #[must_use]
    pub fn sealed_with_sources(&self, width: f64) -> (Self, Vec<usize>) {
        use ogeom_core::FastMap;
        let (mut out, mut sources) = self.without_folds_with_sources();
        if !width.is_finite() || width <= 0.0 {
            return (out, sources);
        }
        // Cracks: each border edge walked the other way from its triangle.
        let uses = edge_uses(&out.triangles);
        let mut onward: FastMap<u32, Vec<u32>> = FastMap::default();
        // The triangle walking each border edge, by its direction there.
        let mut owner: FastMap<(u32, u32), usize> = FastMap::default();
        for (index, t) in out.triangles.iter().enumerate() {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                if uses[&(a.min(b), a.max(b))] == 1 {
                    onward.entry(b).or_default().push(a);
                    owner.insert((a, b), index);
                }
            }
        }
        let mut done: ogeom_core::FastSet<u32> = ogeom_core::FastSet::default();
        let filled_from = out.triangles.len();
        let mut starts: Vec<u32> = onward.keys().copied().collect();
        starts.sort_unstable();
        // The walk each vertex was last put on a ring by, plus one.
        let mut on_walk: Vec<u32> = vec![0; out.positions.len()];
        #[allow(clippy::cast_possible_truncation)]
        for (walk, start) in starts.into_iter().enumerate() {
            let walk = walk as u32 + 1;
            if done.contains(&start) {
                continue;
            }
            let mut ring = vec![start];
            if let Some(mark) = on_walk.get_mut(start as usize) {
                *mark = walk;
            }
            let mut at = start;
            let closed = loop {
                let Some(next) = onward.get(&at) else {
                    break false;
                };
                let [next] = next[..] else {
                    break false;
                };
                if next == start {
                    break true;
                }
                let seen = match on_walk.get(next as usize) {
                    Some(&mark) => mark == walk,
                    None => ring.contains(&next),
                };
                if seen || ring.len() > onward.len() {
                    break false;
                }
                ring.push(next);
                if let Some(mark) = on_walk.get_mut(next as usize) {
                    *mark = walk;
                }
                at = next;
            };
            for &v in &ring {
                done.insert(v);
            }
            if !closed || ring.len() < 3 {
                continue;
            }
            let points: Vec<Point> = ring.iter().map(|&v| out.positions[v as usize]).collect();
            // About the ring's first point: the sum is the same about any
            // point, and about a far one it is a difference of huge
            // products.
            let anchor = points[0];
            let mut normal = Vector::ZERO;
            let mut perimeter = 0.0;
            for (i, p) in points.iter().enumerate() {
                let q = points[(i + 1) % points.len()];
                normal += (*p - anchor).cross(q - anchor);
                perimeter += p.distance(q);
            }
            let area = normal.magnitude() / 2.0;
            if perimeter <= 0.0 || 2.0 * area / perimeter > width {
                continue;
            }
            for i in 1..ring.len() - 1 {
                out.triangles.push([ring[0], ring[i], ring[i + 1]]);
                // The walk steps from `at` to `onward[at]`, so the
                // triangle on side (ring[i], ring[i + 1]) walks it the
                // other way.
                let across = owner
                    .get(&(ring[i + 1], ring[i]))
                    .map_or(0, |&t| sources[t]);
                sources.push(across);
            }
        }
        out.fill_corners_from(filled_from);
        (out, sources)
    }

    /// Give the triangles from `first` on, added without corners, the
    /// corners another triangle has at the same vertices: a triangle
    /// sealing a crack joins vertices of the triangles round it.
    fn fill_corners_from(&mut self, first: usize) {
        let added = first..self.triangles.len();
        if added.is_empty() {
            return;
        }
        if self.corner_parameters.len() == first && first > 0 {
            let mut at: Vec<Option<(f64, f64)>> = vec![None; self.positions.len()];
            for (t, uv) in self.triangles.iter().zip(&self.corner_parameters) {
                for (&v, &c) in t.iter().zip(uv) {
                    at[v as usize].get_or_insert(c);
                }
            }
            for i in added.clone() {
                let t = self.triangles[i];
                self.corner_parameters
                    .push(t.map(|v| at[v as usize].unwrap_or((0.0, 0.0))));
            }
        } else {
            self.corner_parameters.clear();
        }
        if self.corner_normals.len() == first && first > 0 {
            let mut at: Vec<Option<Vector>> = vec![None; self.positions.len()];
            for (t, n) in self.triangles.iter().zip(&self.corner_normals) {
                for (&v, &c) in t.iter().zip(n) {
                    at[v as usize].get_or_insert(c);
                }
            }
            for i in added {
                let t = self.triangles[i];
                self.corner_normals
                    .push(t.map(|v| at[v as usize].unwrap_or(Vector::ZERO)));
            }
        } else {
            self.corner_normals.clear();
        }
    }

    /// Cancel opposite triangles on the same vertices, retaining each survivor's source.
    /// No positions move and no crack is filled.
    #[must_use]
    pub fn without_folds_with_sources(&self) -> (Self, Vec<usize>) {
        let mut keyed: Vec<([u32; 3], usize)> = self
            .triangles
            .iter()
            .enumerate()
            .map(|(i, t)| {
                let mut key = *t;
                key.sort_unstable();
                (key, i)
            })
            .collect();
        keyed.sort_unstable();
        let mut drop = vec![false; self.triangles.len()];
        let mut open: Vec<usize> = Vec::new();
        for group in keyed.chunk_by(|x, y| x.0 == y.0) {
            if group.len() < 2 {
                continue;
            }
            open.clear();
            for &(_, i) in group {
                let t = self.triangles[i];
                let reverse = open.iter().position(|&j| {
                    let u = self.triangles[j];
                    (0..3).any(|k| [u[k], u[(k + 2) % 3], u[(k + 1) % 3]] == t)
                });
                match reverse {
                    Some(at) => {
                        drop[open.remove(at)] = true;
                        drop[i] = true;
                    }
                    None => open.push(i),
                }
            }
        }
        let sources: Vec<_> = (0..self.triangles.len()).filter(|&i| !drop[i]).collect();
        let mut out = self.clone();
        out.triangles = sources.iter().map(|&i| self.triangles[i]).collect();
        (out.corner_parameters, out.corner_normals) = self.corners_of(&sources);
        (out, sources)
    }

    /// Append another mesh, shifting its indices.
    ///
    /// The corner arrays stay one entry per triangle: where only one of the
    /// two meshes has them, the other's are read from its vertices, and
    /// where that has none either they are dropped.
    pub fn append(&mut self, other: &Self) {
        #[allow(clippy::cast_possible_truncation)]
        let offset = self.positions.len() as u32;
        // Kept as they are where both meshes hold them per corner, as every
        // face mesh does, so appending face after face stays linear.
        let joined_p = if self.corner_parameters.len() == self.triangles.len()
            && other.corner_parameters.len() == other.triangles.len()
        {
            let mut a = std::mem::take(&mut self.corner_parameters);
            a.extend_from_slice(&other.corner_parameters);
            a
        } else {
            match (self.all_corner_parameters(), other.all_corner_parameters()) {
                (Some(mut a), Some(b)) => {
                    a.extend(b);
                    a
                }
                _ => Vec::new(),
            }
        };
        let joined_n = if self.corner_normals.len() == self.triangles.len()
            && other.corner_normals.len() == other.triangles.len()
        {
            let mut a = std::mem::take(&mut self.corner_normals);
            a.extend_from_slice(&other.corner_normals);
            a
        } else {
            match (self.all_corner_normals(), other.all_corner_normals()) {
                (Some(mut a), Some(b)) => {
                    a.extend(b);
                    a
                }
                _ => Vec::new(),
            }
        };
        self.positions.extend_from_slice(&other.positions);
        self.normals.extend_from_slice(&other.normals);
        #[allow(deprecated, reason = "kept filled for the release it is deprecated in")]
        self.parameters.extend_from_slice(&other.parameters);
        self.triangles
            .extend(other.triangles.iter().map(|t| t.map(|i| i + offset)));
        self.corner_parameters = joined_p;
        self.corner_normals = joined_n;
        self.deflection_met &= other.deflection_met;
        self.findings.extend_from_slice(&other.findings);
    }

    /// Merge vertices that coincide within `tol`, rewiring the triangles.
    ///
    /// Faces are triangulated independently, so a shared edge produces two
    /// copies of every boundary vertex, at identical positions, since both
    /// came from the same edge discretization, but as separate entries. Merging
    /// them is what turns a pile of face meshes into one closed surface, and
    /// what lets [`Triangulation::is_closed`] answer truthfully.
    #[must_use]
    pub fn welded(&self, tol: Tolerances) -> Self {
        self.welded_with_sources(tol).0
    }

    /// As [`Self::welded`], with the input triangle each output triangle is.
    #[must_use]
    pub fn welded_with_sources(&self, tol: Tolerances) -> (Self, Vec<usize>) {
        use ogeom_core::FastMap;

        // Quantize to a grid a good deal finer than the tolerance, then check
        // the neighbourhood: hashing alone would separate two points that
        // straddle a cell boundary however close they are.
        let cell = tol.confusion().max(f64::MIN_POSITIVE);
        let key = |p: Point| {
            #[allow(clippy::cast_possible_truncation)]
            (
                (p.x / cell).round() as i64,
                (p.y / cell).round() as i64,
                (p.z / cell).round() as i64,
            )
        };

        // Each cell's first and last vertex, and each vertex's successor in
        // its cell: the cells read back in the order they were filled.
        let mut buckets: FastMap<(i64, i64, i64), (u32, u32)> = FastMap::default();
        let mut next_in_cell: Vec<u32> = Vec::with_capacity(self.positions.len());
        let mut remap = vec![0_u32; self.positions.len()];
        let mut out = Self::new();
        out.deflection_met = self.deflection_met;
        out.findings.clone_from(&self.findings);

        for (index, position) in self.positions.iter().enumerate() {
            let (kx, ky, kz) = key(*position);
            let mut found = None;
            'search: for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let Some(&(first, _)) = buckets.get(&(
                            kx.saturating_add(dx),
                            ky.saturating_add(dy),
                            kz.saturating_add(dz),
                        )) else {
                            continue;
                        };
                        let mut candidate = first;
                        while candidate != u32::MAX {
                            if out.positions[candidate as usize].is_equal(*position, tol) {
                                found = Some(candidate);
                                break 'search;
                            }
                            candidate = next_in_cell[candidate as usize];
                        }
                    }
                }
            }

            let target = found.unwrap_or_else(|| {
                #[allow(clippy::cast_possible_truncation)]
                let fresh = out.positions.len() as u32;
                out.positions.push(*position);
                // A mesh read from a file may carry positions alone; the
                // optional per-vertex arrays follow only where they are full.
                if let Some(&normal) = self.normals.get(index)
                    && self.normals.len() == self.positions.len()
                {
                    out.normals.push(normal);
                }
                #[allow(deprecated, reason = "kept filled for the release it is deprecated in")]
                if let Some(&uv) = self.parameters.get(index)
                    && self.parameters.len() == self.positions.len()
                {
                    out.parameters.push(uv);
                }
                next_in_cell.push(u32::MAX);
                if let Some((_, last)) = buckets.get_mut(&(kx, ky, kz)) {
                    next_in_cell[*last as usize] = fresh;
                    *last = fresh;
                } else {
                    buckets.insert((kx, ky, kz), (fresh, fresh));
                }
                fresh
            });
            remap[index] = target;
        }

        let mut sources = Vec::new();
        for (source, t) in self.triangles.iter().enumerate() {
            let mapped = t.map(|i| remap[i as usize]);
            // A triangle whose corners merged is degenerate and contributes
            // nothing but trouble to anything that divides by its area.
            if mapped[0] != mapped[1] && mapped[1] != mapped[2] && mapped[2] != mapped[0] {
                out.triangles.push(mapped);
                sources.push(source);
            }
        }
        (out.corner_parameters, out.corner_normals) = self.corners_of(&sources);
        (out, sources)
    }
}

/// Assign each candidate once per incident side, with stable ties at narrow
/// corners. Hash lookups keep a densely subdivided side linear after sorting.
fn nearest_split_sides<V: Copy + Eq + std::hash::Hash>(
    corners: &[V; 3],
    sides: [&[V]; 3],
    mut distance: impl FnMut(usize, V) -> f64,
) -> ogeom_core::FastMap<V, usize> {
    let mut nearest: ogeom_core::FastMap<V, (usize, f64)> = ogeom_core::FastMap::default();
    for (side, vertices) in sides.into_iter().enumerate() {
        for &v in vertices {
            if corners.contains(&v) {
                continue;
            }
            let d = distance(side, v);
            let entry = nearest.entry(v).or_insert((side, d));
            if d.total_cmp(&entry.1).is_lt() {
                *entry = (side, d);
            }
        }
    }
    nearest
        .into_iter()
        .map(|(v, (side, _))| (v, side))
        .collect()
}

/// A mesh's two per-corner arrays: parameters, then normals.
type Corners = (Vec<[(f64, f64); 3]>, Vec<[Vector; 3]>);

/// The entries of a per-triangle array `all` for the triangles `sources`,
/// or none where `all` does not hold one per triangle of `count`.
fn picked<T: Copy>(all: &[T], count: usize, sources: &[usize]) -> Vec<T> {
    if all.is_empty() || all.len() != count {
        return Vec::new();
    }
    sources.iter().map(|&i| all[i]).collect()
}

/// How many triangles use each undirected edge, keyed (smaller, larger).
fn edge_uses(triangles: &[[u32; 3]]) -> ogeom_core::FastMap<(u32, u32), usize> {
    let mut uses: ogeom_core::FastMap<(u32, u32), usize> = ogeom_core::FastMap::default();
    // A closed mesh has three edges to every two triangles; sized for that,
    // the table is never grown while it is filled.
    uses.reserve(triangles.len() * 3 / 2 + 3);
    for t in triangles {
        for i in 0..3 {
            let (a, b) = (t[i], t[(i + 1) % 3]);
            *uses.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    uses
}

/// Some of a mesh's vertices bucketed by a cubic grid, for the ones near a
/// segment.
struct PointGrid {
    cell: f64,
    /// Each occupied cell's vertices.
    cells: ogeom_core::FastMap<(i64, i64, i64), Vec<u32>>,
    /// The vertices, for a query whose box spans more cells than there are
    /// vertices to look at.
    all: Vec<u32>,
}

impl PointGrid {
    /// `vertices` of `positions` in cells of about `spacing`, never finer
    /// than `reach`; `None` where no finite cell size fits.
    fn over(positions: &[Point], vertices: &[u32], spacing: f64, reach: f64) -> Option<Self> {
        let cell = spacing.max(reach);
        if !cell.is_finite() || cell <= 0.0 {
            return None;
        }
        let mut cells: ogeom_core::FastMap<(i64, i64, i64), Vec<u32>> =
            ogeom_core::FastMap::default();
        for &v in vertices {
            cells
                .entry(Self::key(positions[v as usize], cell))
                .or_default()
                .push(v);
        }
        Some(Self {
            cell,
            cells,
            all: vertices.to_vec(),
        })
    }

    #[allow(clippy::cast_possible_truncation)]
    fn key(p: Point, cell: f64) -> (i64, i64, i64) {
        (
            (p.x / cell).floor() as i64,
            (p.y / cell).floor() as i64,
            (p.z / cell).floor() as i64,
        )
    }

    /// Every vertex that may lie within `reach` of the segment from `a` to
    /// `b`, and some that do not: those in the cells round the segment's
    /// box widened by `reach`, with a cell to spare each way for rounding.
    fn within(&self, a: Point, b: Point, reach: f64, out: &mut Vec<u32>) {
        let low = Point::new(
            a.x.min(b.x) - reach,
            a.y.min(b.y) - reach,
            a.z.min(b.z) - reach,
        );
        let high = Point::new(
            a.x.max(b.x) + reach,
            a.y.max(b.y) + reach,
            a.z.max(b.z) + reach,
        );
        let (lo, hi) = (Self::key(low, self.cell), Self::key(high, self.cell));
        let span = |l: i64, h: i64| h.saturating_sub(l).saturating_add(3).unsigned_abs();
        let boxed = span(lo.0, hi.0)
            .saturating_mul(span(lo.1, hi.1))
            .saturating_mul(span(lo.2, hi.2));
        if boxed > self.cells.len() as u64 {
            out.extend_from_slice(&self.all);
            return;
        }
        for x in lo.0.saturating_sub(1)..=hi.0.saturating_add(1) {
            for y in lo.1.saturating_sub(1)..=hi.1.saturating_add(1) {
                for z in lo.2.saturating_sub(1)..=hi.2.saturating_add(1) {
                    if let Some(found) = self.cells.get(&(x, y, z)) {
                        out.extend_from_slice(found);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn a_mesh_of_positions_alone_welds() {
        let mut mesh = Triangulation::new();
        mesh.positions = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        mesh.triangles = vec![[0, 1, 2], [3, 4, 5]];
        let welded = mesh.welded(T);
        assert_eq!(welded.positions.len(), 4);
        assert!(welded.normals.is_empty() && !welded.has_parameters());
        assert_eq!(welded.triangles, vec![[0, 1, 2], [1, 3, 2]]);
    }

    #[test]
    fn an_empty_mesh_answers_sensibly() {
        let mesh = Triangulation::new();
        assert!(mesh.is_empty());
        assert_eq!(mesh.triangle_count(), 0);
        assert_relative_eq!(mesh.area(), 0.0);
        assert_relative_eq!(mesh.volume(), 0.0);
        assert!(!mesh.is_closed(), "nothing is not closed");
        assert!(mesh.bounds().is_empty());
    }

    /// A closed mesh encloses the same volume wherever it stands. A unit
    /// cube up to 1e8 from the origin has vertices exact in `f64` (spacing
    /// there is 1.5e-8), so any error is the sum's own rounding.
    #[test]
    fn a_closed_mesh_far_from_the_origin_encloses_its_own_volume() {
        for offset in [0.0, 1.0e6, 1.0e7, 1.0e8] {
            let o = Vector::new(offset, -offset, offset);
            let mut mesh = Triangulation::new();
            for k in 0..8 {
                let corner =
                    Point::new(f64::from(k & 1), f64::from((k >> 1) & 1), f64::from(k >> 2));
                mesh.positions.push(corner + o);
            }
            // Corner k sits at (k & 1, k >> 1 & 1, k >> 2); two triangles a
            // side, each wound counter-clockwise seen from outside.
            mesh.triangles = vec![
                [0, 2, 1],
                [1, 2, 3], // z = 0
                [4, 5, 6],
                [5, 7, 6], // z = 1
                [0, 1, 4],
                [1, 5, 4], // y = 0
                [2, 6, 3],
                [3, 6, 7], // y = 1
                [0, 4, 2],
                [2, 4, 6], // x = 0
                [1, 3, 5],
                [3, 7, 5], // x = 1
            ];
            assert!(mesh.is_closed(), "the cube closes at {offset:e}");
            // The cube's own rounding is nil, so 1e-9 leaves room only for
            // summing twelve terms of order one.
            assert_relative_eq!(mesh.volume(), 1.0, epsilon = 1e-9);
        }
    }

    /// A mesh of four corner positions, with whatever triangles are given.
    fn over(points: &[Point], triangles: &[[u32; 3]]) -> Triangulation {
        let mut mesh = Triangulation::new();
        for p in points {
            mesh.positions.push(*p);
            mesh.normals.push(Vector::Z);
        }
        mesh.triangles.extend_from_slice(triangles);
        mesh
    }

    /// Closure is a question about direction, not about how many.
    ///
    /// Counting an edge's uses answers the wrong question twice over: it
    /// calls a mesh with one face inside out closed, and it calls a mesh
    /// with four triangles round one edge open. Neither is what the
    /// divergence theorem asks, which is only that the surface have no
    /// boundary and wind one way.
    #[test]
    fn a_closed_mesh_is_one_crossed_as_often_each_way() {
        let corners = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
        ];
        // A tetrahedron, every face wound outward.
        let solid = over(&corners, &[[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]);
        assert!(solid.is_closed(), "a tetrahedron closes");
        assert!(solid.volume() > 0.0, "and wound outward");

        // One face turned over. Every edge is still used exactly twice, so
        // counting uses calls this closed; two of its edges are now crossed
        // the same way twice, and the volume it gives is wrong.
        let mut flipped = solid.clone();
        flipped.triangles[3] = [1, 3, 2];
        assert!(
            !flipped.is_closed(),
            "a face inside out is not a closed mesh"
        );

        // A hole: one face dropped. Three edges are crossed once.
        let mut holed = solid.clone();
        holed.triangles.pop();
        assert!(!holed.is_closed(), "three edges left dangling");

        // Two tetrahedra sharing the edge 0-1, each closed and outward. The
        // shared edge carries four triangles, two crossing each way. It is
        // not a manifold and it is certainly closed, and its volume is both
        // halves, which is the case a count of uses refuses and the one a
        // real assembly produces where two faces fill a sliver with the same
        // chord.
        let mut pair = solid.clone();
        let mirrored = Point::new(0.0, -1.0, 0.0);
        #[allow(clippy::cast_possible_truncation)]
        let m = pair.positions.len() as u32;
        pair.positions.push(mirrored);
        pair.normals.push(Vector::Z);
        let apex = 3;
        pair.triangles
            .extend_from_slice(&[[0, 1, m], [0, m, apex], [0, apex, 1], [1, apex, m]]);
        let four = pair
            .triangles
            .iter()
            .flat_map(|t| (0..3).map(move |i| (t[i], t[(i + 1) % 3])))
            .filter(|(a, b)| (*a == 0 && *b == 1) || (*a == 1 && *b == 0))
            .count();
        assert_eq!(four, 4, "the shared edge carries four triangles");
        assert!(pair.is_closed(), "and the pair is still closed");
    }

    /// A mesh built with one `(u, v)` per vertex alone reads its corners
    /// from them, and keeps doing so through a weld.
    #[test]
    #[allow(deprecated, reason = "the per-vertex field is what this reads back")]
    fn corners_fall_back_to_per_vertex_parameters() {
        let mut mesh = over(
            &[
                Point::ORIGIN,
                Point::new(1.0, 0.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
            ],
            &[[0, 1, 2]],
        );
        mesh.parameters = vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
        assert!(mesh.has_parameters());
        let expected = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
        assert_eq!(mesh.triangle_parameters(0), Some(expected));
        assert_eq!(mesh.welded(T).triangle_parameters(0), Some(expected));
    }

    /// A triangle split where a neighbour's vertex stands on its side gives
    /// the new corner the parameters part way along that side, and keeps
    /// the old corners' own: the corners at one vertex need not agree.
    #[test]
    fn a_stitched_side_interpolates_its_corners() {
        let points = [
            Point::ORIGIN,
            Point::new(2.0, 0.0, 0.0),
            Point::new(0.0, 2.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, -1.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
        ];
        // The upper triangle's long side runs under the two lower ones,
        // which meet at the middle of it.
        let mut mesh = over(&points, &[[0, 1, 2], [0, 4, 3], [3, 4, 5]]);
        mesh.corner_parameters = vec![
            [(10.0, 0.0), (12.0, 0.0), (10.0, 2.0)],
            [(0.0, 0.0), (1.0, -1.0), (1.0, 0.0)],
            [(1.0, 0.0), (1.0, -1.0), (2.0, 0.0)],
        ];
        mesh.set_corner_normals_from_vertices();
        let welded = mesh.border_welded(1e-6);
        let stitched = welded.border_stitched(1e-6);
        assert_eq!(stitched.triangles.len(), 4);
        assert_eq!(stitched.corner_parameters.len(), 4);
        let mut on_the_side = 0;
        for (t, uv) in stitched.triangles.iter().zip(&stitched.corner_parameters) {
            for (&v, &c) in t.iter().zip(uv) {
                if stitched.positions[v as usize] == Point::new(1.0, 0.0, 0.0) && c.0 > 5.0 {
                    assert!((c.0 - 11.0).abs() < 1e-12 && c.1.abs() < 1e-12, "{c:?}");
                    on_the_side += 1;
                }
            }
        }
        assert_eq!(on_the_side, 2, "both pieces of the upper triangle hold it");
    }

    #[test]
    fn welding_drops_triangles_that_collapse() {
        // Three corners that merge into one describe no area, and anything that
        // divides by a triangle's area would divide by zero.
        let mut mesh = Triangulation::new();
        for _ in 0..3 {
            mesh.positions.push(Point::ORIGIN);
            mesh.normals.push(Vector::Z);
        }
        mesh.triangles.push([0, 1, 2]);
        let welded = mesh.welded(T);
        assert_eq!(welded.vertex_count(), 1);
        assert_eq!(welded.triangle_count(), 0);
    }

    #[test]
    fn appending_shifts_indices_rather_than_overlapping_them() {
        let mut a = Triangulation::new();
        a.positions.push(Point::ORIGIN);
        a.normals.push(Vector::Z);

        let mut b = Triangulation::new();
        b.positions.push(Point::new(1.0, 0.0, 0.0));
        b.normals.push(Vector::Z);
        b.triangles.push([0, 0, 0]);

        a.append(&b);
        assert_eq!(a.vertex_count(), 2);
        assert_eq!(
            a.triangles[0],
            [1, 1, 1],
            "b's index moved past a's vertices"
        );
    }

    #[test]
    fn deflection_failure_propagates_through_the_whole_mesh() {
        // One face that could not meet its tolerance makes the whole mesh's
        // claim untrue, and the flag has to say so rather than being averaged
        // away.
        let mut a = Triangulation::new();
        let mut b = Triangulation::new();
        b.deflection_met = false;
        a.append(&b);
        assert!(!a.deflection_met);
    }

    /// A tetrahedron whose one face is a sliver a hundredth wide, dropped:
    /// sealing at a width past the sliver's closes it again, and at a width
    /// under it leaves it open. A triangle laid twice facing opposite ways
    /// encloses nothing and is cancelled.
    #[test]
    fn sealing_closes_cracks_narrower_than_asked_and_cancels_folds() {
        let corners = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.5, 0.01, 0.0),
            Point::new(0.5, 0.5, 1.0),
        ];
        let solid = over(&corners, &[[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]]);
        assert!(solid.is_closed());
        assert!(solid.volume() > 0.0);
        let mut cracked = solid.clone();
        cracked.triangles.remove(0);
        assert!(!cracked.is_closed());
        let sealed = cracked.sealed(0.05);
        assert!(sealed.is_closed(), "a crack under the width is sealed");
        assert_relative_eq!(sealed.volume(), solid.volume(), epsilon = 1e-12);
        assert!(
            !cracked.sealed(1e-3).is_closed(),
            "a crack wider than asked stays open"
        );
        let mut folded = solid.clone();
        folded.triangles.push([0, 1, 3]);
        folded.triangles.push([0, 3, 1]);
        let unfolded = folded.sealed(0.0);
        assert_eq!(
            unfolded.triangle_count(),
            solid.triangle_count(),
            "the fold is cancelled"
        );
        assert!(unfolded.is_closed());
        let (unfolded, sources) = folded.without_folds_with_sources();
        assert_eq!(sources.len(), unfolded.triangles.len());
        for (triangle, source) in unfolded.triangles.iter().zip(sources) {
            assert_eq!(*triangle, folded.triangles[source]);
        }
    }

    /// A triangle sealing a crack is named after the triangle across the
    /// crack side it closes: that triangle walks one of its sides the other
    /// way. Every other triangle is named after itself.
    #[test]
    fn a_sealing_triangle_is_named_after_its_neighbour_across_the_crack() {
        let corners = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.5, 0.01, 0.0),
            Point::new(0.5, 0.5, 1.0),
        ];
        let mut cracked = over(&corners, &[[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]]);
        cracked.triangles.remove(0);
        let (sealed, sources) = cracked.sealed_with_sources(0.05);
        assert_eq!(sources.len(), sealed.triangles.len());
        assert_eq!(sealed.triangles.len(), cracked.triangles.len() + 1);
        for (i, (triangle, &source)) in sealed.triangles.iter().zip(&sources).enumerate() {
            if i < cracked.triangles.len() {
                assert_eq!(*triangle, cracked.triangles[source]);
                continue;
            }
            let across = cracked.triangles[source];
            let shares_a_side_reversed = (0..3).any(|k| {
                let (a, b) = (triangle[k], triangle[(k + 1) % 3]);
                (0..3).any(|j| (across[j], across[(j + 1) % 3]) == (b, a))
            });
            assert!(
                shares_a_side_reversed,
                "{triangle:?} named after {across:?}"
            );
        }
    }

    #[test]
    fn dense_border_assigns_each_candidate_once_per_side() {
        use std::cell::Cell;
        use std::hash::{Hash, Hasher};

        #[derive(Clone, Copy)]
        struct Vertex<'a>(u32, &'a Cell<usize>);
        impl PartialEq for Vertex<'_> {
            fn eq(&self, other: &Self) -> bool {
                self.1.set(self.1.get() + 1);
                self.0 == other.0
            }
        }
        impl Eq for Vertex<'_> {}
        impl Hash for Vertex<'_> {
            fn hash<H: Hasher>(&self, state: &mut H) {
                self.0.hash(state);
            }
        }

        let comparisons = Cell::new(0);
        let vertices: Vec<_> = (3..10_003).map(|v| Vertex(v, &comparisons)).collect();
        let corners = [
            Vertex(0, &comparisons),
            Vertex(1, &comparisons),
            Vertex(2, &comparisons),
        ];
        let mut visits = 0;
        let nearest =
            nearest_split_sides(&corners, [&vertices, &vertices[..20], &[]], |side, v| {
                visits += 1;
                // Twenty vertices also reach the adjacent side, which is nearer.
                if side == 1 && v.0 < 23 { 0.0 } else { 1.0 }
            });
        // Three corner checks per candidate and hash-table probes fit this
        // linear budget. A split-list scan alone takes over 50 million checks.
        assert!(
            comparisons.get() < 8 * vertices.len(),
            "{} comparisons",
            comparisons.get()
        );
        assert_eq!(visits, vertices.len() + 20);
        assert_eq!(nearest.len(), vertices.len());
        for v in vertices {
            assert_eq!(nearest[&v], usize::from(v.0 < 23));
        }
    }

    /// A split on any of the three sides must create two nonzero triangles,
    /// including the sides incident to the first corner.
    #[test]
    fn stitching_splits_each_side_and_retains_triangle_sources() {
        for [a, b, c] in [[0, 1, 2], [1, 2, 0], [2, 0, 1]] {
            let mut mesh = over(
                &[
                    Point::new(0.0, 0.0, 0.0),
                    Point::new(2.0, 0.0, 0.0),
                    Point::new(0.0, 2.0, 0.0),
                ],
                &[[0, 1, 2]],
            );
            mesh.positions.push(Point::midpoint(
                mesh.positions[a as usize],
                mesh.positions[b as usize],
            ));
            // The neighbour is in another plane, already split at position 3.
            mesh.positions.push(Point::new(0.0, 0.0, 2.0));
            mesh.triangles.extend([[b, 3, 4], [3, a, 4]]);
            let (split, sources) = mesh.border_stitched_with_sources(1e-9);
            assert_eq!(split.triangles.len(), 4);
            assert_eq!(sources.iter().filter(|&&i| i == 0).count(), 2);
            for (t, source) in split.triangles.iter().zip(sources) {
                let [pa, pb, pc] = t.map(|i| split.positions[i as usize]);
                assert!((pb - pa).cross(pc - pa).magnitude() > 0.0);
                if source == 0 {
                    assert!(t.contains(&3) && t.contains(&c));
                    assert!((pb - pa).cross(pc - pa).z > 0.0);
                }
            }
        }
    }

    /// The sliver crack of the tetrahedron above, moved up to 2e7 from the
    /// origin along a tilted direction: a crack's width is its own, so the
    /// same widths seal it and leave it open there as at the origin.
    #[test]
    fn sealing_far_from_the_origin_reads_a_crack_by_its_own_width() {
        for offset in [0.0, 1.0e6, 5.0e6, 1.0e7, 2.0e7] {
            let o = Vector::new(offset, -offset, offset);
            let corners = [
                Point::new(0.0, 0.0, 0.0) + o,
                Point::new(1.0, 0.0, 0.0) + o,
                Point::new(0.5, 0.01, 0.0) + o,
                Point::new(0.5, 0.5, 1.0) + o,
            ];
            let cracked = over(&corners, &[[0, 1, 3], [1, 2, 3], [2, 0, 3]]);
            assert!(
                cracked.sealed(0.05).is_closed(),
                "a crack under the width is sealed at {offset:e}"
            );
            assert!(
                !cracked.sealed(1e-3).is_closed(),
                "a crack wider than asked stays open at {offset:e}"
            );
        }
    }

    /// Two strips meeting along a line, one drawn to twice the other's
    /// spacing there: every segment of the coarse side has a vertex of the
    /// fine side at its middle, and stitching splits each at it.
    #[test]
    fn a_coarse_border_is_split_at_its_neighbours_vertices() {
        const N: u32 = 60;
        let mut mesh = Triangulation::new();
        // The fine strip: a bottom row and a shared row at y = 1.
        for y in [0.0, 1.0] {
            for x in 0..=N {
                mesh.positions.push(Point::new(f64::from(x), y, 0.0));
            }
        }
        let top = |x: u32| N + 1 + x;
        for x in 0..N {
            mesh.triangles.push([x, x + 1, top(x + 1)]);
            mesh.triangles.push([x, top(x + 1), top(x)]);
        }
        // The coarse strip above, on every other vertex of the shared row.
        let apex = top(N) + 1;
        for x in (0..=N).step_by(2) {
            mesh.positions.push(Point::new(f64::from(x), 2.0, 0.0));
        }
        for (i, x) in (0..N).step_by(2).enumerate() {
            let i = u32::try_from(i).unwrap();
            let (a, b) = (apex + i, apex + i + 1);
            mesh.triangles.push([top(x), top(x + 2), b]);
            mesh.triangles.push([top(x), b, a]);
        }

        let stitched = mesh.border_stitched(1e-6);
        let uses = edge_uses(&stitched.triangles);
        for x in 0..N {
            assert_eq!(
                uses.get(&(top(x), top(x + 1))),
                Some(&2),
                "the shared row's segment at {x} is used from both sides"
            );
        }
        assert_eq!(stitched.positions, mesh.positions, "no vertex moves");
    }
}
