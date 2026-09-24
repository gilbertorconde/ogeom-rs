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
/// Vertices carry their parameters as well as their positions, so a caller can
/// ask the exact surface about a triangulated point rather than only the
/// approximation.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Triangulation {
    /// Vertex positions.
    pub positions: Vec<Point>,
    /// Outward unit normals, one per vertex.
    pub normals: Vec<Vector>,
    /// The surface parameters each vertex came from.
    pub parameters: Vec<(f64, f64)>,
    /// Triangles, as indices into the vertex arrays, wound counter-clockwise
    /// about the outward normal.
    pub triangles: Vec<[u32; 3]>,
    /// Whether every face met its requested deflection.
    pub deflection_met: bool,
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
    /// it forms with the origin, and the contributions cancel except over the
    /// enclosed region. An open mesh gives a number with no meaning, and a mesh
    /// wound inward gives the negative, which is why
    /// [`Triangulation::is_closed`] exists to be asked first.
    #[must_use]
    pub fn volume(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| self.positions[i as usize].to_vector());
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
    /// crossing each way. There is no hole there and the volume is right;
    /// demanding exactly two refuses a mesh for being non-manifold when
    /// nothing was asked about manifoldness. Sixty-four bodies of one real
    /// assembly were refused that way, forty-four of them for this alone.
    ///
    /// This also agrees with the topology side at last:
    /// [`is_shell_closed`](../../ogeom_algo/fn.is_shell_closed.html) counts an
    /// edge's uses and accepts any even number, and the two halves of the
    /// kernel should not mean different things by the same word.
    #[must_use]
    pub fn is_closed(&self) -> bool {
        use std::collections::HashMap;
        let mut balance: HashMap<(u32, u32), i64> = HashMap::new();
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
        use std::collections::HashMap;
        if !reach.is_finite() || reach <= 0.0 {
            return self.clone();
        }
        // Border vertices: endpoints of triangle edges used an odd number of
        // times.
        let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
        for t in &self.triangles {
            for i in 0..3 {
                let (a, b) = (t[i], t[(i + 1) % 3]);
                *uses.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        let mut border: Vec<u32> = uses
            .iter()
            .filter(|&(_, &n)| n % 2 == 1)
            .flat_map(|(&(a, b), _)| [a, b])
            .collect();
        border.sort_unstable();
        border.dedup();
        if border.is_empty() {
            return self.clone();
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
        let mut buckets: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::new();
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
                            .get(&(kx + dx, ky + dy, kz + dz))
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
        // Compact: keep every vertex that survives as its own representative
        // or is referenced; simplest is to keep all and let triangles remap.
        out.positions = self.positions.clone();
        out.normals = self.normals.clone();
        out.parameters = self.parameters.clone();
        for t in &self.triangles {
            let mapped = t.map(|i| remap[i as usize]);
            if mapped[0] != mapped[1] && mapped[1] != mapped[2] && mapped[2] != mapped[0] {
                out.triangles.push(mapped);
            }
        }
        out
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
        use std::collections::HashMap;
        if !reach.is_finite() || reach <= 0.0 {
            return self.clone();
        }
        let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
        for t in &self.triangles {
            for i in 0..3 {
                let (a, b) = (t[i], t[(i + 1) % 3]);
                *uses.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        let border_edges: Vec<(u32, u32)> = uses
            .iter()
            .filter(|&(_, &n)| n % 2 == 1)
            .map(|(&e, _)| e)
            .collect();
        if border_edges.is_empty() {
            return self.clone();
        }
        let mut border_vertices: Vec<u32> =
            border_edges.iter().flat_map(|&(a, b)| [a, b]).collect();
        border_vertices.sort_unstable();
        border_vertices.dedup();

        // For every border segment, the border vertices sitting on its
        // interior, ordered along it.
        let mut splits: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
        for &(a, b) in &border_edges {
            let (pa, pb) = (self.positions[a as usize], self.positions[b as usize]);
            let d = pb - pa;
            let l2 = d.dot(d);
            if l2 <= 0.0 {
                continue;
            }
            let mut on: Vec<(f64, u32)> = border_vertices
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
            on.sort_by(|x, y| x.0.total_cmp(&y.0));
            splits.insert((a, b), on.into_iter().map(|(_, v)| v).collect());
        }
        if splits.is_empty() {
            return self.clone();
        }

        let mut out = Self::new();
        out.deflection_met = self.deflection_met;
        out.positions = self.positions.clone();
        out.normals = self.normals.clone();
        out.parameters = self.parameters.clone();
        for t in &self.triangles {
            // The triangle's ring with any split points inserted, fanned from
            // its first corner.
            let mut ring: Vec<u32> = Vec::with_capacity(6);
            let mut any = false;
            for i in 0..3 {
                let (a, b) = (t[i], t[(i + 1) % 3]);
                ring.push(a);
                if let Some(vs) = splits.get(&(a.min(b), a.max(b))) {
                    any = true;
                    if a < b {
                        ring.extend(vs.iter().copied());
                    } else {
                        ring.extend(vs.iter().rev().copied());
                    }
                }
            }
            if !any {
                out.triangles.push(*t);
                continue;
            }
            for i in 1..ring.len() - 1 {
                let tri = [ring[0], ring[i], ring[i + 1]];
                if tri[0] != tri[1] && tri[1] != tri[2] && tri[2] != tri[0] {
                    out.triangles.push(tri);
                }
            }
        }
        out
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
        use std::collections::HashMap;
        let mut out = self.clone();
        // Folds.
        let mut seen: HashMap<[u32; 3], Vec<usize>> = HashMap::new();
        for (i, t) in out.triangles.iter().enumerate() {
            let mut key = *t;
            key.sort_unstable();
            seen.entry(key).or_default().push(i);
        }
        let mut drop = vec![false; out.triangles.len()];
        for list in seen.values() {
            let mut open: Vec<usize> = Vec::new();
            for &i in list {
                let t = out.triangles[i];
                let reverse = open.iter().position(|&j| {
                    let u = out.triangles[j];
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
        let mut index = 0;
        out.triangles.retain(|_| {
            let keep = !drop[index];
            index += 1;
            keep
        });
        if !width.is_finite() || width <= 0.0 {
            return out;
        }
        // Cracks: each border edge walked the other way from its triangle.
        let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
        for t in &out.triangles {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                *uses.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        let mut onward: HashMap<u32, Vec<u32>> = HashMap::new();
        for t in &out.triangles {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                if uses[&(a.min(b), a.max(b))] == 1 {
                    onward.entry(b).or_default().push(a);
                }
            }
        }
        let mut done: std::collections::HashSet<u32> = std::collections::HashSet::new();
        let mut starts: Vec<u32> = onward.keys().copied().collect();
        starts.sort_unstable();
        for start in starts {
            if done.contains(&start) {
                continue;
            }
            let mut ring = vec![start];
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
                if ring.contains(&next) || ring.len() > onward.len() {
                    break false;
                }
                ring.push(next);
                at = next;
            };
            for &v in &ring {
                done.insert(v);
            }
            if !closed || ring.len() < 3 {
                continue;
            }
            let points: Vec<Point> = ring.iter().map(|&v| out.positions[v as usize]).collect();
            let mut normal = Vector::ZERO;
            let mut perimeter = 0.0;
            for (i, p) in points.iter().enumerate() {
                let q = points[(i + 1) % points.len()];
                normal += p.to_vector().cross(q.to_vector());
                perimeter += p.distance(q);
            }
            let area = normal.magnitude() / 2.0;
            if perimeter <= 0.0 || 2.0 * area / perimeter > width {
                continue;
            }
            for i in 1..ring.len() - 1 {
                out.triangles.push([ring[0], ring[i], ring[i + 1]]);
            }
        }
        out
    }

    /// Append another mesh, shifting its indices.
    pub fn append(&mut self, other: &Self) {
        #[allow(clippy::cast_possible_truncation)]
        let offset = self.positions.len() as u32;
        self.positions.extend_from_slice(&other.positions);
        self.normals.extend_from_slice(&other.normals);
        self.parameters.extend_from_slice(&other.parameters);
        self.triangles
            .extend(other.triangles.iter().map(|t| t.map(|i| i + offset)));
        self.deflection_met &= other.deflection_met;
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
        use std::collections::HashMap;

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

        let mut buckets: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::new();
        let mut remap = vec![0_u32; self.positions.len()];
        let mut out = Self::new();
        out.deflection_met = self.deflection_met;

        for (index, position) in self.positions.iter().enumerate() {
            let (kx, ky, kz) = key(*position);
            let mut found = None;
            'search: for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        for &candidate in buckets
                            .get(&(kx + dx, ky + dy, kz + dz))
                            .map_or(&[][..], Vec::as_slice)
                        {
                            if out.positions[candidate as usize].is_equal(*position, tol) {
                                found = Some(candidate);
                                break 'search;
                            }
                        }
                    }
                }
            }

            let target = found.unwrap_or_else(|| {
                #[allow(clippy::cast_possible_truncation)]
                let fresh = out.positions.len() as u32;
                out.positions.push(*position);
                out.normals.push(self.normals[index]);
                out.parameters.push(self.parameters[index]);
                buckets.entry((kx, ky, kz)).or_default().push(fresh);
                fresh
            });
            remap[index] = target;
        }

        for t in &self.triangles {
            let mapped = t.map(|i| remap[i as usize]);
            // A triangle whose corners merged is degenerate and contributes
            // nothing but trouble to anything that divides by its area.
            if mapped[0] != mapped[1] && mapped[1] != mapped[2] && mapped[2] != mapped[0] {
                out.triangles.push(mapped);
            }
        }
        out
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

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

    /// A mesh of four corner positions, with whatever triangles are given.
    fn over(points: &[Point], triangles: &[[u32; 3]]) -> Triangulation {
        let mut mesh = Triangulation::new();
        for p in points {
            mesh.positions.push(*p);
            mesh.normals.push(Vector::Z);
            mesh.parameters.push((0.0, 0.0));
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
        pair.parameters.push((0.0, 0.0));
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

    #[test]
    fn welding_drops_triangles_that_collapse() {
        // Three corners that merge into one describe no area, and anything that
        // divides by a triangle's area would divide by zero.
        let mut mesh = Triangulation::new();
        for _ in 0..3 {
            mesh.positions.push(Point::ORIGIN);
            mesh.normals.push(Vector::Z);
            mesh.parameters.push((0.0, 0.0));
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
        a.parameters.push((0.0, 0.0));

        let mut b = Triangulation::new();
        b.positions.push(Point::new(1.0, 0.0, 0.0));
        b.normals.push(Vector::Z);
        b.parameters.push((1.0, 0.0));
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
    }
}
