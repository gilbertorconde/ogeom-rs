//! Filling a rim of an editable mesh with a patch of new faces.
//!
//! # The rim, as loops
//!
//! A rim is named by its [`RimId`], its lowest-numbered half-edge. A rim
//! that passes through one point more than once, as where two holes touch
//! at a corner, is parted there into loops that each pass through every
//! point once, and each loop is filled on its own. The patches meet at the
//! point the rim passed more than once and close the fans about it into
//! one, so the point becomes a single vertex again.
//!
//! A loop of two edges, and a rim edge that more faces than the rim's own
//! use, cannot be filled into a manifold, and are refused.
//!
//! # The triangulation
//!
//! A loop of up to [`OPTIMAL_UP_TO`] edges is triangulated on its own
//! vertices, by exact dynamic programming over every triangulation. The
//! programme minimises first the *largest bend* and then the area. A bend
//! is the angle between the normals of the two faces that meet across an
//! edge, and the largest is taken over every edge of the patch, including
//! the rim's edges with the faces already on them. A face with collinear
//! corners has no normal, so the triangulation never makes one, and the
//! bend across a rim edge whose face has none counts as zero.
//!
//! A bend `β` is compared as `4 sin²(β/2)`, read from the sum of the unit
//! directions from the edge to the two far corners. That grows with the
//! bend, needs no angle, and keeps its precision for small bends.
//!
//! A bend depends on the faces on both sides of an edge, so the
//! programme's state is a triangle on an edge, not the edge alone. A first
//! pass finds the least largest bend, `θ`. A second pass finds the least
//! area among the triangulations that bend nowhere more than `θ`. Both
//! are exact over all triangulations of the loop: `θ` to the bit, as
//! computed, and the area to the rounding of its sum. Ties go to the
//! lowest-numbered apex. Memory grows as `n³/6` numbers for a loop of `n`
//! edges, about 170 MB at the bound. Time grows as `n³` while the cheapest
//! apexes on each edge fit the faces beyond it, as on smooth holes, and as
//! `n⁴` at worst.
//!
//! The triangulation draws no edge between two points the mesh already
//! joins. Such an edge would have four faces.
//!
//! A longer loop is flattened onto the plane across which it encloses the
//! most area, normal to its vector area, and ear-clipped there by exact
//! orientation tests. A loop whose shadow on that plane crosses or touches
//! itself is refused.
//!
//! # Refinement
//!
//! With [`Fill::refine`], the patch is refined to the spacing of the mesh
//! about it, by Liepa's scheme. Each rim vertex takes the mean length of
//! the mesh's edges at it as its spacing `σ`, and each new vertex the mean
//! of its triangle's corners'. A triangle is split at its centroid `c`,
//! into three, when `c` stands farther than `σ/√2` from each corner `m`
//! (for `σ` both `c`'s and `m`'s own). That is Liepa's density factor: it
//! leaves edges about as long as the rim's. The centroid must also stand
//! that far from every vertex of the patch, so that no two vertices crowd
//! each other. New vertices lie within the rim's convex hull, at least the
//! smallest rim spacing over `√2` apart, so the splitting ends.
//!
//! After each round of splits, each edge inside the patch is flipped while
//! that raises the smallest angle of its two faces. A flip is not made if
//! it would turn a face over against the other, or draw an edge the patch
//! or the mesh already has. Each flip raises the patch's sorted angles
//! lexicographically, so the flipping ends too.
//!
//! # Fairing
//!
//! Only the vertices the refinement made move. Without refinement the
//! patch is the triangulation itself.
//!
//! The Laplacian is the cotangent one, weighted by the patch as refined
//! and, at the rim, by the mesh's faces about it. Weights read from the
//! geometry, rather than one per neighbour, let the result converge as the
//! mesh is refined, however unevenly the mesh about the hole is spaced.
//!
//! - [`Fairing::Flat`] spans the rim tightly: a membrane, with the
//!   Laplacian vanishing at each new vertex and the rim held fixed. It is
//!   exact for linear functions, so a hole in a plane stays in the plane.
//! - [`Fairing::Curved`] carries the surface's curvature across the hole.
//!   It minimises the sum, over the patch's vertices and the rim's, of the
//!   squared Laplacian times each vertex's area (a third of its faces').
//!   The Laplacian at a rim vertex reaches one ring into the surrounding
//!   mesh, which is held fixed. This is the discrete thin-plate energy,
//!   and its stationary condition is a bi-Laplacian vanishing over the
//!   patch.
//!
//! The three coordinates are solved as one system, by conjugate gradients
//! on the normal equations. The fill therefore does not depend on the
//! axes: a hole in a plane is filled in that plane, to rounding. The
//! solve takes about as many steps as there are new vertices, so its time
//! grows as their square.

use super::{Editable, NONE, RimId, Scope, index, prev};
use core::f64::consts::FRAC_1_SQRT_2;
use core::ops::Range;
use ogeom_core::{Exact, FastMap, FastSet, OgeomResult, Predicates, Sign, ogeom_bail};
use ogeom_math::{Point, SparseMatrix, Vector, least_squares_cgnr};

/// The longest loop triangulated optimally, by exact dynamic programming;
/// a longer one is ear-clipped on its plane.
pub const OPTIMAL_UP_TO: usize = 500;

/// How the new vertices of a refined patch are placed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Fairing {
    /// Carrying the surface's curvature across the hole: the thin-plate
    /// patch, with one ring of the mesh about it held fixed.
    Curved,
    /// The tightest patch: a membrane, harmonic over the patch, with the
    /// rim held fixed.
    Flat,
}

/// How a rim is filled.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Fill {
    /// Where the vertices the refinement makes settle.
    pub fairing: Fairing,
    /// Whether to refine the patch to the spacing of the mesh about it.
    /// Without refinement the patch has only the rim's vertices, and
    /// nothing is faired.
    pub refine: bool,
}

/// What a fill made.
#[derive(Debug, Clone, PartialEq)]
pub struct FillReport {
    /// The faces made, numbered after every face there was.
    pub faces: Range<u32>,
    /// The vertices the refinement made, numbered after every vertex there
    /// was. Empty without refinement.
    pub vertices: Range<u32>,
    /// The loops the rim was parted into at the points it passes more than
    /// once, each filled on its own: one for a rim passing every point
    /// once.
    pub loops: usize,
    /// How many of the loops were longer than [`OPTIMAL_UP_TO`] edges, and
    /// were flattened and ear-clipped.
    pub flattened: usize,
    /// The largest bend, in radians, of the triangulation before
    /// refinement. A bend is the angle between the normals of two faces
    /// across an edge of the triangulation, rim edges included. For loops
    /// triangulated optimally, no triangulation of them bends less.
    pub largest_bend: f64,
}

impl Editable {
    /// Fill the rim `rim` with a patch of new faces, triangulated, refined
    /// and faired as `fill` asks. See the [module documentation](self) for
    /// how.
    ///
    /// The patch shares every rim edge with the face beyond it, so the mesh
    /// is closed about it. A mesh whose every rim is filled is closed. No
    /// vertex of the mesh moves. A rim passing through a point more than
    /// once is closed about that point, whose vertices become one: the
    /// lowest-numbered keeps its number. New faces and vertices take
    /// numbers after every one there is, and the new vertices have no
    /// [`source_vertex`]. Nothing changes when the fill is refused.
    ///
    /// [`source_vertex`]: Editable::source_vertex
    ///
    /// # Errors
    ///
    /// - [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if
    ///   `rim` names no rim: a half-edge that is not on a rim, or not the
    ///   lowest-numbered of its rim.
    /// - [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    ///   if the rim cannot be filled into a manifold: a loop of two edges,
    ///   or a rim edge more faces use. Also if a loop has no triangulation
    ///   without a face of collinear corners or an edge the mesh already
    ///   has, if the loops' triangulations would draw one edge twice, or if
    ///   a long loop's shadow on its plane crosses itself.
    /// - [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if ear
    ///   clipping finds no ear it may cut.
    /// - [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if the
    ///   fairing solve does not settle.
    /// - [`OgeomError::Range`](ogeom_core::OgeomError::Range) if the faces
    ///   and vertices made would pass `u32` numbering.
    pub fn fill(&mut self, rim: RimId, fill: &Fill) -> OgeomResult<FillReport> {
        let edges = self.rim_named(rim)?;
        let (mut patch, loops) = self.gather(&edges)?;
        let mut flattened = 0;
        let mut largest_bend = 0.0_f64;
        let mut drawn: FastSet<(u32, u32)> = FastSet::default();
        for lp in &loops {
            let ring = patch.ring(self, lp);
            let free = |i: usize, j: usize| {
                let (a, b) = (ring.ids[i], ring.ids[j]);
                !patch.joined.contains(&(a.min(b), a.max(b)))
            };
            let triangles = if lp.len() <= OPTIMAL_UP_TO {
                let Some(best) = optimal(&ring, &free) else {
                    ogeom_bail!(
                        Construction,
                        "rim {} has a loop of {} edges with no triangulation free of collinear \
                         corners and of edges the mesh already has",
                        rim.0,
                        lp.len()
                    );
                };
                best
            } else {
                flattened += 1;
                ear_clip(&ring, &free)?
            };
            largest_bend = largest_bend.max(radians(measure(&ring, &triangles).0));
            // Each new edge has a face on either side, wound opposite ways:
            // counted once, from the side running it upward.
            let n = ring.ids.len();
            for t in &triangles {
                for k in 0..3 {
                    let (i, j) = (t[k], t[(k + 1) % 3]);
                    if i > j || j - i == 1 || (i == 0 && j == n - 1) {
                        continue;
                    }
                    let (a, b) = (ring.ids[i], ring.ids[j]);
                    if !drawn.insert((a.min(b), a.max(b))) {
                        ogeom_bail!(
                            Construction,
                            "the loops of rim {} would draw one edge twice",
                            rim.0
                        );
                    }
                }
                patch.triangles.push(t.map(|k| ring.ids[k]));
            }
        }
        if fill.refine {
            patch.refine();
        }
        patch.fair(fill.fairing)?;
        let faces = self.corners.len().saturating_add(patch.triangles.len());
        let vertices = self.positions.len().saturating_add(patch.positions.len());
        if faces.saturating_mul(3) >= NONE as usize || vertices >= NONE as usize {
            ogeom_bail!(
                Range,
                "filling rim {} would make {faces} faces and {vertices} vertices, past u32 \
                 numbering",
                rim.0
            );
        }
        let (faces, vertices) = self.commit(&patch);
        Ok(FillReport {
            faces,
            vertices,
            loops: loops.len(),
            flattened,
            largest_bend,
        })
    }

    /// The half-edges of the rim `rim` names, in order from it.
    fn rim_named(&self, rim: RimId) -> OgeomResult<Vec<u32>> {
        let h = rim.0;
        let k = h as usize;
        if k >= self.twin.len() || !self.live[k / 3] || self.twin[k] != NONE {
            ogeom_bail!(Dangling, "half-edge {h} is not on a rim");
        }
        let mut edges = Vec::new();
        let mut g = h;
        loop {
            if g < h {
                ogeom_bail!(
                    Dangling,
                    "half-edge {h} is on the rim of half-edge {g}, which names it"
                );
            }
            edges.push(g);
            g = self.rim_after(g);
            if g == h {
                break;
            }
        }
        Ok(edges)
    }

    /// The patch's rim points, with what the mesh says about each, and the
    /// rim parted into loops at the points it passes more than once.
    fn gather(&self, edges: &[u32]) -> OgeomResult<(Patch, Vec<Vec<u32>>)> {
        let mut local: FastMap<u32, u32> = FastMap::default();
        let mut vertex = Vec::new();
        for &h in edges {
            let v = self.tail(h);
            local.entry(self.point[v as usize]).or_insert_with(|| {
                vertex.push(v);
                index(vertex.len() - 1)
            });
        }
        let r = vertex.len();
        let mut points: Vec<u32> = vertex.iter().map(|&v| self.point[v as usize]).collect();
        points.sort_unstable();
        let faces = self.about(&points);

        // Edges are told apart by their points, so every face at a rim
        // point counts for the edges there. The Laplacian and the spacing
        // read only the fans the rim runs through: a vertex of the rim.
        let on_rim: FastSet<u32> = edges.iter().map(|&h| self.tail(h)).collect();
        let mut uses: FastMap<u64, u32> = FastMap::default();
        let mut joined = FastSet::default();
        let mut next_to: Vec<Vec<(u32, f64)>> = vec![Vec::new(); r];
        let mut mass = vec![0.0; r];
        for &f in &faces {
            let c = self.corners[f as usize];
            let pull = spokes(c.map(|v| self.positions[v as usize]));
            for k in 0..3 {
                *uses.entry(self.edge_key(3 * f + index(k))).or_insert(0) += 1;
                let ends = [c[k], c[(k + 1) % 3]].map(|v| local.get(&self.point[v as usize]));
                if let [Some(&l), Some(&m)] = ends {
                    joined.insert((l.min(m), l.max(m)));
                }
                if !on_rim.contains(&c[k]) {
                    continue;
                }
                if let Some(&l) = local.get(&self.point[c[k] as usize]) {
                    let (a, b) = (c[(k + 1) % 3], c[(k + 2) % 3]);
                    let (wa, wb) = pull.map_or((0.0, 0.0), |(w, third)| {
                        mass[l as usize] += third;
                        (w[(k + 2) % 3], w[(k + 1) % 3])
                    });
                    next_to[l as usize].push((self.point[a as usize], wa));
                    next_to[l as usize].push((self.point[b as usize], wb));
                }
            }
        }
        for &h in edges {
            let n = uses.get(&self.edge_key(h)).copied().unwrap_or(0);
            if n != 1 {
                ogeom_bail!(
                    Construction,
                    "rim edge {h} lies on an edge {n} faces use, which a patch cannot close"
                );
            }
        }

        let positions: Vec<Point> = vertex.iter().map(|&v| self.positions[v as usize]).collect();
        let mut near = Vec::with_capacity(r);
        let mut spacing = Vec::with_capacity(r);
        for (l, list) in next_to.iter_mut().enumerate() {
            merge(list);
            let mut around = Near {
                rim: Vec::new(),
                ring: Vec::new(),
                mass: mass[l],
            };
            let mut length = 0.0;
            for &(q, w) in list.iter() {
                let at = self.positions[q as usize];
                length += at.distance(positions[l]);
                match local.get(&q) {
                    Some(&m) => around.rim.push((m, w)),
                    None => around.ring.push((at, w)),
                }
            }
            spacing.push(length / f64::from(index(list.len().max(1))));
            near.push(around);
        }

        // Parted where a point comes round again: the edges since it last
        // left that point close a loop.
        let mut at = vec![usize::MAX; r];
        let mut stack: Vec<u32> = Vec::new();
        let mut loops = Vec::new();
        let slot = |h: u32| local[&self.point[self.tail(h) as usize]] as usize;
        for &h in edges {
            let l = slot(h);
            if at[l] != usize::MAX {
                let closed: Vec<u32> = stack.drain(at[l]..).collect();
                for &e in &closed {
                    at[slot(e)] = usize::MAX;
                }
                loops.push(closed);
            }
            at[l] = stack.len();
            stack.push(h);
        }
        loops.push(stack);
        if let Some(short) = loops.iter().find(|l| l.len() < 3) {
            ogeom_bail!(
                Construction,
                "rim {} has a loop of {} edges, which no face can fill",
                edges[0],
                short.len()
            );
        }

        Ok((
            Patch {
                vertex,
                local,
                positions,
                spacing,
                triangles: Vec::new(),
                joined,
                near,
            },
            loops,
        ))
    }

    /// Add the patch's vertices and faces, and join them to the mesh.
    /// Returns the faces and vertices made.
    fn commit(&mut self, patch: &Patch) -> (Range<u32>, Range<u32>) {
        let first_face = index(self.corners.len());
        let first_vertex = index(self.positions.len());
        let r = patch.vertex.len();
        let mut global = patch.vertex.clone();
        for &p in &patch.positions[r..] {
            global.push(self.add_vertex(p));
        }
        let mut points: Vec<u32> = patch
            .vertex
            .iter()
            .map(|&v| self.point[v as usize])
            .collect();
        points.sort_unstable();
        let mut faces = self.about(&points);
        points.extend(first_vertex..index(self.positions.len()));
        for t in &patch.triangles {
            faces.push(self.add_face(t.map(|l| global[l as usize])));
        }
        self.relink(&Scope::Points(points), &faces);
        (
            first_face..index(self.corners.len()),
            first_vertex..index(self.positions.len()),
        )
    }
}

/// A patch being built, numbered on its own: the rim's points first, one
/// number each, then the vertices the refinement makes.
struct Patch {
    /// For each rim point, a vertex of the mesh standing for it.
    vertex: Vec<u32>,
    /// The patch number of each rim point.
    local: FastMap<u32, u32>,
    positions: Vec<Point>,
    /// The spacing about each vertex: for a rim point, the mean length of
    /// the mesh's edges at it; for a new vertex, the mean of its
    /// triangle's corners'.
    spacing: Vec<f64>,
    /// The faces, wound as the mesh's.
    triangles: Vec<[u32; 3]>,
    /// The pairs of rim points the mesh joins, smaller first.
    joined: FastSet<(u32, u32)>,
    /// For each rim point, what the mesh's faces about it add to the
    /// Laplacian there.
    near: Vec<Near>,
}

/// What the mesh's faces about a rim point add to the cotangent Laplacian
/// there.
struct Near {
    /// The rim points it shares an edge with, by patch number, each with
    /// the edge's weight.
    rim: Vec<(u32, f64)>,
    /// The other points it shares an edge with, by position, each with the
    /// edge's weight.
    ring: Vec<(Point, f64)>,
    /// A third of the faces' area.
    mass: f64,
}

/// For a triangle, half the cotangent of the angle at each corner, the
/// weight of the edge facing it in the cotangent Laplacian, and a third of
/// the area; `None` for one with collinear corners, which has no angles.
fn spokes(c: [Point; 3]) -> Option<([f64; 3], f64)> {
    let twice = (c[1] - c[0]).cross(c[2] - c[0]).magnitude();
    if twice == 0.0 {
        return None;
    }
    let half_cot = |k: usize| {
        let (a, b) = (c[(k + 1) % 3] - c[k], c[(k + 2) % 3] - c[k]);
        a.dot(b) / twice / 2.0
    };
    Some(([half_cot(0), half_cot(1), half_cot(2)], twice / 6.0))
}

/// Sort `list` by its first members and sum the second over each run of
/// equal first.
fn merge(list: &mut Vec<(u32, f64)>) {
    list.sort_by_key(|e| e.0);
    let mut out: Vec<(u32, f64)> = Vec::with_capacity(list.len());
    for &(k, w) in list.iter() {
        match out.last_mut() {
            Some(last) if last.0 == k => last.1 += w,
            _ => out.push((k, w)),
        }
    }
    *list = out;
}

/// A loop to triangulate.
struct Loop {
    /// The patch numbers of its points, in the rim's order.
    ids: Vec<u32>,
    points: Vec<Point>,
    /// Beyond the edge from point `i` to `i + 1`, the far corner of the
    /// face on it, where that face has a normal.
    beyond: Vec<Option<Point>>,
}

impl Patch {
    /// The loop of the rim's half-edges `edges`.
    fn ring(&self, mesh: &Editable, edges: &[u32]) -> Loop {
        let at = |v: u32| mesh.positions[v as usize];
        let ids: Vec<u32> = edges
            .iter()
            .map(|&h| self.local[&mesh.point[mesh.tail(h) as usize]])
            .collect();
        let points = ids.iter().map(|&l| self.positions[l as usize]).collect();
        let beyond = edges
            .iter()
            .map(|&h| {
                let (a, b, o) = (at(mesh.tail(h)), at(mesh.head(h)), at(mesh.tail(prev(h))));
                ((b - a).cross(o - a) != Vector::ZERO).then_some(o)
            })
            .collect();
        Loop {
            ids,
            points,
            beyond,
        }
    }

    /// Refine to the mesh's spacing: centroid splits and flips, as the
    /// module documentation says.
    fn refine(&mut self) {
        let widest = self.spacing.iter().copied().fold(0.0, f64::max);
        if widest <= 0.0 {
            return;
        }
        // Cells as wide as the farthest a crowding vertex can stand, so the
        // 27 about a point hold every one.
        let width = widest * FRAC_1_SQRT_2;
        #[allow(clippy::cast_possible_truncation, reason = "saturating")]
        let cell = |p: Point| [p.x, p.y, p.z].map(|x| (x / width).floor() as i64);
        let mut grid: FastMap<[i64; 3], Vec<u32>> = FastMap::default();
        for (v, &p) in self.positions.iter().enumerate() {
            grid.entry(cell(p)).or_default().push(index(v));
        }
        loop {
            let mut split = false;
            for t in 0..self.triangles.len() {
                let [a, b, c] = self.triangles[t];
                let corners = [a, b, c].map(|v| self.positions[v as usize]);
                let sum = corners
                    .iter()
                    .fold(Vector::ZERO, |s, &p| s + (p - Point::ORIGIN));
                let g = Point::ORIGIN + sum / 3.0;
                let s = [a, b, c]
                    .iter()
                    .map(|&v| self.spacing[v as usize])
                    .sum::<f64>()
                    / 3.0;
                // Farther than σ/√2: twice the square past σ².
                let far = |d2: f64, sigma: f64| 2.0 * d2 > sigma * sigma;
                let corners_clear = [a, b, c].iter().zip(&corners).all(|(&v, &p)| {
                    let d2 = g.square_distance(p);
                    far(d2, self.spacing[v as usize])
                });
                if !corners_clear {
                    continue;
                }
                let [x, y, z] = cell(g);
                let mut crowded = false;
                'cells: for dx in -1..=1_i64 {
                    for dy in -1..=1_i64 {
                        for dz in -1..=1_i64 {
                            let key = [
                                x.saturating_add(dx),
                                y.saturating_add(dy),
                                z.saturating_add(dz),
                            ];
                            if let Some(list) = grid.get(&key)
                                && list.iter().any(|&v| {
                                    !far(g.square_distance(self.positions[v as usize]), s)
                                })
                            {
                                crowded = true;
                                break 'cells;
                            }
                        }
                    }
                }
                if crowded {
                    continue;
                }
                let v = index(self.positions.len());
                self.positions.push(g);
                self.spacing.push(s);
                grid.entry(cell(g)).or_default().push(v);
                self.triangles[t] = [a, b, v];
                self.triangles.push([b, c, v]);
                self.triangles.push([c, a, v]);
                split = true;
            }
            if !split {
                break;
            }
            self.relax();
        }
    }

    /// Flip each edge inside the patch while that raises the smallest
    /// angle of its two faces, and neither turns over against the other
    /// nor draws an edge there already is.
    fn relax(&mut self) {
        let r = index(self.vertex.len());
        let mut at: FastMap<(u32, u32), u32> = FastMap::default();
        for (t, c) in self.triangles.iter().enumerate() {
            for k in 0..3 {
                at.insert((c[k], c[(k + 1) % 3]), index(t));
            }
        }
        let positions = &self.positions;
        let p = |v: u32| positions[v as usize];
        loop {
            let mut flipped = false;
            'triangles: for t in 0..self.triangles.len() {
                for k in 0..3 {
                    let tri = self.triangles[t];
                    let (u, v, c) = (tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3]);
                    if u > v {
                        continue;
                    }
                    let Some(&t2) = at.get(&(v, u)) else { continue };
                    let other = self.triangles[t2 as usize];
                    let Some(&d) = other.iter().find(|&&w| w != u && w != v) else {
                        continue;
                    };
                    if c == d
                        || at.contains_key(&(c, d))
                        || at.contains_key(&(d, c))
                        || (c < r && d < r && self.joined.contains(&(c.min(d), c.max(d))))
                    {
                        continue;
                    }
                    let n1 = (p(d) - p(u)).cross(p(c) - p(u));
                    let n2 = (p(v) - p(d)).cross(p(c) - p(d));
                    if n1 == Vector::ZERO || n2 == Vector::ZERO || n1.dot(n2) <= 0.0 {
                        continue;
                    }
                    let before =
                        smallest_angle(p(u), p(v), p(c)).min(smallest_angle(p(v), p(u), p(d)));
                    let after =
                        smallest_angle(p(u), p(d), p(c)).min(smallest_angle(p(d), p(v), p(c)));
                    if after <= before {
                        continue;
                    }
                    at.remove(&(u, v));
                    at.remove(&(v, u));
                    self.triangles[t] = [u, d, c];
                    self.triangles[t2 as usize] = [d, v, c];
                    for (e, f) in [((u, d), t), ((d, c), t), ((c, u), t)] {
                        at.insert(e, index(f));
                    }
                    for e in [(d, v), (v, c), (c, d)] {
                        at.insert(e, t2);
                    }
                    flipped = true;
                    continue 'triangles;
                }
            }
            if !flipped {
                break;
            }
        }
    }

    /// Move the new vertices as `fairing` asks.
    fn fair(&mut self, fairing: Fairing) -> OgeomResult<()> {
        let r = self.vertex.len();
        let total = self.positions.len();
        let unknowns = total - r;
        if unknowns == 0 {
            return Ok(());
        }
        // The cotangent Laplacian of the patch as it stands, with the
        // mesh's faces about the rim.
        let mut spoke: Vec<Vec<(u32, f64)>> = vec![Vec::new(); total];
        let mut mass = vec![0.0; total];
        for t in &self.triangles {
            let Some((w, third)) = spokes(t.map(|v| self.positions[v as usize])) else {
                continue;
            };
            for k in 0..3 {
                mass[t[k] as usize] += third;
                let (a, b) = (t[(k + 1) % 3], t[(k + 2) % 3]);
                spoke[a as usize].push((b, w[k]));
                spoke[b as usize].push((a, w[k]));
            }
        }
        for (l, near) in self.near.iter().enumerate() {
            spoke[l].extend(&near.rim);
            mass[l] += near.mass;
        }
        for list in &mut spoke {
            merge(list);
        }
        let column = |l: u32| (l as usize).checked_sub(r);
        let rows: Vec<usize> = (0..total)
            .filter(|&w| {
                w >= r
                    || (fairing == Fairing::Curved
                        && mass[w] > 0.0
                        && spoke[w].iter().any(|&(m, _)| m as usize >= r))
            })
            .collect();
        // The three coordinates as one system, the x rows and columns
        // first, then y, then z. A Curved row is the Laplacian times the
        // root of its vertex's area, so that the squares sum to the
        // thin-plate energy.
        let mut triplets = Vec::new();
        let mut rhs = vec![0.0; 3 * rows.len()];
        let none: Vec<(Point, f64)> = Vec::new();
        for (row, &w) in rows.iter().enumerate() {
            let ring = if w < r { &self.near[w].ring } else { &none };
            // A vertex whose faces all have collinear corners has no
            // weights, and its row is empty however scaled.
            let scale = match fairing {
                Fairing::Curved if mass[w] > 0.0 => 1.0 / mass[w].sqrt(),
                _ => 1.0,
            };
            let here = self.positions[w];
            let mut lap = Vector::ZERO;
            let mut sum = 0.0;
            let mut put = |c: usize, value: f64| {
                for k in 0..3 {
                    triplets.push((k * rows.len() + row, k * unknowns + c, value));
                }
            };
            for &(m, weight) in &spoke[w] {
                lap += (self.positions[m as usize] - here) * weight;
                sum += weight;
                if let Some(c) = column(m) {
                    put(c, scale * weight);
                }
            }
            for &(q, weight) in ring {
                lap += (q - here) * weight;
                sum += weight;
            }
            if let Some(c) = column(index(w)) {
                put(c, -scale * sum);
            }
            let lap = lap * scale;
            rhs[row] = -lap.x;
            rhs[rows.len() + row] = -lap.y;
            rhs[2 * rows.len() + row] = -lap.z;
        }
        let a = SparseMatrix::from_triplets(3 * rows.len(), 3 * unknowns, &triplets);
        // Conjugate gradients end in as many steps as the normal matrix
        // has distinct eigenvalues, at most `unknowns` for its three equal
        // blocks; three times that is allowed for rounding.
        let Some(d) = least_squares_cgnr(&a, &rhs, SETTLED, 3 * unknowns) else {
            ogeom_bail!(
                Numeric,
                "the fairing of {unknowns} new vertices did not settle"
            );
        };
        for c in 0..unknowns {
            self.positions[r + c] += Vector::new(d[c], d[unknowns + c], d[2 * unknowns + c]);
        }
        Ok(())
    }
}

/// Where the fairing solve stops: the residual's gradient this small a
/// part of where it began. On the patches measured, of up to some
/// thousands of vertices, the displacement is then within a few units of
/// rounding of the exact one, for about a tenth more steps than a stop a
/// million times looser.
const SETTLED: f64 = 1e-12;

/// The smallest angle of the triangle `a b c`.
fn smallest_angle(a: Point, b: Point, c: Point) -> f64 {
    let angle = |u: Vector, v: Vector| u.cross(v).magnitude().atan2(u.dot(v));
    angle(b - a, c - a)
        .min(angle(c - b, a - b))
        .min(angle(a - c, b - c))
}

/// A frame about the line from one point to another, in which [`turn`]
/// reads the way a point stands from it.
///
/// [`turn`]: Axis::turn
struct Axis {
    from: Point,
    e1: Vector,
    e2: Vector,
}

impl Axis {
    fn new(a: Point, b: Point) -> Self {
        let u = b - a;
        let length = u.magnitude();
        // The coordinate axis least along the line makes a frame with it.
        let pick = if u.x.abs() <= u.y.abs() && u.x.abs() <= u.z.abs() {
            Vector::X
        } else if u.y.abs() <= u.z.abs() {
            Vector::Y
        } else {
            Vector::Z
        };
        let e1 = u.cross(pick);
        let e1 = e1 / e1.magnitude();
        let e2 = u.cross(e1) / length;
        Self { from: a, e1, e2 }
    }

    /// The unit direction, square to the line, from it to `p`.
    fn turn(&self, p: Point) -> [f64; 2] {
        let w = p - self.from;
        let (x, y) = (w.dot(self.e1), w.dot(self.e2));
        // Scaled to the larger first, so the squares neither overflow nor
        // underflow.
        let m = x.abs().max(y.abs());
        let (x, y) = (x / m, y / m);
        let r = (x * x + y * y).sqrt();
        [x / r, y / r]
    }
}

/// How far two faces sharing an edge are bent, from the turns of their far
/// corners about it: `4 sin²(β / 2)` for the bend `β`, π less the angle
/// between the two half-planes, so 0 for faces flat across the edge and 4
/// for faces folded onto each other. It grows with the bend, and is read
/// from the sum of the two directions, which keeps its precision for small
/// bends.
fn bend(a: [f64; 2], b: [f64; 2]) -> f64 {
    let (x, y) = (a[0] + b[0], a[1] + b[1]);
    x * x + y * y
}

/// The bend `β`, in radians, of the measure [`bend`] gives.
fn radians(fold: f64) -> f64 {
    2.0 * (fold.sqrt() / 2.0).min(1.0).asin()
}

/// The frame of the edge between loop points `i` and `j`, from the lower.
fn axis_of(ring: &Loop, i: usize, j: usize) -> Axis {
    Axis::new(ring.points[i.min(j)], ring.points[i.max(j)])
}

/// The largest bend, as [`bend`] measures it, and the area of a
/// triangulation of `ring`, given as triples of loop indices.
fn measure(ring: &Loop, triangles: &[[usize; 3]]) -> (f64, f64) {
    let n = ring.points.len();
    let p = &ring.points;
    let mut apexes: FastMap<(usize, usize), Vec<usize>> = FastMap::default();
    let mut area = 0.0;
    for t in triangles {
        area += (p[t[1]] - p[t[0]]).cross(p[t[2]] - p[t[0]]).magnitude() / 2.0;
        for k in 0..3 {
            let (i, j) = (t[k], t[(k + 1) % 3]);
            apexes
                .entry((i.min(j), i.max(j)))
                .or_default()
                .push(t[(k + 2) % 3]);
        }
    }
    let mut keys: Vec<(usize, usize)> = apexes.keys().copied().collect();
    keys.sort_unstable();
    let mut largest = 0.0_f64;
    for (i, j) in keys {
        let axis = axis_of(ring, i, j);
        let sides = &apexes[&(i, j)];
        let b = if j - i == 1 || (i == 0 && j == n - 1) {
            let edge = if j - i == 1 { i } else { n - 1 };
            ring.beyond[edge].map_or(0.0, |o| bend(axis.turn(o), axis.turn(p[sides[0]])))
        } else {
            bend(axis.turn(p[sides[0]]), axis.turn(p[sides[1]]))
        };
        largest = largest.max(b);
    }
    (largest, area)
}

/// The triangulation of `ring` with the least largest bend, then the least
/// area, drawing a new edge between loop points `i < j` only where
/// `free(i, j)`. Triangles are wound as the mesh's faces, against the loop;
/// `None` where every triangulation has a face of collinear corners or an
/// edge not free.
fn optimal(ring: &Loop, free: &dyn Fn(usize, usize) -> bool) -> Option<Vec<[usize; 3]>> {
    const NO: f64 = f64::INFINITY;
    let n = ring.points.len();
    let p = &ring.points;
    let mut open = vec![false; n * n];
    for i in 0..n {
        for j in i + 1..n {
            open[i * n + j] = j - i == 1 || (i == 0 && j == n - 1) || free(i, j);
        }
    }
    // The triangles on edge (i, k), one per apex between, from `start`.
    let mut start = vec![0_usize; n * n];
    let mut total = 0;
    for i in 0..n {
        for k in i + 2..n {
            start[i * n + k] = total;
            total += k - i - 1;
        }
    }
    let slot = |i: usize, m: usize, k: usize| start[i * n + k] + (m - i - 1);
    let valid = |i: usize, m: usize, k: usize| {
        open[i * n + m]
            && open[m * n + k]
            && open[i * n + k]
            && (p[m] - p[i]).cross(p[k] - p[i]) != Vector::ZERO
    };
    // The triangle a value about edge (a, b) and a point `c` outside goes
    // to.
    let consumer = |a: usize, b: usize, c: usize| {
        if c < a { slot(c, a, b) } else { slot(a, b, c) }
    };
    let outside = |a: usize, b: usize| (0..a).chain(b + 1..n);

    // First pass: the largest bend within each triangle's side of its
    // edge, least over the triangulations there.
    let mut value = vec![0.0_f64; total];
    let each = |value: &mut Vec<f64>, init: &dyn Fn(usize, usize, usize) -> f64| {
        for i in 0..n {
            for k in i + 2..n {
                for m in i + 1..k {
                    value[slot(i, m, k)] = if valid(i, m, k) { init(i, m, k) } else { NO };
                }
            }
        }
    };
    each(&mut value, &|_, _, _| 0.0);
    for a in 0..n - 1 {
        let axis = axis_of(ring, a, a + 1);
        let far = ring.beyond[a].map(|o| axis.turn(o));
        for c in outside(a, a + 1) {
            let b = far.map_or(0.0, |t| bend(t, axis.turn(p[c])));
            let s = consumer(a, a + 1, c);
            value[s] = value[s].max(b);
        }
    }
    let mut candidates: Vec<(f64, usize, [f64; 2])> = Vec::new();
    let mut theta = NO;
    for len in 2..n {
        for a in 0..n - len {
            let b = a + len;
            let axis = axis_of(ring, a, b);
            gather_candidates(&mut candidates, &value, &axis, p, a, b, &slot);
            let least = |t: Option<[f64; 2]>| {
                let mut best = NO;
                for &(w, _, tx) in &candidates {
                    if w >= best {
                        break;
                    }
                    best = best.min(w.max(t.map_or(0.0, |t| bend(tx, t))));
                }
                best
            };
            if len == n - 1 {
                theta = least(ring.beyond[n - 1].map(|o| axis.turn(o)));
            } else {
                for c in outside(a, b) {
                    let f = least(Some(axis.turn(p[c])));
                    let s = consumer(a, b, c);
                    value[s] = value[s].max(f);
                }
            }
        }
    }
    if theta == NO {
        return None;
    }

    // Second pass: the least area within each triangle's side of its
    // edge, over the triangulations there bending nowhere more than θ.
    each(&mut value, &|i, m, k| {
        (p[m] - p[i]).cross(p[k] - p[i]).magnitude() / 2.0
    });
    let fits = |far: Option<[f64; 2]>, t: [f64; 2]| far.is_none_or(|o| bend(o, t) <= theta);
    for a in 0..n - 1 {
        let axis = axis_of(ring, a, a + 1);
        let far = ring.beyond[a].map(|o| axis.turn(o));
        for c in outside(a, a + 1) {
            if !fits(far, axis.turn(p[c])) {
                value[consumer(a, a + 1, c)] = NO;
            }
        }
    }
    let mut root = None;
    for len in 2..n {
        for a in 0..n - len {
            let b = a + len;
            let axis = axis_of(ring, a, b);
            gather_candidates(&mut candidates, &value, &axis, p, a, b, &slot);
            let first = |t: Option<[f64; 2]>| {
                candidates
                    .iter()
                    .find(|&&(_, _, tx)| fits(t, tx))
                    .map(|&(g, x, _)| (g, x))
            };
            if len == n - 1 {
                root = first(ring.beyond[n - 1].map(|o| axis.turn(o)));
            } else {
                for c in outside(a, b) {
                    let s = consumer(a, b, c);
                    value[s] = match first(Some(axis.turn(p[c]))) {
                        Some((g, _)) => value[s] + g,
                        None => NO,
                    };
                }
            }
        }
    }
    let (_, m) = root?;

    // Each side's best apex given the triangle across its edge, as the
    // second pass chose it: least area, then lowest apex.
    let choose = |a: usize, b: usize, c: usize| {
        let axis = axis_of(ring, a, b);
        let t = axis.turn(p[c]);
        let mut best: Option<(f64, usize)> = None;
        for x in a + 1..b {
            let g = value[slot(a, x, b)];
            if g < NO && bend(axis.turn(p[x]), t) <= theta && best.is_none_or(|(h, _)| g < h) {
                best = Some((g, x));
            }
        }
        best.map(|(_, x)| x)
    };
    let mut triangles = Vec::with_capacity(n - 2);
    let mut stack = vec![(0, m, n - 1)];
    while let Some((i, m, k)) = stack.pop() {
        triangles.push([k, m, i]);
        if m - i >= 2 {
            stack.push((i, choose(i, m, k)?, m));
        }
        if k - m >= 2 {
            stack.push((m, choose(m, k, i)?, k));
        }
    }
    Some(triangles)
}

/// The apexes on edge `(a, b)` with a finite value: the value, the apex,
/// and its turn about the edge, sorted by value and then apex.
fn gather_candidates(
    candidates: &mut Vec<(f64, usize, [f64; 2])>,
    value: &[f64],
    axis: &Axis,
    p: &[Point],
    a: usize,
    b: usize,
    slot: &dyn Fn(usize, usize, usize) -> usize,
) {
    candidates.clear();
    for x in a + 1..b {
        let w = value[slot(a, x, b)];
        if w < f64::INFINITY {
            candidates.push((w, x, axis.turn(p[x])));
        }
    }
    candidates.sort_unstable_by(|u, v| u.0.total_cmp(&v.0).then(u.1.cmp(&v.1)));
}

/// Ear-clip `ring` on the plane across which it encloses the most area,
/// drawing a new edge between loop points `i < j` only where `free(i, j)`.
/// Triangles are wound as the mesh's faces, against the loop.
fn ear_clip(ring: &Loop, free: &dyn Fn(usize, usize) -> bool) -> OgeomResult<Vec<[usize; 3]>> {
    let n = ring.points.len();
    let p = &ring.points;
    // Against the loop, the patch's own way round.
    let order: Vec<usize> = (0..n).map(|j| (n - j) % n).collect();
    let centre = p.iter().fold(Vector::ZERO, |s, &q| s + (q - Point::ORIGIN)) / f64::from(index(n));
    let at = |j: usize| p[order[j]] - Point::ORIGIN - centre;
    let mut normal = Vector::ZERO;
    for j in 0..n {
        normal += at(j).cross(at((j + 1) % n));
    }
    if normal == Vector::ZERO {
        ogeom_bail!(
            Construction,
            "a loop of {n} edges encloses no area across any plane"
        );
    }
    let axis = Axis::new(Point::ORIGIN, Point::ORIGIN + normal);
    // `e1`, `e2` and the normal are right-handed: the loop runs
    // counter-clockwise in those coordinates.
    let q: Vec<[f64; 2]> = (0..n)
        .map(|j| {
            let w = at(j);
            [w.dot(axis.e1), w.dot(axis.e2)]
        })
        .collect();
    let orient = |a: usize, b: usize, c: usize| Exact::orient2d(q[a], q[b], q[c]);
    if !simple(&q) {
        ogeom_bail!(
            Construction,
            "the shadow of a loop of {n} edges on its plane crosses itself"
        );
    }

    let mut next: Vec<usize> = (0..n).map(|j| (j + 1) % n).collect();
    let mut before: Vec<usize> = (0..n).map(|j| (j + n - 1) % n).collect();
    let mut gone = vec![false; n];
    let convex = |before: &[usize], next: &[usize], j: usize| {
        orient(before[j], j, next[j]) == Sign::Positive
    };
    let mut triangles = Vec::with_capacity(n - 2);
    let mut left = n;
    let mut j = 0;
    let mut misses = 0;
    while left > 3 {
        let (a, c) = (before[j], next[j]);
        let ear = convex(&before, &next, j)
            && free(order[a].min(order[c]), order[a].max(order[c]))
            && (0..n).all(|v| {
                gone[v]
                    || v == a
                    || v == j
                    || v == c
                    || convex(&before, &next, v)
                    || orient(a, j, v) == Sign::Negative
                    || orient(j, c, v) == Sign::Negative
                    || orient(c, a, v) == Sign::Negative
            });
        if ear {
            triangles.push([order[a], order[j], order[c]]);
            gone[j] = true;
            next[a] = c;
            before[c] = a;
            left -= 1;
            misses = 0;
            j = a;
        } else {
            j = next[j];
            misses += 1;
            if misses > left {
                ogeom_bail!(
                    NotDone,
                    "ear clipping a loop of {n} edges found no ear to cut with {left} points left"
                );
            }
        }
    }
    let (a, c) = (before[j], next[j]);
    if orient(a, j, c) != Sign::Positive {
        ogeom_bail!(
            NotDone,
            "ear clipping a loop of {n} edges left a last triangle with no area"
        );
    }
    triangles.push([order[a], order[j], order[c]]);
    Ok(triangles)
}

/// Whether the closed polygon `q` is simple: no two of its edges meet but
/// neighbours at their shared corner, and no neighbours run back along
/// each other. Edges are swept in order of their least `x`.
fn simple(q: &[[f64; 2]]) -> bool {
    let n = q.len();
    let edge = |j: usize| (q[j], q[(j + 1) % n]);
    let mut order: Vec<usize> = (0..n).collect();
    let low = |j: usize| edge(j).0[0].min(edge(j).1[0]);
    let high = |j: usize| edge(j).0[0].max(edge(j).1[0]);
    order.sort_by(|&a, &b| low(a).total_cmp(&low(b)).then(a.cmp(&b)));
    let mut active: Vec<usize> = Vec::new();
    for &j in &order {
        active.retain(|&k| high(k) >= low(j));
        for &k in &active {
            let (a, b) = edge(j);
            let (c, d) = edge(k);
            let meets = if (j + 1) % n == k {
                folds_back(a, b, d)
            } else if (k + 1) % n == j {
                folds_back(c, d, b)
            } else {
                segments_meet(a, b, c, d)
            };
            if meets {
                return false;
            }
        }
        active.push(j);
    }
    true
}

/// Whether the edge from `b` to `c` runs back along the edge from `a` to
/// `b`.
fn folds_back(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> bool {
    Exact::orient2d(a, b, c) == Sign::Zero
        && (a[0] - b[0]) * (c[0] - b[0]) + (a[1] - b[1]) * (c[1] - b[1]) > 0.0
}

/// Whether the closed segments `a b` and `c d` share a point, exactly.
fn segments_meet(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    let within = |p: [f64; 2], q: [f64; 2], r: [f64; 2]| {
        r[0] >= p[0].min(q[0])
            && r[0] <= p[0].max(q[0])
            && r[1] >= p[1].min(q[1])
            && r[1] <= p[1].max(q[1])
    };
    let (o1, o2) = (Exact::orient2d(a, b, c), Exact::orient2d(a, b, d));
    let (o3, o4) = (Exact::orient2d(c, d, a), Exact::orient2d(c, d, b));
    if (o1 == Sign::Zero && within(a, b, c))
        || (o2 == Sign::Zero && within(a, b, d))
        || (o3 == Sign::Zero && within(c, d, a))
        || (o4 == Sign::Zero && within(c, d, b))
    {
        return true;
    }
    let opposite = |s: Sign, t: Sign| s != Sign::Zero && t != Sign::Zero && s != t;
    opposite(o1, o2) && opposite(o3, o4)
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    reason = "test code"
)]
mod tests {
    use super::super::tests::{check, grid, mesh, torus};
    use super::*;
    use core::f64::consts::{PI, TAU};
    use ogeom_topo::Triangulation;
    use proptest::prelude::*;

    const CURVED: Fill = Fill {
        fairing: Fairing::Curved,
        refine: true,
    };
    const FLAT: Fill = Fill {
        fairing: Fairing::Flat,
        refine: true,
    };
    const BARE: Fill = Fill {
        fairing: Fairing::Flat,
        refine: false,
    };

    /// A sphere of latitude bands and longitude sides about poles on `z`,
    /// wound outward.
    fn sphere(radius: f64, bands: u32, sides: u32) -> Triangulation {
        let mut points = vec![Point::new(0.0, 0.0, radius)];
        for j in 1..bands {
            let polar = PI * f64::from(j) / f64::from(bands);
            let (s, c) = (radius * polar.sin(), radius * polar.cos());
            for i in 0..sides {
                let a = TAU * f64::from(i) / f64::from(sides);
                points.push(Point::new(s * a.cos(), s * a.sin(), c));
            }
        }
        let south = 1 + (bands - 1) * sides;
        points.push(Point::new(0.0, 0.0, -radius));
        let at = |j: u32, i: u32| 1 + (j - 1) * sides + i % sides;
        let mut triangles = Vec::new();
        for i in 0..sides {
            triangles.push([0, at(1, i), at(1, i + 1)]);
            triangles.push([south, at(bands - 1, i + 1), at(bands - 1, i)]);
        }
        for j in 1..bands - 1 {
            for i in 0..sides {
                let (a, b, c, d) = (at(j, i), at(j + 1, i), at(j + 1, i + 1), at(j, i + 1));
                triangles.extend([[a, b, c], [a, c, d]]);
            }
        }
        mesh(points, triangles)
    }

    /// The faces whose centroids pass `keep`.
    fn faces_where(e: &Editable, keep: impl Fn(Point) -> bool) -> Vec<u32> {
        (0..index(e.face_count()))
            .filter(|&f| {
                e.face(f).is_some_and(|c| {
                    let sum = c.iter().fold(Vector::ZERO, |s, &v| {
                        s + (e.position(v).unwrap() - Point::ORIGIN)
                    });
                    keep(Point::ORIGIN + sum / 3.0)
                })
            })
            .collect()
    }

    /// Fill every rim, each as `fill` asks.
    fn fill_all(e: &mut Editable, fill: &Fill) -> Vec<FillReport> {
        let mut reports = Vec::new();
        while let Some(rim) = e.rims().first() {
            reports.push(e.fill(rim.id(), fill).unwrap());
            check(e);
        }
        reports
    }

    /// The area of the faces `faces`.
    fn area_of(e: &Editable, faces: Range<u32>) -> f64 {
        faces
            .map(|f| {
                let [a, b, c] = e.face(f).unwrap().map(|v| e.position(v).unwrap());
                (b - a).cross(c - a).magnitude() / 2.0
            })
            .sum()
    }

    /// Every triangulation of the polygon `i..=k` on its edge `(i, k)`,
    /// wound as [`optimal`]'s.
    fn every(i: usize, k: usize) -> Vec<Vec<[usize; 3]>> {
        if k - i < 2 {
            return vec![Vec::new()];
        }
        let mut all = Vec::new();
        for m in i + 1..k {
            for left in every(i, m) {
                for right in every(m, k) {
                    let mut t = left.clone();
                    t.extend(&right);
                    t.push([k, m, i]);
                    all.push(t);
                }
            }
        }
        all
    }

    /// A small generator for the loops the programme is checked on.
    struct Draw(u64);

    impl Draw {
        fn next(&mut self) -> f64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            (self.0 >> 11) as f64 / (1_u64 << 53) as f64
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(256))]

        /// The programme's triangulation bends no more than any of a loop's
        /// triangulations, to the bit, and has the least area of those
        /// that bend as little, to the rounding of the sums: checked
        /// against every triangulation of loops of up to nine edges, with
        /// diagonals barred and rim faces with no normal among them.
        #[test]
        fn the_programme_is_optimal_over_every_triangulation(
            n in 3_usize..=9,
            seed in 1_u64..u64::MAX,
        ) {
            let mut draw = Draw(seed);
            let mut points = Vec::new();
            let mut beyond = Vec::new();
            for i in 0..n {
                let a = TAU * i as f64 / n as f64;
                let r = 1.0 + 0.6 * (draw.next() - 0.5);
                points.push(Point::new(r * a.cos(), r * a.sin(), draw.next() - 0.5));
            }
            for i in 0..n {
                let (p, q) = (points[i], points[(i + 1) % n]);
                let out = Vector::new(p.x + q.x, p.y + q.y, 0.0) * 0.5;
                let o = p + (q - p) * 0.5 + out * draw.next()
                    + Vector::new(0.0, 0.0, 2.0 * draw.next() - 1.0);
                beyond.push((draw.next() > 0.1).then_some(o));
            }
            let mut barred = vec![false; n * n];
            for b in &mut barred {
                *b = draw.next() < 0.15;
            }
            let ring = Loop { ids: (0..index(n)).collect(), points, beyond };
            let free = |i: usize, j: usize| !barred[i * n + j];
            let open = |i: usize, j: usize| {
                let (i, j) = (i.min(j), i.max(j));
                j - i == 1 || (i == 0 && j == n - 1) || free(i, j)
            };
            let p = &ring.points;
            let fits = |t: &[[usize; 3]]| {
                t.iter().all(|&[k, m, i]| {
                    open(i, m) && open(m, k) && open(i, k)
                        && (p[m] - p[i]).cross(p[k] - p[i]) != Vector::ZERO
                })
            };
            let best = every(0, n - 1)
                .into_iter()
                .filter(|t| fits(t))
                .map(|t| measure(&ring, &t))
                .min_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
            let found = optimal(&ring, &free);
            match (best, found) {
                (None, None) => {}
                (Some((bend, area)), Some(t)) => {
                    prop_assert_eq!(t.len(), n - 2);
                    prop_assert!(fits(&t));
                    let (b, a) = measure(&ring, &t);
                    prop_assert_eq!(b, bend);
                    prop_assert!((a - area).abs() <= 1e-12 * area, "{} against {}", a, area);
                }
                (best, found) => prop_assert!(false, "{:?} against {:?}", best, found),
            }
        }
    }

    #[test]
    fn a_sphere_cap_filled_curved_comes_back_to_the_sphere() {
        // A cap of polar angle α = π/8 taken off unit spheres of three
        // resolutions. The Laplacians are those of the patch as
        // triangulated, flat across the rim's plane, and of the sphere's
        // faces about it. In the limit the patch is biharmonic over the
        // disc of radius a = sin α, meeting the sphere's position and its
        // slope along the sphere at the rim: height A + B s², with
        // 2 B a = -sin α, so its top stands at cos α + sin² α / 2.
        let alpha = PI / 8.0;
        let top = alpha.cos() + alpha.sin().powi(2) / 2.0;
        let mut worst = Vec::new();
        let mut highest = Vec::new();
        for (bands, sides) in [(16, 24), (32, 48), (64, 96)] {
            let m = sphere(1.0, bands, sides);
            let mut e = Editable::from(&m, 0.0).unwrap();
            let rim_z = alpha.cos();
            e.remove_faces(&faces_where(&e, |c| c.z > rim_z));
            let rims = e.rims();
            assert_eq!(rims.len(), 1);
            assert_eq!(rims[0].edges.len(), sides as usize);
            let report = e.fill(rims[0].id(), &CURVED).unwrap();
            check(&e);
            assert!(e.is_closed());
            assert!(e.non_manifold().is_empty());
            assert_eq!(report.loops, 1);
            assert_eq!(report.flattened, 0);
            assert!(!report.vertices.is_empty());
            for (v, &p) in m.positions.iter().enumerate() {
                if let Some(q) = e.position(index(v)) {
                    assert_eq!(q, p);
                }
            }
            let new = || report.vertices.clone().map(|v| e.position(v).unwrap());
            worst.push(
                new()
                    .map(|p| ((p - Point::ORIGIN).magnitude() - 1.0).abs())
                    .fold(0.0, f64::max),
            );
            highest.push(new().map(|p| p.z).fold(f64::MIN, f64::max));
        }
        // Measured: 0.28 %, 0.38 % and 0.33 % of the radius, and on to
        // 0.31 % and 0.30 % at 128 and 256 bands, tending to the top's
        // 0.29 %.
        assert!(worst.iter().all(|&w| w < 0.01), "{worst:?}");
        let gaps: Vec<f64> = highest.iter().map(|h| (h - top).abs()).collect();
        assert!(
            gaps.windows(2).all(|g| g[1] < g[0]),
            "{highest:?} against {top}"
        );
    }

    #[test]
    fn a_flat_hole_filled_flat_stays_in_its_plane() {
        // A block of three by three cells with an arm of two, out of a grid
        // tilted out of every axis plane.
        let n = 8;
        let u = Vector::new(0.6, 0.3, -0.2);
        let v = Vector::new(-0.1, 0.4, 0.7);
        let normal = u.cross(v) / u.cross(v).magnitude();
        let origin = Point::new(3.0, -2.0, 5.0);
        let mut m = grid(n);
        for p in &mut m.positions {
            *p = origin + u * p.x + v * p.y;
        }
        let cell_area = u.cross(v).magnitude();
        let mut cells = vec![(4, 1), (5, 1)];
        for i in 1..4 {
            for j in 1..4 {
                cells.push((i, j));
            }
        }
        let removed: Vec<u32> = cells
            .iter()
            .flat_map(|&(i, j)| [2 * (n * j + i), 2 * (n * j + i) + 1])
            .collect();
        for fill in [FLAT, BARE] {
            let mut e = Editable::from(&m, 0.0).unwrap();
            e.remove_faces(&removed);
            let rim = e.rims().into_iter().find(|r| r.edges.len() == 16).unwrap();
            let report = e.fill(rim.id(), &fill).unwrap();
            check(&e);
            assert_eq!(e.rims().len(), 1);
            assert!(e.non_manifold().is_empty());
            assert!(
                report.largest_bend < 1e-14,
                "{fill:?} {}",
                report.largest_bend
            );
            assert_eq!(report.vertices.is_empty(), !fill.refine);
            for v in report.vertices.clone() {
                let h = (e.position(v).unwrap() - origin).dot(normal);
                assert!(h.abs() < 1e-14, "{fill:?} {h}");
            }
            let area = area_of(&e, report.faces.clone());
            let hole = cells.len() as f64 * cell_area;
            assert!((area - hole).abs() < 1e-13 * hole, "{fill:?} {area} {hole}");
        }
    }

    #[test]
    fn filling_every_hole_closes_a_torus() {
        let m = torus(16, 12);
        let whole = Editable::from(&m, 0.0).unwrap();
        // A lone face, a quad, a block of three by two quads, and two quads
        // touching at a corner.
        let mut removed = vec![5, 40, 41];
        for (i, j) in [
            (6, 3),
            (6, 4),
            (7, 3),
            (7, 4),
            (8, 3),
            (8, 4),
            (12, 8),
            (13, 9),
        ] {
            removed.extend([2 * (i * 12 + j), 2 * (i * 12 + j) + 1]);
        }
        for fill in [CURVED, FLAT, BARE] {
            let mut e = whole.clone();
            e.remove_faces(&removed);
            assert_eq!(e.rims().len(), 4);
            let reports = fill_all(&mut e, &fill);
            assert_eq!(reports.len(), 4);
            assert!(e.is_closed(), "{fill:?}");
            assert!(e.non_manifold().is_empty());
            assert_eq!(e.components().len(), 1);
            // A torus: V - E + F = 0, each face's three edges each shared.
            let t = e.to_triangulation();
            let (v, f) = (t.positions.len(), t.triangles.len());
            assert_eq!(2 * v, f, "{fill:?}");
            // The block's two inner vertices went with it.
            if !fill.refine {
                assert_eq!(v, m.positions.len() - 2);
            }
        }
    }

    #[test]
    fn holes_touching_at_a_point_are_closed_through_it() {
        // Two cells of a grid meeting at a corner leave one rim through
        // that point twice, at two copies of it.
        let m = grid(4);
        let mut e = Editable::from(&m, 0.0).unwrap();
        let cell = |i: u32, j: u32| [2 * (4 * j + i), 2 * (4 * j + i) + 1];
        e.remove_faces(&[cell(1, 1), cell(2, 2)].concat());
        let pinch = Point::new(2.0, 2.0, 0.0);
        assert_eq!(e.non_manifold().len(), 2);
        let rims = e.rims();
        assert_eq!(rims.len(), 2);
        let inner = rims.iter().find(|r| r.edges.len() == 8).unwrap();
        let at_pinch = inner
            .edges
            .iter()
            .filter(|&&h| e.position(e.edge(h).unwrap()[0]) == Some(pinch))
            .count();
        assert_eq!(at_pinch, 2);
        for fill in [BARE, CURVED] {
            let mut f = e.clone();
            let report = f.fill(inner.id(), &fill).unwrap();
            check(&f);
            assert_eq!(report.loops, 2);
            assert!(f.non_manifold().is_empty(), "{fill:?}");
            assert_eq!(f.rims().len(), 1);
            let back = f.to_triangulation();
            assert_eq!(
                back.positions.iter().filter(|&&p| p == pinch).count(),
                1,
                "{fill:?}"
            );
            if !fill.refine {
                assert_eq!(back.positions, m.positions);
                assert_eq!(back.triangles.len(), m.triangles.len());
            }
        }
    }

    #[test]
    fn an_edge_the_mesh_has_is_not_drawn_again() {
        // A square hole in a grid whose diagonal from (1, 1) to (2, 2) is
        // already an edge, of a fin standing on it. The fill draws the
        // other diagonal, and the fin's edge keeps its one face.
        let mut m = grid(3);
        let fin = index(m.positions.len());
        m.positions.push(Point::new(1.5, 1.5, 1.0));
        m.triangles.push([5, 10, fin]);
        let mut e = Editable::from(&m, 0.0).unwrap();
        e.remove_faces(&[8, 9]);
        let hole = e.rims().into_iter().find(|r| r.edges.len() == 4).unwrap();
        let report = e.fill(hole.id(), &BARE).unwrap();
        check(&e);
        let at = |v: u32| e.position(v).unwrap();
        for f in report.faces {
            let corners = e.face(f).unwrap().map(at);
            assert!(
                !(corners.contains(&Point::new(1.0, 1.0, 0.0))
                    && corners.contains(&Point::new(2.0, 2.0, 0.0))),
                "{corners:?}"
            );
        }
        let rims = e.rims();
        assert_eq!(rims.len(), 2);
        assert!(rims.iter().any(|r| r.edges.len() == 3));
    }

    #[test]
    fn a_long_loop_is_ear_clipped_on_its_plane() {
        // A flat ring of 600 sides between radii 1 and 2: both rims are
        // past the programme's bound.
        let sides = 600;
        let mut positions = Vec::new();
        for r in [1.0, 2.0] {
            for i in 0..sides {
                let a = TAU * f64::from(i) / f64::from(sides);
                positions.push(Point::new(r * a.cos(), r * a.sin(), 0.0));
            }
        }
        let mut triangles = Vec::new();
        for i in 0..sides {
            let j = (i + 1) % sides;
            triangles.extend([[i, sides + i, sides + j], [i, sides + j, j]]);
        }
        let m = mesh(positions, triangles);
        let polygon = |r: f64| f64::from(sides) / 2.0 * r * r * (TAU / f64::from(sides)).sin();
        for fill in [BARE, FLAT] {
            let mut e = Editable::from(&m, 0.0).unwrap();
            let rims = e.rims();
            assert_eq!(rims.len(), 2);
            let inner = rims
                .iter()
                .find(|r| e.position(e.edge(r.edges[0]).unwrap()[0]).unwrap().x.abs() <= 1.0)
                .unwrap();
            let report = e.fill(inner.id(), &fill).unwrap();
            check(&e);
            assert_eq!(report.flattened, 1);
            let area = area_of(&e, report.faces.clone());
            assert!((area - polygon(1.0)).abs() < 1e-12, "{area}");
            for v in report.vertices.clone() {
                assert_eq!(e.position(v).unwrap().z, 0.0);
            }
            let report = e.fill(e.rims()[0].id(), &fill).unwrap();
            assert_eq!(report.flattened, 1);
            check(&e);
            assert!(e.is_closed());
            assert!(e.non_manifold().is_empty());
        }
    }

    #[test]
    fn a_long_loop_whose_shadow_crosses_itself_is_refused() {
        // A ribbon along a trefoil knot, 600 quads long: each rim's shadow
        // crosses itself three times.
        let steps = 600;
        let mut positions = Vec::new();
        for lift in [0.0, 0.05] {
            for i in 0..steps {
                let t = TAU * f64::from(i) / f64::from(steps);
                let r = 3.0 + (3.0 * t).cos();
                positions.push(Point::new(
                    r * (2.0 * t).cos(),
                    r * (2.0 * t).sin(),
                    (3.0 * t).sin() + lift,
                ));
            }
        }
        let mut triangles = Vec::new();
        for i in 0..steps {
            let j = (i + 1) % steps;
            triangles.extend([[i, j, steps + j], [i, steps + j, steps + i]]);
        }
        let mut e = Editable::from(&mesh(positions, triangles), 0.0).unwrap();
        let before = e.clone();
        let rim = e.rims()[0].id();
        let err = e.fill(rim, &BARE).unwrap_err();
        assert!(matches!(err, ogeom_core::OgeomError::Construction(_)));
        assert!(err.to_string().contains("crosses itself"), "{err}");
        assert_eq!(e, before);
    }

    #[test]
    fn a_name_that_is_not_a_rim_is_refused() {
        let mut e = Editable::from(&grid(2), 0.0).unwrap();
        let rim = &e.rims()[0];
        let joined = (0..index(3 * e.face_count()))
            .find(|h| !rim.edges.contains(h))
            .unwrap();
        let later = rim.edges[1];
        let before = e.clone();
        for h in [joined, later, u32::MAX] {
            let err = e.fill(RimId(h), &BARE).unwrap_err();
            assert!(matches!(err, ogeom_core::OgeomError::Dangling(_)), "{err}");
        }
        assert_eq!(e, before);
    }

    #[test]
    fn the_fill_is_deterministic() {
        let m = sphere(2.0, 18, 22);
        let mut a = Editable::from(&m, 0.0).unwrap();
        let removed = faces_where(&a, |c| c.x + 0.4 * c.y > 1.2 || c.z < -1.5);
        a.remove_faces(&removed);
        let mut b = a.clone();
        let ra = fill_all(&mut a, &CURVED);
        let rb = fill_all(&mut b, &CURVED);
        assert_eq!(ra, rb);
        assert_eq!(a, b);
        assert!(a.is_closed());
    }
}
