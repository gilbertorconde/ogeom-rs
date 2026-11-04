//! An indexed triangle mesh held for editing, with the topology queries
//! edits need: which faces meet across each edge, the connected pieces,
//! and the open rims.
//!
//! # Half-edges
//!
//! Side `k` of face `f` is half-edge `3 f + k`, running from the face's
//! corner `k` to corner `k + 1`. Two half-edges are twins when they run one
//! edge opposite ways and the faces are joined across it; a half-edge with
//! no twin lies on a rim.
//!
//! # Points, and the vertices meeting at them
//!
//! Every vertex is a copy of a *point*: the welded vertex of the source it
//! came from. Faces meet across an edge when their corners name the same
//! two points, and are joined across it by this rule, read only from the
//! faces on that edge:
//!
//! - two faces running it opposite ways are joined;
//! - two running it the same way are not, and the edge is open on both;
//! - three or more are taken in turn about the edge, and each one joined
//!   to the next where the space between them lies behind both, as the
//!   material of a closed body does. Two closed bodies sharing an edge each
//!   close on their own; a fin is joined to neither side.
//!
//! The faces about a point, joined across the edges at it, make fans. A
//! point one fan meets is one vertex; a point several fans meet, as at the
//! centre of a bowtie or where two bodies touch, is split into one vertex
//! per fan, each a copy of the point at its position, so that every vertex
//! has one fan and every rim one way on. [`Editable::non_manifold`] reports
//! the copies, and the ends of the edges the rule could not join cleanly.
//!
//! # What edits keep
//!
//! Face and vertex numbers are stable: an edit never renumbers what it
//! leaves. A removed face keeps its number, unused; face `f` below the
//! source's triangle count is source triangle `f`. After an edit, the fans
//! about a point it touched, taken in the order of their lowest-numbered
//! half-edge, each keep the lowest-numbered vertex their corners held that
//! no fan before took, and a fan left without one takes a new copy of the
//! point, numbered after every vertex there is. So a vertex split keeps its
//! number on one side, fans joined again share the lowest number among
//! them, and undoing an edit, flipping a face back, gives back the numbers
//! it had. A vertex an edit leaves with no face is gone.
//! [`Editable::to_triangulation`] is the compaction step: it numbers the
//! vertices left in order, then the faces left in order.
//!
//! # Source vertices
//!
//! Every vertex and every copy of it remembers the source vertices it came
//! from, through every edit: [`Editable::source_vertices`]. A vertex of
//! the source has one, the vertex it was welded onto. A collapse, by
//! [`Editable::simplify`] or by [`Editable::offset`] repairing a fold,
//! leaves the vertex it keeps with the source vertices of both. A vertex
//! an edit makes rather than carries, on an edge it cut or in a hole it
//! filled, has none, and adds none to a vertex it is merged with.
//! [`Editable::source_vertex`] names the source vertex where there is
//! exactly one.
//!
//! A marking of the source's vertices selects the vertices all of whose
//! source vertices it marks, and at least one:
//! [`Editable::region_from_sources`] turns it into the per-vertex region
//! that [`Editable::offset`] and [`Editable::smooth`] take. Provenance is
//! only ever carried; nothing is matched by position.

use ogeom_core::{FastMap, OgeomResult, ogeom_bail};
use ogeom_math::{Aabb, Point, Vector};
use ogeom_topo::Triangulation;

pub mod cut;
pub mod fill;
pub mod offset;
pub mod simplify;
pub mod smooth;

pub use cut::{Cap, CutReport, Side};
pub use fill::{Fairing, Fill, FillReport, OPTIMAL_UP_TO};
pub use offset::{OFFSET_EDGE_ANGLE, OFFSET_FARTHEST, OffsetReport};
pub use simplify::SimplifyReport;
pub use smooth::{PASSES, Smooth, SmoothReport};

/// No half-edge, face or vertex.
const NONE: u32 = u32::MAX;

/// The half-edge after `h` round its face.
const fn next(h: u32) -> u32 {
    if h % 3 == 2 { h - 2 } else { h + 1 }
}

/// The half-edge before `h` round its face.
const fn prev(h: u32) -> u32 {
    if h.is_multiple_of(3) { h + 2 } else { h - 1 }
}

fn index(i: usize) -> u32 {
    u32::try_from(i).unwrap_or(NONE)
}

/// An indexed triangle mesh held as half-edges, for editing.
///
/// Built by [`Editable::from`]; read back by
/// [`Editable::to_triangulation`]. See the [module documentation](self)
/// for how faces are joined, how a point several fans meet is split, and
/// what numbers survive an edit. Edits take `&mut self` and run on one
/// thread; queries take `&self`, and the mesh is `Send` and `Sync`.
#[derive(Debug, Clone, PartialEq)]
pub struct Editable {
    positions: Vec<Point>,
    /// The point each vertex is a copy of: the vertex itself, or for a
    /// copy split off it, the vertex first given that point.
    point: Vec<u32>,
    /// The next copy of the same point, round a ring that holds them all.
    copies: Vec<u32>,
    /// The source vertex each vertex came from, or `NONE` for one an edit
    /// made; for a vertex in `merged`, the lowest of its source vertices.
    source: Vec<u32>,
    /// The source vertices, more than one and in increasing order, of each
    /// vertex a collapse left with several.
    merged: FastMap<u32, Vec<u32>>,
    /// How many vertices the source had.
    source_count: usize,
    /// A half-edge leaving each vertex: the one on a rim where there is
    /// one, so turning from it sweeps the whole fan; `NONE` for a vertex
    /// with no face.
    out: Vec<u32>,
    /// Whether each vertex was left with no face by an edit.
    gone: Vec<bool>,
    corners: Vec<[u32; 3]>,
    /// Whether each face is still there.
    live: Vec<bool>,
    twin: Vec<u32>,
    deflection_met: bool,
}

/// A connected piece of the mesh: faces joined across edges.
#[derive(Debug, Clone, PartialEq)]
pub struct Component {
    /// Its faces, in increasing order.
    pub faces: Vec<u32>,
    /// The box of its faces' corners.
    pub bounds: Aabb,
    /// The sum of its faces' areas.
    pub area: f64,
}

/// A closed loop of open edges.
#[derive(Debug, Clone, PartialEq)]
pub struct Rim {
    /// Its half-edges, each starting where the one before ended, starting
    /// from the lowest-numbered and running with the faces on the left, as
    /// the faces wind.
    pub edges: Vec<u32>,
    /// The sum of the edges' lengths.
    pub length: f64,
    /// The greatest distance between two of its vertices.
    pub span: f64,
}

/// A rim's name: its lowest-numbered half-edge, the first of its
/// [`Rim::edges`]. An edit that leaves the rim's faces as they are leaves
/// its name. One that adds faces on it, or removes or flips its faces, may
/// give it another.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RimId(pub u32);

impl Rim {
    /// The rim's name.
    #[must_use]
    pub fn id(&self) -> RimId {
        RimId(self.edges[0])
    }
}

/// Buffers reused from one point to the next while giving out vertices.
#[derive(Default)]
struct Scratch {
    parent: Vec<usize>,
    fan: Vec<usize>,
    taken: Vec<u32>,
}

/// The points whose vertices an edit relinks: every point, or some.
enum Scope {
    All(usize),
    Points(Vec<u32>),
}

impl Scope {
    fn len(&self) -> usize {
        match self {
            Self::All(n) => *n,
            Self::Points(p) => p.len(),
        }
    }

    fn point(&self, slot: usize) -> u32 {
        match self {
            Self::All(_) => index(slot),
            Self::Points(p) => p[slot],
        }
    }

    fn slot(&self, point: u32) -> Option<usize> {
        match self {
            Self::All(n) => ((point as usize) < *n).then_some(point as usize),
            Self::Points(p) => p.binary_search(&point).ok(),
        }
    }
}

impl Editable {
    /// The mesh as half-edges, its vertices first welded within `weld`.
    ///
    /// With `weld` zero the mesh's own vertex numbers are the points, and a
    /// manifold mesh reads back from [`Editable::to_triangulation`] with the
    /// same vertices and the same triangles in the same order. With `weld`
    /// positive, vertices within `weld` of each other become one: each
    /// vertex in turn joins the lowest-numbered kept vertex within `weld`,
    /// found through a grid of cells `weld` wide, or is kept itself at its
    /// own position. A triangle naming one vertex twice, or two of whose
    /// corners became one, is removed, its face number unused.
    ///
    /// Normals and surface parameters are not carried: an edit leaves
    /// neither true.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if `weld` is
    /// negative or not finite, or a position is not finite;
    /// [`OgeomError::Range`](ogeom_core::OgeomError::Range) if a triangle
    /// names a vertex the mesh does not have, or the mesh has more
    /// half-edges and vertices than `u32` numbers.
    pub fn from(mesh: &Triangulation, weld: f64) -> OgeomResult<Self> {
        if !weld.is_finite() || weld < 0.0 {
            ogeom_bail!(Numeric, "a weld distance of {weld} is not a distance");
        }
        if let Some(k) = mesh.positions.iter().position(|p| !p.is_finite()) {
            ogeom_bail!(Numeric, "mesh vertex {k} is not finite");
        }
        let count = mesh.positions.len();
        if let Some(&k) = mesh
            .triangles
            .iter()
            .flatten()
            .find(|&&k| k as usize >= count)
        {
            ogeom_bail!(Range, "a triangle names vertex {k} of {count}");
        }
        // Every split copies a corner, so the vertices never outnumber
        // the positions and the corners together.
        if count.saturating_add(mesh.triangles.len().saturating_mul(3)) >= NONE as usize {
            ogeom_bail!(
                Range,
                "a mesh of {count} vertices and {} triangles is past u32 numbering",
                mesh.triangles.len()
            );
        }

        let (positions, remap, source) = if weld > 0.0 {
            weld_points(&mesh.positions, weld)
        } else {
            let all: Vec<u32> = (0..count).map(index).collect();
            (mesh.positions.clone(), all.clone(), all)
        };
        let n = positions.len();
        let corners: Vec<[u32; 3]> = mesh
            .triangles
            .iter()
            .map(|t| t.map(|v| remap[v as usize]))
            .collect();
        let live = corners
            .iter()
            .map(|&[a, b, c]| a != b && b != c && c != a)
            .collect();
        let mut mesh = Self {
            positions,
            point: (0..n).map(index).collect(),
            copies: (0..n).map(index).collect(),
            source,
            merged: FastMap::default(),
            source_count: count,
            out: vec![NONE; n],
            gone: vec![false; n],
            twin: vec![NONE; corners.len() * 3],
            corners,
            live,
            deflection_met: mesh.deflection_met,
        };
        let faces: Vec<u32> = (0..mesh.corners.len())
            .filter(|&f| mesh.live[f])
            .map(index)
            .collect();
        mesh.relink(&Scope::All(n), &faces);
        Ok(mesh)
    }

    /// The mesh's vertices and faces, those an edit left, renumbered in
    /// order. A vertex no face of the source used is kept.
    ///
    /// A vertex no edit touched keeps its position; a copy split off a
    /// point stands at the point's position, so two vertices of the result
    /// may coincide.
    #[must_use]
    pub fn to_triangulation(&self) -> Triangulation {
        let mut renumber = vec![NONE; self.positions.len()];
        let mut positions = Vec::with_capacity(self.positions.len());
        for (v, &p) in self.positions.iter().enumerate() {
            if !self.gone[v] {
                renumber[v] = index(positions.len());
                positions.push(p);
            }
        }
        let triangles = self
            .corners
            .iter()
            .zip(&self.live)
            .filter(|&(_, &live)| live)
            .map(|(t, _)| t.map(|v| renumber[v as usize]))
            .collect();
        Triangulation {
            positions,
            triangles,
            deflection_met: self.deflection_met,
            ..Triangulation::new()
        }
    }

    /// How many vertex numbers there are, including those of vertices
    /// gone.
    #[must_use]
    pub fn vertex_count(&self) -> usize {
        self.positions.len()
    }

    /// How many face numbers there are, including those of faces removed.
    #[must_use]
    pub fn face_count(&self) -> usize {
        self.corners.len()
    }

    /// The position of vertex `v`, if it is there.
    #[must_use]
    pub fn position(&self, v: u32) -> Option<Point> {
        let k = v as usize;
        (k < self.positions.len() && !self.gone[k]).then(|| self.positions[k])
    }

    /// The corners of face `f`, counter-clockwise about its outward
    /// normal, if it is there.
    #[must_use]
    pub fn face(&self, f: u32) -> Option<[u32; 3]> {
        let k = f as usize;
        (k < self.corners.len() && self.live[k]).then(|| self.corners[k])
    }

    /// The vertices half-edge `h` runs from and to, if its face is there.
    #[must_use]
    pub fn edge(&self, h: u32) -> Option<[u32; 2]> {
        self.face(h / 3)
            .is_some()
            .then(|| [self.tail(h), self.head(h)])
    }

    /// The source vertex vertex `v` came from, where it came from exactly
    /// one: the one it was welded onto, for every copy split off it and
    /// through every edit. `None` for a vertex an edit made, one a
    /// collapse left with several source vertices, or a number past the
    /// last.
    #[must_use]
    pub fn source_vertex(&self, v: u32) -> Option<u32> {
        match self.source_vertices(v) {
            &[s] => Some(s),
            _ => None,
        }
    }

    /// The source vertices vertex `v` came from, in increasing order: one
    /// for a vertex carried from the source, every one merged into it for
    /// a vertex a collapse kept, and none for a vertex an edit made or a
    /// number past the last. See the [module documentation](self).
    #[must_use]
    pub fn source_vertices(&self, v: u32) -> &[u32] {
        let k = v as usize;
        if k >= self.source.len() || self.source[k] == NONE {
            return &[];
        }
        self.merged
            .get(&v)
            .map_or_else(|| core::slice::from_ref(&self.source[k]), Vec::as_slice)
    }

    /// How many vertices the source had: the length of a marking of them.
    /// [`Editable::simplify`] keeps the count of the mesh first built.
    #[must_use]
    pub fn source_count(&self) -> usize {
        self.source_count
    }

    /// The region, one entry per vertex number, that a marking of the
    /// source's vertices selects: the vertices that have at least one
    /// source vertex, and every one of them marked. A vertex an edit made,
    /// with none, and a vertex gone, are not selected. The result is what
    /// [`Editable::offset`] and [`Smooth::region`] take.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Range`](ogeom_core::OgeomError::Range) if `marked`
    /// does not have one entry per source vertex, as
    /// [`Editable::source_count`] counts them.
    pub fn region_from_sources(&self, marked: &[bool]) -> OgeomResult<Vec<bool>> {
        if marked.len() != self.source_count {
            ogeom_bail!(
                Range,
                "a marking of {} entries for {} source vertices",
                marked.len(),
                self.source_count
            );
        }
        Ok((0..self.positions.len())
            .map(|k| {
                let sources = self.source_vertices(index(k));
                !self.gone[k] && !sources.is_empty() && sources.iter().all(|&s| marked[s as usize])
            })
            .collect())
    }

    /// Give vertex `to` the source vertices of vertex `from`.
    fn copy_sources(&mut self, to: u32, from: u32) {
        self.source[to as usize] = self.source[from as usize];
        match self.merged.get(&from).cloned() {
            Some(set) => self.merged.insert(to, set),
            None => self.merged.remove(&to),
        };
    }

    /// Give vertex `to` the source vertices of vertex `from` as well as
    /// its own.
    fn absorb_sources(&mut self, to: u32, from: u32) {
        let mut set: Vec<u32> = self.source_vertices(to).to_vec();
        set.extend_from_slice(self.source_vertices(from));
        set.sort_unstable();
        set.dedup();
        self.set_sources(to, set);
    }

    /// Set vertex `v`'s source vertices, sorted and once each.
    fn set_sources(&mut self, v: u32, set: Vec<u32>) {
        self.source[v as usize] = set.first().copied().unwrap_or(NONE);
        if set.len() > 1 {
            self.merged.insert(v, set);
        } else {
            self.merged.remove(&v);
        }
    }

    /// Whether there is a face and every edge of every face is joined to
    /// another face.
    #[must_use]
    pub fn is_closed(&self) -> bool {
        let mut any = false;
        for (f, &live) in self.live.iter().enumerate() {
            if live {
                any = true;
                if self.twin[3 * f..3 * f + 3].contains(&NONE) {
                    return false;
                }
            }
        }
        any
    }

    /// The vertices where the faces as they stand are not a manifold, in
    /// increasing order: each copy of a point several fans meet, and both
    /// ends of each edge that more than two faces use, or that two faces
    /// run the same way.
    #[must_use]
    pub fn non_manifold(&self) -> Vec<u32> {
        let mut flagged = vec![false; self.positions.len()];
        for v in 0..self.positions.len() {
            if self.out[v] == NONE || flagged[v] || self.copies[v] as usize == v {
                continue;
            }
            let ring = self.ring(index(v));
            let with_faces: Vec<u32> = ring
                .into_iter()
                .filter(|&w| self.out[w as usize] != NONE)
                .collect();
            if with_faces.len() > 1 {
                for w in with_faces {
                    flagged[w as usize] = true;
                }
            }
        }
        let mut edges: Vec<(u64, u32)> = Vec::new();
        for (f, &live) in self.live.iter().enumerate() {
            if live {
                for h in index(3 * f)..index(3 * f + 3) {
                    edges.push((self.edge_key(h), h));
                }
            }
        }
        edges.sort_unstable();
        for run in edges.chunk_by(|a, b| a.0 == b.0) {
            let clean =
                run.len() == 1 || (run.len() == 2 && self.twin[run[0].1 as usize] == run[1].1);
            if !clean {
                for &(_, h) in run {
                    flagged[self.tail(h) as usize] = true;
                    flagged[self.head(h) as usize] = true;
                }
            }
        }
        (0..self.positions.len())
            .filter(|&v| flagged[v])
            .map(index)
            .collect()
    }

    /// The connected pieces, faces joined across edges, in the order of
    /// their lowest-numbered faces.
    #[must_use]
    pub fn components(&self) -> Vec<Component> {
        let mut seen = vec![false; self.corners.len()];
        let mut pieces = Vec::new();
        let mut stack = Vec::new();
        for seed in 0..self.corners.len() {
            if !self.live[seed] || seen[seed] {
                continue;
            }
            seen[seed] = true;
            stack.push(index(seed));
            let mut faces = Vec::new();
            let (mut lo, mut hi) = (
                Point::new(f64::MAX, f64::MAX, f64::MAX),
                Point::new(f64::MIN, f64::MIN, f64::MIN),
            );
            let mut area = 0.0;
            while let Some(f) = stack.pop() {
                faces.push(f);
                let [a, b, c] = self.corners[f as usize].map(|v| self.positions[v as usize]);
                lo = lo.min(a).min(b).min(c);
                hi = hi.max(a).max(b).max(c);
                area += (b - a).cross(c - a).magnitude() * 0.5;
                for h in 3 * f..3 * f + 3 {
                    let g = self.twin[h as usize];
                    if g != NONE && !seen[(g / 3) as usize] {
                        seen[(g / 3) as usize] = true;
                        stack.push(g / 3);
                    }
                }
            }
            faces.sort_unstable();
            pieces.push(Component {
                faces,
                bounds: Aabb::of_corners(lo, hi),
                area,
            });
        }
        pieces
    }

    /// The rims: the loops of half-edges with no twin, in the order of
    /// their lowest-numbered half-edges.
    #[must_use]
    pub fn rims(&self) -> Vec<Rim> {
        let mut seen = vec![false; self.twin.len()];
        let mut rims = Vec::new();
        for start in 0..self.twin.len() {
            if self.twin[start] != NONE || seen[start] || !self.live[start / 3] {
                continue;
            }
            let mut edges = Vec::new();
            let mut points = Vec::new();
            let mut length = 0.0;
            let mut h = index(start);
            loop {
                seen[h as usize] = true;
                edges.push(h);
                let (a, b) = (self.tail(h), self.head(h));
                let p = self.positions[a as usize];
                points.push(p);
                length += p.distance(self.positions[b as usize]);
                h = self.rim_after(h);
                if h as usize == start {
                    break;
                }
            }
            rims.push(Rim {
                edges,
                length,
                span: diameter(&points),
            });
        }
        rims
    }

    /// Remove faces, and split each vertex whose faces left fall into
    /// several fans into one vertex per fan, so what is left stays a
    /// manifold. Faces not there are passed over.
    pub fn remove_faces(&mut self, faces: &[u32]) {
        let removed = self.listed(faces);
        if removed.is_empty() {
            return;
        }
        let (scope, around) = self.around(&removed);
        for &f in &removed {
            self.live[f as usize] = false;
            for h in 3 * f..3 * f + 3 {
                self.twin[h as usize] = NONE;
            }
        }
        let kept: Vec<u32> = around
            .into_iter()
            .filter(|f| removed.binary_search(f).is_err())
            .collect();
        self.relink(&scope, &kept);
    }

    /// Reverse the winding of faces, and join or part the faces about
    /// them again by the same rule as [`Editable::from`]: an edge two
    /// faces now run the same way opens, and one they now run opposite
    /// ways closes. Each face is flipped once however often it is listed;
    /// faces not there are passed over.
    pub fn flip_faces(&mut self, faces: &[u32]) {
        let flipped = self.listed(faces);
        if flipped.is_empty() {
            return;
        }
        let (scope, around) = self.around(&flipped);
        for &f in &flipped {
            self.corners[f as usize].swap(1, 2);
        }
        self.relink(&scope, &around);
    }

    /// The faces listed that are there, sorted and once each.
    fn listed(&self, faces: &[u32]) -> Vec<u32> {
        let mut listed: Vec<u32> = faces
            .iter()
            .copied()
            .filter(|&f| (f as usize) < self.corners.len() && self.live[f as usize])
            .collect();
        listed.sort_unstable();
        listed.dedup();
        listed
    }

    /// The points at the corners of `faces`, and every face about any of
    /// them, sorted.
    fn around(&self, faces: &[u32]) -> (Scope, Vec<u32>) {
        let mut points: Vec<u32> = faces
            .iter()
            .flat_map(|&f| self.corners[f as usize])
            .map(|v| self.point[v as usize])
            .collect();
        points.sort_unstable();
        points.dedup();
        let around = self.about(&points);
        (Scope::Points(points), around)
    }

    /// Every face with a corner at one of `points`, sorted.
    fn about(&self, points: &[u32]) -> Vec<u32> {
        let mut around = Vec::new();
        for &p in points {
            for v in self.ring(p) {
                let start = self.out[v as usize];
                if start == NONE {
                    continue;
                }
                let mut h = start;
                loop {
                    around.push(h / 3);
                    h = self.twin[prev(h) as usize];
                    if h == NONE || h == start {
                        break;
                    }
                }
            }
        }
        around.sort_unstable();
        around.dedup();
        around
    }

    /// The copies of the point `p`, starting from `p`.
    fn ring(&self, p: u32) -> Vec<u32> {
        let mut ring = vec![p];
        let mut v = self.copies[p as usize];
        while v != p {
            ring.push(v);
            v = self.copies[v as usize];
        }
        ring
    }

    fn tail(&self, h: u32) -> u32 {
        self.corners[(h / 3) as usize][(h % 3) as usize]
    }

    fn head(&self, h: u32) -> u32 {
        self.tail(next(h))
    }

    /// The edge `h` lies on, as its two points, smaller first.
    fn edge_key(&self, h: u32) -> u64 {
        let (a, b) = (
            self.point[self.tail(h) as usize],
            self.point[self.head(h) as usize],
        );
        (u64::from(a.min(b)) << 32) | u64::from(a.max(b))
    }

    /// The rim half-edge leaving the vertex rim half-edge `h` ends at:
    /// turning about that vertex through its fan, from `h`'s face to the
    /// fan's other side.
    fn rim_after(&self, h: u32) -> u32 {
        let mut g = next(h);
        while self.twin[g as usize] != NONE {
            g = next(self.twin[g as usize]);
        }
        g
    }

    /// Join the half-edges on every edge at a point in `scope`, then give
    /// each fan about such a point its vertex. `faces` must be every live
    /// face with a corner at one of those points, sorted.
    fn relink(&mut self, scope: &Scope, faces: &[u32]) {
        // The half-edges leaving each point in scope, in increasing order:
        // one per corner there.
        let mut start = vec![0_usize; scope.len() + 1];
        for &f in faces {
            for v in self.corners[f as usize] {
                if let Some(s) = scope.slot(self.point[v as usize]) {
                    start[s + 1] += 1;
                }
            }
        }
        for s in 0..scope.len() {
            start[s + 1] += start[s];
        }
        let mut fill = start.clone();
        let mut leaving = vec![0_u32; start[scope.len()]];
        for &f in faces {
            for (k, v) in self.corners[f as usize].into_iter().enumerate() {
                if let Some(s) = scope.slot(self.point[v as usize]) {
                    leaving[fill[s]] = 3 * f + index(k);
                    fill[s] += 1;
                }
            }
        }

        let mut on_edges: Vec<(u32, u32)> = Vec::new();
        for s in 0..scope.len() {
            let p = scope.point(s);
            on_edges.clear();
            for &h in &leaving[start[s]..start[s + 1]] {
                let into = prev(h);
                on_edges.push((self.point[self.head(h) as usize], h));
                on_edges.push((self.point[self.tail(into) as usize], into));
            }
            on_edges.sort_unstable();
            for run in on_edges.chunk_by(|a, b| a.0 == b.0) {
                let other = run[0].0;
                // An edge with both ends in scope is joined from its
                // smaller end.
                if other < p && scope.slot(other).is_some() {
                    continue;
                }
                self.join(run, p.min(other));
            }
        }

        let mut scratch = Scratch::default();
        for s in 0..scope.len() {
            let fan = &leaving[start[s]..start[s + 1]];
            self.give_vertices(scope.point(s), fan, &mut scratch);
        }
    }

    /// Join the half-edges of one edge, each the second of a pair in `run`,
    /// `low` the smaller of the edge's two points, by the rule in the
    /// module documentation.
    fn join(&mut self, run: &[(u32, u32)], low: u32) {
        for &(_, h) in run {
            self.twin[h as usize] = NONE;
        }
        let rising = |me: &Self, h: u32| me.point[me.tail(h) as usize] == low;
        match *run {
            [_] => {}
            [(_, a), (_, b)] => {
                if rising(self, a) != rising(self, b) {
                    self.twin[a as usize] = b;
                    self.twin[b as usize] = a;
                }
            }
            _ => {
                // About the edge's axis, from its low end to its high end,
                // each face is a half-plane at an angle. A face running the
                // edge from low to high faces the way the angle grows, so
                // its material lies the way it falls; one running it from
                // high to low has its material the way the angle grows.
                let h0 = run[0].1;
                let (from, to) = if rising(self, h0) {
                    (self.tail(h0), self.head(h0))
                } else {
                    (self.head(h0), self.tail(h0))
                };
                let from = self.positions[from as usize];
                let axis = self.positions[to as usize] - from;
                let length = axis.magnitude();
                let axis = if length > 0.0 { axis / length } else { axis };
                let across: Vec<Vector> = run
                    .iter()
                    .map(|&(_, h)| {
                        let w = self.positions[self.tail(prev(h)) as usize] - from;
                        w - axis * w.dot(axis)
                    })
                    .collect();
                let mut turn: Vec<(f64, u32, bool)> = match across
                    .iter()
                    .find(|w| w.magnitude() > 0.0)
                {
                    Some(&first) => {
                        let e1 = first / first.magnitude();
                        let e2 = axis.cross(e1);
                        run.iter()
                            .zip(&across)
                            .map(|(&(_, h), w)| (w.dot(e2).atan2(w.dot(e1)), h, rising(self, h)))
                            .collect()
                    }
                    None => run
                        .iter()
                        .map(|&(_, h)| (0.0, h, rising(self, h)))
                        .collect(),
                };
                turn.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
                let n = turn.len();
                for i in 0..n {
                    let (a, b) = (turn[i], turn[(i + 1) % n]);
                    if !a.2 && b.2 {
                        self.twin[a.1 as usize] = b.1;
                        self.twin[b.1 as usize] = a.1;
                    }
                }
            }
        }
    }

    /// Group the half-edges `leaving` point `p`, in increasing order, into
    /// fans joined across the edges at `p`, and give each fan a vertex: the
    /// lowest-numbered one its corners held that no earlier fan took, or a
    /// new copy of the point. Copies of `p` no fan took, and that had
    /// faces, are gone.
    fn give_vertices(&mut self, p: u32, leaving: &[u32], scratch: &mut Scratch) {
        fn root(parent: &mut [usize], mut i: usize) -> usize {
            while parent[i] != i {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            i
        }
        let Scratch { parent, fan, taken } = scratch;
        parent.clear();
        parent.extend(0..leaving.len());
        for (i, &h) in leaving.iter().enumerate() {
            let g = self.twin[h as usize];
            if g == NONE {
                continue;
            }
            if let Ok(j) = leaving.binary_search(&next(g)) {
                let (a, b) = (root(parent, i), root(parent, j));
                parent[a.max(b)] = a.min(b);
            }
        }
        // Rooted always at the smaller, a fan's root is its first
        // half-edge, so the fans are numbered in the order of their first.
        fan.clear();
        let mut fans = 0;
        for i in 0..leaving.len() {
            let r = root(parent, i);
            if r == i {
                fan.push(fans);
                fans += 1;
            } else {
                let k = fan[r];
                fan.push(k);
            }
        }

        taken.clear();
        for k in 0..fans {
            let members = || {
                leaving
                    .iter()
                    .zip(fan.iter())
                    .filter(move |&(_, &f)| f == k)
                    .map(|(&h, _)| h)
            };
            let held = members()
                .map(|h| self.tail(h))
                .filter(|v| !taken.contains(v))
                .min();
            let first = members().next().unwrap_or(NONE);
            let v = held.unwrap_or_else(|| self.copy_of(p, self.tail(first)));
            taken.push(v);
            let mut out = NONE;
            for h in members() {
                self.corners[(h / 3) as usize][(h % 3) as usize] = v;
                if out == NONE && self.twin[h as usize] == NONE {
                    out = h;
                }
            }
            self.out[v as usize] = if out == NONE { first } else { out };
            self.gone[v as usize] = false;
        }
        let mut v = p;
        loop {
            let k = v as usize;
            if !taken.contains(&v) && self.out[k] != NONE {
                self.out[k] = NONE;
                self.gone[k] = true;
            }
            v = self.copies[k];
            if v == p {
                break;
            }
        }
    }

    /// A new vertex copying point `p`, standing where vertex `like` does
    /// and remembering its source.
    fn copy_of(&mut self, p: u32, like: u32) -> u32 {
        let v = index(self.positions.len());
        self.positions.push(self.positions[like as usize]);
        self.point.push(p);
        self.source.push(NONE);
        self.copy_sources(v, like);
        self.out.push(NONE);
        self.gone.push(false);
        self.copies.push(self.copies[p as usize]);
        self.copies[p as usize] = v;
        v
    }
}

/// Weld points within `distance`: each in turn joins the lowest-numbered
/// kept point within `distance`, looked for in the cells `distance` wide
/// about it, or is kept. Returns the kept points, the kept point each
/// point became, and the point each kept one was.
fn weld_points(positions: &[Point], distance: f64) -> (Vec<Point>, Vec<u32>, Vec<u32>) {
    // A cast saturates, so a coordinate past `i64`'s cells shares the last
    // one and is still compared by distance.
    #[allow(clippy::cast_possible_truncation, reason = "saturating")]
    let cell = |p: Point| {
        (
            (p.x / distance).floor() as i64,
            (p.y / distance).floor() as i64,
            (p.z / distance).floor() as i64,
        )
    };
    let mut grid: FastMap<(i64, i64, i64), Vec<u32>> =
        FastMap::with_capacity_and_hasher(positions.len(), Default::default());
    let mut kept: Vec<Point> = Vec::with_capacity(positions.len());
    let mut founders = Vec::with_capacity(positions.len());
    let mut remap = Vec::with_capacity(positions.len());
    for (i, p) in positions.iter().enumerate() {
        let (x, y, z) = cell(*p);
        let mut found = NONE;
        for dx in -1..=1_i64 {
            for dy in -1..=1_i64 {
                for dz in -1..=1_i64 {
                    let key = (
                        x.saturating_add(dx),
                        y.saturating_add(dy),
                        z.saturating_add(dz),
                    );
                    if let Some(list) = grid.get(&key)
                        && let Some(&k) = list
                            .iter()
                            .find(|&&k| kept[k as usize].distance(*p) <= distance)
                    {
                        found = found.min(k);
                    }
                }
            }
        }
        if found == NONE {
            found = index(kept.len());
            kept.push(*p);
            founders.push(index(i));
            grid.entry((x, y, z)).or_default().push(found);
        }
        remap.push(found);
    }
    (kept, remap, founders)
}

/// The greatest distance between two of `points`, by branch and bound over
/// a tree of boxes: a pair of boxes is opened only while the farthest their
/// corners stand apart exceeds the best pair found, which, the rounded
/// differences of coordinates growing with the coordinates, bounds every
/// pair of points within them as computed.
fn diameter(points: &[Point]) -> f64 {
    struct Node {
        lo: Point,
        hi: Point,
        range: (usize, usize),
        children: Option<(usize, usize)>,
    }
    fn build(points: &mut [Point], at: usize, nodes: &mut Vec<Node>) -> usize {
        let lo = points.iter().fold(points[0], |a, &p| a.min(p));
        let hi = points.iter().fold(points[0], |a, &p| a.max(p));
        let k = nodes.len();
        nodes.push(Node {
            lo,
            hi,
            range: (at, at + points.len()),
            children: None,
        });
        // Leaves of a few points, compared pairwise.
        if points.len() > 8 {
            let size = hi - lo;
            let axis = if size.x >= size.y && size.x >= size.z {
                0
            } else if size.y >= size.z {
                1
            } else {
                2
            };
            let half = points.len() / 2;
            let along = |p: &Point| [p.x, p.y, p.z][axis];
            points.select_nth_unstable_by(half, |a, b| along(a).total_cmp(&along(b)));
            let (left, right) = points.split_at_mut(half);
            let a = build(left, at, nodes);
            let b = build(right, at + half, nodes);
            nodes[k].children = Some((a, b));
        }
        k
    }
    fn reach(a: &Node, b: &Node) -> f64 {
        let span =
            |lo1: f64, hi1: f64, lo2: f64, hi2: f64| (hi1 - lo2).abs().max((hi2 - lo1).abs());
        Vector::new(
            span(a.lo.x, a.hi.x, b.lo.x, b.hi.x),
            span(a.lo.y, a.hi.y, b.lo.y, b.hi.y),
            span(a.lo.z, a.hi.z, b.lo.z, b.hi.z),
        )
        .magnitude()
    }

    if points.len() < 2 {
        return 0.0;
    }
    let mut points = points.to_vec();
    let mut nodes = Vec::new();
    build(&mut points, 0, &mut nodes);
    // A good first pair: the point farthest from the first, and the point
    // farthest from that.
    let farthest = |from: Point| {
        points
            .iter()
            .copied()
            .max_by(|a, b| a.distance(from).total_cmp(&b.distance(from)))
            .unwrap_or(from)
    };
    let a = farthest(points[0]);
    let mut best = a.distance(farthest(a));
    let mut stack = vec![(0_usize, 0_usize)];
    while let Some((i, j)) = stack.pop() {
        let (a, b) = (&nodes[i], &nodes[j]);
        if reach(a, b) <= best {
            continue;
        }
        match (a.children, b.children) {
            (None, None) => {
                for p in &points[a.range.0..a.range.1] {
                    for q in &points[b.range.0..b.range.1] {
                        best = best.max(p.distance(*q));
                    }
                }
            }
            (Some((l, r)), None) => stack.extend([(l, j), (r, j)]),
            (None, Some((l, r))) => stack.extend([(i, l), (i, r)]),
            (Some((l, r)), Some(_)) if i == j => {
                stack.extend([(l, l), (l, r), (r, r)]);
            }
            (Some((l, r)), Some((m, s))) => {
                if a.range.1 - a.range.0 >= b.range.1 - b.range.0 {
                    stack.extend([(l, j), (r, j)]);
                } else {
                    stack.extend([(i, m), (i, s)]);
                }
            }
        }
    }
    best
}

#[cfg(test)]
mod tests {
    #![allow(
        clippy::unwrap_used,
        clippy::cast_possible_truncation,
        reason = "test code; the f32 casts are the rounding under test"
    )]
    use super::*;
    use core::f64::consts::{SQRT_2, TAU};
    use proptest::prelude::*;

    pub(super) fn mesh(positions: Vec<Point>, triangles: Vec<[u32; 3]>) -> Triangulation {
        Triangulation {
            positions,
            triangles,
            ..Triangulation::new()
        }
    }

    /// A torus about `z`, `rings` by `sides` quads of two triangles each,
    /// wound outward.
    pub(super) fn torus(rings: u32, sides: u32) -> Triangulation {
        let mut positions = Vec::new();
        for i in 0..rings {
            let u = f64::from(i) / f64::from(rings) * TAU;
            for j in 0..sides {
                let v = f64::from(j) / f64::from(sides) * TAU;
                let r = 3.0 + v.cos();
                positions.push(Point::new(r * u.cos(), r * u.sin(), v.sin()));
            }
        }
        let at = |i: u32, j: u32| (i % rings) * sides + (j % sides);
        let mut triangles = Vec::new();
        for i in 0..rings {
            for j in 0..sides {
                let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
                triangles.push([a, b, c]);
                triangles.push([a, c, d]);
            }
        }
        mesh(positions, triangles)
    }

    /// A square of `n` by `n` unit cells in `z = 0`, facing `+z`; cell
    /// `(i, j)` is faces `2 (n j + i)` and the one after.
    pub(super) fn grid(n: u32) -> Triangulation {
        let mut positions = Vec::new();
        for j in 0..=n {
            for i in 0..=n {
                positions.push(Point::new(f64::from(i), f64::from(j), 0.0));
            }
        }
        let at = |i: u32, j: u32| j * (n + 1) + i;
        let mut triangles = Vec::new();
        for j in 0..n {
            for i in 0..n {
                triangles.push([at(i, j), at(i + 1, j), at(i + 1, j + 1)]);
                triangles.push([at(i, j), at(i + 1, j + 1), at(i, j + 1)]);
            }
        }
        mesh(positions, triangles)
    }

    /// The tetrahedron on four points, each face wound away from the
    /// corner opposite it, its corners numbered from `first`.
    fn tetrahedron(points: [Point; 4], first: u32) -> Vec<[u32; 3]> {
        [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
            .into_iter()
            .map(|[a, b, c]| {
                let opposite = 6 - a - b - c;
                let [pa, pb, pc] = [a, b, c].map(|k| points[k]);
                let normal = (pb - pa).cross(pc - pa);
                let t = if normal.dot(points[opposite] - pa) < 0.0 {
                    [a, b, c]
                } else {
                    [a, c, b]
                };
                t.map(|k| first + index(k))
            })
            .collect()
    }

    /// A cube of side `n h` as six separate grids of `n` by `n` cells, each
    /// with its own vertices and wound outward.
    pub(super) fn cube_of_grids(n: u32, h: f64) -> Triangulation {
        let side = f64::from(n) * h;
        let mut out = Triangulation::new();
        // Each face as an origin and two in-plane steps whose cross product
        // points out of the cube.
        let faces = [
            ([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]),
            ([0.0, 0.0, side], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
            ([0.0, side, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
            ([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]),
            ([side, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
        ];
        for (o, u, v) in faces {
            let mut g = grid(n);
            for p in &mut g.positions {
                let (a, b) = (p.x * h, p.y * h);
                *p = Point::new(
                    o[0] + a * u[0] + b * v[0],
                    o[1] + a * u[1] + b * v[1],
                    o[2] + a * u[2] + b * v[2],
                );
            }
            out.append(&g);
        }
        out
    }

    /// Each triangle given its own three vertices.
    fn soup(m: &Triangulation) -> Triangulation {
        let mut out = Triangulation::new();
        for t in &m.triangles {
            let first = index(out.positions.len());
            out.positions.extend(t.map(|v| m.positions[v as usize]));
            out.triangles.push([first, first + 1, first + 2]);
        }
        out
    }

    /// The half-edge structure's invariants: twins run one edge opposite
    /// ways between live faces, every vertex with faces has exactly one fan
    /// and its `out` sweeps it from the rim, and copies share a ring.
    pub(super) fn check(e: &Editable) {
        let mut at: Vec<Vec<u32>> = vec![Vec::new(); e.positions.len()];
        for f in 0..e.corners.len() {
            if !e.live[f] {
                continue;
            }
            for h in index(3 * f)..index(3 * f + 3) {
                let v = e.tail(h);
                assert!(!e.gone[v as usize], "face {f} names gone vertex {v}");
                at[v as usize].push(h);
                let g = e.twin[h as usize];
                if g != NONE {
                    assert!(e.live[(g / 3) as usize]);
                    assert_eq!(e.twin[g as usize], h);
                    assert_eq!((e.tail(g), e.head(g)), (e.head(h), e.tail(h)));
                }
            }
        }
        for (v, leaving) in at.iter().enumerate() {
            let start = e.out[v];
            if leaving.is_empty() {
                assert_eq!(start, NONE, "vertex {v} has no face");
                continue;
            }
            if leaving.iter().any(|&h| e.twin[h as usize] == NONE) {
                assert_eq!(
                    e.twin[start as usize], NONE,
                    "vertex {v} starts off its rim"
                );
            }
            let mut swept = Vec::new();
            let mut h = start;
            loop {
                swept.push(h);
                h = e.twin[prev(h) as usize];
                if h == NONE || h == start {
                    break;
                }
            }
            swept.sort_unstable();
            assert_eq!(&swept, leaving, "vertex {v} has more than one fan");
            assert!(e.ring(e.point[v]).contains(&index(v)));
        }
    }

    fn reads_back(e: &Editable, m: &Triangulation) {
        let back = e.to_triangulation();
        assert_eq!(back.positions, m.positions);
        assert_eq!(back.triangles, m.triangles);
    }

    #[test]
    fn a_manifold_mesh_reads_back_unchanged() {
        let closed = torus(12, 8);
        let e = Editable::from(&closed, 0.0).unwrap();
        check(&e);
        reads_back(&e, &closed);
        assert!(e.is_closed());
        assert!(e.rims().is_empty());
        assert!(e.non_manifold().is_empty());
        assert_eq!(e.components().len(), 1);
        assert!((0..e.vertex_count()).all(|v| e.source_vertex(index(v)) == Some(index(v))));

        let open = grid(5);
        let e = Editable::from(&open, 0.0).unwrap();
        check(&e);
        reads_back(&e, &open);
        assert!(!e.is_closed());
        assert!(e.non_manifold().is_empty());
        let rims = e.rims();
        assert_eq!(rims.len(), 1);
        assert_eq!(rims[0].edges.len(), 20);
        // Unit edges and half-unit triangles sum without rounding.
        assert_eq!(rims[0].length, 20.0);
        assert_eq!(rims[0].span, 50.0_f64.sqrt());
        let pieces = e.components();
        assert_eq!(pieces.len(), 1);
        assert_eq!(pieces[0].area, 25.0);
    }

    #[test]
    fn three_islands_are_three_components_with_their_bounds() {
        let mut positions = Vec::new();
        let mut triangles = Vec::new();
        let corners = [[0.0, 0.0, 0.0], [5.0, 1.0, -2.0], [-4.0, 3.0, 7.0]];
        let sizes = [1.0, 2.0, 0.5];
        for (o, s) in corners.into_iter().zip(sizes) {
            let o = Point::new(o[0], o[1], o[2]);
            let points = [
                o,
                o + Vector::new(s, 0.0, 0.0),
                o + Vector::new(0.0, s, 0.0),
                o + Vector::new(0.0, 0.0, s),
            ];
            triangles.extend(tetrahedron(points, index(positions.len())));
            positions.extend(points);
        }
        let m = mesh(positions, triangles);
        let e = Editable::from(&m, 0.0).unwrap();
        check(&e);
        assert!(e.is_closed());
        let pieces = e.components();
        assert_eq!(pieces.len(), 3);
        for (k, (piece, (o, s))) in pieces
            .iter()
            .zip(corners.into_iter().zip(sizes))
            .enumerate()
        {
            let k = index(k);
            assert_eq!(piece.faces, vec![4 * k, 4 * k + 1, 4 * k + 2, 4 * k + 3]);
            let o = Point::new(o[0], o[1], o[2]);
            assert_eq!(piece.bounds, Aabb::of_corners(o, o + Vector::new(s, s, s)));
            // Three right triangles and an equilateral one of side s√2,
            // to the rounding of four cross products.
            let area = 1.5 * s * s + 3.0_f64.sqrt() / 2.0 * s * s;
            assert!((piece.area - area).abs() < 1e-14 * area);
        }
    }

    #[test]
    fn removing_a_patch_leaves_one_rim_round_it() {
        let (n, h) = (6, 0.5);
        let side = f64::from(n) * h;
        let cube = cube_of_grids(n, h);
        let mut e = Editable::from(&cube, 1e-9).unwrap();
        check(&e);
        assert!(e.is_closed());
        assert_eq!(e.vertex_count(), 6 * 36 + 2);
        assert!(e.non_manifold().is_empty());
        // The top's middle two by two cells.
        let patch: Vec<u32> = (0..e.face_count())
            .map(index)
            .filter(|&f| {
                let ps = e.face(f).unwrap().map(|v| e.position(v).unwrap());
                let c = (ps[0] + (ps[1] - ps[0]) / 3.0) + (ps[2] - ps[0]) / 3.0;
                ps.iter().all(|p| p.z == side)
                    && (c.x - side / 2.0).abs() < h
                    && (c.y - side / 2.0).abs() < h
            })
            .collect();
        assert_eq!(patch.len(), 8);
        e.remove_faces(&patch);
        check(&e);
        assert!(!e.is_closed());
        assert!(e.non_manifold().is_empty());
        let rims = e.rims();
        assert_eq!(rims.len(), 1);
        assert_eq!(rims[0].edges.len(), 8);
        // Every coordinate is a multiple of a half: no rounding anywhere.
        assert_eq!(rims[0].length, 8.0 * h);
        assert_eq!(rims[0].span, 2.0_f64.sqrt());
        // The rim runs with the faces left on its left: round the hole
        // clockwise seen from outside, above it.
        let [a, b] = e.edge(rims[0].edges[0]).unwrap();
        let [_, c] = e.edge(rims[0].edges[1]).unwrap();
        let [a, b, c] = [a, b, c].map(|v| e.position(v).unwrap());
        assert!((b - a).cross(c - b).z <= 0.0);
        let pieces = e.components();
        assert_eq!(pieces.len(), 1);
        assert_eq!(pieces[0].area, 6.0 * side * side - 4.0 * h * h);
        // The patch's middle vertex is gone with its faces.
        let back = e.to_triangulation();
        assert_eq!(back.positions.len(), 6 * 36 + 2 - 1);
        assert_eq!(back.triangles.len(), cube.triangles.len() - 8);
    }

    #[test]
    fn a_bowtie_is_split_into_a_vertex_per_fan() {
        let positions = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(-1.0, -1.0, 0.0),
        ];
        let m = mesh(positions.clone(), vec![[0, 1, 2], [0, 3, 4]]);
        let e = Editable::from(&m, 0.0).unwrap();
        check(&e);
        assert_eq!(e.vertex_count(), 6);
        assert_eq!(e.non_manifold(), vec![0, 5]);
        assert_eq!(e.source_vertex(5), Some(0));
        assert_eq!(e.components().len(), 2);
        let rims = e.rims();
        assert_eq!(rims.len(), 2);
        assert_eq!(rims[0].length, 1.0 + SQRT_2 + 1.0);
        let back = e.to_triangulation();
        assert_eq!(back.triangles, vec![[0, 1, 2], [5, 3, 4]]);
        assert_eq!(back.positions[5], positions[0]);
        assert_eq!(&back.positions[..5], &positions[..]);
    }

    #[test]
    fn removing_faces_splits_the_pinch_they_leave() {
        // A two by two grid less two opposite cells meets itself only at
        // the middle vertex.
        let m = grid(2);
        let mut e = Editable::from(&m, 0.0).unwrap();
        e.remove_faces(&[0, 1, 6, 7]);
        check(&e);
        assert_eq!(e.non_manifold(), vec![4, 9]);
        assert_eq!(e.source_vertex(9), Some(4));
        assert_eq!(e.components().len(), 2);
        let rims = e.rims();
        assert_eq!(rims.len(), 2);
        assert!(rims.iter().all(|r| r.length == 4.0));
        // Putting a cell back is not possible, but flipping both cells
        // over leaves the pinch as it was.
        e.flip_faces(&[2, 3, 4, 5]);
        check(&e);
        assert_eq!(e.non_manifold(), vec![4, 9]);
    }

    #[test]
    fn an_edge_of_three_faces_joins_the_sheet_and_leaves_the_fin() {
        // A sheet across the `x` axis and a fin standing on it.
        let m = mesh(
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
                Point::new(0.5, 1.0, 0.0),
                Point::new(0.5, -1.0, 0.0),
                Point::new(0.5, 0.0, 1.0),
            ],
            vec![[0, 1, 2], [1, 0, 3], [0, 1, 4]],
        );
        let e = Editable::from(&m, 0.0).unwrap();
        check(&e);
        assert_eq!(e.twin[0], 3);
        assert_eq!(e.face(2), Some([5, 6, 4]));
        assert_eq!(e.non_manifold(), vec![0, 1, 5, 6]);
        assert_eq!(e.source_vertex(5), Some(0));
        assert_eq!(e.source_vertex(6), Some(1));
        let pieces = e.components();
        assert_eq!(pieces.len(), 2);
        assert_eq!(pieces[1].faces, vec![2]);
    }

    #[test]
    fn two_bodies_sharing_an_edge_each_close() {
        let points = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.5, 1.0, 0.0),
            Point::new(0.5, 0.5, 1.0),
            Point::new(0.5, -1.0, 0.0),
            Point::new(0.5, -0.5, -1.0),
        ];
        let first = tetrahedron([points[0], points[1], points[2], points[3]], 0);
        let mut second = tetrahedron([points[0], points[1], points[4], points[5]], 0);
        for t in &mut second {
            *t = t.map(|v| if v >= 2 { v + 2 } else { v });
        }
        let m = mesh(points.to_vec(), [first, second].concat());
        let e = Editable::from(&m, 0.0).unwrap();
        check(&e);
        assert!(e.is_closed());
        assert_eq!(e.components().len(), 2);
        assert_eq!(e.non_manifold(), vec![0, 1, 6, 7]);
    }

    #[test]
    fn flipping_faces_reopens_and_recloses_edges() {
        let closed = torus(10, 6);
        let mut e = Editable::from(&closed, 0.0).unwrap();
        let volume = closed.volume();

        let all: Vec<u32> = (0..e.face_count()).map(index).collect();
        e.flip_faces(&all);
        check(&e);
        assert!(e.is_closed());
        assert!(e.non_manifold().is_empty());
        let flipped = e.to_triangulation();
        // The same tetrahedra summed with their signs turned.
        assert_eq!(flipped.volume(), -volume);
        e.flip_faces(&all);
        reads_back(&e, &closed);

        // One face turned against its neighbours parts from all three.
        e.flip_faces(&[17, 17]);
        check(&e);
        assert!(!e.is_closed());
        let rims = e.rims();
        assert_eq!(rims.len(), 2);
        // The same three edges, summed from different starts.
        assert!((rims[0].length - rims[1].length).abs() <= 4.0 * f64::EPSILON * rims[0].length);
        assert_eq!(e.components().len(), 2);
        assert!(!e.non_manifold().is_empty());
        assert!(
            (0..e.vertex_count())
                .filter_map(|v| e.position(index(v)).map(|p| (v, p)))
                .all(|(v, p)| closed.positions[e.source_vertex(index(v)).unwrap() as usize] == p)
        );
        e.flip_faces(&[17]);
        check(&e);
        assert!(e.is_closed());
        assert!(e.non_manifold().is_empty());
        reads_back(&e, &closed);
    }

    #[test]
    fn sources_survive_welding_removal_and_flips() {
        let pieces = soup(&grid(4));
        let mut e = Editable::from(&pieces, 1e-9).unwrap();
        check(&e);
        assert_eq!(e.vertex_count(), 25);
        // Each vertex is the first copy of its position the soup lists.
        for v in 0..25 {
            let s = e.source_vertex(v).unwrap() as usize;
            let p = e.position(v).unwrap();
            assert_eq!(pieces.positions[s], p);
            assert_eq!(pieces.positions.iter().position(|&q| q == p), Some(s));
        }
        e.remove_faces(&[0, 1, 10, 11, 12, 13, 22, 23]);
        e.flip_faces(&[2, 3, 30, 31]);
        check(&e);
        assert!(e.vertex_count() > 25);
        for v in 0..e.vertex_count() {
            if let Some(p) = e.position(index(v)) {
                let s = e.source_vertex(index(v)).unwrap();
                assert_eq!(pieces.positions[s as usize], p);
            }
        }
        assert_eq!(e.source_vertex(index(e.vertex_count())), None);
    }

    #[test]
    fn single_precision_copies_weld_within_their_quantum() {
        // A cube read as a soup of `f32` triangles, each copy of a corner
        // after the first moved by one unit in the last place.
        let cube = cube_of_grids(1, 100.0);
        let mut pieces = soup(&cube);
        let mut seen: Vec<Point> = Vec::new();
        for p in &mut pieces.positions {
            let single = |x: f64| f64::from(x as f32);
            let up = |x: f64| f64::from(f32::from_bits((x as f32).to_bits() + 1));
            let q = Point::new(single(p.x), single(p.y), single(p.z));
            *p = if seen.contains(&q) {
                Point::new(up(q.x), q.y, up(q.z))
            } else {
                seen.push(q);
                q
            };
        }
        // Units in the last place at 100 are 2^-17.
        let quantum = 2.0_f64.powi(-17);
        let e = Editable::from(&pieces, 2.0 * quantum).unwrap();
        check(&e);
        assert!(e.is_closed());
        assert_eq!(e.to_triangulation().positions.len(), 8);
        for v in 0..8 {
            let s = e.source_vertex(v).unwrap() as usize;
            assert_eq!(e.position(v), Some(seen[v as usize]));
            assert_eq!(pieces.positions[s], seen[v as usize]);
        }
        // Closer than the copies stand, they stay apart.
        let apart = Editable::from(&pieces, quantum / 2.0).unwrap();
        assert!(!apart.is_closed());
    }

    #[test]
    fn a_weld_tie_goes_to_the_lowest_kept_point() {
        // The third point is within the distance of both kept ones, and
        // the grid holds the later one in the first cell it looks at.
        let points = [1.5, 0.0, 0.75].map(|x| Point::new(x, 0.0, 0.0));
        let (kept, remap, founders) = weld_points(&points, 1.0);
        assert_eq!(kept.len(), 2);
        assert_eq!(remap, vec![0, 1, 0]);
        assert_eq!(founders, vec![0, 1]);
    }

    #[test]
    fn bad_input_is_refused() {
        let refusal = |m: &Triangulation, weld: f64| Editable::from(m, weld).unwrap_err();
        let m = mesh(vec![Point::new(0.0, 0.0, 0.0); 3], vec![[0, 1, 3]]);
        let e = refusal(&m, 0.0);
        assert!(matches!(e, ogeom_core::OgeomError::Range(_)));
        assert!(e.to_string().contains("names vertex 3 of 3"), "{e}");
        let m = grid(1);
        for weld in [-1.0, f64::NAN] {
            let e = refusal(&m, weld);
            assert!(matches!(e, ogeom_core::OgeomError::Numeric(_)));
            assert!(e.to_string().contains("weld distance"), "{e}");
        }
        let mut m = grid(1);
        m.positions[2].x = f64::NAN;
        let e = refusal(&m, 0.0);
        assert!(e.to_string().contains("vertex 2 is not finite"), "{e}");
    }

    #[test]
    fn an_editable_mesh_crosses_threads() {
        fn send_sync<T: Send + Sync>() {}
        send_sync::<Editable>();
    }

    #[test]
    fn the_span_is_the_farthest_pair() {
        let points: Vec<Point> = (0..200)
            .map(|k| {
                let t = f64::from(k) * 0.37;
                Point::new(t.sin() * 3.0, (t * 1.7).cos() * 2.0, (t * 0.3).sin())
            })
            .collect();
        let brute = points
            .iter()
            .flat_map(|p| points.iter().map(move |q| p.distance(*q)))
            .fold(0.0, f64::max);
        assert_eq!(diameter(&points), brute);
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]

        /// After any removals and flips the structure is sound, and reads
        /// the same as one built afresh from what it writes.
        #[test]
        fn edits_keep_a_sound_manifold(
            edits in prop::collection::vec(
                (any::<bool>(), prop::collection::vec(0_u32..60, 0..8)),
                1..6,
            )
        ) {
            let mut e = Editable::from(&torus(6, 5), 0.0).unwrap();
            for (remove, faces) in &edits {
                if *remove {
                    e.remove_faces(faces);
                } else {
                    e.flip_faces(faces);
                }
                check(&e);
            }
            let written = e.to_triangulation();
            let again = Editable::from(&written, 0.0).unwrap();
            prop_assert_eq!(again.to_triangulation(), written);
            prop_assert_eq!(again.is_closed(), e.is_closed());
            prop_assert_eq!(again.components().len(), e.components().len());
            let lengths = |m: &Editable| {
                let mut l: Vec<usize> = m.rims().iter().map(|r| r.edges.len()).collect();
                l.sort_unstable();
                l
            };
            prop_assert_eq!(lengths(&again), lengths(&e));
        }
    }
}
