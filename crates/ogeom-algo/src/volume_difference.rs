//! The volume between a scanned mesh and the solid it was made to: the
//! material the solid has and the scan lacks, and the material the scan has
//! beyond the solid, as connected regions.
//!
//! # The grid
//!
//! Space over both bodies' boxes is cut into cubic cells `cell` wide, and
//! each cell is read at its centre, its node. A node is in the scan by the
//! parity of the scan's crossings along the line through it parallel to
//! `x` ([`parity`](ogeom_mesh::parity)), exactly, and on the scan counts as
//! in it. A node is in the solid by the same parity against the solid
//! meshed at a sixteenth of a cell. Where the node lies within that chord,
//! plus the confusion distance, of the solid's mesh, the band where the
//! mesh and the solid can disagree, it is read again against the solid's
//! own surfaces by [`SolidBoundary::holds`], and on the solid counts as in
//! it. A node in the solid and out of the scan is missing material; one in
//! the scan and out of the solid is excess material.
//!
//! # Holes in the scan
//!
//! Vertices of the scan at one position are made one first, so a mesh
//! stored as loose triangles reads as the surface they make. The parity
//! needs every edge used an even number of times; a scan with an edge
//! used an odd number of times, a hole's rim, is closed first, each rim
//! filled flat with [`Editable::fill`] on its own vertices. A differing
//! node whose nearest point on the closed scan lies on a fill was not seen
//! by the scanner: it takes the solid's side, so a hole reports no
//! material either way.
//!
//! Such a node is *doubtful*, and so is every node within a cell of a
//! fill, where the fill may stand where the scan would have, and a node
//! the solid's surfaces could not classify, every ray meeting a
//! degeneracy, which keeps the reading of the solid's mesh.
//!
//! # Regions
//!
//! Nodes that differ the same way are joined across faces, edges and
//! corners into regions. A region's volume is its nodes' count times a
//! cell's volume. A node is counted where its centre lies, so the volume
//! is off by at most the cells the region's boundary passes through: about
//! the boundary's area times `√3` cells, and less where the boundary's
//! errors cancel, as they do on a surface not aligned to the grid.
//!
//! A region's depth is how far apart the two surfaces stand across it: the
//! greatest distance, from a vertex of either surface within a cell of the
//! region's cells, to the other surface. A vertex of the scan is measured
//! to the solid's mesh, within its chord of the solid, and a vertex of the
//! solid's mesh to the closed scan. A boss standing off a face reads its
//! height, from its top; a pocket its depth, from its floor. A region no
//! vertex comes near reads the greatest distance from its nodes to either
//! surface.
//!
//! A region shallower than a cell is below the grid's resolution, which
//! puts at least one node into a slab a cell thick, wherever the grid
//! falls and however the slab is turned, but not into a thinner one: it is
//! not reported. A noisy scan
//! lifts single nodes off its surface all over, and each is that shallow.
//!
//! # Layers and lumps
//!
//! A region `d` deep and `V` in volume is as wide as a square slab of that
//! thickness, `w = √(V / d)`. It is a [`Form::Layer`], a face standing
//! off, where the slab's two faces have more area than its four sides,
//! `2 w² ≥ 4 w d`, so where `V ≥ 4 d³`; otherwise a [`Form::Lump`].
//!
//! # Unsure
//!
//! A region is [`Region::unsure`] where a doubtful node lies within a cell
//! of one of its nodes, or where it is shallower than two cells. Then its
//! extent rests on nodes nobody measured, or its volume error, about its
//! area times `√3` cells against its area times its depth, can reach its
//! whole volume.
//!
//! # Cost
//!
//! Each row of nodes along `x` is read in parallel, at a parity line
//! through each body and an exact classification of the nodes in the
//! solid's band, an eighth of those within a cell of its surface.
//! Rows are read and regions numbered in the grid's order, so the answer
//! is the same bit for bit on any thread count.

use ogeom_core::{FastMap, OgeomError, OgeomResult, Tolerances, ogeom_bail, parallel};
use ogeom_math::{Aabb, Point};
use ogeom_mesh::{Bvh, Deflection, Editable, Fairing, Fill, Line, triangulate};
use ogeom_topo::{Model, Shape, Triangulation};

use crate::classify::{SolidBoundary, exact_ring_chord};

/// Which way a region differs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Material {
    /// In the solid and out of the scan: a boss the scan lacks, or a
    /// pocket it has.
    Missing,
    /// In the scan and out of the solid: a boss the solid lacks, or a
    /// pocket the scan fills.
    Excess,
}

/// The shape of a region.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Form {
    /// Wider than twice its depth: a face standing off.
    Layer,
    /// Less than twice as wide as deep: a boss or a pocket.
    Lump,
}

/// A connected volume where the scan and the solid differ the same way.
#[derive(Debug, Clone, PartialEq)]
pub struct Region {
    /// Which way they differ.
    pub side: Material,
    /// Its nodes' count times a cell's volume.
    pub volume: f64,
    /// How far apart the two surfaces stand across it; see the
    /// [module documentation](self).
    pub depth: f64,
    /// A layer or a lump.
    pub kind: Form,
    /// The box of its cells.
    pub bounds: Aabb,
    /// Whether it rests on a hole of the scan, or is too shallow for its
    /// volume to be read; see the [module documentation](self).
    pub unsure: bool,
}

/// How many vertices a worker measures at a time.
const VERTEX_BATCH: usize = 4096;

/// The regions where `mesh`, a scan, and the solid `shape` differ, on a
/// grid of cells `cell` wide; see the [module documentation](self).
///
/// Regions come largest first, equal volumes in the grid's order.
///
/// # Errors
///
/// [`OgeomError::Construction`] if `cell` is not a positive distance, if
/// `shape` is not a solid or a closed shell, or if a hole of the scan
/// cannot be filled, as [`Editable::fill`] refuses.
/// [`OgeomError::NullObject`] if the scan has no triangles.
/// [`OgeomError::NotDone`] if the solid's mesh misses its chord.
/// [`OgeomError::Range`] if the grid would have more than `u32::MAX`
/// nodes along an axis. As [`Bvh::build`] for a scan with a vertex that is
/// not finite.
pub fn volume_difference(
    mesh: &Triangulation,
    model: &Model,
    shape: &Shape,
    cell: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<Region>> {
    if !(cell.is_finite() && cell > 0.0) {
        ogeom_bail!(Construction, "a cell of {cell} is not a distance");
    }
    if mesh.triangles.is_empty() {
        ogeom_bail!(NullObject, "the scan has no triangles");
    }
    let closed = Closed::of(mesh)?;
    let scan = Bvh::build(&closed.mesh)?;

    let chord = cell / 16.0;
    let solid_mesh = triangulate(
        model,
        shape,
        Deflection {
            chord,
            ..Deflection::default()
        },
        tol,
    )?;
    if solid_mesh.is_empty() || !solid_mesh.is_closed() {
        ogeom_bail!(
            Construction,
            "the solid's boundary is not closed, so there is no inside to be in"
        );
    }
    if !solid_mesh.deflection_met {
        ogeom_bail!(NotDone, "the solid's mesh misses its chord of {chord}");
    }
    let solid = Bvh::build(&solid_mesh)?;
    let boundary = SolidBoundary::of(model, shape, exact_ring_chord(tol), tol)?;
    let reader = Reader {
        grid: Grid::over(&scan.bounds().union(&solid.bounds()), cell)?,
        scan: &scan,
        first_fill: closed.first_fill,
        solid: &solid,
        boundary: &boundary,
        model,
        band: tol.confusion() + chord,
        tol,
    };

    let mut read = reader.read()?;
    read.doubtful.extend(closed.near_fills(&reader.grid)?);
    read.doubtful.sort_unstable();
    read.doubtful.dedup();
    let (grid, differ) = (&reader.grid, &read.differ);
    let mut regions = Regions::join(grid, differ);
    regions.near_doubt(grid, differ, &read.doubtful);
    let mut depth = vec![f64::NEG_INFINITY; regions.count];
    regions.deepen(grid, differ, &mesh.positions, &solid, &mut depth);
    regions.deepen(grid, differ, &solid_mesh.positions, &scan, &mut depth);
    regions.deepen_by_nodes(&reader, differ, &mut depth);
    Ok(regions.finish(grid, differ, &depth))
}

/// The scan, closed, and where the fills start among its triangles.
struct Closed {
    mesh: Triangulation,
    /// The first triangle of the fills; every triangle where the scan was
    /// closed.
    first_fill: usize,
}

/// The mesh with its vertices at one position made one, and whether every
/// edge joining two of them is then used an even number of times, as the
/// parity asks.
fn welded(mesh: &Triangulation) -> (Triangulation, bool) {
    let key = |p: Point| [p.x, p.y, p.z].map(|c| (c + 0.0).to_bits());
    let mut first: FastMap<[u64; 3], u32> = FastMap::default();
    let mut out = Triangulation::new();
    let mut to = Vec::with_capacity(mesh.positions.len());
    for &p in &mesh.positions {
        let next = index(out.positions.len());
        let at = *first.entry(key(p)).or_insert(next);
        if at == next {
            out.positions.push(p);
        }
        to.push(at);
    }
    let mut odd: FastMap<(u32, u32), bool> = FastMap::default();
    for t in &mesh.triangles {
        let t = t.map(|v| to.get(v as usize).copied().unwrap_or(v));
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            if a != b {
                let flag = odd.entry((a.min(b), a.max(b))).or_insert(false);
                *flag = !*flag;
            }
        }
        out.triangles.push(t);
    }
    out.deflection_met = mesh.deflection_met;
    let even = odd.values().all(|&o| !o);
    (out, even)
}

impl Closed {
    fn of(mesh: &Triangulation) -> OgeomResult<Self> {
        let (mesh, even) = welded(mesh);
        if even {
            return Ok(Self {
                first_fill: mesh.triangles.len(),
                mesh,
            });
        }
        let mut editable = Editable::from(&mesh, 0.0)?;
        let fill = Fill {
            fairing: Fairing::Flat,
            refine: false,
        };
        let mut made = 0;
        // A fill that joins the vertices where two rims touch can rename
        // the other rim; it is filled on the next pass, under its new name.
        loop {
            let rims = editable.rims();
            if rims.is_empty() {
                break;
            }
            let mut filled = false;
            for rim in rims {
                match editable.fill(rim.id(), &fill) {
                    Ok(report) => {
                        made += report.faces.len();
                        filled = true;
                    }
                    Err(OgeomError::Dangling(_)) => {}
                    Err(e) => ogeom_bail!(
                        Construction,
                        "a hole of the scan, rim {}, cannot be filled: {e}",
                        rim.id().0
                    ),
                }
            }
            if !filled {
                ogeom_bail!(Construction, "the scan's holes cannot be filled");
            }
        }
        let (closed, even) = welded(&editable.to_triangulation());
        if !even {
            ogeom_bail!(
                Construction,
                "the scan is not closed once its holes are filled: an edge is used \
                 an odd number of times"
            );
        }
        Ok(Self {
            first_fill: closed.triangles.len() - made,
            mesh: closed,
        })
    }

    /// The nodes within a cell of a fill, in the grid's order.
    fn near_fills(&self, grid: &Grid) -> OgeomResult<Vec<u64>> {
        let mut fills = Triangulation::new();
        fills.positions.clone_from(&self.mesh.positions);
        fills.triangles = self.mesh.triangles[self.first_fill..].to_vec();
        if fills.triangles.is_empty() {
            return Ok(Vec::new());
        }
        let bvh = Bvh::build(&fills)?;
        let h = grid.cell;
        let found = parallel::map_ordered(&fills.triangles, |_, t| {
            let corners = t.map(|v| fills.positions[v as usize]);
            let reach = Aabb::of_points(&corners).expanded(h);
            let (Some(low), Some(high)) = (reach.low(), reach.high()) else {
                return Vec::new();
            };
            let (from, to) = (grid.cell_of(low), grid.cell_of(high));
            let mut out = Vec::new();
            for k in from[2]..=to[2] {
                for j in from[1]..=to[1] {
                    for i in from[0]..=to[0] {
                        let Some(n) = grid.on([i, j, k]) else {
                            continue;
                        };
                        if bvh.nearest(grid.at(n), h).is_some() {
                            out.push(grid.key(n));
                        }
                    }
                }
            }
            out
        });
        let mut near: Vec<u64> = found.into_iter().flatten().collect();
        near.sort_unstable();
        near.dedup();
        Ok(near)
    }
}

/// The nodes: cell centres, `origin + (i + ½) cell` along each axis.
struct Grid {
    origin: Point,
    cell: f64,
    count: [u32; 3],
}

impl Grid {
    fn over(bounds: &Aabb, cell: f64) -> OgeomResult<Self> {
        let (Some(low), Some(high)) = (bounds.low(), bounds.high()) else {
            ogeom_bail!(NullObject, "neither body has a triangle");
        };
        let size = high - low;
        let along = |extent: f64| -> OgeomResult<u32> {
            let count = (extent / cell).ceil().max(1.0);
            if !(count.is_finite() && count <= f64::from(u32::MAX)) {
                ogeom_bail!(
                    Range,
                    "a grid of {count} nodes along an axis, a cell of {cell} across {extent}"
                );
            }
            #[allow(
                clippy::cast_possible_truncation,
                clippy::cast_sign_loss,
                reason = "a whole number in 1..=u32::MAX, checked above"
            )]
            Ok(count as u32)
        };
        Ok(Self {
            origin: low,
            cell,
            count: [along(size.x)?, along(size.y)?, along(size.z)?],
        })
    }

    fn coordinate(&self, axis: usize, i: u32) -> f64 {
        let o = [self.origin.x, self.origin.y, self.origin.z][axis];
        (f64::from(i) + 0.5).mul_add(self.cell, o)
    }

    fn at(&self, [i, j, k]: [u32; 3]) -> Point {
        Point::new(
            self.coordinate(0, i),
            self.coordinate(1, j),
            self.coordinate(2, k),
        )
    }

    fn key(&self, [i, j, k]: [u32; 3]) -> u64 {
        (u64::from(k) * u64::from(self.count[1]) + u64::from(j)) * u64::from(self.count[0])
            + u64::from(i)
    }

    fn node(&self, key: u64) -> [u32; 3] {
        let nx = u64::from(self.count[0]);
        let ny = u64::from(self.count[1]);
        #[allow(
            clippy::cast_possible_truncation,
            reason = "each index is below its axis's u32 count"
        )]
        [
            (key % nx) as u32,
            ((key / nx) % ny) as u32,
            (key / nx / ny) as u32,
        ]
    }

    /// The node whose cell holds `p` along each axis, as a signed index
    /// that may lie off the grid.
    fn cell_of(&self, p: Point) -> [i64; 3] {
        let o = [self.origin.x, self.origin.y, self.origin.z];
        let q = [p.x, p.y, p.z];
        core::array::from_fn(|a| {
            let at = ((q[a] - o[a]) / self.cell).floor();
            #[allow(
                clippy::cast_possible_truncation,
                reason = "clamped into i64 range first; a node off the grid is dropped"
            )]
            let at = at.clamp(-2.0, f64::from(u32::MAX) + 1.0) as i64;
            at
        })
    }

    /// The node at signed indices `n`, where it lies on the grid.
    fn on(&self, n: [i64; 3]) -> Option<[u32; 3]> {
        let on = (0..3).all(|a| n[a] >= 0 && n[a] < i64::from(self.count[a]));
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "checked to lie in 0..count above"
        )]
        on.then(|| n.map(|v| v as u32))
    }

    /// The nodes within one step of `at` along every axis, `at` among them,
    /// that lie on the grid.
    fn around(&self, at: [i64; 3]) -> impl Iterator<Item = [u32; 3]> + '_ {
        (-1..=1_i64).flat_map(move |dk| {
            (-1..=1_i64).flat_map(move |dj| {
                (-1..=1_i64).filter_map(move |di| self.on([at[0] + di, at[1] + dj, at[2] + dk]))
            })
        })
    }
}

/// What reading the grid needs.
struct Reader<'a> {
    grid: Grid,
    scan: &'a Bvh,
    first_fill: usize,
    solid: &'a Bvh,
    boundary: &'a SolidBoundary,
    model: &'a Model,
    /// How near the solid's mesh a node is read against its surfaces.
    band: f64,
    tol: Tolerances,
}

/// The grid as read: the nodes that differ, in the grid's order, and the
/// doubtful ones, in no order.
struct Read {
    differ: Vec<(u64, Material)>,
    doubtful: Vec<u64>,
}

/// One row's reading.
#[derive(Default)]
struct Row {
    differ: Vec<(u32, Material)>,
    doubtful: Vec<u32>,
}

impl Reader<'_> {
    fn read(&self) -> OgeomResult<Read> {
        let [_, ny, nz] = self.grid.count;
        let rows: Vec<u64> = (0..u64::from(ny) * u64::from(nz)).collect();
        let read = parallel::map_ordered(&rows, |_, &r| {
            #[allow(
                clippy::cast_possible_truncation,
                reason = "each index is below its axis's u32 count"
            )]
            let (j, k) = ((r % u64::from(ny)) as u32, (r / u64::from(ny)) as u32);
            self.row(j, k)
        });
        let nx = u64::from(self.grid.count[0]);
        let mut differ = Vec::new();
        let mut doubtful = Vec::new();
        for (r, row) in rows.iter().zip(read) {
            let row = row?;
            differ.extend(row.differ.iter().map(|&(i, m)| (r * nx + u64::from(i), m)));
            doubtful.extend(row.doubtful.iter().map(|&i| r * nx + u64::from(i)));
        }
        self.unseen(&mut differ, &mut doubtful);
        Ok(Read { differ, doubtful })
    }

    fn row(&self, j: u32, k: u32) -> OgeomResult<Row> {
        let (y, z) = (self.grid.coordinate(1, j), self.grid.coordinate(2, k));
        let nx = self.grid.count[0];
        let scan = Line::through(self.scan, y, z);
        let solid = Line::through(self.solid, y, z);
        let near = self.near_band(y, z);
        let mut row = Row::default();
        for i in 0..nx {
            let x = self.grid.coordinate(0, i);
            let in_scan = scan.lies(x).is_inside_or_on();
            let mut in_solid = solid.lies(x).is_inside_or_on();
            if near[i as usize] {
                let p = Point::new(x, y, z);
                if self.solid.nearest(p, self.band).is_some() {
                    match self.boundary.holds(self.model, p, self.tol) {
                        Ok(c) => in_solid = c.is_inside_or_on(),
                        Err(OgeomError::NotDone(_)) => row.doubtful.push(i),
                        Err(e) => return Err(e),
                    }
                }
            }
            match (in_solid, in_scan) {
                (true, false) => row.differ.push((i, Material::Missing)),
                (false, true) => row.differ.push((i, Material::Excess)),
                _ => {}
            }
        }
        Ok(row)
    }

    /// Which nodes of the row through `(y, z)` lie within the band of a
    /// triangle of the solid's mesh along `x`, a box about it: the only
    /// ones that can be within the band of the mesh.
    fn near_band(&self, y: f64, z: f64) -> Vec<bool> {
        let band = self.band;
        let nx = self.grid.count[0];
        let mut near = vec![false; nx as usize];
        for t in self.solid.triangles_where(|b| match (b.low(), b.high()) {
            (Some(l), Some(h)) => {
                l.y - band <= y && y <= h.y + band && l.z - band <= z && z <= h.z + band
            }
            _ => false,
        }) {
            let Some(c) = self.solid.triangle(t) else {
                continue;
            };
            let low = c[0].x.min(c[1].x).min(c[2].x) - band;
            let high = c[0].x.max(c[1].x).max(c[2].x) + band;
            let first = ((low - self.grid.origin.x) / self.grid.cell - 0.5).ceil();
            let last = ((high - self.grid.origin.x) / self.grid.cell - 0.5).floor();
            let first = first.max(0.0);
            let last = last.min(f64::from(nx) - 1.0);
            if first > last {
                continue;
            }
            #[allow(
                clippy::cast_possible_truncation,
                clippy::cast_sign_loss,
                reason = "clamped into 0..nx above"
            )]
            for slot in &mut near[first as usize..=last as usize] {
                *slot = true;
            }
        }
        near
    }

    /// Give the solid's side to each differing node whose nearest point on
    /// the closed scan lies on a fill, and doubt it.
    fn unseen(&self, differ: &mut Vec<(u64, Material)>, doubtful: &mut Vec<u64>) {
        if self.first_fill == self.scan.triangle_count() || differ.is_empty() {
            return;
        }
        let points: Vec<Point> = differ
            .iter()
            .map(|&(key, _)| self.grid.at(self.grid.node(key)))
            .collect();
        let hits = self.scan.nearest_many(&points, f64::INFINITY);
        let mut kept = Vec::with_capacity(differ.len());
        for (&(key, m), hit) in differ.iter().zip(hits) {
            if hit.is_some_and(|h| h.triangle as usize >= self.first_fill) {
                doubtful.push(key);
            } else {
                kept.push((key, m));
            }
        }
        *differ = kept;
    }
}

/// The differing nodes joined into regions.
struct Regions {
    /// Each differing node's region, in the order of the nodes.
    of: Vec<u32>,
    count: usize,
    /// Whether a doubtful node lies within a cell of each region.
    doubted: Vec<bool>,
}

/// The root of `x`'s set, halving the path to it.
fn find(parent: &mut [u32], mut x: u32) -> u32 {
    while parent[x as usize] != x {
        parent[x as usize] = parent[parent[x as usize] as usize];
        x = parent[x as usize];
    }
    x
}

fn index(i: usize) -> u32 {
    u32::try_from(i).unwrap_or(u32::MAX)
}

/// Where `key` stands among the sorted differing nodes.
fn place(differ: &[(u64, Material)], key: u64) -> Option<usize> {
    differ.binary_search_by_key(&key, |&(k, _)| k).ok()
}

impl Regions {
    fn join(grid: &Grid, differ: &[(u64, Material)]) -> Self {
        let mut parent: Vec<u32> = (0..differ.len()).map(index).collect();
        for (at, &(key, side)) in differ.iter().enumerate() {
            let node = grid.node(key);
            let centre = node.map(i64::from);
            for n in grid.around(centre) {
                let other = grid.key(n);
                if other <= key {
                    continue;
                }
                if let Some(o) = place(differ, other)
                    && differ[o].1 == side
                {
                    let (a, b) = (find(&mut parent, index(at)), find(&mut parent, index(o)));
                    if a != b {
                        parent[a.max(b) as usize] = a.min(b);
                    }
                }
            }
        }
        // Regions numbered by their first node, in the grid's order.
        let mut number = vec![u32::MAX; differ.len()];
        let mut of = Vec::with_capacity(differ.len());
        let mut count = 0;
        for at in 0..differ.len() {
            let root = find(&mut parent, index(at)) as usize;
            if number[root] == u32::MAX {
                number[root] = index(count);
                count += 1;
            }
            of.push(number[root]);
        }
        Self {
            of,
            count,
            doubted: vec![false; count],
        }
    }

    /// The regions of the differing nodes within one step of node `at`.
    fn near(&self, grid: &Grid, differ: &[(u64, Material)], at: [i64; 3], out: &mut Vec<u32>) {
        out.clear();
        for n in grid.around(at) {
            if let Some(o) = place(differ, grid.key(n)) {
                let r = self.of[o];
                if !out.contains(&r) {
                    out.push(r);
                }
            }
        }
    }

    /// Mark the regions within a cell of a doubtful node.
    fn near_doubt(&mut self, grid: &Grid, differ: &[(u64, Material)], doubtful: &[u64]) {
        let mut out = Vec::new();
        for &key in doubtful {
            self.near(grid, differ, grid.node(key).map(i64::from), &mut out);
            for &r in &out {
                self.doubted[r as usize] = true;
            }
        }
    }

    /// Raise each region's depth to the distance, to `other`, of each of
    /// `vertices` within a cell of its cells.
    fn deepen(
        &self,
        grid: &Grid,
        differ: &[(u64, Material)],
        vertices: &[Point],
        other: &Bvh,
        depth: &mut [f64],
    ) {
        if differ.is_empty() {
            return;
        }
        let starts: Vec<usize> = (0..vertices.len()).step_by(VERTEX_BATCH).collect();
        let found = parallel::map_ordered(&starts, |_, &start| {
            let mut out = Vec::new();
            let mut near = Vec::new();
            for &p in &vertices[start..(start + VERTEX_BATCH).min(vertices.len())] {
                self.near(grid, differ, grid.cell_of(p), &mut near);
                if near.is_empty() {
                    continue;
                }
                if let Some(hit) = other.nearest(p, f64::INFINITY) {
                    out.extend(near.iter().map(|&r| (r, hit.distance)));
                }
            }
            out
        });
        for (r, d) in found.into_iter().flatten() {
            depth[r as usize] = depth[r as usize].max(d);
        }
    }

    /// Give each region no vertex came near the greatest distance from its
    /// nodes to either body.
    fn deepen_by_nodes(&self, reader: &Reader<'_>, differ: &[(u64, Material)], depth: &mut [f64]) {
        let (asked, points): (Vec<u32>, Vec<Point>) = differ
            .iter()
            .zip(&self.of)
            .filter(|&(_, &r)| depth[r as usize] == f64::NEG_INFINITY)
            .map(|(&(key, _), &r)| (r, reader.grid.at(reader.grid.node(key))))
            .unzip();
        if points.is_empty() {
            return;
        }
        let to_scan = reader.scan.nearest_many(&points, f64::INFINITY);
        let to_solid = reader.solid.nearest_many(&points, f64::INFINITY);
        for ((r, a), b) in asked.into_iter().zip(to_scan).zip(to_solid) {
            let d = a
                .map_or(0.0, |h| h.distance)
                .max(b.map_or(0.0, |h| h.distance));
            depth[r as usize] = depth[r as usize].max(d);
        }
    }

    /// The regions as reported: deep enough to resolve, largest first.
    fn finish(&self, grid: &Grid, differ: &[(u64, Material)], depth: &[f64]) -> Vec<Region> {
        let mut low = vec![[u32::MAX; 3]; self.count];
        let mut high = vec![[0_u32; 3]; self.count];
        let mut cells = vec![0_u64; self.count];
        let mut side = vec![Material::Missing; self.count];
        for (&(key, m), &r) in differ.iter().zip(&self.of) {
            let r = r as usize;
            let n = grid.node(key);
            for a in 0..3 {
                low[r][a] = low[r][a].min(n[a]);
                high[r][a] = high[r][a].max(n[a]);
            }
            cells[r] += 1;
            side[r] = m;
        }
        let h = grid.cell;
        let mut regions: Vec<(usize, Region)> = (0..self.count)
            .filter(|&r| depth[r] >= h)
            .map(|r| {
                #[allow(
                    clippy::cast_precision_loss,
                    reason = "a count of nodes, exact below 2^53"
                )]
                let volume = cells[r] as f64 * h * h * h;
                let d = depth[r];
                let half = Point::new(h / 2.0, h / 2.0, h / 2.0);
                let (a, b) = (grid.at(low[r]), grid.at(high[r]));
                let region = Region {
                    side: side[r],
                    volume,
                    depth: d,
                    kind: if volume >= 4.0 * d * d * d {
                        Form::Layer
                    } else {
                        Form::Lump
                    },
                    bounds: Aabb::of_corners(
                        Point::new(a.x - half.x, a.y - half.y, a.z - half.z),
                        Point::new(b.x + half.x, b.y + half.y, b.z + half.z),
                    ),
                    unsure: self.doubted[r] || d < 2.0 * h,
                };
                (r, region)
            })
            .collect();
        regions.sort_by(|(ra, a), (rb, b)| b.volume.total_cmp(&a.volume).then(ra.cmp(rb)));
        regions.into_iter().map(|(_, region)| region).collect()
    }
}
