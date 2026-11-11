//! Decimating a triangulation.
//!
//! A mesh fine enough for a mass property is far finer than one needed to draw,
//! and a mesh fine enough to draw a whole assembly is far finer than one needed
//! for the bolt in the corner of it. Decimation is how one tessellation serves
//! both without being computed twice.
//!
//! # The error is measured, not hoped for
//!
//! Collapsing an edge moves the surface. *How far* it moves is what decides
//! whether the collapse is worth making, so every candidate carries the squared
//! distance from the merged vertex to the planes of every face that met there:
//! the quadric error metric of Garland and Heckbert. Summing plane distances
//! this way costs one small symmetric matrix per vertex and makes the choice a
//! comparison rather than a guess.
//!
//! The result reports the worst error it accepted. A decimation that returned
//! only a smaller mesh would be one nothing downstream could decide to trust.
//!
//! # What it will not touch
//!
//! A vertex on a boundary stays. The alternative is a constraint plane that
//! makes boundary collapses expensive but possible, and "expensive but
//! possible" means the outline of a sheet body creeps inward as the mesh
//! coarsens, which is exactly the thing a caller would not think to check.
//! Holding the boundary exactly is a stronger promise and a simpler one.
//!
//! A collapse that would turn a triangle inside out is refused for the same
//! reason: a fold is not a small error, it is a mesh that no longer bounds what
//! it did.
//!
//! # Verified
//!
//! The quadric bounds nothing on its own: it is a sum over planes, and the
//! nearest point of the source need not be on any of them. [`Target::Verified`]
//! makes a collapse only where the points it moves stay within the tolerance
//! of the other mesh and its triangles meet no others but at shared corners,
//! then measures the result against the source both ways, pins what fails
//! and decimates again; see [`Deviation`] for what it reports.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use ogeom_core::{FastMap, FastSet};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::Point;
use ogeom_topo::Triangulation;

mod checked;
mod verified;

/// How far to decimate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Target {
    /// Stop once the mesh is down to this many triangles, whatever the error.
    Triangles(usize),
    /// Collapse only while the error stays within this distance.
    ///
    /// The honest option: the caller says how wrong the mesh may be and gets
    /// however few triangles that allows, rather than naming a count and
    /// discovering the error afterwards.
    Error(f64),
    /// Collapse within `tolerance`, then measure the result against the
    /// source both ways and hold whatever the measurement refuses.
    ///
    /// Each source vertex is measured to the result, and each result vertex
    /// and triangle centre to the source. A collapse is made only where
    /// every point it moves stays within `tolerance`, and only where its
    /// triangles meet no other triangle but at the corners they share, so a
    /// mesh that bounded a solid still does. Where a point is farther than
    /// `tolerance` all the same, which rounding alone can cause, the source
    /// vertices behind it are pinned and the decimation is run again from
    /// the source, until at least `share` of the points on each side are
    /// within `tolerance`. The measurement is reported in
    /// [`Simplified::deviation`].
    ///
    /// With `max_triangles`, a result still larger than that is decimated
    /// further, pins released, and the deviation reported is the one
    /// measured on what came out, however far past `tolerance` it is.
    Verified {
        /// How far, as a distance, a point may be from the other mesh.
        tolerance: f64,
        /// The fraction of points on each side that must be within
        /// `tolerance`, in `(0, 1]`; 0.999 is a scan's usual demand.
        share: f64,
        /// At most this many triangles, whatever the deviation.
        max_triangles: Option<usize>,
    },
}

/// The distance between a decimated mesh and its source, measured both ways.
///
/// The points measured are the source's vertices, against the result, and
/// the result's vertices and triangle centres, against the source; each
/// distance is to the nearest point on the other mesh's triangles.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Deviation {
    /// The smaller of the two sides' fractions of points within the
    /// tolerance.
    pub within_share: f64,
    /// The distance within which the asked share of points lies, on the
    /// worse of the two sides.
    pub distance_at_share: f64,
    /// The point farthest from the other mesh, on either side, and its
    /// distance.
    pub furthest: (Point, f64),
    /// The root mean square distance over every point measured, both sides.
    pub rms: f64,
    /// How many source vertices the measurement pinned, counted before a
    /// `max_triangles` released them.
    pub pinned: usize,
    /// How many decimations were run, each from the source.
    pub rounds: usize,
}

/// What decimation produced.
#[derive(Debug, Clone)]
pub struct Simplified {
    /// The decimated mesh.
    pub mesh: Triangulation,
    /// The worst error accepted, as a distance.
    ///
    /// Zero when nothing was collapsed. Every vertex of the result is within
    /// this of the surface the original described.
    pub error: f64,
    /// How many edge collapses were made.
    pub collapsed: usize,
    /// Whether the target was reached.
    ///
    /// A mesh can run out of *valid* collapses before it runs out of triangles
    /// (every remaining edge is on a boundary or would fold something), and
    /// then the result is as small as it can safely be rather than as small as
    /// was asked for. Reported rather than passed off as success.
    pub target_met: bool,
    /// With [`Target::Verified`], what was measured; `None` otherwise.
    ///
    /// [`Simplified::error`] is then [`Deviation::furthest`]'s distance.
    pub deviation: Option<Deviation>,
    /// For each vertex of the source, the vertex of [`Simplified::mesh`]
    /// it was merged into or kept as; `None` for one no triangle of the
    /// result uses, as a vertex no triangle of the source used. The source
    /// vertices of one result vertex are those its collapses merged.
    pub result_vertex: Vec<Option<u32>>,
}

/// Decimate a mesh.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the target is
/// not a positive count or a positive distance, a share is not in `(0, 1]`, or
/// the mesh names a vertex it does not have;
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) with
/// [`Target::Verified`] if a vertex of a triangle is not finite;
/// [`OgeomError::Range`](ogeom_core::OgeomError::Range) if the mesh has more
/// than `u32::MAX` triangles.
pub fn simplify(mesh: &Triangulation, target: Target, tol: Tolerances) -> OgeomResult<Simplified> {
    match target {
        Target::Verified { tolerance, .. } if !tolerance.is_finite() || tolerance <= 0.0 => {
            ogeom_bail!(Construction, "a tolerance of {tolerance} is not a distance");
        }
        Target::Verified { share, .. } if !(share > 0.0 && share <= 1.0) => {
            ogeom_bail!(
                Construction,
                "a share of {share} is not a fraction of the points"
            );
        }
        Target::Verified {
            max_triangles: Some(0),
            ..
        } => {
            ogeom_bail!(Construction, "a mesh of no triangles describes nothing");
        }
        Target::Triangles(0) => {
            ogeom_bail!(Construction, "a mesh of no triangles describes nothing");
        }
        Target::Error(e) if !e.is_finite() || e <= 0.0 => {
            ogeom_bail!(Construction, "an error budget of {e} is not a distance");
        }
        _ => {}
    }
    if u32::try_from(mesh.triangles.len()).is_err() {
        ogeom_bail!(Range, "a mesh of more than four billion triangles");
    }
    for triangle in &mesh.triangles {
        for index in triangle {
            if *index as usize >= mesh.positions.len() {
                ogeom_bail!(
                    Construction,
                    "a triangle names vertex {index}, and the mesh has {}",
                    mesh.positions.len()
                );
            }
        }
    }

    let (budget, floor) = match target {
        Target::Error(e) => (e * e, 1),
        Target::Triangles(n) => (f64::MAX, n),
        Target::Verified {
            tolerance,
            share,
            max_triangles,
        } => return verified::simplify(mesh, tolerance, share, max_triangles, tol),
    };
    let mut state = State::new(mesh, tol);

    // Within reach the queue makes the collapses the scan makes, sooner.
    let (worst, collapsed) = if mesh.positions.iter().all(|p| within_reach(*p)) {
        state.decimate_by_queue(budget, floor)
    } else {
        state.decimate_by_scan(budget, floor)
    };

    let target_met = match target {
        Target::Triangles(n) => state.live_triangles() <= n,
        // Every collapse made was inside the budget, and the loop stops only
        // when no remaining one is.
        Target::Error(_) | Target::Verified { .. } => true,
    };
    let (mesh, renumbered) = state.harvest();
    Ok(Simplified {
        mesh,
        error: worst.max(0.0).sqrt(),
        collapsed,
        target_met,
        deviation: None,
        result_vertex: state.result_vertex(&renumbered),
    })
}

/// Whether every coordinate of `p` is at most 2^100 in size. A quadric
/// summed over at most 2^32 triangles of such points, weighed at such a
/// point, stays below 2^645 term by term, so no cost overflows to infinity
/// or turns NaN, and costs compare as numbers.
fn within_reach(p: Point) -> bool {
    const REACH: f64 = 1_267_650_600_228_229_401_496_703_205_376.0; // 2^100
    [p.x, p.y, p.z].iter().all(|c| c.abs() <= REACH)
}

/// A collapse [`State::cheapest`] would make, as the queue holds it:
/// cheapest first, ties to the edge the scan of live triangles meets
/// first.
#[derive(Debug, Clone, Copy)]
struct Candidate {
    cost: f64,
    /// Three times the index of the first live triangle holding the edge,
    /// plus the edge's slot in it.
    place: u64,
    /// The edge's ends, as that triangle walks it.
    from: u32,
    to: u32,
}

impl Candidate {
    /// Whether two candidates are the same collapse at the same cost.
    fn same(&self, other: &Self) -> bool {
        self.cost.to_bits() == other.cost.to_bits()
            && self.place == other.place
            && self.from == other.from
            && self.to == other.to
    }
}

impl PartialEq for Candidate {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == core::cmp::Ordering::Equal
    }
}

impl Eq for Candidate {}

impl PartialOrd for Candidate {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Candidate {
    /// By cost, as [`State::cheapest`] compares them (costs within reach
    /// are numbers, so this is a total order), then by place.
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.cost
            .partial_cmp(&other.cost)
            .unwrap_or(core::cmp::Ordering::Equal)
            .then(self.place.cmp(&other.place))
    }
}

/// A 4x4 symmetric quadric, as its upper triangle.
///
/// `v^T Q v` is the sum of squared distances from `v` to a set of planes, which
/// is what makes these addable: the error of merging two vertices is the error
/// of one plus the error of the other.
#[derive(Debug, Clone, Copy, Default)]
struct Quadric([f64; 10]);

impl Quadric {
    /// The quadric of one plane `ax + by + cz + d = 0`, with a unit normal.
    fn of_plane(a: f64, b: f64, c: f64, d: f64) -> Self {
        Self([
            a * a,
            a * b,
            a * c,
            a * d,
            b * b,
            b * c,
            b * d,
            c * c,
            c * d,
            d * d,
        ])
    }

    fn add(&mut self, other: &Self) {
        for (a, b) in self.0.iter_mut().zip(other.0) {
            *a += b;
        }
    }

    /// The total weight of the planes summed in: each plane's unit normal
    /// puts its weight once on the diagonal.
    fn weight(&self) -> f64 {
        self.0[0] + self.0[4] + self.0[7]
    }

    /// The squared distance this quadric assigns to a point.
    fn at(&self, p: Point) -> f64 {
        let [q00, q01, q02, q03, q11, q12, q13, q22, q23, q33] = self.0;
        let (x, y, z) = (p.x, p.y, p.z);
        q00 * x * x
            + 2.0 * q01 * x * y
            + 2.0 * q02 * x * z
            + 2.0 * q03 * x
            + q11 * y * y
            + 2.0 * q12 * y * z
            + 2.0 * q13 * y
            + q22 * z * z
            + 2.0 * q23 * z
            + q33
    }
}

/// The mesh mid-decimation.
struct State {
    positions: Vec<Point>,
    triangles: Vec<[u32; 3]>,
    /// Whether each triangle is still there.
    live: Vec<bool>,
    quadrics: Vec<Quadric>,
    /// Vertices that must not move: the ones on a boundary.
    pinned: Vec<bool>,
    /// Further vertices that must not move: the ones a verified decimation's
    /// measurement held.
    held: Vec<bool>,
    /// The vertex each vertex was merged into; itself while it survives.
    merged_into: Vec<u32>,
    /// Which triangles touch each vertex, live or not; empty for a vertex
    /// no triangle names or one merged away.
    around: Vec<Vec<u32>>,
    tol: Tolerances,
    remaining: usize,
}

impl State {
    fn new(mesh: &Triangulation, tol: Tolerances) -> Self {
        let mut quadrics = vec![Quadric::default(); mesh.positions.len()];
        let mut around: Vec<Vec<u32>> = vec![Vec::new(); mesh.positions.len()];
        let mut uses: FastMap<(u32, u32), usize> = FastMap::default();

        for (i, triangle) in mesh.triangles.iter().enumerate() {
            let [a, b, c] = triangle.map(|v| mesh.positions[v as usize]);
            let normal = (b - a).cross(c - a);
            let length = normal.magnitude();
            if length > tol.confusion() {
                let unit = normal * (1.0 / length);
                // Weighted by area, so a large flat region is not outvoted by a
                // cluster of slivers describing the same plane.
                let plane = Quadric::of_plane(unit.x, unit.y, unit.z, -unit.dot(a.to_vector()));
                let mut weighted = plane;
                for value in &mut weighted.0 {
                    *value *= length;
                }
                for v in triangle {
                    quadrics[*v as usize].add(&weighted);
                }
            }
            for v in triangle {
                #[allow(clippy::cast_possible_truncation)]
                around[*v as usize].push(i as u32);
            }
            for k in 0..3 {
                let (x, y) = (triangle[k], triangle[(k + 1) % 3]);
                *uses.entry((x.min(y), x.max(y))).or_default() += 1;
            }
        }

        // A boundary edge is used once. Both its ends are held.
        let mut pinned = vec![false; mesh.positions.len()];
        for ((a, b), count) in uses {
            if count != 2 {
                pinned[a as usize] = true;
                pinned[b as usize] = true;
            }
        }

        Self {
            positions: mesh.positions.clone(),
            triangles: mesh.triangles.clone(),
            live: vec![true; mesh.triangles.len()],
            quadrics,
            pinned,
            held: vec![false; mesh.positions.len()],
            #[allow(clippy::cast_possible_truncation)]
            merged_into: (0..mesh.positions.len() as u32).collect(),
            around,
            tol,
            remaining: mesh.triangles.len(),
        }
    }

    const fn live_triangles(&self) -> usize {
        self.remaining
    }

    /// The cheapest collapse still worth making, within a squared-error budget,
    /// ties to the edge met first walking the live triangles in order.
    ///
    /// The whole mesh is walked at every call. [`State::decimate_by_queue`]
    /// takes the same collapses in the same order from a queue.
    fn cheapest(&self, budget: f64) -> Option<(f64, u32, u32, Point)> {
        let mut best: Option<(f64, u32, u32, Point)> = None;
        let mut seen: FastSet<(u32, u32)> = FastSet::default();
        for (i, triangle) in self.triangles.iter().enumerate() {
            if !self.live[i] {
                continue;
            }
            for k in 0..3 {
                let (a, b) = (triangle[k], triangle[(k + 1) % 3]);
                let key = (a.min(b), a.max(b));
                if !seen.insert(key) {
                    continue;
                }
                // Either end pinned means the edge is on a boundary, or leads
                // to one. Holding both ends holds the outline exactly.
                if self.pinned[a as usize] || self.pinned[b as usize] {
                    continue;
                }
                let at = Point::from_vector(
                    (self.positions[a as usize].to_vector()
                        + self.positions[b as usize].to_vector())
                        * 0.5,
                );
                let mut merged = self.quadrics[a as usize];
                merged.add(&self.quadrics[b as usize]);
                let cost = merged.at(at).max(0.0);
                if cost > budget {
                    continue;
                }
                if best.is_some_and(|(current, ..)| cost >= current) {
                    continue;
                }
                if self.would_fold(a, b, at) || !self.keeps_manifold(a, b) {
                    continue;
                }
                best = Some((cost, a, b, at));
            }
        }
        best
    }

    /// Collapse the cheapest edge within `budget`, by [`State::cheapest`],
    /// while more than `floor` triangles remain: the worst cost accepted
    /// and how many collapses were made.
    fn decimate_by_scan(&mut self, budget: f64, floor: usize) -> (f64, usize) {
        let mut worst = 0.0_f64;
        let mut collapsed = 0;
        while self.live_triangles() > floor {
            let Some((cost, from, to, at)) = self.cheapest(budget) else {
                break;
            };
            self.collapse(from, to, at);
            worst = worst.max(cost);
            collapsed += 1;
        }
        (worst, collapsed)
    }

    /// [`State::decimate_by_scan`] from a queue: the same collapses in the
    /// same order, each found without walking the whole mesh.
    ///
    /// The queue holds every collapse [`State::cheapest`] would weigh, at
    /// the weight it would give it at this step, and stale entries
    /// besides. A collapse
    /// changes the weighing only of the edges with an end among the
    /// vertices about either end of the collapsed edge, before or after,
    /// or the merged vertex itself: their cost, the first live triangle
    /// holding them, and whether they fold or break the manifold read
    /// nothing else. Those are weighed again and queued. An entry is taken
    /// only where its edge weighs the same when it comes up, so the first
    /// taken is the cheapest edge, ties to the one the scan meets first,
    /// which is what the scan takes. The queue is built afresh once it
    /// holds more than three entries a live triangle, half of them stale
    /// at least.
    fn decimate_by_queue(&mut self, budget: f64, floor: usize) -> (f64, usize) {
        let mut worst = 0.0_f64;
        let mut collapsed = 0;
        let mut queue = self.queue(budget);
        let mut fresh: FastSet<(u32, u32)> = FastSet::default();
        while self.live_triangles() > floor {
            let Some(Reverse(taken)) = queue.pop() else {
                break;
            };
            let Some((now, at)) = self.candidate(taken.from, taken.to, budget) else {
                continue;
            };
            if !now.same(&taken) {
                continue;
            }
            let mut near = self.ring(now.from);
            near.extend(self.ring(now.to));
            self.collapse(now.from, now.to, at);
            near.extend(self.ring(now.to));
            near.push(now.to);
            near.sort_unstable();
            near.dedup();
            fresh.clear();
            for &v in &near {
                for w in self.ring(v) {
                    if fresh.insert((v.min(w), v.max(w)))
                        && let Some((again, _)) = self.candidate(v, w, budget)
                    {
                        queue.push(Reverse(again));
                    }
                }
            }
            worst = worst.max(now.cost);
            collapsed += 1;
            // Each live triangle holds three edge sides and each live edge
            // two: past three entries a triangle, at least half are stale.
            if queue.len() > 3 * self.live_triangles() {
                queue = self.queue(budget);
            }
        }
        (worst, collapsed)
    }

    /// Every collapse [`State::cheapest`] would weigh now, queued.
    fn queue(&self, budget: f64) -> BinaryHeap<Reverse<Candidate>> {
        let mut queue = BinaryHeap::new();
        let mut seen: FastSet<(u32, u32)> = FastSet::default();
        for (i, triangle) in self.triangles.iter().enumerate() {
            if !self.live[i] {
                continue;
            }
            for k in 0..3 {
                let (a, b) = (triangle[k], triangle[(k + 1) % 3]);
                if seen.insert((a.min(b), a.max(b)))
                    && let Some((c, _)) = self.candidate(a, b, budget)
                {
                    queue.push(Reverse(c));
                }
            }
        }
        queue
    }

    /// The collapse of the edge between `a` and `b` as [`State::cheapest`]
    /// weighs it, and the point it merges at: `None` where no live triangle
    /// holds the edge, an end is pinned, the cost is past `budget`, or the
    /// collapse would fold or break the manifold.
    fn candidate(&self, a: u32, b: u32, budget: f64) -> Option<(Candidate, Point)> {
        if self.pinned[a as usize] || self.pinned[b as usize] {
            return None;
        }
        let mut first: Option<(u64, u32, u32)> = None;
        for &index in &self.around[a as usize] {
            if !self.live[index as usize] {
                continue;
            }
            let t = self.triangles[index as usize];
            for k in 0..3 {
                let (x, y) = (t[k], t[(k + 1) % 3]);
                if (x == a && y == b) || (x == b && y == a) {
                    let place = 3 * u64::from(index) + k as u64;
                    if first.is_none_or(|f| place < f.0) {
                        first = Some((place, x, y));
                    }
                }
            }
        }
        let (place, from, to) = first?;
        let at = Point::from_vector(
            (self.positions[from as usize].to_vector() + self.positions[to as usize].to_vector())
                * 0.5,
        );
        let mut merged = self.quadrics[from as usize];
        merged.add(&self.quadrics[to as usize]);
        let cost = merged.at(at).max(0.0);
        if cost > budget || self.would_fold(from, to, at) || !self.keeps_manifold(from, to) {
            return None;
        }
        Some((
            Candidate {
                cost,
                place,
                from,
                to,
            },
            at,
        ))
    }

    /// Whether merging two vertices would turn any surviving triangle over.
    ///
    /// A fold is not a small error. It is a mesh that no longer bounds what it
    /// did, and no error metric measures that: the merged point can sit
    /// exactly on every plane and still put the triangle back to front.
    fn would_fold(&self, from: u32, to: u32, at: Point) -> bool {
        for vertex in [from, to] {
            for &index in &self.around[vertex as usize] {
                if !self.live[index as usize] {
                    continue;
                }
                let triangle = self.triangles[index as usize];
                // Triangles containing both vanish with the collapse.
                if triangle.contains(&from) && triangle.contains(&to) {
                    continue;
                }
                let before = self.normal_of(triangle, None);
                let after = self.normal_of(triangle, Some((from, to, at)));
                let (Some(before), Some(after)) = (before, after) else {
                    return true;
                };
                if before.dot(after) <= 0.0 {
                    return true;
                }
            }
        }
        false
    }

    /// Whether merging two vertices keeps the faces about them a manifold:
    /// the vertices they both neighbour are those across the edge between
    /// them, and no two faces would be left on the same three vertices.
    ///
    /// A third neighbour in common joins two edges into one with four
    /// faces; two faces on the same three vertices, as when a tetrahedron
    /// collapses, close over each other. Neither is a coarser surface.
    fn keeps_manifold(&self, a: u32, b: u32) -> bool {
        let mut across = Vec::new();
        for &index in &self.around[a as usize] {
            let t = self.triangles[index as usize];
            if self.live[index as usize] && t.contains(&b) {
                across.extend(t.iter().filter(|&&w| w != a && w != b));
            }
        }
        across.sort_unstable();
        across.dedup();
        let (ring_a, ring_b) = (self.ring(a), self.ring(b));
        let common = ring_a.iter().filter(|w| ring_b.binary_search(w).is_ok());
        if !common.copied().eq(across.iter().copied()) {
            return false;
        }
        let face = |v: u32, x: u32, y: u32| {
            self.around[v as usize].iter().any(|&index| {
                let t = self.triangles[index as usize];
                self.live[index as usize] && t.contains(&x) && t.contains(&y)
            })
        };
        !(across.len() == 2 && face(a, across[0], across[1]) && face(b, across[0], across[1]))
    }

    /// A triangle's normal, optionally with one collapse applied.
    fn normal_of(
        &self,
        triangle: [u32; 3],
        collapse: Option<(u32, u32, Point)>,
    ) -> Option<ogeom_math::Vector> {
        let at = |v: u32| match collapse {
            Some((from, to, p)) if v == from || v == to => p,
            _ => self.positions[v as usize],
        };
        let (a, b, c) = (at(triangle[0]), at(triangle[1]), at(triangle[2]));
        let normal = (b - a).cross(c - a);
        if normal.magnitude() <= self.tol.confusion() {
            return None;
        }
        Some(normal)
    }

    /// Merge two vertices at a point.
    fn collapse(&mut self, from: u32, to: u32, at: Point) {
        self.positions[to as usize] = at;
        self.merged_into[from as usize] = to;
        let mut merged = self.quadrics[from as usize];
        merged.add(&self.quadrics[to as usize]);
        self.quadrics[to as usize] = merged;

        let touching = std::mem::take(&mut self.around[from as usize]);
        let mut emptied = Vec::new();
        for index in touching {
            if !self.live[index as usize] {
                continue;
            }
            let triangle = &mut self.triangles[index as usize];
            for v in triangle.iter_mut() {
                if *v == from {
                    *v = to;
                }
            }
            // A triangle naming one vertex twice has no area left.
            let [a, b, c] = *triangle;
            if a == b || b == c || c == a {
                self.live[index as usize] = false;
                self.remaining -= 1;
                emptied.extend([a, b, c]);
            } else {
                self.around[to as usize].push(index);
            }
        }
        // The vertices about a triangle that went drop it, the order of the
        // rest kept, so that a vertex that has taken in many others does
        // not walk their dead triangles at every look.
        emptied.sort_unstable();
        emptied.dedup();
        for v in emptied {
            let live = &self.live;
            self.around[v as usize].retain(|&t| live[t as usize]);
        }
    }

    /// The surviving mesh, with unused vertices dropped, and each surviving
    /// vertex's number in it.
    fn harvest(&self) -> (Triangulation, FastMap<u32, u32>) {
        let mut out = Triangulation::new();
        let mut moved: FastMap<u32, u32> = FastMap::default();
        for (index, triangle) in self.triangles.iter().enumerate() {
            if !self.live[index] {
                continue;
            }
            let mut mapped = [0_u32; 3];
            for (slot, v) in mapped.iter_mut().zip(triangle) {
                *slot = *moved.entry(*v).or_insert_with(|| {
                    #[allow(clippy::cast_possible_truncation)]
                    let fresh = out.positions.len() as u32;
                    out.positions.push(self.positions[*v as usize]);
                    fresh
                });
            }
            out.triangles.push(mapped);
        }
        // Normals are recomputed from the surviving geometry rather than
        // carried over: a merged vertex's old normal described a surface that
        // is no longer there.
        out.normals = vec![ogeom_math::Vector::ZERO; out.positions.len()];
        for triangle in &out.triangles {
            let [a, b, c] = triangle.map(|v| out.positions[v as usize]);
            let normal = (b - a).cross(c - a);
            for v in triangle {
                out.normals[*v as usize] += normal;
            }
        }
        for normal in &mut out.normals {
            let length = normal.magnitude();
            if length > self.tol.confusion() {
                *normal *= 1.0 / length;
            }
        }
        out.parameters = vec![(0.0, 0.0); out.positions.len()];
        out.deflection_met = false;
        (out, moved)
    }

    /// The vertex a vertex of the source ended up as.
    fn survivor(&self, mut v: u32) -> u32 {
        while self.merged_into[v as usize] != v {
            v = self.merged_into[v as usize];
        }
        v
    }

    /// For each source vertex, the number in the harvested mesh of the
    /// vertex it ended up as, given each survivor's number there.
    fn result_vertex(&self, renumbered: &FastMap<u32, u32>) -> Vec<Option<u32>> {
        let mut root = self.merged_into.clone();
        for v in 0..root.len() {
            let mut r = root[v];
            while root[r as usize] != r {
                r = root[r as usize];
            }
            let mut w = v;
            while root[w] != r {
                let up = root[w] as usize;
                root[w] = r;
                w = up;
            }
        }
        root.iter().map(|r| renumbered.get(r).copied()).collect()
    }

    /// Whether an edge may collapse at all.
    fn free(&self, a: u32, b: u32) -> bool {
        !(self.pinned[a as usize]
            || self.pinned[b as usize]
            || self.held[a as usize]
            || self.held[b as usize])
    }

    /// The cost of collapsing an edge to its midpoint, as a squared distance:
    /// the mean of the squared distances to the planes merged in, weighted as
    /// they were summed.
    ///
    /// Divided by the weight so that it does not grow with the size of the
    /// triangles, which [`Target::Verified`] needs to set a budget in the
    /// units of its tolerance. A mean can hide one far plane, so the budget
    /// only orders and filters; the measurement afterwards decides.
    fn mean_cost(&self, a: u32, b: u32) -> (f64, Point) {
        let at = Point::from_vector(
            (self.positions[a as usize].to_vector() + self.positions[b as usize].to_vector()) * 0.5,
        );
        let mut merged = self.quadrics[a as usize];
        merged.add(&self.quadrics[b as usize]);
        let cost = merged.at(at).max(0.0);
        let weight = merged.weight();
        (if weight > 0.0 { cost / weight } else { cost }, at)
    }

    /// The vertices sharing a live triangle with `v`, in ascending order.
    fn ring(&self, v: u32) -> Vec<u32> {
        let mut out = Vec::new();
        for &index in &self.around[v as usize] {
            if self.live[index as usize] {
                out.extend(self.triangles[index as usize].iter().filter(|&&w| w != v));
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::{Deflection, triangulate};
    use ogeom_algo::{make_box, make_sphere};
    use ogeom_math::Frame;
    use ogeom_topo::Model;

    const T: Tolerances = Tolerances::millimetres();

    fn sphere(chord: f64) -> Triangulation {
        let mut model = Model::new();
        let built = make_sphere(&mut model, Frame::WORLD, 10.0, T).unwrap();
        triangulate(
            &model,
            &built.shape,
            Deflection {
                chord,
                ..Deflection::default()
            },
            T,
        )
        .unwrap()
    }

    #[test]
    fn decimating_a_sphere_keeps_it_a_sphere_to_the_error_it_reports() {
        // The property that matters: the result says how far it moved, and
        // every vertex of it really is within that of the original surface.
        let mesh = sphere(0.02);
        let before = mesh.triangle_count();
        let done = simplify(&mesh, Target::Triangles(before / 4), T).unwrap();

        assert!(done.collapsed > 0);
        assert!(
            done.mesh.triangle_count() < before,
            "nothing was removed: {} of {before}",
            done.mesh.triangle_count()
        );
        for p in &done.mesh.positions {
            let off = (p.to_vector().magnitude() - 10.0).abs();
            assert!(
                off <= done.error + 1e-9,
                "a vertex is {off} off the sphere, but the reported error is {}",
                done.error
            );
        }
    }

    #[test]
    fn a_tighter_error_budget_removes_less() {
        let mesh = sphere(0.02);
        let loose = simplify(&mesh, Target::Error(0.5), T).unwrap();
        let tight = simplify(&mesh, Target::Error(0.01), T).unwrap();

        assert!(
            tight.mesh.triangle_count() >= loose.mesh.triangle_count(),
            "a tighter budget should keep more: {} against {}",
            tight.mesh.triangle_count(),
            loose.mesh.triangle_count()
        );
        assert!(
            tight.error <= 0.01 + 1e-12,
            "over budget at {}",
            tight.error
        );
        assert!(loose.error <= 0.5 + 1e-12);
        assert!(tight.target_met && loose.target_met);
    }

    #[test]
    fn the_mesh_stays_closed() {
        // A collapse that opened a hole would be a decimation that changed what
        // the mesh bounds, which is a different thing from making it coarser.
        let mesh = sphere(0.05);
        assert!(mesh.is_closed());
        let done = simplify(&mesh, Target::Triangles(mesh.triangle_count() / 2), T).unwrap();
        assert!(
            done.mesh.is_closed(),
            "decimation opened the mesh after {} collapses",
            done.collapsed
        );
        assert!(done.mesh.volume() > 0.0, "and it turned inside out");
    }

    #[test]
    fn a_collapse_that_would_join_two_edges_into_one_is_not_made() {
        // A double pyramid on a triangle with two of its corners close
        // together: merging those two leaves each pole's three faces on two
        // vertices, every edge of the pair under four faces.
        let mut mesh = Triangulation::new();
        mesh.positions = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(0.01, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(0.0, 0.3, 1.0),
            Point::new(0.0, 0.3, -1.0),
        ];
        mesh.triangles = vec![
            [3, 0, 1],
            [3, 1, 2],
            [3, 2, 0],
            [4, 1, 0],
            [4, 2, 1],
            [4, 0, 2],
        ];
        for target in [Target::Triangles(4), Target::Error(1.0)] {
            let done = simplify(&mesh, target, T).unwrap();
            let mut uses: FastMap<(u32, u32), usize> = FastMap::default();
            let mut sets = FastSet::default();
            for t in &done.mesh.triangles {
                let mut set = *t;
                set.sort_unstable();
                assert!(sets.insert(set), "two faces on {set:?} after {target:?}");
                for k in 0..3 {
                    let (a, b) = (t[k], t[(k + 1) % 3]);
                    *uses.entry((a.min(b), a.max(b))).or_default() += 1;
                }
            }
            assert!(
                uses.values().all(|&n| n == 2),
                "an edge without two faces after {target:?}: {uses:?}"
            );
            assert!(done.mesh.volume() > 0.0);
        }
    }

    #[test]
    fn a_boundary_is_held_exactly() {
        // A sheet body's outline must not creep inward as the mesh coarsens,
        // and "expensive but possible" would let it.
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 4.0), T).unwrap();
        let face = ogeom_topo::explore_unique(&model, &solid.shape, ogeom_topo::ShapeType::Face)
            .unwrap()[0]
            .clone();
        let sheet = crate::triangulate_face(
            &model,
            &face,
            Deflection {
                chord: 0.05,
                ..Deflection::default()
            },
            T,
        )
        .unwrap();

        let outline = |m: &Triangulation| {
            let mut low = f64::MAX;
            let mut high = f64::MIN;
            for p in &m.positions {
                low = low.min(p.x);
                high = high.max(p.x);
            }
            (low, high)
        };
        let before = outline(&sheet);
        let done = simplify(&sheet, Target::Triangles(2), T).unwrap();
        let after = outline(&done.mesh);
        assert!(
            (before.0 - after.0).abs() < 1e-12 && (before.1 - after.1).abs() < 1e-12,
            "the outline moved from {before:?} to {after:?}"
        );
    }

    #[test]
    fn a_target_that_cannot_be_reached_is_reported_rather_than_claimed() {
        // Every edge of a single triangle is a boundary edge, so there is
        // nothing to collapse. Saying the target was met would be a lie that
        // costs a caller the chance to try something else.
        let mut mesh = Triangulation::new();
        mesh.positions = vec![
            Point::ORIGIN,
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        mesh.triangles = vec![[0, 1, 2]];
        let done = simplify(&mesh, Target::Triangles(1), T).unwrap();
        assert_eq!(done.collapsed, 0);
        assert_eq!(done.mesh.triangle_count(), 1);
    }

    #[test]
    fn a_target_that_describes_nothing_is_refused() {
        let mesh = sphere(0.2);
        assert!(simplify(&mesh, Target::Triangles(0), T).is_err());
        assert!(simplify(&mesh, Target::Error(0.0), T).is_err());
        assert!(simplify(&mesh, Target::Error(-1.0), T).is_err());
        assert!(simplify(&mesh, Target::Error(f64::NAN), T).is_err());

        let mut broken = Triangulation::new();
        broken.positions = vec![Point::ORIGIN];
        broken.triangles = vec![[0, 1, 2]];
        assert!(simplify(&broken, Target::Triangles(1), T).is_err());
    }

    /// A torus of `n` by `m` grid squares, each cut into two triangles,
    /// its tube's radius waved by `wave`; with `open`, one ring of squares
    /// left out, so the grid has a border to hold.
    fn torus_grid(n: usize, m: usize, wave: f64, open: bool) -> Triangulation {
        use core::f64::consts::TAU;
        let mut grid = Triangulation::new();
        for i in 0..n {
            for j in 0..m {
                #[allow(clippy::cast_precision_loss)]
                let (u, v) = (TAU * i as f64 / n as f64, TAU * j as f64 / m as f64);
                let r = wave.mul_add((3.0 * u).sin() * (2.0 * v).cos(), 5.0);
                let rho = r.mul_add(v.cos(), 20.0);
                grid.positions
                    .push(Point::new(rho * u.cos(), rho * u.sin(), r * v.sin()));
            }
        }
        #[allow(clippy::cast_possible_truncation)]
        let at = |i: usize, j: usize| ((i % n) * m + (j % m)) as u32;
        for i in 0..(if open { n - 1 } else { n }) {
            for j in 0..m {
                let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
                grid.triangles.push([a, b, c]);
                grid.triangles.push([a, c, d]);
            }
        }
        grid
    }

    #[test]
    fn the_queue_makes_the_collapses_the_scan_makes() {
        // The queue only finds the scan's choice sooner: on an even grid,
        // where many edges cost the same to the bit and ties go to the
        // edge met first, on one with a border, and on a mesher's sphere,
        // both make the same collapses to the bit.
        let meshes = [
            torus_grid(40, 20, 0.0, false),
            torus_grid(40, 20, 0.5, false),
            torus_grid(40, 20, 0.5, true),
            sphere(0.05),
        ];
        for mesh in &meshes {
            // A budget the halving reached, so a run held to it collapses
            // some edges and stops on the budget, not the count.
            let halved = State::new(mesh, T).decimate_by_scan(f64::MAX, mesh.triangles.len() / 2);
            for (budget, floor) in [
                (f64::MAX, mesh.triangles.len() / 2),
                (f64::MAX, mesh.triangles.len() / 10),
                (halved.0, 1),
            ] {
                let mut scan = State::new(mesh, T);
                let mut queue = State::new(mesh, T);
                let by_scan = scan.decimate_by_scan(budget, floor);
                let by_queue = queue.decimate_by_queue(budget, floor);
                assert!(by_scan.1 > 0);
                assert_eq!(by_scan.0.to_bits(), by_queue.0.to_bits());
                assert_eq!(by_scan.1, by_queue.1);
                assert_eq!(scan.triangles, queue.triangles);
                assert_eq!(scan.live, queue.live);
                assert_eq!(scan.merged_into, queue.merged_into);
                for (a, b) in scan.positions.iter().zip(&queue.positions) {
                    assert_eq!(
                        [a.x, a.y, a.z].map(f64::to_bits),
                        [b.x, b.y, b.z].map(f64::to_bits)
                    );
                }
            }
        }
    }
}
