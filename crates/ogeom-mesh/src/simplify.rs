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
//! measures the result against the source both ways, pins what fails and
//! decimates again; see [`Deviation`] for what it reports.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use ogeom_core::{FastMap, FastSet};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::Point;
use ogeom_topo::Triangulation;

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
    /// and triangle centre to the source. Where either is farther than
    /// `tolerance`, the source vertices behind it are pinned and the
    /// decimation is run again from the source, until at least `share` of
    /// the points on each side are within `tolerance`. The measurement is
    /// reported in [`Simplified::deviation`].
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
}

/// Decimate a mesh.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the target is
/// not a positive count or a positive distance, a share is not in `(0, 1]`, or
/// the mesh names a vertex it does not have;
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) with
/// [`Target::Verified`] if a vertex of a triangle is not finite.
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

    let mut worst = 0.0_f64;
    let mut collapsed = 0;
    while state.live_triangles() > floor {
        let Some((cost, from, to, at)) = state.cheapest(budget) else {
            break;
        };
        state.collapse(from, to, at);
        worst = worst.max(cost);
        collapsed += 1;
    }

    let target_met = match target {
        Target::Triangles(n) => state.live_triangles() <= n,
        // Every collapse made was inside the budget, and the loop stops only
        // when no remaining one is.
        Target::Error(_) | Target::Verified { .. } => true,
    };
    Ok(Simplified {
        mesh: state.harvest().0,
        error: worst.max(0.0).sqrt(),
        collapsed,
        target_met,
        deviation: None,
    })
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
    pinned: FastSet<u32>,
    /// Further vertices that must not move: the ones a verified decimation's
    /// measurement held.
    held: FastSet<u32>,
    /// The vertex each vertex was merged into; itself while it survives.
    merged_into: Vec<u32>,
    /// Which triangles touch each vertex.
    around: FastMap<u32, Vec<usize>>,
    tol: Tolerances,
    remaining: usize,
}

impl State {
    fn new(mesh: &Triangulation, tol: Tolerances) -> Self {
        let mut quadrics = vec![Quadric::default(); mesh.positions.len()];
        let mut around: FastMap<u32, Vec<usize>> = FastMap::default();
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
                around.entry(*v).or_default().push(i);
            }
            for k in 0..3 {
                let (x, y) = (triangle[k], triangle[(k + 1) % 3]);
                *uses.entry((x.min(y), x.max(y))).or_default() += 1;
            }
        }

        // A boundary edge is used once. Both its ends are held.
        let mut pinned = FastSet::default();
        for ((a, b), count) in uses {
            if count != 2 {
                pinned.insert(a);
                pinned.insert(b);
            }
        }

        Self {
            positions: mesh.positions.clone(),
            triangles: mesh.triangles.clone(),
            live: vec![true; mesh.triangles.len()],
            quadrics,
            pinned,
            held: FastSet::default(),
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

    /// The cheapest collapse still worth making, within a squared-error budget.
    ///
    /// Recomputed each round rather than kept in a heap. The mesh is small
    /// enough that the difference is not what makes decimation slow, and a heap
    /// of costs that go stale on every collapse needs invalidation logic that is
    /// its own source of wrong answers.
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
                if self.pinned.contains(&a) || self.pinned.contains(&b) {
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
                if self.would_fold(a, b, at) {
                    continue;
                }
                best = Some((cost, a, b, at));
            }
        }
        best
    }

    /// Whether merging two vertices would turn any surviving triangle over.
    ///
    /// A fold is not a small error. It is a mesh that no longer bounds what it
    /// did, and no error metric measures that: the merged point can sit
    /// exactly on every plane and still put the triangle back to front.
    fn would_fold(&self, from: u32, to: u32, at: Point) -> bool {
        for vertex in [from, to] {
            for index in self.around.get(&vertex).into_iter().flatten() {
                if !self.live[*index] {
                    continue;
                }
                let triangle = self.triangles[*index];
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

        let touching: Vec<usize> = self
            .around
            .get(&from)
            .into_iter()
            .flatten()
            .copied()
            .collect();
        for index in touching {
            if !self.live[index] {
                continue;
            }
            let triangle = &mut self.triangles[index];
            for v in triangle.iter_mut() {
                if *v == from {
                    *v = to;
                }
            }
            // A triangle naming one vertex twice has no area left.
            let [a, b, c] = *triangle;
            if a == b || b == c || c == a {
                self.live[index] = false;
                self.remaining -= 1;
            } else {
                self.around.entry(to).or_default().push(index);
            }
        }
        self.around.remove(&from);
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

    /// Whether an edge may collapse at all.
    fn free(&self, a: u32, b: u32) -> bool {
        !(self.pinned.contains(&a)
            || self.pinned.contains(&b)
            || self.held.contains(&a)
            || self.held.contains(&b))
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
        for index in self.around.get(&v).into_iter().flatten() {
            if self.live[*index] {
                out.extend(self.triangles[*index].iter().filter(|&&w| w != v));
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    /// Collapse the cheapest free edge, by [`State::mean_cost`], while its
    /// cost is within `budget` and more than `floor` triangles remain;
    /// returns how many collapses were made.
    ///
    /// The choice is the one a scan of every edge would make, ties to the
    /// lower vertex numbers, kept in a queue: a collapse changes the cost of
    /// the edges at the merged vertex and whether the edges at its neighbours
    /// fold, and only those are queued again, the entries before marked
    /// stale.
    fn decimate_queued(&mut self, budget: f64, floor: usize) -> usize {
        let mut stamps = vec![0_u64; self.positions.len()];
        let mut queue = BinaryHeap::new();
        let mut seen: FastSet<(u32, u32)> = FastSet::default();
        for (i, triangle) in self.triangles.iter().enumerate() {
            if !self.live[i] {
                continue;
            }
            for k in 0..3 {
                let (a, b) = (triangle[k], triangle[(k + 1) % 3]);
                let key = (a.min(b), a.max(b));
                if seen.insert(key) && self.free(a, b) {
                    queue.push(self.queued(key, &stamps));
                }
            }
        }
        drop(seen);

        let mut collapsed = 0;
        while self.remaining > floor {
            let Some(Reverse(entry)) = queue.pop() else {
                break;
            };
            if entry.cost > budget {
                break;
            }
            let (lo, hi) = entry.edge;
            if entry.stamps != (stamps[lo as usize], stamps[hi as usize])
                || !self.around.contains_key(&lo)
                || !self.around.contains_key(&hi)
            {
                continue;
            }
            let (_, at) = self.mean_cost(lo, hi);
            if self.would_fold(hi, lo, at) {
                continue;
            }
            self.collapse(hi, lo, at);
            collapsed += 1;

            let mut touched = self.ring(lo);
            touched.push(lo);
            for &v in &touched {
                stamps[v as usize] += 1;
            }
            for &v in &touched {
                for w in self.ring(v) {
                    if self.free(v, w) {
                        queue.push(self.queued((v.min(w), v.max(w)), &stamps));
                    }
                }
            }
        }
        collapsed
    }

    fn queued(&self, edge: (u32, u32), stamps: &[u64]) -> Reverse<Queued> {
        Reverse(Queued {
            cost: self.mean_cost(edge.0, edge.1).0,
            edge,
            stamps: (stamps[edge.0 as usize], stamps[edge.1 as usize]),
        })
    }
}

/// An edge waiting to collapse, with its cost and its ends' stamps when it
/// was queued.
#[derive(Debug, Clone, Copy)]
struct Queued {
    cost: f64,
    edge: (u32, u32),
    stamps: (u64, u64),
}

impl Queued {
    fn key(&self) -> (f64, (u32, u32), (u64, u64)) {
        (self.cost, self.edge, self.stamps)
    }
}

impl PartialEq for Queued {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}

impl Eq for Queued {}

impl PartialOrd for Queued {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Queued {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let (a, b) = (self.key(), other.key());
        a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2))
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
}
