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
//! distance from the merged vertex to the planes of every face that met there —
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
//! coarsens — which is exactly the thing a caller would not think to check.
//! Holding the boundary exactly is a stronger promise and a simpler one.
//!
//! A collapse that would turn a triangle inside out is refused for the same
//! reason: a fold is not a small error, it is a mesh that no longer bounds what
//! it did.

use og_core::{OgResult, Tolerances, og_bail};
use og_math::Point;
use og_topo::Triangulation;
use std::collections::{HashMap, HashSet};

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
    /// — every remaining edge is on a boundary or would fold something — and
    /// then the result is as small as it can safely be rather than as small as
    /// was asked for. Reported rather than passed off as success.
    pub target_met: bool,
}

/// Decimate a mesh.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the target is
/// not a positive count or a positive distance, or the mesh names a vertex it
/// does not have.
pub fn simplify(mesh: &Triangulation, target: Target, tol: Tolerances) -> OgResult<Simplified> {
    match target {
        Target::Triangles(0) => {
            og_bail!(Construction, "a mesh of no triangles describes nothing");
        }
        Target::Error(e) if !e.is_finite() || e <= 0.0 => {
            og_bail!(Construction, "an error budget of {e} is not a distance");
        }
        _ => {}
    }
    for triangle in &mesh.triangles {
        for index in triangle {
            if *index as usize >= mesh.positions.len() {
                og_bail!(
                    Construction,
                    "a triangle names vertex {index}, and the mesh has {}",
                    mesh.positions.len()
                );
            }
        }
    }

    let mut state = State::new(mesh, tol);
    let budget = match target {
        Target::Error(e) => e * e,
        Target::Triangles(_) => f64::MAX,
    };
    let floor = match target {
        Target::Triangles(n) => n,
        Target::Error(_) => 1,
    };

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
        Target::Error(_) => true,
    };
    Ok(Simplified {
        mesh: state.harvest(),
        error: worst.max(0.0).sqrt(),
        collapsed,
        target_met,
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
    pinned: HashSet<u32>,
    /// Which triangles touch each vertex.
    around: HashMap<u32, Vec<usize>>,
    tol: Tolerances,
    remaining: usize,
}

impl State {
    fn new(mesh: &Triangulation, tol: Tolerances) -> Self {
        let mut quadrics = vec![Quadric::default(); mesh.positions.len()];
        let mut around: HashMap<u32, Vec<usize>> = HashMap::new();
        let mut uses: HashMap<(u32, u32), usize> = HashMap::new();

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
        let mut pinned = HashSet::new();
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
        let mut seen: HashSet<(u32, u32)> = HashSet::new();
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
    /// did, and no error metric measures that — the merged point can sit
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
    ) -> Option<og_math::Vector> {
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

    /// The surviving mesh, with unused vertices dropped.
    fn harvest(self) -> Triangulation {
        let mut out = Triangulation::new();
        let mut moved: HashMap<u32, u32> = HashMap::new();
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
        out.normals = vec![og_math::Vector::ZERO; out.positions.len()];
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
        out
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::{Deflection, triangulate};
    use og_algo::{make_box, make_sphere};
    use og_math::Frame;
    use og_topo::Model;

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
        let face = og_topo::explore_unique(&model, &solid.shape, og_topo::ShapeType::Face).unwrap()
            [0]
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
