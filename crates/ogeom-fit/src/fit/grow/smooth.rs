//! The smooth free-form region of a mesh about points clicked on it, for a
//! surface fitted to it later: [`grow_smooth`].

use std::collections::VecDeque;

use ogeom_core::{FastMap, OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Point, Vector};
use ogeom_mesh::Bvh;
use ogeom_topo::Triangulation;

use super::{GrowOptions, NONE, as_index};
use crate::Cutoff;
use crate::orient::check_indices;
use ogeom_mesh::weld_points;

/// A smooth region grown over a mesh, with no surface fitted to it.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct GrownRegion {
    /// The region's triangles, in increasing order.
    pub triangles: Vec<u32>,
    /// The region's vertices, in increasing order. Of vertices that lie
    /// within the confusion distance of each other, only the lowest is
    /// named, so a point on a seam between patches meshed apart is named
    /// once.
    pub vertices: Vec<u32>,
    /// The position of each vertex, in the same order: the points a
    /// free-form surface is fitted to.
    pub points: Vec<Point>,
    /// The triangle nearest each seed, in the seeds' order.
    pub seeds: Vec<u32>,
    /// The passes [`GrowOptions::peel`] took triangles off in.
    pub peeled_rings: u32,
}

/// Grow the smooth region of `mesh` about `seeds`: the triangles reached
/// from the ones nearest the seeds across edges where the surface turns by
/// no more than `crease`, in radians, with no surface fitted.
///
/// `index` is the hierarchy built over `mesh`. The mesh's triangles are
/// read as wound consistently; an edge whose two triangles run it the same
/// way is a fold, and growth does not cross it.
///
/// # Seams
///
/// Vertices within the confusion distance of each other are one, so two
/// patches meshed apart along a shared border, each with its own vertices
/// there, meet across it as though meshed together. Where the patches'
/// borders are divided differently, so that one's vertex lies along the
/// other's edge, they do not meet there. An edge used by one triangle, or
/// by more than two, is a border of the mesh, and growth stops at it.
///
/// # Smoothing
///
/// A scan's noise turns neighbouring triangles' normals by about the noise
/// over the triangles' size, often past any bound that would still stop at
/// a shallow crease. So each triangle's normal is smoothed before it is
/// compared, by bilateral filtering: every round of
/// [`GrowOptions::smoothing_rounds`] averages each triangle's normal with
/// its neighbours' across its edges, each weighing its area times
/// `exp(-(θ/s)²)`, where `θ` is the angle between the two normals the last
/// round gave and `s` is the range. The range is twice the median angle
/// between neighbouring triangles' own normals over the mesh, a measure of
/// the noise that the mesh's creases, being few, do not move, and never
/// more than `crease`. A turn of a few times the noise's is then averaged
/// out, while a crease turning three ranges or more keeps its full turn,
/// its sides weighing each other at `e⁻⁹` or less. A mesh whose median turn is
/// within the angular tolerance, as a noise-free mesh of flat faces is, is
/// not smoothed.
///
/// # Growth
///
/// The region starts at the triangle nearest each seed, all at once: seeds
/// on separate smooth parts give their union. A neighbour across an edge
/// joins where the two triangles' smoothed normals turn by no more than
/// `crease`. That is tangent continuity: across a seam where the surface
/// stays tangent but its curvature jumps, as between a plane and a round
/// that leaves it, the normals turn only by the curvature over a
/// triangle's size and growth runs on. Curvature continuity is not asked
/// for. A crease holds where its turn is past `crease` by clearly more
/// than the smoothed normals' noise. A round turns each edge by its angle
/// over the edges across it, and is crossed where that is no more than
/// `crease`, so whether a round is a crease depends on how finely it is
/// meshed, as it does for any bound on the angle across an edge.
///
/// A triangle without area never joins. Growth stops at
/// [`GrowOptions::max_triangles`], the seeds' own triangles kept past it.
///
/// # Peeling
///
/// A scanner rounds a sharp edge, and where the round's first rows turn by
/// less than `crease` the region climbs onto them. [`GrowOptions::peel`]
/// takes them back, pass by pass. The rim is the region's triangles with a
/// corner on its border, so both halves of a split quad there are on it.
/// Before the first pass, the turns between the smoothed normals of every
/// two triangles off the rim that share a corner give the threshold: the
/// cutoff's threshold over them, and no less than the largest of them, so
/// that the noise's rare large turns inside the region do not set it. Each
/// pass then takes every rim triangle whose smoothed normal turns past the
/// threshold from that of a triangle off the rim sharing a corner with it,
/// or, where it shares a corner with none, as in a corner of the region,
/// from that of any triangle of the region; passes end when one takes
/// nothing, or would take every triangle left.
///
/// Of the other options, [`GrowOptions::cutoff`] only sets peeling's
/// threshold; the normal limit, the crease limit, what is held fixed and
/// the side are a fit's and do not apply.
///
/// # Determinism
///
/// Everything runs on the calling thread in an order fixed by the seeds and
/// the triangle indices: the same mesh and seeds give the same region on
/// any thread count.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) with no seeds, a
/// seed that is not finite, a `crease` not above 0 and below π, a cutoff
/// whose number is not positive, peeling with [`Cutoff::All`], which has
/// no threshold, or [`Cutoff::Distance`], whose threshold is a length and
/// not an angle, a triangle cap of zero, or
/// a triangle naming a vertex the mesh lacks;
/// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) where the
/// hierarchy holds a different number of triangles from the mesh;
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// where the mesh has no triangle near a seed, or the triangle nearest one
/// has no area.
pub fn grow_smooth(
    mesh: &Triangulation,
    index: &Bvh,
    seeds: &[Point],
    crease: f64,
    options: &GrowOptions,
    tol: Tolerances,
) -> OgeomResult<GrownRegion> {
    checked(mesh, index, seeds, crease, options)?;
    let surface = Smooth::of(mesh, crease, options.smoothing_rounds, tol);
    let mut hits = Vec::with_capacity(seeds.len());
    for (i, seed) in seeds.iter().enumerate() {
        let Some(hit) = index.nearest(*seed, f64::INFINITY) else {
            ogeom_bail!(Construction, "no triangle of the mesh lies near seed {i}");
        };
        if surface.area[hit.triangle as usize] <= 0.0 {
            ogeom_bail!(
                Construction,
                "the triangle nearest seed {i} has no area, so no normal to grow from"
            );
        }
        hits.push(hit.triangle);
    }
    let mut in_region = surface.spread(&hits, crease.cos(), options.max_triangles);
    let peeled_rings = if options.peel {
        surface.peel(&mut in_region, options.cutoff, tol.angular())
    } else {
        0
    };
    let triangles: Vec<u32> = (0..mesh.triangles.len())
        .filter(|&t| in_region[t])
        .map(as_index)
        .collect();
    let mut lowest = vec![NONE; surface.weld_count];
    for &t in &triangles {
        for v in mesh.triangles[t as usize] {
            let w = surface.weld[v as usize] as usize;
            lowest[w] = lowest[w].min(v);
        }
    }
    let mut vertices: Vec<u32> = lowest.into_iter().filter(|&v| v != NONE).collect();
    vertices.sort_unstable();
    let points = vertices
        .iter()
        .map(|&v| mesh.positions[v as usize])
        .collect();
    Ok(GrownRegion {
        triangles,
        vertices,
        points,
        seeds: hits,
        peeled_rings,
    })
}

fn checked(
    mesh: &Triangulation,
    index: &Bvh,
    seeds: &[Point],
    crease: f64,
    options: &GrowOptions,
) -> OgeomResult<()> {
    if seeds.is_empty() {
        ogeom_bail!(Domain, "growing a region needs at least one seed");
    }
    if let Some(i) = seeds.iter().position(|p| !p.is_finite()) {
        ogeom_bail!(Domain, "seed {i} is not finite");
    }
    if !(crease > 0.0 && crease < core::f64::consts::PI) {
        ogeom_bail!(
            Domain,
            "the crease angle {crease} must lie above 0 and below pi"
        );
    }
    options.cutoff.check()?;
    if options.peel && matches!(options.cutoff, Cutoff::All | Cutoff::Distance(_)) {
        ogeom_bail!(
            Domain,
            "peeling a smooth region needs a threshold on angles, and {:?} gives none",
            options.cutoff
        );
    }
    if options.max_triangles == Some(0) {
        ogeom_bail!(Domain, "a region capped at no triangles");
    }
    check_indices(mesh)?;
    if index.triangle_count() != mesh.triangles.len() {
        ogeom_bail!(
            Dimension,
            "the hierarchy holds {} triangles and the mesh {}",
            index.triangle_count(),
            mesh.triangles.len()
        );
    }
    Ok(())
}

/// A use of an edge: the triangle, which of its edges, and whether it runs
/// from the lower welded vertex.
type Use = (u32, u8, bool);

/// A mesh read for smooth growth: its vertices welded, each triangle's
/// neighbours across its edges, its own normal and area, and its smoothed
/// normal.
struct Smooth {
    /// The welded vertex each mesh vertex is.
    weld: Vec<u32>,
    weld_count: usize,
    /// Each triangle's corners, welded.
    corners: Vec<[u32; 3]>,
    /// The triangle across each edge `(k, k + 1)` where exactly two use it,
    /// running it opposite ways.
    across: Vec<[u32; 3]>,
    /// Unit, or zero for a triangle without area.
    normal: Vec<Vector>,
    area: Vec<f64>,
    smooth: Vec<Vector>,
}

impl Smooth {
    fn of(mesh: &Triangulation, crease: f64, rounds: u32, tol: Tolerances) -> Self {
        let (kept, weld) = weld_points(&mesh.positions, tol.confusion());
        let weld_count = kept.len();
        let count = mesh.triangles.len();
        let mut normal = Vec::with_capacity(count);
        let mut area = Vec::with_capacity(count);
        for t in &mesh.triangles {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            let doubled = (b - a).cross(c - a);
            let length = doubled.magnitude();
            let collapsed = weld[t[0] as usize] == weld[t[1] as usize]
                || weld[t[1] as usize] == weld[t[2] as usize]
                || weld[t[2] as usize] == weld[t[0] as usize];
            if !collapsed && length.is_finite() && length > 0.0 {
                normal.push(doubled / length);
                area.push(0.5 * length);
            } else {
                normal.push(Vector::ZERO);
                area.push(0.0);
            }
        }
        // Each welded edge's uses, as (triangle, edge, forward).
        let mut uses: FastMap<(u32, u32), Vec<Use>> =
            FastMap::with_capacity_and_hasher(3 * count, Default::default());
        for (t, tri) in mesh.triangles.iter().enumerate() {
            if area[t] <= 0.0 {
                continue;
            }
            for k in 0..3_u8 {
                let a = weld[tri[k as usize] as usize];
                let b = weld[tri[(k as usize + 1) % 3] as usize];
                let key = if a < b { (a, b) } else { (b, a) };
                uses.entry(key).or_default().push((as_index(t), k, a < b));
            }
        }
        let mut across = vec![[NONE; 3]; count];
        for list in uses.values() {
            if let [(t, i, f), (u, j, g)] = list[..]
                && f != g
            {
                across[t as usize][i as usize] = u;
                across[u as usize][j as usize] = t;
            }
        }
        let corners = mesh
            .triangles
            .iter()
            .map(|t| t.map(|v| weld[v as usize]))
            .collect();
        let mut surface = Self {
            weld,
            weld_count,
            corners,
            across,
            normal,
            area,
            smooth: Vec::new(),
        };
        surface.smooth = surface.smoothed(crease, rounds, tol.angular());
        surface
    }

    /// Each triangle's normal averaged with its neighbours' across its
    /// edges, `rounds` times over, every round from the last round's
    /// normals: a neighbour weighs its area times `exp(-(θ/s)²)`, `θ` the
    /// angle between the two triangles' last normals and `s` the range,
    /// twice the median angle between neighbours' own normals over the
    /// mesh and at most the crease. A mesh whose median turn is within the
    /// angular tolerance has no noise to smooth.
    fn smoothed(&self, crease: f64, rounds: u32, angular: f64) -> Vec<Vector> {
        let mut smooth = self.normal.clone();
        let mut turns: Vec<f64> = (0..self.across.len())
            .flat_map(|t| {
                self.across[t]
                    .iter()
                    .filter(move |&&n| n != NONE && (n as usize) > t)
                    .map(move |&n| self.normal[t].dot(self.normal[n as usize]))
            })
            .map(|cos| cos.clamp(-1.0, 1.0).acos())
            .collect();
        if turns.is_empty() {
            return smooth;
        }
        let middle = turns.len() / 2;
        let (_, median, _) = turns.select_nth_unstable_by(middle, f64::total_cmp);
        let range = (2.0 * *median).min(crease);
        if range <= angular {
            return smooth;
        }
        for _ in 0..rounds {
            let last = smooth;
            smooth = (0..last.len())
                .map(|t| {
                    let mut sum = last[t] * self.area[t];
                    for &n in &self.across[t] {
                        if n != NONE {
                            let n = n as usize;
                            let angle = last[t].dot(last[n]).clamp(-1.0, 1.0).acos();
                            let weight = (-(angle / range).powi(2)).exp();
                            sum += last[n] * (self.area[n] * weight);
                        }
                    }
                    let m = sum.magnitude();
                    if m > 0.0 { sum / m } else { last[t] }
                })
                .collect();
        }
        smooth
    }

    /// The triangles reached from the starts across edges whose smoothed
    /// normals turn by no more than the crease, breadth first, until there
    /// are no more or the cap is reached.
    fn spread(&self, starts: &[u32], cos_crease: f64, cap: Option<usize>) -> Vec<bool> {
        let mut in_region = vec![false; self.area.len()];
        let mut queue = VecDeque::new();
        let mut size = 0;
        for &t in starts {
            if !in_region[t as usize] {
                in_region[t as usize] = true;
                size += 1;
                queue.push_back(t);
            }
        }
        while let Some(t) = queue.pop_front() {
            for &n in &self.across[t as usize] {
                if cap.is_some_and(|most| size >= most) {
                    return in_region;
                }
                if n == NONE || in_region[n as usize] {
                    continue;
                }
                if self.smooth[t as usize].dot(self.smooth[n as usize]) >= cos_crease {
                    in_region[n as usize] = true;
                    size += 1;
                    queue.push_back(n);
                }
            }
        }
        in_region
    }

    /// The angle between two triangles' smoothed normals, exact near zero.
    fn turn(&self, t: u32, n: u32) -> f64 {
        let (a, b) = (self.smooth[t as usize], self.smooth[n as usize]);
        a.cross(b).magnitude().atan2(a.dot(b))
    }

    /// The triangles about each welded vertex, as starts into one list.
    fn fans(&self) -> (Vec<usize>, Vec<u32>) {
        let mut start = vec![0_usize; self.weld_count + 1];
        for corners in &self.corners {
            for &w in corners {
                start[w as usize + 1] += 1;
            }
        }
        for w in 0..self.weld_count {
            start[w + 1] += start[w];
        }
        let mut fill = start.clone();
        let mut fan = vec![0_u32; start[self.weld_count]];
        for (t, corners) in self.corners.iter().enumerate() {
            for &w in corners {
                fan[fill[w as usize]] = as_index(t);
                fill[w as usize] += 1;
            }
        }
        (start, fan)
    }

    /// Take the rim off, pass by pass, as [`grow_smooth`] describes; the
    /// number of passes that took a triangle.
    fn peel(&self, in_region: &mut [bool], cutoff: Cutoff, angular: f64) -> u32 {
        let (start, fan) = self.fans();
        let about = |t: u32| {
            self.corners[t as usize].into_iter().flat_map(|w| {
                fan[start[w as usize]..start[w as usize + 1]]
                    .iter()
                    .copied()
            })
        };
        let mut passes = 0;
        let mut threshold = None;
        loop {
            let mut border = vec![false; self.weld_count];
            for t in 0..in_region.len() {
                if !in_region[t] {
                    continue;
                }
                for (k, &n) in self.across[t].iter().enumerate() {
                    if n == NONE || !in_region[n as usize] {
                        border[self.corners[t][k] as usize] = true;
                        border[self.corners[t][(k + 1) % 3] as usize] = true;
                    }
                }
            }
            let on_rim: Vec<bool> = (0..in_region.len())
                .map(|t| in_region[t] && self.corners[t].iter().any(|&w| border[w as usize]))
                .collect();
            let off_rim = |n: u32| in_region[n as usize] && !on_rim[n as usize];
            let inside: Vec<f64> = (0..in_region.len())
                .map(as_index)
                .filter(|&t| off_rim(t))
                .flat_map(|t| about(t).filter(move |&n| n > t).map(move |n| (t, n)))
                .filter(|&(_, n)| off_rim(n))
                .map(|(t, n)| self.turn(t, n))
                .collect();
            if inside.len() < 2 {
                break;
            }
            let threshold = *threshold.get_or_insert_with(|| {
                let most = inside.iter().fold(angular, |m, &a| m.max(a));
                cutoff.threshold(&inside).max(most)
            });
            let rim: Vec<u32> = (0..in_region.len())
                .map(as_index)
                .filter(|&t| on_rim[t as usize])
                .collect();
            let peeled: Vec<u32> = rim
                .iter()
                .copied()
                .filter(|&t| {
                    let inner = about(t).any(off_rim);
                    about(t).any(|n| {
                        n != t
                            && in_region[n as usize]
                            && (off_rim(n) || !inner)
                            && self.turn(t, n) > threshold
                    })
                })
                .collect();
            let left = in_region.iter().filter(|&&r| r).count();
            if peeled.is_empty() || peeled.len() == left {
                break;
            }
            for &t in &peeled {
                in_region[t as usize] = false;
            }
            passes += 1;
        }
        passes
    }
}
