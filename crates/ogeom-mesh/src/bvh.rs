//! A bounding-volume hierarchy over a triangle mesh: the nearest triangle to
//! a point, every triangle within a ball, and the first triangle along a ray.
//!
//! Built once, immutable after, and shared freely across threads. Every
//! answer is the answer brute force gives over the same triangles with the
//! same arithmetic, ties going to the lowest triangle index, so it depends on
//! neither the tree's shape, the order it is walked in, nor the thread count.
//!
//! # Why the answers equal brute force
//!
//! A query visits a node only if the node's box could hold a triangle that
//! passes the query's own test. That is decided with the triangle test's own
//! floating-point operations applied to the box's corners, and rounding is
//! monotone (a larger operand never rounds to a smaller result), so a box's
//! computed bound can never exclude a triangle inside it whose computed value
//! passes:
//!
//! - **nearest and ball**: the closest point on a triangle is clamped into the
//!   triangle's own box before its distance is taken, so each axis gap to it
//!   is at least the computed gap to any box that holds the triangle;
//! - **ray**: the ray is sheared onto the `z` axis and every vertex
//!   transformed by the same operations the box corners are, so the computed
//!   vertices of a triangle lie within the box's computed sheared bounds. A
//!   hit is decided by the exact signs of three orientations of the origin
//!   against those computed vertices, so a hit needs the origin within their
//!   hull, and the hit's depth is clamped into their depth range.
//!
//! # Watertight rays
//!
//! Two triangles that share an edge share its two vertices, so they compute
//! the same sheared coordinates for it, and the exact orientation of the ray
//! against it has opposite signs from the two sides. A ray through the edge
//! (or a vertex) is therefore inside at least one of them: it cannot slip
//! between triangles, however it is aimed.

use ogeom_core::{Exact, OgeomResult, Predicates, Sign, ogeom_bail};
use ogeom_math::{Aabb, Point, Vector};
use ogeom_topo::{Model, Shape, Triangulation};

use crate::discretize::Deflection;

/// What a query found on the mesh.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hit {
    /// The triangle, as its index in the mesh the hierarchy was built from.
    pub triangle: u32,
    /// The point on the triangle.
    pub point: Point,
    /// For a nearest query, the distance from the query point to
    /// [`Hit::point`]; for a ray, the distance along the ray.
    pub distance: f64,
    /// [`Hit::point`]'s weights on the triangle's three corners, in the
    /// triangle's own corner order. They sum to one.
    pub barycentric: [f64; 3],
}

type P = [f64; 3];

/// One node of the tree, its box tight around the triangles below it.
#[derive(Debug, Clone, Copy)]
struct Node {
    low: P,
    high: P,
    /// A leaf's first slot, or an interior node's right child. Its left
    /// child is the next node.
    first: u32,
    /// How many triangles a leaf holds; zero for an interior node.
    count: u32,
}

/// A bounding-volume hierarchy over a mesh's triangles.
///
/// Built with a surface-area heuristic over binned centroids. Holds its own
/// copy of the triangles' corners, so the mesh it was built from may be
/// dropped or changed; queries answer about the mesh as it was.
#[derive(Debug, Clone)]
pub struct Bvh {
    nodes: Vec<Node>,
    /// Each slot's triangle corners, slots in leaf order.
    corners: Vec<[P; 3]>,
    /// The triangle in each slot.
    order: Vec<u32>,
    /// The slot holding each triangle.
    slot: Vec<u32>,
}

/// The most triangles a leaf holds when splitting it would not pay.
const LEAF_MOST: usize = 8;

/// How many bins each axis's centroids are counted into to price a split.
const BINS: usize = 16;

/// How many queries a batch hands a worker at a time.
const BATCH: usize = 256;

impl Bvh {
    /// Build the hierarchy over a mesh's triangles.
    ///
    /// Degenerate triangles are kept: a nearest query can still land on
    /// one, and a ray through one finds it only where it is a proper
    /// triangle as the ray sees it.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Range`](ogeom_core::OgeomError::Range) if a triangle
    /// names a vertex the mesh does not have, and
    /// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if a
    /// triangle's corner is not finite.
    pub fn build(mesh: &Triangulation) -> OgeomResult<Self> {
        let mut corners = Vec::with_capacity(mesh.triangles.len());
        for triangle in &mesh.triangles {
            let mut three = [[0.0; 3]; 3];
            for (corner, &k) in three.iter_mut().zip(triangle) {
                let Some(p) = mesh.positions.get(k as usize) else {
                    ogeom_bail!(
                        Range,
                        "a triangle names vertex {k} of {}",
                        mesh.positions.len()
                    );
                };
                if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
                    ogeom_bail!(Numeric, "a triangle's corner {k} is not finite");
                }
                *corner = [p.x, p.y, p.z];
            }
            corners.push(three);
        }
        if u32::try_from(corners.len()).is_err() {
            ogeom_bail!(Range, "a mesh of more than four billion triangles");
        }
        Ok(Self::over(&corners))
    }

    /// Build the hierarchy over a shape's triangulation at a deflection.
    ///
    /// The triangle indices are those of
    /// [`triangulate`](crate::triangulate::triangulate) for the same
    /// arguments, which is deterministic, and [`Bvh::triangle`] gives each
    /// triangle's corners without it.
    ///
    /// # Errors
    ///
    /// As [`triangulate`](crate::triangulate::triangulate) and
    /// [`Bvh::build`].
    pub fn of_shape(
        model: &Model,
        shape: &Shape,
        deflection: Deflection,
        tol: ogeom_core::Tolerances,
    ) -> OgeomResult<Self> {
        Self::build(&crate::triangulate::triangulate(
            model, shape, deflection, tol,
        )?)
    }

    fn over(corners: &[[P; 3]]) -> Self {
        let mut order: Vec<u32> = (0..corners.len())
            .map(|i| u32::try_from(i).unwrap_or(u32::MAX))
            .collect();
        let boxes: Vec<(P, P)> = corners.iter().map(box_of).collect();
        let nodes = build_nodes(&boxes, &mut order);
        let mut slot = vec![0_u32; order.len()];
        for (at, &t) in order.iter().enumerate() {
            slot[t as usize] = u32::try_from(at).unwrap_or(u32::MAX);
        }
        Self {
            nodes,
            corners: order.iter().map(|&t| corners[t as usize]).collect(),
            order,
            slot,
        }
    }

    /// How many triangles the hierarchy holds.
    #[must_use]
    pub fn triangle_count(&self) -> usize {
        self.order.len()
    }

    /// Whether the hierarchy holds no triangle.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.order.is_empty()
    }

    /// A triangle's three corners, or `None` past the last triangle.
    #[must_use]
    pub fn triangle(&self, triangle: u32) -> Option<[Point; 3]> {
        let slot = *self.slot.get(triangle as usize)?;
        Some(self.corners[slot as usize].map(|[x, y, z]| Point::new(x, y, z)))
    }

    /// The bounds of every triangle.
    #[must_use]
    pub fn bounds(&self) -> Aabb {
        self.nodes.first().map_or(Aabb::EMPTY, |root| {
            Aabb::of_corners(point(root.low), point(root.high))
        })
    }

    /// The point on the mesh nearest `p`, no farther than `max`; the lowest
    /// triangle index among equally near ones. `max` may be infinite.
    #[must_use]
    pub fn nearest(&self, p: Point, max: f64) -> Option<Hit> {
        self.nearest_with(p, max, &mut Vec::new())
    }

    /// [`Bvh::nearest`] for each point, in parallel, answers in point order.
    #[must_use]
    pub fn nearest_many(&self, points: &[Point], max: f64) -> Vec<Option<Hit>> {
        batched(points, |&p| p, |&p, stack| self.nearest_with(p, max, stack))
    }

    fn nearest_with(&self, p: Point, max: f64, stack: &mut Vec<(u32, f64)>) -> Option<Hit> {
        let p = [p.x, p.y, p.z];
        let reach = max * max;
        if self.nodes.is_empty() || !p.iter().all(|c| c.is_finite()) || max.is_nan() || max < 0.0 {
            return None;
        }
        // The best so far as (squared distance, triangle, slot); `reach`
        // stands in until a triangle is found.
        let mut best: Option<(f64, u32, usize)> = None;
        let mut bound = reach;
        stack.clear();
        let root = box_gap(&self.nodes[0], p);
        if root > bound {
            return None;
        }
        stack.push((0, root));
        while let Some((at, near)) = stack.pop() {
            if near > bound {
                continue;
            }
            let node = &self.nodes[at as usize];
            if node.count > 0 {
                let first = node.first as usize;
                for slot in first..first + node.count as usize {
                    let triangle = self.order[slot];
                    let (_, _, d2) = closest(p, &self.corners[slot]);
                    let better = match best {
                        None => d2 <= bound,
                        Some((held, index, _)) => d2 < held || (d2 == held && triangle < index),
                    };
                    if better {
                        best = Some((d2, triangle, slot));
                        bound = d2;
                    }
                }
                continue;
            }
            let (left, right) = (at + 1, node.first);
            let gl = box_gap(&self.nodes[left as usize], p);
            let gr = box_gap(&self.nodes[right as usize], p);
            // The nearer child is popped first.
            let (near, far) = if gl <= gr {
                ((left, gl), (right, gr))
            } else {
                ((right, gr), (left, gl))
            };
            if far.1 <= bound {
                stack.push(far);
            }
            if near.1 <= bound {
                stack.push(near);
            }
        }
        let (_, triangle, slot) = best?;
        let (q, barycentric, d2) = closest(p, &self.corners[slot]);
        Some(Hit {
            triangle,
            point: point(q),
            distance: d2.sqrt(),
            barycentric,
        })
    }

    /// Every triangle with a point within `radius` of `centre`, in index
    /// order.
    #[must_use]
    pub fn within_ball(&self, centre: Point, radius: f64) -> Vec<u32> {
        let c = [centre.x, centre.y, centre.z];
        let reach = radius * radius;
        let mut found = Vec::new();
        if self.nodes.is_empty()
            || !c.iter().all(|v| v.is_finite())
            || radius.is_nan()
            || radius < 0.0
        {
            return found;
        }
        let mut stack = vec![0_u32];
        while let Some(at) = stack.pop() {
            let node = &self.nodes[at as usize];
            if box_gap(node, c) > reach {
                continue;
            }
            if node.count > 0 {
                let first = node.first as usize;
                for slot in first..first + node.count as usize {
                    if closest(c, &self.corners[slot]).2 <= reach {
                        found.push(self.order[slot]);
                    }
                }
            } else {
                stack.push(node.first);
                stack.push(at + 1);
            }
        }
        found.sort_unstable();
        found
    }

    /// Every triangle whose bounding box passes `test`, in index order.
    ///
    /// `test` prunes the walk: a node's box that fails it is not entered,
    /// so `test` must pass every box that contains a box it passes (a
    /// region overlap does, a containment test does not).
    #[must_use]
    pub fn triangles_where(&self, test: impl Fn(&Aabb) -> bool) -> Vec<u32> {
        let mut found = Vec::new();
        if self.nodes.is_empty() {
            return found;
        }
        let mut stack = vec![0_u32];
        while let Some(at) = stack.pop() {
            let node = &self.nodes[at as usize];
            if !test(&Aabb::of_corners(point(node.low), point(node.high))) {
                continue;
            }
            if node.count > 0 {
                let first = node.first as usize;
                for slot in first..first + node.count as usize {
                    let (low, high) = box_of(&self.corners[slot]);
                    if test(&Aabb::of_corners(point(low), point(high))) {
                        found.push(self.order[slot]);
                    }
                }
            } else {
                stack.push(node.first);
                stack.push(at + 1);
            }
        }
        found.sort_unstable();
        found
    }

    /// The first triangle along a ray, no farther than `max` along it, and
    /// the lowest triangle index among equally far ones.
    ///
    /// `direction` need not be unit length; [`Hit::distance`] is measured
    /// along it in model units. A triangle faces the ray when the ray meets
    /// the side its counter-clockwise winding points out of; with
    /// `both_sides` false only those count. A hit at the origin itself
    /// counts, at distance zero. A ray in a triangle's plane does not hit
    /// it, having no area to cross there; in a closed mesh the triangles
    /// across the edges it passes through are hit instead.
    #[must_use]
    pub fn ray(&self, origin: Point, direction: Vector, max: f64, both_sides: bool) -> Option<Hit> {
        let shear = Shear::new(origin, direction)?;
        self.ray_with(&shear, max, both_sides, &mut Vec::new())
    }

    /// [`Bvh::ray`] for each origin and direction, in parallel, answers in
    /// input order.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) if there
    /// are not as many directions as origins.
    pub fn rays(
        &self,
        origins: &[Point],
        directions: &[Vector],
        max: f64,
        both_sides: bool,
    ) -> OgeomResult<Vec<Option<Hit>>> {
        if origins.len() != directions.len() {
            ogeom_bail!(
                Dimension,
                "{} ray origins against {} directions",
                origins.len(),
                directions.len()
            );
        }
        let pairs: Vec<(Point, Vector)> = origins
            .iter()
            .copied()
            .zip(directions.iter().copied())
            .collect();
        Ok(batched(
            &pairs,
            |&(o, _)| o,
            |&(o, d), stack| {
                let shear = Shear::new(o, d)?;
                self.ray_with(&shear, max, both_sides, stack)
            },
        ))
    }

    /// Every triangle along a ray within `max`, nearest first, ties in
    /// index order. Each is the hit [`Bvh::ray`] reports when it is first.
    #[must_use]
    pub fn hits_along(
        &self,
        origin: Point,
        direction: Vector,
        max: f64,
        both_sides: bool,
    ) -> Vec<Hit> {
        let mut hits = Vec::new();
        let Some(shear) = Shear::new(origin, direction) else {
            return hits;
        };
        if self.nodes.is_empty() || max.is_nan() || max < 0.0 {
            return hits;
        }
        let mut stack = vec![0_u32];
        while let Some(at) = stack.pop() {
            let node = &self.nodes[at as usize];
            if shear.box_depth(node).is_none_or(|near| near > max) {
                continue;
            }
            if node.count > 0 {
                let first = node.first as usize;
                for slot in first..first + node.count as usize {
                    if let Some(hit) = shear.triangle(&self.corners[slot], max, both_sides) {
                        hits.push(Hit {
                            triangle: self.order[slot],
                            ..hit
                        });
                    }
                }
            } else {
                stack.push(node.first);
                stack.push(at + 1);
            }
        }
        hits.sort_by(|a, b| {
            a.distance
                .total_cmp(&b.distance)
                .then(a.triangle.cmp(&b.triangle))
        });
        hits
    }

    fn ray_with(
        &self,
        shear: &Shear,
        max: f64,
        both_sides: bool,
        stack: &mut Vec<(u32, f64)>,
    ) -> Option<Hit> {
        if self.nodes.is_empty() || max.is_nan() || max < 0.0 {
            return None;
        }
        let mut best: Option<Hit> = None;
        let mut bound = max;
        stack.clear();
        let near = shear.box_depth(&self.nodes[0])?;
        if near > bound {
            return None;
        }
        stack.push((0, near));
        while let Some((at, near)) = stack.pop() {
            if near > bound {
                continue;
            }
            let node = &self.nodes[at as usize];
            if node.count > 0 {
                let first = node.first as usize;
                for slot in first..first + node.count as usize {
                    let Some(hit) = shear.triangle(&self.corners[slot], bound, both_sides) else {
                        continue;
                    };
                    let triangle = self.order[slot];
                    let better = best.is_none_or(|held| {
                        hit.distance < held.distance
                            || (hit.distance == held.distance && triangle < held.triangle)
                    });
                    if better {
                        bound = hit.distance;
                        best = Some(Hit { triangle, ..hit });
                    }
                }
                continue;
            }
            let (left, right) = (at + 1, node.first);
            let dl = shear.box_depth(&self.nodes[left as usize]);
            let dr = shear.box_depth(&self.nodes[right as usize]);
            let (first, second) = match (dl, dr) {
                (Some(l), Some(r)) if r < l => (Some((right, r)), Some((left, l))),
                (l, r) => (l.map(|l| (left, l)), r.map(|r| (right, r))),
            };
            for (index, depth) in [second, first].into_iter().flatten() {
                if depth <= bound {
                    stack.push((index, depth));
                }
            }
        }
        best
    }
}

/// Answer `queries` on the kernel's workers, answers in query order.
///
/// The queries are taken in the order of their places along a Morton curve
/// over their bounds, so that each worker's run of them walks the same few
/// nodes and triangles while they are in cache. Each answer depends on its
/// query alone, so the order changes only the time.
pub(crate) fn batched<T: Sync, R: Send + Copy + Default>(
    queries: &[T],
    place: impl Fn(&T) -> Point,
    answer: impl Fn(&T, &mut Vec<(u32, f64)>) -> R + Sync,
) -> Vec<R> {
    let places: Vec<Point> = queries.iter().map(place).collect();
    let bounds = places
        .iter()
        .filter(|p| p.x.is_finite() && p.y.is_finite() && p.z.is_finite())
        .fold(Aabb::EMPTY, |acc, &p| acc.with_point(p));
    let mut order: Vec<(u64, usize)> = places
        .iter()
        .enumerate()
        .map(|(i, &p)| (morton(&bounds, p), i))
        .collect();
    order.sort_unstable();
    let chunks: Vec<&[(u64, usize)]> = order.chunks(BATCH).collect();
    let answered = ogeom_core::parallel::map_ordered(&chunks, |_, chunk| {
        let mut stack = Vec::new();
        chunk
            .iter()
            .map(|&(_, i)| answer(&queries[i], &mut stack))
            .collect::<Vec<R>>()
    });
    let mut out = vec![R::default(); queries.len()];
    for (&(_, i), r) in order.iter().zip(answered.into_iter().flatten()) {
        out[i] = r;
    }
    out
}

/// A point's place along a Morton curve over `bounds`: 21 bits per axis,
/// interleaved. Off the bounds (or not finite), the nearest cell.
fn morton(bounds: &Aabb, p: Point) -> u64 {
    let (Some(low), Some(high)) = (bounds.low(), bounds.high()) else {
        return 0;
    };
    let cell = |v: f64, lo: f64, hi: f64| -> u64 {
        let span = hi - lo;
        let unit = if span > 0.0 { (v - lo) / span } else { 0.0 };
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let at = (unit.clamp(0.0, 1.0) * 2_097_151.0) as u64;
        // Spread 21 bits three apart.
        let mut x = at & 0x1f_ffff;
        x = (x | x << 32) & 0x1f_0000_0000_ffff;
        x = (x | x << 16) & 0x1f_0000_ff00_00ff;
        x = (x | x << 8) & 0x100f_00f0_0f00_f00f;
        x = (x | x << 4) & 0x10c3_0c30_c30c_30c3;
        x = (x | x << 2) & 0x1249_2492_4924_9249;
        x
    };
    cell(p.x, low.x, high.x) | cell(p.y, low.y, high.y) << 1 | cell(p.z, low.z, high.z) << 2
}

pub(crate) fn point(p: P) -> Point {
    Point::new(p[0], p[1], p[2])
}

fn box_of(t: &[P; 3]) -> (P, P) {
    let mut low = t[0];
    let mut high = t[0];
    for c in &t[1..] {
        for axis in 0..3 {
            low[axis] = low[axis].min(c[axis]);
            high[axis] = high[axis].max(c[axis]);
        }
    }
    (low, high)
}

/// The squared distance from `p` to a box, summed over the axes in order:
/// the same operations [`closest`] takes its distance with.
pub(crate) fn gap(low: &P, high: &P, p: P) -> f64 {
    let mut sum = 0.0;
    for axis in 0..3 {
        let g = if p[axis] < low[axis] {
            low[axis] - p[axis]
        } else if p[axis] > high[axis] {
            p[axis] - high[axis]
        } else {
            0.0
        };
        sum += g * g;
    }
    sum
}

fn box_gap(node: &Node, p: P) -> f64 {
    gap(&node.low, &node.high, p)
}

/// The squared distance between two points, summed over the axes in order.
pub(crate) fn distance2(a: P, b: P) -> f64 {
    let mut sum = 0.0;
    for axis in 0..3 {
        let d = a[axis] - b[axis];
        sum += d * d;
    }
    sum
}

fn sub(a: P, b: P) -> P {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: P, b: P) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn along(a: P, d: P, s: f64) -> P {
    [a[0] + d[0] * s, a[1] + d[1] * s, a[2] + d[2] * s]
}

/// The point on a triangle nearest `p`, its barycentric weights, and its
/// squared distance to `p`.
///
/// By the Voronoi regions of the corners, edges and face. A triangle too
/// degenerate for its face region is answered by its nearest edge. The point
/// is clamped into the triangle's box, which rounding can step outside of,
/// so that no box holding the triangle is farther from `p` than the answer.
fn closest(p: P, t: &[P; 3]) -> (P, [f64; 3], f64) {
    let (q, w) = closest_unclamped(p, t);
    let (low, high) = box_of(t);
    let q = [0, 1, 2].map(|axis| q[axis].clamp(low[axis], high[axis]));
    (q, w, distance2(q, p))
}

fn closest_unclamped(p: P, [a, b, c]: &[P; 3]) -> (P, [f64; 3]) {
    let (a, b, c) = (*a, *b, *c);
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return (a, [1.0, 0.0, 0.0]);
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return (b, [0.0, 1.0, 0.0]);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 && d1 - d3 > 0.0 {
        let v = d1 / (d1 - d3);
        return (along(a, ab, v), [1.0 - v, v, 0.0]);
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return (c, [0.0, 0.0, 1.0]);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 && d2 - d6 > 0.0 {
        let w = d2 / (d2 - d6);
        return (along(a, ac, w), [1.0 - w, 0.0, w]);
    }
    let va = d3 * d6 - d5 * d4;
    let (e, f) = (d4 - d3, d5 - d6);
    if va <= 0.0 && e >= 0.0 && f >= 0.0 && e + f > 0.0 {
        let w = e / (e + f);
        return (along(b, sub(c, b), w), [0.0, 1.0 - w, w]);
    }
    let sum = va + vb + vc;
    if sum > 0.0 && sum.is_finite() {
        let v = vb / sum;
        let w = vc / sum;
        if v >= 0.0 && w >= 0.0 && v + w <= 1.0 {
            return (along(along(a, ab, v), ac, w), [1.0 - v - w, v, w]);
        }
    }
    nearest_edge(p, a, b, c)
}

/// The nearest point on a triangle's three edges, the first edge winning a
/// tie: the answer for a triangle with no face region of its own.
fn nearest_edge(p: P, a: P, b: P, c: P) -> (P, [f64; 3]) {
    let on = |from: P, to: P| -> (P, f64) {
        let d = sub(to, from);
        let length2 = dot(d, d);
        let s = if length2 > 0.0 {
            (dot(sub(p, from), d) / length2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        (along(from, d, s), s)
    };
    let (qab, sab) = on(a, b);
    let (qbc, sbc) = on(b, c);
    let (qca, sca) = on(c, a);
    let mut best = (qab, [1.0 - sab, sab, 0.0]);
    let mut held = distance2(qab, p);
    for (q, w) in [(qbc, [0.0, 1.0 - sbc, sbc]), (qca, [sca, 0.0, 1.0 - sca])] {
        let d = distance2(q, p);
        if d < held {
            held = d;
            best = (q, w);
        }
    }
    best
}

/// A ray sheared onto the `z` axis: its origin moved to zero, its dominant
/// axis made `z`, and the other two leaned so the ray runs straight up it.
/// The two lateral axes are swapped when the ray runs down its dominant
/// axis, so that a counter-clockwise triangle keeps its sense.
#[derive(Debug, Clone, Copy)]
struct Shear {
    origin: P,
    axes: [usize; 3],
    lean: [f64; 2],
    scale: f64,
}

impl Shear {
    fn new(origin: Point, direction: Vector) -> Option<Self> {
        let length = direction.magnitude();
        let o = [origin.x, origin.y, origin.z];
        if !(length > 0.0 && length.is_finite()) || !o.iter().all(|v| v.is_finite()) {
            return None;
        }
        let d = [
            direction.x / length,
            direction.y / length,
            direction.z / length,
        ];
        let mut kz = 0;
        for axis in 1..3 {
            if d[axis].abs() > d[kz].abs() {
                kz = axis;
            }
        }
        let mut kx = (kz + 1) % 3;
        let mut ky = (kx + 1) % 3;
        if d[kz] < 0.0 {
            core::mem::swap(&mut kx, &mut ky);
        }
        Some(Self {
            origin: o,
            axes: [kx, ky, kz],
            lean: [d[kx] / d[kz], d[ky] / d[kz]],
            scale: 1.0 / d[kz],
        })
    }

    /// A corner in sheared coordinates; `z` is the distance along the ray.
    fn transform(&self, p: P) -> P {
        let [kx, ky, kz] = self.axes;
        let dz = p[kz] - self.origin[kz];
        [
            (p[kx] - self.origin[kx]) - self.lean[0] * dz,
            (p[ky] - self.origin[ky]) - self.lean[1] * dz,
            self.scale * dz,
        ]
    }

    /// The nearest depth a box could hold a hit at, or `None` where the ray
    /// passes it by or it lies wholly behind the origin.
    ///
    /// The bounds of each sheared coordinate over the box, taken with
    /// [`Shear::transform`]'s own operations at the corners that make it
    /// least and greatest.
    fn box_depth(&self, node: &Node) -> Option<f64> {
        let [kx, ky, kz] = self.axes;
        let z0 = node.low[kz] - self.origin[kz];
        let z1 = node.high[kz] - self.origin[kz];
        let lateral = |k: usize, lean: f64| -> bool {
            let (s0, s1) = (lean * z0, lean * z1);
            let least = (node.low[k] - self.origin[k]) - s0.max(s1);
            let most = (node.high[k] - self.origin[k]) - s0.min(s1);
            least <= 0.0 && most >= 0.0
        };
        if !lateral(kx, self.lean[0]) || !lateral(ky, self.lean[1]) {
            return None;
        }
        let (d0, d1) = (self.scale * z0, self.scale * z1);
        let (near, far) = if d0 <= d1 { (d0, d1) } else { (d1, d0) };
        (far >= 0.0).then_some(near.max(0.0))
    }

    /// The ray's hit on one triangle within `max`, its triangle index left
    /// for the caller.
    fn triangle(&self, corners: &[P; 3], max: f64, both_sides: bool) -> Option<Hit> {
        let [a, b, c] = corners.map(|p| self.transform(p));
        let origin = [0.0, 0.0];
        let signs = [
            Exact::orient2d([c[0], c[1]], [b[0], b[1]], origin),
            Exact::orient2d([a[0], a[1]], [c[0], c[1]], origin),
            Exact::orient2d([b[0], b[1]], [a[0], a[1]], origin),
        ];
        let any = |s: Sign| signs.contains(&s);
        if (any(Sign::Negative) && any(Sign::Positive))
            || !(any(Sign::Negative) || any(Sign::Positive))
        {
            return None;
        }
        let facing = any(Sign::Positive);
        if !both_sides && !facing {
            return None;
        }
        // The weights in floating point, each held to its exact sign: one
        // the rounding put on the wrong side is too small to matter, and
        // zero keeps the three of one sign.
        let raw = [
            c[0] * b[1] - c[1] * b[0],
            a[0] * c[1] - a[1] * c[0],
            b[0] * a[1] - b[1] * a[0],
        ];
        let weights: [f64; 3] = [0, 1, 2].map(|i| match signs[i] {
            Sign::Zero => 0.0,
            Sign::Positive if raw[i] > 0.0 => raw[i],
            Sign::Negative if raw[i] < 0.0 => raw[i],
            _ => 0.0,
        });
        let det = weights[0] + weights[1] + weights[2];
        let barycentric = if det == 0.0 || !det.is_finite() {
            [1.0 / 3.0; 3]
        } else {
            weights.map(|w| w / det)
        };
        let depth = barycentric[0] * a[2] + barycentric[1] * b[2] + barycentric[2] * c[2];
        let least = a[2].min(b[2]).min(c[2]);
        let most = a[2].max(b[2]).max(c[2]);
        let depth = depth.clamp(least, most);
        if !(depth >= 0.0 && depth <= max) {
            return None;
        }
        let [pa, pb, pc] = *corners;
        let at = [0, 1, 2]
            .map(|k| barycentric[0] * pa[k] + barycentric[1] * pb[k] + barycentric[2] * pc[k]);
        Some(Hit {
            triangle: 0,
            point: point(at),
            distance: depth,
            barycentric,
        })
    }
}

/// Build the nodes over the triangles' boxes, reordering `order` into leaf
/// order: depth first, each left child right after its parent.
fn build_nodes(boxes: &[(P, P)], order: &mut [u32]) -> Vec<Node> {
    let mut nodes: Vec<Node> = Vec::new();
    if order.is_empty() {
        return nodes;
    }
    let centre = |t: u32| -> P {
        let (low, high) = boxes[t as usize];
        [0, 1, 2].map(|k| 0.5 * (low[k] + high[k]))
    };
    // Each pending range carries the parent whose right child it is.
    let mut pending: Vec<(usize, usize, Option<usize>)> = vec![(0, order.len(), None)];
    while let Some((start, end, parent)) = pending.pop() {
        let index = nodes.len();
        if let Some(parent) = parent {
            nodes[parent].first = u32::try_from(index).unwrap_or(u32::MAX);
        }
        let range = &mut order[start..end];
        let (mut low, mut high) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        let (mut clow, mut chigh) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for &t in range.iter() {
            let (l, h) = boxes[t as usize];
            let c = centre(t);
            for k in 0..3 {
                low[k] = low[k].min(l[k]);
                high[k] = high[k].max(h[k]);
                clow[k] = clow[k].min(c[k]);
                chigh[k] = chigh[k].max(c[k]);
            }
        }
        let count = end - start;
        let leaf = Node {
            low,
            high,
            first: u32::try_from(start).unwrap_or(u32::MAX),
            count: u32::try_from(count).unwrap_or(u32::MAX),
        };
        if count <= 2 {
            nodes.push(leaf);
            continue;
        }
        let mid = match sah_split(boxes, range, &clow, &chigh, area(&low, &high)) {
            Some(split) if split.pays || count > LEAF_MOST => {
                let bin_of = binner(&clow, &chigh, split.axis);
                start + partition(range, |t| bin_of(centre(t)[split.axis]) <= split.bin)
            }
            None if count > LEAF_MOST => start + count / 2,
            _ => {
                nodes.push(leaf);
                continue;
            }
        };
        nodes.push(Node {
            first: 0,
            count: 0,
            ..leaf
        });
        pending.push((mid, end, Some(index)));
        pending.push((start, mid, None));
    }
    nodes
}

/// The cheapest split of a range: after `bin` along `axis`.
struct Split {
    axis: usize,
    bin: usize,
    /// Whether it is cheaper than keeping the range as one leaf.
    pays: bool,
}

fn area(low: &P, high: &P) -> f64 {
    let d = [0, 1, 2].map(|k| (high[k] - low[k]).max(0.0));
    d[0] * d[1] + d[1] * d[2] + d[2] * d[0]
}

/// The bin a centroid coordinate falls in along `axis`.
fn binner(clow: &P, chigh: &P, axis: usize) -> impl Fn(f64) -> usize {
    let low = clow[axis];
    let scale = bin_count_f64() / (chigh[axis] - low);
    move |c: f64| {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let bin = ((c - low) * scale) as usize;
        bin.min(BINS - 1)
    }
}

#[allow(clippy::cast_precision_loss)]
fn bin_count_f64() -> f64 {
    BINS as f64
}

/// The cheapest split by the surface-area heuristic, binning the centroids
/// along each axis they spread over: a split costs one step down plus the
/// triangles of each side weighted by its share of the parent's area,
/// against a leaf's cost of every triangle. `None` where the centroids do
/// not spread along any axis.
fn sah_split(boxes: &[(P, P)], range: &[u32], clow: &P, chigh: &P, whole: f64) -> Option<Split> {
    let mut best: Option<(f64, usize, usize)> = None;
    for axis in 0..3 {
        if chigh[axis] <= clow[axis] {
            continue;
        }
        let bin_of = binner(clow, chigh, axis);
        let mut counts = [0_usize; BINS];
        let mut lows = [[f64::INFINITY; 3]; BINS];
        let mut highs = [[f64::NEG_INFINITY; 3]; BINS];
        for &t in range {
            let (l, h) = boxes[t as usize];
            let b = bin_of(0.5 * (l[axis] + h[axis]));
            counts[b] += 1;
            for k in 0..3 {
                lows[b][k] = lows[b][k].min(l[k]);
                highs[b][k] = highs[b][k].max(h[k]);
            }
        }
        // The right side's area and count for a split after each bin.
        let mut right = [(0.0, 0_usize); BINS];
        let (mut low, mut high, mut n) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3], 0);
        for b in (1..BINS).rev() {
            for k in 0..3 {
                low[k] = low[k].min(lows[b][k]);
                high[k] = high[k].max(highs[b][k]);
            }
            n += counts[b];
            right[b - 1] = (if n > 0 { area(&low, &high) } else { 0.0 }, n);
        }
        let (mut low, mut high, mut n) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3], 0);
        for b in 0..BINS - 1 {
            for k in 0..3 {
                low[k] = low[k].min(lows[b][k]);
                high[k] = high[k].max(highs[b][k]);
            }
            n += counts[b];
            let (ra, rn) = right[b];
            if n == 0 || rn == 0 {
                continue;
            }
            #[allow(clippy::cast_precision_loss)]
            let cost = area(&low, &high) * n as f64 + ra * rn as f64;
            if best.is_none_or(|(held, _, _)| cost < held) {
                best = Some((cost, axis, b));
            }
        }
    }
    let (cost, axis, bin) = best?;
    #[allow(clippy::cast_precision_loss)]
    let leaf = range.len() as f64;
    Some(Split {
        axis,
        bin,
        pays: whole <= 0.0 || 1.0 + cost / whole < leaf,
    })
}

/// Move the items `left` keeps to the front, keeping the order within each
/// side; how many it kept.
fn partition(range: &mut [u32], left: impl Fn(u32) -> bool) -> usize {
    let mut kept = 0;
    let mut rest = Vec::new();
    for i in 0..range.len() {
        let t = range[i];
        if left(t) {
            range[kept] = t;
            kept += 1;
        } else {
            rest.push(t);
        }
    }
    range[kept..].copy_from_slice(&rest);
    kept
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::cast_precision_loss)]
pub(crate) mod tests {
    use super::*;

    /// A small deterministic generator (xorshift), so every run draws the
    /// same meshes and queries.
    pub(crate) struct Seeded(pub(crate) u64);

    impl Seeded {
        pub(crate) fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }

        /// Uniform in `[lo, hi)`.
        pub(crate) fn range(&mut self, lo: f64, hi: f64) -> f64 {
            #[allow(clippy::cast_precision_loss)]
            let unit = (self.next() >> 11) as f64 / (1_u64 << 53) as f64;
            lo + (hi - lo) * unit
        }

        pub(crate) fn below(&mut self, n: usize) -> usize {
            usize::try_from(self.next() % n as u64).unwrap()
        }

        pub(crate) fn point(&mut self, lo: f64, hi: f64) -> Point {
            Point::new(self.range(lo, hi), self.range(lo, hi), self.range(lo, hi))
        }
    }

    fn corners_of(mesh: &Triangulation) -> Vec<[P; 3]> {
        mesh.triangles
            .iter()
            .map(|t| {
                t.map(|k| {
                    let p = mesh.positions[k as usize];
                    [p.x, p.y, p.z]
                })
            })
            .collect()
    }

    fn brute_nearest(corners: &[[P; 3]], p: Point, max: f64) -> Option<Hit> {
        let q = [p.x, p.y, p.z];
        let reach = max * max;
        let mut best: Option<(f64, usize)> = None;
        for (i, t) in corners.iter().enumerate() {
            let d2 = closest(q, t).2;
            if d2 <= reach && best.is_none_or(|(held, _)| d2 < held) {
                best = Some((d2, i));
            }
        }
        let (_, i) = best?;
        let (at, barycentric, d2) = closest(q, &corners[i]);
        Some(Hit {
            triangle: u32::try_from(i).unwrap(),
            point: point(at),
            distance: d2.sqrt(),
            barycentric,
        })
    }

    fn brute_ball(corners: &[[P; 3]], c: Point, r: f64) -> Vec<u32> {
        let q = [c.x, c.y, c.z];
        (0..corners.len())
            .filter(|&i| closest(q, &corners[i]).2 <= r * r)
            .map(|i| u32::try_from(i).unwrap())
            .collect()
    }

    fn brute_ray(corners: &[[P; 3]], o: Point, d: Vector, max: f64, both: bool) -> Option<Hit> {
        let shear = Shear::new(o, d)?;
        let mut best: Option<Hit> = None;
        for (i, t) in corners.iter().enumerate() {
            if let Some(hit) = shear.triangle(t, max, both)
                && best.is_none_or(|held| hit.distance < held.distance)
            {
                best = Some(Hit {
                    triangle: u32::try_from(i).unwrap(),
                    ..hit
                });
            }
        }
        best
    }

    /// A closed surface: a lumpy torus of `around` by `across` quads, two
    /// triangles each, wound outward.
    fn torus(around: u32, across: u32, rng: &mut Seeded) -> Triangulation {
        let mut mesh = Triangulation::new();
        for i in 0..around {
            for j in 0..across {
                let u = std::f64::consts::TAU * f64::from(i) / f64::from(around);
                let v = std::f64::consts::TAU * f64::from(j) / f64::from(across);
                let lump = 1.0 + 0.05 * rng.range(-1.0, 1.0);
                let r = (2.0 + 0.7 * v.cos()) * lump;
                mesh.positions
                    .push(Point::new(r * u.cos(), r * u.sin(), 0.7 * v.sin() * lump));
            }
        }
        let at = |i: u32, j: u32| (i % around) * across + (j % across);
        for i in 0..around {
            for j in 0..across {
                let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
                mesh.triangles.push([a, b, c]);
                mesh.triangles.push([a, c, d]);
            }
        }
        mesh
    }

    /// A soup of random triangles, large and small, with degenerate ones
    /// mixed in: a repeated corner, three collinear corners, one point.
    fn soup(count: usize, rng: &mut Seeded) -> Triangulation {
        let mut mesh = Triangulation::new();
        for i in 0..count {
            let a = rng.point(-1.0, 1.0);
            let size = if i % 3 == 0 { 0.5 } else { 0.05 };
            let b = a + (rng.point(-1.0, 1.0) - Point::ORIGIN) * size;
            let c = a + (rng.point(-1.0, 1.0) - Point::ORIGIN) * size;
            let corners = match i % 17 {
                3 => [a, b, b],
                5 => [a, b, a + (b - a) * 0.25],
                7 => [a, a, a],
                _ => [a, b, c],
            };
            let first = u32::try_from(mesh.positions.len()).unwrap();
            mesh.positions.extend(corners);
            mesh.triangles.push([first, first + 1, first + 2]);
        }
        mesh
    }

    fn direction(rng: &mut Seeded) -> Vector {
        loop {
            let v = rng.point(-1.0, 1.0) - Point::ORIGIN;
            if v.magnitude() > 0.1 {
                return v;
            }
        }
    }

    /// Every query kind on every mesh agrees with brute force, bit for bit:
    /// 1e5 queries of each kind over closed meshes and soups.
    #[test]
    fn every_query_equals_brute_force_on_random_meshes() {
        let mut rng = Seeded(0x05ee_d577);
        let mut meshes = Vec::new();
        for _ in 0..5 {
            meshes.push(torus(
                12 + u32::try_from(rng.below(12)).unwrap(),
                8,
                &mut rng,
            ));
            meshes.push(soup(150 + rng.below(250), &mut rng));
        }
        let per = 100_000 / meshes.len();
        let mut found = [0_usize; 4];
        for mesh in &meshes {
            let bvh = Bvh::build(mesh).unwrap();
            let corners = corners_of(mesh);
            // On a corner, on an edge, near the surface, or anywhere.
            let target = |rng: &mut Seeded, kind: usize| -> Point {
                let t = &corners[rng.below(corners.len())];
                let [a, b, c] = t.map(point);
                match kind {
                    0 => a,
                    1 => a + (b - a) * rng.range(0.0, 1.0),
                    2 => {
                        let (u, v) = (rng.range(0.0, 1.0), rng.range(0.0, 1.0));
                        let (u, v) = if u + v > 1.0 {
                            (1.0 - u, 1.0 - v)
                        } else {
                            (u, v)
                        };
                        a + (b - a) * u + (c - a) * v + (rng.point(-0.05, 0.05) - Point::ORIGIN)
                    }
                    _ => rng.point(-3.5, 3.5),
                }
            };
            for q in 0..per {
                let p = target(&mut rng, q % 4);
                let max = if q % 4 == 0 {
                    f64::INFINITY
                } else {
                    rng.range(0.0, 1.5)
                };
                let near = bvh.nearest(p, max);
                assert_eq!(near, brute_nearest(&corners, p, max), "{p:?} {max}");
                found[0] += usize::from(near.is_some());
                let r = rng.range(0.0, 0.6);
                let ball = bvh.within_ball(p, r);
                assert_eq!(ball, brute_ball(&corners, p, r), "{p:?} {r}");
                found[1] += usize::from(!ball.is_empty());
                // From anywhere, through a corner, an edge or a face.
                let o = rng.point(-3.5, 3.5);
                let aim = target(&mut rng, (q / 4) % 3);
                let d = if aim == o {
                    direction(&mut rng)
                } else {
                    aim - o
                };
                let max = if q % 3 == 0 {
                    f64::INFINITY
                } else {
                    rng.range(0.0, 6.0)
                };
                for both in [false, true] {
                    let hit = bvh.ray(o, d, max, both);
                    assert_eq!(
                        hit,
                        brute_ray(&corners, o, d, max, both),
                        "{o:?} {d:?} {max} {both}"
                    );
                    found[2 + usize::from(both)] += usize::from(hit.is_some());
                }
            }
        }
        // Each kind found something often enough to have been tested.
        assert!(found.iter().all(|&n| n > 25_000), "{found:?}");
    }

    /// The batch queries answer exactly as the single ones, in input order,
    /// at any thread count.
    #[test]
    fn batches_answer_as_single_queries_at_any_thread_count() {
        let mut rng = Seeded(77);
        let mesh = torus(30, 12, &mut rng);
        let bvh = Bvh::build(&mesh).unwrap();
        let points: Vec<Point> = (0..3000).map(|_| rng.point(-3.0, 3.0)).collect();
        let dirs: Vec<Vector> = (0..3000).map(|_| direction(&mut rng)).collect();
        let near: Vec<Option<Hit>> = points.iter().map(|&p| bvh.nearest(p, 1.0)).collect();
        let rays: Vec<Option<Hit>> = points
            .iter()
            .zip(&dirs)
            .map(|(&p, &d)| bvh.ray(p, d, f64::INFINITY, true))
            .collect();
        for threads in [1, 3, 8] {
            ogeom_core::parallel::set_threads(threads);
            assert_eq!(bvh.nearest_many(&points, 1.0), near);
            assert_eq!(bvh.rays(&points, &dirs, f64::INFINITY, true).unwrap(), rays);
        }
        ogeom_core::parallel::set_threads(0);
        assert!(bvh.rays(&points, &dirs[1..], 1.0, true).is_err());
    }

    /// A ray aimed exactly through a vertex or along an edge of a closed
    /// mesh, from inside it, always strikes the mesh: nothing slips between
    /// the triangles that meet there.
    #[test]
    fn rays_through_edges_and_vertices_do_not_slip_through() {
        let mut rng = Seeded(9);
        let mesh = torus(40, 16, &mut rng);
        let bvh = Bvh::build(&mesh).unwrap();
        let mut aimed = 0;
        for t in &mesh.triangles {
            let [a, b, _] = t.map(|k| mesh.positions[k as usize]);
            // From a point on the torus's core circle under the edge, so
            // the ray starts inside the solid.
            let mid = Point::new((a.x + b.x) / 2.0, (a.y + b.y) / 2.0, 0.0);
            let core = (mid - Point::ORIGIN) / (mid - Point::ORIGIN).magnitude() * 2.0;
            let inside = Point::ORIGIN + core;
            for target in [a, b, a + (b - a) * 0.5, a + (b - a) * 0.25] {
                // Outward-wound, so met from inside on its back.
                let hit = bvh.ray(inside, target - inside, f64::INFINITY, true);
                assert!(
                    hit.is_some(),
                    "a ray from {inside:?} through {target:?} slipped"
                );
                aimed += 1;
            }
        }
        assert_eq!(aimed, 4 * mesh.triangles.len());
        // Straight down onto a flat grid's vertices and edges, where
        // every corner is exactly representable and shared by six
        // triangles.
        let mut grid = Triangulation::new();
        for i in 0..5_u32 {
            for j in 0..5_u32 {
                grid.positions
                    .push(Point::new(f64::from(i) * 0.25, f64::from(j) * 0.25, 0.0));
            }
        }
        for i in 0..4 {
            for j in 0..4 {
                let (a, b, c, d) = (
                    i * 5 + j,
                    (i + 1) * 5 + j,
                    (i + 1) * 5 + j + 1,
                    i * 5 + j + 1,
                );
                grid.triangles.push([a, b, c]);
                grid.triangles.push([a, c, d]);
            }
        }
        let flat = Bvh::build(&grid).unwrap();
        for i in 0..=8 {
            for j in 0..=8 {
                let at = Point::new(f64::from(i) * 0.125, f64::from(j) * 0.125, 1.0);
                let hit = flat
                    .ray(at, Vector::new(0.0, 0.0, -1.0), 2.0, false)
                    .expect("a ray onto the grid hits it");
                assert_eq!(hit.distance, 1.0);
                // The lowest index among the triangles that meet there.
                let corners = corners_of(&grid);
                let least = brute_ray(&corners, at, Vector::new(0.0, 0.0, -1.0), 2.0, false);
                assert_eq!(Some(hit), least);
            }
        }
    }

    /// A triangle's front is the side its counter-clockwise winding points
    /// out of. A back-facing hit counts only when both sides are asked for.
    #[test]
    fn a_ray_meets_the_front_or_both_sides_at_the_closed_form_depth() {
        let mut mesh = Triangulation::new();
        mesh.positions = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(4.0, 0.0, 0.0),
            Point::new(0.0, 4.0, 0.0),
        ];
        mesh.triangles = vec![[0, 1, 2]];
        let bvh = Bvh::build(&mesh).unwrap();
        let down = bvh
            .ray(
                Point::new(1.0, 2.0, 3.0),
                Vector::new(0.0, 0.0, -2.0),
                10.0,
                false,
            )
            .unwrap();
        assert_eq!(down.distance, 3.0);
        assert_eq!(down.point, Point::new(1.0, 2.0, 0.0));
        assert_eq!(down.barycentric, [0.25, 0.25, 0.5]);
        let up = Vector::new(0.0, 0.0, 1.0);
        assert!(
            bvh.ray(Point::new(1.0, 2.0, -3.0), up, 10.0, false)
                .is_none()
        );
        assert_eq!(
            bvh.ray(Point::new(1.0, 2.0, -3.0), up, 10.0, true)
                .unwrap()
                .distance,
            3.0
        );
        // Every dominant axis and both senses of it agree on the front.
        for (o, d) in [
            (Point::new(1.0, 1.0, 1.0), Vector::new(0.1, -0.2, -1.0)),
            (Point::new(1.5, 1.0, 1.0), Vector::new(-1.0, 0.0, -0.9)),
            (Point::new(1.0, 2.0, 1.0), Vector::new(0.3, -1.0, -0.8)),
            (Point::new(3.0, 0.5, 1.0), Vector::new(-1.0, 0.1, -0.5)),
            (Point::new(0.5, 3.0, 1.0), Vector::new(0.0, -1.0, -0.5)),
            (Point::new(0.2, 1.0, 1.0), Vector::new(1.0, 0.2, -0.6)),
            (Point::new(1.0, 0.2, 1.0), Vector::new(0.2, 1.0, -0.6)),
        ] {
            assert!(bvh.ray(o, d, 10.0, false).is_some(), "{d:?}");
            let below = Point::new(o.x, o.y, -o.z);
            let rising = Vector::new(d.x, d.y, -d.z);
            assert!(bvh.ray(below, rising, 10.0, false).is_none(), "{rising:?}");
            assert!(bvh.ray(below, rising, 10.0, true).is_some(), "{rising:?}");
        }
        // Past `max`, or behind the origin, is no hit.
        assert!(
            bvh.ray(Point::new(1.0, 2.0, 3.0), -up, 2.9, false)
                .is_none()
        );
        assert!(bvh.ray(Point::new(1.0, 2.0, 3.0), up, 10.0, true).is_none());
        // Every hit along a line through the triangle, from either side.
        assert_eq!(
            bvh.hits_along(Point::new(1.0, 2.0, 3.0), -up, 10.0, true)
                .len(),
            1
        );
    }

    /// Queries exactly on a vertex or an edge answer distance zero there,
    /// and a vertex six triangles share goes to the lowest of them.
    #[test]
    fn queries_on_vertices_and_edges_land_on_them() {
        let mut rng = Seeded(31);
        let mesh = torus(16, 8, &mut rng);
        let bvh = Bvh::build(&mesh).unwrap();
        for (k, &p) in mesh.positions.iter().enumerate() {
            let hit = bvh.nearest(p, f64::INFINITY).unwrap();
            assert_eq!(hit.distance, 0.0);
            assert_eq!(hit.point, p);
            let k = u32::try_from(k).unwrap();
            let least = mesh.triangles.iter().position(|t| t.contains(&k)).unwrap();
            assert_eq!(hit.triangle as usize, least);
            assert!(bvh.within_ball(p, 0.0).contains(&hit.triangle));
        }
        for t in &mesh.triangles {
            let [a, b, _] = t.map(|k| mesh.positions[k as usize]);
            let on = a + (b - a) * 0.5;
            let hit = bvh.nearest(on, f64::INFINITY).unwrap();
            assert!(hit.distance <= 1e-14, "{}", hit.distance);
        }
    }

    /// A degenerate triangle is still a set of points a query can land on:
    /// a repeated corner is its segment, a collapsed one is its point.
    #[test]
    fn degenerate_triangles_answer_as_their_segments_and_points() {
        let mut mesh = Triangulation::new();
        mesh.positions = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
            Point::new(50.0, 50.0, 50.0),
            Point::new(1.0, 0.0, 0.0),
        ];
        // A segment twice over (a repeated corner, then collinear corners)
        // and a single point.
        mesh.triangles = vec![[0, 1, 2], [0, 4, 1], [3, 3, 3]];
        let bvh = Bvh::build(&mesh).unwrap();
        let hit = bvh
            .nearest(Point::new(1.5, 3.0, 4.0), f64::INFINITY)
            .unwrap();
        assert_eq!(hit.triangle, 0);
        assert_eq!(hit.distance, 5.0);
        assert_eq!(hit.point, Point::new(1.5, 0.0, 0.0));
        let hit = bvh
            .nearest(Point::new(50.0, 50.0, 52.0), f64::INFINITY)
            .unwrap();
        assert_eq!((hit.triangle, hit.distance), (2, 2.0));
        assert_eq!(bvh.within_ball(Point::new(1.0, 1.0, 0.0), 1.0), vec![0, 1]);
        // A ray in no triangle's plane meets none of them: they have no
        // area to cross.
        assert!(
            bvh.ray(
                Point::new(1.0, 0.0, 1.0),
                Vector::new(0.0, 0.0, -1.0),
                5.0,
                true
            )
            .is_none()
        );
    }

    #[test]
    fn an_empty_mesh_answers_nothing() {
        let bvh = Bvh::build(&Triangulation::new()).unwrap();
        assert!(bvh.is_empty());
        assert!(bvh.bounds().is_empty());
        let p = Point::new(0.0, 0.0, 0.0);
        assert!(bvh.nearest(p, f64::INFINITY).is_none());
        assert!(bvh.within_ball(p, 1e9).is_empty());
        assert!(
            bvh.ray(p, Vector::new(1.0, 0.0, 0.0), f64::INFINITY, true)
                .is_none()
        );
        assert!(
            bvh.hits_along(p, Vector::new(1.0, 0.0, 0.0), 1.0, true)
                .is_empty()
        );
        assert!(bvh.triangles_where(|_| true).is_empty());
        assert_eq!(bvh.nearest_many(&[p, p], 1.0), vec![None, None]);
        assert!(bvh.triangle(0).is_none());
    }

    /// A negative or undefined reach finds nothing, even where a triangle
    /// is as far as its magnitude.
    #[test]
    fn a_negative_or_undefined_reach_finds_nothing() {
        let mut mesh = Triangulation::new();
        mesh.positions = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        mesh.triangles = vec![[0, 1, 2]];
        let bvh = Bvh::build(&mesh).unwrap();
        let above = Point::new(0.2, 0.2, 1.0);
        assert!(bvh.nearest(above, 1.0).is_some());
        for bad in [-1.0, -2.0, f64::NAN] {
            assert!(bvh.nearest(above, bad).is_none());
            assert!(bvh.within_ball(above, bad).is_empty());
            assert!(
                bvh.ray(above, Vector::new(0.0, 0.0, -1.0), bad, true)
                    .is_none()
            );
        }
    }

    #[test]
    fn a_mesh_naming_a_missing_or_infinite_corner_is_refused() {
        let mut mesh = Triangulation::new();
        mesh.positions = vec![Point::new(0.0, 0.0, 0.0); 2];
        mesh.triangles = vec![[0, 1, 2]];
        assert!(Bvh::build(&mesh).is_err());
        mesh.positions.push(Point::new(f64::NAN, 0.0, 0.0));
        assert!(Bvh::build(&mesh).is_err());
    }

    /// Over a meshed box the nearest point to one outside a face is its foot
    /// on that face, and each triangle's corners come back by index.
    #[test]
    fn a_shape_s_tessellation_indexes_as_its_triangulation() {
        use ogeom_core::Tolerances;
        let tol = Tolerances::millimetres();
        let mut model = Model::new();
        let solid =
            ogeom_algo::make_box(&mut model, ogeom_math::Frame::WORLD, (2.0, 3.0, 4.0), tol)
                .unwrap()
                .shape;
        let deflection = Deflection::default();
        let bvh = Bvh::of_shape(&model, &solid, deflection, tol).unwrap();
        let mesh = crate::triangulate::triangulate(&model, &solid, deflection, tol).unwrap();
        assert_eq!(bvh.triangle_count(), mesh.triangles.len());
        for (i, t) in mesh.triangles.iter().enumerate() {
            let corners = t.map(|k| mesh.positions[k as usize]);
            assert_eq!(bvh.triangle(u32::try_from(i).unwrap()), Some(corners));
        }
        let hit = bvh
            .nearest(Point::new(1.0, 1.5, 9.0), f64::INFINITY)
            .unwrap();
        assert_eq!(hit.distance, 5.0);
        assert_eq!(hit.point, Point::new(1.0, 1.5, 4.0));
        let found = bvh.triangles_where(|b| b.contains(Point::new(1.0, 1.5, 4.0)));
        assert!(found.contains(&hit.triangle));
    }
}
