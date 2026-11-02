//! Regions of the plane bounded by weighted chains of segments, found
//! robustly: the engine under [`shadow`](crate::shadow()) and the outline
//! fit's offsets.
//!
//! A set of directed segments with integer weights, closed as a chain (at
//! every point the weights in equal the weights out), gives every point off
//! the segments a winding count: the weights of the segments a ray from it
//! crosses, from right to left minus from left to right. [`region`] returns
//! the boundary of the points whose counts a rule accepts. Two independent
//! counts are carried, so a rule can take a difference as well as a union.
//!
//! Every decision is exact. The endpoints are first snapped to a grid of
//! `2^-38` of the chain's extent, a power of two, and from then on every
//! coordinate is an integer and every test is integer arithmetic: the
//! orientation of three points, where two segments cross, whether a segment
//! meets a grid cell. The crossings are resolved by snap rounding: every
//! cell holding an endpoint or a crossing is *hot*, and every segment is
//! replaced by the path through the centres of the hot cells it meets, in
//! the order it meets them. The paths cross nowhere but at shared centres,
//! each moved by at most `sqrt(2)/2` of a cell from its segment, and a
//! segment's path starts and ends where the segment does, so the chain
//! stays closed and every count stays defined.
//!
//! The paths are then the constraints of a triangulation. A count changes
//! across a constraint by its weight and across nothing else, so one walk
//! from the unbounded face, whose counts are zero, gives every triangle's
//! counts. The boundary is every triangle edge with an accepted triangle on
//! its left and a rejected one on its right, linked into loops by turning
//! through accepted triangles only, so two loops touching at a corner come
//! out as two loops. Points where a loop runs straight on are dropped.

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::Point2;
use spade::handles::{FixedDirectedEdgeHandle, FixedVertexHandle};
use spade::{ConstrainedDelaunayTriangulation, Triangulation as _};

/// One closed loop of a region's boundary.
///
/// The region lies on the left: an outer boundary runs counter-clockwise
/// and a hole's clockwise.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Loop2 {
    /// The corners in order. The first is not repeated at the end.
    pub points: Vec<Point2>,
}

impl Loop2 {
    /// The signed area the loop encloses: positive for an outer boundary,
    /// negative for a hole.
    #[must_use]
    pub fn area(&self) -> f64 {
        let n = self.points.len();
        let mut twice = 0.0;
        for i in 0..n {
            let (a, b) = (self.points[i], self.points[(i + 1) % n]);
            twice += a.x.mul_add(b.y, -(a.y * b.x));
        }
        twice / 2.0
    }

    /// Whether the loop bounds a hole: whether it runs clockwise.
    #[must_use]
    pub fn is_hole(&self) -> bool {
        self.area() < 0.0
    }
}

/// A directed segment of a chain, with its weight in each of the two
/// counts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Edge2 {
    /// Where the segment starts.
    pub from: Point2,
    /// Where it ends.
    pub to: Point2,
    /// How much it adds to each count of the points on its left, against
    /// the points on its right.
    pub weight: [i32; 2],
}

/// The boundary of the points whose winding counts `inside` accepts.
///
/// The loops are within half a cell's diagonal of the exact boundary, a
/// cell being the power of two at or above `2^-38` of the widest extent of
/// the endpoints: at most `2^-37` of it.
///
/// The edges must form a closed chain in each count: at every endpoint the
/// weights of the edges arriving sum to those of the edges leaving. A
/// chain that is not closed has no counts, and its region is whatever the
/// walk from the unbounded face meets first. `inside` must reject `[0, 0]`,
/// the counts of every point far away.
///
/// # Errors
///
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if an endpoint
/// is not finite, and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `inside` accepts `[0, 0]`.
pub fn region(edges: &[Edge2], inside: impl Fn([i64; 2]) -> bool) -> OgeomResult<Vec<Loop2>> {
    if inside([0, 0]) {
        ogeom_bail!(
            Construction,
            "a region accepting the unbounded face has no boundary"
        );
    }
    if let Some(e) = edges
        .iter()
        .find(|e| !e.from.is_finite() || !e.to.is_finite())
    {
        ogeom_bail!(Numeric, "a chain edge from {:?} is not finite", e.from);
    }
    let Some(grid) = Grid::around(edges.iter().flat_map(|e| [e.from, e.to])) else {
        return Ok(Vec::new());
    };
    let segments = edges
        .iter()
        .map(|e| Seg {
            a: grid.snap(e.from),
            b: grid.snap(e.to),
            w: e.weight,
        })
        .collect();
    Ok(trace(segments, &inside)?
        .into_iter()
        .map(|ring| Loop2 {
            points: ring.into_iter().map(|q| grid.point(q)).collect(),
        })
        .collect())
}

/// A grid coordinate as a double, exactly: grid coordinates are below
/// `2^41`.
#[allow(
    clippy::cast_precision_loss,
    reason = "grid coordinates are below 2^41, exact in a double"
)]
pub(crate) const fn real(v: i64) -> f64 {
    v as f64
}

/// A count as a double, to size buckets by.
#[allow(clippy::cast_precision_loss, reason = "an estimate, rounding harmless")]
pub(crate) const fn estimate(n: usize) -> f64 {
    n as f64
}

/// A double holding an integer of grid size, as one.
#[allow(
    clippy::cast_possible_truncation,
    reason = "callers pass values rounded to integers below 2^42"
)]
pub(crate) const fn whole(v: f64) -> i64 {
    v as i64
}

/// A bucket or point index as the `u32` the lists keep. Bucket counts are
/// capped at `4096^2` and points at the `2^32` the triangulation allows.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "indices are below 2^32"
)]
pub(crate) const fn slot(v: i64) -> u32 {
    v as u32
}

/// A point on the grid, in units of half a cell: cell centres are even.
pub(crate) type Q = [i64; 2];

/// The grid points snap to: square cells, a power of two wide, from an
/// origin on a cell corner's multiple.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Grid {
    origin: [f64; 2],
    cell: f64,
}

/// How many cells the widest extent spans at most.
const SPAN_BITS: i32 = 38;

impl Grid {
    /// The grid for points spread over this box, or `None` for none.
    pub(crate) fn around(points: impl Iterator<Item = Point2>) -> Option<Self> {
        let mut lo = [f64::INFINITY; 2];
        let mut hi = [f64::NEG_INFINITY; 2];
        for p in points {
            lo = [lo[0].min(p.x), lo[1].min(p.y)];
            hi = [hi[0].max(p.x), hi[1].max(p.y)];
        }
        if lo[0] > hi[0] {
            return None;
        }
        let extent = (hi[0] - lo[0]).max(hi[1] - lo[1]);
        // A single point spans nothing; any cell size then serves.
        let cell = if extent > 0.0 {
            2f64.powi(i32::try_from(whole(extent.log2().ceil())).unwrap_or(0) - SPAN_BITS)
        } else {
            1.0
        };
        Some(Self {
            origin: [(lo[0] / cell).floor() * cell, (lo[1] / cell).floor() * cell],
            cell,
        })
    }

    /// The centre of the cell a point falls in.
    pub(crate) fn snap(&self, p: Point2) -> Q {
        let k = |v: f64, o: f64| 2 * whole(((v - o) / self.cell).round());
        [k(p.x, self.origin[0]), k(p.y, self.origin[1])]
    }

    /// A grid point back in the plane.
    pub(crate) fn point(&self, q: Q) -> Point2 {
        let c = self.cell / 2.0;
        Point2::new(
            real(q[0]).mul_add(c, self.origin[0]),
            real(q[1]).mul_add(c, self.origin[1]),
        )
    }
}

/// A directed segment between grid points, and its weights.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Seg {
    pub(crate) a: Q,
    pub(crate) b: Q,
    pub(crate) w: [i32; 2],
}

/// Twice the signed area of `(a, b, c)`, exactly. Grid coordinates are
/// below `2^40`, so the products are below `2^80`.
pub(crate) fn orient(a: Q, b: Q, c: Q) -> i128 {
    let d = |p: Q, q: Q| [i128::from(q[0] - p[0]), i128::from(q[1] - p[1])];
    let (u, v) = (d(a, b), d(a, c));
    u[0] * v[1] - u[1] * v[0]
}

/// Sum the weights of segments running between the same two points either
/// way, and drop those that cancel or have no length.
fn cancel(mut segments: Vec<Seg>) -> Vec<Seg> {
    for s in &mut segments {
        if s.b < s.a {
            *s = Seg {
                a: s.b,
                b: s.a,
                w: [-s.w[0], -s.w[1]],
            };
        }
    }
    segments.sort_unstable_by_key(|s| (s.a, s.b));
    let mut out: Vec<Seg> = Vec::with_capacity(segments.len());
    for s in segments {
        match out.last_mut() {
            Some(t) if t.a == s.a && t.b == s.b => {
                t.w = [t.w[0] + s.w[0], t.w[1] + s.w[1]];
            }
            _ => out.push(s),
        }
    }
    out.retain(|s| s.a != s.b && s.w != [0, 0]);
    out
}

/// Uniform buckets over the grid, each listing what may lie in it.
struct Buckets {
    lo: Q,
    size: i64,
    columns: i64,
    rows: i64,
}

impl Buckets {
    fn over(segments: &[Seg]) -> Self {
        let mut lo = [i64::MAX; 2];
        let mut hi = [i64::MIN; 2];
        for s in segments {
            for p in [s.a, s.b] {
                lo = [lo[0].min(p[0]), lo[1].min(p[1])];
                hi = [hi[0].max(p[0]), hi[1].max(p[1])];
            }
        }
        // A margin of the two units the segment lists are widened by.
        lo = [lo[0] - 2, lo[1] - 2];
        hi = [hi[0] + 2, hi[1] + 2];
        let span = (hi[0] - lo[0]).max(hi[1] - lo[1]) + 1;
        let per_side = whole(estimate(segments.len()).sqrt().ceil()).clamp(1, 4096);
        let size = (span + per_side - 1) / per_side;
        Self {
            lo,
            size,
            columns: (hi[0] - lo[0]) / size + 1,
            rows: (hi[1] - lo[1]) / size + 1,
        }
    }

    fn column(&self, x: i64) -> i64 {
        ((x - self.lo[0]) / self.size).clamp(0, self.columns - 1)
    }

    fn row(&self, y: i64) -> i64 {
        ((y - self.lo[1]) / self.size).clamp(0, self.rows - 1)
    }

    fn of(&self, p: Q) -> u32 {
        slot(self.row(p[1]) * self.columns + self.column(p[0]))
    }

    /// Every bucket holding a point within two units of the segment, and
    /// perhaps a few more: within each row of buckets the segment's span
    /// is found in floating point, then widened by two units, far more than
    /// its rounding.
    fn along(&self, s: &Seg, mut visit: impl FnMut(u32)) {
        let (y_lo, y_hi) = (s.a[1].min(s.b[1]) - 2, s.a[1].max(s.b[1]) + 2);
        let (x_lo, x_hi) = (s.a[0].min(s.b[0]) - 2, s.a[0].max(s.b[0]) + 2);
        let dx = real(s.b[0] - s.a[0]);
        let dy = real(s.b[1] - s.a[1]);
        let x_at = |y: i64| (real(y - s.a[1]) * dx / dy) + real(s.a[0]);
        for r in self.row(y_lo)..=self.row(y_hi) {
            let (from, to) = if dy == 0.0 {
                (x_lo, x_hi)
            } else {
                let top = self.lo[1] + r * self.size;
                let y0 = top.max(y_lo);
                let y1 = (top + self.size - 1).min(y_hi);
                let (xa, xb) = (x_at(y0), x_at(y1));
                (
                    whole((xa.min(xb) - 2.0).floor()).max(x_lo),
                    whole((xa.max(xb) + 2.0).ceil()).min(x_hi),
                )
            };
            for c in self.column(from)..=self.column(to) {
                visit(slot(r * self.columns + c));
            }
        }
    }

    /// Compressed lists of the items each bucket holds.
    fn lists(&self, pairs: &mut [(u32, u32)]) -> (Vec<usize>, Vec<u32>) {
        pairs.sort_unstable();
        let count = slot(self.rows * self.columns) as usize;
        let mut start = vec![0usize; count + 1];
        for &(b, _) in pairs.iter() {
            start[b as usize + 1] += 1;
        }
        for k in 0..count {
            start[k + 1] += start[k];
        }
        (start, pairs.iter().map(|&(_, i)| i).collect())
    }
}

/// The cell of the point where two segments cross at a point inside both,
/// or `None` if they do not. Segments that touch or overlap meet at an
/// endpoint, whose cell is hot already.
fn crossing(s: &Seg, t: &Seg) -> Option<Q> {
    let s1 = orient(s.a, s.b, t.a).signum();
    let s2 = orient(s.a, s.b, t.b).signum();
    if s1 * s2 >= 0 {
        return None;
    }
    let along = orient(t.a, t.b, s.a);
    let other = orient(t.a, t.b, s.b);
    if along.signum() * other.signum() >= 0 {
        return None;
    }
    // The crossing is s.a + (s.b - s.a) * along / (along - other). Every
    // term is below 2^123; the cell's centre is the even integer nearest.
    let (mut num, mut den) = (along, along - other);
    if den < 0 {
        (num, den) = (-num, -den);
    }
    let cell = |a: i64, b: i64| {
        let x = i128::from(a) * den + i128::from(b - a) * num;
        i64::try_from(2 * (x + den).div_euclid(2 * den)).unwrap_or(i64::MAX)
    };
    Some([cell(s.a[0], s.b[0]), cell(s.a[1], s.b[1])])
}

/// A bound on a segment's parameter, `num / den` with `den > 0`, and
/// whether the bound itself is excluded.
#[derive(Debug, Clone, Copy)]
struct Bound {
    num: i128,
    den: i128,
    open: bool,
}

impl Bound {
    /// Compare the values, then a closed bound before an open one.
    fn cmp(self, other: Self) -> core::cmp::Ordering {
        self.value_cmp(other).then(self.open.cmp(&other.open))
    }

    /// Compare the values alone. Numerators are below `2^41` and
    /// denominators below `2^40`.
    fn value_cmp(self, other: Self) -> core::cmp::Ordering {
        (self.num * other.den).cmp(&(other.num * self.den))
    }
}

/// Where along `s` it first meets the half-open cell centred on `c`,
/// `[c - 1, c + 1)` on each axis, or `None` if it does not.
fn entry(s: &Seg, c: Q) -> Option<Bound> {
    let mut lo = Bound {
        num: 0,
        den: 1,
        open: false,
    };
    let mut hi = Bound {
        num: 1,
        den: 1,
        open: false,
    };
    for (k, &centre) in c.iter().enumerate() {
        let (p, d) = (i128::from(s.a[k]), i128::from(s.b[k] - s.a[k]));
        let (left, right) = (i128::from(centre - 1), i128::from(centre + 1));
        if d == 0 {
            if p < left || p >= right {
                return None;
            }
            continue;
        }
        // left <= p + t d < right.
        let (from, to) = if d > 0 {
            (
                Bound {
                    num: left - p,
                    den: d,
                    open: false,
                },
                Bound {
                    num: right - p,
                    den: d,
                    open: true,
                },
            )
        } else {
            (
                Bound {
                    num: p - right,
                    den: -d,
                    open: true,
                },
                Bound {
                    num: p - left,
                    den: -d,
                    open: false,
                },
            )
        };
        // The tighter of each pair: among equal values, the open one.
        if from.cmp(lo).is_gt() {
            lo = from;
        }
        let upper = to.value_cmp(hi);
        if upper.is_lt() || (upper.is_eq() && to.open) {
            hi = to;
        }
    }
    let ordered = lo.value_cmp(hi);
    let meets = ordered.is_lt() || (ordered.is_eq() && !lo.open && !hi.open);
    meets.then_some(lo)
}

/// Snap-round the segments and return the boundary of the accepted region
/// as loops of grid points.
pub(crate) fn trace(
    segments: Vec<Seg>,
    inside: &dyn Fn([i64; 2]) -> bool,
) -> OgeomResult<Vec<Vec<Q>>> {
    let segments = cancel(segments);
    if segments.is_empty() {
        return Ok(Vec::new());
    }
    let buckets = Buckets::over(&segments);
    let mut pairs = Vec::with_capacity(segments.len() * 2);
    for (i, s) in segments.iter().enumerate() {
        buckets.along(s, |b| pairs.push((b, slot(i as i64))));
    }
    let (start, held) = buckets.lists(&mut pairs);

    // Hot cells: every endpoint's, and every crossing's, each crossing
    // found once, in the bucket its cell's centre falls in.
    let mut hot: Vec<Q> = segments.iter().flat_map(|s| [s.a, s.b]).collect();
    for b in 0..start.len() - 1 {
        let here = &held[start[b]..start[b + 1]];
        for (k, &i) in here.iter().enumerate() {
            let s = &segments[i as usize];
            for &j in &here[k + 1..] {
                let t = &segments[j as usize];
                if let Some(c) = crossing(s, t)
                    && buckets.of(c) as usize == b
                {
                    hot.push(c);
                }
            }
        }
    }
    hot.sort_unstable();
    hot.dedup();
    let mut placed: Vec<(u32, u32)> = hot
        .iter()
        .enumerate()
        .map(|(k, &c)| (buckets.of(c), slot(k as i64)))
        .collect();
    let (hot_start, hot_held) = buckets.lists(&mut placed);

    // Each segment as the path through the hot cells it meets.
    let mut rounded = Vec::with_capacity(segments.len());
    let mut met: Vec<(Bound, Q)> = Vec::new();
    for s in &segments {
        met.clear();
        buckets.along(s, |b| {
            let b = b as usize;
            for &k in &hot_held[hot_start[b]..hot_start[b + 1]] {
                let c = hot[k as usize];
                if let Some(at) = entry(s, c) {
                    met.push((at, c));
                }
            }
        });
        met.sort_unstable_by(|x, y| x.0.cmp(y.0));
        for pair in met.windows(2) {
            rounded.push(Seg {
                a: pair[0].1,
                b: pair[1].1,
                w: s.w,
            });
        }
    }
    regions_of(&cancel(rounded), inside)
}

/// Triangulate non-crossing segments, count each triangle, and link the
/// edges between accepted and rejected triangles into loops.
fn regions_of(segments: &[Seg], inside: &dyn Fn([i64; 2]) -> bool) -> OgeomResult<Vec<Vec<Q>>> {
    if segments.is_empty() {
        return Ok(Vec::new());
    }
    let mut points: Vec<Q> = segments.iter().flat_map(|s| [s.a, s.b]).collect();
    points.sort_unstable();
    points.dedup();
    let index = |q: Q| points.binary_search(&q).unwrap_or(0);
    let vertices: Vec<spade::Point2<f64>> = points
        .iter()
        .map(|q| spade::Point2::new(real(q[0]), real(q[1])))
        .collect();
    let Ok(mut cdt) =
        ConstrainedDelaunayTriangulation::<spade::Point2<f64>>::bulk_load_cdt(vertices, Vec::new())
    else {
        ogeom_bail!(Numeric, "the region's points could not be triangulated");
    };
    let mut constraints: Vec<(FixedDirectedEdgeHandle, [i32; 2])> = Vec::new();
    for s in segments {
        let (from, to) = (index(s.a), index(s.b));
        let made = cdt.try_add_constraint(
            FixedVertexHandle::from_index(from),
            FixedVertexHandle::from_index(to),
        );
        if made.is_empty() {
            ogeom_bail!(
                Construction,
                "a snap-rounded segment crosses another; the region is not traced"
            );
        }
        let along = [s.b[0] - s.a[0], s.b[1] - s.a[1]];
        for h in made {
            let e = cdt.directed_edge(h);
            let [p, q] = e.positions();
            let forward = (q.x - p.x).mul_add(real(along[0]), (q.y - p.y) * real(along[1])) > 0.0;
            constraints.push((if forward { h } else { h.rev() }, s.w));
        }
    }
    // Weights in each direction of every edge, by directed edge index.
    let mut weight = vec![[0i64; 2]; cdt.num_directed_edges()];
    for (h, w) in constraints {
        for k in 0..2 {
            weight[h.index()][k] += i64::from(w[k]);
            weight[h.rev().index()][k] -= i64::from(w[k]);
        }
    }

    // Counts by face, walked from the unbounded face. The face on an
    // edge's left counts its weight more than the face on its right.
    let mut count: Vec<Option<[i64; 2]>> = vec![None; cdt.num_all_faces()];
    let outer = cdt.outer_face().fix().index();
    count[outer] = Some([0, 0]);
    let mut sides = vec![[None; 3]; cdt.num_all_faces()];
    for face in cdt.inner_faces() {
        sides[face.fix().index()] = face.adjacent_edges().map(|e| Some(e.fix()));
    }
    let mut stack = vec![outer];
    let mut around = Vec::with_capacity(3);
    while let Some(f) = stack.pop() {
        let Some(here) = count[f] else { continue };
        around.clear();
        if f == outer {
            for e in cdt.convex_hull() {
                around.push(if e.face().is_outer() {
                    e.fix()
                } else {
                    e.rev().fix()
                });
            }
        } else {
            around.extend(sides[f].iter().flatten());
        }
        for &e in &around {
            let w = weight[e.index()];
            let g = cdt.directed_edge(e).rev().face().fix().index();
            if count[g].is_none() {
                count[g] = Some([here[0] - w[0], here[1] - w[1]]);
                stack.push(g);
            }
        }
    }
    let accepted: Vec<bool> = count
        .iter()
        .enumerate()
        .map(|(f, c)| f != outer && c.is_some_and(inside))
        .collect();

    let boundary = |e: FixedDirectedEdgeHandle| {
        let h = cdt.directed_edge(e);
        accepted[h.face().fix().index()] && !accepted[h.rev().face().fix().index()]
    };
    let mut used = vec![false; cdt.num_directed_edges()];
    let mut loops = Vec::new();
    for e in cdt.directed_edges() {
        let first = e.fix();
        if used[first.index()] || !boundary(first) {
            continue;
        }
        let mut ring = Vec::new();
        let mut at = first;
        loop {
            used[at.index()] = true;
            let h = cdt.directed_edge(at);
            ring.push(points[h.from().fix().index()]);
            // Turn about the end through accepted faces to the next edge
            // with a rejected face on its right.
            let mut next = h.next();
            let mut turns = 0usize;
            while !boundary(next.fix()) {
                next = next.rev().next();
                turns += 1;
                if turns > cdt.num_directed_edges() {
                    ogeom_bail!(Construction, "a region boundary does not link up");
                }
            }
            at = next.fix();
            if at == first {
                break;
            }
            if used[at.index()] {
                ogeom_bail!(Construction, "a region boundary links into itself");
            }
        }
        loops.push(straightened(ring));
    }
    Ok(loops)
}

/// The loop without the points where it runs straight on.
fn straightened(mut ring: Vec<Q>) -> Vec<Q> {
    loop {
        let n = ring.len();
        if n < 3 {
            return ring;
        }
        let keep: Vec<bool> = (0..n)
            .map(|i| orient(ring[(i + n - 1) % n], ring[i], ring[(i + 1) % n]) != 0)
            .collect();
        if keep.iter().all(|&k| k) {
            return ring;
        }
        let mut k = 0;
        ring.retain(|_| {
            k += 1;
            keep[k - 1]
        });
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;

    fn square(x: f64, y: f64, s: f64) -> Vec<Edge2> {
        let p = [
            Point2::new(x, y),
            Point2::new(x + s, y),
            Point2::new(x + s, y + s),
            Point2::new(x, y + s),
        ];
        (0..4)
            .map(|i| Edge2 {
                from: p[i],
                to: p[(i + 1) % 4],
                weight: [1, 0],
            })
            .collect()
    }

    /// Two overlapping unit squares unite into one outline of eight
    /// corners, two of them where the sides cross, of area `2 - 1/4`, and
    /// meet in a square of area `1/4`.
    #[test]
    fn two_overlapping_squares_unite() {
        let mut edges = square(0.0, 0.0, 1.0);
        edges.extend(square(0.5, 0.5, 1.0));
        let loops = region(&edges, |w| w[0] > 0).unwrap();
        assert_eq!(loops.len(), 1);
        assert_eq!(loops[0].points.len(), 8);
        assert!((loops[0].area() - 1.75).abs() < 1e-9);
        let both = region(&edges, |w| w[0] > 1).unwrap();
        assert_eq!(both.len(), 1);
        assert!((both[0].area() - 0.25).abs() < 1e-9);
    }

    /// A square with a smaller one taken away in the second count keeps
    /// the smaller as a hole, clockwise.
    #[test]
    fn a_difference_keeps_its_hole() {
        let mut edges = square(0.0, 0.0, 3.0);
        edges.extend(square(1.0, 1.0, 1.0).into_iter().map(|e| Edge2 {
            weight: [0, 1],
            ..e
        }));
        let loops = region(&edges, |w| w[0] > 0 && w[1] == 0).unwrap();
        assert_eq!(loops.len(), 2);
        let mut areas: Vec<f64> = loops.iter().map(Loop2::area).collect();
        areas.sort_by(f64::total_cmp);
        assert!((areas[0] + 1.0).abs() < 1e-9 && (areas[1] - 9.0).abs() < 1e-9);
    }

    /// Squares touching at one corner stay two loops.
    #[test]
    fn squares_touching_at_a_corner_stay_apart() {
        let mut edges = square(0.0, 0.0, 1.0);
        edges.extend(square(1.0, 1.0, 1.0));
        let loops = region(&edges, |w| w[0] > 0).unwrap();
        assert_eq!(loops.len(), 2);
        assert!(loops.iter().all(|l| l.points.len() == 4));
    }

    /// Many thin slivers crossing near one point, the case where rounded
    /// crossings go wrong: the union of a fan of long thin triangles about
    /// a common centre is traced without error, and its area is the sum
    /// of the parts less what they share, which here is under the sum.
    #[test]
    fn a_fan_of_slivers_crossing_near_a_point() {
        let mut edges = Vec::new();
        for k in 0..64 {
            let a = f64::from(k) * 0.049_087_385;
            let (c, s) = (a.cos(), a.sin());
            let p = [
                Point2::new(-c + 1e-9 * s, -s - 1e-9 * c),
                Point2::new(c + 1e-9 * s, s - 1e-9 * c),
                Point2::new(c - 1e-3 * s, s + 1e-3 * c),
            ];
            for i in 0..3 {
                edges.push(Edge2 {
                    from: p[i],
                    to: p[(i + 1) % 3],
                    weight: [1, 0],
                });
            }
        }
        let loops = region(&edges, |w| w[0] > 0).unwrap();
        let area: f64 = loops.iter().map(Loop2::area).sum();
        assert!(area > 0.0 && area < 64.0 * 1e-3);
    }
}
