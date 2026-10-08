//! Surface feet on a B-spline patch by branch and bound over its knot spans.
//!
//! Each knot span of a patch lies inside the box of the `(p+1) x (q+1)`
//! control points it reads, since its basis functions are non-negative and
//! sum to one there, and a rational patch with positive weights keeps the
//! same hull; a block of spans lies inside the box of the control points
//! the block reads. A search pops the block nearest the target, quarters
//! it until it is one span, opens that span into its Bezier patch, and
//! halves patches whose bound still beats the nearest point found so far.
//! A patch's corners are points on the surface, so they give that nearest
//! point as the search goes, and the first patch opened polishes its
//! nearest corner by Newton into a foot. Once no open box can come nearer
//! than the nearest point less [`SETTLE`] of the confusion, its parameters
//! seed the final Newton
//! polish (Ma and Hewitt, CAGD 20, 2003; Selimovic, CAGD 23, 2006).

use super::{SurfaceProjection, refine_foot};
use crate::{BSplineSurface, SurfaceGeometry};
use core::cmp::Ordering;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_math::{Point, Weighted};
use smallvec::SmallVec;
use std::collections::BinaryHeap;

/// A homogeneous control point: the point times its weight, then the weight.
type Homogeneous = [f64; 4];

/// A Bezier patch's control net, `p + 1` rows of `q + 1` points.
type Net = SmallVec<[Homogeneous; 16]>;

/// One row or column of a net.
type Line = SmallVec<[Homogeneous; 8]>;

/// The fraction of the confusion the search settles the least distance
/// to, and within which two distances tie. Finer than the confusion: a
/// point on the patch projects onto itself, not merely within the
/// confusion of it, where Newton would stop.
const SETTLE: f64 = 0.01;

/// The search starts from at most this many blocks of spans per direction.
const BLOCKS: usize = 16;

/// At most this many boxes and patches opened per projection; the nearest
/// corner found by then seeds the polish. A target equidistant from a
/// whole region (the centre of a spherical cap) keeps every patch of it
/// open down to the confusion, and this is what ends that search.
const OPENED: usize = 16_384;

/// An axis-aligned box.
#[derive(Debug, Clone, Copy)]
struct Bounds {
    lo: [f64; 3],
    hi: [f64; 3],
}

impl Bounds {
    const EMPTY: Self = Self {
        lo: [f64::INFINITY; 3],
        hi: [f64::NEG_INFINITY; 3],
    };

    /// Grow to hold `p`, whose coordinates are finite.
    fn add(&mut self, p: [f64; 3]) {
        for ((lo, hi), x) in self.lo.iter_mut().zip(&mut self.hi).zip(p) {
            if x < *lo {
                *lo = x;
            }
            if x > *hi {
                *hi = x;
            }
        }
    }

    /// The distance from `t` to the box, zero inside it.
    fn gap(&self, t: [f64; 3]) -> f64 {
        let mut sum = 0.0;
        for ((lo, hi), x) in self.lo.iter().zip(&self.hi).zip(t) {
            let out = (lo - x).max(x - hi).max(0.0);
            sum += out * out;
        }
        sum.sqrt()
    }
}

/// A rectangle of spans: half-open index ranges along `u` and `v`.
type Block = ((u32, u32), (u32, u32));

/// A B-spline patch's spans and the boxes of its starting blocks: what a
/// search needs besides the patch itself, built once per patch.
#[derive(Debug, Clone)]
pub(super) struct SpanBoxes {
    /// Per direction, the index of the knot each span starts at.
    starts: (Vec<u32>, Vec<u32>),
    degree: (usize, usize),
    /// The reciprocal of the weight every control point shares, where
    /// they share one.
    uniform: Option<f64>,
    /// The starting blocks and their boxes.
    blocks: Vec<(Block, Bounds)>,
}

/// The patch a search reads: knots, degrees and net.
#[derive(Clone, Copy)]
struct Patch<'a> {
    knots: (&'a [f64], &'a [f64]),
    degree: (usize, usize),
    net: &'a [Weighted<Point>],
    v_count: usize,
    /// The reciprocal of the weight every control point shares, where
    /// they share one: then every point of every Bezier net carries that
    /// same weight, and a multiplication takes it out.
    uniform: Option<f64>,
}

impl Patch<'_> {
    fn of(spline: &BSplineSurface, uniform: Option<f64>) -> Patch<'_> {
        Patch {
            knots: (spline.u_knots().knots(), spline.v_knots().knots()),
            degree: (spline.u_knots().degree(), spline.v_knots().degree()),
            net: spline.grid().points(),
            v_count: spline.grid().v_count(),
            uniform,
        }
    }

    fn homogeneous(&self, i: usize, j: usize) -> Homogeneous {
        let w = self.net[i * self.v_count + j];
        [w.scaled.x, w.scaled.y, w.scaled.z, w.weight]
    }

    /// A homogeneous point of this patch or of a Bezier net cut from it,
    /// back in space.
    fn cartesian(&self, h: Homogeneous) -> [f64; 3] {
        match self.uniform {
            Some(k) => [h[0] * k, h[1] * k, h[2] * k],
            None => [h[0] / h[3], h[1] / h[3], h[2] / h[3]],
        }
    }
}

/// The spans of one direction: the knot index each starts at.
fn span_starts(knots: &ogeom_math::KnotVector) -> Vec<u32> {
    let all = knots.knots();
    #[allow(clippy::cast_possible_truncation)]
    (knots.degree()..knots.control_point_count())
        .filter(|&s| all[s] < all[s + 1])
        .map(|s| s as u32)
        .collect()
}

/// `count` indices cut into at most [`BLOCKS`] ranges of near-equal length.
fn cut(count: usize) -> Vec<(u32, u32)> {
    let pieces = count.min(BLOCKS);
    #[allow(clippy::cast_possible_truncation)]
    (0..pieces)
        .map(|k| {
            (
                (k * count / pieces) as u32,
                ((k + 1) * count / pieces) as u32,
            )
        })
        .collect()
}

/// The control indices spans `range` read along one direction.
fn controls(starts: &[u32], degree: usize, range: (u32, u32)) -> core::ops::RangeInclusive<usize> {
    starts[range.0 as usize] as usize - degree..=starts[range.1 as usize - 1] as usize
}

impl SpanBoxes {
    /// The spans of `spline`; `None` where a weight is not positive or a
    /// control point is not finite, so the hull says nothing.
    pub(super) fn of(spline: &BSplineSurface) -> Option<Self> {
        let usable = spline.grid().points().iter().all(|w| {
            w.weight > 0.0
                && w.weight.is_finite()
                && w.scaled.x.is_finite()
                && w.scaled.y.is_finite()
                && w.scaled.z.is_finite()
        });
        let starts = (span_starts(spline.u_knots()), span_starts(spline.v_knots()));
        if !usable || starts.0.is_empty() || starts.1.is_empty() {
            return None;
        }
        let net = spline.grid().points();
        let first = net[0].weight;
        let uniform = net.iter().all(|w| w.weight == first).then(|| 1.0 / first);
        let mut boxes = Self {
            starts,
            degree: (spline.u_knots().degree(), spline.v_knots().degree()),
            uniform,
            blocks: Vec::new(),
        };
        let patch = Patch::of(spline, uniform);
        let (along_u, along_v) = (cut(boxes.starts.0.len()), cut(boxes.starts.1.len()));
        for &u in &along_u {
            for &v in &along_v {
                boxes
                    .blocks
                    .push(((u, v), boxes.net_bounds(&patch, (u, v))));
            }
        }
        Some(boxes)
    }

    /// The box of the control points a block of spans reads.
    fn net_bounds(&self, patch: &Patch<'_>, (u, v): Block) -> Bounds {
        let cols = controls(&self.starts.1, self.degree.1, v);
        let mut bounds = Bounds::EMPTY;
        for i in controls(&self.starts.0, self.degree.0, u) {
            for j in cols.clone() {
                bounds.add(patch.cartesian(patch.homogeneous(i, j)));
            }
        }
        bounds
    }

    /// A lower bound on the distance from `target` to one span, from the
    /// control points it reads, as [`hull_gap`] takes it.
    fn span_gap(&self, patch: &Patch<'_>, (u, v): Block, target: [f64; 3]) -> f64 {
        let cols = controls(&self.starts.1, self.degree.1, v);
        let mut points: SmallVec<[[f64; 3]; 16]> = SmallVec::new();
        for i in controls(&self.starts.0, self.degree.0, u) {
            for j in cols.clone() {
                points.push(patch.cartesian(patch.homogeneous(i, j)));
            }
        }
        hull_gap(&points, self.degree.1 + 1, target)
    }

    /// The nearest point on `surface` to `target`: within [`SETTLE`] of
    /// the confusion of the least distance over the whole patch, then
    /// polished by Newton. `surface` is the B-spline patch this was built
    /// from.
    ///
    /// Two feet equally near to within that are both right; the one
    /// returned is fixed by the inputs. The seed is the nearest point the
    /// search met, one tying with it giving way to the one at the smaller
    /// `u`, then `v`, and the search visits its boxes nearest bound first,
    /// a tie to the one pushed first.
    pub(super) fn project(
        &self,
        surface: &SurfaceGeometry,
        target: Point,
        tol: Tolerances,
    ) -> OgeomResult<SurfaceProjection> {
        let SurfaceGeometry::BSpline(spline) = surface else {
            return refine_foot(surface, target, surface_start(surface), tol);
        };
        let mut search = Search {
            boxes: self,
            patch: Patch::of(spline, self.uniform),
            surface,
            tol,
            foot: None,
            target: [target.x, target.y, target.z],
            eps: tol.confusion() * SETTLE,
            tie: tol.confusion() * SETTLE,
            heap: BinaryHeap::new(),
            quarters: Vec::new(),
            patches: Vec::new(),
            serial: 0,
            best: (f64::INFINITY, (0.0, 0.0)),
        };
        for (k, (_, bounds)) in self.blocks.iter().enumerate() {
            #[allow(clippy::cast_possible_truncation)]
            search.push(bounds.gap(search.target), Item::Block(k as u32));
        }
        search.run(tol.parametric());
        match search.foot {
            Some(foot) if search.best == (foot.distance, foot.parameters) => Ok(foot),
            _ if search.best.0.is_finite() => refine_foot(surface, target, search.best.1, tol),
            _ => refine_foot(surface, target, surface_start(surface), tol),
        }
    }
}

/// The first corner of a surface's domain.
fn surface_start(surface: &SurfaceGeometry) -> (f64, f64) {
    use crate::Surface as _;
    let ((ua, _), (va, _)) = surface.domain();
    (ua, va)
}

/// What a search entry stands for.
#[derive(Debug, Clone, Copy)]
enum Item {
    /// One of the starting blocks.
    Block(u32),
    /// A block found by quartering, an index into the search's quarters.
    Quarter(u32),
    /// A Bezier patch, an index into the search's patches.
    Bezier(u32),
}

/// A search entry, nearest bound first, then first pushed.
#[derive(Debug, Clone, Copy)]
struct Entry {
    bound: f64,
    serial: u32,
    item: Item,
}

impl PartialEq for Entry {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Entry {}

impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Entry {
    /// Reversed, so the standard max-heap pops the nearest bound.
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .bound
            .total_cmp(&self.bound)
            .then(other.serial.cmp(&self.serial))
    }
}

/// A Bezier patch under search: its net and its parameter rectangle.
#[derive(Debug, Clone)]
struct Bezier {
    net: Net,
    u: (f64, f64),
    v: (f64, f64),
}

/// The search's state for one target.
struct Search<'a> {
    boxes: &'a SpanBoxes,
    patch: Patch<'a>,
    surface: &'a SurfaceGeometry,
    tol: Tolerances,
    /// The Newton polish from the nearest corner of the first patch
    /// opened: a foot early makes the nearest point so far exact, so the
    /// search prunes against the foot itself rather than a corner beside
    /// it.
    foot: Option<SurfaceProjection>,
    target: [f64; 3],
    eps: f64,
    /// Distances closer than this are a tie.
    tie: f64,
    heap: BinaryHeap<Entry>,
    /// Blocks found by quartering.
    quarters: Vec<Block>,
    patches: Vec<Bezier>,
    serial: u32,
    /// The nearest corner so far: its distance and parameters.
    best: (f64, (f64, f64)),
}

impl Search<'_> {
    fn push(&mut self, bound: f64, item: Item) {
        self.heap.push(Entry {
            bound,
            serial: self.serial,
            item,
        });
        self.serial = self.serial.wrapping_add(1);
    }

    /// Keep a point at distance `d` and parameters `at` if it is the
    /// nearest so far. Two within the tie of each other are equally near,
    /// and the one at the smaller `u`, then `v`, is kept: a point on a
    /// closed patch's seam is as near at either end of the domain, and
    /// rounding alone would otherwise pick the end.
    fn consider(&mut self, d: f64, at: (f64, f64)) {
        let (best, best_at) = self.best;
        if d < best - self.tie || (d <= best + self.tie && at < best_at) {
            self.best = (d, at);
        }
    }

    /// Whether a bound can still beat the nearest corner by more than the
    /// confusion.
    fn open(&self, bound: f64) -> bool {
        bound < self.best.0 - self.eps
    }

    fn run(&mut self, small: f64) {
        let mut opened = 0;
        while let Some(entry) = self.heap.pop() {
            if !self.open(entry.bound) || opened >= OPENED {
                break;
            }
            opened += 1;
            match entry.item {
                Item::Block(k) => self.quarter(self.boxes.blocks[k as usize].0),
                Item::Quarter(k) => self.quarter(self.quarters[k as usize]),
                Item::Bezier(k) => self.halve(k as usize, small),
            }
        }
    }

    /// A block of one span opens into its Bezier patch; a larger one is
    /// halved along each direction longer than a span.
    fn quarter(&mut self, (u, v): Block) {
        if u.1 - u.0 == 1 && v.1 - v.0 == 1 {
            let bezier = self.bezier(u.0 as usize, v.0 as usize);
            self.offer(bezier);
            return;
        }
        let halves = |(a, b): (u32, u32)| -> SmallVec<[(u32, u32); 2]> {
            if b - a > 1 {
                let mid = a + (b - a) / 2;
                smallvec::smallvec![(a, mid), (mid, b)]
            } else {
                smallvec::smallvec![(a, b)]
            }
        };
        for u in halves(u) {
            for v in halves(v) {
                let bound = if u.1 - u.0 == 1 && v.1 - v.0 == 1 {
                    self.boxes.span_gap(&self.patch, (u, v), self.target)
                } else {
                    self.boxes.net_bounds(&self.patch, (u, v)).gap(self.target)
                };
                if self.open(bound) {
                    #[allow(clippy::cast_possible_truncation)]
                    let k = self.quarters.len() as u32;
                    self.quarters.push((u, v));
                    self.push(bound, Item::Quarter(k));
                }
            }
        }
    }

    /// Span `(a, b)` as a Bezier patch: each column of its control points
    /// converted along `u`, then each row of the result along `v`.
    fn bezier(&self, a: usize, b: usize) -> Bezier {
        let patch = &self.patch;
        let (p, q) = patch.degree;
        let (su, sv) = (
            self.boxes.starts.0[a] as usize,
            self.boxes.starts.1[b] as usize,
        );
        let (along_u, along_v) = (
            blossom_ratios(patch.knots.0, p, su),
            blossom_ratios(patch.knots.1, q, sv),
        );
        let mut net = Net::with_capacity((p + 1) * (q + 1));
        for i in su - p..=su {
            for j in sv - q..=sv {
                net.push(patch.homogeneous(i, j));
            }
        }
        for j in 0..=q {
            to_bezier(&along_u, &mut net, p, |k| k * (q + 1) + j);
        }
        for i in 0..=p {
            to_bezier(&along_v, &mut net, q, |k| i * (q + 1) + k);
        }
        Bezier {
            net,
            u: (patch.knots.0[su], patch.knots.0[su + 1]),
            v: (patch.knots.1[sv], patch.knots.1[sv + 1]),
        }
    }

    /// Halve patch `k` both ways, unless it is already below the
    /// parametric tolerance both ways.
    fn halve(&mut self, k: usize, small: f64) {
        let (p, q) = self.patch.degree;
        let bezier = &self.patches[k];
        let (u, v) = (bezier.u, bezier.v);
        if u.1 - u.0 <= small && v.1 - v.0 <= small {
            return;
        }
        let (low, high) = split(&bezier.net, (p, q), true);
        let um = f64::midpoint(u.0, u.1);
        let vm = f64::midpoint(v.0, v.1);
        for (net, u) in [(low, (u.0, um)), (high, (um, u.1))] {
            let (left, right) = split(&net, (p, q), false);
            self.offer(Bezier {
                net: left,
                u,
                v: (v.0, vm),
            });
            self.offer(Bezier {
                net: right,
                u,
                v: (vm, v.1),
            });
        }
    }

    /// Take a patch's corners as candidates and queue it under its bound,
    /// unless the bound already rules it out.
    fn offer(&mut self, bezier: Bezier) {
        let (p, q) = self.patch.degree;
        let points: SmallVec<[[f64; 3]; 16]> = bezier
            .net
            .iter()
            .map(|h| self.patch.cartesian(*h))
            .collect();
        for (k, at) in [
            (0, (bezier.u.0, bezier.v.0)),
            (q, (bezier.u.0, bezier.v.1)),
            (p * (q + 1), (bezier.u.1, bezier.v.0)),
            (p * (q + 1) + q, (bezier.u.1, bezier.v.1)),
        ] {
            let off = sub(points[k], self.target);
            let d = dot(off, off).sqrt();
            self.consider(d, at);
        }
        if self.foot.is_none() {
            let [x, y, z] = self.target;
            if let Ok(foot) = refine_foot(self.surface, Point::new(x, y, z), self.best.1, self.tol)
            {
                self.consider(foot.distance, foot.parameters);
                self.foot = Some(foot);
            }
        }
        let bound = hull_gap(&points, q + 1, self.target);
        if self.open(bound) {
            #[allow(clippy::cast_possible_truncation)]
            let k = self.patches.len() as u32;
            self.patches.push(bezier);
            self.push(bound, Item::Bezier(k));
        }
    }
}

/// A lower bound on the distance from `target` to a surface piece inside
/// the hull of `points`, rows of `width`: the largest of its distances to
/// their axis-aligned box and to their boxes in two frames fitted to the
/// four corner points, one square to the piece's `u` edges and one to its
/// `v` edges. Any frame bounds the hull. The fitted ones hug a small,
/// nearly flat piece to its bending, sideways too: a foot on a piece's
/// edge is bounded as closely as one inside it.
fn hull_gap(points: &[[f64; 3]], width: usize, target: [f64; 3]) -> f64 {
    let mut aabb = Bounds::EMPTY;
    for &point in points {
        aabb.add(point);
    }
    let mut bound = aabb.gap(target);
    let last = points.len() - 1;
    let (c00, c01) = (points[0], points[width - 1]);
    let (c10, c11) = (points[last + 1 - width], points[last]);
    let du: [f64; 3] = core::array::from_fn(|k| (c10[k] - c00[k]) + (c11[k] - c01[k]));
    let dv: [f64; 3] = core::array::from_fn(|k| (c01[k] - c00[k]) + (c11[k] - c10[k]));
    let normal = cross(du, dv);
    let (nn, uu, vv) = (dot(normal, normal), dot(du, du), dot(dv, dv));
    if uu > 0.0 && vv > 0.0 && nn > 1e-24 * uu * vv {
        let n = scaled(normal, 1.0 / nn.sqrt());
        // Frame one runs along `du`, frame two along `dv`; both share
        // the normal.
        let (a1, v2) = (scaled(du, 1.0 / uu.sqrt()), scaled(dv, 1.0 / vv.sqrt()));
        let (b1, a2) = (cross(n, a1), cross(v2, n));
        let coordinates = |p: [f64; 3]| {
            let off = sub(p, c00);
            (
                [dot(off, a1), dot(off, b1), dot(off, n)],
                [dot(off, a2), dot(off, v2)],
            )
        };
        let (mut one, mut two) = (Bounds::EMPTY, Bounds::EMPTY);
        for &point in points {
            let (first, second) = coordinates(point);
            one.add(first);
            two.add([second[0], second[1], first[2]]);
        }
        let (first, second) = coordinates(target);
        bound = bound
            .max(one.gap(first))
            .max(two.gap([second[0], second[1], first[2]]));
    }
    bound
}

/// The ratios de Boor's recurrence blends by to turn the span starting at
/// knot `s` into Bezier form: Bezier point `j` is the blossom at the
/// span's start repeated `degree - j` times and its end `j` times, a
/// parameter per level of the recurrence. They depend on the knots alone,
/// so one table serves every line of a patch.
fn blossom_ratios(knots: &[f64], degree: usize, s: usize) -> SmallVec<[f64; 32]> {
    let p = degree;
    let (a, b) = (knots[s], knots[s + 1]);
    let mut out = SmallVec::new();
    for j in 0..=p {
        for r in 1..=p {
            let t = if r <= p - j { a } else { b };
            for i in (r..=p).rev() {
                let left = knots[s - p + i];
                let right = knots[s + 1 + i - r];
                out.push((t - left) / (right - left));
            }
        }
    }
    out
}

/// Turn one line of a net, the points at `at(0..=degree)`, from the
/// control points a span reads into its Bezier points, by the recurrence
/// [`blossom_ratios`] tabulates.
fn to_bezier(ratios: &[f64], net: &mut Net, degree: usize, at: impl Fn(usize) -> usize) {
    let p = degree;
    let source: Line = (0..=p).map(&at).map(|k| net[k]).collect();
    let mut work = source.clone();
    let mut ratio = ratios.iter().copied();
    for j in 0..=p {
        work.copy_from_slice(&source);
        for r in 1..=p {
            for i in (r..=p).rev() {
                let t = ratio.next().unwrap_or(0.0);
                work[i] = lerp(work[i - 1], work[i], t);
            }
        }
        net[at(j)] = work[p];
    }
}

/// The point a fraction `t` of the way from `a` to `b`. Plain products,
/// not fused: without a target feature for it a fused multiply-add is a
/// library call, and the search blends hundreds of points per target.
fn lerp(a: Homogeneous, b: Homogeneous, t: f64) -> Homogeneous {
    [
        (b[0] - a[0]) * t + a[0],
        (b[1] - a[1]) * t + a[1],
        (b[2] - a[2]) * t + a[2],
        (b[3] - a[3]) * t + a[3],
    ]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn scaled(a: [f64; 3], k: f64) -> [f64; 3] {
    [a[0] * k, a[1] * k, a[2] * k]
}

/// A Bezier net of `p + 1` rows of `q + 1` points halved at the middle of
/// `u` (`along_u`) or of `v`, by de Casteljau down each line of the net.
fn split(net: &Net, (p, q): (usize, usize), along_u: bool) -> (Net, Net) {
    let (degree, lines) = if along_u { (p, q + 1) } else { (q, p + 1) };
    let at = |line: usize, k: usize| {
        if along_u {
            k * (q + 1) + line
        } else {
            line * (q + 1) + k
        }
    };
    let mut low = net.clone();
    let mut high = net.clone();
    let mut work: Line = SmallVec::from_elem([0.0; 4], degree + 1);
    for line in 0..lines {
        for (k, slot) in work.iter_mut().enumerate() {
            *slot = net[at(line, k)];
        }
        for r in 1..=degree {
            for i in 0..=degree - r {
                work[i] = lerp(work[i], work[i + 1], 0.5);
            }
            // Level `r` gives the low half its point `r` and the high half
            // its point `degree - r`.
            low[at(line, r)] = work[0];
            high[at(line, degree - r)] = work[degree - r];
        }
    }
    (low, high)
}
