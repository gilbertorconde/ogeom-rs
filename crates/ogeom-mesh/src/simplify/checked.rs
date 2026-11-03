//! Decimation that checks each collapse before making it.
//!
//! [`Target::Verified`](super::Target::Verified) measures its result against
//! the source both ways. A collapse is made here only where that
//! measurement, taken of the collapse's own neighbourhood, would still pass,
//! and only where the triangles it leaves meet no other triangle but at the
//! corners they share.
//!
//! # What a collapse can change
//!
//! Merging `from` into `to` at a point moves `to` and replaces `from`, so
//! the triangles that change are exactly those about either, and every other
//! triangle and vertex of the result is unchanged. The measurement points
//! of the result that change are therefore the merged vertex and the
//! centres of the triangles about it afterwards, its fan, and each is
//! measured to the source before the collapse is made.
//!
//! Each source vertex is held by one live triangle of the result it was
//! found within tolerance of, so its distance to the result is at most its
//! distance to that triangle. When the triangle changes, the vertex stays
//! with it if it is still near enough, and otherwise goes to a triangle of
//! the new fan that is, the fan covering what the old triangles about the
//! edge covered; a collapse that would leave it farther than the tolerance
//! from all of them is not made. The distance is the one a nearest query
//! over the result takes for that triangle, so the measurement afterwards
//! finds the vertex no farther.
//!
//! # Crossings
//!
//! Before a collapse, each triangle of the new fan is tested exactly, as
//! [`boolean`](crate::boolean()) tests an operand, against every other
//! triangle whose box meets its own, found through a grid of the current
//! triangles kept up to date collapse by collapse. A mesh that meets
//! itself only at shared corners before the collapse does so after it.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use ogeom_core::{Exact, FastMap, FastSet, Predicates};
use ogeom_math::Point;
use ogeom_topo::Triangulation;

use super::State;
use crate::Bvh;
use crate::boolean::apart_but_at_shared;
use crate::bvh::distance2_to;

/// A collapse of `from` into `to` at `at`, and the triangles it leaves.
struct Collapse {
    from: u32,
    to: u32,
    at: Point,
    /// The triangles about the edge that stay, with their vertices and
    /// corners as they will be.
    fan: Vec<(u32, [u32; 3], [Point; 3])>,
    /// The triangles about the edge that vanish.
    gone: Vec<u32>,
}

impl Collapse {
    fn of(state: &State, from: u32, to: u32, at: Point) -> Self {
        let mut fan = Vec::new();
        let mut gone = Vec::new();
        for v in [to, from] {
            for &index in &state.around[v as usize] {
                if !state.live[index as usize] {
                    continue;
                }
                let triangle = state.triangles[index as usize];
                if triangle.contains(&from) && triangle.contains(&to) {
                    if v == to {
                        gone.push(index);
                    }
                    continue;
                }
                let vertices = triangle.map(|w| if w == from { to } else { w });
                let corners = vertices.map(|w| {
                    if w == to {
                        at
                    } else {
                        state.positions[w as usize]
                    }
                });
                fan.push((index, vertices, corners));
            }
        }
        Self {
            from,
            to,
            at,
            fan,
            gone,
        }
    }

    fn changes(&self, triangle: u32) -> bool {
        self.gone.contains(&triangle) || self.fan.iter().any(|f| f.0 == triangle)
    }
}

/// The centre of a triangle, as the measurement takes it.
pub(super) fn centre_of([a, b, c]: [Point; 3]) -> Point {
    Point::from_vector((a.to_vector() + b.to_vector() + c.to_vector()) * (1.0 / 3.0))
}

/// What keeps each collapse within the measurement's tolerance.
pub(super) struct Watch<'a> {
    source: &'a Bvh,
    points: &'a [Point],
    tolerance: f64,
    /// For each triangle, the source vertices it holds.
    owned: Vec<Vec<u32>>,
}

impl<'a> Watch<'a> {
    /// Each source vertex held by a triangle it is a corner of.
    pub(super) fn new(source: &'a Bvh, mesh: &'a Triangulation, tolerance: f64) -> Self {
        let mut held = vec![false; mesh.positions.len()];
        let mut owned = vec![Vec::new(); mesh.triangles.len()];
        for (index, triangle) in mesh.triangles.iter().enumerate() {
            for &v in triangle {
                if !held[v as usize] {
                    held[v as usize] = true;
                    owned[index].push(v);
                }
            }
        }
        Self {
            source,
            points: &mesh.positions,
            tolerance,
            owned,
        }
    }

    fn near_source(&self, p: Point) -> bool {
        self.source
            .nearest(p, self.tolerance)
            .is_some_and(|hit| hit.distance <= self.tolerance)
    }

    /// Where each source vertex held by a triangle the collapse changes
    /// would go, or `None` if one would be past the tolerance, or the merged
    /// vertex or a centre of the new fan would be.
    fn allows(&self, collapse: &Collapse) -> Option<Vec<(u32, u32)>> {
        let mut moves = Vec::new();
        let within = |p: Point, f: &(u32, [u32; 3], [Point; 3])| {
            distance2_to(p, f.2).sqrt() <= self.tolerance
        };
        for &t in &collapse.gone {
            for &s in &self.owned[t as usize] {
                let p = self.points[s as usize];
                moves.push((s, collapse.fan.iter().find(|f| within(p, f))?.0));
            }
        }
        for (k, f) in collapse.fan.iter().enumerate() {
            for &s in &self.owned[f.0 as usize] {
                let p = self.points[s as usize];
                // Kept by its own triangle where that is near enough, so
                // most stay where they were.
                let at = if within(p, f) {
                    f.0
                } else {
                    collapse
                        .fan
                        .iter()
                        .enumerate()
                        .find(|&(j, g)| j != k && within(p, g))?
                        .1
                        .0
                };
                moves.push((s, at));
            }
        }
        if !self.near_source(collapse.at) {
            return None;
        }
        for f in &collapse.fan {
            if !self.near_source(centre_of(f.2)) {
                return None;
            }
        }
        Some(moves)
    }

    fn apply(&mut self, collapse: &Collapse, moves: Vec<(u32, u32)>) {
        for &t in collapse
            .gone
            .iter()
            .chain(collapse.fan.iter().map(|f| &f.0))
        {
            self.owned[t as usize].clear();
        }
        for (s, t) in moves {
            self.owned[t as usize].push(s);
        }
    }
}

type Box3 = ([f64; 3], [f64; 3]);

fn box_of(corners: &[Point; 3]) -> Box3 {
    let mut low = [corners[0].x, corners[0].y, corners[0].z];
    let mut high = low;
    for p in &corners[1..] {
        for (k, c) in [p.x, p.y, p.z].into_iter().enumerate() {
            low[k] = low[k].min(c);
            high[k] = high[k].max(c);
        }
    }
    (low, high)
}

/// Whether two triangles meet only at the vertices they share, decided by
/// the side tests [`apart_but_at_shared`] begins with where they settle it,
/// and by it where they do not.
///
/// Sharing an edge, the triangles meet only along it when the fourth vertex
/// is off the plane of the other three; sharing a corner, when one
/// triangle's other two corners are strictly on one side of the other's
/// plane; sharing nothing, when one triangle's corners are all strictly on
/// one side of the other's plane.
fn apart((av, ap): ([u32; 3], [Point; 3]), (bv, bp): ([u32; 3], [Point; 3])) -> bool {
    let xyz = |p: Point| [p.x, p.y, p.z];
    let side = |t: &[Point; 3], q: Point| Exact::orient3d(xyz(t[0]), xyz(t[1]), xyz(t[2]), xyz(q));
    let shared = |v: &[u32; 3], k: usize, w: &[u32; 3]| w.contains(&v[k]);
    // The corners of `y` not shared with `x`, all strictly on one side of
    // `x`'s plane.
    let off = |x: &[Point; 3], yv: &[u32; 3], yp: &[Point; 3], xv: &[u32; 3]| {
        let mut sign = None;
        for (v, &q) in yv.iter().zip(yp) {
            if xv.contains(v) {
                continue;
            }
            let s = side(x, q);
            if s.is_zero() || sign.is_some_and(|t| t != s) {
                return false;
            }
            sign = Some(s);
        }
        true
    };
    let count = (0..3).filter(|&k| shared(&av, k, &bv)).count();
    let settled = match count {
        2 => off(&ap, &bv, &bp, &av),
        0 | 1 => off(&ap, &bv, &bp, &av) || off(&bp, &av, &ap, &bv),
        _ => false,
    };
    settled || apart_but_at_shared((av, ap), (bv, bp))
}

fn boxes_meet(a: &Box3, b: &Box3) -> bool {
    (0..3).all(|k| a.0[k] <= b.1[k] && b.0[k] <= a.1[k])
}

/// The current triangles, filed by where their boxes are.
///
/// A loose grid at several levels, each cell a power of two wide and twice
/// as wide as the last. A triangle is filed once, at the finest level whose
/// cells are more than twice as wide as its box, in the cell holding its
/// box's low corner, so its box reaches at most half a cell past that
/// cell's high side. The cell size sets only how many triangles a query
/// looks at; the answer is the same at any.
pub(super) struct Cells {
    base: f64,
    /// Each level's cells, by the slot of [`Cells::filed`] holding each.
    levels: Vec<FastMap<[i64; 3], u32>>,
    /// The triangles filed in each cell, with their boxes.
    filed: Vec<Vec<(u32, Box3)>>,
    /// Slots of [`Cells::filed`] left empty, to reuse.
    spare: Vec<u32>,
    /// Where each live triangle is filed: its level, cell, the cell's slot
    /// and its place there.
    home: Vec<Home>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct Home {
    level: usize,
    at: [i64; 3],
    slot: u32,
    place: u32,
}

impl Cells {
    /// The grid over the state's live triangles, its finest cells the
    /// first power of two more than twice their mean edge, so that a
    /// triangle of that edge is filed at the finest level.
    pub(super) fn over(state: &State) -> Self {
        let mut length = 0.0;
        let mut count = 0_u32;
        for (index, triangle) in state.triangles.iter().enumerate() {
            if state.live[index] {
                for k in 0..3 {
                    let (a, b) = (triangle[k], triangle[(k + 1) % 3]);
                    length +=
                        (state.positions[a as usize] - state.positions[b as usize]).magnitude();
                    count += 1;
                }
            }
        }
        let mean = if count > 0 {
            length / f64::from(count)
        } else {
            0.0
        };
        let mut base = 1.0_f64;
        if mean.is_finite() && mean > 0.0 {
            while base <= 2.0 * mean && base.is_finite() {
                base *= 2.0;
            }
            while base / 2.0 > 2.0 * mean && base / 2.0 > f64::MIN_POSITIVE {
                base /= 2.0;
            }
        }
        let mut cells = Self {
            base,
            levels: Vec::new(),
            filed: Vec::new(),
            spare: Vec::new(),
            home: vec![Home::default(); state.triangles.len()],
        };
        for (index, triangle) in state.triangles.iter().enumerate() {
            if state.live[index] {
                let corners = triangle.map(|v| state.positions[v as usize]);
                let (level, at) = cells.home_of(&box_of(&corners));
                #[allow(clippy::cast_possible_truncation)]
                cells.add(index as u32, level, at, box_of(&corners));
            }
        }
        cells
    }

    /// The width of a level's cells, a power of two, so that a coordinate
    /// divides by it exactly.
    fn size(&self, level: usize) -> f64 {
        #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
        let level = level as i32;
        self.base * 2.0_f64.powi(level)
    }

    fn cell(&self, level: usize, p: [f64; 3]) -> [i64; 3] {
        let size = self.size(level);
        #[allow(clippy::cast_possible_truncation)]
        p.map(|c| (c / size).floor() as i64)
    }

    /// The level and cell a box is filed in. An extent that rounds to less
    /// than half a cell is less than half a cell.
    fn home_of(&self, b: &Box3) -> (usize, [i64; 3]) {
        let extent = (0..3).map(|k| b.1[k] - b.0[k]).fold(0.0, f64::max);
        let mut level = 0;
        while extent >= self.size(level) / 2.0 && self.size(level).is_finite() {
            level += 1;
        }
        (level, self.cell(level, b.0))
    }

    fn add(&mut self, triangle: u32, level: usize, at: [i64; 3], b: Box3) {
        while self.levels.len() <= level {
            self.levels.push(FastMap::default());
        }
        let (filed, spare) = (&mut self.filed, &mut self.spare);
        let slot = *self.levels[level].entry(at).or_insert_with(|| {
            spare.pop().unwrap_or_else(|| {
                filed.push(Vec::new());
                u32::try_from(filed.len() - 1).unwrap_or(u32::MAX)
            })
        });
        let cell = &mut self.filed[slot as usize];
        cell.push((triangle, b));
        self.home[triangle as usize] = Home {
            level,
            at,
            slot,
            place: u32::try_from(cell.len() - 1).unwrap_or(u32::MAX),
        };
    }

    fn unfile(&mut self, triangle: u32) {
        let home = self.home[triangle as usize];
        let cell = &mut self.filed[home.slot as usize];
        cell.swap_remove(home.place as usize);
        if let Some(&(moved, _)) = cell.get(home.place as usize) {
            self.home[moved as usize].place = home.place;
        }
        if cell.is_empty() {
            self.levels[home.level].remove(&home.at);
            self.spare.push(home.slot);
        }
    }

    fn refile_one(&mut self, triangle: u32, corners: &[Point; 3]) {
        let b = box_of(corners);
        let (level, at) = self.home_of(&b);
        let home = self.home[triangle as usize];
        if (home.level, home.at) == (level, at) {
            self.filed[home.slot as usize][home.place as usize].1 = b;
            return;
        }
        self.unfile(triangle);
        self.add(triangle, level, at, b);
    }

    /// The live triangles whose boxes meet `query`.
    fn near(&self, query: &Box3) -> Vec<u32> {
        let mut found = Vec::new();
        for (level, map) in self.levels.iter().enumerate() {
            if map.is_empty() {
                continue;
            }
            // A box meeting the query has its low corner at most the query's
            // high side, and more than half a cell below the query's low
            // side; with exact division, in these cells.
            let low = self.cell(level, query.0).map(|c| c.saturating_sub(1));
            let high = self.cell(level, query.1);
            let cells: i128 = (0..3).map(|k| i128::from(high[k] - low[k]) + 1).product();
            let take = |slot: u32, found: &mut Vec<u32>| {
                let filed = &self.filed[slot as usize];
                found.extend(
                    filed
                        .iter()
                        .filter(|f| boxes_meet(query, &f.1))
                        .map(|f| f.0),
                );
            };
            #[allow(clippy::cast_possible_wrap)]
            if cells > map.len() as i128 {
                for (at, &slot) in map {
                    if (0..3).all(|k| low[k] <= at[k] && at[k] <= high[k]) {
                        take(slot, &mut found);
                    }
                }
            } else {
                for i in low[0]..=high[0] {
                    for j in low[1]..=high[1] {
                        for k in low[2]..=high[2] {
                            if let Some(&slot) = map.get(&[i, j, k]) {
                                take(slot, &mut found);
                            }
                        }
                    }
                }
            }
        }
        found
    }

    /// Whether a triangle of the collapse's fan would meet another triangle
    /// of the state anywhere but at the corners they share.
    fn crosses(&self, state: &State, collapse: &Collapse) -> bool {
        let fan = &collapse.fan;
        for (i, a) in fan.iter().enumerate() {
            for b in &fan[i + 1..] {
                if !apart((a.1, a.2), (b.1, b.2)) {
                    return true;
                }
            }
        }
        let Some(all) = collapse.bounds() else {
            return false;
        };
        let boxes: Vec<Box3> = fan.iter().map(|f| box_of(&f.2)).collect();
        for u in self.near(&all) {
            if collapse.changes(u) {
                continue;
            }
            let vertices = state.triangles[u as usize];
            let corners = vertices.map(|v| state.positions[v as usize]);
            let b = box_of(&corners);
            for (f, fb) in fan.iter().zip(&boxes) {
                if boxes_meet(fb, &b) && !apart((f.1, f.2), (vertices, corners)) {
                    return true;
                }
            }
        }
        false
    }

    /// File the collapse's fan with its new corners, and drop what it
    /// removed.
    fn refile(&mut self, collapse: &Collapse) {
        for &t in &collapse.gone {
            self.unfile(t);
        }
        for f in &collapse.fan {
            self.refile_one(f.0, &f.2);
        }
    }
}

impl Collapse {
    /// The box of the fan, or `None` for no fan.
    fn bounds(&self) -> Option<Box3> {
        let mut boxes = self.fan.iter().map(|f| box_of(&f.2));
        let mut all = boxes.next()?;
        for b in boxes {
            for k in 0..3 {
                all.0[k] = all.0[k].min(b.0[k]);
                all.1[k] = all.1[k].max(b.1[k]);
            }
        }
        Some(all)
    }

    /// Whether the fan meets another collapse's fan anywhere but at the
    /// corners they share.
    fn crosses_fan(&self, other: &Self) -> bool {
        for f in &self.fan {
            let fb = box_of(&f.2);
            for g in &other.fan {
                if boxes_meet(&fb, &box_of(&g.2)) && !apart((f.1, f.2), (g.1, g.2)) {
                    return true;
                }
            }
        }
        false
    }
}

/// How many collapses are weighed at once.
///
/// Collapses whose neighbourhoods share no vertex change disjoint sets of
/// triangles, so each is weighed against the mesh as it was before any of
/// them, on the kernel's workers, and they are then made in order. A fixed
/// number rather than one per worker keeps the result the same on any
/// machine.
const BATCH: usize = 512;

/// A collapse weighed, and where the source vertices it moves would go.
struct Weighed {
    collapse: Collapse,
    moves: Vec<(u32, u32)>,
}

/// Weigh collapsing `hi` into `lo`: `None` if it would fold a triangle,
/// leave the mesh other than a manifold, move a measured point past the
/// watch's tolerance, or make a triangle meet another away from their
/// shared corners.
fn weigh(
    state: &State,
    cells: &Cells,
    watch: Option<&Watch<'_>>,
    (lo, hi): (u32, u32),
) -> Option<Weighed> {
    let (_, at) = state.mean_cost(lo, hi);
    if state.would_fold(hi, lo, at) || !state.keeps_manifold(lo, hi) {
        return None;
    }
    let collapse = Collapse::of(state, hi, lo, at);
    let moves = match watch {
        Some(w) => w.allows(&collapse)?,
        None => Vec::new(),
    };
    if cells.crosses(state, &collapse) {
        return None;
    }
    Some(Weighed { collapse, moves })
}

/// Collapse the cheapest free edge, by [`State::mean_cost`], while its cost
/// is within `budget` and more than `floor` triangles remain, making only
/// collapses that leave the mesh meeting itself only at shared corners and,
/// with a [`Watch`], within its tolerance; returns how many were made.
///
/// Edges are taken cheapest first, ties to the lower vertex numbers, from a
/// queue, up to [`BATCH`] at a time whose neighbourhoods share no vertex;
/// an edge whose neighbourhood meets one already taken waits for the next
/// batch. A collapse changes the cost of the edges at the merged vertex,
/// which are queued again, the entries before marked stale; and whether the
/// edges about its neighbours pass, so an edge refused is set aside at both
/// ends and queued again when a collapse next changes the triangles about
/// either.
pub(super) fn decimate(
    state: &mut State,
    budget: f64,
    floor: usize,
    cells: &mut Cells,
    mut watch: Option<&mut Watch<'_>>,
) -> usize {
    let mut stamps = vec![0_u32; state.positions.len()];
    // The edges refused, each listed at both ends until queued again.
    let mut parked: Vec<Vec<u32>> = vec![Vec::new(); state.positions.len()];
    let mut refused: FastSet<(u32, u32)> = FastSet::default();
    let mut queue = BinaryHeap::new();
    let mut seen: FastSet<(u32, u32)> = FastSet::default();
    for (i, triangle) in state.triangles.iter().enumerate() {
        if !state.live[i] {
            continue;
        }
        for k in 0..3 {
            let (a, b) = (triangle[k], triangle[(k + 1) % 3]);
            let key = (a.min(b), a.max(b));
            if seen.insert(key) && state.free(a, b) {
                queue.push(queued(state, key, &stamps));
            }
        }
    }
    drop(seen);

    let alive = |state: &State, v: u32| !state.around[v as usize].is_empty();
    // The batch each vertex was last claimed for.
    let mut claimed = vec![0_u32; state.positions.len()];
    let mut batch_number = 0_u32;
    let mut collapsed = 0;
    while state.remaining > floor {
        batch_number += 1;
        let mut batch: Vec<(u32, u32)> = Vec::new();
        let mut waiting = Vec::new();
        while batch.len() < BATCH {
            let Some(Reverse(entry)) = queue.pop() else {
                break;
            };
            if entry.cost > budget {
                // The batch may queue cheaper edges again.
                waiting.push(Reverse(entry));
                break;
            }
            let (lo, hi) = entry.edge;
            if entry.stamps != (stamps[lo as usize], stamps[hi as usize])
                || !alive(state, lo)
                || !alive(state, hi)
            {
                continue;
            }
            let near = || {
                [lo, hi]
                    .into_iter()
                    .flat_map(|v| &state.around[v as usize])
                    .filter(|&&t| state.live[t as usize])
                    .flat_map(|&t| state.triangles[t as usize])
            };
            if near().any(|v| claimed[v as usize] == batch_number) {
                waiting.push(Reverse(entry));
                // Once the cheapest left meets the batch, what follows
                // would mostly wait as well.
                if waiting.len() >= BATCH {
                    break;
                }
                continue;
            }
            for v in near() {
                claimed[v as usize] = batch_number;
            }
            batch.push((lo, hi));
        }
        // An edge waits only on one taken, so with none taken every edge
        // left is past the budget.
        if batch.is_empty() {
            break;
        }

        let weighed = {
            let (state, cells, watch) = (&*state, &*cells, watch.as_deref());
            ogeom_core::parallel::map_ordered(&batch, |_, &edge| weigh(state, cells, watch, edge))
        };
        let mut made: Vec<(Box3, Collapse)> = Vec::new();
        for (&(lo, hi), weighed) in batch.iter().zip(weighed) {
            if state.remaining <= floor {
                break;
            }
            // Weighed against the mesh before the batch; the collapses made
            // since changed only their own fans.
            let passed = weighed.filter(|w| {
                let mine = w.collapse.bounds();
                !made.iter().any(|(bounds, other)| {
                    mine.is_some_and(|b| boxes_meet(bounds, &b)) && w.collapse.crosses_fan(other)
                })
            });
            let Some(Weighed { collapse, moves }) = passed else {
                if refused.insert((lo, hi)) {
                    parked[lo as usize].push(hi);
                    parked[hi as usize].push(lo);
                }
                continue;
            };
            state.collapse(collapse.from, collapse.to, collapse.at);
            if let Some(w) = watch.as_deref_mut() {
                w.apply(&collapse, moves);
            }
            cells.refile(&collapse);
            collapsed += 1;

            stamps[lo as usize] = stamps[lo as usize].wrapping_add(1);
            for w in std::mem::take(&mut parked[hi as usize]) {
                refused.remove(&(hi.min(w), hi.max(w)));
            }
            for w in std::mem::take(&mut parked[lo as usize]) {
                refused.remove(&(lo.min(w), lo.max(w)));
            }
            let ring = state.ring(lo);
            for &w in &ring {
                if state.free(lo, w) {
                    queue.push(queued(state, (lo.min(w), lo.max(w)), &stamps));
                }
            }
            for &v in &ring {
                for w in std::mem::take(&mut parked[v as usize]) {
                    let edge = (v.min(w), v.max(w));
                    if refused.remove(&edge) && w != lo && alive(state, w) && state.free(v, w) {
                        queue.push(queued(state, edge, &stamps));
                    }
                }
            }
            if let Some(bounds) = collapse.bounds() {
                made.push((bounds, collapse));
            }
        }
        queue.extend(waiting);
        // A closed mesh has three edges for every two triangles, so a queue
        // past two entries a triangle is at least a quarter stale: those
        // are dropped at once rather than popped one by one.
        if queue.len() > 2 * state.remaining {
            let mut entries = std::mem::take(&mut queue).into_vec();
            entries.retain(|Reverse(entry)| {
                let (lo, hi) = entry.edge;
                entry.stamps == (stamps[lo as usize], stamps[hi as usize])
                    && alive(state, lo)
                    && alive(state, hi)
            });
            queue = BinaryHeap::from(entries);
        }
    }
    collapsed
}

fn queued(state: &State, edge: (u32, u32), stamps: &[u32]) -> Reverse<Queued> {
    Reverse(Queued {
        cost: state.mean_cost(edge.0, edge.1).0,
        edge,
        stamps: (stamps[edge.0 as usize], stamps[edge.1 as usize]),
    })
}

/// An edge waiting to collapse, with its cost and its ends' stamps when it
/// was queued.
#[derive(Debug, Clone, Copy)]
struct Queued {
    cost: f64,
    edge: (u32, u32),
    stamps: (u32, u32),
}

impl Queued {
    fn key(&self) -> (f64, (u32, u32), (u32, u32)) {
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
