//! A k-d tree over a point cloud: the nearest point, the `k` nearest, and
//! every point within a ball.
//!
//! The same contract as [`Bvh`](crate::bvh::Bvh): built once, immutable,
//! shared across threads, and every answer the one brute force gives over
//! the same points with the same arithmetic, ties going to the lowest point
//! index. A cell is entered only if its box's computed squared distance is no
//! more than the best found so far, and since every axis gap to a box is at
//! most the gap to a point in it, rounded the same way, no point a query
//! would accept lies in a cell it skips.

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::Point;

use crate::bvh::{distance2, gap};

type P = [f64; 3];

/// A point a query found, by its index in the cloud.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Neighbour {
    /// The point's index in the cloud the index was built from.
    pub index: u32,
    /// Its distance from the query point.
    pub distance: f64,
}

#[derive(Debug, Clone, Copy)]
struct Cell {
    low: P,
    high: P,
    /// A leaf's first slot, or an interior cell's upper half. Its lower
    /// half is the next cell.
    first: u32,
    /// How many points a leaf holds; zero for an interior cell.
    count: u32,
}

/// A k-d tree over points, split at the median of each cell's widest axis.
///
/// Holds its own copy of the points.
#[derive(Debug, Clone)]
pub struct PointIndex {
    cells: Vec<Cell>,
    /// The points, slots in leaf order.
    points: Vec<P>,
    /// The point in each slot.
    order: Vec<u32>,
}

/// The most points a leaf holds.
const LEAF_MOST: usize = 8;

impl PointIndex {
    /// Build the tree over a cloud.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if a point
    /// is not finite, and
    /// [`OgeomError::Range`](ogeom_core::OgeomError::Range) for a cloud of
    /// more than four billion points.
    pub fn build(points: &[Point]) -> OgeomResult<Self> {
        if u32::try_from(points.len()).is_err() {
            ogeom_bail!(Range, "a cloud of more than four billion points");
        }
        let mut held = Vec::with_capacity(points.len());
        for (i, p) in points.iter().enumerate() {
            if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
                ogeom_bail!(Numeric, "point {i} is not finite");
            }
            held.push([p.x, p.y, p.z]);
        }
        let mut order: Vec<u32> = (0..held.len())
            .map(|i| u32::try_from(i).unwrap_or(u32::MAX))
            .collect();
        let cells = build_cells(&held, &mut order);
        Ok(Self {
            cells,
            points: order.iter().map(|&i| held[i as usize]).collect(),
            order,
        })
    }

    /// How many points the tree holds.
    #[must_use]
    pub fn len(&self) -> usize {
        self.order.len()
    }

    /// Whether the tree holds no point.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.order.is_empty()
    }

    /// The point nearest `p`, no farther than `max`; the lowest index among
    /// equally near ones. `max` may be infinite.
    #[must_use]
    pub fn nearest(&self, p: Point, max: f64) -> Option<Neighbour> {
        self.k_nearest_with(p, 1, max, &mut Vec::new()).pop()
    }

    /// [`PointIndex::nearest`] for each point, in parallel, answers in point
    /// order.
    #[must_use]
    pub fn nearest_many(&self, points: &[Point], max: f64) -> Vec<Option<Neighbour>> {
        crate::bvh::batched(
            points,
            |&p| p,
            |&p, stack| self.k_nearest_with(p, 1, max, stack).pop(),
        )
    }

    /// The `k` points nearest `p` no farther than `max`, nearest first, ties
    /// in index order: fewer where fewer are that near.
    #[must_use]
    pub fn k_nearest(&self, p: Point, k: usize, max: f64) -> Vec<Neighbour> {
        self.k_nearest_with(p, k, max, &mut Vec::new())
    }

    fn k_nearest_with(
        &self,
        p: Point,
        k: usize,
        max: f64,
        stack: &mut Vec<(u32, f64)>,
    ) -> Vec<Neighbour> {
        let p = [p.x, p.y, p.z];
        let reach = max * max;
        if k == 0
            || self.cells.is_empty()
            || !p.iter().all(|c| c.is_finite())
            || max.is_nan()
            || max < 0.0
        {
            return Vec::new();
        }
        // The best so far as (squared distance, index), sorted, at most `k`.
        let mut best: Vec<(f64, u32)> = Vec::with_capacity(k.min(self.order.len()) + 1);
        let bound = |best: &Vec<(f64, u32)>| {
            if best.len() == k {
                best[k - 1].0
            } else {
                reach
            }
        };
        stack.clear();
        stack.push((0, cell_gap(&self.cells[0], p)));
        while let Some((at, g)) = stack.pop() {
            if g > bound(&best) {
                continue;
            }
            let cell = &self.cells[at as usize];
            if cell.count > 0 {
                let first = cell.first as usize;
                for slot in first..first + cell.count as usize {
                    let d2 = distance2(self.points[slot], p);
                    let index = self.order[slot];
                    let full = best.len() == k;
                    if d2 > reach || (full && (d2, index) >= best[k - 1]) {
                        continue;
                    }
                    let at = best.partition_point(|&held| held < (d2, index));
                    best.insert(at, (d2, index));
                    best.truncate(k);
                }
                continue;
            }
            let (lower, upper) = (at + 1, cell.first);
            let gl = cell_gap(&self.cells[lower as usize], p);
            let gu = cell_gap(&self.cells[upper as usize], p);
            let (near, far) = if gl <= gu {
                ((lower, gl), (upper, gu))
            } else {
                ((upper, gu), (lower, gl))
            };
            let limit = bound(&best);
            if far.1 <= limit {
                stack.push(far);
            }
            if near.1 <= limit {
                stack.push(near);
            }
        }
        best.into_iter()
            .map(|(d2, index)| Neighbour {
                index,
                distance: d2.sqrt(),
            })
            .collect()
    }

    /// Every point within `radius` of `centre`, in index order.
    #[must_use]
    pub fn within_ball(&self, centre: Point, radius: f64) -> Vec<u32> {
        let c = [centre.x, centre.y, centre.z];
        let reach = radius * radius;
        let mut found = Vec::new();
        if self.cells.is_empty()
            || !c.iter().all(|v| v.is_finite())
            || radius.is_nan()
            || radius < 0.0
        {
            return found;
        }
        let mut stack = vec![0_u32];
        while let Some(at) = stack.pop() {
            let cell = &self.cells[at as usize];
            if cell_gap(cell, c) > reach {
                continue;
            }
            if cell.count > 0 {
                let first = cell.first as usize;
                for slot in first..first + cell.count as usize {
                    if distance2(self.points[slot], c) <= reach {
                        found.push(self.order[slot]);
                    }
                }
            } else {
                stack.push(cell.first);
                stack.push(at + 1);
            }
        }
        found.sort_unstable();
        found
    }
}

fn cell_gap(cell: &Cell, p: P) -> f64 {
    gap(&cell.low, &cell.high, p)
}

/// Build the cells, reordering `order` into leaf order: depth first, each
/// lower half right after its parent.
fn build_cells(points: &[P], order: &mut [u32]) -> Vec<Cell> {
    let mut cells: Vec<Cell> = Vec::new();
    if order.is_empty() {
        return cells;
    }
    let mut pending: Vec<(usize, usize, Option<usize>)> = vec![(0, order.len(), None)];
    while let Some((start, end, parent)) = pending.pop() {
        let index = cells.len();
        if let Some(parent) = parent {
            cells[parent].first = u32::try_from(index).unwrap_or(u32::MAX);
        }
        let range = &mut order[start..end];
        let (mut low, mut high) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for &i in range.iter() {
            let p = points[i as usize];
            for k in 0..3 {
                low[k] = low[k].min(p[k]);
                high[k] = high[k].max(p[k]);
            }
        }
        let count = end - start;
        let cell = Cell {
            low,
            high,
            first: u32::try_from(start).unwrap_or(u32::MAX),
            count: u32::try_from(count).unwrap_or(u32::MAX),
        };
        if count <= LEAF_MOST {
            cells.push(cell);
            continue;
        }
        let spread = [0, 1, 2].map(|k| high[k] - low[k]);
        let axis = if spread[0] >= spread[1] && spread[0] >= spread[2] {
            0
        } else if spread[1] >= spread[2] {
            1
        } else {
            2
        };
        let half = count / 2;
        range.select_nth_unstable_by(half, |&a, &b| {
            points[a as usize][axis]
                .total_cmp(&points[b as usize][axis])
                .then(a.cmp(&b))
        });
        cells.push(Cell {
            first: 0,
            count: 0,
            ..cell
        });
        pending.push((start + half, end, Some(index)));
        pending.push((start, start + half, None));
    }
    cells
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::bvh::tests::Seeded;

    fn brute(points: &[Point], p: Point, max: f64) -> Vec<Neighbour> {
        let q = [p.x, p.y, p.z];
        let mut all: Vec<(f64, u32)> = points
            .iter()
            .enumerate()
            .map(|(i, x)| (distance2([x.x, x.y, x.z], q), u32::try_from(i).unwrap()))
            .filter(|&(d2, _)| d2 <= max * max)
            .collect();
        all.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        all.into_iter()
            .map(|(d2, index)| Neighbour {
                index,
                distance: d2.sqrt(),
            })
            .collect()
    }

    /// A cloud with repeats and points on a coarse lattice, so ties in
    /// distance are common, plus uniform noise.
    fn cloud(count: usize, rng: &mut Seeded) -> Vec<Point> {
        let mut points = Vec::with_capacity(count);
        for i in 0..count {
            let p = match i % 4 {
                0 if i > 0 => points[rng.below(i)],
                1 => {
                    let at =
                        |rng: &mut Seeded| f64::from(u32::try_from(rng.below(5)).unwrap()) * 0.5;
                    Point::new(at(rng), at(rng), at(rng))
                }
                _ => rng.point(0.0, 2.0),
            };
            points.push(p);
        }
        points
    }

    /// Nearest, k nearest and ball queries agree with brute force, bit for
    /// bit and tie for tie: 1e5 of each.
    #[test]
    fn every_query_equals_brute_force_on_random_clouds() {
        let mut rng = Seeded(0x0009_0177);
        let mut found = [0_usize; 3];
        for round in 0..10 {
            let points = cloud(100 + 40 * round, &mut rng);
            let index = PointIndex::build(&points).unwrap();
            assert_eq!(index.len(), points.len());
            for q in 0..10_000 {
                let p = if q % 3 == 0 {
                    points[rng.below(points.len())]
                } else {
                    rng.point(-0.5, 2.5)
                };
                let max = if q % 2 == 0 {
                    f64::INFINITY
                } else {
                    rng.range(0.0, 0.5)
                };
                let all = brute(&points, p, max);
                let near = index.nearest(p, max);
                assert_eq!(near, all.first().copied(), "{p:?} {max}");
                found[0] += usize::from(near.is_some());
                let k = 1 + rng.below(12);
                let some = index.k_nearest(p, k, max);
                assert_eq!(some, all[..k.min(all.len())], "{p:?} {k} {max}");
                found[1] += usize::from(some.len() == k);
                let r = rng.range(0.0, 0.5);
                let mut ball: Vec<u32> = brute(&points, p, r).iter().map(|n| n.index).collect();
                ball.sort_unstable();
                assert_eq!(index.within_ball(p, r), ball, "{p:?} {r}");
                found[2] += usize::from(!ball.is_empty());
            }
        }
        assert!(found.iter().all(|&n| n > 25_000), "{found:?}");
    }

    /// A query on a repeated point goes to the lowest of its copies.
    #[test]
    fn ties_go_to_the_lowest_index() {
        let p = Point::new(1.0, 2.0, 3.0);
        let mut points = vec![Point::new(9.0, 9.0, 9.0); 20];
        points[7] = p;
        points[13] = p;
        points[3] = Point::new(1.0, 2.0, 4.0);
        points[17] = Point::new(1.0, 2.0, 2.0);
        let index = PointIndex::build(&points).unwrap();
        assert_eq!(
            index.nearest(p, 0.0),
            Some(Neighbour {
                index: 7,
                distance: 0.0
            })
        );
        let three: Vec<u32> = index.k_nearest(p, 3, 1.0).iter().map(|n| n.index).collect();
        assert_eq!(three, vec![7, 13, 3]);
        assert_eq!(index.k_nearest(p, 50, 1.0).len(), 4);
        assert_eq!(index.within_ball(p, 1.0), vec![3, 7, 13, 17]);
        ogeom_core::parallel::set_threads(4);
        let many = index.nearest_many(&vec![p; 600], 0.5);
        ogeom_core::parallel::set_threads(0);
        assert!(many.iter().all(|n| *n
            == Some(Neighbour {
                index: 7,
                distance: 0.0
            })));
    }

    #[test]
    fn an_empty_cloud_answers_nothing_and_a_bad_point_is_refused() {
        let index = PointIndex::build(&[]).unwrap();
        assert!(index.is_empty());
        let p = Point::new(0.0, 0.0, 0.0);
        assert!(index.nearest(p, f64::INFINITY).is_none());
        assert!(index.k_nearest(p, 3, f64::INFINITY).is_empty());
        assert!(index.within_ball(p, 1.0).is_empty());
        assert!(PointIndex::build(&[p, Point::new(0.0, f64::INFINITY, 0.0)]).is_err());
        let one = PointIndex::build(&[p]).unwrap();
        assert!(one.k_nearest(p, 0, 1.0).is_empty());
        assert!(one.nearest(Point::new(2.0, 0.0, 0.0), 1.0).is_none());
        // A negative reach is no reach, whatever its square.
        assert!(one.nearest(Point::new(1.0, 0.0, 0.0), -1.0).is_none());
        assert!(one.within_ball(Point::new(1.0, 0.0, 0.0), -1.0).is_empty());
    }
}
