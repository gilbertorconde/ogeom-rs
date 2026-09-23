//! Points binned on a uniform grid, so the entries near a point are found
//! among the cells around it rather than by walking every entry.

use std::collections::HashMap;

use ogeom_math::Point;

/// The cells a query may search before a walk over every entry is the
/// cheaper answer.
const MOST_CELLS: i64 = 512;

/// Entry indices binned by the cell their point falls in.
#[derive(Debug)]
pub(crate) struct Bins {
    /// The cell's width; a width that is not finite and positive bins
    /// nothing, and every query answers `None`.
    cell: f64,
    bins: HashMap<(i64, i64, i64), Vec<usize>>,
}

impl Bins {
    pub(crate) fn new(cell: f64) -> Self {
        Self {
            cell,
            bins: HashMap::new(),
        }
    }

    fn usable(&self) -> bool {
        self.cell.is_finite() && self.cell > 0.0
    }

    // A cast saturates, so a coordinate past `i64`'s cells shares the last
    // one; the division rounds monotonically, so a point never bins below
    // a lower one.
    #[allow(clippy::cast_possible_truncation, reason = "saturating")]
    fn key(&self, x: f64) -> i64 {
        (x / self.cell).floor() as i64
    }

    pub(crate) fn insert(&mut self, p: Point, index: usize) {
        if self.usable() {
            self.bins
                .entry((self.key(p.x), self.key(p.y), self.key(p.z)))
                .or_default()
                .push(index);
        }
    }

    /// Every entry binned within `radius` of `p` on each axis, in ascending
    /// order — a superset of the entries whose point lies that close — or
    /// `None` when the grid cannot answer: the cell is unusable, the radius
    /// is not finite, or it spans more cells than a walk over every entry
    /// would cost.
    pub(crate) fn near(&self, p: Point, radius: f64) -> Option<Vec<usize>> {
        if !self.usable() || !radius.is_finite() {
            return None;
        }
        // The box's corners are rounded outward by more than the
        // subtraction and addition can lose, so a point on its face is
        // never binned outside it.
        let span = |c: f64| {
            let slack = (c.abs() + radius) * 1e-12;
            (self.key(c - radius - slack), self.key(c + radius + slack))
        };
        let ranges = [span(p.x), span(p.y), span(p.z)];
        let mut cells: i64 = 1;
        for (lo, hi) in ranges {
            cells = cells.saturating_mul(hi.saturating_sub(lo).saturating_add(1));
        }
        if cells > MOST_CELLS {
            return None;
        }
        let mut out = Vec::new();
        for x in ranges[0].0..=ranges[0].1 {
            for y in ranges[1].0..=ranges[1].1 {
                for z in ranges[2].0..=ranges[2].1 {
                    if let Some(list) = self.bins.get(&(x, y, z)) {
                        out.extend_from_slice(list);
                    }
                }
            }
        }
        out.sort_unstable();
        Some(out)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    /// Every point within the radius is found, whatever cell it binned in,
    /// and a query too wide for the grid says so.
    #[test]
    fn near_finds_every_point_within_reach() {
        let mut bins = Bins::new(0.5);
        let points: Vec<Point> = (0..200_u32)
            .map(|i| {
                let t = f64::from(i);
                Point::new((t * 0.37).sin() * 3.0, (t * 0.11).cos() * 3.0, t * 0.01)
            })
            .collect();
        for (i, p) in points.iter().enumerate() {
            bins.insert(*p, i);
        }
        for q in &points {
            let found = bins.near(*q, 0.7).unwrap();
            assert!(found.windows(2).all(|w| w[0] < w[1]));
            for (i, p) in points.iter().enumerate() {
                if p.distance(*q) <= 0.7 {
                    assert!(found.binary_search(&i).is_ok());
                }
            }
        }
        assert!(bins.near(Point::ORIGIN, 100.0).is_none());
        assert!(Bins::new(0.0).near(Point::ORIGIN, 1.0).is_none());
    }
}
