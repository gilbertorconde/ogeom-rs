//! Surfaces sampled into triangles, and where triangles and segments cross.

use super::*;

/// One sampled triangle of a surface, with the parameters it came from.
#[derive(Debug, Clone)]
pub(crate) struct Cell {
    pub(crate) corners: [Point; 3],
    pub(crate) at: (f64, f64),
    pub(crate) low: Point,
    pub(crate) high: Point,
    /// How far the surface bows from the flat cell: its middle's distance
    /// from the middle of the diagonal the two triangles share.
    pub(crate) sag: f64,
    /// The parameters of the three corners, in `corners` order.
    pub(crate) params: [(f64, f64); 3],
}

/// Sample a surface into triangles.
pub(crate) fn sample(surface: &SurfaceGeometry, grid: usize, tol: Tolerances) -> Vec<Cell> {
    sample_by(surface, (grid, grid), tol)
}

/// Sample a surface into triangles, `counts` cells along `u` and along `v`.
pub(crate) fn sample_by(
    surface: &SurfaceGeometry,
    counts: (usize, usize),
    tol: Tolerances,
) -> Vec<Cell> {
    let ((ua, ub), (va, vb)) = surface.domain();
    // An unbounded domain would put the samples a billion units apart and find
    // nothing. Clamped to something a real model lives inside.
    let limit = 1.0e6;
    let (ua, ub) = (ua.max(-limit), ub.min(limit));
    let (va, vb) = (va.max(-limit), vb.min(limit));

    let mut out = Vec::new();
    #[allow(clippy::cast_precision_loss)]
    let (nu, nv) = (counts.0 as f64, counts.1 as f64);
    let at = |s: f64, t: f64| {
        let (u, v) = (ua + (ub - ua) * s, va + (vb - va) * t);
        surface.point_at(u, v, tol).ok().map(|p| ((u, v), p))
    };
    // Each grid corner is evaluated once and shared by the cells around it,
    // a row of corners at a time.
    #[allow(clippy::cast_precision_loss)]
    let row = |i: usize| -> Vec<_> {
        (0..=counts.1)
            .map(|j| at(i as f64 / nu, j as f64 / nv))
            .collect()
    };
    let mut below = row(0);
    for i in 0..counts.0 {
        let above = row(i + 1);
        for j in 0..counts.1 {
            #[allow(clippy::cast_precision_loss)]
            let (s0, s1) = (i as f64 / nu, (i + 1) as f64 / nu);
            #[allow(clippy::cast_precision_loss)]
            let (t0, t1) = (j as f64 / nv, (j + 1) as f64 / nv);
            let (Some((p00, a00)), Some((p10, a10)), Some((p01, a01)), Some((p11, a11))) =
                (below[j], above[j], below[j + 1], above[j + 1])
            else {
                continue;
            };
            let sag = at(f64::midpoint(s0, s1), f64::midpoint(t0, t1))
                .map_or(0.0, |(_, middle)| middle.distance(a00.midpoint(a11)));
            for (corners, params) in [
                ([a00, a10, a11], [p00, p10, p11]),
                ([a00, a11, a01], [p00, p11, p01]),
            ] {
                let low = Point::new(
                    corners.iter().map(|p| p.x).fold(f64::MAX, f64::min),
                    corners.iter().map(|p| p.y).fold(f64::MAX, f64::min),
                    corners.iter().map(|p| p.z).fold(f64::MAX, f64::min),
                );
                let high = Point::new(
                    corners.iter().map(|p| p.x).fold(f64::MIN, f64::max),
                    corners.iter().map(|p| p.y).fold(f64::MIN, f64::max),
                    corners.iter().map(|p| p.z).fold(f64::MIN, f64::max),
                );
                out.push(Cell {
                    corners,
                    at: p00,
                    low,
                    high,
                    sag,
                    params,
                });
            }
        }
        below = above;
    }
    out
}

/// Whether two cells' boxes come within a margin of each other.
/// Sampled cells binned by their boxes on a coarse grid, for finding the
/// cells whose boxes could meet a given one without testing them all.
pub(super) struct CellBins {
    low: Point,
    size: f64,
    counts: [usize; 3],
    bins: Vec<Vec<usize>>,
    margin: f64,
}

impl CellBins {
    /// At most this many bins a side: a grid on a clamped plane's million
    /// units would otherwise be all empty bins.
    const MOST: usize = 48;

    pub(super) fn over(cells: &[Cell], margin: f64) -> Self {
        let mut low = Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut high = Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for c in cells {
            low = Point::new(low.x.min(c.low.x), low.y.min(c.low.y), low.z.min(c.low.z));
            high = Point::new(
                high.x.max(c.high.x),
                high.y.max(c.high.y),
                high.z.max(c.high.z),
            );
        }
        let extent = (high - low).magnitude();
        #[allow(clippy::cast_precision_loss)]
        let size = if extent.is_finite() && extent > 0.0 {
            (extent / Self::MOST as f64).max(margin)
        } else {
            1.0
        };
        let count = |lo: f64, hi: f64| -> usize {
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let n = ((hi - lo) / size).floor() as usize + 1;
            n.clamp(1, Self::MOST + 1)
        };
        let counts = if extent.is_finite() {
            [
                count(low.x, high.x),
                count(low.y, high.y),
                count(low.z, high.z),
            ]
        } else {
            [1, 1, 1]
        };
        let mut bins = vec![Vec::new(); counts[0] * counts[1] * counts[2]];
        let mut this = Self {
            low,
            size,
            counts,
            bins: Vec::new(),
            margin,
        };
        for (i, c) in cells.iter().enumerate() {
            this.each_bin(c.low, c.high, |k| bins[k].push(i));
        }
        this.bins = bins;
        this
    }

    /// Every bin a box grown by the margin touches.
    fn each_bin(&self, low: Point, high: Point, mut visit: impl FnMut(usize)) {
        let index = |x: f64, lo: f64, n: usize| -> usize {
            if !x.is_finite() {
                return if x > 0.0 { n - 1 } else { 0 };
            }
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let k = ((x - lo) / self.size).floor().max(0.0) as usize;
            k.min(n - 1)
        };
        let m = self.margin;
        let [nx, ny, nz] = self.counts;
        let (x0, x1) = (
            index(low.x - m, self.low.x, nx),
            index(high.x + m, self.low.x, nx),
        );
        let (y0, y1) = (
            index(low.y - m, self.low.y, ny),
            index(high.y + m, self.low.y, ny),
        );
        let (z0, z1) = (
            index(low.z - m, self.low.z, nz),
            index(high.z + m, self.low.z, nz),
        );
        for x in x0..=x1 {
            for y in y0..=y1 {
                for z in z0..=z1 {
                    visit((x * ny + y) * nz + z);
                }
            }
        }
    }

    /// The cells whose boxes could meet `cell`'s, lowest index first.
    pub(super) fn candidates(&self, cell: &Cell, out: &mut Vec<usize>) {
        out.clear();
        self.each_bin(cell.low, cell.high, |k| {
            out.extend_from_slice(&self.bins[k])
        });
        out.sort_unstable();
        out.dedup();
    }
}

pub(super) fn overlap(a: &Cell, b: &Cell, margin: f64) -> bool {
    a.low.x <= b.high.x + margin
        && b.low.x <= a.high.x + margin
        && a.low.y <= b.high.y + margin
        && b.low.y <= a.high.y + margin
        && a.low.z <= b.high.z + margin
        && b.low.z <= a.high.z + margin
}

/// A point where two triangles cross, if they do.
///
/// Each triangle's edges are tested against the other's plane and then against
/// the triangle itself. Only an approximate answer is needed: it is a seed, and
/// the Newton correction that follows is what makes it a point on the curve.
pub(super) fn triangles_cross(a: &Cell, b: &Cell) -> Option<Point> {
    for (edges, target) in [(a, b), (b, a)] {
        for k in 0..3 {
            let (from, to) = (edges.corners[k], edges.corners[(k + 1) % 3]);
            if let Some(hit) = segment_meets_triangle(from, to, target.corners) {
                return Some(hit);
            }
        }
    }
    None
}

/// Where a segment crosses a triangle.
pub(crate) fn segment_meets_triangle(from: Point, to: Point, t: [Point; 3]) -> Option<Point> {
    let direction = to - from;
    let (e1, e2) = (t[1] - t[0], t[2] - t[0]);
    let h = direction.cross(e2);
    let determinant = e1.dot(h);
    if determinant.abs() <= f64::MIN_POSITIVE {
        return None;
    }
    let inverse = 1.0 / determinant;
    let s = from - t[0];
    let u = inverse * s.dot(h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = s.cross(e1);
    let v = inverse * direction.dot(q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let along = inverse * e2.dot(q);
    if !(0.0..=1.0).contains(&along) {
        return None;
    }
    Some(from + direction * along)
}
