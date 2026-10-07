//! Where a projected curve crosses the drawing's contours.
//!
//! Visibility along a curve changes only where its projection passes
//! behind a contour (a silhouette or a free boundary) or where it pierces a
//! face, so a curve split at its contour crossings has one visibility per
//! piece, up to piercings (quantitative invisibility). The contours are
//! projected segments binned in a grid over the drawing, so a query walks
//! only the cells its own segment passes through.

use ogeom_core::Tolerances;
use ogeom_math::Point2;

/// Projected contour segments, binned by the cells they cover.
pub(crate) struct Contours {
    segments: Vec<(Point2, Point2)>,
    cells: Vec<Vec<u32>>,
    low: Point2,
    size: f64,
    side: usize,
}

impl Contours {
    /// Bin `segments` in a grid of about one segment per cell.
    pub(crate) fn new(segments: Vec<(Point2, Point2)>) -> Self {
        let (mut low, mut high) = (
            Point2::new(f64::INFINITY, f64::INFINITY),
            Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        );
        for (a, b) in &segments {
            for p in [a, b] {
                low = Point2::new(low.x.min(p.x), low.y.min(p.y));
                high = Point2::new(high.x.max(p.x), high.y.max(p.y));
            }
        }
        #[allow(
            clippy::cast_precision_loss,
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss
        )]
        let side = ((segments.len() as f64).sqrt().ceil() as usize).clamp(1, 512);
        let span = (high.x - low.x).max(high.y - low.y);
        #[allow(clippy::cast_precision_loss)]
        let size = if span.is_finite() && span > 0.0 {
            span / side as f64
        } else {
            1.0
        };
        if !low.x.is_finite() || !low.y.is_finite() {
            low = Point2::new(0.0, 0.0);
        }
        let mut out = Self {
            segments: Vec::new(),
            cells: vec![Vec::new(); side * side],
            low,
            size,
            side,
        };
        for (index, (a, b)) in segments.iter().enumerate() {
            if !(a.x.is_finite() && a.y.is_finite() && b.x.is_finite() && b.y.is_finite()) {
                continue;
            }
            #[allow(clippy::cast_possible_truncation)]
            let index = index as u32;
            // A cell's worth of margin each way, so a crossing found a hair
            // outside either segment's box still meets it in a shared cell.
            let (c0, c1) = (
                out.cell(a.x.min(b.x) - size, out.low.x),
                out.cell(a.x.max(b.x) + size, out.low.x),
            );
            let (r0, r1) = (
                out.cell(a.y.min(b.y) - size, out.low.y),
                out.cell(a.y.max(b.y) + size, out.low.y),
            );
            for r in r0..=r1 {
                for c in c0..=c1 {
                    out.cells[r * side + c].push(index);
                }
            }
        }
        out.segments = segments;
        out
    }

    fn cell(&self, x: f64, lo: f64) -> usize {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let k = ((x - lo) / self.size).floor().max(0.0) as usize;
        k.min(self.side - 1)
    }

    /// The contour segments whose cells the segment `a`-`b` passes
    /// through, each once.
    fn near(&self, a: Point2, b: Point2) -> Vec<u32> {
        let mut out: Vec<u32> = Vec::new();
        let (c0, c1) = (
            self.cell(a.x.min(b.x), self.low.x),
            self.cell(a.x.max(b.x), self.low.x),
        );
        let dx = b.x - a.x;
        for c in c0..=c1 {
            // The stretch of the segment over this column, by its y range.
            #[allow(clippy::cast_precision_loss)]
            let (x0, x1) = (
                (self.low.x + c as f64 * self.size).max(a.x.min(b.x)),
                (self.low.x + (c + 1) as f64 * self.size).min(a.x.max(b.x)),
            );
            let (y0, y1) = if dx.abs() <= f64::EPSILON * (a.x.abs() + b.x.abs() + 1.0) {
                (a.y.min(b.y), a.y.max(b.y))
            } else {
                let at = |x: f64| ((x - a.x) / dx).clamp(0.0, 1.0).mul_add(b.y - a.y, a.y);
                (at(x0).min(at(x1)), at(x0).max(at(x1)))
            };
            let (r0, r1) = (
                self.cell(y0, self.low.y).saturating_sub(1),
                (self.cell(y1, self.low.y) + 1).min(self.side - 1),
            );
            for r in r0..=r1 {
                out.extend_from_slice(&self.cells[r * self.side + c]);
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    /// Where the segment `a`-`b` crosses a contour, as fractions of the
    /// way from `a` to `b` strictly inside it, ascending, no two closer
    /// than the confusion tolerance along the segment. Parallel contours
    /// cross nothing.
    pub(crate) fn crossings(&self, a: Point2, b: Point2, tol: Tolerances) -> Vec<f64> {
        let (dx, dy) = (b.x - a.x, b.y - a.y);
        let length = dx.hypot(dy);
        if length <= tol.confusion() {
            return Vec::new();
        }
        let near = tol.confusion() / length;
        let mut out: Vec<f64> = Vec::new();
        for index in self.near(a, b) {
            let (c, d) = self.segments[index as usize];
            let (ex, ey) = (d.x - c.x, d.y - c.y);
            let across = ex.hypot(ey);
            if across <= tol.confusion() {
                continue;
            }
            let denominator = dx * ey - dy * ex;
            if denominator.abs() <= 1e-12 * length * across {
                continue;
            }
            let (wx, wy) = (c.x - a.x, c.y - a.y);
            let t = (wx * ey - wy * ex) / denominator;
            let s = (wx * dy - wy * dx) / denominator;
            let slack = tol.confusion() / across;
            if t > near && t < 1.0 - near && s >= -slack && s <= 1.0 + slack {
                out.push(t);
            }
        }
        out.sort_by(f64::total_cmp);
        out.dedup_by(|later, earlier| *later - *earlier <= near);
        out
    }

    /// Positions along a polyline (`k + f` is the fraction `f` of the way
    /// along segment `k`) at its own points and where it crosses a contour,
    /// ascending.
    pub(crate) fn split(&self, projected: &[Point2], tol: Tolerances) -> Vec<f64> {
        let mut out = Vec::with_capacity(projected.len());
        for (k, pair) in projected.windows(2).enumerate() {
            #[allow(clippy::cast_precision_loss)]
            let k = k as f64;
            out.push(k);
            out.extend(
                self.crossings(pair[0], pair[1], tol)
                    .into_iter()
                    .map(|f| k + f),
            );
        }
        if let Some(last) = projected.len().checked_sub(1) {
            #[allow(clippy::cast_precision_loss)]
            out.push(last as f64);
        }
        out
    }
}

/// The segment a position along a polyline of `points` points falls in,
/// and the fraction along it.
pub(crate) fn locate(position: f64, points: usize) -> (usize, f64) {
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss
    )]
    let k = (position.floor().max(0.0) as usize).min(points.saturating_sub(2));
    #[allow(clippy::cast_precision_loss)]
    let f = position - k as f64;
    (k, f)
}

#[cfg(test)]
mod tests {
    use super::*;

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn a_segment_is_split_where_it_crosses_and_not_where_it_touches() {
        let contours = Contours::new(vec![
            (Point2::new(3.0, -1.0), Point2::new(3.0, 1.0)),
            // Parallel: no crossing.
            (Point2::new(0.0, 0.5), Point2::new(10.0, 0.5)),
            // Ends on the segment's own end: not inside it.
            (Point2::new(0.0, 0.0), Point2::new(0.0, 5.0)),
            (Point2::new(7.5, 0.0), Point2::new(9.0, 3.0)),
            // Far off to the side.
            (Point2::new(5.0, 4.0), Point2::new(6.0, 5.0)),
        ]);
        let found = contours.crossings(Point2::new(0.0, 0.0), Point2::new(10.0, 0.0), T);
        assert_eq!(found.len(), 2, "{found:?}");
        assert!((found[0] - 0.3).abs() < 1e-12);
        assert!((found[1] - 0.75).abs() < 1e-12);
    }
}
