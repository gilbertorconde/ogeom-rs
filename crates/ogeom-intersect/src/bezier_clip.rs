//! A spline surface cut into its Bézier patches, and the places a straight
//! segment may meet it, found by the patches' control hulls.
//!
//! A rational Bézier patch with positive weights lies in the convex hull of
//! its control points, and so in their box: a segment that misses the box
//! grown by a margin comes no nearer the patch than the margin. Halving a
//! patch by de Casteljau gives two patches whose hulls hug the surface
//! closer, so the boxes a segment meets narrow to where it comes near the
//! surface, and nowhere else is asked.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::BSplineSurface;
use ogeom_math::Point;
use ogeom_math::bspline::{Blend, Weighted, to_bezier_segments};

/// One Bézier patch: its parameter rectangle and its control net, `u`
/// major, `(p + 1) * (q + 1)` points.
#[derive(Debug, Clone)]
struct Patch {
    u: (f64, f64),
    v: (f64, f64),
    net: Vec<Weighted<Point>>,
}

/// A spline surface's Bézier patches, ready for segments.
#[derive(Debug, Clone)]
pub(crate) struct Patches {
    p: usize,
    q: usize,
    patches: Vec<Patch>,
    /// The size below which a patch is a seed rather than halved again.
    leaf: f64,
}

impl Patches {
    /// The patches of `spline`, or `None` where it cannot be cut into them.
    pub(crate) fn of(spline: &BSplineSurface, tol: Tolerances) -> OgeomResult<Option<Self>> {
        let (uk, vk) = (spline.u_knots(), spline.v_knots());
        let grid = spline.grid();
        let (p, q) = (uk.degree(), vk.degree());
        let (nu, nv) = (grid.u_count(), grid.v_count());
        // Each row along `u` cut into its segments; every row shares the
        // knots, so they cut alike.
        let mut rows = Vec::with_capacity(nv);
        for j in 0..nv {
            let mut row = Vec::with_capacity(nu);
            for i in 0..nu {
                let Some(point) = grid.get(i, j) else {
                    return Ok(None);
                };
                row.push(point);
            }
            rows.push(to_bezier_segments(uk, &row, tol)?);
        }
        let Some(first) = rows.first() else {
            return Ok(None);
        };
        let mut patches = Vec::new();
        let (mut low, mut high) = (
            Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
            Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
        );
        for (s, &(u, _)) in first.iter().enumerate() {
            // Column `i` of this `u` segment, along `v`, cut into segments.
            let mut columns = Vec::with_capacity(p + 1);
            for i in 0..=p {
                let column: Vec<Weighted<Point>> = rows.iter().map(|r| r[s].1[i]).collect();
                columns.push(to_bezier_segments(vk, &column, tol)?);
            }
            for t in 0..columns[0].len() {
                let v = columns[0][t].0;
                let mut net = Vec::with_capacity((p + 1) * (q + 1));
                for column in &columns {
                    net.extend_from_slice(&column[t].1);
                }
                for w in &net {
                    if w.weight <= 0.0 {
                        return Ok(None);
                    }
                    let r = w.point();
                    low = Point::new(low.x.min(r.x), low.y.min(r.y), low.z.min(r.z));
                    high = Point::new(high.x.max(r.x), high.y.max(r.y), high.z.max(r.z));
                }
                patches.push(Patch { u, v, net });
            }
        }
        let diagonal = (high - low).magnitude();
        if !diagonal.is_finite() {
            return Ok(None);
        }
        Ok(Some(Self {
            p,
            q,
            patches,
            leaf: diagonal * LEAF,
        }))
    }

    /// Where the segment `a..b` may meet the surface within `margin`: for
    /// each smallest patch whose control box the segment meets, the
    /// fraction along the segment nearest its middle and its middle's
    /// parameters.
    pub(crate) fn seeds(&self, a: Point, b: Point, margin: f64) -> Vec<(f64, (f64, f64))> {
        let mut out = Vec::new();
        let mut stack: Vec<(Patch, u32)> = self.patches.iter().map(|p| (p.clone(), 0)).collect();
        while let Some((patch, depth)) = stack.pop() {
            let (low, high) = net_box(&patch.net);
            if !segment_meets_box(a, b, low, high, margin) {
                continue;
            }
            let size = (high - low).magnitude();
            if size <= self.leaf || depth >= DEEPEST {
                let middle = low.midpoint(high);
                let ab = b - a;
                let length = ab.dot(ab);
                let f = if length > 0.0 {
                    ((middle - a).dot(ab) / length).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                out.push((
                    f,
                    (
                        f64::midpoint(patch.u.0, patch.u.1),
                        f64::midpoint(patch.v.0, patch.v.1),
                    ),
                ));
                continue;
            }
            // Halved across the direction the net runs farther in.
            let (along_u, along_v) = self.extents(&patch.net);
            let (one, two) = if along_u >= along_v {
                self.halve_u(&patch)
            } else {
                self.halve_v(&patch)
            };
            stack.push((one, depth + 1));
            stack.push((two, depth + 1));
        }
        out
    }

    fn at(&self, net: &[Weighted<Point>], i: usize, j: usize) -> Weighted<Point> {
        net[i * (self.q + 1) + j]
    }

    /// How far the net's first and last rows stand apart each way.
    fn extents(&self, net: &[Weighted<Point>]) -> (f64, f64) {
        let (p, q) = (self.p, self.q);
        let mut along_u: f64 = 0.0;
        let mut along_v: f64 = 0.0;
        for j in 0..=q {
            along_u = along_u.max(
                self.at(net, 0, j)
                    .point()
                    .distance(self.at(net, p, j).point()),
            );
        }
        for i in 0..=p {
            along_v = along_v.max(
                self.at(net, i, 0)
                    .point()
                    .distance(self.at(net, i, q).point()),
            );
        }
        (along_u, along_v)
    }

    fn halve_u(&self, patch: &Patch) -> (Patch, Patch) {
        let (p, q) = (self.p, self.q);
        let mut left = vec![Weighted::zero(); patch.net.len()];
        let mut right = left.clone();
        for j in 0..=q {
            let column: Vec<Weighted<Point>> = (0..=p).map(|i| self.at(&patch.net, i, j)).collect();
            let (l, r) = de_casteljau(&column);
            for i in 0..=p {
                left[i * (q + 1) + j] = l[i];
                right[i * (q + 1) + j] = r[i];
            }
        }
        let mid = f64::midpoint(patch.u.0, patch.u.1);
        (
            Patch {
                u: (patch.u.0, mid),
                v: patch.v,
                net: left,
            },
            Patch {
                u: (mid, patch.u.1),
                v: patch.v,
                net: right,
            },
        )
    }

    fn halve_v(&self, patch: &Patch) -> (Patch, Patch) {
        let (p, q) = (self.p, self.q);
        let mut lower = vec![Weighted::zero(); patch.net.len()];
        let mut upper = lower.clone();
        for i in 0..=p {
            let row: Vec<Weighted<Point>> = (0..=q).map(|j| self.at(&patch.net, i, j)).collect();
            let (l, r) = de_casteljau(&row);
            for j in 0..=q {
                lower[i * (q + 1) + j] = l[j];
                upper[i * (q + 1) + j] = r[j];
            }
        }
        let mid = f64::midpoint(patch.v.0, patch.v.1);
        (
            Patch {
                u: patch.u,
                v: (patch.v.0, mid),
                net: lower,
            },
            Patch {
                u: patch.u,
                v: (mid, patch.v.1),
                net: upper,
            },
        )
    }
}

/// A patch this small, against the surface's whole size, is a seed.
const LEAF: f64 = 1.0 / 512.0;

/// The most halvings a patch takes before it is a seed regardless.
const DEEPEST: u32 = 24;

/// A Bézier curve's points halved at its middle: the two halves' points.
fn de_casteljau(points: &[Weighted<Point>]) -> (Vec<Weighted<Point>>, Vec<Weighted<Point>>) {
    let n = points.len();
    let mut work = points.to_vec();
    let mut left = Vec::with_capacity(n);
    let mut right = vec![Weighted::zero(); n];
    left.push(work[0]);
    right[n - 1] = work[n - 1];
    for level in 1..n {
        for k in 0..n - level {
            work[k] = work[k].lerp(work[k + 1], 0.5);
        }
        left.push(work[0]);
        right[n - 1 - level] = work[n - 1 - level];
    }
    (left, right)
}

/// The box of a net's points.
fn net_box(net: &[Weighted<Point>]) -> (Point, Point) {
    let mut low = Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut high = Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for w in net {
        let r = w.point();
        low = Point::new(low.x.min(r.x), low.y.min(r.y), low.z.min(r.z));
        high = Point::new(high.x.max(r.x), high.y.max(r.y), high.z.max(r.z));
    }
    (low, high)
}

/// Whether the segment `a..b` meets the box `low..high` grown by `margin`.
pub(crate) fn segment_meets_box(a: Point, b: Point, low: Point, high: Point, margin: f64) -> bool {
    let (mut t0, mut t1) = (0.0_f64, 1.0_f64);
    for (s, d, lo, hi) in [
        (a.x, b.x - a.x, low.x - margin, high.x + margin),
        (a.y, b.y - a.y, low.y - margin, high.y + margin),
        (a.z, b.z - a.z, low.z - margin, high.z + margin),
    ] {
        if d == 0.0 {
            if s < lo || s > hi {
                return false;
            }
            continue;
        }
        let (e0, e1) = ((lo - s) / d, (hi - s) / d);
        t0 = t0.max(e0.min(e1));
        t1 = t1.min(e0.max(e1));
        if t0 > t1 {
            return false;
        }
    }
    true
}
