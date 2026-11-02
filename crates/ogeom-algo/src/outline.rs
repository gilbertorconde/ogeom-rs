//! An outline of lines, arcs and circles fitted around a silhouette,
//! clear of it by a play and within a tolerance band.
//!
//! The band is every point whose signed distance from the silhouette,
//! positive outside it and negative inside, lies between `play` and
//! `play + 2 tolerance`. The fit has three stages, and only the last one
//! decides what is returned.
//!
//! 1. **The middle of the band.** The silhouette is first traced afresh as
//!    a region, so its loops bound it exactly. Its offset by `play +
//!    tolerance` is then the union of the silhouette with a rectangle along
//!    each edge, on the outside, and a sector at each convex corner, or for
//!    a negative offset the silhouette less such rectangles inside and
//!    sectors at the reflex corners. The sectors' arcs are chords sagging a
//!    quarter of the tolerance, the offset is taken an eighth further out
//!    to centre the sag, and the region engine of
//!    [`ogeom_mesh::planar`] traces the result: loops whose every point is
//!    within an eighth of the tolerance, and the engine's snapping, of the
//!    middle of the band. Holes the band closes over and parts it pinches
//!    off are gone from them, and as the boundary of a region they never
//!    run out and back.
//! 2. **The fit.** Each loop is split into runs a line or a circle follows
//!    within three quarters of the tolerance, less the snapping: the whole
//!    loop as one circle where that holds. A line or an arc spanning a run
//!    lies that close to the loop, so within the band. Consecutive pieces
//!    are joined where their curves meet: a corner turning less than
//!    `round_below` takes the largest tangent arc found to fit the band, a
//!    sharper one stays a corner where that fits and takes an arc where it
//!    does not, and a joint where neither fits is bridged by the short line
//!    between the pieces' ends, which lies near the loop too. No piece is
//!    moved past its middle by a joint, so the joints at its two ends
//!    never meet.
//! 3. **The check.** Every piece is sampled at a quarter of the tolerance
//!    and every sample's signed distance from the silhouette must lie an
//!    eighth of the tolerance inside the band. Distance is 1-Lipschitz, so
//!    the whole piece then lies in the band. A loop with any piece failing
//!    is returned as the polygon of its middle-of-the-band loop, which is
//!    checked the same way.

use ogeom_core::{Exact, OgeomResult, Predicates, Sign, ogeom_bail};
use ogeom_math::{Point2, Vector2};
use ogeom_mesh::{Edge2, Loop2, region};

use crate::fit::segments::{Shape, circle_through, deviation, runs};

/// The corner arcs of the band's middle sag from their circles by at most
/// this fraction of the tolerance.
const SAG: f64 = 32.0;

/// How an outline is fitted around a silhouette.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Outline {
    /// The least clearance from the silhouette. Negative for an outline
    /// inside it, a press fit.
    pub play: f64,
    /// Half the width of the band: every point of the outline lies between
    /// `play` and `play + 2 tolerance` outside the silhouette.
    pub tolerance: f64,
    /// Corners turning by less than this, in radians, are rounded by the
    /// largest tangent arc found to fit the band.
    pub round_below: f64,
}

/// One piece of a fitted outline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Piece2 {
    /// A straight piece.
    Line {
        /// Where it starts.
        from: Point2,
        /// Where it ends.
        to: Point2,
    },
    /// A circular arc, counter-clockwise for a positive sweep.
    Arc {
        /// Its centre.
        centre: Point2,
        /// Its radius.
        radius: f64,
        /// The angle it starts at, from the `x` axis.
        start: f64,
        /// The signed angle it turns through, less than a full turn.
        sweep: f64,
    },
    /// A full circle, a whole outline on its own.
    Circle {
        /// Its centre.
        centre: Point2,
        /// Its radius.
        radius: f64,
        /// Whether it runs counter-clockwise, as an outer boundary does.
        counter_clockwise: bool,
    },
}

impl Piece2 {
    /// The point a fraction `t` of the way along, from `0` at the start to
    /// `1` at the end.
    #[must_use]
    pub fn point_at(&self, t: f64) -> Point2 {
        match *self {
            Self::Line { from, to } => from.lerp(to, t),
            Self::Arc {
                centre,
                radius,
                start,
                sweep,
            } => on_circle(centre, radius, sweep.mul_add(t, start)),
            Self::Circle {
                centre,
                radius,
                counter_clockwise,
            } => {
                let turn = if counter_clockwise { 1.0 } else { -1.0 };
                on_circle(centre, radius, turn * core::f64::consts::TAU * t)
            }
        }
    }

    /// Where it starts.
    #[must_use]
    pub fn start(&self) -> Point2 {
        self.point_at(0.0)
    }

    /// Where it ends: where it starts, for a circle.
    #[must_use]
    pub fn end(&self) -> Point2 {
        self.point_at(1.0)
    }

    /// Its length.
    #[must_use]
    pub fn length(&self) -> f64 {
        match *self {
            Self::Line { from, to } => from.distance(to),
            Self::Arc { radius, sweep, .. } => radius * sweep.abs(),
            Self::Circle { radius, .. } => radius * core::f64::consts::TAU,
        }
    }
}

/// A closed outline: pieces end to end, each starting where the one
/// before ends and the last ending where the first starts. An outer
/// boundary runs counter-clockwise and a hole clockwise.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Wire2d {
    /// The pieces in order.
    pub pieces: Vec<Piece2>,
}

/// A whole number held in a double, as a count.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "callers pass whole numbers from zero to below 2^32"
)]
const fn whole(v: f64) -> usize {
    v as usize
}

/// `k / n` for counts, exactly: counts are below `2^32`.
#[allow(clippy::cast_precision_loss, reason = "counts below 2^32 are exact")]
fn fraction(k: usize, n: usize) -> f64 {
    k as f64 / n as f64
}

fn on_circle(centre: Point2, radius: f64, angle: f64) -> Point2 {
    Point2::new(
        radius.mul_add(angle.cos(), centre.x),
        radius.mul_add(angle.sin(), centre.y),
    )
}

/// Fit an outline of lines, tangent arcs and circles around a silhouette,
/// every point of it between `play` and `play + 2 tolerance` outside the
/// region the loops bound. See the module documentation.
///
/// The loops bound the points they wind around at least once, as
/// [`shadow`](ogeom_mesh::shadow()) gives them: outer boundaries
/// counter-clockwise, holes clockwise. The outline has one wire for each
/// loop of the band's middle: a hole narrower than the band is closed
/// over, and with a negative play a part thinner than the band is gone.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// the tolerance is not positive, the play is not finite, or `round_below`
/// is not an angle from zero to a half turn, and
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if a point is not
/// finite.
pub fn fit_outline(loops: &[Loop2], outline: &Outline) -> OgeomResult<Vec<Wire2d>> {
    let Outline {
        play,
        tolerance,
        round_below,
    } = *outline;
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        ogeom_bail!(
            Construction,
            "an outline tolerance of {tolerance} is no band"
        );
    }
    if !play.is_finite() {
        ogeom_bail!(Construction, "an outline play of {play} is not a distance");
    }
    if !(0.0..=core::f64::consts::PI).contains(&round_below) {
        ogeom_bail!(
            Construction,
            "rounding corners turning less than {round_below} rad is not a turn"
        );
    }
    let silhouette = region(&chain(loops), |w| w[0] > 0)?;
    if silhouette.is_empty() {
        return Ok(Vec::new());
    }
    let snap = snapping(&silhouette);
    let middle = play + tolerance;
    let targets = if middle == 0.0 {
        silhouette.clone()
    } else {
        offset(&silhouette, middle, tolerance / SAG)?
    };
    let spacing = tolerance / 4.0;
    // The margin covers half the spacing, and the silhouette's own
    // snapping off the loops it was traced from.
    let margin = spacing / 2.0 + snap;
    let band = Band::new(
        &silhouette,
        play + margin,
        2.0f64.mul_add(tolerance, play) - margin,
    );
    let fit = Fit {
        band: &band,
        within: tolerance
            - tolerance / (2.0 * SAG)
            - spacing / 2.0
            - 2.0 * (snap + snapping(&targets)),
        round_below,
        spacing,
    };
    let mut wires = Vec::with_capacity(targets.len());
    for target in &targets {
        if target.points.len() < 3 {
            continue;
        }
        wires.push(fit.wire(target)?);
    }
    Ok(wires)
}

/// The loops as a chain counted in the first count.
fn chain(loops: &[Loop2]) -> Vec<Edge2> {
    let mut edges = Vec::new();
    for l in loops {
        let n = l.points.len();
        for i in 0..n {
            edges.push(Edge2 {
                from: l.points[i],
                to: l.points[(i + 1) % n],
                weight: [1, 0],
            });
        }
    }
    edges
}

/// A bound on how far the region engine moves a boundary traced from
/// these loops: a cell's diagonal, the cell at most `2^-37` of the extent.
fn snapping(loops: &[Loop2]) -> f64 {
    let mut lo = [f64::INFINITY; 2];
    let mut hi = [f64::NEG_INFINITY; 2];
    for p in loops.iter().flat_map(|l| &l.points) {
        lo = [lo[0].min(p.x), lo[1].min(p.y)];
        hi = [hi[0].max(p.x), hi[1].max(p.y)];
    }
    let extent = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(0.0);
    extent * 2f64.powi(-36)
}

/// The silhouette offset by `distance`, outward for a positive one, its
/// corner arcs as chords sagging at most `sag`, the offset taken half the
/// sag further so that the boundary lies within half the sag of the exact
/// offset either way.
fn offset(silhouette: &[Loop2], distance: f64, sag: f64) -> OgeomResult<Vec<Loop2>> {
    let outward = distance > 0.0;
    let r = distance.abs() + sag / 2.0;
    // The widest chord sagging `sag` from a circle of radius `r`, and no
    // more than a quarter turn so that no sector is flat.
    let widest = if sag < r {
        (2.0 * (1.0 - sag / r).acos()).min(core::f64::consts::FRAC_PI_2)
    } else {
        core::f64::consts::FRAC_PI_2
    };
    let count = if outward { 0 } else { 1 };
    let mut edges = chain(silhouette);
    let mut piece = |ring: &[Point2]| {
        let n = ring.len();
        for i in 0..n {
            let mut weight = [0, 0];
            weight[count] = 1;
            edges.push(Edge2 {
                from: ring[i],
                to: ring[(i + 1) % n],
                weight,
            });
        }
    };
    for l in silhouette {
        let n = l.points.len();
        for i in 0..n {
            let (a, b, c) = (l.points[i], l.points[(i + 1) % n], l.points[(i + 2) % n]);
            let (u, v) = (b - a, c - b);
            let (lu, lv) = (u.magnitude(), v.magnitude());
            if lu == 0.0 || lv == 0.0 {
                continue;
            }
            let (u, v) = (u / lu, v / lv);
            // The right of travel is outside, the region being on the left.
            let side = |d: Vector2| {
                if outward {
                    -d.perpendicular()
                } else {
                    d.perpendicular()
                }
            };
            let nu = side(u);
            if outward {
                piece(&[a, a + nu * r, b + nu * r, b]);
            } else {
                piece(&[a, b, b + nu * r, a + nu * r]);
            }
            let turn = u.cross(v).atan2(u.dot(v));
            // A sector where the corner opens away from the region's side
            // being offset: convex corners outward, reflex ones inward.
            if (outward && turn > 0.0) || (!outward && turn < 0.0) {
                let steps = whole((turn.abs() / widest).ceil()).max(1);
                let mut ring = vec![b];
                for k in 0..=steps {
                    let t = turn * fraction(k, steps);
                    let (s, co) = t.sin_cos();
                    let d = Vector2::new(nu.x * co - nu.y * s, nu.x * s + nu.y * co);
                    ring.push(b + d * r);
                }
                if !outward {
                    ring[1..].reverse();
                }
                piece(&ring);
            }
        }
    }
    if outward {
        region(&edges, |w| w[0] > 0)
    } else {
        region(&edges, |w| w[0] > 0 && w[1] == 0)
    }
}

/// The band about a silhouette: which points lie between two signed
/// distances from it.
struct Band {
    segments: Vec<[Point2; 2]>,
    lo: f64,
    hi: f64,
    reach: f64,
    origin: Point2,
    cell: f64,
    columns: usize,
    rows: usize,
    start: Vec<usize>,
    held: Vec<usize>,
}

impl Band {
    fn new(silhouette: &[Loop2], lo: f64, hi: f64) -> Self {
        let mut segments = Vec::new();
        for l in silhouette {
            let n = l.points.len();
            for i in 0..n {
                segments.push([l.points[i], l.points[(i + 1) % n]]);
            }
        }
        let mut min = Point2::new(f64::INFINITY, f64::INFINITY);
        let mut max = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
        for p in segments.iter().flatten() {
            min = Point2::new(min.x.min(p.x), min.y.min(p.y));
            max = Point2::new(max.x.max(p.x), max.y.max(p.y));
        }
        let reach = lo.abs().max(hi.abs());
        let extent = (max.x - min.x).max(max.y - min.y);
        let per_side = fraction(segments.len(), 1).sqrt().ceil().max(1.0);
        let cell = (extent / per_side).max(reach / 2.0).max(f64::MIN_POSITIVE);
        let columns = whole((max.x - min.x) / cell) + 1;
        let rows = whole((max.y - min.y) / cell) + 1;
        let mut band = Self {
            segments,
            lo,
            hi,
            reach,
            origin: min,
            cell,
            columns,
            rows,
            start: Vec::new(),
            held: Vec::new(),
        };
        let mut pairs: Vec<(usize, usize)> = Vec::new();
        for (k, s) in band.segments.iter().enumerate() {
            let (c0, r0) = band.cell_of(Point2::new(s[0].x.min(s[1].x), s[0].y.min(s[1].y)));
            let (c1, r1) = band.cell_of(Point2::new(s[0].x.max(s[1].x), s[0].y.max(s[1].y)));
            for r in r0..=r1 {
                for c in c0..=c1 {
                    pairs.push((r * columns + c, k));
                }
            }
        }
        pairs.sort_unstable();
        let mut start = vec![0usize; rows * columns + 1];
        for &(b, _) in &pairs {
            start[b + 1] += 1;
        }
        for k in 0..rows * columns {
            start[k + 1] += start[k];
        }
        band.start = start;
        band.held = pairs.into_iter().map(|(_, k)| k).collect();
        band
    }

    fn cell_of(&self, p: Point2) -> (usize, usize) {
        let at = |v: f64, o: f64, n: usize| {
            let k = ((v - o) / self.cell).floor();
            if k < 0.0 {
                0
            } else {
                whole(k.min(fraction(n - 1, 1)))
            }
        };
        (
            at(p.x, self.origin.x, self.columns),
            at(p.y, self.origin.y, self.rows),
        )
    }

    fn cell(&self, c: usize, r: usize) -> &[usize] {
        let b = r * self.columns + c;
        &self.held[self.start[b]..self.start[b + 1]]
    }

    /// Whether the point's signed distance lies within the band.
    fn holds(&self, p: Point2) -> bool {
        let (c0, r0) = self.cell_of(Point2::new(p.x - self.reach, p.y - self.reach));
        let (c1, r1) = self.cell_of(Point2::new(p.x + self.reach, p.y + self.reach));
        let mut nearest = f64::INFINITY;
        for r in r0..=r1 {
            for c in c0..=c1 {
                for &k in self.cell(c, r) {
                    let [a, b] = self.segments[k];
                    nearest = nearest.min(distance_to_segment(p, a, b));
                }
            }
        }
        if nearest > self.reach {
            return false;
        }
        let signed = if self.inside(p) { -nearest } else { nearest };
        self.lo <= signed && signed <= self.hi
    }

    /// Whether the point is inside the silhouette, by the exact winding of
    /// the loops about it.
    fn inside(&self, p: Point2) -> bool {
        let (c0, r) = self.cell_of(p);
        let mut seen: Vec<usize> = (c0..self.columns)
            .flat_map(|c| self.cell(c, r).iter().copied())
            .collect();
        seen.sort_unstable();
        seen.dedup();
        let at = [p.x, p.y];
        let mut winding = 0i32;
        for k in seen {
            let [a, b] = self.segments[k];
            let side = Exact::orient2d([a.x, a.y], [b.x, b.y], at);
            if a.y <= p.y && p.y < b.y && side == Sign::Positive {
                winding += 1;
            } else if b.y <= p.y && p.y < a.y && side == Sign::Negative {
                winding -= 1;
            }
        }
        winding > 0
    }

    /// Whether every sample of a stretch holds, the samples a sixteenth
    /// closer than `spacing` so that the difference covers the rounding of
    /// the distances.
    fn holds_along(&self, point: impl Fn(f64) -> Point2, length: f64, spacing: f64) -> bool {
        let n = (length / (spacing * 15.0 / 16.0)).ceil().max(1.0);
        if n > 1e8 {
            return false;
        }
        let n = whole(n);
        (0..=n).all(|k| self.holds(point(fraction(k, n))))
    }

    fn holds_piece(&self, piece: &Piece2, spacing: f64) -> bool {
        self.holds_along(|t| piece.point_at(t), piece.length(), spacing)
    }
}

fn distance_to_segment(p: Point2, a: Point2, b: Point2) -> f64 {
    let d = b - a;
    let len2 = d.dot(d);
    let t = if len2 > 0.0 {
        ((p - a).dot(d) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (a + d * t).distance(p)
}

/// A curve a piece runs along, in its direction of travel.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Path {
    /// `at + t along`, `along` a unit vector.
    Line { at: Point2, along: Vector2 },
    /// `centre + radius (cos t, sin t)`, `t` rising if counter-clockwise.
    Circle {
        centre: Point2,
        radius: f64,
        ccw: bool,
    },
}

impl Path {
    fn point(self, t: f64) -> Point2 {
        match self {
            Self::Line { at, along } => at + along * t,
            Self::Circle { centre, radius, .. } => on_circle(centre, radius, t),
        }
    }

    fn tangent(self, t: f64) -> Vector2 {
        match self {
            Self::Line { along, .. } => along,
            Self::Circle { ccw, .. } => {
                let (s, c) = t.sin_cos();
                if ccw {
                    Vector2::new(-s, c)
                } else {
                    Vector2::new(s, -c)
                }
            }
        }
    }

    /// The parameter of a point on the path, an angle taken nearest `near`.
    fn parameter(self, p: Point2, near: f64) -> f64 {
        match self {
            Self::Line { at, along } => (p - at).dot(along),
            Self::Circle { centre, .. } => {
                let a = (p.y - centre.y).atan2(p.x - centre.x);
                let tau = core::f64::consts::TAU;
                a + ((near - a) / tau).round() * tau
            }
        }
    }

    /// The nearest point of the path's curve.
    pub(crate) fn foot(self, p: Point2) -> Option<Point2> {
        match self {
            Self::Line { at, along } => Some(at + along * (p - at).dot(along)),
            Self::Circle { centre, radius, .. } => {
                let d = p - centre;
                let len = d.magnitude();
                (len > 0.0).then(|| centre + d * (radius / len))
            }
        }
    }

    /// The path moved `by` to the left of travel.
    fn moved(self, by: f64) -> Option<Self> {
        match self {
            Self::Line { at, along } => Some(Self::Line {
                at: at + along.perpendicular() * by,
                along,
            }),
            Self::Circle {
                centre,
                radius,
                ccw,
            } => {
                let radius = if ccw { radius - by } else { radius + by };
                (radius > 0.0).then_some(Self::Circle {
                    centre,
                    radius,
                    ccw,
                })
            }
        }
    }
}

/// Where two paths' curves cross.
pub(crate) fn meet(a: Path, b: Path) -> Vec<Point2> {
    match (a, b) {
        (Path::Line { at: p, along: u }, Path::Line { at: q, along: v }) => {
            let den = u.cross(v);
            if den == 0.0 {
                return Vec::new();
            }
            vec![p + u * ((q - p).cross(v) / den)]
        }
        (Path::Line { at, along }, Path::Circle { centre, radius, .. })
        | (Path::Circle { centre, radius, .. }, Path::Line { at, along }) => {
            let foot = at + along * (centre - at).dot(along);
            let d2 = radius.mul_add(radius, -foot.square_distance(centre));
            if d2 < 0.0 {
                return Vec::new();
            }
            let h = d2.sqrt();
            vec![foot + along * h, foot - along * h]
        }
        (
            Path::Circle {
                centre: c0,
                radius: r0,
                ..
            },
            Path::Circle {
                centre: c1,
                radius: r1,
                ..
            },
        ) => {
            let d = c1 - c0;
            let len = d.magnitude();
            if len == 0.0 {
                return Vec::new();
            }
            let along = (len.mul_add(len, r0.mul_add(r0, -r1 * r1))) / (2.0 * len);
            let h2 = r0.mul_add(r0, -along * along);
            if h2 < 0.0 {
                return Vec::new();
            }
            let u = d / len;
            let mid = c0 + u * along;
            let h = h2.sqrt();
            vec![mid + u.perpendicular() * h, mid - u.perpendicular() * h]
        }
    }
}

fn nearest_to(points: &[Point2], p: Point2) -> Option<Point2> {
    points
        .iter()
        .copied()
        .min_by(|a, b| a.square_distance(p).total_cmp(&b.square_distance(p)))
}

/// A run's path and the parameters its fit spans.
#[derive(Debug, Clone, Copy)]
struct Span {
    path: Path,
    from: f64,
    to: f64,
}

impl Span {
    /// How far along the span's own range a parameter lies: 0 at its
    /// start, 1 at its end.
    fn progress(&self, t: f64) -> f64 {
        (t - self.from) / (self.to - self.from)
    }

    fn length(&self) -> f64 {
        match self.path {
            Path::Line { .. } => (self.to - self.from).abs(),
            Path::Circle { radius, .. } => radius * (self.to - self.from).abs(),
        }
    }

    /// The piece between two parameters.
    fn piece(&self, from: f64, to: f64) -> Piece2 {
        match self.path {
            Path::Line { .. } => Piece2::Line {
                from: self.path.point(from),
                to: self.path.point(to),
            },
            Path::Circle { centre, radius, .. } => Piece2::Arc {
                centre,
                radius,
                start: from,
                sweep: to - from,
            },
        }
    }
}

/// How two consecutive pieces meet.
#[derive(Debug, Clone, Copy)]
enum Joint {
    /// The first ends at `a`, the second starts at `b`, both at `at`.
    Corner { a: f64, b: f64, at: Point2 },
    /// The first ends at `a`, the arc follows, the second starts at `b`.
    Round { a: f64, b: f64, arc: Piece2 },
    /// Each keeps its own end and a line joins them.
    Bridge,
}

struct Fit<'a> {
    band: &'a Band,
    within: f64,
    round_below: f64,
    spacing: f64,
}

impl Fit<'_> {
    /// The outline of one loop of the band's middle.
    fn wire(&self, target: &Loop2) -> OgeomResult<Wire2d> {
        if let Some(w) = self.fitted(target)
            && w.pieces
                .iter()
                .all(|p| self.band.holds_piece(p, self.spacing))
        {
            return Ok(w);
        }
        let n = target.points.len();
        let polygon = Wire2d {
            pieces: (0..n)
                .map(|i| Piece2::Line {
                    from: target.points[i],
                    to: target.points[(i + 1) % n],
                })
                .collect(),
        };
        if polygon
            .pieces
            .iter()
            .all(|p| self.band.holds_piece(p, self.spacing))
        {
            return Ok(polygon);
        }
        ogeom_bail!(
            Construction,
            "the middle of the band leaves the band; the outline is not fitted"
        )
    }

    /// The pieces fitted to a loop, before the check.
    fn fitted(&self, target: &Loop2) -> Option<Wire2d> {
        let points = dense(target, self.within);
        if points.len() < 4 || self.within.is_nan() || self.within <= 0.0 {
            return None;
        }
        if let Some(Shape::Circle { centre, radius }) = circle_through(&points[..points.len() - 1])
            && deviation(&points, Shape::Circle { centre, radius }) <= self.within
        {
            let circle = Piece2::Circle {
                centre,
                radius,
                counter_clockwise: target.area() > 0.0,
            };
            if self.band.holds_piece(&circle, self.spacing) {
                return Some(Wire2d {
                    pieces: vec![circle],
                });
            }
        }
        let runs = runs(&points, self.within);
        if runs.len() < 2 {
            return None;
        }
        let mut spans = Vec::with_capacity(runs.len());
        for run in &runs {
            let (p, q) = (points[run.first], points[run.last]);
            let span = match run.shape {
                Shape::Line { at, along } => Span {
                    path: Path::Line { at, along },
                    from: (p - at).dot(along),
                    to: (q - at).dot(along),
                },
                Shape::Circle { centre, radius } => {
                    let angle = |x: Point2| (x.y - centre.y).atan2(x.x - centre.x);
                    let mut sweep = 0.0;
                    for k in run.first..run.last {
                        let (u, v) = (points[k] - centre, points[k + 1] - centre);
                        sweep += u.cross(v).atan2(u.dot(v));
                    }
                    let from = angle(p);
                    Span {
                        path: Path::Circle {
                            centre,
                            radius,
                            ccw: sweep > 0.0,
                        },
                        from,
                        to: from + sweep,
                    }
                }
            };
            let forward = match span.path {
                Path::Line { .. } => span.to > span.from,
                Path::Circle { .. } => span.to != span.from,
            };
            if !forward {
                return None;
            }
            spans.push(span);
        }
        // Each span with the loop point it shares with the next.
        let mut chain: Vec<(Span, Point2)> = spans
            .into_iter()
            .zip(runs.iter().map(|r| points[r.last]))
            .collect();
        let mut joints: Vec<Option<Joint>> = vec![None; chain.len()];
        // A piece between two longer lines goes where the lines themselves
        // meet at a corner or round one off: a corner's own rounding is the
        // largest arc that fits, which the run fit only approximates by an
        // arc or cuts off by a chamfer.
        let mut k = 0;
        while chain.len() > 3 && k < chain.len() {
            let n = chain.len();
            let (before, after) = ((k + n - 1) % n, (k + 1) % n);
            let is_line = |s: &Span| matches!(s.path, Path::Line { .. });
            let shorter =
                chain[k].0.length() < chain[before].0.length().min(chain[after].0.length());
            if shorter && is_line(&chain[before].0) && is_line(&chain[after].0) {
                let arc = &chain[k].0;
                let middle = arc.path.point(0.5 * (arc.from + arc.to));
                let joint = self.joint(&chain[before].0, &chain[after].0, middle);
                if !matches!(joint, Joint::Bridge) {
                    chain.remove(k);
                    joints.remove(k);
                    let before = if k == 0 { chain.len() - 1 } else { k - 1 };
                    joints[before] = Some(joint);
                    continue;
                }
            }
            k += 1;
        }
        let count = chain.len();
        let joints: Vec<Joint> = (0..count)
            .map(|k| {
                joints[k].unwrap_or_else(|| {
                    self.joint(&chain[k].0, &chain[(k + 1) % count].0, chain[k].1)
                })
            })
            .collect();
        let spans: Vec<Span> = chain.into_iter().map(|(s, _)| s).collect();
        let mut pieces = Vec::with_capacity(2 * count);
        for k in 0..count {
            let span = &spans[k];
            let from = match joints[(k + count - 1) % count] {
                Joint::Corner { b, .. } | Joint::Round { b, .. } => b,
                Joint::Bridge => span.from,
            };
            let (to, after) = match joints[k] {
                Joint::Corner { a, .. } => (a, None),
                Joint::Round { a, arc, .. } => (a, Some(arc)),
                Joint::Bridge => {
                    let next = &spans[(k + 1) % count];
                    let line = Piece2::Line {
                        from: span.path.point(span.to),
                        to: next.path.point(next.from),
                    };
                    (span.to, Some(line))
                }
            };
            let mut piece = span.piece(from, to);
            // Corners meet at one computed point, shared by both lines.
            if let Piece2::Line {
                from: start,
                to: end,
            } = &mut piece
            {
                if let Joint::Corner { at, .. } = joints[(k + count - 1) % count] {
                    *start = at;
                }
                if let Joint::Corner { at, .. } = joints[k] {
                    *end = at;
                }
            }
            if piece.length() > 0.0 {
                pieces.push(piece);
            }
            if let Some(p) = after
                && p.length() > 0.0
            {
                pieces.push(p);
            }
        }
        Some(Wire2d { pieces })
    }

    /// How the piece on `a` meets the piece on `b`, both fitted through
    /// the loop point `q`.
    fn joint(&self, a: &Span, b: &Span, q: Point2) -> Joint {
        let corner = nearest_to(&meet(a.path, b.path), q).and_then(|x| {
            let (ta, tb) = (a.path.parameter(x, a.to), b.path.parameter(x, b.from));
            (a.progress(ta) >= 0.5 && b.progress(tb) <= 0.5).then_some((ta, tb, x))
        });
        let (ta, tb) = corner.map_or((a.to, b.from), |(ta, tb, _)| (ta, tb));
        let (u, v) = (a.path.tangent(ta), b.path.tangent(tb));
        let turn = u.cross(v).atan2(u.dot(v));
        let sharp = corner.filter(|&(ta, tb, _)| self.extension_holds(a, b, ta, tb));
        let at_corner = sharp.map(|(a, b, at)| Joint::Corner { a, b, at });
        if turn.abs() >= self.round_below
            && let Some(j) = at_corner
        {
            return j;
        }
        if turn != 0.0
            && let Some(j) = self.round(a, b, turn, corner.map_or(q, |c| c.2))
        {
            return j;
        }
        at_corner.unwrap_or(Joint::Bridge)
    }

    /// Whether the stretches by which a joint lengthens either piece
    /// beyond its fitted range hold.
    fn extension_holds(&self, a: &Span, b: &Span, ta: f64, tb: f64) -> bool {
        let stretch = |s: &Span, t0: f64, t1: f64| {
            let length = s.piece(t0, t1).length();
            self.band
                .holds_along(|t| s.path.point(t0 + (t1 - t0) * t), length, self.spacing)
        };
        (a.progress(ta) <= 1.0 || stretch(a, a.to, ta))
            && (b.progress(tb) >= 0.0 || stretch(b, tb, b.from))
    }

    /// The tangent arc of radius `r` from the piece on `a` to the piece on
    /// `b` about a centre on the inside of the turn, nearest `near`.
    fn arc(&self, a: &Span, b: &Span, r: f64, turn: f64, near: Point2) -> Option<Joint> {
        let side = turn.signum();
        let centre = nearest_to(
            &meet(a.path.moved(side * r)?, b.path.moved(side * r)?),
            near,
        )?;
        let (pa, pb) = (a.path.foot(centre)?, b.path.foot(centre)?);
        let (ta, tb) = (a.path.parameter(pa, a.to), b.path.parameter(pb, b.from));
        if a.progress(ta) < 0.5 || b.progress(tb) > 0.5 {
            return None;
        }
        let (u, v) = (pa - centre, pb - centre);
        let sweep = u.cross(v).atan2(u.dot(v));
        if sweep.signum() != side {
            return None;
        }
        // The arc leaves the first piece and enters the second along them.
        let start = u.y.atan2(u.x);
        let arc = Piece2::Arc {
            centre,
            radius: r,
            start,
            sweep,
        };
        let along = |t: f64| {
            let d = on_circle(Point2::ORIGIN, 1.0, sweep.mul_add(t, start)).to_vector();
            d.perpendicular() * side
        };
        if along(0.0).dot(a.path.tangent(ta)) <= 0.0 || along(1.0).dot(b.path.tangent(tb)) <= 0.0 {
            return None;
        }
        Some(Joint::Round { a: ta, b: tb, arc })
    }

    fn arc_holds(&self, a: &Span, b: &Span, joint: &Joint) -> bool {
        match joint {
            Joint::Round { a: ta, b: tb, arc } => {
                self.band.holds_piece(arc, self.spacing) && self.extension_holds(a, b, *ta, *tb)
            }
            _ => false,
        }
    }

    /// The largest tangent arc found to fit the band: radii shrinking by
    /// a factor from one leaving each piece half its length, then halving
    /// the factor between the last that failed and the first that held.
    fn round(&self, a: &Span, b: &Span, turn: f64, near: Point2) -> Option<Joint> {
        let reach = 0.5 * a.length().min(b.length());
        let largest = reach / (turn.abs() / 2.0).tan().max(1e-9);
        let shrink: f64 = 0.7;
        let mut failed = None;
        for k in 0..80 {
            let r = largest * shrink.powi(k);
            let Some(joint) = self.arc(a, b, r, turn, near) else {
                failed = Some(r);
                continue;
            };
            if !self.arc_holds(a, b, &joint) {
                failed = Some(r);
                continue;
            }
            let Some(mut above) = failed else {
                return Some(joint);
            };
            let (mut below, mut best) = (r, joint);
            for _ in 0..12 {
                let mid = (above * below).sqrt();
                match self.arc(a, b, mid, turn, near) {
                    Some(j) if self.arc_holds(a, b, &j) => {
                        below = mid;
                        best = j;
                    }
                    _ => above = mid,
                }
            }
            return Some(best);
        }
        None
    }
}

/// A loop's points from its sharpest corner round and back to it, with
/// points added so that none are farther apart than `spacing`, as the run
/// fit wants them.
fn dense(target: &Loop2, spacing: f64) -> Vec<Point2> {
    let p = &target.points;
    let n = p.len();
    let turn = |i: usize| {
        let (u, v) = (p[i] - p[(i + n - 1) % n], p[(i + 1) % n] - p[i]);
        u.cross(v).atan2(u.dot(v)).abs()
    };
    let first = (0..n)
        .max_by(|&i, &j| turn(i).total_cmp(&turn(j)))
        .unwrap_or(0);
    let perimeter: f64 = (0..n).map(|i| p[i].distance(p[(i + 1) % n])).sum();
    // Two million points at most.
    let spacing = spacing.max(perimeter / 2e6);
    let mut out = Vec::new();
    for k in 0..n {
        let (a, b) = (p[(first + k) % n], p[(first + k + 1) % n]);
        let steps = whole((a.distance(b) / spacing).ceil()).max(1);
        for s in 0..steps {
            out.push(a.lerp(b, fraction(s, steps)));
        }
    }
    out.push(p[first]);
    out
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;

    fn square(s: f64) -> Loop2 {
        Loop2 {
            points: vec![
                Point2::new(0.0, 0.0),
                Point2::new(s, 0.0),
                Point2::new(s, s),
                Point2::new(0.0, s),
            ],
        }
    }

    /// A square with no play and a tolerance wide enough for its corners
    /// comes back as four lines meeting at corners, one tolerance out.
    #[test]
    fn a_square_keeps_its_corners() {
        let o = Outline {
            play: 0.0,
            tolerance: 0.1,
            round_below: 0.5,
        };
        let wires = fit_outline(&[square(10.0)], &o).unwrap();
        assert_eq!(wires.len(), 1);
        let w = &wires[0];
        assert_eq!(w.pieces.len(), 4, "{w:?}");
        for p in &w.pieces {
            assert!(matches!(p, Piece2::Line { .. }), "{w:?}");
        }
    }

    /// With a play too large for a sharp corner to stay in the band, the
    /// corners round.
    #[test]
    fn a_large_play_rounds_the_corners() {
        let o = Outline {
            play: 1.0,
            tolerance: 0.1,
            round_below: 0.5,
        };
        let wires = fit_outline(&[square(10.0)], &o).unwrap();
        assert_eq!(wires.len(), 1);
        let arcs = wires[0]
            .pieces
            .iter()
            .filter(|p| matches!(p, Piece2::Arc { .. }))
            .count();
        assert_eq!(arcs, 4, "{:?}", wires[0]);
    }
}
