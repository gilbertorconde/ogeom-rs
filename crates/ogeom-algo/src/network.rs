//! A network of B-spline patches fitted to a scan at once.
//!
//! The caller lays out the network: four-sided patches, each bounded by
//! curves drawn on the scan (a polyline, or a curve), a curve between two
//! patches shared by both. [`fit_network`] fits every patch to the points
//! in one least-squares system, in which:
//!
//! - **a shared curve is one set of poles.** Both patches' boundary rows
//!   are the same unknowns, and the corners where curves meet are one
//!   unknown each, so the patches meet exactly, with no gap to report but
//!   rounding.
//! - **a curve that is not a crease is held tangent, pole by pole.** Each
//!   pole of the curve and the poles beside it in the two patches are
//!   held on one line, at a ratio fixed along the curve:
//!   `B − L = λ (R − B)`. The cross derivatives of the two patches along
//!   the curve are then parallel everywhere, the same tangent plane on
//!   either side. The ratio is read off a fit without the constraints, as
//!   the ratio of the pole spacings either side, and is shared by every
//!   curve that runs straight on through a corner of four patches, which
//!   is what lets the constraints of the curves crossing there agree.
//! - **a corner of other than four patches is left near-tangent.** The
//!   triples hold the parameter lines that cross a curve to one direction
//!   on both sides of it. At a corner of three patches (or five, or the
//!   point a three-sided patch's missing side shrinks to) the lines
//!   crossing a curve there lean one way on one side and the other way on
//!   the other, and no triple holds; tangency there needs the cross
//!   derivatives to differ by a multiple of the curve's own tangent, which
//!   is not a pole-by-pole relation on a curve of several spans. So a curve
//!   with such a corner at an end is fitted without the triples: near
//!   tangent as the points are smooth, its angle measured into
//!   [`Network::worst_tangent`], and the widest angle between the
//!   patches' normals at the corner reported in [`Network::corners`].
//! - **spans are raised per strip.** The curves opposite each other in a
//!   patch share knots, and so do all the curves of a strip of patches
//!   joined through their opposite sides. Where a patch misses the points
//!   by more than the tolerance its strips gain spans, the way
//!   [`fit_surface_scattered_over`](ogeom_geom::fit::fit_surface_scattered_over)
//!   refines a single patch.
//!
//! The constraints are solved for exactly, by elimination, before the
//! least squares: what is held tangent is tangent to rounding.

mod build;
mod solve;
mod topology;

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::fit::Ended;
use ogeom_geom::{BSplineSurface, Curve, SurfaceGeometry, SurfaceSeeds, project_on_surface_from};
use ogeom_math::{ControlGrid, KnotVector, Point, Vector};
use ogeom_mesh::PointIndex;
use ogeom_topo::{Model, Shape};

use solve::{Basis, Lower, null_space, solve_spd};
use topology::Topology;

/// The degree of every patch and curve: the least that is curvature
/// continuous within a patch and carries the thin-plate fairing.
const DEGREE: usize = 3;

/// The upper quartile of the standard normal distribution, `Φ⁻¹(3/4)`: the
/// median of the absolute value of a normal deviate, in deviations.
const NORMAL_QUARTILE: f64 = 0.674_489_750_196_081_7;

/// How many points a worker projects at a time.
const POINT_BATCH: usize = 2048;

/// One side of a patch: the curve along it and which way it runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PatchSide {
    /// The curve's index in the network's curves.
    pub curve: usize,
    /// Whether the curve runs against the patch's walk round its sides.
    pub reversed: bool,
}

/// A patch of the network: its four sides, walked round it.
///
/// Side 0 runs along `u` at `v = 0`, side 1 up `v` at `u = 1`, side 2 back
/// along `u` at `v = 1`, side 3 back down `v` at `u = 0`, so the walk turns
/// counter-clockwise seen from the side the patch's normal `S_u × S_v`
/// faces. Two patches sharing a curve walk it in opposite directions, as
/// the faces of a shell do. One side may be missing: a three-sided patch,
/// that side shrunk to the point where its neighbours meet.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PatchDef {
    /// The sides, in walking order.
    pub sides: [Option<PatchSide>; 4],
}

impl PatchDef {
    /// A four-sided patch from its curves in walking order, each with
    /// whether it runs against the walk.
    #[must_use]
    pub fn four(sides: [(usize, bool); 4]) -> Self {
        Self {
            sides: sides.map(|(curve, reversed)| Some(PatchSide { curve, reversed })),
        }
    }
}

/// A curve of the network, as the caller draws it on the scan.
#[derive(Debug, Clone, PartialEq)]
pub enum CurveDef {
    /// Points along the curve, in order, the first and last at its ends.
    Polyline(Vec<Point>),
    /// A curve over a range of its parameter, read within the tolerance.
    Curve {
        /// The curve.
        curve: Curve,
        /// The range it is used over.
        range: (f64, f64),
    },
}

/// How [`fit_network`] fits.
#[derive(Debug, Clone, PartialEq)]
pub struct NetworkOptions {
    /// The distance every point is to come within; also how far the
    /// curves meeting at a corner may stand apart there, and the chord a
    /// [`CurveDef::Curve`] is read to.
    pub tolerance: f64,
    /// The thin-plate fairing weight of each patch, over the unit square
    /// its parameters run over, against the mean squared distance of its
    /// points, as
    /// [`ScatteredOptions::fairing`](ogeom_geom::fit::ScatteredOptions::fairing)
    /// weighs a single patch. It settles the poles no point reaches: a
    /// patch with no points needs it above nought.
    pub smoothing: f64,
    /// The curves left creased: their patches share the curve but not a
    /// tangent plane across it.
    pub creases: Vec<usize>,
    /// The knot spans every strip starts with.
    pub spans: usize,
    /// The most knot spans a strip may reach.
    pub most: usize,
}

/// How one patch fits its points.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PatchFit {
    /// The root mean square distance from its points to it; nought with
    /// no points.
    pub sigma: f64,
    /// The farthest of its points from it.
    pub worst: f64,
    /// How many points it fits.
    pub points: usize,
    /// Its knot spans, `u` then `v`.
    pub spans: (usize, usize),
}

/// A corner inside the network where other than four patches meet, every
/// curve at it held tangent.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IrregularCorner {
    /// Where it is.
    pub point: Point,
    /// How many patches meet there.
    pub valence: usize,
    /// The widest angle between two of their normals there, in radians.
    pub angle: f64,
}

/// What [`fit_network`] produced.
#[derive(Debug, Clone)]
pub struct Network {
    /// A face per patch, in the patches' order, sharing an edge per shared
    /// curve.
    pub faces: Vec<Shape>,
    /// The faces as one shell.
    pub shell: Shape,
    /// An edge per curve, in the curves' order.
    pub edges: Vec<Shape>,
    /// Each patch's surface, over the unit square.
    pub surfaces: Vec<BSplineSurface>,
    /// The widest angle, in radians, between the two patches' normals
    /// across a curve between two patches that is not a crease, read at
    /// both ends and the Gauss nodes of every span (an end at a missing
    /// side's point left out, where a normal is not defined), and at the
    /// irregular corners.
    pub worst_tangent: f64,
    /// The widest distance between a curve's edge and the boundary of a
    /// patch it bounds, read at the same nodes.
    pub worst_gap: f64,
    /// How each patch fits its points.
    pub per_patch: Vec<PatchFit>,
    /// The corners inside the network where other than four patches meet,
    /// no curve there a crease, and their angles.
    pub corners: Vec<IrregularCorner>,
    /// How many points lie past the network's border, nearest a side no
    /// other patch shares, and were left out.
    pub outside: usize,
    /// Why refinement stopped; never [`Ended::Fixed`].
    pub ended: Ended,
}

/// Fit a network of patches to a scan's points, all at once.
///
/// Each patch is first seeded: every curve is drawn as a B-spline with its
/// poles on the caller's curve (at the Greville abscissae of its chord
/// parameters), its spans doubled from [`NetworkOptions::spans`] until it
/// passes within the tolerance of the caller's points or reaches
/// [`NetworkOptions::most`], and each patch is the Coons patch of its four
/// seeds. Each point is taken by the seed nearest it, at the parameters
/// of its nearest point there; a point whose nearest point lies on a side
/// no other patch shares lies past the network and is left out. The points
/// are those of the region the network covers: a point of the scan beyond
/// it is fitted by the nearest patch unless it lies past a free side.
///
/// The patches are then fitted as the [module documentation](self)
/// describes, each weighing the mean squared distance of its points plus
/// [`NetworkOptions::smoothing`] times its thin-plate energy, and every
/// distance reported is from a point to its nearest point on its patch.
/// While a point misses by more than the tolerance, each strip of a patch
/// that misses is tried at twice its spans, and the try Mallows' criterion
/// scores best is kept if it scores below the fit so far, as
/// [`fit_surface_scattered_over`](ogeom_geom::fit::fit_surface_scattered_over)
/// refines, the controls counted after the constraints. Every such strip is
/// tried, not the worst patch's alone: a patch held tangent to a neighbour
/// too coarse for the shape misses for its neighbour. A patch misses only
/// in a cell holding as many points as the sixteen poles that bear on it,
/// and a point farther off than the noise of its cell reaches is a stray
/// rather than a miss.
///
/// The faces are built in `model` on the fitted patches, every curve an
/// edge on the poles the patches share, with the parameter lines of the
/// patches as its pcurves, and sewn by those shared edges into one shell.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// the tolerance is not a distance, the smoothing is not finite and
/// non-negative, the spans are nought or above their bound, there is no
/// patch, a patch has more than one side missing, names a curve that does
/// not exist or is bounded twice by one curve or meets itself at a
/// corner, a curve bounds no patch or more than two, two patches walk a
/// shared curve the same way, a curve closes on itself, a crease names no
/// curve, a curve has fewer than two points or a point that is not
/// finite, or curves meeting at a corner stand farther apart there than
/// the tolerance; [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric)
/// if the points and the fairing cannot settle every pole (raise the
/// smoothing); as [`project_on_surface`](ogeom_geom::project_on_surface)
/// and the face builders otherwise.
pub fn fit_network(
    model: &mut Model,
    patches: &[PatchDef],
    curves: &[CurveDef],
    points: &PointIndex,
    options: &NetworkOptions,
    tol: Tolerances,
) -> OgeomResult<Network> {
    let tolerance = options.tolerance;
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        ogeom_bail!(Construction, "a tolerance of {tolerance} is not a distance");
    }
    if !(options.smoothing.is_finite() && options.smoothing >= 0.0) {
        ogeom_bail!(
            Construction,
            "a smoothing weight of {} is not a weight",
            options.smoothing
        );
    }
    if options.spans == 0 || options.most < options.spans {
        ogeom_bail!(
            Construction,
            "{} starting spans within {} is not a span budget",
            options.spans,
            options.most
        );
    }
    let topo = Topology::of(patches, curves.len(), &options.creases)?;
    let drawn = drawn_curves(curves, tolerance, tol)?;
    let corners = corner_points(&topo, &drawn, tolerance)?;
    let seeds = Seeds::new(&topo, patches, &drawn, &corners, options, tol)?;
    let cloud = points.points();
    let taken = seeds.take(&topo, patches, &cloud, tol)?;

    let fit = Fit {
        topo: &topo,
        patches,
        taken: &taken,
        smoothing: options.smoothing,
        tol,
    };
    let mut spans = vec![options.spans; topo.strips];
    let mut current = fit.round(&spans)?;
    let ended = loop {
        if current.distances.iter().all(|d| *d <= tolerance) {
            break Ended::Met;
        }
        let (missed, any, strays) = fit.missed(&current, tolerance, (DEGREE + 1) * (DEGREE + 1));
        if !any && strays {
            break Ended::Strays;
        }
        // The strips of every patch that misses, each tried alone: a patch
        // can miss for a neighbour too coarse to follow the shape, which
        // its tangency binds it to.
        let mut tries: Vec<usize> = missed
            .iter()
            .flat_map(|&(_, p)| {
                let (su, sv) = topo.patch_strips[p];
                [su, sv]
            })
            .filter(|&s| spans[s] < options.most)
            .collect();
        tries.sort_unstable();
        tries.dedup();
        if tries.is_empty() {
            break Ended::Bounds;
        }
        // Mallows' criterion, as the single-patch refinement scores: the
        // squared distances plus twice the noise's variance per free pole,
        // the variance read off the median distance.
        let variance = (median(&current.distances) / NORMAL_QUARTILE).powi(2);
        let score = |r: &Round| {
            let square: f64 = r.distances.iter().map(|d| d * d).sum();
            2.0f64.mul_add(variance * precise(r.free), square)
        };
        let mut best: Option<(f64, Vec<usize>, Round)> = None;
        for raised in tries {
            let mut next = spans.clone();
            next[raised] = (2 * next[raised]).min(options.most);
            let tried = fit.round(&next)?;
            let scored = score(&tried);
            if best.as_ref().is_none_or(|(s, _, _)| scored < *s) {
                best = Some((scored, next, tried));
            }
        }
        match best {
            Some((scored, next, tried)) if scored < score(&current) => {
                spans = next;
                current = tried;
            }
            _ => break Ended::Noise,
        }
    };

    let built = build::faces(model, &topo, patches, &current, tol)?;
    let measured = build::measure(&topo, patches, &current, tol)?;
    let mut per_patch = vec![
        PatchFit {
            sigma: 0.0,
            worst: 0.0,
            points: 0,
            spans: (0, 0),
        };
        patches.len()
    ];
    for (p, fit) in per_patch.iter_mut().enumerate() {
        let (su, sv) = topo.patch_strips[p];
        fit.spans = (spans[su], spans[sv]);
    }
    for (&p, &d) in taken.patch.iter().zip(&current.distances) {
        let fit = &mut per_patch[p];
        fit.points += 1;
        fit.sigma += d * d;
        fit.worst = fit.worst.max(d);
    }
    for fit in &mut per_patch {
        if fit.points > 0 {
            fit.sigma = (fit.sigma / precise(fit.points)).sqrt();
        }
    }
    Ok(Network {
        faces: built.faces,
        shell: built.shell,
        edges: built.edges,
        surfaces: current.surfaces,
        worst_tangent: measured.worst_tangent,
        worst_gap: measured.worst_gap,
        per_patch,
        corners: measured.corners,
        outside: taken.outside,
        ended,
    })
}

/// A count as a float.
#[allow(clippy::cast_precision_loss, reason = "counts far below 2^52")]
fn precise(n: usize) -> f64 {
    n as f64
}

/// The median of non-empty `values`: the mean of the middle two of an even
/// count.
fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let k = sorted.len() / 2;
    if sorted.len().is_multiple_of(2) {
        f64::midpoint(sorted[k - 1], sorted[k])
    } else {
        sorted[k]
    }
}

/// Each curve as points along it.
fn drawn_curves(
    curves: &[CurveDef],
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<Vec<Point>>> {
    let mut out = Vec::with_capacity(curves.len());
    for (c, def) in curves.iter().enumerate() {
        let points = match def {
            CurveDef::Polyline(points) => points.clone(),
            CurveDef::Curve { curve, range } => {
                let deflection = ogeom_mesh::Deflection {
                    chord: tolerance,
                    ..ogeom_mesh::Deflection::default()
                };
                ogeom_mesh::discretize(curve, *range, deflection, tol)?.points
            }
        };
        if points.len() < 2 {
            ogeom_bail!(Construction, "curve {c} has fewer than two points");
        }
        if let Some(p) = points
            .iter()
            .find(|p| !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()))
        {
            ogeom_bail!(Construction, "curve {c}'s point {p:?} is not finite");
        }
        out.push(points);
    }
    Ok(out)
}

/// Each vertex: the mean of the curve ends that meet there, all within the
/// tolerance of it.
fn corner_points(topo: &Topology, drawn: &[Vec<Point>], tolerance: f64) -> OgeomResult<Vec<Point>> {
    let mut sums = vec![(Vector::ZERO, 0usize); topo.vertices];
    let ends = |c: usize| {
        let d = &drawn[c];
        (d[0], d[d.len() - 1])
    };
    for (c, &(a, b)) in topo.curve_ends.iter().enumerate() {
        let (pa, pb) = ends(c);
        sums[a].0 += pa.to_vector();
        sums[a].1 += 1;
        sums[b].0 += pb.to_vector();
        sums[b].1 += 1;
    }
    let corners: Vec<Point> = sums
        .iter()
        .map(|(s, n)| Point::ORIGIN + *s / precise((*n).max(1)))
        .collect();
    for (c, &(a, b)) in topo.curve_ends.iter().enumerate() {
        let (pa, pb) = ends(c);
        for (v, p) in [(a, pa), (b, pb)] {
            let off = p.distance(corners[v]);
            if off > tolerance {
                ogeom_bail!(
                    Construction,
                    "curve {c} ends {off} from the corner it meets other curves at, \
                     beyond the tolerance"
                );
            }
        }
    }
    Ok(corners)
}

/// A knot vector of `spans` even spans over the unit interval.
fn even_knots(spans: usize) -> OgeomResult<KnotVector> {
    KnotVector::clamped_uniform(DEGREE, spans + DEGREE)
}

/// A knot vector's Greville abscissae: where each pole sits along it.
fn greville(knots: &KnotVector) -> Vec<f64> {
    let t = knots.knots();
    (0..knots.control_point_count())
        .map(|i| t[i + 1..=i + DEGREE].iter().sum::<f64>() / precise(DEGREE))
        .collect()
}

/// Where the poles of the patches' nets stand among the unknowns: a slot
/// per vertex, per inner pole of each curve and per inner pole of each
/// patch.
struct Layout {
    /// The poles along each strip.
    controls: Vec<usize>,
    /// Each curve's first inner pole's slot.
    curve_base: Vec<usize>,
    /// Each patch's first inner pole's slot.
    patch_base: Vec<usize>,
    /// How many slots there are.
    slots: usize,
}

impl Layout {
    fn new(topo: &Topology, spans: &[usize]) -> Self {
        let controls: Vec<usize> = spans.iter().map(|s| s + DEGREE).collect();
        let mut next = topo.vertices;
        let curve_base = topo
            .curve_strip
            .iter()
            .map(|&s| {
                let base = next;
                next += controls[s] - 2;
                base
            })
            .collect();
        let patch_base = topo
            .patch_strips
            .iter()
            .map(|&(su, sv)| {
                let base = next;
                next += (controls[su] - 2) * (controls[sv] - 2);
                base
            })
            .collect();
        Self {
            controls,
            curve_base,
            patch_base,
            slots: next,
        }
    }

    /// A patch's net, `u` poles by `v` poles.
    fn dims(&self, topo: &Topology, p: usize) -> (usize, usize) {
        let (su, sv) = topo.patch_strips[p];
        (self.controls[su], self.controls[sv])
    }

    /// A curve's poles.
    fn along(&self, topo: &Topology, c: usize) -> usize {
        self.controls[topo.curve_strip[c]]
    }

    /// The slot of a curve's pole `q`.
    fn curve_slot(&self, topo: &Topology, c: usize, q: usize) -> usize {
        let n = self.along(topo, c);
        let (a, b) = topo.curve_ends[c];
        if q == 0 {
            a
        } else if q == n - 1 {
            b
        } else {
            self.curve_base[c] + q - 1
        }
    }

    /// The slot of pole `(i, j)` of patch `p`.
    fn slot(
        &self,
        topo: &Topology,
        patches: &[PatchDef],
        p: usize,
        (i, j): (usize, usize),
    ) -> usize {
        let (nu, nv) = self.dims(topo, p);
        let corners = topo.corner_vertex[p];
        let (k, w) = match (i, j) {
            (0, 0) => return corners[0],
            (x, 0) if x == nu - 1 => return corners[1],
            (x, y) if x == nu - 1 && y == nv - 1 => return corners[2],
            (0, y) if y == nv - 1 => return corners[3],
            (_, 0) => (0, i),
            (x, _) if x == nu - 1 => (1, j),
            (_, y) if y == nv - 1 => (2, nu - 1 - i),
            (0, _) => (3, nv - 1 - j),
            _ => return self.patch_base[p] + (i - 1) * (nv - 2) + (j - 1),
        };
        match patches[p].sides[k] {
            None => corners[k],
            Some(s) => {
                let n = if k % 2 == 0 { nu } else { nv };
                self.curve_slot(topo, s.curve, if s.reversed { n - 1 - w } else { w })
            }
        }
    }
}

/// The pole of side `k` of an `nu` by `nv` net at `w` poles along its
/// walk, and the pole one row in from it.
fn side_pole(k: usize, w: usize, (nu, nv): (usize, usize)) -> ((usize, usize), (usize, usize)) {
    match k {
        0 => ((w, 0), (w, 1)),
        1 => ((nu - 1, w), (nu - 2, w)),
        2 => ((nu - 1 - w, nv - 1), (nu - 1 - w, nv - 2)),
        _ => ((0, nv - 1 - w), (1, nv - 1 - w)),
    }
}

/// The chart point of a patch's corner `k`.
fn corner_chart(k: usize) -> (f64, f64) {
    match k % 4 {
        0 => (0.0, 0.0),
        1 => (1.0, 0.0),
        2 => (1.0, 1.0),
        _ => (0.0, 1.0),
    }
}

/// The seed patches the points are parameterised over.
struct Seeds {
    surfaces: Vec<SurfaceGeometry>,
    searches: Vec<SurfaceSeeds>,
    /// The box of each seed's net, which holds the patch.
    boxes: Vec<(Point, Point)>,
}

/// The points each patch takes and where on it.
struct Taken {
    patch: Vec<usize>,
    point: Vec<Point>,
    at: Vec<(f64, f64)>,
    /// Each patch's points, by index.
    of_patch: Vec<Vec<usize>>,
    outside: usize,
}

impl Seeds {
    fn new(
        topo: &Topology,
        patches: &[PatchDef],
        drawn: &[Vec<Point>],
        corners: &[Point],
        options: &NetworkOptions,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        // Each curve's chord parameters over the unit interval.
        let chords: Vec<Vec<f64>> = drawn
            .iter()
            .map(|d| {
                let mut at = vec![0.0];
                for w in d.windows(2) {
                    at.push(at[at.len() - 1] + w[0].distance(w[1]));
                }
                let total = at[at.len() - 1];
                at.iter()
                    .map(|t| if total > 0.0 { t / total } else { 0.0 })
                    .collect()
            })
            .collect();
        let along = |c: usize, t: f64| -> Point {
            let (d, at) = (&drawn[c], &chords[c]);
            let k = at.partition_point(|x| *x <= t).clamp(1, d.len() - 1);
            let (a, b) = (at[k - 1], at[k]);
            let f = if b > a {
                ((t - a) / (b - a)).clamp(0.0, 1.0)
            } else {
                0.0
            };
            d[k - 1] + (d[k] - d[k - 1]) * f
        };
        let poles_of = |c: usize, knots: &KnotVector| -> Vec<Point> {
            let mut poles: Vec<Point> = greville(knots).iter().map(|&t| along(c, t)).collect();
            let (a, b) = topo.curve_ends[c];
            let n = poles.len();
            poles[0] = corners[a];
            poles[n - 1] = corners[b];
            poles
        };
        // Each strip's spans: doubled until its curves pass within the
        // tolerance of their points.
        let mut strip_knots = Vec::with_capacity(topo.strips);
        for s in 0..topo.strips {
            let mut spans = options.spans;
            loop {
                let knots = even_knots(spans)?;
                let mut within = true;
                for c in (0..drawn.len()).filter(|&c| topo.curve_strip[c] == s) {
                    let curve =
                        ogeom_geom::BSplineCurve::new(knots.clone(), poles_of(c, &knots), tol)?;
                    for (p, &t) in drawn[c].iter().zip(&chords[c]) {
                        use ogeom_geom::Curve3d as _;
                        if curve.point_at(t, tol)?.distance(*p) > options.tolerance {
                            within = false;
                        }
                    }
                }
                if within || spans >= options.most {
                    strip_knots.push(knots);
                    break;
                }
                spans = (2 * spans).min(options.most);
            }
        }
        let mut surfaces = Vec::with_capacity(patches.len());
        let mut searches = Vec::with_capacity(patches.len());
        let mut boxes = Vec::with_capacity(patches.len());
        for (p, patch) in patches.iter().enumerate() {
            let (su, sv) = topo.patch_strips[p];
            let (ku, kv) = (strip_knots[su].clone(), strip_knots[sv].clone());
            let (nu, nv) = (ku.control_point_count(), kv.control_point_count());
            let mut net = vec![Point::ORIGIN; nu * nv];
            for (k, side) in patch.sides.iter().enumerate() {
                let n = if k % 2 == 0 { nu } else { nv };
                let poles = match side {
                    Some(s) => {
                        let mut poles = poles_of(s.curve, &strip_knots[topo.curve_strip[s.curve]]);
                        if s.reversed {
                            poles.reverse();
                        }
                        poles
                    }
                    None => vec![corners[topo.corner_vertex[p][k]]; n],
                };
                for (w, pole) in poles.into_iter().enumerate() {
                    let ((i, j), _) = side_pole(k, w, (nu, nv));
                    net[i * nv + j] = pole;
                }
            }
            // The bilinearly blended Coons patch of the four sides: on a
            // net whose inner poles stand at the Greville abscissae's blend,
            // linear precision makes it the Coons surface exactly.
            let (xs, ys) = (greville(&ku), greville(&kv));
            let at = |i: usize, j: usize| net[i * nv + j].to_vector();
            let (c00, c10, c11, c01) = (at(0, 0), at(nu - 1, 0), at(nu - 1, nv - 1), at(0, nv - 1));
            let mut filled = net.clone();
            for i in 1..nu - 1 {
                for j in 1..nv - 1 {
                    let (x, y) = (xs[i], ys[j]);
                    let ruled = at(i, 0) * (1.0 - y)
                        + at(i, nv - 1) * y
                        + at(0, j) * (1.0 - x)
                        + at(nu - 1, j) * x;
                    let bilinear = c00 * ((1.0 - x) * (1.0 - y))
                        + c10 * (x * (1.0 - y))
                        + c11 * (x * y)
                        + c01 * ((1.0 - x) * y);
                    filled[i * nv + j] = Point::ORIGIN + (ruled - bilinear);
                }
            }
            boxes.push(net_box(&filled));
            let grid = ControlGrid::new(filled, nu, nv)?;
            let surface = SurfaceGeometry::BSpline(BSplineSurface::new(ku, kv, &grid, tol)?);
            // A B-spline patch seeds its projections by its span boxes; the
            // sample count is for surfaces without them.
            searches.push(SurfaceSeeds::over(&surface, 0, tol)?);
            surfaces.push(surface);
        }
        Ok(Self {
            surfaces,
            searches,
            boxes,
        })
    }

    /// Give each point to the seed nearest it, at the parameters of its
    /// nearest point there; a point whose nearest point on a seed lies on a
    /// side no other patch shares is past that seed.
    fn take(
        &self,
        topo: &Topology,
        patches: &[PatchDef],
        cloud: &[Point],
        tol: Tolerances,
    ) -> OgeomResult<Taken> {
        let border: Vec<[bool; 4]> = patches
            .iter()
            .map(|patch| {
                patch
                    .sides
                    .map(|s| s.is_some_and(|s| topo.users[s.curve].len() < 2))
            })
            .collect();
        let past = |p: usize, (u, v): (f64, f64)| {
            let b = border[p];
            (b[0] && v <= 0.0) || (b[1] && u >= 1.0) || (b[2] && v >= 1.0) || (b[3] && u <= 0.0)
        };
        let starts: Vec<usize> = (0..cloud.len()).step_by(POINT_BATCH).collect();
        let found = ogeom_core::parallel::map_ordered(&starts, |_, &start| {
            let mut out = Vec::new();
            let mut order: Vec<(f64, usize)> = Vec::with_capacity(self.surfaces.len());
            for &q in &cloud[start..(start + POINT_BATCH).min(cloud.len())] {
                order.clear();
                order.extend(
                    self.boxes
                        .iter()
                        .enumerate()
                        .map(|(p, b)| (box_distance(*b, q), p)),
                );
                order.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
                let mut best: Option<(usize, (f64, f64), f64)> = None;
                for &(reach, p) in &order {
                    if best.is_some_and(|b| reach > b.2) {
                        break;
                    }
                    let foot = self.searches[p].project(&self.surfaces[p], q, tol)?;
                    if best.is_none_or(|b| foot.distance < b.2) {
                        best = Some((p, foot.parameters, foot.distance));
                    }
                }
                out.push(best.filter(|b| !past(b.0, b.1)).map(|(p, at, _)| (p, at)));
            }
            Ok::<_, ogeom_core::OgeomError>(out)
        });
        let mut taken = Taken {
            patch: Vec::new(),
            point: Vec::new(),
            at: Vec::new(),
            of_patch: vec![Vec::new(); patches.len()],
            outside: 0,
        };
        let mut i = 0;
        for batch in found {
            for found in batch? {
                match found {
                    Some((p, at)) => {
                        taken.of_patch[p].push(taken.patch.len());
                        taken.patch.push(p);
                        taken.point.push(cloud[i]);
                        taken.at.push(at);
                    }
                    None => taken.outside += 1,
                }
                i += 1;
            }
        }
        Ok(taken)
    }
}

/// The box of a net's poles.
fn net_box(net: &[Point]) -> (Point, Point) {
    let mut lo = net[0];
    let mut hi = net[0];
    for p in net {
        lo = Point::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
        hi = Point::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
    }
    (lo, hi)
}

/// The distance from a point to a box, nought inside it.
fn box_distance((lo, hi): (Point, Point), p: Point) -> f64 {
    let gap = |x: f64, a: f64, b: f64| (a - x).max(0.0).max(x - b);
    let (dx, dy, dz) = (
        gap(p.x, lo.x, hi.x),
        gap(p.y, lo.y, hi.y),
        gap(p.z, lo.z, hi.z),
    );
    dz.mul_add(dz, dx.mul_add(dx, dy * dy)).sqrt()
}

/// The integrals over a knot vector's domain of the products of its basis
/// functions' derivatives, for orders nought to two: `grams[k][i][d]` is
/// the integral of `N_i^(k) N_(i+d)^(k)`, for `d` up to the degree.
fn derivative_grams(knots: &KnotVector) -> [Vec<Vec<f64>>; 3] {
    let p = knots.degree();
    let n = knots.control_point_count();
    let mut grams: [Vec<Vec<f64>>; 3] = core::array::from_fn(|_| vec![vec![0.0; p + 1]; n]);
    let distinct = knots.distinct();
    for w in distinct.windows(2) {
        for (t, weight) in ogeom_math::gauss_legendre_rule(w[0].0, w[1].0) {
            let span = knots.span_unchecked(t);
            let rows = knots.basis_derivatives(span, t, 2);
            for (k, gram) in grams.iter_mut().enumerate() {
                let Some(row) = rows.get(k) else { continue };
                for i in 0..=p {
                    for j in i..=p {
                        gram[span - p + i][j - i] += weight * row[i] * row[j];
                    }
                }
            }
        }
    }
    grams
}

/// One patch's own normal equations over its net, `u` outer: the lower
/// band of its points' Gram matrix and thin-plate energy, and the
/// right-hand sides.
struct Local {
    nv: usize,
    band: usize,
    lower: Vec<Vec<f64>>,
    rhs: Vec<[f64; 3]>,
}

/// The fitting problem, fixed across rounds.
struct Fit<'a> {
    topo: &'a Topology,
    patches: &'a [PatchDef],
    taken: &'a Taken,
    smoothing: f64,
    tol: Tolerances,
}

/// A fit over given spans, each point at its nearest point on it.
struct Round {
    layout: Layout,
    knots: Vec<KnotVector>,
    /// Each slot's pole.
    poles: Vec<Point>,
    surfaces: Vec<BSplineSurface>,
    /// How many unknowns the constraints leave free.
    free: usize,
    feet: Vec<(f64, f64)>,
    distances: Vec<f64>,
}

impl Fit<'_> {
    /// Fit over `spans` per strip.
    fn round(&self, spans: &[usize]) -> OgeomResult<Round> {
        let topo = self.topo;
        let layout = Layout::new(topo, spans);
        let knots: Vec<KnotVector> = spans
            .iter()
            .map(|&s| even_knots(s))
            .collect::<OgeomResult<_>>()?;
        let grams: Vec<[Vec<Vec<f64>>; 3]> = knots.iter().map(derivative_grams).collect();
        let locals: Vec<Local> = ogeom_core::parallel::map_ordered(
            &(0..self.patches.len()).collect::<Vec<_>>(),
            |_, &p| self.local(p, &knots, &grams),
        );

        // Without the tangency, for the ratios across each curve.
        let loose = self.solve(&layout, &locals, &Basis::identity(layout.slots))?;
        let rows = self.tangency(&layout, &loose)?;
        let basis = null_space(layout.slots, &rows);
        let poles = self.solve(&layout, &locals, &basis)?;

        let mut surfaces = Vec::with_capacity(self.patches.len());
        for p in 0..self.patches.len() {
            let (nu, nv) = layout.dims(topo, p);
            let mut net = Vec::with_capacity(nu * nv);
            for i in 0..nu {
                for j in 0..nv {
                    net.push(poles[layout.slot(topo, self.patches, p, (i, j))]);
                }
            }
            let (su, sv) = topo.patch_strips[p];
            let grid = ControlGrid::new(net, nu, nv)?;
            surfaces.push(BSplineSurface::new(
                knots[su].clone(),
                knots[sv].clone(),
                &grid,
                self.tol,
            )?);
        }
        let geometry: Vec<SurfaceGeometry> = surfaces
            .iter()
            .map(|s| SurfaceGeometry::BSpline(s.clone()))
            .collect();
        let taken = self.taken;
        let starts: Vec<usize> = (0..taken.point.len()).step_by(POINT_BATCH).collect();
        let found = ogeom_core::parallel::map_ordered(&starts, |_, &start| {
            (start..(start + POINT_BATCH).min(taken.point.len()))
                .map(|i| {
                    project_on_surface_from(
                        &geometry[taken.patch[i]],
                        taken.point[i],
                        taken.at[i],
                        self.tol,
                    )
                })
                .collect::<OgeomResult<Vec<_>>>()
        });
        let mut feet = Vec::with_capacity(taken.point.len());
        let mut distances = Vec::with_capacity(taken.point.len());
        for batch in found {
            for foot in batch? {
                feet.push(foot.parameters);
                distances.push(foot.distance);
            }
        }
        Ok(Round {
            layout,
            knots,
            poles,
            surfaces,
            free: basis.free,
            feet,
            distances,
        })
    }

    /// Patch `p`'s normal equations over its own net.
    fn local(&self, p: usize, knots: &[KnotVector], grams: &[[Vec<Vec<f64>>; 3]]) -> Local {
        let (su, sv) = self.topo.patch_strips[p];
        let (ku, kv) = (&knots[su], &knots[sv]);
        let (nu, nv) = (ku.control_point_count(), kv.control_point_count());
        let band = DEGREE * nv + DEGREE;
        let mut lower = vec![vec![0.0; band + 1]; nu * nv];
        let mut rhs = vec![[0.0; 3]; nu * nv];
        let mine = &self.taken.of_patch[p];
        if !mine.is_empty() {
            let share = 1.0 / precise(mine.len());
            let mut terms: Vec<(usize, f64)> = Vec::with_capacity((DEGREE + 1) * (DEGREE + 1));
            for &i in mine {
                let (u, v) = self.taken.at[i];
                let q = self.taken.point[i];
                let (us, vs) = (ku.span_unchecked(u), kv.span_unchecked(v));
                let (bu, bv) = (ku.basis(us, u), kv.basis(vs, v));
                terms.clear();
                for (a, x) in bu.iter().enumerate() {
                    for (b, y) in bv.iter().enumerate() {
                        terms.push(((us - DEGREE + a) * nv + vs - DEGREE + b, x * y));
                    }
                }
                for (k, &(r, x)) in terms.iter().enumerate() {
                    rhs[r][0] += share * x * q.x;
                    rhs[r][1] += share * x * q.y;
                    rhs[r][2] += share * x * q.z;
                    for &(c, y) in &terms[..=k] {
                        let (r, c) = if r >= c { (r, c) } else { (c, r) };
                        lower[r][r - c] += share * x * y;
                    }
                }
            }
        }
        if self.smoothing > 0.0 {
            let (gu, gv) = (&grams[su], &grams[sv]);
            for a in 0..nu {
                for b in 0..nv {
                    let r = a * nv + b;
                    for c in a.saturating_sub(DEGREE)..=a {
                        for d in b.saturating_sub(DEGREE)..(b + DEGREE + 1).min(nv) {
                            let col = c * nv + d;
                            if col > r {
                                continue;
                            }
                            let (f, x) = (a.min(c), a.abs_diff(c));
                            let (e, y) = (b.min(d), b.abs_diff(d));
                            let energy = gu[2][f][x] * gv[0][e][y]
                                + 2.0 * gu[1][f][x] * gv[1][e][y]
                                + gu[0][f][x] * gv[2][e][y];
                            lower[r][r - col] += self.smoothing * energy;
                        }
                    }
                }
            }
        }
        Local {
            nv,
            band,
            lower,
            rhs,
        }
    }

    /// Solve the patches' equations together over the free unknowns of
    /// `basis`; every slot's pole.
    fn solve(&self, layout: &Layout, locals: &[Local], basis: &Basis) -> OgeomResult<Vec<Point>> {
        let mut normal = Lower::default();
        let mut rhs = vec![vec![0.0; basis.free]; 3];
        for (p, local) in locals.iter().enumerate() {
            let slot_of =
                |r: usize| layout.slot(self.topo, self.patches, p, (r / local.nv, r % local.nv));
            for (r, row) in local.lower.iter().enumerate() {
                let sr = &basis.of_slot[slot_of(r)];
                for (y, z) in sr {
                    for (side, value) in rhs.iter_mut().zip(local.rhs[r]) {
                        side[*y] += z * value;
                    }
                }
                for (d, &value) in row.iter().enumerate().take(local.band.min(r) + 1) {
                    if value == 0.0 {
                        continue;
                    }
                    if d == 0 {
                        for &(y1, z1) in sr {
                            for &(y2, z2) in sr {
                                normal.add(y1, y2, z1 * z2 * value);
                            }
                        }
                    } else {
                        let sc = &basis.of_slot[slot_of(r - d)];
                        for &(y1, z1) in sr {
                            for &(y2, z2) in sc {
                                normal.add(y1, y2, z1 * z2 * value);
                                normal.add(y2, y1, z1 * z2 * value);
                            }
                        }
                    }
                }
            }
        }
        let Some(solved) = solve_spd(basis.free, &normal, &rhs) else {
            ogeom_bail!(
                Numeric,
                "the network's system is singular; raise the smoothing so the fairing settles \
                 the poles no point reaches"
            );
        };
        Ok(basis
            .of_slot
            .iter()
            .map(|combo| {
                let mut p = [0.0; 3];
                for &(y, z) in combo {
                    for (axis, value) in p.iter_mut().enumerate() {
                        *value += z * solved[axis][y];
                    }
                }
                Point::new(p[0], p[1], p[2])
            })
            .collect())
    }

    /// The tangency rows: for each pole of each curve held tangent, the
    /// pole and the poles one row into either patch on one line,
    /// `(1 + λ) B − L − λ R = 0`, `L` in the curve's first patch, `λ` its
    /// chain's ratio of pole spacings in the fit `loose`.
    fn tangency(&self, layout: &Layout, loose: &[Point]) -> OgeomResult<Vec<Vec<(usize, f64)>>> {
        let topo = self.topo;
        // Each held triple: (curve, pole, B, L, R) slots.
        let mut triples: Vec<(usize, usize, usize, usize)> = Vec::new();
        for c in 0..topo.users.len() {
            if !topo.held[c] {
                continue;
            }
            let n = layout.along(topo, c);
            let beside = |(p, k): (usize, usize), q: usize| {
                let reversed = self.patches[p].sides[k].is_some_and(|s| s.reversed);
                let w = if reversed { n - 1 - q } else { q };
                let (_, inner) = side_pole(k, w, layout.dims(topo, p));
                layout.slot(topo, self.patches, p, inner)
            };
            for q in 0..n {
                triples.push((
                    c,
                    layout.curve_slot(topo, c, q),
                    beside(topo.users[c][0], q),
                    beside(topo.users[c][1], q),
                ));
            }
        }
        let mut sums = vec![(0.0, 0.0); topo.chains];
        for &(c, b, l, r) in &triples {
            let Some((chain, inverse)) = topo.chain[c] else {
                continue;
            };
            let (near, far) = (loose[b].distance(loose[l]), loose[r].distance(loose[b]));
            let (num, den) = if inverse { (far, near) } else { (near, far) };
            sums[chain].0 += num;
            sums[chain].1 += den;
        }
        let mut ratio = Vec::with_capacity(topo.chains);
        for (chain, &(num, den)) in sums.iter().enumerate() {
            if topo.conflicted[chain] {
                ratio.push(1.0);
                continue;
            }
            let rho = num / den;
            if !(rho > 0.0 && rho.is_finite()) {
                ogeom_bail!(
                    Numeric,
                    "the poles either side of a tangent curve collapse onto it; raise the \
                     smoothing"
                );
            }
            ratio.push(rho);
        }
        Ok(triples
            .iter()
            .filter_map(|&(c, b, l, r)| {
                let (chain, inverse) = topo.chain[c]?;
                let lambda = if inverse {
                    1.0 / ratio[chain]
                } else {
                    ratio[chain]
                };
                Some(vec![(b, 1.0 + lambda), (l, -1.0), (r, -lambda)])
            })
            .collect())
    }

    /// The worst miss of each patch among its cells holding at least
    /// `judged` points, that is not a stray, with the patch; worst first.
    /// Then whether any point misses without being a stray, and whether
    /// any is a stray.
    ///
    /// A stray stands farther from the fit than noise as wide as the rest
    /// of its cell (their root mean square distance) reaches over the
    /// whole cloud, `√(2 ln m)` deviations.
    fn missed(
        &self,
        round: &Round,
        tolerance: f64,
        judged: usize,
    ) -> (Vec<(f64, usize)>, bool, bool) {
        let topo = self.topo;
        let reach = (2.0 * precise(round.distances.len().max(1)).ln()).sqrt();
        let (mut any, mut strays) = (false, false);
        let mut missed = Vec::new();
        for (p, mine) in self.taken.of_patch.iter().enumerate() {
            let (su, sv) = topo.patch_strips[p];
            let (ku, kv) = (&round.knots[su], &round.knots[sv]);
            let (cu, cv) = (
                ku.control_point_count() - DEGREE,
                kv.control_point_count() - DEGREE,
            );
            let cell_of = |i: usize| {
                let (u, v) = round.feet[i];
                (ku.span_unchecked(u) - DEGREE) * cv + kv.span_unchecked(v) - DEGREE
            };
            let mut count = vec![0usize; cu * cv];
            let mut square = vec![0.0; cu * cv];
            for &i in mine {
                let k = cell_of(i);
                count[k] += 1;
                square[k] += round.distances[i] * round.distances[i];
            }
            let mut worst: Option<f64> = None;
            for &i in mine {
                let d = round.distances[i];
                if d <= tolerance {
                    continue;
                }
                let k = cell_of(i);
                let rest = count[k] - 1;
                let stray =
                    rest > 0 && d > reach * ((square[k] - d * d).max(0.0) / precise(rest)).sqrt();
                if stray {
                    strays = true;
                    continue;
                }
                any = true;
                if count[k] >= judged && worst.is_none_or(|w| d > w) {
                    worst = Some(d);
                }
            }
            if let Some(d) = worst {
                missed.push((d, p));
            }
        }
        missed.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        (missed, any, strays)
    }
}
