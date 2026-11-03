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
//! - **a corner of an odd number of patches is held tangent where it
//!   can be.** The triples hold the parameter lines that cross a curve to
//!   one direction on both sides of it. At a corner of three patches (or
//!   five, or the point a three-sided patch's missing side shrinks to) the
//!   lines crossing a curve there lean one way on one side and the other
//!   way on the other, and no triple holds; tangency there needs the cross
//!   derivatives to differ by a multiple of the curve's own tangent,
//!   `D₁ + λ D₂ + γ T = 0`. With `λ` constant and `γ` linear along the
//!   curve that is a cubic spline on the curve's knots, nought where it is
//!   nought at as many points as it has coefficients: an exact linear
//!   condition, though not a pole-by-pole one. At the corner it puts the
//!   curves' first poles in one plane, at ratios read off the fit without
//!   the constraints. A corner of three or five patches whose curves all
//!   run on to the network's border is held so, and every curve into it is
//!   tangent to rounding, as long as the fit still comes within the
//!   tolerance: the one ratio along each curve the condition allows is a
//!   restriction the points at their seed's parameters need not admit.
//!   Where they do not, and at any other corner of other than four patches
//!   (an even count, whose twists close in a cycle that can be singular, or
//!   a curve running on to a corner of four, which closes one there), a
//!   curve with such a corner at an end is fitted without the triples:
//!   near tangent as the points are smooth, its angle measured into
//!   [`Network::worst_tangent`]. Every such corner's widest angle between
//!   the patches' normals is reported in [`Network::corners`].
//! - **a side may be a chain of curves.** A curve may end on the middle of
//!   a patch's side (a T-junction), the side then a chain of curves meeting
//!   there ([`ChainedPatch`]). The patch takes a knot of multiplicity three
//!   at the junction, so its boundary row is the chain's poles, shared
//!   exactly, and the side opposite crosses the same junction: it is
//!   split alike, or is one curve whose knots carry the junction too. The
//!   patch is held `C¹` across the parameter line between the two, pole by
//!   pole, at the ratio of the knot spacings either side. That line runs
//!   straight on from the curve ending at the junction, so that curve is
//!   held tangent at that same ratio, and the chain's curves either side
//!   of the junction share theirs, as at a corner of four patches.
//! - **spans are raised per strip.** The curves opposite each other in a
//!   patch share knots, and so do all the curves of a strip of patches
//!   joined through their opposite sides. Where a patch misses the points
//!   by more than the tolerance its strips gain spans, the way
//!   [`fit_surface_scattered_over`](ogeom_geom::fit::fit_surface_scattered_over)
//!   refines a single patch.
//! - **points may be reparameterised.** Once refinement ends,
//!   [`fit_network_chained`] can move every point's parameters to its
//!   nearest point on its patch and fit the patches again, a fit kept only
//!   where it comes nearer the points.
//!
//! The constraints are solved for exactly, by elimination, before the
//! least squares: what is held tangent is tangent to rounding.

mod build;
mod solve;
mod topology;

use std::collections::BTreeMap;

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

/// A patch of the network whose sides may each be a chain of curves.
///
/// The sides are walked as a [`PatchDef`]'s, each side's curves in walking
/// order, each meeting the next at a junction, where a curve of a
/// neighbouring patch may end (a T-junction). An empty side is missing, as
/// a [`PatchDef`]'s `None`.
///
/// A side split into several curves has its patch cross the same
/// junctions from the side opposite: that side is split into as many
/// curves, the `i`-th along one direction opposite the `i`-th along it on
/// the other, or is one curve, which then carries the junctions to the
/// patch beyond it. A curve of a split side is not itself crossed by
/// junctions elsewhere.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainedPatch {
    /// The sides, in walking order, each its curves in walking order.
    pub sides: [Vec<PatchSide>; 4],
}

impl ChainedPatch {
    /// A four-sided patch from its sides' curves in walking order, each
    /// with whether it runs against the walk.
    #[must_use]
    pub fn four(sides: [Vec<(usize, bool)>; 4]) -> Self {
        Self {
            sides: sides.map(|side| {
                side.into_iter()
                    .map(|(curve, reversed)| PatchSide { curve, reversed })
                    .collect()
            }),
        }
    }
}

impl From<PatchDef> for ChainedPatch {
    fn from(patch: PatchDef) -> Self {
        Self {
            sides: patch.sides.map(|side| side.into_iter().collect()),
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
/// curve at it smooth (a patch whose side runs straight through it counts
/// twice): held tangent to rounding, or near tangent, as the [module
/// documentation](self) describes.
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
/// rather than a miss. Where corners of an odd number of patches are held
/// tangent and the refined fit still misses a point, the network is
/// refined again from the start with those corners near tangent instead.
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
    let chained: Vec<ChainedPatch> = patches.iter().copied().map(ChainedPatch::from).collect();
    fit_network_chained(model, &chained, curves, points, options, 0, tol)
}

/// [`fit_network`] over patches whose sides may be chains of curves, a
/// curve ending on the middle of another's side at a T-junction, with
/// parameter correction.
///
/// Every curve of a chained side is shared pole for pole as a whole side
/// is, and a T-junction is held tangent as the [module
/// documentation](self) describes; a curve that runs between two
/// junctions whose ratios differ is left near tangent, its angle measured
/// into [`Network::worst_tangent`].
///
/// Once refinement ends, every point's parameters move to its nearest
/// point on its patch and the patches are fitted again at the same spans,
/// up to `corrections` times. A fit is kept only where its root mean
/// square distance to the points is below the one before, and, where every
/// point was within the tolerance, every point stays within it; the first
/// that is not ends the corrections. [`fit_network`] corrects none. A
/// correction fits the points, not the shape they sample: where the points
/// leave a pole unsettled (the point a three-sided patch shrinks its
/// missing side to, which only the points round it place), it can move that
/// pole as freely as the seed's parameters did.
///
/// # Errors
///
/// As [`fit_network`], and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if two
/// opposite sides of a patch are split into different numbers of curves,
/// a curve of a split side is itself crossed by junctions elsewhere, or the
/// curves of a split side meet one another, or the patch's corners, other
/// than at their junctions.
pub fn fit_network_chained(
    model: &mut Model,
    patches: &[ChainedPatch],
    curves: &[CurveDef],
    points: &PointIndex,
    options: &NetworkOptions,
    corrections: usize,
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
    let weights = strip_weights(&topo, &drawn);
    let seeds = Seeds::new(&topo, patches, &drawn, &corners, &weights, options, tol)?;
    let cloud = points.points();
    let taken = seeds.take(&topo, patches, &cloud, tol)?;

    let mut fit = Fit {
        topo: &topo,
        patches,
        taken: &taken,
        weights: &weights,
        smoothing: options.smoothing,
        corners: !topo.odd_corners.is_empty(),
        tol,
    };
    let (mut spans, mut current, mut ended) = fit.refine(options)?;
    // A corner held exactly tangent binds its curves' cross derivatives to
    // a relation with one ratio along each, which the points at the seed's
    // parameters need not admit within the tolerance at any spans. Where
    // they do not, the corners are left near tangent, as without them.
    if fit.corners && ended != Ended::Met {
        fit.corners = false;
        (spans, current, ended) = fit.refine(options)?;
    }
    // Parameter correction, at the spans refinement settled on. Refinement
    // fits every try at the seed's parameters, so that tries compare as
    // fits of the same points at the same parameters rather than as
    // corrections that have run further or less far.
    for _ in 0..corrections {
        let next = fit.round(&spans, &current.feet)?;
        let kept_within = ended != Ended::Met || next.distances.iter().all(|d| *d <= tolerance);
        if next.square() < current.square() && kept_within {
            current = next;
        } else {
            break;
        }
    }

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
        let [su, sv] = &topo.patch_strips[p];
        fit.spans = (
            su.iter().map(|&s| spans[s]).sum(),
            sv.iter().map(|&s| spans[s]).sum(),
        );
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

/// Each strip's weight: the mean length of the curves that cross it alone,
/// which sets its share of a run of strips joined at junctions. One for a
/// strip no such curve crosses, which no junction then bounds.
fn strip_weights(topo: &Topology, drawn: &[Vec<Point>]) -> Vec<f64> {
    let mut sums = vec![(0.0, 0usize); topo.strips];
    for (c, strips) in topo.curve_strip.iter().enumerate() {
        if let [s] = strips.as_slice() {
            sums[*s].0 += drawn[c]
                .windows(2)
                .map(|w| w[0].distance(w[1]))
                .sum::<f64>();
            sums[*s].1 += 1;
        }
    }
    sums.iter()
        .map(|&(sum, n)| if n > 0 { sum / precise(n) } else { 1.0 })
        .collect()
}

/// The knots of a run of strips over the unit interval, and the bounds of
/// each strip's share of it: each strip's spans even over its share, the
/// shares in proportion to the strips' weights, and a knot of multiplicity
/// the degree at each junction between two strips, where the run's poles
/// pass through the junction's.
fn run_knots(
    strips: &[usize],
    spans: &[usize],
    weights: &[f64],
) -> OgeomResult<(KnotVector, Vec<f64>)> {
    if let [one] = strips {
        return Ok((even_knots(spans[*one])?, vec![0.0, 1.0]));
    }
    let total: f64 = strips.iter().map(|&s| weights[s]).sum();
    let mut bounds = vec![0.0];
    let mut run = 0.0;
    for (i, &s) in strips.iter().enumerate() {
        run += weights[s];
        bounds.push(if i + 1 == strips.len() {
            1.0
        } else {
            run / total
        });
    }
    let mut knots = vec![0.0; DEGREE + 1];
    for (i, &s) in strips.iter().enumerate() {
        let (a, b) = (bounds[i], bounds[i + 1]);
        for k in 1..spans[s] {
            knots.push((b - a).mul_add(precise(k) / precise(spans[s]), a));
        }
        if i + 1 < strips.len() {
            knots.extend(core::iter::repeat_n(b, DEGREE));
        }
    }
    knots.extend(core::iter::repeat_n(1.0, DEGREE + 1));
    Ok((KnotVector::new(knots, DEGREE)?, bounds))
}

/// The poles of a run of strips of `controls` poles each, a junction's
/// pole shared by the strips either side.
fn run_poles(controls: &[usize], strips: &[usize]) -> usize {
    strips.iter().map(|&s| controls[s]).sum::<usize>() + 1 - strips.len()
}

/// Where the poles of the patches' nets stand among the unknowns: a slot
/// per vertex, per inner pole of each curve and per inner pole of each
/// patch.
struct Layout {
    /// The poles along each strip.
    controls: Vec<usize>,
    /// Each curve's poles.
    curve_poles: Vec<usize>,
    /// Each patch's net, `u` poles by `v` poles.
    dims: Vec<(usize, usize)>,
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
        let curve_poles: Vec<usize> = topo
            .curve_strip
            .iter()
            .map(|strips| run_poles(&controls, strips))
            .collect();
        let dims: Vec<(usize, usize)> = topo
            .patch_strips
            .iter()
            .map(|[u, v]| (run_poles(&controls, u), run_poles(&controls, v)))
            .collect();
        let mut next = topo.vertices;
        let curve_base = curve_poles
            .iter()
            .map(|&n| {
                let base = next;
                next += n - 2;
                base
            })
            .collect();
        let patch_base = dims
            .iter()
            .map(|&(nu, nv)| {
                let base = next;
                next += (nu - 2) * (nv - 2);
                base
            })
            .collect();
        Self {
            controls,
            curve_poles,
            dims,
            curve_base,
            patch_base,
            slots: next,
        }
    }

    /// A patch's net, `u` poles by `v` poles.
    fn dims(&self, p: usize) -> (usize, usize) {
        self.dims[p]
    }

    /// A curve's poles.
    fn along(&self, c: usize) -> usize {
        self.curve_poles[c]
    }

    /// The slot of a curve's pole `q`.
    fn curve_slot(&self, topo: &Topology, c: usize, q: usize) -> usize {
        let n = self.along(c);
        let (a, b) = topo.curve_ends[c];
        if q == 0 {
            a
        } else if q == n - 1 {
            b
        } else {
            self.curve_base[c] + q - 1
        }
    }

    /// How many poles along side `k` of patch `p`'s walk curve `c`'s pole
    /// `q` stands.
    fn walk_index(
        &self,
        patches: &[ChainedPatch],
        (p, k): (usize, usize),
        c: usize,
        q: usize,
    ) -> usize {
        let mut offset = 0;
        for piece in &patches[p].sides[k] {
            let n = self.along(piece.curve);
            if piece.curve == c {
                return offset + if piece.reversed { n - 1 - q } else { q };
            }
            offset += n - 1;
        }
        offset
    }

    /// The slot of pole `(i, j)` of patch `p`.
    fn slot(
        &self,
        topo: &Topology,
        patches: &[ChainedPatch],
        p: usize,
        (i, j): (usize, usize),
    ) -> usize {
        let (nu, nv) = self.dims(p);
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
        // The curve of the side the pole lies on, its curves' poles laid
        // end to end along the walk, a junction's pole shared.
        let mut offset = 0;
        for piece in &patches[p].sides[k] {
            let n = self.along(piece.curve);
            if w < offset + n {
                let q = w - offset;
                return self.curve_slot(
                    topo,
                    piece.curve,
                    if piece.reversed { n - 1 - q } else { q },
                );
            }
            offset += n - 1;
        }
        corners[k]
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

/// The share of side `k`'s walk, from its start, that its `i`-th curve
/// covers, `bounds` the bounds of the shares of the strips along the
/// side's direction.
fn walk_share(k: usize, i: usize, count: usize, bounds: &[f64]) -> (f64, f64) {
    if count < 2 {
        return (0.0, 1.0);
    }
    let m = bounds.len() - 1;
    if k < 2 {
        (bounds[i], bounds[i + 1])
    } else {
        (1.0 - bounds[m - i], 1.0 - bounds[m - 1 - i])
    }
}

/// Which curve of side `k`, of `count`, holds the point at `x` along the
/// side's direction, `bounds` the bounds of the shares of its strips.
fn piece_at(k: usize, x: f64, count: usize, bounds: &[f64]) -> usize {
    if count < 2 {
        return 0;
    }
    let m = bounds.len() - 1;
    let strip = bounds[1..m].iter().filter(|&&b| x > b).count();
    if k < 2 { strip } else { m - 1 - strip }
}

/// The seed patches the points are parameterised over.
struct Seeds {
    surfaces: Vec<SurfaceGeometry>,
    searches: Vec<SurfaceSeeds>,
    /// The box of each seed's net, which holds the patch.
    boxes: Vec<(Point, Point)>,
    /// The bounds of each seed's strips' shares along `u` and `v`.
    bounds: Vec<[Vec<f64>; 2]>,
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
    #[allow(clippy::too_many_lines, reason = "the seed curves, then the patches")]
    fn new(
        topo: &Topology,
        patches: &[ChainedPatch],
        drawn: &[Vec<Point>],
        corners: &[Point],
        weights: &[f64],
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
        // Each strip's spans: doubled until every curve crossing it passes
        // within the tolerance of its points.
        let mut spans = vec![options.spans; topo.strips];
        let mut settled = vec![false; drawn.len()];
        loop {
            let mut raise = vec![false; topo.strips];
            for c in 0..drawn.len() {
                if settled[c] {
                    continue;
                }
                let (knots, _) = run_knots(&topo.curve_strip[c], &spans, weights)?;
                let curve = ogeom_geom::BSplineCurve::new(knots.clone(), poles_of(c, &knots), tol)?;
                let mut within = true;
                for (p, &t) in drawn[c].iter().zip(&chords[c]) {
                    use ogeom_geom::Curve3d as _;
                    if curve.point_at(t, tol)?.distance(*p) > options.tolerance {
                        within = false;
                    }
                }
                let raised: Vec<usize> = topo.curve_strip[c]
                    .iter()
                    .copied()
                    .filter(|&s| spans[s] < options.most)
                    .collect();
                if within || raised.is_empty() {
                    settled[c] = true;
                }
                if !within {
                    for s in raised {
                        raise[s] = true;
                    }
                }
            }
            if !raise.contains(&true) {
                break;
            }
            for (s, spans) in spans.iter_mut().enumerate() {
                if raise[s] {
                    *spans = (2 * *spans).min(options.most);
                }
            }
            // A curve crossing a raised strip is drawn again.
            for (c, strips) in topo.curve_strip.iter().enumerate() {
                if strips.iter().any(|&s| raise[s]) {
                    settled[c] = false;
                }
            }
        }
        let mut surfaces = Vec::with_capacity(patches.len());
        let mut searches = Vec::with_capacity(patches.len());
        let mut boxes = Vec::with_capacity(patches.len());
        let mut all_bounds = Vec::with_capacity(patches.len());
        for (p, patch) in patches.iter().enumerate() {
            let [su, sv] = &topo.patch_strips[p];
            let (ku, bu) = run_knots(su, &spans, weights)?;
            let (kv, bv) = run_knots(sv, &spans, weights)?;
            let (nu, nv) = (ku.control_point_count(), kv.control_point_count());
            let mut net = vec![Point::ORIGIN; nu * nv];
            for (k, side) in patch.sides.iter().enumerate() {
                let n = if k % 2 == 0 { nu } else { nv };
                let poles = if side.is_empty() {
                    vec![corners[topo.corner_vertex[p][k]]; n]
                } else {
                    let mut poles: Vec<Point> = Vec::with_capacity(n);
                    for piece in side {
                        let (knots, _) =
                            run_knots(&topo.curve_strip[piece.curve], &spans, weights)?;
                        let mut own = poles_of(piece.curve, &knots);
                        if piece.reversed {
                            own.reverse();
                        }
                        // A junction's pole ends one curve and starts the
                        // next.
                        let skip = usize::from(!poles.is_empty());
                        poles.extend(own.into_iter().skip(skip));
                    }
                    poles
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
            all_bounds.push([bu, bv]);
        }
        Ok(Self {
            surfaces,
            searches,
            boxes,
            bounds: all_bounds,
        })
    }

    /// Give each point to the seed nearest it, at the parameters of its
    /// nearest point there; a point whose nearest point on a seed lies on a
    /// side no other patch shares is past that seed.
    fn take(
        &self,
        topo: &Topology,
        patches: &[ChainedPatch],
        cloud: &[Point],
        tol: Tolerances,
    ) -> OgeomResult<Taken> {
        // Whether each curve of each side is on the network's border.
        let border: Vec<[Vec<bool>; 4]> = patches
            .iter()
            .map(|patch| {
                core::array::from_fn(|k| {
                    patch.sides[k]
                        .iter()
                        .map(|s| topo.users[s.curve].len() < 2)
                        .collect()
                })
            })
            .collect();
        let past = |p: usize, (u, v): (f64, f64)| {
            let free = |k: usize, x: f64| {
                let sides = &border[p][k];
                let i = piece_at(k, x, sides.len(), &self.bounds[p][k % 2]);
                sides.get(i).copied().unwrap_or(false)
            };
            (v <= 0.0 && free(0, u))
                || (u >= 1.0 && free(1, v))
                || (v >= 1.0 && free(2, u))
                || (u <= 0.0 && free(3, v))
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
    patches: &'a [ChainedPatch],
    taken: &'a Taken,
    weights: &'a [f64],
    smoothing: f64,
    /// Whether the corners of an odd number of patches are held exactly
    /// tangent.
    corners: bool,
    tol: Tolerances,
}

/// A fit over given spans, each point at its nearest point on it.
struct Round {
    layout: Layout,
    /// Each patch's knots along `u` and `v`.
    patch_knots: Vec<[KnotVector; 2]>,
    /// The bounds of each patch's strips' shares along `u` and `v`.
    patch_bounds: Vec<[Vec<f64>; 2]>,
    /// Each curve's knots.
    curve_knots: Vec<KnotVector>,
    /// Each slot's pole.
    poles: Vec<Point>,
    surfaces: Vec<BSplineSurface>,
    /// How many unknowns the constraints leave free.
    free: usize,
    feet: Vec<(f64, f64)>,
    distances: Vec<f64>,
}

impl Round {
    /// The sum of the squared distances.
    fn square(&self) -> f64 {
        self.distances.iter().map(|d| d * d).sum()
    }
}

/// Slots and their weights.
type Weighted = Vec<(usize, f64)>;

/// A constraint row keyed by its slots: the middle one, then the other two
/// in order.
type Keyed = BTreeMap<(usize, usize, usize), Vec<(usize, f64)>>;

impl Fit<'_> {
    /// Fit from [`NetworkOptions::spans`] per strip, raising them while
    /// a point misses, as [`fit_network`] describes; the spans, the fit and
    /// why refinement stopped.
    fn refine(&self, options: &NetworkOptions) -> OgeomResult<(Vec<usize>, Round, Ended)> {
        let (topo, tolerance) = (self.topo, options.tolerance);
        let mut spans = vec![options.spans; topo.strips];
        let mut current = self.round(&spans, &self.taken.at)?;
        let ended = loop {
            if current.distances.iter().all(|d| *d <= tolerance) {
                break Ended::Met;
            }
            let (missed, any, strays) =
                self.missed(&current, tolerance, (DEGREE + 1) * (DEGREE + 1));
            if !any && strays {
                break Ended::Strays;
            }
            // The strips of every patch that misses, each tried alone: a
            // patch can miss for a neighbour too coarse to follow the
            // shape, which its tangency binds it to.
            let mut tries: Vec<usize> = missed
                .iter()
                .flat_map(|&(_, p)| topo.patch_strips[p].iter().flatten().copied())
                .filter(|&s| spans[s] < options.most)
                .collect();
            tries.sort_unstable();
            tries.dedup();
            if tries.is_empty() {
                break Ended::Bounds;
            }
            // Mallows' criterion, as the single-patch refinement scores:
            // the squared distances plus twice the noise's variance per
            // free pole, the variance read off the median distance.
            let variance = (median(&current.distances) / NORMAL_QUARTILE).powi(2);
            let score = |r: &Round| 2.0f64.mul_add(variance * precise(r.free), r.square());
            let mut best: Option<(f64, Vec<usize>, Round)> = None;
            for raised in tries {
                let mut next = spans.clone();
                next[raised] = (2 * next[raised]).min(options.most);
                let tried = self.round(&next, &self.taken.at)?;
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
        Ok((spans, current, ended))
    }

    /// Fit over `spans` per strip, the points at `at`.
    fn round(&self, spans: &[usize], at: &[(f64, f64)]) -> OgeomResult<Round> {
        let topo = self.topo;
        let layout = Layout::new(topo, spans);
        let mut patch_knots = Vec::with_capacity(self.patches.len());
        let mut patch_bounds = Vec::with_capacity(self.patches.len());
        for [su, sv] in &topo.patch_strips {
            let (ku, bu) = run_knots(su, spans, self.weights)?;
            let (kv, bv) = run_knots(sv, spans, self.weights)?;
            patch_knots.push([ku, kv]);
            patch_bounds.push([bu, bv]);
        }
        let curve_knots: Vec<KnotVector> = topo
            .curve_strip
            .iter()
            .map(|strips| run_knots(strips, spans, self.weights).map(|k| k.0))
            .collect::<OgeomResult<_>>()?;
        let locals: Vec<Local> = ogeom_core::parallel::map_ordered(
            &(0..self.patches.len()).collect::<Vec<_>>(),
            |_, &p| self.local(p, at, &patch_knots[p]),
        );

        // Without the tangency, for the ratios across each curve; the
        // patches held smooth across their junctions all the same.
        let junctions = self.junction_rows(&layout, spans);
        let loose_basis = if junctions.is_empty() {
            Basis::identity(layout.slots)
        } else {
            null_space(
                layout.slots,
                &junctions.values().cloned().collect::<Vec<_>>(),
            )
        };
        let loose = self.solve(&layout, &locals, &loose_basis)?;
        let mut rows: Vec<Vec<(usize, f64)>> = junctions.values().cloned().collect();
        rows.extend(self.tangency(&layout, &loose, &junctions)?);
        if self.corners {
            rows.extend(self.corner_rows(&layout, &loose, &curve_knots));
        }
        let basis = null_space(layout.slots, &rows);
        let poles = self.solve(&layout, &locals, &basis)?;

        let mut surfaces = Vec::with_capacity(self.patches.len());
        for (p, [ku, kv]) in patch_knots.iter().enumerate() {
            let (nu, nv) = layout.dims(p);
            let mut net = Vec::with_capacity(nu * nv);
            for i in 0..nu {
                for j in 0..nv {
                    net.push(poles[layout.slot(topo, self.patches, p, (i, j))]);
                }
            }
            let grid = ControlGrid::new(net, nu, nv)?;
            surfaces.push(BSplineSurface::new(
                ku.clone(),
                kv.clone(),
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
                        at[i],
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
            patch_knots,
            patch_bounds,
            curve_knots,
            poles,
            surfaces,
            free: basis.free,
            feet,
            distances,
        })
    }

    /// The rows holding each patch `C¹` across the parameter line through
    /// the junctions of its split sides, along every row of poles across
    /// it: the pole at the junction's knot and its neighbours either side
    /// on one line, `(h₁ + h₂) J − h₂ P − h₁ Q = 0`, `h₁` and `h₂` the knot
    /// spacings either side. A row two patches hold alike (along the curve
    /// carrying the junction between them) is kept once: the spacings are
    /// each strip's weight over its spans, the same numbers from either
    /// patch.
    fn junction_rows(&self, layout: &Layout, spans: &[usize]) -> Keyed {
        let topo = self.topo;
        let spacing = |s: usize| self.weights[s] / precise(spans[s]);
        let mut rows = Keyed::new();
        for p in 0..self.patches.len() {
            let (nu, nv) = layout.dims(p);
            for (dir, strips) in topo.patch_strips[p].iter().enumerate() {
                let across = if dir == 0 { nv } else { nu };
                let mut at = 0;
                for pair in strips.windows(2) {
                    at += layout.controls[pair[0]] - 1;
                    let (h1, h2) = (spacing(pair[0]), spacing(pair[1]));
                    for j in 0..across {
                        let slot = |i: usize| {
                            let pole = if dir == 0 { (i, j) } else { (j, i) };
                            layout.slot(topo, self.patches, p, pole)
                        };
                        let (mid, before, after) = (slot(at), slot(at - 1), slot(at + 1));
                        // A row shrunk to a missing side's point holds
                        // nothing.
                        if mid == before || mid == after || before == after {
                            continue;
                        }
                        rows.entry((mid, before.min(after), before.max(after)))
                            .or_insert_with(|| vec![(mid, h1 + h2), (before, -h2), (after, -h1)]);
                    }
                }
            }
        }
        rows
    }

    /// Patch `p`'s normal equations over its own net.
    fn local(&self, p: usize, at: &[(f64, f64)], [ku, kv]: &[KnotVector; 2]) -> Local {
        let (nu, nv) = (ku.control_point_count(), kv.control_point_count());
        let band = DEGREE * nv + DEGREE;
        let mut lower = vec![vec![0.0; band + 1]; nu * nv];
        let mut rhs = vec![[0.0; 3]; nu * nv];
        let mine = &self.taken.of_patch[p];
        if !mine.is_empty() {
            let share = 1.0 / precise(mine.len());
            let mut terms: Vec<(usize, f64)> = Vec::with_capacity((DEGREE + 1) * (DEGREE + 1));
            for &i in mine {
                let (u, v) = at[i];
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
            let (gu, gv) = (derivative_grams(ku), derivative_grams(kv));
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
    ///
    /// A triple that is one of the `junctions` rows (a curve ending at a
    /// T-junction, the line on through it the patch's own) is held by that
    /// row, and pins its chain's ratio to the row's; a chain pinned to two
    /// ratios is left untied, its curves near tangent.
    fn tangency(
        &self,
        layout: &Layout,
        loose: &[Point],
        junctions: &Keyed,
    ) -> OgeomResult<Vec<Vec<(usize, f64)>>> {
        let topo = self.topo;
        // Each held triple: (curve, B, L, R) slots.
        let mut triples: Vec<(usize, usize, usize, usize)> = Vec::new();
        for c in 0..topo.users.len() {
            if !topo.held[c] {
                continue;
            }
            let n = layout.along(c);
            let beside = |user: (usize, usize), q: usize| self.beside(layout, user, c, q);
            for q in 0..n {
                triples.push((
                    c,
                    layout.curve_slot(topo, c, q),
                    beside(topo.users[c][0], q),
                    beside(topo.users[c][1], q),
                ));
            }
        }
        // Each chain's pin, as the ratio's numerator and denominator, and
        // whether two pins disagree.
        let mut pins: Vec<Option<(f64, f64)>> = vec![None; topo.chains];
        let mut untied = vec![false; topo.chains];
        let mut sums = vec![(0.0, 0.0); topo.chains];
        for &(c, b, l, r) in &triples {
            let Some((chain, inverse)) = topo.chain[c] else {
                continue;
            };
            if let Some(row) = junctions.get(&(b, l.min(r), l.max(r))) {
                // The row reads `a B − x L − y R`, so `λ = y / x`.
                let of = |slot: usize| row.iter().find(|t| t.0 == slot).map_or(0.0, |t| -t.1);
                let (x, y) = (of(l), of(r));
                let pin = if inverse { (x, y) } else { (y, x) };
                match pins[chain] {
                    None => pins[chain] = Some(pin),
                    Some(held) => {
                        if held.0 * pin.1 != pin.0 * held.1 {
                            untied[chain] = true;
                        }
                    }
                }
                continue;
            }
            let (near, far) = (loose[b].distance(loose[l]), loose[r].distance(loose[b]));
            let (num, den) = if inverse { (far, near) } else { (near, far) };
            sums[chain].0 += num;
            sums[chain].1 += den;
        }
        let mut ratio = Vec::with_capacity(topo.chains);
        for (chain, &(num, den)) in sums.iter().enumerate() {
            if let Some((num, den)) = pins[chain] {
                // Opposite parities somewhere are served by a ratio of one
                // alone.
                if topo.conflicted[chain] && num != den {
                    untied[chain] = true;
                }
                ratio.push(num / den);
                continue;
            }
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
                if untied[chain] || junctions.contains_key(&(b, l.min(r), l.max(r))) {
                    return None;
                }
                let lambda = if inverse {
                    1.0 / ratio[chain]
                } else {
                    ratio[chain]
                };
                Some(vec![(b, 1.0 + lambda), (l, -1.0), (r, -lambda)])
            })
            .collect())
    }

    /// The slot of the pole one row into patch side `user` beside curve
    /// `c`'s pole `q`.
    fn beside(&self, layout: &Layout, user: (usize, usize), c: usize, q: usize) -> usize {
        let w = layout.walk_index(self.patches, user, c, q);
        let (_, inner) = side_pole(user.1, w, layout.dims(user.0));
        layout.slot(self.topo, self.patches, user.0, inner)
    }

    /// The rows holding tangent every curve into a corner of an odd number
    /// of patches whose curves all run on to the border.
    ///
    /// Across such a curve the cross derivatives of its two patches are
    /// `D₁(t)` and `D₂(t)` (the differences of the poles one row in from
    /// the curve's, on the curve's knots), and the patches share a tangent
    /// plane wherever `D₁ + λ D₂ + γ T = 0`, `T` the curve's own derivative
    /// and `λ > 0`. At the corner `D₁` and `D₂` run along the neighbouring
    /// curves, so there the condition is a relation among the first poles
    /// of the curves along each of them: they lie in one plane through the
    /// corner. With `λ` and `γ` constant along the curve the left side is
    /// a cubic spline with double knots, `C¹` where the curve's knots are,
    /// so it is nought throughout where it is nought at the corner and at
    /// two points inside each span and the far end, which is as many points
    /// as it has coefficients: an exact, linear condition on the poles,
    /// though no pole-by-pole one.
    ///
    /// The relations at the corner come from the fit `loose`: its first
    /// poles projected on the plane their consecutive cross products face,
    /// each curve's `λ` and `γ` the unique relation among its own first
    /// pole and the two beside it there. Those of `n` curves round a
    /// corner have `n − 2` independent ones, so `n − 2` consecutive ones
    /// are held and the other two follow. The cross derivatives at the
    /// corner differentiated along each curve bring in each patch's twist
    /// pole, once from either side, in a cycle that an odd number of
    /// positive ratios never makes singular; at an even number it can
    /// (the vertex enclosure condition), and a curve running on to a corner
    /// of four patches would close such a cycle there too, so neither is
    /// held. A corner whose relations do not give every curve a positive
    /// ratio is left near tangent.
    fn corner_rows(
        &self,
        layout: &Layout,
        loose: &[Point],
        curve_knots: &[KnotVector],
    ) -> Vec<Vec<(usize, f64)>> {
        let topo = self.topo;
        let mut rows = Vec::new();
        for (v, curves) in &topo.odd_corners {
            let v = *v;
            // Each curve's pole next to the corner, and the poles beside it
            // there in its first and second patch: those of the curves
            // either side.
            let near: Vec<(usize, usize, usize)> = curves
                .iter()
                .map(|&c| {
                    let n = layout.along(c);
                    let (at, next) = if topo.curve_ends[c].0 == v {
                        (0, 1)
                    } else {
                        (n - 1, n - 2)
                    };
                    let [first, second] = [topo.users[c][0], topo.users[c][1]];
                    (
                        self.beside(layout, first, c, at),
                        layout.curve_slot(topo, c, next),
                        self.beside(layout, second, c, at),
                    )
                })
                .collect();
            // The curves in turn round the corner.
            let mut order = vec![0];
            while order.len() < curves.len() {
                let last = near[order[order.len() - 1]];
                let Some(next) = (0..curves.len())
                    .find(|&i| !order.contains(&i) && (near[i].1 == last.0 || near[i].1 == last.2))
                else {
                    break;
                };
                order.push(next);
            }
            let closes = {
                let (first, last) = (near[order[0]], near[order[order.len() - 1]]);
                order.len() == curves.len() && (first.0 == last.1 || first.2 == last.1)
            };
            if !closes {
                continue;
            }
            let e = |slot: usize| loose[slot] - loose[v];
            let mut normal = Vector::ZERO;
            for (k, &i) in order.iter().enumerate() {
                let j = order[(k + 1) % order.len()];
                normal += e(near[i].1).cross(e(near[j].1));
            }
            let first = e(near[order[0]].1);
            let x = first - normal * (first.dot(normal) / normal.dot(normal));
            let y = normal.cross(x);
            let flat = |slot: usize| {
                let d = e(slot);
                (d.dot(x), d.dot(y))
            };
            // Each curve's `λ` and `γ` at the corner: `e₁ + λ e₂ + g e = 0`
            // in the plane.
            let mut weights = Vec::with_capacity(curves.len());
            for &(left, own, right) in &near {
                let ((lx, ly), (ox, oy), (rx, ry)) = (flat(left), flat(own), flat(right));
                let det = rx.mul_add(oy, -(ry * ox));
                let lambda = (ly.mul_add(ox, -(lx * oy))) / det;
                let g = (lx.mul_add(ry, -(ly * rx))) / det;
                if !(lambda > 0.0 && lambda.is_finite() && g.is_finite()) {
                    break;
                }
                weights.push((lambda, g));
            }
            if weights.len() < curves.len() {
                continue;
            }
            for &i in &order[..curves.len() - 2] {
                let (left, own, right) = near[i];
                let (lambda, g) = weights[i];
                rows.push(vec![
                    (left, 1.0),
                    (right, lambda),
                    (own, g),
                    (v, -(1.0 + lambda + g)),
                ]);
            }
            for (i, &c) in curves.iter().enumerate() {
                rows.extend(self.along_rows(layout, loose, c, v, weights[i], &curve_knots[c]));
            }
        }
        rows
    }

    /// The rows holding `D₁ + λ D₂ + γ T` nought along curve `c` past its
    /// end at `v`: at two points inside each span and at its other end.
    /// `weights` are `λ` and the relation's `g` at `v`, `γ T` there being
    /// `g` times the first pole's difference from `v`. `γ` runs linearly
    /// from there to the other end, which keeps `γ T` a cubic spline on the
    /// curve's knots; its value at the other end is the one that leaves the
    /// fit `loose` least out of tangency at the same points.
    fn along_rows(
        &self,
        layout: &Layout,
        loose: &[Point],
        c: usize,
        v: usize,
        (lambda, g): (f64, f64),
        knots: &KnotVector,
    ) -> Vec<Vec<(usize, f64)>> {
        let topo = self.topo;
        let t = knots.knots();
        let n = layout.along(c);
        let starts = topo.curve_ends[c].0 == v;
        // `T` at `v` is `3 / h` times the first pole's difference from
        // `v`, `h` the knot span there, and the other way about at the end.
        let near = if starts {
            g * (t[DEGREE + 1] - t[DEGREE]) / precise(DEGREE)
        } else {
            -g * (t[n] - t[n - 1]) / precise(DEGREE)
        };
        let [first, second] = [topo.users[c][0], topo.users[c][1]];
        let mut at = Vec::new();
        for w in knots.distinct().windows(2) {
            let (a, b) = (w[0].0, w[1].0);
            at.push((b - a).mul_add(1.0 / 3.0, a));
            at.push((b - a).mul_add(2.0 / 3.0, a));
        }
        at.push(if starts { 1.0 } else { 0.0 });
        // At each point, the slots and their weights in `D₁ + λ D₂`, in
        // `T`, and the share of the way from `v`.
        let terms: Vec<(Weighted, Weighted, f64)> = at
            .iter()
            .map(|&x| {
                let span = knots.span_unchecked(x);
                let basis = knots.basis_derivatives(span, x, 1);
                let mut cross = Vec::with_capacity(3 * (DEGREE + 1));
                let mut slope = Vec::with_capacity(DEGREE + 1);
                for i in 0..=DEGREE {
                    let q = span - DEGREE + i;
                    let (value, rate) = (basis[0][i], basis[1][i]);
                    let own = layout.curve_slot(topo, c, q);
                    cross.push((self.beside(layout, first, c, q), value));
                    cross.push((self.beside(layout, second, c, q), lambda * value));
                    cross.push((own, -(1.0 + lambda) * value));
                    slope.push((own, rate));
                }
                (cross, slope, if starts { x } else { 1.0 - x })
            })
            .collect();
        let sum = |row: &[(usize, f64)]| {
            row.iter()
                .fold(Vector::ZERO, |acc, &(s, w)| acc + loose[s].to_vector() * w)
        };
        // `γ = (1 − s) γ₀ + s γ₁`: the least squares `γ₁` over the points.
        let (mut num, mut den) = (0.0, 0.0);
        for (cross, slope, share) in &terms {
            let tangent = sum(slope);
            let rest = sum(cross) + tangent * ((1.0 - share) * near);
            let along = tangent * *share;
            num -= rest.dot(along);
            den += along.dot(along);
        }
        let far = if den > 0.0 { num / den } else { near };
        terms
            .iter()
            .map(|(cross, slope, share)| {
                let gamma = (1.0 - share).mul_add(near, share * far);
                let mut row: BTreeMap<usize, f64> = BTreeMap::new();
                for &(s, w) in cross {
                    *row.entry(s).or_insert(0.0) += w;
                }
                for &(s, w) in slope {
                    *row.entry(s).or_insert(0.0) += gamma * w;
                }
                row.into_iter().filter(|t| t.1 != 0.0).collect()
            })
            .collect()
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
        let reach = (2.0 * precise(round.distances.len().max(1)).ln()).sqrt();
        let (mut any, mut strays) = (false, false);
        let mut missed = Vec::new();
        for (p, mine) in self.taken.of_patch.iter().enumerate() {
            let [ku, kv] = &round.patch_knots[p];
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
