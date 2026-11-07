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
//! - **a corner of other than four patches is held tangent where it
//!   can be.** The triples hold the parameter lines that cross a curve to
//!   one direction on both sides of it. At a corner of three patches (or
//!   five, or six, or the point a three-sided patch's missing side shrinks
//!   to) the lines crossing a curve there lean one way on one side and the
//!   other way on the other, and no triple holds; tangency there needs the
//!   cross derivatives to differ by a multiple of the curve's own tangent,
//!   `D₁ + λ D₂ + γ T = 0`. At the corner that puts the curves' first poles
//!   in one plane, at ratios read off the fit without the constraints. A
//!   corner of an odd number of patches whose curves all run on to the
//!   network's border is first held with `λ` constant and `γ` linear along
//!   each whole curve: a cubic spline on the curve's knots, nought where it
//!   is nought at as many points as it has coefficients, an exact linear
//!   condition though not a pole-by-pole one. It binds each curve to a
//!   single cubic where `γ` is not nought at its knots. Otherwise, or
//!   where that misses the points, a corner is held by its spans: the knot
//!   past each curve's span next to the corner is doubled, that span is a
//!   quadratic, and `γ` is a quadratic over it with a double root at its
//!   far end and nought past it. Past that span the curve is held by the
//!   pole triples at its corner's ratio, so its far end may be on the
//!   border, a corner of four patches (an edge's fillet attached to a
//!   rounded box corner) or a junction, and spans raised there shrink the
//!   quadratic span. A junction pins the curve to its patch's knot
//!   spacings; the line through the junction is then held at the corner's
//!   ratio instead, tangent but not `C¹`. A curve may also run on to
//!   another corner held by its spans, which pins it to a second ratio:
//!   the first poles round a corner may be scaled to give its curves other
//!   ratios, so both corners are scaled to the geometric mean of the two,
//!   moved as little as closes every cycle of a corner's scales all of
//!   whose curves are asked a ratio. Where neither way meets the points, a
//!   corner of an odd number of patches, or an even number above four,
//!   whose curves each run on to the border or to another such corner has
//!   its gluing spread along its curves: `λ` constant and `γ` linear along
//!   each whole curve, from the corner's `g` to a value read off the fit
//!   at the border or to the other corner's, every knot inside the curve's
//!   strip doubled and the curves along it held `C²` across each. The
//!   curve is then as free to refine as any, its patches only `C¹` across
//!   those knots. At an even corner the twists round it close in a cycle
//!   that is always singular, which leaves one condition on the curves'
//!   second poles, held with the rest. This is the gluing
//!   a Coons patch's parameters carry, `γ` falling along the whole curve,
//!   which the points at their seed's parameters admit where `γ` nought
//!   past one span does not. So a corner of three patches or more,
//!   other than four, whose curves each run on to the border, a corner of
//!   four patches, a junction or another such corner is held, every curve
//!   into it tangent to rounding, as long as the fit still comes within
//!   the tolerance: the one ratio along each curve these conditions allow
//!   is a restriction the points at their seed's parameters need not
//!   admit. Where no way of holding it does, and at any other corner of
//!   other than four patches (an apex, or a junction), a curve with such a
//!   corner at an end is fitted without the triples: near tangent as the
//!   points are smooth, its angle measured into [`Network::worst_tangent`].
//!   Every such corner's widest angle between the patches' normals is
//!   reported in [`Network::corners`].
//! - **a fan is held tangent where it can be.** At a point every patch
//!   there shrinks its missing side to, both patches beside each curve
//!   into it shrink, so the parameter lines crossing the curve turn
//!   neither way and its pole triples hold it from end to end. Each
//!   patch's normal at the point is the limit `S_d × S_dw`, `d` across
//!   its missing side and `w` along it: from the point to the row of poles
//!   one in from that side, crossed with the row's own derivative. The
//!   patches share one tangent plane there where every patch's such row
//!   lies in one plane through the point; the rows are held to the layout
//!   they take in the fit without the constraints, projected on a plane
//!   through the point, its orientation, scale and shear left free. A fit
//!   holding a fan is kept only where every patch's normal at the point is
//!   defined all along its missing side and points one way, which the
//!   Bernstein coefficients of `(R − V) × R'` bound, `R` the row and `V`
//!   the point. Each way of holding the corners is first tried with the
//!   fans held, and only where none of them comes within the tolerance are
//!   the fans left near tangent. A point where a three-sided patch meets
//!   four-sided ones is left near tangent: there the three-sided patch's
//!   cross derivative vanishes and the other's does not, which no ratio
//!   fixed along the curve serves.
//! - **a side may be a chain of curves.** A curve may end on the middle of
//!   a patch's side (a T-junction), the side then a chain of curves meeting
//!   there ([`ChainedPatch`]). The patch takes a knot of multiplicity three
//!   at the junction, so its boundary row is the chain's poles, shared
//!   exactly, and the side opposite crosses the same junction: it is
//!   split alike, or is one curve whose knots carry the junction too. The
//!   patch is held `C¹` across the parameter line between the two, pole by
//!   pole, at the ratio of the knot spacings either side, or tangent at a
//!   corner's ratio where a curve from a corner held by its spans ends
//!   there. That line runs straight on from the curve ending at the
//!   junction, so that curve is held tangent at that same ratio, and the
//!   chain's curves either side of the junction share theirs, as at a
//!   corner of four patches.
//! - **spans are raised per strip.** The curves opposite each other in a
//!   patch share knots, and so do all the curves of a strip of patches
//!   joined through their opposite sides. Where a patch misses the points
//!   by more than the tolerance its strips gain spans, the way
//!   [`fit_surface_scattered_over`](ogeom_geom::fit::fit_surface_scattered_over)
//!   refines a single patch.
//! - **points may be reparameterised.** Once refinement ends,
//!   [`fit_network_chained`] can move every point's parameters to its
//!   nearest point on its patch and fit the patches again, a fit kept only
//!   where it comes nearer the points. Where every way of holding the
//!   corners misses, and so does the fit near tangent, each held fit that
//!   missed is finished that way at its spans before the near-tangent one
//!   is kept: a point whose nearest point lies on a side its patch shares
//!   is taken by the patch across that side where that one is nearer, as
//!   the seeds take each point by the nearest. A fit is kept only while it
//!   holds the corners, takes off more of the squared distances than twice
//!   the noise's variance, read once off the fit at the seed's parameters,
//!   and leaves the worst miss that is not a stray no worse, for at most
//!   eight steps; the first that comes within the tolerance is the
//!   network.
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

/// The most steps a held fit that missed is finished by. Re-taking the
/// points fixes the few a seed took by the wrong patch within a step or
/// two; past that, each step is a correction converging linearly, which
/// creeps the points under the tolerance while the patches' corners and
/// free borders, which no point beyond them holds, drift off the shape.
const FINISH_STEPS: usize = 8;

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
/// junctions from the side opposite: that side is split too, or is one
/// curve, which then carries the junctions to the patch beyond it. A curve
/// so carrying junctions may itself be one of the curves of a split side
/// (nested junctions): every side then crosses as many strips as the
/// curves along it carry between them, and the side opposite as many, the
/// `i`-th along one direction opposite the `i`-th along it on the other.
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
/// rather than a miss. Where corners of other than four patches can be
/// held tangent, the network is refined from the start with each way of
/// holding them the [module documentation](self) describes in turn, until
/// one comes within the tolerance; failing that, with those corners near
/// tangent. Where that misses too, each held way is finished at its spans,
/// its points re-taken at their nearest points on it, as the module
/// documentation describes, and the first to come within the tolerance is
/// kept. The points a patch fits ([`PatchFit::points`]) are then those it
/// takes as finished.
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
/// opposite sides of a patch, one of them split, cross different numbers
/// of strips, the split sides carry their junctions round a loop back to
/// themselves, or the curves of a split side meet one another, or the
/// patch's corners, other than at their junctions.
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
        corners: Corners::Loose,
        fans: false,
        tol,
    };
    // A corner held exactly tangent binds its curves' cross derivatives to
    // a relation with one ratio along each, which the points at the seed's
    // parameters need not admit within the tolerance. Each way of holding
    // the corners is refined in turn until one comes within it, the corners
    // left near tangent last. Holding a corner by its spans holds every
    // corner held whole and more besides; where it holds no more, the
    // corners are held whole first.
    let mut tries = Vec::with_capacity(4);
    let whole = !topo.odd_corners.is_empty();
    let spanned = !topo.span_corners.is_empty();
    let more = topo
        .span_corners
        .iter()
        .any(|(v, _)| !topo.odd_corners.iter().any(|(w, _)| w == v));
    if spanned && more {
        tries.push(Corners::Spans);
    }
    if whole {
        tries.push(Corners::Whole);
    }
    if spanned && !more {
        tries.push(Corners::Spans);
    }
    if !topo.spread_corners.is_empty() {
        tries.push(Corners::Spread);
    }
    tries.push(Corners::Loose);
    // Where there are fans, each way is first tried with them held too,
    // the near-tangent corners' way only kept where it comes within the
    // tolerance with them held.
    let mut tries: Vec<(Corners, bool)> = tries.into_iter().map(|c| (c, false)).collect();
    if !topo.fans.is_empty() {
        let fanned: Vec<(Corners, bool)> = tries.iter().map(|&(c, _)| (c, true)).collect();
        tries.splice(0..0, fanned);
    }
    let mut settled = None;
    let mut missed = Vec::new();
    for (corners, fans) in tries {
        fit.corners = corners;
        fit.fans = fans;
        let (spans, current, ended) = fit.refine(options)?;
        if (ended == Ended::Met && current.held) || (corners == Corners::Loose && !fans) {
            settled = Some((spans, current, ended));
            break;
        }
        if current.held {
            missed.push((corners, fans, spans, current));
        }
    }
    let Some((mut spans, mut current, mut ended)) = settled else {
        ogeom_bail!(Numeric, "no way of holding the corners was tried");
    };
    // Where the corners near tangent miss too, each held try that missed
    // is finished at its spans, in the order they were tried: the points
    // re-taken at their nearest points on the fit.
    let mut retaken = None;
    if ended != Ended::Met {
        for (corners, fans, at, round) in missed {
            let held = Fit {
                corners,
                fans,
                ..fit
            };
            if let Some((round, again)) = held.finish(&at, round, tolerance)? {
                fit.corners = corners;
                fit.fans = fans;
                (spans, current, ended) = (at, round, Ended::Met);
                retaken = Some(again);
                break;
            }
        }
    }
    let fit = Fit {
        taken: retaken.as_ref().unwrap_or(&taken),
        ..fit
    };
    // Parameter correction, at the spans refinement settled on. Refinement
    // fits every try at the seed's parameters, so that tries compare as
    // fits of the same points at the same parameters rather than as
    // corrections that have run further or less far.
    for _ in 0..corrections {
        let next = fit.round(&spans, &current.feet)?;
        let kept_within = ended != Ended::Met || next.distances.iter().all(|d| *d <= tolerance);
        if next.square() < current.square() && kept_within && next.held {
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
    for (&p, &d) in fit.taken.patch.iter().zip(&current.distances) {
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

/// Which knots inside a strip are doubled.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct Doubled {
    /// Each end's first knot in, in the strip's own direction.
    ends: [bool; 2],
    /// Every knot inside it.
    every: bool,
}

/// How many knots each strip takes beyond one per span: one for each knot
/// inside it that is doubled.
fn extra_knots(spans: usize, doubled: Doubled) -> usize {
    let ends = doubled.ends;
    match spans {
        1 => 0,
        _ if doubled.every => spans - 1,
        2 => usize::from(ends[0] || ends[1]),
        _ => usize::from(ends[0]) + usize::from(ends[1]),
    }
}

/// The knots of a run of strips over the unit interval, and the bounds of
/// each strip's share of it: each strip's spans even over its share, the
/// shares in proportion to the strips' weights, and a knot of multiplicity
/// the degree at each junction between two strips, where the run's poles
/// pass through the junction's. `flips` says which strips the run crosses
/// against their own direction, and `doubled` which knots inside each
/// strip are doubled (none where it is empty).
fn run_knots(
    strips: &[usize],
    flips: &[bool],
    spans: &[usize],
    weights: &[f64],
    doubled: &[Doubled],
) -> OgeomResult<(KnotVector, Vec<f64>)> {
    let ends = |s: usize| doubled.get(s).copied().unwrap_or_default();
    if let [one] = strips
        && extra_knots(spans[*one], ends(*one)) == 0
    {
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
        let (
            n,
            Doubled {
                ends: [first, last],
                every,
            },
        ) = (spans[s], ends(s));
        let flip = flips.get(i).copied().unwrap_or(false);
        for k in 1..n {
            let x = if strips.len() == 1 {
                precise(k) / precise(n)
            } else {
                (b - a).mul_add(precise(k) / precise(n), a)
            };
            // The knot's place counted in the strip's own direction.
            let own = if flip { n - k } else { k };
            let twice = every || (own == 1 && first) || (own == n - 1 && last);
            knots.extend(core::iter::repeat_n(x, 1 + usize::from(twice)));
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
    /// The layout over `spans` per strip, the knots inside each strip
    /// `doubled` names doubled.
    fn new(topo: &Topology, spans: &[usize], doubled: &[Doubled]) -> Self {
        let controls: Vec<usize> = spans
            .iter()
            .zip(doubled)
            .map(|(&s, &ends)| s + DEGREE + extra_knots(s, ends))
            .collect();
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

/// The strips along the walk of each curve of side `k` of patch `p`.
fn side_parts(topo: &Topology, patches: &[ChainedPatch], (p, k): (usize, usize)) -> Vec<usize> {
    patches[p].sides[k]
        .iter()
        .map(|s| topo.curve_strip[s.curve].len())
        .collect()
}

/// The share of side `k`'s walk, from its start, that its `i`-th curve
/// covers, `parts` the strips each of its curves crosses in walking order
/// and `bounds` the bounds of the shares of the strips along the side's
/// direction.
fn walk_share(k: usize, i: usize, parts: &[usize], bounds: &[f64]) -> (f64, f64) {
    if parts.len() < 2 {
        return (0.0, 1.0);
    }
    let m = bounds.len() - 1;
    let (from, to) = (
        parts[..i].iter().sum::<usize>(),
        parts[..=i].iter().sum::<usize>(),
    );
    if k < 2 {
        (bounds[from], bounds[to])
    } else {
        (1.0 - bounds[m - from], 1.0 - bounds[m - to])
    }
}

/// Which curve of side `k` holds the point at `x` along the side's
/// direction, `parts` the strips each of its curves crosses in walking
/// order and `bounds` the bounds of the shares of its strips.
fn piece_at(k: usize, x: f64, parts: &[usize], bounds: &[f64]) -> usize {
    if parts.len() < 2 {
        return 0;
    }
    let m = bounds.len() - 1;
    let strip = bounds[1..m].iter().filter(|&&b| x > b).count();
    let walked = if k < 2 { strip } else { m - 1 - strip };
    let mut end = 0;
    for (i, n) in parts.iter().enumerate() {
        end += n;
        if walked < end {
            return i;
        }
    }
    parts.len() - 1
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
#[derive(Clone)]
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
                let (knots, _) = run_knots(&topo.curve_strip[c], &[], &spans, weights, &[])?;
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
            let (ku, bu) = run_knots(su, &[], &spans, weights, &[])?;
            let (kv, bv) = run_knots(sv, &[], &spans, weights, &[])?;
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
                            run_knots(&topo.curve_strip[piece.curve], &[], &spans, weights, &[])?;
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
        // Whether each curve of each side is on the network's border, and
        // the strips each crosses.
        let border: Vec<[(Vec<bool>, Vec<usize>); 4]> = (0..patches.len())
            .map(|p| {
                core::array::from_fn(|k| {
                    (
                        patches[p].sides[k]
                            .iter()
                            .map(|s| topo.users[s.curve].len() < 2)
                            .collect(),
                        side_parts(topo, patches, (p, k)),
                    )
                })
            })
            .collect();
        let past = |p: usize, (u, v): (f64, f64)| {
            let free = |k: usize, x: f64| {
                let (sides, parts) = &border[p][k];
                let i = piece_at(k, x, parts, &self.bounds[p][k % 2]);
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

/// Which corners of other than four patches are held exactly tangent.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Corners {
    /// None: every curve into one is fitted near tangent.
    Loose,
    /// The corners of an odd number of patches whose curves all run on to
    /// the border, `γ` linear along each whole curve.
    Whole,
    /// The corners of an odd number of patches, or an even number above
    /// four, held by the span of each curve next to them: that span a
    /// quadratic, `γ` quadratic over it and nought past it, so that past it
    /// the curve is held by the pole triples and may run on to a corner of
    /// four patches, a junction or another such corner.
    Spans,
    /// The corners of an odd number of patches, or an even number above
    /// four, whose curves all run on to the border or to another such
    /// corner, `γ` linear along each whole curve, every knot inside the
    /// curve's strip doubled and the curve held `C²` across each.
    Spread,
}

/// The fitting problem, fixed across rounds.
struct Fit<'a> {
    topo: &'a Topology,
    patches: &'a [ChainedPatch],
    taken: &'a Taken,
    weights: &'a [f64],
    smoothing: f64,
    /// Which corners of other than four patches are held exactly tangent,
    /// and how.
    corners: Corners,
    /// Whether the fans are held exactly tangent.
    fans: bool,
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
    /// Whether every corner the way of holding them names is held.
    held: bool,
}

impl Round {
    /// The sum of the squared distances.
    fn square(&self) -> f64 {
        self.distances.iter().map(|d| d * d).sum()
    }
}

/// Slots and their weights.
type Weighted = Vec<(usize, f64)>;

/// A constraint row: slots and their coefficients, summing to nought.
type Row = Vec<(usize, f64)>;

/// How the curves into a corner of other than four patches relate there.
struct Relation {
    /// The corner's vertex.
    vertex: usize,
    /// Each curve's pole next to the corner and the poles beside it there in
    /// its first and second patch.
    near: Vec<(usize, usize, usize)>,
    /// The curves in turn round the corner.
    order: Vec<usize>,
    /// Each curve's `λ` and `g` there.
    weights: Vec<(f64, f64)>,
    /// Each curve's first pole off the corner, in the plane the relation
    /// is read in.
    flat: Vec<(f64, f64)>,
    /// The curves either side of each: the one whose first pole lies in
    /// its first patch, and the one in its second.
    sides: Vec<(usize, usize)>,
}

/// Each curve's `λ` and `g` at a corner, `e₁ + λ e₂ + g e = 0` among its
/// own first pole `e` and those of the curves either side in the plane;
/// `None` where a ratio is not positive.
fn plane_weights(flat: &[(f64, f64)], sides: &[(usize, usize)]) -> Option<Vec<(f64, f64)>> {
    let mut weights = Vec::with_capacity(flat.len());
    for (i, &(l, r)) in sides.iter().enumerate() {
        let ((lx, ly), (ox, oy), (rx, ry)) = (flat[l], flat[i], flat[r]);
        let det = rx.mul_add(oy, -(ry * ox));
        let lambda = (ly.mul_add(ox, -(lx * oy))) / det;
        let g = (lx.mul_add(ry, -(ly * rx))) / det;
        if !(lambda > 0.0 && lambda.is_finite() && g.is_finite()) {
            return None;
        }
        weights.push((lambda, g));
    }
    Some(weights)
}

impl Relation {
    /// Scale the curves' first poles so each curve of `wanted` takes its
    /// ratio there, then hold it to that ratio to the bit; false, the
    /// relation untouched, where a ratio is then not positive.
    ///
    /// Scaling the first pole of the curve on a curve's first side by `a`
    /// and that on its second by `b` scales its `λ` by `a / b`, and its own
    /// leaves `λ` alone: the ratios wanted are differences of the scales'
    /// logarithms along the curves either side of each, each run of them
    /// settled with the mean logarithm nought. Those pairs close in one
    /// cycle round an odd corner and two round an even one; where every
    /// curve of a cycle is wanted, the ratios wanted have been [closed
    /// round it](close_cycles), and its last pair follows from the others.
    fn rescale(&mut self, wanted: &[(usize, f64)]) -> bool {
        let n = self.flat.len();
        // Each scale's logarithm against its run's first: a union-find
        // whose links carry the difference.
        let mut parent: Vec<usize> = (0..n).collect();
        let mut offset = vec![0.0; n];
        let root = |parent: &[usize], offset: &[f64], mut x: usize| {
            let mut sum = 0.0;
            while parent[x] != x {
                sum += offset[x];
                x = parent[x];
            }
            (x, sum)
        };
        for &(i, target) in wanted {
            let (l, r) = self.sides[i];
            // `log a_l − log a_r = log(target / λ)`.
            let step = (target / self.weights[i].0).ln();
            let ((rl, ol), (rr, or)) = (root(&parent, &offset, l), root(&parent, &offset, r));
            if rl == rr {
                continue;
            }
            parent[rr] = rl;
            offset[rr] = ol - step - or;
        }
        let logs: Vec<(usize, f64)> = (0..n).map(|x| root(&parent, &offset, x)).collect();
        let mut means = vec![(0.0, 0usize); n];
        for &(r, log) in &logs {
            means[r].0 += log;
            means[r].1 += 1;
        }
        let flat: Vec<(f64, f64)> = self
            .flat
            .iter()
            .zip(&logs)
            .map(|(&(x, y), &(r, log))| {
                let scale = (log - means[r].0 / precise(means[r].1)).exp();
                (x * scale, y * scale)
            })
            .collect();
        let Some(mut weights) = plane_weights(&flat, &self.sides) else {
            return false;
        };
        for &(i, target) in wanted {
            weights[i].0 = target;
        }
        self.flat = flat;
        self.weights = weights;
        true
    }

    /// The rows holding the relation: `n − 2` consecutive curves' own, the
    /// other two following from them.
    fn rows(&self) -> Vec<Vec<(usize, f64)>> {
        self.order[..self.order.len() - 2]
            .iter()
            .map(|&i| {
                let (left, own, right) = self.near[i];
                let (lambda, g) = self.weights[i];
                vec![
                    (left, 1.0),
                    (right, lambda),
                    (own, g),
                    (self.vertex, -(1.0 + lambda + g)),
                ]
            })
            .collect()
    }
}

/// A corner held by its spans, with its relation where it has one.
struct Related {
    vertex: usize,
    curves: Vec<usize>,
    relation: Option<Relation>,
}

/// A constraint row keyed by its slots: the middle one, then the other two
/// in order.
type Keyed = BTreeMap<(usize, usize, usize), Vec<(usize, f64)>>;

/// The parameter lines through the junctions: each a patch's rows across
/// one junction, joined with the rows of the patches beyond the curves
/// that carry the junction on.
struct Lines {
    /// Each junction row's line, and whether the row as kept runs against
    /// the line's direction.
    of: BTreeMap<(usize, usize, usize), (usize, bool)>,
    /// The number of lines.
    count: usize,
    /// The lines met both ways, which hold their patches' own ratio alone.
    torn: Vec<bool>,
}

/// What the tangency rows come to: the rows, which chains are left untied,
/// the junction rows held at a ratio of their own in place of their
/// patches' knot spacings, and the triples held.
struct Tangency {
    rows: Vec<Row>,
    untied: Vec<bool>,
    relaxed: Keyed,
    /// Each triple held, by its curve's pole: the poles beside it in the
    /// curve's first and second patch, and `λ`.
    held: BTreeMap<usize, (usize, usize, f64)>,
}

/// The members of a class of chains or of junction lines, each with
/// whether its ratio is the inverse of the class's root's.
type Members = Vec<(usize, bool)>;

/// A class of chains and junction lines that corners pin: each corner's
/// pin (the corner, its curve, and the ratio it pins the root to), and its
/// members.
struct Pinned {
    pins: Vec<(usize, usize, f64)>,
    chains: Members,
    lines: Members,
}

/// A class of chains and junction lines a corner settles on one ratio:
/// the ratio of the class's root, and each chain and line with whether its
/// own is the inverse.
struct Settled {
    ratio: f64,
    chains: Members,
    lines: Members,
}

/// A curve of a corner asked the ratio of a settled class: the curve's
/// place round the corner, the class, whether the chain's ratio is the
/// class's inverse, and whether the curve's is its chain's.
#[derive(Debug, Clone, Copy)]
struct Asked {
    curve: usize,
    class: usize,
    odd: bool,
    inverse: bool,
}

/// A combination of the classes' log ratios with integer coefficients,
/// and a constant.
type Potential = (BTreeMap<usize, i64>, f64);

/// `a + sign b`.
fn potential_add(a: &Potential, b: &Potential, sign: i64) -> Potential {
    let mut out = a.0.clone();
    for (&f, &x) in &b.0 {
        *out.entry(f).or_insert(0) += sign * x;
    }
    out.retain(|_, x| *x != 0);
    (out, precise_signed(sign).mul_add(b.1, a.1))
}

/// Move the log ratios `logs` of the settled classes as little as lets
/// every corner meet every ratio asked of it; the classes that cannot be
/// met so.
///
/// A corner [rescaled](Relation::rescale) meets the ratios asked of its
/// curves through the differences of its scales' logarithms, one per
/// curve between the scales of the curves either side. Those pairs close
/// in one cycle round an odd corner and two round an even one, and where
/// every curve of a cycle is asked a ratio, the differences round it sum
/// to nought: the ratios asked, signed as the cycle meets each curve, sum
/// to the sum of the curves' own log ratios. Round an odd corner that sum
/// is nought, as the ratios multiply to one; round an even one it is the
/// cycle's own. Each such cycle is a linear condition with integer
/// coefficients on the classes' log ratios, so the ratios are moved to
/// the nearest that meet them all, by least squares over independent
/// ones. Conditions that depend on one another are met where all of them
/// are an odd corner's, their sums all nought; otherwise, or where the
/// integers outgrow their type, every class in a condition fails.
fn close_cycles(
    corners: &[Related],
    asked: &BTreeMap<usize, Vec<Asked>>,
    logs: &mut [Option<f64>],
) -> Vec<bool> {
    let mut failed = vec![false; logs.len()];
    // Each cycle's condition: the coefficients, the value they sum to,
    // and whether it is an even corner's.
    let mut conditions: Vec<(BTreeMap<usize, i64>, f64, bool)> = Vec::new();
    for (&k, links) in asked {
        let Some(relation) = &corners[k].relation else {
            continue;
        };
        let n = relation.flat.len();
        let even = n.is_multiple_of(2);
        // Each scale's log against its run's root: a union-find whose
        // links carry the difference.
        let mut parent: Vec<usize> = (0..n).collect();
        let mut offset: Vec<Potential> = vec![(BTreeMap::new(), 0.0); n];
        let root = |parent: &[usize], offset: &[Potential], mut x: usize| {
            let mut sum: Potential = (BTreeMap::new(), 0.0);
            while parent[x] != x {
                sum = potential_add(&sum, &offset[x], 1);
                x = parent[x];
            }
            (x, sum)
        };
        for a in links {
            let (l, r) = relation.sides[a.curve];
            // `log a_l − log a_r = ±x − log λ`, `x` the class's log ratio.
            let sign = if a.odd ^ a.inverse { -1 } else { 1 };
            let step: Potential = (
                BTreeMap::from([(a.class, sign)]),
                -relation.weights[a.curve].0.ln(),
            );
            let ((rl, ol), (rr, or)) = (root(&parent, &offset, l), root(&parent, &offset, r));
            let gap = potential_add(&potential_add(&ol, &step, -1), &or, -1);
            if rl == rr {
                // `ol − or − step` is nought: its coefficients against
                // minus its constant, nought round an odd corner.
                conditions.push((gap.0, if even { -gap.1 } else { 0.0 }, even));
            } else {
                parent[rr] = rl;
                offset[rr] = gap;
            }
        }
    }
    if conditions.is_empty() {
        return failed;
    }
    let fail_all = |failed: &mut Vec<bool>| {
        for (row, _, _) in &conditions {
            for &f in row.keys() {
                failed[f] = true;
            }
        }
    };
    // The independent conditions, by fraction-free elimination over the
    // integers: each kept row with its pivot column, zero in the pivot
    // columns of the rows kept before it.
    let columns = logs.len();
    let mut kept: Vec<(usize, Vec<i128>)> = Vec::new();
    let mut independent = Vec::new();
    let mut dependent = false;
    for (index, (row, _, _)) in conditions.iter().enumerate() {
        let mut dense = vec![0i128; columns];
        for (&f, &c) in row {
            dense[f] = i128::from(c);
        }
        for (pivot, base) in &kept {
            let x = dense[*pivot];
            if x == 0 {
                continue;
            }
            let y = base[*pivot];
            let mut divisor = 0i128;
            for (d, &b) in dense.iter_mut().zip(base) {
                let Some(v) = d
                    .checked_mul(y)
                    .and_then(|v| v.checked_sub(b.checked_mul(x)?))
                else {
                    fail_all(&mut failed);
                    return failed;
                };
                *d = v;
                divisor = gcd(divisor, v);
            }
            if divisor > 1 {
                for d in &mut dense {
                    *d /= divisor;
                }
            }
        }
        match dense.iter().position(|&d| d != 0) {
            Some(pivot) => {
                kept.push((pivot, dense));
                independent.push(index);
            }
            None => dependent = true,
        }
    }
    if dependent && conditions.iter().any(|c| c.2) {
        fail_all(&mut failed);
        return failed;
    }
    // The nearest log ratios meeting the independent conditions `C x = d`:
    // `x − Cᵀ (C Cᵀ)⁻¹ (C x − d)`.
    let rows: Vec<&(BTreeMap<usize, i64>, f64, bool)> =
        independent.iter().map(|&i| &conditions[i]).collect();
    let value = |f: usize| logs[f].unwrap_or(0.0);
    let mut gram = Lower::default();
    for (i, a) in rows.iter().enumerate() {
        for (j, b) in rows.iter().enumerate().take(i + 1) {
            let dot: i64 =
                a.0.iter()
                    .filter_map(|(f, x)| b.0.get(f).map(|y| x * y))
                    .sum();
            gram.add(i, j, precise_signed(dot));
        }
    }
    let residual: Vec<f64> = rows
        .iter()
        .map(|(row, d, _)| {
            row.iter()
                .map(|(&f, &x)| precise_signed(x) * value(f))
                .sum::<f64>()
                - d
        })
        .collect();
    let Some(solved) = solve_spd(rows.len(), &gram, &[residual]) else {
        fail_all(&mut failed);
        return failed;
    };
    let mut shift: BTreeMap<usize, f64> = BTreeMap::new();
    for ((row, _, _), y) in rows.iter().zip(&solved[0]) {
        for (&f, &x) in row {
            *shift.entry(f).or_insert(0.0) += precise_signed(x) * y;
        }
    }
    for (f, moved) in shift {
        if let Some(log) = logs[f].as_mut() {
            *log -= moved;
        }
    }
    failed
}

/// The greatest common divisor of two integers' magnitudes.
fn gcd(a: i128, b: i128) -> i128 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// A small signed count as a float.
#[allow(clippy::cast_precision_loss, reason = "coefficients far below 2^52")]
fn precise_signed(n: i64) -> f64 {
    n as f64
}

/// `γ` at curve `c`'s start, where `starts`, or its end, for a corner
/// there whose relation gives the curve `g`: `T` there is `3 / h` times
/// the first pole's difference from the corner, `h` the knot span there,
/// and the other way about at the end.
fn end_gamma(layout: &Layout, c: usize, g: f64, knots: &KnotVector, starts: bool) -> f64 {
    let t = knots.knots();
    let n = layout.along(c);
    if starts {
        g * (t[DEGREE + 1] - t[DEGREE]) / precise(DEGREE)
    } else {
        -g * (t[n] - t[n - 1]) / precise(DEGREE)
    }
}

/// Two points inside each span of `knots` and the end away from the
/// corner, at its start where `starts`.
fn along_points(knots: &KnotVector, starts: bool) -> Vec<f64> {
    let mut at = Vec::new();
    for w in knots.distinct().windows(2) {
        let (a, b) = (w[0].0, w[1].0);
        at.push((b - a).mul_add(1.0 / 3.0, a));
        at.push((b - a).mul_add(2.0 / 3.0, a));
    }
    at.push(if starts { 1.0 } else { 0.0 });
    at
}

/// `γ₁` in `γ = (1 − s) γ₀ + s γ₁`, `γ₀` being `near`: the least squares
/// value over `terms`, read in the fit `loose`.
fn least_far(loose: &[Point], terms: &[(Weighted, Weighted, f64)], near: f64) -> f64 {
    let sum = |row: &[(usize, f64)]| {
        row.iter()
            .fold(Vector::ZERO, |acc, &(s, w)| acc + loose[s].to_vector() * w)
    };
    let (mut num, mut den) = (0.0, 0.0);
    for (cross, slope, share) in terms {
        let tangent = sum(slope);
        let rest = sum(cross) + tangent * ((1.0 - share) * near);
        let along = tangent * *share;
        num -= rest.dot(along);
        den += along.dot(along);
    }
    if den > 0.0 { num / den } else { near }
}

/// The row `cross + gamma slope`, its slots merged.
fn gamma_row(cross: &[(usize, f64)], slope: &[(usize, f64)], gamma: f64) -> Row {
    let mut row: BTreeMap<usize, f64> = BTreeMap::new();
    for &(s, w) in cross {
        *row.entry(s).or_insert(0.0) += w;
    }
    for &(s, w) in slope {
        *row.entry(s).or_insert(0.0) += gamma * w;
    }
    row.into_iter().filter(|t| t.1 != 0.0).collect()
}

/// The rows holding curve `c` `C²` across each doubled knot inside it,
/// `knots` its knots: its second derivative the same either side, scaled
/// by the square of the width of the span before.
fn smooth_rows(layout: &Layout, topo: &Topology, c: usize, knots: &KnotVector) -> Vec<Row> {
    let t = knots.knots();
    let distinct = knots.distinct();
    let mut rows = Vec::new();
    let mut first = 0;
    for (k, &(x, m)) in distinct.iter().enumerate() {
        if k > 0 && k + 1 < distinct.len() && m == 2 {
            let (left, right) = (first - 1, first + 1);
            let width = x - t[left];
            let scale = width * width;
            let mut row: BTreeMap<usize, f64> = BTreeMap::new();
            for (span, sign) in [(left, scale), (right, -scale)] {
                let bend = knots.basis_derivatives(span, x, 2);
                for (i, d) in bend[2].iter().enumerate() {
                    *row.entry(layout.curve_slot(topo, c, span - DEGREE + i))
                        .or_insert(0.0) += sign * d;
                }
            }
            rows.push(row.into_iter().filter(|e| e.1 != 0.0).collect());
        }
        first += m;
    }
    rows
}

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
        let doubled = self.doubled();
        let layout = Layout::new(topo, spans, &doubled);
        let mut patch_knots = Vec::with_capacity(self.patches.len());
        let mut patch_bounds = Vec::with_capacity(self.patches.len());
        for ([su, sv], [fu, fv]) in topo.patch_strips.iter().zip(&topo.patch_flip) {
            let (ku, bu) = run_knots(su, fu, spans, self.weights, &doubled)?;
            let (kv, bv) = run_knots(sv, fv, spans, self.weights, &doubled)?;
            patch_knots.push([ku, kv]);
            patch_bounds.push([bu, bv]);
        }
        let curve_knots: Vec<KnotVector> = topo
            .curve_strip
            .iter()
            .zip(&topo.curve_flip)
            .map(|(strips, flips)| {
                run_knots(strips, flips, spans, self.weights, &doubled).map(|k| k.0)
            })
            .collect::<OgeomResult<_>>()?;
        let locals: Vec<Local> = ogeom_core::parallel::map_ordered(
            &(0..self.patches.len()).collect::<Vec<_>>(),
            |_, &p| self.local(p, at, &patch_knots[p]),
        );

        // Without the tangency, for the ratios across each curve; the
        // patches held smooth across their junctions all the same.
        let (junctions, lines) = self.junction_rows(&layout, spans);
        let loose_basis = if junctions.is_empty() {
            Basis::identity(layout.slots)
        } else {
            null_space(
                layout.slots,
                &junctions.values().cloned().collect::<Vec<_>>(),
            )
        };
        let loose = self.solve(&layout, &locals, &loose_basis)?;
        let junction_rows = |relaxed: &Keyed| -> Vec<Row> {
            junctions
                .iter()
                .map(|(key, row)| relaxed.get(key).unwrap_or(row).clone())
                .collect()
        };
        let mut rows: Vec<Row>;
        let mut held = true;
        let triples;
        match self.corners {
            Corners::Loose => {
                let held = self.tangency(&layout, &loose, &junctions, &lines, &mut [])?;
                rows = junction_rows(&held.relaxed);
                rows.extend(held.rows);
                triples = held.held;
            }
            Corners::Whole => {
                let held = self.tangency(&layout, &loose, &junctions, &lines, &mut [])?;
                rows = junction_rows(&held.relaxed);
                rows.extend(held.rows);
                triples = held.held;
                rows.extend(self.corner_rows(&layout, &loose, &curve_knots));
            }
            Corners::Spread => {
                let tangent = self.tangency(&layout, &loose, &junctions, &lines, &mut [])?;
                rows = junction_rows(&tangent.relaxed);
                rows.extend(tangent.rows);
                triples = tangent.held;
                match self.spread_rows(&layout, &loose, &curve_knots) {
                    Some(spread) => rows.extend(spread),
                    None => held = false,
                }
            }
            Corners::Spans => {
                let mut related: Vec<Related> = self
                    .topo
                    .span_corners
                    .iter()
                    .map(|(v, curves)| Related {
                        vertex: *v,
                        curves: curves.clone(),
                        relation: self.relate(&layout, &loose, *v, curves),
                    })
                    .collect();
                let held = self.tangency(&layout, &loose, &junctions, &lines, &mut related)?;
                rows = junction_rows(&held.relaxed);
                rows.extend(held.rows);
                triples = held.held;
                let untied = held.untied;
                let holding = &self.topo.spanned;
                for corner in &related {
                    let Some(relation) = &corner.relation else {
                        continue;
                    };
                    let tied = corner
                        .curves
                        .iter()
                        .all(|&c| holding.chain[c].is_some_and(|(chain, _)| !untied[chain]));
                    if !tied {
                        continue;
                    }
                    rows.extend(relation.rows());
                    for (i, &c) in corner.curves.iter().enumerate() {
                        rows.extend(self.span_rows(
                            &layout,
                            c,
                            corner.vertex,
                            relation.weights[i],
                            &curve_knots[c],
                        ));
                    }
                }
            }
        }
        if self.fans {
            match self.fan_rows(&layout, &loose, &triples) {
                Some(fan) => rows.extend(fan),
                None => held = false,
            }
        }
        let basis = null_space(layout.slots, &rows);
        let poles = self.solve(&layout, &locals, &basis)?;
        if self.fans && held {
            held = self.fans_regular(&layout, &patch_knots, &poles)?;
        }

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
            held,
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
    ///
    /// Also the lines the rows lie on: a patch's rows across one junction
    /// are one line, and lines sharing a row are one.
    fn junction_rows(&self, layout: &Layout, spans: &[usize]) -> (Keyed, Lines) {
        let topo = self.topo;
        let spacing = |s: usize| self.weights[s] / precise(spans[s]);
        let mut rows = Keyed::new();
        // Each row's first patch line and the pole before its middle there.
        let mut owner: BTreeMap<(usize, usize, usize), (usize, usize)> = BTreeMap::new();
        let mut joined = Vec::new();
        let mut patch_lines = 0;
        for p in 0..self.patches.len() {
            let (nu, nv) = layout.dims(p);
            for (dir, strips) in topo.patch_strips[p].iter().enumerate() {
                let across = if dir == 0 { nv } else { nu };
                let mut at = 0;
                for pair in strips.windows(2) {
                    at += layout.controls[pair[0]] - 1;
                    let (h1, h2) = (spacing(pair[0]), spacing(pair[1]));
                    let line = patch_lines;
                    patch_lines += 1;
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
                        let key = (mid, before.min(after), before.max(after));
                        match owner.get(&key) {
                            Some(&(first, first_before)) => {
                                joined.push((first, line, first_before != before));
                            }
                            None => {
                                owner.insert(key, (line, before));
                            }
                        }
                        rows.entry(key)
                            .or_insert_with(|| vec![(mid, h1 + h2), (before, -h2), (after, -h1)]);
                    }
                }
            }
        }
        let mut classes = topology::Classes::new(patch_lines);
        let mut torn_lines = Vec::new();
        for (a, b, odd) in joined {
            if !classes.join(a, b, odd) {
                torn_lines.push(a);
            }
        }
        let (named, count) = classes.numbered(0..patch_lines);
        let mut torn = vec![false; count];
        for line in torn_lines {
            if let Some(id) = named[line] {
                torn[id] = true;
            }
        }
        let of = owner
            .into_iter()
            .map(|(key, (line, _))| {
                let flip = classes.find(line).1;
                (key, (named[line].unwrap_or(0), flip))
            })
            .collect();
        (rows, Lines { of, count, torn })
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

    /// The curves held by the pole triples, as the corners and fans are
    /// held.
    fn holding(&self) -> &topology::Holding {
        match (self.corners == Corners::Spans, self.fans) {
            (true, false) => &self.topo.spanned,
            (false, false) => &self.topo.plain,
            (true, true) => &self.topo.spanned_fans,
            (false, true) => &self.topo.plain_fans,
        }
    }

    /// Which knots inside each strip are doubled: the first knot in from
    /// the end at each corner held by its spans, of each curve into it, and
    /// every knot inside the strip of each curve into a corner whose gluing
    /// is spread along its curves.
    fn doubled(&self) -> Vec<Doubled> {
        let topo = self.topo;
        let mut doubled = vec![Doubled::default(); topo.strips];
        if self.corners == Corners::Spans {
            for (v, curves) in &topo.span_corners {
                for &c in curves {
                    let starts = topo.curve_ends[c].0 == *v;
                    doubled[topo.curve_strip[c][0]].ends
                        [usize::from(starts == topo.curve_flip[c][0])] = true;
                }
            }
        }
        if self.corners == Corners::Spread {
            for (_, curves) in &topo.spread_corners {
                for &c in curves {
                    doubled[topo.curve_strip[c][0]].every = true;
                }
            }
        }
        doubled
    }

    /// The tangency rows: for each pole of each curve held tangent, the
    /// pole and the poles one row into either patch on one line,
    /// `(1 + λ) B − L − λ R = 0`, `L` in the curve's first patch, `λ` its
    /// chain's ratio of pole spacings in the fit `loose`; which chains are
    /// left untied; and the junction rows held at a ratio of their own.
    ///
    /// A triple that is one of the `junctions` rows (a curve ending at a
    /// T-junction, the line on through it the patch's own) is held by that
    /// row, and pins its chain's ratio to the row's. A curve into one of the
    /// `corners` pins its chain's ratio to its own `λ` there, and its two
    /// poles nearest the corner are held by the corner's rows instead; a
    /// corner without a relation unties its curves' chains. A chain pinned
    /// to two ratios is left untied, its curves near tangent, but where a
    /// corner pins it:
    ///
    /// - the lines through the junctions it meets are held at the corner's
    ///   ratio rather than their knot spacings', pole by pole as across a
    ///   curve, so their patches are tangent across them though not `C¹`;
    /// - where it meets several corners, it takes the geometric mean of
    ///   their ratios, moved as little as [closes](close_cycles) every
    ///   cycle of a corner's scales asked a ratio at each of its curves, and
    ///   each corner is [rescaled](Relation::rescale) to it.
    ///
    /// Chains and lines that pin one another are settled together, and
    /// where one of them cannot be, all of them are left as before.
    #[allow(
        clippy::too_many_lines,
        reason = "the pins, then the chains they settle"
    )]
    fn tangency(
        &self,
        layout: &Layout,
        loose: &[Point],
        junctions: &Keyed,
        lines: &Lines,
        corners: &mut [Related],
    ) -> OgeomResult<Tangency> {
        let topo = self.topo;
        let holding = self.holding();
        // The poles of each curve the corners hold.
        let mut skipped: Vec<[bool; 2]> = vec![[false; 2]; topo.users.len()];
        for corner in corners.iter() {
            for &c in &corner.curves {
                skipped[c][usize::from(topo.curve_ends[c].0 != corner.vertex)] = true;
            }
        }
        // Each held triple: (curve, B, L, R) slots.
        let mut triples: Vec<(usize, usize, usize, usize)> = Vec::new();
        for (c, [at_start, at_end]) in skipped.into_iter().enumerate() {
            if !holding.held[c] {
                continue;
            }
            let n = layout.along(c);
            let beside = |user: (usize, usize), q: usize| self.beside(layout, user, c, q);
            for q in 0..n {
                if (at_start && q < 2) || (at_end && q + 2 >= n) {
                    continue;
                }
                let (b, l, r) = (
                    layout.curve_slot(topo, c, q),
                    beside(topo.users[c][0], q),
                    beside(topo.users[c][1], q),
                );
                // At a fan's point both rows beside the curve shrink to
                // it, and the triple holds nothing.
                if b == l && b == r {
                    continue;
                }
                triples.push((c, b, l, r));
            }
        }
        let settled = self.settle(&triples, junctions, lines, corners);

        // Each chain's pin, as the ratio's numerator and denominator, and
        // whether two pins disagree.
        let mut pins: Vec<Option<(f64, f64)>> = vec![None; holding.chains];
        let mut untied = vec![false; holding.chains];
        let mut pin = |chain: usize, pinned: (f64, f64)| match pins[chain] {
            None => pins[chain] = Some(pinned),
            Some(held) => {
                if held.0 * pinned.1 != pinned.0 * held.1 {
                    untied[chain] = true;
                }
            }
        };
        let mut unpinned = Vec::new();
        for corner in corners.iter() {
            for (i, &c) in corner.curves.iter().enumerate() {
                let Some((chain, inverse)) = holding.chain[c] else {
                    continue;
                };
                match &corner.relation {
                    Some(relation) => {
                        let lambda = relation.weights[i].0;
                        pin(
                            chain,
                            if inverse {
                                (1.0, lambda)
                            } else {
                                (lambda, 1.0)
                            },
                        );
                    }
                    None => unpinned.push(chain),
                }
            }
        }
        let mut sums = vec![(0.0, 0.0); holding.chains];
        for &(c, b, l, r) in &triples {
            let Some((chain, inverse)) = holding.chain[c] else {
                continue;
            };
            if let Some(row) = junctions.get(&(b, l.min(r), l.max(r))) {
                // The row reads `a B − x L − y R`, so `λ = y / x`.
                let of = |slot: usize| row.iter().find(|t| t.0 == slot).map_or(0.0, |t| -t.1);
                let (x, y) = (of(l), of(r));
                pin(chain, if inverse { (x, y) } else { (y, x) });
                continue;
            }
            let (near, far) = (loose[b].distance(loose[l]), loose[r].distance(loose[b]));
            let (num, den) = if inverse { (far, near) } else { (near, far) };
            sums[chain].0 += num;
            sums[chain].1 += den;
        }
        for chain in unpinned {
            untied[chain] = true;
        }
        let mut ratio = Vec::with_capacity(holding.chains);
        for (chain, &(num, den)) in sums.iter().enumerate() {
            if let Some((num, den)) = pins[chain] {
                // Opposite parities somewhere are served by a ratio of one
                // alone.
                if holding.conflicted[chain] && num != den {
                    untied[chain] = true;
                }
                ratio.push(num / den);
                continue;
            }
            if holding.conflicted[chain] {
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

        // The settled classes override their chains, and a corner with a
        // curve untied is held by none of its rows: its curves are untied
        // together. A settled class that loses a chain so is unsettled, and
        // the corners are gone over again.
        let (pinned_ratio, pinned_untied) = (ratio.clone(), untied.clone());
        let mut kept = vec![true; settled.len()];
        loop {
            ratio.clone_from(&pinned_ratio);
            untied.clone_from(&pinned_untied);
            for (class, keep) in settled.iter().zip(&kept) {
                if !keep {
                    continue;
                }
                for &(chain, odd) in &class.chains {
                    ratio[chain] = if odd { 1.0 / class.ratio } else { class.ratio };
                    untied[chain] = false;
                }
            }
            loop {
                let mut changed = false;
                for corner in corners.iter() {
                    let chains: Vec<usize> = corner
                        .curves
                        .iter()
                        .filter_map(|&c| holding.chain[c].map(|(chain, _)| chain))
                        .collect();
                    if chains.iter().any(|&chain| untied[chain]) {
                        for chain in chains {
                            changed |= !untied[chain];
                            untied[chain] = true;
                        }
                    }
                }
                if !changed {
                    break;
                }
            }
            let mut lost = false;
            for (class, keep) in settled.iter().zip(kept.iter_mut()) {
                if *keep && class.chains.iter().any(|&(chain, _)| untied[chain]) {
                    *keep = false;
                    lost = true;
                }
            }
            if !lost {
                break;
            }
        }
        let mut relaxed = Keyed::new();
        for (class, _) in settled.iter().zip(&kept).filter(|(_, keep)| **keep) {
            for &(line, odd) in &class.lines {
                let mu = if odd { 1.0 / class.ratio } else { class.ratio };
                for (key, &(of, flip)) in &lines.of {
                    if of != line {
                        continue;
                    }
                    let row = &junctions[key];
                    let (mid, before, after) = if flip {
                        (row[0].0, row[2].0, row[1].0)
                    } else {
                        (row[0].0, row[1].0, row[2].0)
                    };
                    relaxed.insert(*key, vec![(mid, 1.0 + mu), (before, -1.0), (after, -mu)]);
                }
            }
        }
        let mut held = BTreeMap::new();
        let rows = triples
            .iter()
            .filter_map(|&(c, b, l, r)| {
                let (chain, inverse) = holding.chain[c]?;
                if untied[chain] || junctions.contains_key(&(b, l.min(r), l.max(r))) {
                    return None;
                }
                let lambda = if inverse {
                    1.0 / ratio[chain]
                } else {
                    ratio[chain]
                };
                held.insert(b, (l, r, lambda));
                Some(vec![(b, 1.0 + lambda), (l, -1.0), (r, -lambda)])
            })
            .collect();
        Ok(Tangency {
            rows,
            untied,
            relaxed,
            held,
        })
    }

    /// The chains and junction lines a corner pins along with something
    /// else, settled on one ratio each, as [`tangency`](Self::tangency)
    /// describes; the corners they meet more than once are rescaled to it.
    fn settle(
        &self,
        triples: &[(usize, usize, usize, usize)],
        junctions: &Keyed,
        lines: &Lines,
        corners: &mut [Related],
    ) -> Vec<Settled> {
        let holding = self.holding();
        let chains = holding.chains;
        // Chains, then lines, each node's ratio that of its class's root or
        // its inverse.
        let mut classes = topology::Classes::new(chains + lines.count);
        let mut torn = Vec::new();
        for &(c, b, l, r) in triples {
            let Some((chain, inverse)) = holding.chain[c] else {
                continue;
            };
            let key = (b, l.min(r), l.max(r));
            let (Some(row), Some(&(line, flip))) = (junctions.get(&key), lines.of.get(&key)) else {
                continue;
            };
            // The line's own direction runs from the pole before its middle
            // to the one after, as the row was first kept, or against it.
            let before = if flip { row[2].0 } else { row[1].0 };
            if !classes.join(chain, chains + line, inverse ^ (l != before)) {
                torn.push(chain);
            }
        }
        for (line, &bad) in lines.torn.iter().enumerate() {
            if bad {
                torn.push(chains + line);
            }
        }
        // Each class's corners' pins: the corner, the curve, and the ratio
        // they pin the class's root to.
        let mut pinned: BTreeMap<usize, Vec<(usize, usize, f64)>> = BTreeMap::new();
        let mut spoiled: Vec<usize> = Vec::new();
        for (k, corner) in corners.iter().enumerate() {
            for (i, &c) in corner.curves.iter().enumerate() {
                let Some((chain, inverse)) = holding.chain[c] else {
                    continue;
                };
                let (root, odd) = classes.find(chain);
                match &corner.relation {
                    Some(relation) => {
                        let lambda = relation.weights[i].0;
                        let value = if inverse ^ odd { 1.0 / lambda } else { lambda };
                        pinned.entry(root).or_default().push((k, i, value));
                    }
                    None => spoiled.push(root),
                }
            }
        }
        for node in torn {
            spoiled.push(classes.find(node).0);
        }
        for chain in 0..chains {
            if holding.conflicted[chain] {
                spoiled.push(classes.find(chain).0);
            }
        }
        let mut members: BTreeMap<usize, (Members, Members)> = BTreeMap::new();
        for node in 0..chains + lines.count {
            let (root, odd) = classes.find(node);
            let entry = members.entry(root).or_default();
            if node < chains {
                entry.0.push((node, odd));
            } else {
                entry.1.push((node - chains, odd));
            }
        }
        // The classes pinned more than once or meeting a junction, with
        // whether their corners were rescaled to them.
        let mut found: Vec<Pinned> = Vec::new();
        for (root, pins) in pinned {
            let Some((chain_members, line_members)) = members.remove(&root) else {
                continue;
            };
            if spoiled.contains(&root) || (pins.len() < 2 && line_members.is_empty()) {
                continue;
            }
            found.push(Pinned {
                pins,
                chains: chain_members,
                lines: line_members,
            });
        }
        // A class pinned by several corners takes the geometric mean of
        // their ratios, read before any corner is rescaled, moved as little
        // as [closes](close_cycles) the corners asked a ratio for every
        // curve of a cycle of their scales; each corner is then rescaled
        // once, to every ratio asked of it.
        let mut logs: Vec<Option<f64>> = vec![None; found.len()];
        let mut asked: BTreeMap<usize, Vec<Asked>> = BTreeMap::new();
        for (f, class) in found.iter().enumerate() {
            let (pins, chain_members) = (&class.pins, &class.chains);
            if pins.len() < 2 {
                continue;
            }
            let sum: f64 = pins.iter().map(|p| p.2.ln()).sum();
            logs[f] = Some(sum / precise(pins.len()));
            for &(k, i, _) in pins {
                let Some((chain, inverse)) = holding.chain[corners[k].curves[i]] else {
                    continue;
                };
                let odd = chain_members
                    .iter()
                    .find(|m| m.0 == chain)
                    .is_some_and(|m| m.1);
                asked.entry(k).or_default().push(Asked {
                    curve: i,
                    class: f,
                    odd,
                    inverse,
                });
            }
        }
        let mut failed = close_cycles(corners, &asked, &mut logs);
        let ratios: Vec<Option<f64>> = logs.iter().map(|l| l.map(f64::exp)).collect();
        for (k, asked) in asked {
            let mut pairs = Vec::with_capacity(asked.len());
            for a in &asked {
                let Some(ratio) = ratios[a.class].filter(|_| !failed[a.class]) else {
                    continue;
                };
                // The ratio in the curve's own direction, as the triples
                // read it.
                let rho = if a.odd { 1.0 / ratio } else { ratio };
                pairs.push((a.curve, if a.inverse { 1.0 / rho } else { rho }));
            }
            let rescaled = pairs.len() == asked.len()
                && corners[k]
                    .relation
                    .as_mut()
                    .is_some_and(|relation| relation.rescale(&pairs));
            if !rescaled {
                for a in asked {
                    failed[a.class] = true;
                }
            }
        }
        // A class pinned by one corner takes its ratio, as rescaled.
        let mut settled = Vec::new();
        for (f, class) in found.into_iter().enumerate() {
            let Pinned {
                pins,
                chains,
                lines,
            } = class;
            if failed[f] {
                continue;
            }
            let ratio = match ratios[f] {
                Some(ratio) => ratio,
                None => {
                    let (k, i, _) = pins[0];
                    let Some(relation) = &corners[k].relation else {
                        continue;
                    };
                    let Some((chain, inverse)) = holding.chain[corners[k].curves[i]] else {
                        continue;
                    };
                    let odd = chains.iter().find(|m| m.0 == chain).is_some_and(|m| m.1);
                    let lambda = relation.weights[i].0;
                    if inverse ^ odd { 1.0 / lambda } else { lambda }
                }
            };
            settled.push(Settled {
                ratio,
                chains,
                lines,
            });
        }
        settled
    }

    /// The slot of the pole one row into patch side `user` beside curve
    /// `c`'s pole `q`.
    /// The slots of the row of poles one in from the missing side of fan
    /// patch `p`, `k` that side, along the patch's knots in that side's
    /// direction: its ends the poles next to the fan's point of the curves
    /// either side.
    fn fan_ring(&self, layout: &Layout, (p, k): (usize, usize)) -> Vec<usize> {
        let (nu, nv) = layout.dims(p);
        let pole = |w: usize| match k {
            0 => (w, 1),
            1 => (nu - 2, w),
            2 => (w, nv - 2),
            _ => (1, w),
        };
        let along = if k % 2 == 0 { nu } else { nv };
        (0..along)
            .map(|w| layout.slot(self.topo, self.patches, p, pole(w)))
            .collect()
    }

    /// The rows holding each fan's patches to one tangent plane at its
    /// point; `None` where a curve into a fan is not held by its triple
    /// next to the point, or the poles round it are not spread over a
    /// plane in the fit `loose`.
    ///
    /// At a fan's point every patch's cross derivative along its missing
    /// side is nought, and its normal there is the limit `S_d × S_dw`, `d`
    /// across the side and `w` along it: from the point to the row of poles
    /// one in from the side, crossed with that row's own derivative. The
    /// patches share a tangent plane there where those rows, every patch's,
    /// lie in one plane through the point. That plane is an unknown
    /// direction, which rows on the poles cannot name, so the rows are laid
    /// out as they lie in the fit without the constraints, projected on
    /// the plane through the point that their consecutive poles' cross
    /// products, each row walked as its patch walks round, sum to the
    /// normal of: each pole is held to two of them, the
    /// anchors, as `P − (1 − α − β) V − α A − β B = 0`, `α` and `β` its
    /// coordinates there against theirs. The plane and the layout's scale
    /// and shear stay free through `A` and `B`. A pole next to the point on
    /// a curve into the fan is held instead by the curve's triple, which
    /// puts it on the line through the poles beside it.
    fn fan_rows(
        &self,
        layout: &Layout,
        loose: &[Point],
        triples: &BTreeMap<usize, (usize, usize, f64)>,
    ) -> Option<Vec<Row>> {
        let topo = self.topo;
        let mut rows = Vec::new();
        for fan in &topo.fans {
            let v = fan.vertex;
            for &c in &fan.curves {
                let n = layout.along(c);
                let next = if topo.curve_ends[c].0 == v { 1 } else { n - 2 };
                triples.get(&layout.curve_slot(topo, c, next))?;
            }
            let rings: Vec<Vec<usize>> = fan
                .patches
                .iter()
                .map(|&side| self.fan_ring(layout, side))
                .collect();
            let e = |slot: usize| loose[slot] - loose[v];
            // Each row walked as its patch walks its missing side, which
            // turns every patch's one way round the point.
            let mut normal = Vector::ZERO;
            for (ring, &(_, k)) in rings.iter().zip(&fan.patches) {
                for w in ring.windows(2) {
                    let turn = e(w[0]).cross(e(w[1]));
                    normal += if k < 2 { turn } else { -turn };
                }
            }
            let inner: Vec<usize> = rings
                .iter()
                .flat_map(|ring| ring[1..ring.len() - 1].iter().copied())
                .collect();
            let first = e(inner[0]);
            let x = first - normal * (first.dot(normal) / normal.dot(normal));
            let x = x / x.magnitude();
            let y = (normal / normal.magnitude()).cross(x);
            let flat = |slot: usize| {
                let d = e(slot);
                (d.dot(x), d.dot(y))
            };
            let flats: Vec<(f64, f64)> = inner.iter().map(|&s| flat(s)).collect();
            let (ax, ay) = flats[0];
            let cross = |(bx, by): (f64, f64)| ax.mul_add(by, -(ay * bx));
            let b = (1..flats.len()).max_by(|&i, &j| {
                let sine = |k: usize| {
                    let (bx, by) = flats[k];
                    cross(flats[k]).abs() / bx.hypot(by)
                };
                sine(i).total_cmp(&sine(j))
            })?;
            let det = cross(flats[b]);
            if !(det.abs() > 0.0 && det.is_finite()) {
                return None;
            }
            let (bx, by) = flats[b];
            for (i, &(px, py)) in flats.iter().enumerate() {
                if i == 0 || i == b {
                    continue;
                }
                let alpha = px.mul_add(by, -(py * bx)) / det;
                let beta = ax.mul_add(py, -(ay * px)) / det;
                if !(alpha.is_finite() && beta.is_finite()) {
                    return None;
                }
                rows.push(vec![
                    (inner[i], 1.0),
                    (inner[0], -alpha),
                    (inner[b], -beta),
                    (v, alpha + beta - 1.0),
                ]);
            }
        }
        Some(rows)
    }

    /// Whether every fan patch's normal at its fan's point, held to one
    /// plane by [`fan_rows`](Self::fan_rows), is defined all along its
    /// missing side and points one way there: the limit normal is
    /// `(R − V) × R'` along the row `R` of poles one in from the side, a
    /// spline whose Bernstein coefficients, span by span, all on one side
    /// of the plane bound it off nought. Across each curve into the fan
    /// the triples hold the two patches' normals alike, and so up to the
    /// point, so every patch's normal there is the same.
    fn fans_regular(
        &self,
        layout: &Layout,
        patch_knots: &[[KnotVector; 2]],
        poles: &[Point],
    ) -> OgeomResult<bool> {
        for fan in &self.topo.fans {
            let at = poles[fan.vertex];
            for &(p, k) in &fan.patches {
                let ring: Vec<Vector> = self
                    .fan_ring(layout, (p, k))
                    .iter()
                    .map(|&s| poles[s] - at)
                    .collect();
                let segments = ogeom_math::bspline::to_bezier_segments(
                    &patch_knots[p][k % 2],
                    &ring,
                    self.tol,
                )?;
                // `(R − V) × R'` over each span as a quintic, from the
                // cubic's poles `b` and its derivative's, `Δb`: each
                // coefficient the sum of `C(3,i) C(2,j) b_i × Δb_j` over
                // `i + j = k`, over `C(5,k)`, dropping the derivative's
                // positive scale.
                const CUBIC: [f64; 4] = [1.0, 3.0, 3.0, 1.0];
                const QUADRATIC: [f64; 3] = [1.0, 2.0, 1.0];
                const QUINTIC: [f64; 6] = [1.0, 5.0, 10.0, 10.0, 5.0, 1.0];
                let mut coefficients = Vec::new();
                for (_, b) in &segments {
                    let d: Vec<Vector> = b.windows(2).map(|w| w[1] - w[0]).collect();
                    for (n, &q) in QUINTIC.iter().enumerate() {
                        let mut sum = Vector::ZERO;
                        for (i, &ci) in CUBIC.iter().enumerate() {
                            if let Some(j) = n.checked_sub(i)
                                && j < 3
                            {
                                sum += b[i].cross(d[j]) * (ci * QUADRATIC[j]);
                            }
                        }
                        coefficients.push(sum / q);
                    }
                }
                let normal = coefficients.iter().fold(Vector::ZERO, |a, &c| a + c);
                if !coefficients.iter().all(|c| c.dot(normal) > 0.0) {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

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
    /// and `λ > 0`. At the corner the condition is the corner's
    /// [relation](Self::relate). With `λ` constant and `γ` linear along the
    /// curve the left side is a cubic spline with double knots, `C¹` where
    /// the curve's knots are, so it is nought throughout where it is nought
    /// at the corner and at two points inside each span and the far end,
    /// which is as many points as it has coefficients: an exact, linear
    /// condition on the poles, though no pole-by-pole one.
    fn corner_rows(
        &self,
        layout: &Layout,
        loose: &[Point],
        curve_knots: &[KnotVector],
    ) -> Vec<Vec<(usize, f64)>> {
        let mut rows = Vec::new();
        for (v, curves) in &self.topo.odd_corners {
            let Some(relation) = self.relate(layout, loose, *v, curves) else {
                continue;
            };
            rows.extend(relation.rows());
            for (i, &c) in curves.iter().enumerate() {
                rows.extend(self.along_rows(
                    layout,
                    loose,
                    c,
                    *v,
                    relation.weights[i],
                    &curve_knots[c],
                ));
            }
        }
        rows
    }

    /// The rows holding tangent every curve into a corner whose gluing is
    /// spread along its curves; `None` where a corner has no relation or the
    /// corners a curve joins cannot be scaled to one ratio for it.
    ///
    /// Each curve is held with `λ` constant and `γ` linear along it, as
    /// [`corner_rows`](Self::corner_rows) holds a curve, but every knot
    /// inside its strip is doubled, so its patches' cross derivatives are
    /// only `C¹` across those knots, and the curve is held `C²` across
    /// each, so `γ T` is `C¹` there too: the left side is a `C¹` cubic
    /// spline on the doubled knots, which [`spread_along_rows`](Self::spread_along_rows)
    /// holds nought. With single knots `γ T` would jump in its second
    /// derivative by `γ` times the jump of `T''` at each knot, binding the
    /// curve to a single cubic; with them doubled it is a `C²` cubic spline,
    /// refined as any. A curve joining two such corners is asked a ratio by
    /// each: they are scaled to the geometric mean of the two, moved as
    /// little as [closes](close_cycles) every cycle of a corner's scales,
    /// and `γ` runs from the one corner's `g` to the other's. The other
    /// curves along the doubled strips carry no gluing, and are held `C²`
    /// across the doubled knots as they are where the knots are single.
    fn spread_rows(
        &self,
        layout: &Layout,
        loose: &[Point],
        curve_knots: &[KnotVector],
    ) -> Option<Vec<Row>> {
        let topo = self.topo;
        let mut corners = Vec::with_capacity(topo.spread_corners.len());
        for (v, curves) in &topo.spread_corners {
            let relation = self.relate(layout, loose, *v, curves)?;
            corners.push(Related {
                vertex: *v,
                curves: curves.clone(),
                relation: Some(relation),
            });
        }
        // Each curve's place round the corners it runs into.
        let mut ends: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
        for (k, corner) in corners.iter().enumerate() {
            for (i, &c) in corner.curves.iter().enumerate() {
                ends.entry(c).or_default().push((k, i));
            }
        }
        let mut logs = Vec::new();
        let mut asked: BTreeMap<usize, Vec<Asked>> = BTreeMap::new();
        for at in ends.values().filter(|at| at.len() == 2) {
            let class = logs.len();
            let mut sum = 0.0;
            for &(k, i) in at {
                sum += corners[k].relation.as_ref()?.weights[i].0.ln();
                asked.entry(k).or_default().push(Asked {
                    curve: i,
                    class,
                    odd: false,
                    inverse: false,
                });
            }
            logs.push(Some(sum / 2.0));
        }
        if close_cycles(&corners, &asked, &mut logs).contains(&true) {
            return None;
        }
        for (&k, asked) in &asked {
            let pairs = asked
                .iter()
                .map(|a| Some((a.curve, logs[a.class]?.exp())))
                .collect::<Option<Vec<_>>>()?;
            if !corners[k].relation.as_mut()?.rescale(&pairs) {
                return None;
            }
        }
        let mut rows = Vec::new();
        for corner in &corners {
            rows.extend(corner.relation.as_ref()?.rows());
        }
        for (&c, at) in &ends {
            let (k, i) = at[0];
            let relation = corners[k].relation.as_ref()?;
            let other = match at.get(1) {
                Some(&(m, j)) => Some(corners[m].relation.as_ref()?.weights[j].1),
                None => None,
            };
            rows.extend(self.spread_along_rows(
                layout,
                loose,
                c,
                corners[k].vertex,
                relation.weights[i],
                &curve_knots[c],
                other,
            ));
            rows.extend(smooth_rows(layout, topo, c, &curve_knots[c]));
        }
        let doubled = self.doubled();
        for (c, strips) in topo.curve_strip.iter().enumerate() {
            if !ends.contains_key(&c) && strips.iter().any(|&s| doubled[s].every) {
                rows.extend(smooth_rows(layout, topo, c, &curve_knots[c]));
            }
        }
        Some(rows)
    }

    /// How the curves into corner `v` of other than four patches relate
    /// there, read off the fit `loose`; `None` where they do not close
    /// round it or a curve's ratio is not positive.
    ///
    /// At the corner `D₁` and `D₂` run along the neighbouring curves, so
    /// there `D₁ + λ D₂ + γ T = 0` is a relation among the first poles of
    /// the curves along each of them: they lie in one plane through the
    /// corner. The relations come from `loose`: its first poles projected
    /// on the plane their consecutive cross products face, each curve's `λ`
    /// and `γ` the unique relation among its own first pole and the two
    /// beside it there. Those of `n` curves round a corner have `n − 2`
    /// independent ones, so `n − 2` consecutive ones are held and the other
    /// two follow.
    ///
    /// The cross derivatives at the corner differentiated along each curve
    /// bring in each patch's twist pole, once from either side, in a cycle.
    /// Each curve's ratio is that of the cross products of its first pole
    /// with those of the curves either side, so the ratios round the corner
    /// multiply to one: the cycle is never singular at an odd number of
    /// curves and always singular at an even one (the vertex enclosure
    /// condition), where it leaves one condition on the curves' own poles.
    /// Corners held whole are odd only; a corner held by its spans or with
    /// its gluing spread may be even, the condition then falling on the
    /// poles of the curves' spans next to it, or on their second poles.
    fn relate(
        &self,
        layout: &Layout,
        loose: &[Point],
        v: usize,
        curves: &[usize],
    ) -> Option<Relation> {
        let topo = self.topo;
        // Each curve's pole next to the corner, and the poles beside it
        // there in its first and second patch: those of the curves either
        // side.
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
            let next = (0..curves.len())
                .find(|&i| !order.contains(&i) && (near[i].1 == last.0 || near[i].1 == last.2))?;
            order.push(next);
        }
        let (first, last) = (near[order[0]], near[order[order.len() - 1]]);
        if !(first.0 == last.1 || first.2 == last.1) {
            return None;
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
        let flat: Vec<(f64, f64)> = near.iter().map(|&(_, own, _)| flat(own)).collect();
        let curve_at = |slot: usize| near.iter().position(|n| n.1 == slot);
        let sides = near
            .iter()
            .map(|&(left, _, right)| Some((curve_at(left)?, curve_at(right)?)))
            .collect::<Option<Vec<_>>>()?;
        let weights = plane_weights(&flat, &sides)?;
        Some(Relation {
            vertex: v,
            near,
            order,
            weights,
            flat,
            sides,
        })
    }

    /// The rows holding `D₁ + λ D₂ + γ T` nought along the span of curve
    /// `c` next to its end at corner `v`, which the curve's other poles'
    /// triples leave free, and that span a quadratic; `weights` are `λ`
    /// and the relation's `g` at `v`.
    ///
    /// The knot past that span is doubled, so over it the two basis
    /// functions nearest the corner are the cubic Bernstein polynomials
    /// `(1 − x)³` and `3 x (1 − x)²`, `x` running across the span from the
    /// corner, and the patches are only `C¹` across its far end. With `γ`
    /// quadratic over the span, `γ(0) (1 − x)²`, and nought past it, and `T`
    /// linear over it, `γ T` is `γ(0) T(0)` times the first and `γ(0) T(1)
    /// / 3` times the second; the curve's triples past those two poles hold
    /// the rest nought. So the triple of the pole next to the corner is
    /// `−γ(0) T(0)`, which is the corner's relation, and that of the pole
    /// after it `−γ(0) T(1) / 3`, `T(1)` the curve's derivative at the far
    /// end of the span: the left side is then nought throughout, `C¹` across
    /// the span's far end on either side, as `γ` has a double root there.
    /// Past the span `γ` is nought, so the curve is held as a curve between
    /// corners of four patches is, and its far end may be one.
    fn span_rows(
        &self,
        layout: &Layout,
        c: usize,
        v: usize,
        (lambda, g): (f64, f64),
        knots: &KnotVector,
    ) -> Vec<Vec<(usize, f64)>> {
        let topo = self.topo;
        let t = knots.knots();
        let n = layout.along(c);
        let starts = topo.curve_ends[c].0 == v;
        // The span next to the corner, its far end, its width, and the pole
        // after the corner's next one.
        let (span, far, width, next) = if starts {
            (DEGREE, t[DEGREE + 1], t[DEGREE + 1] - t[DEGREE], 1)
        } else {
            (n - 1, t[n - 1], t[n] - t[n - 1], n - 2)
        };
        // `γ(0) T` in the poles' differences: `T` at the corner is `3 / h`
        // times the first pole's difference from it, and the other way
        // about at the curve's end.
        let near = if starts {
            g * width / precise(DEGREE)
        } else {
            -g * width / precise(DEGREE)
        };
        let own = |i: usize| layout.curve_slot(topo, c, span - DEGREE + i);
        let [first, second] = [topo.users[c][0], topo.users[c][1]];
        let mut row: BTreeMap<usize, f64> = BTreeMap::new();
        *row.entry(self.beside(layout, first, c, next))
            .or_insert(0.0) += 1.0;
        *row.entry(self.beside(layout, second, c, next))
            .or_insert(0.0) += lambda;
        *row.entry(layout.curve_slot(topo, c, next)).or_insert(0.0) -= 1.0 + lambda;
        let slope = knots.basis_derivatives(span, far, 1);
        for (i, rate) in slope[1].iter().enumerate() {
            *row.entry(own(i)).or_insert(0.0) += near / precise(DEGREE) * rate;
        }
        // The span's third derivative, nought, scaled to its width.
        let middle = f64::midpoint(t[span], t[span + 1]);
        let jerk = knots.basis_derivatives(span, middle, DEGREE);
        let cube = width * width * width;
        vec![
            row.into_iter().filter(|e| e.1 != 0.0).collect(),
            jerk[DEGREE]
                .iter()
                .enumerate()
                .map(|(i, d)| (own(i), d * cube))
                .collect(),
        ]
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
        let starts = self.topo.curve_ends[c].0 == v;
        let near = end_gamma(layout, c, g, knots, starts);
        let terms = self.along_terms(
            layout,
            c,
            lambda,
            knots,
            starts,
            &along_points(knots, starts),
        );
        let far = least_far(loose, &terms, near);
        terms
            .iter()
            .map(|(cross, slope, share)| {
                gamma_row(cross, slope, (1.0 - share).mul_add(near, share * far))
            })
            .collect()
    }

    /// At each of `at` along curve `c`, the slots and their weights in
    /// `D₁ + λ D₂` and in `T`, and the share of the way from the corner
    /// the rows hold it from, at its start where `starts`.
    fn along_terms(
        &self,
        layout: &Layout,
        c: usize,
        lambda: f64,
        knots: &KnotVector,
        starts: bool,
        at: &[f64],
    ) -> Vec<(Weighted, Weighted, f64)> {
        let topo = self.topo;
        let [first, second] = [topo.users[c][0], topo.users[c][1]];
        at.iter()
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
            .collect()
    }

    /// The rows holding `D₁ + λ D₂ + γ T` nought along curve `c` from its
    /// end at corner `v`, every knot inside it doubled and the curve held
    /// `C²` across each by [`smooth_rows`]: the left side is then a `C¹`
    /// cubic spline on the doubled knots, nought where it and its
    /// derivative are nought at every knot. `weights` are `λ` and the
    /// relation's `g` at `v`, and `γ` runs linearly from there to the
    /// other end, its value there as [`along_rows`](Self::along_rows) reads
    /// it, or, where the other end is a corner whose relation gives the
    /// curve `g` there as `other`, that corner's. The value at a corner is
    /// its relation's to hold, so only the derivative is held there.
    ///
    /// Held at the knots, each row bears on the poles of the spans either
    /// side alone, however many spans the curve has; collocated inside the
    /// spans as [`along_rows`](Self::along_rows) holds them, the rows would
    /// tie each span to the last, and the far end's row would follow from
    /// the others ever more nearly as the spans multiply.
    #[allow(
        clippy::too_many_arguments,
        reason = "the curve, its corner's weights, and the other end's"
    )]
    fn spread_along_rows(
        &self,
        layout: &Layout,
        loose: &[Point],
        c: usize,
        v: usize,
        (lambda, g): (f64, f64),
        knots: &KnotVector,
        other: Option<f64>,
    ) -> Vec<Row> {
        let topo = self.topo;
        let starts = topo.curve_ends[c].0 == v;
        let near = end_gamma(layout, c, g, knots, starts);
        let far = match other {
            // The other corner's `g` at the curve's other end, read as
            // `near` reads its own.
            Some(other) => end_gamma(layout, c, other, knots, !starts),
            None => {
                let points = along_points(knots, starts);
                least_far(
                    loose,
                    &self.along_terms(layout, c, lambda, knots, starts, &points),
                    near,
                )
            }
        };
        // `γ` and its derivative along the curve's parameter.
        let gamma = |share: f64| (1.0 - share).mul_add(near, share * far);
        let rate = if starts { far - near } else { near - far };
        let [first, second] = [topo.users[c][0], topo.users[c][1]];
        let mut rows = Vec::new();
        let distinct = knots.distinct();
        let last = distinct.len() - 1;
        for (k, &(x, _)) in distinct.iter().enumerate() {
            let share = if starts { x } else { 1.0 - x };
            let (at_v, at_other) = if starts {
                (k == 0, k == last)
            } else {
                (k == last, k == 0)
            };
            if !(at_v || (at_other && other.is_some())) {
                let terms = self.along_terms(layout, c, lambda, knots, starts, &[x]);
                rows.extend(
                    terms
                        .iter()
                        .map(|(cross, slope, share)| gamma_row(cross, slope, gamma(*share))),
                );
            }
            // The derivative, scaled to the width of the span it is read
            // in, the one after the knot but at the curve's end.
            let span = knots.span_unchecked(x);
            let width = if k < last {
                distinct[k + 1].0 - x
            } else {
                x - distinct[k - 1].0
            };
            let basis = knots.basis_derivatives(span, x, 2);
            let mut row: BTreeMap<usize, f64> = BTreeMap::new();
            for i in 0..=DEGREE {
                let q = span - DEGREE + i;
                let (slope, bend) = (basis[1][i], basis[2][i]);
                let own = layout.curve_slot(topo, c, q);
                *row.entry(self.beside(layout, first, c, q)).or_insert(0.0) += width * slope;
                *row.entry(self.beside(layout, second, c, q)).or_insert(0.0) +=
                    width * lambda * slope;
                *row.entry(own).or_insert(0.0) +=
                    width * (rate.mul_add(slope, gamma(share) * bend) - (1.0 + lambda) * slope);
            }
            rows.push(row.into_iter().filter(|e| e.1 != 0.0).collect());
        }
        rows
    }

    /// Finish a try that missed at `spans`, from its fit `start` at the
    /// seed's parameters: while not every point is within the tolerance,
    /// the points are re-taken at their nearest points on the fit and the
    /// patches fitted again at the same spans. The fit and the points as
    /// re-taken where it comes within the tolerance; `None` where it stops
    /// first.
    ///
    /// The noise's variance is read once, off `start`'s median distance as
    /// refinement reads it, and held: re-read after a correction, whose
    /// distances are all along the normal, it would shrink and let the fit
    /// follow the noise. A fit is kept only where it holds every corner the
    /// try names, takes off more of the squared distances than twice that
    /// variance (what Mallows' criterion charges a free pole), and the
    /// worst miss that is not a stray grows no worse. The corrections stop
    /// where the squared distances are down to the noise left over the free
    /// controls, `σ̂₀²(m − free)`, and after [`FINISH_STEPS`].
    fn finish(
        &self,
        spans: &[usize],
        start: Round,
        tolerance: f64,
    ) -> OgeomResult<Option<(Round, Taken)>> {
        let judged = (DEGREE + 1) * (DEGREE + 1);
        let variance = (median(&start.distances) / NORMAL_QUARTILE).powi(2);
        let floor = variance * precise(start.distances.len().saturating_sub(start.free));
        let worst = |fit: &Fit, round: &Round| {
            fit.missed(round, tolerance, judged)
                .0
                .first()
                .map_or(0.0, |w| w.0)
        };
        let mut taken = self.taken.clone();
        let mut current = start;
        for step in 0..=FINISH_STEPS {
            if current.distances.iter().all(|d| *d <= tolerance) {
                return Ok(Some((current, taken)));
            }
            if step == FINISH_STEPS || variance <= 0.0 || current.square() <= floor {
                return Ok(None);
            }
            let fit = Fit {
                taken: &taken,
                ..*self
            };
            let again = fit.retake(&current)?;
            let refit = Fit {
                taken: &again,
                ..*self
            };
            let next = refit.round(spans, &again.at)?;
            let kept = next.held
                && current.square() - next.square() > 2.0 * variance
                && worst(&refit, &next) <= worst(&fit, &current);
            if !kept {
                return Ok(None);
            }
            (taken, current) = (again, next);
        }
        Ok(None)
    }

    /// The points re-taken at their nearest points on `round`'s patches:
    /// each by its own patch, or, where that foot lies on a side its patch
    /// shares, by the patch across the side where that is nearer. Such a
    /// point lies past its patch's side, on the other's side of the curve,
    /// and the nearest patch takes it, as the seeds took every point.
    fn retake(&self, round: &Round) -> OgeomResult<Taken> {
        let topo = self.topo;
        let geometry: Vec<SurfaceGeometry> = round
            .surfaces
            .iter()
            .map(|s| SurfaceGeometry::BSpline(s.clone()))
            .collect();
        let mut taken = self.taken.clone();
        taken.at.clone_from(&round.feet);
        for (i, &p) in self.taken.patch.iter().enumerate() {
            let (u, v) = round.feet[i];
            let mut best: Option<(usize, (f64, f64), f64)> = None;
            for (k, x, on) in [
                (0, u, v <= 0.0),
                (1, v, u >= 1.0),
                (2, u, v >= 1.0),
                (3, v, u <= 0.0),
            ] {
                if !on {
                    continue;
                }
                let parts = side_parts(topo, self.patches, (p, k));
                let piece = piece_at(k, x, &parts, &round.patch_bounds[p][k % 2]);
                let Some(side) = self.patches[p].sides[k].get(piece) else {
                    continue;
                };
                for &(q, _) in &topo.users[side.curve] {
                    if q == p {
                        continue;
                    }
                    let foot =
                        ogeom_geom::project_on_surface(&geometry[q], taken.point[i], 0, self.tol)?;
                    if foot.distance < best.map_or(round.distances[i], |b| b.2) {
                        best = Some((q, foot.parameters, foot.distance));
                    }
                }
            }
            if let Some((q, at, _)) = best {
                taken.patch[i] = q;
                taken.at[i] = at;
            }
        }
        for mine in &mut taken.of_patch {
            mine.clear();
        }
        for (i, &p) in taken.patch.iter().enumerate() {
            taken.of_patch[p].push(i);
        }
        Ok(taken)
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
