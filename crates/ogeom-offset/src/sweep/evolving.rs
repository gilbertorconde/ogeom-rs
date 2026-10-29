//! Pipe shells whose section scales and twists along the spine.
//!
//! The section is the profile placed in the law's frame at each point of
//! the spine, as for [`make_pipe_shell_with`](super::make_pipe_shell_with),
//! and before that scaled about the spine and turned about it, each by its
//! own law of the share of the spine's length run there: linear, a
//! monotone cubic through given values, or a function the caller gives.
//!
//! Scaled without a twist down a straight spine, every profile point runs
//! along a straight line, at a constant speed, to its scaled image at the
//! end: the solid is the ruled loft between the profile and that image,
//! exact (a pyramid's frustum over a polygon, a cone's over a circle, the
//! pyramid or cone itself where the scale reaches zero). Everywhere else
//! the walls are skinned through the placed sections and held to the
//! tolerance between them against the law itself.

use std::sync::Arc;

use ogeom_algo::{Built, History, curve_length, find_plane};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve;
use ogeom_math::{Axis, Direction, Point, Transform};
use ogeom_topo::{Filter, Model, Shape, ShapeType, VertexData, explore};

use super::corners::law_cornered;
use super::placement::{PipeContact, profile_in_contact};
use super::{
    PipeCorners, PipeLaw, SpineStation, WalkCurve, closed_pipe_shell, closes_on_a_corner,
    densified, evenly_by, has_sharp_corner, law_loft, make_loft, make_pipe_shell_with, moved_ring,
    placed_path, profile_square_at_start, shell_stations, spine_closes, spine_curve_of,
    straight_legs,
};

/// How a pipe's section scales and twists along its spine, each in
/// proportion to the length run along the spine: none at its start, all
/// of it at its end.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PipeEvolution {
    /// The section's scale at the spine's end, about the spine; it is one
    /// at the start. Zero closes the pipe on the spine's end point.
    pub scale: f64,
    /// The angle in radians the section has turned through about the
    /// spine by its end, counter-clockwise about the direction of travel.
    pub twist: f64,
}

/// How a quantity runs along a pipe's spine, against the share of the
/// spine's length run: from its value at the start (one for a scale,
/// naught for a twist) at a share of naught.
#[derive(Debug, Clone, PartialEq)]
pub enum SpineLaw {
    /// Linearly to this value at the spine's end.
    Linear(f64),
    /// Through these `(share, value)` points, the shares rising strictly
    /// and the last of them one, the spine's end. Between two neighbouring
    /// points the law is the cubic that keeps between their values, and it
    /// runs on through each point with one slope (the monotone cubic
    /// Hermite interpolant): a value the points never leave, the law never
    /// leaves either.
    Through(Vec<(f64, f64)>),
    /// The value the function gives at each share, from naught to one. It
    /// gives the start value at naught exactly. The sweep reads it at its
    /// stations and wherever it checks its skin between them, and holds
    /// the skin to the tolerance against what it reads there; it is to be
    /// continuous over the spine, and a scale law above naught all along.
    Function(SpineFunction),
}

/// A law of the share of a pipe's spine run, given as a function from a
/// share in `[0, 1]` to the law's value there.
#[derive(Clone)]
pub struct SpineFunction(Arc<dyn Fn(f64) -> f64 + Send + Sync>);

impl SpineFunction {
    /// The law `f`.
    pub fn new(f: impl Fn(f64) -> f64 + Send + Sync + 'static) -> Self {
        Self(Arc::new(f))
    }

    /// The law's value at share `lambda`.
    pub fn at(&self, lambda: f64) -> f64 {
        (self.0)(lambda)
    }
}

impl core::fmt::Debug for SpineFunction {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str("SpineFunction(..)")
    }
}

/// Two laws are one where they are the same function.
impl PartialEq for SpineFunction {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

/// How a pipe's section scales and twists along its spine, each by its
/// own law of the share of the spine's length run, and how its profile is
/// brought onto the spine first.
#[derive(Debug, Clone, PartialEq)]
pub struct PipeLaws {
    /// The section's scale about the spine; one at the start.
    pub scale: SpineLaw,
    /// The angle in radians the section has turned through about the
    /// spine, counter-clockwise about the direction of travel; naught at
    /// the start.
    pub twist: SpineLaw,
    /// How the profile is moved onto the spine's start before it is
    /// swept; as placed by default.
    pub contact: PipeContact,
}

/// The section neither scaled nor twisted, its profile swept as placed.
impl Default for PipeLaws {
    fn default() -> Self {
        Self {
            scale: SpineLaw::Linear(1.0),
            twist: SpineLaw::Linear(0.0),
            contact: PipeContact::AsPlaced,
        }
    }
}

impl From<PipeEvolution> for PipeLaws {
    fn from(evolution: PipeEvolution) -> Self {
        Self {
            scale: SpineLaw::Linear(evolution.scale),
            twist: SpineLaw::Linear(evolution.twist),
            contact: PipeContact::AsPlaced,
        }
    }
}

/// A [`SpineLaw`] read for evaluation: its values at rising shares from
/// naught to one, and the slope the law runs through each with.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Ramp {
    shares: Vec<f64>,
    values: Vec<f64>,
    slopes: Vec<f64>,
    /// Whether the law runs straight from its start value to its end one.
    linear: bool,
    /// A law given as a function, read in place of the cubics.
    function: Option<FunctionRead>,
}

/// A [`SpineFunction`] read over a stretch of the spine: at a share `λ` of
/// the stretch it gives `gain f(from + (to - from) λ) + offset + tilt λ`.
#[derive(Debug, Clone, PartialEq)]
struct FunctionRead {
    law: SpineFunction,
    from: f64,
    to: f64,
    gain: f64,
    offset: f64,
    tilt: f64,
}

impl FunctionRead {
    fn at(&self, lambda: f64) -> f64 {
        self.gain * self.law.at(self.from + (self.to - self.from) * lambda)
            + self.offset
            + self.tilt * lambda
    }
}

/// How many even shares a law given as a function is read at to bound
/// from below how far it runs, which seeds the stations a skin starts
/// from (the skin then refines itself between them), and to refuse a
/// scale it shows at or below naught.
const FUNCTION_READS: u32 = 64;

impl Ramp {
    /// The law that holds `value` all along.
    pub(super) fn level(value: f64) -> Self {
        Self {
            shares: vec![0.0, 1.0],
            values: vec![value, value],
            slopes: vec![0.0, 0.0],
            linear: true,
            function: None,
        }
    }

    /// The law `f`, which gives `start` at the spine's start; `what` names
    /// it in a refusal.
    fn of_function(f: &SpineFunction, start: f64, what: &str) -> OgeomResult<Self> {
        let (first, last) = (f.at(0.0), f.at(1.0));
        if first != start {
            ogeom_bail!(
                Construction,
                "a {what} law given as a function gives {start} at the spine's start, \
                 not {first}"
            );
        }
        if !last.is_finite() {
            ogeom_bail!(Construction, "a {what} law takes finite values");
        }
        Ok(Self {
            shares: vec![0.0, 1.0],
            values: vec![first, last],
            slopes: vec![0.0, 0.0],
            linear: false,
            function: Some(FunctionRead {
                law: f.clone(),
                from: 0.0,
                to: 1.0,
                gain: 1.0,
                offset: 0.0,
                tilt: 0.0,
            }),
        })
    }

    /// `law` read from `start` at the spine's start; `what` names it in a
    /// refusal.
    fn of(law: &SpineLaw, start: f64, what: &str) -> OgeomResult<Self> {
        let mut shares = vec![0.0];
        let mut values = vec![start];
        match law {
            SpineLaw::Linear(end) => {
                shares.push(1.0);
                values.push(*end);
            }
            SpineLaw::Through(points) => {
                for &(share, value) in points {
                    shares.push(share);
                    values.push(value);
                }
            }
            SpineLaw::Function(f) => return Self::of_function(f, start, what),
        }
        if values.iter().any(|v| !v.is_finite()) {
            ogeom_bail!(Construction, "a {what} law takes finite values");
        }
        if shares
            .windows(2)
            .any(|w| w[1].partial_cmp(&w[0]) != Some(core::cmp::Ordering::Greater))
            || shares[shares.len() - 1] != 1.0
        {
            ogeom_bail!(
                Construction,
                "a {what} law's shares of the spine rise strictly from naught to one, \
                 the last of them one"
            );
        }
        // Fritsch and Carlson's slopes: each span's secant, the mean of two
        // where they agree in sign and none where they do not, each pair
        // then drawn in onto the circle of radius three about the secant,
        // inside which the span's cubic keeps between its ends.
        let n = shares.len();
        let secants: Vec<f64> = (0..n - 1)
            .map(|k| (values[k + 1] - values[k]) / (shares[k + 1] - shares[k]))
            .collect();
        let mut slopes = vec![0.0; n];
        slopes[0] = secants[0];
        slopes[n - 1] = secants[n - 2];
        for k in 1..n - 1 {
            if secants[k - 1] * secants[k] > 0.0 {
                slopes[k] = f64::midpoint(secants[k - 1], secants[k]);
            }
        }
        for k in 0..n - 1 {
            if secants[k] == 0.0 {
                slopes[k] = 0.0;
                slopes[k + 1] = 0.0;
                continue;
            }
            let (a, b) = (slopes[k] / secants[k], slopes[k + 1] / secants[k]);
            let reach = a.hypot(b);
            if reach > 3.0 {
                slopes[k] = 3.0 * a / reach * secants[k];
                slopes[k + 1] = 3.0 * b / reach * secants[k];
            }
        }
        Ok(Self {
            linear: matches!(law, SpineLaw::Linear(_)),
            shares,
            values,
            slopes,
            function: None,
        })
    }

    /// The span of the law `lambda` falls in, and how far along it.
    fn span(&self, lambda: f64) -> (usize, f64) {
        let k = self.shares[1..self.shares.len() - 1].partition_point(|s| *s < lambda);
        let t = (lambda - self.shares[k]) / (self.shares[k + 1] - self.shares[k]);
        (k, t)
    }

    /// The law's value a share `lambda` of the way along.
    pub(super) fn at(&self, lambda: f64) -> f64 {
        let lambda = lambda.clamp(0.0, 1.0);
        if let Some(f) = &self.function {
            return f.at(lambda);
        }
        if self.linear {
            return self.values[0] + (self.values[1] - self.values[0]) * lambda;
        }
        let (k, t) = self.span(lambda);
        if self.values[k] == self.values[k + 1]
            && self.slopes[k] == 0.0
            && self.slopes[k + 1] == 0.0
        {
            return self.values[k];
        }
        let h = self.shares[k + 1] - self.shares[k];
        let (t2, t3) = (t * t, t * t * t);
        (2.0 * t3 - 3.0 * t2 + 1.0) * self.values[k]
            + (t3 - 2.0 * t2 + t) * h * self.slopes[k]
            + (3.0 * t2 - 2.0 * t3) * self.values[k + 1]
            + (t3 - t2) * h * self.slopes[k + 1]
    }

    /// The law's rate of change with the share a share `lambda` of the way
    /// along.
    fn slope_at(&self, lambda: f64) -> f64 {
        let lambda = lambda.clamp(0.0, 1.0);
        if self.linear {
            return self.values[1] - self.values[0];
        }
        let (k, t) = self.span(lambda);
        let h = self.shares[k + 1] - self.shares[k];
        let t2 = t * t;
        ((6.0 * t2 - 6.0 * t) * (self.values[k] - self.values[k + 1])) / h
            + (3.0 * t2 - 4.0 * t + 1.0) * self.slopes[k]
            + (3.0 * t2 - 2.0 * t) * self.slopes[k + 1]
    }

    /// The law over the stretch from share `from` to share `to`, read
    /// against the share of that stretch: the same cubics, cut where the
    /// stretch ends.
    fn window(&self, from: f64, to: f64) -> Self {
        let width = to - from;
        if let Some(f) = &self.function {
            let span = f.to - f.from;
            let read = FunctionRead {
                from: f.from + span * from,
                to: f.from + span * to,
                offset: f.offset + f.tilt * from,
                tilt: f.tilt * width,
                ..f.clone()
            };
            return Self {
                shares: vec![0.0, 1.0],
                values: vec![read.at(0.0), read.at(1.0)],
                slopes: vec![0.0, 0.0],
                linear: false,
                function: Some(read),
            };
        }
        if self.linear {
            return Self {
                shares: vec![0.0, 1.0],
                values: vec![self.at(from), self.at(to)],
                slopes: vec![0.0, 0.0],
                linear: true,
                function: None,
            };
        }
        let mut shares = vec![0.0];
        let mut values = vec![self.at(from)];
        let mut slopes = vec![self.slope_at(from) * width];
        for k in 0..self.shares.len() {
            if self.shares[k] > from && self.shares[k] < to {
                shares.push((self.shares[k] - from) / width);
                values.push(self.values[k]);
                slopes.push(self.slopes[k] * width);
            }
        }
        shares.push(1.0);
        values.push(self.at(to));
        slopes.push(self.slope_at(to) * width);
        Self {
            shares,
            values,
            slopes,
            linear: false,
            function: None,
        }
    }

    /// The law with every value taken through `f`, an affine map whose
    /// rate is `rate`.
    fn mapped(mut self, f: impl Fn(f64) -> f64, rate: f64) -> Self {
        if let Some(read) = &mut self.function {
            read.gain *= rate;
            read.offset = f(read.offset);
            read.tilt *= rate;
        }
        for v in &mut self.values {
            *v = f(*v);
        }
        for m in &mut self.slopes {
            *m *= rate;
        }
        self
    }

    /// The value at the spine's end.
    pub(super) fn end(&self) -> f64 {
        self.values[self.values.len() - 1]
    }

    /// How far from the start the law holds its start value, and from how
    /// far it holds its end value to the end: the shares bounding its
    /// level spans at either end, where its values and slopes stand still.
    fn still_ends(&self) -> (f64, f64) {
        let n = self.shares.len();
        if self.function.is_some() {
            // A function's level stretches are not known: it is read as
            // moving all along.
            return (0.0, 1.0);
        }
        if self.linear {
            return if self.values[0] == self.values[1] {
                (1.0, 0.0)
            } else {
                (0.0, 1.0)
            };
        }
        let flat = |k: usize| {
            self.values[k] == self.values[k + 1]
                && self.slopes[k] == 0.0
                && self.slopes[k + 1] == 0.0
        };
        let mut lead = 0;
        while lead + 1 < n && flat(lead) {
            lead += 1;
        }
        let mut trail = n - 1;
        while trail > 0 && flat(trail - 1) {
            trail -= 1;
        }
        (self.shares[lead], self.shares[trail])
    }

    /// Whether the law holds `value` all along.
    pub(super) fn is_level(&self, value: f64) -> bool {
        self.function.is_none() && self.values.iter().all(|v| *v == value)
    }

    /// Whether the law runs straight from its start value to its end one.
    fn is_linear(&self) -> bool {
        self.linear
    }

    /// The law's values at [`FUNCTION_READS`] even steps from start to
    /// end, for a law given as a function.
    fn function_reads(&self) -> Option<Vec<f64>> {
        self.function.as_ref().map(|f| {
            (0..=FUNCTION_READS)
                .map(|k| f.at(f64::from(k) / f64::from(FUNCTION_READS)))
                .collect()
        })
    }

    /// The least value the law takes: each span keeps between its ends. A
    /// function's is the least it shows where it is read.
    fn least(&self) -> f64 {
        let reads = self.function_reads().unwrap_or_else(|| self.values.clone());
        reads.iter().copied().fold(f64::INFINITY, |a, b| {
            if b.is_nan() {
                f64::NEG_INFINITY
            } else {
                a.min(b)
            }
        })
    }

    /// How far the law runs, up and down, from start to end: each span
    /// runs one way only. A function runs at least as far as it shows
    /// where it is read.
    fn variation(&self) -> f64 {
        let reads = self.function_reads().unwrap_or_else(|| self.values.clone());
        reads.windows(2).map(|w| (w[1] - w[0]).abs()).sum()
    }

    /// The same law ending at `end`.
    fn ending_at(&self, end: f64) -> OgeomResult<Self> {
        if let Some(f) = &self.function {
            let read = FunctionRead {
                tilt: f.tilt + end - f.at(1.0),
                ..f.clone()
            };
            return Ok(Self {
                values: vec![read.at(0.0), read.at(1.0)],
                function: Some(read),
                ..self.clone()
            });
        }
        let points: Vec<(f64, f64)> = self.shares[1..]
            .iter()
            .zip(&self.values[1..])
            .map(|(s, v)| (*s, *v))
            .collect();
        let mut points = points;
        if let Some(last) = points.last_mut() {
            last.1 = end;
        }
        let law = if self.linear {
            SpineLaw::Linear(end)
        } else {
            SpineLaw::Through(points)
        };
        Self::of(&law, self.values[0], "twist")
    }
}

/// A pipe's scale and twist laws, read for evaluation.
#[derive(Debug, Clone)]
pub(super) struct Evolution {
    scale: Ramp,
    twist: Ramp,
}

impl Evolution {
    /// The scale and turn a share `lambda` of the way along, in the start
    /// frame's coordinates: about its origin, the spine, and its `z`, the
    /// spine's tangent.
    pub(super) fn local(&self, lambda: f64, tol: Tolerances) -> OgeomResult<Transform> {
        let (scale, twist) = (self.scale.at(lambda), self.twist.at(lambda));
        // A law given as a function is checked wherever it is read.
        if self.scale.function.is_some() && !(scale.is_finite() && scale > 0.0) {
            ogeom_bail!(
                Construction,
                "a scale law given as a function stays above naught; it gives {scale} \
                 where the sweep reads it"
            );
        }
        if self.twist.function.is_some() && !twist.is_finite() {
            ogeom_bail!(
                Construction,
                "a twist law takes finite values; it gives {twist} where the sweep reads it"
            );
        }
        let turn = Transform::rotation(
            Axis {
                location: Point::ORIGIN,
                direction: Direction::Z,
            },
            twist,
        );
        Ok(turn * Transform::scaling(Point::ORIGIN, scale, tol)?)
    }

    /// The shares of the way along bounding where the laws move the
    /// section: before the first and after the second both stand still.
    pub(super) fn still_ends(&self) -> (f64, f64) {
        let (scale, twist) = (self.scale.still_ends(), self.twist.still_ends());
        (scale.0.min(twist.0), scale.1.max(twist.1))
    }

    /// Whether the laws neither scale nor turn the section anywhere.
    pub(super) fn is_still(&self) -> bool {
        self.scale.is_level(1.0) && self.twist.is_level(0.0)
    }

    /// The laws over the stretch from share `from` to share `to`, against
    /// the share of that stretch and from the section as it stands at
    /// `from`: the scale over its value there, the twist less its value
    /// there. About the spine a scale and a turn commute, so the section
    /// at `from` carried on by these is the profile carried by the whole
    /// laws.
    pub(super) fn window(&self, from: f64, to: f64) -> Self {
        let (s0, t0) = (self.scale.at(from), self.twist.at(from));
        Self {
            scale: self.scale.window(from, to).mapped(|v| v / s0, 1.0 / s0),
            twist: self.twist.window(from, to).mapped(|v| v - t0, 1.0),
        }
    }
}

/// The length of a spine's edges and of the run up to each.
pub(super) struct SpineLengths {
    before: Vec<f64>,
    total: f64,
}

impl SpineLengths {
    /// The lengths along `curves`, a spine's edges in order.
    pub(super) fn of(curves: &[WalkCurve], tol: Tolerances) -> OgeomResult<Self> {
        let mut before = Vec::with_capacity(curves.len());
        let mut total = 0.0;
        for (curve, range, _) in curves {
            before.push(total);
            total += curve_length(curve, *range, tol)?;
        }
        if total <= tol.confusion() {
            ogeom_bail!(Construction, "the spine has no length");
        }
        Ok(Self { before, total })
    }

    /// The share of the spine's length run to parameter `t` on edge `edge`.
    pub(super) fn fraction(
        &self,
        curves: &[WalkCurve],
        edge: usize,
        t: f64,
        tol: Tolerances,
    ) -> OgeomResult<f64> {
        let (curve, range, reversed) = &curves[edge];
        let run = if *reversed {
            curve_length(curve, (t, range.1), tol)?
        } else {
            curve_length(curve, (range.0, t), tol)?
        };
        Ok(((self.before[edge] + run) / self.total).clamp(0.0, 1.0))
    }
}

/// Sweep a planar profile along a spine under a frame law, the section
/// scaling and twisting as it goes.
///
/// The profile (a closed wire, or a face whose holes ride along) stands
/// square to the spine at its start, as for a pipe shell. At each point of
/// the spine the section is the profile scaled about the spine by
/// `1 + (scale - 1) λ` and turned about it by `twist λ`, `λ` the share of
/// the spine's length run there, then carried by `law`. With no scale and
/// no twist this is [`make_pipe_shell_with`] with mitred corners; it is
/// [`make_pipe_shell_laws`] with linear laws.
///
/// Scaled without a twist down one straight edge under the
/// rotation-minimizing, binormal or fixed law (none of which turns the
/// section along a line), the solid is the ruled loft from the profile to
/// its scaled image, exact; a scale of zero closes it on the spine's end,
/// the exact pyramid or cone. Otherwise each loop of the profile is
/// skinned through its placed sections and held to `tolerance` between
/// them against the law, along an open spine with no sharp corner. Round
/// a smooth closed spine the section comes home to itself only unscaled,
/// and twisted whole turns or a part of one that turns it onto itself (a
/// square a quarter turn, a regular hexagon a sixth, a circle about its
/// centre any angle): the skin is closed round the spine, each hole a
/// void tunnel of its own shell. Each loop of the profile must come home
/// onto itself, a faceted one edge onto edge.
///
/// Along an open spine each hole sweeps its own solid under the same law,
/// cut from the outer one.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `scale` is negative or `twist` not finite; if the profile is not a
/// closed planar wire or face standing square to the spine at its start;
/// if a skinned sweep is asked round a sharp corner, of the fixed law, or
/// to a scale of zero; if a section is scaled, or twisted by a part of a
/// turn that does not bring each of its loops home onto itself, round a
/// closed spine; and as
/// [`make_pipe_shell_with`] or the law refuses the spine.
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the skin
/// cannot reach `tolerance`.
pub fn make_pipe_shell_evolving(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    evolution: PipeEvolution,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    make_pipe_shell_laws(
        model,
        profile,
        spine,
        law,
        &evolution.into(),
        PipeCorners::Mitre,
        tolerance,
        tol,
    )
}

/// [`make_pipe_shell_evolving`] with a law of its own for the scale and
/// for the twist: at each point of the spine the section is the profile
/// scaled about the spine by `laws.scale` and turned about it by
/// `laws.twist`, each read at the share of the spine's length run there,
/// then carried by `law`.
///
/// Linear laws build as [`make_pipe_shell_evolving`] does. Any other law
/// skins each loop through its placed sections, held to `tolerance`
/// against the laws; round a smooth closed spine the twist may follow any
/// law that ends where a linear one may, unscaled. A scale law through
/// values is above naught all along: a cubic span keeps between its ends.
/// A law given as a function is read wherever the sweep places a section
/// or checks its skin, and is never taken as standing still: every run it
/// moves along is skinned, and round a closed spine it twists only, as a
/// scale it is refused there.
///
/// Where `laws.contact` asks, the profile is first moved rigidly onto the
/// spine's start as [`PipeContact`] defines, and the moved profile is
/// swept as above; the history records the given profile as generating
/// the result.
///
/// Round sharp corners a face sweeps run by run as
/// [`make_pipe_shell_with`] does, each corner turned as `corners` asks:
/// each smooth run is swept under the laws over its share of the spine's
/// length, from the section it starts with (a run the laws neither scale
/// nor turn is the prism or pipe shell of that section), and the runs
/// meet at each corner on the section the laws give there.
///
/// # Errors
///
/// As [`make_pipe_shell_evolving`] and, round sharp corners, as
/// [`make_pipe_shell_with`]; and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a
/// law's shares do not rise strictly to one, a value is not finite, a
/// scale law through values reaches naught or below, or a law given as a
/// function does not give the start value at the spine's start or, as a
/// scale, gives naught or below where it is read; and, where
/// `laws.contact` moves the profile, if the profile is not planar, or its
/// plane holds the spine's start tangent under a correction.
#[allow(clippy::too_many_arguments)]
pub fn make_pipe_shell_laws(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    laws: &PipeLaws,
    corners: PipeCorners,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if let SpineLaw::Linear(scale) = laws.scale
        && (!scale.is_finite() || scale < 0.0)
    {
        ogeom_bail!(
            Construction,
            "a section scales to {scale} at the spine's end; a scale is naught or more"
        );
    }
    if let SpineLaw::Linear(twist) = laws.twist
        && !twist.is_finite()
    {
        ogeom_bail!(Construction, "a twist of {twist} is not an angle");
    }
    let evolution = Evolution {
        scale: Ramp::of(&laws.scale, 1.0, "scale")?,
        twist: Ramp::of(&laws.twist, 0.0, "twist")?,
    };
    // The spine is swept where it stands.
    let placed_spine = placed_path(model, spine, tol)?;
    if placed_spine != *spine {
        let mut built = make_pipe_shell_laws(
            model,
            profile,
            &placed_spine,
            law,
            laws,
            corners,
            tolerance,
            tol,
        )?;
        built.history.generate(spine, built.shape.clone());
        return Ok(built);
    }
    let placed = profile_in_contact(model, profile, spine, laws.contact, tol)?;
    if placed != *profile {
        let mut built = make_pipe_shell_laws(
            model,
            &placed,
            spine,
            law,
            &PipeLaws {
                contact: PipeContact::AsPlaced,
                ..laws.clone()
            },
            corners,
            tolerance,
            tol,
        )?;
        built.history.generate(profile, built.shape.clone());
        return Ok(built);
    }
    if matches!(laws.scale, SpineLaw::Through(_) | SpineLaw::Function(_))
        && evolution.scale.least() <= 0.0
    {
        ogeom_bail!(
            Construction,
            "a scale law through values or given as a function stays above naught; a \
             section shrinks to a point only linearly"
        );
    }
    if evolution.is_still() {
        return make_pipe_shell_with(model, profile, spine, law, corners, tolerance, tol);
    }
    let rings: Vec<Shape> = match model.kind_of(profile)? {
        ShapeType::Face => explore(model, profile, Filter::OfType(ShapeType::Wire))?,
        ShapeType::Wire => vec![profile.clone()],
        other => ogeom_bail!(
            Construction,
            "a pipe shell sweeps a planar wire or face, not a {other:?}"
        ),
    };
    let Some(plane) = find_plane(model, profile, tol)? else {
        ogeom_bail!(Construction, "a pipe shell sweeps a planar profile");
    };
    let stations = shell_stations(model, spine, tol)?;
    profile_square_at_start(model, &plane, spine, stations[0].at, tol)?;
    if rings.len() > 1 && evolution.scale.end() == 0.0 {
        ogeom_bail!(
            Construction,
            "a profile with holes cannot close on a point; its holes would meet"
        );
    }
    if spine_closes(model, spine, tol)?
        && !has_sharp_corner(model, spine, tol)?
        && !closes_on_a_corner(model, spine, tol)?
        && !matches!(law, PipeLaw::Fixed)
    {
        if !evolution.scale.is_level(1.0) {
            ogeom_bail!(
                Construction,
                "a section scaling round a closed spine would not come home to itself; \
                 it scales along an open one"
            );
        }
        // Whole turns bring the section home to itself; a part of a turn
        // only where its symmetry does, which the closed shell checks.
        let end = evolution.twist.end();
        let turns = end / core::f64::consts::TAU;
        let twist = if (turns - turns.round()).abs() * core::f64::consts::TAU > tol.angular() {
            evolution.twist
        } else {
            evolution
                .twist
                .ending_at(turns.round() * core::f64::consts::TAU)?
        };
        return closed_law_pipe(model, profile, spine, law, &twist, tolerance, tol);
    }
    if !matches!(law, PipeLaw::Fixed)
        && (has_sharp_corner(model, spine, tol)? || closes_on_a_corner(model, spine, tol)?)
    {
        return law_cornered(
            model,
            profile,
            spine,
            law,
            corners,
            Some(&evolution),
            tolerance,
            tol,
        );
    }
    let shape = evolved_run(
        model, profile, spine, law, &evolution, false, tolerance, tol,
    )?;
    let mut history = History::new();
    history.generate(spine, shape.clone());
    history.generate(profile, shape.clone());
    Ok(Built::new(shape, history))
}

/// A planar profile swept along a spine with no sharp corner under the
/// evolving laws: each loop its own solid, the holes' cut from the outer.
/// Where `exact_ends` asks it, a skinned loop of one edge with an exact
/// rational form ends on its end sections exactly.
#[allow(clippy::too_many_arguments)]
pub(super) fn evolved_run(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    evolution: &Evolution,
    exact_ends: bool,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let rings: Vec<Shape> = match model.kind_of(profile)? {
        ShapeType::Face => explore(model, profile, Filter::OfType(ShapeType::Wire))?,
        _ => vec![profile.clone()],
    };
    let mut solids = Vec::with_capacity(rings.len());
    for ring in &rings {
        solids.push(evolving_ring(
            model, ring, spine, law, evolution, exact_ends, tolerance, tol,
        )?);
    }
    let mut shape = solids[0].clone();
    for hole in &solids[1..] {
        shape = ogeom_bool::cut(model, &shape, hole, tol)?.shape;
    }
    Ok(shape)
}

/// One closed ring of the profile swept under the evolving law, as a solid.
#[allow(clippy::too_many_arguments)]
fn evolving_ring(
    model: &mut Model,
    ring: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    evolution: &Evolution,
    exact_ends: bool,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let level = matches!(
        law,
        PipeLaw::RotationMinimizing | PipeLaw::Binormal(_) | PipeLaw::Fixed
    );
    let scale = evolution.scale.end();
    if evolution.twist.is_level(0.0)
        && evolution.scale.is_linear()
        && level
        && let Some(legs) = straight_legs(model, spine, tol)?
        && legs.len() == 1
    {
        let (start, run) = legs[0];
        let end = if scale == 0.0 {
            model.add_vertex(VertexData::new(start + run))
        } else {
            let motion = Transform::translation(run) * Transform::scaling(start, scale, tol)?;
            moved_ring(model, ring, &motion, tol)?
        };
        // A polygon or a circle lofts to its image as planes or a cone;
        // any other ring as the ruled walls between matched edges.
        let shape = if scale == 0.0 || analytic_ring(model, ring, start, tol)? {
            make_loft(model, ring, &end, tol)?.shape
        } else {
            crate::loft::make_loft_solid(model, &[ring.clone(), end], &[], true, tol)?.shape
        };
        return Ok(shape);
    }
    if matches!(law, PipeLaw::Fixed) {
        ogeom_bail!(
            Construction,
            "a fixed section scales or twists only down one straight edge, unturned"
        );
    }
    if scale == 0.0 {
        ogeom_bail!(
            Construction,
            "a section shrinks to a point only down one straight edge, unturned"
        );
    }
    let stations = shell_stations(model, spine, tol)?;
    let ring_spine =
        stations[0].at.distance(stations[stations.len() - 1].at) <= tol.confusion() * 10.0;
    if ring_spine || has_sharp_corner(model, spine, tol)? {
        ogeom_bail!(
            Construction,
            "a scaling or twisting section sweeps an open spine with no sharp corner"
        );
    }
    // As many stations as the law's turning or the twist asks, each edge
    // its share by length.
    let turned = densified(model, spine, stations.clone(), law, tolerance, tol)?.len();
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let twisted = (evolution.twist.variation() / super::MOST_TURN).ceil() as usize + 1;
    let wanted = turned.max(twisted);
    let lengths = edge_lengths(model, spine, tol)?;
    let total: f64 = lengths.iter().sum();
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    let share = |e: usize| -> u32 {
        let length = lengths.get(e).copied().unwrap_or(0.0);
        ((wanted as f64 * length / total).ceil() as u32).max(2)
    };
    let stations = evenly_by(model, spine, stations, &share, false, tol)?;
    Ok(law_loft(
        model,
        ring,
        spine,
        &stations,
        law,
        Some(evolution),
        exact_ends,
        tolerance,
        tol,
    )?
    .shape)
}

/// A planar profile swept round a smooth closed spine under `law`, turned
/// about the spine by the `twist` law of the length run (ending on whole
/// turns, or a part of one the section's symmetry closes): stations as
/// many as the law's turning or the twist asks, each edge its share by
/// length, and the skin closed on itself.
pub(super) fn closed_law_pipe(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    twist: &Ramp,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let stations = shell_stations(model, spine, tol)?;
    let turned = densified(model, spine, stations.clone(), law, tolerance, tol)?.len();
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let twisted = (twist.variation() / super::MOST_TURN).ceil() as usize + 1;
    let wanted = turned.max(twisted);
    let lengths = edge_lengths(model, spine, tol)?;
    let total: f64 = lengths.iter().sum();
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    let share = |e: usize| -> u32 {
        let length = lengths.get(e).copied().unwrap_or(0.0);
        ((wanted as f64 * length / total).ceil() as u32).max(2)
    };
    let stations = evenly_by(model, spine, stations, &share, false, tol)?;
    // Where two edges meet the walk stands twice; the loop takes it once,
    // the frame between it and the next station read off the next edge.
    let mut kept: Vec<SpineStation> = Vec::with_capacity(stations.len());
    for station in stations {
        if kept
            .last()
            .is_some_and(|p| p.at.distance(station.at) <= tol.confusion())
        {
            continue;
        }
        kept.push(station);
    }
    closed_pipe_shell(model, profile, spine, kept, law, twist, tolerance, tol)
}

/// Whether a ring is all straight edges, or one circle about `axis_at`:
/// what [`make_loft`] lofts to its image scaled about `axis_at` exactly.
fn analytic_ring(
    model: &Model,
    ring: &Shape,
    axis_at: Point,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let edges = explore(model, ring, Filter::OfType(ShapeType::Edge))?;
    let mut kinds = Vec::with_capacity(edges.len());
    for edge in &edges {
        kinds.push(spine_curve_of(model, edge)?.0);
    }
    Ok(kinds.iter().all(|c| matches!(c, Curve::Line(_)))
        || (kinds.len() == 1
            && matches!(&kinds[0], Curve::Circle(c)
                if c.circle().centre().distance(axis_at) <= tol.confusion())))
}

/// Each spine edge's length, in the spine's order.
fn edge_lengths(model: &Model, spine: &Shape, tol: Tolerances) -> OgeomResult<Vec<f64>> {
    let edges: Vec<Shape> = match model.kind_of(spine)? {
        ShapeType::Edge => vec![spine.clone()],
        _ => model.ordered_children_of(spine)?,
    };
    let mut out = Vec::with_capacity(edges.len());
    for edge in &edges {
        let (curve, range) = spine_curve_of(model, edge)?;
        out.push(curve_length(&curve, range, tol)?);
    }
    Ok(out)
}
