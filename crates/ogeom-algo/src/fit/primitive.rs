//! Best fits of primitives to point sets, with their residuals.
//!
//! Where [`recognize_points`](crate::recognize_points) decides whether points
//! *are* a surface, these fit the surface regardless and report how well it
//! fits: every point's signed residual, their RMS over the points kept, and
//! the signed extremes. Each fit is orthogonal least squares, the sum of the
//! squared distances from the points to the shape made least.
//!
//! A plane and a line are closed forms: the points' centroid and the
//! covariance's smallest or largest direction. Every other kind starts from
//! a closed-form estimate and is refined by damped Gauss-Newton steps on the
//! distances themselves: a sphere and a circle from the algebraic fit
//! through the `|c|² − r²` substitution, and a cylinder, a cone and a torus
//! from the axis every normal line meets, estimated from the normals when
//! the caller gives them and from each point's neighbours when not.
//!
//! The fits run on the calling thread, in point order, and draw no
//! randomness but a fixed sequence: the same points give the same answer
//! bit for bit on any machine and any thread count.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{
    Axis, Axis2, Circle, Circle2, Cone, Cylinder, Direction, Direction2, Frame, Plane, Point,
    Point2, Sphere, Torus, Vector,
};

use crate::Cutoff;
use crate::recognize::{self, Canonical};

/// Which side of a closed surface or curve the material lies on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
    /// The material is inside: a boss, a shaft, a ball. A residual is
    /// positive outside the shape.
    Inside,
    /// The material is outside: a bore, a socket. A residual is positive
    /// inside the shape.
    Outside,
}

/// What a fit holds fixed instead of fitting.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct FitFixed {
    /// A sphere's, a cylinder's or a circle's radius, or a torus's tube
    /// (minor) radius: a precision ball partly scanned fits its centre only.
    pub radius: Option<f64>,
    /// A plane's normal, a line's direction, or the axis direction of a
    /// cylinder, a cone, a torus or a circle in space; its sense is free.
    pub direction: Option<Direction>,
    /// The whole axis of a cylinder, a cone, a torus or a circle in space,
    /// as a line; its sense is free.
    pub axis: Option<Axis>,
    /// A cone's half angle, in radians.
    pub half_angle: Option<f64>,
}

/// How to fit.
///
/// `S` is the shape fitted, the type of the starting guess.
#[derive(Debug, Clone, PartialEq)]
pub struct FitOptions<'a, S> {
    /// A normal for each point, used for the start and to tell which side
    /// the material is on; never in the least squares itself, which
    /// measures the points alone. Their sense is the material's outside.
    pub normals: Option<&'a [Vector]>,
    /// Which points are kept between rounds.
    pub cutoff: Cutoff,
    /// Start from a consensus over small subsets of the points, scored by
    /// their median residual, rather than from all of them: a start that
    /// many outliers cannot pull away, and whose residuals choose the
    /// first round's points.
    pub robust_start: bool,
    /// What is held instead of fitted.
    pub fixed: FitFixed,
    /// Start from this shape rather than from the estimate; a robust start
    /// then only chooses the first round's points, by their residuals from
    /// it. A closed-form fit (a plane, a line) has no start, and ignores it.
    pub guess: Option<S>,
    /// Which side the material is on. Where it is not given and normals
    /// are, the normals decide.
    pub side: Option<Side>,
    /// The radius a fit must land within, inclusive: the radius
    /// [`FitFixed::radius`] would hold.
    pub radius_range: Option<(f64, f64)>,
}

impl<S> Default for FitOptions<'_, S> {
    fn default() -> Self {
        Self {
            normals: None,
            cutoff: Cutoff::All,
            robust_start: false,
            fixed: FitFixed::default(),
            guess: None,
            side: None,
            radius_range: None,
        }
    }
}

/// A fitted shape and how well it fits.
#[derive(Debug, Clone, PartialEq)]
pub struct Fit<S> {
    /// The shape.
    pub shape: S,
    /// The RMS of the residuals of the points kept.
    pub sigma: f64,
    /// The least and the greatest residual of the points kept; their
    /// difference is the form error, peak to peak.
    pub form: (f64, f64),
    /// Whether each point was kept by the last round.
    pub used: Vec<bool>,
    /// Each point's signed distance to the shape: positive off the
    /// material's side where the side is known, else away from the centre,
    /// the axis or the spine; along the normal for a plane, to the left for
    /// a line in the plane. A line in space has no side, and its residual
    /// is the distance.
    pub residuals: Vec<f64>,
    /// The least-squares steps taken, over every round.
    pub iterations: u32,
    /// The rounds fitted: one, unless the cutoff dropped points.
    pub rounds: u32,
    /// Whether the last round's steps converged and the kept points stopped
    /// changing.
    pub converged: bool,
    /// The side the residuals are signed by, given or read from the
    /// normals; `None` for a kind without one or where neither said.
    pub side: Option<Side>,
}

/// The plane through points in space.
///
/// # Errors
///
/// As every fit: [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// with too few points (three, or one with the normal fixed) or with them
/// all within the confusion distance of one line, and
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for an option this
/// kind has no use for (a side, a radius, an axis) or a point that is not
/// finite.
pub fn plane(
    points: &[Point],
    options: &FitOptions<'_, Plane>,
    tol: Tolerances,
) -> OgeomResult<Fit<Plane>> {
    let guess = options.guess.map(|p| State::about(p.origin(), p.normal()));
    run(Kind::Plane, points, options, guess, tol)?
        .shaped(|s| Ok(Plane::new(Frame::about(s.origin, s.axis))))
}

/// The line through points in space.
///
/// # Errors
///
/// As [`plane`], with two points needed, or one with the direction fixed,
/// and refused where they all lie within the confusion distance of one
/// point.
pub fn line3(
    points: &[Point],
    options: &FitOptions<'_, Axis>,
    tol: Tolerances,
) -> OgeomResult<Fit<Axis>> {
    let guess = options.guess.map(|a| State::about(a.location, a.direction));
    run(Kind::Line3, points, options, guess, tol)?.shaped(|s| Ok(Axis::new(s.origin, s.axis)))
}

/// The sphere through points.
///
/// # Errors
///
/// As [`plane`], with four points needed, or three with the radius fixed;
/// refused where they lie on one plane, and with
/// [`OgeomError::Range`](ogeom_core::OgeomError::Range) where the radius
/// lands outside [`FitOptions::radius_range`].
pub fn sphere(
    points: &[Point],
    options: &FitOptions<'_, Sphere>,
    tol: Tolerances,
) -> OgeomResult<Fit<Sphere>> {
    let guess = options.guess.map(|s| State {
        radius: s.radius(),
        ..State::about(s.centre(), Direction::Z)
    });
    run(Kind::Sphere, points, options, guess, tol)?
        .shaped(|s| Sphere::centred(s.origin, s.radius, tol))
}

/// The cylinder through points.
///
/// # Errors
///
/// As [`sphere`], with five points needed, fewer with parts of it fixed;
/// refused where no axis is determined.
pub fn cylinder(
    points: &[Point],
    options: &FitOptions<'_, Cylinder>,
    tol: Tolerances,
) -> OgeomResult<Fit<Cylinder>> {
    let guess = options.guess.map(|c| State {
        radius: c.radius(),
        ..State::about(c.frame().origin(), c.frame().z())
    });
    run(Kind::Cylinder, points, options, guess, tol)?
        .shaped(|s| Cylinder::new(Frame::about(s.origin, s.axis), s.radius, tol))
}

/// The cone through points.
///
/// The cone's frame sits on its axis level with the points' centroid, and
/// widens along its axis.
///
/// # Errors
///
/// As [`sphere`], with six points needed, fewer with parts of it fixed;
/// refused where no axis is determined or where the points have no taper,
/// lying on a cylinder. A radius range does not apply: a cone's radius
/// changes along it.
pub fn cone(
    points: &[Point],
    options: &FitOptions<'_, Cone>,
    tol: Tolerances,
) -> OgeomResult<Fit<Cone>> {
    let guess = options.guess.map(|c| State {
        radius: c.reference_radius(),
        second: c.half_angle(),
        ..State::about(c.frame().origin(), c.frame().z())
    });
    run(Kind::Cone, points, options, guess, tol)?.shaped(|s| {
        Cone::new(Frame::about(s.origin, s.axis), s.radius, s.second, tol).map_err(|_| {
            ogeom_core::ogeom_err!(
                Construction,
                "the points have no taper a cone can take: half angle {} at radius {}",
                s.second,
                s.radius
            )
        })
    })
}

/// The torus through points.
///
/// [`FitFixed::radius`] and [`FitOptions::radius_range`] are the tube's
/// (minor) radius.
///
/// # Errors
///
/// As [`sphere`], with seven points needed, fewer with parts of it fixed;
/// refused where no axis is determined.
pub fn torus(
    points: &[Point],
    options: &FitOptions<'_, Torus>,
    tol: Tolerances,
) -> OgeomResult<Fit<Torus>> {
    let guess = options.guess.map(|t| State {
        radius: t.minor_radius(),
        second: t.major_radius(),
        ..State::about(t.centre(), t.frame().z())
    });
    run(Kind::Torus, points, options, guess, tol)?
        .shaped(|s| Torus::new(Frame::about(s.origin, s.axis), s.second, s.radius, tol))
}

/// The circle through points in space: the distance a residual measures
/// is to the circle itself, out of its plane as well as across it, and is
/// signed by which side of the cylinder through the circle the point lies.
///
/// # Errors
///
/// As [`sphere`], with three points needed; refused where they lie on one
/// line.
pub fn circle3(
    points: &[Point],
    options: &FitOptions<'_, Circle>,
    tol: Tolerances,
) -> OgeomResult<Fit<Circle>> {
    let guess = options.guess.map(|c| State {
        radius: c.radius(),
        ..State::about(c.centre(), c.frame().z())
    });
    run(Kind::Circle3, points, options, guess, tol)?
        .shaped(|s| Circle::new(Frame::about(s.origin, s.axis), s.radius, tol))
}

/// The circle through points in the plane.
///
/// # Errors
///
/// As [`sphere`], with three points needed; refused where they lie on one
/// line. Normals, a direction and an axis do not apply.
pub fn circle2(
    points: &[Point2],
    options: &FitOptions<'_, Circle2>,
    tol: Tolerances,
) -> OgeomResult<Fit<Circle2>> {
    let lifted: Vec<Point> = points.iter().map(|p| p.to_3d()).collect();
    let guess = options.guess.map(|c| State {
        radius: c.radius(),
        ..State::about(c.centre().to_3d(), Direction::Z)
    });
    run(Kind::Circle2, &lifted, options, guess, tol)?
        .shaped(|s| Circle2::centred(s.origin.xy(), s.radius, tol))
}

/// The line through points in the plane.
///
/// # Errors
///
/// As [`line3`]. Normals, a side and anything fixed do not apply.
pub fn line2(
    points: &[Point2],
    options: &FitOptions<'_, Axis2>,
    tol: Tolerances,
) -> OgeomResult<Fit<Axis2>> {
    let lifted: Vec<Point> = points.iter().map(|p| p.to_3d()).collect();
    run(Kind::Line2, &lifted, options, None, tol)?.shaped(|s| {
        Ok(Axis2 {
            location: s.origin.xy(),
            direction: Direction2::new(s.axis.vector().xy(), tol)?,
        })
    })
}

/// The kinds of shape, as the fitting machinery tells them apart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Kind {
    Plane,
    Line3,
    Line2,
    Sphere,
    Circle2,
    Circle3,
    Cylinder,
    Cone,
    Torus,
}

impl Kind {
    const fn name(self) -> &'static str {
        match self {
            Self::Plane => "plane",
            Self::Line3 => "line",
            Self::Line2 => "line in the plane",
            Self::Sphere => "sphere",
            Self::Circle2 => "circle in the plane",
            Self::Circle3 => "circle",
            Self::Cylinder => "cylinder",
            Self::Cone => "cone",
            Self::Torus => "torus",
        }
    }

    /// Fitted in closed form, with no start and no steps.
    const fn closed(self) -> bool {
        matches!(self, Self::Plane | Self::Line3 | Self::Line2)
    }

    /// A kind in the plane, which has no normals to take.
    const fn planar(self) -> bool {
        matches!(self, Self::Line2 | Self::Circle2)
    }

    /// A kind that encloses a side, and so has residuals signed by it.
    const fn sided(self) -> bool {
        !matches!(self, Self::Plane | Self::Line3 | Self::Line2)
    }

    /// A kind with the one radius [`FitFixed::radius`] holds.
    const fn radial(self) -> bool {
        matches!(
            self,
            Self::Sphere | Self::Circle2 | Self::Circle3 | Self::Cylinder | Self::Torus
        )
    }

    /// A kind with a direction to hold.
    const fn directed(self) -> bool {
        matches!(
            self,
            Self::Plane | Self::Line3 | Self::Circle3 | Self::Cylinder | Self::Cone | Self::Torus
        )
    }

    /// A kind with an axis line to hold.
    const fn axial(self) -> bool {
        matches!(
            self,
            Self::Circle3 | Self::Cylinder | Self::Cone | Self::Torus
        )
    }

    /// How many distances each point contributes: a circle in space is
    /// missed across its plane and within it.
    const fn components(self) -> usize {
        if matches!(self, Self::Circle3) { 2 } else { 1 }
    }

    /// A surface of revolution, whose start is read from normals.
    const fn revolved(self) -> bool {
        matches!(self, Self::Cylinder | Self::Cone | Self::Torus)
    }
}

/// A shape of any kind as the fits move it: a point (a centre, or a point
/// on the axis or the plane), a direction (a normal or an axis), a radius
/// (a cone's at the point, a torus's tube) and a second size (a cone's half
/// angle, a torus's major radius).
#[derive(Debug, Clone, Copy, PartialEq)]
struct State {
    origin: Point,
    axis: Direction,
    radius: f64,
    second: f64,
}

impl State {
    const fn about(origin: Point, axis: Direction) -> Self {
        Self {
            origin,
            axis,
            radius: 0.0,
            second: 0.0,
        }
    }
}

/// The parameters a step moves: two tilts of the direction, two moves of
/// the point across it, one along it, the radius and the second size. A
/// sphere and a circle in the plane move their centre along the world's
/// axes through the same three moves.
const PARAMETERS: usize = 7;
const TILT: [usize; 2] = [0, 1];
const ACROSS: [usize; 2] = [2, 3];
const ALONG: usize = 4;
const RADIUS: usize = 5;
const SECOND: usize = 6;

/// Which parameters a kind fits, with what the caller holds left out.
fn free(kind: Kind, fixed: &FitFixed) -> [bool; PARAMETERS] {
    let mut on = [false; PARAMETERS];
    let tilt = fixed.direction.is_none() && fixed.axis.is_none();
    let across = fixed.axis.is_none();
    let radius = fixed.radius.is_none();
    match kind {
        Kind::Plane | Kind::Line3 | Kind::Line2 => {}
        Kind::Sphere => {
            on[ACROSS[0]] = true;
            on[ACROSS[1]] = true;
            on[ALONG] = true;
            on[RADIUS] = radius;
        }
        Kind::Circle2 => {
            on[ACROSS[0]] = true;
            on[ACROSS[1]] = true;
            on[RADIUS] = radius;
        }
        Kind::Cylinder | Kind::Cone | Kind::Torus | Kind::Circle3 => {
            on[TILT[0]] = tilt;
            on[TILT[1]] = tilt;
            on[ACROSS[0]] = across;
            on[ACROSS[1]] = across;
            on[ALONG] = matches!(kind, Kind::Torus | Kind::Circle3);
            on[RADIUS] = radius;
            on[SECOND] = match kind {
                Kind::Cone => fixed.half_angle.is_none(),
                Kind::Torus => true,
                _ => false,
            };
        }
    }
    on
}

/// The fewest points that determine the fit.
pub(super) fn needed(kind: Kind, fixed: &FitFixed) -> usize {
    let held = fixed.direction.is_some();
    match kind {
        Kind::Plane => {
            if held {
                1
            } else {
                3
            }
        }
        Kind::Line3 | Kind::Line2 => {
            if held {
                1
            } else {
                2
            }
        }
        _ => {
            let count = free(kind, fixed).iter().filter(|on| **on).count();
            count.div_ceil(kind.components()).max(1)
        }
    }
}

/// A fitted surface asked about any point: its residual there, signed as
/// the fit signed its own, and the direction that residual grows in.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Model {
    kind: Kind,
    state: State,
    /// −1 where the material is outside the shape, so the residual and the
    /// normal are the material's; 1 otherwise.
    sense: f64,
    /// Whether the fit knew the material's side, or for a plane, its normal
    /// was turned by the normals the fit was given.
    oriented: bool,
    /// The parameters the fit moved.
    free: [bool; PARAMETERS],
    /// The inverse of the fit's normal matrix over the parameters it moved,
    /// over its kept points: the parameters' covariance in units of the
    /// noise's variance. Zero until [`Model::measured_over`].
    spread: [[f64; PARAMETERS]; PARAMETERS],
}

impl Model {
    fn sided<S>(kind: Kind, state: State, fit: &Fit<S>) -> Self {
        Self {
            kind,
            state,
            sense: if fit.side == Some(Side::Outside) {
                -1.0
            } else {
                1.0
            },
            oriented: fit.side.is_some(),
            free: [false; PARAMETERS],
            spread: [[0.0; PARAMETERS]; PARAMETERS],
        }
    }

    /// A plane fitted with normals: its normal faces them.
    pub(super) fn plane(fit: &Fit<Plane>) -> Self {
        Self {
            kind: Kind::Plane,
            state: State::about(fit.shape.origin(), fit.shape.normal()),
            sense: 1.0,
            oriented: true,
            free: [false; PARAMETERS],
            spread: [[0.0; PARAMETERS]; PARAMETERS],
        }
    }

    pub(super) fn sphere(fit: &Fit<Sphere>) -> Self {
        let s = &fit.shape;
        let state = State {
            radius: s.radius(),
            ..State::about(s.centre(), Direction::Z)
        };
        Self::sided(Kind::Sphere, state, fit)
    }

    pub(super) fn cylinder(fit: &Fit<Cylinder>) -> Self {
        let c = &fit.shape;
        let state = State {
            radius: c.radius(),
            ..State::about(c.frame().origin(), c.frame().z())
        };
        Self::sided(Kind::Cylinder, state, fit)
    }

    pub(super) fn cone(fit: &Fit<Cone>) -> Self {
        let c = &fit.shape;
        let state = State {
            radius: c.reference_radius(),
            second: c.half_angle(),
            ..State::about(c.frame().origin(), c.frame().z())
        };
        Self::sided(Kind::Cone, state, fit)
    }

    pub(super) fn torus(fit: &Fit<Torus>) -> Self {
        let t = &fit.shape;
        let state = State {
            radius: t.minor_radius(),
            second: t.major_radius(),
            ..State::about(t.centre(), t.frame().z())
        };
        Self::sided(Kind::Torus, state, fit)
    }

    /// The point's residual, signed as the fit's residuals are.
    pub(super) fn residual(&self, p: Point) -> f64 {
        self.sense * residual(self.kind, &self.state, p)
    }

    /// The unit direction the point's residual grows in: the material's
    /// outward normal where the model is oriented.
    pub(super) fn normal(&self, p: Point) -> Vector {
        outward(self.kind, &self.state, p) * self.sense
    }

    /// Whether [`Model::normal`] has the material's sense rather than an
    /// arbitrary one.
    pub(super) const fn oriented(&self) -> bool {
        self.oriented
    }

    /// The model with the spread of its parameters measured over the
    /// points its fit kept, `fixed` as the fit held it.
    ///
    /// The normal matrix is inverted over its eigenvalues above the
    /// rounding floor, as the fit's own steps are: a direction the points
    /// do not determine is one the fit did not move either.
    pub(super) fn measured_over(
        mut self,
        points: &[Point],
        used: &[bool],
        fixed: &FitFixed,
    ) -> Self {
        let free = if self.kind == Kind::Plane {
            let mut on = [false; PARAMETERS];
            on[ALONG] = true;
            let tilt = fixed.direction.is_none();
            on[TILT[0]] = tilt;
            on[TILT[1]] = tilt;
            on
        } else {
            free(self.kind, fixed)
        };
        let index: Vec<usize> = (0..PARAMETERS).filter(|k| free[*k]).collect();
        let size = index.len();
        self.free = free;
        self.spread = [[0.0; PARAMETERS]; PARAMETERS];
        if size == 0 {
            return self;
        }
        let frame = basis(self.kind, &self.state);
        let mut normal = nalgebra::DMatrix::<f64>::zeros(size, size);
        for (p, keep) in points.iter().zip(used) {
            if !*keep {
                continue;
            }
            let (rows, count) = self.rows(frame, *p);
            for (_, d) in &rows[..count] {
                for (i, &ki) in index.iter().enumerate() {
                    for (j, &kj) in index.iter().enumerate() {
                        normal[(i, j)] += d[ki] * d[kj];
                    }
                }
            }
        }
        let eigen = nalgebra::SymmetricEigen::new(normal);
        let largest = eigen.eigenvalues.iter().copied().fold(0.0, f64::max);
        #[allow(clippy::cast_precision_loss, reason = "a handful of parameters")]
        let floor = largest * f64::EPSILON * size as f64;
        for (k, value) in eigen.eigenvalues.iter().enumerate() {
            if *value > floor {
                let v = eigen.eigenvectors.column(k);
                for (i, &ki) in index.iter().enumerate() {
                    for (j, &kj) in index.iter().enumerate() {
                        self.spread[ki][kj] += v[i] * v[j] / value;
                    }
                }
            }
        }
        self
    }

    /// How many parameters the fit moved.
    pub(super) fn parameters(&self) -> usize {
        self.free.iter().filter(|on| **on).count()
    }

    /// How much the fit's own uncertainty moves its residual at `p`, in
    /// units of the noise's variance: the residual of a point the fit did
    /// not see varies as the noise's variance times one plus this.
    pub(super) fn leverage(&self, p: Point) -> f64 {
        let (rows, count) = self.rows(basis(self.kind, &self.state), p);
        let mut sum = 0.0;
        for (_, d) in &rows[..count] {
            for i in 0..PARAMETERS {
                if !self.free[i] {
                    continue;
                }
                for j in 0..PARAMETERS {
                    if self.free[j] {
                        sum += d[i] * self.spread[i][j] * d[j];
                    }
                }
            }
        }
        sum
    }

    /// A point's distances and their derivatives, a plane's included: its
    /// point moved along its normal and its normal tilted.
    fn rows(&self, frame: [Vector; 3], p: Point) -> ([(f64, [f64; PARAMETERS]); 2], usize) {
        if self.kind != Kind::Plane {
            return distances(self.kind, &self.state, frame, p);
        }
        let [e1, e2, a] = frame;
        let w = p - self.state.origin;
        let mut d = [0.0; PARAMETERS];
        d[TILT[0]] = w.dot(e1);
        d[TILT[1]] = w.dot(e2);
        d[ALONG] = -1.0;
        ([(w.dot(a), d), (0.0, [0.0; PARAMETERS])], 1)
    }
}

/// A fit's state and the bookkeeping the public result is made from.
struct Outcome {
    state: State,
    sigma: f64,
    form: (f64, f64),
    used: Vec<bool>,
    residuals: Vec<f64>,
    iterations: u32,
    rounds: u32,
    converged: bool,
    side: Option<Side>,
}

impl Outcome {
    fn shaped<S>(self, build: impl FnOnce(&State) -> OgeomResult<S>) -> OgeomResult<Fit<S>> {
        Ok(Fit {
            shape: build(&self.state)?,
            sigma: self.sigma,
            form: self.form,
            used: self.used,
            residuals: self.residuals,
            iterations: self.iterations,
            rounds: self.rounds,
            converged: self.converged,
            side: self.side,
        })
    }
}

/// The most rounds a cutoff fits before giving up on the kept points
/// settling: each round costs a whole fit, and a set still changing after
/// this many is oscillating between two, reported as not converged.
const MOST_ROUNDS: u32 = 50;

/// The most points a start is estimated from, taken evenly through the
/// set: the start needs only to land where the least squares converges
/// from, which a few thousand points fix as well as millions, and this
/// bounds its cost.
const START_POINTS: usize = 4096;

/// The size of a Gaussian's median absolute residual in its standard
/// deviations, `Φ⁻¹(3/4)`: a robust start's spread is its median residual
/// divided by this.
const MEDIAN_IN_SIGMAS: f64 = 0.674_489_750_196_081_7;

fn run<S>(
    kind: Kind,
    points: &[Point],
    options: &FitOptions<'_, S>,
    guess: Option<State>,
    tol: Tolerances,
) -> OgeomResult<Outcome> {
    let normals = checked(kind, points, options)?;
    let fixed = &options.fixed;
    let needed = needed(kind, fixed);
    if points.len() < needed {
        ogeom_bail!(
            Construction,
            "a {} fit needs at least {needed} points, got {}",
            kind.name(),
            points.len()
        );
    }

    let stride = points.len().div_ceil(START_POINTS).max(1);
    let sample: Vec<Point> = points.iter().step_by(stride).copied().collect();
    let sample_normals: Option<Vec<Vector>> = match &normals {
        Some(n) => Some(n.iter().step_by(stride).copied().collect()),
        None if kind.revolved() && start_reads_normals(kind, fixed) => {
            Some(estimated_normals(&sample)?)
        }
        None => None,
    };
    let sample_normals = sample_normals.as_deref();
    let given = normals.is_some();

    let mut state = match guess {
        Some(g) if !kind.closed() => pinned(kind, g, fixed),
        _ if options.robust_start => consensus(kind, &sample, sample_normals, given, fixed, tol)?,
        _ => start(kind, &sample, sample_normals, given, fixed, tol)?,
    };

    let mut used = vec![true; points.len()];
    if options.robust_start && options.cutoff != Cutoff::All {
        let sizes: Vec<f64> = points.iter().map(|p| residual(kind, &state, *p)).collect();
        let threshold = robust_threshold(options.cutoff, &sizes);
        for (keep, r) in used.iter_mut().zip(&sizes) {
            *keep = r.abs() <= threshold;
        }
        refuse_too_few_kept(kind, &used, needed)?;
    }

    let mut iterations = 0_u32;
    let mut rounds = 0_u32;
    let mut converged;
    let mut residuals;
    loop {
        rounds += 1;
        let (next, steps, settled) = if kind.closed() {
            (
                closed_form(kind, &state, points, &used, fixed, tol)?,
                0,
                true,
            )
        } else {
            least_squares(kind, state, &free(kind, fixed), points, &used, tol)
        };
        state = next;
        iterations += steps;
        converged = settled;
        residuals = points
            .iter()
            .map(|p| residual(kind, &state, *p))
            .collect::<Vec<_>>();
        if options.cutoff == Cutoff::All {
            break;
        }
        let kept: Vec<f64> = residuals
            .iter()
            .zip(&used)
            .filter_map(|(r, keep)| keep.then_some(*r))
            .collect();
        let threshold = options.cutoff.threshold(&kept);
        let next_used: Vec<bool> = residuals.iter().map(|r| r.abs() <= threshold).collect();
        if next_used == used {
            break;
        }
        refuse_too_few_kept(kind, &next_used, needed)?;
        if rounds >= MOST_ROUNDS {
            converged = false;
            break;
        }
        used = next_used;
    }
    if kind == Kind::Cone && state.second < 0.0 {
        // The same cone seen from its other end.
        state.axis = -state.axis;
        state.second = -state.second;
    }

    let mut side = None;
    if kind.sided() {
        side = options.side.or_else(|| {
            normals.as_deref().and_then(|n| {
                let vote = facing(kind, &state, points, n, &used);
                if vote > 0.0 {
                    Some(Side::Inside)
                } else if vote < 0.0 {
                    Some(Side::Outside)
                } else {
                    None
                }
            })
        });
        if side == Some(Side::Outside) {
            for r in &mut residuals {
                *r = -*r;
            }
        }
    } else if kind == Kind::Plane
        && let Some(n) = normals.as_deref()
        && facing(kind, &state, points, n, &used) < 0.0
    {
        state.axis = -state.axis;
        for r in &mut residuals {
            *r = -*r;
        }
    }

    let (mut sum, mut count) = (0.0, 0.0);
    let mut form = (f64::INFINITY, f64::NEG_INFINITY);
    for (r, keep) in residuals.iter().zip(&used) {
        if *keep {
            sum += r * r;
            count += 1.0;
            form = (form.0.min(*r), form.1.max(*r));
        }
    }
    let sigma = (sum / count).sqrt();

    if let Some((lo, hi)) = options.radius_range
        && !(lo..=hi).contains(&state.radius)
    {
        ogeom_bail!(
            Range,
            "the fitted {} radius {} lies outside the range [{lo}, {hi}]",
            kind.name(),
            state.radius
        );
    }
    Ok(Outcome {
        state,
        sigma,
        form,
        used,
        residuals,
        iterations,
        rounds,
        converged,
        side,
    })
}

/// The options checked against the kind and the points, and the normals
/// made unit.
fn checked<S>(
    kind: Kind,
    points: &[Point],
    options: &FitOptions<'_, S>,
) -> OgeomResult<Option<Vec<Vector>>> {
    let name = kind.name();
    if let Some(i) = points.iter().position(|p| !p.is_finite()) {
        ogeom_bail!(Domain, "point {i} is not finite");
    }
    options.cutoff.check()?;
    let fixed = &options.fixed;
    let refuse = |what: &str| -> OgeomResult<()> {
        ogeom_bail!(Domain, "{what} does not apply to a {name} fit")
    };
    if fixed.radius.is_some() && !kind.radial() {
        refuse("a fixed radius")?;
    }
    if fixed.direction.is_some() && !kind.directed() {
        refuse("a fixed direction")?;
    }
    if fixed.axis.is_some() && !kind.axial() {
        refuse("a fixed axis")?;
    }
    if fixed.half_angle.is_some() && kind != Kind::Cone {
        refuse("a fixed half angle")?;
    }
    if options.side.is_some() && !kind.sided() {
        refuse("a side")?;
    }
    if options.radius_range.is_some() && !kind.radial() {
        refuse("a radius range")?;
    }
    if options.normals.is_some() && kind.planar() {
        refuse("normals")?;
    }
    if fixed.direction.is_some() && fixed.axis.is_some() {
        ogeom_bail!(
            Domain,
            "a fixed axis holds its direction; give the axis or the direction, not both"
        );
    }
    if let Some(r) = fixed.radius
        && !(r.is_finite() && r > 0.0)
    {
        ogeom_bail!(Domain, "the fixed radius {r} must be finite and positive");
    }
    if let Some(a) = fixed.half_angle
        && !(a.is_finite() && a > 0.0 && a < core::f64::consts::FRAC_PI_2)
    {
        ogeom_bail!(
            Domain,
            "the fixed half angle {a} must lie strictly between 0 and pi/2"
        );
    }
    if let Some((lo, hi)) = options.radius_range
        && !(lo.is_finite() && hi.is_finite() && 0.0 <= lo && lo <= hi)
    {
        ogeom_bail!(
            Domain,
            "the radius range [{lo}, {hi}] must be finite, non-negative and in order"
        );
    }
    let Some(normals) = options.normals else {
        return Ok(None);
    };
    if normals.len() != points.len() {
        ogeom_bail!(
            Dimension,
            "{} normals for {} points",
            normals.len(),
            points.len()
        );
    }
    normals
        .iter()
        .enumerate()
        .map(|(i, n)| {
            let m = n.magnitude();
            if !(m.is_finite() && m > 0.0) {
                ogeom_bail!(Domain, "normal {i} has no direction");
            }
            Ok(*n / m)
        })
        .collect::<OgeomResult<Vec<_>>>()
        .map(Some)
}

fn refuse_too_few_kept(kind: Kind, used: &[bool], needed: usize) -> OgeomResult<()> {
    let kept = used.iter().filter(|u| **u).count();
    if kept < needed {
        ogeom_bail!(
            Construction,
            "the cutoff keeps {kept} points, fewer than the {needed} a {} fit needs",
            kind.name()
        );
    }
    Ok(())
}

/// A cutoff's first threshold over a robust start's residuals: the spread
/// read from their median, which the outliers a consensus start ignores
/// cannot inflate the way they inflate an RMS.
fn robust_threshold(cutoff: Cutoff, residuals: &[f64]) -> f64 {
    let median = Cutoff::Median(1.0).threshold(residuals);
    match cutoff {
        Cutoff::Sigma(k) => k * median / MEDIAN_IN_SIGMAS,
        other => other.threshold(residuals),
    }
}

/// How the normals face the shape: the sum over the kept points of each
/// normal's component along the direction its residual grows.
fn facing(kind: Kind, state: &State, points: &[Point], normals: &[Vector], used: &[bool]) -> f64 {
    let mut sum = 0.0;
    for ((p, n), keep) in points.iter().zip(normals).zip(used) {
        if *keep {
            sum += n.dot(outward(kind, state, *p));
        }
    }
    sum
}

/// The direction a point's residual grows in.
fn outward(kind: Kind, state: &State, p: Point) -> Vector {
    let a = state.axis.vector();
    let w = p - state.origin;
    let h = w.dot(a);
    let q = w - a * h;
    let unit = |v: Vector| {
        let m = v.magnitude();
        if m > 0.0 { v / m } else { Vector::ZERO }
    };
    match kind {
        Kind::Plane => a,
        Kind::Sphere | Kind::Circle2 | Kind::Line3 => unit(w),
        Kind::Line2 => Vector::new(-a.y, a.x, 0.0),
        Kind::Cylinder | Kind::Circle3 => unit(q),
        Kind::Cone => {
            let (s, c) = state.second.sin_cos();
            unit(q) * c - a * s
        }
        Kind::Torus => {
            let rho = q.magnitude();
            let radial = unit(q);
            unit(radial * (rho - state.second) + a * h)
        }
    }
}

/// A point's signed distance to the shape, positive away from its centre,
/// its axis or its spine, along a plane's normal, to the left of a line in
/// the plane; a line in space's distance unsigned.
fn residual(kind: Kind, state: &State, p: Point) -> f64 {
    let a = state.axis.vector();
    let w = p - state.origin;
    let h = w.dot(a);
    let rho = (w - a * h).magnitude();
    match kind {
        Kind::Plane => h,
        Kind::Line3 => rho,
        Kind::Line2 => a.x.mul_add(w.y, -(a.y * w.x)),
        Kind::Sphere | Kind::Circle2 => w.magnitude() - state.radius,
        Kind::Cylinder => rho - state.radius,
        Kind::Cone => {
            let (s, c) = state.second.sin_cos();
            (rho - state.radius).mul_add(c, -(h * s))
        }
        Kind::Torus => (rho - state.second).hypot(h) - state.radius,
        Kind::Circle3 => {
            let across = rho - state.radius;
            across.hypot(h).copysign(across)
        }
    }
}

/// The directions a step moves the point along: two across the axis, and
/// the axis itself; the world's axes for a sphere and a circle in the
/// plane.
fn basis(kind: Kind, state: &State) -> [Vector; 3] {
    if matches!(kind, Kind::Sphere | Kind::Circle2) {
        return [Vector::X, Vector::Y, Vector::Z];
    }
    let a = state.axis;
    let e1 = a.any_perpendicular().vector();
    [e1, a.vector().cross(e1), a.vector()]
}

/// A point's distances to the shape (one, or two for a circle in space),
/// each with its derivative in every parameter.
fn distances(
    kind: Kind,
    state: &State,
    [e1, e2, a]: [Vector; 3],
    p: Point,
) -> ([(f64, [f64; PARAMETERS]); 2], usize) {
    let mut out = [(0.0, [0.0; PARAMETERS]); 2];
    let w = p - state.origin;
    if matches!(kind, Kind::Sphere | Kind::Circle2) {
        let length = w.magnitude();
        let u = if length > 0.0 { w / length } else { Vector::X };
        let d = &mut out[0].1;
        d[ACROSS[0]] = -u.x;
        d[ACROSS[1]] = -u.y;
        d[ALONG] = -u.z;
        d[RADIUS] = -1.0;
        out[0].0 = length - state.radius;
        return (out, 1);
    }
    let h = w.dot(a);
    let q = w - a * h;
    let rho = q.magnitude();
    let n = if rho > 0.0 { q / rho } else { e1 };
    let es = [e1, e2];
    // How the distance from the axis and the height along it move: the
    // point moved across by `e` takes the distance by `-n·e`, and the axis
    // tilted toward `e` takes it by `-h n·e` and the height by `w·e`.
    let rho_across = es.map(|e| -n.dot(e));
    let rho_tilt = es.map(|e| -h * n.dot(e));
    let h_tilt = es.map(|e| w.dot(e));
    match kind {
        Kind::Cylinder => {
            let d = &mut out[0].1;
            for i in 0..2 {
                d[TILT[i]] = rho_tilt[i];
                d[ACROSS[i]] = rho_across[i];
            }
            d[RADIUS] = -1.0;
            out[0].0 = rho - state.radius;
        }
        Kind::Cone => {
            let (s, c) = state.second.sin_cos();
            let d = &mut out[0].1;
            for i in 0..2 {
                d[TILT[i]] = c.mul_add(rho_tilt[i], -(s * h_tilt[i]));
                d[ACROSS[i]] = c * rho_across[i];
            }
            d[RADIUS] = -c;
            d[SECOND] = (state.radius - rho).mul_add(s, -(h * c));
            out[0].0 = (rho - state.radius).mul_add(c, -(h * s));
        }
        Kind::Torus => {
            let off = rho - state.second;
            let tube = off.hypot(h);
            let (m, k) = if tube > 0.0 {
                (off / tube, h / tube)
            } else {
                (1.0, 0.0)
            };
            let d = &mut out[0].1;
            for i in 0..2 {
                d[TILT[i]] = m.mul_add(rho_tilt[i], k * h_tilt[i]);
                d[ACROSS[i]] = m * rho_across[i];
            }
            d[ALONG] = -k;
            d[RADIUS] = -1.0;
            d[SECOND] = -m;
            out[0].0 = tube - state.radius;
        }
        Kind::Circle3 => {
            let (height, across) = out.split_at_mut(1);
            let dh = &mut height[0].1;
            let dr = &mut across[0].1;
            for i in 0..2 {
                dh[TILT[i]] = h_tilt[i];
                dr[TILT[i]] = rho_tilt[i];
                dr[ACROSS[i]] = rho_across[i];
            }
            dh[ALONG] = -1.0;
            dr[RADIUS] = -1.0;
            height[0].0 = h;
            across[0].0 = rho - state.radius;
            return (out, 2);
        }
        Kind::Plane | Kind::Line3 | Kind::Line2 | Kind::Sphere | Kind::Circle2 => {}
    }
    (out, 1)
}

/// The state moved by a step in the parameters, then its point slid along
/// the axis level with the kept points' centroid where nothing pins it
/// there (a cylinder's, a cone's), so the steps cannot wander along it.
fn stepped(
    kind: Kind,
    state: &State,
    [e1, e2, a]: [Vector; 3],
    step: &[f64; PARAMETERS],
    centroid: Point,
    tol: Tolerances,
) -> State {
    let mut next = *state;
    next.origin = state.origin + e1 * step[ACROSS[0]] + e2 * step[ACROSS[1]] + a * step[ALONG];
    if step[TILT[0]] != 0.0 || step[TILT[1]] != 0.0 {
        let tilted = a + e1 * step[TILT[0]] + e2 * step[TILT[1]];
        if let Ok(d) = Direction::new(tilted, tol) {
            next.axis = d;
        }
    }
    next.radius += step[RADIUS];
    next.second += step[SECOND];
    if matches!(kind, Kind::Cylinder | Kind::Cone) {
        let a = next.axis.vector();
        let shift = (centroid - next.origin).dot(a);
        next.origin += a * shift;
        if kind == Kind::Cone {
            // The radius is the cone's at its point.
            next.radius = shift.mul_add(next.second.tan(), next.radius);
        }
    }
    next
}

/// The most Gauss-Newton steps a round takes. Converging steps shrink
/// fast, by orders of magnitude each, from any start in the basin; a fit
/// still stepping after this many is reported as not converged.
const MOST_STEPS: u32 = 100;

/// Orthogonal least squares on the kept points, from `state`: Gauss-Newton
/// steps through the normal equations' eigenvectors, the directions the
/// points determine no better than the equations' rounding left where
/// they are, and each step halved until it lowers the cost. Converged once
/// a step moves the points' distances by less than the confusion distance.
fn least_squares(
    kind: Kind,
    mut state: State,
    free: &[bool; PARAMETERS],
    points: &[Point],
    used: &[bool],
    tol: Tolerances,
) -> (State, u32, bool) {
    let index: Vec<usize> = (0..PARAMETERS).filter(|k| free[*k]).collect();
    let size = index.len();
    if size == 0 {
        return (state, 0, true);
    }
    let kept: Vec<Point> = points
        .iter()
        .zip(used)
        .filter_map(|(p, keep)| keep.then_some(*p))
        .collect();
    let centroid = Point::from_vector(recognize::centroid(&kept));
    #[allow(
        clippy::cast_precision_loss,
        reason = "point counts are far below 2^52"
    )]
    let measured = (kept.len() * kind.components()) as f64;
    let cost = |s: &State| -> f64 {
        let frame = basis(kind, s);
        let mut sum = 0.0;
        for p in &kept {
            let (rows, count) = distances(kind, s, frame, *p);
            for (r, _) in &rows[..count] {
                sum += r * r;
            }
        }
        sum
    };
    let mut current = cost(&state);
    for step in 0..MOST_STEPS {
        let frame = basis(kind, &state);
        let mut normal = nalgebra::DMatrix::<f64>::zeros(size, size);
        let mut gradient = nalgebra::DVector::<f64>::zeros(size);
        for p in &kept {
            let (rows, count) = distances(kind, &state, frame, *p);
            for (r, d) in &rows[..count] {
                for (i, &ki) in index.iter().enumerate() {
                    gradient[i] += d[ki] * r;
                    for (j, &kj) in index.iter().enumerate().skip(i) {
                        normal[(i, j)] += d[ki] * d[kj];
                    }
                }
            }
        }
        for i in 0..size {
            for j in 0..i {
                normal[(i, j)] = normal[(j, i)];
            }
        }
        let eigen = nalgebra::SymmetricEigen::new(normal.clone());
        let largest = eigen.eigenvalues.iter().copied().fold(0.0, f64::max);
        #[allow(clippy::cast_precision_loss, reason = "a handful of parameters")]
        let floor = largest * f64::EPSILON * size as f64;
        let mut delta = nalgebra::DVector::<f64>::zeros(size);
        for (k, value) in eigen.eigenvalues.iter().enumerate() {
            if *value > floor {
                let v = eigen.eigenvectors.column(k);
                delta -= v * (v.dot(&gradient) / value);
            }
        }
        let movement = ((delta.dot(&(&normal * &delta))) / measured).sqrt();
        if !movement.is_finite() {
            return (state, step, false);
        }
        let mut scale = 1.0;
        let mut accepted = false;
        loop {
            let mut full = [0.0; PARAMETERS];
            for (i, &k) in index.iter().enumerate() {
                full[k] = delta[i] * scale;
            }
            let trial = stepped(kind, &state, frame, &full, centroid, tol);
            let c = cost(&trial);
            if c.is_finite() && c < current {
                state = trial;
                current = c;
                accepted = true;
                break;
            }
            if movement * scale <= tol.confusion() {
                break;
            }
            scale *= 0.5;
        }
        // A step too small to move a distance measurably, or one no
        // shortening of which lowers the cost, leaves the fit at its least
        // to rounding.
        if movement <= tol.confusion() || !accepted {
            return (state, step + 1, true);
        }
    }
    (state, MOST_STEPS, false)
}

/// A plane or a line over the kept points, in closed form: the centroid,
/// and the covariance's smallest direction (a plane's normal) or largest
/// (a line's direction), or the held direction. The sense follows the
/// last round's.
fn closed_form(
    kind: Kind,
    last: &State,
    points: &[Point],
    used: &[bool],
    fixed: &FitFixed,
    tol: Tolerances,
) -> OgeomResult<State> {
    let kept: Vec<Point> = points
        .iter()
        .zip(used)
        .filter_map(|(p, keep)| keep.then_some(*p))
        .collect();
    let c = recognize::centroid(&kept);
    let origin = Point::from_vector(c);
    let axis = if let Some(d) = fixed.direction {
        d
    } else {
        let m = recognize::covariance(kept.iter().map(|p| p.to_vector() - c));
        let eigen = nalgebra::SymmetricEigen::new(m);
        let mut order = [0, 1, 2];
        order.sort_by(|i, j| eigen.eigenvalues[*i].total_cmp(&eigen.eigenvalues[*j]));
        #[allow(
            clippy::cast_precision_loss,
            reason = "point counts are far below 2^52"
        )]
        let count = kept.len() as f64;
        let spread = |k: usize| (eigen.eigenvalues[order[k]].max(0.0) / count).sqrt();
        let (pick, flat, what) = if kind == Kind::Plane {
            (order[0], spread(1), "one line")
        } else {
            (order[2], spread(2), "one point")
        };
        if flat <= tol.confusion() {
            ogeom_bail!(
                Construction,
                "the points' spread is degenerate for a {} fit: they lie within {} of {what}",
                kind.name(),
                tol.confusion()
            );
        }
        let v = eigen.eigenvectors.column(pick);
        Direction::new(Vector::new(v[0], v[1], v[2]), tol)?
    };
    let axis = if axis.vector().dot(last.axis.vector()) < 0.0 {
        -axis
    } else {
        axis
    };
    Ok(State::about(origin, axis))
}

/// Whether a surface of revolution's start reads normals: one whose axis
/// is not given whole, and a cylinder only where not even its direction is.
fn start_reads_normals(kind: Kind, fixed: &FitFixed) -> bool {
    fixed.axis.is_none() && !(kind == Kind::Cylinder && fixed.direction.is_some())
}

/// How many neighbours a point's normal is estimated from, by the plane
/// through them: a few more than the three a plane needs, for noise to
/// average, and few enough to stay within a patch's curvature.
const NORMAL_NEIGHBOURS: usize = 12;

/// Each point's normal, in no particular sense, as the smallest direction
/// of its nearest neighbours' covariance.
fn estimated_normals(points: &[Point]) -> OgeomResult<Vec<Vector>> {
    let index = ogeom_mesh::PointIndex::build(points)?;
    Ok(points
        .iter()
        .map(|p| {
            let near: Vec<Point> = index
                .k_nearest(*p, NORMAL_NEIGHBOURS, f64::INFINITY)
                .iter()
                .map(|n| points[n.index as usize])
                .collect();
            let c = recognize::centroid(&near);
            recognize::smallest_direction(recognize::covariance(
                near.iter().map(|q| q.to_vector() - c),
            ))
        })
        .collect())
}

/// The state with what the caller holds put in: the radius, the half
/// angle, and the axis or its direction, in whichever sense the state's
/// own runs.
fn pinned(kind: Kind, mut state: State, fixed: &FitFixed) -> State {
    if let Some(r) = fixed.radius {
        state.radius = r;
    }
    if let Some(a) = fixed.half_angle {
        state.second = a;
    }
    let along = |d: Direction, current: Direction| {
        if d.vector().dot(current.vector()) < 0.0 {
            -d
        } else {
            d
        }
    };
    if let Some(d) = fixed.direction {
        state.axis = along(d, state.axis);
    }
    if let Some(axis) = fixed.axis {
        let d = along(axis.direction, state.axis);
        let foot = axis.location + d.vector() * (state.origin - axis.location).dot(d.vector());
        state.axis = d;
        state.origin = foot;
    }
    if kind == Kind::Cone && state.second < 0.0 {
        state.axis = -state.axis;
        state.second = -state.second;
    }
    state
}

fn degenerate(kind: Kind, why: &str) -> ogeom_core::OgeomError {
    ogeom_core::ogeom_err!(
        Construction,
        "the points' spread is degenerate for a {} fit: {why}",
        kind.name()
    )
}

/// The start a fit refines, from points and (for a surface of revolution)
/// their normals: `given` where the caller's, else estimated.
fn start(
    kind: Kind,
    points: &[Point],
    normals: Option<&[Vector]>,
    given: bool,
    fixed: &FitFixed,
    tol: Tolerances,
) -> OgeomResult<State> {
    let c = Point::from_vector(recognize::centroid(points));
    let state = match kind {
        Kind::Plane | Kind::Line3 | Kind::Line2 => {
            let all = vec![true; points.len()];
            let seed = State::about(c, Direction::Z);
            let seed = closed_form(kind, &seed, points, &all, fixed, tol)?;
            let v = seed.axis.vector();
            // A sense of the covariance's own choosing would leave the
            // answer to the eigen solver; the largest component is made
            // positive instead.
            let largest = [v.x, v.y, v.z]
                .into_iter()
                .fold(0.0_f64, |m, x| if x.abs() > m.abs() { x } else { m });
            if largest < 0.0 {
                State::about(seed.origin, -seed.axis)
            } else {
                seed
            }
        }
        Kind::Sphere => match recognize::fit_sphere(points, tol) {
            Some(Canonical::Sphere(s)) => State {
                radius: s.radius(),
                ..State::about(s.centre(), Direction::Z)
            },
            _ => return Err(degenerate(kind, "they lie on one plane")),
        },
        Kind::Circle2 => {
            let (centre, radius) = circle_in_plane(points, c, Direction::Z, tol)
                .ok_or_else(|| degenerate(kind, "they lie on one line"))?;
            State {
                radius,
                ..State::about(centre, Direction::Z)
            }
        }
        Kind::Circle3 => {
            if let Some(axis) = fixed.axis {
                let foot = on_line(axis, c);
                let radius = mean_distance(points, axis);
                State {
                    radius,
                    ..State::about(foot, axis.direction)
                }
            } else if let Some(d) = fixed.direction {
                let (centre, radius) = circle_in_plane(points, c, d, tol)
                    .ok_or_else(|| degenerate(kind, "they lie on one line"))?;
                State {
                    radius,
                    ..State::about(centre, d)
                }
            } else {
                let (centre, normal, radius, _) = recognize::circle_through(points, tol)
                    .ok_or_else(|| degenerate(kind, "they lie on one line"))?;
                State {
                    radius,
                    ..State::about(centre, normal)
                }
            }
        }
        Kind::Cylinder => {
            if let Some(axis) = fixed.axis {
                State {
                    radius: mean_distance(points, axis),
                    ..State::about(on_line(axis, c), axis.direction)
                }
            } else if let Some(d) = fixed.direction {
                let (centre, radius) = circle_in_plane(points, c, d, tol)
                    .ok_or_else(|| degenerate(kind, "their shadow along the axis is a line"))?;
                State {
                    radius,
                    ..State::about(centre, d)
                }
            } else {
                let normals = normals.ok_or_else(|| degenerate(kind, "no normals"))?;
                match recognize::fit_cylinder(points, normals, tol) {
                    Some(Canonical::Cylinder(cyl)) => State {
                        radius: cyl.radius(),
                        ..State::about(cyl.frame().origin(), cyl.frame().z())
                    },
                    _ => return Err(degenerate(kind, "their normals meet no one axis")),
                }
            }
        }
        Kind::Cone => {
            let (through, axis) = revolution_axis(kind, points, normals, fixed, tol)?;
            match recognize::ruled_about(points, through, axis, 0.0, tol) {
                Some(Canonical::Cone(cone)) => State {
                    radius: cone.reference_radius(),
                    second: cone.half_angle(),
                    ..State::about(cone.frame().origin(), cone.frame().z())
                },
                Some(_) => {
                    return Err(degenerate(
                        kind,
                        "they have no taper, and lie on a cylinder",
                    ));
                }
                None => return Err(degenerate(kind, "they span no height along the axis")),
            }
        }
        Kind::Torus => torus_start(points, normals, given, fixed, tol)?,
    };
    Ok(pinned(kind, state, fixed))
}

/// A cone's or a torus's axis: the one held, or the line every normal line
/// meets, through the held direction where that is held.
fn revolution_axis(
    kind: Kind,
    points: &[Point],
    normals: Option<&[Vector]>,
    fixed: &FitFixed,
    tol: Tolerances,
) -> OgeomResult<(Point, Direction)> {
    if let Some(axis) = fixed.axis {
        return Ok((axis.location, axis.direction));
    }
    let normals = normals.ok_or_else(|| degenerate(kind, "no normals"))?;
    let found = match fixed.direction {
        Some(d) => axis_along(points, normals, d).map(|p| (p, d)),
        None => recognize::revolution_axis(points, normals, tol),
    };
    found.ok_or_else(|| degenerate(kind, "their normals meet no one axis"))
}

/// The point where the line of direction `d` that every normal line meets
/// best crosses the plane through the origin across `d`: the line through
/// `x` meets the normal line through `p` along `n` where
/// `(p − x)·(d × n) = 0`, linear in `x`.
fn axis_along(points: &[Point], normals: &[Vector], d: Direction) -> Option<Point> {
    let reference = *points.first()?;
    let mut m = nalgebra::Matrix3::<f64>::zeros();
    let mut b = nalgebra::Vector3::<f64>::zeros();
    for (p, n) in points.iter().zip(normals) {
        let row = d.vector().cross(*n);
        let r = nalgebra::Vector3::new(row.x, row.y, row.z);
        m += r * r.transpose();
        b += r * (*p - reference).dot(row);
    }
    let scale = m.norm().max(f64::MIN_POSITIVE);
    let x = m.pseudo_inverse(scale * f64::EPSILON * 3.0).ok()? * b;
    Some(reference + Vector::new(x[0], x[1], x[2]))
}

/// A torus's start: its axis, then the circle the points trace in the
/// half-plane through it; and where the caller's own normals are given and
/// the axis is free, also the spine the normals point back to, whichever
/// fits the points closer.
fn torus_start(
    points: &[Point],
    normals: Option<&[Vector]>,
    given: bool,
    fixed: &FitFixed,
    tol: Tolerances,
) -> OgeomResult<State> {
    let kind = Kind::Torus;
    let c = Point::from_vector(recognize::centroid(points));
    let mut best: Option<(f64, State)> = None;
    let mut consider = |state: State| {
        let rms: f64 = points
            .iter()
            .map(|p| residual(kind, &state, *p).powi(2))
            .sum();
        if rms.is_finite() && best.is_none_or(|(held, _)| rms < held) {
            best = Some((rms, state));
        }
    };
    if let Ok((through, axis)) = revolution_axis(kind, points, normals, fixed, tol) {
        let line = Axis::new(through, axis);
        let base = on_line(line, c);
        if let Some((major, height, minor)) = recognize::profile_circle(points, base, line) {
            consider(State {
                radius: minor,
                second: major,
                ..State::about(base + axis.vector() * height, axis)
            });
        }
    }
    if given
        && fixed.axis.is_none()
        && fixed.direction.is_none()
        && let Some(n) = normals
        && let Some(Canonical::Torus(t)) = recognize::fit_torus(points, n, tol)
    {
        consider(State {
            radius: t.minor_radius(),
            second: t.major_radius(),
            ..State::about(t.centre(), t.frame().z())
        });
    }
    best.map(|(_, s)| s)
        .ok_or_else(|| degenerate(kind, "they trace no circle about any axis"))
}

/// The circle through points' shadows on the plane through `through`
/// across `normal`, by the algebraic fit: its centre, on that plane, and
/// its radius.
fn circle_in_plane(
    points: &[Point],
    through: Point,
    normal: Direction,
    tol: Tolerances,
) -> Option<(Point, f64)> {
    let a = normal.vector();
    let e1 = normal.any_perpendicular().vector();
    let e2 = a.cross(e1);
    let mut m = nalgebra::Matrix3::zeros();
    let mut b = nalgebra::Vector3::zeros();
    for p in points {
        let d = *p - through;
        let (x, y) = (d.dot(e1), d.dot(e2));
        let row = nalgebra::Vector3::new(-2.0 * x, -2.0 * y, 1.0);
        m += row * row.transpose();
        b += row * -x.mul_add(x, y * y);
    }
    let solved = m.lu().solve(&b)?;
    let r2 = solved[0].mul_add(solved[0], solved[1] * solved[1]) - solved[2];
    (r2 > tol.confusion() * tol.confusion())
        .then(|| (through + e1 * solved[0] + e2 * solved[1], r2.sqrt()))
}

/// The point of the axis nearest `p`.
fn on_line(axis: Axis, p: Point) -> Point {
    let d = axis.direction.vector();
    axis.location + d * (p - axis.location).dot(d)
}

/// The points' mean distance from an axis.
fn mean_distance(points: &[Point], axis: Axis) -> f64 {
    let d = axis.direction.vector();
    let sum: f64 = points
        .iter()
        .map(|p| {
            let w = *p - axis.location;
            (w - d * w.dot(d)).magnitude()
        })
        .sum();
    #[allow(
        clippy::cast_precision_loss,
        reason = "point counts are far below 2^52"
    )]
    let n = points.len() as f64;
    sum / n
}

/// The probability that a consensus start draws at least one subset free
/// of outliers is at least this, where half the points are outliers.
const CONSENSUS_MISS: f64 = 1e-3;

/// A start by consensus: the kind's start on many small subsets of the
/// points, drawn in a fixed sequence, and the one whose median residual
/// over all of them is least. The draws are as many as make a subset with
/// no outlier all but certain, with half the points outliers.
fn consensus(
    kind: Kind,
    points: &[Point],
    normals: Option<&[Vector]>,
    given: bool,
    fixed: &FitFixed,
    tol: Tolerances,
) -> OgeomResult<State> {
    let size = needed(kind, &FitFixed::default()).max(needed(kind, fixed));
    if points.len() <= size {
        return start(kind, points, normals, given, fixed, tol);
    }
    let clean = 0.5_f64.powi(i32::try_from(size).unwrap_or(i32::MAX));
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "a positive count of a few thousand"
    )]
    let draws = (CONSENSUS_MISS.ln() / (-clean).ln_1p()).ceil() as usize;
    let mut sequence = Draws(0x9E37_79B9_7F4A_7C15);
    let mut best: Option<(f64, State)> = None;
    let mut picked: Vec<usize> = Vec::with_capacity(size);
    let mut sizes = vec![0.0; points.len()];
    for _ in 0..draws {
        picked.clear();
        while picked.len() < size {
            let i = sequence.below(points.len());
            if !picked.contains(&i) {
                picked.push(i);
            }
        }
        let subset: Vec<Point> = picked.iter().map(|i| points[*i]).collect();
        let subset_normals: Option<Vec<Vector>> =
            normals.map(|n| picked.iter().map(|i| n[*i]).collect());
        let Ok(candidate) = start(kind, &subset, subset_normals.as_deref(), given, fixed, tol)
        else {
            continue;
        };
        for (s, p) in sizes.iter_mut().zip(points) {
            *s = residual(kind, &candidate, *p);
        }
        let median = Cutoff::Median(1.0).threshold(&sizes);
        if median.is_finite() && best.is_none_or(|(held, _)| median < held) {
            best = Some((median, candidate));
        }
    }
    match best {
        Some((_, state)) => Ok(state),
        None => start(kind, points, normals, given, fixed, tol),
    }
}

/// A fixed sequence of draws (xorshift64*), so a consensus start is the
/// same on every run.
struct Draws(u64);

impl Draws {
    /// A draw below `n`.
    fn below(&mut self, n: usize) -> usize {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        let wide = u128::from(x.wrapping_mul(0x2545_F491_4F6C_DD1D)) * n as u128;
        usize::try_from(wide >> 64).unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;

    const T: Tolerances = Tolerances::millimetres();

    /// Every kind's derivatives in the parameters it fits are those of its
    /// distances: each against a central difference through the same step the least squares takes,
    /// at points off the shape on every side.
    #[test]
    fn the_derivatives_are_the_distances_own() {
        let axis = Direction::new(Vector::new(0.2, -0.3, 1.0), T).unwrap();
        let state = State {
            origin: Point::new(1.0, 2.0, -1.0),
            axis,
            radius: 3.0,
            second: 0.4,
        };
        let torus = State {
            second: 9.0,
            ..state
        };
        let points = [
            Point::new(4.0, 1.0, 2.0),
            Point::new(-3.0, 5.0, 0.5),
            Point::new(2.5, -6.0, -4.0),
            Point::new(9.0, 3.0, 1.0),
        ];
        for (kind, at) in [
            (Kind::Sphere, state),
            (Kind::Circle3, state),
            (Kind::Cylinder, state),
            (Kind::Cone, state),
            (Kind::Torus, torus),
        ] {
            let frame = basis(kind, &at);
            for p in points {
                let (rows, count) = distances(kind, &at, frame, p);
                let on = free(kind, &FitFixed::default());
                for k in (0..PARAMETERS).filter(|k| on[*k]) {
                    let h = 1e-6;
                    let moved = |sign: f64| {
                        let mut step = [0.0; PARAMETERS];
                        step[k] = sign * h;
                        // Slid along the axis to where the point already
                        // stands, so only the step moves it.
                        let s = stepped(kind, &at, frame, &step, at.origin, T);
                        distances(kind, &s, basis(kind, &s), p).0
                    };
                    let (up, down) = (moved(1.0), moved(-1.0));
                    for c in 0..count {
                        let numeric = (up[c].0 - down[c].0) / (2.0 * h);
                        let analytic = rows[c].1[k];
                        assert!(
                            (numeric - analytic).abs() < 1e-6,
                            "{kind:?} parameter {k} component {c}: {analytic} against {numeric}"
                        );
                    }
                }
            }
        }
    }

    /// A plane's leverage at a point is the closed form of linear least
    /// squares: `1/n` for its offset and `wᵀ S⁻¹ w` for its tilt, `w` the
    /// point's offset from the centroid across the plane and `S` the
    /// points' scatter there.
    #[test]
    #[allow(clippy::cast_precision_loss, reason = "a few points")]
    fn a_planes_leverage_is_the_closed_form() {
        let points: Vec<Point> = (0..7)
            .flat_map(|i| (0..5).map(move |j| (i, j)))
            .map(|(i, j)| {
                let (x, y) = (f64::from(i) * 0.7, f64::from(j) * 1.3 - 0.4);
                Point::new(x, y, 0.01 * f64::from((i * 3 + j * 5) % 7) - 0.03)
            })
            .collect();
        let normals = vec![Vector::Z; points.len()];
        let options = FitOptions {
            normals: Some(&normals),
            ..FitOptions::default()
        };
        let fit = plane(&points, &options, T).unwrap();
        let model = Model::plane(&fit).measured_over(&points, &fit.used, &FitFixed::default());
        assert_eq!(model.parameters(), 3);

        let n = fit.shape.normal().vector();
        let e1 = fit.shape.normal().any_perpendicular().vector();
        let e2 = n.cross(e1);
        let count = points.len() as f64;
        let centroid = points
            .iter()
            .fold(Vector::ZERO, |sum, p| sum + p.to_vector())
            / count;
        let across = |p: Point| {
            let w = p.to_vector() - centroid;
            (w.dot(e1), w.dot(e2))
        };
        let (mut a, mut b, mut c) = (0.0, 0.0, 0.0);
        for p in &points {
            let (u, v) = across(*p);
            a += u * u;
            b += u * v;
            c += v * v;
        }
        let det = a.mul_add(c, -(b * b));
        for p in [Point::new(9.0, -4.0, 0.2), Point::new(2.1, 2.2, 0.0)] {
            let (u, v) = across(p);
            let closed = 1.0 / count + (c * u * u - 2.0 * b * u * v + a * v * v) / det;
            let got = model.leverage(p);
            assert!(
                (got - closed).abs() <= 1e-12 * closed,
                "{got} against {closed}"
            );
        }
    }
}
