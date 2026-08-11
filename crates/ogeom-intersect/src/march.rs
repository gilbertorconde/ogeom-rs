//! The general surface/surface intersector: seed, then walk.
//!
//! Where two surfaces meet has no closed form in general — the curve is
//! transcendental, and `docs/DATA_MODEL.md` §9 is blunt about the consequence:
//! there is no exact answer to be exact about, which is why the topology carries
//! tolerances. What there is instead is a curve that can be *followed*, one
//! corrected step at a time, to a stated accuracy.
//!
//! # Two problems, kept apart on purpose
//!
//! **Finding a branch** and **following one** fail in completely different ways,
//! and lumping them together is how an intersector comes to look better than it
//! is. A tracer that follows one branch beautifully while never noticing the
//! second reports a smooth, accurate, *wrong* answer — and the obvious accuracy
//! measure, "is every point on both surfaces", scores it perfectly.
//!
//! So [`seeds`] and [`trace`] are separate, separately testable, and separately
//! measured. Seeding is polyhedral: both surfaces are sampled into triangles and
//! the triangle pairs that cross give starting points. It finds a branch if the
//! sampling resolves it, and *misses one thinner than the grid* — which is a
//! real limitation with a knob attached rather than a mystery.
//!
//! # Following the curve
//!
//! At a point on both surfaces the intersection runs along the cross product of
//! the two normals: the one direction that stays in both tangent planes. Step
//! along it and you leave both surfaces slightly; a Newton correction brings you
//! back.
//!
//! The correction has four unknowns — two parameters on each surface — and three
//! equations, `A(u1,v1) = B(u2,v2)`. That is deliberately one short, because the
//! solution set *is* the curve and pinning it to a point needs one more
//! condition. The fourth is a plane across the direction of travel: it says how
//! far along to land, and it is what turns "somewhere on the curve" into "the
//! next point".
//!
//! # What it reports about itself
//!
//! Whether the curve closed, and whether it ran out of steps. A polyline that
//! stopped because it hit a limit is not the same answer as one that stopped
//! because the curve ended, and a caller that cannot tell them apart will treat
//! a truncated branch as a complete one.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Surface, SurfaceGeometry};
use ogeom_math::{Point, Vector, solve};

/// How hard to look, and how closely to follow.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Marching {
    /// How far the polyline may sit from the true curve, in space.
    pub chord: f64,
    /// How finely each surface is sampled when looking for branches.
    ///
    /// The limitation with a knob on it: a branch narrower than one cell can be
    /// stepped over entirely. Raising this costs time quadratically and is the
    /// only thing that makes a thin branch findable.
    pub grid: usize,
    /// A ceiling on the points in one branch, so a curve that will not close
    /// cannot run forever.
    pub max_points: usize,
}

impl Default for Marching {
    fn default() -> Self {
        Self {
            chord: 1e-4,
            grid: 24,
            max_points: 20_000,
        }
    }
}

impl Marching {
    /// Check the settings are usable.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the chord is
    /// not positive, the grid is too coarse to hold a triangle, or no points are
    /// allowed.
    pub fn validate(&self) -> OgeomResult<()> {
        if !self.chord.is_finite() || self.chord <= 0.0 {
            ogeom_bail!(Construction, "a chord of {} is not a distance", self.chord);
        }
        if self.grid < 2 {
            ogeom_bail!(Construction, "a sampling grid needs at least two steps");
        }
        if self.max_points < 2 {
            ogeom_bail!(Construction, "a branch needs at least two points");
        }
        Ok(())
    }
}

/// A point that lies on both surfaces, with where it is on each.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Contact {
    /// Parameters on the first surface.
    pub on_a: (f64, f64),
    /// Parameters on the second.
    pub on_b: (f64, f64),
    /// Where that is in space.
    pub point: Point,
}

/// Why a traced branch stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stopped {
    /// It came back to where it started.
    Closed,
    /// It reached the edge of one surface's domain.
    LeftTheDomain,
    /// The correction stopped converging — a tangency or a singular point.
    ///
    /// Reported rather than pushed through. Marching past a point where the two
    /// normals are parallel is how a tracer jumps onto the wrong branch, and a
    /// wrong branch is a plausible answer to a different question.
    Stalled,
    /// It hit [`Marching::max_points`].
    ///
    /// Distinct from every other reason, because this one means the answer is
    /// *incomplete* rather than finished.
    RanOut,
}

/// One traced branch.
#[derive(Debug, Clone, PartialEq)]
pub struct Traced {
    /// The points along it, in order.
    pub points: Vec<Point>,
    /// Where each point is on the first surface.
    pub on_a: Vec<(f64, f64)>,
    /// Where each point is on the second.
    pub on_b: Vec<(f64, f64)>,
    /// Why it stopped.
    pub stopped: Stopped,
}

impl Traced {
    /// Whether the branch is finished rather than truncated.
    #[must_use]
    pub const fn complete(&self) -> bool {
        !matches!(self.stopped, Stopped::RanOut)
    }

    /// Whether it closed on itself.
    #[must_use]
    pub const fn closed(&self) -> bool {
        matches!(self.stopped, Stopped::Closed)
    }
}

/// Starting points on the intersection, one per branch found.
///
/// Polyhedral: both surfaces are sampled into triangles, the pairs that cross
/// give approximate points, and each is corrected onto both surfaces exactly.
/// Points that land on the same spot are merged, so a branch crossing many
/// cells yields one seed rather than dozens.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the settings are
/// unusable.
pub fn seeds(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Vec<Contact>> {
    options.validate()?;
    let (mesh_a, mesh_b) = (sample(a, options.grid, tol), sample(b, options.grid, tol));

    let mut found: Vec<Contact> = Vec::new();
    for cell_a in &mesh_a {
        for cell_b in &mesh_b {
            // Cheap rejection first: most pairs are nowhere near each other,
            // and the segment test below is far from free.
            if !overlap(cell_a, cell_b, options.chord) {
                continue;
            }
            let Some(guess) = triangles_cross(cell_a, cell_b) else {
                continue;
            };
            let start = [cell_a.at.0, cell_a.at.1, cell_b.at.0, cell_b.at.1];
            let Some(contact) = correct(a, b, start, guess, None, tol) else {
                continue;
            };
            // One seed per branch, not one per cell it passes through. The
            // spacing is the grid's, since two distinct branches closer than
            // that were never going to be told apart by this sampling anyway.
            let apart = span(a).max(span(b)) / f64::from(u32::try_from(options.grid).unwrap_or(1));
            if found
                .iter()
                .any(|c| c.point.distance(contact.point) <= apart)
            {
                continue;
            }
            found.push(contact);
        }
    }
    Ok(found)
}

/// Every branch of the intersection: seed, trace each, and keep the distinct
/// ones.
///
/// A branch crossing many sampling cells produces many seeds, and tracing from
/// any of them gives the same curve. So a seed already lying on something
/// traced is dropped rather than followed again — which is what makes the
/// *number* of branches returned meaningful, and it is the number a boolean
/// will act on.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the settings are
/// unusable.
pub fn branches(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Vec<Traced>> {
    let found = seeds(a, b, options, tol)?;
    let mut out: Vec<Traced> = Vec::new();
    for seed in found {
        // Already on something we have followed.
        let reach = options.chord.max(tol.confusion()) * 8.0;
        if out
            .iter()
            .any(|branch| passes_near(branch, seed.point, reach))
        {
            continue;
        }
        // A seed that will not trace — a tangency — is reported by being
        // absent rather than by an error, since the other branches are still
        // real answers.
        if let Ok(branch) = trace(a, b, seed, options, tol)
            && branch.points.len() >= 2
            && !is_fragment(&branch, options)
        {
            out.push(branch);
        }
    }
    Ok(stitch_stalled(out, a, b, options, tol))
}

/// Below this sine the surfaces count as tangent at a point — the
/// branch-point certificate a stall end must carry to participate in
/// stitching.
const BRANCH_POINT_SINE: f64 = 0.05;

/// The sine of the normal angle at a contact — the transversality measure.
fn crossing_sine(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    on_a: (f64, f64),
    on_b: (f64, f64),
    tol: Tolerances,
) -> f64 {
    let Ok(na) = a.normal_at(on_a.0, on_a.1, tol) else {
        return 0.0;
    };
    let Ok(nb) = b.normal_at(on_b.0, on_b.1, tol) else {
        return 0.0;
    };
    na.vector().cross(nb.vector()).magnitude()
}

/// Whether a stalled trace is a fragment rather than a curve.
///
/// Coincident or near-coincident surfaces defeat the tangency check at a seed —
/// rounding in the corrected parameters leaves the two normals a whisker apart,
/// the walk takes a couple of steps, and then stalls where the arithmetic gives
/// out. What comes back lies on both surfaces perfectly and describes nothing:
/// identical spheres yielded six such fragments, each a few points long.
///
/// A stalled branch shorter than a handful of chords carries no information the
/// seed did not, so it is noise from a degenerate configuration and dropped. A
/// *real* stalled branch — one that ran into a genuine tangency — has length
/// behind it and is kept, because a truncated real answer is still an answer.
///
/// The marcher is deliberately not a coincidence detector: for the pairs with
/// closed forms, [`surface_surface`](crate::surface_surface) answers
/// [`Same`](crate::Meeting::Same), and that check belongs before this one.
fn is_fragment(branch: &Traced, options: Marching) -> bool {
    if branch.stopped != Stopped::Stalled {
        return false;
    }
    let length: f64 = branch
        .points
        .windows(2)
        .map(|pair| pair[0].distance(pair[1]))
        .sum();
    length < options.chord * 10.0
}

/// Whether a traced branch passes within a distance of a point.
///
/// Measured against the polyline's *segments*, not its vertices. The vertices
/// are a marching step apart — far more than the chord tolerance — so a seed
/// sitting neatly between two of them looks distant from both, and comparing to
/// vertices alone reported one circle nine times.
fn passes_near(branch: &Traced, p: Point, reach: f64) -> bool {
    branch
        .points
        .windows(2)
        .any(|pair| distance_to_segment(p, pair[0], pair[1]) <= reach)
}

/// Distance from a point to a segment.
fn distance_to_segment(p: Point, a: Point, b: Point) -> f64 {
    let along = b - a;
    let length = along.square_magnitude();
    if length <= f64::MIN_POSITIVE {
        return p.distance(a);
    }
    let t = ((p - a).dot(along) / length).clamp(0.0, 1.0);
    p.distance(a + along * t)
}

/// Follow the intersection from a starting point, in both directions.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the settings are
/// unusable; [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the two surfaces
/// are tangent at the seed, where there is no single direction to follow.
pub fn trace(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    from: Contact,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Traced> {
    options.validate()?;
    if tangent_at(a, b, from, tol).is_none() {
        ogeom_bail!(
            NotDone,
            "the surfaces are tangent here, so the intersection has no single \
             direction to follow; that is a branch point and needs the seed \
             moved off it"
        );
    }

    // Forwards first. If it closes, that is the whole branch and there is
    // nothing behind us.
    let ahead = walk(a, b, from, 1.0, options, tol)?;
    if ahead.stopped == Stopped::Closed {
        return Ok(ahead);
    }
    let behind = walk(a, b, from, -1.0, options, tol)?;

    // Join them, the backward half reversed and its shared first point dropped.
    let mut points = behind.points;
    let mut on_a = behind.on_a;
    let mut on_b = behind.on_b;
    points.reverse();
    on_a.reverse();
    on_b.reverse();
    points.pop();
    on_a.pop();
    on_b.pop();
    points.extend(ahead.points);
    on_a.extend(ahead.on_a);
    on_b.extend(ahead.on_b);

    // The worse of the two reasons: a branch truncated at either end is
    // truncated.
    let stopped = if ahead.stopped == Stopped::RanOut || behind.stopped == Stopped::RanOut {
        Stopped::RanOut
    } else if ahead.stopped == Stopped::Stalled || behind.stopped == Stopped::Stalled {
        Stopped::Stalled
    } else {
        Stopped::LeftTheDomain
    };
    Ok(Traced {
        points,
        on_a,
        on_b,
        stopped,
    })
}

/// Two surfaces, as a condition for the walker: four unknowns and three
/// equations saying the two points coincide.
///
/// The intersector's own walk goes through [`crate::walk`] like everything
/// else, and what is *not* generic lives here — the domain clamps a surface
/// pair needs, and the tangent, which the intersector computes from the two
/// normals rather than from the null space so that it can refuse a crossing
/// too shallow to be more than the correction's own noise.
struct SurfacePair<'s> {
    a: &'s SurfaceGeometry,
    b: &'s SurfaceGeometry,
}

impl crate::walk::Condition for SurfacePair<'_> {
    fn unknowns(&self) -> usize {
        4
    }

    fn position(&self, x: &[f64], tol: Tolerances) -> Option<Point> {
        self.a.point_at(x[0], x[1], tol).ok()
    }

    fn position_gradient(&self, x: &[f64], tol: Tolerances) -> Option<Vec<Vector>> {
        let (au, av) = self.a.d1_at(x[0], x[1], tol).ok()?;
        // The point is taken from the first surface, so it does not move with
        // the second's parameters at all.
        Some(vec![au, av, Vector::ZERO, Vector::ZERO])
    }

    fn system(&self, x: &[f64], tol: Tolerances) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
        let pa = self.a.point_at(x[0], x[1], tol).ok()?;
        let pb = self.b.point_at(x[2], x[3], tol).ok()?;
        let (au, av) = self.a.d1_at(x[0], x[1], tol).ok()?;
        let (bu, bv) = self.b.d1_at(x[2], x[3], tol).ok()?;
        let gap = pa - pb;
        Some((
            vec![gap.x, gap.y, gap.z],
            vec![
                vec![au.x, av.x, -bu.x, -bv.x],
                vec![au.y, av.y, -bu.y, -bv.y],
                vec![au.z, av.z, -bu.z, -bv.z],
            ],
        ))
    }

    fn clamp(&self, x: &mut [f64]) {
        let (ua, va) = clamp(self.a, x[0], x[1]);
        let (ub, vb) = clamp(self.b, x[2], x[3]);
        x[0] = ua;
        x[1] = va;
        x[2] = ub;
        x[3] = vb;
    }

    fn outside(&self, x: &[f64], tol: Tolerances) -> bool {
        outside(self.a, (x[0], x[1]), tol) || outside(self.b, (x[2], x[3]), tol)
    }

    fn near_edge(&self, x: &[f64]) -> bool {
        near_edge(self.a, (x[0], x[1])) || near_edge(self.b, (x[2], x[3]))
    }

    fn extent(&self) -> f64 {
        span(self.a).max(span(self.b))
    }

    fn tangent_is_oriented(&self) -> bool {
        // The cross product of the two normals, whose sign is the surfaces'
        // own and whose flip at a tangency is what stops the march.
        true
    }

    fn tangent(&self, x: &[f64], tol: Tolerances) -> Option<Vector> {
        tangent_at(
            self.a,
            self.b,
            Contact {
                on_a: (x[0], x[1]),
                on_b: (x[2], x[3]),
                point: Point::ORIGIN,
            },
            tol,
        )
    }
}

/// Walk one way from a seed.
fn walk(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    from: Contact,
    sense: f64,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Traced> {
    let pair = SurfacePair { a, b };
    let start = [from.on_a.0, from.on_a.1, from.on_b.0, from.on_b.1];
    let walked = crate::walk::walk_one_way(&pair, &start, sense, options, tol)?;
    Ok(Traced {
        on_a: walked.states.iter().map(|x| (x[0], x[1])).collect(),
        on_b: walked.states.iter().map(|x| (x[2], x[3])).collect(),
        points: walked.points,
        stopped: walked.stopped,
    })
}

/// The sine of the shallowest crossing angle the marcher will follow.
///
/// One microradian, and the number is set by the *correction*, not by taste.
/// `correct` accepts a residual up to the confusion tolerance, so the two
/// parameter points of a contact can disagree by that much in space — and on
/// coincident or near-coincident surfaces, that disagreement shows up as a
/// spurious angle between the two computed normals of about the residual over
/// the local feature size. A gate below that floor reads the correction's own
/// noise as a direction and marches along it: identical spheres came back as
/// six confident little curves that existed nowhere but in rounding.
///
/// So below this angle the marcher cannot tell an ultra-shallow crossing from
/// coincidence, and refuses both rather than guessing. A genuine crossing
/// shallower than a microradian is also one the Newton correction cannot
/// reliably follow — its travel constraint becomes numerically dependent on
/// the surface-gap rows at exactly the same rate — so the gate refuses what
/// could not have been followed anyway.
const SHALLOWEST: f64 = 1e-6;

/// The direction the intersection runs at a contact.
///
/// The cross product of the two normals: the one direction lying in both
/// tangent planes. `None` where the normals are parallel to within
/// [`SHALLOWEST`] — the surfaces are tangent or coincident there, and the
/// intersection has no direction the marcher can trust.
fn tangent_at(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    at: Contact,
    tol: Tolerances,
) -> Option<Vector> {
    let na = normal_at(a, at.on_a, tol)?;
    let nb = normal_at(b, at.on_b, tol)?;
    let cross = na.cross(nb);
    let length = cross.magnitude();
    // The decision is made through intervals rather than a bare compare:
    // each normal component carries the correction's stated residual as an
    // uncertainty, the squared cross magnitude is computed as an enclosure,
    // and only a crossing *certainly* above the floor is followed. A sine
    // inside the enclosure's undecided band is exactly the case the floor
    // exists for — the correction's own noise masquerading as an angle —
    // and it is refused with a certificate instead of a guess.
    let floor = tol.angular().max(SHALLOWEST);
    let widen = |value: f64| ogeom_math::Interval::about(value, tol.confusion());
    let (ax, ay, az) = (widen(na.x), widen(na.y), widen(na.z));
    let (bx, by, bz) = (widen(nb.x), widen(nb.y), widen(nb.z));
    let cx = ay.mul(&bz).sub(&az.mul(&by));
    let cy = az.mul(&bx).sub(&ax.mul(&bz));
    let cz = ax.mul(&by).sub(&ay.mul(&bx));
    let magnitude2 = cx.square().add(&cy.square()).add(&cz.square());
    let above = magnitude2.sub(&ogeom_math::Interval::point(floor * floor));
    if above.certain_sign() != Some(ogeom_core::Sign::Positive) || length <= f64::MIN_POSITIVE {
        return None;
    }
    Some(cross * (1.0 / length))
}

/// A surface's unit normal at a parameter.
fn normal_at(surface: &SurfaceGeometry, at: (f64, f64), tol: Tolerances) -> Option<Vector> {
    let (du, dv) = surface.d1_at(at.0, at.1, tol).ok()?;
    let cross = du.cross(dv);
    let length = cross.magnitude();
    if length <= tol.confusion() {
        return None;
    }
    Some(cross * (1.0 / length))
}

/// Bring a parameter guess onto both surfaces.
///
/// Three equations say the two surface points coincide; the fourth says how far
/// along the direction of travel to land. Without that fourth the system is
/// underdetermined — its solution set *is* the curve — and Newton would wander
/// along it instead of converging to a point.
///
/// `constraint` is `(anchor, direction, distance)`. Without one, the guess
/// itself is used as the anchor and the direction is the intersection tangent,
/// which is what seeding wants: land anywhere on the curve near here.
fn correct(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    start: [f64; 4],
    guess: Point,
    constraint: Option<(Point, Vector, f64)>,
    tol: Tolerances,
) -> Option<Contact> {
    let (anchor, along, reach) = match constraint {
        Some(given) => given,
        None => {
            // No direction to travel: hold the guess still along whichever way
            // the curve runs, so the solve slides onto the curve rather than
            // along it.
            let at = Contact {
                on_a: (start[0], start[1]),
                on_b: (start[2], start[3]),
                point: guess,
            };
            (guess, tangent_at(a, b, at, tol).unwrap_or(Vector::X), 0.0)
        }
    };

    let system = |x: &[f64]| {
        let (ua, va) = clamp(a, x[0], x[1]);
        let (ub, vb) = clamp(b, x[2], x[3]);
        let pa = a.point_at(ua, va, tol).unwrap_or(Point::ORIGIN);
        let pb = b.point_at(ub, vb, tol).unwrap_or(Point::ORIGIN);
        let (au, av) = a.d1_at(ua, va, tol).unwrap_or((Vector::ZERO, Vector::ZERO));
        let (bu, bv) = b.d1_at(ub, vb, tol).unwrap_or((Vector::ZERO, Vector::ZERO));

        let gap = pa - pb;
        let residual = vec![gap.x, gap.y, gap.z, (pa - anchor).dot(along) - reach];
        let jacobian = vec![
            vec![au.x, av.x, -bu.x, -bv.x],
            vec![au.y, av.y, -bu.y, -bv.y],
            vec![au.z, av.z, -bu.z, -bv.z],
            vec![au.dot(along), av.dot(along), 0.0, 0.0],
        ];
        (residual, jacobian)
    };

    let criteria = solve::Criteria {
        residual: tol.confusion() * 0.01,
        step: tol.parametric(),
        max_iterations: 40,
    };
    let found = solve::newton_system(system, &start, criteria).ok()?;
    if found.residual > tol.confusion() {
        return None;
    }
    let (ua, va) = clamp(a, found.value[0], found.value[1]);
    let (ub, vb) = clamp(b, found.value[2], found.value[3]);
    Some(Contact {
        on_a: (ua, va),
        on_b: (ub, vb),
        point: a.point_at(ua, va, tol).ok()?,
    })
}

/// Hold a parameter inside a surface's domain.
///
/// A periodic direction wraps instead, so a curve crossing a cylinder's seam
/// keeps going rather than stopping at a boundary that is not one.
fn clamp(surface: &SurfaceGeometry, u: f64, v: f64) -> (f64, f64) {
    let ((ua, ub), (va, vb)) = surface.domain();
    let fold = |x: f64, lo: f64, hi: f64, periodic: bool| {
        if !periodic {
            return x.clamp(lo, hi);
        }
        let span = hi - lo;
        if span <= 0.0 {
            return x;
        }
        lo + (x - lo).rem_euclid(span)
    };
    (
        fold(u, ua, ub, surface.is_periodic_u()),
        fold(v, va, vb, surface.is_periodic_v()),
    )
}

/// Whether a parameter sits close enough to a non-periodic edge that a stalled
/// walk there means the edge rather than a singularity.
///
/// The band is a fraction of the domain's own span — a walk stalls within a
/// step of the boundary, and the step is far larger than the strict band
/// [`outside`] uses to decide a point has actually crossed.
fn near_edge(surface: &SurfaceGeometry, at: (f64, f64)) -> bool {
    let ((ua, ub), (va, vb)) = surface.domain();
    let close = |x: f64, lo: f64, hi: f64, periodic: bool| {
        !periodic && {
            let band = (hi - lo).abs() * 1e-4;
            x <= lo + band || x >= hi - band
        }
    };
    close(at.0, ua, ub, surface.is_periodic_u()) || close(at.1, va, vb, surface.is_periodic_v())
}

/// Whether a parameter has left a surface's domain, in a direction that has one.
fn outside(surface: &SurfaceGeometry, at: (f64, f64), tol: Tolerances) -> bool {
    let ((ua, ub), (va, vb)) = surface.domain();
    let past = |x: f64, lo: f64, hi: f64, periodic: bool| {
        !periodic && (x <= lo + tol.parametric() || x >= hi - tol.parametric())
    };
    past(at.0, ua, ub, surface.is_periodic_u()) || past(at.1, va, vb, surface.is_periodic_v())
}

/// A surface's rough size, for spacing seeds.
fn span(surface: &SurfaceGeometry) -> f64 {
    let ((ua, ub), (va, vb)) = surface.domain();
    let tol = Tolerances::millimetres();
    let corners = [(ua, va), (ub, va), (ua, vb), (ub, vb)];
    let mut low = Point::new(f64::MAX, f64::MAX, f64::MAX);
    let mut high = Point::new(f64::MIN, f64::MIN, f64::MIN);
    for (u, v) in corners {
        if let Ok(p) = surface.point_at(u, v, tol) {
            low = Point::new(low.x.min(p.x), low.y.min(p.y), low.z.min(p.z));
            high = Point::new(high.x.max(p.x), high.y.max(p.y), high.z.max(p.z));
        }
    }
    let size = (high - low).magnitude();
    if size.is_finite() && size > 0.0 {
        size
    } else {
        1.0
    }
}

/// One sampled triangle of a surface, with the parameters it came from.
pub(crate) struct Cell {
    pub(crate) corners: [Point; 3],
    pub(crate) at: (f64, f64),
    pub(crate) low: Point,
    pub(crate) high: Point,
}

/// Sample a surface into triangles.
pub(crate) fn sample(surface: &SurfaceGeometry, grid: usize, tol: Tolerances) -> Vec<Cell> {
    let ((ua, ub), (va, vb)) = surface.domain();
    // An unbounded domain would put the samples a billion units apart and find
    // nothing. Clamped to something a real model lives inside.
    let limit = 1.0e6;
    let (ua, ub) = (ua.max(-limit), ub.min(limit));
    let (va, vb) = (va.max(-limit), vb.min(limit));

    let mut out = Vec::new();
    #[allow(clippy::cast_precision_loss)]
    let n = grid as f64;
    for i in 0..grid {
        for j in 0..grid {
            #[allow(clippy::cast_precision_loss)]
            let (s0, s1) = (i as f64 / n, (i + 1) as f64 / n);
            #[allow(clippy::cast_precision_loss)]
            let (t0, t1) = (j as f64 / n, (j + 1) as f64 / n);
            let at = |s: f64, t: f64| {
                let (u, v) = (ua + (ub - ua) * s, va + (vb - va) * t);
                surface.point_at(u, v, tol).map(|p| ((u, v), p))
            };
            let (Ok((p00, a00)), Ok((_, a10)), Ok((_, a01)), Ok((_, a11))) =
                (at(s0, t0), at(s1, t0), at(s0, t1), at(s1, t1))
            else {
                continue;
            };
            for corners in [[a00, a10, a11], [a00, a11, a01]] {
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
                });
            }
        }
    }
    out
}

/// Whether two cells' boxes come within a margin of each other.
fn overlap(a: &Cell, b: &Cell, margin: f64) -> bool {
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
/// the triangle itself. Only an approximate answer is needed — it is a seed, and
/// the Newton correction that follows is what makes it a point on the curve.
fn triangles_cross(a: &Cell, b: &Cell) -> Option<Point> {
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

/// One arc between branch points, cut from a stalled fragment.
struct Arc {
    points: Vec<Point>,
    on_a: Vec<(f64, f64)>,
    on_b: Vec<(f64, f64)>,
    /// The branch-point cluster each end attaches to, if any.
    head_bp: Option<usize>,
    tail_bp: Option<usize>,
}

impl Arc {
    fn length(&self) -> f64 {
        self.points
            .windows(2)
            .map(|pair| pair[0].distance(pair[1]))
            .sum()
    }

    /// The direction the arc leaves its end, read across a window deep
    /// enough to stand clear of any residual wander.
    fn outgoing(&self, tail: bool) -> Option<Vector> {
        let n = self.points.len();
        if n < 2 {
            return None;
        }
        let window = (n - 1).min(24);
        let (at, back) = if tail {
            (n - 1, n - 1 - window)
        } else {
            (0, window)
        };
        let out = self.points[at] - self.points[back];
        let m = out.magnitude();
        (m > f64::MIN_POSITIVE).then(|| out / m)
    }
}

/// Whether a stalled branch is transversal *somewhere* in its interior —
/// the certificate that it is a curve passing branch points rather than
/// tangential-contact debris. A plane resting on a torus produces
/// fragments tangent along their whole length; a real curve through a
/// branch point is tangent only in passing.
fn interior_is_transversal(
    branch: &Traced,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> bool {
    let n = branch.points.len();
    if n < 5 {
        return false;
    }
    [n / 4, n / 2, 3 * n / 4]
        .into_iter()
        .any(|i| crossing_sine(a, b, branch.on_a[i], branch.on_b[i], tol) > BRANCH_POINT_SINE)
}

/// Join stalled fragments that meet at branch points into the curves they
/// belong to — the stitching an earlier plan owed, now delivered.
///
/// Where the two normals become parallel the intersection has no single
/// direction: the walk stalls there, wanders in place while the correction
/// gives out, and a curve that passes *through* the singularity comes back
/// as fragments with parked ends. The reassembly is geometric, not
/// bookkeeping: cluster the near-tangent stall ends into branch points;
/// cut every fragment at its branch-point visits, which trims the wander
/// and separates the arcs a walk-through glued together; drop debris and
/// duplicate coverage; then, at each branch point, pair arc ends whose
/// tangents continue one another — a smooth curve crosses the singularity
/// collinearly, and the crossing curve turns through the crossing angle —
/// and chain the pairs into whole curves, closing the loops that close.
///
/// Only fragments transversal somewhere in their interior participate:
/// tangential contact along a whole curve is a different phenomenon and
/// keeps its honest fragments.
fn stitch_stalled(
    found: Vec<Traced>,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> Vec<Traced> {
    let reach = options.chord.max(tol.confusion()) * 60.0;
    const CONTINUES: f64 = 0.5;

    let (candidates, mut out): (Vec<Traced>, Vec<Traced>) = found.into_iter().partition(|branch| {
        branch.stopped == Stopped::Stalled && interior_is_transversal(branch, a, b, tol)
    });
    if candidates.is_empty() {
        return out;
    }

    // Branch points: the near-tangent stall ends, clustered.
    let mut bps: Vec<Point> = Vec::new();
    for branch in &candidates {
        let n = branch.points.len();
        for at in [0, n - 1] {
            if crossing_sine(a, b, branch.on_a[at], branch.on_b[at], tol) < BRANCH_POINT_SINE {
                let p = branch.points[at];
                if !bps.iter().any(|held| held.distance(p) <= reach) {
                    bps.push(p);
                }
            }
        }
    }
    if bps.is_empty() {
        out.extend(candidates);
        return out;
    }
    let bp_of =
        |p: Point| -> Option<usize> { bps.iter().position(|held| held.distance(p) <= reach) };

    // Cut each fragment at its branch-point visits: maximal runs of points
    // clear of every branch point become arcs, attached to the branch
    // points beside them.
    let mut arcs: Vec<Arc> = Vec::new();
    for branch in &candidates {
        let n = branch.points.len();
        let mut run_start: Option<usize> = None;
        for i in 0..=n {
            let near = i < n && bp_of(branch.points[i]).is_some();
            match (run_start, near, i == n) {
                (None, false, false) => run_start = Some(i),
                (Some(s), true, _) | (Some(s), _, true) => {
                    let e = i;
                    if e > s + 1 {
                        let head_bp = if s > 0 {
                            bp_of(branch.points[s - 1])
                        } else {
                            None
                        };
                        let tail_bp = if e < n { bp_of(branch.points[e]) } else { None };
                        arcs.push(Arc {
                            points: branch.points[s..e].to_vec(),
                            on_a: branch.on_a[s..e].to_vec(),
                            on_b: branch.on_b[s..e].to_vec(),
                            head_bp,
                            tail_bp,
                        });
                    }
                    run_start = None;
                }
                _ => {}
            }
        }
    }

    // Debris and duplicate coverage out; longest first so the fuller
    // tracing of a doubly-walked arc is the one kept.
    arcs.retain(|arc| arc.length() > options.chord * 10.0 && arc.points.len() >= 4);
    arcs.sort_by(|x, y| {
        y.length()
            .partial_cmp(&x.length())
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    let mut kept: Vec<Arc> = Vec::new();
    'candidate: for arc in arcs {
        let n = arc.points.len();
        for probe in [n / 4, n / 2, 3 * n / 4] {
            let p = arc.points[probe];
            if kept.iter().any(|held| {
                held.points
                    .windows(2)
                    .any(|pair| distance_to_segment(p, pair[0], pair[1]) <= reach)
            }) {
                continue 'candidate;
            }
        }
        kept.push(arc);
    }

    // Pair arc ends at each branch point by tangent continuation: the
    // smooth curve runs straight through, the crossing one turns.
    let ends: Vec<(usize, bool, usize, Vector)> = kept
        .iter()
        .enumerate()
        .flat_map(|(i, arc)| {
            [(false, arc.head_bp), (true, arc.tail_bp)]
                .into_iter()
                .filter_map(move |(tail, bp)| Some((i, tail, bp?, arc.outgoing(tail)?)))
        })
        .collect();
    let mut partner: Vec<Option<usize>> = vec![None; ends.len()];
    for bp in 0..bps.len() {
        loop {
            let mut best: Option<(usize, usize, f64)> = None;
            for x in 0..ends.len() {
                if partner[x].is_some() || ends[x].2 != bp {
                    continue;
                }
                for y in (x + 1)..ends.len() {
                    if partner[y].is_some() || ends[y].2 != bp {
                        continue;
                    }
                    let score = -ends[x].3.dot(ends[y].3);
                    if score > CONTINUES && best.is_none_or(|(_, _, held)| score > held) {
                        best = Some((x, y, score));
                    }
                }
            }
            let Some((x, y, _)) = best else { break };
            partner[x] = Some(y);
            partner[y] = Some(x);
        }
    }

    // Chain the arcs through the pairings into whole curves.
    let end_index = |arc: usize, tail: bool| -> Option<usize> {
        ends.iter().position(|e| e.0 == arc && e.1 == tail)
    };
    let mut used = vec![false; kept.len()];
    for start in 0..kept.len() {
        if used[start] {
            continue;
        }
        // Walk backwards first to a free entry, unless the chain loops.
        let mut first = start;
        let mut first_reversed = false;
        let mut seen_back = vec![false; kept.len()];
        loop {
            seen_back[first] = true;
            // Forward traversal enters an arc at its head, reversed at its
            // tail — so the entry end is named by the orientation flag.
            let Some(entry) = end_index(first, first_reversed) else {
                break;
            };
            let Some(p) = partner[entry] else { break };
            let (prev, prev_tail, _, _) = ends[p];
            if seen_back[prev] {
                break; // the chain is a loop; any start serves
            }
            first = prev;
            // The previous arc *leaves* through the paired end: leaving at
            // its tail means it runs forward.
            first_reversed = !prev_tail;
        }

        // Now walk forwards from `first`, consuming arcs.
        let mut points: Vec<Point> = Vec::new();
        let mut on_a: Vec<(f64, f64)> = Vec::new();
        let mut on_b: Vec<(f64, f64)> = Vec::new();
        let mut current = first;
        let mut reversed = first_reversed;
        let mut closed = false;
        loop {
            used[current] = true;
            let arc = &kept[current];
            type Run = (Vec<Point>, Vec<(f64, f64)>, Vec<(f64, f64)>);
            let (pts, pa, pb): Run = if reversed {
                (
                    arc.points.iter().rev().copied().collect(),
                    arc.on_a.iter().rev().copied().collect(),
                    arc.on_b.iter().rev().copied().collect(),
                )
            } else {
                (arc.points.clone(), arc.on_a.clone(), arc.on_b.clone())
            };
            // Insert the branch point itself at the junction.
            if !points.is_empty() {
                let joint_bp = if reversed { arc.tail_bp } else { arc.head_bp };
                if let Some(bp) = joint_bp {
                    points.push(bps[bp]);
                    on_a.push(pa[0]);
                    on_b.push(pb[0]);
                }
            }
            points.extend(pts);
            on_a.extend(pa);
            on_b.extend(pb);

            let leaving = end_index(current, !reversed);
            let Some(l) = leaving else { break };
            let Some(p) = partner[l] else { break };
            let (next, next_tail, _, _) = ends[p];
            if used[next] {
                closed = next == first;
                break;
            }
            current = next;
            reversed = next_tail;
        }
        if closed && points.len() > 3 {
            let bridge = points[0];
            let ba = on_a[0];
            let bb = on_b[0];
            points.push(bridge);
            on_a.push(ba);
            on_b.push(bb);
        }
        out.push(Traced {
            points,
            on_a,
            on_b,
            stopped: if closed {
                Stopped::Closed
            } else {
                Stopped::Stalled
            },
        });
    }
    out
}

/// Newton projection of a point onto a surface, warm-started — the local
/// tool the tangential walker corrects with.
fn nearest_on(
    surface: &SurfaceGeometry,
    seed: (f64, f64),
    target: Point,
    tol: Tolerances,
) -> Option<((f64, f64), Point)> {
    let (mut u, mut v) = seed;
    for _ in 0..16 {
        let (u_ok, v_ok) = surface.normalize_parameters(u, v, tol).ok()?;
        u = u_ok;
        v = v_ok;
        let p = surface.point_at(u, v, tol).ok()?;
        let (su, sv) = surface.d1_at(u, v, tol).ok()?;
        let r = p - target;
        let (a11, a12, a22) = (su.dot(su), su.dot(sv), sv.dot(sv));
        let det = a11.mul_add(a22, -(a12 * a12));
        if det.abs() <= f64::MIN_POSITIVE {
            break;
        }
        let (b1, b2) = (-su.dot(r), -sv.dot(r));
        let du = b1.mul_add(a22, -(b2 * a12)) / det;
        let dv = a11.mul_add(b2, -(a12 * b1)) / det;
        u += du;
        v += dv;
        if du.hypot(dv) < 1e-14 {
            break;
        }
    }
    let (u, v) = surface.normalize_parameters(u, v, tol).ok()?;
    Some(((u, v), surface.point_at(u, v, tol).ok()?))
}

/// Trace tangential contact along a curve — the walker an earlier plan
/// owed, following the valley of the gap function rather than a crossing.
///
/// Where two surfaces touch along a whole curve there is no transversal
/// direction to march: the crossing angle is zero along the entire
/// contact, and the crossing walker honestly stalls. But the contact is
/// still a curve, and it is the locus where the *gap* between the surfaces
/// stays zero. This walker steps along the contact and corrects each step
/// transversally: project the candidate onto the first surface, project
/// that onto the second, and slide on the first surface to close the gap —
/// a minimization, not a root-find, because at tangency the gap touches
/// zero without crossing it.
///
/// The seed must be a genuine contact: on both surfaces within tolerance
/// and near-tangent there. A transversal crossing is refused — the
/// ordinary walker owns those.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// settings are unusable or the seed is not a tangential contact.
pub fn trace_tangential(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    from: Contact,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Traced> {
    options.validate()?;
    let accept = tol.confusion() * 100.0;
    let sine = crossing_sine(a, b, from.on_a, from.on_b, tol);
    if sine > BRANCH_POINT_SINE {
        ogeom_bail!(
            Construction,
            "the surfaces cross here at sine {sine}; tangential tracing wants a contact"
        );
    }
    let reach = span(a).max(span(b));
    let step = (options.chord * reach)
        .sqrt()
        .clamp(tol.confusion(), reach / 16.0);

    type Walked = (Vec<Point>, Vec<(f64, f64)>, Vec<(f64, f64)>, Stopped);
    let walk_one = |sense: f64| -> OgeomResult<Walked> {
        let mut points = vec![from.point];
        let mut on_a = vec![from.on_a];
        let mut on_b = vec![from.on_b];
        let mut at = from;
        let mut previous: Option<Vector> = None;
        let mut stopped = Stopped::RanOut;
        while points.len() < options.max_points {
            ogeom_core::progress::checkpoint()?;
            // The contact direction: in the common tangent plane. With the
            // normals parallel, one surface's normal serves for both; the
            // step direction is the previous one projected back into the
            // tangent plane, or any tangent direction to begin with.
            let Some(normal) = normal_at(a, at.on_a, tol) else {
                stopped = Stopped::Stalled;
                break;
            };
            let direction = match previous {
                Some(d) => {
                    let flat = d - normal * d.dot(normal);
                    let m = flat.magnitude();
                    if m <= f64::MIN_POSITIVE {
                        stopped = Stopped::Stalled;
                        break;
                    }
                    flat / m
                }
                None => {
                    // First step: the tangent direction along which the gap
                    // grows least, found by sampling the tangent circle.
                    let (su, _) = a.d1_at(at.on_a.0, at.on_a.1, tol).map_err(|_| {
                        ogeom_core::ogeom_err!(Construction, "the seed cannot be evaluated")
                    })?;
                    let t1 = {
                        let flat = su - normal * su.dot(normal);
                        let m = flat.magnitude();
                        if m <= f64::MIN_POSITIVE {
                            stopped = Stopped::Stalled;
                            break;
                        }
                        flat / m
                    };
                    let t2 = normal.cross(t1);
                    let mut best = (f64::INFINITY, t1);
                    for k in 0..16 {
                        let angle = core::f64::consts::TAU * f64::from(k) / 16.0;
                        let dir = t1 * angle.cos() + t2 * angle.sin();
                        let probe = at.point + dir * step;
                        let Some((_, qa)) = nearest_on(a, at.on_a, probe, tol) else {
                            continue;
                        };
                        let Some((_, qb)) = nearest_on(b, at.on_b, qa, tol) else {
                            continue;
                        };
                        let gap = qa.distance(qb);
                        if gap < best.0 {
                            best = (gap, dir);
                        }
                    }
                    best.1 * sense
                }
            };

            // Step and correct: onto a, gap closed against b by sliding on
            // a a few times.
            let mut candidate = at.point + direction * step;
            let mut pa = at.on_a;
            let mut pb = at.on_b;
            let mut gap = f64::INFINITY;
            for _ in 0..8 {
                let Some((ua, qa)) = nearest_on(a, pa, candidate, tol) else {
                    break;
                };
                let Some((ub, qb)) = nearest_on(b, pb, qa, tol) else {
                    break;
                };
                pa = ua;
                pb = ub;
                gap = qa.distance(qb);
                if gap <= tol.confusion() {
                    candidate = qa;
                    break;
                }
                // Slide the working point toward the midpoint of the gap.
                candidate = qa + (qb - qa) * 0.5;
            }
            if gap > accept {
                stopped = Stopped::Stalled;
                break;
            }
            let next = Contact {
                on_a: pa,
                on_b: pb,
                point: candidate,
            };
            if points.len() > 3 && next.point.distance(from.point) <= step {
                points.push(from.point);
                on_a.push(from.on_a);
                on_b.push(from.on_b);
                stopped = Stopped::Closed;
                break;
            }
            if next.point.distance(at.point) <= step * 1e-3 {
                stopped = Stopped::Stalled;
                break;
            }
            previous = Some(next.point - at.point);
            points.push(next.point);
            on_a.push(next.on_a);
            on_b.push(next.on_b);
            at = next;
        }
        Ok((points, on_a, on_b, stopped))
    };

    let (points, on_a, on_b, stopped) = walk_one(1.0)?;
    if stopped == Stopped::Closed {
        return Ok(Traced {
            points,
            on_a,
            on_b,
            stopped,
        });
    }
    let (mut back_points, mut back_a, mut back_b, back_stopped) = walk_one(-1.0)?;
    back_points.reverse();
    back_a.reverse();
    back_b.reverse();
    back_points.pop();
    back_a.pop();
    back_b.pop();
    back_points.extend(points);
    back_a.extend(on_a);
    back_b.extend(on_b);
    let stopped = if stopped == Stopped::RanOut || back_stopped == Stopped::RanOut {
        Stopped::RanOut
    } else if stopped == Stopped::Stalled || back_stopped == Stopped::Stalled {
        Stopped::Stalled
    } else {
        Stopped::LeftTheDomain
    };
    Ok(Traced {
        points: back_points,
        on_a: back_a,
        on_b: back_b,
        stopped,
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::print_stdout)]
mod tests {
    use super::*;
    use ogeom_geom::{CylinderSurface, PlaneSurface, SphereSurface};
    use ogeom_math::{Cylinder, Direction, Frame, Plane, Sphere};

    const T: Tolerances = Tolerances::millimetres();

    fn cylinder(origin: Point, axis: Vector, radius: f64, height: (f64, f64)) -> SurfaceGeometry {
        let frame = Frame::new(
            origin,
            Direction::new(axis, T).unwrap(),
            Direction::from_cross(axis, Vector::new(0.3, 0.5, 0.9), T).unwrap(),
            T,
        )
        .unwrap();
        CylinderSurface::new(Cylinder::new(frame, radius, T).unwrap(), height)
            .unwrap()
            .into()
    }

    fn sphere(centre: Point, radius: f64) -> SurfaceGeometry {
        SphereSurface::new(Sphere::centred(centre, radius, T).unwrap()).into()
    }

    fn plane(origin: Point, normal: Vector) -> SurfaceGeometry {
        PlaneSurface::over(
            Plane::through(origin, Direction::new(normal, T).unwrap()),
            (-8.0, 8.0),
            (-8.0, 8.0),
        )
        .unwrap()
        .into()
    }

    /// How far a point is from a quadric, in closed form.
    fn off(surface: &SurfaceGeometry, p: Point) -> f64 {
        match surface {
            SurfaceGeometry::Plane(x) => x.plane().distance_to(p),
            SurfaceGeometry::Sphere(x) => x.sphere().distance_to(p),
            SurfaceGeometry::Cylinder(x) => x.cylinder().distance_to(p),
            _ => 0.0,
        }
    }

    /// The worst distance from a traced branch to either surface.
    fn deviation(a: &SurfaceGeometry, b: &SurfaceGeometry, traced: &Traced) -> f64 {
        traced
            .points
            .iter()
            .map(|p| off(a, *p).abs().max(off(b, *p).abs()))
            .fold(0.0_f64, f64::max)
    }

    #[test]
    fn two_crossed_cylinders_are_traced_onto_both_of_them() {
        // The case with no closed form, and the one the analytic module
        // explicitly refuses: two cylinders on perpendicular axes meet in a
        // quartic space curve. Every point of the trace must be on both.
        let a = cylinder(Point::ORIGIN, Vector::Z, 1.0, (-4.0, 4.0));
        let b = cylinder(Point::ORIGIN, Vector::X, 1.0, (-4.0, 4.0));
        let options = Marching {
            chord: 1e-5,
            ..Marching::default()
        };

        let found = branches(&a, &b, options, T).unwrap();
        assert_eq!(
            found.len(),
            2,
            "two equal cylinders crossing at right angles meet in two closed \
             curves — the Steinmetz solid's seams"
        );

        let mut worst = 0.0_f64;
        for branch in &found {
            assert!(branch.closed(), "each seam is a closed loop");
            assert!(
                branch.points.len() > 100,
                "a branch of only {} points",
                branch.points.len()
            );
            worst = worst.max(deviation(&a, &b, branch));
        }
        println!(
            "crossed cylinders: {} branches, worst deviation {worst:e}",
            found.len()
        );
        assert!(worst < 1e-7, "traced off the surfaces by {worst:e}");
    }

    #[test]
    fn unequal_crossed_cylinders_meet_in_two_curves_as_well() {
        // The non-degenerate cousin. Equal radii put the two curves through
        // each other at the tangency points, so getting the count right there
        // says less than getting it right here, where there is no singularity
        // for a tracer to be lucky about.
        let a = cylinder(Point::ORIGIN, Vector::Z, 1.0, (-4.0, 4.0));
        let b = cylinder(Point::ORIGIN, Vector::X, 1.6, (-4.0, 4.0));
        let options = Marching {
            chord: 1e-5,
            ..Marching::default()
        };

        let found = branches(&a, &b, options, T).unwrap();
        assert_eq!(found.len(), 2);
        for branch in &found {
            assert!(branch.closed());
            assert!(deviation(&a, &b, branch) < 1e-7);
        }
    }

    #[test]
    fn a_traced_circle_agrees_with_the_circle_it_should_be() {
        // A sphere cut by a plane through its centre is a circle of known
        // radius, and the marcher does not know that. Tracing it and checking
        // against the closed form is the strongest single check there is: it
        // tests the tracer against an answer derived independently of it.
        let s = sphere(Point::ORIGIN, 3.0);
        let cut = plane(Point::ORIGIN, Vector::Z);
        let options = Marching {
            chord: 1e-6,
            ..Marching::default()
        };

        let found = seeds(&s, &cut, options, T).unwrap();
        assert!(!found.is_empty());
        let branch = trace(&s, &cut, found[0], options, T).unwrap();

        assert!(
            branch.closed(),
            "a plane through a sphere gives a closed loop"
        );
        for p in &branch.points {
            let radius = (p.x * p.x + p.y * p.y).sqrt();
            assert!(
                (radius - 3.0).abs() < 1e-7,
                "a point at radius {radius} on a circle of 3"
            );
            assert!(p.z.abs() < 1e-7, "off the cutting plane by {}", p.z);
        }
    }

    #[test]
    fn a_branch_that_leaves_the_surface_says_so_rather_than_stopping_quietly() {
        // A truncated branch and a finished one are different answers, and a
        // caller that cannot tell them apart treats the first as the second.
        let s = sphere(Point::ORIGIN, 3.0);
        let cut = plane(Point::new(0.0, 0.0, 0.0), Vector::Z);
        let options = Marching {
            chord: 1e-4,
            max_points: 8,
            ..Marching::default()
        };
        let found = seeds(&s, &cut, options, T).unwrap();
        let branch = trace(&s, &cut, found[0], options, T).unwrap();
        assert_eq!(branch.stopped, Stopped::RanOut);
        assert!(!branch.complete(), "a truncated branch is not complete");
    }

    #[test]
    fn tangent_surfaces_are_refused_rather_than_followed_onto_a_guess() {
        // Where the normals are parallel the intersection has no single
        // direction, and marching through such a point is how a tracer changes
        // branch without noticing. A sphere resting on a plane is the case.
        let s = sphere(Point::new(0.0, 0.0, 3.0), 3.0);
        let ground = plane(Point::ORIGIN, Vector::Z);
        let touch = Contact {
            on_a: (0.0, -core::f64::consts::FRAC_PI_2),
            on_b: (0.0, 0.0),
            point: Point::ORIGIN,
        };
        let err = trace(&s, &ground, touch, Marching::default(), T).unwrap_err();
        assert!(err.to_string().contains("tangent"), "unexpected: {err}");
    }

    #[test]
    fn the_number_of_branches_is_the_number_there_are() {
        // The failure the accuracy measure cannot see. Every point of one
        // circle is on both surfaces, so returning one of two scores perfectly
        // — the count is the only thing that catches it.
        let options = Marching {
            chord: 1e-5,
            ..Marching::default()
        };

        // A sphere cut by a plane off its centre: one circle.
        let one = branches(
            &sphere(Point::ORIGIN, 3.0),
            &plane(Point::new(0.0, 0.0, 1.0), Vector::Z),
            options,
            T,
        )
        .unwrap();
        assert_eq!(one.len(), 1, "one plane through a sphere cuts one circle");
        assert!(one[0].closed());

        // A sphere and a coaxial cylinder narrower than it: two circles, one
        // above and one below.
        let two = branches(
            &sphere(Point::ORIGIN, 3.0),
            &cylinder(Point::ORIGIN, Vector::Z, 1.5, (-4.0, 4.0)),
            options,
            T,
        )
        .unwrap();
        assert_eq!(two.len(), 2, "a coaxial cylinder cuts a sphere twice");
        for branch in &two {
            assert!(branch.closed(), "each is a closed circle");
        }
        // And they are on opposite sides, rather than the same one twice.
        let heights: Vec<f64> = two.iter().map(|b| b.points[0].z).collect();
        assert!(
            heights[0] * heights[1] < 0.0,
            "both branches came back on the same side: {heights:?}"
        );
    }

    #[test]
    fn a_branch_thinner_than_the_sampling_is_missed_and_the_knob_finds_it() {
        // The stated limitation of polyhedral seeding, pinned so it is a known
        // boundary rather than a surprise. Two spheres barely overlapping meet
        // in a small circle; a coarse grid steps over it entirely.
        let a = sphere(Point::ORIGIN, 3.0);
        let b = sphere(Point::new(5.98, 0.0, 0.0), 3.0);

        let coarse = seeds(
            &a,
            &b,
            Marching {
                grid: 6,
                ..Marching::default()
            },
            T,
        )
        .unwrap();
        let fine = seeds(
            &a,
            &b,
            Marching {
                grid: 120,
                ..Marching::default()
            },
            T,
        )
        .unwrap();
        assert!(
            coarse.len() < fine.len(),
            "a finer grid should find what a coarse one steps over: {} against {}",
            coarse.len(),
            fine.len()
        );
        assert!(!fine.is_empty(), "the branch is there to be found");
    }

    #[test]
    fn settings_that_could_not_work_are_refused() {
        let a = sphere(Point::ORIGIN, 1.0);
        let b = plane(Point::ORIGIN, Vector::Z);
        for options in [
            Marching {
                chord: 0.0,
                ..Marching::default()
            },
            Marching {
                grid: 1,
                ..Marching::default()
            },
            Marching {
                max_points: 1,
                ..Marching::default()
            },
        ] {
            assert!(seeds(&a, &b, options, T).is_err());
        }
    }
}
