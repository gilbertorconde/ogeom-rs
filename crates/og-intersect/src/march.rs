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

use og_core::{OgResult, Tolerances, og_bail};
use og_geom::{Surface, SurfaceGeometry};
use og_math::{Point, Vector, solve};

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
    /// [`OgError::Construction`](og_core::OgError::Construction) if the chord is
    /// not positive, the grid is too coarse to hold a triangle, or no points are
    /// allowed.
    pub fn validate(&self) -> OgResult<()> {
        if !self.chord.is_finite() || self.chord <= 0.0 {
            og_bail!(Construction, "a chord of {} is not a distance", self.chord);
        }
        if self.grid < 2 {
            og_bail!(Construction, "a sampling grid needs at least two steps");
        }
        if self.max_points < 2 {
            og_bail!(Construction, "a branch needs at least two points");
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
/// [`OgError::Construction`](og_core::OgError::Construction) if the settings are
/// unusable.
pub fn seeds(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> OgResult<Vec<Contact>> {
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
/// [`OgError::Construction`](og_core::OgError::Construction) if the settings are
/// unusable.
pub fn branches(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> OgResult<Vec<Traced>> {
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
        {
            out.push(branch);
        }
    }
    Ok(out)
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
/// [`OgError::Construction`](og_core::OgError::Construction) if the settings are
/// unusable; [`OgError::NotDone`](og_core::OgError::NotDone) if the two surfaces
/// are tangent at the seed, where there is no single direction to follow.
pub fn trace(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    from: Contact,
    options: Marching,
    tol: Tolerances,
) -> OgResult<Traced> {
    options.validate()?;
    if tangent_at(a, b, from, tol).is_none() {
        og_bail!(
            NotDone,
            "the surfaces are tangent here, so the intersection has no single \
             direction to follow; that is a branch point and needs the seed \
             moved off it"
        );
    }

    // Forwards first. If it closes, that is the whole branch and there is
    // nothing behind us.
    let ahead = walk(a, b, from, 1.0, options, tol);
    if ahead.stopped == Stopped::Closed {
        return Ok(ahead);
    }
    let behind = walk(a, b, from, -1.0, options, tol);

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

/// Walk one way from a seed.
fn walk(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    from: Contact,
    sense: f64,
    options: Marching,
    tol: Tolerances,
) -> Traced {
    let mut points = vec![from.point];
    let mut on_a = vec![from.on_a];
    let mut on_b = vec![from.on_b];
    let mut at = from;
    let mut stopped = Stopped::RanOut;

    // The step is set by how far the chord may sag from the arc, and the sag is
    // measured rather than assumed: it is about `h * turn / 8`, where `turn` is
    // the angle between successive tangents. So the step that just meets the
    // tolerance is found by control rather than by a constant — an earlier
    // version capped it at a multiple of the chord, which made a circle of
    // radius three need three hundred thousand points to get round.
    let reach = span(a).max(span(b));
    let ceiling = reach / 8.0;
    let mut step = (options.chord * reach)
        .sqrt()
        .clamp(tol.confusion(), ceiling);

    while points.len() < options.max_points {
        let Some(direction) = tangent_at(a, b, at, tol) else {
            stopped = Stopped::Stalled;
            break;
        };
        let along = direction * sense;

        let mut taken = None;
        for _ in 0..40 {
            let guess = at.point + along * step;
            let start = [at.on_a.0, at.on_a.1, at.on_b.0, at.on_b.1];
            let Some(next) = correct(a, b, start, guess, Some((at.point, along, step)), tol) else {
                step *= 0.5;
                if step <= tol.confusion() {
                    break;
                }
                continue;
            };
            // How far the curve turned over the step, and therefore how far the
            // straight chord sits from it.
            let turn = tangent_at(a, b, next, tol)
                .map_or(0.0, |t| direction.dot(t).clamp(-1.0, 1.0).acos());
            let sag = step * turn / 8.0;
            if sag <= options.chord || step <= tol.confusion() * 8.0 {
                // Aim the next step at exactly the tolerance. Sag grows with
                // the square of the step, so the correction is a square root,
                // and it is damped so one tight corner does not make the whole
                // rest of the curve expensive or one straight stretch overshoot.
                let scale = if sag > 0.0 {
                    (options.chord / sag).sqrt().clamp(0.5, 2.0)
                } else {
                    2.0
                };
                taken = Some((next, (step * scale).clamp(tol.confusion(), ceiling)));
                break;
            }
            step *= (options.chord / sag).sqrt().clamp(0.25, 0.9);
        }
        let Some((next, following)) = taken else {
            stopped = Stopped::Stalled;
            break;
        };

        // Back where we started: a closed loop. Only checked once the walk has
        // gone far enough to have left, or every branch would close at once.
        if points.len() > 3 && next.point.distance(from.point) <= step {
            points.push(from.point);
            on_a.push(from.on_a);
            on_b.push(from.on_b);
            stopped = Stopped::Closed;
            break;
        }
        if outside(a, next.on_a, tol) || outside(b, next.on_b, tol) {
            stopped = Stopped::LeftTheDomain;
            break;
        }

        points.push(next.point);
        on_a.push(next.on_a);
        on_b.push(next.on_b);
        at = next;
        step = following;
    }

    Traced {
        points,
        on_a,
        on_b,
        stopped,
    }
}

/// The direction the intersection runs at a contact.
///
/// The cross product of the two normals: the one direction lying in both
/// tangent planes. `None` where the normals are parallel — the surfaces are
/// tangent, and the intersection has no single direction there.
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
    // Scaled against the normals, which are unit, so this is the sine of the
    // angle between the surfaces rather than an absolute length.
    if length <= tol.angular().max(1e-9) {
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
struct Cell {
    corners: [Point; 3],
    at: (f64, f64),
    low: Point,
    high: Point,
}

/// Sample a surface into triangles.
fn sample(surface: &SurfaceGeometry, grid: usize, tol: Tolerances) -> Vec<Cell> {
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
fn segment_meets_triangle(from: Point, to: Point, t: [Point; 3]) -> Option<Point> {
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

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::print_stdout)]
mod tests {
    use super::*;
    use og_geom::{CylinderSurface, PlaneSurface, SphereSurface};
    use og_math::{Cylinder, Direction, Frame, Plane, Sphere};

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
