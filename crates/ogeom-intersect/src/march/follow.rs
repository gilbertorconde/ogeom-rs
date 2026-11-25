//! Following one branch from a seed: the walk, its ends, and the Newton correction.

use super::*;

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
    // nothing behind. A backward half that closes is the whole branch too,
    // whatever stopped the forward one, turned to run the way that set off.
    let ahead = walk(a, b, from, 1.0, options, tol)?;
    if ahead.stopped == Stopped::Closed {
        return Ok(ahead);
    }
    let behind = walk(a, b, from, -1.0, options, tol)?;
    if behind.stopped == Stopped::Closed {
        return Ok(reversed(behind));
    }
    let last_step = |walked: &[Point]| -> f64 {
        walked
            .windows(2)
            .last()
            .map_or(0.0, |w| w[0].distance(w[1]))
    };
    let steps = last_step(&ahead.points).max(last_step(&behind.points));

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
    let mut stopped = if ahead.stopped == Stopped::RanOut || behind.stopped == Stopped::RanOut {
        Stopped::RanOut
    } else if ahead.stopped == Stopped::Stalled || behind.stopped == Stopped::Stalled {
        Stopped::Stalled
    } else {
        Stopped::LeftTheDomain
    };
    // A loop cut at a seam. A closed patch is clamped, not periodic: a
    // section that runs round it (a rim's circle on a converted drum)
    // is walked from the seed to the seam one way and to the seam the
    // other, each walk stopping a fraction of a step short of it, and the
    // two ends meet where the surface closes on itself. That is the whole
    // loop, and it is closed: left open, the arrangement downstream would
    // hold a circle with two ends at one point and find no face piece to keep.
    // The ends are within a couple of the walks' own last steps of each
    // other, and the loop is closed exactly on its first point.
    if stopped == Stopped::LeftTheDomain && points.len() > 3 {
        let gap = points[0].distance(points[points.len() - 1]);
        if gap <= (steps * 2.0).max(tol.confusion() * 10.0) {
            // Both walks may have landed on the seam itself, one point.
            if gap <= tol.confusion() {
                points.pop();
                on_a.pop();
                on_b.pop();
            }
            points.push(points[0]);
            on_a.push(on_a[0]);
            on_b.push(on_b[0]);
            stopped = Stopped::Closed;
        }
    }
    Ok(Traced {
        points,
        on_a,
        on_b,
        stopped,
    })
}

/// A traced curve run the other way.
pub(crate) fn reversed(mut traced: Traced) -> Traced {
    traced.points.reverse();
    traced.on_a.reverse();
    traced.on_b.reverse();
    traced
}

/// Two surfaces, as a condition for the walker: four unknowns and three
/// equations saying the two points coincide.
///
/// The intersector's own walk goes through [`crate::walk`] like everything
/// else, and what is *not* generic lives here: the domain clamps a surface
/// pair needs, and the tangent, which the intersector computes from the two
/// normals rather than from the null space so that it can refuse a crossing
/// too shallow to be more than the correction's own noise.
pub(super) struct SurfacePair<'s> {
    pub(super) a: &'s SurfaceGeometry,
    pub(super) b: &'s SurfaceGeometry,
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
        Some(self.system_at(x, tol)?.0)
    }

    fn system_at(
        &self,
        x: &[f64],
        tol: Tolerances,
    ) -> Option<((Vec<f64>, Vec<Vec<f64>>), Point, Vec<Vector>)> {
        let (pa, au, av) = self.a.point_d1_at(x[0], x[1], tol).ok()?;
        let (pb, bu, bv) = self.b.point_d1_at(x[2], x[3], tol).ok()?;
        let gap = pa - pb;
        Some((
            (
                vec![gap.x, gap.y, gap.z],
                vec![
                    vec![au.x, av.x, -bu.x, -bv.x],
                    vec![au.y, av.y, -bu.y, -bv.y],
                    vec![au.z, av.z, -bu.z, -bv.z],
                ],
            ),
            pa,
            vec![au, av, Vector::ZERO, Vector::ZERO],
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

    fn tangent_from(
        &self,
        x: &[f64],
        jacobian: &[Vec<f64>],
        gradient: &[Vector],
        tol: Tolerances,
    ) -> Option<Vector> {
        // The system's columns are the two surfaces' first derivatives, the
        // second's negated: the normals the tangent is made of.
        if jacobian.len() != 3 || jacobian.iter().any(|row| row.len() != 4) {
            return self.tangent(x, tol);
        }
        let _ = gradient;
        let column = |c: usize| Vector::new(jacobian[0][c], jacobian[1][c], jacobian[2][c]);
        let na = unit_normal(column(0), column(1), tol)?;
        let nb = unit_normal(-column(2), -column(3), tol)?;
        tangent_of_normals(na, nb, tol)
    }
}

/// Walk one way from a seed.
pub(super) fn walk(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    from: Contact,
    sense: f64,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Traced> {
    let pair = SurfacePair { a, b };
    let start = [from.on_a.0, from.on_a.1, from.on_b.0, from.on_b.1];
    let (mut walked, past) = crate::walk::walk_one_way_past(&pair, &start, sense, options, tol)?;
    if walked.stopped == Stopped::LeftTheDomain
        && !past.is_some_and(|past| end_at_touch(&pair, &mut walked, &past, tol))
    {
        land_on_edge(&pair, &mut walked, tol);
    }
    Ok(Traced {
        on_a: walked.states.iter().map(|x| (x[0], x[1])).collect(),
        on_b: walked.states.iter().map(|x| (x[2], x[3])).collect(),
        points: walked.points,
        stopped: walked.stopped,
    })
}

/// End a walk that left a surface's domain on the edge it left through.
///
/// The walk stops a fraction of a step short of the edge, and a section
/// ending there misses the edge by that much: across a patch boundary the
/// section on the next patch starts on the edge, and the two ends stand
/// apart by more than either curve's tolerance. The last point is carried
/// onto the nearest non-periodic bound of either surface, the pair still
/// meeting there, and kept if it lies ahead within a couple of steps.
fn land_on_edge(pair: &SurfacePair<'_>, walked: &mut crate::walk::Walked, tol: Tolerances) {
    use crate::walk::Condition as _;
    let n = walked.states.len();
    if n < 2 {
        return;
    }
    let (last, before) = (&walked.states[n - 1], &walked.states[n - 2]);
    let step = walked.points[n - 1].distance(walked.points[n - 2]);
    let heading = walked.points[n - 1] - walked.points[n - 2];
    let mut bounds: Vec<(f64, usize, f64)> = Vec::new();
    for (surface, first) in [(pair.a, 0_usize), (pair.b, 2)] {
        let ((ua, ub), (va, vb)) = surface.domain();
        for (k, (lo, hi), periodic) in [
            (first, (ua, ub), surface.is_periodic_u()),
            (first + 1, (va, vb), surface.is_periodic_v()),
        ] {
            let span = hi - lo;
            if periodic || !span.is_finite() || span <= 0.0 {
                continue;
            }
            // Only a bound the walk was heading for, and one that is an
            // edge: a sphere's pole bounds its chart at a point.
            let moving = last[k] - before[k];
            let across = if k == first {
                surface.domain().1
            } else {
                surface.domain().0
            };
            let collapses = |bound: f64| {
                let at = |t: f64| {
                    if k == first {
                        surface.point_at(bound, t, tol)
                    } else {
                        surface.point_at(t, bound, tol)
                    }
                };
                match (at(across.0), at(f64::midpoint(across.0, across.1))) {
                    (Ok(p), Ok(q)) => p.distance(q) <= tol.confusion(),
                    _ => true,
                }
            };
            for (bound, toward) in [(lo, moving < 0.0), (hi, moving > 0.0)] {
                if (toward || (last[k] - bound).abs() <= span * 1e-4) && !collapses(bound) {
                    bounds.push(((last[k] - bound).abs() / span, k, bound));
                }
            }
        }
    }
    bounds.sort_by(|x, y| x.0.total_cmp(&y.0));
    for (_, k, bound) in bounds {
        let system = |x: &[f64; 4]| {
            let mut residual = [f64::INFINITY; 4];
            let mut jacobian = [[0.0; 4]; 4];
            let Some(((rows, matrix), _, _)) = pair.system_at(x, tol) else {
                return (residual, jacobian);
            };
            residual[..3].copy_from_slice(&rows);
            for (to, row) in jacobian.iter_mut().zip(&matrix) {
                to.copy_from_slice(row);
            }
            residual[3] = x[k] - bound;
            jacobian[3][k] = 1.0;
            (residual, jacobian)
        };
        let criteria = solve::Criteria {
            residual: tol.confusion() * 0.01,
            step: tol.parametric(),
            max_iterations: 40,
        };
        let Ok(start) = <[f64; 4]>::try_from(last.as_slice()) else {
            return;
        };
        let Ok((at, norm, _, _)) = solve::newton_system_fixed(system, start, criteria) else {
            continue;
        };
        if norm > tol.confusion() {
            continue;
        }
        let within = [(pair.a, 0_usize), (pair.b, 2)]
            .iter()
            .all(|(surface, first)| {
                let ((ua, ub), (va, vb)) = surface.domain();
                [
                    (*first, ua, ub, surface.is_periodic_u()),
                    (first + 1, va, vb, surface.is_periodic_v()),
                ]
                .iter()
                .all(|&(i, lo, hi, periodic)| {
                    periodic || (at[i] >= lo - tol.parametric() && at[i] <= hi + tol.parametric())
                })
            });
        if !within {
            continue;
        }
        let Ok(point) = pair.a.point_at(at[0], at[1], tol) else {
            continue;
        };
        let ahead = point - walked.points[n - 1];
        if ahead.magnitude() > (step * 2.0).max(tol.confusion())
            || ahead.dot(heading) < -tol.confusion() * step
        {
            continue;
        }
        if ahead.magnitude() <= tol.confusion() {
            walked.states[n - 1] = at.to_vec();
            walked.points[n - 1] = point;
        } else {
            walked.states.push(at.to_vec());
            walked.points.push(point);
        }
        return;
    }
}

/// End a walk whose last step landed on a surface's edge it only touches,
/// at the touch: whether it did.
///
/// A section touching an edge, a drum's rim lying in a plane tangent to the
/// other surface, runs along the edge to second order: within the
/// confusion distance of the other surface along a stretch of the edge
/// about `√(2·ε·ρ)` long, `ρ` the radius their parting curves at. A step
/// landing there is held on the edge with its residual within the
/// correction's, and reads as having left. Solving the section onto the
/// edge from there is singular, its root a double one: the solve stops
/// where the residual, the square of the distance left, is small, short of
/// the touch by the square root of the residual. The touch is solved as one
/// instead, from the landing, by [`touch_on_edge`], and kept where it lies
/// ahead of the walk's last point.
fn end_at_touch(
    pair: &SurfacePair<'_>,
    walked: &mut crate::walk::Walked,
    (landed, beyond): &crate::walk::Past,
    tol: Tolerances,
) -> bool {
    let n = walked.points.len();
    if n == 0 || landed.len() != 4 {
        return false;
    }
    let heading = *beyond - walked.points[n - 1];
    for (surface, first) in [(pair.a, 0_usize), (pair.b, 2)] {
        let ((ua, ub), (va, vb)) = surface.domain();
        for (k, (lo, hi), periodic) in [
            (first, (ua, ub), surface.is_periodic_u()),
            (first + 1, (va, vb), surface.is_periodic_v()),
        ] {
            if periodic {
                continue;
            }
            for bound in [lo, hi] {
                if (landed[k] - bound).abs() > tol.parametric() {
                    continue;
                }
                let Some(at) = touch_on_edge(pair, landed, k, bound, tol) else {
                    continue;
                };
                let Ok(point) = pair.a.point_at(at[0], at[1], tol) else {
                    continue;
                };
                let ahead = point - walked.points[n - 1];
                if ahead.dot(heading) <= 0.0 {
                    continue;
                }
                if ahead.magnitude() <= tol.confusion() {
                    walked.states[n - 1] = at.to_vec();
                    walked.points[n - 1] = point;
                } else {
                    walked.states.push(at.to_vec());
                    walked.points.push(point);
                }
                return true;
            }
        }
    }
    false
}

/// Where the section touches a surface's edge without crossing it, solved
/// from a state beside it: the state there, or `None` if the edge has no
/// such touch there.
///
/// `k` is the unknown the edge holds at `bound`. Along the edge the
/// bordered surface's point runs as a curve of its free parameter, and the
/// section touches the edge where that curve touches the other surface, as
/// [`crate::curve_surface::touch_solved`] solves it. The touch is kept only
/// within the confusion distance, and only where the section leaves the
/// edge into the domain: a section outside it, touching the edge from
/// there, meets the domain at the touch alone, and a walk along the edge
/// beside it is the correction held on the edge, not the section.
pub(super) fn touch_on_edge(
    pair: &SurfacePair<'_>,
    from: &[f64],
    k: usize,
    bound: f64,
    tol: Tolerances,
) -> Option<[f64; 4]> {
    let (edged, other, first, rest) = if k < 2 {
        (pair.a, pair.b, 0_usize, 2_usize)
    } else {
        (pair.b, pair.a, 2, 0)
    };
    let along_u = k == first + 1;
    let free = if along_u { first } else { first + 1 };
    let place = |s: f64| if along_u { (s, bound) } else { (bound, s) };
    let along = |s: f64| {
        let (u, v) = place(s);
        let (p, su, sv) = edged.point_d1_at(u, v, tol).ok()?;
        let (suu, _, svv) = edged.d2_at(u, v, tol).ok()?;
        Some(if along_u { (p, su, suu) } else { (p, sv, svv) })
    };
    let start = [from[free], from[rest], from[rest + 1]];
    let found = crate::curve_surface::touch_solved(along, other, start, tol)?;
    let (u, v) = place(found[0]);
    let within = |surface: &SurfaceGeometry, (u, v): (f64, f64)| {
        let ((ua, ub), (va, vb)) = surface.domain();
        [
            (u, ua, ub, surface.is_periodic_u()),
            (v, va, vb, surface.is_periodic_v()),
        ]
        .iter()
        .all(|&(x, lo, hi, periodic)| {
            periodic || (x >= lo - tol.parametric() && x <= hi + tol.parametric())
        })
    };
    if !within(edged, (u, v)) || !within(other, (found[1], found[2])) {
        return None;
    }
    let (p, su, sv) = edged.point_d1_at(u, v, tol).ok()?;
    let (suu, _, svv) = edged.d2_at(u, v, tol).ok()?;
    let (q, ou, ov) = other.point_d1_at(found[1], found[2], tol).ok()?;
    if p.distance(q) > tol.confusion() {
        return None;
    }
    // Along the edge the gap to the other surface grows as `g·s²/2`, `g`
    // the edge's normal curvature less the other surface's along it, and
    // off the edge as `a·t`, `a` the edge's cross derivative along the
    // normal. So the section stands `−g·s² / (2·a)` across the edge: inside
    // where that has the sign of the way into the domain. Where `g` or `a`
    // vanishes the contact is of higher order, and its side is not read.
    let (c, c2, across) = if along_u {
        (su, suu, sv)
    } else {
        (sv, svv, su)
    };
    let (ouu, ouv, ovv) = other.d2_at(found[1], found[2], tol).ok()?;
    let normal = ou.cross(ov);
    let metric = [[ou.dot(ou), ou.dot(ov)], [ou.dot(ov), ov.dot(ov)]];
    let [cu, cv] = solve::linear_fixed(metric, [c.dot(ou), c.dot(ov)])?;
    let bend = (c2 - (ouu * (cu * cu) + ouv * (2.0 * cu * cv) + ovv * (cv * cv))).dot(normal);
    let ((ua, ub), (va, vb)) = edged.domain();
    let (lo, hi) = if along_u { (va, vb) } else { (ua, ub) };
    let inward = if (bound - hi).abs() <= (bound - lo).abs() {
        -1.0
    } else {
        1.0
    };
    let side = -bend * across.dot(normal) * inward;
    if side <= 0.0 || side.is_nan() {
        return None;
    }
    let mut at = [0.0; 4];
    at[first] = u;
    at[first + 1] = v;
    at[rest] = found[1];
    at[rest + 1] = found[2];
    Some(at)
}

/// The sine of the shallowest crossing angle the marcher will follow.
///
/// One microradian, and the number is set by the *correction*, not by taste.
/// `correct` accepts a residual up to the confusion tolerance, so the two
/// parameter points of a contact can disagree by that much in space, and on
/// coincident or near-coincident surfaces, that disagreement shows up as a
/// spurious angle between the two computed normals of about the residual over
/// the local feature size. A gate below that floor reads the correction's own
/// noise as a direction and marches along it: identical spheres come back as
/// confident little curves that exist nowhere but in rounding.
///
/// So below this angle the marcher cannot tell an ultra-shallow crossing from
/// coincidence, and refuses both rather than guessing. A genuine crossing
/// shallower than a microradian is also one the Newton correction cannot
/// reliably follow (its travel constraint becomes numerically dependent on
/// the surface-gap rows at exactly the same rate), so the gate refuses what
/// could not have been followed anyway.
const SHALLOWEST: f64 = 1e-6;

/// The direction the intersection runs at a contact.
///
/// The cross product of the two normals: the one direction lying in both
/// tangent planes. `None` where the normals are parallel to within
/// [`SHALLOWEST`]: the surfaces are tangent or coincident there, and the
/// intersection has no direction the marcher can trust.
pub(super) fn tangent_at(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    at: Contact,
    tol: Tolerances,
) -> Option<Vector> {
    let na = normal_at(a, at.on_a, tol)?;
    let nb = normal_at(b, at.on_b, tol)?;
    tangent_of_normals(na, nb, tol)
}

/// The intersection tangent from the two unit normals, as [`tangent_at`]
/// decides it.
pub(super) fn tangent_of_normals(na: Vector, nb: Vector, tol: Tolerances) -> Option<Vector> {
    let cross = na.cross(nb);
    let length = cross.magnitude();
    // The decision is made through intervals rather than a bare compare:
    // each normal component carries the correction's stated residual as an
    // uncertainty, the squared cross magnitude is computed as an enclosure,
    // and only a crossing *certainly* above the floor is followed. A sine
    // inside the enclosure's undecided band is exactly the case the floor
    // exists for (the correction's own noise masquerading as an angle),
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
pub(super) fn normal_at(
    surface: &SurfaceGeometry,
    at: (f64, f64),
    tol: Tolerances,
) -> Option<Vector> {
    let (du, dv) = surface.d1_at(at.0, at.1, tol).ok()?;
    unit_normal(du, dv, tol)
}

/// The unit normal of two tangents, `None` where they are parallel.
pub(super) fn unit_normal(du: Vector, dv: Vector, tol: Tolerances) -> Option<Vector> {
    let cross = du.cross(dv);
    let length = cross.magnitude();
    if length <= tol.confusion() {
        return None;
    }
    Some(cross * (1.0 / length))
}

/// A state [`correct`] measured, the system's Jacobian and residual there,
/// and the first surface's two derivatives.
pub(super) type Measured = ([f64; 4], [[f64; 4]; 4], [f64; 4], [Vector; 2]);

/// Bring a parameter guess onto both surfaces.
///
/// Three equations say the two surface points coincide; the fourth says how far
/// along the direction of travel to land. Without that fourth the system is
/// underdetermined (its solution set *is* the curve), and Newton would wander
/// along it instead of converging to a point.
///
/// `constraint` is `(anchor, direction, distance)`. Without one, the guess
/// itself is used as the anchor and the direction is the intersection tangent,
/// which is what seeding wants: land anywhere on the curve near here. The
/// point it lands on stands within `within` of the section as far as the
/// solve can tell, as [`crate::walk::settling_residual`] sets out: where the
/// surfaces cross at a shallow angle the residual is taken further down.
pub(super) fn correct(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    start: [f64; 4],
    guess: Point,
    constraint: Option<(Point, Vector, f64)>,
    within: f64,
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

    // The state and Jacobian last measured, and the first surface's
    // derivatives there: what the point moves by with each unknown.
    let last: core::cell::Cell<Option<Measured>> = core::cell::Cell::new(None);
    let system = |x: &[f64; 4]| {
        last.set(None);
        let (ua, va) = clamp(a, x[0], x[1]);
        let (ub, vb) = clamp(b, x[2], x[3]);
        // Where either surface cannot be evaluated the residual is
        // infinite, so the damped step backs off rather than reading a
        // made-up point as a root.
        let (Ok((pa, au, av)), Ok((pb, bu, bv))) =
            (a.point_d1_at(ua, va, tol), b.point_d1_at(ub, vb, tol))
        else {
            return ([f64::INFINITY; 4], [[0.0; 4]; 4]);
        };

        let gap = pa - pb;
        let residual = [gap.x, gap.y, gap.z, (pa - anchor).dot(along) - reach];
        let jacobian = [
            [au.x, av.x, -bu.x, -bv.x],
            [au.y, av.y, -bu.y, -bv.y],
            [au.z, av.z, -bu.z, -bv.z],
            [au.dot(along), av.dot(along), 0.0, 0.0],
        ];
        last.set(Some((*x, jacobian, residual, [au, av])));
        (residual, jacobian)
    };

    let criteria = solve::Criteria {
        residual: tol.confusion() * 0.01,
        step: tol.parametric(),
        max_iterations: 40,
    };
    let mut found = solve::newton_system_fixed(system, start, criteria).ok()?;
    if found.1 > tol.confusion() {
        return None;
    }
    if last.get().is_none_or(|(x, _, _, _)| x != found.0) {
        system(&found.0);
    }
    if let Some((_, jacobian, residual, [au, av])) = last.get() {
        let gradient = [au, av, Vector::ZERO, Vector::ZERO];
        if let Some(settled) = crate::walk::settle(
            system,
            (found.0, &residual),
            (&jacobian, &gradient),
            within,
            criteria,
        ) {
            found.0 = settled;
        }
    }
    let (ua, va) = clamp(a, found.0[0], found.0[1]);
    let (ub, vb) = clamp(b, found.0[2], found.0[3]);
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
pub(crate) fn clamp(surface: &SurfaceGeometry, u: f64, v: f64) -> (f64, f64) {
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
/// The band is a fraction of the domain's own span: a walk stalls within a
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
pub(super) fn outside(surface: &SurfaceGeometry, at: (f64, f64), tol: Tolerances) -> bool {
    let ((ua, ub), (va, vb)) = surface.domain();
    let past = |x: f64, lo: f64, hi: f64, periodic: bool| {
        !periodic && (x <= lo + tol.parametric() || x >= hi - tol.parametric())
    };
    past(at.0, ua, ub, surface.is_periodic_u()) || past(at.1, va, vb, surface.is_periodic_v())
}

/// A surface's rough size, for spacing seeds.
pub(crate) fn span(surface: &SurfaceGeometry) -> f64 {
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
