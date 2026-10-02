//! Tangential contact: two surfaces that touch along a curve without
//! crossing it.
//!
//! Where two surfaces cross, the curve they share runs along the cross
//! product of their normals, and the crossing walker in [`crate::march`]
//! follows that. Where they touch, the normals are parallel all along the
//! curve, the cross product vanishes, and the crossing walker stalls in
//! fragments. The contact is still a curve, with a description of its own:
//! the locus where the two surfaces meet with parallel normals. This module
//! follows that locus.
//!
//! # Prediction
//!
//! At a contact the two surfaces share a tangent plane, and along the
//! contact their normal curvatures agree: both are the normal curvature of
//! the one curve they share. The contact therefore runs along the null
//! direction of the difference of the two second fundamental forms, taken
//! in the common tangent plane. Where that difference has no distinguished
//! null direction (both its eigenvalues of one size), the surfaces touch at
//! an isolated point or osculate in every direction, and there is no curve
//! to predict along.
//!
//! # Correction
//!
//! A Newton step on four conditions in the four parameters: the second
//! surface's point lies on the first surface's normal line through its own
//! point (two conditions), the two normals agree across the contact (one
//! condition; their disagreement along the contact is the null direction's
//! and carries no information), and the point lies one step ahead along the
//! prediction. The gap along the normal is not a condition: at a contact it
//! touches zero without crossing it, so it has no root a Newton step could
//! find. It is measured instead, and a corrected point whose gap exceeds the
//! acceptance is where the surfaces part, so the trace stops there.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Surface as _, SurfaceGeometry};
use ogeom_math::{Point, Vector, solve};

use crate::march::{BRANCH_POINT_SINE, Contact, Marching, Stopped, Traced, clamp, span};

/// How far apart the two points of a corrected contact may sit, in units of
/// the confusion distance.
const ACCEPT: f64 = 100.0;

/// The ratio of the two eigenvalues of the second forms' difference below
/// which the smaller one names the contact's direction.
const DISTINCT: f64 = 0.25;

/// A surface's point, derivatives and unit normal at one parameter.
#[derive(Debug, Clone, Copy)]
struct Jet {
    p: Point,
    du: Vector,
    dv: Vector,
    duu: Vector,
    duv: Vector,
    dvv: Vector,
    n: Vector,
    /// The length of `du x dv`.
    area: f64,
}

fn jet(surface: &SurfaceGeometry, at: (f64, f64), tol: Tolerances) -> Option<Jet> {
    let p = surface.point_at(at.0, at.1, tol).ok()?;
    let (du, dv) = surface.d1_at(at.0, at.1, tol).ok()?;
    let (duu, duv, dvv) = surface.d2_at(at.0, at.1, tol).ok()?;
    let cross = du.cross(dv);
    let area = cross.magnitude();
    if area.is_nan() || area <= tol.confusion() {
        return None;
    }
    Some(Jet {
        p,
        du,
        dv,
        duu,
        duv,
        dvv,
        n: cross * (1.0 / area),
        area,
    })
}

impl Jet {
    /// The parameter step whose image is the tangent vector `e`.
    fn lift(&self, e: Vector) -> Option<(f64, f64)> {
        let (g11, g12, g22) = (
            self.du.dot(self.du),
            self.du.dot(self.dv),
            self.dv.dot(self.dv),
        );
        let det = g11.mul_add(g22, -(g12 * g12));
        if det <= f64::MIN_POSITIVE {
            return None;
        }
        let (r1, r2) = (e.dot(self.du), e.dot(self.dv));
        Some((
            r1.mul_add(g22, -(r2 * g12)) / det,
            g11.mul_add(r2, -(g12 * r1)) / det,
        ))
    }

    /// The second fundamental form against the normal `n`, in the
    /// orthonormal tangent basis `e1`, `e2`: `[II(e1,e1), II(e1,e2),
    /// II(e2,e2)]`.
    fn second_form(&self, n: Vector, e1: Vector, e2: Vector) -> Option<[f64; 3]> {
        let (a1, b1) = self.lift(e1)?;
        let (a2, b2) = self.lift(e2)?;
        let (l, m, nn) = (n.dot(self.duu), n.dot(self.duv), n.dot(self.dvv));
        let form =
            |a: f64, b: f64, c: f64, d: f64| a * c * l + a.mul_add(d, b * c) * m + b * d * nn;
        Some([
            form(a1, b1, a1, b1),
            form(a1, b1, a2, b2),
            form(a2, b2, a2, b2),
        ])
    }

    /// The derivatives of the unit normal along `u` and `v`.
    fn normal_rates(&self) -> (Vector, Vector) {
        let rate = |dn: Vector| (dn - self.n * self.n.dot(dn)) * (1.0 / self.area);
        (
            rate(self.duu.cross(self.dv) + self.du.cross(self.duv)),
            rate(self.duv.cross(self.dv) + self.du.cross(self.dvv)),
        )
    }
}

/// The second surface's normal turned to agree with the first's.
fn aligned(ja: &Jet, jb: &Jet) -> (Vector, f64) {
    let sign = if ja.n.dot(jb.n) < 0.0 { -1.0 } else { 1.0 };
    (jb.n * sign, sign)
}

/// The direction the contact runs at a point where the surfaces touch, and
/// the first surface's normal curvature along it.
///
/// The null direction of the difference of the two second fundamental
/// forms, signed to continue `previous`. Where the difference has no
/// distinguished null direction, `previous` carried into the tangent plane
/// stands in if there is one, and `None` otherwise.
fn contact_direction(ja: &Jet, jb: &Jet, previous: Option<Vector>) -> Option<(Vector, f64)> {
    let n = ja.n;
    let (nb, _) = aligned(ja, jb);
    let flat_u = ja.du - n * n.dot(ja.du);
    let e1 = flat_u * (1.0 / flat_u.magnitude());
    if !e1.x.is_finite() {
        return None;
    }
    let e2 = n.cross(e1);
    let fa = ja.second_form(n, e1, e2)?;
    let fb = jb.second_form(nb, e1, e2)?;
    let (p, q, r) = (fa[0] - fb[0], fa[1] - fb[1], fa[2] - fb[2]);
    let half_sum = f64::midpoint(p, r);
    let radius = ((p - r) * 0.5).hypot(q);
    let (low, high) = (half_sum - radius, half_sum + radius);
    let (small, large) = if low.abs() <= high.abs() {
        (low, high)
    } else {
        (high, low)
    };
    let carried = previous.and_then(|d| {
        let flat = d - n * n.dot(d);
        let m = flat.magnitude();
        (m > f64::MIN_POSITIVE).then(|| flat * (1.0 / m))
    });
    let direction = if large.abs() > f64::MIN_POSITIVE && small.abs() <= DISTINCT * large.abs() {
        // The eigenvector of `small`, from whichever row of the shifted
        // matrix is the better conditioned.
        let first = (q, small - p);
        let second = (small - r, q);
        let (x, y) = if first.0.hypot(first.1) >= second.0.hypot(second.1) {
            first
        } else {
            second
        };
        let m = x.hypot(y);
        if m <= f64::MIN_POSITIVE {
            return None;
        }
        let t = e1 * (x / m) + e2 * (y / m);
        match carried {
            Some(d) if d.dot(t) < 0.0 => -t,
            _ => t,
        }
    } else {
        carried?
    };
    let (x, y) = (direction.dot(e1), direction.dot(e2));
    let bend = (fa[0] * x).mul_add(x, (2.0 * fa[1] * x).mul_add(y, fa[2] * y * y));
    Some((direction, bend))
}

/// The fourth condition of a correction: a step along the prediction, or a
/// parameter held on a bound.
#[derive(Debug, Clone, Copy)]
enum Pin {
    Step {
        anchor: Point,
        along: Vector,
        reach: f64,
    },
    Bound {
        index: usize,
        value: f64,
    },
}

/// A corrected contact: the raw parameters the solve reached (which may lie
/// past a bound), the contact itself, and its measured gap and normal sine.
#[derive(Debug, Clone, Copy)]
struct Corrected {
    raw: [f64; 4],
    contact: Contact,
    gap: f64,
    sine: f64,
}

/// Bring a guess onto the contact locus. `across` is the tangent direction
/// across the contact, along which the normals are made to agree.
fn correct(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    start: [f64; 4],
    across: Vector,
    pin: Pin,
    tol: Tolerances,
) -> Option<Corrected> {
    let system = |x: &[f64; 4]| {
        let (ua, va) = clamp(a, x[0], x[1]);
        let (ub, vb) = clamp(b, x[2], x[3]);
        let (Some(ja), Some(jb)) = (jet(a, (ua, va), tol), jet(b, (ub, vb), tol)) else {
            return ([f64::INFINITY; 4], [[0.0; 4]; 4]);
        };
        let d = ja.p - jb.p;
        let (lu, lv) = (ja.du.magnitude(), ja.dv.magnitude());
        let (nb, sign) = aligned(&ja, &jb);
        let (na_u, na_v) = ja.normal_rates();
        let (nb_u, nb_v) = jb.normal_rates();
        let (last, last_row) = match pin {
            Pin::Step {
                anchor,
                along,
                reach,
            } => (
                (ja.p - anchor).dot(along) - reach,
                [ja.du.dot(along), ja.dv.dot(along), 0.0, 0.0],
            ),
            Pin::Bound { index, value } => {
                let mut row = [0.0; 4];
                row[index] = 1.0;
                (x[index] - value, row)
            }
        };
        let residual = [
            d.dot(ja.du) / lu,
            d.dot(ja.dv) / lv,
            (nb - ja.n).dot(across),
            last,
        ];
        let jacobian = [
            [
                (ja.du.dot(ja.du) + d.dot(ja.duu)) / lu,
                (ja.dv.dot(ja.du) + d.dot(ja.duv)) / lu,
                -jb.du.dot(ja.du) / lu,
                -jb.dv.dot(ja.du) / lu,
            ],
            [
                (ja.du.dot(ja.dv) + d.dot(ja.duv)) / lv,
                (ja.dv.dot(ja.dv) + d.dot(ja.dvv)) / lv,
                -jb.du.dot(ja.dv) / lv,
                -jb.dv.dot(ja.dv) / lv,
            ],
            [
                -na_u.dot(across),
                -na_v.dot(across),
                nb_u.dot(across) * sign,
                nb_v.dot(across) * sign,
            ],
            last_row,
        ];
        (residual, jacobian)
    };
    let criteria = solve::Criteria {
        residual: tol.confusion() * 1e-3,
        step: tol.parametric() * 1e-3,
        max_iterations: 40,
    };
    let (raw, ..) = solve::newton_system_fixed(system, start, criteria).ok()?;
    let on_a = clamp(a, raw[0], raw[1]);
    let on_b = clamp(b, raw[2], raw[3]);
    let ja = jet(a, on_a, tol)?;
    let jb = jet(b, on_b, tol)?;
    Some(Corrected {
        raw,
        contact: Contact {
            on_a,
            on_b,
            point: ja.p,
        },
        gap: ja.p.distance(jb.p),
        sine: ja.n.cross(jb.n).magnitude(),
    })
}

/// The first non-periodic parameter a solve carried past its bound, as the
/// index and the bound, with how far along the move from `from` the bound
/// lies. `None` while every parameter is inside.
fn crossed_bound(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    from: [f64; 4],
    to: [f64; 4],
    tol: Tolerances,
) -> Option<(usize, f64)> {
    let mut first: Option<(f64, usize, f64)> = None;
    for (surface, offset) in [(a, 0_usize), (b, 2)] {
        let ((ua, ub), (va, vb)) = surface.domain();
        for (k, (lo, hi), periodic) in [
            (offset, (ua, ub), surface.is_periodic_u()),
            (offset + 1, (va, vb), surface.is_periodic_v()),
        ] {
            if periodic || !(hi - lo).is_finite() {
                continue;
            }
            let slack = tol.parametric();
            let bound = if to[k] < lo - slack {
                lo
            } else if to[k] > hi + slack {
                hi
            } else {
                continue;
            };
            let moved = to[k] - from[k];
            let fraction = if moved.abs() > f64::MIN_POSITIVE {
                ((bound - from[k]) / moved).clamp(0.0, 1.0)
            } else {
                0.0
            };
            if first.is_none_or(|(f, ..)| fraction < f) {
                first = Some((fraction, k, bound));
            }
        }
    }
    first.map(|(_, k, bound)| (k, bound))
}

/// A walk one way along the contact.
struct Walked {
    points: Vec<Point>,
    on_a: Vec<(f64, f64)>,
    on_b: Vec<(f64, f64)>,
    stopped: Stopped,
    /// The longest step taken.
    longest: f64,
}

/// The settings one trace walks with.
struct Walk<'s> {
    a: &'s SurfaceGeometry,
    b: &'s SurfaceGeometry,
    options: Marching,
    accept: f64,
    longest: f64,
    shortest: f64,
    tol: Tolerances,
}

impl Walk<'_> {
    /// Whether a corrected point is a contact.
    fn holds(&self, c: &Corrected) -> bool {
        c.gap <= self.accept && c.sine <= BRANCH_POINT_SINE
    }

    /// The step a curve bending at `curvature` affords within the chord.
    fn step_for(&self, curvature: f64) -> f64 {
        let curvature = curvature
            .abs()
            .max(1.0 / self.longest.max(f64::MIN_POSITIVE));
        (8.0 * self.options.chord / curvature)
            .sqrt()
            .clamp(self.shortest, self.longest)
    }

    fn one_way(&self, from: Contact, sense: f64) -> OgeomResult<Walked> {
        let (a, b, tol) = (self.a, self.b, self.tol);
        let mut out = Walked {
            points: vec![from.point],
            on_a: vec![from.on_a],
            on_b: vec![from.on_b],
            stopped: Stopped::RanOut,
            longest: 0.0,
        };
        let mut at = from;
        let mut previous: Option<(Vector, f64)> = None;
        let mut travelled = 0.0;
        while out.points.len() < self.options.max_points {
            ogeom_core::progress::checkpoint()?;
            let (Some(ja), Some(jb)) = (jet(a, at.on_a, tol), jet(b, at.on_b, tol)) else {
                out.stopped = Stopped::Stalled;
                break;
            };
            let Some((mut along, bend)) = contact_direction(&ja, &jb, previous.map(|(d, _)| d))
            else {
                out.stopped = Stopped::Stalled;
                break;
            };
            if previous.is_none() {
                along *= sense;
            }
            // The step: the chord's allowance for the sharper of the
            // surface's own bend along the contact and the contact's
            // turning over the last step, grown at most twofold a step.
            let mut step = match previous {
                Some((d, last)) => {
                    let turning = d.cross(along).magnitude().atan2(d.dot(along)) / last;
                    self.step_for(bend.abs().max(turning)).min(last * 2.0)
                }
                None => self.step_for(bend),
            };
            let across = ja.n.cross(along);
            let state = [at.on_a.0, at.on_a.1, at.on_b.0, at.on_b.1];
            let next = loop {
                let guess = {
                    let ta = ja.lift(along * step);
                    let tb = jb.lift(along * step);
                    match (ta, tb) {
                        (Some((da, ea)), Some((db, eb))) => {
                            [state[0] + da, state[1] + ea, state[2] + db, state[3] + eb]
                        }
                        _ => state,
                    }
                };
                let pin = Pin::Step {
                    anchor: at.point,
                    along,
                    reach: step,
                };
                if let Some(c) = correct(a, b, guess, across, pin, tol) {
                    if let Some((index, value)) = crossed_bound(a, b, state, c.raw, tol) {
                        // The contact runs off a surface's edge within this
                        // step: it ends on the edge, where the corrector
                        // lands it with the parameter held on the bound.
                        // A landing that fails is taken as a step too
                        // long, and the step is halved.
                        let pin = Pin::Bound { index, value };
                        if let Some(edge) = correct(a, b, state, across, pin, tol)
                            && self.holds(&edge)
                            && (edge.contact.point - at.point).dot(along) > -tol.confusion()
                            && edge.contact.point.distance(at.point) <= step * 2.0
                        {
                            break Err(Some(edge.contact));
                        }
                    } else if self.holds(&c) && c.contact.point.distance(at.point) > step * 0.5 {
                        break Ok(c.contact);
                    }
                }
                step *= 0.5;
                if step < self.shortest {
                    break Err(None);
                }
            };
            let next = match next {
                Ok(next) => next,
                Err(Some(edge)) => {
                    let taken = edge.point.distance(at.point);
                    if taken <= tol.confusion() {
                        let last = out.points.len() - 1;
                        out.points[last] = edge.point;
                        out.on_a[last] = edge.on_a;
                        out.on_b[last] = edge.on_b;
                    } else {
                        out.longest = out.longest.max(taken);
                        out.points.push(edge.point);
                        out.on_a.push(edge.on_a);
                        out.on_b.push(edge.on_b);
                    }
                    out.stopped = Stopped::LeftTheDomain;
                    break;
                }
                Err(None) => {
                    out.stopped = Stopped::Stalled;
                    break;
                }
            };
            let taken = next.point.distance(at.point);
            // Back where it began: a loop.
            if out.points.len() > 3
                && travelled > taken * 2.0
                && next.point.distance(from.point) <= taken
                && (from.point - at.point).dot(along) > 0.0
            {
                out.points.push(from.point);
                out.on_a.push(from.on_a);
                out.on_b.push(from.on_b);
                out.stopped = Stopped::Closed;
                break;
            }
            travelled += taken;
            out.longest = out.longest.max(taken);
            previous = Some((next.point - at.point, taken));
            out.points.push(next.point);
            out.on_a.push(next.on_a);
            out.on_b.push(next.on_b);
            at = next;
        }
        Ok(out)
    }
}

/// Trace the tangential contact of two surfaces through a seed.
///
/// The seed need only be near the contact: it is first corrected onto it.
/// The trace walks both ways from there (see the module documentation for
/// the predictor and the corrector) and stops where the contact closes,
/// runs off either surface's domain, or ends: where the surfaces part by
/// more than a hundred confusion distances, where their normals part by
/// more than a branch point's sine, or where the contact has no direction
/// to follow (an isolated touch, a degenerate chart). Every point it
/// returns was measured against both conditions.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// settings are unusable, the surfaces cross at the seed rather than touch,
/// or the seed corrects to no contact: the surfaces come near there without
/// touching, or touch at an isolated point.
pub fn trace_tangential(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    from: Contact,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Traced> {
    options.validate()?;
    let (Some(ja), Some(jb)) = (jet(a, from.on_a, tol), jet(b, from.on_b, tol)) else {
        ogeom_bail!(
            Construction,
            "the seed cannot be evaluated on both surfaces"
        );
    };
    let sine = ja.n.cross(jb.n).magnitude();
    if sine > BRANCH_POINT_SINE {
        ogeom_bail!(
            Construction,
            "the surfaces cross here at sine {sine}; tangential tracing wants a contact"
        );
    }
    let longest = span(a).max(span(b)) / 16.0;
    let walk = Walk {
        a,
        b,
        options,
        accept: tol.confusion() * ACCEPT,
        longest,
        shortest: (longest * 1e-6).max(tol.confusion() * 10.0),
        tol,
    };
    let Some((along, _)) = contact_direction(&ja, &jb, None) else {
        ogeom_bail!(
            Construction,
            "the surfaces touch at an isolated point here, or osculate in every \
             direction: the contact has no direction to follow"
        );
    };
    let start = [from.on_a.0, from.on_a.1, from.on_b.0, from.on_b.1];
    let pin = Pin::Step {
        anchor: from.point,
        along,
        reach: 0.0,
    };
    let seed = match correct(a, b, start, ja.n.cross(along), pin, tol) {
        Some(c) if walk.holds(&c) => c.contact,
        Some(c) => ogeom_bail!(
            Construction,
            "the surfaces come within {} of each other here, at normal sine {}, \
             without touching",
            c.gap,
            c.sine
        ),
        None => ogeom_bail!(Construction, "the seed does not correct onto a contact"),
    };

    let ahead = walk.one_way(seed, 1.0)?;
    if ahead.stopped == Stopped::Closed {
        return Ok(Traced {
            points: ahead.points,
            on_a: ahead.on_a,
            on_b: ahead.on_b,
            stopped: Stopped::Closed,
        });
    }
    let behind = walk.one_way(seed, -1.0)?;
    let steps = ahead.longest.max(behind.longest);
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
    let mut stopped = if ahead.stopped == Stopped::RanOut || behind.stopped == Stopped::RanOut {
        Stopped::RanOut
    } else if ahead.stopped == Stopped::Stalled || behind.stopped == Stopped::Stalled {
        Stopped::Stalled
    } else {
        Stopped::LeftTheDomain
    };
    // A loop cut at the seam of a patch that closes on itself without
    // being periodic: both walks end on the seam, at one point.
    if stopped == Stopped::LeftTheDomain && points.len() > 3 {
        let gap = points[0].distance(points[points.len() - 1]);
        if gap <= (steps * 2.0).max(tol.confusion() * 10.0) {
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

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use ogeom_geom::{CylinderSurface, PlaneSurface, SphereSurface};
    use ogeom_math::{Cylinder, Direction, Frame, Plane, Sphere};

    const T: Tolerances = Tolerances::millimetres();

    fn ball(radius: f64) -> SurfaceGeometry {
        SphereSurface::new(Sphere::centred(Point::ORIGIN, radius, T).unwrap()).into()
    }

    /// The parameters of a point on a surface, by a coarse search refined
    /// by projection.
    fn locate(surface: &SurfaceGeometry, p: Point) -> (f64, f64) {
        let ((ua, ub), (va, vb)) = surface.domain();
        let mut best = ((ua, va), f64::INFINITY);
        for i in 0..=32 {
            for j in 0..=32 {
                let u = ua + (ub - ua) * f64::from(i) / 32.0;
                let v = va + (vb - va) * f64::from(j) / 32.0;
                let d = surface.point_at(u, v, T).unwrap().distance(p);
                if d < best.1 {
                    best = ((u, v), d);
                }
            }
        }
        crate::march::nearest_on(surface, best.0, p, T).unwrap().0
    }

    fn seed(a: &SurfaceGeometry, b: &SurfaceGeometry, p: Point) -> Contact {
        Contact {
            on_a: locate(a, p),
            on_b: locate(b, p),
            point: p,
        }
    }

    #[test]
    fn a_ball_on_a_plane_touches_at_a_point_and_is_refused() {
        // The second forms differ by the same amount in every direction:
        // no direction is the contact's, because the contact is a point.
        let floor: SurfaceGeometry = PlaneSurface::over(
            Plane::through(
                Point::new(0.6, 0.0, 0.8),
                Direction::new(Vector::new(0.6, 0.0, 0.8), T).unwrap(),
            ),
            (-4.0, 4.0),
            (-4.0, 4.0),
        )
        .unwrap()
        .into();
        let ball = ball(1.0);
        let at = Point::new(0.6, 0.0, 0.8);
        let err = trace_tangential(
            &floor,
            &ball,
            seed(&floor, &ball, at),
            Marching::default(),
            T,
        )
        .unwrap_err();
        assert!(err.to_string().contains("isolated point"), "{err}");
    }

    #[test]
    fn a_ball_short_of_its_drum_is_refused_with_the_gap() {
        // A unit ball in a drum a thousandth wider: the closest approach is
        // the equator, a curve, but the surfaces never meet along it.
        let drum: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 1.001, T).unwrap(), (-2.0, 2.0))
                .unwrap()
                .into();
        let ball = ball(1.0);
        let at = Point::new(1.0005, 0.0, 0.0);
        let err = trace_tangential(&drum, &ball, seed(&drum, &ball, at), Marching::default(), T)
            .unwrap_err();
        assert!(err.to_string().contains("without touching"), "{err}");
    }

    #[test]
    fn a_ball_in_its_drum_is_traced_on_the_contact_conditions() {
        // Every point the trace returns meets the conditions it was
        // corrected onto: on both surfaces, normals agreeing.
        let drum: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 1.0, T).unwrap(), (-2.0, 2.0))
                .unwrap()
                .into();
        let ball = ball(1.0);
        let at = Point::new(0.6, 0.8, 0.0);
        let traced =
            trace_tangential(&drum, &ball, seed(&drum, &ball, at), Marching::default(), T).unwrap();
        assert_eq!(traced.stopped, Stopped::Closed);
        for (i, p) in traced.points.iter().enumerate() {
            let (ja, jb) = (
                jet(&drum, traced.on_a[i], T).unwrap(),
                jet(&ball, traced.on_b[i], T).unwrap(),
            );
            assert!(ja.p.distance(*p) <= 1e-12 && jb.p.distance(*p) <= 1e-9);
            assert!(ja.n.cross(jb.n).magnitude() <= 1e-9);
            assert!(p.z.abs() <= 1e-9 && (p.x.hypot(p.y) - 1.0).abs() <= 1e-9);
        }
    }
}
