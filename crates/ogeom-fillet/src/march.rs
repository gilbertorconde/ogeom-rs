//! The marching blend: a rolling ball followed by solving where it touches.
//!
//! The analytic seats — an edge between two planes, a plane and a cylinder in
//! the configurations that give a torus — are built from their own closed
//! forms elsewhere. This is the general one, and the formulation matters more
//! than the marching does.
//!
//! # Four unknowns, four equations
//!
//! The obvious construction intersects the two supports' *offset* surfaces to
//! get the spine, projects back onto each support for the tangency points, and
//! skins the arcs between them. It works on paper and is the wrong shape: the
//! tangency curves arrive by projection, so the legs' pcurves are **fitted**,
//! and a fitted pcurve on a support is exactly what the boolean cannot treat
//! as same-domain later.
//!
//! So the section's two endpoints are solved for directly. The unknowns are
//! `(u₁, v₁)` on the first support and `(u₂, v₂)` on the second — the two
//! points where the ball touches. Three equations say the ball's centre is the
//! same point computed from either side,
//!
//! > `P₁ + r·n₁ = P₂ + r·n₂`
//!
//! and the fourth ties the section to a **guide** — the edge being blended, or
//! any curve running along the seat — by requiring it to lie in the plane
//! through the guide point normal to the guide's tangent.
//!
//! What comes out is worth the change: the tangency curves emerge **in the
//! supports' own parameters**, so the legs' pcurves are exact by construction
//! rather than fitted.
//!
//! # Marched by the shared walker
//!
//! The guide's parameter joins the unknowns as a fifth, which makes the system
//! four equations in five — a curve — and that is what
//! [`ogeom_intersect::walk`] follows. So the step control, the stall reporting
//! and the closure test are the intersector's own, inherited rather than
//! written a second time, and the step is set by the sag of the *tangency
//! curve* it is walking.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_intersect::walk::Condition;
use ogeom_intersect::{Marching, Stopped};
use ogeom_math::{Point, Vector, solve};

/// Which side of each support the ball rolls on.
///
/// A rolling ball sits at `P + s·r·n` for one sign `s` per support: outward
/// for a convex seat, inward for a concave one, and one of each where the
/// blend runs along a step. The pair is not guessed from the normals —
/// normals cannot tell a step from a slot — but tried, and the combination
/// that gives a ball genuinely touching both is the seat.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Sides {
    /// The sign on the first support's normal.
    pub first: i8,
    /// The sign on the second's.
    pub second: i8,
}

impl Sides {
    /// Every combination, for a caller that does not know the seat.
    const ALL: [Self; 4] = [
        Self {
            first: 1,
            second: 1,
        },
        Self {
            first: 1,
            second: -1,
        },
        Self {
            first: -1,
            second: 1,
        },
        Self {
            first: -1,
            second: -1,
        },
    ];
}

/// Why a marched blend stopped.
///
/// The list is the case checklist the formulation gives for free, and each
/// one is a different thing for a caller to do about it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlendStop {
    /// The seat closed on itself — a blend all the way round a rim.
    Closed,
    /// The section reached the boundary of the first support.
    LeftTheFirstSupport,
    /// Of the second.
    LeftTheSecondSupport,
    /// Of both at once, which is a corner rather than a run-out.
    LeftBothSupports,
    /// The guide ran out before either support did.
    RanPastTheGuide,
    /// The two tangency points collapsed onto each other: the radius is too
    /// large for the local geometry and the ball is not seated but wedged.
    SectionCollapsed,
    /// The correction stopped converging — a tangency or a singular point on
    /// one of the supports.
    Stalled,
    /// The step ceiling, which means the answer is *incomplete* rather than
    /// finished.
    RanOut,
}

/// One marched blend: where the ball touched, and where its centre went.
#[derive(Debug, Clone)]
pub struct MarchedBlend {
    /// The ball's centre at each station — the blend's spine.
    pub spine: Vec<Point>,
    /// Where it touched the first support, in that support's own parameters.
    ///
    /// The point of the whole formulation: these are solved for, not
    /// projected, so a pcurve fitted through them is a pcurve of the curve
    /// itself rather than of a projection of it.
    pub on_first: Vec<(f64, f64)>,
    /// And the second's.
    pub on_second: Vec<(f64, f64)>,
    /// The touch points on the first support, in space.
    pub touch_first: Vec<Point>,
    /// And on the second.
    pub touch_second: Vec<Point>,
    /// Which side of each support the ball rolled on.
    pub sides: Sides,
    /// Why it stopped.
    pub stopped: BlendStop,
}

impl MarchedBlend {
    /// How many stations the march produced.
    #[must_use]
    pub fn len(&self) -> usize {
        self.spine.len()
    }

    /// Whether it produced none.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.spine.is_empty()
    }

    /// Whether the march finished rather than being truncated.
    #[must_use]
    pub const fn complete(&self) -> bool {
        !matches!(self.stopped, BlendStop::RanOut)
    }
}

/// Follow a rolling ball of `radius` seated between two supports, guided by
/// `guide`.
///
/// The guide is the curve the sections stand square to — the edge being
/// blended, or any curve running along the seat. It decides *where* the
/// sections are, not what they are: the ball's own contact conditions decide
/// that, and a guide that is merely near the seat gives the same blend as one
/// exactly on it.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// radius or the settings are unusable.
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if no seating of
/// the ball can be found at the guide's start — which is the honest answer for
/// a radius the corner cannot hold, and names that rather than marching a
/// system that is not solved.
pub fn march_blend(
    first: &SurfaceGeometry,
    second: &SurfaceGeometry,
    radius: f64,
    guide: &Curve,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<MarchedBlend> {
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(Construction, "a blend of radius {radius} rounds nothing");
    }
    options.validate()?;
    let (start, sides) = seat(first, second, radius, guide, tol)?;
    let contact = BallContact {
        first,
        second,
        radius,
        guide,
        sides,
    };
    let walked = ogeom_intersect::walk::follow(&contact, &start, options, tol)?;

    let mut blend = MarchedBlend {
        spine: Vec::with_capacity(walked.states.len()),
        on_first: Vec::with_capacity(walked.states.len()),
        on_second: Vec::with_capacity(walked.states.len()),
        touch_first: Vec::with_capacity(walked.states.len()),
        touch_second: Vec::with_capacity(walked.states.len()),
        sides,
        stopped: BlendStop::Stalled,
    };
    for x in &walked.states {
        let (Ok(p1), Ok(p2)) = (
            first.point_at(x[0], x[1], tol),
            second.point_at(x[2], x[3], tol),
        ) else {
            continue;
        };
        let Some(centre) = contact.centre(x, tol) else {
            continue;
        };
        blend.spine.push(centre);
        blend.on_first.push((x[0], x[1]));
        blend.on_second.push((x[2], x[3]));
        blend.touch_first.push(p1);
        blend.touch_second.push(p2);
    }
    blend.stopped = why(&contact, &walked, tol);
    Ok(blend)
}

/// Which of the walker's reasons this is, in the blend's own vocabulary.
fn why(
    contact: &BallContact<'_>,
    walked: &ogeom_intersect::walk::Walked,
    tol: Tolerances,
) -> BlendStop {
    let Some(last) = walked.states.last() else {
        return BlendStop::Stalled;
    };
    // A section whose two ends have met is not a section: the ball is wedged
    // rather than seated, which is what a radius too large for the local
    // geometry looks like from inside the march.
    let collapsed = matches!(
        (
            contact.first.point_at(last[0], last[1], tol),
            contact.second.point_at(last[2], last[3], tol),
        ),
        (Ok(p1), Ok(p2)) if p1.distance(p2) <= contact.radius * 1e-3
    );
    match walked.stopped {
        Stopped::Closed => BlendStop::Closed,
        Stopped::RanOut => BlendStop::RanOut,
        Stopped::Stalled if collapsed => BlendStop::SectionCollapsed,
        Stopped::Stalled => BlendStop::Stalled,
        Stopped::LeftTheDomain => {
            let left_first = at_edge(contact.first, (last[0], last[1]));
            let left_second = at_edge(contact.second, (last[2], last[3]));
            match (left_first, left_second) {
                (true, true) => BlendStop::LeftBothSupports,
                (true, false) => BlendStop::LeftTheFirstSupport,
                (false, true) => BlendStop::LeftTheSecondSupport,
                (false, false) => BlendStop::RanPastTheGuide,
            }
        }
    }
}

/// Where the ball first sits, and which side of each support it sits on.
///
/// Tried rather than assumed: the four sign combinations are each given a
/// Newton solve from the guide's own start, and the one that seats a ball
/// touching two *distinct* points wins. A radius the corner cannot hold seats
/// none of them, and that is what the refusal says.
fn seat(
    first: &SurfaceGeometry,
    second: &SurfaceGeometry,
    radius: f64,
    guide: &Curve,
    tol: Tolerances,
) -> OgeomResult<([f64; 5], Sides)> {
    let (lo, hi) = guide.domain();
    let at = f64::midpoint(lo, hi);
    let anchor = guide.point_at(at, tol)?;
    let near_first = ogeom_algo::project_on_surface(first, anchor, 24, tol)?;
    let near_second = ogeom_algo::project_on_surface(second, anchor, 24, tol)?;
    let start = [
        near_first.parameters.0,
        near_first.parameters.1,
        near_second.parameters.0,
        near_second.parameters.1,
        at,
    ];

    for sides in Sides::ALL {
        let contact = BallContact {
            first,
            second,
            radius,
            guide,
            sides,
        };
        // The guide parameter is held: at the seed there is nothing yet to
        // march along, only a section to find.
        let system = |x: &[f64]| {
            let mut full = [x[0], x[1], x[2], x[3], at];
            contact.clamp(&mut full);
            let (mut residual, jacobian) = contact
                .system(&full, tol)
                .unwrap_or_else(|| (vec![0.0; 4], vec![vec![0.0; 5]; 4]));
            residual.truncate(4);
            let jacobian = jacobian
                .into_iter()
                .take(4)
                .map(|row| row[..4].to_vec())
                .collect();
            (residual, jacobian)
        };
        let criteria = solve::Criteria {
            residual: tol.confusion() * 0.01,
            step: tol.parametric(),
            max_iterations: 80,
        };
        let Ok(found) = solve::newton_system(system, &start[..4], criteria) else {
            continue;
        };
        if found.residual > tol.confusion() {
            continue;
        }
        let mut x = [
            found.value[0],
            found.value[1],
            found.value[2],
            found.value[3],
            at,
        ];
        contact.clamp(&mut x);
        let (Ok(p1), Ok(p2)) = (
            first.point_at(x[0], x[1], tol),
            second.point_at(x[2], x[3], tol),
        ) else {
            continue;
        };
        // Two distinct touches, and a ball genuinely of the stated radius.
        if p1.distance(p2) <= radius * 1e-3 {
            continue;
        }
        let Some(centre) = contact.centre(&x, tol) else {
            continue;
        };
        if (centre.distance(p1) - radius).abs() > tol.confusion() * 10.0
            || (centre.distance(p2) - radius).abs() > tol.confusion() * 10.0
        {
            continue;
        }
        return Ok((x, sides));
    }
    ogeom_bail!(
        NotDone,
        "no ball of radius {radius} seats between these supports at the \
         guide's start; either the corner cannot hold one or the guide does \
         not run along the seat"
    );
}

/// The rolling ball's contact, as a condition the walker can follow.
struct BallContact<'s> {
    first: &'s SurfaceGeometry,
    second: &'s SurfaceGeometry,
    radius: f64,
    guide: &'s Curve,
    sides: Sides,
}

impl BallContact<'_> {
    /// The ball's centre, computed from the first support's side.
    fn centre(&self, x: &[f64], tol: Tolerances) -> Option<Point> {
        let p = self.first.point_at(x[0], x[1], tol).ok()?;
        let n = unit_normal(self.first, x[0], x[1], tol)?;
        Some(p + n * (f64::from(self.sides.first) * self.radius))
    }
}

impl Condition for BallContact<'_> {
    fn unknowns(&self) -> usize {
        5
    }

    fn position(&self, x: &[f64], tol: Tolerances) -> Option<Point> {
        // The tangency point on the first support: the curve being walked is
        // the first leg, so the step control measures *its* sag.
        self.first.point_at(x[0], x[1], tol).ok()
    }

    fn position_gradient(&self, x: &[f64], tol: Tolerances) -> Option<Vec<Vector>> {
        let (du, dv) = self.first.d1_at(x[0], x[1], tol).ok()?;
        Some(vec![du, dv, Vector::ZERO, Vector::ZERO, Vector::ZERO])
    }

    fn system(&self, x: &[f64], tol: Tolerances) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
        let p1 = self.first.point_at(x[0], x[1], tol).ok()?;
        let p2 = self.second.point_at(x[2], x[3], tol).ok()?;
        let (a1, b1) = self.first.d1_at(x[0], x[1], tol).ok()?;
        let (a2, b2) = self.second.d1_at(x[2], x[3], tol).ok()?;
        let (n1, dn1u, dn1v) = normal_and_derivatives(self.first, x[0], x[1], tol)?;
        let (n2, dn2u, dn2v) = normal_and_derivatives(self.second, x[2], x[3], tol)?;
        let (s1, s2) = (
            f64::from(self.sides.first) * self.radius,
            f64::from(self.sides.second) * self.radius,
        );

        // Three: the ball's centre is the same point from either side.
        let gap = (p1 + n1 * s1) - (p2 + n2 * s2);
        let c1u = a1 + dn1u * s1;
        let c1v = b1 + dn1v * s1;
        let c2u = a2 + dn2u * s2;
        let c2v = b2 + dn2v * s2;

        // One more: the section stands in the plane through the guide point
        // square to the guide. Stated with the *unnormalized* tangent, which
        // is the same plane and a simpler derivative.
        let w = x[4];
        let derivatives = self.guide.derivatives_at(w, 2, tol).ok()?;
        let g = Point::ORIGIN + derivatives[0];
        let (gd, gdd) = (derivatives[1], *derivatives.get(2).unwrap_or(&Vector::ZERO));
        let square = (p1 - g).dot(gd);

        Some((
            vec![gap.x, gap.y, gap.z, square],
            vec![
                vec![c1u.x, c1v.x, -c2u.x, -c2v.x, 0.0],
                vec![c1u.y, c1v.y, -c2u.y, -c2v.y, 0.0],
                vec![c1u.z, c1v.z, -c2u.z, -c2v.z, 0.0],
                vec![
                    a1.dot(gd),
                    b1.dot(gd),
                    0.0,
                    0.0,
                    (p1 - g).dot(gdd) - gd.dot(gd),
                ],
            ],
        ))
    }

    fn clamp(&self, x: &mut [f64]) {
        let (u1, v1) = clamp_to(self.first, x[0], x[1]);
        let (u2, v2) = clamp_to(self.second, x[2], x[3]);
        let (lo, hi) = self.guide.domain();
        x[0] = u1;
        x[1] = v1;
        x[2] = u2;
        x[3] = v2;
        // A guide that closes on itself has no end to stop at: the march
        // wraps its parameter and lets the *geometry* say when it is back
        // where it started.
        x[4] = if self.guide.is_periodic() && hi > lo {
            lo + (x[4] - lo).rem_euclid(hi - lo)
        } else {
            x[4].clamp(lo, hi)
        };
    }

    fn outside(&self, x: &[f64], tol: Tolerances) -> bool {
        let (lo, hi) = self.guide.domain();
        let band = tol.parametric();
        beyond(self.first, (x[0], x[1]), tol)
            || beyond(self.second, (x[2], x[3]), tol)
            || (!self.guide.is_periodic() && (x[4] < lo - band || x[4] > hi + band))
    }

    fn near_edge(&self, x: &[f64]) -> bool {
        let (lo, hi) = self.guide.domain();
        let reach = (hi - lo) * 1e-6;
        at_edge(self.first, (x[0], x[1]))
            || at_edge(self.second, (x[2], x[3]))
            || (!self.guide.is_periodic() && (x[4] <= lo + reach || x[4] >= hi - reach))
    }

    fn extent(&self) -> f64 {
        // A *length*, not a parameter span: the guide's parameter may be an
        // angle, and a step control fed an angle where it wanted a distance
        // walks a small circle in enormous steps and a large one in tiny ones.
        let (lo, hi) = self.guide.domain();
        let mut length = 0.0;
        let mut previous = None;
        for i in 0..=16 {
            let t = (hi - lo).mul_add(f64::from(i) / 16.0, lo);
            let Ok(p) = self.guide.point_at(t, Tolerances::millimetres()) else {
                continue;
            };
            if let Some(last) = previous {
                length += p.distance(last);
            }
            previous = Some(p);
        }
        length.max(self.radius * 8.0)
    }
}

/// A surface's unit normal.
fn unit_normal(surface: &SurfaceGeometry, u: f64, v: f64, tol: Tolerances) -> Option<Vector> {
    let (du, dv) = surface.d1_at(u, v, tol).ok()?;
    let cross = du.cross(dv);
    let length = cross.magnitude();
    if length <= tol.confusion() {
        return None;
    }
    Some(cross / length)
}

/// The unit normal and how it turns with each parameter.
///
/// Exactly, from the surface's own second derivatives: with `c = Sᵤ × Sᵥ` the
/// unnormalized normal, `∂n/∂u` is the part of `∂c/∂u` across `n`, over
/// `|c|` — the projection is what keeps a unit vector unit.
fn normal_and_derivatives(
    surface: &SurfaceGeometry,
    u: f64,
    v: f64,
    tol: Tolerances,
) -> Option<(Vector, Vector, Vector)> {
    let (su, sv) = surface.d1_at(u, v, tol).ok()?;
    let (suu, suv, svv) = surface.d2_at(u, v, tol).ok()?;
    let cross = su.cross(sv);
    let length = cross.magnitude();
    if length <= tol.confusion() {
        return None;
    }
    let n = cross / length;
    let dcu = suu.cross(sv) + su.cross(suv);
    let dcv = suv.cross(sv) + su.cross(svv);
    let across = |d: Vector| (d - n * d.dot(n)) / length;
    Some((n, across(dcu), across(dcv)))
}

/// Hold a parameter pair inside a surface's own domain.
fn clamp_to(surface: &SurfaceGeometry, u: f64, v: f64) -> (f64, f64) {
    let ((ua, ub), (va, vb)) = surface.domain();
    let hold = |value: f64, lo: f64, hi: f64, periodic: bool| {
        if periodic {
            let span = hi - lo;
            if span > 0.0 {
                return lo + (value - lo).rem_euclid(span);
            }
        }
        value.clamp(lo, hi)
    };
    (
        hold(u, ua, ub, surface.is_periodic_u()),
        hold(v, va, vb, surface.is_periodic_v()),
    )
}

/// Whether a parameter pair has left a surface's own domain.
fn beyond(surface: &SurfaceGeometry, at: (f64, f64), tol: Tolerances) -> bool {
    let ((ua, ub), (va, vb)) = surface.domain();
    let band = tol.parametric();
    (!surface.is_periodic_u() && (at.0 < ua - band || at.0 > ub + band))
        || (!surface.is_periodic_v() && (at.1 < va - band || at.1 > vb + band))
}

/// Whether it is at the edge of one — which is how a run-out is told from a
/// singularity.
fn at_edge(surface: &SurfaceGeometry, at: (f64, f64)) -> bool {
    let ((ua, ub), (va, vb)) = surface.domain();
    let near = |value: f64, lo: f64, hi: f64| {
        let reach = (hi - lo) * 1e-6;
        value <= lo + reach || value >= hi - reach
    };
    (!surface.is_periodic_u() && near(at.0, ua, ub))
        || (!surface.is_periodic_v() && near(at.1, va, vb))
}
