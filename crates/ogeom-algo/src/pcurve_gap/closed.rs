//! The gap between an analytic curve and an analytic pcurve's image, in
//! closed form.
//!
//! Over the fraction `f` of the way along both ranges alike, each side of
//! the pairs written here is `base + drift·f + cos·cos θ + sin·sin θ` with
//! `θ = phase + rate·f`: a line has no turning part, a circle or an ellipse
//! no drift, and an iso line of a plane, a cylinder, a cone, a sphere or a
//! torus is one of the two. The gap's widest is then read off the two
//! forms rather than sampled.

use ogeom_core::Tolerances;
use ogeom_geom::{Curve, Curve2d as _, Curve3d as _, PlanarCurve, Surface as _, SurfaceGeometry};
use ogeom_math::{Vector, Vector2};

/// How far a value computed from numbers no larger than a form's size may
/// stand off its exact value, per unit of that size. A point on either side
/// is formed, by the sampling and by the forms here alike, in at most a
/// few dozen rounded operations (the parameter's map and wrap, a sine and a
/// cosine, a frame applied), each off by at most half an ulp of the size,
/// so the measured gap and the bound read off the forms are each within a
/// few dozen ulps of the size of the true gap; 64 covers both.
const ROUNDING: f64 = 64.0 * f64::EPSILON;

/// A point moving with the fraction `f` of the way along a stretch.
#[derive(Debug, Clone, Copy)]
struct Form {
    base: Vector,
    drift: Vector,
    cos: Vector,
    sin: Vector,
    phase: f64,
    rate: f64,
    /// The largest magnitude among the numbers its point is formed from.
    size: f64,
}

impl Form {
    fn affine(base: Vector, drift: Vector, size: f64) -> Self {
        Self {
            base,
            drift,
            cos: Vector::ZERO,
            sin: Vector::ZERO,
            phase: 0.0,
            rate: 0.0,
            size,
        }
    }

    fn turning(centre: Vector, (cos, sin): (Vector, Vector), (start, end): (f64, f64)) -> Self {
        let turn = cos.magnitude() + sin.magnitude();
        Self {
            base: centre,
            drift: Vector::ZERO,
            cos,
            sin,
            phase: start,
            rate: end - start,
            size: centre.magnitude() + turn * (1.0 + start.abs().max(end.abs())),
        }
    }

    fn turns(&self) -> bool {
        self.cos != Vector::ZERO || self.sin != Vector::ZERO
    }
}

/// A planar curve's point over the fraction along its range: the same
/// shape of form, in the plane.
#[derive(Debug, Clone, Copy)]
struct Form2 {
    base: Vector2,
    drift: Vector2,
    cos: Vector2,
    sin: Vector2,
    phase: f64,
    rate: f64,
    /// The largest magnitude among the numbers its point is formed from.
    size: f64,
}

/// `(start, end)` signed by a reversal: the angles a reversed circle or
/// ellipse is evaluated at.
fn angles((start, end): (f64, f64), reversed: bool) -> (f64, f64) {
    if reversed {
        (-start, -end)
    } else {
        (start, end)
    }
}

/// Whether `range` lies within `domain`, so evaluating over it clamps
/// nothing.
fn within(range: (f64, f64), domain: (f64, f64)) -> bool {
    range.0.min(range.1) >= domain.0 && range.0.max(range.1) <= domain.1
}

fn curve_form(curve: &Curve, range: (f64, f64)) -> Option<Form> {
    match curve {
        Curve::Line(line) => {
            if !within(range, line.domain()) {
                return None;
            }
            let axis = line.axis();
            let (origin, along) = (axis.location.to_vector(), axis.direction.vector());
            Some(Form::affine(
                origin + along * range.0,
                along * (range.1 - range.0),
                origin.magnitude() + range.0.abs().max(range.1.abs()),
            ))
        }
        Curve::Circle(circle) => {
            let (c, frame) = (circle.circle(), circle.circle().frame());
            let r = c.radius();
            Some(Form::turning(
                c.centre().to_vector(),
                (frame.x().vector() * r, frame.y().vector() * r),
                angles(range, circle.is_reversed()),
            ))
        }
        Curve::Ellipse(ellipse) => {
            let (e, frame) = (ellipse.ellipse(), ellipse.ellipse().frame());
            Some(Form::turning(
                e.centre().to_vector(),
                (
                    frame.x().vector() * e.major_radius(),
                    frame.y().vector() * e.minor_radius(),
                ),
                angles(range, ellipse.is_reversed()),
            ))
        }
        Curve::Trimmed(trimmed) => {
            // Within its domain, a trimmed curve's parameter is its basis's,
            // or mirrored across the domain, either way affine in its own.
            let (a, b) = trimmed.domain();
            if !within(range, (a, b)) {
                return None;
            }
            let basis = if trimmed.is_reversed() {
                (a + b - range.0, a + b - range.1)
            } else {
                range
            };
            let form = curve_form(trimmed.basis(), basis)?;
            let speed = 1.0 + form.cos.magnitude() + form.sin.magnitude();
            Some(Form {
                size: form.size + (a.abs() + b.abs()) * speed,
                ..form
            })
        }
        _ => None,
    }
}

fn planar_form(pcurve: &PlanarCurve, range: (f64, f64)) -> Option<Form2> {
    let turning = |centre: Vector2, cos: Vector2, sin: Vector2, (phase, end): (f64, f64)| Form2 {
        base: centre,
        drift: Vector2::ZERO,
        cos,
        sin,
        phase,
        rate: end - phase,
        size: centre.magnitude()
            + (cos.magnitude() + sin.magnitude()) * (1.0 + phase.abs().max(end.abs())),
    };
    match pcurve {
        PlanarCurve::Line(line) => {
            if !within(range, line.domain()) {
                return None;
            }
            let axis = line.axis();
            let (origin, along) = (axis.location.to_vector(), axis.direction.vector());
            Some(Form2 {
                base: origin + along * range.0,
                drift: along * (range.1 - range.0),
                cos: Vector2::ZERO,
                sin: Vector2::ZERO,
                phase: 0.0,
                rate: 0.0,
                size: origin.magnitude() + range.0.abs().max(range.1.abs()),
            })
        }
        PlanarCurve::Circle(circle) => {
            let (c, frame) = (circle.circle(), circle.circle().frame());
            let r = c.radius();
            Some(turning(
                c.centre().to_vector(),
                frame.x().vector() * r,
                frame.y().vector() * r,
                angles(range, circle.is_reversed()),
            ))
        }
        PlanarCurve::Ellipse(ellipse) => {
            let (e, frame) = (ellipse.ellipse(), ellipse.ellipse().frame());
            Some(turning(
                e.centre().to_vector(),
                frame.x().vector() * e.major_radius(),
                frame.y().vector() * e.minor_radius(),
                angles(range, ellipse.is_reversed()),
            ))
        }
        _ => None,
    }
}

/// A planar form lifted through `surface`: any form on a plane, whose
/// chart is affine; on a surface of revolution only a line along one chart
/// direction, an iso line, which is a ruling or a circle. `None` where a
/// lifted point would be clamped to the surface's domain or refused.
fn lifted_form(pcurve: &PlanarCurve, range: (f64, f64), surface: &SurfaceGeometry) -> Option<Form> {
    let p = planar_form(pcurve, range)?;
    let lift = |w: Vector2, x: Vector, y: Vector| x * w.x + y * w.y;
    if let SurfaceGeometry::Plane(plane) = surface {
        let frame = plane.plane().frame();
        let (x, y) = (frame.x().vector(), frame.y().vector());
        let origin = frame.origin().to_vector();
        return Some(Form {
            base: origin + lift(p.base, x, y),
            drift: lift(p.drift, x, y),
            cos: lift(p.cos, x, y),
            sin: lift(p.sin, x, y),
            phase: p.phase,
            rate: p.rate,
            size: origin.magnitude() + p.size,
        });
    }
    if p.cos != Vector2::ZERO || p.sin != Vector2::ZERO {
        return None;
    }
    let (start, end) = (p.base, p.base + p.drift);
    let ((ua, ub), (va, vb)) = surface.domain();
    // The iso line's fixed chart coordinate, and the moving one's ends.
    let along_u = p.drift.y == 0.0;
    if !along_u && p.drift.x != 0.0 {
        return None;
    }
    let periodic = (surface.is_periodic_u(), surface.is_periodic_v());
    let inside = |(a, b): (f64, f64), periodic: bool, (lo, hi): (f64, f64)| {
        periodic || within((a, b), (lo, hi))
    };
    if !inside((start.x, end.x), periodic.0, (ua, ub))
        || !inside((start.y, end.y), periodic.1, (va, vb))
    {
        return None;
    }
    let reach = start.x.abs().max(end.x.abs()) + start.y.abs().max(end.y.abs());
    // Along `u` a parallel, a circle about the axis at the line's `v`;
    // along `v` a ruling or a meridian at the line's `u`.
    match surface {
        SurfaceGeometry::Cylinder(cylinder) => {
            let c = cylinder.cylinder();
            let f = c.frame();
            let (x, y, z) = (f.x().vector(), f.y().vector(), f.z().vector());
            let (origin, r) = (f.origin().to_vector(), c.radius());
            let size = origin.magnitude() + r * (1.0 + reach) + reach + p.size;
            Some(if along_u {
                Form {
                    size,
                    ..Form::turning(origin + z * start.y, (x * r, y * r), (start.x, end.x))
                }
            } else {
                let (sin, cos) = start.x.sin_cos();
                let out = x * cos + y * sin;
                Form::affine(origin + out * r + z * start.y, z * (end.y - start.y), size)
            })
        }
        SurfaceGeometry::Cone(cone) => {
            let c = cone.cone();
            let f = c.frame();
            let (x, y, z) = (f.x().vector(), f.y().vector(), f.z().vector());
            let origin = f.origin().to_vector();
            let slope = c.half_angle().tan();
            let size = origin.magnitude()
                + (c.reference_radius().abs() + slope.abs() * reach) * (1.0 + reach)
                + reach
                + p.size;
            Some(if along_u {
                let r = c.radius_at(start.y);
                Form {
                    size,
                    ..Form::turning(origin + z * start.y, (x * r, y * r), (start.x, end.x))
                }
            } else {
                let (sin, cos) = start.x.sin_cos();
                let out = x * cos + y * sin;
                Form::affine(
                    origin + out * c.radius_at(start.y) + z * start.y,
                    (out * slope + z) * (end.y - start.y),
                    size,
                )
            })
        }
        SurfaceGeometry::Sphere(sphere) => {
            let s = sphere.sphere();
            let f = s.frame();
            let (x, y, z) = (f.x().vector(), f.y().vector(), f.z().vector());
            let (centre, r) = (s.centre().to_vector(), s.radius());
            let size = centre.magnitude() + r * (1.0 + reach) + p.size;
            Some(Form {
                size,
                ..if along_u {
                    let (sin, cos) = start.y.sin_cos();
                    let ring = r * cos;
                    Form::turning(
                        centre + z * (r * sin),
                        (x * ring, y * ring),
                        (start.x, end.x),
                    )
                } else {
                    let (sin, cos) = start.x.sin_cos();
                    let out = x * cos + y * sin;
                    Form::turning(centre, (out * r, z * r), (start.y, end.y))
                }
            })
        }
        SurfaceGeometry::Torus(torus) => {
            let t = torus.torus();
            let f = t.frame();
            let (x, y, z) = (f.x().vector(), f.y().vector(), f.z().vector());
            let origin = f.origin().to_vector();
            let (major, minor) = (t.major_radius(), t.minor_radius());
            let size = origin.magnitude() + (major.abs() + minor.abs()) * (1.0 + reach) + p.size;
            Some(Form {
                size,
                ..if along_u {
                    let (sin, cos) = start.y.sin_cos();
                    let radius = minor.mul_add(cos, major);
                    Form::turning(
                        origin + z * (minor * sin),
                        (x * radius, y * radius),
                        (start.x, end.x),
                    )
                } else {
                    let (sin, cos) = start.x.sin_cos();
                    let out = x * cos + y * sin;
                    Form::turning(
                        origin + out * major,
                        (out * minor, z * minor),
                        (start.y, end.y),
                    )
                }
            })
        }
        _ => None,
    }
}

/// The widest of `|u·cos θ + v·sin θ|` over every `θ`: the larger singular
/// value of the pair, the root of their Gram matrix's larger eigenvalue.
fn widest_turn(u: Vector, v: Vector) -> f64 {
    let (p, q, w) = (u.dot(u), v.dot(v), u.dot(v));
    let half = (p - q) * 0.5;
    ((p + q) * 0.5 + half.hypot(w)).max(0.0).sqrt()
}

/// A bound on the widest of `|a(f) - b(f)|` over `[0, 1]`, for two forms
/// that both turn.
///
/// Let `κ` be the sign of the two rates' product and `δ = b.phase -
/// κ·a.phase`, so that `θb = κ·θa + δ + ε(f)` with `ε(f) = (b.rate -
/// κ·a.rate)·f`. Then `b`'s turning part at `κ·θa + δ` is
/// `cos θa·(cos δ·b.cos + sin δ·b.sin) + sin θa·κ·(cos δ·b.sin - sin δ·b.cos)`,
/// and moving its angle by `ε` moves it by at most `|ε|` times the root of
/// `|b.cos|² + |b.sin|²`, its speed's largest. The difference is the bases'
/// and drifts' affine part, widest at an end, plus a turning part in `θa`,
/// no wider than [`widest_turn`], plus that move.
fn turning_bound(a: &Form, b: &Form) -> f64 {
    let kappa: f64 = if a.rate * b.rate < 0.0 { -1.0 } else { 1.0 };
    let delta = kappa.mul_add(-a.phase, b.phase);
    let (sin, cos) = delta.sin_cos();
    let along_cos = b.cos * cos + b.sin * sin;
    let along_sin = (b.sin * cos - b.cos * sin) * kappa;
    let shift = a.base - b.base;
    let affine = shift
        .magnitude()
        .max((shift + a.drift - b.drift).magnitude());
    let slip = kappa.mul_add(-a.rate, b.rate).abs();
    let speed = b.cos.magnitude().hypot(b.sin.magnitude());
    affine + widest_turn(a.cos - along_cos, a.sin - along_sin) + slip * speed
}

/// The widest gap between `curve` and `pcurve` lifted through `surface`,
/// read off their forms, with how far it may stand off what sampling the
/// stretch measures; `None` where either has no form here, where the two
/// are not both straight or both turning, where the bound is not within
/// rounding of the gap measured at the ends, or where it could exceed
/// `stated`.
///
/// Within `stated` the sampled gap compares every sample with the same
/// parameter's point, so its answer is the widest of the direct gaps over
/// the samples: no wider than the true widest, which the bound covers, and
/// no narrower than the widest of the ends. Where both are straight the
/// gap is affine in `f` and its length widest at an end, so the ends
/// measured as sampling measures them are the answer.
pub(super) fn closed_gap(
    (curve, range): (&Curve, (f64, f64)),
    (pcurve, pcurve_range): (&PlanarCurve, (f64, f64)),
    surface: &SurfaceGeometry,
    stated: f64,
    tol: Tolerances,
) -> Option<(f64, f64)> {
    let finite = |(a, b): (f64, f64)| a.is_finite() && b.is_finite();
    if !finite(range) || !finite(pcurve_range) {
        return None;
    }
    let a = curve_form(curve, range)?;
    let b = lifted_form(pcurve, pcurve_range, surface)?;
    let end = |f: f64| -> Option<f64> {
        let uv = pcurve
            .point_at(pcurve_range.0 + (pcurve_range.1 - pcurve_range.0) * f, tol)
            .ok()?;
        let lifted = surface.point_at(uv.x, uv.y, tol).ok()?;
        let at = curve
            .point_at(range.0 + (range.1 - range.0) * f, tol)
            .ok()?;
        Some(at.distance(lifted))
    };
    let ends = end(0.0)?.max(end(1.0)?);
    let margin = ROUNDING * (a.size + b.size);
    let gap = match (a.turns(), b.turns()) {
        (false, false) => ends,
        (true, true) => {
            let bound = turning_bound(&a, &b);
            if bound > ends + margin {
                return None;
            }
            bound
        }
        _ => return None,
    };
    (gap.is_finite() && gap + margin <= stated).then_some((gap, margin))
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::primitive::{make_box, make_cone, make_cylinder, make_sphere, make_torus};
    use ogeom_geom::{Circle2d, CircleCurve, Line2d, LineCurve, PlaneSurface, Reversible as _};
    use ogeom_math::{Circle, Circle2, Direction, Frame, Frame2, Plane, Point, Point2};
    use ogeom_topo::{Model, ShapeType};

    const T: Tolerances = Tolerances::millimetres();
    const STATED: f64 = 1e-7;

    /// The closed form's answer and the sampled one, for a pcurve on the
    /// world plane.
    fn on_plane(
        curve: (&Curve, (f64, f64)),
        pcurve: (&PlanarCurve, (f64, f64)),
    ) -> (Option<(f64, f64)>, f64) {
        let plane: SurfaceGeometry = PlaneSurface::new(Plane::new(Frame::WORLD)).into();
        let closed = closed_gap(curve, pcurve, &plane, STATED, T);
        let sampled = super::super::sampled_gap(curve, pcurve, &plane, STATED, T).unwrap();
        (closed, sampled)
    }

    fn circle(r: f64) -> Curve {
        CircleCurve::new(Circle::new(Frame::WORLD, r, T).unwrap()).into()
    }

    fn circle2(centre: (f64, f64), r: f64) -> PlanarCurve {
        let frame = Frame2::new(Point2::new(centre.0, centre.1), ogeom_math::Direction2::X);
        Circle2d::new(Circle2::new(frame, r, T).unwrap()).into()
    }

    /// Every pcurve of a box, a cylinder, a cone, a sphere and a torus,
    /// placed off the origin and tilted, is a line or a circle against a
    /// line or an iso line, and each is answered in closed form within its
    /// margin of what sampling measures.
    #[test]
    fn every_pcurve_of_the_primitives_is_answered_within_rounding_of_sampling() {
        let mut model = Model::new();
        let tilt = Direction::new(Vector::new(0.0, 0.6, 0.8), T).unwrap();
        let frame = Frame::new(Point::new(12.5, -7.25, 3.0), tilt, Direction::X, T).unwrap();
        let solids = [
            make_box(&mut model, frame, (3.0, 4.0, 5.0), T)
                .unwrap()
                .shape,
            make_cylinder(&mut model, frame, 2.0, 5.0, T).unwrap().shape,
            make_cone(&mut model, frame, 3.0, 1.0, 4.0, T)
                .unwrap()
                .shape,
            make_sphere(&mut model, frame, 2.5, T).unwrap().shape,
            make_torus(&mut model, frame, 4.0, 1.5, T).unwrap().shape,
        ];
        let mut answered = 0;
        for solid in &solids {
            for edge in ogeom_topo::explore_unique(&model, solid, ShapeType::Edge).unwrap() {
                super::super::edge_gap_by(
                    &model,
                    &edge,
                    T,
                    |curve, pcurve, surface, stated, tol| {
                        let sampled =
                            super::super::sampled_gap(curve, pcurve, surface, stated, tol)?;
                        let closed = closed_gap(curve, pcurve, surface, stated, tol);
                        let (gap, margin) =
                            closed.expect("an analytic pair answered in closed form");
                        assert!(
                            (gap - sampled).abs() <= margin,
                            "{gap:e} against {sampled:e}"
                        );
                        answered += 1;
                        Ok(sampled)
                    },
                )
                .unwrap();
            }
        }
        // The box's 12 edges on two faces each, the cylinder's and the
        // cone's two rims on two faces and seam, the sphere's seam and the
        // torus's two seams.
        assert_eq!(answered, 24 + 6 + 6 + 2 + 4);
    }

    /// A circle's centre moved in its plane: the gap is the move all
    /// along, and the closed form says so.
    #[test]
    fn a_circle_moved_in_its_plane_is_answered_with_the_move() {
        let range = (0.3, 2.9);
        let (closed, sampled) =
            on_plane((&circle(2.0), range), (&circle2((3e-9, 0.0), 2.0), range));
        let (gap, margin) = closed.unwrap();
        assert!((gap - 3e-9).abs() <= 1e-15, "{gap:e}");
        assert!(
            (gap - sampled).abs() <= margin,
            "{gap:e} against {sampled:e}"
        );
    }

    /// Two lines leaving one point apart: the gap grows along them, and is
    /// the far ends' distance.
    #[test]
    fn two_lines_leaving_one_point_are_answered_with_their_far_ends_distance() {
        let curve: Curve = LineCurve::segment(Point::ORIGIN, Point::new(10.0, 0.0, 0.0), T)
            .unwrap()
            .into();
        let pcurve: PlanarCurve = Line2d::segment(Point2::ORIGIN, Point2::new(10.0, 2e-8), T)
            .unwrap()
            .into();
        let (closed, sampled) = on_plane((&curve, (0.0, 10.0)), (&pcurve, (0.0, 10.0)));
        let (gap, _) = closed.unwrap();
        assert_eq!(gap, sampled);
        assert!((gap - 2e-8).abs() <= 1e-15, "{gap:e}");
    }

    /// A circle about `-z` against a reversed circle on the plane: the two
    /// angles run opposite ways, and the same parameter's points coincide,
    /// so the gap is answered in closed form as nothing.
    #[test]
    fn a_circle_against_an_image_turning_the_other_way_at_the_same_points_is_answered() {
        let below = Frame::new(Point::ORIGIN, -Direction::Z, Direction::X, T).unwrap();
        let curve: Curve = CircleCurve::new(Circle::new(below, 2.0, T).unwrap()).into();
        let range = (0.4, 2.2);
        let back = circle2((0.0, 0.0), 2.0).reversed();
        let (closed, sampled) = on_plane((&curve, range), (&back, range));
        let (gap, margin) = closed.unwrap();
        assert!(gap <= margin, "{gap:e}");
        assert!(sampled <= margin, "{sampled:e}");
    }

    /// A circle's image run the other way round it stands on it, though
    /// the same parameter's points are apart: sampling measures against the
    /// nearest point of the curve, and the closed form, measuring the
    /// direct gap, leaves it to sampling.
    #[test]
    fn a_circle_run_the_other_way_round_is_left_to_sampling() {
        let range = (0.0, core::f64::consts::TAU);
        let back = circle2((0.0, 0.0), 2.0).reversed();
        let (closed, sampled) = on_plane((&circle(2.0), range), (&back, range));
        assert!(closed.is_none());
        assert!(sampled <= STATED, "{sampled:e}");
    }

    /// A gap past the stated tolerance is left to sampling, which measures
    /// it against the nearest point of the curve.
    #[test]
    fn a_gap_past_the_stated_tolerance_is_left_to_sampling() {
        let range = (0.3, 2.9);
        let (closed, sampled) =
            on_plane((&circle(2.0), range), (&circle2((1e-6, 0.0), 2.0), range));
        assert!(closed.is_none());
        assert!(sampled > STATED, "{sampled:e}");
    }

    /// A circle moved across an arc and widened: the bound adds the move
    /// and the widening, which line up only off the arc, so it stands wider
    /// than the ends measure and is left to sampling.
    #[test]
    fn a_bound_wider_than_the_ends_is_left_to_sampling() {
        let half = core::f64::consts::FRAC_PI_2;
        let range = (half - 0.5, half + 0.5);
        let pcurve = circle2((1e-8, 0.0), 2.0 + 1e-8);
        let (closed, sampled) = on_plane((&circle(2.0), range), (&pcurve, range));
        assert!(closed.is_none());
        assert!(sampled < 1.8e-8, "{sampled:e}");
    }
}
