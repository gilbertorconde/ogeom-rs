//! Where a curve pierces a surface.
//!
//! *Elsewhere* this is `GeomAPI_IntCS` and the line/quadric half of `IntAna`.
//! Two consumers drive it: edge/face interference in the boolean's pave
//! filler, and the exact point-in-solid classifier, which is a ray/surface
//! query per face — the very use `docs/SCOPE.md` schedules "exact after §7".
//!
//! # Well-posed, for once
//!
//! `C(t) = S(u, v)` is three equations in three unknowns — unlike the
//! surface/surface system, nothing has to be pinned for Newton to converge to
//! a point. The analytic cases are still answered in closed form first: a line
//! against a plane or a quadric is a linear or quadratic equation, and solving
//! a quadratic by iteration would be slower and less exact than writing down
//! its roots.
//!
//! # A curve lying in the surface
//!
//! A line in a plane crosses it nowhere and everywhere. That is an overlap —
//! the parameter range of the curve that lies in the surface — and it is a
//! different answer from any list of points. Detected where the analytic
//! forms can see it; the general path reports whatever isolated piercings its
//! sampling resolves, and says so.

use og_core::{OgResult, Tolerances, og_bail};
use og_geom::{Curve, Curve3d, Surface, SurfaceGeometry};
use og_math::{Point, solve};

use crate::march::{Cell, sample, segment_meets_triangle};

/// One piercing of a surface by a curve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Piercing {
    /// The parameter on the curve.
    pub on_curve: f64,
    /// The parameters on the surface.
    pub on_surface: (f64, f64),
    /// Where, taken from the curve.
    pub point: Point,
    /// The distance between the two evaluations there.
    pub gap: f64,
}

/// What a curve does to a surface.
#[derive(Debug, Clone, PartialEq)]
pub struct CurveSurfaceIntersection {
    /// Isolated piercings, in order along the curve.
    pub crossings: Vec<Piercing>,
    /// Parameter ranges of the curve that lie *in* the surface.
    ///
    /// Detected for the analytic cases — a line in a plane. The general path
    /// cannot see lying-on and reports whatever isolated piercings its
    /// sampling resolves.
    pub lying: Vec<(f64, f64)>,
}

impl CurveSurfaceIntersection {
    /// No contact found.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.crossings.is_empty() && self.lying.is_empty()
    }

    const fn empty() -> Self {
        Self {
            crossings: Vec::new(),
            lying: Vec::new(),
        }
    }
}

/// How hard the general path looks.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveSurfaceOptions {
    /// How many segments the curve is sampled into for seeding.
    pub samples: usize,
    /// How finely the surface is sampled, per direction.
    pub grid: usize,
    /// The widest gap that still counts as a piercing.
    pub gap: f64,
}

impl Default for CurveSurfaceOptions {
    fn default() -> Self {
        Self {
            samples: 128,
            grid: 24,
            gap: 1e-7,
        }
    }
}

/// Where a curve pierces a surface.
///
/// Analytic line/plane, line/sphere and line/cylinder are answered in closed
/// form; everything else is seeded polyhedrally and polished by Newton on the
/// well-posed three-by-three system.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the options
/// are unusable.
pub fn intersect_curve_surface(
    curve: &Curve,
    surface: &SurfaceGeometry,
    options: CurveSurfaceOptions,
    tol: Tolerances,
) -> OgResult<CurveSurfaceIntersection> {
    if options.samples < 2 || options.grid < 2 {
        og_bail!(Construction, "seeding needs at least two steps each way");
    }
    if !options.gap.is_finite() || options.gap <= 0.0 {
        og_bail!(Construction, "a gap of {} is not a distance", options.gap);
    }

    match (curve, surface) {
        (Curve::Line(line), SurfaceGeometry::Plane(p)) => {
            Ok(line_plane(line, p.plane(), curve, surface, tol))
        }
        (Curve::Line(line), SurfaceGeometry::Sphere(s)) => Ok(line_quadric(
            line,
            curve,
            surface,
            sphere_roots(line, s.sphere()),
            options,
            tol,
        )),
        (Curve::Line(line), SurfaceGeometry::Cylinder(c)) => Ok(line_quadric(
            line,
            curve,
            surface,
            cylinder_roots(line, c.cylinder()),
            options,
            tol,
        )),
        _ => general(curve, surface, options, tol),
    }
}

// --- analytic ----------------------------------------------------------------

fn line_plane(
    line: &og_geom::LineCurve,
    plane: og_math::Plane,
    curve: &Curve,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> CurveSurfaceIntersection {
    let axis = line.axis();
    let along = plane.normal().dot(axis.direction);
    let height = plane.signed_distance_to(axis.location);

    if along.abs() <= tol.angular() {
        // Parallel: in the plane, or never touching it.
        if height.abs() <= tol.confusion() {
            return CurveSurfaceIntersection {
                crossings: Vec::new(),
                lying: vec![line.domain()],
            };
        }
        return CurveSurfaceIntersection::empty();
    }

    let t = -height / along;
    let (lo, hi) = line.domain();
    if t < lo - tol.parametric() || t > hi + tol.parametric() {
        return CurveSurfaceIntersection::empty();
    }
    let point = axis.location + axis.direction.vector() * t;
    let Some(found) = invert(surface, point, curve, t, tol) else {
        return CurveSurfaceIntersection::empty();
    };
    if found.gap > tol.confusion() {
        // The crossing is real on the unbounded plane but outside this
        // surface's stated extents; the clamped polish says so as a gap.
        return CurveSurfaceIntersection::empty();
    }
    CurveSurfaceIntersection {
        crossings: vec![found],
        lying: Vec::new(),
    }
}

/// The line parameters at which a line meets a sphere.
fn sphere_roots(line: &og_geom::LineCurve, sphere: og_math::Sphere) -> Vec<f64> {
    let axis = line.axis();
    let d = axis.direction.vector();
    let m = axis.location - sphere.centre();
    // |m + t d|^2 = r^2, with |d| = 1.
    let b = m.dot(d);
    let c = sphere.radius().mul_add(-sphere.radius(), m.dot(m));
    let discriminant = b.mul_add(b, -c);
    if discriminant < 0.0 {
        return Vec::new();
    }
    let root = discriminant.sqrt();
    if root == 0.0 {
        vec![-b]
    } else {
        vec![-b - root, -b + root]
    }
}

/// The line parameters at which a line meets a cylinder.
fn cylinder_roots(line: &og_geom::LineCurve, cylinder: og_math::Cylinder) -> Vec<f64> {
    let axis = line.axis();
    let w = cylinder.axis().direction.vector();
    // Strip the components along the cylinder's axis; what is left is a 2D
    // circle problem in the perpendicular plane.
    let d = axis.direction.vector();
    let m = axis.location - cylinder.axis().location;
    let d_perp = d - w * d.dot(w);
    let m_perp = m - w * m.dot(w);
    let a = d_perp.dot(d_perp);
    if a <= f64::MIN_POSITIVE {
        // The line runs along the axis direction: on the wall it would lie,
        // not pierce, and lying is not detected here.
        return Vec::new();
    }
    let b = d_perp.dot(m_perp);
    let c = cylinder
        .radius()
        .mul_add(-cylinder.radius(), m_perp.dot(m_perp));
    let discriminant = b.mul_add(b, -(a * c));
    if discriminant < 0.0 {
        return Vec::new();
    }
    let root = discriminant.sqrt();
    if root == 0.0 {
        vec![-b / a]
    } else {
        vec![(-b - root) / a, (-b + root) / a]
    }
}

/// Roots dressed as piercings, filtered by the line's own range and the
/// surface's extents.
fn line_quadric(
    line: &og_geom::LineCurve,
    curve: &Curve,
    surface: &SurfaceGeometry,
    roots: Vec<f64>,
    options: CurveSurfaceOptions,
    tol: Tolerances,
) -> CurveSurfaceIntersection {
    let axis = line.axis();
    let (lo, hi) = line.domain();
    let mut crossings = Vec::new();
    for t in roots {
        if t < lo - tol.parametric() || t > hi + tol.parametric() {
            continue;
        }
        let point = axis.location + axis.direction.vector() * t;
        let Some(found) = invert(surface, point, curve, t, tol) else {
            continue;
        };
        // The extent check is the gap. The polish clamps the surface
        // parameters into the stated domain, so a root beyond the cylinder's
        // height converges to the rim with a gap of exactly how far past it
        // was — a piercing of the unbounded geometry, not of this surface.
        // Discarding the polish's gap and writing zero here was the bug this
        // comment replaces.
        if found.gap > tol.confusion() {
            continue;
        }
        let _ = options;
        crossings.push(Piercing {
            on_curve: t,
            on_surface: found.on_surface,
            point,
            gap: found.gap,
        });
    }
    crossings.sort_by(|a, b| {
        a.on_curve
            .partial_cmp(&b.on_curve)
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    CurveSurfaceIntersection {
        crossings,
        lying: Vec::new(),
    }
}

/// Surface parameters of a point known to lie on an analytic surface.
///
/// Closed-form inversion for the quadrics; refined by one Newton pass so the
/// reported parameters evaluate back onto the point to rounding.
fn invert(
    surface: &SurfaceGeometry,
    point: Point,
    curve: &Curve,
    on_curve: f64,
    tol: Tolerances,
) -> Option<Piercing> {
    let guess = match surface {
        SurfaceGeometry::Plane(p) => {
            let local = p.plane().frame().to_local(point);
            (local.x, local.y)
        }
        SurfaceGeometry::Sphere(s) => {
            let local = s.sphere().frame().to_local(point);
            let latitude = (local.z / s.sphere().radius()).clamp(-1.0, 1.0).asin();
            (
                local.y.atan2(local.x).rem_euclid(core::f64::consts::TAU),
                latitude,
            )
        }
        SurfaceGeometry::Cylinder(c) => {
            let local = c.cylinder().frame().to_local(point);
            (
                local.y.atan2(local.x).rem_euclid(core::f64::consts::TAU),
                local.z,
            )
        }
        _ => return None,
    };
    // One polish step against the curve point, so parameter rounding in the
    // inversion does not survive into the result — and the gap comes with it,
    // because the polish clamps into the surface's extents and the gap is
    // what says whether the clamped answer still touches the curve.
    polish(curve, surface, on_curve, guess, tol)
}

// --- general -----------------------------------------------------------------

fn general(
    curve: &Curve,
    surface: &SurfaceGeometry,
    options: CurveSurfaceOptions,
    tol: Tolerances,
) -> OgResult<CurveSurfaceIntersection> {
    let cells = sample(surface, options.grid, tol);
    let (lo, hi) = curve.domain();

    let mut points = Vec::with_capacity(options.samples + 1);
    for i in 0..=options.samples {
        #[allow(clippy::cast_precision_loss)]
        let t = lo + (hi - lo) * i as f64 / options.samples as f64;
        if let Ok(p) = curve.point_at(t, tol) {
            points.push((t, p));
        }
    }

    let mut crossings: Vec<Piercing> = Vec::new();
    for pair in points.windows(2) {
        let (t0, p0) = pair[0];
        let (t1, p1) = pair[1];
        for cell in &cells {
            if !segment_near_cell(p0, p1, cell, options.gap) {
                continue;
            }
            if segment_meets_triangle(p0, p1, cell.corners).is_none() {
                continue;
            }
            let seed_t = f64::midpoint(t0, t1);
            if let Some(found) = polish(curve, surface, seed_t, cell.at, tol) {
                if found.gap > options.gap {
                    continue;
                }
                let reach = tol.confusion() * 100.0;
                if !crossings
                    .iter()
                    .any(|c| c.point.distance(found.point) <= reach)
                {
                    crossings.push(found);
                }
            }
        }
    }
    crossings.sort_by(|a, b| {
        a.on_curve
            .partial_cmp(&b.on_curve)
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    Ok(CurveSurfaceIntersection {
        crossings,
        lying: Vec::new(),
    })
}

/// Whether a segment's box comes near a cell's.
fn segment_near_cell(a: Point, b: Point, cell: &Cell, margin: f64) -> bool {
    let low = Point::new(a.x.min(b.x), a.y.min(b.y), a.z.min(b.z));
    let high = Point::new(a.x.max(b.x), a.y.max(b.y), a.z.max(b.z));
    low.x <= cell.high.x + margin
        && cell.low.x <= high.x + margin
        && low.y <= cell.high.y + margin
        && cell.low.y <= high.y + margin
        && low.z <= cell.high.z + margin
        && cell.low.z <= high.z + margin
}

/// Newton on the well-posed system `C(t) = S(u, v)`.
fn polish(
    curve: &Curve,
    surface: &SurfaceGeometry,
    seed_t: f64,
    seed_uv: (f64, f64),
    tol: Tolerances,
) -> Option<Piercing> {
    let clamp_t = |t: f64| {
        let (lo, hi) = curve.domain();
        if curve.is_periodic() {
            let span = hi - lo;
            if span > 0.0 {
                return lo + (t - lo).rem_euclid(span);
            }
        }
        t.clamp(lo, hi)
    };
    let clamp_uv = |u: f64, v: f64| {
        let ((ua, ub), (va, vb)) = surface.domain();
        let fold = |x: f64, lo: f64, hi: f64, periodic: bool| {
            if periodic {
                let span = hi - lo;
                if span > 0.0 {
                    return lo + (x - lo).rem_euclid(span);
                }
            }
            x.clamp(lo, hi)
        };
        (
            fold(u, ua, ub, surface.is_periodic_u()),
            fold(v, va, vb, surface.is_periodic_v()),
        )
    };

    let system = |x: &[f64]| {
        let t = clamp_t(x[0]);
        let (u, v) = clamp_uv(x[1], x[2]);
        let pc = curve.point_at(t, tol).unwrap_or(Point::ORIGIN);
        let ps = surface.point_at(u, v, tol).unwrap_or(Point::ORIGIN);
        let dc = curve.d1_at(t, tol).unwrap_or(og_math::Vector::ZERO);
        let (du, dv) = surface
            .d1_at(u, v, tol)
            .unwrap_or((og_math::Vector::ZERO, og_math::Vector::ZERO));
        let gap = pc - ps;
        (
            vec![gap.x, gap.y, gap.z],
            vec![
                vec![dc.x, -du.x, -dv.x],
                vec![dc.y, -du.y, -dv.y],
                vec![dc.z, -du.z, -dv.z],
            ],
        )
    };
    let criteria = solve::Criteria {
        residual: tol.confusion() * 0.01,
        step: tol.parametric(),
        max_iterations: 40,
    };
    let found = solve::newton_system(system, &[seed_t, seed_uv.0, seed_uv.1], criteria).ok()?;
    let t = clamp_t(found.value[0]);
    let (u, v) = clamp_uv(found.value[1], found.value[2]);
    let pc = curve.point_at(t, tol).ok()?;
    let ps = surface.point_at(u, v, tol).ok()?;
    Some(Piercing {
        on_curve: t,
        on_surface: (u, v),
        point: pc,
        gap: pc.distance(ps),
    })
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use og_geom::{
        BSplineCurve, CircleCurve, CylinderSurface, LineCurve, PlaneSurface, SphereSurface,
    };
    use og_math::{Circle, Cylinder, Direction, Frame, KnotVector, Plane, Sphere, Vector};

    const T: Tolerances = Tolerances::millimetres();

    fn sphere(radius: f64) -> SurfaceGeometry {
        SphereSurface::new(Sphere::centred(Point::ORIGIN, radius, T).unwrap()).into()
    }

    fn cylinder(radius: f64, height: (f64, f64)) -> SurfaceGeometry {
        CylinderSurface::new(Cylinder::new(Frame::WORLD, radius, T).unwrap(), height)
            .unwrap()
            .into()
    }

    fn plane(origin: Point, normal: Vector) -> SurfaceGeometry {
        PlaneSurface::over(
            Plane::through(origin, Direction::new(normal, T).unwrap()),
            (-6.0, 6.0),
            (-6.0, 6.0),
        )
        .unwrap()
        .into()
    }

    fn segment(from: Point, to: Point) -> Curve {
        LineCurve::segment(from, to, T).unwrap().into()
    }

    #[test]
    fn a_line_through_a_sphere_pierces_it_where_the_quadratic_says() {
        let ball = sphere(2.0);
        let ray = segment(Point::new(-5.0, 0.0, 0.0), Point::new(5.0, 0.0, 0.0));
        let found =
            intersect_curve_surface(&ray, &ball, CurveSurfaceOptions::default(), T).unwrap();
        assert_eq!(found.crossings.len(), 2);
        assert!(
            found.crossings[0]
                .point
                .is_equal(Point::new(-2.0, 0.0, 0.0), T)
        );
        assert!(
            found.crossings[1]
                .point
                .is_equal(Point::new(2.0, 0.0, 0.0), T)
        );
        for hit in &found.crossings {
            assert!(hit.gap < 1e-12);
            // The surface parameters evaluate back onto the point.
            let lifted = ball
                .point_at(hit.on_surface.0, hit.on_surface.1, T)
                .unwrap();
            assert!(lifted.is_equal(hit.point, T));
        }

        // Tangent: one root. Missing: none.
        let grazing = segment(Point::new(-5.0, 0.0, 2.0), Point::new(5.0, 0.0, 2.0));
        assert_eq!(
            intersect_curve_surface(&grazing, &ball, CurveSurfaceOptions::default(), T)
                .unwrap()
                .crossings
                .len(),
            1
        );
        let missing = segment(Point::new(-5.0, 0.0, 3.0), Point::new(5.0, 0.0, 3.0));
        assert!(
            intersect_curve_surface(&missing, &ball, CurveSurfaceOptions::default(), T)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn a_line_through_a_cylinder_respects_its_height() {
        let drum = cylinder(2.0, (-1.0, 1.0));
        // Crosses the infinite cylinder at z = 0: inside the height, two hits.
        let level = segment(Point::new(-5.0, 0.0, 0.0), Point::new(5.0, 0.0, 0.0));
        assert_eq!(
            intersect_curve_surface(&level, &drum, CurveSurfaceOptions::default(), T)
                .unwrap()
                .crossings
                .len(),
            2
        );
        // Crosses at z = 3: the unbounded geometry meets it, this surface
        // does not reach there.
        let high = segment(Point::new(-5.0, 0.0, 3.0), Point::new(5.0, 0.0, 3.0));
        assert!(
            intersect_curve_surface(&high, &drum, CurveSurfaceOptions::default(), T)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn a_line_lying_in_a_plane_is_an_overlap_not_a_crossing_list() {
        let ground = plane(Point::ORIGIN, Vector::Z);
        let lying = segment(Point::new(-3.0, 1.0, 0.0), Point::new(3.0, 1.0, 0.0));
        let found =
            intersect_curve_surface(&lying, &ground, CurveSurfaceOptions::default(), T).unwrap();
        assert!(found.crossings.is_empty());
        assert_eq!(found.lying.len(), 1);

        let crossing = segment(Point::new(0.0, 0.0, -1.0), Point::new(0.0, 0.0, 1.0));
        let found =
            intersect_curve_surface(&crossing, &ground, CurveSurfaceOptions::default(), T).unwrap();
        assert_eq!(found.crossings.len(), 1);
        assert!(found.crossings[0].point.is_equal(Point::ORIGIN, T));

        let parallel = segment(Point::new(-3.0, 0.0, 1.0), Point::new(3.0, 0.0, 1.0));
        assert!(
            intersect_curve_surface(&parallel, &ground, CurveSurfaceOptions::default(), T)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn a_circle_pierces_a_plane_twice_through_the_general_path() {
        // A circle in the xz plane against the ground: no analytic case
        // handles circle/plane here, so this is the seeded Newton path, and
        // the answer is known exactly anyway.
        let ring: Curve = CircleCurve::new(
            Circle::new(
                Frame::new(Point::new(0.0, 0.0, 0.0), -Direction::Y, Direction::X, T).unwrap(),
                2.0,
                T,
            )
            .unwrap(),
        )
        .into();
        let ground = plane(Point::ORIGIN, Vector::Z);
        let found =
            intersect_curve_surface(&ring, &ground, CurveSurfaceOptions::default(), T).unwrap();
        assert_eq!(found.crossings.len(), 2);
        for hit in &found.crossings {
            assert!(hit.gap < 1e-9);
            assert!(hit.point.z.abs() < 1e-9);
            assert!((hit.point.to_vector().magnitude() - 2.0).abs() < 1e-9);
        }
    }

    #[test]
    fn a_spline_through_a_sphere_is_found_and_polished() {
        // A spline wandering through the ball: piercings with no closed form
        // anywhere, verified implicitly — each reported point is on the
        // sphere to the gap it claims.
        let wander: Curve = BSplineCurve::new(
            KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 3).unwrap(),
            vec![
                Point::new(-4.0, -1.0, -1.0),
                Point::new(-1.0, 2.0, 1.0),
                Point::new(1.0, -2.0, -1.0),
                Point::new(4.0, 1.0, 1.0),
            ],
            T,
        )
        .unwrap()
        .into();
        let ball = sphere(2.0);
        let found =
            intersect_curve_surface(&wander, &ball, CurveSurfaceOptions::default(), T).unwrap();
        assert!(!found.crossings.is_empty(), "the spline passes through");
        for hit in &found.crossings {
            assert!(hit.gap < 1e-9);
            let SurfaceGeometry::Sphere(s) = &ball else {
                unreachable!()
            };
            assert!(s.sphere().distance_to(hit.point).abs() < 1e-9);
        }
    }

    #[test]
    fn unusable_options_are_refused() {
        let ball = sphere(1.0);
        let ray = segment(Point::new(-5.0, 0.0, 0.0), Point::new(5.0, 0.0, 0.0));
        for options in [
            CurveSurfaceOptions {
                samples: 1,
                ..CurveSurfaceOptions::default()
            },
            CurveSurfaceOptions {
                grid: 1,
                ..CurveSurfaceOptions::default()
            },
            CurveSurfaceOptions {
                gap: 0.0,
                ..CurveSurfaceOptions::default()
            },
        ] {
            assert!(intersect_curve_surface(&ray, &ball, options, T).is_err());
        }
    }
}
