//! Measuring whether an intersection is *all* of the intersection.
//!
//! The second instrument the gate needs, and the one that can falsify the first.
//!
//! [`benchmark`](crate::benchmark) measures accuracy: every point of every
//! reported curve lies on both surfaces. That is necessary and it is not
//! sufficient — an intersector that found one of two circles and traced it
//! perfectly scores perfectly, and the answer is still half missing. Accuracy
//! cannot see a branch that was never looked at, because there is nothing wrong
//! with the points it did report.
//!
//! # How completeness is measured without knowing the answer
//!
//! By asking the *surfaces* where the intersection has to be, rather than asking
//! the intersector.
//!
//! Sample one surface into cells. For each cell, take the signed distance from
//! its corners to the other surface. Where that sign changes across a cell, the
//! other surface passes through it — no intersector was consulted to establish
//! that, it follows from the intermediate value theorem and two closed-form
//! distance functions. So the intersection curve crosses that cell, and some
//! reported branch had better cross it too.
//!
//! A cell that the surfaces say is crossed and no branch reaches is a **miss**,
//! and the instrument reports where it is rather than only how many there were.
//!
//! # What it cannot see
//!
//! A branch entirely inside one cell — a tiny loop where two surfaces barely
//! graze — changes no corner's sign and is invisible here, exactly as it is
//! invisible to the seeding it is checking. Both get better with a finer grid
//! and neither is exact, so a clean completeness score at one resolution is a
//! statement about that resolution. The score carries the grid it was measured
//! at for that reason.
//!
//! A tangential contact changes no sign either: the surfaces touch without
//! crossing. Those are counted separately rather than scored as misses, since a
//! curve is not what is there to be found.

use og_core::{OgResult, Tolerances, og_bail};
use og_geom::{Surface, SurfaceGeometry};
use og_math::Point;

use crate::march::Traced;

/// What a completeness measurement found.
#[derive(Debug, Clone, PartialEq)]
pub struct Coverage {
    /// Cells the two surfaces say the intersection crosses.
    pub crossings: usize,
    /// How many of those a reported branch reaches.
    pub covered: usize,
    /// A point in each cell that no branch reached.
    ///
    /// The positions, not just the count: a miss is worth looking at, and a
    /// number sends you hunting for it.
    pub missed: Vec<Point>,
    /// The grid the measurement was taken at.
    ///
    /// Carried with the result because completeness is only ever a statement
    /// about a resolution — a branch narrower than a cell is invisible to this
    /// and to the seeding it is checking.
    pub grid: usize,
}

impl Coverage {
    /// Whether every crossing found a branch.
    #[must_use]
    pub const fn complete(&self) -> bool {
        self.missed.is_empty()
    }

    /// The share of crossings covered, from zero to one.
    ///
    /// One when nothing was crossed, since nothing was missed either — the
    /// alternative is a division by zero dressed up as a failure.
    #[must_use]
    pub fn fraction(&self) -> f64 {
        if self.crossings == 0 {
            return 1.0;
        }
        #[allow(clippy::cast_precision_loss)]
        {
            self.covered as f64 / self.crossings as f64
        }
    }
}

/// Measure how much of the intersection a set of branches actually found.
///
/// `grid` is how finely the first surface is sampled. Both surfaces need a
/// *signed* distance for this to mean anything — the sign change is the whole
/// method — so a pair without one is refused rather than scored.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the grid is too
/// coarse to have cells; [`OgError::NotDone`](og_core::OgError::NotDone) if the
/// second surface has no signed distance, which means this cannot say anything
/// about that pair rather than that the pair is wrong.
pub fn coverage(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    branches: &[Traced],
    grid: usize,
    tol: Tolerances,
) -> OgResult<Coverage> {
    if grid < 2 {
        og_bail!(Construction, "a coverage grid needs at least two steps");
    }
    if signed_distance(b, Point::ORIGIN).is_none() {
        og_bail!(
            NotDone,
            "the second surface has no signed distance, so there is no sign to \
             change and this cannot measure completeness for the pair. That is \
             a limit of the instrument, not a finding about the intersector"
        );
    }

    let ((ua, ub), (va, vb)) = bounded(a);
    #[allow(clippy::cast_precision_loss)]
    let n = grid as f64;
    let mut crossings = 0;
    let mut covered = 0;
    let mut missed = Vec::new();

    for i in 0..grid {
        for j in 0..grid {
            #[allow(clippy::cast_precision_loss)]
            let (s0, s1) = (i as f64 / n, (i + 1) as f64 / n);
            #[allow(clippy::cast_precision_loss)]
            let (t0, t1) = (j as f64 / n, (j + 1) as f64 / n);
            let at = |s: f64, t: f64| a.point_at(ua + (ub - ua) * s, va + (vb - va) * t, tol).ok();
            let corners: Vec<Point> = [(s0, t0), (s1, t0), (s1, t1), (s0, t1)]
                .iter()
                .filter_map(|(s, t)| at(*s, *t))
                .collect();
            if corners.len() < 4 {
                continue;
            }
            let signs: Vec<f64> = corners
                .iter()
                .filter_map(|p| signed_distance(b, *p))
                .collect();
            if signs.len() < 4 {
                continue;
            }
            // The sign changing across the cell means the other surface passes
            // through it. Nothing was asked of the intersector to know that.
            let positive = signs.iter().any(|d| *d > 0.0);
            let negative = signs.iter().any(|d| *d < 0.0);
            if !(positive && negative) {
                continue;
            }
            crossings += 1;

            // A branch passing anywhere through the cell counts. Measured
            // against the polyline's *segments*, not its vertices: the tracer
            // spaces its points by curvature, so on a straight stretch they sit
            // far apart, and comparing to vertices alone marked the cells
            // between two samples of a perfectly traced line as missed.
            let centre = Point::from_vector(
                corners
                    .iter()
                    .fold(og_math::Vector::ZERO, |acc, p| acc + p.to_vector())
                    * 0.25,
            );
            let reach = corners[0]
                .distance(corners[2])
                .max(corners[1].distance(corners[3]));
            let reached = branches.iter().any(|branch| {
                branch
                    .points
                    .windows(2)
                    .any(|pair| segment_distance(centre, pair[0], pair[1]) <= reach)
                    || branch
                        .points
                        .first()
                        .is_some_and(|p| p.distance(centre) <= reach)
            });
            if reached {
                covered += 1;
            } else {
                missed.push(centre);
            }
        }
    }

    Ok(Coverage {
        crossings,
        covered,
        missed,
        grid,
    })
}

/// Distance from a point to a segment.
fn segment_distance(p: Point, a: Point, b: Point) -> f64 {
    let along = b - a;
    let length = along.square_magnitude();
    if length <= f64::MIN_POSITIVE {
        return p.distance(a);
    }
    let t = ((p - a).dot(along) / length).clamp(0.0, 1.0);
    p.distance(a + along * t)
}

/// A surface's domain, clamped to somewhere a real model lives.
///
/// An unbounded plane declares a domain reaching a billion units, and sampling
/// that puts every cell so far from the origin that nothing is ever found near
/// the intersection.
fn bounded(surface: &SurfaceGeometry) -> ((f64, f64), (f64, f64)) {
    let ((ua, ub), (va, vb)) = surface.domain();
    let limit = 1.0e6;
    (
        (ua.max(-limit), ub.min(limit)),
        (va.max(-limit), vb.min(limit)),
    )
}

/// The signed distance to a surface, where that means something.
///
/// `None` for a surface with no inside — a cone's two nappes and a spindle
/// torus's folded branch have no consistent side, and inventing one would put
/// a sign change where there is no crossing.
fn signed_distance(surface: &SurfaceGeometry, p: Point) -> Option<f64> {
    use SurfaceGeometry as S;
    match surface {
        S::Plane(x) => Some(x.plane().signed_distance_to(p)),
        S::Sphere(x) => Some(x.sphere().signed_distance_to(p)),
        S::Cylinder(x) => Some(x.cylinder().signed_distance_to(p)),
        _ => None,
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::print_stdout)]
mod tests {
    use super::*;
    use crate::march::{Marching, branches};
    use og_geom::{CylinderSurface, PlaneSurface, SphereSurface};
    use og_math::{Cylinder, Direction, Frame, Plane, Sphere, Vector};

    const T: Tolerances = Tolerances::millimetres();

    fn sphere(centre: Point, radius: f64) -> SurfaceGeometry {
        SphereSurface::new(Sphere::centred(centre, radius, T).unwrap()).into()
    }

    fn cylinder(axis: Vector, radius: f64) -> SurfaceGeometry {
        let frame = Frame::new(
            Point::ORIGIN,
            Direction::new(axis, T).unwrap(),
            Direction::from_cross(axis, Vector::new(0.3, 0.5, 0.9), T).unwrap(),
            T,
        )
        .unwrap();
        CylinderSurface::new(Cylinder::new(frame, radius, T).unwrap(), (-4.0, 4.0))
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

    fn options() -> Marching {
        Marching {
            chord: 1e-4,
            ..Marching::default()
        }
    }

    #[test]
    fn dropping_a_branch_is_caught() {
        // The test that makes the instrument worth having. A completeness
        // measure that always says "complete" looks exactly like a correct one
        // until something is actually missing, so the negative control is not
        // optional — it is the only evidence the thing works.
        let a = sphere(Point::ORIGIN, 3.0);
        let b = cylinder(Vector::Z, 1.5);

        let found = branches(&a, &b, options(), T).unwrap();
        assert_eq!(found.len(), 2, "a coaxial cylinder cuts a sphere twice");

        let whole = coverage(&a, &b, &found, 40, T).unwrap();
        assert!(
            whole.complete(),
            "{} of {} crossings covered, first miss at {:?}",
            whole.covered,
            whole.crossings,
            whole.missed.first()
        );

        // Now hide one, exactly as an intersector that never seeded it would.
        let half = coverage(&a, &b, &found[..1], 40, T).unwrap();
        assert!(
            !half.complete(),
            "dropping a whole branch went unnoticed, which means this measures \
             nothing"
        );
        println!(
            "coverage with one branch of two: {:.1}% ({} missed)",
            half.fraction() * 100.0,
            half.missed.len()
        );
        assert!(half.fraction() < 0.75, "got {}", half.fraction());

        // And reporting nothing at all is caught most of all.
        let none = coverage(&a, &b, &[], 40, T).unwrap();
        assert_eq!(none.covered, 0);
        assert!(none.crossings > 0);
    }

    #[test]
    fn the_marching_intersector_finds_all_of_what_it_is_asked_for() {
        // The measurement the gate wants, on the cases that have no closed
        // form. Accuracy was already established; this is the other half.
        let cases: Vec<(&str, SurfaceGeometry, SurfaceGeometry)> = vec![
            (
                "sphere/plane",
                sphere(Point::ORIGIN, 3.0),
                plane(Point::new(0.0, 0.0, 1.0), Vector::Z),
            ),
            (
                "sphere/cylinder coaxial",
                sphere(Point::ORIGIN, 3.0),
                cylinder(Vector::Z, 1.5),
            ),
            (
                "sphere/cylinder offset",
                sphere(Point::new(0.6, 0.0, 0.0), 3.0),
                cylinder(Vector::Z, 1.5),
            ),
            (
                "crossed cylinders",
                cylinder(Vector::Z, 1.0),
                cylinder(Vector::X, 1.6),
            ),
        ];
        for (name, a, b) in cases {
            let found = branches(&a, &b, options(), T).unwrap();
            let score = coverage(&a, &b, &found, 40, T).unwrap();
            println!(
                "coverage {name}: {}/{} cells, {} branches",
                score.covered,
                score.crossings,
                found.len()
            );
            assert!(score.crossings > 0, "{name}: nothing to cover");
            assert!(
                score.complete(),
                "{name}: missed {} of {} crossings, first at {:?}",
                score.crossings - score.covered,
                score.crossings,
                score.missed.first()
            );
        }
    }

    #[test]
    fn a_pair_it_cannot_measure_says_so_rather_than_scoring_full_marks() {
        // A cone has no consistent inside, so there is no sign to change. The
        // dangerous answer would be "complete", since a caller reading a
        // hundred percent has no way to tell it apart from a real result.
        let cone: SurfaceGeometry = og_geom::ConeSurface::new(
            og_math::Cone::new(Frame::WORLD, 1.0, 0.4_f64.atan(), T).unwrap(),
            (0.0, 3.0),
        )
        .unwrap()
        .into();
        let err = coverage(&sphere(Point::ORIGIN, 2.0), &cone, &[], 20, T).unwrap_err();
        assert!(err.to_string().contains("signed distance"), "got {err}");
    }

    #[test]
    fn nothing_crossed_is_complete_rather_than_a_division_by_zero() {
        let far = sphere(Point::new(100.0, 0.0, 0.0), 1.0);
        let score = coverage(&sphere(Point::ORIGIN, 1.0), &far, &[], 16, T).unwrap();
        assert_eq!(score.crossings, 0);
        assert!(score.complete());
        assert!((score.fraction() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn a_grid_too_coarse_to_have_cells_is_refused() {
        let a = sphere(Point::ORIGIN, 1.0);
        let b = plane(Point::ORIGIN, Vector::Z);
        assert!(coverage(&a, &b, &[], 1, T).is_err());
        assert!(coverage(&a, &b, &[], 0, T).is_err());
    }
}
