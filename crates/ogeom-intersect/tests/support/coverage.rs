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

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Surface, SurfaceGeometry};
use ogeom_math::Point;

use ogeom_intersect::Traced;

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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the grid is too
/// coarse to have cells; [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the
/// second surface has no signed distance, which means this cannot say anything
/// about that pair rather than that the pair is wrong.
pub fn coverage(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    branches: &[Traced],
    grid: usize,
    tol: Tolerances,
) -> OgeomResult<Coverage> {
    if grid < 2 {
        ogeom_bail!(Construction, "a coverage grid needs at least two steps");
    }
    if signed_distance(b, Point::ORIGIN).is_none() {
        ogeom_bail!(
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
                    .fold(ogeom_math::Vector::ZERO, |acc, p| acc + p.to_vector())
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
