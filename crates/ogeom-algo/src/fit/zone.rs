//! Form as the standards define it: minimum-zone, maximum-inscribed and
//! minimum-circumscribed fits.
//!
//! A least-squares fit ([`plane`](super::plane), [`cylinder`](super::cylinder)
//! and the rest) makes the squared distances least, and its residuals'
//! peak to peak is an upper bound on form. The form a tolerance limits is
//! a Chebyshev (minimax) zone: the narrowest pair of parallel planes,
//! concentric circles, coaxial cylinders or concentric spheres that holds
//! every point. Each minimum zone here starts from the least-squares fit
//! and is never wider than its peak to peak.
//!
//! # Exact and certified
//!
//! - [`flatness`], [`circularity`], [`circumscribed_circle`] and
//!   [`inscribed_circle`] are exact: the global optimum, to rounding.
//!   The slab and the ring are found over a working subset of the points
//!   whose zone, tried in every configuration the optimum can take, is
//!   grown by the points outside it until it holds them all; the
//!   enclosing circle is Welzl's; the empty circle stands on the points'
//!   Voronoi diagram.
//! - [`sphericity`], [`cylindricity`], [`straightness`],
//!   [`circumscribed_cylinder`] and [`inscribed_cylinder`] are not convex
//!   problems. They are settled by sequential linear programming in a trust
//!   region from the least-squares fit, which converges to a point where no
//!   first-order change of the parameters narrows the zone, and
//!   [`Zone::optimal`] says whether the result meets that test: the
//!   contacts' gradients balance, so the zone is least within the
//!   confusion distance among shapes near it.
//!
//! The fits run on the calling thread, in point order, and draw no
//! randomness but a fixed sequence: the same points give the same answer
//! bit for bit.

mod exact;
mod hull;
mod lp;
mod slp;

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Axis, Circle2, Cylinder, Direction, Frame, Plane, Point, Point2, Sphere};

use super::FitOptions;
use exact::Ring;
use slp::{Centre, Goal, Line};

/// A fitted zone.
#[derive(Debug, Clone, PartialEq)]
pub struct Zone<S> {
    /// The zone's middle: the plane halfway between the two planes, the
    /// circle, cylinder or sphere halfway between the two. For an inscribed
    /// or circumscribed fit, the fitted circle or cylinder itself; for
    /// straightness, the axis of the cylinder.
    pub shape: S,
    /// The zone's width: the distance between its two boundaries. For an
    /// inscribed or circumscribed fit, the points' greatest distance from
    /// the centre or axis less their least, their form about that
    /// reference; for straightness, the cylinder's diameter.
    pub width: f64,
    /// The points touching the zone's boundaries, within the confusion
    /// distance, in index order.
    pub contacts: Vec<usize>,
    /// Whether the zone is the optimum: always for the exact fits, and for
    /// the others where the first-order Chebyshev conditions hold.
    pub optimal: bool,
}

/// Flatness: the narrowest slab between two parallel planes holding every
/// point. Exact.
///
/// # Errors
///
/// As [`plane`](super::plane): too few points, points on one line, or a
/// point that is not finite.
pub fn flatness(points: &[Point], tol: Tolerances) -> OgeomResult<Zone<Plane>> {
    let fit = super::plane(points, &FitOptions::default(), tol)?;
    let origin = fit.shape.origin();
    let normal = fit.shape.normal().vector();
    let mut seeds = extremes(&fit.residuals);
    seeds.extend(spanning(points));
    let Some(slab) = exact::narrowest_slab(points, origin, &seeds) else {
        ogeom_bail!(Construction, "the points lie on one line and span no plane");
    };
    let (n, low, high) = if slab.normal.dot(normal) < 0.0 {
        (-slab.normal, -slab.high, -slab.low)
    } else {
        (slab.normal, slab.low, slab.high)
    };
    let slack = tol.confusion();
    let contacts = touching(
        points.iter().map(|p| n.dot(*p - origin)),
        Some(low),
        Some(high),
        slack,
    );
    let middle = origin + n * ((low + high) / 2.0);
    Ok(Zone {
        shape: Plane::new(Frame::about(middle, Direction::new(n, tol)?)),
        width: high - low,
        contacts,
        optimal: true,
    })
}

/// Straightness of a line in space: the thinnest cylinder about a line
/// holding every point, its axis and its diameter.
///
/// # Errors
///
/// As [`line3`](super::line3).
pub fn straightness(points: &[Point], tol: Tolerances) -> OgeomResult<Zone<Axis>> {
    let fit = super::line3(points, &FitOptions::default(), tol)?;
    let start = line(fit.shape.location, fit.shape.direction, points)?;
    let settled = slp::settle(start, points, Goal::Outer, reach(points), tol.confusion());
    let (_, most) = span(&settled.distances);
    let axis = Axis::new(settled.shape.origin, Direction::new(settled.shape.z, tol)?);
    Ok(Zone {
        shape: axis,
        width: 2.0 * most,
        contacts: touching(
            settled.distances.iter().copied(),
            None,
            Some(most),
            tol.confusion(),
        ),
        optimal: settled.optimal,
    })
}

/// Circularity: the narrowest ring between two concentric circles holding
/// every point. Exact.
///
/// # Errors
///
/// As [`circle2`](super::circle2), and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// the points lie nearer a line than any circle: a strip between two
/// parallel lines holds them narrower than every ring, and the narrowest
/// ring is that strip's limit, not reached.
pub fn circularity(points: &[Point2], tol: Tolerances) -> OgeomResult<Zone<Circle2>> {
    let fit = super::circle2(points, &FitOptions::default(), tol)?;
    let mut seeds = extremes(&fit.residuals);
    let lifted: Vec<Point> = points.iter().map(|p| p.to_3d()).collect();
    seeds.extend(spanning(&lifted));
    match exact::narrowest_ring(points, fit.shape.centre(), &seeds) {
        Ring::Circles {
            centre,
            radius,
            deviations,
            inner,
            outer,
        } => Ok(Zone {
            shape: Circle2::centred(centre, radius + (inner + outer) / 2.0, tol)?,
            width: outer - inner,
            contacts: touching(
                deviations.into_iter(),
                Some(inner),
                Some(outer),
                tol.confusion(),
            ),
            optimal: true,
        }),
        Ring::Straight { width } => ogeom_bail!(
            Construction,
            "the points lie between two lines {width} apart, narrower than any ring: \
             they have no narrowest circle"
        ),
    }
}

/// Cylindricity: the narrowest shell between two coaxial cylinders holding
/// every point.
///
/// # Errors
///
/// As [`cylinder`](super::cylinder).
pub fn cylindricity(points: &[Point], tol: Tolerances) -> OgeomResult<Zone<Cylinder>> {
    revolved(points, Goal::Zone, tol)
}

/// The largest cylinder holding no point inside, from the least-squares
/// cylinder: a bore's functional size.
///
/// # Errors
///
/// As [`cylinder`](super::cylinder).
pub fn inscribed_cylinder(points: &[Point], tol: Tolerances) -> OgeomResult<Zone<Cylinder>> {
    revolved(points, Goal::Inner, tol)
}

/// The smallest cylinder holding every point, from the least-squares
/// cylinder: a shaft's functional size.
///
/// # Errors
///
/// As [`cylinder`](super::cylinder).
pub fn circumscribed_cylinder(points: &[Point], tol: Tolerances) -> OgeomResult<Zone<Cylinder>> {
    revolved(points, Goal::Outer, tol)
}

/// Sphericity: the narrowest shell between two concentric spheres holding
/// every point.
///
/// # Errors
///
/// As [`sphere`](super::sphere).
pub fn sphericity(points: &[Point], tol: Tolerances) -> OgeomResult<Zone<Sphere>> {
    let fit = super::sphere(points, &FitOptions::default(), tol)?;
    let slack = tol.confusion();
    let settled = slp::settle(
        Centre(fit.shape.centre()),
        points,
        Goal::Zone,
        reach(points),
        slack,
    );
    let (least, most) = span(&settled.distances);
    Ok(Zone {
        shape: Sphere::centred(settled.shape.0, (least + most) / 2.0, tol)?,
        width: most - least,
        contacts: touching(
            settled.distances.iter().copied(),
            Some(least),
            Some(most),
            slack,
        ),
        optimal: settled.optimal,
    })
}

/// The largest circle whose centre lies in the points' convex hull and
/// which holds none of them inside. Exact.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) with
/// fewer than three points or with them all on one line, and
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for a point that
/// is not finite.
pub fn inscribed_circle(points: &[Point2], tol: Tolerances) -> OgeomResult<Zone<Circle2>> {
    finite(points)?;
    if points.len() < 3 {
        ogeom_bail!(
            Construction,
            "an inscribed circle needs at least 3 points, got {}",
            points.len()
        );
    }
    let Some(centre) = exact::largest_empty_circle(points) else {
        ogeom_bail!(
            Construction,
            "the points lie on one line and enclose no circle"
        );
    };
    let distances: Vec<f64> = points.iter().map(|p| p.distance(centre)).collect();
    let (least, most) = span(&distances);
    Ok(Zone {
        shape: Circle2::centred(centre, least, tol)?,
        width: most - least,
        contacts: touching(
            distances.iter().copied(),
            Some(least),
            None,
            tol.confusion(),
        ),
        optimal: true,
    })
}

/// The smallest circle holding every point. Exact.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) with no
/// points or with them all within the confusion distance of one, and
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for a point that
/// is not finite.
pub fn circumscribed_circle(points: &[Point2], tol: Tolerances) -> OgeomResult<Zone<Circle2>> {
    finite(points)?;
    if points.is_empty() {
        ogeom_bail!(Construction, "a circumscribed circle needs a point");
    }
    let (centre, radius) = exact::enclosing_circle(points);
    let distances: Vec<f64> = points.iter().map(|p| p.distance(centre)).collect();
    let (least, most) = span(&distances);
    Ok(Zone {
        shape: Circle2::centred(centre, radius, tol)?,
        width: most - least,
        contacts: touching(distances.iter().copied(), None, Some(most), tol.confusion()),
        optimal: true,
    })
}

/// A cylinder's zone of `goal`, settled from the least-squares cylinder.
fn revolved(points: &[Point], goal: Goal, tol: Tolerances) -> OgeomResult<Zone<Cylinder>> {
    let fit = super::cylinder(points, &FitOptions::default(), tol)?;
    let frame = fit.shape.frame();
    let start = line(frame.origin(), frame.z(), points)?;
    let slack = tol.confusion();
    let settled = slp::settle(start, points, goal, reach(points), slack);
    let (least, most) = span(&settled.distances);
    let (radius, low, high) = match goal {
        Goal::Zone => ((least + most) / 2.0, Some(least), Some(most)),
        Goal::Outer => (most, None, Some(most)),
        Goal::Inner => (least, Some(least), None),
    };
    let axis = Frame::about(settled.shape.origin, Direction::new(settled.shape.z, tol)?);
    Ok(Zone {
        shape: Cylinder::new(axis, radius, tol)?,
        width: most - least,
        contacts: touching(settled.distances.iter().copied(), low, high, slack),
        optimal: settled.optimal,
    })
}

/// The line through `location` along `direction`, its origin moved level
/// with the points' centroid.
fn line(location: Point, direction: Direction, points: &[Point]) -> OgeomResult<Line> {
    let z = direction.vector();
    let centroid = Point::centroid(points)?;
    let origin = location + z * (centroid - location).dot(z);
    match Line::new(origin, z, points) {
        Some(l) => Ok(l),
        None => ogeom_bail!(Construction, "the fitted axis has no direction"),
    }
}

/// The diagonal of the points' bounding box: the widest step an
/// optimality test need consider.
fn reach(points: &[Point]) -> f64 {
    let (mut lo, mut hi) = (points[0], points[0]);
    for p in points {
        lo = lo.min(*p);
        hi = hi.max(*p);
    }
    lo.distance(hi)
}

fn span(distances: &[f64]) -> (f64, f64) {
    distances
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), d| {
            (lo.min(*d), hi.max(*d))
        })
}

/// The points within `slack` of the boundaries `low` and `high`, or past
/// them.
fn touching(
    values: impl Iterator<Item = f64>,
    low: Option<f64>,
    high: Option<f64>,
    slack: f64,
) -> Vec<usize> {
    values
        .enumerate()
        .filter(|&(_, v)| {
            low.is_some_and(|l| v <= l + slack) || high.is_some_and(|h| v >= h - slack)
        })
        .map(|(i, _)| i)
        .collect()
}

/// The points of the four greatest residuals and the four least: where the
/// zone's contacts are first looked for.
fn extremes(residuals: &[f64]) -> Vec<usize> {
    let mut out = slp::largest(residuals.len(), 4, |i| Some(residuals[i]));
    out.extend(slp::largest(residuals.len(), 4, |i| Some(-residuals[i])));
    out
}

/// Three points spanning the most: the farthest from the first, the
/// farthest from that, and the farthest from the line through those two.
fn spanning(points: &[Point]) -> Vec<usize> {
    let farthest = |key: &dyn Fn(Point) -> f64| {
        slp::largest(points.len(), 1, |i| Some(key(points[i])))
            .first()
            .copied()
            .unwrap_or(0)
    };
    let a = farthest(&|p| p.distance(points[0]));
    let b = farthest(&|p| p.distance(points[a]));
    let along = points[b] - points[a];
    let c = farthest(&|p| (p - points[a]).cross(along).magnitude());
    vec![a, b, c]
}

fn finite(points: &[Point2]) -> OgeomResult<()> {
    if let Some(i) = points
        .iter()
        .position(|p| !(p.x.is_finite() && p.y.is_finite()))
    {
        ogeom_bail!(Domain, "point {i} is not finite");
    }
    Ok(())
}
