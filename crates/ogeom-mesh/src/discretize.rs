//! Turning a curve into a polyline within a stated deflection.
//!
//! The first half of tessellation, and the half that decides whether the second
//! half can succeed: a face's triangulation is built from the polylines of its
//! boundary edges, so two faces meeting along an edge produce a watertight join
//! only if they discretize that edge to the *same* points. That is why
//! discretization is a property of the edge rather than of the face, and why it
//! is a separate step from triangulating the face itself.
//!
//! # What deflection means
//!
//! The chord deflection is the greatest distance between the polyline and the
//! curve it approximates. It is a length, so it scales with the model, and it
//! is the number a caller actually has an opinion about: "no visible error at
//! this zoom", "within machining tolerance".
//!
//! The angular deflection bounds how far the tangent may turn across one
//! segment. Without it a nearly straight curve gets two points and a
//! near-circular one gets far too few near its flattest part — chord error
//! alone does not notice a long, gently curving span.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d};
use ogeom_math::Point;

/// How finely to approximate a curve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Deflection {
    /// The greatest allowed distance from the polyline to the curve.
    pub chord: f64,
    /// The greatest allowed turn in the tangent across one segment, in radians.
    pub angular: f64,
    /// A floor on the number of segments for a curve that bends.
    ///
    /// A closed curve needs at least three to enclose anything, and a nearly
    /// straight arc of one would otherwise collapse to a chord that misses the
    /// bulge entirely — the midpoint test is a sample, and one sample can be
    /// placed exactly where the curve happens to cross its own chord.
    ///
    /// It does *not* apply to a straight curve, which is exactly represented by
    /// its endpoints. Splitting a line adds points that no tolerance asked for,
    /// and every one of them becomes a vertex in every face the edge bounds.
    pub min_segments: usize,
    /// A ceiling, so a pathological curve cannot exhaust memory.
    ///
    /// Reaching it is reported rather than passed off as success — see
    /// [`Polyline::deflection_met`].
    pub max_segments: usize,
}

impl Default for Deflection {
    fn default() -> Self {
        Self {
            // A tenth of a millimetre at unit scale: invisible on screen, and
            // fine enough that a mass property computed from it is accurate to
            // roughly one part in a thousand.
            chord: 1e-1,
            // About 11 degrees, which keeps a circle at 32 segments or so.
            angular: 0.2,
            min_segments: 2,
            max_segments: 4096,
        }
    }
}

impl Deflection {
    /// A deflection with a given chord tolerance and the default angular one.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `chord` is
    /// not finite and positive.
    pub fn with_chord(chord: f64) -> OgeomResult<Self> {
        if !chord.is_finite() || chord <= 0.0 {
            ogeom_bail!(
                Construction,
                "chord deflection {chord} must be finite and positive"
            );
        }
        Ok(Self {
            chord,
            ..Self::default()
        })
    }

    /// A deflection scaled to a model of the given size.
    ///
    /// Relative tolerances are what a caller usually means: "a thousandth of the
    /// part" is a statement that survives the part being modelled in metres
    /// rather than millimetres, and an absolute default does not.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `size` or
    /// `fraction` is not finite and positive.
    pub fn relative(size: f64, fraction: f64) -> OgeomResult<Self> {
        if !size.is_finite() || size <= 0.0 || !fraction.is_finite() || fraction <= 0.0 {
            ogeom_bail!(
                Construction,
                "relative deflection needs a positive size and fraction, got {size} and {fraction}"
            );
        }
        Self::with_chord(size * fraction)
    }

    /// Check that the settings are usable.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a tolerance
    /// is not positive, or the segment bounds are inconsistent.
    pub fn validate(&self) -> OgeomResult<()> {
        if !self.chord.is_finite() || self.chord <= 0.0 {
            ogeom_bail!(
                Construction,
                "chord deflection {} must be positive",
                self.chord
            );
        }
        if !self.angular.is_finite() || self.angular <= 0.0 {
            ogeom_bail!(
                Construction,
                "angular deflection {} must be positive",
                self.angular
            );
        }
        if self.min_segments == 0 {
            ogeom_bail!(Construction, "a polyline needs at least one segment");
        }
        if self.max_segments < self.min_segments {
            ogeom_bail!(
                Construction,
                "segment ceiling {} is below the floor {}",
                self.max_segments,
                self.min_segments
            );
        }
        Ok(())
    }
}

/// A curve approximated by points, with the parameters they came from.
///
/// The parameters are kept, not discarded, because the pcurve of the same edge
/// on an adjacent face has to be sampled at exactly these values for the two
/// faces to meet.
#[derive(Debug, Clone, PartialEq)]
pub struct Polyline {
    /// The points, in order along the curve.
    pub points: Vec<Point>,
    /// The parameter each point came from.
    pub parameters: Vec<f64>,
    /// Whether the requested deflection was actually achieved.
    ///
    /// `false` when the segment ceiling was reached first. Reporting it beats
    /// silently returning a coarser polyline than asked for, which would make
    /// a downstream tolerance claim untrue with nothing to show why.
    pub deflection_met: bool,
}

impl Polyline {
    /// Number of segments.
    #[must_use]
    pub fn segment_count(&self) -> usize {
        self.points.len().saturating_sub(1)
    }

    /// Total length of the polyline.
    ///
    /// An underestimate of the curve's own length, since a chord is shorter
    /// than the arc it spans — and one that improves as the deflection tightens.
    #[must_use]
    pub fn length(&self) -> f64 {
        self.points.windows(2).map(|w| w[0].distance(w[1])).sum()
    }

    /// Whether the polyline returns to where it started.
    #[must_use]
    pub fn is_closed(&self, tol: Tolerances) -> bool {
        match (self.points.first(), self.points.last()) {
            (Some(a), Some(b)) => self.points.len() > 2 && a.is_equal(*b, tol),
            _ => false,
        }
    }
}

/// Whether a curve is a straight line, seeing through any trimming.
///
/// A line is its own polyline, so subdividing it is pure cost: the extra points
/// land on the curve, satisfy every tolerance, and then propagate into the
/// triangulation of every face the edge bounds.
#[must_use]
pub fn is_straight(curve: &Curve) -> bool {
    match curve {
        Curve::Line(_) => true,
        Curve::Trimmed(t) => is_straight(t.basis()),
        _ => false,
    }
}

/// Whether a planar curve is a straight line in parameter space.
#[must_use]
pub fn is_straight_planar(curve: &ogeom_geom::PlanarCurve) -> bool {
    match curve {
        ogeom_geom::PlanarCurve::Line(_) => true,
        ogeom_geom::PlanarCurve::Trimmed(t) => is_straight_planar(t.basis()),
        _ => false,
    }
}

/// Approximate `curve` over `range` within `deflection`.
///
/// Adaptive bisection: split a segment whenever its midpoint is further from the
/// chord than allowed, or the tangent turns too far across it. Uniform sampling
/// is the obvious alternative and wastes points on the straight parts of a curve
/// while still missing the tight ones — the whole difficulty of tessellation is
/// that curvature is not uniform.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the deflection
/// settings are unusable or the range is empty;
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if the range leaves the curve.
pub fn discretize(
    curve: &Curve,
    range: (f64, f64),
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Polyline> {
    deflection.validate()?;
    let (lo, hi) = range;
    if !lo.is_finite() || !hi.is_finite() || hi <= lo + tol.parametric() {
        ogeom_bail!(Construction, "range [{lo}, {hi}] is empty");
    }

    // Start from the floor, so a closed curve is never approximated by a single
    // chord from a point back to itself. A straight curve is exempt: its
    // endpoints already represent it exactly.
    let start = if is_straight(curve) {
        1
    } else {
        deflection.min_segments
    };
    let mut parameters: Vec<f64> = (0..=start)
        .map(|i| {
            #[allow(clippy::cast_precision_loss)]
            let t = i as f64 / start as f64;
            lo + (hi - lo) * t
        })
        .collect();
    let mut points: Vec<Point> = parameters
        .iter()
        .map(|u| curve.point_at(*u, tol))
        .collect::<OgeomResult<_>>()?;

    let mut met = true;
    loop {
        let mut split_at: Option<usize> = None;
        for i in 0..points.len() - 1 {
            if needs_split(
                curve,
                (parameters[i], parameters[i + 1]),
                (points[i], points[i + 1]),
                deflection,
                tol,
            )? {
                split_at = Some(i);
                break;
            }
        }
        let Some(i) = split_at else { break };

        if points.len() > deflection.max_segments {
            met = false;
            break;
        }
        let mid = f64::midpoint(parameters[i], parameters[i + 1]);
        // A split that does not actually divide the interval means the
        // parameters have reached the resolution of f64; refining further would
        // loop without improving anything.
        if mid <= parameters[i] || mid >= parameters[i + 1] {
            met = false;
            break;
        }
        parameters.insert(i + 1, mid);
        points.insert(i + 1, curve.point_at(mid, tol)?);
    }

    Ok(Polyline {
        points,
        parameters,
        deflection_met: met,
    })
}

/// Whether one segment violates either tolerance.
fn needs_split(
    curve: &Curve,
    parameters: (f64, f64),
    ends: (Point, Point),
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let (a, b) = parameters;
    let mid = f64::midpoint(a, b);
    let on_curve = curve.point_at(mid, tol)?;

    // Chord error, measured at the midpoint. Not a bound on the true maximum
    // deviation, which would need the curve's second derivative over the span;
    // it is the standard estimate, and bisection drives it down regardless.
    let chord = ogeom_math::Axis::through(ends.0, ends.1, tol).map_or_else(
        |_| ends.0.distance(on_curve),
        |axis| axis.distance_to(on_curve),
    );
    if chord > deflection.chord {
        return Ok(true);
    }

    // Angular error: how far the tangent turns across the segment. Chord error
    // alone misses a long, gently curving span, which is exactly where a
    // silhouette goes visibly polygonal.
    let (Ok(start), Ok(end)) = (curve.tangent_at(a, tol), curve.tangent_at(b, tol)) else {
        // A cusp has no tangent to compare; the chord test still governs.
        return Ok(false);
    };
    Ok(start.angle(end) > deflection.angular)
}

/// Approximate a planar curve in a surface's parameter space.
///
/// The deflection is measured in parameter units here, not in space, so a caller
/// wanting a spatial tolerance has to convert through the surface's own scale —
/// the two differ by orders of magnitude near a pole. This exists for boundary
/// work in parameter space; for a face's actual boundary, discretize the edge's
/// 3D curve and evaluate the pcurve at those parameters instead, which is what
/// keeps adjacent faces watertight.
///
/// # Errors
///
/// As [`discretize`].
/// Discretize a pcurve with its chord tolerance measured *in space*,
/// through the surface that gives its chart a metric.
///
/// [`discretize_planar`] measures in parameter units because it has no
/// surface to convert through; this is the version that does. Each candidate
/// segment is lifted to the surface and the sagitta measured between world
/// points, so one chord tolerance means one thing whatever the chart's
/// scale — a quarter-turn on a large cylinder refines further than the same
/// quarter-turn on a small one.
///
/// # Errors
///
/// As [`discretize_planar`].
pub fn discretize_on_surface(
    curve: &ogeom_geom::PlanarCurve,
    range: (f64, f64),
    surface: &ogeom_geom::SurfaceGeometry,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<(Vec<ogeom_math::Point2>, Vec<f64>)> {
    use ogeom_geom::Curve2d;
    use ogeom_geom::Surface as _;

    deflection.validate()?;
    let (lo, hi) = range;
    if !lo.is_finite() || !hi.is_finite() || hi <= lo + tol.parametric() {
        ogeom_bail!(Construction, "range [{lo}, {hi}] is empty");
    }
    let lift = |uv: ogeom_math::Point2| -> OgeomResult<ogeom_math::Point> {
        surface.point_at(uv.x, uv.y, tol)
    };

    let start = if is_straight_planar(curve) {
        1
    } else {
        deflection.min_segments
    };
    let mut parameters: Vec<f64> = (0..=start)
        .map(|i| {
            #[allow(clippy::cast_precision_loss)]
            let t = i as f64 / start as f64;
            lo + (hi - lo) * t
        })
        .collect();
    let mut points: Vec<ogeom_math::Point2> = parameters
        .iter()
        .map(|u| curve.point_at(*u, tol))
        .collect::<OgeomResult<_>>()?;
    let mut lifted: Vec<ogeom_math::Point> = points
        .iter()
        .map(|uv| lift(*uv))
        .collect::<OgeomResult<_>>()?;

    while points.len() <= deflection.max_segments {
        let mut split_at = None;
        for i in 0..points.len() - 1 {
            let mid = f64::midpoint(parameters[i], parameters[i + 1]);
            let on_curve = curve.point_at(mid, tol)?;
            let in_space = lift(on_curve)?;
            // Sagitta in world units: the lifted midpoint against the lifted
            // chord.
            let (a, b) = (lifted[i], lifted[i + 1]);
            let chord_vector = b - a;
            let length = chord_vector.magnitude();
            let sagitta = if length <= tol.confusion() {
                in_space.distance(a)
            } else {
                (in_space - a).cross(chord_vector).magnitude() / length
            };
            if sagitta > deflection.chord {
                split_at = Some((i, mid, on_curve, in_space));
                break;
            }
        }
        let Some((i, mid, on_curve, in_space)) = split_at else {
            break;
        };
        parameters.insert(i + 1, mid);
        points.insert(i + 1, on_curve);
        lifted.insert(i + 1, in_space);
    }
    Ok((points, parameters))
}

pub fn discretize_planar(
    curve: &ogeom_geom::PlanarCurve,
    range: (f64, f64),
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<(Vec<ogeom_math::Point2>, Vec<f64>)> {
    use ogeom_geom::Curve2d;

    deflection.validate()?;
    let (lo, hi) = range;
    if !lo.is_finite() || !hi.is_finite() || hi <= lo + tol.parametric() {
        ogeom_bail!(Construction, "range [{lo}, {hi}] is empty");
    }

    let start = if is_straight_planar(curve) {
        1
    } else {
        deflection.min_segments
    };
    let mut parameters: Vec<f64> = (0..=start)
        .map(|i| {
            #[allow(clippy::cast_precision_loss)]
            let t = i as f64 / start as f64;
            lo + (hi - lo) * t
        })
        .collect();
    let mut points: Vec<ogeom_math::Point2> = parameters
        .iter()
        .map(|u| curve.point_at(*u, tol))
        .collect::<OgeomResult<_>>()?;

    while points.len() <= deflection.max_segments {
        let mut split_at = None;
        for i in 0..points.len() - 1 {
            let mid = f64::midpoint(parameters[i], parameters[i + 1]);
            let on_curve = curve.point_at(mid, tol)?;
            let chord = ogeom_math::Axis2::through(points[i], points[i + 1], tol).map_or_else(
                |_| points[i].distance(on_curve),
                |axis| axis.distance_to(on_curve),
            );
            if chord > deflection.chord {
                split_at = Some(i);
                break;
            }
        }
        let Some(i) = split_at else { break };
        let mid = f64::midpoint(parameters[i], parameters[i + 1]);
        if mid <= parameters[i] || mid >= parameters[i + 1] {
            break;
        }
        parameters.insert(i + 1, mid);
        points.insert(i + 1, curve.point_at(mid, tol)?);
    }

    Ok((points, parameters))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ogeom_geom::{BSplineCurve, CircleCurve, LineCurve};
    use ogeom_math::{Circle, Frame, KnotVector};

    const T: Tolerances = Tolerances::millimetres();

    fn circle(radius: f64) -> Curve {
        CircleCurve::new(Circle::new(Frame::WORLD, radius, T).unwrap()).into()
    }

    /// The greatest distance from the curve to the polyline, sampled densely.
    fn worst_error(curve: &Curve, line: &Polyline) -> f64 {
        let mut worst: f64 = 0.0;
        for window in line.parameters.windows(2) {
            for k in 1..16 {
                let t = f64::from(k) / 16.0;
                let u = window[0] + (window[1] - window[0]) * t;
                let on_curve = curve.point_at(u, T).unwrap();
                let a = curve.point_at(window[0], T).unwrap();
                let b = curve.point_at(window[1], T).unwrap();
                let chord = ogeom_math::Axis::through(a, b, T)
                    .map_or(0.0, |axis| axis.distance_to(on_curve));
                worst = worst.max(chord);
            }
        }
        worst
    }

    #[test]
    fn straightness_sees_through_a_trim() {
        // The exemption has to survive trimming, or an edge built from a
        // trimmed line — which is what a solid's edges usually are — pays the
        // floor anyway and the saving evaporates.
        use ogeom_geom::TrimmedCurve;
        let line: Curve = LineCurve::segment(Point::ORIGIN, Point::new(10.0, 0.0, 0.0), T)
            .unwrap()
            .into();
        assert!(is_straight(&line));
        assert!(is_straight(
            &TrimmedCurve::new(line, 2.0, 8.0, T).unwrap().into()
        ));
        assert!(!is_straight(&circle(1.0)));
    }

    #[test]
    fn a_line_needs_no_more_than_the_minimum_segments() {
        // A straight curve has no chord error and no turn, so refinement must
        // stop immediately rather than subdividing to the ceiling.
        let curve: Curve = LineCurve::segment(Point::ORIGIN, Point::new(100.0, 0.0, 0.0), T)
            .unwrap()
            .into();
        let line = discretize(&curve, curve.domain(), Deflection::default(), T).unwrap();
        assert_eq!(line.segment_count(), 1, "a line is its own polyline");
        assert!(line.deflection_met);
        assert_relative_eq!(line.length(), 100.0, epsilon = 1e-9);
    }

    #[test]
    fn a_circle_is_refined_until_the_chord_tolerance_is_met() {
        let curve = circle(10.0);
        for chord in [1.0_f64, 0.1, 0.01, 0.001] {
            let deflection = Deflection {
                chord,
                ..Deflection::default()
            };
            let line = discretize(&curve, curve.domain(), deflection, T).unwrap();
            assert!(line.deflection_met, "gave up at chord {chord}");
            assert!(
                worst_error(&curve, &line) <= chord * 1.5,
                "chord {chord}: worst error {}",
                worst_error(&curve, &line)
            );
        }
    }

    #[test]
    fn a_tighter_tolerance_always_gives_at_least_as_many_segments() {
        let curve = circle(10.0);
        let mut previous = 0;
        for chord in [2.0_f64, 1.0, 0.5, 0.1, 0.01] {
            let line = discretize(
                &curve,
                curve.domain(),
                Deflection {
                    chord,
                    ..Deflection::default()
                },
                T,
            )
            .unwrap();
            assert!(
                line.segment_count() >= previous,
                "chord {chord} gave fewer segments than a looser one"
            );
            previous = line.segment_count();
        }
    }

    #[test]
    fn a_polylines_length_underestimates_the_curve_and_converges_to_it() {
        // Each chord is shorter than the arc it spans, so the polyline is always
        // short — and refining closes the gap.
        let radius = 10.0;
        let curve = circle(radius);
        let exact = core::f64::consts::TAU * radius;

        let coarse = discretize(
            &curve,
            curve.domain(),
            Deflection {
                chord: 1.0,
                ..Deflection::default()
            },
            T,
        )
        .unwrap();
        let fine = discretize(
            &curve,
            curve.domain(),
            Deflection {
                chord: 1e-4,
                ..Deflection::default()
            },
            T,
        )
        .unwrap();

        assert!(coarse.length() < exact);
        assert!(fine.length() < exact);
        assert!(fine.length() > coarse.length());
        assert_relative_eq!(fine.length(), exact, max_relative = 1e-3);
    }

    #[test]
    fn the_angular_tolerance_catches_what_the_chord_one_misses() {
        // A large circle with a loose chord tolerance: the chord error over a
        // whole quadrant is huge in absolute terms but the *ratio* to the radius
        // is what a viewer sees, and the tangent turn is what bounds it.
        let curve = circle(1000.0);
        let chord_only = Deflection {
            chord: 50.0,
            angular: 10.0,
            ..Deflection::default()
        };
        let with_angle = Deflection {
            chord: 50.0,
            angular: 0.1,
            ..Deflection::default()
        };

        let loose = discretize(&curve, curve.domain(), chord_only, T).unwrap();
        let tight = discretize(&curve, curve.domain(), with_angle, T).unwrap();
        assert!(
            tight.segment_count() > loose.segment_count(),
            "the angular limit did nothing: {} vs {}",
            tight.segment_count(),
            loose.segment_count()
        );

        // Every segment really does turn by less than the limit.
        for window in tight.parameters.windows(2) {
            let a = curve.tangent_at(window[0], T).unwrap();
            let b = curve.tangent_at(window[1], T).unwrap();
            assert!(a.angle(b) <= 0.1 + 1e-9);
        }
    }

    #[test]
    fn a_closed_curve_gets_enough_segments_to_enclose_something() {
        // With a minimum of one, a full circle would come out as a single chord
        // from a point back to itself: zero length, zero area, and no error
        // reported anywhere.
        let curve = circle(5.0);
        let line = discretize(
            &curve,
            curve.domain(),
            Deflection {
                chord: 1e6,
                angular: 1e6,
                ..Deflection::default()
            },
            T,
        )
        .unwrap();
        assert!(line.segment_count() >= 2);
        assert!(line.length() > 0.0);
        assert!(line.is_closed(T));
    }

    #[test]
    fn reaching_the_ceiling_is_reported_rather_than_passed_off_as_success() {
        // A downstream tolerance claim built on a polyline that never met its
        // own would be untrue, with nothing to show why.
        let curve = circle(10.0);
        let line = discretize(
            &curve,
            curve.domain(),
            Deflection {
                chord: 1e-12,
                angular: 1e-12,
                min_segments: 2,
                max_segments: 16,
            },
            T,
        )
        .unwrap();
        assert!(!line.deflection_met);
        assert!(line.segment_count() <= 20);
    }

    #[test]
    fn parameters_are_kept_alongside_the_points() {
        // Two faces meeting along an edge must sample its pcurves at exactly
        // these values, or the join is not watertight. Discarding them would
        // make that impossible to arrange.
        let curve = circle(3.0);
        let line = discretize(&curve, curve.domain(), Deflection::default(), T).unwrap();
        assert_eq!(line.points.len(), line.parameters.len());
        for (u, p) in line.parameters.iter().zip(&line.points) {
            assert!(curve.point_at(*u, T).unwrap().is_equal(*p, T));
        }
        // And they increase along the curve.
        assert!(line.parameters.windows(2).all(|w| w[1] > w[0]));
    }

    #[test]
    fn a_spline_is_refined_where_it_curves_and_not_where_it_does_not() {
        // The reason for adaptive rather than uniform sampling: a curve that is
        // straight for half its length and tight for the rest should not pay
        // for the tight part everywhere.
        let control = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(10.0, 0.0, 0.0),
            Point::new(20.0, 0.0, 0.0),
            Point::new(21.0, 8.0, 0.0),
            Point::new(22.0, 0.0, 0.0),
        ];
        let curve: Curve = BSplineCurve::new(
            KnotVector::clamped_uniform(3, control.len()).unwrap(),
            control,
            T,
        )
        .unwrap()
        .into();

        let line = discretize(
            &curve,
            curve.domain(),
            Deflection {
                chord: 0.05,
                ..Deflection::default()
            },
            T,
        )
        .unwrap();
        assert!(line.deflection_met);

        // Segments in the straight first half are longer than in the curved
        // second half.
        let mid = line.points.len() / 2;
        let mean = |points: &[Point]| {
            let gaps: Vec<f64> = points.windows(2).map(|w| w[0].distance(w[1])).collect();
            #[allow(clippy::cast_precision_loss)]
            let count = gaps.len() as f64;
            gaps.iter().sum::<f64>() / count
        };
        let (early, late) = (mean(&line.points[..mid]), mean(&line.points[mid..]));
        assert!(early > late, "uniform spacing: {early} vs {late}");
    }

    #[test]
    fn discretizing_part_of_a_curve_covers_only_that_part() {
        let curve = circle(4.0);
        let line = discretize(&curve, (1.0, 2.0), Deflection::default(), T).unwrap();
        assert_relative_eq!(line.parameters[0], 1.0);
        assert_relative_eq!(line.parameters[line.parameters.len() - 1], 2.0);
        assert!(!line.is_closed(T));
        assert_relative_eq!(line.length(), 4.0, max_relative = 1e-2);
    }

    #[test]
    fn unusable_settings_are_refused() {
        let curve = circle(1.0);
        let bad = [
            Deflection {
                chord: 0.0,
                ..Deflection::default()
            },
            Deflection {
                chord: f64::NAN,
                ..Deflection::default()
            },
            Deflection {
                angular: -1.0,
                ..Deflection::default()
            },
            Deflection {
                min_segments: 0,
                ..Deflection::default()
            },
            Deflection {
                min_segments: 10,
                max_segments: 5,
                ..Deflection::default()
            },
        ];
        for deflection in bad {
            assert!(
                discretize(&curve, curve.domain(), deflection, T).is_err(),
                "accepted {deflection:?}"
            );
        }
        assert!(discretize(&curve, (1.0, 1.0), Deflection::default(), T).is_err());
        assert!(Deflection::with_chord(-1.0).is_err());
        assert!(Deflection::relative(0.0, 0.001).is_err());
        assert!(Deflection::relative(100.0, 0.001).is_ok());
    }

    #[test]
    fn a_relative_deflection_scales_with_the_model() {
        // "A thousandth of the part" survives the part being modelled in metres
        // rather than millimetres; an absolute default does not.
        let small = Deflection::relative(1.0, 0.001).unwrap();
        let large = Deflection::relative(1000.0, 0.001).unwrap();
        assert_relative_eq!(large.chord, small.chord * 1000.0);
    }

    #[test]
    fn a_planar_curve_discretizes_in_parameter_space() {
        let curve: ogeom_geom::PlanarCurve = ogeom_geom::Circle2d::new(
            ogeom_math::Circle2::centred(ogeom_math::Point2::ORIGIN, 5.0, T).unwrap(),
        )
        .into();
        let (points, parameters) = discretize_planar(
            &curve,
            (0.0, core::f64::consts::TAU),
            Deflection {
                chord: 0.05,
                ..Deflection::default()
            },
            T,
        )
        .unwrap();
        assert_eq!(points.len(), parameters.len());
        assert!(points.len() > 20, "only {} points", points.len());
        for p in &points {
            assert_relative_eq!(p.to_vector().magnitude(), 5.0, epsilon = 1e-12);
        }
    }
}
#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod on_surface_tests {
    use super::*;
    use ogeom_core::Tolerances;
    use ogeom_geom::{Circle2d, PlaneSurface, Surface as _};
    use ogeom_math::{Circle2, Frame, Plane, Point2};

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn the_spatial_chord_scales_with_the_surface_not_the_chart() {
        // Two circles in a plane's chart, radii 5 and 100: one chord
        // tolerance, measured in space, refines the big one further and
        // holds both to the same sagitta.
        let plane: ogeom_geom::SurfaceGeometry =
            PlaneSurface::over(Plane::new(Frame::WORLD), (-200.0, 200.0), (-200.0, 200.0))
                .unwrap()
                .into();
        let deflection = Deflection {
            chord: 0.05,
            ..Deflection::default()
        };
        let counts: Vec<usize> = [5.0, 100.0]
            .iter()
            .map(|&radius| {
                let circle: ogeom_geom::PlanarCurve =
                    Circle2d::new(Circle2::centred(Point2::new(0.0, 0.0), radius, T).unwrap())
                        .into();
                let (points, parameters) = discretize_on_surface(
                    &circle,
                    (0.0, core::f64::consts::TAU),
                    &plane,
                    deflection,
                    T,
                )
                .unwrap();
                // Every lifted segment's sagitta is within the chord,
                // measured at the segment's own parameter midpoint.
                for (pair, params) in points.windows(2).zip(parameters.windows(2)) {
                    use ogeom_geom::Curve2d as _;
                    let a = plane.point_at(pair[0].x, pair[0].y, T).unwrap();
                    let b = plane.point_at(pair[1].x, pair[1].y, T).unwrap();
                    let mid = circle
                        .point_at(f64::midpoint(params[0], params[1]), T)
                        .unwrap();
                    let on = plane.point_at(mid.x, mid.y, T).unwrap();
                    let chord = b - a;
                    let sagitta = (on - a).cross(chord).magnitude() / chord.magnitude();
                    assert!(
                        sagitta <= deflection.chord * 1.5,
                        "sagitta {sagitta} at radius {radius}"
                    );
                }
                points.len()
            })
            .collect();
        assert!(
            counts[1] > counts[0] * 2,
            "the larger circle refines further: {counts:?}"
        );
    }
}
