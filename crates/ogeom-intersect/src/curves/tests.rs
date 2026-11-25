use super::general::foot_via_samples;
use super::*;
use ogeom_geom::{BSpline2d, Circle2d, CircleCurve, Line2d, LineCurve};
use ogeom_math::{Circle, Circle2, Direction2, Frame, Frame2, KnotVector, Vector2};

const T: Tolerances = Tolerances::millimetres();

fn line2(from: Point2, to: Point2) -> PlanarCurve {
    Line2d::segment(from, to, T).unwrap().into()
}

fn circle2(centre: Point2, radius: f64) -> PlanarCurve {
    Circle2d::new(
        Circle2::new(
            Frame2::new(centre, Direction2::new(Vector2::new(1.0, 0.0), T).unwrap()),
            radius,
            T,
        )
        .unwrap(),
    )
    .into()
}

#[test]
fn two_lines_cross_where_algebra_says() {
    let a = line2(Point2::new(0.0, 0.0), Point2::new(4.0, 4.0));
    let b = line2(Point2::new(0.0, 4.0), Point2::new(4.0, 0.0));
    let found = intersect_curves_2d(&a, &b, CurveCurveOptions::default(), T).unwrap();
    assert_eq!(found.crossings.len(), 1);
    let hit = &found.crossings[0];
    assert!(hit.point.is_equal(Point2::new(2.0, 2.0), T));
    // Parameters are arc length on a segment.
    approx::assert_relative_eq!(hit.on_a, 8.0_f64.sqrt(), epsilon = 1e-9);

    // Segments that would cross beyond their ends do not.
    let short = line2(Point2::new(0.0, 4.0), Point2::new(1.0, 3.0));
    assert!(
        intersect_curves_2d(&a, &short, CurveCurveOptions::default(), T)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn collinear_lines_overlap_rather_than_crossing_everywhere() {
    let a = line2(Point2::new(0.0, 0.0), Point2::new(10.0, 0.0));
    let b = line2(Point2::new(4.0, 0.0), Point2::new(20.0, 0.0));
    let found = intersect_curves_2d(&a, &b, CurveCurveOptions::default(), T).unwrap();
    assert!(found.crossings.is_empty());
    assert_eq!(found.overlaps.len(), 1);
    let overlap = &found.overlaps[0];
    approx::assert_relative_eq!(overlap.on_a.0, 4.0, epsilon = 1e-9);
    approx::assert_relative_eq!(overlap.on_a.1, 10.0, epsilon = 1e-9);
    approx::assert_relative_eq!(overlap.on_b.0, 0.0, epsilon = 1e-9);
    approx::assert_relative_eq!(overlap.on_b.1, 6.0, epsilon = 1e-9);

    // Parallel but apart: nothing.
    let above = line2(Point2::new(0.0, 1.0), Point2::new(10.0, 1.0));
    assert!(
        intersect_curves_2d(&a, &above, CurveCurveOptions::default(), T)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn a_line_meets_a_circle_in_two_points_one_or_none() {
    let circle = circle2(Point2::new(0.0, 0.0), 2.0);
    let through = line2(Point2::new(-5.0, 0.0), Point2::new(5.0, 0.0));
    let found = intersect_curves_2d(&through, &circle, CurveCurveOptions::default(), T).unwrap();
    assert_eq!(found.crossings.len(), 2);
    for hit in &found.crossings {
        approx::assert_relative_eq!(
            hit.point.distance(Point2::new(0.0, 0.0)),
            2.0,
            epsilon = 1e-9
        );
        // The circle parameter really evaluates to the crossing point.
        let PlanarCurve::Circle(_) = &circle else {
            unreachable!()
        };
        let on_circle = circle.point_at(hit.on_b, T).unwrap();
        assert!(on_circle.is_equal(hit.point, T));
    }

    let tangent = line2(Point2::new(-5.0, 2.0), Point2::new(5.0, 2.0));
    assert_eq!(
        intersect_curves_2d(&tangent, &circle, CurveCurveOptions::default(), T)
            .unwrap()
            .crossings
            .len(),
        1
    );
    let missing = line2(Point2::new(-5.0, 3.0), Point2::new(5.0, 3.0));
    assert!(
        intersect_curves_2d(&missing, &circle, CurveCurveOptions::default(), T)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn two_circles_cross_touch_coincide_or_miss() {
    let a = circle2(Point2::new(0.0, 0.0), 2.0);

    let crossing = circle2(Point2::new(3.0, 0.0), 2.0);
    let found = intersect_curves_2d(&a, &crossing, CurveCurveOptions::default(), T).unwrap();
    assert_eq!(found.crossings.len(), 2);
    for hit in &found.crossings {
        let on_a = a.point_at(hit.on_a, T).unwrap();
        let on_b = crossing.point_at(hit.on_b, T).unwrap();
        assert!(on_a.is_equal(hit.point, T));
        assert!(on_b.is_equal(hit.point, T));
    }

    let touching = circle2(Point2::new(4.0, 0.0), 2.0);
    assert_eq!(
        intersect_curves_2d(&a, &touching, CurveCurveOptions::default(), T)
            .unwrap()
            .crossings
            .len(),
        1
    );

    let same = circle2(Point2::new(0.0, 0.0), 2.0);
    let coincident = intersect_curves_2d(&a, &same, CurveCurveOptions::default(), T).unwrap();
    assert!(coincident.crossings.is_empty());
    assert_eq!(coincident.overlaps.len(), 1);

    let apart = circle2(Point2::new(10.0, 0.0), 2.0);
    assert!(
        intersect_curves_2d(&a, &apart, CurveCurveOptions::default(), T)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn the_general_path_handles_what_has_no_closed_form() {
    // A spline sine-ish wave against a line: three crossings, found by
    // sampling and polished by Newton to rounding.
    let wave: PlanarCurve = BSpline2d::new(
        KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0], 3).unwrap(),
        vec![
            Point2::new(0.0, -1.0),
            Point2::new(1.0, 3.0),
            Point2::new(2.0, -3.0),
            Point2::new(3.0, 3.0),
            Point2::new(4.0, -1.0),
        ],
        T,
    )
    .unwrap()
    .into();
    let axis = line2(Point2::new(-1.0, 0.0), Point2::new(5.0, 0.0));
    let found = intersect_curves_2d(&wave, &axis, CurveCurveOptions::default(), T).unwrap();
    assert_eq!(found.crossings.len(), 3, "a wave crosses its axis thrice");
    for hit in &found.crossings {
        assert!(hit.gap < 1e-9);
        assert!(hit.point.y.abs() < 1e-9);
        let on_wave = wave.point_at(hit.on_a, T).unwrap();
        assert!(on_wave.is_equal(hit.point, T));
    }
}

/// One circle written twice in the plane, a quarter turn apart and
/// wound against each other: the overlap's ranges correspond.
#[test]
fn one_planar_circle_written_twice_states_the_correspondence() {
    let a = circle2(Point2::new(1.0, 2.0), 3.0);
    let quarter = Frame2::new(
        Point2::new(1.0, 2.0),
        Direction2::new(Vector2::new(0.0, 1.0), T).unwrap(),
    );
    let b: PlanarCurve = Circle2d::new(Circle2::new(quarter, 3.0, T).unwrap()).into();
    let flipped: PlanarCurve = ogeom_geom::Reversible::reversed(&b);
    for other in [b, flipped] {
        let found = intersect_curves_2d(&a, &other, CurveCurveOptions::default(), T).unwrap();
        assert_eq!(found.overlaps.len(), 1);
        let overlap = found.overlaps[0];
        for i in 0..=8 {
            let t = f64::from(i) / 8.0;
            let ta = (overlap.on_a.1 - overlap.on_a.0).mul_add(t, overlap.on_a.0);
            let tb = (overlap.on_b.1 - overlap.on_b.0).mul_add(t, overlap.on_b.0);
            let pa = a.point_at(ta, T).unwrap();
            let pb = other
                .point_at(tb.rem_euclid(core::f64::consts::TAU), T)
                .unwrap();
            assert!(pa.distance(pb) < 1e-9, "at {t}: {pa:?} against {pb:?}");
        }
    }
}

/// Collinear segments meeting end to end share their end point, in the
/// plane and in space, as a perpendicular pair meeting there does.
#[test]
fn collinear_segments_meeting_end_to_end_share_a_point() {
    let a = line2(Point2::new(0.0, 0.0), Point2::new(5.0, 0.0));
    let b = line2(Point2::new(5.0, 0.0), Point2::new(9.0, 0.0));
    let found = intersect_curves_2d(&a, &b, CurveCurveOptions::default(), T).unwrap();
    assert!(found.overlaps.is_empty());
    assert_eq!(found.crossings.len(), 1);
    assert!(found.crossings[0].point.distance(Point2::new(5.0, 0.0)) < 1e-12);

    let a: Curve = LineCurve::segment(Point::ORIGIN, Point::new(0.0, 0.0, 5.0), T)
        .unwrap()
        .into();
    let b: Curve = LineCurve::segment(Point::new(0.0, 0.0, 5.0), Point::new(0.0, 0.0, 7.0), T)
        .unwrap()
        .into();
    let found = intersect_curves(&a, &b, CurveCurveOptions::default(), T).unwrap();
    assert!(found.overlaps.is_empty());
    assert_eq!(found.crossings.len(), 1);
    assert!(found.crossings[0].point.distance(Point::new(0.0, 0.0, 5.0)) < 1e-12);
}

/// Two descriptions of one circle overlap over the whole turn, and the
/// overlap's two ranges *correspond*: `on_b`'s ends are where `b` stands
/// at `a`'s own ends. A caller carrying a split from one to the other
/// (the boolean, pairing a hole's arcs against the disc that fills them)
/// reads that correspondence and gets the same point back, whatever phase
/// and winding the two were written with.
#[test]
fn one_circle_written_twice_states_the_correspondence_between_them() {
    use ogeom_math::{Direction, Vector};
    let a: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 3.0, T).unwrap()).into();
    // The same circle seen from underneath, started a third of a turn
    // round: opposite winding, arbitrary phase.
    let third = 2.0 * core::f64::consts::PI / 3.0;
    let flipped = Frame::new(
        Point::ORIGIN,
        -Direction::Z,
        Direction::new(Vector::new(third.cos(), third.sin(), 0.0), T).unwrap(),
        T,
    )
    .unwrap();
    let b: Curve = CircleCurve::new(Circle::new(flipped, 3.0, T).unwrap()).into();
    for pair in [(&a, &b), (&b, &a)] {
        let found = intersect_curves(pair.0, pair.1, CurveCurveOptions::default(), T).unwrap();
        assert!(
            found.crossings.is_empty(),
            "every point is a hit, so none is"
        );
        assert_eq!(found.overlaps.len(), 1);
        let overlap = &found.overlaps[0];
        let span = overlap.on_a.1 - overlap.on_a.0;
        for i in 0..=8 {
            let t = f64::from(i) / 8.0;
            let ta = span.mul_add(t, overlap.on_a.0);
            let tb = (overlap.on_b.1 - overlap.on_b.0).mul_add(t, overlap.on_b.0);
            let pa = pair.0.point_at(ta, T).unwrap();
            let pb = pair
                .1
                .point_at(tb.rem_euclid(core::f64::consts::TAU), T)
                .unwrap();
            assert!(
                pa.distance(pb) < 1e-9,
                "at {t}: {pa:?} against {pb:?} (on_a {:?} on_b {:?})",
                overlap.on_a,
                overlap.on_b
            );
        }
    }
}

#[test]
fn space_curves_cross_within_a_gap_and_report_it() {
    // Two circles that would cross in a shared plane, with one lifted a
    // hair out of it: the crossings become passes with a real, small gap
    // that must be reported, not zeroed. (Not chain links: linked circles
    // never approach, since passing through each other's *disks* is what
    // linked means, and these radii hold the curves a constant two units
    // apart.)
    let a: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 2.0, T).unwrap()).into();
    let lifted = Frame::new(
        Point::new(3.0, 0.0, 0.001),
        ogeom_math::Direction::Z,
        ogeom_math::Direction::X,
        T,
    )
    .unwrap();
    let b: Curve = CircleCurve::new(Circle::new(lifted, 2.0, T).unwrap()).into();

    let options = CurveCurveOptions {
        gap: 1e-2,
        ..CurveCurveOptions::default()
    };
    let found = intersect_curves(&a, &b, options, T).unwrap();
    assert_eq!(found.crossings.len(), 2, "two near-crossings");
    for hit in &found.crossings {
        assert!(hit.gap > 1e-4, "the gap is real and must not be zeroed");
        assert!(hit.gap < 2e-3, "but small: {}", hit.gap);
    }

    // Tighten the gap below the offset and the crossings vanish.
    let strict = CurveCurveOptions {
        gap: 1e-5,
        ..CurveCurveOptions::default()
    };
    assert!(intersect_curves(&a, &b, strict, T).unwrap().is_empty());
}

/// A figure eight passes its double point twice, and a line through
/// that point is crossed once per passage, at parameters half a turn
/// apart. The ring of controls is symmetric under a mirror that shifts
/// the parameter by half a turn and under a point reflection, so the
/// double point is the origin.
#[test]
fn a_line_through_a_figure_eight_s_double_point_is_crossed_twice() {
    let ring = [
        Point::new(2.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
        Point::new(-1.0, -1.0, 0.0),
        Point::new(-2.0, 0.0, 0.0),
        Point::new(-1.0, 1.0, 0.0),
        Point::new(1.0, -1.0, 0.0),
    ];
    let a: Curve = ogeom_geom::BSplineCurve::periodic(&ring, 3, T)
        .unwrap()
        .into();
    let b: Curve = LineCurve::segment(Point::new(0.0, 0.0, -1.0), Point::new(0.0, 0.0, 1.0), T)
        .unwrap()
        .into();
    let found = intersect_curves(&a, &b, CurveCurveOptions::default(), T).unwrap();
    assert_eq!(found.crossings.len(), 2, "{:?}", found.crossings);
    let (lo, hi) = a.domain();
    let apart = (found.crossings[0].on_a - found.crossings[1].on_a).abs();
    assert!(
        (apart - (hi - lo) / 2.0).abs() < 1e-6,
        "half a turn apart: {apart}"
    );
    for hit in &found.crossings {
        assert!(hit.point.distance(Point::ORIGIN) < 1e-9, "{:?}", hit.point);
    }
}

#[test]
fn skew_lines_in_space_miss_and_close_ones_meet() {
    let a: Curve = LineCurve::segment(Point::ORIGIN, Point::new(10.0, 0.0, 0.0), T)
        .unwrap()
        .into();
    let skew: Curve = LineCurve::segment(Point::new(0.0, -5.0, 1.0), Point::new(0.0, 5.0, 1.0), T)
        .unwrap()
        .into();
    assert!(
        intersect_curves(&a, &skew, CurveCurveOptions::default(), T)
            .unwrap()
            .is_empty(),
        "a unit apart is not a crossing"
    );

    let meeting: Curve =
        LineCurve::segment(Point::new(5.0, -5.0, 0.0), Point::new(5.0, 5.0, 0.0), T)
            .unwrap()
            .into();
    let found = intersect_curves(&a, &meeting, CurveCurveOptions::default(), T).unwrap();
    assert_eq!(found.crossings.len(), 1);
    assert!(
        found.crossings[0]
            .point
            .is_equal(Point::new(5.0, 0.0, 0.0), T)
    );
    assert!(found.crossings[0].gap < 1e-12);

    // Collinear 3D lines overlap.
    let collinear: Curve =
        LineCurve::segment(Point::new(4.0, 0.0, 0.0), Point::new(20.0, 0.0, 0.0), T)
            .unwrap()
            .into();
    let shared = intersect_curves(&a, &collinear, CurveCurveOptions::default(), T).unwrap();
    assert_eq!(shared.overlaps.len(), 1);
}

#[test]
fn a_fitted_curve_tracing_an_arc_is_one_overlap_not_a_row_of_crossings() {
    // A spline fitted along a circle's arc sits within its fit budget
    // of the circle everywhere, and "crosses" it at every wobble. The
    // sampling path reports the stretch as one overlap and keeps no
    // crossing inside it. Read as crossings, a section tracing the arc
    // it was cut along would shatter into hundreds of pieces.
    let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 4.0, T).unwrap()).into();
    let points: Vec<Point> = (0..=40)
        .map(|i| {
            let a = 0.2 + 1.0 * f64::from(i) / 40.0;
            Point::new(4.0 * a.cos(), 4.0 * a.sin(), 0.0)
        })
        .collect();
    let fitted: Curve = ogeom_geom::fit::fit_points(&points, 3, 1e-7, T)
        .unwrap()
        .curve
        .into();
    let options = CurveCurveOptions {
        gap: 1e-5,
        ..CurveCurveOptions::default()
    };
    let found = intersect_curves(&fitted, &circle, options, T).unwrap();
    assert_eq!(found.overlaps.len(), 1, "one shared stretch: {found:?}");
    let (lo, hi) = found.overlaps[0].on_a;
    let (fa, fb) = fitted.domain();
    assert!(
        lo - fa < 1e-3 && fb - hi < 1e-3,
        "the whole fit runs along the circle"
    );
    assert!(
        found.crossings.is_empty(),
        "no crossing survives inside the overlap: {:?}",
        found.crossings
    );
}

/// An ellipse's arc from `5.5` round past its seam at `2π` to `7`,
/// spelled as a rational B-spline, runs along the ellipse over that
/// arc: one overlap, read on the ellipse from `5.5` to `7`, the way the
/// arc goes round. Read the short way between the seam's two sides it
/// named the rest of the ellipse, `0.72` to `5.5`.
#[test]
fn an_arc_across_a_periodic_curve_s_seam_overlaps_it_the_way_it_goes() {
    let ellipse: Curve = ogeom_geom::EllipseCurve::new(
        ogeom_math::Ellipse::new(Frame::WORLD, 0.5025, 0.5, T).unwrap(),
    )
    .into();
    let arc: Curve = ellipse.to_bspline_over((5.5, 7.0), T).unwrap().into();
    let found = intersect_curves(&arc, &ellipse, CurveCurveOptions::default(), T).unwrap();
    let [overlap] = found.overlaps.as_slice() else {
        panic!("one shared stretch: {found:?}");
    };
    let (lo, hi) = overlap.on_a;
    assert!(lo.abs() < 1e-9 && (hi - 1.0).abs() < 1e-9);
    let (from, to) = overlap.on_b;
    assert!(
        (from - 5.5).abs() < 1e-9 && (to - 7.0).abs() < 1e-9,
        "{overlap:?}"
    );
}

#[test]
fn an_arc_ending_tangent_to_a_line_is_one_crossing_with_its_reach() {
    // A circle and its tangent line touch at one point, but the
    // stationarity conditions go flat along the touch and every seed
    // converges somewhere in a valley the width of the gap. One contact
    // comes back, at the touch, owning the valley's length as its reach.
    let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 4.0, T).unwrap()).into();
    let line: Curve = LineCurve::segment(Point::new(4.0, -3.0, 0.0), Point::new(4.0, 3.0, 0.0), T)
        .unwrap()
        .into();
    let options = CurveCurveOptions {
        gap: 1e-5,
        ..CurveCurveOptions::default()
    };
    let found = intersect_curves(&circle, &line, options, T).unwrap();
    assert_eq!(found.crossings.len(), 1, "one touch: {:?}", found.crossings);
    let touch = found.crossings[0];
    assert!(
        touch.point.distance(Point::new(4.0, 0.0, 0.0)) < 2e-2,
        "{touch:?}"
    );
    assert!(touch.reach < 5e-2, "the valley is short: {touch:?}");
    assert!(found.overlaps.is_empty());
}

/// Lines touching circles and ellipses in their plane, whole and
/// trimmed, in either order, and passing beside a trim that leaves the
/// touch out: the closed form answers what the sampling path does, the
/// same crossings at the same parameters and gaps owning the same
/// reach. A line a hair beside the touch, missing or crossing, is no
/// double root and stays the sampling path's.
#[test]
fn a_tangent_line_meets_a_conic_as_the_sampling_path_says() {
    use ogeom_geom::{EllipseCurve, TrimmedCurve};
    use ogeom_math::{Direction, Ellipse, Vector};
    use std::f64::consts::{FRAC_PI_2, PI};
    let tilted = Frame::new(
        Point::new(3.0, -2.0, 7.0),
        Direction::new(Vector::new(1.0, 2.0, 2.0), T).unwrap(),
        Direction::new(Vector::new(2.0, -1.0, 0.0), T).unwrap(),
        T,
    )
    .unwrap();
    let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 4.0, T).unwrap()).into();
    let tilted_circle: Curve = CircleCurve::new(Circle::new(tilted, 1.0, T).unwrap()).into();
    let ellipse: Curve = EllipseCurve::new(Ellipse::new(Frame::WORLD, 3.0, 2.0, T).unwrap()).into();
    let arc = |from: f64, to: f64| -> Curve {
        TrimmedCurve::new(circle.clone(), from, to, T)
            .unwrap()
            .into()
    };
    let segment =
        |from: Point, to: Point| -> Curve { LineCurve::segment(from, to, T).unwrap().into() };
    let (ox, oy) = (tilted.x().vector(), tilted.y().vector());
    let (eu, ev) = (3.0 * 2.5_f64.cos(), 2.0 * 2.5_f64.sin());
    let ellipse_tangent = Vector::new(-3.0 * 2.5_f64.sin(), 2.0 * 2.5_f64.cos(), 0.0);
    let along_top = segment(Point::new(-5.0, 4.0, 0.0), Point::new(1.0, 4.0, 0.0));
    let cases: Vec<(&str, Curve, Curve, bool)> = vec![
        ("circle", circle.clone(), along_top.clone(), true),
        (
            "tilted",
            tilted_circle,
            segment(
                tilted.origin() + oy - ox * 2.0,
                tilted.origin() + oy + ox * 5.0,
            ),
            true,
        ),
        (
            "ellipse's end",
            ellipse.clone(),
            segment(Point::new(-3.0, -1.0, 0.0), Point::new(-3.0, 4.0, 0.0)),
            true,
        ),
        (
            "ellipse's side",
            ellipse.clone(),
            segment(Point::new(-6.0, -2.0, 0.0), Point::new(2.0, -2.0, 0.0)),
            true,
        ),
        (
            "ellipse aslant",
            ellipse,
            segment(
                Point::new(eu, ev, 0.0) - ellipse_tangent,
                Point::new(eu, ev, 0.0) + ellipse_tangent * 2.0,
            ),
            true,
        ),
        (
            "arc ending at the touch",
            arc(0.0, FRAC_PI_2),
            along_top.clone(),
            true,
        ),
        (
            "arc starting at the touch",
            arc(FRAC_PI_2, PI),
            along_top.clone(),
            true,
        ),
        (
            "line ending at the touch",
            circle.clone(),
            segment(Point::new(-3.0, 4.0, 0.0), Point::new(0.0, 4.0, 0.0)),
            true,
        ),
        (
            "lifted",
            circle.clone(),
            segment(Point::new(-5.0, 4.0, 2e-8), Point::new(1.0, 4.0, 2e-8)),
            true,
        ),
        (
            "arc away from the touch",
            arc(PI, 2.0 * PI),
            along_top,
            true,
        ),
        (
            "line short of the touch",
            circle.clone(),
            segment(Point::new(1.0, 4.0, 0.0), Point::new(3.0, 4.0, 0.0)),
            true,
        ),
        (
            "missing by a hair",
            circle.clone(),
            segment(
                Point::new(-5.0, 4.0 + 1e-9, 0.0),
                Point::new(1.0, 4.0 + 1e-9, 0.0),
            ),
            false,
        ),
        (
            "crossing a hair inside",
            circle,
            segment(
                Point::new(-5.0, 4.0 - 1e-9, 0.0),
                Point::new(1.0, 4.0 - 1e-9, 0.0),
            ),
            false,
        ),
    ];
    for gap in [1e-7, 1e-5] {
        let options = CurveCurveOptions {
            gap,
            ..CurveCurveOptions::default()
        };
        for (name, conic, line, touches) in &cases {
            for (a, b) in [(conic, line), (line, conic)] {
                let closed = tangent_line_conic_3d(a, b, options, T);
                let Some(closed) = closed else {
                    assert!(!touches, "{name}, gap {gap}: no closed form");
                    continue;
                };
                assert!(touches, "{name}, gap {gap}: {closed:?}");
                let sampled = general_3d(a, b, options, T).unwrap();
                let context = format!("{name}, gap {gap}: {closed:?} against {sampled:?}");
                assert!(closed.overlaps.is_empty() && sampled.overlaps.is_empty());
                assert_eq!(closed.crossings.len(), sampled.crossings.len(), "{context}");
                // The polish lands anywhere in the valley, and the
                // sampling path measures the valley from where it
                // landed: the two crossings stand within its reach, and
                // their reaches differ by no more than they stand apart.
                for (x, y) in closed.crossings.iter().zip(&sampled.crossings) {
                    let apart = x.point.distance(y.point);
                    let (p, q) = (
                        b.point_at(x.on_b, T).unwrap(),
                        a.point_at(x.on_a, T).unwrap(),
                    );
                    assert!(p.distance(q) <= x.gap + 1e-15, "{context}");
                    assert!(apart <= y.reach.max(gap), "{context}");
                    assert!(x.gap <= y.gap + 1e-15, "{context}");
                    assert_eq!(x.reach > 0.0, y.reach > 0.0, "{context}");
                    assert!(
                        (x.reach - y.reach).abs() <= apart + 1e-6 * y.reach,
                        "{context}"
                    );
                }
            }
        }
    }
}

#[test]
fn unusable_options_are_refused() {
    let a = line2(Point2::new(0.0, 0.0), Point2::new(1.0, 0.0));
    for options in [
        CurveCurveOptions {
            samples: 1,
            ..CurveCurveOptions::default()
        },
        CurveCurveOptions {
            gap: 0.0,
            ..CurveCurveOptions::default()
        },
        CurveCurveOptions {
            gap: f64::NAN,
            ..CurveCurveOptions::default()
        },
    ] {
        assert!(intersect_curves_2d(&a, &a.clone(), options, T).is_err());
    }
}

/// Plane sections of one cylinder: a circle and two ellipses tilted
/// different ways. Any two meet where their planes' common line pierces
/// the cylinder, twice, and a circle lifted parallel to another meets
/// it nowhere.
#[test]
fn circles_and_ellipses_in_different_planes_meet_where_the_planes_do() {
    use ogeom_geom::EllipseCurve;
    use ogeom_math::{Direction, Ellipse, Vector};
    let radius = 2.0;
    let tilted = |normal: Vector, major: Vector, lean: f64| -> Curve {
        let frame = Frame::new(
            Point::ORIGIN,
            Direction::new(normal, T).unwrap(),
            Direction::new(major, T).unwrap(),
            T,
        )
        .unwrap();
        EllipseCurve::new(Ellipse::new(frame, radius / lean.cos(), radius, T).unwrap()).into()
    };
    let (p, q) = (0.4_f64, 0.7_f64);
    let about_x = tilted(
        Vector::new(0.0, -p.sin(), p.cos()),
        Vector::new(0.0, p.cos(), p.sin()),
        p,
    );
    let about_y = tilted(
        Vector::new(-q.sin(), 0.0, q.cos()),
        Vector::new(q.cos(), 0.0, q.sin()),
        q,
    );
    let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, radius, T).unwrap()).into();
    for (a, b) in [
        (&about_x, &about_y),
        (&about_y, &about_x),
        (&circle, &about_x),
        (&about_y, &circle),
    ] {
        let found = intersect_curves(a, b, CurveCurveOptions::default(), T).unwrap();
        assert_eq!(found.crossings.len(), 2, "{found:?}");
        for hit in &found.crossings {
            let on_a = a.point_at(hit.on_a, T).unwrap();
            let on_b = b.point_at(hit.on_b, T).unwrap();
            assert!(on_a.distance(on_b) < 1e-9, "{on_a:?} against {on_b:?}");
            assert!((on_a.x.hypot(on_a.y) - radius).abs() < 1e-9);
        }
    }
    let lifted = Frame::new(Point::new(0.0, 0.0, 1.0), Direction::Z, Direction::X, T).unwrap();
    let above: Curve = CircleCurve::new(Circle::new(lifted, radius, T).unwrap()).into();
    assert!(
        intersect_curves(&circle, &above, CurveCurveOptions::default(), T)
            .unwrap()
            .is_empty()
    );
}

/// Two circles in planes a millionth of a radian apart, passing within
/// the gap of each other where their shadows cross: the first crosses
/// the second's plane a twentieth of a unit away from there, where the
/// two are far apart. The near pass is found all the same.
#[test]
fn circles_in_all_but_one_plane_meet_where_they_pass() {
    use ogeom_math::{Direction, Vector};
    let flat: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 2.0, T).unwrap()).into();
    let lean = 1e-6;
    let normal = Direction::new(Vector::new(lean, 0.0, 1.0), T).unwrap();
    let centre = Point::new(1.0, 0.0, lean.mul_add(-0.5, 5e-8));
    let tilted_frame = Frame::new(
        centre,
        normal,
        Direction::new(Vector::new(1.0, 0.0, -lean), T).unwrap(),
        T,
    )
    .unwrap();
    let tilted: Curve = CircleCurve::new(Circle::new(tilted_frame, 2.0, T).unwrap()).into();
    let found = intersect_curves(&flat, &tilted, CurveCurveOptions::default(), T).unwrap();
    assert_eq!(found.crossings.len(), 2, "{found:?}");
    for hit in &found.crossings {
        assert!(hit.gap < 1e-7, "{hit:?}");
        let p = flat.point_at(hit.on_a, T).unwrap();
        assert!((p.x - 0.5).abs() < 1e-4, "{p:?}");
    }
}

/// A line across an ellipse in its plane meets it twice, a line through
/// a circle's plane meets it once where it pierces the circle, and a
/// line through the plane inside the circle misses it. Either order of
/// the pair answers the same, parameters swapped.
#[test]
fn lines_meet_circles_and_ellipses_in_closed_form() {
    use ogeom_geom::EllipseCurve;
    use ogeom_math::Ellipse;
    let ellipse: Curve = EllipseCurve::new(Ellipse::new(Frame::WORLD, 3.0, 2.0, T).unwrap()).into();
    let across: Curve =
        LineCurve::segment(Point::new(-5.0, 1.0, 0.0), Point::new(5.0, 1.0, 0.0), T)
            .unwrap()
            .into();
    for (a, b, line_first) in [(&ellipse, &across, false), (&across, &ellipse, true)] {
        let found = intersect_curves(a, b, CurveCurveOptions::default(), T).unwrap();
        assert_eq!(found.crossings.len(), 2, "{found:?}");
        for hit in &found.crossings {
            let p = a.point_at(hit.on_a, T).unwrap();
            let q = b.point_at(hit.on_b, T).unwrap();
            assert!(p.distance(q) < 1e-9, "{p:?} against {q:?}");
            let on_line = if line_first { p } else { q };
            assert!((on_line.y - 1.0).abs() < 1e-12);
        }
    }
    let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 2.0, T).unwrap()).into();
    let through: Curve =
        LineCurve::segment(Point::new(2.0, 0.0, -1.0), Point::new(2.0, 0.0, 1.0), T)
            .unwrap()
            .into();
    let found = intersect_curves(&circle, &through, CurveCurveOptions::default(), T).unwrap();
    assert_eq!(found.crossings.len(), 1);
    assert!(found.crossings[0].on_a.abs() < 1e-9);
    assert!((found.crossings[0].on_b - 1.0).abs() < 1e-9);
    let inside: Curve =
        LineCurve::segment(Point::new(1.0, 0.0, -1.0), Point::new(1.0, 0.0, 1.0), T)
            .unwrap()
            .into();
    assert!(
        intersect_curves(&circle, &inside, CurveCurveOptions::default(), T)
            .unwrap()
            .is_empty()
    );
}

/// A line in a circle's plane whose chord through it lies a few tenths
/// of a micron deep, past the gap, crosses it twice at a few hundredths
/// of a degree: at the roots `7 +- sqrt(depth (2 r - depth))` along the
/// line, each on both curves to rounding, in either order of the pair.
/// Read as a touch, or not at all, the face the circle bounds and the
/// face across the line keep different splits of the chord.
#[test]
fn a_chord_deeper_than_the_gap_crosses_twice_however_shallowly() {
    use ogeom_math::{Direction, Vector};
    let line: Curve =
        LineCurve::segment(Point::new(0.0, 0.0, 10.0), Point::new(0.0, 14.0, 10.0), T)
            .unwrap()
            .into();
    for depth in [3e-7, 5e-7, 1e-6] {
        for seam in [Vector::new(0.0, 0.0, 1.0), Vector::new(0.0, 1.0, 0.0)] {
            let frame = Frame::new(
                Point::new(0.0, 7.0, 8.0 + depth),
                Direction::X,
                Direction::new(seam, T).unwrap(),
                T,
            )
            .unwrap();
            let circle: Curve = CircleCurve::new(Circle::new(frame, 2.0, T).unwrap()).into();
            // The depth the centre stands at as stored, and how far
            // along the line the coordinates' rounding moves a crossing
            // at the angle `half / r` it crosses at.
            let stored = (8.0 + depth) - 8.0;
            let half = (stored * (4.0 - stored)).sqrt();
            let along = 16.0 * f64::EPSILON * 10.0 * 2.0 / half;
            for (a, b, line_first) in [(&line, &circle, true), (&circle, &line, false)] {
                let found = intersect_curves(a, b, CurveCurveOptions::default(), T).unwrap();
                assert_eq!(found.crossings.len(), 2, "depth {depth}: {found:?}");
                let mut ys: Vec<f64> = found.crossings.iter().map(|c| c.point.y).collect();
                ys.sort_by(f64::total_cmp);
                for (y, want) in ys.iter().zip([7.0 - half, 7.0 + half]) {
                    assert!(
                        (y - want).abs() < along,
                        "depth {depth}: {y} against {want}"
                    );
                }
                for hit in &found.crossings {
                    let (on_line, on_circle) = if line_first {
                        (hit.on_a, hit.on_b)
                    } else {
                        (hit.on_b, hit.on_a)
                    };
                    let p = line.point_at(on_line, T).unwrap();
                    let q = circle.point_at(on_circle, T).unwrap();
                    assert!(p.distance(q) < 1e-12, "{p:?} against {q:?}");
                }
            }
        }
    }
}

/// Two trims of one circle share whatever stretch they share, across the
/// seam or not, walked the same way or the other: every piece comes
/// back inside both windows, and each end of it is the same point on
/// both curves.
#[test]
fn arcs_of_one_circle_overlap_across_the_seam() {
    use ogeom_geom::{Curve3d as _, TrimmedCurve};
    let circle = |x: ogeom_math::Direction, z: ogeom_math::Direction| -> Curve {
        let frame = Frame::new(Point::new(10.0, 20.0, 30.0), z, x, T).unwrap();
        CircleCurve::new(Circle::new(frame, 50.0, T).unwrap()).into()
    };
    let trim = |c: &Curve, a: f64, b: f64| -> Curve {
        TrimmedCurve::new(c.clone(), a, b, T).unwrap().into()
    };
    let (x, y) = (ogeom_math::Direction::X, ogeom_math::Direction::Y);
    let (up, down) = (ogeom_math::Direction::Z, -ogeom_math::Direction::Z);
    let plain = circle(x, up);
    let turned = circle(y, up);
    let backwards = circle(x, down);
    // `a`'s share of its own length that `b` covers.
    for (a, b, share) in [
        (trim(&plain, 5.5, 7.0), trim(&plain, 0.2, 1.0), 0.517 / 1.5),
        (trim(&plain, 0.2, 1.0), trim(&plain, 5.5, 7.0), 0.517 / 0.8),
        (
            trim(&plain, 5.8, 6.8),
            trim(&backwards, 5.0, 6.2),
            0.4336 / 1.0,
        ),
        (trim(&plain, 5.0, 7.5), backwards.clone(), 1.0),
        (trim(&plain, 0.5, 6.0), trim(&turned, 4.0, 6.0), 1.218 / 5.5),
    ] {
        let found = intersect_curves(&a, &b, CurveCurveOptions::default(), T).unwrap();
        let (wa, wb) = (a.domain(), b.domain());
        let mut covered = 0.0;
        for o in &found.overlaps {
            for (t, w) in [
                (o.on_a.0, wa),
                (o.on_a.1, wa),
                (o.on_b.0, wb),
                (o.on_b.1, wb),
            ] {
                assert!(t >= w.0 - 1e-9 && t <= w.1 + 1e-9, "{t} outside {w:?}");
            }
            for (ta, tb) in [(o.on_a.0, o.on_b.0), (o.on_a.1, o.on_b.1)] {
                let gap = a
                    .point_at(ta, T)
                    .unwrap()
                    .distance(b.point_at(tb, T).unwrap());
                assert!(gap < 1e-9, "ends {gap} apart");
            }
            covered += (o.on_a.1 - o.on_a.0).abs();
        }
        let want = share * (wa.1 - wa.0);
        assert!((covered - want).abs() < 2e-3, "{covered} against {want}");
    }
}

/// A trimmed arc answers as its circle does inside its window: a line
/// tangent to it, one crossing it twice a hair below the tangent, and a
/// circle touching it from outside are all found.
#[test]
fn a_trimmed_arc_meets_tangents_as_its_circle_does() {
    use ogeom_geom::Trimmed2d;
    let options = CurveCurveOptions::default();
    let circle = circle2(Point2::new(0.0, 0.0), 10.0);
    let arc: PlanarCurve = Trimmed2d::new(circle, 0.3, 2.9, T).unwrap().into();
    let tangent = line2(Point2::new(-20.0, 10.0), Point2::new(20.0, 10.0));
    let found = intersect_curves_2d(&arc, &tangent, options, T).unwrap();
    assert_eq!(found.crossings.len(), 1);
    assert!(found.crossings[0].point.is_equal(Point2::new(0.0, 10.0), T));
    let grazing = line2(
        Point2::new(-20.0, 10.0 - 1e-4),
        Point2::new(20.0, 10.0 - 1e-4),
    );
    let found = intersect_curves_2d(&arc, &grazing, options, T).unwrap();
    assert_eq!(found.crossings.len(), 2);
    let beside: PlanarCurve = Trimmed2d::new(circle2(Point2::new(20.0, 0.0), 10.0), 2.0, 4.5, T)
        .unwrap()
        .into();
    let whole = circle2(Point2::new(0.0, 0.0), 10.0);
    let found = intersect_curves_2d(&whole, &beside, options, T).unwrap();
    assert_eq!(found.crossings.len(), 1);
    assert!(found.crossings[0].point.is_equal(Point2::new(10.0, 0.0), T));
    // Outside the window, nothing.
    let below = line2(Point2::new(-20.0, -10.0), Point2::new(20.0, -10.0));
    assert!(
        intersect_curves_2d(&arc, &below, options, T)
            .unwrap()
            .is_empty()
    );
}

/// A point's nearest segment of a sampled curve, found by the boxed
/// scan, is the one every segment measured finds, to the bit: from
/// points on the curve, beside it, at its centre and far off, asked in
/// an order that leaves the last segment found far from the next.
#[test]
fn the_boxed_nearest_segment_is_the_one_every_segment_finds() {
    let frame = Frame::new(
        Point::new(120.0, -40.0, 7.5),
        ogeom_math::Direction::new(ogeom_math::Vector::new(0.2, -0.3, 1.0), T).unwrap(),
        ogeom_math::Direction::new(ogeom_math::Vector::new(1.0, 0.0, 0.2), T).unwrap(),
        T,
    )
    .unwrap();
    let curve: Curve = CircleCurve::new(Circle::new(frame, 6.0, T).unwrap()).into();
    let sampled = sample_3d(&curve, 64, T);
    let segments = Segments::of(&sampled);
    let every = |p: Point| {
        let mut seed = (f64::INFINITY, 0.0);
        for j in 1..sampled.points.len() {
            let (_, tb, gap) =
                segments_approach_3d((p, p), (sampled.points[j - 1], sampled.points[j]));
            if gap < seed.0 {
                let (a, b) = (sampled.parameters[j - 1], sampled.parameters[j]);
                seed = (gap, a + (b - a) * tb);
            }
        }
        seed
    };
    for k in 0..97 {
        let t = f64::from(k * 37 % 97) / 97.0 * core::f64::consts::TAU;
        let on = curve.point_at(t, T).unwrap();
        for p in [
            on,
            on + ogeom_math::Vector::new(0.0, 0.0, 1e-9),
            frame.origin(),
            on + ogeom_math::Vector::new(50.0, 3.0, -8.0),
        ] {
            let found = foot_via_samples(&curve, &sampled, &segments, p, T);
            let seed = every(p);
            let expected = foot_on_3d(&curve, p, seed.1, T);
            assert_eq!(
                found.map(|(s, g)| (s.to_bits(), g.to_bits())),
                expected.map(|(s, g)| (s.to_bits(), g.to_bits())),
                "at {p:?}"
            );
        }
    }
}
