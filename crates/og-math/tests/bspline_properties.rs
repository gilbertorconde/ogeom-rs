//! Property tests for the B-spline machinery.
//!
//! The invariant that matters above all others: refinement never moves the
//! curve. Knot insertion, splitting, Bézier decomposition and degree elevation
//! all change the representation and must leave the geometry bit-for-bit
//! indistinguishable within tolerance. If any of them drifts, every algorithm
//! built on top inherits the drift.

#![allow(clippy::unwrap_used, clippy::expect_used)]

use og_core::Tolerances;
use og_math::{
    KnotVector, Point, Weighted,
    bspline::{
        derivatives, elevate_degree, evaluate, evaluate_rational, insert_knot,
        rational_derivatives, reverse, split, to_bezier_segments,
    },
};
use proptest::prelude::*;

const T: Tolerances = Tolerances::millimetres();

fn coord() -> impl Strategy<Value = f64> {
    -50.0f64..50.0
}

fn point() -> impl Strategy<Value = Point> {
    (coord(), coord(), coord()).prop_map(|(x, y, z)| Point::new(x, y, z))
}

/// A clamped uniform B-spline of arbitrary degree and control polygon.
fn curve() -> impl Strategy<Value = (KnotVector, Vec<Point>)> {
    (1usize..5, 2usize..10).prop_flat_map(|(degree, extra)| {
        let count = degree + extra;
        prop::collection::vec(point(), count..=count).prop_map(move |control| {
            (
                KnotVector::clamped_uniform(degree, control.len()).unwrap(),
                control,
            )
        })
    })
}

/// A rational curve: the same, with arbitrary positive weights.
fn rational_curve() -> impl Strategy<Value = (KnotVector, Vec<Weighted<Point>>)> {
    (1usize..4, 2usize..8).prop_flat_map(|(degree, extra)| {
        let count = degree + extra;
        (
            prop::collection::vec(point(), count..=count),
            prop::collection::vec(0.1f64..10.0, count..=count),
        )
            .prop_map(move |(points, weights)| {
                let control: Vec<_> = points
                    .into_iter()
                    .zip(weights)
                    .map(|(p, w)| Weighted::new(p, w, T).unwrap())
                    .collect();
                (
                    KnotVector::clamped_uniform(degree, control.len()).unwrap(),
                    control,
                )
            })
    })
}

/// Sample a curve evenly across its domain.
fn sample(knots: &KnotVector, control: &[Point], n: usize) -> Vec<Point> {
    let (a, b) = knots.domain();
    (0..=n)
        .map(|i| {
            #[allow(clippy::cast_precision_loss)]
            let u = a + (b - a) * (i as f64 / n as f64);
            evaluate(knots, control, u, T).unwrap()
        })
        .collect()
}

/// Curves here span coordinates up to 50 with degrees up to 5, so the tolerance
/// is stated relative to that scale rather than as an absolute distance.
fn same_curve(a: &[Point], b: &[Point]) -> bool {
    a.iter().zip(b).all(|(p, q)| {
        let scale = p.to_vector().magnitude().max(1.0);
        p.distance(*q) <= 1e-10 * scale
    })
}

proptest! {
    #[test]
    fn a_clamped_curve_starts_and_ends_at_its_end_control_points(
        (k, c) in curve(),
    ) {
        let (a, b) = k.domain();
        prop_assert!(evaluate(&k, &c, a, T).unwrap().is_equal(c[0], T));
        prop_assert!(evaluate(&k, &c, b, T).unwrap().is_equal(c[c.len() - 1], T));
    }

    /// A B-spline lies inside the convex hull of its control points. Checked
    /// through the weaker but easy-to-state consequence: it lies inside their
    /// axis-aligned bounding box.
    #[test]
    fn a_curve_stays_within_its_control_polygons_bounds((k, c) in curve()) {
        let low = c.iter().fold(c[0], |acc, p| acc.min(*p));
        let high = c.iter().fold(c[0], |acc, p| acc.max(*p));
        for p in sample(&k, &c, 40) {
            for axis in 0..3 {
                let (v, lo, hi) = (
                    p.coord(axis).unwrap(),
                    low.coord(axis).unwrap(),
                    high.coord(axis).unwrap(),
                );
                prop_assert!(v >= lo - 1e-9 && v <= hi + 1e-9, "left the hull on axis {axis}");
            }
        }
    }

    /// The single most important invariant in the module.
    #[test]
    fn knot_insertion_never_moves_the_curve(
        (k, c) in curve(), position in 0.05f64..0.95, count in 1usize..4,
    ) {
        let (a, b) = k.domain();
        let u = a + (b - a) * position;
        let existing = k.multiplicity_of(u);
        prop_assume!(existing + count <= k.degree());

        let before = sample(&k, &c, 60);
        let (k2, c2) = insert_knot(&k, &c, u, count, T).unwrap();
        prop_assert_eq!(c2.len(), c.len() + count);
        prop_assert_eq!(k2.multiplicity_of(u), existing + count);
        prop_assert!(same_curve(&before, &sample(&k2, &c2, 60)));
    }

    #[test]
    fn inserting_one_at_a_time_matches_inserting_all_at_once(
        (k, c) in curve(), position in 0.05f64..0.95,
    ) {
        prop_assume!(k.degree() >= 2);
        let (a, b) = k.domain();
        let u = a + (b - a) * position;
        prop_assume!(k.multiplicity_of(u) == 0);

        let (ka, ca) = insert_knot(&k, &c, u, 2, T).unwrap();
        let (k1, c1) = insert_knot(&k, &c, u, 1, T).unwrap();
        let (kb, cb) = insert_knot(&k1, &c1, u, 1, T).unwrap();

        prop_assert_eq!(ka.knots(), kb.knots());
        prop_assert!(same_curve(&sample(&ka, &ca, 40), &sample(&kb, &cb, 40)));
    }

    #[test]
    fn insertion_beyond_the_degree_is_refused(
        (k, c) in curve(), position in 0.05f64..0.95,
    ) {
        let (a, b) = k.domain();
        let u = a + (b - a) * position;
        prop_assert!(insert_knot(&k, &c, u, k.degree() + 1, T).is_err());
    }

    #[test]
    fn splitting_reproduces_both_halves((k, c) in curve(), position in 0.1f64..0.9) {
        let (a, b) = k.domain();
        let u = a + (b - a) * position;
        let ((lk, lc), (rk, rc)) = split(&k, &c, u, T).unwrap();

        prop_assert!(lk.is_clamped() && rk.is_clamped());
        prop_assert!((lk.domain().1 - u).abs() <= 1e-12);
        prop_assert!((rk.domain().0 - u).abs() <= 1e-12);

        for i in 0..=30 {
            let t = f64::from(i) / 30.0;
            for (kk, cc, lo, hi) in [(&lk, &lc, a, u), (&rk, &rc, u, b)] {
                let v = lo + (hi - lo) * t;
                let half = evaluate(kk, cc, v, T).unwrap();
                let whole = evaluate(&k, &c, v, T).unwrap();
                let scale = whole.to_vector().magnitude().max(1.0);
                prop_assert!(half.distance(whole) <= 1e-9 * scale, "diverges at {}", v);
            }
        }
    }

    #[test]
    fn the_two_halves_meet_at_the_split_point((k, c) in curve(), position in 0.1f64..0.9) {
        let (a, b) = k.domain();
        let u = a + (b - a) * position;
        let ((lk, lc), (rk, rc)) = split(&k, &c, u, T).unwrap();
        // The left curve's last control point is the right curve's first, so
        // the halves are joined rather than merely adjacent.
        prop_assert!(lc[lc.len() - 1].is_equal(rc[0], T));
        prop_assert!(
            evaluate(&lk, &lc, lk.domain().1, T)
                .unwrap()
                .is_equal(evaluate(&rk, &rc, rk.domain().0, T).unwrap(), T)
        );
    }

    #[test]
    fn bezier_decomposition_reproduces_the_curve((k, c) in curve()) {
        let segments = to_bezier_segments(&k, &c, T).unwrap();
        prop_assert!(!segments.is_empty());
        for (_, points) in &segments {
            prop_assert_eq!(points.len(), k.degree() + 1);
        }
        // Segments tile the domain end to end.
        prop_assert!((segments[0].0.0 - k.domain().0).abs() <= 1e-12);
        prop_assert!((segments[segments.len() - 1].0.1 - k.domain().1).abs() <= 1e-12);
        for pair in segments.windows(2) {
            prop_assert!((pair[0].0.1 - pair[1].0.0).abs() <= 1e-12);
        }

        for ((a, b), points) in &segments {
            let bezier = KnotVector::clamped_uniform(k.degree(), points.len())
                .unwrap()
                .reparameterized(*a, *b)
                .unwrap();
            for i in 0..=15 {
                let u = a + (b - a) * (f64::from(i) / 15.0);
                let piece = evaluate(&bezier, points, u, T).unwrap();
                let whole = evaluate(&k, &c, u, T).unwrap();
                let scale = whole.to_vector().magnitude().max(1.0);
                prop_assert!(piece.distance(whole) <= 1e-9 * scale);
            }
        }
    }

    #[test]
    fn degree_elevation_never_moves_the_curve((k, c) in curve()) {
        let before = sample(&k, &c, 60);
        let (k2, c2) = elevate_degree(&k, &c, T).unwrap();
        prop_assert_eq!(k2.degree(), k.degree() + 1);
        prop_assert!((k2.domain().0 - k.domain().0).abs() <= 1e-12);
        prop_assert!((k2.domain().1 - k.domain().1).abs() <= 1e-12);
        prop_assert!(same_curve(&before, &sample(&k2, &c2, 60)));
    }

    #[test]
    fn reversal_traverses_the_same_points_backwards((k, c) in curve()) {
        let (rk, rc) = reverse(&k, &c);
        let (a, b) = k.domain();
        for i in 0..=30 {
            let t = f64::from(i) / 30.0;
            let forward = evaluate(&k, &c, a + (b - a) * t, T).unwrap();
            let backward = evaluate(&rk, &rc, a + (b - a) * (1.0 - t), T).unwrap();
            let scale = forward.to_vector().magnitude().max(1.0);
            prop_assert!(forward.distance(backward) <= 1e-10 * scale);
        }
    }

    #[test]
    fn reversing_twice_restores_the_curve((k, c) in curve()) {
        let (rk, rc) = reverse(&k, &c);
        let (kk, cc) = reverse(&rk, &rc);
        prop_assert!(same_curve(&sample(&k, &c, 40), &sample(&kk, &cc, 40)));
    }

    #[test]
    fn derivatives_agree_with_central_differences((k, c) in curve(), position in 0.15f64..0.85) {
        let (a, b) = k.domain();
        let u = a + (b - a) * position;
        let h = 1e-6;
        let d = derivatives(&k, &c, u, 1, T).unwrap();

        prop_assert!(d[0].is_equal(evaluate(&k, &c, u, T).unwrap(), T));
        let numeric = (evaluate(&k, &c, u + h, T).unwrap()
            - evaluate(&k, &c, u - h, T).unwrap())
            * (1.0 / (2.0 * h));
        let scale = numeric.magnitude().max(1.0);
        prop_assert!((d[1].to_vector() - numeric).magnitude() <= 1e-4 * scale);
    }

    #[test]
    fn derivatives_above_the_degree_vanish((k, c) in curve(), position in 0.15f64..0.85) {
        let (a, b) = k.domain();
        let u = a + (b - a) * position;
        let n = k.degree() + 2;
        let d = derivatives(&k, &c, u, n, T).unwrap();
        for above in &d[k.degree() + 1..=n] {
            prop_assert!(above.to_vector().magnitude() <= 1e-9);
        }
    }

    #[test]
    fn rational_evaluation_is_unaffected_by_a_uniform_weight_scaling(
        (k, c) in rational_curve(), factor in 0.1f64..10.0,
    ) {
        // Multiplying every weight by the same amount is the same projective
        // point, so the curve must not move.
        let scaled: Vec<_> = c
            .iter()
            .map(|w| Weighted::new(w.point(), w.weight * factor, T).unwrap())
            .collect();
        let (a, b) = k.domain();
        for i in 0..=30 {
            let u = a + (b - a) * (f64::from(i) / 30.0);
            let p = evaluate_rational(&k, &c, u, T).unwrap();
            let q = evaluate_rational(&k, &scaled, u, T).unwrap();
            prop_assert!(p.distance(q) <= 1e-9 * p.to_vector().magnitude().max(1.0));
        }
    }

    #[test]
    fn knot_insertion_never_moves_a_rational_curve(
        (k, c) in rational_curve(), position in 0.05f64..0.95,
    ) {
        let (a, b) = k.domain();
        let u = a + (b - a) * position;
        prop_assume!(k.multiplicity_of(u) == 0);

        let (k2, c2) = insert_knot(&k, &c, u, 1, T).unwrap();
        for i in 0..=30 {
            let v = a + (b - a) * (f64::from(i) / 30.0);
            let p = evaluate_rational(&k, &c, v, T).unwrap();
            let q = evaluate_rational(&k2, &c2, v, T).unwrap();
            prop_assert!(p.distance(q) <= 1e-9 * p.to_vector().magnitude().max(1.0));
        }
    }

    #[test]
    fn rational_derivatives_agree_with_central_differences(
        (k, c) in rational_curve(), position in 0.15f64..0.85,
    ) {
        let (a, b) = k.domain();
        let u = a + (b - a) * position;
        let h = 1e-6;
        let d = rational_derivatives(&k, &c, u, 1, T).unwrap();

        prop_assert!(d[0].is_equal(evaluate_rational(&k, &c, u, T).unwrap(), T));
        let numeric = (evaluate_rational(&k, &c, u + h, T).unwrap()
            - evaluate_rational(&k, &c, u - h, T).unwrap())
            * (1.0 / (2.0 * h));
        let scale = numeric.magnitude().max(1.0);
        prop_assert!(
            (d[1].to_vector() - numeric).magnitude() <= 1e-3 * scale,
            "{:?} vs {numeric:?}",
            d[1]
        );
    }

    /// A rational curve with all weights equal is a non-rational one, so the
    /// two evaluation paths must agree exactly in shape.
    #[test]
    fn unit_weights_reduce_the_rational_case_to_the_polynomial_one((k, c) in curve()) {
        let weighted: Vec<_> = c.iter().map(|p| Weighted::new(*p, 1.0, T).unwrap()).collect();
        let (a, b) = k.domain();
        for i in 0..=30 {
            let u = a + (b - a) * (f64::from(i) / 30.0);
            let plain = evaluate(&k, &c, u, T).unwrap();
            let rational = evaluate_rational(&k, &weighted, u, T).unwrap();
            prop_assert!(plain.distance(rational) <= 1e-11 * plain.to_vector().magnitude().max(1.0));
        }
    }
}
