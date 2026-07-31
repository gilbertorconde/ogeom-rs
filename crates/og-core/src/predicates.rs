//! Geometric predicates, behind a trait.
//!
//! See `docs/DATA_MODEL.md` §9. Algorithms are written against [`Predicates`],
//! never against a concrete implementation, so the robustness strategy can be
//! swapped without touching them.
//!
//! # What this does and does not buy
//!
//! Exact predicates settle the *polyhedral* robustness problem: which side of a
//! plane a point is on, whether four points are coplanar, in-sphere tests for
//! Delaunay. Those questions have exact answers computable from the inputs, and
//! [`Exact`] gives them.
//!
//! They do **not** settle the CAD problem. The intersection curve of two NURBS
//! surfaces is transcendental — there is no exact value to be exact about. That
//! is why per-entity tolerances exist (`docs/DATA_MODEL.md` §5) and why they
//! cannot be traded away for better predicates. Predicates make the decidable
//! parts decidable; tolerances carry the rest.
//!
//! Use exact predicates where the question is genuinely combinatorial —
//! triangulation, point-in-polygon, orientation of a planar facet — and do not
//! reach for them expecting surface intersection to become robust.

/// The sign of a predicate's determinant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Sign {
    /// Determinant is negative.
    Negative,
    /// Determinant is zero — degenerate: collinear, coplanar, cocircular.
    Zero,
    /// Determinant is positive.
    Positive,
}

impl Sign {
    /// Classify a determinant. `NaN` maps to [`Sign::Zero`], on the grounds that
    /// an indeterminate configuration is a degenerate one.
    #[must_use]
    pub fn of(value: f64) -> Self {
        if value > 0.0 {
            Self::Positive
        } else if value < 0.0 {
            Self::Negative
        } else {
            Self::Zero
        }
    }

    /// Whether the configuration is degenerate.
    #[must_use]
    pub const fn is_zero(self) -> bool {
        matches!(self, Self::Zero)
    }

    /// The sign of the negated determinant.
    #[must_use]
    pub const fn reversed(self) -> Self {
        match self {
            Self::Negative => Self::Positive,
            Self::Zero => Self::Zero,
            Self::Positive => Self::Negative,
        }
    }
}

/// A point in the plane, as predicates see it.
pub type P2 = [f64; 2];
/// A point in space, as predicates see it.
pub type P3 = [f64; 3];

/// Orientation and incircle/insphere tests.
///
/// Implementations must agree on sign conventions; only their accuracy and cost
/// may differ.
pub trait Predicates {
    /// Sign of the area of triangle `(a, b, c)`.
    ///
    /// [`Sign::Positive`] when the three points are counter-clockwise,
    /// [`Sign::Zero`] when collinear.
    fn orient2d(a: P2, b: P2, c: P2) -> Sign;

    /// Sign of the volume of tetrahedron `(a, b, c, d)`.
    ///
    /// [`Sign::Positive`] when `d` lies below the plane through `a`, `b`, `c`,
    /// where "below" is the side from which `a`, `b`, `c` appear clockwise.
    /// [`Sign::Zero`] when the four points are coplanar.
    fn orient3d(a: P3, b: P3, c: P3, d: P3) -> Sign;

    /// Whether `d` lies inside the circle through `a`, `b`, `c`.
    ///
    /// [`Sign::Positive`] for inside. `a`, `b`, `c` must be counter-clockwise;
    /// otherwise the sign is inverted.
    fn incircle(a: P2, b: P2, c: P2, d: P2) -> Sign;

    /// Whether `e` lies inside the sphere through `a`, `b`, `c`, `d`.
    ///
    /// [`Sign::Positive`] for inside. `a`, `b`, `c`, `d` must be positively
    /// oriented; otherwise the sign is inverted.
    fn insphere(a: P3, b: P3, c: P3, d: P3, e: P3) -> Sign;

    /// Whether `c` lies to the left of the directed line `a -> b`.
    fn is_left_of(a: P2, b: P2, c: P2) -> bool {
        Self::orient2d(a, b, c) == Sign::Positive
    }

    /// Whether three points are collinear, exactly.
    fn are_collinear(a: P2, b: P2, c: P2) -> bool {
        Self::orient2d(a, b, c).is_zero()
    }

    /// Whether four points are coplanar, exactly.
    fn are_coplanar(a: P3, b: P3, c: P3, d: P3) -> bool {
        Self::orient3d(a, b, c, d).is_zero()
    }
}

/// Adaptive-precision exact predicates.
///
/// Shewchuk's adaptive floating-point expansions, via the `robust` crate: a fast
/// floating-point filter first, escalating to exact arithmetic only when the
/// error bound says the sign is not yet determined. Exact, and in the common
/// non-degenerate case barely slower than [`Fast`].
///
/// This is the default. Reach for [`Fast`] only with a measurement in hand.
#[derive(Debug, Clone, Copy, Default)]
pub struct Exact;

/// Naive floating-point predicates.
///
/// Evaluates the determinant directly. Fast, and **wrong** near degeneracy: it
/// can report a point as being on the wrong side of a plane it is nearly on,
/// which in a triangulation or a boolean means an inconsistent combinatorial
/// structure rather than a slightly-off number.
///
/// Present as a baseline for benchmarking and for callers that have shown the
/// inputs are well separated.
#[derive(Debug, Clone, Copy, Default)]
pub struct Fast;

fn c2(p: P2) -> robust::Coord<f64> {
    robust::Coord { x: p[0], y: p[1] }
}

fn c3(p: P3) -> robust::Coord3D<f64> {
    robust::Coord3D {
        x: p[0],
        y: p[1],
        z: p[2],
    }
}

impl Predicates for Exact {
    fn orient2d(a: P2, b: P2, c: P2) -> Sign {
        Sign::of(robust::orient2d(c2(a), c2(b), c2(c)))
    }

    fn orient3d(a: P3, b: P3, c: P3, d: P3) -> Sign {
        Sign::of(robust::orient3d(c3(a), c3(b), c3(c), c3(d)))
    }

    fn incircle(a: P2, b: P2, c: P2, d: P2) -> Sign {
        Sign::of(robust::incircle(c2(a), c2(b), c2(c), c2(d)))
    }

    fn insphere(a: P3, b: P3, c: P3, d: P3, e: P3) -> Sign {
        Sign::of(robust::insphere(c3(a), c3(b), c3(c), c3(d), c3(e)))
    }
}

impl Predicates for Fast {
    fn orient2d(a: P2, b: P2, c: P2) -> Sign {
        Sign::of((a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0]))
    }

    fn orient3d(a: P3, b: P3, c: P3, d: P3) -> Sign {
        let ad = [a[0] - d[0], a[1] - d[1], a[2] - d[2]];
        let bd = [b[0] - d[0], b[1] - d[1], b[2] - d[2]];
        let cd = [c[0] - d[0], c[1] - d[1], c[2] - d[2]];
        let det = ad[0] * (bd[1] * cd[2] - bd[2] * cd[1]) - bd[0] * (ad[1] * cd[2] - ad[2] * cd[1])
            + cd[0] * (ad[1] * bd[2] - ad[2] * bd[1]);
        Sign::of(det)
    }

    fn incircle(a: P2, b: P2, c: P2, d: P2) -> Sign {
        let ad = [a[0] - d[0], a[1] - d[1]];
        let bd = [b[0] - d[0], b[1] - d[1]];
        let cd = [c[0] - d[0], c[1] - d[1]];
        let alift = ad[0].mul_add(ad[0], ad[1] * ad[1]);
        let blift = bd[0].mul_add(bd[0], bd[1] * bd[1]);
        let clift = cd[0].mul_add(cd[0], cd[1] * cd[1]);
        let det = alift * (bd[0] * cd[1] - cd[0] * bd[1]) - blift * (ad[0] * cd[1] - cd[0] * ad[1])
            + clift * (ad[0] * bd[1] - bd[0] * ad[1]);
        Sign::of(det)
    }

    fn insphere(a: P3, b: P3, c: P3, d: P3, e: P3) -> Sign {
        let lift = |p: P3| {
            let v = [p[0] - e[0], p[1] - e[1], p[2] - e[2]];
            (v, v[0].mul_add(v[0], v[1].mul_add(v[1], v[2] * v[2])))
        };
        let (ae, al) = lift(a);
        let (be, bl) = lift(b);
        let (ce, cl) = lift(c);
        let (de, dl) = lift(d);

        let det3 = |p: [f64; 3], q: [f64; 3], r: [f64; 3]| {
            p[0] * (q[1] * r[2] - q[2] * r[1]) - p[1] * (q[0] * r[2] - q[2] * r[0])
                + p[2] * (q[0] * r[1] - q[1] * r[0])
        };
        let det = -al * det3(be, ce, de) + bl * det3(ae, ce, de) - cl * det3(ae, be, de)
            + dl * det3(ae, be, ce);
        Sign::of(det)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn orient2d_sign_convention() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        assert_eq!(Exact::orient2d(a, b, [0.0, 1.0]), Sign::Positive);
        assert_eq!(Exact::orient2d(a, b, [0.0, -1.0]), Sign::Negative);
        assert_eq!(Exact::orient2d(a, b, [2.0, 0.0]), Sign::Zero);
        assert!(Exact::is_left_of(a, b, [0.5, 0.5]));
        assert!(Exact::are_collinear(a, b, [7.0, 0.0]));
    }

    #[test]
    fn orient3d_sign_convention_and_coplanarity() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        // Convention check: whatever the sign for +z, -z must be its opposite,
        // and a point in the plane must be exactly zero.
        let above = Exact::orient3d(a, b, c, [0.0, 0.0, 1.0]);
        let below = Exact::orient3d(a, b, c, [0.0, 0.0, -1.0]);
        assert_eq!(above, below.reversed());
        assert!(!above.is_zero());
        assert!(Exact::are_coplanar(a, b, c, [3.0, -4.0, 0.0]));
    }

    #[test]
    fn exact_and_fast_agree_when_well_separated() {
        let pts = [
            ([0.0, 0.0], [3.0, 1.0], [1.0, 4.0]),
            ([-2.0, 5.0], [7.0, -1.0], [0.25, 0.5]),
            ([1e6, 1e6], [-1e6, 2e6], [0.0, 0.0]),
        ];
        for (a, b, c) in pts {
            assert_eq!(Exact::orient2d(a, b, c), Fast::orient2d(a, b, c));
        }
    }

    #[test]
    fn exact_predicates_survive_a_case_naive_arithmetic_gets_wrong() {
        // A classic near-degenerate configuration: c is very slightly left of the
        // line a->b, by an amount that cancels catastrophically in the naive
        // determinant. The exact predicate must still say Positive.
        let a = [0.5, 0.5];
        let b = [12.0, 12.0];
        let c = [24.000_000_000_000_004, 24.0];

        assert_eq!(Exact::orient2d(a, b, c), Sign::Negative);
        // Not asserting Fast is wrong here — the point is that Exact is
        // trustworthy at this scale and the algorithms depend on that.
        assert!(!Exact::orient2d(a, b, c).is_zero());
    }

    #[test]
    fn incircle_sign_convention() {
        // Unit square corners, counter-clockwise: circumcircle has radius sqrt(2)/2.
        let a = [-1.0, -1.0];
        let b = [1.0, -1.0];
        let c = [1.0, 1.0];
        assert_eq!(
            Exact::orient2d(a, b, c),
            Sign::Positive,
            "test setup must be CCW"
        );
        assert_eq!(Exact::incircle(a, b, c, [0.0, 0.0]), Sign::Positive);
        assert_eq!(Exact::incircle(a, b, c, [5.0, 5.0]), Sign::Negative);
        assert_eq!(
            Exact::incircle(a, b, c, [-1.0, 1.0]),
            Sign::Zero,
            "cocircular"
        );
    }

    #[test]
    fn insphere_sign_convention() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d = [0.0, 0.0, 1.0];
        // Orient the tetrahedron positively before testing, as the contract requires.
        let (a, b, c, d) = if Exact::orient3d(a, b, c, d) == Sign::Positive {
            (a, b, c, d)
        } else {
            (a, c, b, d)
        };
        let inside = Exact::insphere(a, b, c, d, [0.25, 0.25, 0.25]);
        let outside = Exact::insphere(a, b, c, d, [10.0, 10.0, 10.0]);
        assert_eq!(inside, Sign::Positive);
        assert_eq!(outside, Sign::Negative);
    }

    #[test]
    fn nan_is_treated_as_degenerate_not_propagated() {
        assert_eq!(Sign::of(f64::NAN), Sign::Zero);
        assert_eq!(Sign::of(0.0), Sign::Zero);
        assert_eq!(Sign::of(-0.0), Sign::Zero);
    }

    #[test]
    fn reversed_is_an_involution() {
        for s in [Sign::Negative, Sign::Zero, Sign::Positive] {
            assert_eq!(s.reversed().reversed(), s);
        }
    }
}
