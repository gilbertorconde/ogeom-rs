//! Property tests for the `docs/DATA_MODEL.md` invariants.
//!
//! These check laws rather than examples: an arena key never aliases, tolerances
//! only ever widen, exact predicates are antisymmetric and agree with themselves
//! under relabelling. Laws are what algorithms downstream will actually lean on.

#![allow(clippy::unwrap_used, clippy::expect_used)]

use og_core::{Arena, Exact, OpId, Predicates, ProvenanceTable, Role, Sign, Tolerance, Tolerances};
use proptest::prelude::*;

/// Finite, non-pathological coordinates. Excluding subnormals and huge
/// magnitudes keeps the tests about the predicates rather than about f64.
fn coord() -> impl Strategy<Value = f64> {
    prop_oneof![
        // Small integers: the regime where exact answers are checkable by hand.
        (-20i32..20).prop_map(f64::from),
        -1e6f64..1e6f64,
    ]
}

fn p2() -> impl Strategy<Value = [f64; 2]> {
    (coord(), coord()).prop_map(|(x, y)| [x, y])
}

fn p3() -> impl Strategy<Value = [f64; 3]> {
    (coord(), coord(), coord()).prop_map(|(x, y, z)| [x, y, z])
}

proptest! {
    /// §11 — a key resolves to what was inserted under it, whatever else happened
    /// to the arena in between.
    #[test]
    fn arena_keys_resolve_to_their_own_value(values in prop::collection::vec(any::<u64>(), 1..64)) {
        let mut arena = Arena::new();
        let keys: Vec<_> = values.iter().map(|&v| arena.insert(v)).collect();
        for (key, &value) in keys.iter().zip(&values) {
            prop_assert_eq!(arena.get(*key), Some(&value));
        }
    }

    /// §11 — the central safety claim of generational indices: after a slot is
    /// freed and reused, the old key must not resolve. This is the property that
    /// turns a silent wrong answer into a clean `None`.
    #[test]
    fn removed_keys_never_alias_later_inserts(
        n in 1usize..32,
        remove_at in prop::collection::vec(any::<prop::sample::Index>(), 1..16),
    ) {
        let mut arena = Arena::new();
        let mut live: Vec<_> = (0..n).map(|i| (arena.insert(i), i)).collect();
        let mut dead = Vec::new();

        for idx in remove_at {
            if live.is_empty() {
                break;
            }
            let (key, _) = live.remove(idx.index(live.len()));
            arena.remove(key);
            dead.push(key);

            // Refill, which may well reuse the slot we just freed.
            let value = 1000 + dead.len();
            live.push((arena.insert(value), value));
        }

        for key in dead {
            prop_assert!(arena.get(key).is_none(), "stale key {key:?} resolved");
        }
        for (key, value) in live {
            prop_assert_eq!(arena.get(key), Some(&value));
        }
    }

    /// §5 — widening is a join: commutative, associative, idempotent, and never
    /// decreasing. Boolean operations lean on this when they inflate tolerances.
    #[test]
    fn tolerance_widening_is_a_join(a in 0.0f64..1e3, b in 0.0f64..1e3, c in 0.0f64..1e3) {
        let (a, b, c) = (
            Tolerance::new(a).unwrap(),
            Tolerance::new(b).unwrap(),
            Tolerance::new(c).unwrap(),
        );
        prop_assert_eq!(a.widen(b), b.widen(a), "commutative");
        prop_assert_eq!(a.widen(b).widen(c), a.widen(b.widen(c)), "associative");
        prop_assert_eq!(a.widen(a), a, "idempotent");
        prop_assert!(a.widen(b).get() >= a.get(), "never shrinks");
        prop_assert!(a.widen(b).get() >= b.get(), "never shrinks");
    }

    /// §5 — a tolerance always covers separations up to its own magnitude, and
    /// never covers anything beyond it.
    #[test]
    fn tolerance_covers_exactly_its_radius(t in 1e-6f64..1e3, d in 0.0f64..1e4) {
        let tol = Tolerance::new(t).unwrap();
        prop_assert_eq!(tol.covers(d), d <= tol.get());
        prop_assert_eq!(tol.covers(-d), tol.covers(d), "sign-independent");
    }

    /// §5 — scale conversion is exact enough to round-trip. A model in inches
    /// must not accumulate tolerance drift simply by being read back.
    #[test]
    fn tolerances_scale_consistently(scale in 1e-3f64..1e4) {
        let t = Tolerances::with_scale(scale).unwrap();
        prop_assert_eq!(t.scale(), scale);
        prop_assert!(t.confusion() > 0.0 && t.confusion().is_finite());
        // Angular and parametric tolerances are dimensionless and must not move.
        prop_assert_eq!(t.angular(), Tolerances::millimetres().angular());
        prop_assert_eq!(t.parametric(), Tolerances::millimetres().parametric());
    }

    /// §9 — swapping two arguments of an orientation predicate flips its sign.
    /// Exactly, with no epsilon: this is what makes a triangulation's
    /// combinatorial decisions consistent no matter which way an edge is walked.
    #[test]
    fn orient2d_is_antisymmetric(a in p2(), b in p2(), c in p2()) {
        prop_assert_eq!(Exact::orient2d(a, b, c), Exact::orient2d(b, a, c).reversed());
        prop_assert_eq!(Exact::orient2d(a, b, c), Exact::orient2d(a, c, b).reversed());
        // An even permutation preserves it.
        prop_assert_eq!(Exact::orient2d(a, b, c), Exact::orient2d(b, c, a));
    }

    /// §9 — likewise in 3D.
    #[test]
    fn orient3d_is_antisymmetric(a in p3(), b in p3(), c in p3(), d in p3()) {
        prop_assert_eq!(Exact::orient3d(a, b, c, d), Exact::orient3d(b, a, c, d).reversed());
        prop_assert_eq!(Exact::orient3d(a, b, c, d), Exact::orient3d(a, b, d, c).reversed());
    }

    /// §9 — a repeated point is a degenerate configuration, always, with no
    /// dependence on magnitude or ordering.
    #[test]
    fn degenerate_inputs_are_exactly_zero(a in p2(), b in p2(), x in p3(), y in p3()) {
        prop_assert_eq!(Exact::orient2d(a, a, b), Sign::Zero);
        prop_assert_eq!(Exact::orient2d(a, b, b), Sign::Zero);
        prop_assert_eq!(Exact::orient3d(x, x, y, y), Sign::Zero);
    }

    /// §9 — a point that provably lies on the line through `a` and `b` is
    /// reported as collinear *exactly*, not merely nearly, across ten orders of
    /// magnitude. Naive determinant evaluation loses this at scale.
    ///
    /// The construction has to be exact in binary floating point or the test
    /// measures rounding rather than the predicate: coordinates are small
    /// integers times a power of two, so every product and sum below is
    /// representable and `c` genuinely lies on the line.
    #[test]
    fn exact_collinearity_holds_at_scale(
        ax in -1000i32..1000, ay in -1000i32..1000,
        dx in -1000i32..1000, dy in -1000i32..1000,
        k in -8i32..8,
        exponent in -20i32..20,
    ) {
        prop_assume!(dx != 0 || dy != 0);
        let scale = f64::powi(2.0, exponent);
        let at = |v: i32| f64::from(v) * scale;

        let a = [at(ax), at(ay)];
        let b = [at(ax + dx), at(ay + dy)];
        let c = [at(ax + k * dx), at(ay + k * dy)];

        prop_assert!(
            Exact::are_collinear(a, b, c),
            "a={a:?} b={b:?} c={c:?} reported non-collinear"
        );
    }

    /// §8 — provenance is append-only and every recorded id resolves in the table
    /// that issued it, with ids strictly increasing.
    #[test]
    fn provenance_ids_are_unique_and_resolvable(n in 1usize..128) {
        let mut table = ProvenanceTable::new();
        let mut ids = Vec::new();
        for i in 0..n {
            let id = if ids.is_empty() {
                table.primitive(OpId(0), Role::SOLE)
            } else {
                #[allow(clippy::cast_possible_truncation)]
                table.derived(OpId(i as u32), [ids[i - 1]], Role::op_defined(0))
            };
            ids.push(id);
        }
        prop_assert_eq!(table.len(), n);
        for pair in ids.windows(2) {
            prop_assert!(pair[1] > pair[0], "ids must increase");
        }
        for id in &ids {
            prop_assert!(table.get(*id).is_some());
        }
        // Every entity in one chain traces back to the single primitive root.
        prop_assert_eq!(table.roots(*ids.last().unwrap()), vec![ids[0]]);
    }
}
