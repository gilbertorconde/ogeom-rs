//! Property tests for history composition.
//!
//! Composition is the one operation here with enough structure to get subtly
//! wrong, and the failure mode is nasty: a reference that survives one rebuild
//! and dies on the next is far harder to diagnose than one that dies at once.
//! So the laws are checked over randomized chains rather than the handful of
//! shapes a unit test can enumerate.

#![allow(clippy::unwrap_used, clippy::expect_used)]

use og_algo::History;
use og_math::Point;
use og_topo::{Model, Shape};
use proptest::prelude::*;

/// A pool of distinct shapes to build histories over.
fn pool(n: usize) -> Vec<Shape> {
    let mut model = Model::new();
    #[allow(clippy::cast_precision_loss)]
    (0..n)
        .map(|i| model.add_point(Point::new(i as f64, 0.0, 0.0)))
        .collect()
}

/// What one operation did to one shape.
#[derive(Debug, Clone)]
enum Step {
    Modify(usize, Vec<usize>),
    Generate(usize, usize),
    Delete(usize),
}

fn step(size: usize) -> impl Strategy<Value = Step> {
    prop_oneof![
        (0..size, prop::collection::vec(0..size, 1..3)).prop_map(|(i, out)| Step::Modify(i, out)),
        (0..size, 0..size).prop_map(|(i, o)| Step::Generate(i, o)),
        (0..size).prop_map(Step::Delete),
    ]
}

fn history(size: usize) -> impl Strategy<Value = Vec<Step>> {
    prop::collection::vec(step(size), 0..6)
}

fn build(steps: &[Step], shapes: &[Shape]) -> History {
    let mut h = History::new();
    for s in steps {
        match s {
            Step::Modify(i, outs) => {
                for o in outs {
                    h.modify(&shapes[*i], shapes[*o].clone());
                }
            }
            Step::Generate(i, o) => h.generate(&shapes[*i], shapes[*o].clone()),
            Step::Delete(i) => h.delete(&shapes[*i]),
        }
    }
    h
}

/// Compare two histories by what they claim about every shape in the pool.
fn agrees(a: &History, b: &History, shapes: &[Shape]) -> Result<(), String> {
    for s in shapes {
        if a.is_deleted(s) != b.is_deleted(s) {
            return Err(format!("deletion disagrees for {s:?}"));
        }
        let (am, bm) = (a.modified(s), b.modified(s));
        if am.len() != bm.len() || !am.iter().all(|x| bm.iter().any(|y| x.is_same(y))) {
            return Err(format!("modification disagrees for {s:?}"));
        }
        let (ag, bg) = (a.generated(s), b.generated(s));
        if ag.len() != bg.len() || !ag.iter().all(|x| bg.iter().any(|y| x.is_same(y))) {
            return Err(format!("generation disagrees for {s:?}"));
        }
    }
    Ok(())
}

proptest! {
    /// The empty history changes nothing, on either side.
    #[test]
    fn identity_composes_away(steps in history(5)) {
        let shapes = pool(5);
        let h = build(&steps, &shapes);
        prop_assert!(agrees(&h.then(&History::identity()), &h, &shapes).is_ok());
        prop_assert!(agrees(&History::identity().then(&h), &h, &shapes).is_ok());
    }

    /// Bracketing a chain differently must not change what it says. A caller
    /// that batches three operations as (a,b),c and one that batches them as
    /// a,(b,c) have to end up with the same model.
    #[test]
    fn composition_is_associative(a in history(5), b in history(5), c in history(5)) {
        let shapes = pool(5);
        let (ha, hb, hc) = (build(&a, &shapes), build(&b, &shapes), build(&c, &shapes));
        let left = ha.then(&hb).then(&hc);
        let right = ha.then(&hb.then(&hc));
        if let Err(why) = agrees(&left, &right, &shapes) {
            return Err(TestCaseError::fail(why));
        }
    }

    /// Deletion is absorbing: once a shape is gone from the first step, no
    /// later step can bring it back.
    #[test]
    fn a_deletion_survives_anything_that_follows(later in history(5), victim in 0usize..5) {
        let shapes = pool(5);
        let mut first = History::new();
        first.delete(&shapes[victim]);
        let composed = first.then(&build(&later, &shapes));
        prop_assert!(composed.is_deleted(&shapes[victim]));
        prop_assert!(composed.trace(&shapes[victim]).is_empty());
    }

    /// A shape is never both deleted and modified, however the history was
    /// assembled or composed.
    #[test]
    fn deletion_and_modification_never_coexist(a in history(5), b in history(5)) {
        let shapes = pool(5);
        let composed = build(&a, &shapes).then(&build(&b, &shapes));
        for s in &shapes {
            prop_assert!(
                !(composed.is_deleted(s) && !composed.modified(s).is_empty()),
                "{s:?} is both deleted and modified"
            );
        }
    }

    /// Tracing is total: every shape either has somewhere it went, or is gone.
    /// A caller resolving a stored reference should never get an answer that
    /// means neither.
    #[test]
    fn tracing_answers_for_every_shape(a in history(5), b in history(5)) {
        let shapes = pool(5);
        let composed = build(&a, &shapes).then(&build(&b, &shapes));
        for s in &shapes {
            let trace = composed.trace(s);
            prop_assert_eq!(
                trace.is_empty(),
                composed.is_deleted(s),
                "an empty trace must mean deleted, and nothing else"
            );
        }
    }

    /// An untouched shape traces to itself. Without this every caller would
    /// have to special-case the common path.
    #[test]
    fn an_untouched_shape_traces_to_itself(steps in history(5)) {
        let shapes = pool(6);
        // Shape 5 is outside the range the steps can name, so nothing touches it.
        let h = build(&steps, &shapes);
        let untouched = &shapes[5];
        prop_assert!(!h.is_affected(untouched));
        prop_assert_eq!(h.trace(untouched).len(), 1);
        prop_assert!(h.trace(untouched)[0].is_same(untouched));
    }

    /// Every shape a history mentions appears in `inputs` exactly once.
    #[test]
    fn inputs_are_complete_and_free_of_duplicates(steps in history(5)) {
        let shapes = pool(5);
        let h = build(&steps, &shapes);
        let inputs = h.inputs();

        for s in &shapes {
            prop_assert_eq!(
                h.is_affected(s),
                inputs.iter().any(|i| i.is_same(s)),
                "affected and listed must agree for {:?}", s
            );
        }
        for (i, a) in inputs.iter().enumerate() {
            for b in &inputs[i + 1..] {
                prop_assert!(!a.is_same(b), "duplicate entry in inputs");
            }
        }
    }

    /// Orientation is not part of a history's key, so a caller holding the
    /// reversed handle to an edge finds the same record.
    #[test]
    fn records_are_reachable_through_either_orientation(steps in history(5)) {
        let shapes = pool(5);
        let h = build(&steps, &shapes);
        for s in &shapes {
            let flipped = s.reversed();
            prop_assert_eq!(h.is_deleted(s), h.is_deleted(&flipped));
            prop_assert_eq!(h.modified(s).len(), h.modified(&flipped).len());
            prop_assert_eq!(h.generated(s).len(), h.generated(&flipped).len());
        }
    }

    /// Folding a sequence matches composing it pairwise, left to right.
    #[test]
    fn chaining_matches_repeated_composition(steps in prop::collection::vec(history(4), 0..4)) {
        let shapes = pool(4);
        let histories: Vec<History> = steps.iter().map(|s| build(s, &shapes)).collect();
        let chained = History::chain(&histories);
        let folded = histories
            .iter()
            .fold(History::identity(), |acc, h| acc.then(h));
        prop_assert!(agrees(&chained, &folded, &shapes).is_ok());
    }
}
