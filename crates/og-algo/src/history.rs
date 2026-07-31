//! What an operation did to its inputs.
//!
//! `docs/DATA_MODEL.md` §7. Every operation in this crate and above it reports
//! three things about each shape it was given:
//!
//! - **generated** — new entities made *from* it that did not exist before. A
//!   prism's side faces are generated from the profile's edges.
//! - **modified** — what it *became*. A face split in two is modified into both
//!   halves.
//! - **deleted** — it has no image in the result at all.
//!
//! This is not bookkeeping for its own sake. A parametric application records
//! "fillet *that* edge" and must still find that edge after the model is
//! rebuilt with different dimensions; it does so by walking history. Half-
//! populated history does not error — it reopens the document with the wrong
//! faces filleted — which is why every operation populates it from the commit
//! that introduces it rather than later.
//!
//! # Composition is the hard part
//!
//! Operations chain. If A modifies `x` into `y` and B then modifies `y` into
//! `z`, the composed history must say A-then-B modified `x` into `z`. Getting
//! that wrong is how a reference survives one rebuild and dies on the next,
//! which is far harder to diagnose than dying immediately. See
//! [`History::then`].

use std::collections::{HashMap, HashSet};

use og_topo::{SameKey, Shape};

/// A record of what one operation, or a chain of them, did.
///
/// Keyed by [`SameKey`] — node and placement, ignoring orientation. An edge and
/// its reverse are the same edge, and history about one is history about both;
/// keying on orientation would silently split every record in two.
#[derive(Debug, Clone, Default)]
pub struct History {
    generated: HashMap<SameKey, Vec<Shape>>,
    modified: HashMap<SameKey, Vec<Shape>>,
    deleted: HashSet<SameKey>,
}

impl History {
    /// An empty history.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// A history in which nothing happened to anything.
    ///
    /// The identity for [`History::then`]: composing with it changes nothing.
    #[must_use]
    pub fn identity() -> Self {
        Self::new()
    }

    /// Record that `input` produced `output` as a new entity.
    pub fn generate(&mut self, input: &Shape, output: Shape) {
        self.generated
            .entry(SameKey(input.clone()))
            .or_default()
            .push(output);
    }

    /// Record that `input` became `output`.
    ///
    /// Withdraws any deletion recorded for the same shape. Deletion and
    /// modification are contradictory claims, and the guard has to run both
    /// ways: clearing modifications on delete but not deletions on modify
    /// leaves a history that says both, and two callers reach opposite
    /// conclusions from it.
    pub fn modify(&mut self, input: &Shape, output: Shape) {
        let key = SameKey(input.clone());
        self.deleted.remove(&key);
        self.modified.entry(key).or_default().push(output);
    }

    /// Record that `input` has no image in the result.
    ///
    /// A shape cannot be both deleted and modified: if it became something, it
    /// was not deleted. Recording a deletion drops any modification record for
    /// the same shape, so the two can never disagree.
    ///
    /// Deletion and *generation* are a different matter, and coexist freely: a
    /// swept profile edge is consumed by the sweep — deleted — while generating
    /// the side face that grew from it. Anything that treats a deletion as the
    /// end of the story about a shape loses that face's ancestry.
    pub fn delete(&mut self, input: &Shape) {
        let key = SameKey(input.clone());
        self.modified.remove(&key);
        self.deleted.insert(key);
    }

    /// New entities made from `input`.
    #[must_use]
    pub fn generated(&self, input: &Shape) -> &[Shape] {
        self.generated
            .get(&SameKey(input.clone()))
            .map_or(&[], Vec::as_slice)
    }

    /// What `input` became.
    ///
    /// Empty for a shape the operation left alone: "unchanged" and "modified
    /// into nothing" are different, and the second is [`History::is_deleted`].
    #[must_use]
    pub fn modified(&self, input: &Shape) -> &[Shape] {
        self.modified
            .get(&SameKey(input.clone()))
            .map_or(&[], Vec::as_slice)
    }

    /// Whether `input` has no image in the result.
    #[must_use]
    pub fn is_deleted(&self, input: &Shape) -> bool {
        self.deleted.contains(&SameKey(input.clone()))
    }

    /// Whether the operation touched `input` at all.
    #[must_use]
    pub fn is_affected(&self, input: &Shape) -> bool {
        let key = SameKey(input.clone());
        self.deleted.contains(&key)
            || self.modified.contains_key(&key)
            || self.generated.contains_key(&key)
    }

    /// Where `input` ended up: what it became, or itself if it was untouched.
    ///
    /// The question a caller resolving a stored reference actually has. A shape
    /// an operation ignored is still there, and reporting nothing for it would
    /// make every caller special-case the common path.
    ///
    /// Returns an empty slice only for a shape that was deleted.
    #[must_use]
    pub fn trace<'a>(&'a self, input: &'a Shape) -> &'a [Shape] {
        if self.is_deleted(input) {
            return &[];
        }
        let images = self.modified(input);
        if images.is_empty() {
            core::slice::from_ref(input)
        } else {
            images
        }
    }

    /// Every shape this history has something to say about.
    #[must_use]
    pub fn inputs(&self) -> Vec<Shape> {
        let mut out: Vec<Shape> = Vec::new();
        let mut seen = HashSet::new();
        for key in self
            .generated
            .keys()
            .chain(self.modified.keys())
            .chain(self.deleted.iter())
        {
            if seen.insert(key.clone()) {
                out.push(key.0.clone());
            }
        }
        out
    }

    /// Whether this history records nothing.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.generated.is_empty() && self.modified.is_empty() && self.deleted.is_empty()
    }

    /// This history followed by `later`.
    ///
    /// Composition, and the operation everything chained depends on. For each
    /// shape either history knows about, the result answers where it ended up
    /// after both steps.
    ///
    /// The subtlety is telling *unchanged* from *modified into itself*. A shape
    /// neither step touched must come out with no record at all, not a
    /// modification saying it became itself — otherwise composing with an empty
    /// history would invent records, and composition would not have an
    /// identity. So a modification is recorded only when one of the two steps
    /// actually reported one.
    #[must_use]
    pub fn then(&self, later: &Self) -> Self {
        let mut out = Self::new();

        let mut subjects: Vec<Shape> = self.inputs();
        for input in later.inputs() {
            if !subjects.iter().any(|s| s.is_same(&input)) {
                subjects.push(input);
            }
        }

        for input in subjects {
            let gone_already = self.is_deleted(&input);

            // What this input is after the first step. Nothing, if it was
            // consumed; its images, if it changed; itself, if it was left alone.
            let first_images = self.modified(&input);
            let changed_by_self = !first_images.is_empty();
            let after_first: Vec<Shape> = if gone_already {
                Vec::new()
            } else if changed_by_self {
                first_images.to_vec()
            } else {
                vec![input.clone()]
            };

            if gone_already {
                out.delete(&input);
            } else {
                let mut changed_by_later = false;
                let mut final_images: Vec<Shape> = Vec::new();
                for image in &after_first {
                    if later.is_deleted(image) {
                        // Surviving the first step but not the second is a
                        // change, and one that ends in nothing.
                        changed_by_later = true;
                        continue;
                    }
                    if !later.modified(image).is_empty() {
                        changed_by_later = true;
                    }
                    final_images.extend(later.trace(image).iter().cloned());
                }

                if final_images.is_empty() {
                    out.delete(&input);
                } else if changed_by_self || changed_by_later {
                    for image in final_images {
                        out.modify(&input, image);
                    }
                }
                // Otherwise neither step touched it, and silence is the answer.
            }

            // Generation is reported whether or not the input survived: a
            // consumed edge still explains where the face that replaced it came
            // from, and that ancestry is the whole point of the record.
            for made in self.generated(&input) {
                for image in later.trace(made) {
                    out.generate(&input, image.clone());
                }
            }
            for image in &after_first {
                for made in later.generated(image) {
                    out.generate(&input, made.clone());
                }
            }
        }

        out
    }

    /// Fold a sequence of histories into one, in order.
    #[must_use]
    pub fn chain(steps: &[Self]) -> Self {
        steps
            .iter()
            .fold(Self::identity(), |acc, step| acc.then(step))
    }
}

/// A shape together with the history of the operation that produced it.
///
/// The return type of every operation. Bundling them means an operation cannot
/// return a result without saying what it did to get there.
#[derive(Debug, Clone)]
pub struct Built {
    /// The result.
    pub shape: Shape,
    /// What the operation did to its inputs.
    pub history: History,
}

impl Built {
    /// A result with its history.
    #[must_use]
    pub const fn new(shape: Shape, history: History) -> Self {
        Self { shape, history }
    }

    /// A result of an operation that had no inputs to report on.
    ///
    /// For a primitive built from numbers rather than from existing topology —
    /// there is nothing for the history to say.
    #[must_use]
    pub fn from_nothing(shape: Shape) -> Self {
        Self::new(shape, History::identity())
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use og_math::Point;
    use og_topo::Model;

    fn shapes(n: usize) -> (Model, Vec<Shape>) {
        let mut model = Model::new();
        #[allow(clippy::cast_precision_loss)]
        let shapes = (0..n)
            .map(|i| model.add_point(Point::new(i as f64, 0.0, 0.0)))
            .collect();
        (model, shapes)
    }

    #[test]
    fn an_empty_history_reports_everything_as_untouched() {
        let (_, s) = shapes(1);
        let h = History::new();
        assert!(h.is_empty());
        assert!(!h.is_affected(&s[0]));
        assert!(!h.is_deleted(&s[0]));
        assert!(h.generated(&s[0]).is_empty());
        assert!(h.modified(&s[0]).is_empty());
        // Tracing an untouched shape finds the shape itself, not nothing.
        assert_eq!(h.trace(&s[0]).len(), 1);
        assert!(h.trace(&s[0])[0].is_same(&s[0]));
    }

    #[test]
    fn generated_modified_and_deleted_are_distinct_claims() {
        let (_, s) = shapes(4);
        let mut h = History::new();
        h.generate(&s[0], s[1].clone());
        h.modify(&s[2], s[3].clone());
        h.delete(&s[1]);

        assert_eq!(h.generated(&s[0]).len(), 1);
        assert!(h.modified(&s[0]).is_empty(), "generating is not modifying");
        assert_eq!(h.modified(&s[2]).len(), 1);
        assert!(h.is_deleted(&s[1]));
        assert!(!h.is_deleted(&s[2]));
    }

    #[test]
    fn a_shape_cannot_be_both_deleted_and_modified() {
        // If it became something it was not deleted, and a record claiming both
        // would let two callers reach opposite conclusions.
        let (_, s) = shapes(2);
        let mut h = History::new();
        h.modify(&s[0], s[1].clone());
        assert_eq!(h.modified(&s[0]).len(), 1);

        h.delete(&s[0]);
        assert!(h.is_deleted(&s[0]));
        assert!(
            h.modified(&s[0]).is_empty(),
            "the modification was withdrawn"
        );
        assert!(h.trace(&s[0]).is_empty());
    }

    #[test]
    fn history_is_keyed_ignoring_orientation() {
        // An edge and its reverse are the same edge. Keying on orientation
        // would split every record in two, so a caller holding the reversed
        // handle would find nothing.
        let (_, s) = shapes(2);
        let mut h = History::new();
        h.modify(&s[0], s[1].clone());

        let reversed = s[0].reversed();
        assert_eq!(h.modified(&reversed).len(), 1, "same edge, other way round");
        assert!(h.is_affected(&reversed));
    }

    #[test]
    fn composition_follows_a_shape_through_two_operations() {
        // a -> b in the first, b -> c in the second. The composed history must
        // say a -> c; anything else and a stored reference survives one rebuild
        // and dies on the next.
        let (_, s) = shapes(3);
        let (a, b, c) = (&s[0], &s[1], &s[2]);

        let mut first = History::new();
        first.modify(a, b.clone());
        let mut second = History::new();
        second.modify(b, c.clone());

        let composed = first.then(&second);
        assert_eq!(composed.modified(a).len(), 1);
        assert!(composed.modified(a)[0].is_same(c));
        assert!(!composed.is_deleted(a));
    }

    #[test]
    fn composition_reports_a_shape_deleted_by_the_second_step() {
        let (_, s) = shapes(2);
        let (a, b) = (&s[0], &s[1]);

        let mut first = History::new();
        first.modify(a, b.clone());
        let mut second = History::new();
        second.delete(b);

        let composed = first.then(&second);
        assert!(
            composed.is_deleted(a),
            "a survived the first step but not the second, so it is gone"
        );
        assert!(composed.trace(a).is_empty());
    }

    #[test]
    fn composition_keeps_a_deletion_from_the_first_step() {
        let (_, s) = shapes(2);
        let mut first = History::new();
        first.delete(&s[0]);
        let second = History::new();
        assert!(first.then(&second).is_deleted(&s[0]));
    }

    #[test]
    fn composition_splits_when_the_second_step_splits() {
        // a -> b, then b -> {c, d}. The composed answer is a -> {c, d}, which is
        // what a caller filleting "that edge" needs after two rebuilds.
        let (_, s) = shapes(4);
        let (a, b, c, d) = (&s[0], &s[1], &s[2], &s[3]);

        let mut first = History::new();
        first.modify(a, b.clone());
        let mut second = History::new();
        second.modify(b, c.clone());
        second.modify(b, d.clone());

        let composed = first.then(&second);
        assert_eq!(composed.modified(a).len(), 2);
        assert!(composed.modified(a).iter().any(|s| s.is_same(c)));
        assert!(composed.modified(a).iter().any(|s| s.is_same(d)));
    }

    #[test]
    fn composition_carries_generated_entities_forward() {
        // The first step generates b from a; the second turns b into c. What a
        // generated from the pair is c, not b — the intermediate is gone.
        let (_, s) = shapes(3);
        let (a, b, c) = (&s[0], &s[1], &s[2]);

        let mut first = History::new();
        first.generate(a, b.clone());
        let mut second = History::new();
        second.modify(b, c.clone());

        let composed = first.then(&second);
        assert_eq!(composed.generated(a).len(), 1);
        assert!(composed.generated(a)[0].is_same(c));
    }

    #[test]
    fn composition_passes_through_shapes_only_the_later_step_knows() {
        let (_, s) = shapes(3);
        let (a, b, c) = (&s[0], &s[1], &s[2]);

        let mut first = History::new();
        first.modify(a, a.clone());
        let mut second = History::new();
        second.modify(b, c.clone());

        let composed = first.then(&second);
        assert_eq!(composed.modified(b).len(), 1, "the first step never saw b");
        assert!(composed.modified(b)[0].is_same(c));
    }

    #[test]
    fn the_empty_history_is_an_identity_for_composition() {
        let (_, s) = shapes(3);
        let mut h = History::new();
        h.modify(&s[0], s[1].clone());
        h.generate(&s[0], s[2].clone());
        h.delete(&s[1]);

        for composed in [h.then(&History::identity()), History::identity().then(&h)] {
            assert_eq!(composed.modified(&s[0]).len(), h.modified(&s[0]).len());
            assert_eq!(composed.generated(&s[0]).len(), h.generated(&s[0]).len());
            assert_eq!(composed.is_deleted(&s[1]), h.is_deleted(&s[1]));
        }
    }

    #[test]
    fn composition_is_associative() {
        // Three chained operations must give the same answer however the chain
        // is bracketed, or a caller that batches differently gets a different
        // model.
        let (_, s) = shapes(4);
        let (a, b, c, d) = (&s[0], &s[1], &s[2], &s[3]);

        let mut one = History::new();
        one.modify(a, b.clone());
        let mut two = History::new();
        two.modify(b, c.clone());
        let mut three = History::new();
        three.modify(c, d.clone());

        let left = one.then(&two).then(&three);
        let right = one.then(&two.then(&three));

        assert_eq!(left.modified(a).len(), right.modified(a).len());
        assert!(left.modified(a)[0].is_same(&right.modified(a)[0]));
        assert!(left.modified(a)[0].is_same(d));
    }

    #[test]
    fn chaining_a_sequence_matches_folding_it_by_hand() {
        let (_, s) = shapes(4);
        let mut steps = Vec::new();
        for i in 0..3 {
            let mut h = History::new();
            h.modify(&s[i], s[i + 1].clone());
            steps.push(h);
        }
        let chained = History::chain(&steps);
        assert!(chained.modified(&s[0])[0].is_same(&s[3]));
        assert!(History::chain(&[]).is_empty());
    }

    #[test]
    fn inputs_lists_every_shape_the_history_mentions_once() {
        let (_, s) = shapes(3);
        let mut h = History::new();
        h.modify(&s[0], s[1].clone());
        h.generate(&s[0], s[2].clone());
        h.delete(&s[1]);

        let inputs = h.inputs();
        assert_eq!(inputs.len(), 2, "s0 appears twice but is listed once");
        assert!(inputs.iter().any(|x| x.is_same(&s[0])));
        assert!(inputs.iter().any(|x| x.is_same(&s[1])));
    }

    #[test]
    fn a_consumed_shape_still_reports_what_it_generated() {
        // A sweep consumes its profile edge and grows a side face from it. Both
        // facts are true at once, and losing the second loses that face's
        // ancestry — which is what a rebuild needs to find it again.
        let (_, s) = shapes(2);
        let (edge, face) = (&s[0], &s[1]);

        let mut sweep = History::new();
        sweep.generate(edge, face.clone());
        sweep.delete(edge);

        assert!(sweep.is_deleted(edge));
        assert_eq!(
            sweep.generated(edge).len(),
            1,
            "deletion is not the end of the story"
        );

        // And composition must carry it through.
        let composed = sweep.then(&History::identity());
        assert!(composed.is_deleted(edge));
        assert_eq!(composed.generated(edge).len(), 1);
        assert!(composed.generated(edge)[0].is_same(face));
    }

    #[test]
    fn a_modification_withdraws_an_earlier_deletion() {
        // The guard has to run both ways. Clearing modifications on delete but
        // not deletions on modify leaves a history saying both, and two callers
        // reach opposite conclusions from it.
        let (_, s) = shapes(2);
        let mut h = History::new();
        h.delete(&s[0]);
        h.modify(&s[0], s[1].clone());
        assert!(!h.is_deleted(&s[0]));
        assert_eq!(h.modified(&s[0]).len(), 1);
    }

    #[test]
    fn a_primitive_reports_a_history_with_nothing_in_it() {
        // Built from numbers, not from topology: there is nothing to say, and
        // saying nothing is different from failing to say anything.
        let (_, s) = shapes(1);
        let built = Built::from_nothing(s[0].clone());
        assert!(built.history.is_empty());
        assert!(built.shape.is_same(&s[0]));
    }
}
