//! Placement as a chain of transforms.
//!
//! `docs/DATA_MODEL.md` §2. A [`Location`] is a *sequence* of `(datum, power)`
//! pairs, not a 4×4 matrix, and that is load-bearing rather than stylistic:
//!
//! - **Composition is concatenation.** No matrix product, and no drift from
//!   composing the same placement a thousand times down an assembly tree.
//! - **Identity is structural.** Two shapes are at the same place when their
//!   chains match — decided by comparing a handful of integers rather than
//!   sixteen floats against a tolerance. That is what lets ten thousand
//!   identical fasteners share one piece of geometry *and* be recognisable as
//!   instances of it.
//! - **Inverses are exact.** Negate the powers; no matrix inversion, no
//!   rounding.
//!
//! The composed [`Transform`] is derived on demand. It is deliberately *not*
//! cached inside the location: a location is used as a hash key, and a value
//! with interior mutability has no business being one — the hazard being that
//! a key's hash can change while it sits in a map. Composing a chain is a few
//! transform products, and chains are short; a caller that finds it hot can
//! memoize outside.

use og_core::{Arena, Key, OgResult, Tolerances, og_bail};
use og_math::{Transform, TransformKind};
use smallvec::SmallVec;

/// A rigid or similarity transform that placements are built from.
///
/// Shared: many locations refer to the same datum, and comparing two references
/// to it is what makes placement identity cheap.
pub type Datum = Transform;

/// A handle to a shared [`Datum`].
pub type DatumId = Key<Datum>;

/// The store of transforms that [`Location`] chains refer into.
///
/// One per document. Interning is not an optimisation here — it is what gives
/// placements a stable notion of sameness, since two chains naming the same
/// datum are known to agree without any floating-point comparison.
///
/// # Handles are relative to their store
///
/// A [`DatumId`] means nothing without the store that issued it, and this type
/// cannot tell a handle from another store apart from one of its own — two
/// stores that have each interned one transform both answer to the same first
/// handle. Generations catch a *stale* handle within one store; they cannot
/// catch a foreign one.
///
/// That is the documented trade of arena-based storage
/// (`docs/DATA_MODEL.md` §11), and it is why a document holds exactly one
/// store. Mixing handles between documents is a caller error that will not be
/// reported — it will resolve to whatever transform happens to sit at that
/// index.
#[derive(Debug, Clone, Default)]
pub struct DatumStore {
    arena: Arena<Datum>,
}

impl DatumStore {
    /// An empty store.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            arena: Arena::new(),
        }
    }

    /// Intern a transform, returning a handle to it.
    ///
    /// Identical transforms are *not* deduplicated: recognising two matrices as
    /// equal is a tolerance question, and the whole point of the chain
    /// representation is to avoid asking it. Callers that want sharing hold on
    /// to the handle.
    pub fn insert(&mut self, transform: Datum) -> DatumId {
        self.arena.insert(transform)
    }

    /// The transform behind `id`.
    #[must_use]
    pub fn get(&self, id: DatumId) -> Option<Datum> {
        self.arena.get(id).copied()
    }

    /// Number of interned transforms.
    #[must_use]
    pub fn len(&self) -> usize {
        self.arena.len()
    }

    /// Whether nothing has been interned.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.arena.is_empty()
    }

    /// Every datum, with its handle, in arena order.
    pub fn iter(&self) -> impl Iterator<Item = (DatumId, Datum)> {
        self.arena.iter().map(|(id, t)| (id, *t))
    }
}

/// A placement: a chain of `(datum, power)` pairs.
///
/// Applied left to right, so the first entry is the outermost transform. An
/// empty chain is the identity.
///
/// Powers are integers, so a datum applied twice costs one entry rather than
/// two, and its inverse is the same entry with the sign flipped.
/// Equality and hashing are structural: two locations agree when their chains
/// match entry for entry. Deliberately *not* a comparison of the composed
/// transforms — that would be a tolerance question, and the answer would depend
/// on rounding rather than on what the model says.
#[derive(Debug, Clone, Default, PartialEq, Eq, Hash)]
pub struct Location {
    chain: SmallVec<[(DatumId, i32); 2]>,
}

impl Location {
    /// The identity placement.
    #[must_use]
    pub fn identity() -> Self {
        Self::default()
    }

    /// A placement of one datum, applied once.
    #[must_use]
    pub fn of(datum: DatumId) -> Self {
        let mut chain = SmallVec::new();
        chain.push((datum, 1));
        Self { chain }
    }

    /// A placement of one datum, applied `power` times.
    ///
    /// A power of zero gives the identity; a negative power gives the inverse.
    #[must_use]
    pub fn powered(datum: DatumId, power: i32) -> Self {
        if power == 0 {
            return Self::identity();
        }
        let mut chain = SmallVec::new();
        chain.push((datum, power));
        Self { chain }
    }

    /// Whether this is the identity.
    #[must_use]
    pub fn is_identity(&self) -> bool {
        self.chain.is_empty()
    }

    /// The chain, outermost entry first.
    #[must_use]
    pub fn chain(&self) -> &[(DatumId, i32)] {
        &self.chain
    }

    /// Number of entries.
    #[must_use]
    pub fn depth(&self) -> usize {
        self.chain.len()
    }

    /// This placement followed by `inner`.
    ///
    /// `outer.then(inner)` applies `inner` first, then `outer` — the same order
    /// as transform composition, so a sub-shape's placement composed with its
    /// parent's reads the way the tree does.
    ///
    /// Adjacent entries naming the same datum are merged by adding their
    /// powers, and an entry whose power reaches zero is dropped. That keeps a
    /// chain from growing without bound as a placement is composed and undone
    /// repeatedly, and it is what makes `a.then(a.inverted())` come out exactly
    /// equal to the identity rather than merely close to it.
    #[must_use]
    pub fn then(&self, inner: &Self) -> Self {
        let mut chain = self.chain.clone();
        for &(datum, power) in &inner.chain {
            match chain.last_mut() {
                Some((last, last_power)) if *last == datum => {
                    *last_power += power;
                    if *last_power == 0 {
                        chain.pop();
                    }
                }
                _ => chain.push((datum, power)),
            }
        }
        Self { chain }
    }

    /// The inverse placement.
    ///
    /// Exact: the chain reverses and every power negates. No matrix is
    /// inverted, so `l.then(&l.inverted())` is the identity structurally, not
    /// approximately.
    #[must_use]
    pub fn inverted(&self) -> Self {
        Self {
            chain: self.chain.iter().rev().map(|&(d, p)| (d, -p)).collect(),
        }
    }

    /// The composed transform.
    ///
    /// # Errors
    ///
    /// [`OgError::Dangling`](og_core::OgError::Dangling) if the chain names a
    /// datum the store does not hold, and
    /// [`OgError::Numeric`](og_core::OgError::Numeric) if a negative power
    /// requires inverting a degenerate transform.
    pub fn composed(&self, store: &DatumStore) -> OgResult<Transform> {
        let mut result = Transform::IDENTITY;
        for &(id, power) in &self.chain {
            let Some(datum) = store.get(id) else {
                og_bail!(Dangling, "location refers to a datum not in this store");
            };
            let step = if power >= 0 { datum } else { datum.inverse()? };
            for _ in 0..power.unsigned_abs() {
                result = result * step;
            }
        }
        Ok(result)
    }

    /// Whether the composed transform preserves handedness.
    ///
    /// A shape placed by a handedness-reversing location has to have its
    /// orientation flipped to stay consistent, or a mirrored solid ends up
    /// inside out.
    ///
    /// # Errors
    ///
    /// As [`Location::composed`].
    pub fn preserves_handedness(&self, store: &DatumStore) -> OgResult<bool> {
        Ok(self.composed(store)?.preserves_handedness())
    }

    /// Whether two placements put a shape in the same position.
    ///
    /// Falls back to comparing the composed transforms, which costs more than
    /// [`PartialEq`] and answers a different question: two chains built by
    /// different routes can describe the same placement.
    ///
    /// # Errors
    ///
    /// As [`Location::composed`].
    pub fn is_same_placement(
        &self,
        other: &Self,
        store: &DatumStore,
        tol: Tolerances,
    ) -> OgResult<bool> {
        if self == other {
            return Ok(true);
        }
        Ok(self.composed(store)?.is_equal(&other.composed(store)?, tol))
    }

    /// The kind of the composed transform, for dispatch.
    ///
    /// # Errors
    ///
    /// As [`Location::composed`].
    pub fn kind(&self, store: &DatumStore) -> OgResult<TransformKind> {
        if self.is_identity() {
            return Ok(TransformKind::Identity);
        }
        Ok(self.composed(store)?.kind())
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use og_math::{Axis, Point, Vector};

    const T: Tolerances = Tolerances::millimetres();

    fn store() -> (DatumStore, DatumId, DatumId) {
        let mut s = DatumStore::new();
        let a = s.insert(Transform::translation(Vector::new(1.0, 0.0, 0.0)));
        let b = s.insert(Transform::rotation(Axis::Z, core::f64::consts::FRAC_PI_2));
        (s, a, b)
    }

    #[test]
    fn the_identity_is_empty_and_composes_to_nothing() {
        let (s, _, _) = store();
        let id = Location::identity();
        assert!(id.is_identity());
        assert_eq!(id.depth(), 0);
        assert_eq!(id.composed(&s).unwrap().kind(), TransformKind::Identity);
    }

    #[test]
    fn composition_applies_the_inner_placement_first() {
        let (s, a, b) = store();
        let outer = Location::of(a);
        let inner = Location::of(b);
        let combined = outer.then(&inner);

        let expected = outer.composed(&s).unwrap() * inner.composed(&s).unwrap();
        assert!(combined.composed(&s).unwrap().is_equal(&expected, T));

        // And the two orders differ, as they must.
        assert!(
            !combined
                .composed(&s)
                .unwrap()
                .is_equal(&inner.then(&outer).composed(&s).unwrap(), T)
        );
    }

    #[test]
    fn inversion_is_exact_rather_than_approximate() {
        // The property a matrix representation cannot offer: composing a
        // placement with its inverse gives the identity *structurally*, with no
        // residue to accumulate down an assembly tree.
        let (_, a, b) = store();
        let l = Location::of(a)
            .then(&Location::of(b))
            .then(&Location::of(a));
        let round_trip = l.then(&l.inverted());
        assert!(round_trip.is_identity(), "chain: {:?}", round_trip.chain());
        assert_eq!(round_trip, Location::identity());
    }

    #[test]
    fn repeated_composition_does_not_grow_the_chain() {
        // A placement applied a hundred times is one entry with a power of 100,
        // not a hundred entries — so an assembly that nests deeply stays cheap
        // to compare and to store.
        let (s, a, _) = store();
        let mut l = Location::identity();
        for _ in 0..100 {
            l = l.then(&Location::of(a));
        }
        assert_eq!(l.depth(), 1);
        assert_eq!(l.chain()[0].1, 100);
        assert!(
            l.composed(&s)
                .unwrap()
                .apply(Point::ORIGIN)
                .is_equal(Point::new(100.0, 0.0, 0.0), T)
        );
    }

    #[test]
    fn powers_cancel_exactly() {
        let (s, a, _) = store();
        let forward = Location::powered(a, 5);
        let back = Location::powered(a, -5);
        assert!(forward.then(&back).is_identity());
        assert_eq!(Location::powered(a, 0), Location::identity());
        assert!(
            Location::powered(a, -2)
                .composed(&s)
                .unwrap()
                .apply(Point::ORIGIN)
                .is_equal(Point::new(-2.0, 0.0, 0.0), T)
        );
    }

    #[test]
    fn equality_is_structural_not_numerical() {
        // Two chains that compose to the same transform are still different
        // placements. Asking whether they *are* the same is a comparison of
        // integers; asking whether they *land* in the same place is a separate,
        // costlier question with its own method.
        let mut s = DatumStore::new();
        let a = s.insert(Transform::translation(Vector::new(1.0, 0.0, 0.0)));
        let b = s.insert(Transform::translation(Vector::new(1.0, 0.0, 0.0)));

        let via_a = Location::of(a);
        let via_b = Location::of(b);
        assert_ne!(via_a, via_b, "different datums are different placements");
        assert!(via_a.is_same_placement(&via_b, &s, T).unwrap());

        assert_eq!(via_a, Location::of(a));
        assert!(via_a.is_same_placement(&Location::of(a), &s, T).unwrap());
    }

    #[test]
    fn locations_hash_consistently_with_equality() {
        use std::collections::HashSet;
        let (_, a, b) = store();
        let mut set = HashSet::new();
        set.insert(Location::of(a));
        set.insert(Location::of(a));
        set.insert(Location::of(b));
        set.insert(Location::identity());
        assert_eq!(set.len(), 3);
        assert!(set.contains(&Location::of(a)));
    }

    #[test]
    fn evaluating_a_location_does_not_change_its_identity() {
        // A location is used as a hash key, so nothing it does may alter how it
        // compares or hashes. It holds no interior mutability at all, which is
        // what makes that guarantee rather than a hope.
        let (s, a, _) = store();
        let x = Location::of(a);
        let y = Location::of(a);
        let _ = x.composed(&s).unwrap();
        assert_eq!(x, y);

        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        let hash = |l: &Location| {
            let mut h = DefaultHasher::new();
            l.hash(&mut h);
            h.finish()
        };
        assert_eq!(hash(&x), hash(&y));
    }

    #[test]
    fn a_location_is_send_and_sync() {
        // Needed for the parallel algorithms further up the stack, and easy to
        // lose to a cache tucked inside the type.
        const fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<Location>();
        assert_send_sync::<DatumStore>();
    }

    #[test]
    fn a_handle_the_store_does_not_hold_is_reported_rather_than_ignored() {
        let (s, _, _) = store();
        let mut other = DatumStore::new();
        // Past the end of `s`, so it genuinely does not resolve.
        let mut beyond = other.insert(Transform::translation(Vector::Y));
        for _ in 0..10 {
            beyond = other.insert(Transform::translation(Vector::Y));
        }
        assert!(Location::of(beyond).composed(&s).is_err());
    }

    #[test]
    fn a_handle_from_another_store_at_the_same_index_is_not_detectable() {
        // Not a defect to fix but a boundary to know: an arena handle is
        // meaningful only relative to the arena that issued it. Generations
        // catch a stale handle within one store; nothing here can catch a
        // foreign one, and it will silently resolve to whatever sits at that
        // index. This is why a document holds exactly one store.
        let (s, first, _) = store();
        let mut other = DatumStore::new();
        let foreign = other.insert(Transform::translation(Vector::new(0.0, 99.0, 0.0)));

        let resolved = Location::of(foreign).composed(&s).unwrap();
        let ours = Location::of(first).composed(&s).unwrap();
        assert!(
            resolved.is_equal(&ours, T),
            "the foreign handle resolved to our own first datum"
        );
    }

    #[test]
    fn handedness_follows_the_composed_transform() {
        let mut s = DatumStore::new();
        let mirror = s.insert(Transform::plane_mirror(
            Point::ORIGIN,
            og_math::Direction::Z,
        ));
        let once = Location::of(mirror);
        assert!(!once.preserves_handedness(&s).unwrap());
        // Two mirrors make a rotation.
        assert!(
            Location::powered(mirror, 2)
                .preserves_handedness(&s)
                .unwrap()
        );
        assert!(Location::identity().preserves_handedness(&s).unwrap());
    }

    #[test]
    fn composing_different_datums_keeps_both_entries() {
        let (_, a, b) = store();
        let l = Location::of(a).then(&Location::of(b));
        assert_eq!(l.depth(), 2);
        assert_eq!(l.chain(), &[(a, 1), (b, 1)]);
    }

    #[test]
    fn the_datum_store_does_not_deduplicate() {
        // Deduplicating would mean deciding two matrices are equal, which is a
        // tolerance question — precisely the one the chain exists to avoid.
        let mut s = DatumStore::new();
        let t = Transform::translation(Vector::X);
        let a = s.insert(t);
        let b = s.insert(t);
        assert_ne!(a, b);
        assert_eq!(s.len(), 2);
        assert!(!s.is_empty());
    }
}
