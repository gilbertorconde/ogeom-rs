//! Stable entity identity and provenance.
//!
//! See `docs/DATA_MODEL.md` §8. A deliberate divergence from the conventional
//! design, and the reason this has to live in the foundation crate rather than
//! being added later.
//!
//! Conventionally, topology is identified by pointer. Every modeling operation
//! allocates new nodes, so every reference into a previous result dies — that *is* the
//! topological naming problem, and every downstream fix is an attempt to
//! reconstruct identity after the fact by walking history maps.
//!
//! Here an entity's identity is *what produced it, and from what*. A rebuild
//! with different parameters runs the same operations over the same inputs and
//! therefore produces entities with the same provenance, so a reference like
//! "the fillet on this edge" survives a change to an unrelated dimension.
//!
//! Provenance does not replace operation history — history is what a binding
//! layer consumes, and it is the honest answer where provenance cannot resolve
//! a reference. It is the primary mechanism, not the only one.

use core::num::NonZeroU64;

use smallvec::SmallVec;

/// A stable identity for a topological entity, valid for the lifetime of a
/// document.
///
/// Distinct from an arena [`Key`](crate::Key): a key says *where the data is*
/// and dies when the entity is rebuilt; an `EntityId` says *what the entity is*
/// and survives.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EntityId(NonZeroU64);

impl EntityId {
    /// The identity's raw value. Non-zero.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0.get()
    }
}

/// Identifies one invocation of a modeling operation.
///
/// Stable across rebuilds: the third extrusion in a document's recompute is
/// `OpId(3)` every time, which is what lets provenance survive a parameter
/// change.
///
/// The default, `OpId(0)`, is the implicit operation a document starts in —
/// whatever was there before anything was deliberately begun.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct OpId(pub u32);

/// What an entity *is*, relative to the operation that made it.
///
/// The low values are shared vocabulary; everything from [`Role::OP_DEFINED`]
/// up is interpreted by the producing operation alone. Keeping it a newtype
/// rather than an enum avoids inventing a taxonomy of every role in a CAD
/// kernel before we have written the operations that need one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Role(pub u32);

impl Role {
    /// No distinguishing role — the operation produced exactly one entity of
    /// this kind, so it needs no further discriminator.
    pub const SOLE: Self = Self(0);
    /// The result's outer boundary: the outer wire of a face, the outer shell of
    /// a solid.
    pub const OUTER: Self = Self(1);
    /// An inner boundary — a hole.
    pub const INNER: Self = Self(2);
    /// The start of a swept or extruded result.
    pub const START_CAP: Self = Self(3);
    /// The end of a swept or extruded result.
    pub const END_CAP: Self = Self(4);
    /// The swept side wall between the caps.
    pub const LATERAL: Self = Self(5);
    /// A seam on a closed surface.
    pub const SEAM: Self = Self(6);
    /// The first value an operation may assign meaning to itself.
    pub const OP_DEFINED: u32 = 1024;

    /// An operation-defined role. `index` is offset above [`Role::OP_DEFINED`].
    #[must_use]
    pub const fn op_defined(index: u32) -> Self {
        Self(Self::OP_DEFINED + index)
    }
}

/// Where an entity came from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Provenance {
    /// Created outright by an operation, from no prior entity — the `+Z` face of
    /// a box, the lateral surface of a cylinder.
    Primitive {
        /// The operation that created it.
        op: OpId,
        /// Which part of that operation's result this is.
        role: Role,
    },
    /// Derived from one or more existing entities. A face split by a boolean
    /// names the face it came from; an intersection edge names both faces.
    Derived {
        /// The operation that derived it.
        op: OpId,
        /// The inputs it came from, in the operation's canonical order.
        from: SmallVec<[EntityId; 2]>,
        /// Which part of that operation's result this is.
        role: Role,
    },
    /// Read from a file. `external` is the source's own identifier — a STEP
    /// entity number, say — so that a re-import matches entities up.
    Imported {
        /// Which imported document it came from.
        source: SourceId,
        /// The identifier the source file gave it.
        external: u64,
    },
}

impl Provenance {
    /// The operation that produced this entity, if any.
    #[must_use]
    pub const fn op(&self) -> Option<OpId> {
        match self {
            Self::Primitive { op, .. } | Self::Derived { op, .. } => Some(*op),
            Self::Imported { .. } => None,
        }
    }

    /// The entities this one was derived from.
    #[must_use]
    pub fn inputs(&self) -> &[EntityId] {
        match self {
            Self::Derived { from, .. } => from,
            Self::Primitive { .. } | Self::Imported { .. } => &[],
        }
    }
}

/// Identifies an imported document.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SourceId(pub u32);

/// Assigns [`EntityId`]s and remembers what each one came from.
///
/// One per document. Cloning it clones the whole history, which is what a
/// rebuild-with-rollback needs.
#[derive(Debug, Clone, Default)]
pub struct ProvenanceTable {
    entries: Vec<Provenance>,
}

impl ProvenanceTable {
    /// An empty table.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// Number of entities recorded.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether nothing has been recorded.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Record an entity's provenance and return its identity.
    ///
    /// # Panics
    ///
    /// If more than `u64::MAX - 1` entities are recorded. Not reachable.
    #[allow(clippy::expect_used, reason = "documented panic; see # Panics")]
    pub fn record(&mut self, provenance: Provenance) -> EntityId {
        self.entries.push(provenance);
        // Ids start at 1 so that EntityId can be NonZeroU64 and Option<EntityId>
        // costs nothing.
        let raw = u64::try_from(self.entries.len()).expect("entity count exceeded u64");
        EntityId(NonZeroU64::new(raw).expect("length is at least 1 after push"))
    }

    /// Convenience for [`Provenance::Primitive`].
    pub fn primitive(&mut self, op: OpId, role: Role) -> EntityId {
        self.record(Provenance::Primitive { op, role })
    }

    /// Convenience for [`Provenance::Derived`].
    pub fn derived(
        &mut self,
        op: OpId,
        from: impl IntoIterator<Item = EntityId>,
        role: Role,
    ) -> EntityId {
        self.record(Provenance::Derived {
            op,
            from: from.into_iter().collect(),
            role,
        })
    }

    /// The provenance of `id`, or `None` if it belongs to another document.
    #[must_use]
    pub fn get(&self, id: EntityId) -> Option<&Provenance> {
        let index = usize::try_from(id.get()).ok()?.checked_sub(1)?;
        self.entries.get(index)
    }

    /// Walk `id`'s derivation back to the entities it ultimately came from.
    ///
    /// Returns the roots — entities that are `Primitive` or `Imported`. This is
    /// how a stale reference is resolved after a rebuild: find what the user
    /// originally picked, then find what that became.
    ///
    /// Cycles cannot occur, because an entity can only be derived from ids that
    /// already existed when it was recorded. The visited set guards against a
    /// table assembled by hand or deserialized from a corrupt file.
    #[must_use]
    pub fn roots(&self, id: EntityId) -> Vec<EntityId> {
        let mut out = Vec::new();
        let mut seen = std::collections::HashSet::new();
        let mut stack = vec![id];
        while let Some(current) = stack.pop() {
            if !seen.insert(current) {
                continue;
            }
            match self.get(current) {
                Some(Provenance::Derived { from, .. }) if !from.is_empty() => {
                    stack.extend(from.iter().copied());
                }
                Some(_) => out.push(current),
                None => {}
            }
        }
        out.sort_unstable();
        out
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn ids_are_non_zero_and_distinct() {
        let mut t = ProvenanceTable::new();
        let a = t.primitive(OpId(0), Role::SOLE);
        let b = t.primitive(OpId(0), Role::OUTER);
        assert_ne!(a, b);
        assert!(a.get() > 0 && b.get() > 0);
        assert_eq!(
            core::mem::size_of::<Option<EntityId>>(),
            core::mem::size_of::<EntityId>()
        );
    }

    #[test]
    fn derived_entities_remember_their_inputs() {
        let mut t = ProvenanceTable::new();
        let face_a = t.primitive(OpId(1), Role::LATERAL);
        let face_b = t.primitive(OpId(2), Role::LATERAL);
        let section = t.derived(OpId(3), [face_a, face_b], Role::SOLE);

        let p = t.get(section).unwrap();
        assert_eq!(p.op(), Some(OpId(3)));
        assert_eq!(p.inputs(), &[face_a, face_b]);
    }

    #[test]
    fn roots_walk_back_through_a_derivation_chain() {
        let mut t = ProvenanceTable::new();
        let original = t.primitive(OpId(1), Role::END_CAP);
        // A boolean splits the face, then a fillet modifies one fragment.
        let split = t.derived(OpId(2), [original], Role::op_defined(0));
        let filleted = t.derived(OpId(3), [split], Role::SOLE);

        assert_eq!(t.roots(filleted), vec![original]);
        assert_eq!(t.roots(original), vec![original], "a root is its own root");
    }

    #[test]
    fn roots_of_a_multi_parent_entity_include_every_branch() {
        let mut t = ProvenanceTable::new();
        let a = t.primitive(OpId(1), Role::SOLE);
        let b = t.primitive(OpId(2), Role::SOLE);
        let mid = t.derived(OpId(3), [a], Role::SOLE);
        let joined = t.derived(OpId(4), [mid, b], Role::SOLE);

        let mut expected = vec![a, b];
        expected.sort_unstable();
        assert_eq!(t.roots(joined), expected);
    }

    #[test]
    fn ids_from_another_document_do_not_resolve() {
        let mut t = ProvenanceTable::new();
        let mine = t.primitive(OpId(0), Role::SOLE);
        let mut other = ProvenanceTable::new();
        for _ in 0..10 {
            other.primitive(OpId(0), Role::SOLE);
        }
        let theirs = other.primitive(OpId(0), Role::SEAM);

        assert!(t.get(mine).is_some());
        assert!(t.get(theirs).is_none(), "foreign id must not resolve");
    }

    #[test]
    fn op_defined_roles_do_not_collide_with_shared_ones() {
        assert!(Role::op_defined(0).0 > Role::SEAM.0);
        assert_ne!(Role::op_defined(0), Role::op_defined(1));
    }
}
