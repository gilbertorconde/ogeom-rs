//! Typed generational arenas.
//!
//! Topology lives in arenas rather than behind reference counting — see
//! `docs/DATA_MODEL.md` §11. Keys are small, `Copy`, comparable and hashable,
//! which is what makes stable entity identity possible at all.
//!
//! Slots are generational: freeing a slot bumps its generation, so a stale key
//! fails to resolve instead of silently aliasing whatever was allocated there
//! next. That failure mode is worth eight bytes per key in a kernel where the
//! alternative is a wrong answer rather than a crash.
//!
//! # Keys are scoped to the arena that issued them
//!
//! A generation catches a key that has outlived its slot. It cannot catch a key
//! from a *different* arena, because index 3 generation 0 means something in
//! every arena — so a handle from one document resolved against another comes
//! back with whatever sits at that index, and answers confidently about the
//! wrong entity. Nothing about the result says so.
//!
//! Every arena therefore takes an identifier the first time something is put
//! in it, every key it issues carries that identifier, and every lookup
//! compares it. A foreign key resolves to `None`, exactly as a stale one does.
//! The cost is four bytes per key and one comparison per lookup, against a
//! whole class of silent wrong answers.
//!
//! Cloning an arena keeps its identifier, because a clone is the same document
//! and handles into it should keep working. Identifiers are per-process and are
//! never serialized: a document read back from a file is a new arena with a new
//! identifier, and the reader re-stamps the handles it read.

use core::fmt;
use core::hash::{Hash, Hasher};
use core::marker::PhantomData;
use core::sync::atomic::{AtomicU32, Ordering};

/// Hands out arena identifiers.
///
/// Starts at one so that zero can mean *unscoped* — the state of a key built
/// by a deserializer that does not yet know which arena it will belong to.
static NEXT_SCOPE: AtomicU32 = AtomicU32::new(1);

/// An identifier that no arena in this process shares.
fn next_scope() -> u32 {
    NEXT_SCOPE.fetch_add(1, Ordering::Relaxed)
}

/// The identifier a key carries before it has been bound to an arena.
///
/// A key with this scope resolves in no arena at all. That is deliberate: a
/// handle read from a file is meaningless until the reader says which document
/// it belongs to.
pub const UNSCOPED: u32 = 0;

/// A handle into an [`Arena<T>`].
///
/// Phantom-typed, so a `Key<Face>` cannot be used to index an `Arena<Edge>`.
/// The marker is `fn() -> T` so the key stays `Copy`, `Send` and `Sync`
/// regardless of `T`.
pub struct Key<T> {
    index: u32,
    generation: u32,
    scope: u32,
    marker: PhantomData<fn() -> T>,
}

impl<T> Key<T> {
    const fn new(index: u32, generation: u32, scope: u32) -> Self {
        Self {
            index,
            generation,
            scope,
            marker: PhantomData,
        }
    }

    /// Which arena issued this key.
    ///
    /// [`UNSCOPED`] for a key that has not been bound to one.
    #[must_use]
    pub const fn scope(self) -> u32 {
        self.scope
    }

    /// This key, bound to the arena with the given identifier.
    ///
    /// For a deserializer, which rebuilds handles before it has an arena to
    /// bind them to. Nothing else should need it: a key that came from an arena
    /// already names the right one, and moving a key between arenas is the
    /// mistake the scope exists to catch.
    #[must_use]
    pub const fn with_scope(self, scope: u32) -> Self {
        Self { scope, ..self }
    }

    /// Position of the slot this key refers to.
    #[must_use]
    pub const fn index(self) -> u32 {
        self.index
    }

    /// Generation stamp, used to detect a key outliving its slot.
    #[must_use]
    pub const fn generation(self) -> u32 {
        self.generation
    }

    /// A key naming a given slot, for reading a document back from a file.
    ///
    /// Deliberately narrow. Forging a handle is precisely what generations
    /// exist to prevent, and [`Arena::insert`] is what issues one within a
    /// process. But a file records the handles a document was written with, and
    /// a reader that could not rebuild them would have to renumber everything —
    /// which is to say, hand back a different document.
    ///
    /// A key made this way is not trusted: it resolves through [`Arena::get`]
    /// like any other, so a stale or out-of-range one comes back `None` rather
    /// than aliasing whatever sits at that index.
    #[must_use]
    pub const fn from_parts(index: u32, generation: u32) -> Self {
        Self::new(index, generation, UNSCOPED)
    }
}

// Derived impls would demand `T: Clone` and friends; the key holds no `T`.
impl<T> Clone for Key<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Key<T> {}
impl<T> PartialEq for Key<T> {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index
            && self.generation == other.generation
            && self.scope == other.scope
    }
}
impl<T> Eq for Key<T> {}
impl<T> Hash for Key<T> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.index.hash(state);
        self.generation.hash(state);
        self.scope.hash(state);
    }
}
impl<T> PartialOrd for Key<T> {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl<T> Ord for Key<T> {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        (self.scope, self.index, self.generation).cmp(&(other.scope, other.index, other.generation))
    }
}
impl<T> fmt::Debug for Key<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Key({}v{}@{})", self.index, self.generation, self.scope)
    }
}

#[derive(Debug, Clone)]
enum Slot<T> {
    Occupied {
        generation: u32,
        value: T,
    },
    Vacant {
        generation: u32,
        next_free: Option<u32>,
    },
}

/// A generational arena of `T`.
#[derive(Debug, Clone)]
pub struct Arena<T> {
    slots: Vec<Slot<T>>,
    free_head: Option<u32>,
    len: usize,
    /// Which arena this is. [`UNSCOPED`] until the first insert, because
    /// `new` is `const` and a counter cannot be read from one — and an arena
    /// with nothing in it has issued no keys to disagree with.
    scope: u32,
}

impl<T> Default for Arena<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Arena<T> {
    /// An empty arena.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            slots: Vec::new(),
            free_head: None,
            len: 0,
            scope: UNSCOPED,
        }
    }

    /// An empty arena with room for `capacity` entries.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            slots: Vec::with_capacity(capacity),
            free_head: None,
            len: 0,
            scope: UNSCOPED,
        }
    }

    /// Number of live entries.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.len
    }

    /// Whether there are no live entries.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Which arena this is, for stamping keys that were rebuilt elsewhere.
    ///
    /// [`UNSCOPED`] until the first insert.
    #[must_use]
    pub const fn scope(&self) -> u32 {
        self.scope
    }

    /// Whether a key was issued by this arena.
    ///
    /// Distinct from [`Arena::contains`], which also asks whether the slot is
    /// still live. This asks only whether the key belongs here at all, which is
    /// the question a caller wants when reporting *why* a lookup failed.
    #[must_use]
    pub const fn issued(&self, key: Key<T>) -> bool {
        key.scope == self.scope
    }

    /// Insert a value, returning its key.
    ///
    /// The first insert is what fixes the arena's identity, since [`Arena::new`]
    /// is `const` and cannot read a counter. That is safe because an arena with
    /// nothing in it has issued no keys to disagree with.
    ///
    /// # Panics
    ///
    /// If the arena exceeds `u32::MAX` slots. A single model reaching four
    /// billion topological entities is a bug elsewhere, not a case to handle.
    #[allow(clippy::expect_used, reason = "documented panic; see # Panics")]
    pub fn insert(&mut self, value: T) -> Key<T> {
        if self.scope == UNSCOPED {
            self.scope = next_scope();
        }
        self.len += 1;
        match self.free_head {
            Some(index) => {
                let idx = index as usize;
                let (generation, next_free) = match &self.slots[idx] {
                    Slot::Vacant {
                        generation,
                        next_free,
                    } => (*generation, *next_free),
                    Slot::Occupied { .. } => unreachable!("free list pointed at an occupied slot"),
                };
                self.free_head = next_free;
                self.slots[idx] = Slot::Occupied { generation, value };
                Key::new(index, generation, self.scope)
            }
            None => {
                let index = u32::try_from(self.slots.len()).expect("arena exceeded u32::MAX slots");
                self.slots.push(Slot::Occupied {
                    generation: 0,
                    value,
                });
                Key::new(index, 0, self.scope)
            }
        }
    }

    /// Borrow the value behind `key`, or `None` if the key is stale.
    #[must_use]
    pub fn get(&self, key: Key<T>) -> Option<&T> {
        if key.scope != self.scope {
            return None;
        }
        match self.slots.get(key.index as usize)? {
            Slot::Occupied { generation, value } if *generation == key.generation => Some(value),
            _ => None,
        }
    }

    /// Mutably borrow the value behind `key`, or `None` if the key is stale.
    pub fn get_mut(&mut self, key: Key<T>) -> Option<&mut T> {
        if key.scope != self.scope {
            return None;
        }
        match self.slots.get_mut(key.index as usize)? {
            Slot::Occupied { generation, value } if *generation == key.generation => Some(value),
            _ => None,
        }
    }

    /// Whether `key` resolves to a live entry.
    #[must_use]
    pub fn contains(&self, key: Key<T>) -> bool {
        self.get(key).is_some()
    }

    /// Remove and return the value behind `key`, if it is live.
    ///
    /// The slot's generation is bumped, invalidating every outstanding copy of
    /// `key`.
    pub fn remove(&mut self, key: Key<T>) -> Option<T> {
        if key.scope != self.scope {
            return None;
        }
        let slot = self.slots.get_mut(key.index as usize)?;
        let generation = match slot {
            Slot::Occupied { generation, .. } if *generation == key.generation => *generation,
            _ => return None,
        };
        // Saturating rather than wrapping: a slot recycled 4 billion times stops
        // being reusable, which is strictly better than handing out a generation
        // that collides with a key someone still holds.
        let next = generation.saturating_add(1);
        let replaced = core::mem::replace(
            slot,
            Slot::Vacant {
                generation: next,
                next_free: self.free_head,
            },
        );
        if next != u32::MAX {
            self.free_head = Some(key.index);
        }
        self.len -= 1;
        match replaced {
            Slot::Occupied { value, .. } => Some(value),
            Slot::Vacant { .. } => None,
        }
    }

    /// Iterate over live `(key, &value)` pairs, in slot order.
    pub fn iter(&self) -> impl Iterator<Item = (Key<T>, &T)> {
        let scope = self.scope;
        self.slots
            .iter()
            .enumerate()
            .filter_map(move |(i, slot)| match slot {
                Slot::Occupied { generation, value } => {
                    // `insert` refuses to grow past u32::MAX, so this cannot truncate.
                    #[allow(clippy::cast_possible_truncation)]
                    Some((Key::new(i as u32, *generation, scope), value))
                }
                Slot::Vacant { .. } => None,
            })
    }

    /// Iterate over live `(key, &mut value)` pairs, in slot order.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = (Key<T>, &mut T)> {
        let scope = self.scope;
        self.slots
            .iter_mut()
            .enumerate()
            .filter_map(move |(i, slot)| match slot {
                Slot::Occupied { generation, value } =>
                {
                    #[allow(clippy::cast_possible_truncation)]
                    Some((Key::new(i as u32, *generation, scope), value))
                }
                Slot::Vacant { .. } => None,
            })
    }

    /// Iterate over live values.
    pub fn values(&self) -> impl Iterator<Item = &T> {
        self.iter().map(|(_, v)| v)
    }

    /// Consume the arena, yielding its live values in index order.
    ///
    /// For appending one arena's contents onto another: the receiving arena
    /// hands out its own keys, so the values travel bare.
    pub fn into_values(self) -> impl Iterator<Item = T> {
        self.slots.into_iter().filter_map(|slot| match slot {
            Slot::Occupied { value, .. } => Some(value),
            Slot::Vacant { .. } => None,
        })
    }

    /// Whether the arena has only ever been appended to: every slot occupied,
    /// every generation zero.
    ///
    /// When this holds, [`Arena::len`] is also the next index [`Arena::insert`]
    /// will hand out — the precondition for extending the arena by offset,
    /// where a caller predicts the keys of entries it is about to append.
    #[must_use]
    pub fn is_dense(&self) -> bool {
        self.len == self.slots.len()
            && self
                .slots
                .iter()
                .all(|slot| matches!(slot, Slot::Occupied { generation: 0, .. }))
    }

    /// Remove every entry, bumping all generations so existing keys go stale.
    pub fn clear(&mut self) {
        let keys: Vec<_> = self.iter().map(|(k, _)| k).collect();
        for key in keys {
            self.remove(key);
        }
    }
}

impl<T> core::ops::Index<Key<T>> for Arena<T> {
    type Output = T;

    /// # Panics
    ///
    /// If the key is stale. Use [`Arena::get`] where that is a possibility.
    #[allow(
        clippy::expect_used,
        reason = "Index cannot return Result; see # Panics"
    )]
    fn index(&self, key: Key<T>) -> &T {
        self.get(key).expect("stale arena key")
    }
}

impl<T> core::ops::IndexMut<Key<T>> for Arena<T> {
    /// # Panics
    ///
    /// If the key is stale. Use [`Arena::get_mut`] where that is a possibility.
    #[allow(
        clippy::expect_used,
        reason = "IndexMut cannot return Result; see # Panics"
    )]
    fn index_mut(&mut self, key: Key<T>) -> &mut T {
        self.get_mut(key).expect("stale arena key")
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn insert_and_get() {
        let mut a = Arena::new();
        let k1 = a.insert("one");
        let k2 = a.insert("two");
        assert_eq!(a.get(k1), Some(&"one"));
        assert_eq!(a.get(k2), Some(&"two"));
        assert_eq!(a.len(), 2);
    }

    #[test]
    fn removed_key_goes_stale_and_does_not_alias() {
        let mut a = Arena::new();
        let old = a.insert(1_u32);
        assert_eq!(a.remove(old), Some(1));

        // The slot is reused, but the old key must not resolve to the new value.
        let new = a.insert(2_u32);
        assert_eq!(new.index(), old.index(), "slot should have been reused");
        assert_eq!(a.get(new), Some(&2));
        assert_eq!(a.get(old), None, "stale key aliased a live entry");
        assert!(!a.contains(old));
    }

    #[test]
    fn double_remove_is_none() {
        let mut a = Arena::new();
        let k = a.insert(7_u8);
        assert_eq!(a.remove(k), Some(7));
        assert_eq!(a.remove(k), None);
        assert_eq!(a.len(), 0);
    }

    #[test]
    fn iteration_skips_holes() {
        let mut a = Arena::new();
        let keys: Vec<_> = (0..5_u32).map(|i| a.insert(i)).collect();
        a.remove(keys[1]);
        a.remove(keys[3]);
        let live: Vec<_> = a.values().copied().collect();
        assert_eq!(live, vec![0, 2, 4]);
        assert_eq!(a.len(), 3);
    }

    #[test]
    fn into_values_yields_values_in_index_order() {
        let mut a = Arena::new();
        for i in 0..5_u32 {
            a.insert(i * 10);
        }
        let values: Vec<_> = a.into_values().collect();
        assert_eq!(values, vec![0, 10, 20, 30, 40]);
    }

    #[test]
    fn an_arena_that_never_removed_is_dense() {
        let mut a = Arena::new();
        assert!(a.is_dense(), "an empty arena has no holes");
        for i in 0..4_u32 {
            a.insert(i);
        }
        assert!(a.is_dense());
    }

    #[test]
    fn a_removal_makes_an_arena_not_dense() {
        let mut a = Arena::new();
        let keys: Vec<_> = (0..3_u32).map(|i| a.insert(i)).collect();
        a.remove(keys[1]);
        assert!(!a.is_dense(), "a vacant slot is a hole");

        // Refilling the slot does not restore density either: the recycled
        // entry sits at a bumped generation, so `len` no longer predicts the
        // keys of future appends alone.
        a.insert(9);
        assert!(!a.is_dense(), "a recycled slot is off generation zero");
    }

    #[test]
    fn clear_invalidates_every_key() {
        let mut a = Arena::new();
        let keys: Vec<_> = (0..4_u32).map(|i| a.insert(i)).collect();
        a.clear();
        assert!(a.is_empty());
        assert!(keys.iter().all(|&k| a.get(k).is_none()));
    }

    #[test]
    fn keys_are_hashable_and_distinct() {
        use std::collections::HashSet;
        let mut a = Arena::new();
        let set: HashSet<_> = (0..64_u32).map(|i| a.insert(i)).collect();
        assert_eq!(set.len(), 64);
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod scope_tests {
    use super::*;

    #[test]
    fn a_key_from_one_arena_does_not_resolve_in_another() {
        // The whole reason the scope exists. Index 0 generation 0 means
        // something in every arena, so without it this lookup succeeds and
        // answers about the wrong value — confidently, with nothing about the
        // result to say so.
        let mut a: Arena<&str> = Arena::new();
        let mut b: Arena<&str> = Arena::new();
        let here = a.insert("in a");
        let there = b.insert("in b");

        assert_eq!(a.get(here), Some(&"in a"));
        assert_eq!(b.get(there), Some(&"in b"));
        assert_eq!(here.index(), there.index(), "the same slot in both");
        assert_eq!(a.get(there), None, "a foreign key must not resolve");
        assert_eq!(b.get(here), None);
        assert!(!a.issued(there));
    }

    #[test]
    fn foreign_keys_are_refused_by_every_route_in() {
        let mut a: Arena<u32> = Arena::new();
        let mut b: Arena<u32> = Arena::new();
        let key = a.insert(1);
        b.insert(2);

        assert!(!b.contains(key));
        assert_eq!(b.get_mut(key), None);
        assert_eq!(b.remove(key), None, "and it must not remove something else");
        assert_eq!(b.len(), 1, "nothing was taken out");
    }

    #[test]
    fn keys_from_different_arenas_are_not_equal_and_do_not_collide() {
        // Equality and hashing have to agree with resolution, or a map keyed on
        // handles merges entries from two documents.
        use std::collections::HashSet;
        let mut a: Arena<u32> = Arena::new();
        let mut b: Arena<u32> = Arena::new();
        let here = a.insert(1);
        let there = b.insert(2);

        assert_ne!(here, there);
        let mut set = HashSet::new();
        set.insert(here);
        set.insert(there);
        assert_eq!(set.len(), 2, "two documents' handles collided in a map");
    }

    #[test]
    fn an_unscoped_key_resolves_nowhere_until_it_is_bound() {
        // What a deserializer builds. It names a slot but no arena, and a
        // handle that names no arena is meaningless until someone says which.
        let mut a: Arena<u32> = Arena::new();
        let real = a.insert(7);
        let loose: Key<u32> = Key::from_parts(real.index(), real.generation());

        assert_eq!(loose.scope(), UNSCOPED);
        assert_eq!(a.get(loose), None);
        assert_eq!(a.get(loose.with_scope(a.scope())), Some(&7));
    }

    #[test]
    fn a_clone_answers_to_the_originals_handles() {
        // A clone is the same document — a snapshot — so handles into it keep
        // working. If it took a fresh identifier, every handle a caller held
        // would silently stop resolving after a clone.
        let mut a: Arena<u32> = Arena::new();
        let key = a.insert(5);
        let copy = a.clone();
        assert_eq!(copy.get(key), Some(&5));
    }

    #[test]
    fn an_empty_arena_has_issued_nothing_to_disagree_with() {
        // `new` is const, so the identifier cannot be taken until the first
        // insert. That is safe precisely because an arena with nothing in it
        // has handed out no keys.
        let empty: Arena<u32> = Arena::new();
        assert_eq!(empty.scope(), UNSCOPED);
        let mut used: Arena<u32> = Arena::new();
        used.insert(1);
        assert_ne!(used.scope(), UNSCOPED);
    }
}
