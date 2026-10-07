//! Typed generational arenas.
//!
//! Topology lives in arenas rather than behind reference counting: see
//! `docs/DATA_MODEL.md` §11. Keys are small, `Copy`, comparable and hashable,
//! which is what makes stable entity identity possible at all.
//!
//! Slots are generational: freeing a slot bumps its generation, so a stale key
//! fails to resolve instead of silently aliasing whatever was allocated there
//! next. That failure mode is worth eight bytes per key in a kernel where the
//! alternative is a wrong answer rather than a crash. A slot dropped by
//! [`Arena::retain`] is never handed out again at all.
//!
//! # Keys are scoped to the arena that issued them
//!
//! A generation catches a key that has outlived its slot. It cannot catch a key
//! from a *different* arena, because index 3 generation 0 means something in
//! every arena, so a handle from one document resolved against another comes
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

use hashbrown::HashMap;

/// Hands out arena identifiers.
///
/// Starts at one so that zero can mean *unscoped*: the state of a key built
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
    /// a reader that could not rebuild them would have to renumber everything,
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

// Derived impls would demand `T: Clone` and friends. The key holds no `T`.
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

/// A live value, with the slot it sits in and the generation it was stored
/// at.
#[derive(Debug, Clone)]
struct Entry<T> {
    slot: u32,
    generation: u32,
    value: T,
}

/// A generational arena of `T`.
///
/// Live values are stored packed, in two runs. The *tail* holds every slot
/// from some index on, each one occupied, so its values are found by
/// arithmetic; an arena that has only been appended to is all tail. The
/// values before it are *scattered*: slots whose neighbours were dropped,
/// found through a map from slot to position. An arena holds what is live
/// and nothing for the slots it no longer uses, so cloning one costs what
/// it holds, however many slots it has handed out.
#[derive(Debug, Clone)]
pub struct Arena<T> {
    /// The scattered values, then the tail's, in slot order unless
    /// `ordered` says otherwise.
    entries: Vec<Entry<T>>,
    /// The first slot of the tail.
    tail: u32,
    /// Where the tail starts in `entries`: the number of scattered values.
    tail_at: u32,
    /// Slot to position in `entries`, for the scattered values.
    scattered: Option<HashMap<u32, u32>>,
    /// Whether `entries` runs in slot order. Removing a value and refilling
    /// a freed slot can break the order; iterating mutably restores it.
    ordered: bool,
    /// Freed slots for `insert` to reuse, each with the generation it takes
    /// next. A slot dropped by [`Arena::retain`] is never reused.
    free: Vec<(u32, u32)>,
    /// How many slots have been handed out: the slot a fresh insert takes.
    next_slot: u32,
    /// Which arena this is. [`UNSCOPED`] until the first insert, because
    /// `new` is `const` and a counter cannot be read from one, and an arena
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
            entries: Vec::new(),
            tail: 0,
            tail_at: 0,
            scattered: None,
            ordered: true,
            free: Vec::new(),
            next_slot: 0,
            scope: UNSCOPED,
        }
    }

    /// An empty arena with room for `capacity` entries.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            entries: Vec::with_capacity(capacity),
            ..Self::new()
        }
    }

    /// Number of live entries.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether there are no live entries.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.entries.is_empty()
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

    /// The slot the next [`Arena::insert`] takes, at generation zero, when
    /// it takes a fresh one; `None` while a freed slot waits to be reused.
    ///
    /// The precondition for extending the arena by offset, where a caller
    /// predicts the keys of entries it is about to append: the `k`-th
    /// insert from here lands at this index plus `k`.
    #[must_use]
    pub fn next_index(&self) -> Option<u32> {
        self.free.is_empty().then_some(self.next_slot)
    }

    /// The live key at slot `index`, if the slot holds a value.
    #[must_use]
    pub fn key_at(&self, index: u32) -> Option<Key<T>> {
        let entry = self.entries.get(self.position(index)?)?;
        (entry.slot == index).then(|| Key::new(index, entry.generation, self.scope))
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
        if let Some((slot, generation)) = self.free.pop() {
            self.scatter_tail();
            let at = u32::try_from(self.entries.len()).expect("arena exceeded u32::MAX slots");
            self.ordered &= self.entries.last().is_none_or(|last| last.slot < slot);
            self.entries.push(Entry {
                slot,
                generation,
                value,
            });
            self.scattered
                .get_or_insert_with(HashMap::new)
                .insert(slot, at);
            self.tail_at = at + 1;
            return Key::new(slot, generation, self.scope);
        }
        let slot = self.next_slot;
        self.next_slot = slot.checked_add(1).expect("arena exceeded u32::MAX slots");
        self.entries.push(Entry {
            slot,
            generation: 0,
            value,
        });
        Key::new(slot, 0, self.scope)
    }

    /// Where slot `index`'s value would sit in `entries`: past the end for
    /// a slot never handed out, `None` for a scattered slot with no value.
    #[inline]
    fn position(&self, index: u32) -> Option<usize> {
        if index >= self.tail {
            Some((index - self.tail) as usize + self.tail_at as usize)
        } else {
            self.scattered.as_ref()?.get(&index).map(|&at| at as usize)
        }
    }

    /// The position of `key`'s value, if the key is live here.
    #[inline]
    fn live(&self, key: Key<T>) -> Option<usize> {
        if key.scope != self.scope {
            return None;
        }
        let at = self.position(key.index)?;
        let entry = self.entries.get(at)?;
        (entry.slot == key.index && entry.generation == key.generation).then_some(at)
    }

    /// Borrow the value behind `key`, or `None` if the key is stale.
    #[must_use]
    #[inline]
    pub fn get(&self, key: Key<T>) -> Option<&T> {
        let at = self.live(key)?;
        self.entries.get(at).map(|entry| &entry.value)
    }

    /// Mutably borrow the value behind `key`, or `None` if the key is stale.
    #[inline]
    pub fn get_mut(&mut self, key: Key<T>) -> Option<&mut T> {
        let at = self.live(key)?;
        self.entries.get_mut(at).map(|entry| &mut entry.value)
    }

    /// Whether `key` resolves to a live entry.
    #[must_use]
    pub fn contains(&self, key: Key<T>) -> bool {
        self.live(key).is_some()
    }

    /// Move the tail's values into the scattered ones, leaving the tail
    /// empty and starting at the next fresh slot.
    #[allow(
        clippy::cast_possible_truncation,
        reason = "`insert` keeps positions within u32"
    )]
    fn scatter_tail(&mut self) {
        let from = self.tail_at as usize;
        if from < self.entries.len() {
            let map = self.scattered.get_or_insert_with(HashMap::new);
            for (at, entry) in self.entries.iter().enumerate().skip(from) {
                map.insert(entry.slot, at as u32);
            }
        }
        self.tail = self.next_slot;
        self.tail_at = self.entries.len() as u32;
    }

    /// Remove and return the value behind `key`, if it is live.
    ///
    /// The slot's generation is bumped, invalidating every outstanding copy of
    /// `key`, and a later insert reuses the slot.
    #[allow(
        clippy::cast_possible_truncation,
        reason = "`insert` keeps positions within u32"
    )]
    pub fn remove(&mut self, key: Key<T>) -> Option<T> {
        self.live(key)?;
        self.scatter_tail();
        let map = self.scattered.as_mut()?;
        let at = map.remove(&key.index)? as usize;
        let entry = self.entries.swap_remove(at);
        if let Some(moved) = self.entries.get(at) {
            map.insert(moved.slot, at as u32);
            self.ordered = false;
        }
        self.tail_at -= 1;
        // Saturating rather than wrapping: a slot recycled 4 billion times stops
        // being reusable, which is strictly better than handing out a generation
        // that collides with a key someone still holds.
        let next = entry.generation.saturating_add(1);
        if next != u32::MAX {
            self.free.push((key.index, next));
        }
        Some(entry.value)
    }

    /// Keep the entries `keep` answers `true` for and drop the rest.
    ///
    /// Every kept key still resolves to its value. A dropped key fails to
    /// resolve, and its slot is never handed out again, so no later entry
    /// answers to it. What is dropped stops costing memory: the arena holds
    /// the kept values, packed.
    #[allow(
        clippy::cast_possible_truncation,
        reason = "`insert` keeps positions within u32"
    )]
    pub fn retain(&mut self, mut keep: impl FnMut(Key<T>, &T) -> bool) {
        let scope = self.scope;
        self.entries
            .retain(|entry| keep(Key::new(entry.slot, entry.generation, scope), &entry.value));
        if !self.ordered {
            self.entries.sort_unstable_by_key(|entry| entry.slot);
            self.ordered = true;
        }
        // The longest run of entries ending at the last slot handed out,
        // every slot in it occupied, stays the tail; the rest are scattered.
        let mut tail_len = 0_usize;
        while tail_len < self.entries.len() {
            let entry = &self.entries[self.entries.len() - 1 - tail_len];
            if u64::from(entry.slot) + 1 + tail_len as u64 != u64::from(self.next_slot) {
                break;
            }
            tail_len += 1;
        }
        let tail_at = self.entries.len() - tail_len;
        self.tail = self.next_slot - tail_len as u32;
        self.tail_at = tail_at as u32;
        self.scattered = (tail_at > 0).then(|| {
            self.entries[..tail_at]
                .iter()
                .enumerate()
                .map(|(at, entry)| (entry.slot, at as u32))
                .collect()
        });
        self.entries.shrink_to_fit();
    }

    /// The positions of `entries`, in slot order.
    fn slot_order(&self) -> impl Iterator<Item = usize> + '_ {
        let sorted = (!self.ordered).then(|| {
            let mut order: Vec<usize> = (0..self.entries.len()).collect();
            order.sort_unstable_by_key(|&at| self.entries[at].slot);
            order
        });
        (0..self.entries.len()).map(move |k| sorted.as_ref().map_or(k, |order| order[k]))
    }

    /// Put `entries` in slot order, where removals have disturbed it.
    #[allow(
        clippy::cast_possible_truncation,
        reason = "`insert` keeps positions within u32"
    )]
    fn sort(&mut self) {
        if self.ordered {
            return;
        }
        self.scatter_tail();
        self.entries.sort_unstable_by_key(|entry| entry.slot);
        if let Some(map) = &mut self.scattered {
            for (at, entry) in self.entries.iter().enumerate() {
                map.insert(entry.slot, at as u32);
            }
        }
        self.ordered = true;
    }

    /// Iterate over live `(key, &value)` pairs, in slot order.
    pub fn iter(&self) -> impl Iterator<Item = (Key<T>, &T)> {
        let scope = self.scope;
        self.slot_order().map(move |at| {
            let entry = &self.entries[at];
            (Key::new(entry.slot, entry.generation, scope), &entry.value)
        })
    }

    /// Iterate over live `(key, &mut value)` pairs, in slot order.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = (Key<T>, &mut T)> {
        self.sort();
        let scope = self.scope;
        self.entries.iter_mut().map(move |entry| {
            (
                Key::new(entry.slot, entry.generation, scope),
                &mut entry.value,
            )
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
    pub fn into_values(mut self) -> impl Iterator<Item = T> {
        self.sort();
        self.entries.into_iter().map(|entry| entry.value)
    }

    /// Whether the arena has only ever been appended to: every slot occupied,
    /// every generation zero.
    ///
    /// When this holds, [`Arena::len`] is also the next index [`Arena::insert`]
    /// will hand out.
    #[must_use]
    pub fn is_dense(&self) -> bool {
        self.entries.len() == self.next_slot as usize
            && self.entries.iter().all(|entry| entry.generation == 0)
    }

    /// Remove every entry. Existing keys go stale, and their slots are not
    /// handed out again.
    pub fn clear(&mut self) {
        self.retain(|_, _| false);
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
    fn retained_keys_resolve_and_dropped_slots_are_never_handed_out_again() {
        let mut a = Arena::new();
        let keys: Vec<_> = (0..10_u32).map(|i| a.insert(i)).collect();
        a.retain(|_, v| v % 3 == 0);
        assert_eq!(a.len(), 4);
        for (i, key) in keys.iter().enumerate() {
            let want = (i % 3 == 0).then_some(u32::try_from(i).unwrap());
            assert_eq!(a.get(*key).copied(), want, "slot {i}");
        }
        assert_eq!(a.next_index(), Some(10), "a fresh insert appends");
        let fresh = a.insert(10);
        assert_eq!(fresh.index(), 10);
        assert!(keys.iter().all(|&k| k != fresh));
        assert_eq!(
            a.values().copied().collect::<Vec<_>>(),
            vec![0, 3, 6, 9, 10]
        );
        assert_eq!(a.key_at(3), Some(keys[3]));
        assert_eq!(a.key_at(4), None);

        // A clone answers to the same keys, and a removal from the retained
        // arena still refills its slot under a new generation.
        let copy = a.clone();
        assert!(keys.iter().all(|&k| copy.get(k) == a.get(k)));
        assert_eq!(a.remove(keys[6]), Some(6));
        assert_eq!(a.next_index(), None, "the freed slot is refilled first");
        let refill = a.insert(60);
        assert_eq!(refill.index(), 6);
        assert_eq!(a.get(keys[6]), None);
        assert_eq!(
            a.iter().map(|(k, v)| (k.index(), *v)).collect::<Vec<_>>(),
            vec![(0, 0), (3, 3), (6, 60), (9, 9), (10, 10)]
        );
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
        // answers about the wrong value, confidently, with nothing about the
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
        // A clone is the same document (a snapshot), so handles into it keep
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
