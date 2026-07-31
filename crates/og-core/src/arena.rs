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

use core::fmt;
use core::hash::{Hash, Hasher};
use core::marker::PhantomData;

/// A handle into an [`Arena<T>`].
///
/// Phantom-typed, so a `Key<Face>` cannot be used to index an `Arena<Edge>`.
/// The marker is `fn() -> T` so the key stays `Copy`, `Send` and `Sync`
/// regardless of `T`.
pub struct Key<T> {
    index: u32,
    generation: u32,
    marker: PhantomData<fn() -> T>,
}

impl<T> Key<T> {
    fn new(index: u32, generation: u32) -> Self {
        Self {
            index,
            generation,
            marker: PhantomData,
        }
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
        Self {
            index,
            generation,
            marker: PhantomData,
        }
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
        self.index == other.index && self.generation == other.generation
    }
}
impl<T> Eq for Key<T> {}
impl<T> Hash for Key<T> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.index.hash(state);
        self.generation.hash(state);
    }
}
impl<T> PartialOrd for Key<T> {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl<T> Ord for Key<T> {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        (self.index, self.generation).cmp(&(other.index, other.generation))
    }
}
impl<T> fmt::Debug for Key<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Key({}v{})", self.index, self.generation)
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
        }
    }

    /// An empty arena with room for `capacity` entries.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            slots: Vec::with_capacity(capacity),
            free_head: None,
            len: 0,
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

    /// Insert a value, returning its key.
    ///
    /// # Panics
    ///
    /// If the arena exceeds `u32::MAX` slots. A single model reaching four
    /// billion topological entities is a bug elsewhere, not a case to handle.
    #[allow(clippy::expect_used, reason = "documented panic; see # Panics")]
    pub fn insert(&mut self, value: T) -> Key<T> {
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
                Key::new(index, generation)
            }
            None => {
                let index = u32::try_from(self.slots.len()).expect("arena exceeded u32::MAX slots");
                self.slots.push(Slot::Occupied {
                    generation: 0,
                    value,
                });
                Key::new(index, 0)
            }
        }
    }

    /// Borrow the value behind `key`, or `None` if the key is stale.
    #[must_use]
    pub fn get(&self, key: Key<T>) -> Option<&T> {
        match self.slots.get(key.index as usize)? {
            Slot::Occupied { generation, value } if *generation == key.generation => Some(value),
            _ => None,
        }
    }

    /// Mutably borrow the value behind `key`, or `None` if the key is stale.
    pub fn get_mut(&mut self, key: Key<T>) -> Option<&mut T> {
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
        self.slots
            .iter()
            .enumerate()
            .filter_map(|(i, slot)| match slot {
                Slot::Occupied { generation, value } => {
                    // `insert` refuses to grow past u32::MAX, so this cannot truncate.
                    #[allow(clippy::cast_possible_truncation)]
                    Some((Key::new(i as u32, *generation), value))
                }
                Slot::Vacant { .. } => None,
            })
    }

    /// Iterate over live `(key, &mut value)` pairs, in slot order.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = (Key<T>, &mut T)> {
        self.slots
            .iter_mut()
            .enumerate()
            .filter_map(|(i, slot)| match slot {
                Slot::Occupied { generation, value } =>
                {
                    #[allow(clippy::cast_possible_truncation)]
                    Some((Key::new(i as u32, *generation), value))
                }
                Slot::Vacant { .. } => None,
            })
    }

    /// Iterate over live values.
    pub fn values(&self) -> impl Iterator<Item = &T> {
        self.iter().map(|(_, v)| v)
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
