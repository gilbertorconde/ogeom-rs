//! Each face's box, and what readers keep of a whole shape, kept beside the
//! topology and forgotten when anything they were read from changes.
//!
//! The box is filled through a shared reference, the first time a reader
//! asks, so a model shared across threads fills it as freely as one held
//! alone. It is forgotten through the exclusive reference that changes the
//! face: a node handed out for editing forgets the box of every face above
//! it, and a surface rewritten in place forgets every box.

use std::any::Any;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

use ogeom_math::Aabb;

use crate::shape::Shape;

/// The kept face boxes, one slot per face; each face node names its slot.
#[derive(Default)]
pub(crate) struct FaceBoxes {
    /// Per face, its box in its node's own frame, once found.
    boxes: Vec<OnceLock<Aabb>>,
    /// The geometry store's revision the boxes were found against.
    pub(crate) revision: u64,
    /// Whether any box has been kept since every box was last forgotten.
    any: AtomicBool,
}

impl Clone for FaceBoxes {
    fn clone(&self) -> Self {
        Self {
            boxes: self.boxes.clone(),
            revision: self.revision,
            any: AtomicBool::new(self.any.load(Ordering::Relaxed)),
        }
    }
}

impl core::fmt::Debug for FaceBoxes {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let kept = self.boxes.iter().filter(|b| b.get().is_some()).count();
        f.debug_struct("FaceBoxes")
            .field("faces", &self.boxes.len())
            .field("kept", &kept)
            .finish()
    }
}

impl FaceBoxes {
    /// A slot for a face just added.
    pub(crate) fn add(&mut self) -> u32 {
        let slot = u32::try_from(self.boxes.len()).unwrap_or(u32::MAX);
        self.boxes.push(OnceLock::new());
        slot
    }

    /// The box kept in `slot`, or the place to keep one.
    pub(crate) fn slot(&self, slot: u32) -> Option<&OnceLock<Aabb>> {
        self.boxes.get(slot as usize)
    }

    /// Forget the box kept in `slot`.
    pub(crate) fn forget(&mut self, slot: u32) {
        if let Some(kept) = self.boxes.get_mut(slot as usize) {
            kept.take();
        }
    }

    /// Keep only the slots `kept` names, in that order: the `k`-th becomes
    /// slot `k`, with its box.
    pub(crate) fn keep_slots(&mut self, kept: &[u32]) {
        let boxes = kept
            .iter()
            .map(|&slot| {
                self.boxes
                    .get_mut(slot as usize)
                    .map(core::mem::take)
                    .unwrap_or_default()
            })
            .collect();
        self.boxes = boxes;
    }

    /// Keep `found` in `slot` unless a box is kept there already, and
    /// return the one kept.
    pub(crate) fn keep(&self, slot: &OnceLock<Aabb>, found: Aabb) -> Aabb {
        self.any.store(true, Ordering::Relaxed);
        *slot.get_or_init(|| found)
    }

    /// Forget every box.
    pub(crate) fn forget_all(&mut self) {
        if !*self.any.get_mut() {
            return;
        }
        for kept in &mut self.boxes {
            kept.take();
        }
        *self.any.get_mut() = false;
    }
}

/// How many shapes keep what was read of them: the few large solids an
/// operation runs on in turn, the oldest forgotten first.
const SHAPES: usize = 4;

/// What readers keep of a shape and everything below it, forgotten when
/// any node below it is changed or the geometry is rewritten in place.
///
/// Filled through a shared reference, as the face boxes are. A node is
/// changed only through the exclusive reference that forgets every shape
/// holding it, so what is kept was read from the nodes as they stand.
#[derive(Default)]
pub(crate) struct KeptReads {
    shapes: Mutex<Vec<KeptShape>>,
}

/// One shape and what was read of it.
#[derive(Clone)]
struct KeptShape {
    shape: Shape,
    /// The node index of the shape and of every node below it.
    below: Arc<NodeSet>,
    /// The geometry store's count of rewrites the reads were made against.
    edits: u64,
    reads: Vec<Arc<dyn Any + Send + Sync>>,
}

impl Clone for KeptReads {
    fn clone(&self) -> Self {
        Self {
            shapes: Mutex::new(self.lock().clone()),
        }
    }
}

impl core::fmt::Debug for KeptReads {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let shapes = self.lock();
        f.debug_struct("KeptReads")
            .field("shapes", &shapes.len())
            .field(
                "reads",
                &shapes.iter().map(|s| s.reads.len()).sum::<usize>(),
            )
            .finish()
    }
}

impl KeptReads {
    fn lock(&self) -> std::sync::MutexGuard<'_, Vec<KeptShape>> {
        self.shapes
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    fn shapes_mut(&mut self) -> &mut Vec<KeptShape> {
        self.shapes
            .get_mut()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    /// The read of type `T` kept for `shape` that `same` accepts, made
    /// against the geometry as `edits` counts it.
    pub(crate) fn find<T: Any + Send + Sync>(
        &self,
        shape: &Shape,
        edits: u64,
        same: impl Fn(&T) -> bool,
    ) -> Option<Arc<T>> {
        let shapes = self.lock();
        let kept = shapes
            .iter()
            .find(|k| k.edits == edits && k.shape == *shape)?;
        kept.reads.iter().find_map(|read| {
            let read = Arc::clone(read).downcast::<T>().ok()?;
            same(&read).then_some(read)
        })
    }

    /// The nodes kept as below `shape`, where anything is kept for it.
    pub(crate) fn below(&self, shape: &Shape, edits: u64) -> Option<Arc<NodeSet>> {
        self.lock()
            .iter()
            .find(|k| k.edits == edits && k.shape == *shape)
            .map(|k| Arc::clone(&k.below))
    }

    /// Keep `read` for `shape`, whose nodes are `below`, beside what is
    /// kept for it already; a read of the same type that `same` accepts
    /// is replaced.
    pub(crate) fn keep<T: Any + Send + Sync>(
        &self,
        shape: &Shape,
        edits: u64,
        below: Arc<NodeSet>,
        read: Arc<T>,
        same: impl Fn(&T) -> bool,
    ) {
        let mut shapes = self.lock();
        shapes.retain(|k| k.edits == edits);
        let at = if let Some(at) = shapes.iter().position(|k| k.shape == *shape) {
            at
        } else {
            if shapes.len() >= SHAPES {
                shapes.remove(0);
            }
            shapes.push(KeptShape {
                shape: shape.clone(),
                below,
                edits,
                reads: Vec::new(),
            });
            shapes.len() - 1
        };
        let reads = &mut shapes[at].reads;
        reads.retain(|r| !Arc::clone(r).downcast::<T>().is_ok_and(|r| same(&r)));
        reads.push(read);
    }

    /// Forget what is kept for every shape holding node `index`.
    pub(crate) fn forget(&mut self, index: u32) {
        let shapes = self.shapes_mut();
        if !shapes.is_empty() {
            shapes.retain(|k| !k.below.contains(index));
        }
    }

    /// Forget what is kept for every shape `keep` refuses.
    pub(crate) fn retain(&mut self, keep: impl Fn(&Shape) -> bool) {
        self.shapes_mut().retain(|k| keep(&k.shape));
    }
}

/// A set of node indices, a bit each from the first word holding one.
#[derive(Debug, Default)]
pub(crate) struct NodeSet {
    /// The index of `words[0]`'s first bit, over 64.
    first: usize,
    words: Vec<u64>,
}

impl NodeSet {
    /// Add `index`; whether it was not in the set.
    pub(crate) fn insert(&mut self, index: u32) -> bool {
        let at = index as usize / 64;
        let bit = 1_u64 << (index % 64);
        if self.words.is_empty() {
            self.first = at;
        } else if at < self.first {
            let grow = self.first - at;
            self.words.splice(0..0, core::iter::repeat_n(0, grow));
            self.first = at;
        }
        let word = at - self.first;
        if word >= self.words.len() {
            self.words.resize(word + 1, 0);
        }
        let was = self.words[word] & bit != 0;
        self.words[word] |= bit;
        !was
    }

    /// Whether `index` is in the set.
    pub(crate) fn contains(&self, index: u32) -> bool {
        let at = index as usize / 64;
        at.checked_sub(self.first)
            .and_then(|word| self.words.get(word))
            .is_some_and(|w| w & (1_u64 << (index % 64)) != 0)
    }
}

#[cfg(test)]
mod tests {
    use super::NodeSet;

    #[test]
    fn a_node_set_holds_what_was_put_in_wherever_it_starts() {
        let mut set = NodeSet::default();
        let put = [700_u32, 70, 7, 64, 63, 701, 6400];
        for i in put {
            assert!(set.insert(i));
        }
        assert!(!set.insert(70));
        for i in 0..7000 {
            assert_eq!(set.contains(i), put.contains(&i), "{i}");
        }
    }
}
