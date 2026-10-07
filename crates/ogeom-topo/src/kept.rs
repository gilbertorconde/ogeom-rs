//! Each face's box, kept beside the topology and forgotten when anything
//! the face is made of changes.
//!
//! The box is filled through a shared reference, the first time a reader
//! asks, so a model shared across threads fills it as freely as one held
//! alone. It is forgotten through the exclusive reference that changes the
//! face: a node handed out for editing forgets the box of every face above
//! it, and a surface rewritten in place forgets every box.

use std::sync::OnceLock;
use std::sync::atomic::{AtomicBool, Ordering};

use ogeom_math::Aabb;

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
