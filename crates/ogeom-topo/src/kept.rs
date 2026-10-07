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
use smallvec::SmallVec;

use crate::shape::{ShapeType, TShape, TShapeId};

/// No face slot: the node is not a face.
const NONE: u32 = u32::MAX;

/// The kept face boxes, and the links that say which faces a node bounds.
#[derive(Default)]
pub(crate) struct FaceBoxes {
    /// Per node index, the node's slot in `boxes`, or [`NONE`].
    slot_of: Vec<u32>,
    /// Per face, its box in its node's own frame, once found.
    boxes: Vec<OnceLock<Aabb>>,
    /// Per node index, the edges, wires and faces that hold it: the way up
    /// from a vertex, an edge or a wire to the faces it bounds.
    held_by: Vec<SmallVec<[u32; 2]>>,
    /// The geometry store's revision the boxes were found against.
    pub(crate) revision: u64,
    /// Whether any box has been kept since every box was last forgotten.
    any: AtomicBool,
}

impl Clone for FaceBoxes {
    fn clone(&self) -> Self {
        Self {
            slot_of: self.slot_of.clone(),
            boxes: self.boxes.clone(),
            held_by: self.held_by.clone(),
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
    /// Record a node just added: a slot if it is a face, and its links to
    /// the children it holds if it is part of a face.
    pub(crate) fn note(&mut self, id: TShapeId, node: &TShape) {
        let at = id.index() as usize;
        if self.slot_of.len() <= at {
            self.slot_of.resize(at + 1, NONE);
            self.held_by.resize(at + 1, SmallVec::new());
        }
        if node.kind() == ShapeType::Face {
            self.slot_of[at] = u32::try_from(self.boxes.len()).unwrap_or(NONE);
            self.boxes.push(OnceLock::new());
        }
        if matches!(
            node.kind(),
            ShapeType::Edge | ShapeType::Wire | ShapeType::Face
        ) {
            for child in node.children() {
                let below = child.node().index() as usize;
                if self.held_by.len() <= below {
                    self.slot_of.resize(below + 1, NONE);
                    self.held_by.resize(below + 1, SmallVec::new());
                }
                if !self.held_by[below].contains(&id.index()) {
                    self.held_by[below].push(id.index());
                }
            }
        }
    }

    /// The slot of the face at node index `node`.
    pub(crate) fn slot(&self, node: u32) -> Option<&OnceLock<Aabb>> {
        let slot = *self.slot_of.get(node as usize)?;
        self.boxes.get(slot as usize)
    }

    /// Forget the box of every face that holds node index `node`, the node
    /// itself included where it is a face.
    pub(crate) fn forget_above(&mut self, node: u32) {
        let mut stack: SmallVec<[u32; 8]> = smallvec::smallvec![node];
        while let Some(at) = stack.pop() {
            let Some(&slot) = self.slot_of.get(at as usize) else {
                continue;
            };
            if let Some(kept) = self.boxes.get_mut(slot as usize) {
                kept.take();
            }
            if let Some(up) = self.held_by.get(at as usize) {
                stack.extend(up.iter().copied());
            }
        }
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
