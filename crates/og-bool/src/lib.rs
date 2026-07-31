//! Boolean operations.
//!
//! *Elsewhere:* `BOPDS`, `BOPAlgo`, `BOPTools`, `IntTools` and `BRepAlgoAPI`.
//!
//! The architectural insight worth preserving: **fuse, common, cut, section and
//! split are one algorithm — general fuse — plus a selection predicate.** The
//! pipeline is
//!
//! 1. a data structure of indexed sub-shapes, per-type interference lists
//!    (V/V, V/E, V/F, E/E, E/F, F/F), pave blocks and common blocks;
//! 2. a pave filler running in strictly increasing dimension, ending in face/face
//!    intersection, section-edge construction and pcurve generation;
//! 3. a builder that splits faces in 2D parametric space, unifies same-domain
//!    faces, rebuilds solids from face sets, and repairs tolerances;
//! 4. filters that select from the general-fuse result.
