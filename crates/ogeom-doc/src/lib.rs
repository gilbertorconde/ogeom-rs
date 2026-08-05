//! Document and product structure.
//!
//! Parts, assemblies and instances with placements; colours and materials per shape
//! and per sub-shape; layers; names and user properties; PMI and GD&T — dimensions,
//! tolerances, datums, annotations; validation properties; transactional undo and
//! redo; and the native persistence format.
//!
//! Instancing leans directly on the location chain (`docs/DATA_MODEL.md` §2): an
//! assembly of ten thousand identical fasteners is ten thousand placements over one
//! piece of geometry, and that only works because placement is structural rather
//! than a matrix to be compared with an epsilon.
//!
//! Deliberately not a label-and-attribute tree. The capability is required; that
//! particular design, and compatibility with anyone's file format, are not.

pub mod pmi;
pub mod structure;

pub use pmi::{Datum, Dimension, GeometricTolerance, MeasureKind, Pmi};
pub use structure::{Colour, Document, Instance, Occurrence, Product, ProductId, ProductKind};
