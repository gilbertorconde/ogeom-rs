//! Document and product structure.
//!
//! *Elsewhere:* `XCAFDoc`, `XCAFPrs` and `XCAFDimTolObjects` — the exchange
//! document, which is part of DataExchange and in scope. The application
//! framework those sit on top of is not, which is what the last paragraph here
//! is about.
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

pub mod attributes;
pub mod pmi;
pub mod structure;
pub mod view;

pub use attributes::{
    Layer, LayerId, Material, MaterialId, Property, PropertyValue, Texture, TextureId,
    TextureMapping, ValidationProperties,
};
pub use pmi::{
    Annotated, Callout, Datum, DatumTarget, DatumTargetKind, Dimension, GeometricTolerance,
    MeasureKind, Pmi,
};
pub use structure::{Colour, Document, Instance, Occurrence, Product, ProductId, ProductKind};
pub use view::{Note, View};
