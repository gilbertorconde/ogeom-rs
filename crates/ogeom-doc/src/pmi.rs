//! Semantic PMI: dimensions, geometric tolerances, datums.
//!
//! The machine-readable annotations AP242 calls semantic PMI — the values,
//! not the leader lines. A dimension carries what it measures and its
//! plus/minus bounds; a geometric tolerance carries its kind, magnitude and
//! the datums it references; a datum is a letter on a feature. Everything
//! anchors to topology nodes, the same way colours and names do, so an
//! annotation survives every placement of the shape it describes.

use ogeom_topo::TShapeId;

/// What a dimension's value measures.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeasureKind {
    /// A length, in the document's length unit.
    Length,
    /// An angle, in radians.
    Angle,
}

/// A dimensional characteristic: a size or a location, with its values.
#[derive(Debug, Clone)]
pub struct Dimension {
    /// What the file calls it: `diameter`, `linear distance`, …
    pub name: String,
    /// The stated values — one for a plain dimension, several when the file
    /// states a value with explicit bounds.
    pub values: Vec<f64>,
    /// Length or angle.
    pub kind: MeasureKind,
    /// The upper allowance, when a plus/minus tolerance applies.
    pub plus: Option<f64>,
    /// The lower allowance — typically negative — when one applies.
    pub minus: Option<f64>,
    /// The topology the dimension measures, one group per feature: a size
    /// has one group, a location has one per end.
    pub features: Vec<Vec<TShapeId>>,
    /// Whether the dimension runs *between* two features — a location —
    /// rather than sizing one.
    pub location: bool,
}

impl Dimension {
    /// Every measured node, features flattened.
    pub fn items(&self) -> impl Iterator<Item = TShapeId> + '_ {
        self.features.iter().flatten().copied()
    }
}

/// A geometric tolerance: flatness, position, profile, and their kin.
#[derive(Debug, Clone)]
pub struct GeometricTolerance {
    /// The kind, as a lower-case word: `flatness`, `position`,
    /// `surface_profile`, `perpendicularity`, …
    pub kind: String,
    /// The annotation's own name.
    pub name: String,
    /// The tolerance zone's magnitude.
    pub magnitude: f64,
    /// The datum letters the tolerance references, in precedence order.
    pub datums: Vec<String>,
    /// The topology the tolerance controls.
    pub items: Vec<TShapeId>,
}

/// A datum: a letter naming a feature other annotations reference.
#[derive(Debug, Clone)]
pub struct Datum {
    /// The letter: `A`, `B`, …
    pub label: String,
    /// The feature's topology.
    pub items: Vec<TShapeId>,
}

/// A document's semantic PMI, in file order.
#[derive(Debug, Clone, Default)]
pub struct Pmi {
    /// Dimensional characteristics.
    pub dimensions: Vec<Dimension>,
    /// Geometric tolerances.
    pub tolerances: Vec<GeometricTolerance>,
    /// Datums.
    pub datums: Vec<Datum>,
}

impl Pmi {
    /// No annotations.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Whether anything is annotated.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.dimensions.is_empty() && self.tolerances.is_empty() && self.datums.is_empty()
    }
}
