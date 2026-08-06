//! PMI: dimensions, geometric tolerances, datums, their targets, and the
//! drawing of them.
//!
//! The machine-readable annotations AP242 calls *semantic* PMI — the values,
//! not the leader lines. A dimension carries what it measures and its
//! plus/minus bounds; a geometric tolerance carries its kind, magnitude and
//! the datums it references; a datum is a letter on a feature, and its
//! *targets* are the pads a fixture actually contacts it at. Everything
//! anchors to topology nodes, the same way colours and names do, so an
//! annotation survives every placement of the shape it describes.
//!
//! And the *presentation* kind, which is the leader lines: where each
//! annotation is drawn, in which plane, as which polylines. It is kept
//! separate because it is separate — a drawing places its callouts where a
//! draughtsman put them, and neither half derives from the other. What ties
//! them is [`Callout::annotates`], and a document may carry either alone.

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
    /// Zone and material-condition modifiers, as lower-case words:
    /// `maximum_material_requirement`, `unequally_disposed`, … — the
    /// exchange vocabulary for Ⓜ, Ⓤ and their kin.
    pub modifiers: Vec<String>,
    /// The datum letters the tolerance references, in precedence order. A
    /// composite reference — two datums acting as one, ISO's `A-B` — is a
    /// single entry with its labels hyphen-joined.
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

/// A document's PMI, in file order: the semantic annotations, the targets a
/// datum is actually contacted at, and the drawn presentation.
#[derive(Debug, Clone, Default)]
pub struct Pmi {
    /// Dimensional characteristics.
    pub dimensions: Vec<Dimension>,
    /// Geometric tolerances.
    pub tolerances: Vec<GeometricTolerance>,
    /// Datums.
    pub datums: Vec<Datum>,
    /// The targets datums are established at.
    pub targets: Vec<DatumTarget>,
    /// The drawn annotations.
    pub callouts: Vec<Callout>,
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
        self.dimensions.is_empty()
            && self.tolerances.is_empty()
            && self.datums.is_empty()
            && self.targets.is_empty()
            && self.callouts.is_empty()
    }

    /// The targets belonging to one datum letter, in the order the file gave.
    pub fn targets_of<'a>(&'a self, datum: &'a str) -> impl Iterator<Item = &'a DatumTarget> + 'a {
        self.targets.iter().filter(move |t| t.datum == datum)
    }
}

/// Where a datum is actually contacted: the target a fixture touches it at.
///
/// A datum plane on a casting is not contacted over its whole face — it rests
/// on three pads — and the drawing says so with targets: `A1`, `A2`, `A3`,
/// each a point, a line or an area of stated size. The distinction matters to
/// anything that inspects the part, because the datum it should establish is
/// the one the targets define and not the nominal surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DatumTargetKind {
    /// A point contact.
    Point,
    /// A line contact, of the stated length.
    Line {
        /// How long the contact runs.
        length: f64,
    },
    /// A rectangular area.
    Rectangle {
        /// Along the target's own reference direction.
        length: f64,
        /// Across it.
        width: f64,
    },
    /// A circular area.
    Circle {
        /// The contact patch's diameter.
        diameter: f64,
    },
}

/// One datum target: which datum, which target of it, and where.
#[derive(Debug, Clone)]
pub struct DatumTarget {
    /// The datum's letter — `A` for `A1`.
    pub datum: String,
    /// The target's own number — `1` for `A1`.
    pub index: u32,
    /// Point, line, or area, with its size.
    pub kind: DatumTargetKind,
    /// Where the target sits, in the part's own coordinates.
    pub at: ogeom_math::Point,
    /// The target's own frame: its normal, and the direction its length runs
    /// in. Absent where the file placed the target by point alone, which a
    /// point target does not need.
    pub frame: Option<ogeom_math::Frame>,
    /// The topology the target is placed on.
    pub items: Vec<TShapeId>,
}

impl DatumTarget {
    /// The identifier a drawing shows: the letter and the number, `A1`.
    #[must_use]
    pub fn identifier(&self) -> String {
        format!("{}{}", self.datum, self.index)
    }
}

/// One drawn annotation: the geometry a viewer puts on the screen.
///
/// *Presentation* PMI, as against the semantic kind above. The two are
/// separate on purpose and in the file: the semantic annotation says a
/// tolerance is `0.1` and controls this face, and the presentation says where
/// its frame and leader are drawn. Neither derives from the other — a drawing
/// places its callouts where a draughtsman put them — so a document that
/// wants both carries both, and [`Callout::annotates`] is the link between.
#[derive(Debug, Clone)]
pub struct Callout {
    /// The name the file gave it, which is how a drawing names its own
    /// annotations: `Flatness.1`, `Linear Size.3`.
    pub name: String,
    /// The plane the annotation is drawn in, where the file stated one.
    pub plane: Option<ogeom_math::Frame>,
    /// The drawn geometry: polylines, in the part's own coordinates.
    pub polylines: Vec<Vec<ogeom_math::Point>>,
    /// Which semantic annotation this draws, where the file said.
    pub annotates: Option<Annotated>,
}

/// Which semantic annotation a callout draws.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Annotated {
    /// The dimension at this index of [`Pmi::dimensions`].
    Dimension(usize),
    /// The tolerance at this index of [`Pmi::tolerances`].
    Tolerance(usize),
    /// The datum at this index of [`Pmi::datums`].
    Datum(usize),
}
