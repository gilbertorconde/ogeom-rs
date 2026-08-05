//! The document's attribute layer: properties, materials, layers, and
//! validation values.
//!
//! Everything here is *data about* shapes rather than geometry: the
//! free-form key–value pairs an application pins to a face, the material a
//! body is meant to be cut from, the layers a drawing organizes itself by,
//! and the mass-property check values an exchange partner records so the
//! receiver can verify a translation did not quietly lose a boss. The
//! document stores and round-trips them; computing anything — a volume to
//! compare against a validation record — stays with the code that owns the
//! geometry.

use ogeom_math::Point;

use crate::structure::Colour;

/// A user-defined property's value.
#[derive(Debug, Clone, PartialEq)]
pub enum PropertyValue {
    /// Free text.
    Text(String),
    /// A number, meaning whatever the name says it means.
    Number(f64),
    /// A yes or a no.
    Flag(bool),
}

/// One user-defined property: a name and a value.
#[derive(Debug, Clone, PartialEq)]
pub struct Property {
    /// The key.
    pub name: String,
    /// The value.
    pub value: PropertyValue,
}

/// A material in the document's own list.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct MaterialId(pub(crate) usize);

impl MaterialId {
    /// The position in the document's material list.
    #[must_use]
    pub fn index(&self) -> usize {
        self.0
    }
}

/// A material: what a body is meant to be made of.
#[derive(Debug, Clone, PartialEq)]
pub struct Material {
    /// The name — "AISI 304", "PA12".
    pub name: String,
    /// Density in kilograms per cubic metre, where known.
    pub density: Option<f64>,
    /// A display colour, where one belongs to the material.
    pub colour: Option<Colour>,
}

/// A layer in the document's own list.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LayerId(pub(crate) usize);

impl LayerId {
    /// The position in the document's layer list.
    #[must_use]
    pub fn index(&self) -> usize {
        self.0
    }
}

/// A layer: a named grouping with a visibility flag.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Layer {
    /// The name.
    pub name: String,
    /// Whether the layer is shown.
    pub visible: bool,
}

/// Mass-property check values, recorded so a receiver can verify a
/// translation preserved the body they describe.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ValidationProperties {
    /// Enclosed volume, in the model's cubic length unit.
    pub volume: f64,
    /// Surface area, in the squared unit.
    pub area: f64,
    /// The centroid.
    pub centroid: Point,
}

impl ValidationProperties {
    /// Whether another set of values agrees within a relative tolerance —
    /// the centroid compared against the body's own size, taken from the
    /// cube root of its volume.
    #[must_use]
    pub fn agrees_with(&self, other: &Self, relative: f64) -> bool {
        let close = |a: f64, b: f64, scale: f64| (a - b).abs() <= scale.abs().max(1.0) * relative;
        let size = self.volume.abs().cbrt();
        close(self.volume, other.volume, self.volume)
            && close(self.area, other.area, self.area)
            && close(self.centroid.x, other.centroid.x, size)
            && close(self.centroid.y, other.centroid.y, size)
            && close(self.centroid.z, other.centroid.z, size)
    }
}
