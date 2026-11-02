//! Construction geometry: named points, axes, curves and surfaces a product
//! carries for reference.
//!
//! Construction geometry has no topology and encloses nothing: a datum point,
//! a reference axis, a fitted section curve, a bounded reference plane. It is
//! document data like notes and views: it takes part in undo, persists in the
//! native document format, and goes through STEP as the geometric set it is
//! there.

use ogeom_geom::{Curve, SurfaceGeometry};
use ogeom_math::{Frame, Point};

/// The geometry of one construction item.
#[derive(Debug, Clone, PartialEq)]
pub enum ConstructionGeometry {
    /// A point.
    Point(Point),
    /// An axis placement: an origin, its `z` axis and its `x` axis, with `y`
    /// as `z × x`.
    Axis(Frame),
    /// A curve over a finite stretch of its parameter.
    Curve {
        /// The curve.
        curve: Curve,
        /// The stretch of its parameter the item is.
        range: (f64, f64),
    },
    /// A surface over a finite window of its parameters.
    Surface {
        /// The surface.
        surface: SurfaceGeometry,
        /// The window's `u` stretch.
        u: (f64, f64),
        /// The window's `v` stretch.
        v: (f64, f64),
    },
}

/// A named construction item and the product it belongs to.
#[derive(Debug, Clone, PartialEq)]
pub struct Construction {
    /// The item's name, as the author gave it.
    pub name: String,
    /// What it is.
    pub geometry: ConstructionGeometry,
    /// The product carrying it.
    pub product: crate::structure::ProductId,
}
