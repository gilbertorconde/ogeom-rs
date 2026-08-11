//! Saved views and standalone notes.
//!
//! A saved view is how an annotated part organises its own presentation: a
//! named camera, and the subset of the drawn callouts meant to be visible
//! from it. Real AP242 files use views to sort PMI by sheet, and reading
//! them without views flattens a structure the author meant. A note is
//! smaller still — text a person attached to the document or to a product,
//! with no geometry behind it.
//!
//! Both are document data like colours and layers: they take part in undo,
//! they persist in the native document format, and views round-trip through
//! STEP as the draughting models they are there.

use ogeom_math::{Frame, Plane};

/// A saved view: a camera and what it is meant to show.
#[derive(Debug, Clone, PartialEq)]
pub struct View {
    /// The view's name, as the author gave it.
    pub name: String,
    /// The camera: origin is the eye's target, `z` the viewing direction,
    /// `x` the screen's right.
    pub frame: Frame,
    /// A section plane the view clips at, where one was stated.
    pub clipping: Option<Plane>,
    /// Indices into the document's PMI callouts: the annotations this view
    /// presents. An index rather than a copy, so restyling a callout
    /// restyles it in every view that shows it.
    pub callouts: Vec<usize>,
}

/// A note: text attached to the document or to a product, with an author.
#[derive(Debug, Clone, PartialEq)]
pub struct Note {
    /// Who wrote it.
    pub author: String,
    /// The text.
    pub text: String,
    /// The product it is about, or `None` for the document itself.
    pub product: Option<crate::structure::ProductId>,
}
