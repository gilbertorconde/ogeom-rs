//! Fitting curves and surfaces to points.
//!
//! # Curves through points
//!
//! [`interpolate`], [`approximate`] and [`approximate_within`] fit a
//! B-spline to points. They live in [`ogeom_geom::fit`] and are re-exported
//! here.
//!
//! # Primitives
//!
//! The other half of fitting measured points: the plane, the sphere, the
//! cylinder, the cone, the torus, the circle or the line that comes nearest
//! them, with every point's residual ([`plane`], [`sphere`], [`cylinder`],
//! [`cone`], [`torus`], [`circle3`], [`line3`], [`circle2`], [`line2`]).
//!
//! # Zones
//!
//! Form as a tolerance limits it is not the least-squares fit's peak to
//! peak but the narrowest zone holding every point: [`zone`] fits the
//! minimum zone of a plane, a line, a circle, a cylinder or a sphere, and
//! the largest inscribed and smallest circumscribed circles and
//! cylinders.
//!
//! # Regions
//!
//! On a mesh, which points belong to the surface is itself the question:
//! [`grow`](fn@grow) takes clicks on a surface, grows the region of the mesh a kind
//! of surface covers around them, and fits it. For a circle or a line it
//! traces the crease nearest the click instead, along edges where the
//! surface turns sharply: a bore's rim, a box's edge. For a free-form
//! surface, [`grow_smooth`] grows the region the surface stays smooth over
//! and fits nothing, leaving its points for a free-form fit.

mod grow;
mod primitive;
pub(crate) mod segments;
mod tails;
pub mod zone;

pub use grow::{
    AnyFit, Crease, FitKind, GrowOptions, Grown, GrownRegion, Meeting, SeedReport, grow,
    grow_smooth,
};
pub use ogeom_geom::fit::{Spacing, approximate, approximate_within, interpolate};
pub use primitive::{
    Fit, FitFixed, FitOptions, Side, circle2, circle3, cone, cylinder, line2, line3, plane, sphere,
    torus,
};
