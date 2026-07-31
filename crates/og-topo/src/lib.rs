//! The boundary-representation data model. What this model can express is what
//! the kernel can do, so its invariants are normative — see `docs/DATA_MODEL.md`.
//!
//! *Elsewhere:* `TopoDS`, `TopAbs`, `TopExp`, `TopTools`, `TopLoc` and `BRep`.
//!
//! Non-negotiable invariants:
//!
//! - A shape is a **(TShape, Location, Orientation) triple** — a positionless,
//!   orientationless shared topology node plus a placement plus an orientation.
//!   Cheap to copy; the heavy data lives once in an arena.
//! - **Location is a chain** of `(datum, power)` pairs, not a flat 4x4. Composition
//!   is list concatenation and the composed transform is computed lazily. This is
//!   what lets an assembly of 10,000 identical bolts share one piece of geometry,
//!   and what keeps identical-instance detection a cheap chain comparison.
//! - **Orientation composes multiplicatively** down the tree. An edge's effective
//!   orientation inside a face depends on that face's orientation inside its shell.
//! - **Identity is a trichotomy**: same (tshape + location), equal (+ orientation),
//!   partner (tshape only) — each with a consistent hasher. Conflating these is a
//!   classic source of silent bugs.
//! - **Tolerances are per-entity**, with `tol(vertex) >= tol(edge) >= tol(face)`, and
//!   they grow through operations. This is the kernel's answer to inexact
//!   arithmetic; every production kernel works this way.
//! - **An edge carries a list of representations**, not one curve: a 3D curve, one
//!   pcurve per adjacent face, two pcurves for a seam edge on a closed surface, plus
//!   polygon and triangulation representations — with a `same_parameter` flag and a
//!   repair routine.

pub mod entity;
pub mod location;
pub mod model;
pub mod shape;
pub mod tessellation;

pub use entity::{
    CurveId, EdgeData, EdgeRepr, FaceData, GeometryStore, NodeData, PCurveId, SurfaceId,
    TriangulationId, VertexData, check_containment, enforce_containment,
};
pub use location::{Datum, DatumId, DatumStore, Location};
pub use model::{Filter, Model, ModelParts, ancestors_of, explore, explore_unique};
pub use shape::{Orientation, PartnerKey, SameKey, Shape, ShapeType, TShape, TShapeId};
pub use tessellation::Triangulation;
