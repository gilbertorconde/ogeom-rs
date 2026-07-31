//! Shape construction and query — everything that does not require surface/surface
//! intersection.
//!
//! *Elsewhere:* the `BRepBuilderAPI`, `BRepPrimAPI`, `BRepAdaptor`, `BRepTools`,
//! `BRepLib`, `BRepGProp`, `BRepCheck`, `BRepClass3d`, `BRepExtrema`, `GeomAPI`,
//! `GCPnts` and `GProp` families, plus the classical curve constructors.
//!
//! Every operation in this crate emits history (`generated` / `modified` /
//! `is_deleted`) from the start. That is not optional and cannot be retrofitted:
//! downstream stable naming is built directly on it.
