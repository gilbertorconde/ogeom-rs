//! IGES exchange.
//!
//! *Elsewhere:* the IGES half of the exchange layer. The format is ANSI
//! Y14.26M / IGES 5.3, a fixed-column deck from the 1980s that a great deal
//! of legacy CAD still speaks; implementing it from the published
//! specification is interoperation, the same footing STEP was built on.

mod parse;
mod read;
mod write;

pub use read::{IgesImport, IgesReport, read_iges};
pub use write::write_iges;

/// The entity types the reader and writer translate, by number and by the
/// name the specification gives each. The names are load-bearing twice over:
/// refusal messages quote them, and `tools/parity.py exchange` greps them to
/// keep the parity ledger's IGES coverage figure regenerable from this file
/// rather than asserted.
pub const SUPPORTED_ENTITIES: &[(i64, &str)] = &[
    (100, "CIRCULAR_ARC"),
    (102, "COMPOSITE_CURVE"),
    (104, "CONIC_ARC"),
    (108, "PLANE"),
    (110, "LINE"),
    (112, "SPLINE_CURVE"),
    (116, "POINT"),
    (120, "SURFACE_OF_REVOLUTION"),
    (122, "TABULATED_CYLINDER"),
    (123, "DIRECTION"),
    (124, "TRANSFORMATION_MATRIX"),
    (126, "B_SPLINE_CURVE"),
    (128, "B_SPLINE_SURFACE"),
    (141, "BOUNDARY"),
    (142, "CURVE_ON_SURFACE"),
    (143, "BOUNDED_SURFACE"),
    (144, "TRIMMED_SURFACE"),
    (186, "MANIFOLD_SOLID"),
    (190, "PLANE_SURFACE"),
    (192, "RIGHT_CIRCULAR_CYLINDRICAL_SURFACE"),
    (194, "RIGHT_CIRCULAR_CONICAL_SURFACE"),
    (196, "SPHERICAL_SURFACE"),
    (198, "TOROIDAL_SURFACE"),
    (314, "COLOR"),
    (502, "VERTEX_LIST"),
    (504, "EDGE_LIST"),
    (508, "LOOP"),
    (510, "FACE"),
    (514, "SHELL"),
];

/// The specification's name for an entity type, where this module knows it.
#[must_use]
pub fn entity_name(kind: i64) -> Option<&'static str> {
    SUPPORTED_ENTITIES
        .iter()
        .find(|(k, _)| *k == kind)
        .map(|(_, n)| *n)
}
