//! Pcurves for imported faces, shared by the exchange readers.
//!
//! The machinery itself lives in `ogeom_algo::pcurve_fit` — fitting a trim
//! by projection is geometry, not exchange — and the readers reach it
//! through this shim under the names they always used.

pub(crate) use ogeom_algo::pcurve_fit::fit_projected_pcurve;
