//! The old round of a concave edge whose ball rolls on the far edge of the
//! narrow face beside it, where [`crate::engine::rolled`] reads a run.

use crate::engine::rolled::rolled_run;
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_topo::{Model, Shape};

/// `None` where [`rolled_run`] finds no run of the ball resting on the far
/// edge of the narrow face beside `edge`; where it finds one, the request
/// is the engine's and is refused by name here.
///
/// # Errors
///
/// As [`rolled_run`], and the refusal where it finds a run.
pub(crate) fn rolled_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    // A ball on a narrow face's far edge, a shelf, a rim, carried across a
    // change of support or down a staircase is the engine's or refused by
    // name.
    match rolled_run(model, solid, edge, radius, tol)? {
        None => Ok(None),
        Some(_) => Err(crate::fillet::declined()),
    }
}
