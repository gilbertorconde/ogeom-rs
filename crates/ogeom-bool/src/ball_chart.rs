//! The booleans run on a ball, or a sphere face, restated on a chart whose
//! poles stand clear of the other operand (see `ogeom_algo`'s sphere
//! charts), trying the next best chart where the best refuses.

use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_topo::{Model, Shape};

/// How many charts of a ball whose poles stand clear of the other operand
/// a boolean is tried on, the best first, before it is refused.
const CHARTS_TRIED: usize = 6;

/// `run` on the operands, where either is or holds a whole ball or a sphere
/// face with a pole on or beside the other's surfaces, with that ball or
/// face first charted about an axis whose poles stand clear; the history
/// then runs from the operands as given. Otherwise `run` on the operands
/// as given.
///
/// Where the run on the best such chart refuses, it is run on the next
/// best, up to [`CHARTS_TRIED`] charts: each is the same ball, restated
/// exactly, so any that resolves the boolean resolves it rightly. Where
/// every one refuses, the boolean is refused by name, the last refusal
/// with it. The chart the operands come with is not taken: there a
/// section beside a pole can close where the chart's longitude has
/// collapsed, and give a cut that is valid but short, or a face wound
/// inward.
pub(crate) fn with_poles_clear(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
    mut run: impl FnMut(&mut Model, &Shape, &Shape) -> OgeomResult<Built>,
) -> OgeomResult<Built> {
    #[cfg(feature = "test-hooks")]
    if crate::work::is_unturned() {
        return run(model, a, b);
    }
    let mut last: Option<OgeomError> = None;
    for rank in 0..CHARTS_TRIED {
        let Some((ta, tb, before)) = ogeom_algo::operands_clear_of_poles(model, a, b, rank, tol)?
        else {
            break;
        };
        #[cfg(feature = "test-hooks")]
        let refused = crate::work::is_turn_refused(rank);
        #[cfg(not(feature = "test-hooks"))]
        let refused = false;
        let result = if refused {
            Err(ogeom_core::ogeom_err!(
                NotDone,
                "the turned run is taken as refused"
            ))
        } else {
            run(model, &ta, &tb)
        };
        match result {
            Ok(built) => return Ok(Built::new(built.shape, before.then(&built.history))),
            Err(e @ (OgeomError::Cancelled | OgeomError::Dangling(_))) => return Err(e),
            Err(e) => last = Some(e),
        }
    }
    match last {
        Some(e) => Err(ogeom_core::ogeom_err!(
            NotDone,
            "a ball's pole stands on or beside the other operand, and the boolean on every \
             chart of the ball whose poles stand clear of it is refused (the last: {e}); on \
             its own chart a section beside a pole is not resolved soundly; see \
             docs/CAPABILITIES.md, bool.booleans"
        )),
        None => run(model, a, b),
    }
}
