//! Degenerate edges restored where a face's boundary crosses a pole, by
//! [`ogeom_algo::pole_crossings`] on every face of a shape.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use crate::Reshape;

/// Stage, for every wire whose crossings of poles are read wrongly the
/// shortest way across, a rebuilt wire with a degenerate edge spanning
/// each; returns how many edges were made. The wires are solved by
/// [`ogeom_algo::pole_crossings`]; a face whose surface is not in `model`
/// is left as it is.
pub(crate) fn stage_pole_edges(
    model: &mut Model,
    shape: &Shape,
    reshape: &mut Reshape,
    tol: Tolerances,
) -> OgeomResult<usize> {
    let mut made = 0;
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
            continue;
        };
        if model.geometry().surface(data.surface).is_none() {
            continue;
        }
        for (wire, rebuilt) in ogeom_algo::pole_crossings(model, &face, tol)? {
            let before = model.ordered_children_of(&wire)?.len();
            made += model.ordered_children_of(&rebuilt)?.len() - before;
            reshape.replace(&wire, rebuilt);
        }
    }
    Ok(made)
}
