//! Rounds on seats only the march speaks (ADR 0001, amendment "marched seats"): a
//! crease between faces no closed-form section reads (a spline edge, the
//! seam of two cylinders crossing, a converted solid's every edge). Each
//! edge's strip is the march's band, closed into its wedge by its legs on
//! the two faces and its caps, built on the solid as given with the
//! request's other edges as its mates; the engine applies the wedges in
//! one boolean each, as it applies every tool of the one-boolean fallback.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::fallback::{Sense, Tool};

/// `once` rounded on its marched seats: every edge's wedge (see
/// [`crate::fillet::marched_wedges`]) cut where the edge is convex and
/// fused where it is concave, each wedge's faces named for the edge it
/// rounds (or the edge its builder credits them to). A one-edge request
/// whose band the boolean refuses takes the strip beside it off whole
/// where the march reads that round. `None` where an edge is not a marched
/// seat or its wedge does not build, a boolean refuses, or the result is
/// not valid: the request is the current path's then.
///
/// # Errors
///
/// As the march and the booleans, other than their refusals, which
/// decline; a one-edge request's refusal the current path would raise,
/// by name.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut scratch = model.clone();
    let Some(wedges) = crate::fillet::marched_wedges(&mut scratch, solid, once, tol)? else {
        return Ok(None);
    };
    if let Some(built) = applied(&mut scratch, solid, once, &wedges, tol)? {
        *model = scratch;
        return Ok(Some(built));
    }
    let [(edge, radius)] = once else {
        return Ok(None);
    };
    let mut scratch = model.clone();
    let Some(wedges) =
        crate::fillet::marched_strip_taken_whole(&mut scratch, solid, edge, *radius, tol)?
    else {
        return Ok(None);
    };
    let Some(built) = applied(&mut scratch, solid, once, &wedges, tol)? else {
        return Ok(None);
    };
    *model = scratch;
    Ok(Some(built))
}

/// `wedges` applied to `solid` as the one-boolean fallback's tools, each
/// wedge's faces named for its edge or the edge it credits them to;
/// `None` where a boolean refuses, an edge of `once` is left nothing of,
/// or the result is not valid.
fn applied(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    wedges: &[crate::support::Wedge],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut tools = Vec::with_capacity(wedges.len());
    for wedge in wedges {
        let Some(own) = wedge.edges.first() else {
            return Ok(None);
        };
        let faces = explore_unique(model, &wedge.solid, ShapeType::Face)?
            .into_iter()
            .map(|face| {
                let from = wedge
                    .credit
                    .iter()
                    .find(|(f, _)| f.is_same(&face))
                    .map_or_else(|| own.clone(), |(_, e)| e.clone());
                (face, from)
            })
            .collect();
        tools.push(Tool {
            solid: wedge.solid.clone(),
            sense: if wedge.additive {
                Sense::Fuse
            } else {
                Sense::Cut
            },
            faces,
        });
    }
    let edges: Vec<Shape> = once.iter().map(|(e, _)| e.clone()).collect();
    let Some(built) = super::fallback::applied(model, solid, &edges, &tools, tol)? else {
        return Ok(None);
    };
    if !ogeom_algo::check(model, &built.shape, tol)?.is_valid() {
        return Ok(None);
    }
    Ok(Some(built))
}
