//! Varying rounds on curved seats (ADR 0001, amendment "varying rounds on
//! curved seats"): a law of radius over a chain the planar engine does not
//! read (a rim, an open arc, a drum's ruling, a tangent loop round a pad, a
//! marched crease). Each band is the march under the law, closed into its
//! wedge, built on the solid as given with the chain's other edges as its
//! mates; the engine applies the wedges in one boolean each, as it applies
//! a constant round's on a marched seat.

use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_topo::{Model, Shape};

use crate::law::RadiusLaw;

/// `chain` rounded under `law` on `solid` by the bands' wedges; `None`
/// where a band does not build as a wedge set aside, a boolean refuses, an
/// edge is left nothing of, or the result is not valid: the request is the
/// current path's then.
///
/// # Errors
///
/// [`OgeomError::Cancelled`] where the work is cancelled, and reading the
/// chain.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    chain: &[Shape],
    law: &RadiusLaw,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some(first) = chain.first() else {
        return Ok(None);
    };
    let mut scratch = model.clone();
    let (built, wedges) = crate::support::collecting_wedges(|| {
        crate::fillet::chain_law_bands(&mut scratch, solid, chain, law, tol)
    });
    let collected = match built {
        Ok(built) if built.shape.is_same(solid) && !wedges.is_empty() => built,
        Err(e @ OgeomError::Cancelled) => return Err(e),
        _ => return Ok(None),
    };
    let wedges: Vec<crate::support::Wedge> = wedges
        .into_iter()
        .map(|mut w| {
            if w.edges.is_empty() {
                w.edges.push(first.clone());
            }
            w
        })
        .collect();
    // An open chain of arcs of one curve is one band over them all, named
    // for the first: the others go with it, as the band's own history
    // reads.
    let named = |e: &Shape| {
        wedges.iter().any(|w| {
            w.edges
                .iter()
                .chain(w.credit.iter().map(|(_, c)| c))
                .any(|n| n.is_same(e))
        })
    };
    let (carried, gone): (Vec<&Shape>, Vec<&Shape>) = chain.iter().partition(|e| named(e));
    if gone.iter().any(|e| !collected.history.is_deleted(e)) {
        return Ok(None);
    }
    let largest = crate::law::Profile::new(law, false, tol)?.largest();
    let once: Vec<(Shape, f64)> = carried.iter().map(|e| ((*e).clone(), largest)).collect();
    let Some(mut built) = super::marched_seats::applied(&mut scratch, solid, &once, &wedges, tol)?
    else {
        return Ok(None);
    };
    for e in gone {
        built.history.delete(e);
    }
    *model = scratch;
    Ok(Some(built))
}
