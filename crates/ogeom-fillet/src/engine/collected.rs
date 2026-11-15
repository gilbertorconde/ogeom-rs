//! Rounds whose bands the current path's section builders make, on seats
//! the planar engine does not read (ADR 0001, amendments "varying rounds on
//! curved seats" and "round distances on curved seats"): a varying radius
//! over a rim, an arc, a drum's ruling, a tangent loop round a pad, a
//! marched crease; a round set back its own distance on each face over a
//! rim, a spline crease, a slanted edge. Each band is built on the solid as
//! given, with the request's other edges as its mates, and collected as a
//! wedge rather than applied; the engine applies the wedges in one boolean
//! each, as it applies a constant round's on a marched seat.

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
pub(crate) fn varying(
    model: &mut Model,
    solid: &Shape,
    chain: &[Shape],
    law: &RadiusLaw,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let largest = crate::law::Profile::new(law, false, tol)?.largest();
    rounds(model, solid, chain, largest, tol, |model| {
        crate::fillet::chain_law_bands(model, solid, chain, law, tol)
    })
}

/// `edge` rounded on `solid` set back `on_face` on `face` and `on_other` on
/// its other face, with `profile`, by its band's wedges; `None` as for
/// [`varying`].
///
/// # Errors
///
/// As [`varying`].
#[allow(clippy::too_many_arguments, reason = "fillet_edge_distances' own")]
pub(crate) fn distances(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    face: &Shape,
    on_face: f64,
    on_other: f64,
    profile: crate::RoundProfile,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let edges = std::slice::from_ref(edge);
    rounds(model, solid, edges, on_face.max(on_other), tol, |model| {
        crate::asym::distances_bands(model, solid, edge, face, on_face, on_other, profile, tol)
    })
}

/// `asymmetric`, a request of asymmetric rounds alone (`blend_edges`), each
/// tangent chain of them closed as one, rounded on `solid` by the bands'
/// wedges; `None` as for [`varying`].
///
/// # Errors
///
/// As [`varying`].
pub(crate) fn asymmetric(
    model: &mut Model,
    solid: &Shape,
    asymmetric: &[crate::asym::Asymmetric],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let edges: Vec<Shape> = asymmetric.iter().map(|a| a.edge.clone()).collect();
    let size = asymmetric
        .iter()
        .map(|a| a.on_face.max(a.on_other))
        .fold(0.0, f64::max);
    rounds(model, solid, &edges, size, tol, |model| {
        crate::asym::after_symmetric(
            model,
            &[],
            asymmetric.to_vec(),
            |_| Ok(Built::new(solid.clone(), ogeom_algo::History::new())),
            tol,
        )
    })
}

/// `edges` rounded on `solid` by the wedges `build` sets aside on a copy
/// of `model`, each named for its edge (the first of `edges` where it
/// names none), applied as the one-boolean fallback's tools. An edge no
/// wedge names must be one the build deleted (an arc of a crease another
/// arc's band covers): it goes with that band.
fn rounds(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    size: f64,
    tol: Tolerances,
    build: impl FnOnce(&mut Model) -> OgeomResult<Built>,
) -> OgeomResult<Option<Built>> {
    let Some(first) = edges.first() else {
        return Ok(None);
    };
    let mut scratch = model.clone();
    let (built, wedges) = crate::support::collecting_wedges(|| build(&mut scratch));
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
    let named = |e: &Shape| {
        wedges.iter().any(|w| {
            w.edges
                .iter()
                .chain(w.credit.iter().map(|(_, c)| c))
                .any(|n| n.is_same(e))
        })
    };
    let (carried, gone): (Vec<&Shape>, Vec<&Shape>) = edges.iter().partition(|e| named(e));
    if gone.iter().any(|e| !collected.history.is_deleted(e)) {
        return Ok(None);
    }
    let once: Vec<(Shape, f64)> = carried.iter().map(|e| ((*e).clone(), size)).collect();
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
