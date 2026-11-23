//! Rounds and bevels in one request, closed as one corner where they meet.

use crate::{Chamfer, ChamferCorner, Overflow, RoundProfile};
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_topo::{Model, Shape};

/// How [`blend_edges`] blends one edge.
#[derive(Debug, Clone)]
pub enum Blend {
    /// A round of this radius, as [`fillet_edges_with`](crate::fillet_edges_with).
    Round(f64),
    /// A bevel, as [`chamfer_edges_with`](crate::chamfer_edges_with).
    Chamfer(Chamfer),
    /// A round meeting `face` `on_face` from the edge and the edge's other
    /// face `on_other` from it, its section `profile`, as
    /// [`fillet_edge_distances`](crate::fillet_edge_distances).
    RoundDistances {
        /// The face the first distance runs along.
        face: Shape,
        /// The distance along `face`.
        on_face: f64,
        /// The distance along the edge's other face.
        on_other: f64,
        /// The section between the two contacts.
        profile: RoundProfile,
    },
}

/// Round and bevel the edges of `solid` as `edges` says, in one request.
///
/// A request of rounds alone is
/// [`fillet_edges_with`](crate::fillet_edges_with), and one of bevels
/// alone [`chamfer_edges_cornered`](crate::chamfer_edges_cornered) with
/// `corner`, exactly. Where both are asked, the bevels are cut first and
/// the rounds then, each kind in the order of its edges' nodes, so the
/// order of the request does not matter; at a vertex where only one kind
/// meets, that is what the kind's own call builds.
///
/// Where a round meets a bevel at a vertex, the rounds meet what the
/// bevels left there, as the rounds' own call on the bevelled solid does:
/// two features, as mainstream applications take them (ADR 0001, section
/// 13, question 3).
///
/// Where a blend reaches geometry other than its own two faces, the
/// behaviours of [`Overflow`] are tried in its order, every one on:
/// [`blend_edges_overflow`] chooses among them.
///
/// # Errors
///
/// As the rounds' and the bevels' own calls, and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// an edge is asked twice, or a round's edge does not come through the
/// bevels as one edge.
pub fn blend_edges(
    model: &mut Model,
    solid: &Shape,
    edges: &[(Shape, Blend)],
    corner: ChamferCorner,
    tol: Tolerances,
) -> OgeomResult<Built> {
    blend_edges_overflow(model, solid, edges, corner, Overflow::default(), tol)
}

/// [`blend_edges`], a blend reaching geometry other than its own two faces
/// handled as `overflow` says.
///
/// With [`Overflow::default`] this is [`blend_edges`]. The behaviours that
/// apply:
///
/// - `blend_into_round`, for a request of rounds alone of one radius
///   (circular rounds, and two equal distances of the elliptical profile,
///   which are the circular round): the boss's foot of
///   [`fillet_edges_overflow`](crate::fillet_edges_overflow), bent into
///   the round it reaches. For bevels alone, the foot bevel of
///   [`chamfer_edges_overflow`](crate::chamfer_edges_overflow); in a
///   request mixing bevels and rounds, the bevels take it as there.
/// - `trim_through`: the construction [`blend_edges`] built before this
///   option set. With it off, a result changing any face other than the
///   edges' own and those at their ends is refused by name.
///
/// Hold edge and the narrow passages are not built for this call, and
/// neither is blend into round for two different distances or a conic
/// profile: switching them on changes nothing and switching them off
/// refuses nothing. A request blend into round would take while it and
/// trim through are both off is refused by name.
///
/// # Errors
///
/// As [`blend_edges`], and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction), by
/// name, where a blend reaches other geometry and no behaviour that is on
/// builds it.
pub fn blend_edges_overflow(
    model: &mut Model,
    solid: &Shape,
    edges: &[(Shape, Blend)],
    corner: ChamferCorner,
    overflow: Overflow,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut rounds: Vec<(Shape, f64)> = Vec::new();
    let mut bevels: Vec<(Shape, Chamfer)> = Vec::new();
    let mut asymmetric: Vec<crate::asym::Asymmetric> = Vec::new();
    for (k, (edge, blend)) in edges.iter().enumerate() {
        if edges[..k].iter().any(|(e, _)| e.is_same(edge)) {
            ogeom_bail!(Construction, "an edge is asked to be blended twice");
        }
        match blend {
            Blend::Round(r) => rounds.push((edge.clone(), *r)),
            Blend::Chamfer(c) => bevels.push((edge.clone(), c.clone())),
            Blend::RoundDistances {
                face,
                on_face,
                on_other,
                profile,
            } => {
                // Equal distances, elliptical: the circular round, which
                // meets the request's other rounds and bevels as one.
                if *profile == RoundProfile::Elliptical && on_face == on_other {
                    let r = crate::asym::equal_round(model, solid, edge, face, *on_face, tol)?;
                    rounds.push((edge.clone(), r));
                } else {
                    asymmetric.push(crate::asym::Asymmetric {
                        edge: edge.clone(),
                        face: face.clone(),
                        on_face: *on_face,
                        on_other: *on_other,
                        profile: *profile,
                    });
                }
            }
        }
    }
    if bevels.is_empty()
        && asymmetric.is_empty()
        && let Some(&(_, radius)) = rounds.first()
        && rounds.iter().all(|&(_, r)| r == radius)
    {
        let shapes: Vec<Shape> = rounds.iter().map(|(e, _)| e.clone()).collect();
        if let Some(built) =
            crate::fillet::blend_into_round_step(model, solid, &shapes, radius, overflow, tol)?
        {
            return Ok(built);
        }
    }
    if rounds.is_empty() && asymmetric.is_empty() {
        return crate::chamfer_edges_overflow(model, solid, &bevels, corner, overflow, tol);
    }
    let sized: Vec<(Shape, f64)> = edges
        .iter()
        .map(|(e, blend)| {
            let size = match blend {
                Blend::Round(r) => *r,
                Blend::Chamfer(c) => crate::chamfer_overflow::longer_leg(c),
                Blend::RoundDistances {
                    on_face, on_other, ..
                } => on_face.max(*on_other),
            };
            (e.clone(), size)
        })
        .collect();
    let built = blended_as_today(
        model,
        solid,
        rounds,
        bevels,
        &asymmetric,
        corner,
        overflow,
        tol,
    );
    crate::fillet::trimmed_through(model, solid, &sized, built, overflow, tol)
}

/// [`blend_edges`]' construction of the parsed request, its bevels taking
/// `overflow` as [`chamfer_edges_overflow`](crate::chamfer_edges_overflow)
/// does.
#[allow(
    clippy::too_many_arguments,
    reason = "the request's three kinds, the corner and the overflow"
)]
fn blended_as_today(
    model: &mut Model,
    solid: &Shape,
    rounds: Vec<(Shape, f64)>,
    bevels: Vec<(Shape, Chamfer)>,
    asymmetric: &[crate::asym::Asymmetric],
    corner: ChamferCorner,
    overflow: Overflow,
    tol: Tolerances,
) -> OgeomResult<Built> {
    // Rounds alone, of one radius and of their own distances: the blend
    // engine first, as one request.
    if !asymmetric.is_empty()
        && bevels.is_empty()
        && let Some(request) =
            crate::engine::Request::from_mixed_rounds(model, solid, &rounds, asymmetric, tol)?
        && let Some(built) = crate::engine::run(model, solid, &request, tol)?
    {
        return Ok(built);
    }
    if !asymmetric.is_empty() {
        return crate::engine::mixed_asymmetric(
            model, solid, &rounds, &bevels, asymmetric, corner, tol,
        );
    }
    blend_symmetric(model, solid, rounds, bevels, corner, overflow, tol)
}

/// [`blend_edges`] of the rounds and bevels of one radius each.
fn blend_symmetric(
    model: &mut Model,
    solid: &Shape,
    mut rounds: Vec<(Shape, f64)>,
    mut bevels: Vec<(Shape, Chamfer)>,
    corner: ChamferCorner,
    overflow: Overflow,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if rounds.is_empty() {
        return crate::chamfer_edges_overflow(model, solid, &bevels, corner, overflow, tol);
    }
    if bevels.is_empty() {
        return crate::fillet_edges_with(model, solid, &rounds, tol);
    }
    rounds.sort_by_key(|(e, _)| e.node());
    bevels.sort_by_key(|(e, _)| e.node());
    // The bevels, then the rounds on what they leave (ADR 0001, section 13,
    // question 3): two features, whatever the request's order.
    let beveled = crate::chamfer_edges_overflow(model, solid, &bevels, corner, overflow, tol)?;
    let mut traced: Vec<(Shape, f64)> = Vec::with_capacity(rounds.len());
    for (edge, r) in &rounds {
        let now = beveled.history.trace(edge);
        let [one] = now else {
            ogeom_bail!(
                Construction,
                "a round's edge does not come through the bevels as one edge"
            );
        };
        traced.push((one.clone(), *r));
    }
    let rounded = crate::fillet_edges_with(model, &beveled.shape, &traced, tol)?;
    Ok(Built {
        shape: rounded.shape.clone(),
        history: beveled.history.then(&rounded.history),
    })
}
