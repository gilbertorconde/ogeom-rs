//! Rounds and bevels in one request, closed as one corner where they meet.

use crate::{Chamfer, ChamferCorner, RoundProfile};
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
    // Rounds alone, of one radius and of their own distances: the blend
    // engine first, as one request.
    if !asymmetric.is_empty()
        && bevels.is_empty()
        && let Some(request) =
            crate::engine::Request::from_mixed_rounds(model, solid, &rounds, &asymmetric, tol)?
        && let Some(built) = crate::engine::run(model, solid, &request, tol)?
    {
        return Ok(built);
    }
    if !asymmetric.is_empty() {
        // The corners where an asymmetric round meets the request's other
        // edges, read off the solid as it stands, their edges set apart.
        let corners =
            crate::asym::asymmetric_corners(model, solid, &rounds, &bevels, &asymmetric, tol)?;
        let at_corner = |e: &Shape| corners.edges.iter().any(|c| c.is_same(e));
        let (corner_rounds, rounds): (Vec<_>, Vec<_>) =
            rounds.into_iter().partition(|(e, _)| at_corner(e));
        let (corner_bevels, bevels): (Vec<_>, Vec<_>) =
            bevels.into_iter().partition(|(e, _)| at_corner(e));
        let (corner_asymmetric, asymmetric): (Vec<_>, Vec<_>) =
            asymmetric.into_iter().partition(|a| at_corner(&a.edge));
        let others: Vec<Shape> = rounds
            .iter()
            .map(|(e, _)| e.clone())
            .chain(bevels.iter().map(|(e, _)| e.clone()))
            .chain(corners.edges.iter().cloned())
            .collect();
        let symmetric_only = rounds.is_empty() && bevels.is_empty();
        let built = crate::asym::after_symmetric(
            model,
            &others,
            asymmetric,
            |model| {
                if symmetric_only {
                    Ok(Built::new(solid.clone(), ogeom_algo::History::new()))
                } else {
                    blend_symmetric(model, solid, rounds, bevels, corner, tol)
                }
            },
            tol,
        )?;
        return crate::asym::close_asymmetric_corners(
            model,
            solid,
            built,
            &corners,
            (&corner_rounds, &corner_bevels, &corner_asymmetric),
            tol,
        );
    }
    blend_symmetric(model, solid, rounds, bevels, corner, tol)
}

/// [`blend_edges`] of the rounds and bevels of one radius each.
fn blend_symmetric(
    model: &mut Model,
    solid: &Shape,
    mut rounds: Vec<(Shape, f64)>,
    mut bevels: Vec<(Shape, Chamfer)>,
    corner: ChamferCorner,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if rounds.is_empty() {
        return crate::chamfer_edges_cornered(model, solid, &bevels, corner, tol);
    }
    if bevels.is_empty() {
        return crate::fillet_edges_with(model, solid, &rounds, tol);
    }
    rounds.sort_by_key(|(e, _)| e.node());
    bevels.sort_by_key(|(e, _)| e.node());
    // The bevels, then the rounds on what they leave (ADR 0001, section 13,
    // question 3): two features, whatever the request's order.
    let beveled = crate::chamfer_edges_cornered(model, solid, &bevels, corner, tol)?;
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
