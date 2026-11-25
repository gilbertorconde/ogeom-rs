//! Overflow handling for chamfers: a bevel reaching geometry other than its
//! own two faces (ADR 0001, amendment "overflow handling").

use crate::{Chamfer, ChamferCorner, Overflow};
use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances, ogeom_bail, ogeom_err};
use ogeom_topo::{Model, Shape, ShapeType};

/// [`chamfer_edges_cornered`](crate::chamfer_edges_cornered), a bevel
/// reaching geometry other than its own two faces handled as `overflow`
/// says.
///
/// A flat bevel takes two of the behaviours of [`Overflow`], tried in this
/// order:
///
/// - `blend_into_round`: where a bevel's leg reaches an existing round of
///   the body, the bevel runs on over that round with no step. Built for a
///   boss's foot: the closed rim between a planar floor and a cylindrical
///   wall standing square on it, or an open arc of it whose ends stand on
///   planes parallel to the wall's axis (leaning from the plane through it
///   by at most 60 degrees), every edge of the rim asking the same
///   legs, whose leg on the floor crosses straight edges of the floor
///   along which cylindrical rounds, convex or concave, meet the floor
///   tangentially, over arcs apart. Over the floor the bevel is the cone
///   between its legs; past each tangent line its leg on the floor is
///   where the wall's offset meets the round, and each section is the
///   straight ruling between the legs, the pieces meeting G1. An open
///   arc's ends are capped flush in their planes through the axis, or
///   trimmed by their planes beside it.
/// - `trim_through`: the bevel keeps its own surface, and the faces it
///   reaches are trimmed or extended to meet it: the engine's own
///   construction, what the chamfer built before this option set.
///
/// `hold_edge` and `hold_edge_at_end` do not apply to a bevel, which has
/// no ball to rest on an edge: they are ignored, as are the narrow-passage
/// switches.
///
/// With [`Overflow::default`] this is
/// [`chamfer_edges_cornered`](crate::chamfer_edges_cornered). A behaviour
/// switched off passes the request to the next. With `trim_through` off, a
/// result where a bevel changes a face other than its own edge's and those
/// at its ends (whichever other edge of the request is seated on it), or
/// carries a face on past its side over nothing, is refused by name, as is
/// a request that blend into round would take while it is off.
///
/// A foot bevel whose leg on the floor touches the round's tangent line,
/// or crosses it at under a hundredth of a radian, leaves the floor and the
/// round wedges too sharp to build at the crossings: where neither
/// behaviour builds it, it is refused by that name.
///
/// # Errors
///
/// As [`chamfer_edges_cornered`](crate::chamfer_edges_cornered), and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction), by
/// name, where a bevel reaches other geometry and no behaviour that is on
/// builds it.
#[deprecated(
    since = "0.22.0",
    note = "use `BlendRequest`: one `.edge(e, spec)` per spec, then `.chamfer_corner(c).overflow(o).build(..)`"
)]
pub fn chamfer_edges_overflow(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    corner: ChamferCorner,
    overflow: Overflow,
    tol: Tolerances,
) -> OgeomResult<Built> {
    crate::BlendRequest::bevels(specs)
        .chamfer_corner(corner)
        .overflow(overflow)
        .build(model, solid, tol)
}

/// [`chamfer_edges_overflow`]'s construction, which a
/// [`BlendRequest`](crate::BlendRequest) of bevels alone builds under a
/// point corner or overflow options.
pub(crate) fn bevels_overflow(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    corner: ChamferCorner,
    overflow: Overflow,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if specs.is_empty() {
        ogeom_bail!(Construction, "a chain of no edges bevels nothing");
    }
    if overflow.blend_into_round {
        if let Some(built) = bevelled_into_round(model, solid, specs, tol)? {
            return Ok(built);
        }
    } else if !overflow.trim_through && reaches_round(model, solid, specs, tol)? {
        ogeom_bail!(
            Construction,
            "the chamfer reaches an existing round of the body, and both blend into \
             round and trim through are off; see docs/CAPABILITIES.md, \
             fillet.chamfer-overflow"
        );
    }
    let built = match crate::chamfer_corner::cornered_as_today(model, solid, specs, corner, tol) {
        Ok(built) => built,
        Err(e @ OgeomError::Construction(_)) => {
            return Err(grazing(model, solid, specs, tol)?.unwrap_or(e));
        }
        Err(e) => return Err(e),
    };
    if overflow.trim_through {
        return Ok(built);
    }
    let edges: Vec<(Shape, f64)> = specs
        .iter()
        .map(|(e, spec)| (e.clone(), longer_leg(spec)))
        .collect();
    if let Some(through) = crate::engine::trims_through(model, solid, &edges, &built, tol)? {
        return Err(crate::fillet::refused_through(
            through,
            "chamfer",
            "fillet.chamfer-overflow",
        ));
    }
    Ok(built)
}

/// The longer of the two legs `spec` asks: the reach of its bevel from
/// the edge.
pub(crate) fn longer_leg(spec: &Chamfer) -> f64 {
    match spec {
        Chamfer::Symmetric(d) => *d,
        Chamfer::Distances {
            on_face, on_other, ..
        } => on_face.max(*on_other),
        Chamfer::Angle {
            distance, angle, ..
        } => distance.max(distance * angle.tan().abs()),
    }
}

/// The foot bevel of [`Overflow::blend_into_round`], on a solid or on each
/// solid of a compound; `None` where it does not take the request.
///
/// # Errors
///
/// As reading the solid and the booleans, other than their refusals.
pub(crate) fn bevelled_into_round(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    if let Some(built) = crate::compound::per_solid(
        model,
        solid,
        specs,
        |(edge, _)| edge,
        |model, part, own| on_solid(model, part, own, tol),
    )? {
        return Ok(Some(built));
    }
    if model.kind_of(solid)? != ShapeType::Solid {
        return Ok(None);
    }
    on_solid(model, solid, specs, tol)
}

fn on_solid(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some((edges, size)) = crate::engine::overflow_bevel_size(model, solid, specs, tol)? else {
        return Ok(None);
    };
    crate::engine::overflow_blended(model, solid, &edges, size, tol)
}

/// Whether blend into round would take `specs` on `solid`, or on the
/// solid of a compound holding them.
fn reaches_round(
    model: &Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<bool> {
    for (part, edges, size) in feet(model, solid, specs, tol)? {
        if crate::engine::overflow_reach(model, &part, &edges, size, tol)? {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The refusal of a foot bevel of `specs` whose leg on the floor grazes
/// the tangent line of an existing round, on `solid` or a solid of a
/// compound holding it; `None` where none does.
///
/// # Errors
///
/// As reading the solid.
pub(crate) fn grazing(
    model: &Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<OgeomError>> {
    for (part, edges, size) in feet(model, solid, specs, tol)? {
        let Some(angle) = crate::engine::overflow_grazes(model, &part, &edges, size, tol)? else {
            continue;
        };
        return Ok(Some(ogeom_err!(
            Construction,
            "the chamfer's leg on the floor grazes the tangent line of an existing \
             round, meeting it at {angle:.1e} rad, under a hundredth: the floor \
             and the round are left wedges too sharp to build where the leg meets \
             it; see docs/CAPABILITIES.md, fillet.chamfer-overflow"
        )));
    }
    Ok(None)
}

/// Each solid of `solid` (itself, or each of a compound's) holding edges
/// of `specs` that ask one bevel of a foot rim: the solid, the edges and
/// the bevel.
fn feet(
    model: &Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Vec<(Shape, Vec<Shape>, crate::engine::OverflowSize)>> {
    let parts = if model.kind_of(solid)? == ShapeType::Compound {
        model.children_of(solid)?
    } else {
        vec![solid.clone()]
    };
    let mut out = Vec::new();
    for part in parts {
        let edges = ogeom_topo::explore_unique(model, &part, ShapeType::Edge)?;
        let own: Vec<(Shape, Chamfer)> = specs
            .iter()
            .filter(|(e, _)| edges.iter().any(|o| o.is_same(e)))
            .cloned()
            .collect();
        if own.is_empty() {
            continue;
        }
        if let Some((edges, size)) = crate::engine::overflow_bevel_size(model, &part, &own, tol)? {
            out.push((part, edges, size));
        }
    }
    Ok(out)
}
