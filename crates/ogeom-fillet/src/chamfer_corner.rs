//! Where three bevels meet at a vertex: a point or a facet.

use crate::Chamfer;
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_topo::{Model, Shape};

/// How a chamfer request closes a vertex where three or more of its bevels
/// meet.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChamferCorner {
    /// The bevels run on to the point they share and meet each other from
    /// there, with no facet between them; at a convex vertex of four bevels
    /// or more on planes sharing no point, each runs on until the others
    /// cut it.
    Point,
    /// Each bevel stops at the plane through the points where two bevels'
    /// contact lines cross on the face they share, and the corner is the
    /// polygon there: a triangle at a corner of three bevels. The
    /// default: what [`chamfer_edges`](crate::chamfer_edges),
    /// [`chamfer_edges_with`](crate::chamfer_edges_with) and
    /// [`chamfer_edges_checked`](crate::chamfer_edges_checked) build. At a
    /// box's corner of bevels `d` it takes `d³/12` more than
    /// [`ChamferCorner::Point`].
    Facet,
}

/// [`chamfer_edges_with`](crate::chamfer_edges_with), the vertices where
/// three or more of its bevels meet closed as `corner` says.
///
/// [`ChamferCorner::Facet`] is exactly
/// [`chamfer_edges_with`](crate::chamfer_edges_with): the blend engine
/// closes each such vertex, convex, concave or mixed, with the polygon
/// through its faces' setback corners, and a vertex of four bevels or more
/// on planes whose corners leave one plane with the patch through them
/// (the bilinear patch at four, the regular polygon's Wachspress
/// coordinates at five or more). A box's corner of three bevels `d` gains
/// the triangle cutting each edge at `2d` from the corner.
///
/// [`ChamferCorner::Point`] is the engine's mitre: on planes the bevels run
/// on to the point their planes share, or at a convex vertex of four or
/// more sharing none until the others cut them; on curved faces to the
/// point their surfaces share, meeting along their sections from there.
/// The two differ only at a vertex where three or more bevels meet: a
/// request with no such vertex builds with either exactly what
/// [`chamfer_edges_with`](crate::chamfer_edges_with) builds.
///
/// # Errors
///
/// As [`chamfer_edges_with`](crate::chamfer_edges_with): a request the
/// engine does not build is refused by name.
pub fn chamfer_edges_cornered(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    corner: ChamferCorner,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if specs.is_empty() {
        ogeom_bail!(Construction, "a chain of no edges bevels nothing");
    }
    if let Some(built) = crate::compound::each_solid(
        model,
        solid,
        specs,
        |(edge, _)| edge,
        |model, part, own| chamfer_edges_cornered(model, part, own, corner, tol),
    )? {
        return Ok(built);
    }
    let (specs, _) = crate::chamfer::without_seams(model, solid, specs, tol)?;
    if specs.is_empty() {
        return Ok(crate::chamfer::unchanged(solid));
    }
    let specs = &specs[..];
    if crate::chamfer::below_a_floor(model, solid, specs, tol)?.is_some() {
        return Err(crate::chamfer::below_a_floor_refused());
    }
    if corner == ChamferCorner::Point && three_meet(model, specs)? {
        // The blend engine mitres in one pass, on planes or on curved
        // faces.
        return crate::chamfer::engine_chamfer_mitred(model, solid, specs, tol)?
            .ok_or_else(crate::chamfer::declined);
    }
    // Every vertex its facet. Where no three bevels meet, the mitre and
    // the facet are one corner, and `Point` builds what the default does.
    if let Some(built) = crate::chamfer::engine_chamfer(model, solid, specs, tol)? {
        return Ok(built);
    }
    Err(crate::chamfer::declined())
}

/// Whether three or more of `specs`' edges meet at one vertex: the only
/// place [`ChamferCorner::Point`] and [`ChamferCorner::Facet`] differ.
fn three_meet(model: &Model, specs: &[(Shape, Chamfer)]) -> OgeomResult<bool> {
    let mut count: ogeom_core::FastMap<ogeom_topo::TShapeId, usize> =
        ogeom_core::FastMap::default();
    for (edge, _) in specs {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
            continue;
        };
        let ends = if a.node() == b.node() {
            vec![a.node()]
        } else {
            vec![a.node(), b.node()]
        };
        for v in ends {
            *count.entry(v).or_default() += 1;
        }
    }
    Ok(count.values().any(|&n| n >= 3))
}
