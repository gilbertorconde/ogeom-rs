//! Arcs of one closed crease asked together (a pipe tee's junction loop the
//! boolean split into several edges): one round over the whole loop, as
//! tangent edges asked in one request are one feature. The march rounds the
//! loop from its first arc; the arcs after it lie on the crease that band
//! rounded and go with it.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of, explore_unique};

/// `once` rounded on `solid` where its edges, two or more at one radius,
/// run end to end round one closed crease between the same two faces: the
/// first edge's marched band, which runs the loop, every other edge gone
/// from the result and deleted in the history. `None` otherwise, or where
/// the band leaves any of the other edges standing.
///
/// # Errors
///
/// As reading the edges and the march.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let [(first, radius), rest @ ..] = once else {
        return Ok(None);
    };
    if rest.is_empty() || rest.iter().any(|(_, r)| r != radius) {
        return Ok(None);
    }
    // One pair of faces along every edge.
    let faces_of = |edge: &Shape| -> OgeomResult<Vec<Shape>> {
        let mut out: Vec<Shape> = Vec::new();
        for f in ancestors_of(model, solid, edge, ShapeType::Face)? {
            if !out.iter().any(|g| g.is_same(&f)) {
                out.push(f);
            }
        }
        Ok(out)
    };
    let pair = faces_of(first)?;
    if pair.len() != 2 {
        return Ok(None);
    }
    for (edge, _) in rest {
        let faces = faces_of(edge)?;
        if faces.len() != 2 || !faces.iter().all(|f| pair.iter().any(|g| g.is_same(f))) {
            return Ok(None);
        }
    }
    // A closed loop: every end met by exactly two of the edges.
    let mut ends: Vec<Shape> = Vec::new();
    for (edge, _) in once {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
            return Ok(None);
        };
        ends.push(a);
        ends.push(b);
    }
    if ends
        .iter()
        .any(|v| ends.iter().filter(|w| w.is_same(v)).count() != 2)
    {
        return Ok(None);
    }
    let mut scratch = model.clone();
    let Some(mut built) =
        super::marched_seats::rounds(&mut scratch, solid, &[(first.clone(), *radius)], false, tol)?
    else {
        return Ok(None);
    };
    let left = explore_unique(&scratch, &built.shape, ShapeType::Edge)?;
    for (edge, _) in rest {
        if !built.history.trace(edge).is_empty() || left.iter().any(|e| e.node() == edge.node()) {
            return Ok(None);
        }
        built.history.delete(edge);
    }
    *model = scratch;
    Ok(Some(built))
}
