//! The fit and overlap guards a request's balls pass before any band is
//! built: each ball sits on both faces of its edge, a ball resting on a
//! narrow face's far edge is read as such, and no band overlaps another
//! across a face they share. The engine's marched seats, the blame and
//! largest-fit searches and the one-edge entries all ask these.

use ogeom_core::{OgeomError, Tolerances};
use ogeom_topo::{Model, Shape, ShapeType};

/// The checks a request's bands pass before any is built: each ball sits
/// on both faces of its edge, and no band overlaps another across a face
/// they share. Which edges, sharing nothing with the rest of the request,
/// had their ball rest on the far edge of a narrow face beside them, to
/// be rounded as alone, comes back.
///
/// # Errors
///
/// [`OgeomError::Construction`] where a ball does not sit on its faces or
/// two bands overlap, naming the edges refused.
pub(crate) fn request_fits(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    radii: &[f64],
    tol: Tolerances,
) -> Result<Vec<bool>, crate::support::Blamed> {
    use crate::support::Blamed;
    let blame = |i: usize| {
        move |error: OgeomError| Blamed {
            edges: vec![i],
            error,
        }
    };
    // The ball touches each face a setback back from the edge; set back
    // past a face's far side, the band would cut through the face, and set
    // back into another edge's band, the two balls overlap.
    let mut contacts = Vec::with_capacity(edges.len());
    // An edge sharing no face and no vertex with another of the request
    // stands alone: no corner or junction takes its ends, and its band
    // is built as it is alone.
    let mut touched: Vec<(Vec<ogeom_topo::TShapeId>, Vec<ogeom_topo::TShapeId>)> =
        Vec::with_capacity(edges.len());
    for edge in edges {
        let faces = ogeom_topo::ancestors_of(model, solid, edge, ShapeType::Face)?
            .iter()
            .map(Shape::node)
            .collect();
        let ends = ogeom_algo::edge_vertices(model, edge)?
            .map(|(a, b)| vec![a.node(), b.node()])
            .unwrap_or_default();
        touched.push((faces, ends));
    }
    let alone = |i: usize| {
        let (faces, ends) = &touched[i];
        touched.iter().enumerate().all(|(j, (f, e))| {
            j == i || (!f.iter().any(|x| faces.contains(x)) && !e.iter().any(|x| ends.contains(x)))
        })
    };
    // Edges sharing nothing with the rest of the request whose ball rests
    // on the far edge of a narrow face beside them, rounded as alone.
    let mut rolled = vec![false; edges.len()];
    crate::support::keeping_face_meshes(|| -> Result<(), Blamed> {
        // A band checked against no other is held to what one edge alone
        // is.
        for (i, edge) in edges.iter().enumerate() {
            if let Some(run) =
                crate::shelf::run(model, solid, edge, radii[i], tol).map_err(blame(i))?
            {
                contacts.push(run.contacts(tol).map_err(blame(i))?);
                rolled[i] = true;
                continue;
            }
            if let Some(run) =
                crate::rolled_rim::rim_run(model, solid, edge, radii[i], tol).map_err(blame(i))?
            {
                contacts.push(run.contacts());
                rolled[i] = true;
                continue;
            }
            match crate::support::ball_fits(model, solid, edge, radii[i], !alone(i), tol) {
                Ok(found) => contacts.push(found),
                Err(refusal @ OgeomError::Construction(_)) if edges.len() > 1 && alone(i) => {
                    if crate::rolled::rolled_run(model, solid, edge, radii[i], tol)
                        .map_err(blame(i))?
                        .is_none()
                    {
                        return Err(blame(i)(refusal));
                    }
                    rolled[i] = true;
                    contacts.push(Vec::new());
                }
                Err(other) => return Err(blame(i)(other)),
            }
        }
        // A band lies within its setback of its edge, its contacts its
        // farthest points. A band standing alone whose ball runs off a face
        // is cut off on the face across that face's side, which may be one
        // another band stands on; where the two edges stand within their
        // setbacks' sum, it is checked against the others as a neighbour's
        // is.
        let reach =
            |found: &[crate::support::Contact]| found.iter().map(|c| c.setback).reduce(f64::max);
        let mut boxes = Vec::with_capacity(edges.len());
        for edge in edges {
            boxes.push(ogeom_algo::shape_bounds(model, edge, tol)?);
        }
        let reaches: Vec<Option<f64>> = contacts.iter().map(|found| reach(found)).collect();
        for (i, edge) in edges.iter().enumerate() {
            // A carried run already supplied its contacts on every support,
            // including the face across the narrow strip. Asking the ordinary
            // two-face ball again would discard that transfer. The same
            // bands_clear check below still compares all its contacts.
            if rolled[i] {
                continue;
            }
            let Some(own) = reaches[i].filter(|_| alone(i)) else {
                continue;
            };
            let near = reaches.iter().enumerate().any(|(j, other)| {
                j != i
                    && other.is_some_and(|other| boxes[i].distance_to_box(&boxes[j]) < own + other)
            });
            if near {
                contacts[i] = crate::support::ball_fits(model, solid, edge, radii[i], true, tol)
                    .map_err(blame(i))?;
            }
        }
        Ok(())
    })?;
    crate::support::bands_clear(model, edges, &contacts, tol)?;
    Ok(rolled)
}
