//! Rounds and bevels in one request, closed as one corner where they meet.

use std::collections::BTreeMap;

use crate::setback::CornerEnd;
use crate::{Chamfer, ChamferCorner};
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Continuity;
use ogeom_topo::{Model, Shape, ShapeType, TShapeId, explore_unique};

/// How [`blend_edges`] blends one edge.
#[derive(Debug, Clone)]
pub enum Blend {
    /// A round of this radius, as [`fillet_edges_with`](crate::fillet_edges_with).
    Round(f64),
    /// A bevel, as [`chamfer_edges_with`](crate::chamfer_edges_with).
    Chamfer(Chamfer),
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
/// At a vertex where a round meets a bevel, every edge there must be in
/// the request, and the vertex a convex corner of planes: each band is
/// stopped at the least setback where its section clears every other band
/// and one of its own widths past it, the cell short of the sections cut
/// away (all that the order of the build leaves there), and the opening
/// closed by one patch. On each face between two bands it is trimmed by a
/// conic tangent to both bands' lines there. The patch is tangent to each
/// round's band and to each face, and meets each bevel's end in position
/// only: the bevel meets its faces in a crease, so where its end meets a
/// face the patch takes the face's tangent plane.
///
/// A vertex of three edges where one round meets one bevel and the third
/// edge is left sharp closes as
/// [`round_vertex_setback_sharp`](crate::round_vertex_setback_sharp)
/// closes a corner: the sharp edge is cut as far out as the farther band,
/// at `Q`, where a three-sided patch held across both of its faces turns
/// its normal from one face's to the other's, and one patch over the rest
/// meets the round's band and the faces between rounds tangent, and the
/// bevel's end and the trims from it in position only.
///
/// # Errors
///
/// As the rounds' and the bevels' own calls, and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// an edge is asked twice, or a round meets a bevel at a vertex where an
/// edge is left sharp other than one beside one round and one bevel; or
/// where that vertex is not a convex corner of planes.
pub fn blend_edges(
    model: &mut Model,
    solid: &Shape,
    edges: &[(Shape, Blend)],
    corner: ChamferCorner,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut rounds: Vec<(Shape, f64)> = Vec::new();
    let mut bevels: Vec<(Shape, Chamfer)> = Vec::new();
    for (k, (edge, blend)) in edges.iter().enumerate() {
        if edges[..k].iter().any(|(e, _)| e.is_same(edge)) {
            ogeom_bail!(Construction, "an edge is asked to be blended twice");
        }
        match blend {
            Blend::Round(r) => rounds.push((edge.clone(), *r)),
            Blend::Chamfer(c) => bevels.push((edge.clone(), c.clone())),
        }
    }
    if rounds.is_empty() {
        return crate::chamfer_edges_cornered(model, solid, &bevels, corner, tol);
    }
    if bevels.is_empty() {
        return crate::fillet_edges_with(model, solid, &rounds, tol);
    }
    rounds.sort_by_key(|(e, _)| e.node());
    bevels.sort_by_key(|(e, _)| e.node());

    // The vertices where a round meets a bevel, each with the request's
    // edges there.
    let mut at: BTreeMap<TShapeId, (Shape, Vec<usize>, Vec<usize>)> = BTreeMap::new();
    for (i, (edge, _)) in rounds.iter().enumerate() {
        for v in ends(model, edge)? {
            at.entry(v.node())
                .or_insert_with(|| (v.clone(), Vec::new(), Vec::new()))
                .1
                .push(i);
        }
    }
    for (i, (edge, _)) in bevels.iter().enumerate() {
        for v in ends(model, edge)? {
            at.entry(v.node())
                .or_insert_with(|| (v.clone(), Vec::new(), Vec::new()))
                .2
                .push(i);
        }
    }
    let solid_edges = explore_unique(model, solid, ShapeType::Edge)?;
    let mut mixed: Vec<(Shape, Vec<usize>, Vec<usize>)> = Vec::new();
    for (vertex, r, b) in at.into_values() {
        if r.is_empty() || b.is_empty() {
            continue;
        }
        let mut meeting = 0;
        for e in &solid_edges {
            if ends(model, e)?.iter().any(|v| v.is_same(&vertex)) {
                meeting += 1;
            }
        }
        // One edge left sharp beside one round and one bevel: the opening
        // crosses it at a point, where the patch's normal turns.
        let one_sharp = meeting == 3 && r.len() == 1 && b.len() == 1;
        if meeting != r.len() + b.len() && !one_sharp {
            ogeom_bail!(
                Construction,
                "a round meets a bevel at a vertex where an edge is left sharp: the \
                 patch closing the round's opening would cross the bevel's sharp end \
                 line on the face beside it and be tangent there to two planes; blend \
                 every edge at the vertex, or leave one sharp beside one round and one \
                 bevel"
            );
        }
        mixed.push((vertex, r, b));
    }

    // Each mixed corner's frame and setbacks, read off the solid as it
    // stands.
    let mut corners = Vec::with_capacity(mixed.len());
    for (vertex, r, b) in &mixed {
        let largest = r.iter().map(|&i| rounds[i].1).fold(0.0, f64::max);
        let frame = crate::corner::corner_frame(model, solid, vertex, largest, tol)?;
        let n = frame.d.len();
        let mut widths = Vec::with_capacity(n);
        let mut kinds = Vec::with_capacity(n);
        let mut sharp: Option<usize> = None;
        for k in 0..n {
            let (ma, mb) = (frame.m[frame.ring[k]], frame.m[frame.ring[(k + 1) % n]]);
            let dk = frame.d[k].vector();
            let along = |edge: &Shape| -> OgeomResult<bool> {
                let Some((a, z)) = ogeom_algo::edge_vertices(model, edge)? else {
                    return Ok(false);
                };
                let far = if a.is_same(vertex) { z } else { a };
                let p = point_of(model, &far)? - frame.corner;
                let length = p.magnitude();
                Ok(length > tol.confusion() && p.dot(dk) / length >= 1.0 - tol.angular().sqrt())
            };
            let mut found = None;
            for &i in r {
                if along(&rounds[i].0)? {
                    let cos = -ma.dot(mb);
                    let sin = (1.0 - cos * cos).max(0.0).sqrt();
                    if sin <= tol.angular() {
                        ogeom_bail!(
                            Construction,
                            "two faces at a corner meet flat; there is no edge to round"
                        );
                    }
                    let w = rounds[i].1 * (1.0 + cos) / sin;
                    found = Some(((w, w), CornerEnd::Round));
                }
            }
            for &i in b {
                if along(&bevels[i].0)? {
                    let (edge, spec) = &bevels[i];
                    let crate::chamfer::Legs::Straight(seat, d) =
                        crate::chamfer::legs(model, solid, edge, spec, tol)?
                    else {
                        ogeom_bail!(
                            Construction,
                            "a corner of rounds and bevels takes bevels on straight edges \
                             between planes"
                        );
                    };
                    let on = |m: ogeom_math::Vector| {
                        (0..2)
                            .find(|&j| seat.normals[j].cross(m).magnitude() <= tol.angular() * 10.0)
                    };
                    let (Some(ja), Some(jb)) = (on(ma), on(mb)) else {
                        ogeom_bail!(
                            Construction,
                            "a bevel at a corner does not lie on the corner's faces"
                        );
                    };
                    found = Some(((d[ja], d[jb]), CornerEnd::Bevel));
                }
            }
            let Some((w, kind)) = found.or_else(|| {
                // The one edge left sharp: a band of no width.
                (sharp.is_none() && n == 3 && r.len() + b.len() == 2).then(|| {
                    sharp = Some(k);
                    ((0.0, 0.0), CornerEnd::Round)
                })
            }) else {
                ogeom_bail!(
                    Construction,
                    "an edge at a corner of rounds and bevels is not one of the request's, \
                     or is curved"
                );
            };
            widths.push(w);
            kinds.push(kind);
        }
        let mut setback = crate::setback::least_setbacks_widths(&frame, &widths, tol)?;
        if let Some(ks) = sharp {
            // The sharp edge cut as far out as the farther of the bands
            // beside it, clear of their corner.
            setback[ks] = (0..n)
                .filter(|&k| k != ks)
                .map(|k| setback[k])
                .fold(0.0, f64::max);
        }
        corners.push((vertex.clone(), frame, setback, kinds, sharp));
    }

    // The bevels, then the rounds on what they leave.
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
    let mut built = Built {
        shape: rounded.shape.clone(),
        history: beveled.history.then(&rounded.history),
    };
    for (vertex, frame, setback, kinds, sharp) in &corners {
        built = match sharp {
            Some(ks) => crate::sharp::close_sharp_corner(
                model, built, vertex, frame, setback, *ks, kinds, tol,
            )?,
            None => crate::setback::close_corner_ends(
                model,
                built,
                vertex,
                frame,
                setback,
                kinds,
                Continuity::G1,
                tol,
            )?,
        };
    }
    Ok(built)
}

/// An edge's two vertices.
fn ends(model: &Model, edge: &Shape) -> OgeomResult<Vec<Shape>> {
    Ok(match ogeom_algo::edge_vertices(model, edge)? {
        Some((a, b)) => vec![a, b],
        None => Vec::new(),
    })
}

/// A vertex's point, placed.
fn point_of(model: &Model, vertex: &Shape) -> OgeomResult<ogeom_math::Point> {
    let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
        ogeom_bail!(Dangling, "a vertex is not in this model");
    };
    Ok(vertex.transform(model.datums())?.apply(data.point))
}
