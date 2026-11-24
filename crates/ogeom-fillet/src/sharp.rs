//! A setback corner with one edge left sharp: the bands of the vertex's
//! other edges stopped at their setbacks, and the opening between them
//! closed by a patch that runs into a point on the sharp edge.
//!
//! The opening crosses the sharp edge, and where it does, at `Q`, the
//! patch is to be tangent to both of that edge's faces: two tangent planes
//! at one point, which no patch with a tangent plane there has. The patch
//! at `Q` is three-sided, its fourth side shrunk to `Q`, and its normal
//! turns along that shrunk side from one face's to the other's, as a
//! cone's turns round its apex: tangent to both faces up to `Q` and,
//! locally, a cone at `Q` fanning between them.
//!
//! The sharp edge is read as a band of no width: its end on its section
//! plane is the point `Q`, and the trims on its two faces run from `Q`
//! into the bands' rails, each heading from `Q` for its band's rail half
//! way back to the vertex, so the two trims at `Q` part and the cone there
//! opens.
//!
//! The patch at `Q` runs out to the middle of each of those trims, where a
//! curve between them bounds it; it is fitted by
//! [`ogeom_fit::fit_network`], its halves of the trims held to them and
//! held across the faces, so it is tangent to both faces exactly up to the
//! held curves' own projection. The rest of the opening, six sides with
//! that curve, is filled by one patch G1 to every side, the patch at `Q`
//! among its supports.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Continuity;
use ogeom_topo::{Model, Shape};

use crate::setback::{Setbacks, corner_edges};

/// Round the edges of `vertex` other than `sharp` with bands at `radius`,
/// each stopped at its setback, and close the opening between them with a
/// patch that runs into the point `sharp`'s setback names on it.
///
/// The vertex is one of three planar faces and three straight edges, a
/// convex corner one ball touches every face of (a box's corner). The two
/// bands are cut square to their edges at their setbacks; the sharp edge
/// is cut at its own setback, at `Q`. On the face between the two bands
/// the trim is the conic tangent to both rails, as
/// [`round_vertex_setback`](crate::round_vertex_setback) trims; on each face
/// of the sharp edge the trim is the conic from `Q` tangent to the band's
/// rail, heading from `Q` for that rail half way back to the vertex.
///
/// The opening is closed by two patches. At `Q`, a three-sided one bounded
/// by the halves of the two trims from `Q` and a curve between their
/// middles, its fourth side shrunk to `Q`: fitted by
/// [`ogeom_fit::fit_network`] to arcs swept from `Q` that leave each trim
/// square across it in its face, its halves of the trims held to them and
/// held across the two faces, so it is tangent to each face exactly up to
/// the held curves' own projection, its normal turning at `Q` from one
/// face's to the other's. Over the rest, one patch from
/// [`make_filling_n`](crate::make_filling_n) meeting the two bands' end
/// arcs, the trims and the patch at `Q` at G1, to the bands' fit target.
/// Where it misses, or the patch at `Q` folds its normal back at `Q`, the
/// arcs are drawn fuller or flatter and the two patches fitted again.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// `radius` or a setback is not a positive length, a setback is shorter
/// than the ball's own along a rounded edge or reaches past an edge's far
/// end, the sharp edge's setback lies inside the bands' own corner, the
/// continuity is not G1, `sharp` does not leave the vertex, or the vertex
/// is not a convex planar corner of three edges one ball touches every
/// face of; [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where,
/// however its arcs are drawn, the rest misses the bands' fit target or
/// the patch at `Q` folds its normal back.
pub fn round_vertex_setback_sharp(
    model: &mut Model,
    solid: &Shape,
    vertex: &Shape,
    sharp: &Shape,
    radius: f64,
    setbacks: &Setbacks,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut own: Vec<Shape> = vec![vertex.clone(), sharp.clone()];
    own.extend(setbacks.per_edge.iter().map(|(e, _)| e.clone()));
    ogeom_algo::on_own_nodes(model, solid, &own, |model, solid, parts| {
        let named: Vec<(Shape, f64)> = parts[2..]
            .iter()
            .zip(&setbacks.per_edge)
            .map(|(e, (_, s))| (e.clone(), *s))
            .collect();
        let own = Setbacks {
            per_edge: named,
            default: setbacks.default,
            continuity: setbacks.continuity,
        };
        sharp_own(model, solid, &parts[0], &parts[1], radius, &own, tol)
    })
}

/// [`round_vertex_setback_sharp`] on a shape every node below which is its
/// own.
fn sharp_own(
    model: &mut Model,
    solid: &Shape,
    vertex: &Shape,
    sharp: &Shape,
    radius: f64,
    setbacks: &Setbacks,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(Construction, "a blend radius must be a positive distance");
    }
    if setbacks.continuity != Continuity::G1 {
        ogeom_bail!(
            Construction,
            "a setback corner with an edge left sharp meets its bands and faces \
             at G1, not {:?}",
            setbacks.continuity
        );
    }
    let frame = crate::engine::vertex_read::corner_frame(model, solid, vertex, radius, tol)?;
    let n = frame.d.len();
    let (edges, reach) = corner_edges(model, solid, vertex, &frame, tol)?;
    let Some(ks) = edges.iter().position(|e| e.is_same(sharp)) else {
        ogeom_bail!(
            Construction,
            "the edge left sharp does not leave the vertex"
        );
    };
    if n != 3 {
        ogeom_bail!(
            Construction,
            "a setback corner with an edge left sharp closes a vertex of three \
             edges; one of {n} is still owed; see docs/CAPABILITIES.md, \
             fillet.setback-corners"
        );
    }
    let setback: Vec<f64> = edges
        .iter()
        .map(|e| {
            setbacks
                .per_edge
                .iter()
                .find(|(named, _)| named.is_same(e))
                .map_or(setbacks.default, |(_, s)| *s)
        })
        .collect();
    for k in 0..n {
        let s = setback[k];
        if !s.is_finite() || s <= tol.confusion() {
            ogeom_bail!(Construction, "a setback of {s} is not a positive distance");
        }
        let own = (frame.far - frame.corner).dot(frame.d[k].vector());
        if k == ks && s <= own + tol.confusion() {
            ogeom_bail!(
                Construction,
                "the edge left sharp is cut at {s}, inside the bands' own \
                 corner, which reaches {own} along it; cut it farther from \
                 the vertex"
            );
        }
        if k != ks && s < own - tol.confusion() {
            ogeom_bail!(
                Construction,
                "a setback of {s} stops the band inside the ball's own corner, \
                 which reaches {own} along its edge"
            );
        }
        if s >= reach[k] - tol.confusion() {
            ogeom_bail!(
                Construction,
                "a setback of {s} reaches past its edge, {} long",
                reach[k]
            );
        }
    }
    crate::engine::vertex_setbacks::guard(model, solid, vertex, &frame, &setback, tol)?;
    crate::engine::vertex_setbacks::sharp(model, solid, vertex, &edges, &setback, ks, radius, tol)
}
