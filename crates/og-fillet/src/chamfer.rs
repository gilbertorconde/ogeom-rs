//! Chamfers: the bevel that replaces an edge.
//!
//! P4's opening stone, built deliberately on M3's shoulders: a chamfer along
//! a straight edge between two planar faces is a wedge subtracted, and the
//! wedge's own faces lie *exactly* on the solid's — coplanar, materials
//! aligned — which is the same-domain case the boolean learned to resolve.
//! The blend machinery proper (rolling-ball fillets, variable radii, vertex
//! blends) comes after; a chamfer is the member of the family that needs no
//! new surface, only the machinery already proven.

use crate::support::{planar_face, planar_seat, subtract_wedge};
use og_algo::Built;
use og_core::{OgResult, Tolerances, og_bail};
use og_topo::{Model, Shape};

/// Bevel a straight edge of a solid, cutting `distance` back along each of
/// its two faces.
///
/// The edge must be straight, convex, and shared by exactly two planar faces;
/// the distances are equal (the symmetric chamfer). The result is the boolean
/// difference with a wedge whose legs run along the two faces — so the
/// history reads as a cut: the two faces are modified into their trimmed
/// pieces, the edge's neighbourhood gains the bevel face.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the edge is
/// not straight, not convex, not shared by exactly two planar faces of
/// `solid`, or `distance` is not a usable length.
pub fn chamfer_edge(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    distance: f64,
    tol: Tolerances,
) -> OgResult<Built> {
    if !distance.is_finite() || distance <= tol.confusion() {
        og_bail!(Construction, "a chamfer of {distance} cuts nothing");
    }
    let seat = planar_seat(model, solid, edge, tol)?;
    let a = seat.leg(0, tol)?;
    let b = seat.leg(1, tol)?;

    let travel = seat.end - seat.start;
    let apex0 = seat.start;
    let apex1 = seat.end;
    let a0 = apex0 + a * distance;
    let b0 = apex0 + b * distance;
    let a1 = a0 + travel;
    let b1 = b0 + travel;

    // The bevel's outward normal is the leg bisector exactly: the bevel plane
    // holds `a - b` and the edge direction, and the bisector is perpendicular
    // to both, pointing from the apex toward the cut.
    let bevel_out = {
        let u = a + b;
        let m = u.magnitude();
        if m <= tol.angular() {
            og_bail!(Construction, "the faces meet too sharply to seat a chamfer");
        }
        u / m
    };

    // The wedge: a triangular prism whose apex line is the edge and whose
    // legs run `distance` along each face. Built from five explicit planar
    // faces rather than swept, because a sweep's walls are extrusion
    // surfaces even when they are geometrically planes, and the boolean's
    // same-domain resolution — which is what makes the coplanar legs melt
    // into the solid's own faces — recognises coincidence between *planes*.
    let faces = [
        planar_face(model, &[apex0, a0, b0], -seat.along, tol)?,
        planar_face(model, &[apex1, a1, b1], seat.along, tol)?,
        planar_face(model, &[apex0, a0, a1, apex1], seat.normals[0], tol)?,
        planar_face(model, &[apex0, b0, b1, apex1], seat.normals[1], tol)?,
        planar_face(model, &[a0, b0, b1, a1], bevel_out, tol)?,
    ];
    subtract_wedge(model, solid, edge, &faces, tol)
}
