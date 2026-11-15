//! A round wider than a narrow face beside its edge takes that face off
//! whole (ADR 0001, amendment "M5d").
//!
//! At a straight edge `E` between planar faces `F` and `N`, where `N` is
//! narrower than the ball's setback on it, no ball touches both of `E`'s
//! faces: the ball touching `F` and the face `G` beyond `N`'s far edge
//! holds `N` in the corner it takes off (a convex edge) or under the fill
//! it adds (a concave one). The round is then the round of the crease `F`
//! and `G` close in once `N` is gone: `N` comes off a copy of the solid
//! ([`ogeom_bool::remove_faces`], its neighbours closing over it), that
//! crease is rounded there by the engine, and the copy less its rounding
//! and the rounding less the copy are applied to the solid as given, as
//! the general non-local route applies them (amendment "M5c"). Built only
//! where `N` is gone from the result, wholly taken; anything else (a ball
//! that would still touch `N`) is handed back.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, SurfaceGeometry};
use ogeom_math::Point;
use ogeom_topo::{Model, Shape};

use super::curved::read;
use super::defeatured::Rounder;
use super::far_edge::{far_edge, into_face};
use super::trimmed::{faces_of, outward};
use super::{Id, id};

/// The rounds of `asked` on `solid` where at least one edge's ball is wider
/// than a narrow planar face beside it, that face taken off whole and the
/// edge rounded as the crease its neighbours close in, the rest as asked,
/// all by `round` on a copy of the solid without those faces. `None` where
/// no edge is of that kind, or a narrow face is not wholly taken.
///
/// # Errors
///
/// As rounding the copy, other than refusals, which decline.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    asked: &[(Shape, f64)],
    round: Rounder,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    let mut narrow: Vec<usize> = Vec::new();
    let mut across: Vec<(Shape, [Shape; 2])> = Vec::new();
    for (edge, radius) in asked {
        let e = id(edge);
        if !curved.edges.contains_key(&e) {
            return Ok(None);
        }
        if let Some((n, [f, g])) = narrow_beside(model, &curved, e, *radius, tol)? {
            if !narrow.contains(&n) {
                narrow.push(n);
            }
            across.push((
                edge.clone(),
                [curved.faces[f].shape.clone(), curved.faces[g].shape.clone()],
            ));
        }
    }
    if narrow.is_empty() {
        return Ok(None);
    }
    // The faces the narrow ones meet stay: only the narrow faces come off.
    let far: Vec<Shape> = narrow
        .iter()
        .map(|n| curved.faces[*n].shape.clone())
        .collect();
    let mut near: Vec<Shape> = Vec::new();
    for n in &narrow {
        for (x, _) in curved.faces[*n].loops.iter().flatten() {
            for f in faces_of(&curved, *x) {
                if !narrow.contains(&f) && !near.iter().any(|s| s.is_same(&curved.faces[f].shape)) {
                    near.push(curved.faces[f].shape.clone());
                }
            }
        }
    }
    let mut scratch = model.clone();
    let Some(tools) = super::defeatured::tools_across(
        &mut scratch,
        solid,
        asked,
        &far,
        &near,
        &across,
        round,
        tol,
    )?
    else {
        return Ok(None);
    };
    let edges: Vec<Shape> = asked.iter().map(|(e, _)| e.clone()).collect();
    let Some(built) = super::fallback::applied(&mut scratch, solid, &edges, &tools, tol)? else {
        return Ok(None);
    };
    // Each narrow face wholly taken: nothing of it is left.
    if far.iter().any(|n| !built.history.trace(n).is_empty()) {
        return Ok(None);
    }
    match super::guards::solid_fits(&scratch, &built.shape, tol) {
        Ok(()) => {
            *model = scratch;
            Ok(Some(built))
        }
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(None),
    }
}

/// The narrow face beside straight `edge` that its ball at `radius` is
/// wider than, with the edge's other face and the face beyond the narrow
/// one's far edge: `None` where neither of its planar faces is such a
/// strip (its first edge met crossing it square from the edge's middle
/// straight, nearer than the setback).
fn narrow_beside(
    model: &Model,
    curved: &super::curved::Curved,
    edge: Id,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(usize, [usize; 2])>> {
    let eps = tol.confusion();
    let c = &curved.edges[&edge];
    if c.ends.0 == c.ends.1 || !matches!(c.curve, Curve::Line(_)) {
        return Ok(None);
    }
    let (a, b) = (curved.points[&c.ends.0], curved.points[&c.ends.1]);
    let length = a.distance(b);
    if length <= eps {
        return Ok(None);
    }
    let along = (b - a) / length;
    let mid: Point = a.midpoint(b);
    let [fa, fb] = faces_of(curved, edge)[..] else {
        return Ok(None);
    };
    if [fa, fb]
        .iter()
        .any(|f| !matches!(curved.faces[*f].surface, SurfaceGeometry::Plane(_)))
    {
        return Ok(None);
    }
    let (Some(na), Some(nb)) = (outward(curved, fa, mid), outward(curved, fb, mid)) else {
        return Ok(None);
    };
    let turn = na.dot(nb).clamp(-1.0, 1.0).acos();
    let setback = radius * (0.5 * turn).tan();
    for (n, f) in [(fa, fb), (fb, fa)] {
        let Some(into) = into_face(model, curved, n, (mid, along, length), tol)? else {
            continue;
        };
        let Some((far, at, _)) = far_edge(curved, n, edge, mid, into, tol)? else {
            continue;
        };
        if (at - mid).magnitude() >= setback - eps {
            continue;
        }
        let beyond = faces_of(curved, far);
        let [ga, gb] = beyond[..] else {
            continue;
        };
        let g = if ga == n { gb } else { ga };
        if g == f {
            continue;
        }
        return Ok(Some((n, [f, g])));
    }
    Ok(None)
}
