//! A round of an edge on a marched seat that takes a narrow strip beside it
//! off whole: the marched counterpart of [`crate::consumed_rim`].
//!
//! Where a strip lies between the edge's face and a face beyond it (a
//! bevelled or small-rounded rim where a bore leaves a drum, rounded
//! larger), no ball touches both of the edge's faces. The ball that rounds
//! the corner then touches the two faces either side of the strip, rolling
//! along the crossing of their surfaces nearest the strip, as
//! [`crate::pair_marched`] rolls it between two faces that share no edge.
//! The strip goes with the corner where it lies wholly in it, over the
//! whole loop of a closed seat or over the runs of an open one, each
//! capped in its end sections' planes:
//!
//! - the ball touches both outer faces, past each strip edge (the march
//!   refuses a touch off its face);
//! - no face of the solid reaches into the ball's lune at any station, the
//!   strip included ([`crate::pair_marched::face_inside`]);
//! - every point of the strip's mesh stands behind both outer surfaces, in
//!   the corner they make, and between some run's two caps;
//! - no edge of the solid but the strip's own and the seams stands within
//!   the corner's reach of the crossing, strictly between some run's caps:
//!   the farthest a touch stands from its station's crossing point.

use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_geom::Curve3d as _;
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

/// The faces of `solid` holding `edge`.
fn faces_at(model: &Model, solid: &Shape, edge: &Shape) -> OgeomResult<Vec<Shape>> {
    let mut found = Vec::new();
    for face in explore_unique(model, solid, ShapeType::Face)? {
        if explore_unique(model, &face, ShapeType::Edge)?
            .iter()
            .any(|e| e.is_same(edge))
        {
            found.push(face);
        }
    }
    Ok(found)
}

/// The round of `edge` of `solid` at `radius` taking a strip beside it off
/// whole, where one of the edge's faces is such a strip between the other
/// and a face beyond it on a marched seat; `None` otherwise, the caller's
/// refusal standing.
pub(crate) fn consumed_marched(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    crate::engine::old_marched_constructor("consumed_marched");
    let beside = faces_at(model, solid, edge)?;
    let [one, two] = beside.as_slice() else {
        return Ok(None);
    };
    for (strip, outer) in [(one, two), (two, one)] {
        let Some(beyond) = strip_beyond(model, solid, strip, outer)? else {
            continue;
        };
        match round_over(model, solid, edge, strip, [outer, &beyond], radius, tol) {
            Ok(Some(built)) => return Ok(Some(built)),
            Ok(None) | Err(OgeomError::Construction(_) | OgeomError::NotDone(_)) => {}
            Err(e) => return Err(e),
        }
    }
    Ok(None)
}

/// The face across `strip` from `outer`: every edge of the strip but its
/// seams is shared with `outer`, with that face, or with a face at one of
/// the strip's ends, and the face across is the one sharing the most of
/// the strip's border after `outer`. Both sides hold some.
fn strip_beyond(
    model: &Model,
    solid: &Shape,
    strip: &Shape,
    outer: &Shape,
) -> OgeomResult<Option<Shape>> {
    let mut shares: Vec<(Shape, f64)> = Vec::new();
    let mut on_outer = false;
    for e in explore_unique(model, strip, ShapeType::Edge)? {
        let holders = faces_at(model, solid, &e)?;
        let others: Vec<&Shape> = holders.iter().filter(|f| !f.is_same(strip)).collect();
        match others.as_slice() {
            // A seam of the strip itself.
            [] => {}
            [f] if f.is_same(outer) => on_outer = true,
            [f] => {
                let length = edge_length(model, &e)?;
                match shares.iter_mut().find(|(g, _)| g.is_same(f)) {
                    Some((_, total)) => *total += length,
                    None => shares.push(((*f).clone(), length)),
                }
            }
            _ => return Ok(None),
        }
    }
    if !on_outer {
        return Ok(None);
    }
    Ok(shares
        .into_iter()
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(f, _)| f))
}

/// An edge's length, summed over chords.
fn edge_length(model: &Model, edge: &Shape) -> OgeomResult<f64> {
    let Some((curve, (a, b))) = ogeom_algo::placed_edge_curve(model, edge, Tolerances::default())?
    else {
        return Ok(0.0);
    };
    let tol = Tolerances::default();
    let mut length = 0.0;
    let mut last = curve.point_at(a, tol)?;
    for k in 1..=64 {
        let p = curve.point_at((b - a).mul_add(f64::from(k) / 64.0, a), tol)?;
        length += p.distance(last);
        last = p;
    }
    Ok(length)
}

/// The round of the corner of `outer`, the two faces either side of
/// `strip`, along their crossing nearest `edge`, applied where the strip
/// lies wholly in the corner it cuts off.
fn round_over(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    _strip: &Shape,
    outer: [&Shape; 2],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let (curve, range) = crate::support::edge_curve(model, edge, tol)?;
    let near = curve.point_at(f64::midpoint(range.0, range.1), tol)?;
    // The round is checked whole before it is taken, so its lines of contact
    // may be held to the round's own fit where the tighter one misses.
    let hosts = crate::pair_marched::hosts_of(model, outer, radius, tol)?.loosened();
    // The ball rolls behind both faces: a convex corner, cut off.
    let seat = crate::pair_marched::pair_seat_near(
        model,
        outer,
        &hosts,
        radius,
        true,
        Some(solid),
        Some(near),
        tol,
    )?;
    if crate::pair_marched::face_inside(model, solid, outer, &seat, radius, tol)?.is_some() {
        return Ok(None);
    }
    // The ball rolls in the material: beside its arc's middle, the solid.
    for &[_, lune] in seat.probes() {
        if !matches!(
            ogeom_algo::classify_in_solid_exact(model, solid, lune, tol)?,
            ogeom_algo::Containment::In
        ) {
            return Ok(None);
        }
    }
    // A strip the ball's lunes keep clear is taken off whole by the
    // engine's marched seats, or refused by name.
    Err(crate::fillet::declined())
}
