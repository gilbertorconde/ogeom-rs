//! A round of an edge on a marched seat that takes a narrow strip beside it
//! off whole: the marched counterpart of [`crate::consumed_rim`].
//!
//! Where a strip lies between the edge's face and a face beyond it (a
//! bevelled or small-rounded rim where a bore leaves a drum, rounded
//! larger), no ball touches both of the edge's faces. The ball that rounds
//! the corner then touches the two faces either side of the strip, rolling
//! along the crossing of their surfaces nearest the strip, as
//! [`crate::pair_marched`] rolls it between two faces that share no edge.
//! The strip goes with the corner where it lies wholly in it:
//!
//! - the ball touches both outer faces, past each strip edge (the march
//!   refuses a touch off its face);
//! - no face of the solid reaches into the ball's lune at any station, the
//!   strip included ([`crate::pair_marched::face_inside`]);
//! - every point of the strip's mesh stands behind both outer surfaces, in
//!   the corner they make;
//! - no edge of the solid but the strip's own and the seams stands within
//!   the corner's reach of the crossing: the farthest a touch stands from
//!   its station's crossing point.

use ogeom_algo::Built;
use ogeom_core::{OgeomError, OgeomResult, Tolerances};
use ogeom_geom::{Curve3d as _, Surface as _};
use ogeom_math::Point;
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
/// seams is shared with `outer` or with that one face, and both sides hold
/// some.
fn strip_beyond(
    model: &Model,
    solid: &Shape,
    strip: &Shape,
    outer: &Shape,
) -> OgeomResult<Option<Shape>> {
    let mut beyond: Option<Shape> = None;
    let mut on_outer = false;
    for e in explore_unique(model, strip, ShapeType::Edge)? {
        let holders = faces_at(model, solid, &e)?;
        let others: Vec<&Shape> = holders.iter().filter(|f| !f.is_same(strip)).collect();
        match others.as_slice() {
            // A seam of the strip itself.
            [] => {}
            [f] if f.is_same(outer) => on_outer = true,
            [f] => match &beyond {
                Some(b) if !b.is_same(f) => return Ok(None),
                Some(_) => {}
                None => beyond = Some((*f).clone()),
            },
            _ => return Ok(None),
        }
    }
    Ok(beyond.filter(|b| on_outer && !b.is_same(outer)))
}

/// The round of the corner of `outer`, the two faces either side of
/// `strip`, along their crossing nearest `edge`, applied where the strip
/// lies wholly in the corner it cuts off.
fn round_over(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    strip: &Shape,
    outer: [&Shape; 2],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let (curve, range) = crate::support::edge_curve(model, edge, tol)?;
    let near = curve.point_at(f64::midpoint(range.0, range.1), tol)?;
    let hosts = crate::pair_marched::hosts_of(model, outer, radius, tol)?;
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
    if !seat.closed() {
        return Ok(None);
    }
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
    // The corner's reach: the farthest a touch stands from its station's
    // point on the crossing.
    let guide = seat.guide().clone();
    let mut reach: f64 = 0.0;
    for blend in seat.runs() {
        for i in 0..blend.len() {
            let crease = guide.point_at(blend.along[i], tol)?;
            reach = reach
                .max(crease.distance(blend.touch_first[i]))
                .max(crease.distance(blend.touch_second[i]));
        }
    }
    if !strip_behind(model, strip, outer, reach, tol)? {
        return Ok(None);
    }
    if other_edge_within(model, solid, strip, &guide, reach, tol)? {
        return Ok(None);
    }
    let mut built = crate::pair_marched::apply_to_solid(model, solid, seat, radius, true, tol)?;
    built.history.delete(edge);
    crate::support::credit_new_faces(model, solid, edge, &mut built)?;
    Ok(Some(built))
}

/// Whether every point of the strip's mesh stands behind both outer faces'
/// surfaces, within its edges' tolerances: in the corner they make.
fn strip_behind(
    model: &Model,
    strip: &Shape,
    outer: [&Shape; 2],
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let deflection = ogeom_mesh::Deflection {
        chord: reach * 1e-3,
        ..ogeom_mesh::Deflection::default()
    };
    let mesh = ogeom_mesh::triangulate_face(model, strip, deflection, tol)?;
    // The strip's edges are as far off its neighbours' surfaces as their
    // tolerances say, and so is its border; the strip stands behind a
    // surface where it is no further in front of it than that.
    let mut slack = tol.confusion();
    for e in explore_unique(model, strip, ShapeType::Edge)? {
        if let Some(t) = model.tolerance_of(&e)? {
            slack = slack.max(t.get());
        }
    }
    for face in outer {
        let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
            return Ok(false);
        };
        let Some(stored) = model.geometry().surface(data.surface).cloned() else {
            return Ok(false);
        };
        let surface = {
            use ogeom_geom::Transformable as _;
            stored.transformed(&face.transform(model.datums())?, tol)?
        };
        let sign = if face.orientation() == ogeom_topo::Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
        for &p in &mesh.positions {
            let foot = ogeom_algo::project_on_surface(&surface, p, 32, tol)?;
            let (u, v) = foot.parameters;
            let normal = surface.normal_at(u, v, tol)?.vector() * sign;
            let at = surface.point_at(u, v, tol)?;
            if (p - at).dot(normal) > slack {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Whether an edge of `solid` other than the strip's own and the seams
/// stands within `reach` of the crossing `guide` anywhere, sampled along
/// each edge.
fn other_edge_within(
    model: &Model,
    solid: &Shape,
    strip: &Shape,
    guide: &ogeom_geom::Curve,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let own = explore_unique(model, strip, ShapeType::Edge)?;
    for e in explore_unique(model, solid, ShapeType::Edge)? {
        // The strip's own edges, and seams, which one face holds on both
        // sides and which bound nothing.
        if own.iter().any(|o| o.is_same(&e)) || faces_at(model, solid, &e)?.len() < 2 {
            continue;
        }
        let Some((curve, (a, b))) = ogeom_algo::placed_edge_curve(model, &e, tol)? else {
            continue;
        };
        for k in 0..=64 {
            let t = (b - a).mul_add(f64::from(k) / 64.0, a);
            let p: Point = curve.point_at(t, tol)?;
            if ogeom_algo::project_on_curve(guide, p, 64, tol)?.distance <= reach {
                return Ok(true);
            }
        }
    }
    Ok(false)
}
