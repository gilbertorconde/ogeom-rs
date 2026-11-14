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
    strip: &Shape,
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
    // Over an open seat the round is runs capped in their end sections'
    // planes: the strip must stand between some run's two caps, and no
    // other edge stands within reach strictly between them. Round a closed
    // seat the round has no caps.
    let spans = if seat.closed() {
        Vec::new()
    } else {
        let mut spans = Vec::with_capacity(seat.runs().len());
        for run in seat.runs() {
            spans.push(run_caps(run, tol)?);
        }
        spans
    };
    let within = |p: Point, slack: f64| -> bool {
        spans.is_empty()
            || spans.iter().any(|caps| {
                caps.iter()
                    .all(|(on, normal)| (p - *on).dot(*normal) >= -slack)
            })
    };
    if !strip_behind(model, strip, outer, reach, &within, tol)? {
        return Ok(None);
    }
    if other_edge_within(model, solid, strip, &guide, reach, &within, tol)? {
        return Ok(None);
    }
    let mut built = crate::pair_marched::apply_to_solid(model, solid, seat, radius, true, tol)?;
    if !ogeom_algo::check(model, &built.shape, tol)?.is_valid() {
        return Ok(None);
    }
    built.history.delete(edge);
    crate::support::credit_new_faces(model, solid, edge, &mut built)?;
    Ok(Some(built))
}

/// Whether every point of the strip's mesh stands behind both outer faces'
/// surfaces, within its edges' tolerances: in the corner they make; and
/// `within` the round's span, the same slack allowed.
fn strip_behind(
    model: &Model,
    strip: &Shape,
    outer: [&Shape; 2],
    reach: f64,
    within: &dyn Fn(Point, f64) -> bool,
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
    // The caps stand where the march finds the ball leaving the faces,
    // which the round holds only to the bands' fit target: the strip is
    // between them within that, whatever its own edges hold.
    let between = slack.max(crate::marched::band_fit_target(tol));
    if mesh.positions.iter().any(|p| !within(*p, between)) {
        return Ok(false);
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
/// stands within `reach` of the crossing `guide` anywhere strictly
/// `within` the round's span, sampled along each edge.
fn other_edge_within(
    model: &Model,
    solid: &Shape,
    strip: &Shape,
    guide: &ogeom_geom::Curve,
    reach: f64,
    within: &dyn Fn(Point, f64) -> bool,
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
            // A point on or past a cap (a wall the round ends on) is the
            // cap's to meet, not the corner's.
            if !within(p, -tol.confusion()) {
                continue;
            }
            if ogeom_algo::project_on_curve(guide, p, 64, tol)?.distance <= reach {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// A capped run's two end planes, each as a point on it and its unit
/// normal into the run: the plane of the ball's arc at the run's first
/// station and at its last.
fn run_caps(
    run: &crate::march::MarchedBlend,
    tol: Tolerances,
) -> OgeomResult<[(Point, ogeom_math::Vector); 2]> {
    let n = run.len();
    if n < 2 {
        ogeom_core::ogeom_bail!(Construction, "a capped run of fewer than two stations");
    }
    let cap = |at: usize, toward: usize| -> OgeomResult<(Point, ogeom_math::Vector)> {
        let centre = run.spine[at];
        let normal = (run.touch_first[at] - centre).cross(run.touch_second[at] - centre);
        let m = normal.magnitude();
        if m <= tol.confusion() * tol.confusion() {
            ogeom_core::ogeom_bail!(Construction, "a run's end section has no plane");
        }
        let normal = normal / m;
        let inward = if (run.spine[toward] - centre).dot(normal) >= 0.0 {
            normal
        } else {
            -normal
        };
        Ok((centre, inward))
    };
    Ok([cap(0, n - 1)?, cap(n - 1, 0)?])
}
