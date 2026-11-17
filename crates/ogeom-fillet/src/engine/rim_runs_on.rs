//! A rim's round whose band runs on past an end of its arc (ADR 0001,
//! amendment "corners where a round's rail misses").
//!
//! At a vertex of three edges of one convexity where a rim's arc ends (a
//! boss's rim where a rib's wall leaves it, turned towards the material),
//! the round's rail on the cap may cross the edge beside it well back from
//! the vertex: the ball still finds material past the arc's end. The
//! round's surface then runs on past the vertex round its own circle until
//! the end face trims it, and the end face is trimmed back to it, as the
//! mainstream applications end a fillet there; it is not capped in the
//! meridian plane through the vertex. Its tool is the rim's section turned
//! about the axis from where, behind one end, the section has left the
//! solid to where, past the other, it has; the solid is cut by it. A band
//! that does not run on at an end stops there inside the material and ends
//! in the plane through the axis at the vertex where that vertex is a
//! corner convex in plan, or a seam where the rim's circle runs on
//! unasked, or where a face beside the cap or the wall carries it on; at an
//! inside corner such a free end is amendment "M5e"'s.
//!
//! Other edges asked with such a rim are rounded alone by the engine and
//! the rim's tool is cut from what that leaves: past where the bands meet
//! each runs out as it does alone, and the request takes the union of the
//! cuts, whatever the order.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{CircleCurve, Curve, Curve3d as _, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Axis, Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::curved::{Curved, read};
use super::defeatured::Rounder;
use super::fallback::{Sense, Tool};
use super::flush::{EndSurface, cut_far, cut_material_far, rejoined};
use super::trimmed::{faces_of, outward};
use super::{Id, id};

/// A point a band ends at and the surface it ends flush on there.
type EndCut = (Point, EndSurface);

/// A rim to build, with its span and its end cuts on planes and past walls.
type Planned = (Rim, (f64, f64), Vec<EndCut>, Vec<EndCut>);

/// A rim asked: its edge, its circle's arc, the cap and the wall, and its
/// radius.
struct Rim {
    shape: Shape,
    circle: ogeom_math::Circle,
    range: (f64, f64),
    top: usize,
    wall: usize,
    radius: f64,
    /// The arc's vertices, at its range's start and end.
    ends: (Id, Id),
}

/// How many stations along a rim its ball's contacts are read at, ends
/// apart, and how many points across the way to each.
const ALONG: u32 = 16;
const ACROSS: u32 = 8;

/// The rounds of `asked` on `solid` where some asked edge is a convex rim
/// between a plane square to its axis and a coaxial cylinder or cone,
/// open, its band run on past an end of its arc in the solid where it runs
/// on and else, at a corner convex in plan or a seam, ended in the plane
/// through the axis at the vertex, each
/// such rim's ball sitting on its faces along the arc and its ends' faces
/// meeting its own at an angle; the other edges asked rounded alone by
/// `round`. `None` otherwise, or where a tool or a boolean will not build.
///
/// # Errors
///
/// As reading the solid, `round` and the booleans, other than refusals,
/// which decline.
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
    let mut rims: Vec<Planned> = Vec::new();
    let mut others: Vec<(Shape, f64)> = Vec::new();
    for (edge, radius) in asked {
        let Some(rim) = rim_of(model, solid, &curved, id(edge), edge, *radius, tol)? else {
            others.push((edge.clone(), *radius));
            continue;
        };
        if !seated(model, solid, &curved, &rim, tol)? || !ends_apart(&curved, &rim, id(edge), tol) {
            return Ok(None);
        }
        let Some((span, longer)) = span(model, solid, &curved, &rim, tol)? else {
            return Ok(None);
        };
        // Where it does not run on, the band ends at its vertex in the
        // plane through the axis at a corner convex in plan or where the
        // cap or wall runs on; at an inside corner it ends flush on the
        // plane of the face it meets there, as amendment "M5e" ends a free
        // end (SMLib's mixed 3x1 corner: "the end face will be extended out
        // to the fillet surface").
        let Some(planes) = end_cuts(&curved, &rim, tol)? else {
            return Ok(None);
        };
        // Another asked crease leaving that inside corner is no end face:
        // the two rounds meet there, each carried on past the other
        // (amendment "M5e"), which this route does not build.
        if !longer && meets_asked(&curved, edge, asked, &planes, tol) {
            return Ok(None);
        }
        let Some(clips) = wall_clips(model, solid, &curved, &rim, tol)? else {
            return Ok(None);
        };
        rims.push((rim, span, planes, clips));
    }
    if rims.is_empty() {
        return Ok(None);
    }
    // A rim's vertex of three edges or more, every one asked, is a corner
    // of rounds, the corner patch's (one call rounds a cube's corner with
    // one ball); this route rounds the rim on past a vertex holding an edge
    // left sharp (amendment "corners where a round's rail misses"), or
    // through a vertex where two asked arcs continue one another.
    let asked_ids: Vec<Id> = asked.iter().map(|(e, _)| id(e)).collect();
    for (rim, ..) in &rims {
        for v in [rim.ends.0, rim.ends.1] {
            let at: Vec<Id> = curved
                .edges
                .iter()
                .filter(|(_, g)| g.ends.0 == v || g.ends.1 == v)
                .map(|(x, _)| *x)
                .collect();
            if at.len() >= 3 && at.iter().all(|x| asked_ids.contains(x)) {
                return Ok(None);
            }
        }
    }
    let mut scratch = model.clone();
    // The other edges on marched seats: their wedges as each builds alone,
    // applied with the rims' tools, the rims' first.
    if !others.is_empty()
        && let Some(wedges) = super::marched_seats::wedges(&mut scratch, solid, &others, tol)?
    {
        let mut tools = Vec::with_capacity(rims.len() + wedges.len());
        for (rim, span, planes, clips) in &rims {
            let Some(made) = rim_tool(
                &mut scratch,
                &curved,
                rim,
                *span,
                (planes, clips),
                solid,
                tol,
            )?
            else {
                return Ok(None);
            };
            tools.extend(made);
        }
        for wedge in &wedges {
            let Some(own) = wedge.edges.first() else {
                return Ok(None);
            };
            let faces = explore_unique(&scratch, &wedge.solid, ShapeType::Face)?
                .into_iter()
                .map(|face| {
                    let from = wedge
                        .credit
                        .iter()
                        .find(|(f, _)| f.is_same(&face))
                        .map_or_else(|| own.clone(), |(_, e)| e.clone());
                    (face, from)
                })
                .collect();
            tools.push(Tool {
                solid: wedge.solid.clone(),
                sense: if wedge.additive {
                    Sense::Fuse
                } else {
                    Sense::Cut
                },
                faces,
            });
        }
        let edges: Vec<Shape> = asked.iter().map(|(e, _)| e.clone()).collect();
        let Some(built) = super::fallback::applied(&mut scratch, solid, &edges, &tools, tol)?
        else {
            return Ok(None);
        };
        return match super::guards::solid_fits(&scratch, &built.shape, tol) {
            Ok(()) => {
                *model = scratch;
                Ok(Some(built))
            }
            Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
            Err(_) => Ok(None),
        };
    }
    let mut scratch = model.clone();
    // The other edges, rounded alone.
    let base = if others.is_empty() {
        Built::from_nothing(solid.clone())
    } else {
        let Some(built) = round(&mut scratch, solid, &others, tol)? else {
            return Ok(None);
        };
        built
    };
    let mut tools = Vec::with_capacity(rims.len());
    for (rim, span, planes, clips) in &rims {
        let Some(made) = rim_tool(
            &mut scratch,
            &curved,
            rim,
            *span,
            (planes, clips),
            solid,
            tol,
        )?
        else {
            return Ok(None);
        };
        tools.extend(made);
    }
    let Some(cut) = super::fallback::apply(&mut scratch, &base.shape, &tools, tol)? else {
        return Ok(None);
    };
    // A band ending flush on a plane leaves a face on that plane beside
    // the face it ended on: one face again, so a crease running on along
    // both stays one edge to its real end.
    let kept: Vec<Shape> = explore_unique(&scratch, &base.shape, ShapeType::Face)?
        .iter()
        .flat_map(|f| cut.history.trace(f).to_vec())
        .collect();
    let around: Vec<Shape> = explore_unique(&scratch, &cut.shape, ShapeType::Face)?
        .into_iter()
        .filter(|f| {
            !kept.iter().any(|k| k.is_same(f))
                && rims
                    .iter()
                    .flat_map(|(_, _, planes, _)| planes)
                    .any(|(_, end)| on_end_plane(&scratch, f, end, tol))
        })
        .collect();
    let unified = rejoined(
        &mut scratch,
        Built::from_nothing(cut.shape.clone()),
        &around,
        tol,
    )?;
    let mut history = base.history.then(&cut.history).then(&unified.history);
    for (rim, _, _, _) in &rims {
        let made: Vec<Shape> = cut
            .history
            .generated(&rim.shape)
            .iter()
            .flat_map(|f| unified.history.trace(f).to_vec())
            .collect();
        if made.is_empty() {
            return Ok(None);
        }
        for f in made {
            history.generate(&rim.shape, f);
        }
        history.delete(&rim.shape);
    }
    for (edge, _) in asked {
        if history.generated(edge).is_empty() {
            return Ok(None);
        }
    }
    let fits = super::guards::solid_fits(&scratch, &unified.shape, tol);
    match fits {
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        // Where the rule's own figure touches itself along an edge (two
        // bands' wall contacts at one depth, the gap between them closing
        // to the edge), it is zero thickness: refused by name, as
        // SOLIDWORKS refuses it (ADR 0001, amendment "M5e", bullet "not
        // built"). Read whether or not the guard passes: a torus and a cone
        // back to back pass it as a cusp edge, and only the side the void
        // stands on tells a crack from a knife of material.
        _ if super::guards::touches_itself(&scratch, &unified.shape, tol)? => {
            ogeom_core::ogeom_bail!(
                Construction,
                "a rim's round would leave the solid touching itself along an edge \
                 (zero thickness): its band and an earlier round's meet their walls \
                 at one depth there; see docs/CAPABILITIES.md, fillet.edge-blends"
            )
        }
        Ok(()) => {
            *model = scratch;
            Ok(Some(Built::new(unified.shape, history)))
        }
        Err(_) => Ok(None),
    }
}

/// Whether `face` lies on the plane `end` a band ends flush on.
fn on_end_plane(model: &Model, face: &Shape, end: &EndSurface, tol: Tolerances) -> bool {
    let EndSurface::Plane { at, keep } = end else {
        return false;
    };
    model
        .node(face)
        .and_then(|n| n.data().as_face())
        .and_then(|d| model.geometry().surface(d.surface))
        .is_some_and(|s| match s {
            SurfaceGeometry::Plane(p) => {
                let plane = p.plane();
                plane.normal().vector().cross(*keep).magnitude() <= tol.angular()
                    && keep.dot(plane.origin() - *at).abs() <= tol.confusion()
            }
            _ => false,
        })
}

/// Whether another of the `asked` creases leaves `edge` at a vertex where
/// its band ends flush, at one of `ends`.
fn meets_asked(
    curved: &Curved,
    edge: &Shape,
    asked: &[(Shape, f64)],
    ends: &[EndCut],
    tol: Tolerances,
) -> bool {
    asked.iter().any(|(other, _)| {
        id(other) != id(edge)
            && curved.edges.get(&id(other)).is_some_and(|g| {
                [g.ends.0, g.ends.1].iter().any(|v| {
                    curved.points.get(v).is_some_and(|p| {
                        ends.iter()
                            .any(|(at, _)| p.distance(*at) <= tol.confusion())
                    })
                })
            })
    })
}

/// The rim's band over `span` as tools, its band's faces named for the
/// rim: over its own arc the section standing off the wall (no tool face
/// on the solid's own wall), past each end the section bounded by the
/// wall's surface carried on, which is no face of the solid there and
/// takes nothing past it ([`span`] keeps the wall's contact outside the
/// solid there). `None` where one will not build.
fn rim_tool(
    model: &mut Model,
    curved: &Curved,
    rim: &Rim,
    span: (f64, f64),
    (planes, clips): (&[EndCut], &[EndCut]),
    solid: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Tool>>> {
    let pieces = [
        ((rim.range.0, rim.range.1), true),
        ((span.0, rim.range.0), false),
        ((rim.range.1, span.1), false),
    ];
    let mut out = Vec::with_capacity(3);
    // A piece shorter than a twentieth of the run's step is the reading's
    // own tolerance at a face the band ends square on.
    let least = 0.0125 * rim.radius / rim.circle.radius();
    for ((t0, t1), standoff) in pieces {
        if t1 - t0 <= least.max(tol.angular()) {
            continue;
        }
        let Some((mut made, mut band)) = tool(model, curved, rim, (t0, t1), standoff, tol)? else {
            return Ok(None);
        };
        // Cut flush on each inside corner's end plane, the crease's side
        // kept, the band's faces traced through the cut.
        for (at, end) in planes {
            let reach = 3.0 * rim.radius;
            let Some((kept, traced)) = cut_far(model, &made, &band, end, *at, reach, None, tol)?
            else {
                return Ok(None);
            };
            made = kept;
            band = traced;
        }
        // Where the solid holds material past the wall near an end, the
        // band's standoff past the wall takes none of it.
        for (at, wall) in clips {
            let reach = 3.0 * rim.radius;
            let Some((kept, traced)) =
                cut_material_far(model, &made, &band, wall, *at, reach, solid, tol)?
            else {
                return Ok(None);
            };
            made = kept;
            band = traced;
        }
        out.push(Tool {
            solid: made,
            sense: Sense::Cut,
            faces: band.into_iter().map(|f| (f, rim.shape.clone())).collect(),
        });
    }
    Ok(Some(out))
}

/// `edge` as a rim: an open arc of a circle between a planar cap square to
/// its axis and a coaxial cylinder or cone, convex; `None` otherwise.
fn rim_of(
    model: &Model,
    solid: &Shape,
    curved: &Curved,
    edge: Id,
    shape: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Rim>> {
    let Some(c) = curved.edges.get(&edge) else {
        return Ok(None);
    };
    let Curve::Circle(circle) = &c.curve else {
        return Ok(None);
    };
    if c.ends.0 == c.ends.1 {
        return Ok(None);
    }
    let circle = circle.circle();
    let axis = circle.frame().z().vector();
    let [fa, fb] = faces_of(curved, edge)[..] else {
        return Ok(None);
    };
    let (top, wall) = match (&curved.faces[fa].surface, &curved.faces[fb].surface) {
        (SurfaceGeometry::Plane(_), SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_)) => {
            (fa, fb)
        }
        (SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_), SurfaceGeometry::Plane(_)) => {
            (fb, fa)
        }
        _ => return Ok(None),
    };
    let mid = c.curve.point_at(0.5 * (c.range.0 + c.range.1), tol)?;
    let (Some(nt), Some(nw)) = (outward(curved, top, mid), outward(curved, wall, mid)) else {
        return Ok(None);
    };
    if nt.cross(axis).magnitude() > tol.angular() {
        return Ok(None);
    }
    // Convex: just in from the rim, behind both faces, the material.
    let inside = mid - (nt + nw) * (0.1 * radius);
    if ogeom_algo::classify_in_solid_exact(model, solid, inside, tol)?
        != ogeom_algo::Containment::In
    {
        return Ok(None);
    }
    Ok(Some(Rim {
        shape: shape.clone(),
        circle,
        range: c.range,
        top,
        wall,
        radius,
        ends: c.ends,
    }))
}

/// The rim's ball at parameter `t` of its circle, the faces' surfaces
/// carried on past the arc: the rim's point, the centre, and the contacts
/// on the cap and on the wall.
fn ball(curved: &Curved, rim: &Rim, t: f64, tol: Tolerances) -> OgeomResult<Option<[Point; 4]>> {
    let p = CircleCurve::new(rim.circle).point_at(t, tol)?;
    let (Some(nt), Some(nw)) = (outward(curved, rim.top, p), outward(curved, rim.wall, p)) else {
        return Ok(None);
    };
    let g = nt.dot(nw);
    if (1.0 + g).abs() <= tol.angular() {
        return Ok(None);
    }
    let centre = p - (nt + nw) * (rim.radius / (1.0 + g));
    Ok(Some([
        p,
        centre,
        centre + nt * rim.radius,
        centre + nw * rim.radius,
    ]))
}

/// Whether the rim's ball sits on both its faces along its arc: each
/// contact, and the way to it from the rim across its face, on the
/// solid's boundary.
fn seated(
    model: &Model,
    solid: &Shape,
    curved: &Curved,
    rim: &Rim,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let (a, b) = rim.range;
    for i in 1..ALONG {
        let t = a + (b - a) * f64::from(i) / f64::from(ALONG);
        let Some([p, _, on_top, on_wall]) = ball(curved, rim, t, tol)? else {
            return Ok(false);
        };
        for contact in [on_top, on_wall] {
            for j in 1..=ACROSS {
                let q = p + (contact - p) * (f64::from(j) / f64::from(ACROSS));
                if ogeom_algo::classify_in_solid_exact(model, solid, q, tol)?
                    != ogeom_algo::Containment::On
                {
                    return Ok(false);
                }
            }
        }
    }
    Ok(true)
}

/// Whether each end of the rim's arc stands on the solid's own faces, each
/// meeting the rim's at an angle: one running on tangent from them (a band
/// an earlier round left) is that round's to meet.
fn ends_apart(curved: &Curved, rim: &Rim, edge: Id, tol: Tolerances) -> bool {
    let c = &curved.edges[&edge];
    for v in [c.ends.0, c.ends.1] {
        let at = curved.points[&v];
        let (Some(nt), Some(nw)) = (outward(curved, rim.top, at), outward(curved, rim.wall, at))
        else {
            return false;
        };
        for (f, face) in curved.faces.iter().enumerate() {
            if f == rim.top
                || f == rim.wall
                || !face.loops.iter().flatten().any(|(x, _)| {
                    let g = &curved.edges[x];
                    g.ends.0 == v || g.ends.1 == v
                })
            {
                continue;
            }
            // The cap or the wall carried on past a boolean's split (a
            // neck fused on flush with the drum's top) is that face itself.
            if continues(curved, rim, f, tol) {
                continue;
            }
            let Some(n) = outward(curved, f, at) else {
                return false;
            };
            // A curved face tangent to the cap or wall there (an earlier
            // round's band) is the end face the band ends flush on.
            let curved_face = !matches!(face.surface, SurfaceGeometry::Plane(_));
            if !curved_face
                && (n.cross(nt).magnitude() <= tol.angular()
                    || n.cross(nw).magnitude() <= tol.angular())
            {
                return false;
            }
        }
    }
    true
}

/// The surfaces the rim's band ends flush on, one per end at an inside
/// corner, each with its vertex: the surface of the one face at the vertex,
/// past the cap and the wall and the faces carrying them on, that the arc
/// carried on runs into (the end face; an earlier round's band among them).
/// An end where the arc carried on leaves
/// through every such face (a corner convex in plan, or a seam the circle
/// runs on through) has none. `None` where an end meets more than one such
/// face, one [`EndSurface`] does not read, or the arc runs tangent into an
/// earlier round's band.
fn end_cuts(curved: &Curved, rim: &Rim, tol: Tolerances) -> OgeomResult<Option<Vec<EndCut>>> {
    let mut cuts = Vec::new();
    for (k, sign) in [(0, -1.0), (1, 1.0)] {
        let v = if k == 0 { rim.ends.0 } else { rim.ends.1 };
        let from = if k == 0 { rim.range.0 } else { rim.range.1 };
        let at = curved.points[&v];
        let arc = CircleCurve::new(rim.circle);
        let ahead = arc.d1_at(from, tol)? * sign;
        let ahead = ahead / ahead.magnitude();
        let mut blocking = Vec::new();
        for (f, face) in curved.faces.iter().enumerate() {
            if f == rim.top
                || f == rim.wall
                || !face.loops.iter().flatten().any(|(x, _)| {
                    let g = &curved.edges[x];
                    g.ends.0 == v || g.ends.1 == v
                })
                || continues(curved, rim, f, tol)
            {
                continue;
            }
            let Some(n) = outward(curved, f, at) else {
                return Ok(None);
            };
            // A curved face the arc carried on runs tangent to (an earlier
            // round's band beside the cap) is no end face: the band runs on
            // into it and is trimmed by it (amendments "a later round
            // reaching an earlier band at its end" and "a later round
            // reaching a band seated on a curved face").
            if !matches!(curved.faces[f].surface, SurfaceGeometry::Plane(_))
                && n.dot(ahead).abs() <= 1e-6
            {
                continue;
            }
            if n.dot(ahead) > tol.angular() {
                continue;
            }
            blocking.push((f, n));
        }
        // The crease's side: where the arc stands just back from its end.
        let back = from - sign * 1e-3 * (rim.range.1 - rim.range.0);
        let q = arc.point_at(back, tol)?;
        match blocking.as_slice() {
            [] => {}
            [(f, n)] => {
                let Some(end) = EndSurface::of(&curved.faces[*f].surface, at, *n, q) else {
                    return Ok(None);
                };
                cuts.push((at, end));
            }
            _ => return Ok(None),
        }
    }
    Ok(Some(cuts))
}

/// The rim's wall surface, its outside the far side, at each end of its
/// arc where the solid holds material just past the wall within the arc
/// (a block fused against the drum, at an inside corner), with that end's
/// vertex. `None` where the wall's surface is not one [`EndSurface`] reads.
fn wall_clips(
    model: &Model,
    solid: &Shape,
    curved: &Curved,
    rim: &Rim,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<EndCut>>> {
    let mut clips = Vec::new();
    let arc = CircleCurve::new(rim.circle);
    for (k, sign) in [(0, 1.0), (1, -1.0)] {
        let from = if k == 0 { rim.range.0 } else { rim.range.1 };
        let v = if k == 0 { rim.ends.0 } else { rim.ends.1 };
        let at = curved.points[&v];
        // Stations into the arc from the end, the nearest a few
        // thousandths of the radius on, each the wall's contact carried a
        // hundredth of the radius out past the wall: the material a block
        // leaning against the drum holds there is a thin wedge at the end.
        let mut held = false;
        for step in [0.002, 0.01, 0.05, 0.25] {
            let t = from + sign * step * rim.radius / rim.circle.radius();
            if (t - rim.range.0) * (rim.range.1 - t) < 0.0 {
                break;
            }
            let Some([_, _, _, on_wall]) = ball(curved, rim, t, tol)? else {
                return Ok(None);
            };
            let Some(nw) = outward(curved, rim.wall, arc.point_at(t, tol)?) else {
                return Ok(None);
            };
            let probe = on_wall + nw * (0.01 * rim.radius);
            if ogeom_algo::classify_in_solid_exact(model, solid, probe, tol)?
                == ogeom_algo::Containment::In
            {
                held = true;
                break;
            }
        }
        if !held {
            continue;
        }
        let Some(nw) = outward(curved, rim.wall, at) else {
            return Ok(None);
        };
        // The drum's side: a point just inside the wall at the end.
        let Some(wall) = EndSurface::of(
            &curved.faces[rim.wall].surface,
            at,
            nw,
            at - nw * (0.1 * rim.radius),
        ) else {
            return Ok(None);
        };
        clips.push((at, wall));
    }
    Ok(Some(clips))
}

/// Whether face `f` carries on the rim's cap or wall: a plane the cap's
/// own, or the wall's own cylinder, faced alike.
fn continues(curved: &Curved, rim: &Rim, f: usize, tol: Tolerances) -> bool {
    let same = |a: &SurfaceGeometry, b: &SurfaceGeometry| match (a, b) {
        (SurfaceGeometry::Plane(p), SurfaceGeometry::Plane(q)) => {
            let (p, q) = (p.plane(), q.plane());
            p.normal().vector().cross(q.normal().vector()).magnitude() <= tol.angular()
                && p.distance_to(q.origin()) <= tol.confusion()
        }
        (SurfaceGeometry::Cone(p), SurfaceGeometry::Cone(q)) => {
            let (p, q) = (p.cone(), q.cone());
            (p.half_angle() - q.half_angle()).abs() <= tol.angular()
                && p.apex().distance(q.apex()) <= tol.confusion()
                && p.axis()
                    .direction
                    .vector()
                    .cross(q.axis().direction.vector())
                    .magnitude()
                    <= tol.angular()
        }
        (SurfaceGeometry::Cylinder(p), SurfaceGeometry::Cylinder(q)) => {
            let (p, q) = (p.cylinder(), q.cylinder());
            (p.radius() - q.radius()).abs() <= tol.confusion()
                && p.axis()
                    .direction
                    .vector()
                    .cross(q.axis().direction.vector())
                    .magnitude()
                    <= tol.angular()
                && p.axis().distance_to(q.axis().location) <= tol.confusion()
        }
        _ => false,
    };
    let faced = |g: usize| curved.faces[g].shape.orientation();
    [rim.top, rim.wall].iter().any(|&g| {
        g != f && faced(g) == faced(f) && same(&curved.faces[g].surface, &curved.faces[f].surface)
    })
}

/// The stations of the rim's circle its band runs between (ADR 0001,
/// amendment "corners where a round's rail misses"): from each end of its
/// arc on round its circle while the end face trims the band, its contact
/// on the cap standing in the solid and its contact on the wall's surface
/// carried on standing outside it, to the first station where that trim
/// reaches either contact, found by halving to rounding, and never past
/// where the cap's contact first meets the cap's boundary, crossing or
/// touching it (amendment "a band's run past its end"); and whether either
/// lies past its end of the arc. `None` where the band runs on half a
/// turn.
fn span(
    model: &Model,
    solid: &Shape,
    curved: &Curved,
    rim: &Rim,
    tol: Tolerances,
) -> OgeomResult<Option<((f64, f64), bool)>> {
    let step = 0.25 * rim.radius / rim.circle.radius();
    let steps = (core::f64::consts::PI / step).ceil();
    let trims = |t: f64, at_band: bool| -> OgeomResult<Option<bool>> {
        let Some([_, _, on_top, on_wall]) = ball(curved, rim, t, tol)? else {
            return Ok(None);
        };
        let on = ogeom_algo::classify_in_solid_exact(model, solid, on_top, tol)?;
        let off = ogeom_algo::classify_in_solid_exact(model, solid, on_wall, tol)?;
        if !at_band {
            return Ok(Some(
                on != ogeom_algo::Containment::Out && off == ogeom_algo::Containment::Out,
            ));
        }
        // Past an end where an earlier round's band took the edge, not its
        // faces (amendment "a later round reaching a band seated on a curved
        // face"), the ball runs on seated on its own faces, its contact on
        // the cap in air only where that band cut the cap back, until a
        // contact is in the material.
        let on_cap = on != ogeom_algo::Containment::Out
            || outward(curved, rim.top, on_top).is_some_and(|n| {
                super::ruling_out::cut_back(
                    model,
                    solid,
                    curved,
                    (on_top, n),
                    4.0 * rim.radius,
                    tol,
                )
                .unwrap_or(false)
            });
        Ok(Some(on_cap && off != ogeom_algo::Containment::In))
    };
    // Whether the contact on the cap at `t` stands in the solid, and
    // whether the section there holds any of it, read at points across the
    // corner it takes, outside the ball.
    let cap_held = |t: f64| -> OgeomResult<bool> {
        let Some([_, _, on_top, _]) = ball(curved, rim, t, tol)? else {
            return Ok(true);
        };
        Ok(
            ogeom_algo::classify_in_solid_exact(model, solid, on_top, tol)?
                != ogeom_algo::Containment::Out,
        )
    };
    let holds = |t: f64| -> OgeomResult<bool> {
        let Some([p, centre, on_top, on_wall]) = ball(curved, rim, t, tol)? else {
            return Ok(true);
        };
        for i in 0..=8 {
            for j in 0..=(8 - i) {
                let (a, b) = (f64::from(i) / 8.0, f64::from(j) / 8.0);
                let q = p + (on_top - p) * a + (on_wall - p) * b;
                if q.distance(centre) < rim.radius {
                    continue;
                }
                if ogeom_algo::classify_in_solid_exact(model, solid, q, tol)?
                    == ogeom_algo::Containment::In
                {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    };
    let mut ends = [rim.range.0, rim.range.1];
    let mut longer = false;
    for (k, sign) in [(0, -1.0), (1, 1.0)] {
        let from = ends[k];
        // Only at a corner convex in plan, the arc carried on leaving
        // through every face at its vertex but its own (all its edges of
        // one convexity): at a re-entrant corner it runs into the material.
        let v = if k == 0 { rim.ends.0 } else { rim.ends.1 };
        let at = curved.points[&v];
        let ahead = CircleCurve::new(rim.circle).d1_at(from, tol)? * sign;
        let at_band = curved.faces.iter().enumerate().any(|(f, face)| {
            f != rim.top
                && f != rim.wall
                && !matches!(face.surface, SurfaceGeometry::Plane(_))
                && face.loops.iter().flatten().any(|(x, _)| {
                    let g = &curved.edges[x];
                    g.ends.0 == v || g.ends.1 == v
                })
                && outward(curved, f, at)
                    .is_some_and(|n| n.dot(ahead / ahead.magnitude()).abs() <= 1e-6)
        });
        let leaves = curved.faces.iter().enumerate().all(|(f, face)| {
            f == rim.top
                || f == rim.wall
                || !face.loops.iter().flatten().any(|(x, _)| {
                    let g = &curved.edges[x];
                    g.ends.0 == v || g.ends.1 == v
                })
                || outward(curved, f, at).is_some_and(|n| n.dot(ahead) > tol.angular())
                // An earlier round's band the arc runs tangent into: the
                // band runs on past, trimmed by it.
                || (!matches!(face.surface, SurfaceGeometry::Plane(_))
                    && outward(curved, f, at)
                        .is_some_and(|n| n.dot(ahead / ahead.magnitude()).abs() <= 1e-6))
        });
        if !leaves {
            continue;
        }
        // The first point past the end where the contact on the cap meets
        // the cap's boundary, crossing it or touching it: past there the
        // section holds no material connected to the arc's (ADR 0001,
        // amendment "a band's run past its end").
        let limit = cap_boundary(curved, rim, from, sign, tol)?;
        let within = |t: f64| limit.is_none_or(|l| sign * (l - t) > 0.0);
        let mut good = from;
        let mut bad = None;
        let mut i = 1.0;
        while i <= steps {
            let t = from + sign * step * i;
            if !within(t) {
                bad = limit;
                break;
            }
            match trims(t, at_band)? {
                None => return Ok(None),
                Some(true) => good = t,
                Some(false) => {
                    bad = Some(t);
                    break;
                }
            }
            i += 1.0;
        }
        let Some(mut bad) = bad else {
            return Ok(None);
        };
        for _ in 0..60 {
            let mid = f64::midpoint(good, bad);
            if !within(mid) {
                bad = mid;
                continue;
            }
            match trims(mid, at_band)? {
                None => return Ok(None),
                Some(true) => good = mid,
                Some(false) => bad = mid,
            }
        }
        // Past an end it stops square on, the reading's own tolerance at the
        // face it ends on: no running on.
        longer |= (good - from).abs() > 0.05 * step;
        // Where the band stops because its contact on the cap has left the
        // solid, it holds nothing past there: the tool runs half a step on,
        // its end clear of where that contact crosses the cap's edge. Where
        // it stops because its contact on the wall's surface enters the
        // solid, it ends there in the plane through the axis, set a few
        // dozen confusions back along the wall so that plane's edge on the
        // wall's surface does not lie on the end face it reaches.
        let past = good + sign * 0.5 * step;
        ends[k] = if !cap_held(past)? && !holds(past)? {
            past
        } else {
            good - sign * (32.0 * tol.confusion() / rim.circle.radius())
        };
    }
    Ok(Some(((ends[0], ends[1]), longer)))
}

/// Where the rim's contact on its cap, run on from `from` the way `sign`
/// takes it round the circle, first meets an edge of the cap other than
/// the rim's own, crossing it or touching it: a parameter on the rim's
/// circle, or `None` where it meets none within a turn.
///
/// # Errors
///
/// As intersecting the curves.
fn cap_boundary(
    curved: &Curved,
    rim: &Rim,
    from: f64,
    sign: f64,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let Some([_, _, on_top, _]) = ball(curved, rim, from, tol)? else {
        return Ok(None);
    };
    let frame = rim.circle.frame();
    let axis = frame.z().vector();
    let off = on_top - frame.origin();
    let radius = (off - axis * axis.dot(off)).magnitude();
    let height = axis.dot(off);
    let at = Frame::new(frame.origin() + axis * height, frame.z(), frame.x(), tol)?;
    let contact_circle = ogeom_math::Circle::new(at, radius, tol)?;
    let contact: Curve = CircleCurve::new(contact_circle).into();
    let period = core::f64::consts::TAU;
    let mut nearest: Option<f64> = None;
    for (x, _) in curved.faces[rim.top].loops.iter().flatten() {
        if *x == id(&rim.shape) {
            continue;
        }
        // Where the cap runs on tangent into an earlier round's band, the
        // band is no boundary: the cap's contact runs on past it to where
        // it meets the band's other seat, a plane, as it would have met
        // that face's edge with the cap had the band not been there.
        if let Some(seats) = band_seats(curved, rim, *x, tol)? {
            for plane in seats {
                for t in circle_meets_plane(&contact_circle, plane, tol) {
                    let ahead = (sign * (t - from)).rem_euclid(core::f64::consts::TAU);
                    if ahead > tol.parametric() && nearest.is_none_or(|n| ahead < n) {
                        nearest = Some(ahead);
                    }
                }
            }
            continue;
        }
        let e = &curved.edges[x];
        let (dlo, dhi) = e.curve.domain();
        let boundary: Curve = if (e.range.0 - dlo).abs() <= tol.parametric()
            && (e.range.1 - dhi).abs() <= tol.parametric()
        {
            e.curve.clone()
        } else {
            ogeom_geom::TrimmedCurve::new(e.curve.clone(), e.range.0, e.range.1, tol)?.into()
        };
        let found = ogeom_intersect::curves::intersect_curves(
            &contact,
            &boundary,
            ogeom_intersect::curves::CurveCurveOptions {
                gap: tol.confusion(),
                ..ogeom_intersect::curves::CurveCurveOptions::default()
            },
            tol,
        )?;
        for c in found.crossings {
            let ahead = (sign * (c.on_a - from)).rem_euclid(period);
            if ahead > tol.parametric() && nearest.is_none_or(|n| ahead < n) {
                nearest = Some(ahead);
            }
        }
    }
    Ok(nearest.map(|d| from + sign * d))
}

/// Where the cap meets an earlier round's band tangentially along edge
/// `x`: that band's other seats that are planes, each by a point and its
/// unit normal; `None` where `x` is no such edge.
fn band_seats(
    curved: &Curved,
    rim: &Rim,
    x: Id,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<(Point, ogeom_math::Vector)>>> {
    let e = &curved.edges[&x];
    let mid = e.curve.point_at(0.5 * (e.range.0 + e.range.1), tol)?;
    let Some(nt) = outward(curved, rim.top, mid) else {
        return Ok(None);
    };
    let tangent =
        |f: usize, g: usize, at: Point| match (outward(curved, f, at), outward(curved, g, at)) {
            (Some(a), Some(b)) => a.cross(b).magnitude() <= 1e-6 && a.dot(b) > 0.0,
            _ => false,
        };
    let Some(band) = faces_of(curved, x).into_iter().find(|g| {
        *g != rim.top
            && !matches!(curved.faces[*g].surface, SurfaceGeometry::Plane(_))
            && outward(curved, *g, mid)
                .is_some_and(|n| n.cross(nt).magnitude() <= 1e-6 && n.dot(nt) > 0.0)
    }) else {
        return Ok(None);
    };
    let mut seats = Vec::new();
    for (y, _) in curved.faces[band].loops.iter().flatten() {
        let edge = &curved.edges[y];
        let at = edge
            .curve
            .point_at(0.5 * (edge.range.0 + edge.range.1), tol)?;
        for g in faces_of(curved, *y) {
            if g == band || g == rim.top || !tangent(band, g, at) {
                continue;
            }
            if let SurfaceGeometry::Plane(p) = &curved.faces[g].surface {
                let plane = p.plane();
                let seat = (plane.origin(), plane.normal().vector());
                if !seats.iter().any(|(o, n): &(Point, ogeom_math::Vector)| {
                    n.cross(seat.1).magnitude() <= 1e-9
                        && (seat.0 - *o).dot(*n).abs() <= tol.confusion()
                }) {
                    seats.push(seat);
                }
            }
        }
    }
    Ok(Some(seats))
}

/// The parameters at which `circle` meets the plane through `plane.0` with
/// normal `plane.1`, a touch counted once.
fn circle_meets_plane(
    circle: &ogeom_math::Circle,
    (on, n): (Point, ogeom_math::Vector),
    tol: Tolerances,
) -> Vec<f64> {
    let frame = circle.frame();
    let (e1, e2) = (frame.x().vector(), frame.y().vector());
    let r = circle.radius();
    let (a, b, d) = (r * n.dot(e1), r * n.dot(e2), n.dot(frame.origin() - on));
    let m = a.hypot(b);
    if m <= tol.confusion() {
        return Vec::new();
    }
    let c = -d / m;
    if c.abs() > 1.0 + tol.confusion() / r {
        return Vec::new();
    }
    let base = b.atan2(a);
    let spread = c.clamp(-1.0, 1.0).acos();
    if spread <= 1e-9 {
        vec![base]
    } else {
        vec![base + spread, base - spread]
    }
}

/// The rim's band over `span` of its circle: its section at the span's
/// start, from the rim's point by its two contacts and the ball's arc,
/// turned about the axis to the span's end; with the band's faces in it.
/// The section's other sides stand off the faces' surfaces, a quarter
/// radius out, so no tool face lies on the solid's own.
fn tool(
    model: &mut Model,
    curved: &Curved,
    rim: &Rim,
    (t0, t1): (f64, f64),
    standoff: bool,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
    let built = (|| -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
        let Some([p, centre, on_top, on_wall]) = ball(curved, rim, t0, tol)? else {
            return Ok(None);
        };
        let r = rim.radius;
        let (x, y) = (on_top - centre, on_wall - centre);
        let normal = x.cross(y);
        if normal.magnitude() <= tol.confusion() * r {
            return Ok(None);
        }
        let Some(nw) = outward(curved, rim.wall, p) else {
            return Ok(None);
        };
        let Some(nt) = outward(curved, rim.top, p) else {
            return Ok(None);
        };
        let off = r * 0.25;
        let past = p + (p - on_top) / p.distance(on_top) * off;
        let out = on_wall + nw * off;
        // Past the arc the cap's side stands off the cap, above it (no tool
        // face on the solid's own cap there either), and the wall's side on
        // the wall's surface carried on, outside the solid.
        let (lift_p, lift_t) = (p + nt * off, on_top + nt * off);
        let [vt, vq, vo, vw, vp, vlp, vlt] = [on_top, past, out, on_wall, p, lift_p, lift_t]
            .map(|q| ogeom_algo::make_vertex(model, q).shape);
        let mut line = |from: Point, to: Point, a: &Shape, b: &Shape| -> OgeomResult<Shape> {
            Ok(ogeom_algo::make_edge_between(
                model,
                LineCurve::segment(from, to, tol)?.into(),
                (0.0, from.distance(to)),
                a,
                b,
                tol,
            )?
            .shape)
        };
        // Without the standoff, the rim's point closes the section, the wall's
        // surface and the cap's its two other sides.
        let (back, across, along_top) = if standoff {
            (
                line(on_wall, out, &vw, &vo)?,
                line(out, past, &vo, &vq)?,
                line(past, on_top, &vq, &vt)?,
            )
        } else {
            (
                line(on_wall, p, &vw, &vp)?,
                line(p, lift_p, &vp, &vlp)?,
                line(lift_p, lift_t, &vlp, &vlt)?,
            )
        };
        let down = if standoff {
            None
        } else {
            Some(line(lift_t, on_top, &vlt, &vt)?)
        };
        let frame = Frame::new(
            centre,
            Direction::new(normal, tol)?,
            Direction::new(x, tol)?,
            tol,
        )?;
        let arc_curve: Curve = CircleCurve::new(ogeom_math::Circle::new(frame, r, tol)?).into();
        let arc = ogeom_algo::make_edge_between(
            model,
            arc_curve,
            (0.0, x.angle(y, tol)?),
            &vt,
            &vw,
            tol,
        )?
        .shape;
        let mut sides = vec![arc.clone(), back, across, along_top];
        sides.extend(down);
        let wire = ogeom_algo::make_wire(model, &sides, tol)?.shape;
        let tangent = CircleCurve::new(rim.circle).d1_at(t0, tol)?;
        let plane = PlaneSurface::new(ogeom_math::Plane::through(p, Direction::new(tangent, tol)?));
        let face = ogeom_algo::make_face(model, plane.into(), &[wire], tol)?.shape;
        let turned = ogeom_algo::make_revolution(
            model,
            &face,
            Axis {
                location: rim.circle.centre(),
                direction: rim.circle.frame().z(),
            },
            t1 - t0,
            tol,
        )?;
        let band = turned.history.generated(&arc).to_vec();
        Ok((!band.is_empty()).then_some((turned.shape, band)))
    })();
    match built {
        Ok(found) => Ok(found),
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(None),
    }
}
