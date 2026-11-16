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
//! solid to where, past the other, it has; the solid is cut by it.
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
use super::trimmed::{faces_of, outward};
use super::{Id, id};

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
/// between a plane square to its axis and a coaxial cylinder or cone, open,
/// whose band runs on past an end of its arc in the solid, each such rim's
/// ball sitting on its faces along the arc and its ends' faces meeting its
/// own at an angle; the other edges asked rounded alone by `round`. `None`
/// otherwise, or where a tool or a boolean will not build.
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
    let mut rims: Vec<(Rim, (f64, f64))> = Vec::new();
    let mut others: Vec<(Shape, f64)> = Vec::new();
    let mut runs_on = false;
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
        runs_on |= longer;
        rims.push((rim, span));
    }
    if !runs_on {
        return Ok(None);
    }
    let mut scratch = model.clone();
    // The other edges on marched seats: their wedges as each builds alone,
    // applied with the rims' tools, the rims' first.
    if !others.is_empty()
        && let Some(wedges) = super::marched_seats::wedges(&mut scratch, solid, &others, tol)?
    {
        let mut tools = Vec::with_capacity(rims.len() + wedges.len());
        for (rim, span) in &rims {
            let Some(made) = rim_tool(&mut scratch, &curved, rim, *span, tol)? else {
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
    for (rim, span) in &rims {
        let Some(made) = rim_tool(&mut scratch, &curved, rim, *span, tol)? else {
            return Ok(None);
        };
        tools.extend(made);
    }
    let Some(cut) = super::fallback::apply(&mut scratch, &base.shape, &tools, tol)? else {
        return Ok(None);
    };
    let mut history = base.history.then(&cut.history);
    for (rim, _) in &rims {
        if cut.history.generated(&rim.shape).is_empty() {
            return Ok(None);
        }
        for f in cut.history.generated(&rim.shape) {
            history.generate(&rim.shape, f.clone());
        }
        history.delete(&rim.shape);
    }
    for (edge, _) in asked {
        if history.generated(edge).is_empty() {
            return Ok(None);
        }
    }
    match super::guards::solid_fits(&scratch, &cut.shape, tol) {
        Ok(()) => {
            *model = scratch;
            Ok(Some(Built::new(cut.shape, history)))
        }
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(None),
    }
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
        let Some((made, band)) = tool(model, curved, rim, (t0, t1), standoff, tol)? else {
            return Ok(None);
        };
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
            let Some(n) = outward(curved, f, at) else {
                return false;
            };
            if n.cross(nt).magnitude() <= tol.angular() || n.cross(nw).magnitude() <= tol.angular()
            {
                return false;
            }
        }
    }
    true
}

/// The stations of the rim's circle its band runs between (ADR 0001,
/// amendment "corners where a round's rail misses"): from each end of its
/// arc on round its circle while the end face trims the band, its contact
/// on the cap standing in the solid and its contact on the wall's surface
/// carried on standing outside it, to the first station where that trim
/// reaches either contact, found by halving to rounding; and whether
/// either lies past its end of the arc. `None` where the band runs on half
/// a turn.
fn span(
    model: &Model,
    solid: &Shape,
    curved: &Curved,
    rim: &Rim,
    tol: Tolerances,
) -> OgeomResult<Option<((f64, f64), bool)>> {
    let step = 0.25 * rim.radius / rim.circle.radius();
    let steps = (core::f64::consts::PI / step).ceil();
    let trims = |t: f64| -> OgeomResult<Option<bool>> {
        let Some([_, _, on_top, on_wall]) = ball(curved, rim, t, tol)? else {
            return Ok(None);
        };
        let on = ogeom_algo::classify_in_solid_exact(model, solid, on_top, tol)?;
        let off = ogeom_algo::classify_in_solid_exact(model, solid, on_wall, tol)?;
        Ok(Some(
            on != ogeom_algo::Containment::Out && off == ogeom_algo::Containment::Out,
        ))
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
        let leaves = curved.faces.iter().enumerate().all(|(f, face)| {
            f == rim.top
                || f == rim.wall
                || !face.loops.iter().flatten().any(|(x, _)| {
                    let g = &curved.edges[x];
                    g.ends.0 == v || g.ends.1 == v
                })
                || outward(curved, f, at).is_some_and(|n| n.dot(ahead) > tol.angular())
        });
        if !leaves {
            continue;
        }
        let mut good = from;
        let mut bad = None;
        let mut i = 1.0;
        while i <= steps {
            let t = from + sign * step * i;
            match trims(t)? {
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
            match trims(mid)? {
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
