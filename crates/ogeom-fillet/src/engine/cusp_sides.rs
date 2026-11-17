//! Rounds either side of a cusp (ADR 0001, amendment "rounds either side
//! of a cusp"): where a flat side runs tangent into a drum's wall, the flat
//! top between them narrows to a cusp at the tangent point. A round of the
//! top's rim on the wall (the drum's foot or head) and a round of the top's
//! edge with the side meet there tangentially, the two edges continuing
//! each other through the cusp.
//!
//! Where no asked rim and asked straight edge leave the cusp on one side of
//! it (that pair turns through an S, the trimmed assembly's cusp), each
//! round is its own closed tool, applied by the one boolean:
//!
//! - a rim's round is its section (the corner between the cap and the wall
//!   less the ball's disc, or the ball's fill) turned about the drum's axis
//!   over its arc, ending in the meridian planes through the arc's ends;
//!   where it runs past the side at a cusp, the side's plane cuts it off,
//!   its fill ending flush with the side;
//! - a straight edge's round is its section swept along the edge, ending in
//!   the plane square to the edge through the cusp, as a band ends where a
//!   round of another kind continues it tangentially, and at its far end in
//!   its end face's plane.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{CircleCurve, Curve, Curve3d as _, LineCurve, PlaneSurface};
use ogeom_math::{Axis, Circle, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::fallback::{Sense, Tool};

/// A cusp of the solid: the tangent point, and the side's plane by a point
/// on it and its outward normal.
struct Cusp {
    at: Point,
    side: (Point, Vector),
    /// The top, the rim's cap.
    top: Shape,
}

/// An asked edge, read.
enum Asked {
    /// A rim between a planar cap and a coaxial cylinder.
    Rim(Box<crate::support::RevolvedSeat>, CircleCurve, (f64, f64)),
    /// A straight edge between two planes, by its two faces' outward
    /// normals.
    Straight([Shape; 2]),
}

/// A face's outward normal where it is a plane.
fn plane_of(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Option<(Point, Vector)>> {
    let planar = model
        .node(face)
        .and_then(|n| n.data().as_face())
        .and_then(|d| model.geometry().surface(d.surface))
        .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Plane(_)));
    if !planar {
        return Ok(None);
    }
    Ok(Some(ogeom_algo::face_normal(model, face, tol)?))
}

/// The faces of `solid` holding `edge`, each once.
fn faces_at(model: &Model, solid: &Shape, edge: &Shape) -> OgeomResult<Vec<Shape>> {
    let mut found: Vec<Shape> = Vec::new();
    for f in ogeom_topo::ancestors_of(model, solid, edge, ShapeType::Face)? {
        if !found.iter().any(|g| g.is_same(&f)) {
            found.push(f);
        }
    }
    Ok(found)
}

/// A vertex's point.
fn point_of(model: &Model, vertex: &Shape) -> OgeomResult<Point> {
    let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
        ogeom_core::ogeom_bail!(Construction, "a vertex holds no point");
    };
    Ok(vertex.transform(model.datums())?.apply(data.point))
}

/// `edge` read as a rim or a straight edge between planes; `None` for any
/// other.
fn read(model: &Model, solid: &Shape, edge: &Shape, tol: Tolerances) -> OgeomResult<Option<Asked>> {
    let (curve, range) = crate::support::edge_curve(model, edge, tol)?;
    match curve {
        Curve::Circle(c) => match crate::support::revolved_seat(model, solid, edge, &c, tol) {
            Ok(seat) if seat.slope == 0.0 => Ok(Some(Asked::Rim(Box::new(seat), c, range))),
            Ok(_) | Err(ogeom_core::OgeomError::Construction(_)) => Ok(None),
            Err(e) => Err(e),
        },
        Curve::Line(_) => {
            let faces = faces_at(model, solid, edge)?;
            let [a, b] = faces.as_slice() else {
                return Ok(None);
            };
            if plane_of(model, a, tol)?.is_none() || plane_of(model, b, tol)?.is_none() {
                return Ok(None);
            }
            Ok(Some(Asked::Straight([a.clone(), b.clone()])))
        }
        _ => Ok(None),
    }
}

/// The cusps on a rim's circle: a vertex of the cap where a straight edge
/// of the cap with a planar side leaves, the side tangent to the wall
/// there, the wall's outward normal the side's.
fn cusps_of(
    model: &Model,
    solid: &Shape,
    seat: &crate::support::RevolvedSeat,
    tol: Tolerances,
) -> OgeomResult<Vec<Cusp>> {
    let mut found = Vec::new();
    for e in explore_unique(model, &seat.cap_face, ShapeType::Edge)? {
        let (curve, _) = crate::support::edge_curve(model, &e, tol)?;
        if !matches!(curve, Curve::Line(_)) {
            continue;
        }
        let Some(side) = faces_at(model, solid, &e)?
            .into_iter()
            .find(|f| !f.is_same(&seat.cap_face))
        else {
            continue;
        };
        let Some((on, out)) = plane_of(model, &side, tol)? else {
            continue;
        };
        let Some((a, b)) = ogeom_algo::edge_vertices(model, &e)? else {
            continue;
        };
        for v in [&a, &b] {
            let p = point_of(model, v)?;
            let d = p - seat.centre;
            let radial = d - seat.up * d.dot(seat.up);
            if (radial.magnitude() - seat.radius).abs() > tol.confusion() * 10.0
                || d.dot(seat.up).abs() > tol.confusion() * 10.0
            {
                continue;
            }
            let wall_out = radial / radial.magnitude() * seat.sigma;
            let other = if v.is_same(&a) { &b } else { &a };
            let along = point_of(model, other)? - p;
            if wall_out.dot(out) >= 1.0 - 1e-9
                && pinched(model, &seat.cap_face, p, along, tol)?
                && !found
                    .iter()
                    .any(|c: &Cusp| c.at.distance(p) <= tol.confusion())
            {
                found.push(Cusp {
                    at: p,
                    side: (on, out),
                    top: seat.cap_face.clone(),
                });
            }
        }
    }
    Ok(found)
}

/// Whether `cap` narrows to nothing at `p`: an arc of its boundary ending
/// at `p` leaves it the way its straight edge does, `along`. Where the arc
/// leaves the other way, the two edges continue each other through `p`
/// and the cap runs on straight across it (a concave band's end face
/// joined to the face it ends on): a smooth turn of the boundary, no cusp.
/// No arc ending at `p` reads as pinched.
fn pinched(
    model: &Model,
    cap: &Shape,
    p: Point,
    along: Vector,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let near = tol.confusion() * 10.0;
    for e in explore_unique(model, cap, ShapeType::Edge)? {
        let (curve, (t0, t1)) = crate::support::edge_curve(model, &e, tol)?;
        if !matches!(curve, Curve::Circle(_)) {
            continue;
        }
        let step = (t1 - t0) * 1e-3;
        for (end, inward) in [(t0, t0 + step), (t1, t1 - step)] {
            if curve.point_at(end, tol)?.distance(p) <= near {
                let leaves = curve.point_at(inward, tol)? - p;
                return Ok(leaves.dot(along) > 0.0);
            }
        }
    }
    Ok(true)
}

/// The section of a round: the corner `c`, its contacts `a` and `b`, and
/// the arc between them about `o`, as a planar face in the plane `normal`
/// stands square to.
fn section(
    model: &mut Model,
    [c, a, b, o]: [Point; 4],
    normal: Vector,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let x = Direction::new(a - o, tol)?;
    let mut z = Direction::new(normal, tol)?;
    let y = z.vector().cross(x.vector());
    let turn = |y: Vector| (b - o).dot(y).atan2((b - o).dot(x.vector()));
    if turn(y) < 0.0 {
        z = Direction::new(-normal, tol)?;
    }
    let frame = Frame::new(o, z, x, tol)?;
    let y = z.vector().cross(x.vector());
    let arc = Curve::Circle(CircleCurve::new(Circle::new(frame, radius, tol)?));
    let [vc, va, vb] = [c, a, b].map(|p| ogeom_algo::make_vertex(model, p).shape);
    let line = |model: &mut Model, (p, vp): (Point, &Shape), (q, vq): (Point, &Shape)| {
        let curve = Curve::Line(LineCurve::segment(p, q, tol)?);
        let domain = curve.domain();
        Ok::<_, ogeom_core::OgeomError>(
            ogeom_algo::make_edge_between(model, curve, domain, vp, vq, tol)?.shape,
        )
    };
    let ca = line(model, (c, &vc), (a, &va))?;
    let ab = ogeom_algo::make_edge_between(model, arc, (0.0, turn(y)), &va, &vb, tol)?.shape;
    let bc = line(model, (b, &vb), (c, &vc))?;
    let wire = ogeom_algo::make_wire(model, &[ca, ab, bc], tol)?.shape;
    let plane = Plane::new(Frame::new(c, Direction::new(normal, tol)?, x, tol)?);
    Ok(ogeom_algo::make_face(model, PlaneSurface::new(plane).into(), &[wire], tol)?.shape)
}

/// `tool` less what lies past the plane through `on` on the side `out`
/// points to.
fn cut_past(
    model: &mut Model,
    tool: &Shape,
    (on, out): (Point, Vector),
    span: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let outward = Direction::new(out, tol)?;
    let across = {
        let x = if out.x.abs() < 0.9 {
            Vector::new(1.0, 0.0, 0.0)
        } else {
            Vector::new(0.0, 1.0, 0.0)
        };
        Direction::new(x - out * x.dot(out), tol)?
    };
    let frame = Frame::new(on, outward, across, tol)?;
    let corner = on - (frame.x().vector() + frame.y().vector()) * span;
    let beyond = ogeom_algo::make_box(
        model,
        Frame::new(corner, outward, across, tol)?,
        (span * 2.0, span * 2.0, span * 2.0),
        tol,
    )?
    .shape;
    let kept = ogeom_bool::cut(model, tool, &beyond, tol)?.shape;
    let solids = explore_unique(model, &kept, ShapeType::Solid)?;
    Ok(match solids.as_slice() {
        [one] => Some(one.clone()),
        _ => None,
    })
}

/// A rim's tool over its arc, cut off at each cusp's side it reaches.
fn rim_tool(
    model: &mut Model,
    (seat, circle, range): (&crate::support::RevolvedSeat, &CircleCurve, (f64, f64)),
    radius: f64,
    (cuts, whole): (&[(Point, Vector)], bool),
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let ball = seat.ball(radius);
    let up = seat.up;
    let start = circle.point_at(range.0, tol)?;
    let end = circle.point_at(range.1, tol)?;
    let radial_of = |p: Point| {
        let d = p - seat.centre;
        let r = d - up * d.dot(up);
        r / r.magnitude()
    };
    let x0 = radial_of(start);
    // The arc's turn about `up`, the way the circle runs.
    let runs_up = circle.circle().frame().z().vector().dot(up) > 0.0;
    let closed = whole || start.distance(end) <= tol.confusion();
    let sweep = if closed {
        core::f64::consts::TAU
    } else {
        let x1 = radial_of(end);
        let mut a = x0.cross(x1).dot(up).atan2(x0.dot(x1));
        if !runs_up {
            a = -a;
        }
        if a <= 0.0 {
            a += core::f64::consts::TAU;
        }
        a
    };
    // The section in the meridian through the arc's start; where the arc
    // runs against `up`, the turn starts from its end.
    let from = if runs_up || closed {
        x0
    } else {
        radial_of(end)
    };
    let at = |rho: f64, h: f64| seat.centre + from * rho + up * h;
    let corner = at(seat.radius, 0.0);
    let on_cap = at(ball.rho, 0.0);
    let on_wall = at(ball.wall_rho, ball.wall_height);
    let centre = at(ball.rho, ball.height);
    let face = section(
        model,
        [corner, on_cap, on_wall, centre],
        up.cross(from),
        radius,
        tol,
    )?;
    let axis = Axis::new(seat.centre, Direction::new(up, tol)?);
    let mut tool = ogeom_algo::make_revolution(model, &face, axis, sweep, tol)?.shape;
    let span = (seat.radius + 2.0 * radius) * 4.0;
    for cut in cuts {
        let Some(kept) = cut_past(model, &tool, *cut, span, tol)? else {
            return Ok(None);
        };
        tool = kept;
    }
    Ok(Some(tool))
}

/// A straight edge's tool: its section swept along it, ending square to it
/// at a cusp and in its end face's plane elsewhere.
fn straight_tool(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    faces: &[Shape; 2],
    radius: f64,
    cusps: &[Cusp],
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, bool)>> {
    let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
        return Ok(None);
    };
    let (pa, pb) = (point_of(model, &a)?, point_of(model, &b)?);
    let length = pa.distance(pb);
    let along = (pb - pa) / length;
    let (Some((_, n0)), Some((_, n1))) = (
        plane_of(model, &faces[0], tol)?,
        plane_of(model, &faces[1], tol)?,
    ) else {
        return Ok(None);
    };
    // Convex where the material lies behind both faces next to the edge.
    let probe = pa.midpoint(pb) - (n0 + n1) * (radius * 1e-3);
    let convex = matches!(
        ogeom_algo::classify_in_solid_exact(model, solid, probe, tol)?,
        ogeom_algo::Containment::In
    );
    let sign = if convex { -1.0 } else { 1.0 };
    // The ball's centre in the section through `pa`: `radius` from both
    // planes on the round's side, square to the edge.
    let (u0, u1) = (n0 - along * n0.dot(along), n1 - along * n1.dot(along));
    let g = u0.dot(u1);
    let det = 1.0 - g * g;
    if det <= 1e-12 {
        return Ok(None);
    }
    let (k0, k1) = ((1.0 - g) / det, (1.0 - g) / det);
    let reach = 2.0 * radius;
    let start = pa - along * reach;
    let offset = (u0 * k0 + u1 * k1) * (sign * radius);
    let centre = start + offset;
    let on0 = centre - u0 * (sign * radius);
    let on1 = centre - u1 * (sign * radius);
    let face = section(model, [start, on0, on1, centre], along, radius, tol)?;
    let mut tool = ogeom_algo::make_prism(model, &face, along * (length + 2.0 * reach), tol)?.shape;
    let span = (length + 4.0 * reach) * 2.0;
    for (v, p, leaving) in [(&a, pa, -along), (&b, pb, along)] {
        let cut = if cusps
            .iter()
            .any(|c| c.at.distance(p) <= tol.confusion() * 10.0)
        {
            (p, leaving)
        } else {
            let mut end = None;
            for f in ogeom_topo::ancestors_of(model, solid, v, ShapeType::Face)? {
                if faces.iter().any(|g| g.is_same(&f)) {
                    continue;
                }
                if let Some(plane) = plane_of(model, &f, tol)? {
                    end = Some(plane);
                }
            }
            match end {
                Some(plane) => plane,
                None => return Ok(None),
            }
        };
        let Some(kept) = cut_past(model, &tool, cut, span, tol)? else {
            return Ok(None);
        };
        tool = kept;
    }
    Ok(Some((tool, !convex)))
}

/// `once` rounded either side of the cusps its rims reach, each round its
/// own tool; `None` where an asked edge is not a rim on a drum or a
/// straight edge between planes, no rim reaches a cusp, a rim and a
/// straight edge leave a cusp on one side of it, a tool does not build, a
/// boolean refuses, or the result is not valid.
///
/// # Errors
///
/// As reading the edges and the booleans, other than their refusals.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut asked = Vec::with_capacity(once.len());
    let mut cusps: Vec<Cusp> = Vec::new();
    for (edge, _) in once {
        let Some(read) = read(model, solid, edge, tol)? else {
            return Ok(None);
        };
        if let Asked::Rim(seat, ..) = &read {
            for c in cusps_of(model, solid, seat, tol)? {
                if !cusps.iter().any(|d| d.at.distance(c.at) <= tol.confusion()) {
                    cusps.push(c);
                }
            }
        }
        asked.push(read);
    }
    if cusps.is_empty() {
        return Ok(None);
    }
    // A rim and a straight edge leaving a cusp the same way turn through
    // an S, which these tools do not build.
    let leaving = |model: &Model, edge: &Shape, at: Point| -> OgeomResult<Option<Vector>> {
        let (curve, (lo, hi)) = crate::support::edge_curve(model, edge, tol)?;
        for (t, s) in [(lo, 1.0), (hi, -1.0)] {
            if curve.point_at(t, tol)?.distance(at) <= tol.confusion() * 10.0 {
                return Ok(Some(curve.d1_at(t, tol)? * s));
            }
        }
        Ok(None)
    };
    for cusp in &cusps {
        let mut rims = Vec::new();
        let mut lines = Vec::new();
        for ((edge, _), read) in once.iter().zip(&asked) {
            if let Some(d) = leaving(model, edge, cusp.at)? {
                match read {
                    Asked::Rim(..) => rims.push(d),
                    Asked::Straight(_) => lines.push(d),
                }
            }
        }
        if rims.iter().any(|r| lines.iter().any(|l| r.dot(*l) > 0.0)) {
            // Where the top runs on past the cusp, its edges there left
            // sharp, an S on one side of the cusp has nothing to close
            // against: refused by name.
            let mut at_cusp = 0_usize;
            for e in explore_unique(model, &cusp.top, ShapeType::Edge)? {
                if let Some((a, b)) = ogeom_algo::edge_vertices(model, &e)?
                    && (point_of(model, &a)?.distance(cusp.at) <= tol.confusion() * 10.0
                        || point_of(model, &b)?.distance(cusp.at) <= tol.confusion() * 10.0)
                {
                    at_cusp += 1;
                }
            }
            if at_cusp > rims.len() + lines.len() {
                ogeom_core::ogeom_bail!(
                    Construction,
                    "a round turning through an S on one side of a cusp did not close: \
                     the top runs on past the tangent point with its edges there left \
                     sharp; round the chain on both sides of the tangent point; see \
                     docs/CAPABILITIES.md, fillet.edge-blends"
                );
            }
            return Ok(None);
        }
    }
    let mut tools = Vec::with_capacity(once.len());
    let mut scratch = model.clone();
    for ((edge, radius), read) in once.iter().zip(&asked) {
        let (solid_tool, additive) = match read {
            Asked::Rim(seat, circle, range) => {
                let own: Vec<(Point, Vector)> = cusps_of(&scratch, solid, seat, tol)?
                    .into_iter()
                    .map(|c| c.side)
                    .collect();
                let Some(t) = rim_tool(
                    &mut scratch,
                    (seat, circle, *range),
                    *radius,
                    (&own, false),
                    tol,
                )?
                else {
                    return Ok(None);
                };
                (t, seat.additive())
            }
            Asked::Straight(faces) => {
                match straight_tool(&mut scratch, solid, edge, faces, *radius, &cusps, tol)? {
                    Some(t) => t,
                    None => return Ok(None),
                }
            }
        };
        // The round's own face is named for its edge: the torus a rim's
        // ball sweeps, the cylinder a straight edge's does.
        let mut faces = Vec::new();
        for f in explore_unique(&scratch, &solid_tool, ShapeType::Face)? {
            let surface = scratch
                .node(&f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| scratch.geometry().surface(d.surface));
            let own = match read {
                Asked::Rim(..) => matches!(surface, Some(ogeom_geom::SurfaceGeometry::Torus(_))),
                Asked::Straight(_) => {
                    matches!(surface, Some(ogeom_geom::SurfaceGeometry::Cylinder(_)))
                }
            };
            if own {
                faces.push((f, edge.clone()));
            }
        }
        tools.push(Tool {
            solid: solid_tool,
            sense: if additive { Sense::Fuse } else { Sense::Cut },
            faces,
        });
    }
    // Applied in the solid's own order of the edges, whatever the
    // request's.
    let order = explore_unique(&scratch, solid, ShapeType::Edge)?;
    let rank = |e: &Shape| {
        order
            .iter()
            .position(|o| o.is_same(e))
            .unwrap_or(usize::MAX)
    };
    let mut ranked: Vec<(usize, Tool)> = once.iter().map(|(e, _)| rank(e)).zip(tools).collect();
    ranked.sort_by_key(|(k, _)| *k);
    let tools: Vec<Tool> = ranked.into_iter().map(|(_, t)| t).collect();
    let edges: Vec<Shape> = once.iter().map(|(e, _)| e.clone()).collect();
    let Some(built) = super::fallback::applied(&mut scratch, solid, &edges, &tools, tol)? else {
        return Ok(None);
    };
    if !ogeom_algo::check(&scratch, &built.shape, tol)?.is_valid() {
        return Ok(None);
    }
    *model = scratch;
    Ok(Some(built))
}

/// One rim arc of a drum ending obliquely on flats (a drum cut flat across
/// its rim), rounded alone (ADR 0001, amendment "a rim arc ending obliquely
/// on a flat"): the ball rolls on past each end through the corner it still
/// finds in the solid, while its contact on the cap stands on the cap. The
/// round's tool is the rim's section turned a whole turn about the drum's
/// axis, cut by the plane of each face the arc ends on. The arc may come in
/// pieces of one circle (a wall's seam), and be concave (a drum's foot on a
/// block's top, filled); an end's plane may stand square to the arc (a
/// meridian), where the round ends flush in it. `None` where the request is
/// not one such arc, the cap and the wall share an edge not asked, an end's
/// face is not one plane, or holds the rim's way on, the tool does not
/// build, the boolean refuses, or the result is not valid.
///
/// # Errors
///
/// As reading the edge and the booleans, other than their refusals.
pub(crate) fn oblique_rim_rounds(
    model: &mut Model,
    solid: &Shape,
    once: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some((edge, radius)) = once.first() else {
        return Ok(None);
    };
    let Some(Asked::Rim(seat, circle, range)) = read(model, solid, edge, tol)? else {
        return Ok(None);
    };
    // The arc may come in pieces of one circle (a wall's seam splits a
    // drum's foot), asked with one radius: the arc ends where a piece's
    // vertex is no other piece's.
    let mut vertices: Vec<Shape> = Vec::new();
    for (piece, r) in once {
        if (r - radius).abs() > tol.confusion() {
            return Ok(None);
        }
        if !piece.is_same(edge) {
            let Some(Asked::Rim(other, c, _)) = read(model, solid, piece, tol)? else {
                return Ok(None);
            };
            let same = |p: Point, q: Point| p.distance(q) <= tol.confusion();
            if !(same(other.centre, seat.centre)
                && (other.radius - seat.radius).abs() <= tol.confusion()
                && other.up.dot(seat.up) >= 1.0 - tol.angular()
                && (c.circle().radius() - circle.circle().radius()).abs() <= tol.confusion()
                && other.additive() == seat.additive())
            {
                return Ok(None);
            }
        }
        let Some((a, b)) = ogeom_algo::edge_vertices(model, piece)? else {
            return Ok(None);
        };
        if a.is_same(&b) {
            return Ok(None);
        }
        for v in [a, b] {
            match vertices.iter().position(|w| w.is_same(&v)) {
                Some(i) => {
                    vertices.remove(i);
                }
                None => vertices.push(v),
            }
        }
    }
    let [a, b] = vertices.as_slice() else {
        return Ok(None);
    };
    // The section turned a whole turn rounds every edge the cap and the
    // wall share: those must be the pieces asked (a drum cut across its
    // rim leaves it one arc; a drum cutting a block's corners leaves it
    // an arc at each).
    let cap_edges = explore_unique(model, &seat.cap_face, ShapeType::Edge)?;
    for shared in explore_unique(model, &seat.wall_face, ShapeType::Edge)? {
        if cap_edges.iter().any(|e| e.is_same(&shared))
            && !once.iter().any(|(piece, _)| piece.is_same(&shared))
        {
            return Ok(None);
        }
    }
    let mut cuts: Vec<(Point, Vector)> = Vec::new();
    for v in [a, b] {
        let p = point_of(model, v)?;
        // The way the circle runs on past the arc's end at `v`: out of
        // the arc, away from the pieces.
        let radial = {
            let d = p - seat.centre;
            let r = d - seat.up * d.dot(seat.up);
            r / r.magnitude()
        };
        let tangent = seat.up.cross(radial);
        let into_arc = once.iter().find_map(|(piece, _)| {
            let (c, (lo, hi)) = crate::support::edge_curve(model, piece, tol).ok()?;
            let Curve::Circle(c) = c else { return None };
            let mid = c.point_at(0.5 * (lo + hi), tol).ok()?;
            let ends = ogeom_algo::edge_vertices(model, piece).ok()??;
            (ends.0.is_same(v) || ends.1.is_same(v)).then_some(mid - p)
        });
        let Some(into_arc) = into_arc else {
            return Ok(None);
        };
        let ahead = if tangent.dot(into_arc) > 0.0 {
            -tangent
        } else {
            tangent
        };
        let mut ends = Vec::new();
        for f in ogeom_topo::ancestors_of(model, solid, v, ShapeType::Face)? {
            if once.iter().any(|(piece, _)| {
                faces_at(model, solid, piece).is_ok_and(|fs| fs.iter().any(|g| g.is_same(&f)))
            }) {
                continue;
            }
            if !ends.iter().any(|g: &Shape| g.is_same(&f)) {
                ends.push(f);
            }
        }
        let [end] = ends.as_slice() else {
            return Ok(None);
        };
        let Some((on, out)) = plane_of(model, end, tol)? else {
            return Ok(None);
        };
        // The face's plane must not hold the rim's way on (square to the
        // meridian there): oblique, or square to the rim (a meridian
        // plane, where the round ends flush in it).
        let lean = out.dot(ahead) / ahead.magnitude();
        if lean.abs() <= tol.angular() {
            return Ok(None);
        }
        if (p - on).dot(out).abs() > tol.confusion() * 10.0 {
            return Ok(None);
        }
        if !cuts.iter().any(|(o, n)| {
            n.dot(out) >= 1.0 - tol.angular() && (*o - on).dot(out).abs() <= tol.confusion()
        }) {
            cuts.push((on, out));
        }
    }
    let mut scratch = model.clone();
    let Some(tool) = rim_tool(
        &mut scratch,
        (&seat, &circle, range),
        *radius,
        (&cuts, true),
        tol,
    )?
    else {
        return Ok(None);
    };
    let mut faces = Vec::new();
    for f in explore_unique(&scratch, &tool, ShapeType::Face)? {
        let torus = scratch
            .node(&f)
            .and_then(|n| n.data().as_face())
            .and_then(|d| scratch.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Torus(_)));
        if torus {
            for (piece, _) in once {
                faces.push((f.clone(), piece.clone()));
            }
        }
    }
    let tools = [Tool {
        solid: tool,
        sense: if seat.additive() {
            Sense::Fuse
        } else {
            Sense::Cut
        },
        faces,
    }];
    let edges: Vec<Shape> = once.iter().map(|(e, _)| e.clone()).collect();
    let Some(built) = super::fallback::applied(&mut scratch, solid, &edges, &tools, tol)? else {
        return Ok(None);
    };
    if !ogeom_algo::check(&scratch, &built.shape, tol)?.is_valid() {
        return Ok(None);
    }
    *model = scratch;
    Ok(Some(built))
}
