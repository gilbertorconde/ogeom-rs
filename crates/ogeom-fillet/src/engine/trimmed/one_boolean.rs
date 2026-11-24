//! Chamfers on curved solids closed into tools and applied in one boolean.

use super::*;

/// Curved chamfers the rings do not place (a bevel rolling over its face's
/// far side, or reaching a face it does not end at) built by the one
/// boolean (ADR 0001, amendments "M5a" and "M5b"): each asked edge's
/// bevel, read as the curved path reads it, closed into its tool by
/// [`close_tool`] with its own ends, the tools applied together. Only
/// where no two asked edges meet and each edge's ends are each on one
/// other face, a plane square to the edge; `None` otherwise, or where a
/// tool will not close, a boolean refuses or a bevel is left nothing of.
///
/// # Errors
///
/// As reading the solid and building the tools and the booleans, other
/// than their refusals.
pub(crate) fn creases_one_boolean(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    // Each asked edge's ends: a joint where it runs on into another asked
    // piece of its curve between the same two faces (capped there on the
    // plane square to the curve, as two features on the pieces would end),
    // else its own end. Asked edges meeting any other way decline.
    let asked: Vec<Id> = specs.iter().map(|(e, _)| id(e)).collect();
    let mut caps: Vec<[EndCap; 2]> = Vec::with_capacity(specs.len());
    // Whether asked edges run on into each other tangentially across
    // different faces (a pocket floor's loop of lines and arcs): there a
    // bevel may be concave, its tool fused in.
    let mut either_sense = false;
    for e in &asked {
        let Some(c) = curved.edges.get(e) else {
            return Ok(None);
        };
        if c.ends.0 == c.ends.1 {
            return Ok(None);
        }
        // Two faces on one carrier, facing one way: a seam between them.
        let one = |a: usize, b: usize| -> OgeomResult<bool> {
            ogeom_heal::same_domain(model, &curved.faces[a].shape, &curved.faces[b].shape, tol)
        };
        let own = faces_of(&curved, *e);
        let mut pair = [EndCap::Own, EndCap::Own];
        for (k, v) in [c.ends.0, c.ends.1].into_iter().enumerate() {
            let others: Vec<Id> = asked
                .iter()
                .copied()
                .filter(|x| {
                    x != e && {
                        let o = &curved.edges[x];
                        o.ends.0 == v || o.ends.1 == v
                    }
                })
                .collect();
            match others[..] {
                [] => {}
                [o] => {
                    // The same two surfaces either side, and nothing else at
                    // the joint but seams (edges between faces on one
                    // surface, which a boolean left).
                    let beside = faces_of(&curved, o);
                    let [a0, a1] = own[..] else {
                        return Ok(None);
                    };
                    let [b0, b1] = beside[..] else {
                        return Ok(None);
                    };
                    let same = (one(a0, b0)? && one(a1, b1)?) || (one(a0, b1)? && one(a1, b0)?);
                    let mut seams_only = true;
                    for (x, y) in &curved.edges {
                        if *x == *e || *x == o || !(y.ends.0 == v || y.ends.1 == v) {
                            continue;
                        }
                        let [f, g] = faces_of(&curved, *x)[..] else {
                            seams_only = false;
                            break;
                        };
                        if !one(f, g)? {
                            seams_only = false;
                            break;
                        }
                    }
                    if !(same && seams_only) {
                        if !tangent_joint(model, &curved, (*e, o), v, tol)? {
                            return Ok(None);
                        }
                        either_sense = true;
                    }
                    let t = if k == 0 { c.range.0 } else { c.range.1 };
                    let h = (c.range.1 - c.range.0) * 1e-6;
                    let p = c.curve.point_at(t, tol)?;
                    let q = c.curve.point_at(if k == 0 { t + h } else { t - h }, tol)?;
                    pair[k] = EndCap::Plane {
                        at: p,
                        normal: q - p,
                    };
                }
                _ => return Ok(None),
            }
        }
        caps.push(pair);
    }
    let mut tools = Vec::with_capacity(specs.len());
    for (edge, spec) in specs {
        let e = id(edge);
        let strip = match bevel_of(model, solid, &curved, e, edge, spec, either_sense, tol)? {
            Some(strip) => strip,
            None => {
                let [f0, f1] = faces_of(&curved, e)[..] else {
                    return Ok(None);
                };
                let Some((strip, _)) =
                    general_crease(model, solid, &curved, e, edge, [f0, f1], spec, None, tol)?
                else {
                    return Ok(None);
                };
                strip
            }
        };
        let Some(tool) = close_tool(model, &curved, &strip, caps[tools.len()], tol)? else {
            return Ok(None);
        };
        tools.push(tool);
    }
    let edges: Vec<Shape> = specs.iter().map(|(e, _)| e.clone()).collect();
    super::fallback::applied(model, solid, &edges, &tools, tol)
}

/// Whether asked edges `a` and `b` run on into each other at `v` as one
/// tangent chain across different faces: their curves leave `v` back to
/// back (within [`super::setbacks::TANGENT`]), they share one face (or two
/// on one carrier facing one way), their other faces stand tangent at `v`,
/// and every other edge at `v` is a seam or the smooth edge between those
/// two. The bevels' sections then meet on the plane square to the chain
/// there, and each tool is capped on it.
fn tangent_joint(
    model: &Model,
    curved: &Curved,
    (a, b): (Id, Id),
    v: Id,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let leaving = |e: Id| -> OgeomResult<Vector> {
        let c = &curved.edges[&e];
        let (t, sign) = if c.ends.0 == v {
            (c.range.0, 1.0)
        } else {
            (c.range.1, -1.0)
        };
        let d = c.curve.d1_at(t, tol)? * sign;
        Ok(d / d.magnitude())
    };
    let (ta, tb) = (leaving(a)?, leaving(b)?);
    if ta.dot(tb) >= 0.0 || ta.cross(tb).magnitude() > super::setbacks::TANGENT {
        return Ok(false);
    }
    let one = |x: usize, y: usize| -> OgeomResult<bool> {
        Ok(x == y
            || ogeom_heal::same_domain(model, &curved.faces[x].shape, &curved.faces[y].shape, tol)?)
    };
    let (&[a0, a1], &[b0, b1]) = (&faces_of(curved, a)[..], &faces_of(curved, b)[..]) else {
        return Ok(false);
    };
    // The faces the two edges do not share, one each side.
    let mut others = None;
    for (x, xo) in [(a0, a1), (a1, a0)] {
        for (y, yo) in [(b0, b1), (b1, b0)] {
            if others.is_none() && one(x, y)? && !one(xo, yo)? {
                others = Some((xo, yo));
            }
        }
    }
    let Some((fa, fb)) = others else {
        return Ok(false);
    };
    let p = curved.points[&v];
    let (Some(na), Some(nb)) = (outward(curved, fa, p), outward(curved, fb, p)) else {
        return Ok(false);
    };
    if na.dot(nb) <= 0.0
        || na.cross(nb).magnitude() > super::setbacks::TANGENT * na.magnitude() * nb.magnitude()
    {
        return Ok(false);
    }
    for (x, c) in &curved.edges {
        if *x == a || *x == b || !(c.ends.0 == v || c.ends.1 == v) {
            continue;
        }
        let [f, g] = faces_of(curved, *x)[..] else {
            return Ok(false);
        };
        let between = (f == fa && g == fb) || (f == fb && g == fa);
        if !between && !one(f, g)? {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Rim chamfers between two surfaces of revolution about the rim's own
/// axis (a plane square to it, a coaxial cylinder, cone, sphere or torus)
/// that no other route builds, by the one boolean (ADR 0001, amendment
/// "M5b"): a rim rolling over or reaching a face, or pieces of one rim
/// asked apart. Each tool is the bevel's section, in the half-plane
/// through the axis at the arc's start (the rim, each face's trace to its
/// contact, the bevel's line between the contacts), turned about the axis
/// over the arc's sweep, its ends on the planes through the axis (the
/// step at a joint); arcs of one circle asked alike that close it turn
/// once, whole, and unless pieces are asked apart an arc asked alone whose
/// circle's pieces between the same two carriers close the turn stands for
/// the whole rim. Cut where the rim is convex, fused where concave; the
/// bevel's face is named for every edge of its turn. `None` where an
/// asked edge is not such a rim, its strip's contacts are not circles
/// about its axis, a tool will not build, a boolean refuses or a bevel is
/// left nothing of.
///
/// # Errors
///
/// As reading the solid and its strips, other than their refusals.
pub(crate) fn rims_revolved(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    apart: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let (eps, ang) = (tol.confusion(), tol.angular());
    if specs.is_empty() {
        return Ok(None);
    }
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    // Each rim: its axis, its start and sweep about it, and its section's
    // points (the rim, the contact on each face) at the start.
    struct Rim {
        origin: Point,
        axis: Vector,
        radius: f64,
        start: Point,
        sweep: f64,
        faces: [usize; 2],
        contacts: [(Point, f64); 2],
    }
    let mut rims: Vec<Rim> = Vec::with_capacity(specs.len());
    for (edge, spec) in specs {
        let e = id(edge);
        let Some(c) = curved.edges.get(&e) else {
            return Ok(None);
        };
        let Curve::Circle(circle) = &c.curve else {
            return Ok(None);
        };
        let Some(strip) = bevel_of(model, solid, &curved, e, edge, spec, false, tol)? else {
            return Ok(None);
        };
        let circle = circle.circle();
        let (origin, axis) = (circle.centre(), circle.frame().z().vector());
        let off_axis = |p: Point| {
            let q = p - origin;
            (q - axis * q.dot(axis)).magnitude()
        };
        // Each contact a circle about the rim's axis: its centre's height
        // along the axis and its radius.
        let mut contacts = [(Point::ORIGIN, 0.0); 2];
        for (k, contact) in strip.contacts.iter().enumerate() {
            let Curve::Circle(ring) = contact else {
                return Ok(None);
            };
            let ring = ring.circle();
            if ring.frame().z().vector().cross(axis).magnitude() > ang
                || off_axis(ring.centre()) > eps
            {
                return Ok(None);
            }
            contacts[k] = (ring.centre(), ring.radius());
        }
        let closed = c.ends.0 == c.ends.1
            || !apart
                && !asked_apart(&curved, specs, e, tol)
                && pieces_close(model, &curved, e, tol)?;
        let (t0, t1) = c.range;
        rims.push(Rim {
            origin,
            axis,
            radius: circle.radius(),
            start: c.curve.point_at(t0, tol)?,
            sweep: if closed {
                core::f64::consts::TAU
            } else {
                t1 - t0
            },
            faces: strip.faces,
            contacts,
        });
    }
    // Arcs of one circle asked alike turn once, whole, where they close it.
    let same = |i: &Rim, j: &Rim| {
        i.origin.distance(j.origin) <= eps
            && (i.radius - j.radius).abs() <= eps
            && i.axis.dot(j.axis) >= 1.0 - ang
            && (0..2).all(|k| {
                i.contacts[k].0.distance(j.contacts[k].0) <= eps
                    && (i.contacts[k].1 - j.contacts[k].1).abs() <= eps
            })
    };
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for i in 0..rims.len() {
        match groups.iter_mut().find(|g| same(&rims[g[0]], &rims[i])) {
            Some(g) => g.push(i),
            None => groups.push(vec![i]),
        }
    }
    let mut runs: Vec<(Vec<usize>, f64)> = Vec::new();
    for group in groups {
        let turned: f64 = group.iter().map(|k| rims[*k].sweep).sum();
        let whole = |k: &usize| (rims[*k].sweep - core::f64::consts::TAU).abs() <= ang;
        if group.iter().any(whole)
            || group.len() > 1 && (turned - core::f64::consts::TAU).abs() <= ang
        {
            runs.push((group, core::f64::consts::TAU));
        } else {
            runs.extend(group.into_iter().map(|k| (vec![k], rims[k].sweep)));
        }
    }
    let mut tools = Vec::with_capacity(runs.len());
    for (run, sweep) in runs {
        let rim = &rims[run[0]];
        let Some(tool) = revolved_tool(model, &curved, rim_section(rim), sweep, tol)? else {
            return Ok(None);
        };
        let mut faces = Vec::new();
        for k in &run {
            faces.push((tool.1.clone(), specs[*k].0.clone()));
        }
        tools.push(super::fallback::Tool {
            solid: tool.0,
            sense: tool.2,
            faces,
        });
    }
    let edges: Vec<Shape> = specs.iter().map(|(e, _)| e.clone()).collect();
    return super::fallback::applied(model, solid, &edges, &tools, tol);

    /// A rim's section at its start: the axis, the rim's point, its faces
    /// and the contact on each in the same half-plane.
    fn rim_section(rim: &Rim) -> RimSection {
        let q = rim.start - rim.origin;
        let radial = q - rim.axis * q.dot(rim.axis);
        let u = radial / radial.magnitude();
        let at = |(centre, radius): (Point, f64)| centre + u * radius;
        RimSection {
            origin: rim.origin,
            axis: rim.axis,
            radial: u,
            rim: rim.start,
            faces: rim.faces,
            contacts: [at(rim.contacts[0]), at(rim.contacts[1])],
        }
    }
}

/// Chamfers on a convex solid whose bevels meet and roll over or reach a
/// face, by booleans (ADR 0001, amendment "M5b"): on a convex solid a
/// bevel's notch takes everything past its surface, so what is left is
/// the solid in common with each bevel's material side, a half-space for
/// a plane and the inside of a cone or cylinder about a rim's axis. Only
/// where every face is a plane or a cylinder, cone or sphere with the
/// material inside, every edge is convex at its samples (a connected
/// solid convex near every point is convex), each asked edge's bevel is
/// such a plane, cone or cylinder with the edge outside it, and no vertex
/// has three asked edges (those close with a facet). Each bevel's face is
/// named for its edge. `None` otherwise or where a boolean refuses.
///
/// # Errors
///
/// As reading the solid and its strips, other than their refusals.
pub(crate) fn convex_one_boolean(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let ang = tol.angular();
    if specs.is_empty() {
        return Ok(None);
    }
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    if !convex_solid(model, &curved, tol)? {
        return Ok(None);
    }
    // No vertex of three asked edges.
    let mut at_vertex: FastMap<Id, usize> = FastMap::default();
    for (edge, _) in specs {
        let Some(c) = curved.edges.get(&id(edge)) else {
            return Ok(None);
        };
        for v in [c.ends.0, c.ends.1] {
            *at_vertex.entry(v).or_default() += 1;
        }
        if c.ends.0 == c.ends.1 {
            *at_vertex.entry(c.ends.0).or_default() -= 1;
        }
    }
    if at_vertex.values().any(|n| *n >= 3) {
        return Ok(None);
    }
    let bounds = ogeom_algo::shape_bounds(model, solid, tol)?;
    let (Some(centre), reach) = (bounds.centre(), 2.0 * bounds.diagonal()) else {
        return Ok(None);
    };
    let margin = 0.01 * bounds.diagonal() + tol.confusion();
    let corners = bounds.corners();
    // Each bevel's material side as a solid holding the whole solid but
    // what the bevel takes, and its face on the bevel's surface.
    let mut sides: Vec<(Shape, Vec<Shape>, Shape)> = Vec::with_capacity(specs.len());
    for (edge, spec) in specs {
        let e = id(edge);
        let Some(strip) = bevel_of(model, solid, &curved, e, edge, spec, false, tol)? else {
            return Ok(None);
        };
        let c = &curved.edges[&e];
        let mid = c.curve.point_at(f64::midpoint(c.range.0, c.range.1), tol)?;
        let side = match &strip.surface {
            SurfaceGeometry::Plane(plane) => {
                let plane = plane.plane();
                let n = plane.normal().vector();
                let off = (mid - plane.origin()).dot(n);
                if off.abs() <= tol.confusion() {
                    return Ok(None);
                }
                // Into the material: away from the edge.
                let inward = n * -off.signum();
                let foot = centre - n * (centre - plane.origin()).dot(n);
                let z = Direction::new(inward, tol)?;
                let x = Direction::new(
                    if inward.x.abs() < 0.9 {
                        Vector::X
                    } else {
                        Vector::Y
                    },
                    tol,
                )?;
                let frame = Frame::new(foot, z, x, tol)?;
                let origin = foot - frame.x().vector() * reach - frame.y().vector() * reach;
                let frame = Frame::new(origin, z, frame.x(), tol)?;
                let block = ogeom_algo::make_box(
                    model,
                    frame,
                    (2.0 * reach, 2.0 * reach, 2.0 * reach),
                    tol,
                )?
                .shape;
                let mut on = Vec::new();
                for face in ogeom_topo::explore_unique(model, &block, ogeom_topo::ShapeType::Face)?
                {
                    let (at, out) = ogeom_algo::face_normal(model, &face, tol)?;
                    if out.dot(inward) <= -1.0 + ang
                        && (at - foot).dot(inward).abs() <= tol.confusion()
                    {
                        on.push(face);
                    }
                }
                (block, on)
            }
            SurfaceGeometry::Cone(cone) => {
                let cone = cone.cone();
                let Some(side) = revolved_side(
                    model,
                    cone.frame(),
                    |z| cone.radius_at(z),
                    (&corners, margin, mid),
                    tol,
                )?
                else {
                    return Ok(None);
                };
                side
            }
            SurfaceGeometry::Cylinder(cylinder) => {
                let cylinder = cylinder.cylinder();
                let radius = cylinder.radius();
                let Some(side) = revolved_side(
                    model,
                    cylinder.frame(),
                    |_| radius,
                    (&corners, margin, mid),
                    tol,
                )?
                else {
                    return Ok(None);
                };
                side
            }
            _ => return Ok(None),
        };
        sides.push((side.0, side.1, edge.clone()));
    }
    let mut built = Built::from_nothing(solid.clone());
    for (side, _, _) in &sides {
        let Ok(step) = ogeom_bool::common(model, &built.shape, side, tol) else {
            return Ok(None);
        };
        built = Built {
            shape: step.shape.clone(),
            history: built.history.then(&step.history),
        };
    }
    let kept = ogeom_topo::explore_unique(model, &built.shape, ogeom_topo::ShapeType::Face)?;
    let mut history: History = built.history.clone();
    for (_, faces, from) in &sides {
        for face in faces {
            for made in built
                .history
                .modified(face)
                .iter()
                .chain(built.history.generated(face))
            {
                if kept.iter().any(|k| k.is_same(made)) {
                    history.generate(from, made.clone());
                }
            }
        }
    }
    for (edge, _) in specs {
        if history.generated(edge).is_empty() {
            return Ok(None);
        }
        history.delete(edge);
    }
    Ok(Some(Built::new(built.shape, history)))
}

/// The inside of a cone or cylinder about `frame`'s `z`, radius `radius`
/// at height `z` along it, as a solid over the heights the solid's bounds
/// `corners` reach (and `margin` either side), with its face on the
/// surface; `None` where the edge's middle `mid` is not outside it or the
/// apex falls among those heights.
#[allow(clippy::type_complexity, reason = "a solid and its named faces")]
fn revolved_side(
    model: &mut Model,
    frame: Frame,
    radius: impl Fn(f64) -> f64,
    (corners, margin, mid): (&[Point], f64, Point),
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
    let (o, a) = (frame.origin(), frame.z().vector());
    let heights: Vec<f64> = corners.iter().map(|p| (*p - o).dot(a)).collect();
    let low = heights.iter().copied().fold(f64::INFINITY, f64::min) - margin;
    let high = heights.iter().copied().fold(f64::NEG_INFINITY, f64::max) + margin;
    let (r_low, r_high) = (radius(low), radius(high));
    if !(r_low > tol.confusion() && r_high > tol.confusion()) {
        return Ok(None);
    }
    let q = mid - o;
    let h = q.dot(a);
    if (q - a * h).magnitude() <= radius(h) + tol.confusion() {
        return Ok(None);
    }
    let base = Frame::new(o + a * low, frame.z(), frame.x(), tol)?;
    let block = if (r_low - r_high).abs() <= tol.confusion() {
        ogeom_algo::make_cylinder(model, base, r_low, high - low, tol)?.shape
    } else {
        ogeom_algo::make_cone(model, base, r_low, r_high, high - low, tol)?.shape
    };
    let mut on = Vec::new();
    for face in ogeom_topo::explore_unique(model, &block, ogeom_topo::ShapeType::Face)? {
        let (_, out) = ogeom_algo::face_normal(model, &face, tol)?;
        if out.dot(a).abs() < 1.0 - tol.angular() {
            on.push(face);
        }
    }
    Ok(Some((block, on)))
}

/// Whether `curved` is convex: every face a plane, or a cylinder, cone or
/// sphere with the material inside, and every edge between two faces
/// convex at a quarter, the middle and three quarters of the way along it
/// (tangent counts as convex). A connected solid convex near every point
/// is convex.
fn convex_solid(model: &Model, curved: &Curved, tol: Tolerances) -> OgeomResult<bool> {
    let ang = tol.angular();
    for (fi, f) in curved.faces.iter().enumerate() {
        let Some((e, _)) = f.loops.iter().flatten().next() else {
            return Ok(false);
        };
        let c = &curved.edges[e];
        let p = c.curve.point_at(f64::midpoint(c.range.0, c.range.1), tol)?;
        let Some(n) = outward(curved, fi, p) else {
            return Ok(false);
        };
        let radial = match &f.surface {
            SurfaceGeometry::Plane(_) => continue,
            SurfaceGeometry::Cylinder(s) => {
                let axis = s.cylinder().axis();
                let q = p - axis.location;
                q - axis.direction.vector() * q.dot(axis.direction.vector())
            }
            SurfaceGeometry::Cone(s) => {
                let axis = s.cone().axis();
                let q = p - axis.location;
                q - axis.direction.vector() * q.dot(axis.direction.vector())
            }
            SurfaceGeometry::Sphere(s) => p - s.sphere().centre(),
            _ => return Ok(false),
        };
        if n.dot(radial) <= 0.0 {
            return Ok(false);
        }
    }
    for (e, c) in &curved.edges {
        let [f0, f1] = faces_of(curved, *e)[..] else {
            continue;
        };
        for share in [0.25, 0.5, 0.75] {
            let t = c.range.0 + (c.range.1 - c.range.0) * share;
            let p = c.curve.point_at(t, tol)?;
            let (Some(into), Some(n1)) = (
                into_face(model, curved, (f0, *e), (p, t), tol)?,
                outward(curved, f1, p),
            ) else {
                return Ok(false);
            };
            if n1.dot(into) > ang {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// The unit direction into face `face` from its edge `edge` at `p`, the
/// edge's point at `t`: square to the edge in the face's tangent plane, the
/// way a short step (a hundredth of the edge's span) lands inside the face
/// (read on the face itself, not from its loops' orientation). `None` where
/// neither way, or both, does.
pub(super) fn into_face(
    model: &Model,
    curved: &Curved,
    (face, edge): (usize, Id),
    (p, t): (Point, f64),
    tol: Tolerances,
) -> OgeomResult<Option<Vector>> {
    let c = &curved.edges[&edge];
    let d = c.curve.d1_at(t, tol)?;
    let Some(n) = outward(curved, face, p) else {
        return Ok(None);
    };
    let across = n.cross(d);
    let m = across.magnitude();
    if m <= 0.0 {
        return Ok(None);
    }
    let across = across / m;
    let span = c
        .curve
        .point_at(c.range.0, tol)?
        .distance(c.curve.point_at(c.range.1, tol)?)
        .max(1.0);
    // Clear of the face's boundary band, still close to the edge.
    let step = (span * 1e-2).max(100.0 * tol.confusion());
    let shape = &curved.faces[face].shape;
    let lands = |dir: Vector| -> OgeomResult<bool> {
        let q = p + dir * step;
        let foot =
            ogeom_geom::project::project_on_surface(&curved.faces[face].surface, q, 16, tol)?;
        // The trim drawn finer than the step, so its uncertain band (the
        // tolerance and the chord's sag) leaves the step clear.
        let fine = ogeom_mesh::Deflection {
            chord: step * 0.1,
            ..ogeom_mesh::Deflection::default()
        };
        Ok(
            ogeom_algo::classify_on_face(model, shape, foot.point, fine, tol)?
                == ogeom_algo::Containment::In,
        )
    };
    match (lands(across)?, lands(-across)?) {
        (true, false) => Ok(Some(across)),
        (false, true) => Ok(Some(-across)),
        _ => Ok(None),
    }
}
