//! A convex straight edge between a plane and a cylinder along its ruling
//! (a cube's side meeting the drum that shaves it), rounded where an end
//! meets an earlier round's band tangent to one of its faces: its band runs
//! on past that end and the earlier band trims it (ADR 0001, amendments "a
//! later round reaching an earlier band at its end" and "a later round
//! reaching a band seated on a curved face"). The band's section is the
//! same all along the ruling, so its spandrel is a prism: carried on past
//! both ends, cut from the solid as it stands.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{CircleCurve, Curve, Curve3d as _, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::curved::{Curved, read};
use super::fallback::{Sense, Tool};
use super::trimmed::{faces_of, outward};
use super::{Id, id};

/// How many stations along the ruling its ball's seat is read at, and how
/// many points across the way to each contact.
const ALONG: u32 = 16;
const ACROSS: u32 = 8;

/// How far apart two faces' normals may stand and the faces still touch
/// tangentially at a vertex a boolean left a hair off the contact.
const TANGENT: f64 = 1e-6;

/// The ruling read: its line, its plane face and cylinder face, and the
/// ball's section at a point of it.
struct Ruling {
    shape: Shape,
    ends: (Id, Id),
    a: Point,
    u: Vector,
    length: f64,
    plane: usize,
    wall: usize,
    radius: f64,
}

/// The rounds of `asked` (one edge) where it is such a ruling with an end
/// at an earlier band; `None` otherwise, or where a tool or a boolean will
/// not build.
///
/// # Errors
///
/// As reading the solid, classifying points and building the tool.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    asked: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let [(edge, radius)] = asked else {
        {
            return Ok(None);
        }
    };
    let Some(curved) = read(model, solid, tol)? else {
        {
            return Ok(None);
        }
    };
    let Some(ruling) = ruling_of(&curved, edge, *radius, tol)? else {
        {
            return Ok(None);
        }
    };
    if !at_earlier_band(&curved, &ruling, tol)? {
        {
            return Ok(None);
        }
    }
    if !seated(model, solid, &curved, &ruling, tol)? {
        {
            return Ok(None);
        }
    }
    // Carried on along the ruling past both ends as far as the solid
    // reaches: the cut takes only what the solid holds.
    let (mut lo, mut hi) = (
        Vector::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
        Vector::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for p in curved.points.values() {
        let v = *p - Point::ORIGIN;
        lo = Vector::new(lo.x.min(v.x), lo.y.min(v.y), lo.z.min(v.z));
        hi = Vector::new(hi.x.max(v.x), hi.y.max(v.y), hi.z.max(v.z));
    }
    let reach = (hi - lo).magnitude() + 4.0 * ruling.radius;
    let mut scratch = model.clone();
    let Some((tool, band)) = prism(&mut scratch, &curved, &ruling, reach, tol)? else {
        {
            return Ok(None);
        }
    };
    let Ok(taken) = ogeom_bool::common(&mut scratch, &tool, solid, tol) else {
        {
            return Ok(None);
        }
    };
    if explore_unique(&scratch, &taken.shape, ShapeType::Solid)?.len() != 1 {
        {
            return Ok(None);
        }
    }
    let tools = [Tool {
        solid: tool,
        sense: Sense::Cut,
        faces: vec![(band, ruling.shape.clone())],
    }];
    let Some(built) = super::fallback::applied(
        &mut scratch,
        solid,
        std::slice::from_ref(&ruling.shape),
        &tools,
        tol,
    )?
    else {
        {
            return Ok(None);
        }
    };
    match super::guards::solid_fits(&scratch, &built.shape, tol) {
        Ok(()) => {
            *model = scratch;
            Ok(Some(built))
        }
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(None),
    }
}

/// `edge` as a convex ruling between a plane and a cylinder whose axis
/// runs along it, read on its own line.
fn ruling_of(
    curved: &Curved,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Ruling>> {
    let Some(c) = curved.edges.get(&id(edge)) else {
        {
            return Ok(None);
        }
    };
    if c.ends.0 == c.ends.1 || !matches!(c.curve, Curve::Line(_)) {
        {
            return Ok(None);
        }
    }
    if !(radius.is_finite() && radius > tol.confusion()) {
        {
            return Ok(None);
        }
    }
    let (a, b) = (
        c.curve.point_at(c.range.0, tol)?,
        c.curve.point_at(c.range.1, tol)?,
    );
    let length = a.distance(b);
    if length <= tol.confusion() {
        {
            return Ok(None);
        }
    }
    let u = (b - a) / length;
    let [fa, fb] = faces_of(curved, id(edge))[..] else {
        {
            return Ok(None);
        }
    };
    let (plane, wall) = match (&curved.faces[fa].surface, &curved.faces[fb].surface) {
        (SurfaceGeometry::Plane(_), SurfaceGeometry::Cylinder(_)) => (fa, fb),
        (SurfaceGeometry::Cylinder(_), SurfaceGeometry::Plane(_)) => (fb, fa),
        _ => return Ok(None),
    };
    let SurfaceGeometry::Cylinder(cyl) = &curved.faces[wall].surface else {
        {
            return Ok(None);
        }
    };
    if cyl
        .cylinder()
        .axis()
        .direction
        .vector()
        .cross(u)
        .magnitude()
        > tol.angular()
    {
        {
            return Ok(None);
        }
    }
    let ruling = Ruling {
        shape: edge.clone(),
        ends: c.ends,
        a,
        u,
        length,
        plane,
        wall,
        radius,
    };
    // Convex, the material inside the drum and behind the plane: the ball
    // sits there, a radius off each.
    let mid = a + u * (length * 0.5);
    let Some([_, _, _, _]) = ball(curved, &ruling, mid, tol)? else {
        {
            return Ok(None);
        }
    };
    Ok(Some(ruling))
}

/// At `at` on the ruling: the edge's point, the ball's centre, its contact
/// on the plane and its contact on the cylinder; `None` where no ball sits
/// in the convex corner there.
fn ball(
    curved: &Curved,
    ruling: &Ruling,
    at: Point,
    tol: Tolerances,
) -> OgeomResult<Option<[Point; 4]>> {
    let (Some(np), Some(nc)) = (
        outward(curved, ruling.plane, at),
        outward(curved, ruling.wall, at),
    ) else {
        return Ok(None);
    };
    let SurfaceGeometry::Cylinder(cyl) = &curved.faces[ruling.wall].surface else {
        return Ok(None);
    };
    let cyl = cyl.cylinder();
    let axis = cyl.axis();
    let d = axis.direction.vector();
    let o = axis.location + d * d.dot(at - axis.location);
    let big = cyl.radius();
    let r = ruling.radius;
    // Convex: the cylinder's outward normal at the edge points away from
    // its axis (material inside), and the plane's from the material.
    if (at - o).dot(nc) <= 0.0 || r >= big {
        return Ok(None);
    }
    // The centre a radius behind the plane, `big - r` from the axis: along
    // the plane's trace in the section, the root nearer the edge.
    let e2 = ruling.u.cross(np);
    let base = at - np * r;
    let w = base - o;
    let (qb, qc) = (2.0 * w.dot(e2), w.dot(w) - (big - r) * (big - r));
    let disc = qb.mul_add(qb, -4.0 * qc);
    if disc < 0.0 {
        return Ok(None);
    }
    let root = disc.sqrt();
    let t = [(-qb - root) / 2.0, (-qb + root) / 2.0]
        .into_iter()
        .min_by(|x, y| x.abs().total_cmp(&y.abs()))
        .unwrap_or(0.0);
    let centre = base + e2 * t;
    let on_plane = centre + np * r;
    let radial = centre - o;
    let on_wall = o + radial * (big / radial.magnitude());
    // The contacts on the corner's sides of the edge: the plane's contact
    // inside the cylinder, the cylinder's behind the plane.
    if (on_plane - o).magnitude() > big + tol.confusion()
        || (on_wall - at).dot(np) > tol.confusion()
    {
        return Ok(None);
    }
    Ok(Some([at, centre, on_plane, on_wall]))
}

/// Whether one of the ruling's ends meets an earlier round's band (a
/// cylinder, sphere or torus meeting two faces or more tangentially along
/// its edges) tangent there to one of the ruling's faces.
fn at_earlier_band(curved: &Curved, ruling: &Ruling, tol: Tolerances) -> OgeomResult<bool> {
    for v in [ruling.ends.0, ruling.ends.1] {
        let at = curved.points[&v];
        for (f, face) in curved.faces.iter().enumerate() {
            if f == ruling.plane
                || f == ruling.wall
                || !face.loops.iter().flatten().any(|(x, _)| {
                    let g = &curved.edges[x];
                    g.ends.0 == v || g.ends.1 == v
                })
            {
                continue;
            }
            if !matches!(
                face.surface,
                SurfaceGeometry::Cylinder(_)
                    | SurfaceGeometry::Torus(_)
                    | SurfaceGeometry::Sphere(_)
            ) {
                continue;
            }
            let Some(n) = outward(curved, f, at) else {
                continue;
            };
            let touches = [ruling.plane, ruling.wall].iter().any(|g| {
                outward(curved, *g, at)
                    .is_some_and(|m| m.cross(n).magnitude() <= TANGENT && m.dot(n) > 0.0)
            });
            if touches && tangent_seats(curved, f, tol)? >= 2 {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// How many faces face `f` meets tangentially along its edges.
fn tangent_seats(curved: &Curved, f: usize, tol: Tolerances) -> OgeomResult<usize> {
    let mut seats: Vec<usize> = Vec::new();
    for (x, _) in curved.faces[f].loops.iter().flatten() {
        let edge = &curved.edges[x];
        let mid = edge
            .curve
            .point_at(0.5 * (edge.range.0 + edge.range.1), tol)?;
        for g in faces_of(curved, *x) {
            if g == f || seats.contains(&g) {
                continue;
            }
            if let (Some(a), Some(b)) = (outward(curved, f, mid), outward(curved, g, mid))
                && a.cross(b).magnitude() <= tol.angular() * 100.0
                && a.dot(b) > 0.0
            {
                seats.push(g);
            }
        }
    }
    Ok(seats.len())
}

/// Whether the ruling's ball sits on both its faces along it: each
/// contact, and the way to it from the edge across its face, on the solid's
/// boundary, or in air where an earlier band cut the face back (amendment
/// "a crease's seat past an earlier blend").
fn seated(
    model: &Model,
    solid: &Shape,
    curved: &Curved,
    ruling: &Ruling,
    tol: Tolerances,
) -> OgeomResult<bool> {
    for i in 1..ALONG {
        let at = ruling.a + ruling.u * (ruling.length * f64::from(i) / f64::from(ALONG));
        let Some([p, _, on_plane, on_wall]) = ball(curved, ruling, at, tol)? else {
            {
                return Ok(false);
            }
        };
        for (contact, face) in [(on_plane, ruling.plane), (on_wall, ruling.wall)] {
            for j in 1..=ACROSS {
                let mut q = p + (contact - p) * (f64::from(j) / f64::from(ACROSS));
                // On the cylinder, the way to the contact over its surface.
                if face == ruling.wall
                    && let SurfaceGeometry::Cylinder(cyl) = &curved.faces[face].surface
                {
                    let cyl = cyl.cylinder();
                    let axis = cyl.axis();
                    let d = axis.direction.vector();
                    let o = axis.location + d * d.dot(q - axis.location);
                    q = o + (q - o) * (cyl.radius() / (q - o).magnitude());
                }
                let held = ogeom_algo::classify_in_solid_exact(model, solid, q, tol)?;
                if held == ogeom_algo::Containment::On {
                    continue;
                }
                if held == ogeom_algo::Containment::Out
                    && let Some(n) = outward(curved, face, q)
                    && cut_back(model, solid, curved, (q, n), 4.0 * ruling.radius, tol)?
                {
                    continue;
                }
                {
                    return Ok(false);
                }
            }
        }
    }
    Ok(true)
}

/// Whether `p`, a point of a face's surface (outward normal `n`) in air
/// beside the ruling, stands where an earlier band cut that face back:
/// straight in from it the solid is first met on such a band (a face
/// meeting two others or more tangentially), within `reach`, not over a
/// groove.
pub(crate) fn cut_back(
    model: &Model,
    solid: &Shape,
    curved: &Curved,
    (p, n): (Point, Vector),
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    use ogeom_algo::Containment;
    let steps = 256_u32;
    let mut air = 0.0;
    for k in 1..=steps {
        let t = reach * f64::from(k) / f64::from(steps);
        let Ok(at) = ogeom_algo::classify_in_solid_exact(model, solid, p - n * t, tol) else {
            return Ok(false);
        };
        if at != Containment::Out {
            break;
        }
        air = t;
        if k == steps {
            return Ok(false);
        }
    }
    // The first point met, halved to rounding.
    let (mut lo, mut hi) = (air, air + reach / f64::from(steps));
    for _ in 0..50 {
        let mid = f64::midpoint(lo, hi);
        let Ok(at) = ogeom_algo::classify_in_solid_exact(model, solid, p - n * mid, tol) else {
            return Ok(false);
        };
        if at == Containment::Out {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let met = p - n * hi;
    let deflection = ogeom_mesh::Deflection::default();
    for (f, face) in curved.faces.iter().enumerate() {
        if !matches!(
            face.surface,
            SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Torus(_) | SurfaceGeometry::Sphere(_)
        ) || tangent_seats(curved, f, tol)? < 2
        {
            continue;
        }
        if ogeom_algo::classify_on_face(model, &face.shape, met, deflection, tol)
            .is_ok_and(|c| c != Containment::Out)
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The ruling's spandrel as a prism `reach` past each end: its section the
/// corner between the edge and the ball's arc, standing a quarter radius off
/// the plane and the cylinder past the contacts so no tool face lies on the
/// solid's own past them; with its band face (the ball's arc swept).
fn prism(
    model: &mut Model,
    curved: &Curved,
    ruling: &Ruling,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Shape)>> {
    let start = ruling.a - ruling.u * reach;
    let Some([p, centre, on_plane, on_wall]) = ball(curved, ruling, start, tol)? else {
        {
            return Ok(None);
        }
    };
    let Some(nc) = outward(curved, ruling.wall, ruling.a) else {
        {
            return Ok(None);
        }
    };
    let r = ruling.radius;
    let off = r * 0.25;
    let past = p + (p - on_plane) / p.distance(on_plane) * off;
    let out = on_wall + nc * off;
    let [vp, vq, vo, vw] =
        [on_plane, past, out, on_wall].map(|q| ogeom_algo::make_vertex(model, q).shape);
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
    let back = line(on_wall, out, &vw, &vo)?;
    let across = line(out, past, &vo, &vq)?;
    let along = line(past, on_plane, &vq, &vp)?;
    let (x, y) = (on_plane - centre, on_wall - centre);
    let normal = x.cross(y);
    if normal.magnitude() <= tol.confusion() * r {
        {
            return Ok(None);
        }
    }
    let frame = Frame::new(
        centre,
        Direction::new(normal, tol)?,
        Direction::new(x, tol)?,
        tol,
    )?;
    let arc_curve: Curve = CircleCurve::new(ogeom_math::Circle::new(frame, r, tol)?).into();
    let arc =
        ogeom_algo::make_edge_between(model, arc_curve, (0.0, x.angle(y, tol)?), &vp, &vw, tol)?
            .shape;
    let wire = ogeom_algo::make_wire(model, &[arc.clone(), back, across, along], tol)?.shape;
    let plane = PlaneSurface::new(ogeom_math::Plane::through(
        start,
        Direction::new(ruling.u, tol)?,
    ));
    let face = ogeom_algo::make_face(model, plane.into(), &[wire], tol)?.shape;
    let swept =
        ogeom_algo::make_prism(model, &face, ruling.u * (ruling.length + 2.0 * reach), tol)?;
    let Some(band) = swept.history.generated(&arc).first().cloned() else {
        {
            return Ok(None);
        }
    };
    let solids = explore_unique(model, &swept.shape, ShapeType::Solid)?;
    let [solid] = solids.as_slice() else {
        {
            return Ok(None);
        }
    };
    Ok(Some((solid.clone(), band)))
}
