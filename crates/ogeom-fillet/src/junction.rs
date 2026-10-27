//! Two blends of one request meeting at a vertex where their creases turn.
//!
//! Two convex creases of one `fillet_edges` call that share a vertex and a
//! host face, and leave the vertex other than tangentially, are two bands
//! meeting at a corner. Each band carried on past the vertex is the prism
//! of its end section along its crease's end tangent: a straight band of
//! the same radius, whose ball centres run along a line a radius in from
//! both its hosts' tangent planes. Two such lines lie a radius below the
//! shared host, so they meet, and two straight bands of one radius whose
//! centre lines meet cross exactly in the plane through that point
//! bisecting the lines.
//!
//! Where the turn is re-entrant (the one crease carried on runs into the
//! material behind the other's second face), the solid goes on past the
//! vertex beside both bands, and the corner both rounds take is the cell
//! past both caps that lies in both continuations: each continuation up to
//! that bisecting plane. Each band is carried on past its cap to the plane,
//! and the two are joined across it into one wedge.
//!
//! Where the turn is convex (the crease carried on leaves the solid through
//! the other's second face), each band's continuation lies inside the
//! other's wedge, and the two wedges as they stand remove their union; the
//! bands themselves meet along their own intersection. Where the two bands'
//! sections in the bisecting plane agree to the fit a marched band is held
//! to, the turn is shallow enough that the bands continue each other: their
//! intersection is that plane to the fit, and the two nearly coincident
//! band surfaces are not asked to cross. Each band stops in that plane, and
//! the two are joined across it into one wedge.
//!
//! The plane is the two bands' meeting only where the reflection in it
//! carries one band's end section onto the other's: the shared host a
//! plane, the reflection taking the one wall's tangent plane at the vertex
//! onto the other's, and the plane through the vertex, each to the fit a
//! marched band is held to.
//!
//! Where the sections differ (a drafted wall beside a square one) and both
//! creases are straight, the junction is crossed: no plane bounds it. Across
//! a re-entrant turn the corner takes what lies above both bands carried
//! on, what stands past a band's wall's plane lying above that band as much
//! as what its ball rolls off, so each band is carried on with its section
//! widened past its wall and the two are intersected in one boolean, which
//! bounds the corner by the bands' own intersection ([`corner`]). Across a
//! convex turn each band runs out past the vertex as it does alone, and
//! the two cuts take their union; so too where one crease turns about an
//! axis along the shared host's normal on a drum (a D-shaped block's arc
//! beside its chord), whose band's curvature keeps the two from meeting in
//! a plane.
//!
//! One crease of a re-entrant crossed junction may turn about an axis along
//! the shared host's normal, on a drum about that axis, beside a straight
//! one, the material within the drum (a boss) or outside it (a bore). Its
//! band is a torus band whose section is the same at every turn, and
//! before its cap what lies above it is its band or what stands off the
//! drum. The corner is that band's part from its cap to where its section
//! stands wholly past the other's wall, widened by what stands off the
//! drum there and intersected as above, less the open corner off the drum
//! past the other's wall; the rest of the band is joined back on. Where
//! the straight band carried on stands off the drum before the cap beyond
//! that part, no junction is taken, nor elsewhere (a curved shared host, a
//! curved crease whose sections differ otherwise): the bands meet as they
//! did without one.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Orientation, Shape, ShapeType, explore_unique};

use crate::fillet::Mate;

/// Two edges of a request meeting at a vertex where their creases turn.
#[derive(Clone)]
pub(crate) struct Junction {
    /// Where the creases meet: their shared vertex, or where their end
    /// tangents meet within that vertex's tolerance of it.
    pub at: Point,
    /// The two edges, by their index in the request.
    pub pair: [usize; 2],
    /// Each crease's direction past the vertex, away from its own edge.
    pub outward: [Vector; 2],
    /// The turn runs into the material rather than out of it.
    pub reentrant: bool,
    /// The plane the bands meet flush in, through a point on it and by its
    /// unit normal towards the first edge's side.
    pub flush: (Point, Vector),
    /// The bands' end sections are no reflection of each other, so no
    /// plane bounds the meeting: across a re-entrant turn the boolean's own
    /// intersection of the bands carried on bounds the corner, and across
    /// a convex one each band runs out as it does alone.
    pub crossed: bool,
    /// The shared host's outward unit normal at the vertex.
    pub top: Vector,
    /// Each edge's wall's outward unit normal at the vertex.
    pub sides: [Vector; 2],
    /// Of a crossed re-entrant junction whose one crease is a circle about
    /// an axis along `top`, on the cylinder through it about that axis, and
    /// whose other is straight: that band, a torus band, its section the
    /// same at every turn about the axis.
    pub turned: Option<Turned>,
}

/// The turned band of a crossed junction: its crease a circle about an
/// axis along the shared host's normal, on the cylinder through it.
#[derive(Clone, Copy)]
pub(crate) struct Turned {
    /// The band's slot in the junction.
    pub slot: usize,
    /// The circle's centre.
    pub centre: Point,
    /// The circle's radius, the cylinder's.
    pub big: f64,
    /// The material lies within the cylinder (a boss), not outside it (a
    /// bore).
    pub within: bool,
}

impl Junction {
    /// The slot `index` holds in this junction, if it is one of the pair.
    pub(crate) fn slot(&self, index: usize) -> Option<usize> {
        self.pair.iter().position(|&i| i == index)
    }

    /// The bisecting plane as the half-space the band in `slot` keeps: a
    /// point on it and the outward normal of the side it keeps.
    pub(crate) fn keeps(&self, slot: usize) -> (Point, Vector) {
        let (on, normal) = self.flush;
        (on, if slot == 0 { normal } else { -normal })
    }
}

/// The junctions among `edges` of `solid`, given the faces each lies
/// between, the request's mates, and the radius: every pair of edges
/// sharing an end vertex no third edge of the request ends at (a corner the
/// corner tool speaks to), each leaving it along its mate's direction at
/// its end nearer the vertex, where the meeting is exact. A pair whose
/// junction cannot be read is no junction.
pub(crate) fn junctions(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    hosts: &[Vec<Shape>],
    mates: &[Mate],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<Junction>> {
    let ends: Vec<Vec<Shape>> = edges
        .iter()
        .map(|e| -> OgeomResult<Vec<Shape>> {
            Ok(ogeom_algo::edge_vertices(model, e)?
                .map(|(a, b)| if a.is_same(&b) { vec![a] } else { vec![a, b] })
                .unwrap_or_default())
        })
        .collect::<OgeomResult<_>>()?;
    let point_of = |v: &Shape| -> Option<Point> {
        let data = model.node(v)?.data().as_vertex()?;
        Some(v.transform(model.datums()).ok()?.apply(data.point))
    };
    // Each edge's own end nearest the vertex, and its direction leaving it.
    let end_at = |mate: &Mate, at: Point| -> (Point, Vector) {
        let [(p, a), (q, b)] = mate.ends;
        if p.distance(at) <= q.distance(at) {
            (p, a)
        } else {
            (q, b)
        }
    };
    // The faces holding each edge as it stands: a prism's top edge and its
    // bottom one may be one edge placed twice, whose faces the hosts list
    // holds together.
    let mut held_by: Vec<Vec<Shape>> = Vec::with_capacity(edges.len());
    for (edge, faces) in edges.iter().zip(hosts) {
        let mut own = Vec::new();
        for face in faces {
            let holds = explore_unique(model, face, ShapeType::Edge)?
                .iter()
                .any(|e| crate::support::same_occurrence(model, e, edge, tol));
            if holds {
                own.push(face.clone());
            }
        }
        held_by.push(own);
    }
    let hosts = held_by.as_slice();
    let mut found = Vec::new();
    for i in 0..edges.len() {
        for j in i + 1..edges.len() {
            for v in &ends[i] {
                if !ends[j].iter().any(|w| w.is_same(v)) {
                    continue;
                }
                if ends
                    .iter()
                    .filter(|e| e.iter().any(|w| w.is_same(v)))
                    .count()
                    > 2
                {
                    continue;
                }
                let Some(at) = point_of(v) else {
                    continue;
                };
                // The vertex may stand off the faces it bounds within its
                // own tolerance; their normals are read that far.
                let held = model
                    .tolerance_of(v)?
                    .map_or(tol.confusion(), ogeom_core::Tolerance::get)
                    .max(tol.confusion() * 1e3);
                let (end_i, leaving_i) = end_at(&mates[i], at);
                let (end_j, leaving_j) = end_at(&mates[j], at);
                // The creases meet where their end tangents do: the vertex,
                // or where a vertex stands off the walls within its
                // tolerance, the corner the edges' lines make.
                let at = closest_on_lines(end_i, leaving_i, end_j, leaving_j)
                    .filter(|meet| meet.distance(at) <= held)
                    .unwrap_or(at);
                if let Ok(Some(junction)) = junction(
                    model,
                    solid,
                    [i, j],
                    [&hosts[i], &hosts[j]],
                    [edges[i].clone(), edges[j].clone()],
                    (at, held),
                    [leaving_i, leaving_j],
                    radius,
                    tol,
                ) {
                    found.push(junction);
                }
            }
        }
    }
    Ok(found)
}

/// The junction of two edges of `solid` at `at`, each leaving it along its
/// own `leaving` direction, or `None` where they are no pair of convex
/// creases turning on one shared planar host whose bands meet exactly in
/// the plane bisecting them.
#[allow(clippy::too_many_arguments, reason = "one question, all its data")]
#[allow(clippy::too_many_lines, reason = "one junction, read end to end")]
fn junction(
    model: &Model,
    solid: &Shape,
    pair: [usize; 2],
    hosts: [&Vec<Shape>; 2],
    edges: [Shape; 2],
    (at, held): (Point, f64),
    leaving: [Vector; 2],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Junction>> {
    // A tangent pair is a chain, whose caps stand flush in one section.
    if leaving[0].cross(-leaving[1]).magnitude() <= 1e-2 && leaving[0].dot(-leaving[1]) > 0.0 {
        return Ok(None);
    }
    // The host the two creases share: one face, or two faces on one
    // surface (a top face a boolean split along a line through the vertex).
    let mut shared: Option<(&Shape, &Shape)> = None;
    for a in hosts[0] {
        for b in hosts[1] {
            if one_host(model, a, b, at, radius, tol)? {
                if shared.is_some() {
                    return Ok(None);
                }
                shared = Some((a, b));
            }
        }
    }
    let Some((top_0, top_1)) = shared else {
        return Ok(None);
    };
    if plane_of(model, top_0, tol)?.is_none() || plane_of(model, top_1, tol)?.is_none() {
        return Ok(None);
    }
    let other = |k: usize, top: &Shape| hosts[k].iter().find(|f| !f.is_same(top)).cloned();
    let (Some(wall_0), Some(wall_1)) = (other(0, top_0), other(1, top_1)) else {
        return Ok(None);
    };
    // Both creases convex: a fill in a re-entrant edge is material, and two
    // bands meet as described only where both cut. Read off the wall, the
    // smaller face: the shared host may hold a great many others' trims.
    for (k, edge) in edges.iter().enumerate() {
        let top = if k == 0 { top_0 } else { top_1 };
        let wall = if k == 0 { &wall_0 } else { &wall_1 };
        if !crease_is_convex(model, edge, [wall, top], radius, tol)? {
            return Ok(None);
        }
    }
    let outward = [-leaving[0], -leaving[1]];
    let normal_of = |face: &Shape| -> OgeomResult<Option<Vector>> {
        let Some(normal) = crate::marched::face_normal_within(model, face, at, held, tol)? else {
            return Ok(None);
        };
        let sign = if face.orientation() == Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
        Ok(Some(normal * sign))
    };
    let (Some(top), Some(side_0), Some(side_1)) =
        (normal_of(top_0)?, normal_of(&wall_0)?, normal_of(&wall_1)?)
    else {
        return Ok(None);
    };
    // Each crease carried on past the vertex against the other's second
    // face: out of the solid through it, or into the material behind it.
    // A crease a radius on stands that far off the face's tangent plane,
    // read against the face's own tolerance.
    let lean = |face: &Shape, normal: Vector, along: Vector| -> OgeomResult<Option<f64>> {
        let held = model
            .tolerance_of(face)?
            .map_or(tol.confusion(), ogeom_core::Tolerance::get);
        let off = normal.dot(along) * radius;
        Ok((off.abs() > held).then_some(off))
    };
    let (Some(lean_0), Some(lean_1)) = (
        lean(&wall_1, side_1, outward[0])?,
        lean(&wall_0, side_0, outward[1])?,
    ) else {
        return Ok(None);
    };
    if (lean_0 < 0.0) != (lean_1 < 0.0) {
        return Ok(None);
    }
    let reentrant = lean_0 < 0.0;
    // Each band's ball centre at the vertex, a radius in from the tangent
    // planes of its two hosts there, and the point where the two centre
    // lines, along the creases' end tangents, meet.
    let centre = |side: Vector| -> Point { at - (top + side) * (radius / (1.0 + top.dot(side))) };
    // The plane they meet in holds the shared host's normal, so only the
    // lines' shadows on that host place it: taken there, two creases at a
    // shallow turn do not turn a hair of skew between the lines into a
    // shift along them.
    let flat = |v: Vector| v - top * v.dot(top);
    let (c0, c1) = (centre(side_0), centre(side_1));
    let (c0, c1) = (at + flat(c0 - at), at + flat(c1 - at));
    let (u0, u1) = (flat(outward[0]), flat(outward[1]));
    let (m0, m1) = (u0.magnitude(), u1.magnitude());
    if m0 <= tol.confusion() || m1 <= tol.confusion() {
        return Ok(None);
    }
    let Some(meet) = closest_on_lines(c0, u0 / m0, c1, u1 / m1) else {
        return Ok(None);
    };
    let bisector = {
        let u = leaving[0] - leaving[1];
        let u = u - top * u.dot(top);
        let m = u.magnitude();
        if m <= tol.confusion() {
            return Ok(None);
        }
        u / m
    };
    // The plane is the bands' meeting where the reflection in it carries
    // the one band's end section onto the other's: the shared host's normal
    // lies in the plane, so the host goes to itself, and the one wall's
    // tangent plane must go to the other's. A section point a radius from
    // the vertex moves by the radius times the walls' normals' miss, and by
    // twice the vertex's distance off the plane.
    let reflect = |v: Vector| v - bisector * (2.0 * v.dot(bisector));
    let miss = radius.mul_add(
        (reflect(side_0) - side_1).magnitude(),
        2.0 * (at - meet).dot(bisector).abs(),
    );
    let fit = crate::marched::band_fit_target(tol);
    let turn = outward[0].dot(leaving[1]).clamp(-1.0, 1.0).acos();
    // Where the end sections are no reflection of each other, two straight
    // creases still meet exactly: each band carried on is the straight band
    // of its end section. Across a re-entrant turn the boolean bounds the
    // corner both take along the bands' own intersection, whether the
    // material holds it asked of the corner itself; across a convex one
    // each band runs out past the vertex as it does alone, built from the
    // solid as it stands, and the two cuts take their union.
    let straight = || -> OgeomResult<bool> {
        for edge in &edges {
            let (curve, _) = crate::support::edge_curve(model, edge, tol)?;
            if !matches!(curve, ogeom_geom::Curve::Line(_)) {
                return Ok(false);
            }
        }
        Ok(true)
    };
    // One crease a circle about an axis along the shared host's normal, on
    // the cylinder through it, the other straight.
    let turned_beside_straight = || -> OgeomResult<bool> {
        for (k, edge) in edges.iter().enumerate() {
            let wall = if k == 0 { &wall_0 } else { &wall_1 };
            let (other, _) = crate::support::edge_curve(model, &edges[1 - k], tol)?;
            if matches!(other, ogeom_geom::Curve::Line(_))
                && turned_about(model, edge, wall, top, tol)?.is_some()
            {
                return Ok(true);
            }
        }
        Ok(false)
    };
    let mut turned = None;
    let crossed = if reentrant {
        if miss > fit {
            if !straight()? {
                // One band turning about an axis along the shared host's
                // normal beside a straight one: the turned band's section
                // is its end section at every turn, and the corner is built
                // from that band as it stands before its cap.
                for (k, edge) in edges.iter().enumerate() {
                    let wall = if k == 0 { &wall_0 } else { &wall_1 };
                    let (other, _) = crate::support::edge_curve(model, &edges[1 - k], tol)?;
                    if matches!(other, ogeom_geom::Curve::Line(_)) {
                        turned = turned_about(model, edge, wall, top, tol)?.map(|(centre, big)| {
                            // The wall's outward normal points off the axis
                            // where the material lies within the cylinder.
                            let from_axis = at - centre;
                            let side = if k == 0 { side_0 } else { side_1 };
                            Turned {
                                slot: k,
                                centre,
                                big,
                                within: side.dot(from_axis - top * from_axis.dot(top)) > 0.0,
                            }
                        });
                    }
                }
                if turned.is_none() {
                    return Ok(None);
                }
            }
            true
        } else {
            // Carried on, each band is a straight band of its end section,
            // so the reflection alone decides the meeting. The cell it
            // closes reaches from the vertex to where the bands' contacts on
            // the shared host meet, a radius under it: the material must
            // hold that reach, or the continuations run out of the solid (a
            // turn folding back on itself) and the cell is not the corner's.
            let probe = crate::support::material_probe(model, solid, tol)?;
            // The vertex itself stands on the solid's edge; its
            // neighbourhood past it is what the cell takes.
            let under = top * (radius * 0.5);
            for k in 1..=8 {
                let p = at + (meet - at) * (f64::from(k) / 8.0) - under;
                if probe.holds(model, p, tol)? != ogeom_algo::Containment::In {
                    return Ok(None);
                }
            }
            false
        }
    } else {
        // Across a convex turn the bands themselves meet. Each sweeps its
        // section along its centre line; where that line curves, the
        // section a distance s along stands kappa s^2 / 2 off the straight
        // continuation. In the bisecting plane a section point a distance d
        // from the vertex's edge stands d tan(turn/2) along, and d reaches
        // the radius.
        let mut bend = 0.0;
        for edge in &edges {
            let kappa = crease_curvature(model, edge, at, tol)?;
            if kappa * radius >= 1.0 {
                return Ok(None);
            }
            bend += kappa / (1.0 - kappa * radius);
        }
        let along = radius * (turn / 2.0).tan();
        if (bend * along).mul_add(along / 2.0, miss) > fit {
            if !straight()? && !turned_beside_straight()? {
                return Ok(None);
            }
            true
        } else {
            false
        }
    };
    Ok(Some(Junction {
        at,
        pair,
        outward,
        reentrant,
        flush: (meet, bisector),
        crossed,
        top,
        sides: [side_0, side_1],
        turned,
    }))
}

/// Where `edge` is a circle about an axis along the unit `top`, on a `wall`
/// that is the cylinder through it about that axis, each within the wall's
/// own tolerance: the circle's centre and radius.
fn turned_about(
    model: &Model,
    edge: &Shape,
    wall: &Shape,
    top: Vector,
    tol: Tolerances,
) -> OgeomResult<Option<(Point, f64)>> {
    use ogeom_geom::Transformable as _;
    let (curve, _) = crate::support::edge_curve(model, edge, tol)?;
    let ogeom_geom::Curve::Circle(circle) = curve else {
        return Ok(None);
    };
    let circle = circle.circle();
    if circle.frame().z().vector().cross(top).magnitude() > tol.angular() {
        return Ok(None);
    }
    let Some(surface) = model
        .node(wall)
        .and_then(|n| n.data().as_face())
        .and_then(|d| model.geometry().surface(d.surface))
        .cloned()
    else {
        return Ok(None);
    };
    let ogeom_geom::SurfaceGeometry::Cylinder(placed) =
        surface.transformed(&wall.transform(model.datums())?, tol)?
    else {
        return Ok(None);
    };
    let cylinder = placed.cylinder();
    let axis = cylinder.axis();
    let along = axis.direction.vector();
    if along.cross(top).magnitude() > tol.angular() {
        return Ok(None);
    }
    let held = model
        .tolerance_of(wall)?
        .map_or(tol.confusion(), ogeom_core::Tolerance::get)
        .max(tol.confusion());
    let off = circle.centre() - axis.location;
    if (off - along * off.dot(along)).magnitude() > held
        || (cylinder.radius() - circle.radius()).abs() > held
    {
        return Ok(None);
    }
    Ok(Some((circle.centre(), circle.radius())))
}

/// Whether `edge` between `hosts` is a convex crease: the material lies
/// within the angle the faces make there.
fn crease_is_convex(
    model: &Model,
    edge: &Shape,
    hosts: [&Shape; 2],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    use ogeom_geom::Transformable as _;
    let (guide, range) = crate::support::edge_curve(model, edge, tol)?;
    let mut placed = Vec::with_capacity(2);
    for face in hosts {
        let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
            ogeom_core::ogeom_bail!(Dangling, "face is not in this model");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            ogeom_core::ogeom_bail!(Dangling, "face refers to a surface not in this model");
        };
        let sign = if face.orientation() == Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
        placed.push((
            surface.transformed(&face.transform(model.datums())?, tol)?,
            sign,
        ));
    }
    crate::marched::crease_convexity(
        model,
        hosts[0],
        [(&placed[0].0, placed[0].1), (&placed[1].0, placed[1].1)],
        &guide,
        range,
        radius,
        tol,
    )
}

/// The point midway between the nearest points of two lines, each through
/// a point along a unit direction; `None` for parallel lines.
fn closest_on_lines(p: Point, u: Vector, q: Point, v: Vector) -> Option<Point> {
    let w = p - q;
    let b = u.dot(v);
    let denominator = 1.0 - b * b;
    if denominator <= 1e-12 {
        return None;
    }
    let (d, e) = (u.dot(w), v.dot(w));
    let s = (b * e - d) / denominator;
    let t = (e - b * d) / denominator;
    Some((p + u * s).midpoint(q + v * t))
}

/// Whether faces `a` and `b` are one host where they meet `at`: the same
/// face, two faces on one surface, or two planes that agree, within the
/// faces' own tolerances, over the `radius` round `at` a band spans.
fn one_host(
    model: &Model,
    a: &Shape,
    b: &Shape,
    at: Point,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    if a.is_same(b) {
        return Ok(true);
    }
    let surface_of = |f: &Shape| {
        model
            .node(f)
            .and_then(|n| n.data().as_face())
            .map(|d| d.surface)
    };
    let (Some(sa), Some(sb)) = (surface_of(a), surface_of(b)) else {
        return Ok(false);
    };
    if sa == sb && a.location() == b.location() {
        return Ok(true);
    }
    let (Some((pa, na)), Some((pb, nb))) = (plane_of(model, a, tol)?, plane_of(model, b, tol)?)
    else {
        return Ok(false);
    };
    let held = [a, b]
        .iter()
        .map(|f| {
            model
                .tolerance_of(f)
                .map(|t| t.map_or(tol.confusion(), ogeom_core::Tolerance::get))
        })
        .collect::<OgeomResult<Vec<f64>>>()?
        .into_iter()
        .fold(tol.confusion(), f64::max);
    // Two in-plane directions round `at` on `a`, a radius out each way.
    let (u, v) = {
        let seed = if na.x.abs() < 0.9 {
            Vector::new(1.0, 0.0, 0.0)
        } else {
            Vector::new(0.0, 1.0, 0.0)
        };
        let u = seed - na * seed.dot(na);
        let u = u / u.magnitude();
        (u, na.cross(u))
    };
    let on_a = |p: Point| p - na * (p - pa).dot(na);
    for p in [
        at,
        at + u * radius,
        at - u * radius,
        at + v * radius,
        at - v * radius,
    ] {
        if (on_a(p) - pb).dot(nb).abs() > held {
            return Ok(false);
        }
    }
    // And facing the same way: the material on one side of both.
    let facing = |f: &Shape| {
        if f.orientation() == Orientation::Reversed {
            -1.0
        } else {
            1.0
        }
    };
    Ok(na.dot(nb) * facing(a) * facing(b) > 0.0)
}

/// The curvature of an edge's curve at its end nearest `at`.
fn crease_curvature(model: &Model, edge: &Shape, at: Point, tol: Tolerances) -> OgeomResult<f64> {
    let (curve, (lo, hi)) = crate::support::edge_curve(model, edge, tol)?;
    let t = if curve.point_at(lo, tol)?.distance(at) <= curve.point_at(hi, tol)?.distance(at) {
        lo
    } else {
        hi
    };
    let h = (hi - lo) * 1e-5;
    let (a, b) = if t == lo { (t, t + h) } else { (t - h, t) };
    let d1 = curve.d1_at(t, tol)?;
    let d2 = (curve.d1_at(b, tol)? - curve.d1_at(a, tol)?) / (b - a);
    let speed = d1.magnitude();
    if speed <= tol.confusion() {
        return Ok(0.0);
    }
    Ok(d1.cross(d2).magnitude() / speed.powi(3))
}

/// The planar face of `wedge` that caps it at `at`, square to the crease
/// leaving there along `outward`: the planar face through `at` whose normal
/// lies nearest `outward`.
pub(crate) fn cap_at(
    model: &Model,
    wedge: &Shape,
    at: Point,
    outward: Vector,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let near = crate::marched::band_fit_target(tol);
    let mut best: Option<(f64, Shape)> = None;
    for face in explore_unique(model, wedge, ShapeType::Face)? {
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
            continue;
        };
        let Some(ogeom_geom::SurfaceGeometry::Plane(plane)) =
            model.geometry().surface(data.surface)
        else {
            continue;
        };
        use ogeom_geom::Transformable as _;
        let placed = ogeom_geom::SurfaceGeometry::Plane(*plane)
            .transformed(&face.transform(model.datums())?, tol)?;
        let ogeom_geom::SurfaceGeometry::Plane(placed) = placed else {
            continue;
        };
        let plane = placed.plane();
        if plane.signed_distance_to(at).abs() > near {
            continue;
        }
        // A band along a circle has both caps in one plane through its axis
        // where it turns half way round: the cap is the one reaching `at`.
        if ogeom_algo::shape_bounds(model, &face, tol)?.distance_to(at) > near {
            continue;
        }
        let along = plane.normal().vector().dot(outward).abs();
        if best.as_ref().is_none_or(|(b, _)| along > *b) {
            best = Some((along, face));
        }
    }
    Ok(best.filter(|(along, _)| *along > 0.5).map(|(_, f)| f))
}

/// Two wedges each ending in the plane through `on` with unit `normal`
/// at the junction's vertex `at`, joined across it into one: every face of
/// either but the two in that plane at `at`, sewn across `gap`, the most
/// their sections there part by. With it, each face of either wedge kept
/// and the face it became.
pub(crate) fn joined(
    model: &mut Model,
    wedges: [&Shape; 2],
    (on, normal): (Point, Vector),
    at: Point,
    gap: f64,
    tol: Tolerances,
) -> OgeomResult<(Shape, Vec<(Shape, Shape)>)> {
    // The sewing keeps the first of two twin edges. A fitted band is held
    // only to its fit and its edges are widened to say so, while a face on
    // an exact surface (a torus beside it, say, met again by its own rim's
    // next piece) reads its edges in its own chart: the wedge with fewer
    // fitted faces goes first, so its edges are the ones kept.
    let fitted = |model: &Model, wedge: &Shape| -> OgeomResult<usize> {
        Ok(explore_unique(model, wedge, ShapeType::Face)?
            .iter()
            .filter(|f| {
                model
                    .node(f)
                    .and_then(|n| n.data().as_face())
                    .and_then(|d| model.geometry().surface(d.surface))
                    .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::BSpline(_)))
            })
            .count())
    };
    let wedges = if fitted(model, wedges[1])? < fitted(model, wedges[0])? {
        [wedges[1], wedges[0]]
    } else {
        wedges
    };
    let mut faces = Vec::new();
    // A loop closing on itself is one wedge, both its ends in the plane.
    let count = if wedges[0].is_same(wedges[1]) { 1 } else { 2 };
    for wedge in &wedges[..count] {
        for face in explore_unique(model, wedge, ShapeType::Face)? {
            if !end_face(model, &face, (on, normal), at, gap, tol)? {
                faces.push(face);
            }
        }
    }
    let sewn = ogeom_algo::sew_within(model, &faces, gap, tol)?;
    let traced: Vec<(Shape, Shape)> = faces
        .iter()
        .flat_map(|f| {
            sewn.history
                .trace(f)
                .iter()
                .map(|image| (f.clone(), image.clone()))
                .collect::<Vec<_>>()
        })
        .collect();
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_core::ogeom_bail!(
            Construction,
            "two bands meeting flush at a shallow turn did not join into one wedge"
        );
    }
    for face in explore_unique(model, &sewn.shells[0], ShapeType::Face)? {
        trimmed_on_its_surface(model, &face, gap, tol)?;
    }
    ogeom_algo::state_pcurve_gaps(model, &sewn.shells[0], tol)?;
    Ok((
        model.add_solid(std::slice::from_ref(&sewn.shells[0]))?,
        traced,
    ))
}

/// Whether `wedge` ends in the plane through `on` with unit `normal` at
/// the junction's point `at`: one of its faces lies in that plane there.
pub(crate) fn ends_in(
    model: &Model,
    wedge: &Shape,
    plane: (Point, Vector),
    at: Point,
    gap: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    for face in explore_unique(model, wedge, ShapeType::Face)? {
        if end_face(model, &face, plane, at, gap, tol)? {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether `face` is a wedge's end in the plane through `on` with unit
/// `normal` at `at`: a plane lying in it within `gap`, and reaching `at`.
/// Two corners of one solid may share a bisecting plane (opposite corners
/// of a square), so the face must stand at this one.
fn end_face(
    model: &Model,
    face: &Shape,
    (on, normal): (Point, Vector),
    at: Point,
    gap: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let here = ogeom_algo::shape_bounds(model, face, tol)?.distance_to(at) <= gap;
    Ok(here && in_plane(model, face, on, normal, gap)?)
}

/// Whether `face` is a plane lying in the plane through `at` with unit
/// `normal`, within `gap`.
fn in_plane(model: &Model, face: &Shape, at: Point, normal: Vector, gap: f64) -> OgeomResult<bool> {
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
        return Ok(false);
    };
    if !matches!(
        model.geometry().surface(data.surface),
        Some(ogeom_geom::SurfaceGeometry::Plane(_))
    ) {
        return Ok(false);
    }
    for vertex in explore_unique(model, face, ShapeType::Vertex)? {
        let Some(v) = model.node(&vertex).and_then(|n| n.data().as_vertex()) else {
            return Ok(false);
        };
        let p = vertex.transform(model.datums())?.apply(v.point);
        if (p - at).dot(normal).abs() > gap {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The two bands of a crossed junction and the corner they take across
/// their re-entrant turn, in one wedge, with the faces of it that came from
/// the second band.
///
/// Carried on past the vertex, each band is the straight band of its end
/// section, and what it leaves is what lies under both its hosts' planes
/// outside its ball: what stands past its wall's plane lies above the band
/// as much as what its ball rolls off. Near the corner the material taken
/// is what lies above both bands, bounded by the boolean's own intersection
/// of the two. So each band is taken as its wedge carried on through the
/// corner, its section widened past its wall over the reach where the
/// other's wedge stands, from as far back into its own wedge as the other's
/// wedge, or the other band carried on past this one's wall, reaches, and
/// the two are intersected in one boolean: each band's wedge is kept to
/// what stands above the other carried on, and past both caps the two
/// continuations keep the cell above both.
///
/// What stands past both walls' planes is the open corner, taken out of
/// the result. The corner is taken only where the material holds all of
/// it but the wedges' own run-outs (a turn folding back past the solid's
/// far side leaves it), and where it holds all that each band carried on
/// keeps past the other's wall, which lies above both. Elsewhere this
/// refuses, and the request rounds without crossed junctions.
pub(crate) fn corner(
    model: &mut Model,
    wedges: [&Shape; 2],
    junction: &Junction,
    solid: &Shape,
    tol: Tolerances,
) -> OgeomResult<(Shape, Vec<Shape>)> {
    let at = junction.at;
    let near = tol.confusion() * 1e3;
    // A turned band takes part from its cap back to where its section
    // stands wholly past the other's wall; the rest stands as built.
    let mut parts = [wedges[0].clone(), wedges[1].clone()];
    let mut rest = None;
    let mut sector = None;
    if let Some(turned) = junction.turned {
        let slot = turned.slot;
        let (part, far, plane, quadrant) =
            turned_split(model, wedges[slot], junction, slot, turned.centre, tol)?;
        parts[slot] = part;
        rest = Some((slot, far, plane));
        sector = Some(quadrant);
    }
    let wedges = [&parts[0], &parts[1]];
    let mut caps = Vec::with_capacity(2);
    let mut sections = Vec::with_capacity(2);
    for (slot, wedge) in wedges.iter().enumerate() {
        let Some(cap) = cap_at(model, wedge, at, junction.outward[slot], tol)? else {
            ogeom_core::ogeom_bail!(
                Construction,
                "a band meeting another across a re-entrant turn has no cap at the vertex"
            );
        };
        sections.push(Section::of(model, &cap, near, tol)?);
        caps.push(cap);
    }
    // Each wedge is the straight band of its end section from its cap back
    // to its far end, trimmed there by planes or run out, its faces planes
    // and its band a cylinder. Along the crease every face reaches no
    // farther than its edges do, each within the hull of its control
    // points, so the wedge lies within the hull of its section's hull
    // points placed at the cap and as far back as its edges reach. A
    // turned band's part is a torus band between a plane and a cylinder,
    // within the corners of its own guaranteed bound.
    let mut hulls: Vec<Vec<Point>> = Vec::with_capacity(2);
    for (slot, wedge) in wedges.iter().enumerate() {
        if junction.turned.is_some_and(|t| t.slot == slot) {
            for face in explore_unique(model, wedge, ShapeType::Face)? {
                let surface = model
                    .node(&face)
                    .and_then(|n| n.data().as_face())
                    .and_then(|d| model.geometry().surface(d.surface));
                if !matches!(
                    surface,
                    Some(
                        ogeom_geom::SurfaceGeometry::Plane(_)
                            | ogeom_geom::SurfaceGeometry::Cylinder(_)
                            | ogeom_geom::SurfaceGeometry::Torus(_)
                    )
                ) {
                    ogeom_core::ogeom_bail!(
                        Construction,
                        "a crossed junction's turned band is not a torus band"
                    );
                }
            }
            hulls.push(ogeom_algo::shape_bounds(model, wedge, tol)?.corners());
            continue;
        }
        for face in explore_unique(model, wedge, ShapeType::Face)? {
            let surface = model
                .node(&face)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface));
            let along = junction.outward[slot];
            let straight = match surface {
                Some(
                    ogeom_geom::SurfaceGeometry::Plane(_)
                    | ogeom_geom::SurfaceGeometry::Cylinder(_),
                ) => true,
                Some(ogeom_geom::SurfaceGeometry::Extrusion(e)) => {
                    e.direction().vector().cross(along).magnitude() <= tol.angular()
                }
                _ => false,
            };
            if !straight {
                ogeom_core::ogeom_bail!(
                    Construction,
                    "a crossed junction's band is not a straight band between planes"
                );
            }
        }
        let mut back = 0.0_f64;
        for edge in explore_unique(model, wedge, ShapeType::Edge)? {
            let (curve, range) = crate::support::edge_curve(model, &edge, tol)?;
            for point in curve.to_bspline_over(range, tol)?.control_points() {
                back = back.max(-(point.point() - at).dot(junction.outward[slot]));
            }
        }
        let shift = -junction.outward[slot] * back;
        let hull = sections[slot].hull();
        let far: Vec<Point> = hull.iter().map(|&q| q + shift).collect();
        hulls.push(hull.into_iter().chain(far).collect());
    }
    // How far the cell can reach: a point of one continuation a distance s
    // along stands s sin(angle) off the other's crease line, less its own
    // section's reach, so the two hold no common point past
    // (w0 + w1) / sin(angle) along either.
    let reach = [sections[0].reach(at), sections[1].reach(at)];
    let sine = junction.outward[0].cross(junction.outward[1]).magnitude();
    if sine <= tol.angular() {
        ogeom_core::ogeom_bail!(Construction, "two bands meeting across a turn run parallel");
    }
    let cell = (reach[0] + reach[1]) / sine + reach[0] + reach[1];
    // Where the walls' meeting line leans, a band carried on stands past
    // the other's wall's plane for a stretch past its cap: the material
    // there lies above both bands, beside the band and before the other's
    // cap. A section point q a distance s = (q - v).n past the other's
    // wall, the band leaving along u with u.n < 0 (the turn re-entrant),
    // stays past it for s / -(u.n) along. Its hull reaches the far end of
    // that stretch.
    let mut crossing: Vec<Vec<Point>> = Vec::with_capacity(2);
    for (section, (along, wall)) in sections.iter().zip(
        junction
            .outward
            .into_iter()
            .zip([junction.sides[1], junction.sides[0]]),
    ) {
        let rate = along.dot(wall);
        if rate >= -tol.angular() {
            ogeom_core::ogeom_bail!(
                Construction,
                "a band carried on across a re-entrant turn does not run behind the other's wall"
            );
        }
        crossing.push(
            section
                .hull()
                .into_iter()
                .filter_map(|q| {
                    let s = (q - at).dot(wall);
                    (s > 0.0).then(|| q + along * (s / -rate))
                })
                .collect(),
        );
    }
    let mut carried: [Option<Shape>; 2] = [None, None];
    let mut ahead = [0.0; 2];
    let mut spans = [0.0; 2];
    // A turned band first: past its cap its widening holds material past
    // the other's wall's plane (between its wall's tangent plane and the
    // cylinder, a bore's material), which the other band carried on must
    // reach. What stands off the cylinder past the other's wall is the
    // open corner, which the corner leaves.
    let order = junction.turned.map_or([0, 1], |t| [t.slot, 1 - t.slot]);
    for slot in order {
        let along = junction.outward[slot];
        let mut others: Vec<Point> = hulls[1 - slot]
            .iter()
            .chain(&crossing[1 - slot])
            .copied()
            .collect();
        if let Some(turned) = junction.turned.filter(|t| t.slot != slot)
            && let Some(widened) = &carried[turned.slot]
        {
            let bounds = ogeom_algo::shape_bounds(model, widened, tol)?;
            let size = bounds
                .corners()
                .iter()
                .map(|p| p.distance(at))
                .fold(0.0_f64, f64::max);
            let past = half_space(model, (at, junction.sides[slot]), size, tol)?;
            let beyond = half_space(model, (at, junction.outward[turned.slot]), size, tol)?;
            let past = ogeom_bool::common(model, &past, &beyond, tol)?.shape;
            let past = ogeom_bool::common(model, widened, &past, tol)?.shape;
            let mut held = ogeom_bool::common(model, &past, solid, tol)?.shape;
            if !explore_unique(model, &held, ShapeType::Face)?.is_empty() {
                let off = off_turned_wall(model, junction, turned, size, tol)?;
                held = ogeom_bool::cut(model, &held, &off, tol)?.shape;
            }
            if !explore_unique(model, &held, ShapeType::Face)?.is_empty() {
                others.extend(ogeom_algo::shape_bounds(model, &held, tol)?.corners());
            }
        }
        let own = crossing[slot]
            .iter()
            .map(|p| (*p - at).dot(along))
            .fold(0.0_f64, f64::max);
        // Back into its own wedge as far as the other's reaches, on past its
        // cap as far as the other's wedge or the cell does, and across past
        // its wall as far as the other's wedge stands from the vertex.
        let behind = others
            .iter()
            .map(|p| -(*p - at).dot(along))
            .fold(0.0_f64, f64::max);

        ahead[slot] = others
            .iter()
            .map(|p| (*p - at).dot(along))
            .fold(cell.max(own), f64::max)
            + reach[slot];
        let span = others
            .iter()
            .map(|p| p.distance(at))
            .fold(0.0_f64, f64::max)
            + reach[slot];
        spans[slot] = span;
        // A turned band is carried on from its cap, and widened before it
        // over its part below.
        let turned = junction.turned.is_some_and(|t| t.slot == slot);
        let (wedge, cap) = if behind > near && !turned {
            let start = at - along * behind;
            let trimmed = crate::support::trimmed_wedge(
                model,
                wedges[slot],
                (start, -along),
                behind + reach[slot],
                tol,
            )?;
            let Some(cap) = cap_at(model, &trimmed, start, along, tol)? else {
                ogeom_core::ogeom_bail!(
                    Construction,
                    "a band cut back for its junction's corner has no end there"
                );
            };
            (trimmed, cap)
        } else {
            (wedges[slot].clone(), caps[slot].clone())
        };
        carried[slot] = Some(carried_widened(
            model,
            &wedge,
            &cap,
            along,
            (at + along * ahead[slot], along),
            (junction.top, span),
            tol,
        )?);
    }
    let [Some(first_carried), Some(second_carried)] = carried else {
        ogeom_core::ogeom_bail!(
            Construction,
            "a band of a crossed junction was not carried on"
        );
    };
    let mut carried = [first_carried, second_carried];
    // Before its cap, what lies above a turned band is the band or what
    // stands off its wall, the cylinder: the band widened past its wall
    // turns with it. Over its part that widening is what stands off the
    // cylinder within the part's sector, under the shared host and before
    // the cap; past the part's sector the straight band carried on must
    // stand nowhere off the cylinder before the cap ([`beside_turned`]),
    // so nothing there lies above both.
    let size = ahead[0] + ahead[1] + spans[0] + spans[1];
    if let (Some(turned), Some(sector)) = (junction.turned, &sector) {
        let other = 1 - turned.slot;
        let along = junction.outward[other];
        let back = hulls[other]
            .iter()
            .chain(&crossing[other])
            .map(|p| (*p - at).dot(along))
            .fold(0.0_f64, f64::min);
        if !beside_turned(
            junction,
            turned,
            &sections[other].hull(),
            (back, ahead[other]),
            sector,
            near,
        ) {
            ogeom_core::ogeom_bail!(
                Construction,
                "a straight band carried on stands off a turned wall past its part's sector"
            );
        }
        let off = off_turned_wall(model, junction, turned, size, tol)?;
        let before = half_space(model, (at, -junction.outward[turned.slot]), size, tol)?;
        let under = half_space(model, (at, -junction.top), size, tol)?;
        let held = ogeom_bool::common(model, &sector.solid, &before, tol)?.shape;
        let held = ogeom_bool::common(model, &held, &under, tol)?.shape;
        let widening = ogeom_bool::common(model, &held, &off, tol)?.shape;
        if !explore_unique(model, &widening, ShapeType::Face)?.is_empty() {
            let slot = turned.slot;
            carried[slot] = ogeom_bool::fuse(model, &carried[slot], &widening, tol)?.shape;
        }
    }
    // Taken in an order of the bands' own directions, not the order the
    // request named the edges in, so either naming builds the same corner.
    let [u, v] = junction.outward;
    let first = usize::from(
        (u.x, u.y, u.z)
            .partial_cmp(&(v.x, v.y, v.z))
            .is_some_and(core::cmp::Ordering::is_gt),
    );
    let common = ogeom_bool::common(model, &carried[first], &carried[1 - first], tol)?;
    // What stands past both walls' planes is the open corner: the material
    // lies behind one wall or the other there. The widenings, each reaching
    // past its own wall, meet in it, and are taken out of it, so neither
    // need stop at the other's wall. Beside a turned band it is what stands
    // off the cylinder past the other's wall's plane: the cylinder bends
    // off its tangent plane, in towards the material of a bore.
    let past_first = half_space(model, (at, junction.sides[0]), size, tol)?;
    let past_second = half_space(model, (at, junction.sides[1]), size, tol)?;
    let open = match junction.turned {
        Some(turned) => {
            let off = off_turned_wall(model, junction, turned, size, tol)?;
            let other = if turned.slot == 0 {
                &past_second
            } else {
                &past_first
            };
            ogeom_bool::common(model, &off, other, tol)?.shape
        }
        None => ogeom_bool::common(model, &past_first, &past_second, tol)?.shape,
    };
    let cleared = ogeom_bool::cut(model, &common.shape, &open, tol)?;
    let history = common.history.then(&cleared.history);
    let taken = cleared.shape;
    if explore_unique(model, &taken, ShapeType::Face)?.is_empty() {
        ogeom_core::ogeom_bail!(Construction, "two bands carried on take no corner");
    }
    // The corner must lie in the material: nothing of it outside the solid
    // but what the two wedges themselves take there (a band run out past an
    // oblique end face at its far end).
    let mut outside = ogeom_bool::cut(model, &taken, solid, tol)?;
    for wedge in wedges {
        if explore_unique(model, &outside.shape, ShapeType::Face)?.is_empty() {
            break;
        }
        outside = ogeom_bool::cut(model, &outside.shape, wedge, tol)?;
    }
    if !explore_unique(model, &outside.shape, ShapeType::Face)?.is_empty() {
        ogeom_core::ogeom_bail!(
            Construction,
            "the corner two bands carried on take leaves the material"
        );
    }
    // Nor may it reach either continuation's far end: no face of it lies
    // in either far end's plane.
    for face in explore_unique(model, &taken, ShapeType::Face)? {
        for (along, far) in junction.outward.into_iter().zip(ahead) {
            if in_plane(model, &face, at + along * far, along, near)? {
                ogeom_core::ogeom_bail!(
                    Construction,
                    "the corner reaches a band continuation's far end"
                );
            }
        }
    }
    // Nor may it leave behind what lies above both bands for certain: each
    // band carried on, its wedge with it, where it stands past the other's
    // wall's plane. That is above the band, and past the other's wall above
    // the other, so the corner must hold it all but the open corner.
    for slot in 0..2 {
        let along = junction.outward[slot];
        let band = carried_to(
            model,
            wedges[slot],
            &caps[slot],
            along,
            (at + along * ahead[slot], along),
            tol,
        )?;
        let past = half_space(
            model,
            (at, junction.sides[1 - slot]),
            ahead[0] + ahead[1],
            tol,
        )?;
        let above = ogeom_bool::common(model, &band, &past, tol)?.shape;
        let mut left = ogeom_bool::cut(model, &above, &taken, tol)?.shape;
        if junction.turned.is_some() && !explore_unique(model, &left, ShapeType::Face)?.is_empty() {
            left = ogeom_bool::cut(model, &left, &open, tol)?.shape;
        }
        if !explore_unique(model, &left, ShapeType::Face)?.is_empty() {
            ogeom_core::ogeom_bail!(
                Construction,
                "the corner two bands carried on take leaves material above both"
            );
        }
    }
    // What both bands carried on hold needs no check of its own: each band
    // carried on lies within its widening, which reaches as far back and
    // as far on, and behind its own wall, out of the open corner past both,
    // so the corner holds it.
    let mut of_second = Vec::new();
    for face in explore_unique(model, &carried[1], ShapeType::Face)? {
        of_second.extend(history.trace(&face).iter().cloned());
    }
    let Some((slot, far, (plane, on))) = rest else {
        return Ok((taken, of_second));
    };
    // The turned band's rest joined back on across the half-plane it was
    // cut in, where the corner holds the band's own section.
    if slot == 1 {
        of_second.extend(explore_unique(model, &far, ShapeType::Face)?);
    }
    let (joined, traced) = joined(model, [&taken, &far], plane, on, near, tol)?;
    let of_second = traced
        .into_iter()
        .filter(|(face, _)| of_second.iter().any(|f| f.is_same(face)))
        .map(|(_, image)| image)
        .collect();
    Ok((joined, of_second))
}

/// What of a cube of side `2 size` round the vertex of `junction` stands
/// off its turned wall, the cylinder about the turned band's axis: outside
/// it round a boss, within it in a bore.
fn off_turned_wall(
    model: &mut Model,
    junction: &Junction,
    turned: Turned,
    size: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let axis = junction.top;
    let seed = if axis.x.abs() < 0.9 {
        Vector::new(1.0, 0.0, 0.0)
    } else {
        Vector::new(0.0, 1.0, 0.0)
    };
    let across = seed - axis * seed.dot(axis);
    let reach = size + (junction.at - turned.centre).magnitude();
    let frame = ogeom_math::Frame::new(
        turned.centre - axis * reach,
        ogeom_math::Direction::new(axis, tol)?,
        ogeom_math::Direction::new(across, tol)?,
        tol,
    )?;
    let drum = ogeom_algo::make_cylinder(model, frame, turned.big, reach * 2.0, tol)?.shape;
    let normal = junction.sides[turned.slot];
    let cube = half_space(model, (junction.at - normal * size, normal), size, tol)?;
    Ok(if turned.within {
        ogeom_bool::cut(model, &cube, &drum, tol)?.shape
    } else {
        ogeom_bool::common(model, &cube, &drum, tol)?.shape
    })
}

/// The quarter turn about a turned band's axis its part stands in: a box
/// on the half-plane the part was cut in, as a solid, the foot of the axis
/// there, and the unit normals, square to the axis, of the box's two faces
/// through it, each into the box.
pub(crate) struct Sector {
    solid: Shape,
    foot: Point,
    faces: [Vector; 2],
}

/// A turned band's wedge cut in two: its part at the vertex, the rest, and
/// the plane of the cut, as a point and a unit normal, with a point of the
/// cut on it, and the quarter turn the part stands in.
type Split = (Shape, Shape, ((Point, Vector), Point), Sector);

/// The turned band's wedge in `slot` of `junction`, its crease a circle
/// about the axis through `centre` along the shared host's normal, cut in
/// two in a half-plane from the axis: the part from its cap at the vertex
/// back to that half-plane, the rest, the plane with a point of the cut
/// on it, and the quarter turn the part stands in.
///
/// The band's section turned back about the axis lies within the hull of
/// its hull points turned alike. The cut stands at twice the first turn,
/// in steps of a sixty-fourth of a quarter turn up to half of one, at
/// which those points stand wholly past the other band's wall's plane, and
/// still do: any such turn cuts the corner off the same way, and twice it
/// keeps the cut clear of the plane. Past it the band stands above the
/// other band carried on, and the corner leaves it as built.
fn turned_split(
    model: &mut Model,
    wedge: &Shape,
    junction: &Junction,
    slot: usize,
    centre: Point,
    tol: Tolerances,
) -> OgeomResult<Split> {
    let at = junction.at;
    let near = tol.confusion() * 1e3;
    let axis = junction.top;
    let Some(cap) = cap_at(model, wedge, at, junction.outward[slot], tol)? else {
        ogeom_core::ogeom_bail!(
            Construction,
            "a band meeting another across a re-entrant turn has no cap at the vertex"
        );
    };
    let hull = Section::of(model, &cap, near, tol)?.hull();
    // The way round the axis the crease leaves the vertex.
    let sense = axis.cross(at - centre).dot(-junction.outward[slot]);
    if sense.abs() <= near {
        ogeom_core::ogeom_bail!(
            Construction,
            "a turned band's crease does not turn about its axis"
        );
    }
    let turn = |p: Point, angle: f64| -> Point {
        let v = p - centre;
        let (s, c) = (angle * sense.signum()).sin_cos();
        centre + v * c + axis.cross(v) * s + axis * (axis.dot(v) * (1.0 - c))
    };
    let (on, wall) = (at, junction.sides[1 - slot]);
    let quarter = core::f64::consts::FRAC_PI_2;
    let past = |angle: f64| hull.iter().all(|&q| (turn(q, angle) - on).dot(wall) > near);
    let Some(first) = (1..=32)
        .map(|k| quarter * f64::from(k) / 64.0)
        .find(|&angle| past(angle))
    else {
        ogeom_core::ogeom_bail!(
            Construction,
            "a turned band does not stand past the other's wall within a quarter turn"
        );
    };
    let angle = first * 2.0;
    if !past(angle) {
        ogeom_core::ogeom_bail!(
            Construction,
            "a turned band comes back behind the other's wall within a quarter turn"
        );
    }
    // A box on the half-plane from the axis, the vertex's side of it.
    let foot = centre + axis * (at - centre).dot(axis);
    let outward = {
        let v = turn(at, angle) - foot;
        v / v.magnitude()
    };
    let mut up = axis;
    if (at - foot).dot(up.cross(outward)) < 0.0 {
        up = -up;
    }
    let size = ogeom_algo::shape_bounds(model, wedge, tol)?.diagonal() + (at - foot).magnitude();
    let frame = ogeom_math::Frame::new(
        foot - up * size,
        ogeom_math::Direction::new(up, tol)?,
        ogeom_math::Direction::new(outward, tol)?,
        tol,
    )?;
    let sector = ogeom_algo::make_box(model, frame, (size, size, size * 2.0), tol)?.shape;
    let one = |model: &Model, shape: &Shape| -> OgeomResult<Shape> {
        let solids = explore_unique(model, shape, ShapeType::Solid)?;
        let [solid] = solids.as_slice() else {
            ogeom_core::ogeom_bail!(Construction, "a turned band cut in two fell apart");
        };
        Ok(solid.clone())
    };
    let part = ogeom_bool::common(model, wedge, &sector, tol)?.shape;
    let part = one(model, &part)?;
    let rest = ogeom_bool::cut(model, wedge, &sector, tol)?.shape;
    let rest = one(model, &rest)?;
    let normal = up.cross(outward);
    let quadrant = Sector {
        solid: sector,
        foot,
        faces: [outward, normal],
    };
    Ok((part, rest, ((foot, normal), turn(at, angle)), quadrant))
}

/// Whether the straight band of `junction`, beside its turned one, carried
/// on along its crease from `lo` to `hi` past the vertex, its section's
/// `hull` points across it, stands nowhere before the turned band's cap
/// and outside the part's `sector` off the turned wall, within `held`.
///
/// Seen along the axis the band lies within a strip, and before the cap
/// lies a half-plane, less `held`. Outside the sector lies the union of two
/// half-planes, each less `held`. Each meeting of the strip with the
/// half-plane before the cap and one of those is a convex polygon. A
/// boss's cylinder, the material within it, holds a convex polygon wholly
/// where it holds its corners; a bore's, the material outside it, where the
/// polygon's nearest point to the axis lies outside it.
fn beside_turned(
    junction: &Junction,
    turned: Turned,
    hull: &[Point],
    (lo, hi): (f64, f64),
    sector: &Sector,
    held: f64,
) -> bool {
    let axis = junction.top;
    let flat = |v: Vector| v - axis * v.dot(axis);
    let unit = |v: Vector| {
        let m = v.magnitude();
        (m > held).then(|| v / m)
    };
    let slot = turned.slot;
    let (Some(d), Some(t)) = (
        unit(flat(junction.outward[1 - slot])),
        unit(flat(junction.outward[slot])),
    ) else {
        return false;
    };
    let m = axis.cross(d);
    let at = junction.at;
    let (w_lo, w_hi) = hull.iter().fold((0.0_f64, 0.0_f64), |(a, b), q| {
        let w = (*q - at).dot(m);
        (a.min(w), b.max(w))
    });
    // Plan positions from the vertex: the strip's rectangle, clipped to the
    // half-plane before the cap.
    let strip = [
        d * lo + m * w_lo,
        d * hi + m * w_lo,
        d * hi + m * w_hi,
        d * lo + m * w_hi,
    ];
    let before = clipped(&strip, |v| v.dot(t) + held);
    let from_axis = flat(at - turned.centre);
    let from_foot = flat(at - sector.foot);
    sector.faces.iter().all(|&inward| {
        let inward = flat(inward);
        let polygon = clipped(&before, |v| (from_foot + v).dot(inward) + held);
        let placed: Vec<Vector> = polygon.iter().map(|&v| from_axis + v).collect();
        if turned.within {
            placed.iter().all(|v| v.magnitude() <= turned.big + held)
        } else {
            nearest_to_origin(&placed, axis) >= turned.big - held
        }
    })
}

/// What of the convex `polygon` stands where `g` is at most zero, `g`
/// affine: the polygon clipped to that half-plane.
fn clipped(polygon: &[Vector], g: impl Fn(Vector) -> f64) -> Vec<Vector> {
    let mut out = Vec::with_capacity(polygon.len() + 1);
    for k in 0..polygon.len() {
        let (a, b) = (polygon[k], polygon[(k + 1) % polygon.len()]);
        let (ga, gb) = (g(a), g(b));
        if ga <= 0.0 {
            out.push(a);
        }
        if (ga < 0.0) != (gb < 0.0) && ga != gb {
            out.push(a + (b - a) * (ga / (ga - gb)));
        }
    }
    out
}

/// How far the convex `polygon`, its corners in order in a plane square to
/// `axis` through the origin, stands from the origin: nothing if it holds
/// the origin, infinitely far if it is empty.
fn nearest_to_origin(polygon: &[Vector], axis: Vector) -> f64 {
    let n = polygon.len();
    if n == 0 {
        return f64::INFINITY;
    }
    let mut nearest = f64::INFINITY;
    let (mut left, mut right) = (false, false);
    for k in 0..n {
        let (a, b) = (polygon[k], polygon[(k + 1) % n]);
        let edge = b - a;
        let length = edge.dot(edge);
        let f = if length > 0.0 {
            (-a.dot(edge) / length).clamp(0.0, 1.0)
        } else {
            0.0
        };
        nearest = nearest.min((a + edge * f).magnitude());
        // The origin's side of each edge: on one side of all, it is inside.
        let side = edge.cross(-a).dot(axis);
        left |= side > 0.0;
        right |= side < 0.0;
    }
    if n >= 3 && !(left && right) {
        return 0.0;
    }
    nearest
}

/// The half-space past the plane through `on` with unit `normal`, as a
/// cube of side `2 size` standing on the plane with `on` at its middle.
fn half_space(
    model: &mut Model,
    (on, normal): (Point, Vector),
    size: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let z = ogeom_math::Direction::new(normal, tol)?;
    let seed = if normal.x.abs() < 0.9 {
        Vector::new(1.0, 0.0, 0.0)
    } else {
        Vector::new(0.0, 1.0, 0.0)
    };
    let x = ogeom_math::Direction::new(seed - normal * seed.dot(normal), tol)?;
    let y = normal.cross(x.vector());
    let origin = on - (x.vector() + y) * size;
    Ok(ogeom_algo::make_box(
        model,
        ogeom_math::Frame::new(origin, z, x, tol)?,
        (size * 2.0, size * 2.0, size * 2.0),
        tol,
    )?
    .shape)
}

/// A curved side of a cap: its curve, its range, its ends and its middle.
type ArcPiece = (ogeom_geom::Curve, (f64, f64), (Point, Point), Point);

/// A band's end section on its planar cap: its two straight legs, the
/// ball's arc between their far ends, and the cap's plane.
struct Section {
    lines: Vec<(Point, Point)>,
    arc: ogeom_geom::BSplineCurve,
    plane: (Point, Vector),
}

impl Section {
    /// The section on `cap`, its legs told from its arc within `near`.
    fn of(model: &Model, cap: &Shape, near: f64, tol: Tolerances) -> OgeomResult<Self> {
        let mut lines: Vec<(Point, Point)> = Vec::new();
        let mut pieces = Vec::new();
        for edge in explore_unique(model, cap, ShapeType::Edge)? {
            let (curve, range) = crate::support::edge_curve(model, &edge, tol)?;
            let ends = (curve.point_at(range.0, tol)?, curve.point_at(range.1, tol)?);
            let middle = curve.point_at(f64::midpoint(range.0, range.1), tol)?;
            if middle.distance(ends.0.midpoint(ends.1)) <= near {
                lines.push(ends);
            } else {
                pieces.push((curve, range, ends, middle));
            }
        }
        let arc = match pieces.as_slice() {
            [] => ogeom_core::ogeom_bail!(Construction, "a band's end cap has no arc"),
            [(curve, range, _, _)] => curve.to_bspline_over(*range, tol)?,
            // The ball's arc split in pieces (at its circle's own seam, which
            // follows the part's placement): one circle through them all,
            // from the end no other piece meets to the other.
            _ => Self::one_arc(&pieces, near, tol)?,
        };
        let Some(plane) = plane_of(model, cap, tol)? else {
            ogeom_core::ogeom_bail!(Construction, "a band's end cap is not planar");
        };
        Ok(Self { lines, arc, plane })
    }

    /// The section's corners and seventeen points along its arc, where it
    /// is placed against a plane.
    fn points(&self, tol: Tolerances) -> OgeomResult<Vec<Point>> {
        use ogeom_geom::Curve3d as _;
        let arc: ogeom_geom::Curve = self.arc.clone().into();
        let (lo, hi) = arc.domain();
        let mut points: Vec<Point> = self.lines.iter().flat_map(|&(p, q)| [p, q]).collect();
        for k in 0..=16 {
            points.push(arc.point_at((hi - lo).mul_add(f64::from(k) / 16.0, lo), tol)?);
        }
        Ok(points)
    }

    /// The arc the circular `pieces` (each a curve, its range, its ends and
    /// its middle) make end to end, as one circle's arc.
    fn one_arc(
        pieces: &[ArcPiece],
        near: f64,
        tol: Tolerances,
    ) -> OgeomResult<ogeom_geom::BSplineCurve> {
        let ends: Vec<Point> = pieces
            .iter()
            .flat_map(|(_, _, (a, b), _)| [*a, *b])
            .collect();
        let lone: Vec<Point> = ends
            .iter()
            .copied()
            .filter(|p| ends.iter().filter(|q| q.distance(*p) <= near).count() == 1)
            .collect();
        let [start, end] = lone.as_slice() else {
            ogeom_core::ogeom_bail!(Construction, "a band's end cap's arc pieces do not chain");
        };
        let circle = ogeom_math::Circle::through(*start, pieces[0].3, *end, tol)?;
        for (curve, range, _, middle) in pieces {
            if !matches!(curve, ogeom_geom::Curve::Circle(_)) {
                ogeom_core::ogeom_bail!(Construction, "a band's end cap has two curved sides");
            }
            for p in [*middle, curve.point_at(range.0, tol)?] {
                if (p.distance(circle.centre()) - circle.radius()).abs() > near {
                    ogeom_core::ogeom_bail!(
                        Construction,
                        "a band's end cap's arc pieces lie on no one circle"
                    );
                }
            }
        }
        let frame = circle.frame();
        let angle = |p: Point| {
            let d = p - frame.origin();
            d.dot(frame.y().vector()).atan2(d.dot(frame.x().vector()))
        };
        let turn = core::f64::consts::TAU;
        let a0 = angle(*start);
        let lift = |a: f64| a0 + (a - a0).rem_euclid(turn);
        let (am, a1) = (lift(angle(pieces[0].3)), lift(angle(*end)));
        // The arc through the pieces, one way round or the other.
        let range = if am <= a1 { (a0, a1) } else { (a1 - turn, a0) };
        ogeom_geom::Curve::Circle(ogeom_geom::CircleCurve::new(circle)).to_bspline_over(range, tol)
    }

    /// Points whose hull holds the section: its legs' ends and its arc's
    /// control points, the arc's weights all positive.
    fn hull(&self) -> Vec<Point> {
        self.lines
            .iter()
            .flat_map(|&(p, q)| [p, q])
            .chain(self.arc.control_points().iter().map(|w| w.point()))
            .collect()
    }

    /// The farthest the section reaches from `at`.
    fn reach(&self, at: Point) -> f64 {
        self.hull()
            .iter()
            .map(|p| p.distance(at))
            .fold(0.0, f64::max)
    }

    /// The section widened to bound, besides itself, what stands past the
    /// wall's plane under the shared host with outward unit normal `top`
    /// within `span`: the wall's leg carried on down past its contact, out
    /// past the wall, back up to the host's plane past the crease, and in
    /// along the host to the crease, whose own leg then runs on to its
    /// contact. With it, the sides that close the part it gains in the
    /// section's plane.
    fn widened(
        &self,
        top: Vector,
        span: f64,
        near: f64,
    ) -> OgeomResult<(Self, Vec<(Point, Point)>)> {
        let [(a, b), (c, d)] = self.lines.as_slice() else {
            ogeom_core::ogeom_bail!(Construction, "a band's end cap has no two straight legs");
        };
        let (crease, [far_0, far_1]) = if a.distance(*c) <= near {
            (*a, [*b, *d])
        } else if a.distance(*d) <= near {
            (*a, [*b, *c])
        } else if b.distance(*c) <= near {
            (*b, [*a, *d])
        } else if b.distance(*d) <= near {
            (*b, [*a, *c])
        } else {
            ogeom_core::ogeom_bail!(Construction, "a band's end cap's legs do not meet");
        };
        let unit = |v: Vector| -> OgeomResult<Vector> {
            let m = v.magnitude();
            if m <= near {
                ogeom_core::ogeom_bail!(Construction, "a band's end cap has a leg of no length");
            }
            Ok(v / m)
        };
        let (u0, u1) = (unit(far_0 - crease)?, unit(far_1 - crease)?);
        let (on_top, on_wall) = if u0.dot(top).abs() <= u1.dot(top).abs() {
            (far_0, far_1)
        } else {
            (far_1, far_0)
        };
        let down = unit(on_wall - crease)?;
        let mut out = unit(self.plane.1.cross(down))?;
        if (on_top - crease).dot(out) > 0.0 {
            out = -out;
        }
        let past = unit(crease - on_top)?;
        let deep = on_wall + down * span;
        let wide = deep + out * span;
        let level = crease + past * span;
        let gained = vec![
            (on_wall, deep),
            (deep, wide),
            (wide, level),
            (level, crease),
        ];
        let mut lines = gained.clone();
        lines.push((crease, on_top));
        let mut closing = gained;
        closing.push((crease, on_wall));
        Ok((
            Self {
                lines,
                arc: self.arc.clone(),
                plane: self.plane,
            },
            closing,
        ))
    }
}

/// `wedge` carried on past its planar end face `cap` along `along`, the
/// crease's direction there, to the plane through `on` with unit `normal`:
/// each point of the cap slides along `along` until it meets the plane, so
/// the band's far end is the cap's arc mapped obliquely onto the plane, an
/// affine image, and the continuation is a straight band of the cap's own
/// section. Its faces replace the cap; a side whose ends already stand in
/// the plane (the wall's leg on a vertical edge the plane holds) has
/// nothing to sweep and is left out.
pub(crate) fn carried_to(
    model: &mut Model,
    wedge: &Shape,
    cap: &Shape,
    along: Vector,
    plane: (Point, Vector),
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let mut rings = Rings {
        pool: Vec::new(),
        near: tol.confusion() * 1e3,
        tol,
    };
    let section = Section::of(model, cap, rings.near, tol)?;
    let mut faces: Vec<Shape> = explore_unique(model, wedge, ShapeType::Face)?
        .into_iter()
        .filter(|f| !f.is_same(cap))
        .collect();
    faces.extend(swept_past(model, &section, along, plane, &mut rings, tol)?);
    closed_wedge(model, &faces, rings.near, tol)
}

/// `wedge` carried on past its planar end face `cap` along `along` to the
/// plane through `on` with unit `normal`, the band kept on the normal's
/// side: [`carried_to`] where the whole cap stands before the plane.
///
/// A wall drafted so its foot runs out under the corner (two walls leaning
/// alike at a re-entrant turn) puts part of the cap past the plane, which
/// no point of it slides forwards onto. Along a straight `edge` between
/// two planar `hosts` the band is the straight band of its cap's section
/// on either side of the cap, so it is carried on past every point's
/// meeting with the plane, to a plane square to the crease, and cut back
/// to the plane: what stands of it on the normal's side, however far past
/// the plane it was carried. Every point of it past the plane lies within
/// the reach of the section's points from `on`, where the cut takes it.
/// Elsewhere the band before its cap is no straight band, and this refuses
/// as [`carried_to`] does.
pub(crate) fn carried_into(
    model: &mut Model,
    wedge: &Shape,
    cap: &Shape,
    (edge, hosts): (&Shape, &[Shape]),
    along: Vector,
    (on, normal): (Point, Vector),
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let near = tol.confusion() * 1e3;
    let section = Section::of(model, cap, near, tol)?;
    let rate = along.dot(normal);
    let points = section.points(tol)?;
    let slide = |q: Point| -> f64 { (on - q).dot(normal) / rate };
    let behind = rate.abs() > tol.angular() && points.iter().any(|&q| slide(q) < -near);
    let (curve, _) = crate::support::edge_curve(model, edge, tol)?;
    let mut straight = matches!(curve, ogeom_geom::Curve::Line(_)) && hosts.len() == 2;
    for host in hosts {
        straight = straight && plane_of(model, host, tol)?.is_some();
    }
    if !behind || !straight {
        return carried_to(model, wedge, cap, along, (on, normal), tol);
    }
    let back = points.iter().map(|&q| slide(q)).fold(0.0_f64, f64::min);
    let reach = points.iter().map(|&q| slide(q)).fold(0.0_f64, f64::max);
    if reach <= near {
        ogeom_core::ogeom_bail!(
            Construction,
            "a band's end section stands wholly past the plane its junction meets in"
        );
    }
    // Any reach past the farthest meeting gives the same band once cut
    // back; twice it keeps the band's far end clear of the plane.
    let far = section.plane.0 + along * (reach * 2.0);
    let carried = carried_to(model, wedge, cap, along, (far, along), tol)?;
    // What of a straight band stands past the plane is what the cut is
    // to take, so any box holding it takes the same: twice the farthest
    // the section's points stand from `on`, wherever along the band, holds
    // the arc's bulge between them too.
    let span = points
        .iter()
        .flat_map(|&q| [back, 0.0, reach * 2.0].map(|s| (q + along * s).distance(on)))
        .fold(near, f64::max);
    crate::support::trimmed_wedge(model, &carried, (on, normal), span * 2.0, tol)
}

/// `wedge` carried on past its planar end face `cap` as [`carried_to`]
/// carries it, its section widened past the wall under the shared host
/// with outward unit normal `top` within `span`: the wedge and what stands
/// above its band carried on, in one solid.
fn carried_widened(
    model: &mut Model,
    wedge: &Shape,
    cap: &Shape,
    along: Vector,
    plane: (Point, Vector),
    (top, span): (Vector, f64),
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let mut rings = Rings {
        pool: Vec::new(),
        near: tol.confusion() * 1e3,
        tol,
    };
    let (section, closing) =
        Section::of(model, cap, rings.near, tol)?.widened(top, span, rings.near)?;
    let mut faces: Vec<Shape> = explore_unique(model, wedge, ShapeType::Face)?
        .into_iter()
        .filter(|f| !f.is_same(cap))
        .collect();
    let mut edges = Vec::new();
    for (p, q) in closing {
        if let Some(e) = rings.segment(model, p, q)? {
            edges.push(e);
        }
    }
    faces.push(face_on(
        model,
        plane_through(section.plane, tol)?,
        edges,
        tol,
    )?);
    faces.extend(swept_past(model, &section, along, plane, &mut rings, tol)?);
    closed_wedge(model, &faces, rings.near, tol)
}

/// A face on `surface` bounded by `edges`, in any order.
fn face_on(
    model: &mut Model,
    surface: ogeom_geom::SurfaceGeometry,
    edges: Vec<Shape>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let wire = ogeom_algo::make_wire_unordered(model, &edges, tol)?.shape;
    Ok(ogeom_algo::make_face(model, surface, &[wire], tol)?.shape)
}

/// The plane through `on` with unit `normal`, as a surface.
fn plane_through(
    (on, normal): (Point, Vector),
    tol: Tolerances,
) -> OgeomResult<ogeom_geom::SurfaceGeometry> {
    Ok(ogeom_geom::PlaneSurface::over(
        ogeom_math::Plane::through(on, ogeom_math::Direction::new(normal, tol)?),
        (-1e3, 1e3),
        (-1e3, 1e3),
    )?
    .into())
}

/// The faces of `section` swept along `along` to the plane through `on`
/// with unit `normal`: each leg's sweep, the arc's, and the far end in the
/// plane. Every point of the section must slide forwards onto the plane, or
/// stand in it: one behind it would carry the band back over its own end.
fn swept_past(
    model: &mut Model,
    section: &Section,
    along: Vector,
    (on, normal): (Point, Vector),
    rings: &mut Rings,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    use ogeom_geom::{Curve, ExtrusionSurface};
    let rate = along.dot(normal);
    if rate.abs() <= tol.angular() {
        ogeom_core::ogeom_bail!(Construction, "a band runs parallel to the plane it meets");
    }
    let slide = |q: Point| -> f64 { (on - q).dot(normal) / rate };
    let image = |q: Point| -> Point { q + along * slide(q) };
    let near = rings.near;
    let arc: Curve = section.arc.clone().into();
    let (lo, hi) = arc.domain();
    if section.points(tol)?.iter().any(|&q| slide(q) < -near) {
        ogeom_core::ogeom_bail!(
            Construction,
            "a band's end section stands past the plane its junction meets in"
        );
    }
    let moved: Vec<ogeom_math::Weighted<Point>> = section
        .arc
        .control_points()
        .iter()
        .map(|w| ogeom_math::Weighted::new(image(w.point()), w.weight, tol))
        .collect::<OgeomResult<_>>()?;
    let far: Curve = ogeom_geom::BSplineCurve::rational(section.arc.knots().clone(), moved)?.into();
    let (a0, a1) = (arc.point_at(lo, tol)?, arc.point_at(hi, tol)?);
    let mut faces = Vec::new();
    // Each leg swept to the plane: a quad, or a triangle where one end
    // already stands in it, or nothing where both do.
    for &(p, q) in &section.lines {
        let (p2, q2) = (image(p), image(q));
        let mut edges = Vec::new();
        for (x, y) in [(p, q), (q, q2), (q2, p2), (p2, p)] {
            if let Some(e) = rings.segment(model, x, y)? {
                edges.push(e);
            }
        }
        if edges.len() < 3 {
            continue;
        }
        let plane_normal = (q - p).cross(along);
        let m = plane_normal.magnitude();
        if m <= tol.angular() {
            continue;
        }
        faces.push(face_on(
            model,
            plane_through((p, plane_normal / m), tol)?,
            edges,
            tol,
        )?);
    }
    // The band swept to the plane, on the extrusion of the arc.
    {
        let reach = [a0, a1]
            .iter()
            .map(|q| slide(*q))
            .fold(0.0_f64, f64::max)
            .max(near);
        let mut edges = vec![
            rings.edge(model, arc.clone())?,
            rings.edge(model, far.clone())?,
        ];
        for q in [a0, a1] {
            if let Some(e) = rings.segment(model, q, image(q))? {
                edges.push(e);
            }
        }
        let surface = ExtrusionSurface::over(
            arc.clone(),
            ogeom_math::Direction::new(along, tol)?,
            (-reach, reach * 2.0),
        )?;
        faces.push(face_on(model, surface.into(), edges, tol)?);
    }
    // The far end in the plane.
    {
        let mut edges = vec![rings.edge(model, far)?];
        for (p, q) in &section.lines {
            if let Some(e) = rings.segment(model, image(*p), image(*q))? {
                edges.push(e);
            }
        }
        faces.push(face_on(
            model,
            plane_through((on, normal), tol)?,
            edges,
            tol,
        )?);
    }
    Ok(faces)
}

/// `faces` sewn within `near` into one closed shell, every edge trimmed on
/// each face's surface, and made a solid.
fn closed_wedge(
    model: &mut Model,
    faces: &[Shape],
    near: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let sewn = ogeom_algo::sew_within(model, faces, near, tol)?;
    let [shell] = sewn.shells.as_slice() else {
        ogeom_core::ogeom_bail!(
            Construction,
            "a band carried on past its cap did not close into one wedge"
        );
    };
    if !ogeom_algo::is_shell_closed(model, shell)? {
        ogeom_core::ogeom_bail!(
            Construction,
            "a band carried on past its cap did not close into one wedge"
        );
    }
    for face in explore_unique(model, shell, ShapeType::Face)? {
        trimmed_on_its_surface(model, &face, near, tol)?;
    }
    ogeom_algo::state_pcurve_gaps(model, shell, tol)?;
    // A vertex is held at least as loosely as every edge it bounds, which
    // the gaps just stated may have widened.
    for edge in explore_unique(model, shell, ShapeType::Edge)? {
        let Some(held) = model.tolerance_of(&edge)? else {
            continue;
        };
        if let Some((a, b)) = ogeom_algo::edge_vertices(model, &edge)? {
            for vertex in [a, b] {
                model.widen(&vertex, held)?;
            }
        }
    }
    Ok(ogeom_algo::make_solid(model, std::slice::from_ref(shell))?.shape)
}

/// Give every edge of `face` a trim on the face's surface: the closed form
/// where there is one, otherwise a fit of the curve's projection, which
/// lies on the surface by construction and within `near` of it as built.
fn trimmed_on_its_surface(
    model: &mut Model,
    face: &Shape,
    near: f64,
    tol: Tolerances,
) -> OgeomResult<()> {
    use ogeom_geom::Transformable as _;
    let Some(data) = model.node(face).and_then(|n| n.data().as_face().cloned()) else {
        ogeom_core::ogeom_bail!(Dangling, "face is not in this model");
    };
    let Some(surface) = model.geometry().surface(data.surface).cloned() else {
        ogeom_core::ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };
    let unplace = face.transform(model.datums())?.inverse()?;
    for edge in explore_unique(model, face, ShapeType::Edge)? {
        let Some(edge_data) = model.node(&edge).and_then(|n| n.data().as_edge().cloned()) else {
            continue;
        };
        if edge_data
            .pcurve_for(data.surface, edge.location())
            .is_some()
        {
            continue;
        }
        let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = edge_data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve).cloned() else {
            continue;
        };
        let range = *range;
        let placed = geometry.transformed(&(unplace * edge.transform(model.datums())?), tol)?;
        let pcurve = match ogeom_intersect::exact_pcurve_of(&placed, &surface, tol) {
            Some(exact) => exact,
            None => {
                ogeom_algo::pcurve_fit::fit_projected_pcurve_capped(
                    &placed, range, &surface, near, tol,
                )?
                .0
            }
        };
        ogeom_algo::attach_pcurve(
            model,
            &edge,
            pcurve,
            data.surface,
            edge.location().clone(),
            range,
        )?;
    }
    Ok(())
}

/// The rings of a continuation's faces, their vertices one per point so
/// the faces share them.
struct Rings {
    pool: Vec<(Point, Shape)>,
    near: f64,
    tol: Tolerances,
}

impl Rings {
    /// The vertex at `p`: one already made within reach, widened to reach
    /// `p` too, or a new one.
    fn vertex(&mut self, model: &mut Model, p: Point) -> OgeomResult<Shape> {
        if let Some((q, v)) = self.pool.iter().find(|(q, _)| q.distance(p) <= self.near) {
            let gap = q.distance(p);
            if gap > self.tol.confusion() {
                model.widen(v, ogeom_core::Tolerance::new(gap * 2.0)?)?;
            }
            return Ok(v.clone());
        }
        let v = ogeom_algo::make_vertex(model, p).shape;
        self.pool.push((p, v.clone()));
        Ok(v)
    }

    fn edge(&mut self, model: &mut Model, curve: ogeom_geom::Curve) -> OgeomResult<Shape> {
        let domain = curve.domain();
        let from = self.vertex(model, curve.point_at(domain.0, self.tol)?)?;
        let to = self.vertex(model, curve.point_at(domain.1, self.tol)?)?;
        Ok(ogeom_algo::make_edge_between(model, curve, domain, &from, &to, self.tol)?.shape)
    }

    fn segment(&mut self, model: &mut Model, p: Point, q: Point) -> OgeomResult<Option<Shape>> {
        if p.distance(q) <= self.near {
            return Ok(None);
        }
        let line = ogeom_geom::LineCurve::segment(p, q, self.tol)?;
        Ok(Some(self.edge(model, line.into())?))
    }
}

/// A planar face's plane, placed: a point on it and its unit normal.
fn plane_of(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Option<(Point, Vector)>> {
    use ogeom_geom::Transformable as _;
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
        return Ok(None);
    };
    let Some(surface @ ogeom_geom::SurfaceGeometry::Plane(_)) =
        model.geometry().surface(data.surface).cloned()
    else {
        return Ok(None);
    };
    let ogeom_geom::SurfaceGeometry::Plane(placed) =
        surface.transformed(&face.transform(model.datums())?, tol)?
    else {
        return Ok(None);
    };
    let plane = placed.plane();
    Ok(Some((plane.origin(), plane.normal().vector())))
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;

    /// A box's top edge's wedge carried on past its cap: every vertex of
    /// the band carried on is held at least as loosely as the edges it
    /// bounds, so the solid passes `check`.
    #[test]
    fn a_band_carried_on_past_its_cap_checks_valid() {
        let tol = Tolerances::millimetres();
        let mut model = Model::new();
        let block =
            ogeom_algo::make_box(&mut model, ogeom_math::Frame::WORLD, (4.0, 4.0, 2.0), tol)
                .unwrap()
                .shape;
        let edge = explore_unique(&model, &block, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let b = ogeom_algo::shape_bounds(&model, e, tol).unwrap();
                let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
                (lo.z - 2.0).abs() < 1e-6
                    && (hi.z - 2.0).abs() < 1e-6
                    && lo.y.abs() < 1e-6
                    && hi.y.abs() < 1e-6
            })
            .unwrap();
        let (built, wedges) = crate::support::collecting_wedges(|| {
            crate::fillet::fillet_edge(&mut model, &block, &edge, 0.5, tol)
        });
        built.unwrap();
        let [wedge] = wedges.as_slice() else {
            panic!("one wedge");
        };
        let at = Point::new(4.0, 0.0, 2.0);
        let along = Vector::new(1.0, 0.0, 0.0);
        let cap = cap_at(&model, &wedge.solid, at, along, tol)
            .unwrap()
            .unwrap();
        let carried = carried_to(
            &mut model,
            &wedge.solid,
            &cap,
            along,
            (at + along * 1.5, Vector::new(0.6, 0.8, 0.0)),
            tol,
        )
        .unwrap();
        let diagnosis = ogeom_algo::check(&model, &carried, tol).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    }
}
