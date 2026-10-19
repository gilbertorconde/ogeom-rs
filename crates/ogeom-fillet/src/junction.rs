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
//! marched band is held to. Elsewhere (a drafted wall beside a square one,
//! a curved shared host) the sections differ and no junction is taken: the
//! bands meet as they did without one.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Orientation, Shape, ShapeType, explore_unique};

use crate::fillet::Mate;

/// Two edges of a request meeting at a vertex where their creases turn.
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
    if reentrant {
        // Carried on, each band is a straight band of its end section, so
        // the reflection alone decides the meeting. The cell it closes
        // reaches from the vertex to where the bands' contacts on the
        // shared host meet, a radius under it: the material must hold that
        // reach, or the continuations run out of the solid (a turn folding
        // back on itself) and the cell is not the corner's.
        if miss > fit {
            return Ok(None);
        }
        let probe = crate::support::material_probe(model, solid, tol)?;
        // The vertex itself stands on the solid's edge; its neighbourhood
        // past it is what the cell takes.
        let under = top * (radius * 0.5);
        for k in 1..=8 {
            let p = at + (meet - at) * (f64::from(k) / 8.0) - under;
            if probe.holds(model, p, tol)? != ogeom_algo::Containment::In {
                return Ok(None);
            }
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
            return Ok(None);
        }
    }
    Ok(Some(Junction {
        at,
        pair,
        outward,
        reentrant,
        flush: (meet, bisector),
    }))
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
    for wedge in wedges {
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
    (on, normal): (Point, Vector),
    tol: Tolerances,
) -> OgeomResult<Shape> {
    use ogeom_geom::{Curve, ExtrusionSurface, PlaneSurface};
    let rate = along.dot(normal);
    if rate.abs() <= tol.angular() {
        ogeom_core::ogeom_bail!(Construction, "a band runs parallel to the plane it meets");
    }
    let slide = |q: Point| -> f64 { (on - q).dot(normal) / rate };
    let image = |q: Point| -> Point { q + along * slide(q) };
    let near = tol.confusion() * 1e3;
    // The cap's edges: its two legs are straight, and the third is the
    // ball's arc.
    let mut lines: Vec<(Point, Point)> = Vec::new();
    let mut arc: Option<(Curve, (f64, f64))> = None;
    for edge in explore_unique(model, cap, ShapeType::Edge)? {
        let (curve, range) = crate::support::edge_curve(model, &edge, tol)?;
        let ends = (curve.point_at(range.0, tol)?, curve.point_at(range.1, tol)?);
        let middle = curve.point_at(f64::midpoint(range.0, range.1), tol)?;
        if middle.distance(ends.0.midpoint(ends.1)) <= near {
            lines.push(ends);
        } else if arc.is_none() {
            arc = Some((curve, range));
        } else {
            ogeom_core::ogeom_bail!(Construction, "a band's end cap has two curved sides");
        }
    }
    let Some((arc, range)) = arc else {
        ogeom_core::ogeom_bail!(Construction, "a band's end cap has no arc");
    };
    // Every point of the cap slides forwards onto the plane, or stands in
    // it: one behind it would carry the band back over its own end.
    let mut corners: Vec<Point> = lines.iter().flat_map(|&(p, q)| [p, q]).collect();
    for k in 0..=16 {
        corners.push(arc.point_at(
            (range.1 - range.0).mul_add(f64::from(k) / 16.0, range.0),
            tol,
        )?);
    }
    if corners.iter().any(|&q| slide(q) < -near) {
        ogeom_core::ogeom_bail!(
            Construction,
            "a band's end section stands past the plane its junction meets in"
        );
    }
    let arc = arc.to_bspline_over(range, tol)?;
    let moved: Vec<ogeom_math::Weighted<Point>> = arc
        .control_points()
        .iter()
        .map(|w| ogeom_math::Weighted::new(image(w.point()), w.weight, tol))
        .collect::<OgeomResult<_>>()?;
    let far: Curve = ogeom_geom::BSplineCurve::rational(arc.knots().clone(), moved)?.into();
    let arc: Curve = arc.into();
    let (a0, a1) = (
        arc.point_at(arc.domain().0, tol)?,
        arc.point_at(arc.domain().1, tol)?,
    );
    let mut rings = Rings {
        pool: Vec::new(),
        near,
        tol,
    };
    let mut faces: Vec<Shape> = explore_unique(model, wedge, ShapeType::Face)?
        .into_iter()
        .filter(|f| !f.is_same(cap))
        .collect();
    let face_on = |model: &mut Model,
                   surface: ogeom_geom::SurfaceGeometry,
                   edges: Vec<Shape>|
     -> OgeomResult<Shape> {
        let wire = ogeom_algo::make_wire_unordered(model, &edges, tol)?.shape;
        Ok(ogeom_algo::make_face(model, surface, &[wire], tol)?.shape)
    };
    // Each leg swept to the plane: a quad, or a triangle where one end
    // already stands in it, or nothing where both do.
    for &(p, q) in &lines {
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
        let surface = PlaneSurface::over(
            ogeom_math::Plane::through(p, ogeom_math::Direction::new(plane_normal / m, tol)?),
            (-1e3, 1e3),
            (-1e3, 1e3),
        )?;
        faces.push(face_on(model, surface.into(), edges)?);
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
        faces.push(face_on(model, surface.into(), edges)?);
    }
    // The far end in the plane.
    {
        let mut edges = vec![rings.edge(model, far)?];
        for (p, q) in &lines {
            if let Some(e) = rings.segment(model, image(*p), image(*q))? {
                edges.push(e);
            }
        }
        let surface = PlaneSurface::over(
            ogeom_math::Plane::through(on, ogeom_math::Direction::new(normal, tol)?),
            (-1e3, 1e3),
            (-1e3, 1e3),
        )?;
        faces.push(face_on(model, surface.into(), edges)?);
    }
    let sewn = ogeom_algo::sew_within(model, &faces, near, tol)?;
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
