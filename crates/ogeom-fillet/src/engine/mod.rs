//! The blend engine: one model for fillets and chamfers, built in fixed
//! phases from the solid as given and the whole request
//! (docs/adr/0001-blend-architecture.md).
//!
//! - P0 [`classify`]: the solid read once, every element keyed by its node,
//!   every edge's seat, every vertex's class.
//! - P1 [`sections`]: each blend's [`SectionLaw`](sections::SectionLaw):
//!   its contact's setback on each face, and its strip.
//! - P2 [`setbacks`]: each face's setback corners, where its contact lines
//!   meet: blends end by intersection, not by case.
//! - P3 [`vertices`]: each vertex closed by one
//!   [`VertexMethod`](vertices::VertexMethod).
//! - P4 [`assemble`]: the result's faces built once over shared points and
//!   edges, each original face keeping its surface.
//! - P5 [`guards`]: checks independent of the construction.
//!
//! The backbone is the face-setback structure: every face shrinks to its
//! setback corners, every blended edge becomes a strip, every vertex a
//! vertex blend. Nothing depends on the order of the request.
//!
//! The polyhedral classes (a solid bounded by planes) go through here.

mod assemble;
mod classify;
mod consumed;
#[cfg_attr(
    not(test),
    allow(
        dead_code,
        reason = "rounds of curved solids go through it from the rim layer on"
    )
)]
mod curved;
mod guards;
mod request;
mod rims;
mod round_strips;
mod rounds;
mod seams;
mod sections;
mod setbacks;
mod trimmed;
mod vertices;

pub(crate) use request::{Request, Section};
pub(crate) use rims::assemble as round_rims;
pub(crate) use round_strips::round as round_curved;
pub(crate) use seams::{Dissolved, dissolve, dissolve_curved, is_seam};
pub(crate) use trimmed::{
    chamfer as chamfer_curved, chamfer_mitred as chamfer_curved_mitred, chamfer_patched,
};

/// Whether `solid` is one the polyhedral classes read: one solid bounded by
/// planes, each face its rings of straight edges.
///
/// # Errors
///
/// As reading the solid.
pub(crate) fn bounded_by_planes(
    model: &Model,
    solid: &Shape,
    tol: Tolerances,
) -> OgeomResult<bool> {
    Ok(classify::read_planar(model, solid, tol)?.is_some())
}

use ogeom_algo::{Built, History};
use ogeom_core::{FastMap, FastSet, OgeomResult, Tolerances};
use ogeom_topo::{Model, Shape};

use assemble::Rings;
use sections::{Strip, StripCorners};

/// A node of the solid, by its arena index.
pub(crate) type Id = u32;

pub(crate) fn id(shape: &Shape) -> Id {
    shape.node().index()
}

/// Where a ring of the result comes from, for its history.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Origin {
    /// A face of the solid, shrunk to its setback corners.
    Face(usize),
    /// A blended edge's strip.
    Edge(Id),
    /// A vertex's blend.
    Vertex(Id),
}

/// The blend of `request` on `solid`, or `None` where the solid or the
/// request lies outside the engine's classes (for now: a solid bounded by
/// planes, each face one ring of straight edges and a ring per hole (rounds:
/// no holes), every edge between two
/// faces, every vertex's faces closing round it once; and the request
/// inside the setback domain: a face its setbacks consume convex and every
/// bevel and vertex round it convex ([`consumed`]), no two faces crossing,
/// every
/// chamfer corner's setback corners on one plane).
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// a vertex needs a method not built yet, or the built solid fails a guard.
pub(crate) fn run(
    model: &mut Model,
    solid: &Shape,
    request: &Request,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    // P0.
    let Some(planar) = classify::read_planar(model, solid, tol)? else {
        return Ok(None);
    };
    let Some(seats) = planar.seats() else {
        return Ok(None);
    };
    let mut by_edge: FastMap<Id, &Section> = FastMap::default();
    for blend in &request.blends {
        by_edge.insert(id(&blend.edge), &blend.section);
    }
    let blended: FastSet<Id> = by_edge.keys().copied().collect();
    let classes = planar.classes(&seats, &blended);
    // P1: each blended edge's setbacks on its two faces.
    let mut setback: FastMap<(Id, Id), f64> = FastMap::default();
    let mut keys: Vec<(Id, Id)> = seats.keys().copied().collect();
    keys.sort_unstable();
    for key in &keys {
        let seat = &seats[key];
        let Some(section) = by_edge.get(&seat.edge) else {
            continue;
        };
        let d = sections::law(section).setbacks(seat, &planar)?;
        for (k, face) in seat.faces.iter().enumerate() {
            setback.insert((id(&planar.faces[*face].shape), seat.edge), d[k]);
        }
    }
    // P2.
    let Some(corners) = setbacks::setback_corners(&planar, &setback, tol)? else {
        return Ok(None);
    };
    // Rounds are curved: their strips and vertex blends go into a sheet.
    if request
        .blends
        .iter()
        .any(|b| matches!(b.section, Section::Round { .. }))
    {
        if planar.faces.iter().any(|f| !f.holes.is_empty()) {
            return Ok(None);
        }
        return rounds::assemble(model, &planar, &seats, &by_edge, &corners, tol);
    }
    // The mitred vertices, each of three chamfers or more, and the point
    // their planes share; where four or more share no one point, the
    // vertex is cut back by every plane round it instead (`scattered`),
    // each bevel running on until the others cut it.
    let mut mitres: FastMap<Id, ogeom_math::Point> = FastMap::default();
    let mut scattered: FastSet<Id> = FastSet::default();
    if request.mitre {
        let mut vertices: Vec<Id> = planar.points.keys().copied().collect();
        vertices.sort_unstable();
        for v in vertices {
            let at: Vec<&classify::Seat> = keys
                .iter()
                .filter(|(a, b)| *a == v || *b == v)
                .map(|key| &seats[key])
                .filter(|seat| blended.contains(&seat.edge))
                .collect();
            if at.len() < 3 {
                continue;
            }
            let planes: Vec<(ogeom_math::Vector, f64)> = at
                .iter()
                .map(|seat| {
                    let (a, b) = seat.ends;
                    let (p, q, r) = (
                        corners[&(seat.faces[0], a)],
                        corners[&(seat.faces[0], b)],
                        corners[&(seat.faces[1], a)],
                    );
                    let n = (q - p).cross(r - p);
                    let n = n / n.magnitude();
                    (n, n.dot(p.to_vector()))
                })
                .collect();
            let ((a, da), (b, db), (c, dc)) = (planes[0], planes[1], planes[2]);
            let det = a.dot(b.cross(c));
            if det.abs() <= 1e-12 {
                return Ok(None);
            }
            let x = (b.cross(c) * da + c.cross(a) * db + a.cross(b) * dc) / det;
            // Past three bevels, the point only where every plane passes
            // through it (a pyramid's apex bevelled alike).
            if planes[3..]
                .iter()
                .any(|(n, d)| (n.dot(x) - d).abs() > tol.confusion())
            {
                scattered.insert(v);
            } else {
                mitres.insert(v, ogeom_math::Point::from_vector(x));
            }
        }
    }
    let cluster = match consumed::consumed(
        &planar, &seats, &blended, &corners, &classes, &scattered, tol,
    ) {
        consumed::Consumed::Nothing => None,
        consumed::Consumed::Cluster(cluster) => Some(cluster),
        consumed::Consumed::Outside => return Ok(None),
    };
    // A mitre to one point inside a cut-back cluster is not built.
    if cluster.is_some() && !mitres.is_empty() {
        return Ok(None);
    }
    let eaten = |fi: usize| cluster.as_ref().is_some_and(|c| c.faces.contains(&fi));
    let mut rings = Rings::new(tol);
    let mut sorted: Vec<((usize, Id), ogeom_math::Point)> =
        corners.iter().map(|(k, p)| (*k, *p)).collect();
    sorted.sort_unstable_by_key(|(k, _)| *k);
    let at: FastMap<(usize, Id), usize> = sorted
        .into_iter()
        .map(|(key, p)| (key, rings.point(p)))
        .collect();
    // A consumed face starts from its own outline, which the clip below
    // cuts back to what survives of it.
    let mut outlines: FastMap<usize, Option<Id>> = FastMap::default();
    for (fi, face) in planar.faces.iter().enumerate() {
        if eaten(fi) {
            let ring: Vec<usize> = face
                .corners
                .iter()
                .map(|v| {
                    let i = rings.point(planar.points[v]);
                    outlines.insert(i, Some(*v));
                    i
                })
                .collect();
            rings.push(ring, Origin::Face(fi), Some(face.normal));
            continue;
        }
        let mut shrunk: Vec<Vec<usize>> = face
            .loops()
            .map(|corners| corners.iter().map(|v| at[&(fi, *v)]).collect())
            .collect();
        let outer = shrunk.remove(0);
        rings.push_holed(outer, shrunk, Origin::Face(fi), Some(face.normal));
    }
    // How far a bevel runs on past a corner the clip cuts back: across the
    // whole solid, so no plane round the corner trims it short.
    let reach = {
        let (mut lo, mut hi) = (
            ogeom_math::Vector::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
            ogeom_math::Vector::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
        );
        for p in planar.points.values() {
            let v = p.to_vector();
            lo = ogeom_math::Vector::new(lo.x.min(v.x), lo.y.min(v.y), lo.z.min(v.z));
            hi = ogeom_math::Vector::new(hi.x.max(v.x), hi.y.max(v.y), hi.z.max(v.z));
        }
        (hi - lo).magnitude()
    };
    let mut run_on: FastMap<usize, Option<Id>> = FastMap::default();
    // The strips, each by its law.
    for key in &keys {
        let seat = &seats[key];
        let Some(section) = by_edge.get(&seat.edge) else {
            continue;
        };
        let (a, b) = seat.ends;
        // At a corner the clip cuts back, the bevel runs on past it, far
        // enough for the planes there to trim it to what it keeps; the
        // run-on lies past the solid, so any of them may take it.
        let mut end = |face: usize, v: Id, sign: f64| -> usize {
            if cluster
                .as_ref()
                .is_some_and(|c| c.vertices.contains(&v) && !c.grouped.contains(&v))
            {
                let p = corners[&(face, v)] + seat.along * (sign * reach);
                let i = rings.point(p);
                run_on.insert(i, None);
                i
            } else {
                at[&(face, v)]
            }
        };
        let corners = StripCorners {
            first: [end(seat.faces[0], a, -1.0), end(seat.faces[0], b, 1.0)],
            second: [end(seat.faces[1], a, -1.0), end(seat.faces[1], b, 1.0)],
        };
        match sections::law(section).strip(seat, corners)? {
            Strip::Planar(mut ring) => {
                // At a mitred vertex the bevel runs on to the point its
                // plane shares with the other two there: after its first
                // face's corner at the start, after its second's at the
                // end.
                if let Some(m) = mitres.get(&b) {
                    let i = rings.point(*m);
                    ring.push(i);
                }
                if let Some(m) = mitres.get(&a) {
                    let i = rings.point(*m);
                    ring.insert(2, i);
                }
                let facing = seat.normals[0] + seat.normals[1];
                rings.push(ring, Origin::Edge(seat.edge), Some(facing));
            }
        }
    }
    // P3.
    let mut vertices: Vec<Id> = planar.points.keys().copied().collect();
    vertices.sort_unstable();
    for v in vertices {
        let Some(order) = classify::faces_round(&planar.faces, v) else {
            return Ok(None);
        };
        let met: Vec<&Section> = keys
            .iter()
            .filter(|(a, b)| *a == v || *b == v)
            .filter_map(|key| by_edge.get(&seats[key].edge).copied())
            .collect();
        if cluster.as_ref().is_some_and(|c| c.unfaceted.contains(&v)) || mitres.contains_key(&v) {
            continue;
        }
        match vertices::method(classes[&v], &met)? {
            vertices::VertexMethod::Facet => {
                let ring = order.iter().map(|fi| at[&(*fi, v)]).collect();
                rings.push(ring, Origin::Vertex(v), None);
            }
        }
    }
    // Round consumed faces every ring is cut back by the planes there.
    if let Some(cluster) = &cluster {
        let mut allowed: FastMap<usize, Option<Id>> = outlines;
        allowed.extend(run_on.iter().map(|(i, v)| (*i, *v)));
        for ((_, v), i) in &at {
            if cluster.vertices.contains(v) {
                allowed.insert(*i, Some(*v));
            }
        }
        if !rings.clip_local(
            &cluster.cuts,
            &cluster.groups,
            &mut allowed,
            |ring, at, cut| cluster.may_take(ring, at, cut),
        ) {
            return Ok(None);
        }
    }
    // P4.
    let Some(assembled) = rings.assemble()? else {
        return Ok(None);
    };
    let planes: Vec<Option<ogeom_math::Plane>> = assembled
        .origins
        .iter()
        .map(|from| match from.first() {
            Some(Origin::Face(fi)) => planar.faces[*fi].plane,
            _ => None,
        })
        .collect();
    let (built, faces) = if assembled.holes.iter().all(Vec::is_empty) {
        ogeom_algo::make_faceted_solid(model, &assembled.points, &assembled.rings, &planes, tol)?
    } else {
        let holed: Vec<Vec<Vec<usize>>> = assembled
            .rings
            .iter()
            .zip(&assembled.holes)
            .map(|(ring, inner)| {
                std::iter::once(ring.clone())
                    .chain(inner.iter().cloned())
                    .collect()
            })
            .collect();
        ogeom_algo::make_faceted_solid_with_holes(model, &assembled.points, &holed, &planes, tol)?
    };
    let mut history = History::new();
    let mut gone: Vec<Id> = blended.iter().copied().collect();
    gone.sort_unstable();
    for e in gone {
        history.delete(&planar.shapes[&e]);
    }
    for (face, from) in faces.iter().zip(&assembled.origins) {
        for origin in from {
            match origin {
                Origin::Face(fi) => history.modify(&planar.faces[*fi].shape, face.clone()),
                Origin::Edge(e) => history.generate(&planar.shapes[e], face.clone()),
                Origin::Vertex(v) => history.generate(&planar.shapes[v], face.clone()),
            }
        }
    }
    // An edge left as it was is the new edge between the faces its two faces
    // became; a vertex left where it was is the new vertex there.
    let into: FastMap<usize, Vec<Shape>> = {
        let mut into: FastMap<usize, Vec<Shape>> = FastMap::default();
        for (face, from) in faces.iter().zip(&assembled.origins) {
            for origin in from {
                if let Origin::Face(fi) = origin {
                    into.entry(*fi).or_default().push(face.clone());
                }
            }
        }
        into
    };
    for key in &keys {
        let seat = &seats[key];
        if blended.contains(&seat.edge) {
            continue;
        }
        let (Some(first), Some(second)) = (into.get(&seat.faces[0]), into.get(&seat.faces[1]))
        else {
            continue;
        };
        for e in ogeom_topo::explore_unique(model, &built.shape, ogeom_topo::ShapeType::Edge)? {
            let around =
                ogeom_topo::ancestors_of(model, &built.shape, &e, ogeom_topo::ShapeType::Face)?;
            let on = |held: &[Shape]| around.iter().any(|f| held.iter().any(|g| g.is_same(f)));
            if on(first) && on(second) {
                history.modify(&planar.shapes[&seat.edge], e);
            }
        }
    }
    let mut kept: Vec<Id> = planar.points.keys().copied().collect();
    kept.sort_unstable();
    let new_vertices =
        ogeom_topo::explore_unique(model, &built.shape, ogeom_topo::ShapeType::Vertex)?;
    for v in kept {
        let at = planar.points[&v];
        for w in &new_vertices {
            if crate::setback::point_of(model, w)?.distance(at) <= tol.confusion() {
                history.modify(&planar.shapes[&v], w.clone());
            }
        }
    }
    // What a consumed face leaves nothing of (the face, its sharp edges and
    // its vertices), and a mitred corner's vertex, where nothing of the
    // result came of them.
    if let Some(cluster) = &cluster {
        let mut gone: Vec<&Shape> = Vec::new();
        let mut eaten: Vec<usize> = cluster.faces.iter().copied().collect();
        eaten.sort_unstable();
        for fi in eaten {
            let face = &planar.faces[fi];
            gone.push(&face.shape);
            let n = face.corners.len();
            for k in 0..n {
                let (a, b) = (face.corners[k], face.corners[(k + 1) % n]);
                gone.push(&planar.shapes[&planar.edges[&(a.min(b), a.max(b))]]);
                gone.push(&planar.shapes[&a]);
            }
        }
        let mut mitred: Vec<Id> = cluster.vertices.iter().copied().collect();
        mitred.sort_unstable();
        gone.extend(mitred.iter().map(|v| &planar.shapes[v]));
        for shape in gone {
            if !history.is_affected(shape) {
                history.delete(shape);
            }
        }
    }
    // P5.
    guards::solid_fits(model, &built.shape, tol)?;
    Ok(Some(Built::new(built.shape, history)))
}
