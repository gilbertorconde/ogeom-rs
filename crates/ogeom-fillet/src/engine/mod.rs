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
mod guards;
mod request;
mod rounds;
mod seams;
mod sections;
mod setbacks;
mod vertices;

pub(crate) use request::{Request, Section};
pub(crate) use seams::{Dissolved, dissolve, is_seam};

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
/// planes, each face one ring of straight edges, every edge between two
/// faces, every vertex's faces closing round it once; and the request
/// inside the setback domain: no face consumed by its setbacks but a band
/// between two parallel bevels ([`consumed`]), no two faces crossing, every
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
    let Some(mut corners) = setbacks::setback_corners(&planar, &setback, tol)? else {
        return Ok(None);
    };
    // Rounds are curved: their strips and vertex blends go into a sheet.
    if request
        .blends
        .iter()
        .any(|b| matches!(b.section, Section::Round { .. }))
    {
        return rounds::assemble(model, &planar, &seats, &by_edge, &corners, tol);
    }
    let Some(consumed) = consumed::consumed(&planar, &seats, &setback, &mut corners, tol) else {
        return Ok(None);
    };
    let mut rings = Rings::new(tol);
    let mut sorted: Vec<((usize, Id), ogeom_math::Point)> =
        corners.iter().map(|(k, p)| (*k, *p)).collect();
    sorted.sort_unstable_by_key(|(k, _)| *k);
    let at: FastMap<(usize, Id), usize> = sorted
        .into_iter()
        .map(|(key, p)| (key, rings.point(p)))
        .collect();
    for (fi, face) in planar.faces.iter().enumerate() {
        if consumed.faces.contains(&fi) {
            continue;
        }
        let n = face.corners.len();
        let mut ring = Vec::with_capacity(n);
        for k in 0..n {
            let (v, w) = (face.corners[k], face.corners[(k + 1) % n]);
            ring.push(at[&(fi, v)]);
            if let Some(ridge) = consumed.ridges.get(&(fi, v, w)) {
                ring.push(rings.point(*ridge));
            }
        }
        rings.push(ring, Origin::Face(fi), Some(face.normal));
    }
    // The strips, each by its law.
    for key in &keys {
        let seat = &seats[key];
        let Some(section) = by_edge.get(&seat.edge) else {
            continue;
        };
        let (a, b) = seat.ends;
        let corners = StripCorners {
            first: [at[&(seat.faces[0], a)], at[&(seat.faces[0], b)]],
            second: [at[&(seat.faces[1], a)], at[&(seat.faces[1], b)]],
        };
        match sections::law(section).strip(seat, corners)? {
            Strip::Planar(ring) => {
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
        match vertices::method(classes[&v], &met)? {
            vertices::VertexMethod::Facet => {
                let ring = order.iter().map(|fi| at[&(*fi, v)]).collect();
                rings.push(ring, Origin::Vertex(v), None);
            }
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
    let (built, faces) =
        ogeom_algo::make_faceted_solid(model, &assembled.points, &assembled.rings, &planes, tol)?;
    let mut history = History::new();
    let mut gone: Vec<Id> = blended.iter().copied().collect();
    gone.sort_unstable();
    for e in gone {
        history.delete(&planar.shapes[&e]);
    }
    let mut eaten: Vec<usize> = consumed.faces.iter().copied().collect();
    eaten.sort_unstable();
    for fi in eaten {
        history.delete(&planar.faces[fi].shape);
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
    // What a consumed band leaves nothing of: its sharp edges, and its
    // vertices where no facet survives them.
    for e in &consumed.gone {
        let shape = &planar.shapes[e];
        if history.generated(shape).is_empty() {
            history.delete(shape);
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
    // P5.
    guards::solid_fits(model, &built.shape, tol)?;
    Ok(Some(Built::new(built.shape, history)))
}
