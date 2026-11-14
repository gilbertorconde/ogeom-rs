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
mod fallback;
mod guards;
mod request;
mod rims;
mod round_strips;
mod rounds;
mod seams;
mod seats;
mod sections;
mod setbacks;
mod trimmed;
mod vertices;

pub(crate) use fallback::{Sense, Tool, applied as apply_tools, straight_non_local};
pub(crate) use request::{Request, Section};
pub(crate) use rims::assemble as round_rims;
pub(crate) use round_strips::{non_local_tools as round_non_local_tools, round as round_curved};
pub(crate) use seams::{Dissolved, Seamless, dissolve, dissolve_curved, is_seam};
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

/// What makes a solid bounded by planes degenerate, where it is: an edge
/// bounded by other than two faces (two blocks touching along an edge), an
/// edge closing on itself, two edges between one pair of vertices, or a
/// vertex whose faces close round it more than once (two blocks sharing
/// only a corner in one shell). The engine refuses these by name (ADR
/// 0001, amendment "what M6 refuses by name"); `None` for a solid with a
/// curved face, or with none of them.
///
/// # Errors
///
/// As reading the solid.
pub(crate) fn degenerate_planar(model: &Model, solid: &Shape) -> OgeomResult<Option<String>> {
    use ogeom_topo::{Filter, ShapeType, explore, explore_unique};
    let faces = explore_unique(model, solid, ShapeType::Face)?;
    for f in &faces {
        let planar = model
            .node(f)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Plane(_)));
        if !planar {
            return Ok(None);
        }
    }
    let mut uses: FastMap<Id, Vec<usize>> = FastMap::default();
    for (fi, f) in faces.iter().enumerate() {
        for e in explore(model, f, Filter::OfType(ShapeType::Edge))? {
            uses.entry(id(&e)).or_default().push(fi);
        }
    }
    let mut pairs: FastSet<(Id, Id)> = FastSet::default();
    let mut at: FastMap<Id, Vec<Id>> = FastMap::default();
    for e in explore_unique(model, solid, ShapeType::Edge)? {
        let n = uses.get(&id(&e)).map_or(0, Vec::len);
        if n != 2 {
            return Ok(Some(format!("an edge is bounded by {n} faces, not two")));
        }
        let Some((a, b)) = ogeom_algo::edge_vertices(model, &e)? else {
            continue;
        };
        if a.is_same(&b) {
            return Ok(Some("an edge closes on itself".into()));
        }
        let (x, y) = (id(&a), id(&b));
        if !pairs.insert((x.min(y), x.max(y))) {
            return Ok(Some("two edges join one pair of vertices".into()));
        }
        at.entry(x).or_default().push(id(&e));
        at.entry(y).or_default().push(id(&e));
    }
    // Round each vertex its faces close once: joined face to face across
    // its edges, they make one fan, not two cones touching at a point.
    let mut vertices: Vec<(&Id, &Vec<Id>)> = at.iter().collect();
    vertices.sort_unstable_by_key(|(v, _)| **v);
    for (_, edges) in vertices {
        let mut root: FastMap<usize, usize> = FastMap::default();
        fn find(root: &mut FastMap<usize, usize>, f: usize) -> usize {
            let up = *root.entry(f).or_insert(f);
            if up == f {
                return f;
            }
            let top = find(root, up);
            root.insert(f, top);
            top
        }
        for e in edges {
            if let [a, b] = uses[e].as_slice() {
                let (ra, rb) = (find(&mut root, *a), find(&mut root, *b));
                root.insert(ra, rb);
            }
        }
        let faces: Vec<usize> = root.keys().copied().collect();
        let mut tops: Vec<usize> = faces.into_iter().map(|f| find(&mut root, f)).collect();
        tops.sort_unstable();
        tops.dedup();
        if tops.len() > 1 {
            return Ok(Some(
                "a vertex's faces close round it more than once".into(),
            ));
        }
    }
    Ok(None)
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
    let Some(mut planar) = classify::read_planar(model, solid, tol)? else {
        return Ok(None);
    };
    for split in &request.splits {
        if !planar.split(model, split)? {
            return Ok(None);
        }
    }
    let planar = planar;
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
    let mut steps = FastMap::default();
    let Some(corners) = setbacks::setback_corners(&planar, &setback, &mut steps, tol)? else {
        return Ok(None);
    };
    // A round reaching or touching a face not round its ends: the non-local
    // test, then the one-boolean fallback (ADR 0001, M5c); where that
    // declines, the sheet as before.
    if request
        .blends
        .iter()
        .all(|b| matches!(b.section, Section::Round { .. }))
        && fallback::non_local(&planar, &seats, &blended, &corners, &steps, tol).is_some()
        && let Some(built) = fallback::rounds_one_boolean(
            model, solid, &planar, &seats, &by_edge, &blended, &corners, &steps, tol,
        )?
    {
        return Ok(Some(built));
    }
    // Rounds are curved: their strips and vertex blends go into a sheet.
    if request
        .blends
        .iter()
        .any(|b| matches!(b.section, Section::Round { .. }))
    {
        return rounds::assemble(model, &planar, &seats, &by_edge, &corners, &steps, tol);
    }
    // A bevel reaching or touching a face not round its ends changes a face
    // the request does not name, which the rings do not trim: the
    // non-local test, then the one-boolean fallback (ADR 0001, M5a).
    // A bevel rolling over its face's far side notches the faces beyond:
    // that is the rings' (amendment "M5b"), not a non-local blend.
    if !request.mitre
        && !consumed::rolls_over(&planar, &seats, &blended, &corners, tol)
        && fallback::non_local(&planar, &seats, &blended, &corners, &steps, tol).is_some()
    {
        return fallback::one_boolean(
            model, solid, &planar, &seats, &blended, &corners, &steps, tol,
        );
    }
    // Chamfers set back differently along one line step there: only where
    // the line's two faces alone meet, at a joint between two pieces of one
    // edge. Elsewhere (a third face at the vertex) they are the current
    // path's.
    let mut joints: Vec<Id> = steps.keys().map(|(_, v)| *v).collect();
    joints.sort_unstable();
    joints.dedup();
    for v in &joints {
        if classify::faces_round(&planar.faces, *v).is_none_or(|round| round.len() != 2) {
            return fallback::one_boolean(
                model, solid, &planar, &seats, &blended, &corners, &steps, tol,
            );
        }
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
        consumed::Consumed::Boolean if !request.mitre => {
            return fallback::one_boolean(
                model, solid, &planar, &seats, &blended, &corners, &steps, tol,
            );
        }
        consumed::Consumed::Boolean => return Ok(None),
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
            .map(|corners| {
                corners
                    .iter()
                    .flat_map(|v| match steps.get(&(fi, *v)) {
                        Some([p, q]) => vec![rings.point(*p), rings.point(*q)],
                        None => vec![at[&(fi, *v)]],
                    })
                    .collect()
            })
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
        // At a step, the corner on this edge's own line: the face's first
        // where it arrives along the edge, its second where it leaves.
        let mut end = |face: usize, v: Id, sign: f64| -> usize {
            if let Some([p, q]) = steps.get(&(face, v)) {
                let arriving = (face == seat.faces[0]) == (v == b);
                return rings.point(if arriving { *p } else { *q });
            }
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
    // Each joint's step: what the deeper bevel takes there past the
    // shallower, in the plane through the joint square to the line. Each
    // bevel's end section is the triangle of the joint and its corners on
    // the two faces, so the step is the larger less the smaller, facing the
    // deeper side.
    for v in &joints {
        let Some(round) = classify::faces_round(&planar.faces, *v) else {
            return Ok(None);
        };
        let (fa, fb) = (round[0], round[1]);
        let pair = |fi: usize| {
            steps
                .get(&(fi, *v))
                .copied()
                .unwrap_or([corners[&(fi, *v)]; 2])
        };
        let (pa, pb) = (pair(fa), pair(fb));
        let at = planar.points[v];
        // `fa` arrives along the line from one side and leaves to the
        // other; `fb` walks it the other way round.
        let (back, on) = ((pa[0] - at).magnitude(), (pa[1] - at).magnitude());
        let (on_b, back_b) = ((pb[0] - at).magnitude(), (pb[1] - at).magnitude());
        let eps = tol.confusion();
        let ahead = on >= back - eps && on_b >= back_b - eps;
        let behind = on <= back + eps && on_b <= back_b + eps;
        if ahead && behind {
            continue;
        }
        let Some((_, next)) = planar.faces[fa].beside(*v) else {
            return Ok(None);
        };
        let toward = planar.points[&next] - at;
        // Each piece of the step wound to face its way.
        let mut put = |mut ring: Vec<ogeom_math::Point>, facing: ogeom_math::Vector| {
            let turn = (1..ring.len() - 1).fold(ogeom_math::Vector::ZERO, |n, k| {
                n + (ring[k] - ring[0]).cross(ring[k + 1] - ring[0])
            });
            if turn.dot(facing) < 0.0 {
                ring.reverse();
            }
            let ring: Vec<usize> = ring.iter().map(|p| rings.point(*p)).collect();
            rings.push(ring, Origin::Vertex(*v), Some(facing));
        };
        if ahead || behind {
            put(
                vec![pa[0], pa[1], pb[0], pb[1]],
                if ahead { toward } else { -toward },
            );
            continue;
        }
        // The two end sections cross: one bevel the deeper on `fa`, the
        // other on `fb`. Their far sides cross at one point, and the step
        // is a triangle on each face's side of it, facing the side whose
        // bevel is deeper there.
        let (u, w) = (pa[1] - pb[0], pa[0] - pb[1]);
        let n = u.cross(w);
        if n.magnitude() <= eps * eps {
            return Ok(None);
        }
        let t = (pa[0] - pb[0]).cross(w).dot(n) / n.dot(n);
        let cross = pb[0] + u * t;
        put(
            vec![pa[0], pa[1], cross],
            if on > back { toward } else { -toward },
        );
        put(
            vec![pb[0], pb[1], cross],
            if on_b > back_b { toward } else { -toward },
        );
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
        let mut eaten: Vec<usize> = cluster.faces.iter().copied().collect();
        eaten.sort_unstable();
        for fi in &eaten {
            let face = &planar.faces[*fi];
            if face.holes.is_empty() {
                continue;
            }
            let holes: Vec<Vec<ogeom_math::Point>> = face
                .holes
                .iter()
                .map(|h| h.iter().map(|v| planar.points[v]).collect())
                .collect();
            if !rings.drop_in(*fi, face.normal, &holes) {
                return Ok(None);
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
