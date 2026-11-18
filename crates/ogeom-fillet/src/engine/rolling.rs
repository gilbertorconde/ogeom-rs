//! Rolling corners (ADR 0001, amendment "rolling corners"): rounds of one
//! radius whose vertices close by the ball rolling through them, every
//! corner face a piece of the ball's sweep (a sphere or a torus), never a
//! patch, and a vertex that is not so refused by name.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::SurfaceGeometry;
use ogeom_math::Point;
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of};

use super::Request;

/// `edges` of `solid` rounded at `radius`, each vertex where two or more
/// meet closed by the ball rolling through it.
///
/// # Errors
///
/// As [`crate::fillet_edges`] where the rounds at every vertex are of one
/// sense. Elsewhere
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction),
/// naming the vertex, where one mixing convex and concave rounds touches a
/// curved face or joins other than three edges, or where the engine's
/// rounds close it on no sphere or torus; and where the request lies
/// outside the engine's rounds between planes.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !(radius.is_finite() && radius > tol.confusion()) {
        ogeom_bail!(Construction, "a fillet of radius {radius} rounds nothing");
    }
    let mut once: Vec<Shape> = Vec::with_capacity(edges.len());
    for edge in edges {
        if !once.iter().any(|e| e.is_same(edge)) {
            once.push(edge.clone());
        }
    }
    let baked = ogeom_algo::baked_where_placed(model, solid, tol)?;
    let placed = !baked.shape.is_same(solid);
    let mut traced = Vec::with_capacity(once.len());
    for edge in &once {
        if !placed {
            traced.push(edge.clone());
            continue;
        }
        let [only] = baked.history.trace(edge) else {
            ogeom_bail!(
                Construction,
                "an edge of the request is not one of the solid's"
            );
        };
        traced.push(only.clone());
    }
    let work = baked.shape.clone();
    let meetings = meetings(model, &traced)?;
    // No vertex where rounds meet: no corner, and the rounds are the
    // default's.
    if meetings.is_empty() {
        return crate::fillet_edges(model, solid, edges, radius, tol);
    }
    // Where the rounds meeting at a vertex are all of one sense, the
    // default's corner is already the ball's: the band running on through
    // a tangent junction, the sphere where one ball touches every face,
    // the bands crossing beside an edge left sharp, the apex's envelope of
    // balls. Only a vertex mixing convex and concave rounds needs the
    // ball's sweep round the lone band's end, built here between planes.
    let mut mixed: Vec<&(Point, Shape)> = Vec::new();
    for meeting in &meetings {
        let (mut convex, mut concave) = (false, false);
        for edge in &traced {
            let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
                continue;
            };
            if !(a.is_same(&meeting.1) || b.is_same(&meeting.1)) {
                continue;
            }
            match crate::setback::edge_is_convex(model, &work, edge, radius, tol)? {
                Some(true) => convex = true,
                Some(false) => concave = true,
                None => (convex, concave) = (true, true),
            }
        }
        if convex && concave {
            mixed.push(meeting);
        }
    }
    let between_planes = meetings
        .iter()
        .map(|(_, vertex)| between_planes(model, &work, vertex))
        .collect::<OgeomResult<Vec<bool>>>()?
        .into_iter()
        .all(|planar| planar);
    if mixed.is_empty() && !between_planes {
        return crate::fillet_edges(model, solid, edges, radius, tol);
    }
    for (p, vertex) in mixed {
        let named = |why: String| -> ogeom_core::OgeomError {
            ogeom_core::ogeom_err!(
                Construction,
                "a rolling corner does not close the vertex at ({:.6}, {:.6}, {:.6}): {why}",
                p.x,
                p.y,
                p.z
            )
        };
        for face in ancestors_of(model, &work, vertex, ShapeType::Face)? {
            let planar = model
                .node(&face)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(|s| matches!(s, SurfaceGeometry::Plane(_)));
            if !planar {
                return Err(named(
                    "a face there is curved, and a rolling corner closes only between planes"
                        .into(),
                ));
            }
        }
        let mut around: Vec<Shape> = Vec::new();
        for e in ancestors_of(model, &work, vertex, ShapeType::Edge)? {
            let degenerate = model
                .node(&e)
                .and_then(|n| n.data().as_edge())
                .is_some_and(|d| d.degenerate);
            if !degenerate && !around.iter().any(|a| a.is_same(&e)) {
                around.push(e);
            }
        }
        if around.len() != 3 {
            return Err(named(format!(
                "it joins {} edges, and a rolling corner closes a vertex of three",
                around.len()
            )));
        }
    }
    let rounds: Vec<(Shape, f64)> = traced.iter().map(|e| (e.clone(), radius)).collect();
    let mut request = Request::from_rounds(&rounds, tol)?;
    request.rolling = true;
    let Some(built) = super::run(model, &work, &request, tol)? else {
        ogeom_bail!(
            Construction,
            "a rolling corner is built where every round of the request stands between \
             planes and the rounds close among themselves; this request rounds a solid \
             with a curved face or a face with a hole, or a round that reaches a face \
             beyond its ends, and is refused"
        );
    };
    Ok(if placed {
        Built::new(built.shape, baked.history.then(&built.history))
    } else {
        built
    })
}

/// Whether every face at `vertex` of `solid` is a plane and the vertex
/// joins three edges: the corners the planar rounds close.
fn between_planes(model: &Model, solid: &Shape, vertex: &Shape) -> OgeomResult<bool> {
    for face in ancestors_of(model, solid, vertex, ShapeType::Face)? {
        let planar = model
            .node(&face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, SurfaceGeometry::Plane(_)));
        if !planar {
            return Ok(false);
        }
    }
    let mut around: Vec<Shape> = Vec::new();
    for e in ancestors_of(model, solid, vertex, ShapeType::Edge)? {
        let degenerate = model
            .node(&e)
            .and_then(|n| n.data().as_edge())
            .is_some_and(|d| d.degenerate);
        if !degenerate && !around.iter().any(|a| a.is_same(&e)) {
            around.push(e);
        }
    }
    Ok(around.len() == 3)
}

/// The vertices where two or more of `edges` meet, with their
/// points, in the order of their points.
fn meetings(model: &Model, edges: &[Shape]) -> OgeomResult<Vec<(Point, Shape)>> {
    let mut ends: Vec<(Shape, usize)> = Vec::new();
    for edge in edges {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
            continue;
        };
        let pair = if a.is_same(&b) { vec![a] } else { vec![a, b] };
        for v in pair {
            match ends.iter_mut().find(|(w, _)| w.is_same(&v)) {
                Some((_, n)) => *n += 1,
                None => ends.push((v, 1)),
            }
        }
    }
    let mut out: Vec<(Point, Shape)> = Vec::new();
    for (v, n) in ends {
        if n < 2 {
            continue;
        }
        let Some(p) = model
            .node(&v)
            .and_then(|node| node.data().as_vertex().map(|d| d.point))
        else {
            ogeom_bail!(Construction, "a vertex of the request has no point");
        };
        out.push((p, v));
    }
    out.sort_by(|(p, _), (q, _)| {
        p.x.total_cmp(&q.x)
            .then(p.y.total_cmp(&q.y))
            .then(p.z.total_cmp(&q.z))
    });
    Ok(out)
}
