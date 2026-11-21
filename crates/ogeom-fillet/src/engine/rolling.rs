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
/// meet closed by the ball rolling through it, a round reaching other
/// geometry handled as `overflow` says.
///
/// # Errors
///
/// As [`crate::fillet_edges_overflow`] where the rounds at every vertex
/// are of one sense, and by name with trim through off where the mixed
/// corners' rounds change another face. Elsewhere
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction),
/// naming the vertex, where one mixing convex and concave rounds joins
/// other than three edges, touches a curved face where its lone band's
/// sweep is no torus ([`lone_bands_first`]), or where the engine's rounds
/// close it on no sphere or torus; and where the request lies outside the
/// engine's rounds between planes.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    overflow: crate::Overflow,
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
        return crate::fillet_edges_overflow(model, solid, edges, radius, overflow, tol);
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
        return crate::fillet_edges_overflow(model, solid, edges, radius, overflow, tol);
    }
    // A mixed vertex on a curved face: the ball's sweep round the lone
    // band's end is the torus the lone band's foot is rounded by, so the
    // lone bands go first and the rest roll round their feet.
    let mut curved = false;
    for (_, vertex) in &mixed {
        for face in ancestors_of(model, &work, vertex, ShapeType::Face)? {
            curved |= !is_plane(model, &face);
        }
    }
    if curved {
        let built = lone_bands_first(model, &work, &traced, &mixed, radius, overflow, tol)?;
        return Ok(if placed {
            Built::new(built.shape, baked.history.then(&built.history))
        } else {
            built
        });
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
    let built = crate::fillet::trimmed_through(model, &work, &traced, Ok(built), overflow, tol)?;
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

/// Whether `face` lies on a plane.
fn is_plane(model: &Model, face: &Shape) -> bool {
    model
        .node(face)
        .and_then(|n| n.data().as_face())
        .and_then(|d| model.geometry().surface(d.surface))
        .is_some_and(|s| matches!(s, SurfaceGeometry::Plane(_)))
}

/// The ends of `edge`: none for an edge without vertices, one for a loop.
fn ends(model: &Model, edge: &Shape) -> OgeomResult<Vec<Shape>> {
    Ok(match ogeom_algo::edge_vertices(model, edge)? {
        Some((a, b)) if a.is_same(&b) => vec![a],
        Some((a, b)) => vec![a, b],
        None => Vec::new(),
    })
}

/// `edges` of `solid` rounded at `radius` where a vertex mixing convex and
/// concave rounds touches a curved face (ADR 0001, amendment "mixed
/// rolling corners on curved faces").
///
/// At such a vertex of three edges the lone round's edge is a line, and
/// each face beside it a plane, a cylinder ruled along it or a cone it is
/// a ruling of. Square to a plane cap (the face across from it), its band
/// is a cylinder of radius `radius` about an axis square to the cap, and
/// its ball's sweep round the band's end is the torus about that axis,
/// major radius twice the radius: the round of the band's foot on the
/// cap, which runs on tangent into the other two rounds. Leaning to the
/// cap, the band runs past the cap and is clipped there, and the ball
/// rolling round its foot changes walls where its contact crosses the
/// band's rails: the round between is the foot's pipe (amendment "a ball
/// changing walls off the vertex"). On a curved cap (a rib under a drum's
/// top) the band runs into the cap, its foot the curve they meet on, and
/// the ball changes walls the same way. So the lone bands are built first,
/// with every round reached from them through vertices of one sense, and
/// the rest of the request then rounds with each lone band's foot on its
/// cap as one tangent chain.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction),
/// naming the vertex, where a mixed vertex is not so: it joins other than
/// three edges, a convex and a concave round meet beside an edge left
/// sharp, the lone edge is curved, a face beside it is not ruled along
/// it, or a round meeting the lone round there is reached from another
/// lone band first. As [`crate::fillet_edges_overflow`] for either pass.
#[allow(
    clippy::too_many_lines,
    reason = "one rule, read and applied end to end"
)]
fn lone_bands_first(
    model: &mut Model,
    work: &Shape,
    traced: &[Shape],
    mixed: &[&(Point, Shape)],
    radius: f64,
    overflow: crate::Overflow,
    tol: Tolerances,
) -> OgeomResult<Built> {
    use ogeom_geom::Curve3d as _;
    let conf = tol.confusion();
    let named = |p: Point, why: &str| -> ogeom_core::OgeomError {
        ogeom_core::ogeom_err!(
            Construction,
            "a rolling corner does not close the vertex at ({:.6}, {:.6}, {:.6}): {why}",
            p.x,
            p.y,
            p.z
        )
    };
    let mut senses = Vec::with_capacity(traced.len());
    let mut traced_ends = Vec::with_capacity(traced.len());
    for edge in traced {
        senses.push(crate::setback::edge_is_convex(
            model, work, edge, radius, tol,
        )?);
        traced_ends.push(ends(model, edge)?);
    }
    let at = |vertex: &Shape| -> Vec<usize> {
        (0..traced.len())
            .filter(|&i| traced_ends[i].iter().any(|v| v.is_same(vertex)))
            .collect()
    };
    // Each mixed vertex: its lone edge, the request's edges there, its
    // point and the cap's surface.
    let mut lones: Vec<(usize, Vec<usize>, Point, SurfaceGeometry)> = Vec::new();
    for (p, vertex) in mixed.iter().map(|m| (m.0, &m.1)) {
        let mut around: Vec<Shape> = Vec::new();
        for e in ancestors_of(model, work, vertex, ShapeType::Edge)? {
            let degenerate = model
                .node(&e)
                .and_then(|n| n.data().as_edge())
                .is_some_and(|d| d.degenerate);
            if !degenerate && !around.iter().any(|a| a.is_same(&e)) {
                around.push(e);
            }
        }
        if around.len() != 3 {
            return Err(named(
                p,
                &format!(
                    "it joins {} edges, and a rolling corner closes a vertex of three",
                    around.len()
                ),
            ));
        }
        let asked = at(vertex);
        if asked.len() != 3 {
            return Err(named(
                p,
                "a convex and a concave round meet there beside an edge left sharp",
            ));
        }
        let mut convex = Vec::new();
        let mut concave = Vec::new();
        for &i in &asked {
            match senses[i] {
                Some(true) => convex.push(i),
                Some(false) => concave.push(i),
                None => return Err(named(p, "a round's side there is not read")),
            }
        }
        let lone = if convex.len() == 1 {
            convex[0]
        } else {
            concave[0]
        };
        let mut faces: Vec<Shape> = Vec::new();
        for f in ancestors_of(model, work, vertex, ShapeType::Face)? {
            if !faces.iter().any(|g| g.is_same(&f)) {
                faces.push(f);
            }
        }
        let mut beside = Vec::new();
        let mut caps = Vec::new();
        for f in faces {
            let holds = ogeom_topo::explore_unique(model, &f, ShapeType::Edge)?
                .iter()
                .any(|e| e.is_same(&traced[lone]));
            if holds {
                beside.push(f);
            } else {
                caps.push(f);
            }
        }
        let ([first, second], [cap]) = (beside.as_slice(), caps.as_slice()) else {
            return Err(named(
                p,
                "its faces do not stand two beside the lone round and one across",
            ));
        };
        let surface = |f: &Shape| -> OgeomResult<Option<SurfaceGeometry>> {
            Ok(super::rolled::host(model, f, tol)?.map(|(s, _)| s))
        };
        // A curved cap's corner is the foot's pipe too.
        let Some(cap_surface) = surface(cap)? else {
            return Err(named(
                p,
                "the face across from the lone round has no surface",
            ));
        };
        let curved_edge = || {
            named(
                p,
                "the lone round's edge is curved, and the ball's sweep round its band's end \
                 is no torus",
            )
        };
        let (curve, (t0, t1)) =
            ogeom_algo::placed_edge_curve(model, &traced[lone], tol)?.ok_or_else(curved_edge)?;
        if !matches!(curve, ogeom_geom::Curve::Line(_)) {
            return Err(curved_edge());
        }
        let (a, b) = (curve.point_at(t0, tol)?, curve.point_at(t1, tol)?);
        let along = (b - a) / (b - a).magnitude();
        // A lone edge leaning to the cap is rounded too: its foot's ball
        // changes walls off the vertex, where the foot's pipe takes it
        // (amendment "a ball changing walls off the vertex").
        for f in [first, second] {
            let ruled = match surface(f)? {
                Some(SurfaceGeometry::Plane(s)) => {
                    s.plane().normal().vector().dot(along).abs() <= tol.angular()
                }
                Some(SurfaceGeometry::Cylinder(s)) => {
                    s.cylinder()
                        .axis()
                        .direction
                        .vector()
                        .cross(along)
                        .magnitude()
                        <= tol.angular()
                }
                // A cone along one of its rulings.
                Some(cone @ SurfaceGeometry::Cone(_)) => {
                    let mut on = true;
                    for k in 0..=2 {
                        let q = a + (b - a) * (f64::from(k) / 2.0);
                        on &= ogeom_geom::project_on_surface(&cone, q, 16, tol)?.distance
                            <= conf * 10.0;
                    }
                    on
                }
                _ => false,
            };
            if !ruled {
                return Err(named(
                    p,
                    "a face beside the lone round is neither a plane, a cylinder along it nor \
                     a cone it is a ruling of, so its band is not one the engine rounds there",
                ));
            }
        }
        lones.push((lone, asked, p, cap_surface));
    }
    // The first pass: the lone bands, and every round reached from them
    // through vertices of one sense, whose corners are the default's.
    let is_mixed = |v: &Shape| mixed.iter().any(|m| m.1.is_same(v));
    let mut first = vec![false; traced.len()];
    for (lone, ..) in &lones {
        first[*lone] = true;
    }
    loop {
        let mut grew = false;
        for i in 0..traced.len() {
            if !first[i] {
                continue;
            }
            for v in &traced_ends[i] {
                if is_mixed(v) {
                    continue;
                }
                for j in at(v) {
                    if !first[j] {
                        first[j] = true;
                        grew = true;
                    }
                }
            }
        }
        if !grew {
            break;
        }
    }
    for (lone, asked, p, _) in &lones {
        if asked.iter().any(|&j| j != *lone && first[j]) {
            return Err(named(
                *p,
                "a round meeting the lone round there is reached from another lone band \
                 through corners of one sense, and would be rounded before it",
            ));
        }
    }
    let ones: Vec<Shape> = (0..traced.len())
        .filter(|&i| first[i])
        .map(|i| traced[i].clone())
        .collect();
    let built = crate::fillet_edges_overflow(model, work, &ones, radius, overflow, tol)?;
    let left = ogeom_topo::explore_unique(model, &built.shape, ShapeType::Edge)?;
    let kept = |e: &Shape| left.iter().any(|l| l.is_same(e));
    let mut second: Vec<Shape> = Vec::new();
    for i in (0..traced.len()).filter(|&i| !first[i]) {
        let images: Vec<Shape> = built
            .history
            .trace(&traced[i])
            .iter()
            .filter(|e| kept(e))
            .cloned()
            .collect();
        if images.is_empty() {
            ogeom_bail!(
                Construction,
                "an edge of the request is gone once the lone bands are rounded"
            );
        }
        for e in images {
            if !second.iter().any(|s| s.is_same(&e)) {
                second.push(e);
            }
        }
    }
    // Each lone band's foot: its edge on the cap.
    for (lone, _, p, cap) in &lones {
        let mut found = false;
        for face in built.history.generated(&traced[*lone]) {
            if model.kind_of(face)? != ShapeType::Face {
                continue;
            }
            for e in ogeom_topo::explore_unique(model, face, ShapeType::Edge)? {
                if !kept(&e) {
                    continue;
                }
                let Some((curve, (a, b))) = ogeom_algo::placed_edge_curve(model, &e, tol)? else {
                    continue;
                };
                let mut on_cap = true;
                for k in 0..=4 {
                    let q = curve.point_at(a + (b - a) * f64::from(k) / 4.0, tol)?;
                    on_cap &= ogeom_geom::project_on_surface(cap, q, 16, tol)?.distance <= conf;
                }
                if on_cap {
                    found = true;
                    if !second.iter().any(|s| s.is_same(&e)) {
                        second.push(e);
                    }
                }
            }
        }
        if !found {
            return Err(named(
                *p,
                "the lone band leaves no foot on the face across from it",
            ));
        }
    }
    let rolled = crate::fillet_edges_overflow(model, &built.shape, &second, radius, overflow, tol)?;
    Ok(Built::new(
        rolled.shape,
        built.history.then(&rolled.history),
    ))
}
