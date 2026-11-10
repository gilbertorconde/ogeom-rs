//! Setback corners: the bands meeting at a vertex stopped a chosen
//! distance along their edges, and one N-sided patch over the opening.
//!
//! The bands are rounded as [`fillet_edges`](crate::fillet_edges) rounds
//! them, meeting in the ball's own corner. A cell is then cut from the
//! corner: the polyhedron bounded by the planes through the vertex and,
//! for each edge, the plane square to it at its setback. What is left of
//! each band ends on its section plane in an exact arc. On each face the
//! opening is closed by a conic from one band's rail end to the next one's,
//! tangent to both rails, so the face's trim turns smoothly where the
//! bands leave it. The opening is then a loop of arcs and conics, filled
//! by one patch meeting each band and each face at the continuity asked.

use ogeom_algo::{Built, History};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Continuity, Curve, Curve3d as _};
use ogeom_math::{Direction, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use crate::fill_n::{FillBoundary, make_filling_n, make_filling_n_within};

/// How far each band meeting a vertex stops from it, and how the patch
/// between them meets the bands and the faces.
#[derive(Debug, Clone)]
pub struct Setbacks {
    /// Distance from the vertex along each edge at which that edge's band
    /// stops; edges not listed take `default`.
    pub per_edge: Vec<(Shape, f64)>,
    /// The setback of every edge not listed.
    pub default: f64,
    /// How the patch meets each band and face: [`Continuity::G1`]
    /// (tangent). [`Continuity::G2`] is refused while bands are circular,
    /// G1 to their faces.
    pub continuity: Continuity,
}

/// Round `vertex` of `solid` with the bands of its edges at `radius`, each
/// stopped at its setback, and one patch over the opening between them.
///
/// Each band is cut square to its edge at its setback, measured from the
/// vertex along the edge. On each face between two bands the trim runs
/// from one band's rail end to the next one's on a conic tangent to both
/// rails. The opening those arcs and conics bound is filled by one patch
/// meeting each band end and each face at `setbacks.continuity`, fitted
/// to the bands' own fit target. A setback equal to the ball's own (the
/// distance along the edge to the foot of the ball's centre, the radius at
/// a square corner) on every edge is [`fillet_edges`](crate::fillet_edges)'
/// corner.
///
/// What `solid` shares with other shapes is copied before it is worked on
/// ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// `radius` or a setback is not a positive length, a setback is shorter
/// than the ball's own along its edge or reaches past the edge's far end,
/// the continuity is neither G1 nor G2, the vertex is not a convex planar
/// corner one ball touches every face of (a curved face meeting it, or
/// planes sharing no tangent ball, are still owed), or the opening does not
/// close as a loop of the bands' arcs and the faces' conics;
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where the
/// patch misses its target.
pub fn round_vertex_setback(
    model: &mut Model,
    solid: &Shape,
    vertex: &Shape,
    radius: f64,
    setbacks: &Setbacks,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut own: Vec<Shape> = vec![vertex.clone()];
    own.extend(setbacks.per_edge.iter().map(|(e, _)| e.clone()));
    ogeom_algo::on_own_nodes(model, solid, &own, |model, solid, parts| {
        let named: Vec<(Shape, f64)> = parts[1..]
            .iter()
            .zip(&setbacks.per_edge)
            .map(|(e, (_, s))| (e.clone(), *s))
            .collect();
        let own = Setbacks {
            per_edge: named,
            default: setbacks.default,
            continuity: setbacks.continuity,
        };
        setback_own(model, solid, &parts[0], radius, &own, tol)
    })
}

/// The point a vertex stands at.
pub(crate) fn point_of(model: &Model, vertex: &Shape) -> OgeomResult<Point> {
    let Some(raw) = model
        .node(vertex)
        .and_then(|n| n.data().as_vertex().map(|d| d.point))
    else {
        ogeom_bail!(Construction, "the vertex holds no point");
    };
    Ok(vertex.transform(model.datums())?.apply(raw))
}

/// A setback corner read off the solid before its edges are rounded: the
/// corner's frame, the solid's edge along each of its edge directions with
/// that edge's length, and each edge's setback and the ball's own setback
/// along it.
struct SetbackCorner {
    frame: crate::corner::CornerFrame,
    edges: Vec<Shape>,
    reach: Vec<f64>,
    setback: Vec<f64>,
    own_setback: Vec<f64>,
}

impl SetbackCorner {
    /// Whether every setback is the ball's own, the corner `fillet_edges`
    /// closes.
    fn is_the_balls(&self, tol: Tolerances) -> bool {
        self.setback
            .iter()
            .zip(&self.own_setback)
            .all(|(s, own)| *s <= own + tol.confusion())
    }
}

/// The setback corner at `vertex` of `solid`, its edges not yet rounded,
/// checked as [`round_vertex_setback`] states.
fn read_corner(
    model: &Model,
    solid: &Shape,
    vertex: &Shape,
    radius: f64,
    setbacks: &Setbacks,
    tol: Tolerances,
) -> OgeomResult<SetbackCorner> {
    if !radius.is_finite() || radius <= tol.confusion() {
        ogeom_bail!(Construction, "a blend radius must be a positive distance");
    }
    if setbacks.continuity == Continuity::G2 {
        // Where a band's end arc meets a face, the arc (radius r, tangent
        // to the face) bends the patch at 1/r along itself, and the face
        // asks it not to bend across the face's trim, which runs square to
        // the arc there: no patch is G2 to both while the band itself is
        // only G1 to the face.
        ogeom_bail!(
            Construction,
            "a G2 setback corner needs bands curvature continuous with their \
             faces; a circular band meets its face at G1, so its end arc \
             bends the patch at 1/radius where the face asks none; see \
             docs/CAPABILITIES.md, fillet.setback-corners"
        );
    }
    if setbacks.continuity != Continuity::G1 {
        ogeom_bail!(
            Construction,
            "a setback corner meets its bands and faces at G1, not {:?}",
            setbacks.continuity
        );
    }
    let frame = crate::corner::corner_frame(model, solid, vertex, radius, tol)?;
    let n = frame.d.len();

    let (edges, reach) = corner_edges(model, solid, vertex, &frame, tol)?;
    let setback: Vec<f64> = edges
        .iter()
        .map(|e| {
            setbacks
                .per_edge
                .iter()
                .find(|(named, _)| named.is_same(e))
                .map_or(setbacks.default, |(_, s)| *s)
        })
        .collect();
    let own_setback: Vec<f64> = frame
        .d
        .iter()
        .map(|dk| (frame.far - frame.corner).dot(dk.vector()))
        .collect();
    for k in 0..n {
        let s = setback[k];
        if !s.is_finite() || s <= tol.confusion() {
            ogeom_bail!(Construction, "a setback of {s} is not a positive distance");
        }
        if s < own_setback[k] - tol.confusion() {
            ogeom_bail!(
                Construction,
                "a setback of {s} stops the band inside the ball's own corner, \
                 which reaches {} along its edge",
                own_setback[k]
            );
        }
        if s >= reach[k] - tol.confusion() {
            ogeom_bail!(
                Construction,
                "a setback of {s} reaches past its edge, {} long",
                reach[k]
            );
        }
    }
    Ok(SetbackCorner {
        frame,
        edges,
        reach,
        setback,
        own_setback,
    })
}

/// [`round_vertex_setback`] on a shape every node below which is its own.
fn setback_own(
    model: &mut Model,
    solid: &Shape,
    vertex: &Shape,
    radius: f64,
    setbacks: &Setbacks,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let corner = read_corner(model, solid, vertex, radius, setbacks, tol)?;

    // The bands, meeting in the ball's corner: all that a setback equal to
    // the ball's own asks for.
    let banded = crate::fillet::fillet_edges(model, solid, &corner.edges, radius, tol)?;
    if corner.is_the_balls(tol) {
        return Ok(banded);
    }
    close_corner(
        model,
        solid,
        banded,
        vertex,
        &corner.frame,
        &corner.setback,
        setbacks.continuity,
        tol,
    )
}

/// Round `edges` of `solid` at `radius`, as [`fillet_edges`](crate::fillet_edges)
/// does, with each vertex of `corners` closed as a setback corner, as
/// [`round_vertex_setback`] closes it.
///
/// Every corner is read off the solid before its edges are rounded, the
/// edges are rounded once, and each corner is then cut back to its
/// setbacks and its opening filled. Vertices not named close in the
/// ball's own corner. A setback corner rounds every edge that meets its
/// vertex, so each of them must be in `edges`. Along an edge between two
/// named corners, the two setbacks must stop short of each other.
///
/// What `solid` shares with other shapes is copied before it is worked on
/// ([`ogeom_algo::on_own_nodes`]), which leaves them as they are.
///
/// # Errors
///
/// As [`round_vertex_setback`] for each corner and
/// [`fillet_edges`](crate::fillet_edges) for the edges;
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// where a vertex is named twice, an edge meeting a named vertex is not in
/// `edges`, or the setbacks at both ends of an edge meet.
pub fn fillet_edges_cornered(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    corners: &[(Shape, Setbacks)],
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut own: Vec<Shape> = edges.to_vec();
    own.extend(corners.iter().map(|(v, _)| v.clone()));
    for (_, setbacks) in corners {
        own.extend(setbacks.per_edge.iter().map(|(e, _)| e.clone()));
    }
    ogeom_algo::on_own_nodes(model, solid, &own, |model, solid, parts| {
        let (edges, rest) = parts.split_at(edges.len());
        let (vertices, mut named) = rest.split_at(corners.len());
        let mut own_corners = Vec::with_capacity(corners.len());
        for (vertex, (_, setbacks)) in vertices.iter().zip(corners) {
            let (these, others) = named.split_at(setbacks.per_edge.len());
            named = others;
            own_corners.push((
                vertex.clone(),
                Setbacks {
                    per_edge: these
                        .iter()
                        .zip(&setbacks.per_edge)
                        .map(|(e, (_, s))| (e.clone(), *s))
                        .collect(),
                    default: setbacks.default,
                    continuity: setbacks.continuity,
                },
            ));
        }
        cornered_own(model, solid, edges, radius, &own_corners, tol)
    })
}

/// [`fillet_edges_cornered`] on a shape every node below which is its own.
fn cornered_own(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    corners: &[(Shape, Setbacks)],
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut read: Vec<SetbackCorner> = Vec::with_capacity(corners.len());
    for (i, (vertex, setbacks)) in corners.iter().enumerate() {
        if corners[..i].iter().any(|(v, _)| v.is_same(vertex)) {
            ogeom_bail!(Construction, "a vertex is named twice as a setback corner");
        }
        let corner = read_corner(model, solid, vertex, radius, setbacks, tol)?;
        if corner
            .edges
            .iter()
            .any(|e| !edges.iter().any(|asked| asked.is_same(e)))
        {
            ogeom_bail!(
                Construction,
                "a setback corner rounds every edge meeting its vertex, and an edge \
                 meeting a named corner is not among the edges to round"
            );
        }
        read.push(corner);
    }
    // Two corners on one edge cut it back from both ends: the cells must
    // stop short of each other along it, each band's middle left whole.
    for (i, a) in read.iter().enumerate() {
        for b in &read[i + 1..] {
            for (ka, ea) in a.edges.iter().enumerate() {
                for (kb, eb) in b.edges.iter().enumerate() {
                    if ea.is_same(eb)
                        && a.setback[ka] + b.setback[kb] >= a.reach[ka] - tol.confusion()
                    {
                        ogeom_bail!(
                            Construction,
                            "the setbacks of {} and {} at both ends of an edge {} long \
                             meet; a band's middle must be left between two setback \
                             corners",
                            a.setback[ka],
                            b.setback[kb],
                            a.reach[ka]
                        );
                    }
                }
            }
        }
    }
    let mut built = crate::fillet::fillet_edges(model, solid, edges, radius, tol)?;
    for ((vertex, setbacks), corner) in corners.iter().zip(&read) {
        if corner.is_the_balls(tol) {
            continue;
        }
        built = close_corner(
            model,
            solid,
            built,
            vertex,
            &corner.frame,
            &corner.setback,
            setbacks.continuity,
            tol,
        )?;
    }
    Ok(built)
}

/// The solid's edge along each of `frame`'s edge directions from `vertex`,
/// and each one's length; refused by name where none does, or where the
/// nearest is curved.
pub(crate) fn corner_edges(
    model: &Model,
    solid: &Shape,
    vertex: &Shape,
    frame: &crate::corner::CornerFrame,
    tol: Tolerances,
) -> OgeomResult<(Vec<Shape>, Vec<f64>)> {
    let n = frame.d.len();
    let at_vertex: Vec<Shape> = explore_unique(model, solid, ShapeType::Edge)?
        .into_iter()
        .filter(|e| {
            ogeom_algo::edge_vertices(model, e)
                .ok()
                .flatten()
                .is_some_and(|(a, b)| a.is_same(vertex) || b.is_same(vertex))
        })
        .collect();
    let mut edges: Vec<Shape> = Vec::with_capacity(n);
    let mut reach: Vec<f64> = Vec::with_capacity(n);
    for dk in &frame.d {
        let mut best: Option<(f64, Shape, f64)> = None;
        for e in &at_vertex {
            let Some((a, b)) = ogeom_algo::edge_vertices(model, e)? else {
                continue;
            };
            let far = if a.is_same(vertex) { b } else { a };
            let along = point_of(model, &far)? - frame.corner;
            let length = along.magnitude();
            if length <= tol.confusion() {
                continue;
            }
            let cos = along.dot(dk.vector()) / length;
            if best.as_ref().is_none_or(|(c, _, _)| cos > *c) {
                best = Some((cos, e.clone(), length));
            }
        }
        let Some((cos, edge, length)) = best else {
            ogeom_bail!(
                Construction,
                "no edge of the solid leaves the vertex along a corner edge"
            );
        };
        if cos < 1.0 - tol.angular().sqrt() {
            ogeom_bail!(
                Construction,
                "a setback corner rounds straight edges; an edge leaving the \
                 vertex is curved, which is still owed; see \
                 docs/CAPABILITIES.md, fillet.setback-corners"
            );
        }
        edges.push(edge);
        reach.push(length);
    }
    Ok((edges, reach))
}

/// How a band at a setback corner ends where the patch meets it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CornerEnd {
    /// A round: its end arc on its section plane, which the patch meets at
    /// the corner's continuity.
    Round,
    /// A bevel: the straight end of its face on its section plane, which
    /// the patch meets in position only, as it meets the trims from that
    /// end on the faces either side: the bevel meets those faces in a
    /// crease, and the patch, running along the bevel's end, leaves each
    /// face at the bevel's angle there.
    Bevel,
}

/// Refused by name unless the cell cut from `banded` at `vertex` holds
/// nothing of `solid`, the solid the corner was asked on, but the corner
/// itself. What is foreign is a face of `solid` that does not meet the
/// vertex and lies on none of the corner's own planes through it
/// (`planes`, each a normal and a point): the bands, the bevels and the
/// ball's own corner are made by the rounding, on surfaces of their own.
/// Refused where the cut takes a face on a foreign surface whole (every
/// corner of it inside the cell, `walls`, each an inward normal and a
/// point), where a face on a foreign surface that the cut leaves has an
/// edge on a wall inside the cell (the cut passed through it), or where the
/// cut leaves fewer shells. Any of these (a void, a hole, another feature)
/// would be left out of what is sewn back, and the solid would close
/// without it.
#[allow(
    clippy::too_many_arguments,
    reason = "one cell and what it is read against"
)]
pub(crate) fn guard_cell(
    model: &Model,
    solid: &Shape,
    banded: &Shape,
    cut: &Shape,
    vertex: &Shape,
    planes: &[(Vector, Point)],
    walls: &[(Vector, Point)],
    tol: Tolerances,
) -> OgeomResult<()> {
    let near = tol.confusion() * 100.0;
    let surface_of = |face: &Shape| {
        model
            .node(face)
            .and_then(|n| n.data().as_face())
            .map(|d| d.surface)
    };
    let on_a_plane = |face: &Shape| -> OgeomResult<bool> {
        let planar = model
            .node(face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Plane(_)));
        if !planar {
            return Ok(false);
        }
        let (origin, outward) = ogeom_algo::face_normal(model, face, tol)?;
        Ok(planes.iter().any(|(normal, at)| {
            outward.cross(*normal).magnitude() <= tol.angular() * 10.0
                && normal.dot(origin - *at).abs() <= near
        }))
    };
    let mut foreign = Vec::new();
    for face in explore_unique(model, solid, ShapeType::Face)? {
        let meets = explore_unique(model, &face, ShapeType::Vertex)?
            .iter()
            .any(|v| v.is_same(vertex));
        if !meets
            && !on_a_plane(&face)?
            && let Some(surface) = surface_of(&face)
        {
            foreign.push(surface);
        }
    }
    let is_foreign = |face: &Shape| surface_of(face).is_some_and(|s| foreign.contains(&s));
    let inside = |p: Point| walls.iter().all(|(n, at)| n.dot(p - *at) >= -near);
    let refused = || {
        ogeom_core::ogeom_err!(
            Construction,
            "another feature of the solid lies inside the setback corner's cell (a face \
             that is neither one of the corner's faces nor made by its bands, or a shell \
             of its own), and cutting the cell would remove it; see \
             docs/CAPABILITIES.md, fillet.setback-corners"
        )
    };
    // A face the cut takes whole.
    for face in explore_unique(model, banded, ShapeType::Face)? {
        if !is_foreign(&face) {
            continue;
        }
        let mut all_in = true;
        for v in explore_unique(model, &face, ShapeType::Vertex)? {
            if !inside(point_of(model, &v)?) {
                all_in = false;
                break;
            }
        }
        if all_in {
            return Err(refused());
        }
    }
    // A face the cut passes through.
    for face in explore_unique(model, cut, ShapeType::Face)? {
        if !is_foreign(&face) {
            continue;
        }
        for e in explore_unique(model, &face, ShapeType::Edge)? {
            let Some((curve, (lo, hi))) = ogeom_algo::placed_edge_curve(model, &e, tol)? else {
                continue;
            };
            let mid = curve.point_at(f64::midpoint(lo, hi), tol)?;
            if walls.iter().any(|(n, at)| n.dot(mid - *at).abs() <= near) && inside(mid) {
                return Err(refused());
            }
        }
    }
    let shells = |shape: &Shape| explore_unique(model, shape, ShapeType::Shell).map(|s| s.len());
    if shells(cut)? < shells(banded)? {
        return Err(refused());
    }
    Ok(())
}

/// The bands of `banded` meeting at `vertex`, whose corner reads as
/// `frame`, stopped at `setback` along their edges and the opening between
/// them filled by one patch at `continuity`.
#[allow(clippy::too_many_arguments, reason = "one corner and its reading")]
fn close_corner(
    model: &mut Model,
    solid: &Shape,
    banded: Built,
    vertex: &Shape,
    frame: &crate::corner::CornerFrame,
    setback: &[f64],
    continuity: Continuity,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let ends = vec![CornerEnd::Round; frame.d.len()];
    close_corner_ends(
        model, solid, banded, vertex, frame, setback, &ends, continuity, tol,
    )
}

/// As [`close_corner`], each band ending as `ends` says.
#[allow(
    clippy::too_many_arguments,
    clippy::too_many_lines,
    reason = "one corner, assembled end to end"
)]
pub(crate) fn close_corner_ends(
    model: &mut Model,
    solid: &Shape,
    banded: Built,
    vertex: &Shape,
    frame: &crate::corner::CornerFrame,
    setback: &[f64],
    ends: &[CornerEnd],
    continuity: Continuity,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let n = frame.d.len();
    // The cell the patch replaces: in from every plane through the vertex,
    // and short of each band's section.
    let mut walls: Vec<(Vector, Point)> = frame.m.iter().map(|mk| (*mk, frame.corner)).collect();
    let sections: Vec<(Vector, Point)> = (0..n)
        .map(|k| {
            let dk = frame.d[k].vector();
            (-dk, frame.corner + dk * setback[k])
        })
        .collect();
    walls.extend(sections.iter().copied());
    let cell = crate::corner::convex_block(model, &walls, tol)?;
    let cut = ogeom_bool::cut(model, &banded.shape, &cell, tol)?;
    let planes: Vec<(Vector, Point)> = frame.m.iter().map(|mk| (*mk, frame.corner)).collect();
    guard_cell(
        model,
        solid,
        &banded.shape,
        &cut.shape,
        vertex,
        &planes,
        &walls,
        tol,
    )?;
    let history = banded.history.then(&cut.history);

    let on_plane = |p: Point, (normal, at): (Vector, Point)| -> bool {
        normal.dot(p - at).abs() <= tol.confusion() * 100.0
    };
    let section_of = |face: &Shape, model: &Model| -> OgeomResult<Option<usize>> {
        let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
            return Ok(None);
        };
        if !matches!(
            model.geometry().surface(data.surface),
            Some(ogeom_geom::SurfaceGeometry::Plane(_))
        ) {
            return Ok(None);
        }
        let (origin, outward) = ogeom_algo::face_normal(model, face, tol)?;
        let Some(k) = (0..n).find(|&k| {
            outward.cross(sections[k].0).magnitude() <= tol.angular() * 10.0
                && on_plane(origin, sections[k])
        }) else {
            return Ok(None);
        };
        for v in explore_unique(model, face, ShapeType::Vertex)? {
            if !in_cell(&walls, point_of(model, &v)?, tol) {
                return Ok(None);
            }
        }
        Ok(Some(k))
    };

    // Each band's end arc, on its section plane, and the band face it ends.
    let faces = explore_unique(model, &cut.shape, ShapeType::Face)?;
    let mut kept: Vec<Shape> = Vec::with_capacity(faces.len());
    for face in &faces {
        if section_of(face, model)?.is_none() {
            kept.push(face.clone());
        }
    }
    let mut arcs: Vec<(Shape, Shape)> = Vec::with_capacity(n);
    for (k, section) in sections.iter().enumerate() {
        if ends[k] == CornerEnd::Bevel {
            arcs.push(bevel_end(model, &kept, frame, *section, &walls, tol)?);
            continue;
        }
        arcs.push(round_end(model, &kept, k, *section, &walls, tol)?);
    }

    // On host plane `ring[k + 1]`, between band `k` and band `k + 1`: each
    // band's arc end there, and where the two section planes meet it.
    let host = |t: usize| -> (Vector, Point) { (frame.m[frame.ring[t]], frame.corner) };
    let end_on =
        |model: &Model, arc: &Shape, plane: (Vector, Point)| -> OgeomResult<(Shape, Point)> {
            let Some((a, b)) = ogeom_algo::edge_vertices(model, arc)? else {
                ogeom_bail!(Construction, "a band's end arc has no vertices");
            };
            for v in [a, b] {
                let p = point_of(model, &v)?;
                if on_plane(p, plane) {
                    return Ok((v, p));
                }
            }
            ogeom_bail!(Construction, "a band's end arc does not reach its face");
        };
    let fit = crate::marched::band_fit_target(tol);
    let mut boundary: Vec<FillBoundary> = Vec::with_capacity(2 * n);
    let mut fillers: Vec<Shape> = Vec::with_capacity(n);
    for k in 0..n {
        let next = (k + 1) % n;
        let plane = host(next);
        let (v0, p0) = end_on(model, &arcs[k].0, plane)?;
        let (v1, p1) = end_on(model, &arcs[next].0, plane)?;
        let Some(q) = crate::corner::planes_meet(&plane, &sections[k], &sections[next]) else {
            ogeom_bail!(
                Construction,
                "two section planes do not meet on the face between them"
            );
        };
        // The rails run along the edges: each conic leaves its band's rail
        // along that edge's direction, towards the vertex.
        let (t0, t1) = (-frame.d[k].vector(), -frame.d[next].vector());
        let (conic, apex) = trim_conic((p0, t0), (p1, t1), tol)?;
        let (trim, filler) =
            trim_and_filler(model, conic, apex, (&v0, p0), (&v1, p1), q, plane, tol)?;
        boundary.push(FillBoundary {
            edge: arcs[k].0.clone(),
            support: Some(arcs[k].1.clone()),
            continuity: match ends[k] {
                CornerEnd::Round => continuity,
                CornerEnd::Bevel => Continuity::C0,
            },
        });
        // A trim from a bevel's end meets its face in position only: the
        // patch runs along the bevel's end there, which leaves the face at
        // the bevel's angle, so it cannot take the face's tangent plane.
        // Between two rounds both ends lie tangent to the face.
        boundary.push(FillBoundary {
            edge: trim,
            support: Some(filler.clone()),
            continuity: if ends[k] == CornerEnd::Round && ends[next] == CornerEnd::Round {
                continuity
            } else {
                Continuity::C0
            },
        });
        fillers.push(filler);
    }

    // The vertex is convex, so what of the patch stands past one of its
    // planes lies outside the solid: the patch is held within them.
    let within: Vec<(Vector, Point)> = frame.m.iter().map(|m| (-*m, frame.corner)).collect();
    let filled = make_filling_n_within(model, &boundary, &[], fit, &within, tol)?;
    let patch = filled.built.shape.clone();
    let mut faces = kept;
    faces.extend(fillers);
    faces.push(patch.clone());
    let sewn = ogeom_algo::sew(model, &faces, tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the setback corner's faces did not close");
    }
    let rounded = model.add_solid(std::slice::from_ref(&sewn.shells[0]))?;
    // Sewing may rebuild a face whose edges it joins: an earlier corner's
    // patch, as well as this one, is followed to the face it became.
    let mut history: History = history.then(&sewn.history);
    for image in sewn.history.trace(&patch) {
        history.generate(vertex, image.clone());
    }
    Ok(Built::new(rounded, history))
}

/// The edge on `conic` from `v0` to `v1`, and the face's filler inside the
/// corner's cell: the planar face on `plane` bounded by the trim and the
/// two section lines from its ends to `q`, where they meet.
#[allow(clippy::too_many_arguments, reason = "one face's trim and its filler")]
pub(crate) fn trim_and_filler(
    model: &mut Model,
    conic: Curve,
    apex: Point,
    (v0, p0): (&Shape, Point),
    (v1, p1): (&Shape, Point),
    q: Point,
    plane: (Vector, Point),
    tol: Tolerances,
) -> OgeomResult<(Shape, Shape)> {
    let trim = ogeom_algo::make_edge_between(model, conic, (0.0, 1.0), v0, v1, tol)?.shape;
    let filler = filler_inside(
        model,
        vec![trim.clone()],
        apex,
        (v0, p0),
        (v1, p1),
        q,
        plane,
        tol,
    )?;
    Ok((trim, filler))
}

/// The face's filler inside the corner's cell: the planar face on `plane`
/// bounded by `trim`, the edges in order from `v0` to `v1` along a trim
/// whose middle control point is `apex`, and the two section lines from
/// its ends to `q`, where they meet.
#[allow(clippy::too_many_arguments, reason = "one face's filler")]
pub(crate) fn filler_inside(
    model: &mut Model,
    trim: Vec<Shape>,
    apex: Point,
    (v0, p0): (&Shape, Point),
    (v1, p1): (&Shape, Point),
    q: Point,
    plane: (Vector, Point),
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let vq = ogeom_algo::make_vertex(model, q).shape;
    let mut ring = trim;
    ring.push(crate::support::segment_between(
        model,
        (v1, p1),
        (&vq, q),
        tol,
    )?);
    ring.push(crate::support::segment_between(
        model,
        (&vq, q),
        (v0, p0),
        tol,
    )?);
    let outward = -plane.0;
    let reach = (p0.distance(q) + p1.distance(q) + apex.distance(q)) * 2.0 + 1.0;
    let surface: ogeom_geom::SurfaceGeometry = ogeom_geom::PlaneSurface::over(
        Plane::through(q, Direction::new(outward, tol)?),
        (-reach, reach),
        (-reach, reach),
    )?
    .into();
    let turn = (p0 - q).cross(p1 - q);
    let ring: Vec<Shape> = if turn.dot(outward) > 0.0 {
        ring
    } else {
        ring.iter().rev().map(Shape::reversed).collect()
    };
    Ok(ogeom_algo::make_face_with_pcurves(model, surface, &[ring], tol)?.shape)
}

/// The conic from `p0` to `p1` on a face, tangent at each end to the line
/// from it along its `t` (pointing into the conic), and its middle control
/// point, where those lines meet; weighted by the cosine of half the turn
/// between them, refused by name where the lines are parallel or do not run
/// towards each other.
pub(crate) fn trim_conic(
    (p0, t0): (Point, Vector),
    (p1, t1): (Point, Vector),
    tol: Tolerances,
) -> OgeomResult<(Curve, Point)> {
    let apex = {
        // p0 + a t0 = p1 + b t1, in the plane.
        let w = p1 - p0;
        let c = t0.cross(t1);
        let cc = c.dot(c);
        if cc <= tol.angular() {
            ogeom_bail!(Construction, "two rails on a face are parallel");
        }
        let a = w.cross(t1).dot(c) / cc;
        p0 + t0 * a
    };
    if (apex - p0).dot(t0) <= 0.0 || (apex - p1).dot(t1) <= 0.0 {
        ogeom_bail!(
            Construction,
            "the rails on a face do not run towards each other past their setbacks"
        );
    }
    let (u0, u1) = ((apex - p0).normalized(tol)?, (p1 - apex).normalized(tol)?);
    let weight = ((1.0 + u0.dot(u1)) / 2.0).max(0.0).sqrt();
    let knots = ogeom_math::KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2)?;
    let control = vec![
        ogeom_math::Weighted::new(p0, 1.0, tol)?,
        ogeom_math::Weighted::new(apex, weight, tol)?,
        ogeom_math::Weighted::new(p1, 1.0, tol)?,
    ];
    Ok((
        Curve::BSpline(ogeom_geom::BSplineCurve::rational(knots, control)?),
        apex,
    ))
}

/// Whether `p` lies inside the convex cell bounded by `walls`, each an
/// inward normal and a point on it, to the slack the corner's planes are
/// read to: this corner's, not another's on the same planes.
pub(crate) fn in_cell(walls: &[(Vector, Point)], p: Point, tol: Tolerances) -> bool {
    walls
        .iter()
        .all(|(normal, at)| normal.dot(p - *at) >= -tol.confusion() * 100.0)
}

/// Band `k`'s end arc on its section plane `section`, and the band face it
/// ends: the curved edge of a curved face of `kept` lying on the plane
/// inside the cell `walls` bound.
pub(crate) fn round_end(
    model: &Model,
    kept: &[Shape],
    k: usize,
    section: (Vector, Point),
    walls: &[(Vector, Point)],
    tol: Tolerances,
) -> OgeomResult<(Shape, Shape)> {
    let mut found: Option<(Shape, Shape)> = None;
    for face in kept {
        let is_plane = model
            .node(face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Plane(_)));
        if is_plane {
            continue;
        }
        for e in explore_unique(model, face, ShapeType::Edge)? {
            let Some((c, (lo, hi))) = ogeom_algo::placed_edge_curve(model, &e, tol)? else {
                continue;
            };
            if matches!(c, Curve::Line(_)) {
                continue;
            }
            let on = [lo, f64::midpoint(lo, hi), hi]
                .iter()
                .map(|t| c.point_at(*t, tol))
                .collect::<OgeomResult<Vec<Point>>>()?
                .into_iter()
                .all(|p| {
                    section.0.dot(p - section.1).abs() <= tol.confusion() * 100.0
                        && in_cell(walls, p, tol)
                });
            if on {
                found = Some((e, face.clone()));
            }
        }
    }
    let Some(found) = found else {
        ogeom_bail!(
            Construction,
            "band {k} does not end in an arc on its section plane"
        );
    };
    Ok(found)
}

/// A bevel's straight end on its section plane `section`, and the bevel's
/// face: the planar face of `kept`, on none of the corner's host planes,
/// with a straight edge lying on the section plane inside the cell `walls`
/// bound.
pub(crate) fn bevel_end(
    model: &Model,
    kept: &[Shape],
    frame: &crate::corner::CornerFrame,
    section: (Vector, Point),
    walls: &[(Vector, Point)],
    tol: Tolerances,
) -> OgeomResult<(Shape, Shape)> {
    let on_plane = |p: Point| {
        section.0.dot(p - section.1).abs() <= tol.confusion() * 100.0 && in_cell(walls, p, tol)
    };
    for face in kept {
        let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
            continue;
        };
        if !matches!(
            model.geometry().surface(data.surface),
            Some(ogeom_geom::SurfaceGeometry::Plane(_))
        ) {
            continue;
        }
        let (_, outward) = ogeom_algo::face_normal(model, face, tol)?;
        let host = frame
            .m
            .iter()
            .any(|mk| mk.cross(outward).magnitude() <= tol.angular() * 10.0);
        if host || outward.cross(section.0).magnitude() <= tol.angular() * 10.0 {
            continue;
        }
        for e in explore_unique(model, face, ShapeType::Edge)? {
            let Some((c, (lo, hi))) = ogeom_algo::placed_edge_curve(model, &e, tol)? else {
                continue;
            };
            if !matches!(c, Curve::Line(_)) {
                continue;
            }
            let on = [lo, hi]
                .iter()
                .map(|t| c.point_at(*t, tol))
                .collect::<OgeomResult<Vec<Point>>>()?
                .into_iter()
                .all(on_plane);
            if on {
                return Ok((e, face.clone()));
            }
        }
    }
    ogeom_bail!(
        Construction,
        "a bevel does not end in a straight edge on its section plane"
    );
}

/// Close the corner at `vertex` of `banded`, where the bands of `request`
/// (each edge with its radius) meet at different radii and no ball rolls
/// round: each band stopped one of its own widths past the least setback
/// where its section clears every other band, and the opening filled G1
/// as [`round_vertex_setback`] fills it.
///
/// `banded` holds the bands as they were built up to one another, the
/// later run on through the earlier; all of that lies inside the cell the
/// patch replaces, so what is left does not depend on the order they were
/// built in.
pub(crate) fn close_mixed_corner(
    model: &mut Model,
    solid: &Shape,
    banded: Built,
    vertex: &Shape,
    request: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Built> {
    let largest = request.iter().map(|(_, r)| *r).fold(0.0, f64::max);
    let frame = crate::corner::corner_frame(model, &banded.shape, vertex, largest, tol)?;
    let n = frame.d.len();
    // Each corner edge's radius: the request edge leaving the vertex along
    // it.
    let mut radius: Vec<f64> = Vec::with_capacity(n);
    for dk in &frame.d {
        let mut best: Option<(f64, f64)> = None;
        for (edge, r) in request {
            let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
                continue;
            };
            let far = if a.is_same(vertex) {
                b
            } else if b.is_same(vertex) {
                a
            } else {
                continue;
            };
            let along = point_of(model, &far)? - frame.corner;
            let length = along.magnitude();
            if length <= tol.confusion() {
                continue;
            }
            let cos = along.dot(dk.vector()) / length;
            if best.is_none_or(|(c, _)| cos > c) {
                best = Some((cos, *r));
            }
        }
        match best {
            Some((cos, r)) if cos >= 1.0 - tol.angular().sqrt() => radius.push(r),
            _ => ogeom_bail!(
                Construction,
                "a corner of different radii where an edge of the vertex is not \
                 rounded, or is curved, is still owed; see docs/CAPABILITIES.md, \
                 fillet.setback-corners"
            ),
        }
    }
    let setback = least_setbacks(&frame, &radius, tol)?;
    close_corner(
        model,
        solid,
        banded,
        vertex,
        &frame,
        &setback,
        Continuity::G1,
        tol,
    )
}

/// For each edge of the corner, the setback its band stops at: the least
/// distance along the edge at which the band's section clears every other
/// band of the corner, and past that the band's own width on its faces.
/// At the least distance a rail's end stands on its neighbour's rail line
/// and the face's trim would leave it square; one width further, the trim
/// leaves each rail along a leg as long as its band is wide.
///
/// Band `k` of radius `r_k` between planes at an interior angle `phi_k`
/// touches each plane `w_k = r_k (1 + cos phi_k) / sin phi_k` from its
/// edge, and everything it takes away lies within `w_k` of the edge's
/// line. Its section at distance `s` is the triangle of the edge's point
/// and the two touch points; the distance from another edge's line is
/// convex over it, so its corners decide. Each corner's distance to line
/// `j` is `|A s + B|`, increasing past `s* = -A·B / A·A`; the band clears
/// band `j` from the larger root of `|A s + B| = w_j`, or from `s*` where
/// it is never nearer.
fn least_setbacks(
    frame: &crate::corner::CornerFrame,
    radius: &[f64],
    tol: Tolerances,
) -> OgeomResult<Vec<f64>> {
    let n = frame.d.len();
    let planes = |k: usize| (frame.m[frame.ring[k]], frame.m[frame.ring[(k + 1) % n]]);
    let reach: Vec<f64> = (0..n)
        .map(|k| {
            let (ma, mb) = planes(k);
            let cos = -ma.dot(mb);
            let sin = (1.0 - cos * cos).max(0.0).sqrt();
            if sin <= tol.angular() {
                f64::INFINITY
            } else {
                radius[k] * (1.0 + cos) / sin
            }
        })
        .collect();
    let widths: Vec<(f64, f64)> = reach.iter().map(|&w| (w, w)).collect();
    least_setbacks_widths(frame, &widths, tol)
}

/// As [`least_setbacks`], each band taking `widths[k]` off the two faces of
/// its edge, in the order of the corner's ring: the first on the plane
/// `ring[k]`, the second on `ring[k + 1]`. Everything a band takes away
/// lies within the larger of its two widths of its edge's line.
pub(crate) fn least_setbacks_widths(
    frame: &crate::corner::CornerFrame,
    widths: &[(f64, f64)],
    tol: Tolerances,
) -> OgeomResult<Vec<f64>> {
    let n = frame.d.len();
    let planes = |k: usize| (frame.m[frame.ring[k]], frame.m[frame.ring[(k + 1) % n]]);
    let reach: Vec<f64> = widths.iter().map(|&(a, b)| a.max(b)).collect();
    let mut setback = vec![0.0_f64; n];
    let mut corners: Vec<[Vector; 3]> = Vec::with_capacity(n);
    for k in 0..n {
        if !reach[k].is_finite() {
            ogeom_bail!(
                Construction,
                "two faces at a corner meet flat; there is no edge to round"
            );
        }
        let dk = frame.d[k].vector();
        let (ma, mb) = planes(k);
        // Into each face, square to the edge.
        let leg = |own: Vector, other: Vector| -> OgeomResult<Vector> {
            let u = own.cross(dk).normalized(tol)?;
            Ok(if u.dot(other) < 0.0 { -u } else { u })
        };
        let offsets = [
            Vector::ZERO,
            leg(ma, mb)? * widths[k].0,
            leg(mb, ma)? * widths[k].1,
        ];
        let mut least = reach[k];
        for j in (0..n).filter(|&j| j != k) {
            let dj = frame.d[j].vector();
            let a = dk - dj * dk.dot(dj);
            let aa = a.dot(a);
            if aa <= tol.angular() {
                ogeom_bail!(Construction, "two edges at a corner run together");
            }
            for o in offsets {
                let b = o - dj * o.dot(dj);
                let ab = a.dot(b);
                let turn = -ab / aa;
                let disc = ab * ab - aa * (b.dot(b) - reach[j] * reach[j]);
                let from = if disc >= 0.0 {
                    ((-ab + disc.sqrt()) / aa).max(turn)
                } else {
                    turn
                };
                least = least.max(from);
            }
        }
        setback[k] = least + reach[k];
        corners.push(offsets);
    }
    // And each band, up to its own section, short of every other band's
    // section plane: band `j` takes away the sweep of its section triangle
    // from the vertex to `s_j`, which reaches along `d_k` as far as
    // `s_j max(0, d_j·d_k) + o·d_k` at its farthest corner `o`, a touch
    // point on its rail. Where that reaches band `k`'s section, band `k`
    // stops one of its own widths past it. The bounds tie the setbacks
    // together, each a share below one of the others, so they settle by
    // going round until none moves.
    let ties = |setback: &[f64], k: usize| -> f64 {
        let dk = frame.d[k].vector();
        let mut far = f64::NEG_INFINITY;
        for (j, offsets) in corners.iter().enumerate().filter(|&(j, _)| j != k) {
            let along = frame.d[j].vector().dot(dk).max(0.0);
            for o in offsets {
                far = far.max(setback[j].mul_add(along, o.dot(dk)));
            }
        }
        far
    };
    for _ in 0..=64 * n {
        let mut moved = false;
        for k in 0..n {
            let far = ties(&setback, k);
            // An edge of no width (one left sharp) moves only where it is
            // passed, or it would be moved by nothing for ever.
            if far > setback[k] || (far >= setback[k] && reach[k] > 0.0) {
                setback[k] = far + reach[k];
                moved = true;
            }
        }
        if !moved {
            return Ok(setback);
        }
    }
    ogeom_bail!(
        Construction,
        "the bands at this corner reach past each other's sections however far \
         they are set back"
    )
}

/// A vertex where a convex and a concave band of one request meet: three
/// straight edges, every one of them rounded, and the three planar faces
/// between them.
pub(crate) struct SignedCorner {
    /// The vertex itself, a node of the solid the request was given on.
    pub vertex: Shape,
    corner: Point,
    /// Each edge's unit direction away from the vertex, in ring order.
    d: Vec<Vector>,
    /// Each edge's radius.
    radius: Vec<f64>,
    /// Face `i`, between edge `i` and edge `i + 1`: its outward normal.
    normal: Vec<Vector>,
    /// Edge `i`'s legs: the unit directions square to it into face `i - 1`
    /// and face `i`, its two faces.
    legs: Vec<[Vector; 2]>,
}

/// A request edge with its radius and whether it is convex.
type Kinded = (Shape, f64, bool);

/// The vertices of `solid` where the request's bands of both kinds meet,
/// each read as a [`SignedCorner`]; refused by name where such a vertex is
/// not one this closes.
///
/// No one patch closes such a vertex unless every edge of it is rounded:
/// with an edge left sharp, the opening between the bands' ends crosses
/// it, and a patch tangent to both faces there would need two tangent
/// planes at one point.
pub(crate) fn signed_corners(
    model: &Model,
    solid: &Shape,
    request: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Vec<SignedCorner>> {
    // The vertices two or more of the request's edges meet at.
    let mut by_vertex: Vec<(Shape, Vec<(Shape, f64)>)> = Vec::new();
    for (edge, radius) in request {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
            continue;
        };
        for v in [a, b] {
            if let Some(slot) = by_vertex.iter_mut().find(|(held, _)| held.is_same(&v)) {
                if !slot.1.iter().any(|(e, _)| e.is_same(edge)) {
                    slot.1.push((edge.clone(), *radius));
                }
            } else {
                by_vertex.push((v, vec![(edge.clone(), *radius)]));
            }
        }
    }
    by_vertex.retain(|(_, at)| at.len() >= 2);
    // Each such edge's kind, where it reads.
    let mut kinds: Vec<(Shape, Option<bool>)> = Vec::new();
    let mut read = |model: &Model, edge: &Shape, radius: f64| -> OgeomResult<Option<bool>> {
        if let Some((_, kind)) = kinds.iter().find(|(e, _)| e.is_same(edge)) {
            return Ok(*kind);
        }
        let kind = match edge_is_convex(model, solid, edge, radius, tol) {
            Ok(kind) => kind,
            Err(ogeom_core::OgeomError::Construction(_) | ogeom_core::OgeomError::NotDone(_)) => {
                None
            }
            Err(e) => return Err(e),
        };
        kinds.push((edge.clone(), kind));
        Ok(kind)
    };
    let mut signed: Vec<(Shape, Vec<Kinded>)> = Vec::new();
    for (vertex, at) in by_vertex {
        let mut known = Vec::with_capacity(at.len());
        for (edge, radius) in at {
            if let Some(convex) = read(model, &edge, radius)? {
                known.push((edge, radius, convex));
            }
        }
        signed.push((vertex, known));
    }
    let by_vertex = signed;
    let mut corners = Vec::new();
    for (vertex, at) in by_vertex {
        if !(at.iter().any(|e| e.2) && at.iter().any(|e| !e.2)) {
            continue;
        }
        corners.push(signed_corner(model, solid, &vertex, &at, tol)?);
    }
    Ok(corners)
}

/// Whether `edge` of `solid` is convex: from its planar seat, else read on
/// its two faces at its middle as a marched seat reads it; `None` where it
/// does not part two faces.
fn edge_is_convex(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<bool>> {
    if let Some(seat) = crate::consumed::seat_of(model, solid, edge, tol)? {
        return Ok(Some(seat.convex));
    }
    let mut faces: Vec<Shape> = Vec::with_capacity(2);
    for f in explore_unique(model, solid, ShapeType::Face)? {
        if explore_unique(model, &f, ShapeType::Edge)?
            .iter()
            .any(|e| e.is_same(edge))
        {
            faces.push(f);
        }
    }
    let [first, second] = faces.as_slice() else {
        return Ok(None);
    };
    let (Some(a), Some(b)) = (
        crate::rolled::host(model, first, tol)?,
        crate::rolled::host(model, second, tol)?,
    ) else {
        return Ok(None);
    };
    let Some((curve, range)) = ogeom_algo::placed_edge_curve(model, edge, tol)? else {
        return Ok(None);
    };
    crate::marched::crease_convexity(
        model,
        first,
        [(&a.0, a.1), (&b.0, b.1)],
        &curve,
        range,
        radius,
        tol,
    )
    .map(Some)
}

/// The refusal for a signed corner this does not close.
fn signed_refusal(why: &str) -> ogeom_core::OgeomError {
    ogeom_core::ogeom_err!(
        Construction,
        "a convex and a concave round meeting at one vertex close as one corner \
         only where {why}; no one patch is tangent to both faces where the opening \
         between the bands crosses an edge left sharp, and the two rounds alone \
         would take the corner in the order asked; see docs/CAPABILITIES.md, \
         fillet.setback-corners"
    )
}

fn signed_corner(
    model: &Model,
    solid: &Shape,
    vertex: &Shape,
    request: &[Kinded],
    tol: Tolerances,
) -> OgeomResult<SignedCorner> {
    let corner = point_of(model, vertex)?;
    // The faces at the vertex, coplanar neighbours (a fused solid's seam)
    // taken as one: each group's planar faces and its outward normal.
    let mut groups: Vec<(Vec<Shape>, crate::rolled::Flat)> = Vec::new();
    for f in explore_unique(model, solid, ShapeType::Face)? {
        if !explore_unique(model, &f, ShapeType::Vertex)?
            .iter()
            .any(|v| v.is_same(vertex))
        {
            continue;
        }
        let Some(flat) = crate::rolled::flat(model, &f, tol)? else {
            return Err(signed_refusal("every face of the vertex is planar"));
        };
        let same = groups.iter_mut().find(|(_, g)| {
            g.normal.dot(flat.normal) >= 1.0 - tol.angular()
                && g.normal.dot(flat.origin - g.origin).abs() <= tol.confusion() * 100.0
        });
        match same {
            Some((members, g)) => {
                members.push(f.clone());
                g.sides.extend(flat.sides.iter().copied());
            }
            None => groups.push((vec![f.clone()], flat)),
        }
    }
    // The edges at the vertex between two groups; a seam inside one is none.
    let mut at_vertex: Vec<Shape> = Vec::new();
    let mut edge_groups: Vec<[usize; 2]> = Vec::new();
    for e in explore_unique(model, solid, ShapeType::Edge)? {
        if !ogeom_algo::edge_vertices(model, &e)?
            .is_some_and(|(a, b)| a.is_same(vertex) || b.is_same(vertex))
        {
            continue;
        }
        let mut own: Vec<usize> = Vec::new();
        for (gi, (members, _)) in groups.iter().enumerate() {
            for f in members {
                if explore_unique(model, f, ShapeType::Edge)?
                    .iter()
                    .any(|x| x.is_same(&e))
                    && !own.contains(&gi)
                {
                    own.push(gi);
                }
            }
        }
        match own.as_slice() {
            [_] => {}
            [a, b] => {
                at_vertex.push(e);
                edge_groups.push([*a, *b]);
            }
            _ => {
                return Err(signed_refusal(
                    "each edge of the vertex parts two of its faces",
                ));
            }
        }
    }
    if at_vertex.len() != 3 || groups.len() != 3 {
        return Err(signed_refusal("three edges meet at the vertex"));
    }
    let mut radius_of = Vec::with_capacity(3);
    let mut dir = Vec::with_capacity(3);
    let mut reach = Vec::with_capacity(3);
    for e in &at_vertex {
        let Some((_, r, _)) = request.iter().find(|(asked, _, _)| asked.is_same(e)) else {
            return Err(signed_refusal(
                "the request rounds every edge of the vertex",
            ));
        };
        let Some([p, q]) = crate::consumed::segment(model, e, tol)? else {
            return Err(signed_refusal("every edge of the vertex is straight"));
        };
        let far = if p.distance(corner) < q.distance(corner) {
            q
        } else {
            p
        };
        let along = far - corner;
        reach.push(along.magnitude());
        dir.push(along.normalized(tol)?);
        radius_of.push(*r);
    }
    // Each face's two edges at the vertex.
    let flats: Vec<crate::rolled::Flat> = groups.iter().map(|(_, g)| g.clone()).collect();
    let mut pairs: Vec<[usize; 2]> = Vec::with_capacity(3);
    for gi in 0..3 {
        let own: Vec<usize> = (0..3).filter(|&k| edge_groups[k].contains(&gi)).collect();
        let [a, b] = own.as_slice() else {
            return Err(signed_refusal(
                "each face of the vertex holds two of its edges",
            ));
        };
        pairs.push([*a, *b]);
    }
    // The ring: edge 0, the face after it, the edge after that, and so on.
    let mut ring_edges = vec![0_usize];
    let mut ring_faces: Vec<usize> = Vec::with_capacity(3);
    let mut used = [false; 3];
    while ring_faces.len() < 3 {
        let last = *ring_edges.last().unwrap_or(&0);
        let Some(fi) = (0..3).find(|&i| !used[i] && pairs[i].contains(&last)) else {
            return Err(signed_refusal("the faces of the vertex close round it"));
        };
        used[fi] = true;
        ring_faces.push(fi);
        let next = if pairs[fi][0] == last {
            pairs[fi][1]
        } else {
            pairs[fi][0]
        };
        if ring_faces.len() < 3 {
            ring_edges.push(next);
        }
    }
    let d: Vec<Vector> = ring_edges.iter().map(|&k| dir[k]).collect();
    let radius: Vec<f64> = ring_edges.iter().map(|&k| radius_of[k]).collect();
    let normal: Vec<Vector> = ring_faces.iter().map(|&f| flats[f].normal).collect();
    // Into a face square to an edge: the side a point just off the edge's
    // middle lies inside the face on.
    let into = |fi: usize, k: usize| -> OgeomResult<Vector> {
        let flat = &flats[fi];
        let u = flat.normal.cross(dir[k]).normalized(tol)?;
        let base = corner + dir[k] * (reach[k] * 0.5);
        let step = reach[k] * 1e-3;
        let frame = crate::rolled::basis(flat.normal);
        match (
            crate::rolled::inside(&flat.sides, base + u * step, frame),
            crate::rolled::inside(&flat.sides, base - u * step, frame),
        ) {
            (true, false) => Ok(u),
            (false, true) => Ok(-u),
            _ => Err(signed_refusal(
                "each face lies on one side of each of its edges",
            )),
        }
    };
    let mut legs = Vec::with_capacity(3);
    for i in 0..3 {
        let k = ring_edges[i];
        let before = ring_faces[(i + 2) % 3];
        let after = ring_faces[i];
        legs.push([into(before, k)?, into(after, k)?]);
    }
    Ok(SignedCorner {
        vertex: vertex.clone(),
        corner,
        d,
        radius,
        normal,
        legs,
    })
}

/// Each band's width on its faces: `r tan(a / 2)`, `a` the turn between
/// its faces' normals, the same for a convex band and a concave one.
fn signed_widths(c: &SignedCorner, tol: Tolerances) -> OgeomResult<Vec<f64>> {
    (0..3)
        .map(|i| {
            let (a, b) = (c.normal[(i + 2) % 3], c.normal[i]);
            let turn = a.dot(b).clamp(-1.0, 1.0).acos();
            if turn <= tol.angular() || turn >= core::f64::consts::PI - tol.angular() {
                return Err(signed_refusal(
                    "the faces of each edge turn by less than a half turn",
                ));
            }
            Ok(c.radius[i] * (turn / 2.0).tan())
        })
        .collect()
}

/// For each edge, the setback its band stops at: as [`least_setbacks`]
/// reads it from each band's section triangle (the edge's point and its
/// touch on each face), which holds what the band takes away or adds,
/// whichever its kind.
fn signed_setbacks(c: &SignedCorner, width: &[f64], tol: Tolerances) -> OgeomResult<Vec<f64>> {
    let mut setback = vec![0.0_f64; 3];
    for k in 0..3 {
        let dk = c.d[k];
        let offsets = [
            Vector::ZERO,
            c.legs[k][0] * width[k],
            c.legs[k][1] * width[k],
        ];
        let mut least = width[k];
        for j in (0..3).filter(|&j| j != k) {
            let dj = c.d[j];
            let a = dk - dj * dk.dot(dj);
            let aa = a.dot(a);
            if aa <= tol.angular() {
                return Err(signed_refusal("no two edges of the vertex run together"));
            }
            for o in offsets {
                let b = o - dj * o.dot(dj);
                let ab = a.dot(b);
                let turn = -ab / aa;
                let disc = ab * ab - aa * (b.dot(b) - width[j] * width[j]);
                let from = if disc >= 0.0 {
                    ((-ab + disc.sqrt()) / aa).max(turn)
                } else {
                    turn
                };
                least = least.max(from);
            }
        }
        setback[k] = least + width[k];
    }
    Ok(setback)
}

/// Close the signed corner `c` of `banded`: each band stopped at its
/// setback and the opening between them filled by one patch G1 to every
/// band end and face.
///
/// The cell cut away is the box square to the three edges, from each
/// band's section back past the vertex by twice the widest band, so it
/// holds every meeting of the bands; whatever order they were built in,
/// all of that goes. Its walls behind the vertex stand in the material,
/// where the solid before and after agree, so the faces the cut leaves on
/// every wall are dropped and nothing is put back there. On each face the
/// trim runs on a conic from one band's rail end to the next one's,
/// tangent to both; the face's part kept inside the box is the side of the
/// conic away from the vertex, out to the box's walls.
#[allow(clippy::too_many_lines, reason = "one corner, assembled end to end")]
pub(crate) fn close_signed_corner(
    model: &mut Model,
    solid: &Shape,
    banded: Built,
    c: &SignedCorner,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let width = signed_widths(c, tol)?;
    let setback = signed_setbacks(c, &width, tol)?;
    let back = 2.0 * width.iter().copied().fold(0.0, f64::max);
    let sections: Vec<(Vector, Point)> = (0..3)
        .map(|k| (-c.d[k], c.corner + c.d[k] * setback[k]))
        .collect();
    let mut walls = sections.clone();
    walls.extend((0..3).map(|k| (c.d[k], c.corner - c.d[k] * back)));
    let cell = crate::corner::convex_block(model, &walls, tol)?;
    let cut = ogeom_bool::cut(model, &banded.shape, &cell, tol)?;
    let planes: Vec<(Vector, Point)> = c.normal.iter().map(|n| (*n, c.corner)).collect();
    guard_cell(
        model,
        solid,
        &banded.shape,
        &cut.shape,
        &c.vertex,
        &planes,
        &walls,
        tol,
    )?;
    let history = banded.history.then(&cut.history);
    let slack = tol.confusion() * 100.0;
    let on_plane =
        |p: Point, (normal, at): (Vector, Point)| -> bool { normal.dot(p - at).abs() <= slack };
    let is_plane = |face: &Shape, model: &Model| {
        model
            .node(face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Plane(_)))
    };
    // Inside the cell, to the slack: this corner's, not another's on the
    // same planes.
    let in_cell = |p: Point| walls.iter().all(|(n, at)| n.dot(p - *at) >= -slack);
    let on_wall = |face: &Shape, model: &Model| -> OgeomResult<bool> {
        if !is_plane(face, model) {
            return Ok(false);
        }
        let (origin, outward) = ogeom_algo::face_normal(model, face, tol)?;
        if !walls
            .iter()
            .any(|w| outward.cross(w.0).magnitude() <= tol.angular() * 10.0 && on_plane(origin, *w))
        {
            return Ok(false);
        }
        for v in explore_unique(model, face, ShapeType::Vertex)? {
            if !in_cell(point_of(model, &v)?) {
                return Ok(false);
            }
        }
        Ok(true)
    };
    let faces = explore_unique(model, &cut.shape, ShapeType::Face)?;
    let mut kept: Vec<Shape> = Vec::with_capacity(faces.len());
    for face in &faces {
        if !on_wall(face, model)? {
            kept.push(face.clone());
        }
    }
    // Every corner of the kept faces.
    let mut held: Vec<Point> = Vec::new();
    for face in &kept {
        for v in explore_unique(model, face, ShapeType::Vertex)? {
            let p = point_of(model, &v)?;
            if !held.iter().any(|q| q.distance(p) <= slack) {
                held.push(p);
            }
        }
    }
    // Each band's end arc on its section plane, and the band face it ends.
    let mut arcs: Vec<(Shape, Shape)> = Vec::with_capacity(3);
    for (k, section) in sections.iter().enumerate() {
        let mut found: Option<(Shape, Shape)> = None;
        for face in &kept {
            if is_plane(face, model) {
                continue;
            }
            for e in explore_unique(model, face, ShapeType::Edge)? {
                let Some((curve, (lo, hi))) = ogeom_algo::placed_edge_curve(model, &e, tol)? else {
                    continue;
                };
                if matches!(curve, Curve::Line(_)) {
                    continue;
                }
                let on = [lo, f64::midpoint(lo, hi), hi]
                    .iter()
                    .map(|t| curve.point_at(*t, tol))
                    .collect::<OgeomResult<Vec<Point>>>()?
                    .into_iter()
                    .all(|p| on_plane(p, *section) && in_cell(p));
                if on {
                    found = Some((e, face.clone()));
                }
            }
        }
        let Some(found) = found else {
            ogeom_bail!(
                Construction,
                "band {k} of a signed corner does not end in an arc on its section plane"
            );
        };
        arcs.push(found);
    }
    let end_on =
        |model: &Model, arc: &Shape, plane: (Vector, Point)| -> OgeomResult<(Shape, Point)> {
            let Some((a, b)) = ogeom_algo::edge_vertices(model, arc)? else {
                ogeom_bail!(Construction, "a band's end arc has no vertices");
            };
            for v in [a, b] {
                let p = point_of(model, &v)?;
                if on_plane(p, plane) {
                    return Ok((v, p));
                }
            }
            ogeom_bail!(Construction, "a band's end arc does not reach its face");
        };
    let fit = crate::marched::band_fit_target(tol);
    let mut boundary: Vec<FillBoundary> = Vec::with_capacity(6);
    let mut fillers: Vec<Shape> = Vec::with_capacity(3);
    for k in 0..3 {
        let next = (k + 1) % 3;
        let plane = (c.normal[k], c.corner);
        let (v0, p0) = end_on(model, &arcs[k].0, plane)?;
        let (v1, p1) = end_on(model, &arcs[next].0, plane)?;
        let (t0, t1) = (-c.d[k], -c.d[next]);
        let apex = {
            let w = p1 - p0;
            let cr = t0.cross(t1);
            let cc = cr.dot(cr);
            if cc <= tol.angular() {
                ogeom_bail!(Construction, "two rails on a face are parallel");
            }
            p0 + t0 * (w.cross(t1).dot(cr) / cc)
        };
        if (apex - p0).dot(t0) <= 0.0 || (apex - p1).dot(t1) <= 0.0 {
            ogeom_bail!(
                Construction,
                "the rails on a face do not run towards each other past their setbacks"
            );
        }
        let (u0, u1) = ((apex - p0).normalized(tol)?, (p1 - apex).normalized(tol)?);
        let weight = ((1.0 + u0.dot(u1)) / 2.0).max(0.0).sqrt();
        let knots = ogeom_math::KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2)?;
        let control = vec![
            ogeom_math::Weighted::new(p0, 1.0, tol)?,
            ogeom_math::Weighted::new(apex, weight, tol)?,
            ogeom_math::Weighted::new(p1, 1.0, tol)?,
        ];
        let conic = ogeom_geom::BSplineCurve::rational(knots, control)?;
        let samples: Vec<Point> = (0..=16)
            .map(|i| conic.point_at(f64::from(i) / 16.0, tol))
            .collect::<OgeomResult<_>>()?;
        let trim =
            ogeom_algo::make_edge_between(model, Curve::BSpline(conic), (0.0, 1.0), &v0, &v1, tol)?
                .shape;
        // The face's plane inside the box: a convex polygon, ordered round
        // its middle.
        let outward = c.normal[k];
        let mut corners: Vec<Point> = Vec::new();
        for i in 0..walls.len() {
            for j in i + 1..walls.len() {
                let Some(p) = crate::corner::planes_meet(&plane, &walls[i], &walls[j]) else {
                    continue;
                };
                if walls.iter().any(|(n, at)| n.dot(p - *at) < -slack)
                    || corners.iter().any(|q| q.distance(p) <= slack)
                {
                    continue;
                }
                corners.push(p);
            }
        }
        if corners.len() < 3 {
            ogeom_bail!(Construction, "a face of a signed corner misses its cell");
        }
        let (e1, e2) = crate::rolled::basis(outward);
        let middle = corners
            .iter()
            .fold(Vector::ZERO, |acc, p| acc + (*p - Point::ORIGIN))
            * (1.0 / 3.0_f64.max(f64::from(u32::try_from(corners.len()).unwrap_or(3))));
        let middle = Point::ORIGIN + middle;
        let angle = |p: &Point| (*p - middle).dot(e2).atan2((*p - middle).dot(e1));
        corners.sort_by(|a, b| angle(a).total_cmp(&angle(b)));
        // Where each rail end sits on the polygon's boundary.
        let place = |p: Point| -> Option<(usize, f64)> {
            let m = corners.len();
            (0..m).find_map(|i| {
                let (a, b) = (corners[i], corners[(i + 1) % m]);
                let ab = b - a;
                let t = (p - a).dot(ab) / ab.dot(ab);
                ((-1e-9..=1.0 + 1e-9).contains(&t) && (a + ab * t).distance(p) <= slack)
                    .then_some((i, t))
            })
        };
        let (Some(at0), Some(at1)) = (place(p0), place(p1)) else {
            ogeom_bail!(Construction, "a band's rail end misses its section wall");
        };
        // The two ways round the polygon from p1 back to p0.
        let walk = |forward: bool| -> Vec<Point> {
            let m = corners.len();
            let mut path = Vec::new();
            if forward {
                let mut i = at1.0;
                loop {
                    i = (i + 1) % m;
                    path.push(corners[i]);
                    if i == at0.0 {
                        break;
                    }
                }
            } else {
                let mut i = at1.0;
                loop {
                    path.push(corners[i]);
                    if i == (at0.0 + 1) % m {
                        break;
                    }
                    i = (i + m - 1) % m;
                }
            }
            path
        };
        let flat2 = |p: &Point| ((*p - c.corner).dot(e1), (*p - c.corner).dot(e2));
        let holds_vertex = |path: &[Point]| -> bool {
            let mut ring: Vec<(f64, f64)> = samples.iter().map(flat2).collect();
            ring.extend(path.iter().map(flat2));
            let mut odd = false;
            for i in 0..ring.len() {
                let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
                if (a.1 > 0.0) != (b.1 > 0.0) && 0.0 < a.0 + (b.0 - a.0) * (0.0 - a.1) / (b.1 - a.1)
                {
                    odd = !odd;
                }
            }
            odd
        };
        let same_side = at0.0 == at1.0;
        let path = if same_side {
            Vec::new()
        } else {
            let (a, b) = (walk(true), walk(false));
            match (holds_vertex(&a), holds_vertex(&b)) {
                (false, true) => a,
                (true, false) => b,
                _ => ogeom_bail!(
                    Construction,
                    "a face's trim does not part the vertex from the cell's walls"
                ),
            }
        };
        let path: Vec<Point> = path
            .into_iter()
            .filter(|p| p.distance(p0) > slack && p.distance(p1) > slack)
            .collect();
        // The straight run back from p1 to p0, broken wherever a kept
        // face's corner stands on it (a seam of coplanar faces running
        // into the cell), so each piece meets one kept edge.
        let mut stops: Vec<Point> = Vec::new();
        let mut from = p1;
        for p in path.iter().copied().chain(core::iter::once(p0)) {
            let run = p - from;
            let length = run.magnitude();
            let mut between: Vec<(f64, Point)> = held
                .iter()
                .filter_map(|q| {
                    let t = (*q - from).dot(run) / (length * length);
                    (t > 0.0
                        && t < 1.0
                        && (from + run * t).distance(*q) <= slack
                        && q.distance(from) > slack
                        && q.distance(p) > slack)
                        .then_some((t, *q))
                })
                .collect();
            between.sort_by(|a, b| a.0.total_cmp(&b.0));
            stops.extend(between.into_iter().map(|(_, q)| q));
            stops.push(p);
            from = p;
        }
        stops.pop();
        let mut ring = vec![trim.clone()];
        let mut prev = (v1.clone(), p1);
        for p in &stops {
            let v = ogeom_algo::make_vertex(model, *p).shape;
            ring.push(crate::support::segment_between(
                model,
                (&prev.0, prev.1),
                (&v, *p),
                tol,
            )?);
            prev = (v, *p);
        }
        ring.push(crate::support::segment_between(
            model,
            (&prev.0, prev.1),
            (&v0, p0),
            tol,
        )?);
        let mut loop2: Vec<(f64, f64)> = samples.iter().map(flat2).collect();
        loop2.extend(path.iter().map(flat2));
        let area: f64 = (0..loop2.len())
            .map(|i| {
                let (a, b) = (loop2[i], loop2[(i + 1) % loop2.len()]);
                a.0 * b.1 - b.0 * a.1
            })
            .sum();
        let ring: Vec<Shape> = if area > 0.0 {
            ring
        } else {
            ring.iter().rev().map(Shape::reversed).collect()
        };
        let reach = corners
            .iter()
            .map(|p| p.distance(c.corner))
            .fold(0.0, f64::max)
            * 2.0
            + 1.0;
        let surface: ogeom_geom::SurfaceGeometry = ogeom_geom::PlaneSurface::over(
            Plane::through(c.corner, Direction::new(outward, tol)?),
            (-reach, reach),
            (-reach, reach),
        )?
        .into();
        let filler = ogeom_algo::make_face_with_pcurves(model, surface, &[ring], tol)?.shape;
        boundary.push(FillBoundary {
            edge: arcs[k].0.clone(),
            support: Some(arcs[k].1.clone()),
            continuity: Continuity::G1,
        });
        boundary.push(FillBoundary {
            edge: trim,
            support: Some(filler.clone()),
            continuity: Continuity::G1,
        });
        fillers.push(filler);
    }
    let filled = make_filling_n(model, &boundary, &[], fit, tol)?;
    let patch = filled.built.shape.clone();
    let mut faces = kept;
    faces.extend(fillers);
    faces.push(patch.clone());
    let sewn = ogeom_algo::sew(model, &faces, tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the signed corner's faces did not close");
    }
    let rounded = model.add_solid(std::slice::from_ref(&sewn.shells[0]))?;
    let mut history: History = history;
    history.generate(&c.vertex, patch);
    Ok(Built::new(rounded, history))
}
