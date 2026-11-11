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

use crate::fill_n::{FillBoundary, make_filling_n_region, make_filling_n_within};

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
    /// The edge the request leaves sharp, by its place in the ring: a band
    /// of no width, its radius nought.
    sharp: Option<usize>,
    /// Each edge's length, from the vertex to its far end.
    reach: Vec<f64>,
}

/// A request edge with its radius and whether it is convex.
type Kinded = (Shape, f64, bool);

/// The vertices of `solid` where the request's bands of both kinds meet,
/// each read as a [`SignedCorner`]; refused by name where such a vertex is
/// not one this closes.
///
/// An edge of such a vertex the request leaves sharp is read as a band of
/// no width: the opening between the bands' ends crosses it, and is closed
/// through a point on it where a patch turns its normal from one of its
/// faces to the other's ([`crate::sharp::close_through_q`]). Such an edge
/// is to be convex.
pub(crate) fn signed_corners(
    model: &Model,
    solid: &Shape,
    request: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Vec<SignedCorner>> {
    let mut corners = Vec::new();
    for (vertex, at) in signed_vertices(model, solid, request, tol)? {
        if runs_through(model, solid, &vertex, &at, tol)?
            || crate::cusp::cusp_corner(model, solid, &vertex, &at, tol)?.is_some()
        {
            continue;
        }
        corners.push(signed_corner(model, solid, &vertex, &at, tol)?);
    }
    Ok(corners)
}

/// The cusp corners of the request: vertices where its convex and concave
/// bands meet that [`crate::cusp::cusp_corner`] reads as one.
pub(crate) fn cusp_corners(
    model: &Model,
    solid: &Shape,
    request: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Vec<crate::cusp::CuspCorner>> {
    // Read in node order, so the corner and its sides are the same in any
    // order asked.
    let mut sorted = request.to_vec();
    sorted.sort_by_key(|(edge, _)| edge.node());
    let mut corners = Vec::new();
    for (vertex, at) in signed_vertices(model, solid, &sorted, tol)? {
        if let Some(c) = crate::cusp::cusp_corner(model, solid, &vertex, &at, tol)? {
            corners.push(c);
        }
    }
    Ok(corners)
}

/// The vertices two or more of the request's edges meet at, bands of both
/// kinds among them, each with its edges, radii and kinds.
fn signed_vertices(
    model: &Model,
    solid: &Shape,
    request: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Vec<(Shape, Vec<Kinded>)>> {
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
    Ok(signed
        .into_iter()
        .filter(|(_, at)| at.iter().any(|e| e.2) && at.iter().any(|e| !e.2))
        .collect())
}

/// Whether the request's two edges at `vertex`, one convex and one
/// concave, run on through it into each other, tangent there and leaving
/// it opposite ways, where a face of the vertex is curved: a tangent chain
/// turning from convex to concave where a wall runs tangent into a curved
/// one (a plate's side into a boss). The two bands then lie either side of
/// the vertex and meet only across their end sections there, so there is
/// no corner to close: each is built as the request's own bands are, and
/// refused as they are where its face narrows to a cusp there and would
/// leave the band's fill over nothing.
fn runs_through(
    model: &Model,
    solid: &Shape,
    vertex: &Shape,
    at: &[Kinded],
    tol: Tolerances,
) -> OgeomResult<bool> {
    let [(first, _, _), (second, _, _)] = at else {
        return Ok(false);
    };
    let corner = point_of(model, vertex)?;
    // Each edge's unit tangent at the vertex, pointing away from it.
    let leaving = |edge: &Shape| -> OgeomResult<Option<Vector>> {
        let Some((curve, (lo, hi))) = ogeom_algo::placed_edge_curve(model, edge, tol)? else {
            return Ok(None);
        };
        let (at_lo, at_hi) = (curve.point_at(lo, tol)?, curve.point_at(hi, tol)?);
        let (t, sense) = if at_lo.distance(corner) <= at_hi.distance(corner) {
            (lo, 1.0)
        } else {
            (hi, -1.0)
        };
        Ok(Some(curve.d1_at(t, tol)?.normalized(tol)? * sense))
    };
    let (Some(a), Some(b)) = (leaving(first)?, leaving(second)?) else {
        return Ok(false);
    };
    if a.dot(b) > -1.0 + tol.angular() * 10.0 {
        return Ok(false);
    }
    for f in explore_unique(model, solid, ShapeType::Face)? {
        if explore_unique(model, &f, ShapeType::Vertex)?
            .iter()
            .any(|v| v.is_same(vertex))
            && crate::rolled::flat(model, &f, tol)?.is_none()
        {
            return Ok(true);
        }
    }
    Ok(false)
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
         only where {why}; the two rounds alone would take the corner in the \
         order asked; see docs/CAPABILITIES.md, fillet.setback-corners"
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
    let mut sharp_edge: Option<usize> = None;
    for (k, e) in at_vertex.iter().enumerate() {
        let r = match request.iter().find(|(asked, _, _)| asked.is_same(e)) {
            Some((_, r, _)) => *r,
            None if sharp_edge.is_none() => {
                sharp_edge = Some(k);
                0.0
            }
            None => {
                return Err(signed_refusal(
                    "the request rounds all but at most one edge of the vertex",
                ));
            }
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
        radius_of.push(r);
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
    let length: Vec<f64> = ring_edges.iter().map(|&k| reach[k]).collect();
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
    let sharp = sharp_edge.and_then(|k| ring_edges.iter().position(|&e| e == k));
    // A sharp edge's patch is a cone at its point, which turns outward
    // across a convex edge: the face after it leans in from the face
    // before's outward side.
    if let Some(i) = sharp
        && legs[i][1].dot(normal[(i + 2) % 3]) >= 0.0
    {
        return Err(signed_refusal("an edge the request leaves sharp is convex"));
    }
    Ok(SignedCorner {
        vertex: vertex.clone(),
        corner,
        d,
        radius,
        normal,
        legs,
        sharp,
        reach: length,
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
        for j in (0..3).filter(|&j| j != k && c.sharp != Some(j)) {
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
    // The sharp edge is cut past where the bands reach along it by the
    // widest band, so the trims from its point have room to turn into the
    // bands' rails.
    if let Some(k) = c.sharp {
        setback[k] += width.iter().copied().fold(0.0, f64::max);
    }
    // Each edge is cut short of its far end, where its band's section
    // would fall on the face there and leave no end arc.
    for (s, length) in setback.iter().zip(&c.reach) {
        if *s >= length - tol.confusion() * 100.0 {
            return Err(signed_refusal(
                "each edge runs on past its band's setback, its own width beyond \
                 the least at which its section clears the others",
            ));
        }
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
///
/// With an edge left sharp, that edge is cut at its setback, at `Q`, and
/// the trims beside it run from `Q` into the bands' rails, split at their
/// middles; on a face reflex at the vertex the trim is a cubic swinging
/// round the vertex inside the face. The opening is then closed as
/// [`crate::round_vertex_setback_sharp`] closes it: a patch at `Q`, its
/// fourth side shrunk there, and one over the rest.
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
    let mut arcs: Vec<Option<(Shape, Shape)>> = Vec::with_capacity(3);
    for (k, section) in sections.iter().enumerate() {
        if c.sharp == Some(k) {
            arcs.push(None);
            continue;
        }
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
        arcs.push(Some(found));
    }
    // Where the cut leaves the sharp edge, at its setback: `Q`.
    let at_q: Option<(Shape, Point)> = match c.sharp {
        None => None,
        Some(k) => {
            let q = sections[k].1;
            let Some(v) = explore_unique(model, &cut.shape, ShapeType::Vertex)?
                .into_iter()
                .find(|v| point_of(model, v).is_ok_and(|p| p.distance(q) <= slack))
            else {
                ogeom_bail!(
                    Construction,
                    "the cut does not leave the sharp edge of a signed corner at its setback"
                );
            };
            Some((v, q))
        }
    };
    let end_on = |model: &Model, k: usize, plane: (Vector, Point)| -> OgeomResult<(Shape, Point)> {
        let Some((arc, _)) = &arcs[k] else {
            return at_q
                .clone()
                .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "a band has no end arc"));
        };
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
    let mut trims: Vec<crate::sharp::Trim> = Vec::with_capacity(3);
    #[allow(
        clippy::needless_range_loop,
        reason = "face k lies between edge k and edge k + 1"
    )]
    for k in 0..3 {
        let next = (k + 1) % 3;
        let plane = (c.normal[k], c.corner);
        let (v0, p0) = end_on(model, k, plane)?;
        let (v1, p1) = end_on(model, next, plane)?;
        // A rail's trim leaves it along its edge, towards the vertex; one
        // from the sharp edge's point heads for the band's rail half way
        // back to the vertex.
        let into = |from: (usize, Point), to: (usize, Point)| -> OgeomResult<Vector> {
            if c.sharp == Some(from.0) {
                let aim = to.1 - c.d[to.0] * (setback[to.0] / 2.0);
                Ok((aim - from.1).normalized(tol)?)
            } else {
                Ok(-c.d[from.0])
            }
        };
        // On a face reflex at the vertex, beside the sharp edge, the trim
        // from its point swings round the vertex inside the face: a cubic
        // leaving the point into the face half way between the sharp edge,
        // back towards the vertex, and square across it, and coming onto
        // the band's rail from a band's width past the sharp edge's line.
        let reflex = c.sharp.and_then(|ks| {
            let (band, leg) = if k == ks {
                (next, c.legs[ks][1])
            } else if next == ks {
                (k, c.legs[ks][0])
            } else {
                return None;
            };
            (c.d[band].dot(leg) < 0.0).then_some((ks, band, leg))
        });
        let conic = if let Some((ks, band, leg)) = reflex {
            let (q, rail) = if k == ks { (p0, p1) } else { (p1, p0) };
            let w = width[band];
            let back = -c.d[band];
            let near_q = q + (leg - c.d[ks]) * (w / 2.0);
            let past = rail + back * ((rail - q).dot(-back) + w);
            let poles = if k == ks {
                vec![q, near_q, past, rail]
            } else {
                vec![rail, past, near_q, q]
            };
            let knots =
                ogeom_math::KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 3)?;
            ogeom_geom::BSplineCurve::new(knots, poles, tol)?
        } else {
            let (t0, t1) = (into((k, p0), (next, p1))?, into((next, p1), (k, p0))?);
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
            ogeom_geom::BSplineCurve::rational(knots, control)?
        };
        let samples: Vec<Point> = (0..=16)
            .map(|i| conic.point_at(f64::from(i) / 16.0, tol))
            .collect::<OgeomResult<_>>()?;
        let conic = Curve::BSpline(conic);
        let at_sharp = c.sharp == Some(k) || c.sharp == Some(next);
        let trim = crate::sharp::trim_edges(model, &conic, (&v0, &v1), at_sharp, tol)?;
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
        let mut ring = trim.clone();
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
        if let Some((arc, face)) = &arcs[k] {
            boundary.push(FillBoundary {
                edge: arc.clone(),
                support: Some(face.clone()),
                continuity: Continuity::G1,
            });
        }
        if let [whole] = trim.as_slice() {
            boundary.push(FillBoundary {
                edge: whole.clone(),
                support: Some(filler.clone()),
                continuity: Continuity::G1,
            });
        }
        trims.push(crate::sharp::Trim {
            curve: conic,
            edges: trim,
            normal: outward,
            filler: filler.clone(),
        });
        fillers.push(filler);
    }
    let mut faces = kept;
    faces.extend(fillers);
    let mut history: History = history;
    // The patches are held inside what the solid and its concave rounds
    // hold near the vertex: out of each concave round's ball.
    let region = |p: Point| signed_region(c, p);
    let Some(ks) = c.sharp else {
        let filled = make_filling_n_region(model, &boundary, &[], fit, &region, tol)?;
        let patch = filled.built.shape.clone();
        signed_patch_within(model, c, &[&patch], tol)?;
        faces.push(patch.clone());
        let sewn = ogeom_algo::sew(model, &faces, tol)?;
        if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
            ogeom_bail!(Construction, "the signed corner's faces did not close");
        }
        let rounded = model.add_solid(std::slice::from_ref(&sewn.shells[0]))?;
        history.generate(&c.vertex, patch);
        return Ok(Built::new(rounded, history));
    };
    // The opening runs into the sharp edge's point: a patch there turning
    // its normal from one of the edge's faces to the other's, and one over
    // the rest.
    let (r1, r2) = ((ks + 1) % 3, (ks + 2) % 3);
    let (Some(arc1), Some(arc2)) = (&arcs[r1], &arcs[r2]) else {
        ogeom_bail!(Construction, "a rounded edge has no end arc");
    };
    let patches = crate::sharp::close_through_q(
        model,
        [&trims[ks], &trims[r1], &trims[r2]],
        [arc1, arc2],
        [true; 5],
        &crate::sharp::Held {
            // A signed corner is not convex (the bracket's front is reflex
            // at the vertex), so its planes do not bound the solid there.
            within: &[],
            // The sharp edge is convex: the solid holds the piece within its
            // two faces' planes there, and the arcs of the patch at `Q` are
            // kept inside it.
            seam: &[(c.normal[(ks + 2) % 3], c.corner), (c.normal[ks], c.corner)],
            region: Some(&region),
        },
        tol,
    )?;
    signed_patch_within(model, c, &[&patches.0, &patches.1], tol)?;
    crate::sharp::sew_corner(model, faces, patches, &c.vertex, history, tol)
}

/// Refuses by name where a signed corner's patches stand outside what the
/// solid and its concave rounds hold near the vertex, by more than the
/// bands' fit target, proved over every point by [`certify_within`]. A convex
/// round only takes material off and a concave one adds no more than the
/// fill between its faces and its ball, so a patch past both stands where
/// no round of the request puts material.
fn signed_patch_within(
    model: &Model,
    c: &SignedCorner,
    patches: &[&Shape],
    tol: Tolerances,
) -> OgeomResult<()> {
    let fit = crate::marched::band_fit_target(tol);
    for patch in patches {
        let certified = certify_within(model, c, patch, fit, tol)?;
        if let Err((out, at)) = certified {
            ogeom_bail!(
                NotDone,
                "the signed corner's patch stands {out} outside what the solid and its \
                 concave rounds hold, at {at:?}"
            );
        }
    }
    Ok(())
}

/// The worst any vertex of `patch`'s mesh, drawn at a quarter of the
/// band's fit target, stands outside [`signed_region`], and where: the sampled check the
/// certified one replaced, kept to cross-check it.
fn mesh_outside(
    model: &Model,
    c: &SignedCorner,
    patch: &Shape,
    tol: Tolerances,
) -> OgeomResult<(f64, Point)> {
    let deflection = ogeom_mesh::Deflection {
        chord: crate::marched::band_fit_target(tol) * 0.25,
        ..ogeom_mesh::Deflection::default()
    };
    let mesh = ogeom_mesh::triangulate_face(model, patch, deflection, tol)?;
    let mut worst = (f64::NEG_INFINITY, Point::ORIGIN);
    for p in &mesh.positions {
        let out = signed_outside(c, *p);
        if out > worst.0 {
            worst = (out, *p);
        }
    }
    Ok(worst)
}

/// How finely a span of a patch is split before a part of it no bound
/// certifies is taken as standing outside.
const CERTIFY_DEPTH: u32 = 12;

/// Whether every point of `patch` lies within `fit` of [`signed_region`]:
/// `Ok` when proved, or a point of the patch that stands farther out and
/// by how much (a point the proof could not clear, where none is found).
///
/// The patch's surface is a B-spline whose weights are all positive, so
/// over each span it lies in the convex hull of that span's control
/// points, and over any rectangle of its parameters in the hull of the
/// control points of its segment there. Each piece of the region bounds
/// itself over a hull: a plane's height by its largest at the hull's
/// points, a concave round's ball by the hull's distance from its axis;
/// and the region's maxima and minima of pieces bound by the same of the
/// pieces' bounds. A rectangle whose bound is within `fit` is cleared; one
/// wholly outside the face's trim is passed over; the rest is split in
/// four, down to [`CERTIFY_DEPTH`] halvings of a span, a point inside the
/// trim read at each split.
fn certify_within(
    model: &Model,
    c: &SignedCorner,
    patch: &Shape,
    fit: f64,
    tol: Tolerances,
) -> OgeomResult<Result<(), (f64, Point)>> {
    let surface = model
        .node(patch)
        .and_then(|n| n.data().as_face())
        .and_then(|d| model.geometry().surface(d.surface));
    let Some(ogeom_geom::SurfaceGeometry::BSpline(surface)) = surface else {
        // Not a B-spline: read on the mesh as before.
        let (out, at) = mesh_outside(model, c, patch, tol)?;
        return Ok(if out > fit { Err((out, at)) } else { Ok(()) });
    };
    let grid = surface.grid();
    if grid.points().iter().any(|w| w.weight <= 0.0) {
        let (out, at) = mesh_outside(model, c, patch, tol)?;
        return Ok(if out > fit { Err((out, at)) } else { Ok(()) });
    }
    let place = patch.transform(model.datums())?;
    // The trim in the surface's parameters, and how far its polylines may
    // stand from the trim itself.
    let chord = crate::marched::band_fit_target(tol) * 0.25;
    let rings = ogeom_mesh::face_boundary(
        model,
        patch,
        ogeom_mesh::Deflection {
            chord,
            ..ogeom_mesh::Deflection::default()
        },
        tol,
    )?;
    let slack = chord * 16.0;
    let index = TrimIndex::of(&rings);
    let near_trim = |lo: ogeom_math::Point2, hi: ogeom_math::Point2| -> bool {
        let centre = ogeom_math::Point2::new(f64::midpoint(lo.x, hi.x), f64::midpoint(lo.y, hi.y));
        let half = (hi.x - lo.x).hypot(hi.y - lo.y) * 0.5;
        index.near(centre, half + slack) || index.inside(centre)
    };
    let (u_knots, v_knots) = (surface.u_knots(), surface.v_knots());
    let (pu, pv) = (u_knots.degree(), v_knots.degree());
    let (un, vn) = (grid.u_count(), grid.v_count());
    let points: Vec<Point> = grid
        .points()
        .iter()
        .map(|w| place.apply(w.point()))
        .collect();
    let rational = grid.points().iter().any(|w| w.weight != 1.0);
    if rational {
        // Weighted: each rectangle's segment, its hull that of its
        // projected control points, the weights all positive.
        let us: Vec<f64> = u_knots.distinct().iter().map(|k| k.0).collect();
        let vs: Vec<f64> = v_knots.distinct().iter().map(|k| k.0).collect();
        for w in us.windows(2) {
            for z in vs.windows(2) {
                let rect = (
                    ogeom_math::Point2::new(w[0], z[0]),
                    ogeom_math::Point2::new(w[1], z[1]),
                );
                if let Err(found) =
                    certify_segment(surface, c, &place, rect, fit, &near_trim, &rings, 0, tol)?
                {
                    return Ok(Err(found));
                }
            }
        }
        return Ok(Ok(()));
    }
    // Polynomial: the patch's Bezier pieces, along u for each column of
    // control points and then along v for each row of those.
    let mut columns: Vec<Vec<ogeom_math::bspline::BezierSegment<Point>>> = Vec::with_capacity(vn);
    for j in 0..vn {
        let column: Vec<Point> = (0..un).map(|i| points[i * vn + j]).collect();
        columns.push(ogeom_math::bspline::to_bezier_segments(
            u_knots, &column, tol,
        )?);
    }
    // Each column of pieces along u on its own worker; the first piece
    // found standing outside, in the pieces' order, whatever the count.
    let spans: Vec<usize> = (0..columns[0].len()).collect();
    let found = ogeom_core::parallel::map_ordered(&spans, |_, &s_u| {
        let (ua, ub) = columns[0][s_u].0;
        let mut rows: Vec<Vec<ogeom_math::bspline::BezierSegment<Point>>> =
            Vec::with_capacity(pu + 1);
        for k in 0..=pu {
            let row: Vec<Point> = columns.iter().map(|col| col[s_u].1[k]).collect();
            rows.push(ogeom_math::bspline::to_bezier_segments(v_knots, &row, tol)?);
        }
        for s_v in 0..rows[0].len() {
            let (va, vb) = rows[0][s_v].0;
            let net: Vec<Point> = rows.iter().flat_map(|r| r[s_v].1.iter().copied()).collect();
            let rect = (
                ogeom_math::Point2::new(ua, va),
                ogeom_math::Point2::new(ub, vb),
            );
            let piece = Bezier { net, pu, pv };
            if let Err(found) = certify_bezier(
                surface, c, &place, &piece, rect, fit, &near_trim, &rings, 0, tol,
            )? {
                return Ok(Err(found));
            }
        }
        Ok(Ok(()))
    });
    for result in found {
        if let Err(found) = result? {
            return Ok(Err(found));
        }
    }
    Ok(Ok(()))
}

/// A Bezier patch's control net, `pv + 1` points along v for each of its
/// `pu + 1` rows along u.
struct Bezier {
    net: Vec<Point>,
    pu: usize,
    pv: usize,
}

impl Bezier {
    /// The two halves split at the middle of u, then of v: four pieces, in
    /// the order of their rectangles' corners (low u low v, high u low v,
    /// low u high v, high u high v).
    fn quarters(&self) -> [Self; 4] {
        let (pu, pv) = (self.pu, self.pv);
        let at = |net: &[Point], i: usize, j: usize| net[i * (pv + 1) + j];
        // de Casteljau at a half along u, for each column.
        let mut low = vec![Point::ORIGIN; self.net.len()];
        let mut high = vec![Point::ORIGIN; self.net.len()];
        for j in 0..=pv {
            let mut line: Vec<Point> = (0..=pu).map(|i| at(&self.net, i, j)).collect();
            for r in 0..=pu {
                low[r * (pv + 1) + j] = line[0];
                high[(pu - r) * (pv + 1) + j] = line[pu - r];
                for k in 0..pu - r {
                    line[k] = line[k].midpoint(line[k + 1]);
                }
            }
        }
        let split_v = |net: &[Point]| -> (Vec<Point>, Vec<Point>) {
            let mut a = vec![Point::ORIGIN; net.len()];
            let mut b = vec![Point::ORIGIN; net.len()];
            for i in 0..=pu {
                let mut line: Vec<Point> = (0..=pv).map(|j| net[i * (pv + 1) + j]).collect();
                for r in 0..=pv {
                    a[i * (pv + 1) + r] = line[0];
                    b[i * (pv + 1) + pv - r] = line[pv - r];
                    for k in 0..pv - r {
                        line[k] = line[k].midpoint(line[k + 1]);
                    }
                }
            }
            (a, b)
        };
        let (ll, lh) = split_v(&low);
        let (hl, hh) = split_v(&high);
        [ll, hl, lh, hh].map(|net| Self { net, pu, pv })
    }
}

/// [`certify_within`] over one Bezier piece of a polynomial patch, over
/// rectangle `(lo, hi)` of its parameters: cleared by its net's hull, or
/// split in four by de Casteljau.
#[allow(clippy::too_many_arguments, reason = "one piece of a recursion")]
fn certify_bezier(
    surface: &ogeom_geom::BSplineSurface,
    c: &SignedCorner,
    place: &ogeom_math::Transform,
    piece: &Bezier,
    (lo, hi): (ogeom_math::Point2, ogeom_math::Point2),
    fit: f64,
    near_trim: &dyn Fn(ogeom_math::Point2, ogeom_math::Point2) -> bool,
    rings: &[Vec<ogeom_math::Point2>],
    depth: u32,
    tol: Tolerances,
) -> OgeomResult<Result<(), (f64, Point)>> {
    if !near_trim(lo, hi) {
        return Ok(Ok(()));
    }
    let bound = signed_bound(c, &piece.net);
    if bound <= fit {
        return Ok(Ok(()));
    }
    // A point of the face here that stands past `fit` settles it.
    let mid = ogeom_math::Point2::new(f64::midpoint(lo.x, hi.x), f64::midpoint(lo.y, hi.y));
    let p = place.apply(ogeom_geom::Surface::point_at(surface, mid.x, mid.y, tol)?);
    if ogeom_mesh::inside_boundary(rings, mid) {
        let out = signed_outside(c, p);
        if out > fit {
            return Ok(Err((out, p)));
        }
    }
    if depth >= CERTIFY_DEPTH {
        return Ok(Err((bound, p)));
    }
    let quarters = piece.quarters();
    let rects = [
        (lo, mid),
        (
            ogeom_math::Point2::new(mid.x, lo.y),
            ogeom_math::Point2::new(hi.x, mid.y),
        ),
        (
            ogeom_math::Point2::new(lo.x, mid.y),
            ogeom_math::Point2::new(mid.x, hi.y),
        ),
        (mid, hi),
    ];
    for (quarter, rect) in quarters.iter().zip(rects) {
        if let Err(found) = certify_bezier(
            surface,
            c,
            place,
            quarter,
            rect,
            fit,
            near_trim,
            rings,
            depth + 1,
            tol,
        )? {
            return Ok(Err(found));
        }
    }
    Ok(Ok(()))
}

/// [`certify_within`] over one rectangle of a weighted patch, split in four
/// while its segment's hull does not clear it.
#[allow(clippy::too_many_arguments, reason = "one rectangle of a recursion")]
fn certify_segment(
    surface: &ogeom_geom::BSplineSurface,
    c: &SignedCorner,
    place: &ogeom_math::Transform,
    (lo, hi): (ogeom_math::Point2, ogeom_math::Point2),
    fit: f64,
    near_trim: &dyn Fn(ogeom_math::Point2, ogeom_math::Point2) -> bool,
    rings: &[Vec<ogeom_math::Point2>],
    depth: u32,
    tol: Tolerances,
) -> OgeomResult<Result<(), (f64, Point)>> {
    if !near_trim(lo, hi) {
        return Ok(Ok(()));
    }
    let piece = surface.segment((lo.x, hi.x), (lo.y, hi.y), tol)?;
    let hull: Vec<Point> = piece
        .grid()
        .points()
        .iter()
        .map(|w| place.apply(w.point()))
        .collect();
    let bound = signed_bound(c, &hull);
    if bound <= fit {
        return Ok(Ok(()));
    }
    let mid = ogeom_math::Point2::new(f64::midpoint(lo.x, hi.x), f64::midpoint(lo.y, hi.y));
    let p = place.apply(ogeom_geom::Surface::point_at(surface, mid.x, mid.y, tol)?);
    if ogeom_mesh::inside_boundary(rings, mid) {
        let out = signed_outside(c, p);
        if out > fit {
            return Ok(Err((out, p)));
        }
    }
    if depth >= CERTIFY_DEPTH {
        return Ok(Err((bound, p)));
    }
    for rect in [
        (lo, mid),
        (
            ogeom_math::Point2::new(mid.x, lo.y),
            ogeom_math::Point2::new(hi.x, mid.y),
        ),
        (
            ogeom_math::Point2::new(lo.x, mid.y),
            ogeom_math::Point2::new(mid.x, hi.y),
        ),
        (mid, hi),
    ] {
        if let Err(found) = certify_segment(
            surface,
            c,
            place,
            rect,
            fit,
            near_trim,
            rings,
            depth + 1,
            tol,
        )? {
            return Ok(Err(found));
        }
    }
    Ok(Ok(()))
}

/// A face's trim rings, bucketed: their segments on a grid of cells over
/// their bounds, for the nearest within a reach, and by bands of height,
/// for the crossings a point's ray to `+x` makes.
struct TrimIndex {
    segments: Vec<(ogeom_math::Point2, ogeom_math::Point2)>,
    lo: ogeom_math::Point2,
    hi: ogeom_math::Point2,
    step: (f64, f64),
    cells: Vec<Vec<usize>>,
    bands: Vec<Vec<usize>>,
}

/// The cells of a [`TrimIndex`] each way.
const TRIM_CELLS: usize = 64;

impl TrimIndex {
    fn of(rings: &[Vec<ogeom_math::Point2>]) -> Self {
        let mut segments = Vec::new();
        for ring in rings {
            for k in 0..ring.len() {
                segments.push((ring[k], ring[(k + 1) % ring.len()]));
            }
        }
        let (mut lo, mut hi) = (
            ogeom_math::Point2::new(f64::INFINITY, f64::INFINITY),
            ogeom_math::Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        );
        for (a, b) in &segments {
            for p in [a, b] {
                lo = ogeom_math::Point2::new(lo.x.min(p.x), lo.y.min(p.y));
                hi = ogeom_math::Point2::new(hi.x.max(p.x), hi.y.max(p.y));
            }
        }
        #[expect(clippy::cast_precision_loss, reason = "a cell count of sixty-four")]
        let count = TRIM_CELLS as f64;
        let step = (
            ((hi.x - lo.x) / count).max(f64::MIN_POSITIVE),
            ((hi.y - lo.y) / count).max(f64::MIN_POSITIVE),
        );
        let mut index = Self {
            segments,
            lo,
            hi,
            step,
            cells: vec![Vec::new(); TRIM_CELLS * TRIM_CELLS],
            bands: vec![Vec::new(); TRIM_CELLS],
        };
        for (k, (a, b)) in index.segments.iter().enumerate() {
            let (i0, j0) = index.cell(ogeom_math::Point2::new(a.x.min(b.x), a.y.min(b.y)));
            let (i1, j1) = index.cell(ogeom_math::Point2::new(a.x.max(b.x), a.y.max(b.y)));
            for i in i0..=i1 {
                for j in j0..=j1 {
                    index.cells[i * TRIM_CELLS + j].push(k);
                }
            }
            for j in j0..=j1 {
                index.bands[j].push(k);
            }
        }
        index
    }

    /// The cell holding `p`, clamped onto the grid.
    fn cell(&self, p: ogeom_math::Point2) -> (usize, usize) {
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "a cell index, clamped into the grid"
        )]
        let at = |t: f64, lo: f64, step: f64| -> usize {
            (((t - lo) / step).floor().max(0.0) as usize).min(TRIM_CELLS - 1)
        };
        (
            at(p.x, self.lo.x, self.step.0),
            at(p.y, self.lo.y, self.step.1),
        )
    }

    /// Whether a segment comes within `reach` of `p`: those of the cells
    /// the square of side `2 reach` about `p` meets, every segment
    /// within `reach` lying in one of them.
    fn near(&self, p: ogeom_math::Point2, reach: f64) -> bool {
        let (i0, j0) = self.cell(ogeom_math::Point2::new(p.x - reach, p.y - reach));
        let (i1, j1) = self.cell(ogeom_math::Point2::new(p.x + reach, p.y + reach));
        (i0..=i1).any(|i| {
            (j0..=j1).any(|j| {
                self.cells[i * TRIM_CELLS + j].iter().any(|&k| {
                    let (a, b) = self.segments[k];
                    segment_distance(p, a, b) <= reach
                })
            })
        })
    }

    /// Whether `p` lies inside the rings, by the parity of the crossings
    /// its ray to `+x` makes: only segments of its band of height can
    /// cross it.
    fn inside(&self, p: ogeom_math::Point2) -> bool {
        if p.y < self.lo.y || p.y > self.hi.y {
            return false;
        }
        let (_, j) = self.cell(p);
        let mut crossings = 0usize;
        for &k in &self.bands[j] {
            let (a, b) = self.segments[k];
            if (a.y > p.y) != (b.y > p.y) {
                let x = a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x);
                if x > p.x {
                    crossings += 1;
                }
            }
        }
        crossings % 2 == 1
    }
}

/// The distance from `p` to the segment from `a` to `b`.
fn segment_distance(p: ogeom_math::Point2, a: ogeom_math::Point2, b: ogeom_math::Point2) -> f64 {
    let (ux, uy) = (b.x - a.x, b.y - a.y);
    let length = ux.mul_add(ux, uy * uy);
    let t = if length > 0.0 {
        (((p.x - a.x) * ux + (p.y - a.y) * uy) / length).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (p.x - (a.x + ux * t)).hypot(p.y - (a.y + uy * t))
}

/// An upper bound of [`signed_outside`] over the convex hull of `points`:
/// each piece of the region bounded over the hull, a plane's height by its
/// largest at the points (it is affine), a concave round's ball by the
/// hull's least distance from its axis, and combined as the region
/// combines the pieces, a maximum's bound the maximum of bounds and a
/// minimum's the minimum.
fn signed_bound(c: &SignedCorner, points: &[Point]) -> f64 {
    let xs: Vec<Vector> = points.iter().map(|p| *p - c.corner).collect();
    let top = |n: Vector| {
        xs.iter()
            .map(|x| n.dot(*x))
            .fold(f64::NEG_INFINITY, f64::max)
    };
    let h = |f: usize| top(c.normal[f]);
    let concave: Vec<usize> = (0..3)
        .filter(|&k| c.legs[k][1].dot(c.normal[(k + 2) % 3]) >= 0.0)
        .collect();
    let solid = match concave.as_slice() {
        [k] => h((k + 1) % 3).max(h((k + 2) % 3).min(h(*k))),
        [k1, k2] => {
            let faces = |k: usize| [(k + 2) % 3, k];
            let shared = (0..3)
                .find(|f| faces(*k1).contains(f) && faces(*k2).contains(f))
                .unwrap_or(0);
            let others: Vec<usize> = (0..3).filter(|&f| f != shared).collect();
            h(shared).min(h(others[0]).max(h(others[1])))
        }
        _ => h(0).max(h(1)).max(h(2)),
    };
    let mut out = solid;
    for &k in &concave {
        let r = c.radius[k];
        if r <= 0.0 {
            continue;
        }
        let (na, nb, d) = (c.normal[(k + 2) % 3], c.normal[k], c.d[k]);
        let centre = (na + nb) * (r / (1.0 + na.dot(nb)));
        let contact = centre - na * r;
        let towards = na + nb;
        let towards = towards / towards.magnitude();
        // The hull seen along the edge: its least distance from the ball's
        // axis, nought where the axis passes through it.
        let (e1, e2) = crate::rolled::basis(d);
        let flat: Vec<(f64, f64)> = xs.iter().map(|x| (e1.dot(*x), e2.dot(*x))).collect();
        let axis = (e1.dot(centre), e2.dot(centre));
        let ball = r - hull_distance(&flat, axis);
        let fill = top(-na)
            .max(top(-nb))
            .max(ball)
            .max(top(towards) - towards.dot(contact))
            .max(top(-d));
        out = out.min(fill);
    }
    out
}

/// The distance from `p` to the convex hull of `points` in the plane:
/// nought where `p` lies inside it.
fn hull_distance(points: &[(f64, f64)], p: (f64, f64)) -> f64 {
    let mut sorted: Vec<(f64, f64)> = points.to_vec();
    sorted.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    sorted.dedup();
    if sorted.len() == 1 {
        return (sorted[0].0 - p.0).hypot(sorted[0].1 - p.1);
    }
    let cross = |o: (f64, f64), a: (f64, f64), b: (f64, f64)| {
        (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
    };
    let mut hull: Vec<(f64, f64)> = Vec::with_capacity(2 * sorted.len());
    for &q in sorted.iter().chain(sorted.iter().rev().skip(1)) {
        while hull.len() >= 2 && cross(hull[hull.len() - 2], hull[hull.len() - 1], q) <= 0.0 {
            hull.pop();
        }
        hull.push(q);
    }
    hull.pop();
    let n = hull.len();
    let edge = |a: (f64, f64), b: (f64, f64)| {
        segment_distance(
            ogeom_math::Point2::new(p.0, p.1),
            ogeom_math::Point2::new(a.0, a.1),
            ogeom_math::Point2::new(b.0, b.1),
        )
    };
    if n < 3 {
        return (0..n)
            .map(|k| edge(hull[k], hull[(k + 1) % n]))
            .fold(f64::INFINITY, f64::min);
    }
    // Counter-clockwise: inside where `p` is left of every edge.
    if (0..n).all(|k| cross(hull[k], hull[(k + 1) % n], p) >= 0.0) {
        return 0.0;
    }
    (0..n)
        .map(|k| edge(hull[k], hull[(k + 1) % n]))
        .fold(f64::INFINITY, f64::min)
}

/// How far `p` stands outside the solid near the corner's vertex together
/// with each concave round's fill there (negative inside): the vertex's
/// three face half-spaces, met across a convex edge and joined across a
/// concave one, and each concave edge's fill between its two faces, off its
/// ball and short of the chord between the ball's contacts, along the edge
/// from the vertex.
fn signed_outside(c: &SignedCorner, p: Point) -> f64 {
    signed_region(c, p).0
}

/// [`signed_outside`] with the outward unit normal of the boundary piece
/// that decides it at `p`: a face's plane, a concave round's ball, the
/// chord between the ball's contacts, or the plane square to the edge
/// through the vertex.
fn signed_region(c: &SignedCorner, p: Point) -> (f64, Vector) {
    let x = p - c.corner;
    let h = |f: usize| (c.normal[f].dot(x), c.normal[f]);
    let most = |a: (f64, Vector), b: (f64, Vector)| if b.0 > a.0 { b } else { a };
    let least = |a: (f64, Vector), b: (f64, Vector)| if b.0 < a.0 { b } else { a };
    // Edge `k` lies between face `k - 1` and face `k`; it is convex where
    // its leg into face `k` runs inside face `k - 1`.
    let concave: Vec<usize> = (0..3)
        .filter(|&k| c.legs[k][1].dot(c.normal[(k + 2) % 3]) >= 0.0)
        .collect();
    let solid = match concave.as_slice() {
        [k] => most(h((k + 1) % 3), least(h((k + 2) % 3), h(*k))),
        [k1, k2] => {
            let faces = |k: usize| [(k + 2) % 3, k];
            let shared = (0..3)
                .find(|f| faces(*k1).contains(f) && faces(*k2).contains(f))
                .unwrap_or(0);
            let others: Vec<usize> = (0..3).filter(|&f| f != shared).collect();
            least(h(shared), most(h(others[0]), h(others[1])))
        }
        _ => most(most(h(0), h(1)), h(2)),
    };
    let mut out = solid;
    for &k in &concave {
        let r = c.radius[k];
        if r <= 0.0 {
            continue;
        }
        let (na, nb, d) = (c.normal[(k + 2) % 3], c.normal[k], c.d[k]);
        let centre = (na + nb) * (r / (1.0 + na.dot(nb)));
        let across = x - d * d.dot(x);
        let contact = centre - na * r;
        let towards = na + nb;
        let towards = towards / towards.magnitude();
        let off = across - centre;
        let reach = off.magnitude();
        let ball = if reach > 0.0 {
            (r - reach, off * (-1.0 / reach))
        } else {
            (r, -towards)
        };
        let fill = most(
            most(
                most((-na.dot(x), -na), (-nb.dot(x), -nb)),
                most(ball, (towards.dot(across - contact), towards)),
            ),
            (-d.dot(x), -d),
        );
        out = least(out, fill);
    }
    out
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

    use super::{SignedCorner, certify_within, mesh_outside, signed_corners};
    use ogeom_core::Tolerances;
    use ogeom_geom::{Curve3d as _, SurfaceGeometry};
    use ogeom_math::{Direction, Frame, Point};
    use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    fn edge_near(m: &Model, s: &Shape, q: Point) -> Shape {
        let d = |e: &Shape| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            (0..=64)
                .map(|i| {
                    c.point_at(a + (b - a) * f64::from(i) / 64.0, T)
                        .unwrap()
                        .distance(q)
                })
                .fold(f64::INFINITY, f64::min)
        };
        explore_unique(m, s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .min_by(|a, b| d(a).total_cmp(&d(b)))
            .unwrap()
    }

    /// The L-bracket's corner rounded at `radius` on `edges` (of the front
    /// edge, the re-entrant edge and the upright's front edge), its signed
    /// corner, and the fitted faces the rounding closed it with.
    fn corner(radius: f64, three: bool) -> (Model, SignedCorner, Vec<Shape>) {
        let mut m = Model::new();
        let frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
        let a = ogeom_algo::make_box(&mut m, frame, (20.0, 20.0, 5.0), T)
            .unwrap()
            .shape;
        let b = ogeom_algo::make_box(&mut m, frame, (5.0, 20.0, 20.0), T)
            .unwrap()
            .shape;
        let l = ogeom_bool::fuse(&mut m, &a, &b, T).unwrap().shape;
        let mut edges = vec![
            edge_near(&m, &l, Point::new(12.5, 0.0, 5.0)),
            edge_near(&m, &l, Point::new(5.0, 10.0, 5.0)),
        ];
        if three {
            edges.push(edge_near(&m, &l, Point::new(5.0, 0.0, 12.5)));
        }
        let request: Vec<(Shape, f64)> = edges.iter().map(|e| (e.clone(), radius)).collect();
        let c = signed_corners(&m, &l, &request, T)
            .unwrap()
            .into_iter()
            .next()
            .expect("the bracket's signed corner");
        let rounded = crate::fillet_edges(&mut m, &l, &edges, radius, T)
            .unwrap()
            .shape;
        let patches: Vec<Shape> = explore_unique(&m, &rounded, ShapeType::Face)
            .unwrap()
            .into_iter()
            .filter(|f| {
                let data = m.node(f).unwrap().data().as_face().unwrap();
                matches!(
                    m.geometry().surface(data.surface),
                    Some(SurfaceGeometry::BSpline(_))
                ) && ogeom_algo::shape_bounds(&m, f, T)
                    .unwrap()
                    .distance_to(Point::new(5.0, 0.0, 5.0))
                    < 2.0 * radius
            })
            .collect();
        assert!(!patches.is_empty(), "the corner's patches");
        (m, c, patches)
    }

    /// The certified check proves what the fine mesh it replaced samples:
    /// on the bracket's corners, with and without the upright's edge left
    /// sharp, both clear every patch at the bands' fit target; asked to
    /// keep the patch a tenth of that target inside the region, which the
    /// patch's sides on the faces cannot, the mesh finds a point past it
    /// and the proof clears nothing the mesh does not.
    #[test]
    fn the_certified_check_agrees_with_the_fine_mesh() {
        let fit = crate::marched::band_fit_target(T);
        for (radius, three) in [(2.0, false), (2.0, true)] {
            let (m, c, patches) = corner(radius, three);
            for patch in &patches {
                let sampled = mesh_outside(&m, &c, patch, T).unwrap();
                let proved = certify_within(&m, &c, patch, fit, T).unwrap();
                assert!(sampled.0 <= fit, "the mesh finds {} outside", sampled.0);
                assert!(proved.is_ok(), "the proof refuses {proved:?}");
                let strict = -0.1 * fit;
                let proved = certify_within(&m, &c, patch, strict, T).unwrap();
                assert!(sampled.0 > strict, "the mesh finds every point well inside");
                assert!(proved.is_err(), "the proof clears what the mesh does not");
            }
        }
    }
}
