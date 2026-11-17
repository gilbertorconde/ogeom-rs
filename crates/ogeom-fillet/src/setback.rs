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

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Continuity, Curve, Curve3d as _};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

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
    #[cfg(feature = "old-path-log")]
    let _old_path_request = crate::old_path::public_request("round_vertex_setback");
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
    crate::engine::vertex_setbacks::guard(model, solid, vertex, &frame, &setback, tol)?;
    Ok(SetbackCorner {
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

    cornered_by_engine(
        model,
        solid,
        &corner.edges,
        radius,
        &[(vertex.clone(), setbacks.clone())],
        std::slice::from_ref(&corner),
        tol,
    )?
    .ok_or_else(|| {
        ogeom_core::OgeomError::Construction(
            "the blend engine could not close this setback corner".into(),
        )
    })
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
    #[cfg(feature = "old-path-log")]
    let _old_path_request = crate::old_path::public_request("fillet_edges_cornered");
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
    cornered_by_engine(model, solid, edges, radius, corners, &read, tol)?.ok_or_else(|| {
        ogeom_core::OgeomError::Construction(
            "the blend engine could not close these setback corners".into(),
        )
    })
}

/// [`cornered_own`] by the blend engine: the rounds of `edges` at `radius`
/// with each corner of `corners` not the ball's own set back (`read`, the
/// corners as read off the solid); `None` hands the request back.
fn cornered_by_engine(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
    corners: &[(Shape, Setbacks)],
    read: &[SetbackCorner],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut once: Vec<(Shape, f64)> = Vec::with_capacity(edges.len());
    for edge in edges {
        if !once.iter().any(|(e, _)| e.is_same(edge)) {
            once.push((edge.clone(), radius));
        }
    }
    if once.is_empty() || !(radius.is_finite() && radius > tol.confusion()) {
        return Ok(None);
    }
    let mut request = crate::engine::Request::from_rounds(&once, tol)?;
    request.corners = corners
        .iter()
        .zip(read)
        .filter(|(_, corner)| !corner.is_the_balls(tol))
        .map(|((vertex, setbacks), _)| (vertex.clone(), setbacks.clone()))
        .collect();
    let baked = ogeom_algo::baked_where_placed(model, solid, tol)?;
    if !baked.shape.is_same(solid) {
        let trace = |part: &Shape| -> OgeomResult<Shape> {
            let [only] = baked.history.trace(part) else {
                ogeom_bail!(Construction, "a placed setback request has no unique image");
            };
            Ok(only.clone())
        };
        for blend in &mut request.blends {
            blend.edge = trace(&blend.edge)?;
        }
        for (vertex, setbacks) in &mut request.corners {
            *vertex = trace(vertex)?;
            for (edge, _) in &mut setbacks.per_edge {
                *edge = trace(edge)?;
            }
        }
    }
    Ok(
        crate::engine::run(model, &baked.shape, &request, tol)?.map(|built| Built {
            shape: built.shape,
            history: baked.history.then(&built.history),
        }),
    )
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

pub(crate) use crate::engine::vertex_setbacks::guard_cell;

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
) -> OgeomResult<Vec<Shape>> {
    let mut corners = Vec::new();
    for (vertex, at) in signed_vertices(model, solid, request, tol)? {
        if runs_through(model, solid, &vertex, &at, tol)?
            || crate::cusp::cusp_corner(model, solid, &vertex, &at, tol)?
        {
            continue;
        }
        signed_corner(model, solid, &vertex, &at, tol)?;
        corners.push(vertex);
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
) -> OgeomResult<Vec<Shape>> {
    // Read in node order, so the corner and its sides are the same in any
    // order asked.
    let mut sorted = request.to_vec();
    sorted.sort_by_key(|(edge, _)| edge.node());
    let mut corners = Vec::new();
    for (vertex, at) in signed_vertices(model, solid, &sorted, tol)? {
        if crate::cusp::cusp_corner(model, solid, &vertex, &at, tol)? {
            corners.push(vertex);
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
) -> OgeomResult<()> {
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
    Ok(())
}
