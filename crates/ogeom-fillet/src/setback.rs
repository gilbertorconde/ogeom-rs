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

use crate::fill_n::{FillBoundary, make_filling_n};

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
fn point_of(model: &Model, vertex: &Shape) -> OgeomResult<Point> {
    let Some(raw) = model
        .node(vertex)
        .and_then(|n| n.data().as_vertex().map(|d| d.point))
    else {
        ogeom_bail!(Construction, "the vertex holds no point");
    };
    Ok(vertex.transform(model.datums())?.apply(raw))
}

/// [`round_vertex_setback`] on a shape every node below which is its own.
#[allow(clippy::too_many_lines, reason = "one corner, assembled end to end")]
fn setback_own(
    model: &mut Model,
    solid: &Shape,
    vertex: &Shape,
    radius: f64,
    setbacks: &Setbacks,
    tol: Tolerances,
) -> OgeomResult<Built> {
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

    // The solid's edge along each of the corner's edge directions.
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

    // The bands, meeting in the ball's corner: all that a setback equal to
    // the ball's own asks for.
    let banded = crate::fillet::fillet_edges(model, solid, &edges, radius, tol)?;
    if (0..n).all(|k| setback[k] <= own_setback[k] + tol.confusion()) {
        return Ok(banded);
    }
    close_corner(
        model,
        banded,
        vertex,
        &frame,
        &setback,
        setbacks.continuity,
        tol,
    )
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

/// The bands of `banded` meeting at `vertex`, whose corner reads as
/// `frame`, stopped at `setback` along their edges and the opening between
/// them filled by one patch at `continuity`.
fn close_corner(
    model: &mut Model,
    banded: Built,
    vertex: &Shape,
    frame: &crate::corner::CornerFrame,
    setback: &[f64],
    continuity: Continuity,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let ends = vec![CornerEnd::Round; frame.d.len()];
    close_corner_ends(
        model, banded, vertex, frame, setback, &ends, continuity, tol,
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
        Ok((0..n).find(|&k| {
            outward.cross(sections[k].0).magnitude() <= tol.angular() * 10.0
                && on_plane(origin, sections[k])
        }))
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
            arcs.push(bevel_end(model, &kept, frame, *section, tol)?);
            continue;
        }
        let mut found: Option<(Shape, Shape)> = None;
        for face in &kept {
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
                    .all(|p| on_plane(p, *section));
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
        arcs.push(found);
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
        let conic = Curve::BSpline(ogeom_geom::BSplineCurve::rational(knots, control)?);
        let trim = ogeom_algo::make_edge_between(model, conic, (0.0, 1.0), &v0, &v1, tol)?.shape;
        let vq = ogeom_algo::make_vertex(model, q).shape;
        let ring = vec![
            trim.clone(),
            crate::support::segment_between(model, (&v1, p1), (&vq, q), tol)?,
            crate::support::segment_between(model, (&vq, q), (&v0, p0), tol)?,
        ];
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
        let filler = ogeom_algo::make_face_with_pcurves(model, surface, &[ring], tol)?.shape;
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

    let filled = make_filling_n(model, &boundary, &[], fit, tol)?;
    let patch = filled.built.shape.clone();
    let mut faces = kept;
    faces.extend(fillers);
    faces.push(patch.clone());
    let sewn = ogeom_algo::sew(model, &faces, tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the setback corner's faces did not close");
    }
    let rounded = model.add_solid(std::slice::from_ref(&sewn.shells[0]))?;
    let mut history: History = history;
    history.generate(vertex, patch);
    Ok(Built::new(rounded, history))
}

/// A bevel's straight end on its section plane `section`, and the bevel's
/// face: the planar face of `kept`, on none of the corner's host planes,
/// with a straight edge lying on the section plane.
fn bevel_end(
    model: &Model,
    kept: &[Shape],
    frame: &crate::corner::CornerFrame,
    section: (Vector, Point),
    tol: Tolerances,
) -> OgeomResult<(Shape, Shape)> {
    let on_plane = |p: Point| section.0.dot(p - section.1).abs() <= tol.confusion() * 100.0;
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
    close_corner(model, banded, vertex, &frame, &setback, Continuity::G1, tol)
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
            if far >= setback[k] {
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
