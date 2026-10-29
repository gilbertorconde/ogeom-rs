//! Trimmed runs of a pipe glued on the planes they meet on.
//!
//! Two runs trimmed on the plane they meet on each leave a face there: the
//! region of the plane its own solid covers. Their union has no face where
//! both regions cover the plane, a face looking out of the first run where
//! only its region does, and one looking out of the second where only the
//! second's does. Those faces are built here from the two regions' own
//! boundary curves, split where the curves cross, and the runs' other
//! faces sewn to them, so no boolean is asked to reconcile two boundaries
//! drawn apart.
//!
//! A boundary curve of one region that runs within `gap` of the other
//! region's boundary all along is taken for that boundary: the strip
//! between them is no wider than the runs' skins may stand off the solids
//! they are held to, and the two curves are sewn as one edge whose
//! tolerance spans them.

use ogeom_algo::{
    make_edge_between, make_face_with_pcurves, make_solid, project_on_curve, sew_within,
    shape_bounds,
};
use ogeom_core::{OgeomResult, Tolerance, Tolerances, ogeom_bail};
use ogeom_geom::{
    Curve, Curve3d as _, PlaneSurface, SurfaceGeometry, Transformable as _, TrimmedCurve,
};
use ogeom_intersect::{CurveCurveOptions, intersect_curves};
use ogeom_math::{Direction, Frame, Plane, Point, Point2, Vector};
use ogeom_topo::{Filter, Model, NodeData, Shape, ShapeType, VertexData, explore};

use super::spine_curve_of;

/// Where two trimmed runs meet: a point of the plane and its normal,
/// looking out of run `a` into run `b`.
pub(super) struct Meeting {
    pub(super) a: usize,
    pub(super) b: usize,
    pub(super) at: Point,
    pub(super) normal: Vector,
    /// How far apart two boundary curves on the plane may run and be
    /// taken for one.
    pub(super) gap: f64,
}

/// Samples along a boundary curve where it is measured against the other
/// region's boundary.
const MEASURED: usize = 64;

/// Segments each boundary curve is seeded with when two are crossed.
const SEEDS: usize = 512;

/// A boundary curve of one region on the meeting plane, walked with the
/// region on its left about the plane's normal.
#[derive(Clone)]
struct Bound {
    curve: Curve,
    from: f64,
    to: f64,
    /// Whether it bounds the first run's region.
    first: bool,
    /// The curve's polyline in the plane's coordinates, from `from` to `to`.
    dense: Vec<Point2>,
}

impl Bound {
    fn at(&self, t: f64, tol: Tolerances) -> OgeomResult<Point> {
        self.curve.point_at(t, tol)
    }

    fn start(&self, tol: Tolerances) -> OgeomResult<Point> {
        self.at(self.from, tol)
    }

    fn end(&self, tol: Tolerances) -> OgeomResult<Point> {
        self.at(self.to, tol)
    }

    fn reversed(&self) -> Self {
        let mut dense = self.dense.clone();
        dense.reverse();
        Self {
            curve: self.curve.clone(),
            from: self.to,
            to: self.from,
            first: self.first,
            dense,
        }
    }

    /// The stretch from `from` to `to` of this bound's curve.
    fn piece(&self, from: f64, to: f64, frame: &Frame, tol: Tolerances) -> OgeomResult<Self> {
        Ok(Self {
            curve: self.curve.clone(),
            from,
            to,
            first: self.first,
            dense: dense_between(&self.curve, from, to, frame, tol)?,
        })
    }

    /// The distance from `p` to this bound's curve between its ends.
    fn distance(&self, p: Point, tol: Tolerances) -> OgeomResult<f64> {
        let found = project_on_curve(&self.curve, p, MEASURED, tol)?;
        let (lo, hi) = (self.from.min(self.to), self.from.max(self.to));
        let inside =
            found.parameter >= lo - tol.parametric() && found.parameter <= hi + tol.parametric();
        let ends = p.distance(self.start(tol)?).min(p.distance(self.end(tol)?));
        Ok(if inside {
            found.distance.min(ends)
        } else {
            ends
        })
    }
}

/// The trimmed runs `solids` united where they meet on each of `meetings`,
/// each run touching the next only on those planes. Boundaries of the two
/// regions on a plane that run within the meeting's `gap` of each other
/// are taken for one.
///
/// # Errors
///
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where a run has
/// no face on a plane it meets another on, where the two regions'
/// boundaries neither keep within `gap` of each other nor cross cleanly,
/// or where the faces do not sew into one closed shell.
pub(super) fn glued(
    model: &mut Model,
    solids: &[Shape],
    meetings: &[Meeting],
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let mut on_planes: Vec<Shape> = Vec::new();
    let mut built: Vec<Shape> = Vec::new();
    let mut sew_gap = tol.confusion();
    for meeting in meetings {
        let normal = Direction::new(meeting.normal, tol)?;
        let frame = Frame::about(meeting.at, normal);
        // Two meetings of the same runs on one plane (the two corners of
        // a symmetric lens) are one arrangement, taken whole at the first.
        let unused = |faces: Vec<Shape>| -> Vec<Shape> {
            faces
                .into_iter()
                .filter(|f| !on_planes.iter().any(|g| g.is_same(f)))
                .collect()
        };
        let first = unused(faces_on(model, &solids[meeting.a], &frame, tol)?);
        let second = unused(faces_on(model, &solids[meeting.b], &frame, tol)?);
        if first.is_empty() && second.is_empty() {
            continue;
        }
        if first.is_empty() || second.is_empty() {
            ogeom_bail!(
                NotDone,
                "a run has no face on the plane it meets the next on"
            );
        }
        let mut bounds = bounds_of(model, &first, true, &frame, tol)?;
        bounds.extend(bounds_of(model, &second, false, &frame, tol)?);
        let (faces, spanned) = between(model, &bounds, &frame, meeting.gap, tol)?;
        sew_gap = sew_gap.max(spanned);
        on_planes.extend(first);
        on_planes.extend(second);
        built.extend(faces);
    }
    let mut faces = built;
    for solid in solids {
        for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
            if !on_planes.iter().any(|f| f.is_same(&face)) {
                faces.push(face);
            }
        }
    }
    let sewn = sew_within(model, &faces, sew_gap, tol)?;
    if sewn.shells.is_empty() || !sewn.free_edges.is_empty() {
        ogeom_bail!(
            NotDone,
            "the runs glued on their meeting planes did not close"
        );
    }
    // A section with a hole swept round a closed spine bounds a void as
    // well: the outer shell is the one reaching farthest, and the others
    // stand inside it.
    let mut shells = sewn.shells;
    let reach = |model: &Model, shell: &Shape| -> OgeomResult<f64> {
        let b = shape_bounds(model, shell, tol)?;
        Ok(match (b.low(), b.high()) {
            (Some(low), Some(high)) => low.distance(high),
            _ => 0.0,
        })
    };
    let mut widest = 0;
    for k in 1..shells.len() {
        if reach(model, &shells[k])? > reach(model, &shells[widest])? {
            widest = k;
        }
    }
    shells.swap(0, widest);
    Ok(make_solid(model, &shells)?.shape)
}

/// The faces of `solid` on the plane of `frame`.
fn faces_on(
    model: &Model,
    solid: &Shape,
    frame: &Frame,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let mut found = Vec::new();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data()) else {
            continue;
        };
        let Some(SurfaceGeometry::Plane(p)) = model.geometry().surface(data.surface) else {
            continue;
        };
        let plane = p.plane();
        if plane
            .normal()
            .vector()
            .cross(frame.z().vector())
            .magnitude()
            <= tol.angular()
            && plane.distance_to(frame.origin()) <= tol.confusion() * 100.0
            && !found.iter().any(|f: &Shape| f.is_same(&face))
        {
            found.push(face);
        }
    }
    Ok(found)
}

/// The boundary curves of `faces`, each walked with its face on its left
/// about `frame`'s `z`: the ring of a face enclosing the most area
/// anticlockwise, its others clockwise.
fn bounds_of(
    model: &Model,
    faces: &[Shape],
    first: bool,
    frame: &Frame,
    tol: Tolerances,
) -> OgeomResult<Vec<Bound>> {
    let mut out = Vec::new();
    for face in faces {
        let mut rings: Vec<(f64, Vec<Bound>)> = Vec::new();
        for wire in model.ordered_children_of(face)? {
            let mut ring = Vec::new();
            for edge in explore(model, &wire, Filter::OfType(ShapeType::Edge))? {
                let (curve, range) = spine_curve_of(model, &edge)?;
                let curve = curve.transformed(&edge.transform(model.datums())?, tol)?;
                let (from, to) = if edge.orientation() == ogeom_topo::Orientation::Reversed {
                    (range.1, range.0)
                } else {
                    range
                };
                ring.push(Bound {
                    dense: dense_between(&curve, from, to, frame, tol)?,
                    curve,
                    from,
                    to,
                    first,
                });
            }
            let area = signed_area(ring.iter().map(|b| b.dense.as_slice()));
            rings.push((area, ring));
        }
        let Some(outer) =
            (0..rings.len()).max_by(|&i, &j| rings[i].0.abs().total_cmp(&rings[j].0.abs()))
        else {
            continue;
        };
        for (k, (area, ring)) in rings.into_iter().enumerate() {
            let flip = (k == outer) == (area < 0.0);
            if flip {
                out.extend(ring.iter().rev().map(Bound::reversed));
            } else {
                out.extend(ring);
            }
        }
    }
    Ok(out)
}

/// The faces of the meeting plane that bound the union of the two regions
/// `bounds` encloses, and the widest gap between two boundary curves taken
/// for one.
fn between(
    model: &mut Model,
    bounds: &[Bound],
    frame: &Frame,
    gap: f64,
    tol: Tolerances,
) -> OgeomResult<(Vec<Shape>, f64)> {
    // A bound is shared where every sample along it lies within `gap` of
    // the other region's boundary and the two run the same way there; a
    // shared bound of one region needs the other's to lie within `gap` of
    // its own shared ones, so the shared bounds are pared until both sides
    // agree.
    let mut shared: Vec<bool> = vec![true; bounds.len()];
    let mut spanned;
    loop {
        let mut changed = false;
        spanned = 0.0_f64;
        for i in 0..bounds.len() {
            if !shared[i] {
                continue;
            }
            let others: Vec<&Bound> = bounds
                .iter()
                .zip(&shared)
                .filter(|(b, s)| **s && b.first != bounds[i].first)
                .map(|(b, _)| b)
                .collect();
            match off_by(&bounds[i], &others, gap, tol)? {
                Some(off) => spanned = spanned.max(off),
                None => {
                    shared[i] = false;
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    // Shared bounds meet free ones at ends the two regions place apart by
    // up to the gap between their shared bounds.
    let join = spanned.max(tol.confusion() * 10.0);

    // Each free bound split where it crosses a free bound of the other
    // region, or where one of those ends on it. A shared bound of the other
    // region it may meet only at that bound's ends.
    let free: Vec<usize> = (0..bounds.len()).filter(|&i| !shared[i]).collect();
    let mut cuts: Vec<Vec<f64>> = vec![Vec::new(); bounds.len()];
    for &i in &free {
        for j in 0..bounds.len() {
            let (a, b) = (&bounds[i], &bounds[j]);
            if a.first == b.first || (!shared[j] && j < i) {
                continue;
            }
            let found = intersect_curves(
                &trimmed(a, tol)?,
                &trimmed(b, tol)?,
                CurveCurveOptions {
                    samples: SEEDS,
                    gap: tol.confusion(),
                },
                tol,
            )?;
            if !found.overlaps.is_empty() {
                ogeom_bail!(
                    NotDone,
                    "two regions' boundaries on a meeting plane share a stretch"
                );
            }
            for crossing in found.crossings {
                if shared[j] {
                    let p = crossing.point;
                    if p.distance(b.start(tol)?) > join && p.distance(b.end(tol)?) > join {
                        ogeom_bail!(
                            NotDone,
                            "a region's boundary on a meeting plane crosses the other's \
                             where the two are taken for one"
                        );
                    }
                    continue;
                }
                cuts[i].push(crossing.on_a);
                cuts[j].push(crossing.on_b);
            }
            if !shared[j] {
                for (on, end) in [
                    (i, b.start(tol)?),
                    (i, b.end(tol)?),
                    (j, a.start(tol)?),
                    (j, a.end(tol)?),
                ] {
                    let found = project_on_curve(&bounds[on].curve, end, MEASURED, tol)?;
                    if found.distance <= join {
                        cuts[on].push(found.parameter);
                    }
                }
            }
        }
    }

    // Each piece of a free bound lies inside the other region or out of
    // it. The face looking out of the first run is walked by the first
    // region's pieces outside the second and the second's inside the
    // first, reversed; the face looking out of the second alike.
    let mut up: Vec<Bound> = Vec::new();
    let mut down: Vec<Bound> = Vec::new();
    for &i in &free {
        let bound = &bounds[i];
        let ts = cuts_of(bound, &cuts[i], join, tol)?;
        let mut stops = vec![bound.from];
        stops.extend(ts);
        stops.push(bound.to);
        let other: Vec<&Bound> = bounds.iter().filter(|b| b.first != bound.first).collect();
        for w in stops.windows(2) {
            let piece = bound.piece(w[0], w[1], frame, tol)?;
            let mid = bound.at(f64::midpoint(w[0], w[1]), tol)?;
            let clear = other
                .iter()
                .map(|b| b.distance(mid, tol))
                .try_fold(f64::INFINITY, |m, d| d.map(|d| m.min(d)))?;
            // Ten times what the polylines the winding is read on stray.
            if clear <= tol.confusion() * 100.0 {
                ogeom_bail!(
                    NotDone,
                    "a piece of a region's boundary on a meeting plane runs along the other's"
                );
            }
            let inside = winding(&other, to_plane(frame, mid)) != 0;
            match (bound.first, inside) {
                (true, false) => up.push(piece),
                (false, true) => up.push(piece.reversed()),
                (false, false) => down.push(piece),
                (true, true) => down.push(piece.reversed()),
            }
        }
    }
    let mut corners = Vec::new();
    let mut faces = faces_of(model, &up, frame.z(), frame, join, &mut corners, tol)?;
    faces.extend(faces_of(
        model,
        &down,
        -frame.z(),
        frame,
        join,
        &mut corners,
        tol,
    )?);
    Ok((faces, spanned))
}

/// How far `bound` lies off `others` where every sample along it lies
/// within `gap` of them and the nearest runs its way, or `None`.
fn off_by(bound: &Bound, others: &[&Bound], gap: f64, tol: Tolerances) -> OgeomResult<Option<f64>> {
    if others.is_empty() {
        return Ok(None);
    }
    let mut widest = 0.0_f64;
    for s in 0..=MEASURED {
        #[allow(clippy::cast_precision_loss, reason = "a small count")]
        let f = s as f64 / MEASURED as f64;
        let t = bound.from + (bound.to - bound.from) * f;
        let p = bound.at(t, tol)?;
        let mut best = (f64::INFINITY, 0usize);
        for (k, other) in others.iter().enumerate() {
            let d = other.distance(p, tol)?;
            if d < best.0 {
                best = (d, k);
            }
        }
        if best.0 > gap {
            return Ok(None);
        }
        widest = widest.max(best.0);
        if s == MEASURED / 2
            && direction(bound, t, tol)?.dot(heading_near(others[best.1], p, tol)?) <= 0.0
        {
            return Ok(None);
        }
    }
    // The ends, which sewing joins, lie within the gap of the other
    // region's ends.
    for p in [bound.start(tol)?, bound.end(tol)?] {
        let mut nearest = f64::INFINITY;
        for other in others {
            nearest = nearest
                .min(p.distance(other.start(tol)?))
                .min(p.distance(other.end(tol)?));
        }
        if nearest > gap {
            return Ok(None);
        }
        widest = widest.max(nearest);
    }
    Ok(Some(widest))
}

/// The unit tangent of `bound` at `t`, the way it is walked.
fn direction(bound: &Bound, t: f64, tol: Tolerances) -> OgeomResult<Vector> {
    let d = bound.curve.d1_at(t, tol)?;
    let sense = if bound.to >= bound.from { 1.0 } else { -1.0 };
    Ok(d * sense)
}

/// The tangent of `bound` the way it is walked, at its point nearest `p`.
fn heading_near(bound: &Bound, p: Point, tol: Tolerances) -> OgeomResult<Vector> {
    let found = project_on_curve(&bound.curve, p, MEASURED, tol)?;
    let (lo, hi) = (bound.from.min(bound.to), bound.from.max(bound.to));
    direction(bound, found.parameter.clamp(lo, hi), tol)
}

/// The parameters in `cuts` strictly inside `bound`, more than `join` from
/// its ends and from each other, in its walking order.
fn cuts_of(bound: &Bound, cuts: &[f64], join: f64, tol: Tolerances) -> OgeomResult<Vec<f64>> {
    let (lo, hi) = (bound.from.min(bound.to), bound.from.max(bound.to));
    let mut ts: Vec<f64> = cuts
        .iter()
        .copied()
        .filter(|t| *t > lo && *t < hi)
        .collect();
    ts.sort_by(f64::total_cmp);
    if bound.to < bound.from {
        ts.reverse();
    }
    let mut kept: Vec<f64> = Vec::new();
    let mut last = bound.start(tol)?;
    for t in ts {
        let p = bound.at(t, tol)?;
        if p.distance(last) > join && p.distance(bound.end(tol)?) > join {
            kept.push(t);
            last = p;
        }
    }
    Ok(kept)
}

/// `bound`'s curve over its range, as a curve of its own.
fn trimmed(bound: &Bound, tol: Tolerances) -> OgeomResult<Curve> {
    let (lo, hi) = (bound.from.min(bound.to), bound.from.max(bound.to));
    Ok(TrimmedCurve::new(bound.curve.clone(), lo, hi, tol)?.into())
}

/// The planar faces normal to `normal` that the pieces `pieces` walk round,
/// each ring anticlockwise about the plane's `z` round the region it
/// encloses and clockwise round a hole.
fn faces_of(
    model: &mut Model,
    pieces: &[Bound],
    normal: Direction,
    frame: &Frame,
    join: f64,
    corners: &mut Vec<(Point, Shape)>,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let rings = rings_of(pieces, join, tol)?;
    let mut outers: Vec<(f64, Vec<Bound>)> = Vec::new();
    let mut holes: Vec<Vec<Bound>> = Vec::new();
    for ring in rings {
        let area = signed_area(ring.iter().map(|b| b.dense.as_slice()));
        if area > 0.0 {
            outers.push((area, ring));
        } else {
            holes.push(ring);
        }
    }
    let mut held: Vec<Vec<Vec<Bound>>> = vec![Vec::new(); outers.len()];
    for hole in holes {
        let probe = hole[0].dense[0];
        let Some(k) = (0..outers.len())
            .filter(|&k| winding(&outers[k].1.iter().collect::<Vec<_>>(), probe) != 0)
            .min_by(|&i, &j| outers[i].0.total_cmp(&outers[j].0))
        else {
            ogeom_bail!(NotDone, "a hole on a meeting plane lies in no face");
        };
        held[k].push(hole);
    }
    // Facing against the plane's `z`, each ring is walked back.
    let back = normal.vector().dot(frame.z().vector()) < 0.0;
    let walked = |ring: Vec<Bound>| -> Vec<Bound> {
        if back {
            ring.iter().rev().map(Bound::reversed).collect()
        } else {
            ring
        }
    };
    let surface: SurfaceGeometry =
        PlaneSurface::new(Plane::new(Frame::about(frame.origin(), normal))).into();
    let mut faces = Vec::new();
    for ((_, outer), holes) in outers.into_iter().zip(held) {
        let mut rings = vec![edges_of(model, &walked(outer), join, corners, tol)?];
        for hole in holes {
            rings.push(edges_of(model, &walked(hole), join, corners, tol)?);
        }
        faces.push(make_face_with_pcurves(model, surface.clone(), &rings, tol)?.shape);
    }
    Ok(faces)
}

/// The pieces chained end to start into closed rings, each piece followed
/// by the one starting nearest its end, within `join`.
fn rings_of(pieces: &[Bound], join: f64, tol: Tolerances) -> OgeomResult<Vec<Vec<Bound>>> {
    let starts: Vec<Point> = pieces
        .iter()
        .map(|p| p.start(tol))
        .collect::<OgeomResult<_>>()?;
    let ends: Vec<Point> = pieces
        .iter()
        .map(|p| p.end(tol))
        .collect::<OgeomResult<_>>()?;
    let mut used = vec![false; pieces.len()];
    let mut rings = Vec::new();
    for first in 0..pieces.len() {
        if used[first] {
            continue;
        }
        used[first] = true;
        let mut ring = vec![first];
        loop {
            let last = ring[ring.len() - 1];
            if ends[last].distance(starts[first]) <= join {
                break;
            }
            let next = (0..pieces.len())
                .filter(|&k| !used[k] && ends[last].distance(starts[k]) <= join)
                .min_by(|&i, &j| {
                    ends[last]
                        .distance(starts[i])
                        .total_cmp(&ends[last].distance(starts[j]))
                });
            let Some(next) = next else {
                ogeom_bail!(
                    NotDone,
                    "the boundary pieces on a meeting plane do not close"
                );
            };
            used[next] = true;
            ring.push(next);
        }
        rings.push(ring.into_iter().map(|k| pieces[k].clone()).collect());
    }
    Ok(rings)
}

/// The edges of a closed ring through `ring`'s pieces, one vertex where
/// each meets the next, its tolerance spanning both pieces' ends there.
/// Rings meeting at a point (two faces either side of a crossing) meet on
/// one vertex: `corners` holds those made so far.
fn edges_of(
    model: &mut Model,
    ring: &[Bound],
    join: f64,
    corners: &mut Vec<(Point, Shape)>,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let n = ring.len();
    let mut vertices = Vec::with_capacity(n);
    for k in 0..n {
        let (arrive, leave) = (ring[(k + n - 1) % n].end(tol)?, ring[k].start(tol)?);
        let at = arrive + (leave - arrive) * 0.5;
        if let Some((_, v)) = corners.iter().find(|(p, _)| p.distance(at) <= join) {
            let v = v.clone();
            let reach = arrive
                .distance(placed(model, &v)?)
                .max(leave.distance(placed(model, &v)?));
            if reach > tol.confusion() {
                model.widen(&v, Tolerance::new(reach)?)?;
            }
            vertices.push(v);
            continue;
        }
        let reach = (arrive.distance(leave) * 0.5).max(tol.confusion());
        let data = VertexData::with_tolerance(at, reach)?;
        let v = model.add_vertex(data);
        corners.push((at, v.clone()));
        vertices.push(v);
    }
    let mut edges = Vec::with_capacity(n);
    for k in 0..n {
        let piece = &ring[k];
        let (from, to) = (&vertices[k], &vertices[(k + 1) % n]);
        let edge = if piece.to > piece.from {
            make_edge_between(
                model,
                piece.curve.clone(),
                (piece.from, piece.to),
                from,
                to,
                tol,
            )?
            .shape
        } else {
            make_edge_between(
                model,
                piece.curve.clone(),
                (piece.to, piece.from),
                to,
                from,
                tol,
            )?
            .shape
            .reversed()
        };
        edges.push(edge);
    }
    Ok(edges)
}

/// `curve`'s polyline from `from` to `to` in the plane's coordinates, its
/// chords within ten confusions of the curve.
fn dense_between(
    curve: &Curve,
    from: f64,
    to: f64,
    frame: &Frame,
    tol: Tolerances,
) -> OgeomResult<Vec<Point2>> {
    const FIRST: u32 = 64;
    let at = |f: f64| curve.point_at(from + (to - from) * f, tol);
    let mut dense = Vec::new();
    for i in 0..FIRST {
        let (f0, f1) = (
            f64::from(i) / f64::from(FIRST),
            f64::from(i + 1) / f64::from(FIRST),
        );
        super::chords_within(
            &at,
            (f0, at(f0)?),
            (f1, at(f1)?),
            tol.confusion() * 10.0,
            0,
            &mut dense,
        )?;
    }
    dense.push(at(1.0)?);
    Ok(dense.into_iter().map(|p| to_plane(frame, p)).collect())
}

fn to_plane(frame: &Frame, p: Point) -> Point2 {
    let local = frame.to_local(p);
    Point2::new(local.x, local.y)
}

/// Twice the area enclosed by the polylines walked end to end, positive
/// anticlockwise.
fn signed_area<'a>(lines: impl Iterator<Item = &'a [Point2]>) -> f64 {
    let mut sum = 0.0;
    for line in lines {
        for w in line.windows(2) {
            sum += w[0].x * w[1].y - w[1].x * w[0].y;
        }
    }
    sum
}

/// How many times the bounds wind round `p`, anticlockwise counted up.
fn winding(bounds: &[&Bound], p: Point2) -> i32 {
    let mut turned = 0.0;
    for b in bounds {
        for w in b.dense.windows(2) {
            let (u, v) = ((w[0].x - p.x, w[0].y - p.y), (w[1].x - p.x, w[1].y - p.y));
            turned += (u.0 * v.1 - u.1 * v.0).atan2(u.0 * v.0 + u.1 * v.1);
        }
    }
    #[allow(clippy::cast_possible_truncation, reason = "a whole count of turns")]
    let turns = (turned / core::f64::consts::TAU).round() as i32;
    turns
}

/// Where a vertex stands.
fn placed(model: &Model, vertex: &Shape) -> OgeomResult<Point> {
    let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
        ogeom_bail!(Construction, "a vertex holds no point");
    };
    Ok(vertex.transform(model.datums())?.apply(data.point))
}
