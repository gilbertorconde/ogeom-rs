//! The trimmed-face assembly of the strips on a curved solid.

use super::*;

/// P3, P4 and P5 for strips on a solid with curved faces, whatever section
/// law made them: the sheet of the trimmed faces and the strips, its shell
/// and solid with its seams joined, each strip proved within the material
/// its two faces bound, and the history; `None` at a vertex outside the
/// cases built, or a strip found past its material.
///
/// # Errors
///
/// As building the sheet.
#[cfg_attr(
    not(test),
    allow(
        dead_code,
        reason = "the round sections on curved seats assemble through it from M3 on"
    )
)]
pub(crate) fn assemble(
    model: &mut Model,
    curved: &Curved,
    bevels: &[Strip],
    blended: &FastSet<Id>,
    corners: &FastMap<(usize, Id), Corner>,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    assemble_consuming(
        model,
        curved,
        bevels,
        blended,
        corners,
        &FastSet::default(),
        tol,
    )
}

/// [`assemble`] where the strips consume the faces of `consumed` whole:
/// each such face is not rebuilt, its strips meet each other along their
/// contacts on it (given as one curve, the ridge), its corners are where
/// that ridge stands at each vertex's station, and its sharp edges go.
///
/// # Errors
///
/// As [`assemble`].
pub(crate) fn assemble_consuming(
    model: &mut Model,
    curved: &Curved,
    bevels: &[Strip],
    blended: &FastSet<Id>,
    corners: &FastMap<(usize, Id), Corner>,
    consumed: &FastSet<usize>,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    assemble_widened(
        model,
        curved,
        bevels,
        blended,
        corners,
        consumed,
        &FastMap::default(),
        false,
        &[],
        tol,
    )
}

/// One run of a rolled round's ball on an edge of the face its contact
/// left, as the assembly takes it: the strip's index, the side whose
/// contact leaves (`bevels[strip].faces[side]`), the crossings where it
/// leaves and comes back, and the run [`super::seats::edge_run`] read.
#[derive(Debug, Clone)]
pub(crate) struct RolledRun {
    pub strip: usize,
    pub side: usize,
    pub enter: super::setbacks::RollOver,
    pub leave: super::setbacks::RollOver,
    pub run: super::seats::EdgeRun,
}

/// [`assemble`] with rolled rounds: each run in `rolled` splits its strip
/// into the pieces either side, each between its faces, and the piece
/// between, the ball resting on the edge (ADR 0001, section 13, question
/// 8). The face the contact leaves loses the edge's piece under the round,
/// and splits where that cuts it apart.
///
/// # Errors
///
/// As [`assemble`].
pub(crate) fn assemble_rolled(
    model: &mut Model,
    curved: &Curved,
    bevels: &[Strip],
    blended: &FastSet<Id>,
    corners: &FastMap<(usize, Id), Corner>,
    rolled: &[RolledRun],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    assemble_widened(
        model,
        curved,
        bevels,
        blended,
        corners,
        &FastSet::default(),
        &FastMap::default(),
        false,
        rolled,
        tol,
    )
}

/// [`assemble_consuming`], each contact keyed by face and edge in `widths`
/// a fitted curve standing within that of its face, its edges widened to
/// it; with `mitre`, each vertex of three bevels or more mitred instead of
/// closed by its facet.
///
/// # Errors
///
/// As [`assemble_consuming`].
#[allow(clippy::too_many_arguments, reason = "the assembly's whole input")]
pub(super) fn assemble_widened(
    model: &mut Model,
    curved: &Curved,
    bevels: &[Strip],
    blended: &FastSet<Id>,
    corners: &FastMap<(usize, Id), Corner>,
    consumed: &FastSet<usize>,
    widths: &FastMap<(usize, Id), f64>,
    mitre: bool,
    rolled: &[RolledRun],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut points = Points {
        points: Vec::new(),
        merge: widths
            .values()
            .fold(tol.confusion(), |m, w| m.max(*w))
            .min(MERGE_CAP),
    };
    let mut edges: Vec<SheetCurve> = Vec::new();
    let mut made: FastMap<Key, usize> = FastMap::default();
    let contact_of = |fi: usize, e: Id| -> Option<&Curve> {
        bevels.iter().find(|b| b.edge == e).and_then(|b| {
            b.faces
                .iter()
                .position(|f| *f == fi)
                .map(|k| &b.contacts[k])
        })
    };
    // What each rolled run cuts, by face and the edge whose contact (or
    // which edge itself) the walk takes: the contact on the face it leaves
    // gone between the crossings, the other face's taken by the run's
    // contact between its feet there, and each edge the ball rests on gone
    // from the face it leaves under the round and cut at the crossings on
    // its others.
    let mut spans: FastMap<(usize, Id), Vec<Span>> = FastMap::default();
    for r in rolled {
        let b = &bevels[r.strip];
        let (left, other) = (b.faces[r.side], b.faces[1 - r.side]);
        let off_left = ogeom_algo::PointDistance::of(model, &curved.faces[left].shape, tol)?;
        let c = &b.contacts[r.side];
        let gone = span_between(
            c,
            (r.enter.at, r.enter.point),
            (r.leave.at, r.leave.point),
            |lo, hi| {
                let p = c.point_at(0.5 * (lo + hi), tol)?;
                Ok(off_left.nearest(p, tol)?.0 > tol.confusion())
            },
            Fate::Drop,
            tol,
        )?;
        spans.entry((left, b.edge)).or_default().push(gone);
        let c = &b.contacts[1 - r.side];
        let (Some(&q0), Some(&q1)) = (r.run.feet.first(), r.run.feet.last()) else {
            return Ok(None);
        };
        // The way round holding the run's middle foot, on a closed contact.
        let middle = r.run.feet[r.run.feet.len() / 2];
        let (t0, t1) = (
            param_on(c, c.domain(), q0, tol)?,
            param_on(c, c.domain(), q1, tol)?,
        );
        let tm = param_on(c, c.domain(), middle, tol)?;
        let period = period_of(c, tol);
        let taken = span_between(
            c,
            (t0, q0),
            (t1, q1),
            |lo, hi| Ok(wrapped(tm, lo, period) <= hi),
            Fate::Replace(r.run.contact.clone()),
            tol,
        )?;
        spans.entry((other, b.edge)).or_default().push(taken);
        let n = r.run.over.len();
        let (Some(first), Some(last)) = (r.run.on_over.first(), r.run.on_over.last()) else {
            return Ok(None);
        };
        for (k, e) in r.run.over.iter().enumerate() {
            let edge = &curved.edges[e];
            // The end of this edge it shares with the chain's edge `j`.
            let toward = |j: usize| -> Option<f64> {
                let next = &curved.edges[&r.run.over[j]];
                if edge.ends.1 == next.ends.0 || edge.ends.1 == next.ends.1 {
                    Some(edge.range.1)
                } else if edge.ends.0 == next.ends.0 || edge.ends.0 == next.ends.1 {
                    Some(edge.range.0)
                } else {
                    None
                }
            };
            let at =
                |t: f64| -> OgeomResult<(f64, Point)> { Ok((t, edge.curve.point_at(t, tol)?)) };
            let (a, z) = if n == 1 {
                ((first.1, r.enter.point), (last.1, r.leave.point))
            } else if k == 0 {
                let Some(t) = toward(1) else {
                    return Ok(None);
                };
                ((first.1, r.enter.point), at(t)?)
            } else if k == n - 1 {
                let Some(t) = toward(n - 2) else {
                    return Ok(None);
                };
                (at(t)?, (last.1, r.leave.point))
            } else {
                (at(edge.range.0)?, at(edge.range.1)?)
            };
            let held = r.run.on_over[r.run.on_over.len() / 2];
            let middle = (held.0 == k).then_some(held.1);
            for fi in faces_of(curved, *e) {
                let fate = if fi == left { Fate::Drop } else { Fate::Keep };
                let period = period_of(&edge.curve, tol);
                let span = span_between(
                    &edge.curve,
                    a,
                    z,
                    |lo, hi| Ok(middle.is_none_or(|m| wrapped(m, lo, period) <= hi)),
                    fate,
                    tol,
                )?;
                spans.entry((fi, *e)).or_default().push(span);
            }
        }
    }
    // Each sharp edge's end at each of its vertices: the corner a blend
    // moved there on it, or the vertex itself. Where its two faces' corners
    // stand apart on it, the end is the one further from the vertex, and
    // the piece on to the other is the vertex's tail on it.
    let mut sharp_end: FastMap<(Id, Id), Point> = FastMap::default();
    let mut tails: FastMap<Id, Vec<(Id, Point, Point)>> = FastMap::default();
    for (e, edge) in &curved.edges {
        let kept: Vec<usize> = faces_of(curved, *e)
            .into_iter()
            .filter(|fi| !consumed.contains(fi))
            .collect();
        if blended.contains(e) || kept.is_empty() {
            continue;
        }
        for (end, v) in [edge.ends.0, edge.ends.1].into_iter().enumerate() {
            let at = curved.points[&v];
            let mut moved: Vec<Point> = Vec::new();
            for fi in kept.iter().copied() {
                if let Some(c) = corner_for(corners, curved, fi, v, *e, end == 0, tol)?
                    && c.point.distance(at) > tol.confusion()
                    && let Some(t) = c.on(*e)
                    && edge.curve.point_at(t, tol)?.distance(c.point) <= tol.confusion()
                    && !moved.iter().any(|m| m.distance(c.point) <= tol.confusion())
                {
                    moved.push(c.point);
                }
            }
            // A face whose corner stays at the vertex runs on to it.
            let stays = kept.len() == 2
                && kept.iter().any(|fi| {
                    corners
                        .get(&(*fi, v))
                        .is_some_and(|c| c.point.distance(at) <= tol.confusion())
                });
            match moved.as_slice() {
                [] => {
                    sharp_end.insert((*e, v), at);
                }
                [p] if stays && edge.ends.0 != edge.ends.1 => {
                    sharp_end.insert((*e, v), *p);
                    tails.entry(v).or_default().push((*e, *p, at));
                }
                [p] => {
                    sharp_end.insert((*e, v), *p);
                }
                [p, q] if edge.ends.0 != edge.ends.1 => {
                    let (far, near) = if p.distance(at) >= q.distance(at) {
                        (*p, *q)
                    } else {
                        (*q, *p)
                    };
                    sharp_end.insert((*e, v), far);
                    tails.entry(v).or_default().push((*e, far, near));
                }
                _ => return Ok(None),
            }
        }
    }
    // A vertex's blends: one, or two running on into each other; each
    // bevel's end there, the segment between its two faces' corners.
    // Each cusp's corners on the side and the wall moved to its lobes'.
    let cusps = cusps_of(curved, bevels, corners, tol)?;
    // A one-sided cusp's sharp edges on the side or the wall (the edge the
    // two share below a top ending there, the wall's seam, the top's edges
    // with them where it runs on past the vertex) run on to the vertex: the
    // S's feet meet them there. Two-sided, the patch spans the vertex and
    // they end on its sides.
    for c in cusps.iter().filter(|c| c.lobes.len() == 1) {
        let [top, side, wall] = c.faces;
        for (e, edge) in &curved.edges {
            let faces = faces_of(curved, *e);
            if (edge.ends.0 == c.v || edge.ends.1 == c.v)
                && !blended.contains(e)
                && faces.iter().any(|f| *f == side || *f == wall)
                && faces.iter().all(|f| *f == top || *f == side || *f == wall)
            {
                sharp_end.insert((*e, c.v), curved.points[&c.v]);
                if let Some(t) = tails.get_mut(&c.v) {
                    t.retain(|(te, _, _)| te != e);
                }
            }
        }
    }
    let reset;
    let corners = if cusps.is_empty() {
        corners
    } else {
        let mut map = corners.clone();
        for c in &cusps {
            let [_, side, wall] = c.faces;
            for (k, lobe) in c.lobes.iter().enumerate() {
                let rim = &bevels[lobe.rim];
                let ring = &rim.contacts[usize::from(rim.faces[1] == wall)];
                let at = if k == 0 { c.v } else { pinched(c.v) };
                let line = &bevels[lobe.line];
                let contact = &line.contacts[usize::from(line.faces[1] == side)];
                map.insert(
                    (side, at),
                    Corner {
                        point: lobe.side,
                        on: vec![(
                            line.edge,
                            param_on(contact, contact.domain(), lobe.side, tol)?,
                        )],
                    },
                );
                map.insert(
                    (wall, at),
                    Corner {
                        point: lobe.wall,
                        on: vec![(rim.edge, param_on(ring, ring.domain(), lobe.wall, tol)?)],
                    },
                );
            }
        }
        reset = map;
        &reset
    };
    let mut sections: FastMap<Id, Vec<(Point, Point, usize)>> = FastMap::default();
    for (bi, b) in bevels.iter().enumerate() {
        let e = &curved.edges[&b.edge];
        for (end, v) in [e.ends.0, e.ends.1].into_iter().enumerate() {
            let (Some(c0), Some(c1)) = (
                corner_for(corners, curved, b.faces[0], v, b.edge, end == 0, tol)?,
                corner_for(corners, curved, b.faces[1], v, b.edge, end == 0, tol)?,
            ) else {
                return Ok(None);
            };
            if c0.point.distance(c1.point) > tol.confusion() {
                sections
                    .entry(v)
                    .or_default()
                    .push((c0.point, c1.point, bi));
            }
        }
    }
    let mut facets: Vec<(Id, Vec<Point>)> = Vec::new();
    // A capped bevel's end: the way its cap faces, where the faces at the
    // vertex all stand square to it (a rim's arc running on into the rest
    // of its circle left sharp).
    let mut facing: FastMap<Id, Vector> = FastMap::default();
    // Each mitred vertex's point, where every bevel there passes.
    let mut mitred: FastMap<Id, Point> = FastMap::default();
    // Each vertex where rounds meet on one ball: its corners in turn, the
    // ball's centre and radius.
    let mut balls: Vec<(Id, Vec<Point>, Point, f64)> = Vec::new();
    // Each vertex where rounds meet on no one ball: its corners in turn and
    // the rounds' end arcs' centres.
    let mut vertex_patches: Vec<(Id, Vec<Point>, Vec<Point>)> = Vec::new();
    let mut keyed: Vec<&Id> = curved.points.keys().collect();
    keyed.sort_unstable();
    for v in keyed {
        if cusps.iter().any(|c| c.v == *v) {
            continue;
        }
        let mut at_v: Vec<Id> = curved
            .edges
            .iter()
            .filter(|(e, c)| blended.contains(*e) && (c.ends.0 == *v || c.ends.1 == *v))
            .map(|(e, _)| *e)
            .collect();
        at_v.sort_unstable();
        at_v.dedup();
        match at_v.len() {
            0 => {}
            1 => {
                // One bevel running on into an edge left sharp along it: its
                // end is capped by the plane square to the edge there,
                // through the bevel's end, the tails of the sharp edges whose
                // faces end apart on it, and the vertex.
                let at = curved.points[v];
                let blended_edge = &curved.edges[&at_v[0]];
                let tangent = |c: &CurvedEdge| -> OgeomResult<Vector> {
                    let t = if c.ends.0 == *v { c.range.0 } else { c.range.1 };
                    c.curve.d1_at(t, tol)
                };
                let ta = tangent(blended_edge)?;
                let mut runs_on = false;
                for (e, c) in &curved.edges {
                    if *e != at_v[0] && (c.ends.0 == *v || c.ends.1 == *v) {
                        let tb = tangent(c)?;
                        runs_on |= ta.cross(tb).magnitude()
                            <= super::setbacks::TANGENT * ta.magnitude() * tb.magnitude();
                    }
                }
                if !runs_on {
                    continue;
                }
                let Some(ends) = sections.get(v) else {
                    continue;
                };
                let [(p, q, bi)] = ends[..] else {
                    continue;
                };
                // The cap faces into the bevel at a convex edge, the corner
                // left past it standing behind the cap; out of it at a
                // concave one, closing the fill's end.
                let leaving = if blended_edge.ends.0 == *v { ta } else { -ta };
                let across: Vector = bevels[bi]
                    .faces
                    .iter()
                    .map(|fi| outward_at(curved, *fi, at))
                    .fold(Vector::ZERO, |sum, n| sum + n);
                let chord = p.midpoint(q) - at;
                facing.insert(
                    *v,
                    if chord.dot(across) < 0.0 {
                        leaving
                    } else {
                        -leaving
                    },
                );
                // From the bevel's end on along the tails to the vertex,
                // the cap's last side back to the bevel's end.
                let close = |a: Point, b: Point| a.distance(b) <= tol.confusion();
                let mut cycle = vec![p, q];
                let mut last = q;
                while !close(last, at) {
                    let next = tails
                        .get(v)
                        .into_iter()
                        .flatten()
                        .find_map(|(_, far, near)| {
                            if close(*far, last) {
                                Some(*near)
                            } else if close(*near, last) && !cycle.iter().any(|c| close(*c, *far)) {
                                Some(*far)
                            } else {
                                None
                            }
                        });
                    match next {
                        Some(n) if !cycle.iter().any(|c| close(*c, n)) => {
                            cycle.push(n);
                            last = n;
                        }
                        _ => break,
                    }
                }
                if !close(last, at) {
                    cycle.push(at);
                }
                facets.push((*v, cycle));
            }
            2 => {
                // Two bevels running on into each other share their end.
                if let Some(ends) = sections.get(v) {
                    let first = ends[0];
                    let same = |(p, q, _): &(Point, Point, usize)| {
                        (p.distance(first.0) <= tol.confusion()
                            && q.distance(first.1) <= tol.confusion())
                            || (p.distance(first.1) <= tol.confusion()
                                && q.distance(first.0) <= tol.confusion())
                    };
                    if !ends.iter().all(same) {
                        return Ok(None);
                    }
                }
            }
            _ => {
                // Three bevels or more, every edge there bevelled by a
                // segment: the facet through the faces' corners, its sides
                // the bevels' ends, where those corners lie on one plane. A
                // face's own seam there (a boss's, from its foot) is no edge
                // to bevel.
                let degree = curved
                    .edges
                    .iter()
                    .filter(|(e, c)| {
                        (c.ends.0 == *v || c.ends.1 == *v) && faces_of(curved, **e).len() != 1
                    })
                    .count();
                let Some(ends) = sections.get(v) else {
                    return Ok(None);
                };
                let lines = ends
                    .iter()
                    .all(|(_, _, bi)| matches!(bevels[*bi].section, Section::Line));
                let arcs = ends
                    .iter()
                    .all(|(_, _, bi)| matches!(bevels[*bi].section, Section::Arc { .. }));
                // Rounds only where every edge there is rounded: no tails.
                if degree < at_v.len() || !(lines || arcs) || (arcs && degree != at_v.len()) {
                    return Ok(None);
                }
                // Round the corner, the bevels' ends and the tails of the
                // sharp edges there.
                let mut round = ends.clone();
                round.extend(
                    tails
                        .get(v)
                        .into_iter()
                        .flatten()
                        .map(|(_, far, near)| (*far, *near, usize::MAX)),
                );
                let Some(cycle) = chained(&round, tol) else {
                    return Ok(None);
                };
                if lines && mitre {
                    let surfaces: Vec<&SurfaceGeometry> =
                        ends.iter().map(|(_, _, bi)| &bevels[*bi].surface).collect();
                    #[expect(clippy::cast_precision_loss, reason = "a few corners")]
                    let start = cycle.iter().fold(Vector::ZERO, |s, p| s + p.to_vector())
                        / cycle.len() as f64;
                    let Some(x) = shared_point(&surfaces, Point::from_vector(start), tol)? else {
                        return Ok(None);
                    };
                    mitred.insert(*v, x);
                    continue;
                }
                if lines {
                    facets.push((*v, cycle));
                    continue;
                }
                // Rounds every one: the sphere where their end arcs are
                // arcs of one ball, each about its spine's point there; the
                // patch through them where they are not.
                let mut centres: Vec<Point> = Vec::with_capacity(ends.len());
                for (p, _, bi) in ends {
                    let Section::Arc { spine } = &bevels[*bi].section else {
                        return Ok(None);
                    };
                    centres.push(spine.point_at(param_on(spine, spine.domain(), *p, tol)?, tol)?);
                }
                let centre = centres[0];
                let radius = cycle[0].distance(centre);
                let one_ball = centres
                    .iter()
                    .all(|c| c.distance(centre) <= tol.confusion())
                    && cycle
                        .iter()
                        .all(|p| (p.distance(centre) - radius).abs() <= tol.confusion());
                if one_ball {
                    balls.push((*v, cycle, centre, radius));
                } else {
                    vertex_patches.push((*v, cycle, centres));
                }
            }
        }
    }
    // Each facet's plane through its corners, facing out, and its corners
    // in turn about that side.
    let mut facet_planes: FastMap<Id, SurfaceGeometry> = FastMap::default();
    let mut facet_rings: Vec<(Id, Vec<Point>)> = Vec::with_capacity(facets.len());
    let mut patches: FastMap<Id, SurfaceGeometry> = FastMap::default();
    let mut polygons: Vec<(Id, Vec<Point>)> = Vec::new();
    let mut middles: FastMap<Id, Point> = FastMap::default();
    for (v, cycle) in &facets {
        let mut out = Vector::ZERO;
        for (e, c) in &curved.edges {
            if c.ends.0 == *v || c.ends.1 == *v {
                for fi in faces_of(curved, *e) {
                    out += outward_at(curved, fi, curved.points[v]);
                }
            }
        }
        let mut newell = Vector::ZERO;
        for k in 0..cycle.len() {
            let (p, q) = (
                cycle[k].to_vector(),
                cycle[(k + 1) % cycle.len()].to_vector(),
            );
            newell += p.cross(q);
        }
        let m = newell.magnitude();
        if m <= tol.confusion() * tol.confusion() {
            return Ok(None);
        }
        if let Some(cap) = facing.get(v)
            && newell.dot(out).abs() <= tol.angular() * m * out.magnitude()
        {
            out = *cap;
        }
        let mut cycle = cycle.clone();
        if newell.dot(out) < 0.0 {
            cycle.reverse();
        }
        let normal = if newell.dot(out) < 0.0 {
            -newell
        } else {
            newell
        } / m;
        if cycle
            .iter()
            .any(|p| (*p - cycle[0]).dot(normal).abs() > tol.confusion())
        {
            // Off one plane: the patch through the corners, its sides the
            // segments between them, where each bevel there is a plane and
            // its end that segment; at four corners the bilinear patch.
            let curved_here = sections[v]
                .iter()
                .any(|(_, _, bi)| !matches!(bevels[*bi].surface, SurfaceGeometry::Plane(_)));
            if curved_here {
                // A curved bevel there: its end is where its surface meets
                // the plane through its two corners and the hole's middle,
                // and the patch takes that curve for its side.
                #[expect(clippy::cast_precision_loss, reason = "a handful of corners")]
                let middle = Point::from_vector(
                    cycle
                        .iter()
                        .fold(Vector::ZERO, |sum, p| sum + p.to_vector())
                        / cycle.len() as f64,
                );
                middles.insert(*v, middle);
                polygons.push((*v, cycle));
            } else if cycle.len() == 4 {
                patches.insert(*v, bilinear(&cycle, tol)?);
                facet_rings.push((*v, cycle));
            } else {
                // Five corners or more: their Wachspress blend, built on
                // the sheet's edges once it stands.
                polygons.push((*v, cycle));
            }
            continue;
        }
        let frame = Frame::new(
            cycle[0],
            Direction::new(normal, tol)?,
            Direction::new(cycle[1] - cycle[0], tol)?,
            tol,
        )?;
        facet_planes.insert(*v, PlaneSurface::new(Plane::new(frame)).into());
        facet_rings.push((*v, cycle));
    }
    // What a strip's end at `v`, from `p` to `q`, lies on: the vertex's
    // facet, or a plane face there holding both, where the strip ends flush
    // on it, or the other strip sharing the end where either is a plane,
    // or the plane through both and along an edge left sharp there, where
    // two strips meet across it mirrored. Two curved strips that are not
    // mirrored there meet on the other strip itself, the section marched.
    let cut_at = |strip: &Strip, v: Id, p: Point, q: Point| -> Option<SurfaceGeometry> {
        if let Some(plane) = facet_planes.get(&v) {
            return Some(plane.clone());
        }
        if let Some(middle) = middles.get(&v) {
            let normal = Direction::new((q - p).cross(*middle - p), tol).ok()?;
            let plane =
                Plane::new(Frame::new(p, normal, Direction::new(q - p, tol).ok()?, tol).ok()?);
            return Some(PlaneSurface::new(plane).into());
        }
        let at_v = || {
            curved
                .edges
                .iter()
                .filter(move |(_, c)| c.ends.0 == v || c.ends.1 == v)
        };
        let flush_plane = at_v()
            .flat_map(|(e, _)| faces_of(curved, *e))
            .filter(|fi| !strip.faces.contains(fi))
            .find_map(|fi| match &curved.faces[fi].surface {
                SurfaceGeometry::Plane(plane) => {
                    let pl = plane.plane();
                    (pl.signed_distance_to(p).abs() <= tol.confusion()
                        && pl.signed_distance_to(q).abs() <= tol.confusion())
                    .then(|| curved.faces[fi].surface.clone())
                }
                _ => None,
            });
        if flush_plane.is_some() {
            return flush_plane;
        }
        // Another strip sharing this end: where either of the two is a
        // plane, they meet along the other cut by it.
        let close = |a: Point, b: Point| a.distance(b) <= tol.confusion();
        let other = sections.get(&v).and_then(|ends| {
            ends.iter().find_map(|(x, y, bi)| {
                let same = (close(*x, p) && close(*y, q)) || (close(*x, q) && close(*y, p));
                (same && bevels[*bi].edge != strip.edge).then_some(&bevels[*bi])
            })
        });
        if let Some(other) = other
            && (matches!(other.surface, SurfaceGeometry::Plane(_))
                || matches!(strip.surface, SurfaceGeometry::Plane(_)))
        {
            return Some(other.surface.clone());
        }
        let at = curved.points[&v];
        let along_edge: Option<SurfaceGeometry> = at_v()
            .filter(|(e, _)| !blended.contains(*e))
            .find_map(|(_, c)| {
                let t = if c.ends.0 == v { c.range.0 } else { c.range.1 };
                let along = c.curve.d1_at(t, tol).ok()?;
                let normal = Direction::new(along.cross(q - p), tol).ok()?;
                let plane =
                    Plane::new(Frame::new(p, normal, Direction::new(q - p, tol).ok()?, tol).ok()?);
                (plane.signed_distance_to(at).abs() <= tol.confusion())
                    .then(|| PlaneSurface::new(plane).into())
            });
        if let Some(other) = other
            && other.surface != strip.surface
        {
            let mirror = along_edge.as_ref().is_some_and(|plane| {
                mirrored(&strip.surface, &other.surface, plane, p, q, tol).unwrap_or(false)
            });
            if !mirror {
                return Some(other.surface.clone());
            }
        }
        along_edge
    };
    // The width each marched section's edge is built to, by its index.
    let marched_widths: core::cell::RefCell<FastMap<usize, f64>> =
        core::cell::RefCell::new(FastMap::default());
    let section = |points: &mut Points,
                   edges: &mut Vec<SheetCurve>,
                   made: &mut FastMap<Key, usize>,
                   strip: &Strip,
                   cut: Option<&SurfaceGeometry>,
                   p: Point,
                   q: Point|
     -> OgeomResult<(usize, bool)> {
        let (i, j) = (points.index(p), points.index(q));
        let key = Key::Section(i.min(j), i.max(j));
        if let Some(k) = made.get(&key) {
            return Ok((*k, edges[*k].ends.0 == i));
        }
        let (curve, range, ends): (Curve, (f64, f64), (usize, usize)) = match &strip.section {
            Section::Line => {
                // Two planes meet in the segment; a plane and a curved
                // surface in the curved one cut by the plane.
                let met = match (&strip.surface, cut) {
                    (_, None) | (SurfaceGeometry::Plane(_), Some(SurfaceGeometry::Plane(_))) => {
                        None
                    }
                    (SurfaceGeometry::Plane(_), Some(cut)) => {
                        met_between(cut, &strip.surface, p, q, tol)?
                    }
                    (surface, Some(cut)) => met_between(surface, cut, p, q, tol)?,
                };
                // With no closed form, the segment where it lies on both
                // (a ruling), else the section marched.
                let met = match (met, cut) {
                    (None, Some(cut))
                        if !(matches!(strip.surface, SurfaceGeometry::Plane(_))
                            && matches!(cut, SurfaceGeometry::Plane(_)))
                            && !segment_on(&[&strip.surface, cut], p, q, tol)? =>
                    {
                        marched_between(&strip.surface, cut, p, q, tol)?.map(|(met, width)| {
                            marched_widths.borrow_mut().insert(edges.len(), width);
                            met
                        })
                    }
                    (met, _) => met,
                };
                match met {
                    Some((curve, range, p_low)) => {
                        (curve, range, if p_low { (i, j) } else { (j, i) })
                    }
                    None => {
                        let (a, b) = (points.points[i.min(j)], points.points[i.max(j)]);
                        let line = LineCurve::segment(a, b, tol)?;
                        let range = line.domain();
                        (line.into(), range, (i.min(j), i.max(j)))
                    }
                }
            }
            Section::Arc { spine } => {
                let (a, b) = (points.points[i], points.points[j]);
                let centre = spine.point_at(param_on(spine, spine.domain(), a, tol)?, tol)?;
                let far = spine.point_at(param_on(spine, spine.domain(), b, tol)?, tol)?;
                if centre.distance(far) > tol.confusion() {
                    // Its corners at two stations: on a torus, the curve
                    // running evenly in both its angles; on a straight
                    // strip, the strip cut by the plane through both (the
                    // face it ends flush on, or the one holding square to
                    // the axis the chord of the section at the first).
                    if let SurfaceGeometry::Torus(torus) = &strip.surface {
                        let (curve, range) = torus_skew(&torus.torus(), a, b, tol)?;
                        edges.push(SheetCurve {
                            ends: (i, j),
                            curve,
                            range,
                        });
                        made.insert(key, edges.len() - 1);
                        return Ok((edges.len() - 1, true));
                    }
                    // A fitted strip (no closed form) ending on a face's
                    // plane: the strip cut by it, marched, its edge built to
                    // how far the march stands off.
                    if matches!(strip.surface, SurfaceGeometry::BSpline(_))
                        && let Some(cut @ SurfaceGeometry::Plane(_)) = cut
                        && let Some(((curve, range, a_low), width)) =
                            marched_between(&strip.surface, cut, a, b, tol)?
                    {
                        marched_widths.borrow_mut().insert(edges.len(), width);
                        let ends = if a_low { (i, j) } else { (j, i) };
                        edges.push(SheetCurve { ends, curve, range });
                        made.insert(key, edges.len() - 1);
                        return Ok((edges.len() - 1, ends.0 == i));
                    }
                    let Curve::Line(_) = spine else {
                        ogeom_core::ogeom_bail!(
                            Construction,
                            "a strip along a curved spine ending skew"
                        );
                    };
                    let normal = match cut {
                        Some(SurfaceGeometry::Plane(plane)) => {
                            Some(plane.plane().normal().vector())
                        }
                        _ => None,
                    };
                    let control = skew_control(centre, far - centre, a, b, normal, tol)?;
                    let (curve, _) = rational_quadratics(&[control], tol)?;
                    edges.push(SheetCurve {
                        ends: (i, j),
                        curve,
                        range: (0.0, 1.0),
                    });
                    made.insert(key, edges.len() - 1);
                    return Ok((edges.len() - 1, true));
                }
                let (x, w) = (a - centre, b - centre);
                let normal = x.cross(w);
                let frame = Frame::new(
                    centre,
                    Direction::new(normal, tol)?,
                    Direction::new(x, tol)?,
                    tol,
                )?;
                let y = normal.cross(x);
                let angle = w.dot(y).atan2(w.dot(x) * normal.magnitude());
                let circle = CircleCurve::new(Circle::new(frame, x.magnitude(), tol)?);
                (circle.into(), (0.0, angle), (i, j))
            }
        };
        edges.push(SheetCurve { ends, curve, range });
        made.insert(key, edges.len() - 1);
        Ok((edges.len() - 1, ends.0 == i))
    };
    // The tail of a sharp edge at `v` from `p` to `q`, where it has one
    // between those two: its edge and whether it runs from `p`.
    let tail_between = |points: &mut Points,
                        edges: &mut Vec<SheetCurve>,
                        made: &mut FastMap<Key, usize>,
                        v: Id,
                        p: Point,
                        q: Point|
     -> OgeomResult<Option<(usize, bool)>> {
        let close = |x: Point, y: Point| x.distance(y) <= tol.confusion();
        let found = tails.get(&v).and_then(|ts| {
            ts.iter()
                .find(|(_, far, near)| {
                    (close(*far, p) && close(*near, q)) || (close(*far, q) && close(*near, p))
                })
                .copied()
        });
        let Some((e, far, near)) = found else {
            return Ok(None);
        };
        let i = points.index(p);
        if let Some(k) = made.get(&Key::Tail(e, v)) {
            return Ok(Some((*k, edges[*k].ends.0 == i)));
        }
        let edge = &curved.edges[&e];
        let (tf, tn) = (
            param_on(&edge.curve, edge.range, far, tol)?,
            param_on(&edge.curve, edge.range, near, tol)?,
        );
        let (f, n) = (points.index(far), points.index(near));
        let (range, ends) = if tf <= tn {
            ((tf, tn), (f, n))
        } else {
            ((tn, tf), (n, f))
        };
        edges.push(SheetCurve {
            ends,
            curve: edge.curve.clone(),
            range,
        });
        made.insert(Key::Tail(e, v), edges.len() - 1);
        Ok(Some((edges.len() - 1, ends.0 == i)))
    };
    // A strip's end at a mitred vertex, from `a` to `z`: on to the
    // vertex's point along its section with the bevel sharing `a`, then on
    // to `z` along its section with the one sharing `z`. `false` where the
    // end stands at no mitred vertex.
    let mitred_end = |points: &mut Points,
                      edges: &mut Vec<SheetCurve>,
                      made: &mut FastMap<Key, usize>,
                      ring: &mut Vec<(usize, bool)>,
                      strip: &Strip,
                      a: Point,
                      z: Point|
     -> OgeomResult<bool> {
        let close = |x: Point, y: Point| x.distance(y) <= tol.confusion();
        let found = mitred.iter().find(|(v, _)| {
            sections.get(*v).is_some_and(|ends| {
                ends.iter().any(|(p, q, bi)| {
                    bevels[*bi].edge == strip.edge
                        && ((close(*p, a) && close(*q, z)) || (close(*p, z) && close(*q, a)))
                })
            })
        });
        let Some((v, x)) = found else {
            return Ok(false);
        };
        let beside = |at: Point| {
            sections[v]
                .iter()
                .find(|(p, q, bi)| {
                    bevels[*bi].edge != strip.edge && (close(*p, at) || close(*q, at))
                })
                .map(|(_, _, bi)| &bevels[*bi].surface)
        };
        let (Some(first), Some(second)) = (beside(a), beside(z)) else {
            ogeom_core::ogeom_bail!(
                Construction,
                "a mitred bevel's end shares no corner with another bevel"
            );
        };
        ring.push(section(points, edges, made, strip, Some(first), a, *x)?);
        ring.push(section(points, edges, made, strip, Some(second), *x, z)?);
        Ok(true)
    };
    // A face's side of the cap at a vertex, from `a` to `b`: the face cut
    // by the cap's plane, shared by the face and the cap.
    let cap_side = |points: &mut Points,
                    edges: &mut Vec<SheetCurve>,
                    made: &mut FastMap<Key, usize>,
                    surface: &SurfaceGeometry,
                    plane: &SurfaceGeometry,
                    a: Point,
                    b: Point|
     -> OgeomResult<(usize, bool)> {
        let (i, j) = (points.index(a), points.index(b));
        let key = Key::Section(i.min(j), i.max(j));
        if let Some(k) = made.get(&key) {
            return Ok((*k, edges[*k].ends.0 == i));
        }
        let met = match surface {
            SurfaceGeometry::Plane(_) => None,
            _ => met_between(surface, plane, a, b, tol)?,
        };
        let (curve, range, ends): (Curve, (f64, f64), (usize, usize)) = match met {
            Some((curve, range, a_low)) => (curve, range, if a_low { (i, j) } else { (j, i) }),
            None => {
                let line = LineCurve::segment(a, b, tol)?;
                let range = line.domain();
                (line.into(), range, (i, j))
            }
        };
        edges.push(SheetCurve { ends, curve, range });
        made.insert(key, edges.len() - 1);
        Ok((edges.len() - 1, ends.0 == i))
    };
    // Where a face's ring stopped short of its next piece, from `p` to
    // `q`, at the vertices `at`: one bevel's end flush on the face, or two
    // meeting where a face between them was consumed.
    let bridge = |points: &mut Points,
                  edges: &mut Vec<SheetCurve>,
                  made: &mut FastMap<Key, usize>,
                  face: &super::curved::CurvedFace,
                  at: &[Id],
                  p: Point,
                  q: Point|
     -> OgeomResult<Option<Vec<(usize, bool)>>> {
        let close = |x: Point, y: Point| x.distance(y) <= tol.confusion();
        let ends: Vec<&(Point, Point, usize)> = at
            .iter()
            .filter_map(|v| sections.get(v))
            .flatten()
            .collect();
        let joins = |a: Point, b: Point| {
            ends.iter()
                .copied()
                .find(|(x, y, _)| (close(*x, a) && close(*y, b)) || (close(*y, a) && close(*x, b)))
        };
        let mut hops: Vec<(Point, Point, usize)> = Vec::new();
        if let Some((_, _, bi)) = joins(p, q) {
            hops.push((p, q, *bi));
        } else {
            let found = ends.iter().find_map(|(x, y, bi)| {
                let middle = if close(*x, p) {
                    *y
                } else if close(*y, p) {
                    *x
                } else {
                    return None;
                };
                let (_, _, other) = joins(middle, q)?;
                (other != bi).then_some((middle, *bi, *other))
            });
            let Some((middle, first, second)) = found else {
                // A capped vertex: the face's side of the cap, from its
                // corner to the vertex.
                let capped = at.iter().find_map(|v| {
                    let plane = facet_planes.get(v)?;
                    let point = curved.points[v];
                    (close(point, p) || close(point, q)).then_some(plane)
                });
                let Some(plane) = capped else {
                    return Ok(None);
                };
                return Ok(Some(vec![cap_side(
                    points,
                    edges,
                    made,
                    &face.surface,
                    plane,
                    p,
                    q,
                )?]));
            };
            hops.push((p, middle, first));
            hops.push((middle, q, second));
        }
        let mut pieces = Vec::with_capacity(hops.len());
        for (a, b, bi) in hops {
            let strip = &bevels[bi];
            if !flush(face, &strip.section, a, b, tol)? {
                return Ok(None);
            }
            pieces.push(section(
                points,
                edges,
                made,
                strip,
                Some(&face.surface),
                a,
                b,
            )?);
        }
        Ok(Some(pieces))
    };
    // Each cusp's hole: four sides, each one curve the patch meets, walked
    // from the side face over the first lobe to the wall, along the wall,
    // back over the second lobe and along the side. A lobe's side is its
    // straight strip's end, cut by its meridian plane, running on into the
    // rim's end in that plane, one rational quadratic; the wall's is one arc
    // of the rim's ring, cut at the seam's foot; the side's a straight. Each
    // edge of a side is a range of its curve, so the patch takes it exactly;
    // they are the strips' ends there and the faces' pieces where they stop
    // short of the patch.
    let mut bounds: FastMap<Id, Vec<(Point, Point, usize)>> = FastMap::default();
    let mut holes: Vec<(usize, Vec<HoleSpec>)> = Vec::with_capacity(cusps.len());
    // Each S's foot's edge, with how far its fitted curve stands off its
    // face, measured finely: it bends faster near the vertex than the
    // edges' common sampling reads.
    let mut feet_widths: Vec<(usize, f64)> = Vec::new();
    for (ci, c) in cusps.iter().enumerate() {
        let [_, side, wall] = c.faces;
        let mut lobes = Vec::with_capacity(2);
        let mut ring_at = Point::ORIGIN;
        for lobe in &c.lobes {
            let rim = &bevels[lobe.rim];
            let SurfaceGeometry::Torus(torus) = &rim.surface else {
                return Ok(None);
            };
            let torus = torus.torus();
            let (o, z) = (torus.frame().origin(), torus.frame().z().vector());
            let Curve::Circle(ring) = &rim.contacts[usize::from(rim.faces[1] == wall)] else {
                return Ok(None);
            };
            ring_at = ring.circle().frame().origin();
            let line = &bevels[lobe.line];
            let Section::Arc { spine } = &line.section else {
                return Ok(None);
            };
            let centre = spine.point_at(param_on(spine, spine.domain(), lobe.side, tol)?, tol)?;
            let far = spine.point_at(param_on(spine, spine.domain(), lobe.top, tol)?, tol)?;
            let SurfaceGeometry::Plane(plane) = &lobe.plane else {
                return Ok(None);
            };
            let normal = plane.plane().normal().vector();
            let straight =
                skew_control(centre, far - centre, lobe.side, lobe.top, Some(normal), tol)?;
            let radial = {
                let d = lobe.top - o;
                (d - z * d.dot(z)).normalized(tol)?
            };
            let core = o + radial * torus.major_radius();
            let round = arc_control(core, lobe.top, lobe.wall, tol)?;
            let (curve, ends) = rational_quadratics(&[straight, round], tol)?;
            let (i, j, k) = (
                points.index(lobe.side),
                points.index(lobe.top),
                points.index(lobe.wall),
            );
            let low = piece_of(&mut edges, &mut made, &curve, (i, j), (0.0, ends[0]));
            let high = piece_of(&mut edges, &mut made, &curve, (j, k), (ends[0], ends[1]));
            lobes.push(HoleSpec {
                curve,
                range: (0.0, ends[1]),
                forward: true,
                pieces: vec![
                    (low, (0.0, ends[0]), line.surface.clone()),
                    (high, (ends[0], ends[1]), rim.surface.clone()),
                ],
            });
        }
        // Two lobes: one hole over the vertex, the lobes and, between them,
        // the rim's contact on the wall and the straight strip's on the
        // side, so the patch stands clear of the vertex, where the side and
        // the wall are tangent.
        if let [a, b] = c.lobes.as_slice() {
            let straight: Curve = LineCurve::segment(b.side, a.side, tol)?.into();
            let range = straight.domain();
            edges.push(SheetCurve {
                ends: (points.index(b.side), points.index(a.side)),
                curve: straight.clone(),
                range,
            });
            bounds
                .entry(c.v)
                .or_default()
                .push((b.side, a.side, edges.len() - 1));
            let along_side = HoleSpec {
                curve: straight,
                range,
                forward: true,
                pieces: vec![(edges.len() - 1, range, curved.faces[side].surface.clone())],
            };
            let (arc, _) = rational_quadratics(&[arc_control(ring_at, a.wall, b.wall, tol)?], tol)?;
            let mut stops = vec![(a.wall, 0.0)];
            if let Some(f) = c.foot {
                stops.push((f, param_on(&arc, (0.0, 1.0), f, tol)?));
            }
            stops.push((b.wall, 1.0));
            let mut pieces = Vec::with_capacity(2);
            for w in stops.windows(2) {
                let ((p, t0), (q, t1)) = (w[0], w[1]);
                edges.push(SheetCurve {
                    ends: (points.index(p), points.index(q)),
                    curve: arc.clone(),
                    range: (t0, t1),
                });
                bounds.entry(c.v).or_default().push((p, q, edges.len() - 1));
                pieces.push((
                    edges.len() - 1,
                    (t0, t1),
                    curved.faces[wall].surface.clone(),
                ));
            }
            let along_wall = HoleSpec {
                curve: arc,
                range: (0.0, 1.0),
                forward: true,
                pieces,
            };
            let [first, second] = <[HoleSpec; 2]>::try_from(lobes).map_err(|_| {
                ogeom_core::OgeomError::Construction("a cusp of other than two lobes".into())
            })?;
            holes.push((ci, vec![first, along_wall, second.back(), along_side]));
            continue;
        }
        // One lobe: its hole is the lobe, its foot on the wall down to the
        // vertex, and its foot on the side back to the lobe.
        let at = curved.points[&c.v];
        for (spec, (lobe, (wall_foot, side_foot))) in
            lobes.into_iter().zip(c.lobes.iter().zip(&c.feet))
        {
            let mut hole = vec![spec];
            for (curve, p, q, face) in [
                (wall_foot, lobe.wall, at, wall),
                (side_foot, at, lobe.side, side),
            ] {
                let range = curve.domain();
                let surface = &curved.faces[face].surface;
                let mut off: f64 = 0.0;
                for i in 0..=256 {
                    let t = (range.1 - range.0).mul_add(f64::from(i) / 256.0, range.0);
                    let at = curve.point_at(t, tol)?;
                    off = off.max(
                        ogeom_geom::project::project_on_surface(surface, at, 16, tol)?.distance,
                    );
                }
                edges.push(SheetCurve {
                    ends: (points.index(p), points.index(q)),
                    curve: curve.clone(),
                    range,
                });
                feet_widths.push((edges.len() - 1, 1.25 * off + tol.confusion()));
                bounds.entry(c.v).or_default().push((p, q, edges.len() - 1));
                hole.push(HoleSpec {
                    curve: curve.clone(),
                    range,
                    forward: true,
                    pieces: vec![(edges.len() - 1, range, curved.faces[face].surface.clone())],
                });
            }
            holes.push((ci, hole));
        }
    }
    // The pieces of a face's boundary a cusp's patch meets, joining `p` to
    // `q` at `v`: one piece, or two meeting at a point between (a wall's
    // ring arc and the tangent line down to the side's contact), each its
    // edge and whether it runs from the end nearer `p`.
    let bound_between = |v: Id, p: Point, q: Point| -> Option<Vec<(usize, bool)>> {
        let near = |a: Point, b: Point| a.distance(b) <= tol.confusion();
        let b = bounds.get(&v)?;
        // The pieces from `p` to `q`, each taken once: searched depth first,
        // the shortest found first since a direct piece is tried first.
        let mut stack: Vec<(Point, Vec<usize>)> = vec![(p, Vec::new())];
        let mut best: Option<Vec<usize>> = None;
        while let Some((at, path)) = stack.pop() {
            if near(at, q) && !path.is_empty() {
                if best.as_ref().is_none_or(|b| path.len() < b.len()) {
                    best = Some(path);
                }
                continue;
            }
            for (i, (x, y, _)) in b.iter().enumerate() {
                if path.contains(&i) {
                    continue;
                }
                let next = if near(*x, at) {
                    *y
                } else if near(*y, at) {
                    *x
                } else {
                    continue;
                };
                let mut longer = path.clone();
                longer.push(i);
                stack.push((next, longer));
            }
        }
        let path = best?;
        let mut at = p;
        let mut walk = Vec::with_capacity(path.len());
        for i in path {
            let (x, y, k) = b[i];
            let from_at = near(x, at);
            walk.push((k, from_at));
            at = if from_at { y } else { x };
        }
        Some(walk)
    };
    // The faces a rolled round cut a piece from, their rings rejoined, and
    // the faces it split off them, each with the face it split from.
    let mut relink: FastSet<usize> = FastSet::default();
    let mut split_faces: Vec<(usize, SheetRings)> = Vec::new();
    // P4: each face rebuilt, its loops in the stored order, each edge's
    // piece between the corners at its ends, walked the surface's way.
    let mut rings_of: Vec<SheetRings> = Vec::with_capacity(curved.faces.len() + bevels.len());
    for (fi, face) in curved.faces.iter().enumerate() {
        if consumed.contains(&fi) {
            // Not rebuilt; a placeholder keeps the faces' indices.
            rings_of.push(SheetRings {
                surface: face.surface.clone(),
                rings: Vec::new(),
                reversed: face.reversed,
            });
            continue;
        }
        let mut rings = Vec::with_capacity(face.loops.len());
        for walk in &face.loops {
            let mut ring: Vec<(usize, bool)> = Vec::new();
            let mut last: Option<Point> = None;
            let mut first: Option<Point> = None;
            // The vertices of the edges skipped since the last piece, and
            // before the first, which the loop's closing gap takes.
            let mut skipped: Vec<Id> = Vec::new();
            let mut leading: Vec<Id> = Vec::new();
            for (e, w) in walk {
                let forward = *w != face.reversed;
                let edge = &curved.edges[e];
                // An edge left sharp against a consumed face goes with it;
                // the strips' ends there close the gap.
                if !blended.contains(e) && faces_of(curved, *e).iter().any(|f| consumed.contains(f))
                {
                    skipped.extend([edge.ends.0, edge.ends.1]);
                    continue;
                }
                let (vs, ve) = if forward {
                    (edge.ends.0, edge.ends.1)
                } else {
                    (edge.ends.1, edge.ends.0)
                };
                // This edge's piece on this face: its contact, or itself.
                let (curve, start, end, key) = if blended.contains(e) {
                    let (Some(c), Some(a), Some(b)) = (
                        contact_of(fi, *e),
                        corner_for(corners, curved, fi, vs, *e, forward, tol)?,
                        corner_for(corners, curved, fi, ve, *e, !forward, tol)?,
                    ) else {
                        return Ok(None);
                    };
                    let (Some(ta), Some(tb)) = (a.on(*e), b.on(*e)) else {
                        return Ok(None);
                    };
                    (
                        c.clone(),
                        (a.point, ta),
                        (b.point, tb),
                        Key::Contact(fi, *e),
                    )
                } else {
                    let (pa, pb) = (sharp_end[&(*e, vs)], sharp_end[&(*e, ve)]);
                    let ta = param_on(&edge.curve, edge.range, pa, tol)?;
                    let tb = param_on(&edge.curve, edge.range, pb, tol)?;
                    (edge.curve.clone(), (pa, ta), (pb, tb), Key::Sharp(*e))
                };
                // Where the last piece stopped short of this one, a
                // bevel's end at the vertex runs between them.
                if let Some(p) = last
                    && p.distance(start.0) > tol.confusion()
                    && let Some(piece) =
                        tail_between(&mut points, &mut edges, &mut made, vs, p, start.0)?
                {
                    ring.push(piece);
                    last = Some(start.0);
                }
                // A piece of the face's boundary a cusp's patch meets.
                if let Some(p) = last
                    && p.distance(start.0) > tol.confusion()
                    && let Some(pieces) = bound_between(vs, p, start.0)
                {
                    ring.extend(pieces);
                    last = Some(start.0);
                }
                if let Some(p) = last
                    && p.distance(start.0) > tol.confusion()
                {
                    let mut at = vec![vs];
                    at.append(&mut skipped);
                    let Some(pieces) =
                        bridge(&mut points, &mut edges, &mut made, face, &at, p, start.0)?
                    else {
                        return Ok(None);
                    };
                    ring.extend(pieces);
                }
                if last.is_none() {
                    leading.append(&mut skipped);
                }
                skipped.clear();
                let along_edge = {
                    let t = param_on(&edge.curve, edge.range, start.0, tol)?;
                    edge.curve.d1_at(t, tol)? * if forward { 1.0 } else { -1.0 }
                };
                let rising = runs(&curve, start.1, along_edge, tol)?;
                let (range, against) = piece(&curve, start.1, end.1, rising, tol);
                // A piece a rolled round cuts: its parts in the walk's order,
                // each kept, dropped (the ring rejoined below) or taken by
                // the run's contact.
                if let Some(cuts) = spans.get(&(fi, *e)) {
                    let period = period_of(&curve, tol);
                    let (lo_point, hi_point) = if against {
                        (end.0, start.0)
                    } else {
                        (start.0, end.0)
                    };
                    let mut stops: Vec<(f64, Point)> =
                        vec![(range.0, lo_point), (range.1, hi_point)];
                    for c in cuts {
                        for (t, p) in [c.lo, c.hi] {
                            let t = wrapped(t, range.0, period);
                            if t > range.0 + tol.parametric() && t < range.1 - tol.parametric() {
                                stops.push((t, p));
                            }
                        }
                    }
                    stops.sort_by(|a, b| a.0.total_cmp(&b.0));
                    stops.dedup_by(|a, b| (a.0 - b.0).abs() <= tol.parametric());
                    let mut parts: Vec<Part<'_>> = Vec::with_capacity(stops.len() - 1);
                    for w in stops.windows(2) {
                        let mid = 0.5 * (w[0].0 + w[1].0);
                        let fate = cuts
                            .iter()
                            .find(|c| {
                                let t = wrapped(mid, c.lo.0, period);
                                t >= c.lo.0 && t < c.hi.0
                            })
                            .map(|c| &c.fate);
                        parts.push((w[0], w[1], fate));
                    }
                    if against {
                        parts.reverse();
                    }
                    let owner = if blended.contains(e) { fi } else { usize::MAX };
                    for (a, z, fate) in parts {
                        let (ia, iz) = (points.index(a.1), points.index(z.1));
                        let key = Key::Piece(owner, *e, ia.min(iz), ia.max(iz));
                        let (k, forward) = match fate {
                            Some(Fate::Drop) => {
                                relink.insert(fi);
                                continue;
                            }
                            Some(Fate::Replace(c)) => {
                                let k = if let Some(k) = made.get(&key) {
                                    *k
                                } else {
                                    let (lo, hi) = c.domain();
                                    let starts_at_a = c.point_at(lo, tol)?.distance(a.1)
                                        <= c.point_at(lo, tol)?.distance(z.1);
                                    edges.push(SheetCurve {
                                        ends: if starts_at_a { (ia, iz) } else { (iz, ia) },
                                        curve: c.clone(),
                                        range: (lo, hi),
                                    });
                                    made.insert(key, edges.len() - 1);
                                    edges.len() - 1
                                };
                                // Walked from `a` to `z`, or back.
                                let from = if against { iz } else { ia };
                                (k, edges[k].ends.0 == from)
                            }
                            Some(Fate::Keep) | None => {
                                let k = if let Some(k) = made.get(&key) {
                                    *k
                                } else {
                                    edges.push(SheetCurve {
                                        ends: (ia, iz),
                                        curve: curve.clone(),
                                        range: (a.0, z.0),
                                    });
                                    made.insert(key, edges.len() - 1);
                                    edges.len() - 1
                                };
                                (k, !against)
                            }
                        };
                        ring.push((k, forward));
                    }
                    if first.is_none() {
                        first = Some(start.0);
                    }
                    last = Some(end.0);
                    continue;
                }
                let (i, j) = (points.index(start.0), points.index(end.0));
                let k = if let Some(k) = made.get(&key) {
                    *k
                } else {
                    let (lo_end, hi_end) = if against { (j, i) } else { (i, j) };
                    edges.push(SheetCurve {
                        ends: (lo_end, hi_end),
                        curve: curve.clone(),
                        range,
                    });
                    made.insert(key, edges.len() - 1);
                    edges.len() - 1
                };
                ring.push((k, !against));
                if first.is_none() {
                    first = Some(start.0);
                }
                last = Some(end.0);
            }
            // The loop closes: a bevel's end where it stopped short.
            if let (Some(p), Some(q)) = (last, first)
                && p.distance(q) > tol.confusion()
                && let Some(piece) = tail_between(
                    &mut points,
                    &mut edges,
                    &mut made,
                    loop_start_vertex(curved, face, walk),
                    p,
                    q,
                )?
            {
                ring.push(piece);
                last = Some(q);
            }
            if let (Some(p), Some(q)) = (last, first)
                && p.distance(q) > tol.confusion()
                && let Some(pieces) = bound_between(loop_start_vertex(curved, face, walk), p, q)
            {
                ring.extend(pieces);
                last = Some(q);
            }
            if let (Some(p), Some(q)) = (last, first)
                && p.distance(q) > tol.confusion()
            {
                let mut at = vec![loop_start_vertex(curved, face, walk)];
                at.append(&mut skipped);
                at.append(&mut leading);
                let Some(pieces) = bridge(&mut points, &mut edges, &mut made, face, &at, p, q)?
                else {
                    return Ok(None);
                };
                ring.extend(pieces);
            }
            rings.push(ring);
        }
        // A face a rolled round cut a piece from: its pieces rejoined end to
        // end into rings, and each ring walking the surface's way round with
        // the rings inside it a face of its own.
        let rings = if relink.contains(&fi) {
            let Some(rejoined) = rejoined(&edges, &rings) else {
                return Ok(None);
            };
            let Some(mut groups) = grouped(&face.surface, &edges, rejoined, tol)? else {
                return Ok(None);
            };
            let first = groups.remove(0);
            for group in groups {
                split_faces.push((
                    fi,
                    SheetRings {
                        surface: face.surface.clone(),
                        rings: group,
                        reversed: face.reversed,
                    },
                ));
            }
            first
        } else {
            rings
        };
        rings_of.push(SheetRings {
            surface: face.surface.clone(),
            rings,
            reversed: face.reversed,
        });
    }
    // Each bevel: its contact on its first face, walked against the way
    // that face walks it (faces beside each other walk a shared edge
    // opposite ways, each as it faces out), its end, its contact on its
    // second face the other way round, and its other end.
    // The middles of the contacts made on consumed faces.
    let middles_made = core::cell::RefCell::new(Points {
        points: Vec::new(),
        merge: tol.confusion(),
    });
    // The faces each rolled strip's pieces after its first take, appended
    // after every other face, with the strip's index.
    let mut rolled_faces: Vec<(usize, SheetRings)> = Vec::new();
    for b in bevels {
        let e = &curved.edges[&b.edge];
        // A rolled strip: its pieces either side of its run, each between
        // its faces, and the run's piece, the ball on the edge.
        if let Some(r) = rolled.iter().find(|r| bevels[r.strip].edge == b.edge) {
            if rolled.iter().filter(|q| q.strip == r.strip).count() != 1 {
                return Ok(None);
            }
            let (left, other) = (b.faces[r.side], b.faces[1 - r.side]);
            let (Some(&q0), Some(&q1)) = (r.run.feet.first(), r.run.feet.last()) else {
                return Ok(None);
            };
            let walked = |fi: usize, k: usize| -> Option<bool> {
                face_walks(&rings_of[fi], k).or_else(|| {
                    split_faces
                        .iter()
                        .filter(|(f, _)| *f == fi)
                        .find_map(|(_, rings)| face_walks(rings, k))
                })
            };
            let piece_key = |owner: usize, edge: Id, a: Point, z: Point, points: &mut Points| {
                let (i, j) = (points.index(a), points.index(z));
                Key::Piece(owner, edge, i.min(j), i.max(j))
            };
            let mut rings: Vec<Vec<(usize, bool)>> = Vec::with_capacity(3);
            // The pieces between the faces, at each end of the edge.
            for (end, v) in [e.ends.0, e.ends.1].into_iter().enumerate() {
                let (Some(cl), Some(co)) = (
                    corner_for(corners, curved, left, v, b.edge, end == 0, tol)?,
                    corner_for(corners, curved, other, v, b.edge, end == 0, tol)?,
                ) else {
                    return Ok(None);
                };
                // The crossing on this end's side: the one nearer its corner
                // along the contact.
                let c = &b.contacts[r.side];
                let tc = param_on(c, c.domain(), cl.point, tol)?;
                let period = period_of(c, tol);
                let gap = |t: f64| {
                    let d = t - tc;
                    period.map_or(d.abs(), |p| {
                        let d = d.rem_euclid(p);
                        d.min(p - d)
                    })
                };
                let (p, q) = if gap(r.enter.at) <= gap(r.leave.at) {
                    (r.enter.point, q0)
                } else {
                    (r.leave.point, q1)
                };
                let (a, z) = (cl.point, co.point);
                let key_l = piece_key(left, b.edge, a, p, &mut points);
                let key_o = piece_key(other, b.edge, z, q, &mut points);
                let (Some(&kl), Some(&ko)) = (made.get(&key_l), made.get(&key_o)) else {
                    return Ok(None);
                };
                let switch = section(&mut points, &mut edges, &mut made, b, None, p, q)?;
                let cut = cut_at(b, v, a, z);
                let capped = section(&mut points, &mut edges, &mut made, b, cut.as_ref(), z, a)?;
                let Some(wl) = walked(left, kl) else {
                    return Ok(None);
                };
                let Some(ring) = chained_ring(&edges, (kl, !wl), &[ko, switch.0, capped.0]) else {
                    return Ok(None);
                };
                rings.push(ring);
            }
            // The run's piece: the other face's contact along it, walked
            // against that face's way, the switches and the edges between
            // the crossings.
            let key_o = piece_key(other, b.edge, q0, q1, &mut points);
            let Some(&ko) = made.get(&key_o) else {
                return Ok(None);
            };
            let mut stops = vec![r.enter.point];
            for w in r.run.over.windows(2) {
                let (a, z) = (&curved.edges[&w[0]], &curved.edges[&w[1]]);
                let v = if a.ends.0 == z.ends.0 || a.ends.0 == z.ends.1 {
                    a.ends.0
                } else {
                    a.ends.1
                };
                stops.push(curved.points[&v]);
            }
            stops.push(r.leave.point);
            let mut others: Vec<usize> = Vec::with_capacity(stops.len() + 1);
            for (k, w) in stops.windows(2).enumerate() {
                let key = piece_key(usize::MAX, r.run.over[k], w[0], w[1], &mut points);
                let Some(&kk) = made.get(&key) else {
                    return Ok(None);
                };
                others.push(kk);
            }
            for (p, q) in [(r.enter.point, q0), (r.leave.point, q1)] {
                others.push(section(&mut points, &mut edges, &mut made, b, None, p, q)?.0);
            }
            let Some(wo) = walked(other, ko) else {
                return Ok(None);
            };
            let Some(run_ring) = chained_ring(&edges, (ko, !wo), &others) else {
                return Ok(None);
            };
            // Each faces out, away from both faces there; a ring runs the
            // surface's way.
            let oriented = |surface: &SurfaceGeometry,
                            ring: Vec<(usize, bool)>,
                            at: Point|
             -> OgeomResult<SheetRings> {
                // Read where the piece stands: round an arc the faces'
                // normals turn.
                let out = outward_at(curved, left, at) + outward_at(curved, other, at);
                let reversed = natural_normal(surface, at, tol)?.dot(out) < 0.0;
                let ring = if reversed {
                    ring.iter().rev().map(|(k, f)| (*k, !f)).collect()
                } else {
                    ring
                };
                Ok(SheetRings {
                    surface: surface.clone(),
                    rings: vec![ring],
                    reversed,
                })
            };
            let mut made_rings = Vec::with_capacity(3);
            for ring in rings {
                let at = points.points[edges[ring[0].0].ends.0];
                made_rings.push(oriented(&b.surface, ring, at)?);
            }
            let middle = r.run.centres.len() / 2;
            let run_at = r.run.feet[middle].midpoint(r.run.on_edge[middle]);
            let run_at =
                ogeom_geom::project::project_on_surface(&r.run.surface, run_at, 16, tol)?.point;
            made_rings.push(oriented(&r.run.surface, run_ring, run_at)?);
            let mut made_rings = made_rings.into_iter();
            let Some(first) = made_rings.next() else {
                return Ok(None);
            };
            rings_of.push(first);
            let strip = rolled.iter().position(|q| q.strip == r.strip).unwrap_or(0);
            rolled_faces.extend(made_rings.map(|m| (rolled[strip].strip, m)));
            continue;
        }
        // A contact's piece: the face's own, or, on a consumed face, made
        // here between the corners at the edge's ends, shared by every
        // strip meeting along it there.
        let piece_of = |k: usize,
                        points: &mut Points,
                        edges: &mut Vec<SheetCurve>,
                        made: &mut FastMap<Key, usize>|
         -> OgeomResult<Option<usize>> {
            let fi = b.faces[k];
            if let Some(found) = made.get(&Key::Contact(fi, b.edge)) {
                return Ok(Some(*found));
            }
            let (Some(a), Some(z)) = (corners.get(&(fi, e.ends.0)), corners.get(&(fi, e.ends.1)))
            else {
                return Ok(None);
            };
            let contact = &b.contacts[k];
            let (ta, tz) = (
                param_on(contact, contact.domain(), a.point, tol)?,
                param_on(contact, contact.domain(), z.point, tol)?,
            );
            let (i, j) = (points.index(a.point), points.index(z.point));
            // Round a closed contact, the way that stands beside the edge:
            // its middle nearer the edge's own.
            let beside = {
                let middle = e.curve.point_at(0.5 * (e.range.0 + e.range.1), tol)?;
                let off = |(lo, hi): (f64, f64)| -> OgeomResult<f64> {
                    Ok(contact.point_at(0.5 * (lo + hi), tol)?.distance(middle))
                };
                let (up, down) = (
                    piece(contact, ta, tz, true, tol),
                    piece(contact, ta, tz, false, tol),
                );
                if !(contact.is_periodic() || contact.is_closed(tol)) {
                    // Along an open one, the way its ends run.
                    piece(contact, ta, tz, ta <= tz, tol)
                } else if off(down.0)? < off(up.0)? {
                    down
                } else {
                    up
                }
            };
            let (range, against) = beside;
            let m = middles_made
                .borrow_mut()
                .index(contact.point_at(0.5 * (range.0 + range.1), tol)?);
            let key = Key::Shared(fi, i.min(j), i.max(j), m);
            if let Some(found) = made.get(&key).copied() {
                made.insert(Key::Contact(fi, b.edge), found);
                return Ok(Some(found));
            }
            let ends = if against { (j, i) } else { (i, j) };
            edges.push(SheetCurve {
                ends,
                curve: contact.clone(),
                range,
            });
            made.insert(key, edges.len() - 1);
            made.insert(Key::Contact(fi, b.edge), edges.len() - 1);
            Ok(Some(edges.len() - 1))
        };
        let (Some(k0), Some(k1)) = (
            piece_of(0, &mut points, &mut edges, &mut made)?,
            piece_of(1, &mut points, &mut edges, &mut made)?,
        ) else {
            return Ok(None);
        };
        // Walk the contact on a face rebuilt against the way that face walks
        // it, as it faces out; then on to the other contact.
        let r = usize::from(consumed.contains(&b.faces[0]));
        if consumed.contains(&b.faces[r]) {
            return Ok(None);
        }
        let (kr, ko) = if r == 0 { (k0, k1) } else { (k1, k0) };
        let walked = face_walks(&rings_of[b.faces[r]], kr).or_else(|| {
            split_faces
                .iter()
                .filter(|(f, _)| *f == b.faces[r])
                .find_map(|(_, rings)| face_walks(rings, kr))
        });
        let Some(walked) = walked else {
            return Ok(None);
        };
        let sr = !(walked != rings_of[b.faces[r]].reversed);
        let ends = |edges: &Vec<SheetCurve>, k: usize, forward: bool| {
            let (a, z) = edges[k].ends;
            if forward { (a, z) } else { (z, a) }
        };
        let (pr, qr) = ends(&edges, kr, sr);
        // The vertex an end of the strip stands at: the one whose corners
        // it joins.
        let merge = points.merge;
        let vertex_of = |x: Point| -> Id {
            let near = |v: Id| {
                [b.faces[0], b.faces[1]].iter().any(|f| {
                    [v, pinched(v)].iter().any(|w| {
                        corners
                            .get(&(*f, *w))
                            .is_some_and(|c| c.point.distance(x) <= merge)
                    })
                })
            };
            if near(e.ends.0) { e.ends.0 } else { e.ends.1 }
        };
        let other = b.faces[1 - r];
        // The other contact: from its corner where the first ended; a
        // closed one the other way round the strip's axis.
        // The vertex's corner on the other face, or its second where the
        // face passes the vertex twice: the one ending the other contact.
        let v = vertex_of(points.points[qr]);
        let start_o = [v, pinched(v)]
            .iter()
            .filter_map(|w| corners.get(&(other, *w)))
            .map(|c| points.index(c.point))
            .find(|i| *i == edges[ko].ends.0 || *i == edges[ko].ends.1);
        let Some(start_o) = start_o else {
            return Ok(None);
        };
        let so = if edges[ko].ends.0 != edges[ko].ends.1 {
            edges[ko].ends.0 == start_o
        } else {
            match (
                turning(&b.surface, &edges[kr], tol)?,
                turning(&b.surface, &edges[ko], tol)?,
            ) {
                (Some(turn_r), Some(turn_o)) => {
                    (turn_r * if sr { 1.0 } else { -1.0 }) * turn_o < 0.0
                }
                // A closed ruled bevel's contacts come from the same stations
                // made compatible, so they run one way along it: the other
                // is walked against the first.
                _ if matches!(b.surface, SurfaceGeometry::BSpline(_)) => !sr,
                _ => return Ok(None),
            }
        };
        let (po, qo) = ends(&edges, ko, so);
        let mut ring: Vec<(usize, bool)> = vec![(kr, sr)];
        if qr != po {
            let (a, z) = (points.points[qr], points.points[po]);
            if !mitred_end(&mut points, &mut edges, &mut made, &mut ring, b, a, z)? {
                let cut = cut_at(b, vertex_of(a), a, z);
                ring.push(section(
                    &mut points,
                    &mut edges,
                    &mut made,
                    b,
                    cut.as_ref(),
                    a,
                    z,
                )?);
            }
        }
        if consumed.contains(&other) {
            // On a consumed face the strips meet along the ridge, each
            // other strip's corner on it splitting it: walk its pieces.
            let curve = edges[ko].curve.clone();
            // The piece's range runs from its first end to its second,
            // either way along an open curve.
            let (r0, r1) = edges[ko].range;
            let (lo, hi) = (r0.min(r1), r0.max(r1));
            let mut cuts: Vec<(f64, usize)> = Vec::new();
            for ((fi, _), c) in corners.iter() {
                if *fi != other {
                    continue;
                }
                let t = param_on(&curve, (lo, hi), c.point, tol)?;
                let i = points.index(c.point);
                if t > lo + tol.parametric()
                    && t < hi - tol.parametric()
                    && curve.point_at(t, tol)?.distance(c.point) <= tol.confusion()
                    && i != edges[ko].ends.0
                    && i != edges[ko].ends.1
                    && !cuts.iter().any(|(_, j)| *j == i)
                {
                    cuts.push((t, i));
                }
            }
            if cuts.is_empty() {
                ring.push((ko, so));
            } else {
                cuts.sort_by(|a, b| a.0.total_cmp(&b.0));
                if r0 > r1 {
                    cuts.reverse();
                }
                let mut stops = vec![(r0, edges[ko].ends.0)];
                stops.extend(cuts);
                stops.push((r1, edges[ko].ends.1));
                let mut pieces: Vec<(usize, bool)> = Vec::with_capacity(stops.len() - 1);
                for w in stops.windows(2) {
                    let ((ta, ia), (tb, ib)) = (w[0], w[1]);
                    let m = middles_made
                        .borrow_mut()
                        .index(curve.point_at(0.5 * (ta + tb), tol)?);
                    let key = Key::Shared(other, ia.min(ib), ia.max(ib), m);
                    let k = if let Some(k) = made.get(&key) {
                        *k
                    } else {
                        edges.push(SheetCurve {
                            ends: (ia, ib),
                            curve: curve.clone(),
                            range: (ta, tb),
                        });
                        made.insert(key, edges.len() - 1);
                        if let Some(w) = widths.get(&(other, b.edge)) {
                            marched_widths.borrow_mut().insert(edges.len() - 1, *w);
                        }
                        edges.len() - 1
                    };
                    pieces.push((k, edges[k].ends.0 == ia));
                }
                if !so {
                    pieces.reverse();
                    for piece in &mut pieces {
                        piece.1 = !piece.1;
                    }
                }
                ring.extend(pieces);
            }
        } else {
            ring.push((ko, so));
        }
        if qo != pr {
            let (a, z) = (points.points[qo], points.points[pr]);
            if !mitred_end(&mut points, &mut edges, &mut made, &mut ring, b, a, z)? {
                let cut = cut_at(b, vertex_of(a), a, z);
                ring.push(section(
                    &mut points,
                    &mut edges,
                    &mut made,
                    b,
                    cut.as_ref(),
                    a,
                    z,
                )?);
            }
        }
        // The bevel faces out: away from both faces at the edge, read at
        // its middle (an arc's chord's middle stands off it, at a half
        // circle's centre).
        let mid = e.curve.point_at(0.5 * (e.range.0 + e.range.1), tol)?;
        let out = outward_at(curved, b.faces[0], mid) + outward_at(curved, b.faces[1], mid);
        let natural = natural_normal(&b.surface, mid, tol)?;
        let reversed = natural.dot(out) < 0.0;
        // The ring above runs as the bevel faces out; the sheet wants it
        // the surface's way.
        let ring = if reversed {
            ring.iter().rev().map(|(k, f)| (*k, !f)).collect()
        } else {
            ring
        };
        rings_of.push(SheetRings {
            surface: b.surface.clone(),
            rings: vec![ring],
            reversed,
        });
    }
    // Two strips on one surface running on into each other (pieces of one
    // edge the boolean split) are one face: their rings join across the end
    // they share, which goes, but at a fitted band's seam. `joined[i]` is
    // the strip ring `i` was joined into.
    let first_strip = curved.faces.len();
    let mut joined: Vec<usize> = (0..bevels.len()).collect();
    let find = |joined: &Vec<usize>, mut i: usize| {
        while joined[i] != i {
            i = joined[i];
        }
        i
    };
    let mut keyed_v: Vec<&Id> = curved.points.keys().collect();
    keyed_v.sort_unstable();
    for v in keyed_v {
        let Some(ends) = sections.get(v) else {
            continue;
        };
        let [(p, q, a), (_, _, b)] = ends[..] else {
            continue;
        };
        if bevels[a].surface != bevels[b].surface {
            continue;
        }
        // Not where a closed band's contacts start, the seam of its chart,
        // on a surface that does not repeat (a fitted band); one that
        // repeats (a cone) takes its seam where its faces are joined.
        let repeats = ogeom_geom::Surface::is_periodic_u(&bevels[a].surface)
            || ogeom_geom::Surface::is_periodic_v(&bevels[a].surface);
        let seam = !repeats
            && bevels[a].contacts.first().is_some_and(|c| {
                let (lo, hi) = c.domain();
                lo.is_finite()
                    && hi.is_finite()
                    && c.point_at(lo, tol).is_ok_and(|s| {
                        c.point_at(hi, tol)
                            .is_ok_and(|e| s.distance(e) <= points.merge)
                            && (s.distance(p) <= points.merge || s.distance(q) <= points.merge)
                    })
            });
        if seam {
            continue;
        }
        let (i, j) = (points.index(p), points.index(q));
        let Some(&k) = made.get(&Key::Section(i.min(j), i.max(j))) else {
            continue;
        };
        let (ra, rb) = (find(&joined, a), find(&joined, b));
        if ra == rb || rings_of[first_strip + ra].reversed != rings_of[first_strip + rb].reversed {
            continue;
        }
        // The one starting lower on its curve keeps the joined ring.
        let start = |r: usize| curved.edges[&bevels[r].edge].range.0;
        let (ra, rb) = if start(rb) < start(ra) {
            (rb, ra)
        } else {
            (ra, rb)
        };
        let ring_a = rings_of[first_strip + ra].rings[0].clone();
        let ring_b = rings_of[first_strip + rb].rings[0].clone();
        let (Some(ia), Some(ib)) = (
            ring_a.iter().position(|(x, _)| *x == k),
            ring_b.iter().position(|(x, _)| *x == k),
        ) else {
            continue;
        };
        if ring_a[ia].1 == ring_b[ib].1 {
            continue;
        }
        let mut ring: Vec<(usize, bool)> = ring_a[ia + 1..].to_vec();
        ring.extend_from_slice(&ring_a[..ia]);
        ring.extend_from_slice(&ring_b[ib + 1..]);
        ring.extend_from_slice(&ring_b[..ib]);
        rings_of[first_strip + ra].rings = vec![ring];
        joined[rb] = ra;
    }
    // Each facet: its ring the bevels' ends in turn about the side it
    // faces.
    for (v, cycle) in &facet_rings {
        let mut ring = Vec::with_capacity(cycle.len());
        for k in 0..cycle.len() {
            let (p, q) = (cycle[k], cycle[(k + 1) % cycle.len()]);
            let strip = sections.get(v).and_then(|ends| {
                ends.iter()
                    .find(|(x, y, _)| {
                        (x.distance(p) <= tol.confusion() && y.distance(q) <= tol.confusion())
                            || (y.distance(p) <= tol.confusion()
                                && x.distance(q) <= tol.confusion())
                    })
                    .map(|(_, _, bi)| &bevels[*bi])
            });
            let Some(strip) = strip else {
                if let Some(piece) = tail_between(&mut points, &mut edges, &mut made, *v, p, q)? {
                    ring.push(piece);
                    continue;
                }
                // A cap's side on a face: from the face's corner at the
                // vertex to the vertex itself.
                let at = curved.points[v];
                let corner = if p.distance(at) <= tol.confusion() {
                    q
                } else {
                    p
                };
                // A face passing the vertex twice (a wall's seam leaving it)
                // keeps its second corner there under the pinched key.
                let face = (0..curved.faces.len()).find(|fi| {
                    [*v, pinched(*v)].iter().any(|w| {
                        corners
                            .get(&(*fi, *w))
                            .is_some_and(|c| c.point.distance(corner) <= points.merge)
                    })
                });
                let (Some(fi), Some(plane)) = (face, facet_planes.get(v)) else {
                    return Ok(None);
                };
                let surface = curved.faces[fi].surface.clone();
                ring.push(cap_side(
                    &mut points,
                    &mut edges,
                    &mut made,
                    &surface,
                    plane,
                    p,
                    q,
                )?);
                continue;
            };
            ring.push(section(
                &mut points,
                &mut edges,
                &mut made,
                strip,
                facet_planes.get(v),
                p,
                q,
            )?);
        }
        let Some(surface) = facet_planes.get(v).or_else(|| patches.get(v)) else {
            return Ok(None);
        };
        rings_of.push(SheetRings {
            surface: surface.clone(),
            rings: vec![ring],
            reversed: false,
        });
    }
    // Each ball: its ring the rounds' end arcs in turn, wound about its
    // surface's own normal (off the centre), the face turned over where the
    // corner is concave. Convex corners only: there the sphere stands inside
    // the material all the faces round the vertex bound.
    let mut ball_regions: Vec<ogeom_algo::HalfSpaces> = Vec::with_capacity(balls.len());
    for (v, cycle, centre, radius) in &balls {
        let mut out = Vector::ZERO;
        let mut regions = Vec::new();
        let mut round: Vec<usize> = Vec::new();
        for (e, c) in &curved.edges {
            if c.ends.0 == *v || c.ends.1 == *v {
                for fi in faces_of(curved, *e) {
                    out += outward_at(curved, fi, curved.points[v]);
                    if !round.contains(&fi) {
                        round.push(fi);
                    }
                }
            }
        }
        if (*centre - curved.points[v]).dot(out) >= 0.0 {
            return Ok(None);
        }
        round.sort_unstable();
        for fi in round {
            let Some(m) = material(curved, fi) else {
                return Ok(None);
            };
            regions.push(m);
        }
        ball_regions.push(ogeom_algo::HalfSpaces::Every(regions));
        let mut cycle = cycle.clone();
        let mut newell = Vector::ZERO;
        let mut mean = Vector::ZERO;
        for k in 0..cycle.len() {
            let (p, q) = (
                cycle[k].to_vector(),
                cycle[(k + 1) % cycle.len()].to_vector(),
            );
            newell += p.cross(q);
            mean += p;
        }
        #[expect(clippy::cast_precision_loss, reason = "a vertex's few rounds")]
        let mean = Point::from_vector(mean / cycle.len() as f64);
        if newell.dot(mean - *centre) < 0.0 {
            cycle.reverse();
        }
        let surface = crate::round::sphere::corner_sphere(*centre, *radius, &cycle, tol)?;
        let mut ring = Vec::with_capacity(cycle.len());
        for k in 0..cycle.len() {
            let (p, q) = (cycle[k], cycle[(k + 1) % cycle.len()]);
            let strip = sections.get(v).and_then(|ends| {
                ends.iter()
                    .find(|(x, y, _)| {
                        (x.distance(p) <= tol.confusion() && y.distance(q) <= tol.confusion())
                            || (y.distance(p) <= tol.confusion()
                                && x.distance(q) <= tol.confusion())
                    })
                    .map(|(_, _, bi)| &bevels[*bi])
            });
            let Some(strip) = strip else {
                return Ok(None);
            };
            ring.push(section(
                &mut points,
                &mut edges,
                &mut made,
                strip,
                None,
                p,
                q,
            )?);
        }
        rings_of.push(SheetRings {
            surface,
            rings: vec![ring],
            reversed: false,
        });
    }
    // Each patch: the rounds' end arcs in turn, wound counter-clockwise seen
    // from outside, each with the strip it ends; its region the material of
    // every face round the vertex. Convex corners only: every end arc's
    // centre inside the material.
    let mut patch_rings: Vec<PatchRing> = Vec::with_capacity(vertex_patches.len());
    for (v, cycle, centres) in &vertex_patches {
        let mut out = Vector::ZERO;
        let mut round: Vec<usize> = Vec::new();
        for (e, c) in &curved.edges {
            if c.ends.0 == *v || c.ends.1 == *v {
                for fi in faces_of(curved, *e) {
                    out += outward_at(curved, fi, curved.points[v]);
                    if !round.contains(&fi) {
                        round.push(fi);
                    }
                }
            }
        }
        if centres
            .iter()
            .any(|c| (*c - curved.points[v]).dot(out) >= 0.0)
        {
            return Ok(None);
        }
        round.sort_unstable();
        let mut regions = Vec::with_capacity(round.len());
        for fi in round {
            let Some(m) = material(curved, fi) else {
                return Ok(None);
            };
            regions.push(m);
        }
        let mut cycle = cycle.clone();
        let mut newell = Vector::ZERO;
        for k in 0..cycle.len() {
            newell += cycle[k]
                .to_vector()
                .cross(cycle[(k + 1) % cycle.len()].to_vector());
        }
        if newell.dot(out) < 0.0 {
            cycle.reverse();
        }
        let mut ring = Vec::with_capacity(cycle.len());
        for k in 0..cycle.len() {
            let (p, q) = (cycle[k], cycle[(k + 1) % cycle.len()]);
            let found = sections.get(v).and_then(|ends| {
                ends.iter().find(|(x, y, _)| {
                    (x.distance(p) <= tol.confusion() && y.distance(q) <= tol.confusion())
                        || (y.distance(p) <= tol.confusion() && x.distance(q) <= tol.confusion())
                })
            });
            let Some(&(_, _, bi)) = found else {
                return Ok(None);
            };
            let (edge, forward) =
                section(&mut points, &mut edges, &mut made, &bevels[bi], None, p, q)?;
            ring.push((edge, forward, bi));
        }
        patch_rings.push((*v, ring, ogeom_algo::HalfSpaces::Every(regions)));
    }
    // The faces a rolled round split off, then its strips' further pieces,
    // after every other face.
    let first_split = rings_of.len();
    rings_of.extend(split_faces.iter().map(|(_, r)| r.clone()));
    let first_rolled = rings_of.len();
    rings_of.extend(rolled_faces.iter().map(|(_, r)| r.clone()));
    // A consumed face is not built, nor a strip joined into another: the
    // sheet is every other.
    let away = |i: usize| {
        consumed.contains(&i)
            || (i >= first_strip
                && i < first_strip + bevels.len()
                && joined[i - first_strip] != i - first_strip)
    };
    let built_rings: Vec<SheetRings> = rings_of
        .iter()
        .enumerate()
        .filter(|(i, _)| !away(*i))
        .map(|(_, r)| r.clone())
        .collect();
    let mut edge_widths = vec![0.0; edges.len()];
    for (k, w) in marched_widths.borrow().iter() {
        edge_widths[*k] = *w;
    }
    for (k, w) in &feet_widths {
        if *w > tol.confusion() {
            edge_widths[*k] = edge_widths[*k].max(*w);
        }
    }
    // A bevel's end meeting a fitted contact or ridge stands on it within
    // its width: so is its edge built.
    let mut at_ridge: FastMap<usize, f64> = FastMap::default();
    for b in bevels {
        let e = &curved.edges[&b.edge];
        for fi in b.faces {
            if let Some(w) = widths.get(&(fi, b.edge)) {
                for v in [e.ends.0, e.ends.1] {
                    if let Some(c) = corners.get(&(fi, v)) {
                        let i = points.index(c.point);
                        let held = at_ridge.entry(i).or_insert(0.0);
                        *held = held.max(*w);
                    }
                }
            }
        }
    }
    for (key, k) in &made {
        if let Key::Section(i, j) = key {
            let w = at_ridge
                .get(i)
                .copied()
                .unwrap_or(0.0)
                .max(at_ridge.get(j).copied().unwrap_or(0.0));
            edge_widths[*k] = edge_widths[*k].max(w);
        }
    }
    // An edge left sharp keeps the tolerance it was read at, and its
    // vertices': an inexact one (a boolean's fitted section) passes its ends
    // no closer than that.
    let held = |shape: &Shape| -> OgeomResult<f64> {
        Ok(model
            .tolerance_of(shape)?
            .map_or(0.0, ogeom_core::Tolerance::get))
    };
    for (key, k) in &made {
        if let Key::Sharp(e) | Key::Tail(e, _) = key {
            let c = &curved.edges[e];
            let mut w = held(&c.shape)?;
            for v in [c.ends.0, c.ends.1] {
                if let Some(shape) = curved.shapes.get(&v) {
                    w = w.max(held(shape)?);
                }
            }
            if w > tol.confusion() {
                edge_widths[*k] = edge_widths[*k].max(w);
            }
        }
    }
    for (key, k) in &made {
        if let Key::Contact(fi, e) = key
            && let Some(w) = widths.get(&(*fi, *e))
        {
            edge_widths[*k] = *w;
        }
    }
    // A fitted edge reaches its corners within how far its curve's ends
    // stand from them (two contacts meeting midway where they nearly
    // touch), as far as the bands' fit target.
    let fit = crate::marched::band_fit_target(tol);
    for (k, e) in edges.iter().enumerate() {
        if edge_widths[k] <= 0.0 {
            continue;
        }
        let (a, b) = (
            e.curve.point_at(e.range.0, tol)?,
            e.curve.point_at(e.range.1, tol)?,
        );
        let off = a
            .distance(points.points[e.ends.0])
            .max(b.distance(points.points[e.ends.1]));
        if off <= fit {
            edge_widths[k] = edge_widths[k].max(off + tol.confusion());
        }
    }
    // Every edge is built to how far it stands off the faces it bounds,
    // measured, as far as the bands' fit target: a fitted contact's,
    // a marched section's or a bevel's end between two fitted contacts'
    // corners where they nearly touch. Sampled, so a quarter more is
    // allowed for what the samples miss between them.
    {
        const SAMPLES: u32 = 16;
        let mut used: Vec<Vec<&SurfaceGeometry>> = vec![Vec::new(); edges.len()];
        for r in &built_rings {
            for ring in &r.rings {
                for (k, _) in ring {
                    used[*k].push(&r.surface);
                }
            }
        }
        for (k, e) in edges.iter().enumerate() {
            let mut off: f64 = 0.0;
            for i in 0..=SAMPLES {
                let t =
                    (e.range.1 - e.range.0).mul_add(f64::from(i) / f64::from(SAMPLES), e.range.0);
                let p = e.curve.point_at(t, tol)?;
                for surface in &used[k] {
                    let foot = ogeom_geom::project::project_on_surface(surface, p, 16, tol)?;
                    off = off.max(foot.distance);
                }
            }
            if off > tol.confusion() && off <= fit {
                edge_widths[k] = edge_widths[k].max(1.25 * off + tol.confusion());
            }
        }
    }
    // Every edge of a round's fitted piece (a rolled run's, or a strip with
    // no closed form) is held to the width its chart's curve is fitted to
    // there.
    let fitted_rounds = !rolled.is_empty()
        || bevels.iter().any(|b| {
            matches!(b.section, Section::Arc { .. })
                && matches!(b.surface, SurfaceGeometry::BSpline(_))
        });
    for rings in built_rings
        .iter()
        .filter(|r| matches!(r.surface, SurfaceGeometry::BSpline(_)))
        .filter(|_| fitted_rounds)
    {
        for (k, _) in rings.rings.iter().flatten() {
            let off = rolled.iter().fold(0.0_f64, |m, r| m.max(r.run.off));
            edge_widths[*k] = edge_widths[*k].max(ROLLED_EDGE_WIDTH).max(2.0 * off);
        }
    }
    // A fitted edge running to an end of its curve is moved to end on its
    // corner there exactly (its end control point, which a clamped curve
    // passes through, shifted by the little it stands off), so its curve on
    // each face meets the next edge's at the corner, not beside it.
    for (k, e) in edges.iter_mut().enumerate() {
        if edge_widths[k] <= 0.0 {
            continue;
        }
        let Curve::BSpline(spline) = &e.curve else {
            continue;
        };
        if spline.is_rational() || spline.is_periodic() {
            continue;
        }
        let (lo, hi) = spline.domain();
        let mut control: Vec<Point> = spline.control_points().iter().map(|w| w.point()).collect();
        let n = control.len();
        // A closed curve keeps closed: both its ends move together, onto
        // the one corner at its start and end.
        let closed = spline
            .point_at(lo, tol)?
            .distance(spline.point_at(hi, tol)?)
            <= tol.confusion();
        let mut moved = false;
        for (t, at) in [(e.range.0, e.ends.0), (e.range.1, e.ends.1)] {
            if (t - lo).abs() > tol.parametric() && (t - hi).abs() > tol.parametric() {
                continue;
            }
            let end = if (t - lo).abs() <= tol.parametric() {
                0
            } else {
                n - 1
            };
            let off = points.points[at] - spline.point_at(t, tol)?;
            if off.magnitude() > tol.confusion() * 1e-3 && off.magnitude() <= edge_widths[k] {
                if closed {
                    control[0] = spline.control_points()[0].point() + off;
                    control[n - 1] = spline.control_points()[n - 1].point() + off;
                } else {
                    control[end] += off;
                }
                moved = true;
            }
        }
        if moved {
            e.curve = ogeom_geom::BSplineCurve::new(spline.knots().clone(), control, tol)?.into();
        }
    }
    let mut sheet = super::curved::build_with_widths(
        model,
        &points.points,
        &edges,
        &built_rings,
        &edge_widths,
        tol,
    )?;
    // Each corner of five bevels or more off one plane: the patch through
    // its corners on the bevels' end segments, held to the bands' fit
    // target where it is interpolated.
    let fit = crate::marched::band_fit_target(tol);
    let mut polygon_faces: Vec<(Id, Vec<Shape>)> = Vec::with_capacity(polygons.len());
    for (v, cycle) in &polygons {
        let index = |p: Point| {
            points
                .points
                .iter()
                .position(|x| x.distance(p) <= tol.confusion())
        };
        let mut sides = Vec::with_capacity(cycle.len());
        for k in 0..cycle.len() {
            let (p, q) = (cycle[k], cycle[(k + 1) % cycle.len()]);
            let (Some(i), Some(j)) = (index(p), index(q)) else {
                return Ok(None);
            };
            let Some(&e) = made.get(&Key::Section(i.min(j), i.max(j))) else {
                return Ok(None);
            };
            let strip = sections[v].iter().find_map(|(x, y, bi)| {
                let at = |a: Point, b: Point| {
                    a.distance(p) <= tol.confusion() && b.distance(q) <= tol.confusion()
                };
                (at(*x, *y) || at(*y, *x)).then_some(&bevels[*bi])
            });
            let Some(strip) = strip else {
                return Ok(None);
            };
            sides.push(crate::round::transfinite::HoleSide {
                edge: sheet.edges[e].clone(),
                curve: edges[e].curve.clone(),
                range: edges[e].range,
                forward: edges[e].ends.0 == i,
                support: strip.surface.clone(),
                pieces: Vec::new(),
            });
        }
        let patch =
            match crate::round::transfinite::corner_patch_through_corners(model, &sides, tol) {
                Ok(patch) => patch,
                Err(e @ ogeom_core::OgeomError::Cancelled) => return Err(e),
                Err(_) => return Ok(None),
            };
        if !(patch.gap <= fit && patch.off <= fit) {
            return Ok(None);
        }
        polygon_faces.push((*v, patch.faces));
    }
    // The patches, each through its rounds' ends as built, tangent to them,
    // held to the bands' fit target, or the request is handed back.
    let mut vertex_faces: Vec<(Id, Vec<Shape>, ogeom_algo::HalfSpaces)> =
        Vec::with_capacity(patch_rings.len());
    for (v, ring, region) in patch_rings {
        let sides: Vec<crate::round::transfinite::HoleSide> = ring
            .iter()
            .map(|&(e, forward, bi)| crate::round::transfinite::HoleSide {
                edge: sheet.edges[e].clone(),
                curve: edges[e].curve.clone(),
                range: edges[e].range,
                forward,
                support: bevels[bi].surface.clone(),
                pieces: Vec::new(),
            })
            .collect();
        let patch = crate::round::transfinite::corner_patch(model, &sides, tol)?;
        if !(patch.gap <= fit && patch.angle <= fit) {
            return Ok(None);
        }
        vertex_faces.push((v, patch.faces, region));
    }
    // Each cusp's patch, through its four sides wound counter-clockwise seen
    // from outside, proved where the plate's corner, the boss or the foot's
    // fill lies: within the top and the side, within the wall above the
    // top, or outside both the top and the wall.
    for (ci, mut sides) in holes {
        let c = &cusps[ci];
        let [top, side, wall] = c.faces;
        let at = curved.points[&c.v];
        let (Some(n0), Some(n1)) = (outward(curved, side, at), outward(curved, wall, at)) else {
            return Ok(None);
        };
        let mut ring = Vec::with_capacity(4);
        for spec in &sides {
            ring.push(spec.start(tol)?);
        }
        let mut newell = Vector::ZERO;
        for k in 0..ring.len() {
            newell += ring[k]
                .to_vector()
                .cross(ring[(k + 1) % ring.len()].to_vector());
        }
        if newell.dot(n0 + n1) < 0.0 {
            sides.reverse();
            sides = sides.into_iter().map(HoleSpec::back).collect();
            // A lobe's hole keeps its lobe first, so its feet's corner, the
            // cusp, is the patch's pole.
            if sides.len() == 3 {
                sides.rotate_right(1);
            }
        }
        let hole: Vec<crate::round::transfinite::HoleSide> =
            sides.iter().map(|spec| spec.side(&sheet.edges)).collect();
        let patch = crate::round::transfinite::corner_patch(model, &hole, tol)?;
        if !(patch.gap <= fit && patch.angle <= fit) {
            return Ok(None);
        }
        let (Some(mt), Some(ms), Some(mw)) = (
            material(curved, top),
            material(curved, side),
            material(curved, wall),
        ) else {
            return Ok(None);
        };
        use ogeom_algo::HalfSpaces;
        let region = HalfSpaces::Any(vec![
            HalfSpaces::Every(vec![mt.clone(), ms]),
            HalfSpaces::Every(vec![mw.clone(), outside(&mt)]),
            HalfSpaces::Every(vec![outside(&mt), outside(&mw)]),
        ]);
        vertex_faces.push((c.v, patch.faces, region));
    }
    let shell_faces: Vec<Shape> = sheet
        .faces
        .iter()
        .chain(polygon_faces.iter().flat_map(|(_, f)| f))
        .chain(vertex_faces.iter().flat_map(|(_, f, _)| f))
        .cloned()
        .collect();
    let shell = ogeom_algo::make_shell(model, &shell_faces)?.shape;
    // The sheet's faces again by the rings' indices, a consumed face's slot
    // holding the face's own shape, never used.
    let mut by_index: Vec<Shape> = Vec::with_capacity(rings_of.len());
    let mut next = sheet.faces.iter();
    for i in 0..rings_of.len() {
        match curved.faces.get(i) {
            Some(face) if consumed.contains(&i) => by_index.push(face.shape.clone()),
            _ if away(i) => by_index.push(curved.faces[0].shape.clone()),
            _ => by_index.push(
                next.next()
                    .cloned()
                    .unwrap_or_else(|| curved.faces[0].shape.clone()),
            ),
        }
    }
    // A strip joined into another stands for that one's face.
    for bi in 0..bevels.len() {
        let root = find(&joined, bi);
        if root != bi {
            by_index[first_strip + bi] = by_index[first_strip + root].clone();
        }
    }
    sheet.faces = by_index;
    let solid = ogeom_algo::make_solid(model, &[shell])?.shape;
    // An edge a face walks twice, once each way (a wall's or a cone
    // bevel's seam), needs both columns of the seam: the builder fits it
    // one pcurve, and a later walk of the face in its chart would not
    // close. Nothing else holds these nodes, so nothing is copied.
    ogeom_algo::join_seam_columns(model, &solid, tol)?;
    // P5: each bevel proved within the material its two faces bound there,
    // to the bands' fit target: under a plane, inside a wall (outside a
    // bore's).
    for (bi, b) in bevels.iter().enumerate() {
        // An extruded wall has no region the certify reads: its bevel's
        // rulings stand in their sections on its material side by
        // construction ([`extruded_crease`]), and only the cap is proved.
        let by_section = |k: usize| {
            matches!(
                curved.faces[b.faces[k]].surface,
                SurfaceGeometry::Extrusion(_)
            ) && matches!(b.surface, SurfaceGeometry::BSpline(_))
        };
        let face = &sheet.faces[curved.faces.len() + bi];
        // A face with no region the certify reads: an extruded wall holds
        // the bevel by construction; any other is measured, the bevel's
        // mesh standing on the face's material side at the foot of each
        // point, within the fit target.
        let mut pieces = Vec::with_capacity(2);
        for k in 0..2 {
            match material(curved, b.faces[k]) {
                Some(region) => pieces.push(region),
                None if by_section(k) => {}
                None if matches!(b.surface, SurfaceGeometry::BSpline(_)) => {
                    if !measured_within(model, face, &curved.faces[b.faces[k]], fit, tol)? {
                        return Ok(None);
                    }
                }
                None => return Ok(None),
            }
        }
        if pieces.is_empty() {
            continue;
        }
        // A round whose ball stands outside both faces' material (a concave
        // edge's) fills the corner: it stands outside each of them.
        let concave = pieces.len() == 2
            && match &b.section {
                Section::Arc { spine } => {
                    let (lo, hi) = spine.domain();
                    let at = spine.point_at(0.5 * (lo + hi), tol)?;
                    pieces[0].height(at).0 > 0.0 && pieces[1].height(at).0 > 0.0
                }
                Section::Line => false,
            };
        let region = if concave {
            ogeom_algo::HalfSpaces::Every(pieces.iter().map(outside).collect())
        } else {
            ogeom_algo::HalfSpaces::Every(pieces)
        };
        if ogeom_algo::certify_face_within(model, face, &region, fit, tol)?.is_err() {
            return Ok(None);
        }
        // A rolled strip's further pieces: those between its faces in the
        // same region; the run's piece within the other face's material and
        // the material the edge it rests on bounds (each face's, both at a
        // convex edge, either at a concave one), outside both where the
        // round fills.
        for (n, (strip, rings)) in rolled_faces.iter().enumerate() {
            if bevels[*strip].edge != b.edge {
                continue;
            }
            let piece = &sheet.faces[first_rolled + n];
            let within = if !matches!(rings.surface, SurfaceGeometry::BSpline(_)) {
                region.clone()
            } else {
                let Some(r) = rolled.iter().find(|r| r.strip == *strip) else {
                    return Ok(None);
                };
                let Some(edge_region) = resting_region(curved, r, tol)? else {
                    return Ok(None);
                };
                let Some(o) = material(curved, b.faces[1 - r.side]) else {
                    return Ok(None);
                };
                if concave {
                    ogeom_algo::HalfSpaces::Every(vec![outside(&o), outside(&edge_region)])
                } else {
                    ogeom_algo::HalfSpaces::Every(vec![o, edge_region])
                }
            };
            if ogeom_algo::certify_face_within(model, piece, &within, fit, tol)?.is_err() {
                return Ok(None);
            }
        }
    }
    // Each ball proved within the material of the faces round its vertex.
    let first_ball = curved.faces.len() + bevels.len() + facet_rings.len();
    for (k, region) in ball_regions.iter().enumerate() {
        if ogeom_algo::certify_face_within(model, &sheet.faces[first_ball + k], region, fit, tol)?
            .is_err()
        {
            return Ok(None);
        }
    }
    // Each patch within the material every face at its vertex bounds.
    let patch_faces: Vec<(Id, Shape)> = facet_rings
        .iter()
        .enumerate()
        .filter(|(_, (v, _))| patches.contains_key(v))
        .map(|(fi, (v, _))| {
            (
                *v,
                sheet.faces[curved.faces.len() + bevels.len() + fi].clone(),
            )
        })
        .chain(
            polygon_faces
                .iter()
                .flat_map(|(v, faces)| faces.iter().map(move |f| (*v, f.clone()))),
        )
        .collect();
    for (v, face) in &patch_faces {
        let mut round: Vec<usize> = curved
            .edges
            .iter()
            .filter(|(_, c)| c.ends.0 == *v || c.ends.1 == *v)
            .flat_map(|(e, _)| faces_of(curved, *e))
            .collect();
        round.sort_unstable();
        round.dedup();
        let Some(pieces) = round
            .iter()
            .map(|fi| material(curved, *fi))
            .collect::<Option<Vec<_>>>()
        else {
            return Ok(None);
        };
        let region = ogeom_algo::HalfSpaces::Every(pieces);
        if ogeom_algo::certify_face_within(model, face, &region, fit, tol)?.is_err() {
            return Ok(None);
        }
    }
    // History: faces modified, bevelled edges deleted into their bevels,
    // sharp edges modified into their pieces, vertices kept where nothing
    // moved off them.
    let mut history = History::new();
    for (fi, face) in curved.faces.iter().enumerate() {
        if consumed.contains(&fi) {
            history.delete(&face.shape);
        } else {
            history.modify(&face.shape, sheet.faces[fi].clone());
        }
    }
    for (bi, b) in bevels.iter().enumerate() {
        let edge = &curved.shapes[&b.edge];
        history.delete(edge);
        history.generate(edge, sheet.faces[curved.faces.len() + bi].clone());
    }
    for (fi, (v, _)) in facet_rings.iter().enumerate() {
        history.generate(
            &curved.shapes[v],
            sheet.faces[curved.faces.len() + bevels.len() + fi].clone(),
        );
    }
    for (v, faces) in &polygon_faces {
        for face in faces {
            history.generate(&curved.shapes[v], face.clone());
        }
    }
    for (k, (v, _, _, _)) in balls.iter().enumerate() {
        history.generate(&curved.shapes[v], sheet.faces[first_ball + k].clone());
    }
    // Each patch proved within the material of the faces round its vertex,
    // and made from it.
    for (v, faces, region) in &vertex_faces {
        for face in faces {
            if ogeom_algo::certify_face_within(model, face, region, fit, tol)?.is_err() {
                return Ok(None);
            }
            history.generate(&curved.shapes[v], face.clone());
        }
    }
    for (e, _) in &curved.edges {
        if let Some(k) = made.get(&Key::Sharp(*e)) {
            history.modify(&curved.shapes[e], sheet.edges[*k].clone());
        }
    }
    // A rolled round's faces split off theirs, its strips' further pieces
    // made from their edges, and each edge it cut modified into its pieces.
    for (n, (fi, _)) in split_faces.iter().enumerate() {
        history.modify(
            &curved.faces[*fi].shape,
            sheet.faces[first_split + n].clone(),
        );
    }
    for (n, (strip, _)) in rolled_faces.iter().enumerate() {
        history.generate(
            &curved.shapes[&bevels[*strip].edge],
            sheet.faces[first_rolled + n].clone(),
        );
    }
    let mut cut: Vec<(Id, usize)> = made
        .iter()
        .filter_map(|(key, k)| match key {
            Key::Piece(usize::MAX, e, _, _) => Some((*e, *k)),
            _ => None,
        })
        .collect();
    cut.sort_unstable();
    for (e, k) in cut {
        history.modify(&curved.shapes[&e], sheet.edges[k].clone());
    }
    for (v, p) in &curved.points {
        let kept = curved.faces.iter().enumerate().all(|(fi, _)| {
            [*v, pinched(*v)].iter().all(|w| {
                corners
                    .get(&(fi, *w))
                    .is_none_or(|c| c.point.distance(*p) <= tol.confusion())
            })
        });
        let shape = &curved.shapes[v];
        if kept
            && let Some(i) = points
                .points
                .iter()
                .position(|q| q.distance(*p) <= tol.confusion())
        {
            history.modify(shape, sheet.vertices[i].clone());
        } else {
            history.delete(shape);
        }
    }
    Ok(Some(Built::new(solid, history)))
}

/// The bevels' ends at a vertex chained into one cycle of corners, each
/// end joining the corners of two faces; `None` where they do not close one
/// cycle through every corner once.
fn chained(ends: &[(Point, Point, usize)], tol: Tolerances) -> Option<Vec<Point>> {
    let same = |a: Point, b: Point| a.distance(b) <= tol.confusion();
    let mut cycle = vec![ends[0].0, ends[0].1];
    let mut used = vec![false; ends.len()];
    used[0] = true;
    loop {
        let last = *cycle.last()?;
        let next = ends
            .iter()
            .enumerate()
            .find(|(i, (p, q, _))| !used[*i] && (same(*p, last) || same(*q, last)))
            .map(|(i, (p, q, _))| (i, if same(*p, last) { *q } else { *p }));
        let Some((i, far)) = next else {
            break;
        };
        used[i] = true;
        if same(far, cycle[0]) {
            return used.iter().all(|u| *u).then_some(cycle);
        }
        cycle.push(far);
    }
    None
}

/// Which way a sheet edge, walked as stored, turns about a revolved strip's
/// axis: the sign of its tangent's turn at its start, `None` off a cone, a
/// cylinder or a torus.
fn turning(
    surface: &SurfaceGeometry,
    edge: &SheetCurve,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let (o, a) = match surface {
        SurfaceGeometry::Cone(c) => {
            let f = c.cone().frame();
            (f.origin(), f.z().vector())
        }
        SurfaceGeometry::Cylinder(c) => {
            let axis = c.cylinder().axis();
            (axis.location, axis.direction.vector())
        }
        SurfaceGeometry::Torus(t) => {
            let f = t.torus().frame();
            (f.origin(), f.z().vector())
        }
        _ => return Ok(None),
    };
    let t = edge.range.0;
    let (p, d) = (edge.curve.point_at(t, tol)?, edge.curve.d1_at(t, tol)?);
    Ok(Some((p - o).cross(d).dot(a).signum()))
}
