//! The filler's records: sections, contacts, tangencies, grazes, paves, and
//! the chart folding they share.

use super::*;

// --- the filler --------------------------------------------------------------

/// One section curve between a face of `a` and a face of `b`.
pub(crate) struct SectionRec {
    pub(crate) curve: Curve,
    /// The pcurve on the A face's surface, sharing the curve's parameter.
    pub(crate) pc_a: PlanarCurve,
    /// The same on the B face.
    pub(crate) pc_b: PlanarCurve,
    pub(crate) face_a: usize,
    pub(crate) face_b: usize,
    pub(crate) closed: bool,
    /// How far the curve may sit from the true intersection: zero for an
    /// exact section, the trace-plus-fit budget for a marched one. Crossing
    /// filters widen their acceptance by this, or a fitted section would
    /// never register against the edges it genuinely meets.
    pub(crate) tolerance: f64,
}

/// A closed circular section whose seam stands a hair from a vertex of
/// either face on it, turned to start at that vertex.
///
/// The section's own start is wherever its frame happened to point (a
/// world axis projected into its plane), and a piece reaching across it is
/// split there to stay within one period. Where a face's vertex sits on
/// the circle a sliver of a turn from that start, the split leaves an arc
/// the paving takes for dust and collapses, dragging the vertex along the
/// circle. Starting the circle at the vertex leaves nothing to split. A
/// seam already on a vertex, or well clear of every one, is left alone.
pub(crate) fn seamed_at_a_stop(
    section: &SectionRec,
    fa: &GFace,
    fb: &GFace,
    tol: Tolerances,
) -> Option<SectionRec> {
    use core::f64::consts::TAU;
    if !section.closed || section.tolerance > 0.0 {
        return None;
    }
    let Curve::Circle(c) = &section.curve else {
        return None;
    };
    let circle = c.circle();
    let (centre, radius) = (circle.centre(), circle.radius());
    let frame = circle.frame();
    // The angle a point stands at along the curve's own parameter.
    let sweep = if c.is_reversed() { -1.0 } else { 1.0 };
    let angle_of = |p: Point| {
        let local = frame.to_local(p);
        (sweep * local.y.atan2(local.x)).rem_euclid(TAU)
    };
    let mut nearest: Option<f64> = None;
    for face in [fa, fb] {
        for e in &face.edges {
            for end in [e.crange.0, e.crange.1] {
                let Ok(p) = e.curve.point_at(end, tol) else {
                    continue;
                };
                let local = frame.to_local(p);
                let off = local.z.hypot(local.x.hypot(local.y) - radius);
                if off > tol.confusion() * 10.0 {
                    continue;
                }
                let a = angle_of(p);
                let shift = if a > TAU / 2.0 { a - TAU } else { a };
                if nearest.is_none_or(|n| shift.abs() < n.abs()) {
                    nearest = Some(shift);
                }
            }
        }
    }
    let shift = nearest?;
    if shift.abs() <= tol.parametric() || shift.abs() > 1e-2 {
        return None;
    }
    let start = c.point_at(shift, tol).ok()?;
    let x = ogeom_math::Direction::new(start - centre, tol).ok()?;
    let z = if c.is_reversed() {
        -frame.z()
    } else {
        frame.z()
    };
    let turned = ogeom_math::Frame::new(centre, z, x, tol).ok()?;
    let curve: Curve =
        ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(turned, radius, tol).ok()?).into();
    let pc_a = ogeom_intersect::exact_pcurve_of(&curve, &fa.surface, tol)?;
    let pc_b = ogeom_intersect::exact_pcurve_of(&curve, &fb.surface, tol)?;
    Some(SectionRec {
        curve,
        pc_a,
        pc_b,
        face_a: section.face_a,
        face_b: section.face_b,
        closed: true,
        tolerance: section.tolerance,
    })
}

/// A closed marched section with its seam moved to `u`, its curve and both
/// images the same from the new start, each parameter naming the point it
/// named: the domain begins at `u`, and past the old end each image is the
/// old one a lap earlier, carried by the whole periods it climbs over the
/// lap.
///
/// A closed section's kept piece that runs across its seam is split there
/// to stay within one period. Moved to the piece's start, the seam splits
/// nothing.
///
/// `None` (the piece split at the seam, as before) unless the curve is a
/// clamped spline closed to the confusion and each image is a spline whose
/// end stands on its start, carried by whole periods of its chart, to the
/// parametric tolerance: then nothing moves but the seam.
pub(crate) fn reseamed_section(
    section: &SectionRec,
    u: f64,
    surfaces: [&SurfaceGeometry; 2],
    tol: Tolerances,
) -> Option<SectionRec> {
    let Curve::BSpline(spline) = &section.curve else {
        return None;
    };
    if spline.is_periodic() {
        return None;
    }
    let (start, end) = spline.domain();
    if spline
        .point_at(start, tol)
        .ok()?
        .distance(spline.point_at(end, tol).ok()?)
        > tol.confusion()
    {
        return None;
    }
    let image = |pcurve: &PlanarCurve, surface: &SurfaceGeometry| -> Option<PlanarCurve> {
        let PlanarCurve::BSpline(image) = pcurve else {
            return None;
        };
        let domain = image.domain();
        if (domain.0 - start).abs() > tol.parametric() || (domain.1 - end).abs() > tol.parametric()
        {
            return None;
        }
        let (head, tail) = (
            image.point_at(start, tol).ok()?,
            image.point_at(end, tol).ok()?,
        );
        // The lap is the whole periods nearest the image's climb; the
        // reseam itself refuses an image that misses it.
        let climb = [tail.x - head.x, tail.y - head.y];
        let mut lap = [0.0; 2];
        for (k, span) in chart_spans(surface, tol).into_iter().enumerate() {
            if let Some(span) = span {
                lap[k] = (climb[k] / span).round() * span;
            }
        }
        image
            .reseamed_at(u, ogeom_math::Vector2::new(lap[0], lap[1]), tol)
            .ok()
            .map(PlanarCurve::BSpline)
    };
    Some(SectionRec {
        curve: Curve::BSpline(spline.reseamed_at(u, tol).ok()?),
        pc_a: image(&section.pc_a, surfaces[0])?,
        pc_b: image(&section.pc_b, surfaces[1])?,
        face_a: section.face_a,
        face_b: section.face_b,
        closed: true,
        tolerance: section.tolerance,
    })
}

/// A face's welded chart outline, borrowed as the trim tests take it.
pub(crate) fn outline_refs(face: &GFace, tol: Tolerances) -> OgeomResult<Vec<&[Point2]>> {
    let Some(lines) = face.outline(tol) else {
        ogeom_bail!(
            Construction,
            "a face's boundary would not polyline in its chart"
        );
    };
    Ok(lines.iter().map(Vec::as_slice).collect())
}

/// The parameter intervals of a contact curve over which *both* faces
/// actually reach it.
///
/// Two surfaces touch along the whole of their contact. Two faces touch
/// along whatever part of it their trims both hold. That part is found by
/// sampling (sixty-four stations along the curve, each asked of both
/// charts) rather than by intersecting the contact with the boundary
/// edges, because a contact meets those boundaries tangentially too and the
/// crossing finder is the wrong instrument for it. The cost of sampling is
/// the usual one: a stretch shorter than a station can be missed, and an
/// endpoint is placed within a station of the truth.
pub(crate) fn contact_intervals(
    fused: &GeneralFused,
    contact: &TangentRec,
    tol: Tolerances,
) -> OgeomResult<Vec<(f64, f64)>> {
    // The same welded outline, both seam columns, that every other trim
    // test reads: with one column a full drum's trim is open, and every
    // station of a contact on it reads outside.
    let refs_a = outline_refs(&fused.a.faces[contact.face_a], tol)?;
    let refs_b = outline_refs(&fused.b.faces[contact.face_b], tol)?;

    const STATIONS: usize = 64;
    let domain = contact.curve.domain();
    let span = domain.1 - domain.0;
    let mut runs: Vec<(f64, f64)> = Vec::new();
    let mut open: Option<f64> = None;
    for k in 0..=STATIONS {
        #[expect(
            clippy::cast_precision_loss,
            reason = "a station index, far below the mantissa"
        )]
        let t = span.mul_add(k as f64 / STATIONS as f64, domain.0);
        // A periodic chart's parameters run out past its own window. The
        // trim test is only meaningful once they are folded back into it.
        let held = matches!(
            (
                contact.pc_a.point_at(t, tol),
                contact.pc_b.point_at(t, tol)
            ),
            (Ok(pa), Ok(pb))
                if arrange::inside_many(
                    &refs_a,
                    fold_point_into_chart(pa, &fused.a.faces[contact.face_a].surface, tol),
                ) && arrange::inside_many(
                    &refs_b,
                    fold_point_into_chart(pb, &fused.b.faces[contact.face_b].surface, tol),
                )
        );
        match (held, open) {
            (true, None) => open = Some(t),
            (false, Some(from)) => {
                if t - from > tol.parametric() {
                    runs.push((from, t));
                }
                open = None;
            }
            _ => {}
        }
    }
    if let Some(from) = open
        && domain.1 - from > tol.parametric()
    {
        runs.push((from, domain.1));
    }
    Ok(runs)
}

/// One curve along which two faces *touch* without crossing.
///
/// A contact carries no boundary parity (neither face passes through the
/// other), so it takes no part in the arrangement or the classification.
/// It is still a curve that exists on the result, and a section view
/// through a tangency has to show it, so it is carried alongside.
pub(crate) struct TangentRec {
    pub(crate) curve: Curve,
    pub(crate) pc_a: PlanarCurve,
    pub(crate) pc_b: PlanarCurve,
    pub(crate) face_a: usize,
    pub(crate) face_b: usize,
}

/// A marched section of faces `fa` and `fb` with each end that stopped
/// short of a vertex of either face moved onto it, and how far it moved.
///
/// Where a section runs into a point both surfaces pass through nearly
/// tangent to each other (a transfinite patch's point side lying on the
/// plane it meets), the surfaces stay within the marcher's tolerance of
/// each other for a stretch, and the trace stops somewhere along it, short
/// of the point: its end then stands off the vertex the other curves end
/// at, and the wire through it does not close. An end within
/// [`SNAP_REACH`] of the faces' chord scale (and past [`SNAP_FLOOR`]
/// confusions, which the arrangement welds) of a vertex that is a pole of
/// either surface ([`SNAP_POLE`]) lying on both
/// surfaces (to within [`SNAP_ON`] confusions), the surfaces nearly
/// touching at the end ([`SNAP_TOUCH`]), is moved there: the curve's
/// end control point, and each pcurve's, to the vertex's parameters found
/// from the end's own (a pole keeps its side's parameter). The section's
/// tolerance widens by the move. Curves or pcurves other than clamped
/// B-splines are left as they are.
pub(crate) fn snapped_to_vertices(
    fa: &GFace,
    fb: &GFace,
    curve: Curve,
    pa: PlanarCurve,
    pb: PlanarCurve,
    tol: Tolerances,
) -> OgeomResult<(Curve, PlanarCurve, PlanarCurve, f64)> {
    let (Curve::BSpline(c), PlanarCurve::BSpline(qa), PlanarCurve::BSpline(qb)) =
        (&curve, &pa, &pb)
    else {
        return Ok((curve, pa, pb, 0.0));
    };
    let reach = SNAP_REACH * fa.chord_scale.min(fb.chord_scale);
    let on = SNAP_ON * tol.confusion();
    let mut vertices: Vec<Point> = Vec::new();
    for e in fa.edges.iter().chain(&fb.edges) {
        for t in [e.crange.0, e.crange.1] {
            vertices.push(e.curve.point_at(t, tol)?);
        }
    }
    let (lo, hi) = curve.domain();
    let mut control: Vec<Weighted<Point>> = c.control_points().to_vec();
    let mut control_a: Vec<Weighted<Point2>> = qa.control_points().to_vec();
    let mut control_b: Vec<Weighted<Point2>> = qb.control_points().to_vec();
    let mut moved = 0.0_f64;
    for (end, t) in [(0usize, lo), (1, hi)] {
        let at = curve.point_at(t, tol)?;
        // Only where the surfaces nearly touch at the end: a section
        // crossing an edge beside a vertex ends there truly.
        let mut normals = [Vector::ZERO; 2];
        for (k, (face, pcurve)) in [(fa, &pa), (fb, &pb)].into_iter().enumerate() {
            let (pl, ph) = pcurve.domain();
            let uv = pcurve.point_at(if end == 0 { pl } else { ph }, tol)?;
            let Ok(n) = face.surface.normal_at(uv.x, uv.y, tol) else {
                continue;
            };
            normals[k] = n.vector();
        }
        if normals[0].cross(normals[1]).magnitude() > SNAP_TOUCH {
            continue;
        }
        let Some(&vertex) = vertices
            .iter()
            .filter(|v| {
                let d = v.distance(at);
                d > SNAP_FLOOR * tol.confusion() && d <= reach
            })
            .min_by(|a, b| a.distance(at).total_cmp(&b.distance(at)))
        else {
            continue;
        };
        // The vertex on both surfaces, its parameters on each found from the
        // end's own.
        let mut images = [Point2::ORIGIN; 2];
        let mut ok = true;
        for (k, (face, pcurve)) in [(fa, &pa), (fb, &pb)].into_iter().enumerate() {
            let (pl, ph) = pcurve.domain();
            let guess = pcurve.point_at(if end == 0 { pl } else { ph }, tol)?;
            let foot = ogeom_algo::project_on_surface_from(
                &face.surface,
                vertex,
                (guess.x, guess.y),
                tol,
            )?;
            if foot.distance > on {
                ok = false;
                break;
            }
            images[k] = Point2::new(foot.parameters.0, foot.parameters.1);
        }
        if !ok {
            continue;
        }
        // The vertex a pole of either surface: a side of its chart
        // collapsed to the point, one derivative vanishing there.
        let mut pole = false;
        for (k, face) in [fa, fb].into_iter().enumerate() {
            let Ok((du, dv)) = face.surface.d1_at(images[k].x, images[k].y, tol) else {
                continue;
            };
            let (mu, mv) = (du.magnitude(), dv.magnitude());
            if mu.min(mv) <= SNAP_POLE * mu.max(mv) {
                pole = true;
            }
        }
        if !pole {
            continue;
        }
        let move_to = |points: &mut [Weighted<Point>], p: Point| {
            let k = if end == 0 { 0 } else { points.len() - 1 };
            let w = points[k].weight;
            points[k] = Weighted {
                scaled: Point::from_vector(p.to_vector() * w),
                weight: w,
            };
        };
        let move_to_2d = |points: &mut [Weighted<Point2>], p: Point2| {
            let k = if end == 0 { 0 } else { points.len() - 1 };
            let w = points[k].weight;
            points[k] = Weighted {
                scaled: Point2::new(p.x * w, p.y * w),
                weight: w,
            };
        };
        move_to(&mut control, vertex);
        move_to_2d(&mut control_a, images[0]);
        move_to_2d(&mut control_b, images[1]);
        moved = moved.max(vertex.distance(at));
        ogeom_core::diag!(Channel::Wire, "SNAP end {end} {at:?} -> {vertex:?}");
    }
    if moved == 0.0 {
        return Ok((curve, pa, pb, 0.0));
    }
    let (Ok(c), Ok(qa), Ok(qb)) = (
        ogeom_geom::BSplineCurve::rational(c.knots().clone(), control),
        ogeom_geom::BSpline2d::rational(qa.knots().clone(), control_a),
        ogeom_geom::BSpline2d::rational(qb.knots().clone(), control_b),
    ) else {
        return Ok((curve, pa, pb, 0.0));
    };
    Ok((
        Curve::BSpline(c),
        PlanarCurve::BSpline(qa),
        PlanarCurve::BSpline(qb),
        moved,
    ))
}

/// How far a marched section's end may stop short of a vertex both
/// surfaces pass through and still be moved onto it
/// ([`snapped_to_vertices`]), as a share of the faces' chord scale.
pub(crate) const SNAP_REACH: f64 = 1e-3;

/// How far off each surface, in confusions, a vertex a marched section's
/// end is moved onto may lie.
pub(crate) const SNAP_ON: f64 = 100.0;

/// How small one chart derivative must be against the other, at a vertex,
/// for the vertex to be a pole of the surface (a side collapsed to it).
pub(crate) const SNAP_POLE: f64 = 1e-6;

/// How far, in confusions, a marched section's end may stand off a vertex
/// and be left where it is: the arrangement welds the wire through such a
/// gap itself, other curves meeting the end where it stands.
pub(crate) const SNAP_FLOOR: f64 = 1e3;

/// The sine of the angle between the two surfaces' normals at a marched
/// section's end under which they nearly touch there, and the end may be
/// moved onto a vertex.
pub(crate) const SNAP_TOUCH: f64 = 1e-2;

/// Whether a line edge of one face, its ends at `ends`, runs along a line
/// edge of `target` end to end, each end within `weld` of that edge.
pub(crate) fn shared_line(
    ends: &[Point],
    target: &GFace,
    weld: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    for t in &target.edges {
        if !matches!(&*t.curve, Curve::Line(_)) {
            continue;
        }
        let (a, b) = (
            t.curve.point_at(t.crange.0, tol)?,
            t.curve.point_at(t.crange.1, tol)?,
        );
        let along = (b.x - a.x, b.y - a.y, b.z - a.z);
        let length2 = along.0 * along.0 + along.1 * along.1 + along.2 * along.2;
        if length2 <= 0.0 {
            continue;
        }
        let beside = |p: &Point| {
            let s =
                ((p.x - a.x) * along.0 + (p.y - a.y) * along.1 + (p.z - a.z) * along.2) / length2;
            let s = s.clamp(0.0, 1.0);
            p.distance(Point::new(
                a.x + s * along.0,
                a.y + s * along.1,
                a.z + s * along.2,
            ))
        };
        if ends.iter().all(|p| beside(p) <= weld) {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether a face bounded by straight edges lies on one side of `plane`, its
/// corners past it by rounding at most: `None` where it crosses, has no
/// edges, or has an edge of another kind or a pole, whose extent off the
/// plane its corners do not bound; otherwise whether an edge of it lies in
/// the plane, both ends within the confusion distance.
pub(crate) fn one_side_of(
    face: &GFace,
    plane: ogeom_math::quadric::Plane,
    tol: Tolerances,
) -> OgeomResult<Option<bool>> {
    if face.edges.is_empty() || !face.poles.is_empty() {
        return Ok(None);
    }
    let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut edge_in = false;
    for e in &face.edges {
        if !matches!(&*e.curve, Curve::Line(_)) {
            return Ok(None);
        }
        let mut ends_in = true;
        for t in [e.crange.0, e.crange.1] {
            let d = plane.signed_distance_to(e.curve.point_at(t, tol)?);
            (low, high) = (low.min(d), high.max(d));
            ends_in &= d.abs() <= tol.confusion();
        }
        edge_in |= ends_in;
    }
    // A corner the face shares with the plane's own face stands off it by
    // rounding alone. A corner past it by more, however little, is a
    // crossing: two blocks hinged a hair apart cross a few nanometres inside
    // a side face, and that line bounds the sliver between them.
    let past = tol.confusion() * 1e-3;
    Ok((low >= -past || high <= past).then_some(edge_in))
}

/// A strand's tolerance as a junction may trust it: a fitted section whose
/// trace failed reports a budget of metres, and a weld that believed it
/// would join every vertex of the model. Nothing this pipeline fits is
/// honestly looser than ten thousand confusions.
pub(crate) fn honest(tolerance: f64, tol: Tolerances) -> f64 {
    tolerance.min(tol.confusion() * 1e4)
}

/// How far off the other surface a boundary edge may lie and still be
/// where a grazing section runs: ten confusions, fixed.
pub(crate) const GRAZE_REACH: f64 = 10.0;

/// The stretches of boundary edge a grazing section runs along, carried
/// onto the other face as contacts.
///
/// Where a face's boundary edge lies in the other surface's tangent plane
/// and touches it at one point (a drill whose wall runs along a plane
/// through a rim the rim's circle touches, or along a round's tangent line
/// into the corner where the round's end arc touches it), their section
/// meets the edge there to the fourth order, and stays within the
/// confusion distance of the edge for a stretch round it. The marcher stops
/// on the stretch, and the section ends on the edge short of where it
/// meets it: both ends on one edge, the loop handed back open, or one end
/// on an edge short of the edge's vertex where the section truly ends. The
/// face holding the edge is split there by its own boundary, but the other
/// face is not split along the stretch, and its piece the section bounds
/// never closes. The edge's stretch from the end to the other end or to
/// the vertex, where it lies on the other surface to within
/// [`GRAZE_REACH`] confusions, is the rest of the section: carried as a
/// contact it splits the other face as the section would, and both sides
/// share it as the edge it is.
///
/// With `crossings`, on the retries run again once every other has refused
/// a piece read on the other solid (see [`settled`]), a section line of two
/// planes also ends where it crosses an edge of either face: a line edge,
/// and on the retries after those, any edge ([`section_crossings`]). Two
/// walls leaning a hair apart cross a few microns inside one's edge, and
/// the strip between stands within rounding of the other wall: the edges
/// bounding the strip lie on that wall from the crossing on (an arc from
/// one crossing to the other), and carried there they split the other
/// wall's matching strip off, so the two are read as one contact and one
/// of them kept. Every earlier attempt leaves the edges where they are, so
/// whatever those resolve is resolved as before.
pub(crate) fn graze_contacts(
    ga: &GSolid,
    gb: &GSolid,
    sections: &[SectionRec],
    contacts: &[ContactRec],
    crossings: Crossings,
    tol: Tolerances,
) -> OgeomResult<Vec<ContactRec>> {
    let reach = tol.confusion() * GRAZE_REACH;
    // A stretch shorter than the weld is the junction itself, which the
    // crossing at the edge's end already snaps onto the vertex.
    let weld = tol.confusion() * 1e2;
    let mut out: Vec<ContactRec> = Vec::new();
    for section in sections.iter().filter(|s| !s.closed) {
        let (lo, hi) = section.curve.domain();
        let ends = [
            section.curve.point_at(lo, tol)?,
            section.curve.point_at(hi, tol)?,
        ];
        let middle = section.curve.point_at(f64::midpoint(lo, hi), tol)?;
        let near = tol.confusion().max(honest(section.tolerance, tol));
        for (owner, target, target_from_a, target_face) in [
            (
                &ga.faces[section.face_a],
                &gb.faces[section.face_b],
                false,
                section.face_b,
            ),
            (
                &gb.faces[section.face_b],
                &ga.faces[section.face_a],
                true,
                section.face_a,
            ),
        ] {
            for e in &owner.edges {
                let held = |c: &ContactRec| {
                    c.node == e.node
                        && c.target_from_a == target_from_a
                        && c.target_face == target_face
                };
                let reached = near.max(e.standoff);
                // A section line whose ends stand off the edge may still
                // leave the owner face across it. Crossed once, the
                // crossing is its one end there, the other staying out of
                // reach, and its middle is read an edge's length along from
                // the crossing; crossed twice, the two crossings are its
                // ends and its middle is halfway between them.
                let crossed = if crossings != Crossings::Off
                    && ends.iter().all(|p| e.bound.distance_to(*p) > reached)
                {
                    section_crossings(&section.curve, e, &owner.surface, crossings, reached, tol)?
                } else {
                    None
                };
                let (ends, middle) = match crossed {
                    Some(Crossed::Once(at, along)) => {
                        let length = e
                            .curve
                            .point_at(e.crange.0, tol)?
                            .distance(e.curve.point_at(e.crange.1, tol)?);
                        ([at, ends[1]], at + along * length)
                    }
                    Some(Crossed::Twice(first, second)) => {
                        ([first, second], first.midpoint(second))
                    }
                    None => (ends, middle),
                };
                let crossed = crossed.is_some() && !matches!(&*e.curve, Curve::Line(_));
                if ends.iter().all(|p| e.bound.distance_to(*p) > reached)
                    || contacts.iter().chain(&out).any(held)
                {
                    continue;
                }
                let mut feet: [Option<f64>; 2] = [None, None];
                for (k, end) in ends.iter().enumerate() {
                    if e.bound.distance_to(*end) > reached {
                        continue;
                    }
                    let foot = ogeom_algo::project_on_curve(&e.curve, *end, 64, tol)?;
                    if foot.distance <= reached
                        && foot.parameter >= e.crange.0 - tol.parametric()
                        && foot.parameter <= e.crange.1 + tol.parametric()
                    {
                        feet[k] = Some(foot.parameter.clamp(e.crange.0, e.crange.1));
                    }
                }
                // A section running along the edge (a wall's seam lying in
                // a plane it is tangent to) is that edge, not one that
                // leaves it.
                if feet.iter().all(Option::is_none)
                    || ogeom_algo::project_on_curve(&e.curve, middle, 64, tol)?.distance <= near
                {
                    continue;
                }
                let closed = e
                    .curve
                    .point_at(e.crange.0, tol)?
                    .distance(e.curve.point_at(e.crange.1, tol)?)
                    <= tol.confusion();
                // The ways the stretch may run, each a list of ranges on the
                // edge: between two ends (on a closed edge either way round
                // it, across its vertex), or from one end to a vertex.
                let mut ways: Vec<Vec<(f64, f64)>> = Vec::new();
                if let [Some(f0), Some(f1)] = feet {
                    let (t0, t1) = (f0.min(f1), f0.max(f1));
                    if t1 - t0 > tol.parametric() {
                        ways.push(vec![(t0, t1)]);
                        if closed {
                            ways.push(
                                [(t1, e.crange.1), (e.crange.0, t0)]
                                    .into_iter()
                                    .filter(|(a, b)| b - a > tol.parametric())
                                    .collect(),
                            );
                        }
                    }
                } else if !closed {
                    for t in feet.into_iter().flatten() {
                        ways.push(vec![(e.crange.0, t)]);
                        ways.push(vec![(t, e.crange.1)]);
                    }
                }
                let mut stretches = None;
                for way in ways {
                    let mut lies = !way.is_empty();
                    'way: for &(a, b) in &way {
                        let span = e
                            .curve
                            .point_at(a, tol)?
                            .distance(e.curve.point_at(b, tol)?);
                        if b - a <= tol.parametric() || span <= weld {
                            lies = false;
                            break;
                        }
                        for i in 0..=16 {
                            let t = a + (b - a) * f64::from(i) / 16.0;
                            let p = e.curve.point_at(t, tol)?;
                            if ogeom_algo::project_on_surface(&target.surface, p, 16, tol)?.distance
                                > reach
                            {
                                lies = false;
                                break 'way;
                            }
                        }
                    }
                    if lies {
                        stretches = Some(way);
                        break;
                    }
                }
                let Some(stretches) = stretches else {
                    continue;
                };
                for crange in stretches {
                    let Ok((pcurve, _, true, ..)) = ogeom_algo::pcurve_fit::fit_projected_pcurve(
                        &e.curve,
                        crange,
                        &target.surface,
                        tol,
                    ) else {
                        continue;
                    };
                    ogeom_core::diag!(
                        Channel::Wire,
                        "GRAZE edge {} over {crange:?} onto face {target_face} (from_a {target_from_a})",
                        e.node.index()
                    );
                    out.push(ContactRec {
                        curve: (*e.curve).clone(),
                        crange,
                        pcurve,
                        prange: crange,
                        node: e.node,
                        tolerance: e.tolerance.max(owner.tolerance + target.tolerance),
                        target_from_a,
                        target_face,
                        bound: e.bound,
                        crossed,
                    });
                }
            }
        }
    }
    Ok(out)
}

/// Which edges [`graze_contacts`] carries onto the other face from where a
/// section line crosses them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Crossings {
    /// None: an edge is carried only from a section's own ends.
    Off,
    /// Line edges.
    Lines,
    /// Edges of any curve.
    Curves,
}

/// Where a section line crosses an edge of the face it lies on, as
/// [`section_crossings`] finds it.
pub(crate) enum Crossed {
    /// One crossing, on the line, with the line's unit direction.
    Once(Point, ogeom_math::Vector),
    /// Two crossings, on the line.
    Twice(Point, Point),
}

/// Where a section line crosses the boundary edge `e` of the face whose
/// surface is `owner`, within the edge's range: `None` unless the section
/// is a line, or where the edge is crossed other than once or twice.
///
/// A line edge is crossed where the two lines come nearest each other
/// ([`lines_crossing`]). Any other edge needs the face to be a plane, which
/// holds both it and the line. A circle is crossed where the line, laid in
/// the circle's plane, meets it: the two roots of a quadratic. Any other
/// curve is crossed where its signed distance from the line, square to the
/// line in the plane, changes sign: read at 64 steps along the range and
/// halved down to rounding at each change. Two crossings within one step
/// are not seen; what is found still carries a stretch only where the
/// stretch lies on the other face end to end. Each crossing is given on
/// the line, and kept only where the edge's nearest point to it is within
/// `reach` and inside the edge's range.
pub(crate) fn section_crossings(
    section: &Curve,
    e: &BoundaryEdge,
    owner: &SurfaceGeometry,
    crossings: Crossings,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Crossed>> {
    if matches!(&*e.curve, Curve::Line(_)) {
        return Ok(
            lines_crossing(section, &e.curve, tol)?.map(|(at, along)| Crossed::Once(at, along))
        );
    }
    if crossings != Crossings::Curves {
        return Ok(None);
    }
    let (Curve::Line(_), SurfaceGeometry::Plane(plane)) = (section, owner) else {
        return Ok(None);
    };
    let p = section.point_at(0.0, tol)?;
    let u = section.point_at(1.0, tol)? - p;
    let length = u.magnitude();
    if length <= 0.0 {
        return Ok(None);
    }
    let u = u / length;
    let mut found: Vec<Point> = Vec::new();
    if let Curve::Circle(c) = &*e.curve {
        let circle = c.circle();
        let axis = circle.frame().z().vector();
        let w = p - circle.centre();
        let w = w - axis * w.dot(axis);
        let v = u - axis * u.dot(axis);
        let (square, half, rest) = (
            v.dot(v),
            w.dot(v),
            w.dot(w) - circle.radius() * circle.radius(),
        );
        let discriminant = half * half - square * rest;
        if square <= 0.0 || discriminant <= 0.0 {
            return Ok(None);
        }
        let q = -(half + discriminant.sqrt().copysign(half));
        if q == 0.0 {
            return Ok(None);
        }
        found.extend([p + u * (q / square), p + u * (rest / q)]);
    } else {
        let across = u.cross(plane.plane().normal().vector());
        let size = across.magnitude();
        if size <= 0.0 {
            return Ok(None);
        }
        let across = across / size;
        let side = |q: Point| (q - p).dot(across);
        let corners = e.bound.expanded(reach).corners();
        if corners.iter().all(|q| side(*q) > 0.0) || corners.iter().all(|q| side(*q) < 0.0) {
            return Ok(None);
        }
        let (lo, hi) = e.crange;
        let at = |t: f64| -> OgeomResult<f64> { Ok(side(e.curve.point_at(t, tol)?)) };
        let mut before = (lo, at(lo)?);
        for i in 1..=64 {
            let t = lo + (hi - lo) * f64::from(i) / 64.0;
            let now = (t, at(t)?);
            if (before.1 < 0.0) != (now.1 < 0.0) {
                let (mut a, mut b, mut fa) = (before.0, now.0, before.1);
                for _ in 0..64 {
                    let m = f64::midpoint(a, b);
                    if m <= a || m >= b {
                        break;
                    }
                    let fm = at(m)?;
                    if (fa < 0.0) == (fm < 0.0) {
                        (a, fa) = (m, fm);
                    } else {
                        b = m;
                    }
                }
                let q = e.curve.point_at(f64::midpoint(a, b), tol)?;
                found.push(p + u * (q - p).dot(u));
            }
            before = now;
        }
    }
    let mut kept: Vec<Point> = Vec::new();
    for q in found {
        let foot = ogeom_algo::project_on_curve(&e.curve, q, 64, tol)?;
        if foot.distance <= reach
            && foot.parameter >= e.crange.0 - tol.parametric()
            && foot.parameter <= e.crange.1 + tol.parametric()
        {
            kept.push(q);
        }
    }
    Ok(match kept[..] {
        [at] => Some(Crossed::Once(at, u)),
        [first, second] => Some(Crossed::Twice(first, second)),
        _ => None,
    })
}

/// Where two lines come nearest each other, on the first, with the first's
/// unit direction: `None` unless both are lines and they are not parallel.
pub(crate) fn lines_crossing(
    first: &Curve,
    second: &Curve,
    tol: Tolerances,
) -> OgeomResult<Option<(Point, ogeom_math::Vector)>> {
    if !matches!(first, Curve::Line(_)) || !matches!(second, Curve::Line(_)) {
        return Ok(None);
    }
    let (p, q) = (first.point_at(0.0, tol)?, second.point_at(0.0, tol)?);
    let (u, v) = (
        first.point_at(1.0, tol)? - p,
        second.point_at(1.0, tol)? - q,
    );
    let across = u.cross(v);
    let square = across.dot(across);
    if square <= 0.0 {
        return Ok(None);
    }
    let s = (q - p).cross(v).dot(across) / square;
    Ok(Some((p + u * s, u / u.magnitude())))
}

/// The boundary edges a tangential contact runs along, carried onto the
/// other face as contacts.
///
/// Where two faces touch along a curve, neither passes through the other
/// there and the curve splits nothing: it is set aside from the
/// arrangement. But where the curve runs along a boundary edge of one
/// face, the solid's boundary turns at that edge onto a neighbour, and the
/// other face can pass from outside the solid into it across the line it
/// touches: a drill tangent to a wall along the edge where a concave round
/// leaves the wall, the drill flatter than the round, cuts a sliver of
/// material from the round's side on, and its wall is lined from the
/// tangent line to where it leaves the solid. The edge is carried onto the
/// other face as a contact, so the face is split along it and its pieces
/// either side are classified on their own; where both lie the same way,
/// the split costs an edge both faces already share.
pub(crate) fn tangent_contacts(
    ga: &GSolid,
    gb: &GSolid,
    tangents: &[TangentRec],
    contacts: &[ContactRec],
    tol: Tolerances,
) -> OgeomResult<Vec<ContactRec>> {
    let reach = tol.confusion() * GRAZE_REACH;
    let weld = tol.confusion() * 1e2;
    let mut out: Vec<ContactRec> = Vec::new();
    for contact in tangents {
        let (lo, hi) = contact.curve.domain();
        let mut contact_bound = ogeom_math::Aabb::EMPTY;
        for i in 0..=16 {
            let t = lo + (hi - lo) * f64::from(i) / 16.0;
            contact_bound = contact_bound.with_point(contact.curve.point_at(t, tol)?);
        }
        for (solid, owner, target, target_from_a, target_face) in [
            (
                ga,
                &ga.faces[contact.face_a],
                &gb.faces[contact.face_b],
                false,
                contact.face_b,
            ),
            (
                gb,
                &gb.faces[contact.face_b],
                &ga.faces[contact.face_a],
                true,
                contact.face_a,
            ),
        ] {
            for e in &owner.edges {
                let held = |c: &ContactRec| {
                    c.node == e.node
                        && c.target_from_a == target_from_a
                        && c.target_face == target_face
                };
                let near = reach.max(e.tolerance);
                if e.other_side.is_some()
                    || !e.bound.expanded(near).intersects(&contact_bound)
                    || contacts.iter().chain(&out).any(held)
                {
                    continue;
                }
                let span = e
                    .curve
                    .point_at(e.crange.0, tol)?
                    .distance(e.curve.point_at(e.crange.1, tol)?);
                if span <= weld {
                    continue;
                }
                // The edge lies along the contact curve and on the other
                // surface, end to end.
                let mut lies = true;
                for i in 0..=16 {
                    let t = e.crange.0 + (e.crange.1 - e.crange.0) * f64::from(i) / 16.0;
                    let p = e.curve.point_at(t, tol)?;
                    if ogeom_algo::project_on_curve(&contact.curve, p, 64, tol)?.distance > near
                        || ogeom_algo::project_on_surface(&target.surface, p, 16, tol)?.distance
                            > near
                    {
                        lies = false;
                        break;
                    }
                }
                if !lies {
                    continue;
                }
                // The solid's boundary turns across the other surface at
                // the edge: the face across the edge leaves that surface to
                // the other side from this one. A seam, or an edge whose
                // two faces both stay outside (a drill tangent to a box's
                // top along the top's edge), splits nothing.
                let Some(across) = solid.faces.iter().find(|f| {
                    !core::ptr::eq(*f, owner) && f.edges.iter().any(|x| x.node == e.node)
                }) else {
                    continue;
                };
                let Some(across_edge) = across.edges.iter().find(|x| x.node == e.node) else {
                    continue;
                };
                let (Some(own_side), Some(other_side)) = (
                    side_off_edge(owner, e, &target.surface, near, tol)?,
                    side_off_edge(across, across_edge, &target.surface, near, tol)?,
                ) else {
                    continue;
                };
                if own_side.signum() == other_side.signum() {
                    continue;
                }
                let Ok((pcurve, _, true, ..)) = ogeom_algo::pcurve_fit::fit_projected_pcurve(
                    &e.curve,
                    e.crange,
                    &target.surface,
                    tol,
                ) else {
                    continue;
                };
                ogeom_core::diag!(
                    Channel::Wire,
                    "TANGENT edge {} over {:?} onto face {target_face} (from_a {target_from_a})",
                    e.node.index(),
                    e.crange
                );
                out.push(ContactRec {
                    curve: (*e.curve).clone(),
                    crange: e.crange,
                    pcurve,
                    prange: e.crange,
                    node: e.node,
                    tolerance: e.tolerance.max(owner.tolerance + target.tolerance),
                    target_from_a,
                    target_face,
                    bound: e.bound,
                    crossed: false,
                });
            }
        }
    }
    Ok(out)
}

/// Which side of `other` the face lies on just off its edge `e`: the
/// signed distance along `other`'s normal of a point of the face a short
/// way in from the middle of the edge, or `None` where no such point is
/// decided (the face is narrower than any step tried, or lies on `other`
/// there to within the confusion distance).
///
/// The step is taken in the face's chart, square to the edge's pcurve,
/// starting at a twentieth of the edge's span and halving until the point
/// lands inside the face's trim.
pub(crate) fn side_off_edge(
    face: &GFace,
    e: &BoundaryEdge,
    other: &SurfaceGeometry,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let (p0, p1) = e.prange;
    let tm = f64::midpoint(p0, p1);
    let h = (p1 - p0).abs() * 1e-3;
    if h <= 0.0 {
        return Ok(None);
    }
    let at = e.pcurve.point_at(tm, tol)?;
    let (before, after) = (
        e.pcurve.point_at(tm - h, tol)?,
        e.pcurve.point_at(tm + h, tol)?,
    );
    let along = (after.x - before.x, after.y - before.y);
    let length = along.0.hypot(along.1);
    if length <= f64::MIN_POSITIVE {
        return Ok(None);
    }
    let square = (-along.1 / length, along.0 / length);
    let (su, sv) = face.surface.d1_at(at.x, at.y, tol)?;
    let per_param = su.magnitude() * square.0.abs() + sv.magnitude() * square.1.abs();
    if per_param <= f64::MIN_POSITIVE {
        return Ok(None);
    }
    let span = e
        .curve
        .point_at(e.crange.0, tol)?
        .distance(e.curve.point_at(e.crange.1, tol)?);
    let mut step = span * 0.05;
    for _ in 0..8 {
        for sign in [1.0, -1.0] {
            let k = sign * step / per_param;
            let q2 = Point2::new(square.0.mul_add(k, at.x), square.1.mul_add(k, at.y));
            let Ok(q) = face.surface.point_at(q2.x, q2.y, tol) else {
                continue;
            };
            if chart_point_within(face, q, reach, tol).is_none() {
                continue;
            }
            let foot = ogeom_geom::project_on_surface(other, q, 16, tol)?;
            let normal = other.normal_at(foot.parameters.0, foot.parameters.1, tol)?;
            let side = (q - foot.point).dot(normal.vector());
            if side.abs() > tol.confusion() {
                return Ok(Some(side));
            }
            return Ok(None);
        }
        step *= 0.5;
    }
    Ok(None)
}

/// One split an edge is asked for, with how honestly it can be placed.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Pave {
    /// The parameter on the edge's curve.
    pub(crate) t: f64,
    /// How far the strand that asked for it may honestly sit from the exact
    /// junction: its own tolerance, and its tangential doubt where the
    /// crossing was a touch.
    pub(crate) honesty: f64,
    /// The face of the other solid whose section crossed the edge here
    /// (whether that face is of `a`, and its index), where the pave is
    /// such a crossing.
    pub(crate) across: Option<(bool, usize)>,
    /// How far along the edge the crossing may stand from the one the
    /// face across the edge finds: a loose edge's tolerance over the sine
    /// of its lean out of the face that crosses it. Paves that near are
    /// one junction, which reaches only as far as they stand apart.
    pub(crate) spread: f64,
}

/// How far along a loose edge its crossing by `surface` at `at` may
/// stand from the crossing the face across the edge finds: the edge's
/// tolerance over the sine of its lean out of the surface there, held to
/// [`LEAN_SPREAD`] tolerances. None for a tight edge, and none where the
/// lean cannot be read (a cone's apex, a curve with no tangent there).
pub(crate) fn lean_spread(
    e: &BoundaryEdge,
    t: f64,
    surface: &SurfaceGeometry,
    at: Point,
    tol: Tolerances,
) -> f64 {
    if e.tolerance <= tol.confusion() * DRAWN_APART {
        return 0.0;
    }
    let Ok(along) = e.curve.d1_at(t, tol) else {
        return 0.0;
    };
    let length = along.magnitude();
    let Ok(foot) = ogeom_geom::project_on_surface(surface, at, 16, tol) else {
        return 0.0;
    };
    let Ok(normal) = surface.normal_at(foot.parameters.0, foot.parameters.1, tol) else {
        return 0.0;
    };
    if length <= 0.0 {
        return 0.0;
    }
    let sine = along.dot(normal.vector()).abs() / length;
    let most = e.tolerance * LEAN_SPREAD;
    if sine * most > e.tolerance {
        e.tolerance / sine
    } else {
        most
    }
}

/// How far along an edge its crossing by a fitted section at `at` may
/// stand from the true one where the edge touches the section's other
/// surface there: zero where it crosses that surface at an angle.
///
/// The section lies within `honesty` of the surface, so it meets the edge
/// anywhere along the stretch where the edge does. Along the edge, its
/// signed distance off the surface runs `g(s) = g1 s + g2 s^2 / 2` from the
/// crossing: `g1` the edge's lean out of the surface's tangent plane, `g2`
/// the edge's curvature towards the surface's normal less the surface's own
/// normal curvature along the edge. Where `g1^2 <= 2 |g2| honesty` the
/// valley's floor lies within `honesty` of the surface: the edge touches
/// the surface rather than crossing it twice, and the stretch runs to the
/// far side of the valley, `(|g1| + sqrt(g1^2 + 2 |g2| honesty)) / |g2|` from
/// the crossing. Two faces whose surfaces are tangent where a section
/// leaves them (two rounds of one top face whose rails cross) meet it so:
/// the face across splits its own edge at the rails' crossing, as much as
/// a hundredth of a millimetre from where the section meets this one.
/// Elsewhere the lean bounds the stretch at `honesty / |g1|`, which the
/// crossing's own honesty already covers.
pub(crate) fn touch_spread(
    (curve, t): (&Curve, f64),
    surface: &SurfaceGeometry,
    at: Point,
    honesty: f64,
    tol: Tolerances,
) -> f64 {
    if honesty <= 0.0 {
        return 0.0;
    }
    let Ok(d) = curve.derivatives_at(t, 2, tol) else {
        return 0.0;
    };
    let (Some(&d1), Some(&d2)) = (d.get(1), d.get(2)) else {
        return 0.0;
    };
    let speed = d1.magnitude();
    let Ok(foot) = ogeom_geom::project_on_surface(surface, at, 16, tol) else {
        return 0.0;
    };
    let (u, v) = foot.parameters;
    let (Ok((_, su, sv)), Ok((suu, suv, svv))) =
        (surface.point_d1_at(u, v, tol), surface.d2_at(u, v, tol))
    else {
        return 0.0;
    };
    let normal = su.cross(sv);
    let area = normal.magnitude();
    if speed <= 0.0 || area <= 0.0 {
        return 0.0;
    }
    let normal = normal * (1.0 / area);
    let tangent = d1 * (1.0 / speed);
    // The edge's unit tangent in the surface's parameters, solved on the
    // first fundamental form (whose determinant is `area^2`), and the
    // surface's normal curvature along it.
    let (e, f, g) = (su.dot(su), su.dot(sv), sv.dot(sv));
    let (tu, tv) = (tangent.dot(su), tangent.dot(sv));
    let du = g.mul_add(tu, -f * tv) / (area * area);
    let dv = e.mul_add(tv, -f * tu) / (area * area);
    let bends = (du * du).mul_add(
        suu.dot(normal),
        (2.0 * du * dv).mul_add(suv.dot(normal), dv * dv * svv.dot(normal)),
    );
    let lean = tangent.dot(normal).abs();
    let curving = (d2 - tangent * d2.dot(tangent)) * (1.0 / (speed * speed));
    let parting = (curving.dot(normal) - bends).abs();
    if parting <= 0.0 || lean * lean > 2.0 * parting * honesty {
        return 0.0;
    }
    (lean + lean.mul_add(lean, 2.0 * parting * honesty).sqrt()) / parting
}

/// Paves the edge cannot tell apart, as one junction each.
///
/// A tolerant rail meeting a wedge's faces near one corner collects a
/// cluster of crossings inside its own stated radius, and a fitted section
/// and a fitted contact asking for the same corner land a few microns
/// apart. Split at each, the edge shatters into dust no weld downstream can
/// rejoin. So consecutive paves whose gap along the edge in space is within
/// the edge's honesty or either pave's own are one cluster, its first pave
/// speaking for it, and so are the two ends of a graze (see [`grazed`]),
/// the graze's nearest approach speaking for them. A pave known to within
/// the edge's own floor (a crossing found in closed form) speaks for its
/// cluster over looser ones: where a marched section ends a few microns off
/// a point two edges cross at, the edges split at that crossing, not at the
/// march's end. A looser pave never takes over from the first.
#[derive(Debug, Clone, Copy)]
pub(crate) struct PaveCluster {
    /// The representative parameter: the cluster's first pave, or its most
    /// honest one known to within the edge's floor, or the nearest approach
    /// of a graze it holds.
    pub(crate) t: f64,
    /// Where it sits.
    pub(crate) at: Point,
    /// How far the cluster's members reach from the representative.
    pub(crate) span: f64,
    /// The loosest honesty among the members.
    pub(crate) honesty: f64,
    /// How many paves the cluster holds.
    pub(crate) members: usize,
}

pub(crate) fn cluster_paves(
    curve: &Curve,
    crange: (f64, f64),
    edge_tolerance: f64,
    paves: &[Pave],
    off: Option<OffFace<'_>>,
    along: Along<'_>,
    tol: Tolerances,
) -> OgeomResult<Vec<PaveCluster>> {
    let mut ts: Vec<Pave> = paves
        .iter()
        .copied()
        .filter(|p| p.t > crange.0 + tol.parametric() && p.t < crange.1 - tol.parametric())
        .collect();
    ts.sort_by(|x, y| x.t.partial_cmp(&y.t).unwrap_or(core::cmp::Ordering::Equal));
    ts.dedup_by(|x, y| {
        if (x.t - y.t).abs() <= tol.parametric() {
            y.honesty = y.honesty.max(x.honesty);
            true
        } else {
            false
        }
    });
    let floor = edge_tolerance.max(tol.confusion() * 10.0);
    // A pave within the edge's honesty of its own end *is* the end: split
    // there, the sliver between them is dust one face keeps and the
    // face across the edge drops, and the sew finds it used once.
    let ends = [
        curve.point_at(crange.0, tol)?,
        curve.point_at(crange.1, tol)?,
    ];
    let mut clusters: Vec<PaveCluster> = Vec::new();
    let mut prev: Option<(Point, Pave)> = None;
    // Where the current cluster's members stand.
    let mut held_at: Vec<Point> = Vec::new();
    // The current cluster's representative: the honesty it was chosen for,
    // or none once a graze's nearest approach speaks for the cluster.
    let mut speaking: Option<f64> = None;
    for pave in ts {
        let at = curve.point_at(pave.t, tol)?;
        if ends
            .iter()
            .any(|e| e.distance(at) <= floor.max(pave.honesty))
        {
            continue;
        }
        // Two crossings each known to within its honesty are one point
        // where they stand no further apart than both together. Two
        // crossings of a loose edge within their spreads are one junction
        // too, while every member stands within those spreads of the new
        // one: a series of shallow crossings, each within its spread of the
        // last, does not chain into a junction longer than the spreads.
        let joined = prev.is_some_and(|(held, last): (Point, Pave)| {
            let apart = held.distance(at);
            apart <= floor.max(last.honesty + pave.honesty)
                || (apart <= last.spread + pave.spread
                    && held_at
                        .iter()
                        .all(|p| p.distance(at) <= last.spread + pave.spread))
        });
        let touch = match (prev, off) {
            (Some((held, last)), Some(off)) if !joined => grazed(
                curve,
                (held, last),
                (at, pave),
                edge_tolerance,
                off,
                along,
                tol,
            )?,
            _ => None,
        };
        if (joined || touch.is_some())
            && let Some(cluster) = clusters.last_mut()
        {
            held_at.push(at);
            if let Some((t, point)) = touch {
                cluster.t = t;
                cluster.at = point;
                speaking = None;
            } else if pave.honesty <= floor && speaking.is_some_and(|best| pave.honesty < best) {
                cluster.t = pave.t;
                cluster.at = at;
                speaking = Some(pave.honesty);
            }
            cluster.span = held_at
                .iter()
                .map(|p| p.distance(cluster.at))
                .fold(0.0, f64::max);
            cluster.honesty = cluster.honesty.max(pave.honesty);
            cluster.members += 1;
        } else {
            held_at.clear();
            held_at.push(at);
            speaking = Some(pave.honesty);
            clusters.push(PaveCluster {
                t: pave.t,
                at,
                span: 0.0,
                honesty: pave.honesty,
                members: 1,
            });
        }
        prev = Some((at, pave));
    }
    Ok(clusters)
}

/// How far a point stands off a face of the other solid, the face named as
/// [`Pave::across`] names it.
pub(crate) type OffFace<'f> = &'f dyn Fn((bool, usize), Point) -> OgeomResult<f64>;

/// A stretch of an edge carried onto a face of the other solid from where
/// section lines cross it ([`ContactRec::crossed`]): its range on the edge,
/// and that face as [`Pave::across`] names it.
pub(crate) type Stretch = ((f64, f64), (bool, usize));

/// The [`Stretch`]es of one edge.
pub(crate) type Along<'h> = &'h [Stretch];

/// The [`Stretch`]es of every edge that has any.
pub(crate) type AlongEdges = ogeom_core::FastMap<EdgeKey, Vec<Stretch>>;

/// [`AlongEdges`] for the contacts carried.
pub(crate) fn crossed_along(contacts: &[ContactRec]) -> AlongEdges {
    let mut along = AlongEdges::default();
    for c in contacts.iter().filter(|c| c.crossed) {
        along
            .entry(c.node)
            .or_default()
            .push((c.crange, (c.target_from_a, c.target_face)));
    }
    along
}

/// How far `p` stands off the surface of a face of `a` (the flag set) or
/// of `b`.
pub(crate) fn off_face(
    ga: &GSolid,
    gb: &GSolid,
    (from_a, fi): (bool, usize),
    p: Point,
    tol: Tolerances,
) -> OgeomResult<f64> {
    let face = if from_a { &ga.faces[fi] } else { &gb.faces[fi] };
    Ok(ogeom_algo::project_on_surface(&face.surface, p, 16, tol)?.distance)
}

/// Where two consecutive crossings of an edge by sections on one face of
/// the other solid are the two ends of a graze: the edge between them never
/// leaving that face by more than the crossings' honesty or its own
/// tolerance. An exact edge dipping deeper into the face than its
/// tolerance crosses it twice, as the surfaces either side of it do, and
/// the sections ending there end apart. A face tangent
/// to the other solid's face along a line meets a third face of its own
/// solid a hair off the tangency (a round beside the face, its tangent edge
/// a few microns from the other solid). The sections either side of the
/// touch then cross the round's edge twice, the root of the hair times the
/// curvature apart, on an edge the other solid cannot be told from between
/// them. Split twice, the edge keeps a stub no face across it has. The two
/// are one junction at the edge's nearest approach to the face, its
/// parameter and point, while the stretch stays short: a graze and not an
/// edge running along the face. The honesty asked of the edge is held to
/// the fit-slop bound [`honest`] keeps (a micron at millimetre
/// tolerances), however loosely the sections state themselves. The two
/// crossings may stand at most a hundred times that bound apart, whatever
/// the sections state: how far apart a graze's crossings land is set by
/// the gap and the curvature at the touch, not by how closely the
/// sections were fitted, and a window read off their fit error would
/// split a graze or not as a fit happens to land. So the stub a junction
/// takes in is at most 0.1 long at millimetre tolerances, and the
/// junction's reach never grows with a section's stated error. A stretch
/// `along` carries onto the face between the two is an edge running along
/// it, and no graze.
pub(crate) fn grazed(
    curve: &Curve,
    (from, first): (Point, Pave),
    (to, second): (Point, Pave),
    edge_tolerance: f64,
    off: OffFace<'_>,
    along: Along<'_>,
    tol: Tolerances,
) -> OgeomResult<Option<(f64, Point)>> {
    let (Some(face), Some(other)) = (first.across, second.across) else {
        return Ok(None);
    };
    let honesty = honest(first.honesty.max(second.honesty).max(edge_tolerance), tol);
    if face != other || from.distance(to) > honest(f64::INFINITY, tol) * 1e2 {
        return Ok(None);
    }
    let (lo, hi) = (first.t.min(second.t), first.t.max(second.t));
    if along.iter().any(|&((a, b), onto)| {
        onto == face && a <= lo + tol.parametric() && b >= hi - tol.parametric()
    }) {
        return Ok(None);
    }
    let mut nearest: Option<(f64, f64, Point)> = None;
    for k in 1..8 {
        let t = (second.t - first.t).mul_add(f64::from(k) / 8.0, first.t);
        let p = curve.point_at(t, tol)?;
        let d = off(face, p)?;
        if d > honesty {
            return Ok(None);
        }
        if nearest.is_none_or(|(best, ..)| d < best) {
            nearest = Some((d, t, p));
        }
    }
    Ok(nearest.map(|(_, t, p)| (t, p)))
}

/// One kept sub-range of one section.
#[derive(Clone)]
pub(crate) struct SectionPiece {
    pub(crate) section: usize,
    pub(crate) range: (f64, f64),
    /// Whether this sub-range already *is* boundary on the A face, the B
    /// face, or neither. A piece that is boundary on one side is still the
    /// splitting curve on the other, and only the face it duplicates leaves
    /// it out.
    pub(crate) hugs: [bool; 2],
    /// Set when the piece was admitted by a hug: a split of the other face
    /// along this edge of the hugging face: (edge node, target from A,
    /// target face). Two sections hugging one edge onto one face are the
    /// same split, kept once.
    pub(crate) hug_key: Option<(EdgeKey, bool, usize)>,
}

/// One boundary edge of one argument's face, lying in a face of the other
/// argument's own surface: the splitting curve same-domain contact
/// contributes, since coincident surfaces have no section curve to offer.
pub(crate) struct ContactRec {
    /// The owner's world curve and the sub-range its edge covers.
    pub(crate) curve: Curve,
    pub(crate) crange: (f64, f64),
    /// The curve spoken in the *target* face's chart, over its own window,
    /// mapped proportionally from `crange` exactly as an edge's own stored
    /// pcurve is.
    pub(crate) pcurve: PlanarCurve,
    pub(crate) prange: (f64, f64),
    /// The owner edge's node, whose paves this record shares.
    pub(crate) node: EdgeKey,
    /// The owner edge's own tolerance: how far its curve may honestly sit
    /// from the exact geometry it meets, which is how far a crossing filter
    /// must reach to see a fitted rail cross an exact seam.
    pub(crate) tolerance: f64,
    /// The target: which argument's face list, and which face.
    pub(crate) target_from_a: bool,
    pub(crate) target_face: usize,
    /// The owner edge's conservative box: a target edge whose own box
    /// stands further off than a crossing can reach is not intersected.
    pub(crate) bound: ogeom_math::Aabb,
    /// Whether the stretch runs between two crossings of the edge by a
    /// section line, or from one to a vertex, where the edge is a curve
    /// (see [`graze_contacts`]): the edge runs along the target face there,
    /// however short the stretch.
    pub(crate) crossed: bool,
}

/// A same-domain contact edge's image in the shared surface's chart, where
/// no closed form exists: fitted by projection into the target's chart, or
/// (where the edge overhangs the target's own window, a blend's leg on a
/// spline host continued past the face) into the owner's, which is the
/// same chart wherever an extension kept its parent's parameters. The
/// owner's image is trusted only where lifting it through the *target*
/// lands back on the edge, sampled over the run the target's window holds.
///
/// `None` for an edge lying wholly outside the target's window, which the
/// target face cannot meet.
///
/// # Errors
///
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if neither
/// chart holds the edge, or the owner's image lies off the target's window
/// while the edge lies on it.
pub(crate) fn projected_into_shared_chart(
    curve: &Curve,
    range: (f64, f64),
    owner: &SurfaceGeometry,
    target: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<Option<PlanarCurve>> {
    let against_target = ogeom_algo::pcurve_fit::fit_projected_pcurve(curve, range, target, tol);
    if let Ok(fitted) = &against_target
        && fitted.2
    {
        return Ok(Some(fitted.0.clone()));
    }
    let Ok(on_owner) = ogeom_algo::pcurve_fit::fit_projected_pcurve(curve, range, owner, tol)
    else {
        ogeom_bail!(
            NotDone,
            "same-domain contact whose edge could not be projected into the \
             shared surface's chart: {}",
            against_target
                .err()
                .map_or_else(String::new, |e| e.to_string())
        );
    };
    let ((ua, ub), (va, vb)) = target.domain();
    let reach = tol.confusion() * 1e4;
    let at = |i: u32| range.0 + (range.1 - range.0) * f64::from(i) / 16.0;
    let mut checked = 0usize;
    for i in 0..=16 {
        let t = at(i);
        let uv = on_owner.0.point_at(t, tol)?;
        if uv.x < ua || uv.x > ub || uv.y < va || uv.y > vb {
            continue;
        }
        let lifted = target.point_at(uv.x, uv.y, tol)?;
        if lifted.distance(curve.point_at(t, tol)?) > reach {
            ogeom_bail!(
                NotDone,
                "same-domain contact whose edge overhangs the shared surface's \
                 window, on a chart the two surfaces do not share"
            );
        }
        checked += 1;
    }
    if checked == 0 {
        // The owner's image clear of the target's window says the edge is
        // clear of it only on a shared chart. Where the charts differ (the
        // axis reversed, so the owner's heights read negated), the edge's
        // own feet on the target say whether it reaches the window.
        for i in 0..=16 {
            if on_window(curve.point_at(at(i), tol)?, target, reach, tol)? {
                ogeom_bail!(
                    NotDone,
                    "same-domain contact whose edge lies on the shared surface's \
                     window, on a chart the two surfaces do not share"
                );
            }
        }
        return Ok(None);
    }
    Ok(Some(on_owner.0))
}

/// Whether a point's foot on a surface lies within `reach` of it and clear
/// of the rim of its window. A foot the search pinned to the rim is past
/// the window, not on it.
pub(crate) fn on_window(
    p: Point,
    surface: &SurfaceGeometry,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let ((ua, ub), (va, vb)) = surface.domain();
    let (mu, mv) = ((ub - ua).abs() * 1e-9, (vb - va).abs() * 1e-9);
    let Ok(foot) = ogeom_algo::project_on_surface(surface, p, 16, tol) else {
        return Ok(false);
    };
    let (u, v) = foot.parameters;
    Ok(foot.distance <= reach
        && (!u.is_finite() || !ua.is_finite() || (u > ua + mu && u < ub - mu))
        && (!v.is_finite() || !va.is_finite() || (v > va + mv && v < vb - mv)))
}

/// The stretches of an edge that lie on a coincident surface's window, and
/// their images in that surface's chart, where the edge as a whole cannot
/// be imaged there: it runs on past the window (a converted part's wall
/// edge against a pad's shorter side swept from the part's own outline),
/// and only what lies within it can split the face on that surface.
///
/// The edge is sampled, and the samples whose foot on the surface is
/// within `reach` and clear of the window's rim make runs, one for each
/// time the edge comes onto the window. Each run's ends are narrowed onto
/// the rim by bisection and its image fitted over the run alone, so an
/// edge that leaves the window and comes back has one image per stretch.
/// Empty where no sample lands on the window. An error where a stretch's
/// fit fails.
pub(crate) fn clipped_into_chart(
    curve: &Curve,
    range: (f64, f64),
    target: &SurfaceGeometry,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<(PlanarCurve, (f64, f64))>> {
    const SAMPLES: u32 = 64;
    let on =
        |t: f64| -> OgeomResult<bool> { on_window(curve.point_at(t, tol)?, target, reach, tol) };
    let at = |i: u32| range.0 + (range.1 - range.0) * f64::from(i) / f64::from(SAMPLES);
    let mut inside = Vec::with_capacity(SAMPLES as usize + 1);
    for i in 0..=SAMPLES {
        inside.push(on(at(i))?);
    }
    // Each end narrowed between the last sample on the window and the
    // first past it.
    let narrow = |mut good: f64, mut bad: f64| -> OgeomResult<f64> {
        for _ in 0..40 {
            let mid = f64::midpoint(good, bad);
            if on(mid)? {
                good = mid;
            } else {
                bad = mid;
            }
        }
        Ok(good)
    };
    let mut stretches = Vec::new();
    let mut next = 0usize;
    while let Some(offset) = inside[next..].iter().position(|&x| x) {
        let first = next + offset;
        let last = inside[first..]
            .iter()
            .position(|&x| !x)
            .map_or(inside.len() - 1, |gap| first + gap - 1);
        next = last + 1;
        #[allow(clippy::cast_possible_truncation, reason = "sample indices")]
        let lo = if first == 0 {
            range.0
        } else {
            narrow(at(first as u32), at(first as u32 - 1))?
        };
        #[allow(clippy::cast_possible_truncation, reason = "sample indices")]
        let hi = if last == SAMPLES as usize {
            range.1
        } else {
            narrow(at(last as u32), at(last as u32 + 1))?
        };
        if hi - lo <= tol.parametric() {
            continue;
        }
        let fitted = ogeom_algo::pcurve_fit::fit_projected_pcurve(curve, (lo, hi), target, tol)?;
        if !fitted.2 {
            ogeom_bail!(
                NotDone,
                "same-domain contact whose edge's stretch on the shared surface's \
                 window has no image in its chart"
            );
        }
        stretches.push((fitted.0, (lo, hi)));
    }
    Ok(stretches)
}

/// The map from one cylinder's, cone's or torus's chart to another's where
/// both are one surface on two frames: the same axis line, within `slop`,
/// and the same radii, and for a torus the same centre on it. The frames
/// differ by a turn about the axis, a slide along it (a cylinder's or a
/// cone's), the axis reversed, or the handedness, so the map is exact:
/// `(u, v)` goes to `(σ·u + φ, s·v + c)`, `s` the sign between the axes,
/// `σ` the sign between the axes the two angles turn about (`x × y`, which
/// is `z` signed by the frame's handedness), `φ` the owner's angle zero read
/// in the target's frame, and `c` the owner's origin's height on the
/// target's axis, zero on a torus. On a torus `v` is the tube's angle,
/// turned over with the axis. A cone's `v` is a height along its axis too,
/// its radius `r + v tan(a)` on either frame: one cone where the half
/// angles agree up to the axes' sign and the owner's reference radius is
/// the target's radius at the owner's origin. Where `σ` and `s` agree the
/// map is a rigid motion of the plane, a half turn where both are negative;
/// where they differ it is a mirror.
pub(crate) fn coaxial_chart_map(
    owner: &SurfaceGeometry,
    target: &SurfaceGeometry,
    slop: f64,
    tol: Tolerances,
) -> Option<ogeom_math::Transform2> {
    let reach = slop.max(tol.confusion());
    let mut cones = None;
    let (fa, fb) = match (owner, target) {
        (SurfaceGeometry::Cylinder(a), SurfaceGeometry::Cylinder(b)) => {
            let (a, b) = (a.cylinder(), b.cylinder());
            if (a.radius() - b.radius()).abs() > reach {
                return None;
            }
            (a.frame(), b.frame())
        }
        (SurfaceGeometry::Cone(a), SurfaceGeometry::Cone(b)) => {
            let (a, b) = (a.cone(), b.cone());
            cones = Some((a, b));
            (a.frame(), b.frame())
        }
        (SurfaceGeometry::Torus(a), SurfaceGeometry::Torus(b)) => {
            let (a, b) = (a.torus(), b.torus());
            if (a.major_radius() - b.major_radius()).abs() > reach
                || (a.minor_radius() - b.minor_radius()).abs() > reach
                || a.frame().origin().distance(b.frame().origin()) > reach
            {
                return None;
            }
            (a.frame(), b.frame())
        }
        _ => return None,
    };
    let (za, zb) = (fa.z().vector(), fb.z().vector());
    let axis = ogeom_math::Axis::new(fb.origin(), fb.z());
    // The axes are one line where the owner's stands within `reach` of the
    // target's wherever the owner's face can be: at its origin, and tilted
    // by their angle over the farthest the face reaches from it, which is
    // the ends of a cylinder's or a cone's height window and a torus's
    // outer radius. An unbounded window is held to the angle alone.
    let lever = match owner {
        SurfaceGeometry::Torus(t) => t.torus().major_radius() + t.torus().minor_radius(),
        _ => {
            let (_, (va, vb)) = owner.domain();
            va.abs().max(vb.abs())
        }
    };
    let tilt = za.cross(zb).magnitude();
    let tilted = if lever.is_finite() {
        tilt * lever > reach
    } else {
        tilt > tol.angular()
    };
    if tilted || axis.distance_to(fa.origin()) > reach {
        return None;
    }
    let s = if za.dot(zb) >= 0.0 { 1.0 } else { -1.0 };
    let sigma = s * fa.handedness().sign() * fb.handedness().sign();
    let x = fa.x().vector();
    let phi = x.dot(fb.y().vector()).atan2(x.dot(fb.x().vector()));
    let c = if matches!(target, SurfaceGeometry::Torus(_)) {
        0.0
    } else {
        (fa.origin() - fb.origin()).dot(zb)
    };
    if let Some((a, b)) = cones
        && ((a.half_angle() - s * b.half_angle()).abs() > tol.angular()
            || (a.reference_radius() - b.radius_at(c)).abs() > reach)
    {
        return None;
    }
    let slide = ogeom_math::Transform2::translation(ogeom_math::Vector2::new(phi, c));
    Some(match (sigma > 0.0, s > 0.0) {
        (true, true) => slide,
        (false, false) => {
            ogeom_math::Transform2::rotation(Point2::new(phi / 2.0, c / 2.0), core::f64::consts::PI)
        }
        // `u` to `φ - u`: the mirror in the column `u = φ / 2`, slid by `c`.
        (false, true) => {
            ogeom_math::Transform2::translation(ogeom_math::Vector2::new(0.0, c))
                * ogeom_math::Transform2::line_mirror(
                    Point2::new(phi / 2.0, 0.0),
                    ogeom_math::Direction2::X,
                )
        }
        // `v` to `c - v`: the mirror in the row `v = c / 2`, slid by `φ`.
        (true, false) => {
            ogeom_math::Transform2::translation(ogeom_math::Vector2::new(phi, 0.0))
                * ogeom_math::Transform2::line_mirror(
                    Point2::new(0.0, c / 2.0),
                    ogeom_math::Direction2::Y,
                )
        }
    })
}

/// An owner face's pcurve carried into the target's chart by
/// [`coaxial_chart_map`], on the branch of the target's `u` window its middle
/// falls in.
pub(crate) fn carried_across(
    pcurve: &PlanarCurve,
    prange: (f64, f64),
    owner: &SurfaceGeometry,
    target: &SurfaceGeometry,
    slop: f64,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    let map = coaxial_chart_map(owner, target, slop, tol)?;
    let carried = pcurve.transformed(&map, tol).ok()?;
    let middle = carried.point_at(0.5 * (prange.0 + prange.1), tol).ok()?;
    let ((ua, ub), (va, vb)) = target.domain();
    let turn = core::f64::consts::TAU;
    let turns = |at: f64, lo: f64, hi: f64| -> f64 {
        if at < lo {
            ((lo - at) / turn).ceil()
        } else if at > hi {
            -((at - hi) / turn).ceil()
        } else {
            0.0
        }
    };
    // A torus's tube angle is periodic too: carried over with a reversed
    // axis it may land a turn off the target's `v` window.
    let along_v = if target.is_periodic_v() {
        turns(middle.y, va, vb)
    } else {
        0.0
    };
    let shift = ogeom_math::Vector2::new(turns(middle.x, ua, ub) * turn, along_v * turn);
    if shift.x == 0.0 && shift.y == 0.0 {
        return Some(carried);
    }
    carried
        .transformed(&ogeom_math::Transform2::translation(shift), tol)
        .ok()
}

/// Whether two surfaces are the *identical chart*: the same
/// parameterization, frame and all, not merely the same point set.
///
/// The five analytics only: a fitted surface never qualifies, because two
/// independent fits of one geometry agree nowhere in parameter space. This
/// is what lets a stored pcurve stand in for a closed-form projection in the
/// same-domain melt: on the identical chart, the owner's pcurve already is
/// the projection.
pub(crate) fn same_chart(a: &SurfaceGeometry, b: &SurfaceGeometry, tol: Tolerances) -> bool {
    use SurfaceGeometry as S;
    // The handedness too: on frames sharing `x` and `z`, a left-handed
    // chart's `y` is the right-handed one's reversed, and its angle runs
    // the other way round.
    let frames = |fa: ogeom_math::Frame, fb: ogeom_math::Frame| -> bool {
        fa.handedness() == fb.handedness()
            && fa.origin().distance(fb.origin()) <= tol.confusion()
            && fa.z().vector().dot(fb.z().vector()) >= 1.0 - tol.angular()
            && fa.x().vector().dot(fb.x().vector()) >= 1.0 - tol.angular()
    };
    match (a, b) {
        (S::Plane(x), S::Plane(y)) => frames(x.plane().frame(), y.plane().frame()),
        (S::Cylinder(x), S::Cylinder(y)) => {
            frames(x.cylinder().frame(), y.cylinder().frame())
                && (x.cylinder().radius() - y.cylinder().radius()).abs() <= tol.confusion()
        }
        (S::Cone(x), S::Cone(y)) => {
            frames(x.cone().frame(), y.cone().frame())
                && (x.cone().reference_radius() - y.cone().reference_radius()).abs()
                    <= tol.confusion()
                && (x.cone().half_angle() - y.cone().half_angle()).abs() <= tol.angular()
        }
        (S::Sphere(x), S::Sphere(y)) => {
            frames(x.sphere().frame(), y.sphere().frame())
                && (x.sphere().radius() - y.sphere().radius()).abs() <= tol.confusion()
        }
        (S::Torus(x), S::Torus(y)) => {
            frames(x.torus().frame(), y.torus().frame())
                && (x.torus().major_radius() - y.torus().major_radius()).abs() <= tol.confusion()
                && (x.torus().minor_radius() - y.torus().minor_radius()).abs() <= tol.confusion()
        }
        // Two patches are one chart only as one patch: the same knots and
        // the same net, to the bit: a blend's leg built on the host's own
        // patch, widened in place and shared. Any other pair of patches
        // that merely coincides is not.
        (S::BSpline(x), S::BSpline(y)) => x == y,
        _ => false,
    }
}

/// What a face's arrangement strand stands for.
#[derive(Clone)]
pub(crate) enum Tag {
    /// A sub-range of boundary edge `edge` (index into the face's edges).
    Boundary { edge: usize, range: (f64, f64) },
    /// A sub-range of a global section.
    Section { section: usize, range: (f64, f64) },
    /// A sub-range of a global contact edge: another face's boundary edge
    /// lying in this face's own surface, splitting it.
    Contact { contact: usize, range: (f64, f64) },
    /// A sub-range of a pole of the face's own chart. The edge is a point in
    /// space whatever the range says. The range is where it runs in the
    /// chart, which is the only place it has length.
    Pole { pole: usize, range: (f64, f64) },
}

/// Where a piece stands relative to the other solid.
///
/// `In` and `Out` are the classifier's words. The two `On` states split the
/// case the classifier cannot decide alone: a piece lying on the other
/// boundary bounds material on one side here and one side there, and whether
/// those sides agree (outward normals aligned) or oppose is what every
/// operation's filter turns on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PieceState {
    In,
    Out,
    /// On the other boundary, outward normals aligned: both solids' material
    /// on the same side. One copy of the piece bounds the union and the
    /// intersection alike.
    OnAligned,
    /// On the other boundary, outward normals opposed: material on both
    /// sides. The contact is interior to the union and vanishes from it.
    OnOpposed,
    /// Not classified: the other argument is a sheet, which bounds nothing.
    Unread,
}

/// A point in a polyline's *interior*: its half-way member, except that a
/// straight strand is two points and its half-way member is an endpoint,
/// which may sit exactly on a seam or a trim. The chord midpoint is on the
/// curve for a straight image and strictly inside either way.
pub(crate) fn interior_of(line: &[Point2]) -> Point2 {
    if line.len() == 2 {
        Point2::new(
            f64::midpoint(line[0].x, line[1].x),
            f64::midpoint(line[0].y, line[1].y),
        )
    } else {
        line[line.len() / 2]
    }
}

/// Fold a param-space polyline into a periodic surface's chart, by one
/// constant offset per axis.
///
/// A section's pcurve is *unwrapped* across a seam: continuous, and allowed
/// to leave the stated domain, because that is what crossing a seam is. The
/// arrangement lives in one chart, and the filler has already split every
/// section at its seam crossings, so each strand spans at most one chart
/// width and a single period shift per axis brings it home. The shift is
/// chosen by the strand's midpoint, so endpoints sitting exactly on the
/// chart's edge stay on whichever side the strand's body is.
pub(crate) fn fold_into_chart(line: &mut [Point2], surface: &SurfaceGeometry, tol: Tolerances) {
    let ((ua, _), (va, _)) = surface.domain();
    if line.is_empty() {
        return;
    }
    let mid = interior_of(line);
    let [span_u, span_v] = chart_spans(surface, tol);
    if let Some(span) = span_u {
        let shift = (ua + (mid.x - ua).rem_euclid(span)) - mid.x;
        for p in line.iter_mut() {
            p.x += shift;
        }
    }
    if let Some(span) = span_v {
        let shift = (va + (mid.y - va).rem_euclid(span)) - mid.y;
        for p in line.iter_mut() {
            p.y += shift;
        }
    }
}

/// The width a chart comes round by along each direction: its period, or
/// its span where it closes on itself without repeating (a drum spelt as a
/// clamped spline). A section's pcurve is unwrapped across either alike,
/// so a point of it lies a whole span off the chart's window wherever the
/// section crosses the closure, and folds back by that span.
pub(crate) fn chart_spans(surface: &SurfaceGeometry, tol: Tolerances) -> [Option<f64>; 2] {
    let ((ua, ub), (va, vb)) = surface.domain();
    [
        ((surface.is_periodic_u() || surface.is_closed_u(tol)) && ub > ua).then_some(ub - ua),
        ((surface.is_periodic_v() || surface.is_closed_v(tol)) && vb > va).then_some(vb - va),
    ]
}

/// The shifts by whole periods worth trying on a chart point after its fold,
/// nearest first: none, then one period along each periodic direction.
pub(crate) fn period_shifts(surface: &SurfaceGeometry, tol: Tolerances) -> Vec<(f64, f64)> {
    let mut out = vec![(0.0, 0.0)];
    let [span_u, span_v] = chart_spans(surface, tol);
    if let Some(span) = span_u {
        out.push((span, 0.0));
        out.push((-span, 0.0));
    }
    if let Some(span) = span_v {
        out.push((0.0, span));
        out.push((0.0, -span));
    }
    out
}

/// [`period_shifts`], then the shift by whole periods along both directions
/// at once that carries `from` nearest `centre`, the middle of the face's
/// trim, where it moves the point along both.
///
/// A torus face's trim can sit a period off the window in both directions
/// (its rims' circles run past a turn, its tube angle below zero), and a
/// point folded into the window reaches it only by that diagonal shift. It
/// is the one diagonal tried: asked of every diagonal, a trim still open
/// where a pole is yet to close it reads a point a period off both ways as
/// inside.
pub(crate) fn shifts_toward(
    surface: &SurfaceGeometry,
    from: Point2,
    centre: Point2,
    tol: Tolerances,
) -> Vec<(f64, f64)> {
    let mut out = period_shifts(surface, tol);
    if let [Some(du), Some(dv)] = chart_spans(surface, tol) {
        let toward = (
            ((centre.x - from.x) / du).round() * du,
            ((centre.y - from.y) / dv).round() * dv,
        );
        if toward.0 != 0.0 && toward.1 != 0.0 {
            out.push(toward);
        }
    }
    out
}

/// The middle of the box round chart lines.
pub(crate) fn middle_of<'a>(lines: impl IntoIterator<Item = &'a [Point2]>) -> Point2 {
    let (low, high) = lines.into_iter().flatten().fold(
        (
            Point2::new(f64::INFINITY, f64::INFINITY),
            Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        ),
        |(low, high), p| {
            (
                Point2::new(low.x.min(p.x), low.y.min(p.y)),
                Point2::new(high.x.max(p.x), high.y.max(p.y)),
            )
        },
    );
    Point2::new(f64::midpoint(low.x, high.x), f64::midpoint(low.y, high.y))
}

/// A chart point folded to the side of the face's seam its trim is on.
///
/// The surface's fold puts a point in the period starting at the surface's
/// own zero, which is where a face's trim is only when its seam runs there.
/// A patch whose seam was cut elsewhere (a band opened along the widest gap
/// its mesh left, running diagonally round) covers a period starting at the
/// seam, and a point just short of it belongs a period on. The fold is kept
/// where it lands inside the trim, or where no shift does.
///
/// The trim is the face's strands, which close only with the pieces still
/// being placed, so their even-odd count can name a copy that lies wholly
/// outside the box round them (a contact's middle level with the end of
/// the strand it continues, the ray running through that end). No copy
/// outside that box is inside. Where none inside is found and the fold
/// itself lies outside the box, the copy nearest the trim is taken: a
/// point on the face's boundary is in the column the boundary is in.
pub(crate) fn fold_inside(
    p: Point2,
    surface: &SurfaceGeometry,
    trim: &[&[Point2]],
    tol: Tolerances,
) -> Point2 {
    let folded = fold_point_into_chart(p, surface, tol);
    // A chart that does not come round has the one fold, inside or not.
    if trim.is_empty() || chart_spans(surface, tol).iter().all(Option::is_none) {
        return folded;
    }
    let (low, high) = trim.iter().flat_map(|line| line.iter()).fold(
        (
            Point2::new(f64::INFINITY, f64::INFINITY),
            Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        ),
        |(low, high), q| {
            (
                Point2::new(low.x.min(q.x), low.y.min(q.y)),
                Point2::new(high.x.max(q.x), high.y.max(q.y)),
            )
        },
    );
    let boxed = |q: &Point2| {
        q.x >= low.x - PARAM_SNAP
            && q.x <= high.x + PARAM_SNAP
            && q.y >= low.y - PARAM_SNAP
            && q.y <= high.y + PARAM_SNAP
    };
    let copies: Vec<Point2> = shifts_toward(surface, folded, middle_of(trim.iter().copied()), tol)
        .into_iter()
        .map(|(du, dv)| Point2::new(folded.x + du, folded.y + dv))
        .collect();
    if let Some(q) = copies.iter().find(|q| boxed(q) && inside_many(trim, **q)) {
        return *q;
    }
    if boxed(&folded) {
        return folded;
    }
    copies
        .into_iter()
        .map(|q| (distance_to_lines(trim, q), q))
        .min_by(|a, b| a.0.total_cmp(&b.0))
        .map_or(folded, |(_, q)| q)
}

/// The distance in the chart from `p` to the nearest of `lines`.
pub(crate) fn distance_to_lines(lines: &[&[Point2]], p: Point2) -> f64 {
    let mut best = f64::INFINITY;
    for line in lines {
        if let [only] = line {
            best = best.min(p.distance(*only));
        }
        for w in line.windows(2) {
            let (a, ab) = (w[0], w[1] - w[0]);
            let length = ab.dot(ab);
            let t = if length > 0.0 {
                ((p - a).dot(ab) / length).clamp(0.0, 1.0)
            } else {
                0.0
            };
            best = best.min(p.distance(a + ab * t));
        }
    }
    best
}

/// A strand folded into the chart as a whole, to the side of the face's
/// seam its interior lies inside the trim on, as [`fold_inside`] does for
/// a point.
pub(crate) fn fold_line_inside(
    line: &mut [Point2],
    surface: &SurfaceGeometry,
    trim: &arrange::Trim,
    tol: Tolerances,
) {
    fold_into_chart(line, surface, tol);
    if line.is_empty() || trim.is_empty() {
        return;
    }
    // The strand's middle can sit level with the boundary's own vertices
    // (a contact running across the chart ends where it paved the trim),
    // where the level ray passes through them. The leaning ray cannot.
    let mid = interior_of(line);
    if let Some((du, dv)) = shifts_toward(surface, mid, trim.middle(), tol)
        .into_iter()
        .find(|(du, dv)| trim.inside_slanted(Point2::new(mid.x + du, mid.y + dv)))
    {
        for p in line.iter_mut() {
            p.x += du;
            p.y += dv;
        }
    }
}

/// A face's boundary in its chart, each edge's image sampled, both sides of
/// a seam: enough to tell which side of the seam a point is on.
pub(crate) fn face_trim_lines(face: &GFace, tol: Tolerances) -> Vec<Vec<Point2>> {
    let sample = |pcurve: &PlanarCurve, range: (f64, f64)| -> Vec<Point2> {
        (0..=16)
            .filter_map(|k| {
                pcurve
                    .point_at(range.0 + (range.1 - range.0) * f64::from(k) / 16.0, tol)
                    .ok()
            })
            .collect()
    };
    let mut out = Vec::new();
    for e in &face.edges {
        out.push(sample(&e.pcurve, e.prange));
        if let Some((other, range)) = &e.other_side {
            out.push(sample(other, *range));
        }
    }
    out
}

/// Remove period tears from a sampled polyline, axis by axis.
pub(crate) fn unwrap_polyline(line: &mut [Point2], surface: &SurfaceGeometry, tol: Tolerances) {
    // A surface closed without being periodic (a sphere spelt as a spline)
    // tears the same way: a closed section's sample at its end folds back to
    // its start, one column across the chart from where the strand arrives.
    let [span_u, span_v] = chart_spans(surface, tol);
    let spans = (span_u.unwrap_or(0.0), span_v.unwrap_or(0.0));
    for i in 1..line.len() {
        if spans.0 > 0.0 {
            while line[i].x - line[i - 1].x > spans.0 * 0.5 {
                line[i].x -= spans.0;
            }
            while line[i].x - line[i - 1].x < -spans.0 * 0.5 {
                line[i].x += spans.0;
            }
        }
        if spans.1 > 0.0 {
            while line[i].y - line[i - 1].y > spans.1 * 0.5 {
                line[i].y -= spans.1;
            }
            while line[i].y - line[i - 1].y < -spans.1 * 0.5 {
                line[i].y += spans.1;
            }
        }
    }
}

/// Fold a point into a chart that comes round, as [`fold_into_chart`].
pub(crate) fn fold_point_into_chart(
    p: Point2,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> Point2 {
    let mut one = [p];
    fold_into_chart(&mut one, surface, tol);
    one[0]
}

/// A sub-range of a closed curve, brought into its domain.
///
/// The filler split every wrap interval at the domain end, so a piece fits
/// within one period. The fold of its start may still land the end a hair
/// past the domain, which clamps. A piece starting a hair short of the
/// domain end is not split there (the split would leave dust), so its
/// start folds to just below the end with all its length past it: it is
/// read from the domain start instead, or clamping would shrink it to the
/// hair and weld both its ends into one vertex.
pub(crate) fn folded_range(range: (f64, f64), domain: (f64, f64), closed: bool) -> (f64, f64) {
    if !closed {
        return range;
    }
    let length = range.1 - range.0;
    let mut f0 = fold(range.0, domain);
    if f0 + length - domain.1 > domain.1 - f0 {
        f0 = domain.0.max(f0 - (domain.1 - domain.0));
    }
    let f1 = (f0 + length).min(domain.1);
    (f0, f1)
}

/// Fold a parameter into a closed curve's domain.
pub(crate) fn fold(t: f64, domain: (f64, f64)) -> f64 {
    let span = domain.1 - domain.0;
    if span <= 0.0 {
        return domain.0;
    }
    domain.0 + (t - domain.0).rem_euclid(span)
}

/// The parameter to evaluate a section at: folded where the curve is closed
/// and its pieces may run past the domain end, left alone where it is not.
///
/// An open section's own end *is* the domain end, and folding it lands on the
/// domain start instead: the far end of the curve. On a plane's section
/// through a ball's poles, where the marcher hands back two open half circles
/// each running pole to pole, that is the difference between an edge bounded
/// by the pole it reaches and one bounded by the pole on the other side.
pub(crate) fn at_param(t: f64, domain: (f64, f64), closed: bool) -> f64 {
    if closed { fold(t, domain) } else { t }
}

/// How far apart a face's chart holds two of its edges' ends as one point:
/// its loosest edge, the vertices that edge ends at, and the holes set
/// aside from it.
pub(crate) fn doubt_of(face: &GFace) -> f64 {
    face.edges
        .iter()
        .fold(0.0_f64, |acc, e| {
            acc.max(e.standoff * 2.0)
                .max(e.ends_tolerance + e.tolerance)
        })
        .max(
            face.holes_aside
                .iter()
                .fold(0.0_f64, |acc, hole| acc.max(hole.doubt)),
        )
}

/// How far apart a face's chart holds two points as one at `p`: the doubt
/// of each edge `p` lies within reach of (twice what it stands off the face
/// along it, its tolerance and its vertex's at either end), the face's own
/// tolerance, and the holes set aside from it. A loose edge welds what
/// crosses or ends near it, and leaves a crossing of a tight edge elsewhere
/// on the face where it is. Never more than [`doubt_of`], the face's
/// loosest.
pub(crate) fn doubt_at(face: &GFace, p: Point, tol: Tolerances) -> OgeomResult<f64> {
    let mut doubt = face
        .holes_aside
        .iter()
        .fold(face.tolerance, |acc, hole| acc.max(hole.doubt));
    for e in &face.edges {
        for (end, radius) in e.ends {
            let reach = radius + e.tolerance;
            if p.distance(end) <= reach {
                doubt = doubt.max(reach);
            }
        }
        let along = e.standoff * 2.0;
        if along > doubt
            && e.bound.distance_to(p) <= along
            && distance_to_edge_curve(&e.curve, e.crange, p, tol)? <= along
        {
            doubt = along;
        }
    }
    Ok(doubt.min(doubt_of(face)))
}

/// One hug decision for each stretch of an edge, shared by every section
/// running along it.
///
/// Where a face of one solid passes along an edge of the other, each face
/// of the edge meets it in its own section, and the two sections run
/// between the same crossings of the edge. Each section's hug is asked of
/// its own face's edges at its own width (a fitted section's is wider than
/// an exact one's), and the face both sections cross welds its chart at the
/// loosest doubt of anything it meets. Decided apart, one face of the edge
/// keeps the edge while the other keeps its section, or both keep their
/// sections while the crossed face welds the sliver between them to
/// nothing, and the kept pieces do not close.
///
/// Two such pieces are decided together. Where both stay along the edge
/// and along each other within the widest of their hug widths and the
/// weld at each point along them, the stretch is the edge's on both faces:
/// both pieces hug it, admitted onto the crossed face as a hugging piece
/// is. The weld at a point is the doubt of what lies near it on the three
/// faces ([`doubt_at`]): a loose edge at one end of the stretch does not
/// weld the rest of it.
pub(crate) fn share_hugs(
    pieces: &mut [SectionPiece],
    sections: &[SectionRec],
    solids: [&GSolid; 2],
    tol: Tolerances,
) -> OgeomResult<()> {
    let point_of = |section: &SectionRec, t: f64| -> OgeomResult<Point> {
        let domain = section.curve.domain();
        section
            .curve
            .point_at(at_param(t, domain, section.closed), tol)
    };
    let face_of = |s: &SectionRec, side: usize| if side == 0 { s.face_a } else { s.face_b };
    // Each face's weld: its own doubt and that of every face its sections
    // end on.
    let mut weld: [Vec<f64>; 2] = [
        solids[0].faces.iter().map(doubt_of).collect(),
        solids[1].faces.iter().map(doubt_of).collect(),
    ];
    for s in sections {
        let (a, b) = (
            doubt_of(&solids[0].faces[s.face_a]),
            doubt_of(&solids[1].faces[s.face_b]),
        );
        weld[0][s.face_a] = weld[0][s.face_a].max(b);
        weld[1][s.face_b] = weld[1][s.face_b].max(a);
    }
    // The width a piece already hugs at, or none.
    let hug_width = |p: &SectionPiece, side: usize| -> f64 {
        if p.hugs[side] {
            (tol.confusion() * 1e3).max(sections[p.section].tolerance * 3.0)
        } else {
            0.0
        }
    };
    // Whether the piece stays along the edge within the weld `near` sets
    // at each point of its length.
    let along_edge = |p: &SectionPiece,
                      e: &BoundaryEdge,
                      near: &dyn Fn(Point) -> OgeomResult<f64>|
     -> OgeomResult<bool> {
        let s = &sections[p.section];
        for i in 0..=4 {
            let t = p.range.0 + (p.range.1 - p.range.0) * f64::from(i) / 4.0;
            let at = point_of(s, t)?;
            let (drawn, drange, off) = e.drawn();
            let reach = near(at)?.max(off * 2.0);
            if distance_to_edge_curve(drawn, drange, at, tol)? > reach {
                return Ok(false);
            }
        }
        Ok(true)
    };
    let stretch =
        |s: &SectionRec, p: &SectionPiece| folded_range(p.range, s.curve.domain(), s.closed);
    let middles: Vec<Point> = pieces
        .iter()
        .map(|p| point_of(&sections[p.section], f64::midpoint(p.range.0, p.range.1)))
        .collect::<OgeomResult<_>>()?;
    let mut adopted: Vec<(usize, usize, EdgeKey)> = Vec::new();
    for side in 0..2 {
        let own = solids[side];
        for (x, px) in pieces.iter().enumerate() {
            if px.hugs[1 - side] {
                continue;
            }
            let sx = &sections[px.section];
            let (fx, crossed) = (face_of(sx, side), face_of(sx, 1 - side));
            for (y, py) in pieces.iter().enumerate().skip(x + 1) {
                if py.hugs[1 - side] || py.section == px.section || (px.hugs[side] && py.hugs[side])
                {
                    continue;
                }
                let sy = &sections[py.section];
                let fy = face_of(sy, side);
                if fy == fx || face_of(sy, 1 - side) != crossed {
                    continue;
                }
                let hugged = hug_width(px, side).max(hug_width(py, side));
                let widest = weld[1 - side][crossed].max(hugged);
                let (mx, my) = (middles[x], middles[y]);
                if mx.distance(my) > 2.0 * widest
                    || distance_to_edge_curve(&sx.curve, stretch(sx, px), my, tol)? > widest
                    || distance_to_edge_curve(&sy.curve, stretch(sy, py), mx, tol)? > widest
                {
                    continue;
                }
                let faces = [
                    &own.faces[fx],
                    &own.faces[fy],
                    &solids[1 - side].faces[crossed],
                ];
                let near = |at: Point| -> OgeomResult<f64> {
                    let mut near = hugged;
                    for face in faces {
                        near = near.max(doubt_at(face, at, tol)?);
                    }
                    Ok(near)
                };
                let middle = near(mx)?.max(near(my)?);
                if mx.distance(my) > 2.0 * middle
                    || distance_to_edge_curve(&sx.curve, stretch(sx, px), my, tol)? > middle
                    || distance_to_edge_curve(&sy.curve, stretch(sy, py), mx, tol)? > middle
                {
                    continue;
                }
                for ex in &own.faces[fx].edges {
                    let Some(ey) = own.faces[fy].edges.iter().find(|e| e.node == ex.node) else {
                        continue;
                    };
                    if along_edge(px, ex, &near)? && along_edge(py, ey, &near)? {
                        if !px.hugs[side] {
                            adopted.push((x, side, ex.node));
                        }
                        if !py.hugs[side] {
                            adopted.push((y, side, ey.node));
                        }
                        break;
                    }
                }
            }
        }
    }
    for (y, side, key) in adopted {
        let piece = &mut pieces[y];
        if piece.hugs[0] || piece.hugs[1] {
            continue;
        }
        let section = &sections[piece.section];
        ogeom_core::diag!(
            Channel::Wire,
            "PAVE s{}: ({:.6}, {:.6}) hugs the edge with the section beside it",
            piece.section,
            piece.range.0,
            piece.range.1
        );
        piece.hugs[side] = true;
        piece.hug_key = Some(if side == 0 {
            (key, false, section.face_b)
        } else {
            (key, true, section.face_a)
        });
    }
    Ok(())
}

/// The move of each face of `own` standing on a drum that touches the
/// planes of `other`'s faces it meets along lines, each to within the
/// confusion distance as
/// [`drum_onto_plane`](ogeom_intersect::drum_onto_plane) reads the pair:
/// the one move across the axis that lays the drum on all of them. The
/// faces on one surface move together, the planes any of them meets asked
/// of all. `None` for a face no plane touches so, or whose planes ask for
/// moves no one move makes that is no longer than the longest of them.
///
/// The touching line is the pair's whole section, and every other section
/// of the drum (against a round the plane runs into, against the faces the
/// line ends on) must end where that line does, as it would on a drum
/// touching the plane exactly. Moved by no more than the confusion
/// distance, the drum stands within its edges' tolerance of where it was,
/// its chart unchanged, so the sections found on it are the drum's own and
/// their pcurves read on it as they stand.
pub(crate) fn drums_onto_planes(
    own: &GSolid,
    other: &GSolid,
    tol: Tolerances,
) -> Vec<Option<ogeom_math::Vector>> {
    let mut numbered: ogeom_core::FastMap<&(ogeom_topo::SurfaceId, Location), usize> =
        ogeom_core::FastMap::default();
    let mut asked: Vec<Vec<(ogeom_math::Vector, f64)>> = Vec::new();
    let mut of_face = Vec::with_capacity(own.faces.len());
    for face in &own.faces {
        let next = numbered.len();
        let at = *numbered.entry(&face.carrier).or_insert(next);
        if at == asked.len() {
            asked.push(Vec::new());
        }
        of_face.push(at);
        let SurfaceGeometry::Cylinder(drum) = &face.surface else {
            continue;
        };
        let drum = drum.cylinder();
        for theirs in &other.faces {
            let SurfaceGeometry::Plane(plane) = &theirs.surface else {
                continue;
            };
            if !face.filter.intersects(&theirs.filter) {
                continue;
            }
            let plane = plane.plane();
            if let Some(shift) = ogeom_intersect::drum_onto_plane(&plane, &drum, tol) {
                let normal = plane.normal().vector();
                asked[at].push((normal, shift.dot(normal)));
            }
        }
    }
    let moves: Vec<Option<ogeom_math::Vector>> = asked
        .iter()
        .map(|asked| {
            let &(n0, s0) = asked.first()?;
            // Every normal lies across the axis. The move is along the
            // first normal, and along the first normal not parallel to it
            // as well.
            let mut shift = n0 * s0;
            if let Some(&(n1, s1)) = asked
                .iter()
                .find(|(n, _)| n.cross(n0).magnitude() > tol.angular())
            {
                let (a, b, c) = (n0.dot(n0), n0.dot(n1), n1.dot(n1));
                let det = a.mul_add(c, -(b * b));
                shift = n0 * ((s0 * c - s1 * b) / det) + n1 * ((s1 * a - s0 * b) / det);
            }
            let agrees = asked
                .iter()
                .all(|&(n, s)| (n.dot(shift) - s).abs() <= tol.intersection());
            let longest = asked.iter().map(|(_, s)| s.abs()).fold(0.0, f64::max);
            (agrees && shift.magnitude() <= longest + tol.intersection()).then_some(shift)
        })
        .collect();
    of_face.into_iter().map(|at| moves[at]).collect()
}

/// Each face of `solid` moved as `moves` says (see [`drums_onto_planes`]),
/// with the edges that move along with it: an edge moves where every face
/// holding it moves by the same shift, or holds it on a plane the shift
/// runs along, which the edge moved by it stays on. A seam, a drum's rim
/// shared with its own flat cap, and the vertex where they meet move with
/// the drum, and the solid stays closed where it stood. An edge a face
/// set aside holds, or a face that does not move and takes the edge off
/// its surface, stays where it stands, which the moved faces then reach
/// within the move.
pub(crate) fn move_drums(
    solid: &mut GSolid,
    moves: &[Option<ogeom_math::Vector>],
    tol: Tolerances,
) -> OgeomResult<()> {
    if moves.iter().all(Option::is_none) {
        return Ok(());
    }
    let mut holders: ogeom_core::FastMap<EdgeKey, Vec<(usize, usize)>> =
        ogeom_core::FastMap::default();
    for (f, face) in solid.faces.iter().enumerate() {
        for (k, e) in face.edges.iter().enumerate() {
            holders.entry(e.node).or_default().push((f, k));
        }
    }
    let along_plane = |face: &GFace, shift: ogeom_math::Vector| match &face.surface {
        SurfaceGeometry::Plane(plane) => {
            plane.plane().normal().vector().dot(shift).abs() <= tol.angular() * shift.magnitude()
        }
        _ => false,
    };
    let mut shift_of: ogeom_core::FastMap<EdgeKey, ogeom_math::Vector> =
        ogeom_core::FastMap::default();
    for (key, held) in &holders {
        if solid.shared.contains_key(key) {
            continue;
        }
        let Some(shift) = held.iter().find_map(|&(f, _)| moves[f]) else {
            continue;
        };
        if held.iter().all(|&(f, _)| match moves[f] {
            Some(other) => other == shift,
            None => along_plane(&solid.faces[f], shift),
        }) {
            shift_of.insert(*key, shift);
        }
    }
    // A vertex moves where only moving edges end on it.
    let staying: ogeom_core::FastSet<[u64; 3]> = solid
        .faces
        .iter()
        .flat_map(|face| &face.edges)
        .filter(|e| !shift_of.contains_key(&e.node))
        .flat_map(|e| e.vertices)
        .map(|p| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()])
        .collect();
    let mut curves: ogeom_core::FastMap<EdgeKey, std::sync::Arc<Curve>> =
        ogeom_core::FastMap::default();
    for (face, own) in solid.faces.iter_mut().zip(moves) {
        let far = own
            .or_else(|| {
                face.edges
                    .iter()
                    .find_map(|e| shift_of.get(&e.node).copied())
            })
            .map_or(0.0, |shift| shift.magnitude());
        if far == 0.0 {
            continue;
        }
        if let Some(shift) = own {
            face.surface = face
                .surface
                .transformed(&ogeom_math::Transform::translation(*shift), tol)?;
        }
        face.bound = face.bound.expanded(far);
        face.filter = face.filter.expanded(far);
        let chart_shift = |shift: ogeom_math::Vector| match &face.surface {
            SurfaceGeometry::Plane(plane) => {
                let frame = plane.plane().frame();
                Some(ogeom_math::Transform2::translation(
                    ogeom_math::Vector2::new(
                        shift.dot(frame.x().vector()),
                        shift.dot(frame.y().vector()),
                    ),
                ))
            }
            _ => None,
        };
        for e in &mut face.edges {
            let moved = shift_of.get(&e.node).copied();
            // The face's own image of an edge is its pcurve on its surface,
            // which moves with the surface.
            if let Some(image) = &e.image
                && let Some(shift) = moved.or(*own)
            {
                e.image = Some(std::sync::Arc::new(
                    image.transformed(&ogeom_math::Transform::translation(shift), tol)?,
                ));
            }
            let Some(shift) = moved else {
                continue;
            };
            let moving = ogeom_math::Transform::translation(shift);
            e.curve = match curves.get(&e.node) {
                Some(curve) => curve.clone(),
                None => {
                    let curve = std::sync::Arc::new(e.curve.transformed(&moving, tol)?);
                    curves.insert(e.node, curve.clone());
                    curve
                }
            };
            // On a plane that stays, the edge moves in the plane's chart.
            if own.is_none()
                && let Some(by) = chart_shift(shift)
            {
                e.pcurve = e.pcurve.transformed(&by, tol)?;
            }
            for end in &mut e.ends {
                end.0 = moving.apply(end.0);
            }
            for at in &mut e.vertices {
                if !staying.contains(&[at.x.to_bits(), at.y.to_bits(), at.z.to_bits()]) {
                    *at = moving.apply(*at);
                }
            }
            e.bound = e.bound.expanded(far);
        }
    }
    Ok(())
}
