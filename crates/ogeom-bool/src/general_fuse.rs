//! The general fuse: every face of each argument split by the other's and
//! its pieces classified.

use super::*;

pub(crate) fn general_fuse(
    model: &Model,
    a: &Shape,
    b: &Shape,
    aside: Option<&[aside::Aside; 2]>,
    tol: Tolerances,
) -> OgeomResult<GeneralFused> {
    general_fuse_as(model, a, b, Operands::Solids, aside, tol)
}

/// [`general_fuse`] for an operation that keeps pieces by their state,
/// refusing a classification that cannot hold: the boundaries meet along
/// sections, every piece of both reads outside the other, and yet the
/// solids overlap.
///
/// Two solids sharing volume share it within a region bounded by pieces of
/// their boundaries, and each such piece lies inside the other solid or on
/// its boundary. So where no piece is in or on, the solids are apart, and a
/// section that leaves them overlapping was not arranged into the faces it
/// splits (a loop through a point where the surfaces touch, read as one
/// ring round both its lobes). Kept as read, the cut would keep the whole
/// solid and the union both solids whole. The overlap is asked of sampled
/// boundaries, only where every piece reads out and sections were found.
pub(crate) fn general_fuse_classified(
    model: &Model,
    a: &Shape,
    b: &Shape,
    aside: Option<&[aside::Aside; 2]>,
    tol: Tolerances,
) -> OgeomResult<GeneralFused> {
    let fused = general_fuse(model, a, b, aside, tol)?;
    #[cfg(feature = "test-hooks")]
    work::add(|w| w.pieces += fused.pieces.len());
    if !fused.sections.is_empty()
        && fused.pieces.iter().all(|p| p.state == PieceState::Out)
        && !apart(model, a, b, tol)?
    {
        ogeom_bail!(
            NotDone,
            "the boundaries cross yet every piece of both reads outside the other, \
             though the solids overlap; a section was not arranged into the faces it \
             splits"
        );
    }
    Ok(fused)
}

/// For each gathered face, the faces of the same solid it shares a vertex
/// with, itself among them, by their places in `solid.faces`.
pub(crate) fn neighbouring_faces(
    model: &Model,
    solid: &GSolid,
) -> OgeomResult<Vec<ogeom_core::FastSet<usize>>> {
    // A vertex is its node where it is placed, as an edge is (see
    // [`EdgeKey`]): one node placed twice is two corners in space.
    let mut by_vertex: ogeom_core::FastMap<EdgeKey, Vec<usize>> = ogeom_core::FastMap::default();
    for (index, face) in solid.faces.iter().enumerate() {
        for vertex in explore_unique(model, &face.face, ShapeType::Vertex)? {
            by_vertex
                .entry(EdgeKey::of(&vertex))
                .or_default()
                .push(index);
        }
    }
    let mut neighbours: Vec<ogeom_core::FastSet<usize>> = (0..solid.faces.len())
        .map(|index| std::iter::once(index).collect())
        .collect();
    for faces in by_vertex.values() {
        for &i in faces {
            neighbours[i].extend(faces.iter().copied());
        }
    }
    Ok(neighbours)
}

/// The general fuse of `a` and `b` taken as `operands` says. Where `a` is a
/// sheet only its pieces are made, and the pieces of `b` are not. The
/// faces `aside` names are left out.
pub(crate) fn general_fuse_as(
    model: &Model,
    a: &Shape,
    b: &Shape,
    operands: Operands,
    aside: Option<&[aside::Aside; 2]>,
    tol: Tolerances,
) -> OgeomResult<GeneralFused> {
    // Read here, on the caller's thread: the pieces are classified on
    // worker threads, which do not share it.
    let settle = SETTLE_FROM_BOTH_SIDES.get();
    let read_at_centre = PROBE_AT_CENTRE.get();
    let read_thin = AT_STATED.get();
    ogeom_core::progress::stage("boolean: gather");
    let mut ga = gather(
        model,
        a,
        matches!(operands, Operands::SheetBySolid | Operands::SheetBySheet),
        aside.map(|aside| &aside[0]),
        tol,
    )?;
    let mut gb = gather(
        model,
        b,
        operands == Operands::SheetBySheet,
        aside.map(|aside| &aside[1]),
        tol,
    )?;
    let (onto_a, onto_b) = (
        drums_onto_planes(&ga, &gb, tol),
        drums_onto_planes(&gb, &ga, tol),
    );
    move_drums(&mut ga, &onto_a, tol)?;
    move_drums(&mut gb, &onto_b, tol)?;
    #[cfg(feature = "test-hooks")]
    work::add(|w| w.faces_gathered += ga.faces.len() + gb.faces.len());
    let neighbours = if operands == Operands::Itself {
        Some(neighbouring_faces(model, &ga)?)
    } else {
        None
    };
    ogeom_core::progress::stage("boolean: intersect");
    let (
        sections,
        section_pieces,
        contacts,
        tangents,
        contact_along,
        mut paves,
        same_a,
        same_b,
        hugs,
    ) = fill(&ga, &gb, false, neighbours.as_deref(), tol)?;
    if operands == Operands::SheetBySheet && same_a.iter().any(|partners| !partners.is_empty()) {
        ogeom_bail!(
            NotDone,
            "the sheets lie on one surface where they overlap; a sheet split \
             by another lying along it is not resolved"
        );
    }
    for section in &sections {
        for (pcurve, face) in [
            (&section.pc_a, &ga.faces[section.face_a]),
            (&section.pc_b, &gb.faces[section.face_b]),
        ] {
            past_a_pole(section, pcurve, face, tol)?;
        }
    }
    let along = crossed_along(&contacts);
    let mut junctions = pave_junctions(&ga, &gb, &paves, &along, tol)?;
    if diag::enabled(Channel::Wire) {
        for j in &junctions {
            ogeom_core::diag!(
                Channel::Wire,
                "JUNCTION from paves at {:?} reach {:.3e}",
                j.at,
                j.reach
            );
        }
        for j in &hugs {
            ogeom_core::diag!(
                Channel::Wire,
                "JUNCTION from hugs at {:?} reach {:.3e}",
                j.at,
                j.reach
            );
        }
    }
    junctions.extend(hugs);
    // A vertex of an edge shared with a face set aside is that vertex in
    // the rebuild, which every strand ending within its span names.
    let mut kept_bins = bins::Bins::new(tol.confusion() * 1e2);
    let kept_points: Vec<Point> = ga
        .shared_vertices
        .iter()
        .chain(&gb.shared_vertices)
        .map(|(at, _)| *at)
        .collect();
    for (index, at) in kept_points.iter().enumerate() {
        kept_bins.insert(*at, index);
    }
    let mut near_kept = Vec::new();
    for face in ga.faces.iter().chain(gb.faces.iter()) {
        // An input vertex that owns a span (the corner an earlier boolean
        // welded three rims into, each ending a fraction of a micron from
        // the others) is a junction here as well: the rebuild's positional
        // weld reaches only a few honesties, and minted separately those
        // ends stand as three vertices with hairlines between them.
        for e in &face.edges {
            for (at, radius) in e.ends {
                let mut kept = || {
                    !kept_points.is_empty()
                        && if kept_bins.near_into(at, radius, &mut near_kept) {
                            near_kept
                                .iter()
                                .any(|&i| kept_points[i].distance(at) <= radius)
                        } else {
                            kept_points.iter().any(|p| p.distance(at) <= radius)
                        }
                };
                if radius > tol.confusion() * 10.0 && !kept() {
                    ogeom_core::diag!(
                        Channel::Wire,
                        "JUNCTION from vertex at {at:?} radius {radius:.3e}"
                    );
                    junctions.push(Junction {
                        at,
                        reach: radius,
                        onto_vertex: false,
                    });
                }
            }
        }
    }
    // Junctions whose balls overlap are one junction: two pave clusters a
    // fraction of a micron apart at a corner where rims meet, each owning
    // a span that reaches into the other's, would weld a strand end to
    // whichever covers it first, and the corner would stand as two
    // vertices with a hairline between them, kept on one face, dropped
    // on another. Merged, the corner is one vertex whose reach covers
    // every member's span, on every face alike.
    junctions = merge_junctions(junctions);
    for (from_a, gathered) in [(true, &ga), (false, &gb)] {
        seams_through_junctions(
            gathered,
            from_a,
            (&sections, &section_pieces),
            &junctions,
            &mut paves,
            tol,
        )?;
    }
    if diag::enabled(Channel::Wire) {
        for j in &junctions {
            ogeom_core::diag!(
                Channel::Wire,
                "JUNCTION at {:?} reach {:.3e}",
                j.at,
                j.reach
            );
        }
        let kind = |sf: &SurfaceGeometry| match sf {
            SurfaceGeometry::Plane(_) => "plane",
            SurfaceGeometry::Cylinder(_) => "cyl",
            SurfaceGeometry::BSpline(_) => "bspline",
            _ => "other",
        };
        for (si, sec) in sections.iter().enumerate() {
            let fa = &ga.faces[sec.face_a];
            let fb = &gb.faces[sec.face_b];
            let d = sec.curve.domain();
            let mut worst_a: f64 = 0.0;
            let mut worst_b: f64 = 0.0;
            let mut len = 0.0;
            let mut prev: Option<Point> = None;
            for i in 0..=32 {
                let t = if i == 32 {
                    d.1
                } else {
                    d.0 + (d.1 - d.0) * f64::from(i) / 32.0
                };
                let at = sec.curve.point_at(t, tol)?;
                if let Some(p) = prev {
                    len += p.distance(at);
                }
                prev = Some(at);
                let near = |f: &GFace| -> OgeomResult<f64> {
                    let mut best = f64::INFINITY;
                    for e in &f.edges {
                        best = best.min(distance_to_edge_curve(&e.curve, e.crange, at, tol)?);
                    }
                    Ok(best)
                };
                worst_a = worst_a.max(near(fa)?);
                worst_b = worst_b.max(near(fb)?);
            }
            let pieces = section_pieces.iter().filter(|p| p.section == si).count();
            ogeom_core::diag!(
                Channel::Wire,
                "SECTION s{si} a{}({}) x b{}({}) tol {:.2e} closed {} len {:.4} pieces {pieces} hug-dist a {:.3e} b {:.3e}",
                sec.face_a,
                kind(&fa.surface),
                sec.face_b,
                kind(&fb.surface),
                sec.tolerance,
                sec.closed,
                len,
                worst_a,
                worst_b
            );
        }
        ogeom_core::diag!(
            Channel::Wire,
            "TANGENTS {}  CONTACTS {}",
            tangents.len(),
            contacts.len()
        );
    }
    // The strict audit: fill again with every pair admitted and demand the
    // same kept material. The production result above is always the
    // filtered run (under audit too), so the audit compares rather than
    // substitutes, and a filtered run *is* replayed exactly. Zero cost with
    // the variable unset.
    if diag::enabled(Channel::BoolAudit) {
        // A pair only the unfiltered fill tries may refuse (a march beside a
        // pole the filter never sends it), which says nothing of the
        // filter. The audit says so and compares nothing.
        match fill(&ga, &gb, true, neighbours.as_deref(), tol) {
            Ok((audit_sections, audit_pieces, ..)) => audit_fill_equivalence(
                (&sections, &section_pieces),
                (&audit_sections, &audit_pieces),
                tol,
            ),
            Err(e @ ogeom_core::OgeomError::Cancelled) => return Err(e),
            Err(e) => ogeom_core::diag!(
                Channel::BoolAudit,
                "bound-filter audit skipped: the unfiltered fill refused: {e}"
            ),
        }
    }

    ogeom_core::progress::stage("boolean: split");
    let mut pieces: Vec<FacePiece> = Vec::new();
    // Every face's strands, before any face is arranged: the dust decision
    // (a piece shorter than its face's snap collapses to one node) must
    // be one decision per shared edge piece, and a sub-piece of one edge
    // lies in several charts, each with its own snap and its own metric.
    // A band's rail hugging a wedge's cylinder for a quarter of a micron
    // is dust in the wall's chart and a strand in the cylinder's, and the
    // sew then finds the cylinder's piece used once.
    let chart_length =
        |line: &[Point2]| -> f64 { line.windows(2).map(|w| w[0].distance(w[1])).sum::<f64>() };
    // Each edge's end vertices where they stand, each with its radius and
    // the edge's own: a contact strand on the other solid's face ends at
    // one of them.
    let mut edge_vertices: ogeom_core::FastMap<EdgeKey, [(Point, f64); 2]> =
        ogeom_core::FastMap::default();
    for e in ga
        .faces
        .iter()
        .chain(gb.faces.iter())
        .flat_map(|f| f.edges.iter())
    {
        edge_vertices.entry(e.node).or_insert([
            (e.vertices[0], e.ends[0].1 + e.tolerance),
            (e.vertices[1], e.ends[1].1 + e.tolerance),
        ]);
    }
    let strands_of = |from_a: bool,
                      fi: usize,
                      face: &GFace|
     -> OgeomResult<(Vec<Strand<Tag>>, f64, Vec<Junction>)> {
        let mut strands: Vec<Strand<Tag>> = Vec::new();
        // How far apart the paves one junction stands for lie. A section
        // meeting the edge all but tangentially (its fitted end sliding
        // along it) paves it a few microns from where an exact section
        // crosses. The first pave speaks for both, and the other's strand
        // must still reach it.
        // Each such junction where it stands in space and in this chart,
        // and how far its members reach: the strands ending at a member are
        // carried onto it below.
        let mut gathered: Vec<(Point, Point2, f64)> = Vec::new();
        for (ei, e) in face.edges.iter().enumerate() {
            let mut stops = vec![e.crange.0];
            if let Some(ts) = paves.get(&e.node) {
                // Paves the edge itself cannot tell apart are one
                // junction, and the cluster's representative speaks for it.
                let off = |face, p| off_face(&ga, &gb, face, p, tol);
                let along = along.get(&e.node).map_or(&[][..], Vec::as_slice);
                for c in cluster_paves(&e.curve, e.crange, e.tolerance, ts, Some(&off), along, tol)?
                {
                    if c.members > 1 {
                        let chart = e.pcurve.point_at(rescale(c.t, e.crange, e.prange), tol)?;
                        gathered.push((c.at, chart, c.span + c.honesty));
                    }
                    stops.push(c.t);
                }
            }
            stops.push(e.crange.1);
            // A closed boundary edge (a cap's full circle) needs two
            // distinct endpoints per strand.
            let closed_edge = e
                .curve
                .point_at(e.crange.0, tol)?
                .distance(e.curve.point_at(e.crange.1, tol)?)
                <= tol.confusion();
            if closed_edge && stops.len() == 2 {
                stops.insert(1, f64::midpoint(e.crange.0, e.crange.1));
            }
            for pair in stops.windows(2) {
                let sub = (pair[0], pair[1]);
                if sub.1 - sub.0 <= tol.parametric() {
                    continue;
                }
                strands.push(Strand {
                    polyline: pcurve_polyline(
                        &e.pcurve,
                        e.prange,
                        e.crange,
                        sub,
                        &face.surface,
                        tol,
                    )?,
                    tag: Tag::Boundary {
                        edge: ei,
                        range: sub,
                    },
                    boundary: true,
                    sag: f64::INFINITY,
                });
                if let Some((other_pc, orange)) = &e.other_side {
                    strands.push(Strand {
                        polyline: pcurve_polyline(
                            other_pc,
                            *orange,
                            e.crange,
                            sub,
                            &face.surface,
                            tol,
                        )?,
                        tag: Tag::Boundary {
                            edge: ei,
                            range: sub,
                        },
                        boundary: true,
                        sag: f64::INFINITY,
                    });
                }
            }
        }
        // The trim the strands are folded inside: the boundary strands as
        // they stand and the poles whole, with their boxes, drawn once for
        // the face. The poles become strands only after the sections that
        // cut them, and a trim without them is open where it reaches one: a
        // point inside reads outside along a ray leaving through the gap.
        let pole_lines = face
            .poles
            .iter()
            .map(|pole| {
                pcurve_polyline(
                    &pole.pcurve,
                    pole.prange,
                    pole.prange,
                    pole.prange,
                    &face.surface,
                    tol,
                )
            })
            .collect::<OgeomResult<Vec<_>>>()?;
        let mut section_trim: Option<arrange::Trim> = None;
        for sp in &section_pieces {
            let section = &sections[sp.section];
            let (belongs, pcurve) = if from_a {
                (section.face_a == fi && !sp.hugs[0], &section.pc_a)
            } else {
                (section.face_b == fi && !sp.hugs[1], &section.pc_b)
            };
            if !belongs {
                continue;
            }
            let domain = section.curve.domain();
            let sub = if section.closed {
                (
                    fold(sp.range.0, domain),
                    sp.range.1 - sp.range.0 + fold(sp.range.0, domain),
                )
            } else {
                sp.range
            };
            // The pcurve shares the curve's parameterization. Sampling
            // uses folded parameters for periodic curves.
            let fold_at = |t: f64| if section.closed { fold(t, domain) } else { t };
            // A section's sampling on a plane holds its chords to a sag,
            // which the arrangement can read it to.
            let (count, sag) = if !matches!(face.surface, SurfaceGeometry::Plane(_)) {
                (32, f64::INFINITY)
            } else if matches!(pcurve, PlanarCurve::Line(_)) {
                (32, 0.0)
            } else {
                match circle_sag_count(pcurve, sub.0, sub.1, 32) {
                    Some(count) => (
                        count,
                        circle_sag(pcurve, sub.0, sub.1, count).unwrap_or(f64::INFINITY),
                    ),
                    None => sag_count(&|t| pcurve.point_at(fold_at(t), tol), sub.0, sub.1, 32)?,
                }
            };
            let mut line = Vec::with_capacity(count + 1);
            for i in 0..=count {
                #[allow(clippy::cast_precision_loss)]
                let t = sub.0 + (sub.1 - sub.0) * i as f64 / count as f64;
                let tf = if section.closed { fold(t, domain) } else { t };
                line.push(pcurve.point_at(tf, tol)?);
            }
            // Folding the parameter can tear the sampled polyline at the
            // period. Unwrap it pointwise, then bring the whole strand
            // into the chart with one shift.
            unwrap_polyline(&mut line, &face.surface, tol);
            let trim = section_trim.get_or_insert_with(|| {
                arrange::Trim::new(
                    strands
                        .iter()
                        .filter(|st| st.boundary)
                        .map(|st| st.polyline.clone())
                        .chain(pole_lines.iter().cloned())
                        .collect(),
                )
            });
            fold_line_inside(&mut line, &face.surface, trim, tol);
            strands.push(Strand {
                polyline: line,
                tag: Tag::Section {
                    section: sp.section,
                    range: sp.range,
                },
                boundary: false,
                sag,
            });
        }
        // A section's end the paving read at a pole (its crossing within
        // the window that reads a crossing as an edge's end vertex, at the
        // pole vertex of a seam) stands beside the pole in space by at most
        // that window, but in the chart that can be anywhere along the
        // pole's line, the chart's free parameter shrinking to nothing
        // there: off the pole row by the window over the radius, and along
        // it by any amount, past the line's end where it crosses the seam
        // beside the pole. It is the pole, and is carried onto the pole's
        // line at its foot, held to the line (every point of which is the
        // pole), where the pole is cut for it below.
        for strand in &mut strands {
            let Tag::Section { section, .. } = strand.tag else {
                continue;
            };
            let s = &sections[section];
            let reach = tol.confusion().max(s.tolerance * 2.0);
            let window = (reach + honest(s.tolerance, tol)).max(tol.confusion() * 1e2);
            let last = strand.polyline.len() - 1;
            for end in [0, last] {
                let at = strand.polyline[end];
                let Ok(q) = face.surface.point_at(at.x, at.y, tol) else {
                    continue;
                };
                for pole in &face.poles {
                    let PlanarCurve::Line(line) = &pole.pcurve else {
                        continue;
                    };
                    if q.distance(pole.point) > window {
                        continue;
                    }
                    let axis = line.axis();
                    let along = (at - axis.location)
                        .dot(axis.direction.vector())
                        .clamp(pole.prange.0, pole.prange.1);
                    strand.polyline[end] = axis.point_at(along);
                }
            }
        }
        // A boundary edge ending at a pole (a patch's side running into its
        // collapsed corner) ends there in space, but its pcurve may end off
        // the pole's line by the edge's tolerance over the chart's shrinking
        // scale, and anywhere along it: it is carried onto the line at its
        // foot, as a section's end is above.
        for strand in &mut strands {
            let Tag::Boundary { edge, .. } = strand.tag else {
                continue;
            };
            let window = (face.edges[edge].tolerance * 2.0).max(tol.confusion() * 1e2);
            let last = strand.polyline.len() - 1;
            for end in [0, last] {
                let at = strand.polyline[end];
                let Ok(q) = face.surface.point_at(at.x, at.y, tol) else {
                    continue;
                };
                for pole in &face.poles {
                    let PlanarCurve::Line(line) = &pole.pcurve else {
                        continue;
                    };
                    if q.distance(pole.point) > window {
                        continue;
                    }
                    let axis = line.axis();
                    let along = (at - axis.location)
                        .dot(axis.direction.vector())
                        .clamp(pole.prange.0, pole.prange.1);
                    strand.polyline[end] = axis.point_at(along);
                }
            }
        }
        // The poles, after the sections, because a section can end *on*
        // a pole (a plane through a ball's axis cuts it exactly there),
        // and the pole has to be cut where that happens or the two meet
        // at no shared node and the arrangement sees a dangling section.
        // A boundary edge ending on it cuts it too. Where any strand ends
        // on the pole, a piece of it with an end no strand meets lies past
        // the face's own boundary (the face spans only part of the pole's
        // line) and would only dangle: it is left out.
        for (pi, pole) in face.poles.iter().enumerate() {
            let mut stops = vec![pole.prange.0, pole.prange.1];
            let mut met: Vec<f64> = Vec::new();
            if let PlanarCurve::Line(line) = &pole.pcurve {
                let axis = line.axis();
                for strand in &strands {
                    if matches!(strand.tag, Tag::Pole { .. }) {
                        continue;
                    }
                    for end in [strand.polyline.first(), strand.polyline.last()]
                        .into_iter()
                        .flatten()
                    {
                        let along = (*end - axis.location).dot(axis.direction.vector());
                        let foot = axis.point_at(along);
                        if foot.distance(*end) > PARAM_SNAP {
                            continue;
                        }
                        met.push(along);
                        if along > pole.prange.0 + PARAM_SNAP && along < pole.prange.1 - PARAM_SNAP
                        {
                            stops.push(along);
                        }
                    }
                }
            }
            stops.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
            stops.dedup_by(|a, b| (*a - *b).abs() <= PARAM_SNAP);
            let meets = |t: f64| met.iter().any(|m| (m - t).abs() <= PARAM_SNAP);
            for pair in stops.windows(2) {
                if !met.is_empty() && !(meets(pair[0]) && meets(pair[1])) {
                    continue;
                }
                let sub = (pair[0], pair[1]);
                strands.push(Strand {
                    polyline: pcurve_polyline(
                        &pole.pcurve,
                        pole.prange,
                        pole.prange,
                        sub,
                        &face.surface,
                        tol,
                    )?,
                    tag: Tag::Pole {
                        pole: pi,
                        range: sub,
                    },
                    boundary: true,
                    sag: f64::INFINITY,
                });
            }
        }
        let mut contact_trim: Option<arrange::Trim> = None;
        for (ci, contact) in contacts.iter().enumerate() {
            if contact.target_from_a != from_a || contact.target_face != fi {
                continue;
            }
            // The owner's paves split its edge. The contact strands split
            // at the same parameters, so the sub-edges rebuilt from both
            // sides are the same edges and sew shared.
            let mut stops = vec![contact.crange.0];
            if let Some(ts) = paves.get(&contact.node) {
                stops.extend(
                    cluster_paves(
                        &contact.curve,
                        contact.crange,
                        contact.tolerance,
                        ts,
                        None,
                        &[],
                        tol,
                    )?
                    .iter()
                    .map(|c| c.t),
                );
            }
            stops.push(contact.crange.1);
            let closed_contact = contact
                .curve
                .point_at(contact.crange.0, tol)?
                .distance(contact.curve.point_at(contact.crange.1, tol)?)
                <= tol.confusion();
            if closed_contact && stops.len() == 2 {
                stops.insert(1, f64::midpoint(contact.crange.0, contact.crange.1));
            }
            for pair in stops.windows(2) {
                let sub = (pair[0], pair[1]);
                if sub.1 - sub.0 <= tol.parametric() {
                    continue;
                }
                let mid_t = f64::midpoint(sub.0, sub.1);
                if contact_along[ci]
                    .iter()
                    .any(|(lo, hi)| mid_t >= *lo && mid_t <= *hi)
                {
                    // Already boundary on both sides.
                    ogeom_core::diag!(
                        Channel::Strands,
                        "CONTACT c{ci} {sub:?} on fi={fi}: along boundary, skipped"
                    );
                    continue;
                }
                // Keep only what lies inside this face's trim. The rest
                // of the owner's boundary splits nothing here.
                let mut line = pcurve_polyline(
                    &contact.pcurve,
                    contact.prange,
                    contact.crange,
                    sub,
                    &face.surface,
                    tol,
                )?;
                unwrap_polyline(&mut line, &face.surface, tol);
                let trim = contact_trim.get_or_insert_with(|| {
                    arrange::Trim::new(
                        strands
                            .iter()
                            .filter(|st| st.boundary)
                            .map(|st| st.polyline.clone())
                            .collect(),
                    )
                });
                fold_line_inside(&mut line, &face.surface, trim, tol);
                let mid = interior_of(&line);
                // On the last attempt, where the stretch stands nearest an
                // edge's inside, it is inside the face just where it stands
                // on the side the face keeps its material on: the edge's
                // own curve says so where the outline's chords sag past it.
                let beside = if read_thin {
                    inside_beside_edge(face, contact.curve.point_at(mid_t, tol)?, tol)?
                } else {
                    None
                };
                if !beside.unwrap_or_else(|| trim.inside_slanted(mid)) {
                    ogeom_core::diag!(
                        Channel::Strands,
                        "CONTACT c{ci} {sub:?} on fi={fi}: outside the trim at {mid:?}, skipped"
                    );
                    continue;
                }
                strands.push(Strand {
                    polyline: line,
                    tag: Tag::Contact {
                        contact: ci,
                        range: sub,
                    },
                    boundary: false,
                    sag: f64::INFINITY,
                });
            }
        }

        // A tolerant contact's chart image meets the boundary it paved
        // only as closely as its own slop allows. The arrangement's node
        // weld reaches that far on this face, or the strand dangles a
        // few microns from the junction it belongs to.
        // A fitted section that hugs one of this face's edges leaves it
        // at the hug's far end by up to the hug's own width. Its strand
        // must still find the edge's node there.
        // A strand ending at another member of a junction than the one that
        // speaks for it ends where its curve all but runs along the edge (a
        // section crossing it at a shallow angle, whose crossing is only
        // known along the edge to within the run). An exact section moved
        // along to the junction stays within tolerance of both curves, and
        // carried there it meets the edge's split exactly, so the face's
        // weld need not widen to the cluster's span, which would swallow a
        // sliver beside the edge whole. A fitted one wanders by more than
        // that, and the weld widens for it.
        // On the last attempt that weld is each such end's alone: the face's
        // node weld also decides which regions are too thin to hold a point,
        // and a junction's run along one edge says nothing of a sliver
        // elsewhere on the face.
        let mut spread = 0.0_f64;
        let mut spread_at: Vec<[f64; 2]> = vec![[0.0; 2]; strands.len()];
        for (si, strand) in strands
            .iter_mut()
            .enumerate()
            .filter(|(_, st)| !st.boundary)
        {
            let exact = matches!(strand.tag, Tag::Section { section, .. }
                if sections.get(section).is_some_and(|s| s.tolerance <= 0.0));
            let last = strand.polyline.len().saturating_sub(1);
            for (side, end) in [0, last].into_iter().enumerate() {
                let Some(&at) = strand.polyline.get(end) else {
                    continue;
                };
                let Ok(q) = face.surface.point_at(at.x, at.y, tol) else {
                    continue;
                };
                let landed = gathered
                    .iter()
                    .filter(|(junction, _, reach)| junction.distance(q) <= *reach)
                    .min_by(|x, y| x.0.distance(q).total_cmp(&y.0.distance(q)));
                match landed {
                    // A junction on a seam stands on both sides of it, and
                    // the strand is carried to the side it arrives from.
                    Some(&(_, chart, _)) if exact => {
                        strand.polyline[end] = period_shifts(&face.surface, tol)
                            .into_iter()
                            .map(|(du, dv)| Point2::new(chart.x + du, chart.y + dv))
                            .min_by(|p, r| p.distance(at).total_cmp(&r.distance(at)))
                            .unwrap_or(chart);
                    }
                    Some(&(_, _, reach)) => {
                        spread = spread.max(reach);
                        spread_at[si][side] = reach;
                    }
                    None => {}
                }
            }
        }
        // What the face holds everywhere in its chart: its own tolerance
        // and its holes set aside, and twice the run's confusion, which two
        // ends each within it of one point keep between them (a fuzzy
        // run's confusion is its fuzz). What its edges and vertices add
        // holds only near them, and each strand end is welded at the doubt
        // where it stands, below.
        let everywhere = |face: &GFace| -> f64 {
            face.holes_aside
                .iter()
                .fold(face.tolerance, |acc, hole| acc.max(hole.doubt))
                .min(doubt_of(face))
        };
        let near = contacts
            .iter()
            .filter(|c| c.target_from_a == from_a && c.target_face == fi)
            .fold(PARAM_SNAP, |acc, c| acc.max(c.tolerance * 2.0))
            .max(tol.confusion() * 2.0)
            .max(everywhere(face))
            .max(if read_thin { 0.0 } else { spread })
            .max(
                if !read_thin
                    && sections.iter().any(|s| {
                        s.tolerance > 0.0
                            && if from_a {
                                s.face_a == fi
                            } else {
                                s.face_b == fi
                            }
                    })
                {
                    tol.confusion() * 1e3
                } else {
                    0.0
                },
            );
        // The face a section shares with this one.
        let other = |section: usize| -> &GFace {
            let s = &sections[section];
            if from_a {
                &gb.faces[s.face_b]
            } else {
                &ga.faces[s.face_a]
            }
        };
        let face_snap = near.max(
            (0..sections.len())
                .filter(|&k| {
                    if from_a {
                        sections[k].face_a == fi
                    } else {
                        sections[k].face_b == fi
                    }
                })
                .map(|k| everywhere(other(k)))
                .fold(0.0_f64, f64::max),
        );
        // The weld at a strand's end, in space. An edge's end lands in the
        // chart within its vertex's own recorded doubt plus its own image's:
        // a projected pcurve is honest to the edge's tolerance, and the
        // vertex it ends at was welded to some earlier gap. A section ends
        // where it crosses the other face's boundary, and lands there within
        // that boundary's own doubt: two sections through neighbouring
        // facets of a converted mesh stop on the edge they share a few
        // microns apart, each on its own facet's plane, and meet on this
        // face only if the weld reaches that far. Both doubts are those of
        // the edges and vertices near the end ([`doubt_at`]): a loose edge
        // welds what crosses or ends near it, and a section crossing a
        // tight edge elsewhere on the face keeps its crossing, as it does on
        // the face across that edge.
        let mut ends: Vec<(usize, usize, Point2, Point, f64)> = Vec::new();
        for (si, st) in strands.iter().enumerate() {
            if st.polyline.len() < 2 || matches!(st.tag, Tag::Pole { .. }) {
                continue;
            }
            let last = st.polyline.len() - 1;
            for end in [0, last] {
                let at = st.polyline[end];
                let Ok(q) = face.surface.point_at(at.x, at.y, tol) else {
                    continue;
                };
                let mut weld = doubt_at(face, q, tol)?.max(near);
                if let Tag::Section { section, .. } = st.tag {
                    weld = weld.max(doubt_at(other(section), q, tol)?);
                    // An exact section's stop read at an exact edge's end
                    // vertex (see the crossings in [`fill`]) stands beside
                    // the vertex by as far as the section passes it, up to
                    // the weld that read it there, and is that vertex.
                    let s = &sections[section];
                    if s.tolerance <= 0.0 && exact(&s.curve) {
                        let reach = tol.confusion();
                        for e in face
                            .edges
                            .iter()
                            .filter(|e| e.image.is_none() && exact(&e.curve))
                        {
                            let read_at = (reach + e.tolerance.max(e.ends_tolerance))
                                .max(tol.confusion() * 10.0);
                            for end_t in [e.crange.0, e.crange.1] {
                                let gap = e.curve.point_at(end_t, tol)?.distance(q);
                                if gap <= read_at {
                                    weld = weld.max(gap + tol.confusion());
                                }
                            }
                        }
                    }
                }
                // A contact's end is its own edge's end vertex, within that
                // vertex's radius, as an edge of this face ends at its own.
                if let Tag::Contact { contact, .. } = st.tag
                    && let Some(vertices) = edge_vertices.get(&contacts[contact].node)
                {
                    for &(v, radius) in vertices {
                        if q.distance(v) <= radius {
                            weld = weld.max(radius);
                        }
                    }
                }
                // On the last attempt the face weld leaves out a junction's
                // run and a fitted section's hug width, and they hold at the
                // ends they concern: a fitted section's end, hugging an edge
                // within the two curves' radii or crossing it, meets the
                // edge's node within its own reach, twice its tolerance.
                if read_thin {
                    weld = weld.max(spread_at[si][usize::from(end != 0)]);
                    if let Tag::Section { section, .. } = st.tag {
                        weld = weld.max(honest(sections[section].tolerance, tol) * 2.0);
                    }
                }
                ends.push((si, end, at, q, weld));
            }
        }
        // The doubts above are lengths in space, and the weld runs in the
        // chart. Where the chart stretches every direction (a sphere of a
        // few millimetres runs a radian over its whole radius), a gap in
        // space is that many times narrower in the chart, and a weld taken
        // at the space length would swallow whole edges of a small patch.
        // Where some direction shrinks instead (a pole, a thin cylinder's
        // turn), the length stands.
        let stretch = strands
            .iter()
            .flat_map(|st| st.polyline.iter())
            .step_by(8)
            .filter_map(|p| face.surface.d1_at(p.x, p.y, tol).ok())
            .map(|(du, dv)| du.magnitude().min(dv.magnitude()))
            .fold(f64::INFINITY, f64::min);
        let in_chart = |length: f64| -> f64 {
            if stretch.is_finite() && stretch > 1.0 {
                (length / stretch).max(PARAM_SNAP)
            } else {
                length
            }
        };
        let face_snap = in_chart(face_snap);
        let mut reach: Vec<f64> = ends.iter().map(|e| in_chart(e.4)).collect();
        // Where only one direction stretches (a cylinder's turn against its
        // straight length) no single scale holds, and a strand can be short
        // in the chart yet long in space. A strand longer in space than its
        // own doubt is no dust, so no weld near it reaches as far as its
        // chart length, or it collapses and takes its neighbours' ends with
        // it.
        let own_doubt = |si: usize| -> f64 {
            let at_ends = ends
                .iter()
                .filter(|e| e.0 == si)
                .fold(near, |acc, e| acc.max(e.4));
            match strands[si].tag {
                Tag::Section { section, .. } => at_ends.max(everywhere(other(section))),
                Tag::Pole { .. } => f64::INFINITY,
                Tag::Boundary { .. } | Tag::Contact { .. } => at_ends,
            }
        };
        let loosest = reach.iter().copied().fold(face_snap, f64::max);
        let kept: Vec<(Point2, Point2, f64)> = strands
            .iter()
            .enumerate()
            .filter(|(_, st)| st.polyline.len() >= 2 && chart_length(&st.polyline) <= loosest)
            .filter_map(|(si, st)| {
                let points: Vec<Point> = st
                    .polyline
                    .iter()
                    .filter_map(|p| face.surface.point_at(p.x, p.y, tol).ok())
                    .collect();
                let space: f64 = points.windows(2).map(|w| w[0].distance(w[1])).sum();
                (space > own_doubt(si)).then(|| {
                    (
                        st.polyline[0],
                        st.polyline[st.polyline.len() - 1],
                        chart_length(&st.polyline) / 2.0,
                    )
                })
            })
            .collect();
        let face_snap = kept
            .iter()
            .filter(|k| k.2 * 2.0 <= face_snap)
            .fold(face_snap, |acc, k| acc.min(k.2))
            .max(PARAM_SNAP);
        for (e, r) in ends.iter().zip(reach.iter_mut()) {
            for &(a, b, half) in &kept {
                if half * 2.0 <= *r && (e.2.distance(a) <= *r || e.2.distance(b) <= *r) {
                    *r = r.min(half);
                }
            }
            *r = r.max(face_snap);
        }
        // Each end whose own weld reaches past the face's is welded here,
        // onto one node with every end within its reach: a node of the
        // face's boundary where one stands among them, the nearest to a
        // vertex of the face first. A strand whose ends both land on the
        // node and which is no longer than the weld there is dust.
        let mut class: Vec<usize> = (0..ends.len()).collect();
        fn root(class: &mut [usize], mut i: usize) -> usize {
            while class[i] != i {
                class[i] = class[class[i]];
                i = class[i];
            }
            i
        }
        let mut welds: Vec<Junction> = Vec::new();
        let mut welded = false;
        // Two ends that each stand within a vertex's radius (its own and
        // its edge's) of the same vertex are that vertex: the curves a file
        // closes on one vertex can stop either side of it, each within the
        // radius. The vertex is one of this face's edges', or one ending
        // the other solid's edge a contact strand runs along. In the chart they stand as far apart as both their
        // reaches, or further where the chart stretches more than the
        // reaches allow for (a cylinder narrower than the radius turns a
        // quarter within it); there the stretch between them, read as an
        // edge is at seventeen points, lies within the vertex too, on a
        // chart with no pole's line: about a pole the whole chart may lie
        // within one vertex, and the pole's line is what joins the ends
        // there. Two ends at one point in space are no such gap: apart in
        // the chart, they stand on a pole's line, whose ends the pole keeps
        // apart.
        let poled = strands.iter().any(|st| matches!(st.tag, Tag::Pole { .. }));
        let one_vertex = |i: usize, j: usize| -> bool {
            let (a, b) = (ends[i].3, ends[j].3);
            let (from, to) = (ends[i].2, ends[j].2);
            if a.distance(b) <= tol.confusion() {
                return false;
            }
            let close = from.distance(to) <= reach[i] + reach[j];
            let contact_vertices = [i, j]
                .into_iter()
                .filter_map(|k| match strands[ends[k].0].tag {
                    Tag::Contact { contact, .. } => edge_vertices.get(&contacts[contact].node),
                    _ => None,
                })
                .flatten()
                .copied();
            face.edges
                .iter()
                .flat_map(|e| {
                    e.vertices
                        .iter()
                        .zip(e.ends)
                        .map(|(&v, (_, radius))| (v, radius + e.tolerance))
                })
                .chain(contact_vertices)
                .any(|(v, within)| {
                    a.distance(v) <= within
                        && b.distance(v) <= within
                        && (close
                            || !poled
                                && (1..16).all(|k| {
                                    let t = f64::from(k) / 16.0;
                                    let at = Point2::new(
                                        from.x + (to.x - from.x) * t,
                                        from.y + (to.y - from.y) * t,
                                    );
                                    face.surface
                                        .point_at(at.x, at.y, tol)
                                        .is_ok_and(|q| q.distance(v) <= within)
                                }))
                })
        };
        for i in 0..ends.len() {
            if reach[i] <= face_snap {
                continue;
            }
            for j in 0..ends.len() {
                if j != i && (ends[i].2.distance(ends[j].2) <= reach[i] || one_vertex(i, j)) {
                    let (ri, rj) = (root(&mut class, i), root(&mut class, j));
                    if ri != rj {
                        class[ri.max(rj)] = ri.min(rj);
                        welded = true;
                    }
                }
            }
        }
        if welded {
            let off_vertex = |q: Point| -> f64 {
                face.edges
                    .iter()
                    .flat_map(|e| e.ends.iter())
                    .map(|(v, _)| v.distance(q))
                    .fold(f64::INFINITY, f64::min)
            };
            let rank = |m: usize| {
                (
                    u8::from(!strands[ends[m].0].boundary),
                    off_vertex(ends[m].3),
                )
            };
            let mut node: Vec<Option<usize>> = vec![None; ends.len()];
            let mut weld_of = vec![0.0_f64; ends.len()];
            for (i, &weld) in reach.iter().enumerate() {
                let r = root(&mut class, i);
                weld_of[r] = weld_of[r].max(weld);
                let better = node[r].is_none_or(|k| {
                    let (x, y) = (rank(i), rank(k));
                    x.0.cmp(&y.0).then(x.1.total_cmp(&y.1)).is_lt()
                });
                if better {
                    node[r] = Some(i);
                }
            }
            // A node holding ends that stand apart in space is one vertex
            // there, on every face: a junction reaching each of them, unless
            // one the paving found holds them all already.
            let mut spread_of = vec![0.0_f64; ends.len()];
            let mut held = vec![true; ends.len()];
            for i in 0..ends.len() {
                let r = root(&mut class, i);
                if let Some(k) = node[r] {
                    spread_of[r] = spread_of[r].max(ends[k].3.distance(ends[i].3));
                    held[r] &= junctions.iter().any(|j| {
                        j.at.distance(ends[k].3) <= j.reach && j.at.distance(ends[i].3) <= j.reach
                    });
                }
            }
            for (r, k) in node.iter().enumerate() {
                if let Some(k) = *k
                    && spread_of[r] > tol.confusion() * 1e2
                    && !held[r]
                {
                    welds.push(Junction {
                        at: ends[k].3,
                        reach: spread_of[r] + tol.confusion() * 1e2,
                        onto_vertex: false,
                    });
                }
            }
            let lengths: Vec<f64> = strands
                .iter()
                .map(|st| chart_length(&st.polyline))
                .collect();
            let mut landed: Vec<[Option<usize>; 2]> = vec![[None, None]; strands.len()];
            for i in 0..ends.len() {
                let r = root(&mut class, i);
                let Some(k) = node[r] else {
                    continue;
                };
                let (si, end, at, _, _) = ends[i];
                let target = ends[k].2;
                // A node on a seam stands on both sides of it, and the end
                // is carried to the side it arrives from.
                strands[si].polyline[end] = period_shifts(&face.surface, tol)
                    .into_iter()
                    .map(|(du, dv)| Point2::new(target.x + du, target.y + dv))
                    .min_by(|p, q| p.distance(at).total_cmp(&q.distance(at)))
                    .unwrap_or(target);
                landed[si][usize::from(end != 0)] = Some(r);
            }
            for (si, st) in strands.iter_mut().enumerate() {
                if let [Some(a), Some(b)] = landed[si]
                    && a == b
                    && lengths[si] <= weld_of[a]
                {
                    ogeom_core::diag!(
                        Channel::Strands,
                        "WELD from_a={from_a} fi={fi}: strand {si} is dust at its ends' weld {:.3e}",
                        weld_of[a]
                    );
                    st.polyline = vec![st.polyline[0], st.polyline[0]];
                }
            }
        }
        // A section that runs along one of this face's edges into its end
        // stops within the hug's width of the node the edge was split at,
        // and two sections meeting a rail that grazes this face stop on it
        // as far apart as the rail stays within tolerance of the face. The
        // junction there says each pair is one point. Where a chart unit is
        // shorter than a millimetre, or the graze long, that gap is wider in
        // the chart than the weld, and the section dangles beside the node it
        // ends at. Its end moves onto a boundary node the junction holds, or
        // onto the first section end already there.
        let mut anchors: Vec<(Option<usize>, Point2, Point)> = strands
            .iter()
            .filter(|st| st.boundary && st.polyline.len() >= 2)
            .flat_map(|st| [st.polyline[0], st.polyline[st.polyline.len() - 1]])
            .filter_map(|p| {
                face.surface
                    .point_at(p.x, p.y, tol)
                    .ok()
                    .map(|q| (None, p, q))
            })
            .collect();
        for st in strands
            .iter_mut()
            .filter(|st| matches!(st.tag, Tag::Section { .. }) && st.polyline.len() >= 2)
        {
            let last = st.polyline.len() - 1;
            let mut moved_into: Option<usize> = None;
            for at in [0, last] {
                let p = st.polyline[at];
                if anchors
                    .iter()
                    .any(|(held, n, _)| held.is_none() && n.distance(p) <= face_snap)
                {
                    continue;
                }
                let Ok(q) = face.surface.point_at(p.x, p.y, tol) else {
                    continue;
                };
                let Some((ji, junction)) = junctions
                    .iter()
                    .enumerate()
                    .filter(|(_, j)| j.at.distance(q) <= j.reach)
                    .min_by(|(a, x), (b, y)| {
                        x.at.distance(q).total_cmp(&y.at.distance(q)).then(a.cmp(b))
                    })
                else {
                    continue;
                };
                // The nearest in space the junction also holds, and reached
                // across the chart without leaving it: the chart's midpoint
                // of the two lies midway in space too, where a node in
                // another period's copy would not.
                let anchor = anchors
                    .iter()
                    .filter(|(held, _, s)| {
                        held.is_none_or(|h| h == ji) && junction.at.distance(*s) <= junction.reach
                    })
                    .filter(|(_, n, s)| {
                        face.surface
                            .point_at((n.x + p.x) / 2.0, (n.y + p.y) / 2.0, tol)
                            .is_ok_and(|m| m.distance(q.midpoint(*s)) <= junction.reach)
                    })
                    .min_by(|a, b| {
                        a.2.distance(q)
                            .partial_cmp(&b.2.distance(q))
                            .unwrap_or(core::cmp::Ordering::Equal)
                    })
                    .map(|(_, n, _)| *n);
                match anchor {
                    Some(n) => {
                        st.polyline[at] = n;
                        moved_into = Some(ji);
                    }
                    None => anchors.push((Some(ji), p, q)),
                }
            }
            // A piece of section lying wholly inside one junction, its ends
            // carried onto one node, is that junction: two crossings a
            // shallow angle leaves a hair apart along an edge bound a stub
            // of the section between them, and kept as a strand from the
            // node back to itself it encloses a sliver with no inside.
            // Carried to a point, it is dust, and the faces either side of
            // it drop it alike.
            if let Some(ji) = moved_into
                && st.polyline[0].distance(st.polyline[last]) <= face_snap
            {
                let junction = &junctions[ji];
                let mut inside = true;
                for p in &st.polyline {
                    inside &= face
                        .surface
                        .point_at(p.x, p.y, tol)
                        .is_ok_and(|q| junction.at.distance(q) <= junction.reach);
                }
                if inside {
                    st.polyline = vec![st.polyline[0], st.polyline[0]];
                }
            }
        }
        if diag::enabled(Channel::Strands) {
            ogeom_core::diag!(
                Channel::Strands,
                "FACE-SNAP from_a={from_a} fi={fi}: snap {face_snap:.3e} edges {:?}",
                face.edges
                    .iter()
                    .map(|e| (
                        format!("{:.2e}", e.tolerance),
                        format!("{:.2e}", e.ends_tolerance)
                    ))
                    .collect::<Vec<_>>()
            );
            for (si, st) in strands.iter().enumerate() {
                let tag = match st.tag {
                    Tag::Boundary { edge, range } => format!("Boundary e{edge} {range:?}"),
                    Tag::Contact { contact, range } => format!("Contact c{contact} {range:?}"),
                    Tag::Section { section, range } => format!("Section s{section} {range:?}"),
                    Tag::Pole { pole, range } => format!("Pole p{pole} {range:?}"),
                };
                let (a, b) = (st.polyline[0], st.polyline[st.polyline.len() - 1]);
                ogeom_core::diag!(
                    Channel::Strands,
                    "STRAND from_a={from_a} fi={fi} {si}: boundary={} {tag} pts={} {a:?} .. {b:?}",
                    st.boundary,
                    st.polyline.len()
                );
            }
        }
        Ok((strands, face_snap, welds))
    };
    // The identity of a piece across faces: a boundary or contact piece by
    // its edge node, a section piece by its section: the same section is
    // split the same way on both faces it cuts, and a piece that is dust
    // in a sphere's chart at the pole must be dust on the plane it also
    // lies in.
    let node_of = |face: &GFace, tag: &Tag| -> Option<(usize, usize, (f64, f64))> {
        match tag {
            Tag::Boundary { edge, range } => {
                Some((0, face.edges[*edge].node.index() as usize, *range))
            }
            Tag::Contact { contact, range } => {
                Some((1, contacts[*contact].node.index() as usize, *range))
            }
            Tag::Section { section, range } => Some((2, *section, *range)),
            Tag::Pole { .. } => None,
        }
    };
    /// A face's strands, its snap and the junctions its welds make,
    /// waiting to be arranged.
    type Prepared = Option<(Vec<Strand<Tag>>, f64, Vec<Junction>)>;
    // A face the other solid leaves alone (no section on it, no contact
    // into it, no edge of it split, no coincident partner, no pole, no
    // edge short enough to be anyone's dust) is one piece, itself, and is
    // neither split nor arranged: its rings are its own wires and its state
    // is read once. On a large solid and a small tool that is nearly every
    // face.
    let mut quiet: [Vec<bool>; 2] = [vec![true; ga.faces.len()], vec![true; gb.faces.len()]];
    for s in &sections {
        quiet[0][s.face_a] = false;
        quiet[1][s.face_b] = false;
    }
    for c in &contacts {
        quiet[usize::from(!c.target_from_a)][c.target_face] = false;
    }
    for t in &tangents {
        quiet[0][t.face_a] = false;
        quiet[1][t.face_b] = false;
    }
    for (side, solid, same) in [(0_usize, &ga, &same_a), (1, &gb, &same_b)] {
        for (fi, face) in solid.faces.iter().enumerate() {
            if !same[fi].is_empty()
                || !face.poles.is_empty()
                || face.edges.iter().any(|e| {
                    paves.contains_key(&e.node)
                        || e.ends[0].0.distance(e.ends[1].0) <= QUIET_EDGE
                        || e.bound.diagonal() <= QUIET_EDGE
                })
            {
                quiet[side][fi] = false;
            }
        }
    }
    let mut prepared: [Vec<Prepared>; 2] = [Vec::new(), Vec::new()];
    let mut dust: Vec<(usize, usize, (f64, f64))> = Vec::new();
    for (side, solid, from_a) in [(0_usize, &ga, true), (1, &gb, false)] {
        // A sheet's tool only cuts: its own faces are not split.
        if side == 1 && operands != Operands::Solids {
            break;
        }
        for (fi, face) in solid.faces.iter().enumerate() {
            ogeom_core::progress::checkpoint()?;
            if quiet[side][fi] {
                prepared[side].push(Some((Vec::new(), 0.0, Vec::new())));
                continue;
            }
            let (strands, snap, welds) = strands_of(from_a, fi, face)?;
            for st in &strands {
                if st.polyline.len() >= 2
                    && chart_length(&st.polyline) <= snap
                    && let Some(key) = node_of(face, &st.tag)
                {
                    dust.push(key);
                }
            }
            prepared[side].push(Some((strands, snap, welds)));
        }
    }
    ogeom_core::progress::stage("boolean: split arrange");
    let same_key = |a: &(usize, usize, (f64, f64)), b: &(usize, usize, (f64, f64))| {
        a.0 == b.0
            && a.1 == b.1
            && (a.2.0 - b.2.0).abs() <= tol.parametric()
            && (a.2.1 - b.2.1).abs() <= tol.parametric()
    };
    // The other solid's boundary, prepared once per side. It is asked once
    // per face piece, and what it costs to prepare (every face's trimming
    // rings, polylined) does not depend on the point being asked about.
    // Rebuilt per question, the preparation costs hundreds of times the
    // question.
    // A sheet bounds nothing, so nothing is classified against one.
    let boundaries = [
        if matches!(operands, Operands::SheetBySheet | Operands::Itself) {
            None
        } else {
            Some(ogeom_algo::SolidBoundary::of_kept(
                model,
                &gb.solid,
                tol.confusion() * 1e4,
                gb.closed,
                tol,
            )?)
        },
        // A solid of one shell whose edges setting faces aside read closed
        // is not asked again.
        if operands == Operands::Solids {
            Some(ogeom_algo::SolidBoundary::of_kept(
                model,
                &ga.solid,
                tol.confusion() * 1e4,
                ga.closed,
                tol,
            )?)
        } else {
            None
        },
    ];
    // Each face is split and its pieces classified on its own: it reads
    // the two solids, the sections, the contacts and the junctions found
    // so far, and writes only its own pieces and the junctions it adds. So
    // the faces run in parallel, and their findings are joined in face
    // order, the order one thread would have made them in.
    type FaceWork = std::sync::Mutex<Prepared>;
    let mut work: Vec<(bool, usize, FaceWork)> = Vec::new();
    for (side, slots) in prepared.iter_mut().enumerate() {
        for (fi, slot) in slots.iter_mut().enumerate() {
            let Some(prepared_face) = slot.take() else {
                ogeom_bail!(Construction, "a face was prepared twice");
            };
            work.push((side == 0, fi, std::sync::Mutex::new(Some(prepared_face))));
        }
    }
    let found = ogeom_core::parallel::map_ordered(
        &work,
        |_, (from_a, fi, prepared_face)| -> OgeomResult<(Vec<FacePiece>, Vec<Junction>)> {
            ogeom_core::progress::checkpoint()?;
            let (from_a, fi) = (*from_a, *fi);
            let (own, other, boundary) = if from_a {
                (&ga, &gb.solid, boundaries[0].as_ref())
            } else {
                (&gb, &ga.solid, boundaries[1].as_ref())
            };
            let face = &own.faces[fi];
            let mut junctions: Vec<Junction> = Vec::new();
            let Some((mut strands, mut face_snap, welds)) = prepared_face
                .lock()
                .map_err(|_| ogeom_core::ogeom_err!(Construction, "a face's work was poisoned"))?
                .take()
            else {
                ogeom_bail!(Construction, "a face was prepared twice");
            };
            junctions.extend(welds);
            if quiet[usize::from(!from_a)][fi] {
                if let Some(piece) = quiet_piece(model, from_a, fi, face, boundary, tol)? {
                    return Ok((vec![piece], Vec::new()));
                }
                // Read on the other solid after all: split it as any face.
                let welds;
                (strands, face_snap, welds) = strands_of(from_a, fi, face)?;
                junctions.extend(welds);
            }
            let mut pieces: Vec<FacePiece> = Vec::new();
            // A strand's ends in space, or none for a pole.
            let space_ends = |face: &GFace, tag: &Tag| -> OgeomResult<Option<(Point, Point)>> {
                Ok(Some(match tag {
                    Tag::Boundary { edge, range } => {
                        let e = &face.edges[*edge];
                        (
                            e.curve.point_at(range.0, tol)?,
                            e.curve.point_at(range.1, tol)?,
                        )
                    }
                    Tag::Contact { contact, range } => {
                        let c = &contacts[*contact];
                        (
                            c.curve.point_at(range.0, tol)?,
                            c.curve.point_at(range.1, tol)?,
                        )
                    }
                    Tag::Section { section, range } => {
                        let sec = &sections[*section];
                        let domain = sec.curve.domain();
                        (
                            sec.curve
                                .point_at(at_param(range.0, domain, sec.closed), tol)?,
                            sec.curve
                                .point_at(at_param(range.1, domain, sec.closed), tol)?,
                        )
                    }
                    Tag::Pole { .. } => return Ok(None),
                }))
            };
            // A strand this face's own weld collapses is one point in space
            // as well: its ends become one junction, or the pieces either
            // side of it end on two vertices the chart says are one, and a
            // loop cut at a seam a few microns from its own start closes in
            // the chart and stays open in space.
            for st in &strands {
                if st.polyline.len() >= 2
                    && chart_length(&st.polyline) <= face_snap
                    && let Some((from, to)) = space_ends(face, &st.tag)?
                    && from.distance(to) > tol.confusion() * 1e2
                {
                    junctions.push(Junction {
                        at: from.midpoint(to),
                        reach: from.distance(to) / 2.0 + tol.confusion() * 1e2,
                        onto_vertex: false,
                    });
                }
            }
            // A piece some other chart already collapsed collapses here too:
            // its ends become one node, every neighbour meeting them moves
            // onto it, and one junction owns the span in space so the
            // rebuilt vertices agree on every face.
            let forced: Vec<usize> = strands
                .iter()
                .enumerate()
                .filter(|(_, st)| {
                    st.polyline.len() >= 2
                        && chart_length(&st.polyline) > face_snap
                        && node_of(face, &st.tag)
                            .is_some_and(|k| dust.iter().any(|d| same_key(d, &k)))
                })
                .map(|(i, _)| i)
                .collect();
            if !forced.is_empty() {
                let mut reps: Vec<(Point2, usize)> = Vec::new();
                let canon = |p: Point2, reps: &mut Vec<(Point2, usize)>| -> usize {
                    match reps.iter().position(|(q, _)| q.distance(p) <= face_snap) {
                        Some(i) => i,
                        None => {
                            reps.push((p, reps.len()));
                            reps.len() - 1
                        }
                    }
                };
                fn root(reps: &mut [(Point2, usize)], mut i: usize) -> usize {
                    while reps[i].1 != i {
                        reps[i].1 = reps[reps[i].1].1;
                        i = reps[i].1;
                    }
                    i
                }
                for &i in &forced {
                    let st = &strands[i];
                    let a = canon(st.polyline[0], &mut reps);
                    let b = canon(st.polyline[st.polyline.len() - 1], &mut reps);
                    let (ra, rb) = (root(&mut reps, a), root(&mut reps, b));
                    if ra != rb {
                        reps[rb].1 = ra;
                    }
                    let Some((from, to)) = space_ends(face, &st.tag)? else {
                        continue;
                    };
                    ogeom_core::diag!(
                        Channel::Strands,
                        "DUST from_a={from_a} fi={fi}: strand {i} collapses with its partners, {from:?} -> {to:?}"
                    );
                    junctions.push(Junction {
                        at: from.midpoint(to),
                        reach: from.distance(to) / 2.0 + tol.confusion() * 1e2,
                        onto_vertex: false,
                    });
                }
                for (i, st) in strands.iter_mut().enumerate() {
                    if forced.contains(&i) {
                        continue;
                    }
                    let last = st.polyline.len() - 1;
                    for end in [0, last] {
                        if let Some(k) = reps
                            .iter()
                            .position(|(q, _)| q.distance(st.polyline[end]) <= face_snap)
                        {
                            let r = root(&mut reps, k);
                            st.polyline[end] = reps[r].0;
                        }
                    }
                }
                let mut index = 0_usize;
                strands.retain(|_| {
                    let keep = !forced.contains(&index);
                    index += 1;
                    keep
                });
            }
            // A hole nothing comes near is walked on its own, and the face's
            // walk is spared it: a face crossed by a small tool can hold
            // hundreds of holes the tool never reaches.
            let lone = lone_holes(model, face, &strands, face_snap)?;
            let every: Vec<usize>;
            let (walked, places, lone): (&[Strand<Tag>], &[usize], &[arrange::Lone<Tag>]) =
                match &lone {
                    Some((rest, places, lone)) => (rest, places, lone),
                    None => {
                        every = (0..strands.len()).collect();
                        (&strands, &every, &[])
                    }
                };
            let split = match arrange_pieces(walked, places, face_snap, lone) {
                Ok(split) => split,
                Err(err) => {
                    if diag::enabled(Channel::Arrange) {
                        ogeom_core::diag!(
                            Channel::Arrange,
                            "ARRANGE from_a={from_a} fi={fi} failed: {err}"
                        );
                        for (si, st) in strands.iter().enumerate() {
                            let (a, b) = (st.polyline[0], st.polyline[st.polyline.len() - 1]);
                            let describe = |curve: &Curve, range: (f64, f64)| -> String {
                                let mid = curve
                                    .point_at(f64::midpoint(range.0, range.1), tol)
                                    .map(|p| format!("({:.4},{:.4},{:.4})", p.x, p.y, p.z))
                                    .unwrap_or_default();
                                match curve {
                                    Curve::Circle(c) => format!(
                                        "circle centre {:?} r {:.4} mid {mid}",
                                        c.circle().centre(),
                                        c.circle().radius()
                                    ),
                                    other => {
                                        format!("{:?} mid {mid}", core::mem::discriminant(other))
                                    }
                                }
                            };
                            let geometry = match st.tag {
                                Tag::Boundary { edge, range } => {
                                    describe(&face.edges[edge].curve, range)
                                }
                                Tag::Contact { contact, range } => {
                                    format!(
                                        "{} window {:?} along {:?} node {}",
                                        describe(&contacts[contact].curve, range),
                                        contacts[contact].crange,
                                        contact_along[contact],
                                        contacts[contact].node.index()
                                    )
                                }
                                _ => String::new(),
                            };
                            ogeom_core::diag!(Channel::Arrange, "    {geometry}");
                            ogeom_core::diag!(
                                Channel::Arrange,
                                "  strand {si} boundary={} {:?} {a:?} .. {b:?}",
                                st.boundary,
                                match st.tag {
                                    Tag::Boundary { edge, range } =>
                                        format!("Boundary e{edge} {range:?}"),
                                    Tag::Contact { contact, range } =>
                                        format!("Contact c{contact} {range:?}"),
                                    Tag::Section { section, range } =>
                                        format!("Section s{section} {range:?}"),
                                    Tag::Pole { pole, range } => format!("Pole p{pole} {range:?}"),
                                }
                            );
                        }
                    }
                    return Err(err);
                }
            };
            for mut piece in split {
                let Some(boundary) = boundary else {
                    let probe = face.surface.point_at(piece.first.x, piece.first.y, tol)?;
                    pieces.push(FacePiece {
                        from_a,
                        face: fi,
                        rings: piece.rings,
                        outlines: piece.outlines,
                        probe,
                        state: PieceState::Unread,
                        covered: false,
                        whole: false,
                    });
                    continue;
                };
                // Where a piece stands is asked at its interior probes in
                // turn. The first is the roomiest, and usually the only one
                // needed. The rest are for the piece that merely *touches*
                // the other solid, whose roomiest probe can land on the
                // contact and read neither in nor out.
                //
                // A partner face is asked before the classifier, not after.
                // Where a piece sits on a surface the other solid also
                // carries, the partner *is* the answer, and getting there
                // through the classifier means every ray grazing the shared
                // face, its whole fan of directions exhausted, before it
                // reports the On the partner list already knew, which on a
                // part whose bore is refilled by its own cylinder is orders
                // of magnitude slower.
                //
                // A probe every ray from which is ambiguous (one a hair off
                // the other solid's faces, where they run along this one)
                // says nothing either way, and the next is asked. Only a
                // piece none of whose probes can be read is refused.
                let partners = if from_a { &same_a[fi] } else { &same_b[fi] };
                let mut chosen = None;
                let mut unread = None;
                // What a probe reads: on a coincident partner, or what the
                // other solid's boundary says.
                let ask = |candidate: Point2| -> OgeomResult<(Point, OgeomResult<Containment>)> {
                    let at = face.surface.point_at(candidate.x, candidate.y, tol)?;
                    let shared = !partners.is_empty()
                        && partners.iter().any(|&pi| {
                            let partner = if from_a { &gb.faces[pi] } else { &ga.faces[pi] };
                            chart_point_within(partner, at, partner_reach(face, partner, tol), tol)
                                .is_some()
                        });
                    if shared {
                        return Ok((at, Ok(Containment::On)));
                    }
                    match boundary.holds(model, at, tol) {
                        Err(e @ ogeom_core::OgeomError::NotDone(_)) => Ok((at, Err(e))),
                        Err(e) => Err(e),
                        says => Ok((at, says)),
                    }
                };
                // The roomiest probe settles most pieces on its own, so it is
                // asked before the rest are found, wherever the rest could not
                // change whether it is asked first: on a face that is not a
                // plane every probe is kept, and on a plane bounded by lines
                // and circles its exact rings' count for it does not depend on
                // the rings drawn for the rest. Where it reads neither in nor
                // out, every probe is asked in turn below, as if it had not
                // been.
                let first_kept = match face.surface {
                    SurfaceGeometry::Plane(_) => {
                        let drawn = piece.rings.iter().flatten().all(|t| {
                            fine_strand(face, &sections, from_a, &t.tag).is_some_and(|s| {
                                matches!(s.0, PlanarCurve::Line(_) | PlanarCurve::Circle(_))
                            })
                        });
                        drawn
                            .then(|| {
                                FineRings::of(face, &piece.rings, &sections, from_a, None, tol)
                            })
                            .flatten()
                            .and_then(|every| {
                                let reached = rings_near(&piece.outlines, &[piece.first]);
                                every.inside_beside(piece.first, &reached)
                            })
                    }
                    // On a curved face, a piece on a partner's surface has
                    // its roomiest probe held to its rings drawn finely: a
                    // thin piece's outline chords bow past the piece, and a
                    // probe beyond them reads the other solid off the
                    // partner, as if the piece lay there.
                    // Where the rings cannot be drawn finely (a contact
                    // strand among them), a probe off the face itself is
                    // not the piece's either.
                    _ if !partners.is_empty() => {
                        match FineRings::of_curved(face, &piece.rings, &sections, from_a, tol) {
                            Some(every) => Some(every.inside(piece.first)),
                            None => Some(on_own_face(face, piece.first, tol)?),
                        }
                    }
                    _ => Some(true),
                };
                if first_kept == Some(true)
                    && let (at, Ok(says @ (Containment::In | Containment::Out))) = ask(piece.first)?
                {
                    chosen = Some((piece.first, at, says));
                }
                // A planar piece's probes are held to its exact rings, and
                // where none of them is inside those, the rings offer their
                // own: a sliver narrower than the outline's chords bow has
                // every coarse probe in a neighbour.
                // Only the rings a probe's ray can cross are drawn finely: on a
                // face with hundreds of holes the rest are far from every probe.
                let mut fine = None;
                if chosen.is_none() {
                    let mut interiors: Vec<Point2> = piece.interiors().to_vec();
                    let reached = rings_near(&piece.outlines, &interiors);
                    fine =
                        FineRings::of(face, &piece.rings, &sections, from_a, Some(&reached), tol)
                            .or_else(|| {
                                (!partners.is_empty())
                                    .then(|| {
                                        FineRings::of_curved(
                                            face,
                                            &piece.rings,
                                            &sections,
                                            from_a,
                                            tol,
                                        )
                                    })
                                    .flatten()
                            });
                    if fine.is_none() && !partners.is_empty() {
                        let mut kept = Vec::with_capacity(interiors.len());
                        for p in &interiors {
                            if on_own_face(face, *p, tol)? {
                                kept.push(*p);
                            }
                        }
                        // Where none stands on the face as read, the reading
                        // is not trusted over the piece's own outline.
                        if !kept.is_empty() {
                            interiors = kept;
                        }
                    }
                    if let Some(fine) = &fine {
                        interiors.retain(|p| fine.inside(*p));
                        if interiors.is_empty()
                            && let Some(every) =
                                FineRings::of(face, &piece.rings, &sections, from_a, None, tol)
                                    .or_else(|| {
                                        FineRings::of_curved(
                                            face,
                                            &piece.rings,
                                            &sections,
                                            from_a,
                                            tol,
                                        )
                                    })
                        {
                            interiors =
                                arrange::interior_points_of(&every.points(), tol.parametric());
                        }
                    }
                    for candidate in &interiors {
                        let says = match ask(*candidate)? {
                            (at, Ok(says)) => (at, says),
                            (_, Err(e)) => {
                                unread = Some(e);
                                continue;
                            }
                        };
                        if chosen.is_none() || !matches!(says.1, Containment::On) {
                            chosen = Some((*candidate, says.0, says.1));
                        }
                        if !matches!(says.1, Containment::On) {
                            break;
                        }
                    }
                }
                // On the attempt after kept pieces failed to close, a piece
                // is read at the centre of its area where that lies inside
                // it: the roomiest probe can stand by a corner, and a face
                // leaning microradians off the other solid's, the two
                // sharing an edge, reads within rounding of it there.
                if read_at_centre
                    && matches!(chosen, Some((_, _, Containment::In | Containment::Out)))
                    && let Some(centre) = arrange::centroid_of(&piece.outlines)
                    && arrange::inside_rings(&piece.outlines, centre)
                {
                    if fine.is_none() {
                        let reached = rings_near(&piece.outlines, &[centre]);
                        fine = FineRings::of(
                            face,
                            &piece.rings,
                            &sections,
                            from_a,
                            Some(&reached),
                            tol,
                        );
                    }
                    if fine.as_ref().is_none_or(|f| f.inside(centre)) {
                        let at = face.surface.point_at(centre.x, centre.y, tol)?;
                        if let Ok(says @ (Containment::In | Containment::Out)) =
                            boundary.holds(model, at, tol)
                        {
                            chosen = Some((centre, at, says));
                        }
                    }
                }
                let Some((interior, probe, said)) = chosen else {
                    if let Some(e) = unread {
                        return Err(e);
                    }
                    ogeom_bail!(
                        Construction,
                        "a piece of a face has no interior point to classify at"
                    );
                };
                // On the last attempt, a piece thinner than the band the
                // other solid is read within, every probe of it on that
                // solid's boundary, is read by the strands bounding it.
                let thin = if read_thin && matches!(said, Containment::On) {
                    let theirs = if from_a { &gb } else { &ga };
                    thin_by_strands(face, &piece.rings, &piece.outlines, &contacts, theirs, tol)?
                } else {
                    Thin::Wide
                };
                match thin {
                    Thin::Reads(state) => {
                        pieces.push(FacePiece {
                            from_a,
                            face: fi,
                            rings: piece.rings,
                            outlines: piece.outlines,
                            probe,
                            state,
                            covered: false,
                            whole: false,
                        });
                        continue;
                    }
                    Thin::Disagree => ogeom_bail!(
                        NotDone,
                        "a piece thinner than the boolean reads the other solid \
                         within lies beside that solid's edges on sides that \
                         disagree on whether it is inside: it is dust the \
                         boolean does not resolve"
                    ),
                    Thin::Wide | Thin::Undecided => {}
                }
                let state = match said {
                    Containment::In => PieceState::In,
                    Containment::Out => PieceState::Out,
                    Containment::On => {
                        // On the other boundary: same-domain contact. The
                        // partner face on the shared surface decides whether
                        // the two materials lie on the same side or oppose.
                        let own_normal = outward_normal(face, interior, tol)?;
                        let mut resolved = None;
                        for &pi in partners {
                            let partner = if from_a { &gb.faces[pi] } else { &ga.faces[pi] };
                            let Some(uv) = chart_point_within(
                                partner,
                                probe,
                                partner_reach(face, partner, tol),
                                tol,
                            ) else {
                                continue;
                            };
                            let theirs = outward_normal(partner, uv, tol)?;
                            resolved = Some(if own_normal.dot(theirs) > 0.0 {
                                PieceState::OnAligned
                            } else {
                                PieceState::OnOpposed
                            });
                            break;
                        }
                        let Some(state) = resolved else {
                            // On with no partner containing the probe: the
                            // band's generosity may have read proximity as
                            // coincidence. Ask again at a width where it
                            // cannot. What stays on the boundary there is a
                            // genuine contact, settled below.
                            match ogeom_algo::classify_in_solid_exact_banded(
                                model,
                                other,
                                probe,
                                tol.confusion() * 10.0,
                                tol,
                            )? {
                                Containment::In => {
                                    pieces.push(FacePiece {
                                        from_a,
                                        face: fi,
                                        rings: piece.rings,
                                        outlines: piece.outlines,
                                        probe,
                                        state: PieceState::In,
                                        covered: false,
                                        whole: false,
                                    });
                                    continue;
                                }
                                Containment::Out => {
                                    pieces.push(FacePiece {
                                        from_a,
                                        face: fi,
                                        rings: piece.rings,
                                        outlines: piece.outlines,
                                        probe,
                                        state: PieceState::Out,
                                        covered: false,
                                        whole: false,
                                    });
                                    continue;
                                }
                                Containment::On => {}
                            }
                            // Every probe the piece offered stands on the
                            // other solid's boundary, and none on a coincident
                            // face. Where the two touch only along lines or at
                            // points the piece does not cross the other
                            // boundary there (a crossing would have split it),
                            // so it lies wholly on one side, and any point of
                            // it off the contact says which. Points at golden
                            // fractions of the piece's extent miss a contact
                            // that met its halves and quarters. Where they
                            // read on the boundary too, the contact covers a
                            // region, which is not this case.
                            let mut off_contact = None;
                            for candidate in
                                arrange::off_contact_points(&piece.outlines, tol.parametric())
                            {
                                if !arrange::inside_rings(&piece.outlines, candidate)
                                    || fine.as_ref().is_some_and(|f| !f.inside(candidate))
                                {
                                    continue;
                                }
                                let at = face.surface.point_at(candidate.x, candidate.y, tol)?;
                                match boundary.holds(model, at, tol) {
                                    Ok(Containment::In) => {
                                        off_contact = Some(PieceState::In);
                                        break;
                                    }
                                    Ok(Containment::Out) => {
                                        off_contact = Some(PieceState::Out);
                                        break;
                                    }
                                    Ok(Containment::On)
                                    | Err(ogeom_core::OgeomError::NotDone(_)) => {}
                                    Err(e) => return Err(e),
                                }
                            }
                            if let Some(state) = off_contact {
                                pieces.push(FacePiece {
                                    from_a,
                                    face: fi,
                                    rings: piece.rings,
                                    outlines: piece.outlines,
                                    probe,
                                    state,
                                    covered: false,
                                    whole: false,
                                });
                                continue;
                            }
                            // Still on it, on the boolean's second attempt
                            // (see [`settled`]): a face all but coplanar with
                            // the other's, crossing it at a hair's angle, keeps
                            // a sliver within tolerance of it with no
                            // coincident partner. Just off the piece on both
                            // sides the other solid answers plainly: the
                            // same on both, the piece is that. Inside only
                            // on the piece's own material side, the other's
                            // face runs along it backed the same way, the
                            // aligned contact, which the first solid's piece
                            // speaks for. Inside only beyond it, opposed.
                            if settle {
                                let off = own_normal * (tol.confusion() * 1e3);
                                let side = |q: Point| {
                                    ogeom_algo::classify_in_solid_exact_banded(
                                        model,
                                        other,
                                        q,
                                        tol.confusion() * 10.0,
                                        tol,
                                    )
                                };
                                let (beyond, within) = (side(probe + off)?, side(probe - off)?);
                                let read = match (beyond, within) {
                                    (Containment::In, Containment::In) => {
                                        Some((PieceState::In, false))
                                    }
                                    (Containment::Out, Containment::Out) => {
                                        Some((PieceState::Out, false))
                                    }
                                    (Containment::Out, Containment::In) => {
                                        Some((PieceState::OnAligned, !from_a))
                                    }
                                    (Containment::In, Containment::Out) => {
                                        Some((PieceState::OnOpposed, false))
                                    }
                                    _ => None,
                                };
                                if let Some((state, covered)) = read {
                                    pieces.push(FacePiece {
                                        from_a,
                                        face: fi,
                                        rings: piece.rings,
                                        outlines: piece.outlines,
                                        probe,
                                        state,
                                        covered,
                                        whole: false,
                                    });
                                    continue;
                                }
                            }
                            ogeom_core::diag!(
                                Channel::Wire,
                                "ON-CONTACT piece of {} face {fi} probe {probe:?} partners {partners:?}",
                                if from_a { "A" } else { "B" }
                            );
                            ogeom_bail!(
                                NotDone,
                                "a piece lies on the other solid's boundary \
                             at every point it was asked at, with no \
                             coincident partner face to compare sides \
                             against; a contact over a region of faces not \
                             recognised as coincident is refused; see the \
                             open items in docs/PLAN.md"
                            );
                        };
                        state
                    }
                };
                pieces.push(FacePiece {
                    from_a,
                    face: fi,
                    rings: piece.rings,
                    outlines: piece.outlines,
                    probe,
                    state,
                    covered: false,
                    whole: false,
                });
            }
            // A face nothing crosses, its boundary unsplit, is one piece:
            // itself.
            let untouched = strands.iter().all(|st| match &st.tag {
                Tag::Boundary { edge, range } => {
                    let full = face.edges[*edge].crange;
                    (range.0.min(range.1) - full.0.min(full.1)).abs() <= tol.parametric()
                        && (range.0.max(range.1) - full.0.max(full.1)).abs() <= tol.parametric()
                }
                Tag::Pole { .. } => true,
                Tag::Section { .. } | Tag::Contact { .. } => false,
            });
            if untouched
                && let [piece] = pieces.as_mut_slice()
                && matches!(
                    piece.state,
                    PieceState::In | PieceState::Out | PieceState::Unread
                )
            {
                piece.whole = true;
            }

            Ok((pieces, junctions))
        },
    );
    for face in found {
        let (face_pieces, face_junctions) = face?;
        pieces.extend(face_pieces);
        junctions.extend(face_junctions);
    }
    let coincident = mark_covered_coincidences(&ga, &gb, &mut pieces, tol);
    if diag::enabled(Channel::Arrange) {
        for (fi, partners) in same_a.iter().enumerate() {
            if !partners.is_empty() {
                ogeom_core::diag!(Channel::Arrange, "SAME a{fi} with b{partners:?}");
            }
        }
        for (ci, contact) in contacts.iter().enumerate() {
            ogeom_core::diag!(
                Channel::Arrange,
                "CONTACT c{ci} onto {}{} over {:?} prange {:?} node {} pcurve {:?}",
                if contact.target_from_a { "a" } else { "b" },
                contact.target_face,
                contact.crange,
                contact.prange,
                contact.node.index(),
                contact.pcurve
            );
        }
        for (i, p) in pieces.iter().enumerate() {
            let own = if p.from_a {
                &ga.faces[p.face]
            } else {
                &gb.faces[p.face]
            };
            let kind = match &own.surface {
                SurfaceGeometry::Plane(_) => "plane",
                SurfaceGeometry::Cylinder(_) => "cyl",
                SurfaceGeometry::Torus(_) => "torus",
                SurfaceGeometry::BSpline(_) => "bspline",
                _ => "other",
            };
            ogeom_core::diag!(
                Channel::Arrange,
                "PIECE {i} from_a={} face={} {kind} state={:?} covered={} probe={:?} rings {:?}",
                p.from_a,
                p.face,
                p.state,
                p.covered,
                p.probe,
                p.rings.iter().map(Vec::len).collect::<Vec<_>>()
            );
        }
    }
    ogeom_core::progress::stage("boolean: classified");
    // Merged again: the arrangement adds junctions of its own (a strand
    // that collapsed with its partners becomes one) after the first
    // merge, and two of those at one triple point, each within the other's
    // reach, weld the ends of one band to two vertices with a hairline
    // between them.
    let junctions = merge_junctions(junctions);
    Ok(GeneralFused {
        a: ga,
        b: gb,
        sections,
        contacts,
        tangents,
        pieces,
        coincident,
        junctions,
    })
}
