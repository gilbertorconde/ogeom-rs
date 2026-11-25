//! The general fuse's pieces: junctions, outlines, thin pieces, hugs and
//! coincidences, and the settling of pieces on the other boundary.

use super::*;

// --- the general fuse --------------------------------------------------------

/// One piece of one argument face, classified against the other argument.
pub(crate) struct FacePiece {
    /// Which argument's face list, and which face.
    pub(crate) from_a: bool,
    pub(crate) face: usize,
    pub(crate) rings: Vec<Vec<Traversal<Tag>>>,
    /// The rings as chart polylines, for the coincidence pairing below.
    pub(crate) outlines: Vec<Vec<Point2>>,
    /// The interior point the state was decided at, in space.
    pub(crate) probe: Point,
    pub(crate) state: PieceState,
    /// Whether the *other* argument already contributes this same patch of
    /// this same surface.
    ///
    /// Two faces lying on one surface with their material on the same side
    /// bound the union and the intersection once between them, not twice,
    /// but "once" is a statement about the *patch*, not about which argument
    /// it came from. Where a tool's cap fills a hole its own bore left, the
    /// part has no piece there at all, and dropping the tool's by argument
    /// identity leaves the result with a hole in it. So the duplicate is
    /// identified by containment, and only a piece genuinely stood in for is
    /// dropped.
    pub(crate) covered: bool,
    /// Whether the piece is its whole face, untouched: nothing crosses it
    /// and none of its edges is split.
    pub(crate) whole: bool,
}

pub(crate) struct GeneralFused {
    pub(crate) a: GSolid,
    pub(crate) b: GSolid,
    pub(crate) sections: Vec<SectionRec>,
    pub(crate) contacts: Vec<ContactRec>,
    /// Curves the two boundaries touch along without crossing. They take no
    /// part in the classification below (that is what tangency means) and
    /// are carried for the consumers that want the contact itself.
    pub(crate) tangents: Vec<TangentRec>,
    pub(crate) pieces: Vec<FacePiece>,
    /// Pieces of the two arguments standing for one patch of one surface,
    /// as indices into `pieces`: the first argument's piece, then the
    /// second's.
    pub(crate) coincident: Vec<(usize, usize)>,
    /// Junctions several paves describe, each resolved once for every
    /// strand that ends in it.
    pub(crate) junctions: Vec<Junction>,
}

/// One junction the paving found several times over.
///
/// A tolerant rail meeting a wedge's faces near one corner collects a
/// cluster of crossings inside its own stated radius: each face's section
/// stops at its own crossing with the rail, and the crossings sit a few
/// tenths of a micron to a few hundred apart along it. The rail's strands
/// split once per cluster, at its first pave. The sections' ends still name
/// their own crossings. And the two descriptions of the junction can sit
/// apart by the cluster's whole span, which no single strand's honesty
/// covers. So the junction is one vertex standing at the first pave, owning
/// the span the paves disagree by, and every strand end inside that span
/// names it, whichever pave its own trim stopped at.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Junction {
    /// Where the cluster's first pave sits.
    pub(crate) at: Point,
    /// How far a strand end may sit from `at` and still be this junction:
    /// the cluster's span plus the rail's own reach.
    pub(crate) reach: f64,
    /// Whether `at` is an edge's own end vertex that a section tangent to
    /// the edge was welded onto: the section's chart image is bent onto it
    /// rather than left where the section stopped.
    pub(crate) onto_vertex: bool,
}

/// Paves on each seam of `solid`'s faces where a junction the seam passes
/// through holds section ends of the face on both sides of it.
///
/// The ends of one junction are one point, and a strand end is carried onto
/// a node the junction holds, or onto another end, across the face's chart.
/// Ends either side of a seam stand near its two columns, a period apart in
/// the chart, and no path inside the chart joins them: with no node on the
/// seam they dangle on either side of it. Split where it passes nearest the
/// junction, the seam puts a node on each column, and each end is carried
/// onto the one on its side, as ends a junction holds on one side of the
/// seam are carried onto one another.
pub(crate) fn seams_through_junctions(
    solid: &GSolid,
    from_a: bool,
    (sections, pieces): (&[SectionRec], &[SectionPiece]),
    junctions: &[Junction],
    paves: &mut ogeom_core::FastMap<EdgeKey, Vec<Pave>>,
    tol: Tolerances,
) -> OgeomResult<()> {
    for (fi, face) in solid.faces.iter().enumerate() {
        let seams: Vec<&BoundaryEdge> = face
            .edges
            .iter()
            .filter(|e| e.other_side.is_some())
            .collect();
        if seams.is_empty() {
            continue;
        }
        // Each end of the face's section pieces, in space and in the chart.
        let mut ends: Vec<(Point, Point2)> = Vec::new();
        for sp in pieces {
            let section = &sections[sp.section];
            let (face_of, pcurve, hugs) = if from_a {
                (section.face_a, &section.pc_a, sp.hugs[0])
            } else {
                (section.face_b, &section.pc_b, sp.hugs[1])
            };
            if face_of != fi || hugs || section.closed {
                continue;
            }
            for t in [sp.range.0, sp.range.1] {
                ends.push((section.curve.point_at(t, tol)?, pcurve.point_at(t, tol)?));
            }
        }
        for e in seams {
            let Some((other, orange)) = &e.other_side else {
                continue;
            };
            for junction in junctions {
                if !e.bound.expanded(junction.reach).contains(junction.at) {
                    continue;
                }
                let (gap, t) = nearest_on_edge_curve(&e.curve, e.crange, junction.at, tol)?;
                if gap > junction.reach {
                    continue;
                }
                let paved = paves.get(&e.node).is_some_and(|list| {
                    list.iter().any(|p| {
                        e.curve
                            .point_at(p.t, tol)
                            .is_ok_and(|q| junction.at.distance(q) <= junction.reach)
                    })
                });
                if paved {
                    continue;
                }
                let columns = [
                    e.pcurve.point_at(rescale(t, e.crange, e.prange), tol)?,
                    other.point_at(rescale(t, e.crange, *orange), tol)?,
                ];
                let mut sides = [false; 2];
                for (at, chart) in &ends {
                    if junction.at.distance(*at) <= junction.reach {
                        let nearer =
                            usize::from(chart.distance(columns[1]) < chart.distance(columns[0]));
                        sides[nearer] = true;
                    }
                }
                if sides == [true; 2] {
                    ogeom_core::diag!(
                        Channel::Wire,
                        "PAVE seam e{} at {t:.6}: a junction at {:?} holds ends on both sides",
                        e.node.index(),
                        junction.at
                    );
                    paves.entry(e.node).or_default().push(Pave {
                        t,
                        honesty: gap,
                        across: None,
                        spread: 0.0,
                    });
                }
            }
        }
    }
    Ok(())
}

/// Junctions whose balls overlap, merged transitively into one junction
/// each: at the members' centroid, reaching as far as the farthest member's
/// own reach extends from it.
pub(crate) fn merge_junctions(junctions: Vec<Junction>) -> Vec<Junction> {
    let n = junctions.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    // Swept along x: two balls overlap only if their centres are within
    // both reaches along every axis, so each junction is asked only of those
    // ahead of it by no more than its reach and the widest. The groups are
    // the same as asking every pair, which for the tens of thousands of
    // junctions a converted part carries costs seconds.
    let widest = junctions.iter().map(|j| j.reach).fold(0.0_f64, f64::max);
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&x, &y| junctions[x].at.x.total_cmp(&junctions[y].at.x));
    for (k, &i) in order.iter().enumerate() {
        let limit = junctions[i].at.x + junctions[i].reach + widest;
        for &j in &order[k + 1..] {
            if junctions[j].at.x > limit {
                break;
            }
            if junctions[i].at.distance(junctions[j].at) <= junctions[i].reach + junctions[j].reach
            {
                let (ri, rj) = (root(&mut parent, i), root(&mut parent, j));
                if ri != rj {
                    parent[rj] = ri;
                }
            }
        }
    }
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut of_root: Vec<Option<usize>> = vec![None; n];
    for i in 0..n {
        let r = root(&mut parent, i);
        match of_root[r] {
            Some(g) => groups[g].push(i),
            None => {
                of_root[r] = Some(groups.len());
                groups.push(vec![i]);
            }
        }
    }
    // A group is one junction only while it stays a junction's size: a
    // chain of overlapping balls along an edge would otherwise merge into
    // one reaching the chain's length, and the rebuild would weld the edge
    // into a point. A group reaching beyond four of its widest member stays
    // as it was.
    let mut out = Vec::with_capacity(groups.len());
    for members in groups {
        if members.len() == 1 {
            out.push(junctions[members[0]]);
            continue;
        }
        let mut sum = ogeom_math::Vector::ZERO;
        for &m in &members {
            sum += junctions[m].at - Point::ORIGIN;
        }
        #[allow(clippy::cast_precision_loss)]
        let at = Point::ORIGIN + sum / (members.len() as f64);
        let reach = members
            .iter()
            .map(|&m| at.distance(junctions[m].at) + junctions[m].reach)
            .fold(0.0_f64, f64::max);
        let widest = members
            .iter()
            .map(|&m| junctions[m].reach)
            .fold(0.0_f64, f64::max);
        if reach > widest * 4.0 {
            out.extend(members.iter().map(|&m| junctions[m]));
        } else {
            out.push(Junction {
                at,
                reach,
                onto_vertex: members.iter().any(|&m| junctions[m].onto_vertex),
            });
        }
    }
    out
}

/// The junctions the paves describe more than once, or more loosely than
/// the edge itself, per edge in face order.
pub(crate) fn pave_junctions(
    ga: &GSolid,
    gb: &GSolid,
    paves: &ogeom_core::FastMap<EdgeKey, Vec<Pave>>,
    along: &AlongEdges,
    tol: Tolerances,
) -> OgeomResult<Vec<Junction>> {
    let mut seen: ogeom_core::FastSet<EdgeKey> = ogeom_core::FastSet::default();
    let mut junctions = Vec::new();
    for e in ga
        .faces
        .iter()
        .chain(gb.faces.iter())
        .flat_map(|f| f.edges.iter())
    {
        if !seen.insert(e.node) {
            continue;
        }
        let Some(ts) = paves.get(&e.node) else {
            continue;
        };
        // A tolerant edge (a merged facet group's boundary, a fitted rail)
        // is met by a section up to its own radius off it, so the section's
        // end and the edge's split point are one junction that far apart,
        // on either side of the edge.
        let tolerant = e.tolerance > tol.confusion() * 1e2;
        let floor = if tolerant {
            e.tolerance * 2.0
        } else {
            e.tolerance
        }
        .max(tol.confusion() * 10.0);
        let off = |face, p| off_face(ga, gb, face, p, tol);
        let along = along.get(&e.node).map_or(&[][..], Vec::as_slice);
        for cluster in cluster_paves(&e.curve, e.crange, e.tolerance, ts, Some(&off), along, tol)? {
            if cluster.members > 1 || cluster.honesty > tol.confusion() * 1e2 || tolerant {
                junctions.push(Junction {
                    at: cluster.at,
                    reach: cluster.span + cluster.honesty.max(floor),
                    onto_vertex: false,
                });
            }
        }
    }
    Ok(junctions)
}

/// How long an edge of a face left alone must be: a shorter one may be
/// dust some other chart collapses, which only the face's own split finds.
pub(crate) const QUIET_EDGE: f64 = 1e-2;

/// The one piece of a face the other solid leaves alone: its own wires as
/// rings, its state read at the middle of one of its edges, which no
/// section touches and which therefore stands inside the other solid or
/// out of it with the whole face, or [`PieceState::Unread`] where there is no
/// `boundary` to read it against. `None` where that edge reads on the other
/// boundary after all, or a wire's edge is not one the face was gathered
/// with.
pub(crate) fn quiet_piece(
    model: &Model,
    from_a: bool,
    fi: usize,
    face: &GFace,
    boundary: Option<&ogeom_algo::SolidBoundary>,
    tol: Tolerances,
) -> OgeomResult<Option<FacePiece>> {
    let forward = face.face.oriented(ogeom_topo::Orientation::Forward);
    // A face with many holes is looked up by key, not walked once per edge.
    let indexed: Option<ogeom_core::FastMap<EdgeKey, usize>> = (face.edges.len() > 16).then(|| {
        let mut first =
            ogeom_core::FastMap::with_capacity_and_hasher(face.edges.len(), Default::default());
        for (index, e) in face.edges.iter().enumerate() {
            first.entry(e.node).or_insert(index);
        }
        first
    });
    let aside: ogeom_core::FastSet<ogeom_topo::TShapeId> = face
        .holes_aside
        .iter()
        .map(|hole| hole.wire.node())
        .collect();
    let mut rings = Vec::new();
    for wire in model.ordered_children_of(&forward)? {
        if aside.contains(&wire.node()) {
            continue;
        }
        let mut ring = Vec::new();
        for edge in model.ordered_children_of(&wire)? {
            let key = EdgeKey::of(&edge);
            if face.slits.contains(&key) {
                continue;
            }
            let found = match &indexed {
                Some(first) => first.get(&key).copied(),
                None => face.edges.iter().position(|e| e.node == key),
            };
            let Some(index) = found else {
                return Ok(None);
            };
            ring.push(Traversal {
                tag: Tag::Boundary {
                    edge: index,
                    range: face.edges[index].crange,
                },
                reversed: edge.orientation() == ogeom_topo::Orientation::Reversed,
                nodes: None,
            });
        }
        rings.push(ring);
    }
    let Some(e) = face.edges.first() else {
        return Ok(None);
    };
    let probe = e
        .curve
        .point_at(f64::midpoint(e.crange.0, e.crange.1), tol)?;
    let Some(boundary) = boundary else {
        return Ok(Some(FacePiece {
            from_a,
            face: fi,
            rings,
            outlines: Vec::new(),
            probe,
            state: PieceState::Unread,
            covered: false,
            whole: true,
        }));
    };
    let state = match boundary.holds(model, probe, tol) {
        Ok(Containment::In) => PieceState::In,
        Ok(Containment::Out) => PieceState::Out,
        Ok(Containment::On) | Err(ogeom_core::OgeomError::NotDone(_)) => return Ok(None),
        Err(e) => return Err(e),
    };
    Ok(Some(FacePiece {
        from_a,
        face: fi,
        rings,
        outlines: Vec::new(),
        probe,
        state,
        covered: false,
        whole: true,
    }))
}

/// A face's strands with its lone holes taken out: the rest, the place of
/// each among all of them, and the holes, each walked alone.
pub(crate) type Apart = (Vec<Strand<Tag>>, Vec<usize>, Vec<arrange::Lone<Tag>>);

/// The holes of a planar face whose strands stand apart from every other
/// strand, each walked alone (see [`arrange::Lone`]), with the rest of the
/// strands and their places among all of them; `None` where no hole does.
///
/// A hole stands apart when no strand of another wire ends within a few
/// snaps of its box, no section, contact or pole strand's box comes that
/// near it, and its box lies inside the outer wire's.
pub(crate) fn lone_holes(
    model: &Model,
    face: &GFace,
    strands: &[Strand<Tag>],
    snap: f64,
) -> OgeomResult<Option<Apart>> {
    // Fewer strands than this are walked together for no cost worth
    // sparing.
    const FEW: usize = 64;
    if strands.len() < FEW
        || !matches!(face.surface, SurfaceGeometry::Plane(_))
        || !face.poles.is_empty()
    {
        return Ok(None);
    }
    // Each edge's wire, `None` for an edge in no wire or in two.
    let index: ogeom_core::FastMap<EdgeKey, usize> = face
        .edges
        .iter()
        .enumerate()
        .map(|(i, e)| (e.node, i))
        .collect();
    let mut wire_of: Vec<Option<Option<usize>>> = vec![None; face.edges.len()];
    let wires = model.ordered_children_of(&face.face.oriented(ogeom_topo::Orientation::Forward))?;
    for (w, wire) in wires.iter().enumerate() {
        for edge in model.ordered_children_of(wire)? {
            if let Some(&i) = index.get(&EdgeKey::of(&edge)) {
                wire_of[i] = match wire_of[i] {
                    None => Some(Some(w)),
                    Some(_) => Some(None),
                };
            }
        }
    }
    let strand_wire: Vec<Option<usize>> = strands
        .iter()
        .map(|st| match st.tag {
            Tag::Boundary { edge, .. } if st.boundary => {
                wire_of.get(edge).copied().flatten().flatten()
            }
            _ => None,
        })
        .collect();
    let mut by_wire: Vec<Vec<usize>> = vec![Vec::new(); wires.len()];
    for (k, w) in strand_wire.iter().enumerate() {
        if let Some(w) = w {
            by_wire[*w].push(k);
        }
    }
    let boxed = |ks: &[usize]| -> (Point2, Point2) {
        let mut lo = Point2::new(f64::INFINITY, f64::INFINITY);
        let mut hi = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
        for &k in ks {
            for p in &strands[k].polyline {
                lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
                hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
            }
        }
        (lo, hi)
    };
    let Some(outer) = by_wire.first().filter(|ks| !ks.is_empty()) else {
        return Ok(None);
    };
    let (outer_lo, outer_hi) = boxed(outer);
    // Every strand end, by abscissa; and the box of every strand no wire
    // owns.
    let mut ends: Vec<(f64, f64, usize)> = Vec::with_capacity(strands.len() * 2);
    let mut loose: Vec<(Point2, Point2)> = Vec::new();
    for (k, st) in strands.iter().enumerate() {
        for p in [st.polyline.first(), st.polyline.last()]
            .into_iter()
            .flatten()
        {
            ends.push((p.x, p.y, k));
        }
        if strand_wire[k].is_none() {
            loose.push(boxed(&[k]));
        }
    }
    ends.sort_by(|a, b| a.0.total_cmp(&b.0));
    let margin = snap * 4.0;
    let mut lone = Vec::new();
    let mut taken = vec![false; strands.len()];
    for (w, ks) in by_wire.iter().enumerate().skip(1) {
        if ks.is_empty() {
            continue;
        }
        let (lo, hi) = boxed(ks);
        let (lo, hi) = (
            Point2::new(lo.x - margin, lo.y - margin),
            Point2::new(hi.x + margin, hi.y + margin),
        );
        if !(lo.x > outer_lo.x && lo.y > outer_lo.y && hi.x < outer_hi.x && hi.y < outer_hi.y) {
            continue;
        }
        let first = ends.partition_point(|e| e.0 < lo.x);
        let crowded = ends[first..]
            .iter()
            .take_while(|e| e.0 <= hi.x)
            .any(|&(_, y, k)| y >= lo.y && y <= hi.y && strand_wire[k] != Some(w));
        let crossed = loose
            .iter()
            .any(|(a, b)| a.x <= hi.x && b.x >= lo.x && a.y <= hi.y && b.y >= lo.y);
        if crowded || crossed {
            continue;
        }
        let mine: Vec<Strand<Tag>> = ks.iter().map(|&k| strands[k].clone()).collect();
        let Some(ring) = arrange::lone_ring(&mine, ks, snap) else {
            continue;
        };
        for &k in ks {
            taken[k] = true;
        }
        lone.push(ring);
    }
    if lone.is_empty() {
        return Ok(None);
    }
    let (rest, places): (Vec<Strand<Tag>>, Vec<usize>) = strands
        .iter()
        .enumerate()
        .filter(|(k, _)| !taken[*k])
        .map(|(k, st)| (st.clone(), k))
        .unzip();
    Ok(Some((rest, places, lone)))
}

/// The face's outward normal at a chart point: the surface's, flipped when
/// the face presents its other side.
pub(crate) fn outward_normal(
    face: &GFace,
    at: Point2,
    tol: Tolerances,
) -> OgeomResult<ogeom_math::Vector> {
    let n = face.surface.normal_at(at.x, at.y, tol)?.vector();
    Ok(
        if face.face.orientation() == ogeom_topo::Orientation::Reversed {
            -n
        } else {
            n
        },
    )
}

/// How far apart two of a face's edges may honestly end in its chart: the
/// loosest edge's tolerance, or the hug width where an edge is a fitted
/// section that once ran along another before parting from it.
pub(crate) fn outline_snap(face: &GFace, tol: Tolerances) -> f64 {
    face.edges
        .iter()
        .fold(tol.confusion() * 1e2, |acc, e| acc.max(e.tolerance * 2.0))
        .max(
            if face.edges.iter().any(|e| e.tolerance > tol.confusion()) {
                tol.confusion() * 1e3
            } else {
                0.0
            },
        )
}

/// Polyline ends sorted by their chart `x`, for finding the ends within a
/// reach of a point without comparing every pair.
pub(crate) struct EndsByX {
    /// Indices into the ends, by ascending `x`.
    order: Vec<usize>,
    xs: Vec<f64>,
}

impl EndsByX {
    fn new(ends: &[(usize, bool, Point2)]) -> Self {
        let mut order: Vec<usize> = (0..ends.len()).collect();
        order.sort_by(|&a, &b| ends[a].2.x.total_cmp(&ends[b].2.x));
        let xs = order.iter().map(|&k| ends[k].2.x).collect();
        Self { order, xs }
    }

    /// The ends whose `x` lies within `reach` of `p`'s, in ascending order
    /// of their index: a superset of the ends within `reach` of `p`.
    fn near(&self, p: Point2, reach: f64) -> Vec<usize> {
        let lo = self.xs.partition_point(|&x| x < p.x - reach);
        let hi = self.xs.partition_point(|&x| x <= p.x + reach);
        let mut out = self.order[lo..hi.max(lo)].to_vec();
        out.sort_unstable();
        out
    }
}

/// Join the outline's ends the weld left apart, each to the unmatched end
/// nearest it where that end is nearest it too and the opening is small
/// against the outline.
///
/// The topology says a face's boundary closes. Where two edges' pcurves end
/// further apart in the chart than their tolerances reach (a fitted blend
/// edge meeting a column on a narrow torus wedge, a pole's row with no edge
/// along it), a ray cast for containment passes through the opening, counts
/// one crossing, and a point level with it reads inside. A short segment
/// across each opening closes the outline as the topology does.
pub(crate) fn bridge_outline_gaps(lines: &mut Vec<Vec<Point2>>) {
    let ends: Vec<(usize, bool, Point2)> = lines
        .iter()
        .enumerate()
        .filter(|(_, line)| line.len() >= 2)
        .flat_map(|(i, line)| [(i, false, line[0]), (i, true, line[line.len() - 1])])
        .collect();
    // Ends that meet within the parametric snap but not exactly are joined
    // where they stand: the weld leaves them, since moving one by rounding
    // noise takes a probe on a seam off it, but a ray level with the sliver
    // between them still slips through.
    let by_x = EndsByX::new(&ends);
    let mut joins: Vec<Vec<Point2>> = Vec::new();
    for (k, &(i, end_i, p)) in ends.iter().enumerate() {
        for m in by_x.near(p, PARAM_SNAP) {
            let (j, end_j, q) = ends[m];
            let d = p.distance(q);
            if m > k && (j, end_j) != (i, end_i) && d > 0.0 && d <= PARAM_SNAP {
                joins.push(vec![p, q]);
            }
        }
    }
    let open: Vec<Point2> = ends
        .iter()
        .filter(|&&(i, end_i, p)| {
            !by_x.near(p, PARAM_SNAP).into_iter().any(|m| {
                let (j, end_j, q) = ends[m];
                (j, end_j) != (i, end_i) && p.distance(q) <= PARAM_SNAP
            })
        })
        .map(|&(_, _, p)| p)
        .collect();
    // Only an opening small against the outline itself: a pole's row
    // spans the chart and is the pole's to bound, not a gap.
    let (mut lo, mut hi) = (
        Point2::new(f64::INFINITY, f64::INFINITY),
        Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for p in lines.iter().flatten() {
        lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
        hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
    }
    let reach = lo.distance(hi) * 1e-3;
    let nearest = |k: usize| {
        (0..open.len()).filter(|&m| m != k).min_by(|&a, &b| {
            open[k]
                .distance(open[a])
                .total_cmp(&open[k].distance(open[b]))
        })
    };
    for k in 0..open.len() {
        if let Some(m) = nearest(k)
            && k < m
            && nearest(m) == Some(k)
            && open[k].distance(open[m]) <= reach
        {
            lines.push(vec![open[k], open[m]]);
        }
    }
    lines.extend(joins);
}

/// Close the gaps between a face's outline polylines.
///
/// The edges are polylined one by one, and where two of them meet at a
/// vertex that owns some tolerance (a fitted section's end welded to the
/// edge it hugged, a few dozen microns off), their polylines stop that far
/// apart. A ray cast for containment slips through such a gap, and a probe
/// standing inside the face within a hug's width of the seam reads as
/// outside. Ends within the face's snap ([`outline_snap`]) of another
/// edge's end are made one point, so the outline is closed exactly as the
/// topology says it is.
///
/// The snap is a length in space, and so is every gap it is held against:
/// two ends meet only where their points on the surface do too. A chart
/// stretched along one direction (a wall a micron tall charted over a
/// tenth of a unit) holds a whole edge within a chart distance the snap
/// would read as a gap. Nor does the weld reach as far as half the face's
/// shortest edge, measured the same way: the two ends of that edge would
/// then read as one point, and the outline would close over nothing. Two
/// ends of one edge are never joined.
pub(crate) fn weld_outline_ends(lines: &mut [Vec<Point2>], face: &GFace, tol: Tolerances) {
    let snap = outline_snap(face, tol);
    if snap <= 0.0 {
        return;
    }
    // Which edge drew each line: a seam draws two, one per column.
    let owner: Vec<usize> = face
        .edges
        .iter()
        .enumerate()
        .flat_map(|(i, e)| std::iter::repeat_n(i, 1 + usize::from(e.other_side.is_some())))
        .collect();
    if owner.len() != lines.len() {
        return;
    }
    let lift = |p: Point2| face.surface.point_at(p.x, p.y, tol).ok();
    // How far apart in space two ends may be welded: the snap, and no
    // further than half the face's shortest edge there. A line with no
    // length in space (a pole's or an apex's row) has no ends to tell
    // apart.
    let mut within = snap;
    for line in lines.iter() {
        let points: Vec<Point> = line.iter().filter_map(|p| lift(*p)).collect();
        let length: f64 = points.windows(2).map(|w| w[0].distance(w[1])).sum();
        if points.len() == line.len() && length > tol.confusion() {
            within = within.min(length / 2.0);
        }
    }
    // Every polyline end, as (line, is its last point, where), and its
    // point on the surface.
    let ends: Vec<(usize, bool, Point2)> = lines
        .iter()
        .enumerate()
        .flat_map(|(i, line)| {
            let first = line.first().map(|p| (i, false, *p));
            let last = line.last().map(|p| (i, true, *p));
            first.into_iter().chain(last)
        })
        .collect();
    let lifted: Vec<Option<Point>> = ends.iter().map(|&(_, _, p)| lift(p)).collect();
    // Each end moves onto the nearest other end within reach that sorts
    // before it, so a matched pair lands on one point rather than trading
    // places. Within reach in the chart, as a gap the snap is read across
    // there, and on the surface.
    let by_x = EndsByX::new(&ends);
    let mut moves: Vec<(usize, bool, Point2)> = Vec::new();
    for (k, &(i, end_i, p)) in ends.iter().enumerate() {
        let Some(at) = lifted[k] else {
            continue;
        };
        let nearest = by_x
            .near(p, snap)
            .into_iter()
            .filter_map(|m| {
                let (j, end_j, q) = ends[m];
                let d = p.distance(q);
                // Ends closer than the parametric snap already meet for
                // every purpose here, and moving one by rounding noise
                // takes a probe that stands exactly on a chart's seam off
                // it.
                if (j, end_j) >= (i, end_i) || owner[j] == owner[i] || d <= PARAM_SNAP || d > snap {
                    return None;
                }
                (lifted[m]?.distance(at) <= within).then_some((d, q))
            })
            .min_by(|a, b| a.0.total_cmp(&b.0));
        if let Some((_, q)) = nearest {
            moves.push((i, end_i, q));
        }
    }
    for (i, end_i, q) in moves {
        let slot = if end_i {
            lines[i].last_mut()
        } else {
            lines[i].first_mut()
        };
        if let Some(p) = slot {
            *p = q;
        }
    }
}

/// A curve's length and the diagonal of the box its points span, both from
/// the same 64 samples: a trace lapping one place many times is far longer
/// than its box, and a loop is about a turn round it.
pub(crate) fn length_and_reach(curve: &Curve, tol: Tolerances) -> OgeomResult<(f64, f64)> {
    use ogeom_geom::Curve3d as _;
    let (lo, hi) = curve.domain();
    let mut length = 0.0_f64;
    let mut last: Option<Point> = None;
    let mut span = ogeom_math::Aabb::default();
    for k in 0..=64 {
        let t = if k == 64 {
            hi
        } else {
            lo + (hi - lo) * f64::from(k) / 64.0
        };
        let p = curve.point_at(t, tol)?;
        if let Some(q) = last {
            length += q.distance(p);
        }
        last = Some(p);
        span = span.with_point(p);
    }
    Ok((length, span.diagonal()))
}

/// Whether a point a face's surface holds to within `margin` lies inside
/// the face's trim, or within `margin` of its boundary as the face draws it.
pub(crate) fn within_face(
    face: &GFace,
    p: Point,
    margin: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    if chart_point_within(face, p, margin.max(tol.confusion() * 1e2), tol).is_some() {
        return Ok(true);
    }
    for e in &face.edges {
        if !e.bound.expanded(margin).contains(p) {
            continue;
        }
        let (drawn, drange, _) = e.drawn();
        if distance_to_edge_curve(drawn, drange, p, tol)? <= margin {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Where a point sits in a face's chart, for a point up to `reach` off the
/// face's surface and inside its trim. `None` elsewhere.
pub(crate) fn chart_point_within(
    face: &GFace,
    p: Point,
    reach: f64,
    tol: Tolerances,
) -> Option<Point2> {
    // Closed-form inversion for the analytic surfaces: the same-domain
    // resolution asks "where does this probe sit in the partner's chart", and
    // the partner may be any surface a face melts along: a plane against a
    // plane, but equally a wall band against the cylinder it copies. The
    // reach check keeps the answer honest: a point off the surface has no
    // chart position, whatever the inversion returns.
    use ogeom_math::elementary;
    let raw = match &face.surface {
        SurfaceGeometry::Plane(x) => {
            let local = x.plane().frame().to_local(p);
            if local.z.abs() > reach {
                return None;
            }
            Point2::new(local.x, local.y)
        }
        SurfaceGeometry::Cylinder(x) => {
            let cylinder = x.cylinder();
            if cylinder.distance_to(p) > reach {
                return None;
            }
            let (u, v) = elementary::cylinder_parameters(&cylinder, p, tol).ok()?;
            Point2::new(u, v)
        }
        SurfaceGeometry::Cone(x) => {
            let cone = x.cone();
            if cone.distance_to(p) > reach {
                return None;
            }
            let (u, v) = elementary::cone_parameters(&cone, p, tol).ok()?;
            Point2::new(u, v)
        }
        SurfaceGeometry::Sphere(x) => {
            let sphere = x.sphere();
            if sphere.distance_to(p) > reach {
                return None;
            }
            let (u, v) = elementary::sphere_parameters(&sphere, p, tol).ok()?;
            Point2::new(u, v)
        }
        SurfaceGeometry::Torus(x) => {
            let torus = x.torus();
            if torus.distance_to(p) > reach {
                return None;
            }
            let (u, v) = elementary::torus_parameters(&torus, p, tol).ok()?;
            Point2::new(u, v)
        }
        // No closed form (a fitted patch, a swept or revolved surface)
        // inverts by projection, held to the same reach: a blend's leg on
        // a spline host is a partner like any other.
        other => {
            let foot = ogeom_algo::project_on_surface(other, p, 24, tol).ok()?;
            if foot.distance > reach {
                return None;
            }
            Point2::new(foot.parameters.0, foot.parameters.1)
        }
    };
    let at = fold_point_into_chart(raw, &face.surface, tol);
    let lines = face.outline(tol)?;
    let trim = face.outline_trim(tol)?;
    // The face's boundary polylines are unwrapped (a winding ring may span
    // any one period's window, not necessarily the chart's canonical one,
    // and a wire chained onto one branch of the chart may sit whole periods
    // away from it), so the probe is carried to the outline's own branch
    // first, and asked at the neighbouring images as well.
    let mut shifts = vec![0.0];
    if face.surface.is_periodic_u() {
        let ((ua, ub), _) = face.surface.domain();
        if ub > ua {
            let period = ub - ua;
            let (sum, count) = lines
                .iter()
                .flatten()
                .fold((0.0_f64, 0_usize), |(s, n), q| (s + q.x, n + 1));
            #[allow(clippy::cast_precision_loss)]
            let centre = if count > 0 { sum / count as f64 } else { at.x };
            let home = ((centre - at.x) / period).round() * period;
            shifts = vec![home, home + period, home - period];
        }
    }
    for shift in shifts {
        let shifted = Point2::new(at.x + shift, at.y);
        if trim.inside(shifted) {
            return Some(shifted);
        }
    }
    None
}

/// Whether `p`, on the face's surface, stands inside the face, read at
/// the face's edge nearest it: where the nearest point of the boundary is
/// inside an edge, `p` stands off it square to it, and is inside just
/// where it stands on the side the face keeps its material on, the left
/// of the edge as the face walks it about its outward normal. `None`
/// where the nearest point is an edge's end (a corner reads both sides),
/// a seam's, or within the edge's own radius of it.
pub(crate) fn inside_beside_edge(
    face: &GFace,
    p: Point,
    tol: Tolerances,
) -> OgeomResult<Option<bool>> {
    let mut nearest: Option<(f64, usize, f64)> = None;
    for (k, e) in face.edges.iter().enumerate() {
        let (off, t) = nearest_on_edge_curve(&e.curve, e.crange, p, tol)?;
        if nearest.is_none_or(|n| off < n.0) {
            nearest = Some((off, k, t));
        }
    }
    let Some((off, k, t)) = nearest else {
        return Ok(None);
    };
    let e = &face.edges[k];
    let (lo, hi) = (e.crange.0.min(e.crange.1), e.crange.0.max(e.crange.1));
    let Some(forward) = e.walked else {
        return Ok(None);
    };
    if off <= e.tolerance.max(tol.confusion())
        || t <= lo + tol.parametric()
        || t >= hi - tol.parametric()
    {
        return Ok(None);
    }
    let at = e.pcurve.point_at(rescale(t, e.crange, e.prange), tol)?;
    let tangent = e.curve.d1_at(t, tol)?;
    let walked = if forward { tangent } else { -tangent };
    let left = outward_normal(face, at, tol)?.cross(walked);
    Ok(Some((p - e.curve.point_at(t, tol)?).dot(left) > 0.0))
}

/// Whether the point of `face`'s chart at `uv` lies on the face, read as a
/// partner's probe is read on it: a probe a piece's coarse outline puts
/// outside its face (beyond the chords' bow) is not the piece's.
pub(crate) fn on_own_face(face: &GFace, uv: Point2, tol: Tolerances) -> OgeomResult<bool> {
    let at = face.surface.point_at(uv.x, uv.y, tol)?;
    Ok(chart_point_within(face, at, partner_reach(face, face, tol), tol).is_some())
}

/// What the strands bounding a piece say of it, where its probes cannot.
pub(crate) enum Thin {
    /// The piece is no thinner than the band, or not planar: its probes
    /// stand.
    Wide,
    /// Thin, but some strand bounding it says nothing.
    Undecided,
    /// Thin, and every strand the other solid draws round it agrees.
    Reads(PieceState),
    /// Thin, and the strands disagree.
    Disagree,
}

/// A thin planar piece read by the other solid's edges it is bounded by.
///
/// A piece whose area is at most the classification band times half its
/// perimeter is on average narrower than the band, and no probe of it can
/// stand outside the band: each reads it on the other solid's boundary. A
/// strand that is the image of the other solid's edge `E`, between its
/// faces `g` on this face's surface and `h`, separates "on `g`" from
/// "beyond `g`". Which side the piece lies on is the arrangement's: its
/// rings keep it on their left about the chart's normal, as `g` keeps its
/// material on the left of its own rings about its outward normal. Beyond
/// `g` the other solid's wedge at `E` answers: past a convex edge (where
/// `h` turns away behind `g`) the plane beyond `g` is outside it, past a
/// concave one inside. Neither answer depends on how thin the piece is.
///
/// Only contact strands speak: a piece a section bounds is left to its
/// probes, as is one beside an edge whose two faces run on tangentially,
/// one it lies on `g`'s side of (a coincidence, which the partner reads),
/// or one whose faces are not both gathered.
pub(crate) fn thin_by_strands(
    face: &GFace,
    rings: &[Vec<Traversal<Tag>>],
    outlines: &[Vec<Point2>],
    contacts: &[ContactRec],
    theirs: &GSolid,
    tol: Tolerances,
) -> OgeomResult<Thin> {
    if !matches!(face.surface, SurfaceGeometry::Plane(_)) {
        return Ok(Thin::Wide);
    }
    // A plane's chart is its own frame, so lengths and areas in it are
    // the plane's.
    let (mut area, mut perimeter) = (0.0, 0.0);
    for (k, ring) in outlines.iter().enumerate() {
        let n = ring.len();
        let mut twice = 0.0;
        for i in 0..n {
            let (a, b) = (ring[i], ring[(i + 1) % n]);
            twice += a.x * b.y - b.x * a.y;
            perimeter += (b - a).magnitude();
        }
        area += if k == 0 { twice.abs() } else { -twice.abs() } / 2.0;
    }
    if area > tol.confusion() * 1e4 * perimeter / 2.0 {
        return Ok(Thin::Wide);
    }
    let mut read: Option<PieceState> = None;
    let mut spoke = false;
    for traversal in rings.iter().flatten() {
        let (contact, range) = match &traversal.tag {
            Tag::Contact { contact, range } => (&contacts[*contact], *range),
            Tag::Section { .. } => return Ok(Thin::Undecided),
            Tag::Boundary { .. } | Tag::Pole { .. } => continue,
        };
        let Some(state) = beyond_edge(face, contact, range, traversal.reversed, theirs, tol)?
        else {
            return Ok(Thin::Undecided);
        };
        spoke = true;
        match read {
            Some(said) if said != state => return Ok(Thin::Disagree),
            _ => read = Some(state),
        }
    }
    Ok(match read {
        Some(state) if spoke => Thin::Reads(state),
        _ => Thin::Undecided,
    })
}

/// Where the piece on the left of the contact strand over `range`, walked
/// back where `reversed`, stands against the solid whose edge the contact
/// is: `None` where the strand says nothing (see [`thin_by_strands`]).
pub(crate) fn beyond_edge(
    face: &GFace,
    contact: &ContactRec,
    range: (f64, f64),
    reversed: bool,
    theirs: &GSolid,
    tol: Tolerances,
) -> OgeomResult<Option<PieceState>> {
    let t = f64::midpoint(range.0, range.1);
    let p = contact.curve.point_at(t, tol)?;
    let along = contact.curve.d1_at(t, tol)?;
    let walked = if (range.1 < range.0) != reversed {
        -along
    } else {
        along
    };
    let uv = contact
        .pcurve
        .point_at(rescale(t, contact.crange, contact.prange), tol)?;
    let chart = face.surface.normal_at(uv.x, uv.y, tol)?.vector();
    let Ok(left) = chart.cross(walked).normalized(tol) else {
        return Ok(None);
    };
    // The edge's two faces, each with the way into its material from the
    // edge and its outward normal there.
    let mut sides = Vec::with_capacity(2);
    for other in &theirs.faces {
        for e in other.edges.iter().filter(|e| e.node == contact.node) {
            let Some(forward) = e.walked else {
                return Ok(None);
            };
            let (_, te) = nearest_on_edge_curve(&e.curve, e.crange, p, tol)?;
            let at = e.pcurve.point_at(rescale(te, e.crange, e.prange), tol)?;
            let normal = outward_normal(other, at, tol)?;
            let tangent = e.curve.d1_at(te, tol)?;
            let tangent = if forward { tangent } else { -tangent };
            let Ok(into) = normal.cross(tangent).normalized(tol) else {
                return Ok(None);
            };
            let Ok(normal) = normal.normalized(tol) else {
                return Ok(None);
            };
            sides.push((into, normal));
        }
    }
    let [(into_a, normal_a), (into_b, normal_b)] = sides[..] else {
        return Ok(None);
    };
    // `g` is the face lying along this one: its way in runs along the
    // piece's, one way or the other.
    let along_this = |into: ogeom_math::Vector| into.cross(left).magnitude() < TANGENT_TURN;
    let ((into_g, normal_g), into_h) = match (along_this(into_a), along_this(into_b)) {
        (true, false) => ((into_a, normal_a), into_b),
        (false, true) => ((into_b, normal_b), into_a),
        _ => return Ok(None),
    };
    if into_g.dot(left) > 0.0 {
        return Ok(None);
    }
    let turn = into_h.dot(normal_g);
    if turn.abs() < TANGENT_TURN {
        return Ok(None);
    }
    Ok(Some(if turn < 0.0 {
        PieceState::Out
    } else {
        PieceState::In
    }))
}

/// The sine of the turn under which two curves meet tangentially at a
/// crossing: the parallel-noise gate for pave placement. One degree is far
/// below any deliberate crossing and far above the smear a tangency leaves
/// in the general intersector.
pub(crate) const TANGENT_TURN: f64 = 2e-2;

/// Whether two curves meet tangentially at a crossing: their tangents turn
/// by less than [`TANGENT_TURN`].
pub(crate) fn tangential(
    a: &Curve,
    ta: f64,
    b: &Curve,
    tb: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let da = a.d1_at(ta, tol)?;
    let db = b.d1_at(tb, tol)?;
    let (ma, mb) = (da.magnitude(), db.magnitude());
    if ma <= tol.confusion() || mb <= tol.confusion() {
        return Ok(true);
    }
    Ok(da.cross(db).magnitude() / (ma * mb) < TANGENT_TURN)
}

/// The chart image of an exact curve on a surface, fitted at the curve's
/// own parameters from its points inverted on the surface, with the
/// distance the image honestly sits off the curve, or nothing, when the
/// fit misses its budget by more than a chart's swing beside a pole allows.
pub(crate) fn fitted_image(
    curve: &Curve,
    surface: &SurfaceGeometry,
    budget: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(PlanarCurve, f64)>> {
    use ogeom_geom::Surface as _;
    const SAMPLES: usize = 96;
    let (lo, hi) = curve.domain();
    if !hi.is_finite() || !lo.is_finite() || hi <= lo {
        return Ok(None);
    }
    let ((ua, ub), (va, vb)) = surface.domain();
    let periods = (
        if surface.is_periodic_u() {
            ub - ua
        } else {
            0.0
        },
        if surface.is_periodic_v() {
            vb - va
        } else {
            0.0
        },
    );
    // The image sampled where it moves: beside a pole a circle's image
    // swings half a turn within a hair of arc, and a uniform sampling puts
    // one point in the swing. Between any two samples whose images sit
    // further apart than a small step of the chart a sample is added, until
    // the image is walked at that step or the budget of samples is spent.
    let invert = |t: f64, guess: Option<(f64, f64)>| -> OgeomResult<Option<((f64, f64), f64)>> {
        let p = curve.point_at(t, tol)?;
        // A plane's image is its frame's own coordinates, exact and
        // unbounded: the face's chart window is a trim, not a limit on
        // where the curve may be read, and a projection clamped to it puts
        // the circle's far side off the surface by the window's shortfall.
        if let SurfaceGeometry::Plane(plane) = surface {
            let local = plane.plane().frame().to_local(p);
            let back = plane
                .plane()
                .frame()
                .to_world(ogeom_math::Point::new(local.x, local.y, 0.0));
            return Ok(Some(((local.x, local.y), back.distance(p))));
        }
        // A sphere's image is its frame's longitude and latitude, exact
        // where the projection's Newton cannot refine beside a pole and
        // hands back a coarse scan sample instead. The reading is checked
        // round the trip and left to the projection if the chart's
        // convention is not this one.
        if let SurfaceGeometry::Sphere(ball) = surface {
            let frame = ball.sphere().frame();
            let local = frame.to_local(p);
            let flat = local.x.hypot(local.y);
            // Kept a hair inside the poles: a spline fitted through points
            // on the pole itself overshoots it between them, and the chart
            // refuses a latitude past its end.
            let half = core::f64::consts::FRAC_PI_2 - 1e-6;
            let v = local.z.atan2(flat).clamp(-half, half);
            let u = local.y.atan2(local.x).rem_euclid(core::f64::consts::TAU);
            let back = surface.point_at(u, v, tol)?;
            if back.distance(p) <= budget.max(tol.confusion() * 1e2) {
                return Ok(Some(((u, v), back.distance(p))));
            }
        }
        // A cone's image is read in closed form too, past its window's ends:
        // an open conic it cuts is stated over a span that climbs beyond the
        // window, which the face trims, and a projection clamped to the
        // window puts those points off the surface.
        if let SurfaceGeometry::Cone(cone) = surface
            && let Ok((u, v)) = ogeom_math::elementary::cone_parameters(&cone.cone(), p, tol)
        {
            let back = ogeom_math::elementary::cone_at(&cone.cone(), u, v).point;
            if back.distance(p) <= budget.max(tol.confusion() * 1e2) {
                return Ok(Some(((u, v), back.distance(p))));
            }
        }
        let mut foot = match guess {
            Some(g) => ogeom_algo::project_on_surface_from(surface, p, g, tol)?,
            None => ogeom_algo::project_on_surface(surface, p, 24, tol)?,
        };
        if foot.distance > budget.max(tol.confusion() * 1e2) && guess.is_some() {
            foot = ogeom_algo::project_on_surface(surface, p, 48, tol)?;
        }
        if foot.distance > budget.max(tol.confusion() * 1e2) {
            ogeom_core::diag!(
                Channel::Wire,
                "IMAGE FIT: a point sits {:.2e} off the surface; no image",
                foot.distance
            );
            return Ok(None);
        }
        Ok(Some((foot.parameters, foot.distance)))
    };
    let unwrap = |uv: (f64, f64), before: Option<Point2>| -> Point2 {
        let mut uv = Point2::new(uv.0, uv.1);
        if let Some(prev) = before {
            for (coord, period, before) in [
                (&mut uv.x, periods.0, prev.x),
                (&mut uv.y, periods.1, prev.y),
            ] {
                if period > 0.0 {
                    while *coord - before > period / 2.0 {
                        *coord -= period;
                    }
                    while before - *coord > period / 2.0 {
                        *coord += period;
                    }
                }
            }
        }
        uv
    };
    let mut samples: Vec<(f64, Point2)> = Vec::with_capacity(SAMPLES + 1);
    let mut guess: Option<(f64, f64)> = None;
    for k in 0..=SAMPLES {
        #[allow(clippy::cast_precision_loss)]
        let t = lo + (hi - lo) * (k as f64) / (SAMPLES as f64);
        let Some((uv, _)) = invert(t, guess)? else {
            return Ok(None);
        };
        guess = Some(uv);
        let before = samples.last().map(|(_, p)| *p);
        samples.push((t, unwrap(uv, before)));
    }
    let step = {
        let (su, sv) = (
            if periods.0 > 0.0 { periods.0 } else { ub - ua },
            if periods.1 > 0.0 { periods.1 } else { vb - va },
        );
        su.min(sv) * 0.01
    };
    const CAP: usize = 4096;
    let mut refined = true;
    while refined && samples.len() < CAP {
        refined = false;
        let mut next: Vec<(f64, Point2)> = Vec::with_capacity(samples.len() * 2);
        for pair in samples.windows(2) {
            let ((t0, p0), (t1, p1)) = (pair[0], pair[1]);
            next.push((t0, p0));
            if p0.distance(p1) > step && t1 - t0 > (hi - lo) * 1e-7 && next.len() < CAP {
                let tm = f64::midpoint(t0, t1);
                let Some((uv, _)) = invert(tm, Some((p0.x, p0.y)))? else {
                    return Ok(None);
                };
                next.push((tm, unwrap(uv, Some(p0))));
                refined = true;
            }
        }
        next.push(samples[samples.len() - 1]);
        samples = next;
    }
    let params: Vec<f64> = samples.iter().map(|(t, _)| *t).collect();
    let image: Vec<Point2> = samples.iter().map(|(_, p)| *p).collect();
    let fitted =
        ogeom_geom::fit::fit_points_2d_at(&params, &image, 3, tol.confusion() * 10.0, tol)?;
    let planar: PlanarCurve = fitted.curve.into();
    // The honest error is in space: the surface read through the image
    // against the curve itself, between the samples as well as at them.
    // Checked at the samples and between each pair of them: where the
    // image swings the samples are dense, and a uniform check would step
    // over the swing and read the fit as honest.
    // An image that leaves the chart between its samples (a spline
    // overshooting the pole it was fitted up to) is no image either.
    let mut off = 0.0_f64;
    let mut checks: Vec<f64> = Vec::with_capacity(params.len() * 2);
    for pair in params.windows(2) {
        checks.push(pair[0]);
        checks.push(f64::midpoint(pair[0], pair[1]));
    }
    if let Some(&t) = params.last() {
        checks.push(t);
    }
    // A cone is read past its window's ends, as its points were inverted.
    // So is a cylinder where its window refuses the point: a circle a
    // rounding off square to the axis climbs a hair past the wall's end
    // round half its turn, and the fit through points held at that end
    // follows it there.
    let lift = |u: f64, v: f64| match surface {
        SurfaceGeometry::Cone(cone) => {
            Ok(ogeom_math::elementary::cone_at(&cone.cone(), u, v).point)
        }
        SurfaceGeometry::Cylinder(drum) => surface
            .point_at(u, v, tol)
            .or_else(|_| Ok(ogeom_math::elementary::cylinder_at(&drum.cylinder(), u, v).point)),
        _ => surface.point_at(u, v, tol),
    };
    for t in checks {
        let uv = planar.point_at(t, tol)?;
        let Ok(on_surface) = lift(uv.x, uv.y) else {
            ogeom_core::diag!(
                Channel::Wire,
                "IMAGE FIT: the image leaves the chart at {uv:?}; no image"
            );
            return Ok(None);
        };
        off = off.max(on_surface.distance(curve.point_at(t, tol)?));
    }
    ogeom_core::diag!(
        Channel::Wire,
        "IMAGE FIT: image off the curve by {off:.2e} against {:.2e} (fit error {:.2e} met {})",
        budget.max(tol.confusion() * 1e3),
        fitted.error,
        fitted.met
    );
    if off > budget.max(tol.confusion() * 1e3) {
        return Ok(None);
    }
    Ok(Some((planar, off)))
}

/// The distance from a point to a bounded edge curve, through a sampling
/// fine enough for the along-boundary question it answers.
pub(crate) fn distance_to_edge_curve(
    curve: &Curve,
    crange: (f64, f64),
    p: Point,
    tol: Tolerances,
) -> OgeomResult<f64> {
    Ok(nearest_on_edge_curve(curve, crange, p, tol)?.0)
}

/// Absolute distance for the hug's strict radius test. Other sampling
/// queries compare the growth of their sampled distances, or use their
/// feet to seed paving; they retain their own sampling parameterization.
pub(crate) fn hug_distance(
    curve: &Curve,
    crange: (f64, f64),
    p: Point,
    tol: Tolerances,
) -> OgeomResult<f64> {
    let (_, seed) = nearest_on_edge_curve(curve, crange, p, tol)?;
    // A chord's foot is only a seed: its distance can be above or below
    // the distance to the curve, even for a point lying exactly on it.
    // Refine on the curve itself before asking the micron-wide hug filters.
    let fine = (crange.1 - crange.0) / (48.0 * 8.0);
    let distance_at = |t| curve.point_at(t, tol).map(|at| at.square_distance(p));
    let refined = ogeom_math::solve::minimize(
        |t| distance_at(t).unwrap_or(f64::INFINITY),
        (seed - fine).max(crange.0),
        (seed + fine).min(crange.1),
        ogeom_math::solve::Criteria {
            residual: 0.0,
            step: tol.parametric(),
            max_iterations: 100,
        },
    )?;
    let mut best = distance_at(seed)?;
    for t in [crange.0, crange.1, refined.value] {
        best = best.min(distance_at(t)?);
    }
    Ok(best.sqrt())
}

/// [`distance_to_edge_curve`], with the parameter of the nearest point.
pub(crate) fn nearest_on_edge_curve(
    curve: &Curve,
    crange: (f64, f64),
    p: Point,
    tol: Tolerances,
) -> OgeomResult<(f64, f64)> {
    // Coarse bracket, then two rounds of local refinement: the answer feeds
    // the hug filters, whose widths are fractions of a millimetre, and a
    // long arc's 48-segment polyline sags by more than that on its own.
    let scan = |lo: f64, hi: f64, steps: u32| -> OgeomResult<(f64, f64)> {
        let mut best = f64::INFINITY;
        let mut best_t = lo;
        let mut previous: Option<(f64, Point)> = None;
        for i in 0..=steps {
            let t = lo + (hi - lo) * f64::from(i) / f64::from(steps);
            let at = curve.point_at(t, tol)?;
            if let Some((t0, last)) = previous {
                let d = at - last;
                let len2 = d.dot(d);
                let s = if len2 > 0.0 {
                    ((p - last).dot(d) / len2).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let dist = p.distance(last + d * s);
                if dist < best {
                    best = dist;
                    best_t = (t - t0).mul_add(s, t0);
                }
            }
            previous = Some((t, at));
        }
        Ok((best, best_t))
    };
    let span = crange.1 - crange.0;
    let (_, t1) = scan(crange.0, crange.1, 48)?;
    let step = span / 48.0;
    let (_, t2) = scan((t1 - step).max(crange.0), (t1 + step).min(crange.1), 16)?;
    let fine = span / (48.0 * 8.0);
    scan((t2 - fine).max(crange.0), (t2 + fine).min(crange.1), 16)
}

/// Pair up the pieces two coincident faces contribute for the same patch of
/// one surface, and mark the second argument's copy as stood in for where
/// the material is on the same side.
///
/// Returns each pair, first argument's piece first, the opposed ones
/// included: whichever of the two the operation drops, its face is recorded
/// as having become the one it keeps.
///
/// The substitution is by *region*, not by argument: a piece of `b` is a
/// duplicate exactly where a piece of `a`, on the same side and on a surface
/// they share, already covers the point `b`'s piece was classified at. Where
/// `a` has nothing there (a bore refilled by the cylinder that cut it, whose
/// caps fill holes the part no longer has faces for), nothing stands in, and
/// `b`'s piece is the only description of that patch there is.
pub(crate) fn mark_covered_coincidences(
    ga: &GSolid,
    gb: &GSolid,
    pieces: &mut [FacePiece],
    tol: Tolerances,
) -> Vec<(usize, usize)> {
    let on = |state: PieceState| matches!(state, PieceState::OnAligned | PieceState::OnOpposed);
    // Which pieces of the first argument stand on a shared surface.
    let from_a: Vec<(usize, usize, PieceState)> = pieces
        .iter()
        .enumerate()
        .filter(|(_, p)| p.from_a && on(p.state))
        .map(|(i, p)| (i, p.face, p.state))
        .collect();
    if from_a.is_empty() {
        return Vec::new();
    }
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for (index, piece) in pieces.iter().enumerate() {
        if piece.from_a || !on(piece.state) {
            continue;
        }
        for &(other, face_a, state) in &from_a {
            if state != piece.state {
                continue;
            }
            let host = &ga.faces[face_a];
            // The point `b`'s piece stands at, read in `a`'s face's chart,
            // as far off it as the two faces say they may stand. A point
            // that is not on that surface at all has no chart position.
            let reach = partner_reach(&gb.faces[piece.face], host, tol);
            let Some(at) = chart_point_within(host, piece.probe, reach, tol) else {
                continue;
            };
            // The host piece's rings, not its boundary strands: a ring is a
            // closed loop whose last point does not repeat its first, so the
            // test has to close it. Asked as though the rings were strands
            // that jointly close, the segment from the ring's end back to its
            // start goes uncounted, and a point the ring plainly encloses
            // comes back outside whenever that missing segment would have
            // been crossed.
            if inside_rings(&pieces[other].outlines, at) {
                pairs.push((other, index));
                break;
            }
        }
    }
    for &(_, index) in &pairs {
        if pieces[index].state == PieceState::OnAligned {
            pieces[index].covered = true;
        }
    }
    pairs
}

thread_local! {
    /// Whether a piece read on the other solid's boundary, with no coincident
    /// partner, is settled by asking just off it on both sides. Off on the
    /// first attempt, so a configuration the nested fallback answers exactly
    /// keeps that answer. On for the second. Read once per boolean on the
    /// calling thread.
    pub(crate) static SETTLE_FROM_BOTH_SIDES: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// Whether a piece that reads in or out is read again at the centre of
    /// its area. Off but for the attempt after kept pieces fail to close.
    /// Read once per boolean on the calling thread.
    pub(crate) static PROBE_AT_CENTRE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// The farthest apart a contact and a boundary edge it was read along
    /// stand, where that is past what either edge states: the paving takes
    /// the two for one edge and splits neither face along them, so only the
    /// sew can join them. Zero where every such pair lies within its stated
    /// tolerances. Gathered by the attempts since [`settled`] began.
    pub(crate) static ALONG_APART: std::cell::Cell<f64> = const { std::cell::Cell::new(0.0) };
    /// The gap the kept pieces are sewn across: zero but for the attempt
    /// [`settled`] makes last, across the distance [`ALONG_APART`] measured.
    /// Read once per assembly on the calling thread.
    pub(crate) static SEW_ACROSS: std::cell::Cell<f64> = const { std::cell::Cell::new(0.0) };
    /// Whether contacts are read at the radii their curves hold to (see
    /// [`settled`]). Off but for the last attempt. Read once per boolean
    /// on the calling thread.
    pub(crate) static AT_STATED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// Which edges a section line also ends where it crosses, for the edges
    /// it carries onto the other face (see [`graze_contacts`]). Off but for
    /// the retries run again once every other has refused a piece read on
    /// the other solid's boundary. Read once per boolean on the calling
    /// thread.
    pub(crate) static FROM_CROSSINGS: std::cell::Cell<Crossings> = const { std::cell::Cell::new(Crossings::Off) };
    /// The farthest apart two faces may stand and still be dropped as a
    /// membrane (see [`without_membranes`]), whatever weld the attempt
    /// itself runs at: during the retries at a widened confusion, the weld
    /// their result is held to, and unbounded otherwise. A shell dropped
    /// whole leaves no edge or vertex behind to state how far apart its
    /// faces stood. Read once per assembly on the calling thread.
    pub(crate) static MEMBRANE_WITHIN: std::cell::Cell<f64> = const { std::cell::Cell::new(f64::INFINITY) };
    /// Whether a plane's section with an extrusion of a B-spline profile
    /// crossing a rim of it is split at the rim in closed form
    /// ([`ogeom_intersect::plane_extrusion_at_rims`]) rather than marched.
    /// Off but for the retries [`settled`] runs once every other has
    /// refused. Read once per assembly on the calling thread.
    pub(crate) static AT_RIMS: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// Whether an assembly since [`settled`] began met such a pair with
    /// [`AT_RIMS`] off.
    pub(crate) static RIM_CROSSED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// Whether two faces are asked whether they coincide within the
    /// tolerances they state ([`coincide_as_stated`]) where those are the
    /// default ones too. Off but for the attempt [`settled`] makes last.
    /// Read once per assembly on the calling thread.
    pub(crate) static AS_STATED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// An operation run once, and where it refuses, again in three steps. Where
/// the refusal is a piece on the other solid's boundary
/// with no coincident partner, once more with such pieces settled from both
/// sides of them. Where the kept pieces did not close, once more with each
/// piece read at the centre of its area: a face
/// leaning off the other's by microradians, the two meeting along an edge
/// and nothing splitting it, stands within rounding of the other near that
/// edge, and a probe there can read the wrong side. And where it still
/// refuses, at a confusion distance ten
/// and then a hundred times wider: a sliver a few microns thick (two walls
/// a few microns apart, a face tilted a few microradians) is then the kind
/// the boolean welds, and the result is kept only if nothing in it states
/// a tolerance past the weld distance, a hundred confusion distances (and
/// the ten the welding run states its own vertices to). A shell those
/// retries drop as a membrane leaves nothing behind to state that, so it
/// is dropped only where its faces stand within that same weld of each
/// other, and a sliver any thicker is never lost to the retry.
///
/// Then, at the given confusion, with contacts read at the radii their
/// curves hold to rather than at the loose width a fitted curve is
/// measured within: a fitted curve leaving an exact one tangentially
/// leaves a lune between them wider than those radii, which is geometry.
/// Each stretch the two share ends where they part by those radii,
/// measured once for the pair so both sides split at one vertex; a shallow
/// crossing splits both curves, at a vertex of either it cannot be told
/// from; a contact stretch near an edge is read inside or out at that
/// edge; a fitted section's hug width welds only a face it leaves a strand
/// on; and a piece thinner than the classification band is read by the
/// other solid's edges bounding it ([`thin_by_strands`]).
///
/// Last, where the refusal was a piece on the other solid's boundary and
/// every retry refused too, the same retries again with the line edges a
/// section line crosses carried onto the other face from the crossing on
/// ([`graze_contacts`]): a strip of one wall between such a crossing and
/// its edge, standing on another wall a hair's angle off it, is then split
/// off both and read as one contact. Where those refuse too, once more
/// with edges of any curve carried so, the result kept only where nothing
/// in it states a tolerance past the weld distance: an arc's strip is
/// shorter than a line's, and a junction welding its two ends together
/// would weld the faces beside it into one. What the retries without them
/// resolve is resolved as before.
///
/// Where all of that refuses and a plane crossed an extrusion of a
/// B-spline profile across a rim of it, whose marched section stops short
/// of the rim, the run once more from the first attempt with that section
/// split where it crosses the rim in closed form ([`AT_RIMS`]), its result
/// kept only where nothing in it states a tolerance past the weld
/// distance.
///
/// And where the kept pieces still did not close and the first attempt
/// read a contact along a boundary edge it stands further from than either
/// states (a floor edge lying its tolerance off a mirror plane, against its
/// mirror image as far off the other way), once more with the pieces sewn
/// across that distance ([`SEW_ACROSS`]), where it is within the weld
/// distance the retries above hold to: the paving took the two for one
/// edge and split neither face along them, so the edge kept widens to
/// reach where the other ran. The attempt runs at the given confusion, and
/// the sew widens nothing further than that distance beyond what the
/// pieces already state, so its result is not held to the weld distance a
/// second time: a part may itself state wider tolerances.
///
/// After all of that, where the refusal was a piece on the other solid's
/// boundary with no coincident partner and every retry above refused too,
/// once more with every pair of faces asked whether they coincide within
/// the tolerances they state, the default ones included ([`AS_STATED`]). A
/// face leaning a rounding off another it meets along an edge (a pocket's
/// floor sketched a tenth of a micron off the plane it stands on) parts
/// from it by a hair past the confusion distance at its far edge: the
/// closed form reads the two as crossing along the edge they share, which
/// splits neither, and the leaning face reads on the other's boundary
/// everywhere with no partner. Within what the two faces state they are one
/// surface, and each face is split by the other's boundary as a coincident
/// face is.
pub(crate) fn settled(
    model: &mut Model,
    tol: Tolerances,
    mut attempt: impl FnMut(&mut Model, Tolerances) -> OgeomResult<Built>,
) -> OgeomResult<Built> {
    RIM_CROSSED.set(false);
    ALONG_APART.set(0.0);
    let refusal = match attempt(model, tol) {
        Ok(built) => return Ok(built),
        Err(e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_))) => {
            return Err(e);
        }
        Err(e) => e,
    };
    let apart = ALONG_APART.get();
    if let Some(built) = resettled(model, tol, &refusal, &mut attempt)? {
        return Ok(built);
    }
    if apart > 0.0
        && apart + tol.confusion() <= weld_limit(tol)
        && refusal.to_string().contains("did not close into a shell")
    {
        struct Plain;
        impl Drop for Plain {
            fn drop(&mut self) {
                SEW_ACROSS.set(0.0);
            }
        }
        SEW_ACROSS.set(apart + tol.confusion());
        let _plain = Plain;
        match attempt(model, tol) {
            Ok(built) => return Ok(built),
            Err(e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_))) => {
                return Err(e);
            }
            Err(_) => {}
        }
    }
    if refusal.to_string().contains("no coincident partner face") {
        struct Plain;
        impl Drop for Plain {
            fn drop(&mut self) {
                AS_STATED.set(false);
            }
        }
        AS_STATED.set(true);
        let _plain = Plain;
        match attempt(model, tol) {
            Ok(built) => return Ok(built),
            Err(e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_))) => {
                return Err(e);
            }
            Err(_) => {}
        }
    }
    Err(refusal)
}

/// The retries [`settled`] runs after the first attempt refused with
/// `refusal`, short of the last, and the first result one of them keeps.
pub(crate) fn resettled(
    model: &mut Model,
    tol: Tolerances,
    refusal: &ogeom_core::OgeomError,
    attempt: &mut impl FnMut(&mut Model, Tolerances) -> OgeomResult<Built>,
) -> OgeomResult<Option<Built>> {
    if let Some(built) = settled_after(model, tol, refusal, attempt)? {
        return Ok(Some(built));
    }
    if !RIM_CROSSED.get() {
        return Ok(None);
    }
    struct Rims;
    impl Drop for Rims {
        fn drop(&mut self) {
            AT_RIMS.set(false);
        }
    }
    AT_RIMS.set(true);
    let _rims = Rims;
    let limit = weld_limit(tol);
    let at_rims = match attempt(model, tol) {
        Ok(built) if welded_within(model, &built.shape, limit)? => return Ok(Some(built)),
        Err(e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_))) => {
            return Err(e);
        }
        Ok(_) => return Ok(None),
        Err(e) => e,
    };
    if let Some(built) = settled_after(model, tol, &at_rims, attempt)?
        && welded_within(model, &built.shape, limit)?
    {
        return Ok(Some(built));
    }
    Ok(None)
}

/// The retries [`settled`] runs after an attempt refused with `refusal`,
/// pieces on the other solid's boundary settled from both sides of them,
/// and the first result one of them keeps.
pub(crate) fn settled_after(
    model: &mut Model,
    tol: Tolerances,
    refusal: &ogeom_core::OgeomError,
    attempt: &mut impl FnMut(&mut Model, Tolerances) -> OgeomResult<Built>,
) -> OgeomResult<Option<Built>> {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            SETTLE_FROM_BOTH_SIDES.set(false);
        }
    }
    SETTLE_FROM_BOTH_SIDES.set(true);
    let _reset = Reset;
    if let Some(built) = retried(model, tol, refusal, attempt)? {
        return Ok(Some(built));
    }
    if refusal.to_string().contains("no coincident partner face") {
        struct Plain;
        impl Drop for Plain {
            fn drop(&mut self) {
                FROM_CROSSINGS.set(Crossings::Off);
            }
        }
        let _plain = Plain;
        FROM_CROSSINGS.set(Crossings::Lines);
        if let Some(built) = retried(model, tol, refusal, attempt)? {
            return Ok(Some(built));
        }
        FROM_CROSSINGS.set(Crossings::Curves);
        if let Some(built) = retried(model, tol, refusal, attempt)?
            && welded_within(model, &built.shape, weld_limit(tol))?
        {
            return Ok(Some(built));
        }
    }
    Ok(None)
}

/// The retries [`settled`] runs after a first attempt refused with
/// `refusal`, each in turn, and the first result one of them keeps.
pub(crate) fn retried(
    model: &mut Model,
    tol: Tolerances,
    refusal: &ogeom_core::OgeomError,
    attempt: &mut impl FnMut(&mut Model, Tolerances) -> OgeomResult<Built>,
) -> OgeomResult<Option<Built>> {
    if refusal.to_string().contains("no coincident partner face") {
        match attempt(model, tol) {
            Ok(built) => return Ok(Some(built)),
            Err(ogeom_core::OgeomError::Cancelled) => {
                return Err(ogeom_core::OgeomError::Cancelled);
            }
            Err(_) => {}
        }
    }
    if refusal.to_string().contains("did not close into a shell") {
        struct Clear;
        impl Drop for Clear {
            fn drop(&mut self) {
                PROBE_AT_CENTRE.set(false);
            }
        }
        PROBE_AT_CENTRE.set(true);
        let _clear = Clear;
        match attempt(model, tol) {
            Ok(built) => return Ok(Some(built)),
            Err(ogeom_core::OgeomError::Cancelled) => {
                return Err(ogeom_core::OgeomError::Cancelled);
            }
            Err(_) => {}
        }
    }
    let limit = weld_limit(tol);
    struct Restore(f64);
    impl Drop for Restore {
        fn drop(&mut self) {
            MEMBRANE_WITHIN.set(self.0);
        }
    }
    let restore = Restore(MEMBRANE_WITHIN.replace(limit.min(MEMBRANE_WITHIN.get())));
    for widen in [10.0, 100.0] {
        let Ok(wider) = Tolerances::with_scale(tol.scale() / widen) else {
            continue;
        };
        match attempt(model, wider) {
            Ok(built) if welded_within(model, &built.shape, limit)? => return Ok(Some(built)),
            Err(ogeom_core::OgeomError::Cancelled) => {
                return Err(ogeom_core::OgeomError::Cancelled);
            }
            _ => {}
        }
    }
    drop(restore);
    {
        struct Off;
        impl Drop for Off {
            fn drop(&mut self) {
                AT_STATED.set(false);
            }
        }
        AT_STATED.set(true);
        let _off = Off;
        #[cfg(feature = "test-hooks")]
        work::add(|w| w.stated_retries += 1);
        match attempt(model, tol) {
            Ok(built) => return Ok(Some(built)),
            Err(ogeom_core::OgeomError::Cancelled) => {
                return Err(ogeom_core::OgeomError::Cancelled);
            }
            Err(_) => {}
        }
    }
    Ok(None)
}

/// A pole's pcurve over `range` as the straight line it runs along in the
/// chart, where it does (a line trimmed to the pole's stretch, as a
/// boolean writes a pole edge back): the poles are cut where strands end
/// on them by reading their line, and a pole read as some other curve was
/// never cut. Any other curve is kept as it is.
pub(crate) fn straight_pole(
    planar: &PlanarCurve,
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<(PlanarCurve, (f64, f64))> {
    use ogeom_geom::Curve2d as _;
    if matches!(planar, PlanarCurve::Line(_)) {
        return Ok((planar.clone(), range));
    }
    let (a, b) = (
        planar.point_at(range.0, tol)?,
        planar.point_at(range.1, tol)?,
    );
    let length = a.distance(b);
    if length <= tol.parametric() {
        return Ok((planar.clone(), range));
    }
    for k in 1..8 {
        let t = range.0 + (range.1 - range.0) * f64::from(k) / 8.0;
        let (run, off) = (b - a, planar.point_at(t, tol)? - a);
        if (run.x * off.y - run.y * off.x).abs() / length > tol.parametric() {
            return Ok((planar.clone(), range));
        }
    }
    let line = ogeom_geom::Line2d::segment(a, b, tol)?;
    let domain = line.domain();
    Ok((line.into(), domain))
}

/// The widest tolerance a retried result may state: the weld distance, and
/// the wider confusion the welding run itself states its vertices to.
pub(crate) fn weld_limit(tol: Tolerances) -> f64 {
    (tol.confusion() * 1e2 + tol.confusion() * 10.0) * (1.0 + 1e-9)
}

/// Whether every edge and vertex of `shape` states a tolerance within
/// `weld`.
pub(crate) fn welded_within(model: &Model, shape: &Shape, weld: f64) -> OgeomResult<bool> {
    for kind in [ShapeType::Edge, ShapeType::Vertex] {
        for sub in explore_unique(model, shape, kind)? {
            if model.tolerance_of(&sub)?.is_some_and(|t| t.get() > weld) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// The strict audit's verdict: the filtered fill and the unfiltered fill
/// keep the same material.
///
/// Pieces are compared as merged parameter intervals per section: the
/// unfiltered run sees more curves, and a dropped pair's dead section can
/// still cross a live one and split its pieces differently, so individual
/// piece boundaries are noise and the kept *union* is the signal. Sections
/// are keyed by their face pair and a sampled point, because indices differ
/// between the runs.
pub(crate) fn audit_fill_equivalence(
    filtered: (&[SectionRec], &[SectionPiece]),
    unfiltered: (&[SectionRec], &[SectionPiece]),
    tol: Tolerances,
) {
    type Key = (usize, usize, [i64; 3]);
    let key_of = |section: &SectionRec, tol: Tolerances| -> Key {
        let (lo, hi) = section.curve.domain();
        let p = section
            .curve
            .point_at(f64::midpoint(lo, hi), tol)
            .unwrap_or(ogeom_math::Point::ORIGIN);
        let grid = (tol.confusion() * 100.0).max(1e-9);
        #[allow(clippy::cast_possible_truncation)]
        let q = |x: f64| (x / grid).round() as i64;
        (section.face_a, section.face_b, [q(p.x), q(p.y), q(p.z)])
    };
    let merged = |sections: &[SectionRec], pieces: &[SectionPiece], tol: Tolerances| {
        let mut kept: ogeom_core::FastMap<Key, Vec<(f64, f64)>> = ogeom_core::FastMap::default();
        for piece in pieces {
            kept.entry(key_of(&sections[piece.section], tol))
                .or_default()
                .push(piece.range);
        }
        let slop = tol.parametric() * 10.0;
        for ranges in kept.values_mut() {
            ranges.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(core::cmp::Ordering::Equal));
            let mut out: Vec<(f64, f64)> = Vec::with_capacity(ranges.len());
            for &(lo, hi) in ranges.iter() {
                match out.last_mut() {
                    Some(last) if lo <= last.1 + slop => last.1 = last.1.max(hi),
                    _ => out.push((lo, hi)),
                }
            }
            *ranges = out;
        }
        kept
    };
    let a = merged(filtered.0, filtered.1, tol);
    let b = merged(unfiltered.0, unfiltered.1, tol);
    let slop = tol.parametric() * 20.0;
    let matches = |x: &Vec<(f64, f64)>, y: &Vec<(f64, f64)>| {
        x.len() == y.len()
            && x.iter()
                .zip(y)
                .all(|(p, q)| (p.0 - q.0).abs() <= slop && (p.1 - q.1).abs() <= slop)
    };
    for (key, ranges) in &b {
        let held = a.get(key);
        assert!(
            held.is_some_and(|r| matches(r, ranges)),
            "strict bound-filter audit: the unfiltered fill keeps {ranges:?} \
             on section {key:?}, the filtered fill keeps {held:?}; the \
             filter drops material"
        );
    }
    for (key, ranges) in &a {
        assert!(
            b.contains_key(key),
            "strict bound-filter audit: the filtered fill keeps {ranges:?} \
             on section {key:?} the unfiltered fill never made"
        );
    }
}
