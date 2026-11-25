//! Rebuilding: the kept pieces made into faces, edges and shells.

use super::*;

// --- rebuilding --------------------------------------------------------------

/// Everything the rebuild shares across pieces.
pub(crate) struct Rebuild<'m> {
    model: &'m mut Model,
    /// World surface ids, minted once per source face.
    surfaces_a: Vec<Option<ogeom_topo::SurfaceId>>,
    surfaces_b: Vec<Option<ogeom_topo::SurfaceId>>,
    /// Vertices shared by position: a wire's connectivity is checked by node
    /// identity, so two sub-edges meeting at a point must *name* the same
    /// vertex, not merely coincide there.
    vertices: Vec<(Point, Shape)>,
    /// `vertices` binned by position.
    vertex_bins: bins::Bins,
    /// The entries a bin query found, kept between queries.
    near: Vec<usize>,
    /// The widest tolerance any vertex in `vertices` holds: every widening
    /// of one goes through [`Rebuild::vertex`] or [`Rebuild::widen`], so an
    /// end is compared only with the vertices this far and a weld from it.
    widest: f64,
    /// How far two honest descriptions of one junction may sit apart: a
    /// hundred confusions as the floor, widened to three times the loosest
    /// contact edge's or fitted section's own tolerance when one took part:
    /// the strand's own, the crossing's, and the edge's it stopped on.
    weld: f64,
    /// The paving's junctions, each minted as a vertex the first time a
    /// strand end lands inside it.
    junctions: Vec<(Junction, Option<Shape>)>,
    /// `junctions` binned by position.
    junction_bins: bins::Bins,
    /// The widest reach of any junction.
    junction_reach: f64,
    /// Vertices minted from junctions at an edge's own end vertex, which a
    /// tangent section's chart image is bent onto.
    onto_vertex: ogeom_core::FastSet<ogeom_topo::TShapeId>,
    /// Boundary edges rebuilt whole, by the edge they were: the face on the
    /// other side of one finds it built and adds its own pcurve, where a
    /// second build would leave the sew a twin to find.
    whole_edges: ogeom_core::FastMap<EdgeKey, Shape>,
    /// The vertices of edges shared with faces set aside, which the rebuild
    /// names as they are and never changes, each where it stands.
    kept: ogeom_core::FastMap<ogeom_topo::TShapeId, Point>,
    /// Whether a strand end asked a kept vertex to reach further than it
    /// does.
    strained: bool,
}

impl Rebuild<'_> {
    fn vertex(&mut self, p: Point, tol: Tolerances) -> Shape {
        // A junction several paves described is one vertex owning the span
        // they disagree by, and every strand end inside it names that
        // vertex, whichever pave its own trim stopped at.
        // The rim is as generous as the positional weld below: two strand
        // ends a weld apart must not fall on opposite sides of it.
        let slack = self.weld.max(tol.confusion() * 1e2);
        ogeom_core::diag!(Channel::Wire, "VERTEX ask {p:?}");
        let reach = self.junction_reach;
        let inside = |slot: &usize| {
            let j = &self.junctions[*slot].0;
            j.at.distance(p) <= j.reach + slack
        };
        let mut near = core::mem::take(&mut self.near);
        let first = if self.junction_bins.near_into(p, reach + slack, &mut near) {
            near.iter().copied().find(inside)
        } else {
            (0..self.junctions.len()).find(inside)
        };
        if let Some(slot) = first {
            let junction = self.junctions[slot].0;
            let shape = match &self.junctions[slot].1 {
                Some(shape) => shape.clone(),
                None => {
                    let shape = make_vertex(self.model, junction.at).shape;
                    self.junctions[slot].1 = Some(shape.clone());
                    self.remember(junction.at, &shape);
                    shape
                }
            };
            self.near = near;
            if junction.onto_vertex {
                self.onto_vertex.insert(shape.node());
            }
            // The junction owns its span, and an end welded in from the rim
            // widens it by what it actually sat off by, as any weld does.
            let gap = junction.at.distance(p);
            if let Some(node) = self.model.node_mut(&shape)
                && let ogeom_topo::NodeData::Vertex(data) = node.data_mut()
            {
                // A confusion of slack past the measured gap: an edge built
                // against this vertex measures the same gap through its own
                // rounding, and a tolerance equal to it fails by an ulp.
                data.tolerance = data
                    .tolerance
                    .widen_to(junction.reach.max(gap) + tol.confusion());
            }
            self.note_tolerance(&shape);
            return shape;
        }
        // The weld reach covers what the inputs may honestly disagree by: a
        // boundary curve that is itself a fitted intersection from an
        // earlier boolean carries a couple of microns of slop, and two
        // descriptions of one junction arrive that far apart. A hundred
        // confusions at millimetre tolerances is ten microns, below any
        // feature this pipeline can resolve, and every weld wider than
        // confusion is recorded on the vertex, not papered over.
        // And as far again as the vertex already owns: a vertex widened by
        // the descriptions it has taken in claims its point within that
        // tolerance, an end arriving within that and the weld of it may be
        // the same point, and a boolean cannot show it is not. Several
        // descriptions of one junction (three bands meeting at an apex)
        // can arrive in turn, each within the last's reach and none within
        // the first's, and measured from the first alone the wire round
        // the band's end has two vertices where it needs one. So every end
        // taken in is remembered where it arrived, and the next end is
        // measured from the description nearest it.
        let floor = self.weld.max(tol.confusion() * 1e2);
        if !self
            .vertex_bins
            .near_into(p, floor + self.widest, &mut near)
        {
            near.extend(0..self.vertices.len());
        }
        let found = near
            .iter()
            .filter_map(|&index| {
                let (q, shape) = &self.vertices[index];
                let own = self
                    .model
                    .node(shape)
                    .and_then(|n| n.data().as_vertex())
                    .map_or(0.0, |d| d.tolerance.get());
                let gap = q.distance(p);
                (gap <= floor.max(own + floor)).then_some((gap, index))
            })
            .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(core::cmp::Ordering::Equal))
            .map(|(gap, index)| (gap, self.vertices[index].1.clone()));
        self.near = near;
        if let Some((gap, shape)) = found {
            // Two descriptions of one junction may disagree by a general
            // crossing's residual. The vertex's tolerance is where that
            // disagreement is recorded, so the sub-edges built against
            // either description still reach it honestly.
            let (at, own) = self
                .model
                .node(&shape)
                .and_then(|n| n.data().as_vertex())
                .map_or((p, 0.0), |d| (d.point, d.tolerance.get()));
            let off = at.distance(p);
            if let Some(stands) = self.kept.get(&shape.node()) {
                // A kept vertex is never widened: an end it does not reach
                // as it stands strains the rebuild.
                let off = stands.distance(p);
                self.strained |= off > tol.confusion() && off + tol.confusion() > own;
            } else if off > tol.confusion()
                && let Some(node) = self.model.node_mut(&shape)
                && let ogeom_topo::NodeData::Vertex(data) = node.data_mut()
            {
                data.tolerance = data.tolerance.widen_to(off + tol.confusion());
            }
            self.note_tolerance(&shape);
            if gap > tol.confusion() {
                self.remember(p, &shape);
            }
            return shape;
        }
        let shape = make_vertex(self.model, p).shape;
        self.remember(p, &shape);
        shape
    }

    /// Take in `shape` as described at `p`.
    fn remember(&mut self, p: Point, shape: &Shape) {
        self.vertex_bins.insert(p, self.vertices.len());
        self.vertices.push((p, shape.clone()));
        self.note_tolerance(shape);
    }

    /// Keep `widest` covering `shape`'s tolerance.
    fn note_tolerance(&mut self, shape: &Shape) {
        if let Some(data) = self.model.node(shape).and_then(|n| n.data().as_vertex()) {
            self.widest = self.widest.max(data.tolerance.get());
        }
    }

    /// Widen a vertex this rebuild handed out. A kept vertex is not
    /// widened: one that would need it strains the rebuild.
    fn widen(&mut self, vertex: &Shape, to: f64, tol: Tolerances) -> OgeomResult<()> {
        if self.kept.contains_key(&vertex.node()) {
            let own = self
                .model
                .node(vertex)
                .and_then(|n| n.data().as_vertex())
                .map_or(0.0, |d| d.tolerance.get());
            self.strained |= own < to;
            return Ok(());
        }
        self.model
            .widen(vertex, ogeom_core::Tolerance::new_in(to, tol)?)?;
        self.note_tolerance(vertex);
        Ok(())
    }

    fn surface_id(
        &mut self,
        fused: &GeneralFused,
        from_a: bool,
        face: usize,
    ) -> ogeom_topo::SurfaceId {
        let slot = if from_a {
            &mut self.surfaces_a[face]
        } else {
            &mut self.surfaces_b[face]
        };
        if let Some(id) = slot {
            return *id;
        }
        let own = if from_a {
            fused.a.faces[face].own_surface
        } else {
            fused.b.faces[face].own_surface
        };
        if let Some(id) = own {
            *slot = Some(id);
            return id;
        }
        let surface = if from_a {
            fused.a.faces[face].surface.clone()
        } else {
            fused.b.faces[face].surface.clone()
        };
        let id = self.model.geometry_mut().add_surface(surface);
        *slot = Some(id);
        id
    }
}

/// The sub-edges built of each strand of a piece, keyed by strand, each
/// with the parameter range it covers, in the order they were built.
pub(crate) type StrandEdges = ogeom_core::FastMap<(usize, u8), Vec<((f64, f64), Shape)>>;

/// Build one piece as a face, orientation matching its source face's side.
pub(crate) fn build_piece(
    rebuild: &mut Rebuild,
    fused: &GeneralFused,
    piece: &FacePiece,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let own = if piece.from_a { &fused.a } else { &fused.b };
    let face = &own.faces[piece.face];
    let surface_id = rebuild.surface_id(fused, piece.from_a, piece.face);

    // Sub-edges cached within the piece, so a seam used from both sides is
    // one edge appearing twice. Keyed by strand, each with the ranges built
    // of it in order.
    let mut cache: StrandEdges = ogeom_core::FastMap::default();
    let mut wires = Vec::new();
    for ring in &piece.rings {
        let mut edges = Vec::with_capacity(ring.len());
        // The edges shared with faces set aside this ring has taken, each
        // the way it runs.
        let mut taken: Vec<(EdgeKey, bool)> = Vec::new();
        for traversal in ring {
            // An edge shared with a face set aside is itself, whole: nothing
            // reaches it to split it, and it carries its pcurve on the
            // face's own surface already. Its strands (a closed edge's two
            // halves) stand for it once each way it is walked.
            if let Tag::Boundary { edge, .. } = &traversal.tag
                && let Some(kept) = own.shared.get(&face.edges[*edge].node)
            {
                if face.own_surface != Some(surface_id) {
                    rebuild.strained = true;
                }
                let key = (face.edges[*edge].node, traversal.reversed);
                if !taken.contains(&key) {
                    taken.push(key);
                    edges.push(if traversal.reversed {
                        kept.reversed()
                    } else {
                        kept.clone()
                    });
                }
                continue;
            }
            let (key_edge, key_kind, range) = match &traversal.tag {
                Tag::Boundary { edge, range } => (*edge, 0_u8, *range),
                Tag::Section { section, range } => (*section, 1, *range),
                Tag::Contact { contact, range } => (*contact, 2, *range),
                Tag::Pole { pole, range } => (*pole, 3, *range),
            };
            let near = |a: (f64, f64), b: (f64, f64)| {
                (a.0 - b.0).abs() <= tol.parametric() && (a.1 - b.1).abs() <= tol.parametric()
            };
            if diag::enabled(Channel::Wire) && (range.1 - range.0).abs() < 1e-4 {
                let kind = match &traversal.tag {
                    Tag::Boundary { .. } => "boundary",
                    Tag::Section { .. } => "section",
                    Tag::Contact { .. } => "contact",
                    Tag::Pole { .. } => "pole",
                };
                ogeom_core::diag!(
                    Channel::Wire,
                    "DUST {kind} {key_edge} range {range:?} in piece from_a={} face={}",
                    piece.from_a,
                    piece.face
                );
            }
            let built = if let Some((_, shape)) = cache
                .get(&(key_edge, key_kind))
                .and_then(|built| built.iter().find(|(r, _)| near(*r, range)))
            {
                shape.clone()
            } else {
                let shape = build_sub_edge(
                    rebuild,
                    fused,
                    piece.from_a,
                    face,
                    surface_id,
                    (&traversal.tag, traversal.nodes),
                    tol,
                )?;
                cache
                    .entry((key_edge, key_kind))
                    .or_default()
                    .push((range, shape.clone()));
                shape
            };
            // A piece whose two ends welded to one vertex and whose whole
            // length lies within the vertex's reach is that vertex's dust
            // (a rim's last fraction of a micron before a pole corner) and
            // no edge of the wire. A closed edge on one vertex, a full
            // circle, reaches far from it and stays.
            if !matches!(traversal.tag, Tag::Pole { .. })
                && let Some((v0, v1)) = ogeom_algo::edge_vertices(rebuild.model, &built)?
                && v0.node() == v1.node()
                && let Some(vd) = rebuild
                    .model
                    .node(&v0)
                    .and_then(|n| n.data().as_vertex())
                    .map(|d| (d.point, d.tolerance.get()))
                && let Some(ed) = rebuild.model.node(&built).and_then(|n| n.data().as_edge())
                && let Some(ogeom_topo::EdgeRepr::Curve3d {
                    curve, range: r3, ..
                }) = ed.curve3d()
                && let Some(g) = rebuild.model.geometry().curve(*curve)
                && g.point_at(f64::midpoint(r3.0, r3.1), tol)?.distance(vd.0)
                    <= vd.1.max(tol.confusion() * 1e4)
            {
                ogeom_core::diag!(
                    Channel::Wire,
                    "DUST piece closing on one vertex at {:?} dropped from piece from_a={} face={}",
                    vd.0,
                    piece.from_a,
                    piece.face
                );
                continue;
            }
            edges.push(if traversal.reversed {
                built.reversed()
            } else {
                built
            });
        }
        // A ring every piece of which is dust stands within its vertices'
        // reach of one point: a hole of it is no hole, and a boundary of
        // it bounds no face.
        if edges.is_empty() {
            if wires.is_empty() {
                ogeom_core::diag!(
                    Channel::Wire,
                    "DUST piece from_a={} face={}: its boundary closes on one point",
                    piece.from_a,
                    piece.face
                );
                return Ok(None);
            }
            continue;
        }
        let wire = match make_wire(rebuild.model, &edges, tol) {
            Ok(w) => w.shape,
            Err(e) => {
                if diag::enabled(Channel::Wire) {
                    ogeom_core::diag!(
                        Channel::Wire,
                        "WIRE FAIL piece from_a={} face={}",
                        piece.from_a,
                        piece.face
                    );
                    for built in &edges {
                        if let Some((a, b)) = ogeom_algo::edge_vertices(rebuild.model, built)? {
                            let at = |v: &Shape| {
                                rebuild
                                    .model
                                    .node(v)
                                    .and_then(|n| n.data().as_vertex())
                                    .map(|d| (d.point, d.tolerance.get()))
                            };
                            ogeom_core::diag!(
                                Channel::Wire,
                                "   built {:?}{}: v{} {:?} -> v{} {:?}",
                                built.node(),
                                if built.orientation() == ogeom_topo::Orientation::Reversed {
                                    " rev"
                                } else {
                                    ""
                                },
                                a.node().index(),
                                at(&a),
                                b.node().index(),
                                at(&b)
                            );
                        }
                    }
                    for t in ring {
                        let tag = match &t.tag {
                            Tag::Boundary { edge, range } => format!(
                                "Boundary e{edge} {range:?} tol {:.2e}",
                                face.edges[*edge].tolerance
                            ),
                            Tag::Contact { contact, range } => format!(
                                "Contact c{contact} {range:?} tol {:.2e}",
                                fused.contacts[*contact].tolerance
                            ),
                            Tag::Section { section, range } => format!(
                                "Section s{section} {range:?} tol {:.2e}",
                                fused.sections[*section].tolerance
                            ),
                            Tag::Pole { pole, range } => format!("Pole p{pole} {range:?}"),
                        };
                        ogeom_core::diag!(Channel::Wire, "   {tag} reversed={}", t.reversed);
                    }
                }
                return Err(e);
            }
        };
        wires.push(wire);
    }
    let built = match make_face_on(rebuild.model, surface_id, &wires, tol) {
        Ok(b) => b.shape,
        Err(e) => {
            if diag::enabled(Channel::Wire) {
                ogeom_core::diag!(
                    Channel::Wire,
                    "FACE FAIL piece from_a={} face={}: {e}",
                    piece.from_a,
                    piece.face
                );
                for (wi, wire) in wires.iter().enumerate() {
                    for edge in rebuild.model.ordered_children_of(wire)? {
                        if let Some((a, b)) = ogeom_algo::edge_vertices(rebuild.model, &edge)? {
                            let at = |v: &Shape| {
                                rebuild
                                    .model
                                    .node(v)
                                    .and_then(|n| n.data().as_vertex())
                                    .map(|d| d.point)
                            };
                            ogeom_core::diag!(
                                Channel::Wire,
                                "   wire {wi} edge {:?}{}: v{} {:?} -> v{} {:?}",
                                edge.node(),
                                if edge.orientation() == ogeom_topo::Orientation::Reversed {
                                    " rev"
                                } else {
                                    ""
                                },
                                a.node().index(),
                                at(&a),
                                b.node().index(),
                                at(&b)
                            );
                        }
                    }
                }
            }
            return Err(e);
        }
    };
    Ok(Some(
        if face.face.orientation() == ogeom_topo::Orientation::Reversed {
            built.reversed()
        } else {
            built
        },
    ))
}

/// Build the exact sub-edge a tag names, pcurves attached, with the nodes
/// the face's arrangement closed its ring on at the strand's two ends.
pub(crate) fn build_sub_edge(
    rebuild: &mut Rebuild,
    fused: &GeneralFused,
    from_a: bool,
    face: &GFace,
    surface_id: ogeom_topo::SurfaceId,
    (tag, nodes): (&Tag, Option<[Point2; 2]>),
    tol: Tolerances,
) -> OgeomResult<Shape> {
    match tag {
        Tag::Pole { pole, range } => {
            // A pole rebuilds as what it was: one vertex, an edge with no
            // curve bounded by it twice, and the chart line that says where
            // it runs in this face's parameters.
            let p = &face.poles[*pole];
            let at = rebuild.vertex(p.point, tol);
            let model = &mut *rebuild.model;
            let mut data = ogeom_topo::EdgeData::new_in(tol);
            data.degenerate = true;
            let built = model.add_edge(data, &[at.clone(), at])?;
            ogeom_algo::attach_pcurve(
                model,
                &built,
                p.pcurve.clone(),
                surface_id,
                Location::identity(),
                *range,
            )?;
            Ok(built)
        }
        Tag::Boundary { edge, range } => {
            let e = &face.edges[*edge];
            let whole = e.other_side.is_none()
                && (range.0 - e.crange.0).abs() <= tol.parametric()
                && (range.1 - e.crange.1).abs() <= tol.parametric();
            // On the face's own surface, which other faces may stand on
            // too, an edge its node states as a seam stays one: the faces
            // either side each walk their own column of it, which two
            // plain pcurves on the one surface could not tell apart.
            let seam = match &e.seam {
                Some(pair) if e.other_side.is_none() && face.own_surface == Some(surface_id) => {
                    Some(pair)
                }
                _ => None,
            };
            if whole && let Some(built) = rebuild.whole_edges.get(&e.node).cloned() {
                let attached = rebuild
                    .model
                    .node(&built)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(|d| d.pcurve_on(surface_id).is_some());
                match seam {
                    // The face on the other side, on the same surface,
                    // attached the seam with the edge's first build.
                    Some(_) if attached => {}
                    Some((forward, reversed)) => ogeom_algo::attach_seam(
                        rebuild.model,
                        &built,
                        forward.clone(),
                        reversed.clone(),
                        surface_id,
                        Location::identity(),
                        e.prange,
                    )?,
                    None => ogeom_algo::attach_pcurve(
                        rebuild.model,
                        &built,
                        e.pcurve.clone(),
                        surface_id,
                        Location::identity(),
                        e.prange,
                    )?,
                }
                return Ok(built);
            }
            let from = e.curve.point_at(range.0, tol)?;
            let to = e.curve.point_at(range.1, tol)?;
            let v0 = rebuild.vertex(from, tol);
            let v1 = rebuild.vertex(to, tol);
            let model = &mut *rebuild.model;
            let built = make_edge_between(model, (*e.curve).clone(), *range, &v0, &v1, tol)?.shape;
            // A piece of a tolerant edge is the same curve with the same
            // honest radius: a fitted rail's stated slop must survive the
            // split, or the next boolean over this solid measures the rail
            // against a tolerance it never had. The ends own it too: a
            // junction on a tolerant curve is a junction to the curve's own
            // resolution, and the wires rebuilt through it meet within that.
            // A fuzzy boolean's confusion is its fuzz, which is no floor
            // here: a piece is born at the smallest tolerance an edge
            // states, and keeps whatever its edge stated above that.
            if e.tolerance > tol.least_tolerance() {
                if let Some(node) = model.node_mut(&built)
                    && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen_to(e.tolerance);
                }
                for v in [&v0, &v1] {
                    rebuild.widen(v, e.tolerance, tol)?;
                }
            }
            let model = &mut *rebuild.model;
            let sub_p = (
                rescale(range.0, e.crange, e.prange),
                rescale(range.1, e.crange, e.prange),
            );
            match (&e.other_side, seam) {
                (None, Some((forward, reversed))) => ogeom_algo::attach_seam(
                    model,
                    &built,
                    forward.clone(),
                    reversed.clone(),
                    surface_id,
                    Location::identity(),
                    sub_p,
                )?,
                (None, None) => ogeom_algo::attach_pcurve(
                    model,
                    &built,
                    e.pcurve.clone(),
                    surface_id,
                    Location::identity(),
                    sub_p,
                )?,
                (Some((other, _)), _) => {
                    // A seam: both sides attach over the one range a seam
                    // carries, and an occurrence picks its side by its
                    // orientation.
                    ogeom_algo::attach_seam(
                        model,
                        &built,
                        e.pcurve.clone(),
                        other.clone(),
                        surface_id,
                        Location::identity(),
                        sub_p,
                    )?;
                }
            }
            // An edge's pcurve need only lie on its curve within the
            // tolerance, not keep pace with it: lifted, it may run a little
            // ahead or behind. Cut short, the piece's end is read at one
            // parameter on both, and that lag stands between the piece's
            // vertex and its image in the chart. The piece states it.
            if !whole
                && let Some(gap) = ogeom_algo::edge_pcurve_gap(model, &built, tol)?
                && let Some(data) = model.node(&built).and_then(|n| n.data().as_edge())
                && gap > data.tolerance.get()
            {
                let widened = ogeom_core::Tolerance::new_in(gap * (1.0 + 1e-6), tol)?;
                if let Some(node) = model.node_mut(&built)
                    && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen(widened);
                }
                for v in [&v0, &v1] {
                    rebuild.widen(v, widened.get(), tol)?;
                }
            }
            if whole {
                rebuild.whole_edges.insert(e.node, built.clone());
            }
            Ok(built)
        }
        Tag::Contact { contact, range } => {
            let c = &fused.contacts[*contact];
            let from = c.curve.point_at(range.0, tol)?;
            let to = c.curve.point_at(range.1, tol)?;
            let v0 = rebuild.vertex(from, tol);
            let v1 = rebuild.vertex(to, tol);
            let model = &mut *rebuild.model;
            let built = make_edge_between(model, c.curve.clone(), *range, &v0, &v1, tol)?.shape;
            // A tolerant contact's pieces and ends own its stated slop, as a
            // boundary's and a section's do: three kinds of strand close one
            // ring, and every junction meets at the strands' own honesty.
            if c.tolerance > tol.least_tolerance() {
                if let Some(node) = model.node_mut(&built)
                    && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen_to(c.tolerance);
                }
                for v in [&v0, &v1] {
                    rebuild.widen(v, c.tolerance, tol)?;
                }
            }
            let model = &mut *rebuild.model;
            // The stored image keeps its own window. The attached copy names
            // the sub-window this piece covers under the proportional map.
            let sub_p = (
                rescale(range.0, c.crange, c.prange),
                rescale(range.1, c.crange, c.prange),
            );
            let mid = c.pcurve.point_at(f64::midpoint(sub_p.0, sub_p.1), tol)?;
            let trim: Vec<&[Point2]> = face.trim_lines(tol).iter().map(Vec::as_slice).collect();
            let folded = fold_inside(mid, &face.surface, &trim, tol);
            let shifted = c
                .pcurve
                .transformed(&ogeom_math::Transform2::translation(folded - mid), tol)?;
            ogeom_algo::attach_pcurve(
                model,
                &built,
                shifted,
                surface_id,
                Location::identity(),
                sub_p,
            )?;
            // A contact's curve lies on the target's surface only to the
            // weld or the fuzz (a rim on a drum a micron off its own), and
            // its image there, lifted, stands off the curve by as much. The
            // piece states that, and its ends hold it.
            if let Some(gap) = ogeom_algo::edge_pcurve_gap(model, &built, tol)?
                && let Some(data) = model.node(&built).and_then(|n| n.data().as_edge())
                && gap > data.tolerance.get()
            {
                let widened = ogeom_core::Tolerance::new_in(gap * (1.0 + 1e-6), tol)?;
                if let Some(node) = model.node_mut(&built)
                    && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen(widened);
                }
                for v in [&v0, &v1] {
                    rebuild.widen(v, widened.get(), tol)?;
                }
            }
            Ok(built)
        }
        Tag::Section { section, range } => {
            let s = &fused.sections[*section];
            let domain = s.curve.domain();
            let (f0, f1) = folded_range(*range, domain, s.closed);
            let from = s.curve.point_at(at_param(f0, domain, s.closed), tol)?;
            let to = s.curve.point_at(at_param(f1, domain, s.closed), tol)?;
            let v0 = rebuild.vertex(from, tol);
            let v1 = rebuild.vertex(to, tol);
            let model = &mut *rebuild.model;
            let built = make_edge_between(model, s.curve.clone(), (f0, f1), &v0, &v1, tol)?.shape;
            // A piece welded into a junction some way off its own end (a
            // section tangent to an edge at the edge's vertex, stopped
            // where it began to hug the edge) owns that gap as the vertex
            // does: the face it bounds meets its neighbour there only that
            // closely, in its chart as in space.
            let at = |v: &Shape| {
                model
                    .node(v)
                    .and_then(|n| n.data().as_vertex())
                    .map(|d| d.point)
            };
            let ends = [at(&v0), at(&v1)];
            let welded = [(ends[0], from), (ends[1], to)]
                .iter()
                .filter_map(|(v, p)| v.map(|v| v.distance(*p)))
                .fold(0.0, f64::max);
            let own = s.tolerance.max(welded);
            // A fitted section's pieces and their ends own the section's
            // stated slop, exactly as a tolerant boundary's do: the ring
            // they close alternates between the two, and both sides must
            // meet at the junction's own resolution.
            if own > tol.least_tolerance() {
                if let Some(node) = model.node_mut(&built)
                    && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen_to(own);
                }
                for v in [&v0, &v1] {
                    rebuild.widen(v, own, tol)?;
                }
            }
            let model = &mut *rebuild.model;
            // The section's pcurve is unwrapped across any seam. The face's
            // triangulator lives in one chart, so the attached copy is folded
            // home by the same period shift the arrangement gave this
            // strand's polyline (to the side of the face's seam its trim is
            // on), decided by the sub-range's midpoint, so an
            // endpoint sitting exactly on the chart's edge stays on the side
            // the arc's body is.
            let pcurve = if from_a { &s.pc_a } else { &s.pc_b };
            let mid = pcurve.point_at(f64::midpoint(f0, f1), tol)?;
            let trim: Vec<&[Point2]> = face.trim_lines(tol).iter().map(Vec::as_slice).collect();
            let folded = fold_inside(mid, &face.surface, &trim, tol);
            let shifted =
                pcurve.transformed(&ogeom_math::Transform2::translation(folded - mid), tol)?;
            // The chart image ends where the vertex is, as the edge does in
            // space to its tolerance: a welded end left where the section
            // stopped leaves the face's outline open by the weld in its
            // chart, wider there than any later arrangement's snap.
            let loose = s.tolerance.max(tol.confusion() * 1e2);
            let bent = [
                rebuild.onto_vertex.contains(&v0.node()),
                rebuild.onto_vertex.contains(&v1.node()),
            ];
            let targets = [
                ends[0].filter(|v| bent[0] && v.distance(from) > loose),
                ends[1].filter(|v| bent[1] && v.distance(to) > loose),
            ];
            // Elsewhere the image ends on the node the arrangement closed
            // the ring on, where its own end stands further off it than any
            // later arrangement's snap and the node lifts nearer the piece's
            // end in space than the image's end does. A fitted image can lag
            // its curve by a few microns. Beside a pole, where a micron spans
            // a hundred times as much of the chart, that leaves it a
            // ten-thousandth of a radian off the seam column its neighbour
            // ends on, and the face's ring open there to every later
            // boolean. A piece whose own end stops short of its node (a loop
            // whose last hair was welded away) keeps the image of that end.
            let mut onto = [None, None];
            for (k, ((target, t), at)) in targets.iter().zip([f0, f1]).zip([from, to]).enumerate() {
                let end = shifted.point_at(t, tol)?;
                if let Some(target) = target {
                    onto[k] = Some(chart_foot(&face.surface, *target, end, tol)?);
                    continue;
                }
                let Some(node) = nodes.map(|n| n[k]) else {
                    continue;
                };
                let node = period_shifts(&face.surface, tol)
                    .into_iter()
                    .map(|(du, dv)| Point2::new(node.x + du, node.y + dv))
                    .min_by(|p, q| p.distance(end).total_cmp(&q.distance(end)))
                    .unwrap_or(node);
                let vertex = if k == 0 { &v0 } else { &v1 };
                let radius = model.tolerance_of(vertex)?.map_or(0.0, |r| r.get());
                if node.distance(end) > param_snap(&face.surface, tol)
                    && !on_collapsed_row(&face.surface, node, radius, tol)
                    && face.surface.point_at(node.x, node.y, tol)?.distance(at)
                        < face.surface.point_at(end.x, end.y, tol)?.distance(at)
                {
                    onto[k] = Some(node);
                }
            }
            let shifted = if onto.iter().any(Option::is_some) {
                pcurve_onto_ends(&shifted, (f0, f1), &face.surface, onto, tol)?
            } else {
                shifted
            };
            ogeom_algo::attach_pcurve(
                model,
                &built,
                shifted,
                surface_id,
                Location::identity(),
                (f0, f1),
            )?;
            // The section's image was fitted through samples and may have
            // been bent onto its welded ends; lifted, it stands off the
            // curve by what is measured along the whole piece, and the
            // piece and its ends state that.
            if let Some(gap) = ogeom_algo::edge_pcurve_gap(model, &built, tol)?
                && let Some(data) = model.node(&built).and_then(|n| n.data().as_edge())
                && gap > data.tolerance.get()
            {
                let widened = ogeom_core::Tolerance::new_in(gap * (1.0 + 1e-6), tol)?;
                if let Some(node) = model.node_mut(&built)
                    && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen(widened);
                }
                for v in [&v0, &v1] {
                    rebuild.widen(v, widened.get(), tol)?;
                }
            }
            Ok(built)
        }
    }
}

/// The point of `surface`'s chart nearest `end` that lifts to the foot of
/// `target` on the surface: the foot may land a turn away from `end` on a
/// periodic chart, whether or not the surface's wrapper says it is
/// periodic, and a shifted candidate counts where the surface agrees it is
/// the foot.
pub(crate) fn chart_foot(
    surface: &SurfaceGeometry,
    target: Point,
    end: Point2,
    tol: Tolerances,
) -> OgeomResult<Point2> {
    let foot = ogeom_algo::project_on_surface_from(surface, target, (end.x, end.y), tol)?;
    let raw = Point2::new(foot.parameters.0, foot.parameters.1);
    let turn = core::f64::consts::TAU;
    let mut shifts = period_shifts(surface, tol);
    for du in [-turn, 0.0, turn] {
        for dv in [-turn, 0.0, turn] {
            shifts.push((du, dv));
        }
    }
    Ok(shifts
        .into_iter()
        .map(|(du, dv)| Point2::new(raw.x + du, raw.y + dv))
        .filter(|c| {
            surface
                .point_at(c.x, c.y, tol)
                .is_ok_and(|q| q.distance(foot.point) <= tol.confusion() * 1e2)
        })
        .min_by(|a, b| {
            a.distance(end)
                .partial_cmp(&b.distance(end))
                .unwrap_or(core::cmp::Ordering::Equal)
        })
        .unwrap_or(raw))
}

/// Whether the chart row or column through `at` lies, across the chart's
/// whole width that way, within `radius` of one point in space: a pole's
/// row or an apex's, along which a point moves nowhere, so a pcurve's end
/// is placed on it by the way it arrives and not by a node there.
pub(crate) fn on_collapsed_row(
    surface: &SurfaceGeometry,
    at: Point2,
    radius: f64,
    tol: Tolerances,
) -> bool {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let Ok((du, dv)) = surface.d1_at(at.x, at.y, tol) else {
        return true;
    };
    du.magnitude() * (u1 - u0).abs() <= radius || dv.magnitude() * (v1 - v0).abs() <= radius
}

/// A pcurve over `range` bent at its ends onto the chart points `targets`,
/// where given, and refitted at the same parameters, at every span of the
/// image, so the image stays same-parameter with its edge to the
/// correction's size.
///
/// Each end's correction fades linearly to nothing at the other end, and,
/// read on the surface's first derivatives, nowhere moves the image
/// further in space than it moves the end: where a chart step costs more
/// in space than it does at the end (away from a pole the end stands
/// beside), the step is scaled down by as much. A ten-thousandth of a
/// radian round a ball's axis is a few microns beside its pole and a
/// ten-thousandth of its radius at its equator.
pub(crate) fn pcurve_onto_ends(
    pcurve: &PlanarCurve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    targets: [Option<Point2>; 2],
    tol: Tolerances,
) -> OgeomResult<PlanarCurve> {
    const SAMPLES: usize = 32;
    let ends = [
        pcurve.point_at(range.0, tol)?,
        pcurve.point_at(range.1, tol)?,
    ];
    // The length in space of a chart step `by` taken at `at`.
    let cost = |at: Point2, by: ogeom_math::Vector2| -> Option<f64> {
        let (du, dv) = surface.d1_at(at.x, at.y, tol).ok()?;
        Some((du * by.x + dv * by.y).magnitude())
    };
    let mut deltas = [ogeom_math::Vector2::new(0.0, 0.0); 2];
    let mut at_end = [None; 2];
    for (k, target) in targets.iter().enumerate() {
        let Some(target) = target else {
            continue;
        };
        deltas[k] = *target - ends[k];
        at_end[k] = cost(ends[k], deltas[k]);
    }
    // How much of an end's correction a point of the image takes, its
    // linear share scaled so the point moves no further in space than the
    // end does.
    let share = |k: usize, p: Point2, linear: f64| -> f64 {
        let (Some(end), Some(here)) = (at_end[k], cost(p, deltas[k])) else {
            return linear;
        };
        if here > end {
            linear * end / here
        } else {
            linear
        }
    };
    // Refitted at every span of the image and held to the target between
    // the samples as well as at them; the edge states what the bend and
    // the fit leave once the image is attached.
    let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
    let stations = ogeom_algo::traced::stations_2d(pcurve, (lo, hi), SAMPLES);
    let fitted = ogeom_algo::traced::fit_traced_2d(
        |t| {
            let along = if hi > lo {
                (t - range.0) / (range.1 - range.0)
            } else {
                0.0
            };
            let p = pcurve.point_at(t, tol)?;
            Ok(p + deltas[0] * share(0, p, 1.0 - along) + deltas[1] * share(1, p, along))
        },
        &stations,
        3,
        tol.confusion() * 10.0,
        tol,
    )?;
    Ok(PlanarCurve::BSpline(fitted.curve))
}

/// Sew kept pieces, demand closure, and nest shells into solids and voids.
///
/// The faces each solid set aside come through as they are where
/// `keep_aside` keeps that solid's outside, and are deleted where it does
/// not.
pub(crate) fn assemble_result(
    model: &mut Model,
    fused: &GeneralFused,
    kept: &[(usize, bool)],
    keep_aside: [bool; 2],
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
) -> OgeomResult<Built> {
    ogeom_core::progress::stage("boolean: assemble");
    let local = !fused.a.aside.is_empty() || !fused.b.aside.is_empty();
    let mut history = History::new();
    let source_face = |piece: &FacePiece| -> Shape {
        if piece.from_a {
            fused.a.faces[piece.face].face.clone()
        } else {
            fused.b.faces[piece.face].face.clone()
        }
    };

    if kept.is_empty() {
        if (keep_aside[0] && !fused.a.aside.is_empty())
            || (keep_aside[1] && !fused.b.aside.is_empty())
        {
            ogeom_bail!(NotDone, "no piece is kept beside the faces set aside");
        }
        // A legitimate answer: cutting a solid away entirely leaves nothing.
        let empty = model.add_compound(&[])?;
        for piece in &fused.pieces {
            history.delete(&source_face(piece));
        }
        for (face, _) in fused.a.aside.iter().chain(&fused.b.aside) {
            history.delete(face);
        }
        history.modify(a, empty.clone());
        history.modify(b, empty.clone());
        return Ok(Built::new(empty, history));
    }

    let Rebuilt {
        faces,
        settled,
        weld: floor,
        holes,
        kept: built_from,
    } = rebuilt_pieces(model, fused, kept, &mut history, tol)?;
    let sewn = sew_around_within(model, &faces, &settled, SEW_ACROSS.get(), tol)?;
    // Sewing rebuilds the pieces onto shared edges: the result's faces are
    // its faces, reached from the inputs through both steps.
    history = history.followed_by(&sewn.history);
    // What the sew and the passes after it made of each rebuilt face, where
    // faces set aside are put back beside them.
    let mut after = if local {
        sewn.history.clone()
    } else {
        History::new()
    };
    let mut sewn = sewn;
    let dropped = without_membranes(model, &mut sewn.shells, floor, tol)?;
    history = history.followed_by(&dropped);
    let held: ogeom_core::FastSet<ogeom_topo::TShapeId> = fused
        .a
        .shared
        .values()
        .chain(fused.b.shared.values())
        .map(Shape::node)
        .collect();
    let joined = seam_join::join_across_seams(model, &mut sewn.shells, &held, tol)?;
    history = history.followed_by(&joined).without_repeated_images();
    let mut closed = None;
    let mut made = None;
    if local {
        after = after.followed_by(&dropped).followed_by(&joined);
        let (shells, shut, new) = with_aside(
            model,
            fused,
            &built_from,
            &faces,
            &holes,
            &after,
            &sewn.shells,
            keep_aside,
            &mut history,
            tol,
        )?;
        sewn.shells = shells;
        closed = Some(shut);
        made = Some(new);
    }
    if sewn.shells.is_empty() {
        // Membranes all through: what was kept encloses nothing.
        let empty = model.add_compound(&[])?;
        history.modify(a, empty.clone());
        history.modify(b, empty.clone());
        return Ok(Built::new(empty, history));
    }
    for (i, shell) in sewn.shells.iter().enumerate() {
        let shut = match &closed {
            Some(closed) => closed[i],
            None => is_shell_closed(model, shell)?,
        };
        if !shut {
            // Env-gated forensics: the open shell's unshared edges, the
            // question every failure here starts from.
            if diag::enabled(Channel::Arrange) {
                use ogeom_geom::Curve3d as _;
                for edge in ogeom_topo::explore_unique(model, shell, ShapeType::Edge)? {
                    let users = ogeom_topo::explore(model, shell, Filter::OfType(ShapeType::Face))?
                        .iter()
                        .filter(|f| {
                            ogeom_topo::explore_unique(model, f, ShapeType::Edge)
                                .map(|es| es.iter().any(|e2| e2.node() == edge.node()))
                                .unwrap_or(false)
                        })
                        .count();
                    let mut occurrences = 0_usize;
                    for f in ogeom_topo::explore(model, shell, Filter::OfType(ShapeType::Face))? {
                        for wire in model.children_of(&f)? {
                            for e2 in model.children_of(&wire)? {
                                if e2.node() == edge.node() {
                                    occurrences += 1;
                                }
                            }
                        }
                    }
                    if occurrences % 2 == 1
                        && let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge())
                        && let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) =
                            data.curve3d()
                        && let Some(g) = model.geometry().curve(*curve)
                    {
                        let a = g.point_at(range.0, tol)?;
                        let b = g.point_at(range.1, tol)?;
                        ogeom_core::diag!(
                            Channel::Arrange,
                            "  open edge uses={occurrences} faces={users} {:?} range={range:?}: {a:?} -> {b:?}",
                            core::mem::discriminant(g)
                        );
                        // Twins: other edges of the shell with the same ends,
                        // which the sew should have found to be this edge.
                        for other in ogeom_topo::explore_unique(model, shell, ShapeType::Edge)? {
                            if other.node() == edge.node() {
                                continue;
                            }
                            let Some(od) = model.node(&other).and_then(|n| n.data().as_edge())
                            else {
                                continue;
                            };
                            let Some(ogeom_topo::EdgeRepr::Curve3d {
                                curve: oc,
                                range: orange,
                                ..
                            }) = od.curve3d()
                            else {
                                continue;
                            };
                            let Some(og) = model.geometry().curve(*oc) else {
                                continue;
                            };
                            let oa = og.point_at(orange.0, tol)?;
                            let ob = og.point_at(orange.1, tol)?;
                            let same = (oa.distance(a) <= 1e-3 && ob.distance(b) <= 1e-3)
                                || (oa.distance(b) <= 1e-3 && ob.distance(a) <= 1e-3);
                            if same {
                                let om = og.point_at(f64::midpoint(orange.0, orange.1), tol)?;
                                let m = g.point_at(f64::midpoint(range.0, range.1), tol)?;
                                ogeom_core::diag!(
                                    Channel::Arrange,
                                    "    twin {:?} range={orange:?} mid gap {:.2e} tol {:.2e} vs {:.2e}",
                                    core::mem::discriminant(og),
                                    om.distance(m),
                                    od.tolerance.get(),
                                    data.tolerance.get()
                                );
                            }
                        }
                        for f in ogeom_topo::explore(model, shell, Filter::OfType(ShapeType::Face))?
                        {
                            let uses_it = ogeom_topo::explore_unique(model, &f, ShapeType::Edge)
                                .map(|es| es.iter().any(|e2| e2.node() == edge.node()))
                                .unwrap_or(false);
                            if !uses_it {
                                continue;
                            }
                            if let Some(ogeom_topo::NodeData::Face(fd)) =
                                model.node(&f).map(|n| n.data())
                                && let Some(sg) = model.geometry().surface(fd.surface)
                            {
                                let bound = shape_bounds(model, &f, tol)?;
                                ogeom_core::diag!(
                                    Channel::Arrange,
                                    "    used by face on {:?} bound {:?}",
                                    core::mem::discriminant(sg),
                                    bound
                                );
                            }
                        }
                    }
                }
            }
            if diag::enabled(Channel::Arrange) {
                use ogeom_geom::Curve3d as _;
                for f in ogeom_topo::explore(model, shell, Filter::OfType(ShapeType::Face))? {
                    let Some(ogeom_topo::NodeData::Face(fd)) = model.node(&f).map(|n| n.data())
                    else {
                        continue;
                    };
                    let Some(sg) = model.geometry().surface(fd.surface) else {
                        continue;
                    };
                    let edges = ogeom_topo::explore_unique(model, &f, ShapeType::Edge)?;
                    ogeom_core::diag!(
                        Channel::Arrange,
                        "  face {:?} with {} edges",
                        core::mem::discriminant(sg),
                        edges.len()
                    );
                    if matches!(sg, SurfaceGeometry::Sphere(_)) {
                        for e in &edges {
                            if let Some(ed) = model.node(e).and_then(|n| n.data().as_edge())
                                && let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) =
                                    ed.curve3d()
                                && let Some(g) = model.geometry().curve(*curve)
                            {
                                let a = g.point_at(range.0, tol)?;
                                let b = g.point_at(range.1, tol)?;
                                ogeom_core::diag!(
                                    Channel::Arrange,
                                    "    sphere edge {:?} ({:.4},{:.4},{:.4}) -> ({:.4},{:.4},{:.4})",
                                    core::mem::discriminant(g),
                                    a.x,
                                    a.y,
                                    a.z,
                                    b.x,
                                    b.y,
                                    b.z
                                );
                            }
                        }
                    }
                }
            }
            ogeom_bail!(
                NotDone,
                "the kept pieces did not close into a shell; the configuration \
                 is beyond what the boolean currently resolves"
            );
        }
    }

    // Nest: a shell inside another is that solid's void, and a shell inside
    // a void is a solid of its own. The kept pieces are already oriented
    // as the material's boundary, so a void's faces face into it as they
    // are.
    let mut solids = Vec::new();
    for (outer, voids) in nest_shells(model, &sewn.shells, tol)? {
        let mut group = vec![sewn.shells[outer].clone()];
        group.extend(voids.iter().map(|&j| sewn.shells[j].clone()));
        solids.push(model.add_solid(&group)?);
    }
    let result = if solids.len() == 1 {
        solids.remove(0)
    } else {
        model.add_compound(&solids)?
    };
    // A face made here whose boundary does not close in its surface's chart
    // (a section stopped short of where its neighbour starts, at a contact
    // the marcher did not follow through) holds a strand that dangles, and
    // the result would check broken. It is refused instead, for the retries
    // to take another way.
    // Set aside, the faces made are known; otherwise they are the result's
    // faces that are not the operands'.
    let built_faces = match &made {
        Some(made) => made.clone(),
        None => {
            let given: ogeom_core::FastSet<ogeom_topo::TShapeId> = [a, b]
                .into_iter()
                .map(|operand| explore_unique(model, operand, ShapeType::Face))
                .collect::<OgeomResult<Vec<_>>>()?
                .into_iter()
                .flatten()
                .map(|face| face.node())
                .collect();
            explore_unique(model, &result, ShapeType::Face)?
                .into_iter()
                .filter(|f| !given.contains(&f.node()))
                .collect()
        }
    };
    // A wire of a face the boolean gathered that comes through as the same
    // wire (a hole clear of the other solid, put back into the piece that
    // holds it) closes as it did.
    let given_wires: ogeom_core::FastSet<ogeom_topo::TShapeId> = fused
        .a
        .faces
        .iter()
        .chain(&fused.b.faces)
        .map(|gathered| model.ordered_children_of(&gathered.face))
        .collect::<OgeomResult<Vec<_>>>()?
        .into_iter()
        .flatten()
        .map(|wire| wire.node())
        .collect();
    let taken = |wire: &Shape| given_wires.contains(&wire.node());
    for face in &built_faces {
        if !ogeom_algo::chart_breaks(model, face, &taken, tol)?.is_empty() {
            ogeom_bail!(
                NotDone,
                "a face the boolean made does not close in its surface's chart; the \
                 configuration is beyond what the boolean currently resolves"
            );
        }
    }
    // The faces set aside keep the boxes they had: only the faces made
    // are looked at.
    match &made {
        Some(made) => keep_face_boxes(model, made)?,
        None => keep_face_boxes(model, &explore_unique(model, &result, ShapeType::Face)?)?,
    }
    history.modify(a, result.clone());
    history.modify(b, result.clone());
    Ok(Built::new(result, history))
}

/// Each of `faces` on a ruled surface (a plane, a drum, a cone) keeps its
/// box, found from its edges alone: the next boolean on the result reads
/// every face's box to find the faces it can set aside.
pub(crate) fn keep_face_boxes(model: &Model, faces: &[Shape]) -> OgeomResult<()> {
    for face in faces {
        if model.face_bounds(face).is_some() {
            continue;
        }
        let ruled = model
            .node(face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| {
                matches!(
                    s,
                    SurfaceGeometry::Plane(_)
                        | SurfaceGeometry::Cylinder(_)
                        | SurfaceGeometry::Cone(_)
                )
            });
        if ruled {
            ogeom_algo::face_bounds(model, face)?;
        }
    }
    Ok(())
}

/// The kept pieces of a sheet sewn along the edges they share: one shell,
/// or a compound of shells where they fall apart. Nothing is asked of the
/// shells' closure. A sheet that is one face has its history in its
/// pieces; a shell or a compound also becomes the result.
pub(crate) fn assemble_sheet(
    model: &mut Model,
    fused: &GeneralFused,
    kept: &[usize],
    sheet: &Shape,
    tol: Tolerances,
) -> OgeomResult<Built> {
    ogeom_core::progress::stage("boolean: assemble");
    let mut history = History::new();
    let whole = model.kind_of(sheet)? != ShapeType::Face;
    if kept.is_empty() {
        let empty = model.add_compound(&[])?;
        for piece in &fused.pieces {
            history.delete(&fused.a.faces[piece.face].face);
        }
        if whole {
            history.modify(sheet, empty.clone());
        }
        return Ok(Built::new(empty, history));
    }
    let kept: Vec<(usize, bool)> = kept.iter().map(|&i| (i, false)).collect();
    let Rebuilt { faces, settled, .. } = rebuilt_pieces(model, fused, &kept, &mut history, tol)?;
    let mut sewn = sew_around(model, &faces, &settled, tol)?;
    history = history.followed_by(&sewn.history);
    let joined = seam_join::join_across_seams(
        model,
        &mut sewn.shells,
        &ogeom_core::FastSet::default(),
        tol,
    )?;
    history = history.followed_by(&joined).without_repeated_images();
    let result = if sewn.shells.len() == 1 {
        sewn.shells.remove(0)
    } else {
        model.add_compound(&sewn.shells)?
    };
    if whole {
        history.modify(sheet, result.clone());
    }
    Ok(Built::new(result, history))
}

/// The kept pieces rebuilt as faces, each recorded in `history` against the
/// face it came from and every face no piece is kept of deleted, whether
/// each is a copy already sewn to all its neighbours, and the distance the
/// rebuild welds their ends within.
pub(crate) fn rebuilt_pieces(
    model: &mut Model,
    fused: &GeneralFused,
    kept: &[(usize, bool)],
    history: &mut History,
    tol: Tolerances,
) -> OgeomResult<Rebuilt> {
    let source_face = |piece: &FacePiece| -> Shape {
        if piece.from_a {
            fused.a.faces[piece.face].face.clone()
        } else {
            fused.b.faces[piece.face].face.clone()
        }
    };
    let weld = fused
        .contacts
        .iter()
        .map(|c| c.tolerance)
        .chain(fused.sections.iter().map(|s| s.tolerance))
        .fold(0.0_f64, |acc, t| acc.max(honest(t, tol) * 3.0));
    let floor = weld.max(tol.confusion() * 1e2);
    let junction_reach = fused
        .junctions
        .iter()
        .fold(0.0_f64, |acc, j| acc.max(j.reach));
    let mut junction_bins = bins::Bins::new(junction_reach + floor);
    for (index, j) in fused.junctions.iter().enumerate() {
        junction_bins.insert(j.at, index);
    }
    // The vertices' cells are sized to the loosest tolerance a strand
    // brings to its ends, so an end is compared across a few cells, not
    // every vertex.
    let loosest = fused
        .a
        .faces
        .iter()
        .chain(fused.b.faces.iter())
        .flat_map(|f| f.edges.iter().map(|e| e.tolerance))
        .chain(fused.contacts.iter().map(|c| c.tolerance))
        .chain(fused.sections.iter().map(|s| s.tolerance))
        .fold(junction_reach, f64::max);
    let mut rebuild = Rebuild {
        model,
        surfaces_a: vec![None; fused.a.faces.len()],
        surfaces_b: vec![None; fused.b.faces.len()],
        vertices: Vec::new(),
        vertex_bins: bins::Bins::new(floor + loosest),
        near: Vec::new(),
        widest: 0.0,
        weld,
        junctions: fused.junctions.iter().map(|j| (*j, None)).collect(),
        junction_bins,
        junction_reach,
        onto_vertex: ogeom_core::FastSet::default(),
        whole_edges: ogeom_core::FastMap::default(),
        kept: ogeom_core::FastMap::default(),
        strained: false,
    };
    for (at, vertex) in fused
        .a
        .shared_vertices
        .iter()
        .chain(&fused.b.shared_vertices)
    {
        rebuild.remember(*at, vertex);
        rebuild.kept.insert(vertex.node(), *at);
    }
    // A face the operation never touched (nothing crosses it, no edge of it
    // is split, no junction lands at its corners) is rebuilt as an exact
    // copy of itself. Its neighbours across its unsplit edges build those
    // edges from the same curves, so whichever twin the sew keeps is its
    // own boundary.
    let clean = |piece: &FacePiece, model: &Model| -> OgeomResult<bool> {
        if !piece.whole {
            return Ok(false);
        }
        let face = &source_face(piece);
        for vertex in explore_unique(model, face, ShapeType::Vertex)? {
            let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) else {
                continue;
            };
            let at = vertex.transform(model.datums())?.apply(data.point);
            if fused
                .junctions
                .iter()
                .any(|j| j.at.distance(at) <= j.reach + floor)
            {
                return Ok(false);
            }
        }
        Ok(true)
    };
    // Each kept piece built first: one whose boundary closes on a point is
    // dust, and none of what follows sees it.
    let mut pieces_built: Vec<(usize, bool, Shape)> = Vec::with_capacity(kept.len());
    for &(index, flip) in kept {
        if let Some(built) = build_piece(&mut rebuild, fused, &fused.pieces[index], tol)? {
            pieces_built.push((index, flip, built));
        }
    }
    let kept: Vec<(usize, bool)> = pieces_built.iter().map(|&(i, f, _)| (i, f)).collect();
    let kept = kept.as_slice();
    // One kept turned over (a tool's face lining the cavity a cut leaves)
    // faces the other way, which no copy does.
    let mut copies = Vec::with_capacity(kept.len());
    for &(index, flip) in kept {
        copies.push(!flip && clean(&fused.pieces[index], rebuild.model)?);
    }
    let mut faces = Vec::new();
    let mut kept_sources: ogeom_core::FastSet<Shape> = ogeom_core::FastSet::default();
    let mut built_of: ogeom_core::FastMap<usize, usize> = ogeom_core::FastMap::default();
    for (slot, (index, flip, built)) in pieces_built.into_iter().enumerate() {
        built_of.insert(index, slot);
        let piece = &fused.pieces[index];
        let mut built = built;
        if flip {
            built = built.reversed();
        }
        if copies[slot] {
            history.copy(&source_face(piece), built.clone());
        } else {
            history.modify(&source_face(piece), built.clone());
        }
        kept_sources.insert(source_face(piece));
        faces.push(built);
    }
    if rebuild.strained {
        ogeom_bail!(
            NotDone,
            "a strand ends further from a vertex shared with a face set aside \
             than the vertex reaches"
        );
    }
    let holes = holes_put_back(rebuild.model, fused, kept, tol)?;
    // A piece dropped because the other argument's coincident piece is
    // kept for the same patch: its face became that piece's face.
    let mut stood_for: ogeom_core::FastSet<Shape> = ogeom_core::FastSet::default();
    for &(index_a, index_b) in &fused.coincident {
        let (dropped, standing) = match (built_of.get(&index_a), built_of.get(&index_b)) {
            (Some(&slot), None) => (index_b, slot),
            (None, Some(&slot)) => (index_a, slot),
            _ => continue,
        };
        let source = source_face(&fused.pieces[dropped]);
        if !history
            .modified(&source)
            .iter()
            .any(|f| f.is_same(&faces[standing]))
        {
            history.modify(&source, faces[standing].clone());
        }
        stood_for.insert(source);
    }
    for piece in &fused.pieces {
        let source = source_face(piece);
        if !kept_sources.contains(&source) && !stood_for.contains(&source) {
            history.delete(&source);
        }
    }

    // A copy shares each edge with the face of its own solid across it,
    // built once for both: no edge a copy holds has a twin in its own
    // solid. Where it meets a face of the other solid, its edges may have
    // twins only the sew can find, so only a copy all of whose vertices are
    // held by faces of its own solid alone comes through the sew as it
    // stands.
    let class = |slot: usize| -> u8 {
        match (copies[slot], fused.pieces[kept[slot].0].from_a) {
            (true, true) => 1,
            (true, false) => 2,
            (false, true) => 4,
            (false, false) => 8,
        }
    };
    let mut held: ogeom_core::FastMap<ogeom_topo::TShapeId, u8> = ogeom_core::FastMap::default();
    let mut corners = Vec::with_capacity(faces.len());
    for (slot, face) in faces.iter().enumerate() {
        let vertices = explore_unique(rebuild.model, face, ShapeType::Vertex)?;
        for vertex in &vertices {
            *held.entry(vertex.node()).or_default() |= class(slot);
        }
        corners.push(vertices);
    }
    let settled = corners
        .iter()
        .enumerate()
        .map(|(slot, vertices)| {
            let own = if fused.pieces[kept[slot].0].from_a {
                1 | 4
            } else {
                2 | 8
            };
            copies[slot]
                && vertices
                    .iter()
                    .all(|v| held.get(&v.node()).is_some_and(|by| by & !own == 0))
        })
        .collect();

    Ok(Rebuilt {
        faces,
        settled,
        weld: floor,
        holes,
        kept: kept.to_vec(),
    })
}

/// The shells the result is bounded by, with the faces set aside put back:
/// each face taking back holes left out of its arrangement is made again
/// with them, the faces set aside join the faces sewn where `keep_aside`
/// keeps them, and the shells are found again from the edges they all
/// share. `history` records the faces set aside and the faces made again.
#[allow(
    clippy::too_many_arguments,
    reason = "the assembly's own state, read once"
)]
pub(crate) fn with_aside(
    model: &mut Model,
    fused: &GeneralFused,
    kept: &[(usize, bool)],
    faces: &[Shape],
    holes: &[Vec<usize>],
    after: &History,
    shells: &[Shape],
    keep_aside: [bool; 2],
    history: &mut History,
    tol: Tolerances,
) -> OgeomResult<(Vec<Shape>, Vec<bool>, Vec<Shape>)> {
    if kept.len() != faces.len() {
        ogeom_bail!(
            Construction,
            "{} kept pieces for {} rebuilt faces; each face is read with its own piece",
            kept.len(),
            faces.len()
        );
    }
    // Each face made again with its holes, with the wires it walks itself
    // and the faces set aside beside its holes.
    let mut remade: ogeom_core::FastMap<ogeom_topo::TShapeId, Remade> =
        ogeom_core::FastMap::default();
    let mut put_back = History::new();
    // Each face's place: its solid, the place there of the face it came
    // from, and its piece's.
    let mut rank: ogeom_core::FastMap<ogeom_topo::TShapeId, (bool, usize, usize)> =
        ogeom_core::FastMap::default();
    for (slot, built) in faces.iter().enumerate() {
        let piece = &fused.pieces[kept[slot].0];
        let own = if piece.from_a { &fused.a } else { &fused.b };
        let source = &own.faces[piece.face].face;
        let place = (
            !piece.from_a,
            own.order
                .get(&ogeom_topo::SameKey(source.clone()))
                .copied()
                .unwrap_or(usize::MAX),
            slot,
        );
        rank.insert(built.node(), place);
        for image in after.modified(built) {
            rank.insert(image.node(), place);
        }
    }
    for (slot, (built, holes)) in faces.iter().zip(holes).enumerate() {
        if holes.is_empty() {
            continue;
        }
        let image = match after.modified(built) {
            [] if !after.is_deleted(built) => built.clone(),
            [image] => image.clone(),
            _ => ogeom_bail!(NotDone, "a face taking back its holes was split or dropped"),
        };
        let piece = &fused.pieces[kept[slot].0];
        let own = if piece.from_a { &fused.a } else { &fused.b };
        let source = &own.faces[piece.face];
        let side = usize::from(!piece.from_a);
        let wires: Vec<Shape> = holes
            .iter()
            .map(|&h| source.holes_aside[h].wire.clone())
            .collect();
        let joins: Vec<(usize, u32)> = holes
            .iter()
            .flat_map(|&h| source.holes_aside[h].beside.iter().map(move |&j| (side, j)))
            .collect();
        let walk = model.node(&image).map_or(0, |n| n.children().len());
        let made = with_holes(model, &image, &wires, tol)?;
        // A plane's holes lie inside its outer boundary, which alone sets
        // its box: the face with its holes back keeps the box of the face
        // without them.
        ogeom_algo::face_bounds(model, &image)?;
        let found = model.kept_face_box(&image, |_| {
            Err(ogeom_core::ogeom_err!(
                NotDone,
                "the face's box was not kept"
            ))
        });
        if let Ok(found) = found {
            model.kept_face_box(&made, |_| Ok(found))?;
        }
        let copied = after.copy_of(built).is_some() || after.modified(built).is_empty();
        if copied {
            put_back.copy(&image, made.clone());
        } else {
            put_back.modify(&image, made.clone());
        }
        remade.insert(image.node(), (made, walk, joins));
    }
    // Each face with its place: its solid, the place there of the face it
    // came from, and its piece's.
    let mut all = Vec::new();
    for shell in shells {
        for face in model.children_of(shell)? {
            let place = rank
                .get(&face.node())
                .copied()
                .unwrap_or((true, usize::MAX, usize::MAX));
            all.push(match remade.remove(&face.node()) {
                Some((made, walk, joins)) => (
                    place,
                    aside::Member {
                        face: made,
                        aside: None,
                        walk,
                        joins,
                    },
                ),
                None => (
                    place,
                    aside::Member {
                        face,
                        aside: None,
                        walk: usize::MAX,
                        joins: Vec::new(),
                    },
                ),
            });
        }
    }
    *history = core::mem::take(history).followed_by(&put_back);
    // The faces set aside are no image of anything the history holds so
    // far, so they are recorded as they stand.
    for (g, keep, from_b) in [
        (&fused.a, keep_aside[0], false),
        (&fused.b, keep_aside[1], true),
    ] {
        for (face, place) in &g.aside {
            if keep {
                history.copy(face, face.clone());
                all.push((
                    (from_b, *place, 0),
                    aside::Member {
                        face: face.clone(),
                        aside: Some((usize::from(from_b), *place)),
                        walk: 0,
                        joins: Vec::new(),
                    },
                ));
            } else {
                history.delete(face);
            }
        }
    }
    // In the order the solids list the faces they came from, as the
    // general fuse of every face lists them.
    all.sort_by_key(|(place, _)| *place);
    let all: Vec<aside::Member> = all.into_iter().map(|(_, member)| member).collect();
    if all.is_empty() {
        return Ok((Vec::new(), Vec::new(), Vec::new()));
    }
    let made = all
        .iter()
        .filter(|member| member.aside.is_none())
        .map(|member| member.face.clone())
        .collect();
    let (shells, closed) = aside::shells_of(model, &all, [&fused.a.beside, &fused.b.beside])?;
    Ok((shells, closed, made))
}

/// A face made again with its holes, how many of its wires it walks
/// itself, and the faces set aside beside its holes.
pub(crate) type Remade = (Shape, usize, Vec<(usize, u32)>);

/// What [`rebuilt_pieces`] makes of the kept pieces.
pub(crate) struct Rebuilt {
    /// One face per kept piece, in order.
    faces: Vec<Shape>,
    /// Whether each comes through the sew as it stands.
    settled: Vec<bool>,
    /// The distance the rebuild welds ends within.
    weld: f64,
    /// The holes left out of the arrangement each face takes back, by
    /// their place among its source face's.
    holes: Vec<Vec<usize>>,
    /// The kept pieces each face is built from, in order: those given, less
    /// any whose boundary closed on one point.
    kept: Vec<(usize, bool)>,
}

/// The holes left out of each face's arrangement, given to the kept piece
/// holding each.
///
/// A hole left out stands clear of the other solid, so the piece holding it
/// lies outside the other solid. Where the operation keeps what lies
/// outside, the face's one kept piece holds it, or of several the one whose
/// outline it lies within. Where it keeps what lies inside, no kept piece
/// holds one.
pub(crate) fn holes_put_back(
    model: &Model,
    fused: &GeneralFused,
    kept: &[(usize, bool)],
    tol: Tolerances,
) -> OgeomResult<Vec<Vec<usize>>> {
    let mut holes = vec![Vec::new(); kept.len()];
    let mut by_face: std::collections::BTreeMap<(bool, usize), Vec<usize>> =
        std::collections::BTreeMap::new();
    for (slot, &(index, _)) in kept.iter().enumerate() {
        let piece = &fused.pieces[index];
        let own = if piece.from_a { &fused.a } else { &fused.b };
        if !own.faces[piece.face].holes_aside.is_empty() {
            by_face
                .entry((piece.from_a, piece.face))
                .or_default()
                .push(slot);
        }
    }
    for ((from_a, fi), slots) in by_face {
        let face = if from_a {
            &fused.a.faces[fi]
        } else {
            &fused.b.faces[fi]
        };
        let outside = slots
            .iter()
            .any(|&slot| fused.pieces[kept[slot].0].state == PieceState::Out);
        if let [slot] = slots.as_slice()
            && outside
        {
            holes[*slot] = (0..face.holes_aside.len()).collect();
            continue;
        }
        for (h, hole) in face.holes_aside.iter().enumerate() {
            let Some(at) = hole_point(model, face, &hole.wire, tol)? else {
                ogeom_bail!(NotDone, "a hole set aside has no point in its face's chart");
            };
            let mut within = slots
                .iter()
                .filter(|&&slot| inside_rings(&fused.pieces[kept[slot].0].outlines, at));
            match (within.next(), within.next()) {
                (Some(&slot), None) if fused.pieces[kept[slot].0].state == PieceState::Out => {
                    holes[slot].push(h);
                }
                (None, _) if !outside => {}
                _ => ogeom_bail!(
                    NotDone,
                    "a hole set aside lies in no one piece standing outside the other solid"
                ),
            }
        }
    }
    Ok(holes)
}

/// A point of a hole's first edge in its face's chart.
pub(crate) fn hole_point(
    model: &Model,
    face: &GFace,
    hole: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Point2>> {
    let Some(surface) = face.own_surface else {
        return Ok(None);
    };
    let Some(edge) = model.children_of(hole)?.into_iter().next() else {
        return Ok(None);
    };
    let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
        return Ok(None);
    };
    let (pc, range) = match data.pcurve_for(surface, edge.location()) {
        Some(EdgeRepr::PCurve { curve, range, .. }) => (*curve, *range),
        Some(EdgeRepr::Seam { forward, range, .. }) => (*forward, *range),
        _ => return Ok(None),
    };
    let Some(pcurve) = model.geometry().pcurve(pc) else {
        return Ok(None);
    };
    Ok(Some(pcurve.point_at(f64::midpoint(range.0, range.1), tol)?))
}

/// `image` with the holes `holes` added: a new face on the same surface,
/// facing the same way.
pub(crate) fn with_holes(
    model: &mut Model,
    image: &Shape,
    holes: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Shape> {
    if !image.location().is_identity() {
        ogeom_bail!(
            NotDone,
            "a face taking back its holes stands under a placement"
        );
    }
    let Some(NodeData::Face(data)) = model.node(image).map(|n| n.data().clone()) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let mut data = *data;
    data.triangulation = None;
    // The holes' edges are the solid's own and are not widened: a face
    // looser than one of them takes none back. No edge is tighter than the
    // least tolerance.
    if data.tolerance.get() > tol.least_tolerance() {
        for hole in holes {
            for edge in model.children_of(hole)? {
                if model
                    .node(&edge)
                    .and_then(|n| n.data().as_edge())
                    .is_none_or(|e| e.tolerance.get() < data.tolerance.get())
                {
                    ogeom_bail!(NotDone, "a face is looser than a hole it takes back");
                }
            }
        }
    }
    let mut wires = model.children_of(&Shape::of(image.node()))?;
    wires.extend(holes.iter().cloned());
    Ok(model.add_face(data, &wires)?.oriented(image.orientation()))
}

/// Each shell without its membranes, and the history of their removal.
///
/// Two faces standing apart by less than the weld (two boxes stacked a few
/// microns apart) have their corners welded into one vertex each, and their
/// edges sewn into one: the two faces then bound one shell along the same
/// edges, lying on each other and facing apart, a wall of no thickness with
/// material on both sides. That is the opposed pair the boolean drops where
/// the contact is exact, and it is dropped here the same way. A pair bounded
/// by the same edges but standing apart (a lens's two caps on their one
/// rim) is a solid's own boundary and stays. A shell that is membranes
/// through (a sliver thinner than the weld) encloses nothing and goes.
/// The faces of a pair stand within [`MEMBRANE_WITHIN`] of each other too.
pub(crate) fn without_membranes(
    model: &mut Model,
    shells: &mut Vec<Shape>,
    weld: f64,
    tol: Tolerances,
) -> OgeomResult<History> {
    use ogeom_geom::Transformable as _;
    let mut dropped = History::new();
    let within = MEMBRANE_WITHIN.get();
    let bounding = |model: &Model, face: &Shape| -> OgeomResult<Vec<ogeom_topo::TShapeId>> {
        let mut edges: Vec<ogeom_topo::TShapeId> = explore_unique(model, face, ShapeType::Edge)?
            .iter()
            .filter(|e| {
                model
                    .node(e)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(|d| !d.degenerate)
            })
            .map(Shape::node)
            .collect();
        edges.sort_unstable();
        Ok(edges)
    };
    // Whether `second` lies on `first`, facing away from it.
    let membrane = |model: &Model, first: &Shape, second: &Shape| -> OgeomResult<bool> {
        let (at, n) = ogeom_algo::face_normal(model, first, tol)?;
        let (_, m) = ogeom_algo::face_normal(model, second, tol)?;
        if n.dot(m) > -0.99 * n.magnitude() * m.magnitude() {
            return Ok(false);
        }
        // The first's inner points, every one on the second's surface: a
        // thin lens (a tube's tip poking through a cap, cut off by it) is
        // two faces on the same loop facing apart too, as thin as a weld
        // only near its rim.
        let mesh =
            ogeom_mesh::triangulate_face(model, first, ogeom_mesh::Deflection::default(), tol)?;
        let mut probes: Vec<Point> = mesh
            .triangles
            .iter()
            .map(|t| {
                let [p, q, r] = t.map(|i| mesh.positions[i as usize]);
                Point::from_vector((p.to_vector() + q.to_vector() + r.to_vector()) / 3.0)
            })
            .collect();
        if probes.is_empty() {
            probes.push(at);
        }
        let Some(NodeData::Face(data)) = model.node(second).map(|n| n.data()) else {
            return Ok(false);
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            return Ok(false);
        };
        let placed = surface.transformed(&second.transform(model.datums())?, tol)?;
        // As far as the weld that joined their corners reached: its floor
        // past what the corners themselves own.
        let owned = explore_unique(model, first, ShapeType::Vertex)?
            .iter()
            .filter_map(|v| model.node(v).and_then(|n| n.data().as_vertex()))
            .fold(0.0_f64, |acc, d| acc.max(d.tolerance.get()));
        let reach = (weld + owned).min(within);
        for probe in probes {
            let foot = ogeom_algo::project_on_surface(&placed, probe, 16, tol)?;
            if foot.point.distance(probe) > reach {
                return Ok(false);
            }
        }
        Ok(true)
    };
    let mut emptied = Vec::new();
    for (index, shell) in shells.iter_mut().enumerate() {
        let faces = ogeom_topo::explore(model, shell, Filter::OfType(ShapeType::Face))?;
        let edges = faces
            .iter()
            .map(|f| bounding(model, f))
            .collect::<OgeomResult<Vec<_>>>()?;
        // Only faces bounded by the same edges can be a membrane's pair, so
        // the pairs are taken within each such group, in face order.
        let mut alike: ogeom_core::FastMap<&[ogeom_topo::TShapeId], Vec<usize>> =
            ogeom_core::FastMap::default();
        for (i, bounded_by) in edges.iter().enumerate() {
            if !bounded_by.is_empty() {
                alike.entry(bounded_by.as_slice()).or_default().push(i);
            }
        }
        let mut groups: Vec<Vec<usize>> = alike.into_values().filter(|g| g.len() > 1).collect();
        groups.sort_unstable();
        let mut gone = vec![false; faces.len()];
        for group in &groups {
            for (k, &i) in group.iter().enumerate() {
                for &j in &group[k + 1..] {
                    if gone[i] || gone[j] {
                        continue;
                    }
                    if membrane(model, &faces[i], &faces[j])? {
                        gone[i] = true;
                        gone[j] = true;
                    }
                }
            }
        }
        if !gone.contains(&true) {
            continue;
        }
        let mut kept = Vec::new();
        for (face, gone) in faces.iter().zip(&gone) {
            if *gone {
                dropped.delete(face);
            } else {
                kept.push(face.clone());
            }
        }
        if kept.is_empty() {
            emptied.push(index);
        } else {
            *shell = model.add_shell(&kept)?;
        }
    }
    for index in emptied.into_iter().rev() {
        shells.remove(index);
    }
    Ok(dropped)
}

/// Closed shells grouped into solids: each outer shell with the voids
/// directly inside it, as indices into `shells`.
///
/// A shell's depth is how many others enclose it. Even depth is a solid's
/// outer boundary (a free shell, or an island standing in another's
/// cavity), odd depth a void, belonging to the enclosing shell one level
/// out. Enclosure is decided by classifying the inner shell's vertices,
/// not by bounding boxes: a part standing in another's notch sits inside
/// its box and outside its material. A pair whose vertices all lie on the
/// other (shells touching everywhere they could be tested) falls back to
/// the boxes.
pub(crate) fn nest_shells(
    model: &Model,
    shells: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Vec<(usize, Vec<usize>)>> {
    if shells.len() == 1 {
        return Ok(vec![(0, Vec::new())]);
    }
    let mut bounds = Vec::with_capacity(shells.len());
    for shell in shells {
        bounds.push(shape_bounds(model, shell, tol)?);
    }
    // A container's box widened by the tolerance a touch is read within: a
    // void meeting its outer shell at a point (a ball resting against the
    // inside of a wall) may reach a hair past it, and is still within it.
    let reach = tol.confusion() * 1e2;
    let mut boundaries: Vec<Option<ogeom_algo::SolidBoundary>> =
        (0..shells.len()).map(|_| None).collect();
    let mut encloses = vec![vec![false; shells.len()]; shells.len()];
    for i in 0..shells.len() {
        for j in 0..shells.len() {
            if i == j || !bounds[i].expanded(reach).contains_box(&bounds[j]) {
                continue;
            }
            if boundaries[i].is_none() {
                boundaries[i] = Some(ogeom_algo::SolidBoundary::of(
                    model,
                    &shells[i],
                    tol.confusion() * 1e4,
                    tol,
                )?);
            }
            let Some(boundary) = boundaries[i].as_ref() else {
                continue;
            };
            let mut verdict = None;
            for vertex in explore_unique(model, &shells[j], ShapeType::Vertex)? {
                let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) else {
                    continue;
                };
                let at = vertex.transform(model.datums())?.apply(data.point);
                match boundary.holds(model, at, tol)? {
                    Containment::In => verdict = Some(true),
                    Containment::Out => verdict = Some(false),
                    Containment::On => continue,
                }
                break;
            }
            encloses[i][j] = verdict.unwrap_or(true);
        }
    }
    let depth: Vec<usize> = (0..shells.len())
        .map(|j| (0..shells.len()).filter(|&i| encloses[i][j]).count())
        .collect();
    let mut out: Vec<(usize, Vec<usize>)> = (0..shells.len())
        .filter(|&i| depth[i].is_multiple_of(2))
        .map(|i| (i, Vec::new()))
        .collect();
    for j in (0..shells.len()).filter(|&j| !depth[j].is_multiple_of(2)) {
        let parent = (0..shells.len()).find(|&i| encloses[i][j] && depth[i] + 1 == depth[j]);
        if let Some(parent) = parent
            && let Some(slot) = out.iter_mut().find(|(i, _)| *i == parent)
        {
            slot.1.push(j);
        }
    }
    Ok(out)
}
