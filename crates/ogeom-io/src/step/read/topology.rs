//! Topology: vertices, edges, faces, shells and solids, and the pcurves their bounds carry.

use super::*;

impl<'a> Reader<'a> {
    // --- topology ------------------------------------------------------------

    fn vertex(&mut self, id: u64) -> OgeomResult<Shape> {
        if let Some(shape) = self.vertices.get(&id) {
            return Ok(shape.clone());
        }
        let args = self.args(id, "VERTEX_POINT")?;
        let point = self.point(ref_at(args, 1))?;
        let shape = make_vertex(&mut self.model, point).shape;
        self.vertices.insert(id, shape.clone());
        Ok(shape)
    }

    /// An `EDGE_CURVE`, built once, curve-forward.
    ///
    /// Returns the edge, its curve, its range, and whether STEP's edge
    /// direction runs *against* the curve, which each use folds into its
    /// own orientation.
    fn edge(&mut self, id: u64) -> OgeomResult<Option<BuiltEdge>> {
        if let Some(found) = self.edges.get(&id) {
            return Ok(Some(found.clone()));
        }
        let args = self.args(id, "EDGE_CURVE")?;
        let v1 = ref_at(args, 1);
        let v2 = ref_at(args, 2);
        let Some(mut curve) = self.curve(ref_at(args, 3))? else {
            return Ok(None);
        };
        let same_sense = !args.get(4).is_some_and(|a| a.is_enum("F"));

        let p1 = {
            let vargs = self.args(v1, "VERTEX_POINT")?;
            self.point(ref_at(vargs, 1))?
        };
        let p2 = {
            let vargs = self.args(v2, "VERTEX_POINT")?;
            self.point(ref_at(vargs, 1))?
        };

        // The edge is built along the curve's own parameter; a STEP edge
        // running the other way is flagged, and every use composes the flag
        // into its orientation.
        let (start, end, mut flipped) = if same_sense {
            (p1, p2, false)
        } else {
            (p2, p1, true)
        };
        let (range, closed) = match (
            self.parameter_of(&curve, start),
            self.parameter_of(&curve, end),
        ) {
            (Some(a), Some(b)) => {
                let period = if curve.is_periodic() {
                    let (lo, hi) = curve.domain();
                    hi - lo
                } else {
                    0.0
                };
                if v1 == v2 {
                    ((a, a + if period > 0.0 { period } else { 0.0 }), true)
                } else if period > 0.0 && b <= a + self.tol.parametric() {
                    ((a, b + period), false)
                } else if b < a - self.tol.parametric() {
                    // The vertices stand at descending parameters on an
                    // open curve whatever the sense flag says: a
                    // mesh-to-STEP converter writes every edge forward and
                    // lets the line run the other way. The edge is built
                    // along the curve's own parameter and runs against it.
                    flipped = !flipped;
                    self.report.warnings.push(format!(
                        "#{id}: the edge's vertices run against its curve's \
                         parameter despite its sense flag; the edge was reversed"
                    ));
                    self.tally("edge-against-curve", a - b, id);
                    ((b, a), false)
                } else {
                    ((a, b), false)
                }
            }
            _ => {
                // No closed-form inversion: take the curve's own domain and
                // hold the endpoints to it.
                //
                // A closed edge whose one vertex sits on the curve but away
                // from its seam (a fitted loop written with its start
                // wherever the fit began, the vertex millimetres along it)
                // would otherwise be held to a seam the vertex misses by
                // that much, and the vertex's tolerance widened to say so.
                // The seam is moved to the vertex instead: the same curve,
                // begun where the edge does.
                if v1 == v2
                    && let Curve::BSpline(spline) = &curve
                {
                    let (lo, hi) = curve.domain();
                    let gap = curve.point_at(lo, self.tol)?.distance(start);
                    if gap > self.tol.confusion() * 10.0 {
                        let samples = (spline.control_points().len() * 8).max(64);
                        if let Ok(found) = project_on_curve(&curve, start, samples, self.tol)
                            && found.distance <= self.tol.confusion() * 1e4
                            && found.parameter > lo + self.tol.parametric()
                            && found.parameter < hi - self.tol.parametric()
                            && let Ok(moved) = spline.reseamed_at(found.parameter, self.tol)
                        {
                            curve = Curve::BSpline(moved);
                            self.report.warnings.push(format!(
                                "#{id}: a closed edge's vertex sits {gap:.2e} along its \
                                 curve from the seam; the seam was moved to the vertex"
                            ));
                            self.tally("reseamed", gap, id);
                        }
                    }
                }
                let (lo, hi) = curve.domain();
                let head = curve.point_at(lo, self.tol)?;
                let tail = curve.point_at(hi, self.tol)?;
                let (head_miss, tail_miss) = (head.distance(start), tail.distance(end));
                if head_miss > self.tol.confusion() * 10.0
                    || tail_miss > self.tol.confusion() * 10.0
                {
                    // An open edge whose vertices stand *on* the curve but
                    // short of its ends (a file that writes the whole
                    // spline and lets the vertices say where the edge
                    // stops, millimetres in) takes the window between the
                    // vertices' own feet. Held to the whole curve, the edge
                    // overshoots its neighbours and the face it bounds draws
                    // as nothing at all.
                    let window = if v1 == v2 {
                        None
                    } else {
                        window_between_feet(&curve, start, end, self.tol)
                    };
                    if let Some(window) = window {
                        self.report.warnings.push(format!(
                            "#{id}: edge endpoints sit {head_miss:.2e} and {tail_miss:.2e} \
                             from its curve's ends; the window between the vertices' \
                             feet on the curve was taken"
                        ));
                        self.tally("vertex-window", head_miss.max(tail_miss), id);
                        (window, false)
                    } else {
                        self.report.warnings.push(format!(
                            "#{id}: edge endpoints sit {head_miss:.2e} and {tail_miss:.2e} \
                             from its curve's ends; the curve's own domain was taken"
                        ));
                        ((lo, hi), v1 == v2)
                    }
                } else {
                    ((lo, hi), v1 == v2)
                }
            }
        };
        let _ = closed;

        let (vlo, vhi) = if flipped {
            (self.vertex(v2)?, self.vertex(v1)?)
        } else {
            (self.vertex(v1)?, self.vertex(v2)?)
        };
        // Real files are imprecise. Where the curve's end misses its vertex
        // by more than the default
        // tolerance, the vertex's tolerance *grows* to state the gap; the
        // data model's growing tolerances are exactly for this, and the
        // warning keeps the healing visible.
        for (vertex, t) in [(&vlo, range.0), (&vhi, range.1)] {
            let end = curve.point_at(t, self.tol)?;
            let stated = {
                let Some(node) = self.model.node(vertex) else {
                    continue;
                };
                let Some(data) = node.data().as_vertex() else {
                    continue;
                };
                (data.point, data.tolerance)
            };
            let gap = end.distance(stated.0);
            if gap > stated.1.get() {
                if let Some(node) = self.model.node_mut(vertex)
                    && let ogeom_topo::NodeData::Vertex(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen_to(gap + self.tol.confusion());
                }
                self.warn_vertex_miss(id, gap);
            }
        }
        let shape =
            make_edge_between(&mut self.model, curve.clone(), range, &vlo, &vhi, self.tol)?.shape;
        let entry = (shape, curve, range, flipped);
        self.edges.insert(id, entry.clone());
        Ok(Some(entry))
    }

    /// The edge a bound names as the bound walks it: its pieces where it
    /// was cut at a pole, in order along its curve, or the edge itself
    /// under its own id.
    fn edge_pieces(&mut self, id: u64) -> OgeomResult<Option<Vec<(u64, BuiltEdge)>>> {
        if let Some(pieces) = self.pieces.get(&id) {
            return Ok(Some(pieces.clone()));
        }
        Ok(self.edge(id)?.map(|built| vec![(id, built)]))
    }

    /// The keys a bound's use of edge `id` stands under.
    fn piece_keys(&self, id: u64) -> Vec<u64> {
        self.pieces.get(&id).map_or_else(
            || vec![id],
            |pieces| pieces.iter().map(|(k, _)| *k).collect(),
        )
    }

    /// What a file's edge id stands for in the read model: the pieces it
    /// was cut into at a pole, in order along its curve, or the edge
    /// itself; nothing when it never built. Everything that names an edge
    /// by id (a style, a shape aspect) resolves through here, so it lands
    /// on the edges the faces use.
    pub(super) fn edge_shapes(&self, id: u64) -> Vec<Shape> {
        if let Some(pieces) = self.pieces.get(&id) {
            return pieces
                .iter()
                .map(|(_, (piece, ..))| piece.clone())
                .collect();
        }
        self.edges
            .get(&id)
            .map(|(shape, ..)| shape.clone())
            .into_iter()
            .collect()
    }

    /// Cut every edge of these faces whose curve runs through a pole of a
    /// face it bounds (a sphere's pole, a cone's apex, a patch's collapsed
    /// side) at that pole, before any face is built.
    ///
    /// No single image of such a curve follows it through the chart: at the
    /// pole it leaves along another column than the one it arrived on, and
    /// a fit across the jump stands millimetres off the edge. Cut there,
    /// each piece has an image of its own (a meridian of a sphere is
    /// exactly a straight column), and the face's wire closes across the
    /// pole row as it does at any pole vertex. The cut is made on the edge,
    /// so every face it bounds walks the same pieces.
    ///
    /// Best-effort: anything unreadable here is left whole, and the face
    /// that uses it reports it as it would have.
    fn cut_at_poles(&mut self, face_ids: &[u64]) {
        let mut crossings: BTreeMap<u64, Vec<crate::pcurves::Crossing>> = BTreeMap::new();
        for &fid in face_ids {
            let Ok(args) = self.face_args(fid) else {
                continue;
            };
            let Some(surface) = args
                .get(2)
                .and_then(Arg::reference)
                .and_then(|sid| self.surface(sid).ok().flatten())
            else {
                continue;
            };
            let poles = crate::pcurves::chart_poles(&surface, self.tol);
            if poles.is_empty() {
                continue;
            }
            let bounds: Vec<u64> = args
                .get(1)
                .and_then(Arg::list)
                .unwrap_or(&[])
                .iter()
                .filter_map(Arg::reference)
                .collect();
            for bound in bounds {
                let Ok((loop_id, _)) = self.bound_args(bound) else {
                    continue;
                };
                let Ok(loop_args) = self.args(loop_id, "EDGE_LOOP") else {
                    continue;
                };
                let edge_ids: Vec<u64> = loop_args
                    .get(1)
                    .and_then(Arg::list)
                    .unwrap_or(&[])
                    .iter()
                    .filter_map(Arg::reference)
                    .filter_map(|oe| {
                        self.args(oe, "ORIENTED_EDGE")
                            .ok()
                            .and_then(|a| a.get(3).and_then(Arg::reference))
                    })
                    .collect();
                for edge_id in edge_ids {
                    if self.pieces.contains_key(&edge_id) {
                        continue;
                    }
                    let Ok(Some((_, curve, range, _))) = self.edge(edge_id) else {
                        continue;
                    };
                    let found = crate::pcurves::pole_crossings(&curve, range, &poles, self.tol);
                    if found.is_empty() {
                        continue;
                    }
                    let entry = crossings.entry(edge_id).or_default();
                    for crossing in found {
                        if entry.iter().all(|(t, ..)| (t - crossing.0).abs() > 1e-9) {
                            entry.push(crossing);
                        }
                    }
                }
            }
        }
        for (edge_id, mut at) in crossings {
            at.sort_by(|x, y| x.0.total_cmp(&y.0));
            let Some((shape, curve, range, flipped)) = self.edges.get(&edge_id).cloned() else {
                continue;
            };
            let Ok(cut) = crate::pcurves::split_at_poles(
                &mut self.model,
                &shape,
                &curve,
                range,
                &at,
                &mut self.pole_vertices,
                self.tol,
            ) else {
                continue;
            };
            let pieces = cut
                .into_iter()
                .map(|(piece, window)| {
                    self.piece_key -= 1;
                    (self.piece_key, (piece, curve.clone(), window, flipped))
                })
                .collect();
            self.pieces.insert(edge_id, pieces);
        }
    }

    /// Split one rim of every band whose two rims start on different
    /// columns, before any face is built.
    ///
    /// A face on a periodic surface bounded by two whole circles that are
    /// parallels of it needs a seam joining a vertex of one rim to a vertex
    /// of the other along one column. Where the rims' vertices stand on
    /// different columns, the second rim (or the first, where the second
    /// is already in pieces) is split on the other's column. The cut is
    /// made on the edge, so every face it bounds walks both pieces, and the
    /// band is then built through the common column.
    ///
    /// A split rim changes every face that holds it, so the faces are
    /// passed over again until a pass splits and marks nothing: a band on
    /// the far side of the rim, seen whole before the split, is marked as a
    /// band of rings then, or it would be built with a rim in pieces and no
    /// seam.
    ///
    /// Best-effort: a face whose rims resist is left as it was, and the
    /// face reports it as it would have.
    fn split_band_rims(&mut self, face_ids: &[u64]) {
        loop {
            let before = (self.pieces.len(), self.band_faces.len());
            for &fid in face_ids {
                let _ = self.split_band_rim(fid);
            }
            if (self.pieces.len(), self.band_faces.len()) == before {
                break;
            }
        }
    }

    /// [`Self::split_band_rims`] for one face.
    fn split_band_rim(&mut self, fid: u64) -> OgeomResult<()> {
        let args = self.face_args(fid)?.to_vec();
        let Some(surface) = args
            .get(2)
            .and_then(Arg::reference)
            .and_then(|sid| self.surface(sid).ok().flatten())
        else {
            return Ok(());
        };
        if !surface.is_periodic_u() {
            return Ok(());
        }
        let bounds: Vec<u64> = args
            .get(1)
            .and_then(Arg::list)
            .unwrap_or(&[])
            .iter()
            .filter_map(Arg::reference)
            .collect();
        let [first, second] = bounds[..] else {
            return Ok(());
        };
        // Each bound one oriented edge: the rim's file id.
        let mut rims = Vec::with_capacity(2);
        for bound in [first, second] {
            let (loop_id, _) = self.bound_args(bound)?;
            let Ok(loop_args) = self.args(loop_id, "EDGE_LOOP") else {
                return Ok(());
            };
            let edge_ids: Vec<u64> = loop_args
                .get(1)
                .and_then(Arg::list)
                .unwrap_or(&[])
                .iter()
                .filter_map(Arg::reference)
                .filter_map(|oe| {
                    self.args(oe, "ORIENTED_EDGE")
                        .ok()
                        .and_then(|a| a.get(3).and_then(Arg::reference))
                })
                .collect();
            let [edge_id] = edge_ids[..] else {
                return Ok(());
            };
            rims.push(edge_id);
        }
        // Each rim's edges as the faces will walk them, and its circle.
        let mut rings: Vec<(Vec<Shape>, ogeom_math::Circle)> = Vec::with_capacity(2);
        for &edge_id in &rims {
            let Some(pieces) = self.edge_pieces(edge_id)? else {
                return Ok(());
            };
            let Curve::Circle(c) = &pieces[0].1.1 else {
                return Ok(());
            };
            let circle = c.circle();
            let shapes: Vec<Shape> = pieces.into_iter().map(|(_, (s, ..))| s).collect();
            if shapes.len() == 1 {
                let Some((a, b)) = ogeom_algo::edge_vertices(&self.model, &shapes[0])? else {
                    return Ok(());
                };
                if !a.is_same(&b) {
                    return Ok(());
                }
            }
            rings.push((shapes, circle));
        }
        let all: Vec<&Shape> = rings.iter().flat_map(|(s, _)| s.iter()).collect();
        if !ogeom_algo::rings_are_parallels(&self.model, &surface, &all, self.tol)? {
            return Ok(());
        }
        let starts = |model: &Model, ring: &[Shape]| -> OgeomResult<Vec<(Shape, Point)>> {
            let mut out = Vec::with_capacity(ring.len());
            for edge in ring {
                if let Some((v, _)) = ogeom_algo::edge_vertices(model, edge)?
                    && let Some(data) = model.node(&v).and_then(|n| n.data().as_vertex())
                {
                    let at = data.point;
                    out.push((v, at));
                }
            }
            Ok(out)
        };
        let starts_a = starts(&self.model, &rings[0].0)?;
        let starts_b = starts(&self.model, &rings[1].0)?;
        let shared = starts_a.iter().any(|(_, p)| {
            starts_b
                .iter()
                .any(|(_, q)| same_column(&rings[1].1, *p, *q, self.tol))
        });
        let split_any = rings.iter().any(|(s, _)| s.len() > 1);
        if shared {
            if split_any {
                self.band_faces.insert(fid);
            }
            return Ok(());
        }
        // Split a rim still whole on the other's first vertex column.
        let (target, anchor) = if rings[1].0.len() == 1 {
            (1, &starts_a[0].0)
        } else if rings[0].0.len() == 1 {
            (0, &starts_b[0].0)
        } else {
            return Ok(());
        };
        let edge_id = rims[target];
        let Some((shape, curve, _, flipped)) = self.edges.get(&edge_id).cloned() else {
            return Ok(());
        };
        let Some(split) =
            ogeom_algo::split_ring_on_column(&mut self.model, &surface, &shape, anchor, self.tol)?
        else {
            return Ok(());
        };
        let mut pieces = Vec::with_capacity(2);
        for piece in split.history.modified(&shape) {
            let Some(ogeom_topo::EdgeRepr::Curve3d { range, .. }) = self
                .model
                .node(piece)
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.curve3d())
            else {
                return Ok(());
            };
            self.piece_key -= 1;
            pieces.push((
                self.piece_key,
                (piece.clone(), curve.clone(), *range, flipped),
            ));
        }
        self.pieces.insert(edge_id, pieces);
        self.band_faces.insert(fid);
        Ok(())
    }

    /// A band face whose rim was split on the other rim's column: the band
    /// of rings through a column both rims have a vertex on. Each ring's
    /// occurrences are turned to start there and keep the direction the
    /// file walks them, which says which side the face faces.
    fn band_of_split_rims(
        &mut self,
        surface: &SurfaceGeometry,
        wires: &[Shape],
    ) -> OgeomResult<Shape> {
        let [wire_a, wire_b] = wires else {
            ogeom_bail!(Construction, "a band has two rims");
        };
        let ring_a = self.model.ordered_children_of(wire_a)?;
        let ring_b = self.model.ordered_children_of(wire_b)?;
        let circle = {
            let Some(ogeom_topo::EdgeRepr::Curve3d { curve, .. }) = self
                .model
                .node(&ring_b[0])
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.curve3d())
            else {
                ogeom_bail!(Construction, "a band rim has no curve");
            };
            let Some(Curve::Circle(c)) = self.model.geometry().curve(*curve) else {
                ogeom_bail!(Construction, "a band rim is not a circle");
            };
            c.circle()
        };
        let start = |model: &Model, edge: &Shape| -> OgeomResult<Point> {
            let Some((v, _)) = ogeom_algo::edge_vertices(model, edge)? else {
                ogeom_bail!(Construction, "a band rim has no vertex");
            };
            let Some(data) = model.node(&v).and_then(|n| n.data().as_vertex()) else {
                ogeom_bail!(Construction, "a band rim's vertex holds no point");
            };
            Ok(data.point)
        };
        for i in 0..ring_a.len() {
            let p = start(&self.model, &ring_a[i])?;
            for j in 0..ring_b.len() {
                let q = start(&self.model, &ring_b[j])?;
                if !same_column(&circle, p, q, self.tol) {
                    continue;
                }
                let mut a = ring_a.clone();
                a.rotate_left(i);
                let mut b = ring_b.clone();
                b.rotate_left(j);
                return ogeom_algo::make_band_of_rings(&mut self.model, surface, &a, &b, self.tol);
            }
        }
        ogeom_bail!(Construction, "the band's rims share no vertex column")
    }

    fn face(&mut self, id: u64) -> OgeomResult<Option<Shape>> {
        if let Some(shape) = self.faces.get(&id) {
            return Ok(Some(shape.clone()));
        }
        let args = self.face_args(id)?;
        let bounds: Vec<u64> = args
            .get(1)
            .and_then(Arg::list)
            .unwrap_or(&[])
            .iter()
            .filter_map(Arg::reference)
            .collect();
        let Some(surface) = self.surface(ref_at(args, 2))? else {
            return Ok(None);
        };
        let face_forward = !args.get(3).is_some_and(|a| a.is_enum("F"));
        let surface_id = self.model.geometry_mut().add_surface(surface.clone());

        // Outer bound first, so the face's first wire is its outer ring.
        let mut ordered = bounds.clone();
        ordered.sort_by_key(|b| {
            self.instances
                .get(*b)
                .map_or(1, |i| i32::from(i.part("FACE_OUTER_BOUND").is_none()))
        });

        // Which edges this face uses twice: those are seams, and get both
        // sides' pcurves.
        let mut edge_uses: FastMap<u64, usize> = FastMap::default();
        for &bound in &ordered {
            let bargs = self.bound_args(bound)?;
            if self.instance(bargs.0)?.part("VERTEX_LOOP").is_some() {
                continue;
            }
            let loop_args = self.args(bargs.0, "EDGE_LOOP")?;
            for oe in loop_args.get(1).and_then(Arg::list).unwrap_or(&[]) {
                if let Some(oe_id) = oe.reference() {
                    let oargs = self.args(oe_id, "ORIENTED_EDGE")?;
                    if let Some(e) = oargs.get(3).and_then(Arg::reference) {
                        for key in self.piece_keys(e) {
                            *edge_uses.entry(key).or_default() += 1;
                        }
                    }
                }
            }
        }

        let mut wires = Vec::new();
        let mut annotated: FastSet<u64> = FastSet::default();
        for &bound in &ordered {
            let (loop_id, bound_forward) = self.bound_args(bound)?;
            if let Some(vertex_loop) = {
                let instance = self.instance(loop_id)?;
                instance.part("VERTEX_LOOP")
            } {
                // A loop of one vertex: a pole or an apex. It has no edges,
                // but it still bounds the face in parameter space, as a
                // degenerate edge running across the chart at the row the
                // point collapses to, exactly as native cones and spheres
                // are built.
                let vertex_id = vertex_loop.get(1).and_then(Arg::reference).unwrap_or(0);
                let vargs = self.args(vertex_id, "VERTEX_POINT")?;
                let at = self.point(ref_at(vargs, 1))?;
                let vertex = self.vertex(vertex_id)?;
                let projection = ogeom_algo::project_on_surface(&surface, at, 32, self.tol)?;
                let ((ua, ub), _) = surface.domain();
                let row = projection.parameters.1;
                let mut data = ogeom_topo::EdgeData::new();
                data.degenerate = true;
                let edge = self
                    .model
                    .add_edge(data, &[vertex.clone(), vertex.clone()])?;
                let pcurve: PlanarCurve = ogeom_geom::Line2d::segment(
                    ogeom_math::Point2::new(ua, row),
                    ogeom_math::Point2::new(ub, row),
                    self.tol,
                )?
                .into();
                ogeom_algo::attach_pcurve(
                    &mut self.model,
                    &edge,
                    pcurve,
                    surface_id,
                    Location::identity(),
                    (0.0, ub - ua),
                )?;
                wires.push(make_wire(&mut self.model, &[edge], self.tol)?.shape);
                continue;
            }
            let loop_args = self.args(loop_id, "EDGE_LOOP")?;
            let mut uses: Vec<BoundUse> = Vec::new();
            for oe in loop_args.get(1).and_then(Arg::list).unwrap_or(&[]) {
                let Some(oe_id) = oe.reference() else {
                    continue;
                };
                let oargs = self.args(oe_id, "ORIENTED_EDGE")?;
                let Some(edge_id) = oargs.get(3).and_then(Arg::reference) else {
                    continue;
                };
                let forward = !oargs.get(4).is_some_and(|a| a.is_enum("F"));
                let Some(pieces) = self.edge_pieces(edge_id)? else {
                    self.report.warnings.push(format!(
                        "#{id}: a bound references unreadable edge #{edge_id}; \
                         the face is skipped"
                    ));
                    return Ok(None);
                };
                let flipped = pieces.first().is_some_and(|(_, built)| built.3);
                // The use's direction composes the loop's, the bound's and
                // the edge-against-curve flag. The pieces of a cut edge go
                // in the order the loop as written walks them; reversing a
                // reversed bound below turns them with the rest.
                let mut use_forward = forward == bound_forward;
                if flipped {
                    use_forward = !use_forward;
                }
                let mut walked: Vec<BoundUse> = pieces
                    .into_iter()
                    .map(|(key, (shape, curve, range, _))| {
                        let placed = if use_forward {
                            shape.clone()
                        } else {
                            shape.reversed()
                        };
                        (placed, shape, key, curve, range)
                    })
                    .collect();
                if forward == flipped {
                    walked.reverse();
                }
                uses.extend(walked);
            }
            if !bound_forward {
                uses.reverse();
            }
            // The window may have grown to hold this bound, and the images
            // are derived on the surface as it stands.
            self.widen_window(surface_id, &uses)?;
            let widened = self
                .model
                .geometry()
                .surface(surface_id)
                .cloned()
                .unwrap_or_else(|| surface.clone());
            self.chart_bound(id, &uses, &widened, surface_id, &edge_uses, &mut annotated)?;
            let edges: Vec<Shape> = uses.into_iter().map(|(placed, ..)| placed).collect();
            if edges.is_empty() {
                continue;
            }
            wires.push(make_wire(&mut self.model, &edges, self.tol)?.shape);
        }
        if wires.is_empty() {
            self.report
                .warnings
                .push(format!("#{id}: a face with no readable bounds is skipped"));
            return Ok(None);
        }
        // As the bounds left it: a window stretched to hold them is the
        // surface every path below builds on, including the band, which
        // registers the value it is handed rather than the one on record.
        let surface = self
            .model
            .geometry()
            .surface(surface_id)
            .cloned()
            .unwrap_or(surface);

        // A band whose rim was split on the other rim's column takes its
        // seam there. The rings are walked as the file walks them about the
        // face's own normal, which already says which side the face faces.
        if wires.len() == 2 && self.band_faces.contains(&id) {
            match self.band_of_split_rims(&surface, &wires) {
                Ok(shape) => {
                    self.faces.insert(id, shape.clone());
                    return Ok(Some(shape));
                }
                Err(e) => {
                    self.report.warnings.push(format!(
                        "#{id}: no seam could be synthesised ({e}); the \
                         face may not triangulate"
                    ));
                }
            }
        }

        // A periodic face bound only by closed rings (a cylinder band
        // between two circles) arrives without a seam edge, which is a
        // legitimate STEP shape and an open rectangle in this kernel's
        // chart. The seam is synthesised the way native cylinders build it:
        // one edge at the period join, appearing in the wire twice.
        if wires.len() == 2
            && surface.is_periodic_u()
            && let [(e_lo, _v_lo), (e_hi, _v_hi)] =
                closed_ring_edges(&self.model, &wires)?.as_slice()
            // Only rings that are parallels of the surface make a band.
            // Two closed circles that merely lie on it (a button head's
            // rims, square to the screw on a sphere whose chart runs along
            // z) bound a legitimate face on their own, nested loops in the
            // chart, and take the ordinary path below without a word.
            && ogeom_algo::rings_are_parallels(&self.model, &surface, &[e_lo, e_hi], self.tol)?
        {
            {
                // The band construction is ogeom-algo's make_revolution_band:
                // one authority shared with the healer. Anything it refuses
                // (a spline ring, ring vertices at different angles, a
                // surface with no closed-form iso-curve) becomes a warning
                // and the raw bounds, not an error.
                match ogeom_algo::make_revolution_band(
                    &mut self.model,
                    &surface,
                    e_lo,
                    e_hi,
                    self.tol,
                ) {
                    Ok(built) => {
                        let shape = if face_forward {
                            built
                        } else {
                            built.reversed()
                        };
                        self.faces.insert(id, shape.clone());
                        return Ok(Some(shape));
                    }
                    Err(e) => {
                        self.report.warnings.push(format!(
                            "#{id}: no seam could be synthesised ({e}); the \
                             face may not triangulate"
                        ));
                    }
                }
            }
        }

        // A cone face bounded by a *single* closed ring: the apex is the
        // other boundary, and some files simply never write it: no vertex
        // loop, nothing. The geometry leaves one choice of what the face
        // means, so the apex is synthesised and the band built as if the
        // file had said so.
        if wires.len() == 1
            && matches!(surface, SurfaceGeometry::Cone(_))
            && let [(ring, _)] = closed_ring_edges(&self.model, &wires)?.as_slice()
        {
            match ogeom_algo::make_apex_band(&mut self.model, &surface, ring, self.tol) {
                Ok(built) => {
                    let shape = if face_forward {
                        built
                    } else {
                        built.reversed()
                    };
                    self.faces.insert(id, shape.clone());
                    return Ok(Some(shape));
                }
                Err(e) => {
                    self.report.warnings.push(format!(
                        "#{id}: no apex could be synthesised ({e}); the face \
                         may not triangulate"
                    ));
                }
            }
        }

        // A sphere face bounded by a single parallel: the pole on the side
        // the face's material lies is the other boundary, which files leave
        // unwritten as they leave a cone's apex. The ring is walked with the
        // material on its left about the face's own normal, so the pole
        // that left points toward is the one the face reaches.
        if wires.len() == 1
            && let SurfaceGeometry::Sphere(ball) = &surface
            && let [(ring, at)] = closed_ring_edges(&self.model, &wires)?.as_slice()
            && ogeom_algo::rings_are_parallels(&self.model, &surface, &[ring], self.tol)?
            && let Some(pole) = pole_on_the_left(
                &self.model,
                &surface,
                ball.sphere(),
                ring,
                at,
                face_forward,
                self.tol,
            )?
        {
            let vertex = self.model.add_vertex(ogeom_topo::VertexData::new(pole));
            let mut data = ogeom_topo::EdgeData::new();
            data.degenerate = true;
            let collapsed = self.model.add_edge(data, &[vertex.clone(), vertex])?;
            match ogeom_algo::make_revolution_band(
                &mut self.model,
                &surface,
                ring,
                &collapsed,
                self.tol,
            ) {
                Ok(built) => {
                    let shape = if face_forward {
                        built
                    } else {
                        built.reversed()
                    };
                    self.faces.insert(id, shape.clone());
                    return Ok(Some(shape));
                }
                Err(e) => {
                    self.report.warnings.push(format!(
                        "#{id}: no pole could be synthesised ({e}); the face \
                         may not triangulate"
                    ));
                }
            }
        }

        // A spline closed round with no seam written: the rings close into
        // one loop up and down the period's join.
        if let Some(belt) = ogeom_algo::belted_wire(&mut self.model, surface_id, &wires, self.tol)?
        {
            wires = vec![belt];
        }
        // Every wire's images on one branch of a periodic chart: a file's
        // pcurves, or fits of them, answer in whatever phase they like, and
        // a hole loop straddling a drum's seam arrives half on each branch.
        ogeom_algo::chain_wire_branches(&mut self.model, surface_id, &wires, self.tol)?;
        // With every image on one branch, a boundary meeting a pole from two
        // sides of the chart is closed across the pole row.
        for wire in &mut wires {
            let edges = self.model.ordered_children_of(wire)?;
            let closed =
                ogeom_algo::closed_at_poles(&mut self.model, surface_id, &edges, self.tol)?;
            if closed.len() != edges.len() {
                *wire = make_wire(&mut self.model, &closed, self.tol)?.shape;
            }
        }
        // A file walks each loop with the face's material on its left about
        // the face's own normal; a face turned against its surface is stored
        // with its loops walked back, so that turning it walks them as
        // written.
        if !face_forward {
            wires = walked_back(&mut self.model, &wires, self.tol)?;
        }
        let mut built = make_face_on(&mut self.model, surface_id, &wires, self.tol)?.shape;
        // A pole crossed a whole turn round, beside a joint a whole period
        // apart where `closed_at_poles` lays no run: the crossing solved
        // for the face as built, and its wires rebuilt with the runs in.
        let crossed = ogeom_algo::pole_crossings(&mut self.model, &built, self.tol)?;
        if !crossed.is_empty() {
            for wire in &mut wires {
                if let Some((_, new)) = crossed.iter().find(|(old, _)| old.node() == wire.node()) {
                    *wire = new.clone();
                }
            }
            built = make_face_on(&mut self.model, surface_id, &wires, self.tol)?.shape;
            self.report.warnings.push(format!(
                "#{id}: the boundary crosses a pole of its surface a whole turn \
                 round; a degenerate edge was laid across it"
            ));
            self.tally("pole-crossed", 0.0, id);
        }
        let shape = if face_forward {
            built
        } else {
            built.reversed()
        };
        self.faces.insert(id, shape.clone());
        Ok(Some(shape))
    }

    /// A bound's loop and orientation, whichever of the two bound kinds it is.
    fn bound_args(&mut self, id: u64) -> OgeomResult<(u64, bool)> {
        let instance = self.instance(id)?;
        let args = instance
            .part("FACE_OUTER_BOUND")
            .or_else(|| instance.part("FACE_BOUND"))
            .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "#{id} is not a face bound"))?
            .to_vec();
        let loop_id = args.get(1).and_then(Arg::reference).unwrap_or(0);
        let forward = !args.get(2).is_some_and(|a| a.is_enum("F"));
        Ok((loop_id, forward))
    }

    /// Stretch a surface's parameter window to hold the edges that bound it.
    ///
    /// The window a reader gives a plane, a cylinder or a cone is a
    /// convention: those surfaces are unbounded, and [`SURFACE_EXTENT`] is
    /// a guess at how far past its own geometry a file will reach. A file
    /// can falsify the guess: a surface whose own origin sits far from the
    /// face it carries has trim parameters past any fixed window, the
    /// surface refuses to be evaluated where its own face lies, and the
    /// face draws as a hole.
    ///
    /// So the window is measured rather than guessed: the edges' own
    /// points, projected onto the surface's chart by the geometry that
    /// defines it, with the guess kept as a floor. Widening is safe for a
    /// surface other faces share, since a window only ever grows.
    fn widen_window(
        &mut self,
        surface_id: ogeom_topo::SurfaceId,
        uses: &[BoundUse],
    ) -> OgeomResult<()> {
        use ogeom_geom::Curve3d as _;
        let Some(surface) = self.model.geometry().surface(surface_id).cloned() else {
            return Ok(());
        };
        // Where each surface keeps its unbounded directions, and how a point
        // of space reads on them.
        let along: Box<ChartReading> = match &surface {
            SurfaceGeometry::Plane(p) => {
                let frame = p.plane().frame();
                let (origin, x, y) = (frame.origin(), frame.x().vector(), frame.y().vector());
                Box::new(move |at: Point| (Some((at - origin).dot(x)), Some((at - origin).dot(y))))
            }
            SurfaceGeometry::Cylinder(c) => {
                let axis = c.cylinder().axis();
                let (origin, direction) = (axis.location, axis.direction.vector());
                Box::new(move |at: Point| (None, Some((at - origin).dot(direction))))
            }
            SurfaceGeometry::Cone(c) => {
                let axis = c.cone().axis();
                let (origin, direction) = (axis.location, axis.direction.vector());
                Box::new(move |at: Point| (None, Some((at - origin).dot(direction))))
            }
            // A swept curve: a point's sweep parameter is its reach along
            // the direction less the swept curve's own, read against the
            // middle of the curve's reach; the curve's half-reach either
            // way is inside the margin the window grows by.
            SurfaceGeometry::Extrusion(e) => {
                let direction = e.direction().vector();
                let (lo, hi) = e.curve().domain();
                let (mut least, mut most) = (f64::INFINITY, f64::NEG_INFINITY);
                for k in 0..=32 {
                    let t = lo + (hi - lo) * f64::from(k) / 32.0;
                    if let Ok(p) = e.curve().point_at(t, self.tol) {
                        let reach = p.to_vector().dot(direction);
                        least = least.min(reach);
                        most = most.max(reach);
                    }
                }
                let middle = f64::midpoint(least, most);
                Box::new(move |at: Point| (None, Some(at.to_vector().dot(direction) - middle)))
            }
            _ => return Ok(()),
        };
        const STATIONS: usize = 4;
        let (mut ua, mut ub) = (-SURFACE_EXTENT, SURFACE_EXTENT);
        let (mut va, mut vb) = (-SURFACE_EXTENT, SURFACE_EXTENT);
        for (_, _, _, curve, range) in uses {
            for step in 0..=STATIONS {
                #[allow(clippy::cast_precision_loss)]
                let t = range.0 + (range.1 - range.0) * (step as f64 / STATIONS as f64);
                let Ok(at) = curve.point_at(t, self.tol) else {
                    continue;
                };
                let (u, v) = along(at);
                if let Some(u) = u {
                    ua = ua.min(u);
                    ub = ub.max(u);
                }
                if let Some(v) = v {
                    va = va.min(v);
                    vb = vb.max(v);
                }
            }
        }
        let ((was_u, was_v), grown) = (surface.domain(), (ub - ua) * 0.05);
        if ua >= was_u.0 && ub <= was_u.1 && va >= was_v.0 && vb <= was_v.1 {
            return Ok(());
        }
        let (u, v) = (
            (ua.min(was_u.0) - grown, ub.max(was_u.1) + grown),
            (va.min(was_v.0) - grown, vb.max(was_v.1) + grown),
        );
        let wider = match surface {
            SurfaceGeometry::Plane(p) => PlaneSurface::over(p.plane(), u, v)?.into(),
            SurfaceGeometry::Cylinder(c) => CylinderSurface::new(c.cylinder(), v)?.into(),
            SurfaceGeometry::Cone(c) => ConeSurface::new(c.cone(), v)?.into(),
            SurfaceGeometry::Extrusion(e) => {
                ExtrusionSurface::over(e.curve().clone(), e.direction(), v)?.into()
            }
            other => other,
        };
        if let Some(held) = self.model.geometry_mut().surface_mut(surface_id) {
            *held = wider;
        }
        Ok(())
    }

    /// One bound's images, chained around the face's chart.
    ///
    /// Each is derived on its own (the exact projection where the pair has
    /// a closed form, the fitted one where it does not), and a periodic
    /// chart then leaves a branch to choose. Chosen edge by edge the choice
    /// is arbitrary and the wire comes apart: a drilled block's bore wall
    /// reads back with one rim two whole turns from the other, both the
    /// right circle and neither meeting the seam the wire closes on. So
    /// each image after the first is slid by whole periods until its start
    /// meets where the last one ended.
    ///
    /// A seam falls out of the same walk. The wire uses it twice, up one
    /// column of the chart and down the other, and those columns are one
    /// image a period apart, so the walk finds both, and the use that runs
    /// forward is the forward side. Where it cannot, because the chart
    /// closes without being periodic, the other column goes a chart's width
    /// over, toward the middle.
    fn chart_bound(
        &mut self,
        face_id: u64,
        uses: &[BoundUse],
        surface: &SurfaceGeometry,
        surface_id: ogeom_topo::SurfaceId,
        edge_uses: &FastMap<u64, usize>,
        annotated: &mut FastSet<u64>,
    ) -> OgeomResult<()> {
        use ogeom_geom::Curve2d as _;
        let mut images: FastMap<u64, PlanarCurve> = FastMap::default();
        let mut sides: FastMap<u64, (Option<PlanarCurve>, Option<PlanarCurve>)> =
            FastMap::default();
        let mut order: Vec<u64> = Vec::new();
        let mut previous: Option<ogeom_math::Point2> = None;
        let mut opening: Option<(u64, bool, f64)> = None;
        for (placed, shape, edge_id, curve, range) in uses {
            if !images.contains_key(edge_id) {
                let Some(image) =
                    self.image_for(face_id, *edge_id, shape, curve, *range, surface)?
                else {
                    continue;
                };
                images.insert(*edge_id, image);
                order.push(*edge_id);
            }
            let Some(image) = images.get(edge_id).cloned() else {
                continue;
            };
            let backwards = placed.orientation() == ogeom_topo::Orientation::Reversed;
            let (start, end) = if backwards {
                (range.1, range.0)
            } else {
                (range.0, range.1)
            };
            if previous.is_none() && edge_uses.get(edge_id).copied().unwrap_or(1) == 1 {
                opening = Some((*edge_id, backwards, start));
            }
            let mut image =
                crate::pcurves::shifted_to_meet(&image, start, previous, surface, self.tol)?;
            let walked = sides.entry(*edge_id).or_default();
            let first = if backwards { &walked.0 } else { &walked.1 };
            if edge_uses.get(edge_id).copied().unwrap_or(1) > 1
                && let Some(first) = first
            {
                image = crate::pcurves::across_from(image, first, *range, surface, self.tol)?;
            }
            previous = Some(image.point_at(end, self.tol)?);
            if backwards {
                walked.1 = Some(image);
            } else {
                walked.0 = Some(image);
            }
        }
        // The first image, which no edge before it placed, meets the last
        // across a closing line it lies on.
        if let (Some((edge_id, backwards, start)), Some(end)) = (opening, previous)
            && let Some(walked) = sides.get_mut(&edge_id)
        {
            let side = if backwards {
                &mut walked.1
            } else {
                &mut walked.0
            };
            if let Some(image) = side.as_ref() {
                *side = Some(crate::pcurves::onto_closing_side(
                    image, start, end, surface, self.tol,
                )?);
            }
        }
        for edge_id in order {
            if !annotated.insert(edge_id) {
                continue;
            }
            let Some((shape, range)) = uses
                .iter()
                .find(|(_, _, id, _, _)| *id == edge_id)
                .map(|(_, shape, _, _, range)| (shape.clone(), *range))
            else {
                continue;
            };
            let Some(image) = images.get(&edge_id).cloned() else {
                continue;
            };
            let walked = sides.remove(&edge_id).unwrap_or((None, None));
            let seam = edge_uses.get(&edge_id).copied().unwrap_or(1) > 1;
            let columns = match &walked {
                (Some(forward), Some(reversed)) => {
                    let at = forward.point_at(range.0, self.tol)?;
                    let other = reversed.point_at(range.0, self.tol)?;
                    (at.distance(other) > self.tol.confusion())
                        .then(|| (forward.clone(), reversed.clone()))
                }
                _ => None,
            };
            match (seam, columns) {
                (true, Some((forward, reversed))) => {
                    self.record_columns(&shape, forward, reversed, surface_id, range)?;
                }
                (true, None) => {
                    self.record_pcurve(&shape, image, surface_id, true, range)?;
                }
                (false, _) => {
                    let (forward, reversed) = walked;
                    let placed = forward.or(reversed).unwrap_or(image);
                    self.record_pcurve(&shape, placed, surface_id, false, range)?;
                }
            }
        }
        Ok(())
    }

    fn image_for(
        &mut self,
        face_id: u64,
        edge_id: u64,
        edge: &Shape,
        curve: &Curve,
        range: (f64, f64),
        surface: &SurfaceGeometry,
    ) -> OgeomResult<Option<PlanarCurve>> {
        // Taken from the table where the pass at the head of the solid
        // already derived it. The fallbacks below serve the faces no pass
        // covered: a face reached outside a solid walk, or one whose
        // preparation refused.
        let (prepared, offered) = match self.pcurves.remove(&(face_id, edge_id)) {
            Some(prepared) => prepared,
            None => {
                let stated = self.stated_pcurves(face_id, edge_id);
                let offered = !stated.is_empty();
                (
                    derive_pcurve(curve, range, surface, &stated, self.tol),
                    offered,
                )
            }
        };
        if matches!(prepared, PreparedPcurve::Stated { .. }) {
            self.report.pcurves_kept += 1;
        } else if offered {
            self.report.pcurves_derived += 1;
        }
        let (pcurve, off) = match prepared {
            PreparedPcurve::Exact { curve, off } | PreparedPcurve::Stated { curve, off } => {
                (crate::pcurves::over_range(&curve, range), off)
            }
            PreparedPcurve::Fitted {
                curve: fitted,
                error,
                met,
                worst_off,
                off,
                warning,
            } => {
                if let Some(w) = warning {
                    self.warn_slop(w, worst_off, face_id);
                }
                if !met {
                    self.warn_fit_short(face_id, error);
                }
                (fitted, off.max(worst_off))
            }
            PreparedPcurve::Refused(why) => {
                self.report.warnings.push(format!(
                    "face #{face_id}: no pcurve for an edge on this surface ({why}); \
                     the face may not triangulate"
                ));
                self.note_untrimmed(face_id);
                return Ok(None);
            }
        };
        // The edge provably stands `off` from the face it bounds, whether
        // the curve lies off the surface or a fit strays between its
        // samples; its tolerance grows to cover that, and its vertices
        // follow once the solid is whole. Past the millimetre a file's
        // slop stays under, the pcurve is a fit gone astray (a curve
        // through a pole left whole, whose image leaps across the chart
        // there), and
        // an edge that wide would swallow its neighbours; the face is
        // reported instead.
        if off > self.tol.confusion() * 1e7 {
            self.warn_fit_short(face_id, off);
        } else if off > self.tol.confusion()
            && let Some(node) = self.model.node_mut(edge)
            && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
        {
            data.tolerance = data.tolerance.widen_to(off + self.tol.confusion());
        }
        Ok(Some(pcurve))
    }

    /// The pcurves the file states for an edge on a face's surface. A key
    /// standing for a piece of an edge cut at a pole names no edge in the
    /// file, and the file's pcurves run over the whole edge in any case.
    fn stated_pcurves(&mut self, face_id: u64, edge_id: u64) -> Vec<PlanarCurve> {
        if self.pieces.contains_key(&edge_id) || !self.instances.contains(edge_id) {
            return Vec::new();
        }
        let Some(surface_id) = self
            .face_args(face_id)
            .ok()
            .and_then(|args| args.get(2).and_then(Arg::reference))
        else {
            return Vec::new();
        };
        self.file_pcurves(edge_id, surface_id)
    }

    /// Count one occurrence of a warning kind toward the summary.
    fn tally(&mut self, kind: &'static str, measured: f64, exemplar: u64) {
        let entry = self.tallies.entry(kind).or_insert((0, 0.0, exemplar));
        entry.0 += 1;
        if measured > entry.1 {
            entry.1 = measured;
            entry.2 = exemplar;
        }
    }

    /// A curve end missing its vertex: the prose, and the count.
    fn warn_vertex_miss(&mut self, id: u64, gap: f64) {
        self.report.warnings.push(format!(
            "#{id}: a curve end misses its vertex by {gap:.2e}; the vertex \
             tolerance grew to say so"
        ));
        self.tally("vertex-miss", gap, id);
    }

    /// A fitted trim that stopped short of its target: prose and count.
    fn warn_fit_short(&mut self, face_id: u64, error: f64) {
        self.report.warnings.push(format!(
            "face #{face_id}: a projected pcurve fit stopped at {error:.2e}; \
             the face's mesh may sit that far off along this edge"
        ));
        self.tally("fit-short", error, face_id);
    }

    /// The file's own boundary slop, carried and counted.
    fn warn_slop(&mut self, prose: String, worst: f64, exemplar: u64) {
        self.report.warnings.push(prose);
        self.tally("boundary-slop", worst, exemplar);
    }

    /// Note a face left without a complete trim, once.
    fn note_untrimmed(&mut self, face_id: u64) {
        if self.untrimmed_ids.last() != Some(&face_id) {
            self.untrimmed_ids.push(face_id);
            self.tally("untrimmed", 0.0, face_id);
        }
    }

    /// Attach a seam whose two columns the wire's own walk already found.
    fn record_columns(
        &mut self,
        edge: &Shape,
        forward: PlanarCurve,
        reversed: PlanarCurve,
        surface_id: ogeom_topo::SurfaceId,
        range: (f64, f64),
    ) -> OgeomResult<()> {
        ogeom_algo::attach_seam(
            &mut self.model,
            edge,
            forward,
            reversed,
            surface_id,
            Location::identity(),
            range,
        )
    }

    /// Attach a derived pcurve, seaming it where the edge bounds the chart
    /// twice and the walk could not say where its other column stands.
    fn record_pcurve(
        &mut self,
        edge: &Shape,
        pcurve: PlanarCurve,
        surface_id: ogeom_topo::SurfaceId,
        seam: bool,
        range: (f64, f64),
    ) -> OgeomResult<()> {
        if seam {
            let Some(surface) = self.model.geometry().surface(surface_id).cloned() else {
                ogeom_bail!(Dangling, "the surface is not in this model");
            };
            // One side is where the projection landed; the other is one
            // period over, along the axis the surface actually closes on.
            let other = crate::pcurves::seam_other_side(&pcurve, range, &surface, self.tol)?;
            ogeom_algo::attach_seam(
                &mut self.model,
                edge,
                pcurve,
                other,
                surface_id,
                Location::identity(),
                range,
            )?;
        } else {
            ogeom_algo::attach_pcurve(
                &mut self.model,
                edge,
                pcurve,
                surface_id,
                Location::identity(),
                range,
            )?;
        }
        Ok(())
    }

    /// Derive every pcurve this solid's faces will want, in parallel.
    ///
    /// Deriving one is a pure function of a curve, its range and a surface:
    /// it reads no model and depends on no order, and it is most of the time
    /// spent building solids. The walk that follows attaches them, in file
    /// order.
    ///
    /// Resolving *what* to derive still runs on the walk's own thread, since
    /// it builds edges and vertices into the model. That part is cheap; it is
    /// the projection and the fitting that are not.
    ///
    /// Best-effort by design: anything this cannot resolve is simply left out
    /// of the table, and `chart_bound` derives it on the walk. So a face
    /// shape this does not anticipate costs time, never correctness.
    fn prepare_pcurves(&mut self, face_ids: &[u64]) {
        struct Job {
            face: u64,
            edge: u64,
            curve: Curve,
            range: (f64, f64),
            /// The pcurves the file states for it, tried first.
            stated: Vec<PlanarCurve>,
            /// Index into `surfaces`. The surface is held once per face, not
            /// once per edge: a B-spline patch owns its whole control grid,
            /// and cloning that per edge costs more than the projection it
            /// was cloned for.
            surface: usize,
        }
        let mut surfaces: Vec<SurfaceGeometry> = Vec::new();
        let mut jobs: Vec<Job> = Vec::new();
        let mut seen: FastSet<(u64, u64)> = FastSet::default();
        for &fid in face_ids {
            let Ok(args) = self.face_args(fid) else {
                continue;
            };
            let Some(surface) = args
                .get(2)
                .and_then(Arg::reference)
                .and_then(|sid| self.surface(sid).ok().flatten())
            else {
                continue;
            };
            surfaces.push(surface);
            let at = surfaces.len() - 1;
            let bounds: Vec<u64> = args
                .get(1)
                .and_then(Arg::list)
                .unwrap_or(&[])
                .iter()
                .filter_map(Arg::reference)
                .collect();
            for bound in bounds {
                let Ok((loop_id, _)) = self.bound_args(bound) else {
                    continue;
                };
                let Ok(loop_args) = self.args(loop_id, "EDGE_LOOP") else {
                    continue;
                };
                let uses: Vec<u64> = loop_args
                    .get(1)
                    .and_then(Arg::list)
                    .unwrap_or(&[])
                    .iter()
                    .filter_map(Arg::reference)
                    .collect();
                for oe_id in uses {
                    let Ok(oargs) = self.args(oe_id, "ORIENTED_EDGE") else {
                        continue;
                    };
                    let Some(edge_id) = oargs.get(3).and_then(Arg::reference) else {
                        continue;
                    };
                    if !seen.insert((fid, edge_id)) {
                        continue;
                    }
                    if let Ok(Some(pieces)) = self.edge_pieces(edge_id) {
                        for (key, (_, curve, range, _)) in pieces {
                            let stated = self.stated_pcurves(fid, key);
                            jobs.push(Job {
                                face: fid,
                                edge: key,
                                curve,
                                range,
                                stated,
                                surface: at,
                            });
                        }
                    }
                }
            }
        }
        // Below a handful of edges the threads cost more than the work; the
        // sequential path through `chart_bound` is already correct, so the
        // table is simply left empty and the walk derives them itself.
        if jobs.len() < 16 {
            return;
        }
        let tol = self.tol;
        let derived = ogeom_core::parallel::map_ordered(&jobs, |_, job| {
            derive_pcurve(
                &job.curve,
                job.range,
                &surfaces[job.surface],
                &job.stated,
                tol,
            )
        });
        for (job, pcurve) in jobs.iter().zip(derived) {
            self.pcurves
                .insert((job.face, job.edge), (pcurve, !job.stated.is_empty()));
        }
    }

    /// The solid a `MANIFOLD_SOLID_BREP` names: its shell's faces, sewn.
    ///
    /// `BREP_WITH_VOIDS` is the same entity with cavities: a subtype of
    /// `MANIFOLD_SOLID_BREP`, so its first two attributes are the name and
    /// the outer shell, and a third names the shells that bound the voids.
    /// A reader matching on the leading keyword alone does not see it, and
    /// the part vanishes from the import. The voids join the solid as
    /// shells of their own, every normal pointing away from the material:
    /// out of the body on the outside, into the cavity within. A void is
    /// oriented by its geometry where its flags
    /// would have it add material.
    pub(super) fn solid(&mut self, id: u64) -> OgeomResult<Shape> {
        let instance = self.instance(id)?;
        let args = instance
            .part("MANIFOLD_SOLID_BREP")
            .or_else(|| instance.part("BREP_WITH_VOIDS"))
            .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "#{id} is not a solid"))?
            .to_vec();
        let shell_id = args.get(1).and_then(Arg::reference).unwrap_or(0);
        let Some(shell) = self.shell(shell_id)? else {
            ogeom_bail!(Construction, "#{id}: a solid with no readable faces");
        };
        if !ogeom_algo::is_shell_closed(&self.model, &shell)? {
            self.report.warnings.push(format!(
                "#{id}: the shell does not close as read; measures needing an \
                 inside will refuse it"
            ));
        }
        let mut shells = vec![shell];
        for void_id in args
            .get(2)
            .and_then(Arg::list)
            .unwrap_or(&[])
            .iter()
            .filter_map(Arg::reference)
            .collect::<Vec<u64>>()
        {
            match self.shell(void_id) {
                Ok(Some(void)) => {
                    // A void is a cavity whatever its flags say: its faces
                    // point into it, so it encloses negative volume. One a
                    // writer turned twice reads as added material, and is
                    // turned back as it is built.
                    let encloses = ogeom_mesh::triangulate(
                        &self.model,
                        &void,
                        ogeom_mesh::Deflection::for_tolerances(self.tol),
                        self.tol,
                    )
                    .ok()
                    .filter(ogeom_topo::Triangulation::is_closed)
                    .map(|mesh| mesh.volume());
                    if encloses.is_some_and(|v| v > 0.0) {
                        self.report.warnings.push(format!(
                            "#{id}: void shell #{void_id} was written the wrong way \
                             out; it is read as the cavity it bounds"
                        ));
                        shells.push(void.reversed());
                    } else {
                        shells.push(void);
                    }
                }
                Ok(None) => self.report.warnings.push(format!(
                    "#{id}: void shell #{void_id} has no readable faces; the \
                     cavity is missing from the solid"
                )),
                Err(refusal) => self.report.warnings.push(format!(
                    "#{id}: void shell #{void_id} could not be read ({refusal}); \
                     the cavity is missing from the solid"
                )),
            }
        }
        Ok(make_solid(&mut self.model, &shells)?.shape)
    }

    /// The shape a `SHELL_BASED_SURFACE_MODEL` names: each of its shells
    /// sewn from its faces, a compound of them when there are several.
    ///
    /// A surface body is what a modeller exports for a part built from
    /// faces rather than from a solid, often one single-face body per face,
    /// and a reader that walks only `MANIFOLD_SOLID_BREP` leaves such a
    /// part invisible. The
    /// shells stay shells: the file did not call them solids, and a closed
    /// one is still the file's surface model, not this reader's promotion.
    pub(super) fn surface_model(&mut self, id: u64) -> OgeomResult<Shape> {
        let args = self.args(id, "SHELL_BASED_SURFACE_MODEL")?;
        let shell_ids: Vec<u64> = args
            .get(1)
            .and_then(Arg::list)
            .unwrap_or(&[])
            .iter()
            .filter_map(Arg::reference)
            .collect();
        let mut shells = Vec::new();
        for shell_id in shell_ids {
            if let Some(shell) = self.shell(shell_id)? {
                shells.push(shell);
            }
        }
        match shells.len() {
            0 => ogeom_bail!(
                Construction,
                "#{id}: a surface model with no readable faces"
            ),
            1 => Ok(shells.remove(0)),
            _ => Ok(ogeom_algo::build::make_compound(&mut self.model, &shells)?.shape),
        }
    }

    /// The faces of a `CLOSED_SHELL` or `OPEN_SHELL`, sewn; `None` when not
    /// one of them could be read.
    ///
    /// An `ORIENTED_CLOSED_SHELL` is a use of another shell the other way
    /// round (how a solid's voids are named) and resolves to that shell,
    /// reversed when the use says so.
    fn shell(&mut self, shell_id: u64) -> OgeomResult<Option<Shape>> {
        self.nested(|reader| reader.shell_at(shell_id))
    }

    fn shell_at(&mut self, shell_id: u64) -> OgeomResult<Option<Shape>> {
        let shell_instance = self.instance(shell_id)?;
        if let Some(oriented) = shell_instance.part("ORIENTED_CLOSED_SHELL") {
            let Some(base) = oriented.get(2).and_then(Arg::reference) else {
                ogeom_bail!(Construction, "#{shell_id}: an oriented shell names none");
            };
            let forward = !oriented.get(3).is_some_and(|a| a.is_enum("F"));
            return Ok(self
                .shell(base)?
                .map(|shell| if forward { shell } else { shell.reversed() }));
        }
        let shell_args = shell_instance
            .part("CLOSED_SHELL")
            .or_else(|| shell_instance.part("OPEN_SHELL"))
            .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "#{shell_id} is not a shell"))?
            .to_vec();
        let face_ids: Vec<u64> = shell_args
            .get(1)
            .and_then(Arg::list)
            .unwrap_or(&[])
            .iter()
            .filter_map(Arg::reference)
            .collect();
        self.cut_at_poles(&face_ids);
        self.split_band_rims(&face_ids);
        self.prepare_pcurves(&face_ids);
        let mut read = Vec::new();
        for fid in face_ids {
            if let Some(face) = self.face(fid)? {
                read.push((fid, face));
            }
        }
        if read.is_empty() {
            return Ok(None);
        }
        self.turn_backward_loops(shell_id, &mut read)?;
        let faces: Vec<Shape> = read.into_iter().map(|(_, face)| face).collect();
        Ok(Some(make_shell(&mut self.model, &faces)?.shape))
    }

    /// Walk back the loops of the faces a shell's own edges say run
    /// against their sense.
    ///
    /// Each edge between two faces is walked once each way, every face
    /// keeping its material on the left of its loops. A file whose faces
    /// turned against their surfaces were written with their loops walked
    /// backward has such faces walking every edge the same way as the face
    /// beside them. Where the shell's walks pick out a set of faces whose
    /// loops, walked back, leave every edge walked once each way, and every
    /// face in the set is turned against its surface, those faces are
    /// rebuilt with their loops walked back, and the report says how many.
    /// A shell whose walks pick out no such set is left as read.
    fn turn_backward_loops(&mut self, shell_id: u64, read: &mut [(u64, Shape)]) -> OgeomResult<()> {
        // Each edge's uses: the face, and whether it walks the edge forward.
        let mut uses: FastMap<(ogeom_topo::TShapeId, Location), Vec<(usize, bool)>> =
            FastMap::default();
        for (i, (_, face)) in read.iter().enumerate() {
            for wire in self.model.children_of(face)? {
                for edge in self.model.children_of(&wire)? {
                    if self
                        .model
                        .node(&edge)
                        .and_then(|n| n.data().as_edge())
                        .is_some_and(|d| d.degenerate)
                    {
                        continue;
                    }
                    uses.entry((edge.node(), edge.location().clone()))
                        .or_default()
                        .push((i, edge.orientation() == ogeom_topo::Orientation::Forward));
                }
            }
        }
        // Whether each face's loops must be walked back relative to a
        // reference face of its group: walked the same way by two faces,
        // an edge asks one of the two to turn; once each way, both or
        // neither.
        let pairs: Vec<(usize, usize, bool)> = uses
            .values()
            .filter_map(|u| match u.as_slice() {
                [(a, x), (b, y)] if a != b => Some((*a, *b, x == y)),
                _ => None,
            })
            .collect();
        if !pairs.iter().any(|&(_, _, same)| same) {
            return Ok(());
        }
        let mut beside: Vec<Vec<(usize, bool)>> = vec![Vec::new(); read.len()];
        for &(a, b, same) in &pairs {
            beside[a].push((b, same));
            beside[b].push((a, same));
        }
        let mut turn: Vec<Option<bool>> = vec![None; read.len()];
        let mut groups: Vec<Vec<usize>> = Vec::new();
        for start in 0..read.len() {
            if turn[start].is_some() {
                continue;
            }
            turn[start] = Some(false);
            let mut group = vec![start];
            let mut k = 0;
            while k < group.len() {
                let f = group[k];
                k += 1;
                let mine = turn[f] == Some(true);
                for &(g, same) in &beside[f] {
                    let want = mine ^ same;
                    match turn[g] {
                        None => {
                            turn[g] = Some(want);
                            group.push(g);
                        }
                        Some(had) if had != want => return Ok(()),
                        Some(_) => {}
                    }
                }
            }
            groups.push(group);
        }
        let backward = |i: usize| read[i].1.orientation() == ogeom_topo::Orientation::Reversed;
        let mut walk_back = Vec::new();
        for group in groups {
            let (marked, unmarked): (Vec<usize>, Vec<usize>) =
                group.iter().partition(|&&f| turn[f] == Some(true));
            let candidates = [marked, unmarked];
            let Some(chosen) = candidates
                .iter()
                .filter(|set| set.iter().all(|&f| backward(f)))
                .min_by_key(|set| set.len())
            else {
                return Ok(());
            };
            walk_back.extend(chosen.iter().copied());
        }
        if walk_back.is_empty() {
            return Ok(());
        }
        for &i in &walk_back {
            let (fid, face) = read[i].clone();
            let built = face.oriented(ogeom_topo::Orientation::Forward);
            let Some(surface_id) = self
                .model
                .node(&built)
                .and_then(|n| n.data().as_face())
                .map(|d| d.surface)
            else {
                ogeom_bail!(Construction, "#{fid}: a read face holds no face data");
            };
            let wires = self.model.ordered_children_of(&built)?;
            let wires = walked_back(&mut self.model, &wires, self.tol)?;
            let rebuilt = make_face_on(&mut self.model, surface_id, &wires, self.tol)?.shape;
            let turned = rebuilt.oriented(face.orientation());
            self.faces.insert(fid, turned.clone());
            read[i].1 = turned;
        }
        self.report.warnings.push(format!(
            "#{shell_id}: {} face(s) turned against their surfaces walked their \
             loops backward; their loops were walked back",
            walk_back.len()
        ));
        Ok(())
    }
}
