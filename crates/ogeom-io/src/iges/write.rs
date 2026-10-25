//! From a document to an IGES deck.
//!
//! The writer emits manifold solid B-rep objects (entity 186 over shells
//! (514), faces (510), loops (508), one edge list (504) and one vertex list
//! (502) per solid), and sheets (shells no solid owns, faces no shell owns)
//! as trimmed surfaces (144), one per face, each bounded by curves on the
//! surface (142) over its model-space edges, and wireframe (wires and edges
//! outside every face) as independent curves. Surfaces go out in their
//! analytic spellings where IGES has one and as rational B-splines (128)
//! where it does not. Curves defined in a
//! plane of their own, arcs and ellipses, are written in definition space
//! with a transformation matrix (124) carrying them to model space, which is
//! how the format wants them.
//!
//! Everything is written in millimetres, model space, with geometry baked:
//! the same decision the STEP writer made: a file carries positions, not this
//! kernel's location chains.

use ogeom_core::FastMap;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::Reversible as _;
use ogeom_geom::Surface as _;
use ogeom_geom::Transformable as _;
use ogeom_geom::{BSplineCurve, Curve, SurfaceGeometry};
use ogeom_math::{Frame, Point, Transform};
use ogeom_topo::{EdgeRepr, Filter, NodeData, Shape, ShapeType, explore};

/// Write a document's solids and sheets as an IGES file.
///
/// A solid becomes a manifold solid B-rep object (186). A sheet's faces
/// become independent trimmed surfaces (144), which a reader sews back into
/// the sheet they were; a sheet that closes reads back as a solid.
/// A wire or edge outside every face is wireframe, written as its curve
/// (an edge's own, a composite curve (102) of a wire's edges in the order
/// the wire walks them) standing independent.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// document holds no part, a part holds a vertex outside every edge or
/// nothing at all, a wire holds only degenerate edges, or a body carries
/// geometry with no IGES spelling here; the error names the entity and the
/// capability row.
pub fn write_iges(document: &ogeom_doc::Document, tol: Tolerances) -> OgeomResult<String> {
    let mut writer = Writer {
        model: document.model(),
        entities: Vec::new(),
        vertices: FastMap::default(),
        edges: FastMap::default(),
        vertex_coords: Vec::new(),
        edge_records: Vec::new(),
        tol,
    };

    let mut parts = Vec::new();
    for (_, product) in document.products() {
        let ogeom_doc::ProductKind::Part { shape } = &product.kind else {
            continue;
        };
        let bodies = crate::bodies::bodies_of(writer.model, shape, "IGES", true)?;
        parts.push((product.name.clone(), bodies));
    }
    if parts.is_empty() {
        ogeom_bail!(Construction, "the document holds no part to write as IGES");
    }
    for (label, bodies) in parts {
        for solid in &bodies.solids {
            writer.solid(solid, &label)?;
        }
        for shell in &bodies.shells {
            for face in writer.model.ordered_children_of(shell)? {
                writer.trimmed_surface(&face, &label)?;
            }
        }
        for face in &bodies.faces {
            writer.trimmed_surface(face, &label)?;
        }
        for wire in &bodies.wires {
            writer.wireframe(wire, &label)?;
        }
    }
    Ok(writer.serialize())
}

/// One pending entity: directory fields and parameter text.
struct Pending {
    kind: i64,
    form: i64,
    /// Index into `entities` of a 124 transform, or `None`.
    transform: Option<usize>,
    /// Independent (top-level) or physically subordinate.
    independent: bool,
    label: String,
    params: String,
}

struct Writer<'a> {
    model: &'a ogeom_topo::Model,
    entities: Vec<Pending>,
    /// Vertex index (1-based) in the current solid's 502, keyed by node
    /// *and position*: an instanced node placed twice is two vertices in the
    /// file, and keying by node alone would weld a prism's top to its bottom.
    vertices: FastMap<(ogeom_topo::TShapeId, [u64; 3]), usize>,
    /// Edge index (1-based) in the current solid's 504, keyed by node and
    /// placement for the same reason.
    edges: FastMap<(ogeom_topo::TShapeId, [u64; 3]), usize>,
    /// The current solid's vertex coordinates, in 502 order.
    vertex_coords: Vec<Point>,
    /// The current solid's edge records: (curve entity, start index, end index).
    edge_records: Vec<(usize, usize, usize)>,
    tol: Tolerances,
}

impl Writer<'_> {
    fn push(&mut self, p: Pending) -> usize {
        self.entities.push(p);
        self.entities.len() - 1
    }

    /// The directory pointer an entity index will serialize as.
    #[allow(clippy::cast_possible_wrap, reason = "entity counts are small")]
    fn de(&self, index: usize) -> i64 {
        2 * index as i64 + 1
    }

    fn solid(&mut self, solid: &Shape, label: &str) -> OgeomResult<()> {
        // Fresh per-solid lists; their entities are created after the faces
        // so their contents are complete, and patched into the loops by a
        // placeholder scheme below.
        self.vertices.clear();
        self.edges.clear();
        self.vertex_coords.clear();
        self.edge_records.clear();

        // The first shell is the outer boundary and each further one bounds
        // a cavity. A void is written as the standard names it, a shell
        // facing out of the cavity used the other way round, which leaves
        // every face of the solid facing away from the material.
        let shells = explore(self.model, solid, Filter::OfType(ShapeType::Shell))?;
        if shells.is_empty() {
            ogeom_bail!(Construction, "a solid with no shell cannot be written");
        }
        let mut shell_faces = Vec::with_capacity(shells.len());
        for (k, shell) in shells.iter().enumerate() {
            let mut face_entities = Vec::new();
            for face in self.model.ordered_children_of(shell)? {
                let face = if k == 0 { face } else { face.reversed() };
                let (entity, turned) = self.face(&face)?;
                let forward = face.orientation() != ogeom_topo::Orientation::Reversed;
                face_entities.push((entity, forward != turned));
            }
            shell_faces.push(face_entities);
        }

        // The lists exist in full.
        let vertex_list = {
            let mut params = format!("{}", self.vertex_coords.len());
            for p in &self.vertex_coords {
                params.push_str(&format!(",{},{},{}", fmt(p.x), fmt(p.y), fmt(p.z)));
            }
            self.push(Pending {
                kind: 502,
                form: 1,
                transform: None,
                independent: false,
                label: String::new(),
                params,
            })
        };
        let edge_list = {
            let records = std::mem::take(&mut self.edge_records);
            let mut params = format!("{}", records.len());
            let mut curve_des = Vec::new();
            for (curve_entity, sv, tv) in &records {
                curve_des.push(self.de(*curve_entity));
                params.push_str(&format!(",{},@V,{sv},@V,{tv}", self.de(*curve_entity)));
            }
            self.push(Pending {
                kind: 504,
                form: 1,
                transform: None,
                independent: false,
                label: String::new(),
                params,
            })
        };
        // Loops referred to the lists before the lists existed; the
        // placeholders resolve here.
        let vlist_de = self.de(vertex_list);
        let elist_de = self.de(edge_list);
        for e in &mut self.entities {
            if e.kind == 508 || e.kind == 504 {
                e.params = e.params.replace("@E", &elist_de.to_string());
                e.params = e.params.replace("@V", &vlist_de.to_string());
            }
        }

        let mut shell_des = Vec::with_capacity(shell_faces.len());
        for face_entities in &shell_faces {
            let mut params = format!("{}", face_entities.len());
            for (face, forward) in face_entities {
                let flag = i32::from(*forward);
                params.push_str(&format!(",{},{flag}", self.de(*face)));
            }
            let shell_entity = self.push(Pending {
                kind: 514,
                form: 1,
                transform: None,
                independent: false,
                label: String::new(),
                params,
            });
            shell_des.push(self.de(shell_entity));
        }
        let mut params = format!("{},1,{}", shell_des[0], shell_des.len() - 1);
        for void in &shell_des[1..] {
            params.push_str(&format!(",{void},0"));
        }
        // The part's whole name, past the label's eight columns, as the
        // solid's Name property (406 form 15): no back pointers, one
        // property.
        if let Some(name) = self.name_property(label) {
            params.push_str(&format!(",0,1,{}", self.de(name)));
        }
        self.push(Pending {
            kind: 186,
            form: 0,
            transform: None,
            independent: true,
            label: label.chars().take(8).collect(),
            params,
        });
        Ok(())
    }

    /// A Name property (406 form 15) holding `name`, or `None` for an empty
    /// name or one a string parameter cannot carry whole here: past ASCII,
    /// or holding a delimiter, which a record break could fall beside.
    fn name_property(&mut self, name: &str) -> Option<usize> {
        let plain = !name.is_empty()
            && name
                .chars()
                .all(|c| (c.is_ascii_graphic() || c == ' ') && c != ',' && c != ';');
        plain.then(|| {
            self.push(Pending {
                kind: 406,
                form: 15,
                transform: None,
                independent: false,
                label: String::new(),
                params: format!("1,{}H{name}", name.len()),
            })
        })
    }

    /// A face (510) over its surface written in model space, its loops
    /// running about the file's surface normal, and whether the face faces
    /// against that normal (see [`crate::bodies::turned_in_file`]): the shell's flag for the
    /// face is its orientation, turned once more where it does.
    fn face(&mut self, face: &Shape) -> OgeomResult<(usize, bool)> {
        let placement = face.transform(self.model.datums())?;
        let surface = {
            let Some(node) = self.model.node(face) else {
                ogeom_bail!(Dangling, "face is not in this model");
            };
            let NodeData::Face(data) = node.data() else {
                ogeom_bail!(Construction, "face node holds no face data");
            };
            let Some(surface) = self.model.geometry().surface(data.surface) else {
                ogeom_bail!(Dangling, "face refers to a surface not in this model");
            };
            surface.clone().transformed(&placement, self.tol)?
        };
        let restated = self.restates(&surface)?;
        let turned = crate::bodies::turned_in_file(&placement, restated);
        let surface_entity = self.surface(&surface)?;
        // A face's loops keep it on their left about its surface's normal,
        // as the face stores them, and run the other way about a normal the
        // file restates; the shell's flag says which side of the surface the
        // face faces. The outer loop stays first.
        let walked = if restated {
            ogeom_topo::Orientation::Reversed
        } else {
            ogeom_topo::Orientation::Forward
        };
        let wires = self.model.children_of(&face.oriented(walked))?;
        let mut loop_entities = Vec::new();
        for wire in &wires {
            if let Some(entity) = self.wire(wire)? {
                loop_entities.push(entity);
            }
        }
        if loop_entities.is_empty() {
            ogeom_bail!(
                Construction,
                "a face whose every boundary is degenerate cannot be written \
                 as IGES; see docs/CAPABILITIES.md, io.iges"
            );
        }
        let mut params = format!("{},{},1", self.de(surface_entity), loop_entities.len());
        for entity in &loop_entities {
            params.push_str(&format!(",{}", self.de(*entity)));
        }
        let entity = self.push(Pending {
            kind: 510,
            form: 1,
            transform: None,
            independent: false,
            label: String::new(),
            params,
        });
        Ok((entity, turned))
    }

    /// A face as an independent trimmed surface (144): its surface, and one
    /// curve on the surface (142) per boundary over the model-space curves
    /// of the boundary's edges.
    ///
    /// The format gives a trimmed surface no sense of its own: its normal is
    /// its surface's. A reversed face therefore goes out on its surface
    /// turned over, a plane with its normal negated and anything else as its
    /// exact B-spline patch with the two parameter directions exchanged.
    fn trimmed_surface(&mut self, face: &Shape, label: &str) -> OgeomResult<()> {
        let placement = face.transform(self.model.datums())?;
        let surface = {
            let Some(node) = self.model.node(face) else {
                ogeom_bail!(Dangling, "face is not in this model");
            };
            let NodeData::Face(data) = node.data() else {
                ogeom_bail!(Construction, "face node holds no face data");
            };
            let Some(surface) = self.model.geometry().surface(data.surface) else {
                ogeom_bail!(Dangling, "face refers to a surface not in this model");
            };
            surface.clone().transformed(&placement, self.tol)?
        };
        let surface_entity = if face.orientation() == ogeom_topo::Orientation::Reversed {
            self.turned_surface(&surface)?
        } else {
            self.surface(&surface)?
        };
        // Stored order, the outer boundary first; the walk of a reversed
        // face lists its holes first.
        let mut boundaries = Vec::new();
        for wire in self.model.children_of(face)? {
            if let Some(curve) = self.boundary_curve(&wire)? {
                // Unspecified creation, no parameter-space curve, the
                // model-space curve preferred.
                let params = format!("0,{},0,{},2", self.de(surface_entity), self.de(curve));
                boundaries.push(self.push(Pending {
                    kind: 142,
                    form: 0,
                    transform: None,
                    independent: false,
                    label: String::new(),
                    params,
                }));
            }
        }
        // A face with no wire at all is its whole surface: the trimmed
        // surface whose outer boundary is the surface's own.
        if self.model.children_of(face)?.is_empty() {
            self.push(Pending {
                kind: 144,
                form: 0,
                transform: None,
                independent: true,
                label: label.chars().take(8).collect(),
                params: format!("{},0,0,0", self.de(surface_entity)),
            });
            return Ok(());
        }
        let Some((outer, inner)) = boundaries.split_first() else {
            ogeom_bail!(
                Construction,
                "a face whose every boundary is degenerate cannot be written \
                 as IGES; see docs/CAPABILITIES.md, io.iges"
            );
        };
        let mut params = format!(
            "{},1,{},{}",
            self.de(surface_entity),
            inner.len(),
            self.de(*outer)
        );
        for boundary in inner {
            params.push_str(&format!(",{}", self.de(*boundary)));
        }
        self.push(Pending {
            kind: 144,
            form: 0,
            transform: None,
            independent: true,
            label: label.chars().take(8).collect(),
            params,
        });
        Ok(())
    }

    /// A wire's edges as one model-space curve: the edge's own curve when
    /// there is one, a composite curve (102) of them when there are more,
    /// `None` when every edge is degenerate. Each piece runs the way the wire
    /// walks it, in the wire's order: the composite curve has no sense flag
    /// for its pieces, so an edge the wire uses reversed goes out as its
    /// curve reversed.
    fn boundary_curve(&mut self, wire: &Shape) -> OgeomResult<Option<usize>> {
        let edges = self.model.ordered_children_of(wire)?;
        self.chained_curve(&edges)
    }

    /// Wireframe, a free wire or edge, as its curve standing independent.
    fn wireframe(&mut self, shape: &Shape, label: &str) -> OgeomResult<()> {
        let edges = if self.model.kind_of(shape)? == ShapeType::Edge {
            vec![shape.clone()]
        } else {
            self.model.ordered_children_of(shape)?
        };
        let Some(curve) = self.chained_curve(&edges)? else {
            ogeom_bail!(
                Construction,
                "a free wire whose every edge is degenerate cannot be written as IGES"
            );
        };
        let entity = &mut self.entities[curve];
        entity.independent = true;
        entity.label = label.chars().take(8).collect();
        Ok(())
    }

    /// Edges as one model-space curve, as [`Self::boundary_curve`] writes a
    /// wire's.
    fn chained_curve(&mut self, edges: &[Shape]) -> OgeomResult<Option<usize>> {
        let mut pieces = Vec::new();
        for edge in edges {
            let forward = edge.orientation() != ogeom_topo::Orientation::Reversed;
            let (curve, range) = {
                let Some(data) = self.model.node(edge).and_then(|n| n.data().as_edge()) else {
                    ogeom_bail!(Construction, "edge node holds no edge data");
                };
                if data.degenerate {
                    continue;
                }
                let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                    ogeom_bail!(Construction, "an edge with no curve cannot be written");
                };
                let Some(geometry) = self.model.geometry().curve(*curve) else {
                    ogeom_bail!(Dangling, "curve is not in this model");
                };
                let placement = edge.transform(self.model.datums())?;
                (geometry.clone().transformed(&placement, self.tol)?, *range)
            };
            pieces.push(self.curve(&curve, range, forward)?);
        }
        match pieces.as_slice() {
            [] => Ok(None),
            [only] => Ok(Some(*only)),
            _ => {
                let mut params = format!("{}", pieces.len());
                for piece in &pieces {
                    params.push_str(&format!(",{}", self.de(*piece)));
                }
                Ok(Some(self.push(Pending {
                    kind: 102,
                    form: 0,
                    transform: None,
                    independent: false,
                    label: String::new(),
                    params,
                })))
            }
        }
    }

    /// A surface with its normal turned over.
    fn turned_surface(&mut self, surface: &SurfaceGeometry) -> OgeomResult<usize> {
        match surface {
            SurfaceGeometry::Plane(p) => {
                let frame = p.plane().frame();
                let point = self.point_entity(frame.origin());
                let normal = self.direction_entity(-frame.z());
                let params = format!("{},{}", self.de(point), self.de(normal));
                Ok(self.push(Pending {
                    kind: 190,
                    form: 0,
                    transform: None,
                    independent: false,
                    label: String::new(),
                    params,
                }))
            }
            SurfaceGeometry::BSpline(b) => Ok(self.nurbs_surface(b, true)),
            other => {
                if let Some(entity) = self.swept_surface(other, false)? {
                    return Ok(entity);
                }
                let bspline = other.to_bspline(self.tol)?;
                Ok(self.nurbs_surface(&bspline, true))
            }
        }
    }

    /// Whether this writer restates `surface`'s normal (see
    /// [`crate::bodies::restates`]): an offset goes out as its analytic
    /// form where it has one, or over its plane or spline basis, and either
    /// says for it; any other offset goes out as a spline patch.
    fn restates(&self, surface: &SurfaceGeometry) -> OgeomResult<bool> {
        if let SurfaceGeometry::Offset(o) = surface {
            if let Some(analytic) = o.analytic(self.tol)? {
                return Ok(crate::bodies::restates(&analytic, false));
            }
            if matches!(
                o.basis(),
                SurfaceGeometry::Plane(_) | SurfaceGeometry::BSpline(_)
            ) {
                return Ok(crate::bodies::restates(o.basis(), false));
            }
        }
        Ok(crate::bodies::restates(surface, false))
    }

    /// A surface of revolution (120), tabulated cylinder (122) or offset
    /// surface (140) in its own spelling, its normal turned over where not
    /// `forward`: a revolution's generatrix and a tabulated cylinder's
    /// directrix written the other way round, an offset over its basis
    /// turned over. `None` for an offset whose basis is not a plane or a
    /// B-spline patch and has no analytic form, and for any other surface.
    fn swept_surface(
        &mut self,
        surface: &SurfaceGeometry,
        forward: bool,
    ) -> OgeomResult<Option<usize>> {
        use SurfaceGeometry as S;
        let (kind, params) = match surface {
            // The axis as a line, the generatrix over its whole range, and
            // the turn as the start and terminate angles.
            S::Revolution(r) => {
                let ((a0, a1), (c0, c1)) = r.domain();
                let axis = r.axis();
                let (from, to) = (axis.location, axis.location + axis.direction.vector());
                let line = self.push(Pending {
                    kind: 110,
                    form: 0,
                    transform: None,
                    independent: false,
                    label: String::new(),
                    params: format!(
                        "{},{},{},{},{},{}",
                        fmt(from.x),
                        fmt(from.y),
                        fmt(from.z),
                        fmt(to.x),
                        fmt(to.y),
                        fmt(to.z)
                    ),
                });
                let generatrix = self.curve(r.curve(), (c0, c1), forward)?;
                (
                    120,
                    format!(
                        "{},{},{},{}",
                        self.de(line),
                        self.de(generatrix),
                        fmt(a0),
                        fmt(a1)
                    ),
                )
            }
            // The directrix at the start of the sweep's window, and the far
            // end of the generator through the directrix's first point.
            S::Extrusion(e) => {
                let ((u0, u1), (v0, v1)) = e.domain();
                let along = e.direction().vector();
                let directrix = e
                    .curve()
                    .transformed(&Transform::translation(along * v0), self.tol)?;
                let start = directrix.point_at(if forward { u0 } else { u1 }, self.tol)?;
                let far = start + along * (v1 - v0);
                let curve = self.curve(&directrix, (u0, u1), forward)?;
                (
                    122,
                    format!(
                        "{},{},{},{}",
                        self.de(curve),
                        fmt(far.x),
                        fmt(far.y),
                        fmt(far.z)
                    ),
                )
            }
            S::Offset(o) => {
                if let Some(analytic) = o.analytic(self.tol)? {
                    let entity = if forward {
                        self.surface(&analytic)?
                    } else {
                        self.turned_surface(&analytic)?
                    };
                    return Ok(Some(entity));
                }
                let basis = o.basis();
                if !matches!(basis, S::Plane(_) | S::BSpline(_)) {
                    return Ok(None);
                }
                // The reader takes the offset's side from the basis normal
                // at the middle of its domain, which a plane and a written
                // B-spline patch keep (a patch written turned over swaps its
                // parameters, the middle staying where it was).
                let ((ua, ub), (va, vb)) = basis.domain();
                let normal = basis
                    .normal_at(f64::midpoint(ua, ub), f64::midpoint(va, vb), self.tol)?
                    .vector();
                let (written, normal, distance) = if forward {
                    (self.surface(basis)?, normal, o.distance())
                } else {
                    (self.turned_surface(basis)?, -normal, -o.distance())
                };
                let indicator = if distance < 0.0 { -normal } else { normal };
                (
                    140,
                    format!(
                        "{},{},{},{},{}",
                        fmt(indicator.x),
                        fmt(indicator.y),
                        fmt(indicator.z),
                        fmt(distance.abs()),
                        self.de(written)
                    ),
                )
            }
            _ => return Ok(None),
        };
        Ok(Some(self.push(Pending {
            kind,
            form: 0,
            transform: None,
            independent: false,
            label: String::new(),
            params,
        })))
    }

    /// A wire as a loop (508), or `None` when every edge in it is degenerate:
    /// a pole has no curve, and the reader rebuilds chart degeneracies
    /// from the surface itself.
    fn wire(&mut self, wire: &Shape) -> OgeomResult<Option<usize>> {
        let children = self.model.ordered_children_of(wire)?;
        let degenerate = |edge: &Shape| {
            self.model
                .node(edge)
                .and_then(|n| n.data().as_edge())
                .is_some_and(|d| d.degenerate)
        };
        let live: Vec<&Shape> = children.iter().filter(|e| !degenerate(e)).collect();
        if live.is_empty() {
            return Ok(None);
        }
        let mut entries = Vec::new();
        for edge in live {
            let index = self.edge(edge)?;
            let flag = i32::from(edge.orientation() != ogeom_topo::Orientation::Reversed);
            entries.push(format!("0,@E,{index},{flag},0"));
        }
        let params = format!("{},{}", entries.len(), entries.join(","));
        Ok(Some(self.push(Pending {
            kind: 508,
            form: 1,
            transform: None,
            independent: false,
            label: String::new(),
            params,
        })))
    }

    /// The edge's index in the current solid's 504 list, creating it once
    /// per placement.
    fn edge(&mut self, edge: &Shape) -> OgeomResult<usize> {
        let placement = edge.transform(self.model.datums())?;
        let key = (edge.node(), transform_bits(&placement));
        if let Some(&index) = self.edges.get(&key) {
            return Ok(index);
        }
        let (curve, range) = {
            let Some(data) = self.model.node(edge).and_then(|n| n.data().as_edge()) else {
                ogeom_bail!(Construction, "edge node holds no edge data");
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                ogeom_bail!(Construction, "an edge with no curve cannot be written");
            };
            let Some(geometry) = self.model.geometry().curve(*curve) else {
                ogeom_bail!(Dangling, "curve is not in this model");
            };
            (geometry.clone().transformed(&placement, self.tol)?, *range)
        };
        let curve_entity = self.curve(&curve, range, true)?;
        let vertices = self.model.children_of(edge)?;
        let (from, to) = match vertices.len() {
            0 => ogeom_bail!(Construction, "an edge with no vertices cannot be written"),
            1 => (vertices[0].clone(), vertices[0].clone()),
            _ => (vertices[0].clone(), vertices[vertices.len() - 1].clone()),
        };
        // Each vertex carries its own composed placement; children_of has
        // already folded the edge's chain in, and an instanced vertex (a
        // prism's top corner is its bottom corner, moved) adds a hop of its
        // own that the edge's placement alone would drop.
        let sv = {
            let at = from.transform(self.model.datums())?;
            self.vertex(&from, &at)?
        };
        let tv = {
            let at = to.transform(self.model.datums())?;
            self.vertex(&to, &at)?
        };
        self.edge_records.push((curve_entity, sv, tv));
        let index = self.edge_records.len();
        self.edges.insert(key, index);
        Ok(index)
    }

    /// The vertex's index in the current solid's 502 list, creating it once.
    fn vertex(&mut self, vertex: &Shape, placement: &Transform) -> OgeomResult<usize> {
        let Some(data) = self.model.node(vertex).and_then(|n| n.data().as_vertex()) else {
            ogeom_bail!(Construction, "vertex node holds no vertex data");
        };
        let at = placement.apply(data.point);
        let key = (
            vertex.node(),
            [at.x.to_bits(), at.y.to_bits(), at.z.to_bits()],
        );
        if let Some(&index) = self.vertices.get(&key) {
            return Ok(index);
        }
        self.vertex_coords.push(at);
        let index = self.vertex_coords.len();
        self.vertices.insert(key, index);
        Ok(index)
    }

    /// A curve over a range, in its IGES spelling, running from `range.0`
    /// to `range.1` when `forward` and the other way when not.
    fn curve(&mut self, curve: &Curve, range: (f64, f64), forward: bool) -> OgeomResult<usize> {
        // The angle on the underlying conic at each end of the walk.
        let walk = |reversed: bool| {
            let angle = |u: f64| if reversed { -u } else { u };
            let (s, e) = if forward { range } else { (range.1, range.0) };
            // A conic arc runs counter-clockwise about its frame's `z`; one
            // walked clockwise goes out about the frame turned over `x`,
            // where each angle negates.
            let counter_clockwise = reversed != forward;
            let sign = if counter_clockwise { 1.0 } else { -1.0 };
            (sign * angle(s), sign * angle(e), counter_clockwise)
        };
        match curve {
            Curve::Line(line) => {
                let (from, to) = if forward { range } else { (range.1, range.0) };
                let a = line.point_at(from, self.tol)?;
                let b = line.point_at(to, self.tol)?;
                Ok(self.push(Pending {
                    kind: 110,
                    form: 0,
                    transform: None,
                    independent: false,
                    label: String::new(),
                    params: format!(
                        "{},{},{},{},{},{}",
                        fmt(a.x),
                        fmt(a.y),
                        fmt(a.z),
                        fmt(b.x),
                        fmt(b.y),
                        fmt(b.z)
                    ),
                }))
            }
            Curve::Circle(c) => {
                let circle = c.circle();
                let r = circle.radius();
                let (s, e, counter_clockwise) = walk(c.is_reversed());
                let frame = if counter_clockwise {
                    circle.frame()
                } else {
                    circle.frame().with_z_reversed()
                };
                let transform = self.transform_entity(&frame);
                let params = format!(
                    "0.,0.,0.,{},{},{},{}",
                    fmt(r * s.cos()),
                    fmt(r * s.sin()),
                    fmt(r * e.cos()),
                    fmt(r * e.sin()),
                );
                Ok(self.push(Pending {
                    kind: 100,
                    form: 0,
                    transform: Some(transform),
                    independent: false,
                    label: String::new(),
                    params,
                }))
            }
            Curve::Ellipse(el) => {
                let ellipse = el.ellipse();
                let (a, b) = (ellipse.major_radius(), ellipse.minor_radius());
                let (s, e, counter_clockwise) = walk(el.is_reversed());
                let frame = if counter_clockwise {
                    ellipse.frame()
                } else {
                    ellipse.frame().with_z_reversed()
                };
                let transform = self.transform_entity(&frame);
                let at = |t: f64| (a * t.cos(), b * t.sin());
                let (sx, sy) = at(s);
                let (ex, ey) = at(e);
                // x²/a² + y²/b² − 1 = 0, spelt in the general coefficients.
                let params = format!(
                    "{},0.,{},0.,0.,-1.,0.,{},{},{},{}",
                    fmt(1.0 / (a * a)),
                    fmt(1.0 / (b * b)),
                    fmt(sx),
                    fmt(sy),
                    fmt(ex),
                    fmt(ey),
                );
                Ok(self.push(Pending {
                    kind: 104,
                    form: 1,
                    transform: Some(transform),
                    independent: false,
                    label: String::new(),
                    params,
                }))
            }
            // A reversed trim runs its basis backwards, each parameter
            // mirrored in its own domain.
            Curve::Trimmed(t) if t.is_reversed() => {
                let (lo, hi) = ogeom_geom::Curve3d::domain(curve);
                self.curve(t.basis(), (lo + hi - range.1, lo + hi - range.0), !forward)
            }
            Curve::Trimmed(t) => self.curve(t.basis(), range, forward),
            Curve::BSpline(b) if forward => self.nurbs_curve(b, range),
            Curve::BSpline(b) => {
                // Reversal keeps the domain and mirrors each parameter in it.
                let (lo, hi) = ogeom_geom::Curve3d::domain(b);
                let back = reversed_bspline(b)?;
                self.nurbs_curve(&back, (lo + hi - range.1, lo + hi - range.0))
            }
            other => {
                // The exact conversion carries anything with a closed NURBS
                // form; what has none (a helix) is refused by name there.
                let mut bspline = other.to_bspline_over(range, self.tol)?;
                if !forward {
                    bspline = reversed_bspline(&bspline)?;
                }
                self.nurbs_curve(&bspline, ogeom_geom::Curve3d::domain(&bspline))
            }
        }
    }

    fn nurbs_curve(&mut self, curve: &BSplineCurve, range: (f64, f64)) -> OgeomResult<usize> {
        let knots = curve.knots();
        let control = curve.control_points();
        let degree = knots.degree();
        let k = control.len() - 1;
        let rational = curve.is_rational();
        let closed = {
            let (lo, hi) = ogeom_geom::Curve3d::domain(curve);
            let a = curve.point_at(lo, self.tol)?;
            let b = curve.point_at(hi, self.tol)?;
            i32::from(a.distance(b) < self.tol.confusion())
        };
        let mut params = format!("{k},{degree},0,{closed},{},0", i32::from(!rational));
        for t in knots.knots() {
            params.push_str(&format!(",{}", fmt(*t)));
        }
        for w in control {
            params.push_str(&format!(",{}", fmt(w.weight)));
        }
        for w in control {
            let p = (*w).point();
            params.push_str(&format!(",{},{},{}", fmt(p.x), fmt(p.y), fmt(p.z)));
        }
        params.push_str(&format!(",{},{}", fmt(range.0), fmt(range.1)));
        Ok(self.push(Pending {
            kind: 126,
            form: 0,
            transform: None,
            independent: false,
            label: String::new(),
            params,
        }))
    }

    /// A definition-space→model transform (124) for a frame.
    fn transform_entity(&mut self, frame: &Frame) -> usize {
        let (x, y, z) = (frame.x().vector(), frame.y().vector(), frame.z().vector());
        let o = frame.origin();
        let params = format!(
            "{},{},{},{},{},{},{},{},{},{},{},{}",
            fmt(x.x),
            fmt(y.x),
            fmt(z.x),
            fmt(o.x),
            fmt(x.y),
            fmt(y.y),
            fmt(z.y),
            fmt(o.y),
            fmt(x.z),
            fmt(y.z),
            fmt(z.z),
            fmt(o.z),
        );
        self.push(Pending {
            kind: 124,
            form: 0,
            transform: None,
            independent: false,
            label: String::new(),
            params,
        })
    }

    /// A point entity (116).
    fn point_entity(&mut self, p: Point) -> usize {
        self.push(Pending {
            kind: 116,
            form: 0,
            transform: None,
            independent: false,
            label: String::new(),
            params: format!("{},{},{},0", fmt(p.x), fmt(p.y), fmt(p.z)),
        })
    }

    /// A direction entity (123).
    fn direction_entity(&mut self, d: ogeom_math::Direction) -> usize {
        let v = d.vector();
        self.push(Pending {
            kind: 123,
            form: 0,
            transform: None,
            independent: false,
            label: String::new(),
            params: format!("{},{},{}", fmt(v.x), fmt(v.y), fmt(v.z)),
        })
    }

    /// A surface in its IGES spelling.
    fn surface(&mut self, surface: &SurfaceGeometry) -> OgeomResult<usize> {
        use SurfaceGeometry as S;
        match surface {
            S::Plane(p) => {
                let frame = p.plane().frame();
                let point = self.point_entity(frame.origin());
                let normal = self.direction_entity(frame.z());
                let params = format!("{},{}", self.de(point), self.de(normal));
                Ok(self.push(Pending {
                    kind: 190,
                    form: 0,
                    transform: None,
                    independent: false,
                    label: String::new(),
                    params,
                }))
            }
            S::Cylinder(c) => {
                let cyl = c.cylinder();
                let frame = cyl.frame();
                let point = self.point_entity(frame.origin());
                let axis = self.direction_entity(frame.z());
                let params = format!("{},{},{}", self.de(point), self.de(axis), fmt(cyl.radius()));
                Ok(self.push(Pending {
                    kind: 192,
                    form: 0,
                    transform: None,
                    independent: false,
                    label: String::new(),
                    params,
                }))
            }
            S::Cone(c) => {
                let cone = c.cone();
                let frame = cone.frame();
                let point = self.point_entity(frame.origin());
                let axis = self.direction_entity(frame.z());
                let params = format!(
                    "{},{},{},{}",
                    self.de(point),
                    self.de(axis),
                    fmt(cone.reference_radius()),
                    fmt(cone.half_angle().to_degrees()),
                );
                Ok(self.push(Pending {
                    kind: 194,
                    form: 0,
                    transform: None,
                    independent: false,
                    label: String::new(),
                    params,
                }))
            }
            // Form 1, with the axis and the direction `u` starts from: a
            // sphere's poles are on its axis, and the face's poles and seam
            // with them.
            S::Sphere(s) => {
                let sphere = s.sphere();
                let frame = sphere.frame();
                let point = self.point_entity(sphere.centre());
                let axis = self.direction_entity(frame.z());
                let reference = self.direction_entity(frame.x());
                let params = format!(
                    "{},{},{},{}",
                    self.de(point),
                    fmt(sphere.radius()),
                    self.de(axis),
                    self.de(reference)
                );
                Ok(self.push(Pending {
                    kind: 196,
                    form: 1,
                    transform: None,
                    independent: false,
                    label: String::new(),
                    params,
                }))
            }
            S::Torus(t) => {
                let torus = t.torus();
                let frame = torus.frame();
                let point = self.point_entity(frame.origin());
                let axis = self.direction_entity(frame.z());
                let params = format!(
                    "{},{},{},{}",
                    self.de(point),
                    self.de(axis),
                    fmt(torus.major_radius()),
                    fmt(torus.minor_radius()),
                );
                Ok(self.push(Pending {
                    kind: 198,
                    form: 0,
                    transform: None,
                    independent: false,
                    label: String::new(),
                    params,
                }))
            }
            S::Revolution(_) | S::Extrusion(_) | S::Offset(_) => {
                match self.swept_surface(surface, true)? {
                    Some(entity) => Ok(entity),
                    None => {
                        let bspline = surface.to_bspline(self.tol)?;
                        self.surface(&SurfaceGeometry::BSpline(bspline))
                    }
                }
            }
            S::BSpline(b) => Ok(self.nurbs_surface(b, false)),
            other => {
                // Swept, offset and trimmed carriers convert exactly where a
                // closed NURBS form exists; the conversion refuses by name
                // where it does not.
                let bspline = other.to_bspline(self.tol)?;
                self.surface(&SurfaceGeometry::BSpline(bspline))
            }
        }
    }

    /// A B-spline surface (128). `exchanged` writes it with its two
    /// parameter directions swapped: the same points, the normal turned over.
    fn nurbs_surface(&mut self, b: &ogeom_geom::BSplineSurface, exchanged: bool) -> usize {
        let (mut uk, mut vk) = (b.u_knots(), b.v_knots());
        let grid = b.grid();
        let (mut nu, mut nv) = (grid.u_count(), grid.v_count());
        let ((mut u0, mut u1), (mut v0, mut v1)) = b.domain();
        if exchanged {
            std::mem::swap(&mut uk, &mut vk);
            std::mem::swap(&mut nu, &mut nv);
            std::mem::swap(&mut u0, &mut v0);
            std::mem::swap(&mut u1, &mut v1);
        }
        // The control net at the file's (u, v), whichever way it is written.
        let at = |u: usize, v: usize| {
            if exchanged {
                grid.get(v, u)
            } else {
                grid.get(u, v)
            }
        };
        let (k1, k2) = (nu - 1, nv - 1);
        let (m1, m2) = (uk.degree(), vk.degree());
        let mut params = format!(
            "{k1},{k2},{m1},{m2},0,0,{},0,0",
            i32::from(!b.is_rational())
        );
        for t in uk.knots() {
            params.push_str(&format!(",{}", fmt(*t)));
        }
        for t in vk.knots() {
            params.push_str(&format!(",{}", fmt(*t)));
        }
        // The file wants the first (u) index varying fastest.
        for v in 0..nv {
            for u in 0..nu {
                let w = at(u, v).map_or(1.0, |w| w.weight);
                params.push_str(&format!(",{}", fmt(w)));
            }
        }
        for v in 0..nv {
            for u in 0..nu {
                let p = at(u, v).map_or(Point::ORIGIN, ogeom_math::Weighted::point);
                params.push_str(&format!(",{},{},{}", fmt(p.x), fmt(p.y), fmt(p.z)));
            }
        }
        params.push_str(&format!(",{},{},{},{}", fmt(u0), fmt(u1), fmt(v0), fmt(v1)));
        self.push(Pending {
            kind: 128,
            form: 0,
            transform: None,
            independent: false,
            label: String::new(),
            params,
        })
    }

    /// The deck: start, global, directory, parameters, terminate.
    fn serialize(&self) -> String {
        let mut s = String::new();
        fn push_record(s: &mut String, body: &str, section: char, seq: usize) {
            s.push_str(&format!("{body:<72}{section}{seq:>7}\n"));
        }
        push_record(&mut s, "ogeom IGES writer", 'S', 1);

        // Parameter text per entity, then directory and parameter sections
        // interleaved by the format's mutual pointers.
        let mut param_lines: Vec<Vec<String>> = Vec::with_capacity(self.entities.len());
        for e in &self.entities {
            let full = format!("{},{};", e.kind, e.params);
            param_lines.push(wrap_params(&full));
        }
        let mut param_starts = Vec::with_capacity(self.entities.len());
        let mut next_param = 1usize;
        for lines in &param_lines {
            param_starts.push(next_param);
            next_param += lines.len();
        }

        let globals = [
            "1H,".to_string(),
            "1H;".to_string(),
            "5Hogeom".to_string(),
            "9Hmodel.igs".to_string(),
            "5Hogeom".to_string(),
            "5Hogeom".to_string(),
            "32".to_string(),
            "308".to_string(),
            "15".to_string(),
            "308".to_string(),
            "15".to_string(),
            "5Hogeom".to_string(),
            "1.".to_string(),
            // Millimetres, stated twice as the format wants.
            "2".to_string(),
            "2HMM".to_string(),
            "1".to_string(),
            "0.01".to_string(),
            "15H20260807.000000".to_string(),
            fmt(1e-7),
            "0.".to_string(),
            "5Hogeom".to_string(),
            "5Hogeom".to_string(),
            "11".to_string(),
            "0".to_string(),
            "15H20260807.000000".to_string(),
        ];
        let global_text = globals.join(",") + ";";
        let global_lines = wrap_params(&global_text);
        for (i, line) in global_lines.iter().enumerate() {
            push_record(&mut s, line, 'G', i + 1);
        }

        let mut d_seq = 1usize;
        for (i, e) in self.entities.iter().enumerate() {
            let transform_de = e.transform.map_or(0, |t| self.de(t));
            let status = if e.independent {
                "00000000"
            } else {
                "00010000"
            };
            let line1 = format!(
                "{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{status}",
                e.kind, param_starts[i], 0, 0, 0, 0, transform_de, 0
            );
            push_record(&mut s, &line1, 'D', d_seq);
            let line2 = format!(
                "{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
                e.kind,
                0,
                0,
                param_lines[i].len(),
                e.form,
                "",
                "",
                e.label,
                0
            );
            push_record(&mut s, &line2, 'D', d_seq + 1);
            d_seq += 2;
        }

        let mut p_seq = 1usize;
        for (i, lines) in param_lines.iter().enumerate() {
            let back = self.de(i);
            for line in lines {
                s.push_str(&format!("{line:<64}{back:>8}P{p_seq:>7}\n"));
                p_seq += 1;
            }
        }
        let tail = format!(
            "S{:>7}G{:>7}D{:>7}P{:>7}",
            1,
            global_lines.len(),
            d_seq - 1,
            p_seq - 1
        );
        push_record(&mut s, &tail, 'T', 1);
        s
    }
}

/// The same B-spline traversed the other way, over the same domain.
fn reversed_bspline(b: &BSplineCurve) -> OgeomResult<BSplineCurve> {
    match Curve::BSpline(b.clone()).reversed() {
        Curve::BSpline(back) => Ok(back),
        _ => ogeom_bail!(Construction, "a reversed B-spline is a B-spline"),
    }
}

/// A rigid transform quantized to bits, for per-placement deduplication,
/// the same probe the STEP writer uses.
fn transform_bits(t: &Transform) -> [u64; 3] {
    let p = t.apply(Point::new(0.123_456_789, 9.87, -3.21));
    [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
}

/// A real in the file, the shorter of its positional and exponent
/// spellings, each the shortest that reads back to the same value, with a
/// decimal point kept so a reader that types by spelling reads a real.
///
/// Positional notation alone spells a coefficient of 1.5e-51 in
/// sixty-nine characters, and a record holds sixty-four.
fn fmt(v: f64) -> String {
    let positional = {
        let s = format!("{v}");
        if s.contains('.') { s } else { format!("{s}.") }
    };
    let exponent = {
        let s = format!("{v:E}");
        match s.split_once('E') {
            Some((mantissa, power)) if !mantissa.contains('.') => {
                format!("{mantissa}.E{power}")
            }
            _ => s,
        }
    };
    if exponent.len() < positional.len() {
        exponent
    } else {
        positional
    }
}

/// Parameter text into records of at most 64 data columns, split between
/// parameters.
///
/// The reader joins a parameter record's data columns as they stand,
/// padding included, so a parameter is never split where padding would
/// land inside it: one that does not fit starts a new record. One longer
/// than a whole record (only a long name can be) fills every record it
/// crosses to exactly 64 columns, so the pieces rejoin with nothing
/// between them.
fn wrap_params(text: &str) -> Vec<String> {
    const COLUMNS: usize = 64;
    let mut lines = Vec::new();
    let mut current = String::new();
    for piece in text.split_inclusive(',') {
        if current.len() + piece.len() > COLUMNS && piece.len() <= COLUMNS {
            lines.push(std::mem::take(&mut current));
        }
        let mut rest = piece;
        while current.len() + rest.len() > COLUMNS {
            let mut room = COLUMNS - current.len();
            while !rest.is_char_boundary(room) {
                room -= 1;
            }
            current.push_str(&rest[..room]);
            rest = &rest[room..];
            lines.push(std::mem::take(&mut current));
        }
        current.push_str(rest);
    }
    if !current.is_empty() {
        lines.push(current);
    }
    lines
}

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code")]
mod tests {
    use super::{fmt, wrap_params};

    #[test]
    fn a_real_takes_its_shorter_spelling_and_reads_back() {
        for v in [
            0.0,
            1.0,
            -2.5,
            1e-51,
            1.48776941870366e-51,
            2.76e-76,
            1e20,
            123_456.789,
        ] {
            let s = fmt(v);
            assert!(s.len() <= 24, "{v} spelt {s}");
            assert!(s.contains('.'), "{s} reads as a real");
            let back: f64 = s.replace('D', "E").parse().unwrap();
            assert_eq!(back, v, "{s}");
        }
        assert_eq!(fmt(1.0), "1.");
        assert_eq!(fmt(1e-51), "1.E-51");
    }

    #[test]
    fn a_parameter_longer_than_a_record_fills_the_records_it_crosses() {
        let name = format!("100H{},", "x".repeat(100));
        let text = format!("128,{name}1.5,2.;");
        let lines = wrap_params(&text);
        assert!(lines.iter().all(|l| l.len() <= 64));
        // The pieces rejoin as written: every record the long name crosses
        // is full, so padding never lands inside it.
        let long: Vec<&String> = lines.iter().filter(|l| l.contains('x')).collect();
        for l in &long[..long.len() - 1] {
            assert_eq!(l.len(), 64);
        }
        assert_eq!(lines.concat(), text);
    }
}
