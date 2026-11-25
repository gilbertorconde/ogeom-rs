//! Product structure, names, colours and construction geometry.

use super::*;

impl<'a> Reader<'a> {
    // --- product structure, names, colours -----------------------------------

    /// Every `GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION` and
    /// `GEOMETRICALLY_BOUNDED_SURFACE_SHAPE_REPRESENTATION`, in file order,
    /// with the construction geometry it holds: the points, bounded curves
    /// and bounded surfaces of its geometric sets, and the axis placements
    /// that are items of its own. An element the reader cannot read is a
    /// warning, and the rest of the set still reads.
    pub(super) fn construction_sets(&mut self) -> OgeomResult<ConstructionSets> {
        let mut reps: Vec<u64> = self
            .instances
            .iter()
            .filter(|(_, inst)| {
                inst.part("GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION")
                    .is_some()
                    || inst
                        .part("GEOMETRICALLY_BOUNDED_SURFACE_SHAPE_REPRESENTATION")
                        .is_some()
            })
            .map(|(id, _)| id)
            .collect();
        reps.sort_unstable();
        let mut out = Vec::new();
        for rep in reps {
            let items: Vec<u64> = {
                let Some(inst) = self.instances.get(rep) else {
                    continue;
                };
                let args = inst
                    .part("GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION")
                    .or_else(|| inst.part("GEOMETRICALLY_BOUNDED_SURFACE_SHAPE_REPRESENTATION"))
                    .unwrap_or(&[]);
                args.get(1)
                    .and_then(Arg::list)
                    .unwrap_or(&[])
                    .iter()
                    .filter_map(Arg::reference)
                    .collect()
            };
            let mut read = Vec::new();
            for item in items {
                ogeom_core::progress::checkpoint()?;
                let Some(inst) = self.instances.get(item) else {
                    continue;
                };
                let set = inst
                    .part("GEOMETRIC_CURVE_SET")
                    .or_else(|| inst.part("GEOMETRIC_SET"));
                if let Some(args) = set {
                    let elements: Vec<u64> = args
                        .get(1)
                        .and_then(Arg::list)
                        .unwrap_or(&[])
                        .iter()
                        .filter_map(Arg::reference)
                        .collect();
                    self.visited.mark(item);
                    for element in elements {
                        if let Some(found) = self.construction_item(element) {
                            read.push(found);
                        }
                    }
                } else if inst.part("AXIS2_PLACEMENT_3D").is_some() {
                    let name = item_name(inst);
                    match self.frame(item) {
                        Ok(frame) => read.push((name, ConstructionGeometry::Axis(frame))),
                        Err(e) => self
                            .report
                            .warnings
                            .push(format!("#{item}: a construction axis is not read: {e}")),
                    }
                }
            }
            if !read.is_empty() {
                self.visited.mark(rep);
                out.push((rep, read));
            }
        }
        Ok(out)
    }

    /// One element of a geometric set, named: a Cartesian point, a bounded
    /// curve over its stretch (a trimmed curve's between its trims), or a
    /// bounded surface over its window (a rectangularly trimmed one's in
    /// its basis's parameters). `None`, with a warning, for what does not
    /// read or is not bounded.
    fn construction_item(&mut self, id: u64) -> Option<(String, ConstructionGeometry)> {
        let inst = self.instances.get(id)?;
        let name = item_name(inst);
        let keyword = inst.keyword().to_owned();
        let surface = inst.parts.iter().any(|(k, _)| {
            k == "PLANE"
                || k.ends_with("_SURFACE")
                || (k.starts_with("SURFACE_") && !k.ends_with("_CURVE"))
        });
        let bounded = |(lo, hi): (f64, f64)| {
            lo.is_finite() && hi.is_finite() && lo < hi && hi - lo <= SURFACE_EXTENT
        };
        let read: OgeomResult<Option<ConstructionGeometry>> =
            if inst.part("CARTESIAN_POINT").is_some() {
                self.point(id).map(|p| Some(ConstructionGeometry::Point(p)))
            } else if surface {
                self.construction_surface(id, &keyword).map(|found| {
                    found
                        .filter(|(_, u, v)| bounded(*u) && bounded(*v))
                        .map(|(surface, u, v)| ConstructionGeometry::Surface { surface, u, v })
                })
            } else if keyword == "TRIMMED_CURVE" {
                self.construction_trim(id).map(|found| {
                    found.map(|(curve, range)| ConstructionGeometry::Curve { curve, range })
                })
            } else {
                self.curve(id).map(|found| {
                    found
                        .map(|curve| (curve.domain(), curve))
                        .filter(|(range, _)| bounded(*range))
                        .map(|(range, curve)| ConstructionGeometry::Curve { curve, range })
                })
            };
        match read {
            Ok(Some(geometry)) => Some((name, geometry)),
            Ok(None) => {
                self.report.warnings.push(format!(
                    "#{id}: a geometric set's {keyword} is not read: it is not bounded, \
                     or not a kind read yet"
                ));
                None
            }
            Err(e) => {
                self.report.warnings.push(format!(
                    "#{id}: a geometric set's {keyword} is not read: {e}"
                ));
                None
            }
        }
    }

    /// A trimmed curve of a geometric set over its stretch: between its
    /// trims' parameters where the file names the parameters its master
    /// trim and the basis is one this kernel parameterizes as the file
    /// does (a line, a conic of its own frame, a spline), running forward;
    /// else as [`Self::bounded_curve`] reads it, between its trims' points.
    fn construction_trim(&mut self, id: u64) -> OgeomResult<Option<(Curve, (f64, f64))>> {
        let args = self.args(id, "TRIMMED_CURVE")?;
        let basis_id = args.get(1).and_then(Arg::reference).unwrap_or(0);
        let forward = !args.get(4).is_some_and(|a| a.is_enum("F"));
        let by_parameter = args.get(5).is_some_and(|a| a.is_enum("PARAMETER"));
        let parameter = |trim: Option<&Arg>| {
            trim.and_then(Arg::list)?
                .iter()
                .find_map(|item| match item {
                    Arg::Typed(typed) if typed.0 == "PARAMETER_VALUE" => {
                        typed.1.first().and_then(Arg::number)
                    }
                    _ => None,
                })
        };
        let (Some(v0), Some(v1)) = (parameter(args.get(2)), parameter(args.get(3))) else {
            return self.bounded_curve(id);
        };
        if !forward || !by_parameter {
            return self.bounded_curve(id);
        }
        let basis_keyword = self.instance(basis_id)?.keyword().to_owned();
        let scale = match basis_keyword.as_str() {
            "LINE" => {
                let line = self.args(basis_id, "LINE")?;
                let vector = self.args(ref_at(line, 2), "VECTOR")?;
                vector.get(2).and_then(Arg::number).unwrap_or(1.0) * self.report.scale_mm
            }
            "CIRCLE" | "ELLIPSE" => self.angle_scale,
            "B_SPLINE_CURVE_WITH_KNOTS" => 1.0,
            _ if self
                .instance(basis_id)?
                .part("B_SPLINE_CURVE_WITH_KNOTS")
                .is_some() =>
            {
                1.0
            }
            _ => return self.bounded_curve(id),
        };
        let Some(basis) = self.curve(basis_id)? else {
            return Ok(None);
        };
        let (a, mut b) = (v0 * scale, v1 * scale);
        if basis.is_periodic() && b <= a {
            let (lo, hi) = basis.domain();
            b += hi - lo;
        }
        if a >= b {
            return self.bounded_curve(id);
        }
        Ok(Some((basis, (a, b))))
    }

    /// A surface of a geometric set with its window: a rectangularly
    /// trimmed surface's basis over the trim's parameters, read in the
    /// file's units for the basis's kind; any other surface over its own
    /// domain.
    #[allow(clippy::type_complexity)]
    fn construction_surface(
        &mut self,
        id: u64,
        keyword: &str,
    ) -> OgeomResult<Option<(SurfaceGeometry, (f64, f64), (f64, f64))>> {
        if keyword != "RECTANGULAR_TRIMMED_SURFACE" {
            return Ok(self.surface(id)?.map(|surface| {
                let (u, v) = surface.domain();
                (surface, u, v)
            }));
        }
        let args = self.args(id, "RECTANGULAR_TRIMMED_SURFACE")?;
        let basis_id = args.get(1).and_then(Arg::reference).unwrap_or(0);
        let Some(basis) = self.surface(basis_id)? else {
            return Ok(None);
        };
        let (length, angle) = (self.report.scale_mm, self.angle_scale);
        // How the file measures each parameter of the basis: a plane's in
        // lengths, an angle round each axis of revolution, a spline's and a
        // swept curve's own parameter as it is.
        let (u_scale, v_scale) = match self.instance(basis_id)?.keyword() {
            "PLANE" => (length, length),
            "CYLINDRICAL_SURFACE" | "CONICAL_SURFACE" => (angle, length),
            "SPHERICAL_SURFACE" | "TOROIDAL_SURFACE" | "DEGENERATE_TOROIDAL_SURFACE" => {
                (angle, angle)
            }
            "SURFACE_OF_LINEAR_EXTRUSION" => (1.0, length),
            "SURFACE_OF_REVOLUTION" => (1.0, angle),
            _ => (1.0, 1.0),
        };
        let n = |i: usize| args.get(i).and_then(Arg::number).unwrap_or(0.0);
        let (u1, u2) = (n(2) * u_scale, n(3) * u_scale);
        let (v1, v2) = (n(4) * v_scale, n(5) * v_scale);
        Ok(Some((
            basis,
            (u1.min(u2), u1.max(u2)),
            (v1.min(v2), v1.max(v2)),
        )))
    }

    /// Construction items into the document under `product`; one the
    /// document refuses is a warning.
    fn add_construction(
        &mut self,
        document: &mut ogeom_doc::Document,
        product: ogeom_doc::ProductId,
        items: Vec<(String, ConstructionGeometry)>,
    ) {
        for (name, geometry) in items {
            if let Err(e) = document.add_construction(product, &name, geometry) {
                self.report
                    .warnings
                    .push(format!("construction item '{name}' is not kept: {e}"));
            }
        }
    }

    /// Assemble the document: the model, plus everything the file says about
    /// products, assemblies, placements and appearance.
    ///
    /// Takes the model out of the reader; geometry reading is over by the
    /// time structure is read. Structure that resists becomes a warning and a
    /// flat document, never an error: the geometry is already good, and a
    /// mangled product tree should not take it down.
    pub(super) fn document(
        &mut self,
        by_item: &FastMap<u64, Shape>,
        solids: &[Shape],
        shells: &[Shape],
        construction: ConstructionSets,
    ) -> OgeomResult<ogeom_doc::Document> {
        // The graph is walked before the model moves, because frames scale
        // through the reader's own unit handling.
        let mut sets: FastMap<u64, Vec<(String, ConstructionGeometry)>> = FastMap::default();
        let order: Vec<u64> = construction.iter().map(|(rep, _)| *rep).collect();
        sets.extend(construction);
        let structure = self.product_structure(by_item, &mut sets);
        let colours = self.colours(by_item);
        let pmi = self.pmi_of();

        let mut document = ogeom_doc::Document::over(std::mem::take(&mut self.model));
        match structure {
            Some(products) => self.build_products(&mut document, products),
            None => {
                for (i, solid) in solids.iter().enumerate() {
                    document.add_part(format!("solid-{i}"), solid.clone());
                }
                for (i, shell) in shells.iter().enumerate() {
                    document.add_part(format!("shell-{i}"), shell.clone());
                }
            }
        }
        // A set no product reaches is a part of its own, named as its
        // representation is.
        for (i, rep) in order.into_iter().enumerate() {
            let Some(items) = sets.remove(&rep) else {
                continue;
            };
            let name = self
                .instances
                .get(rep)
                .and_then(|inst| inst.parts.first())
                .and_then(|(_, args)| match args.first() {
                    Some(Arg::Str(name)) if !name.is_empty() => Some(name.clone()),
                    _ => None,
                })
                .unwrap_or_else(|| format!("construction-{i}"));
            let empty = ogeom_algo::build::make_compound(document.model_mut(), &[])?.shape;
            let part = document.add_part(name, empty);
            self.add_construction(&mut document, part, items);
        }
        for (shape, colour) in colours {
            document.set_colour(&shape, colour);
        }
        *document.pmi_mut() = pmi;
        for view in self.views() {
            document.add_view(view);
        }
        Ok(document)
    }

    /// The file's product graph, or `None` when it has none worth the name.
    fn product_structure(
        &mut self,
        by_item: &FastMap<u64, Shape>,
        sets: &mut FastMap<u64, Vec<(String, ConstructionGeometry)>>,
    ) -> Option<Vec<PdEntry>> {
        // PRODUCT_DEFINITION -> name, via formation and product.
        let mut pds: Vec<u64> = self
            .instances
            .iter()
            .filter(|(_, inst)| inst.part("PRODUCT_DEFINITION").is_some())
            .filter(|(_, inst)| inst.part("PRODUCT_DEFINITION_RELATIONSHIP").is_none())
            .map(|(id, _)| id)
            .collect();
        pds.sort_unstable();
        if pds.is_empty() {
            return None;
        }

        // SHAPE_DEFINITION_REPRESENTATION: definition (a PRODUCT_DEFINITION_SHAPE
        // over a PD or a usage) -> shape representation.
        let mut sr_of_pd: FastMap<u64, u64> = FastMap::default();
        let sdrs: Vec<(u64, u64)> = self
            .instances
            .iter()
            .filter(|(_, inst)| inst.part("SHAPE_DEFINITION_REPRESENTATION").is_some())
            .filter_map(|(id, _)| {
                let args = self.args(id, "SHAPE_DEFINITION_REPRESENTATION").ok()?;
                Some((args.first()?.reference()?, args.get(1)?.reference()?))
            })
            .collect();
        for (pds_id, sr) in sdrs {
            if let Some(definition) = self.definition_of_shape(pds_id) {
                sr_of_pd.insert(definition, sr);
            }
        }

        // A product's solids may live one representation over: AP203 files
        // routinely tie the product to a bare axis representation and hang
        // the B-rep off it through a plain SHAPE_REPRESENTATION_RELATIONSHIP.
        // Only the plain ones are followed; the transformation-carrying kind
        // is an assembly edge, and following it would leak one product's
        // geometry into another.
        let mut linked: FastMap<u64, Vec<u64>> = FastMap::default();
        let mut srrs: Vec<u64> = self
            .instances
            .iter()
            .filter(|(_, inst)| {
                inst.part("SHAPE_REPRESENTATION_RELATIONSHIP").is_some()
                    && inst
                        .part("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION")
                        .is_none()
            })
            .map(|(id, _)| id)
            .collect();
        srrs.sort_unstable();
        for srr in srrs {
            let args = {
                let Ok(instance) = self.instance(srr) else {
                    continue;
                };
                let Some(args) = instance
                    .part("SHAPE_REPRESENTATION_RELATIONSHIP")
                    .or_else(|| instance.part("REPRESENTATION_RELATIONSHIP"))
                else {
                    continue;
                };
                args.to_vec()
            };
            if let (Some(a), Some(b)) = (
                args.get(2).and_then(Arg::reference),
                args.get(3).and_then(Arg::reference),
            ) {
                linked.entry(a).or_default().push(b);
                linked.entry(b).or_default().push(a);
            }
        }

        let mut entries = Vec::new();
        for pd in pds {
            let name = self.product_name(pd).unwrap_or_else(|| format!("#{pd}"));
            let mut shapes: Vec<Shape> = Vec::new();
            let mut construction = Vec::new();
            if let Some(&sr) = sr_of_pd.get(&pd) {
                let mut reps = vec![sr];
                reps.extend(linked.get(&sr).into_iter().flatten().copied());
                for rep in reps {
                    if let Some(items) = self.representation_items(rep) {
                        shapes.extend(items.iter().filter_map(|item| by_item.get(item).cloned()));
                    }
                    construction.extend(sets.remove(&rep).into_iter().flatten());
                }
            }
            entries.push(PdEntry {
                pd,
                name,
                shapes,
                construction,
                children: Vec::new(),
            });
        }

        // NEXT_ASSEMBLY_USAGE_OCCURRENCE: parent -> child, with the placement
        // recovered from the CONTEXT_DEPENDENT_SHAPE_REPRESENTATION over it.
        let mut nauos: Vec<u64> = self
            .instances
            .iter()
            .filter(|(_, inst)| inst.part("NEXT_ASSEMBLY_USAGE_OCCURRENCE").is_some())
            .map(|(id, _)| id)
            .collect();
        nauos.sort_unstable();
        let entry_of_pd: FastMap<u64, usize> =
            entries.iter().enumerate().map(|(i, e)| (e.pd, i)).collect();
        for nauo in nauos {
            let Ok(args) = self.args(nauo, "NEXT_ASSEMBLY_USAGE_OCCURRENCE") else {
                continue;
            };
            let (Some(parent), Some(child)) = (
                args.get(3).and_then(Arg::reference),
                args.get(4).and_then(Arg::reference),
            ) else {
                continue;
            };
            let name = args
                .get(5)
                .and_then(|a| match a {
                    Arg::Str(s) if !s.is_empty() => Some(s.clone()),
                    _ => None,
                })
                .or_else(|| {
                    self.args(nauo, "NEXT_ASSEMBLY_USAGE_OCCURRENCE")
                        .ok()
                        .and_then(|args| match args.get(1) {
                            Some(Arg::Str(s)) if !s.is_empty() => Some(s.clone()),
                            _ => None,
                        })
                });
            let child_sr = sr_of_pd.get(&child).copied();
            let at = self
                .usage_transform(nauo, child_sr)
                .unwrap_or(Transform::IDENTITY);
            if let Some(&at_index) = entry_of_pd.get(&parent) {
                entries[at_index].children.push((child, at, name));
            }
        }
        Some(entries)
    }

    /// The `PRODUCT_DEFINITION` (or usage) a `PRODUCT_DEFINITION_SHAPE` is over.
    fn definition_of_shape(&mut self, pds_id: u64) -> Option<u64> {
        let args = self.args(pds_id, "PRODUCT_DEFINITION_SHAPE").ok()?;
        args.get(2).and_then(Arg::reference)
    }

    /// A product definition's product name.
    fn product_name(&mut self, pd: u64) -> Option<String> {
        let formation = self
            .args(pd, "PRODUCT_DEFINITION")
            .ok()?
            .get(2)
            .and_then(Arg::reference)?;
        // The formation, plain or with its source named (a mesh converter's
        // habit, and a modeller's for an assembly's parts) carries the
        // product in its third slot either way, whether it stands alone or
        // as one part of a complex instance.
        let product = [
            "PRODUCT_DEFINITION_FORMATION",
            "PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE",
        ]
        .into_iter()
        .find_map(|keyword| {
            self.instance(formation)
                .ok()?
                .part(keyword)
                .or_else(|| self.args(formation, keyword).ok())?
                .get(2)
                .and_then(Arg::reference)
        })?;
        // The product's name, or its id where the name is blank: a mesh
        // converter fills the id and leaves the name empty.
        let args = self.args(product, "PRODUCT").ok()?;
        [args.get(1), args.first()]
            .into_iter()
            .find_map(|arg| match arg {
                Some(Arg::Str(name)) if !name.is_empty() => Some(name.clone()),
                _ => None,
            })
    }

    /// A representation's item references.
    pub(super) fn representation_items(&mut self, sr: u64) -> Option<Vec<u64>> {
        let instance = self.instance(sr).ok()?;
        // Every subtype a file names a shape representation by: the B-rep
        // kinds, the surface-model kinds a mesh converter writes, and the
        // plain one.
        let args = [
            "SHAPE_REPRESENTATION",
            "ADVANCED_BREP_SHAPE_REPRESENTATION",
            "MANIFOLD_SURFACE_SHAPE_REPRESENTATION",
            "FACETED_BREP_SHAPE_REPRESENTATION",
            "GEOMETRICALLY_BOUNDED_SURFACE_SHAPE_REPRESENTATION",
            "GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION",
            "REPRESENTATION",
        ]
        .into_iter()
        .find_map(|keyword| instance.part(keyword))?
        .to_vec();
        Some(
            args.get(1)?
                .list()?
                .iter()
                .filter_map(Arg::reference)
                .collect(),
        )
    }

    /// The placement a usage's `CONTEXT_DEPENDENT_SHAPE_REPRESENTATION` states.
    ///
    /// The transformation aligns an axis placement in the child's space with
    /// one in the parent's; which item is which follows from which side of
    /// the representation relationship is the child's own shape
    /// representation, not from argument order, because real files disagree
    /// about the order.
    fn usage_transform(&mut self, nauo: u64, child_sr: Option<u64>) -> Option<Transform> {
        if self.cdsr_of_nauo.is_none() {
            // One pass over the CDSRs, each resolved to the usage it
            // describes, in ascending id order, so a usage described twice
            // keeps the lowest-id one.
            let mut cdsrs: Vec<u64> = self
                .instances
                .iter()
                .filter(|(_, inst)| {
                    inst.part("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION")
                        .is_some()
                })
                .map(|(id, _)| id)
                .collect();
            cdsrs.sort_unstable();
            let mut index: FastMap<u64, u64> = FastMap::default();
            for id in cdsrs {
                let Some(owner) = self
                    .args(id, "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION")
                    .ok()
                    .and_then(|args| args.get(1).and_then(Arg::reference))
                    .and_then(|pds| self.definition_of_shape(pds))
                else {
                    continue;
                };
                index.entry(owner).or_insert(id);
            }
            self.cdsr_of_nauo = Some(index);
        }
        let cdsr = *self.cdsr_of_nauo.as_ref()?.get(&nauo)?;
        let rr = self
            .args(cdsr, "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION")
            .ok()?
            .first()
            .and_then(Arg::reference)?;
        let (rep_1, rep_2) = {
            let args = self.args(rr, "REPRESENTATION_RELATIONSHIP").ok()?;
            (
                args.get(2).and_then(Arg::reference)?,
                args.get(3).and_then(Arg::reference)?,
            )
        };
        let idt = self
            .args(rr, "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION")
            .ok()?
            .first()
            .and_then(Arg::reference)?;
        let (item_1, item_2) = {
            let args = self.args(idt, "ITEM_DEFINED_TRANSFORMATION").ok()?;
            (
                args.get(2).and_then(Arg::reference)?,
                args.get(3).and_then(Arg::reference)?,
            )
        };
        let frame_1 = self.frame(item_1).ok()?;
        let frame_2 = self.frame(item_2).ok()?;
        // item_1 pairs with rep_1. When rep_1 is the child's representation,
        // the child's frame_1 lands on the parent's frame_2.
        let child_first = match child_sr {
            Some(sr) => rep_1 == sr || rep_2 != sr,
            None => true,
        };
        Some(if child_first {
            Transform::from_frame(&frame_2) * Transform::to_frame(&frame_1)
        } else {
            Transform::from_frame(&frame_1) * Transform::to_frame(&frame_2)
        })
    }

    /// Products into the document: assemblies for the parents, parts for the
    /// shaped, instances for the usage edges.
    fn build_products(&mut self, document: &mut ogeom_doc::Document, mut entries: Vec<PdEntry>) {
        let mut ids: FastMap<u64, ogeom_doc::ProductId> = FastMap::default();
        for entry in &mut entries {
            let shape = match entry.shapes.len() {
                // A product whose content is all construction geometry is
                // a part over nothing, carrying it.
                0 if !entry.construction.is_empty() && entry.children.is_empty() => {
                    match ogeom_algo::build::make_compound(document.model_mut(), &[]) {
                        Ok(built) => Some(built.shape),
                        Err(_) => None,
                    }
                }
                0 => None,
                1 => Some(entry.shapes[0].clone()),
                _ => match ogeom_algo::build::make_compound(document.model_mut(), &entry.shapes) {
                    Ok(built) => Some(built.shape),
                    Err(_) => Some(entry.shapes[0].clone()),
                },
            };
            if entry.children.is_empty() {
                if let Some(shape) = shape {
                    ids.insert(entry.pd, document.add_part(&entry.name, shape));
                }
                // A product with neither shape nor children holds nothing a
                // document can say; it is left out.
            } else {
                let assembly = document.add_assembly(&entry.name);
                ids.insert(entry.pd, assembly);
                // An assembly with its own geometry keeps it as a body part
                // placed at identity: rare, but files do it.
                if let Some(shape) = shape {
                    let body = document.add_part(format!("{}-body", entry.name), shape);
                    let _ = document.add_instance(assembly, body, Transform::IDENTITY, None);
                }
            }
            if let Some(&id) = ids.get(&entry.pd) {
                let items = std::mem::take(&mut entry.construction);
                self.add_construction(document, id, items);
            }
        }
        for entry in &entries {
            let Some(&parent) = ids.get(&entry.pd) else {
                continue;
            };
            for (child, at, name) in &entry.children {
                let Some(&child_id) = ids.get(child) else {
                    self.report.warnings.push(format!(
                        "#{child}: an assembly child holds nothing readable"
                    ));
                    continue;
                };
                if let Err(e) = document.add_instance(parent, child_id, *at, name.clone()) {
                    self.report
                        .warnings
                        .push(format!("assembly edge #{} -> #{child}: {e}", entry.pd));
                }
            }
        }
    }

    /// Colours from styled items, keyed to the shapes they style.
    fn colours(&mut self, by_item: &FastMap<u64, Shape>) -> Vec<(Shape, ogeom_doc::Colour)> {
        let mut styled: Vec<u64> = self
            .instances
            .iter()
            .filter(|(_, inst)| {
                inst.part("STYLED_ITEM").is_some() || inst.part("OVER_RIDING_STYLED_ITEM").is_some()
            })
            .map(|(id, _)| id)
            .collect();
        // Sorted so an item styled twice resolves the same way every run;
        // overriding styles carry higher instance numbers in practice, and a
        // later same-shape entry wins the map insertion downstream.
        styled.sort_unstable();
        let mut out = Vec::new();
        for id in styled {
            let args = {
                let Ok(instance) = self.instance(id) else {
                    continue;
                };
                let Some(args) = instance
                    .part("STYLED_ITEM")
                    .or_else(|| instance.part("OVER_RIDING_STYLED_ITEM"))
                else {
                    continue;
                };
                args.to_vec()
            };
            let Some(item) = args.get(2).and_then(Arg::reference) else {
                continue;
            };
            let shapes = match by_item.get(&item).or_else(|| self.faces.get(&item)) {
                Some(shape) => vec![shape.clone()],
                None => self.edge_shapes(item),
            };
            if shapes.is_empty() {
                continue;
            }
            let styles: Vec<u64> = args
                .get(1)
                .and_then(Arg::list)
                .map(|list| list.iter().filter_map(Arg::reference).collect())
                .unwrap_or_default();
            if let Some(colour) = styles.iter().find_map(|&style| self.colour_in(style, 0)) {
                out.extend(shapes.into_iter().map(|shape| (shape, colour)));
            }
        }
        out
    }

    /// The first `COLOUR_RGB` reachable from a presentation style, depth-bounded.
    ///
    /// The styled-item chain has five links and real files rearrange them, so
    /// the walk follows references rather than the textbook path.
    fn colour_in(&mut self, id: u64, depth: usize) -> Option<ogeom_doc::Colour> {
        if depth > 6 {
            return None;
        }
        let (keyword, args) = {
            let instance = self.instance(id).ok()?;
            let all: Vec<Arg> = instance
                .parts()
                .flat_map(|(_, args)| args.iter().cloned())
                .collect();
            (instance.keyword().to_owned(), all)
        };
        if keyword == "COLOUR_RGB" {
            let channel = |i: usize| args.get(i).and_then(Arg::number);
            return Some(ogeom_doc::Colour::rgb(
                channel(1)?,
                channel(2)?,
                channel(3)?,
            ));
        }
        // A colour named rather than mixed: one of the eight the standard
        // (ISO 10303-46) defines, each channel fully on or off, read in any
        // case and with any spaces around it. Any other name is no colour,
        // said once for the entity.
        if keyword == "DRAUGHTING_PRE_DEFINED_COLOUR" {
            let name = match args.first() {
                Some(Arg::Str(name)) => name.as_str(),
                _ => "",
            };
            let channels = match name.trim().to_ascii_lowercase().as_str() {
                "black" => Some((0.0, 0.0, 0.0)),
                "red" => Some((1.0, 0.0, 0.0)),
                "green" => Some((0.0, 1.0, 0.0)),
                "blue" => Some((0.0, 0.0, 1.0)),
                "yellow" => Some((1.0, 1.0, 0.0)),
                "magenta" => Some((1.0, 0.0, 1.0)),
                "cyan" => Some((0.0, 1.0, 1.0)),
                "white" => Some((1.0, 1.0, 1.0)),
                _ => None,
            };
            if channels.is_none() && self.unnamed_colours.insert(id) {
                self.report.warnings.push(format!(
                    "#{id}: DRAUGHTING_PRE_DEFINED_COLOUR '{name}' names none of the eight \
                     predefined colours; the item is left uncoloured"
                ));
            }
            return channels.map(|(r, g, b)| ogeom_doc::Colour::rgb(r, g, b));
        }
        let mut refs: Vec<u64> = Vec::new();
        collect_refs(&args, &mut refs);
        refs.into_iter()
            .find_map(|next| self.colour_in(next, depth + 1))
    }
}

/// A product definition gathered from the file: its shapes and its usage
/// edges, before anything is committed to the document.
struct PdEntry {
    pd: u64,
    name: String,
    shapes: Vec<Shape>,
    /// The construction geometry its representations carry.
    construction: Vec<(String, ConstructionGeometry)>,
    children: Vec<(u64, Transform, Option<String>)>,
}
