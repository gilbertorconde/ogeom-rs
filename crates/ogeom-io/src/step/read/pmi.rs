//! Semantic PMI: dimensions, tolerances, datums, views and callouts.

use super::*;

impl<'a> Reader<'a> {
    // --- semantic PMI --------------------------------------------------------

    /// The file's semantic PMI: dimensions, geometric tolerances, datums.
    ///
    /// Annotations that resist stay out with a warning; PMI never takes the
    /// geometry down.
    pub(super) fn pmi_of(&mut self) -> ogeom_doc::Pmi {
        let mut pmi = ogeom_doc::Pmi::new();
        // Which STEP instance each annotation came from, so the presentation
        // pass can name the annotation a callout draws exactly rather than by
        // matching a string two annotations may share.
        let mut annotation_ids: FastMap<u64, ogeom_doc::Annotated> = FastMap::default();

        // Which topology each shape aspect describes: directly through
        // GEOMETRIC_ITEM_SPECIFIC_USAGE, and one relationship step outward,
        // because composite aspects hold their pieces through relationships.
        let mut aspect_items: FastMap<u64, Vec<ogeom_topo::TShapeId>> = FastMap::default();
        let mut gisus = self.ids_with("GEOMETRIC_ITEM_SPECIFIC_USAGE");
        gisus.sort_unstable();
        for id in gisus {
            let Ok(args) = self.args(id, "GEOMETRIC_ITEM_SPECIFIC_USAGE") else {
                continue;
            };
            let (Some(aspect), Some(item)) = (
                args.get(2).and_then(Arg::reference),
                args.get(4).and_then(Arg::reference),
            ) else {
                continue;
            };
            let nodes = match self.faces.get(&item) {
                Some(face) => vec![face.node()],
                None => self.edge_shapes(item).iter().map(Shape::node).collect(),
            };
            if !nodes.is_empty() {
                aspect_items.entry(aspect).or_default().extend(nodes);
            }
        }
        let mut adjacency: FastMap<u64, Vec<u64>> = FastMap::default();
        for id in self.ids_with("SHAPE_ASPECT_RELATIONSHIP") {
            let Ok(args) = self.args(id, "SHAPE_ASPECT_RELATIONSHIP") else {
                continue;
            };
            if let (Some(a), Some(b)) = (
                args.get(2).and_then(Arg::reference),
                args.get(3).and_then(Arg::reference),
            ) {
                adjacency.entry(a).or_default().push(b);
                adjacency.entry(b).or_default().push(a);
            }
        }
        let items_for = |aspect: u64| -> Vec<ogeom_topo::TShapeId> {
            // Three relationship steps: a composite aspect holds components,
            // a derived aspect sits behind a composite, and a datum one link
            // behind its features: the deepest chain real files exhibit.
            let mut reach = vec![aspect];
            for _ in 0..3 {
                let mut next = reach.clone();
                for a in &reach {
                    next.extend(adjacency.get(a).into_iter().flatten().copied());
                }
                next.sort_unstable();
                next.dedup();
                reach = next;
            }
            let mut out: Vec<ogeom_topo::TShapeId> = Vec::new();
            for a in reach {
                out.extend(aspect_items.get(&a).into_iter().flatten().copied());
            }
            out.sort_unstable();
            out.dedup();
            out
        };

        // Dimensions: characteristic -> representation, values from the
        // measure items, bounds from any plus/minus tolerance over the same
        // characteristic.
        let mut plus_minus: FastMap<u64, (Option<f64>, Option<f64>)> = FastMap::default();
        for id in self.ids_with("PLUS_MINUS_TOLERANCE") {
            let Ok(args) = self.args(id, "PLUS_MINUS_TOLERANCE") else {
                continue;
            };
            let (Some(tv), Some(dim)) = (
                args.first().and_then(Arg::reference),
                args.get(1).and_then(Arg::reference),
            ) else {
                continue;
            };
            let Ok(tv_args) = self.args(tv, "TOLERANCE_VALUE") else {
                continue;
            };
            let lower = tv_args
                .first()
                .and_then(Arg::reference)
                .and_then(|r| self.measure_value(r))
                .map(|(v, _)| v);
            let upper = tv_args
                .get(1)
                .and_then(Arg::reference)
                .and_then(|r| self.measure_value(r))
                .map(|(v, _)| v);
            plus_minus.insert(dim, (lower, upper));
        }
        let mut dcrs = self.ids_with("DIMENSIONAL_CHARACTERISTIC_REPRESENTATION");
        dcrs.sort_unstable();
        for id in dcrs {
            let Ok(args) = self.args(id, "DIMENSIONAL_CHARACTERISTIC_REPRESENTATION") else {
                continue;
            };
            let (Some(dim), Some(sdr)) = (
                args.first().and_then(Arg::reference),
                args.get(1).and_then(Arg::reference),
            ) else {
                continue;
            };
            let mut values = Vec::new();
            let mut kind = ogeom_doc::MeasureKind::Length;
            if let Ok(sdr_args) = self.args(sdr, "SHAPE_DIMENSION_REPRESENTATION") {
                for item in sdr_args.get(1).and_then(Arg::list).unwrap_or(&[]) {
                    if let Some(r) = item.reference()
                        && let Some((v, k)) = self.measure_value(r)
                    {
                        values.push(v);
                        kind = k;
                    }
                }
            }
            let (name, location, aspects) = self.dimension_shape(dim);
            let features: Vec<Vec<ogeom_topo::TShapeId>> =
                aspects.into_iter().map(&items_for).collect();
            let (minus, plus) = plus_minus.get(&dim).copied().unwrap_or((None, None));
            annotation_ids.insert(dim, ogeom_doc::Annotated::Dimension(pmi.dimensions.len()));
            pmi.dimensions.push(ogeom_doc::Dimension {
                name,
                values,
                kind,
                plus,
                minus,
                features,
                location,
            });
        }

        // Geometric tolerances: any instance whose subtype names one. The
        // complex form keeps its attributes on the GEOMETRIC_TOLERANCE part;
        // the simple form flattens them into the subtype's own list.
        let is_subtype = |k: &str| {
            k.ends_with("_TOLERANCE")
                && !matches!(
                    k,
                    "GEOMETRIC_TOLERANCE"
                        | "GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE"
                        | "GEOMETRIC_TOLERANCE_WITH_DEFINED_UNIT"
                        | "GEOMETRIC_TOLERANCE_WITH_MODIFIERS"
                        | "GEOMETRIC_TOLERANCE_WITH_MAXIMUM_TOLERANCE"
                        | "PLUS_MINUS_TOLERANCE"
                        | "TOLERANCE_VALUE"
                )
        };
        let mut gts: Vec<u64> = self
            .instances
            .iter()
            .filter(|(_, inst)| inst.parts().any(|(k, _)| is_subtype(k)))
            .map(|(id, _)| id)
            .collect();
        gts.sort_unstable();
        for id in gts {
            let (name, magnitude_ref, aspect, kind, modifiers, datum_refs) = {
                let Ok(instance) = self.instance(id) else {
                    continue;
                };
                let subtype = instance
                    .parts()
                    .map(|(k, _)| k.to_owned())
                    .find(|k| is_subtype(k));
                let Some(subtype) = subtype else {
                    continue;
                };
                let base = instance
                    .part("GEOMETRIC_TOLERANCE")
                    .or_else(|| instance.part(&subtype))
                    .unwrap_or(&[])
                    .to_vec();
                let name = match base.first() {
                    Some(Arg::Str(s)) => s.clone(),
                    _ => String::new(),
                };
                let magnitude_ref = base.get(2).and_then(Arg::reference);
                let aspect = base.get(3).and_then(Arg::reference);
                let kind = Some(subtype.trim_end_matches("_TOLERANCE").to_lowercase());
                // Modifiers ride on their own part in the complex form, as
                // a list of enumeration words.
                let modifiers: Vec<String> = instance
                    .part("GEOMETRIC_TOLERANCE_WITH_MODIFIERS")
                    .and_then(|args| args.first())
                    .and_then(Arg::list)
                    .map(|list| {
                        list.iter()
                            .filter_map(|arg| match arg {
                                Arg::Enum(word) => Some(word.to_lowercase()),
                                _ => None,
                            })
                            .collect()
                    })
                    .unwrap_or_default();
                // The complex form keeps the datum list on its own part;
                // the simple form appends it as the subtype's fifth argument.
                let datum_refs: Vec<u64> = instance
                    .part("GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE")
                    .and_then(|args| args.first())
                    .and_then(Arg::list)
                    .or_else(|| base.get(4).and_then(Arg::list))
                    .map(|list| list.iter().filter_map(Arg::reference).collect())
                    .unwrap_or_default();
                (name, magnitude_ref, aspect, kind, modifiers, datum_refs)
            };
            let Some(kind) = kind else {
                continue;
            };
            let magnitude = magnitude_ref
                .and_then(|r| self.measure_value(r))
                .map_or(0.0, |(v, _)| v);
            // A datum system states the whole frame, in precedence order
            // and with each reference's modifiers; any other reference is
            // one datum, or one common datum.
            let mut datums: Vec<ogeom_doc::DatumReference> = Vec::new();
            let mut datum_system = false;
            for &r in &datum_refs {
                if let Some(frame) = self.datum_system(r) {
                    datums.extend(frame);
                    datum_system = true;
                } else if let Some(letter) = self.datum_letter(r, 0) {
                    datums.push(letter.into());
                }
            }
            let items = aspect.map(items_for).unwrap_or_default();
            annotation_ids.insert(id, ogeom_doc::Annotated::Tolerance(pmi.tolerances.len()));
            pmi.tolerances.push(ogeom_doc::GeometricTolerance {
                kind,
                name,
                magnitude,
                modifiers,
                datums,
                datum_system,
                items,
            });
        }

        // Datums: the letters, with the features they mark reached through
        // the aspect graph.
        let mut datums = self.ids_with("DATUM");
        datums.sort_unstable();
        for id in datums {
            let letter = {
                let Ok(instance) = self.instance(id) else {
                    continue;
                };
                if instance.part("DATUM_FEATURE").is_some()
                    || instance.part("DATUM_REFERENCE").is_some()
                    || instance.part("DATUM_REFERENCE_COMPARTMENT").is_some()
                    || instance.part("DATUM_SYSTEM").is_some()
                {
                    continue;
                }
                let Some(args) = instance.part("DATUM") else {
                    continue;
                };
                match args.get(4) {
                    Some(Arg::Str(s)) if !s.is_empty() => s.clone(),
                    _ => continue,
                }
            };
            annotation_ids.insert(id, ogeom_doc::Annotated::Datum(pmi.datums.len()));
            pmi.datums.push(ogeom_doc::Datum {
                label: letter,
                items: items_for(id),
            });
        }

        // Datum targets: the pads a datum is actually established at. The
        // target's identifier is the letter's number (`A1` is target 1 of
        // datum A), and its placement and size come through the shape
        // representation the feature is associated with.
        let mut targets = self.ids_with("PLACED_DATUM_TARGET_FEATURE");
        targets.sort_unstable();
        for id in targets {
            if let Some(target) = self.datum_target(id, &items_for) {
                pmi.targets.push(target);
            }
        }

        // Presentation: what a viewer draws. A callout holds tessellated
        // annotation occurrences, each an indexed set of polylines over a
        // coordinates list; an annotation plane says which plane they are
        // drawn in and which callouts it holds; and a model item association
        // says which semantic annotation a callout is the picture of.
        let (callouts, callout_index) = self.callouts(&annotation_ids, &aspect_items);
        pmi.callouts = callouts;
        self.callout_index = callout_index;
        pmi
    }

    /// Saved views: every *named* draughting model is one; the unnamed one
    /// is the annotation-plane container this writer emits itself. The
    /// camera item gives the frame; the callout items give the subset.
    pub(super) fn views(&mut self) -> Vec<ogeom_doc::View> {
        let mut out = Vec::new();
        for id in self.ids_with("DRAUGHTING_MODEL") {
            let Ok(args) = self.args(id, "DRAUGHTING_MODEL") else {
                continue;
            };
            let name = match args.first() {
                Some(Arg::Str(text)) => text.clone(),
                _ => String::new(),
            };
            if name.is_empty() {
                continue;
            }
            let mut frame = None;
            let mut callouts = Vec::new();
            for item in args.get(1).and_then(Arg::list).unwrap_or(&[]) {
                let Some(item) = item.reference() else {
                    continue;
                };
                if let Ok(cam) = self.args(item, "CAMERA_MODEL_D3") {
                    if let Some(placement) = cam.get(1).and_then(Arg::reference) {
                        frame = self.frame(placement).ok();
                    }
                    continue;
                }
                if let Some(&index) = self.callout_index.get(&item) {
                    callouts.push(index);
                }
            }
            let Some(frame) = frame else { continue };
            out.push(ogeom_doc::View {
                name,
                frame,
                clipping: None,
                callouts,
            });
        }
        out
    }

    /// One placed datum target: which datum, which number, where and how big.
    fn datum_target(
        &mut self,
        id: u64,
        items_for: &dyn Fn(u64) -> Vec<ogeom_topo::TShapeId>,
    ) -> Option<ogeom_doc::DatumTarget> {
        let (target_id, description) = {
            let instance = self.instance(id).ok()?;
            let args = instance.part("PLACED_DATUM_TARGET_FEATURE")?;
            let target = match args.get(4) {
                Some(Arg::Str(s)) if !s.is_empty() => s.clone(),
                _ => return None,
            };
            let description = match args.get(1) {
                Some(Arg::Str(s)) => s.to_ascii_lowercase(),
                _ => String::new(),
            };
            (target, description)
        };
        // `A1`: the letter is the datum, the digits its number.
        let split = target_id
            .find(|c: char| c.is_ascii_digit())
            .unwrap_or(target_id.len());
        let (letter, number) = target_id.split_at(split);
        let index = number.parse::<u32>().unwrap_or(1);
        let datum = if letter.is_empty() {
            self.datum_letter(id, 0).unwrap_or_default()
        } else {
            letter.to_owned()
        };

        // The target's placement and size live in a shape representation the
        // target's own property definition names. The lengths come in the
        // file's own unit, as every length does. Both hops go through the
        // indexes rather than a scan of every entity per property *per
        // target*.
        self.ensure_property_indexes();
        let mut frame = None;
        let mut lengths: Vec<f64> = Vec::new();
        let properties: Vec<u64> = self
            .properties_of_definition
            .as_ref()
            .and_then(|m| m.get(&id).cloned())
            .unwrap_or_default();
        for property in properties {
            let sdrs: Vec<u64> = self
                .sdrs_of_property
                .as_ref()
                .and_then(|m| m.get(&property).cloned())
                .unwrap_or_default();
            for sdr in sdrs {
                let Ok(args) = self.args(sdr, "SHAPE_DEFINITION_REPRESENTATION") else {
                    continue;
                };
                let Some(rep) = args.get(1).and_then(Arg::reference) else {
                    continue;
                };
                for item in self.representation_items(rep).unwrap_or_default() {
                    let keyword = self
                        .instance(item)
                        .map(|i| i.keyword().to_owned())
                        .unwrap_or_default();
                    match keyword.as_str() {
                        "AXIS2_PLACEMENT_3D" => frame = self.frame(item).ok(),
                        "LENGTH_MEASURE_WITH_UNIT" | "MEASURE_REPRESENTATION_ITEM" => {
                            if let Some((value, _)) = self.measure_value(item) {
                                lengths.push(value);
                            }
                        }
                        _ => {}
                    }
                }
            }
        }

        // What the target *is* follows from the description the file gives
        // and how many sizes it carries: an area needs two, a circle one, a
        // line one, a point none.
        let kind = if description.contains("circle") || description.contains("circular") {
            ogeom_doc::DatumTargetKind::Circle {
                diameter: lengths.first().copied().unwrap_or(0.0),
            }
        } else if description.contains("rectangle") || lengths.len() >= 2 {
            ogeom_doc::DatumTargetKind::Rectangle {
                length: lengths.first().copied().unwrap_or(0.0),
                width: lengths.get(1).copied().unwrap_or(0.0),
            }
        } else if description.contains("line") || lengths.len() == 1 {
            ogeom_doc::DatumTargetKind::Line {
                length: lengths.first().copied().unwrap_or(0.0),
            }
        } else {
            ogeom_doc::DatumTargetKind::Point
        };
        Some(ogeom_doc::DatumTarget {
            datum,
            index,
            kind,
            at: frame.map_or(Point::ORIGIN, |f: Frame| f.origin()),
            frame,
            items: items_for(id),
        })
    }

    /// The drawn annotations: callouts, their polylines, their planes, and
    /// which semantic annotation each one draws.
    fn callouts(
        &mut self,
        annotation_ids: &FastMap<u64, ogeom_doc::Annotated>,
        aspect_items: &FastMap<u64, Vec<ogeom_topo::TShapeId>>,
    ) -> (Vec<ogeom_doc::Callout>, FastMap<u64, usize>) {
        // Which callout each annotation plane holds, and the plane's frame.
        let mut plane_of: FastMap<u64, Frame> = FastMap::default();
        for id in self.ids_with("ANNOTATION_PLANE") {
            let Ok(args) = self.args(id, "ANNOTATION_PLANE") else {
                continue;
            };
            let frame = args
                .get(2)
                .and_then(Arg::reference)
                .and_then(|plane| {
                    let inner = self.args(plane, "PLANE").ok()?;
                    inner.get(1).and_then(Arg::reference)
                })
                .and_then(|placement| self.frame(placement).ok());
            let Some(frame) = frame else { continue };
            for element in args.get(3).and_then(Arg::list).unwrap_or(&[]) {
                if let Some(callout) = element.reference() {
                    plane_of.insert(callout, frame);
                }
            }
        }

        // Which semantic annotation each callout draws. The association names
        // the annotation's own STEP id, so the link is made by matching that
        // against the ids the semantic pass already resolved.
        // A callout may be associated more than once (once with the shape
        // aspect the annotation is *about*, once with the annotation itself),
        // so every association is kept and the first that resolves to a
        // semantic annotation is the answer.
        let mut draws: FastMap<u64, Vec<u64>> = FastMap::default();
        let mut associations = self.ids_with("DRAUGHTING_MODEL_ITEM_ASSOCIATION");
        associations.sort_unstable();
        for id in associations {
            let Ok(args) = self.args(id, "DRAUGHTING_MODEL_ITEM_ASSOCIATION") else {
                continue;
            };
            if let (Some(definition), Some(item)) = (
                args.get(2).and_then(Arg::reference),
                args.get(4).and_then(Arg::reference),
            ) {
                draws.entry(item).or_default().push(definition);
            }
        }

        let mut out = Vec::new();
        let mut index_of: FastMap<u64, usize> = FastMap::default();
        let mut callouts = self.ids_with("DRAUGHTING_CALLOUT");
        callouts.sort_unstable();
        for id in callouts {
            let Ok(args) = self.args(id, "DRAUGHTING_CALLOUT") else {
                continue;
            };
            let name = match args.first() {
                Some(Arg::Str(s)) => s.clone(),
                _ => String::new(),
            };
            let mut polylines = Vec::new();
            let mut text: Option<String> = None;
            for element in args.get(1).and_then(Arg::list).unwrap_or(&[]).to_vec() {
                let Some(occurrence) = element.reference() else {
                    continue;
                };
                if let Some(words) = self.annotation_text(occurrence) {
                    text.get_or_insert(words);
                    continue;
                }
                polylines.extend(self.annotation_polylines(occurrence));
            }
            if polylines.is_empty() && text.is_none() && !plane_of.contains_key(&id) {
                continue;
            }
            index_of.insert(id, out.len());
            let definitions = draws.get(&id).cloned().unwrap_or_default();
            let annotates = definitions
                .iter()
                .find_map(|d| annotation_ids.get(d).copied());
            // A note draws no semantic annotation; what it is about is the
            // shape aspect it is associated with.
            let mut items: Vec<ogeom_topo::TShapeId> = Vec::new();
            if annotates.is_none() {
                for d in &definitions {
                    for &node in aspect_items.get(d).into_iter().flatten() {
                        if !items.contains(&node) {
                            items.push(node);
                        }
                    }
                }
            }
            out.push(ogeom_doc::Callout {
                name,
                plane: plane_of.get(&id).copied(),
                polylines,
                annotates,
                items,
                text,
            });
        }
        (out, index_of)
    }

    /// The words an annotation text occurrence carries: its text literal's.
    fn annotation_text(&mut self, occurrence: u64) -> Option<String> {
        let literal = {
            let instance = self.instance(occurrence).ok()?;
            let args = instance.part("ANNOTATION_TEXT_OCCURRENCE")?;
            args.get(2).and_then(Arg::reference)?
        };
        let args = self.args(literal, "TEXT_LITERAL").ok()?;
        match args.get(1) {
            Some(Arg::Str(words)) => Some(words.clone()),
            _ => None,
        }
    }

    /// The polylines one annotation occurrence draws.
    ///
    /// A tessellated occurrence names a curve set, which names a coordinates
    /// list and gives each polyline as *one-based* indices into it. That is
    /// the whole of the drawn geometry, and it is read as it is written
    /// rather than resampled.
    fn annotation_polylines(&mut self, occurrence: u64) -> Vec<Vec<Point>> {
        let item = {
            let Ok(instance) = self.instance(occurrence) else {
                return Vec::new();
            };
            let args = instance
                .part("TESSELLATED_ANNOTATION_OCCURRENCE")
                .or_else(|| instance.part("ANNOTATION_OCCURRENCE"))
                .or_else(|| instance.part("STYLED_ITEM"));
            match args.and_then(|a| a.get(2).and_then(Arg::reference)) {
                Some(item) => item,
                None => return Vec::new(),
            }
        };
        self.tessellated_polylines(item, 0)
    }

    /// The polylines a tessellated item holds, however deep it nests them.
    ///
    /// The drawn geometry of one annotation is a *set*: a frame's box, its
    /// leader, its text strokes, each its own curve set over its own
    /// coordinates list, gathered under one item, which may itself be
    /// repositioned by a placement, and that placement is applied here rather
    /// than left for a consumer to discover.
    fn tessellated_polylines(&mut self, item: u64, depth: usize) -> Vec<Vec<Point>> {
        if depth > 4 {
            return Vec::new();
        }
        let (children, curve_set, placement) = {
            let Ok(instance) = self.instance(item) else {
                return Vec::new();
            };
            let children: Vec<u64> = instance
                .part("TESSELLATED_GEOMETRIC_SET")
                .and_then(|args| args.first().and_then(Arg::list))
                .map(|list| list.iter().filter_map(Arg::reference).collect())
                .unwrap_or_default();
            let curve_set = instance.part("TESSELLATED_CURVE_SET");
            let placement = instance
                .part("REPOSITIONED_TESSELLATED_ITEM")
                .and_then(|args| args.first().and_then(Arg::reference));
            (children, curve_set, placement)
        };

        let mut out = Vec::new();
        for child in children {
            out.extend(self.tessellated_polylines(child, depth + 1));
        }
        if let Some(args) = curve_set
            && let Some(list) = args.get(1).and_then(Arg::reference)
        {
            let points = self.coordinates_list(list);
            for line in args.get(2).and_then(Arg::list).unwrap_or(&[]) {
                let Some(indices) = line.list() else { continue };
                let mut polyline = Vec::with_capacity(indices.len());
                for index in indices {
                    // One-based, as the format states them.
                    let Some(k) = index.number() else { continue };
                    #[expect(
                        clippy::cast_possible_truncation,
                        clippy::cast_sign_loss,
                        reason = "an index into a list the file itself sized"
                    )]
                    let k = k as usize;
                    if k >= 1 && k <= points.len() {
                        polyline.push(points[k - 1]);
                    }
                }
                if polyline.len() >= 2 {
                    out.push(polyline);
                }
            }
        }
        if let Some(placement) = placement
            && let Ok(frame) = self.frame(placement)
        {
            for polyline in &mut out {
                for p in polyline.iter_mut() {
                    *p = frame.to_world(*p);
                }
            }
        }
        out
    }

    /// A `COORDINATES_LIST`'s points, in the document's own length unit.
    fn coordinates_list(&mut self, id: u64) -> Vec<Point> {
        let Ok(args) = self.args(id, "COORDINATES_LIST") else {
            return Vec::new();
        };
        let scale = self.report.scale_mm;
        args.get(2)
            .and_then(Arg::list)
            .unwrap_or(&[])
            .iter()
            .filter_map(|entry| {
                let coords = entry.list()?;
                let value = |i: usize| coords.get(i).and_then(Arg::number).unwrap_or(0.0) * scale;
                Some(Point::new(value(0), value(1), value(2)))
            })
            .collect()
    }

    /// Every instance id carrying a part with this keyword.
    fn ids_with(&self, keyword: &str) -> Vec<u64> {
        self.instances
            .iter()
            .filter(|(_, inst)| inst.part(keyword).is_some())
            .map(|(id, _)| id)
            .collect()
    }

    /// One pass over the exchange builds both property-chain indexes:
    /// definition → its `PROPERTY_DEFINITION`s (by the definition argument),
    /// property → its `SHAPE_DEFINITION_REPRESENTATION`s (by the definition
    /// they represent). Ascending ids inside each list, so a lookup visits
    /// candidates in file order, lowest id first.
    fn ensure_property_indexes(&mut self) {
        if self.properties_of_definition.is_some() {
            return;
        }
        let mut properties: FastMap<u64, Vec<u64>> = FastMap::default();
        let mut sdrs: FastMap<u64, Vec<u64>> = FastMap::default();
        let instances = self.instances;
        for (id, instance) in instances.iter() {
            if let Some(args) = instance.part("PROPERTY_DEFINITION") {
                // Read for the index is read: the skipped table should not
                // claim the reader never looked.
                self.visited.mark(id);
                if let Some(definition) = args.get(2).and_then(Arg::reference) {
                    properties.entry(definition).or_default().push(id);
                }
            }
            if let Some(args) = instance.part("SHAPE_DEFINITION_REPRESENTATION") {
                self.visited.mark(id);
                if let Some(property) = args.first().and_then(Arg::reference) {
                    sdrs.entry(property).or_default().push(id);
                }
            }
        }
        for list in properties.values_mut().chain(sdrs.values_mut()) {
            list.sort_unstable();
        }
        self.properties_of_definition = Some(properties);
        self.sdrs_of_property = Some(sdrs);
    }

    /// A measure item's value and kind, scaled into the document's units.
    fn measure_value(&mut self, id: u64) -> Option<(f64, ogeom_doc::MeasureKind)> {
        let args = {
            let instance = self.instance(id).ok()?;
            instance.part("MEASURE_WITH_UNIT")?.to_vec()
        };
        match args.first() {
            Some(Arg::Typed(typed)) => {
                let (kind, inner) = (&typed.0, &typed.1);
                let value = inner.first().and_then(Arg::number)?;
                if kind.contains("ANGLE") {
                    Some((value * self.angle_scale, ogeom_doc::MeasureKind::Angle))
                } else {
                    Some((value * self.report.scale_mm, ogeom_doc::MeasureKind::Length))
                }
            }
            _ => None,
        }
    }

    /// A dimensional characteristic's name, kind and aspects.
    ///
    /// Sizes apply to one feature; locations (linear or angular) run
    /// between two. `ANGULAR_SIZE` and `ANGULAR_LOCATION` are the same
    /// shapes with an extra angle-selection argument at the end.
    fn dimension_shape(&mut self, dim: u64) -> (String, bool, Vec<u64>) {
        let size = {
            let Ok(instance) = self.instance(dim) else {
                return (String::new(), false, Vec::new());
            };
            instance
                .part("DIMENSIONAL_SIZE")
                .or_else(|| instance.part("ANGULAR_SIZE"))
        };
        if let Some(args) = size {
            let name = match args.get(1) {
                Some(Arg::Str(s)) => s.clone(),
                _ => String::new(),
            };
            let aspects = args.first().and_then(Arg::reference).into_iter().collect();
            return (name, false, aspects);
        }
        let location = {
            let Ok(instance) = self.instance(dim) else {
                return (String::new(), false, Vec::new());
            };
            instance
                .part("DIMENSIONAL_LOCATION")
                .or_else(|| instance.part("ANGULAR_LOCATION"))
        };
        if let Some(args) = location {
            let name = match args.first() {
                Some(Arg::Str(s)) => s.clone(),
                _ => String::new(),
            };
            let aspects = [args.get(2), args.get(3)]
                .into_iter()
                .flatten()
                .filter_map(Arg::reference)
                .collect();
            return (name, true, aspects);
        }
        (String::new(), false, Vec::new())
    }

    /// The references of a datum system, primary first: each compartment's
    /// datum (a common datum's letters hyphen-joined) and its modifiers.
    /// `None` when `id` is no datum system, or names a datum no letter
    /// resolves, the caller then reading it as a plain reference.
    fn datum_system(&mut self, id: u64) -> Option<Vec<ogeom_doc::DatumReference>> {
        let compartments: Vec<u64> = {
            let instance = self.instance(id).ok()?;
            instance
                .part("DATUM_SYSTEM")?
                .get(4)
                .and_then(Arg::list)?
                .iter()
                .filter_map(Arg::reference)
                .collect()
        };
        let mut out = Vec::with_capacity(compartments.len());
        for compartment in compartments {
            let args = self.args(compartment, "DATUM_REFERENCE_COMPARTMENT").ok()?;
            let label = match args.get(4)? {
                Arg::Ref(base) => self.datum_letter(*base, 0)?,
                Arg::List(elements) => {
                    // A common datum: one element per constituent, each
                    // naming its datum as a compartment does.
                    let mut letters = Vec::with_capacity(elements.len());
                    for element in elements {
                        let element = element.reference()?;
                        let datum = match self.args(element, "DATUM_REFERENCE_ELEMENT") {
                            Ok(args) => args.get(4).and_then(Arg::reference)?,
                            Err(_) => element,
                        };
                        letters.push(self.datum_letter(datum, 0)?);
                    }
                    letters.join("-")
                }
                _ => return None,
            };
            let mut modifiers = Vec::new();
            for modifier in args.get(5).and_then(Arg::list).unwrap_or(&[]) {
                match modifier {
                    Arg::Typed(typed) => {
                        if let Some(Arg::Enum(word)) = typed.1.first() {
                            modifiers.push(ogeom_doc::DatumModifier::Simple(word.to_lowercase()));
                        }
                    }
                    Arg::Enum(word) => {
                        modifiers.push(ogeom_doc::DatumModifier::Simple(word.to_lowercase()));
                    }
                    Arg::Ref(valued) => {
                        if let Some(modifier) = self.datum_modifier_with_value(*valued) {
                            modifiers.push(modifier);
                        }
                    }
                    _ => {}
                }
            }
            out.push(ogeom_doc::DatumReference { label, modifiers });
        }
        Some(out)
    }

    /// A `DATUM_REFERENCE_MODIFIER_WITH_VALUE`: its kind and its length.
    fn datum_modifier_with_value(&mut self, id: u64) -> Option<ogeom_doc::DatumModifier> {
        let args = self.args(id, "DATUM_REFERENCE_MODIFIER_WITH_VALUE").ok()?;
        let kind = match args.first()? {
            Arg::Enum(word) => ogeom_doc::DatumModifierKind::from_word(word)?,
            _ => return None,
        };
        let value = self.length_value(args.get(1)?.reference()?)?;
        Some(ogeom_doc::DatumModifier::WithValue { kind, value })
    }

    /// A length measure with unit, in millimetres, whether written as the
    /// complex measure item or as the plain entity.
    fn length_value(&mut self, id: u64) -> Option<f64> {
        if let Some((value, ogeom_doc::MeasureKind::Length)) = self.measure_value(id) {
            return Some(value);
        }
        let args = self.args(id, "LENGTH_MEASURE_WITH_UNIT").ok()?;
        match args.first()? {
            Arg::Typed(typed) => {
                Some(typed.1.first().and_then(Arg::number)? * self.report.scale_mm)
            }
            _ => None,
        }
    }

    /// The datum letter reachable from a datum reference, depth-bounded: the
    /// reference chain runs through compartments and systems, and files
    /// arrange it differently.
    fn datum_letter(&mut self, id: u64, depth: usize) -> Option<String> {
        if depth > 4 {
            return None;
        }
        let (letter, members, refs) = {
            let instance = self.instance(id).ok()?;
            let letter = instance.part("DATUM").and_then(|args| match args.get(4) {
                Some(Arg::Str(s)) if !s.is_empty() => Some(s.clone()),
                _ => None,
            });
            // A compartment or common datum carrying a list of constituent
            // references is a composite: every constituent contributes a
            // letter, and the letters act as one datum.
            let members: Vec<u64> = ["DATUM_REFERENCE_COMPARTMENT", "COMMON_DATUM"]
                .iter()
                .find_map(|k| instance.part(k))
                .into_iter()
                .flatten()
                .find_map(Arg::list)
                .map(|list| list.iter().filter_map(Arg::reference).collect())
                .unwrap_or_default();
            let mut refs = Vec::new();
            for (_, args) in instance.parts() {
                collect_refs(args, &mut refs);
            }
            (letter, members, refs)
        };
        if let Some(letter) = letter {
            return Some(letter);
        }
        if !members.is_empty() {
            let letters: Vec<String> = members
                .into_iter()
                .filter_map(|r| self.datum_letter(r, depth + 1))
                .collect();
            if !letters.is_empty() {
                return Some(letters.join("-"));
            }
        }
        refs.into_iter()
            .find_map(|r| self.datum_letter(r, depth + 1))
    }
}
