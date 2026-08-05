//! Writing STEP: the document's products, assemblies, colours and B-rep.
//!
//! The mirror of the reader, and deliberately written against the same
//! vocabulary: every entity this writer emits is one the reader parses, so
//! writing what was read and reading it back is the honest round-trip test.
//! AP214's schema name goes in the header — the entities used are the common
//! AP203/AP214/AP242 core.
//!
//! Geometry is written in world coordinates: every face, edge and vertex is
//! transformed through its occurrence's own placement chain before it is
//! serialized, and shared nodes are deduplicated *per placement* — a prism's
//! bottom and top edge are one node at two locations, and the file needs
//! both. Surfaces the format has no analytic name for — extrusions,
//! revolutions — go out as their exact rational B-spline patches (§3's
//! conversion), so nothing is fitted on the way out.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_doc::{Document, ProductId, ProductKind};
use ogeom_geom::Transformable as _;
use ogeom_geom::{Curve, SurfaceGeometry};
use ogeom_math::{Frame, Point, Transform, Vector};
use ogeom_topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, explore};
use std::collections::HashMap;
use std::fmt::Write as _;

/// Write a document as a STEP exchange file.
///
/// Products become `PRODUCT` trees; parts carry their solids as
/// `MANIFOLD_SOLID_BREP`s; assemblies become usage occurrences with their
/// placements; colours become styled items over the written solids and
/// faces.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a
/// shape's structure cannot be expressed — a non-rigid instance placement, a
/// solid with no shell.
pub fn write_step(document: &Document, tol: Tolerances) -> OgeomResult<String> {
    let mut writer = Writer {
        model: document.model(),
        entities: Vec::new(),
        points: HashMap::new(),
        directions: HashMap::new(),
        vertices: HashMap::new(),
        edges: HashMap::new(),
        written_nodes: Vec::new(),
        tol,
    };

    let app = writer.entity("APPLICATION_CONTEXT('automotive design')".into());
    let _proto = writer.entity(format!(
        "APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2010,#{app})"
    ));
    let lu = writer.entity("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))".into());
    let au = writer.entity("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))".into());
    let su = writer.entity("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())".into());
    let unc = writer.entity(format!(
        "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.0E-06),#{lu},'distance_accuracy_value','')"
    ));
    let gctx = writer.entity(format!(
        "(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{unc}))GLOBAL_UNIT_ASSIGNED_CONTEXT((#{lu},#{au},#{su}))REPRESENTATION_CONTEXT('Context','3D'))"
    ));
    let pctx = writer.entity(format!("PRODUCT_CONTEXT('',#{app},'mechanical')"));
    let pdctx = writer.entity(format!(
        "PRODUCT_DEFINITION_CONTEXT('part definition',#{app},'design')"
    ));

    // Products, in document order; every product gets its definition and its
    // shape representation before the assembly edges tie them together.
    let mut pd_of: HashMap<ProductId, u64> = HashMap::new();
    let mut sr_of: HashMap<ProductId, u64> = HashMap::new();
    let mut anchor_pds: Option<(u64, u64)> = None;
    for (id, product) in document.products() {
        let name = escape(&product.name);
        let p = writer.entity(format!("PRODUCT('{name}','{name}','',(#{pctx}))"));
        let formation = writer.entity(format!("PRODUCT_DEFINITION_FORMATION('','',#{p})"));
        let pd = writer.entity(format!(
            "PRODUCT_DEFINITION('design','',#{formation},#{pdctx})"
        ));
        pd_of.insert(id, pd);

        let world = writer.frame(&Frame::WORLD);
        let sr = match &product.kind {
            ProductKind::Part { shape } => {
                let mut items = vec![world];
                for solid in writer.solids_of(shape)? {
                    items.push(solid);
                }
                let list = items
                    .iter()
                    .map(|i| format!("#{i}"))
                    .collect::<Vec<_>>()
                    .join(",");
                writer.entity(format!(
                    "ADVANCED_BREP_SHAPE_REPRESENTATION('{name}',({list}),#{gctx})"
                ))
            }
            ProductKind::Assembly { .. } => {
                writer.entity(format!("SHAPE_REPRESENTATION('{name}',(#{world}),#{gctx})"))
            }
        };
        sr_of.insert(id, sr);
        let pds = writer.entity(format!("PRODUCT_DEFINITION_SHAPE('','',#{pd})"));
        writer.entity(format!("SHAPE_DEFINITION_REPRESENTATION(#{pds},#{sr})"));
        if anchor_pds.is_none() && matches!(product.kind, ProductKind::Part { .. }) {
            anchor_pds = Some((pds, sr));
        }
    }

    // Assembly edges: one usage occurrence per instance, its placement said
    // through the transformation between the parent's world frame and the
    // child's placement frame.
    let mut usage = 0_usize;
    for (id, product) in document.products() {
        let ProductKind::Assembly { children } = &product.kind else {
            continue;
        };
        for instance in children {
            usage += 1;
            let designator = instance
                .name
                .clone()
                .unwrap_or_else(|| format!("occurrence-{usage}"));
            let designator = escape(&designator);
            let (parent_pd, child_pd) = (pd_of[&id], pd_of[&instance.product]);
            let (parent_sr, child_sr) = (sr_of[&id], sr_of[&instance.product]);
            let nauo = writer.entity(format!(
                "NEXT_ASSEMBLY_USAGE_OCCURRENCE('{designator}','{designator}','',#{parent_pd},#{child_pd},$)"
            ));
            // The location resolved through the model's own datum store: a
            // placed dummy shape shares the resolution path every traversal
            // uses.
            let at = location_transform(&instance.location, document.model())?;
            let placed = writer.placement_frame(&at)?;
            let world = writer.frame(&Frame::WORLD);
            let idt = writer.entity(format!(
                "ITEM_DEFINED_TRANSFORMATION('','',#{world},#{placed})"
            ));
            let rr = writer.entity(format!(
                "(REPRESENTATION_RELATIONSHIP('','',#{child_sr},#{parent_sr})REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#{idt})SHAPE_REPRESENTATION_RELATIONSHIP())"
            ));
            let pds = writer.entity(format!("PRODUCT_DEFINITION_SHAPE('','',#{nauo})"));
            writer.entity(format!(
                "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#{rr},#{pds})"
            ));
        }
    }

    // Colours: a styled item over every written entity whose node the
    // document colours, plus product colours carried by their solids.
    let mut styled = Vec::new();
    let node_colours: HashMap<_, _> = document.colours().collect();
    for (node, step_id) in writer.written_nodes.clone() {
        if let Some(colour) = node_colours.get(&node) {
            styled.push(writer.styled_item(step_id, *colour));
        }
    }
    for (id, product) in document.products() {
        let Some(colour) = product.colour else {
            continue;
        };
        let ProductKind::Part { shape } = &product.kind else {
            continue;
        };
        let _ = id;
        for (node, step_id) in writer.written_nodes.clone() {
            if node == shape.node() && !node_colours.contains_key(&node) {
                styled.push(writer.styled_item(step_id, colour));
            }
        }
    }
    if !styled.is_empty() {
        let list = styled
            .iter()
            .map(|i| format!("#{i}"))
            .collect::<Vec<_>>()
            .join(",");
        writer.entity(format!(
            "MECHANICAL_DESIGN_GEOMETRIC_PRESENTATION_REPRESENTATION('',({list}),#{gctx})"
        ));
    }

    // Semantic PMI: datums first, so tolerances can reference their letters.
    let pmi = document.pmi();
    if !pmi.is_empty() {
        let Some((pds, absr)) = anchor_pds else {
            ogeom_bail!(
                Construction,
                "PMI needs at least one part to anchor its aspects to"
            );
        };
        writer.pmi(pmi, pds, absr, lu, au, gctx)?;
    }

    let mut out = String::new();
    out.push_str("ISO-10303-21;\nHEADER;\n");
    out.push_str("FILE_DESCRIPTION(('written by ogeom'),'2;1');\n");
    out.push_str("FILE_NAME('','',('ogeom'),('ogeom'),'ogeom','ogeom','');\n");
    out.push_str("FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\n");
    out.push_str("ENDSEC;\nDATA;\n");
    for (i, entity) in writer.entities.iter().enumerate() {
        let _ = writeln!(out, "#{}={entity};", i + 1);
    }
    out.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
    Ok(out)
}

/// The state of one write: the entity buffer and the per-placement caches.
struct Writer<'a> {
    model: &'a Model,
    entities: Vec<String>,
    points: HashMap<[u64; 3], u64>,
    directions: HashMap<[u64; 3], u64>,
    /// Vertex occurrences by node and world position bits.
    vertices: HashMap<(ogeom_topo::TShapeId, [u64; 3]), u64>,
    /// Edge occurrences by node and placement bits.
    edges: HashMap<(ogeom_topo::TShapeId, [u64; 3]), u64>,
    /// Every solid and face written, with its entity id — the hooks colours
    /// attach to.
    written_nodes: Vec<(ogeom_topo::TShapeId, u64)>,
    tol: Tolerances,
}

impl Writer<'_> {
    fn entity(&mut self, text: String) -> u64 {
        self.entities.push(text);
        self.entities.len() as u64
    }

    fn point(&mut self, p: Point) -> u64 {
        let key = [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
        if let Some(&id) = self.points.get(&key) {
            return id;
        }
        let id = self.entity(format!(
            "CARTESIAN_POINT('',({},{},{}))",
            real(p.x),
            real(p.y),
            real(p.z)
        ));
        self.points.insert(key, id);
        id
    }

    fn direction(&mut self, v: Vector) -> u64 {
        let key = [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()];
        if let Some(&id) = self.directions.get(&key) {
            return id;
        }
        let id = self.entity(format!(
            "DIRECTION('',({},{},{}))",
            real(v.x),
            real(v.y),
            real(v.z)
        ));
        self.directions.insert(key, id);
        id
    }

    fn frame(&mut self, frame: &Frame) -> u64 {
        let origin = self.point(frame.origin());
        let z = self.direction(frame.z().vector());
        let x = self.direction(frame.x().vector());
        self.entity(format!("AXIS2_PLACEMENT_3D('',#{origin},#{z},#{x})"))
    }

    /// A rigid transform as the frame it carries the world onto.
    fn placement_frame(&mut self, at: &Transform) -> OgeomResult<u64> {
        let origin = at.apply(Point::ORIGIN);
        let z = at.apply_vector(Vector::new(0.0, 0.0, 1.0));
        let x = at.apply_vector(Vector::new(1.0, 0.0, 0.0));
        let frame = Frame::new(
            origin,
            ogeom_math::Direction::new(z, self.tol)?,
            ogeom_math::Direction::new(x, self.tol)?,
            self.tol,
        )
        .map_err(|_| {
            ogeom_core::ogeom_err!(
                Construction,
                "an instance placement is not rigid; STEP cannot state it"
            )
        })?;
        Ok(self.frame(&frame))
    }

    /// Every solid under a part's shape, written.
    fn solids_of(&mut self, shape: &Shape) -> OgeomResult<Vec<u64>> {
        let mut out = Vec::new();
        for solid in explore(self.model, shape, Filter::OfType(ShapeType::Solid))? {
            out.push(self.solid(&solid)?);
        }
        if out.is_empty() {
            ogeom_bail!(Construction, "a part's shape holds no solid to write");
        }
        Ok(out)
    }

    fn solid(&mut self, solid: &Shape) -> OgeomResult<u64> {
        let shells = explore(self.model, solid, Filter::OfType(ShapeType::Shell))?;
        let Some(shell) = shells.first() else {
            ogeom_bail!(Construction, "a solid with no shell cannot be written");
        };
        let mut faces = Vec::new();
        for face in self.model.ordered_children_of(shell)? {
            faces.push(self.face(&face)?);
        }
        let list = faces
            .iter()
            .map(|i| format!("#{i}"))
            .collect::<Vec<_>>()
            .join(",");
        let shell_id = self.entity(format!("CLOSED_SHELL('',({list}))"));
        let msb = self.entity(format!("MANIFOLD_SOLID_BREP('',#{shell_id})"));
        self.written_nodes.push((solid.node(), msb));
        Ok(msb)
    }

    fn face(&mut self, face: &Shape) -> OgeomResult<u64> {
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
        let surface_id = self.surface(&surface)?;

        let mut bounds = Vec::new();
        for (index, wire) in self.model.ordered_children_of(face)?.iter().enumerate() {
            let keyword = if index == 0 {
                "FACE_OUTER_BOUND"
            } else {
                "FACE_BOUND"
            };
            let loop_id = self.wire(wire)?;
            bounds.push(self.entity(format!("{keyword}('',#{loop_id},.T.)")));
        }
        let list = bounds
            .iter()
            .map(|i| format!("#{i}"))
            .collect::<Vec<_>>()
            .join(",");
        let sense = if face.orientation() == ogeom_topo::Orientation::Reversed {
            ".F."
        } else {
            ".T."
        };
        let id = self.entity(format!("ADVANCED_FACE('',({list}),#{surface_id},{sense})"));
        self.written_nodes.push((face.node(), id));
        Ok(id)
    }

    /// A wire as an `EDGE_LOOP`, or a `VERTEX_LOOP` when every edge in it is
    /// degenerate — a pole or an apex has no curve to serialize, and STEP's
    /// own spelling for it is the loop of one vertex.
    fn wire(&mut self, wire: &Shape) -> OgeomResult<u64> {
        let children = self.model.ordered_children_of(wire)?;
        let degenerate = |edge: &Shape| {
            self.model
                .node(edge)
                .and_then(|n| n.data().as_edge())
                .is_some_and(|d| d.degenerate)
        };
        if !children.is_empty() && children.iter().all(degenerate) {
            let edge = &children[0];
            let Some(vertex) = self.model.children_of(edge)?.first().cloned() else {
                ogeom_bail!(Construction, "a degenerate edge has no vertex");
            };
            let vertex_id = self.vertex(&vertex, &edge.transform(self.model.datums())?)?;
            return Ok(self.entity(format!("VERTEX_LOOP('',#{vertex_id})")));
        }
        let mut oriented = Vec::new();
        for edge in &children {
            if degenerate(edge) {
                continue;
            }
            let edge_id = self.edge(edge)?;
            let sense = if edge.orientation() == ogeom_topo::Orientation::Reversed {
                ".F."
            } else {
                ".T."
            };
            oriented.push(self.entity(format!("ORIENTED_EDGE('',*,*,#{edge_id},{sense})")));
        }
        let list = oriented
            .iter()
            .map(|i| format!("#{i}"))
            .collect::<Vec<_>>()
            .join(",");
        Ok(self.entity(format!("EDGE_LOOP('',({list}))")))
    }

    fn edge(&mut self, edge: &Shape) -> OgeomResult<u64> {
        let placement = edge.transform(self.model.datums())?;
        let key = (edge.node(), transform_bits(&placement));
        if let Some(&id) = self.edges.get(&key) {
            return Ok(id);
        }
        let (curve, range) = {
            let Some(node) = self.model.node(edge) else {
                ogeom_bail!(Dangling, "edge is not in this model");
            };
            let Some(data) = node.data().as_edge() else {
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
        let curve_id = self.curve(&curve, range)?;
        let vertices = self.model.children_of(edge)?;
        let (from, to) = match vertices.len() {
            0 => ogeom_bail!(Construction, "an edge with no vertices cannot be written"),
            1 => (vertices[0].clone(), vertices[0].clone()),
            _ => (vertices[0].clone(), vertices[vertices.len() - 1].clone()),
        };
        let from_id = self.vertex(&from, &placement)?;
        let to_id = self.vertex(&to, &placement)?;
        let id = self.entity(format!(
            "EDGE_CURVE('',#{from_id},#{to_id},#{curve_id},.T.)"
        ));
        self.edges.insert(key, id);
        Ok(id)
    }

    fn vertex(&mut self, vertex: &Shape, placement: &Transform) -> OgeomResult<u64> {
        let Some(data) = self.model.node(vertex).and_then(|n| n.data().as_vertex()) else {
            ogeom_bail!(Construction, "vertex node holds no vertex data");
        };
        let at = placement.apply(data.point);
        let key = (
            vertex.node(),
            [at.x.to_bits(), at.y.to_bits(), at.z.to_bits()],
        );
        if let Some(&id) = self.vertices.get(&key) {
            return Ok(id);
        }
        let point = self.point(at);
        let id = self.entity(format!("VERTEX_POINT('',#{point})"));
        self.vertices.insert(key, id);
        Ok(id)
    }

    /// A curve for an `EDGE_CURVE`: the analytic spelling where STEP has
    /// one, the exact B-spline conversion where it does not.
    fn curve(&mut self, curve: &Curve, range: (f64, f64)) -> OgeomResult<u64> {
        match curve {
            Curve::Line(line) => {
                let axis = line.axis();
                let origin = self.point(axis.location);
                let d = self.direction(axis.direction.vector());
                let vector = self.entity(format!("VECTOR('',#{d},1.0)"));
                Ok(self.entity(format!("LINE('',#{origin},#{vector})")))
            }
            Curve::Circle(c) => {
                let circle = c.circle();
                let frame = self.frame(&circle.frame());
                Ok(self.entity(format!("CIRCLE('',#{frame},{})", real(circle.radius()))))
            }
            Curve::Ellipse(el) => {
                let ellipse = el.ellipse();
                let frame = self.frame(&ellipse.frame());
                Ok(self.entity(format!(
                    "ELLIPSE('',#{frame},{},{})",
                    real(ellipse.major_radius()),
                    real(ellipse.minor_radius())
                )))
            }
            Curve::BSpline(b) => self.bspline_curve(b),
            other => {
                // Exact for every remaining kind this kernel has: the
                // conversion is the §3 machinery, not a fit.
                let spline = other.to_bspline_over(range, self.tol)?;
                self.bspline_curve(&spline)
            }
        }
    }

    fn bspline_curve(&mut self, spline: &ogeom_geom::BSplineCurve) -> OgeomResult<u64> {
        let control: Vec<String> = spline
            .control_points()
            .iter()
            .map(|c| {
                let p = Point::from_vector(c.scaled.to_vector() / c.weight);
                format!("#{}", self.point(p))
            })
            .collect();
        let (mults, knots) = compress_knots(spline.knots().knots());
        let degree = spline.knots().degree();
        let rational = spline
            .control_points()
            .iter()
            .any(|c| (c.weight - 1.0).abs() > 1e-12);
        let control = control.join(",");
        let mults = mults
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(",");
        let knots = knots.iter().map(|k| real(*k)).collect::<Vec<_>>().join(",");
        if rational {
            let weights = spline
                .control_points()
                .iter()
                .map(|c| real(c.weight))
                .collect::<Vec<_>>()
                .join(",");
            Ok(self.entity(format!(
                "(BOUNDED_CURVE()B_SPLINE_CURVE({degree},({control}),.UNSPECIFIED.,.F.,.F.)B_SPLINE_CURVE_WITH_KNOTS(({mults}),({knots}),.UNSPECIFIED.)CURVE()GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE(({weights}))REPRESENTATION_ITEM(''))"
            )))
        } else {
            Ok(self.entity(format!(
                "B_SPLINE_CURVE_WITH_KNOTS('',{degree},({control}),.UNSPECIFIED.,.F.,.F.,({mults}),({knots}),.UNSPECIFIED.)"
            )))
        }
    }

    /// A surface: analytic where STEP has the word, the exact B-spline patch
    /// where it does not. Trimmed surfaces write their basis — the trim is
    /// the face's own topology.
    fn surface(&mut self, surface: &SurfaceGeometry) -> OgeomResult<u64> {
        match surface {
            SurfaceGeometry::Plane(s) => {
                let frame = self.frame(&s.plane().frame());
                Ok(self.entity(format!("PLANE('',#{frame})")))
            }
            SurfaceGeometry::Cylinder(s) => {
                let cylinder = s.cylinder();
                let frame = self.frame(&cylinder.frame());
                Ok(self.entity(format!(
                    "CYLINDRICAL_SURFACE('',#{frame},{})",
                    real(cylinder.radius())
                )))
            }
            SurfaceGeometry::Cone(s) => {
                let cone = s.cone();
                let frame = self.frame(&cone.frame());
                Ok(self.entity(format!(
                    "CONICAL_SURFACE('',#{frame},{},{})",
                    real(cone.radius_at(0.0)),
                    real(cone.half_angle())
                )))
            }
            SurfaceGeometry::Sphere(s) => {
                let sphere = s.sphere();
                let frame = self.frame(&sphere.frame());
                Ok(self.entity(format!(
                    "SPHERICAL_SURFACE('',#{frame},{})",
                    real(sphere.radius())
                )))
            }
            SurfaceGeometry::Torus(s) => {
                let torus = s.torus();
                let frame = self.frame(&torus.frame());
                Ok(self.entity(format!(
                    "TOROIDAL_SURFACE('',#{frame},{},{})",
                    real(torus.major_radius()),
                    real(torus.minor_radius())
                )))
            }
            SurfaceGeometry::BSpline(b) => self.bspline_surface(b),
            SurfaceGeometry::Trimmed(t) => self.surface(t.basis()),
            other => {
                let patch = other.to_bspline(self.tol)?;
                self.bspline_surface(&patch)
            }
        }
    }

    fn bspline_surface(&mut self, patch: &ogeom_geom::BSplineSurface) -> OgeomResult<u64> {
        let grid = patch.grid();
        let (nu, nv) = (grid.u_count(), grid.v_count());
        let mut rows = Vec::with_capacity(nu);
        let mut weights_rows = Vec::with_capacity(nu);
        let mut rational = false;
        for u in 0..nu {
            let mut row = Vec::with_capacity(nv);
            let mut wrow = Vec::with_capacity(nv);
            for v in 0..nv {
                let Some(c) = grid.get(u, v) else {
                    ogeom_bail!(Construction, "a control grid cell is missing");
                };
                let p = Point::from_vector(c.scaled.to_vector() / c.weight);
                row.push(format!("#{}", self.point(p)));
                wrow.push(real(c.weight));
                rational |= (c.weight - 1.0).abs() > 1e-12;
            }
            rows.push(format!("({})", row.join(",")));
            weights_rows.push(format!("({})", wrow.join(",")));
        }
        let grid_text = rows.join(",");
        let (u_deg, v_deg) = (patch.u_knots().degree(), patch.v_knots().degree());
        let (um, uk) = compress_knots(patch.u_knots().knots());
        let (vm, vk) = compress_knots(patch.v_knots().knots());
        let fmt_m = |m: &[usize]| {
            m.iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",")
        };
        let fmt_k = |k: &[f64]| k.iter().map(|v| real(*v)).collect::<Vec<_>>().join(",");
        let (um, uk, vm, vk) = (fmt_m(&um), fmt_k(&uk), fmt_m(&vm), fmt_k(&vk));
        if rational {
            let weights = weights_rows.join(",");
            Ok(self.entity(format!(
                "(BOUNDED_SURFACE()B_SPLINE_SURFACE({u_deg},{v_deg},({grid_text}),.UNSPECIFIED.,.F.,.F.,.F.)B_SPLINE_SURFACE_WITH_KNOTS(({um}),({vm}),({uk}),({vk}),.UNSPECIFIED.)GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_SURFACE(({weights}))REPRESENTATION_ITEM('')SURFACE())"
            )))
        } else {
            Ok(self.entity(format!(
                "B_SPLINE_SURFACE_WITH_KNOTS('',{u_deg},{v_deg},({grid_text}),.UNSPECIFIED.,.F.,.F.,.F.,({um}),({vm}),({uk}),({vk}),.UNSPECIFIED.)"
            )))
        }
    }

    /// The document's PMI, written over the aspects of the anchor part.
    #[allow(clippy::many_single_char_names)]
    fn pmi(
        &mut self,
        pmi: &ogeom_doc::Pmi,
        pds: u64,
        absr: u64,
        lu: u64,
        au: u64,
        gctx: u64,
    ) -> OgeomResult<()> {
        let by_node: HashMap<ogeom_topo::TShapeId, u64> =
            self.written_nodes.iter().copied().collect();
        let aspect_for = |w: &mut Self, items: &[ogeom_topo::TShapeId]| -> u64 {
            let aspect = w.entity(format!("SHAPE_ASPECT('','',#{pds},.T.)"));
            for item in items {
                if let Some(&step_id) = by_node.get(item) {
                    w.entity(format!(
                        "GEOMETRIC_ITEM_SPECIFIC_USAGE('','',#{aspect},#{absr},#{step_id})"
                    ));
                }
            }
            aspect
        };

        let mut datum_ids: HashMap<&str, u64> = HashMap::new();
        for datum in &pmi.datums {
            let label = escape(&datum.label);
            let id = self.entity(format!("DATUM('','',#{pds},.F.,'{label}')"));
            for item in &datum.items {
                if let Some(&step_id) = by_node.get(item) {
                    self.entity(format!(
                        "GEOMETRIC_ITEM_SPECIFIC_USAGE('','',#{id},#{absr},#{step_id})"
                    ));
                }
            }
            datum_ids.insert(datum.label.as_str(), id);
        }

        for dimension in &pmi.dimensions {
            let name = escape(&dimension.name);
            let angular = dimension.kind == ogeom_doc::MeasureKind::Angle;
            let dim = if dimension.location {
                // A location runs between two features; each keeps its own
                // aspect, so what was read as two ends writes as two ends.
                let empty = Vec::new();
                let first = dimension.features.first().unwrap_or(&empty);
                let second = dimension.features.get(1).unwrap_or(first);
                let a = aspect_for(self, first);
                let b = aspect_for(self, second);
                if angular {
                    self.entity(format!("ANGULAR_LOCATION('{name}','',#{a},#{b},.EQUAL.)"))
                } else {
                    self.entity(format!("DIMENSIONAL_LOCATION('{name}','',#{a},#{b})"))
                }
            } else {
                let empty = Vec::new();
                let items = dimension.features.first().unwrap_or(&empty);
                let aspect = aspect_for(self, items);
                if angular {
                    self.entity(format!("ANGULAR_SIZE(#{aspect},'{name}',.EQUAL.)"))
                } else {
                    self.entity(format!("DIMENSIONAL_SIZE(#{aspect},'{name}')"))
                }
            };
            let measures: Vec<String> = dimension
                .values
                .iter()
                .map(|&v| format!("#{}", self.measure(v, dimension.kind, lu, au)))
                .collect();
            let list = measures.join(",");
            let sdr = self.entity(format!(
                "SHAPE_DIMENSION_REPRESENTATION('',({list}),#{gctx})"
            ));
            self.entity(format!(
                "DIMENSIONAL_CHARACTERISTIC_REPRESENTATION(#{dim},#{sdr})"
            ));
            if dimension.plus.is_some() || dimension.minus.is_some() {
                let lower = self.measure(dimension.minus.unwrap_or(0.0), dimension.kind, lu, au);
                let upper = self.measure(dimension.plus.unwrap_or(0.0), dimension.kind, lu, au);
                let tv = self.entity(format!("TOLERANCE_VALUE(#{lower},#{upper})"));
                self.entity(format!("PLUS_MINUS_TOLERANCE(#{tv},#{dim})"));
            }
        }

        for tolerance in &pmi.tolerances {
            let aspect = aspect_for(self, &tolerance.items);
            let name = escape(&tolerance.name);
            let magnitude =
                self.measure(tolerance.magnitude, ogeom_doc::MeasureKind::Length, lu, au);
            let keyword = format!("{}_TOLERANCE", tolerance.kind.to_uppercase());
            if tolerance.datums.is_empty() {
                self.entity(format!("{keyword}('{name}','',#{magnitude},#{aspect})"));
            } else {
                let refs: Vec<String> = tolerance
                    .datums
                    .iter()
                    .filter_map(|d| datum_ids.get(d.as_str()))
                    .map(|i| format!("#{i}"))
                    .collect();
                let refs = refs.join(",");
                self.entity(format!(
                    "{keyword}('{name}','',#{magnitude},#{aspect},({refs}))"
                ));
            }
        }
        Ok(())
    }

    /// A measure representation item, in the document's own units.
    fn measure(&mut self, value: f64, kind: ogeom_doc::MeasureKind, lu: u64, au: u64) -> u64 {
        let v = real(value);
        match kind {
            ogeom_doc::MeasureKind::Length => self.entity(format!(
                "(LENGTH_MEASURE_WITH_UNIT()MEASURE_REPRESENTATION_ITEM()MEASURE_WITH_UNIT(LENGTH_MEASURE({v}),#{lu})REPRESENTATION_ITEM(''))"
            )),
            ogeom_doc::MeasureKind::Angle => self.entity(format!(
                "(MEASURE_REPRESENTATION_ITEM()MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE({v}),#{au})PLANE_ANGLE_MEASURE_WITH_UNIT()REPRESENTATION_ITEM(''))"
            )),
        }
    }

    /// A styled item colouring one written entity.
    fn styled_item(&mut self, item: u64, colour: ogeom_doc::Colour) -> u64 {
        let rgb = self.entity(format!(
            "COLOUR_RGB('',{},{},{})",
            real(colour.r),
            real(colour.g),
            real(colour.b)
        ));
        let fasc = self.entity(format!("FILL_AREA_STYLE_COLOUR('',#{rgb})"));
        let fas = self.entity(format!("FILL_AREA_STYLE('',(#{fasc}))"));
        let ssfa = self.entity(format!("SURFACE_STYLE_FILL_AREA(#{fas})"));
        let sss = self.entity(format!("SURFACE_SIDE_STYLE('',(#{ssfa}))"));
        let ssu = self.entity(format!("SURFACE_STYLE_USAGE(.BOTH.,#{sss})"));
        let psa = self.entity(format!("PRESENTATION_STYLE_ASSIGNMENT((#{ssu}))"));
        self.entity(format!("STYLED_ITEM('',(#{psa}),#{item})"))
    }
}

/// A Part 21 real: shortest round-trip form, decimal point guaranteed,
/// exponent uppercased.
fn real(v: f64) -> String {
    let mut s = format!("{v:?}");
    if let Some(e) = s.find(['e', 'E']) {
        let (mantissa, exponent) = s.split_at(e);
        let mut m = mantissa.to_string();
        if !m.contains('.') {
            m.push_str(".0");
        }
        s = format!("{m}E{}", &exponent[1..]);
    } else if !s.contains('.') {
        s.push_str(".0");
    }
    s
}

/// A string literal's body, quotes doubled per Part 21.
fn escape(s: &str) -> String {
    s.replace('\'', "''")
}

/// Knots as STEP states them: distinct values with multiplicities.
fn compress_knots(knots: &[f64]) -> (Vec<usize>, Vec<f64>) {
    let mut mults = Vec::new();
    let mut values = Vec::new();
    for &k in knots {
        match values.last() {
            Some(&last) if k == last => {
                if let Some(m) = mults.last_mut() {
                    *m += 1;
                }
            }
            _ => {
                values.push(k);
                mults.push(1);
            }
        }
    }
    (mults, values)
}

/// A rigid transform quantized to bits, for per-placement deduplication.
fn transform_bits(t: &Transform) -> [u64; 3] {
    let p = t.apply(Point::new(0.123_456_789, 9.87, -3.21));
    [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
}

/// A location resolved to the rigid transform it composes to.
fn location_transform(location: &ogeom_topo::Location, model: &Model) -> OgeomResult<Transform> {
    let mut out = Transform::IDENTITY;
    for &(datum, power) in location.chain() {
        let Some(t) = model.datums().get(datum) else {
            ogeom_bail!(Dangling, "an instance placement names a missing datum");
        };
        let step = if power >= 0 { t } else { t.inverse()? };
        for _ in 0..power.unsigned_abs() {
            out = out * step;
        }
    }
    Ok(out)
}
