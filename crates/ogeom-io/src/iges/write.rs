//! From a document to an IGES deck.
//!
//! The writer emits manifold solid B-rep objects — entity 186 over shells
//! (514), faces (510), loops (508), one edge list (504) and one vertex list
//! (502) per solid — with surfaces in their analytic spellings where IGES has
//! one and as rational B-splines (128) where it does not. Curves defined in a
//! plane of their own, arcs and ellipses, are written in definition space
//! with a transformation matrix (124) carrying them to model space, which is
//! how the format wants them.
//!
//! Everything is written in millimetres, model space, with geometry baked —
//! the same decision the STEP writer made: a file carries positions, not this
//! kernel's location chains.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::Surface as _;
use ogeom_geom::Transformable as _;
use ogeom_geom::{BSplineCurve, Curve, SurfaceGeometry};
use ogeom_math::{Frame, Point, Transform};
use ogeom_topo::{EdgeRepr, Filter, NodeData, Shape, ShapeType, explore};
use std::collections::HashMap;

/// Write a document's solids as an IGES file.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// document holds no solid, or a solid carries geometry with no IGES
/// spelling here yet — the error names the entity and the parity row.
pub fn write_iges(document: &ogeom_doc::Document, tol: Tolerances) -> OgeomResult<String> {
    let mut writer = Writer {
        model: document.model(),
        entities: Vec::new(),
        vertices: HashMap::new(),
        edges: HashMap::new(),
        vertex_coords: Vec::new(),
        edge_records: Vec::new(),
        tol,
    };

    let mut solids = Vec::new();
    for (_, product) in document.products() {
        let ogeom_doc::ProductKind::Part { shape } = &product.kind else {
            continue;
        };
        for solid in explore(writer.model, shape, Filter::OfType(ShapeType::Solid))? {
            solids.push((product.name.clone(), solid));
        }
    }
    if solids.is_empty() {
        ogeom_bail!(Construction, "the document holds no solid to write as IGES");
    }
    for (label, solid) in solids {
        writer.solid(&solid, &label)?;
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
    vertices: HashMap<(ogeom_topo::TShapeId, [u64; 3]), usize>,
    /// Edge index (1-based) in the current solid's 504, keyed by node and
    /// placement for the same reason.
    edges: HashMap<(ogeom_topo::TShapeId, [u64; 3]), usize>,
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

        let shells = explore(self.model, solid, Filter::OfType(ShapeType::Shell))?;
        let Some(shell) = shells.first() else {
            ogeom_bail!(Construction, "a solid with no shell cannot be written");
        };
        let mut face_entities = Vec::new();
        for face in self.model.ordered_children_of(shell)? {
            face_entities.push((self.face(&face)?, face.orientation()));
        }

        // Now the lists exist in full.
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
        // placeholders resolve now.
        let vlist_de = self.de(vertex_list);
        let elist_de = self.de(edge_list);
        for e in &mut self.entities {
            if e.kind == 508 || e.kind == 504 {
                e.params = e.params.replace("@E", &elist_de.to_string());
                e.params = e.params.replace("@V", &vlist_de.to_string());
            }
        }

        let shell_entity = {
            let mut params = format!("{}", face_entities.len());
            for (face, orientation) in &face_entities {
                let flag = i32::from(*orientation != ogeom_topo::Orientation::Reversed);
                params.push_str(&format!(",{},{flag}", self.de(*face)));
            }
            self.push(Pending {
                kind: 514,
                form: 1,
                transform: None,
                independent: false,
                label: String::new(),
                params,
            })
        };
        let shell_de = self.de(shell_entity);
        self.push(Pending {
            kind: 186,
            form: 0,
            transform: None,
            independent: true,
            label: label.chars().take(8).collect(),
            params: format!("{shell_de},1,0"),
        });
        Ok(())
    }

    fn face(&mut self, face: &Shape) -> OgeomResult<usize> {
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
        let surface_entity = self.surface(&surface)?;
        let wires = self.model.ordered_children_of(face)?;
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
                 as IGES — docs/PARITY.md, io.iges"
            );
        }
        let mut params = format!("{},{},1", self.de(surface_entity), loop_entities.len());
        for entity in &loop_entities {
            params.push_str(&format!(",{}", self.de(*entity)));
        }
        Ok(self.push(Pending {
            kind: 510,
            form: 1,
            transform: None,
            independent: false,
            label: String::new(),
            params,
        }))
    }

    /// A wire as a loop (508), or `None` when every edge in it is degenerate
    /// — a pole has no curve, and the reader rebuilds chart degeneracies
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
        let curve_entity = self.curve(&curve, range)?;
        let vertices = self.model.children_of(edge)?;
        let (from, to) = match vertices.len() {
            0 => ogeom_bail!(Construction, "an edge with no vertices cannot be written"),
            1 => (vertices[0].clone(), vertices[0].clone()),
            _ => (vertices[0].clone(), vertices[vertices.len() - 1].clone()),
        };
        // Each vertex carries its own composed placement — children_of has
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

    /// A curve over a range, in its IGES spelling.
    fn curve(&mut self, curve: &Curve, range: (f64, f64)) -> OgeomResult<usize> {
        match curve {
            Curve::Line(line) => {
                let a = line.point_at(range.0, self.tol)?;
                let b = line.point_at(range.1, self.tol)?;
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
                let transform = self.transform_entity(&circle.frame());
                let (s, e) = (range.0, range.1);
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
                let transform = self.transform_entity(&ellipse.frame());
                let at = |t: f64| (a * t.cos(), b * t.sin());
                let (sx, sy) = at(range.0);
                let (ex, ey) = at(range.1);
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
            Curve::Trimmed(t) => self.curve(t.basis(), range),
            Curve::BSpline(b) => self.nurbs_curve(b, range),
            other => {
                // The exact conversion carries anything with a closed NURBS
                // form; what has none — a helix — is refused by name there.
                let bspline = other.to_bspline_over(range, self.tol)?;
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
            S::Sphere(s) => {
                let sphere = s.sphere();
                let point = self.point_entity(sphere.centre());
                let params = format!("{},{}", self.de(point), fmt(sphere.radius()));
                Ok(self.push(Pending {
                    kind: 196,
                    form: 0,
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
            S::BSpline(b) => {
                let (uk, vk) = (b.u_knots(), b.v_knots());
                let grid = b.grid();
                let (nu, nv) = (grid.u_count(), grid.v_count());
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
                        let w = grid.get(u, v).map_or(1.0, |w| w.weight);
                        params.push_str(&format!(",{}", fmt(w)));
                    }
                }
                for v in 0..nv {
                    for u in 0..nu {
                        let p = grid
                            .get(u, v)
                            .map_or(Point::ORIGIN, ogeom_math::Weighted::point);
                        params.push_str(&format!(",{},{},{}", fmt(p.x), fmt(p.y), fmt(p.z)));
                    }
                }
                let ((u0, u1), (v0, v1)) = b.domain();
                params.push_str(&format!(",{},{},{},{}", fmt(u0), fmt(u1), fmt(v0), fmt(v1)));
                Ok(self.push(Pending {
                    kind: 128,
                    form: 0,
                    transform: None,
                    independent: false,
                    label: String::new(),
                    params,
                }))
            }
            other => {
                // Swept, offset and trimmed carriers convert exactly where a
                // closed NURBS form exists; the conversion refuses by name
                // where it does not.
                let bspline = other.to_bspline(self.tol)?;
                self.surface(&SurfaceGeometry::BSpline(bspline))
            }
        }
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

/// A rigid transform quantized to bits, for per-placement deduplication —
/// the same probe the STEP writer uses.
fn transform_bits(t: &Transform) -> [u64; 3] {
    let p = t.apply(Point::new(0.123_456_789, 9.87, -3.21));
    [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
}

/// A real in the file: decimal point kept, so a reader that types by
/// spelling reads a real.
fn fmt(v: f64) -> String {
    let s = format!("{v}");
    if s.contains('.') || s.contains('e') || s.contains('E') {
        s
    } else {
        format!("{s}.")
    }
}

/// Parameter text into 64-column lines, split at delimiters.
fn wrap_params(text: &str) -> Vec<String> {
    let mut lines = Vec::new();
    let mut current = String::new();
    for piece in text.split_inclusive(',') {
        if current.len() + piece.len() > 64 {
            lines.push(std::mem::take(&mut current));
        }
        current.push_str(piece);
    }
    if !current.is_empty() {
        lines.push(current);
    }
    lines
}
