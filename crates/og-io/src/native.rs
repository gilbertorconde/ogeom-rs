//! The native `.og` format: a whole document, written so it can be read back.
//!
//! Topology, geometry, placements, per-entity tolerances, provenance and the
//! cached tessellation — everything a [`Model`] holds. This is what the
//! differential harness round-trips through, which sets the bar: a format that
//! quietly dropped what it did not handle would make every comparison through
//! it pass, including the ones that should not.
//!
//! # Text, not binary
//!
//! The point of a round trip is to *disagree usefully* when it goes wrong, and
//! a disagreement you can only see through a hex dump is one nobody will read.
//! So: one record per line, tagged, in dependency order — datums, geometry,
//! entities, then topology. A model written, read and written again produces
//! the same bytes, so `diff` is the whole comparison tool.
//!
//! Floats are printed in Rust's shortest round-tripping form, which is exact:
//! parsing what this writes gives the same `f64` back, bit for bit. It is not a
//! lossy text approximation of a binary truth.
//!
//! # Handles are preserved, not renumbered
//!
//! A document's entity identities are the thing references are recorded
//! against (`docs/DATA_MODEL.md` §8). Rebuilding a model through the builders
//! would mint fresh ones and quietly invalidate every reference in every
//! document that pointed into it, so reading goes through
//! [`Model::from_parts`], which reproduces the arenas exactly and checks that
//! every handle resolves.
//!
//! # What it refuses
//!
//! A representation, curve, surface or node kind this writer does not know how
//! to write is an error, not a silent omission. That is the whole discipline of
//! the format: `write` either produces a file that reads back as the same
//! model, or it says it cannot.

use og_core::{
    EntityId, Key, OgResult, OpId, Provenance, Role, SourceId, Tolerance, Tolerances, og_bail,
};
use og_geom::{
    BSpline2d, BSplineCurve, BSplineSurface, Circle2d, CircleCurve, ConeSurface, Curve, Curve2d,
    Curve3d, CylinderSurface, Ellipse2d, EllipseCurve, ExtrusionSurface, HyperbolaCurve, Line2d,
    LineCurve, ParabolaCurve, PlanarCurve, PlaneSurface, RevolutionSurface, SphereSurface, Surface,
    SurfaceGeometry, TorusSurface, Trimmed2d, TrimmedCurve, TrimmedSurface,
};
use og_math::{
    Axis, Axis2, Circle, Circle2, Cone, ControlGrid, Cylinder, Direction, Direction2, Ellipse,
    Ellipse2, Frame, Frame2, Hyperbola, KnotVector, Parabola, Plane, Point, Point2, Sphere, Torus,
    Transform, Vector, Weighted,
};
use og_topo::{
    DatumId, EdgeData, EdgeRepr, FaceData, GeometryStore, Location, Model, ModelParts, NodeData,
    Orientation, Shape, ShapeType, TShape, TShapeId, Triangulation, VertexData,
};

/// The format version this reads and writes.
///
/// Bumped when the grammar changes in a way an older reader could not follow.
/// A file naming a version this does not know is refused rather than guessed
/// at.
pub const VERSION: u32 = 1;

/// The word every file starts with.
const MAGIC: &str = "openGeometry";

/// What to include when writing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WriteOptions {
    /// Whether to write the cached triangulations.
    ///
    /// They are derived data — recomputable from the topology at any deflection
    /// — and far larger than the topology that produced them, so a file meant
    /// for reading by a human is better without. Left on by default, because
    /// the alternative is a round trip that is not the identity for a model
    /// that has meshes attached, and this is the format the harness compares
    /// through.
    pub triangulations: bool,
}

impl Default for WriteOptions {
    fn default() -> Self {
        Self {
            triangulations: true,
        }
    }
}

/// Write a model, and the shapes it is a document *about*.
///
/// `roots` are recorded so a reader gets back the same handles the writer had.
/// A model with no roots is legal — it is a document of loose geometry — but a
/// reader then has nothing to start a traversal from.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the model holds
/// something this does not know how to write, which is refused rather than
/// dropped, or if a root does not resolve in the model.
pub fn write(model: &Model, roots: &[Shape], options: WriteOptions) -> OgResult<String> {
    for root in roots {
        if model.node(root).is_none() {
            og_bail!(Construction, "a root shape is not in this model");
        }
    }

    let mut out = String::new();
    out.push_str(&format!("{MAGIC} {VERSION}\n"));

    for (id, datum) in model.datums().iter() {
        let mut t = Vec::new();
        w(&mut t, "datum");
        key(&mut t, id);
        transform(&mut t, &datum)?;
        emit(&mut out, &t);
    }

    let geometry = model.geometry();
    for (id, c) in geometry.curves() {
        let mut t = Vec::new();
        w(&mut t, "curve");
        key(&mut t, id);
        curve(&mut t, c)?;
        emit(&mut out, &t);
    }
    for (id, c) in geometry.pcurves() {
        let mut t = Vec::new();
        w(&mut t, "pcurve");
        key(&mut t, id);
        pcurve(&mut t, c)?;
        emit(&mut out, &t);
    }
    for (id, s) in geometry.surfaces() {
        let mut t = Vec::new();
        w(&mut t, "surface");
        key(&mut t, id);
        surface(&mut t, s)?;
        emit(&mut out, &t);
    }
    if options.triangulations {
        for (id, mesh) in geometry.triangulations() {
            write_mesh(&mut out, id, mesh);
        }
    } else if geometry.triangulation_count() > 0 {
        // Said out loud in the file rather than left to be inferred from an
        // absence: a reader comparing two documents needs to know the meshes
        // were left out on purpose.
        out.push_str(&format!(
            "# {} triangulation(s) omitted by request\n",
            geometry.triangulation_count()
        ));
    }

    for (id, entry) in model.provenance().iter() {
        let mut t = Vec::new();
        w(&mut t, "entity");
        u(&mut t, id.get());
        provenance(&mut t, entry);
        emit(&mut out, &t);
    }

    for (id, node) in model.nodes() {
        write_node(&mut out, id, node, options)?;
    }
    for (node, entity) in model.identities() {
        let mut t = Vec::new();
        w(&mut t, "identity");
        key(&mut t, node);
        u(&mut t, entity.get());
        emit(&mut out, &t);
    }

    let mut t = Vec::new();
    w(&mut t, "operation");
    u(&mut t, u64::from(model.current_operation().0));
    emit(&mut out, &t);

    for root in roots {
        let mut t = Vec::new();
        w(&mut t, "root");
        shape(&mut t, root);
        emit(&mut out, &t);
    }
    Ok(out)
}

/// Read a model and its root shapes back.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the text is not
/// this format, names a version this does not know, or is malformed;
/// [`OgError::Dangling`](og_core::OgError::Dangling) if it names a handle that
/// is not in the file.
pub fn read(text: &str) -> OgResult<(Model, Vec<Shape>)> {
    let tol = Tolerances::millimetres();
    let mut cursor = Cursor::new(text);

    if cursor.word()? != MAGIC {
        og_bail!(Construction, "not an openGeometry document");
    }
    let version = cursor.count()?;
    if version != VERSION as usize {
        og_bail!(
            Construction,
            "document is version {version}; this reads version {VERSION}"
        );
    }

    let mut parts = ModelParts::default();
    let mut geometry = GeometryStore::new();
    let mut roots = Vec::new();

    while !cursor.done() {
        let tag = cursor.word()?.to_string();
        match tag.as_str() {
            "datum" => {
                let key = cursor.key()?;
                expect_datum(key, parts.datums.len())?;
                parts.datums.push(cursor.transform(tol)?);
            }
            "curve" => {
                let (index, _) = cursor.key()?;
                expect_index(index, geometry.counts().0)?;
                geometry.add_curve(cursor.curve(tol)?);
            }
            "pcurve" => {
                let (index, _) = cursor.key()?;
                expect_index(index, geometry.counts().1)?;
                geometry.add_pcurve(cursor.pcurve(tol)?);
            }
            "surface" => {
                let (index, _) = cursor.key()?;
                expect_index(index, geometry.counts().2)?;
                geometry.add_surface(cursor.surface(tol)?);
            }
            "mesh" => {
                let (index, _) = cursor.key()?;
                expect_index(index, geometry.triangulation_count())?;
                geometry.add_triangulation(cursor.mesh()?);
            }
            "entity" => {
                let id = cursor.count()?;
                if id != parts.provenance.len() + 1 {
                    og_bail!(
                        Construction,
                        "entities must be written in the order their identities \
                         were issued; expected {}, got {id}",
                        parts.provenance.len() + 1
                    );
                }
                parts.provenance.push(cursor.provenance()?);
            }
            "node" => {
                let (index, _) = cursor.key()?;
                expect_index(index, parts.nodes.len())?;
                parts.nodes.push(cursor.node()?);
            }
            "identity" => {
                let node = cursor.shape_key()?;
                let raw = cursor.count()?;
                let Some(entity) = EntityId::from_raw(raw as u64) else {
                    og_bail!(Construction, "identity 0 was never issued");
                };
                parts.identity.push((node, entity));
            }
            "operation" => parts.current_op = OpId(cursor.small()?),
            "root" => roots.push(cursor.shape()?),
            other => og_bail!(Construction, "unknown record `{other}`"),
        }
    }

    parts.geometry = geometry;
    let model = Model::from_parts(parts)?;
    // The roots were rebuilt alongside the rest but travel outside `ModelParts`,
    // so they are still unbound: they name no arena and resolve nowhere. Binding
    // them checks them too, which is why a file naming a root that is not there
    // fails here rather than in whatever first tries to use it.
    let roots = roots
        .iter()
        .map(|root| model.bind(root))
        .collect::<OgResult<Vec<_>>>()?;
    Ok((model, roots))
}

/// Refuse a record written out of arena order.
///
/// The order *is* the numbering — a reader replays inserts to reproduce the
/// handles — so a gap or a repeat would silently shift every later reference.
fn expect_index(index: u32, next: usize) -> OgResult<()> {
    if index as usize != next {
        og_bail!(
            Construction,
            "records must run in arena order; expected index {next}, got {index}"
        );
    }
    Ok(())
}

/// As [`expect_index`], for a datum, whose generation must be a fresh arena's.
fn expect_datum(key: (u32, u32), next: usize) -> OgResult<()> {
    expect_index(key.0, next)?;
    if key.1 != 0 {
        og_bail!(
            Construction,
            "a datum with generation {} cannot be rebuilt: a fresh arena hands \
             out generation 0, so the handle would not match what the file says",
            key.1
        );
    }
    Ok(())
}

// --- writing -----------------------------------------------------------------

/// Append a word.
fn w(t: &mut Vec<String>, s: &str) {
    t.push(s.to_string());
}

/// Append a number, in the shortest form that parses back to it exactly.
fn n(t: &mut Vec<String>, v: f64) {
    t.push(format!("{v:?}"));
}

/// Append an unsigned integer.
fn u(t: &mut Vec<String>, v: u64) {
    t.push(v.to_string());
}

/// Append an arena handle, as `index:generation`.
fn key<T>(t: &mut Vec<String>, id: Key<T>) {
    t.push(format!("{}:{}", id.index(), id.generation()));
}

/// Append a flag.
fn flag(t: &mut Vec<String>, v: bool) {
    t.push(if v { "1" } else { "0" }.to_string());
}

/// Finish a record.
fn emit(out: &mut String, t: &[String]) {
    out.push_str(&t.join(" "));
    out.push('\n');
}

fn point(t: &mut Vec<String>, p: Point) {
    n(t, p.x);
    n(t, p.y);
    n(t, p.z);
}

fn point2(t: &mut Vec<String>, p: Point2) {
    n(t, p.x);
    n(t, p.y);
}

fn vector(t: &mut Vec<String>, v: Vector) {
    n(t, v.x);
    n(t, v.y);
    n(t, v.z);
}

fn direction(t: &mut Vec<String>, d: Direction) {
    vector(t, d.vector());
}

fn direction2(t: &mut Vec<String>, d: Direction2) {
    n(t, d.vector().x);
    n(t, d.vector().y);
}

/// A frame, as origin then all three axes.
///
/// `y` is written even though a right-handed frame derives it, because a
/// mirrored frame does not: writing only `z` and `x` would quietly turn every
/// left-handed frame right-handed on the way back.
fn frame(t: &mut Vec<String>, f: &Frame) {
    point(t, f.origin());
    direction(t, f.x());
    direction(t, f.y());
    direction(t, f.z());
}

fn frame2(t: &mut Vec<String>, f: &Frame2) {
    point2(t, f.origin());
    direction2(t, f.x());
    direction2(t, f.y());
}

fn axis(t: &mut Vec<String>, a: Axis) {
    point(t, a.location);
    direction(t, a.direction);
}

fn axis2(t: &mut Vec<String>, a: Axis2) {
    point2(t, a.location);
    direction2(t, a.direction);
}

fn range(t: &mut Vec<String>, r: (f64, f64)) {
    n(t, r.0);
    n(t, r.1);
}

/// A similarity transform: its unit linear part, its scale, its translation.
fn transform(t: &mut Vec<String>, x: &Transform) -> OgResult<()> {
    let m = x.linear();
    for row in 0..3 {
        for column in 0..3 {
            n(t, m.get(row, column)?);
        }
    }
    n(t, x.scale_factor());
    vector(t, x.translation_vector());
    Ok(())
}

/// A placement chain, or `-` for the identity.
fn location(t: &mut Vec<String>, l: &Location) {
    if l.is_identity() {
        w(t, "-");
        return;
    }
    let chain: Vec<String> = l
        .chain()
        .iter()
        .map(|(datum, power)| format!("{}:{}^{power}", datum.index(), datum.generation()))
        .collect();
    w(t, &chain.join(","));
}

/// A shape triple, as one token: node, orientation, placement.
fn shape(t: &mut Vec<String>, s: &Shape) {
    let orientation = match s.orientation() {
        Orientation::Forward => "F",
        Orientation::Reversed => "R",
        Orientation::Internal => "I",
        Orientation::External => "E",
    };
    let mut token = format!(
        "{}:{}/{orientation}",
        s.node().index(),
        s.node().generation()
    );
    if !s.location().is_identity() {
        let chain: Vec<String> = s
            .location()
            .chain()
            .iter()
            .map(|(datum, power)| format!("{}:{}^{power}", datum.index(), datum.generation()))
            .collect();
        token.push('/');
        token.push_str(&chain.join(","));
    }
    w(t, &token);
}

fn knots(t: &mut Vec<String>, k: &KnotVector) {
    u(t, k.degree() as u64);
    u(t, k.knots().len() as u64);
    for value in k.knots() {
        n(t, *value);
    }
}

fn weighted(t: &mut Vec<String>, c: Weighted<Point>) {
    // Written in the homogeneous form it is stored in, so no multiply and
    // divide stands between what was held and what comes back.
    point(t, c.scaled);
    n(t, c.weight);
}

fn weighted2(t: &mut Vec<String>, c: Weighted<Point2>) {
    point2(t, c.scaled);
    n(t, c.weight);
}

fn curve(t: &mut Vec<String>, c: &Curve) -> OgResult<()> {
    match c {
        Curve::Line(l) => {
            w(t, "line");
            axis(t, l.axis());
            range(t, l.domain());
        }
        Curve::Circle(c) => {
            w(t, "circle");
            frame(t, &c.circle().frame());
            n(t, c.circle().radius());
            flag(t, c.is_reversed());
        }
        Curve::Ellipse(e) => {
            w(t, "ellipse");
            frame(t, &e.ellipse().frame());
            n(t, e.ellipse().major_radius());
            n(t, e.ellipse().minor_radius());
            flag(t, e.is_reversed());
        }
        Curve::Hyperbola(h) => {
            w(t, "hyperbola");
            frame(t, &h.hyperbola().frame());
            n(t, h.hyperbola().major_radius());
            n(t, h.hyperbola().minor_radius());
            range(t, h.domain());
            flag(t, h.is_reversed());
        }
        Curve::Parabola(p) => {
            w(t, "parabola");
            frame(t, &p.parabola().frame());
            n(t, p.parabola().focal());
            range(t, p.domain());
            flag(t, p.is_reversed());
        }
        Curve::BSpline(b) => {
            if b.is_periodic() {
                og_bail!(
                    Construction,
                    "a periodic B-spline curve cannot be written: nothing \
                     builds one today, so the format has no way to say so, and \
                     dropping the flag would give back a different curve"
                );
            }
            w(t, "bspline");
            knots(t, b.knots());
            u(t, b.control_points().len() as u64);
            for c in b.control_points() {
                weighted(t, *c);
            }
        }
        Curve::Trimmed(x) => {
            w(t, "trimmed");
            range(t, x.domain());
            flag(t, x.is_reversed());
            curve(t, x.basis())?;
        }
    }
    Ok(())
}

fn pcurve(t: &mut Vec<String>, c: &PlanarCurve) -> OgResult<()> {
    match c {
        PlanarCurve::Line(l) => {
            w(t, "line2");
            axis2(t, l.axis());
            range(t, l.domain());
        }
        PlanarCurve::Circle(c) => {
            w(t, "circle2");
            frame2(t, &c.circle().frame());
            n(t, c.circle().radius());
            flag(t, c.is_reversed());
        }
        PlanarCurve::Ellipse(e) => {
            w(t, "ellipse2");
            frame2(t, &e.ellipse().frame());
            n(t, e.ellipse().major_radius());
            n(t, e.ellipse().minor_radius());
            flag(t, e.is_reversed());
        }
        PlanarCurve::BSpline(b) => {
            w(t, "bspline2");
            knots(t, b.knots());
            u(t, b.control_points().len() as u64);
            for c in b.control_points() {
                weighted2(t, *c);
            }
        }
        PlanarCurve::Trimmed(x) => {
            w(t, "trimmed2");
            range(t, x.domain());
            flag(t, x.is_reversed());
            pcurve(t, x.basis())?;
        }
    }
    Ok(())
}

fn surface(t: &mut Vec<String>, s: &SurfaceGeometry) -> OgResult<()> {
    let (u_domain, v_domain) = s.domain();
    match s {
        SurfaceGeometry::Plane(p) => {
            w(t, "plane");
            frame(t, &p.plane().frame());
            range(t, u_domain);
            range(t, v_domain);
        }
        SurfaceGeometry::Cylinder(c) => {
            w(t, "cylinder");
            frame(t, &c.cylinder().frame());
            n(t, c.cylinder().radius());
            range(t, v_domain);
        }
        SurfaceGeometry::Cone(c) => {
            w(t, "cone");
            frame(t, &c.cone().frame());
            n(t, c.cone().reference_radius());
            n(t, c.cone().half_angle());
            range(t, v_domain);
        }
        SurfaceGeometry::Sphere(s) => {
            w(t, "sphere");
            frame(t, &s.sphere().frame());
            n(t, s.sphere().radius());
        }
        SurfaceGeometry::Torus(x) => {
            w(t, "torus");
            frame(t, &x.torus().frame());
            n(t, x.torus().major_radius());
            n(t, x.torus().minor_radius());
        }
        SurfaceGeometry::BSpline(b) => {
            w(t, "bsurface");
            knots(t, b.u_knots());
            knots(t, b.v_knots());
            u(t, b.grid().u_count() as u64);
            u(t, b.grid().v_count() as u64);
            for c in b.grid().points() {
                weighted(t, *c);
            }
        }
        SurfaceGeometry::Revolution(r) => {
            w(t, "revolution");
            axis(t, r.axis());
            range(t, u_domain);
            curve(t, r.curve())?;
        }
        SurfaceGeometry::Extrusion(e) => {
            w(t, "extrusion");
            direction(t, e.direction());
            range(t, v_domain);
            curve(t, e.curve())?;
        }
        SurfaceGeometry::Trimmed(x) => {
            w(t, "tsurface");
            range(t, u_domain);
            range(t, v_domain);
            surface(t, x.basis())?;
        }
    }
    Ok(())
}

fn provenance(t: &mut Vec<String>, p: &Provenance) {
    match p {
        Provenance::Primitive { op, role } => {
            w(t, "primitive");
            u(t, u64::from(op.0));
            u(t, u64::from(role.0));
        }
        Provenance::Derived { op, from, role } => {
            w(t, "derived");
            u(t, u64::from(op.0));
            u(t, u64::from(role.0));
            u(t, from.len() as u64);
            for source in from {
                u(t, source.get());
            }
        }
        Provenance::Imported { source, external } => {
            w(t, "imported");
            u(t, u64::from(source.0));
            u(t, *external);
        }
    }
}

/// A cached mesh: a header line, then a line per vertex and a line per triangle.
///
/// Broken across lines rather than run together, because a mesh is the largest
/// thing in the file by far and a one-line-per-vertex diff is the difference
/// between reading a disagreement and scrolling past it.
fn write_mesh(out: &mut String, id: og_topo::TriangulationId, mesh: &Triangulation) {
    let mut t = Vec::new();
    w(&mut t, "mesh");
    key(&mut t, id);
    u(&mut t, mesh.positions.len() as u64);
    u(&mut t, mesh.triangles.len() as u64);
    flag(&mut t, mesh.deflection_met);
    emit(out, &t);

    for i in 0..mesh.positions.len() {
        let mut t = Vec::new();
        w(&mut t, "v");
        point(&mut t, mesh.positions[i]);
        // A mesh may carry fewer normals or parameters than positions if it was
        // assembled by hand; the zeros keep the record rectangular and the
        // reader keeps whatever length it is told.
        vector(&mut t, mesh.normals.get(i).copied().unwrap_or(Vector::ZERO));
        let (u_at, v_at) = mesh.parameters.get(i).copied().unwrap_or((0.0, 0.0));
        n(&mut t, u_at);
        n(&mut t, v_at);
        emit(out, &t);
    }
    for triangle in &mesh.triangles {
        let mut t = Vec::new();
        w(&mut t, "f");
        for index in triangle {
            u(&mut t, u64::from(*index));
        }
        emit(out, &t);
    }
}

fn write_node(
    out: &mut String,
    id: TShapeId,
    node: &TShape,
    options: WriteOptions,
) -> OgResult<()> {
    let mut t = Vec::new();
    w(&mut t, "node");
    key(&mut t, id);
    w(
        &mut t,
        match node.kind() {
            ShapeType::Vertex => "vertex",
            ShapeType::Edge => "edge",
            ShapeType::Wire => "wire",
            ShapeType::Face => "face",
            ShapeType::Shell => "shell",
            ShapeType::Solid => "solid",
            ShapeType::CompSolid => "compsolid",
            ShapeType::Compound => "compound",
        },
    );

    let mut representations = Vec::new();
    match node.data() {
        NodeData::Vertex(v) => {
            n(&mut t, v.tolerance.get());
            point(&mut t, v.point);
        }
        NodeData::Edge(e) => {
            n(&mut t, e.tolerance.get());
            flag(&mut t, e.same_parameter());
            flag(&mut t, e.degenerate);
            u(&mut t, e.representations.len() as u64);
            for repr in &e.representations {
                let mut line = Vec::new();
                w(&mut line, "r");
                write_repr(&mut line, repr)?;
                representations.push(line);
            }
        }
        NodeData::Face(f) => {
            n(&mut t, f.tolerance.get());
            key(&mut t, f.surface);
            location(&mut t, &f.location);
            flag(&mut t, f.natural_restriction);
            match f.triangulation.filter(|_| options.triangulations) {
                Some(mesh) => key(&mut t, mesh),
                None => w(&mut t, "-"),
            }
        }
        NodeData::Container => {}
    }

    u(&mut t, node.children().len() as u64);
    for child in node.children() {
        shape(&mut t, child);
    }
    emit(out, &t);
    for line in &representations {
        emit(out, line);
    }
    Ok(())
}

fn write_repr(t: &mut Vec<String>, repr: &EdgeRepr) -> OgResult<()> {
    match repr {
        EdgeRepr::Curve3d {
            curve,
            location: at,
            range: r,
        } => {
            w(t, "curve3d");
            key(t, *curve);
            location(t, at);
            range(t, *r);
        }
        EdgeRepr::PCurve {
            curve,
            surface,
            location: at,
            range: r,
        } => {
            w(t, "pcurve");
            key(t, *curve);
            key(t, *surface);
            location(t, at);
            range(t, *r);
        }
        EdgeRepr::Seam {
            forward,
            reversed,
            surface,
            location: at,
            range: r,
        } => {
            w(t, "seam");
            key(t, *forward);
            key(t, *reversed);
            key(t, *surface);
            location(t, at);
            range(t, *r);
        }
        EdgeRepr::Polyline {
            points,
            parameters,
            location: at,
            deflection,
        } => {
            w(t, "polyline");
            location(t, at);
            n(t, *deflection);
            u(t, points.len() as u64);
            for p in points {
                point(t, *p);
            }
            u(t, parameters.len() as u64);
            for at in parameters {
                n(t, *at);
            }
        }
        // `EdgeRepr` is non-exhaustive, so a variant added later lands here.
        // Refused, not skipped: a file missing one of an edge's descriptions of
        // itself reads back as an edge that has lost a face's pcurve, and
        // nothing downstream would say where it went.
        other => og_bail!(
            Construction,
            "this version writes no edge representation of that kind: {other:?}"
        ),
    }
    Ok(())
}

// --- reading -----------------------------------------------------------------

/// A position in the token stream.
///
/// The file is written a record to a line, but read as a flat stream of
/// whitespace-separated words: layout is for the person reading a diff, and
/// making the parser depend on it would only add a second thing to get wrong.
struct Cursor<'a> {
    items: Vec<&'a str>,
    at: usize,
}

impl<'a> Cursor<'a> {
    fn new(text: &'a str) -> Self {
        let items = text
            .lines()
            .map(|line| line.split('#').next().unwrap_or(""))
            .flat_map(str::split_whitespace)
            .collect();
        Self { items, at: 0 }
    }

    fn done(&self) -> bool {
        self.at >= self.items.len()
    }

    fn word(&mut self) -> OgResult<&'a str> {
        let Some(item) = self.items.get(self.at) else {
            og_bail!(Construction, "the document ends part-way through a record");
        };
        self.at += 1;
        Ok(item)
    }

    fn number(&mut self) -> OgResult<f64> {
        let word = self.word()?;
        word.parse::<f64>()
            .map_err(|_| og_core::og_err!(Construction, "`{word}` is not a number"))
    }

    fn count(&mut self) -> OgResult<usize> {
        let word = self.word()?;
        word.parse::<usize>()
            .map_err(|_| og_core::og_err!(Construction, "`{word}` is not a count"))
    }

    fn small(&mut self) -> OgResult<u32> {
        let word = self.word()?;
        word.parse::<u32>()
            .map_err(|_| og_core::og_err!(Construction, "`{word}` is not an identifier"))
    }

    fn flag(&mut self) -> OgResult<bool> {
        match self.word()? {
            "0" => Ok(false),
            "1" => Ok(true),
            other => og_bail!(Construction, "`{other}` is not a flag"),
        }
    }

    /// An `index:generation` pair.
    fn key(&mut self) -> OgResult<(u32, u32)> {
        parse_key(self.word()?)
    }

    fn handle<T>(&mut self) -> OgResult<Key<T>> {
        let (index, generation) = self.key()?;
        Ok(Key::from_parts(index, generation))
    }

    fn shape_key(&mut self) -> OgResult<TShapeId> {
        self.handle()
    }

    fn point(&mut self) -> OgResult<Point> {
        Ok(Point::new(self.number()?, self.number()?, self.number()?))
    }

    fn point2(&mut self) -> OgResult<Point2> {
        Ok(Point2::new(self.number()?, self.number()?))
    }

    fn vector(&mut self) -> OgResult<Vector> {
        Ok(Vector::new(self.number()?, self.number()?, self.number()?))
    }

    /// A direction, checked rather than renormalized.
    ///
    /// Normalizing what is already a unit vector moves it by a bit or two, and
    /// that is enough to make writing what was read produce a different file.
    fn direction(&mut self, tol: Tolerances) -> OgResult<Direction> {
        Direction::unit(self.vector()?, tol)
    }

    fn direction2(&mut self, tol: Tolerances) -> OgResult<Direction2> {
        Direction2::unit(og_math::Vector2::new(self.number()?, self.number()?), tol)
    }

    fn frame(&mut self, tol: Tolerances) -> OgResult<Frame> {
        let origin = self.point()?;
        let x = self.direction(tol)?;
        let y = self.direction(tol)?;
        let z = self.direction(tol)?;
        Frame::from_axes(origin, x, y, z, tol)
    }

    fn frame2(&mut self, tol: Tolerances) -> OgResult<Frame2> {
        let origin = self.point2()?;
        let x = self.direction2(tol)?;
        let y = self.direction2(tol)?;
        Frame2::from_axes(origin, x, y, tol)
    }

    fn axis(&mut self, tol: Tolerances) -> OgResult<Axis> {
        Ok(Axis::new(self.point()?, self.direction(tol)?))
    }

    fn axis2(&mut self, tol: Tolerances) -> OgResult<Axis2> {
        Ok(Axis2::new(self.point2()?, self.direction2(tol)?))
    }

    fn range(&mut self) -> OgResult<(f64, f64)> {
        Ok((self.number()?, self.number()?))
    }

    fn transform(&mut self, tol: Tolerances) -> OgResult<Transform> {
        let mut rows = [[0.0_f64; 3]; 3];
        for row in &mut rows {
            for cell in row.iter_mut() {
                *cell = self.number()?;
            }
        }
        let scale = self.number()?;
        let translation = self.vector()?;
        // Put back together from the parts a transform is stored as, rather
        // than multiplied into one matrix and factored apart again: the round
        // trip through a scale that divides out is not exact, and a placement
        // that drifts a little on every save is worse than one that fails.
        // The kind is re-derived, since it is a function of the other three.
        Transform::from_parts(
            og_math::Matrix3::new(rows),
            scale,
            translation,
            tol.angular().max(1e-9),
        )
    }

    fn location(&mut self) -> OgResult<Location> {
        let word = self.word()?;
        parse_location(word)
    }

    fn shape(&mut self) -> OgResult<Shape> {
        let word = self.word()?;
        let mut parts = word.split('/');
        let (Some(node), Some(orientation)) = (parts.next(), parts.next()) else {
            og_bail!(Construction, "`{word}` is not a shape");
        };
        let (index, generation) = parse_key(node)?;
        let orientation = match orientation {
            "F" => Orientation::Forward,
            "R" => Orientation::Reversed,
            "I" => Orientation::Internal,
            "E" => Orientation::External,
            other => og_bail!(Construction, "`{other}` is not an orientation"),
        };
        let location = match parts.next() {
            Some(chain) => parse_location(chain)?,
            None => Location::identity(),
        };
        if parts.next().is_some() {
            og_bail!(Construction, "`{word}` has more parts than a shape has");
        }
        Ok(Shape::new(
            Key::from_parts(index, generation),
            location,
            orientation,
        ))
    }

    fn knots(&mut self) -> OgResult<KnotVector> {
        let degree = self.count()?;
        let n = self.count()?;
        let mut values = Vec::with_capacity(n);
        for _ in 0..n {
            values.push(self.number()?);
        }
        KnotVector::new(values, degree)
    }

    fn weighted(&mut self) -> OgResult<Weighted<Point>> {
        Ok(Weighted {
            scaled: self.point()?,
            weight: self.number()?,
        })
    }

    fn weighted2(&mut self) -> OgResult<Weighted<Point2>> {
        Ok(Weighted {
            scaled: self.point2()?,
            weight: self.number()?,
        })
    }

    fn curve(&mut self, tol: Tolerances) -> OgResult<Curve> {
        Ok(match self.word()? {
            "line" => {
                let axis = self.axis(tol)?;
                let (lo, hi) = self.range()?;
                LineCurve::over(axis, lo, hi)?.into()
            }
            "circle" => {
                let circle = Circle::new(self.frame(tol)?, self.number()?, tol)?;
                reverse_if(CircleCurve::new(circle).into(), self.flag()?)
            }
            "ellipse" => {
                let frame = self.frame(tol)?;
                let ellipse = Ellipse::new(frame, self.number()?, self.number()?, tol)?;
                reverse_if(EllipseCurve::new(ellipse).into(), self.flag()?)
            }
            "hyperbola" => {
                let frame = self.frame(tol)?;
                let h = Hyperbola::new(frame, self.number()?, self.number()?, tol)?;
                let (lo, hi) = self.range()?;
                let curve = HyperbolaCurve::new(h, extent(lo, hi)?)?;
                reverse_if(curve.into(), self.flag()?)
            }
            "parabola" => {
                let frame = self.frame(tol)?;
                let p = Parabola::new(frame, self.number()?, tol)?;
                let (lo, hi) = self.range()?;
                let curve = ParabolaCurve::new(p, extent(lo, hi)?)?;
                reverse_if(curve.into(), self.flag()?)
            }
            "bspline" => {
                let knots = self.knots()?;
                let n = self.count()?;
                let mut control = Vec::with_capacity(n);
                for _ in 0..n {
                    control.push(self.weighted()?);
                }
                BSplineCurve::rational(knots, control)?.into()
            }
            "trimmed" => {
                let (lo, hi) = self.range()?;
                let reversed = self.flag()?;
                let basis = self.curve(tol)?;
                let curve = TrimmedCurve::new(basis, lo, hi, tol)?;
                reverse_if(Curve::Trimmed(Box::new(curve)), reversed)
            }
            other => og_bail!(Construction, "`{other}` is not a curve this reads"),
        })
    }

    fn pcurve(&mut self, tol: Tolerances) -> OgResult<PlanarCurve> {
        Ok(match self.word()? {
            "line2" => {
                let axis = self.axis2(tol)?;
                let (lo, hi) = self.range()?;
                Line2d::over(axis, lo, hi)?.into()
            }
            "circle2" => {
                let circle = Circle2::new(self.frame2(tol)?, self.number()?, tol)?;
                reverse_if(Circle2d::new(circle).into(), self.flag()?)
            }
            "ellipse2" => {
                let frame = self.frame2(tol)?;
                let ellipse = Ellipse2::new(frame, self.number()?, self.number()?, tol)?;
                reverse_if(Ellipse2d::new(ellipse).into(), self.flag()?)
            }
            "bspline2" => {
                let knots = self.knots()?;
                let n = self.count()?;
                let mut control = Vec::with_capacity(n);
                for _ in 0..n {
                    control.push(self.weighted2()?);
                }
                BSpline2d::rational(knots, control)?.into()
            }
            "trimmed2" => {
                let (lo, hi) = self.range()?;
                let reversed = self.flag()?;
                let basis = self.pcurve(tol)?;
                let curve = Trimmed2d::new(basis, lo, hi, tol)?;
                reverse_if(PlanarCurve::Trimmed(Box::new(curve)), reversed)
            }
            other => og_bail!(Construction, "`{other}` is not a planar curve this reads"),
        })
    }

    fn surface(&mut self, tol: Tolerances) -> OgResult<SurfaceGeometry> {
        Ok(match self.word()? {
            "plane" => {
                let plane = Plane::new(self.frame(tol)?);
                PlaneSurface::over(plane, self.range()?, self.range()?)?.into()
            }
            "cylinder" => {
                let cylinder = Cylinder::new(self.frame(tol)?, self.number()?, tol)?;
                CylinderSurface::new(cylinder, self.range()?)?.into()
            }
            "cone" => {
                let frame = self.frame(tol)?;
                let cone = Cone::new(frame, self.number()?, self.number()?, tol)?;
                ConeSurface::new(cone, self.range()?)?.into()
            }
            "sphere" => {
                let sphere = Sphere::new(self.frame(tol)?, self.number()?, tol)?;
                SphereSurface::new(sphere).into()
            }
            "torus" => {
                let frame = self.frame(tol)?;
                let torus = Torus::new(frame, self.number()?, self.number()?, tol)?;
                TorusSurface::new(torus).into()
            }
            "bsurface" => {
                let u_knots = self.knots()?;
                let v_knots = self.knots()?;
                let u_count = self.count()?;
                let v_count = self.count()?;
                let mut points = Vec::with_capacity(u_count * v_count);
                for _ in 0..u_count * v_count {
                    points.push(self.weighted()?);
                }
                let grid = ControlGrid::new(points, u_count, v_count)?;
                BSplineSurface::rational(u_knots, v_knots, grid)?.into()
            }
            "revolution" => {
                let axis = self.axis(tol)?;
                let (start, end) = self.range()?;
                if start != 0.0 {
                    og_bail!(
                        Construction,
                        "a revolution's angle starts at zero; this one starts at \
                         {start}"
                    );
                }
                let basis = self.curve(tol)?;
                RevolutionSurface::new(basis, axis, end)?.into()
            }
            "extrusion" => {
                let direction = self.direction(tol)?;
                let (start, end) = self.range()?;
                if start != 0.0 {
                    og_bail!(
                        Construction,
                        "an extrusion's extent starts at zero; this one starts \
                         at {start}"
                    );
                }
                let basis = self.curve(tol)?;
                ExtrusionSurface::new(basis, direction, end)?.into()
            }
            "tsurface" => {
                let u_range = self.range()?;
                let v_range = self.range()?;
                let basis = self.surface(tol)?;
                SurfaceGeometry::Trimmed(Box::new(TrimmedSurface::new(
                    basis, u_range, v_range, tol,
                )?))
            }
            other => og_bail!(Construction, "`{other}` is not a surface this reads"),
        })
    }

    fn mesh(&mut self) -> OgResult<Triangulation> {
        let vertices = self.count()?;
        let triangles = self.count()?;
        let mut mesh = Triangulation::new();
        mesh.deflection_met = self.flag()?;
        for _ in 0..vertices {
            if self.word()? != "v" {
                og_bail!(Construction, "expected a mesh vertex");
            }
            mesh.positions.push(self.point()?);
            mesh.normals.push(self.vector()?);
            mesh.parameters.push((self.number()?, self.number()?));
        }
        for _ in 0..triangles {
            if self.word()? != "f" {
                og_bail!(Construction, "expected a mesh triangle");
            }
            let mut corners = [0_u32; 3];
            for corner in &mut corners {
                *corner = self.small()?;
                if *corner as usize >= vertices {
                    og_bail!(
                        Dangling,
                        "a triangle names vertex {corner}, and the mesh has \
                         {vertices}"
                    );
                }
            }
            mesh.triangles.push(corners);
        }
        Ok(mesh)
    }

    fn provenance(&mut self) -> OgResult<Provenance> {
        Ok(match self.word()? {
            "primitive" => Provenance::Primitive {
                op: OpId(self.small()?),
                role: Role(self.small()?),
            },
            "derived" => {
                let op = OpId(self.small()?);
                let role = Role(self.small()?);
                let n = self.count()?;
                let mut from = Vec::with_capacity(n);
                for _ in 0..n {
                    let raw = self.count()?;
                    let Some(id) = EntityId::from_raw(raw as u64) else {
                        og_bail!(Construction, "identity 0 was never issued");
                    };
                    from.push(id);
                }
                Provenance::Derived {
                    op,
                    from: from.into_iter().collect(),
                    role,
                }
            }
            "imported" => Provenance::Imported {
                source: SourceId(self.small()?),
                external: self.count()? as u64,
            },
            other => og_bail!(Construction, "`{other}` is not a provenance"),
        })
    }

    fn node(&mut self) -> OgResult<TShape> {
        let kind = match self.word()? {
            "vertex" => ShapeType::Vertex,
            "edge" => ShapeType::Edge,
            "wire" => ShapeType::Wire,
            "face" => ShapeType::Face,
            "shell" => ShapeType::Shell,
            "solid" => ShapeType::Solid,
            "compsolid" => ShapeType::CompSolid,
            "compound" => ShapeType::Compound,
            other => og_bail!(Construction, "`{other}` is not a kind of shape"),
        };

        // An edge's representations are records of their own, written after the
        // node line — so its data cannot be finished until the children on that
        // line have been read past.
        let mut pending = None;
        let mut data = match kind {
            ShapeType::Vertex => {
                let tolerance = self.number()?;
                NodeData::Vertex(VertexData::with_tolerance(self.point()?, tolerance)?)
            }
            ShapeType::Edge => {
                let mut edge = EdgeData::new();
                edge.tolerance = Tolerance::new(self.number()?)?;
                let agrees = self.flag()?;
                edge.degenerate = self.flag()?;
                pending = Some((agrees, self.count()?));
                NodeData::Edge(Box::new(edge))
            }
            ShapeType::Face => {
                let tolerance = Tolerance::new(self.number()?)?;
                let surface = self.handle()?;
                let at = self.location()?;
                let natural = self.flag()?;
                let mut face = if natural {
                    FaceData::natural(surface, at)
                } else {
                    FaceData::new(surface, at)
                };
                face.tolerance = tolerance;
                face.triangulation = match self.word()? {
                    "-" => None,
                    word => Some(Key::from_parts_of(parse_key(word)?)),
                };
                NodeData::Face(Box::new(face))
            }
            _ => NodeData::Container,
        };

        let n = self.count()?;
        let mut children = Vec::with_capacity(n);
        for _ in 0..n {
            children.push(self.shape()?);
        }

        if let (Some((agrees, count)), NodeData::Edge(edge)) = (pending, &mut data) {
            for _ in 0..count {
                if self.word()? != "r" {
                    og_bail!(Construction, "expected an edge representation");
                }
                let repr = self.repr()?;
                edge.add(repr);
            }
            // Set last: `add` clears it, which is the right default for a
            // caller attaching a representation and the wrong one for a reader
            // restoring a claim the document already makes.
            edge.assert_same_parameter(agrees);
        }
        Ok(TShape::new(kind, data, children))
    }

    fn repr(&mut self) -> OgResult<EdgeRepr> {
        Ok(match self.word()? {
            "curve3d" => EdgeRepr::Curve3d {
                curve: self.handle()?,
                location: self.location()?,
                range: self.range()?,
            },
            "pcurve" => EdgeRepr::PCurve {
                curve: self.handle()?,
                surface: self.handle()?,
                location: self.location()?,
                range: self.range()?,
            },
            "seam" => EdgeRepr::Seam {
                forward: self.handle()?,
                reversed: self.handle()?,
                surface: self.handle()?,
                location: self.location()?,
                range: self.range()?,
            },
            "polyline" => {
                let location = self.location()?;
                let deflection = self.number()?;
                let n = self.count()?;
                let mut points = Vec::with_capacity(n);
                for _ in 0..n {
                    points.push(self.point()?);
                }
                let n = self.count()?;
                let mut parameters = Vec::with_capacity(n);
                for _ in 0..n {
                    parameters.push(self.number()?);
                }
                EdgeRepr::Polyline {
                    points,
                    parameters,
                    location,
                    deflection,
                }
            }
            other => og_bail!(
                Construction,
                "`{other}` is not an edge representation this reads"
            ),
        })
    }
}

/// A trimmed curve's half-width, from the symmetric domain it was written as.
fn extent(lo: f64, hi: f64) -> OgResult<f64> {
    if (lo + hi).abs() > 1e-12 * hi.abs().max(1.0) {
        og_bail!(
            Construction,
            "a hyperbola or parabola spans [-extent, extent]; [{lo}, {hi}] is \
             not symmetric, and this version has no way to express it"
        );
    }
    Ok(hi)
}

/// Reverse a curve if the file says it runs backwards.
fn reverse_if<T: og_geom::Reversible>(curve: T, reversed: bool) -> T {
    if reversed { curve.reversed() } else { curve }
}

fn parse_key(word: &str) -> OgResult<(u32, u32)> {
    let Some((index, generation)) = word.split_once(':') else {
        og_bail!(Construction, "`{word}` is not a handle");
    };
    let (Ok(index), Ok(generation)) = (index.parse::<u32>(), generation.parse::<u32>()) else {
        og_bail!(Construction, "`{word}` is not a handle");
    };
    Ok((index, generation))
}

fn parse_location(word: &str) -> OgResult<Location> {
    if word == "-" {
        return Ok(Location::identity());
    }
    let mut location = Location::identity();
    for step in word.split(',') {
        let Some((handle, power)) = step.rsplit_once('^') else {
            og_bail!(Construction, "`{step}` is not a placement step");
        };
        let (index, generation) = parse_key(handle)?;
        let Ok(power) = power.parse::<i32>() else {
            og_bail!(Construction, "`{power}` is not a power");
        };
        let datum: DatumId = Key::from_parts(index, generation);
        location = location.then(&Location::powered(datum, power));
    }
    Ok(location)
}

/// A handle from an already-parsed `(index, generation)` pair.
trait FromParsedKey: Sized {
    fn from_parts_of(parts: (u32, u32)) -> Self;
}

impl<T> FromParsedKey for Key<T> {
    fn from_parts_of((index, generation): (u32, u32)) -> Self {
        Self::from_parts(index, generation)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use og_algo::{
        check, check_tessellation, make_box, make_cone, make_cylinder, make_prism, make_revolution,
        make_sphere, make_torus, make_wedge,
    };
    use og_math::{Frame, Vector};
    use og_mesh::{Deflection, triangulate};
    use og_topo::explore_unique;

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> Deflection {
        Deflection {
            chord: 0.02,
            ..Deflection::default()
        }
    }

    /// One model holding every kind of shape this can build, and the roots that
    /// name them.
    fn everything() -> (Model, Vec<Shape>) {
        let mut model = Model::new();
        let mut roots = vec![
            make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T)
                .unwrap()
                .shape,
            make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T)
                .unwrap()
                .shape,
            make_sphere(&mut model, Frame::WORLD, 3.0, T).unwrap().shape,
            make_cone(&mut model, Frame::WORLD, 3.0, 1.0, 4.0, T)
                .unwrap()
                .shape,
            make_cone(&mut model, Frame::WORLD, 3.0, 0.0, 4.0, T)
                .unwrap()
                .shape,
            make_torus(&mut model, Frame::WORLD, 5.0, 2.0, T)
                .unwrap()
                .shape,
            make_wedge(&mut model, Frame::WORLD, (4.0, 4.0, 6.0), (2.0, 2.0), T)
                .unwrap()
                .shape,
        ];

        // A prism, so an extrusion surface and a placement datum are in there.
        let face = explore_unique(&model, &roots[0], ShapeType::Face).unwrap()[0].clone();
        roots.push(
            make_prism(&mut model, &face, Vector::new(0.0, 0.0, 2.0), T)
                .unwrap()
                .shape,
        );
        // And a revolution, so a surface of revolution and a seam are too. The
        // profile has to lie in a plane the axis runs through and stay clear of
        // it, or the swept solid would pass through itself — so it is the box's
        // `-y` face, spanning x and z, turned about a z axis five units away.
        let side = explore_unique(&model, &roots[0], ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                model.provenance_of(f).and_then(og_core::Provenance::role)
                    == Some(og_algo::primitive::roles::FACE_MIN_Y)
            })
            .expect("the box has a -y face");
        let axis = og_math::Axis::new(og_math::Point::new(-5.0, 0.0, 0.0), og_math::Direction::Z);
        roots.push(
            make_revolution(&mut model, &side, axis, 1.0, T)
                .unwrap()
                .shape,
        );
        (model, roots)
    }

    #[test]
    fn a_file_that_is_not_one_is_refused() {
        assert!(read("").is_err());
        assert!(read("something else 1").is_err());
        assert!(read(&format!("{MAGIC} 99")).is_err());
        assert!(read(&format!("{MAGIC} 1\nnonsense\n")).is_err());
    }

    #[test]
    fn writing_what_was_read_gives_the_same_bytes() {
        // The property the whole format exists for. If a round trip through it
        // is not the identity, `diff` says exactly where — which is why this is
        // text and not a binary blob.
        let (model, roots) = everything();
        let first = write(&model, &roots, WriteOptions::default()).unwrap();
        let (restored, restored_roots) = read(&first).unwrap();
        let second = write(&restored, &restored_roots, WriteOptions::default()).unwrap();

        if first != second {
            let mismatch = first
                .lines()
                .zip(second.lines())
                .enumerate()
                .find(|(_, (a, b))| a != b);
            panic!("the round trip changed the document: {mismatch:?}");
        }
        assert_eq!(restored_roots.len(), roots.len());
    }

    #[test]
    fn the_restored_model_is_the_same_model() {
        // Byte equality of the file is necessary and not sufficient: what
        // matters is that the shapes still answer the same way. Same counts,
        // same validity, same measurements, same handles.
        let (model, roots) = everything();
        let text = write(&model, &roots, WriteOptions::default()).unwrap();
        let (restored, restored_roots) = read(&text).unwrap();

        assert_eq!(restored.node_count(), model.node_count());
        assert_eq!(restored.geometry().counts(), model.geometry().counts());
        assert_eq!(restored.current_operation(), model.current_operation());

        for (before, after) in roots.iter().zip(&restored_roots) {
            // The *slot* came back, not merely something like it: same index,
            // same generation. Not the same handle, though — a restored
            // document is a new set of arenas and its handles say so, which is
            // what stops one document's shape from resolving against another.
            assert_eq!(before.node().index(), after.node().index());
            assert_eq!(before.node().generation(), after.node().generation());
            assert_ne!(
                before.node().scope(),
                after.node().scope(),
                "a restored document should not answer to the original's handles"
            );
            assert!(model.node(after).is_none(), "and not the other way round");
            // A reference kept across a save is re-found by its *identity*,
            // not by its handle. `bind` deliberately refuses a shape from
            // another document — relabelling one would hand back something that
            // resolves and answers about a different entity — so the way
            // through is the thing §8 exists for.
            assert!(
                restored.bind(before).is_err(),
                "a foreign handle should not be re-homed"
            );

            for kind in [
                ShapeType::Face,
                ShapeType::Edge,
                ShapeType::Vertex,
                ShapeType::Shell,
            ] {
                assert_eq!(
                    explore_unique(&restored, after, kind).unwrap().len(),
                    explore_unique(&model, before, kind).unwrap().len(),
                    "{kind:?} count changed"
                );
            }

            assert!(
                check(&restored, after, T).unwrap().is_valid(),
                "restored shape is invalid: {}",
                check(&restored, after, T).unwrap()
            );
            assert!(
                check_tessellation(&restored, after, fine(), T)
                    .unwrap()
                    .is_valid(),
                "restored shape's mesh came apart"
            );

            let before_mesh = triangulate(&model, before, fine(), T).unwrap();
            let after_mesh = triangulate(&restored, after, fine(), T).unwrap();
            assert_eq!(after_mesh.triangle_count(), before_mesh.triangle_count());
            assert_relative_eq!(
                after_mesh.volume(),
                before_mesh.volume(),
                max_relative = 0.0
            );
        }
    }

    #[test]
    fn provenance_and_identity_survive() {
        // The reason reading does not go through the builders. A rebuild
        // through them would mint fresh identities, and every reference
        // recorded against the old ones — "the top face of that box" — would
        // resolve to nothing or, worse, to something else.
        let (model, roots) = everything();
        let text = write(&model, &roots, WriteOptions::default()).unwrap();
        let (restored, restored_roots) = read(&text).unwrap();

        for (before, after) in roots.iter().zip(&restored_roots) {
            let faces = explore_unique(&model, before, ShapeType::Face).unwrap();
            let restored_faces = explore_unique(&restored, after, ShapeType::Face).unwrap();
            assert_eq!(faces.len(), restored_faces.len());
            for (a, b) in faces.iter().zip(&restored_faces) {
                assert_eq!(
                    restored.identity_of(b),
                    model.identity_of(a),
                    "an identity was renumbered"
                );
                // And the identity is enough to find the entity again in the
                // reloaded document, which is what makes a saved reference
                // survive at all.
                if let Some(id) = model.identity_of(a) {
                    let found = restored.shape_of(id).expect("the entity is still there");
                    assert!(found.is_partner(b), "identity found the wrong node");
                }
                assert_eq!(
                    restored.provenance_of(b),
                    model.provenance_of(a),
                    "a provenance record changed"
                );
                assert_eq!(restored.roots_of(b), model.roots_of(a));
            }
        }
        assert_eq!(
            restored.provenance().len(),
            model.provenance().len(),
            "the table changed length"
        );
    }

    #[test]
    fn tolerances_survive_entity_by_entity() {
        // Per-entity and they only grow (`docs/DATA_MODEL.md` §5). A format
        // that wrote one number for the document would quietly tighten every
        // edge a repair had widened.
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let edge = explore_unique(&model, &solid, ShapeType::Edge).unwrap()[0].clone();
        model.widen(&edge, Tolerance::new(1e-3).unwrap()).unwrap();

        let text = write(
            &model,
            std::slice::from_ref(&solid),
            WriteOptions::default(),
        )
        .unwrap();
        let (restored, roots) = read(&text).unwrap();
        for kind in [ShapeType::Vertex, ShapeType::Edge, ShapeType::Face] {
            let before = explore_unique(&model, &solid, kind).unwrap();
            let after = explore_unique(&restored, &roots[0], kind).unwrap();
            for (a, b) in before.iter().zip(&after) {
                assert_eq!(
                    restored.tolerance_of(b).unwrap().map(Tolerance::get),
                    model.tolerance_of(a).unwrap().map(Tolerance::get),
                    "a {kind:?} tolerance changed"
                );
            }
        }
        assert!(
            restored
                .tolerance_of(&roots[0])
                .into_iter()
                .flatten()
                .count()
                <= 1
        );
    }

    #[test]
    fn a_cached_mesh_survives_and_omitting_it_says_so() {
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T)
            .unwrap()
            .shape;
        og_mesh::tessellate(&mut model, &solid, fine(), T).unwrap();
        assert!(model.geometry().triangulation_count() > 0);

        let text = write(
            &model,
            std::slice::from_ref(&solid),
            WriteOptions::default(),
        )
        .unwrap();
        let (restored, _) = read(&text).unwrap();
        assert_eq!(
            restored.geometry().triangulation_count(),
            model.geometry().triangulation_count()
        );
        for (before, after) in model
            .geometry()
            .triangulations()
            .zip(restored.geometry().triangulations())
        {
            assert_eq!(before.1, after.1, "a cached mesh changed");
        }

        // Omitted on request, and the file says so rather than leaving it to be
        // inferred from an absence.
        let without = write(
            &model,
            std::slice::from_ref(&solid),
            WriteOptions {
                triangulations: false,
            },
        )
        .unwrap();
        assert!(without.contains("triangulation(s) omitted by request"));
        assert!(
            without.len() < text.len(),
            "omitting the mesh saved nothing"
        );
        let (bare, bare_roots) = read(&without).unwrap();
        assert_eq!(bare.geometry().triangulation_count(), 0);
        // And the shape is still a shape: the mesh was a cache, not the model.
        assert!(check(&bare, &bare_roots[0], T).unwrap().is_valid());
    }

    #[test]
    fn a_document_naming_a_handle_that_is_not_there_is_refused() {
        // A corrupt file is an error, not a model that answers about the wrong
        // entity.
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let text = write(
            &model,
            std::slice::from_ref(&solid),
            WriteOptions::default(),
        )
        .unwrap();

        let broken = text.replace("curve3d 0:0", "curve3d 999:0");
        assert_ne!(broken, text, "the substitution found nothing to do");
        assert!(read(&broken).is_err(), "a dangling curve was accepted");

        // And an identity nobody issued.
        let broken = format!("{text}identity 0:0 99999\n");
        assert!(read(&broken).is_err(), "a dangling identity was accepted");

        // And a node bound to an identity that is fine, but to a node that is
        // not there.
        let broken = format!("{text}identity 9999:0 1\n");
        assert!(
            read(&broken).is_err(),
            "a dangling node binding was accepted"
        );
    }

    #[test]
    fn records_out_of_arena_order_are_refused() {
        // The order *is* the numbering: a reader replays inserts to reproduce
        // the handles. A gap would shift every later reference by one and
        // nothing would say so.
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let text = write(
            &model,
            std::slice::from_ref(&solid),
            WriteOptions::default(),
        )
        .unwrap();
        let shuffled: String = {
            let mut lines: Vec<&str> = text.lines().collect();
            let first = lines
                .iter()
                .position(|l| l.starts_with("curve 0:0"))
                .expect("a curve");
            lines.swap(first, first + 1);
            lines.join("\n")
        };
        assert!(
            read(&shuffled).is_err(),
            "an out-of-order curve was accepted"
        );
    }

    #[test]
    fn comments_and_blank_lines_are_ignored() {
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let text = write(
            &model,
            std::slice::from_ref(&solid),
            WriteOptions::default(),
        )
        .unwrap();
        let annotated = format!("# a note\n\n{text}\n\n# and another\n");
        let (restored, roots) = read(&annotated).unwrap();
        assert_eq!(restored.node_count(), model.node_count());
        assert_eq!(roots.len(), 1);
    }

    #[test]
    fn a_root_that_is_not_in_the_model_is_refused() {
        let mut model = Model::new();
        make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let beyond = Shape::of(Key::from_parts(9999, 0));
        assert!(
            write(
                &model,
                std::slice::from_ref(&beyond),
                WriteOptions::default()
            )
            .is_err()
        );

        // A handle from a different model is caught by the same check: arena
        // keys carry the identifier of the arena that issued them, so a foreign
        // one resolves to nothing here rather than to whatever node sits at
        // that index.
        let mut elsewhere = Model::new();
        let foreign = og_algo::make_sphere(&mut elsewhere, Frame::WORLD, 1.0, T)
            .unwrap()
            .shape;
        assert!(
            write(
                &model,
                std::slice::from_ref(&foreign),
                WriteOptions::default()
            )
            .is_err()
        );
    }

    #[test]
    fn every_float_comes_back_bit_for_bit() {
        // The format is exact, not merely close. A text format that printed six
        // decimals would give a round trip that always disagrees a little, and
        // then a real disagreement has nothing to stand out against.
        let mut model = Model::new();
        let awkward = [
            0.1_f64,
            -0.0,
            1e-300,
            1.0 / 3.0,
            f64::MAX / 4.0,
            f64::MIN_POSITIVE,
            std::f64::consts::PI,
        ];
        let mut points = Vec::new();
        for value in awkward {
            points.push(model.add_point(og_math::Point::new(value, -value, value * 2.0)));
        }
        let text = write(&model, &points, WriteOptions::default()).unwrap();
        let (restored, roots) = read(&text).unwrap();
        for (before, after) in points.iter().zip(&roots) {
            let a = model
                .node(before)
                .unwrap()
                .data()
                .as_vertex()
                .unwrap()
                .point;
            let b = restored
                .node(after)
                .unwrap()
                .data()
                .as_vertex()
                .unwrap()
                .point;
            assert_eq!(a.x.to_bits(), b.x.to_bits());
            assert_eq!(a.y.to_bits(), b.y.to_bits());
            assert_eq!(a.z.to_bits(), b.z.to_bits());
        }
    }
}
