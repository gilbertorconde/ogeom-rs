//! The `.brep` interchange text format, both directions.
//!
//! An interchange format for boundary representation, implemented from its
//! published specification. It carries exact geometry and full topology,
//! which is what makes it a bridge for an application moving onto this
//! kernel that a mesh cannot be.
//!
//! The file is a header followed by tables (placements, geometry, topology),
//! and the last of those is the interesting one. The header's first line,
//! the content type, is left out by a writer saving a shape on its own (as
//! a project archive holds its parts), so the version line alone is enough.
//! Topology records are written leaves first and *numbered backwards*: the
//! final record is number one, and a record refers to its children by how
//! far above it they sit. So a parent can only name children already
//! written, the file needs no forward references, and a reader can build
//! the model in a single pass with nothing left dangling.
//!
//! ## What crosses
//!
//! Everything this kernel's own geometry can say: lines, circles, ellipses,
//! parabolas, hyperbolas, B-splines, and the trimmed and offset forms over
//! them, in two dimensions and three; planes, cylinders, cones, spheres,
//! tori, extrusions, revolutions, B-spline surfaces, and trimmed and offset
//! forms over those. All of the topology, with every occurrence's
//! orientation and placement, and edges carrying their curve, their pcurves
//! and their seams.
//!
//! ## What does not
//!
//! One thing the reader deliberately does *not* take on faith: whether an
//! edge's representations agree on parameterization. That is a claim, the
//! writing kernel's claim about its own data, and everything downstream
//! relies on it, so it is measured here instead, each representation
//! evaluated against the others at matched parameters and the claim
//! re-established only where they actually land together. A file can
//! therefore come back with the claim *established* where its writer never
//! made it, which is the right direction to be wrong in.
//!
//! Cached triangulations and the polygons that go with them are parsed and
//! skipped rather than half-honoured: this kernel tessellates on demand and
//! keeps its own cache, and a mesh read from a file would be a mesh nobody
//! could say the deflection of. The bookkeeping flags each record carries
//! are read and dropped for the same reason: they describe the state of the
//! writer's own session, not the shape.

use ogeom_core::{FastMap, FastSet};
use std::fmt::Write as _;

use ogeom_core::{OgeomResult, Tolerance, Tolerances, ogeom_bail};
use ogeom_geom::{
    BSplineCurve, BSplineSurface, Circle2d, CircleCurve, ConeSurface, Curve, Curve2d as _,
    Curve3d as _, CylinderSurface, Ellipse2d, EllipseCurve, ExtrusionSurface, HyperbolaCurve,
    Line2d, LineCurve, ParabolaCurve, PlanarCurve, PlaneSurface, RevolutionSurface, SphereSurface,
    Surface as _, SurfaceGeometry, TorusSurface, Transformable as _, TrimmedCurve, TrimmedSurface,
};
use ogeom_math::{
    Axis, Axis2, Circle, Circle2, Cone, ControlGrid, Cylinder, Direction, Direction2, Ellipse,
    Ellipse2, Frame, Frame2, Handedness, Hyperbola, Hyperbola2, KnotVector, Matrix3, Parabola,
    Plane, Point, Point2, Sphere, Torus, Transform, Transform2, TransformKind, Vector, Vector2,
    Weighted,
};
use ogeom_topo::{
    CurveId, EdgeData, EdgeRepr, FaceData, Location, Model, NodeData, Orientation, PCurveId, Shape,
    ShapeType, SurfaceId, TShape, TShapeId, VertexData,
};

/// What the format puts at the top of every file, and what a reader looks
/// for to know it is one. Data, not prose: these bytes are the format.
const CONTENT_TYPE: &str = "DBRep_DrawableShape";
const VERSION: &str = "CASCADE Topology V1, (c) Matra-Datavision";
/// What every version line opens with, whichever version it names.
const VERSION_FAMILY: &str = "CASCADE Topology V";

// --- writing -----------------------------------------------------------------

/// Write a shape as interchange text.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// shape carries geometry the format has no record for (a helical curve, a
/// curve living on a surface, the sinusoidal pcurve an oblique section
/// leaves on a cylinder), each refused by name rather than approximated into
/// something the file would claim was exact.
pub fn write(model: &Model, root: &Shape, tol: Tolerances) -> OgeomResult<String> {
    let mut tables = Tables::default();
    tables.gather(model, root)?;

    let mut out = String::new();
    let _ = writeln!(out, "{CONTENT_TYPE}\n");
    let _ = writeln!(out, "{VERSION}");

    // Placements. Every location is written as its own matrix rather than as
    // a composition of others: the composition form exists to let a file
    // share a chain between many shapes, and nothing here is repeated often
    // enough for that to pay.
    let _ = writeln!(out, "Locations {}", tables.locations.len());
    for location in &tables.locations {
        let transform = location.composed(model.datums())?;
        let linear = transform.linear();
        let scale = transform.scale_factor();
        let translation = transform.translation_vector();
        // The record's kind on its own line, then the three rows.
        let _ = writeln!(out, "1");
        for row in 0..3 {
            let mut line = String::new();
            for column in 0..3 {
                let _ = write!(line, " {}", real(linear.get(row, column)? * scale));
            }
            let _ = write!(
                line,
                " {}",
                real(match row {
                    0 => translation.x,
                    1 => translation.y,
                    _ => translation.z,
                })
            );
            let _ = writeln!(out, "{}", line.trim_start());
        }
    }

    let _ = writeln!(out, "Curve2ds {}", tables.pcurves.len());
    for id in &tables.pcurves {
        let Some(geometry) = model.geometry().pcurve(*id) else {
            ogeom_bail!(Dangling, "a pcurve is not in this model");
        };
        // A cone's chart runs along its axis here and along its ruling in
        // the format.
        match tables.ruled.get(id) {
            Some(cos) => write_pcurve(&mut out, &geometry.stretched(1.0, 1.0 / cos, tol)?)?,
            None => write_pcurve(&mut out, geometry)?,
        }
    }

    let _ = writeln!(out, "Curves {}", tables.curves.len());
    for id in &tables.curves {
        let Some(geometry) = model.geometry().curve(*id) else {
            ogeom_bail!(Dangling, "a curve is not in this model");
        };
        write_curve(&mut out, geometry)?;
    }

    // The two cached-mesh sections, empty. They are part of the file's
    // shape whether or not anything is in them.
    let _ = writeln!(out, "Polygon3D 0");
    let _ = writeln!(out, "PolygonOnTriangulations 0");

    let _ = writeln!(out, "Surfaces {}", tables.surfaces.len());
    for id in &tables.surfaces {
        let Some(geometry) = model.geometry().surface(*id) else {
            ogeom_bail!(Dangling, "a surface is not in this model");
        };
        write_surface(&mut out, geometry)?;
    }
    let _ = writeln!(out, "Triangulations 0");

    // Topology, leaves first. Record numbering counts backwards from the
    // last, so a child written earlier has the larger number.
    let total = tables.order.len();
    let _ = writeln!(out, "\nTShapes {total}");
    let number_of = |node: TShapeId| -> OgeomResult<usize> {
        let Some(position) = tables.index_of.get(&node) else {
            ogeom_bail!(Construction, "a subshape was never written");
        };
        Ok(total - position)
    };

    for node_id in &tables.order {
        let shape = Shape::of(*node_id);
        let Some(node) = model.node(&shape) else {
            ogeom_bail!(Dangling, "a shape is not in this model");
        };
        write_record(&mut out, model, &shape, &tables, tol)?;

        // Flags, then the children. Seven flags the format asks for and this
        // kernel does not keep: they record the writing session's own
        // bookkeeping (whether a shape was visited, modified, checked), and
        // reading them back as geometry would be reading somebody else's
        // scratch paper. Orientable is the one that is always true here.
        let _ = writeln!(out, "\n0001000");
        let mut line = String::new();
        for child in node.children() {
            let _ = write!(
                line,
                "{}{} {} ",
                match child.orientation() {
                    Orientation::Forward => "+",
                    Orientation::Reversed => "-",
                    Orientation::Internal => "i",
                    Orientation::External => "e",
                },
                number_of(child.node())?,
                tables.location_number(child.location()),
            );
        }
        let _ = writeln!(out, "{line}*");
    }

    // The final record says how the whole model is oriented and placed.
    let _ = writeln!(
        out,
        "\n{}{} {}",
        match root.orientation() {
            Orientation::Forward => "+",
            Orientation::Reversed => "-",
            Orientation::Internal => "i",
            Orientation::External => "e",
        },
        number_of(root.node())?,
        tables.location_number(root.location()),
    );
    Ok(out)
}

/// The tables a file's sections are: geometry by the order it is written,
/// topology leaves first, and every placement anything refers to.
#[derive(Default)]
struct Tables {
    curves: Vec<CurveId>,
    pcurves: Vec<PCurveId>,
    surfaces: Vec<SurfaceId>,
    locations: Vec<Location>,
    order: Vec<TShapeId>,
    index_of: FastMap<TShapeId, usize>,
    curve_at: FastMap<CurveId, usize>,
    pcurve_at: FastMap<PCurveId, usize>,
    surface_at: FastMap<SurfaceId, usize>,
    location_at: FastMap<Location, usize>,
    /// Charts on a cone, with the cosine of its half-angle.
    ruled: FastMap<PCurveId, f64>,
    /// Where each of an edge's charts is written as used: the placement of
    /// the face's surface as seen from the edge, once per use the faces
    /// make of it, by the chart's index among the edge's representations.
    chart_uses: FastMap<TShapeId, Vec<(usize, Location)>>,
}

impl Tables {
    /// Walk the shape leaves-first, numbering everything it refers to.
    fn gather(&mut self, model: &Model, root: &Shape) -> OgeomResult<()> {
        self.visit(model, root)?;
        // The file's own shape is placed by a record like any other.
        self.take_location(root.location());
        Ok(())
    }

    fn visit(&mut self, model: &Model, shape: &Shape) -> OgeomResult<()> {
        if self.index_of.contains_key(&shape.node()) {
            return Ok(());
        }
        let Some(node) = model.node(shape) else {
            ogeom_bail!(Dangling, "shape is not in this model");
        };
        // Children first, so a parent's record can name them by how far
        // above it they were written.
        for child in node.children() {
            self.visit(model, child)?;
            self.take_location(child.location());
        }
        match node.data() {
            NodeData::Edge(data) => {
                for representation in &data.representations {
                    match representation {
                        EdgeRepr::Curve3d {
                            curve, location, ..
                        } => {
                            self.take_curve(*curve);
                            self.take_location(location);
                        }
                        EdgeRepr::PCurve {
                            curve,
                            surface,
                            location,
                            ..
                        } => {
                            self.take_pcurve(*curve);
                            self.take_surface(*surface);
                            self.take_location(&location.inverted());
                            self.note_ruled(model, *surface, &[*curve]);
                        }
                        EdgeRepr::Seam {
                            forward,
                            reversed,
                            surface,
                            location,
                            ..
                        } => {
                            self.take_pcurve(*forward);
                            self.take_pcurve(*reversed);
                            self.take_surface(*surface);
                            self.take_location(&location.inverted());
                            self.note_ruled(model, *surface, &[*forward, *reversed]);
                        }
                        _ => {}
                    }
                }
            }
            NodeData::Face(data) => {
                self.take_surface(data.surface);
                self.take_location(&data.location);
                // The format names an edge's chart on a face by where the
                // face's surface sits as seen from the edge: the face's own
                // placement taken back through the edge's inside the face.
                for wire in model.children_of(&Shape::of(shape.node()))? {
                    for edge in model.children_of(&wire)? {
                        let Some(edge_data) = model.node(&edge).and_then(|n| n.data().as_edge())
                        else {
                            continue;
                        };
                        let Some(index) = edge_data.representations.iter().position(|r| {
                            Some(r) == edge_data.pcurve_for(data.surface, edge.location())
                        }) else {
                            continue;
                        };
                        let seen = edge.location().inverted().then(&data.location);
                        self.take_location(&seen);
                        let uses = self.chart_uses.entry(edge.node()).or_default();
                        if !uses.contains(&(index, seen.clone())) {
                            uses.push((index, seen));
                        }
                    }
                }
            }
            NodeData::Vertex(_) | NodeData::Container => {}
        }
        self.index_of.insert(shape.node(), self.order.len());
        self.order.push(shape.node());
        Ok(())
    }

    fn note_ruled(&mut self, model: &Model, surface: SurfaceId, curves: &[PCurveId]) {
        if let Some(SurfaceGeometry::Cone(cone)) = model.geometry().surface(surface) {
            for curve in curves {
                self.ruled.insert(*curve, cone.cone().half_angle().cos());
            }
        }
    }

    fn take_curve(&mut self, id: CurveId) {
        if !self.curve_at.contains_key(&id) {
            self.curve_at.insert(id, self.curves.len());
            self.curves.push(id);
        }
    }

    fn take_pcurve(&mut self, id: PCurveId) {
        if !self.pcurve_at.contains_key(&id) {
            self.pcurve_at.insert(id, self.pcurves.len());
            self.pcurves.push(id);
        }
    }

    fn take_surface(&mut self, id: SurfaceId) {
        if !self.surface_at.contains_key(&id) {
            self.surface_at.insert(id, self.surfaces.len());
            self.surfaces.push(id);
        }
    }

    fn take_location(&mut self, location: &Location) {
        if location.is_identity() || self.location_at.contains_key(location) {
            return;
        }
        self.location_at
            .insert(location.clone(), self.locations.len());
        self.locations.push(location.clone());
    }

    /// The file's number for a placement: one-based, and zero for identity,
    /// which the format spells as "no placement at all".
    fn location_number(&self, location: &Location) -> usize {
        if location.is_identity() {
            return 0;
        }
        self.location_at.get(location).map_or(0, |at| at + 1)
    }

    fn curve_number(&self, id: CurveId) -> usize {
        self.curve_at.get(&id).map_or(0, |at| at + 1)
    }

    fn pcurve_number(&self, id: PCurveId) -> usize {
        self.pcurve_at.get(&id).map_or(0, |at| at + 1)
    }

    fn surface_number(&self, id: SurfaceId) -> usize {
        self.surface_at.get(&id).map_or(0, |at| at + 1)
    }
}

/// One topology record's own data: everything before the flag word.
fn write_record(
    out: &mut String,
    model: &Model,
    shape: &Shape,
    tables: &Tables,
    tol: Tolerances,
) -> OgeomResult<()> {
    let Some(node) = model.node(shape) else {
        ogeom_bail!(Dangling, "a shape is not in this model");
    };
    match node.data() {
        NodeData::Vertex(data) => {
            let _ = writeln!(out, "Ve");
            let _ = writeln!(out, "{}", real(data.tolerance.get()));
            let _ = writeln!(
                out,
                "{} {} {}",
                real(data.point.x),
                real(data.point.y),
                real(data.point.z)
            );
            // A vertex's parameters on the curves through it are recoverable
            // by projection and are not written, so the list is empty.
            let _ = writeln!(out, "0 0");
        }
        NodeData::Edge(data) => {
            let _ = writeln!(out, "Ed");
            let _ = writeln!(
                out,
                " {} {} {} {}",
                real(data.tolerance.get()),
                usize::from(data.same_parameter()),
                1,
                usize::from(data.degenerate)
            );
            for (index, representation) in data.representations.iter().enumerate() {
                // Each use the faces make of a chart is its own line, at
                // the placement that use names; a chart no face uses is
                // written as if its face sat where the edge does.
                let uses: Vec<&Location> = tables
                    .chart_uses
                    .get(&shape.node())
                    .map(|uses| {
                        uses.iter()
                            .filter(|(at, _)| *at == index)
                            .map(|(_, seen)| seen)
                            .collect()
                    })
                    .unwrap_or_default();
                let inverted = representation.location().map(Location::inverted);
                let seen: Vec<&Location> = if uses.is_empty() {
                    inverted.iter().collect()
                } else {
                    uses
                };
                match representation {
                    EdgeRepr::Curve3d {
                        curve,
                        range,
                        location,
                    } => {
                        let _ = writeln!(
                            out,
                            "1 {} {} {} {}",
                            tables.curve_number(*curve),
                            tables.location_number(location),
                            real(range.0),
                            real(range.1)
                        );
                    }
                    EdgeRepr::PCurve {
                        curve,
                        range,
                        surface,
                        ..
                    } => {
                        for seen in seen {
                            let _ = writeln!(
                                out,
                                "2 {} {} {} {} {}",
                                tables.pcurve_number(*curve),
                                tables.surface_number(*surface),
                                tables.location_number(seen),
                                real(range.0),
                                real(range.1)
                            );
                        }
                    }
                    EdgeRepr::Seam {
                        forward,
                        reversed,
                        range,
                        surface,
                        ..
                    } => {
                        for seen in seen {
                            let _ = writeln!(
                                out,
                                "3 {} {} C0 {} {} {} {}",
                                tables.pcurve_number(*forward),
                                tables.pcurve_number(*reversed),
                                tables.surface_number(*surface),
                                tables.location_number(seen),
                                real(range.0),
                                real(range.1)
                            );
                        }
                    }
                    other => ogeom_bail!(
                        Construction,
                        "an edge carries a representation this format has no \
                         record for: {other:?}"
                    ),
                }
            }
            let _ = writeln!(out, "0");
        }
        NodeData::Face(data) => {
            let _ = writeln!(out, "Fa");
            let _ = writeln!(
                out,
                "0 {} {} {}",
                real(data.tolerance.get()),
                tables.surface_number(data.surface),
                tables.location_number(&data.location)
            );
        }
        NodeData::Container => {
            let kind = match model.kind_of(shape)? {
                ShapeType::Wire => "Wi",
                ShapeType::Shell => "Sh",
                ShapeType::Solid => "So",
                ShapeType::CompSolid => "CS",
                ShapeType::Compound => "Co",
                other => ogeom_bail!(Construction, "a {other:?} has no record in this format"),
            };
            let _ = writeln!(out, "{kind}\n");
        }
    }
    let _ = tol;
    Ok(())
}

fn write_curve(out: &mut String, curve: &Curve) -> OgeomResult<()> {
    match curve {
        Curve::Line(l) => {
            let axis = l.axis();
            let _ = writeln!(
                out,
                "1 {} {}",
                point(axis.location),
                direction(axis.direction.vector())
            );
        }
        Curve::Circle(c) => {
            let circle = c.circle();
            let frame = circle.frame();
            let _ = writeln!(
                out,
                "2 {} {} {} {} {}",
                point(frame.origin()),
                direction(frame.z().vector()),
                direction(frame.x().vector()),
                direction(frame.y().vector()),
                real(circle.radius())
            );
        }
        Curve::Ellipse(e) => {
            let ellipse = e.ellipse();
            let frame = ellipse.frame();
            let _ = writeln!(
                out,
                "3 {} {} {} {} {} {}",
                point(frame.origin()),
                direction(frame.z().vector()),
                direction(frame.x().vector()),
                direction(frame.y().vector()),
                real(ellipse.major_radius()),
                real(ellipse.minor_radius())
            );
        }
        Curve::Parabola(p) => {
            let parabola = p.parabola();
            let frame = parabola.frame();
            let _ = writeln!(
                out,
                "4 {} {} {} {} {}",
                point(frame.origin()),
                direction(frame.z().vector()),
                direction(frame.x().vector()),
                direction(frame.y().vector()),
                real(parabola.focal())
            );
        }
        Curve::Hyperbola(h) => {
            let hyperbola = h.hyperbola();
            let frame = hyperbola.frame();
            let _ = writeln!(
                out,
                "5 {} {} {} {} {} {}",
                point(frame.origin()),
                direction(frame.z().vector()),
                direction(frame.x().vector()),
                direction(frame.y().vector()),
                real(hyperbola.major_radius()),
                real(hyperbola.minor_radius())
            );
        }
        Curve::BSpline(b) => {
            let rational = b
                .control_points()
                .iter()
                .any(|c| (c.weight - 1.0).abs() > 0.0);
            let distinct = b.knots().distinct();
            let mut line = format!(
                "7 {} 0  {} {} {}",
                usize::from(rational),
                b.knots().degree(),
                b.control_points().len(),
                distinct.len()
            );
            for control in b.control_points() {
                let at = control.point();
                let _ = write!(line, " {}", point(at));
                if rational {
                    let _ = write!(line, " {}", real(control.weight));
                }
            }
            let _ = writeln!(out, "{line}");
            let mut knot_line = String::new();
            for (value, multiplicity) in distinct {
                let _ = write!(knot_line, " {} {multiplicity}", real(value));
            }
            let _ = writeln!(out, "{knot_line}");
        }
        Curve::Trimmed(t) => {
            let (lo, hi) = t.domain();
            let _ = writeln!(out, "8 {} {}", real(lo), real(hi));
            write_curve(out, t.basis())?;
        }
        Curve::Offset(o) => {
            let _ = writeln!(
                out,
                "9 {} {}",
                real(o.distance()),
                direction(o.reference().vector())
            );
            write_curve(out, o.basis())?;
        }
        Curve::Helix(_) | Curve::OnSurface(_) => ogeom_bail!(
            Construction,
            "the format has no record for a helix or a curve carried on a \
             surface, and writing one as anything else would be writing a \
             different curve"
        ),
    }
    Ok(())
}

fn write_pcurve(out: &mut String, curve: &PlanarCurve) -> OgeomResult<()> {
    match curve {
        PlanarCurve::Line(l) => {
            let axis = l.axis();
            let _ = writeln!(
                out,
                "1 {} {}",
                point2(axis.location),
                direction2(axis.direction.vector())
            );
        }
        PlanarCurve::Circle(c) => {
            let circle = c.circle();
            let frame = circle.frame();
            let _ = writeln!(
                out,
                "2 {} {} {} {}",
                point2(frame.origin()),
                direction2(frame.x().vector()),
                direction2(frame.y().vector()),
                real(circle.radius())
            );
        }
        PlanarCurve::Ellipse(e) => {
            let ellipse = e.ellipse();
            let frame = ellipse.frame();
            let _ = writeln!(
                out,
                "3 {} {} {} {} {}",
                point2(frame.origin()),
                direction2(frame.x().vector()),
                direction2(frame.y().vector()),
                real(ellipse.major_radius()),
                real(ellipse.minor_radius())
            );
        }
        PlanarCurve::Hyperbola(h) => {
            // The record states the branch over its natural parameter,
            // forward; a branch run backwards has no record.
            if h.is_reversed() {
                ogeom_bail!(
                    Construction,
                    "the format has no record for a hyperbola run backwards"
                );
            }
            let hyperbola = h.hyperbola();
            let frame = hyperbola.frame();
            let _ = writeln!(
                out,
                "5 {} {} {} {} {}",
                point2(frame.origin()),
                direction2(frame.x().vector()),
                direction2(frame.y().vector()),
                real(hyperbola.major_radius()),
                real(hyperbola.minor_radius())
            );
        }
        PlanarCurve::BSpline(b) => {
            let rational = b
                .control_points()
                .iter()
                .any(|c| (c.weight - 1.0).abs() > 0.0);
            let distinct = b.knots().distinct();
            let mut line = format!(
                "7 {} 0  {} {} {}",
                usize::from(rational),
                b.knots().degree(),
                b.control_points().len(),
                distinct.len()
            );
            for control in b.control_points() {
                let at = control.point();
                let _ = write!(line, " {}", point2(at));
                if rational {
                    let _ = write!(line, " {}", real(control.weight));
                }
            }
            let _ = writeln!(out, "{line}");
            let mut knot_line = String::new();
            for (value, multiplicity) in distinct {
                let _ = write!(knot_line, " {} {multiplicity}", real(value));
            }
            let _ = writeln!(out, "{knot_line}");
        }
        PlanarCurve::Trimmed(t) => {
            let (lo, hi) = t.domain();
            let _ = writeln!(out, "8 {} {}", real(lo), real(hi));
            write_pcurve(out, t.basis())?;
        }
        PlanarCurve::Offset(o) => {
            let _ = writeln!(out, "9 {}", real(o.distance()));
            write_pcurve(out, o.basis())?;
        }
        PlanarCurve::Trig(_) => ogeom_bail!(
            Construction,
            "the format has no record for the sinusoidal chart curve an \
             oblique section leaves on a cylinder; writing it as a spline \
             would be writing a fit and calling it exact"
        ),
    }
    Ok(())
}

fn write_surface(out: &mut String, surface: &SurfaceGeometry) -> OgeomResult<()> {
    // The elementary surfaces all begin the same way: an origin and the
    // frame's three directions, normal first.
    let seat = |frame: &Frame| {
        format!(
            "{} {} {} {}",
            point(frame.origin()),
            direction(frame.z().vector()),
            direction(frame.x().vector()),
            direction(frame.y().vector())
        )
    };
    match surface {
        SurfaceGeometry::Plane(p) => {
            let _ = writeln!(out, "1 {}", seat(&p.plane().frame()));
        }
        SurfaceGeometry::Cylinder(c) => {
            let _ = writeln!(
                out,
                "2 {} {}",
                seat(&c.cylinder().frame()),
                real(c.cylinder().radius())
            );
        }
        SurfaceGeometry::Cone(c) => {
            let cone = c.cone();
            let _ = writeln!(
                out,
                "3 {} {} {}",
                seat(&cone.frame()),
                real(cone.reference_radius()),
                real(cone.half_angle())
            );
        }
        SurfaceGeometry::Sphere(s) => {
            let _ = writeln!(
                out,
                "4 {} {}",
                seat(&s.sphere().frame()),
                real(s.sphere().radius())
            );
        }
        SurfaceGeometry::Torus(t) => {
            let torus = t.torus();
            let _ = writeln!(
                out,
                "5 {} {} {}",
                seat(&torus.frame()),
                real(torus.major_radius()),
                real(torus.minor_radius())
            );
        }
        SurfaceGeometry::Extrusion(e) => {
            let _ = writeln!(out, "6 {}", direction(e.direction().vector()));
            write_curve(out, e.curve())?;
        }
        SurfaceGeometry::Revolution(r) => {
            let axis = r.axis();
            let _ = writeln!(
                out,
                "7 {} {}",
                point(axis.location),
                direction(axis.direction.vector())
            );
            write_curve(out, r.curve())?;
        }
        SurfaceGeometry::BSpline(b) => {
            let grid = b.grid();
            let rational = grid.points().iter().any(|c| (c.weight - 1.0).abs() > 0.0);
            let u_distinct = b.u_knots().distinct();
            let v_distinct = b.v_knots().distinct();
            let _ = writeln!(
                out,
                "9 {} {} 0 0 {} {} {} {} {} {}",
                usize::from(rational),
                usize::from(rational),
                b.u_knots().degree(),
                b.v_knots().degree(),
                grid.u_count(),
                grid.v_count(),
                u_distinct.len(),
                v_distinct.len()
            );
            // Poles run u-major: all of one u row's v values, then the next.
            for i in 0..grid.u_count() {
                let mut line = String::new();
                for j in 0..grid.v_count() {
                    let control = grid.points()[i * grid.v_count() + j];
                    let at = control.point();
                    let _ = write!(line, " {}", point(at));
                    if rational {
                        let _ = write!(line, " {}", real(control.weight));
                    }
                }
                let _ = writeln!(out, "{line}");
            }
            let mut line = String::new();
            for (value, multiplicity) in u_distinct {
                let _ = write!(line, " {} {multiplicity}", real(value));
            }
            let _ = writeln!(out, "{line}");
            let mut line = String::new();
            for (value, multiplicity) in v_distinct {
                let _ = write!(line, " {} {multiplicity}", real(value));
            }
            let _ = writeln!(out, "{line}");
        }
        SurfaceGeometry::Trimmed(t) => {
            let ((u0, u1), (v0, v1)) = t.domain();
            let _ = writeln!(
                out,
                "10 {} {} {} {}",
                real(u0),
                real(u1),
                real(v0),
                real(v1)
            );
            write_surface(out, t.basis())?;
        }
        SurfaceGeometry::Offset(o) => {
            let _ = writeln!(out, "11 {}", real(o.distance()));
            write_surface(out, o.basis())?;
        }
    }
    Ok(())
}

/// A number, written so it reads back as itself.
fn real(v: f64) -> String {
    let mut text = format!("{v:?}");
    if text.ends_with(".0") {
        text.truncate(text.len() - 2);
    }
    text
}

fn point(p: Point) -> String {
    format!("{} {} {}", real(p.x), real(p.y), real(p.z))
}

fn point2(p: Point2) -> String {
    format!("{} {}", real(p.x), real(p.y))
}

fn direction(v: Vector) -> String {
    format!("{} {} {}", real(v.x), real(v.y), real(v.z))
}

fn direction2(v: Vector2) -> String {
    format!("{} {}", real(v.x), real(v.y))
}

// --- reading -----------------------------------------------------------------

/// Read a shape from interchange text, into a model of its own.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// text is not this format, or carries a record this does not read, which
/// it says by number rather than by shrugging.
pub fn read(text: &str, tol: Tolerances) -> OgeomResult<(Model, Shape)> {
    let (model, shape, _) = read_numbered(text, tol)?;
    Ok((model, shape))
}

/// A shape's edges and faces in the order the file numbers them.
///
/// An application that names subshapes by number (the fifth edge, the
/// second face) counts them the way the file is written: depth first from
/// the file's own shape, each record's references in the order written,
/// every edge or face counted the first time it is met at a placement.
/// Faces come out in that order from exploring the shape read too, but
/// edges need not: the file keeps a wire's edges as a set, in whatever
/// order its writer met them, and the wire read is walked end to end.
#[derive(Debug, Clone, Default)]
pub struct Numbering {
    /// The edges, first numbered first.
    pub edges: Vec<Shape>,
    /// The faces, first numbered first.
    pub faces: Vec<Shape>,
}

/// [`read`], with the edges and faces as the file numbers them.
///
/// # Errors
///
/// As [`read`].
pub fn read_numbered(text: &str, tol: Tolerances) -> OgeomResult<(Model, Shape, Numbering)> {
    // A file written for the drawing tool opens with the content-type line;
    // one written straight from a shape (what project archives hold) opens
    // with the version line, often after a blank one.
    let mut lines = text.lines().skip_while(|line| line.trim().is_empty());
    let Some(first) = lines.next() else {
        ogeom_bail!(Construction, "the file is empty");
    };
    let first = first.trim();
    if first != CONTENT_TYPE && !first.starts_with(VERSION_FAMILY) {
        ogeom_bail!(
            Construction,
            "this is not the interchange format: it opens with `{first}`"
        );
    }
    // The version line is whatever the writer put there; every version this
    // reads differs only in sections this skips.
    let rest: String = lines.collect::<Vec<_>>().join("\n");
    let Some(at) = rest.find("Locations") else {
        ogeom_bail!(Construction, "the file has no placement table");
    };
    let mut cursor = Cursor::new(&rest[at..]);
    let mut model = Model::new();
    let mut built = Built::default();

    cursor.expect("Locations")?;
    let count = cursor.count()?;
    for _ in 0..count {
        let kind = cursor.count()?;
        match kind {
            1 => {
                let mut rows = [[0.0_f64; 3]; 3];
                let mut translation = [0.0_f64; 3];
                for (row, offset) in rows.iter_mut().zip(&mut translation) {
                    for cell in row.iter_mut() {
                        *cell = cursor.number()?;
                    }
                    *offset = cursor.number()?;
                }
                let transform = transform_of(rows, translation, tol)?;
                // One datum per placement that moves anything, so two
                // shapes placed by it are at the same place structurally;
                // one that moves nothing is left out of every chain.
                built
                    .chains
                    .push(if transform.kind() == TransformKind::Identity {
                        Location::identity()
                    } else {
                        Location::of(model.add_datum(transform))
                    });
                built.locations.push(transform);
            }
            2 => {
                // A composition of placements already read, each raised to a
                // power, ending at a zero. Each factor applies after the ones
                // listed before it: the first is the innermost. The chain
                // keeps the factors, as the file does, so placements built
                // from the same factors are the same placement.
                let mut composed = Transform::IDENTITY;
                let mut chain = Location::identity();
                loop {
                    let index = cursor.count()?;
                    if index == 0 {
                        break;
                    }
                    let power: i64 = cursor.integer()?;
                    let (Some(base), Some(factor)) = (
                        built.locations.get(index - 1).copied(),
                        built.chains.get(index - 1).cloned(),
                    ) else {
                        ogeom_bail!(
                            Construction,
                            "a composed placement names placement {index}, which is not above it"
                        );
                    };
                    let (step, factor) = if power < 0 {
                        (base.inverse()?, factor.inverted())
                    } else {
                        (base, factor)
                    };
                    for _ in 0..power.abs() {
                        composed = step * composed;
                        chain = factor.then(&chain);
                    }
                }
                built.locations.push(composed);
                built.chains.push(chain);
            }
            other => ogeom_bail!(
                Construction,
                "placement record {other} is not one this reads"
            ),
        }
    }

    cursor.expect("Curve2ds")?;
    let count = cursor.count()?;
    for _ in 0..count {
        let curve = cursor.pcurve(tol)?;
        built.pcurves.push(model.geometry_mut().add_pcurve(curve));
    }

    cursor.expect("Curves")?;
    let count = cursor.count()?;
    for _ in 0..count {
        let curve = cursor.curve(tol)?;
        built.curves.push(model.geometry_mut().add_curve(curve));
    }

    // The cached-mesh sections, counted so they can be stepped over.
    cursor.expect("Polygon3D")?;
    let count = cursor.count()?;
    for _ in 0..count {
        cursor.skip_polygon3d()?;
    }
    cursor.expect("PolygonOnTriangulations")?;
    let count = cursor.count()?;
    for _ in 0..count {
        cursor.skip_polygon_on_triangulation()?;
    }

    cursor.expect("Surfaces")?;
    let count = cursor.count()?;
    for _ in 0..count {
        let surface = cursor.surface(tol)?;
        built
            .surfaces
            .push(model.geometry_mut().add_surface(surface));
    }

    cursor.expect("Triangulations")?;
    let count = cursor.count()?;
    for _ in 0..count {
        cursor.skip_triangulation()?;
    }

    cursor.expect("TShapes")?;
    let total = cursor.count()?;
    // Records are read in file order and numbered backwards, so the record
    // just read is number `total - written`, and every child it names has
    // already been built.
    let mut shapes: Vec<Shape> = Vec::with_capacity(cursor.room(total));
    for _ in 0..total {
        let shape = cursor.record(&mut model, &mut built, &shapes, total, tol)?;
        shapes.push(shape);
    }
    // The faces have given every edge its charts: whether those agree with
    // the edge's curve can be measured here, and not before.
    for edge in &shapes {
        let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let agrees = agree_on_parameter(&model, data, tol);
        if let Some(NodeData::Edge(data)) = model.node_mut(edge).map(TShape::data_mut) {
            data.assert_same_parameter(agrees);
        }
    }
    // The last reference is the file's own shape: how the whole model is
    // oriented and placed.
    let (orientation, number, location) = cursor.reference()?;
    let Some(root) = pick(&shapes, number, total) else {
        ogeom_bail!(
            Construction,
            "the file's own shape is numbered {number}, which is not there"
        );
    };
    let sense = built.sense(&root, orientation);
    let root = placed(root, &built, location)?.composed(sense);
    let numbering = numbered(&model, &built, &shapes, (number, location))?;
    Ok((model, root, numbering))
}

/// How far an unbounded conic runs when it is read.
///
/// The format's parabola and hyperbola records carry a focus and a frame and
/// no parameter window; the window is the *edge's* business, and every edge
/// states its own range. So the curve is built over a window wide enough for
/// any edge a file of ordinary size could hold, and the edge trims it.
const OPEN_EXTENT: f64 = 1e6;

/// Everything read so far that later records refer to by number.
#[derive(Default)]
struct Built {
    locations: Vec<Transform>,
    /// Each placement as a chain of the file's elementary placements.
    chains: Vec<Location>,
    curves: Vec<CurveId>,
    pcurves: Vec<PCurveId>,
    surfaces: Vec<SurfaceId>,
    /// Curves and surfaces carried to where a record places them, by the
    /// file's curve or surface number and placement number.
    placed_curves: FastMap<(usize, usize), CurveId>,
    placed_surfaces: FastMap<(usize, usize), SurfaceId>,
    /// Each edge's charts as the file states them, held until the faces
    /// that use them are read: a chart is bound to a face occurrence, which
    /// the edge record alone cannot name.
    charts: FastMap<TShapeId, Vec<StatedChart>>,
    /// How the file's charts read on a surface whose chart is not the
    /// file's own, and the charts carried onto them.
    chart_maps: FastMap<SurfaceId, ChartMap>,
    mapped_charts: FastMap<(PCurveId, SurfaceId), PCurveId>,
    /// Faces on such a surface: each built with its wires walked the other
    /// way and used turned over, so it bounds what the file says it does.
    turned_faces: FastSet<TShapeId>,
    /// Each record's references as the file writes them, in its order: a
    /// wire here is a walk, and the file's own order is what numbers its
    /// edges.
    written: Vec<Vec<(usize, usize)>>,
    /// Faces rebuilt where a shell places them, by the face and that
    /// placement.
    settled: FastMap<(TShapeId, Location), TShapeId>,
}

/// How a chart curve in the file reads on a surface as it is built here.
#[derive(Clone, Copy)]
struct ChartMap {
    /// The surface was rebuilt right-handed and its chart turned over.
    turned: bool,
    /// What `v` is multiplied by: a cone's chart runs along the ruling in
    /// the file and along the axis here.
    v_scale: f64,
}

/// One chart representation of an edge, in the file's own numbers.
struct StatedChart {
    /// One curve, or a seam's two (forward use first).
    curves: (PCurveId, Option<PCurveId>),
    surface: usize,
    location: usize,
    range: (f64, f64),
}

impl Built {
    /// How a reference uses a shape: as the file says, turned over for a
    /// face built turned over.
    fn sense(&self, shape: &Shape, stated: Orientation) -> Orientation {
        if self.turned_faces.contains(&shape.node()) {
            stated.compose(Orientation::Reversed)
        } else {
            stated
        }
    }

    /// The transform the file's placement number names; zero is none.
    fn transform(&self, number: usize) -> OgeomResult<Transform> {
        if number == 0 {
            return Ok(Transform::IDENTITY);
        }
        let Some(transform) = self.locations.get(number - 1) else {
            ogeom_bail!(
                Construction,
                "a record names placement {number}, which is not there"
            );
        };
        Ok(*transform)
    }

    /// The placement the file's number names, as a location.
    fn location(&self, number: usize) -> OgeomResult<Location> {
        if number == 0 {
            return Ok(Location::identity());
        }
        let Some(chain) = self.chains.get(number - 1) else {
            ogeom_bail!(
                Construction,
                "a record names placement {number}, which is not there"
            );
        };
        Ok(chain.clone())
    }

    /// Curve `number` where placement `location` puts it.
    fn curve(
        &mut self,
        model: &mut Model,
        number: usize,
        location: usize,
        tol: Tolerances,
    ) -> OgeomResult<CurveId> {
        let Some(id) = number
            .checked_sub(1)
            .and_then(|at| self.curves.get(at))
            .copied()
        else {
            ogeom_bail!(Construction, "an edge names curve {number}");
        };
        let transform = self.transform(location)?;
        if transform.kind() == TransformKind::Identity {
            return Ok(id);
        }
        if let Some(placed) = self.placed_curves.get(&(number, location)) {
            return Ok(*placed);
        }
        let Some(curve) = model.geometry().curve(id) else {
            ogeom_bail!(Dangling, "a curve read from the file is not in the model");
        };
        let moved = curve.transformed(&transform, tol)?;
        let placed = model.geometry_mut().add_curve(moved);
        self.placed_curves.insert((number, location), placed);
        Ok(placed)
    }

    /// Surface `number` where placement `location` puts it, on a
    /// right-handed seat.
    fn surface(
        &mut self,
        model: &mut Model,
        number: usize,
        location: usize,
        tol: Tolerances,
    ) -> OgeomResult<SurfaceId> {
        let Some(id) = number
            .checked_sub(1)
            .and_then(|at| self.surfaces.get(at))
            .copied()
        else {
            ogeom_bail!(Construction, "a face names surface {number}");
        };
        if let Some(placed) = self.placed_surfaces.get(&(number, location)) {
            return Ok(*placed);
        }
        let Some(surface) = model.geometry().surface(id) else {
            ogeom_bail!(Dangling, "a surface read from the file is not in the model");
        };
        let transform = self.transform(location)?;
        let moved = transform.kind() != TransformKind::Identity;
        let surface = if moved {
            surface.transformed(&transform, tol)?
        } else {
            surface.clone()
        };
        let (surface, rebuilt, turned) = match right_handed(&surface, tol)? {
            Some((upright, turned)) => (upright, true, turned),
            None => (surface, false, false),
        };
        let v_scale = match &surface {
            SurfaceGeometry::Cone(c) => c.cone().half_angle().cos(),
            _ => 1.0,
        };
        let placed = if moved || rebuilt {
            model.geometry_mut().add_surface(surface)
        } else {
            id
        };
        if turned || v_scale != 1.0 {
            self.chart_maps.insert(placed, ChartMap { turned, v_scale });
        }
        self.placed_surfaces.insert((number, location), placed);
        Ok(placed)
    }

    /// Whether faces on surface `on` are built turned over.
    fn turned(&self, on: SurfaceId) -> bool {
        self.chart_maps.get(&on).is_some_and(|map| map.turned)
    }

    /// A chart curve as it reads on surface `on`, carried by the surface's
    /// chart map where it has one.
    fn chart_on(
        &mut self,
        model: &mut Model,
        curve: PCurveId,
        on: SurfaceId,
        tol: Tolerances,
    ) -> OgeomResult<PCurveId> {
        let Some(map) = self.chart_maps.get(&on).copied() else {
            return Ok(curve);
        };
        if let Some(mapped) = self.mapped_charts.get(&(curve, on)) {
            return Ok(*mapped);
        }
        let Some(geometry) = model.geometry().pcurve(curve) else {
            ogeom_bail!(
                Dangling,
                "a chart curve read from the file is not in the model"
            );
        };
        let mut mapped = geometry.stretched(1.0, map.v_scale, tol)?;
        if map.turned {
            mapped = mapped.transformed(&turned_chart(), tol)?;
        }
        let mapped = model.geometry_mut().add_pcurve(mapped);
        self.mapped_charts.insert((curve, on), mapped);
        Ok(mapped)
    }

    /// Give each edge of a face the chart the file binds to its use there.
    ///
    /// The file keys a chart by where the face's surface sits as seen from
    /// the edge: the face's own placement, taken back through the edge's
    /// placement inside the face. Here a chart is keyed the other way round,
    /// by the edge's placement inside the face, on the face's surface
    /// already carried to where the face puts it.
    fn bind_charts(
        &mut self,
        model: &mut Model,
        (surface, location): (usize, usize),
        on: SurfaceId,
        wires: &[Shape],
        tol: Tolerances,
    ) -> OgeomResult<()> {
        let face_placement = self.transform(location)?;
        for wire in wires {
            for edge in model.children_of(wire)? {
                let within = edge.location().clone();
                let seen = within.composed(model.datums())?.inverse()? * face_placement;
                let mut found = None;
                for chart in self.charts.get(&edge.node()).map_or(&[][..], Vec::as_slice) {
                    if chart.surface == surface
                        && self.transform(chart.location)?.is_equal(&seen, tol)
                    {
                        found = Some((chart.curves, chart.range));
                        break;
                    }
                }
                let Some((curves, range)) = found else {
                    plane_chart(model, &edge, on, tol)?;
                    continue;
                };
                let representation = match curves {
                    (curve, None) => EdgeRepr::PCurve {
                        curve: self.chart_on(model, curve, on, tol)?,
                        range,
                        surface: on,
                        location: within.clone(),
                    },
                    // A seam's sides are named by the edge's sense in the
                    // face, which a turned face walks the other way.
                    (first, Some(second)) => {
                        let (first, second) = (
                            self.chart_on(model, first, on, tol)?,
                            self.chart_on(model, second, on, tol)?,
                        );
                        let (forward, reversed) = if self.turned(on) {
                            (second, first)
                        } else {
                            (first, second)
                        };
                        EdgeRepr::Seam {
                            forward,
                            reversed,
                            range,
                            surface: on,
                            location: within.clone(),
                        }
                    }
                };
                let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(TShape::data_mut) else {
                    ogeom_bail!(Construction, "a wire holds something that is not an edge");
                };
                let held = data
                    .representations
                    .iter()
                    .any(|r| r.surface() == Some(on) && r.location() == Some(&within));
                if !held {
                    data.add(representation);
                }
            }
        }
        Ok(())
    }
}

/// The shape a backward record number names.
fn pick(shapes: &[Shape], number: usize, total: usize) -> Option<Shape> {
    // Record number `n` counts back from the last: it is the one at position
    // `total - n` in file order.
    total
        .checked_sub(number)
        .and_then(|at| shapes.get(at))
        .cloned()
}

/// A transform from the format's three rows and translation column.
fn transform_of(
    rows: [[f64; 3]; 3],
    translation: [f64; 3],
    tol: Tolerances,
) -> OgeomResult<Transform> {
    // The three-by-three carries the scale inside it: its determinant is the
    // cube of the uniform factor, and dividing it out leaves the rotation.
    let m = Matrix3::new(rows);
    let determinant = m.determinant();
    if determinant.abs() <= f64::MIN_POSITIVE {
        ogeom_bail!(
            Construction,
            "a placement with a singular matrix places nothing"
        );
    }
    let scale = determinant.cbrt();
    let mut unit = [[0.0_f64; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            unit[row][column] = rows[row][column] / scale;
        }
    }
    Transform::from_parts(
        Matrix3::new(unit),
        scale,
        Vector::new(translation[0], translation[1], translation[2]),
        tol.angular().max(1e-9),
    )
}

/// A surface on a left-handed seat, rebuilt on a right-handed one; `None`
/// for one already right-handed, or of a kind with no seat.
///
/// A plane keeps its chart, its normal taken as `x` crossed with `y`, which
/// is the side the file's faces read as theirs; so does a trimmed one. A surface of revolution
/// keeps its axis, its turn running the other way: the chart is turned
/// over (`u` read as a full turn less `u`), and with it the side its faces
/// call outward, which the second value says.
fn right_handed(
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<Option<(SurfaceGeometry, bool)>> {
    use ogeom_geom::Surface as _;
    let upright = |frame: Frame| Frame::new(frame.origin(), frame.z(), frame.x(), tol);
    let left = |frame: Frame| frame.handedness() == Handedness::Left;
    let (_, v) = surface.domain();
    Ok(Some(match surface {
        SurfaceGeometry::Plane(s) if left(s.plane().frame()) => {
            let frame = s.plane().frame();
            let (u, v) = s.domain();
            let normal = Direction::new(-frame.z().vector(), tol)?;
            let frame = Frame::new(frame.origin(), normal, frame.x(), tol)?;
            (PlaneSurface::over(Plane::new(frame), u, v)?.into(), false)
        }
        SurfaceGeometry::Cylinder(s) if left(s.cylinder().frame()) => {
            let c = s.cylinder();
            let c = Cylinder::new(upright(c.frame())?, c.radius(), tol)?;
            (CylinderSurface::new(c, v)?.into(), true)
        }
        SurfaceGeometry::Cone(s) if left(s.cone().frame()) => {
            let c = s.cone();
            let c = Cone::new(
                upright(c.frame())?,
                c.reference_radius(),
                c.half_angle(),
                tol,
            )?;
            (ConeSurface::new(c, v)?.into(), true)
        }
        SurfaceGeometry::Sphere(s) if left(s.sphere().frame()) => {
            let b = s.sphere();
            (
                SphereSurface::new(Sphere::new(upright(b.frame())?, b.radius(), tol)?).into(),
                true,
            )
        }
        SurfaceGeometry::Torus(s) if left(s.torus().frame()) => {
            let t = s.torus();
            let t = Torus::new(upright(t.frame())?, t.major_radius(), t.minor_radius(), tol)?;
            (TorusSurface::new(t).into(), true)
        }
        // A trimmed plane keeps its window, which the rebuilt plane's chart
        // shares.
        SurfaceGeometry::Trimmed(t) => match right_handed(t.basis(), tol)? {
            Some((basis, false)) => {
                let (u, v) = surface.domain();
                (
                    SurfaceGeometry::Trimmed(Box::new(TrimmedSurface::new(basis, u, v, tol)?)),
                    false,
                )
            }
            _ => return Ok(None),
        },
        _ => return Ok(None),
    }))
}

/// The chart of a surface of revolution turned over: `u` read as a full
/// turn less `u`, `v` as it was.
fn turned_chart() -> Transform2 {
    Transform2::line_mirror(Point2::new(core::f64::consts::PI, 0.0), Direction2::X)
}

/// Give an edge its chart on a plane the file states none for.
///
/// The format leaves an edge on a plane without one where its writer can
/// project the curve whenever it is asked, so a reader has to do the same:
/// the image of a curve in a plane's chart is exact for every curve the
/// format carries. Anything but a plane, trimmed or not, the file must
/// chart itself.
fn plane_chart(model: &mut Model, edge: &Shape, on: SurfaceId, tol: Tolerances) -> OgeomResult<()> {
    let Some(surface) = model.geometry().surface(on).cloned() else {
        return Ok(());
    };
    let flat = match &surface {
        SurfaceGeometry::Plane(_) => true,
        SurfaceGeometry::Trimmed(t) => matches!(t.basis(), SurfaceGeometry::Plane(_)),
        _ => false,
    };
    if !flat {
        return Ok(());
    }
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        return Ok(());
    };
    let within = edge.location().clone();
    if data.degenerate
        || data
            .pcurve_for(on, &within)
            .is_some_and(|r| r.location() == Some(&within))
    {
        return Ok(());
    }
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d().cloned() else {
        return Ok(());
    };
    let Some(geometry) = model.geometry().curve(curve) else {
        ogeom_bail!(Dangling, "an edge's curve is not in the model");
    };
    // The edge's curve where the face uses it.
    let geometry = geometry.transformed(&within.composed(model.datums())?, tol)?;
    let image = match ogeom_intersect::exact_pcurve_over(&geometry, range, &surface, tol) {
        Some(exact) => crate::pcurves::over_range(&exact, range),
        None => crate::pcurves::fit_projected_pcurve(&geometry, range, &surface, tol)?.0,
    };
    let image = model.geometry_mut().add_pcurve(image);
    let Some(NodeData::Edge(data)) = model.node_mut(edge).map(TShape::data_mut) else {
        return Ok(());
    };
    data.add(EdgeRepr::PCurve {
        curve: image,
        range,
        surface: on,
        location: within,
    });
    Ok(())
}

/// The file's numbering of the edges and faces below its own shape.
fn numbered(
    model: &Model,
    built: &Built,
    shapes: &[Shape],
    root: (usize, usize),
) -> OgeomResult<Numbering> {
    let total = shapes.len();
    let mut numbering = Numbering::default();
    let mut met: FastSet<(TShapeId, Location)> = FastSet::default();
    // Each entry: a record, where its parent stands, and its own placement
    // there.
    let mut stack = vec![(root.0, Location::identity(), built.location(root.1)?)];
    while let Some((number, outer, own)) = stack.pop() {
        let (Some(index), Some(shape)) = (total.checked_sub(number), pick(shapes, number, total))
        else {
            ogeom_bail!(
                Construction,
                "a record names subshape {number}, which is not above it"
            );
        };
        let at = outer.then(&own);
        let kind = model.kind_of(&shape)?;
        if matches!(kind, ShapeType::Edge | ShapeType::Face) {
            if !met.insert((shape.node(), at.clone())) {
                continue;
            }
            if kind == ShapeType::Edge {
                numbering.edges.push(shape.moved(&at));
                continue;
            }
            // A face rebuilt where its shell placed it is the face there.
            numbering
                .faces
                .push(match built.settled.get(&(shape.node(), own)) {
                    Some(rebuilt) => Shape::of(*rebuilt).moved(&outer),
                    None => shape.moved(&at),
                });
        }
        if kind == ShapeType::Vertex {
            continue;
        }
        let Some(written) = built.written.get(index) else {
            continue;
        };
        for &(child, location) in written.iter().rev() {
            stack.push((child, at.clone(), built.location(location)?));
        }
    }
    Ok(numbering)
}

/// A face used at a placement, rebuilt where that puts it when the
/// placement undoes part of its edges' placements inside it.
///
/// An edge's chart on a face is found from the edge's placement in the
/// whole shape, whose tail is its placement inside the face. Where the
/// face's placement cancels against that (a face turned one way, its wire
/// turned back) the tail is gone from the chain; the face is then built
/// anew in place, its surface carried there and its edges' charts bound to
/// the placements they have there.
fn settled(
    model: &mut Model,
    built: &mut Built,
    face: Shape,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    if face.location().is_identity() || model.kind_of(&face)? != ShapeType::Face {
        return Ok(face);
    }
    let Some(data) = model.node(&face).and_then(|n| n.data().as_face()).cloned() else {
        return Ok(face);
    };
    let bare = Shape::of(face.node());
    let wires = model.children_of(&bare)?;
    let mut inside = Vec::new();
    for wire in &wires {
        inside.push(model.children_of(&wire.oriented(Orientation::Forward))?);
    }
    let outer = face.location().clone();
    if inside
        .iter()
        .flatten()
        .all(|edge| outer.then(edge.location()).ends_with(edge.location()))
    {
        return Ok(face);
    }
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "a face's surface is not in the model");
    };
    let surface = surface.transformed(&outer.composed(model.datums())?, tol)?;
    let surface = model.geometry_mut().add_surface(surface);
    let mut placed_wires = Vec::with_capacity(wires.len());
    for (wire, edges) in wires.iter().zip(&inside) {
        let mut placed_edges = Vec::with_capacity(edges.len());
        for edge in edges {
            let at = outer.then(edge.location());
            let chart = model
                .node(edge)
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.pcurve_for(data.surface, edge.location()))
                .cloned();
            let bound = match chart {
                Some(EdgeRepr::PCurve { curve, range, .. }) => Some(EdgeRepr::PCurve {
                    curve,
                    range,
                    surface,
                    location: at.clone(),
                }),
                Some(EdgeRepr::Seam {
                    forward,
                    reversed,
                    range,
                    ..
                }) => Some(EdgeRepr::Seam {
                    forward,
                    reversed,
                    range,
                    surface,
                    location: at.clone(),
                }),
                _ => None,
            };
            if let Some(bound) = bound
                && let Some(NodeData::Edge(edge_data)) = model.node_mut(edge).map(TShape::data_mut)
            {
                edge_data.add(bound);
            }
            placed_edges.push(
                Shape::of(edge.node())
                    .moved(&at)
                    .oriented(edge.orientation()),
            );
        }
        placed_wires.push(model.add_wire(&placed_edges)?.oriented(wire.orientation()));
    }
    let rebuilt = model.add_face(
        FaceData {
            surface,
            location: Location::identity(),
            ..data
        },
        &placed_wires,
    )?;
    built.settled.insert((face.node(), outer), rebuilt.node());
    Ok(rebuilt.oriented(face.orientation()))
}

/// A wire's edges in the order they run, each starting where the one
/// before it ends; `None` where they do not make one walk.
///
/// The format keeps a wire's edges as a set, in whatever order its writer
/// met them (a prism's side lists its two rims before the ruling between
/// them), while a wire here is a walk.
fn walked(model: &Model, edges: &[Shape]) -> OgeomResult<Option<Vec<Shape>>> {
    // Where each use of an edge starts and ends: its vertices as the use
    // meets them.
    let mut ends = Vec::with_capacity(edges.len());
    for edge in edges {
        let mut start = None;
        let mut end = None;
        for vertex in model.children_of(edge)? {
            let at = Some((vertex.node(), vertex.location().clone()));
            match vertex.orientation() {
                Orientation::Forward => start = at,
                Orientation::Reversed => end = at,
                Orientation::Internal | Orientation::External => {}
            }
        }
        let (Some(start), Some(end)) = (start, end) else {
            return Ok(None);
        };
        ends.push((start, end));
    }
    // A walk through every edge once, closed loops included wherever they
    // sit: each step leaves by an edge not yet taken from the vertex it is
    // at, and a step that finds none settles its edge into the walk
    // (Hierholzer's way), so a loop met mid-walk is spliced in rather
    // than left over.
    let mut leaving: FastMap<&(TShapeId, Location), Vec<usize>> = FastMap::default();
    for (i, (start, _)) in ends.iter().enumerate().rev() {
        leaving.entry(start).or_default().push(i);
    }
    let Some((first, _)) = ends.first() else {
        return Ok(None);
    };
    let mut stack: Vec<(&(TShapeId, Location), Option<usize>)> = vec![(first, None)];
    let mut order = Vec::with_capacity(edges.len());
    while let Some(&(at, by)) = stack.last() {
        if let Some(next) = leaving.get_mut(at).and_then(Vec::pop) {
            stack.push((&ends[next].1, Some(next)));
        } else {
            stack.pop();
            order.extend(by);
        }
    }
    order.reverse();
    if order.len() != edges.len() {
        return Ok(None);
    }
    if order.iter().enumerate().all(|(at, &i)| at == i) {
        return Ok(None);
    }
    Ok(Some(order.into_iter().map(|i| edges[i].clone()).collect()))
}

/// Place a shape by the file's placement number.
fn placed(shape: Shape, built: &Built, number: usize) -> OgeomResult<Shape> {
    let location = built.location(number)?;
    if location.is_identity() {
        return Ok(shape);
    }
    Ok(shape.moved(&location))
}

/// A word, or a number and the continuity written straight after it.
///
/// A seam's record runs its second chart number into its continuity with no
/// space between (`72CN` for chart 72, continuity `CN`), so the two are
/// told apart by the continuity names, the only words that can follow.
fn split_continuity(word: &str) -> impl Iterator<Item = &str> {
    let split = ["C0", "C1", "C2", "C3", "CN", "G1", "G2"]
        .iter()
        .find_map(|name| {
            let number = word.strip_suffix(name)?;
            (!number.is_empty() && number.bytes().all(|b| b.is_ascii_digit()))
                .then(|| word.split_at(number.len()))
        });
    let (first, second) = split.map_or((word, None), |(a, b)| (a, Some(b)));
    std::iter::once(first).chain(second)
}

/// A whitespace-delimited walk over the file's body.
struct Cursor<'a> {
    tokens: Vec<&'a str>,
    at: usize,
}

impl<'a> Cursor<'a> {
    fn new(text: &'a str) -> Self {
        Self {
            tokens: text.split_whitespace().flat_map(split_continuity).collect(),
            at: 0,
        }
    }

    fn word(&mut self) -> OgeomResult<&'a str> {
        let Some(token) = self.tokens.get(self.at) else {
            ogeom_bail!(Construction, "the file ends in the middle of a record");
        };
        self.at += 1;
        Ok(token)
    }

    fn peek(&self) -> Option<&'a str> {
        self.tokens.get(self.at).copied()
    }

    fn expect(&mut self, what: &str) -> OgeomResult<()> {
        let token = self.word()?;
        if token != what {
            ogeom_bail!(Construction, "expected the {what} table, found `{token}`");
        }
        Ok(())
    }

    fn number(&mut self) -> OgeomResult<f64> {
        let token = self.word()?;
        let Ok(value) = token.parse::<f64>() else {
            ogeom_bail!(Construction, "`{token}` is not a number");
        };
        Ok(value)
    }

    fn integer(&mut self) -> OgeomResult<i64> {
        let token = self.word()?;
        let Ok(value) = token.parse::<i64>() else {
            ogeom_bail!(Construction, "`{token}` is not a whole number");
        };
        Ok(value)
    }

    /// Room to reserve for `n` items read from here on: never more than the
    /// words left, since every item takes at least one, so a count in the
    /// file cannot size an allocation past what the file holds.
    fn room(&self, n: usize) -> usize {
        n.min(self.tokens.len().saturating_sub(self.at))
    }

    fn count(&mut self) -> OgeomResult<usize> {
        let value = self.integer()?;
        let Ok(count) = usize::try_from(value) else {
            ogeom_bail!(Construction, "{value} is not a count");
        };
        Ok(count)
    }

    fn flag(&mut self) -> OgeomResult<bool> {
        Ok(self.count()? != 0)
    }

    fn point(&mut self) -> OgeomResult<Point> {
        Ok(Point::new(self.number()?, self.number()?, self.number()?))
    }

    fn point2(&mut self) -> OgeomResult<Point2> {
        Ok(Point2::new(self.number()?, self.number()?))
    }

    fn direction(&mut self, tol: Tolerances) -> OgeomResult<Direction> {
        let v = Vector::new(self.number()?, self.number()?, self.number()?);
        Direction::new(v, tol)
    }

    fn direction2(&mut self, tol: Tolerances) -> OgeomResult<Direction2> {
        let v = Vector2::new(self.number()?, self.number()?);
        Direction2::new(v, tol)
    }

    /// The origin-normal-x-y seat every elementary record opens with.
    ///
    /// A surface's seat may be left-handed, its `y` against the normal
    /// crossed with `x`: the chart then runs the other way across the
    /// surface, and the side its faces call outward turns with it, so the
    /// `y` the file states is kept rather than derived.
    fn seat(&mut self, tol: Tolerances) -> OgeomResult<Frame> {
        let origin = self.point()?;
        let normal = self.direction(tol)?;
        let x = self.direction(tol)?;
        let y = self.direction(tol)?;
        Frame::from_axes(origin, x, y, normal, tol)
    }

    fn seat2(&mut self, tol: Tolerances) -> OgeomResult<Frame2> {
        let origin = self.point2()?;
        let x = self.direction2(tol)?;
        let y = self.direction2(tol)?;
        Frame2::from_axes(origin, x, y, tol)
    }

    /// Poles and knots, as both spline records spell them.
    fn spline_knots(&mut self, degree: usize, count: usize) -> OgeomResult<KnotVector> {
        let mut flat = Vec::new();
        for _ in 0..count {
            let value = self.number()?;
            let multiplicity = self.count()?;
            // No valid spline repeats a knot more often than its order, and
            // no real one is of a degree past 64.
            if multiplicity > degree.min(64) + 1 {
                ogeom_bail!(
                    Construction,
                    "a knot of multiplicity {multiplicity} on a degree {degree} spline"
                );
            }
            flat.extend(std::iter::repeat_n(value, multiplicity));
        }
        KnotVector::new(flat, degree)
    }

    fn curve(&mut self, tol: Tolerances) -> OgeomResult<Curve> {
        let kind = self.count()?;
        Ok(match kind {
            1 => {
                let origin = self.point()?;
                let direction = self.direction(tol)?;
                // The record carries no window; the edge's own range is what
                // bounds it, so the curve is built wide and trimmed there.
                LineCurve::over(Axis::new(origin, direction), -OPEN_EXTENT, OPEN_EXTENT)?.into()
            }
            2 => {
                let frame = self.seat(tol)?;
                CircleCurve::new(Circle::new(frame, self.number()?, tol)?).into()
            }
            3 => {
                let frame = self.seat(tol)?;
                let ellipse = Ellipse::new(frame, self.number()?, self.number()?, tol)?;
                EllipseCurve::new(ellipse).into()
            }
            4 => {
                let frame = self.seat(tol)?;
                let parabola = Parabola::new(frame, self.number()?, tol)?;
                ParabolaCurve::new(parabola, OPEN_EXTENT)?.into()
            }
            5 => {
                let frame = self.seat(tol)?;
                let hyperbola = Hyperbola::new(frame, self.number()?, self.number()?, tol)?;
                HyperbolaCurve::new(hyperbola, OPEN_EXTENT)?.into()
            }
            7 => {
                let rational = self.flag()?;
                let _periodic = self.flag()?;
                let degree = self.count()?;
                let poles = self.count()?;
                let knot_count = self.count()?;
                let mut control = Vec::with_capacity(self.room(poles));
                for _ in 0..poles {
                    let at = self.point()?;
                    let weight = if rational { self.number()? } else { 1.0 };
                    control.push(Weighted::new(at, weight, tol)?);
                }
                let knots = self.spline_knots(degree, knot_count)?;
                BSplineCurve::rational(knots, control)?.into()
            }
            8 => {
                let lo = self.number()?;
                let hi = self.number()?;
                let basis = self.curve(tol)?;
                Curve::Trimmed(Box::new(TrimmedCurve::new(basis, lo, hi, tol)?))
            }
            9 => {
                let distance = self.number()?;
                let reference = self.direction(tol)?;
                let basis = self.curve(tol)?;
                Curve::Offset(Box::new(ogeom_geom::OffsetCurve::new(
                    basis, distance, reference,
                )?))
            }
            other => ogeom_bail!(
                Construction,
                "curve record {other} is one this does not read; a Bezier \
                 curve is record 6, and nothing here builds one"
            ),
        })
    }

    fn pcurve(&mut self, tol: Tolerances) -> OgeomResult<PlanarCurve> {
        let kind = self.count()?;
        Ok(match kind {
            1 => {
                let origin = self.point2()?;
                let direction = self.direction2(tol)?;
                Line2d::over(Axis2::new(origin, direction), -OPEN_EXTENT, OPEN_EXTENT)?.into()
            }
            2 => {
                let frame = self.seat2(tol)?;
                Circle2d::new(Circle2::new(frame, self.number()?, tol)?).into()
            }
            3 => {
                let frame = self.seat2(tol)?;
                let ellipse = Ellipse2::new(frame, self.number()?, self.number()?, tol)?;
                Ellipse2d::new(ellipse).into()
            }
            5 => {
                let frame = self.seat2(tol)?;
                let hyperbola = Hyperbola2::new(frame, self.number()?, self.number()?, tol)?;
                ogeom_geom::Hyperbola2d::over(hyperbola, -OPEN_EXTENT, OPEN_EXTENT)?.into()
            }
            7 => {
                let rational = self.flag()?;
                let _periodic = self.flag()?;
                let degree = self.count()?;
                let poles = self.count()?;
                let knot_count = self.count()?;
                let mut control = Vec::with_capacity(self.room(poles));
                for _ in 0..poles {
                    let at = self.point2()?;
                    let weight = if rational { self.number()? } else { 1.0 };
                    control.push(Weighted::new(at, weight, tol)?);
                }
                let knots = self.spline_knots(degree, knot_count)?;
                ogeom_geom::BSpline2d::rational(knots, control)?.into()
            }
            8 => {
                let lo = self.number()?;
                let hi = self.number()?;
                let basis = self.pcurve(tol)?;
                PlanarCurve::Trimmed(Box::new(ogeom_geom::Trimmed2d::new(basis, lo, hi, tol)?))
            }
            9 => {
                let distance = self.number()?;
                let basis = self.pcurve(tol)?;
                PlanarCurve::Offset(Box::new(ogeom_geom::Offset2d::new(basis, distance)?))
            }
            other => ogeom_bail!(
                Construction,
                "chart curve record {other} is one this does not read"
            ),
        })
    }

    fn surface(&mut self, tol: Tolerances) -> OgeomResult<SurfaceGeometry> {
        let kind = self.count()?;
        // An unbounded carrier needs a window to live in; the format gives
        // none, so one wide enough for anything the file can hold is used,
        // and the face's own trim is what bounds it in practice.
        const REACH: (f64, f64) = (-1e9, 1e9);
        Ok(match kind {
            1 => {
                let frame = self.seat(tol)?;
                PlaneSurface::over(Plane::new(frame), REACH, REACH)?.into()
            }
            2 => {
                let frame = self.seat(tol)?;
                let cylinder = Cylinder::new(frame, self.number()?, tol)?;
                CylinderSurface::new(cylinder, REACH)?.into()
            }
            3 => {
                let frame = self.seat(tol)?;
                let cone = Cone::new(frame, self.number()?, self.number()?, tol)?;
                ConeSurface::new(cone, REACH)?.into()
            }
            4 => {
                let frame = self.seat(tol)?;
                SphereSurface::new(Sphere::new(frame, self.number()?, tol)?).into()
            }
            5 => {
                let frame = self.seat(tol)?;
                let torus = Torus::new(frame, self.number()?, self.number()?, tol)?;
                TorusSurface::new(torus).into()
            }
            6 => {
                let direction = self.direction(tol)?;
                let basis = self.curve(tol)?;
                ExtrusionSurface::over(basis, direction, REACH)?.into()
            }
            7 => {
                let origin = self.point()?;
                let direction = self.direction(tol)?;
                let basis = self.curve(tol)?;
                RevolutionSurface::new(basis, Axis::new(origin, direction), core::f64::consts::TAU)?
                    .into()
            }
            9 => {
                let u_rational = self.flag()?;
                let v_rational = self.flag()?;
                let _u_periodic = self.flag()?;
                let _v_periodic = self.flag()?;
                let u_degree = self.count()?;
                let v_degree = self.count()?;
                let u_poles = self.count()?;
                let v_poles = self.count()?;
                let u_knot_count = self.count()?;
                let v_knot_count = self.count()?;
                let rational = u_rational || v_rational;
                let Some(cells) = u_poles.checked_mul(v_poles) else {
                    ogeom_bail!(Construction, "a {u_poles} by {v_poles} pole grid");
                };
                let mut points = Vec::with_capacity(self.room(cells));
                for _ in 0..cells {
                    let at = self.point()?;
                    let weight = if rational { self.number()? } else { 1.0 };
                    points.push(Weighted::new(at, weight, tol)?);
                }
                let u_knots = self.spline_knots(u_degree, u_knot_count)?;
                let v_knots = self.spline_knots(v_degree, v_knot_count)?;
                let grid = ControlGrid::new(points, u_poles, v_poles)?;
                BSplineSurface::rational(u_knots, v_knots, grid)?.into()
            }
            10 => {
                let u = (self.number()?, self.number()?);
                let v = (self.number()?, self.number()?);
                let basis = self.surface(tol)?;
                SurfaceGeometry::Trimmed(Box::new(TrimmedSurface::new(basis, u, v, tol)?))
            }
            11 => {
                let distance = self.number()?;
                let basis = self.surface(tol)?;
                SurfaceGeometry::Offset(Box::new(ogeom_geom::OffsetSurface::new(basis, distance)?))
            }
            other => ogeom_bail!(
                Construction,
                "surface record {other} is one this does not read; a Bezier \
                 surface is record 8, and nothing here builds one"
            ),
        })
    }

    /// A subshape reference: its orientation, its backward record number and
    /// its placement.
    fn reference(&mut self) -> OgeomResult<(Orientation, usize, usize)> {
        let token = self.word()?;
        let (orientation, digits) = match token.split_at(1) {
            ("+", rest) => (Orientation::Forward, rest),
            ("-", rest) => (Orientation::Reversed, rest),
            ("i", rest) => (Orientation::Internal, rest),
            ("e", rest) => (Orientation::External, rest),
            _ => ogeom_bail!(Construction, "`{token}` does not name a subshape"),
        };
        let Ok(number) = digits.parse::<usize>() else {
            ogeom_bail!(Construction, "`{token}` does not name a subshape number");
        };
        Ok((orientation, number, self.count()?))
    }

    fn skip_polygon3d(&mut self) -> OgeomResult<()> {
        let nodes = self.count()?;
        let has_parameters = self.flag()?;
        let _deflection = self.number()?;
        for _ in 0..nodes * 3 {
            let _ = self.number()?;
        }
        if has_parameters {
            for _ in 0..nodes {
                let _ = self.number()?;
            }
        }
        Ok(())
    }

    fn skip_polygon_on_triangulation(&mut self) -> OgeomResult<()> {
        let nodes = self.count()?;
        for _ in 0..nodes {
            let _ = self.integer()?;
        }
        // An optional parameter list follows, introduced by `p`.
        if self.peek() == Some("p") {
            let _ = self.word()?;
            let _deflection = self.number()?;
            let _flag = self.count()?;
            let count = self.count()?;
            for _ in 0..count {
                let _ = self.number()?;
            }
        }
        Ok(())
    }

    fn skip_triangulation(&mut self) -> OgeomResult<()> {
        let nodes = self.count()?;
        let triangles = self.count()?;
        let has_parameters = self.flag()?;
        let _deflection = self.number()?;
        for _ in 0..nodes * 3 {
            let _ = self.number()?;
        }
        if has_parameters {
            for _ in 0..nodes * 2 {
                let _ = self.number()?;
            }
        }
        for _ in 0..triangles * 3 {
            let _ = self.integer()?;
        }
        Ok(())
    }

    /// One topology record, built into the model against what came before.
    fn record(
        &mut self,
        model: &mut Model,
        built: &mut Built,
        shapes: &[Shape],
        total: usize,
        tol: Tolerances,
    ) -> OgeomResult<Shape> {
        let kind = self.word()?;
        let data = match kind {
            "Ve" => {
                let tolerance = self.number()?;
                let at = self.point()?;
                // The parameter list, which this recovers by projection where
                // it needs it rather than trusting a file's copy.
                loop {
                    let first = self.word()?;
                    let second = self.word()?;
                    if first == "0" && second == "0" {
                        break;
                    }
                    // A representation: its data depends on the kind in
                    // `second`, and every one ends with a placement number.
                    match second {
                        "1" => {
                            let _ = self.count()?;
                        }
                        "2" | "3" => {
                            let _ = self.count()?;
                            let _ = self.count()?;
                        }
                        other => ogeom_bail!(
                            Construction,
                            "a vertex representation of kind {other} is not one this reads"
                        ),
                    }
                    let _ = self.count()?;
                }
                Record::Vertex(VertexData {
                    point: at,
                    tolerance: Tolerance::new(tolerance.max(tol.confusion()))?,
                })
            }
            "Ed" => {
                let tolerance = self.number()?;
                let _same_parameter = self.flag()?;
                let _same_range = self.flag()?;
                let degenerate = self.flag()?;
                let mut data = EdgeData::new();
                data.tolerance = Tolerance::new(tolerance.max(tol.confusion()))?;
                data.degenerate = degenerate;
                let mut charts = Vec::new();
                loop {
                    let kind = self.count()?;
                    match kind {
                        0 => break,
                        1 => {
                            let curve = self.count()?;
                            let location = self.count()?;
                            let lo = self.number()?;
                            let hi = self.number()?;
                            let id = built.curve(model, curve, location, tol)?;
                            data.add(EdgeRepr::Curve3d {
                                curve: id,
                                range: (lo, hi),
                                location: Location::identity(),
                            });
                        }
                        2 | 3 => {
                            let first = self.count()?;
                            let second = if kind == 3 {
                                let reversed = self.count()?;
                                let _continuity = self.word()?;
                                Some(reversed)
                            } else {
                                None
                            };
                            let surface = self.count()?;
                            let location = self.count()?;
                            let lo = self.number()?;
                            let hi = self.number()?;
                            let chart = |number: usize| {
                                number
                                    .checked_sub(1)
                                    .and_then(|at| built.pcurves.get(at))
                                    .copied()
                            };
                            let (Some(curve), true, true) = (
                                chart(first),
                                second.is_none_or(|n| chart(n).is_some()),
                                (1..=built.surfaces.len()).contains(&surface),
                            ) else {
                                ogeom_bail!(
                                    Construction,
                                    "an edge names chart curve {first} on surface {surface}"
                                );
                            };
                            charts.push(StatedChart {
                                curves: (curve, second.and_then(chart)),
                                surface,
                                location,
                                range: (lo, hi),
                            });
                        }
                        4 => {
                            let _continuity = self.word()?;
                            for _ in 0..4 {
                                let _ = self.count()?;
                            }
                        }
                        5 => {
                            let _ = self.count()?;
                            let _ = self.count()?;
                        }
                        6 => {
                            for _ in 0..3 {
                                let _ = self.count()?;
                            }
                        }
                        7 => {
                            for _ in 0..4 {
                                let _ = self.count()?;
                            }
                        }
                        other => ogeom_bail!(
                            Construction,
                            "an edge representation of kind {other} is not one this reads"
                        ),
                    }
                }
                // The file states whether its representations agree on
                // parameterization. That is the writer's claim about its own
                // data, and it is checked here rather than believed: once the
                // faces have bound the edge's charts, the representations are
                // evaluated against each other, and the claim is only
                // re-established when they actually agree.
                data.assert_same_parameter(agree_on_parameter(model, &data, tol));
                Record::Edge(Box::new(data), charts)
            }
            "Fa" => {
                let _natural = self.flag()?;
                let tolerance = self.number()?;
                let surface = self.count()?;
                let location = self.count()?;
                if self.peek() == Some("2") {
                    let _ = self.word()?;
                    let _triangulation = self.count()?;
                }
                let on = built.surface(model, surface.max(1), location, tol)?;
                Record::Face(
                    Box::new(FaceData {
                        surface: on,
                        location: Location::identity(),
                        tolerance: Tolerance::new(tolerance.max(tol.confusion()))?,
                        natural_restriction: false,
                        triangulation: None,
                    }),
                    (surface.max(1), location),
                )
            }
            "Wi" => Record::Container(ShapeType::Wire),
            "Sh" => Record::Container(ShapeType::Shell),
            "So" => Record::Container(ShapeType::Solid),
            "CS" => Record::Container(ShapeType::CompSolid),
            "Co" => Record::Container(ShapeType::Compound),
            other => ogeom_bail!(Construction, "`{other}` is not a shape record this reads"),
        };

        // The flag word, which says nothing about the shape, then the
        // children.
        let _flags = self.word()?;
        let mut children = Vec::new();
        let mut written = Vec::new();
        while self.peek() != Some("*") {
            let (orientation, number, location) = self.reference()?;
            written.push((number, location));
            let Some(child) = pick(shapes, number, total) else {
                ogeom_bail!(
                    Construction,
                    "a record names subshape {number}, which is not above it"
                );
            };
            let sense = built.sense(&child, orientation);
            children.push(placed(child, built, location)?.composed(sense));
        }
        let _ = self.word()?;
        built.written.push(written);

        Ok(match data {
            Record::Vertex(data) => model.add_vertex(data),
            Record::Edge(data, charts) => {
                // The file may list an edge's end before its start; an edge
                // here holds its start first.
                children.sort_by_key(|v| v.orientation() != Orientation::Forward);
                let edge = model.add_edge(*data, &children)?;
                if !charts.is_empty() {
                    built.charts.insert(edge.node(), charts);
                }
                edge
            }
            Record::Face(data, stated) => {
                // A wire the file places inside its face is rebuilt in the
                // face's own frame, its edges carrying the placement: what
                // reads a face's boundary takes the wire's frame for the
                // face's.
                let mut wires = Vec::with_capacity(children.len());
                for wire in children {
                    if wire.location().is_identity() {
                        wires.push(wire);
                        continue;
                    }
                    let edges = model.children_of(&wire.oriented(Orientation::Forward))?;
                    wires.push(model.add_wire(&edges)?.oriented(wire.orientation()));
                }
                built.bind_charts(model, stated, data.surface, &wires, tol)?;
                if built.turned(data.surface) {
                    let wires: Vec<Shape> = wires.iter().map(Shape::reversed).collect();
                    let face = model.add_face(*data, &wires)?;
                    built.turned_faces.insert(face.node());
                    face
                } else {
                    model.add_face(*data, &wires)?
                }
            }
            Record::Container(ShapeType::Wire) => {
                let walk = walked(model, &children)?;
                model.add_wire(walk.as_deref().unwrap_or(&children))?
            }
            Record::Container(ShapeType::Shell) => {
                let faces = children
                    .into_iter()
                    .map(|face| settled(model, built, face, tol))
                    .collect::<OgeomResult<Vec<_>>>()?;
                model.add_shell(&faces)?
            }
            Record::Container(ShapeType::Solid) => model.add_solid(&children)?,
            Record::Container(ShapeType::CompSolid) => model.add_compsolid(&children)?,
            Record::Container(_) => model.add_compound(&children)?,
        })
    }
}

/// What a record turned out to be, before its children are known.
enum Record {
    Vertex(VertexData),
    /// The edge, with the charts it states for the faces to bind.
    Edge(Box<EdgeData>, Vec<StatedChart>),
    /// The face, with its surface and placement as the file numbers them.
    Face(Box<FaceData>, (usize, usize)),
    Container(ShapeType),
}

/// Whether an edge's representations land on the same point at the same
/// parameter, which is what the same-parameter claim means.
///
/// Sampled: nine stations across the edge's own range, each read through
/// every representation it carries, all held to the edge's stated tolerance.
/// A file whose pcurve was fitted independently of its curve fails this, and
/// it should: the claim is what every algorithm downstream relies on.
fn agree_on_parameter(model: &Model, data: &EdgeData, tol: Tolerances) -> bool {
    use ogeom_geom::Surface as _;
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        return false;
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        return false;
    };
    let allowed = data.tolerance.get().max(tol.confusion());
    for representation in &data.representations {
        // A chart bound to the edge's use at a placement follows the curve
        // carried there.
        let Ok(moved) = representation
            .location()
            .map_or(Ok(Transform::IDENTITY), |l| l.composed(model.datums()))
        else {
            return false;
        };
        let (pcurve, prange, surface) = match representation {
            EdgeRepr::PCurve {
                curve,
                range,
                surface,
                ..
            } => (*curve, *range, *surface),
            EdgeRepr::Seam {
                forward,
                range,
                surface,
                ..
            } => (*forward, *range, *surface),
            _ => continue,
        };
        let (Some(chart), Some(on)) = (
            model.geometry().pcurve(pcurve),
            model.geometry().surface(surface),
        ) else {
            return false;
        };
        for k in 0..=8 {
            let f = f64::from(k) / 8.0;
            let t = (range.1 - range.0).mul_add(f, range.0);
            let pt = (prange.1 - prange.0).mul_add(f, prange.0);
            let (Ok(at), Ok(uv)) = (geometry.point_at(t, tol), chart.point_at(pt, tol)) else {
                return false;
            };
            let Ok(lifted) = on.point_at(uv.x, uv.y, tol) else {
                return false;
            };
            if lifted.distance(moved.apply(at)) > allowed {
                return false;
            }
        }
    }
    true
}
