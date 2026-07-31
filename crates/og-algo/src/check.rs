//! Checking a shape against the model's invariants.
//!
//! Everything in `docs/DATA_MODEL.md` that can be checked with geometry in
//! hand, checked in one place. The builders enforce what they can at the moment
//! of construction; this catches what only becomes wrong later — an edge whose
//! tolerance was widened past its face's, a shell left open by an operation
//! that dropped a face, a pcurve that has stopped agreeing with its curve.
//!
//! # It reports, it does not judge
//!
//! The result is a list of what is wrong and where, not a boolean. A boolean
//! answers "should I panic", which is never the question: an imported shape is
//! usually invalid in some specific, fixable way, and healing it needs to know
//! which way. A caller that only wants the boolean asks
//! [`Diagnosis::is_valid`].
//!
//! # Severity is not a comment
//!
//! [`Severity::Broken`] means an algorithm reading this shape will get a wrong
//! answer rather than an error — an open shell has no inside, so every
//! containment test against it is a coin toss. [`Severity::Suspect`] means
//! something is out of order but every operation will still behave: a tolerance
//! larger than the feature it describes is alarming and not yet wrong.
//!
//! The distinction is what lets a pipeline decide. Booleans refuse `Broken`
//! input because they would produce nonsense from it; they proceed on
//! `Suspect` because refusing would reject most real imported geometry.

use std::collections::HashMap;
use std::fmt;

use og_core::{OgResult, Tolerances, og_bail};
use og_geom::{Curve2d, Curve3d, Surface};
use og_mesh::Deflection;
use og_topo::{EdgeRepr, Filter, Model, Shape, ShapeType, TShapeId, explore, explore_unique};

/// How badly a problem breaks the shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Severity {
    /// Out of order, but every operation will still behave.
    Suspect,
    /// An algorithm reading this shape will get a wrong answer, not an error.
    Broken,
}

/// One thing wrong with a shape.
#[derive(Debug, Clone, PartialEq)]
pub struct Problem {
    /// How badly it breaks things.
    pub severity: Severity,
    /// The sub-shape it is about.
    pub at: Shape,
    /// What kind of sub-shape that is, so a report reads without a lookup.
    pub kind: ShapeType,
    /// What is wrong, in a sentence.
    pub what: String,
}

impl fmt::Display for Problem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mark = match self.severity {
            Severity::Broken => "broken",
            Severity::Suspect => "suspect",
        };
        write!(f, "[{mark}] {:?}: {}", self.kind, self.what)
    }
}

/// Everything wrong with a shape.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Diagnosis {
    /// The problems found, in the order they were found.
    pub problems: Vec<Problem>,
}

impl Diagnosis {
    /// Whether nothing is wrong at all.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.problems.is_empty()
    }

    /// Whether anything would make an algorithm answer wrongly.
    ///
    /// The question a boolean or a mass property should ask before starting.
    #[must_use]
    pub fn is_usable(&self) -> bool {
        !self.problems.iter().any(|p| p.severity == Severity::Broken)
    }

    /// The worst severity found.
    #[must_use]
    pub fn worst(&self) -> Option<Severity> {
        self.problems.iter().map(|p| p.severity).max()
    }

    /// The problems of one severity.
    #[must_use]
    pub fn of(&self, severity: Severity) -> Vec<&Problem> {
        self.problems
            .iter()
            .filter(|p| p.severity == severity)
            .collect()
    }

    fn note(&mut self, severity: Severity, at: &Shape, kind: ShapeType, what: String) {
        self.problems.push(Problem {
            severity,
            at: at.clone(),
            kind,
            what,
        });
    }
}

impl fmt::Display for Diagnosis {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.problems.is_empty() {
            return write!(f, "valid");
        }
        for (i, problem) in self.problems.iter().enumerate() {
            if i > 0 {
                writeln!(f)?;
            }
            write!(f, "{problem}")?;
        }
        Ok(())
    }
}

/// Check a shape and everything below it.
///
/// # Errors
///
/// [`OgError::Dangling`](og_core::OgError::Dangling) if a handle does not
/// resolve. A dangling handle is not a *finding* — it means the shape and the
/// model do not belong together, and every other answer would be about
/// something that is not there.
pub fn check(model: &Model, shape: &Shape, tol: Tolerances) -> OgResult<Diagnosis> {
    if model.node(shape).is_none() {
        og_bail!(Dangling, "shape refers to a node not in this model");
    }
    let mut found = Diagnosis::default();

    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        check_edge(model, &edge, tol, &mut found)?;
    }
    for wire in explore_unique(model, shape, ShapeType::Wire)? {
        check_wire(model, &wire, tol, &mut found)?;
    }
    for face in explore_unique(model, shape, ShapeType::Face)? {
        check_face(model, &face, tol, &mut found)?;
    }
    for shell in explore_unique(model, shape, ShapeType::Shell)? {
        check_shell(model, &shell, &mut found)?;
    }
    check_containment(model, shape, &mut found)?;
    Ok(found)
}

/// Check that a shape's *tessellation* agrees with its topology.
///
/// Separate from [`check`] because it needs a deflection, and because it asks a
/// different question: not "is this shape well formed" but "do its two
/// descriptions of itself agree". A shell whose edges are all used twice is
/// closed as far as the topology knows. If the mesh built from it still has a
/// boundary, then some face's pcurves do not cover the region its edges claim
/// to bound — and the topology cannot see that, because the defect is entirely
/// in parameter space.
///
/// That failure is worth its own function because it is *invisible* to every
/// other check. Face counts look right, the shell closes, each face
/// triangulates without error, and the solid still has a slit down it. The
/// first thing to notice is usually a volume that is quietly wrong.
///
/// Reports the position of the unshared edges, not just their number: where the
/// mesh comes apart is the whole diagnosis, and a count sends you looking.
///
/// # Errors
///
/// As [`og_mesh::triangulate()`].
pub fn check_tessellation(
    model: &Model,
    shape: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgResult<Diagnosis> {
    let mut found = Diagnosis::default();

    for shell in explore_unique(model, shape, ShapeType::Shell)? {
        // An open shell is *meant* to have a boundary, so a mesh with one is
        // agreement, not disagreement. Only a shell the topology calls closed
        // makes a claim the mesh can contradict.
        if !crate::build::is_shell_closed(model, &shell)? {
            continue;
        }
        let mesh = og_mesh::triangulate(model, &shell, deflection, tol)?;
        if mesh.is_empty() {
            found.note(
                Severity::Broken,
                &shell,
                ShapeType::Shell,
                "the topology says this shell is closed and it tessellates to \
                 nothing at all"
                    .into(),
            );
            continue;
        }
        if let Some(report) = open_edges(&mesh) {
            found.note(Severity::Broken, &shell, ShapeType::Shell, report);
        }
    }
    Ok(found)
}

/// Describe a mesh's unshared edges, or `None` if every edge is shared twice.
fn open_edges(mesh: &og_topo::Triangulation) -> Option<String> {
    let mut uses: HashMap<(u32, u32), usize> = HashMap::new();
    for triangle in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (triangle[i], triangle[(i + 1) % 3]);
            *uses.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }

    let mut loose: Vec<&(u32, u32)> = uses
        .iter()
        .filter(|(_, n)| **n != 2)
        .map(|(e, _)| e)
        .collect();
    if loose.is_empty() {
        return None;
    }
    // Deterministic: a diagnosis that names a different edge each run is one
    // nobody can act on.
    loose.sort_unstable();

    let sample: Vec<String> = loose
        .iter()
        .take(3)
        .map(|(a, b)| {
            let (p, q) = (mesh.positions[*a as usize], mesh.positions[*b as usize]);
            format!(
                "({:.6}, {:.6}, {:.6})-({:.6}, {:.6}, {:.6})",
                p.x, p.y, p.z, q.x, q.y, q.z
            )
        })
        .collect();

    Some(format!(
        "the topology says this shell is closed, but its mesh has {} triangle \
         edge(s) not shared by two triangles, so the tessellated solid has a \
         slit in it. The first are at {}. This is a parameter-space defect — \
         some face's pcurves do not cover the region its edges bound — and no \
         topological check can see it",
        loose.len(),
        sample.join(", ")
    ))
}

/// An edge's curve must reach the vertices it claims to join, and its
/// representations must agree with each other.
fn check_edge(model: &Model, edge: &Shape, tol: Tolerances, found: &mut Diagnosis) -> OgResult<()> {
    let Some(node) = model.node(edge) else {
        og_bail!(Dangling, "edge is not in this model");
    };
    let Some(data) = node.data().as_edge() else {
        return Ok(());
    };
    let reach = data.tolerance.get().max(tol.confusion());

    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        // A degenerate edge is *supposed* to have no curve. Any other edge
        // without one has nowhere in space it runs, and every algorithm that
        // walks a boundary will skip it silently.
        if !data.degenerate {
            found.note(
                Severity::Broken,
                edge,
                ShapeType::Edge,
                "no curve in space, and not marked degenerate; a boundary walk \
                 will step over it without noticing"
                    .into(),
            );
        }
        return Ok(());
    };
    if data.degenerate {
        found.note(
            Severity::Suspect,
            edge,
            ShapeType::Edge,
            "marked degenerate but carries a curve in space".into(),
        );
    }
    let Some(geometry) = model.geometry().curve(*curve) else {
        og_bail!(Dangling, "curve is not in this model");
    };

    let placement = edge.transform(model.datums())?;
    let bounds = model.children_of(edge)?;
    for (parameter, vertex) in [(range.0, bounds.first()), (range.1, bounds.last())] {
        let Some(vertex) = vertex else { continue };
        let Some(point) = model
            .node(vertex)
            .and_then(|n| n.data().as_vertex())
            .map(|v| v.point)
        else {
            continue;
        };
        let placed = vertex.transform(model.datums())?.apply(point);
        let on_curve = placement.apply(geometry.point_at(parameter, tol)?);
        let gap = on_curve.distance(placed);
        if gap > reach {
            found.note(
                Severity::Broken,
                edge,
                ShapeType::Edge,
                format!(
                    "curve stops {gap} from the vertex it should meet, outside \
                     its tolerance of {reach}; the boundary has a gap there"
                ),
            );
        }
    }

    // `same_parameter` is a claim, and a false one is worse than no claim: every
    // algorithm evaluates whichever representation is cheapest and assumes the
    // answer is interchangeable.
    if data.same_parameter() {
        check_same_parameter(model, edge, data, geometry, *range, reach, tol, found)?;
    }
    Ok(())
}

/// Verify that every pcurve lands where the 3D curve does.
#[allow(clippy::too_many_arguments)]
fn check_same_parameter(
    model: &Model,
    edge: &Shape,
    data: &og_topo::EdgeData,
    curve: &og_geom::Curve,
    range: (f64, f64),
    reach: f64,
    tol: Tolerances,
    found: &mut Diagnosis,
) -> OgResult<()> {
    const SAMPLES: usize = 8;
    for repr in &data.representations {
        let (pcurve_id, pcurve_range, surface_id) = match repr {
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
        let (Some(pcurve), Some(surface)) = (
            model.geometry().pcurve(pcurve_id),
            model.geometry().surface(surface_id),
        ) else {
            og_bail!(Dangling, "an edge names geometry not in this model");
        };

        for i in 0..=SAMPLES {
            #[allow(clippy::cast_precision_loss)]
            let t = i as f64 / SAMPLES as f64;
            let on_curve = curve.point_at(range.0 + (range.1 - range.0) * t, tol)?;
            let at =
                pcurve.point_at(pcurve_range.0 + (pcurve_range.1 - pcurve_range.0) * t, tol)?;
            let Ok(on_surface) = surface.point_at(at.x, at.y, tol) else {
                continue;
            };
            let gap = on_curve.distance(on_surface);
            if gap > reach {
                found.note(
                    Severity::Broken,
                    edge,
                    ShapeType::Edge,
                    format!(
                        "claims same_parameter but its pcurve is {gap} from its \
                         curve at parameter {t} of the range, outside the edge's \
                         tolerance of {reach}"
                    ),
                );
                break;
            }
        }
    }
    Ok(())
}

/// A wire's edges must meet end to end.
fn check_wire(model: &Model, wire: &Shape, tol: Tolerances, found: &mut Diagnosis) -> OgResult<()> {
    let edges = model.ordered_children_of(wire)?;
    if edges.is_empty() {
        found.note(
            Severity::Broken,
            wire,
            ShapeType::Wire,
            "has no edges, so it bounds nothing".into(),
        );
        return Ok(());
    }
    for i in 0..edges.len() {
        let (Some((_, end)), Some((next, _))) = (
            crate::build::edge_vertices(model, &edges[i])?,
            crate::build::edge_vertices(model, &edges[(i + 1) % edges.len()])?,
        ) else {
            found.note(
                Severity::Broken,
                wire,
                ShapeType::Wire,
                format!("edge {i} has no bounding vertices, so it joins nothing"),
            );
            continue;
        };
        if !end.is_same(&next) && !model.same_position(&end, &next, tol)? {
            found.note(
                Severity::Broken,
                wire,
                ShapeType::Wire,
                format!(
                    "edge {i} ends where edge {} does not begin; a face built on \
                     this has a gap in its boundary",
                    (i + 1) % edges.len()
                ),
            );
        }
    }
    Ok(())
}

/// Every edge of a face needs a pcurve on that face's surface.
fn check_face(
    model: &Model,
    face: &Shape,
    _tol: Tolerances,
    found: &mut Diagnosis,
) -> OgResult<()> {
    let Some(node) = model.node(face) else {
        og_bail!(Dangling, "face is not in this model");
    };
    let Some(data) = node.data().as_face() else {
        return Ok(());
    };
    if model.geometry().surface(data.surface).is_none() {
        og_bail!(Dangling, "face names a surface not in this model");
    }

    let wires = model.children_of(face)?;
    if wires.is_empty() && !data.natural_restriction {
        found.note(
            Severity::Broken,
            face,
            ShapeType::Face,
            "has no wires and is not marked as covering its whole surface, so \
             what it is a face *of* is undefined"
                .into(),
        );
    }

    for wire in &wires {
        for edge in model.children_of(wire)? {
            let Some(edge_data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                continue;
            };
            if edge_data
                .pcurve_for(data.surface, edge.location())
                .is_none()
            {
                found.note(
                    Severity::Broken,
                    &edge,
                    ShapeType::Edge,
                    "bounds a face it has no pcurve on; the face cannot be split \
                     or triangulated in its own parameter space"
                        .into(),
                );
            }
        }
    }
    Ok(())
}

/// A shell is closed when every edge is used an even number of times.
///
/// Reported as `Suspect` rather than `Broken` on its own: an open shell is a
/// perfectly good surface and plenty of operations want one. It becomes
/// `Broken` only when something asks it to bound a volume, which is a question
/// this function is not being asked.
fn check_shell(model: &Model, shell: &Shape, found: &mut Diagnosis) -> OgResult<()> {
    let mut uses: HashMap<TShapeId, usize> = HashMap::new();
    for face in explore(model, shell, Filter::OfType(ShapeType::Face))? {
        for wire in model.children_of(&face)? {
            for edge in model.children_of(&wire)? {
                if model
                    .node(&edge)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(|d| d.degenerate)
                {
                    continue;
                }
                *uses.entry(edge.node()).or_default() += 1;
            }
        }
    }
    let odd = uses.values().filter(|n| *n % 2 == 1).count();
    if odd > 0 {
        found.note(
            Severity::Suspect,
            shell,
            ShapeType::Shell,
            format!(
                "{odd} edge(s) used an odd number of times, so the shell is open \
                 along them; it encloses no volume"
            ),
        );
    }
    Ok(())
}

/// Tolerance containment: a face is no looser than its edges, an edge no looser
/// than its vertices.
///
/// The rule is transitive and the check has to be too. Checking one level would
/// pass a face whose edge is fine and whose *vertex* is tighter than the face —
/// and the containment claim is about the face reaching the vertex.
fn check_containment(model: &Model, shape: &Shape, found: &mut Diagnosis) -> OgResult<()> {
    for face in explore_unique(model, shape, ShapeType::Face)? {
        compare(model, &face, ShapeType::Face, found)?;
    }
    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        compare(model, &edge, ShapeType::Edge, found)?;
    }
    Ok(())
}

/// Compare one shape's tolerance against everything below it.
fn compare(model: &Model, shape: &Shape, kind: ShapeType, found: &mut Diagnosis) -> OgResult<()> {
    let Some(bounding) = model.tolerance_of(shape)? else {
        return Ok(());
    };
    for below in explore(model, shape, Filter::All)? {
        if below.is_same(shape) {
            continue;
        }
        let Some(bounded) = model.tolerance_of(&below)? else {
            continue;
        };
        if bounded.get() < bounding.get() {
            found.note(
                Severity::Broken,
                &below,
                model.kind_of(&below)?,
                format!(
                    "tolerance {} is tighter than the {kind:?} that bounds it \
                     ({}); the bound does not reliably contain what it bounds",
                    bounded.get(),
                    bounding.get()
                ),
            );
        }
    }
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::{make_box, make_cylinder, make_sphere, make_torus};
    use og_core::Tolerance;
    use og_math::{Frame, Point};
    use og_topo::NodeData;
    use og_topo::VertexData;

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn every_primitive_is_valid() {
        // The check earns its keep only if the things known to be right pass
        // it. A checker that flags a correct box is worse than none, because
        // every real finding then reads as noise.
        let mut model = Model::new();
        let shapes = [
            make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T)
                .unwrap()
                .shape,
            make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T)
                .unwrap()
                .shape,
            make_sphere(&mut model, Frame::WORLD, 3.0, T).unwrap().shape,
            make_torus(&mut model, Frame::WORLD, 5.0, 2.0, T)
                .unwrap()
                .shape,
        ];
        for shape in &shapes {
            let found = check(&model, shape, T).unwrap();
            assert!(found.is_valid(), "a primitive was flagged: {found}");
        }
    }

    #[test]
    fn a_prism_is_valid() {
        use og_math::Vector;
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let face = explore_unique(&model, &solid, ShapeType::Face).unwrap()[0].clone();
        let prism = crate::make_prism(&mut model, &face, Vector::new(0.0, 0.0, 2.0), T)
            .unwrap()
            .shape;

        let found = check(&model, &prism, T).unwrap();
        assert!(found.is_valid(), "the prism was flagged: {found}");
    }

    #[test]
    fn a_single_face_is_reported_open_but_still_usable() {
        // An open shell is a perfectly good surface, and plenty of operations
        // want one. Calling it broken would make the checker useless for every
        // sheet body.
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let face = explore_unique(&model, &solid, ShapeType::Face).unwrap()[0].clone();
        let shell = crate::build::make_shell(&mut model, std::slice::from_ref(&face))
            .unwrap()
            .shape;

        let found = check(&model, &shell, T).unwrap();
        assert!(!found.is_valid(), "an open shell is worth reporting");
        assert!(found.is_usable(), "but nothing here answers wrongly");
        assert_eq!(found.worst(), Some(Severity::Suspect));
        assert_eq!(found.of(Severity::Suspect).len(), 1);
    }

    #[test]
    fn a_vertex_tighter_than_its_edge_is_caught() {
        // The containment rule runs the other way from intuition: the *bound*
        // is looser, and a vertex tighter than the edge that caps it means the
        // edge does not reliably reach it.
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let edge = explore_unique(&model, &solid, ShapeType::Edge).unwrap()[0].clone();

        // Widen the edge alone, bypassing the cascading repair that exists to
        // stop exactly this.
        let loose = Tolerance::new(1e-3).unwrap();
        if let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(og_topo::TShape::data_mut) {
            data.tolerance = loose;
        }

        let found = check(&model, &edge, T).unwrap();
        assert!(!found.is_usable(), "a broken containment is not usable");
        let broken = found.of(Severity::Broken);
        assert!(!broken.is_empty());
        assert!(broken.iter().all(|p| p.kind == ShapeType::Vertex));
    }

    #[test]
    fn an_edge_with_no_curve_and_no_excuse_is_caught() {
        let mut model = Model::new();
        let v = model.add_vertex(VertexData::new(Point::ORIGIN));
        let edge = model
            .add_edge(og_topo::EdgeData::new(), &[v.clone(), v])
            .unwrap();

        let found = check(&model, &edge, T).unwrap();
        assert!(!found.is_usable());
        assert!(found.problems[0].what.contains("not marked degenerate"));
    }

    #[test]
    fn an_edge_whose_curve_misses_its_vertex_is_caught() {
        // The failure that leaves a gap in every face built on the wire, and
        // which nothing notices until something walks the boundary.
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let vertex = explore_unique(&model, &solid, ShapeType::Vertex).unwrap()[0].clone();
        if let Some(NodeData::Vertex(data)) = model.node_mut(&vertex).map(og_topo::TShape::data_mut)
        {
            data.point = Point::new(50.0, 50.0, 50.0);
        }

        let found = check(&model, &solid, T).unwrap();
        assert!(!found.is_usable());
        assert!(
            found
                .of(Severity::Broken)
                .iter()
                .any(|p| p.what.contains("from the vertex it should meet")),
            "got {found}"
        );
    }

    #[test]
    fn a_face_whose_edge_has_no_pcurve_is_caught() {
        // Without a pcurve the face cannot be split in a boolean or
        // triangulated at all, and the failure surfaces far from its cause.
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let face = explore_unique(&model, &solid, ShapeType::Face).unwrap()[0].clone();
        let edge = model
            .children_of(&model.children_of(&face).unwrap()[0])
            .unwrap()[0]
            .clone();
        if let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(og_topo::TShape::data_mut) {
            data.representations.retain(|r| r.is_curve3d());
        }

        let found = check(&model, &face, T).unwrap();
        assert!(!found.is_usable());
        assert!(
            found
                .of(Severity::Broken)
                .iter()
                .any(|p| p.what.contains("no pcurve on")),
            "got {found}"
        );
    }

    #[test]
    fn a_diagnosis_reads_as_a_report_rather_than_a_debug_dump() {
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        assert_eq!(check(&model, &solid, T).unwrap().to_string(), "valid");

        let face = explore_unique(&model, &solid, ShapeType::Face).unwrap()[0].clone();
        let shell = crate::build::make_shell(&mut model, std::slice::from_ref(&face))
            .unwrap()
            .shape;
        let text = check(&model, &shell, T).unwrap().to_string();
        assert!(text.starts_with("[suspect] Shell:"), "got {text}");
        assert!(text.contains("open"), "got {text}");
    }

    #[test]
    fn a_handle_that_does_not_resolve_is_an_error_not_a_finding() {
        // A dangling handle means the shape and the model do not belong
        // together. Every finding would then be about something that is not
        // there, which is worse than no finding.
        //
        // Note what this does *not* catch: a handle from a different model
        // whose arena happens to have a node at the same index resolves
        // silently, because arena keys are not scoped to a model. Detecting
        // that needs an identifier on the model itself — see the deferred list
        // in `docs/SCOPE.md`.
        let mut other = Model::new();
        for _ in 0..4 {
            other.add_vertex(VertexData::new(Point::ORIGIN));
        }
        let beyond = other.add_vertex(VertexData::new(Point::ORIGIN));

        let empty = Model::new();
        assert!(check(&empty, &beyond, T).is_err());
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tessellation_tests {
    use super::*;
    use crate::{make_box, make_cone, make_cylinder, make_sphere, make_torus};
    use og_math::{Frame, Point};
    use og_topo::{NodeData, VertexData};

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> Deflection {
        Deflection {
            chord: 0.02,
            ..Deflection::default()
        }
    }

    #[test]
    fn every_primitive_tessellates_into_a_mesh_that_agrees_with_its_topology() {
        // The regression net for every seam and every pole. A primitive whose
        // shell closes but whose mesh does not is the failure this exists to
        // name, and it is invisible to every other check.
        let mut model = Model::new();
        let shapes = [
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
        ];
        for shape in &shapes {
            let found = check_tessellation(&model, shape, fine(), T).unwrap();
            assert!(found.is_valid(), "a primitive's mesh came apart: {found}");
        }
    }

    #[test]
    fn a_prism_tessellates_into_an_agreeing_mesh_whichever_face_it_swept() {
        // This check found the defect that made the distinction matter:
        // sweeping the *downward* face of a box produced four lateral faces
        // that every one of them failed to triangulate, while the shell still
        // closed and every other check passed. Both directions are covered
        // here now, so a regression cannot hide behind the one that worked.
        use og_math::Vector;
        for role in [
            crate::primitive::roles::FACE_MAX_Z,
            crate::primitive::roles::FACE_MIN_Z,
        ] {
            let mut model = Model::new();
            let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
                .unwrap()
                .shape;
            let face = explore_unique(&model, &solid, ShapeType::Face)
                .unwrap()
                .into_iter()
                .find(|f| model.provenance_of(f).and_then(og_core::Provenance::role) == Some(role))
                .expect("the box has a face with that role");
            let prism = crate::make_prism(&mut model, &face, Vector::new(0.0, 0.0, 2.0), T)
                .unwrap()
                .shape;
            assert!(
                check_tessellation(&model, &prism, fine(), T)
                    .unwrap()
                    .is_valid(),
                "{role:?}"
            );
            assert!(check(&model, &prism, T).unwrap().is_valid(), "{role:?}");
        }
    }

    #[test]
    fn moving_a_vertex_does_not_move_the_mesh() {
        // Worth pinning, because it is unintuitive and it invalidated an
        // earlier attempt at a test here. Tessellation reads curves and
        // pcurves, never vertex positions — so a vertex moved off its edges is
        // caught by `check` (the curve no longer reaches it) and is invisible
        // to `check_tessellation`. The two checks genuinely see different
        // things, which is why both exist.
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T)
            .unwrap()
            .shape;
        let before = og_mesh::triangulate(&model, &solid, fine(), T).unwrap();

        let vertex = explore_unique(&model, &solid, ShapeType::Vertex).unwrap()[0].clone();
        if let Some(NodeData::Vertex(data)) = model.node_mut(&vertex).map(og_topo::TShape::data_mut)
        {
            data.point = Point::new(0.5, 0.5, 0.5);
        }

        let after = og_mesh::triangulate(&model, &solid, fine(), T).unwrap();
        assert_eq!(before.positions, after.positions);
        assert!(
            check_tessellation(&model, &solid, fine(), T)
                .unwrap()
                .is_valid()
        );
        assert!(
            !check(&model, &solid, T).unwrap().is_usable(),
            "check sees it"
        );
    }

    #[test]
    fn an_open_shell_is_not_reported_because_it_never_claimed_to_close() {
        // A mesh with a boundary is agreement here, not disagreement. Flagging
        // it would make the check useless for every sheet body.
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let face = explore_unique(&model, &solid, ShapeType::Face).unwrap()[0].clone();
        let shell = crate::build::make_shell(&mut model, std::slice::from_ref(&face))
            .unwrap()
            .shape;
        assert!(
            check_tessellation(&model, &shell, fine(), T)
                .unwrap()
                .is_valid()
        );
    }

    #[test]
    fn a_shape_with_no_shell_has_nothing_to_disagree_about() {
        let mut model = Model::new();
        let vertex = model.add_vertex(VertexData::new(Point::ORIGIN));
        assert!(
            check_tessellation(&model, &vertex, fine(), T)
                .unwrap()
                .is_valid()
        );
    }
}
