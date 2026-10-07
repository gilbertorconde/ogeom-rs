//! Checking a shape against the model's invariants.
//!
//! Everything in `docs/DATA_MODEL.md` that can be checked with geometry in
//! hand, checked in one place. The builders enforce what they can at the moment
//! of construction; this catches what only becomes wrong later: an edge whose
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
//! answer rather than an error: an open shell has no inside, so every
//! containment test against it is a coin toss. [`Severity::Suspect`] means
//! something is out of order but every operation will still behave: a tolerance
//! larger than the feature it describes is alarming and not yet wrong.
//!
//! The distinction is what lets a pipeline decide. Booleans refuse `Broken`
//! input because they would produce nonsense from it; they proceed on
//! `Suspect` because refusing would reject most real imported geometry.

use ogeom_core::{FastMap, FastSet};
use std::fmt;

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve2d, Curve3d, Surface};
use ogeom_math::Point2;
use ogeom_mesh::Deflection;
use ogeom_topo::{EdgeRepr, Filter, Model, Shape, ShapeType, TShapeId, explore, explore_unique};

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
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle does not
/// resolve. A dangling handle is not a *finding*; it means the shape and the
/// model do not belong together, and every other answer would be about
/// something that is not there.
pub fn check(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Diagnosis> {
    if model.node(shape).is_none() {
        ogeom_bail!(Dangling, "shape refers to a node not in this model");
    }
    let mut found = Diagnosis::default();

    // Every distinct sub-shape from one walk, by type. The walk is
    // pre-order, so each type's shapes come in the order a walk for that
    // type alone would give.
    let mut distinct: FastMap<ShapeType, Vec<Shape>> = FastMap::default();
    let mut seen = FastSet::default();
    for sub in explore(model, shape, Filter::All)? {
        if seen.insert(ogeom_topo::SameKey(sub.clone())) {
            distinct.entry(model.kind_of(&sub)?).or_default().push(sub);
        }
    }
    let of = |kind: ShapeType| distinct.get(&kind).map_or(&[][..], Vec::as_slice);

    for edge in of(ShapeType::Edge) {
        check_edge(model, edge, tol, &mut found)?;
    }
    for wire in of(ShapeType::Wire) {
        check_wire(model, wire, tol, &mut found)?;
    }
    for face in of(ShapeType::Face) {
        check_face(model, face, tol, &mut found)?;
        check_face_chart(model, face, tol, &mut found)?;
    }
    for shell in of(ShapeType::Shell) {
        check_shell(model, shell, &mut found)?;
    }
    for solid in of(ShapeType::Solid) {
        check_orientation(model, solid, tol, &mut found)?;
    }
    // Tolerance containment: a face is no looser than its edges, an edge no
    // looser than its vertices, checked through every level below each.
    for face in of(ShapeType::Face) {
        compare(model, face, ShapeType::Face, &mut found)?;
    }
    for edge in of(ShapeType::Edge) {
        compare(model, edge, ShapeType::Edge, &mut found)?;
    }
    Ok(found)
}

/// Report every face of a solid that faces into its material.
fn check_orientation(
    model: &Model,
    solid: &Shape,
    tol: Tolerances,
    found: &mut Diagnosis,
) -> OgeomResult<()> {
    for face in inside_out_faces(model, solid, tol)? {
        found.note(
            Severity::Broken,
            &face,
            ShapeType::Face,
            "the face's outward normal points into the material: its \
             orientation is reversed"
                .into(),
        );
    }
    Ok(())
}

/// The faces of a solid that face into its material.
///
/// A face's flag is the only thing that says which side of its surface the
/// material is on, and a flipped one leaves every edge used twice and the
/// shell closed: nothing topological sees it. Nor do the flags asked about
/// each other along shared edges: that walk reads which side of an edge a
/// face lies on from its chart, and where the chart cannot say (an apex, a
/// boundary whose middle is not inside it) it has to take the flags as they
/// come. So every face is probed: a point a little off it along its outward
/// normal must lie outside the material and the point as far the other way
/// inside it. A face whose outside is inside is named; a probe that cannot
/// tell names nothing, so a face is never named on a guess. A face whose
/// chart cannot be walked is not asked at all.
///
/// # Errors
///
/// [`OgeomError::Cancelled`](ogeom_core::OgeomError::Cancelled) and
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) only; a
/// question the geometry cannot answer names no face.
pub fn inside_out_faces(model: &Model, solid: &Shape, tol: Tolerances) -> OgeomResult<Vec<Shape>> {
    let boundary = match probe_boundary(model, solid, tol) {
        Ok(boundary) => boundary,
        // A question the geometry cannot answer (an open shell, a chart
        // walk off its domain) is no finding; only cancellation and a
        // broken model are errors.
        Err(e @ (ogeom_core::OgeomError::Cancelled | ogeom_core::OgeomError::Dangling(_))) => {
            return Err(e);
        }
        Err(_) => return Ok(Vec::new()),
    };
    let faces = explore_unique(model, solid, ShapeType::Face)?;
    turned_in_against(model, &boundary, &faces, tol)
}

/// The boundary of `solid` its faces are probed against, as
/// [`inside_out_faces`] builds it.
pub(crate) fn probe_boundary(
    model: &Model,
    solid: &Shape,
    tol: Tolerances,
) -> OgeomResult<crate::SolidBoundary> {
    crate::SolidBoundary::of(model, solid, tol.confusion() * 1e4, tol)
}

/// The faces of `faces` whose probe against `boundary` finds material on
/// their outside, probed in parallel.
fn turned_in_against(
    model: &Model,
    boundary: &crate::SolidBoundary,
    faces: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let inward = ogeom_core::parallel::map_ordered(faces, |_, face| {
        faces_inward(model, face, boundary, tol)
    });
    let mut out = Vec::new();
    for (face, inward) in faces.iter().zip(inward) {
        if inward? == Some(true) {
            out.push(face.clone());
        }
    }
    Ok(out)
}

/// Whether a face's outward normal points into the solid's material; `None`
/// where the probe cannot tell.
///
/// Probed at the middle of the face's largest mesh triangle, or of the next
/// largest where that cannot tell. Each probe steps off the face both ways,
/// first by a fraction of the triangle's size and then by shorter steps: a
/// step longer than the solid is thick lands past its far side, and one
/// across a gap narrower than the step lands in the material beyond, and
/// either leaves both points on the same side. The shortest step stays well
/// clear of the band within which a point counts as on the face.
pub(crate) fn faces_inward(
    model: &Model,
    face: &Shape,
    boundary: &crate::SolidBoundary,
    tol: Tolerances,
) -> OgeomResult<Option<bool>> {
    /// How many of the largest triangles are asked before giving up.
    const PROBES: usize = 3;
    let Ok(mesh) = ogeom_mesh::triangulate_face(model, face, Deflection::default(), tol) else {
        return Ok(None);
    };
    if mesh.parameters.len() != mesh.positions.len() {
        return Ok(None);
    }
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let area = |t: &[u32; 3]| {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        (b - a).cross(c - a).magnitude()
    };
    let mut largest: Vec<(f64, &[u32; 3])> = mesh.triangles.iter().map(|t| (area(t), t)).collect();
    largest.sort_by(|x, y| y.0.total_cmp(&x.0));
    let placement = face.transform(model.datums())?;
    let flip = face.orientation() == ogeom_topo::Orientation::Reversed;
    // Twenty times the band a point counts as on the face within, and
    // never below a thousand confusion distances.
    let shortest = data.tolerance.get().max(tol.confusion() * 1e3) * 20.0;
    for (size, triangle) in largest.into_iter().take(PROBES) {
        let uv = triangle.map(|i| mesh.parameters[i as usize]);
        let (u, v) = (
            (uv[0].0 + uv[1].0 + uv[2].0) / 3.0,
            (uv[0].1 + uv[1].1 + uv[2].1) / 3.0,
        );
        let (Ok(point), Ok(normal)) = (surface.point_at(u, v, tol), surface.normal_at(u, v, tol))
        else {
            continue;
        };
        let point = placement.apply(point);
        let mut outward = placement.apply_vector(normal.vector());
        if flip {
            outward = -outward;
        }
        let mut step = (size.sqrt() * 0.05).max(tol.confusion() * 1e5);
        loop {
            ogeom_core::progress::checkpoint()?;
            if let (Ok(outside), Ok(inside)) = (
                boundary.holds(model, point + outward * step, tol),
                boundary.holds(model, point - outward * step, tol),
            ) {
                match (outside, inside) {
                    (crate::Containment::In, crate::Containment::Out) => return Ok(Some(true)),
                    (crate::Containment::Out, crate::Containment::In) => return Ok(Some(false)),
                    _ => {}
                }
            }
            if step <= shortest {
                break;
            }
            step = (step * 0.1).max(shortest);
        }
    }
    Ok(None)
}

/// Check that a shape's *tessellation* agrees with its topology.
///
/// Separate from [`check`] because it needs a deflection, and because it asks a
/// different question: not "is this shape well formed" but "do its two
/// descriptions of itself agree". A shell whose edges are all used twice is
/// closed as far as the topology knows. If the mesh built from it still has a
/// boundary, then some face's pcurves do not cover the region its edges claim
/// to bound, and the topology cannot see that, because the defect is entirely
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
/// As [`ogeom_mesh::triangulate()`].
pub fn check_tessellation(
    model: &Model,
    shape: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Diagnosis> {
    let mut found = Diagnosis::default();

    for shell in explore_unique(model, shape, ShapeType::Shell)? {
        // An open shell is *meant* to have a boundary, so a mesh with one is
        // agreement, not disagreement. Only a shell the topology calls closed
        // makes a claim the mesh can contradict.
        if !crate::build::is_shell_closed(model, &shell)? {
            continue;
        }
        let mesh = ogeom_mesh::triangulate(model, &shell, deflection, tol)?;
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
fn open_edges(mesh: &ogeom_topo::Triangulation) -> Option<String> {
    let mut uses: FastMap<(u32, u32), usize> = FastMap::default();
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
         slit in it. The first are at {}. This is a parameter-space defect \
         (some face's pcurves do not cover the region its edges bound), and no \
         topological check can see it",
        loose.len(),
        sample.join(", ")
    ))
}

/// An edge's curve must reach the vertices it claims to join, and its
/// representations must agree with each other.
fn check_edge(
    model: &Model,
    edge: &Shape,
    tol: Tolerances,
    found: &mut Diagnosis,
) -> OgeomResult<()> {
    let Some(node) = model.node(edge) else {
        ogeom_bail!(Dangling, "edge is not in this model");
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
        ogeom_bail!(Dangling, "curve is not in this model");
    };

    let placement = edge.transform(model.datums())?;
    let bounds = model.children_of(edge)?;
    for (parameter, vertex) in [(range.0, bounds.first()), (range.1, bounds.last())] {
        let Some(vertex) = vertex else { continue };
        let Some((point, vertex_reach)) = model
            .node(vertex)
            .and_then(|n| n.data().as_vertex())
            .map(|v| (v.point, v.tolerance.get()))
        else {
            continue;
        };
        let placed = vertex.transform(model.datums())?.apply(point);
        let on_curve = placement.apply(geometry.point_at(parameter, tol)?);
        let gap = on_curve.distance(placed);
        // The junction's own stated tolerance is the radius within which
        // things meeting it may stray; the same acceptance construction
        // applies. A checker stricter than the builder would condemn what
        // the builder rightly admitted and honestly recorded.
        let reach = reach.max(vertex_reach);
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

    let curve_location = match data.curve3d() {
        Some(EdgeRepr::Curve3d { location, .. }) => location.clone(),
        _ => ogeom_topo::Location::identity(),
    };
    check_pcurves(
        model,
        edge,
        data,
        (geometry, *range, &curve_location),
        reach,
        tol,
        found,
    )?;
    Ok(())
}

/// Verify that every pcurve lands where the 3D curve does, to the edge's
/// tolerance.
///
/// The tolerance is the radius about the curve within which every
/// description of the edge lies, so each pcurve lifted through its surface
/// stays inside it all along the edge, not only at its ends. An edge
/// claiming `same_parameter` is held to the claim: the lifted point at each
/// parameter within the tolerance of the curve's point at the same
/// parameter, since every algorithm evaluates whichever representation is
/// cheapest and assumes the answer interchangeable. An edge making no claim
/// may pace its pcurves differently, so a lifted point is measured against
/// the nearest point of the curve's stretch. A pcurve at another placement
/// than the curve describes the edge where that occurrence stands, and is
/// not compared.
fn check_pcurves(
    model: &Model,
    edge: &Shape,
    data: &ogeom_topo::EdgeData,
    (curve, range, location): (&ogeom_geom::Curve, (f64, f64), &ogeom_topo::Location),
    reach: f64,
    tol: Tolerances,
    found: &mut Diagnosis,
) -> OgeomResult<()> {
    check_plain_pcurves_apart(model, edge, data, tol, found)?;
    let claimed = data.same_parameter();
    for repr in &data.representations {
        let (sides, pcurve_range, surface_id, at) = match repr {
            EdgeRepr::PCurve {
                curve,
                range,
                surface,
                location,
            } => ([Some(*curve), None], *range, *surface, location),
            EdgeRepr::Seam {
                forward,
                reversed,
                range,
                surface,
                location,
            } => (
                [Some(*forward), Some(*reversed)],
                *range,
                *surface,
                location,
            ),
            _ => continue,
        };
        if at != location {
            continue;
        }
        let Some(surface) = model.geometry().surface(surface_id) else {
            ogeom_bail!(Dangling, "an edge names geometry not in this model");
        };
        for pcurve_id in sides.into_iter().flatten() {
            let Some(pcurve) = model.geometry().pcurve(pcurve_id) else {
                ogeom_bail!(Dangling, "an edge names geometry not in this model");
            };
            let Some((gap, t)) = pcurve_off_curve(
                (curve, range),
                (pcurve, pcurve_range),
                surface,
                reach,
                claimed,
                tol,
            )?
            else {
                continue;
            };
            let what = if claimed {
                format!(
                    "claims same_parameter but its pcurve is {gap} from its curve at \
                     parameter {t} of the range, outside the edge's tolerance of {reach}"
                )
            } else {
                format!(
                    "has a pcurve that leaves its curve by {gap} at parameter {t} of \
                     the range, outside the edge's tolerance of {reach}; the face it \
                     bounds and the curve its neighbour follows part there"
                )
            };
            found.note(Severity::Broken, edge, ShapeType::Edge, what);
        }
    }
    Ok(())
}

/// How far a pcurve, lifted through its surface, leaves its edge's curve
/// beyond `reach`: the widest gap and where along the range it stands, or
/// `None` where it keeps within `reach` everywhere sampled. Held to the
/// same parameter where `paced`; otherwise a lifted point is measured
/// against the nearest point of the curve's stretch.
pub(crate) fn pcurve_off_curve(
    (curve, range): (&ogeom_geom::Curve, (f64, f64)),
    (pcurve, pcurve_range): (&ogeom_geom::PlanarCurve, (f64, f64)),
    surface: &ogeom_geom::SurfaceGeometry,
    reach: f64,
    paced: bool,
    tol: Tolerances,
) -> OgeomResult<Option<(f64, f64)>> {
    // A pcurve fitted through its edge's points strays most between them,
    // so the samples are many more than a fit takes per span.
    const SAMPLES: u32 = 32;
    let mut widest: Option<(f64, f64)> = None;
    for i in 0..=SAMPLES {
        let t = f64::from(i) / f64::from(SAMPLES);
        let on_curve = curve.point_at(range.0 + (range.1 - range.0) * t, tol)?;
        let uv = pcurve.point_at(pcurve_range.0 + (pcurve_range.1 - pcurve_range.0) * t, tol)?;
        let Ok(lifted) = surface.point_at(uv.x, uv.y, tol) else {
            continue;
        };
        let mut gap = on_curve.distance(lifted);
        if gap > reach && !paced {
            gap = gap.min(crate::pcurve_gap::nearest_on_stretch(
                curve, range, lifted, tol,
            )?);
        }
        if gap > reach && widest.is_none_or(|(g, _)| gap > g) {
            widest = Some((gap, t));
        }
    }
    Ok(widest)
}

/// A wire's edges must meet end to end.
fn check_wire(
    model: &Model,
    wire: &Shape,
    tol: Tolerances,
    found: &mut Diagnosis,
) -> OgeomResult<()> {
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
        if !end.is_same(&next)
            && !model.same_position(&end, &next, tol)?
            && !crate::build::one_point(model, &end, &next, tol)?
        {
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
) -> OgeomResult<()> {
    let Some(node) = model.node(face) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let Some(data) = node.data().as_face() else {
        return Ok(());
    };
    if model.geometry().surface(data.surface).is_none() {
        ogeom_bail!(Dangling, "face names a surface not in this model");
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

/// An edge's plain pcurves on one surface, at one placement, must agree.
///
/// A face finds its edge's pcurve by surface and placement, and the lookup
/// answers with the first that matches. Where two plain pcurves stand on
/// one surface and part (the two columns of a closed surface's seam,
/// attached as two plain pcurves for two faces on that surface), every face
/// reads the first, and a face that walks the other column reads its
/// boundary a period away. Apart is measured along the surface between
/// them: two fits of one edge stay within its tolerance of each other,
/// two columns of a seam stand the surface's girth apart.
fn check_plain_pcurves_apart(
    model: &Model,
    edge: &Shape,
    data: &ogeom_topo::EdgeData,
    tol: Tolerances,
    found: &mut Diagnosis,
) -> OgeomResult<()> {
    let reach = CHART_GAP_SLACK * data.tolerance.get().max(tol.approximation());
    let plain: Vec<_> = data
        .representations
        .iter()
        .filter_map(|r| match r {
            EdgeRepr::PCurve {
                curve,
                surface,
                location,
                range,
            } => Some((*curve, *surface, location, *range)),
            _ => None,
        })
        .collect();
    for (i, a) in plain.iter().enumerate() {
        for b in &plain[i + 1..] {
            if a.1 != b.1 || a.2 != b.2 {
                continue;
            }
            let (Some(surface), Some(pa), Some(pb)) = (
                model.geometry().surface(a.1),
                model.geometry().pcurve(a.0),
                model.geometry().pcurve(b.0),
            ) else {
                ogeom_bail!(Dangling, "an edge names geometry not in this model");
            };
            let apart = pcurves_apart(surface, (pa, a.3), (pb, b.3), tol)?;
            if apart > reach {
                found.note(
                    Severity::Broken,
                    edge,
                    ShapeType::Edge,
                    format!(
                        "holds two plain pcurves on one surface that stand {apart} \
                         apart along it; a lookup by surface finds the first for \
                         every face, so a face walking the other reads its boundary \
                         in the wrong place"
                    ),
                );
                return Ok(());
            }
        }
    }
    Ok(())
}

/// The widest distance between two pcurves, each sampled evenly over its
/// own range, measured along the surface between them.
fn pcurves_apart(
    surface: &ogeom_geom::SurfaceGeometry,
    (a, ra): (&ogeom_geom::PlanarCurve, (f64, f64)),
    (b, rb): (&ogeom_geom::PlanarCurve, (f64, f64)),
    tol: Tolerances,
) -> OgeomResult<f64> {
    const SAMPLES: u32 = 4;
    let mut widest = 0.0_f64;
    for i in 0..=SAMPLES {
        let t = f64::from(i) / f64::from(SAMPLES);
        let pa = a.point_at(ra.0 + (ra.1 - ra.0) * t, tol)?;
        let pb = b.point_at(rb.0 + (rb.1 - rb.0) * t, tol)?;
        if let Some(apart) = along_surface(surface, pa, pb, tol) {
            widest = widest.max(apart);
        }
    }
    Ok(widest)
}

/// A face's boundary, read edge by edge through the pcurve lookup the
/// boolean and the mesher use, must close in its surface's chart.
///
/// Each edge is read as its face walks it: a plain pcurve from its start to
/// its end, or reversed; a seam on the side whose start continues the point
/// walked to, as the mesher picks it. At each joint, one edge's end and the
/// next edge's start stand within the joint's tolerance of each other,
/// measured along the surface between them, or a whole period apart along
/// a direction the surface closes in. Otherwise the boundary does not close
/// there.
///
/// A ring that goes once round the surface with no seam (a band's rim, one
/// circle) closes by a single period's jump, and each such jump counts
/// toward the ring's turns round the surface. Jumps that cancel (a period
/// one way at one joint, back at another) at an edge that holds another
/// plain pcurve on the surface close nothing: the lookup found the first,
/// and the face walks the other, on the far column of the surface's seam.
/// Cancelling jumps at edges with one pcurve each are a ring folding across
/// the chart's join, which the mesher folds back, and are not reported.
///
/// Either way a boolean through the face finds a strand that dangles. A
/// wire with an edge that has no pcurve on the face is reported by
/// [`check_face`] and not walked.
fn check_face_chart(
    model: &Model,
    face: &Shape,
    tol: Tolerances,
    found: &mut Diagnosis,
) -> OgeomResult<()> {
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
        return Ok(());
    };
    let surface_id = data.surface;
    let Some(surface) = model.geometry().surface(surface_id) else {
        ogeom_bail!(Dangling, "face names a surface not in this model");
    };
    let ((u0, u1), (v0, v1)) = surface.domain();
    let periods = [
        (surface.is_periodic_u() || surface.is_closed_u(tol)).then_some(u1 - u0),
        (surface.is_periodic_v() || surface.is_closed_v(tol)).then_some(v1 - v0),
    ];
    let stored = face.oriented(ogeom_topo::Orientation::Forward);
    'wires: for wire in model.ordered_children_of(&stored)? {
        let edges = model.ordered_children_of(&wire)?;
        let Some(ends) = chart_ends(model, surface_id, &edges, tol)? else {
            continue;
        };
        // The walk starts off a seam: a seam's side is picked by the point
        // already walked to.
        let Some(start) = ends.iter().position(|(e, _)| e.len() == 1) else {
            continue;
        };
        let n = edges.len();
        let mut walked: Vec<(Point2, Point2)> = Vec::with_capacity(n);
        for k in 0..n {
            let sides = &ends[(start + k) % n].0;
            let picked = match walked.last() {
                Some(&(_, last)) => sides
                    .iter()
                    .min_by(|a, b| last.distance(a.0).total_cmp(&last.distance(b.0))),
                None => sides.first(),
            };
            let Some(&picked) = picked else { break };
            walked.push(picked);
        }
        // Whole periods jumped at the joints, signed and in all, per
        // direction.
        let mut turns = [0_i64; 2];
        let mut jumps = [0_i64; 2];
        let mut misread = false;
        for k in 0..walked.len() {
            let (i, j) = ((start + k) % n, (start + k + 1) % n);
            let (end, next) = (walked[k].1, walked[(k + 1) % walked.len()].0);
            let joint = crate::build::edge_vertices(model, &edges[i])?
                .and_then(|(_, v)| model.tolerance_of(&v).ok().flatten())
                .map_or(0.0, |t| t.get());
            let reach = CHART_GAP_SLACK
                * [&edges[i], &edges[j]]
                    .iter()
                    .filter_map(|e| model.tolerance_of(e).ok().flatten())
                    .map(|t| t.get())
                    .fold(joint.max(tol.approximation()), f64::max);
            let Some(apart) = along_surface(surface, end, next, tol) else {
                continue;
            };
            if apart <= reach {
                continue;
            }
            // The next edge's start brought to the end's copy of the chart.
            let mut whole = [0_i64; 2];
            let mut moved = next;
            for (axis, period) in periods.iter().enumerate() {
                let Some(period) = period.filter(|p| *p > 0.0 && p.is_finite()) else {
                    continue;
                };
                let (from, to) = if axis == 0 {
                    (end.x, &mut moved.x)
                } else {
                    (end.y, &mut moved.y)
                };
                #[allow(clippy::cast_possible_truncation, reason = "a few periods")]
                let shift = ((from - *to) / period).round() as i64;
                #[allow(clippy::cast_precision_loss, reason = "a few periods")]
                {
                    *to += shift as f64 * period;
                }
                whole[axis] = -shift;
            }
            let closes = whole != [0, 0]
                && along_surface(surface, end, moved, tol).is_some_and(|gap| gap <= reach);
            if !closes {
                found.note(
                    Severity::Broken,
                    face,
                    ShapeType::Face,
                    format!(
                        "its boundary does not close in its surface's chart: edge \
                         {i} ends {apart} along the surface from where edge {j} \
                         begins, outside the joint's tolerance of {reach}; a \
                         boolean through the face finds a strand that dangles"
                    ),
                );
                continue 'wires;
            }
            for axis in 0..2 {
                turns[axis] += whole[axis];
                jumps[axis] += whole[axis].abs();
            }
            misread |= ends[i].1 || ends[j].1;
        }
        if misread && (0..2).any(|axis| jumps[axis] > turns[axis].abs()) {
            found.note(
                Severity::Broken,
                face,
                ShapeType::Face,
                "its boundary does not close in its surface's chart: it jumps a \
                 whole period at joints whose jumps cancel, so an edge between \
                 them is read on the far column of the surface's seam; a boolean \
                 through the face finds a strand that dangles"
                    .into(),
            );
        }
    }
    Ok(())
}

/// How many times its joint's tolerance a gap in a face's chart may span
/// before the boundary is open there. A fitted trim's ends wander from its
/// vertex by about the tolerance its edge states, and never less than the
/// fitting accuracy, which is where the tolerance taken starts.
const CHART_GAP_SLACK: f64 = 10.0;

/// An edge's ends in a face's chart as the face walks it (one pair for a
/// pcurve, one per side for a seam), and whether the edge holds another
/// plain pcurve on the surface at the same placement.
type ChartEnds = (Vec<(Point2, Point2)>, bool);

/// Each edge's [`ChartEnds`]; `None` where an edge has no pcurve on the
/// surface.
fn chart_ends(
    model: &Model,
    surface: ogeom_topo::SurfaceId,
    edges: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Option<Vec<ChartEnds>>> {
    let mut ends = Vec::with_capacity(edges.len());
    for edge in edges {
        let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
            return Ok(None);
        };
        let (ids, range, another) = match data.pcurve_for(surface, edge.location()) {
            Some(EdgeRepr::PCurve {
                curve,
                range,
                location,
                ..
            }) => {
                let plain = data
                    .representations
                    .iter()
                    .filter(|r| {
                        matches!(r, EdgeRepr::PCurve { .. })
                            && r.surface() == Some(surface)
                            && r.location() == Some(location)
                    })
                    .count();
                (vec![*curve], *range, plain > 1)
            }
            Some(EdgeRepr::Seam {
                forward,
                reversed,
                range,
                ..
            }) => (vec![*forward, *reversed], *range, false),
            _ => return Ok(None),
        };
        let mut sides = Vec::with_capacity(ids.len());
        for id in ids {
            let Some(pcurve) = model.geometry().pcurve(id) else {
                ogeom_bail!(Dangling, "an edge names geometry not in this model");
            };
            let (a, b) = (
                pcurve.point_at(range.0, tol)?,
                pcurve.point_at(range.1, tol)?,
            );
            sides.push(if edge.orientation() == ogeom_topo::Orientation::Reversed {
                (b, a)
            } else {
                (a, b)
            });
        }
        ends.push((sides, another));
    }
    Ok(Some(ends))
}

/// The length of the straight chart segment from `a` to `b`, lifted onto
/// the surface: nothing where the two meet at a pole, a period's worth of
/// surface where they stand on two columns of a seam. `None` where the
/// surface cannot be evaluated along it.
fn along_surface(
    surface: &ogeom_geom::SurfaceGeometry,
    a: Point2,
    b: Point2,
    tol: Tolerances,
) -> Option<f64> {
    const PIECES: u32 = 8;
    if a.distance(b) <= tol.parametric() {
        return Some(0.0);
    }
    let at = |t: f64| surface.point_at(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, tol);
    let mut last = at(0.0).ok()?;
    let mut length = 0.0;
    for i in 1..=PIECES {
        let next = at(f64::from(i) / f64::from(PIECES)).ok()?;
        length += last.distance(next);
        last = next;
    }
    Some(length)
}

/// A shell is closed when every edge is used an even number of times.
///
/// Reported as `Suspect` rather than `Broken` on its own: an open shell is a
/// perfectly good surface and plenty of operations want one. It becomes
/// `Broken` only when something asks it to bound a volume, which is a question
/// this function is not being asked.
fn check_shell(model: &Model, shell: &Shape, found: &mut Diagnosis) -> OgeomResult<()> {
    let mut uses: FastMap<TShapeId, usize> = FastMap::default();
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
    for edge in walked_one_way(model, shell)? {
        found.note(
            Severity::Suspect,
            &edge,
            ShapeType::Edge,
            "walked the same way by both its faces: one of them does not keep \
             its material on the left of its rings"
                .into(),
        );
    }
    Ok(())
}

/// The edges between two faces of `shape` that both faces walk the same
/// way, in the order the walk first meets them.
///
/// Each face keeps its material on the left of its rings about its outward
/// normal, so an edge two faces share is walked once each way, senses
/// composed. That holds for any edge with two faces, closed shell or sheet;
/// an edge with one face or more than two says nothing. A degenerate edge
/// (a pole, an apex) bounds nothing and is not asked, and a seam, walked
/// both ways by its one face, passes.
pub(crate) fn walked_one_way(model: &Model, shape: &Shape) -> OgeomResult<Vec<Shape>> {
    type Walks = (Shape, usize, usize);
    let mut walks: FastMap<(TShapeId, ogeom_topo::Location), Walks> = FastMap::default();
    let mut order = Vec::new();
    for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
        for edge in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
            if model
                .node(&edge)
                .and_then(|n| n.data().as_edge())
                .is_some_and(|d| d.degenerate)
            {
                continue;
            }
            let key = (edge.node(), edge.location().clone());
            let walk = walks.entry(key.clone()).or_insert_with(|| {
                order.push(key);
                (edge.clone(), 0, 0)
            });
            walk.1 += 1;
            if edge.orientation() == ogeom_topo::Orientation::Forward {
                walk.2 += 1;
            }
        }
    }
    Ok(order
        .into_iter()
        .filter_map(|key| {
            let (edge, uses, forward) = walks.remove(&key)?;
            (uses == 2 && forward != 1).then_some(edge)
        })
        .collect())
}

/// Restore tolerance containment below `shape`: every edge widened to at
/// least the faces it bounds, every vertex to at least the edges it bounds.
///
/// The rule [`check`] enforces, established the only way the data model
/// allows: by raising what is bounded, never lowering what bounds. Each
/// face and then each edge is widened to its own tolerance through
/// [`Model::widen`], which cascades to everything below it and leaves
/// anything already looser as it is. Returns how many entities grew.
///
/// An operation that widens an edge's tolerance by writing it directly (a
/// reader recording how far a pcurve sits from its curve) leaves the
/// edge's vertices behind; this is the pass that brings them along. What
/// `shape` shares with other shapes is copied first ([`Model::unshare`]),
/// so their tolerances stay as they are.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the shape,
/// or anything below it, does not resolve in this model.
pub fn restore_containment(model: &mut Model, shape: &Shape) -> OgeomResult<usize> {
    model.unshare(shape)?;
    let bounded: Vec<Shape> = explore_unique(model, shape, ShapeType::Edge)?
        .into_iter()
        .chain(explore_unique(model, shape, ShapeType::Vertex)?)
        .collect();
    let before: Vec<f64> = bounded
        .iter()
        .map(|s| model.tolerance_of(s).map(|t| t.map_or(0.0, |t| t.get())))
        .collect::<OgeomResult<_>>()?;
    for kind in [ShapeType::Face, ShapeType::Edge] {
        for bounding in explore_unique(model, shape, kind)? {
            if let Some(own) = model.tolerance_of(&bounding)? {
                model.widen(&bounding, own)?;
            }
        }
    }
    let mut grown = 0;
    for (s, was) in bounded.iter().zip(before) {
        if model.tolerance_of(s)?.is_some_and(|t| t.get() > was) {
            grown += 1;
        }
    }
    Ok(grown)
}

/// Compare one shape's tolerance against everything below it.
fn compare(
    model: &Model,
    shape: &Shape,
    kind: ShapeType,
    found: &mut Diagnosis,
) -> OgeomResult<()> {
    let Some(bounding) = model.tolerance_of(shape)? else {
        return Ok(());
    };
    // Each occurrence once: a vertex reached through both its edges is one
    // containment claim, not two.
    let mut seen = FastSet::default();
    for below in explore(model, shape, Filter::All)? {
        if below.is_same(shape) || !seen.insert(ogeom_topo::SameKey(below.clone())) {
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

/// Faces of one shape that reach each other without sharing topology:
/// self-intersection, detected as the interference it is.
///
/// Every unordered pair of distinct faces that share no edge and no vertex
/// node is put through the exact minimum-distance machinery; a pair within
/// the confusion tolerance of touching is reported. Adjacent faces meet at
/// their shared boundary by construction and are not interference; a valid
/// solid therefore reports nothing, and a sheet folded through itself names
/// the faces that cross.
///
/// # Errors
///
/// As [`crate::distance_between_shapes`].
pub fn check_self_intersection(
    model: &Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgeomResult<Vec<(Shape, Shape)>> {
    crossings_among(model, shape, None, tol)
}

/// [`check_self_intersection`] over the pairs holding at least one of
/// `near`: after an edit that moved some faces, the pairs of faces it left
/// alone cannot have begun to cross.
///
/// # Errors
///
/// As [`check_self_intersection`].
pub fn check_self_intersection_near(
    model: &Model,
    shape: &Shape,
    near: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Vec<(Shape, Shape)>> {
    let near: ogeom_core::FastSet<TShapeId> = near.iter().map(Shape::node).collect();
    crossings_among(model, shape, Some(&near), tol)
}

fn crossings_among(
    model: &Model,
    shape: &Shape,
    near: Option<&ogeom_core::FastSet<TShapeId>>,
    tol: Tolerances,
) -> OgeomResult<Vec<(Shape, Shape)>> {
    use ogeom_topo::explore_unique;
    let faces = explore_unique(model, shape, ShapeType::Face)?;
    let asked = |i: usize| near.is_none_or(|n| n.contains(&faces[i].node()));
    // The topology below each face, for the adjacency exclusion; each face
    // gathered and bounded once, not once per pair it is in.
    let mut below: Vec<ogeom_core::FastSet<TShapeId>> = Vec::with_capacity(faces.len());
    let mut gathered = Vec::with_capacity(faces.len());
    let mut bounds = Vec::with_capacity(faces.len());
    for face in &faces {
        let mut set = ogeom_core::FastSet::default();
        for kind in [ShapeType::Edge, ShapeType::Vertex] {
            for sub in explore_unique(model, face, kind)? {
                set.insert(sub.node());
            }
        }
        below.push(set);
        gathered.push(crate::proximity::Elements::of(model, face, tol)?);
        bounds.push(crate::measure::shape_bounds(model, face, tol)?.expanded(tol.confusion()));
    }

    let mut crossings = Vec::new();
    for i in 0..faces.len() {
        for j in i + 1..faces.len() {
            ogeom_core::progress::checkpoint()?;
            if !(asked(i) || asked(j))
                || !bounds[i].intersects(&bounds[j])
                || !below[i].is_disjoint(&below[j])
            {
                continue;
            }
            let reach = crate::proximity::distance_between_prepared(
                &gathered[i],
                &gathered[j],
                ogeom_intersect::ExtremaOptions::default(),
                tol,
            )?;
            if reach.distance <= tol.confusion() {
                crossings.push((faces[i].clone(), faces[j].clone()));
            }
        }
    }
    Ok(crossings)
}

/// How many edges between two faces of `shape` both faces walk the same
/// way: none where every face keeps its material on the left of its rings.
#[cfg(test)]
pub(crate) fn edges_walked_one_way(model: &Model, shape: &Shape) -> usize {
    walked_one_way(model, shape).map_or(0, |edges| edges.len())
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::{make_box, make_cylinder, make_sphere, make_torus};
    use ogeom_core::Tolerance;
    use ogeom_math::{Frame, Point};
    use ogeom_topo::NodeData;
    use ogeom_topo::VertexData;

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
        use ogeom_math::Vector;
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
        if let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(ogeom_topo::TShape::data_mut)
        {
            data.tolerance = loose;
        }

        let found = check(&model, &edge, T).unwrap();
        assert!(!found.is_usable(), "a broken containment is not usable");
        let broken = found.of(Severity::Broken);
        assert!(!broken.is_empty());
        assert!(broken.iter().all(|p| p.kind == ShapeType::Vertex));
    }

    #[test]
    fn a_face_too_loose_for_its_boundary_names_each_part_once() {
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let face = explore_unique(&model, &solid, ShapeType::Face).unwrap()[0].clone();
        if let Some(NodeData::Face(data)) = model.node_mut(&face).map(ogeom_topo::TShape::data_mut)
        {
            data.tolerance = Tolerance::new(1e-3).unwrap();
        }
        // Four edges and four corners, each reached through two edges and
        // named once.
        let broken = check(&model, &face, T).unwrap().of(Severity::Broken).len();
        assert_eq!(broken, 8);
    }

    #[test]
    fn an_edge_with_no_curve_and_no_excuse_is_caught() {
        let mut model = Model::new();
        let v = model.add_vertex(VertexData::new(Point::ORIGIN));
        let edge = model
            .add_edge(ogeom_topo::EdgeData::new(), &[v.clone(), v])
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
        if let Some(NodeData::Vertex(data)) =
            model.node_mut(&vertex).map(ogeom_topo::TShape::data_mut)
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
        if let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(ogeom_topo::TShape::data_mut)
        {
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
    fn an_edge_whose_pcurve_leaves_its_curve_is_caught() {
        // A pcurve standing off its edge's curve by more than the edge's
        // tolerance parts the face it bounds from the neighbour following
        // the curve, whether or not the edge claims the two share a
        // parameter. A box edge's pcurve on one face is moved a hundredth
        // sideways in that face's chart.
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let face = explore_unique(&model, &solid, ShapeType::Face).unwrap()[0].clone();
        let face_node = model.node(&face).unwrap().clone();
        let NodeData::Face(face_data) = face_node.data() else {
            panic!("not a face")
        };
        let surface = face_data.surface;
        let edge = model
            .children_of(&model.children_of(&face).unwrap()[0])
            .unwrap()[0]
            .clone();
        assert!(check(&model, &solid, T).unwrap().is_valid());

        let NodeData::Edge(data) = model.node(&edge).unwrap().data().clone() else {
            panic!("not an edge")
        };
        let original = data
            .representations
            .iter()
            .find_map(|r| match r {
                EdgeRepr::PCurve {
                    curve, surface: s, ..
                } if *s == surface => Some(*curve),
                _ => None,
            })
            .expect("the edge has a pcurve on the face");
        let basis = model.geometry().pcurve(original).unwrap().clone();
        let moved = model
            .geometry_mut()
            .add_pcurve(ogeom_geom::PlanarCurve::Offset(Box::new(
                ogeom_geom::Offset2d::new(basis, 0.01).unwrap(),
            )));
        if let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(ogeom_topo::TShape::data_mut)
        {
            for repr in &mut data.representations {
                if let EdgeRepr::PCurve {
                    curve, surface: s, ..
                } = repr
                    && *s == surface
                {
                    *curve = moved;
                }
            }
        }

        let found = check(&model, &solid, T).unwrap();
        assert!(!found.is_usable(), "got {found}");
        assert!(
            found
                .of(Severity::Broken)
                .iter()
                .any(|p| p.kind == ShapeType::Edge && p.what.contains("leaves its curve by")),
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
        // A handle from a *different* model is caught too, and by the same
        // route: arena keys carry the identifier of the arena that issued them,
        // so a foreign one resolves to nothing rather than to whatever sits at
        // that index.
        let mut other = Model::new();
        for _ in 0..4 {
            other.add_vertex(VertexData::new(Point::ORIGIN));
        }
        let beyond = other.add_vertex(VertexData::new(Point::ORIGIN));

        let empty = Model::new();
        assert!(check(&empty, &beyond, T).is_err());
    }

    /// The edges of `diagnosis` named as walked the same way by both their
    /// faces.
    fn walked_one_way_named(diagnosis: &Diagnosis) -> Vec<TShapeId> {
        diagnosis
            .problems
            .iter()
            .filter(|p| {
                p.severity == Severity::Suspect
                    && p.kind == ShapeType::Edge
                    && p.what.starts_with("walked the same way by both its faces")
            })
            .map(|p| p.at.node())
            .collect()
    }

    /// A face of a box rebuilt with its ring walked back, its sense kept,
    /// walks each of its four edges the same way as the face beside it:
    /// `check` names those four edges, in the closed box and in a sheet of
    /// that face and one neighbour, where it names the one edge they share.
    #[test]
    fn an_edge_both_faces_walk_the_same_way_is_named() {
        let mut model = Model::new();
        let block = make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T)
            .unwrap()
            .shape;
        let faces = explore_unique(&model, &block, ShapeType::Face).unwrap();
        let face = faces[2].clone();
        let data = model.node(&face).unwrap().data().as_face().unwrap().clone();
        let stored = face.oriented(ogeom_topo::Orientation::Forward);
        let [wire] = &model.ordered_children_of(&stored).unwrap()[..] else {
            panic!("a box face has one ring");
        };
        let ring = model.ordered_children_of(wire).unwrap();
        let back: Vec<Shape> = ring.iter().rev().map(Shape::reversed).collect();
        let back = crate::build::make_wire(&mut model, &back, T).unwrap().shape;
        let walked_back = model
            .add_face(data, &[back])
            .unwrap()
            .oriented(face.orientation());

        let mut held = faces.clone();
        held[2] = walked_back.clone();
        let shell = crate::build::make_shell(&mut model, &held).unwrap().shape;
        let solid = crate::build::make_solid(&mut model, &[shell])
            .unwrap()
            .shape;
        let diagnosis = check(&model, &solid, T).unwrap();
        let mut named = walked_one_way_named(&diagnosis);
        named.sort();
        let mut ring_nodes: Vec<TShapeId> = ring.iter().map(Shape::node).collect();
        ring_nodes.sort();
        assert_eq!(named, ring_nodes, "{diagnosis}");
        assert!(diagnosis.is_usable(), "{diagnosis}");

        let beside = faces
            .iter()
            .find(|f| {
                let theirs = explore_unique(&model, f, ShapeType::Edge).unwrap();
                !f.is_same(&face) && theirs.iter().any(|e| ring_nodes.contains(&e.node()))
            })
            .unwrap()
            .clone();
        let sheet = crate::build::make_shell(&mut model, &[walked_back, beside])
            .unwrap()
            .shape;
        assert_eq!(
            walked_one_way_named(&check(&model, &sheet, T).unwrap()).len(),
            1
        );
        assert!(walked_one_way_named(&check(&model, &block, T).unwrap()).is_empty());
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tessellation_tests {
    use super::*;
    use crate::{make_box, make_cone, make_cylinder, make_sphere, make_torus};
    use ogeom_math::{Frame, Point};
    use ogeom_topo::{NodeData, VertexData};

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
        // A sweep of a box's downward face must triangulate as well as one of
        // its upward face: lateral walls wound the wrong way fail to
        // triangulate while the shell still closes and every other check
        // passes. Both directions are covered, so a regression cannot hide
        // behind the one that works.
        use ogeom_math::Vector;
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
                .find(|f| {
                    model
                        .provenance_of(f)
                        .and_then(ogeom_core::Provenance::role)
                        == Some(role)
                })
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
        // Tessellation reads curves and pcurves, never vertex positions, so a
        // vertex moved off its edges is caught by `check` (the curve no longer
        // reaches it) and is invisible to `check_tessellation`. The two checks
        // see different things, which is why both exist.
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T)
            .unwrap()
            .shape;
        let before = ogeom_mesh::triangulate(&model, &solid, fine(), T).unwrap();

        let vertex = explore_unique(&model, &solid, ShapeType::Vertex).unwrap()[0].clone();
        if let Some(NodeData::Vertex(data)) =
            model.node_mut(&vertex).map(ogeom_topo::TShape::data_mut)
        {
            data.point = Point::new(0.5, 0.5, 0.5);
        }

        let after = ogeom_mesh::triangulate(&model, &solid, fine(), T).unwrap();
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
