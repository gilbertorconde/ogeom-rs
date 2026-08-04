//! Building topology from geometry.
//!
//! One level above [`Model`]'s raw builders: these take geometry, derive the
//! topology it implies, check the invariants that need geometry to check, and
//! report history.
//!
//! # What "checked" means here
//!
//! [`Model::add_wire`] verifies that its children are edges. That is all it can
//! do, because it has no geometry. [`make_wire`] additionally verifies that
//! consecutive edges actually *meet* — which is the property that makes a wire
//! a connected path rather than a bag of edges, and which every algorithm
//! downstream assumes without checking. A wire whose edges do not join produces
//! a face with a gap in its boundary, and the first thing to notice is usually
//! a boolean, several operations later.

use std::collections::HashMap;

use og_core::{OgResult, Tolerances, og_bail};
use og_geom::{Curve, Curve3d, SurfaceGeometry};
use og_math::Point;
use og_topo::{EdgeData, EdgeRepr, FaceData, Location, Model, Shape, ShapeType, VertexData};

use crate::history::{Built, History};

/// Roles a builder assigns, so a rebuild can match entities up.
pub mod roles {
    use og_core::Role;

    /// The vertex an edge starts at.
    pub const EDGE_START: Role = Role::op_defined(0);
    /// The vertex an edge ends at.
    pub const EDGE_END: Role = Role::op_defined(1);
}

/// Add a vertex at `point`.
pub fn make_vertex(model: &mut Model, point: Point) -> Built {
    Built::from_nothing(model.add_vertex(VertexData::new(point)))
}

/// Build an edge on `curve`, bounded by its own endpoints.
///
/// Creates the two bounding vertices from the curve's ends, so the edge's
/// topology and its geometry cannot disagree about where it starts and stops. A
/// closed curve gets one vertex named twice, which is what keeps "walk to the
/// end" meaningful for a full circle.
///
/// # Errors
///
/// [`OgError::Domain`](og_core::OgError::Domain) if `range` leaves the curve's
/// domain, and [`OgError::Construction`](og_core::OgError::Construction) if it
/// is empty.
pub fn make_edge(
    model: &mut Model,
    curve: Curve,
    range: (f64, f64),
    tol: Tolerances,
) -> OgResult<Built> {
    let (lo, hi) = range;
    if !lo.is_finite() || !hi.is_finite() || hi <= lo + tol.parametric() {
        og_bail!(Construction, "edge range [{lo}, {hi}] is empty");
    }
    let start = curve.point_at(lo, tol)?;
    let end = curve.point_at(hi, tol)?;

    let start_vertex = model.add_vertex(VertexData::new(start));
    // A closed edge names one vertex twice rather than two coincident ones:
    // two would leave the wire looking open at the join.
    let end_vertex = if start.is_equal(end, tol) {
        start_vertex.clone()
    } else {
        model.add_vertex(VertexData::new(end))
    };

    let id = model.geometry_mut().add_curve(curve);
    let data = EdgeData::on_curve(id, Location::identity(), range);
    let edge = model.add_edge(data, &[start_vertex.clone(), end_vertex.clone()])?;

    model.set_derived(
        &start_vertex,
        std::slice::from_ref(&edge),
        roles::EDGE_START,
    )?;
    if !end_vertex.is_same(&start_vertex) {
        model.set_derived(&end_vertex, std::slice::from_ref(&edge), roles::EDGE_END)?;
    }

    let mut history = History::new();
    history.generate(&edge, start_vertex);
    if !end_vertex.is_same(&edge) {
        history.generate(&edge, end_vertex);
    }
    Ok(Built::new(edge, history))
}

/// Build an edge between two existing vertices, on `curve`.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if either shape is
/// not a vertex, or the curve's ends do not reach them within tolerance.
pub fn make_edge_between(
    model: &mut Model,
    curve: Curve,
    range: (f64, f64),
    from: &Shape,
    to: &Shape,
    tol: Tolerances,
) -> OgResult<Built> {
    for v in [from, to] {
        if model.kind_of(v)? != ShapeType::Vertex {
            og_bail!(Construction, "an edge is bounded by vertices");
        }
    }
    // The geometry has to actually reach the vertices it claims to join. An
    // edge whose curve stops short leaves a gap that only shows up later, in
    // whatever first tries to walk the boundary.
    let ends = [(range.0, from), (range.1, to)];
    for (parameter, vertex) in ends {
        let Some(node) = model.node(vertex) else {
            og_bail!(Dangling, "vertex is not in this model");
        };
        let Some(data) = node.data().as_vertex() else {
            og_bail!(Construction, "vertex node holds no point");
        };
        // Through the vertex's own placement, not against its stored point. A
        // vertex is a triple like any other shape, and the same node appears at
        // different places — the two ends of a prism are one vertex twice.
        // Comparing against the stored point would put both ends at the origin
        // of the placement and reject every edge that joins them.
        let placed = vertex.transform(model.datums())?.apply(data.point);
        let on_curve = curve.point_at(parameter, tol)?;
        let reach = data.tolerance.get().max(tol.confusion());
        if !on_curve.is_within(placed, reach) {
            og_bail!(
                Construction,
                "curve at {parameter} is {} from the vertex it should meet, \
                 outside its tolerance of {reach}",
                on_curve.distance(placed)
            );
        }
    }

    let id = model.geometry_mut().add_curve(curve);
    let data = EdgeData::on_curve(id, Location::identity(), range);
    let edge = model.add_edge(data, &[from.clone(), to.clone()])?;
    Ok(Built::from_nothing(edge))
}

/// The vertices an edge runs between, in the direction it is traversed.
///
/// A reversed edge runs the other way, so its start is its stored second
/// bound. Every caller that walks a boundary needs this, and getting it from
/// the raw children instead is how a wire ends up appearing disconnected.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if `edge` is not
/// an edge.
pub fn edge_vertices(model: &Model, edge: &Shape) -> OgResult<Option<(Shape, Shape)>> {
    if model.kind_of(edge)? != ShapeType::Edge {
        og_bail!(Construction, "expected an edge");
    }
    let bounds = model.children_of(edge)?;
    let (first, last) = match bounds.len() {
        0 => return Ok(None),
        1 => (bounds[0].clone(), bounds[0].clone()),
        _ => (bounds[0].clone(), bounds[bounds.len() - 1].clone()),
    };
    Ok(Some(
        if edge.orientation() == og_topo::Orientation::Reversed {
            (last, first)
        } else {
            (first, last)
        },
    ))
}

/// Build a wire of straight segments through a sequence of points.
///
/// `closed` adds a final segment back to the first point — and does it by
/// naming the *first vertex again* rather than making a coincident second one,
/// which is what keeps the wire closed under [`is_wire_closed`] rather than
/// merely looking closed.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if there are fewer
/// than two points, fewer than three for a closed polygon, if consecutive
/// points coincide within tolerance — a zero-length edge is not an edge — or if
/// a closed polygon's ends do not meet.
pub fn make_polygon(
    model: &mut Model,
    points: &[Point],
    closed: bool,
    tol: Tolerances,
) -> OgResult<Built> {
    let least = if closed { 3 } else { 2 };
    if points.len() < least {
        og_bail!(
            Construction,
            "a {} polygon needs at least {least} points, got {}",
            if closed { "closed" } else { "open" },
            points.len()
        );
    }
    for i in 1..points.len() {
        if points[i].is_equal(points[i - 1], tol) {
            og_bail!(
                Construction,
                "points {} and {i} coincide, so the edge between them has no \
                 length",
                i - 1
            );
        }
    }
    // Whether or not `closed` was asked for. A caller that repeated the first
    // point at the end wants a loop, and building it as written would give the
    // wire two vertices in the same place — which every later boundary walk
    // treats as a gap that happens to be zero wide. `closed` produces the loop
    // by naming the first vertex again, which is the thing that actually
    // closes.
    if points[0].is_equal(points[points.len() - 1], tol) {
        og_bail!(
            Construction,
            "the first and last points coincide; pass `closed` and drop the \
             repeat, or the wire gets two vertices in one place rather than a \
             closed loop"
        );
    }

    let mut vertices: Vec<Shape> = Vec::with_capacity(points.len());
    for point in points {
        vertices.push(model.add_vertex(VertexData::new(*point)));
    }

    let mut edges = Vec::with_capacity(points.len());
    let segments = if closed {
        points.len()
    } else {
        points.len() - 1
    };
    for i in 0..segments {
        let (from, to) = (points[i], points[(i + 1) % points.len()]);
        let curve: Curve = og_geom::LineCurve::segment(from, to, tol)?.into();
        edges.push(
            make_edge_between(
                model,
                curve,
                (0.0, from.distance(to)),
                &vertices[i],
                // The closing segment names the first vertex again rather than
                // a second one at the same place.
                &vertices[(i + 1) % points.len()],
                tol,
            )?
            .shape,
        );
    }
    make_wire(model, &edges, tol)
}

/// The plane a shape lies in, if it lies in one.
///
/// Fitted to the shape's geometry and then *checked*: the fit always produces a
/// plane, and the answer is only useful if every sample is within tolerance of
/// it. `None` means the shape is not planar, which is a fact about the shape
/// rather than a failure.
///
/// Curved edges are sampled along their length, not only at their ends. A
/// circular arc's endpoints lie in a great many planes that the arc does not.
///
/// # Errors
///
/// [`OgError::Dangling`](og_core::OgError::Dangling) if a handle fails to
/// resolve.
pub fn find_plane(
    model: &Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgResult<Option<og_math::Plane>> {
    let points = sample_shape(model, shape, tol)?;
    if points.len() < 3 {
        return Ok(None);
    }
    let Some((centroid, normal)) = crate::measure::least_squares_plane(&points, tol) else {
        return Ok(None);
    };
    // Fitting always yields a plane. Whether the shape is *in* it is the
    // question, and it is answered by measuring, not by having fitted.
    let reach = tol.confusion().max(
        points
            .iter()
            .map(|p| normal.dot_vector(*p - centroid).abs())
            .fold(0.0_f64, f64::max),
    );
    if reach > tol.confusion() {
        return Ok(None);
    }
    Ok(Some(og_math::Plane::new(og_math::Frame::about(
        centroid, normal,
    ))))
}

/// Points along a shape's geometry: every vertex, and every curved edge sampled
/// along its length.
fn sample_shape(model: &Model, shape: &Shape, tol: Tolerances) -> OgResult<Vec<Point>> {
    /// Enough to catch a curve leaving a candidate plane, without the cost of a
    /// real discretization — this is a yes-or-no question, not a mesh.
    const ALONG_EDGE: usize = 8;

    let mut points = Vec::new();
    for vertex in og_topo::explore_unique(model, shape, ShapeType::Vertex)? {
        if let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) {
            points.push(vertex.transform(model.datums())?.apply(data.point));
        }
    }
    for edge in og_topo::explore_unique(model, shape, ShapeType::Edge)? {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            continue;
        };
        if geometry.kind() == og_geom::CurveKind::Line {
            continue;
        }
        let placement = edge.transform(model.datums())?;
        for i in 1..ALONG_EDGE {
            #[allow(clippy::cast_precision_loss)]
            let t = i as f64 / ALONG_EDGE as f64;
            let at = range.0 + (range.1 - range.0) * t;
            points.push(placement.apply(geometry.point_at(at, tol)?));
        }
    }
    Ok(points)
}

/// Build a wire from edges that meet end to end.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if a child is not
/// an edge, the list is empty, or consecutive edges do not share a vertex.
pub fn make_wire(model: &mut Model, edges: &[Shape], tol: Tolerances) -> OgResult<Built> {
    if edges.is_empty() {
        og_bail!(Construction, "a wire needs at least one edge");
    }
    check_connected(model, edges, tol)?;

    let wire = model.add_wire(edges)?;
    let mut history = History::new();
    for edge in edges {
        // The edges are not consumed — they are shared, and an edge between two
        // faces belongs to both wires. Reporting them modified into the wire
        // would claim the edge ceased to exist.
        history.generate(edge, wire.clone());
    }
    Ok(Built::new(wire, history))
}

/// Verify that consecutive edges share a vertex.
fn check_connected(model: &Model, edges: &[Shape], tol: Tolerances) -> OgResult<()> {
    let mut ends: Vec<Option<(Shape, Shape)>> = Vec::with_capacity(edges.len());
    for edge in edges {
        ends.push(edge_vertices(model, edge)?);
    }
    for i in 0..edges.len().saturating_sub(1) {
        let (Some((_, end)), Some((next_start, _))) = (&ends[i], &ends[i + 1]) else {
            // An unbounded edge cannot be shown to join anything, and pretending
            // otherwise is worse than saying so.
            og_bail!(
                Construction,
                "edge {i} or {} has no bounding vertices, so the wire cannot be \
                 shown to connect",
                i + 1
            );
        };
        if !end.is_same(next_start) && !model.same_position(end, next_start, tol)? {
            og_bail!(
                Construction,
                "edge {i} ends where edge {} does not begin; a wire whose edges \
                 do not meet leaves a gap in every face built on it",
                i + 1
            );
        }
    }
    Ok(())
}

/// Whether a wire's last edge returns to its first edge's start.
///
/// # Errors
///
/// As [`edge_vertices`].
pub fn is_wire_closed(model: &Model, wire: &Shape, tol: Tolerances) -> OgResult<bool> {
    let edges = model.children_of(wire)?;
    let (Some(first), Some(last)) = (edges.first(), edges.last()) else {
        return Ok(false);
    };
    let (Some((start, _)), Some((_, end))) =
        (edge_vertices(model, first)?, edge_vertices(model, last)?)
    else {
        return Ok(false);
    };
    Ok(start.is_same(&end) || model.same_position(&start, &end, tol)?)
}

/// Build a face on `surface`, bounded by `wires`.
///
/// The first wire is the outer boundary; any others are holes. Every wire must
/// be closed, since an open boundary encloses nothing.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if a wire is open
/// or is not a wire.
pub fn make_face(
    model: &mut Model,
    surface: SurfaceGeometry,
    wires: &[Shape],
    tol: Tolerances,
) -> OgResult<Built> {
    for (i, wire) in wires.iter().enumerate() {
        if model.kind_of(wire)? != ShapeType::Wire {
            og_bail!(Construction, "a face is bounded by wires");
        }
        if !is_wire_closed(model, wire, tol)? {
            og_bail!(
                Construction,
                "wire {i} is open; an open boundary encloses no area"
            );
        }
    }

    let id = model.geometry_mut().add_surface(surface);
    make_face_on(model, id, wires, tol)
}

/// Build a face on a surface the model already holds.
///
/// The distinction from [`make_face`] matters more than it looks: a pcurve
/// names the surface it is drawn on by id, so an edge's pcurve and the face
/// bounded by that edge have to name the *same* id. Registering the surface
/// twice gives two ids for one surface, and every lookup of "the pcurve on this
/// face" comes back empty even though the pcurve is right there.
///
/// # Errors
///
/// As [`make_face`].
pub fn make_face_on(
    model: &mut Model,
    surface: og_topo::SurfaceId,
    wires: &[Shape],
    tol: Tolerances,
) -> OgResult<Built> {
    for (i, wire) in wires.iter().enumerate() {
        if model.kind_of(wire)? != ShapeType::Wire {
            og_bail!(Construction, "a face is bounded by wires");
        }
        if !is_wire_closed(model, wire, tol)? {
            og_bail!(
                Construction,
                "wire {i} is open; an open boundary encloses no area"
            );
        }
    }

    let data = FaceData::new(surface, Location::identity());
    let face = model.add_face(data, wires)?;

    let mut history = History::new();
    for wire in wires {
        history.generate(wire, face.clone());
    }
    Ok(Built::new(face, history))
}

/// Build a face on `surface` from per-wire edge lists, attaching an exact
/// same-parameter pcurve to every edge.
///
/// The construction path for faces whose curves were *chosen* to have
/// closed-form charts — blend wedges, offset rebuilds. Every edge's curve
/// must lie on the surface in a configuration
/// [`og_intersect::exact_pcurve_of`] recognises; a fitted pcurve here would
/// manufacture disagreement where none exists, so an edge with no closed
/// form is refused instead.
///
/// # Errors
///
/// As [`make_wire`] and [`make_face`], and
/// [`OgError::Construction`](og_core::OgError::Construction) if an edge has
/// no 3D curve or no closed-form pcurve on `surface`.
pub fn make_face_with_pcurves(
    model: &mut Model,
    surface: SurfaceGeometry,
    wires: &[Vec<Shape>],
    tol: Tolerances,
) -> OgResult<Built> {
    let mut rings: Vec<Shape> = Vec::with_capacity(wires.len());
    for edges in wires {
        rings.push(make_wire(model, edges, tol)?.shape);
    }
    let built = make_face(model, surface.clone(), &rings, tol)?;
    let surface_id = {
        let Some(node) = model.node(&built.shape) else {
            og_bail!(Dangling, "the face just built is not in this model");
        };
        let og_topo::NodeData::Face(data) = node.data() else {
            og_bail!(Construction, "the face holds no face data");
        };
        data.surface
    };
    for edge in og_topo::explore(
        model,
        &built.shape,
        og_topo::Filter::OfType(ShapeType::Edge),
    )? {
        let (curve, prange) = {
            let Some(node) = model.node(&edge) else {
                og_bail!(Dangling, "edge is not in this model");
            };
            let Some(data) = node.data().as_edge() else {
                og_bail!(Construction, "edge node holds no edge data");
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                og_bail!(Construction, "a face edge has no 3D curve");
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                og_bail!(Dangling, "curve is not in this model");
            };
            (geometry.clone(), *range)
        };
        let Some(pcurve) = og_intersect::exact_pcurve_of(&curve, &surface, tol) else {
            og_bail!(
                Construction,
                "a face edge has no closed-form pcurve on its surface"
            );
        };
        attach_pcurve(
            model,
            &edge,
            pcurve,
            surface_id,
            Location::identity(),
            prange,
        )?;
    }
    Ok(built)
}

/// Build a face covering the whole of `surface`, with no trimming.
///
/// # Errors
///
/// Never fails for a well-formed surface; the signature matches [`make_face`]
/// so the two are interchangeable at a call site.
pub fn make_natural_face(model: &mut Model, surface: SurfaceGeometry) -> OgResult<Built> {
    let id = model.geometry_mut().add_surface(surface);
    let face = model.add_face(FaceData::natural(id, Location::identity()), &[])?;
    Ok(Built::from_nothing(face))
}

/// Build a shell from faces.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if a child is not
/// a face, or the list is empty.
/// Build a compound from any shapes, with history.
///
/// The grouping container: a compound holds anything, orders nothing, and
/// claims nothing about closure. What this adds over the raw model call is
/// the same thing every builder adds — a history that says what went in.
///
/// # Errors
///
/// As [`Model::add_compound`].
pub fn make_compound(model: &mut Model, shapes: &[Shape]) -> OgResult<Built> {
    let compound = model.add_compound(shapes)?;
    let mut history = History::new();
    for shape in shapes {
        history.generate(shape, compound.clone());
    }
    Ok(Built::new(compound, history))
}

/// Build a shell from faces.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if a child is
/// not a face, or the list is empty.
pub fn make_shell(model: &mut Model, faces: &[Shape]) -> OgResult<Built> {
    let shell = model.add_shell(faces)?;
    let mut history = History::new();
    for face in faces {
        history.generate(face, shell.clone());
    }
    Ok(Built::new(shell, history))
}

/// Build a solid from shells.
///
/// The first shell is the outer boundary; any others are voids inside it.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if a child is not
/// a shell, or the list is empty.
pub fn make_solid(model: &mut Model, shells: &[Shape]) -> OgResult<Built> {
    let solid = model.add_solid(shells)?;
    let mut history = History::new();
    for shell in shells {
        history.generate(shell, solid.clone());
    }
    Ok(Built::new(solid, history))
}

/// Whether every edge in a shell is shared by exactly two faces.
///
/// The defining property of a closed shell, and the one that decides whether it
/// can bound a solid. An edge used once is a free boundary — the shell has a
/// hole. An edge used three or more times is non-manifold, which is legitimate
/// topology but not something that encloses a volume.
///
/// # Errors
///
/// As [`og_topo::explore`].
pub fn is_shell_closed(model: &Model, shell: &Shape) -> OgResult<bool> {
    // Counted by *use*, not by how many distinct faces an edge belongs to. A
    // seam edge bounds one face twice — up one side of its parameter rectangle
    // and down the other — so counting faces would call every cylinder, sphere
    // and torus open, which is precisely backwards.
    let mut uses: HashMap<og_topo::TShapeId, usize> = HashMap::new();
    for face in og_topo::explore(model, shell, og_topo::Filter::OfType(ShapeType::Face))? {
        for wire in model.children_of(&face)? {
            for edge in model.children_of(&wire)? {
                // A degenerate edge — a sphere's pole, a cone's apex — has no
                // length, so there is no gap along it for a second face to
                // close. Counting it would call every sphere and every true
                // cone open, and the thing that is actually open, isn't.
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
    Ok(!uses.is_empty() && uses.values().all(|n| n % 2 == 0))
}

/// Attach a pcurve to an edge, describing it in a surface's parameter space.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if `edge` is not
/// an edge; [`OgError::Dangling`](og_core::OgError::Dangling) if it does not
/// resolve.
pub fn attach_pcurve(
    model: &mut Model,
    edge: &Shape,
    pcurve: og_geom::PlanarCurve,
    surface: og_topo::SurfaceId,
    location: Location,
    range: (f64, f64),
) -> OgResult<()> {
    if model.kind_of(edge)? != ShapeType::Edge {
        og_bail!(Construction, "pcurves attach to edges");
    }
    let id = model.geometry_mut().add_pcurve(pcurve);
    let Some(node) = model.node_mut(edge) else {
        og_bail!(Dangling, "edge is not in this model");
    };
    let og_topo::NodeData::Edge(data) = node.data_mut() else {
        og_bail!(Construction, "edge node holds no edge data");
    };
    data.add(EdgeRepr::PCurve {
        curve: id,
        surface,
        location,
        range,
    });
    Ok(())
}

/// Attach a seam representation to an edge: one pcurve per side of a closed
/// surface's join.
///
/// The counterpart of [`attach_pcurve`] for the edge that bounds one face
/// twice — up one side of the parameter rectangle and down the other. Which
/// pcurve applies to an occurrence is decided by that occurrence's
/// orientation, which is the only thing distinguishing the two.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if `edge` is not
/// an edge; [`OgError::Dangling`](og_core::OgError::Dangling) if it is not in
/// this model.
pub fn attach_seam(
    model: &mut Model,
    edge: &Shape,
    forward: og_geom::PlanarCurve,
    reversed: og_geom::PlanarCurve,
    surface: og_topo::SurfaceId,
    location: Location,
    range: (f64, f64),
) -> OgResult<()> {
    if model.kind_of(edge)? != ShapeType::Edge {
        og_bail!(Construction, "seams attach to edges");
    }
    let forward = model.geometry_mut().add_pcurve(forward);
    let reversed = model.geometry_mut().add_pcurve(reversed);
    let Some(node) = model.node_mut(edge) else {
        og_bail!(Dangling, "edge is not in this model");
    };
    let og_topo::NodeData::Edge(data) = node.data_mut() else {
        og_bail!(Construction, "edge node holds no edge data");
    };
    data.add(EdgeRepr::Seam {
        forward,
        reversed,
        surface,
        location,
        range,
    });
    Ok(())
}

/// Build the face of a revolution band: two closed circle rings joined by a
/// synthesised seam, pcurves attached window-coherently.
///
/// The one authority for a job three call sites got subtly wrong three
/// different ways. The chart walk decides everything: the bottom ring is
/// traversed forward, and where its chart line *ends* is where the first
/// seam column stands; the top ring's occurrence direction is chosen so its
/// walk starts there — rings winding the same way traverse opposite, rings
/// winding opposite traverse alike — and the seam's two pcurves are assigned
/// to match which occurrence the triangulator will hand them to. The face is
/// built on a fresh copy of the surface, so no stale annotation from another
/// phase of the same rings can apply.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if a ring is
/// not a circle, the rings are not rings of this surface, or the surface's
/// iso-curve has no closed form.
pub fn make_revolution_band(
    model: &mut Model,
    surface: &og_geom::SurfaceGeometry,
    ring_lo: &Shape,
    ring_hi: &Shape,
    tol: Tolerances,
) -> OgResult<Shape> {
    use og_geom::Surface as _;

    let surface_id = model.geometry_mut().add_surface(surface.clone());
    let ((ua_dom, ub_dom), _) = surface.domain();
    let span = ub_dom - ua_dom;
    let axis_z = surface_iso_axis(surface)
        .ok_or_else(|| og_core::og_err!(Construction, "the surface has no revolution axis"))?;

    // Each ring's curve, range, winding, vertex, and chart row.
    struct Ring {
        edge: Shape,
        vertex: Shape,
        crange: (f64, f64),
        winding: f64,
        row: f64,
    }
    let mut rings = Vec::new();
    for used in [ring_lo, ring_hi] {
        // The caller's edges may carry orientation from their old wire uses;
        // the band is built on the curves' own directions, and the walk
        // chooses each occurrence's orientation itself.
        let edge = &if used.orientation() == og_topo::Orientation::Reversed {
            used.reversed()
        } else {
            used.clone()
        };
        let Some((vertex, other)) = edge_vertices(model, edge)? else {
            og_bail!(Construction, "a band ring has no vertex");
        };
        if !vertex.is_same(&other) {
            og_bail!(Construction, "a band ring is not closed");
        }
        let (curve, crange) = {
            let Some(node) = model.node(edge) else {
                og_bail!(Dangling, "edge is not in this model");
            };
            let Some(data) = node.data().as_edge() else {
                og_bail!(Construction, "edge node holds no edge data");
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                og_bail!(Construction, "a band ring has no curve");
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                og_bail!(Dangling, "curve is not in this model");
            };
            (geometry.clone(), *range)
        };
        let og_geom::Curve::Circle(c) = &curve else {
            og_bail!(Construction, "a band ring is not a circle");
        };
        let winding = c.circle().frame().z().vector().dot(axis_z).signum();
        let at = {
            let Some(node) = model.node(&vertex) else {
                og_bail!(Dangling, "vertex is not in this model");
            };
            let Some(data) = node.data().as_vertex() else {
                og_bail!(Construction, "vertex node holds no vertex data");
            };
            data.point
        };
        let (_, row) = crate::measure::project_on_surface(surface, at, 32, tol)?.parameters;
        rings.push(Ring {
            edge: edge.clone(),
            vertex,
            crange,
            winding,
            row,
        });
    }
    let anchor = {
        let Some(node) = model.node(&rings[0].vertex) else {
            og_bail!(Dangling, "vertex is not in this model");
        };
        let Some(data) = node.data().as_vertex() else {
            og_bail!(Construction, "vertex node holds no vertex data");
        };
        data.point
    };
    let (ua, _) = crate::measure::project_on_surface(surface, anchor, 32, tol)?.parameters;

    // Window-coherent ring pcurves: u(t) spans [ua, ua + span] whichever way
    // each ring winds.
    for ring in &rings {
        let u_start = if ring.winding > 0.0 { ua } else { ua + span };
        let origin = og_math::Point2::new(ring.winding.mul_add(-ring.crange.0, u_start), ring.row);
        let pcurve: og_geom::PlanarCurve = og_geom::Line2d::over(
            og_math::Axis2::new(
                origin,
                og_math::Direction2::new(og_math::Vector2::new(ring.winding, 0.0), tol)?,
            ),
            ring.crange.0,
            ring.crange.1,
        )?
        .into();
        attach_pcurve(
            model,
            &ring.edge,
            pcurve,
            surface_id,
            Location::identity(),
            ring.crange,
        )?;
    }

    // The seam runs along the surface's own iso-curve at the anchor angle,
    // parameterized by `v`, built along increasing `v`.
    let (va, vb) = (rings[0].row, rings[1].row);
    let Some(seam_curve) = surface_iso_u_curve(surface, ua, tol) else {
        og_bail!(
            Construction,
            "the surface's iso-curve has no closed form; no seam can be built"
        );
    };
    let (range, from, to, downward) = if va <= vb {
        (
            (va, vb),
            rings[0].vertex.clone(),
            rings[1].vertex.clone(),
            false,
        )
    } else {
        (
            (vb, va),
            rings[1].vertex.clone(),
            rings[0].vertex.clone(),
            true,
        )
    };
    let seam = make_edge_between(model, seam_curve, range, &from, &to, tol)?.shape;

    // The walk closes only if the top ring's traversal starts where the
    // bottom's ends, and the seam sides sit at the columns the walk visits.
    let bottom_end = if rings[0].winding > 0.0 {
        ua + span
    } else {
        ua
    };
    let other_col = if rings[0].winding > 0.0 {
        ua
    } else {
        ua + span
    };
    let hi_reversed = (rings[1].winding - rings[0].winding).abs() < 0.5;
    let column = |u: f64| -> OgResult<og_geom::PlanarCurve> {
        Ok(og_geom::Line2d::over(
            og_math::Axis2::new(
                og_math::Point2::new(u, 0.0),
                og_math::Direction2::new(og_math::Vector2::new(0.0, 1.0), tol)?,
            ),
            range.0 - 1.0,
            range.1 + 1.0,
        )?
        .into())
    };
    // The first seam occurrence in the wire is `up`; the triangulator hands
    // a Forward occurrence the `forward` pcurve. `up` is Forward exactly
    // when the seam was built upward.
    let (forward_col, reversed_col) = if downward {
        (other_col, bottom_end)
    } else {
        (bottom_end, other_col)
    };
    attach_seam(
        model,
        &seam,
        column(forward_col)?,
        column(reversed_col)?,
        surface_id,
        Location::identity(),
        range,
    )?;

    let up = if downward {
        seam.reversed()
    } else {
        seam.clone()
    };
    let top = if hi_reversed {
        rings[1].edge.reversed()
    } else {
        rings[1].edge.clone()
    };
    let ring = vec![rings[0].edge.clone(), up.clone(), top, up.reversed()];
    let wire = make_wire(model, &ring, tol)?.shape;
    Ok(make_face_on(model, surface_id, &[wire], tol)?.shape)
}

/// The revolution axis direction of a periodic analytic surface.
fn surface_iso_axis(surface: &og_geom::SurfaceGeometry) -> Option<og_math::Vector> {
    match surface {
        og_geom::SurfaceGeometry::Cylinder(c) => Some(c.cylinder().frame().z().vector()),
        og_geom::SurfaceGeometry::Cone(c) => Some(c.cone().frame().z().vector()),
        og_geom::SurfaceGeometry::Torus(t) => Some(t.torus().frame().z().vector()),
        og_geom::SurfaceGeometry::Sphere(s) => Some(s.sphere().frame().z().vector()),
        _ => None,
    }
}

/// The surface's `u = at` iso-curve, parameterized by `v` exactly.
///
/// A ruling on a cylinder, a tube circle on a torus — the curve a seam runs
/// along. `None` where no closed form exists. Shared by the STEP reader's
/// seam synthesis and the healer's ring re-anchoring, so the two cannot
/// disagree about what a seam is.
pub fn surface_iso_u_curve(
    surface: &SurfaceGeometry,
    at: f64,
    tol: Tolerances,
) -> Option<og_geom::Curve> {
    match surface {
        og_geom::SurfaceGeometry::Cylinder(c) => {
            let cylinder = c.cylinder();
            let frame = cylinder.frame();
            let radial = frame.x().vector() * at.cos() + frame.y().vector() * at.sin();
            let location = frame.origin() + radial * cylinder.radius();
            let axis = og_math::Axis {
                location,
                direction: frame.z(),
            };
            Some(og_geom::LineCurve::new(axis).into())
        }
        og_geom::SurfaceGeometry::Torus(t) => {
            let torus = t.torus();
            let frame = torus.frame();
            let radial = frame.x().vector() * at.cos() + frame.y().vector() * at.sin();
            let centre = frame.origin() + radial * torus.major_radius();
            // The tube circle framed so its own angle *is* the surface's v:
            // x toward the outer equator, y along the axis, which makes
            // z = x cross y the tangential direction.
            let circle_frame = og_math::Frame::new(
                centre,
                og_math::Direction::new(radial.cross(frame.z().vector()), tol).ok()?,
                og_math::Direction::new(radial, tol).ok()?,
                tol,
            )
            .ok()?;
            let circle = og_math::Circle::new(circle_frame, torus.minor_radius(), tol).ok()?;
            Some(og_geom::CircleCurve::new(circle).into())
        }
        _ => None,
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use og_geom::{CircleCurve, LineCurve, PlaneSurface};
    use og_math::{Circle, Frame, Plane};
    use og_topo::explore_unique;

    const T: Tolerances = Tolerances::millimetres();

    fn segment(a: Point, b: Point) -> Curve {
        LineCurve::segment(a, b, T).unwrap().into()
    }

    /// Four edges forming a closed square in the xy plane.
    fn square_edges(model: &mut Model) -> Vec<Shape> {
        let corners = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let mut vertices: Vec<Shape> = corners
            .iter()
            .map(|p| model.add_vertex(VertexData::new(*p)))
            .collect();
        vertices.push(vertices[0].clone());

        (0..4)
            .map(|i| {
                let curve = segment(corners[i], corners[(i + 1) % 4]);
                let length = corners[i].distance(corners[(i + 1) % 4]);
                make_edge_between(
                    model,
                    curve,
                    (0.0, length),
                    &vertices[i],
                    &vertices[i + 1],
                    T,
                )
                .unwrap()
                .shape
            })
            .collect()
    }

    #[test]
    fn an_edge_takes_its_vertices_from_its_own_geometry() {
        // Deriving them means the topology and the geometry cannot disagree
        // about where the edge starts and stops.
        let mut model = Model::new();
        let built = make_edge(
            &mut model,
            segment(Point::ORIGIN, Point::new(3.0, 4.0, 0.0)),
            (0.0, 5.0),
            T,
        )
        .unwrap();

        let (start, end) = edge_vertices(&model, &built.shape).unwrap().unwrap();
        let point_of = |v: &Shape| model.node(v).unwrap().data().as_vertex().unwrap().point;
        assert!(point_of(&start).is_equal(Point::ORIGIN, T));
        assert!(point_of(&end).is_equal(Point::new(3.0, 4.0, 0.0), T));
        assert_eq!(built.history.generated(&built.shape).len(), 2);
    }

    #[test]
    fn a_closed_edge_names_one_vertex_twice() {
        // Two coincident vertices would leave the wire looking open at the join,
        // and the gap only surfaces when something tries to walk the boundary.
        let mut model = Model::new();
        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 2.0, T).unwrap()).into();
        let built = make_edge(&mut model, circle, (0.0, core::f64::consts::TAU), T).unwrap();

        let (start, end) = edge_vertices(&model, &built.shape).unwrap().unwrap();
        assert!(start.is_same(&end), "a closed edge starts where it ends");
        assert_eq!(
            explore_unique(&model, &built.shape, ShapeType::Vertex)
                .unwrap()
                .len(),
            1
        );
    }

    #[test]
    fn an_edge_whose_curve_does_not_reach_its_vertices_is_refused() {
        // The gap would show up much later, in whatever first walks the
        // boundary, with nothing left to say where it came from.
        let mut model = Model::new();
        let a = model.add_vertex(VertexData::new(Point::ORIGIN));
        let b = model.add_vertex(VertexData::new(Point::new(10.0, 0.0, 0.0)));
        let curve = segment(Point::ORIGIN, Point::new(5.0, 0.0, 0.0));

        assert!(
            make_edge_between(&mut model, curve.clone(), (0.0, 5.0), &a, &b, T).is_err(),
            "the curve stops half way"
        );
        let c = model.add_vertex(VertexData::new(Point::new(5.0, 0.0, 0.0)));
        assert!(make_edge_between(&mut model, curve, (0.0, 5.0), &a, &c, T).is_ok());
    }

    #[test]
    fn a_vertex_with_a_loose_tolerance_admits_a_curve_that_stops_short() {
        // The reach is the vertex's own tolerance, not a fixed epsilon: a
        // vertex that has been widened by an earlier repair genuinely does
        // occupy that much space.
        let mut model = Model::new();
        let a = model.add_vertex(VertexData::new(Point::ORIGIN));
        let b = model
            .add_vertex(VertexData::with_tolerance(Point::new(5.001, 0.0, 0.0), 1e-2).unwrap());
        let curve = segment(Point::ORIGIN, Point::new(5.0, 0.0, 0.0));
        assert!(make_edge_between(&mut model, curve, (0.0, 5.0), &a, &b, T).is_ok());
    }

    #[test]
    fn edge_vertices_follow_the_direction_of_travel() {
        let mut model = Model::new();
        let built = make_edge(
            &mut model,
            segment(Point::ORIGIN, Point::new(1.0, 0.0, 0.0)),
            (0.0, 1.0),
            T,
        )
        .unwrap();

        let (forward_start, forward_end) = edge_vertices(&model, &built.shape).unwrap().unwrap();
        let (back_start, back_end) = edge_vertices(&model, &built.shape.reversed())
            .unwrap()
            .unwrap();

        assert!(back_start.is_same(&forward_end));
        assert!(back_end.is_same(&forward_start));
    }

    #[test]
    fn a_wire_whose_edges_do_not_meet_is_refused() {
        // Model::add_wire cannot catch this — it has no geometry. This is the
        // check that keeps a face from being built on a boundary with a gap.
        let mut model = Model::new();
        let joined = make_edge(
            &mut model,
            segment(Point::ORIGIN, Point::new(1.0, 0.0, 0.0)),
            (0.0, 1.0),
            T,
        )
        .unwrap()
        .shape;
        let apart = make_edge(
            &mut model,
            segment(Point::new(5.0, 0.0, 0.0), Point::new(6.0, 0.0, 0.0)),
            (0.0, 1.0),
            T,
        )
        .unwrap()
        .shape;

        let err = make_wire(&mut model, &[joined.clone(), apart], T).unwrap_err();
        assert!(err.to_string().contains("gap"), "unexpected message: {err}");
        assert!(make_wire(&mut model, &[joined], T).is_ok());
    }

    #[test]
    fn a_wire_of_edges_that_meet_is_accepted_and_reported_closed() {
        let mut model = Model::new();
        let edges = square_edges(&mut model);
        let wire = make_wire(&mut model, &edges, T).unwrap();
        assert!(is_wire_closed(&model, &wire.shape, T).unwrap());

        // The edges are generated into the wire, not consumed by it: an edge
        // between two faces belongs to both wires.
        for edge in &edges {
            assert!(!wire.history.is_deleted(edge));
            assert_eq!(wire.history.generated(edge).len(), 1);
        }
    }

    #[test]
    fn an_open_wire_is_reported_open() {
        let mut model = Model::new();
        let edges = square_edges(&mut model);
        let open = make_wire(&mut model, &edges[..3], T).unwrap();
        assert!(!is_wire_closed(&model, &open.shape, T).unwrap());
    }

    #[test]
    fn a_face_refuses_an_open_boundary() {
        let mut model = Model::new();
        let edges = square_edges(&mut model);
        let open = make_wire(&mut model, &edges[..3], T).unwrap().shape;
        let plane: SurfaceGeometry = PlaneSurface::new(Plane::new(Frame::WORLD)).into();

        let err = make_face(&mut model, plane.clone(), &[open], T).unwrap_err();
        assert!(
            err.to_string().contains("open"),
            "unexpected message: {err}"
        );

        let closed = make_wire(&mut model, &edges, T).unwrap().shape;
        assert!(make_face(&mut model, plane, &[closed], T).is_ok());
    }

    #[test]
    fn a_natural_face_needs_no_wires_at_all() {
        let mut model = Model::new();
        let built = make_natural_face(
            &mut model,
            PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        )
        .unwrap();
        let node = model.node(&built.shape).unwrap();
        assert!(node.data().as_face().unwrap().natural_restriction);
        assert!(built.history.is_empty());
    }

    #[test]
    fn a_shell_of_one_face_is_not_closed() {
        // A single face has free edges all round, so it bounds nothing.
        let mut model = Model::new();
        let edges = square_edges(&mut model);
        let wire = make_wire(&mut model, &edges, T).unwrap().shape;
        let face = make_face(
            &mut model,
            PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
            &[wire],
            T,
        )
        .unwrap()
        .shape;
        let shell = make_shell(&mut model, &[face]).unwrap();
        assert!(!is_shell_closed(&model, &shell.shape).unwrap());
    }

    #[test]
    fn two_faces_sharing_every_edge_form_a_closed_shell() {
        // The degenerate closed shell: the same square from both sides. Every
        // edge is used exactly twice, which is what closure means.
        let mut model = Model::new();
        let edges = square_edges(&mut model);
        let plane: SurfaceGeometry = PlaneSurface::new(Plane::new(Frame::WORLD)).into();

        let front_wire = make_wire(&mut model, &edges, T).unwrap().shape;
        let reversed: Vec<Shape> = edges.iter().rev().map(Shape::reversed).collect();
        let back_wire = make_wire(&mut model, &reversed, T).unwrap().shape;

        let front = make_face(&mut model, plane.clone(), &[front_wire], T)
            .unwrap()
            .shape;
        let back = make_face(&mut model, plane, &[back_wire], T).unwrap().shape;
        let shell = make_shell(&mut model, &[front, back]).unwrap().shape;

        assert!(is_shell_closed(&model, &shell).unwrap());
        assert_eq!(
            explore_unique(&model, &shell, ShapeType::Edge)
                .unwrap()
                .len(),
            4
        );
        assert!(make_solid(&mut model, &[shell]).is_ok());
    }

    #[test]
    fn a_pcurve_can_be_attached_after_the_edge_exists() {
        // An edge learns about a face's parameter space when it joins that
        // face, not when it is created — it may join several.
        let mut model = Model::new();
        let edge = make_edge(
            &mut model,
            segment(Point::ORIGIN, Point::new(1.0, 0.0, 0.0)),
            (0.0, 1.0),
            T,
        )
        .unwrap()
        .shape;

        let surface = model
            .geometry_mut()
            .add_surface(PlaneSurface::new(Plane::new(Frame::WORLD)).into());
        let pcurve =
            og_geom::Line2d::segment(og_math::Point2::ORIGIN, og_math::Point2::new(1.0, 0.0), T)
                .unwrap()
                .into();
        attach_pcurve(
            &mut model,
            &edge,
            pcurve,
            surface,
            Location::identity(),
            (0.0, 1.0),
        )
        .unwrap();

        let data = model.node(&edge).unwrap().data().as_edge().unwrap();
        assert_eq!(data.representations.len(), 2, "the curve and now a pcurve");
        assert!(data.pcurve_on(surface).is_some());
        assert!(
            !data.same_parameter(),
            "the new representation has not been shown to agree with the curve"
        );
    }

    #[test]
    fn degenerate_edge_ranges_are_refused() {
        let mut model = Model::new();
        let curve = segment(Point::ORIGIN, Point::new(1.0, 0.0, 0.0));
        assert!(make_edge(&mut model, curve.clone(), (1.0, 0.0), T).is_err());
        assert!(make_edge(&mut model, curve.clone(), (0.5, 0.5), T).is_err());
        assert!(make_edge(&mut model, curve.clone(), (0.0, f64::NAN), T).is_err());
        assert!(make_edge(&mut model, curve, (0.0, 1.0), T).is_ok());
    }

    #[test]
    fn builders_reject_children_of_the_wrong_type() {
        let mut model = Model::new();
        let vertex = model.add_vertex(VertexData::new(Point::ORIGIN));
        assert!(make_wire(&mut model, std::slice::from_ref(&vertex), T).is_err());
        assert!(make_shell(&mut model, std::slice::from_ref(&vertex)).is_err());
        assert!(make_solid(&mut model, &[vertex]).is_err());
        assert!(make_wire(&mut model, &[], T).is_err());
    }

    #[test]
    fn an_edges_vertices_record_where_they_came_from() {
        // Provenance is what a rebuild matches against; a vertex that does not
        // say which edge produced it cannot be found again.
        let mut model = Model::new();
        model.begin_operation();
        let built = make_edge(
            &mut model,
            segment(Point::ORIGIN, Point::new(1.0, 0.0, 0.0)),
            (0.0, 1.0),
            T,
        )
        .unwrap();

        let (start, end) = edge_vertices(&model, &built.shape).unwrap().unwrap();
        for (vertex, role) in [(start, roles::EDGE_START), (end, roles::EDGE_END)] {
            let provenance = model.provenance_of(&vertex).unwrap();
            assert!(
                matches!(provenance, og_core::Provenance::Derived { role: r, .. } if *r == role),
                "expected a derived vertex with role {role:?}, got {provenance:?}"
            );
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod polygon_tests {
    use super::*;
    use og_geom::{CircleCurve, PlaneSurface};
    use og_math::{Circle, Frame, Plane};
    use og_topo::explore_unique;

    const T: Tolerances = Tolerances::millimetres();

    fn square() -> Vec<Point> {
        vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
            Point::new(2.0, 2.0, 0.0),
            Point::new(0.0, 2.0, 0.0),
        ]
    }

    #[test]
    fn a_closed_polygon_names_its_first_vertex_again_rather_than_a_second_one() {
        // Two coincident vertices would leave the wire looking open at the
        // join, and the gap only surfaces when something walks the boundary.
        let mut model = Model::new();
        let built = make_polygon(&mut model, &square(), true, T).unwrap();

        assert_eq!(
            explore_unique(&model, &built.shape, ShapeType::Edge)
                .unwrap()
                .len(),
            4
        );
        assert_eq!(
            explore_unique(&model, &built.shape, ShapeType::Vertex)
                .unwrap()
                .len(),
            4,
            "four corners, not five"
        );
        assert!(is_wire_closed(&model, &built.shape, T).unwrap());
    }

    #[test]
    fn an_open_polygon_is_one_edge_short_and_says_it_is_open() {
        let mut model = Model::new();
        let built = make_polygon(&mut model, &square(), false, T).unwrap();
        assert_eq!(
            explore_unique(&model, &built.shape, ShapeType::Edge)
                .unwrap()
                .len(),
            3
        );
        assert!(!is_wire_closed(&model, &built.shape, T).unwrap());
    }

    #[test]
    fn a_polygon_that_describes_nothing_is_refused() {
        let mut model = Model::new();
        assert!(make_polygon(&mut model, &[], false, T).is_err());
        assert!(make_polygon(&mut model, &square()[..1], false, T).is_err());
        assert!(
            make_polygon(&mut model, &square()[..2], true, T).is_err(),
            "two points do not enclose anything"
        );

        // A repeated point is a zero-length edge, which is not an edge.
        let mut doubled = square();
        doubled.insert(2, doubled[1]);
        assert!(make_polygon(&mut model, &doubled, true, T).is_err());

        // Repeating the first point at the end is refused either way. Closed,
        // it would add a zero-length segment; open, it would leave the wire
        // with two vertices in one place, which every later boundary walk
        // reads as a gap that happens to be zero wide.
        let mut wrapped = square();
        wrapped.push(wrapped[0]);
        for closed in [true, false] {
            let err = make_polygon(&mut model, &wrapped, closed, T).unwrap_err();
            assert!(
                err.to_string().contains("first and last points coincide"),
                "unexpected message: {err}"
            );
        }
    }

    #[test]
    fn a_planar_shape_reports_the_plane_it_lies_in() {
        let mut model = Model::new();
        let wire = make_polygon(&mut model, &square(), true, T).unwrap().shape;
        let plane = find_plane(&model, &wire, T)
            .unwrap()
            .expect("a flat square");
        assert!(
            plane.normal().is_parallel(og_math::Direction::Z, T),
            "got {:?}",
            plane.normal()
        );
        for p in square() {
            assert!(plane.distance_to(p) < 1e-9);
        }
    }

    #[test]
    fn a_shape_that_is_not_planar_says_so_rather_than_fitting_one_anyway() {
        // Every set of three or more points has a best-fit plane, including a
        // set nowhere near one. Returning it would be a confident wrong answer.
        let mut model = Model::new();
        let mut skew = square();
        skew[2] = Point::new(2.0, 2.0, 1.0);
        let wire = make_polygon(&mut model, &skew, true, T).unwrap().shape;
        assert!(find_plane(&model, &wire, T).unwrap().is_none());
    }

    #[test]
    fn a_curved_edge_is_sampled_along_its_length_not_only_at_its_ends() {
        // An arc's endpoints lie in a great many planes the arc itself does
        // not. Checking only the vertices would call this shape planar.
        let mut model = Model::new();
        let circle = Circle::new(
            Frame::new(
                Point::ORIGIN,
                og_math::Direction::Z,
                og_math::Direction::X,
                T,
            )
            .unwrap(),
            2.0,
            T,
        )
        .unwrap();
        let arc = make_edge(
            &mut model,
            CircleCurve::new(circle).into(),
            (0.0, std::f64::consts::PI),
            T,
        )
        .unwrap()
        .shape;
        // In its own plane, it is planar.
        assert!(find_plane(&model, &arc, T).unwrap().is_some());

        // The two endpoints alone would admit the plane through them and the
        // z axis; the arc does not lie in it, and sampling catches that.
        let plane = find_plane(&model, &arc, T).unwrap().unwrap();
        assert!(plane.normal().is_parallel(og_math::Direction::Z, T));
    }

    #[test]
    fn a_face_is_planar_when_its_surface_is() {
        let mut model = Model::new();
        let wire = make_polygon(&mut model, &square(), true, T).unwrap().shape;
        let face = make_face(
            &mut model,
            PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
            std::slice::from_ref(&wire),
            T,
        )
        .unwrap()
        .shape;
        assert!(find_plane(&model, &face, T).unwrap().is_some());

        let solid = crate::make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        assert!(
            find_plane(&model, &solid, T).unwrap().is_none(),
            "a box is not planar"
        );
    }
}
