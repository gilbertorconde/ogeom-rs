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

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d, SurfaceGeometry};
use ogeom_math::Point;
use ogeom_topo::{EdgeData, EdgeRepr, FaceData, Location, Model, Shape, ShapeType, VertexData};

use crate::history::{Built, History};

/// Roles a builder assigns, so a rebuild can match entities up.
pub mod roles {
    use ogeom_core::Role;

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
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if `range` leaves the curve's
/// domain, and [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if it
/// is empty.
pub fn make_edge(
    model: &mut Model,
    curve: Curve,
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<Built> {
    let (lo, hi) = range;
    if !lo.is_finite() || !hi.is_finite() || hi <= lo + tol.parametric() {
        ogeom_bail!(Construction, "edge range [{lo}, {hi}] is empty");
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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if either shape is
/// not a vertex, or the curve's ends do not reach them within tolerance.
pub fn make_edge_between(
    model: &mut Model,
    curve: Curve,
    range: (f64, f64),
    from: &Shape,
    to: &Shape,
    tol: Tolerances,
) -> OgeomResult<Built> {
    for v in [from, to] {
        if model.kind_of(v)? != ShapeType::Vertex {
            ogeom_bail!(Construction, "an edge is bounded by vertices");
        }
    }
    // The geometry has to actually reach the vertices it claims to join. An
    // edge whose curve stops short leaves a gap that only shows up later, in
    // whatever first tries to walk the boundary.
    let ends = [(range.0, from), (range.1, to)];
    for (parameter, vertex) in ends {
        let Some(node) = model.node(vertex) else {
            ogeom_bail!(Dangling, "vertex is not in this model");
        };
        let Some(data) = node.data().as_vertex() else {
            ogeom_bail!(Construction, "vertex node holds no point");
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
            ogeom_bail!(
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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `edge` is not
/// an edge.
pub fn edge_vertices(model: &Model, edge: &Shape) -> OgeomResult<Option<(Shape, Shape)>> {
    if model.kind_of(edge)? != ShapeType::Edge {
        ogeom_bail!(Construction, "expected an edge");
    }
    let bounds = model.children_of(edge)?;
    let (first, last) = match bounds.len() {
        0 => return Ok(None),
        1 => (bounds[0].clone(), bounds[0].clone()),
        _ => (bounds[0].clone(), bounds[bounds.len() - 1].clone()),
    };
    Ok(Some(
        if edge.orientation() == ogeom_topo::Orientation::Reversed {
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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if there are fewer
/// than two points, fewer than three for a closed polygon, if consecutive
/// points coincide within tolerance — a zero-length edge is not an edge — or if
/// a closed polygon's ends do not meet.
pub fn make_polygon(
    model: &mut Model,
    points: &[Point],
    closed: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let least = if closed { 3 } else { 2 };
    if points.len() < least {
        ogeom_bail!(
            Construction,
            "a {} polygon needs at least {least} points, got {}",
            if closed { "closed" } else { "open" },
            points.len()
        );
    }
    for i in 1..points.len() {
        if points[i].is_equal(points[i - 1], tol) {
            ogeom_bail!(
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
        ogeom_bail!(
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
        let curve: Curve = ogeom_geom::LineCurve::segment(from, to, tol)?.into();
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
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle fails to
/// resolve.
pub fn find_plane(
    model: &Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<ogeom_math::Plane>> {
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
    Ok(Some(ogeom_math::Plane::new(ogeom_math::Frame::about(
        centroid, normal,
    ))))
}

/// Points along a shape's geometry: every vertex, and every curved edge sampled
/// along its length.
fn sample_shape(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Vec<Point>> {
    /// Enough to catch a curve leaving a candidate plane, without the cost of a
    /// real discretization — this is a yes-or-no question, not a mesh.
    const ALONG_EDGE: usize = 8;

    let mut points = Vec::new();
    for vertex in ogeom_topo::explore_unique(model, shape, ShapeType::Vertex)? {
        if let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) {
            points.push(vertex.transform(model.datums())?.apply(data.point));
        }
    }
    for edge in ogeom_topo::explore_unique(model, shape, ShapeType::Edge)? {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            continue;
        };
        if geometry.kind() == ogeom_geom::CurveKind::Line {
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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a child is not
/// an edge, the list is empty, or consecutive edges do not share a vertex.
pub fn make_wire(model: &mut Model, edges: &[Shape], tol: Tolerances) -> OgeomResult<Built> {
    if edges.is_empty() {
        ogeom_bail!(Construction, "a wire needs at least one edge");
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
fn check_connected(model: &Model, edges: &[Shape], tol: Tolerances) -> OgeomResult<()> {
    let mut ends: Vec<Option<(Shape, Shape)>> = Vec::with_capacity(edges.len());
    for edge in edges {
        ends.push(edge_vertices(model, edge)?);
    }
    for i in 0..edges.len().saturating_sub(1) {
        let (Some((_, end)), Some((next_start, _))) = (&ends[i], &ends[i + 1]) else {
            // An unbounded edge cannot be shown to join anything, and pretending
            // otherwise is worse than saying so.
            ogeom_bail!(
                Construction,
                "edge {i} or {} has no bounding vertices, so the wire cannot be \
                 shown to connect",
                i + 1
            );
        };
        if !end.is_same(next_start) && !model.same_position(end, next_start, tol)? {
            ogeom_bail!(
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
pub fn is_wire_closed(model: &Model, wire: &Shape, tol: Tolerances) -> OgeomResult<bool> {
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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a wire is open
/// or is not a wire.
pub fn make_face(
    model: &mut Model,
    surface: SurfaceGeometry,
    wires: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Built> {
    for (i, wire) in wires.iter().enumerate() {
        if model.kind_of(wire)? != ShapeType::Wire {
            ogeom_bail!(Construction, "a face is bounded by wires");
        }
        if !is_wire_closed(model, wire, tol)? {
            ogeom_bail!(
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
    surface: ogeom_topo::SurfaceId,
    wires: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Built> {
    for (i, wire) in wires.iter().enumerate() {
        if model.kind_of(wire)? != ShapeType::Wire {
            ogeom_bail!(Construction, "a face is bounded by wires");
        }
        if !is_wire_closed(model, wire, tol)? {
            ogeom_bail!(
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
/// [`ogeom_intersect::exact_pcurve_of`] recognises; a fitted pcurve here would
/// manufacture disagreement where none exists, so an edge with no closed
/// form is refused instead.
///
/// # Errors
///
/// As [`make_wire`] and [`make_face`], and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if an edge has
/// no 3D curve or no closed-form pcurve on `surface`.
pub fn make_face_with_pcurves(
    model: &mut Model,
    surface: SurfaceGeometry,
    wires: &[Vec<Shape>],
    tol: Tolerances,
) -> OgeomResult<Built> {
    let mut rings: Vec<Shape> = Vec::with_capacity(wires.len());
    for edges in wires {
        rings.push(make_wire(model, edges, tol)?.shape);
    }
    let built = make_face(model, surface.clone(), &rings, tol)?;
    let surface_id = {
        let Some(node) = model.node(&built.shape) else {
            ogeom_bail!(Dangling, "the face just built is not in this model");
        };
        let ogeom_topo::NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "the face holds no face data");
        };
        data.surface
    };
    // Walked in each wire's own traversal order, so that on a periodic
    // chart every pcurve lands on the *same branch* as its neighbour: an
    // exact inversion answers with whatever phase its construction likes,
    // and two images a full period apart describe the same points while
    // tearing the boundary the arrangement walks. Each pcurve after the
    // first is shifted by whole periods until its start meets the previous
    // traversal's end.
    use ogeom_geom::Curve2d as _;
    use ogeom_geom::Surface as _;
    let ((ua, ub), (va, vb)) = surface.domain();
    let u_period = surface.is_periodic_u().then_some(ub - ua);
    let v_period = surface.is_periodic_v().then_some(vb - va);
    let mut done: Vec<Shape> = Vec::new();
    for edges in wires {
        let mut prev_end: Option<ogeom_math::Point2> = None;
        for edge in edges {
            let (curve, prange) = {
                let Some(node) = model.node(edge) else {
                    ogeom_bail!(Dangling, "edge is not in this model");
                };
                let Some(data) = node.data().as_edge() else {
                    ogeom_bail!(Construction, "edge node holds no edge data");
                };
                let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                    ogeom_bail!(Construction, "a face edge has no 3D curve");
                };
                let Some(geometry) = model.geometry().curve(*curve) else {
                    ogeom_bail!(Dangling, "curve is not in this model");
                };
                (geometry.clone(), *range)
            };
            let mut pcurve = match ogeom_intersect::exact_pcurve_of(&curve, &surface, tol) {
                Some(exact) => exact,
                // No closed form — a fitted surface, mostly. The edge lies on
                // the surface by construction here (a rebuilt boundary, a
                // recovered intersection), so the projected fit speaks it: the
                // same machinery the exchange readers trust, at the same cap,
                // and the measured offset widens the edge honestly.
                None => {
                    let (fitted, _, _, worst_off, _) =
                        crate::pcurve_fit::fit_projected_pcurve(&curve, prange, &surface, tol)?;
                    if worst_off > tol.confusion()
                        && let Some(node) = model.node_mut(edge)
                        && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
                    {
                        data.tolerance = data.tolerance.widen_to(worst_off + tol.confusion());
                    }
                    fitted
                }
            };
            let (t_start, t_end) = if edge.orientation() == ogeom_topo::Orientation::Reversed {
                (prange.1, prange.0)
            } else {
                (prange.0, prange.1)
            };
            if u_period.is_some() || v_period.is_some() {
                let start = pcurve.point_at(t_start, tol)?;
                if let Some(prev) = prev_end {
                    let shift = ogeom_math::Vector2::new(
                        u_period.map_or(0.0, |p| ((prev.x - start.x) / p).round() * p),
                        v_period.map_or(0.0, |p| ((prev.y - start.y) / p).round() * p),
                    );
                    if shift.x != 0.0 || shift.y != 0.0 {
                        pcurve =
                            pcurve.transformed(&ogeom_math::Transform2::translation(shift), tol)?;
                    }
                }
                prev_end = Some(pcurve.point_at(t_end, tol)?);
            }
            if done.iter().any(|e| e.is_same(edge)) {
                continue;
            }
            done.push(edge.clone());
            attach_pcurve(
                model,
                edge,
                pcurve,
                surface_id,
                Location::identity(),
                prange,
            )?;
        }
    }
    Ok(built)
}

/// Build a face covering the whole of `surface`, with no trimming.
///
/// # Errors
///
/// Never fails for a well-formed surface; the signature matches [`make_face`]
/// so the two are interchangeable at a call site.
pub fn make_natural_face(model: &mut Model, surface: SurfaceGeometry) -> OgeomResult<Built> {
    let id = model.geometry_mut().add_surface(surface);
    let face = model.add_face(FaceData::natural(id, Location::identity()), &[])?;
    Ok(Built::from_nothing(face))
}

/// Build a shell from faces.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a child is not
/// a face, or the list is empty.
/// Build a compsolid: solids gluing along shared faces.
///
/// A compsolid is more than a bag of solids — its members must actually
/// glue. Every pair connected through the whole is connected through shared
/// *face nodes*: the same face entity bounding two solids, once from each
/// side, which is the sharing a boolean or a sew produces. A set of solids
/// that merely touch, faces coincident but distinct, is a compound, not a
/// compsolid, and is refused as one.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// fewer than two solids are given, a member is not a solid, or the members
/// do not glue into one connected whole through shared faces.
pub fn make_compsolid(model: &mut Model, solids: &[Shape]) -> OgeomResult<Built> {
    if solids.len() < 2 {
        ogeom_bail!(
            Construction,
            "a compsolid glues at least two solids; one solid is just a solid"
        );
    }
    // Which face nodes bound each solid.
    let mut faces_of: Vec<std::collections::HashSet<ogeom_topo::TShapeId>> = Vec::new();
    for solid in solids {
        if model.kind_of(solid)? != ShapeType::Solid {
            ogeom_bail!(Construction, "a compsolid's members are solids");
        }
        let set = ogeom_topo::explore(model, solid, ogeom_topo::Filter::OfType(ShapeType::Face))?
            .iter()
            .map(Shape::node)
            .collect();
        faces_of.push(set);
    }
    // Connectivity over shared face nodes: a flood fill from the first.
    let mut joined = vec![false; solids.len()];
    joined[0] = true;
    let mut frontier = vec![0_usize];
    while let Some(current) = frontier.pop() {
        for (other, other_faces) in faces_of.iter().enumerate() {
            if joined[other] {
                continue;
            }
            if !faces_of[current].is_disjoint(other_faces) {
                joined[other] = true;
                frontier.push(other);
            }
        }
    }
    if joined.iter().any(|j| !j) {
        ogeom_bail!(
            Construction,
            "the solids do not glue into one connected whole: at least one              shares no face with the rest. Coincident-but-distinct faces are              a compound's arrangement, not a compsolid's; sew or fuse first."
        );
    }

    let compsolid = model.add_compsolid(solids)?;
    let mut history = History::new();
    for solid in solids {
        history.generate(solid, compsolid.clone());
    }
    Ok(Built::new(compsolid, history))
}

/// Build a compound from any shapes, with history.
///
/// The grouping container: a compound holds anything, orders nothing, and
/// claims nothing about closure. What this adds over the raw model call is
/// the same thing every builder adds — a history that says what went in.
///
/// # Errors
///
/// As [`Model::add_compound`].
pub fn make_compound(model: &mut Model, shapes: &[Shape]) -> OgeomResult<Built> {
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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a child is
/// not a face, or the list is empty.
pub fn make_shell(model: &mut Model, faces: &[Shape]) -> OgeomResult<Built> {
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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a child is not
/// a shell, or the list is empty.
pub fn make_solid(model: &mut Model, shells: &[Shape]) -> OgeomResult<Built> {
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
/// As [`ogeom_topo::explore`].
pub fn is_shell_closed(model: &Model, shell: &Shape) -> OgeomResult<bool> {
    // Counted by *use*, not by how many distinct faces an edge belongs to. A
    // seam edge bounds one face twice — up one side of its parameter rectangle
    // and down the other — so counting faces would call every cylinder, sphere
    // and torus open, which is precisely backwards.
    let mut uses: HashMap<ogeom_topo::TShapeId, usize> = HashMap::new();
    for face in ogeom_topo::explore(model, shell, ogeom_topo::Filter::OfType(ShapeType::Face))? {
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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `edge` is not
/// an edge; [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if it does not
/// resolve.
pub fn attach_pcurve(
    model: &mut Model,
    edge: &Shape,
    pcurve: ogeom_geom::PlanarCurve,
    surface: ogeom_topo::SurfaceId,
    location: Location,
    range: (f64, f64),
) -> OgeomResult<()> {
    if model.kind_of(edge)? != ShapeType::Edge {
        ogeom_bail!(Construction, "pcurves attach to edges");
    }
    let id = model.geometry_mut().add_pcurve(pcurve);
    let Some(node) = model.node_mut(edge) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let ogeom_topo::NodeData::Edge(data) = node.data_mut() else {
        ogeom_bail!(Construction, "edge node holds no edge data");
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
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `edge` is not
/// an edge; [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if it is not in
/// this model.
pub fn attach_seam(
    model: &mut Model,
    edge: &Shape,
    forward: ogeom_geom::PlanarCurve,
    reversed: ogeom_geom::PlanarCurve,
    surface: ogeom_topo::SurfaceId,
    location: Location,
    range: (f64, f64),
) -> OgeomResult<()> {
    if model.kind_of(edge)? != ShapeType::Edge {
        ogeom_bail!(Construction, "seams attach to edges");
    }
    let forward = model.geometry_mut().add_pcurve(forward);
    let reversed = model.geometry_mut().add_pcurve(reversed);
    let Some(node) = model.node_mut(edge) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let ogeom_topo::NodeData::Edge(data) = node.data_mut() else {
        ogeom_bail!(Construction, "edge node holds no edge data");
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

/// Build the face of a revolution band: two closed rings joined by a
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
/// One ring may be *degenerate* — an edge with no curve, both ends the same
/// vertex, flagged as such — standing for an apex or a pole: a rim of no
/// length that still bounds the face in parameter space. It takes the row
/// the collapsed point sits on and traverses the chart opposite the real
/// ring, exactly as native cones and spheres bound their tips.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a ring is
/// neither a circle nor degenerate, both rings are degenerate, the rings are
/// not rings of this surface, or the surface's iso-curve has no closed form.
pub fn make_revolution_band(
    model: &mut Model,
    surface: &ogeom_geom::SurfaceGeometry,
    ring_lo: &Shape,
    ring_hi: &Shape,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    use ogeom_geom::Surface as _;

    let surface_id = model.geometry_mut().add_surface(surface.clone());
    let ((ua_dom, ub_dom), _) = surface.domain();
    let span = ub_dom - ua_dom;
    let axis_z = surface_iso_axis(surface).ok_or_else(|| {
        ogeom_core::ogeom_err!(Construction, "the surface has no revolution axis")
    })?;

    // Each ring's curve, range, winding, vertex, and chart row.
    struct Ring {
        edge: Shape,
        vertex: Shape,
        crange: (f64, f64),
        winding: f64,
        row: f64,
        degenerate: bool,
    }
    let mut rings = Vec::new();
    for used in [ring_lo, ring_hi] {
        // The caller's edges may carry orientation from their old wire uses;
        // the band is built on the curves' own directions, and the walk
        // chooses each occurrence's orientation itself.
        let edge = &if used.orientation() == ogeom_topo::Orientation::Reversed {
            used.reversed()
        } else {
            used.clone()
        };
        let Some((vertex, other)) = edge_vertices(model, edge)? else {
            ogeom_bail!(Construction, "a band ring has no vertex");
        };
        if !vertex.is_same(&other) {
            ogeom_bail!(Construction, "a band ring is not closed");
        }
        let at = {
            let Some(node) = model.node(&vertex) else {
                ogeom_bail!(Dangling, "vertex is not in this model");
            };
            let Some(data) = node.data().as_vertex() else {
                ogeom_bail!(Construction, "vertex node holds no vertex data");
            };
            data.point
        };
        let degenerate = {
            let Some(node) = model.node(edge) else {
                ogeom_bail!(Dangling, "edge is not in this model");
            };
            let Some(data) = node.data().as_edge() else {
                ogeom_bail!(Construction, "edge node holds no edge data");
            };
            data.degenerate
        };
        if degenerate {
            // An apex or a pole: no curve to read, no winding of its own.
            // Its row comes from the collapsed point — which sits on the
            // axis, where iterative projection has no nearest angle, so only
            // the closed-form inversion can place it.
            let Some(uv) = analytic_chart_of(surface, at) else {
                ogeom_bail!(
                    Construction,
                    "a degenerate ring's row cannot be found on this surface"
                );
            };
            rings.push(Ring {
                edge: edge.clone(),
                vertex,
                crange: (0.0, span),
                winding: 0.0,
                row: uv.y,
                degenerate: true,
            });
            continue;
        }
        let (curve, crange) = {
            let Some(node) = model.node(edge) else {
                ogeom_bail!(Dangling, "edge is not in this model");
            };
            let Some(data) = node.data().as_edge() else {
                ogeom_bail!(Construction, "edge node holds no edge data");
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                ogeom_bail!(Construction, "a band ring has no curve");
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                ogeom_bail!(Dangling, "curve is not in this model");
            };
            (geometry.clone(), *range)
        };
        let ogeom_geom::Curve::Circle(c) = &curve else {
            ogeom_bail!(Construction, "a band ring is not a circle");
        };
        let winding = c.circle().frame().z().vector().dot(axis_z).signum();
        let row = match analytic_chart_of(surface, at) {
            Some(uv) => uv.y,
            None => {
                crate::measure::project_on_surface(surface, at, 32, tol)?
                    .parameters
                    .1
            }
        };
        rings.push(Ring {
            edge: edge.clone(),
            vertex,
            crange,
            winding,
            row,
            degenerate: false,
        });
    }
    // A degenerate ring bounds nothing by itself, and the chart anchor must
    // come from a rim that has an angle: the real ring goes first, and the
    // degenerate one traverses opposite it.
    if rings[0].degenerate && rings[1].degenerate {
        ogeom_bail!(Construction, "both band rings are degenerate");
    }
    if rings[0].degenerate {
        rings.swap(0, 1);
    }
    if rings[1].degenerate {
        rings[1].winding = -rings[0].winding;
    }
    let anchor = {
        let Some(node) = model.node(&rings[0].vertex) else {
            ogeom_bail!(Dangling, "vertex is not in this model");
        };
        let Some(data) = node.data().as_vertex() else {
            ogeom_bail!(Construction, "vertex node holds no vertex data");
        };
        data.point
    };
    let ua = match analytic_chart_of(surface, anchor) {
        Some(uv) => uv.x,
        None => {
            crate::measure::project_on_surface(surface, anchor, 32, tol)?
                .parameters
                .0
        }
    };

    // Window-coherent ring pcurves: u(t) spans [ua, ua + span] whichever way
    // each ring winds.
    for ring in &rings {
        let u_start = if ring.winding > 0.0 { ua } else { ua + span };
        let origin =
            ogeom_math::Point2::new(ring.winding.mul_add(-ring.crange.0, u_start), ring.row);
        let pcurve: ogeom_geom::PlanarCurve = ogeom_geom::Line2d::over(
            ogeom_math::Axis2::new(
                origin,
                ogeom_math::Direction2::new(ogeom_math::Vector2::new(ring.winding, 0.0), tol)?,
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
        ogeom_bail!(
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
    // The edge lives on the curve's own parameterization; the chart rows map
    // onto it linearly, which is exactly the rescale a pcurve range states.
    let curve_range = (
        iso_curve_parameter_at(surface, range.0),
        iso_curve_parameter_at(surface, range.1),
    );
    let seam = make_edge_between(model, seam_curve, curve_range, &from, &to, tol)?.shape;

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
    let column = |u: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
        Ok(ogeom_geom::Line2d::over(
            ogeom_math::Axis2::new(
                ogeom_math::Point2::new(u, 0.0),
                ogeom_math::Direction2::new(ogeom_math::Vector2::new(0.0, 1.0), tol)?,
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

/// Build the face of a band between two closed rings that need not be
/// circles, each with a caller-supplied chart image.
///
/// [`make_revolution_band`]'s general sibling: where the band derives each
/// circle's row and pcurve itself, this takes the rings as they come — a
/// fitted tangency curve, a wavy trim — with the pcurve the caller already
/// knows for each, same-parameter over the edge's own range. The seam stays
/// an exact iso-column of the surface, which is what routes around the
/// closed-form refusal in [`make_face_with_pcurves`].
///
/// The caller's contract: both rings are closed, both chart images run the
/// full period *forward* in `u`, and both start on the same column — the
/// seam runs there, between the two start vertices.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a
/// ring is not closed, the images disagree about the start column or do not
/// run the period forward, or the surface's iso-curve has no closed form.
pub fn make_band_between(
    model: &mut Model,
    surface: &SurfaceGeometry,
    rings: [(&Shape, ogeom_geom::PlanarCurve); 2],
    tol: Tolerances,
) -> OgeomResult<Shape> {
    use ogeom_geom::Curve2d as _;
    use ogeom_geom::Surface as _;

    let surface_id = model.geometry_mut().add_surface(surface.clone());
    let ((ua_dom, ub_dom), _) = surface.domain();
    let span = ub_dom - ua_dom;

    struct Ring {
        edge: Shape,
        vertex: Shape,
        crange: (f64, f64),
        pcurve: ogeom_geom::PlanarCurve,
        start: ogeom_math::Point2,
    }
    let mut prepared: Vec<Ring> = Vec::with_capacity(2);
    for (used, pcurve) in rings {
        let edge = if used.orientation() == ogeom_topo::Orientation::Reversed {
            used.reversed()
        } else {
            used.clone()
        };
        let Some((vertex, other)) = edge_vertices(model, &edge)? else {
            ogeom_bail!(Construction, "a band ring has no vertex");
        };
        if !vertex.is_same(&other) {
            ogeom_bail!(Construction, "a band ring is not closed");
        }
        let crange = {
            let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                ogeom_bail!(Construction, "a band ring holds no edge data");
            };
            let Some(EdgeRepr::Curve3d { range, .. }) = data.curve3d() else {
                ogeom_bail!(Construction, "a band ring has no curve");
            };
            *range
        };
        let start = pcurve.point_at(crange.0, tol)?;
        let end = pcurve.point_at(crange.1, tol)?;
        if (end.x - start.x - span).abs() > 1e-3 {
            ogeom_bail!(
                Construction,
                "a band ring's chart image must run the period forward in u"
            );
        }
        prepared.push(Ring {
            edge,
            vertex,
            crange,
            pcurve,
            start,
        });
    }
    // The connector between the rings' starts. On one column it is the
    // surface's own iso; on different columns it is the straight chart
    // segment — a ruling, a parallel arc, or a helix, all exact on a
    // cylinder, and refused by name elsewhere.
    let (start0, start1) = (prepared[0].start, prepared[1].start);
    let dcol = start1.x - start0.x;
    // The iso is for rings genuinely sharing a column; starts even a whisker
    // apart take the chart segment, which meets both start vertices exactly
    // where an off-column iso would miss one by the offset times the radius.
    let (seam, up_is_forward, chart_from, chart_to) = if dcol.abs() <= 1e-12 {
        let column = start0.x;
        let (va, vb) = (start0.y, start1.y);
        let Some(seam_curve) = surface_iso_u_curve(surface, column, tol) else {
            ogeom_bail!(
                Construction,
                "the surface's iso-curve has no closed form; no seam can be built"
            );
        };
        let (range, from, to, downward) = if va <= vb {
            (
                (va, vb),
                prepared[0].vertex.clone(),
                prepared[1].vertex.clone(),
                false,
            )
        } else {
            (
                (vb, va),
                prepared[1].vertex.clone(),
                prepared[0].vertex.clone(),
                true,
            )
        };
        let curve_range = (
            iso_curve_parameter_at(surface, range.0),
            iso_curve_parameter_at(surface, range.1),
        );
        let seam = make_edge_between(model, seam_curve, curve_range, &from, &to, tol)?.shape;
        let (a, b) = (
            ogeom_math::Point2::new(column, range.0),
            ogeom_math::Point2::new(column, range.1),
        );
        (seam, !downward, a, b)
    } else {
        let SurfaceGeometry::Cylinder(c) = surface else {
            ogeom_bail!(
                Construction,
                "a band whose rings start on different columns needs a chart \
                 connector only a cylinder speaks exactly"
            );
        };
        let cylinder = c.cylinder();
        let frame = cylinder.frame();
        // The chart segment, parameterized by the angle so the pcurve's
        // linear map onto it is exact.
        let (a, b, from, to, forward) = if dcol > 0.0 {
            (
                start0,
                start1,
                prepared[0].vertex.clone(),
                prepared[1].vertex.clone(),
                true,
            )
        } else {
            (
                start1,
                start0,
                prepared[1].vertex.clone(),
                prepared[0].vertex.clone(),
                false,
            )
        };
        let drow = b.y - a.y;
        let connector: ogeom_geom::Curve = if drow.abs() <= 1e-9 {
            // A parallel: an arc at one height.
            let lifted = ogeom_math::Frame::new(
                frame.origin() + frame.z().vector() * a.y,
                frame.z(),
                frame.x(),
                tol,
            )?;
            ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(lifted, cylinder.radius(), tol)?)
                .into()
        } else {
            // A helix through both chart points: the pitch is the slope, and
            // the frame rides down so height zero lands where it must.
            let pitch = core::f64::consts::TAU * drow / (b.x - a.x);
            let dropped = ogeom_math::Frame::new(
                frame.origin()
                    + frame.z().vector() * (pitch / core::f64::consts::TAU).mul_add(-a.x, a.y),
                frame.z(),
                frame.x(),
                tol,
            )?;
            ogeom_geom::HelixCurve::over(dropped, cylinder.radius(), pitch, a.x, b.x)?.into()
        };
        let seam = make_edge_between(model, connector, (a.x, b.x), &from, &to, tol)?.shape;
        (seam, forward, a, b)
    };

    // Both rings run forward, so the walk closes as [lo, up, hi rev, down]:
    // the bottom ends at the far column, and the seam's two chart images sit
    // one period apart, assigned to whichever occurrence the walk meets
    // first.
    let segment_line = |shift: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
        let a = ogeom_math::Point2::new(chart_from.x + shift, chart_from.y);
        let b = ogeom_math::Point2::new(chart_to.x + shift, chart_to.y);
        let d = ogeom_math::Vector2::new(b.x - a.x, b.y - a.y);
        let m = d.x.hypot(d.y);
        if m <= 1e-12 {
            ogeom_bail!(Construction, "the seam's chart image has no length");
        }
        Ok(ogeom_geom::Line2d::over(
            ogeom_math::Axis2::new(a, ogeom_math::Direction2::new(d, tol)?),
            0.0,
            m,
        )?
        .into())
    };
    let (forward_shift, reversed_shift) = if up_is_forward {
        (span, 0.0)
    } else {
        (0.0, span)
    };
    // The pcurves run over their own arc length, start-for-start with the
    // seam's curve range, so that is the window the attachment states.
    let seam_length = (chart_to.x - chart_from.x).hypot(chart_to.y - chart_from.y);
    attach_seam(
        model,
        &seam,
        segment_line(forward_shift)?,
        segment_line(reversed_shift)?,
        surface_id,
        Location::identity(),
        (0.0, seam_length),
    )?;
    for ring in &prepared {
        attach_pcurve(
            model,
            &ring.edge,
            ring.pcurve.clone(),
            surface_id,
            Location::identity(),
            ring.crange,
        )?;
    }

    let up = if up_is_forward {
        seam.clone()
    } else {
        seam.reversed()
    };
    let walk = vec![
        prepared[0].edge.clone(),
        up.clone(),
        prepared[1].edge.reversed(),
        up.reversed(),
    ];
    let wire = make_wire(model, &walk, tol)?.shape;
    Ok(make_face_on(model, surface_id, &[wire], tol)?.shape)
}

/// Build the face of a revolution cap: one closed ring belting a cone, with
/// the apex the file never wrote synthesised as the degenerate ring the band
/// needs.
///
/// A cone face bounded by a single circle has exactly one other boundary the
/// geometry permits — the apex — because the region away from the apex is
/// unbounded. The apex becomes a vertex and a degenerate edge, and the rest
/// is [`make_revolution_band`], one authority for the seam either way.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the surface
/// is not a cone, or the ring resists the band construction.
pub fn make_apex_band(
    model: &mut Model,
    surface: &SurfaceGeometry,
    ring: &Shape,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let SurfaceGeometry::Cone(c) = surface else {
        ogeom_bail!(
            Construction,
            "only a cone has the single apex a one-ring cap implies"
        );
    };
    let apex = model.add_vertex(VertexData::new(c.cone().apex()));
    let mut data = EdgeData::new();
    data.degenerate = true;
    let edge = model.add_edge(data, &[apex.clone(), apex])?;
    make_revolution_band(model, surface, ring, &edge, tol)
}

/// The chart coordinates of a point on a periodic analytic surface, by
/// closed-form inversion, folded into the surface's own stated window.
///
/// The band synthesis needs the *row* each ring stands on and the *column*
/// the seam anchors to; iterative projection over an imported surface's
/// enormous stated extents can converge to a clamped boundary and place a
/// seam whole units from the vertex it must meet, so the analytic kinds
/// invert exactly instead.
fn analytic_chart_of(
    surface: &SurfaceGeometry,
    p: ogeom_math::Point,
) -> Option<ogeom_math::Point2> {
    use ogeom_geom::Surface as _;
    let tau = core::f64::consts::TAU;
    let raw = match surface {
        SurfaceGeometry::Cylinder(s) => {
            let l = s.cylinder().frame().to_local(p);
            ogeom_math::Point2::new(l.y.atan2(l.x), l.z)
        }
        SurfaceGeometry::Cone(s) => {
            let l = s.cone().frame().to_local(p);
            ogeom_math::Point2::new(l.y.atan2(l.x), l.z)
        }
        SurfaceGeometry::Sphere(s) => {
            let sphere = s.sphere();
            let l = sphere.frame().to_local(p);
            let lat = (l.z / sphere.radius()).clamp(-1.0, 1.0).asin();
            ogeom_math::Point2::new(l.y.atan2(l.x), lat)
        }
        SurfaceGeometry::Torus(s) => {
            let torus = s.torus();
            let l = torus.frame().to_local(p);
            let radial = l.x.hypot(l.y) - torus.major_radius();
            ogeom_math::Point2::new(l.y.atan2(l.x), l.z.atan2(radial))
        }
        _ => return None,
    };
    let ((ua, _), (va, vb)) = surface.domain();
    let u = ua + (raw.x - ua).rem_euclid(tau);
    let v = if surface.is_periodic_v() {
        va + (raw.y - va).rem_euclid(vb - va)
    } else {
        raw.y
    };
    Some(ogeom_math::Point2::new(u, v))
}

/// The revolution axis direction of a periodic analytic surface.
fn surface_iso_axis(surface: &ogeom_geom::SurfaceGeometry) -> Option<ogeom_math::Vector> {
    match surface {
        ogeom_geom::SurfaceGeometry::Cylinder(c) => Some(c.cylinder().frame().z().vector()),
        ogeom_geom::SurfaceGeometry::Cone(c) => Some(c.cone().frame().z().vector()),
        ogeom_geom::SurfaceGeometry::Torus(t) => Some(t.torus().frame().z().vector()),
        ogeom_geom::SurfaceGeometry::Sphere(s) => Some(s.sphere().frame().z().vector()),
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
) -> Option<ogeom_geom::Curve> {
    match surface {
        ogeom_geom::SurfaceGeometry::Cylinder(c) => {
            let cylinder = c.cylinder();
            let frame = cylinder.frame();
            let radial = frame.x().vector() * at.cos() + frame.y().vector() * at.sin();
            let location = frame.origin() + radial * cylinder.radius();
            let axis = ogeom_math::Axis {
                location,
                direction: frame.z(),
            };
            Some(ogeom_geom::LineCurve::new(axis).into())
        }
        ogeom_geom::SurfaceGeometry::Torus(t) => {
            let torus = t.torus();
            let frame = torus.frame();
            let radial = frame.x().vector() * at.cos() + frame.y().vector() * at.sin();
            let centre = frame.origin() + radial * torus.major_radius();
            // The tube circle framed so its own angle *is* the surface's v:
            // x toward the outer equator, y along the axis, which makes
            // z = x cross y the tangential direction.
            let circle_frame = ogeom_math::Frame::new(
                centre,
                ogeom_math::Direction::new(radial.cross(frame.z().vector()), tol).ok()?,
                ogeom_math::Direction::new(radial, tol).ok()?,
                tol,
            )
            .ok()?;
            let circle = ogeom_math::Circle::new(circle_frame, torus.minor_radius(), tol).ok()?;
            Some(ogeom_geom::CircleCurve::new(circle).into())
        }
        ogeom_geom::SurfaceGeometry::Cone(c) => {
            let cone = c.cone();
            let frame = cone.frame();
            let radial = frame.x().vector() * at.cos() + frame.y().vector() * at.sin();
            // The ruling through v = 0, arc-length parameterized: the chart's
            // v maps onto it linearly, by t = v / cos(half angle), which the
            // seam construction recovers through `iso_curve_parameter_at`.
            let location = frame.origin() + radial * cone.radius_at(0.0);
            let direction = ogeom_math::Direction::new(
                frame.z().vector() + radial * cone.half_angle().tan(),
                tol,
            )
            .ok()?;
            Some(
                ogeom_geom::LineCurve::new(ogeom_math::Axis {
                    location,
                    direction,
                })
                .into(),
            )
        }
        ogeom_geom::SurfaceGeometry::Sphere(sp) => {
            let sphere = sp.sphere();
            let frame = sphere.frame();
            let radial = frame.x().vector() * at.cos() + frame.y().vector() * at.sin();
            // The meridian framed so its own angle is the latitude exactly:
            // x one radius out along the parallel, y toward the north pole.
            let circle_frame = ogeom_math::Frame::new(
                frame.origin(),
                ogeom_math::Direction::new(radial.cross(frame.z().vector()), tol).ok()?,
                ogeom_math::Direction::new(radial, tol).ok()?,
                tol,
            )
            .ok()?;
            let circle = ogeom_math::Circle::new(circle_frame, sphere.radius(), tol).ok()?;
            Some(ogeom_geom::CircleCurve::new(circle).into())
        }
        _ => None,
    }
}

/// The parameter on a surface's iso-curve that lands at chart row `v`.
///
/// Identity for the kinds whose iso-curve is parameterized by `v` itself —
/// a cylinder ruling, a torus tube circle, a sphere meridian — and the slant
/// rescale for a cone, whose ruling is arc-length parameterized while the
/// chart's `v` is the height.
fn iso_curve_parameter_at(surface: &SurfaceGeometry, v: f64) -> f64 {
    match surface {
        ogeom_geom::SurfaceGeometry::Cone(c) => v / c.cone().half_angle().cos(),
        _ => v,
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;

    /// A face on a fitted surface builds with fitted trims: the closed-form
    /// refusal falls back to the projected fit, which is what lets a
    /// defeaturing rebuild stand a face on a spline neighbour (issue #19).
    #[test]
    fn a_face_on_a_fitted_surface_gains_projected_pcurves() {
        let tol = ogeom_core::Tolerances::millimetres();
        let mut model = Model::new();
        // A gently wavy fitted patch, nothing analytic recognises.
        let rows: Vec<Vec<ogeom_math::Point>> = (0..5)
            .map(|j| {
                (0..5)
                    .map(|i| {
                        let (x, y) = (f64::from(i) * 2.5, f64::from(j) * 2.5);
                        ogeom_math::Point::new(x, y, 0.4 * (x * 0.7).sin() * (y * 0.5).cos())
                    })
                    .collect()
            })
            .collect();
        let surface: ogeom_geom::SurfaceGeometry =
            ogeom_geom::fit::fit_surface_grid(&rows, 3, 1e-6, tol)
                .unwrap()
                .curve
                .into();
        // The border iso-curves, straight off the fitted chart.
        use ogeom_geom::Surface as _;
        let ((ua, ub), (va, vb)) = surface.domain();
        let iso = |fixed_u: Option<f64>, fixed_v: Option<f64>| -> Curve {
            const N: usize = 33;
            let pts: Vec<ogeom_math::Point> = (0..N)
                .map(|k| {
                    let t = f64::from(u32::try_from(k).unwrap())
                        / f64::from(u32::try_from(N - 1).unwrap());
                    let (u, v) = match (fixed_u, fixed_v) {
                        (Some(u), None) => (u, va + (vb - va) * t),
                        (None, Some(v)) => (ua + (ub - ua) * t, v),
                        _ => unreachable!(),
                    };
                    surface.point_at(u, v, tol).unwrap()
                })
                .collect();
            Curve::BSpline(
                ogeom_geom::fit::fit_points(&pts, 3, 1e-9, tol)
                    .unwrap()
                    .curve,
            )
        };
        // Shared corner vertices, widened to absorb the border fits' own
        // slack, so the ring connects by node rather than by luck.
        let corners: Vec<Shape> = [(ua, va), (ub, va), (ub, vb), (ua, vb)]
            .into_iter()
            .map(|(u, v)| {
                let p = surface.point_at(u, v, tol).unwrap();
                let vertex = model.add_vertex(ogeom_topo::VertexData::new(p));
                if let Some(node) = model.node_mut(&vertex)
                    && let ogeom_topo::NodeData::Vertex(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen_to(1e-4);
                }
                vertex
            })
            .collect();
        let ring = [
            (iso(None, Some(va)), 0usize, 1usize),
            (iso(Some(ub), None), 1, 2),
            (iso(None, Some(vb)), 3, 2),
            (iso(Some(ua), None), 0, 3),
        ];
        let edges: Vec<Shape> = ring
            .into_iter()
            .enumerate()
            .map(|(k, (c, from, to))| {
                use ogeom_geom::Curve3d as _;
                let domain = c.domain();
                let edge =
                    make_edge_between(&mut model, c, domain, &corners[from], &corners[to], tol)
                        .unwrap()
                        .shape;
                if k >= 2 { edge.reversed() } else { edge }
            })
            .collect();
        let built =
            make_face_with_pcurves(&mut model, surface, std::slice::from_ref(&edges), tol).unwrap();
        // Every edge carries a pcurve for the face's surface now.
        let data = model
            .node(&built.shape)
            .unwrap()
            .data()
            .as_face()
            .unwrap()
            .clone();
        for edge in &edges {
            let e = model.node(edge).unwrap().data().as_edge().unwrap();
            assert!(
                e.pcurve_for(data.surface, edge.location()).is_some(),
                "a fitted trim was attached"
            );
        }
    }

    #[test]
    fn a_band_whose_rings_start_apart_gets_a_helical_connector() {
        // Two circles at different heights whose parameter origins differ:
        // the seam has no single column to run down, so it runs the straight
        // chart segment — a helix on the cylinder, exactly.
        let mut model = Model::new();
        let radius = 2.0;
        let cylinder = ogeom_geom::CylinderSurface::new(
            ogeom_math::Cylinder::new(Frame::WORLD, radius, T).unwrap(),
            (-1.0, 4.0),
        )
        .unwrap();
        let surface: SurfaceGeometry = cylinder.into();

        let ring_at = |model: &mut Model, height: f64, phase: f64| {
            let x = ogeom_math::Direction::new(
                ogeom_math::Vector::new(phase.cos(), phase.sin(), 0.0),
                T,
            )
            .unwrap();
            let frame =
                Frame::new(Point::new(0.0, 0.0, height), ogeom_math::Direction::Z, x, T).unwrap();
            let curve: Curve = CircleCurve::new(Circle::new(frame, radius, T).unwrap()).into();
            let domain = curve.domain();
            let edge = make_edge(model, curve, domain, T).unwrap().shape;
            let pcurve: ogeom_geom::PlanarCurve = ogeom_geom::Line2d::over(
                ogeom_math::Axis2::new(
                    ogeom_math::Point2::new(phase, height),
                    ogeom_math::Direction2::new(ogeom_math::Vector2::new(1.0, 0.0), T).unwrap(),
                ),
                0.0,
                core::f64::consts::TAU,
            )
            .unwrap()
            .into();
            (edge, pcurve)
        };
        let (lo, lo_pcurve) = ring_at(&mut model, 0.0, 0.0);
        let (hi, hi_pcurve) = ring_at(&mut model, 2.0, 0.4);

        let band = make_band_between(
            &mut model,
            &surface,
            [(&lo, lo_pcurve), (&hi, hi_pcurve)],
            T,
        )
        .unwrap();

        assert_eq!(
            explore_unique(&model, &band, ShapeType::Edge)
                .unwrap()
                .len(),
            3
        );
        // The connector's curve is a genuine helix through both starts.
        let helix = explore_unique(&model, &band, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find_map(|e| {
                let data = model.node(&e)?.data().as_edge()?;
                let ogeom_topo::EdgeRepr::Curve3d { curve, .. } = data.curve3d()? else {
                    return None;
                };
                match model.geometry().curve(*curve)? {
                    Curve::Helix(h) => Some(*h),
                    _ => None,
                }
            })
            .expect("the connector is a helix");
        use ogeom_geom::Curve3d as _;
        let at_start = Curve::Helix(helix).point_at(0.0, T).unwrap();
        let at_end = Curve::Helix(helix).point_at(0.4, T).unwrap();
        assert!(at_start.distance(Point::new(2.0, 0.0, 0.0)) < 1e-9);
        assert!(at_end.distance(Point::new(2.0 * 0.4_f64.cos(), 2.0 * 0.4_f64.sin(), 2.0)) < 1e-9);
    }

    #[test]
    fn a_band_between_a_circle_and_a_fitted_wavy_ring_closes_with_its_seam() {
        use ogeom_geom::Curve2d as _;
        use ogeom_geom::Curve3d as _;
        use ogeom_geom::Surface as _;

        let mut model = Model::new();
        let radius = 2.0;
        let cylinder = ogeom_geom::CylinderSurface::new(
            ogeom_math::Cylinder::new(Frame::WORLD, radius, T).unwrap(),
            (-1.0, 4.0),
        )
        .unwrap();
        let surface: SurfaceGeometry = cylinder.into();

        // The lower ring: the exact circle at z = 0, its chart image the row.
        let lo_curve: Curve =
            CircleCurve::new(Circle::new(Frame::WORLD, radius, T).unwrap()).into();
        let lo_domain = lo_curve.domain();
        let lo = make_edge(&mut model, lo_curve, lo_domain, T).unwrap().shape;
        let lo_pcurve: ogeom_geom::PlanarCurve = ogeom_geom::Line2d::over(
            ogeom_math::Axis2::new(
                ogeom_math::Point2::new(0.0, 0.0),
                ogeom_math::Direction2::new(ogeom_math::Vector2::new(1.0, 0.0), T).unwrap(),
            ),
            lo_domain.0,
            lo_domain.1,
        )
        .unwrap()
        .into();

        // The upper ring: a wavy loop on the cylinder, fitted jointly with
        // its own chart image so curve and pcurve stay same-parameter.
        let tau = core::f64::consts::TAU;
        let samples = 96;
        let mut points = Vec::with_capacity(samples + 1);
        let mut chart = Vec::with_capacity(samples + 1);
        for i in 0..=samples {
            #[allow(clippy::cast_precision_loss)]
            let u = tau * (i as f64) / (samples as f64);
            let v = (3.0 * u).sin().mul_add(0.2, 1.5);
            points.push(Point::new(radius * u.cos(), radius * u.sin(), v));
            chart.push(ogeom_math::Point2::new(u, v));
        }
        let (fitted, hi_chart, _) =
            ogeom_geom::fit::fit_points_joint_closed(&points, &chart, &chart, 3, 1e-4, T).unwrap();
        assert!(fitted.met, "the wavy ring fit should meet its tolerance");
        let hi_curve: Curve = fitted.curve.into();
        let hi_domain = hi_curve.domain();
        let hi = make_edge(&mut model, hi_curve.clone(), hi_domain, T)
            .unwrap()
            .shape;

        let band = make_band_between(
            &mut model,
            &surface,
            [(&lo, lo_pcurve), (&hi, hi_chart.into())],
            T,
        )
        .unwrap();

        // Four traversals, three distinct edges: two rings and the seam,
        // used once per side.
        let occurrences =
            ogeom_topo::explore(&model, &band, ogeom_topo::Filter::OfType(ShapeType::Edge))
                .unwrap();
        assert_eq!(occurrences.len(), 4);
        assert_eq!(
            explore_unique(&model, &band, ShapeType::Edge)
                .unwrap()
                .len(),
            3
        );

        // The wavy ring's chart image lands on the curve itself: the fit is
        // same-parameter, and the surface evaluates the image back onto it.
        let hi_pcurve = {
            let data = model.node(&hi).unwrap().data().as_edge().unwrap();
            let repr = data
                .representations
                .iter()
                .find_map(|r| match r {
                    ogeom_topo::EdgeRepr::PCurve { curve, .. } => Some(*curve),
                    _ => None,
                })
                .expect("the ring carries its pcurve");
            model.geometry().pcurve(repr).unwrap().clone()
        };
        for i in 0..8 {
            #[allow(clippy::cast_precision_loss)]
            let t = hi_domain.0 + (hi_domain.1 - hi_domain.0) * (i as f64) / 8.0;
            let uv = hi_pcurve.point_at(t, T).unwrap();
            let through_chart = surface.point_at(uv.x, uv.y, T).unwrap();
            let direct = hi_curve.point_at(t, T).unwrap();
            assert!(
                through_chart.distance(direct) < 1e-4,
                "same-parameter drift at t = {t}"
            );
        }
    }

    use ogeom_geom::{CircleCurve, LineCurve, PlaneSurface};
    use ogeom_math::{Circle, Frame, Plane};
    use ogeom_topo::explore_unique;

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
        let pcurve = ogeom_geom::Line2d::segment(
            ogeom_math::Point2::ORIGIN,
            ogeom_math::Point2::new(1.0, 0.0),
            T,
        )
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
                matches!(provenance, ogeom_core::Provenance::Derived { role: r, .. } if *r == role),
                "expected a derived vertex with role {role:?}, got {provenance:?}"
            );
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod polygon_tests {
    use super::*;
    use ogeom_geom::{CircleCurve, PlaneSurface};
    use ogeom_math::{Circle, Frame, Plane};
    use ogeom_topo::explore_unique;

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
            plane.normal().is_parallel(ogeom_math::Direction::Z, T),
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
                ogeom_math::Direction::Z,
                ogeom_math::Direction::X,
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
        assert!(plane.normal().is_parallel(ogeom_math::Direction::Z, T));
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
#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod compsolid_tests {
    use super::*;
    use ogeom_geom::PlaneSurface;
    use ogeom_math::{Frame, Plane, Point};

    const T: Tolerances = Tolerances::millimetres();

    /// Two unit cubes glued along one shared face node at `x = 1`.
    fn glued_cubes(model: &mut Model) -> (Shape, Shape) {
        let a = crate::make_box(model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        // The face of `a` at x = 1: the one whose every vertex has x = 1.
        let shared = ogeom_topo::explore(model, &a, ogeom_topo::Filter::OfType(ShapeType::Face))
            .unwrap()
            .into_iter()
            .find(|f| {
                ogeom_topo::explore(model, f, ogeom_topo::Filter::OfType(ShapeType::Vertex))
                    .unwrap()
                    .iter()
                    .all(|v| {
                        let p = model.node(v).unwrap().data().as_vertex().unwrap().point;
                        let world = v.transform(model.datums()).unwrap().apply(p);
                        (world.x - 1.0).abs() < 1e-9
                    })
            })
            .expect("the box has a face at x = 1");

        // The second cube: the shared face seen from the other side, plus
        // five new faces over the shared boundary vertices.
        let corner = |x: f64, y: f64, z: f64| Point::new(x, y, z);
        let quad = |model: &mut Model, points: [Point; 4]| -> Shape {
            let wire = crate::make_polygon(model, &points, true, T).unwrap().shape;
            let plane = crate::find_plane(model, &wire, T).unwrap().unwrap();
            make_face(
                model,
                PlaneSurface::over(Plane::new(plane.frame()), (-4.0, 4.0), (-4.0, 4.0))
                    .unwrap()
                    .into(),
                std::slice::from_ref(&wire),
                T,
            )
            .unwrap()
            .shape
        };
        let faces = [
            quad(
                model,
                [
                    corner(2.0, 0.0, 0.0),
                    corner(2.0, 1.0, 0.0),
                    corner(2.0, 1.0, 1.0),
                    corner(2.0, 0.0, 1.0),
                ],
            ),
            quad(
                model,
                [
                    corner(1.0, 0.0, 0.0),
                    corner(2.0, 0.0, 0.0),
                    corner(2.0, 0.0, 1.0),
                    corner(1.0, 0.0, 1.0),
                ],
            ),
            quad(
                model,
                [
                    corner(1.0, 1.0, 0.0),
                    corner(1.0, 1.0, 1.0),
                    corner(2.0, 1.0, 1.0),
                    corner(2.0, 1.0, 0.0),
                ],
            ),
            quad(
                model,
                [
                    corner(1.0, 0.0, 0.0),
                    corner(1.0, 1.0, 0.0),
                    corner(2.0, 1.0, 0.0),
                    corner(2.0, 0.0, 0.0),
                ],
            ),
            quad(
                model,
                [
                    corner(1.0, 0.0, 1.0),
                    corner(2.0, 0.0, 1.0),
                    corner(2.0, 1.0, 1.0),
                    corner(1.0, 1.0, 1.0),
                ],
            ),
        ];
        let mut shell_faces = vec![shared.reversed()];
        shell_faces.extend(faces);
        let shell = make_shell(model, &shell_faces).unwrap().shape;
        let b = make_solid(model, std::slice::from_ref(&shell))
            .unwrap()
            .shape;
        (a, b)
    }

    #[test]
    fn glued_solids_build_a_compsolid_and_loose_ones_are_refused() {
        let mut model = Model::new();
        let (a, b) = glued_cubes(&mut model);
        let built = make_compsolid(&mut model, &[a.clone(), b.clone()]).unwrap();
        assert_eq!(model.kind_of(&built.shape).unwrap(), ShapeType::CompSolid);
        assert_eq!(
            built.history.generated(&a),
            std::slice::from_ref(&built.shape)
        );

        // Two boxes merely sitting apart share nothing and are refused.
        let mut model = Model::new();
        let a = crate::make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let far = Frame::new(
            Point::new(5.0, 0.0, 0.0),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
            T,
        )
        .unwrap();
        let b = crate::make_box(&mut model, far, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        assert!(make_compsolid(&mut model, &[a.clone(), b]).is_err());
        assert!(make_compsolid(&mut model, std::slice::from_ref(&a)).is_err());
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod band_tests {
    use super::*;
    use approx::assert_relative_eq;
    use ogeom_core::Tolerances;
    use ogeom_geom::CircleCurve;
    use ogeom_math::{Circle, Cone, Frame, Sphere};

    const T: Tolerances = Tolerances::millimetres();
    const TAU: f64 = core::f64::consts::TAU;

    fn mesh_area(model: &Model, face: &Shape) -> f64 {
        let fine = ogeom_mesh::Deflection {
            chord: 1e-3,
            ..ogeom_mesh::Deflection::default()
        };
        let mesh = ogeom_mesh::triangulate(model, face, fine, T).unwrap();
        mesh.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
                (b - a).cross(c - a).magnitude() / 2.0
            })
            .sum()
    }

    /// A full circle edge at `frame`'s origin, radius `r`, single vertex.
    fn ring(model: &mut Model, frame: Frame, r: f64) -> Shape {
        let circle = Circle::new(frame, r, T).unwrap();
        let vertex = make_vertex(model, frame.origin() + frame.x() * r).shape;
        make_edge_between(
            model,
            CircleCurve::new(circle).into(),
            (0.0, TAU),
            &vertex,
            &vertex,
            T,
        )
        .unwrap()
        .shape
    }

    #[test]
    fn a_cone_cap_takes_its_apex_as_a_degenerate_ring() {
        // Half angle 45 degrees, reference radius 2: apex at z = -2, slant
        // length 2*sqrt(2), lateral area pi * r * slant.
        let mut model = Model::new();
        let cone = Cone::new(Frame::WORLD, 2.0, core::f64::consts::FRAC_PI_4, T).unwrap();
        let surface: SurfaceGeometry = ogeom_geom::ConeSurface::new(cone, (-3.0, 3.0))
            .unwrap()
            .into();
        let rim = ring(&mut model, Frame::WORLD, 2.0);
        let face = crate::make_apex_band(&mut model, &surface, &rim, T).unwrap();
        assert_relative_eq!(
            mesh_area(&model, &face),
            core::f64::consts::PI * 2.0 * 2.0 * core::f64::consts::SQRT_2,
            max_relative = 1e-2
        );
    }

    #[test]
    fn a_degenerate_ring_is_accepted_whichever_side_it_is_passed_on() {
        let mut model = Model::new();
        let cone = Cone::new(Frame::WORLD, 2.0, core::f64::consts::FRAC_PI_4, T).unwrap();
        let surface: SurfaceGeometry = ogeom_geom::ConeSurface::new(cone, (-3.0, 3.0))
            .unwrap()
            .into();
        let rim = ring(&mut model, Frame::WORLD, 2.0);
        let apex = make_vertex(&mut model, cone.apex()).shape;
        let mut data = EdgeData::new();
        data.degenerate = true;
        let degenerate = model.add_edge(data, &[apex.clone(), apex]).unwrap();
        // Degenerate first: the band swaps it into place rather than asking
        // the caller to know which ring is real.
        let face = make_revolution_band(&mut model, &surface, &degenerate, &rim, T).unwrap();
        assert_relative_eq!(
            mesh_area(&model, &face),
            core::f64::consts::PI * 2.0 * 2.0 * core::f64::consts::SQRT_2,
            max_relative = 1e-2
        );
    }

    #[test]
    fn a_sphere_cap_closes_against_its_pole() {
        // The equator ring and the north pole: a hemisphere, area 2 pi r^2.
        let mut model = Model::new();
        let sphere = Sphere::new(Frame::WORLD, 3.0, T).unwrap();
        let surface: SurfaceGeometry = ogeom_geom::SphereSurface::new(sphere).into();
        let rim = ring(&mut model, Frame::WORLD, 3.0);
        let pole = make_vertex(&mut model, Point::new(0.0, 0.0, 3.0)).shape;
        let mut data = EdgeData::new();
        data.degenerate = true;
        let degenerate = model.add_edge(data, &[pole.clone(), pole]).unwrap();
        let face = make_revolution_band(&mut model, &surface, &rim, &degenerate, T).unwrap();
        assert_relative_eq!(
            mesh_area(&model, &face),
            2.0 * core::f64::consts::PI * 9.0,
            max_relative = 1e-2
        );
    }

    #[test]
    fn two_degenerate_rings_bound_nothing_and_are_refused() {
        let mut model = Model::new();
        let cone = Cone::new(Frame::WORLD, 2.0, core::f64::consts::FRAC_PI_4, T).unwrap();
        let surface: SurfaceGeometry = ogeom_geom::ConeSurface::new(cone, (-3.0, 3.0))
            .unwrap()
            .into();
        let apex = make_vertex(&mut model, cone.apex()).shape;
        let mut data = EdgeData::new();
        data.degenerate = true;
        let a = model
            .add_edge(data.clone(), &[apex.clone(), apex.clone()])
            .unwrap();
        let b = model.add_edge(data, &[apex.clone(), apex]).unwrap();
        assert!(make_revolution_band(&mut model, &surface, &a, &b, T).is_err());
    }

    #[test]
    fn an_apex_band_needs_a_cone() {
        let mut model = Model::new();
        let sphere = Sphere::new(Frame::WORLD, 3.0, T).unwrap();
        let surface: SurfaceGeometry = ogeom_geom::SphereSurface::new(sphere).into();
        let rim = ring(&mut model, Frame::WORLD, 3.0);
        assert!(crate::make_apex_band(&mut model, &surface, &rim, T).is_err());
    }
}
