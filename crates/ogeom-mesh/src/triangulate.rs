//! Triangulating a face.
//!
//! The second half of tessellation. A face is a trimmed region of a surface, so
//! the triangulation is built in the surface's `(u, v)` parameter space — where
//! the region is an ordinary polygon with holes — and then lifted back into
//! space by evaluating the surface at each vertex.
//!
//! # Why parameter space
//!
//! Triangulating in 3D would mean deciding which side of a curved boundary a
//! point falls on, in space, which is the point-in-solid problem. In parameter
//! space the boundary is a closed 2D polygon and the question is a winding
//! count. The surface does the rest.
//!
//! The cost is that parameter space is distorted: equal steps in `(u, v)` cover
//! very different distances near a sphere's pole than near its equator. So the
//! interior points are chosen by measuring deflection *in space* and the
//! triangulation is done in parameter space — measuring where the answer
//! matters, connecting where it is easy.
//!
//! # Watertightness
//!
//! A face's boundary points come from discretizing the *edge's* 3D curve and
//! evaluating the pcurve at those same parameters. Two faces sharing an edge
//! therefore place their boundary vertices at identical spatial positions, and
//! the join has no gap. Discretizing each face's pcurve independently would
//! give each face its own idea of where the edge runs, and the seams would show.

use ogeom_core::{Exact, OgeomResult, Predicates, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::{Curve2d, Surface, SurfaceGeometry};
use ogeom_math::{Direction, Point, Point2, Vector};
use ogeom_topo::{EdgeRepr, Model, NodeData, Orientation, Shape, ShapeType, Triangulation};
use spade::{ConstrainedDelaunayTriangulation, Point2 as SpadePoint, Triangulation as _};

use crate::discretize::{Deflection, discretize};

/// Triangulate one face.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `face` is not a
/// face or its geometry is missing;
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle fails to
/// resolve; [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the boundary
/// cannot be triangulated.
pub fn triangulate_face(
    model: &Model,
    face: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Triangulation> {
    deflection.validate()?;
    if model.kind_of(face)? != ShapeType::Face {
        ogeom_bail!(Construction, "expected a face");
    }
    let Some(node) = model.node(face) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let NodeData::Face(data) = node.data() else {
        ogeom_bail!(Construction, "face node holds no face data");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };
    let placement = face.transform(model.datums())?;

    let (uv, anchors, met) = trimming_rings(model, face, data.surface, surface, deflection, tol)?;
    let planar = triangulate_region(&uv, surface, deflection, tol)?;

    // Boundary vertices take their positions from their edges' own curves —
    // the shared authority — keyed by their exact parameter-space bits.
    let mut anchored: std::collections::HashMap<(u64, u64), Point> =
        std::collections::HashMap::new();
    for (ring, ring_anchor) in uv.iter().zip(&anchors) {
        for (p, a) in ring.iter().zip(ring_anchor) {
            if let Some(point) = a {
                anchored.insert((p.x.to_bits(), p.y.to_bits()), *point);
            }
        }
    }

    // Lift into space. The normal follows the face's orientation, not the
    // surface's: a reversed face presents the other side, and a renderer or a
    // volume computation that ignored that would have the solid inside out.
    let flip = face.orientation() == Orientation::Reversed;
    let mut mesh = Triangulation::new();
    mesh.deflection_met = met;
    let mut hits = 0_usize;
    for (u, v) in planar.parameters {
        // Anchors are already in world coordinates — their edges' own
        // placements applied — where surface lifts still need the face's.
        let point = match anchored.get(&(u.to_bits(), v.to_bits())) {
            Some(anchor) => {
                hits += 1;
                *anchor
            }
            None => placement.apply(surface.point_at(u, v, tol)?),
        };
        let normal = surface
            .normal_at(u, v, tol)
            .map_or(Vector::ZERO, |n| placement.apply_vector(n.vector()));
        mesh.positions.push(point);
        mesh.normals.push(if flip { -normal } else { normal });
        mesh.parameters.push((u, v));
    }
    if std::env::var("OGEOM_MESH_DEBUG").is_ok() {
        eprintln!(
            "DBG anchors map={} hits={} verts={}",
            anchored.len(),
            hits,
            mesh.positions.len()
        );
    }
    mesh.triangles = planar
        .triangles
        .into_iter()
        .map(|t| if flip { [t[0], t[2], t[1]] } else { t })
        .collect();
    Ok(mesh)
}

/// Triangulate every face below a shape, welded into one mesh.
///
/// # Errors
///
/// As [`triangulate_face`].
pub fn triangulate(
    model: &Model,
    shape: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Triangulation> {
    let mut mesh = Triangulation::new();
    for face in ogeom_topo::explore(model, shape, ogeom_topo::Filter::OfType(ShapeType::Face))? {
        mesh.append(&triangulate_face(model, &face, deflection, tol)?);
    }
    let mesh = mesh.welded(tol);

    // A second, border-only pass at the tolerance the model itself recorded:
    // imported edges carry the file's slop in their widened tolerances, and
    // two faces lifting the same edge through disagreeing geometry land that
    // far apart. Interior edges are already manifold and are not touched.
    let mut reach = 0.0_f64;
    for kind in [ShapeType::Edge, ShapeType::Vertex] {
        for shape in ogeom_topo::explore(model, shape, ogeom_topo::Filter::OfType(kind))? {
            let recorded = model.node(&shape).map_or(0.0, |n| match n.data() {
                NodeData::Edge(d) => d.tolerance.get(),
                NodeData::Vertex(d) => d.tolerance.get(),
                _ => 0.0,
            });
            reach = reach.max(recorded);
        }
    }
    // Floored at a tenth of the chord the caller asked for: a border gap
    // smaller than that is below the resolution of the mesh they accepted,
    // whether or not the model recorded the slop that caused it.
    let reach = reach.max(deflection.chord * 0.1);
    if reach > tol.confusion() {
        let reach = reach + tol.confusion();
        Ok(mesh.border_welded(reach).border_stitched(reach))
    } else {
        Ok(mesh)
    }
}

/// A face's trimming boundary, in its surface's parameter space.
///
/// The outer wire first, then any holes, each as a closed ring of `(u, v)`
/// points with no repeated closing point. A face with no wires covers its
/// surface's whole domain, and gets that rectangle as its boundary.
///
/// Public because parameter-space trimming is not only the triangulator's
/// concern: classifying a point against a face, splitting a face in a boolean,
/// and hidden-line removal all ask the same question of the same rings, and
/// each deriving them separately would be three chances to disagree.
///
/// # Errors
///
/// As [`triangulate_face`].
pub fn face_boundary(
    model: &Model,
    face: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Vec<Vec<Point2>>> {
    deflection.validate()?;
    if model.kind_of(face)? != ShapeType::Face {
        ogeom_bail!(Construction, "expected a face");
    }
    let Some(node) = model.node(face) else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let NodeData::Face(data) = node.data() else {
        ogeom_bail!(Construction, "face node holds no face data");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "face refers to a surface not in this model");
    };

    Ok(trimming_rings(model, face, data.surface, surface, deflection, tol)?.0)
}

/// The rings bounding a face in parameter space, and whether every boundary
/// edge met its deflection.
/// Boundary rings with, per ring vertex, the 3D anchor its edge's own curve
/// provides — `None` where an edge has no 3D curve to defer to.
type RingsWithAnchors = (Vec<Vec<Point2>>, Vec<Vec<Option<Point>>>, bool);

fn trimming_rings(
    model: &Model,
    face: &Shape,
    id: ogeom_topo::SurfaceId,
    surface: &SurfaceGeometry,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<RingsWithAnchors> {
    let mut rings = Vec::new();
    let mut ring_anchors = Vec::new();
    let mut met = true;
    for wire in model.ordered_children_of(face)? {
        let (ring, anchors, ring_met) = boundary_ring(model, &wire, id, deflection, tol)?;
        met &= ring_met;
        if ring.len() >= 3 {
            rings.push(ring);
            ring_anchors.push(anchors);
        }
    }
    if rings.is_empty() {
        // A face with no wires covers its surface's whole domain, so the domain
        // rectangle is the boundary.
        let ring = domain_ring(surface, deflection, tol);
        ring_anchors.push(vec![None; ring.len()]);
        rings.push(ring);
    }
    Ok((rings, ring_anchors, met))
}

/// Whether a point in parameter space lies inside the region `rings` bound.
///
/// Even-odd winding: inside the outer ring and outside every hole. The rings
/// come from wires whose direction already encodes outer from inner, but
/// counting crossings does not depend on that being right, which makes it
/// robust to a wire that was built the wrong way round.
///
/// Says nothing about a point *on* a ring — the crossing count of a boundary
/// point is whichever side rounding puts it. A caller that cares has to measure
/// its distance to the boundary and decide, which is what classification does.
///
/// Decided with exact predicates. See [`inside_boundary_with`].
#[must_use]
pub fn inside_boundary(rings: &[Vec<Point2>], p: Point2) -> bool {
    inside_boundary_with::<Exact>(rings, p)
}

/// As [`inside_boundary`], with the predicate implementation named.
///
/// This is the seam `docs/DATA_MODEL.md` §9 describes, and it is here rather
/// than anywhere else because this is where the *combinatorial* decision is.
/// Whether a point is inside a ring is not a measurement that can be a little
/// wrong: it decides whether a triangle is kept or dropped, so an error near a
/// boundary is a hole in the mesh rather than a slightly misplaced one.
///
/// The naive form divides to find where an edge crosses the sampling ray, and
/// that division cancels catastrophically for a point nearly on the edge.
/// `orient2d` answers the same question with no division at all, and
/// [`Exact`] answers it correctly however close the point is.
#[must_use]
pub fn inside_boundary_with<P: Predicates>(rings: &[Vec<Point2>], p: Point2) -> bool {
    let mut inside = false;
    for ring in rings {
        if crosses_odd_times::<P>(ring, p) {
            inside = !inside;
        }
    }
    inside
}

/// The boundary of one wire, in the face's parameter space.
///
/// Each edge is discretized in *space* and its pcurve evaluated at the resulting
/// parameters, so two faces sharing the edge agree on where its points are.
fn boundary_ring(
    model: &Model,
    wire: &Shape,
    surface: ogeom_topo::SurfaceId,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<(Vec<Point2>, Vec<Option<Point>>, bool)> {
    let mut ring: Vec<Point2> = Vec::new();
    let mut anchors: Vec<Option<Point>> = Vec::new();
    let mut met = true;

    // Start the walk off a seam if the wire allows it: a seam's side is
    // chosen by continuity with the point already walked to, and continuity
    // needs something to continue from. The ring is cyclic, so rotating the
    // walk changes nothing it reports.
    let mut children = model.ordered_children_of(wire)?;
    let is_seam = |model: &Model, e: &Shape| -> bool {
        model
            .node(e)
            .and_then(|n| n.data().as_edge())
            .and_then(|d| d.pcurve_for(surface, e.location()))
            .is_some_and(|r| matches!(r, EdgeRepr::Seam { .. }))
    };
    if let Some(start) = children.iter().position(|e| !is_seam(model, e)) {
        children.rotate_left(start);
    }

    for edge in children {
        let Some(node) = model.node(&edge) else {
            ogeom_bail!(Dangling, "edge is not in this model");
        };
        let NodeData::Edge(data) = node.data() else {
            ogeom_bail!(Construction, "edge node holds no edge data");
        };
        let (pcurve_id, pcurve_range) = match data.pcurve_for(surface, edge.location()) {
            Some(EdgeRepr::PCurve { curve, range, .. }) => (*curve, *range),
            // A seam edge runs along a closed surface's join and bounds its
            // face twice — up one side of the parameter rectangle and down
            // the other. Which side this occurrence takes is decided by the
            // ring itself: the side whose oriented start continues the point
            // already walked to. Orientation flags cannot answer it — a
            // reversed face flips every occurrence while the chart columns
            // stay where they were built — but the chart can.
            Some(EdgeRepr::Seam {
                forward,
                reversed,
                range,
                ..
            }) => {
                let (f, r) = (*forward, *reversed);
                let picked = if let Some(last) = ring.last().copied() {
                    let start_of = |id: ogeom_topo::PCurveId| -> Option<Point2> {
                        let pc = model.geometry().pcurve(id)?;
                        let t = if edge.orientation() == Orientation::Reversed {
                            range.1
                        } else {
                            range.0
                        };
                        pc.point_at(t, tol).ok()
                    };
                    match (start_of(f), start_of(r)) {
                        (Some(a), Some(b)) => {
                            if last.distance(a) <= last.distance(b) {
                                f
                            } else {
                                r
                            }
                        }
                        _ => f,
                    }
                } else if edge.orientation() == Orientation::Reversed {
                    r
                } else {
                    f
                };
                (picked, *range)
            }
            _ => ogeom_bail!(
                Construction,
                "edge has no pcurve on this face, so the face cannot be \
                 triangulated in its own parameter space"
            ),
        };
        let Some(pcurve) = model.geometry().pcurve(pcurve_id) else {
            ogeom_bail!(Dangling, "pcurve is not in this model");
        };

        // Sample where the *3D* curve says to, so an adjacent face lands on
        // the same points — and *anchor* the boundary vertices to that curve
        // too: two faces sharing an edge lift the same parameters through
        // different surfaces, and on an imported file those surfaces
        // disagree by the file's own slop. The edge is the shared authority,
        // so its points are the positions both faces use, and the weld is a
        // matter of identity rather than luck.
        let mut edge_anchors: Vec<Option<Point>> = Vec::new();
        let samples = match sample_parameters(model, data, deflection, tol)? {
            Some((parameters, edge_met)) => {
                met &= edge_met;
                if let Some(EdgeRepr::Curve3d { curve, .. }) = data.curve3d()
                    && let Some(geometry) = model.geometry().curve(*curve)
                    && let Ok(edge_placement) = edge.transform(model.datums())
                {
                    for t in &parameters {
                        edge_anchors.push(
                            geometry
                                .point_at(*t, tol)
                                .ok()
                                .map(|p| edge_placement.apply(p)),
                        );
                    }
                }
                map_to_pcurve(&parameters, data, pcurve_range)
            }
            // No 3D curve to defer to — the pcurve's own shape, measured in
            // space through the surface, so the chord tolerance means the
            // same thing it means everywhere else.
            None => {
                let Some(geometry) = model.geometry().surface(surface) else {
                    ogeom_bail!(Dangling, "face refers to a surface not in this model");
                };
                let (_, parameters) = crate::discretize::discretize_on_surface(
                    pcurve,
                    pcurve_range,
                    geometry,
                    deflection,
                    tol,
                )?;
                parameters
            }
        };

        let mut points: Vec<Point2> = samples
            .iter()
            .map(|u| pcurve.point_at(*u, tol))
            .collect::<OgeomResult<_>>()?;
        if edge.orientation() == Orientation::Reversed {
            points.reverse();
            edge_anchors.reverse();
        }
        if edge_anchors.len() != points.len() {
            edge_anchors = vec![None; points.len()];
        }
        // The ends of an edge belong to its *vertices* — the one authority
        // every face and every neighbouring edge shares — but only within the
        // tolerance the vertex itself records. An imported curve ends within
        // the vertex's widened tolerance of it, and lifting the ends through
        // the curve alone would leave each corner split as many ways as there
        // are curves meeting there; a vertex that sits *beyond* its stated
        // tolerance from the curve is not describing the curve's end at all,
        // and the curve stays the authority.
        if !points.is_empty()
            && let Ok(vs) = model.children_of(&edge)
            && vs.len() >= 2
            && let Ok(edge_placement) = edge.transform(model.datums())
        {
            let point_of = |v: &Shape| -> Option<(Point, f64)> {
                let data = model.node(v)?.data().as_vertex()?;
                Some((edge_placement.apply(data.point), data.tolerance.get()))
            };
            let (from, to) = if edge.orientation() == Orientation::Reversed {
                (&vs[vs.len() - 1], &vs[0])
            } else {
                (&vs[0], &vs[vs.len() - 1])
            };
            if let Some((p, within)) = point_of(from)
                && let Some(a) = edge_anchors.first_mut()
                && a.is_none_or(|end| end.distance(p) <= within + tol.confusion())
            {
                *a = Some(p);
            }
            if let Some((p, within)) = point_of(to)
                && let Some(a) = edge_anchors.last_mut()
                && a.is_none_or(|end| end.distance(p) <= within + tol.confusion())
            {
                *a = Some(p);
            }
        }
        // Fold onto the branch that continues the ring. Two faces sharing an
        // edge share its pcurve, and on a periodic surface the pcurve sits in
        // *one* face's window: a cylinder split into two halves has a ruling
        // at u = 0 that the other half needs at u = 2pi. The chart cannot
        // store both; continuity with the ring being walked recovers the
        // right branch, exactly as the seam sides are chosen.
        if let Some(last) = ring.last().copied()
            && let Some(first) = points.first().copied()
            && let Some(geometry) = model.geometry().surface(surface)
        {
            use ogeom_geom::Surface as _;
            let ((ua, ub), (va, vb)) = geometry.domain();
            let mut shift = Point2::new(0.0, 0.0);
            if geometry.is_periodic_u() && (ub - ua) > 0.0 {
                shift.x = ((last.x - first.x) / (ub - ua)).round() * (ub - ua);
            }
            if geometry.is_periodic_v() && (vb - va) > 0.0 {
                shift.y = ((last.y - first.y) / (vb - va)).round() * (vb - va);
            }
            if shift.x != 0.0 || shift.y != 0.0 {
                for p in &mut points {
                    p.x += shift.x;
                    p.y += shift.y;
                }
            }
        }
        // The previous edge already contributed the shared vertex.
        if !ring.is_empty() && !points.is_empty() {
            points.remove(0);
            edge_anchors.remove(0);
        }
        ring.extend(points);
        anchors.extend(edge_anchors);
    }

    // A closed ring repeats its first point at the end; the triangulator wants
    // it named once.
    if ring.len() > 2
        && let (Some(first), Some(last)) = (ring.first().copied(), ring.last().copied())
        && first.is_equal(last, tol)
    {
        ring.pop();
        anchors.pop();
    }
    Ok((ring, anchors, met))
}

/// Parameters at which to sample an edge, taken from its 3D curve.
fn sample_parameters(
    model: &Model,
    data: &ogeom_topo::EdgeData,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Option<(Vec<f64>, bool)>> {
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        return Ok(None);
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    let line = discretize(geometry, *range, deflection, tol)?;
    Ok(Some((line.parameters, line.deflection_met)))
}

/// Map parameters on the 3D curve onto the pcurve's own range.
///
/// The two are parameterized over their own intervals; `same_parameter` means
/// they agree *proportionally*, which is what this converts.
fn map_to_pcurve(
    parameters: &[f64],
    data: &ogeom_topo::EdgeData,
    pcurve_range: (f64, f64),
) -> Vec<f64> {
    let Some(EdgeRepr::Curve3d { range, .. }) = data.curve3d() else {
        return parameters.to_vec();
    };
    let (ca, cb) = *range;
    let (pa, pb) = pcurve_range;
    if (cb - ca).abs() <= f64::MIN_POSITIVE {
        return parameters.to_vec();
    }
    parameters
        .iter()
        .map(|u| pa + (pb - pa) * (u - ca) / (cb - ca))
        .collect()
}

/// The edge of a surface's domain, as a boundary ring.
///
/// Refined, not just the four corners. A triangulation only ever connects the
/// points it is given, so a boundary named by its corners alone forces long
/// triangles reaching right across the domain to find one — on a sphere, a
/// sliver from the equator to the pole. The interior refinement cannot fix that;
/// the missing points are on the boundary.
fn domain_ring(surface: &SurfaceGeometry, deflection: Deflection, tol: Tolerances) -> Vec<Point2> {
    let ((ua, ub), (va, vb)) = surface.domain();
    if ![ua, ub, va, vb].iter().all(|x| x.is_finite()) {
        return Vec::new();
    }

    let along_u = |v: f64| {
        refine_direction(ua, ub, deflection.chord, |a, b| {
            sag_between(surface, (a, v), (b, v), tol)
        })
    };
    let along_v = |u: f64| {
        refine_direction(va, vb, deflection.chord, |a, b| {
            sag_between(surface, (u, a), (u, b), tol)
        })
    };

    // Counter-clockwise around the rectangle. Each side drops its final point,
    // which the next side contributes: a repeated vertex would be a
    // zero-length boundary edge, and a constraint of zero length is not one.
    let (bottom, top) = (along_u(va), along_u(vb));
    let (left, right) = (along_v(ua), along_v(ub));
    let mut ring = Vec::new();
    ring.extend(
        bottom[..bottom.len() - 1]
            .iter()
            .map(|u| Point2::new(*u, va)),
    );
    ring.extend(right[..right.len() - 1].iter().map(|v| Point2::new(ub, *v)));
    ring.extend(top[1..].iter().rev().map(|u| Point2::new(*u, vb)));
    ring.extend(left[1..].iter().rev().map(|v| Point2::new(ua, *v)));
    ring
}

/// A triangulation still in parameter space, before it is lifted onto the
/// surface.
struct PlanarMesh {
    /// The `(u, v)` of each vertex.
    parameters: Vec<(f64, f64)>,
    /// Triangles as indices into `parameters`.
    triangles: Vec<[u32; 3]>,
}

/// Triangulate a region in parameter space, given its boundary rings.
///
/// The first ring is the outer boundary; the rest are holes.
fn triangulate_region(
    rings: &[Vec<Point2>],
    surface: &SurfaceGeometry,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<PlanarMesh> {
    let mut cdt: ConstrainedDelaunayTriangulation<SpadePoint<f64>> =
        ConstrainedDelaunayTriangulation::new();

    // The boundary edges are constraints, so the triangulation respects the
    // trimming rather than spanning across a hole.
    for ring in rings {
        let mut ring_handles = Vec::with_capacity(ring.len());
        for p in ring {
            let handle = cdt
                .insert(SpadePoint::new(p.x, p.y))
                .map_err(|e| ogeom_core::ogeom_err!(NotDone, "boundary insertion failed: {e}"))?;
            ring_handles.push(handle);
        }
        for i in 0..ring_handles.len() {
            let (a, b) = (ring_handles[i], ring_handles[(i + 1) % ring_handles.len()]);
            if a != b && cdt.can_add_constraint(a, b) {
                cdt.add_constraint(a, b);
            }
        }
    }

    // Interior points where the surface bends away from the flat triangle. A
    // planar face needs none, which is why this is driven by measured
    // deflection rather than by a fixed grid.
    add_interior_points(&mut cdt, rings, surface, deflection, tol)?;

    // The grid rows guarantee the deflection along their own lines, but a
    // hole in the face punches a gap through a row, and where the surface is
    // flat in one direction — a cylinder along its axis — there may be no
    // other row for the mesher to reach. The band around the hole then fans
    // from the rim to the far side of the gap, in triangles that sag through
    // the solid by far more than the deflection while every one of their
    // vertices sits exactly on the surface. The boolean caught this as a
    // fused solid whose faces all had the right area and the wrong volume.
    //
    // The repair measures the truth: any kept triangle whose midpoints sag
    // beyond the chord gets its centre inserted, and the loop runs until the
    // mesh is honest or the cap says the surface is being unreasonable.
    if !matches!(surface.kind(), ogeom_geom::SurfaceKind::Plane) {
        for _ in 0..REFINEMENT_ROUNDS {
            let mut worst: Vec<SpadePoint<f64>> = Vec::new();
            for triangle in cdt.inner_faces() {
                let vertices = triangle.vertices();
                let centre = triangle.center();
                let at = Point2::new(centre.x, centre.y);
                if !inside_region(rings, at) {
                    continue;
                }
                let corners: Vec<(f64, f64)> = vertices
                    .iter()
                    .map(|v| {
                        let p = v.position();
                        (p.x, p.y)
                    })
                    .collect();
                // The grid already bounds sag along rows and columns, and a
                // grid triangle's diagonal spanning one cell each way may
                // legitimately sag up to the sum — twice the chord — which
                // was the guarantee before this loop existed. The threshold
                // sits clear above that band so the repair fires only on the
                // fan triangles it exists for, which sag through a hole's
                // gap by tens of chords, and an honest grid — including a
                // perfectly symmetric one, whose mesh must stay symmetric —
                // is left untouched.
                let sagged = (0..3).any(|i| {
                    sag_between(surface, corners[i], corners[(i + 1) % 3], tol)
                        > deflection.chord * 3.0
                });
                if sagged {
                    worst.push(SpadePoint::new(centre.x, centre.y));
                }
            }
            if worst.is_empty() {
                break;
            }
            for point in worst {
                cdt.insert(point).map_err(|e| {
                    ogeom_core::ogeom_err!(NotDone, "refinement insertion failed: {e}")
                })?;
            }
        }
    }

    let mut parameters = Vec::new();
    let mut index_of = std::collections::HashMap::new();
    for (i, vertex) in cdt.vertices().enumerate() {
        let p = vertex.position();
        index_of.insert(vertex.fix(), i);
        parameters.push((p.x, p.y));
    }

    let mut triangles = Vec::new();
    for triangle in cdt.inner_faces() {
        let vertices = triangle.vertices();
        let centre = triangle.center();
        // A constrained Delaunay covers the convex hull of its input, so
        // triangles outside the trimmed region — across a concavity, or inside
        // a hole — have to be discarded. Winding tells them apart.
        if !inside_region(rings, Point2::new(centre.x, centre.y)) {
            continue;
        }
        #[allow(clippy::cast_possible_truncation)]
        let indices = [
            index_of[&vertices[0].fix()] as u32,
            index_of[&vertices[1].fix()] as u32,
            index_of[&vertices[2].fix()] as u32,
        ];
        triangles.push(indices);
    }

    if triangles.is_empty() {
        ogeom_bail!(
            NotDone,
            "the face's boundary enclosed no triangulable region"
        );
    }
    Ok(PlanarMesh {
        parameters,
        triangles,
    })
}

/// Add interior points wherever the surface deviates from flat by more than the
/// deflection allows.
///
/// The two parameter directions are refined *independently*, and that is not an
/// optimization — it is the difference between converging and not. A uniform
/// grid on a sphere puts as many meridians through the pole as through the
/// equator, so the triangles there become arbitrarily thin slivers in space.
/// The summed area of such a mesh does not approach the sphere's; it grows
/// without bound as the grid tightens. (Schwarz's lantern is the standard
/// example: an inscribed polyhedron whose area diverges under refinement.)
///
/// Refining each direction by its own measured sag fixes it at the source. Near
/// a pole the circle of latitude has almost no radius, so a chord right across
/// it sags by almost nothing and the direction stops subdividing after one or
/// two steps — while the meridian direction, whose curvature does not change,
/// keeps refining. The mesh degenerates into a fan, which is the right shape.
fn add_interior_points(
    cdt: &mut ConstrainedDelaunayTriangulation<SpadePoint<f64>>,
    rings: &[Vec<Point2>],
    surface: &SurfaceGeometry,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<()> {
    // A plane is flat everywhere; sampling it would add points that buy nothing.
    if matches!(surface.kind(), ogeom_geom::SurfaceKind::Plane) {
        return Ok(());
    }

    let bound = rings
        .iter()
        .flatten()
        .fold(ogeom_math::Aabb::EMPTY, |acc, p| {
            acc.with_point(Point::new(p.x, p.y, 0.0))
        });
    let (Some(low), Some(high)) = (bound.low(), bound.high()) else {
        return Ok(());
    };

    // The v resolution has to hold everywhere the region reaches, so its sag is
    // the worst over a spread of u probes rather than the sag along one line.
    // A surface of revolution is the same at every u and a lofted one is not.
    #[allow(clippy::cast_precision_loss)]
    let probes: Vec<f64> = (0..=U_PROBES)
        .map(|i| low.x + (high.x - low.x) * i as f64 / U_PROBES as f64)
        .collect();
    let rows = refine_direction(low.y, high.y, deflection.chord, |a, b| {
        probes
            .iter()
            .map(|u| sag_between(surface, (*u, a), (*u, b), tol))
            .fold(0.0_f64, f64::max)
    });

    for &v in &rows[1..rows.len().saturating_sub(1)] {
        // Each row gets its own u resolution, measured at that row.
        let columns = refine_direction(low.x, high.x, deflection.chord, |a, b| {
            sag_between(surface, (a, v), (b, v), tol)
        });
        for &u in &columns[1..columns.len().saturating_sub(1)] {
            // Interior points only: the boundary is already constrained, and a
            // point landing just off a constraint would split it.
            if !inside_region(rings, Point2::new(u, v)) {
                continue;
            }
            cdt.insert(SpadePoint::new(u, v))
                .map_err(|e| ogeom_core::ogeom_err!(NotDone, "interior insertion failed: {e}"))?;
        }
    }
    Ok(())
}

/// How many places across the domain the v resolution is measured at.
const U_PROBES: usize = 8;

/// How many rounds of sag-driven refinement a region may take.
///
/// Each round halves the worst sag roughly; a surface not honest after this
/// many is degenerate, and the cap makes that a coarse mesh rather than an
/// exhausted allocator.
const REFINEMENT_ROUNDS: usize = 12;

/// The most subdivisions one parameter direction may take.
///
/// A surface that has not converged by 512 has a singularity, not a resolution
/// problem, and the ceiling is what makes that a coarse mesh rather than an
/// exhausted allocator.
const MAX_DIRECTION_STEPS: usize = 512;

/// Subdivide `[lo, hi]` until no sub-interval sags further than `chord`.
///
/// The same adaptive bisection [`discretize`] uses on a curve, applied to a
/// line through parameter space. Returns the parameters in increasing order,
/// endpoints included.
fn refine_direction<F: Fn(f64, f64) -> f64>(lo: f64, hi: f64, chord: f64, sag: F) -> Vec<f64> {
    let mut values = vec![lo, f64::midpoint(lo, hi), hi];
    while values.len() < MAX_DIRECTION_STEPS {
        let Some(i) = (0..values.len() - 1).find(|&i| sag(values[i], values[i + 1]) > chord) else {
            break;
        };
        let mid = f64::midpoint(values[i], values[i + 1]);
        // A split that does not divide the interval means the parameters have
        // reached the resolution of f64, and refining further would loop.
        if mid <= values[i] || mid >= values[i + 1] {
            break;
        }
        values.insert(i + 1, mid);
    }
    values
}

/// How far the surface departs from the chord joining two parameter points.
///
/// Measured in space, which is the only place the number means anything: the
/// same step in `u` covers a metre at a sphere's equator and a millimetre near
/// its pole.
fn sag_between(
    surface: &SurfaceGeometry,
    from: (f64, f64),
    to: (f64, f64),
    tol: Tolerances,
) -> f64 {
    let mid = (f64::midpoint(from.0, to.0), f64::midpoint(from.1, to.1));
    let (Ok(a), Ok(b), Ok(m)) = (
        surface.point_at(from.0, from.1, tol),
        surface.point_at(to.0, to.1, tol),
        surface.point_at(mid.0, mid.1, tol),
    ) else {
        // Off the surface's domain; nothing to refine towards.
        return 0.0;
    };
    ogeom_math::Axis::through(a, b, tol).map_or_else(|_| a.distance(m), |axis| axis.distance_to(m))
}

/// Whether a point is inside the region the rings bound.
fn inside_region(rings: &[Vec<Point2>], p: Point2) -> bool {
    inside_boundary_with::<Exact>(rings, p)
}

/// Ray-crossing count for one ring, decided by orientation rather than by
/// arithmetic.
///
/// A horizontal ray in `+u`. The half-open `y` comparison counts a vertex lying
/// exactly on the ray once rather than twice or not at all; which side of the
/// edge the point falls on is then an `orient2d` sign.
///
/// Deliberately *not* "solve for where the edge crosses the ray, then compare".
/// That form divides by the edge's `y` extent, which is near zero for a nearly
/// horizontal edge, and subtracts two nearly equal numbers to compare — so for a
/// point close to the boundary it can answer either way. Here there is no
/// division and the comparison is a determinant's sign, which an exact predicate
/// gets right at any separation.
fn crosses_odd_times<P: Predicates>(ring: &[Point2], p: Point2) -> bool {
    let mut inside = false;
    let n = ring.len();
    let at = |q: Point2| [q.x, q.y];
    for i in 0..n {
        let (a, b) = (ring[i], ring[(i + 1) % n]);
        if (a.y > p.y) == (b.y > p.y) {
            continue;
        }
        // The edge crosses the ray's line. Whether it crosses the ray *itself*
        // — to the right of `p` — is which side of the directed edge `p` is on,
        // read the right way round for the edge's direction in `y`.
        let side = P::orient2d(at(a), at(b), [p.x, p.y]);
        let rightwards = if b.y > a.y {
            side == ogeom_core::Sign::Positive
        } else {
            side == ogeom_core::Sign::Negative
        };
        if rightwards {
            inside = !inside;
        }
    }
    inside
}

/// The unit normal of a triangle, or `None` if it is degenerate.
#[must_use]
pub fn triangle_normal(a: Point, b: Point, c: Point, tol: Tolerances) -> Option<Direction> {
    Direction::from_cross(b - a, c - a, tol).ok()
}

/// Discretize an edge into a polyline in space, for display or coarse queries.
///
/// # Errors
///
/// As [`discretize`].
pub fn polyline_of_edge(
    model: &Model,
    edge: &Shape,
    deflection: Deflection,
    tol: Tolerances,
) -> OgeomResult<Vec<Point>> {
    if model.kind_of(edge)? != ShapeType::Edge {
        ogeom_bail!(Construction, "expected an edge");
    }
    let Some(node) = model.node(edge) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let NodeData::Edge(data) = node.data() else {
        ogeom_bail!(Construction, "edge node holds no edge data");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        return Ok(Vec::new());
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    let placement = edge.transform(model.datums())?;
    let line = discretize(geometry, *range, deflection, tol)?;
    let mut points: Vec<Point> = line.points.iter().map(|p| placement.apply(*p)).collect();
    if edge.orientation() == Orientation::Reversed {
        points.reverse();
    }
    Ok(points)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ogeom_algo::make_box;
    use ogeom_math::Frame;
    use ogeom_topo::explore_unique;

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> Deflection {
        Deflection {
            chord: 1e-3,
            angular: 0.05,
            ..Deflection::default()
        }
    }

    #[test]
    fn a_box_face_triangulates_into_two_triangles() {
        // A planar square needs no interior points at all, which is what makes
        // deflection-driven refinement worth having: a fixed grid would add
        // dozens that buy nothing.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T).unwrap();
        let faces = explore_unique(&model, &built.shape, ShapeType::Face).unwrap();

        for face in &faces {
            let mesh = triangulate_face(&model, face, fine(), T).unwrap();
            assert_eq!(mesh.triangle_count(), 2, "a rectangle is two triangles");
            assert_eq!(mesh.vertex_count(), 4);
            assert!(mesh.deflection_met);
        }
    }

    #[test]
    fn a_boxs_mesh_is_closed_and_reports_the_right_volume() {
        // The end-to-end check: triangulate, weld, and ask the mesh what it
        // encloses. A volume that comes out negative would mean the faces are
        // wound inward; one that is wrong in magnitude would mean the
        // triangulation is not covering the boundary.
        let mut model = Model::new();
        let size = (2.0, 3.0, 4.0);
        let built = make_box(&mut model, Frame::WORLD, size, T).unwrap();
        let mesh = triangulate(&model, &built.shape, fine(), T).unwrap();

        assert_eq!(mesh.triangle_count(), 12, "six faces, two triangles each");
        assert_eq!(mesh.vertex_count(), 8, "welding merged the shared corners");
        assert!(
            mesh.is_closed(),
            "every triangle edge should be shared by two"
        );
        assert_relative_eq!(mesh.volume(), size.0 * size.1 * size.2, epsilon = 1e-9);
        assert_relative_eq!(
            mesh.area(),
            2.0 * (size.0 * size.1 + size.1 * size.2 + size.2 * size.0),
            epsilon = 1e-9
        );
    }

    #[test]
    fn welding_is_what_closes_the_mesh() {
        // Without it each face brings its own copy of every boundary vertex, so
        // no triangle edge is shared and the surface is a pile of loose squares.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();

        let mut loose = Triangulation::new();
        for face in ogeom_topo::explore(
            &model,
            &built.shape,
            ogeom_topo::Filter::OfType(ShapeType::Face),
        )
        .unwrap()
        {
            loose.append(&triangulate_face(&model, &face, fine(), T).unwrap());
        }
        assert_eq!(loose.vertex_count(), 24, "four corners per face, unmerged");
        assert!(!loose.is_closed());

        let welded = loose.welded(T);
        assert_eq!(welded.vertex_count(), 8);
        assert!(welded.is_closed());
    }

    #[test]
    fn a_reversed_face_presents_the_other_side() {
        // A renderer or a volume computation that ignored orientation would
        // have the solid inside out, and nothing about the positions says so.
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let face = explore_unique(&model, &built.shape, ShapeType::Face).unwrap()[0].clone();

        let forward = triangulate_face(&model, &face, fine(), T).unwrap();
        let backward = triangulate_face(&model, &face.reversed(), fine(), T).unwrap();

        assert_eq!(forward.triangle_count(), backward.triangle_count());
        for (a, b) in forward.normals.iter().zip(&backward.normals) {
            assert!(a.is_equal(-*b, T), "normals did not flip: {a:?} vs {b:?}");
        }
        // And the winding flipped with them, so the two agree.
        let winding = |m: &Triangulation, i: usize| {
            let [a, b, c] = m.triangles[i].map(|k| m.positions[k as usize]);
            (b - a).cross(c - a)
        };
        assert!(winding(&forward, 0).dot(winding(&backward, 0)) < 0.0);
    }

    #[test]
    fn every_triangle_vertex_lies_on_the_surface_it_came_from() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (2.0, 1.0, 3.0), T).unwrap();

        for face in explore_unique(&model, &built.shape, ShapeType::Face).unwrap() {
            let mesh = triangulate_face(&model, &face, fine(), T).unwrap();
            let data = model.node(&face).unwrap().data().as_face().unwrap().clone();
            let surface = model.geometry().surface(data.surface).unwrap();
            for (position, (u, v)) in mesh.positions.iter().zip(&mesh.parameters) {
                let exact = surface.point_at(*u, *v, T).unwrap();
                assert!(
                    position.is_equal(exact, T),
                    "{position:?} is not on its surface"
                );
            }
        }
    }

    #[test]
    fn a_finer_deflection_never_gives_fewer_triangles() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let coarse = triangulate(&model, &built.shape, Deflection::default(), T).unwrap();
        let detailed = triangulate(&model, &built.shape, fine(), T).unwrap();
        assert!(detailed.triangle_count() >= coarse.triangle_count());
        // Both enclose the same volume, since the faces are flat.
        assert_relative_eq!(coarse.volume(), 1.0, epsilon = 1e-9);
        assert_relative_eq!(detailed.volume(), 1.0, epsilon = 1e-9);
    }

    #[test]
    fn a_sphere_gets_interior_points_and_converges_on_its_true_area() {
        // The whole point of measuring deflection in space rather than in
        // parameter space. A sphere's parameter rectangle is uniform; the
        // surface it maps to is not, and a fixed grid would be dense at the
        // poles and coarse at the equator.
        use ogeom_algo::make_natural_face;
        use ogeom_geom::SphereSurface;
        use ogeom_math::Sphere;

        let radius = 10.0;
        let exact = 4.0 * std::f64::consts::PI * radius * radius;
        let mut previous = 0.0;

        for chord in [1.0_f64, 0.25, 0.05] {
            let mut model = Model::new();
            let surface = SphereSurface::new(Sphere::new(Frame::WORLD, radius, T).unwrap());
            let face = make_natural_face(&mut model, surface.into()).unwrap().shape;
            let deflection = Deflection {
                chord,
                ..Deflection::default()
            };
            let mesh = triangulate_face(&model, &face, deflection, T).unwrap();

            assert!(
                mesh.triangle_count() > 2,
                "a curved face needs interior points, got {} triangles",
                mesh.triangle_count()
            );
            // Every triangle chord-cuts the sphere, so the area comes in under
            // the truth and climbs as the tolerance tightens.
            let area = mesh.area();
            assert!(area < exact, "a chord-cut area cannot exceed the surface's");
            assert!(
                area > previous,
                "tightening the chord from the previous step lost area: \
                 {area} after {previous}"
            );
            previous = area;
        }
        assert!(
            previous > exact * 0.99,
            "at a chord of 0.05 on a radius of 10 the area should be within a \
             percent, got {previous} against {exact}"
        );
    }

    #[test]
    fn every_sphere_vertex_sits_at_the_right_radius() {
        // Lifting through the surface is what makes the mesh curved at all; a
        // vertex left in parameter space, or lifted with the wrong parameters,
        // would land nowhere near the sphere.
        use ogeom_algo::make_natural_face;
        use ogeom_geom::SphereSurface;
        use ogeom_math::Sphere;

        let mut model = Model::new();
        let surface = SphereSurface::new(Sphere::new(Frame::WORLD, 3.0, T).unwrap());
        let face = make_natural_face(&mut model, surface.into()).unwrap().shape;
        let mesh = triangulate_face(&model, &face, Deflection::default(), T).unwrap();

        for p in &mesh.positions {
            assert_relative_eq!(p.to_vector().magnitude(), 3.0, epsilon = 1e-9);
        }
    }

    #[test]
    fn a_curved_domains_boundary_is_refined_not_just_its_corners() {
        // The corners alone would leave the triangulation nothing to connect to
        // along a side, so it reaches right across the domain for one — and a
        // sliver from a sphere's equator to its pole makes the summed area
        // diverge under refinement rather than converge.
        use ogeom_geom::{PlaneSurface, SphereSurface};
        use ogeom_math::{Plane, Sphere};

        let sphere: SurfaceGeometry =
            SphereSurface::new(Sphere::new(Frame::WORLD, 10.0, T).unwrap()).into();
        let coarse = domain_ring(&sphere, Deflection::default(), T);
        let fine_ring = domain_ring(&sphere, fine(), T);
        assert!(coarse.len() > 4, "a sphere's domain edge is curved");
        assert!(
            fine_ring.len() > coarse.len(),
            "a tighter chord should place more boundary points"
        );

        // No side may repeat a corner: a zero-length constraint is not one.
        for w in fine_ring.windows(2) {
            assert!(!w[0].is_equal(w[1], T), "the ring repeats a point");
        }

        // A plane is flat, so its domain needs only what closes the rectangle.
        let plane: SurfaceGeometry = PlaneSurface::new(Plane::new(Frame::WORLD)).into();
        assert!(domain_ring(&plane, fine(), T).len() >= 4);
    }

    #[test]
    fn winding_detects_points_inside_and_outside_a_ring() {
        let square = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ];
        assert!(inside_region(
            std::slice::from_ref(&square),
            Point2::new(0.5, 0.5)
        ));
        assert!(!inside_region(
            std::slice::from_ref(&square),
            Point2::new(1.5, 0.5)
        ));
        assert!(!inside_region(
            std::slice::from_ref(&square),
            Point2::new(0.5, -0.5)
        ));

        // With a hole, the middle is outside again.
        let hole = vec![
            Point2::new(0.4, 0.4),
            Point2::new(0.6, 0.4),
            Point2::new(0.6, 0.6),
            Point2::new(0.4, 0.6),
        ];
        let with_hole = vec![square, hole];
        assert!(!inside_region(&with_hole, Point2::new(0.5, 0.5)));
        assert!(inside_region(&with_hole, Point2::new(0.2, 0.2)));
    }

    #[test]
    fn a_polyline_of_an_edge_runs_in_the_edges_direction() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        let edge = explore_unique(&model, &built.shape, ShapeType::Edge).unwrap()[0].clone();

        let forward = polyline_of_edge(&model, &edge, fine(), T).unwrap();
        let backward = polyline_of_edge(&model, &edge.reversed(), fine(), T).unwrap();
        assert!(forward.len() >= 2);
        assert!(forward[0].is_equal(backward[backward.len() - 1], T));
        assert!(forward[forward.len() - 1].is_equal(backward[0], T));
    }

    #[test]
    fn triangulating_something_that_is_not_a_face_is_refused() {
        let mut model = Model::new();
        let built = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
        assert!(triangulate_face(&model, &built.shape, fine(), T).is_err());

        let vertex = explore_unique(&model, &built.shape, ShapeType::Vertex).unwrap()[0].clone();
        assert!(triangulate_face(&model, &vertex, fine(), T).is_err());
        assert!(polyline_of_edge(&model, &vertex, fine(), T).is_err());
    }

    #[test]
    fn a_triangles_normal_is_none_when_it_is_degenerate() {
        let a = Point::ORIGIN;
        let b = Point::new(1.0, 0.0, 0.0);
        assert!(triangle_normal(a, b, Point::new(0.0, 1.0, 0.0), T).is_some());
        assert!(triangle_normal(a, b, Point::new(2.0, 0.0, 0.0), T).is_none());
        assert!(triangle_normal(a, a, a, T).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod predicate_tests {
    use super::*;
    use ogeom_core::Fast;

    /// A triangle with one very long, very nearly diagonal edge.
    ///
    /// The configuration where a naive crossing test goes wrong. Subtracting
    /// the edge's ends from a query point close to the line cancels almost
    /// every significant digit, and what is left is rounding rather than
    /// geometry.
    fn sliver() -> Vec<Vec<Point2>> {
        vec![vec![
            Point2::new(0.5, 0.5),
            Point2::new(1000.0, 1000.0),
            Point2::new(1000.0, 0.5),
        ]]
    }

    #[test]
    fn the_two_implementations_are_a_real_choice_and_not_a_decoration() {
        // The point of the seam. In a band of near-degenerate queries the two
        // answer differently — and if they never did, routing through the trait
        // would be ceremony rather than a design.
        //
        // The query points straddle the long edge at a spacing far below what
        // the subtraction can resolve, which is exactly the case a mesh hits
        // when a face's boundary passes close to a triangulation vertex.
        let rings = sliver();
        let step = 2.0_f64.powi(-48);
        let mut disagreements = 0;
        for i in 0..256i32 {
            for j in 0..256i32 {
                let p = Point2::new(
                    250.0 + f64::from(i - 128) * step,
                    250.0 + f64::from(j - 128) * step,
                );
                if inside_boundary_with::<Exact>(&rings, p)
                    != inside_boundary_with::<Fast>(&rings, p)
                {
                    disagreements += 1;
                }
            }
        }
        assert!(
            disagreements > 0,
            "the exact and fast predicates never disagreed, so the seam is not \
             carrying anything"
        );
    }

    #[test]
    fn the_exact_predicate_is_the_one_that_is_right_near_the_edge() {
        // Disagreeing is not enough — the exact one has to be the *correct*
        // one, and that needs points whose side is known without asking either
        // implementation.
        //
        // Stepped in units of the last place rather than by a small distance.
        // A distance below the spacing of `f64` at 250 rounds away entirely,
        // leaving two points that are both exactly *on* the diagonal — where
        // either answer is defensible and the test would be asserting nothing.
        let rings = sliver();
        let base = 250.0_f64;
        let mut off = base;
        for k in 1..64 {
            off = off.next_up();
            assert_ne!(off, base, "step {k} did not move the point at all");
            // Inside this triangle is below the diagonal `y = x`: larger x.
            assert!(
                inside_boundary_with::<Exact>(&rings, Point2::new(off, base)),
                "a point {k} ulps below the diagonal was reported outside"
            );
            assert!(
                !inside_boundary_with::<Exact>(&rings, Point2::new(base, off)),
                "a point {k} ulps above the diagonal was reported inside"
            );
        }
    }

    #[test]
    fn the_exact_answer_is_the_one_the_geometry_supports() {
        // Points placed by construction, so the right answer is known without
        // asking either implementation.
        let square = vec![vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ]];
        assert!(inside_boundary_with::<Exact>(
            &square,
            Point2::new(0.5, 0.5)
        ));
        assert!(!inside_boundary_with::<Exact>(
            &square,
            Point2::new(1.5, 0.5)
        ));
        assert!(!inside_boundary_with::<Exact>(
            &square,
            Point2::new(-0.5, 0.5)
        ));
        assert!(!inside_boundary_with::<Exact>(
            &square,
            Point2::new(0.5, 1.5)
        ));

        // A vertex exactly on the sampling ray is counted once, not twice or
        // not at all — which is what the half-open comparison is for.
        let diamond = vec![vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 0.0),
            Point2::new(1.0, -1.0),
        ]];
        assert!(inside_boundary_with::<Exact>(
            &diamond,
            Point2::new(1.0, 0.0)
        ));
        assert!(!inside_boundary_with::<Exact>(
            &diamond,
            Point2::new(3.0, 0.0)
        ));
        assert!(!inside_boundary_with::<Exact>(
            &diamond,
            Point2::new(-1.0, 0.0)
        ));
    }

    #[test]
    fn a_hole_is_outside_however_either_ring_is_wound() {
        // Even-odd does not depend on the wires being wound consistently, which
        // is what makes it survive imported geometry.
        let outer = vec![
            Point2::new(0.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(4.0, 4.0),
            Point2::new(0.0, 4.0),
        ];
        let hole: Vec<Point2> = vec![
            Point2::new(1.0, 1.0),
            Point2::new(3.0, 1.0),
            Point2::new(3.0, 3.0),
            Point2::new(1.0, 3.0),
        ];
        let backwards: Vec<Point2> = hole.iter().rev().copied().collect();
        for inner in [hole, backwards] {
            let rings = vec![outer.clone(), inner];
            assert!(!inside_boundary_with::<Exact>(
                &rings,
                Point2::new(2.0, 2.0)
            ));
            assert!(inside_boundary_with::<Exact>(&rings, Point2::new(0.5, 0.5)));
        }
    }
}
