//! Offsetting and thickening sheets: a face or an open shell moved along
//! its own normals, and the solid between a sheet and its offset.
//!
//! Every face moves on its own chart: the moved surface answers, at the old
//! surface's `(u, v)`, the old point moved the distance along the normal
//! there. An analytic surface has such a parallel in its own family (a
//! plane translated, a cylinder, cone, sphere or torus with its radius
//! changed), and it is used once it is measured to agree with the moved
//! points. Any other surface (a B-spline, a revolution, an extrusion) is
//! fitted to the moved points at their own parameters, refined until the
//! fit, measured *between* the fitted points against the exact offset, is
//! within the approximation tolerance; the measured deviation widens the
//! moved face's tolerance, which is where it is reported. Because the chart
//! is kept, every pcurve of the sheet carries over unchanged, so a trimmed
//! face, holes and all, is trimmed along the same chart curves as its original.
//!
//! Edges and vertices move the same way: each point along the normal of
//! the faces it bounds. A line or a circle whose move is a translation or
//! a similarity stays a line or a circle (checked at samples, not assumed);
//! any other edge is fitted at its own parameters to the same tolerance.
//! Where two faces meet, both must move the shared boundary to the same
//! place: faces meeting tangentially do, faces meeting at a crease do not,
//! and a crease is refused by name rather than torn or patched.
//!
//! Thickening builds the sheet moved to both of its sides (or the sheet
//! itself and one side) and closes the gap along every free edge with the
//! ruled face between the edge and its offset: a plane along a straight
//! edge of a constant normal, a cylinder along a circle whose normal runs
//! along its axis, and otherwise a B-spline fitted to the rulings on the
//! edge's own parameter, its deviation recorded the same way. The pieces
//! share their edges by construction, so the shell closes without sewing.
//!
//! Refused by name: a distance beyond a face's smallest radius of curvature
//! on the side it moves to (checked over the face's chart window), a
//! free-form face whose normal is undefined somewhere the offset needs it,
//! a crease, an edge shared by more than two faces, and a thickened sheet
//! whose layers run into each other.

use ogeom_algo::{
    Built, History, attach_pcurve, attach_seam, edge_vertices, make_edge_between, make_face_on,
    make_shell, make_solid, make_vertex, make_wire,
};
use ogeom_core::{OgeomResult, Tolerance, Tolerances, ogeom_bail};
use ogeom_geom::{
    Curve, Curve2d as _, Curve3d as _, CylinderSurface, Line2d, LineCurve, OffsetSurface,
    PlanarCurve, PlaneSurface, Surface as _, SurfaceGeometry, Transformable as _,
};
use ogeom_math::{Axis2, Cylinder, Direction, Direction2, Frame, Point, Point2, Transform, Vector};
use ogeom_topo::{
    EdgeData, EdgeRepr, FaceData, Filter, Location, Model, NodeData, Orientation, Shape, ShapeType,
    SurfaceId, TShapeId, explore,
};
use std::collections::HashMap;

/// Offset a sheet (a face or a shell, open or closed) by `distance` along
/// its normals: each face's surface moves on its own chart, its trim and
/// its joins to its neighbours carried over one for one, so the history
/// maps every face, edge and vertex to the one it became.
///
/// A face with no parallel in its own family is fitted within the
/// approximation tolerance, and its tolerance (and so its edges' and
/// vertices') covers the deviation the fit measured.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// the shape is not a face or a shell, `distance` moves nothing, the
/// distance reaches past a face's smallest radius of curvature on the side
/// it moves to, the faces meet at a crease (the offsets of the two sides
/// part or cross there), an edge is shared by more than two faces, or a
/// free-form face has no normal where the offset needs one;
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if a fit does
/// not reach the tolerance.
pub fn offset_sheet(
    model: &mut Model,
    sheet: &Shape,
    distance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !distance.is_finite() || distance.abs() <= tol.confusion() {
        ogeom_bail!(Construction, "an offset of {distance} moves nothing");
    }
    let read = read_sheet(model, sheet, tol)?;
    let layer = moved_layer(model, &read, distance, tol)?;
    let shape = if read.root_is_face {
        layer.faces[0].clone()
    } else {
        make_shell(model, &layer.faces)?.shape
    };
    let mut history = layer.history;
    if !read.root_is_face {
        history.modify(&read.root, shape.clone());
    }
    if layer.faces.len() > 1 && !ogeom_algo::check_self_intersection(model, &shape, tol)?.is_empty()
    {
        ogeom_bail!(
            Construction,
            "the offset sheet runs into itself; the distance is larger than \
             the room between its faces"
        );
    }
    Ok(Built::new(shape, with_prefix(read.prefix, history)))
}

/// Thicken a sheet (a face or an open shell) into a solid: the sheet
/// moved `thickness` along its normals (a negative thickness against
/// them), or by half of it to each side when `both_sides`, the two layers
/// joined by side faces along the sheet's free edges.
///
/// The layers are built as [`offset_sheet`] builds them, so a free-form
/// face is fitted within the approximation tolerance and reports its
/// deviation in its tolerance. Each side face is the ruled face between a
/// free edge's two images: a plane or a cylinder where the rulings make
/// one on the edge's own parameter, and otherwise a B-spline fitted to the
/// rulings, its deviation recorded the same way.
///
/// # Errors
///
/// As [`offset_sheet`], and if the sheet is closed or has no free edge (a
/// closed shell bounds a solid already; hollow it with
/// [`make_thick_solid`](crate::make_thick_solid)), a free edge runs into a
/// point where its face has no normal, or the layers run into each other.
pub fn make_thick_sheet(
    model: &mut Model,
    sheet: &Shape,
    thickness: f64,
    both_sides: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !thickness.is_finite() || thickness.abs() <= tol.confusion() {
        ogeom_bail!(Construction, "a thickness of {thickness} holds nothing");
    }
    let read = read_sheet(model, sheet, tol)?;
    if read.faces.iter().any(|f| f.natural) {
        ogeom_bail!(
            Construction,
            "a face with no boundary edges has nothing to close the sides \
             along; bound it with edges first"
        );
    }
    let free: Vec<usize> = (0..read.edges.len())
        .filter(|&i| read.edges[i].uses.len() == 1 && !read.edges[i].data.degenerate)
        .collect();
    if free.is_empty() {
        ogeom_bail!(
            Construction,
            "the sheet has no free edge; a closed shell is hollowed with \
             make_thick_solid, not thickened"
        );
    }
    let (lo, hi) = if both_sides {
        (-thickness.abs() / 2.0, thickness.abs() / 2.0)
    } else if thickness > 0.0 {
        (0.0, thickness)
    } else {
        (thickness, 0.0)
    };
    let lower = moved_layer(model, &read, lo, tol)?;
    let upper = moved_layer(model, &read, hi, tol)?;

    let mut history = History::new();
    let mut faces: Vec<Shape> = Vec::new();
    for (fi, face) in read.faces.iter().enumerate() {
        // The upper layer looks the way the sheet does; the lower one looks
        // back into the material.
        faces.push(upper.faces[fi].clone());
        faces.push(lower.faces[fi].reversed());
        history.generate(&face.occurrence, upper.faces[fi].clone());
        history.generate(&face.occurrence, lower.faces[fi].reversed());
    }
    let mut risers: HashMap<TShapeId, Shape> = HashMap::new();
    for ei in free {
        let side = side_face(
            model,
            &read,
            ei,
            (lo, hi),
            (&lower, &upper),
            &mut risers,
            tol,
        )?;
        history.generate(&Shape::of(read.edges[ei].node), side.clone());
        faces.push(side);
    }
    let shell = make_shell(model, &faces)?.shape;
    if !ogeom_algo::is_shell_closed(model, &shell)? {
        ogeom_bail!(
            Construction,
            "the thickened sheet does not close: an edge bounds the sheet \
             more than once in a way the side faces cannot follow"
        );
    }
    let solid = make_solid(model, std::slice::from_ref(&shell))?.shape;
    if !ogeom_algo::inside_out_faces(model, &solid, tol)?.is_empty() {
        ogeom_bail!(
            Construction,
            "the thickened sheet folds through itself and turns faces inside out"
        );
    }
    if !ogeom_algo::check_self_intersection(model, &solid, tol)?.is_empty() {
        ogeom_bail!(
            Construction,
            "the thickened sheet runs into itself; the thickness is larger \
             than the room between its faces"
        );
    }
    if !read.root_is_face {
        history.generate(&read.root, solid.clone());
    }
    Ok(Built::new(solid, with_prefix(read.prefix, history)))
}

fn with_prefix(prefix: Option<History>, history: History) -> History {
    match prefix {
        Some(prefix) => prefix.then(&history),
        None => history,
    }
}

/// One face of the sheet, as the sheet holds it.
struct SheetFace {
    occurrence: Shape,
    surface_id: SurfaceId,
    surface: SurfaceGeometry,
    /// `+1` where the sheet's normal is the surface's own, `-1` against it.
    sign: f64,
    natural: bool,
    /// The stored wires: each wire's sense and its edges' nodes and senses.
    wires: Vec<(Orientation, Vec<(TShapeId, Orientation)>)>,
    /// The region of the chart the face covers.
    window: ((f64, f64), (f64, f64)),
}

/// One use of an edge by a face of the sheet.
struct EdgeUse {
    face: usize,
    /// The edge's sense in the bare face.
    sense: Orientation,
    /// The pcurve on the face's surface and the range it runs over, if the
    /// edge carries one.
    pcurve: Option<(PlanarCurve, (f64, f64))>,
}

/// One edge of the sheet.
struct SheetEdge {
    node: TShapeId,
    data: EdgeData,
    /// The curve in space and the edge's range on it; none for a pole.
    curve: Option<(Curve, (f64, f64))>,
    /// The start and end vertex nodes.
    ends: (TShapeId, TShapeId),
    uses: Vec<EdgeUse>,
}

/// A sheet read into the pieces the offset walks.
struct Sheet {
    root: Shape,
    root_is_face: bool,
    faces: Vec<SheetFace>,
    edges: Vec<SheetEdge>,
    /// Each vertex node, its point and where it sits in each face's chart.
    vertices: Vec<SheetVertex>,
    prefix: Option<History>,
}

/// A vertex node, its point, and where it sits in the chart of each face
/// it bounds.
type SheetVertex = (TShapeId, Point, Seats);

/// Places in face charts: the face's index and the parameters there.
type Seats = Vec<(usize, Point2)>;

/// The sheet as faces, edges and vertices, every placement baked in.
fn read_sheet(model: &mut Model, sheet: &Shape, tol: Tolerances) -> OgeomResult<Sheet> {
    let kind = model.kind_of(sheet)?;
    match kind {
        ShapeType::Face | ShapeType::Shell => {}
        ShapeType::Solid | ShapeType::CompSolid => ogeom_bail!(
            Construction,
            "a solid is offset with offset_shape and hollowed with \
             make_thick_solid; a sheet is a face or a shell"
        ),
        other => ogeom_bail!(
            Construction,
            "a sheet is a face or a shell, not a {other:?}"
        ),
    }
    let (root, prefix) = if placed(model, sheet)? {
        baked_sheet(model, sheet, kind, tol)?
    } else {
        (sheet.clone(), None)
    };

    let mut faces: Vec<SheetFace> = Vec::new();
    let mut seen: Vec<TShapeId> = Vec::new();
    for occurrence in explore(model, &root, Filter::OfType(ShapeType::Face))? {
        if seen.contains(&occurrence.node()) {
            ogeom_bail!(Construction, "the sheet holds one face twice");
        }
        seen.push(occurrence.node());
        let Some(node) = model.node(&occurrence) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface).cloned() else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
        };
        let mut wires = Vec::new();
        for wire in node.children() {
            let Some(wire_node) = model.node(wire) else {
                ogeom_bail!(Dangling, "wire is not in this model");
            };
            let edges = wire_node
                .children()
                .iter()
                .map(|e| (e.node(), e.orientation()))
                .collect();
            wires.push((wire.orientation(), edges));
        }
        faces.push(SheetFace {
            occurrence: occurrence.clone(),
            surface_id: data.surface,
            sign: if occurrence.orientation() == Orientation::Reversed {
                -1.0
            } else {
                1.0
            },
            natural: data.natural_restriction || wires.is_empty(),
            surface,
            wires,
            window: ((0.0, 0.0), (0.0, 0.0)),
        });
    }
    if faces.is_empty() {
        ogeom_bail!(Construction, "the sheet has no faces");
    }

    let mut edges: Vec<SheetEdge> = Vec::new();
    let mut edge_index: HashMap<TShapeId, usize> = HashMap::new();
    for (fi, face) in faces.iter().enumerate() {
        for (wire_sense, list) in &face.wires {
            for (edge_node, edge_sense) in list {
                let sense = wire_sense.compose(*edge_sense);
                let index = if let Some(i) = edge_index.get(edge_node) {
                    *i
                } else {
                    let bare = Shape::of(*edge_node);
                    let Some(data) = model.node(&bare).and_then(|n| n.data().as_edge()).cloned()
                    else {
                        ogeom_bail!(Construction, "edge node holds no edge data");
                    };
                    let curve = match data.curve3d() {
                        Some(EdgeRepr::Curve3d { curve, range, .. }) => {
                            let Some(geometry) = model.geometry().curve(*curve) else {
                                ogeom_bail!(Dangling, "curve is not in this model");
                            };
                            Some((geometry.clone(), *range))
                        }
                        _ if data.degenerate => None,
                        _ => ogeom_bail!(Construction, "a sheet edge has no curve in space"),
                    };
                    let Some((start, end)) = edge_vertices(model, &bare)? else {
                        ogeom_bail!(Construction, "a sheet edge has no vertices");
                    };
                    edges.push(SheetEdge {
                        node: *edge_node,
                        data,
                        curve,
                        ends: (start.node(), end.node()),
                        uses: Vec::new(),
                    });
                    edge_index.insert(*edge_node, edges.len() - 1);
                    edges.len() - 1
                };
                let edge = &mut edges[index];
                let pcurve = match edge.data.pcurve_on(face.surface_id) {
                    Some(EdgeRepr::PCurve { curve, range, .. }) => model
                        .geometry()
                        .pcurve(*curve)
                        .cloned()
                        .map(|p| (p, *range)),
                    Some(EdgeRepr::Seam {
                        forward,
                        reversed,
                        range,
                        ..
                    }) => {
                        let id = if sense == Orientation::Reversed {
                            *reversed
                        } else {
                            *forward
                        };
                        model.geometry().pcurve(id).cloned().map(|p| (p, *range))
                    }
                    _ => None,
                };
                if pcurve.is_none() && edge.data.degenerate {
                    ogeom_bail!(Construction, "a pole edge has no pcurve to place it");
                }
                edge.uses.push(EdgeUse {
                    face: fi,
                    sense,
                    pcurve,
                });
            }
        }
    }
    for edge in &edges {
        let faces_met: std::collections::HashSet<usize> =
            edge.uses.iter().map(|u| u.face).collect();
        if edge.uses.len() > 2 {
            ogeom_bail!(
                Construction,
                "an edge is shared by more than two faces; the sheet is not a \
                 surface the offset can move"
            );
        }
        if faces_met.len() == 1 && edge.uses.len() == 2 && edge.uses[0].sense == edge.uses[1].sense
        {
            ogeom_bail!(Construction, "a seam is walked the same way twice");
        }
    }

    // Each face's window, and each vertex's seat in every chart it is on.
    let mut vertices: Vec<SheetVertex> = Vec::new();
    let mut vertex_index: HashMap<TShapeId, usize> = HashMap::new();
    let mut seat = |model: &Model,
                    vertex: TShapeId,
                    face: usize,
                    uv: Option<Point2>,
                    surface: &SurfaceGeometry|
     -> OgeomResult<()> {
        let index = if let Some(i) = vertex_index.get(&vertex) {
            *i
        } else {
            let Some(data) = model
                .node(&Shape::of(vertex))
                .and_then(|n| n.data().as_vertex())
            else {
                ogeom_bail!(Construction, "vertex node holds no point");
            };
            vertices.push((vertex, data.point, Vec::new()));
            vertex_index.insert(vertex, vertices.len() - 1);
            vertices.len() - 1
        };
        let uv = match uv {
            Some(uv) => uv,
            None => {
                let found = ogeom_algo::project_on_surface(surface, vertices[index].1, 32, tol)?;
                Point2::new(found.parameters.0, found.parameters.1)
            }
        };
        vertices[index].2.push((face, uv));
        Ok(())
    };
    let mut boxes: Vec<Option<(Point2, Point2)>> = vec![None; faces.len()];
    let grow = |boxes: &mut Vec<Option<(Point2, Point2)>>, face: usize, uv: Point2| {
        boxes[face] = Some(match boxes[face] {
            None => (uv, uv),
            Some((lo, hi)) => (
                Point2::new(lo.x.min(uv.x), lo.y.min(uv.y)),
                Point2::new(hi.x.max(uv.x), hi.y.max(uv.y)),
            ),
        });
    };
    for edge in &edges {
        for used in &edge.uses {
            let surface = &faces[used.face].surface;
            let (start_uv, end_uv) = if let Some((p, range)) = &used.pcurve {
                for k in 0..=24 {
                    let t = range.0 + (range.1 - range.0) * f64::from(k) / 24.0;
                    grow(&mut boxes, used.face, p.point_at(t, tol)?);
                }
                (
                    Some(p.point_at(range.0, tol)?),
                    Some(p.point_at(range.1, tol)?),
                )
            } else if let Some((curve, range)) = &edge.curve {
                for k in 0..=24 {
                    let t = range.0 + (range.1 - range.0) * f64::from(k) / 24.0;
                    let found =
                        ogeom_algo::project_on_surface(surface, curve.point_at(t, tol)?, 32, tol)?;
                    grow(
                        &mut boxes,
                        used.face,
                        Point2::new(found.parameters.0, found.parameters.1),
                    );
                }
                (None, None)
            } else {
                (None, None)
            };
            seat(model, edge.ends.0, used.face, start_uv, surface)?;
            seat(model, edge.ends.1, used.face, end_uv, surface)?;
        }
    }
    for (fi, face) in faces.iter_mut().enumerate() {
        let ((ua, ub), (va, vb)) = face.surface.domain();
        face.window = match boxes[fi] {
            Some((lo, hi)) if !face.natural => {
                // A pcurve bulges between its samples by a little; the
                // margin holds it, the surface's own domain caps it where
                // the surface is not periodic.
                let mu = (hi.x - lo.x).max(tol.parametric()) * 0.02;
                let mv = (hi.y - lo.y).max(tol.parametric()) * 0.02;
                let (mut u0, mut u1, mut v0, mut v1) = (lo.x - mu, hi.x + mu, lo.y - mv, hi.y + mv);
                if !face.surface.is_periodic_u() {
                    u0 = u0.max(ua);
                    u1 = u1.min(ub);
                }
                if !face.surface.is_periodic_v() {
                    v0 = v0.max(va);
                    v1 = v1.min(vb);
                }
                ((u0, u1), (v0, v1))
            }
            _ => ((ua, ub), (va, vb)),
        };
        let ((u0, u1), (v0, v1)) = face.window;
        if ![u0, u1, v0, v1].iter().all(|x| x.is_finite()) || u1 <= u0 || v1 <= v0 {
            ogeom_bail!(
                Construction,
                "a face covers no finite region of its surface's chart"
            );
        }
    }
    Ok(Sheet {
        root_is_face: kind == ShapeType::Face,
        root,
        faces,
        edges,
        vertices,
        prefix,
    })
}

/// Whether anything in the sheet stands away from its node: a placed
/// occurrence, or geometry stored under a placement.
fn placed(model: &Model, sheet: &Shape) -> OgeomResult<bool> {
    for kind in [ShapeType::Vertex, ShapeType::Edge, ShapeType::Face] {
        for occurrence in explore(model, sheet, Filter::OfType(kind))? {
            if occurrence.transform(model.datums())?.kind() != ogeom_math::TransformKind::Identity {
                return Ok(true);
            }
            let Some(node) = model.node(&occurrence) else {
                ogeom_bail!(Dangling, "shape is not in this model");
            };
            match node.data() {
                NodeData::Face(data) if !data.location.is_identity() => return Ok(true),
                NodeData::Edge(data)
                    if data
                        .representations
                        .iter()
                        .any(|r| r.location().is_some_and(|l| !l.is_identity())) =>
                {
                    return Ok(true);
                }
                _ => {}
            }
        }
    }
    Ok(false)
}

/// The sheet restated in world coordinates, with the history from the
/// placed sheet to the restated one.
fn baked_sheet(
    model: &mut Model,
    sheet: &Shape,
    kind: ShapeType,
    tol: Tolerances,
) -> OgeomResult<(Shape, Option<History>)> {
    let shell = if kind == ShapeType::Face {
        model.add_shell(std::slice::from_ref(sheet))?
    } else {
        sheet.clone()
    };
    let holder = model.add_solid(std::slice::from_ref(&shell))?;
    let baked = ogeom_algo::baked_shape(model, &holder, tol)?;
    let root = if kind == ShapeType::Face {
        let faces = explore(model, &baked.shape, Filter::OfType(ShapeType::Face))?;
        let [face] = faces.as_slice() else {
            ogeom_bail!(Construction, "a placed face did not restate as one face");
        };
        face.clone()
    } else {
        let shells = explore(model, &baked.shape, Filter::OfType(ShapeType::Shell))?;
        let [one] = shells.as_slice() else {
            ogeom_bail!(Construction, "a placed shell did not restate as one shell");
        };
        one.clone()
    };
    let mut history = baked.history;
    if history.modified(sheet).is_empty() {
        history.modify(sheet, root.clone());
    }
    Ok((root, Some(history)))
}

/// A face's surface moved by the offset.
struct MovedSurface {
    id: SurfaceId,
    geometry: SurfaceGeometry,
    /// Whether the geometry is the exact parallel (rather than a fit).
    exact: bool,
    /// The fit's measured deviation from the exact offset; zero when exact.
    deviation: f64,
}

/// The sheet moved by one distance.
struct Layer {
    /// Faces as the sheet holds them, index for index.
    faces: Vec<Shape>,
    /// Each sheet edge's image, bare and forward.
    edges: HashMap<TShapeId, Shape>,
    /// Each sheet vertex's image.
    vertices: HashMap<TShapeId, (Shape, Point)>,
    history: History,
}

/// The sheet moved `distance` along its normals; at zero, a copy on the
/// same geometry.
fn moved_layer(
    model: &mut Model,
    sheet: &Sheet,
    distance: f64,
    tol: Tolerances,
) -> OgeomResult<Layer> {
    let target = tol.approximation();
    let mut history = History::new();

    let mut surfaces: Vec<MovedSurface> = Vec::with_capacity(sheet.faces.len());
    for face in &sheet.faces {
        let (geometry, exact, deviation) = if distance == 0.0 {
            (face.surface.clone(), true, 0.0)
        } else {
            moved_surface(face, distance, target, tol)?
        };
        let id = model.geometry_mut().add_surface(geometry.clone());
        surfaces.push(MovedSurface {
            id,
            geometry,
            exact,
            deviation,
        });
    }

    // Vertices: every face a vertex bounds must move it to one place.
    let mut vertices: HashMap<TShapeId, (Shape, Point)> = HashMap::new();
    let mut vertex_slop: HashMap<TShapeId, f64> = HashMap::new();
    for (node, point, seats) in &sheet.vertices {
        let mut candidates = Vec::with_capacity(seats.len());
        for (fi, uv) in seats {
            candidates.push(moved_point(
                sheet, &surfaces, *fi, *uv, *point, distance, tol,
            )?);
        }
        let (moved, spread) = agreed(&candidates, target, "a vertex")?;
        let shape = make_vertex(model, moved).shape;
        history.modify(&Shape::of(*node), shape.clone());
        vertex_slop.insert(*node, spread);
        vertices.insert(*node, (shape, moved));
    }

    // Edges.
    let mut edges: HashMap<TShapeId, Shape> = HashMap::new();
    for edge in &sheet.edges {
        let (Some((from, from_at)), Some((to, to_at))) = (
            vertices.get(&edge.ends.0).cloned(),
            vertices.get(&edge.ends.1).cloned(),
        ) else {
            ogeom_bail!(Construction, "an edge end has no moved vertex");
        };
        let built = if let Some((curve, range)) = &edge.curve {
            let off = |t: f64| -> OgeomResult<Point> {
                let at = curve.point_at(t, tol)?;
                let mut candidates = Vec::with_capacity(edge.uses.len());
                for used in &edge.uses {
                    let uv = seat_on(sheet, edge, used, (at, t, *range), tol)?;
                    candidates.push(moved_point(
                        sheet, &surfaces, used.face, uv, at, distance, tol,
                    )?);
                }
                Ok(agreed(&candidates, target, "an edge")?.0)
            };
            let (moved, slop) = if distance == 0.0 {
                (curve.clone(), 0.0)
            } else if let Some(exact) = exact_moved_curve(curve, *range, &off, tol)? {
                (exact, 0.0)
            } else {
                fitted_moved_curve(*range, &off, target, tol)?
            };
            // The curve's ends must reach the vertices: within the slop of
            // the fit and the spread the vertex's faces agreed within.
            for (t, (vertex, at)) in [(range.0, (&from, from_at)), (range.1, (&to, to_at))] {
                let gap = moved.point_at(t, tol)?.distance(at);
                if gap > tol.confusion() {
                    model.widen(vertex, Tolerance::new(gap + tol.confusion())?)?;
                }
            }
            let shape = make_edge_between(model, moved, *range, &from, &to, tol)?.shape;
            if slop > 0.0 {
                model.widen(&shape, Tolerance::new(slop + tol.confusion())?)?;
            }
            shape
        } else {
            let mut data = EdgeData::new();
            data.degenerate = true;
            model.add_edge(data, &[from.clone(), to.clone()])?
        };
        for end in [edge.ends.0, edge.ends.1] {
            let slop = vertex_slop.get(&end).copied().unwrap_or(0.0);
            if slop > tol.confusion() {
                model.widen(&vertices[&end].0, Tolerance::new(slop + tol.confusion())?)?;
            }
        }
        // The pcurves carry over: the moved surface keeps the chart.
        let mut done: Vec<SurfaceId> = Vec::new();
        for used in &edge.uses {
            let sid = surfaces[used.face].id;
            if done.contains(&sid) {
                continue;
            }
            done.push(sid);
            let old = sheet.faces[used.face].surface_id;
            match edge.data.pcurve_on(old) {
                Some(EdgeRepr::PCurve { curve, range, .. }) => {
                    let Some(p) = model.geometry().pcurve(*curve).cloned() else {
                        ogeom_bail!(Dangling, "pcurve is not in this model");
                    };
                    attach_pcurve(model, &built, p, sid, Location::identity(), *range)?;
                }
                Some(EdgeRepr::Seam {
                    forward,
                    reversed,
                    range,
                    ..
                }) => {
                    let (Some(f), Some(r)) = (
                        model.geometry().pcurve(*forward).cloned(),
                        model.geometry().pcurve(*reversed).cloned(),
                    ) else {
                        ogeom_bail!(Dangling, "pcurve is not in this model");
                    };
                    attach_seam(model, &built, f, r, sid, Location::identity(), *range)?;
                }
                _ => {}
            }
        }
        // The sheet's own slop moves with it.
        model.widen(&built, edge.data.tolerance)?;
        same_parameter(model, &built, tol)?;
        history.modify(&Shape::of(edge.node), built.clone());
        edges.insert(edge.node, built);
    }

    // Faces, mirroring the sheet's own structure.
    let mut faces: Vec<Shape> = Vec::with_capacity(sheet.faces.len());
    for (fi, face) in sheet.faces.iter().enumerate() {
        let mut wires = Vec::with_capacity(face.wires.len());
        for (wire_sense, list) in &face.wires {
            let ring: Vec<Shape> = list
                .iter()
                .map(|(node, sense)| edges[node].oriented(*sense))
                .collect();
            wires.push(model.add_wire(&ring)?.oriented(*wire_sense));
        }
        let data = if face.natural {
            FaceData::natural(surfaces[fi].id, Location::identity())
        } else {
            FaceData::new(surfaces[fi].id, Location::identity())
        };
        let bare = model.add_face(data, &wires)?;
        if surfaces[fi].deviation > 0.0 {
            model.widen(
                &bare,
                Tolerance::new(surfaces[fi].deviation + tol.confusion())?,
            )?;
        }
        let shape = bare.oriented(face.occurrence.orientation());
        history.modify(&face.occurrence, shape.clone());
        faces.push(shape);
    }
    Ok(Layer {
        faces,
        edges,
        vertices,
        history,
    })
}

/// Where the point `at`, at `t` on the edge's curve over `range`, sits
/// in a face's chart: read off the pcurve, its own range matched to the
/// curve's end for end, where that lands on the point, and projected
/// where it does not (a pcurve parameterized apart from its curve).
fn seat_on(
    sheet: &Sheet,
    edge: &SheetEdge,
    used: &EdgeUse,
    (at, t, range): (Point, f64, (f64, f64)),
    tol: Tolerances,
) -> OgeomResult<Point2> {
    let surface = &sheet.faces[used.face].surface;
    if let Some((p, own)) = &used.pcurve {
        let s = if range == *own {
            t
        } else {
            own.0 + (own.1 - own.0) * (t - range.0) / (range.1 - range.0)
        };
        let uv = p.point_at(s, tol)?;
        let reach = edge.data.tolerance.get().max(tol.confusion()) * 10.0;
        if surface.point_at(uv.x, uv.y, tol)?.distance(at) <= reach {
            return Ok(uv);
        }
    }
    let found = ogeom_algo::project_on_surface(surface, at, 32, tol)?;
    Ok(Point2::new(found.parameters.0, found.parameters.1))
}

/// A point of a face, at `uv` in its chart, moved `distance` along the
/// sheet's normal there. Where the normal is undefined (a pole, an apex)
/// the exact moved surface still answers at the same parameters; a fitted
/// one has nothing to answer with.
fn moved_point(
    sheet: &Sheet,
    surfaces: &[MovedSurface],
    face: usize,
    uv: Point2,
    at: Point,
    distance: f64,
    tol: Tolerances,
) -> OgeomResult<Point> {
    if distance == 0.0 {
        return Ok(at);
    }
    if let Some(n) = sheet_normal(&sheet.faces[face], uv, tol) {
        return Ok(at + n * distance);
    }
    if surfaces[face].exact {
        return surfaces[face].geometry.point_at(uv.x, uv.y, tol);
    }
    ogeom_bail!(
        Construction,
        "a free-form face has no normal at ({}, {}); its offset there is \
         not determined",
        uv.x,
        uv.y
    )
}

/// The sheet's unit normal on a face, or `None` where the chart is
/// degenerate and determines none.
fn sheet_normal(face: &SheetFace, uv: Point2, tol: Tolerances) -> Option<Vector> {
    face.surface
        .normal_at(uv.x, uv.y, tol)
        .ok()
        .map(|n| n.vector() * face.sign)
}

/// The one place several faces move a shared point to, and how far the
/// farthest of them stands from it; refused as a crease when they part by
/// more than `target`.
fn agreed(candidates: &[Point], target: f64, what: &str) -> OgeomResult<(Point, f64)> {
    let Some(first) = candidates.first() else {
        ogeom_bail!(Construction, "{what} bounds no face of the sheet");
    };
    #[allow(clippy::cast_precision_loss, reason = "a handful of faces")]
    let count = candidates.len() as f64;
    let sum = candidates
        .iter()
        .fold(Vector::ZERO, |acc, p| acc + (*p - *first));
    let mean = *first + sum / count;
    let spread = candidates
        .iter()
        .fold(0.0_f64, |acc, p| acc.max(p.distance(mean)));
    if spread > target {
        ogeom_bail!(
            Construction,
            "the sheet's faces meet at a crease: {what} they share moves \
             {} apart on the two sides, which the offset cannot join",
            2.0 * spread
        );
    }
    Ok((mean, spread))
}

/// The sample count along a window or a range when measuring.
const PROBES: i32 = 8;

/// A face's surface moved `distance` along the sheet's normal: the exact
/// parallel where one is measured to keep the chart, else a fit. Returns
/// the geometry, whether it is exact, and the fit's measured deviation.
fn moved_surface(
    face: &SheetFace,
    distance: f64,
    target: f64,
    tol: Tolerances,
) -> OgeomResult<(SurfaceGeometry, bool, f64)> {
    let along = face.sign * distance;
    // An offset surface moves on as its basis offset by the sum.
    let (base, carried) = match &face.surface {
        SurfaceGeometry::Offset(o) => (o.basis(), o.distance()),
        other => (other, 0.0),
    };
    let total = carried + along;
    let ((u0, u1), (v0, v1)) = face.window;
    let probe = |i: i32, j: i32| {
        (
            u0 + (u1 - u0) * f64::from(i) / f64::from(PROBES),
            v0 + (v1 - v0) * f64::from(j) / f64::from(PROBES),
        )
    };

    // The fold: moved toward the side a surface curves to, by its radius
    // of curvature or more, the parallel turns back on itself.
    for i in 0..=PROBES {
        for j in 0..=PROBES {
            let (u, v) = probe(i, j);
            let Some(curvatures) = principal_curvatures(base, u, v, tol)? else {
                continue;
            };
            for k in curvatures {
                if total * k >= 1.0 - 1e-9 {
                    ogeom_bail!(
                        Construction,
                        "an offset of {distance} reaches past the face's \
                         radius of curvature of {} at ({u}, {v}); the moved \
                         face would fold over itself",
                        1.0 / k.abs()
                    );
                }
            }
        }
    }

    let moved_at = |u: f64, v: f64| -> OgeomResult<Option<Point>> {
        let p = face.surface.point_at(u, v, tol)?;
        Ok(face
            .surface
            .normal_at(u, v, tol)
            .ok()
            .map(|n| p + n.vector() * along))
    };

    // The exact parallel, measured against the moved points before use.
    if let Some(parallel) = OffsetSurface::new(base.clone(), total)?.analytic(tol)? {
        let mut worst = 0.0_f64;
        for i in 0..=PROBES {
            for j in 0..=PROBES {
                let (u, v) = probe(i, j);
                if let Some(want) = moved_at(u, v)? {
                    worst = worst.max(parallel.point_at(u, v, tol)?.distance(want));
                }
            }
        }
        if worst <= tol.confusion() * 10.0 {
            return Ok((parallel, true, 0.0));
        }
    }

    // A fit at the chart's own parameters, refined until the points it was
    // not fitted to (the cell centres and edge middles) agree too.
    let mut spans: u32 = 8;
    loop {
        let n = spans;
        let at = |k: u32, lo: f64, hi: f64| lo + (hi - lo) * f64::from(k) / f64::from(2 * n);
        let mut fine: Vec<Vec<Point>> = Vec::with_capacity((2 * n + 1) as usize);
        for j in 0..=2 * n {
            let v = at(j, v0, v1);
            let mut row = Vec::with_capacity((2 * n + 1) as usize);
            for i in 0..=2 * n {
                let u = at(i, u0, u1);
                let Some(p) = moved_at(u, v)? else {
                    ogeom_bail!(
                        Construction,
                        "a free-form face has no normal at ({u}, {v}); its \
                         offset there is not determined"
                    );
                };
                row.push(p);
            }
            fine.push(row);
        }
        let us: Vec<f64> = (0..=n).map(|i| at(2 * i, u0, u1)).collect();
        let vs: Vec<f64> = (0..=n).map(|j| at(2 * j, v0, v1)).collect();
        let rows: Vec<Vec<Point>> = (0..=n as usize)
            .map(|j| (0..=n as usize).map(|i| fine[2 * j][2 * i]).collect())
            .collect();
        let fitted = ogeom_geom::fit::fit_surface_grid_at(&us, &vs, &rows, 3, target * 0.5, tol)?;
        let surface: SurfaceGeometry = fitted.curve.into();
        let mut worst = fitted.error;
        for j in 0..=2 * n {
            for i in 0..=2 * n {
                let want = fine[j as usize][i as usize];
                let (u, v) = (at(i, u0, u1), at(j, v0, v1));
                worst = worst.max(surface.point_at(u, v, tol)?.distance(want));
            }
        }
        if worst <= target {
            // The fit must face the way the exact offset does everywhere it
            // was measured; a fold between the probes shows up here.
            for j in (0..=2 * n).step_by(2) {
                for i in (0..=2 * n).step_by(2) {
                    let (u, v) = (at(i, u0, u1), at(j, v0, v1));
                    let (Ok(a), Ok(b)) = (
                        face.surface.normal_at(u, v, tol),
                        surface.normal_at(u, v, tol),
                    ) else {
                        continue;
                    };
                    if a.vector().dot(b.vector()) <= 0.0 {
                        ogeom_bail!(
                            Construction,
                            "an offset of {distance} turns the face over at \
                             ({u}, {v}); it reaches past the face's radius of \
                             curvature there"
                        );
                    }
                }
            }
            return Ok((surface, false, worst));
        }
        if spans >= 128 {
            ogeom_bail!(
                NotDone,
                "the offset of a free-form face did not fit within {target}; \
                 the closest fit was {worst} off"
            );
        }
        spans *= 2;
    }
}

/// The two principal curvatures at `(u, v)`, signed positive where the
/// surface curves toward its own normal; `None` at a point where the chart
/// determines no normal. From the two fundamental forms: the curvatures
/// are the roots of `k^2 - 2Hk + K`, so no principal direction is needed,
/// which at an umbilic (everywhere on a sphere) is not determined.
fn principal_curvatures(
    surface: &SurfaceGeometry,
    u: f64,
    v: f64,
    tol: Tolerances,
) -> OgeomResult<Option<[f64; 2]>> {
    if surface.is_degenerate_at(u, v, tol)? {
        return Ok(None);
    }
    let jet = surface.jet_at(u, v, tol)?;
    let c = jet.du.cross(jet.dv);
    let scale = jet.du.magnitude().max(jet.dv.magnitude());
    if c.magnitude() <= tol.angular() * scale * scale {
        return Ok(None);
    }
    let n = c / c.magnitude();
    let (e, f, g) = (jet.du.dot(jet.du), jet.du.dot(jet.dv), jet.dv.dot(jet.dv));
    let (l, m, nn) = (jet.d2u.dot(n), jet.duv.dot(n), jet.d2v.dot(n));
    let det = e.mul_add(g, -f * f);
    let gaussian = l.mul_add(nn, -m * m) / det;
    let mean = (e * nn - 2.0 * f * m + g * l) / (2.0 * det);
    let spread = mean.mul_add(mean, -gaussian).max(0.0).sqrt();
    Ok(Some([mean + spread, mean - spread]))
}

/// A curve moved by `off` that is still a line or a circle: the move a
/// translation, or for a circle a similarity about its centre, measured
/// against `off` before it is trusted.
fn exact_moved_curve(
    curve: &Curve,
    range: (f64, f64),
    off: &dyn Fn(f64) -> OgeomResult<Point>,
    tol: Tolerances,
) -> OgeomResult<Option<Curve>> {
    let start = curve.point_at(range.0, tol)?;
    let shift = off(range.0)? - start;
    let mut candidates = vec![Transform::translation(shift)];
    if let Curve::Circle(c) = curve {
        let circle = c.circle();
        let radius = circle.radius();
        let radial = (start - circle.centre()) / radius;
        let axis = circle.frame().z().vector();
        let grown = radius + shift.dot(radial);
        if grown > tol.confusion() {
            candidates.push(
                Transform::translation(axis * shift.dot(axis))
                    * Transform::scaling(circle.centre(), grown / radius, tol)?,
            );
        }
    }
    for candidate in candidates {
        let moved = curve.transformed(&candidate, tol)?;
        let mut worst = 0.0_f64;
        for k in 0..=4 * PROBES {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / f64::from(4 * PROBES);
            worst = worst.max(moved.point_at(t, tol)?.distance(off(t)?));
        }
        if worst <= tol.confusion() * 10.0 {
            return Ok(Some(moved));
        }
    }
    Ok(None)
}

/// A curve fitted to `off` at its own parameters, refined until the
/// midpoints between the fitted samples agree; the curve and the
/// deviation measured.
fn fitted_moved_curve(
    range: (f64, f64),
    off: &dyn Fn(f64) -> OgeomResult<Point>,
    target: f64,
    tol: Tolerances,
) -> OgeomResult<(Curve, f64)> {
    let mut spans: u32 = 16;
    loop {
        let at = |k: u32| range.0 + (range.1 - range.0) * f64::from(k) / f64::from(2 * spans);
        let mut fine = Vec::with_capacity((2 * spans + 1) as usize);
        for k in 0..=2 * spans {
            fine.push(off(at(k))?);
        }
        let params: Vec<f64> = (0..=spans).map(|k| at(2 * k)).collect();
        let points: Vec<Point> = (0..=spans as usize).map(|k| fine[2 * k]).collect();
        let fitted = ogeom_geom::fit::fit_points_at(&params, &points, 3, target * 0.5, tol)?;
        let curve: Curve = fitted.curve.into();
        let mut worst = fitted.error;
        for k in 0..=2 * spans {
            let want = fine[k as usize];
            worst = worst.max(curve.point_at(at(k), tol)?.distance(want));
        }
        if worst <= target {
            return Ok((curve, worst));
        }
        if spans >= 1024 {
            ogeom_bail!(
                NotDone,
                "the offset of an edge did not fit within {target}; the \
                 closest fit was {worst} off"
            );
        }
        spans *= 2;
    }
}

/// Measure how far each pcurve of `edge`, read on its surface, strays from
/// the edge's curve, widen the edge to hold it, and record the agreement.
fn same_parameter(model: &mut Model, edge: &Shape, tol: Tolerances) -> OgeomResult<()> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()).cloned() else {
        ogeom_bail!(Construction, "edge node holds no edge data");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        return Ok(());
    };
    let Some(curve) = model.geometry().curve(*curve).cloned() else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    let mut worst = 0.0_f64;
    for repr in &data.representations {
        let (ids, surface, own) = match repr {
            EdgeRepr::PCurve {
                curve,
                surface,
                range,
                ..
            } => (vec![*curve], *surface, *range),
            EdgeRepr::Seam {
                forward,
                reversed,
                surface,
                range,
                ..
            } => (vec![*forward, *reversed], *surface, *range),
            _ => continue,
        };
        // A pcurve over its own range carries the sheet's own parameter
        // disagreement over as it was; there is no claim to measure.
        if own != *range {
            return Ok(());
        }
        let Some(surface) = model.geometry().surface(surface).cloned() else {
            ogeom_bail!(Dangling, "surface is not in this model");
        };
        for id in ids {
            let Some(p) = model.geometry().pcurve(id).cloned() else {
                ogeom_bail!(Dangling, "pcurve is not in this model");
            };
            for k in 0..=4 * PROBES {
                let t = range.0 + (range.1 - range.0) * f64::from(k) / f64::from(4 * PROBES);
                let uv = p.point_at(t, tol)?;
                let on = surface.point_at(uv.x, uv.y, tol)?;
                worst = worst.max(on.distance(curve.point_at(t, tol)?));
            }
        }
    }
    if worst > tol.confusion() {
        model.widen(edge, Tolerance::new(worst + tol.confusion())?)?;
    }
    if let Some(node) = model.node_mut(edge)
        && let NodeData::Edge(data) = node.data_mut()
    {
        data.assert_same_parameter(true);
    }
    Ok(())
}

/// The ruled face between a free edge's images in the two layers, oriented
/// out of the material, with the risers at its ends shared with the
/// neighbouring side faces through `risers`.
#[allow(
    clippy::too_many_lines,
    reason = "one construction: chart, edges, pcurves, orientation"
)]
fn side_face(
    model: &mut Model,
    sheet: &Sheet,
    ei: usize,
    (lo, hi): (f64, f64),
    (lower, upper): (&Layer, &Layer),
    risers: &mut HashMap<TShapeId, Shape>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let target = tol.approximation();
    let edge = &sheet.edges[ei];
    let used = &edge.uses[0];
    let face = &sheet.faces[used.face];
    let Some((curve, range)) = edge.curve.clone() else {
        ogeom_bail!(Construction, "a free edge has no curve in space");
    };
    let (t0, t1) = range;
    let height = hi - lo;
    // The sheet's normal along the edge, from the face it bounds.
    let normal_at = |t: f64| -> OgeomResult<(Point, Vector)> {
        let at = curve.point_at(t, tol)?;
        let uv = seat_on(sheet, edge, used, (at, t, range), tol)?;
        let Some(n) = sheet_normal(face, uv, tol) else {
            ogeom_bail!(
                Construction,
                "a free edge runs into a point where its face has no \
                 normal; the side there is not determined"
            );
        };
        Ok((at, n))
    };
    let ruled = |t: f64, s: f64| -> OgeomResult<Point> {
        let (at, n) = normal_at(t)?;
        Ok(at + n * (lo + s))
    };

    // The chart: u is the edge's own parameter, v = v0 + sense * s for the
    // height s above the lower layer.
    let measure = |surface: &SurfaceGeometry, v0: f64, sense: f64| -> OgeomResult<f64> {
        let mut worst = 0.0_f64;
        for i in 0..=4 * PROBES {
            let t = t0 + (t1 - t0) * f64::from(i) / f64::from(4 * PROBES);
            for s in [0.0, 0.5 * height, height] {
                let p = surface.point_at(t, sense.mul_add(s, v0), tol)?;
                worst = worst.max(p.distance(ruled(t, s)?));
            }
        }
        Ok(worst)
    };
    let mut chart: Option<(SurfaceGeometry, f64, f64, f64)> = None;
    let (_, n0) = normal_at(t0)?;
    match &curve {
        Curve::Line(line) => {
            let axis = line.axis();
            let across = Direction::new(axis.direction.vector().cross(n0), tol)?;
            let frame = Frame::new(axis.location + n0 * lo, across, axis.direction, tol)?;
            let pad = (t1 - t0).abs() + height;
            let plane: SurfaceGeometry = PlaneSurface::over(
                ogeom_math::Plane::new(frame),
                (t0 - pad, t1 + pad),
                (-pad, height + pad),
            )?
            .into();
            if measure(&plane, 0.0, 1.0)? <= tol.confusion() * 10.0 {
                chart = Some((plane, 0.0, 1.0, 0.0));
            }
        }
        Curve::Circle(c) => {
            let circle = c.circle();
            let axis = circle.frame().z().vector();
            let sense = if n0.dot(axis) >= 0.0 { 1.0 } else { -1.0 };
            let (a, b) = (sense * lo, sense * hi);
            let cylinder: SurfaceGeometry = CylinderSurface::new(
                Cylinder::new(circle.frame(), circle.radius(), tol)?,
                (a.min(b) - height, a.max(b) + height),
            )?
            .into();
            if measure(&cylinder, sense * lo, sense)? <= tol.confusion() * 10.0 {
                chart = Some((cylinder, sense * lo, sense, 0.0));
            }
        }
        _ => {}
    }
    let (surface, v0, sense, deviation) = match chart {
        Some(found) => found,
        None => {
            let mut spans: u32 = 16;
            loop {
                let at = |k: u32| t0 + (t1 - t0) * f64::from(k) / f64::from(2 * spans);
                let us: Vec<f64> = (0..=spans).map(|k| at(2 * k)).collect();
                let rows: Vec<Vec<Point>> = [0.0, height]
                    .iter()
                    .map(|s| us.iter().map(|t| ruled(*t, *s)).collect())
                    .collect::<OgeomResult<_>>()?;
                let fitted = ogeom_geom::fit::fit_surface_grid_at(
                    &us,
                    &[0.0, height],
                    &rows,
                    3,
                    target * 0.5,
                    tol,
                )?;
                let fit: SurfaceGeometry = fitted.curve.into();
                let worst = fitted.error.max(measure(&fit, 0.0, 1.0)?);
                if worst <= target {
                    break (fit, 0.0, 1.0, worst);
                }
                if spans >= 1024 {
                    ogeom_bail!(
                        NotDone,
                        "the side face along a free edge did not fit within \
                         {target}; the closest fit was {worst} off"
                    );
                }
                spans *= 2;
            }
        }
    };
    let side_id = model.geometry_mut().add_surface(surface.clone());

    // The risers: one straight edge per free vertex, lower to upper.
    let mut riser = |model: &mut Model, vertex: TShapeId| -> OgeomResult<Shape> {
        if let Some(found) = risers.get(&vertex) {
            return Ok(found.clone());
        }
        let (from, a) = lower.vertices[&vertex].clone();
        let (to, b) = upper.vertices[&vertex].clone();
        let segment = LineCurve::segment(a, b, tol)?;
        let span = a.distance(b);
        let shape = make_edge_between(model, segment.into(), (0.0, span), &from, &to, tol)?.shape;
        risers.insert(vertex, shape.clone());
        Ok(shape)
    };
    let rise_start = riser(model, edge.ends.0)?;
    let rise_end = riser(model, edge.ends.1)?;
    let low = lower.edges[&edge.node].clone();
    let high = upper.edges[&edge.node].clone();

    // The rectangle (t0, 0) .. (t1, height) walked counter-clockwise in
    // (t, s): along the lower edge, up the end riser, back along the upper
    // edge, down the start riser. Where v runs against s the chart sees it
    // clockwise, and the walk turns round.
    let forward = Orientation::Forward;
    let backward = Orientation::Reversed;
    let walk: Vec<(Shape, Orientation)> = if sense > 0.0 {
        vec![
            (low.clone(), forward),
            (rise_end.clone(), forward),
            (high.clone(), backward),
            (rise_start.clone(), backward),
        ]
    } else {
        vec![
            (rise_start.clone(), forward),
            (high.clone(), forward),
            (rise_end.clone(), backward),
            (low.clone(), backward),
        ]
    };
    let ring: Vec<Shape> = walk.iter().map(|(e, o)| e.oriented(*o)).collect();
    let wire = make_wire(model, &ring, tol)?.shape;
    let bare = make_face_on(model, side_id, std::slice::from_ref(&wire), tol)?.shape;

    let row = |v: f64| -> OgeomResult<PlanarCurve> {
        Ok(Line2d::over(Axis2::new(Point2::new(0.0, v), Direction2::X), t0, t1)?.into())
    };
    let up = Direction2::new(ogeom_math::Vector2::new(0.0, sense), tol)?;
    let column = |t: f64, length: f64| -> OgeomResult<PlanarCurve> {
        Ok(Line2d::over(Axis2::new(Point2::new(t, v0), up), 0.0, length)?.into())
    };
    let identity = Location::identity();
    attach_pcurve(model, &low, row(v0)?, side_id, identity.clone(), range)?;
    attach_pcurve(
        model,
        &high,
        row(sense.mul_add(height, v0))?,
        side_id,
        identity.clone(),
        range,
    )?;
    let length_of = |model: &Model, e: &Shape| -> OgeomResult<f64> {
        match model
            .node(e)
            .and_then(|n| n.data().as_edge())
            .and_then(EdgeData::curve3d)
        {
            Some(EdgeRepr::Curve3d { range, .. }) => Ok(range.1),
            _ => ogeom_bail!(Construction, "a riser has no curve"),
        }
    };
    if rise_start.node() == rise_end.node() {
        // A closed free edge: one riser bounds the side twice, a seam. Its
        // forward use is whichever end the walk climbs.
        let length = length_of(model, &rise_start)?;
        let (climbed, descended) = if sense > 0.0 { (t1, t0) } else { (t0, t1) };
        attach_seam(
            model,
            &rise_start,
            column(climbed, length)?,
            column(descended, length)?,
            side_id,
            identity,
            (0.0, length),
        )?;
    } else {
        let length = length_of(model, &rise_start)?;
        attach_pcurve(
            model,
            &rise_start,
            column(t0, length)?,
            side_id,
            identity.clone(),
            (0.0, length),
        )?;
        let length = length_of(model, &rise_end)?;
        attach_pcurve(
            model,
            &rise_end,
            column(t1, length)?,
            side_id,
            identity,
            (0.0, length),
        )?;
    }
    for e in [&low, &high, &rise_start, &rise_end] {
        same_parameter(model, e, tol)?;
    }
    if deviation > 0.0 {
        model.widen(&bare, Tolerance::new(deviation + tol.confusion())?)?;
    }

    // Out of the material: the face lies to the left of the edge as the
    // bare face walks it, so outward is the walk's tangent across the
    // surface's own normal.
    let mid = f64::midpoint(t0, t1);
    let uv = seat_on(
        sheet,
        edge,
        used,
        (curve.point_at(mid, tol)?, mid, range),
        tol,
    )?;
    let natural = face.surface.normal_at(uv.x, uv.y, tol)?.vector();
    let walked = if used.sense == Orientation::Reversed {
        -curve.d1_at(mid, tol)?
    } else {
        curve.d1_at(mid, tol)?
    };
    let outward = walked.cross(natural);
    let side_normal = surface
        .normal_at(mid, sense.mul_add(0.5 * height, v0), tol)?
        .vector();
    Ok(if side_normal.dot(outward) >= 0.0 {
        bare
    } else {
        bare.reversed()
    })
}
