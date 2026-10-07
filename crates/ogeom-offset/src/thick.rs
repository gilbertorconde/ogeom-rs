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
//! within the approximation tolerance. Where the finest fit still misses
//! it (the offset of a fitted surface is only as smooth as that surface's
//! normal), the closest fit is kept if it is within the fit's bound: the
//! larger of the face's own tolerance and a ten-thousandth of the
//! distance. The measured deviation widens the moved face's tolerance,
//! which is where it is reported. Because the chart
//! is kept, every pcurve of the sheet carries over unchanged, so a trimmed
//! face, holes and all, is trimmed along the same chart curves as its original.
//!
//! Edges and vertices move the same way: each point along the normal of
//! the faces it bounds. A line or a circle whose move is a translation or
//! a similarity stays a line or a circle (checked at samples, not assumed);
//! any other edge is fitted at its own parameters the same way, its bound
//! the largest of its faces'.
//! Where two faces meet, both must move the shared boundary to the same
//! place: faces meeting tangentially do, faces meeting at a crease do not,
//! and an offset sheet refuses a crease by name rather than tearing or
//! patching it.
//!
//! Thickening builds the sheet moved to both of its sides (or the sheet
//! itself and one side) and closes the gap along every free edge with the
//! ruled face between the edge and its offset: a plane along a straight
//! edge of a constant normal, a cylinder along a circle whose normal runs
//! along its axis, and otherwise a B-spline fitted to the rulings on the
//! edge's own parameter, its deviation recorded the same way. The pieces
//! share their edges by construction, so the shell closes without sewing.
//!
//! A thickened sheet joins its faces across a crease with a mitre, as a
//! solid's walls meet when it is hollowed: in each layer the two faces'
//! moved surfaces are cut back or run on to where they cross, found beside
//! each point of the crease in the plane square to it there. A crease is
//! joined where it runs from one border of the sheet to another, each end
//! a vertex where it meets one free edge of each of its two faces, and
//! where those borders and their offsets lie square to it (as the walls of
//! a profile swept square to its plane do), so the side along each such
//! border is flat and reaches the mitre. The mitred edges are lines and
//! circles on analytic faces, each charted exactly on its moved face.
//!
//! Refused by name: a distance beyond a face's smallest radius of curvature
//! on the side it moves to (checked over the face's chart window), a
//! free-form face whose normal is undefined somewhere the offset needs it,
//! a crease in an offset sheet and a crease a thickened sheet cannot mitre,
//! an edge shared by more than two faces, and a thickened sheet whose
//! layers run into each other.

use ogeom_algo::{
    Built, History, attach_pcurve, attach_seam, edge_vertices, make_edge_between, make_face_on,
    make_face_with_pcurves, make_shell, make_solid, make_vertex, make_wire,
};
use ogeom_core::{OgeomResult, Tolerance, Tolerances, ogeom_bail};
use ogeom_geom::{
    Curve, Curve2d as _, Curve3d as _, CylinderSurface, Line2d, LineCurve, OffsetSurface,
    PlanarCurve, PlaneSurface, Surface as _, SurfaceGeometry, Transformable as _,
};
use ogeom_math::{
    Axis2, Cylinder, Direction, Direction2, Frame, Point, Point2, Transform, Transform2, Vector,
    Vector2,
};
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
/// approximation tolerance where refining the fit reaches it, and
/// otherwise within the larger of the face's own tolerance and a
/// ten-thousandth of `distance`; its tolerance (and so its edges' and
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
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) if the finest
/// fit misses that bound.
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
    let layer = moved_layer(model, &read, distance, None, tol)?;
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
/// face is fitted to the same bound and reports its deviation in its
/// tolerance. Where two faces meet at a crease, their offsets meet on the
/// mitre: each layer's faces are cut back or run on to where their moved
/// surfaces cross, and the borders ending at the crease with them. Each
/// side face is the ruled face between a free edge's two images: a plane
/// or a cylinder where the rulings make one on the edge's own parameter,
/// and otherwise a B-spline fitted to the rulings to the same bound (the
/// thickness standing for the distance), its deviation recorded the same
/// way; along a border ending at a crease it is the flat face between the
/// images.
///
/// # Errors
///
/// As [`offset_sheet`] but for creases, and if the sheet is closed or has
/// no free edge (a closed shell bounds a solid already; hollow it with
/// [`make_thick_solid`](crate::make_thick_solid)), a free edge runs into a
/// point where its face has no normal, or the layers run into each other.
/// A crease is refused where it closes on itself, ends where more than
/// one border of each of its faces meets it, is walked the same way by
/// both faces, folds its faces back onto each other, or where a border
/// ending at it is not flat with its offset or does not reach the mitre,
/// or its mitred edges have no exact chart image on the moved faces.
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
    let creases = creases(&read, lo.abs().max(hi.abs()), tol)?;
    let lower = moved_layer(model, &read, lo, Some(&creases), tol)?;
    let upper = moved_layer(model, &read, hi, Some(&creases), tol)?;

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
        let side = if creases.touches(&read.edges[ei]) {
            flat_side_face(model, &read, ei, (&lower, &upper), &mut risers, tol)?
        } else {
            side_face(
                model,
                &read,
                ei,
                (lo, hi),
                (&lower, &upper),
                &mut risers,
                tol,
            )?
        };
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
    // Crossing layers are asked about first: where they cross, material
    // lies on both sides of a face, so a solid that runs into itself also
    // shows faces turned in, and the crossing is the cause.
    if !ogeom_algo::check_self_intersection(model, &solid, tol)?.is_empty() {
        ogeom_bail!(
            Construction,
            "the thickened sheet runs into itself; the thickness is larger \
             than the room between its faces"
        );
    }
    if !ogeom_algo::inside_out_faces(model, &solid, tol)?.is_empty() {
        ogeom_bail!(
            Construction,
            "the thickened sheet folds through itself and turns faces inside out"
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
    /// The face's own tolerance.
    tolerance: f64,
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
            tolerance: data.tolerance.get(),
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
                // Not past a pole, a side of the box that is one point: a
                // chart carried on through the axis of a revolution turns
                // its normal over there, and no offset follows it.
                let pole =
                    |u: Option<f64>, v: Option<f64>| collapsed(&face.surface, u, v, (lo, hi), tol);
                let margin = |m: f64, side: bool| if side { 0.0 } else { m };
                let (mut u0, mut u1, mut v0, mut v1) = (
                    lo.x - margin(mu, pole(Some(lo.x), None)?),
                    hi.x + margin(mu, pole(Some(hi.x), None)?),
                    lo.y - margin(mv, pole(None, Some(lo.y))?),
                    hi.y + margin(mv, pole(None, Some(hi.y))?),
                );
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
///
/// With `creases`, the faces either side of each crease meet on the
/// mitre: the crease's image is where their moved surfaces cross, beside
/// each point of the crease in the plane square to it there, and the
/// borders ending at the crease are cut back or carried on to that line.
/// Without, a crease is refused.
fn moved_layer(
    model: &mut Model,
    sheet: &Sheet,
    distance: f64,
    creases: Option<&Creases>,
    tol: Tolerances,
) -> OgeomResult<Layer> {
    let creases = creases.filter(|c| distance != 0.0 && !c.edges.is_empty());
    let target = tol.approximation();
    let mut history = History::new();

    let mut surfaces: Vec<MovedSurface> = Vec::with_capacity(sheet.faces.len());
    for face in &sheet.faces {
        let (geometry, exact, deviation) = if distance == 0.0 {
            (face.surface.clone(), true, 0.0)
        } else {
            let bound = fit_bound(face.tolerance, distance, tol);
            moved_surface(face, distance, (target, bound), tol)?
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
        if let Some(&(ci, at_end)) = creases.and_then(|c| c.ends.get(node)) {
            let edge = &sheet.edges[ci];
            let Some((_, range)) = &edge.curve else {
                ogeom_bail!(Construction, "a crease has no curve in space");
            };
            let t = if at_end { range.1 } else { range.0 };
            let moved = mitre_at(sheet, &surfaces, edge, t, tol)?;
            let shape = make_vertex(model, moved).shape;
            history.modify(&Shape::of(*node), shape.clone());
            vertex_slop.insert(*node, 0.0);
            vertices.insert(*node, (shape, moved));
            continue;
        }
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
    for (ei, edge) in sheet.edges.iter().enumerate() {
        let is_crease = creases.is_some_and(|c| c.edges.contains(&ei));
        let cut = creases.map_or([false, false], |c| {
            [
                c.ends.contains_key(&edge.ends.0),
                c.ends.contains_key(&edge.ends.1),
            ]
        });
        let mut recharted: Option<(Curve, (f64, f64))> = None;
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
            let mitre = |t: f64| mitre_at(sheet, &surfaces, edge, t, tol);
            let off: &dyn Fn(f64) -> OgeomResult<Point> = if is_crease { &mitre } else { &off };
            let (moved, slop) = if distance == 0.0 {
                (curve.clone(), 0.0)
            } else if let Some(exact) = exact_moved_curve(curve, *range, off, tol)? {
                (exact, 0.0)
            } else {
                let bound = edge
                    .uses
                    .iter()
                    .map(|u| fit_bound(sheet.faces[u.face].tolerance, distance, tol))
                    .fold(target, f64::max);
                fitted_moved_curve(*range, off, (target, bound), tol)?
            };
            // A border ending at a crease runs to the crease's image.
            let (moved, range) = if !is_crease && cut.contains(&true) {
                let ends = [cut[0].then_some(from_at), cut[1].then_some(to_at)];
                retrimmed(&moved, *range, ends, tol)?
            } else {
                (moved, *range)
            };
            if is_crease || cut.contains(&true) {
                recharted = Some((moved.clone(), range));
            }
            // The curve's ends must reach the vertices: within the slop of
            // the fit and the spread the vertex's faces agreed within.
            for (t, (vertex, at)) in [(range.0, (&from, from_at)), (range.1, (&to, to_at))] {
                let gap = moved.point_at(t, tol)?.distance(at);
                if gap > tol.confusion() {
                    model.widen(vertex, Tolerance::new(gap + tol.confusion())?)?;
                }
            }
            let shape = make_edge_between(model, moved, range, &from, &to, tol)?.shape;
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
        // An edge the mitre moved is charted afresh on each moved surface.
        if let Some((curve, range)) = &recharted {
            for used in &edge.uses {
                let surface = &surfaces[used.face].geometry;
                let Some(image) = ogeom_intersect::exact_pcurve_over(curve, *range, surface, tol)
                else {
                    ogeom_bail!(
                        Construction,
                        "the sheet's faces meet at a crease whose mitre has no exact \
                         image on a moved face; a crease is joined between planes, \
                         cylinders, cones, spheres and tori bounded by lines and circles"
                    );
                };
                let near = sheet
                    .vertices
                    .iter()
                    .find(|(node, _, _)| *node == edge.ends.0)
                    .and_then(|(_, _, seats)| seats.iter().find(|(f, _)| *f == used.face))
                    .map(|(_, uv)| *uv);
                let image = match near {
                    Some(near) => on_branch(image, range.0, near, surface, tol)?,
                    None => image,
                };
                let sid = surfaces[used.face].id;
                attach_pcurve(model, &built, image, sid, Location::identity(), *range)?;
            }
            model.widen(&built, edge.data.tolerance)?;
            same_parameter(model, &built, tol)?;
            history.modify(&Shape::of(edge.node), built.clone());
            edges.insert(edge.node, built);
            continue;
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
        let moved = at + n * distance;
        // On a parallel collapsed to a point the chart's normal is any
        // direction, and may point the moved point away from the exact
        // moved surface by as much as twice the distance; that surface
        // answers there instead.
        if surfaces[face].exact {
            let exact = surfaces[face].geometry.point_at(uv.x, uv.y, tol)?;
            if exact.distance(moved) > distance.abs() {
                return Ok(exact);
            }
        }
        return Ok(moved);
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

/// Whether the line of `surface`'s chart at `u` (or at `v`) across the box
/// from `lo` to `hi` is one point within the confusion distance: a pole,
/// or a parallel so near the axis of a revolution that it is one.
fn collapsed(
    surface: &SurfaceGeometry,
    u: Option<f64>,
    v: Option<f64>,
    (lo, hi): (Point2, Point2),
    tol: Tolerances,
) -> OgeomResult<bool> {
    let mut points = Vec::with_capacity(9);
    for k in 0..=8 {
        let t = f64::from(k) / 8.0;
        points.push(surface.point_at(
            u.unwrap_or(lo.x + (hi.x - lo.x) * t),
            v.unwrap_or(lo.y + (hi.y - lo.y) * t),
            tol,
        )?);
    }
    let first = points[0];
    let sum = points
        .iter()
        .fold(Vector::ZERO, |acc, p| acc + (*p - first));
    let centre = first + sum / 9.0;
    Ok(points.iter().all(|p| p.distance(centre) <= tol.confusion()))
}

/// The sample count along a window or a range when measuring.
const PROBES: i32 = 8;

/// The share of the offset distance a fit that cannot reach the
/// approximation tolerance may miss by.
const FIT_SHARE: f64 = 1e-4;

/// The most a fitted offset of a face may miss the exact offset by when
/// refining cannot bring it within the approximation tolerance: the face's
/// own tolerance or a share of the distance, whichever is larger.
fn fit_bound(face_tolerance: f64, distance: f64, tol: Tolerances) -> f64 {
    tol.approximation()
        .max(face_tolerance)
        .max(FIT_SHARE * distance.abs())
}

/// A face's surface moved `distance` along the sheet's normal: the exact
/// parallel where one is measured to keep the chart, else a fit refined
/// toward `target` and kept within `bound`. Returns the geometry, whether
/// it is exact, and the fit's measured deviation.
fn moved_surface(
    face: &SheetFace,
    distance: f64,
    (target, bound): (f64, f64),
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
        // On a pole the basis normal is any direction the chart happens to
        // give, and says nothing of the offset: the window's rows and
        // columns that are one point are not measured.
        let window = (Point2::new(u0, v0), Point2::new(u1, v1));
        let mut pole_columns = Vec::new();
        let mut pole_rows = Vec::new();
        for k in 0..=PROBES {
            let (u, v) = probe(k, k);
            pole_columns.push(collapsed(base, Some(u), None, window, tol)?);
            pole_rows.push(collapsed(base, None, Some(v), window, tol)?);
        }
        let mut worst = 0.0_f64;
        for i in 0..=PROBES {
            for j in 0..=PROBES {
                let (u, v) = probe(i, j);
                if pole_columns[i as usize] || pole_rows[j as usize] {
                    continue;
                }
                if let Some(want) = moved_at(u, v)? {
                    worst = worst.max(parallel.point_at(u, v, tol)?.distance(want));
                }
            }
        }
        if worst <= tol.confusion() * 10.0 {
            return Ok((parallel, true, 0.0));
        }
    }

    let corners = chart_corners(base, face.window);

    // A fit at the chart's own parameters, refined until the points it was
    // not fitted to (the cell centres and edge middles) agree too, or the
    // finest net is reached and the closest fit stands against the bound.
    let mut spans: u32 = 8;
    let mut best: Option<(SurfaceGeometry, f64)> = None;
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
        // The chart's corner lines stand in the fit as knots it may turn
        // at, once the net is fine enough to hold them.
        let held = |kept: &[(f64, usize)]| 4 + 3 * kept.len() <= n as usize + 1;
        let kept_u: &[(f64, usize)] = if held(&corners.0) { &corners.0 } else { &[] };
        let kept_v: &[(f64, usize)] = if held(&corners.1) { &corners.1 } else { &[] };
        let fitted = ogeom_geom::fit::fit_surface_grid_at_with_knots(
            &us,
            &vs,
            &rows,
            3,
            target * 0.5,
            (kept_u, kept_v),
            tol,
        )?;
        let mut surface: SurfaceGeometry = fitted.curve.into();
        let mut worst = fitted.error;
        for j in 0..=2 * n {
            for i in 0..=2 * n {
                let want = fine[j as usize][i as usize];
                let (u, v) = (at(i, u0, u1), at(j, v0, v1));
                worst = worst.max(surface.point_at(u, v, tol)?.distance(want));
            }
        }
        if best.as_ref().is_none_or(|(_, b)| worst < *b) {
            best = Some((surface.clone(), worst));
        }
        let finest = spans >= 128;
        if worst > target && finest {
            let Some((closest, reached)) = best.take() else {
                ogeom_bail!(NotDone, "the offset of a free-form face was not fitted");
            };
            if reached > bound {
                ogeom_bail!(
                    NotDone,
                    "the offset of a free-form face did not fit within {bound}; \
                     the closest fit was {reached} off"
                );
            }
            surface = closest;
            worst = reached;
        }
        if worst <= target || finest {
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
        spans *= 2;
    }
}

/// Interior knots, `(parameter, multiplicity)`.
type Knots = Vec<(f64, usize)>;

/// Where a B-spline chart turns: the interior knots, in `u` and in `v`,
/// across which the moved points are only C0, each as a cubic fit's knot
/// of full multiplicity. A point of the surface is as smooth as its basis
/// is at a knot and its normal one order less, so the moved point is C0
/// wherever the surface is C1 or less.
fn chart_corners(
    surface: &SurfaceGeometry,
    ((u0, u1), (v0, v1)): ((f64, f64), (f64, f64)),
) -> (Knots, Knots) {
    let SurfaceGeometry::BSpline(b) = surface else {
        return (Vec::new(), Vec::new());
    };
    let corners = |knots: &ogeom_math::KnotVector, (lo, hi): (f64, f64)| {
        let degree = knots.degree();
        knots
            .distinct()
            .into_iter()
            .filter(|&(at, multiplicity)| at > lo && at < hi && multiplicity + 1 >= degree)
            .map(|(at, _)| (at, 3))
            .collect()
    };
    (
        corners(b.u_knots(), (u0, u1)),
        corners(b.v_knots(), (v0, v1)),
    )
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
/// midpoints between the fitted samples agree within `target`, or the
/// closest fit kept within `bound` once the finest is reached; the curve
/// and the deviation measured.
fn fitted_moved_curve(
    range: (f64, f64),
    off: &dyn Fn(f64) -> OgeomResult<Point>,
    (target, bound): (f64, f64),
    tol: Tolerances,
) -> OgeomResult<(Curve, f64)> {
    let mut spans: u32 = 16;
    let mut best: Option<(Curve, f64)> = None;
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
        if best.as_ref().is_none_or(|(_, b)| worst < *b) {
            best = Some((curve, worst));
        }
        if spans >= 1024 {
            return match best {
                Some((closest, reached)) if reached <= bound => Ok((closest, reached)),
                Some((_, reached)) => ogeom_bail!(
                    NotDone,
                    "the offset of an edge did not fit within {bound}; the \
                     closest fit was {reached} off"
                ),
                None => ogeom_bail!(NotDone, "the offset of an edge was not fitted"),
            };
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

/// Where the sheet's faces meet at an angle, so that their offsets part
/// on one side and cross on the other.
#[derive(Default)]
struct Creases {
    /// The creased edges, by index into the sheet's edges.
    edges: Vec<usize>,
    /// Each vertex a crease ends at: the crease, and whether the vertex is
    /// its end rather than its start.
    ends: HashMap<TShapeId, (usize, bool)>,
}

impl Creases {
    /// Whether an edge ends at a crease.
    fn touches(&self, edge: &SheetEdge) -> bool {
        self.ends.contains_key(&edge.ends.0) || self.ends.contains_key(&edge.ends.1)
    }
}

/// The sheet's creases at an offset reaching `reach`: the edges two faces
/// share whose normals part there by more than the faces could agree on.
/// Each crease must run between two of the sheet's borders: its ends are
/// vertices where it meets one free edge of each of its two faces and
/// nothing else, and its faces do not fold back onto each other.
fn creases(sheet: &Sheet, reach: f64, tol: Tolerances) -> OgeomResult<Creases> {
    let target = tol.approximation();
    let mut found = Creases::default();
    for (ei, edge) in sheet.edges.iter().enumerate() {
        let [a, b] = edge.uses.as_slice() else {
            continue;
        };
        let Some((curve, range)) = &edge.curve else {
            continue;
        };
        if a.face == b.face {
            continue;
        }
        let mut parting = 0.0_f64;
        for k in 0..=PROBES {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / f64::from(PROBES);
            let at = curve.point_at(t, tol)?;
            let na = sheet_normal(
                &sheet.faces[a.face],
                seat_on(sheet, edge, a, (at, t, *range), tol)?,
                tol,
            );
            let nb = sheet_normal(
                &sheet.faces[b.face],
                seat_on(sheet, edge, b, (at, t, *range), tol)?,
                tol,
            );
            if let (Some(na), Some(nb)) = (na, nb) {
                parting = parting.max((na - nb).magnitude());
            }
        }
        if parting * reach / 2.0 > target {
            found.edges.push(ei);
        }
    }
    for &ei in &found.edges {
        let edge = &sheet.edges[ei];
        if edge.ends.0 == edge.ends.1 {
            ogeom_bail!(
                Construction,
                "the sheet's faces meet at a crease that closes on itself; a crease \
                 is joined where it runs from one border of the sheet to another"
            );
        }
        let walked = |u: &EdgeUse| {
            sheet.faces[u.face]
                .occurrence
                .orientation()
                .compose(u.sense)
        };
        if walked(&edge.uses[0]) == walked(&edge.uses[1]) {
            ogeom_bail!(
                Construction,
                "the sheet's faces walk a crease the same way, so they disagree on \
                 which side of the sheet they face; orient the sheet first"
            );
        }
        let mut sides: Vec<usize> = edge.uses.iter().map(|u| u.face).collect();
        sides.sort_unstable();
        for (vertex, at_end) in [(edge.ends.0, false), (edge.ends.1, true)] {
            let meeting: Vec<&SheetEdge> = sheet
                .edges
                .iter()
                .enumerate()
                .filter(|(i, e)| *i != ei && (e.ends.0 == vertex || e.ends.1 == vertex))
                .map(|(_, e)| e)
                .collect();
            let mut faces: Vec<usize> = meeting
                .iter()
                .filter(|e| e.uses.len() == 1)
                .map(|e| e.uses[0].face)
                .collect();
            faces.sort_unstable();
            if meeting.len() != 2 || faces != sides {
                ogeom_bail!(
                    Construction,
                    "the sheet's faces meet at a crease that ends where more than \
                     its two faces' borders meet; a crease is joined where it runs \
                     from one border of the sheet to another"
                );
            }
            found.ends.insert(vertex, (ei, at_end));
        }
    }
    Ok(found)
}

/// The mitre beside the point at `t` on a crease: where the moved
/// surfaces of the crease's two faces cross in the plane square to the
/// crease there.
fn mitre_at(
    sheet: &Sheet,
    surfaces: &[MovedSurface],
    edge: &SheetEdge,
    t: f64,
    tol: Tolerances,
) -> OgeomResult<Point> {
    let Some((curve, range)) = &edge.curve else {
        ogeom_bail!(Construction, "a crease has no curve in space");
    };
    let [a, b] = edge.uses.as_slice() else {
        ogeom_bail!(Construction, "a crease is shared by two faces");
    };
    let at = curve.point_at(t, tol)?;
    let seat_a = seat_on(sheet, edge, a, (at, t, *range), tol)?;
    let seat_b = seat_on(sheet, edge, b, (at, t, *range), tol)?;
    let (Some(na), Some(nb)) = (
        sheet_normal(&sheet.faces[a.face], seat_a, tol),
        sheet_normal(&sheet.faces[b.face], seat_b, tol),
    ) else {
        ogeom_bail!(
            Construction,
            "a crease runs into a point where a face has no normal; the mitre \
             there is not determined"
        );
    };
    if 1.0 + na.dot(nb) <= 1e-6 {
        ogeom_bail!(
            Construction,
            "the sheet folds back onto itself at a crease; the offsets of its \
             two faces never meet"
        );
    }
    let along = curve.d1_at(t, tol)?;
    let speed = along.magnitude();
    if speed <= f64::MIN_POSITIVE {
        ogeom_bail!(Construction, "a crease stands still at {at:?}");
    }
    mitre_point(
        (&surfaces[a.face].geometry, &surfaces[b.face].geometry),
        (seat_a, seat_b),
        at,
        along / speed,
        tol,
    )
}

/// The point on both surfaces in the plane through `at` square to
/// `tangent`, by Newton's method in the two charts from `seats`.
fn mitre_point(
    (first, second): (&SurfaceGeometry, &SurfaceGeometry),
    (seat_a, seat_b): (Point2, Point2),
    at: Point,
    tangent: Vector,
    tol: Tolerances,
) -> OgeomResult<Point> {
    let mut x = [seat_a.x, seat_a.y, seat_b.x, seat_b.y];
    for _ in 0..32 {
        let (pa, au, av) = first.point_d1_at(x[0], x[1], tol)?;
        let (pb, bu, bv) = second.point_d1_at(x[2], x[3], tol)?;
        let gap = pa - pb;
        let across = (pa - at).dot(tangent);
        if gap.magnitude() <= tol.confusion() * 1e-3 && across.abs() <= tol.confusion() * 1e-3 {
            return Ok(pa + (pb - pa) * 0.5);
        }
        let mut rows = [
            [au.x, av.x, -bu.x, -bv.x, -gap.x],
            [au.y, av.y, -bu.y, -bv.y, -gap.y],
            [au.z, av.z, -bu.z, -bv.z, -gap.z],
            [au.dot(tangent), av.dot(tangent), 0.0, 0.0, -across],
        ];
        let Some(step) = solved(&mut rows) else {
            break;
        };
        for (xi, si) in x.iter_mut().zip(step) {
            *xi += si;
        }
    }
    ogeom_bail!(
        NotDone,
        "the offsets of the faces at a crease were not found to meet beside {at:?}"
    )
}

/// The solution of four linear equations, each row its coefficients and
/// right-hand side, by elimination with partial pivoting; `None` where
/// they are singular.
fn solved(rows: &mut [[f64; 5]; 4]) -> Option<[f64; 4]> {
    let scale = rows
        .iter()
        .flat_map(|r| r[..4].iter())
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    if scale <= f64::MIN_POSITIVE {
        return None;
    }
    for col in 0..4 {
        let pivot = (col..4).max_by(|&i, &j| rows[i][col].abs().total_cmp(&rows[j][col].abs()))?;
        if rows[pivot][col].abs() <= 1e-14 * scale {
            return None;
        }
        rows.swap(col, pivot);
        let lead = rows[col];
        for row in rows.iter_mut().skip(col + 1) {
            let f = row[col] / lead[col];
            for (v, l) in row.iter_mut().zip(lead).skip(col) {
                *v -= f * l;
            }
        }
    }
    let mut out = [0.0; 4];
    for r in (0..4).rev() {
        let mut v = rows[r][4];
        for (c, known) in out.iter().enumerate().skip(r + 1) {
            v -= rows[r][c] * known;
        }
        out[r] = v / rows[r][r];
    }
    Some(out)
}

/// A border's moved curve run to new ends: each end given a point is moved
/// to the curve's foot of that point, a line or a circle carried on past
/// its old end where the point lies beyond it.
fn retrimmed(
    curve: &Curve,
    range: (f64, f64),
    ends: [Option<Point>; 2],
    tol: Tolerances,
) -> OgeomResult<(Curve, (f64, f64))> {
    let mut new = [range.0, range.1];
    for (k, end) in ends.iter().enumerate() {
        let Some(p) = end else {
            continue;
        };
        let t = match curve {
            Curve::Line(line) => {
                let axis = line.axis();
                (*p - axis.location).dot(axis.direction.vector())
            }
            _ => {
                let mut t = new[k];
                for _ in 0..64 {
                    let c = curve.point_at(t, tol)?;
                    let d = curve.d1_at(t, tol)?;
                    let speed = d.dot(d);
                    if speed <= f64::MIN_POSITIVE {
                        break;
                    }
                    let step = (*p - c).dot(d) / speed;
                    t += step;
                    if step.abs() <= 1e-15 * (1.0 + t.abs()) {
                        break;
                    }
                }
                t
            }
        };
        let reached = match curve {
            Curve::Line(line) => line.axis().point_at(t),
            _ => curve.point_at(t, tol)?,
        };
        if reached.distance(*p) > tol.confusion() * 10.0 {
            ogeom_bail!(
                Construction,
                "a border ending at a crease does not reach the crease's mitre \
                 ({:.3e} away); a crease is joined where the borders at its ends \
                 lie square to it",
                reached.distance(*p)
            );
        }
        new[k] = t;
    }
    if new[1] - new[0] <= tol.parametric() {
        ogeom_bail!(
            Construction,
            "the thickness reaches past a border's length where it meets a \
             crease; the faces' offsets cross beyond it"
        );
    }
    Ok(match curve {
        Curve::Line(line) => {
            let (d0, d1) = curve.domain();
            let carried = LineCurve::over(line.axis(), new[0].min(d0), new[1].max(d1))?;
            (carried.into(), (new[0], new[1]))
        }
        // Restarted at the new start, so the range runs from zero within
        // one turn, where every chart image of the circle is defined.
        Curve::Circle(c) => {
            let circle = c.circle();
            let start = curve.point_at(new[0], tol)?;
            let frame = Frame::new(
                circle.centre(),
                circle.frame().z(),
                Direction::new(start - circle.centre(), tol)?,
                tol,
            )?;
            let restarted = ogeom_math::Circle::new(frame, circle.radius(), tol)?;
            (
                ogeom_geom::CircleCurve::new(restarted).into(),
                (0.0, new[1] - new[0]),
            )
        }
        _ => (curve.clone(), (new[0], new[1])),
    })
}

/// A pcurve shifted by whole periods of `surface`'s chart so that its
/// point at `t` lies closest to `near`.
fn on_branch(
    image: PlanarCurve,
    t: f64,
    near: Point2,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<PlanarCurve> {
    let ((ua, ub), (va, vb)) = surface.domain();
    let at = image.point_at(t, tol)?;
    let whole = |periodic: bool, period: f64, gap: f64| {
        if periodic && period > 0.0 {
            (gap / period).round() * period
        } else {
            0.0
        }
    };
    let shift = Vector2::new(
        whole(surface.is_periodic_u(), ub - ua, near.x - at.x),
        whole(surface.is_periodic_v(), vb - va, near.y - at.y),
    );
    if shift.x == 0.0 && shift.y == 0.0 {
        return Ok(image);
    }
    image.transformed(&Transform2::translation(shift), tol)
}

/// The straight riser from a free vertex's image in the lower layer to its
/// image in the upper one, built once and shared through `risers`.
fn riser(
    model: &mut Model,
    vertex: TShapeId,
    (lower, upper): (&Layer, &Layer),
    risers: &mut HashMap<TShapeId, Shape>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
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
}

/// The direction out of the material across a free edge, in the sheet:
/// the face lies to the left of the edge as the bare face walks it.
fn outward_of(sheet: &Sheet, ei: usize, tol: Tolerances) -> OgeomResult<Vector> {
    let edge = &sheet.edges[ei];
    let used = &edge.uses[0];
    let face = &sheet.faces[used.face];
    let Some((curve, range)) = &edge.curve else {
        ogeom_bail!(Construction, "a free edge has no curve in space");
    };
    let mid = f64::midpoint(range.0, range.1);
    let uv = seat_on(
        sheet,
        edge,
        used,
        (curve.point_at(mid, tol)?, mid, *range),
        tol,
    )?;
    let natural = face.surface.normal_at(uv.x, uv.y, tol)?.vector();
    let walked = if used.sense == Orientation::Reversed {
        -curve.d1_at(mid, tol)?
    } else {
        curve.d1_at(mid, tol)?
    };
    Ok(walked.cross(natural))
}

/// Points along an edge's curve over its range, start to end.
fn edge_samples(
    model: &Model,
    edge: &Shape,
    count: i32,
    tol: Tolerances,
) -> OgeomResult<Vec<Point>> {
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .and_then(EdgeData::curve3d)
    else {
        ogeom_bail!(Construction, "a layer edge has no curve");
    };
    let Some(curve) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    (0..=count)
        .map(|k| {
            curve.point_at(
                range.0 + (range.1 - range.0) * f64::from(k) / f64::from(count),
                tol,
            )
        })
        .collect()
}

/// The side face along a free edge that ends at a crease: the flat face
/// bounded by the edge's images in the two layers and the risers at its
/// ends, oriented out of the material.
fn flat_side_face(
    model: &mut Model,
    sheet: &Sheet,
    ei: usize,
    (lower, upper): (&Layer, &Layer),
    risers: &mut HashMap<TShapeId, Shape>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let edge = &sheet.edges[ei];
    let rise_start = riser(model, edge.ends.0, (lower, upper), risers, tol)?;
    let rise_end = riser(model, edge.ends.1, (lower, upper), risers, tol)?;
    let low = lower.edges[&edge.node].clone();
    let high = upper.edges[&edge.node].clone();
    // The walk: along the lower image, up the end riser, back along the
    // upper image and down the start riser; its plane is the one the
    // walk turns counter-clockwise about.
    let mut walk = edge_samples(model, &low, 4 * PROBES, tol)?;
    let mut back = edge_samples(model, &high, 4 * PROBES, tol)?;
    back.reverse();
    walk.extend(back);
    let mut turning = Vector::ZERO;
    for (k, p) in walk.iter().enumerate() {
        let q = walk[(k + 1) % walk.len()];
        turning += Vector::new(
            (p.y - q.y) * (p.z + q.z),
            (p.z - q.z) * (p.x + q.x),
            (p.x - q.x) * (p.y + q.y),
        );
    }
    #[allow(clippy::cast_precision_loss, reason = "a few dozen samples")]
    let count = walk.len() as f64;
    let centre = walk
        .iter()
        .fold(Point::ORIGIN, |acc, p| acc + (*p - Point::ORIGIN) / count);
    let normal = Direction::new(turning, tol)?;
    let flat = walk
        .iter()
        .map(|p| (*p - centre).dot(normal.vector()).abs())
        .fold(0.0_f64, f64::max);
    if flat > tol.confusion() * 10.0 {
        ogeom_bail!(
            Construction,
            "the side along a border ending at a crease is not flat ({flat:.3e} \
             out of its plane); a crease is joined where the borders at its ends \
             and their offsets lie in one plane"
        );
    }
    let plane =
        ogeom_math::Plane::new(Frame::new(centre, normal, normal.any_perpendicular(), tol)?);
    let bare = make_face_with_pcurves(
        model,
        PlaneSurface::new(plane).into(),
        &[vec![
            low.clone(),
            rise_end.clone(),
            high.reversed(),
            rise_start.reversed(),
        ]],
        tol,
    )?
    .shape;
    for e in [&low, &high, &rise_start, &rise_end] {
        same_parameter(model, e, tol)?;
    }
    Ok(if normal.vector().dot(outward_of(sheet, ei, tol)?) >= 0.0 {
        bare
    } else {
        bare.reversed()
    })
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
    let bound = fit_bound(face.tolerance, hi - lo, tol);
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
            let mut best: Option<(SurfaceGeometry, f64)> = None;
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
                if best.as_ref().is_none_or(|(_, b)| worst < *b) {
                    best = Some((fit, worst));
                }
                if spans >= 1024 {
                    match best {
                        Some((closest, reached)) if reached <= bound => {
                            break (closest, 0.0, 1.0, reached);
                        }
                        Some((_, reached)) => ogeom_bail!(
                            NotDone,
                            "the side face along a free edge did not fit within \
                             {bound}; the closest fit was {reached} off"
                        ),
                        None => ogeom_bail!(NotDone, "the side face was not fitted"),
                    }
                }
                spans *= 2;
            }
        }
    };
    let side_id = model.geometry_mut().add_surface(surface.clone());

    // The risers: one straight edge per free vertex, lower to upper.
    let rise_start = riser(model, edge.ends.0, (lower, upper), risers, tol)?;
    let rise_end = riser(model, edge.ends.1, (lower, upper), risers, tol)?;
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

    let mid = f64::midpoint(t0, t1);
    let outward = outward_of(sheet, ei, tol)?;
    let side_normal = surface
        .normal_at(mid, sense.mul_add(0.5 * height, v0), tol)?
        .vector();
    Ok(if side_normal.dot(outward) >= 0.0 {
        bare
    } else {
        bare.reversed()
    })
}
