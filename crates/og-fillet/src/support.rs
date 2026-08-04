//! Shared scaffolding for the subtractive blends.
//!
//! A chamfer and a constant-radius fillet on a straight edge between planar
//! faces differ only in the face that replaces the edge — a bevel plane for
//! one, a tangent cylinder for the other. Everything around that face — the
//! seat on the solid, the legs running along the adjacent faces, faces
//! assembled from explicit curves with exact pcurves — is one piece of
//! scaffolding, kept here so the two operations cannot drift apart.

use og_algo::{attach_pcurve, make_edge_between, make_face, make_vertex, make_wire};
use og_core::{OgResult, Tolerances, og_bail};
use og_geom::Curve3d as _;
use og_geom::{Curve, LineCurve, PlaneSurface, SurfaceGeometry};
use og_math::{Direction, Plane, Point, Vector};
use og_topo::{EdgeRepr, Filter, Location, Model, NodeData, Shape, ShapeType, explore};

/// Where a blend sits on a solid: a straight edge and the two planar faces
/// meeting there, reduced to the numbers the wedge construction runs on.
pub(crate) struct Seat {
    /// The edge's start, at its lower parameter.
    pub start: Point,
    /// The edge's end.
    pub end: Point,
    /// Unit direction from start to end.
    pub along: Vector,
    /// Outward unit normals of the two faces, in discovery order.
    pub normals: [Vector; 2],
    /// The two faces themselves, in the same order — so a caller naming a
    /// face can find which leg it owns.
    pub faces: [Shape; 2],
}

impl Seat {
    /// On face `i`, the unit direction perpendicular to the edge that walks
    /// away from the *other* face — into the material the blend cuts back
    /// along.
    pub fn leg(&self, i: usize, tol: Tolerances) -> OgResult<Vector> {
        let own = self.normals[i];
        let other = self.normals[1 - i];
        let mut t = own.cross(self.along);
        if t.dot(other) > 0.0 {
            t = -t;
        }
        let m = t.magnitude();
        if m <= tol.confusion() {
            og_bail!(Construction, "a face is tangent to its own edge");
        }
        Ok(t / m)
    }
}

/// An edge's 3D curve and range, cloned out of the model.
pub(crate) fn edge_curve(model: &Model, edge: &Shape) -> OgResult<(Curve, (f64, f64))> {
    let Some(node) = model.node(edge) else {
        og_bail!(Dangling, "edge is not in this model");
    };
    let NodeData::Edge(data) = node.data() else {
        og_bail!(Construction, "expected an edge");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        og_bail!(Construction, "the edge has no curve to blend along");
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        og_bail!(Dangling, "curve is not in this model");
    };
    Ok((geometry.clone(), *range))
}

/// Find the seat of a blend: the straight edge's ends and direction, and the
/// outward normals of the exactly two planar faces of `solid` meeting there.
///
/// Refuses concave and tangent edges: the wedge these blends subtract lies in
/// the material only when the edge is convex. A concave blend *adds* material
/// and is a different construction, recorded in the deferred table.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the edge is
/// not straight, is not shared by exactly two planar faces of `solid`, or the
/// edge is concave or tangent.
pub(crate) fn planar_seat(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> OgResult<Seat> {
    let (curve, range) = edge_curve(model, edge)?;
    let Curve::Line(_) = &curve else {
        og_bail!(
            Construction,
            "blending a curved edge needs the marching blend machinery; this \
             is the straight-edge form"
        );
    };
    let start = curve.point_at(range.0, tol)?;
    let end = curve.point_at(range.1, tol)?;
    let along = (end - start) / start.distance(end);

    let mut normals: Vec<Vector> = Vec::new();
    let mut faces: Vec<Shape> = Vec::new();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let touches = explore(model, &face, Filter::OfType(ShapeType::Edge))?
            .iter()
            .any(|e| e.node() == edge.node());
        if !touches {
            continue;
        }
        let Some(node) = model.node(&face) else {
            og_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            og_bail!(Construction, "face node holds no face data");
        };
        let Some(SurfaceGeometry::Plane(plane)) = model.geometry().surface(data.surface) else {
            og_bail!(
                Construction,
                "blending an edge of a curved face needs the marching blend \
                 machinery; this is the planar form"
            );
        };
        let placement = face.transform(model.datums())?;
        let mut normal = placement.apply_vector(plane.plane().normal().vector());
        if face.orientation() == og_topo::Orientation::Reversed {
            normal = -normal;
        }
        normals.push(normal);
        faces.push(face);
    }
    if normals.len() != 2 {
        og_bail!(
            Construction,
            "a blend needs an edge shared by exactly two faces, found {}",
            normals.len()
        );
    }

    let seat = Seat {
        start,
        end,
        along,
        normals: [normals[0], normals[1]],
        faces: [faces[0].clone(), faces[1].clone()],
    };
    // Convexity: on a convex edge, walking along one face away from the other
    // moves *behind* the other face's plane. On a concave edge it moves in
    // front, and the wedge these blends subtract would sit in empty space.
    if seat.leg(0, tol)?.dot(seat.normals[1]) > -tol.angular() {
        og_bail!(
            Construction,
            "the edge is concave or its faces are tangent; the subtractive \
             blend needs a convex edge"
        );
    }
    Ok(seat)
}

/// An edge along the segment from `from` to `to`, parameterized by arc
/// length, joining two existing vertices.
///
/// A wire chains through shared vertex *objects*, not through coincident
/// coordinates — which is why this takes the vertices and not just the
/// points.
pub(crate) fn segment_between(
    model: &mut Model,
    from: (&Shape, Point),
    to: (&Shape, Point),
    tol: Tolerances,
) -> OgResult<Shape> {
    let line = LineCurve::segment(from.1, to.1, tol)?;
    let curve = Curve::Line(line);
    let domain = curve.domain();
    Ok(make_edge_between(model, curve, domain, from.0, to.0, tol)?.shape)
}

/// A face on `surface` bounded by `edges` in traversal order, with an exact
/// same-parameter pcurve attached to every edge.
///
/// The edges' curves must lie on the surface in a configuration
/// [`og_intersect::exact_pcurve_of`] recognises. That is the point: a blend's
/// faces are built from curves *chosen* to have closed-form charts, and a
/// fitted pcurve here would manufacture disagreement where none exists.
pub(crate) fn face_from_edges(
    model: &mut Model,
    surface: SurfaceGeometry,
    edges: &[Shape],
    tol: Tolerances,
) -> OgResult<Shape> {
    let wire = make_wire(model, edges, tol)?.shape;
    let face = make_face(model, surface.clone(), std::slice::from_ref(&wire), tol)?.shape;
    let surface_id = {
        let Some(node) = model.node(&face) else {
            og_bail!(Dangling, "the face just built is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            og_bail!(Construction, "the face holds no face data");
        };
        data.surface
    };
    for pedge in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
        let (curve, prange) = {
            let Some(node) = model.node(&pedge) else {
                og_bail!(Dangling, "edge is not in this model");
            };
            let Some(data) = node.data().as_edge() else {
                og_bail!(Construction, "edge node holds no edge data");
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                og_bail!(Construction, "a blend face edge has no 3D curve");
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                og_bail!(Dangling, "curve is not in this model");
            };
            (geometry.clone(), *range)
        };
        let Some(pcurve) = og_intersect::exact_pcurve_of(&curve, &surface, tol) else {
            og_bail!(
                Construction,
                "a blend face edge has no closed-form pcurve on its surface"
            );
        };
        attach_pcurve(
            model,
            &pedge,
            pcurve,
            surface_id,
            Location::identity(),
            prange,
        )?;
    }
    Ok(face)
}

/// Sew the wedge's faces, demand a closed shell, subtract it from the solid,
/// and report the history: the blend is a cut, and the edge it replaces is
/// gone.
pub(crate) fn subtract_wedge(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    faces: &[Shape],
    tol: Tolerances,
) -> OgResult<og_algo::Built> {
    let sewn = og_algo::sew(model, faces, tol)?;
    if sewn.shells.len() != 1 || !og_algo::is_shell_closed(model, &sewn.shells[0])? {
        og_bail!(Construction, "the blend wedge did not close");
    }
    let wedge = og_algo::make_solid(model, std::slice::from_ref(&sewn.shells[0]))?;
    let mut result = og_bool::cut(model, solid, &wedge.shape, tol)?;
    result.history.delete(edge);
    Ok(result)
}

/// A planar face over explicit corners, with `outward` as its plane normal.
///
/// The corners must be coplanar and `outward` perpendicular to them — the
/// callers know both by construction, which is why this takes the normal
/// instead of rediscovering it.
pub(crate) fn planar_face(
    model: &mut Model,
    corners: &[Point],
    outward: Vector,
    tol: Tolerances,
) -> OgResult<Shape> {
    let normal = Direction::new(outward, tol)?;
    let plane = Plane::through(corners[0], normal);
    let mut reach = 1.0_f64;
    for p in corners {
        reach = reach.max(p.distance(corners[0]) * 2.0);
    }
    let surface = PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?;
    let vertices: Vec<Shape> = corners
        .iter()
        .map(|p| make_vertex(model, *p).shape)
        .collect();
    let mut edges = Vec::with_capacity(corners.len());
    for i in 0..corners.len() {
        let j = (i + 1) % corners.len();
        edges.push(segment_between(
            model,
            (&vertices[i], corners[i]),
            (&vertices[j], corners[j]),
            tol,
        )?);
    }
    face_from_edges(model, surface.into(), &edges, tol)
}
