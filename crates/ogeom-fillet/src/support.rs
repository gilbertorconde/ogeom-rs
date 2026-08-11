//! Shared scaffolding for the subtractive blends.
//!
//! A chamfer and a constant-radius fillet on a straight edge between planar
//! faces differ only in the face that replaces the edge — a bevel plane for
//! one, a tangent cylinder for the other. Everything around that face — the
//! seat on the solid, the legs running along the adjacent faces, faces
//! assembled from explicit curves with exact pcurves — is one piece of
//! scaffolding, kept here so the two operations cannot drift apart.

use ogeom_algo::{make_edge_between, make_vertex};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d as _;
use ogeom_geom::{Curve, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Direction, Plane, Point, Vector};
use ogeom_topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, explore};

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
    /// Whether the edge is convex: material inside the dihedral, so a blend
    /// subtracts. A concave edge's blend adds, with every sign mirrored.
    pub convex: bool,
}

impl Seat {
    /// On face `i`, the unit direction perpendicular to the edge that walks
    /// away from the *other* face — into the material the blend cuts back
    /// along.
    pub fn leg(&self, i: usize, tol: Tolerances) -> OgeomResult<Vector> {
        let own = self.normals[i];
        let other = self.normals[1 - i];
        let mut t = own.cross(self.along);
        if t.dot(other) > 0.0 {
            t = -t;
        }
        let m = t.magnitude();
        if m <= tol.confusion() {
            ogeom_bail!(Construction, "a face is tangent to its own edge");
        }
        Ok(t / m)
    }
}

/// An edge's 3D curve and range, cloned out of the model.
pub(crate) fn edge_curve(model: &Model, edge: &Shape) -> OgeomResult<(Curve, (f64, f64))> {
    let Some(node) = model.node(edge) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let NodeData::Edge(data) = node.data() else {
        ogeom_bail!(Construction, "expected an edge");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        ogeom_bail!(Construction, "the edge has no curve to blend along");
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    Ok((geometry.clone(), *range))
}

/// Find the seat of a blend: the straight edge's ends and direction, and the
/// outward normals of the exactly two planar faces of `solid` meeting there.
///
/// Refuses concave and tangent edges: the wedge these blends subtract lies in
/// the material only when the edge is convex. A concave blend *adds* material
/// and is a different construction — docs/PARITY.md, fillet.edge-blends.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the edge is
/// not straight, is not shared by exactly two planar faces of `solid`, or the
/// edge is concave or tangent.
pub(crate) fn planar_seat(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Seat> {
    let (curve, range) = edge_curve(model, edge)?;
    let Curve::Line(_) = &curve else {
        ogeom_bail!(
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
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(SurfaceGeometry::Plane(plane)) = model.geometry().surface(data.surface) else {
            ogeom_bail!(
                Construction,
                "blending an edge of a curved face needs the marching blend \
                 machinery; this is the planar form"
            );
        };
        let placement = face.transform(model.datums())?;
        let mut normal = placement.apply_vector(plane.plane().normal().vector());
        if face.orientation() == ogeom_topo::Orientation::Reversed {
            normal = -normal;
        }
        normals.push(normal);
        faces.push(face);
    }
    if normals.len() != 2 {
        ogeom_bail!(
            Construction,
            "a blend needs an edge shared by exactly two faces, found {}",
            normals.len()
        );
    }

    let mut seat = Seat {
        start,
        end,
        along,
        normals: [normals[0], normals[1]],
        faces: [faces[0].clone(), faces[1].clone()],
        convex: true,
    };
    // Convexity is read from the face itself, not derived from the normals:
    // the leg construction cannot answer it, because it *chooses* its side.
    // Which way the first face actually extends from the edge — sampled
    // against its own trim — leans behind the other face's plane on a convex
    // edge and in front of it on a concave one.
    let raw = {
        let t = normals[0].cross(along);
        let m = t.magnitude();
        if m <= tol.angular() {
            ogeom_bail!(Construction, "a face is tangent to its own edge");
        }
        t / m
    };
    let mid = curve.point_at(f64::midpoint(range.0, range.1), tol)?;
    let mut face_side: Option<Vector> = None;
    'scales: for scale in [1e-3, 1e-2, 5e-2] {
        let eps = start.distance(end) * scale;
        let deflection = ogeom_mesh::Deflection {
            chord: eps * 0.1,
            ..ogeom_mesh::Deflection::default()
        };
        for dir in [raw, -raw] {
            if crate::support::on_face_side(
                model,
                &seat.faces[0],
                mid + dir * eps,
                deflection,
                tol,
            )? {
                face_side = Some(dir);
                break 'scales;
            }
        }
    }
    let Some(extends) = face_side else {
        ogeom_bail!(
            Construction,
            "cannot read which way the edge's face extends; the face is \
             thinner than the probe can resolve"
        );
    };
    let lean = extends.dot(seat.normals[1]);
    if lean.abs() <= tol.angular() {
        ogeom_bail!(
            Construction,
            "the edge's faces are tangent; there is no corner to blend"
        );
    }
    seat.convex = lean < 0.0;
    Ok(seat)
}

/// Whether `probe` sits strictly inside the face's trim.
pub(crate) fn on_face_side(
    model: &Model,
    face: &Shape,
    probe: ogeom_math::Point,
    deflection: ogeom_mesh::Deflection,
    tol: Tolerances,
) -> OgeomResult<bool> {
    Ok(
        ogeom_algo::classify_on_face(model, face, probe, deflection, tol)?
            == ogeom_algo::Containment::In,
    )
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
) -> OgeomResult<Shape> {
    let line = LineCurve::segment(from.1, to.1, tol)?;
    let curve = Curve::Line(line);
    let domain = curve.domain();
    Ok(make_edge_between(model, curve, domain, from.0, to.0, tol)?.shape)
}

/// A face on `surface` bounded by `edges` in traversal order, with an exact
/// same-parameter pcurve attached to every edge.
///
/// [`ogeom_algo::make_face_with_pcurves`] with one wire: the blend keeps this
/// thin name because every wedge face is a single loop.
pub(crate) fn face_from_edges(
    model: &mut Model,
    surface: SurfaceGeometry,
    edges: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Shape> {
    Ok(ogeom_algo::make_face_with_pcurves(model, surface, &[edges.to_vec()], tol)?.shape)
}

/// Sew the wedge's faces, demand a closed shell, and apply it to the solid —
/// subtracted on a convex edge, fused on a concave one. Either way the
/// history reads the same truth: the edge the blend replaces is gone.
pub(crate) fn apply_wedge(
    model: &mut Model,
    solid: &Shape,
    edge: Option<&Shape>,
    faces: &[Shape],
    additive: bool,
    tol: Tolerances,
) -> OgeomResult<ogeom_algo::Built> {
    let sewn = ogeom_algo::sew(model, faces, tol)?;
    if sewn.shells.len() != 1 || !ogeom_algo::is_shell_closed(model, &sewn.shells[0])? {
        ogeom_bail!(Construction, "the blend wedge did not close");
    }
    let wedge = ogeom_algo::make_solid(model, std::slice::from_ref(&sewn.shells[0]))?;
    let mut result = if additive {
        ogeom_bool::fuse(model, solid, &wedge.shape, tol)?
    } else {
        ogeom_bool::cut(model, solid, &wedge.shape, tol)?
    };
    if let Some(edge) = edge {
        result.history.delete(edge);
    }
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
) -> OgeomResult<Shape> {
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
