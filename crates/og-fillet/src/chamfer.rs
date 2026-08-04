//! Chamfers: the bevel that replaces an edge.
//!
//! P4's opening stone, built deliberately on M3's shoulders: a chamfer along
//! a straight edge between two planar faces is a wedge subtracted, and the
//! wedge's own faces lie *exactly* on the solid's — coplanar, materials
//! aligned — which is the same-domain case the boolean learned to resolve.
//! The blend machinery proper (rolling-ball fillets, variable radii, vertex
//! blends) comes after; a chamfer is the member of the family that needs no
//! new surface, only the machinery already proven.

use og_algo::{Built, find_plane, make_face, make_polygon};
use og_bool::cut;
use og_core::{OgResult, Tolerances, og_bail};
use og_geom::Curve3d as _;
use og_geom::{Curve, SurfaceGeometry};
use og_math::{Point, Vector};
use og_topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, explore};

/// Bevel a straight edge of a solid, cutting `distance` back along each of
/// its two faces.
///
/// The edge must be straight and shared by exactly two planar faces; the
/// distances are equal (the symmetric chamfer). The result is the boolean
/// difference with a wedge whose legs run along the two faces — so the
/// history reads as a cut: the two faces are modified into their trimmed
/// pieces, the edge's neighbourhood gains the bevel face.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the edge is
/// not straight, is not shared by exactly two planar faces of `solid`, or
/// `distance` is not a usable length.
pub fn chamfer_edge(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    distance: f64,
    tol: Tolerances,
) -> OgResult<Built> {
    if !distance.is_finite() || distance <= tol.confusion() {
        og_bail!(Construction, "a chamfer of {distance} cuts nothing");
    }

    // The edge's line, in world space.
    let (curve, range) = {
        let Some(node) = model.node(edge) else {
            og_bail!(Dangling, "edge is not in this model");
        };
        let NodeData::Edge(data) = node.data() else {
            og_bail!(Construction, "expected an edge");
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            og_bail!(Construction, "the edge has no curve to bevel along");
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            og_bail!(Dangling, "curve is not in this model");
        };
        (geometry.clone(), *range)
    };
    let Curve::Line(_) = &curve else {
        og_bail!(
            Construction,
            "chamfering a curved edge needs the blend machinery; this is the \
             straight-edge chamfer"
        );
    };
    let start = curve.point_at(range.0, tol)?;
    let end = curve.point_at(range.1, tol)?;
    let along = (end - start) / start.distance(end);

    // The two faces meeting at the edge, with their outward normals.
    let mut adjacent: Vec<Vector> = Vec::new();
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
                "chamfering an edge of a curved face needs the blend \
                 machinery; this is the planar chamfer"
            );
        };
        let placement = face.transform(model.datums())?;
        let mut normal = placement.apply_vector(plane.plane().normal().vector());
        if face.orientation() == og_topo::Orientation::Reversed {
            normal = -normal;
        }
        adjacent.push(normal);
    }
    if adjacent.len() != 2 {
        og_bail!(
            Construction,
            "a chamfer needs an edge shared by exactly two faces, found {}",
            adjacent.len()
        );
    }

    // On each face, the direction perpendicular to the edge that walks away
    // from the *other* face — into the material each leg cuts back along.
    let leg = |own: Vector, other: Vector| -> OgResult<Vector> {
        let mut t = own.cross(along);
        if t.dot(other) > 0.0 {
            t = -t;
        }
        let m = t.magnitude();
        if m <= tol.confusion() {
            og_bail!(Construction, "a face is tangent to its own edge");
        }
        Ok(t / m)
    };
    let a = leg(adjacent[0], adjacent[1])?;
    let b = leg(adjacent[1], adjacent[0])?;

    // The wedge: a triangular prism whose apex line is the edge and whose
    // legs run `distance` along each face. Built from five explicit planar
    // faces rather than swept, because a sweep's walls are extrusion
    // surfaces even when they are geometrically planes, and the boolean's
    // same-domain resolution — which is what makes the coplanar legs melt
    // into the solid's own faces — recognises coincidence between *planes*.
    let travel = end - start;
    let apex0 = start;
    let a0 = start + a * distance;
    let b0 = start + b * distance;
    let quad = |p0: Point, p1: Point| [p0, p1, p1 + travel, p0 + travel];
    let centroid = Point::new(
        (apex0.x + a0.x + b0.x) / 3.0 + travel.x * 0.5,
        (apex0.y + a0.y + b0.y) / 3.0 + travel.y * 0.5,
        (apex0.z + a0.z + b0.z) / 3.0 + travel.z * 0.5,
    );
    let mut faces = Vec::new();
    let rings: [Vec<Point>; 5] = [
        vec![apex0, a0, b0],
        vec![apex0 + travel, a0 + travel, b0 + travel],
        quad(apex0, a0).to_vec(),
        quad(a0, b0).to_vec(),
        quad(b0, apex0).to_vec(),
    ];
    for ring in rings {
        faces.push(planar_face(model, &ring, centroid, tol)?);
    }
    let sewn = og_algo::sew(model, &faces, tol)?;
    if sewn.shells.len() != 1 || !og_algo::is_shell_closed(model, &sewn.shells[0])? {
        og_bail!(Construction, "the chamfer wedge did not close");
    }
    let wedge = og_algo::make_solid(model, std::slice::from_ref(&sewn.shells[0]))?;

    let mut result = cut(model, solid, &wedge.shape, tol)?;
    result.history.delete(edge);
    Ok(result)
}

/// A planar face over explicit corners, oriented outward from `inside`, with
/// same-parameter pcurves on its own plane.
fn planar_face(
    model: &mut Model,
    corners: &[Point],
    inside: Point,
    tol: Tolerances,
) -> OgResult<Shape> {
    let wire = make_polygon(model, corners, true, tol)?.shape;
    let Some(mut plane) = find_plane(model, &wire, tol)? else {
        og_bail!(Construction, "a wedge face is degenerate");
    };
    if plane.normal().vector().dot(corners[0] - inside) < 0.0 {
        plane = plane.reversed();
    }
    let mut reach = 1.0_f64;
    for p in corners {
        reach = reach.max(p.distance(corners[0]) * 2.0);
    }
    let surface = og_geom::PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?;
    let face = make_face(model, surface.into(), std::slice::from_ref(&wire), tol)?.shape;
    let (surface_id, frame) = {
        let Some(node) = model.node(&face) else {
            og_bail!(Dangling, "the face just built is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            og_bail!(Construction, "the face holds no face data");
        };
        (data.surface, plane.frame())
    };
    for pedge in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
        let (pcurve, prange) = {
            let Some(node) = model.node(&pedge) else {
                og_bail!(Dangling, "edge is not in this model");
            };
            let Some(data) = node.data().as_edge() else {
                og_bail!(Construction, "edge node holds no edge data");
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                continue;
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                og_bail!(Dangling, "curve is not in this model");
            };
            let flat = |p: Point| {
                let l = frame.to_local(p);
                og_math::Point2::new(l.x, l.y)
            };
            let from = flat(geometry.point_at(range.0, tol)?);
            let to = flat(geometry.point_at(range.1, tol)?);
            (
                og_geom::PlanarCurve::from(og_geom::Line2d::segment(from, to, tol)?),
                *range,
            )
        };
        og_algo::attach_pcurve(
            model,
            &pedge,
            pcurve,
            surface_id,
            og_topo::Location::identity(),
            prange,
        )?;
    }
    Ok(face)
}
