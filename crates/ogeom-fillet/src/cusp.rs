//! A cusp corner: where a flat side runs tangent into a cylindrical wall,
//! the flat top between them narrows to a cusp at the tangent point, and a
//! request rounds the top's edge with the side (convex) and its edge with
//! the wall (concave) together.
//!
//! Within `w = 2r` of the cusp the two bands no longer fit side by side on
//! the top, and the sound round there is an S: off the side in a convex
//! turn, onto the wall in a concave one, flattening into the walls' common
//! tangent plane at the tangent point. The bands are built, a cell is cut
//! round the tangent point out to where the cusp is `4r` wide, and the cell
//! is closed again: the side and the wall each down to (or up from) a trim
//! along the S's own tangent line on it, and on each side of the tangent
//! point one patch, G1 to both bands' ends, the top, the side and the wall,
//! held inside the solid and the wall's concave fill.
//!
//! On the side the trim is the S's own foot, section by section along the
//! side. On the wall it leaves the tangent point up the wall's line, along
//! which the wall's normal does not turn, and reaches the concave band's
//! rail along it: at the tangent point the side's normal holds still along
//! the side's trim, and a patch tangent to both there can only be one whose
//! normal also holds still along the wall's.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve3d as _, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use crate::setback::point_of;

/// A face's surface where it stands.
fn placed(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Option<SurfaceGeometry>> {
    use ogeom_geom::Transformable as _;
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    Ok(Some(
        surface
            .clone()
            .transformed(&face.transform(model.datums())?, tol)?,
    ))
}

/// An edge's unit tangent at `at`, an end of it, pointing away from it.
fn leaving(model: &Model, edge: &Shape, at: Point, tol: Tolerances) -> OgeomResult<Option<Vector>> {
    let Some((curve, (lo, hi))) = ogeom_algo::placed_edge_curve(model, edge, tol)? else {
        return Ok(None);
    };
    let (a, b) = (curve.point_at(lo, tol)?, curve.point_at(hi, tol)?);
    let (t, sense) = if a.distance(at) <= b.distance(at) {
        (lo, 1.0)
    } else {
        (hi, -1.0)
    };
    Ok(Some(curve.d1_at(t, tol)?.normalized(tol)? * sense))
}

/// The cusp corner at `vertex`, where `at` (each request edge there with
/// its radius and kind) meet, or `None` where it is not one: the vertex's
/// faces a planar top, a planar side square to it and a cylindrical wall
/// about an axis square to the top, the side tangent to the wall at the
/// vertex; each convex edge a line between the top and the side, each
/// concave one an arc between the top and the wall, every one leaving the
/// vertex along with one of the other kind; one radius for all.
pub(crate) fn cusp_corner(
    model: &Model,
    solid: &Shape,
    vertex: &Shape,
    at: &[(Shape, f64, bool)],
    tol: Tolerances,
) -> OgeomResult<bool> {
    let Some(blend) = at.first().map(|e| e.1) else {
        return Ok(false);
    };
    if at.iter().any(|e| (e.1 - blend).abs() > tol.confusion()) {
        return Ok(false);
    }
    let corner = point_of(model, vertex)?;
    let mut planes: Vec<(Shape, Vector)> = Vec::new();
    let mut walls: Vec<(Shape, ogeom_math::Cylinder)> = Vec::new();
    for f in explore_unique(model, solid, ShapeType::Face)? {
        if !explore_unique(model, &f, ShapeType::Vertex)?
            .iter()
            .any(|v| v.is_same(vertex))
        {
            continue;
        }
        match placed(model, &f, tol)? {
            Some(SurfaceGeometry::Plane(_)) => {
                let (_, outward) = ogeom_algo::face_normal(model, &f, tol)?;
                planes.push((f, outward));
            }
            Some(SurfaceGeometry::Cylinder(c)) => walls.push((f, c.cylinder())),
            _ => return Ok(false),
        }
    }
    let ([(top_a, na), (top_b, nb)], [(wall, cylinder)]) = (planes.as_slice(), walls.as_slice())
    else {
        return Ok(false);
    };
    // The top is the plane a concave edge of the request runs on.
    let on = |face: &Shape, edge: &Shape| -> OgeomResult<bool> {
        Ok(explore_unique(model, face, ShapeType::Edge)?
            .iter()
            .any(|e| e.is_same(edge)))
    };
    let Some(concave) = at.iter().find(|e| !e.2) else {
        return Ok(false);
    };
    let (up, out) = if on(top_a, &concave.0)? {
        (*na, *nb)
    } else if on(top_b, &concave.0)? {
        (*nb, *na)
    } else {
        return Ok(false);
    };
    let axis = cylinder.axis();
    let radius = cylinder.radius();
    let square = tol.angular() * 10.0;
    if up.dot(out).abs() > square || axis.direction.vector().cross(up).magnitude() > square {
        return Ok(false);
    }
    // The side tangent to the wall at the vertex: the axis `radius` behind
    // it, square to it, through the vertex.
    let behind = corner - out * radius;
    let off = (behind - axis.location) - up * (behind - axis.location).dot(up);
    if off.magnitude() > tol.confusion() * 100.0 {
        return Ok(false);
    }
    let _ = wall;
    // Each convex edge a line on the side, leaving along with a concave arc.
    let mut sides: Vec<Vector> = Vec::new();
    let mut paired = vec![false; at.len()];
    for (i, (edge, _, convex)) in at.iter().enumerate() {
        if !convex {
            continue;
        }
        let Some((curve, _)) = ogeom_algo::placed_edge_curve(model, edge, tol)? else {
            return Ok(false);
        };
        if !matches!(curve, Curve::Line(_)) {
            return Ok(false);
        }
        let Some(t) = leaving(model, edge, corner, tol)? else {
            return Ok(false);
        };
        let mut found = false;
        for (j, (other, _, other_convex)) in at.iter().enumerate() {
            if *other_convex || paired[j] {
                continue;
            }
            if leaving(model, other, corner, tol)?.is_some_and(|u| u.dot(t) > 1.0 - square) {
                paired[j] = true;
                found = true;
                break;
            }
        }
        if !found {
            return Ok(false);
        }
        paired[i] = true;
        sides.push(t);
    }
    if sides.is_empty() || paired.iter().any(|p| !p) {
        return Ok(false);
    }
    Ok(true)
}
