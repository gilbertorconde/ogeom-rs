//! Continuation of a band stopped short by another band at the same vertex.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_math::{Direction, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType};

#[allow(
    clippy::too_many_arguments,
    reason = "the corner and the band it closes"
)]
pub(super) fn carry(
    model: &mut Model,
    solid: &Shape,
    corner: Point,
    far: Point,
    axes: &[Direction],
    k: usize,
    planes: [Vector; 2],
    radius: f64,
    from_vertex: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let along = axes[k].vector();
    let on_line = tol.confusion() * 100.0;
    let point_of = |model: &Model, v: &Shape| -> OgeomResult<Option<Point>> {
        let Some(p) = model
            .node(v)
            .and_then(|n| n.data().as_vertex().map(|d| d.point))
        else {
            return Ok(None);
        };
        Ok(Some(v.transform(model.datums())?.apply(p)))
    };
    // Where the band starts along the edge, from the block's plane: the
    // nearest of its vertices.
    let mut ends: Vec<Point> = Vec::new();
    for face in ogeom_topo::explore_unique(model, solid, ShapeType::Face)? {
        if !on_band(model, &face, far, along, radius, tol) {
            continue;
        }
        for v in ogeom_topo::explore_unique(model, &face, ShapeType::Vertex)? {
            if let Some(p) = point_of(model, &v)? {
                ends.push(p);
            }
        }
    }
    let Some(start) = ends
        .iter()
        .map(|p| (*p - far).dot(along))
        .min_by(f64::total_cmp)
    else {
        return Ok(None);
    };
    if start <= tol.confusion() {
        return Ok(None);
    }
    // The band's near end stands against another edge's band of this
    // corner: where its section there meets the crease, a point on that
    // band's cylinder.
    let from = (corner - far).dot(along);
    let gap = |p: Point, axis: Vector| {
        let w = p - far;
        ((w - axis * w.dot(axis)).magnitude() - radius).abs()
    };
    let crease_end = corner + along * (start - from);
    if !axes
        .iter()
        .enumerate()
        .any(|(j, axis)| j != k && gap(crease_end, axis.vector()) <= on_line)
    {
        return Ok(None);
    }
    // No stretch of the crease between the vertex and the band is still a
    // sharp edge of the solid. The crease leaves the other band at its rail,
    // where the band is tangent to the host plane the crease lies in, so
    // the two stand apart only quadratically there and the boolean places
    // the crossing within its vertices' tolerances. A piece each of whose
    // ends lies on one other band's cylinder within the end's own
    // tolerance lies on it within those tolerances all along (a line's
    // distance from a cylinder it stays outside is convex along it): it is
    // that band's boundary, not a stretch of sharp crease.
    for edge in ogeom_topo::explore_unique(model, solid, ShapeType::Edge)? {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, &edge)? else {
            continue;
        };
        let (Some(pa), Some(pb)) = (point_of(model, &a)?, point_of(model, &b)?) else {
            continue;
        };
        let crease = |p: Point| {
            let w = p - corner;
            (w - along * w.dot(along)).magnitude() <= on_line
        };
        if !crease(pa) || !crease(pb) {
            continue;
        }
        let reach = |v: &Shape| {
            model
                .node(v)
                .and_then(|n| n.data().as_vertex())
                .map_or(0.0, |d| d.tolerance.get())
        };
        let (reach_a, reach_b) = (reach(&a), reach(&b));
        let on_band_seam = axes.iter().enumerate().any(|(j, axis)| {
            j != k && gap(pa, axis.vector()) <= reach_a && gap(pb, axis.vector()) <= reach_b
        });
        if on_band_seam {
            continue;
        }
        let (ta, tb) = ((pa - far).dot(along), (pb - far).dot(along));
        if ta.max(tb) > from + on_line && ta.min(tb) < start - on_line {
            return Ok(None);
        }
    }
    let face_on = |model: &Model, inward: Vector| -> OgeomResult<Option<Shape>> {
        for face in ogeom_topo::explore_unique(model, solid, ShapeType::Face)? {
            let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
                continue;
            };
            if !matches!(
                model.geometry().surface(data.surface),
                Some(ogeom_geom::SurfaceGeometry::Plane(_))
            ) {
                continue;
            }
            let (origin, outward) = ogeom_algo::face_normal(model, &face, tol)?;
            if (corner - origin).dot(outward).abs() <= on_line
                && outward.cross(inward).magnitude() < tol.angular() * 10.0
                && outward.dot(inward) < 0.0
            {
                return Ok(Some(face));
            }
        }
        Ok(None)
    };
    let (Some(first), Some(second)) = (face_on(model, planes[0])?, face_on(model, planes[1])?)
    else {
        return Ok(None);
    };
    let foot = corner + along * (far - corner).dot(along);
    let seat = crate::support::Seat {
        start: if from_vertex { corner } else { foot },
        end: foot + along * start,
        along,
        normals: [-planes[0], -planes[1]],
        faces: [first, second],
        convex: true,
    };
    let cut = super::vertex_tools::cut_strip(model, solid, &seat, radius, tol)?;
    // The carried band and the band it reaches lie on one cylinder: one
    // face, as the band the corner closes in one call is.
    let on_axis: Vec<Shape> = ogeom_topo::explore_unique(model, &cut.shape, ShapeType::Face)?
        .into_iter()
        .filter(|face| on_band(model, face, far, along, radius, tol))
        .collect();
    let (unified, _) = ogeom_heal::unify_same_domain_around(model, &cut.shape, &on_axis, tol)?;
    Ok(Some(Built {
        shape: unified.shape,
        history: cut.history.then(&unified.history),
    }))
}

/// Whether `face` lies on the cylinder of `radius` about the line through
/// `far` along `along`.
fn on_band(
    model: &Model,
    face: &Shape,
    far: Point,
    along: Vector,
    radius: f64,
    tol: Tolerances,
) -> bool {
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
        return false;
    };
    let Some(ogeom_geom::SurfaceGeometry::Cylinder(c)) = model.geometry().surface(data.surface)
    else {
        return false;
    };
    let frame = c.cylinder().frame();
    let axis = frame.z().vector();
    let off = far - frame.origin();
    (c.cylinder().radius() - radius).abs() <= tol.confusion()
        && axis.cross(along).magnitude() <= tol.angular()
        && (off - axis * off.dot(axis)).magnitude() <= tol.confusion()
}
