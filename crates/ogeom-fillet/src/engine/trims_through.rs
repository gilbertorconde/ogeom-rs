//! Trim through off (ADR 0001, amendment "overflow handling", "trim
//! through holds per blend"): reading, from a built result, whether a
//! blend of the request reached a face it is not seated on, or carried a
//! face on past its side over nothing, which only trim through builds.
//!
//! The check holds per blend. Each point where the result's boundary
//! newly crosses a face of the solid is attributed to the blend whose
//! reach holds it most tightly: the request edge nearest it, its distance
//! read in that blend's own size. The face must be the attributed blend's
//! own, one of its edge's two faces or a face at its ends, so a blend
//! trimming a face is refused even where another blend of the request is
//! seated on that face.

use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve3d as _, SurfaceGeometry};
use ogeom_math::{Aabb, Point};
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of, explore_unique};

/// What a result built with trim through holds that trim through off
/// refuses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Through {
    /// A blend changed a face it is not seated on.
    Reaches,
    /// The result lays a face of the solid on past its side, over
    /// nothing: a face carried flat as a ledge past where the solid turns
    /// away from it.
    Carries,
}

/// The interior points read along each edge of the result.
const ALONG: u32 = 8;

/// Whether `built`, the blends of `edges` on `solid` (each edge with its
/// blend's size: a round's largest radius, a bevel's longer leg), holds
/// what only trim through builds ([`Through`]); `None` where it does not.
///
/// A face of the solid other than the edges' own and those at their ends
/// that the result changed is [`Through::Reaches`], as is a point where
/// the result's boundary newly crosses a face whose attributed blend is
/// not seated on it. An edge of the result lying along a face's surface,
/// at a point no face of the solid holds, is [`Through::Carries`], unless
/// the face is one at the ends of the attributed blend's edge (not one of
/// its two faces), which a concave band's end extends to cap its fill.
///
/// # Errors
///
/// As exploring the shapes and measuring the faces.
pub(crate) fn through(
    model: &Model,
    solid: &Shape,
    edges: &[(Shape, f64)],
    built: &Built,
    tol: Tolerances,
) -> OgeomResult<Option<Through>> {
    let shapes: Vec<Shape> = edges.iter().map(|(e, _)| e.clone()).collect();
    if changes_others(model, solid, &shapes, built, tol)? {
        return Ok(Some(Through::Reaches));
    }
    crossings(model, solid, edges, built, tol)
}

/// A face of the solid, its surface and its bounds.
struct Seat {
    face: Shape,
    surface: SurfaceGeometry,
    bounds: Aabb,
}

/// An edge of the request: its curve, its blend's size, its two faces
/// and the other faces at its ends.
struct Asked {
    curve: Curve,
    range: (f64, f64),
    size: f64,
    sides: Vec<Shape>,
    ends: Vec<Shape>,
}

/// An edge of the solid with its curve and bounds.
struct Side {
    curve: Curve,
    range: (f64, f64),
    bounds: Aabb,
}

/// The per-blend reading of [`through`]: every edge of the result lying
/// along a face's surface, read at interior points away from the solid's
/// own edges.
fn crossings(
    model: &Model,
    solid: &Shape,
    edges: &[(Shape, f64)],
    built: &Built,
    tol: Tolerances,
) -> OgeomResult<Option<Through>> {
    let eps = tol.confusion() * 10.0;
    let deflection = ogeom_mesh::Deflection::for_tolerances(tol);
    let mut asked = Vec::new();
    for (edge, size) in edges {
        let Some((curve, range)) = ogeom_algo::placed_edge_curve(model, edge, tol)? else {
            return Ok(None);
        };
        let sides = ancestors_of(model, solid, edge, ShapeType::Face)?;
        let mut ends = Vec::new();
        for vertex in explore_unique(model, edge, ShapeType::Vertex)? {
            for face in ancestors_of(model, solid, &vertex, ShapeType::Face)? {
                if !sides.iter().any(|f| f.is_same(&face)) {
                    ends.push(face);
                }
            }
        }
        asked.push(Asked {
            curve,
            range,
            size: size.max(tol.confusion()),
            sides,
            ends,
        });
    }
    let mut seats = Vec::new();
    for face in explore_unique(model, solid, ShapeType::Face)? {
        let Some(surface) = model
            .node(&face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .cloned()
        else {
            continue;
        };
        let bounds = ogeom_algo::shape_bounds(model, &face, tol)?.expanded(eps);
        seats.push(Seat {
            face,
            surface,
            bounds,
        });
    }
    let mut sides = Vec::new();
    for edge in explore_unique(model, solid, ShapeType::Edge)? {
        if let Some((curve, range)) = ogeom_algo::placed_edge_curve(model, &edge, tol)? {
            let bounds = ogeom_algo::shape_bounds(model, &edge, tol)?.expanded(eps);
            sides.push(Side {
                curve,
                range,
                bounds,
            });
        }
    }
    for edge in explore_unique(model, &built.shape, ShapeType::Edge)? {
        let Some((curve, (a, b))) = ogeom_algo::placed_edge_curve(model, &edge, tol)? else {
            continue;
        };
        let mut points = Vec::new();
        for k in 1..ALONG {
            let t = a + (b - a) * f64::from(k) / f64::from(ALONG);
            points.push(curve.point_at(t, tol)?);
        }
        // The faces whose surface the edge runs along, not merely crosses.
        let mut along = Vec::new();
        for seat in &seats {
            let mut on = true;
            for &p in &points {
                let off = ogeom_geom::project::project_on_surface(&seat.surface, p, 16, tol)?;
                if off.distance > eps {
                    on = false;
                    break;
                }
            }
            if on {
                along.push(seat);
            }
        }
        if along.is_empty() {
            continue;
        }
        for &p in &points {
            if on_side(&sides, p, eps, tol)? {
                continue;
            }
            let mut holding = Vec::new();
            for seat in &along {
                if seat.bounds.contains(p)
                    && ogeom_algo::classify_on_face(model, &seat.face, p, deflection, tol)?
                        != ogeom_algo::Containment::Out
                {
                    holding.push(&seat.face);
                }
            }
            if holding.is_empty() {
                // A band's end closes against a face at its ends, which
                // is extended over a concave band's fill: only a face the
                // band is seated on, or none of its own, is carried.
                let faces: Vec<&Shape> = along.iter().map(|s| &s.face).collect();
                if !tightest(&asked, p, eps, tol)?
                    .iter()
                    .any(|e| e.ends.iter().any(|o| faces.iter().any(|f| f.is_same(o))))
                {
                    return Ok(Some(Through::Carries));
                }
                continue;
            }
            if !tightest(&asked, p, eps, tol)?.iter().any(|e| {
                e.sides
                    .iter()
                    .chain(&e.ends)
                    .any(|o| holding.iter().any(|h| h.is_same(o)))
            }) {
                return Ok(Some(Through::Reaches));
            }
        }
    }
    Ok(None)
}

/// Whether `p` lies on an edge of the solid.
fn on_side(sides: &[Side], p: Point, eps: f64, tol: Tolerances) -> OgeomResult<bool> {
    for side in sides {
        if side.bounds.contains(p)
            && crate::support::distance_to_edge(&side.curve, side.range, p, tol)? <= eps
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The blends whose reach holds `p` most tightly: the distance to the
/// edge over the blend's size the least of the request's, to within `eps`.
fn tightest(asked: &[Asked], p: Point, eps: f64, tol: Tolerances) -> OgeomResult<Vec<&Asked>> {
    let mut reads = Vec::with_capacity(asked.len());
    for edge in asked {
        let d = crate::support::distance_to_edge(&edge.curve, edge.range, p, tol)?;
        reads.push((d, edge));
    }
    let least = reads
        .iter()
        .map(|(d, e)| d / e.size)
        .fold(f64::INFINITY, f64::min);
    Ok(reads
        .into_iter()
        .filter(|(d, e)| (d - eps) / e.size <= least)
        .map(|(_, e)| e)
        .collect())
}

/// Whether `built`, the blends of `edges` on `solid`, changed a face of
/// `solid` other than the edges' own and those at their ends: a face the
/// blends reached and trimmed through. A face is unchanged where the
/// result holds a face of its area and its bounds on whose surface a point
/// of it lies.
///
/// # Errors
///
/// As exploring the shapes and measuring the faces.
fn changes_others(
    model: &Model,
    solid: &Shape,
    edges: &[Shape],
    built: &Built,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let eps = tol.confusion() * 10.0;
    let deflection = ogeom_mesh::Deflection::for_tolerances(tol);
    let mut own: Vec<Shape> = Vec::new();
    for edge in edges {
        own.extend(ancestors_of(model, solid, edge, ShapeType::Face)?);
        for vertex in explore_unique(model, edge, ShapeType::Vertex)? {
            own.extend(ancestors_of(model, solid, &vertex, ShapeType::Face)?);
        }
    }
    let read = |face: &Shape| -> OgeomResult<(f64, Aabb)> {
        Ok((
            ogeom_algo::surface_properties(model, face, deflection, tol)?.mass,
            ogeom_algo::shape_bounds(model, face, tol)?,
        ))
    };
    let same_box = |a: &Aabb, b: &Aabb| match (a.low(), a.high(), b.low(), b.high()) {
        (Some(al), Some(ah), Some(bl), Some(bh)) => al.distance(bl).max(ah.distance(bh)) <= eps,
        _ => false,
    };
    let made: Vec<(Shape, f64, Aabb)> = explore_unique(model, &built.shape, ShapeType::Face)?
        .into_iter()
        .map(|f| read(&f).map(|(area, bounds)| (f, area, bounds)))
        .collect::<OgeomResult<_>>()?;
    for face in explore_unique(model, solid, ShapeType::Face)? {
        if own.iter().any(|o| o.is_same(&face)) {
            continue;
        }
        let (area, bounds) = read(&face)?;
        let (p, _) = ogeom_algo::face_normal(model, &face, tol)?;
        let mut kept = false;
        for (g, now, around) in &made {
            if (area - now).abs() > 1e-9 * (1.0 + area) || !same_box(&bounds, around) {
                continue;
            }
            let Some(surface) = model
                .node(g)
                .and_then(|x| x.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
            else {
                continue;
            };
            let on = ogeom_geom::project::project_on_surface(surface, p, 16, tol)?.distance;
            if on <= eps {
                kept = true;
                break;
            }
        }
        if !kept {
            return Ok(true);
        }
    }
    Ok(false)
}
