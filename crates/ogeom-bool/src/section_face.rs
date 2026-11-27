//! Sectioning faces by a plane, with no solid behind them.
//!
//! The boolean's section needs two solids, because it runs the whole
//! arrangement and asks each piece which side of the other it lies on. A
//! face on its own (a sheet read from a file, one face picked from a part)
//! has no inside to ask about, and what a section of it means is simpler:
//! the curves where its surface crosses the plane, kept where they run on
//! the face.

use ogeom_algo::{
    Built, Containment, History, classify_on_face, make_edge_between, make_vertex,
    project_on_curve, shape_bounds,
};
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve3d as _, PlaneSurface, SurfaceGeometry, Transformable as _};
use ogeom_intersect::{CurveCurveOptions, IntersectOptions, SurfaceIntersection};
use ogeom_math::{Plane, Point};
use ogeom_mesh::Deflection;
use ogeom_topo::{EdgeRepr, Filter, Model, Shape, ShapeType, explore, explore_unique};

/// The edges where the faces of `shape` cross `plane`, each trimmed to the
/// face it lies on.
///
/// `shape` is a face, a shell, a solid or a compound of them. A curve is
/// exact where the face's surface and the plane have a closed-form section
/// (a line, a circle, an ellipse) and fitted otherwise. A face lying in the
/// plane gives its own boundary. A curve running along a face's boundary (a
/// seam, or an edge the plane holds) is kept once, however many faces or
/// sides of a seam it bounds. A plane that misses every face gives an
/// empty compound.
///
/// # Errors
///
/// As the surface and curve intersectors, and
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `shape` has no faces.
pub fn section_face(
    model: &mut Model,
    shape: &Shape,
    plane: &Plane,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let faces = explore(model, shape, Filter::OfType(ShapeType::Face))?;
    if faces.is_empty() {
        ogeom_core::ogeom_bail!(Construction, "a section by a plane needs faces to cut");
    }
    // The plane, windowed to reach past the shape on every side, so its
    // sections with the faces are bounded curves.
    let bounds = shape_bounds(model, shape, tol)?;
    let (Some(low), Some(high)) = (bounds.low(), bounds.high()) else {
        return empty(model, shape);
    };
    let reach = low.distance(high).max(tol.confusion() * 1e3);
    let centre = plane.frame().to_local(low.midpoint(high));
    let cutter: SurfaceGeometry = PlaneSurface::over(
        *plane,
        (centre.x - reach, centre.x + reach),
        (centre.y - reach, centre.y + reach),
    )?
    .into();

    let mut edges: Vec<Shape> = Vec::new();
    let mut vertices: Vec<(Point, Shape)> = Vec::new();
    let mut seen_boundary: Vec<ogeom_topo::TShapeId> = Vec::new();
    // What has been emitted along a face's boundary, by a stretch on it or
    // as a coplanar face's edge: the face across that edge, or the other
    // side of a seam, finds the same stretch, and it is kept once.
    let mut along: Vec<Curve> = Vec::new();
    for face in &faces {
        // A face whose box stands wholly on one side of the plane cannot
        // meet it: its corners' signed distances say so without asking the
        // intersector.
        let corners = shape_bounds(model, face, tol)?.corners();
        let side = |c: &Point| plane.signed_distance_to(*c);
        if !corners.is_empty()
            && (corners.iter().all(|c| side(c) > tol.confusion())
                || corners.iter().all(|c| side(c) < -tol.confusion()))
        {
            continue;
        }
        let Some(surface_id) = model
            .node(face)
            .and_then(|n| n.data().as_face())
            .map(|d| d.surface)
        else {
            continue;
        };
        let Some(surface) = model.geometry().surface(surface_id) else {
            continue;
        };
        let placed = surface.transformed(&face.transform(model.datums())?, tol)?;
        let found = ogeom_intersect::intersect_surfaces(
            &placed,
            &cutter,
            IntersectOptions {
                tolerance: tol.confusion() * 0.5,
                ..IntersectOptions::for_tolerances(tol)
            },
            tol,
        )?;
        let curves = match found {
            SurfaceIntersection::Along(curves) => curves,
            // A face in the plane: its boundary is where the two meet.
            SurfaceIntersection::Same => {
                for (edge, curve) in boundary_curves(model, face, tol)? {
                    if seen_boundary.contains(&edge.node()) {
                        continue;
                    }
                    seen_boundary.push(edge.node());
                    let (lo, hi) = curve.domain();
                    let middle = curve.point_at(f64::midpoint(lo, hi), tol)?;
                    if !emitted(&along, middle, tol)? {
                        edges.push(edge);
                        along.push(curve);
                    }
                }
                continue;
            }
            SurfaceIntersection::Apart | SurfaceIntersection::Touching(_) => continue,
        };
        let boundary: Vec<Curve> = boundary_curves(model, face, tol)?
            .into_iter()
            .map(|(_, c)| c)
            .collect();
        for section in curves {
            for (lo, hi, on_boundary) in on_face(model, face, &section.curve, &boundary, tol)? {
                if on_boundary {
                    let middle = section.curve.point_at(f64::midpoint(lo, hi), tol)?;
                    if emitted(&along, middle, tol)? {
                        continue;
                    }
                    along.push(
                        ogeom_geom::TrimmedCurve::new(section.curve.clone(), lo, hi, tol)?.into(),
                    );
                }
                let from = section.curve.point_at(lo, tol)?;
                let to = section.curve.point_at(hi, tol)?;
                let v0 = vertex(model, &mut vertices, from, tol)?;
                let v1 = vertex(model, &mut vertices, to, tol)?;
                edges.push(
                    make_edge_between(model, section.curve.clone(), (lo, hi), &v0, &v1, tol)?.shape,
                );
            }
        }
    }
    let result = model.add_compound(&edges)?;
    let mut history = History::new();
    history.modify(shape, result.clone());
    Ok(Built::new(result, history))
}

fn empty(model: &mut Model, shape: &Shape) -> OgeomResult<Built> {
    let result = model.add_compound(&[])?;
    let mut history = History::new();
    history.modify(shape, result.clone());
    Ok(Built::new(result, history))
}

/// One vertex per place: neighbouring faces' sections meet at the edge they
/// share, and there their ends are one vertex.
pub(crate) fn vertex(
    model: &mut Model,
    vertices: &mut Vec<(Point, Shape)>,
    at: Point,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    if let Some((p, v)) = vertices
        .iter()
        .find(|(p, _)| p.distance(at) <= tol.confusion() * 10.0)
    {
        // The shared vertex stands for this end too: it must reach it, or
        // the edge ending here is refused as missing its vertex.
        let gap = p.distance(at) * 1.01;
        if gap > tol.confusion()
            && let Ok(reach) = ogeom_core::Tolerance::new(gap)
        {
            model.widen(v, reach)?;
        }
        return Ok(v.clone());
    }
    let v = make_vertex(model, at).shape;
    vertices.push((at, v.clone()));
    Ok(v)
}

/// Whether `p` lies on one of the stretches already emitted.
fn emitted(along: &[Curve], p: Point, tol: Tolerances) -> OgeomResult<bool> {
    for curve in along {
        if project_on_curve(curve, p, 32, tol)?.distance <= tol.confusion() * 10.0 {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The face's boundary edges, each with its curve over its own range,
/// placed.
fn boundary_curves(
    model: &Model,
    face: &Shape,
    tol: Tolerances,
) -> OgeomResult<Vec<(Shape, Curve)>> {
    let mut out = Vec::new();
    for edge in explore_unique(model, face, ShapeType::Edge)? {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            continue;
        };
        let placed = geometry
            .clone()
            .transformed(&edge.transform(model.datums())?, tol)?;
        if let Ok(trimmed) = ogeom_geom::TrimmedCurve::new(placed, range.0, range.1, tol) {
            out.push((edge, trimmed.into()));
        }
    }
    Ok(out)
}

/// The stretches of `curve` that run on the face: split wherever it
/// crosses or leaves the face's boundary, each stretch kept if its middle
/// lies on the face, and marked when that middle is on the face's boundary
/// (a seam, or an edge the plane holds) rather than inside it.
fn on_face(
    model: &Model,
    face: &Shape,
    curve: &Curve,
    boundary: &[Curve],
    tol: Tolerances,
) -> OgeomResult<Vec<(f64, f64, bool)>> {
    let (lo, hi) = curve.domain();
    let mut cuts = vec![lo, hi];
    let options = CurveCurveOptions {
        gap: tol.confusion() * 10.0,
        ..CurveCurveOptions::for_tolerances(tol)
    };
    for edge in boundary {
        let found = ogeom_intersect::intersect_curves(curve, edge, options, tol)?;
        cuts.extend(found.crossings.iter().map(|c| c.on_a));
        for overlap in &found.overlaps {
            cuts.push(overlap.on_a.0);
            cuts.push(overlap.on_a.1);
        }
    }
    cuts.retain(|t| *t >= lo && *t <= hi);
    cuts.sort_by(f64::total_cmp);
    cuts.dedup_by(|a, b| (*a - *b).abs() <= tol.parametric());
    let mut kept: Vec<(f64, f64, bool)> = Vec::new();
    for pair in cuts.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        if b - a <= tol.parametric() {
            continue;
        }
        let middle = curve.point_at(f64::midpoint(a, b), tol)?;
        let on_boundary =
            match classify_on_face(model, face, middle, Deflection::for_tolerances(tol), tol)? {
                Containment::In => false,
                Containment::On => true,
                Containment::Out => continue,
            };
        // Stretches that meet end to end on the face, both inside it or
        // both on its boundary, are one edge.
        match kept.last_mut() {
            Some((_, end, on)) if *on == on_boundary && (*end - a).abs() <= tol.parametric() => {
                *end = b;
            }
            _ => kept.push((a, b, on_boundary)),
        }
    }
    Ok(kept)
}
