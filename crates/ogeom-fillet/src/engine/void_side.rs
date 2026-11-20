//! Concave rounds taken as the void's convex rounds (ADR 0001, amendment
//! "rolling corners on arcs, and concave corners from the void's side"):
//! a box about the request less the solid is the void beside it, where
//! every asked edge is convex; what rounding them there takes off the
//! void is the fill the solid gains.

use ogeom_algo::{Built, History};
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

/// How far, in radii, the box stands out from every asked edge: past the
/// run of a band beyond its ends, so no round of the void reaches a face
/// the box made. The odd fraction keeps the box's planes off the solid's
/// own round-numbered faces.
const REACH: f64 = 8.371;

/// `edges` of `solid`, every one concave, rounded at their radii from the
/// void's side; `None` where an edge is not concave, the void's rounds do
/// not build, a fill reaches a face of the box, or the solid gained is not
/// one sound solid.
///
/// # Errors
///
/// As the model's queries; a boolean or round that fails is `None`.
pub(crate) fn rounds(
    model: &mut Model,
    solid: &Shape,
    edges: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    if edges.is_empty()
        || model.kind_of(solid)? != ShapeType::Solid
        || model.children_of(solid)?.len() != 1
    {
        return Ok(None);
    }
    let mut reach = 0.0_f64;
    let mut low: Option<Point> = None;
    let mut high: Option<Point> = None;
    for (edge, radius) in edges {
        if crate::setback::edge_is_convex(model, solid, edge, *radius, tol)? != Some(false) {
            return Ok(None);
        }
        reach = reach.max(*radius * REACH);
        let bounds = ogeom_algo::shape_bounds(model, edge, tol)?;
        let (Some(a), Some(b)) = (bounds.low(), bounds.high()) else {
            return Ok(None);
        };
        low = Some(low.map_or(a, |l| l.min(a)));
        high = Some(high.map_or(b, |h| h.max(b)));
    }
    let (Some(low), Some(high)) = (low, high) else {
        return Ok(None);
    };
    let out = Vector::new(reach, reach, reach);
    let (corner, size) = (low - out, (high - low) + out * 2.0);
    let frame = Frame::new(corner, Direction::Z, Direction::X, tol)?;
    let Ok(boxed) = ogeom_algo::make_box(model, frame, (size.x, size.y, size.z), tol) else {
        return Ok(None);
    };
    let Ok(void) = ogeom_bool::cut(model, &boxed.shape, solid, tol) else {
        return Ok(None);
    };
    // Each asked edge's pieces in the void, grouped by the void's solids.
    let pieces = explore_unique(model, &void.shape, ShapeType::Solid)?;
    let mut asked: Vec<Vec<(Shape, f64)>> = vec![Vec::new(); pieces.len()];
    let mut image_of: Vec<(Shape, Vec<Shape>)> = Vec::with_capacity(edges.len());
    for (edge, radius) in edges {
        let mut images = Vec::new();
        for image in void.history.trace(edge) {
            let Some(k) = (0..pieces.len()).find(|&k| {
                explore_unique(model, &pieces[k], ShapeType::Edge)
                    .is_ok_and(|held| held.iter().any(|e| e.is_same(image)))
            }) else {
                continue;
            };
            asked[k].push((image.clone(), *radius));
            images.push(image.clone());
        }
        if images.is_empty() {
            return Ok(None);
        }
        image_of.push((edge.clone(), images));
    }
    let inner_low = corner + Vector::new(1.0, 1.0, 1.0) * (tol.confusion() * 100.0);
    let inner_high = corner + size - Vector::new(1.0, 1.0, 1.0) * (tol.confusion() * 100.0);
    let mut gained = solid.clone();
    let mut history = History::new();
    // Each void round's band faces, followed to the fills.
    let mut bands: Vec<(Shape, Vec<Shape>)> = Vec::new();
    for (k, piece) in pieces.iter().enumerate() {
        if asked[k].is_empty() {
            continue;
        }
        let Some(rounded) = crate::fillet::rounds_by_engine(model, piece, &asked[k], tol)? else {
            return Ok(None);
        };
        let Ok(fill) = ogeom_bool::cut(model, piece, &rounded.shape, tol) else {
            return Ok(None);
        };
        let fills = explore_unique(model, &fill.shape, ShapeType::Solid)?;
        if fills.is_empty() {
            return Ok(None);
        }
        // A fill touching a face the box made would be the box's corner,
        // not the solid's.
        let bounds = ogeom_algo::shape_bounds(model, &fill.shape, tol)?;
        let (Some(a), Some(b)) = (bounds.low(), bounds.high()) else {
            return Ok(None);
        };
        let inside = |p: Point| {
            p.x > inner_low.x
                && p.y > inner_low.y
                && p.z > inner_low.z
                && p.x < inner_high.x
                && p.y < inner_high.y
                && p.z < inner_high.z
        };
        if !(inside(a) && inside(b)) {
            return Ok(None);
        }
        for (edge, images) in &image_of {
            let mut faces = Vec::new();
            for image in images {
                if !asked[k].iter().any(|(e, _)| e.is_same(image)) {
                    continue;
                }
                for face in rounded.history.generated(image) {
                    faces.extend(fill.history.trace(face).iter().cloned());
                }
            }
            if !faces.is_empty() {
                bands.push((edge.clone(), faces));
            }
        }
        let Ok(fused) = ogeom_bool::fuse(model, &gained, &fill.shape, tol) else {
            return Ok(None);
        };
        // Band faces of earlier fills follow this fuse too.
        for (_, faces) in &mut bands {
            *faces = faces
                .iter()
                .flat_map(|f| fused.history.trace(f).to_vec())
                .collect();
        }
        history = history.then(&fused.history);
        gained = fused.shape;
    }
    let solids = explore_unique(model, &gained, ShapeType::Solid)?;
    let [one] = solids.as_slice() else {
        return Ok(None);
    };
    if model.children_of(one)?.len() != 1 || !ogeom_algo::check(model, one, tol)?.is_valid() {
        return Ok(None);
    }
    let kept = explore_unique(model, one, ShapeType::Face)?;
    for (edge, faces) in bands {
        for face in faces {
            if kept.iter().any(|f| f.is_same(&face)) {
                history.generate(&edge, face);
            }
        }
    }
    Ok(Some(Built::new(one.clone(), history)))
}
