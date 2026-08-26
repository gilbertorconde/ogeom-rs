//! Fixing faces: recomputing the trims a face is missing.
//!
//! *Elsewhere:* the `ShapeFix_Face` / `ShapeFix_Edge` corner of the fixing
//! family.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Transformable as _;
use ogeom_topo::{EdgeRepr, Model, NodeData, Shape, ShapeType};

/// What [`fix_face_pcurves`] did, edge by edge.
#[derive(Debug, Default)]
pub struct FixedTrims {
    /// Edges that gained a fitted pcurve.
    pub fitted: usize,
    /// Edges that already carried one and were left alone.
    pub already: usize,
    /// The worst measured edge-to-surface offset among the fitted, now
    /// recorded in those edges' widened tolerances.
    pub worst: f64,
    /// Edges refused — farther from the surface than the cap — with the
    /// offset each was measured at.
    pub refused: Vec<(Shape, f64)>,
}

/// Give a face's pcurve-less edges the trims projection can honestly fit.
///
/// The reader heals boundary slop up to a millimetre and hands what it
/// refuses over in `untrimmed_faces`, face shape included; this is the
/// instructed follow-up — the
/// same projection fit, at the cap the caller chooses. Each fitted edge's
/// tolerance widens to the offset actually measured, so the model says
/// what it now knows; an edge past the cap is reported, not touched.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `face` is not a face, holds no surface, or an edge carries no space
/// curve to project.
pub fn fix_face_pcurves(
    model: &mut Model,
    face: &Shape,
    cap: f64,
    tol: Tolerances,
) -> OgeomResult<FixedTrims> {
    if model.kind_of(face)? != ShapeType::Face {
        ogeom_bail!(Construction, "fix_face_pcurves fixes a face");
    }
    let (surface_id, surface) = {
        let Some(node) = model.node(face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        let Some(stored) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
        };
        let placement = face.transform(model.datums())?;
        (data.surface, stored.transformed(&placement, tol)?)
    };

    let mut report = FixedTrims::default();
    for wire in model.ordered_children_of(face)? {
        for edge in model.ordered_children_of(&wire)? {
            let (curve, range) = {
                let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                    continue;
                };
                if data.pcurve_for(surface_id, edge.location()).is_some() {
                    report.already += 1;
                    continue;
                }
                let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                    ogeom_bail!(
                        Construction,
                        "an edge has no space curve; nothing can be projected"
                    );
                };
                let Some(geometry) = model.geometry().curve(*curve) else {
                    ogeom_bail!(Dangling, "an edge names a curve not in this model");
                };
                let placed = edge.transform(model.datums())?;
                (geometry.clone().transformed(&placed, tol)?, *range)
            };
            match ogeom_algo::pcurve_fit::fit_projected_pcurve_capped(
                &curve, range, &surface, cap, tol,
            ) {
                Ok((pcurve, _, _, worst_off, _)) => {
                    report.fitted += 1;
                    report.worst = report.worst.max(worst_off);
                    if worst_off > tol.confusion()
                        && let Some(node) = model.node_mut(&edge)
                        && let NodeData::Edge(data) = node.data_mut()
                    {
                        data.tolerance = data.tolerance.widen_to(worst_off + tol.confusion());
                    }
                    ogeom_algo::attach_pcurve(
                        model,
                        &edge,
                        pcurve,
                        surface_id,
                        ogeom_topo::Location::identity(),
                        range,
                    )?;
                }
                Err(refusal) => {
                    // The measured offset travels in the message; the report
                    // carries the number a consumer acts on.
                    let off = refusal
                        .to_string()
                        .split_whitespace()
                        .find_map(|w| w.parse::<f64>().ok())
                        .unwrap_or(f64::INFINITY);
                    report.refused.push((edge.clone(), off));
                }
            }
        }
    }
    Ok(report)
}

/// What [`reanchor_boundaries`] did.
#[derive(Debug, Default)]
pub struct ReanchoredBoundaries {
    /// Edges whose space curves moved onto their face's surface.
    pub moved: usize,
    /// The worst edge-to-surface offset found before moving.
    pub worst_before: f64,
    /// The worst residual after — the fit's honest distance from the
    /// projected samples.
    pub worst_after: f64,
    /// Edges refused — farther out than the cap — with their offsets.
    pub refused: Vec<(Shape, f64)>,
}

/// Move boundary curves onto the surfaces they are supposed to bound.
///
/// The stronger fix behind [`fix_face_pcurves`]: where that fits a *chart*
/// through whatever offset the boundary carries, this moves the boundary
/// itself — each off-surface edge's curve is projected, refitted at its own
/// parameters (so every chart already speaking the old curve keeps its
/// same-parameter law), and replaced throughout the shape. The displacement
/// is not hidden: the edge's and its vertices' tolerances widen to cover
/// where the boundary *was*, because the neighbouring faces still stand on
/// the unmoved geometry and honesty about the gap is what keeps them sewn.
///
/// An edge shared by several faces moves once, onto the first face that
/// claims it in face order; the recorded tolerance covers the rest.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// shape's structure resists rebuilding; refusals past the cap are reported,
/// not thrown.
pub fn reanchor_boundaries(
    model: &mut Model,
    shape: &Shape,
    cap: f64,
    tol: Tolerances,
) -> OgeomResult<(ogeom_algo::Built, ReanchoredBoundaries)> {
    use ogeom_topo::{Filter, explore};
    let mut report = ReanchoredBoundaries::default();
    let mut reshape = crate::reshape::Reshape::new();
    let mut done: std::collections::HashSet<ogeom_topo::TShapeId> =
        std::collections::HashSet::new();

    const SAMPLES: usize = 33;
    for face in explore(model, shape, Filter::OfType(ShapeType::Face))? {
        let surface = {
            let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
                continue;
            };
            let Some(stored) = model.geometry().surface(data.surface) else {
                continue;
            };
            let placement = face.transform(model.datums())?;
            stored.clone().transformed(&placement, tol)?
        };
        for edge in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
            if !done.insert(edge.node()) {
                continue;
            }
            let (curve, range, reprs) = {
                let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                    continue;
                };
                let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                    continue;
                };
                let Some(geometry) = model.geometry().curve(*curve) else {
                    continue;
                };
                let placed = edge.transform(model.datums())?;
                (
                    geometry.clone().transformed(&placed, tol)?,
                    *range,
                    data.representations.clone(),
                )
            };
            // Measure, then move only what is honestly off and under the cap.
            use ogeom_geom::Curve3d as _;
            use ogeom_geom::Surface as _;
            let mut params = Vec::with_capacity(SAMPLES);
            let mut projected = Vec::with_capacity(SAMPLES);
            let mut worst = 0.0_f64;
            let mut seed: Option<(f64, f64)> = None;
            for i in 0..SAMPLES {
                #[allow(clippy::cast_precision_loss, reason = "a sample index")]
                let t = range.0 + (range.1 - range.0) * i as f64 / (SAMPLES - 1) as f64;
                let p = curve.point_at(t, tol)?;
                let hit = match seed {
                    Some(uv) => ogeom_algo::project_on_surface_from(&surface, p, uv, tol)
                        .or_else(|_| ogeom_algo::project_on_surface(&surface, p, 24, tol))?,
                    None => ogeom_algo::project_on_surface(&surface, p, 24, tol)?,
                };
                seed = Some(hit.parameters);
                worst = worst.max(hit.distance);
                params.push(t);
                projected.push(surface.point_at(hit.parameters.0, hit.parameters.1, tol)?);
            }
            if worst <= tol.confusion() * 1e3 {
                continue; // Already on the surface, to the reader's own bar.
            }
            report.worst_before = report.worst_before.max(worst);
            if worst > cap {
                report.refused.push((edge.clone(), worst));
                continue;
            }

            let fitted = ogeom_geom::fit::fit_points_at(
                &params,
                &projected,
                3,
                (tol.confusion() * 1e3).max(worst * 1e-3),
                tol,
            )?;
            report.worst_after = report.worst_after.max(fitted.error);

            // The move is recorded before it is made: ends and edge widen to
            // cover where the boundary was, so every neighbour still meets
            // it within stated tolerance.
            // Stored order, not traversal order: the curve's range runs the
            // stored way, and the rebuilt edge's ends must match it however
            // this occurrence happens to be oriented.
            let bounds = model.children_of(&edge)?;
            let (Some(va), Some(vb)) = (bounds.first().cloned(), bounds.last().cloned()) else {
                continue;
            };
            for v in [&va, &vb] {
                if let Some(node) = model.node_mut(v)
                    && let NodeData::Vertex(data) = node.data_mut()
                {
                    data.tolerance = data.tolerance.widen_to(worst + tol.confusion());
                }
            }
            let rebuilt = ogeom_algo::make_edge_between(
                model,
                ogeom_geom::Curve::BSpline(fitted.curve),
                (range.0, range.1),
                &va,
                &vb,
                tol,
            )?
            .shape;
            if let Some(node) = model.node_mut(&rebuilt)
                && let NodeData::Edge(data) = node.data_mut()
            {
                data.tolerance = data.tolerance.widen_to(worst + tol.confusion());
                // The charts riding the old curve stay: each pcurve speaks
                // its own surface, whose geometry did not move, and the fit
                // at the old parameters keeps the same-parameter law.
                for repr in &reprs {
                    if !matches!(repr, EdgeRepr::Curve3d { .. }) {
                        data.add(repr.clone());
                    }
                }
            }
            reshape.replace(&edge, rebuilt);
            report.moved += 1;
        }
    }
    if reshape.is_empty() {
        return Ok((ogeom_algo::Built::from_nothing(shape.clone()), report));
    }
    let built = reshape.apply(model, shape)?;
    Ok((built, report))
}
