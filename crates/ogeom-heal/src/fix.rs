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
/// The reader heals boundary slop up to a millimetre and names what it
/// refuses in `untrimmed_faces`; this is the instructed follow-up — the
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
