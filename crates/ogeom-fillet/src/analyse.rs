//! What a blend actually achieved, measured rather than asserted.
//!
//! A blend claims two things: that it meets each of its supports along a
//! curve, and that it meets them *smoothly* — the two surfaces sharing a
//! normal there. Both are claims about geometry that construction can get
//! subtly wrong: a fitted section drifts, a marched spine carries its
//! chord budget, a rebuilt face lands on a neighbour a hair off. This
//! module reports the two numbers instead of trusting them.
//!
//! The report is per shared edge, because that is where the claim lives.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve3d as _, Surface as _};
use ogeom_topo::{EdgeRepr, Model, NodeData, Shape, ShapeType, explore_unique};

/// How well a blend face meets one neighbour along one edge.
#[derive(Debug, Clone, PartialEq)]
pub struct BlendContact {
    /// The shared edge.
    pub edge: Shape,
    /// The face on the other side of it.
    pub neighbour: Shape,
    /// The largest angle, in radians, between the two surfaces' normals at
    /// the sampled stations — zero for a tangent join.
    pub tangency_error: f64,
    /// The largest distance, in model units, between the edge's curve and
    /// the two surfaces it is supposed to lie on.
    pub gap: f64,
    /// How many stations were sampled along the edge.
    pub stations: usize,
}

/// Measure a blend face against every face it shares an edge with.
///
/// Both numbers are sampled: `stations` points spread over each shared
/// edge's range, the surfaces read at the parameters the edge's own pcurves
/// give, so nothing is inverted and nothing is guessed. A blend that is
/// tangent everywhere except between two stations reports zero, which is
/// the honest limit of sampling and the reason the count is in the report.
///
/// An edge whose pcurve on either face is missing cannot be measured this
/// way at all: it is reported with its `gap` and `tangency_error` set to
/// infinity rather than quietly skipped, because a blend nobody can measure
/// is not a blend anybody should trust.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `blend` is not a face of `shape`, or `stations` is less than two.
pub fn analyse_blend(
    model: &Model,
    shape: &Shape,
    blend: &Shape,
    stations: usize,
    tol: Tolerances,
) -> OgeomResult<Vec<BlendContact>> {
    if stations < 2 {
        ogeom_bail!(
            Construction,
            "a blend measured at {stations} stations is not measured"
        );
    }
    let faces = explore_unique(model, shape, ShapeType::Face)?;
    if !faces.iter().any(|f| f.node() == blend.node()) {
        ogeom_bail!(Construction, "that face is not part of this shape");
    }
    let surface_of =
        |face: &Shape| -> Option<(ogeom_topo::SurfaceId, ogeom_geom::SurfaceGeometry)> {
            let NodeData::Face(data) = model.node(face)?.data() else {
                return None;
            };
            let geometry = model.geometry().surface(data.surface)?.clone();
            Some((data.surface, geometry))
        };
    let Some((blend_id, blend_surface)) = surface_of(blend) else {
        ogeom_bail!(Construction, "the blend face carries no surface");
    };
    let blend_edges = explore_unique(model, blend, ShapeType::Edge)?;

    let mut out = Vec::new();
    for neighbour in &faces {
        if neighbour.node() == blend.node() {
            continue;
        }
        let Some((other_id, other_surface)) = surface_of(neighbour) else {
            continue;
        };
        for edge in explore_unique(model, neighbour, ShapeType::Edge)? {
            if !blend_edges.iter().any(|e| e.node() == edge.node()) {
                continue;
            }
            let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                continue;
            };
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                continue;
            };
            let Some(geometry) = model.geometry().curve(*curve).cloned() else {
                continue;
            };
            let chart =
                |id: ogeom_topo::SurfaceId| -> Option<(ogeom_geom::PlanarCurve, (f64, f64))> {
                    match data.pcurve_for(id, edge.location())? {
                        EdgeRepr::PCurve { curve, range, .. } => {
                            Some((model.geometry().pcurve(*curve)?.clone(), *range))
                        }
                        EdgeRepr::Seam { forward, range, .. } => {
                            Some((model.geometry().pcurve(*forward)?.clone(), *range))
                        }
                        _ => None,
                    }
                };
            let (Some((pc_blend, pr_blend)), Some((pc_other, pr_other))) =
                (chart(blend_id), chart(other_id))
            else {
                out.push(BlendContact {
                    edge: edge.clone(),
                    neighbour: neighbour.clone(),
                    tangency_error: f64::INFINITY,
                    gap: f64::INFINITY,
                    stations: 0,
                });
                continue;
            };
            let (mut worst_angle, mut worst_gap) = (0.0f64, 0.0f64);
            for k in 0..stations {
                #[expect(
                    clippy::cast_precision_loss,
                    reason = "a station index, far below the mantissa"
                )]
                let f = k as f64 / (stations - 1) as f64;
                let t = (range.1 - range.0).mul_add(f, range.0);
                let on_curve = geometry.point_at(t, tol)?;
                let mut normals = Vec::with_capacity(2);
                for ((pcurve, prange), surface) in [
                    ((&pc_blend, pr_blend), &blend_surface),
                    ((&pc_other, pr_other), &other_surface),
                ] {
                    let pt = (prange.1 - prange.0).mul_add(f, prange.0);
                    let uv = ogeom_geom::Curve2d::point_at(pcurve, pt, tol)?;
                    worst_gap =
                        worst_gap.max(surface.point_at(uv.x, uv.y, tol)?.distance(on_curve));
                    normals.push(surface.normal_at(uv.x, uv.y, tol)?.vector());
                }
                // Orientation is the topology's business, not the join's:
                // two faces meeting smoothly may still be wound opposite
                // ways, so the angle is taken to the nearer sense.
                let dot = normals[0].dot(normals[1]).abs().min(1.0);
                worst_angle = worst_angle.max(dot.acos());
            }
            out.push(BlendContact {
                edge: edge.clone(),
                neighbour: neighbour.clone(),
                tangency_error: worst_angle,
                gap: worst_gap,
                stations,
            });
        }
    }
    Ok(out)
}
