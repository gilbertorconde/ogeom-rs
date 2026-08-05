//! Verifying — and where needed, widening into truth — the `same_parameter`
//! claim.
//!
//! An edge carries several representations of one curve, and nearly every
//! algorithm evaluates whichever is convenient, assuming the answers are
//! interchangeable within the edge's tolerance. The flag that records this is
//! set false whenever a representation is added — honest, but pessimistic:
//! every primitive's edges claim a disagreement they do not have. This is the
//! repair the flag's own documentation demands: measure the actual
//! disagreement, and either confirm the claim or widen the edge's tolerance
//! until the claim is true. Either way, afterwards the flag *means* something.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve2d as _, Curve3d as _, Surface as _};
use ogeom_topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, TShapeId, explore};

/// What one repair pass did.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SameParameterReport {
    /// Edges examined.
    pub checked: usize,
    /// Edges whose representations already agreed within tolerance.
    pub agreed: usize,
    /// Edges whose tolerance had to widen to make the claim true.
    pub widened: usize,
    /// Edges with no pcurves to disagree with — trivially true.
    pub trivial: usize,
}

/// Verify every edge under `shape` and make its `same_parameter` flag true.
///
/// Each pcurve is sampled against the edge's own curve at matched parameters
/// — the same linear range mapping the triangulator uses — and the worst gap
/// decides: within the edge's tolerance, the claim is confirmed; beyond it,
/// the tolerance widens to cover what was measured, which makes the claim
/// true by making the tolerance honest. Degenerate edges and edges with no
/// pcurves are trivially true.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the shape
/// or a representation resolves to nothing.
pub fn repair_same_parameter(
    model: &mut Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgeomResult<SameParameterReport> {
    const SAMPLES: usize = 24;
    let mut report = SameParameterReport::default();
    let mut done: Vec<TShapeId> = Vec::new();
    for edge in explore(model, shape, Filter::OfType(ShapeType::Edge))? {
        if done.contains(&edge.node()) {
            continue;
        }
        done.push(edge.node());
        report.checked += 1;

        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let data = data.clone();
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            // A degenerate edge's pcurve is its whole story; there is nothing
            // for it to disagree with.
            report.trivial += 1;
            set_flag(model, &edge, true);
            continue;
        };
        let Some(curve) = model.geometry().curve(*curve).cloned() else {
            report.trivial += 1;
            continue;
        };
        let (ca, cb) = *range;

        // Every pcurve representation, seam sides included.
        let mut pairs: Vec<(
            ogeom_geom::PlanarCurve,
            (f64, f64),
            ogeom_geom::SurfaceGeometry,
        )> = Vec::new();
        for repr in &data.representations {
            match repr {
                EdgeRepr::PCurve {
                    curve: pc,
                    surface,
                    range,
                    ..
                } => {
                    if let (Some(p), Some(s)) = (
                        model.geometry().pcurve(*pc).cloned(),
                        model.geometry().surface(*surface).cloned(),
                    ) {
                        pairs.push((p, *range, s));
                    }
                }
                EdgeRepr::Seam {
                    forward,
                    reversed,
                    surface,
                    range,
                    ..
                } => {
                    for pc in [forward, reversed] {
                        if let (Some(p), Some(s)) = (
                            model.geometry().pcurve(*pc).cloned(),
                            model.geometry().surface(*surface).cloned(),
                        ) {
                            pairs.push((p, *range, s));
                        }
                    }
                }
                EdgeRepr::Curve3d { .. } => {}
                _ => {}
            }
        }
        if pairs.is_empty() {
            report.trivial += 1;
            set_flag(model, &edge, true);
            continue;
        }

        let mut worst = 0.0_f64;
        for (pcurve, prange, surface) in &pairs {
            for i in 0..=SAMPLES {
                #[allow(clippy::cast_precision_loss)]
                let t = ca + (cb - ca) * i as f64 / SAMPLES as f64;
                // The linear range mapping every consumer uses.
                let u = if (cb - ca).abs() <= f64::MIN_POSITIVE {
                    prange.0
                } else {
                    prange.0 + (prange.1 - prange.0) * (t - ca) / (cb - ca)
                };
                let on_curve = curve.point_at(t, tol)?;
                let uv = pcurve.point_at(u, tol)?;
                let lifted = surface.point_at(uv.x, uv.y, tol)?;
                worst = worst.max(on_curve.distance(lifted));
            }
        }

        let within = {
            let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                continue;
            };
            data.tolerance.get()
        };
        if worst <= within {
            report.agreed += 1;
        } else {
            report.widened += 1;
            if let Some(node) = model.node_mut(&edge)
                && let NodeData::Edge(data) = node.data_mut()
            {
                data.tolerance = data.tolerance.widen_to(worst + tol.confusion());
            }
        }
        set_flag(model, &edge, true);
    }
    Ok(report)
}

/// Set an edge's `same_parameter` claim.
fn set_flag(model: &mut Model, edge: &Shape, agrees: bool) {
    if let Some(node) = model.node_mut(edge)
        && let NodeData::Edge(data) = node.data_mut()
    {
        data.assert_same_parameter(agrees);
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use ogeom_math::Frame;

    const T: Tolerances = Tolerances::millimetres();

    fn all_flags_true(model: &Model, shape: &Shape) -> bool {
        explore(model, shape, Filter::OfType(ShapeType::Edge))
            .unwrap()
            .iter()
            .all(|e| {
                model
                    .node(e)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(ogeom_topo::EdgeData::same_parameter)
            })
    }

    #[test]
    fn a_primitives_edges_agree_and_the_flag_finally_says_so() {
        let mut model = Model::new();
        let solid = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T).unwrap();
        assert!(
            !all_flags_true(&model, &solid.shape),
            "the builder is honest: unverified means false"
        );
        let report = repair_same_parameter(&mut model, &solid.shape, T).unwrap();
        assert_eq!(report.widened, 0, "a native primitive has nothing to widen");
        assert!(report.agreed > 0);
        assert!(all_flags_true(&model, &solid.shape));
    }

    #[test]
    fn a_disagreeing_pcurve_widens_the_tolerance_into_truth() {
        use ogeom_geom::{Line2d, LineCurve, PlaneSurface};
        use ogeom_math::{Plane, Point, Point2};
        let mut model = Model::new();
        let curve: ogeom_geom::Curve =
            LineCurve::segment(Point::new(0.0, 0.0, 0.0), Point::new(10.0, 0.0, 0.0), T)
                .unwrap()
                .into();
        let edge = ogeom_algo::make_edge(&mut model, curve, (0.0, 10.0), T)
            .unwrap()
            .shape;
        let surface = model.geometry_mut().add_surface(
            PlaneSurface::over(Plane::new(Frame::WORLD), (-20.0, 20.0), (-20.0, 20.0))
                .unwrap()
                .into(),
        );
        // A pcurve half a unit off the curve it claims to follow.
        let off = Line2d::segment(Point2::new(0.0, 0.5), Point2::new(10.0, 0.5), T).unwrap();
        ogeom_algo::attach_pcurve(
            &mut model,
            &edge,
            off.into(),
            surface,
            ogeom_topo::Location::identity(),
            (0.0, 10.0),
        )
        .unwrap();

        let report = repair_same_parameter(&mut model, &edge, T).unwrap();
        assert_eq!(report.widened, 1);
        let data = model.node(&edge).unwrap().data().as_edge().unwrap();
        assert!(data.same_parameter());
        assert!(
            data.tolerance.get() >= 0.5,
            "the tolerance covers the measured gap, got {}",
            data.tolerance.get()
        );
    }
}
