//! Verifying (and where needed, widening into truth) the `same_parameter`
//! claim.
//!
//! An edge carries several representations of one curve, and nearly every
//! algorithm evaluates whichever is convenient, assuming the answers are
//! interchangeable within the edge's tolerance. The flag that records this is
//! set false whenever a representation is added: honest, but pessimistic:
//! every primitive's edges claim a disagreement they do not have. This is the
//! repair the flag's own documentation demands: measure the actual
//! disagreement, and either confirm the claim or widen the edge's tolerance
//! until the claim is true. Either way, afterwards the flag *means* something.

use ogeom_core::{FastSet, OgeomResult, Tolerance, Tolerances};
use ogeom_topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, explore};

/// What one repair pass did.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SameParameterReport {
    /// Edges examined.
    pub checked: usize,
    /// Edges whose representations already agreed within tolerance.
    pub agreed: usize,
    /// Edges whose tolerance had to widen to make the claim true.
    pub widened: usize,
    /// Edges with no pcurves to disagree with; trivially true.
    pub trivial: usize,
}

/// Verify every edge under `shape` and make its `same_parameter` flag true.
///
/// Each pcurve, on its surface where the representation places it, is
/// lifted against the edge's own curve at matched parameters (the same
/// linear range mapping the triangulator uses), and the worst gap decides:
/// within the edge's tolerance, the claim is confirmed; beyond it, the
/// tolerance widens to cover what was measured, which makes the claim true
/// by making the tolerance honest. Degenerate edges and edges with no
/// pcurves are trivially true. An edge whose pcurves cannot be measured
/// keeps its flag as it was. The edges are measured in parallel and
/// widened in order.
///
/// Nodes below `shape` that other shapes hold as well are copied first
/// ([`Model::unshare`]), so the repair reaches no other shape.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the shape
/// resolves to nothing.
pub fn repair_same_parameter(
    model: &mut Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgeomResult<SameParameterReport> {
    model.unshare(shape)?;
    let mut seen = FastSet::default();
    let edges: Vec<Shape> = explore(model, shape, Filter::OfType(ShapeType::Edge))?
        .into_iter()
        .filter(|edge| seen.insert(edge.node()))
        .collect();
    let found = {
        let model = &*model;
        ogeom_core::parallel::map_ordered(&edges, |_, edge| measure(model, edge, tol))
    };
    let mut report = SameParameterReport::default();
    for (edge, found) in edges.iter().zip(found) {
        report.checked += 1;
        match found {
            Found::Absent | Found::Unmeasured => {}
            Found::Trivial { flag } => {
                report.trivial += 1;
                if flag {
                    set_flag(model, edge, true);
                }
            }
            Found::Gap(worst) => {
                let within = model
                    .node(edge)
                    .and_then(|n| n.data().as_edge())
                    .map_or(0.0, |data| data.tolerance.get());
                if worst <= within {
                    report.agreed += 1;
                } else {
                    report.widened += 1;
                    model.widen(edge, Tolerance::new(worst + tol.confusion())?)?;
                }
                set_flag(model, edge, true);
            }
        }
    }
    Ok(report)
}

/// What measuring one edge found.
enum Found {
    /// Not an edge in this model.
    Absent,
    /// Nothing to disagree: no curve in space, or no pcurve. `flag` where
    /// the claim is true for it.
    Trivial { flag: bool },
    /// Pcurves that could not be measured against the curve.
    Unmeasured,
    /// The widest gap between a lifted pcurve and the curve.
    Gap(f64),
}

fn measure(model: &Model, edge: &Shape, tol: Tolerances) -> Found {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        return Found::Absent;
    };
    let Some(EdgeRepr::Curve3d { curve, .. }) = data.curve3d() else {
        // A degenerate edge's pcurve is its whole story. There is nothing
        // for it to disagree with.
        return Found::Trivial { flag: true };
    };
    if model.geometry().curve(*curve).is_none() {
        return Found::Trivial { flag: false };
    }
    let any_pcurve = data.representations.iter().any(|repr| match repr {
        EdgeRepr::PCurve { curve, surface, .. } => {
            model.geometry().pcurve(*curve).is_some()
                && model.geometry().surface(*surface).is_some()
        }
        EdgeRepr::Seam {
            forward,
            reversed,
            surface,
            ..
        } => {
            model.geometry().surface(*surface).is_some()
                && [forward, reversed]
                    .iter()
                    .any(|pc| model.geometry().pcurve(**pc).is_some())
        }
        _ => false,
    });
    if !any_pcurve {
        return Found::Trivial { flag: true };
    }
    match crate::upgrade::measure_edge(model, edge, true, tol) {
        Some(worst) if worst.is_finite() => Found::Gap(worst),
        _ => Found::Unmeasured,
    }
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
