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
    /// Edges refused (farther from the surface than the cap), with the
    /// offset each was measured at.
    pub refused: Vec<(Shape, f64)>,
}

/// Give a face's pcurve-less edges the trims projection can honestly fit.
///
/// The reader heals boundary slop up to a millimetre and hands what it
/// refuses over in `untrimmed_faces`, face shape included. This is the
/// follow-up: the same projection fit, at the cap the caller chooses. Each
/// fitted edge's tolerance widens to the offset actually measured and to
/// how far the fitted trim, lifted, stands from the curve anywhere along
/// it, so the model says what it knows. An edge past the cap is reported, not touched.
/// What the face shares with other shapes is copied first
/// ([`Model::unshare`]), so their edges keep the trims they have.
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
    let plan = plan_trims(model, face, tol)?;
    let fits = fit_trims(std::slice::from_ref(&plan), cap, tol);
    let mut report = FixedTrims::default();
    for (plan, fits) in [plan].into_iter().zip(fits) {
        attach_trims(model, plan, fits, tol, &mut report)?;
    }
    Ok(report)
}

/// A face's edge uses as the trim fix finds them, in wire order: each
/// with the curve to project where it holds no pcurve on the face.
pub(crate) struct TrimPlan {
    surface_id: ogeom_topo::SurfaceId,
    /// The face's surface where the face stands, made only when some edge
    /// needs a fit.
    surface: Option<ogeom_geom::SurfaceGeometry>,
    uses: Vec<(Shape, Option<Projected>)>,
}

/// A curve where its edge stands, over the edge's range, to project.
type Projected = (ogeom_geom::Curve, (f64, f64));

/// What one edge's fit came to: the pcurve with the worst sample offset
/// and the offset the edge must state, or a refusal at the offset
/// measured.
pub(crate) enum TrimFit {
    Fitted {
        pcurve: ogeom_geom::PlanarCurve,
        worst_off: f64,
        off: f64,
    },
    Refused(f64),
}

/// Find what [`fix_face_pcurves`] would fit on `face`, without fitting:
/// the face is unshared, and each edge use without a pcurve on it is
/// given its curve where the face stands.
///
/// # Errors
///
/// As [`fix_face_pcurves`].
pub(crate) fn plan_trims(
    model: &mut Model,
    face: &Shape,
    tol: Tolerances,
) -> OgeomResult<TrimPlan> {
    if model.kind_of(face)? != ShapeType::Face {
        ogeom_bail!(Construction, "fix_face_pcurves fixes a face");
    }
    model.unshare(face)?;
    let surface_id = {
        let Some(node) = model.node(face) else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let NodeData::Face(data) = node.data() else {
            ogeom_bail!(Construction, "face node holds no face data");
        };
        if model.geometry().surface(data.surface).is_none() {
            ogeom_bail!(Dangling, "face refers to a surface not in this model");
        }
        data.surface
    };
    let mut uses = Vec::new();
    for wire in model.ordered_children_of(face)? {
        for edge in model.ordered_children_of(&wire)? {
            let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                continue;
            };
            if data.pcurve_for(surface_id, edge.location()).is_some() {
                uses.push((edge, None));
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
            let job = (geometry.clone().transformed(&placed, tol)?, *range);
            uses.push((edge, Some(job)));
        }
    }
    let surface = if uses.iter().any(|(_, job)| job.is_some()) {
        let placement = face.transform(model.datums())?;
        model
            .geometry()
            .surface(surface_id)
            .map(|stored| stored.transformed(&placement, tol))
            .transpose()?
    } else {
        None
    };
    Ok(TrimPlan {
        surface_id,
        surface,
        uses,
    })
}

/// Fit every planned edge of `plans` in parallel, each use's fit in its
/// plan's place: `None` for a use that needs none.
pub(crate) fn fit_trims(
    plans: &[TrimPlan],
    cap: f64,
    tol: Tolerances,
) -> Vec<Vec<Option<OgeomResult<TrimFit>>>> {
    let jobs: Vec<(usize, usize)> = plans
        .iter()
        .enumerate()
        .flat_map(|(p, plan)| {
            plan.uses
                .iter()
                .enumerate()
                .filter(|(_, (_, job))| job.is_some())
                .map(move |(u, _)| (p, u))
        })
        .collect();
    let fitted = ogeom_core::parallel::map_ordered(&jobs, |_, &(p, u)| {
        let plan = &plans[p];
        let (Some((curve, range)), Some(surface)) = (&plan.uses[u].1, &plan.surface) else {
            ogeom_bail!(Construction, "a planned trim has nothing to fit");
        };
        fit_trim(curve, *range, surface, cap, tol)
    });
    let mut out: Vec<Vec<Option<OgeomResult<TrimFit>>>> = plans
        .iter()
        .map(|plan| plan.uses.iter().map(|_| None).collect())
        .collect();
    for ((p, u), fit) in jobs.into_iter().zip(fitted) {
        out[p][u] = Some(fit);
    }
    out
}

/// One edge's projected trim on a surface, held to `cap`.
fn fit_trim(
    curve: &ogeom_geom::Curve,
    range: (f64, f64),
    surface: &ogeom_geom::SurfaceGeometry,
    cap: f64,
    tol: Tolerances,
) -> OgeomResult<TrimFit> {
    use ogeom_algo::pcurve_fit::CappedFit;
    match ogeom_algo::pcurve_fit::fit_projected_pcurve_within(curve, range, surface, cap, tol) {
        Ok(CappedFit::Fitted((pcurve, _, _, worst_off, _))) => {
            // The edge states where the fitted chart lies, lifted and
            // measured densely against the curve, as well as the offset
            // its samples sat at.
            let gap = ogeom_algo::pcurve_fit::lifted_gap(
                (curve, range),
                (&pcurve, range),
                surface,
                false,
                tol,
            )?;
            Ok(TrimFit::Fitted {
                pcurve,
                worst_off,
                off: gap.max(worst_off),
            })
        }
        Ok(CappedFit::TooFar(off)) => Ok(TrimFit::Refused(off)),
        // A projection that cannot converge measured no offset.
        Err(_) => Ok(TrimFit::Refused(f64::INFINITY)),
    }
}

/// Attach a face's fitted trims in wire order, counting into `report`. A
/// use that has gained its pcurve since it was planned (an edge the wire
/// walks twice, or that a face on the same surface fixed first) is left
/// alone.
///
/// # Errors
///
/// As [`fix_face_pcurves`].
pub(crate) fn attach_trims(
    model: &mut Model,
    plan: TrimPlan,
    fits: Vec<Option<OgeomResult<TrimFit>>>,
    tol: Tolerances,
    report: &mut FixedTrims,
) -> OgeomResult<()> {
    for ((edge, job), fit) in plan.uses.into_iter().zip(fits) {
        let attached = model
            .node(&edge)
            .and_then(|n| n.data().as_edge())
            .is_some_and(|data| data.pcurve_for(plan.surface_id, edge.location()).is_some());
        let (Some((_, range)), Some(fit), false) = (job, fit, attached) else {
            report.already += 1;
            continue;
        };
        match fit? {
            TrimFit::Fitted {
                pcurve,
                worst_off,
                off,
            } => {
                report.fitted += 1;
                report.worst = report.worst.max(worst_off);
                // Its vertices widen with it: what bounds the edge is never
                // tighter than the edge.
                if off > tol.confusion() {
                    model.widen(&edge, ogeom_core::Tolerance::new(off + tol.confusion())?)?;
                }
                ogeom_algo::attach_pcurve(
                    model,
                    &edge,
                    pcurve,
                    plan.surface_id,
                    ogeom_topo::Location::identity(),
                    range,
                )?;
            }
            TrimFit::Refused(off) => report.refused.push((edge, off)),
        }
    }
    Ok(())
}

/// What [`reanchor_boundaries`] did.
#[derive(Debug, Default)]
pub struct ReanchoredBoundaries {
    /// Edges whose space curves moved onto their face's surface.
    pub moved: usize,
    /// The worst edge-to-surface offset found before moving.
    pub worst_before: f64,
    /// The worst residual after: the moved curve's distance from the old
    /// curve's projection, measured between the samples as well as at them.
    pub worst_after: f64,
    /// Edges refused (farther out than the cap), with their offsets.
    pub refused: Vec<(Shape, f64)>,
}

/// Move boundary curves onto the surfaces they are supposed to bound.
///
/// The stronger fix behind [`fix_face_pcurves`]: where that fits a *chart*
/// through whatever offset the boundary carries, this moves the boundary
/// itself: each off-surface edge's curve is projected, refitted at its own
/// parameters (so every chart already speaking the old curve keeps its
/// same-parameter law), and replaced throughout the shape. The displacement
/// is not hidden: the edge's and its vertices' tolerances widen to cover
/// where the boundary *was*, because the neighbouring faces still stand on
/// the unmoved geometry and honesty about the gap is what keeps them sewn.
///
/// An edge shared by several faces moves once, onto the first face that
/// claims it in face order; the recorded tolerance covers the rest.
/// What `shape` shares with other shapes is copied first
/// ([`ogeom_algo::on_own_nodes`]), so their vertices keep their tolerances.
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
    ogeom_algo::on_own_nodes_with(
        model,
        shape,
        &[],
        |model, shape, _| reanchor_boundaries_own(model, shape, cap, tol),
        |(built, _)| &mut built.history,
    )
}

/// [`reanchor_boundaries`] on a shape every node below which is its own.
fn reanchor_boundaries_own(
    model: &mut Model,
    shape: &Shape,
    cap: f64,
    tol: Tolerances,
) -> OgeomResult<(ogeom_algo::Built, ReanchoredBoundaries)> {
    use ogeom_topo::{Filter, explore};
    let mut report = ReanchoredBoundaries::default();
    let mut reshape = crate::reshape::Reshape::new();
    let mut done: ogeom_core::FastSet<ogeom_topo::TShapeId> = ogeom_core::FastSet::default();

    /// The fewest samples an edge is measured at.
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
            // Sampled at every span of the curve, so a boundary with more
            // detail than a fixed count of samples is measured whole.
            use ogeom_geom::Curve3d as _;
            use ogeom_geom::Surface as _;
            let stations = ogeom_algo::traced::stations(&curve, range, SAMPLES - 1);
            // Each projection seeded from where the nearest one before it
            // landed: consecutive points of a curve are neighbours on the
            // surface.
            let mut landed: Vec<(f64, (f64, f64), f64)> = Vec::with_capacity(stations.len());
            let mut project = |t: f64| -> OgeomResult<(ogeom_math::Point, f64)> {
                let p = curve.point_at(t, tol)?;
                let at = landed.partition_point(|l| l.0 < t);
                let seed = at.checked_sub(1).or((!landed.is_empty()).then_some(0));
                let hit = match seed.map(|i| landed[i].1) {
                    Some(uv) => ogeom_algo::project_on_surface_from(&surface, p, uv, tol)
                        .or_else(|_| ogeom_algo::project_on_surface(&surface, p, 24, tol))?,
                    None => ogeom_algo::project_on_surface(&surface, p, 24, tol)?,
                };
                landed.insert(at, (t, hit.parameters, hit.distance));
                Ok((
                    surface.point_at(hit.parameters.0, hit.parameters.1, tol)?,
                    hit.distance,
                ))
            };
            let mut worst = 0.0_f64;
            let mut shadow: ogeom_core::FastMap<u64, ogeom_math::Point> =
                ogeom_core::FastMap::with_capacity_and_hasher(stations.len(), Default::default());
            for &t in &stations {
                let (q, off) = project(t)?;
                worst = worst.max(off);
                shadow.insert(t.to_bits(), q);
            }
            if worst <= tol.confusion() * 1e3 {
                continue; // Already on the surface, to the reader's own bar.
            }
            report.worst_before = report.worst_before.max(worst);
            if worst > cap {
                report.refused.push((edge.clone(), worst));
                continue;
            }

            // Fitted at the curve's own parameters and held to the target
            // between the samples too; the error is the worst measured
            // anywhere, and the displacement the edge states covers it.
            let target = (tol.confusion() * 1e3).max(worst * 1e-3);
            let fitted = ogeom_algo::traced::fit_traced(
                |t| match shadow.get(&t.to_bits()) {
                    Some(q) => Ok(*q),
                    None => {
                        let (q, off) = project(t)?;
                        worst = worst.max(off);
                        Ok(q)
                    }
                },
                &stations,
                3,
                target,
                tol,
            )?;
            report.worst_before = report.worst_before.max(worst);
            report.worst_after = report.worst_after.max(fitted.error);
            let moved_by = worst + fitted.error + tol.confusion();

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
                    data.tolerance = data.tolerance.widen_to(moved_by);
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
                data.tolerance = data.tolerance.widen_to(moved_by);
                // The charts riding the old curve stay: each pcurve speaks
                // its own surface, whose geometry did not move, and the fit
                // at the old parameters keeps the same-parameter law.
                for repr in &reprs {
                    if !matches!(repr, EdgeRepr::Curve3d { .. }) {
                        data.add(repr.clone());
                    }
                }
            }
            // The kept charts are measured against the moved curve, and
            // the edge states how far they stand from it.
            ogeom_algo::state_pcurve_gaps_of(model, core::slice::from_ref(&rebuilt), tol)?;
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
