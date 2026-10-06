//! How far an edge's pcurves stand off its curve, and edges' tolerances
//! raised to cover it.
//!
//! An edge's tolerance is the radius about its curve within which every
//! description of the edge lies, each pcurve lifted through its surface
//! included. A producer that fits or projects a pcurve knows its error
//! only where it sampled; these measure it densely, after the fact.

use ogeom_core::{OgeomResult, Tolerance, Tolerances};
use ogeom_geom::{Curve, Curve2d as _, Curve3d as _, PlanarCurve, Surface as _, SurfaceGeometry};
use ogeom_math::Point;
use ogeom_topo::{EdgeRepr, Model, Shape, ShapeType};

/// How much wider than the gap it was measured from an edge's tolerance is
/// recorded: a millionth.
const MARGIN: f64 = 1e-6;

/// How far a pcurve, lifted through its surface, stands off an edge's
/// curve: at each of many samples along both ranges alike, the distance
/// between the two points, or where that exceeds `stated`, the distance
/// from the lifted point to the nearest point of the curve's stretch, and
/// the widest peaks between samples searched for their tops, so a stretch
/// of the edge sampled anywhere stays within the answer. Within `stated`,
/// the answer is only known to be within it.
///
/// Unlike [`lifted_gap`](crate::pcurve_fit::lifted_gap), which charges a
/// pcurve for a pace differing from the curve's by more than a sample's
/// step, this measures where the pcurve lies, as `check` does: an iso line
/// run the other way round a closed rim stands on the rim.
///
/// # Errors
///
/// If the curve or the pcurve cannot be evaluated over its range.
pub fn pcurve_gap(
    (curve, range): (&Curve, (f64, f64)),
    (pcurve, pcurve_range): (&PlanarCurve, (f64, f64)),
    surface: &SurfaceGeometry,
    stated: f64,
    tol: Tolerances,
) -> OgeomResult<f64> {
    // The checker's sample counts divide this one, so every point it
    // samples is one of these.
    const SAMPLES: u32 = 256;
    // A sample this close to the widest may stand beside a peak wider
    // than every sample: the stretch either side of it is searched, for
    // at most this many of the widest.
    const NEAR_PEAK: f64 = 0.9;
    const PEAKS: usize = 4;
    // The gap at `f` of the way along, the nearest point of the curve
    // standing in where the same parameter's is past `within`.
    let gap_beyond = |f: f64, within: f64| -> OgeomResult<Option<f64>> {
        let uv = pcurve.point_at(pcurve_range.0 + (pcurve_range.1 - pcurve_range.0) * f, tol)?;
        let Ok(lifted) = surface.point_at(uv.x, uv.y, tol) else {
            return Ok(None);
        };
        let gap = curve
            .point_at(range.0 + (range.1 - range.0) * f, tol)?
            .distance(lifted);
        Ok(Some(if gap > within {
            gap.min(nearest_on_stretch(curve, range, lifted, tol)?)
        } else {
            gap
        }))
    };
    let step = 1.0 / f64::from(SAMPLES);
    let mut widest: f64 = 0.0;
    let mut gaps = Vec::with_capacity(SAMPLES as usize + 1);
    for i in 0..=SAMPLES {
        // A gap no wider than the widest so far changes nothing, and is
        // kept as found.
        let gap = gap_beyond(f64::from(i) * step, widest.max(stated))?;
        if let Some(gap) = gap {
            widest = widest.max(gap);
        }
        gaps.push(gap);
    }
    // A piece of the edge is sampled between these points, so the peaks
    // between them count: the widest local peaks near the widest sample are
    // searched for their tops.
    let sampled = widest;
    let gap_of = |j: Option<usize>| j.and_then(|j| gaps.get(j).copied().flatten());
    let mut peaks: Vec<(f64, usize)> = gaps
        .iter()
        .enumerate()
        .filter_map(|(i, gap)| {
            let gap = (*gap)?;
            let peak = gap >= sampled * NEAR_PEAK
                && gap_of(i.checked_sub(1)).is_none_or(|g| g <= gap)
                && gap_of(Some(i + 1)).is_none_or(|g| g <= gap);
            peak.then_some((gap, i))
        })
        .collect();
    peaks.sort_by(|a, b| b.0.total_cmp(&a.0));
    let gap_at = |f: f64| gap_beyond(f, stated);
    for &(_, i) in peaks.iter().take(PEAKS) {
        #[allow(clippy::cast_precision_loss, reason = "a sample index")]
        let at = i as f64 * step;
        widest = widest.max(peak_about(&gap_at, at, step)?);
    }
    Ok(widest)
}

/// The widest of `gap_at` within a step either side of `at`, by a
/// golden-section search, never less than at `at`.
fn peak_about(
    gap_at: &impl Fn(f64) -> OgeomResult<Option<f64>>,
    at: f64,
    step: f64,
) -> OgeomResult<f64> {
    const ROUNDS: u32 = 40;
    let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
    let mut widest = gap_at(at)?.unwrap_or(0.0);
    let (mut a, mut b) = ((at - step).max(0.0), (at + step).min(1.0));
    for _ in 0..ROUNDS {
        let (c, d) = (b - (b - a) * ratio, a + (b - a) * ratio);
        let (gc, gd) = (gap_at(c)?.unwrap_or(0.0), gap_at(d)?.unwrap_or(0.0));
        widest = widest.max(gc).max(gd);
        if gc > gd {
            b = d;
        } else {
            a = c;
        }
    }
    Ok(widest)
}

/// The distance from `target` to the nearest point of `curve` over `range`:
/// a scan, then a golden-section search about the nearest sample.
fn nearest_on_stretch(
    curve: &Curve,
    range: (f64, f64),
    target: Point,
    tol: Tolerances,
) -> OgeomResult<f64> {
    const SCAN: u32 = 64;
    let at = |t: f64| -> OgeomResult<f64> { Ok(curve.point_at(t, tol)?.distance(target)) };
    let step = (range.1 - range.0) / f64::from(SCAN);
    let mut best = (range.0, at(range.0)?);
    for k in 1..=SCAN {
        let t = range.0 + step * f64::from(k);
        let d = at(t)?;
        if d < best.1 {
            best = (t, d);
        }
    }
    let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
    let (mut a, mut b) = ((best.0 - step.abs()).max(lo), (best.0 + step.abs()).min(hi));
    let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
    for _ in 0..80 {
        let (c, d) = (b - (b - a) * ratio, a + (b - a) * ratio);
        if at(c)? < at(d)? {
            b = d;
        } else {
            a = c;
        }
    }
    Ok(best.1.min(at(f64::midpoint(a, b))?))
}

/// The widest gap between `edge`'s curve and any of its pcurves at the
/// edge's own placement, measured by [`pcurve_gap`] against the edge's
/// stated tolerance; `None` for an edge with no curve.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the edge
/// or its geometry is not in `model`; as [`pcurve_gap`].
pub fn edge_pcurve_gap(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<Option<f64>> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        ogeom_core::ogeom_bail!(Dangling, "an edge is not in this model");
    };
    let Some(EdgeRepr::Curve3d {
        curve,
        location,
        range,
    }) = data.curve3d()
    else {
        return Ok(None);
    };
    let Some(curve) = model.geometry().curve(*curve) else {
        ogeom_core::ogeom_bail!(Dangling, "an edge names a curve not in this model");
    };
    let stated = data.tolerance.get();
    let mut widest: f64 = 0.0;
    for repr in &data.representations {
        let (sides, pcurve_range, surface, at) = match repr {
            EdgeRepr::PCurve {
                curve,
                range,
                surface,
                location,
            } => ([Some(*curve), None], *range, *surface, location),
            EdgeRepr::Seam {
                forward,
                reversed,
                range,
                surface,
                location,
            } => (
                [Some(*forward), Some(*reversed)],
                *range,
                *surface,
                location,
            ),
            _ => continue,
        };
        // A pcurve at another placement describes the edge where that
        // occurrence stands, not where the curve does.
        if at != location {
            continue;
        }
        let Some(surface) = model.geometry().surface(surface) else {
            ogeom_core::ogeom_bail!(Dangling, "an edge names a surface not in this model");
        };
        for id in sides.into_iter().flatten() {
            let Some(pcurve) = model.geometry().pcurve(id) else {
                ogeom_core::ogeom_bail!(Dangling, "an edge names a pcurve not in this model");
            };
            widest = widest.max(pcurve_gap(
                (curve, *range),
                (pcurve, pcurve_range),
                surface,
                stated,
                tol,
            )?);
        }
    }
    Ok(Some(widest))
}

/// Raise each of `edges`' tolerance, and its vertices', to how far its
/// pcurves stand off its curve, where that is wider than it states.
///
/// # Errors
///
/// As [`edge_pcurve_gap`].
pub fn state_pcurve_gaps_of(
    model: &mut Model,
    edges: &[Shape],
    tol: Tolerances,
) -> OgeomResult<()> {
    let mut wider = Vec::new();
    for edge in edges {
        let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let stated = data.tolerance.get();
        if let Some(gap) = edge_pcurve_gap(model, edge, tol)?
            && gap > stated
        {
            wider.push((edge.clone(), gap));
        }
    }
    for (edge, gap) in wider {
        model.widen(&edge, Tolerance::new(gap * (1.0 + MARGIN))?)?;
    }
    Ok(())
}

/// [`state_pcurve_gaps_of`] over every edge of `shape`.
///
/// # Errors
///
/// As [`edge_pcurve_gap`].
pub fn state_pcurve_gaps(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<()> {
    let edges = ogeom_topo::explore_unique(model, shape, ShapeType::Edge)?;
    state_pcurve_gaps_of(model, &edges, tol)
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use ogeom_geom::{BSpline2d, LineCurve, PlaneSurface};
    use ogeom_math::{KnotVector, Plane, Point2};

    const T: Tolerances = Tolerances::millimetres();

    /// A pcurve rising to a corner half way between two samples, along a
    /// straight edge on the plane: the samples either side stand short of
    /// the corner by a quarter of a percent, and the gap is the corner's
    /// height, where a piece of the edge sampled at its corner finds it.
    #[test]
    fn a_peak_between_samples_is_measured_at_its_top() {
        let curve: Curve = LineCurve::segment(Point::ORIGIN, Point::new(1.0, 0.0, 0.0), T)
            .unwrap()
            .into();
        let surface: SurfaceGeometry =
            PlaneSurface::new(Plane::new(ogeom_math::Frame::WORLD)).into();
        let (corner, height) = (128.5 / 256.0, 1e-3);
        let knots = KnotVector::new(vec![0.0, 0.0, corner, 1.0, 1.0], 1).unwrap();
        let control = vec![
            Point2::new(0.0, 0.0),
            Point2::new(corner, height),
            Point2::new(1.0, 0.0),
        ];
        let pcurve: PlanarCurve = BSpline2d::new(knots, control, T).unwrap().into();
        let gap = pcurve_gap(
            (&curve, (0.0, 1.0)),
            (&pcurve, (0.0, 1.0)),
            &surface,
            1e-7,
            T,
        )
        .unwrap();
        assert!(
            (gap - height).abs() < height * 1e-6,
            "{gap} against {height}"
        );
    }
}
