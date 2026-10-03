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
/// from the lifted point to the nearest point of the curve's stretch.
/// Within `stated`, the answer is only known to be within it.
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
    let mut widest: f64 = 0.0;
    for i in 0..=SAMPLES {
        let f = f64::from(i) / f64::from(SAMPLES);
        let uv = pcurve.point_at(pcurve_range.0 + (pcurve_range.1 - pcurve_range.0) * f, tol)?;
        let Ok(lifted) = surface.point_at(uv.x, uv.y, tol) else {
            continue;
        };
        let mut gap = curve
            .point_at(range.0 + (range.1 - range.0) * f, tol)?
            .distance(lifted);
        if gap > widest.max(stated) {
            gap = gap.min(nearest_on_stretch(curve, range, lifted, tol)?);
        }
        widest = widest.max(gap);
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
