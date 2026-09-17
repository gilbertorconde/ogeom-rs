//! Pcurves for imported faces, shared by the exchange readers.
//!
//! The machinery itself lives in `ogeom_algo::pcurve_fit` — fitting a trim
//! by projection is geometry, not exchange — and the readers reach it
//! through this shim under the names they always used.
//!
//! What is written here is the part neither reader can do edge by edge:
//! *where in the chart* an image goes. A periodic chart offers a branch per
//! turn, every one describing the same points, and only continuity with the
//! edge before chooses between them. Both readers derive each image alone,
//! so both need the same walk afterwards.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{PlanarCurve, SurfaceGeometry};

pub(crate) use ogeom_algo::pcurve_fit::fit_projected_pcurve;

/// An image slid by whole periods until its start meets `previous`.
///
/// The only freedom a periodic chart leaves: every branch describes the
/// same points, and only continuity with the edge before it chooses.
pub(crate) fn shifted_to_meet(
    image: &PlanarCurve,
    start: f64,
    previous: Option<ogeom_math::Point2>,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<PlanarCurve> {
    use ogeom_geom::Curve2d as _;
    use ogeom_geom::Surface as _;
    let Some(previous) = previous else {
        return Ok(image.clone());
    };
    let ((ua, ub), (va, vb)) = surface.domain();
    let u_period = surface.is_periodic_u().then_some(ub - ua);
    let v_period = surface.is_periodic_v().then_some(vb - va);
    if u_period.is_none() && v_period.is_none() {
        return Ok(image.clone());
    }
    // Whole turns only, and only where the gap is plainly some number of
    // them. A branch chosen wrongly leaves a gap within rounding of a whole
    // period; anything else is the file saying something this does not
    // understand — a cone's wire that already stood half a turn open, say —
    // and half a period rounds to one, which would open it further. Left
    // alone, such a wire is exactly as it was.
    let whole = |gap: f64, period: f64| -> f64 {
        let turns = (gap / period).round();
        let left = period.mul_add(-turns, gap);
        if gap.abs() > period * 0.75 && left.abs() < period * 0.25 {
            turns * period
        } else {
            0.0
        }
    };
    let at = image.point_at(start, tol)?;
    let shift = ogeom_math::Vector2::new(
        u_period.map_or(0.0, |p| whole(previous.x - at.x, p)),
        v_period.map_or(0.0, |p| whole(previous.y - at.y, p)),
    );
    if shift.x == 0.0 && shift.y == 0.0 {
        return Ok(image.clone());
    }
    image.transformed(&ogeom_math::Transform2::translation(shift), tol)
}

/// The other column of a seam the wire walked only one way: a period over,
/// toward the middle of the chart.
pub(crate) fn seam_other_side(
    image: &PlanarCurve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<PlanarCurve> {
    use ogeom_geom::Curve2d as _;
    use ogeom_geom::Surface as _;
    let ((ua, ub), (va, vb)) = surface.domain();
    let a = image.point_at(range.0, tol)?;
    let b = image.point_at(range.1, tol)?;
    let runs_in_u = (b.x - a.x).abs() > (b.y - a.y).abs();
    let mid = image.point_at(f64::midpoint(range.0, range.1), tol)?;
    let shift = if runs_in_u && surface.is_periodic_v() {
        let span = vb - va;
        let d = if mid.y - va < span * 0.5 { span } else { -span };
        ogeom_math::Vector2::new(0.0, d)
    } else {
        let span = ub - ua;
        let d = if mid.x - ua < span * 0.5 { span } else { -span };
        ogeom_math::Vector2::new(d, 0.0)
    };
    image.transformed(&ogeom_math::Transform2::translation(shift), tol)
}
