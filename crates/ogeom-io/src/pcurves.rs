//! Pcurves for imported faces, shared by the exchange readers.
//!
//! The machinery itself lives in `ogeom_algo::pcurve_fit` (fitting a trim
//! by projection is geometry, not exchange), and the readers reach it
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
    // understand (a cone's wire that already stood half a turn open, say),
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
///
/// Which axis to step along is the surface's business, not the curve's. A
/// seam lies *on* the join, so it runs along the direction the surface does
/// not close and stands still in the one it does, and closure, not
/// periodicity, is the test, the same distinction a skinned wall forced
/// everywhere else in this module. A patch closed in `v` over `(-π, π)`,
/// whose `u` is a knot range that closes on nothing, had its second column
/// stepped a `u` span sideways and off the chart entirely: the face then
/// drew as one flat triangle across itself, over whatever it was a boss for.
///
/// Where the surface closes on neither axis the edge is a slit rather than
/// a seam (one curve, walked twice), and the same column serves both ways.
pub(crate) fn seam_other_side(
    image: &PlanarCurve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> OgeomResult<PlanarCurve> {
    use ogeom_geom::Curve2d as _;
    use ogeom_geom::Surface as _;
    let ((ua, ub), (va, vb)) = surface.domain();
    let closed_u = surface.is_periodic_u() || surface.is_closed_u(tol);
    let closed_v = surface.is_periodic_v() || surface.is_closed_v(tol);
    let a = image.point_at(range.0, tol)?;
    let b = image.point_at(range.1, tol)?;
    let runs_in_u = (b.x - a.x).abs() >= (b.y - a.y).abs();
    let over_v = match (closed_u, closed_v) {
        (false, false) => return Ok(image.clone()),
        (true, false) => false,
        (false, true) => true,
        // Closed both ways: a torus. The seam stands still in the axis it
        // is a seam of, which is the one it does not run along.
        (true, true) => runs_in_u,
    };
    let mid = image.point_at(f64::midpoint(range.0, range.1), tol)?;
    let shift = if over_v {
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

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;
    use ogeom_geom::Curve2d as _;
    use ogeom_geom::{BSplineSurface, CylinderSurface, Line2d, PlaneSurface};
    use ogeom_math::{ControlGrid, Cylinder, Frame, KnotVector, Plane, Point, Point2};

    const T: Tolerances = Tolerances::millimetres();

    /// A seam steps across the join the surface actually has.
    ///
    /// A boss around a bolt hole is a patch closed in `v` whose `u` is a
    /// knot range closing on nothing. Stepping its second column a `u` span
    /// sideways puts it off the chart, the wire never closes, and the face
    /// collapses to a single triangle laid across whatever it was a boss
    /// for. Closure decides the axis, not periodicity and not `u` by
    /// default: the same distinction a skinned wall forces everywhere.
    #[test]
    fn a_seam_steps_across_the_axis_its_surface_closes_on() {
        use ogeom_geom::Surface as _;
        let pi = core::f64::consts::PI;

        // A tube: one straight run in `u`, a closed square ring in `v`.
        let ring = [
            Point::new(0.0, 1.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
            Point::new(0.0, -1.0, 0.0),
            Point::new(0.0, 0.0, -1.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let mut control = Vec::new();
        for i in 0..2 {
            for p in &ring {
                control.push(Point::new(f64::from(i), p.y, p.z));
            }
        }
        let tube: SurfaceGeometry = BSplineSurface::new(
            KnotVector::clamped_uniform(1, 2).unwrap(),
            KnotVector::new(vec![-pi, -pi, -pi / 2.0, 0.0, pi / 2.0, pi, pi], 1).unwrap(),
            &ControlGrid::new(control, 2, 5).unwrap(),
            T,
        )
        .unwrap()
        .into();
        assert!(!tube.is_closed_u(T) && tube.is_closed_v(T), "closed in v");

        // The seam runs along `v = -pi`, the whole way across `u`.
        let along: PlanarCurve = Line2d::segment(Point2::new(0.0, -pi), Point2::new(1.0, -pi), T)
            .unwrap()
            .into();
        let other = seam_other_side(&along, (0.0, 1.0), &tube, T).unwrap();
        let ends = (
            other.point_at(0.0, T).unwrap(),
            other.point_at(1.0, T).unwrap(),
        );
        assert!(
            (ends.0.y - pi).abs() < 1e-12 && (ends.1.y - pi).abs() < 1e-12,
            "stepped to the other v: {ends:?}"
        );
        assert!(
            (ends.0.x - 0.0).abs() < 1e-12 && (ends.1.x - 1.0).abs() < 1e-12,
            "and stayed where it was in u: {ends:?}"
        );

        // A cylinder closes the other way round, and the seam runs along
        // `u = 0` up its height; the step is a turn in `u`.
        let drum: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 3.0, T).unwrap(), (0.0, 10.0))
                .unwrap()
                .into();
        let up: PlanarCurve = Line2d::segment(Point2::new(0.0, 0.0), Point2::new(0.0, 10.0), T)
            .unwrap()
            .into();
        let over = seam_other_side(&up, (0.0, 1.0), &drum, T).unwrap();
        let at = over.point_at(0.0, T).unwrap();
        assert!(
            (at.x - core::f64::consts::TAU).abs() < 1e-12 && at.y.abs() < 1e-12,
            "a turn round in u: {at:?}"
        );

        // A surface that closes on neither axis has no other side: the
        // edge is a slit, one curve walked twice, and the column stands.
        let flat: SurfaceGeometry = PlaneSurface::over(
            Plane::through(Point::ORIGIN, ogeom_math::Direction::Z),
            (-5.0, 5.0),
            (-5.0, 5.0),
        )
        .unwrap()
        .into();
        let slit = seam_other_side(&along, (0.0, 1.0), &flat, T).unwrap();
        assert_eq!(
            slit.point_at(0.0, T).unwrap(),
            along.point_at(0.0, T).unwrap()
        );
    }
}
