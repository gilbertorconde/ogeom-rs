//! A B-spline surface's own boundary iso-curve has the chart's straight
//! side as its exact pcurve, also where the patch closes to a point along
//! another side and a projection near that point has no one foot.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::{
    BSplineSurface, Curve, Curve2d as _, Curve3d as _, Surface as _, SurfaceGeometry,
};
use ogeom_intersect::exact_pcurve_of;
use ogeom_math::{ControlGrid, KnotVector, Point};

const T: Tolerances = Tolerances::millimetres();

/// A bicubic patch whose last row of controls is one point, so one side of
/// its chart maps to that point.
fn closing_patch() -> BSplineSurface {
    let mut net = Vec::with_capacity(16);
    for i in 0..4_u8 {
        for j in 0..4_u8 {
            let (fi, fj) = (f64::from(i), f64::from(j));
            net.push(if i == 3 {
                Point::new(6.0, 1.0, 2.0)
            } else {
                Point::new(fi * 2.0, fj * 1.5, (fi * fj).sin())
            });
        }
    }
    let cubic = || KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 3).unwrap();
    BSplineSurface::new(cubic(), cubic(), &ControlGrid::new(net, 4, 4).unwrap(), T).unwrap()
}

/// Each boundary iso-curve's pcurve lifts onto the curve at its own
/// parameter, to rounding, the collapsed side among them.
#[test]
fn a_boundary_iso_curve_reads_as_the_chart_s_side() {
    let patch = closing_patch();
    let surface = SurfaceGeometry::BSpline(patch.clone());
    let ((u0, u1), (v0, v1)) = patch.domain();
    let isos = [
        patch.iso_u_curve(u0, T).unwrap(),
        patch.iso_u_curve(u1, T).unwrap(),
        patch.iso_v_curve(v0, T).unwrap(),
        patch.iso_v_curve(v1, T).unwrap(),
    ];
    for iso in isos {
        let curve = Curve::BSpline(iso);
        let pcurve = exact_pcurve_of(&curve, &surface, T).expect("an exact pcurve");
        let (a, b) = curve.domain();
        for k in 0..=16 {
            let t = a + (b - a) * f64::from(k) / 16.0;
            let uv = pcurve.point_at(t, T).unwrap();
            let lifted = patch.point_at(uv.x, uv.y, T).unwrap();
            let on = curve.point_at(t, T).unwrap();
            assert!(lifted.distance(on) < 1e-12, "{lifted:?} against {on:?}");
        }
    }
}

/// An iso-curve inside the domain, or any other curve, is not read so.
#[test]
fn an_inner_iso_curve_is_not_read_as_a_side() {
    let patch = closing_patch();
    let surface = SurfaceGeometry::BSpline(patch.clone());
    let inner = Curve::BSpline(patch.iso_u_curve(0.5, T).unwrap());
    assert!(exact_pcurve_of(&inner, &surface, T).is_none());
}
