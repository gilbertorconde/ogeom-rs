use super::pcurves::paced;
use super::*;

const T: Tolerances = Tolerances::millimetres();

/// A pcurve read at an affine map of the parameter is the same curve
/// at the mapped parameter, for each kind and each way round: a line
/// at a pace of one and at another, a circle reversed and turned, an
/// ellipse turned a half turn, a spline reversed and slowed. Rounding
/// in the turned frame's sines is all that stands between the two.
#[test]
fn a_pcurve_paced_is_the_same_curve_at_the_mapped_parameter() {
    use ogeom_math::{Direction2, Frame2, Point2};
    let frame = Frame2::new(Point2::new(1.0, -2.0), Direction2::from_angle(0.3));
    let curves: Vec<PlanarCurve> = vec![
        ogeom_geom::Line2d::over(
            ogeom_math::Axis2::new(Point2::new(0.5, 1.0), Direction2::from_angle(1.1)),
            -10.0,
            10.0,
        )
        .unwrap()
        .into(),
        ogeom_geom::Circle2d::new(ogeom_math::Circle2::new(frame, 2.0, T).unwrap()).into(),
        ogeom_geom::Ellipse2d::new(ogeom_math::Ellipse2::new(frame, 3.0, 1.0, T).unwrap()).into(),
        ogeom_geom::BSpline2d::new(
            KnotVector::new(vec![0.0, 0.0, 0.0, 2.0, 5.0, 5.0, 5.0], 2).unwrap(),
            vec![
                Point2::new(0.0, 0.0),
                Point2::new(1.0, 3.0),
                Point2::new(4.0, 2.0),
                Point2::new(6.0, -1.0),
            ],
            T,
        )
        .unwrap()
        .into(),
    ];
    let pi = core::f64::consts::PI;
    let maps: [&[(f64, f64)]; 4] = [
        &[(1.0, 2.5), (-1.0, 0.75), (0.4, -1.0)],
        &[(1.0, 0.7), (-1.0, 0.7), (-1.0, -2.0)],
        &[(1.0, pi), (-1.0, pi), (-1.0, 0.0)],
        &[(0.5, 1.0), (-2.0, 4.0)],
    ];
    for (curve, maps) in curves.iter().zip(maps) {
        for &(a, b) in maps {
            let paced = paced(curve, a, b, T).unwrap();
            let (lo, hi) = curve.domain();
            let ((from, to), at) = if curve.is_periodic() {
                ((0.0, 6.0), 6)
            } else {
                let (x, y) = ((lo - b) / a, (hi - b) / a);
                ((x.min(y), x.max(y)), 8)
            };
            for i in 0..=at {
                let t = from + (to - from) * f64::from(i) / f64::from(at);
                let want = curve.point_at(a.mul_add(t, b), T).unwrap();
                let got = paced.point_at(t, T).unwrap();
                assert!(
                    got.distance(want) <= 1e-12,
                    "{curve:?} at {a}·t + {b}, t = {t}: {got:?} against {want:?}"
                );
            }
        }
    }
    // An ellipse turns only by whole half turns, a conic only at a pace
    // of one.
    assert!(paced(&curves[2], 1.0, 0.5, T).is_none());
    assert!(paced(&curves[1], 2.0, 0.0, T).is_none());
}

/// An edge's window between the feet of vertices standing on the
/// curve short of its ends: taken in the curve's order, refused where
/// the vertices are the ends already, run against the curve, or off it.
#[test]
fn an_edges_window_is_taken_between_its_vertices_feet() {
    let spline: Curve = BSplineCurve::new(
        KnotVector::clamped_uniform(3, 6).unwrap(),
        vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 2.0, 0.5),
            Point::new(2.5, 1.0, -0.5),
            Point::new(4.0, 3.0, 1.0),
            Point::new(5.0, 0.5, 0.0),
            Point::new(6.0, 2.0, 2.0),
        ],
        T,
    )
    .unwrap()
    .into();
    let at = |t: f64| spline.point_at(t, T).unwrap();
    let (a, b) = window_between_feet(&spline, at(0.3), at(0.7), T).unwrap();
    assert!(
        (a - 0.3).abs() < 1e-9 && (b - 0.7).abs() < 1e-9,
        "{a} .. {b}"
    );
    let (a, b) = window_between_feet(&spline, at(0.0), at(0.7), T).unwrap();
    assert!(a.abs() < 1e-9 && (b - 0.7).abs() < 1e-9, "{a} .. {b}");
    assert!(window_between_feet(&spline, at(0.0), at(1.0), T).is_none());
    assert!(window_between_feet(&spline, at(0.7), at(0.3), T).is_none());
    let off = at(0.3) + Vector::new(0.0, 0.0, 0.5);
    assert!(window_between_feet(&spline, off, at(0.7), T).is_none());
}

/// A parallel of a cone read as a circle a few tenths of a micron wider
/// than the cone there, its range running past one turn: the exact
/// pcurve is stated over the whole range, and the gap it stands off
/// the circle is measured along all of it.
#[test]
fn a_wrapped_parallel_states_its_gap() {
    let cone: SurfaceGeometry = ConeSurface::new(
        Cone::new(Frame::WORLD, 2.0, core::f64::consts::FRAC_PI_4, T).unwrap(),
        (-10.0, 10.0),
    )
    .unwrap()
    .into();
    let wider = 3e-7;
    let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 2.0 + wider, T).unwrap()).into();
    let range = (3.0, 8.5);
    let PreparedPcurve::Exact { curve, off } = derive_pcurve(&circle, range, &cone, &[], T) else {
        panic!("a parallel has an exact pcurve");
    };
    assert!(curve.point_at(range.1, T).is_ok(), "{curve:?}");
    assert!((off - wider).abs() < wider * 1e-3, "{off}");
}
