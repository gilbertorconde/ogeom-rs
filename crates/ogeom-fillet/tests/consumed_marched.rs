//! A round of a bore's rim on a drum taking a narrow bevel beside it off
//! whole: the marched counterpart of a rim's round over a chamfer.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{chamfer_edges, fillet_edges};
use ogeom_geom::{Curve, Curve3d as _, LineCurve};
use ogeom_math::{Axis, Direction, Frame, KnotVector, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

/// A drum of radius 10 and height 30 with a bore of radius 3 across it
/// along y at height 15.
fn holed_drum(m: &mut Model) -> Shape {
    let drum = ogeom_algo::make_cylinder(m, at(0.0, 0.0, 0.0), 10.0, 30.0, T)
        .unwrap()
        .shape;
    let f = Frame::new(Point::new(0.0, -20.0, 15.0), Direction::Y, Direction::X, T).unwrap();
    let bore = ogeom_algo::make_cylinder(m, f, 3.0, 40.0, T).unwrap().shape;
    ogeom_bool::cut(m, &drum, &bore, T).unwrap().shape
}

/// The edge of `s` passing nearest `p`.
fn edge_near(m: &Model, s: &Shape, p: [f64; 3]) -> Shape {
    let q = Point::new(p[0], p[1], p[2]);
    let d = |e: &Shape| {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, e, T).unwrap() else {
            return f64::INFINITY;
        };
        (0..=64)
            .map(|i| {
                c.point_at(a + (b - a) * f64::from(i) / 64.0, T)
                    .unwrap()
                    .distance(q)
            })
            .fold(f64::INFINITY, f64::min)
    };
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| d(a).total_cmp(&d(b)))
        .unwrap()
}

fn vol(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// The faces of `s` holding `e`.
fn faces_at(m: &Model, s: &Shape, e: &Shape) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            explore_unique(m, f, ShapeType::Edge)
                .unwrap()
                .iter()
                .any(|x| x.is_same(e))
        })
        .collect()
}

/// The bore's rim at y > 0 rounded at `radius` on the plain drum, and on
/// the drum with that rim first bevelled 0.5, rounded from the strip's
/// edge nearest `from`: both volumes, the second checked, the bevel gone.
fn rounded(radius: f64, from: [f64; 3]) -> (f64, f64) {
    let mut m = Model::new();
    let s = holed_drum(&mut m);
    let rim = edge_near(&m, &s, [0.0, 9.5, 18.0]);
    let plain = fillet_edges(&mut m, &s, std::slice::from_ref(&rim), radius, T).unwrap();
    let plain = vol(&m, &plain.shape);
    let bevelled = chamfer_edges(&mut m, &s, &[rim], 0.5, T).unwrap().shape;
    let edge = edge_near(&m, &bevelled, from);
    let strip = faces_at(&m, &bevelled, &edge)
        .into_iter()
        .find(|f| {
            matches!(
                m.node(f)
                    .and_then(|n| n.data().as_face())
                    .and_then(|d| m.geometry().surface(d.surface)),
                Some(ogeom_geom::SurfaceGeometry::BSpline(_))
            )
        })
        .expect("the bevel beside the edge");
    let built = fillet_edges(&mut m, &bevelled, &[edge], radius, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!(
        !explore_unique(&m, &built.shape, ShapeType::Face)
            .unwrap()
            .iter()
            .any(|f| f.is_same(&strip)),
        "the bevel is gone"
    );
    (plain, vol(&m, &built.shape))
}

fn close(v: f64, expected: f64) {
    assert!(
        (v - expected).abs() <= 1e-6 * expected,
        "{v} against {expected}"
    );
}

#[test]
fn a_bevelled_bore_rim_rounds_from_the_drum_edge_as_the_plain_rim() {
    // The drum / bevel edge, half a millimetre down the drum from the rim.
    let (plain, consumed) = rounded(1.0, [0.0, 9.9, 18.6]);
    close(consumed, plain);
}

#[test]
fn a_bevelled_bore_rim_rounds_from_the_bore_edge_as_the_plain_rim() {
    // The bevel / bore edge, half a millimetre into the bore from the rim.
    let (plain, consumed) = rounded(1.0, [0.0, 9.0, 18.0]);
    close(consumed, plain);
}

#[test]
fn a_bevel_wider_than_the_round_s_setback_is_still_refused() {
    // A bevel of 1.5 each way against a ball of 1: the ball touching the
    // drum and the bore stands within the bevel, off both faces.
    let mut m = Model::new();
    let s = holed_drum(&mut m);
    let rim = edge_near(&m, &s, [0.0, 9.5, 18.0]);
    let bevelled = chamfer_edges(&mut m, &s, &[rim], 1.5, T).unwrap().shape;
    let edge = edge_near(&m, &bevelled, [0.0, 9.9, 19.5]);
    assert!(fillet_edges(&mut m, &bevelled, &[edge], 1.0, T).is_err());
}

/// The S-shaped crease's profile in `y = 0`: from the origin along `x`,
/// down by 1 and on along `x` to `x = 10`, its tangent along `x` at both
/// ends.
fn spline() -> Curve {
    let at = |x: f64, z: f64| Point::new(x, 0.0, z);
    ogeom_geom::BSplineCurve::new(
        KnotVector::clamped_uniform(3, 4).unwrap(),
        vec![at(0.0, 0.0), at(3.0, 0.0), at(7.0, -1.0), at(10.0, -1.0)],
        T,
    )
    .unwrap()
    .into()
}

/// A slab 6 thick along `y`: its profile in `y = 0` under the spline,
/// down to `z = -8`, between the end walls `x = 0` and `x = 10`, each
/// square to the spline's end. Returns the slab, the crease along the
/// spline's top, its top face and its wall.
fn slab(m: &mut Model) -> (Shape, Shape, Shape, Shape) {
    let at = |x: f64, z: f64| Point::new(x, 0.0, z);
    let corners: Vec<Shape> = [at(0.0, 0.0), at(10.0, -1.0), at(10.0, -8.0), at(0.0, -8.0)]
        .into_iter()
        .map(|p| ogeom_algo::make_vertex(m, p).shape)
        .collect();
    let line = |m: &mut Model, from: (&Shape, Point), to: (&Shape, Point)| {
        let curve: Curve =
            LineCurve::new(Axis::new(from.1, Direction::new(to.1 - from.1, T).unwrap())).into();
        ogeom_algo::make_edge_between(m, curve, (0.0, from.1.distance(to.1)), from.0, to.0, T)
            .unwrap()
            .shape
    };
    let edges = vec![
        ogeom_algo::make_edge_between(m, spline(), (0.0, 1.0), &corners[0], &corners[1], T)
            .unwrap()
            .shape,
        line(
            m,
            (&corners[1], at(10.0, -1.0)),
            (&corners[2], at(10.0, -8.0)),
        ),
        line(
            m,
            (&corners[2], at(10.0, -8.0)),
            (&corners[3], at(0.0, -8.0)),
        ),
        line(m, (&corners[3], at(0.0, -8.0)), (&corners[0], at(0.0, 0.0))),
    ];
    let wire = ogeom_algo::make_wire(m, &edges, T).unwrap().shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Y),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(m, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let solid = ogeom_algo::make_prism(m, &face, Vector::new(0.0, 6.0, 0.0), T)
        .unwrap()
        .shape;
    let crease = explore_unique(m, &solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let b = ogeom_algo::shape_bounds(m, e, T).unwrap();
            let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
            (lo.y - 6.0).abs() < 1e-6
                && (hi.y - 6.0).abs() < 1e-6
                && lo.x.abs() < 1e-6
                && (hi.x - 10.0).abs() < 1e-6
                && hi.z - lo.z > 0.5
        })
        .expect("the crease along the spline's top");
    let beside: Vec<Shape> = explore_unique(m, &solid, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            explore_unique(m, f, ShapeType::Edge)
                .unwrap()
                .iter()
                .any(|x| x.is_same(&crease))
        })
        .collect();
    let is_plane = |f: &Shape| {
        matches!(
            m.node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| m.geometry().surface(d.surface)),
            Some(ogeom_geom::SurfaceGeometry::Plane(_))
        )
    };
    let top = beside
        .iter()
        .find(|f| is_plane(f))
        .expect("the top")
        .clone();
    let wall = beside
        .iter()
        .find(|f| !is_plane(f))
        .expect("the wall")
        .clone();
    (solid, crease, top, wall)
}

/// The slab's spline crease turns through nothing between its ends, so a
/// section swept along it takes off its area times the crease's length:
/// the plain crease rounded at `radius` takes off `r² (1 - π/4) L`.
fn plain_round(before: f64, radius: f64) -> f64 {
    use ogeom_geom::Curve3d as _;
    let curve = spline();
    let nodes = [
        (-0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
        (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
    ];
    let mut length = 0.0;
    for k in 0..256 {
        let (a, b) = (f64::from(k) / 256.0, f64::from(k + 1) / 256.0);
        for (x, w) in nodes {
            let t = (b - a).mul_add(x, a + b) * 0.5;
            length += w * (b - a) * 0.5 * curve.d1_at(t, T).unwrap().magnitude();
        }
    }
    before - radius * radius * (1.0 - core::f64::consts::FRAC_PI_4) * length
}

/// The slab with its crease bevelled 0.5, rounded at `radius` from the
/// bevel's edge nearest `from`: checked, the bevel gone; the volume and the
/// slab's before the bevel.
fn slab_rounded(radius: f64, from: [f64; 3]) -> (f64, f64) {
    let mut m = Model::new();
    let (solid, crease, _, _) = slab(&mut m);
    let before = vol(&m, &solid);
    let bevelled = chamfer_edges(&mut m, &solid, &[crease], 0.5, T)
        .unwrap()
        .shape;
    let edge = edge_near(&m, &bevelled, from);
    let strip = faces_at(&m, &bevelled, &edge)
        .into_iter()
        .find(|f| {
            matches!(
                m.node(f)
                    .and_then(|n| n.data().as_face())
                    .and_then(|d| m.geometry().surface(d.surface)),
                Some(ogeom_geom::SurfaceGeometry::BSpline(_))
            )
        })
        .expect("the bevel beside the edge");
    let built = fillet_edges(&mut m, &bevelled, &[edge], radius, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!(
        !explore_unique(&m, &built.shape, ShapeType::Face)
            .unwrap()
            .iter()
            .any(|f| f.is_same(&strip)),
        "the bevel is gone"
    );
    (before, vol(&m, &built.shape))
}

#[test]
fn a_bevelled_open_crease_rounds_as_the_plain_crease_from_either_bevel_edge() {
    // A ball of 2 or 3 sets back `r tan(π/8)`, 0.83 or 1.24, across the top
    // and the bevel, past the bevel's own width of 0.71: no round of the
    // bevel's edge seats, and the round takes the bevel off whole.
    for radius in [2.0, 3.0] {
        // The top / bevel edge over the spline's middle, and the bevel /
        // wall edge.
        for from in [[5.0, 6.0, -1.0], [5.0, 5.5, -0.5]] {
            let (before, v) = slab_rounded(radius, from);
            let expected = plain_round(before, radius);
            assert!(
                (v - expected).abs() <= 1e-6 * expected,
                "r {radius} from {from:?}: {v} against {expected}"
            );
        }
    }
}
