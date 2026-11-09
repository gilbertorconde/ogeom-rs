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

/// The slab's crease rounded at `radius`, and the slab with the crease
/// first bevelled 0.5 rounded from the bevel's edge on the top face: both
/// volumes, the second checked, the bevel gone.
fn slab_rounded(radius: f64) -> (f64, f64) {
    let mut m = Model::new();
    let (solid, crease, _, _) = slab(&mut m);
    let plain = fillet_edges(&mut m, &solid, std::slice::from_ref(&crease), radius, T).unwrap();
    let plain = vol(&m, &plain.shape);
    let bevelled = chamfer_edges(&mut m, &solid, &[crease], 0.5, T).unwrap().shape;
    // The top / bevel edge, half a millimetre in from the crease over the
    // spline's middle.
    let edge = edge_near(&m, &bevelled, [5.0, 6.0, -1.0]);
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

#[test]
fn a_bevelled_open_crease_rounds_as_the_plain_crease() {
    // A ball of 2 sets back 2 tan(22.5 degrees), 0.83, across the top and
    // the bevel, past the bevel's own width of 0.71: no round of the
    // bevel's edge seats, and the round takes the bevel off whole.
    let (plain, consumed) = slab_rounded(2.0);
    close(consumed, plain);
}

#[test]
fn zz_open_slab() {
    for r in [2.0, 3.0] {
        let mut m = Model::new();
        let (solid, crease, _, _) = slab(&mut m);
        let plain = fillet_edges(&mut m, &solid, std::slice::from_ref(&crease), r, T).map(|b| vol(&m, &b.shape));
        let bevelled = chamfer_edges(&mut m, &solid, &[crease], 0.5, T).unwrap().shape;
        for (name, p) in [("top", [5.0, 6.0, -1.0]), ("wall", [5.0, 5.5, -0.5])] {
            let e = edge_near(&m, &bevelled, p);
            let got = fillet_edges(&mut m, &bevelled, &[e], r, T).map(|b| vol(&m, &b.shape));
            eprintln!("r {r} {name}: plain {:?} got {:?}", plain.as_ref().map_err(|e| e.to_string()), got.map_err(|e| e.to_string()));
        }
    }
}

/// A slab 6 thick along `y` whose profile in `y = 0` is the half ellipse
/// of semi-axes 4 along `x` and 3 along `z` over the base `z = 0`: its top
/// crease the elliptic arc, square to the base at both ends.
fn elliptic_slab(m: &mut Model) -> (Shape, Shape) {
    use ogeom_geom::EllipseCurve;
    use ogeom_math::Ellipse;
    let frame = Frame::new(Point::ORIGIN, Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(), Direction::X, T).unwrap();
    let arc: Curve = EllipseCurve::new(Ellipse::new(frame, 4.0, 3.0, T).unwrap()).into();
    let a = arc.point_at(0.0, T).unwrap();
    let b = arc.point_at(core::f64::consts::PI, T).unwrap();
    let va = ogeom_algo::make_vertex(m, a).shape;
    let vb = ogeom_algo::make_vertex(m, b).shape;
    let top = ogeom_algo::make_edge_between(m, arc, (0.0, core::f64::consts::PI), &va, &vb, T).unwrap().shape;
    let base: Curve = LineCurve::new(Axis::new(b, Direction::new(a - b, T).unwrap())).into();
    let base = ogeom_algo::make_edge_between(m, base, (0.0, a.distance(b)), &vb, &va, T).unwrap().shape;
    let wire = ogeom_algo::make_wire(m, &[top, base], T).unwrap().shape;
    let plane = ogeom_geom::PlaneSurface::over(Plane::through(Point::ORIGIN, Direction::Y), (-100.0, 100.0), (-100.0, 100.0)).unwrap();
    let face = ogeom_algo::make_face(m, plane.into(), &[wire], T).unwrap().shape;
    let solid = ogeom_algo::make_prism(m, &face, Vector::new(0.0, 6.0, 0.0), T).unwrap().shape;
    let crease = edge_near(m, &solid, [0.0, 6.0, 3.0]);
    (solid, crease)
}

#[test]
fn zz_ellipse() {
    let r = 1.5;
    let mut m = Model::new();
    let (solid, crease) = elliptic_slab(&mut m);
    eprintln!("solid vol {}", vol(&m, &solid));
    let plain = fillet_edges(&mut m, &solid, std::slice::from_ref(&crease), r, T).map(|b| vol(&m, &b.shape));
    let ch = chamfer_edges(&mut m, &solid, std::slice::from_ref(&crease), 0.3, T).map(|b| b.shape);
    eprintln!("plain {:?} chamfer {:?}", plain.as_ref().map_err(|e| e.to_string()), ch.as_ref().map(|_| ()).map_err(|e| e.to_string()));
    if let Ok(ch) = ch {
        for (name, p) in [("top", [0.0, 6.0, 2.7]), ("wall", [0.0, 5.7, 3.0])] {
            let e = edge_near(&m, &ch, p);
            let got = fillet_edges(&mut m, &ch, &[e], r, T).map(|b| (vol(&m, &b.shape), ogeom_algo::check(&m, &b.shape, T).unwrap().is_valid()));
            eprintln!("{name}: {:?}", got.map_err(|e| e.to_string()));
        }
    }
}
