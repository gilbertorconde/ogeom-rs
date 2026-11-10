//! Rounding the edge beside a narrow bevel on a marched seat.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{chamfer_edges, fillet_edges};
use ogeom_geom::{Curve, Curve3d as _, LineCurve};
use ogeom_math::{Axis, Direction, KnotVector, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn vol(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(
        m,
        s,
        ogeom_mesh::Deflection {
            chord: 1e-4,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
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

/// The spline's length and its tangent's turn between its ends, by
/// Gauss–Legendre over 256 pieces.
fn length_and_turn() -> (f64, f64) {
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
    let heading = |t: f64| {
        let d = curve.d1_at(t, T).unwrap();
        d.z.atan2(d.x)
    };
    (length, heading(1.0) - heading(0.0))
}

/// A ball across the corner between the slab's top and its bevel, 135
/// degrees inside, sets back `r tan(π/8)` on each: within the bevel's width
/// of `0.5 √2` the round is the ball's, and the corner it takes off is the
/// section `r² (tan(π/8) - π/8)` swept along the top's crease, which, a
/// constant offset of a crease turning through nothing, is as long as it.
/// Past that width no ball touches the bevel itself: the round is not
/// seated on the bevel's surface continued past it, and takes the bevel off
/// whole instead, as the plain crease's round (`consumed_marched.rs`).
#[test]
fn a_round_beside_a_bevel_seats_on_the_bevel_or_takes_it_off_whole() {
    let (length, turn) = length_and_turn();
    assert!(turn.abs() < 1e-12, "{turn}");
    let eighth = core::f64::consts::FRAC_PI_8;
    for (radius, fits) in [(1.0_f64, true), (3.0, false)] {
        assert_eq!(radius * eighth.tan() < 0.5 * 2.0_f64.sqrt(), fits);
        let mut m = Model::new();
        let (solid, crease, _, _) = slab(&mut m);
        let before = vol(&m, &solid);
        let bevelled = chamfer_edges(&mut m, &solid, &[crease], 0.5, T)
            .unwrap()
            .shape;
        // The top / bevel edge, half a millimetre in from the crease over
        // the spline's middle.
        let edge = edge_near(&m, &bevelled, [5.0, 6.0, -1.0]);
        let built = fillet_edges(&mut m, &bevelled, &[edge], radius, T).unwrap();
        let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        let expected = if fits {
            before - 0.5 * 0.5 * 0.5 * length - radius * radius * (eighth.tan() - eighth) * length
        } else {
            before - radius * radius * (1.0 - core::f64::consts::FRAC_PI_4) * length
        };
        let v = vol(&m, &built.shape);
        assert!(
            (v - expected).abs() <= 1e-6 * expected,
            "{v} against {expected}"
        );
    }
}
