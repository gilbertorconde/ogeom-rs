//! Chamfers on an open crease only the marching machinery seats: a
//! spline's top edge along its wall, between two end walls.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{Chamfer, chamfer_edges_with};
use ogeom_geom::{Curve, Curve3d as _, LineCurve};
use ogeom_math::{Axis, Direction, KnotVector, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

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

/// A slab 2 thick along `y`: its profile in `y = 0` under the spline,
/// down to `z = -5`, between the end walls `x = 0` and `x = 10`, each
/// square to the spline's end. Returns the slab, the crease along the
/// spline's top, its top face and its wall.
fn slab(m: &mut Model) -> (Shape, Shape, Shape, Shape) {
    let at = |x: f64, z: f64| Point::new(x, 0.0, z);
    let corners: Vec<Shape> = [at(0.0, 0.0), at(10.0, -1.0), at(10.0, -5.0), at(0.0, -5.0)]
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
            (&corners[2], at(10.0, -5.0)),
        ),
        line(
            m,
            (&corners[2], at(10.0, -5.0)),
            (&corners[3], at(0.0, -5.0)),
        ),
        line(m, (&corners[3], at(0.0, -5.0)), (&corners[0], at(0.0, 0.0))),
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
    let solid = ogeom_algo::make_prism(m, &face, Vector::new(0.0, 2.0, 0.0), T)
        .unwrap()
        .shape;
    let crease = explore_unique(m, &solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let b = ogeom_algo::shape_bounds(m, e, T).unwrap();
            let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
            (lo.y - 2.0).abs() < 1e-6
                && (hi.y - 2.0).abs() < 1e-6
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

/// The slab's crease chamfered by `spec`, checked, against the volume the
/// section's right triangle of legs `on_top` and `on_wall` sweeps along
/// the crease: its area times the crease's length, the crease turning
/// through nothing between its ends, so the triangle's centroid stands
/// off it for no change in the length swept.
fn chamfered(spec: impl Fn(&Shape, &Shape) -> Chamfer, on_top: f64, on_wall: f64) {
    let (length, turn) = length_and_turn();
    assert!(turn.abs() < 1e-12, "{turn}");
    let mut m = Model::new();
    let (solid, crease, top, wall) = slab(&mut m);
    let before = vol(&m, &solid);
    let built =
        chamfer_edges_with(&mut m, &solid, &[(crease.clone(), spec(&top, &wall))], T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!(built.history.is_deleted(&crease));
    let expected = before - 0.5 * on_top * on_wall * length;
    let v = vol(&m, &built.shape);
    assert!(
        (v - expected).abs() <= 1e-6 * expected,
        "{v} against {expected}"
    );
}

#[test]
fn an_open_spline_crease_chamfers_symmetric() {
    chamfered(|_, _| Chamfer::Symmetric(0.5), 0.5, 0.5);
}

#[test]
fn an_open_spline_crease_chamfers_by_two_distances() {
    chamfered(
        |top, _| Chamfer::Distances {
            face: top.clone(),
            on_face: 0.8,
            on_other: 0.4,
        },
        0.8,
        0.4,
    );
}

#[test]
fn an_open_spline_crease_chamfers_by_distance_and_angle() {
    // The corner is square: the leg down the wall is the tangent's.
    let angle = 35.0_f64.to_radians();
    chamfered(
        |top, _| Chamfer::Angle {
            face: top.clone(),
            distance: 0.6,
            angle,
        },
        0.6,
        0.6 * angle.tan(),
    );
}
