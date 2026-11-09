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

/// The slab with its spline's crease in two edges over the spline's two
/// halves, meeting tangentially at its middle: the slab and the two.
fn split_slab(m: &mut Model) -> (Shape, Vec<Shape>) {
    let at = |x: f64, z: f64| Point::new(x, 0.0, z);
    let middle = spline().point_at(0.5, T).unwrap();
    let corners: Vec<Shape> = [
        at(0.0, 0.0),
        at(10.0, -1.0),
        at(10.0, -5.0),
        at(0.0, -5.0),
        middle,
    ]
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
        ogeom_algo::make_edge_between(m, spline(), (0.0, 0.5), &corners[0], &corners[4], T)
            .unwrap()
            .shape,
        ogeom_algo::make_edge_between(m, spline(), (0.5, 1.0), &corners[4], &corners[1], T)
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
    let creases = edges_where(m, &solid, |a, mid, b| {
        (mid.y - 2.0).abs() < 1e-9
            && (a.y - 2.0).abs() < 1e-9
            && (b.y - 2.0).abs() < 1e-9
            && (a.z - b.z).abs() > 0.1
            && (a.x - b.x).abs() > 1.0
    });
    (solid, creases)
}

/// The edges of `s` whose midpoint `keep` picks, by their own curve.
fn edges_where(m: &Model, s: &Shape, keep: impl Fn(Point, Point, Point) -> bool) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, e, T).unwrap() else {
                return false;
            };
            keep(
                c.point_at(a, T).unwrap(),
                c.point_at(0.5 * (a + b), T).unwrap(),
                c.point_at(b, T).unwrap(),
            )
        })
        .collect()
}

/// The slab's crease bevelled in two halves meeting tangentially at the
/// spline's middle, against the whole crease's swept section: the first
/// half's bevel runs on through the vertex where the crease carries on,
/// and the second, asked alike, is bevelled with it.
fn halves_chamfered(spec: impl Fn(&Shape) -> Chamfer, on_top: f64, on_wall: f64) {
    let (length, _) = length_and_turn();
    let mut m = Model::new();
    let (solid, creases) = split_slab(&mut m);
    assert_eq!(creases.len(), 2, "the crease's two halves");
    let top = {
        let beside: Vec<Shape> = explore_unique(&m, &solid, ShapeType::Face)
            .unwrap()
            .into_iter()
            .filter(|f| {
                matches!(
                    m.node(f)
                        .and_then(|n| n.data().as_face())
                        .and_then(|d| m.geometry().surface(d.surface)),
                    Some(ogeom_geom::SurfaceGeometry::Plane(_))
                ) && explore_unique(&m, f, ShapeType::Edge)
                    .unwrap()
                    .iter()
                    .any(|x| x.is_same(&creases[0]))
            })
            .collect();
        beside[0].clone()
    };
    let before = vol(&m, &solid);
    let specs: Vec<(Shape, Chamfer)> = creases.iter().map(|e| (e.clone(), spec(&top))).collect();
    let built = chamfer_edges_with(&mut m, &solid, &specs, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    for e in &creases {
        assert!(built.history.is_deleted(e));
    }
    let expected = before - 0.5 * on_top * on_wall * length;
    let v = vol(&m, &built.shape);
    assert!(
        (v - expected).abs() <= 1e-6 * expected,
        "{v} against {expected}"
    );
}

#[test]
fn two_halves_of_an_open_crease_bevel_as_the_whole_crease() {
    halves_chamfered(|_| Chamfer::Symmetric(0.5), 0.5, 0.5);
    halves_chamfered(
        |top| Chamfer::Distances {
            face: top.clone(),
            on_face: 0.8,
            on_other: 0.4,
        },
        0.8,
        0.4,
    );
    let angle = 35.0_f64.to_radians();
    halves_chamfered(
        |top| Chamfer::Angle {
            face: top.clone(),
            distance: 0.6,
            angle,
        },
        0.6,
        0.6 * angle.tan(),
    );
}

#[test]
fn an_open_crease_and_a_straight_edge_mitre_at_their_corner_in_either_order() {
    // The spline's crease meets the end wall's top edge at a convex corner,
    // square to it. Were the crease straight there, the two bevels of 0.5
    // would overlap in the corner by the solid a + b <= 1/2, a + c <= 1/2
    // over the positive octant: 1/24. The crease bends off its end tangent
    // by at most `kappa s² / 2` over the overlap's length `s = 1/2`, kappa
    // its curvature there, which moves the overlap's quarter-square
    // section by no more than that.
    let (length, _) = length_and_turn();
    let curve = spline();
    let kappa = (0..=64)
        .map(|i| {
            let t = 0.9 + 0.1 * f64::from(i) / 64.0;
            let d = curve.d1_at(t, T).unwrap();
            let dd = curve.derivatives_at(t, 2, T).unwrap()[2];
            d.cross(dd).magnitude() / d.magnitude().powi(3)
        })
        .fold(0.0_f64, f64::max);
    let bend = kappa * 0.25 / 2.0 * 0.25;
    let mut m = Model::new();
    let (solid, crease, _, _) = slab(&mut m);
    let before = vol(&m, &solid);
    let wall_top = edges_where(&m, &solid, |a, _, b| {
        (a.x - 10.0).abs() < 1e-9
            && (b.x - 10.0).abs() < 1e-9
            && (a.y - 2.0).abs() < 1e-9
            && (b.y - 2.0).abs() < 1e-9
    })[0]
        .clone();
    let mut volumes = Vec::new();
    for order in [[0, 1], [1, 0]] {
        let edges = [crease.clone(), wall_top.clone()];
        let specs: Vec<(Shape, Chamfer)> = order
            .iter()
            .map(|&i| (edges[i].clone(), Chamfer::Symmetric(0.5)))
            .collect();
        let built = chamfer_edges_with(&mut m, &solid, &specs, T).unwrap();
        let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        volumes.push(vol(&m, &built.shape));
    }
    assert!(
        (volumes[0] - volumes[1]).abs() <= 1e-9 * volumes[0],
        "{volumes:?}"
    );
    let straight = before - 0.125 * length - 0.125 * 4.0 + 1.0 / 24.0;
    assert!(
        (volumes[0] - straight).abs() <= bend,
        "{} against {straight} within {bend}",
        volumes[0]
    );
}
