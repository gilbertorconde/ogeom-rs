//! Chamfers on a cone's rim, in all three spellings, against Pappus.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{Chamfer, chamfer_edges_with};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const TAU: f64 = core::f64::consts::TAU;

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

/// A frustum of base radius 10 and top radius 6, 10 high: its top rim
/// and the cap above it.
fn frustum(m: &mut Model) -> (Shape, Shape, Shape) {
    let frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
    let cone = ogeom_algo::make_cone(m, frame, 10.0, 6.0, 10.0, T)
        .unwrap()
        .shape;
    let rim = explore_unique(m, &cone, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (c, (a, _)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            let p = c.point_at(a, T).unwrap();
            (p.z - 10.0).abs() < 1e-9 && (p.x.hypot(p.y) - 6.0).abs() < 1e-9
        })
        .expect("the top rim");
    let cap = explore_unique(m, &cone, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let b = ogeom_algo::shape_bounds(m, f, T).unwrap();
            (b.low().unwrap().z - 10.0).abs() < 1e-3 && (b.high().unwrap().z - 10.0).abs() < 1e-3
        })
        .expect("the top cap");
    (cone, rim, cap)
}

/// The volume the triangle `a b c` of the rim's section sweeps about the
/// axis, by Pappus: its area times the turn of its centroid.
fn pappus(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    let area = 0.5 * ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs();
    area * TAU * (a.0 + b.0 + c.0) / 3.0
}

/// The section's wall leg from the rim, `s` along the slant: down and out.
fn on_wall(s: f64) -> (f64, f64) {
    let k = 116.0_f64.sqrt();
    (6.0 + 4.0 * s / k, 10.0 - 10.0 * s / k)
}

fn chamfered(spec: impl Fn(&Shape) -> Chamfer) -> f64 {
    let mut m = Model::new();
    let (cone, rim, cap) = frustum(&mut m);
    let built = chamfer_edges_with(&mut m, &cone, &[(rim.clone(), spec(&cap))], T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!(built.history.is_deleted(&rim));
    vol(&m, &built.shape)
}

fn whole() -> f64 {
    core::f64::consts::PI * 10.0 / 3.0 * (100.0 + 60.0 + 36.0)
}

fn close(v: f64, expected: f64) {
    assert!(
        (v - expected).abs() <= 1e-6 * expected,
        "{v} against {expected}"
    );
}

#[test]
fn a_cone_rim_chamfers_symmetric_to_the_revolved_triangle() {
    // 2052.507 less 16.891: 2035.62.
    let expected = whole() - pappus((6.0, 10.0), (5.0, 10.0), on_wall(1.0));
    close(chamfered(|_| Chamfer::Symmetric(1.0)), expected);
}

#[test]
fn a_cone_rim_chamfers_by_two_distances() {
    let expected = whole() - pappus((6.0, 10.0), (4.5, 10.0), on_wall(0.8));
    close(
        chamfered(|cap| Chamfer::Distances {
            face: cap.clone(),
            on_face: 1.5,
            on_other: 0.8,
        }),
        expected,
    );
}

#[test]
fn a_cone_rim_chamfers_by_distance_and_angle() {
    // From the cap's contact at 1.2 in, the bevel leaves the cap at 40
    // degrees down into the material; where that line meets the wall's
    // slant is the other contact.
    let (d, angle) = (1.2, 40.0_f64.to_radians());
    let c = (6.0 - d, 10.0);
    let (dir_x, dir_z) = (angle.cos(), -angle.sin());
    // Solve c + t dir = (6, 10) + s (4, -10)/sqrt(116).
    let k = 116.0_f64.sqrt();
    let (wx, wz) = (4.0 / k, -10.0 / k);
    let det = dir_x * (-wz) - (-wx) * dir_z;
    let (rx, rz) = (6.0 - c.0, 10.0 - c.1);
    let s = (dir_x * rz - dir_z * rx) / det;
    let expected = whole() - pappus((6.0, 10.0), c, on_wall(s));
    close(
        chamfered(|cap| Chamfer::Angle {
            face: cap.clone(),
            distance: d,
            angle,
        }),
        expected,
    );
}
