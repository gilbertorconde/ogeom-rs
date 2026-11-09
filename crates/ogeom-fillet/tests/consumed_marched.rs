//! A round of a bore's rim on a drum taking a narrow bevel beside it off
//! whole: the marched counterpart of a rim's round over a chamfer.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{chamfer_edges, fillet_edges};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
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
