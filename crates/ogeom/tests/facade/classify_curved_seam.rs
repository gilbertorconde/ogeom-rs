//! A face classified about an edge that bends more in its chart than in
//! space: the parabola where a drafted drum (a cone) meets a wall drafted
//! as the drum is.
//!
//! The face's trim is a polyline in the cone's chart. Drawn from the seam's
//! points in space alone, it stood up to 6.9e-4 off the seam at a chord of
//! 5e-4, read in space through the chart's scales, and a point on the seam
//! read as off the face.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::Containment;
use ogeom::core::Tolerances;
use ogeom::geom::{Curve3d as _, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, explore, explore_unique};

const T: Tolerances = Tolerances::millimetres();
/// The drum's radius at its top, two high.
const BIG: f64 = 2.5;
/// Where the wall cuts the top, `y = C`.
const C: f64 = 0.5;

/// The wall's unit normal, out of the material, for a draft of `k`.
fn wall_normal(k: f64) -> Vector {
    Vector::new(0.0, 1.0, k) * (1.0 / k.hypot(1.0))
}

/// A drum drafted by `k` over its height (its foot standing out by `2k`),
/// cut by the wall through `y = C` at the top drafted alike: the seam the
/// two leave is a parabola, the wall lying parallel to one of the cone's
/// rulings.
fn drafted_d(model: &mut Model, k: f64) -> Shape {
    let disc = ogeom::algo::make_cone(model, Frame::WORLD, 2.0_f64.mul_add(k, BIG), BIG, 2.0, T)
        .unwrap()
        .shape;
    let normal = wall_normal(k);
    let across = normal.cross(Vector::X);
    let frame = Frame::new(
        Point::new(0.0, C, 2.0) - Vector::X * 20.0 - across * 20.0,
        Direction::new(normal, T).unwrap(),
        Direction::X,
        T,
    )
    .unwrap();
    let past = ogeom::algo::make_box(model, frame, (40.0, 40.0, 40.0), T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(model, &disc, &past, T).unwrap().shape;
    explore(model, &cut, Filter::OfType(ShapeType::Solid)).unwrap()[0].clone()
}

/// Points on the cone face's parabolic seam, at the inner nine tenths of
/// its range, each with the unit direction along the cone across the seam
/// out of the face.
fn seam_points(model: &Model, cone: &Shape, k: f64) -> Vec<(Point, Vector)> {
    let mut out = Vec::new();
    for edge in explore_unique(model, cone, ShapeType::Edge).unwrap() {
        let NodeData::Edge(data) = model.node(&edge).unwrap().data() else {
            unreachable!()
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let curve = model.geometry().curve(*curve).unwrap();
        if !matches!(curve, ogeom::geom::Curve::Parabola(_)) {
            continue;
        }
        for i in 1..=19 {
            let t = (range.1 - range.0).mul_add(f64::from(i) / 20.0, range.0);
            let p = curve.point_at(t, T).unwrap();
            // The cone `|xy| = BIG + k (2 - z)`: its normal is the gradient
            // `(x / r, y / r, k)`.
            let r = p.x.hypot(p.y);
            let normal = Vector::new(p.x / r, p.y / r, k);
            let along = curve.tangent_at(t, T).unwrap().vector();
            let mut across = normal.cross(along);
            across = across * (1.0 / across.magnitude());
            if across.dot(wall_normal(k)) < 0.0 {
                across = -across;
            }
            out.push((p, across));
        }
    }
    out
}

/// A point on the seam lies on the cone face, and points along the cone
/// twice the chord in from it and out past it lie in it and off it: the
/// trim keeps to the seam within the chord in the cone's chart. Before,
/// at the fillet run-out's chord of 5e-4 and finer, about half the seam's
/// points read as off the face, and a point twice the chord in read as on
/// its trim.
#[test]
fn a_point_on_a_parabolic_seam_lies_on_the_cone() {
    for k in [0.1_f64, -0.1] {
        let mut model = Model::new();
        let solid = drafted_d(&mut model, k);
        let cone = explore_unique(&model, &solid, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                let NodeData::Face(data) = model.node(f).unwrap().data() else {
                    return false;
                };
                matches!(
                    model.geometry().surface(data.surface),
                    Some(SurfaceGeometry::Cone(_))
                )
            })
            .unwrap();
        let seam = seam_points(&model, &cone, k);
        assert_eq!(seam.len(), 38, "draft {k}: the seam's two parabolas");
        for chord in [1e-3, 5e-4, 1e-4] {
            let deflection = Deflection::with_chord(chord).unwrap();
            let classify =
                |p: Point| ogeom::algo::classify_on_face(&model, &cone, p, deflection, T).unwrap();
            for &(p, across) in &seam {
                assert_eq!(
                    classify(p),
                    Containment::On,
                    "draft {k}, chord {chord}: {p:?}"
                );
                // A step along the tangent plane, back onto the cone: the
                // cone curves by under 1e-7 over it.
                let step = |s: f64| {
                    let q = p + across * s;
                    let r = q.x.hypot(q.y);
                    let want = k.mul_add(2.0 - q.z, BIG);
                    Point::new(q.x * want / r, q.y * want / r, q.z)
                };
                assert_eq!(
                    classify(step(-2.0 * chord)),
                    Containment::In,
                    "draft {k}, chord {chord}: in from {p:?}"
                );
                assert_eq!(
                    classify(step(2.0 * chord)),
                    Containment::Out,
                    "draft {k}, chord {chord}: out from {p:?}"
                );
            }
        }
    }
}
