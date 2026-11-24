//! A round of varying radius through the blend engine: a box's edge, the
//! ball's radius following a law along it, flush in the end faces.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::FRAC_PI_4;

use ogeom_core::Tolerances;
use ogeom_fillet::RadiusLaw;
use ogeom_math::{Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A 20 by 10 by 6 box's top front edge (along `x` at `y = 0`, `z = 6`)
/// rounded from 1 at its first vertex to 2 at its last. Through the engine
/// with no boolean, valid, the strip generated; the box less the section
/// `(1 - pi/4) r^2` integrated along the edge, `20 (1 - pi/4) 7/3`; and the
/// contact on the top 1 in at the first vertex, 2 at the last.
#[test]
fn a_box_edge_rounds_along_a_linear_law() {
    let mut m = Model::new();
    let s = ogeom_algo::make_box(&mut m, Frame::WORLD, (20.0, 10.0, 6.0), T)
        .unwrap()
        .shape;
    let point = |m: &Model, v: &Shape| m.node(v).unwrap().data().as_vertex().unwrap().point;
    let edge = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
            let (a, b) = (point(&m, &a), point(&m, &b));
            a.y.abs() < 1e-9 && b.y.abs() < 1e-9 && a.z > 5.9 && b.z > 5.9
        })
        .unwrap();
    let (first, last) = ogeom_algo::edge_vertices(&m, &edge).unwrap().unwrap();
    let (first, last) = (point(&m, &first), point(&m, &last));
    let law = RadiusLaw::Linear {
        start: 1.0,
        end: 2.0,
    };
    let (built, work) = ogeom_bool::work::counted(|| {
        ogeom_fillet::fillet_edges_variable(&mut m, &s, std::slice::from_ref(&edge), &law, T)
            .unwrap()
    });
    assert_eq!(work.face_pairs, 0, "built with no boolean");
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    assert!(!built.history.generated(&edge).is_empty(), "the strip");
    let got = ogeom_algo::volume_properties(&m, &built.shape, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass;
    let want = 1200.0 - 20.0 * (1.0 - FRAC_PI_4) * 7.0 / 3.0;
    assert!((got - want).abs() < 1e-6, "{got} against {want}");
    let vertices: Vec<Point> = explore_unique(&m, &built.shape, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| point(&m, v))
        .collect();
    for (end, r) in [(first, 1.0), (last, 2.0)] {
        let contact = Point::new(end.x, r, 6.0);
        assert!(
            vertices.iter().any(|p| p.distance(contact) < 1e-9),
            "the top's contact {r} in at {end:?}"
        );
    }
}
