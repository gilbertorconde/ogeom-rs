//! Booleans against a solid whose rounded corner closes on a patch with a
//! point side lying on a plane of the other solid.
#![allow(clippy::unwrap_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::Point;
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A prism of the polygon `profile` in (y, z), swept from x = 0 to
/// `length`.
fn prism(m: &mut Model, profile: &[(f64, f64)], length: f64) -> Shape {
    let n = profile.len();
    let mut points: Vec<Point> = profile
        .iter()
        .map(|&(y, z)| Point::new(0.0, y, z))
        .collect();
    points.extend(profile.iter().map(|&(y, z)| Point::new(length, y, z)));
    let mut rings: Vec<Vec<usize>> = vec![(0..n).rev().collect(), (n..2 * n).collect()];
    for i in 0..n {
        let j = (i + 1) % n;
        rings.push(vec![i, j, n + j, n + i]);
    }
    ogeom_algo::make_faceted_solid(m, &points, &rings, &[], T)
        .unwrap()
        .0
        .shape
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(
        m,
        s,
        ogeom_mesh::Deflection {
            chord: 1e-5,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

/// An L prism (the profile (0, 0), (10, 0), (10, 10), (8, 10), (8, 2),
/// (0, 2) swept 20 along x), its three edges at (0, 8, 2) rounded 1: the
/// wall's concave foot and two convex edges, their corner closed on one
/// patch whose point side lies at (1, 7, 2) on the floor plane z = 2, the
/// patch nearly tangent to the plane there. The section of the L's floor
/// with the patch runs into that point, and the marcher stopped 4e-4 short
/// of it, so the L less the rounded L and their common part refused ("the
/// kept pieces did not close into a shell"). Each boolean now builds, valid,
/// and the pieces add up: the cut and the common part make the L in volume,
/// and the rounded L less the L with the common part make the rounded L.
#[test]
fn booleans_against_a_corner_patch_with_a_point_side_on_a_plane_close() {
    let mut m = Model::new();
    let l = prism(
        &mut m,
        &[
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (8.0, 10.0),
            (8.0, 2.0),
            (0.0, 2.0),
        ],
        20.0,
    );
    let at = Point::new(0.0, 8.0, 2.0);
    let edges: Vec<Shape> = explore_unique(&m, &l, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
            [a, b].iter().any(|v| {
                m.node(v)
                    .unwrap()
                    .data()
                    .as_vertex()
                    .unwrap()
                    .point
                    .distance(at)
                    < 1e-9
            })
        })
        .collect();
    assert_eq!(edges.len(), 3, "the three edges at the corner");
    let rounded = ogeom_fillet::fillet_edges(&mut m, &l, &edges, 1.0, T)
        .unwrap()
        .shape;
    let less = ogeom_bool::cut(&mut m, &l, &rounded, T).unwrap().shape;
    let both = ogeom_bool::common(&mut m, &l, &rounded, T).unwrap().shape;
    let back = ogeom_bool::cut(&mut m, &rounded, &l, T).unwrap().shape;
    for s in [&less, &both, &back] {
        let found = ogeom_algo::check(&m, s, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
    }
    let (vl, vr) = (volume(&m, &l), volume(&m, &rounded));
    let (vless, vboth, vback) = (volume(&m, &less), volume(&m, &both), volume(&m, &back));
    assert!(vless > 0.0 && vback > 0.0, "{vless} {vback}");
    assert!(
        (vless + vboth - vl).abs() <= 1e-6 * vl,
        "{vless} + {vboth} against {vl}"
    );
    assert!(
        (vback + vboth - vr).abs() <= 1e-6 * vr,
        "{vback} + {vboth} against {vr}"
    );
}
