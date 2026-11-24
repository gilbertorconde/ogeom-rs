//! Rounds of an elliptical or conic section through the blend engine: a
//! box's edge, each face its own distance, the strip the section swept
//! along the edge, flush in the end faces.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::FRAC_PI_4;

use ogeom_core::Tolerances;
use ogeom_fillet::RoundProfile;
use ogeom_math::{Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// A 20 by 10 by 6 box, its top front edge (along `x` at `y = 0`, `z = 6`)
/// and its top and front faces.
fn boxed(m: &mut Model) -> (Shape, Shape, Shape, Shape) {
    let s = ogeom_algo::make_box(m, Frame::WORLD, (20.0, 10.0, 6.0), T)
        .unwrap()
        .shape;
    let edge = explore_unique(m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (a, b) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
            let q = |v: &Shape| m.node(v).unwrap().data().as_vertex().unwrap().point;
            let (a, b) = (q(&a), q(&b));
            a.y.abs() < 1e-9
                && b.y.abs() < 1e-9
                && (a.z - 6.0).abs() < 1e-9
                && (b.z - 6.0).abs() < 1e-9
        })
        .unwrap();
    // The face whose outward normal is `n`.
    let facing = |m: &Model, n: Vector| {
        explore_unique(m, &s, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| ogeom_algo::face_normal(m, f, T).unwrap().1.dot(n) > 0.99)
            .unwrap()
    };
    let top = facing(m, Vector::Z);
    let front = facing(m, -Vector::Y);
    (s, edge, top, front)
}

/// The edge rounded meeting the top 2 in and the front 1 down, an ellipse
/// and a parabola (`rho = 1/2`), named either way round. Through the
/// engine with no boolean, valid, the box less the section's spandrel along
/// the 20: the ellipse's `(1 - pi/4) d1 d2`, the parabola's a third of the
/// triangle `d1 d2 / 2` (a parabola's arc cuts two thirds of its chord's
/// triangle off).
#[test]
fn a_box_edge_s_elliptical_and_conic_rounds_take_their_exact_spandrels() {
    for (profile, spandrel) in [
        (RoundProfile::Elliptical, (1.0 - FRAC_PI_4) * 2.0),
        (RoundProfile::Conic(0.5), 2.0 / 6.0),
    ] {
        let mut volumes = Vec::new();
        for swap in [false, true] {
            let mut m = Model::new();
            let (s, edge, top, front) = boxed(&mut m);
            let (face, on_face, on_other) = if swap {
                (front, 1.0, 2.0)
            } else {
                (top, 2.0, 1.0)
            };
            let (built, work) = ogeom_bool::work::counted(|| {
                ogeom_fillet::fillet_edge_distances(
                    &mut m, &s, &edge, &face, on_face, on_other, profile, T,
                )
                .unwrap()
            });
            assert_eq!(work.face_pairs, 0, "{profile:?}: built with no boolean");
            let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
            assert!(found.is_valid(), "{:?}", found.problems);
            assert!(!built.history.generated(&edge).is_empty(), "the strip");
            let got = volume(&m, &built.shape);
            let want = 1200.0 - spandrel * 20.0;
            assert!(
                (got - want).abs() < 1e-6,
                "{profile:?}: {got} against {want}"
            );
            volumes.push(got);
        }
        assert!((volumes[0] - volumes[1]).abs() < 1e-9, "either face named");
    }
}

/// The box's corner (20, 0, 6) through `blend_edges`: the top front edge
/// elliptical (2 in on the top, 1 down the front), the two others round at
/// 1, closing on one patch. Through the engine with no boolean, valid,
/// every edge generating its strip, the same solid in either order. Each
/// strip is exact clear of the corner's box (2 each way, its setbacks'
/// reach): the volume lies between the box less those strips and that
/// less the corner box too.
#[test]
fn an_elliptical_round_meeting_two_rounds_at_a_corner_closes_on_one_patch() {
    let mut volumes = Vec::new();
    for flip in [false, true] {
        let mut m = Model::new();
        let (s, edge, top, _) = boxed(&mut m);
        let at = |e: &Shape, m: &Model| {
            let (a, b) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
            let q = |v: &Shape| m.node(v).unwrap().data().as_vertex().unwrap().point;
            (q(&a), q(&b))
        };
        let others: Vec<Shape> = explore_unique(&m, &s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                let (a, b) = at(e, &m);
                let corner = Point::new(20.0, 0.0, 6.0);
                !e.is_same(&edge) && (a.distance(corner) < 1e-9 || b.distance(corner) < 1e-9)
            })
            .collect();
        assert_eq!(others.len(), 2);
        let mut asked: Vec<(Shape, ogeom_fillet::Blend)> = vec![(
            edge.clone(),
            ogeom_fillet::Blend::RoundDistances {
                face: top.clone(),
                on_face: 2.0,
                on_other: 1.0,
                profile: RoundProfile::Elliptical,
            },
        )];
        asked.extend(
            others
                .iter()
                .map(|e| (e.clone(), ogeom_fillet::Blend::Round(1.0))),
        );
        if flip {
            asked.reverse();
        }
        let (built, work) = ogeom_bool::work::counted(|| {
            ogeom_fillet::blend_edges(&mut m, &s, &asked, ogeom_fillet::ChamferCorner::Facet, T)
                .unwrap()
        });
        assert_eq!(work.face_pairs, 0, "built with no boolean");
        let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        for (e, _) in &asked {
            assert!(!built.history.generated(e).is_empty(), "each strip");
        }
        volumes.push(volume(&m, &built.shape));
    }
    let round = 1.0 - FRAC_PI_4;
    let strips = round * 2.0 * 18.0 + round * 8.0 + round * 4.0;
    let (hi, lo) = (1200.0 - strips, 1200.0 - strips - 8.0);
    for v in &volumes {
        assert!(*v < hi && *v > lo, "{v} outside ({lo}, {hi})");
    }
    assert!((volumes[0] - volumes[1]).abs() < 1e-9, "either order");
}
