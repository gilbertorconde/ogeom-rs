//! Rounds and bevels in one request, closed as one corner where they meet.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{Blend, Chamfer, ChamferCorner, blend_edges};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn vol(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// A 20 × 20 × 10 box and its edges at the origin along x, y and z.
fn corner(m: &mut Model) -> (Shape, [Shape; 3]) {
    let frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
    let b = ogeom_algo::make_box(m, frame, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let along = |axis: usize| {
        explore_unique(m, &b, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let (c, (s, t)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
                let (p, q) = (c.point_at(s, T).unwrap(), c.point_at(t, T).unwrap());
                let at = |p: Point| p.distance(Point::ORIGIN) < 1e-9;
                let other = if at(p) {
                    q
                } else if at(q) {
                    p
                } else {
                    return false;
                };
                let v = [other.x, other.y, other.z];
                (0..3).all(|i| (i == axis) == (v[i].abs() > 1e-9))
            })
            .unwrap()
    };
    let e = [along(0), along(1), along(2)];
    (b, e)
}

/// The band fit target the corner's patch is fitted to, in length and in
/// angle where it is held tangent.
const FIT: f64 = 2e-4;

/// The request on the box corner, `kinds[i]` for the edge along axis `i`:
/// a bevel `Symmetric(2)`, a round of 2, or none; in `order`.
fn request(e: &[Shape; 3], kinds: [usize; 3], order: [usize; 3]) -> Vec<(Shape, Blend)> {
    order
        .iter()
        .filter(|&&i| kinds[i] != 2)
        .map(|&i| {
            let blend = if kinds[i] == 0 {
                Blend::Chamfer(Chamfer::Symmetric(2.0))
            } else {
                Blend::Round(2.0)
            };
            (e[i].clone(), blend)
        })
        .collect()
}

fn surface(m: &Model, f: &Shape) -> ogeom_geom::SurfaceGeometry {
    let d = m.node(f).unwrap().data().as_face().unwrap();
    m.geometry().surface(d.surface).unwrap().clone()
}

/// Every vertex's point, sorted, as bits.
fn points(m: &Model, s: &Shape) -> Vec<[u64; 3]> {
    let mut out: Vec<[u64; 3]> = explore_unique(m, s, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| {
            let p = m.node(v).unwrap().data().as_vertex().unwrap().point;
            [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
        })
        .collect();
    out.sort_unstable();
    out
}

/// The corner `kinds` asks, checked in every order and against the same
/// request built separately; whether a trim between two rounds is held.
fn closes_as_one_patch(kinds: [usize; 3], trims_held: bool) {
    use ogeom_geom::{Surface as _, SurfaceGeometry};
    let mut first: Option<(u64, Vec<[u64; 3]>)> = None;
    for order in [
        [0usize, 1, 2],
        [2, 1, 0],
        [1, 0, 2],
        [0, 2, 1],
        [1, 2, 0],
        [2, 0, 1],
    ] {
        let mut m = Model::new();
        let (b, e) = corner(&mut m);
        let r = blend_edges(
            &mut m,
            &b,
            &request(&e, kinds, order),
            ChamferCorner::Point,
            T,
        )
        .unwrap();
        let key = (vol(&m, &r.shape).to_bits(), points(&m, &r.shape));
        match &first {
            None => first = Some(key),
            Some(held) => assert_eq!(held, &key, "order {order:?}"),
        }
    }
    let mut m = Model::new();
    let (b, e) = corner(&mut m);
    let r = blend_edges(
        &mut m,
        &b,
        &request(&e, kinds, [0, 1, 2]),
        ChamferCorner::Point,
        T,
    )
    .unwrap();
    assert!(ogeom_algo::check(&m, &r.shape, T).unwrap().is_usable());
    let all = explore_unique(&m, &r.shape, ShapeType::Face).unwrap();
    let patches: Vec<&Shape> = all
        .iter()
        .filter(|f| matches!(surface(&m, f), SurfaceGeometry::BSpline(_)))
        .collect();
    assert_eq!(patches.len(), 1);
    let patch = patches[0];
    let mut held_trims = 0;
    for edge in explore_unique(&m, patch, ShapeType::Edge).unwrap() {
        let other = all
            .iter()
            .find(|f| {
                !f.is_same(patch)
                    && explore_unique(&m, f, ShapeType::Edge)
                        .unwrap()
                        .iter()
                        .any(|g| g.is_same(&edge))
            })
            .unwrap();
        let (sp, so) = (surface(&m, patch), surface(&m, other));
        let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(&m, &edge, T)
            .unwrap()
            .unwrap();
        let (mut angle, mut gap) = (0.0_f64, 0.0_f64);
        for i in 0..=20 {
            let p = c
                .point_at((hi - lo).mul_add(f64::from(i) / 20.0, lo), T)
                .unwrap();
            let fp = ogeom_geom::project_on_surface(&sp, p, 16, T).unwrap();
            let fo = ogeom_geom::project_on_surface(&so, p, 16, T).unwrap();
            gap = gap.max(fp.distance).max(fo.distance);
            let np = sp
                .normal_at(fp.parameters.0, fp.parameters.1, T)
                .unwrap()
                .vector();
            let no = so
                .normal_at(fo.parameters.0, fo.parameters.1, T)
                .unwrap()
                .vector();
            angle = angle.max(np.cross(no).magnitude().atan2(np.dot(no).abs()));
        }
        assert!(gap <= FIT, "the patch stands {gap} off a face it meets");
        match so {
            // A round's band: tangent.
            SurfaceGeometry::Cylinder(_) => assert!(angle <= FIT, "{angle} off a round"),
            SurfaceGeometry::Plane(_) if angle <= FIT => held_trims += 1,
            _ => {}
        }
    }
    assert_eq!(held_trims > 0, trims_held);

    // Away from the corner: the same request built separately, the bevels
    // and then the rounds, has every face on the very same surface.
    let mut s = Model::new();
    let (sb, se) = corner(&mut s);
    let bevels: Vec<(Shape, Chamfer)> = (0..3)
        .filter(|&i| kinds[i] == 0)
        .map(|i| (se[i].clone(), Chamfer::Symmetric(2.0)))
        .collect();
    let beveled =
        ogeom_fillet::chamfer_edges_cornered(&mut s, &sb, &bevels, ChamferCorner::Point, T)
            .unwrap();
    let rounds: Vec<(Shape, f64)> = (0..3)
        .filter(|&i| kinds[i] == 1)
        .map(|i| (beveled.history.trace(&se[i])[0].clone(), 2.0))
        .collect();
    let separate = ogeom_fillet::fillet_edges_with(&mut s, &beveled.shape, &rounds, T).unwrap();
    let mut far = 0;
    for f in explore_unique(&s, &separate.shape, ShapeType::Face).unwrap() {
        let reaches_out = explore_unique(&s, &f, ShapeType::Vertex)
            .unwrap()
            .iter()
            .any(|v| {
                s.node(v)
                    .unwrap()
                    .data()
                    .as_vertex()
                    .unwrap()
                    .point
                    .distance(Point::ORIGIN)
                    > 8.0
            });
        if !reaches_out {
            continue;
        }
        far += 1;
        assert!(
            all.iter().any(|g| surface(&m, g) == surface(&s, &f)),
            "a face away from the corner moved"
        );
    }
    assert!(far >= 8, "{far}");
}

#[test]
fn a_bevel_and_two_rounds_meeting_close_as_one_patch_whatever_the_order() {
    closes_as_one_patch([0, 1, 1], true);
}

#[test]
fn two_bevels_and_a_round_meeting_close_as_one_patch_whatever_the_order() {
    closes_as_one_patch([0, 0, 1], false);
}

/// A round of 2 meeting a bevel `Symmetric(2)` at the box's corner, the
/// third edge left sharp, each edge taking each part in turn and the
/// request asked both ways round. The opening crosses the sharp edge at a
/// point, where a patch held across both of its faces turns its normal
/// from one face's to the other's, and one patch closes the rest: the two
/// orders build the same solid to the bit, valid, its patches tangent to
/// the round's band and, at the sharp edge, to both its faces, and meeting
/// every face they border in position.
#[test]
fn a_round_meeting_a_bevel_beside_a_sharp_edge_closes_whatever_the_order() {
    use ogeom_geom::SurfaceGeometry;
    for kinds in [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        let mut first: Option<(u64, Vec<[u64; 3]>)> = None;
        for order in [[0usize, 1, 2], [2, 1, 0]] {
            let mut m = Model::new();
            let (b, e) = corner(&mut m);
            let r = blend_edges(
                &mut m,
                &b,
                &request(&e, kinds, order),
                ChamferCorner::Point,
                T,
            )
            .unwrap_or_else(|why| panic!("{kinds:?} {order:?}: {why}"));
            let diagnosis = ogeom_algo::check(&m, &r.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{kinds:?}: {:?}", diagnosis.problems);
            let key = (vol(&m, &r.shape).to_bits(), points(&m, &r.shape));
            match &first {
                None => first = Some(key),
                Some(held) => assert_eq!(held, &key, "{kinds:?}: order {order:?}"),
            }
            if order != [0, 1, 2] {
                continue;
            }
            let all = explore_unique(&m, &r.shape, ShapeType::Face).unwrap();
            let patches: Vec<&Shape> = all
                .iter()
                .filter(|f| matches!(surface(&m, f), SurfaceGeometry::BSpline(_)))
                .collect();
            assert_eq!(
                patches.len(),
                2,
                "{kinds:?}: the patch at the sharp edge and the rest"
            );
            let sharp = (0..3).find(|&i| kinds[i] == 2).unwrap();
            // The faces of the sharp edge: the planes through the corner
            // square to the other two axes.
            let on_sharp_face = |n: ogeom_math::Vector| {
                (0..3).any(|i| i != sharp && ([n.x, n.y, n.z][i].abs() - 1.0).abs() < 1e-9)
            };
            let mut worst_round = 0.0_f64;
            let mut at_q_faces = 0;
            for patch in &patches {
                for c in ogeom_fillet::analyse_blend(&m, &r.shape, patch, 17, T).unwrap() {
                    assert!(
                        c.gap <= FIT,
                        "{kinds:?}: a patch stands {} off a face",
                        c.gap
                    );
                    match surface(&m, &c.neighbour) {
                        SurfaceGeometry::Cylinder(_) => {
                            worst_round = worst_round.max(c.tangency_error);
                        }
                        SurfaceGeometry::Plane(_) => {
                            let (_, n) = ogeom_algo::face_normal(&m, &c.neighbour, T).unwrap();
                            if on_sharp_face(n) && c.tangency_error <= FIT {
                                at_q_faces += 1;
                            }
                        }
                        _ => {}
                    }
                }
            }
            eprintln!("{kinds:?}: the patches meet the round within {worst_round} rad");
            assert!(
                worst_round <= FIT,
                "{kinds:?}: the patches meet the round {worst_round} rad off tangent"
            );
            assert!(
                at_q_faces >= 2,
                "{kinds:?}: the patches are tangent to {at_q_faces} trims on the sharp edge's faces"
            );
        }
    }
}

#[test]
fn a_request_of_one_kind_is_that_kind_s_own_call() {
    let mut m = Model::new();
    let (b, e) = corner(&mut m);
    let mut plain = m.clone();
    let rounds: Vec<(Shape, f64)> = e.iter().map(|x| (x.clone(), 2.0)).collect();
    let own = ogeom_fillet::fillet_edges_with(&mut plain, &b, &rounds, T).unwrap();
    let blended = blend_edges(
        &mut m,
        &b,
        &request(&e, [1, 1, 1], [0, 1, 2]),
        ChamferCorner::Point,
        T,
    )
    .unwrap();
    assert_eq!(
        vol(&plain, &own.shape).to_bits(),
        vol(&m, &blended.shape).to_bits()
    );
    assert_eq!(points(&plain, &own.shape), points(&m, &blended.shape));

    let mut m = Model::new();
    let (b, e) = corner(&mut m);
    let mut plain = m.clone();
    let bevels: Vec<(Shape, Chamfer)> = e
        .iter()
        .map(|x| (x.clone(), Chamfer::Symmetric(2.0)))
        .collect();
    let own =
        ogeom_fillet::chamfer_edges_cornered(&mut plain, &b, &bevels, ChamferCorner::Facet, T)
            .unwrap();
    let blended = blend_edges(
        &mut m,
        &b,
        &request(&e, [0, 0, 0], [0, 1, 2]),
        ChamferCorner::Facet,
        T,
    )
    .unwrap();
    assert_eq!(
        vol(&plain, &own.shape).to_bits(),
        vol(&m, &blended.shape).to_bits()
    );
    assert_eq!(points(&plain, &own.shape), points(&m, &blended.shape));
}

#[test]
fn a_round_and_a_bevel_that_do_not_meet_are_their_own_calls_in_turn() {
    // The bevel along x at the origin and a round on the top edge along x
    // at the far side: no vertex shared, so the bevel and then the round,
    // each as its own call builds it.
    let edge_at = |m: &Model, b: &Shape, p: Point| {
        explore_unique(m, b, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let (c, (s, t)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
                c.point_at(f64::midpoint(s, t), T).unwrap().distance(p) < 1e-9
            })
            .unwrap()
    };
    let mut m = Model::new();
    let (b, e) = corner(&mut m);
    let far = edge_at(&m, &b, Point::new(10.0, 20.0, 10.0));
    let mut plain = m.clone();
    let beveled = ogeom_fillet::chamfer_edges_with(
        &mut plain,
        &b,
        &[(e[0].clone(), Chamfer::Symmetric(2.0))],
        T,
    )
    .unwrap();
    let traced = beveled.history.trace(&far)[0].clone();
    let own =
        ogeom_fillet::fillet_edges_with(&mut plain, &beveled.shape, &[(traced, 2.0)], T).unwrap();
    let blended = blend_edges(
        &mut m,
        &b,
        &[
            (far.clone(), Blend::Round(2.0)),
            (e[0].clone(), Blend::Chamfer(Chamfer::Symmetric(2.0))),
        ],
        ChamferCorner::Point,
        T,
    )
    .unwrap();
    assert_eq!(
        vol(&plain, &own.shape).to_bits(),
        vol(&m, &blended.shape).to_bits()
    );
    assert_eq!(points(&plain, &own.shape), points(&m, &blended.shape));
}
