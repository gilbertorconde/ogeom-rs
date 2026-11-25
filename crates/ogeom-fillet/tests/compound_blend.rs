//! A blend on a compound of solids: each solid blended alone, the solids
//! with no edge asked handed back as they are, the compound kept.
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom_core::Tolerances;
use ogeom_fillet::Chamfer;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom_algo::volume_properties(model, shape, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// Two boxes 80 by 30 by 6 side by side in one compound, touching along a
/// face, their tops on one plane: `near` at the origin, `far` beside it
/// at y = 30.
fn touching_pair(model: &mut Model) -> (Shape, Shape, Shape) {
    let near = ogeom_algo::make_box(model, Frame::WORLD, (80.0, 30.0, 6.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(0.0, 30.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let far = ogeom_algo::make_box(model, at, (80.0, 30.0, 6.0), T)
        .unwrap()
        .shape;
    let both = ogeom_algo::make_compound(model, &[near.clone(), far.clone()])
        .unwrap()
        .shape;
    (near, far, both)
}

/// The result holds two solids, `far` the very solid given (its faces
/// untouched, 14400 of volume), and the other `near` blended to `want`;
/// both valid.
fn assert_far_kept(model: &Model, built: &Shape, far: &Shape, want: f64) {
    assert_eq!(model.kind_of(built).unwrap(), ShapeType::Compound);
    let parts = model.children_of(built).unwrap();
    assert_eq!(parts.len(), 2, "two solids back");
    assert_eq!(
        explore_unique(model, built, ShapeType::Solid)
            .unwrap()
            .len(),
        2
    );
    let kept: Vec<&Shape> = parts.iter().filter(|p| p.is_same(far)).collect();
    assert_eq!(kept.len(), 1, "far handed back as it was");
    let faces = |s: &Shape| explore_unique(model, s, ShapeType::Face).unwrap();
    let (was, now) = (faces(far), faces(kept[0]));
    assert_eq!(now.len(), 6, "far's faces");
    assert!(was.iter().all(|f| now.iter().any(|g| g.is_same(f))));
    let got = volume(model, kept[0]);
    assert!((got - 14_400.0).abs() < 1e-6, "far {got}");
    let blended = parts.iter().find(|p| !p.is_same(far)).unwrap();
    let got = volume(model, blended);
    assert!(
        (got - want).abs() <= 1e-8 * want,
        "near {got} against {want}"
    );
    for p in &parts {
        let found = ogeom_algo::check(model, p, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
    }
}

/// Every edge of `near` chamfered 1: `near` loses `2d²(a + b + c) −
/// 16d³/3` (the facet at each corner), `far` is untouched.
#[test]
fn chamfering_one_solid_of_a_compound_leaves_a_touching_one_apart() {
    let mut m = Model::new();
    let (near, far, both) = touching_pair(&mut m);
    let specs: Vec<(Shape, Chamfer)> = explore_unique(&m, &near, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .map(|e| (e, Chamfer::Symmetric(1.0)))
        .collect();
    let built = ogeom_fillet::chamfer_edges_with(&mut m, &both, &specs, T).unwrap();
    let d = 1.0_f64;
    let want = 14_400.0 - (2.0 * d * d * (80.0 + 30.0 + 6.0) - 16.0 * d.powi(3) / 3.0);
    assert_far_kept(&m, &built.shape, &far, want);
}

/// Every edge of `near` rounded 1: `near` is its inner box swept by the
/// ball, `p q w + 2r(pq + qw + wp) + πr²(p + q + w) + 4/3 πr³` with the
/// primes the sides less 2r; `far` is untouched.
#[test]
fn rounding_one_solid_of_a_compound_leaves_a_touching_one_apart() {
    let mut m = Model::new();
    let (near, far, both) = touching_pair(&mut m);
    let edges = explore_unique(&m, &near, ShapeType::Edge).unwrap();
    let r = 1.0_f64;
    let built = ogeom_fillet::fillet_edges(&mut m, &both, &edges, r, T).unwrap();
    let (p, q, w) = (80.0 - 2.0 * r, 30.0 - 2.0 * r, 6.0 - 2.0 * r);
    let want = p * q * w
        + 2.0 * r * (p * q + q * w + w * p)
        + PI * r * r * (p + q + w)
        + 4.0 / 3.0 * PI * r.powi(3);
    assert_far_kept(&m, &built.shape, &far, want);
    let asked: Vec<(Shape, f64)> = edges.iter().map(|e| (e.clone(), r)).collect();
    let mut m2 = m.clone();
    let built = ogeom_fillet::fillet_edges_with(&mut m2, &both, &asked, T).unwrap();
    assert_far_kept(&m2, &built.shape, &far, want);
}

/// An edge held by two solids of one compound (here one box standing in
/// it twice) is refused by name, for a chamfer and a round alike: each
/// body is blended on its own, so no edge belongs to two.
#[test]
fn an_edge_two_solids_of_a_compound_share_is_refused_by_name() {
    let mut m = Model::new();
    let a = ogeom_algo::make_box(&mut m, Frame::WORLD, (10.0, 8.0, 6.0), T)
        .unwrap()
        .shape;
    let twice = ogeom_algo::make_compound(&mut m, &[a.clone(), a.clone()])
        .unwrap()
        .shape;
    let edge = explore_unique(&m, &a, ShapeType::Edge).unwrap()[0].clone();
    let refused = ogeom_fillet::chamfer_edges_with(
        &mut m,
        &twice,
        &[(edge.clone(), Chamfer::Symmetric(1.0))],
        T,
    )
    .unwrap_err();
    assert!(refused.to_string().contains("two solids"), "{refused}");
    let refused = ogeom_fillet::fillet_edges(&mut m, &twice, &[edge], 1.0, T).unwrap_err();
    assert!(refused.to_string().contains("two solids"), "{refused}");
}

/// The top edge of `near` along x on the side it shares with `far`
/// (y = 30, z = 6), and `near`'s top face.
fn shared_side_edge_and_top(model: &Model, near: &Shape) -> (Shape, Shape) {
    let bounds = |s: &Shape| ogeom_algo::shape_bounds(model, s, T).unwrap();
    let edge = explore_unique(model, near, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let b = bounds(e);
            let (l, h) = (b.low().unwrap(), b.high().unwrap());
            (l.y - 30.0).abs() < 1e-6
                && (h.y - 30.0).abs() < 1e-6
                && (l.z - 6.0).abs() < 1e-6
                && h.x - l.x > 50.0
        })
        .unwrap();
    let top = explore_unique(model, near, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| (bounds(f).low().unwrap().z - 6.0).abs() < 1e-6)
        .unwrap();
    (edge, top)
}

/// A blend entry asked on its second shape, given `near`, `near`'s edge
/// on the shared side and its top face.
type Entry =
    fn(&mut Model, &Shape, &Shape, &Shape, &Shape) -> ogeom_core::OgeomResult<ogeom_algo::Built>;

/// The other blend entries, each asked on the touching pair for `near`
/// alone, build `near` as they build it standing by itself and hand `far`
/// back as it was; before, each read the compound as one shape and refused
/// it. Every edge of `near` mitred at its corners loses `2d²(a + b + c) −
/// 6d³`; rolling corners round it to the inner box swept by the ball; the
/// rest are held to the same call on `near` alone.
#[test]
fn every_blend_entry_works_on_one_solid_of_a_compound_alone() {
    let mitred = 14_400.0 - (2.0 * (80.0 + 30.0 + 6.0) - 6.0);
    let (p, q, w) = (78.0, 28.0, 4.0);
    let swept = p * q * w + 2.0 * (p * q + q * w + w * p) + PI * (p + q + w) + 4.0 / 3.0 * PI;
    let entries: Vec<(&str, Entry, Option<f64>)> = vec![
        (
            "chamfer_edges_cornered, point",
            |m, near, target, _, _| {
                let specs: Vec<(Shape, Chamfer)> = explore_unique(m, near, ShapeType::Edge)
                    .unwrap()
                    .into_iter()
                    .map(|e| (e, Chamfer::Symmetric(1.0)))
                    .collect();
                ogeom_fillet::chamfer_edges_cornered(
                    m,
                    target,
                    &specs,
                    ogeom_fillet::ChamferCorner::Point,
                    T,
                )
            },
            Some(mitred),
        ),
        (
            "fillet_edges_corner, rolling",
            |m, near, target, _, _| {
                let edges = explore_unique(m, near, ShapeType::Edge).unwrap();
                ogeom_fillet::fillet_edges_corner(
                    m,
                    target,
                    &edges,
                    1.0,
                    ogeom_fillet::FilletCorner::Rolling,
                    T,
                )
            },
            Some(swept),
        ),
        (
            "fillet_edge_variable",
            |m, _, target, edge, _| {
                ogeom_fillet::fillet_edge_variable(m, target, edge, 1.0, 2.0, T)
            },
            None,
        ),
        (
            "fillet_edges_variable",
            |m, _, target, edge, _| {
                ogeom_fillet::fillet_edges_variable(
                    m,
                    target,
                    std::slice::from_ref(edge),
                    &ogeom_fillet::RadiusLaw::Linear {
                        start: 1.0,
                        end: 2.0,
                    },
                    T,
                )
            },
            None,
        ),
        (
            "fillet_edges_cornered",
            |m, _, target, edge, _| {
                ogeom_fillet::fillet_edges_cornered(
                    m,
                    target,
                    std::slice::from_ref(edge),
                    1.0,
                    &[],
                    T,
                )
            },
            None,
        ),
        (
            "fillet_edge_distances",
            |m, _, target, edge, top| {
                ogeom_fillet::fillet_edge_distances(
                    m,
                    target,
                    edge,
                    top,
                    1.0,
                    2.0,
                    ogeom_fillet::RoundProfile::Elliptical,
                    T,
                )
            },
            None,
        ),
    ];
    for (name, entry, exact) in entries {
        let mut m = Model::new();
        let (near, far, both) = touching_pair(&mut m);
        let (edge, top) = shared_side_edge_and_top(&m, &near);
        let want = exact.unwrap_or_else(|| {
            let mut alone = m.clone();
            let built = entry(&mut alone, &near, &near, &edge, &top).unwrap();
            volume(&alone, &built.shape)
        });
        let built =
            entry(&mut m, &near, &both, &edge, &top).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_far_kept(&m, &built.shape, &far, want);
    }
}
