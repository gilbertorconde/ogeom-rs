//! Refusals naming the request's edges they concern, and partial builds
//! leaving the refused edges sharp.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{
    Chamfer, chamfer_edges_checked, chamfer_edges_partial, chamfer_edges_with, fillet_edges,
    fillet_edges_checked, fillet_edges_partial,
};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn bx(m: &mut Model, o: [f64; 3], s: [f64; 3]) -> Shape {
    ogeom_algo::make_box(m, at(o[0], o[1], o[2]), (s[0], s[1], s[2]), T)
        .unwrap()
        .shape
}

fn vol(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// The edge's two ends.
fn ends(m: &Model, e: &Shape) -> [Point; 2] {
    let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
    [c.point_at(a, T).unwrap(), c.point_at(b, T).unwrap()]
}

/// Whether an edge of the 2 mm plate is horizontal (on the top or bottom).
fn horizontal(m: &Model, e: &Shape) -> bool {
    let [a, b] = ends(m, e);
    (a.z - b.z).abs() < 1e-9
}

/// Whether two horizontal plate edges are the same edge across the wall:
/// one on top, one on the bottom, over each other.
fn across_the_wall(m: &Model, e: &Shape, f: &Shape) -> bool {
    let ([a, b], [c, d]) = (ends(m, e), ends(m, f));
    let flat = |p: Point| (p.x, p.y);
    let same = |p: Point, q: Point| (flat(p).0 - flat(q).0).hypot(flat(p).1 - flat(q).1) < 1e-9;
    (a.z - c.z).abs() > 1.0 && ((same(a, c) && same(b, d)) || (same(a, d) && same(b, c)))
}

#[test]
fn a_plate_s_overlap_names_an_edge_and_the_one_across_the_wall() {
    let mut m = Model::new();
    let plate = bx(&mut m, [0.0; 3], [20.0, 20.0, 2.0]);
    let all = explore_unique(&m, &plate, ShapeType::Edge).unwrap();
    let refusal = fillet_edges_checked(&mut m, &plate, &all, 1.5, T).unwrap_err();
    assert!(refusal.edges.len() >= 2, "{refusal}");
    let (e, f) = (&all[refusal.edges[0]], &all[refusal.edges[1]]);
    assert!(horizontal(&m, e) && across_the_wall(&m, e, f), "{refusal}");
    // The same rounding, chamfered: every wedge seats, and none applies
    // past the face it runs past.
    let specs: Vec<(Shape, Chamfer)> = all
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(2.5)))
        .collect();
    let refusal = chamfer_edges_checked(&mut m, &plate, &specs, T).unwrap_err();
    assert!(!refusal.edges.is_empty(), "{refusal}");
}

#[test]
fn a_notched_box_s_refusal_names_the_notch_edge() {
    let mut m = Model::new();
    let b = bx(&mut m, [0.0; 3], [20.0, 20.0, 10.0]);
    let notch = bx(&mut m, [9.75, -1.0, 8.0], [0.5, 5.0, 5.0]);
    let s = ogeom_bool::cut(&mut m, &b, &notch, T).unwrap().shape;
    let all = explore_unique(&m, &s, ShapeType::Edge).unwrap();
    let refusal = fillet_edges_checked(&mut m, &s, &all, 1.0, T).unwrap_err();
    let named = &all[refusal.edges[0]];
    // The notch's floor edge along y at x = 9.75, z = 8.
    let [a, b] = ends(&m, named);
    assert!(
        (a.x - 9.75).abs() < 1e-9
            && (b.x - 9.75).abs() < 1e-9
            && (a.z - 8.0).abs() < 1e-9
            && (b.z - 8.0).abs() < 1e-9,
        "{refusal}: {a:?} {b:?}"
    );
}

#[test]
fn a_request_that_builds_builds_the_same_checked() {
    let mut m = Model::new();
    let plate = bx(&mut m, [0.0; 3], [20.0, 20.0, 2.0]);
    let all = explore_unique(&m, &plate, ShapeType::Edge).unwrap();
    let mut plain = m.clone();
    let a = fillet_edges(&mut plain, &plate, &all, 0.5, T).unwrap();
    let b = fillet_edges_checked(&mut m, &plate, &all, 0.5, T).unwrap();
    assert!((vol(&plain, &a.shape) - vol(&m, &b.shape)).abs() <= 1e-12);
}

#[test]
fn a_plate_rounds_partially_leaving_one_edge_of_each_pair_sharp() {
    let mut m = Model::new();
    let plate = bx(&mut m, [0.0; 3], [20.0, 20.0, 2.0]);
    let all = explore_unique(&m, &plate, ShapeType::Edge).unwrap();
    let source = m.clone();
    let (built, refused) = fillet_edges_partial(&mut m, &plate, &all, 1.5, T).unwrap();
    assert!(!refused.is_empty());
    let refused_ix: Vec<usize> = refused.iter().map(|r| r.edges[0]).collect();
    let accepted: Vec<usize> = (0..all.len()).filter(|i| !refused_ix.contains(i)).collect();
    for r in &refused {
        let e = &all[r.edges[0]];
        assert!(horizontal(&source, e), "{r}");
        assert!(
            accepted
                .iter()
                .any(|&j| across_the_wall(&source, e, &all[j])),
            "{r}: its partner across the wall is not accepted"
        );
    }
    // The build is fillet_edges on exactly the accepted edges.
    let mut again = source.clone();
    let picked: Vec<Shape> = accepted.iter().map(|&i| all[i].clone()).collect();
    let direct = fillet_edges(&mut again, &plate, &picked, 1.5, T).unwrap();
    let (v, w) = (vol(&m, &built.shape), vol(&again, &direct.shape));
    assert!((v - w).abs() <= 1e-9 * w, "{v} against {w}");
    assert!(ogeom_algo::check(&m, &built.shape, T).unwrap().is_valid());
}

#[test]
fn a_plate_chamfers_partially_leaving_what_runs_past_its_face_sharp() {
    // A chamfer of 2.5 runs past the 2 mm sides: every horizontal edge is
    // refused, and every vertical one, 20 mm faces either side, bevels.
    let mut m = Model::new();
    let plate = bx(&mut m, [0.0; 3], [20.0, 20.0, 2.0]);
    let all = explore_unique(&m, &plate, ShapeType::Edge).unwrap();
    let specs: Vec<(Shape, Chamfer)> = all
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(2.5)))
        .collect();
    let source = m.clone();
    let (built, refused) = chamfer_edges_partial(&mut m, &plate, &specs, T).unwrap();
    let refused_ix: Vec<usize> = refused.iter().map(|r| r.edges[0]).collect();
    for (i, e) in all.iter().enumerate() {
        assert_eq!(refused_ix.contains(&i), horizontal(&source, e), "edge {i}");
    }
    let accepted: Vec<(Shape, Chamfer)> = (0..all.len())
        .filter(|i| !refused_ix.contains(i))
        .map(|i| specs[i].clone())
        .collect();
    let mut again = source.clone();
    let direct = chamfer_edges_with(&mut again, &plate, &accepted, T).unwrap();
    let (v, w) = (vol(&m, &built.shape), vol(&again, &direct.shape));
    assert!((v - w).abs() <= 1e-9 * w, "{v} against {w}");
}
