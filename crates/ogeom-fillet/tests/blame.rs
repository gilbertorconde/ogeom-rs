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
    // The same plate chamfered 2.5 is no refusal: each horizontal bevel runs
    // past the 2 mm side and rolls over onto the face beyond as a notch, its
    // plane carrying on (ADR 0001, amendment "M5b"), so the plate is the
    // polytope of its planes, its bevels' and its corners' facets.
    let specs: Vec<(Shape, Chamfer)> = all
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(2.5)))
        .collect();
    let built = chamfer_edges_checked(&mut m, &plate, &specs, T).unwrap();
    assert!(ogeom_algo::check(&m, &built.shape, T).unwrap().is_valid());
    let (want, faces) = notched_plate(2.5);
    let got = vol(&m, &built.shape);
    assert!((got - want).abs() <= 1e-9 * want, "{got} against {want}");
    assert_eq!(
        explore_unique(&m, &built.shape, ShapeType::Face)
            .unwrap()
            .len(),
        faces
    );
}

/// The 20 × 20 × 2 plate with every edge chamfered `d`, as the convex
/// polytope of the box's planes, each edge's bevel plane (`d` in along both
/// its faces) and each corner's facet plane (through its three faces'
/// setback corners): its volume and its face count. Independent of the
/// engine.
fn notched_plate(d: f64) -> (f64, usize) {
    use ogeom_math::Vector;
    let size = Vector::new(20.0, 20.0, 2.0);
    let (x, y, z) = (Vector::X, Vector::Y, Vector::Z);
    let mut planes: Vec<(Vector, f64)> = vec![
        (-x, 0.0),
        (x, size.x),
        (-y, 0.0),
        (y, size.y),
        (-z, 0.0),
        (z, size.z),
    ];
    // The corner a sign per axis names: the high end where the sign is one.
    let corner = |s: Vector| {
        Vector::new(
            if s.x > 0.0 { size.x } else { 0.0 },
            if s.y > 0.0 { size.y } else { 0.0 },
            if s.z > 0.0 { size.z } else { 0.0 },
        )
    };
    for (a, b) in [(x, y), (x, z), (y, z)] {
        for sa in [-1.0, 1.0] {
            for sb in [-1.0, 1.0] {
                let (na, nb) = (a * sa, b * sb);
                let n = (na + nb) / (na + nb).magnitude();
                planes.push((n, n.dot(corner(na + nb) - nb * d)));
            }
        }
    }
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            for sz in [-1.0, 1.0] {
                let c = corner(Vector::new(sx, sy, sz));
                let (ix, iy, iz) = (x * (sx * d), y * (sy * d), z * (sz * d));
                let (p, q, r) = (c - ix - iy, c - ix - iz, c - iy - iz);
                let mut n = (q - p).cross(r - p);
                n = n / n.magnitude();
                if n.dot(Vector::new(sx, sy, sz)) < 0.0 {
                    n = -n;
                }
                planes.push((n, n.dot(p)));
            }
        }
    }
    convex_polytope(&planes)
}

/// The convex polytope `planes` bound (each an outward unit normal and its
/// offset): its volume, by the divergence theorem over its faces, and its
/// face count.
fn convex_polytope(planes: &[(ogeom_math::Vector, f64)]) -> (f64, usize) {
    use ogeom_math::Vector;
    let mut corners: Vec<Vector> = Vec::new();
    for i in 0..planes.len() {
        for j in i + 1..planes.len() {
            for k in j + 1..planes.len() {
                let (a, b, c) = (planes[i], planes[j], planes[k]);
                let det = a.0.dot(b.0.cross(c.0));
                if det.abs() < 1e-12 {
                    continue;
                }
                let p = (b.0.cross(c.0) * a.1 + c.0.cross(a.0) * b.1 + a.0.cross(b.0) * c.1) / det;
                if planes.iter().all(|(n, o)| n.dot(p) <= o + 1e-9)
                    && !corners.iter().any(|q| (*q - p).magnitude() < 1e-9)
                {
                    corners.push(p);
                }
            }
        }
    }
    let (mut volume, mut faces) = (0.0, 0);
    for (n, o) in planes {
        let on: Vec<Vector> = corners
            .iter()
            .copied()
            .filter(|p| (n.dot(*p) - o).abs() < 1e-9)
            .collect();
        if on.len() < 3 {
            continue;
        }
        #[expect(clippy::cast_precision_loss, reason = "a face's few corners")]
        let centre = on.iter().fold(Vector::ZERO, |s, p| s + *p) / on.len() as f64;
        let u = (on[0] - centre) / (on[0] - centre).magnitude();
        let w = n.cross(u);
        let mut turn: Vec<(f64, Vector)> = on
            .iter()
            .map(|p| ((*p - centre).dot(w).atan2((*p - centre).dot(u)), *p))
            .collect();
        turn.sort_by(|a, b| a.0.total_cmp(&b.0));
        let area: f64 = (0..turn.len())
            .map(|k| {
                (turn[k].1 - centre)
                    .cross(turn[(k + 1) % turn.len()].1 - centre)
                    .dot(*n)
            })
            .sum::<f64>()
            / 2.0;
        if area > 1e-12 {
            volume += area * o / 3.0;
            faces += 1;
        }
    }
    (volume, faces)
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
fn a_plate_chamfers_partially_refusing_nothing_where_it_notches() {
    // A chamfer of 2.5 runs past the 2 mm sides. That used to refuse every
    // horizontal edge; each now rolls over as a notch (ADR 0001, amendment
    // "M5b"), so the partial build refuses nothing and is the whole
    // request: the polytope of the plate's, bevels' and facets' planes.
    let mut m = Model::new();
    let plate = bx(&mut m, [0.0; 3], [20.0, 20.0, 2.0]);
    let all = explore_unique(&m, &plate, ShapeType::Edge).unwrap();
    let specs: Vec<(Shape, Chamfer)> = all
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(2.5)))
        .collect();
    let source = m.clone();
    let (built, refused) = chamfer_edges_partial(&mut m, &plate, &specs, T).unwrap();
    assert!(refused.is_empty(), "{refused:?}");
    let mut again = source.clone();
    let direct = chamfer_edges_with(&mut again, &plate, &specs, T).unwrap();
    let (v, w) = (vol(&m, &built.shape), vol(&again, &direct.shape));
    assert!((v - w).abs() <= 1e-9 * w, "{v} against {w}");
    let (want, _) = notched_plate(2.5);
    assert!((v - want).abs() <= 1e-9 * want, "{v} against {want}");
    // A request still partly refused: an edge of another solid is named,
    // refused by itself, and the rest is the same notched plate.
    let other = bx(&mut m, [40.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let stray = explore_unique(&m, &other, ShapeType::Edge).unwrap()[0].clone();
    let mut asked = specs.clone();
    asked.push((stray, Chamfer::Symmetric(0.2)));
    let (built, refused) = chamfer_edges_partial(&mut m, &plate, &asked, T).unwrap();
    let refused_ix: Vec<usize> = refused.iter().map(|r| r.edges[0]).collect();
    assert_eq!(refused_ix, vec![specs.len()], "only the stray edge");
    let v = vol(&m, &built.shape);
    assert!((v - want).abs() <= 1e-9 * want, "{v} against {want}");
}
