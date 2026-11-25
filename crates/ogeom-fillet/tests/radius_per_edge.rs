//! A radius per edge in one call.
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{fillet_edges, fillet_edges_with};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn bx(m: &mut Model, s: [f64; 3]) -> Shape {
    ogeom_algo::make_box(m, at(0.0, 0.0, 0.0), (s[0], s[1], s[2]), T)
        .unwrap()
        .shape
}

fn vol(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(
        m,
        s,
        ogeom_mesh::Deflection {
            chord: 1e-3,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

/// The edge of `s` whose middle stands at `p`.
fn edge_at(m: &Model, s: &Shape, p: [f64; 3]) -> Shape {
    let q = Point::new(p[0], p[1], p[2]);
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            c.point_at(0.5 * (a + b), T).unwrap().distance(q) < 1e-9
        })
        .expect("an edge there")
}

/// The section a band of radius `r` takes off a right-angled crease.
fn section(r: f64) -> f64 {
    (1.0 - PI / 4.0) * r * r
}

#[test]
fn one_radius_on_every_edge_is_fillet_edges() {
    let mut m = Model::new();
    let b = bx(&mut m, [20.0, 20.0, 10.0]);
    let all = explore_unique(&m, &b, ShapeType::Edge).unwrap();
    let mut plain = m.clone();
    let a = fillet_edges(&mut plain, &b, &all, 2.0, T).unwrap();
    let with: Vec<(Shape, f64)> = all.iter().map(|e| (e.clone(), 2.0)).collect();
    let w = fillet_edges_with(&mut m, &b, &with, T).unwrap();
    let faces = |m: &Model, s: &Shape| explore_unique(m, s, ShapeType::Face).unwrap().len();
    assert_eq!(faces(&plain, &a.shape), faces(&m, &w.shape));
    assert_eq!(vol(&plain, &a.shape).to_bits(), vol(&m, &w.shape).to_bits());
}

/// Two creases of a box's top turning at a vertex, at radii 1 and 3: the
/// cut takes the union of the two bands, whichever order they are asked
/// in. Where the bands overlap at the vertex, the region under the
/// small band's spandrel within the large band's reach,
/// `∫ (1 − √(1 − (z − 9)²)) (3 − √(9 − (z − 7)²)) dz` over the top unit.
#[test]
fn two_radii_turning_at_a_vertex_take_the_union_of_their_bands() {
    let mut m = Model::new();
    let b = bx(&mut m, [20.0, 20.0, 10.0]);
    let front = edge_at(&m, &b, [10.0, 0.0, 10.0]);
    let left = edge_at(&m, &b, [0.0, 10.0, 10.0]);
    let n = 200_000;
    let h = 1.0 / f64::from(n);
    let overlap: f64 = (0..n)
        .map(|k| {
            let z = 9.0 + (f64::from(k) + 0.5) * h;
            let small = 1.0 - (1.0 - (z - 9.0).powi(2)).sqrt();
            let large = 3.0 - (9.0 - (z - 7.0).powi(2)).sqrt();
            small * large * h
        })
        .sum();
    let expected = 4000.0 - 20.0 * section(1.0) - 20.0 * section(3.0) + overlap;
    let mut volumes = Vec::new();
    for asked in [
        [(front.clone(), 1.0), (left.clone(), 3.0)],
        [(left.clone(), 3.0), (front.clone(), 1.0)],
    ] {
        let mut mm = m.clone();
        let built = fillet_edges_with(&mut mm, &b, &asked, T).unwrap();
        let diagnosis = ogeom_algo::check(&mm, &built.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        assert!(built.history.is_deleted(&front) && built.history.is_deleted(&left));
        volumes.push(vol(&mm, &built.shape));
    }
    for v in &volumes {
        assert!((v - expected).abs() < 1e-3, "{v} against {expected}");
    }
}

/// A box's three edges at one corner at three radii: no ball rolls round
/// it, so each strip ends through its faces' setback corners and one patch
/// closes the corner, tangent to every strip's end, the same whatever the
/// order asked.
#[test]
fn three_radii_at_a_corner_close_with_one_patch_in_any_order() {
    let mut m = Model::new();
    let b = bx(&mut m, [20.0, 20.0, 10.0]);
    let edges = [
        edge_at(&m, &b, [10.0, 0.0, 10.0]),
        edge_at(&m, &b, [0.0, 10.0, 10.0]),
        edge_at(&m, &b, [0.0, 0.0, 5.0]),
    ];
    let radii = [2.0, 3.0, 2.5];
    let asked: Vec<(Shape, f64)> = edges.iter().cloned().zip(radii).collect();
    let mut volumes = Vec::new();
    let mut faces = Vec::new();
    for order in [[0, 1, 2], [2, 0, 1], [1, 2, 0]] {
        let request: Vec<(Shape, f64)> = order.iter().map(|&i| asked[i].clone()).collect();
        let rounded = fillet_edges_with(&mut m, &b, &request, T).unwrap();
        let found = ogeom_algo::check(&m, &rounded.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        // The patch: the face bounded by the three strips' ends, which meet
        // at the faces' setback corners, tangent across each.
        let patch = corner_patch(&m, &b, &rounded, [0.0, 0.0, 10.0]);
        let contacts = ogeom_fillet::analyse_blend(&m, &rounded.shape, &patch, 17, T).unwrap();
        assert_eq!(contacts.len(), 3, "{contacts:?}");
        for c in &contacts {
            assert!(c.tangency_error < 1e-3, "a crease across the patch: {c:?}");
        }
        volumes.push(vol(&m, &rounded.shape));
        faces.push(
            explore_unique(&m, &rounded.shape, ShapeType::Face)
                .unwrap()
                .len(),
        );
    }
    for (v, f) in volumes.iter().zip(&faces) {
        assert!((v - volumes[0]).abs() <= volumes[0] * 1e-9, "{volumes:?}");
        assert_eq!(*f, faces[0], "{faces:?}");
    }
    let mut at = |r: f64| {
        let request: Vec<(Shape, f64)> = edges.iter().cloned().map(|e| (e, r)).collect();
        let s = fillet_edges_with(&mut m, &b, &request, T).unwrap().shape;
        vol(&m, &s)
    };
    let (two, three) = (at(2.0), at(3.0));
    assert!(
        three < volumes[0] && volumes[0] < two,
        "{three} < {} < {two}",
        volumes[0]
    );
}

/// The issue's corner at 1, 3 and 2: the patch closes it tangent to every
/// strip's end and within the block, the same solid whatever the order
/// asked, its volume between the corner rounded at 1 and at 3.
#[test]
fn the_issue_s_three_radii_close_within_the_block_in_any_order() {
    let mut m = Model::new();
    let b = bx(&mut m, [20.0, 20.0, 10.0]);
    let edges = [
        edge_at(&m, &b, [10.0, 0.0, 10.0]),
        edge_at(&m, &b, [0.0, 10.0, 10.0]),
        edge_at(&m, &b, [0.0, 0.0, 5.0]),
    ];
    let asked: Vec<(Shape, f64)> = edges.iter().cloned().zip([1.0, 3.0, 2.0]).collect();
    let mut volumes = Vec::new();
    for order in [[0, 1, 2], [2, 0, 1], [1, 2, 0], [2, 1, 0]] {
        let request: Vec<(Shape, f64)> = order.iter().map(|&i| asked[i].clone()).collect();
        let rounded =
            fillet_edges_with(&mut m, &b, &request, T).unwrap_or_else(|e| panic!("{order:?}: {e}"));
        let found = ogeom_algo::check(&m, &rounded.shape, T).unwrap();
        assert!(found.is_valid(), "{order:?}: {:?}", found.problems);
        let patch = corner_patch(&m, &b, &rounded, [0.0, 0.0, 10.0]);
        let contacts = ogeom_fillet::analyse_blend(&m, &rounded.shape, &patch, 17, T).unwrap();
        assert_eq!(contacts.len(), 3, "{contacts:?}");
        let worst = contacts
            .iter()
            .map(|c| c.tangency_error)
            .fold(0.0, f64::max);
        assert!(worst < 2e-4, "{order:?}: tangent within {worst} rad");
        // Within the block: read on the patch's own mesh, finer than the
        // bands' fit target.
        let mesh = ogeom_mesh::triangulate_face(
            &m,
            &patch,
            ogeom_mesh::Deflection {
                chord: 2e-5,
                ..ogeom_mesh::Deflection::default()
            },
            T,
        )
        .unwrap();
        let past = mesh
            .positions
            .iter()
            .map(|p| (-p.x).max(-p.y).max(p.z - 10.0))
            .fold(f64::NEG_INFINITY, f64::max);
        eprintln!(
            "{order:?}: tangent within {worst} rad, {past} past the block's faces, volume {}",
            vol(&m, &rounded.shape)
        );
        assert!(
            past <= 2e-4,
            "{order:?}: the patch stands {past} outside the block"
        );
        volumes.push(vol(&m, &rounded.shape));
    }
    for v in &volumes {
        assert!((v - volumes[0]).abs() <= volumes[0] * 1e-9, "{volumes:?}");
    }
    let mut at = |r: f64| {
        let request: Vec<(Shape, f64)> = edges.iter().cloned().map(|e| (e, r)).collect();
        let s = fillet_edges_with(&mut m, &b, &request, T).unwrap().shape;
        vol(&m, &s)
    };
    let (one, three) = (at(1.0), at(3.0));
    assert!(
        three < volumes[0] && volumes[0] < one,
        "{three} < {} < {one}",
        volumes[0]
    );
}

/// The face `rounded`'s history makes of `b`'s vertex at `at`: the corner's
/// patch.
fn corner_patch(m: &Model, b: &Shape, rounded: &ogeom_algo::Built, at: [f64; 3]) -> Shape {
    let vertex = explore_unique(m, b, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .find(|v| {
            let p = m
                .node(v)
                .and_then(|n| n.data().as_vertex().map(|d| d.point))
                .unwrap();
            p.distance(ogeom_math::Point::new(at[0], at[1], at[2])) < 1e-9
        })
        .expect("the corner");
    let made = rounded.history.generated(&vertex);
    assert_eq!(made.len(), 1, "one patch at the corner");
    made[0].clone()
}

/// A square pyramid, its apex at (0, 0, 10).
fn pyramid(m: &mut Model) -> Shape {
    let mut solid = ogeom_algo::make_box(m, at(-10.0, -10.0, -10.0), (20.0, 20.0, 20.0), T)
        .unwrap()
        .shape;
    // Four planes leaning in from the box's top edges to the apex, each
    // cutting a half-space off.
    for (nx, ny) in [(1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)] {
        let normal = Vector::new(nx, ny, 1.0);
        let across = Vector::new(-ny, nx, 0.0);
        let side = normal.cross(across) / normal.cross(across).magnitude();
        let corner = Point::new(0.0, 0.0, 10.0) - across * 40.0 - side * 40.0;
        let frame = Frame::new(
            corner,
            Direction::new(normal, T).unwrap(),
            Direction::new(across, T).unwrap(),
            T,
        )
        .unwrap();
        let away = ogeom_algo::make_box(m, frame, (80.0, 80.0, 40.0), T)
            .unwrap()
            .shape;
        solid = ogeom_bool::cut(m, &solid, &away, T).unwrap().shape;
    }
    solid
}

/// A pyramid's apex at alternating radii 0.5 and 1.5: each band stops
/// past where the others reach its section plane, not only past where its
/// own section clears them, so a wide band beside a narrow one does not
/// cross the narrow one's section inside the cell. Closed, valid, tangent,
/// the same in either order (the patch's net to the bit), between the
/// rounds at either radius alone, and the corner's patch proved within the
/// pyramid by the hull of its stored net, under the four slants and within
/// the box's walls and floor. Built by the engine's assembly with no
/// boolean: interpolated, the patch leans 4.7e-3 past a face it is tangent
/// to along its side, and is proved within the material only once its net
/// is held within the corner's planes.
#[test]
fn a_pyramid_apex_at_alternating_radii_closes_with_one_patch() {
    for (small, large) in [(0.5, 1.5)] {
        let mut m = Model::new();
        let solid = pyramid(&mut m);
        let apex = Point::new(0.0, 0.0, 10.0);
        let at_apex: Vec<Shape> = explore_unique(&m, &solid, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                let (c, (a, b)) = ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap();
                c.point_at(a, T).unwrap().distance(apex) < 1e-6
                    || c.point_at(b, T).unwrap().distance(apex) < 1e-6
            })
            .collect();
        assert_eq!(at_apex.len(), 4);
        // Round the apex: alternate edges at the two radii, by direction so
        // the alternation is round the apex whatever order they came in.
        let radius = |e: &Shape| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap();
            let far = [c.point_at(a, T).unwrap(), c.point_at(b, T).unwrap()]
                .into_iter()
                .find(|p| p.distance(apex) > 1e-6)
                .unwrap();
            if (far.x > 0.0) == (far.y > 0.0) {
                small
            } else {
                large
            }
        };
        let asked: Vec<(Shape, f64)> = at_apex.iter().map(|e| (e.clone(), radius(e))).collect();
        let mut volumes = Vec::new();
        let mut nets = Vec::new();
        for order in [[0, 1, 2, 3], [3, 1, 0, 2]] {
            let request: Vec<(Shape, f64)> = order.iter().map(|&i| asked[i].clone()).collect();
            let (rounded, work) = ogeom_bool::work::counted(|| {
                fillet_edges_with(&mut m, &solid, &request, T).unwrap()
            });
            assert_eq!(work.face_pairs, 0, "built by the engine's assembly");
            let found = ogeom_algo::check(&m, &rounded.shape, T).unwrap();
            assert!(found.is_valid(), "{small}/{large}: {:?}", found.problems);
            let patch = explore_unique(&m, &rounded.shape, ShapeType::Face)
                .unwrap()
                .into_iter()
                .find(|f| {
                    m.node(f)
                        .and_then(|n| n.data().as_face())
                        .and_then(|d| m.geometry().surface(d.surface))
                        .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::BSpline(_)))
                })
                .expect("the corner patch");
            // Within the pyramid, proved on the stored surface by the hull
            // of its net (subdivided where the whole net does not show it):
            // under the four slants through the apex, within the box's
            // walls and over its floor.
            let slant = |nx: f64, ny: f64| ogeom_algo::HalfSpaces::Below {
                normal: Vector::new(nx, ny, 1.0) / 2.0_f64.sqrt(),
                at: apex,
            };
            let wall = |nx: f64, ny: f64, nz: f64| ogeom_algo::HalfSpaces::Below {
                normal: Vector::new(nx, ny, nz),
                at: Point::new(nx * 10.0, ny * 10.0, nz * 10.0),
            };
            let pyramid_region = ogeom_algo::HalfSpaces::Every(vec![
                slant(1.0, 0.0),
                slant(-1.0, 0.0),
                slant(0.0, 1.0),
                slant(0.0, -1.0),
                wall(1.0, 0.0, 0.0),
                wall(-1.0, 0.0, 0.0),
                wall(0.0, 1.0, 0.0),
                wall(0.0, -1.0, 0.0),
                wall(0.0, 0.0, -1.0),
            ]);
            let proved =
                ogeom_algo::certify_face_within(&m, &patch, &pyramid_region, 1e-6, T).unwrap();
            assert!(
                proved.is_ok(),
                "{small}/{large}: not proved within: {proved:?}"
            );
            if let Some(ogeom_geom::SurfaceGeometry::BSpline(s)) = m
                .node(&patch)
                .and_then(|n| n.data().as_face())
                .and_then(|d| m.geometry().surface(d.surface))
            {
                nets.push(s.grid().points().to_vec());
            }
            // A counterexample finder beside the proof: the fine mesh.
            let fine = ogeom_mesh::Deflection {
                chord: 1e-5,
                ..ogeom_mesh::Deflection::default()
            };
            let mesh = ogeom_mesh::triangulate_face(&m, &patch, fine, T).unwrap();
            let past = mesh
                .positions
                .iter()
                .map(|q| {
                    [(1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)]
                        .iter()
                        .map(|(nx, ny)| (nx * q.x + ny * q.y + q.z - 10.0) / 2.0_f64.sqrt())
                        .fold(f64::NEG_INFINITY, f64::max)
                })
                .fold(f64::NEG_INFINITY, f64::max);
            assert!(
                past <= 1e-6,
                "{small}/{large}: the patch stands {past} past a face"
            );
            for c in ogeom_fillet::analyse_blend(&m, &rounded.shape, &patch, 17, T).unwrap() {
                assert!(c.tangency_error < 1e-3, "a crease across the patch: {c:?}");
            }
            volumes.push(vol(&m, &rounded.shape));
        }
        assert!(
            (volumes[0] - volumes[1]).abs() <= volumes[0] * 1e-9,
            "{volumes:?}"
        );
        assert_eq!(nets.len(), 2);
        assert!(
            nets[0] == nets[1],
            "the patch's net is the same in either order"
        );
        let mut at = |r: f64| {
            let request: Vec<(Shape, f64)> = at_apex.iter().cloned().map(|e| (e, r)).collect();
            let s = fillet_edges_with(&mut m, &solid, &request, T)
                .unwrap()
                .shape;
            vol(&m, &s)
        };
        let (narrow, wide) = (at(small), at(large));
        assert!(
            wide < volumes[0] && volumes[0] < narrow,
            "{small}/{large}: {wide} < {} < {narrow}",
            volumes[0]
        );
    }
}
