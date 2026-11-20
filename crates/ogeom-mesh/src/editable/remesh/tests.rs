#![allow(
    clippy::unwrap_used,
    clippy::cast_precision_loss,
    reason = "test code; counts are far below 2^52"
)]

use super::*;
use crate::editable::tests::{check, grid};
use ogeom_core::OgeomError;

/// A cube of side `side` from the origin as 12 triangles wound outward,
/// vertex `i` at the corner with `x`, `y`, `z` its bits 0, 1 and 2.
fn cube(side: f64) -> Triangulation {
    let positions = (0..8)
        .map(|i| {
            let bit = |k: u32| f64::from((i >> k) & 1) * side;
            Point::new(bit(0), bit(1), bit(2))
        })
        .collect();
    let triangles = vec![
        [0, 2, 3],
        [0, 3, 1],
        [4, 5, 7],
        [4, 7, 6],
        [0, 1, 5],
        [0, 5, 4],
        [2, 6, 7],
        [2, 7, 3],
        [0, 4, 6],
        [0, 6, 2],
        [1, 3, 7],
        [1, 7, 5],
    ];
    Triangulation {
        positions,
        triangles,
        ..Triangulation::new()
    }
}

/// A unit sphere of 80 triangles: the poles and five rings of eight at
/// uneven heights, each ring turned against the one above, so the
/// triangles range from slivers to broad ones.
fn uneven_sphere() -> Triangulation {
    let heights = [0.2, 0.55, 1.3, 2.1, 2.75_f64];
    let mut positions = vec![Point::new(0.0, 0.0, 1.0)];
    for (k, theta) in heights.iter().enumerate() {
        for j in 0..8 {
            let phi = (f64::from(j) + 0.37 * k as f64) * core::f64::consts::TAU / 8.0;
            positions.push(Point::new(
                theta.sin() * phi.cos(),
                theta.sin() * phi.sin(),
                theta.cos(),
            ));
        }
    }
    positions.push(Point::new(0.0, 0.0, -1.0));
    let at = |k: u32, j: u32| 1 + 8 * k + j % 8;
    let mut triangles = Vec::new();
    for j in 0..8 {
        triangles.push([0, at(0, j), at(0, j + 1)]);
        for k in 0..4 {
            triangles.push([at(k, j), at(k + 1, j), at(k + 1, j + 1)]);
            triangles.push([at(k, j), at(k + 1, j + 1), at(k, j + 1)]);
        }
        triangles.push([41, at(4, j + 1), at(4, j)]);
    }
    Triangulation {
        positions,
        triangles,
        ..Triangulation::new()
    }
}

/// Each edge's length, once per edge.
fn edge_lengths(m: &Triangulation) -> Vec<f64> {
    let mut edges: Vec<(u32, u32)> = m
        .triangles
        .iter()
        .flat_map(|t| (0..3).map(move |k| (t[k].min(t[(k + 1) % 3]), t[k].max(t[(k + 1) % 3]))))
        .collect();
    edges.sort_unstable();
    edges.dedup();
    edges
        .iter()
        .map(|&(a, b)| m.positions[a as usize].distance(m.positions[b as usize]))
        .collect()
}

/// The share of `lengths` from `low` to `high`.
fn share_within(lengths: &[f64], low: f64, high: f64) -> f64 {
    lengths
        .iter()
        .filter(|&&l| (low..=high).contains(&l))
        .count() as f64
        / lengths.len() as f64
}

fn remesh(edge: f64) -> Remesh<'static> {
    Remesh {
        edge,
        passes: 5,
        keep_edges: None,
        region: None,
        max_faces: None,
    }
}

/// The distance from a point inside a box from the origin to `side` to the
/// box's surface; `None` for a point outside it.
fn to_cube(p: Point, side: f64) -> Option<f64> {
    let c = [p.x, p.y, p.z];
    c.iter()
        .all(|x| (0.0..=side).contains(x))
        .then(|| c.iter().map(|&x| x.min(side - x)).fold(f64::MAX, f64::min))
}

#[test]
fn a_cube_remeshed_keeps_its_edges_corners_faces_and_volume() {
    let side = 10.0;
    let input = cube(side);
    let mut e = Editable::from(&input, 0.0).unwrap();
    let report = e
        .remesh(&Remesh {
            keep_edges: Some(30_f64.to_radians()),
            ..remesh(1.0)
        })
        .unwrap();
    check(&e);
    assert!(report.split > 0 && report.collapsed > 0 && report.flipped > 0);
    assert!(e.is_closed() && e.non_manifold().is_empty());
    let out = e.to_triangulation();
    assert!(out.is_closed());

    let lengths = edge_lengths(&out);
    let share = share_within(&lengths, 0.5, 1.5);
    assert!(share >= 0.95, "{share} of the edges within 0.5..1.5");

    // Every vertex on the surface, and every face on one of its planes.
    for p in &out.positions {
        let d = to_cube(*p, side).unwrap();
        assert!(d <= 1e-9, "{p:?} {d} off the surface");
    }
    for t in &out.triangles {
        let [a, b, c] = t.map(|v| out.positions[v as usize]);
        let on_one = (0..3).any(|axis| {
            let x = |p: Point| [p.x, p.y, p.z][axis];
            [0.0, side]
                .iter()
                .any(|&s| [a, b, c].iter().all(|&p| (x(p) - s).abs() <= 1e-9))
        });
        assert!(on_one, "{t:?} leaves the cube's planes");
    }

    // The corners stay, with their sources, and the twelve edges are
    // covered by edges of the mesh, end to end.
    for (s, corner) in input.positions.iter().enumerate() {
        let at: Vec<u32> = (0..index(out.positions.len()))
            .filter(|&v| out.positions[v as usize] == *corner)
            .collect();
        assert_eq!(at.len(), 1, "corner {corner:?}");
        assert_eq!(e.source_vertices(at[0]), &[index(s)]);
    }
    let mut cube_edges = Vec::new();
    for a in 0..8_u32 {
        for k in 0..3 {
            let b = a | (1 << k);
            if b != a {
                cube_edges.push((input.positions[a as usize], input.positions[b as usize]));
            }
        }
    }
    assert_eq!(cube_edges.len(), 12);
    let mut edges: Vec<(u32, u32)> = out
        .triangles
        .iter()
        .flat_map(|t| (0..3).map(move |k| (t[k].min(t[(k + 1) % 3]), t[k].max(t[(k + 1) % 3]))))
        .collect();
    edges.sort_unstable();
    edges.dedup();
    for (from, to) in cube_edges {
        let along = (to - from) / side;
        let on = |p: Point| {
            let d = p - from;
            (d - along * d.dot(along)).magnitude() <= 1e-9
        };
        let covered: f64 = edges
            .iter()
            .map(|&(a, b)| (out.positions[a as usize], out.positions[b as usize]))
            .filter(|&(p, q)| on(p) && on(q))
            .map(|(p, q)| p.distance(q))
            .sum();
        assert!(
            (covered - side).abs() <= 1e-9,
            "{from:?}..{to:?} covered {covered}"
        );
    }

    let volume = out.volume();
    assert!((volume - 1000.0).abs() <= 1e-6 * 1000.0, "volume {volume}");
    // On the cube's planes throughout, the new surface is the old.
    assert!(report.max <= 1e-9, "{report:?}");
}

#[test]
fn an_uneven_sphere_remeshed_is_even_closed_and_on_its_surface() {
    let input = uneven_sphere();
    assert!(input.is_closed() && input.triangles.len() == 80);
    let mut e = Editable::from(&input, 0.0).unwrap();
    let report = e.remesh(&remesh(0.1)).unwrap();
    check(&e);
    assert!(e.is_closed() && e.non_manifold().is_empty());
    let out = e.to_triangulation();
    assert!(out.triangles.len() > 1000);

    let lengths = edge_lengths(&out);
    let share = share_within(&lengths, 0.05, 0.15);
    assert!(share >= 0.95, "{share} of the edges within 0.05..0.15");
    let within = share_within(&lengths, 0.07, 0.13);
    assert!(within >= 0.9, "{within} of the edges within 30 % of 0.1");

    let bvh = Bvh::build(&input).unwrap();
    let mut farthest = 0.0_f64;
    for p in &out.positions {
        farthest = farthest.max(bvh.nearest(*p, f64::INFINITY).unwrap().distance);
    }
    assert!(farthest <= 1e-3, "a vertex {farthest} off the input");
    // The report measures at the vertices and at the faces' centroids. A
    // face across a fold of the input, its corners on either side, lies
    // inside by up to its length times the sine of half the fold.
    let fold = input
        .triangles
        .iter()
        .enumerate()
        .flat_map(|(i, s)| {
            input.triangles.iter().skip(i + 1).filter_map(move |t| {
                let shared = s.iter().filter(|v| t.contains(v)).count();
                (shared == 2).then_some((*s, *t))
            })
        })
        .map(|(s, t)| {
            let n = |t: [u32; 3]| {
                let [a, b, c] = t.map(|v| input.positions[v as usize]);
                let n = (b - a).cross(c - a);
                n / n.magnitude()
            };
            n(s).dot(n(t)).clamp(-1.0, 1.0).acos()
        })
        .fold(0.0, f64::max);
    let longest = lengths.iter().copied().fold(0.0, f64::max);
    assert!(report.max >= farthest, "{report:?}");
    assert!(
        report.max <= longest * (fold / 2.0).sin(),
        "{report:?} {fold}"
    );
    assert!(report.rms <= report.max);

    // The input is convex and the result's vertices lie on it, so the
    // result lies inside it; its chords across the input's edges cut off
    // a little, no more than the deviation over the area.
    let (before, after) = (input.volume(), out.volume());
    assert!(after < before);
    let area: f64 = out
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|v| out.positions[v as usize]);
            (b - a).cross(c - a).magnitude() / 2.0
        })
        .sum();
    assert!(before - after <= report.max * area, "{before} to {after}");
}

#[test]
fn remeshing_is_deterministic_and_carries_sources() {
    let input = uneven_sphere();
    let run = || {
        let mut e = Editable::from(&input, 0.0).unwrap();
        let report = e.remesh(&remesh(0.15)).unwrap();
        (e, report)
    };
    let (first, report) = run();
    let (second, again) = run();
    assert_eq!(first, second);
    assert_eq!(report, again);

    // A collapse merges source sets and a split makes a vertex with none,
    // so each source vertex is behind at most one vertex.
    let mut behind = vec![0; input.positions.len()];
    for v in 0..index(first.vertex_count()) {
        if first.position(v).is_some() {
            for &s in first.source_vertices(v) {
                behind[s as usize] += 1;
            }
        }
    }
    assert!(behind.iter().all(|&n| n <= 1));
    assert!(behind.contains(&1));
}

#[test]
fn a_cube_subdivided_flat_twice() {
    let input = cube(10.0);
    let mut e = Editable::from(&input, 0.0).unwrap();
    let report = e.subdivide(2, Subdivision::Flat, None).unwrap();
    check(&e);
    assert_eq!(report.faces, 192);
    assert_eq!(report.moved, 0.0);
    let out = e.to_triangulation();
    assert_eq!(out.triangles.len(), 192);
    assert_eq!(report.added, out.positions.len() - 8);
    assert!(e.is_closed() && e.non_manifold().is_empty());
    // Every coordinate is a multiple of 2.5, so six times the volume sums
    // exactly.
    let six = out
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|v| out.positions[v as usize].to_vector());
            a.dot(b.cross(c))
        })
        .sum::<f64>();
    assert_eq!(six, 6000.0);
    for v in 0..8_u32 {
        assert_eq!(e.position(v), Some(input.positions[v as usize]));
        assert_eq!(e.source_vertices(v), &[v]);
    }
    for v in 8..index(out.positions.len()) {
        assert!(e.source_vertices(v).is_empty());
        assert_eq!(to_cube(out.positions[v as usize], 10.0), Some(0.0));
    }
}

#[test]
fn an_open_sheet_keeps_its_rim() {
    // Ten by ten unit cells, the inner vertices pushed about in the plane.
    let mut input = grid(10);
    for (k, p) in input.positions.iter_mut().enumerate() {
        if p.x > 0.0 && p.x < 10.0 && p.y > 0.0 && p.y < 10.0 {
            let t = k as f64;
            p.x += 0.3 * (t * 1.7).sin();
            p.y += 0.3 * (t * 2.3).cos();
        }
    }
    let mut e = Editable::from(&input, 0.0).unwrap();
    let before = e.rims();
    e.remesh(&remesh(0.4)).unwrap();
    check(&e);
    let rims = e.rims();
    assert_eq!((before.len(), rims.len()), (1, 1));
    assert!(rims[0].edges.len() > before[0].edges.len());
    let (was, is) = (before[0].length, rims[0].length);
    assert!((is - was).abs() <= 1e-9 * was, "{was} to {is}");
    let out = e.to_triangulation();
    assert!(out.positions.iter().all(|p| p.z == 0.0));
    // The rim's vertices are all on the square's sides.
    for &h in &rims[0].edges {
        let p = e.position(e.edge(h).unwrap()[0]).unwrap();
        assert!(
            [p.x, p.y]
                .iter()
                .any(|&c| c.abs() <= 1e-12 || (c - 10.0).abs() <= 1e-12)
        );
    }
    let share = share_within(&edge_lengths(&out), 0.2, 0.6);
    assert!(share >= 0.95, "{share}");
}

#[test]
fn a_region_leaves_the_faces_with_a_corner_outside_it() {
    let mut e = Editable::from(&cube(10.0), 0.0).unwrap();
    e.subdivide(3, Subdivision::Flat, None).unwrap();
    let before = e.to_triangulation();
    assert_eq!(before.triangles.len(), 768);
    let region: Vec<bool> = (0..index(e.vertex_count()))
        .map(|v| e.position(v).is_some_and(|p| p.x <= 5.0))
        .collect();
    let bits = |m: &Triangulation, t: &[u32; 3]| {
        let mut corners = t.map(|v| {
            let p = m.positions[v as usize];
            [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
        });
        // The same face from its lowest corner.
        let first = (0..3).min_by_key(|&k| corners[k]).unwrap();
        corners.rotate_left(first);
        corners
    };
    let outside: Vec<[[u64; 3]; 3]> = before
        .triangles
        .iter()
        .filter(|t| t.iter().any(|&v| before.positions[v as usize].x > 5.0))
        .map(|t| bits(&before, t))
        .collect();
    assert!(!outside.is_empty() && outside.len() < 768);

    let report = e
        .remesh(&Remesh {
            region: Some(&region),
            keep_edges: Some(30_f64.to_radians()),
            ..remesh(0.6)
        })
        .unwrap();
    check(&e);
    assert!(report.split > 0 && report.collapsed > 0);
    assert!(e.is_closed() && e.non_manifold().is_empty());
    let after = e.to_triangulation();
    let mut kept: Vec<[[u64; 3]; 3]> = after
        .triangles
        .iter()
        .filter(|t| t.iter().any(|&v| after.positions[v as usize].x > 5.0))
        .map(|t| bits(&after, t))
        .collect();
    let mut expected = outside;
    kept.sort_unstable();
    expected.sort_unstable();
    assert_eq!(kept, expected);
    assert_ne!(after.triangles.len(), 768);
    assert!((after.volume() - 1000.0).abs() <= 1e-9 * 1000.0);
}

#[test]
fn too_many_faces_are_refused_before_any_work() {
    let input = cube(10.0);
    let mut e = Editable::from(&input, 0.0).unwrap();
    let untouched = e.clone();
    // 600 mm² in triangles of side 0.1 is about 138 564 faces.
    let err = e
        .remesh(&Remesh {
            max_faces: Some(100_000),
            ..remesh(0.1)
        })
        .unwrap_err();
    assert!(matches!(err, OgeomError::Range(_)));
    assert!(
        err.to_string().contains("138564 faces, past the 100000"),
        "{err}"
    );
    assert_eq!(e, untouched);
    // Within the limit it runs and stays within it.
    let done = e
        .remesh(&Remesh {
            max_faces: Some(2_000),
            ..remesh(1.0)
        })
        .unwrap();
    assert!(done.split > 0);
    assert!(e.to_triangulation().triangles.len() <= 2_000);

    let mut e = untouched.clone();
    let err = e.subdivide(3, Subdivision::Flat, Some(767)).unwrap_err();
    assert!(err.to_string().contains("768 faces, past the 767"), "{err}");
    assert_eq!(e, untouched);
    let err = e.subdivide(40, Subdivision::Flat, None).unwrap_err();
    assert!(err.to_string().contains("past u32 numbering"), "{err}");
    assert_eq!(e, untouched);
}

#[test]
fn bad_requests_are_refused() {
    let mut e = Editable::from(&cube(10.0), 0.0).unwrap();
    let untouched = e.clone();
    for edge in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let err = e.remesh(&remesh(edge)).unwrap_err();
        assert!(matches!(err, OgeomError::Numeric(_)), "{err}");
    }
    for angle in [-0.1, 3.2, f64::NAN] {
        let err = e
            .remesh(&Remesh {
                keep_edges: Some(angle),
                ..remesh(1.0)
            })
            .unwrap_err();
        assert!(err.to_string().contains("crease angle"), "{err}");
        let err = e
            .subdivide(
                1,
                Subdivision::Smooth {
                    keep_edges: Some(angle),
                },
                None,
            )
            .unwrap_err();
        assert!(err.to_string().contains("crease angle"), "{err}");
    }
    let short = [true; 3];
    let err = e
        .remesh(&Remesh {
            region: Some(&short),
            ..remesh(1.0)
        })
        .unwrap_err();
    assert!(matches!(err, OgeomError::Range(_)));
    assert!(
        err.to_string().contains("3 entries for 8 vertices"),
        "{err}"
    );
    assert_eq!(e, untouched);
}

/// The regular octahedron about the origin, wound outward.
fn octahedron() -> Triangulation {
    let positions = vec![
        Point::new(1.0, 0.0, 0.0),
        Point::new(-1.0, 0.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
        Point::new(0.0, -1.0, 0.0),
        Point::new(0.0, 0.0, 1.0),
        Point::new(0.0, 0.0, -1.0),
    ];
    let triangles = vec![
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5],
    ];
    Triangulation {
        positions,
        triangles,
        ..Triangulation::new()
    }
}

#[test]
fn loop_subdivision_tends_to_the_limit_surface() {
    let input = octahedron();
    assert!(input.is_closed() && input.volume() > 0.0);
    let smooth = Subdivision::Smooth { keep_edges: None };
    let mut e = Editable::from(&input, 0.0).unwrap();
    let report = e.subdivide(1, smooth, None).unwrap();
    check(&e);
    assert_eq!(report.faces, 32);
    // Valence four: β = (5/8 - (3/8)²)/4 = 31/256, so the apex goes to
    // 1 - 4β = 0.515625; the midpoint of an edge from it to the equator
    // is 3/8 of each end and 1/8 of each corner across, which cancel.
    assert_eq!(e.position(4), Some(Point::new(0.0, 0.0, 0.515_625)));
    let out = e.to_triangulation();
    assert!(out.positions.contains(&Point::new(0.375, 0.0, 0.375)));
    assert_eq!(report.moved, 1.0 - 0.515_625);

    // Loop's limit mask: the apex's limit point is (1 - nω) itself and ω
    // each neighbour, ω = 1/(3/(8β) + n) = 31/220.
    let limit = 1.0 - 4.0 * 31.0 / 220.0;
    let mut gaps = vec![(e.position(4).unwrap().z - limit).abs()];
    for _ in 0..4 {
        e.subdivide(1, smooth, None).unwrap();
        let p = e.position(4).unwrap();
        assert_eq!((p.x, p.y), (0.0, 0.0));
        gaps.push((p.z - limit).abs());
    }
    check(&e);
    assert!(e.is_closed());
    assert!(gaps.windows(2).all(|w| w[1] < 0.75 * w[0]), "{gaps:?}");
    assert!(gaps[4] < 0.1 * gaps[0], "{gaps:?}");
}

#[test]
fn smooth_subdivision_keeps_creases_and_rims_sharp() {
    // The cube's edges are creases at 30°: its corners stay, its edges stay
    // straight, and only its faces round.
    let side = 10.0;
    let input = cube(side);
    let mut e = Editable::from(&input, 0.0).unwrap();
    e.subdivide(
        2,
        Subdivision::Smooth {
            keep_edges: Some(30_f64.to_radians()),
        },
        None,
    )
    .unwrap();
    check(&e);
    assert!(e.is_closed());
    for v in 0..8_u32 {
        assert_eq!(e.position(v), Some(input.positions[v as usize]));
    }
    let out = e.to_triangulation();
    let on_edges = out
        .positions
        .iter()
        .filter(|p| {
            [p.x, p.y, p.z]
                .iter()
                .filter(|&&c| c == 0.0 || c == side)
                .count()
                >= 2
        })
        .count();
    // Each edge cut into four: three new vertices on each, and the corners.
    assert_eq!(on_edges, 12 * 3 + 8);
    assert!(out.volume() < 1000.0);

    // A flat sheet stays flat, its corners stay, and its rim stays on its
    // sides.
    let mut e = Editable::from(&grid(4), 0.0).unwrap();
    e.subdivide(2, Subdivision::Smooth { keep_edges: None }, None)
        .unwrap();
    check(&e);
    let out = e.to_triangulation();
    assert!(out.positions.iter().all(|p| p.z == 0.0));
    let rims = e.rims();
    assert_eq!(rims.len(), 1);
    assert!((rims[0].length - 16.0).abs() <= 1e-12);
    for corner in [(0.0, 0.0), (4.0, 0.0), (0.0, 4.0), (4.0, 4.0)] {
        assert!(out.positions.contains(&Point::new(corner.0, corner.1, 0.0)));
    }
}
