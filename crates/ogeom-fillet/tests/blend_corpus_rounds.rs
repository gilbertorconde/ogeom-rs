//! The acceptance corpus for rounds on planar faces (ADR 0001, M2, §12).
//!
//! Each case asks the public calls the engine will answer, and holds the
//! answer to an exact oracle where one exists, to bounds where the result
//! has a general vertex patch, and on every case to:
//!
//! - every permutation of the request giving the bit-identical solid;
//! - validity;
//! - the region guard: a convex round's material stays inside the original
//!   solid, a concave one's within the round's reach of its edge;
//! - no overhang: no blend face looks down past the original faces beside
//!   it by more than the fill's wobble.
//!
//! The cases run only with `OGEOM_M2_CORPUS` set, until rounds on planar
//! faces route through the engine: today's path fails them on the order
//! of the request and refuses some, and the class switches only once they
//! pass on the engine (ADR 0001, §11). Without it each returns at once.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom_core::Tolerances;
use ogeom_fillet::{Blend, Chamfer, ChamferCorner};
use ogeom_math::{Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// Whether the corpus is asked for: until planar rounds route through the
/// engine, only where `OGEOM_M2_CORPUS` is set.
fn asked() -> bool {
    std::env::var_os("OGEOM_M2_CORPUS").is_some()
}

fn fine() -> ogeom_mesh::Deflection {
    ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    }
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, fine(), T).unwrap().mass
}

/// A vertex's point, placed.
fn point_of(m: &Model, v: &Shape) -> Point {
    let p = m
        .node(v)
        .and_then(|n| n.data().as_vertex().map(|d| d.point))
        .unwrap();
    v.transform(m.datums()).unwrap().apply(p)
}

/// Every vertex of `s`, as bits, sorted.
fn vertex_bits(m: &Model, s: &Shape) -> Vec<[u64; 3]> {
    let mut out: Vec<[u64; 3]> = explore_unique(m, s, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| {
            let p = point_of(m, v);
            [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
        })
        .collect();
    out.sort_unstable();
    out
}

/// A prism: the profile, given as (y, z), in the plane x = 0, swept `length`
/// along +x, stated in place (its end no copy of its start moved along).
fn prism(m: &mut Model, profile: &[(f64, f64)], length: f64) -> Shape {
    let points: Vec<Point> = profile
        .iter()
        .map(|&(y, z)| Point::new(0.0, y, z))
        .collect();
    let wire = ogeom_algo::make_polygon(m, &points, true, T).unwrap().shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::X),
        (-1000.0, 1000.0),
        (-1000.0, 1000.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(m, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let swept = ogeom_algo::make_prism(m, &face, Vector::new(length, 0.0, 0.0), T)
        .unwrap()
        .shape;
    ogeom_algo::baked_shape(m, &swept, T).unwrap().shape
}

/// The edges of `s`, and those along the sweep (both ends `length` apart
/// along x).
fn edges(m: &Model, s: &Shape, length: f64) -> (Vec<Shape>, Vec<Shape>) {
    let all = explore_unique(m, s, ShapeType::Edge).unwrap();
    let along = all
        .iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
            let (p, q) = (point_of(m, &a), point_of(m, &b));
            ((p.x - q.x).abs() - length).abs() < 1e-9
        })
        .cloned()
        .collect();
    (all, along)
}

/// The area a round of radius `r` takes off a convex square corner, or
/// adds to a concave one.
fn square_corner(r: f64) -> f64 {
    (1.0 - PI / 4.0) * r * r
}

/// The area of a simple polygon given as (y, z).
fn polygon_area(profile: &[(f64, f64)]) -> f64 {
    let n = profile.len();
    0.5 * (0..n)
        .map(|i| {
            let (a, b) = (profile[i], profile[(i + 1) % n]);
            a.0 * b.1 - a.1 * b.0
        })
        .sum::<f64>()
        .abs()
}

/// The convex and the concave corners of a rectilinear profile listed
/// anticlockwise: (convex, concave).
fn corners(profile: &[(f64, f64)]) -> (usize, usize) {
    let n = profile.len();
    let mut convex = 0;
    let mut concave = 0;
    for i in 0..n {
        let (p, q, r) = (profile[(i + n - 1) % n], profile[i], profile[(i + 1) % n]);
        let turn = (q.0 - p.0) * (r.1 - q.1) - (q.1 - p.1) * (r.0 - q.0);
        if turn > 0.0 {
            convex += 1;
        } else {
            concave += 1;
        }
    }
    (convex, concave)
}

/// Orders of a request of `n` edges: as given, reversed, and a fixed
/// shuffle.
fn orders(n: usize) -> Vec<Vec<usize>> {
    let given: Vec<usize> = (0..n).collect();
    let reversed: Vec<usize> = (0..n).rev().collect();
    let shuffled: Vec<usize> = (0..n).map(|i| (i * 7 + 3) % n).collect();
    let mut out = vec![given, reversed];
    let mut seen = vec![false; n];
    if shuffled
        .iter()
        .all(|&i| !std::mem::replace(&mut seen[i], true))
    {
        out.push(shuffled);
    }
    out
}

/// The guards every case is held to: validity, the region guard (each
/// vertex of the result inside the original solid, or within `reach` of
/// one of its concave edges, given as segments) and no overhang (each
/// vertex of the result within the original's box).
fn guards(m: &Model, original: &Shape, result: &Shape, concave: &[(Point, Point)], reach: f64) {
    let diagnosis = ogeom_algo::check(m, result, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let bounds = ogeom_algo::shape_bounds(m, original, T)
        .unwrap()
        .expanded(1e-7);
    let probe = ogeom_algo::SolidBoundary::of(m, original, 1e-6, T).unwrap();
    for v in explore_unique(m, result, ShapeType::Vertex).unwrap() {
        let p = point_of(m, &v);
        assert!(
            bounds.contains(p),
            "a vertex at {p:?} overhangs the original's box"
        );
        if probe.holds(m, p, T).unwrap() != ogeom_algo::Containment::Out {
            continue;
        }
        let near = concave.iter().any(|&(a, b)| {
            let d = b - a;
            let t = ((p - a).dot(d) / d.dot(d)).clamp(0.0, 1.0);
            (a + d * t).distance(p) <= reach * 1.5 + 1e-6
        });
        assert!(
            near,
            "a vertex at {p:?} stands outside the original, away from its concave edges"
        );
    }
}

/// Rounds `picked` of `s` at `radius` in every order of [`orders`], checks
/// each against the guards and the others to the bit; the volume.
fn rounded_every_order(
    m: &mut Model,
    s: &Shape,
    picked: &[Shape],
    radius: f64,
    concave: &[(Point, Point)],
) -> f64 {
    let mut first: Option<(u64, Vec<[u64; 3]>)> = None;
    let mut v0 = 0.0;
    for order in orders(picked.len()) {
        let asked: Vec<Shape> = order.iter().map(|&i| picked[i].clone()).collect();
        let built = ogeom_fillet::fillet_edges(m, s, &asked, radius, T).unwrap();
        guards(m, s, &built.shape, concave, radius);
        let v = volume(m, &built.shape);
        let key = (v.to_bits(), vertex_bits(m, &built.shape));
        match &first {
            None => {
                first = Some(key);
                v0 = v;
            }
            Some(k) => assert_eq!(*k, key, "the request's order changes the solid"),
        }
    }
    v0
}

/// The concave edges of a prism along its sweep, at the profile's concave
/// corners.
fn concave_edges(profile: &[(f64, f64)], length: f64) -> Vec<(Point, Point)> {
    let n = profile.len();
    let mut out = Vec::new();
    for i in 0..n {
        let (p, q, r) = (profile[(i + n - 1) % n], profile[i], profile[(i + 1) % n]);
        let turn = (q.0 - p.0) * (r.1 - q.1) - (q.1 - p.1) * (r.0 - q.0);
        if turn < 0.0 {
            out.push((Point::new(0.0, q.0, q.1), Point::new(length, q.0, q.1)));
        }
    }
    out
}

/// Every edge of an a x b x c box at radius r: the inner box swept by the
/// ball, its volume `(a')(b')(c') + 2r(a'b' + b'c' + c'a') + pi r^2 (a' +
/// b' + c') + 4/3 pi r^3`, the primes the box's sides less 2r.
#[test]
#[ignore = "heavy"]
fn every_edge_of_a_box_is_the_inner_box_swept_by_the_ball() {
    if !asked() {
        return;
    }
    let (a, b, c) = (20.0, 14.0, 9.0);
    for r in [0.5, 1.0, 2.0, 4.0] {
        let mut m = Model::new();
        let s = ogeom_algo::make_box(&mut m, Frame::WORLD, (a, b, c), T)
            .unwrap()
            .shape;
        let all = explore_unique(&m, &s, ShapeType::Edge).unwrap();
        let v = rounded_every_order(&mut m, &s, &all, r, &[]);
        let (p, q, w) = (a - 2.0 * r, b - 2.0 * r, c - 2.0 * r);
        let exact = p * q * w
            + 2.0 * r * (p * q + q * w + w * p)
            + PI * r * r * (p + q + w)
            + 4.0 / 3.0 * PI * r.powi(3);
        assert!(
            (v - exact).abs() <= 1e-8 * exact,
            "r {r}: {v} against {exact}"
        );
    }
}

/// A prismatic part: its name, profile, sweep and the radii that fit it.
type Part = (&'static str, Vec<(f64, f64)>, f64, Vec<f64>);

/// The profiles of ADR 0001 §12's prismatic parts, given as (y, z)
/// anticlockwise, their sweep, and the radii that fit them.
fn prismatic() -> Vec<Part> {
    vec![
        (
            "L channel",
            vec![
                (0.0, 0.0),
                (20.0, 0.0),
                (20.0, 5.0),
                (5.0, 5.0),
                (5.0, 20.0),
                (0.0, 20.0),
            ],
            30.0,
            vec![0.5, 1.0, 2.0],
        ),
        (
            "T channel",
            vec![
                (0.0, 0.0),
                (30.0, 0.0),
                (30.0, 5.0),
                (17.5, 5.0),
                (17.5, 20.0),
                (12.5, 20.0),
                (12.5, 5.0),
                (0.0, 5.0),
            ],
            30.0,
            vec![0.5, 1.0, 2.0],
        ),
        (
            "U channel",
            vec![
                (0.0, 0.0),
                (30.0, 0.0),
                (30.0, 20.0),
                (25.0, 20.0),
                (25.0, 5.0),
                (5.0, 5.0),
                (5.0, 20.0),
                (0.0, 20.0),
            ],
            30.0,
            vec![0.5, 1.0, 2.0],
        ),
        (
            "stepped block",
            vec![
                (0.0, 0.0),
                (30.0, 0.0),
                (30.0, 10.0),
                (15.0, 10.0),
                (15.0, 20.0),
                (0.0, 20.0),
            ],
            25.0,
            vec![1.0, 2.0, 4.0],
        ),
        (
            "the #902 bracket",
            vec![
                (0.0, 0.0),
                (60.0, 0.0),
                (60.0, 6.0),
                (30.0, 6.0),
                (30.0, 16.0),
                (15.0, 16.0),
                (15.0, 30.0),
                (0.0, 30.0),
            ],
            80.0,
            vec![1.0, 2.0],
        ),
    ]
}

/// The edges along each prismatic part's sweep, rounded: the rounded
/// profile swept along, its area the profile's less the square corner at
/// each convex corner, plus it at each concave one.
#[test]
#[ignore = "heavy"]
fn the_sweep_edges_of_each_prismatic_part_take_their_rounded_profile_along() {
    if !asked() {
        return;
    }
    for (name, profile, length, radii) in prismatic() {
        let (convex, concave) = corners(&profile);
        let reaches = concave_edges(&profile, length);
        for r in radii {
            let mut m = Model::new();
            let s = prism(&mut m, &profile, length);
            let (_, along) = edges(&m, &s, length);
            assert_eq!(along.len(), profile.len(), "{name}: its sweep edges");
            let v = rounded_every_order(&mut m, &s, &along, r, &reaches);
            #[allow(clippy::cast_precision_loss, reason = "a few corners")]
            let area = polygon_area(&profile) - square_corner(r) * (convex as f64)
                + square_corner(r) * (concave as f64);
            let exact = area * length;
            assert!(
                (v - exact).abs() <= 1e-8 * exact,
                "{name}, r {r}: {v} against {exact}"
            );
        }
    }
}

/// Every edge of each prismatic part, its ends' edges included: spheres at
/// the convex corners of three edges, a patch at the mixed corners where
/// a concave sweep edge meets its end. No closed form; the solid takes
/// less than its convex strips over their whole edges and a corner's
/// worth per vertex, and adds no more than its concave strips do.
#[test]
#[ignore = "heavy"]
fn every_edge_of_each_prismatic_part_rounds_whatever_the_order() {
    if !asked() {
        return;
    }
    for (name, profile, length, radii) in prismatic() {
        let (convex, concave) = corners(&profile);
        let reaches = concave_edges(&profile, length);
        let perimeter: f64 = (0..profile.len())
            .map(|i| {
                let (a, b) = (profile[i], profile[(i + 1) % profile.len()]);
                (b.0 - a.0).hypot(b.1 - a.1)
            })
            .sum();
        for r in radii {
            let mut m = Model::new();
            let s = prism(&mut m, &profile, length);
            let before = volume(&m, &s);
            let (all, _) = edges(&m, &s, length);
            let v = rounded_every_order(&mut m, &s, &all, r, &reaches);
            let a = square_corner(r);
            #[allow(clippy::cast_precision_loss, reason = "a few corners")]
            let (cx, cc, vertices) = (convex as f64, concave as f64, (2 * profile.len()) as f64);
            let taken = a * (cx * length + 2.0 * perimeter) + vertices * r.powi(3);
            let added = a * cc * length;
            assert!(
                v > before - taken && v < before + added,
                "{name}, r {r}: {v} outside ({}, {})",
                before - taken,
                before + added
            );
        }
    }
}

/// Rounds and bevels in one request: the bevels first, as one set, then
/// the rounds on what they leave (ADR 0001 §13, question 3). The box's top
/// edges bevelled and its uprights rounded: whatever the order, one solid,
/// taking more than the bevels alone and no more than the bevels and the
/// rounds over their whole edges.
#[test]
#[ignore = "heavy"]
fn rounds_and_bevels_in_one_request_take_the_bevels_first_whatever_the_order() {
    if !asked() {
        return;
    }
    let (a, b, c) = (20.0, 14.0, 9.0);
    let (d, r) = (1.0, 2.0);
    let mut first: Option<(u64, Vec<[u64; 3]>)> = None;
    let mut m = Model::new();
    let s = ogeom_algo::make_box(&mut m, Frame::WORLD, (a, b, c), T)
        .unwrap()
        .shape;
    let all = explore_unique(&m, &s, ShapeType::Edge).unwrap();
    let kind = |e: &Shape| -> Blend {
        let (p, q) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
        let (p, q) = (point_of(&m, &p), point_of(&m, &q));
        if (p.z - q.z).abs() > 1e-9 {
            Blend::Round(r)
        } else if (p.z - c).abs() < 1e-9 {
            Blend::Chamfer(Chamfer::Symmetric(d))
        } else {
            Blend::Round(r)
        }
    };
    let request: Vec<(Shape, Blend)> = all.iter().map(|e| (e.clone(), kind(e))).collect();
    let bevel_volume = 0.5 * d * d * 2.0 * (a + b);
    for order in orders(request.len()) {
        let asked: Vec<(Shape, Blend)> = order.iter().map(|&i| request[i].clone()).collect();
        let built = ogeom_fillet::blend_edges(&mut m, &s, &asked, ChamferCorner::Facet, T).unwrap();
        guards(&m, &s, &built.shape, &[], r);
        let v = volume(&m, &built.shape);
        let rounds_whole = square_corner(r) * (4.0 * c + 2.0 * (a + b));
        assert!(
            v < a * b * c - bevel_volume
                && v > a * b * c - bevel_volume - rounds_whole - 8.0 * r.powi(3),
            "{v}"
        );
        let key = (v.to_bits(), vertex_bits(&m, &built.shape));
        match &first {
            None => first = Some(key),
            Some(k) => assert_eq!(*k, key, "the request's order changes the solid"),
        }
    }
}

/// Unequal radii at a box's corners: no ball touches all three faces, so
/// each corner is the general patch. The solid is the same in every order,
/// and takes between the strips' sections over their edges short of the
/// corners and over the whole edges with a corner's worth per vertex.
#[test]
#[ignore = "heavy"]
fn unequal_radii_at_a_box_s_corners_patch_whatever_the_order() {
    if !asked() {
        return;
    }
    let (a, b, c) = (20.0, 14.0, 9.0);
    let mut m = Model::new();
    let s = ogeom_algo::make_box(&mut m, Frame::WORLD, (a, b, c), T)
        .unwrap()
        .shape;
    let all = explore_unique(&m, &s, ShapeType::Edge).unwrap();
    // Each edge's radius by its direction, and its length.
    let sized: Vec<(Shape, f64, f64)> = all
        .iter()
        .map(|e| {
            let (p, q) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
            let d = point_of(&m, &q) - point_of(&m, &p);
            let (r, l) = if d.x.abs() > 1e-9 {
                (1.0, d.x.abs())
            } else if d.y.abs() > 1e-9 {
                (1.5, d.y.abs())
            } else {
                (2.0, d.z.abs())
            };
            (e.clone(), r, l)
        })
        .collect();
    let request: Vec<(Shape, f64)> = sized.iter().map(|(e, r, _)| (e.clone(), *r)).collect();
    let whole: f64 = sized.iter().map(|(_, r, l)| square_corner(*r) * l).sum();
    let short: f64 = sized
        .iter()
        .map(|(_, r, l)| square_corner(*r) * (l - 4.0))
        .sum();
    let mut first: Option<(u64, Vec<[u64; 3]>)> = None;
    for order in orders(request.len()) {
        let asked: Vec<(Shape, f64)> = order.iter().map(|&i| request[i].clone()).collect();
        let built = ogeom_fillet::fillet_edges_with(&mut m, &s, &asked, T).unwrap();
        guards(&m, &s, &built.shape, &[], 2.0);
        let v = volume(&m, &built.shape);
        assert!(
            v < a * b * c - short && v > a * b * c - whole - 8.0 * 8.0,
            "{v}"
        );
        let key = (v.to_bits(), vertex_bits(&m, &built.shape));
        match &first {
            None => first = Some(key),
            Some(k) => assert_eq!(*k, key, "the request's order changes the solid"),
        }
    }
}

/// ADR 0001 §7: every edge of a 50-edge prismatic part rounded in under a
/// second, release build. A comb of four teeth on a base: 20 sweep edges
/// and 40 on its ends.
#[test]
#[ignore = "heavy"]
fn every_edge_of_a_sixty_edge_comb_rounds_in_under_a_second() {
    if !asked() {
        return;
    }
    let mut profile = vec![(0.0, 0.0), (18.0, 0.0), (18.0, 6.0)];
    for k in (0..4).rev() {
        #[allow(clippy::cast_precision_loss, reason = "four teeth")]
        let y = 4.0 * k as f64;
        profile.extend([
            (y + 3.0, 6.0),
            (y + 3.0, 12.0),
            (y + 1.0, 12.0),
            (y + 1.0, 6.0),
        ]);
    }
    profile.push((0.0, 6.0));
    let mut m = Model::new();
    let s = prism(&mut m, &profile, 20.0);
    let (all, _) = edges(&m, &s, 20.0);
    assert_eq!(all.len(), 60, "the comb's edges");
    let started = ogeom_core::clock::Instant::now();
    let built = ogeom_fillet::fillet_edges(&mut m, &s, &all, 0.4, T).unwrap();
    let took = started.elapsed();
    let reaches = concave_edges(&profile, 20.0);
    guards(&m, &s, &built.shape, &reaches, 0.4);
    if !cfg!(debug_assertions) {
        assert!(took.as_secs_f64() < 1.0, "every edge took {took:?}");
    }
}
