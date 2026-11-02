//! Shadows of meshes on planes, and outlines fitted around them in a band.
//!
//! Every outline is checked against the shadow it was fitted to: each
//! piece sampled at a hundredth of the tolerance, every sample's signed
//! distance from the shadow, found by brute force over all its edges,
//! between the play and the play plus twice the tolerance.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    reason = "test code"
)]

#[path = "support/cpu_time.rs"]
mod cpu_time;
use cpu_time::cpu_time;

use core::f64::consts::TAU;
use std::time::Duration;

use ogeom::algo::{Outline, Piece2, Wire2d, fit_outline};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Plane, Point, Point2, Vector};
use ogeom::mesh::{Loop2, shadow};
use ogeom::topo::Triangulation;

/// A torus about the `z` axis, `rings` by `sides` quads, two triangles each.
fn torus(major: f64, minor: f64, rings: u32, sides: u32) -> Triangulation {
    let mut mesh = Triangulation::new();
    for i in 0..rings {
        let u = f64::from(i) / f64::from(rings) * TAU;
        for j in 0..sides {
            let v = f64::from(j) / f64::from(sides) * TAU;
            let r = major + minor * v.cos();
            mesh.positions
                .push(Point::new(r * u.cos(), r * u.sin(), minor * v.sin()));
        }
    }
    let at = |i: u32, j: u32| (i % rings) * sides + (j % sides);
    for i in 0..rings {
        for j in 0..sides {
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            mesh.triangles.push([a, b, c]);
            mesh.triangles.push([a, c, d]);
        }
    }
    mesh
}

/// A prism over a simple polygon given counter-clockwise, `height` tall
/// along `z`: the polygon's ear-clipped caps and a quad per side.
fn prism(outline: &[(f64, f64)], height: f64) -> Triangulation {
    let n = outline.len() as u32;
    let mut mesh = Triangulation::new();
    for z in [0.0, height] {
        for &(x, y) in outline {
            mesh.positions.push(Point::new(x, y, z));
        }
    }
    let caps = ear_clip(outline);
    for &[a, b, c] in &caps {
        mesh.triangles.push([a, c, b]);
        mesh.triangles.push([a + n, b + n, c + n]);
    }
    for i in 0..n {
        let j = (i + 1) % n;
        mesh.triangles.push([i, j, j + n]);
        mesh.triangles.push([i, j + n, i + n]);
    }
    mesh
}

fn ear_clip(polygon: &[(f64, f64)]) -> Vec<[u32; 3]> {
    let mut left: Vec<u32> = (0..polygon.len() as u32).collect();
    let at = |k: u32| polygon[k as usize];
    let cross = |a: (f64, f64), b: (f64, f64), c: (f64, f64)| {
        (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
    };
    let mut out = Vec::new();
    while left.len() > 3 {
        let n = left.len();
        let ear = (0..n)
            .find(|&i| {
                let (a, b, c) = (left[(i + n - 1) % n], left[i], left[(i + 1) % n]);
                cross(at(a), at(b), at(c)) > 0.0
                    && left.iter().all(|&k| {
                        k == a
                            || k == b
                            || k == c
                            || cross(at(a), at(b), at(k)) < 0.0
                            || cross(at(b), at(c), at(k)) < 0.0
                            || cross(at(c), at(a), at(k)) < 0.0
                    })
            })
            .unwrap();
        out.push([left[(ear + n - 1) % n], left[ear], left[(ear + 1) % n]]);
        left.remove(ear);
    }
    out.push([left[0], left[1], left[2]]);
    out
}

/// A closed cylinder about the `z` axis, `sides` around and `rows` up.
fn cylinder(radius: f64, height: f64, sides: u32, rows: u32) -> Triangulation {
    let mut mesh = Triangulation::new();
    for j in 0..=rows {
        let z = height * f64::from(j) / f64::from(rows);
        for i in 0..sides {
            let a = f64::from(i) / f64::from(sides) * TAU;
            mesh.positions
                .push(Point::new(radius * a.cos(), radius * a.sin(), z));
        }
    }
    let at = |i: u32, j: u32| j * sides + i % sides;
    for j in 0..rows {
        for i in 0..sides {
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            mesh.triangles.push([a, b, c]);
            mesh.triangles.push([a, c, d]);
        }
    }
    let (bottom, top) = (mesh.positions.len() as u32, mesh.positions.len() as u32 + 1);
    mesh.positions.push(Point::new(0.0, 0.0, 0.0));
    mesh.positions.push(Point::new(0.0, 0.0, height));
    for i in 0..sides {
        mesh.triangles.push([bottom, at(i + 1, 0), at(i, 0)]);
        mesh.triangles.push([top, at(i, rows), at(i + 1, rows)]);
    }
    mesh
}

/// A square plate `side` wide and `thick` thick with a round hole of
/// `radius` through its middle, `4 k` points round the hole.
fn holed_plate(side: f64, thick: f64, radius: f64, k: u32) -> Triangulation {
    let n = 4 * k;
    let mut mesh = Triangulation::new();
    for z in [0.0, thick] {
        for i in 0..n {
            let a = f64::from(i) / f64::from(n) * TAU;
            let (c, s) = (a.cos(), a.sin());
            mesh.positions.push(Point::new(radius * c, radius * s, z));
            // Out along the same ray to the square's edge.
            let t = side / 2.0 / c.abs().max(s.abs());
            mesh.positions.push(Point::new(t * c, t * s, z));
        }
    }
    let hole = |i: u32, top: u32| 2 * (i % n) + top * 2 * n;
    let rim = |i: u32, top: u32| 2 * (i % n) + 1 + top * 2 * n;
    for i in 0..n {
        // Top faces up, bottom down, the rim's wall out, the hole's in.
        mesh.triangles.push([hole(i, 1), rim(i, 1), rim(i + 1, 1)]);
        mesh.triangles
            .push([hole(i, 1), rim(i + 1, 1), hole(i + 1, 1)]);
        mesh.triangles.push([hole(i, 0), rim(i + 1, 0), rim(i, 0)]);
        mesh.triangles
            .push([hole(i, 0), hole(i + 1, 0), rim(i + 1, 0)]);
        mesh.triangles
            .push([rim(i, 0), rim(i + 1, 0), rim(i + 1, 1)]);
        mesh.triangles.push([rim(i, 0), rim(i + 1, 1), rim(i, 1)]);
        mesh.triangles
            .push([hole(i, 0), hole(i, 1), hole(i + 1, 1)]);
        mesh.triangles
            .push([hole(i, 0), hole(i + 1, 1), hole(i + 1, 0)]);
    }
    mesh
}

/// Signed distance from the region the loops bound, by brute force:
/// positive outside.
fn signed_distance(loops: &[Loop2], p: Point2) -> f64 {
    let mut nearest = f64::INFINITY;
    let mut winding = 0;
    for l in loops {
        let n = l.points.len();
        for i in 0..n {
            let (a, b) = (l.points[i], l.points[(i + 1) % n]);
            let d = b - a;
            let t = ((p - a).dot(d) / d.dot(d)).clamp(0.0, 1.0);
            nearest = nearest.min((a + d * t).distance(p));
            let side = d.cross(p - a);
            if a.y <= p.y && p.y < b.y && side > 0.0 {
                winding += 1;
            } else if b.y <= p.y && p.y < a.y && side < 0.0 {
                winding -= 1;
            }
        }
    }
    if winding > 0 { -nearest } else { nearest }
}

/// Every point of every wire, sampled at a hundredth of the tolerance,
/// lies in the band; the wires close up.
fn assert_in_band(loops: &[Loop2], wires: &[Wire2d], o: &Outline) {
    let step = o.tolerance / 100.0;
    for w in wires {
        let n = w.pieces.len();
        for (k, p) in w.pieces.iter().enumerate() {
            let next = w.pieces[(k + 1) % n].start();
            assert!(p.end().distance(next) < 1e-9, "a gap after piece {k}");
            let samples = (p.length() / step).ceil().max(1.0) as usize;
            for s in 0..=samples {
                let q = p.point_at(s as f64 / samples as f64);
                let d = signed_distance(loops, q);
                assert!(
                    o.play <= d && d <= o.play + 2.0 * o.tolerance,
                    "{q:?} at {d}, outside the band [{}, {}]",
                    o.play,
                    o.play + 2.0 * o.tolerance
                );
            }
        }
    }
}

fn along_z() -> Plane {
    Plane::XY
}

/// A cylinder seen along its axis casts its polygon, and the outline is
/// one circle, of radius the cylinder's plus the play and within the
/// tolerance.
#[test]
fn a_cylinder_seen_along_its_axis_gives_one_circle() {
    let (r, sides) = (10.0, 256);
    let mesh = cylinder(r, 20.0, sides, 8);
    let loops = shadow(&mesh, &along_z()).unwrap();
    assert_eq!(loops.len(), 1);
    assert_eq!(loops[0].points.len(), sides as usize);
    let o = Outline {
        play: 0.2,
        tolerance: 0.05,
        round_below: 0.5,
    };
    let wires = fit_outline(&loops, &o).unwrap();
    assert_eq!(wires.len(), 1);
    let [
        Piece2::Circle {
            centre,
            radius,
            counter_clockwise,
        },
    ] = wires[0].pieces[..]
    else {
        panic!("{:?}", wires[0]);
    };
    assert!(counter_clockwise);
    assert!(centre.distance(Point2::ORIGIN) < 1e-3);
    // The polygon's sides sit inside the circle through its corners.
    let inside = r * (core::f64::consts::PI / f64::from(sides)).cos();
    assert!(
        inside + o.play <= radius && radius <= r + o.play + 2.0 * o.tolerance,
        "{radius}"
    );
    assert_in_band(&loops, &wires, &o);
}

const L: [(f64, f64); 6] = [
    (0.0, 0.0),
    (40.0, 0.0),
    (40.0, 5.0),
    (5.0, 5.0),
    (5.0, 30.0),
    (0.0, 30.0),
];

/// An L bracket seen end on gives its six sides as six lines meeting at
/// corners, each corner where the sides moved out to the band's middle
/// meet, within the tolerance.
#[test]
fn an_l_bracket_gives_its_outline_with_corners() {
    let mesh = prism(&L, 12.0);
    let loops = shadow(&mesh, &along_z()).unwrap();
    assert_eq!(loops.len(), 1);
    assert_eq!(loops[0].points.len(), 6);
    assert!((loops[0].area() - (40.0 * 5.0 + 5.0 * 25.0)).abs() < 1e-9);
    let o = Outline {
        play: 0.05,
        tolerance: 0.05,
        round_below: 0.5,
    };
    let wires = fit_outline(&loops, &o).unwrap();
    assert_eq!(wires.len(), 1);
    let pieces = &wires[0].pieces;
    assert_eq!(pieces.len(), 6, "{pieces:?}");
    assert!(pieces.iter().all(|p| matches!(p, Piece2::Line { .. })));
    let m = o.play + o.tolerance;
    let expected = [
        (-m, -m),
        (40.0 + m, -m),
        (40.0 + m, 5.0 + m),
        (5.0 + m, 5.0 + m),
        (5.0 + m, 30.0 + m),
        (-m, 30.0 + m),
    ];
    for &(x, y) in &expected {
        let near = pieces
            .iter()
            .map(|p| p.start().distance(Point2::new(x, y)))
            .fold(f64::INFINITY, f64::min);
        assert!(near < o.tolerance, "no corner near ({x}, {y}): {pieces:?}");
    }
    assert_in_band(&loops, &wires, &o);
}

/// A press fit: with a negative play the outline runs inside the
/// bracket. The five convex corners stay sharp; at the reflex one the
/// sides moved in meet too deep, and an arc rounds it.
#[test]
fn a_negative_play_fits_inside() {
    let mesh = prism(&L, 12.0);
    let loops = shadow(&mesh, &along_z()).unwrap();
    let o = Outline {
        play: -0.2,
        tolerance: 0.05,
        round_below: 0.5,
    };
    let wires = fit_outline(&loops, &o).unwrap();
    assert_eq!(wires.len(), 1);
    let pieces = &wires[0].pieces;
    let lines = pieces
        .iter()
        .filter(|p| matches!(p, Piece2::Line { .. }))
        .count();
    let arcs: Vec<f64> = pieces
        .iter()
        .filter_map(|p| match p {
            Piece2::Arc { radius, .. } => Some(*radius),
            _ => None,
        })
        .collect();
    assert_eq!((lines, arcs.len()), (6, 1), "{pieces:?}");
    // At least the depth of the band's middle, the exact inner offset's
    // own rounding there.
    assert!(arcs[0] >= -(o.play + o.tolerance), "{arcs:?}");
    assert_in_band(&loops, &wires, &o);
}

/// A plate with a hole through it keeps the hole in its shadow, and the
/// outline is the square's four corners and the hole's circle, run
/// clockwise and shrunk by the play.
#[test]
fn a_through_hole_is_kept() {
    let (side, radius) = (40.0, 6.0);
    let mesh = holed_plate(side, 4.0, radius, 32);
    let loops = shadow(&mesh, &along_z()).unwrap();
    assert_eq!(loops.len(), 2);
    let holes: Vec<&Loop2> = loops.iter().filter(|l| l.is_hole()).collect();
    assert_eq!(holes.len(), 1);
    assert!((holes[0].area() + core::f64::consts::PI * radius * radius).abs() < 0.2);
    let o = Outline {
        play: 0.05,
        tolerance: 0.05,
        round_below: 0.5,
    };
    let wires = fit_outline(&loops, &o).unwrap();
    assert_eq!(wires.len(), 2);
    let circle = wires
        .iter()
        .find_map(|w| match w.pieces[..] {
            [
                Piece2::Circle {
                    radius: r,
                    counter_clockwise,
                    ..
                },
            ] => Some((r, counter_clockwise)),
            _ => None,
        })
        .expect("the hole is a circle");
    assert!(!circle.1);
    assert!(
        radius - o.play - 2.0 * o.tolerance <= circle.0 && circle.0 <= radius - o.play,
        "{circle:?}"
    );
    let square = wires
        .iter()
        .find(|w| w.pieces.len() == 4)
        .expect("four sides");
    assert!(
        square
            .pieces
            .iter()
            .all(|p| matches!(p, Piece2::Line { .. }))
    );
    assert_in_band(&loops, &wires, &o);
}

/// Gaussian noise from a fixed seed.
struct Noise(u64);

impl Noise {
    fn uniform(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn gauss(&mut self, sigma: f64) -> f64 {
        let (u, v) = (self.uniform().max(1e-300), self.uniform());
        sigma * (-2.0 * u.ln()).sqrt() * (TAU * v).cos()
    }
}

/// A scanned cylinder, every vertex moved by noise of `sigma = 0.02`:
/// the shadow is jagged, and the outline is still one circle in the band
/// about it.
#[test]
fn a_noisy_cylinder_still_gives_one_circle() {
    let mut mesh = cylinder(10.0, 20.0, 512, 40);
    let mut noise = Noise(0x9e37_79b9_7f4a_7c15);
    for p in &mut mesh.positions {
        *p = Point::new(
            p.x + noise.gauss(0.02),
            p.y + noise.gauss(0.02),
            p.z + noise.gauss(0.02),
        );
    }
    let loops = shadow(&mesh, &along_z()).unwrap();
    assert_eq!(loops.iter().filter(|l| !l.is_hole()).count(), 1);
    let o = Outline {
        play: 0.1,
        tolerance: 0.1,
        round_below: 0.5,
    };
    let (wires, took) = cpu_time(|| fit_outline(&loops, &o).unwrap());
    eprintln!(
        "noisy cylinder: {} shadow points, fitted in {took:?}",
        loops.iter().map(|l| l.points.len()).sum::<usize>()
    );
    assert_eq!(wires.len(), 1);
    let [Piece2::Circle { radius, .. }] = wires[0].pieces[..] else {
        panic!("{:?}", wires[0]);
    };
    // Clear of the noise's outermost points, and not by much more.
    assert!(10.0 + o.play < radius && radius < 10.0 + o.play + 2.0 * o.tolerance + 0.1);
    assert_in_band(&loops, &wires, &o);
}

/// A bracket stored in single precision, far from the origin where the
/// rounding is coarse: the shadow is exact to the rounded points and the
/// outline keeps to the band about it.
#[test]
fn a_single_precision_mesh_far_from_the_origin() {
    let mut mesh = prism(&L, 12.0);
    for p in &mut mesh.positions {
        let moved = |v: f64, by: f64| f64::from((v + by) as f32);
        *p = Point::new(moved(p.x, 1234.567), moved(p.y, -987.654), moved(p.z, 55.5));
    }
    let loops = shadow(&mesh, &along_z()).unwrap();
    assert_eq!(loops.len(), 1);
    assert_eq!(loops[0].points.len(), 6);
    assert!((loops[0].area() - 325.0).abs() < 1e-2);
    let o = Outline {
        play: 0.05,
        tolerance: 0.05,
        round_below: 0.5,
    };
    let wires = fit_outline(&loops, &o).unwrap();
    assert_eq!(wires.len(), 1);
    assert_eq!(wires[0].pieces.len(), 6, "{:?}", wires[0]);
    assert_in_band(&loops, &wires, &o);
}

/// A million triangles, a torus seen along its axis and tilted: its shadow
/// is a ring, the hole kept, of the area the closed form gives, in well
/// under a second of processor time.
#[test]
#[ignore = "heavy"]
fn a_million_triangle_torus_casts_a_ring() {
    let mesh = torus(10.0, 3.0, 1000, 500);
    assert_eq!(mesh.triangles.len(), 1_000_000);
    let plane = Plane::through(Point::ORIGIN, Direction::Z);
    let (loops, took) = cpu_time(|| shadow(&mesh, &plane).unwrap());
    eprintln!("shadow of 1M triangles: {took:?}");
    assert_eq!(loops.len(), 2);
    let area: f64 = loops.iter().map(Loop2::area).sum();
    let ring = core::f64::consts::PI * (13.0f64.powi(2) - 7.0f64.powi(2));
    assert!((area - ring).abs() < 1e-3 * ring, "{area} against {ring}");
    // Under 0.2 s on the development machine; the bound leaves room for a
    // shared CI runner several times slower.
    assert!(took < Duration::from_secs(2), "{took:?} of processor time");

    let tilted = Plane::through(
        Point::ORIGIN,
        Direction::new(Vector::new(0.3, 0.03, 1.0), Tolerances::millimetres()).unwrap(),
    );
    let loops = shadow(&mesh, &tilted).unwrap();
    assert_eq!(loops.len(), 2);
    assert_eq!(loops.iter().filter(|l| l.is_hole()).count(), 1);

    // The same torus scanned: every vertex moved by noise of 0.02.
    let mut noisy = mesh;
    let mut noise = Noise(0x2545_f491_4f6c_dd1d);
    for p in &mut noisy.positions {
        *p = Point::new(
            p.x + noise.gauss(0.02),
            p.y + noise.gauss(0.02),
            p.z + noise.gauss(0.02),
        );
    }
    let (loops, took) = cpu_time(|| shadow(&noisy, &plane).unwrap());
    let o = Outline {
        play: 0.1,
        tolerance: 0.1,
        round_below: 0.5,
    };
    let (wires, fitted) = cpu_time(|| fit_outline(&loops, &o).unwrap());
    eprintln!(
        "noisy 1M triangles: shadow of {} points in {took:?}, outline of {:?} pieces in {fitted:?}",
        loops.iter().map(|l| l.points.len()).sum::<usize>(),
        wires.iter().map(|w| w.pieces.len()).collect::<Vec<_>>()
    );
}
