#![allow(clippy::unwrap_used, clippy::cast_precision_loss)]

use super::*;
use crate::Editable;
use crate::boolean::check_solid;
use std::f64::consts::{PI, TAU};

const T: Tolerances = Tolerances::millimetres();

/// An axis-aligned box from `lo` to `hi`, wound outward.
fn cuboid(lo: P3, hi: P3) -> Triangulation {
    let mut mesh = Triangulation::new();
    for k in 0..8_u32 {
        let pick = |b: u32, axis: usize| {
            if (k >> b) & 1 == 1 {
                hi[axis]
            } else {
                lo[axis]
            }
        };
        mesh.positions
            .push(Point::new(pick(0, 0), pick(1, 1), pick(2, 2)));
    }
    for [a, b, c, d] in [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
    ] {
        mesh.triangles.push([a, b, c]);
        mesh.triangles.push([a, c, d]);
    }
    mesh
}

/// A sphere about the origin of latitude rings and longitude slices, wound
/// outward.
fn sphere(radius: f64, rings: u32, slices: u32) -> Triangulation {
    let mut mesh = Triangulation::new();
    let at = |theta: f64, phi: f64| {
        Point::new(
            radius * theta.sin() * phi.cos(),
            radius * theta.sin() * phi.sin(),
            radius * theta.cos(),
        )
    };
    mesh.positions.push(at(0.0, 0.0));
    for r in 1..rings {
        let theta = PI * f64::from(r) / f64::from(rings);
        for s in 0..slices {
            mesh.positions
                .push(at(theta, TAU * f64::from(s) / f64::from(slices)));
        }
    }
    mesh.positions.push(at(PI, 0.0));
    let ring = |r: u32, s: u32| 1 + (r - 1) * slices + s % slices;
    let last = u32::try_from(mesh.positions.len() - 1).unwrap();
    for s in 0..slices {
        mesh.triangles.push([0, ring(1, s), ring(1, s + 1)]);
        for r in 1..rings - 1 {
            let (a, b) = (ring(r, s), ring(r, s + 1));
            let (c, d) = (ring(r + 1, s), ring(r + 1, s + 1));
            mesh.triangles.push([a, c, d]);
            mesh.triangles.push([a, d, b]);
        }
        mesh.triangles
            .push([ring(rings - 1, s), last, ring(rings - 1, s + 1)]);
    }
    mesh
}

/// The cavity's triangles alone, on the hollowed mesh's vertices.
fn cavity(done: &Hollowed) -> Triangulation {
    let mut out = done.mesh.clone();
    out.triangles.drain(..done.first_cavity_triangle);
    out
}

/// Each cavity vertex's distance to the part.
fn depths(part: &Triangulation, done: &Hollowed) -> Vec<f64> {
    let bvh = Bvh::build(part).unwrap();
    let base = part.positions.len();
    done.mesh.positions[base..]
        .iter()
        .map(|&p| bvh.nearest(p, f64::INFINITY).unwrap().distance)
        .collect()
}

/// The result is one closed 2-manifold mesh.
fn sealed(mesh: &Triangulation) {
    let e = Editable::from(mesh, 0.0).unwrap();
    assert!(e.is_closed(), "the hollowed mesh is open");
    assert!(
        e.non_manifold().is_empty(),
        "the hollowed mesh is not 2-manifold"
    );
    check_solid(mesh, "hollowed", 0).unwrap();
}

/// How far a cavity point may be from the wall: the interpolation error
/// `interpolation`, the clamp's `√3 h CLAMP`, and the decimation's
/// measured deviation.
fn slack(done: &Hollowed, interpolation: f64) -> f64 {
    let decimation = done.deviation.map_or(0.0, |d| d.furthest.1);
    3.0f64.sqrt().mul_add(done.cell * CLAMP, interpolation) + decimation
}

/// The hollowed volume of a cube of side `a` with a wall `w`, and the bound
/// on its error at cell `h`: about each of the inner cube's twelve edges,
/// where the field is the least of two faces' distances, the interpolation
/// falls short by at most `L / √2` within `L = √3 h` of the plane between
/// them, a band of `√2 L²` across; and the faces move by the clamp and the
/// decimation.
fn cube_volume(a: f64, w: f64, done: &Hollowed) -> (f64, f64) {
    let s = a - 2.0 * w;
    let l2 = 3.0 * done.cell * done.cell;
    let bound = (12.0 * s * 2.0f64.sqrt()).mul_add(l2, 6.0 * s * s * slack(done, 0.0));
    (a.powi(3) - s.powi(3), bound)
}

#[test]
fn a_hollowed_cube_keeps_its_wall_and_the_volume_of_a_box_less_a_box() {
    // The points of a cube of side `a` at least `w` from its faces are the
    // cube of side `a - 2w` about its centre, sharp edges and all.
    let (a, w) = (15.0, 2.5);
    let part = cuboid([0.0; 3], [a; 3]);
    let done = hollow(&part, w, &Grid::default(), T).unwrap();
    sealed(&done.mesh);
    assert_eq!(done.shells, 1);
    assert_eq!(done.mesh.positions[..8], part.positions[..]);
    assert_eq!(done.mesh.triangles[..12], part.triangles[..]);
    // The field is 1-Lipschitz, so its interpolation over a tetrahedron
    // errs by at most the tetrahedron's diameter.
    let e = slack(&done, 3.0f64.sqrt() * done.cell);
    for d in depths(&part, &done) {
        assert!((d - w).abs() < e, "a cavity vertex {d} deep");
    }
    let (want, bound) = cube_volume(a, w, &done);
    let got = done.mesh.volume();
    assert!(
        (got - want).abs() <= bound,
        "volume {got} against {want}, bound {bound}"
    );
}

#[test]
#[ignore = "heavy"]
fn a_hollowed_cube_converges_as_the_square_of_the_cell() {
    // The volume lost along the inner cube's sharp edges goes as `h²`:
    // halving the cell quarters it, near enough.
    let (a, w) = (30.0, 2.5);
    let part = cuboid([0.0; 3], [a; 3]);
    let errors: Vec<f64> = [w / 3.0, w / 6.0]
        .iter()
        .map(|&h| {
            let done = hollow(&part, w, &Grid { cell: Some(h) }, T).unwrap();
            sealed(&done.mesh);
            let (want, bound) = cube_volume(a, w, &done);
            let got = done.mesh.volume();
            assert!(
                (got - want).abs() <= bound,
                "volume {got} against {want}, bound {bound}"
            );
            got - want
        })
        .collect();
    let ratio = errors[0] / errors[1];
    assert!((3.0..5.0).contains(&ratio), "errors {errors:?}");
}

#[test]
fn a_hollowed_sphere_has_the_inner_radius() {
    let (r, w) = (6.0, 1.0);
    let part = sphere(r, 48, 96);
    let done = hollow(&part, w, &Grid::default(), T).unwrap();
    sealed(&done.mesh);
    assert_eq!(done.shells, 1);
    assert!(
        done.deviation.is_some(),
        "the decimated cavity was refused and the marched one kept"
    );
    let h = done.cell;
    // The facets lie between the sphere and `sag` inside it, so a point
    // `ρ` from the centre is between `R - sag - ρ` and `R - ρ` from them.
    let bvh = Bvh::build(&part).unwrap();
    let sag = r - bvh.nearest(Point::ORIGIN, f64::INFINITY).unwrap().distance;
    // Where the field is smooth, interpolation over a tetrahedron of
    // diameter `√3 h` errs by at most half its square times the curvature,
    // `1 / ρ`; the facets add `2 sag` either way.
    let inner = r - w;
    let rho = inner - sag - h;
    let e = slack(&done, 3.0 * h * h / (2.0 * rho) + 2.0 * sag);
    let base = part.positions.len();
    for p in &done.mesh.positions[base..] {
        let rho = p.to_vector().magnitude();
        assert!(
            rho > inner - sag - e && rho < inner + e,
            "a cavity vertex {rho} from the centre, {inner} wanted within {e} and the sag {sag}"
        );
    }
    let want = 4.0 / 3.0 * PI * inner.powi(3);
    let got = -cavity(&done).volume();
    let bound = 4.0 * PI * inner * inner * (e + sag);
    assert!((got - want).abs() <= bound, "cavity {got} against {want}");
}

#[test]
#[ignore = "heavy"]
fn a_ball_of_a_million_triangles_hollows_in_time() {
    // A ball 100 across with a 2 wall, of 998 000 triangles: the cavity
    // within 0.2 % of the closed form.
    let (r, w) = (50.0, 2.0);
    let part = sphere(r, 500, 1000);
    let start = ogeom_core::clock::Instant::now();
    let done = hollow(&part, w, &Grid::default(), T).unwrap();
    let took = start.elapsed();
    assert_eq!(done.shells, 1);
    let got = -cavity(&done).volume();
    let want = 4.0 / 3.0 * PI * (r - w).powi(3);
    eprintln!(
        "ball: {} triangles hollowed in {took:?}, cavity {} triangles, volume off by {:.3e}, {:?}",
        part.triangle_count(),
        done.mesh.triangle_count() - done.first_cavity_triangle,
        (got - want) / want,
        done.deviation
    );
    assert!(
        ((got - want) / want).abs() < 0.002,
        "cavity {got} against {want}"
    );
    // About 5 s on eight threads of the machine it was measured on.
    assert!(took < std::time::Duration::from_secs(50), "took {took:?}");
}

#[test]
fn a_part_with_a_cavity_keeps_the_wall_clear_of_it() {
    // A cube of 20 with a cube of 4 cut out of its middle. The points at
    // least 2 from both are the cube of 16 less the inner cube grown by 2,
    // whose edges and corners round: by Steiner's formula its volume is
    // `a³ + 6a²r + 3πar² + 4πr³/3`.
    let (a, inner, w) = (20.0, 4.0, 2.0);
    let mut part = cuboid([0.0; 3], [a; 3]);
    let hole = cuboid([8.0; 3], [12.0; 3]);
    part.positions.extend_from_slice(&hole.positions);
    part.triangles.extend(
        hole.triangles
            .iter()
            .map(|&[x, y, z]| [x + 8, z + 8, y + 8]),
    );
    let done = hollow(&part, w, &Grid::default(), T).unwrap();
    sealed(&done.mesh);
    assert_eq!(done.shells, 2);
    let e = slack(&done, 3.0f64.sqrt() * done.cell);
    for d in depths(&part, &done) {
        assert!((d - w).abs() < e, "a cavity vertex {d} deep");
    }
    let s = a - 2.0 * w;
    let grown = (4.0 / 3.0 * PI).mul_add(
        w.powi(3),
        (3.0 * PI * inner).mul_add(w * w, (6.0 * inner * inner).mul_add(w, inner.powi(3))),
    );
    let want = a.powi(3) - inner.powi(3) - (s.powi(3) - grown);
    // The outer cube's edges as in the cube's bound; the grown cube's
    // rounded edges and corners, curvature at most `1 / w`, as the sphere's.
    let l2 = 3.0 * done.cell * done.cell;
    let area = 6.0f64.mul_add(
        s * s,
        6.0 * inner * inner + 4.0 * PI * w * w + 3.0 * PI * inner * w,
    );
    let bound = (12.0 * s * 2.0f64.sqrt()).mul_add(l2, area * slack(&done, l2 / (2.0 * w)));
    let got = done.mesh.volume();
    assert!(
        (got - want).abs() <= bound,
        "volume {got} against {want}, bound {bound}"
    );
}

#[test]
fn a_part_thinner_than_two_walls_is_left_as_it_is() {
    let part = cuboid([0.0; 3], [10.0, 10.0, 1.5]);
    let done = hollow(&part, 1.0, &Grid::default(), T).unwrap();
    assert_eq!(done.shells, 0);
    assert_eq!(done.mesh.positions, part.positions);
    assert_eq!(done.mesh.triangles, part.triangles);
}

#[test]
fn the_same_part_hollows_the_same_way() {
    let part = sphere(4.0, 24, 48);
    let grid = Grid { cell: Some(0.3) };
    let first = hollow(&part, 0.9, &grid, T).unwrap();
    assert_eq!(first, hollow(&part, 0.9, &grid, T).unwrap());
}

#[test]
fn what_is_not_a_solid_or_a_distance_is_refused() {
    let part = cuboid([0.0; 3], [10.0; 3]);
    for wall in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(hollow(&part, wall, &Grid::default(), T).is_err());
    }
    for cell in [0.0, -1.0, f64::NAN] {
        assert!(hollow(&part, 1.0, &Grid { cell: Some(cell) }, T).is_err());
    }
    let mut open = part.clone();
    open.triangles.pop();
    let message = hollow(&open, 1.0, &Grid::default(), T)
        .unwrap_err()
        .to_string();
    assert!(message.contains("part mesh is open"), "{message}");
    let mut inside_out = part;
    for t in &mut inside_out.triangles {
        t.swap(1, 2);
    }
    let message = hollow(&inside_out, 1.0, &Grid::default(), T)
        .unwrap_err()
        .to_string();
    assert!(message.contains("not wound outward"), "{message}");
}

#[test]
fn a_decimated_cavity_meets_itself_only_at_shared_corners() {
    // A ball's cavity, curved to a radius of a few cells: collapses on one
    // side of it reach across to the other unless each new fan is tested
    // against the triangles near it.
    let (r, w) = (6.0, 1.0);
    let part = sphere(r, 48, 96);
    let checked = crate::boolean::check_solid_kept(&part, "part").unwrap();
    let cell = w / 3.0;
    let marched = Lattice::over(checked.bvh(), cell, w)
        .unwrap()
        .march(checked.bvh())
        .unwrap();
    let target = Target::Verified {
        tolerance: cell / 16.0,
        share: 1.0,
        max_triangles: None,
    };
    let done = simplify(&marched, target, T).unwrap();
    assert!(done.target_met);
    assert!(done.mesh.triangle_count() < marched.triangle_count());
    // The cavity faces into itself, so it is a solid turned inside out.
    let mut outward = done.mesh.clone();
    for t in &mut outward.triangles {
        t.swap(1, 2);
    }
    check_solid(&outward, "decimated", 0).unwrap();
}

#[test]
fn a_cavity_is_checked_as_the_part_joined_to_it_is() {
    // Against the part, a cube of 20: a cavity well inside it, the same
    // wound outward, one poking through the part's face, one crossing
    // itself, and one walling in a hole the part already has.
    let part = cuboid([0.0; 3], [20.0; 3]);
    let inward = |lo: P3, hi: P3| {
        let mut c = cuboid(lo, hi);
        for t in &mut c.triangles {
            t.swap(1, 2);
        }
        c
    };
    let mut crossing = inward([4.0; 3], [8.0; 3]);
    let other = inward([6.0; 3], [10.0; 3]);
    crossing.positions.extend_from_slice(&other.positions);
    crossing
        .triangles
        .extend(other.triangles.iter().map(|t| t.map(|v| v + 8)));
    let cases = [
        (inward([5.0; 3], [15.0; 3]), true),
        (cuboid([5.0; 3], [15.0; 3]), false),
        (inward([5.0; 3], [25.0; 3]), false),
        (crossing, false),
    ];
    let checked = crate::boolean::check_solid_kept(&part, "part").unwrap();
    for (cavity, sound) in &cases {
        let alone = checked.check_cavity(cavity, "hollowed");
        let joined = check_solid(&join(&part, cavity).unwrap(), "hollowed", 12);
        assert_eq!(alone.is_ok(), *sound, "{alone:?}");
        assert_eq!(
            alone.as_ref().err().map(ToString::to_string),
            joined.as_ref().err().map(ToString::to_string)
        );
    }

    // A part with a hole of its own, and a cavity that walls it in: the
    // hole's faces then face into the cavity's material.
    let mut holed = cuboid([0.0; 3], [20.0; 3]);
    let hole = inward([9.0; 3], [11.0; 3]);
    holed.positions.extend_from_slice(&hole.positions);
    holed
        .triangles
        .extend(hole.triangles.iter().map(|t| t.map(|v| v + 8)));
    let checked = crate::boolean::check_solid_kept(&holed, "part").unwrap();
    let around = inward([5.0; 3], [15.0; 3]);
    let alone = checked.check_cavity(&around, "hollowed");
    let joined = check_solid(&join(&holed, &around).unwrap(), "hollowed", 24);
    assert!(alone.is_err());
    assert_eq!(
        alone.err().map(|e| e.to_string()),
        joined.err().map(|e| e.to_string())
    );
}
