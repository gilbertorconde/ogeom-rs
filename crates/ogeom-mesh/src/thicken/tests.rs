#![allow(clippy::unwrap_used, clippy::cast_precision_loss)]

use super::*;
use ogeom_math::Point;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

/// A square sheet of side `size` in the plane `z = 0`, `cells` squares a
/// side, facing `+z`.
fn square(size: f64, cells: u32) -> Triangulation {
    let mut mesh = Triangulation::new();
    let step = size / f64::from(cells);
    for j in 0..=cells {
        for i in 0..=cells {
            mesh.positions
                .push(Point::new(f64::from(i) * step, f64::from(j) * step, 0.0));
        }
    }
    let at = |i: u32, j: u32| j * (cells + 1) + i;
    for j in 0..cells {
        for i in 0..cells {
            mesh.triangles
                .push([at(i, j), at(i + 1, j), at(i + 1, j + 1)]);
            mesh.triangles
                .push([at(i, j), at(i + 1, j + 1), at(i, j + 1)]);
        }
    }
    mesh
}

/// The half of a sphere about the origin above `z = 0`, of latitude rings
/// and longitude slices, facing out, its rim on the equator.
fn hemisphere(radius: f64, rings: u32, slices: u32) -> Triangulation {
    let mut mesh = Triangulation::new();
    let at = |theta: f64, phi: f64| {
        Point::new(
            radius * theta.sin() * phi.cos(),
            radius * theta.sin() * phi.sin(),
            radius * theta.cos(),
        )
    };
    mesh.positions.push(at(0.0, 0.0));
    for r in 1..=rings {
        let theta = FRAC_PI_2 * f64::from(r) / f64::from(rings);
        for s in 0..slices {
            mesh.positions
                .push(at(theta, TAU * f64::from(s) / f64::from(slices)));
        }
    }
    let ring = |r: u32, s: u32| 1 + (r - 1) * slices + s % slices;
    for s in 0..slices {
        mesh.triangles.push([0, ring(1, s), ring(1, s + 1)]);
        for r in 1..rings {
            let (a, b) = (ring(r, s), ring(r, s + 1));
            let (c, d) = (ring(r + 1, s), ring(r + 1, s + 1));
            mesh.triangles.push([a, c, d]);
            mesh.triangles.push([a, d, b]);
        }
    }
    mesh
}

/// The result is one closed 2-manifold mesh, wound outward.
fn sealed(mesh: &Triangulation) {
    let e = Editable::from(mesh, 0.0).unwrap();
    assert!(e.is_closed(), "the thickened mesh is open");
    assert!(
        e.non_manifold().is_empty(),
        "the thickened mesh is not 2-manifold"
    );
    check_solid(mesh, "thickened", 0).unwrap();
}

#[test]
fn a_thickened_hemisphere_holds_the_shell_between_its_spheres() {
    let (r, d, rings, slices) = (10.0, 1.0, 32, 128);
    let sheet = hemisphere(r, rings, slices);
    // The facets lie between the sphere and `sag` inside it, and each layer
    // between its sphere and `sag` inside that, so each layer's volume is
    // off by at most its area times the sag. The rim of an offset layer
    // moves along the normal of the last ring's faces, tilted from the
    // equator by at most a ring's angle, which lifts the wall by at most
    // `|offset| sin(π / 2 rings)` round the rim.
    let bvh = crate::Bvh::build(&sheet).unwrap();
    let sag = r - bvh.nearest(Point::ORIGIN, f64::INFINITY).unwrap().distance;
    let tilt = (FRAC_PI_2 / f64::from(rings)).sin();
    for side in [Thicken::Outward, Thicken::Inward, Thicken::Both] {
        let done = thicken(&sheet, d, side).unwrap();
        sealed(&done.mesh);
        let (inner, outer, moved) = match side {
            Thicken::Outward => (r, r + d, d),
            Thicken::Inward => (r - d, r, d),
            Thicken::Both => (r - d / 2.0, r + d / 2.0, d / 2.0),
        };
        let want = 2.0 / 3.0 * PI * (outer.powi(3) - inner.powi(3));
        let bound = (2.0 * PI).mul_add(
            (outer * outer + inner * inner) * sag,
            2.0 * TAU * outer * moved * moved * tilt,
        );
        let got = done.mesh.volume();
        assert!(
            (got - want).abs() <= bound,
            "{side:?}: volume {got} against {want}, bound {bound}"
        );
        assert_eq!(done.front.folded + done.back.folded, 0);
    }
}

#[test]
fn a_thickened_square_is_a_slab_on_the_side_asked() {
    let sheet = square(10.0, 10);
    for (side, low, high) in [
        (Thicken::Outward, 0.0, 0.5),
        (Thicken::Inward, -0.5, 0.0),
        (Thicken::Both, -0.25, 0.25),
    ] {
        let done = thicken(&sheet, 0.5, side).unwrap();
        sealed(&done.mesh);
        assert!((done.mesh.volume() - 50.0).abs() < 1e-12, "{side:?}");
        let (lo, hi) = done
            .mesh
            .positions
            .iter()
            .fold((f64::MAX, f64::MIN), |(l, h), p| (l.min(p.z), h.max(p.z)));
        assert!(
            (lo - low).abs() < 1e-12 && (hi - high).abs() < 1e-12,
            "{side:?}: from {lo} to {hi}"
        );
        // The two layers and a wall of two triangles along each rim edge.
        assert_eq!(done.mesh.triangles.len(), 2 * 200 + 2 * 40);
    }
}

#[test]
fn the_same_sheet_thickens_the_same_way() {
    let sheet = hemisphere(3.0, 8, 24);
    let first = thicken(&sheet, 0.4, Thicken::Both).unwrap();
    assert_eq!(first, thicken(&sheet, 0.4, Thicken::Both).unwrap());
}

#[test]
fn a_thickness_past_the_sheet_s_bend_is_refused() {
    // Inward by more than the radius, the inner layer passes through the
    // centre and out through itself.
    let sheet = hemisphere(1.0, 8, 24);
    assert!(thicken(&sheet, 0.5, Thicken::Inward).is_ok());
    let message = thicken(&sheet, 3.0, Thicken::Inward)
        .unwrap_err()
        .to_string();
    assert!(message.contains("thickened"), "{message}");
}

#[test]
fn what_is_not_an_open_sheet_or_a_distance_is_refused() {
    let sheet = square(1.0, 2);
    for d in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(thicken(&sheet, d, Thicken::Both).is_err());
    }
    assert!(thicken(&Triangulation::new(), 1.0, Thicken::Both).is_err());
    // A hemisphere closed by a fan over its rim.
    let mut closed = hemisphere(1.0, 4, 8);
    let centre = u32::try_from(closed.positions.len()).unwrap();
    closed.positions.push(Point::ORIGIN);
    let rim = |s: u32| 1 + 3 * 8 + s % 8;
    for s in 0..8 {
        closed.triangles.push([centre, rim(s + 1), rim(s)]);
    }
    let message = thicken(&closed, 0.1, Thicken::Outward)
        .unwrap_err()
        .to_string();
    assert!(message.contains("closed"), "{message}");
    // Two squares sharing one corner.
    let mut bowtie = square(1.0, 1);
    bowtie.positions.extend([
        Point::new(2.0, 1.0, 0.0),
        Point::new(2.0, 2.0, 0.0),
        Point::new(1.0, 2.0, 0.0),
    ]);
    bowtie.triangles.extend([[3, 4, 5], [3, 5, 6]]);
    let message = thicken(&bowtie, 0.1, Thicken::Outward)
        .unwrap_err()
        .to_string();
    assert!(message.contains("2-manifold"), "{message}");
}
