//! Periodic faces imported without the seam or pole edge their chart needs:
//! bands between rings that each go round once, caps bounded by one ring,
//! faces reaching a pole with no edge along the pole's row. Healed, they
//! close in the chart, face the way their loops say, and a boolean cuts
//! them.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

#[path = "support/walks.rs"]
mod walks;

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::with_chord(0.05).unwrap(), T)
        .unwrap()
        .mass
}

/// A part whose bands, caps and pole corners all arrived seamless: healed,
/// it is valid and measures the same, and a drill through a spherical
/// pocket and the bands around it cuts, the cut and the common adding up to
/// the part.
#[test]
fn seamless_faces_heal_and_a_drill_cuts_through_them() {
    let mut import = ogeom::io::read_step(&corpus("nist_ftc_06_asme1_rd.stp"), T).unwrap();
    let solid = import.solids[0].clone();
    let model = import.document.model_mut();
    let before = volume(model, &solid);
    let (healed, count) = ogeom::heal::reanchor_periodic_rings(model, &solid, T).unwrap();
    assert!(count > 0);
    let part = healed.shape;
    assert!(ogeom::algo::check(model, &part, T).unwrap().is_valid());
    assert_eq!(walks::edges_walked_one_way(model, &part), 0);
    let whole = volume(model, &part);
    assert!(
        (whole - before).abs() < before * 1e-4,
        "{whole} against {before}"
    );

    // The pocket's floor faces into the sphere: the material is outside it.
    let beside_pole = Point::new(114.3, 46.625, -107.95);
    let within_pocket = Point::new(114.3, 48.625, -107.95);
    let fine = Deflection::with_chord(0.01).unwrap();
    for (at, want) in [
        (beside_pole, ogeom::algo::Containment::In),
        (within_pocket, ogeom::algo::Containment::Out),
    ] {
        let found = ogeom::algo::classify_in_solid(model, &part, at, fine, T).unwrap();
        assert_eq!(found, want, "at {at:?}");
    }

    let frame = Frame::new(
        Point::new(
            81.848_489_699_818_77,
            80.088_274_165_849_1,
            -453.028_484_344_653,
        ),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    let drill = ogeom::algo::make_cylinder(model, frame, 42.991_498, 493.500_258, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(model, &part, &drill, T).unwrap().shape;
    let common = ogeom::boolean::common(model, &part, &drill, T)
        .unwrap()
        .shape;
    for result in [&cut, &common] {
        assert!(ogeom::algo::check(model, result, T).unwrap().is_valid());
        assert_eq!(walks::edges_walked_one_way(model, result), 0);
    }
    let (a, b) = (volume(model, &cut), volume(model, &common));
    assert!(
        (a + b - whole).abs() < whole * 1e-4,
        "{a} + {b} against {whole}"
    );
}

/// How much of a closed triangle mesh lies inside a cylinder along `y`:
/// the length the mesh encloses along each ray parallel to the axis,
/// integrated over the cylinder's disc by the midpoint rule in polar
/// coordinates. Each crossing adds its height signed by which way the
/// triangle faces, so the sum is the length inside, in any order.
#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "grid indices"
)]
fn mesh_inside_cylinder_along_y(
    mesh: &ogeom::topo::Triangulation,
    centre: (f64, f64),
    radius: f64,
    (rings, spokes): (usize, usize),
) -> f64 {
    // Triangles bucketed by their extent in (x, z), so a ray meets few.
    let cell = 2.0;
    let side = (2.0 * radius / cell).ceil() as usize + 1;
    let (x0, z0) = (centre.0 - radius, centre.1 - radius);
    let bucket = |v: f64, lo: f64| ((v - lo) / cell).floor().clamp(0.0, (side - 1) as f64) as usize;
    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); side * side];
    for (t, tri) in mesh.triangles.iter().enumerate() {
        let p = tri.map(|i| mesh.positions[i as usize]);
        let (xl, xh) = (
            p[0].x.min(p[1].x).min(p[2].x),
            p[0].x.max(p[1].x).max(p[2].x),
        );
        let (zl, zh) = (
            p[0].z.min(p[1].z).min(p[2].z),
            p[0].z.max(p[1].z).max(p[2].z),
        );
        if xh < x0 || zh < z0 || xl > x0 + 2.0 * radius || zl > z0 + 2.0 * radius {
            continue;
        }
        for i in bucket(xl, x0)..=bucket(xh, x0) {
            for k in bucket(zl, z0)..=bucket(zh, z0) {
                buckets[i * side + k].push(t);
            }
        }
    }
    let mut volume = 0.0;
    let (dr, dt) = (
        radius / rings as f64,
        core::f64::consts::TAU / spokes as f64,
    );
    for i in 0..rings {
        let r = (i as f64 + 0.5) * dr;
        for j in 0..spokes {
            // An odd offset keeps the rays off the lines the part's edges
            // run along.
            let t = (j as f64 + 0.25 * (i % 4) as f64 + 0.123) * dt;
            let (x, z) = (centre.0 + r * t.cos(), centre.1 + r * t.sin());
            let mut length = 0.0;
            for &tri in &buckets[bucket(x, x0) * side + bucket(z, z0)] {
                let [a, b, c] = mesh.triangles[tri].map(|i| mesh.positions[i as usize]);
                // Barycentric coordinates of (x, z) in the projected
                // triangle.
                let det = (b.x - a.x) * (c.z - a.z) - (c.x - a.x) * (b.z - a.z);
                if det.abs() < 1e-300 {
                    continue;
                }
                let u = ((x - a.x) * (c.z - a.z) - (c.x - a.x) * (z - a.z)) / det;
                let v = ((b.x - a.x) * (z - a.z) - (x - a.x) * (b.z - a.z)) / det;
                if u < 0.0 || v < 0.0 || u + v > 1.0 {
                    continue;
                }
                let y = a.y + u * (b.y - a.y) + v * (c.y - a.y);
                // The outward normal's y component: positive where the ray
                // leaves the solid.
                let ny = (b.z - a.z) * (c.x - a.x) - (b.x - a.x) * (c.z - a.z);
                length += y * ny.signum();
            }
            volume += length * r * dr * dt;
        }
    }
    volume
}

/// A hemisphere bounded by one great circle through both poles: the
/// circle's chart image jumps half a turn at each pole, so healing splits
/// it there into two exact half meridians with a pole edge between them,
/// and the plane holding the circle takes the halves too. The healed part
/// measures what its fine mesh as read measures, less that mesh's sag,
/// and a drill whose wall crosses the hemisphere cuts it: both results
/// valid, the common matching the drill's share of the mesh as read, the
/// cut and the common adding up to the part.
#[test]
fn a_meridian_circle_splits_at_the_poles_and_a_drill_crosses_the_hemisphere() {
    let mut import = ogeom::io::read_step(&corpus("nist_ftc_06_asme1_rd.stp"), T).unwrap();
    let solid = import.solids[0].clone();
    let model = import.document.model_mut();
    let fine = Deflection::with_chord(0.01).unwrap();
    let mesh = ogeom::mesh::triangulate(model, &solid, fine, T).unwrap();
    let read = ogeom::algo::volume_properties(model, &solid, fine, T)
        .unwrap()
        .mass;
    let part = ogeom::heal::reanchor_periodic_rings(model, &solid, T)
        .unwrap()
        .0
        .shape;
    assert!(ogeom::algo::check(model, &part, T).unwrap().is_valid());
    assert_eq!(walks::edges_walked_one_way(model, &part), 0);
    let whole = ogeom::algo::volume_properties(model, &part, fine, T)
        .unwrap()
        .mass;
    // The mesh as read stands off the part's curved faces by at most its
    // chord, over the whole boundary.
    assert!(
        (whole - read).abs() <= 0.01 * mesh.area(),
        "{whole} against {read}"
    );

    // Along y, through the corner at (61.39, -31.75), wide enough that its
    // wall runs through the hemisphere centred at (0, 97.79, -15.875).
    let (x, z, radius) = (61.394_731_853_578_74, -31.75, 56.427_189_749_833_62);
    let frame = Frame::new(
        Point::new(x, -40.471_773_826_894_63, z),
        Direction::Y,
        Direction::new(ogeom::math::Vector::new(0.0, 0.0, -1.0), T).unwrap(),
        T,
    )
    .unwrap();
    let drill = ogeom::algo::make_cylinder(model, frame, radius, 276.523_547_653_789_25, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(model, &part, &drill, T).unwrap().shape;
    let common = ogeom::boolean::common(model, &part, &drill, T)
        .unwrap()
        .shape;
    for result in [&cut, &common] {
        assert!(ogeom::algo::check(model, result, T).unwrap().is_valid());
        assert_eq!(walks::edges_walked_one_way(model, result), 0);
    }
    let (a, b) = (volume(model, &cut), volume(model, &common));
    let inside = mesh_inside_cylinder_along_y(&mesh, (x, z), radius, (300, 1200));
    assert!(
        (a + b - whole).abs() < whole * 1e-4,
        "{a} + {b} against {whole}"
    );
    // Sampled half as finely each way, the estimate moves by 7e-5 of
    // itself.
    assert!(
        (b - inside).abs() < b * 1e-4,
        "common {b} against {inside} measured on the mesh as read"
    );
}
