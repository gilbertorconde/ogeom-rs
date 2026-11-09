//! Chamfers on the rims where a plane meets a sphere, in all three
//! spellings, against Pappus.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{Chamfer, chamfer_edges_with};
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn vol(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(
        m,
        s,
        ogeom_mesh::Deflection {
            chord: 1e-4,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn ball(m: &mut Model) -> Shape {
    ogeom_algo::make_sphere(m, at(0.0, 0.0, 0.0), 5.0, T)
        .unwrap()
        .shape
}

/// A sphere of radius 5 with everything above z = 2 cut away: one
/// circular rim of radius √21 between the flat and the sphere.
fn flat_side(m: &mut Model) -> Shape {
    let s = ball(m);
    let lid = ogeom_algo::make_box(m, at(-10.0, -10.0, 2.0), (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    ogeom_bool::cut(m, &s, &lid, T).unwrap().shape
}

/// A sphere of radius 5 with a 4 by 4 square hole through it along z.
fn square_hole(m: &mut Model) -> Shape {
    let s = ball(m);
    let hole = ogeom_algo::make_box(m, at(-2.0, -2.0, -10.0), (4.0, 4.0, 20.0), T)
        .unwrap()
        .shape;
    ogeom_bool::cut(m, &s, &hole, T).unwrap().shape
}

fn surface_of(m: &Model, f: &Shape) -> ogeom_geom::SurfaceGeometry {
    m.node(f)
        .and_then(|n| n.data().as_face())
        .and_then(|d| m.geometry().surface(d.surface))
        .cloned()
        .expect("a face's surface")
}

/// Every edge of `s` between a plane and a sphere, with the plane's face
/// and the sphere's.
fn rims(m: &Model, s: &Shape) -> Vec<(Shape, Shape, Shape)> {
    let faces = explore_unique(m, s, ShapeType::Face).unwrap();
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter_map(|e| {
            let beside: Vec<Shape> = faces
                .iter()
                .filter(|f| {
                    explore_unique(m, f, ShapeType::Edge)
                        .unwrap()
                        .iter()
                        .any(|x| x.is_same(&e))
                })
                .cloned()
                .collect();
            let flat = beside
                .iter()
                .find(|f| matches!(surface_of(m, f), ogeom_geom::SurfaceGeometry::Plane(_)))?;
            let wall = beside
                .iter()
                .find(|f| matches!(surface_of(m, f), ogeom_geom::SurfaceGeometry::Sphere(_)))?;
            Some((e.clone(), flat.clone(), wall.clone()))
        })
        .collect()
}

/// The rim's section in the meridian half-plane, `(rho, h)` with the
/// sphere's centre at the origin and the flat at `h = 2`: the rim, the
/// flat's contact `on_flat` in, and the sphere's `on_sphere` down its arc.
fn section(on_flat: f64, on_sphere: f64) -> [(f64, f64); 3] {
    let r = 21.0_f64.sqrt();
    let rim_lat = 2.0_f64.atan2(r);
    let lat = rim_lat - on_sphere / 5.0;
    [
        (r, 2.0),
        (r - on_flat, 2.0),
        (5.0 * lat.cos(), 5.0 * lat.sin()),
    ]
}

/// The area of the section the cut takes off and its centroid's distance
/// from the axis: the triangle rim, flat contact, sphere contact, and the
/// circular segment between the chord from rim to sphere contact and the
/// sphere's own arc.
fn removed(on_flat: f64, on_sphere: f64) -> (f64, f64) {
    let [a, b, c] = section(on_flat, on_sphere);
    let tri = 0.5 * ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs();
    let tri_rho = (a.0 + b.0 + c.0) / 3.0;
    let turn = on_sphere / 5.0;
    let seg = 12.5 * (turn - turn.sin());
    let mid = 2.0_f64.atan2(21.0_f64.sqrt()) - 0.5 * turn;
    let seg_rho = 4.0 * 5.0 * (0.5 * turn).sin().powi(3) / (3.0 * (turn - turn.sin())) * mid.cos();
    let area = tri + seg;
    (area, (tri * tri_rho + seg * seg_rho) / area)
}

/// The sphere of radius 5 above the plane `h = -5` less its cap above
/// `h = 2`.
fn flat_whole() -> f64 {
    let pi = core::f64::consts::PI;
    let cap = pi * 9.0 * (15.0 - 3.0) / 3.0;
    4.0 / 3.0 * pi * 125.0 - cap
}

fn close(v: f64, expected: f64) {
    assert!(
        (v - expected).abs() <= 1e-6 * expected,
        "{v} against {expected}"
    );
}

/// The flat-sided sphere chamfered by `spec`, checked and measured.
fn flat_chamfered(spec: impl Fn(&Shape, &Shape) -> Chamfer) -> f64 {
    let mut m = Model::new();
    let s = flat_side(&mut m);
    // The boolean leaves the rim as arcs; each is chamfered.
    let found = rims(&m, &s);
    assert!(!found.is_empty(), "the rim");
    close(vol(&m, &s), flat_whole());
    let specs: Vec<(Shape, Chamfer)> = found
        .iter()
        .map(|(e, flat, wall)| (e.clone(), spec(flat, wall)))
        .collect();
    let built = chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    for (rim, _, _) in &found {
        assert!(built.history.is_deleted(rim));
    }
    assert!(rims(&m, &built.shape).is_empty(), "no rim left");
    vol(&m, &built.shape)
}

fn flat_expected(on_flat: f64, on_sphere: f64) -> f64 {
    let (area, rho) = removed(on_flat, on_sphere);
    flat_whole() - area * core::f64::consts::TAU * rho
}

#[test]
fn a_sphere_rim_chamfers_symmetric_to_the_revolved_section() {
    close(
        flat_chamfered(|_, _| Chamfer::Symmetric(0.5)),
        flat_expected(0.5, 0.5),
    );
}

#[test]
fn a_sphere_rim_chamfers_by_two_distances() {
    close(
        flat_chamfered(|flat, _| Chamfer::Distances {
            face: flat.clone(),
            on_face: 0.8,
            on_other: 0.4,
        }),
        flat_expected(0.8, 0.4),
    );
    close(
        flat_chamfered(|_, wall| Chamfer::Distances {
            face: wall.clone(),
            on_face: 0.8,
            on_other: 0.4,
        }),
        flat_expected(0.4, 0.8),
    );
}

/// Where the line from the flat's contact `d` in, leaving the flat at
/// `angle` down into the material, meets the sphere's meridian circle:
/// the arc length from the rim down to it.
fn flat_angle_contact(d: f64, angle: f64) -> f64 {
    let r = 21.0_f64.sqrt();
    let c = (r - d, 2.0);
    let dir = (angle.cos(), -angle.sin());
    // |c + t dir| = 5, the root ahead of c.
    let b = c.0 * dir.0 + c.1 * dir.1;
    let cc = c.0 * c.0 + c.1 * c.1 - 25.0;
    let t = -b + (b * b - cc).sqrt();
    let p = (c.0 + t * dir.0, c.1 + t * dir.1);
    5.0 * (2.0_f64.atan2(r) - p.1.atan2(p.0))
}

/// Where the line from the sphere's contact `d` down its arc, leaving the
/// sphere's tangent at `angle` into the material, meets the flat: the
/// distance in from the rim.
fn sphere_angle_contact(d: f64, angle: f64) -> f64 {
    let lat = 2.0_f64.atan2(21.0_f64.sqrt()) - d / 5.0;
    let w = (5.0 * lat.cos(), 5.0 * lat.sin());
    // The tangent up the meridian toward the rim, turned by `angle` toward
    // the centre.
    let up = (-lat.sin(), lat.cos());
    let inward = (-lat.cos(), -lat.sin());
    let dir = (
        up.0 * angle.cos() + inward.0 * angle.sin(),
        up.1 * angle.cos() + inward.1 * angle.sin(),
    );
    let t = (2.0 - w.1) / dir.1;
    21.0_f64.sqrt() - (w.0 + t * dir.0)
}

#[test]
fn a_sphere_rim_chamfers_by_distance_and_angle() {
    let (d, angle) = (0.6, 40.0_f64.to_radians());
    close(
        flat_chamfered(|flat, _| Chamfer::Angle {
            face: flat.clone(),
            distance: d,
            angle,
        }),
        flat_expected(d, flat_angle_contact(d, angle)),
    );
    close(
        flat_chamfered(|_, wall| Chamfer::Angle {
            face: wall.clone(),
            distance: d,
            angle,
        }),
        flat_expected(sphere_angle_contact(d, angle), d),
    );
}

/// A 20 by 20 by 10 block, its top at z = 0, with a sphere of radius 5
/// centred 2 above the top taken out: a pocket whose rim of radius √21
/// has the material outside the sphere.
fn pocket(m: &mut Model) -> Shape {
    let block = ogeom_algo::make_box(m, at(-10.0, -10.0, -10.0), (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let s = ogeom_algo::make_sphere(m, at(0.0, 0.0, 2.0), 5.0, T)
        .unwrap()
        .shape;
    ogeom_bool::cut(m, &block, &s, T).unwrap().shape
}

#[test]
fn a_pocket_rim_chamfers_to_the_revolved_section() {
    // In the meridian, sphere centre at (0, 2): the rim (√21, 0), the
    // top's contact d out, the sphere's d down its arc. The cut takes the
    // triangle less the circular segment its chord cuts off the sphere's
    // side.
    let d = 0.5;
    let r = 21.0_f64.sqrt();
    let rim_lat = (-2.0_f64).atan2(r);
    let lat = rim_lat - d / 5.0;
    let (a, b, c) = (
        (r, 0.0),
        (r + d, 0.0),
        (5.0 * lat.cos(), 2.0 + 5.0 * lat.sin()),
    );
    let tri = 0.5 * ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs();
    let tri_rho = (a.0 + b.0 + c.0) / 3.0;
    let turn = d / 5.0;
    let seg = 12.5 * (turn - turn.sin());
    let seg_rho = 4.0 * 5.0 * (0.5 * turn).sin().powi(3) / (3.0 * (turn - turn.sin()))
        * (rim_lat - 0.5 * turn).cos();
    let tau = core::f64::consts::TAU;
    let whole = 4000.0 - 36.0 * core::f64::consts::PI;
    let expected = whole - tau * (tri * tri_rho - seg * seg_rho);

    let mut m = Model::new();
    let s = pocket(&mut m);
    close(vol(&m, &s), whole);
    let found = rims(&m, &s);
    assert!(!found.is_empty(), "the rim");
    let specs: Vec<(Shape, Chamfer)> = found
        .iter()
        .map(|(e, _, _)| (e.clone(), Chamfer::Symmetric(d)))
        .collect();
    let built = chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!(rims(&m, &built.shape).is_empty(), "no rim left");
    close(vol(&m, &built.shape), expected);
}

/// A rim's section in its meridian half-plane, `(rho, h)` from the rim's
/// centre with `h` along the cap's outward normal, the sphere's centre at
/// height `centre` on the axis: the rim, the cap's contact `on_cap` in,
/// and the sphere's `on_sphere` down its arc away from the cap; then the
/// area the cut takes off (the triangle and the circular segment between
/// the chord to the sphere's contact and the sphere's arc) and its
/// centroid's distance from the axis.
fn section_about(on_cap: f64, on_sphere: f64, centre: f64) -> (f64, f64) {
    let r = (25.0 - centre * centre).sqrt();
    let rim_lat = (-centre).atan2(r);
    let lat = rim_lat - on_sphere / 5.0;
    let (a, b, c) = (
        (r, 0.0),
        (r - on_cap, 0.0),
        (5.0 * lat.cos(), 5.0f64.mul_add(lat.sin(), centre)),
    );
    let tri = 0.5 * ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs();
    let tri_rho = (a.0 + b.0 + c.0) / 3.0;
    let turn = on_sphere / 5.0;
    let seg = 12.5 * (turn - turn.sin());
    let seg_rho = 4.0 * 5.0 * (0.5 * turn).sin().powi(3) / (3.0 * (turn - turn.sin()))
        * (rim_lat - 0.5 * turn).cos();
    let area = tri + seg;
    (area, (tri * tri_rho + seg * seg_rho) / area)
}

/// Where the line from the cap's contact `d` in, leaving the cap at `angle`
/// away from it toward the rim, meets the sphere about `(0, centre)`: the
/// arc length from the rim to it.
fn cap_angle_contact_about(d: f64, angle: f64, centre: f64) -> f64 {
    let r = (25.0 - centre * centre).sqrt();
    let c = (r - d, 0.0 - centre);
    let dir = (angle.cos(), -angle.sin());
    let b = c.0.mul_add(dir.0, c.1 * dir.1);
    let k = c.0.mul_add(c.0, c.1 * c.1) - 25.0;
    let t = -b + b.mul_add(b, -k).sqrt();
    let p = (t.mul_add(dir.0, c.0), t.mul_add(dir.1, c.1));
    5.0 * ((-centre).atan2(r) - p.1.atan2(p.0))
}

/// The square hole with every rim arc chamfered by `spec`, checked; its
/// volume and the hole's before.
fn hole_chamfered(spec: &dyn Fn(&Shape, &Shape) -> Chamfer) -> (f64, f64) {
    let mut m = Model::new();
    let s = square_hole(&mut m);
    let before = vol(&m, &s);
    // Each opening's four arcs, as the pieces the boolean leaves.
    let found = rims(&m, &s);
    assert!(found.len() >= 8, "the hole's rim arcs");
    let specs: Vec<(Shape, Chamfer)> = found
        .iter()
        .map(|(e, flat, wall)| (e.clone(), spec(flat, wall)))
        .collect();
    let built = chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    for (rim, _, _) in &found {
        assert!(built.history.is_deleted(rim), "a rim kept");
    }
    (before, vol(&m, &built.shape))
}

/// The eight arcs' worth of the section on a wall of the hole: each wall
/// stands 2 off the sphere's centre with the centre beyond it, so past
/// the rim the sphere closes toward the axis, and each arc turns through
/// 2 asin(2/√21) about its own axis. Every contact stands nearer the axis
/// than the rim, so the wedges on two walls meeting at a corner of the
/// hole stand apart.
fn eight_arcs(on_cap: f64, on_sphere: f64) -> f64 {
    let sweep = 2.0 * (2.0 / 21.0_f64.sqrt()).asin();
    let (area, rho) = section_about(on_cap, on_sphere, 2.0);
    8.0 * area * sweep * rho
}

#[test]
fn a_square_hole_through_a_sphere_chamfers_its_rim_arcs_symmetric() {
    let (before, v) = hole_chamfered(&|_, _| Chamfer::Symmetric(0.5));
    close(v, before - eight_arcs(0.5, 0.5));
}

#[test]
fn a_square_hole_through_a_sphere_chamfers_its_rim_arcs_by_two_distances() {
    let (before, v) = hole_chamfered(&|flat, _| Chamfer::Distances {
        face: flat.clone(),
        on_face: 0.8,
        on_other: 0.4,
    });
    close(v, before - eight_arcs(0.8, 0.4));
}

#[test]
fn a_square_hole_through_a_sphere_chamfers_its_rim_arcs_by_distance_and_angle() {
    let (d, angle) = (0.6, 40.0_f64.to_radians());
    let (before, v) = hole_chamfered(&|flat, _| Chamfer::Angle {
        face: flat.clone(),
        distance: d,
        angle,
    });
    close(
        v,
        before - eight_arcs(d, cap_angle_contact_about(d, angle, 2.0)),
    );
}
