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

/// What chamfering every rim arc of the square hole takes off, by the
/// boolean: each wall's bevel is the cone about the wall's normal through
/// its two contact rings (`on_cap` in along the wall, `on_sphere` down the
/// sphere's arc), and the bevels of two walls meeting at a corner of the
/// hole meet on the plane bisecting it, so each wall's share is the hole
/// beyond its wall, between the two bisecting planes, outside its cone. The
/// solid is the same about every bisecting plane, so the four walls take
/// off the same. The pieces past each wall's contacts end on the corner's
/// bisecting planes: the mitre a planar pocket's corner takes, without
/// which a notch of solid stands uncut where the arcs end.
fn mitred_arcs(on_cap: f64, on_sphere: f64) -> f64 {
    let mut m = Model::new();
    let s = square_hole(&mut m);
    // The wall at x = 2, its axis the x axis: the cap contact and the
    // sphere's, `(x, rho)`.
    let r = 21.0_f64.sqrt();
    let lat = (-2.0_f64).atan2(r) - on_sphere / 5.0;
    let cap = (2.0, r - on_cap);
    let wall = (2.0 - 5.0f64.mul_add(lat.sin(), 2.0), 5.0 * lat.cos());
    let slope = (wall.1 - cap.1) / (wall.0 - cap.0);
    let rho = |x: f64| slope.mul_add(x - cap.0, cap.1);
    let along_x = Frame::new(Point::new(1.5, 0.0, 0.0), Direction::X, Direction::Y, T).unwrap();
    let core = ogeom_algo::make_cone(&mut m, along_x, rho(1.5), rho(6.5), 5.0, T)
        .unwrap()
        .shape;
    // The hole's quarter beyond the wall: 2 <= x <= 6, |y| <= x, |z| <= 6.
    let points: Vec<Point> = [(2.0, -2.0), (2.0, 2.0), (6.0, 6.0), (6.0, -6.0)]
        .iter()
        .flat_map(|(x, y)| [Point::new(*x, *y, -6.0), Point::new(*x, *y, 6.0)])
        .collect();
    let rings = [
        vec![0, 2, 4, 6],
        vec![1, 7, 5, 3],
        vec![0, 1, 3, 2],
        vec![2, 3, 5, 4],
        vec![4, 5, 7, 6],
        vec![6, 7, 1, 0],
    ];
    let quarter = ogeom_algo::make_polyhedron(&mut m, &points, &rings, T)
        .unwrap()
        .shape;
    let cutter = ogeom_bool::cut(&mut m, &quarter, &core, T).unwrap().shape;
    let piece = ogeom_bool::common(&mut m, &s, &cutter, T).unwrap().shape;
    4.0 * vol(&m, &piece)
}

#[test]
fn a_square_hole_through_a_sphere_chamfers_its_rim_arcs_symmetric() {
    let (before, v) = hole_chamfered(&|_, _| Chamfer::Symmetric(0.5));
    close(v, before - mitred_arcs(0.5, 0.5));
}

#[test]
fn a_square_hole_through_a_sphere_chamfers_its_rim_arcs_by_two_distances() {
    let (before, v) = hole_chamfered(&|flat, _| Chamfer::Distances {
        face: flat.clone(),
        on_face: 0.8,
        on_other: 0.4,
    });
    close(v, before - mitred_arcs(0.8, 0.4));
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
        before - mitred_arcs(d, cap_angle_contact_about(d, angle, 2.0)),
    );
}

/// A sphere of radius 5 between the planes z = ±1: a zone of the sphere
/// two rims of radius √24 bound.
fn zone(m: &mut Model) -> Shape {
    let s = ball(m);
    let slab = ogeom_algo::make_box(m, at(-10.0, -10.0, -1.0), (20.0, 20.0, 2.0), T)
        .unwrap()
        .shape;
    ogeom_bool::common(m, &s, &slab, T).unwrap().shape
}

/// The zone's two rims chamfered 1.2: each sphere contact 1.2 down the
/// arc from its rim, and the arc between the rims is about 2.01, so the
/// two contacts cross on the sphere, which goes; the two bevel cones meet
/// in their ridge, the circle at z = 0 where their radii agree. The
/// boolean leaves each rim in arcs split at different places top and
/// bottom, so the ridge is split at every one of their corners. Valid, the
/// sphere deleted, the flats and a cone on each rim arc left, the same in
/// both orders. The solid is the zone in common with the inside of each
/// bevel's cone, through its flat's contact and its sphere's: by booleans,
/// which share nothing with the engine.
#[test]
fn a_thin_zone_s_rim_bevels_consume_the_sphere_and_meet_in_a_ridge() {
    let d = 1.2_f64;
    let mut m = Model::new();
    let s = zone(&mut m);
    let found = rims(&m, &s);
    assert!(
        found.len() >= 2,
        "the zone's two rims, as the pieces the boolean leaves"
    );
    let spheres: Vec<Shape> = found.iter().map(|(_, _, f)| f.clone()).collect();
    let specs: Vec<(Shape, Chamfer)> = found
        .iter()
        .map(|(e, _, _)| (e.clone(), Chamfer::Symmetric(d)))
        .collect();
    let back: Vec<(Shape, Chamfer)> = specs.iter().rev().cloned().collect();
    // Each bevel's cone, by its radius at two heights: the flat's contact
    // at the rim's height, the sphere's 1.2 down the arc.
    let rim = 24.0_f64.sqrt();
    let below = (1.0_f64 / 5.0).asin() - d / 5.0;
    let (r_sphere, z_sphere) = (5.0 * below.cos(), 5.0 * below.sin());
    let slope = (rim - d - r_sphere) / (1.0 - z_sphere);
    let radius = |z: f64| r_sphere + slope * (z - z_sphere);
    let mut o = Model::new();
    let mut want = zone(&mut o);
    for up in [1.0_f64, -1.0] {
        let frame = Frame::new(
            Point::new(0.0, 0.0, -2.0 * up),
            if up > 0.0 {
                Direction::Z
            } else {
                Direction::new(ogeom_math::Vector::new(0.0, 0.0, -1.0), T).unwrap()
            },
            Direction::X,
            T,
        )
        .unwrap();
        let cone = ogeom_algo::make_cone(&mut o, frame, radius(-2.0), radius(1.5), 3.5, T).unwrap();
        want = ogeom_bool::common(&mut o, &want, &cone.shape, T)
            .unwrap()
            .shape;
    }
    let want = vol(&o, &want);
    for asked in [&specs, &back] {
        let built = chamfer_edges_with(&mut m, &s, asked, T).unwrap();
        let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        for sphere in &spheres {
            assert!(built.history.is_deleted(sphere), "the sphere stays");
        }
        assert_eq!(
            explore_unique(&m, &built.shape, ShapeType::Face)
                .unwrap()
                .len(),
            2 + specs.len(),
            "the flats and a cone on each rim piece"
        );
        let got = vol(&m, &built.shape);
        assert!(
            (got - want).abs() < 1e-6,
            "volume {got} against the booleans' {want}"
        );
    }
}
