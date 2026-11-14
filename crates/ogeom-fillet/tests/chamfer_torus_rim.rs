//! Chamfers on the rims where a plane square to a torus's axis meets the
//! torus, in all three spellings, against Pappus.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{PI, TAU};
use ogeom_core::Tolerances;
use ogeom_fillet::{Chamfer, chamfer_edges_with};
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
/// The torus: major radius 5, minor 2, about z; everything above z = 1
/// cut away.
const BIG: f64 = 5.0;
const SMALL: f64 = 2.0;
const LID: f64 = 1.0;

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

/// The torus less everything above z = 1: two rims on z = 1, of radius
/// 5 + √3 and 5 − √3, each between the flat annulus and the tube.
fn sliced(m: &mut Model) -> Shape {
    let torus = ogeom_algo::make_torus(m, Frame::WORLD, BIG, SMALL, T)
        .unwrap()
        .shape;
    let lid = ogeom_algo::make_box(
        m,
        Frame::new(Point::new(-10.0, -10.0, LID), Direction::Z, Direction::X, T).unwrap(),
        (20.0, 20.0, 5.0),
        T,
    )
    .unwrap()
    .shape;
    ogeom_bool::cut(m, &torus, &lid, T).unwrap().shape
}

fn surface_kind(m: &Model, f: &Shape) -> ogeom_geom::SurfaceKind {
    use ogeom_geom::Surface as _;
    let data = m.node(f).unwrap().data().as_face().unwrap();
    m.geometry().surface(data.surface).unwrap().kind()
}

/// Every edge of `s` on z = 1 between the flat and the tube, with the flat's
/// face and the tube's, and whether it is the outer rim.
fn rims(m: &Model, s: &Shape) -> Vec<(Shape, Shape, Shape, bool)> {
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
                .find(|f| surface_kind(m, f) == ogeom_geom::SurfaceKind::Plane)?;
            let tube = beside
                .iter()
                .find(|f| surface_kind(m, f) == ogeom_geom::SurfaceKind::Torus)?;
            let (c, (a, _)) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap()?;
            let p = ogeom_geom::Curve3d::point_at(&c, a, T).unwrap();
            Some((e.clone(), flat.clone(), tube.clone(), p.x.hypot(p.y) > BIG))
        })
        .collect()
}

/// A rim's section in the meridian, `(rho, z)`: the rim at the tube's angle
/// `alpha` (30° outside, 150° inside), the flat's contact `on_flat` along
/// z = 1 away from the rim (toward the axis outside, away from it inside),
/// the tube's `on_tube` round the tube away from the flat; the area the cut
/// takes off (the triangle and the circular segment between the chord to
/// the tube's contact and the tube) and its centroid's distance from the
/// axis.
fn removed(outer: bool, on_flat: f64, on_tube: f64) -> (f64, f64) {
    let (alpha, side, turn): (f64, f64, f64) = if outer {
        (PI / 6.0, -1.0, -1.0)
    } else {
        (5.0 * PI / 6.0, 1.0, 1.0)
    };
    let at = |b: f64| (SMALL.mul_add(b.cos(), BIG), SMALL * b.sin());
    let a = at(alpha);
    let b = (side.mul_add(on_flat, a.0), LID);
    let turned = on_tube / SMALL;
    let c = at(turn.mul_add(turned, alpha));
    let tri = 0.5 * ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs();
    let tri_rho = (a.0 + b.0 + c.0) / 3.0;
    let seg = 0.5 * SMALL * SMALL * (turned - turned.sin());
    let reach = 4.0 * SMALL * (0.5 * turned).sin().powi(3) / (3.0 * (turned - turned.sin()));
    let seg_rho = reach.mul_add(turn.mul_add(0.5 * turned, alpha).cos(), BIG);
    let area = tri + seg;
    (area, (tri * tri_rho + seg * seg_rho) / area)
}

/// The arc round the tube from the rim to where the line from the flat's
/// contact `d` in, leaving the flat at `angle` toward the rim and down,
/// meets the tube.
fn angle_contact(outer: bool, d: f64, angle: f64) -> f64 {
    let (alpha, side, turn): (f64, f64, f64) = if outer {
        (PI / 6.0, -1.0, -1.0)
    } else {
        (5.0 * PI / 6.0, 1.0, 1.0)
    };
    let rim = SMALL.mul_add(alpha.cos(), BIG);
    let c = (side.mul_add(d, rim) - BIG, LID);
    let dir = (-side * angle.cos(), -angle.sin());
    let b = c.0.mul_add(dir.0, c.1 * dir.1);
    let k = c.0.mul_add(c.0, c.1 * c.1) - SMALL * SMALL;
    let t = -b + b.mul_add(b, -k).sqrt();
    let p = (t.mul_add(dir.0, c.0), t.mul_add(dir.1, c.1));
    let mut swept = turn * (p.1.atan2(p.0) - alpha);
    if swept < 0.0 {
        swept += TAU;
    }
    swept * SMALL
}

fn close(v: f64, expected: f64) {
    assert!(
        (v - expected).abs() <= 1e-6 * expected,
        "{v} against {expected}"
    );
}

/// The sliced torus with both rims chamfered by `spec`, checked: valid,
/// no rim left, each rim's arcs deleted into their bevels; its volume
/// lost.
fn lost(spec: impl Fn(&Shape, &Shape) -> Chamfer) -> f64 {
    let mut m = Model::new();
    let s = sliced(&mut m);
    let before = vol(&m, &s);
    let found = rims(&m, &s);
    assert!(
        found.iter().any(|r| r.3) && found.iter().any(|r| !r.3),
        "both rims"
    );
    let specs: Vec<(Shape, Chamfer)> = found
        .iter()
        .map(|(e, flat, tube, _)| (e.clone(), spec(flat, tube)))
        .collect();
    let built = chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!(rims(&m, &built.shape).is_empty(), "a rim left");
    for (e, _, _, _) in &found {
        assert!(built.history.is_deleted(e), "a rim kept");
        assert!(
            !built.history.generated(e).is_empty(),
            "a rim without its bevel"
        );
    }
    before - vol(&m, &built.shape)
}

/// Both rims' revolved sections at legs `on_flat` and `on_tube`.
fn pappus(on_flat: [f64; 2], on_tube: [f64; 2]) -> f64 {
    let (a, r) = removed(true, on_flat[0], on_tube[0]);
    let (b, q) = removed(false, on_flat[1], on_tube[1]);
    TAU * a.mul_add(r, b * q)
}

#[test]
fn a_torus_rim_chamfers_symmetric_to_its_revolved_section() {
    close(
        lost(|_, _| Chamfer::Symmetric(0.3)),
        pappus([0.3; 2], [0.3; 2]),
    );
}

#[test]
fn a_torus_rim_chamfers_by_two_distances() {
    let got = lost(|flat, _| Chamfer::Distances {
        face: flat.clone(),
        on_face: 0.4,
        on_other: 0.25,
    });
    close(got, pappus([0.4; 2], [0.25; 2]));
}

#[test]
fn a_torus_rim_chamfers_by_distance_and_angle() {
    let (d, angle) = (0.3, 40.0_f64.to_radians());
    let got = lost(|flat, _| Chamfer::Angle {
        face: flat.clone(),
        distance: d,
        angle,
    });
    close(
        got,
        pappus(
            [d; 2],
            [
                angle_contact(true, d, angle),
                angle_contact(false, d, angle),
            ],
        ),
    );
}

/// The angle an arc of a rim turns through: its circle's parameter range.
fn arc_turn(m: &Model, e: &Shape) -> f64 {
    let data = m.node(e).and_then(|n| n.data().as_edge()).expect("an edge");
    let Some(ogeom_topo::EdgeRepr::Curve3d { range, .. }) = data.curve3d() else {
        panic!("an edge in space");
    };
    (range.1 - range.0).abs()
}

/// The sliced torus's outer rim, its arcs asked apart: the first 0.5, the
/// rest 0.8. Each arc's bevel is its section turned over the arc's own
/// sweep, the two stepping on the plane through the axis at each joint
/// (ADR 0001, amendment "M5b"): the volume lost is each arc's share of
/// its section's ring, exactly (Pappus).
#[test]
fn a_torus_rim_s_arcs_asked_apart_step_at_their_joints() {
    let mut m = Model::new();
    let s = sliced(&mut m);
    let before = vol(&m, &s);
    let outer: Vec<Shape> = rims(&m, &s)
        .into_iter()
        .filter(|r| r.3)
        .map(|r| r.0)
        .collect();
    assert!(outer.len() >= 2, "the outer rim's arcs");
    let turns: f64 = outer.iter().map(|e| arc_turn(&m, e)).sum();
    assert!((turns - TAU).abs() < 1e-9, "the arcs close the rim");
    let ds: Vec<f64> = (0..outer.len())
        .map(|k| if k == 0 { 0.5 } else { 0.8 })
        .collect();
    let specs: Vec<(Shape, Chamfer)> = outer
        .iter()
        .zip(&ds)
        .map(|(e, d)| (e.clone(), Chamfer::Symmetric(*d)))
        .collect();
    let built = chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let want: f64 = outer
        .iter()
        .zip(&ds)
        .map(|(e, d)| {
            let (area, rho) = removed(true, *d, *d);
            area * rho * arc_turn(&m, e)
        })
        .sum();
    close(before - vol(&m, &built.shape), want);
    for e in &outer {
        assert_eq!(built.history.generated(e).len(), 1, "its bevel");
    }
}
