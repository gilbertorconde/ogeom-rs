//! Chamfers on the rims of a toroidal groove cut into a flat, the material
//! outside the torus's tube, in all three spellings, against Pappus.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{PI, TAU};
use ogeom_core::Tolerances;
use ogeom_fillet::{Chamfer, chamfer_edges_with};
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
/// The groove's torus: major radius 5, minor 2, about z, its centre 1 above
/// the disc's top at z = 4.
const BIG: f64 = 5.0;
const SMALL: f64 = 2.0;
const TOP: f64 = 4.0;
const CENTRE: f64 = 5.0;

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

/// A disc of radius 8, 4 high, less the torus: a groove in its top with two
/// rims on z = 4, of radius 5 + √3 and 5 − √3, each between the flat and the
/// groove.
fn grooved(m: &mut Model) -> Shape {
    let disc = ogeom_algo::make_cylinder(m, Frame::WORLD, 8.0, TOP, T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(0.0, 0.0, CENTRE), Direction::Z, Direction::X, T).unwrap();
    let torus = ogeom_algo::make_torus(m, at, BIG, SMALL, T).unwrap().shape;
    ogeom_bool::cut(m, &disc, &torus, T).unwrap().shape
}

fn surface_kind(m: &Model, f: &Shape) -> ogeom_geom::SurfaceKind {
    use ogeom_geom::Surface as _;
    let data = m.node(f).unwrap().data().as_face().unwrap();
    m.geometry().surface(data.surface).unwrap().kind()
}

/// Every edge of `s` between the flat and the groove, with the flat's face.
fn rims(m: &Model, s: &Shape) -> Vec<(Shape, Shape)> {
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
            beside
                .iter()
                .find(|f| surface_kind(m, f) == ogeom_geom::SurfaceKind::Torus)?;
            Some((e.clone(), flat.clone()))
        })
        .collect()
}

/// A rim's meridian: its angle round the tube from the tube's centre, the
/// way along the flat into the flat, and the way round the tube into the
/// groove.
fn rim(outer: bool) -> (f64, f64, f64) {
    if outer {
        (-PI / 6.0, 1.0, -1.0)
    } else {
        (-5.0 * PI / 6.0, -1.0, 1.0)
    }
}

/// The area a rim's chamfer takes off in the meridian and its centroid's
/// distance from the axis: the triangle of the rim, the flat's contact
/// `on_flat` along the flat and the groove's `on_tube` round the tube,
/// less the circular segment between the chord to the groove's contact and
/// the tube, which the groove holds already.
fn removed(outer: bool, on_flat: f64, on_tube: f64) -> (f64, f64) {
    let (alpha, side, turn) = rim(outer);
    let at = |b: f64| (SMALL.mul_add(b.cos(), BIG), SMALL.mul_add(b.sin(), CENTRE));
    let a = at(alpha);
    let b = (side.mul_add(on_flat, a.0), TOP);
    let turned = on_tube / SMALL;
    let c = at(turn.mul_add(turned, alpha));
    let tri = 0.5 * ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs();
    let tri_rho = (a.0 + b.0 + c.0) / 3.0;
    let seg = 0.5 * SMALL * SMALL * (turned - turned.sin());
    let reach = 4.0 * SMALL * (0.5 * turned).sin().powi(3) / (3.0 * (turned - turned.sin()));
    let seg_rho = reach.mul_add(turn.mul_add(0.5 * turned, alpha).cos(), BIG);
    let area = tri - seg;
    (area, tri.mul_add(tri_rho, -(seg * seg_rho)) / area)
}

/// The arc round the tube from the rim to where the line from the flat's
/// contact `d` in, leaving the flat at `angle` toward the rim and down,
/// first meets the tube.
fn angle_contact(outer: bool, d: f64, angle: f64) -> f64 {
    let (alpha, side, turn) = rim(outer);
    let rim_rho = SMALL.mul_add(alpha.cos(), BIG);
    let c = (side.mul_add(d, rim_rho) - BIG, TOP - CENTRE);
    let dir = (-side * angle.cos(), -angle.sin());
    let b = c.0.mul_add(dir.0, c.1 * dir.1);
    let k = c.0.mul_add(c.0, c.1 * c.1) - SMALL * SMALL;
    let t = -b - b.mul_add(b, -k).sqrt();
    let p = (t.mul_add(dir.0, c.0), t.mul_add(dir.1, c.1));
    let mut swept = turn * (p.1.atan2(p.0) - alpha);
    if swept < 0.0 {
        swept += TAU;
    }
    swept * SMALL
}

/// The grooved disc with both rims chamfered by `spec`, checked: valid, no
/// rim left, every rim deleted into its bevel; the volume it lost.
fn lost(spec: impl Fn(&Shape) -> Chamfer) -> f64 {
    let mut m = Model::new();
    let s = grooved(&mut m);
    let before = vol(&m, &s);
    let found = rims(&m, &s);
    assert!(found.len() >= 2, "both rims");
    let specs: Vec<(Shape, Chamfer)> = found
        .iter()
        .map(|(e, flat)| (e.clone(), spec(flat)))
        .collect();
    let built = chamfer_edges_with(&mut m, &s, &specs, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!(rims(&m, &built.shape).is_empty(), "a rim left");
    for (e, _) in &found {
        assert!(built.history.is_deleted(e), "a rim kept");
        assert!(
            !built.history.generated(e).is_empty(),
            "a rim without its bevel"
        );
    }
    before - vol(&m, &built.shape)
}

fn pappus(on_flat: [f64; 2], on_tube: [f64; 2]) -> f64 {
    let (a, r) = removed(true, on_flat[0], on_tube[0]);
    let (b, q) = removed(false, on_flat[1], on_tube[1]);
    TAU * a.mul_add(r, b * q)
}

fn close(v: f64, expected: f64) {
    assert!(
        (v - expected).abs() <= 1e-6 * expected,
        "{v} against {expected}"
    );
}

#[test]
fn a_groove_s_rims_chamfer_symmetric_to_their_revolved_sections() {
    close(
        lost(|_| Chamfer::Symmetric(0.3)),
        pappus([0.3; 2], [0.3; 2]),
    );
}

#[test]
fn a_groove_s_rims_chamfer_by_two_distances() {
    let got = lost(|flat| Chamfer::Distances {
        face: flat.clone(),
        on_face: 0.4,
        on_other: 0.25,
    });
    close(got, pappus([0.4; 2], [0.25; 2]));
}

/// At 20° off the flat the line from the flat's contact meets the groove's
/// wall; much steeper, the wall bends back over it first and the line runs
/// under the groove.
#[test]
fn a_groove_s_rims_chamfer_by_distance_and_angle() {
    let (d, angle) = (0.3, 20.0_f64.to_radians());
    let got = lost(|flat| Chamfer::Angle {
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
