//! Options built from the model's tolerances, not the millimetre defaults:
//! the same part in micrometres, millimetres, inches and metres gives the
//! same result, read in millimetres.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::{Curve3d as _, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Plane, Point, Vector};
use ogeom::topo::{EdgeRepr, Model, ShapeType, explore_unique};

/// What a section came to, in millimetres.
struct Section {
    edges: usize,
    length: f64,
    /// The greatest distance of a sampled point from the torus or the plane.
    stray: f64,
}

/// A torus of radii 8 and 2 mm cut by a plane tilted off its axis: a pair
/// with no closed form, marched at the default marching chord.
fn tilted_torus_section(mm_per_unit: f64) -> Section {
    let tol = Tolerances::with_scale(mm_per_unit).unwrap();
    let u = |mm: f64| mm / mm_per_unit;
    let mut model = Model::new();
    let ring = ogeom::algo::make_torus(&mut model, Frame::WORLD, u(8.0), u(2.0), tol)
        .unwrap()
        .shape;
    let skin = explore_unique(&model, &ring, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            !matches!(
                model.geometry().surface(data.surface),
                Some(SurfaceGeometry::Plane(_))
            )
        })
        .unwrap();
    let normal = Direction::new(Vector::new(0.0, 0.3, 1.0), tol).unwrap();
    let plane = Plane::through(Point::new(0.0, 0.0, u(0.5)), normal);
    let section = ogeom::boolean::section_face(&mut model, &skin, &plane, tol)
        .unwrap_or_else(|e| panic!("{mm_per_unit} mm per unit: {e}"))
        .shape;
    let edges = explore_unique(&model, &section, ShapeType::Edge).unwrap();
    let (mut length, mut stray) = (0.0_f64, 0.0_f64);
    for edge in &edges {
        let data = model.node(edge).unwrap().data().as_edge().unwrap().clone();
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            panic!("an edge with no curve")
        };
        let c = model.geometry().curve(*curve).unwrap();
        let n = 20_000_u32;
        let points: Vec<Point> = (0..=n)
            .map(|k| {
                let t = range.0 + (range.1 - range.0) * f64::from(k) / f64::from(n);
                let p = c.point_at(t, tol).unwrap();
                Point::new(p.x * mm_per_unit, p.y * mm_per_unit, p.z * mm_per_unit)
            })
            .collect();
        length += points.windows(2).map(|w| w[0].distance(w[1])).sum::<f64>();
        for p in &points {
            let off_torus = ((p.x.hypot(p.y) - 8.0).hypot(p.z) - 2.0).abs();
            let off_plane = 0.3_f64.mul_add(p.y, p.z - 0.5).abs() / 0.3_f64.hypot(1.0);
            stray = stray.max(off_torus).max(off_plane);
        }
    }
    Section {
        edges: edges.len(),
        length,
        stray,
    }
}

/// The section's length in millimetres, as sampled in every unit.
const LENGTH_MM: f64 = 74.976_475;

/// The section in `mm_per_unit` is one curve on both surfaces to 1e-7 mm,
/// as long as in millimetres.
fn same_as_in_millimetres(mm_per_unit: f64) {
    let s = tilted_torus_section(mm_per_unit);
    assert_eq!(s.edges, 1, "{mm_per_unit} mm per unit");
    assert!(
        (s.length - LENGTH_MM).abs() < 1e-7 * LENGTH_MM,
        "{mm_per_unit} mm per unit: {} mm against {LENGTH_MM} mm",
        s.length
    );
    assert!(
        s.stray < 1e-7,
        "{mm_per_unit} mm per unit: strays {} mm",
        s.stray
    );
}

#[test]
fn a_marched_section_in_millimetres() {
    same_as_in_millimetres(1.0);
}

// Before the marching chord followed the unit, the section in micrometres
// ran out of points, and in inches it came in two pieces straying 1.2e-5 mm
// from the true curve.

#[test]
fn a_marched_section_is_the_same_in_micrometres() {
    same_as_in_millimetres(1e-3);
}

#[test]
fn a_marched_section_is_the_same_in_inches() {
    same_as_in_millimetres(25.4);
}

/// Slow: in metres the section's fit takes about ten times as long as in
/// the other units.
#[test]
#[ignore = "heavy"]
fn a_marched_section_is_the_same_in_metres() {
    same_as_in_millimetres(1e3);
}
