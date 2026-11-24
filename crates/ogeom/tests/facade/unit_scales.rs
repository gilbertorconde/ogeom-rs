//! One part, stated in micrometres, millimetres, inches, metres and
//! kilometres, builds, cuts, checks, measures and reads back the same: the
//! unit scale moves every linear tolerance with the numbers.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_box, make_cylinder, volume_properties};
use ogeom::core::Tolerances;
use ogeom::io::native::{WriteOptions, read, write};
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::Model;

/// A 20 mm cube with a 4 mm square notch through one corner and a bore of
/// radius 3 mm through its middle, every length divided by `mm_per_unit`.
/// Returns the volume read back from a native round trip, over the closed
/// form, both in cubic model units.
fn drilled_cube(mm_per_unit: f64) -> (f64, f64) {
    let tol = Tolerances::with_scale(mm_per_unit).unwrap();
    let u = |mm: f64| mm / mm_per_unit;
    let mut model = Model::with_tolerances(tol);
    let cube = make_box(&mut model, Frame::WORLD, (u(20.0), u(20.0), u(20.0)), tol)
        .unwrap()
        .shape;
    let at = Frame::new(
        Point::new(u(16.0), u(16.0), u(-1.0)),
        Direction::Z,
        Direction::X,
        tol,
    )
    .unwrap();
    let notch = make_box(&mut model, at, (u(8.0), u(8.0), u(22.0)), tol)
        .unwrap()
        .shape;
    let at = Frame::new(
        Point::new(u(10.0), u(10.0), u(-1.0)),
        Direction::Z,
        Direction::X,
        tol,
    )
    .unwrap();
    let drill = make_cylinder(&mut model, at, u(3.0), u(22.0), tol)
        .unwrap()
        .shape;
    let notched = ogeom::boolean::cut(&mut model, &cube, &notch, tol)
        .unwrap()
        .shape;
    let part = ogeom::boolean::cut(&mut model, &notched, &drill, tol)
        .unwrap()
        .shape;
    let diagnosis = check(&model, &part, tol).unwrap();
    assert!(
        diagnosis.is_valid(),
        "{mm_per_unit} mm per unit: {diagnosis:?}"
    );

    let text = write(&model, &[part], WriteOptions::default()).unwrap();
    let (back, roots) = read(&text).unwrap();
    assert_eq!(back.tolerances(), tol);
    let diagnosis = check(&back, &roots[0], tol).unwrap();
    assert!(
        diagnosis.is_valid(),
        "{mm_per_unit} mm per unit, read back: {diagnosis:?}"
    );
    let chord = Deflection::with_chord(u(1e-3)).unwrap();
    let volume = volume_properties(&back, &roots[0], chord, tol)
        .unwrap()
        .mass;
    let want = u(20.0).powi(3)
        - u(4.0) * u(4.0) * u(20.0)
        - core::f64::consts::PI * u(3.0) * u(3.0) * u(20.0);
    (volume, want)
}

#[test]
fn one_part_measures_the_same_in_every_unit() {
    for mm_per_unit in [1e-3, 1.0, 25.4, 1e3, 1e6] {
        let (volume, want) = drilled_cube(mm_per_unit);
        assert!(
            ((volume - want) / want).abs() < 1e-12,
            "{mm_per_unit} mm per unit: {volume} against {want}"
        );
    }
}
