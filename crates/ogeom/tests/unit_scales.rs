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

/// Two drums 4 mm tall, one of radius 2 mm and one 2 mm higher, either a
/// thousandth of a millimetre wider or with its axis that far off: the
/// wider cap leaves an annulus or a crescent that thin. Fuse, common and
/// cut keep their volumes to the closed forms, in cubic millimetres, in
/// every unit.
#[test]
fn drums_a_thousandth_apart_keep_their_volumes_in_every_unit() {
    use core::f64::consts::PI;
    for (radius, offset) in [(2.001, 0.0), (2.0001, 0.0), (2.0, 1e-3)] {
        for mm_per_unit in [1e-3, 1.0, 25.4, 1e3] {
            let tol = Tolerances::with_scale(mm_per_unit).unwrap();
            let u = |mm: f64| mm / mm_per_unit;
            let mut model = Model::with_tolerances(tol);
            let at = |x: f64, z: f64| {
                Frame::new(Point::new(u(x), 0.0, u(z)), Direction::Z, Direction::X, tol).unwrap()
            };
            let low = make_cylinder(&mut model, at(0.0, 0.0), u(2.0), u(4.0), tol)
                .unwrap()
                .shape;
            let high = make_cylinder(&mut model, at(offset, 2.0), u(radius), u(4.0), tol)
                .unwrap()
                .shape;
            let chord = Deflection::with_chord(u(1e-3)).unwrap();
            let mm3 = |model: &Model, shape: &ogeom::topo::Shape| {
                volume_properties(model, shape, chord, tol).unwrap().mass * mm_per_unit.powi(3)
            };
            let case = format!("radius {radius}, offset {offset}, {mm_per_unit} mm per unit");
            let fuse = ogeom::boolean::fuse(&mut model, &low, &high, tol)
                .unwrap_or_else(|e| panic!("{case}: fuse: {e}"));
            let common = ogeom::boolean::common(&mut model, &low, &high, tol)
                .unwrap_or_else(|e| panic!("{case}: common: {e}"));
            let cut = ogeom::boolean::cut(&mut model, &low, &high, tol)
                .unwrap_or_else(|e| panic!("{case}: cut: {e}"));
            let (f, c, k) = (
                mm3(&model, &fuse.shape),
                mm3(&model, &common.shape),
                mm3(&model, &cut.shape),
            );
            let (va, vb) = (16.0 * PI, 4.0 * PI * radius * radius);
            let lens = 8.0 * (offset / 4.0).acos() - offset / 2.0 * (16.0 - offset * offset).sqrt();
            assert!((f + c - va - vb).abs() < 1e-6, "{case}: {f} + {c}");
            assert!((k + c - va).abs() < 1e-6, "{case}: {k} + {c}");
            assert!((c - 2.0 * lens.min(4.0 * PI)).abs() < 1e-6, "{case}: {c}");
        }
    }
}
