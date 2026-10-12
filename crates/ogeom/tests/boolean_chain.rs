//! Booleans on a boolean's result: its faces bounded by sections, sewn
//! within their tolerances, and cut again.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

fn frame(at: [f64; 3], z: [f64; 3], x: [f64; 3]) -> Frame {
    let direction = |v: [f64; 3]| Direction::new(Vector::new(v[0], v[1], v[2]), T).unwrap();
    Frame::new(
        Point::new(at[0], at[1], at[2]),
        direction(z),
        direction(x),
        T,
    )
    .unwrap()
}

/// A small torus fused with a drum, the drum's wall crossing the tube in a
/// meridian loop: on the fused torus face the loop's pcurve stops 7.6e-5
/// short of the column of the tube's seam it closes on, two hundredths of a
/// micron in space, inside the loop's own radius. A box then cuts the
/// fused solid. Read at a fixed parameter snap, the face's far seam column
/// met nothing and was dropped as a column the face does not use, and the
/// strand that ends on it dangled: the cut and the common were refused. Ends
/// that meet within their radii in space meet, and the cut and the common
/// add up to the fused solid.
#[test]
fn a_fused_torus_face_whose_loop_stops_a_hair_short_of_its_seam_is_cut() {
    let mut model = Model::new();
    let ring = ogeom::algo::make_torus(
        &mut model,
        frame(
            [1.603_794_412_022_94, -2.5, 1.701_797_781_993_951_4],
            [
                -0.397_598_504_970_574_9,
                -0.754_395_506_943_598_8,
                0.522_305_320_620_491_5,
            ],
            [
                -0.917_559_496_079_226_4,
                0.326_895_996_389_312_16,
                -0.226_326_266_039_704_52,
            ],
        ),
        1.0,
        0.3,
        T,
    )
    .unwrap()
    .shape;
    let drum = ogeom::algo::make_cylinder(
        &mut model,
        frame(
            [
                1.358_332_931_288_014_6,
                -1.274_486_087_653_093,
                -0.071_769_504_399_047_61,
            ],
            [
                -0.353_600_360_302_786_74,
                0.700_039_055_334_248_5,
                0.620_412_851_414_662_7,
            ],
            [
                -0.935_396_592_464_255,
                -0.264_630_066_205_497_7,
                -0.234_529_620_445_551_62,
            ],
        ),
        3.760_353_843_229_520_3,
        2.167_250_261_850_533_4,
        T,
    )
    .unwrap()
    .shape;
    let block = ogeom::algo::make_box(
        &mut model,
        frame(
            [2.5, 1.5, -1.5],
            [
                -0.514_798_599_584_068_4,
                -0.607_003_965_430_051_1,
                0.605_416_045_227_144_1,
            ],
            [0.0, 0.706_180_069_765_889_5, 0.708_032_279_677_588_9],
        ),
        (2.065_523_355_786_593_7, 7.788_207_790_478_866, 7.0),
        T,
    )
    .unwrap()
    .shape;
    let fused = ogeom::boolean::fuse(&mut model, &ring, &drum, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(&mut model, &fused, &block, T)
        .unwrap()
        .shape;
    let common = ogeom::boolean::common(&mut model, &fused, &block, T)
        .unwrap()
        .shape;
    for (name, shape) in [("cut", &cut), ("common", &common)] {
        let diagnosis = ogeom::algo::check(&model, shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
    }
    let (whole, kept, shared) = (
        volume(&model, &fused),
        volume(&model, &cut),
        volume(&model, &common),
    );
    assert!(
        (kept + shared - whole).abs() <= whole * 1e-7,
        "cut {kept} + common {shared} against {whole}"
    );
}
