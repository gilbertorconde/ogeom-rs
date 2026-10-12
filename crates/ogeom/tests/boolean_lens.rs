//! A torus's tube poking a little way through a drum's or a cone's flat
//! end: what of the tube stands past the end is a thin lens, bounded by a
//! patch of the torus and a patch of the end on one loop.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

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

/// The cut, the common, both ways round, each valid, adding up to the two.
fn holds_the_identities(model: &mut Model, a: &Shape, b: &Shape) -> Shape {
    let (va, vb) = (volume(model, a), volume(model, b));
    let common = ogeom::boolean::common(model, a, b, T).unwrap().shape;
    let cut = ogeom::boolean::cut(model, a, b, T).unwrap().shape;
    let back = ogeom::boolean::cut(model, b, a, T).unwrap().shape;
    for (name, shape) in [("common", &common), ("cut", &cut), ("back", &back)] {
        let diagnosis = ogeom::algo::check(model, shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
    }
    let (vc, vk, vback) = (
        volume(model, &common),
        volume(model, &cut),
        volume(model, &back),
    );
    let slack = (va + vb) * 1e-8;
    assert!(
        (vk + vc - va).abs() <= slack,
        "cut {vk} + common {vc} against {va}"
    );
    assert!(
        (vback + vc - vb).abs() <= slack,
        "back {vback} + common {vc} against {vb}"
    );
    cut
}

/// Leaning, its tube's tip poking a little way past a drum's foot. The cut keeps
/// the tip as a lump of its own, bounded by the patch of the torus past the
/// foot and the patch of the foot inside the tube: two faces on one loop,
/// facing apart, and as thin as a weld only near the loop. Read there
/// alone they looked like a membrane, and the tip was dropped: 0.009 cubic
/// millimetres short.
#[test]
fn a_tube_tip_past_a_drums_foot_is_kept_as_a_lump() {
    let mut model = Model::new();
    let ring = ogeom::algo::make_torus(
        &mut model,
        frame(
            [0.5, 1.5, 1.381_881_667_504_546_7],
            [
                -0.469_591_847_504_161,
                -0.559_580_976_308_394_7,
                0.682_900_159_402_070_2,
            ],
            [0.0, 0.773_488_305_166_557, 0.633_810_572_466_701],
        ),
        2.75,
        1.058_694_984_756_591_4,
        T,
    )
    .unwrap()
    .shape;
    let drum = ogeom::algo::make_cylinder(
        &mut model,
        frame(
            [-1.874_384_946_373_919_7, 0.262_001_941_559_123_7, -1.5],
            [
                0.656_664_884_685_214_8,
                0.603_475_753_123_094_4,
                0.452_336_428_572_658_4,
            ],
            [0.0, 0.599_770_525_268_392_8, -0.800_172_054_635_299_2],
        ),
        2.584_210_061_937_432,
        4.314_983_524_289_919,
        T,
    )
    .unwrap()
    .shape;
    let cut = holds_the_identities(&mut model, &ring, &drum);
    assert_eq!(
        explore_unique(&model, &cut, ShapeType::Solid)
            .unwrap()
            .len(),
        2
    );
}

/// The same past a cone's narrow end.
#[test]
fn a_tube_tip_past_a_cones_end_is_kept_as_a_lump() {
    let mut model = Model::new();
    let ring = ogeom::algo::make_torus(
        &mut model,
        frame(
            [
                -0.199_200_721_145_933,
                0.569_263_811_986_817_6,
                -1.707_653_885_488_181_3,
            ],
            [
                0.617_244_040_520_762_8,
                -0.243_807_884_407_904_68,
                -0.748_042_451_965_224_3,
            ],
            [
                -0.786_771_754_984_635_4,
                -0.191_273_978_417_912_33,
                -0.586_859_839_091_567_9,
            ],
        ),
        2.25,
        0.534_582_073_281_834_8,
        T,
    )
    .unwrap()
    .shape;
    let cone = ogeom::algo::make_cone(
        &mut model,
        frame(
            [
                -1.254_377_920_460_736_8,
                -2.658_596_167_693_578_3,
                -2.441_556_532_810_962,
            ],
            [
                0.793_379_254_183_003_7,
                0.521_952_241_005_792_1,
                -0.313_233_486_621_485_04,
            ],
            [
                -0.608_727_655_879_064_3,
                0.680_281_363_412_529_3,
                -0.408_249_810_240_702_6,
            ],
        ),
        3.397_537_904_225_782_5,
        1.159_425_485_439_682,
        2.5,
        T,
    )
    .unwrap()
    .shape;
    let cut = holds_the_identities(&mut model, &ring, &cone);
    assert_eq!(
        explore_unique(&model, &cut, ShapeType::Solid)
            .unwrap()
            .len(),
        2
    );
}
