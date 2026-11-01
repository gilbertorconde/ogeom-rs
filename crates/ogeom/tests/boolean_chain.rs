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

/// The chart gap, widest over `shape`'s faces, at each joint of a ring
/// where a seam meets an edge of another kind, read across whole periods:
/// how far the ring stands open on its seam column in the chart.
fn widest_seam_joint(model: &Model, shape: &Shape) -> f64 {
    use ogeom::geom::{Curve2d, Surface};
    use ogeom::math::Point2;
    use ogeom::topo::{EdgeRepr, Orientation, ShapeType};
    let mut widest = 0.0_f64;
    for face in ogeom::topo::explore_unique(model, shape, ShapeType::Face).unwrap() {
        let surface_id = model.node(&face).unwrap().data().as_face().unwrap().surface;
        let surface = model.geometry().surface(surface_id).unwrap();
        let ((u0, u1), (v0, v1)) = surface.domain();
        let turns = |periodic: bool, span: f64| {
            if periodic {
                vec![-span, 0.0, span]
            } else {
                vec![0.0]
            }
        };
        let shifts: Vec<(f64, f64)> = turns(surface.is_periodic_u(), u1 - u0)
            .into_iter()
            .flat_map(|du| {
                turns(surface.is_periodic_v(), v1 - v0)
                    .into_iter()
                    .map(move |dv| (du, dv))
            })
            .collect();
        for wire in model
            .ordered_children_of(&face.oriented(Orientation::Forward))
            .unwrap()
        {
            let edges = model.ordered_children_of(&wire).unwrap();
            // Each edge's ends as the ring walks it, one pair per side, and
            // whether it is a seam.
            let ends: Vec<(bool, Vec<(Point2, Point2)>)> = edges
                .iter()
                .map(|edge| {
                    let data = model.node(edge).unwrap().data().as_edge().unwrap();
                    let back = edge.orientation() == Orientation::Reversed;
                    let walked = |id, range: (f64, f64)| {
                        let pcurve = model.geometry().pcurve(id).unwrap();
                        let (a, b) = (
                            pcurve.point_at(range.0, T).unwrap(),
                            pcurve.point_at(range.1, T).unwrap(),
                        );
                        if back { (b, a) } else { (a, b) }
                    };
                    match data.pcurve_for(surface_id, edge.location()) {
                        Some(EdgeRepr::PCurve { curve, range, .. }) => {
                            (false, vec![walked(*curve, *range)])
                        }
                        Some(EdgeRepr::Seam {
                            forward,
                            reversed,
                            range,
                            ..
                        }) => (
                            true,
                            vec![walked(*forward, *range), walked(*reversed, *range)],
                        ),
                        _ => (false, Vec::new()),
                    }
                })
                .collect();
            for i in 0..edges.len() {
                let j = (i + 1) % edges.len();
                if ends[i].0 == ends[j].0 {
                    continue;
                }
                let gap = ends[i]
                    .1
                    .iter()
                    .flat_map(|a| ends[j].1.iter().map(move |b| (a.1, b.0)))
                    .flat_map(|(end, start)| {
                        shifts.iter().map(move |(du, dv)| {
                            end.distance(Point2::new(start.x + du, start.y + dv))
                        })
                    })
                    .fold(f64::INFINITY, f64::min);
                if gap.is_finite() {
                    widest = widest.max(gap);
                }
            }
        }
    }
    widest
}

/// Two balls fused, their section circle crossing the second ball's seam
/// about a hundredth of a radian below its pole, then cut by a cone. Beside
/// the pole a fitted image's lag of a few microns along the circle is a
/// ten-thousandth of a radian across the chart, past a later boolean's
/// arrangement snap: the fused face's ring closes on its seam column in the
/// chart, or the strand ending on the column dangles. The fused volume is
/// the two balls less their lens, and the cut and the common add up to the
/// fused solid.
#[test]
fn a_fused_ball_face_whose_section_crosses_its_seam_beside_the_pole_is_cut() {
    let mut model = Model::new();
    let (r1, r2) = (3.75, 3.219_816_700_387_659_7);
    let (c1, c2) = (
        Point::new(1.0, 1.5, 2.827_136_898_230_423),
        Point::new(1.5, 1.0, -1.218_746_613_467_500_2),
    );
    let first = ogeom::algo::make_sphere(
        &mut model,
        frame(
            [c1.x, c1.y, c1.z],
            [
                0.039_962_428_829_553_68,
                -0.641_765_914_547_37,
                -0.765_858_678_352_097,
            ],
            [
                -0.999_201_183_086_791_1,
                -0.025_667_027_941_364_065,
                -0.030_630_040_721_722_51,
            ],
        ),
        r1,
        T,
    )
    .unwrap()
    .shape;
    let second = ogeom::algo::make_sphere(
        &mut model,
        frame(
            [c2.x, c2.y, c2.z],
            [
                -0.581_169_508_139_111_3,
                -0.641_184_792_228_283_4,
                0.501_122_804_334_942_3,
            ],
            [-0.0, 0.615_794_499_290_288_7, 0.787_906_805_811_336_2],
        ),
        r2,
        T,
    )
    .unwrap()
    .shape;
    let cone = ogeom::algo::make_cone(
        &mut model,
        frame(
            [-2.5, -1.5, -0.877_100_382_536_820_9],
            [
                0.567_663_994_121_792_6,
                0.394_737_864_230_853_44,
                -0.722_453_879_718_392_7,
            ],
            [
                -0.823_260_341_433_797_9,
                0.272_184_218_482_011_4,
                -0.498_154_756_508_098_8,
            ],
        ),
        3.027_096_920_196_425_5,
        2.156_288_217_038_640_4,
        4.0,
        T,
    )
    .unwrap()
    .shape;
    let fused = ogeom::boolean::fuse(&mut model, &first, &second, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &fused, T).unwrap();
    assert!(diagnosis.is_valid(), "fused: {diagnosis}");
    // Read at a millionth of a radian, the snap below which every
    // arrangement's nodes meet.
    let open = widest_seam_joint(&model, &fused);
    assert!(
        open <= 1e-6,
        "the fused ring stands {open} off its seam column"
    );
    // The union of two balls: both, less the lens they share.
    let pi = core::f64::consts::PI;
    let d = c1.distance(c2);
    let lens = pi
        * (r1 + r2 - d).powi(2)
        * (d * d + 2.0 * d * r2 - 3.0 * r2 * r2 + 2.0 * d * r1 + 6.0 * r1 * r2 - 3.0 * r1 * r1)
        / (12.0 * d);
    let whole = 4.0 / 3.0 * pi * (r1.powi(3) + r2.powi(3)) - lens;
    let measured = volume(&model, &fused);
    assert!(
        (measured - whole).abs() <= whole * 1e-9,
        "fused {measured} against {whole}"
    );
    let cut = ogeom::boolean::cut(&mut model, &fused, &cone, T)
        .unwrap()
        .shape;
    let common = ogeom::boolean::common(&mut model, &fused, &cone, T)
        .unwrap()
        .shape;
    for (name, shape) in [("cut", &cut), ("common", &common)] {
        let diagnosis = ogeom::algo::check(&model, shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
    }
    let (kept, shared) = (volume(&model, &cut), volume(&model, &common));
    assert!(
        (kept + shared - measured).abs() <= measured * 1e-9,
        "cut {kept} + common {shared} against {measured}"
    );
}

/// A torus and a ball fused, then cut by a second torus whose sections
/// cross the fused seam twice a hundredth of a millimetre apart, where the
/// seam bows between the crossings by more than that. The arc of each
/// section between the crossings lies inside the other operand and on
/// neither fused face; read against a face outline whose chords cut across
/// the bow, both arcs read inside their faces, and the kept pieces did not
/// close. The cut and the common add up to the fused solid, and the fused
/// solid and the first two's common to the torus and the ball.
#[test]
fn a_section_crossing_a_bowed_fused_seam_twice_close_together_is_trimmed_by_it() {
    let mut model = Model::new();
    let (major, minor) = (2.667_455_290_041_892_6, 0.944_581_982_922_587_7);
    let ring = ogeom::algo::make_torus(
        &mut model,
        frame(
            [
                2.190_991_776_198_397_5,
                -1.672_874_258_614_358_6,
                -1.301_708_966_940_314_8,
            ],
            [
                -0.342_447_767_989_901,
                -0.718_629_797_426_728_3,
                0.605_227_841_766_350_5,
            ],
            [
                -0.939_536_867_929_479,
                0.261_930_296_234_303_2,
                -0.220_596_903_232_740_87,
            ],
        ),
        major,
        minor,
        T,
    )
    .unwrap()
    .shape;
    let radius = 1.887_274_072_171_815_3;
    let ball = ogeom::algo::make_sphere(
        &mut model,
        frame(
            [
                1.749_738_635_071_303,
                0.296_105_230_908_845_66,
                -1.327_061_519_051_863_8,
            ],
            [
                -0.780_490_395_450_219_7,
                -0.597_201_210_969_75,
                0.184_893_094_046_868_05,
            ],
            [-0.0, 0.295_749_561_085_824, 0.955_265_511_320_042_7],
        ),
        radius,
        T,
    )
    .unwrap()
    .shape;
    let cutter = ogeom::algo::make_torus(
        &mut model,
        frame(
            [
                1.073_495_874_185_432_8,
                -0.432_176_061_383_354_07,
                -0.508_354_819_639_134_1,
            ],
            [
                0.487_571_345_070_652_1,
                -0.581_791_432_783_320_9,
                -0.650_993_788_146_957_5,
            ],
            [
                -0.873_083_148_082_698_6,
                -0.324_900_133_573_391_2,
                -0.363_546_035_239_001_7,
            ],
        ),
        1.831_766_660_605_446_6,
        0.681_031_088_147_804,
        T,
    )
    .unwrap()
    .shape;
    let fused = ogeom::boolean::fuse(&mut model, &ring, &ball, T)
        .unwrap()
        .shape;
    let lens = ogeom::boolean::common(&mut model, &ring, &ball, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(&mut model, &fused, &cutter, T)
        .unwrap()
        .shape;
    let common = ogeom::boolean::common(&mut model, &fused, &cutter, T)
        .unwrap()
        .shape;
    for (name, shape) in [("fused", &fused), ("cut", &cut), ("common", &common)] {
        let diagnosis = ogeom::algo::check(&model, shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
    }
    let pi = core::f64::consts::PI;
    let both = 2.0 * pi * pi * major * minor * minor + 4.0 / 3.0 * pi * radius.powi(3);
    let (whole, shared) = (volume(&model, &fused), volume(&model, &lens));
    assert!(
        (whole + shared - both).abs() <= both * 1e-7,
        "fused {whole} + common {shared} against {both}"
    );
    let (kept, inside) = (volume(&model, &cut), volume(&model, &common));
    assert!(
        (kept + inside - whole).abs() <= whole * 1e-7,
        "cut {kept} + common {inside} against {whole}"
    );
}

/// Two primitives fused and the fused solid cut by a third: whether the
/// cut and the common hold and add up to the fused solid, and the fused
/// solid and the first two's common to `first` + `second`, the two's
/// volumes by closed form.
fn fused_then_cut(
    model: &mut Model,
    (first, second, tool): (Shape, Shape, Shape),
    apart: f64,
) -> (f64, f64, f64) {
    let fused = ogeom::boolean::fuse(model, &first, &second, T)
        .unwrap()
        .shape;
    let both = ogeom::boolean::common(model, &first, &second, T)
        .unwrap()
        .shape;
    let whole = volume(model, &fused);
    let shared = volume(model, &both);
    assert!(
        (whole + shared - apart).abs() <= apart * 1e-7,
        "fused {whole} + common {shared} against {apart}"
    );
    let cut = ogeom::boolean::cut(model, &fused, &tool, T).unwrap().shape;
    let common = ogeom::boolean::common(model, &fused, &tool, T)
        .unwrap()
        .shape;
    for (name, shape) in [("cut", &cut), ("common", &common)] {
        let diagnosis = ogeom::algo::check(model, shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
    }
    let (kept, taken) = (volume(model, &cut), volume(model, &common));
    assert!(
        (kept + taken - whole).abs() <= whole * 1e-7,
        "cut {kept} + common {taken} against {whole}"
    );
    (whole, kept, taken)
}

/// A torus and a ball that do not meet, fused into a compound of the two,
/// then cut by a second ball whose section circle with the first passes
/// 1.7e-4 from the first ball's pole. The section's image swings half a
/// turn of longitude within that hair of arc, and the fused compound's
/// ball was left on its own chart, a compound being no ball: no image fits
/// and the marched pair misses its trace, so the cut was refused. Each
/// member of a compound is charted as a solid operand is. The cut and the
/// common add up to the two primitives, and what the tool takes from the
/// ball is the lens of the two balls.
#[test]
fn a_fused_compound_s_ball_whose_section_passes_beside_its_pole_is_cut() {
    let mut model = Model::new();
    let (major, minor) = (3.75, 1.359_517_775_997_370_7);
    let ring = ogeom::algo::make_torus(
        &mut model,
        frame(
            [-1.894_339_489_500_740_6, -2.5, -2.0],
            [
                0.482_773_091_203_115_96,
                0.142_371_002_235_688,
                0.864_095_272_601_692_6,
            ],
            [
                -0.875_745_478_098_624_2,
                0.078_485_005_707_638_17,
                0.476_350_670_692_174_2,
            ],
        ),
        major,
        minor,
        T,
    )
    .unwrap()
    .shape;
    let (r1, r2) = (1.456_050_108_057_675, 4.0);
    let (c1, c2) = (
        Point::new(-0.5, -0.078_066_247_923_747_5, 2.789_663_657_902_101_7),
        Point::new(
            1.344_172_899_361_099_7,
            1.436_321_666_937_523_4,
            0.697_185_813_805_528_9,
        ),
    );
    let ball = ogeom::algo::make_sphere(
        &mut model,
        frame(
            [c1.x, c1.y, c1.z],
            [
                -0.814_119_917_054_905_2,
                0.510_449_050_795_048_7,
                0.276_858_316_105_816_74,
            ],
            [
                -0.580_696_788_913_555_5,
                -0.715_634_642_429_343,
                -0.388_147_263_162_465,
            ],
        ),
        r1,
        T,
    )
    .unwrap()
    .shape;
    let tool = ogeom::algo::make_sphere(
        &mut model,
        frame(
            [c2.x, c2.y, c2.z],
            [
                -0.993_281_840_366_193_7,
                -0.013_779_370_014_499_021,
                -0.114_896_973_679_687_94,
            ],
            [
                0.013_688_097_555_423_431,
                -0.999_905_059_974_197_7,
                0.001_583_358_237_950_541,
            ],
        ),
        r2,
        T,
    )
    .unwrap()
    .shape;
    let pi = core::f64::consts::PI;
    let apart = 2.0 * pi * pi * major * minor * minor + 4.0 / 3.0 * pi * r1.powi(3);
    let tube = ogeom::boolean::common(&mut model, &ring, &tool, T)
        .unwrap()
        .shape;
    let from_tube = volume(&model, &tube);
    let (_, _, taken) = fused_then_cut(&mut model, (ring, ball, tool), apart);
    let d = c1.distance(c2);
    let lens = pi
        * (r1 + r2 - d).powi(2)
        * (d * d + 2.0 * d * r2 - 3.0 * r2 * r2 + 2.0 * d * r1 + 6.0 * r1 * r2 - 3.0 * r1 * r1)
        / (12.0 * d);
    assert!(
        (taken - from_tube - lens).abs() <= lens * 1e-7,
        "the tool takes {} from the ball against the lens {lens}",
        taken - from_tube
    );
}

/// A ball and a cone fused, then cut by a box one of whose planes passes
/// 1.4e-4 from a pole of the fused ball face: it cuts a small circle round
/// the pole, whose image the chart cannot hold. The fused face was charted
/// clear only of curved walls, so it stayed on its own chart, no image
/// fitted and the marched pair missed its trace: the cut was refused. A
/// plane counts against a sphere face's pole as it does against a whole
/// ball's. The cut and the common add up to the fused solid, and the fused
/// solid and the two's common to the ball and the cone.
#[test]
fn a_fused_ball_face_whose_pole_a_plane_passes_beside_is_cut() {
    let mut model = Model::new();
    let radius = 2.75;
    let ball = ogeom::algo::make_sphere(
        &mut model,
        frame(
            [
                0.665_796_095_297_923_6,
                1.923_745_322_398_347_4,
                -1.320_731_024_509_681,
            ],
            [
                -0.592_965_217_962_411_5,
                -0.286_372_249_796_456_15,
                0.752_584_337_355_824_8,
            ],
            [-0.0, 0.934_622_577_331_156_5, 0.355_641_164_578_662_57],
        ),
        radius,
        T,
    )
    .unwrap()
    .shape;
    let (r0, r1, h) = (
        1.289_984_701_409_45,
        0.494_460_021_109_275_7,
        7.635_607_607_293_029_5,
    );
    let cone = ogeom::algo::make_cone(
        &mut model,
        frame(
            [1.783_218_488_520_477_2, 1.0, -0.036_682_021_305_314_28],
            [
                -0.764_175_441_321_044_6,
                -0.227_742_590_617_548_42,
                -0.603_464_338_052_046_4,
            ],
            [
                -0.645_008_445_589_503,
                0.269_818_629_326_169_85,
                0.714_955_951_361_173_3,
            ],
        ),
        r0,
        r1,
        h,
        T,
    )
    .unwrap()
    .shape;
    let block = ogeom::algo::make_box(
        &mut model,
        frame(
            [-2.0, 1.5, 2.0],
            [
                -0.596_315_614_081_461_8,
                0.490_113_572_718_455_57,
                -0.635_764_401_519_777_2,
            ],
            [0.0, -0.791_982_982_762_753_8, -0.610_543_163_923_904_4],
        ),
        (3.662_395_082_978_046, 7.5, 7.0),
        T,
    )
    .unwrap()
    .shape;
    let pi = core::f64::consts::PI;
    let apart = 4.0 / 3.0 * pi * radius.powi(3) + pi * h / 3.0 * (r0 * r0 + r0 * r1 + r1 * r1);
    fused_then_cut(&mut model, (ball, cone, block), apart);
}
