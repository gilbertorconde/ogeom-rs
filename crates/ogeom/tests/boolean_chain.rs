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
