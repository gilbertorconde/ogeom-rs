//! Holes drilled beside the corner of a box rounded on every edge. A hole
//! crosses a fillet's tangent line on a flat face, where the section on
//! the flat and the section on the fillet meet tangentially, and near the
//! corner it crosses fillets and the corner ball together. Oblique tools
//! through a slanted box rounded on every edge take its volume, Steiner's,
//! in two shares.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_box, make_cone, make_cylinder, volume_properties};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

/// Cut and common with each hole are valid and share the box's volume.
#[test]
fn holes_beside_a_rounded_corner_cut_and_fill_valid() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let rounded = ogeom::fillet::fillet_edges(&mut model, &block, &edges, 2.0, T)
        .unwrap()
        .shape;
    let whole = volume(&model, &rounded);
    for (x, y, r) in [
        (4.3, 2.3, 1.4),
        (1.1, 4.45, 1.0),
        (0.65, 4.5, 0.75),
        (1.1, 0.6, 1.4),
        (1.35, 1.9, 1.4),
        (4.6, 3.0, 1.3),
    ] {
        let at = Frame::new(Point::new(x, y, -5.0), Direction::Z, Direction::X, T).unwrap();
        let drill = make_cylinder(&mut model, at, r, 20.0, T).unwrap().shape;
        let mut shares = 0.0;
        for (name, made) in [
            ("cut", ogeom::boolean::cut(&mut model, &rounded, &drill, T)),
            (
                "common",
                ogeom::boolean::common(&mut model, &rounded, &drill, T),
            ),
        ] {
            let made = made.unwrap_or_else(|e| panic!("{name} at ({x}, {y}): {e}"));
            let diagnosis = check(&model, &made.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{name} at ({x}, {y}): {diagnosis}");
            shares += volume(&model, &made.shape);
        }
        assert!(
            (shares - whole).abs() < whole * 1e-4,
            "at ({x}, {y}): {shares} against {whole}"
        );
    }
}

/// A box of sides `dims` on `frame` rounded at `r` on all twelve edges.
fn rounded_box(model: &mut Model, frame: Frame, dims: (f64, f64, f64), r: f64) -> Shape {
    let block = make_box(model, frame, dims, T).unwrap().shape;
    let edges = explore_unique(model, &block, ShapeType::Edge).unwrap();
    ogeom::fillet::fillet_edges(model, &block, &edges, r, T)
        .unwrap()
        .shape
}

/// Steiner's volume of a box of sides `dims` rounded at `r`: the box
/// shrunk by `r`, its faces swept out by `r`, quarter drums along its
/// edges and eighth balls at its corners.
fn steiner(dims: (f64, f64, f64), r: f64) -> f64 {
    use core::f64::consts::PI;
    let (a, b, c) = (dims.0 - 2.0 * r, dims.1 - 2.0 * r, dims.2 - 2.0 * r);
    a * b * c
        + 2.0 * r * (a * b + b * c + a * c)
        + PI * r * r * (a + b + c)
        + 4.0 / 3.0 * PI * r * r * r
}

fn frame(at: (f64, f64, f64), z: (f64, f64, f64), x: (f64, f64, f64)) -> Frame {
    let dir = |v: (f64, f64, f64)| Direction::new(Vector::new(v.0, v.1, v.2), T).unwrap();
    Frame::new(Point::new(at.0, at.1, at.2), dir(z), dir(x), T).unwrap()
}

/// Cut and common of `rounded` with `tool` are valid and add up to `whole`.
/// Their volumes integrate on the exact surfaces and land within 1e-9 of
/// it; a piece kept on the wrong side moves the sum by far more than 1e-8.
fn shares_add_up(model: &mut Model, rounded: &Shape, tool: &Shape, whole: f64) {
    let mut shares = 0.0;
    for (name, made) in [
        ("cut", ogeom::boolean::cut(model, rounded, tool, T)),
        ("common", ogeom::boolean::common(model, rounded, tool, T)),
    ] {
        let made = made.unwrap_or_else(|e| panic!("{name}: {e}"));
        let diagnosis = check(model, &made.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
        shares += volume(model, &made.shape);
    }
    assert!(
        (shares - whole).abs() < whole * 1e-8,
        "{shares} against Steiner's {whole}, off {:.3e}",
        shares - whole
    );
}

/// A cone whose sections with a flat face and the band beside it reach
/// their tangent seam about 1e-4 apart: the sections meet there, and the
/// cut's wires close.
#[test]
fn a_cone_across_a_band_seam_of_a_rounded_box_cuts_and_fills() {
    let mut model = Model::new();
    let dims = (
        2.076_536_253_795_200_3,
        3.666_203_699_331_365,
        5.245_926_180_077_518,
    );
    let r = 0.622_537_196_054_406_9;
    let box_at = frame(
        (
            0.962_314_625_792_049_8,
            0.446_895_910_791_789_6,
            -1.887_915_774_776_346_2,
        ),
        (
            0.816_217_070_618_454_5,
            0.003_544_934_433_293_659_2,
            0.577_734_478_000_831,
        ),
        (
            -0.577_745_353_621_324_3,
            0.005_008_151_048_799_099,
            0.816_201_705_947_766_1,
        ),
    );
    let rounded = rounded_box(&mut model, box_at, dims, r);
    let cone_at = frame(
        (2.0, -0.194_198_616_922_097_06, -1.0),
        (
            -0.682_823_087_852_392_9,
            0.716_101_356_864_785_5,
            -0.144_746_251_737_779_47,
        ),
        (0.0, -0.198_124_102_422_413_02, -0.980_176_943_229_799_3),
    );
    let cone = make_cone(&mut model, cone_at, 1.25, 0.362_179_447_369_744_4, 4.5, T)
        .unwrap()
        .shape;
    shares_add_up(&mut model, &rounded, &cone, steiner(dims, r));
}

/// A drum of radius 3 through a corner ball and its three bands: the
/// pieces kept round the corner close into a shell.
#[test]
fn a_drum_through_a_rounded_corner_cuts_and_fills() {
    let mut model = Model::new();
    let dims = (
        4.331_777_386_478_203,
        5.495_657_236_386_767_5,
        4.455_939_974_181_635,
    );
    let r = 1.391_219_955_259_229_3;
    let box_at = frame(
        (
            0.116_477_004_812_492_49,
            -2.785_948_483_551_330_4,
            0.633_106_713_016_488_7,
        ),
        (
            -0.585_918_674_241_867_6,
            -0.519_831_238_160_654,
            -0.621_670_966_834_281_5,
        ),
        (0.0, -0.767_144_730_889_939_4, 0.641_474_053_931_881_7),
    );
    let rounded = rounded_box(&mut model, box_at, dims, r);
    let drum_at = frame(
        (1.0, -1.0, 0.096_822_504_618_094_38),
        (
            0.467_693_076_303_335_56,
            -0.786_658_129_279_117_3,
            0.403_028_750_856_563_5,
        ),
        (
            -0.883_890_935_793_507_5,
            -0.416_244_296_193_946_7,
            0.213_254_541_588_411_87,
        ),
    );
    let drum = make_cylinder(&mut model, drum_at, 3.0, 4.579_880_319_455_027, T)
        .unwrap()
        .shape;
    shares_add_up(&mut model, &rounded, &drum, steiner(dims, r));
}
