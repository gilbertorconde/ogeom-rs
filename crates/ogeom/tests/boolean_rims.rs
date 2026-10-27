//! Two solids whose walls cross in one loop that their rims cut into arcs,
//! one of them short and tucked into a corner of both walls' charts.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

fn frame(at: Point, axis: Direction) -> Frame {
    let x = if axis == Direction::X {
        Direction::Y
    } else {
        Direction::X
    };
    Frame::new(at, axis, x, T).unwrap()
}

/// The fuse and the common add up to the two, the cut is the first less
/// the common, and every result checks valid.
fn holds_the_identities(model: &mut Model, a: &Shape, b: &Shape) {
    let (va, vb) = (volume(model, a), volume(model, b));
    let common = ogeom::boolean::common(model, a, b, T).unwrap().shape;
    let fuse = ogeom::boolean::fuse(model, a, b, T).unwrap().shape;
    let cut = ogeom::boolean::cut(model, a, b, T).unwrap().shape;
    for (name, shape) in [("common", &common), ("fuse", &fuse), ("cut", &cut)] {
        let diagnosis = ogeom::algo::check(model, shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
    }
    let (vc, vf, vk) = (
        volume(model, &common),
        volume(model, &fuse),
        volume(model, &cut),
    );
    let slack = (va + vb) * 1e-6;
    assert!(vc > 0.0, "the two overlap");
    assert!(
        (vf + vc - va - vb).abs() <= slack,
        "fuse {vf} + common {vc} against {}",
        va + vb
    );
    assert!(
        (vk + vc - va).abs() <= slack,
        "cut {vk} + common {vc} against {va}"
    );
}

/// Two drums across each other: the second stands on a foot a little
/// below the first's axis and the first ends a little past the second's
/// axis. Their walls cross in one loop, which the second's foot and the
/// first's end cut into a long arc and one 0.09 long between the two rims.
/// The short arc lies inside one pair of the marcher's sampling cells and
/// crosses neither, and it is found where the rims pierce the other wall.
/// Without it the common, the fuse and the cut were refused: the kept
/// pieces did not close.
#[test]
fn drums_whose_rims_cut_a_short_arc_from_their_crossing() {
    let mut model = Model::new();
    let a = ogeom::algo::make_cylinder(
        &mut model,
        frame(
            Point::new(-1.419_880_124_370_016_3, -1.0, 3.0),
            Direction::X,
        ),
        1.773_938_116_486_465_8,
        2.5,
        T,
    )
    .unwrap()
    .shape;
    let b = ogeom::algo::make_cylinder(
        &mut model,
        frame(
            Point::new(0.510_316_713_993_857_4, 1.5, 2.230_205_316_817_85),
            Direction::Z,
        ),
        1.040_566_819_682_834_4,
        4.0,
        T,
    )
    .unwrap()
    .shape;
    holds_the_identities(&mut model, &a, &b);
}

/// A drum lying across a cone near the cone's foot: their crossing is cut
/// by the rims into arcs as the two drums' is, and without the rims'
/// piercings to seed from the common was refused the same way.
#[test]
fn a_drum_and_a_cone_whose_rims_cut_a_short_arc_from_their_crossing() {
    let mut model = Model::new();
    let a = ogeom::algo::make_cylinder(
        &mut model,
        frame(
            Point::new(-0.579_165_590_402_497_3, -1.0, -2.159_607_595_130_682),
            Direction::X,
        ),
        3.421_186_347_144_059,
        5.5,
        T,
    )
    .unwrap()
    .shape;
    let b = ogeom::algo::make_cone(
        &mut model,
        frame(
            Point::new(1.0, 1.169_069_147_829_731, -2.295_292_000_395_93),
            Direction::Z,
        ),
        1.25,
        0.345_033_661_471_566_33,
        7.071_819_750_037_607,
        T,
    )
    .unwrap()
    .shape;
    holds_the_identities(&mut model, &a, &b);
}

/// A quarter of a disc, `x >= 0` and `y <= 0`, and half a drum along `x`
/// under its axis, `y <= -0.5`, ending `reach` past where the bore stands
/// at the drum's lowest line `y = -1`. Only a strip of the drum's height
/// `reach` wide stands past the bore, and the drum's wall meets the bore in
/// one short arc inside it, between two points of the drum's end rim. Each
/// wall is a part of a turn, and the boolean marches each over a window of
/// its face, less than a turn: the arc lies inside one pair of sampling
/// cells and is found only where the end rim pierces the bore. Without the
/// rims of a window the arc was never traced, and the cut and the common
/// were refused: the kept pieces did not close.
///
/// The drum stands past the bore only where `y < -sqrt(2.5^2 - end^2)`,
/// `end` the drum's end, and by less than `reach` along `x`: under `reach`
/// times the area of the drum's section below that line, which the
/// assertions allow as well as the measure's own nanometre.
#[test]
fn a_half_drum_ending_a_hair_past_a_quarter_bore() {
    let pi = core::f64::consts::PI;
    for reach in [1e-5, 1e-4, 1e-3] {
        let mut model = Model::new();
        let disc = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 2.5, 3.0, T)
            .unwrap()
            .shape;
        let cut_away = |model: &mut Model, from: &Shape, corner: Point, size| {
            let block = ogeom::algo::make_box(model, frame(corner, Direction::Z), size, T)
                .unwrap()
                .shape;
            ogeom::boolean::cut(model, from, &block, T).unwrap().shape
        };
        let quarter = cut_away(
            &mut model,
            &disc,
            Point::new(-3.0, 0.0, -1.0),
            (6.0, 3.0, 5.0),
        );
        let quarter = cut_away(
            &mut model,
            &quarter,
            Point::new(-3.0, -3.0, -1.0),
            (3.0, 6.0, 5.0),
        );
        let end = 5.25_f64.sqrt() + reach;
        let drum = ogeom::algo::make_cylinder(
            &mut model,
            frame(Point::new(1.0, -0.5, 1.5), Direction::X),
            0.5,
            end - 1.0,
            T,
        )
        .unwrap()
        .shape;
        let half = cut_away(
            &mut model,
            &drum,
            Point::new(0.0, -0.5, 0.0),
            (4.0, 1.0, 3.0),
        );
        let cut = ogeom::boolean::cut(&mut model, &quarter, &half, T)
            .unwrap()
            .shape;
        let common = ogeom::boolean::common(&mut model, &quarter, &half, T)
            .unwrap()
            .shape;
        let inside = (end - 1.0) * pi / 8.0;
        let below = -(6.25 - end * end).sqrt();
        let past = reach * (below + 1.0) * 2.0 * (0.25 - (below + 0.5).powi(2)).sqrt();
        for (name, shape, low) in [
            ("cut", &cut, 6.25 * pi * 3.0 / 4.0 - inside),
            ("common", &common, inside - past),
        ] {
            let diagnosis = ogeom::algo::check(&model, shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{reach}, {name}: {diagnosis}");
            let got = volume(&model, shape);
            assert!(
                got >= low - 1e-9 && got <= low + past + 1e-9,
                "{reach}, {name}: {got} against {low} to {}",
                low + past
            );
        }
    }
}
