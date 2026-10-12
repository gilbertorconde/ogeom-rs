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
