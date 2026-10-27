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

/// The complete elliptic integrals of the first and second kinds, `K(k)`
/// and `E(k)`, by the arithmetic-geometric mean.
fn complete_elliptic(k: f64) -> (f64, f64) {
    let (mut a, mut g) = (1.0_f64, (1.0 - k * k).sqrt());
    let (mut sum, mut weight) = (k * k / 2.0, 0.5);
    while (a - g).abs() > f64::EPSILON * a {
        let c = (a - g) / 2.0;
        (a, g) = (f64::midpoint(a, g), (a * g).sqrt());
        weight *= 2.0;
        sum += weight * c * c;
    }
    let first = core::f64::consts::FRAC_PI_2 / a;
    (first, first * (1.0 - sum))
}

/// A disc of radius `R = 2.5` and height 3, and a drum of radius
/// `r = 0.5` along x through the disc's axis at z = 1.5, from x = 1 to the
/// plane x = 2.5 tangent to the disc's wall. The drum's end cap touches the
/// wall along a line, and the drum's wall meets the disc's in a loop that
/// touches the drum's end rim at z = 1 and z = 2 and turns back from it.
/// The drum's seam stands 0.7 rad round from its top. A walk ending short
/// of a touch, or running over it, leaves the section's arcs apart from
/// the rim's points, and the kept pieces do not close into a shell.
///
/// The drum within the disc is `∫∫ (√(R² − y²) − 1) dy dz` over its
/// section, which with `y = r·sin θ` is
/// `4·r²·R·((1 + k²)·E(k) − (1 − k²)·K(k)) / (3·k²) − π·r²`, `k = r / R`.
/// The checks run with either solid first, the rim at either end of the
/// drum's height.
#[test]
fn a_drum_capped_in_a_plane_tangent_to_a_bore() {
    let pi = core::f64::consts::PI;
    let (big, small) = (2.5_f64, 0.5_f64);
    let k = small / big;
    let (first, second) = complete_elliptic(k);
    let inside = 4.0 * small * small * big * ((1.0 + k * k) * second - (1.0 - k * k) * first)
        / (3.0 * k * k)
        - pi * small * small;
    let seam = Direction::new(
        ogeom::math::Vector::new(0.0, 0.7_f64.sin(), 0.7_f64.cos()),
        T,
    )
    .unwrap();
    for (foot, axis) in [
        (Point::new(1.0, 0.0, 1.5), Direction::X),
        (
            Point::new(big, 0.0, 1.5),
            Direction::new(-ogeom::math::Vector::X, T).unwrap(),
        ),
    ] {
        let mut model = Model::new();
        let disc = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, big, 3.0, T)
            .unwrap()
            .shape;
        let drum = ogeom::algo::make_cylinder(
            &mut model,
            Frame::new(foot, axis, seam, T).unwrap(),
            small,
            big - 1.0,
            T,
        )
        .unwrap()
        .shape;
        for (name, a, b, whole) in [
            ("disc, drum", &disc, &drum, pi * big * big * 3.0),
            ("drum, disc", &drum, &disc, pi * small * small * (big - 1.0)),
        ] {
            let cut = ogeom::boolean::cut(&mut model, a, b, T).unwrap().shape;
            let common = ogeom::boolean::common(&mut model, a, b, T).unwrap().shape;
            for (what, shape, want) in [("cut", &cut, whole - inside), ("common", &common, inside)]
            {
                let diagnosis = ogeom::algo::check(&model, shape, T).unwrap();
                assert!(
                    diagnosis.is_valid(),
                    "{axis:?}, {name}, {what}: {diagnosis}"
                );
                let got = volume(&model, shape);
                assert!(
                    (got - want).abs() <= 1e-9,
                    "{axis:?}, {name}, {what}: {got} against {want}"
                );
            }
        }
    }
}
