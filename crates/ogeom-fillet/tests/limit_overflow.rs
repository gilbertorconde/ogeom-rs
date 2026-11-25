//! Largest-size readings under the request's overflow options: each probe
//! builds with them, so a reading answers the request as it will be built,
//! for setback and rolling corners, radius laws and two distances.
//! The part is a block whose top front edge is rounded, with a boss fused
//! on its top near that round: the boss's foot reaches the round, which
//! blend into round takes and the other behaviours do not.
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{
    BlendSizing, Chamfer, ChamferCorner, FilletCorner, LIMIT_TOLERANCE, LawSize, LawStation,
    Overflow, RadiusLaw, RoundProfile, blend_limit, chamfer_edges_overflow,
    fillet_edge_distances_overflow, fillet_edges, fillet_edges_corner_overflow,
    fillet_edges_overflow, fillet_edges_variable_overflow, largest_chamfer_distance_all,
    largest_chamfer_distance_all_overflow, largest_chamfer_distance_overflow,
    largest_fillet_radius_overflow,
};
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn edges_where(m: &Model, s: &Shape, on: impl Fn(Point) -> bool) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            [0.0, 0.3, 0.7, 1.0]
                .iter()
                .all(|t| on(c.point_at(a + (b - a) * t, T).unwrap()))
        })
        .collect()
}

/// The block 40 x 30 x 15, its top front edge (along `y = 0`) rounded at
/// 5, and a boss of radius 6, 15 tall, fused on its top about `(20, 12)`;
/// with the boss's foot rim.
fn bossed() -> (Model, Shape, Vec<Shape>) {
    let mut m = Model::new();
    let block = ogeom_algo::make_box(&mut m, at(0.0, 0.0, 0.0), (40.0, 30.0, 15.0), T)
        .unwrap()
        .shape;
    let front = edges_where(&m, &block, |p| {
        p.y.abs() < 1e-9 && (p.z - 15.0).abs() < 1e-9
    });
    let rounded = fillet_edges(&mut m, &block, &front, 5.0, T).unwrap().shape;
    let boss = ogeom_algo::make_cylinder(&mut m, at(20.0, 12.0, 15.0), 6.0, 15.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(&mut m, &rounded, &boss, T).unwrap().shape;
    let foot = edges_where(&m, &s, |p| {
        (p.z - 15.0).abs() < 1e-9 && ((p.x - 20.0).hypot(p.y - 12.0) - 6.0).abs() < 1e-9
    });
    assert!(!foot.is_empty(), "the foot rim");
    (m, s, foot)
}

fn blend_off() -> Overflow {
    let mut off = Overflow::default();
    off.blend_into_round = false;
    off
}

fn rounds(m: &Model, s: &Shape, foot: &[Shape], r: f64, overflow: Overflow) -> bool {
    let mut scratch = m.clone();
    fillet_edges_overflow(&mut scratch, s, foot, r, overflow, T).is_ok()
}

fn bevels(m: &Model, s: &Shape, foot: &[Shape], d: f64, overflow: Overflow) -> bool {
    let specs: Vec<_> = foot
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(d)))
        .collect();
    let mut scratch = m.clone();
    chamfer_edges_overflow(&mut scratch, s, &specs, ChamferCorner::Facet, overflow, T).is_ok()
}

/// With blend into round off the foot round builds only up to about
/// 1.01, just past where its contact first reaches the edge round's
/// tangent line; the default reading, 8.62, is refused under those
/// options. The reading under them answers 1.01, builds there and is
/// refused a tolerance above. Under the default options it answers what
/// `largest_fillet_radius_all` does.
#[test]
#[ignore = "heavy"]
fn the_foot_round_reads_its_limit_with_blend_into_round_off() {
    let (m, s, foot) = bossed();
    let default = largest_fillet_radius_overflow(&m, &s, &foot, 20.0, Overflow::default(), T)
        .unwrap()
        .unwrap();
    assert!((default - 8.6195).abs() < 1e-3, "{default}");
    assert!(!rounds(&m, &s, &foot, default, blend_off()));

    let found = blend_limit(
        &m,
        &s,
        &BlendSizing::RoundsOverflow {
            edges: &foot,
            overflow: blend_off(),
        },
        20.0,
        T,
    )
    .unwrap();
    let size = found.size.unwrap();
    assert!(
        (size - 1.0100).abs() <= 1.0100 * LIMIT_TOLERANCE,
        "{found:?}"
    );
    assert_eq!(
        largest_fillet_radius_overflow(&m, &s, &foot, 20.0, blend_off(), T).unwrap(),
        Some(size)
    );
    assert!(rounds(&m, &s, &foot, size, blend_off()), "{found:?}");
    assert!(
        !rounds(&m, &s, &foot, size * (1.0 + LIMIT_TOLERANCE), blend_off()),
        "{found:?}"
    );
}

/// The foot bevel of the same part with blend into round off: the reading
/// builds under those options and is refused a tolerance above. Under the
/// default options it answers what `largest_chamfer_distance_all` does.
#[test]
#[ignore = "heavy"]
fn the_foot_bevel_reads_its_limit_with_blend_into_round_off() {
    let (m, s, foot) = bossed();
    let specs: Vec<_> = foot
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(1.0)))
        .collect();
    let off = largest_chamfer_distance_all_overflow(
        &m,
        &s,
        &specs,
        20.0,
        ChamferCorner::Facet,
        blend_off(),
        T,
    )
    .unwrap()
    .unwrap();
    assert!(bevels(&m, &s, &foot, off, blend_off()), "{off}");
    assert!(
        !bevels(&m, &s, &foot, off * (1.0 + LIMIT_TOLERANCE), blend_off()),
        "{off}"
    );
    let default = largest_chamfer_distance_all(&m, &s, &specs, 20.0, T)
        .unwrap()
        .unwrap();
    assert_eq!(
        largest_chamfer_distance_all_overflow(
            &m,
            &s,
            &specs,
            20.0,
            ChamferCorner::Facet,
            Overflow::default(),
            T,
        )
        .unwrap(),
        Some(default)
    );
}

/// A 2 mm plate's top front edge bevels past its seat as a notch over
/// the front's foot, to 22, which only trim through builds. With trim
/// through off the notch is refused, and the reading stops at the seat,
/// the plate's thickness.
#[test]
fn a_plate_s_bevel_stops_at_its_seat_with_trim_through_off() {
    let mut m = Model::new();
    let s = ogeom_algo::make_box(&mut m, at(0.0, 0.0, 0.0), (20.0, 20.0, 2.0), T)
        .unwrap()
        .shape;
    let edge = edges_where(&m, &s, |p| p.y.abs() < 1e-9 && (p.z - 2.0).abs() < 1e-9);
    let mut trim_off = Overflow::default();
    trim_off.trim_through = false;
    let read = |overflow: Overflow| {
        largest_chamfer_distance_overflow(
            &m,
            &s,
            &edge[0],
            &Chamfer::Symmetric(1.0),
            30.0,
            ChamferCorner::Facet,
            overflow,
            T,
        )
        .unwrap()
        .unwrap()
    };
    let notched = read(Overflow::default());
    assert!(
        (22.0 * (1.0 - 2.0 * LIMIT_TOLERANCE)..=22.0).contains(&notched),
        "{notched}"
    );
    let seated = read(trim_off);
    assert!(
        (2.0 * (1.0 - 2.0 * LIMIT_TOLERANCE)..=2.0).contains(&seated),
        "{seated}"
    );
    assert!(bevels(&m, &s, &edge, seated, trim_off));
    assert!(!bevels(
        &m,
        &s,
        &edge,
        seated * (1.0 + LIMIT_TOLERANCE),
        trim_off
    ));
}

/// Whether `size` builds and the next size up, a [`LIMIT_TOLERANCE`]
/// past it, and 1.01 times it are refused, by `builds`.
fn stops_at(size: f64, builds: impl Fn(f64) -> bool) -> bool {
    builds(size) && !builds(size * (1.0 + LIMIT_TOLERANCE)) && !builds(size * 1.01)
}

/// The foot round with rolling corners: the rim has no vertex where
/// rounds meet, so with blend into round off it stops where the setback
/// round does, about 1.0100, just past the seat 1.0 where its contact on
/// the top reaches the edge round's tangent line (the rim is 6 from the
/// boss's axis, the tangent line 7). The default reading, 8.62, is refused
/// under those options; under the default options the reading is the
/// plain one, probe for probe.
#[test]
#[ignore = "heavy"]
fn rolling_corners_read_their_limit_with_blend_into_round_off() {
    let (m, s, foot) = bossed();
    let read = |overflow: Overflow| {
        let sizing = BlendSizing::RoundsCornerOverflow {
            edges: &foot,
            corner: FilletCorner::Rolling,
            overflow,
        };
        blend_limit(&m, &s, &sizing, 20.0, T).unwrap()
    };
    let rolls = |r: f64| {
        let mut scratch = m.clone();
        fillet_edges_corner_overflow(
            &mut scratch,
            &s,
            &foot,
            r,
            FilletCorner::Rolling,
            blend_off(),
            T,
        )
        .is_ok()
    };
    let plain = BlendSizing::Rounds {
        edges: &foot,
        corner: FilletCorner::Rolling,
    };
    let plain = blend_limit(&m, &s, &plain, 20.0, T).unwrap();
    assert_eq!(read(Overflow::default()), plain);
    assert!(!rolls(plain.size.unwrap()), "{plain:?}");

    let found = read(blend_off());
    let size = found.size.unwrap();
    assert!(
        (1.0..=1.0100 * (1.0 + LIMIT_TOLERANCE)).contains(&size),
        "{found:?}"
    );
    assert!(
        (size - 1.0100).abs() <= 1.0100 * LIMIT_TOLERANCE,
        "{found:?}"
    );
    assert!(stops_at(size, rolls), "{found:?}");
}

/// A law holding one radius on the closed foot rim, its start freed (and
/// so its end, the same point): a level law is the constant round, so with
/// blend into round off it stops where that round does, about 1.0100.
/// The default reading is refused under those options; under the default
/// options the reading is the plain one, probe for probe.
#[test]
#[ignore = "heavy"]
fn a_level_law_reads_its_limit_with_blend_into_round_off() {
    let (m, s, foot) = bossed();
    let law = RadiusLaw::Linear {
        start: 0.5,
        end: 0.5,
    };
    let size = LawSize::At(LawStation::Start);
    let read = |overflow: Overflow| {
        let sizing = BlendSizing::RoundLawOverflow {
            chain: &foot,
            law: &law,
            size,
            overflow,
        };
        blend_limit(&m, &s, &sizing, 20.0, T).unwrap()
    };
    let rounds = |r: f64| {
        let at = RadiusLaw::Linear { start: r, end: r };
        let mut scratch = m.clone();
        fillet_edges_variable_overflow(&mut scratch, &s, &foot, &at, blend_off(), T).is_ok()
    };
    let plain = BlendSizing::RoundLaw {
        chain: &foot,
        law: &law,
        size,
    };
    let plain = blend_limit(&m, &s, &plain, 20.0, T).unwrap();
    assert_eq!(read(Overflow::default()), plain);
    assert!(!rounds(plain.size.unwrap()), "{plain:?}");

    let found = read(blend_off());
    let size = found.size.unwrap();
    assert!(
        (size - 1.0100).abs() <= 1.0100 * LIMIT_TOLERANCE,
        "{found:?}"
    );
    assert!(stops_at(size, rounds), "{found:?}");
}

/// A boss of radius 5 on a plate, a pocket 0.5 square cut into its top
/// beside the rim, its far wall 4.8 from the axis. The rim's round of
/// distances `k` along the top and `0.7 k` down the wall trims through the
/// pocket by default, up to the bound. With trim through off it is refused
/// once its contact on the top, the circle of radius `5 - k`, reaches the
/// pocket's far corner: `k = 5 - hypot(4.8, 0.25)`.
#[test]
#[ignore = "heavy"]
fn two_distances_stop_at_a_pocket_with_trim_through_off() {
    let mut m = Model::new();
    let plate = ogeom_algo::make_box(&mut m, at(-20.0, -20.0, 0.0), (40.0, 40.0, 4.0), T)
        .unwrap()
        .shape;
    let boss = ogeom_algo::make_cylinder(&mut m, at(0.0, 0.0, 4.0), 5.0, 5.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(&mut m, &plate, &boss, T).unwrap().shape;
    let pocket = ogeom_algo::make_box(&mut m, at(4.3, -0.25, 6.0), (0.5, 0.5, 4.0), T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(&mut m, &s, &pocket, T).unwrap().shape;
    let rim = edges_where(&m, &s, |p| {
        (p.z - 9.0).abs() < 1e-6 && (p.x.hypot(p.y) - 5.0).abs() < 1e-6
    });
    let top = explore_unique(&m, &s, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let b = ogeom_algo::shape_bounds(&m, f, T).unwrap();
            [b.low().unwrap(), b.high().unwrap()]
                .iter()
                .all(|p| (p.z - 9.0).abs() < 1e-6)
        })
        .unwrap();
    let mut trim_off = Overflow::default();
    trim_off.trim_through = false;
    let read = |overflow: Overflow| {
        let sizing = BlendSizing::RoundDistancesOverflow {
            edge: &rim[0],
            face: &top,
            on_face: 1.0,
            on_other: 0.7,
            profile: RoundProfile::Elliptical,
            overflow,
        };
        blend_limit(&m, &s, &sizing, 4.0, T).unwrap()
    };
    let plain = BlendSizing::RoundDistances {
        edge: &rim[0],
        face: &top,
        on_face: 1.0,
        on_other: 0.7,
        profile: RoundProfile::Elliptical,
    };
    let plain = blend_limit(&m, &s, &plain, 4.0, T).unwrap();
    assert_eq!(read(Overflow::default()), plain);
    assert_eq!(plain.size, Some(4.0));

    let found = read(trim_off);
    let size = found.size.unwrap();
    let corner = 5.0 - 4.8_f64.hypot(0.25);
    assert!(
        size <= corner && corner <= found.refused.unwrap(),
        "{found:?} against {corner}"
    );
    let rounds = |k: f64| {
        let mut scratch = m.clone();
        fillet_edge_distances_overflow(
            &mut scratch,
            &s,
            &rim[0],
            &top,
            k,
            0.7 * k,
            RoundProfile::Elliptical,
            trim_off,
            T,
        )
        .is_ok()
    };
    assert!(stops_at(size, rounds), "{found:?}");
}
