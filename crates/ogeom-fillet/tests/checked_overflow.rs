//! Checked blends under overflow options
//! (`fillet_edges_corner_overflow_checked`, `chamfer_edges_overflow_checked`):
//! for each option set that refuses a feature, a part holding that feature
//! and edges well clear of it. The refusal names exactly the feature's
//! edges, by their indices in the request; the request without them builds
//! with the same options, and its volume matches the closed form of the
//! clear edges' blends; the feature's edges alone are refused.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom_core::Tolerances;
use ogeom_fillet::{
    BlendRefusal, Chamfer, ChamferCorner, FilletCorner, Overflow, chamfer_edges_overflow,
    chamfer_edges_overflow_checked, fillet_edges, fillet_edges_corner_overflow,
    fillet_edges_corner_overflow_checked,
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

fn edge_at(m: &Model, s: &Shape, on: impl Fn(Point) -> bool) -> Shape {
    let found = edges_where(m, s, on);
    assert_eq!(found.len(), 1, "one edge");
    found[0].clone()
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

fn with(f: impl FnOnce(&mut Overflow)) -> Overflow {
    let mut o = Overflow::default();
    f(&mut o);
    o
}

/// What a round of radius `r` takes off a convex edge between square
/// planes, per unit length.
fn round_section(r: f64) -> f64 {
    r * r * (1.0 - PI / 4.0)
}

/// What rounds of radius `r` take off three convex edges of lengths
/// `lengths` meeting square at a box's corner, closed by the ball's
/// sphere, their other ends flush on end faces: each band past the corner
/// cube, and the cube less the sphere's octant.
fn rounded_corner(lengths: [f64; 3], r: f64) -> f64 {
    lengths
        .iter()
        .map(|l| (l - r) * round_section(r))
        .sum::<f64>()
        + r.powi(3) * (1.0 - PI / 6.0)
}

/// What bevels `d` take off three convex edges of lengths `lengths`
/// meeting square at a box's corner, by inclusion and exclusion of the
/// three prisms `x + y < d` near the corner: each prism `l d²/2`, each pair
/// overlapping in `d³/3`, all three in `d³/4`. The facet takes `d³/12`
/// more: the simplex below the plane `x + y + z = 2d` that the point
/// corner keeps.
fn bevelled_corner(lengths: [f64; 3], d: f64, corner: ChamferCorner) -> f64 {
    let point = lengths.iter().sum::<f64>() * d * d / 2.0 - d.powi(3) + d.powi(3) / 4.0;
    match corner {
        ChamferCorner::Point => point,
        ChamferCorner::Facet => point + d.powi(3) / 12.0,
    }
}

/// The three edges of a box meeting at its corner `(x, y, z0)` with the
/// box running towards `-x`, `-y` and `+z`: along `x`, along `y`, along
/// `z`.
fn corner_edges(m: &Model, s: &Shape, (x, y, z0): (f64, f64, f64)) -> [Shape; 3] {
    [
        edge_at(m, s, |p| (p.y - y).abs() < 1e-9 && (p.z - z0).abs() < 1e-9),
        edge_at(m, s, |p| (p.x - x).abs() < 1e-9 && (p.z - z0).abs() < 1e-9),
        edge_at(m, s, |p| (p.x - x).abs() < 1e-9 && (p.y - y).abs() < 1e-9),
    ]
}

/// One request: the solid, the request's edges with the feature's edges
/// standing among the clear ones, the indices the feature's edges take,
/// and what the clear edges' blends take off.
struct Case {
    name: String,
    solid: Shape,
    edges: Vec<Shape>,
    refused: Vec<usize>,
    /// The edges the refused edges are refused together with: none where
    /// they are refused alone.
    refused_with: Vec<usize>,
    taken: f64,
}

/// The request `clear` with `feature` put after its first edge.
fn case(name: String, solid: Shape, clear: Vec<Shape>, feature: Vec<Shape>, taken: f64) -> Case {
    let mut edges = vec![clear[0].clone()];
    let refused = (1..=feature.len()).collect();
    edges.extend(feature);
    edges.extend(clear.into_iter().skip(1));
    Case {
        name,
        solid,
        edges,
        refused,
        refused_with: Vec::new(),
        taken,
    }
}

/// A request's build under its options, by indices.
type Build<'a> = dyn Fn(&mut Model, &[usize]) -> ogeom_core::OgeomResult<ogeom_algo::Built> + 'a;

/// The checks every case shares, `refusal` being the checked call's:
/// exactly the feature's edges named and the unchecked call's error; the
/// request without them building to the closed form within `1e-9` of the
/// solid's volume; the feature's edges refused alone, or with the edges
/// they are refused together with. The measured error.
fn holds(m: &mut Model, case: &Case, refusal: &BlendRefusal, build: &Build<'_>) -> f64 {
    let all: Vec<usize> = (0..case.edges.len()).collect();
    let unchecked = build(&mut m.clone(), &all).unwrap_err();
    assert_eq!(refusal.edges, case.refused, "{refusal}");
    assert_eq!(refusal.error.to_string(), unchecked.to_string());
    let kept: Vec<usize> = all
        .iter()
        .copied()
        .filter(|i| !refusal.edges.contains(i))
        .collect();
    let before = volume(m, &case.solid);
    let built = build(m, &kept).unwrap();
    assert!(ogeom_algo::check(m, &built.shape, T).unwrap().is_valid());
    let got = before - volume(m, &built.shape);
    let error = (got - case.taken).abs();
    assert!(
        error <= 1e-9 * before,
        "took {got} against the closed form {}",
        case.taken
    );
    let mut together = case.refused_with.clone();
    together.extend(&case.refused);
    assert!(build(&mut m.clone(), &together).is_err());
    eprintln!(
        "{}: named {:?}, taken {got:.9}, closed form {:.9}, error {error:.1e}",
        case.name, refusal.edges, case.taken
    );
    error
}

fn round_holds(m: &mut Model, case: &Case, r: f64, corner: FilletCorner, overflow: Overflow) {
    let refusal =
        fillet_edges_corner_overflow_checked(m, &case.solid, &case.edges, r, corner, overflow, T)
            .unwrap_err();
    let build = |m: &mut Model, set: &[usize]| {
        let picked: Vec<Shape> = set.iter().map(|&i| case.edges[i].clone()).collect();
        fillet_edges_corner_overflow(m, &case.solid, &picked, r, corner, overflow, T)
    };
    holds(m, case, &refusal, &build);
}

fn bevel_holds(
    m: &mut Model,
    case: &Case,
    specs: &[(Shape, Chamfer)],
    corner: ChamferCorner,
    overflow: Overflow,
) {
    let refusal =
        chamfer_edges_overflow_checked(m, &case.solid, specs, corner, overflow, T).unwrap_err();
    let build = |m: &mut Model, set: &[usize]| {
        let picked: Vec<(Shape, Chamfer)> = set.iter().map(|&i| specs[i].clone()).collect();
        chamfer_edges_overflow(m, &case.solid, &picked, corner, overflow, T)
    };
    holds(m, case, &refusal, &build);
}

/// The solid between the plane sections `a` at `x = 0` and `b` at
/// `x = 20`, each a ring of `(y, z)` points, side faces planes.
fn prism(m: &mut Model, a: &[(f64, f64)], b: &[(f64, f64)]) -> Shape {
    let n = a.len();
    let mut points: Vec<Point> = a.iter().map(|(y, z)| Point::new(0.0, *y, *z)).collect();
    points.extend(b.iter().map(|(y, z)| Point::new(20.0, *y, *z)));
    for flip in [false, true] {
        let mut ends = [(0..n).collect::<Vec<_>>(), (n..2 * n).rev().collect()];
        let mut sides: Vec<Vec<usize>> = (0..n)
            .map(|i| vec![i, n + i, n + (i + 1) % n, (i + 1) % n])
            .collect();
        if flip {
            ends.iter_mut().for_each(|r| r.reverse());
            sides.iter_mut().for_each(|r| r.reverse());
        }
        let rings: Vec<Vec<usize>> = ends.into_iter().chain(sides).collect();
        if let Ok((built, _)) = ogeom_algo::make_faceted_solid(m, &points, &rings, &[], T) {
            return built.shape;
        }
    }
    panic!("the rings enclose no volume either way round");
}

/// The block 40 by 22 by 15, its top edges along `y = 0` and `y = 22`
/// rounded at 5, a boss of radius 6 and 15 tall standing on the top about
/// `(20, 11)`; with the boss's foot rim, and the block's corner edges at
/// `(40, 22, 0)`, the upright one rising 10 to the back round's tangent
/// line. The rim is two arcs, front and back, and a round of 3 or a bevel
/// of 3 on either reaches 9 from the axis, past the edge round's tangent
/// line 6 away on its side: each arc reaches an existing round, as the
/// foot of `overflow.rs` and `chamfer_overflow.rs` does.
fn bossed(m: &mut Model) -> (Shape, Vec<Shape>, [Shape; 3]) {
    boss_on(m, 22.0, true)
}

/// The block 40 by `deep` by 15, its top front edge rounded at 5 (and its
/// top back edge with `back`), the boss standing midway across it or at
/// `y = 12` on the deep block; the foot rim and the corner edges at
/// `(40, deep, 0)`.
fn boss_on(m: &mut Model, deep: f64, back: bool) -> (Shape, Vec<Shape>, [Shape; 3]) {
    let block = ogeom_algo::make_box(m, at(0.0, 0.0, 0.0), (40.0, deep, 15.0), T)
        .unwrap()
        .shape;
    let top = edges_where(m, &block, |p| {
        (p.y.abs() < 1e-9 || (back && (p.y - deep).abs() < 1e-9)) && (p.z - 15.0).abs() < 1e-9
    });
    let rounded = fillet_edges(m, &block, &top, 5.0, T).unwrap().shape;
    let y = if back { deep / 2.0 } else { 12.0 };
    let boss = ogeom_algo::make_cylinder(m, at(20.0, y, 15.0), 6.0, 15.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(m, &rounded, &boss, T).unwrap().shape;
    let foot = edges_where(m, &s, |p| {
        (p.z - 15.0).abs() < 1e-9 && ((p.x - 20.0).hypot(p.y - y) - 6.0).abs() < 1e-9
    });
    assert_eq!(foot.len(), 2, "the foot rim's two arcs");
    let corner = corner_edges(m, &s, (40.0, deep, 0.0));
    (s, foot, corner)
}

const BOSSED_CORNER: [f64; 3] = [40.0, 22.0, 10.0];

/// The foot round with blend into round off (trim through's own refusal),
/// with trim through off too, and with every option off; the corner's
/// three rounds close on the ball's sphere. (A rolling corner is refused
/// on a solid with a curved face: the planar parts below take it.)
#[test]
fn a_foot_round_refused_by_its_options_is_named() {
    for overflow in [
        with(|o| o.blend_into_round = false),
        with(|o| {
            o.blend_into_round = false;
            o.trim_through = false;
        }),
        Overflow::none(),
    ] {
        let mut m = Model::new();
        let (s, foot, clear) = bossed(&mut m);
        let taken = rounded_corner(BOSSED_CORNER, 3.0);
        let name = format!("foot round, {overflow:?}");
        let case = case(name, s, clear.to_vec(), foot, taken);
        round_holds(&mut m, &case, 3.0, FilletCorner::Setback, overflow);
    }
}

/// The foot bevel with every option off and with blend into round and
/// trim through off, under either corner: the corner's three bevels close
/// on the point or the facet. The foot of `chamfer_overflow.rs`, one edge
/// round: a whole rim reaching it is refused by name with these options,
/// where a rim reaching two rounds, or one arc of it, is not yet.
#[test]
fn a_foot_bevel_refused_by_its_options_is_named() {
    for corner in [ChamferCorner::Point, ChamferCorner::Facet] {
        for overflow in [
            Overflow::none(),
            with(|o| {
                o.blend_into_round = false;
                o.trim_through = false;
            }),
        ] {
            let mut m = Model::new();
            let (s, foot, clear) = boss_on(&mut m, 30.0, false);
            let taken = bevelled_corner([40.0, 30.0, 15.0], 2.0, corner);
            let name = format!("foot bevel, {corner:?}, {overflow:?}");
            let case = case(name, s, clear.to_vec(), foot, taken);
            let specs: Vec<(Shape, Chamfer)> = case
                .edges
                .iter()
                .enumerate()
                .map(|(i, e)| {
                    let d = if case.refused.contains(&i) { 3.0 } else { 2.0 };
                    (e.clone(), Chamfer::Symmetric(d))
                })
                .collect();
            bevel_holds(&mut m, &case, &specs, corner, overflow);
        }
    }
}

/// A block 30 by 20 by 15 under `z = 0`, a pocket of radius 5 and depth 8
/// sunk into its top about `(15, 7)`: the front top edge's round of 4
/// rolls over the pocket's rim, which hold edge keeps (as `hold_edge.rs`
/// builds it). With the corner edges at `(30, 20, -15)`.
fn pocketed(m: &mut Model) -> (Shape, Shape, [Shape; 3]) {
    let block = ogeom_algo::make_box(m, at(0.0, 0.0, -15.0), (30.0, 20.0, 15.0), T)
        .unwrap()
        .shape;
    let pocket = ogeom_algo::make_cylinder(m, at(15.0, 7.0, -8.0), 5.0, 9.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(m, &block, &pocket, T).unwrap().shape;
    let front = edge_at(m, &s, |p| p.y.abs() < 1e-9 && p.z.abs() < 1e-9);
    let corner = corner_edges(m, &s, (30.0, 20.0, -15.0));
    (s, front, corner)
}

/// Hold edge off with trim through off, and every option off.
#[test]
fn a_round_over_a_pocket_rim_refused_by_its_options_is_named() {
    for overflow in [
        with(|o| {
            o.hold_edge = false;
            o.trim_through = false;
        }),
        Overflow::none(),
    ] {
        let mut m = Model::new();
        let (s, front, clear) = pocketed(&mut m);
        let taken = rounded_corner([30.0, 20.0, 15.0], 4.0);
        let name = format!("pocket rim, {overflow:?}");
        let case = case(name, s, clear.to_vec(), vec![front], taken);
        round_holds(&mut m, &case, 4.0, FilletCorner::Setback, overflow);
    }
}

/// The loft of `hold_edge.rs`: a block 20 long whose top narrows from 6
/// deep to 2, a cliff 6 high dropping from its back edge to a ledge, 12
/// deep and 15 tall below. The front top edge's band end runs over the
/// cliff's edge, which hold edge at the end keeps. With the three edges at
/// its corner `(20, 12, -15)`, square between planes: along `x` 20 long,
/// along `y` 12, and up to the ledge 9. Every face is a plane, so a
/// rolling corner closes that corner.
fn cliffed(m: &mut Model) -> (Shape, Shape, Vec<Shape>) {
    let end = |back: f64| {
        [
            (0.0, 0.0),
            (0.0, -15.0),
            (12.0, -15.0),
            (12.0, -6.0),
            (back, -6.0),
            (back, 0.0),
        ]
    };
    let s = prism(m, &end(6.0), &end(2.0));
    let front = edge_at(m, &s, |p| p.y.abs() < 1e-9 && p.z.abs() < 1e-9);
    let clear = corner_edges(m, &s, (20.0, 12.0, -15.0)).to_vec();
    (s, front, clear)
}

/// Hold edge at the end off with trim through off, hold edge off with trim
/// through off, and every option off; setback and rolling corners, the
/// corner the ball's sphere under both. Under rolling corners a request
/// with a vertex where rounds meet is the engine's, which refuses a round
/// reaching a face beyond its ends whatever the options; the front edge
/// is refused with the corner's first edge by the options, as under
/// setback corners.
#[test]
fn a_band_end_over_a_cliff_refused_by_its_options_is_named() {
    for corner in [FilletCorner::Setback, FilletCorner::Rolling] {
        for overflow in [
            with(|o| {
                o.hold_edge_at_end = false;
                o.trim_through = false;
            }),
            with(|o| {
                o.hold_edge = false;
                o.trim_through = false;
            }),
            Overflow::none(),
        ] {
            let mut m = Model::new();
            let (s, front, clear) = cliffed(&mut m);
            let taken = rounded_corner([20.0, 12.0, 9.0], 4.0);
            let name = format!("cliff, {corner:?}, {overflow:?}");
            let case = case(name, s, clear, vec![front], taken);
            round_holds(&mut m, &case, 4.0, corner, overflow);
        }
    }
}

/// A plate 40 by 40 by 4 with a boss of radius 5 and 5 tall, a pocket
/// sunk into the boss beside its top rim: the rim's round of 1 reaches the
/// pocket's walls and trims through them (as `overflow.rs` builds it).
/// With the plate's corner edges at `(20, 20, 0)`.
fn plate_rim(m: &mut Model) -> (Shape, Vec<Shape>, [Shape; 3]) {
    let plate = ogeom_algo::make_box(m, at(-20.0, -20.0, 0.0), (40.0, 40.0, 4.0), T)
        .unwrap()
        .shape;
    let boss = ogeom_algo::make_cylinder(m, at(0.0, 0.0, 4.0), 5.0, 5.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::fuse(m, &plate, &boss, T).unwrap().shape;
    let pocket = ogeom_algo::make_box(m, at(4.3, -0.25, 6.0), (0.5, 0.5, 4.0), T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(m, &s, &pocket, T).unwrap().shape;
    let rim = edges_where(m, &s, |p| {
        (p.z - 9.0).abs() < 1e-6 && (p.x.hypot(p.y) - 5.0).abs() < 1e-6
    });
    assert!(!rim.is_empty(), "the rim");
    let corner = corner_edges(m, &s, (20.0, 20.0, 0.0));
    (s, rim, corner)
}

/// Trim through off, and every option off.
#[test]
fn a_rim_round_trimming_through_refused_by_its_options_is_named() {
    for overflow in [with(|o| o.trim_through = false), Overflow::none()] {
        let mut m = Model::new();
        let (s, rim, clear) = plate_rim(&mut m);
        let taken = rounded_corner([40.0, 40.0, 4.0], 1.0);
        let name = format!("boss rim, {overflow:?}");
        let case = case(name, s, clear.to_vec(), rim, taken);
        round_holds(&mut m, &case, 1.0, FilletCorner::Setback, overflow);
    }
}

/// A block 20 long with a through slot 6 wide and 8 deep from `y = 10`,
/// its floor at 7 (as `narrow_passage.rs` builds it): the floor's rim
/// rounded at 4 is a narrow passage. Either rim edge's round alone is a
/// plain concave round, but the two bands, each 4 wide on a floor 6 wide,
/// overlap: only the later is refused. With the bottom front and bottom
/// back edges, each 20 long between square planes.
fn slot(m: &mut Model) -> (Shape, Vec<Shape>, Vec<Shape>) {
    let section = [
        (0.0, 0.0),
        (30.0, 0.0),
        (30.0, 15.0),
        (16.0, 15.0),
        (16.0, 7.0),
        (10.0, 7.0),
        (10.0, 15.0),
        (0.0, 15.0),
    ];
    let s = prism(m, &section, &section);
    let rim = edges_where(m, &s, |p| {
        (p.z - 7.0).abs() < 1e-9 && ((p.y - 10.0).abs() < 1e-9 || (p.y - 16.0).abs() < 1e-9)
    });
    assert_eq!(rim.len(), 2, "the floor's rim along its walls");
    let clear = vec![
        edge_at(m, &s, |p| p.y.abs() < 1e-9 && p.z.abs() < 1e-9),
        edge_at(m, &s, |p| (p.y - 30.0).abs() < 1e-9 && p.z.abs() < 1e-9),
    ];
    (s, rim, clear)
}

/// Every narrow-passage strategy off, and every option off.
#[test]
fn a_narrow_passage_refused_by_its_options_is_named() {
    let narrow_off = with(|o| {
        o.narrow_blend_made = false;
        o.narrow_whole_face = false;
        o.narrow_partial_face = false;
    });
    for overflow in [narrow_off, Overflow::none()] {
        let mut m = Model::new();
        let (s, rim, clear) = slot(&mut m);
        // The two clear rounds, less the first rim edge's concave fill.
        let taken = 2.0 * 20.0 * round_section(4.0) - 20.0 * round_section(4.0);
        let name = format!("slot, {overflow:?}");
        let mut case = case(name, s, clear, rim, taken);
        case.refused = vec![2];
        case.refused_with = vec![1];
        round_holds(&mut m, &case, 4.0, FilletCorner::Setback, overflow);
    }
}

/// A plate 10 by 10 by 2, its top front edge bevelled 1 on the top and 3
/// down the 2 mm front, a notch trimming through the floor (as
/// `chamfer_overflow.rs` builds it), refused with trim through off; the
/// plate's top back edge bevelled 0.5, under either corner. The clear edge
/// owns none of the faces the notch trims (ADR 0001, amendment "checked
/// blends under overflow", on a face another edge of the request owns).
#[test]
fn a_notch_refused_with_trim_through_off_is_named() {
    for corner in [ChamferCorner::Point, ChamferCorner::Facet] {
        for overflow in [with(|o| o.trim_through = false), Overflow::none()] {
            let mut m = Model::new();
            let s = ogeom_algo::make_box(&mut m, at(0.0, 0.0, 0.0), (10.0, 10.0, 2.0), T)
                .unwrap()
                .shape;
            let top = explore_unique(&m, &s, ShapeType::Face)
                .unwrap()
                .into_iter()
                .find(|f| ogeom_algo::face_normal(&m, f, T).unwrap().1.z > 0.5)
                .unwrap();
            let notch = edge_at(&m, &s, |p| p.y.abs() < 1e-9 && (p.z - 2.0).abs() < 1e-9);
            let back = edge_at(&m, &s, |p| {
                (p.y - 10.0).abs() < 1e-9 && (p.z - 2.0).abs() < 1e-9
            });
            let name = format!("notch, {corner:?}, {overflow:?}");
            let case = case(name, s, vec![back], vec![notch], 10.0 * 0.5 * 0.5 / 2.0);
            let specs: Vec<(Shape, Chamfer)> = case
                .edges
                .iter()
                .enumerate()
                .map(|(i, e)| {
                    let spec = if case.refused.contains(&i) {
                        Chamfer::Distances {
                            face: top.clone(),
                            on_face: 1.0,
                            on_other: 3.0,
                        }
                    } else {
                        Chamfer::Symmetric(0.5)
                    };
                    (e.clone(), spec)
                })
                .collect();
            bevel_holds(&mut m, &case, &specs, corner, overflow);
        }
    }
}

/// A request that builds is the unchecked call's result, face for face
/// and to the volume, under options and either corner; an empty request
/// names nothing.
#[test]
fn a_request_that_builds_is_the_unchecked_result() {
    let faces = |m: &Model, s: &Shape| explore_unique(m, s, ShapeType::Face).unwrap().len();
    let same = |m: &Model, a: &ogeom_algo::Built, b: &ogeom_algo::Built| {
        assert_eq!(faces(m, &a.shape), faces(m, &b.shape));
        assert!((volume(m, &a.shape) - volume(m, &b.shape)).abs() < 1e-9);
    };
    let mut m = Model::new();
    // The foot of `overflow.rs`, one edge round: blend into round builds
    // it with trim through off (it takes a request of the foot alone).
    let (s, foot, clear) = boss_on(&mut m, 30.0, false);
    let edges = foot;
    let trim_off = with(|o| o.trim_through = false);
    let setback = FilletCorner::Setback;
    let a = fillet_edges_corner_overflow_checked(&mut m, &s, &edges, 3.0, setback, trim_off, T)
        .unwrap();
    let b = fillet_edges_corner_overflow(&mut m, &s, &edges, 3.0, setback, trim_off, T).unwrap();
    same(&m, &a, &b);
    let specs: Vec<(Shape, Chamfer)> = clear
        .iter()
        .map(|e| (e.clone(), Chamfer::Symmetric(2.0)))
        .collect();
    for corner in [ChamferCorner::Point, ChamferCorner::Facet] {
        let none = Overflow::none();
        let a = chamfer_edges_overflow_checked(&mut m, &s, &specs, corner, none, T).unwrap();
        let b = chamfer_edges_overflow(&mut m, &s, &specs, corner, none, T).unwrap();
        same(&m, &a, &b);
    }
    // A rolling corner on the planar loft, every option off.
    let (loft, _, edges) = cliffed(&mut m);
    let (rolling, none) = (FilletCorner::Rolling, Overflow::none());
    let a =
        fillet_edges_corner_overflow_checked(&mut m, &loft, &edges, 4.0, rolling, none, T).unwrap();
    let b = fillet_edges_corner_overflow(&mut m, &loft, &edges, 4.0, rolling, none, T).unwrap();
    same(&m, &a, &b);

    let empty =
        fillet_edges_corner_overflow_checked(&mut m, &s, &[], 3.0, rolling, none, T).unwrap_err();
    assert!(empty.edges.is_empty());
    let empty =
        chamfer_edges_overflow_checked(&mut m, &s, &[], ChamferCorner::Point, none, T).unwrap_err();
    assert!(empty.edges.is_empty());
}
