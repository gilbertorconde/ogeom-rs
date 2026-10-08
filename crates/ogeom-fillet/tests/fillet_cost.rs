//! What one fillet call costs on a solid with many faces away from it, and
//! that rounding many edges of a holed plate at once gives the closed-form
//! solid.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_core::clock::Instant;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};
use std::time::Duration;

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// A `side` x `side` x 10 plate with an `n` x `n` grid of through holes of
/// radius 1.5 at `pitch`, the grid centred on the plate; the plain plate
/// for `n` = 0.
fn plate(model: &mut Model, side: f64, n: u32, pitch: f64) -> Shape {
    let first = (side - pitch * f64::from(n.max(1) - 1)) / 2.0;
    let mut centres = Vec::new();
    for i in 0..n {
        for j in 0..n {
            centres.push((first + pitch * f64::from(i), first + pitch * f64::from(j)));
        }
    }
    drilled(model, side, &centres)
}

/// A `side` x `side` x 10 plate with through holes of radius 1.5 at
/// `centres`.
fn drilled(model: &mut Model, side: f64, centres: &[(f64, f64)]) -> Shape {
    let plate = ogeom_algo::make_box(model, Frame::WORLD, (side, side, 10.0), T)
        .unwrap()
        .shape;
    if centres.is_empty() {
        return plate;
    }
    let mut drums = Vec::new();
    for &(x, y) in centres {
        let at = Frame::new(Point::new(x, y, -1.0), Direction::Z, Direction::X, T).unwrap();
        drums.push(
            ogeom_algo::make_cylinder(model, at, 1.5, 12.0, T)
                .unwrap()
                .shape,
        );
    }
    let grid = model.add_compound(&drums).unwrap();
    ogeom_bool::cut(model, &plate, &grid, T).unwrap().shape
}

/// The edges of the top face: its outer edges, or its hole rims.
fn top_edges(model: &Model, solid: &Shape, rims: bool) -> Vec<Shape> {
    use ogeom_geom::Curve3d as _;
    explore_unique(model, solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let data = model.node(e).unwrap().data().as_edge().unwrap();
            let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                return false;
            };
            let curve = model.geometry().curve(*curve).unwrap();
            let mid = curve.point_at(f64::midpoint(range.0, range.1), T).unwrap();
            let circle = matches!(curve, ogeom_geom::Curve::Circle(_));
            mid.z > 9.99 && circle == rims
        })
        .collect()
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom_algo::volume_properties(
        model,
        shape,
        ogeom_mesh::Deflection {
            chord: 1e-3,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

/// The corner a fillet of radius `r` takes off a right-angled edge: its
/// area, and how far its centroid stands from either face.
fn corner(r: f64) -> (f64, f64) {
    (
        r * r * (1.0 - PI / 4.0),
        r * (10.0 - 3.0 * PI) / (3.0 * (4.0 - PI)),
    )
}

/// Fillet `edges` of `solid` at `radius`: the time it took, the volume it
/// removed, and the result's face count. The result must be valid.
fn fillet(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    radius: f64,
) -> (Duration, f64, usize) {
    let start = Instant::now();
    let rounded = ogeom_fillet::fillet_edges(model, solid, edges, radius, T).unwrap();
    let took = start.elapsed();
    assert!(
        ogeom_algo::check(model, &rounded.shape, T)
            .unwrap()
            .is_valid()
    );
    let faces = explore_unique(model, &rounded.shape, ShapeType::Face)
        .unwrap()
        .len();
    let removed = volume(model, solid) - volume(model, &rounded.shape);
    (took, removed, faces)
}

/// Rims on a checkerboard, each blend's box reaching into its diagonal
/// neighbours' though the blends stay apart, round in one call: the
/// blends whose boxes meet go in different booleans. Each blend removes
/// the corner between the plate and its bore turned about the bore's axis
/// (Pappus).
#[test]
fn rims_whose_blend_boxes_meet_round_together() {
    let mut model = Model::new();
    // Blends of radius 0.5 on holes of radius 1.5 reach 2 out: diagonal
    // neighbours 3.5 apart each way have meeting boxes and blends about
    // 0.95 apart.
    let mut centres = Vec::new();
    for i in 0..4 {
        for j in 0..4 {
            if (i + j) % 2 == 0 {
                centres.push((
                    3.5f64.mul_add(f64::from(i), 2.25),
                    3.5f64.mul_add(f64::from(j), 2.25),
                ));
            }
        }
    }
    let holed = drilled(&mut model, 15.0, &centres);
    let rims = top_edges(&model, &holed, true);
    let (_, removed, faces) = fillet(&mut model, &holed, &rims, 0.5);
    let (area, centroid) = corner(0.5);
    let want = 8.0 * area * 2.0 * PI * (1.5 + centroid);
    assert!(
        (removed - want).abs() < want * 1e-4,
        "{removed} against {want}"
    );
    // The plate's six faces, and per hole its wall and its blend.
    assert_eq!(faces, 6 + 8 * 2);
}

/// The outer edges of a small holed plate's top face meet at mitres, the
/// later blends at each corner running on through the earlier ones until
/// the ball leaves the plate.
#[test]
fn the_outer_edges_of_a_holed_plate_meet_at_mitres() {
    let mut model = Model::new();
    let holed = plate(&mut model, 20.0, 2, 10.0);
    let edges = top_edges(&model, &holed, false);
    assert_eq!(edges.len(), 4);
    let (_, removed, faces) = fillet(&mut model, &holed, &edges, 1.0);
    let (area, centroid) = corner(1.0);
    // Each blend's centroid line is its edge less the centroid's offset at
    // either mitre.
    let want = 4.0 * area * 2.0f64.mul_add(-centroid, 20.0);
    assert!(
        (removed - want).abs() < want * 1e-6,
        "{removed} against {want}"
    );
    assert_eq!(faces, 10 + 4);
}

/// The outer edges of a holed plate's top face fillet within twice the
/// time they take on the plain plate: the run-out probes at the corners
/// read the faces near the corners, and each pass's boolean sets aside
/// every face clear of each of its wedges, the holes among them.
#[test]
#[ignore = "heavy"]
fn a_plate_with_holes_fillets_its_outer_edges_at_about_the_plain_plates_cost() {
    let (area, centroid) = corner(1.0);
    // Each blend runs to the mitres with its neighbours: its centroid's
    // line is the edge less the centroid's offset at either end.
    let want = 4.0 * area * 2.0f64.mul_add(-centroid, 100.0);
    let mut plates = Vec::new();
    for n in [0, 4, 8, 16] {
        let mut model = Model::new();
        let solid = plate(&mut model, 100.0, n, 80.0 / f64::from(n.max(1)));
        let edges = top_edges(&model, &solid, false);
        assert_eq!(edges.len(), 4);
        let (took, removed, faces) = fillet(&mut model, &solid, &edges, 1.0);
        eprintln!("{n} x {n} holes: {took:?}, removed {removed}, {faces} faces");
        assert!(
            (removed - want).abs() < want * 1e-4,
            "{removed} against {want}"
        );
        assert_eq!(faces, 10 + (n * n) as usize);
        plates.push((model, solid, edges));
    }
    // The least of several runs of each side, taken in turns, reads past
    // a machine busy with other work; a round is measured again, up to
    // four, while the bound is not met.
    let least = |(model, solid, edges): &mut (Model, Shape, Vec<Shape>)| {
        (0..3).fold(Duration::MAX, |least, _| {
            let start = Instant::now();
            ogeom_fillet::fillet_edges(model, solid, edges, 1.0, T).unwrap();
            least.min(start.elapsed())
        })
    };
    let mut rounds = Vec::new();
    for _ in 0..4 {
        let (mut plain, mut holed) = (Duration::MAX, Duration::MAX);
        for _ in 0..4 {
            plain = plain.min(least(&mut plates[0]));
            holed = holed.min(least(&mut plates[3]));
        }
        rounds.push((holed, plain));
        if holed <= plain * 2 {
            return;
        }
    }
    panic!("256 holes against the plain plate, (holed, plain) per round: {rounds:?}");
}

/// Every hole rim on a plate's top face fillets in time that grows about
/// linearly with the number of rims. The one boolean that cuts all the
/// blends still grows a little faster than its tool's lumps, so the bound
/// per hole is loose.
#[test]
#[ignore = "heavy"]
fn every_rim_of_a_holed_plate_fillets_in_about_linear_time() {
    let (area, centroid) = corner(0.5);
    let mut per_hole = Vec::new();
    for n in [4_u32, 8, 16] {
        let mut model = Model::new();
        let holed = plate(&mut model, 100.0, n, 80.0 / f64::from(n));
        let rims = top_edges(&model, &holed, true);
        let (took, removed, faces) = fillet(&mut model, &holed, &rims, 0.5);
        let holes = f64::from(n * n);
        eprintln!(
            "{} rims: {took:?}, removed {removed}, {faces} faces",
            rims.len()
        );
        let want = holes * area * 2.0 * PI * (1.5 + centroid);
        assert!(
            (removed - want).abs() < want * 1e-4,
            "{removed} against {want}"
        );
        assert_eq!(faces, 6 + 2 * (n * n) as usize);
        per_hole.push(took.as_secs_f64() / holes);
    }
    assert!(
        per_hole[2] < per_hole[0] * 6.0,
        "per hole {} s at 16 x 16 against {} s at 4 x 4",
        per_hole[2],
        per_hole[0]
    );
}
