//! A face split along a curve that crosses onto it from a sibling face of
//! the same surface, measured against the closed-form sectors of a drum's
//! wall.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    SheetMethod, SheetPlacement, SurfaceSheet, check, face_normal, make_box, make_cylinder,
    surface_properties,
};
use ogeom::boolean::cut;
use ogeom::core::Tolerances;
use ogeom::geom::{Curve2d as _, Line2d, PlanarCurve};
use ogeom::heal::{Projection, split_face};
use ogeom::math::{Direction, Frame, Point, Point2, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn faces(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Face).unwrap()
}

fn area(model: &Model, face: &Shape) -> f64 {
    let fine = Deflection {
        chord: 1e-4,
        angular: 0.01,
        ..Deflection::default()
    };
    surface_properties(model, face, fine, T).unwrap().mass
}

#[test]
fn a_line_crossing_onto_the_lower_half_of_a_slab_cut_drum_splits_it() {
    // A drum of radius 10 and height 10 cut through by a slab 1 thick
    // across y = 0: its wall comes back as two halves of the one cylinder.
    // A sheet laid at the top of the +y half, +x round the drum toward -x,
    // carries the line from (-4, -5) to (20, 5): from the bottom rim at
    // 90 - 22.9 degrees up and round, across the slab near 180 degrees, to
    // the top rim at 90 + 114.6 degrees on the -y half.
    let mut model = Model::new();
    let drum = make_cylinder(&mut model, Frame::WORLD, 10.0, 10.0, T)
        .unwrap()
        .shape;
    let slab = make_box(
        &mut model,
        Frame::new(Point::new(-15.0, -0.5, -1.0), Direction::Z, Direction::X, T).unwrap(),
        (30.0, 1.0, 12.0),
        T,
    )
    .unwrap()
    .shape;
    let halves = cut(&mut model, &drum, &slab, T).unwrap().shape;
    let walls: Vec<Shape> = faces(&model, &halves)
        .into_iter()
        .filter(|f| {
            face_normal(&model, f, T)
                .is_ok_and(|(p, n)| n.z.abs() < 0.5 && (p.x.hypot(p.y) - 10.0).abs() < 1e-6)
        })
        .collect();
    assert_eq!(walls.len(), 2, "the wall in two halves");
    let lower = walls
        .iter()
        .find(|f| face_normal(&model, f, T).unwrap().0.y < 0.0)
        .unwrap()
        .clone();
    let upper = walls
        .iter()
        .find(|f| f.node() != lower.node())
        .unwrap()
        .clone();
    let sheet = SurfaceSheet::new(
        &model,
        &upper,
        &SheetPlacement {
            at: Point::new(0.0, 10.0, 5.0),
            x_direction: Direction::new(Vector::new(-1.0, 0.0, 0.0), T).unwrap(),
            method: SheetMethod::Geodesic,
        },
        30.0,
        T,
    )
    .unwrap();
    let (from, to) = (Point2::new(-4.0, -5.0), Point2::new(20.0, 5.0));
    let line: PlanarCurve = Line2d::segment(from, to, T).unwrap().into();
    let range = line.domain();
    let edge = sheet
        .chart_edge(&mut model, &line, range, 1e-6)
        .unwrap()
        .shape;

    // On the development, arc length round the wall from the sheet's origin
    // is s.x, and the height is s.y + 5. The lower half starts at
    // s0 = 10 (pi / 2 + asin(1 / 20)); above the line, between s0 and the
    // line's end at s.x = 20, lies a triangle of legs 20 - s0 along the top
    // rim and (20 - s0) 10 / 24 down the slab's edge.
    let s0 = 10.0 * (core::f64::consts::FRAC_PI_2 + (0.05_f64).asin());
    let triangle = (to.x - s0).powi(2) * (to.y - from.y) / (to.x - from.x) / 2.0;
    let whole = area(&model, &lower);
    let wall = 10.0 * 10.0 * (core::f64::consts::PI - 2.0 * (0.05_f64).asin());
    assert!((whole - wall).abs() <= 1e-6 * wall, "{whole} {wall}");

    for projection in [Projection::AlongNormals, Projection::OnFace] {
        let mut model = model.clone();
        let built = split_face(
            &mut model,
            &halves,
            &lower,
            core::slice::from_ref(&edge),
            projection,
            T,
        )
        .unwrap_or_else(|e| panic!("{projection:?}: {e}"));
        let pieces = built.history.trace(&lower).to_vec();
        assert_eq!(pieces.len(), 2, "{projection:?}: the half in two");
        let mut areas: Vec<f64> = pieces.iter().map(|p| area(&model, p)).collect();
        areas.sort_by(f64::total_cmp);
        eprintln!(
            "{projection:?}: {areas:?} against {triangle} and {}",
            wall - triangle
        );
        assert!(
            (areas[0] - triangle).abs() <= 1e-5 * triangle,
            "{projection:?}: {} {triangle}",
            areas[0]
        );
        assert!(
            (areas[1] - (wall - triangle)).abs() <= 1e-6 * wall,
            "{projection:?}: {} {}",
            areas[1],
            wall - triangle
        );
        let report = check(&model, &built.shape, T).unwrap();
        assert!(
            report.problems.is_empty(),
            "{projection:?}: {:?}",
            report.problems
        );
    }
}
