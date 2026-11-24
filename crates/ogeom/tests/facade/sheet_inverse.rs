//! A sheet's inverse map, from a surface point back to the sheet point
//! that lands on it, measured against the developments of a cylinder and a
//! cone, a sphere's great circles and a torus's geodesic equator and
//! meridians, and against the forward map both ways round; and a carried
//! edge's pcurves on every face of the sheet's surface.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom::algo::{
    SheetMethod, SheetPlacement, SurfaceSheet, check, make_box, make_cone, make_cylinder,
    make_natural_face, make_sphere, make_torus,
};
use ogeom::boolean::cut;
use ogeom::core::Tolerances;
use ogeom::geom::{Curve2d as _, Curve3d as _, Line2d, PlanarCurve, Surface as _};
use ogeom::heal::{Projection, split_face};
use ogeom::math::{Direction, Frame, Point, Point2, Vector};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The round trip asked of every sheet, on the sheet, in millimetres.
const ROUND_TRIP: f64 = 1e-9;

fn faces(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Face).unwrap()
}

/// The face of `shape` whose normal's `z` is under `below` in size: a
/// drum's or a cone's wall.
fn wall(model: &Model, shape: &Shape, below: f64) -> Shape {
    faces(model, shape)
        .into_iter()
        .find(|f| ogeom::algo::face_normal(model, f, T).is_ok_and(|(_, n)| n.z.abs() < below))
        .unwrap()
}

fn laid(
    model: &Model,
    face: &Shape,
    at: Point,
    x: Direction,
    method: SheetMethod,
    reach: f64,
) -> SurfaceSheet {
    SurfaceSheet::new(
        model,
        face,
        &SheetPlacement {
            at,
            x_direction: x,
            method,
        },
        reach,
        T,
    )
    .unwrap()
}

fn lift(sheet: &SurfaceSheet, s: Point2) -> Point {
    let uv = sheet.chart_point(s).unwrap();
    sheet.surface().point_at(uv.x, uv.y, T).unwrap()
}

/// Sheet points over the disc of radius `reach`: the origin and 24 rays
/// at eight radii out to the rim.
fn disc(reach: f64) -> Vec<Point2> {
    let mut out = vec![Point2::new(0.0, 0.0)];
    for k in 0..24 {
        let angle = f64::from(k) * PI / 12.0 + 0.1;
        for i in 1..=8 {
            let r = reach * f64::from(i) / 8.0;
            out.push(Point2::new(r * angle.cos(), r * angle.sin()));
        }
    }
    out
}

/// The worst miss, on the sheet, of `sheet_point(lift(s))` against `s`
/// over the disc of the sheet's reach, and the worst miss, in space, of
/// `lift(sheet_point(p))` against `p` over the same points landed.
fn round_trips(sheet: &SurfaceSheet, reach: f64) -> (f64, f64) {
    let (mut there, mut back): (f64, f64) = (0.0, 0.0);
    for s in disc(reach) {
        let p = lift(sheet, s);
        let found = sheet
            .sheet_point(p)
            .unwrap()
            .unwrap_or_else(|| panic!("{s:?} comes back past the reach"));
        there = there.max((found - s).magnitude());
        back = back.max(lift(sheet, found).distance(p));
    }
    (there, back)
}

#[test]
fn a_plane_comes_back_as_its_own_coordinates() {
    // The top of a 10 x 10 x 10 box, the sheet at its middle with +x along
    // x: the sheet's +y is the outward normal +z turned against +x, so the
    // point (x, y, 10) is the sheet point (x - 5, y - 5).
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let top = faces(&model, &block)
        .into_iter()
        .find(|f| ogeom::algo::face_normal(&model, f, T).is_ok_and(|(_, n)| n.z > 0.5))
        .unwrap();
    let sheet = laid(
        &model,
        &top,
        Point::new(5.0, 5.0, 10.0),
        Direction::X,
        SheetMethod::Geodesic,
        4.5,
    );
    let mut worst: f64 = 0.0;
    for (x, y) in [(1.0, 3.0), (5.0, 5.0), (8.0, 6.5), (6.0, 1.5)] {
        let s = sheet.sheet_point(Point::new(x, y, 10.0)).unwrap().unwrap();
        worst = worst.max((s - Point2::new(x - 5.0, y - 5.0)).magnitude());
    }
    assert!(worst <= 1e-12, "{worst}");
    assert_eq!(sheet.sheet_point(Point::new(9.9, 9.9, 10.0)).unwrap(), None);
    let (there, back) = round_trips(&sheet, 4.5);
    eprintln!("plane: {there:.3e} on the sheet, {back:.3e} in space");
    assert!(there <= ROUND_TRIP && back <= ROUND_TRIP, "{there} {back}");
}

#[test]
fn a_cylinder_unrolls_back_to_its_development() {
    // Radius 10, height 40, the sheet at (0, -10, 20) with +x round the
    // drum: the point at turn `a` from the origin and height `z` is the
    // sheet point (10 a, z - 20), its development.
    let mut model = Model::new();
    let drum = make_cylinder(&mut model, Frame::WORLD, 10.0, 40.0, T)
        .unwrap()
        .shape;
    let face = wall(&model, &drum, 0.5);
    let sheet = laid(
        &model,
        &face,
        Point::new(0.0, -10.0, 20.0),
        Direction::X,
        SheetMethod::Geodesic,
        25.0,
    );
    let mut worst: f64 = 0.0;
    for a in [-1.9, -1.0, -0.2, 0.0, 0.3, 1.2, 1.85] {
        for z in [10.5, 13.0, 20.0, 29.0] {
            let p = Point::new(10.0 * f64::sin(a), -10.0 * f64::cos(a), z);
            let s = sheet.sheet_point(p).unwrap().unwrap();
            worst = worst.max((s - Point2::new(10.0 * a, z - 20.0)).magnitude());
        }
    }
    assert!(worst <= 1e-12, "{worst}");
    // A point off the surface comes back as its foot's sheet point.
    let off = sheet
        .sheet_point(Point::new(0.0, -12.0, 21.0))
        .unwrap()
        .unwrap();
    assert!(
        (off - Point2::new(0.0, 1.0)).magnitude() <= 1e-12,
        "{off:?}"
    );
    // Half a turn away is past a reach of 25 either way round.
    assert_eq!(
        sheet.sheet_point(Point::new(0.0, 10.0, 20.0)).unwrap(),
        None
    );
    let (there, back) = round_trips(&sheet, 19.9);
    eprintln!("cylinder: {there:.3e} on the sheet, {back:.3e} in space");
    assert!(there <= ROUND_TRIP && back <= ROUND_TRIP, "{there} {back}");
}

#[test]
fn a_cone_unrolls_back_to_its_development() {
    // Base radius 10 at z = 0, 5 at z = 10: the apex at z = 20, the half
    // angle atan(1 / 2). Developed about the apex, the point at slant
    // `rho` and turn `a` from the origin's generator stands at developed
    // angle `a sin(half)` from it; the sheet's +y runs up the slant toward
    // the apex, which stands `rho0` along it.
    let mut model = Model::new();
    let cone = make_cone(&mut model, Frame::WORLD, 10.0, 5.0, 10.0, T)
        .unwrap()
        .shape;
    let face = wall(&model, &cone, 0.9);
    let sheet = laid(
        &model,
        &face,
        Point::new(0.0, -7.5, 5.0),
        Direction::X,
        SheetMethod::Geodesic,
        10.0,
    );
    let half = 0.5_f64.atan();
    let slant = |z: f64| (20.0 - z) / half.cos();
    let rho0 = slant(5.0);
    let mut worst: f64 = 0.0;
    for a in [-0.6, -0.3, 0.0, 0.25, 0.55] {
        for z in [0.5, 4.0, 5.0, 9.5] {
            let radius = (20.0 - z) / 2.0;
            let azimuth = a - PI / 2.0;
            let p = Point::new(radius * azimuth.cos(), radius * azimuth.sin(), z);
            let (rho, psi) = (slant(z), a * half.sin());
            let expected = Point2::new(rho * psi.sin(), rho.mul_add(-psi.cos(), rho0));
            let s = sheet.sheet_point(p).unwrap().unwrap();
            worst = worst.max((s - expected).magnitude());
        }
    }
    assert!(worst <= 1e-12, "{worst}");
    // Over the disc the cone's wall holds about the origin, its height
    // bounded at 0 and 10.
    let (there, back) = round_trips(&sheet, 4.9);
    eprintln!("cone: {there:.3e} on the sheet, {back:.3e} in space");
    assert!(there <= ROUND_TRIP && back <= ROUND_TRIP, "{there} {back}");
}

#[test]
fn a_sphere_comes_back_along_its_great_circles() {
    // Radius 10, the sheet at (0, -10, 0) with +x along the equator and
    // +y up the meridian: the equator's point at turn `a` is (10 a, 0) and
    // the meridian's at latitude `b` is (0, 10 b).
    let mut model = Model::new();
    let ball = make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &ball)[0].clone();
    let at = Point::new(0.0, -10.0, 0.0);
    let sheet = laid(&model, &face, at, Direction::X, SheetMethod::Geodesic, 20.0);
    let mut worst: f64 = 0.0;
    for a in [-1.5, -0.4, 0.7, 1.9] {
        let p = Point::new(10.0 * f64::sin(a), -10.0 * f64::cos(a), 0.0);
        let s = sheet.sheet_point(p).unwrap().unwrap();
        worst = worst.max((s - Point2::new(10.0 * a, 0.0)).magnitude());
        let p = Point::new(0.0, -10.0 * f64::cos(a / 2.0), 10.0 * f64::sin(a / 2.0));
        let s = sheet.sheet_point(p).unwrap().unwrap();
        worst = worst.max((s - Point2::new(0.0, 5.0 * a)).magnitude());
    }
    assert!(worst <= 1e-12, "{worst}");
    let (there, back) = round_trips(&sheet, 20.0);
    eprintln!("sphere: {there:.3e} on the sheet, {back:.3e} in space");
    assert!(there <= ROUND_TRIP && back <= ROUND_TRIP, "{there} {back}");
    // The antipode, where every geodesic from the origin meets, and the
    // chart's pole are refused by name; a point past the reach is none.
    let antipode = sheet
        .sheet_point(Point::new(0.0, 10.0, 0.0))
        .unwrap_err()
        .to_string();
    assert!(antipode.contains("opposite"), "{antipode}");
    let pole = sheet
        .sheet_point(Point::new(0.0, 0.0, 10.0))
        .unwrap_err()
        .to_string();
    assert!(pole.contains("pole"), "{pole}");
    let near = sheet_point_reach(&model, &face, at, 5.0, Point::new(10.0, 0.0, 0.0));
    assert_eq!(near, None);
}

/// `p`'s sheet point on a geodesic sheet of reach `reach` laid at `at`.
fn sheet_point_reach(
    model: &Model,
    face: &Shape,
    at: Point,
    reach: f64,
    p: Point,
) -> Option<Point2> {
    laid(model, face, at, Direction::X, SheetMethod::Geodesic, reach)
        .sheet_point(p)
        .unwrap()
}

#[test]
fn a_torus_comes_back_along_its_equator_and_meridians() {
    // Major radius 20, minor 5, the sheet at (0, -25, 0) on the outer
    // equator with +x along it: the equator and the meridian through the
    // origin are geodesics, so the equator's point at turn `a` is
    // (25 a, 0) and the meridian's at angle `b` up the tube is (0, 5 b).
    // The torus is reached by shooting geodesics.
    let mut model = Model::new();
    let torus = make_torus(&mut model, Frame::WORLD, 20.0, 5.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &torus)[0].clone();
    let sheet = laid(
        &model,
        &face,
        Point::new(0.0, -25.0, 0.0),
        Direction::X,
        SheetMethod::Geodesic,
        7.0,
    );
    let mut worst: f64 = 0.0;
    for a in [-0.25, -0.1, 0.05, 0.2] {
        let p = Point::new(25.0 * f64::sin(a), -25.0 * f64::cos(a), 0.0);
        let s = sheet.sheet_point(p).unwrap().unwrap();
        worst = worst.max((s - Point2::new(25.0 * a, 0.0)).magnitude());
    }
    for b in [-1.3, -0.5, 0.4, 1.35] {
        let p = Point::new(0.0, -(20.0 + 5.0 * f64::cos(b)), 5.0 * f64::sin(b));
        let s = sheet.sheet_point(p).unwrap().unwrap();
        worst = worst.max((s - Point2::new(0.0, 5.0 * b)).magnitude());
    }
    eprintln!("torus against its geodesics: {worst:.3e}");
    assert!(worst <= 1e-9, "{worst}");
    let (there, back) = round_trips(&sheet, 7.0);
    eprintln!("torus: {there:.3e} on the sheet, {back:.3e} in space");
    assert!(there <= ROUND_TRIP && back <= ROUND_TRIP, "{there} {back}");
}

#[test]
fn a_fitted_free_form_patch_comes_back_by_shooting() {
    // A saddle with a twist, fitted through a grid as a bicubic spline,
    // laid on at its middle; and a patch fitted to a cylinder of radius
    // 10, which must come back as the cylinder's development to the
    // fit's own reach.
    let saddle = |x: f64, y: f64| Point::new(x, y, (x * x - y * y) / 40.0 + x * y / 60.0);
    let drum = |x: f64, y: f64| {
        let a = x / 10.0;
        Point::new(10.0 * a.sin(), -10.0 * a.cos(), y)
    };
    let fitted = |shape: &dyn Fn(f64, f64) -> Point| {
        let n = 21;
        let rows: Vec<Vec<Point>> = (0..n)
            .map(|j| {
                let y = -10.0 + 20.0 * f64::from(j) / f64::from(n - 1);
                (0..n)
                    .map(|i| shape(-10.0 + 20.0 * f64::from(i) / f64::from(n - 1), y))
                    .collect::<Vec<_>>()
            })
            .collect();
        ogeom::geom::fit::fit_surface_grid(&rows, 3, 1e-8, T)
            .unwrap()
            .curve
    };

    let mut model = Model::new();
    let face = make_natural_face(&mut model, fitted(&saddle).into())
        .unwrap()
        .shape;
    let sheet = laid(
        &model,
        &face,
        saddle(0.0, 0.0),
        Direction::X,
        SheetMethod::Geodesic,
        6.0,
    );
    let (there, back) = round_trips(&sheet, 6.0);
    eprintln!("saddle: {there:.3e} on the sheet, {back:.3e} in space");
    assert!(there <= ROUND_TRIP && back <= ROUND_TRIP, "{there} {back}");

    let mut model = Model::new();
    let face = make_natural_face(&mut model, fitted(&drum).into())
        .unwrap()
        .shape;
    let at = drum(0.0, 0.0);
    let sheet = laid(&model, &face, at, Direction::X, SheetMethod::Geodesic, 6.0);
    let mut worst: f64 = 0.0;
    for x in [-5.0, -2.0, 0.0, 1.5, 4.0] {
        for y in [-3.0, 0.0, 2.5] {
            let s = sheet.sheet_point(drum(x, y)).unwrap().unwrap();
            // The patch faces out of the drum, so the sheet's +y, the
            // outward normal turned against +x, runs up it.
            worst = worst.max((s - Point2::new(x, y)).magnitude());
        }
    }
    eprintln!("fitted drum against its development: {worst:.3e}");
    assert!(worst <= 1e-5, "{worst}");
}

#[test]
fn a_flattened_sheet_comes_back_through_its_flattening() {
    // The least-distortion and the as-rigid-as-possible sheets on a
    // sphere, and the least-distortion sheet on a torus, each taken back
    // through its geodesic sheet and then its flattening's polynomial.
    let mut model = Model::new();
    let ball = make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &ball)[0].clone();
    let torus = make_torus(&mut model, Frame::WORLD, 20.0, 5.0, T)
        .unwrap()
        .shape;
    let ring = faces(&model, &torus)[0].clone();
    for (face, at, method) in [
        (
            &face,
            Point::new(0.0, -10.0, 0.0),
            SheetMethod::LeastDistortion { resolution: 1e-3 },
        ),
        (
            &face,
            Point::new(0.0, -10.0, 0.0),
            SheetMethod::AsRigidAsPossible {
                resolution: 1e-3,
                stretch: 0.05,
            },
        ),
        (
            &ring,
            Point::new(0.0, -25.0, 0.0),
            SheetMethod::LeastDistortion { resolution: 1e-3 },
        ),
    ] {
        let sheet = laid(&model, face, at, Direction::X, method, 5.0);
        let (there, back) = round_trips(&sheet, 5.0);
        eprintln!("{method:?} at {at:?}: {there:.3e} on the sheet, {back:.3e} in space");
        assert!(there <= ROUND_TRIP && back <= ROUND_TRIP, "{there} {back}");
    }
}

#[test]
fn an_edge_crossing_onto_a_sibling_face_carries_a_pcurve_there() {
    // A drum of radius 10 cut through by a slab across y = 0: its wall
    // comes back as two faces of the one cylinder, one each side of the
    // slab. A sheet line laid on the +y half runs round the drum across the
    // slab onto the -y half, climbing as it goes from the bottom rim to the
    // top one.
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
        .filter(|f| ogeom::algo::face_normal(&model, f, T).is_ok_and(|(_, n)| n.z.abs() < 0.5))
        .filter(|f| {
            ogeom::algo::face_normal(&model, f, T)
                .is_ok_and(|(p, _)| (p.x.hypot(p.y) - 10.0).abs() < 1e-6)
        })
        .collect();
    assert_eq!(walls.len(), 2, "the wall in two halves");
    let near = |f: &Shape| ogeom::algo::face_normal(&model, f, T).unwrap().0.y;
    let (upper, lower) = if near(&walls[0]) > 0.0 {
        (walls[0].clone(), walls[1].clone())
    } else {
        (walls[1].clone(), walls[0].clone())
    };
    let sheet = laid(
        &model,
        &upper,
        Point::new(0.0, 10.0, 5.0),
        Direction::new(Vector::new(-1.0, 0.0, 0.0), T).unwrap(),
        SheetMethod::Geodesic,
        30.0,
    );
    // +x round the drum toward -x, +y up it: from the bottom rim at
    // s.x = 10 to the top rim at s.x = 20, across the slab at a quarter
    // turn from the origin.
    let line: PlanarCurve = Line2d::segment(Point2::new(10.0, -5.0), Point2::new(20.0, 5.0), T)
        .unwrap()
        .into();
    let range = line.domain();
    let edge = sheet
        .chart_edge(&mut model, &line, range, 1e-6)
        .unwrap()
        .shape;

    // The edge's pcurve on each half's surface lifts onto its curve.
    let surface_of =
        |model: &Model, f: &Shape| model.node(f).unwrap().data().as_face().unwrap().surface;
    let curve = {
        let data = model.node(&edge).unwrap().data().as_edge().unwrap();
        match data.curve3d().unwrap() {
            ogeom::topo::EdgeRepr::Curve3d { curve, .. } => {
                model.geometry().curve(*curve).unwrap().clone()
            }
            _ => unreachable!(),
        }
    };
    let mut worst: f64 = 0.0;
    for half in [&upper, &lower] {
        let id = surface_of(&model, half);
        let data = model.node(&edge).unwrap().data().as_edge().unwrap();
        let Some(ogeom::topo::EdgeRepr::PCurve { curve: pc, .. }) =
            data.pcurve_for(id, &ogeom::topo::Location::identity())
        else {
            panic!("no pcurve on {half:?}'s surface");
        };
        let pcurve = model.geometry().pcurve(*pc).unwrap();
        let surface = model.geometry().surface(id).unwrap();
        for k in 0..=200 {
            let t = (range.1 - range.0).mul_add(f64::from(k) / 200.0, range.0);
            let uv = pcurve.point_at(t, T).unwrap();
            let on = surface.point_at(uv.x, uv.y, T).unwrap();
            worst = worst.max(on.distance(curve.point_at(t, T).unwrap()));
        }
    }
    eprintln!("pcurves on both halves against the edge: {worst:.3e}");
    assert!(worst <= T.confusion(), "{worst}");

    // Each half splits along it, and the solid checks clean.
    let mut shape = halves;
    for half in [upper, lower] {
        let at = ogeom::algo::face_normal(&model, &half, T).unwrap().0;
        let face = faces(&model, &shape)
            .into_iter()
            .find(|f| {
                ogeom::algo::face_normal(&model, f, T).is_ok_and(|(p, _)| p.distance(at) < 1e-9)
            })
            .unwrap();
        let built = split_face(
            &mut model,
            &shape,
            &face,
            core::slice::from_ref(&edge),
            Projection::OnFace,
            T,
        )
        .unwrap();
        assert_eq!(built.history.trace(&face).len(), 2, "the half in two");
        shape = built.shape;
        let report = check(&model, &shape, T).unwrap();
        assert!(report.problems.is_empty(), "{:?}", report.problems);
    }
}

#[test]
fn an_edge_carries_its_pcurve_onto_every_face_of_an_equal_surface() {
    // Four drums, each its own surface: one the sheet is laid on, one
    // standing in the same place with the same chart, one in the same
    // place with its chart turned a quarter turn about the axis, and one
    // of radius 12. The carried edge has its pcurve on the first two only:
    // the turned chart lands the same (u, v) elsewhere, and the wider drum
    // is another surface. The second drum's wall then splits along it.
    let mut model = Model::new();
    let turned = Frame::new(Point::ORIGIN, Direction::Z, Direction::Y, T).unwrap();
    let drums: Vec<Shape> = [
        (Frame::WORLD, 10.0),
        (Frame::WORLD, 10.0),
        (turned, 10.0),
        (Frame::WORLD, 12.0),
    ]
    .into_iter()
    .map(|(frame, radius)| {
        make_cylinder(&mut model, frame, radius, 10.0, T)
            .unwrap()
            .shape
    })
    .collect();
    let walls: Vec<Shape> = drums.iter().map(|d| wall(&model, d, 0.5)).collect();
    let surface_of =
        |model: &Model, f: &Shape| model.node(f).unwrap().data().as_face().unwrap().surface;
    let sheet = laid(
        &model,
        &walls[0],
        Point::new(0.0, -10.0, 5.0),
        Direction::X,
        SheetMethod::Geodesic,
        8.0,
    );
    // From the bottom rim to the top one.
    let line: PlanarCurve = Line2d::segment(Point2::new(-4.0, -5.0), Point2::new(4.0, 5.0), T)
        .unwrap()
        .into();
    let range = line.domain();
    let edge = sheet
        .chart_edge(&mut model, &line, range, 1e-6)
        .unwrap()
        .shape;
    let data = model.node(&edge).unwrap().data().as_edge().unwrap().clone();
    let has: Vec<bool> = walls
        .iter()
        .map(|w| {
            data.pcurve_for(surface_of(&model, w), &ogeom::topo::Location::identity())
                .is_some()
        })
        .collect();
    assert_eq!(has, [true, true, false, false]);
    assert!(data.same_parameter());

    let built = split_face(
        &mut model,
        &drums[1],
        &walls[1],
        &[edge],
        Projection::OnFace,
        T,
    )
    .unwrap();
    assert_eq!(built.history.trace(&walls[1]).len(), 2, "the wall in two");
    let report = check(&model, &built.shape, T).unwrap();
    assert!(report.problems.is_empty(), "{:?}", report.problems);
}
