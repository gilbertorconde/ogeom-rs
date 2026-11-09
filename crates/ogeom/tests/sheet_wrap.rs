//! Flat sketches laid onto curved faces by geodesic polar coordinates,
//! measured against the exponential map's closed forms: great circles and
//! geodesic circles on a sphere, the outer equator and the meridians of a
//! torus, and the developments of a cylinder and a cone.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{PI, TAU};

use ogeom::algo::{
    SheetMethod, SheetPlacement, SurfaceSheet, check, make_cone, make_cylinder, make_sphere,
    make_torus, surface_properties, volume_properties,
};
use ogeom::boolean::fuse;
use ogeom::core::Tolerances;
use ogeom::geom::{Circle2d, Curve2d as _, Line2d, PlanarCurve, Surface as _};
use ogeom::heal::{Projection, split_face};
use ogeom::math::{Circle2, Direction, Frame, Point, Point2, Vector};
use ogeom::mesh::Deflection;
use ogeom::offset::make_thick_sheet;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The fit tolerance asked of every carried curve.
const FIT: f64 = 1e-5;

fn fine() -> Deflection {
    Deflection {
        chord: 1e-4,
        angular: 0.01,
        ..Deflection::default()
    }
}

fn faces(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Face).unwrap()
}

fn area(model: &Model, face: &Shape) -> f64 {
    surface_properties(model, face, fine(), T).unwrap().mass
}

/// The 6 x 4 rectangle centred on the sheet's origin.
const RECTANGLE: [(f64, f64); 4] = [(-3.0, -2.0), (3.0, -2.0), (3.0, 2.0), (-3.0, 2.0)];

fn sheet_on(model: &Model, face: &Shape, at: Point, reach: f64) -> SurfaceSheet {
    SurfaceSheet::new(
        model,
        face,
        &SheetPlacement {
            at,
            x_direction: Direction::X,
            method: SheetMethod::Geodesic,
        },
        reach,
        T,
    )
    .unwrap()
}

/// The rectangle's sides carried onto the face as edges.
fn rectangle_edges(model: &mut Model, sheet: &SurfaceSheet) -> Vec<Shape> {
    (0..4)
        .map(|k| {
            let (a, b) = (RECTANGLE[k], RECTANGLE[(k + 1) % 4]);
            let line: PlanarCurve =
                Line2d::segment(Point2::new(a.0, a.1), Point2::new(b.0, b.1), T)
                    .unwrap()
                    .into();
            let range = line.domain();
            sheet.chart_edge(model, &line, range, FIT).unwrap().shape
        })
        .collect()
}

/// The face of `shape` split along the carried rectangle: the pieces
/// inside it (more than one where a seam of the face runs through it),
/// their area, and the centroid of their vertices.
struct Inside {
    shape: Shape,
    pieces: Vec<Shape>,
    area: f64,
    centroid: Point,
}

fn rectangle_inside(
    model: &mut Model,
    shape: &Shape,
    face: &Shape,
    sheet: &SurfaceSheet,
) -> Inside {
    let edges = rectangle_edges(model, sheet);
    let whole = area(model, face);
    let built = split_face(model, shape, face, &edges, Projection::OnFace, T).unwrap();
    assert!(check(model, &built.shape, T).unwrap().is_usable());
    let mut by_area: Vec<(f64, Shape)> = built
        .history
        .trace(face)
        .iter()
        .map(|p| (area(model, p), p.clone()))
        .collect();
    assert!(by_area.len() >= 2, "the rectangle cuts the face");
    by_area.sort_by(|a, b| a.0.total_cmp(&b.0));
    let sum: f64 = by_area.iter().map(|p| p.0).sum();
    assert!((sum - whole).abs() <= whole * 1e-6, "{sum} of {whole}");
    // The outside is the rest of the face, larger than the rectangle.
    by_area.pop();
    let compound = model
        .add_compound(&by_area.iter().map(|p| p.1.clone()).collect::<Vec<_>>())
        .unwrap();
    let vertices = explore_unique(model, &compound, ShapeType::Vertex).unwrap();
    let mut total = Vector::new(0.0, 0.0, 0.0);
    for v in &vertices {
        total += model
            .node(v)
            .unwrap()
            .data()
            .as_vertex()
            .unwrap()
            .point
            .to_vector();
    }
    let count = f64::from(u32::try_from(vertices.len()).unwrap());
    Inside {
        shape: built.shape.clone(),
        area: by_area.iter().map(|p| p.0).sum(),
        pieces: by_area.into_iter().map(|p| p.1).collect(),
        centroid: Point::from_vector(total * (1.0 / count)),
    }
}

/// The angle about `z` between where two sheet points land.
fn turn_about_z(sheet: &SurfaceSheet, a: Point2, b: Point2) -> f64 {
    let lift = |p: Point2| {
        let uv = sheet.chart_point(p).unwrap();
        sheet.surface().point_at(uv.x, uv.y, T).unwrap()
    };
    let (pa, pb) = (lift(a), lift(b));
    let turn = pb.y.atan2(pb.x) - pa.y.atan2(pa.x);
    (turn + PI).rem_euclid(TAU) - PI
}

#[test]
fn a_rectangle_laid_on_a_sphere_keeps_its_area_its_middle_row_and_its_place() {
    let mut model = Model::new();
    let sphere = make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &sphere)[0].clone();
    let sheet = sheet_on(&model, &face, Point::new(0.0, -10.0, 0.0), 5.0);
    let inside = rectangle_inside(&mut model, &sphere, &face, &sheet);
    let a = inside.area;
    assert!((a - 24.0).abs() <= 24.0 * 0.03, "{a}");
    let centroid = inside.centroid;
    assert!(
        centroid.x.abs() <= 0.05 && centroid.z.abs() <= 0.05,
        "{centroid:?}"
    );
    // The equator is a great circle: the middle row runs along it, 6 / R.
    let turn = turn_about_z(&sheet, Point2::new(-3.0, 0.0), Point2::new(3.0, 0.0));
    assert!((turn - 0.6).abs() <= 0.6 * 1e-9, "{turn}");
}

#[test]
fn a_rectangle_laid_on_a_torus_keeps_its_area_its_middle_row_and_its_place() {
    let mut model = Model::new();
    let torus = make_torus(&mut model, Frame::WORLD, 20.0, 5.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &torus)[0].clone();
    let sheet = sheet_on(&model, &face, Point::new(0.0, -25.0, 0.0), 5.0);
    let inside = rectangle_inside(&mut model, &torus, &face, &sheet);
    let a = inside.area;
    assert!((a - 24.0).abs() <= 24.0 * 0.03, "{a}");
    let centroid = inside.centroid;
    assert!(
        centroid.x.abs() <= 0.05 && centroid.z.abs() <= 0.05,
        "{centroid:?}"
    );
    // The outer equator is a geodesic of radius 25.
    let turn = turn_about_z(&sheet, Point2::new(-3.0, 0.0), Point2::new(3.0, 0.0));
    assert!((turn - 6.0 / 25.0).abs() <= 6.0 / 25.0 * 1e-6, "{turn}");
}

#[test]
fn a_line_across_the_torus_equator_lands_on_its_meridian() {
    // At the outer equator the sheet's +y runs up the tube: the meridian
    // through the origin is a geodesic, and length r along it is the tube
    // angle r / 5.
    let mut model = Model::new();
    let torus = make_torus(&mut model, Frame::WORLD, 20.0, 5.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &torus)[0].clone();
    let sheet = sheet_on(&model, &face, Point::new(0.0, -25.0, 0.0), 8.0);
    for k in 1..=8 {
        let r = f64::from(k);
        let uv = sheet.chart_point(Point2::new(0.0, r)).unwrap();
        let p = sheet.surface().point_at(uv.x, uv.y, T).unwrap();
        let tube = r / 5.0;
        let expected = Point::new(0.0, -(20.0 + 5.0 * tube.cos()), 5.0 * tube.sin());
        assert!(p.distance(expected) <= 1e-6, "{r}: {p:?} {expected:?}");
    }
}

#[test]
fn a_ray_on_a_torus_keeps_clairauts_constant() {
    // Every geodesic of a surface of revolution keeps `rho sin(psi)`, its
    // distance from the axis times the sine of its angle to the meridian:
    // at the origin on the outer equator, leaving at 40 degrees from the
    // parallel, it is 25 cos(40 degrees).
    let mut model = Model::new();
    let torus = make_torus(&mut model, Frame::WORLD, 20.0, 5.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &torus)[0].clone();
    let sheet = sheet_on(&model, &face, Point::new(0.0, -25.0, 0.0), 12.0);
    let (s, c) = 40.0_f64.to_radians().sin_cos();
    let lift = |r: f64| {
        let uv = sheet.chart_point(Point2::new(r * c, r * s)).unwrap();
        sheet.surface().point_at(uv.x, uv.y, T).unwrap()
    };
    let constant = 25.0 * c;
    for k in 1..=10 {
        let r = f64::from(k);
        let h = 1e-3;
        let (p, ahead, behind) = (lift(r), lift(r + h), lift(r - h));
        let tangent = (ahead - behind) * (1.0 / (2.0 * h));
        let rho = p.x.hypot(p.y);
        let around = Vector::new(-p.y / rho, p.x / rho, 0.0);
        let held = rho * tangent.dot(around).abs() / tangent.magnitude();
        assert!(
            (held - constant).abs() <= 1e-6,
            "{r}: {held} against {constant}"
        );
        // And the geodesic is run at unit speed: lengths along the ray kept.
        assert!((tangent.magnitude() - 1.0).abs() <= 1e-6, "{r}");
    }
}

#[test]
fn a_circle_round_the_origin_lands_on_the_spheres_geodesic_circle() {
    // The sheet circle of radius 4 about the origin lands where every point
    // stands 4 / R radians from the origin's point, a small circle; its
    // stretch across the rays is sin(0.4) / 0.4.
    let mut model = Model::new();
    let sphere = make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &sphere)[0].clone();
    let sheet = sheet_on(&model, &face, Point::new(0.0, -10.0, 0.0), 4.0);
    let circle: PlanarCurve =
        Circle2d::new(Circle2::centred(Point2::ORIGIN, 4.0, T).unwrap()).into();
    let pcurve = sheet.chart_curve(&circle, (0.0, TAU), FIT).unwrap();
    let expected = 10.0 * 0.4_f64.cos();
    for k in 0..=200 {
        let t = TAU * f64::from(k) / 200.0;
        let uv = pcurve.point_at(t, T).unwrap();
        let p = sheet.surface().point_at(uv.x, uv.y, T).unwrap();
        // Distance along the origin's radius, and the sheet point it carries.
        assert!((-p.y - expected).abs() <= FIT, "{t}: {p:?}");
        let (s, c) = t.sin_cos();
        let sheet_point = Point::new(
            10.0 * 0.4_f64.sin() * c,
            -expected,
            10.0 * 0.4_f64.sin() * s,
        );
        assert!(p.distance(sheet_point) <= FIT, "{t}: {p:?} {sheet_point:?}");
    }
    let distortion = sheet.distortion();
    let across = 0.4_f64.sin() / 0.4;
    assert!(
        (distortion.area.0 - across).abs() <= 1e-12,
        "{distortion:?}"
    );
    assert!(distortion.area.1 <= 1.0 && distortion.length == (distortion.area.0, 1.0));
}

#[test]
fn a_rectangle_on_a_cylinder_or_a_cone_is_their_exact_development() {
    let mut model = Model::new();
    let drum = make_cylinder(&mut model, Frame::WORLD, 10.0, 10.0, T)
        .unwrap()
        .shape;
    let wall = faces(&model, &drum)
        .into_iter()
        .find(|f| ogeom::algo::face_normal(&model, f, T).is_ok_and(|(_, n)| n.z.abs() < 0.5))
        .unwrap();
    let sheet = sheet_on(&model, &wall, Point::new(0.0, -10.0, 5.0), 4.0);
    let inside = rectangle_inside(&mut model, &drum, &wall, &sheet);
    let a = inside.area;
    assert!((a - 24.0).abs() <= 24.0 * 1e-6, "{a}");
    assert_eq!(sheet.distortion().area, (1.0, 1.0));

    // A cone of half angle atan(1 / 2): base radius 10, top 5, height 10.
    let mut model = Model::new();
    let cone = make_cone(&mut model, Frame::WORLD, 10.0, 5.0, 10.0, T)
        .unwrap()
        .shape;
    let wall = faces(&model, &cone)
        .into_iter()
        .find(|f| ogeom::algo::face_normal(&model, f, T).is_ok_and(|(_, n)| n.z.abs() < 0.9))
        .unwrap();
    let sheet = sheet_on(&model, &wall, Point::new(0.0, -7.5, 5.0), 4.0);
    let inside = rectangle_inside(&mut model, &cone, &wall, &sheet);
    let a = inside.area;
    assert!((a - 24.0).abs() <= 24.0 * 1e-6, "{a}");
    // Lengths in the development: the sheet's x axis is a straight line
    // in it, 6 long between the rectangle's middle sides.
    let lift = |p: Point2| {
        let uv = sheet.chart_point(p).unwrap();
        sheet.surface().point_at(uv.x, uv.y, T).unwrap()
    };
    let up = lift(Point2::new(0.0, 2.0)).distance(lift(Point2::new(0.0, -2.0)));
    assert!((up - 4.0).abs() <= 1e-9, "{up}");
}

#[test]
fn a_sheet_is_refused_where_its_map_folds_or_its_origin_is_a_pole() {
    let mut model = Model::new();
    let sphere = make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &sphere)[0].clone();
    let placement = |at| SheetPlacement {
        at,
        x_direction: Direction::X,
        method: SheetMethod::Geodesic,
    };
    let antipode = SurfaceSheet::new(
        &model,
        &face,
        &placement(Point::new(0.0, -10.0, 0.0)),
        32.0,
        T,
    );
    assert!(
        antipode
            .unwrap_err()
            .to_string()
            .contains("opposite its origin")
    );
    let pole = SurfaceSheet::new(
        &model,
        &face,
        &placement(Point::new(0.0, 0.0, 10.0)),
        2.0,
        T,
    );
    assert!(pole.unwrap_err().to_string().contains("pole"));
    let sheet = sheet_on(&model, &face, Point::new(0.0, -10.0, 0.0), 2.0);
    let far = sheet.chart_point(Point2::new(3.0, 0.0));
    assert!(far.unwrap_err().to_string().contains("reach"));
    // The sheet point a quarter turn up lands on the north pole.
    let sheet = sheet_on(&model, &face, Point::new(0.0, -10.0, 0.0), 20.0);
    let up = sheet.chart_point(Point2::new(0.0, 5.0 * PI));
    assert!(up.unwrap_err().to_string().contains("pole"));
}

#[test]
fn a_rectangle_laid_on_a_sphere_thickens_outward_and_fuses_onto_it() {
    let mut model = Model::new();
    let sphere = make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &sphere)[0].clone();
    let sheet = sheet_on(&model, &face, Point::new(0.0, -10.0, 0.0), 5.0);
    let inside = rectangle_inside(&mut model, &sphere, &face, &sheet);
    let piece = model.add_shell(&inside.pieces).unwrap();
    let raised = make_thick_sheet(&mut model, &piece, 1.0, false, T)
        .unwrap()
        .shape;
    let raised_volume = volume_properties(&model, &raised, fine(), T).unwrap().mass;
    // The shell between radii 10 and 11 over the patch: area times
    // (1 + 1/R + 1/(3 R^2)) for a sphere's patch thickened outward.
    let expected = inside.area * (1.0 + 0.1 + 1.0 / 300.0);
    assert!(
        (raised_volume - expected).abs() <= expected * 1e-3,
        "{raised_volume} {expected}"
    );
    let fused = fuse(&mut model, &inside.shape, &raised, T).unwrap().shape;
    assert!(check(&model, &fused, T).unwrap().is_usable());
    let volume = volume_properties(&model, &fused, fine(), T).unwrap().mass;
    let ball = 4.0 / 3.0 * PI * 1000.0;
    assert!(
        (volume - ball - expected).abs() <= expected * 1e-3,
        "{volume}"
    );
}

#[test]
fn a_circle_round_the_spheres_pole_cuts_off_its_geodesic_cap() {
    // Laid ten degrees from the north pole, the sheet circle of radius 3
    // about the origin goes round the pole: its chart curve runs a whole
    // turn of longitude, and it cuts off the cap of geodesic radius 3,
    // 2 pi R^2 (1 - cos(3 / R)).
    let mut model = Model::new();
    let sphere = make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &sphere)[0].clone();
    let latitude = 80.0_f64.to_radians();
    let at = Point::new(0.0, -10.0 * latitude.cos(), 10.0 * latitude.sin());
    let sheet = sheet_on(&model, &face, at, 3.0);
    let circle: PlanarCurve =
        Circle2d::new(Circle2::centred(Point2::ORIGIN, 3.0, T).unwrap()).into();
    let pcurve = sheet.chart_curve(&circle, (0.0, TAU), FIT).unwrap();
    let (start, end) = (
        pcurve.point_at(0.0, T).unwrap(),
        pcurve.point_at(TAU, T).unwrap(),
    );
    assert!(
        ((end.x - start.x).abs() - TAU).abs() <= 1e-9,
        "{start:?} {end:?}"
    );
    assert!((end.y - start.y).abs() <= 1e-9, "{start:?} {end:?}");
    for k in 0..=400 {
        let t = TAU * f64::from(k) / 400.0;
        let uv = pcurve.point_at(t, T).unwrap();
        let p = sheet.surface().point_at(uv.x, uv.y, T).unwrap();
        // Every point stands 3 / R radians from the origin's point.
        let apart = (p.to_vector().dot(at.to_vector()) / 100.0)
            .clamp(-1.0, 1.0)
            .acos();
        assert!((apart * 10.0 - 3.0).abs() <= FIT, "{t}: {apart}");
    }
    let edge = sheet
        .chart_edge(&mut model, &circle, (0.0, TAU), FIT)
        .unwrap()
        .shape;
    let whole = area(&model, &face);
    let built = split_face(&mut model, &sphere, &face, &[edge], Projection::OnFace, T).unwrap();
    assert!(check(&model, &built.shape, T).unwrap().is_usable());
    let mut areas: Vec<f64> = built
        .history
        .trace(&face)
        .iter()
        .map(|p| area(&model, p))
        .collect();
    areas.sort_by(f64::total_cmp);
    assert_eq!(areas.len(), 2, "{areas:?}");
    let cap = TAU * 100.0 * (1.0 - 0.3_f64.cos());
    assert!(
        (areas[0] - cap).abs() <= cap * 1e-5,
        "{areas:?} against {cap}"
    );
    assert!(
        (areas[0] + areas[1] - whole).abs() <= whole * 1e-6,
        "{areas:?}"
    );
}

/// The resolution asked of every least-distortion sheet.
const RESOLUTION: f64 = 1e-3;

fn flat_sheet_on(model: &Model, face: &Shape, at: Point, reach: f64) -> SurfaceSheet {
    SurfaceSheet::new(
        model,
        face,
        &SheetPlacement {
            at,
            x_direction: Direction::X,
            method: SheetMethod::LeastDistortion {
                resolution: RESOLUTION,
            },
        },
        reach,
        T,
    )
    .unwrap()
}

#[test]
fn a_least_distortion_sheet_on_a_sphere_is_its_stereographic_projection() {
    // A geodesic disc of a sphere flattened conformally about its centre
    // is, by its symmetry, the stereographic projection from the antipode,
    // scaled: the sheet point at radius r lands at the geodesic angle
    // theta from the origin where r = c 2R tan(theta / 2), c one constant,
    // on the ray at the sheet point's own angle.
    let mut model = Model::new();
    let sphere = make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &sphere)[0].clone();
    let at = Point::new(0.0, -10.0, 0.0);
    let sheet = flat_sheet_on(&model, &face, at, 5.0);
    assert!(
        sheet.flattening_error() <= RESOLUTION,
        "{}",
        sheet.flattening_error()
    );
    let mut scales: Vec<f64> = Vec::new();
    for k in 0..12 {
        let angle = f64::from(k) * TAU / 12.0;
        for r in [1.0, 2.5, 4.0, 5.0] {
            let (s, c) = angle.sin_cos();
            let uv = sheet.chart_point(Point2::new(r * c, r * s)).unwrap();
            let p = sheet.surface().point_at(uv.x, uv.y, T).unwrap();
            let theta = (-p.y / 10.0).clamp(-1.0, 1.0).acos();
            scales.push(r / (20.0 * (theta / 2.0).tan()));
            // On the ray at the sheet point's angle: x along the sheet's +x,
            // z along its +y (the sphere's outward normal at the origin is
            // -Y, and +X turned about it is +Z).
            let turned = p.z.atan2(p.x);
            let off = (turned - angle + PI).rem_euclid(TAU) - PI;
            assert!(
                off.abs() * 10.0 * theta.sin() <= RESOLUTION,
                "{r} {angle}: {off}"
            );
        }
    }
    let (least, most) = scales
        .iter()
        .fold((f64::INFINITY, 0.0_f64), |(a, b), c| (a.min(*c), b.max(*c)));
    // A miss of the resolution at radius 1 moves the constant by at most
    // that fraction.
    assert!(most - least <= RESOLUTION * most, "{least} {most}");
    let distortion = sheet.distortion();
    // Conformal: at each point every direction stretches alike, so the
    // area ratio's range is the length ratio's squared.
    let (l, a) = (distortion.length, distortion.area);
    assert!((a.0 - l.0 * l.0).abs() <= 1e-2 * a.0, "{distortion:?}");
    assert!((a.1 - l.1 * l.1).abs() <= 1e-2 * a.1, "{distortion:?}");
}

#[test]
fn a_rectangle_laid_least_distorted_on_a_sphere_cuts_out_the_stereographic_area() {
    // The sheet is the stereographic projection scaled, r = k tan(theta / 2),
    // so a sheet region's area on the sphere is the integral over it of the
    // squared stretch (2R / k)^2 / (1 + r^2 / k^2)^2; k is read off the
    // sheet's own rim of reach.
    let mut model = Model::new();
    let sphere = make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &sphere)[0].clone();
    let sheet = flat_sheet_on(&model, &face, Point::new(0.0, -10.0, 0.0), 5.0);
    let uv = sheet.chart_point(Point2::new(5.0, 0.0)).unwrap();
    let p = sheet.surface().point_at(uv.x, uv.y, T).unwrap();
    let theta = (-p.y / 10.0).clamp(-1.0, 1.0).acos();
    let k = 5.0 / (theta / 2.0).tan();
    let inside = rectangle_inside(&mut model, &sphere, &face, &sheet);
    // Gauss-Legendre on five points a side, over a grid of 6 x 4 unit cells.
    let nodes = [
        (0.0, 0.568_888_888_888_888_9),
        (-0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
        (0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
        (-0.906_179_845_938_664, 0.236_926_885_056_189_1),
        (0.906_179_845_938_664, 0.236_926_885_056_189_1),
    ];
    let mut expected = 0.0;
    for i in 0..6 {
        for j in 0..4 {
            for (sx, wx) in nodes {
                for (sy, wy) in nodes {
                    let x = f64::from(i) - 3.0 + 0.5 + sx / 2.0;
                    let y = f64::from(j) - 2.0 + 0.5 + sy / 2.0;
                    let ratio = (x * x + y * y) / (k * k);
                    expected += wx * wy / 4.0 * (20.0 / k).powi(2) / (1.0 + ratio).powi(2);
                }
            }
        }
    }
    // Each side may stand off the smooth map's by the stated error times
    // the most stretch, and off its chart curve by the fit, over the
    // rectangle's perimeter of 20.
    let stretch = sheet.distortion().length.1;
    let allowed = 20.0 * (sheet.flattening_error() * stretch + FIT);
    assert!(
        (inside.area - expected).abs() <= allowed,
        "{} against {expected}, allowed {allowed}",
        inside.area
    );
    assert!(inside.centroid.x.abs() <= 0.05 && inside.centroid.z.abs() <= 0.05);
}

#[test]
fn a_rectangle_laid_least_distorted_on_a_torus_keeps_its_angles_and_place() {
    let mut model = Model::new();
    let torus = make_torus(&mut model, Frame::WORLD, 20.0, 5.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &torus)[0].clone();
    let sheet = flat_sheet_on(&model, &face, Point::new(0.0, -25.0, 0.0), 5.0);
    assert!(
        sheet.flattening_error() <= RESOLUTION,
        "{}",
        sheet.flattening_error()
    );
    let inside = rectangle_inside(&mut model, &torus, &face, &sheet);
    assert!(inside.centroid.x.abs() <= 0.05 && inside.centroid.z.abs() <= 0.05);
    // Conformal: the area ratio's range is the length ratio's squared.
    let distortion = sheet.distortion();
    let (l, a) = (distortion.length, distortion.area);
    assert!((a.0 - l.0 * l.0).abs() <= 1e-2 * a.0, "{distortion:?}");
    assert!((a.1 - l.1 * l.1).abs() <= 1e-2 * a.1, "{distortion:?}");
    // The torus curves less than the sphere: the stretch stays within the
    // sphere's own over the same reach.
    assert!(
        l.0 > 0.9 && l.1 < 1.1,
        "{distortion:?}, area {}",
        inside.area
    );
}

#[test]
fn a_least_distortion_sheet_on_a_cylinder_is_its_development() {
    let mut model = Model::new();
    let drum = make_cylinder(&mut model, Frame::WORLD, 10.0, 10.0, T)
        .unwrap()
        .shape;
    let wall = faces(&model, &drum)
        .into_iter()
        .find(|f| ogeom::algo::face_normal(&model, f, T).is_ok_and(|(_, n)| n.z.abs() < 0.5))
        .unwrap();
    let sheet = flat_sheet_on(&model, &wall, Point::new(0.0, -10.0, 5.0), 4.0);
    assert_eq!(sheet.flattening_error(), 0.0);
    let inside = rectangle_inside(&mut model, &drum, &wall, &sheet);
    assert!((inside.area - 24.0).abs() <= 24.0 * 1e-6, "{}", inside.area);
    assert_eq!(sheet.distortion().area, (1.0, 1.0));
}

#[test]
fn a_sheet_circle_through_the_spheres_pole_cuts_off_its_cap() {
    // Laid ten degrees from the north pole, the sheet circle of radius ten
    // degrees about the origin runs through the pole: carried as two
    // pieces meeting at the pole, it cuts off the cap of geodesic radius
    // ten degrees, 2 pi R^2 (1 - cos 10 degrees).
    let mut model = Model::new();
    let sphere = make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &sphere)[0].clone();
    let angle = 10.0_f64.to_radians();
    let at = Point::new(0.0, -10.0 * angle.sin(), 10.0 * angle.cos());
    let sheet = sheet_on(&model, &face, at, 2.0);
    let circle: PlanarCurve =
        Circle2d::new(Circle2::centred(Point2::ORIGIN, 10.0 * angle, T).unwrap()).into();
    assert!(sheet.chart_curve(&circle, (0.0, TAU), FIT).is_err());
    let edges = sheet
        .chart_edges(&mut model, &circle, (0.0, TAU), FIT)
        .unwrap();
    assert_eq!(edges.len(), 2);
    let whole = area(&model, &face);
    let built = split_face(&mut model, &sphere, &face, &edges, Projection::OnFace, T).unwrap();
    assert!(check(&model, &built.shape, T).unwrap().is_usable());
    let mut areas: Vec<f64> = built
        .history
        .trace(&face)
        .iter()
        .map(|p| area(&model, p))
        .collect();
    areas.sort_by(f64::total_cmp);
    assert_eq!(areas.len(), 2, "{areas:?}");
    let cap = TAU * 100.0 * (1.0 - angle.cos());
    assert!(
        (areas[0] - cap).abs() <= cap * 1e-6,
        "{areas:?} against {cap}"
    );
    assert!(
        (areas[0] + areas[1] - whole).abs() <= whole * 1e-6,
        "{areas:?}"
    );
}

#[test]
fn a_sheet_line_through_the_pole_halves_a_hemisphere() {
    // The sphere cut at its equator, and its upper half cut by a sheet
    // line laid at 45 degrees north along the meridian: the line runs from
    // the equator over the pole to the equator opposite, carried as two
    // pieces meeting at the pole, and halves the hemisphere into quarter
    // spheres of pi R^2 each.
    let mut model = Model::new();
    let sphere = make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    let face = faces(&model, &sphere)[0].clone();
    let equator = {
        let circle = ogeom::math::Circle::new(Frame::WORLD, 10.0, T).unwrap();
        ogeom::algo::make_edge(
            &mut model,
            ogeom::geom::Curve::from(ogeom::geom::CircleCurve::new(circle)),
            (0.0, TAU),
            T,
        )
        .unwrap()
        .shape
    };
    let halves = split_face(
        &mut model,
        &sphere,
        &face,
        &[equator],
        Projection::OnFace,
        T,
    )
    .unwrap();
    let upper = halves
        .history
        .trace(&face)
        .iter()
        .find(|f| ogeom::algo::face_normal(&model, f, T).is_ok_and(|(p, _)| p.z > 0.0))
        .cloned()
        .unwrap();
    let north = 45.0_f64.to_radians();
    let at = Point::new(0.0, -10.0 * north.cos(), 10.0 * north.sin());
    let sheet = sheet_on(&model, &upper, at, 24.0);
    // From just below the equator behind the origin to just below it on
    // the far side of the pole, 45 and 135 degrees of arc away.
    let line: PlanarCurve = Line2d::segment(Point2::new(0.0, -8.0), Point2::new(0.0, 23.7), T)
        .unwrap()
        .into();
    let range = line.domain();
    let edges = sheet.chart_edges(&mut model, &line, range, FIT).unwrap();
    assert_eq!(edges.len(), 2);
    let built = split_face(
        &mut model,
        &halves.shape,
        &upper,
        &edges,
        Projection::OnFace,
        T,
    )
    .unwrap();
    assert!(check(&model, &built.shape, T).unwrap().is_usable());
    // The sphere's seam meridian runs through the quarter on the +x side
    // and stays an edge of it, so that quarter comes back in two pieces.
    let mut areas: Vec<f64> = built
        .history
        .trace(&upper)
        .iter()
        .map(|p| area(&model, p))
        .collect();
    areas.sort_by(f64::total_cmp);
    let quarter = PI * 100.0;
    let largest = areas.pop().unwrap();
    let rest: f64 = areas.iter().sum();
    assert!(
        (largest - quarter).abs() <= quarter * 1e-6,
        "{largest} against {quarter}"
    );
    assert!(
        (rest - quarter).abs() <= quarter * 1e-6,
        "{areas:?} against {quarter}"
    );
}
