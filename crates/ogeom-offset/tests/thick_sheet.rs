//! Sheets offset and thickened: boxes, annular sectors and fitted
//! saddles against their closed forms, and the refusals by name.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_algo::{
    check, edge_vertices, make_edge, make_edge_between, make_face_with_pcurves, make_polygon,
    make_vertex, tight_bounds,
};
use ogeom_core::Tolerances;
use ogeom_geom::{
    BSplineSurface, CircleCurve, Curve, Curve3d as _, CylinderSurface, LineCurve, PlanarCurve,
    PlaneSurface, Surface as _, SurfaceGeometry,
};
use ogeom_math::{Circle, Cylinder, Direction, Frame, Plane, Point, Point2, Transform, Vector};
use ogeom_mesh::Deflection;
use ogeom_offset::{make_thick_sheet, offset_sheet};
use ogeom_topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, explore};
use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_2, PI};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, solid: &Shape) -> f64 {
    let deflection = Deflection {
        chord: 1e-4,
        ..Deflection::default()
    };
    ogeom_algo::volume_properties(model, solid, deflection, T)
        .unwrap()
        .mass
}

fn assert_valid(model: &Model, shape: &Shape) {
    let diagnosis = check(model, shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
}

fn z_range(model: &Model, shape: &Shape) -> (f64, f64) {
    let bounds = tight_bounds(model, shape, T).unwrap();
    (bounds.low().unwrap().z, bounds.high().unwrap().z)
}

/// A planar face in z = 0 over a polygon, with a hole for each inner
/// edge loop given.
fn planar_face(model: &mut Model, corners: &[Point], holes: &[Shape]) -> Shape {
    let wire = make_polygon(model, corners, true, T).unwrap().shape;
    let outer = explore(model, &wire, Filter::OfType(ShapeType::Edge)).unwrap();
    let mut wires = vec![outer];
    for hole in holes {
        wires.push(vec![hole.clone()]);
    }
    make_face_with_pcurves(
        model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &wires,
        T,
    )
    .unwrap()
    .shape
}

fn square(model: &mut Model, x0: f64, x1: f64) -> Shape {
    planar_face(
        model,
        &[
            Point::new(x0, 0.0, 0.0),
            Point::new(x1, 0.0, 0.0),
            Point::new(x1, 10.0, 0.0),
            Point::new(x0, 10.0, 0.0),
        ],
        &[],
    )
}

/// A quarter of the cylinder of radius 10 about z, from z = 0 to 5,
/// its normal pointing away from the axis.
fn quarter_cylinder(model: &mut Model) -> Shape {
    cylinder_arc(model, 10.0, (0.0, FRAC_PI_2), 5.0)
}

/// The part of the cylinder of `radius` about z between the angles
/// `from` and `to`, from z = 0 to `height`, its normal pointing away from
/// the axis.
fn cylinder_arc(model: &mut Model, radius: f64, (from, to): (f64, f64), height: f64) -> Shape {
    let corner = |angle: f64, z: f64| Point::new(radius * angle.cos(), radius * angle.sin(), z);
    let v = [
        make_vertex(model, corner(from, 0.0)).shape,
        make_vertex(model, corner(to, 0.0)).shape,
        make_vertex(model, corner(to, height)).shape,
        make_vertex(model, corner(from, height)).shape,
    ];
    let arc = |model: &mut Model, z: f64, a: &Shape, b: &Shape| {
        let frame = Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X, T).unwrap();
        let circle: Curve = CircleCurve::new(Circle::new(frame, radius, T).unwrap()).into();
        make_edge_between(model, circle, (from, to), a, b, T)
            .unwrap()
            .shape
    };
    let line = |model: &mut Model, a: &Shape, b: &Shape, p: Point, q: Point| {
        let segment: Curve = LineCurve::segment(p, q, T).unwrap().into();
        make_edge_between(model, segment, (0.0, p.distance(q)), a, b, T)
            .unwrap()
            .shape
    };
    let bottom = arc(model, 0.0, &v[0], &v[1]);
    let top = arc(model, height, &v[3], &v[2]);
    let left = line(model, &v[1], &v[2], corner(to, 0.0), corner(to, height));
    let right = line(model, &v[0], &v[3], corner(from, 0.0), corner(from, height));
    let cylinder = CylinderSurface::new(
        Cylinder::new(Frame::WORLD, radius, T).unwrap(),
        (-1.0, height + 1.0),
    )
    .unwrap();
    make_face_with_pcurves(
        model,
        cylinder.into(),
        &[vec![bottom, left, top.reversed(), right.reversed()]],
        T,
    )
    .unwrap()
    .shape
}

/// The saddle z = (x^2 - y^2) / 20 over [-5, 5]^2, fitted as a cubic
/// B-spline, bounded by its four border iso-curves.
fn saddle(model: &mut Model) -> (Shape, SurfaceGeometry) {
    saddle_with(model, None)
}

/// The saddle, with a hole where the chart circle of the given centre
/// and radius lies.
fn saddle_with(model: &mut Model, hole: Option<(Point2, f64)>) -> (Shape, SurfaceGeometry) {
    let n = 21;
    let rows: Vec<Vec<Point>> = (0..n)
        .map(|j| {
            let y = -5.0 + 10.0 * f64::from(j) / f64::from(n - 1);
            (0..n)
                .map(|i| {
                    let x = -5.0 + 10.0 * f64::from(i) / f64::from(n - 1);
                    Point::new(x, y, (x * x - y * y) / 20.0)
                })
                .collect()
        })
        .collect();
    let fitted = ogeom_geom::fit::fit_surface_grid(&rows, 3, 1e-6, T).unwrap();
    let surface: BSplineSurface = fitted.curve;
    let ((u0, u1), (v0, v1)) = surface.domain();
    let corners = [(u0, v0), (u1, v0), (u1, v1), (u0, v1)];
    let v: Vec<Shape> = corners
        .iter()
        .map(|(u, w)| make_vertex(model, surface.point_at(*u, *w, T).unwrap()).shape)
        .collect();
    let iso = |model: &mut Model, curve: ogeom_geom::BSplineCurve, from: &Shape, to: &Shape| {
        let range = curve.domain();
        make_edge_between(model, curve.into(), range, from, to, T)
            .unwrap()
            .shape
    };
    let south = iso(model, surface.iso_v_curve(v0, T).unwrap(), &v[0], &v[1]);
    let east = iso(model, surface.iso_u_curve(u1, T).unwrap(), &v[1], &v[2]);
    let north = iso(model, surface.iso_v_curve(v1, T).unwrap(), &v[3], &v[2]);
    let west = iso(model, surface.iso_u_curve(u0, T).unwrap(), &v[0], &v[3]);
    let geometry: SurfaceGeometry = surface.into();
    let mut wires = vec![vec![south, east, north.reversed(), west.reversed()]];
    if let Some((centre, radius)) = hole {
        let ring: PlanarCurve =
            ogeom_geom::Circle2d::new(ogeom_math::Circle2::centred(centre, radius, T).unwrap())
                .into();
        let traced: Curve = Curve::OnSurface(Box::new(ogeom_geom::CurveOnSurface::new(
            ring,
            geometry.clone(),
        )));
        let edge = make_edge(model, traced, (0.0, 2.0 * PI), T).unwrap().shape;
        wires.push(vec![edge.reversed()]);
    }
    let face = make_face_with_pcurves(model, geometry.clone(), &wires, T)
        .unwrap()
        .shape;
    (face, geometry)
}

#[test]
fn a_square_thickened_by_two_is_a_box_of_volume_200() {
    let mut model = Model::new();
    let face = square(&mut model, 0.0, 10.0);
    let solid = make_thick_sheet(&mut model, &face, 2.0, false, T)
        .unwrap()
        .shape;
    assert_valid(&model, &solid);
    assert!((volume(&model, &solid) - 200.0).abs() < 1e-9);
    let (lo, hi) = z_range(&model, &solid);
    assert!(
        lo.abs() < 1e-12 && (hi - 2.0).abs() < 1e-12,
        "z in [{lo}, {hi}]"
    );

    let both = make_thick_sheet(&mut model, &face, 2.0, true, T)
        .unwrap()
        .shape;
    assert_valid(&model, &both);
    assert!((volume(&model, &both) - 200.0).abs() < 1e-9);
    let (lo, hi) = z_range(&model, &both);
    assert!(
        (lo + 1.0).abs() < 1e-12 && (hi - 1.0).abs() < 1e-12,
        "z in [{lo}, {hi}]"
    );

    // Against the normal, the box hangs below the sheet.
    let under = make_thick_sheet(&mut model, &face, -2.0, false, T)
        .unwrap()
        .shape;
    assert_valid(&model, &under);
    assert!((volume(&model, &under) - 200.0).abs() < 1e-9);
    let (lo, hi) = z_range(&model, &under);
    assert!(
        (lo + 2.0).abs() < 1e-12 && hi.abs() < 1e-12,
        "z in [{lo}, {hi}]"
    );
}

#[test]
fn a_quarter_cylinder_thickens_to_its_annular_sector() {
    let mut model = Model::new();
    let face = quarter_cylinder(&mut model);
    let outward = make_thick_sheet(&mut model, &face, 1.0, false, T)
        .unwrap()
        .shape;
    assert_valid(&model, &outward);
    let want = (11.0_f64.powi(2) - 10.0_f64.powi(2)) * PI / 4.0 * 5.0;
    let got = volume(&model, &outward);
    assert!((got - want).abs() < want * 1e-6, "{got} against {want}");

    let inward = make_thick_sheet(&mut model, &face, -1.0, false, T)
        .unwrap()
        .shape;
    assert_valid(&model, &inward);
    let want = (10.0_f64.powi(2) - 9.0_f64.powi(2)) * PI / 4.0 * 5.0;
    let got = volume(&model, &inward);
    assert!((got - want).abs() < want * 1e-6, "{got} against {want}");
}

#[test]
fn an_offset_past_the_radius_of_curvature_is_refused() {
    let mut model = Model::new();
    let face = quarter_cylinder(&mut model);
    for distance in [-10.0, -12.0] {
        let refused = offset_sheet(&mut model, &face, distance, T).unwrap_err();
        assert!(
            refused.to_string().contains("radius of curvature"),
            "{refused}"
        );
    }
    // Short of the axis the moved face is the same quarter of the
    // cylinder of radius 1.
    let moved = offset_sheet(&mut model, &face, -9.0, T).unwrap().shape;
    assert_valid(&model, &moved);
    let area = ogeom_algo::surface_properties(&model, &moved, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((area - FRAC_PI_2 * 5.0).abs() < 1e-6, "area {area}");
}

#[test]
fn a_fitted_saddle_offset_by_one_stands_one_away_and_keeps_its_border() {
    let mut model = Model::new();
    let (face, original) = saddle(&mut model);
    let built = offset_sheet(&mut model, &face, 1.0, T).unwrap();
    let moved = built.shape.clone();
    assert_valid(&model, &moved);
    assert_eq!(built.history.modified(&face), std::slice::from_ref(&moved));

    // Everywhere on the moved face, one from the original surface.
    let Some(NodeData::Face(data)) = model.node(&moved).map(|n| n.data().clone()) else {
        panic!("the offset is not a face");
    };
    let surface = model.geometry().surface(data.surface).unwrap().clone();
    assert!(matches!(surface, SurfaceGeometry::BSpline(_)));
    let ((u0, u1), (v0, v1)) = original.domain();
    let mut deviation = 0.0_f64;
    for i in 0..=20 {
        for j in 0..=20 {
            let u = u0 + (u1 - u0) * f64::from(i) / 20.0;
            let v = v0 + (v1 - v0) * f64::from(j) / 20.0;
            let p = surface.point_at(u, v, T).unwrap();
            let foot = ogeom_algo::project_on_surface(&original, p, 32, T).unwrap();
            assert!(
                (foot.distance - 1.0).abs() < 1e-3,
                "{} from the saddle at ({u}, {v})",
                foot.distance
            );
            // And it is the point the original's normal reaches.
            let exact =
                original.point_at(u, v, T).unwrap() + original.normal_at(u, v, T).unwrap().vector();
            deviation = deviation.max(p.distance(exact));
        }
    }
    // The deviation is the face's to report, and it reports it.
    assert!(
        deviation <= data.tolerance.get(),
        "measured {deviation} beyond the face's tolerance {}",
        data.tolerance.get()
    );
    assert!(data.tolerance.get() <= 1e-5, "{}", data.tolerance.get());

    // Each border edge is its original moved one along the normal.
    let old_edges = ogeom_topo::explore_unique(&model, &face, ShapeType::Edge).unwrap();
    assert_eq!(old_edges.len(), 4);
    for old in &old_edges {
        let [new] = built.history.modified(old) else {
            panic!("an edge did not map to one edge");
        };
        let read = |shape: &Shape| {
            let data = model.node(shape).unwrap().data().as_edge().unwrap().clone();
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                panic!("no curve");
            };
            (model.geometry().curve(*curve).unwrap().clone(), *range)
        };
        let (before, _) = read(old);
        let (after, range) = read(new);
        for k in 0..=40 {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 40.0;
            let p = after.point_at(t, T).unwrap();
            let foot = ogeom_algo::project_on_surface(&original, p, 32, T).unwrap();
            assert!((foot.distance - 1.0).abs() < 1e-3, "{}", foot.distance);
            let on_border = ogeom_algo::project_on_curve(&before, foot.point, 64, T).unwrap();
            assert!(
                on_border.distance < 1e-3,
                "the border's offset lands {} off the border",
                on_border.distance
            );
        }
    }
}

#[test]
fn a_thickened_saddle_holds_the_volume_its_curvature_says() {
    // Both ways by half: the solid swept by the normal segment of
    // length t has volume, over the sheet, the integral of
    // t + t^3 K / 12 with K the Gaussian curvature.
    let mut model = Model::new();
    let (face, original) = saddle(&mut model);
    let t = 0.8;
    let solid = make_thick_sheet(&mut model, &face, t, true, T)
        .unwrap()
        .shape;
    assert_valid(&model, &solid);
    let want = over_saddle(&original, None, &|u, v| {
        let k = original.curvature_at(u, v, T).unwrap().gaussian();
        t + t.powi(3) * k / 12.0
    });
    let got = volume(&model, &solid);
    assert!((got - want).abs() < want * 1e-5, "{got} against {want}");
}

/// Over the saddle's chart, the integral of `f(u, v)` times the area
/// element, by Simpson's rule on an `n` by `n` grid; over the chart
/// disc of `hole` instead when one is given, in polar coordinates.
fn over_saddle(
    surface: &SurfaceGeometry,
    hole: Option<(Point2, f64)>,
    f: &dyn Fn(f64, f64) -> f64,
) -> f64 {
    let n = 80;
    let simpson = |k: i32| {
        if k == 0 || k == n {
            1.0
        } else if k % 2 == 1 {
            4.0
        } else {
            2.0
        }
    };
    let ((u0, u1), (v0, v1)) = surface.domain();
    let (a0, a1, b0, b1) = match hole {
        None => (u0, u1, v0, v1),
        Some((_, radius)) => (0.0, radius, 0.0, 2.0 * PI),
    };
    let (ha, hb) = ((a1 - a0) / f64::from(n), (b1 - b0) / f64::from(n));
    let mut sum = 0.0;
    for i in 0..=n {
        for j in 0..=n {
            let (a, b) = (a0 + ha * f64::from(i), b0 + hb * f64::from(j));
            let (u, v, jacobian) = match hole {
                None => (a, b, 1.0),
                Some((centre, _)) => (centre.x + a * b.cos(), centre.y + a * b.sin(), a),
            };
            let (su, sv) = surface.d1_at(u, v, T).unwrap();
            sum += simpson(i) * simpson(j) * su.cross(sv).magnitude() * jacobian * f(u, v);
        }
    }
    sum * ha * hb / 9.0
}

#[test]
fn a_trimmed_saddle_thickens_around_its_hole() {
    let hole = (Point2::new(0.45, 0.55), 0.2);
    let mut model = Model::new();
    let (face, original) = saddle_with(&mut model, Some(hole));
    let t = 0.6;
    let solid = make_thick_sheet(&mut model, &face, t, true, T)
        .unwrap()
        .shape;
    assert_valid(&model, &solid);
    let swept = |u: f64, v: f64| {
        let k = original.curvature_at(u, v, T).unwrap().gaussian();
        t + t.powi(3) * k / 12.0
    };
    let want = over_saddle(&original, None, &swept) - over_saddle(&original, Some(hole), &swept);
    let got = volume(&model, &solid);
    assert!((got - want).abs() < want * 1e-5, "{got} against {want}");

    // The hole's offset is the hole's edge moved one along the normal.
    let built = offset_sheet(&mut model, &face, 1.0, T).unwrap();
    assert_valid(&model, &built.shape);
    let ring = ogeom_topo::explore_unique(&model, &face, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| matches!(edge_vertices(&model, e).unwrap(), Some((a, b)) if a.node() == b.node()))
        .unwrap();
    let [moved] = built.history.modified(&ring) else {
        panic!("the hole's edge did not map to one edge");
    };
    let data = model.node(moved).unwrap().data().as_edge().unwrap().clone();
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        panic!("no curve");
    };
    let curve = model.geometry().curve(*curve).unwrap().clone();
    for k in 0..=40 {
        let s = range.0 + (range.1 - range.0) * f64::from(k) / 40.0;
        let p = curve.point_at(s, T).unwrap();
        let foot = ogeom_algo::project_on_surface(&original, p, 32, T).unwrap();
        assert!((foot.distance - 1.0).abs() < 1e-3, "{}", foot.distance);
        let (u, v) = foot.parameters;
        let off_ring = (Point2::new(u, v).distance(hole.0) - hole.1).abs();
        assert!(
            off_ring < 1e-4,
            "the hole's offset lands {off_ring} off the ring"
        );
    }
}

#[test]
fn a_square_with_a_hole_thickens_around_the_hole() {
    let mut model = Model::new();
    let frame = Frame::new(Point::new(5.0, 5.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let ring: Curve = CircleCurve::new(Circle::new(frame, 2.0, T).unwrap()).into();
    let hole = make_edge(&mut model, ring, (0.0, 2.0 * PI), T)
        .unwrap()
        .shape;
    let face = planar_face(
        &mut model,
        &[
            Point::new(0.0, 0.0, 0.0),
            Point::new(10.0, 0.0, 0.0),
            Point::new(10.0, 10.0, 0.0),
            Point::new(0.0, 10.0, 0.0),
        ],
        &[hole.reversed()],
    );
    let solid = make_thick_sheet(&mut model, &face, 2.0, false, T)
        .unwrap()
        .shape;
    assert_valid(&model, &solid);
    let want = (100.0 - PI * 4.0) * 2.0;
    let got = volume(&model, &solid);
    assert!((got - want).abs() < want * 1e-6, "{got} against {want}");
}

#[test]
fn a_shell_of_coplanar_faces_thickens_and_offsets_as_one() {
    let mut model = Model::new();
    let a = square(&mut model, 0.0, 10.0);
    let b = square(&mut model, 10.0, 25.0);
    let sewn = ogeom_algo::sew(&mut model, &[a, b], T).unwrap();
    let [shell] = sewn.shells.as_slice() else {
        panic!("two shells");
    };
    let moved = offset_sheet(&mut model, shell, 3.0, T).unwrap().shape;
    assert_eq!(model.kind_of(&moved).unwrap(), ShapeType::Shell);
    // Open along the same edges as the sheet, and nothing else.
    let before = check(&model, shell, T).unwrap();
    let after = check(&model, &moved, T).unwrap();
    assert!(after.is_usable(), "{after}");
    assert_eq!(
        after
            .problems
            .iter()
            .map(|p| p.what.clone())
            .collect::<Vec<_>>(),
        before
            .problems
            .iter()
            .map(|p| p.what.clone())
            .collect::<Vec<_>>()
    );
    // The shared edge stays shared.
    let edges = ogeom_topo::explore_unique(&model, &moved, ShapeType::Edge).unwrap();
    assert_eq!(edges.len(), 7);
    let (lo, hi) = z_range(&model, &moved);
    assert!((lo - 3.0).abs() < 1e-12 && (hi - 3.0).abs() < 1e-12);

    let solid = make_thick_sheet(&mut model, shell, 1.5, false, T)
        .unwrap()
        .shape;
    assert_valid(&model, &solid);
    assert!((volume(&model, &solid) - 250.0 * 1.5).abs() < 1e-9);
}

/// The square floor [0, 10]^2 in z = 0, facing up, sewn to a square wall
/// rising 10 from its edge at x = 10, bounded counter-clockwise about
/// its normal: back over the floor (-X) when `back`, else away (+X).
fn floor_and_wall(model: &mut Model, back: bool) -> Shape {
    let floor = square(model, 0.0, 10.0);
    let wall = {
        let mut corners = vec![
            Point::new(10.0, 0.0, 0.0),
            Point::new(10.0, 0.0, 10.0),
            Point::new(10.0, 10.0, 10.0),
            Point::new(10.0, 10.0, 0.0),
        ];
        let normal = if back {
            Direction::new(Vector::new(-1.0, 0.0, 0.0), T).unwrap()
        } else {
            corners.reverse();
            Direction::X
        };
        let wire = make_polygon(model, &corners, true, T).unwrap().shape;
        let edges = explore(model, &wire, Filter::OfType(ShapeType::Edge)).unwrap();
        let plane =
            Plane::new(Frame::new(Point::new(10.0, 0.0, 0.0), normal, Direction::Y, T).unwrap());
        make_face_with_pcurves(model, PlaneSurface::new(plane).into(), &[edges], T)
            .unwrap()
            .shape
    };
    let sewn = ogeom_algo::sew(model, &[floor, wall], T).unwrap();
    let [shell] = sewn.shells.as_slice() else {
        panic!("two shells");
    };
    shell.clone()
}

#[test]
fn a_floor_and_wall_thicken_with_a_mitre_but_do_not_offset() {
    let mut model = Model::new();
    // Walking their shared edge opposite ways, the floor and the wall face
    // into the corner between them.
    let shell = floor_and_wall(&mut model, true);
    let refused = offset_sheet(&mut model, &shell, 1.0, T).unwrap_err();
    assert!(refused.to_string().contains("crease"), "{refused}");
    // Into the corner the slabs meet on the mitre and their 1 x 1 x 10
    // overlap is counted once; away from it the mitre adds the 1 x 1 x 10
    // block at the corner. Half to each side, the two cancel.
    for (t, want) in [(1.0, 190.0), (-1.0, 210.0), (2.0, 400.0)] {
        let solid = make_thick_sheet(&mut model, &shell, t, t == 2.0, T)
            .unwrap()
            .shape;
        assert_valid(&model, &solid);
        let got = volume(&model, &solid);
        assert!((got - want).abs() < 1e-9, "{t}: {got} against {want}");
    }
}

#[test]
fn creases_the_mitre_does_not_join_are_refused() {
    let mut model = Model::new();
    // Facing away from the floor, the wall walks the edge they share the
    // way the floor does: they disagree on which side of the sheet is
    // which.
    let shell = floor_and_wall(&mut model, false);
    let refused = make_thick_sheet(&mut model, &shell, 1.0, false, T).unwrap_err();
    assert!(
        refused.to_string().contains("walk a crease the same way"),
        "{refused}"
    );
    // Swept aslant, the folded sheet's borders do not lie square to its
    // crease, and their offsets miss the mitre.
    let wire = make_polygon(
        &mut model,
        &[
            Point::ORIGIN,
            Point::new(10.0, 0.0, 0.0),
            Point::new(10.0, 10.0, 0.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    let aslant = ogeom_algo::make_prism(&mut model, &wire, Vector::new(2.0, 1.0, 5.0), T)
        .unwrap()
        .shape;
    let refused = make_thick_sheet(&mut model, &aslant, 1.0, false, T).unwrap_err();
    assert!(refused.to_string().contains("mitre"), "{refused}");
}

#[test]
fn solids_closed_shells_and_empty_moves_are_refused() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T)
        .unwrap()
        .shape;
    let refused = offset_sheet(&mut model, &solid, 1.0, T).unwrap_err();
    assert!(refused.to_string().contains("offset_shape"), "{refused}");
    let shell = explore(&model, &solid, Filter::OfType(ShapeType::Shell)).unwrap()[0].clone();
    let refused = make_thick_sheet(&mut model, &shell, 1.0, false, T).unwrap_err();
    assert!(refused.to_string().contains("no free edge"), "{refused}");

    let face = square(&mut model, 0.0, 10.0);
    assert!(offset_sheet(&mut model, &face, 0.0, T).is_err());
    assert!(make_thick_sheet(&mut model, &face, f64::NAN, true, T).is_err());
}

/// How many faces of a shape stand on each family of surface.
fn families(model: &Model, shape: &Shape) -> HashMap<&'static str, usize> {
    let mut out = HashMap::new();
    for face in ogeom_topo::explore_unique(model, shape, ShapeType::Face).unwrap() {
        let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data()) else {
            panic!("not a face");
        };
        let name = match model.geometry().surface(data.surface).unwrap() {
            SurfaceGeometry::Plane(_) => "plane",
            SurfaceGeometry::Cylinder(_) => "cylinder",
            SurfaceGeometry::Sphere(_) => "sphere",
            SurfaceGeometry::BSpline(_) => "bspline",
            _ => "other",
        };
        *out.entry(name).or_insert(0) += 1;
    }
    out
}

#[test]
fn side_faces_are_exact_where_the_rulings_allow() {
    let mut model = Model::new();
    let face = square(&mut model, 0.0, 10.0);
    let solid = make_thick_sheet(&mut model, &face, 2.0, false, T)
        .unwrap()
        .shape;
    assert_eq!(families(&model, &solid), HashMap::from([("plane", 6)]));

    // Along the quarter cylinder's straight edges the sides are planes;
    // along its arcs the rulings run radially, an annular sector the
    // edge's angle does not chart as a plane, so they are fitted.
    let face = quarter_cylinder(&mut model);
    let solid = make_thick_sheet(&mut model, &face, 1.0, false, T)
        .unwrap()
        .shape;
    assert_eq!(
        families(&model, &solid),
        HashMap::from([("cylinder", 2), ("plane", 2), ("bspline", 2)])
    );
}

#[test]
fn a_drum_wall_offsets_and_thickens_across_its_seam() {
    let mut model = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 4.0, 6.0, T)
        .unwrap()
        .shape;
    let wall = ogeom_topo::explore_unique(&model, &drum, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let Some(NodeData::Face(data)) = model.node(f).map(|n| n.data()) else {
                return false;
            };
            matches!(
                model.geometry().surface(data.surface),
                Some(SurfaceGeometry::Cylinder(_))
            )
        })
        .unwrap();
    let moved = offset_sheet(&mut model, &wall, 1.5, T).unwrap().shape;
    let area = ogeom_algo::surface_properties(&model, &moved, Deflection::default(), T)
        .unwrap()
        .mass;
    let want = 2.0 * PI * 5.5 * 6.0;
    assert!((area - want).abs() < want * 1e-6, "{area} against {want}");

    let tube = make_thick_sheet(&mut model, &wall, 1.0, false, T)
        .unwrap()
        .shape;
    assert_valid(&model, &tube);
    let want = PI * (25.0 - 16.0) * 6.0;
    let got = volume(&model, &tube);
    assert!((got - want).abs() < want * 1e-6, "{got} against {want}");
}

#[test]
fn a_ball_offsets_through_its_poles() {
    let mut model = Model::new();
    let ball = ogeom_algo::make_sphere(&mut model, Frame::WORLD, 3.0, T)
        .unwrap()
        .shape;
    let shell = explore(&model, &ball, Filter::OfType(ShapeType::Shell)).unwrap()[0].clone();
    let grown = offset_sheet(&mut model, &shell, 2.0, T).unwrap().shape;
    assert_eq!(families(&model, &grown), HashMap::from([("sphere", 1)]));
    let area = ogeom_algo::surface_properties(&model, &grown, Deflection::default(), T)
        .unwrap()
        .mass;
    let want = 4.0 * PI * 25.0;
    assert!((area - want).abs() < want * 1e-6, "{area} against {want}");
    let refused = offset_sheet(&mut model, &shell, -3.0, T).unwrap_err();
    assert!(
        refused.to_string().contains("radius of curvature"),
        "{refused}"
    );
}

#[test]
fn a_placed_sheet_thickens_where_it_stands() {
    let mut model = Model::new();
    let face = square(&mut model, 0.0, 10.0);
    let lifted = ogeom_algo::transformed(
        &mut model,
        &face,
        Transform::translation(Vector::new(0.0, 0.0, 7.0)),
    )
    .unwrap()
    .shape;
    let solid = make_thick_sheet(&mut model, &lifted, 2.0, false, T)
        .unwrap()
        .shape;
    assert_valid(&model, &solid);
    assert!((volume(&model, &solid) - 200.0).abs() < 1e-9);
    let (lo, hi) = z_range(&model, &solid);
    assert!(
        (lo - 7.0).abs() < 1e-12 && (hi - 9.0).abs() < 1e-12,
        "z in [{lo}, {hi}]"
    );
    // The history runs from the placed face to its offset.
    let moved = offset_sheet(&mut model, &lifted, 1.0, T).unwrap();
    assert_eq!(
        moved.history.modified(&lifted),
        std::slice::from_ref(&moved.shape)
    );
    let (lo, hi) = z_range(&model, &moved.shape);
    assert!((lo - 8.0).abs() < 1e-12 && (hi - 8.0).abs() < 1e-12);
}

/// A plane face through the given corners, its own frame read off them.
fn face_through(model: &mut Model, corners: [Point; 4]) -> Shape {
    let wire = make_polygon(model, &corners, true, T).unwrap().shape;
    let edges = explore(model, &wire, Filter::OfType(ShapeType::Edge)).unwrap();
    let x = Direction::new(corners[1] - corners[0], T).unwrap();
    let z = Direction::new((corners[1] - corners[0]).cross(corners[3] - corners[0]), T).unwrap();
    let plane = Plane::new(Frame::new(corners[0], z, x, T).unwrap());
    make_face_with_pcurves(model, PlaneSurface::new(plane).into(), &[edges], T)
        .unwrap()
        .shape
}

#[test]
fn layers_that_run_into_each_other_are_refused() {
    // Three quarters of a drum of radius 5, with a flat flap running on
    // tangentially from each end: the flaps converge, x = 5 and y = -5,
    // and stop 4 along, short of each other.
    let mut model = Model::new();
    let drum = cylinder_arc(&mut model, 5.0, (0.0, 1.5 * PI), 3.0);
    let east = face_through(
        &mut model,
        [
            Point::new(5.0, -4.0, 0.0),
            Point::new(5.0, 0.0, 0.0),
            Point::new(5.0, 0.0, 3.0),
            Point::new(5.0, -4.0, 3.0),
        ],
    );
    let south = face_through(
        &mut model,
        [
            Point::new(0.0, -5.0, 0.0),
            Point::new(4.0, -5.0, 0.0),
            Point::new(4.0, -5.0, 3.0),
            Point::new(0.0, -5.0, 3.0),
        ],
    );
    let sewn = ogeom_algo::sew(&mut model, &[drum, east, south], T).unwrap();
    let [sheet] = sewn.shells.as_slice() else {
        panic!("the pieces did not sew into one sheet");
    };
    // Half a unit toward the axis the walls stand apart: two slabs and a
    // sector of the annulus between radii 4.5 and 5.
    let thin = make_thick_sheet(&mut model, sheet, -0.5, false, T).unwrap();
    assert_valid(&model, &thin.shape);
    let want = 2.0 * 4.0 * 0.5 * 3.0 + (25.0 - 4.5 * 4.5) / 2.0 * 1.5 * PI * 3.0;
    let got = volume(&model, &thin.shape);
    assert!((got - want).abs() < want * 1e-6, "{got} against {want}");
    // Two units toward it the flaps' layers cross.
    let refused = make_thick_sheet(&mut model, sheet, -2.0, false, T).unwrap_err();
    assert!(
        refused.to_string().contains("runs into itself"),
        "{refused}"
    );
}

#[test]
fn an_edge_shared_by_three_faces_is_refused() {
    let mut model = Model::new();
    let corner = |x: f64, y: f64, z: f64| Point::new(x, y, z);
    let a = make_vertex(&mut model, corner(0.0, 0.0, 0.0)).shape;
    let b = make_vertex(&mut model, corner(0.0, 10.0, 0.0)).shape;
    let spine: Curve = LineCurve::segment(corner(0.0, 0.0, 0.0), corner(0.0, 10.0, 0.0), T)
        .unwrap()
        .into();
    let shared = make_edge_between(&mut model, spine, (0.0, 10.0), &a, &b, T)
        .unwrap()
        .shape;
    // Three fins on the one edge: two floors and a wall.
    let mut fins = Vec::new();
    for far in [
        Vector::new(10.0, 0.0, 0.0),
        Vector::new(-10.0, 0.0, 0.0),
        Vector::new(0.0, 0.0, 10.0),
    ] {
        let p = corner(0.0, 0.0, 0.0) + far;
        let q = corner(0.0, 10.0, 0.0) + far;
        let c = make_vertex(&mut model, q).shape;
        let d = make_vertex(&mut model, p).shape;
        let line = |model: &mut Model, from: &Shape, to: &Shape, s: Point, e: Point| {
            let segment: Curve = LineCurve::segment(s, e, T).unwrap().into();
            make_edge_between(model, segment, (0.0, s.distance(e)), from, to, T)
                .unwrap()
                .shape
        };
        let bc = line(&mut model, &b, &c, corner(0.0, 10.0, 0.0), q);
        let cd = line(&mut model, &c, &d, q, p);
        let da = line(&mut model, &d, &a, p, corner(0.0, 0.0, 0.0));
        let z = Direction::new(Vector::new(0.0, 10.0, 0.0).cross(far), T).unwrap();
        let plane = Plane::new(Frame::new(Point::ORIGIN, z, Direction::Y, T).unwrap());
        fins.push(
            make_face_with_pcurves(
                &mut model,
                PlaneSurface::new(plane).into(),
                &[vec![shared.clone(), bc, cd, da]],
                T,
            )
            .unwrap()
            .shape,
        );
    }
    let shell = ogeom_algo::make_shell(&mut model, &fins).unwrap().shape;
    let refused = offset_sheet(&mut model, &shell, 1.0, T).unwrap_err();
    assert!(
        refused.to_string().contains("more than two faces"),
        "{refused}"
    );
}

#[test]
fn a_bare_patch_offsets_but_does_not_thicken() {
    let mut model = Model::new();
    let plane = PlaneSurface::over(Plane::new(Frame::WORLD), (0.0, 4.0), (0.0, 3.0)).unwrap();
    let patch = ogeom_algo::make_natural_face(&mut model, plane.into())
        .unwrap()
        .shape;
    let moved = offset_sheet(&mut model, &patch, 2.0, T).unwrap().shape;
    let area = ogeom_algo::surface_properties(&model, &moved, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((area - 12.0).abs() < 1e-9, "area {area}");
    let refused = make_thick_sheet(&mut model, &patch, 1.0, false, T).unwrap_err();
    assert!(
        refused.to_string().contains("no boundary edges"),
        "{refused}"
    );
}

#[test]
fn a_free_form_patch_without_a_normal_is_refused() {
    // A bilinear patch with one side collapsed to a point: its normal is
    // undefined along that side, and so is its offset there.
    let mut model = Model::new();
    let knots = || ogeom_math::KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
    let grid = ogeom_math::ControlGrid::new(
        vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(0.0, 0.0, 0.0),
            Point::new(10.0, 0.0, 1.0),
            Point::new(10.0, 10.0, -1.0),
        ],
        2,
        2,
    )
    .unwrap();
    let patch = BSplineSurface::new(knots(), knots(), &grid, T).unwrap();
    let face = ogeom_algo::make_natural_face(&mut model, patch.into())
        .unwrap()
        .shape;
    let refused = offset_sheet(&mut model, &face, 1.0, T).unwrap_err();
    assert!(refused.to_string().contains("no normal"), "{refused}");
}

/// The quarter of a ball of radius 3 between the meridians at longitude 0
/// and a quarter turn, pole to pole: two meridians and two poles.
fn lune(model: &mut Model) -> Shape {
    use ogeom_algo::attach_pcurve;
    use ogeom_geom::{Line2d, SphereSurface};
    use ogeom_math::{Axis2, Direction2, Sphere};
    use ogeom_topo::{EdgeData, Location};
    let south = make_vertex(model, Point::new(0.0, 0.0, -3.0)).shape;
    let north = make_vertex(model, Point::new(0.0, 0.0, 3.0)).shape;
    let meridian = |model: &mut Model, longitude: f64| {
        let x = Direction::new(Vector::new(longitude.cos(), longitude.sin(), 0.0), T).unwrap();
        let y = Direction::new(Vector::new(-longitude.sin(), longitude.cos(), 0.0), T).unwrap();
        let circle = Circle::new(
            Frame::new(Point::ORIGIN, y.reversed(), x, T).unwrap(),
            3.0,
            T,
        )
        .unwrap();
        make_edge_between(
            model,
            CircleCurve::new(circle).into(),
            (-FRAC_PI_2, FRAC_PI_2),
            &south,
            &north,
            T,
        )
        .unwrap()
        .shape
    };
    let west = meridian(model, 0.0);
    let east = meridian(model, FRAC_PI_2);
    let pole = |model: &mut Model, at: &Shape| {
        let mut data = EdgeData::new();
        data.degenerate = true;
        model.add_edge(data, &[at.clone(), at.clone()]).unwrap()
    };
    let bottom = pole(model, &south);
    let top = pole(model, &north);
    let surface = model
        .geometry_mut()
        .add_surface(SphereSurface::new(Sphere::new(Frame::WORLD, 3.0, T).unwrap()).into());
    let wire = ogeom_algo::make_wire(
        model,
        &[
            bottom.clone(),
            east.clone(),
            top.reversed(),
            west.reversed(),
        ],
        T,
    )
    .unwrap()
    .shape;
    let face = ogeom_algo::make_face_on(model, surface, &[wire], T)
        .unwrap()
        .shape;
    let line = |from: Point2, along: Direction2, length: f64| -> PlanarCurve {
        Line2d::over(Axis2::new(from, along), 0.0, length)
            .unwrap()
            .into()
    };
    let up = |longitude: f64| -> PlanarCurve {
        Line2d::over(
            Axis2::new(Point2::new(longitude, 0.0), Direction2::Y),
            -FRAC_PI_2,
            FRAC_PI_2,
        )
        .unwrap()
        .into()
    };
    let identity = Location::identity();
    let across = |v: f64| line(Point2::new(0.0, v), Direction2::X, FRAC_PI_2);
    attach_pcurve(
        model,
        &bottom,
        across(-FRAC_PI_2),
        surface,
        identity.clone(),
        (0.0, FRAC_PI_2),
    )
    .unwrap();
    attach_pcurve(
        model,
        &top,
        across(FRAC_PI_2),
        surface,
        identity.clone(),
        (0.0, FRAC_PI_2),
    )
    .unwrap();
    attach_pcurve(
        model,
        &west,
        up(0.0),
        surface,
        identity.clone(),
        (-FRAC_PI_2, FRAC_PI_2),
    )
    .unwrap();
    attach_pcurve(
        model,
        &east,
        up(FRAC_PI_2),
        surface,
        identity,
        (-FRAC_PI_2, FRAC_PI_2),
    )
    .unwrap();
    face
}

#[test]
fn a_lune_offsets_through_its_poles_but_its_sides_do_not_close() {
    let mut model = Model::new();
    let face = lune(&mut model);
    // A lune of dihedral angle a on a ball of radius r has area 2 a r^2.
    let grown = offset_sheet(&mut model, &face, 1.0, T).unwrap().shape;
    assert_valid(&model, &grown);
    let area = ogeom_algo::surface_properties(&model, &grown, Deflection::default(), T)
        .unwrap()
        .mass;
    let want = 2.0 * FRAC_PI_2 * 16.0;
    assert!((area - want).abs() < want * 1e-6, "{area} against {want}");
    // Its free meridians end at the poles, where the sides have no normal
    // to rule along.
    let refused = make_thick_sheet(&mut model, &face, 1.0, false, T).unwrap_err();
    assert!(refused.to_string().contains("no normal"), "{refused}");
}

/// The sheet the open polyline `corners` in z = 0 sweeps rising 5.
fn folded(model: &mut Model, corners: &[Point]) -> Shape {
    let wire = make_polygon(model, corners, false, T).unwrap().shape;
    ogeom_algo::make_prism(model, &wire, Vector::new(0.0, 0.0, 5.0), T)
        .unwrap()
        .shape
}

fn corners(model: &Model, shape: &Shape) -> Vec<Point> {
    explore(model, shape, Filter::OfType(ShapeType::Vertex))
        .unwrap()
        .iter()
        .map(|v| match model.node(v).unwrap().data() {
            NodeData::Vertex(data) => data.point,
            _ => panic!("not a vertex"),
        })
        .collect()
}

/// Thickened by `t` (by half of it to each side when `both`), the sheet is
/// one valid solid of volume `want`.
fn thickens_to(model: &mut Model, sheet: &Shape, t: f64, both: bool, want: f64) -> Shape {
    let solid = make_thick_sheet(model, sheet, t, both, T).unwrap().shape;
    assert_eq!(model.kind_of(&solid).unwrap(), ShapeType::Solid);
    assert_valid(model, &solid);
    let got = volume(model, &solid);
    assert!(
        (got - want).abs() < want * 1e-9,
        "{t}: {got} against {want}"
    );
    solid
}

#[test]
fn a_sheet_folded_square_thickens_with_a_mitre() {
    let mut model = Model::new();
    let sheet = folded(
        &mut model,
        &[
            Point::ORIGIN,
            Point::new(10.0, 0.0, 0.0),
            Point::new(10.0, 10.0, 0.0),
        ],
    );
    // The faces face away from the corner, right of the polyline's
    // travel: that way the slabs part and the mitre fills the 1 x 1 square
    // at the corner, the other way they overlap on it.
    let out = thickens_to(&mut model, &sheet, 1.0, false, 5.0 * 21.0);
    let found = corners(&model, &out);
    for z in [0.0, 5.0] {
        let mitre = Point::new(11.0, -1.0, z);
        assert!(found.iter().any(|c| c.distance(mitre) < 1e-12), "{mitre:?}");
    }
    // The sheet's two faces, their offsets and a side along each of the
    // six free edges.
    assert_eq!(
        explore(&model, &out, Filter::OfType(ShapeType::Face))
            .unwrap()
            .len(),
        10
    );
    let inside = thickens_to(&mut model, &sheet, -1.0, false, 5.0 * 19.0);
    let found = corners(&model, &inside);
    for z in [0.0, 5.0] {
        let mitre = Point::new(9.0, 1.0, z);
        assert!(found.iter().any(|c| c.distance(mitre) < 1e-12), "{mitre:?}");
    }
    // Half to each side: the square gained outside is the one lost inside.
    thickens_to(&mut model, &sheet, 2.0, true, 5.0 * 40.0);
}

#[test]
fn a_sheet_folded_at_sixty_degrees_thickens_with_a_mitre() {
    let mut model = Model::new();
    let turn = PI / 3.0;
    let sheet = folded(
        &mut model,
        &[
            Point::ORIGIN,
            Point::new(10.0, 0.0, 0.0),
            Point::new(10.0 + 10.0 * turn.cos(), 10.0 * turn.sin(), 0.0),
        ],
    );
    // A band of width t along a polyline of length 20 turning by 60
    // degrees gains t^2 tan 30 at the mitre outside the turn and loses as
    // much inside it.
    let corner = (turn / 2.0).tan();
    thickens_to(&mut model, &sheet, 1.0, false, 5.0 * (20.0 + corner));
    thickens_to(&mut model, &sheet, -1.0, false, 5.0 * (20.0 - corner));
    thickens_to(&mut model, &sheet, 2.0, true, 5.0 * 40.0);
}

#[test]
fn a_line_turning_into_an_arc_thickens_with_a_mitre() {
    let mut model = Model::new();
    // The line along +X to (10, 0, 0), then a quarter of the circle of
    // radius 5 about (5, 0, 0) from there to (5, 5, 0): a corner of 90
    // degrees between a plane and a cylinder.
    let v = [
        make_vertex(&mut model, Point::ORIGIN).shape,
        make_vertex(&mut model, Point::new(10.0, 0.0, 0.0)).shape,
        make_vertex(&mut model, Point::new(5.0, 5.0, 0.0)).shape,
    ];
    let segment: Curve = LineCurve::segment(Point::ORIGIN, Point::new(10.0, 0.0, 0.0), T)
        .unwrap()
        .into();
    let line = make_edge_between(&mut model, segment, (0.0, 10.0), &v[0], &v[1], T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(5.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let circle: Curve = CircleCurve::new(Circle::new(frame, 5.0, T).unwrap()).into();
    let arc = make_edge_between(&mut model, circle, (0.0, FRAC_PI_2), &v[1], &v[2], T)
        .unwrap()
        .shape;
    let wire = ogeom_algo::make_wire(&mut model, &[line, arc], T)
        .unwrap()
        .shape;
    let sheet = ogeom_algo::make_prism(&mut model, &wire, Vector::new(0.0, 0.0, 5.0), T)
        .unwrap()
        .shape;
    // Outside the turn the line's offset y = -1 meets the circle of radius
    // 6 at M = (5 + sqrt 35, -1), at the angle -asin(1/6) about the
    // centre. The band is the trapezoid under the line, out to M, and the
    // annular sector between radii 5 and 6 with the sliver of the radius-6
    // disc between the angles -asin(1/6) and 0, less the triangle from the
    // centre to the corner and M.
    let r35 = 35.0_f64.sqrt();
    let outside = (15.0 + r35) / 2.0 + 11.0 * PI / 4.0 + 18.0 * (1.0 / 6.0_f64).asin() - 2.5;
    let solid = thickens_to(&mut model, &sheet, 1.0, false, 5.0 * outside);
    let mitre = Point::new(5.0 + r35, -1.0, 0.0);
    assert!(
        corners(&model, &solid)
            .iter()
            .any(|c| c.distance(mitre) < 1e-9)
    );
    // Inside, y = 1 meets the circle of radius 4 at (5 + sqrt 15, 1), at
    // the angle a = asin(1/4): the trapezoid, the annular sector between
    // radii 4 and 5 from a to 90 degrees, and the radius-5 sector from 0
    // to a less the triangle from the centre to the corner and the mitre.
    let a = 0.25_f64.asin();
    let inside = (15.0 + 15.0_f64.sqrt()) / 2.0 + 4.5 * (FRAC_PI_2 - a) + 12.5 * a - 2.5;
    thickens_to(&mut model, &sheet, -1.0, false, 5.0 * inside);
}
