//! A tool meeting a face of the solid along a curve lying in it: an edge
//! of the tool in the face, or a section through a tool face whose chart
//! closes to one point along a whole side.

#![allow(clippy::unwrap_used)]

use core::f64::consts::PI;

use ogeom_core::Tolerances;
use ogeom_geom::{BSplineCurve, BSplineSurface, Curve, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{ControlGrid, Direction, Frame, KnotVector, Plane, Point};
use ogeom_topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    ogeom_algo::volume_properties(model, shape, fine, T)
        .unwrap()
        .mass
}

/// `cut` and `common` of `a` by `b` each check valid and hold their closed
/// form within `within`.
fn both_hold(model: &mut Model, a: &Shape, b: &Shape, inside: f64, within: f64) {
    let whole = volume(model, a);
    let cut = ogeom_bool::cut(model, a, b, T).unwrap().shape;
    let common = ogeom_bool::common(model, a, b, T).unwrap().shape;
    for (what, shape, want) in [("cut", &cut, whole - inside), ("common", &common, inside)] {
        let found = ogeom_algo::check(model, shape, T).unwrap();
        assert!(found.is_valid(), "{what}: {found}");
        let got = volume(model, shape);
        eprintln!("{what}: {got:.12} against {want:.12}");
        assert!((got - want).abs() < within, "{what}: {got} against {want}");
    }
}

/// A solid with its apex at `p`, over the base `z = 0` bounded by the
/// quadratic arc `b0 m b1` and the lines from its ends to `q`. Its side
/// over the arc is the ruled patch from the arc (`u = 0`) to `p` (`u = 1`):
/// the whole side `u = 1` of its chart is the one point `p`, closed by a
/// degenerate edge running along `v`.
fn closing_patch_solid(model: &mut Model, [b0, m, b1, q, p]: [Point; 5]) -> Shape {
    let k2 = KnotVector::clamped_uniform(2, 3).unwrap();
    let k1 = KnotVector::clamped_uniform(1, 2).unwrap();
    let grid = ControlGrid::new(vec![b0, m, b1, p, p, p], 2, 3).unwrap();
    let patch: SurfaceGeometry = BSplineSurface::new(k1, k2.clone(), &grid, T)
        .unwrap()
        .into();
    let arc: Curve = BSplineCurve::new(k2, vec![b0, m, b1], T).unwrap().into();
    let v: Vec<Shape> = [b0, b1, q, p]
        .iter()
        .map(|x| ogeom_algo::make_vertex(model, *x).shape)
        .collect();
    let line = |model: &mut Model, i: usize, j: usize| {
        let (a, b) = ([b0, b1, q, p][i], [b0, b1, q, p][j]);
        let c: Curve = LineCurve::segment(a, b, T).unwrap().into();
        ogeom_algo::make_edge_between(model, c, (0.0, a.distance(b)), &v[i], &v[j], T)
            .unwrap()
            .shape
    };
    let e_arc = ogeom_algo::make_edge_between(model, arc, (0.0, 1.0), &v[0], &v[1], T)
        .unwrap()
        .shape;
    let (e_b0p, e_b1p, e_b0q, e_b1q, e_qp) = (
        line(model, 0, 3),
        line(model, 1, 3),
        line(model, 0, 2),
        line(model, 1, 2),
        line(model, 2, 3),
    );
    let mut faces = Vec::new();
    let ring = [e_b0p.clone(), e_b1p.reversed(), e_arc.reversed()];
    let wire = ogeom_algo::make_wire(model, &ring, T).unwrap().shape;
    let bare = ogeom_algo::make_face(model, patch, core::slice::from_ref(&wire), T)
        .unwrap()
        .shape;
    ogeom_algo::trimmed_or_fitted_where_bare(model, &bare, T).unwrap();
    let surface = model.node(&bare).unwrap().data().as_face().unwrap().surface;
    let uses = model.children_of(&wire).unwrap();
    let closed = ogeom_algo::closed_at_poles(model, surface, &uses, T).unwrap();
    assert_eq!(closed.len(), uses.len() + 1, "the apex's degenerate edge");
    let wire = ogeom_algo::make_wire(model, &closed, T).unwrap().shape;
    let side = ogeom_algo::make_face_on(model, surface, &[wire], T)
        .unwrap()
        .shape;
    faces.push(side.reversed());
    let centre = Point::new(
        (b0.x + b1.x + q.x + p.x) / 4.0,
        (b0.y + b1.y + q.y + p.y) / 4.0,
        (b0.z + b1.z + q.z + p.z) / 4.0,
    );
    // Each plane faces away from the solid.
    let plane = |a: Point, b: Point, c: Point| -> SurfaceGeometry {
        let (b, c) = if (b - a).cross(c - a).dot(a - centre) > 0.0 {
            (b, c)
        } else {
            (c, b)
        };
        PlaneSurface::new(Plane::through_points(a, b, c, T).unwrap()).into()
    };
    for (ring, surface) in [
        (
            vec![e_arc, e_b1q.clone(), e_b0q.reversed()],
            plane(b0, b1, q),
        ),
        (
            vec![e_b0p, e_qp.reversed(), e_b0q.reversed()],
            plane(b0, p, q),
        ),
        (
            vec![e_b1p, e_qp.reversed(), e_b1q.reversed()],
            plane(b1, p, q),
        ),
    ] {
        let wire = ogeom_algo::make_wire(model, &ring, T).unwrap().shape;
        let face = ogeom_algo::make_face(model, surface, &[wire], T)
            .unwrap()
            .shape;
        ogeom_algo::trimmed_or_fitted_where_bare(model, &face, T).unwrap();
        faces.push(face);
    }
    let shell = ogeom_algo::make_shell(model, &faces).unwrap().shape;
    let solid = ogeom_algo::make_solid(model, &[shell]).unwrap().shape;
    let found = ogeom_algo::check(model, &solid, T).unwrap();
    assert!(found.is_valid(), "the closing patch's solid: {found}");
    solid
}

/// The plane `z = h` crosses the closing patch along the line `u = h / 4`
/// of its chart, beside the side `u = 1` that closes to the apex. A
/// containment ray from that line runs across the closed side, which the
/// face's outline left out: the section read as outside the patch, split
/// nothing, and the kept pieces did not close. The solid is a cone over
/// its base `B` (area `32 / 3`: the triangle `b0 b1 q` of 8 and the arc's
/// segment of `2/3` of the triangle `b0 m b1`), apex 4 above it, and the
/// part above `z = h` is it scaled by `s = 1 - h / 4` about the apex.
#[test]
fn a_box_crosses_a_face_whose_chart_closes_to_a_point_along_a_side() {
    let corners = [
        Point::new(-2.0, -1.0, 0.0),
        Point::new(0.0, -3.0, 0.0),
        Point::new(2.0, -1.0, 0.0),
        Point::new(0.0, 3.0, 0.0),
        Point::new(0.0, 0.0, 4.0),
    ];
    let whole = 32.0 / 3.0 * 4.0 / 3.0;
    for h in [1.0, 2.0, 3.0] {
        let mut model = Model::new();
        let tool = closing_patch_solid(&mut model, corners);
        assert!((volume(&model, &tool) - whole).abs() < 1e-9);
        let block = ogeom_algo::make_box(&mut model, at(-5.0, -5.0, h), (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let s: f64 = 1.0 - h / 4.0;
        both_hold(&mut model, &block, &tool, whole * s.powi(3), 1e-9);
    }
}

/// A capsule (a cylinder of radius 2 along `x` closed by a half ball) whose
/// joining circle lies in the wall `x = 10` of a block 20 by 10 by 10
/// notched over `x > 10, y > 5`: half the circle lies in the wall, half in
/// the material below the notch, and it crosses the wall's edge. The tool
/// edge in the face splits it, and the capsule keeps its cylinder over
/// `0 < x < 10` and the quarter ball under the notch.
#[test]
fn a_tool_edge_lying_in_a_face_and_across_its_edge_splits_it() {
    let mut model = Model::new();
    let block = ogeom_algo::make_box(&mut model, at(0.0, 0.0, 0.0), (20.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let notch = ogeom_algo::make_box(&mut model, at(10.0, 5.0, -1.0), (11.0, 6.0, 12.0), T)
        .unwrap()
        .shape;
    let block = ogeom_bool::cut(&mut model, &block, &notch, T)
        .unwrap()
        .shape;
    let axis = Frame::about(Point::new(-5.0, 5.0, 5.0), Direction::X);
    let drum = ogeom_algo::make_cylinder(&mut model, axis, 2.0, 15.0, T)
        .unwrap()
        .shape;
    let ball = Frame::about(Point::new(10.0, 5.0, 5.0), Direction::X);
    let ball = ogeom_algo::make_sphere(&mut model, ball, 2.0, T)
        .unwrap()
        .shape;
    let capsule = ogeom_bool::fuse(&mut model, &drum, &ball, T).unwrap().shape;
    both_hold(
        &mut model,
        &block,
        &capsule,
        40.0 * PI + 8.0 * PI / 3.0,
        1e-6,
    );
}
