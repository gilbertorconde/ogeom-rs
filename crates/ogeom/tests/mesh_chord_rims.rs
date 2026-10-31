//! A face whose rim runs along straight chords of its curved surface, as a
//! solid converted back from a mesh has them: each chord an edge whose
//! tolerance covers how far it stands inside the surface. The face's mesh
//! says whether it met the chord it was asked for, and says so truly.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{attach_pcurve, make_edge_between, make_face_on, make_vertex, make_wire};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve3d as _, CylinderSurface, Line2d, LineCurve, Surface as _};
use ogeom::math::{Circle, Cylinder, Direction, Frame, Point, Point2};
use ogeom::mesh::Deflection;
use ogeom::topo::{Location, Model, NodeData, Orientation, Shape};

const T: Tolerances = Tolerances::millimetres();

const RADIUS: f64 = 2.0;

/// A cylinder face of radius 2 shaped as one of Voron 2.4 solid 990's
/// converted back from its mesh: a quarter turn high three, its lower left
/// corner cut off by two straight chords a fifth of a radian each across
/// the turn, which stand up to 0.038 inside the surface and carry that as
/// their tolerance. Every other edge lies on the cylinder.
fn chorded_face(model: &mut Model) -> Shape {
    let pi = core::f64::consts::PI;
    let surface =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, RADIUS, T).unwrap(), (-1.0, 4.0)).unwrap();
    let on = |u: f64, v: f64| surface.point_at(u, v, T).unwrap();
    let cylinder = model.geometry_mut().add_surface(surface.into());
    // The chart corners, counter-clockwise: up the cut, along the bottom,
    // up the right side, back along the top and down the left side.
    let (left, right) = (0.75 * pi, 1.25 * pi);
    let corners = [
        (left, 1.0),
        (0.875 * pi, 0.6173),
        (pi, 0.2929),
        (pi, 0.0),
        (right, 0.0),
        (right, 3.0),
        (left, 3.0),
    ];
    let vertices: Vec<Shape> = corners
        .iter()
        .map(|&(u, v)| make_vertex(model, on(u, v)).shape)
        .collect();
    let chart = |model: &mut Model, edge: &Shape, from: (f64, f64), to: (f64, f64)| {
        let (a, b) = (Point2::new(from.0, from.1), Point2::new(to.0, to.1));
        let line = Line2d::segment(a, b, T).unwrap();
        attach_pcurve(
            model,
            edge,
            line.into(),
            cylinder,
            Location::identity(),
            (0.0, a.distance(b)),
        )
        .unwrap();
    };
    let mut edges = Vec::new();
    for k in 0..corners.len() {
        let (from, to) = (corners[k], corners[(k + 1) % corners.len()]);
        let (a, b) = (on(from.0, from.1), on(to.0, to.1));
        let (vf, vt) = (&vertices[k], &vertices[(k + 1) % corners.len()]);
        let edge = if (from.1 - to.1).abs() < 1e-12 {
            // Along a parallel: the circle at that height.
            let frame =
                Frame::new(Point::new(0.0, 0.0, from.1), Direction::Z, Direction::X, T).unwrap();
            let circle = CircleCurve::new(Circle::new(frame, RADIUS, T).unwrap());
            let range = (from.0.min(to.0), from.0.max(to.0));
            let (lo, hi) = if from.0 < to.0 { (vf, vt) } else { (vt, vf) };
            let edge = make_edge_between(model, circle.into(), range, lo, hi, T)
                .unwrap()
                .shape;
            let (start, end) = if from.0 < to.0 {
                (from, to)
            } else {
                (to, from)
            };
            chart(model, &edge, start, end);
            if from.0 < to.0 {
                edge
            } else {
                edge.oriented(Orientation::Reversed)
            }
        } else {
            let line = LineCurve::segment(a, b, T).unwrap();
            let length = a.distance(b);
            let edge = make_edge_between(model, line.into(), (0.0, length), vf, vt, T)
                .unwrap()
                .shape;
            chart(model, &edge, from, to);
            // How far the straight edge stands from the surface along its
            // pcurve: nothing on a ruling, the chord's sag across the cut.
            let lifted =
                Line2d::segment(Point2::new(from.0, from.1), Point2::new(to.0, to.1), T).unwrap();
            let span = Point2::new(from.0, from.1).distance(Point2::new(to.0, to.1));
            let mut apart = 0.0_f64;
            for i in 0..=256 {
                let s = f64::from(i) / 256.0;
                let p = line.point_at(length * s, T).unwrap();
                let q = ogeom::geom::Curve2d::point_at(&lifted, span * s, T).unwrap();
                apart = apart.max(p.distance(on(q.x, q.y)));
            }
            if let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(|n| n.data_mut()) {
                data.tolerance = data.tolerance.widen_to(apart);
            }
            edge
        };
        edges.push(edge);
    }
    let wire = make_wire(model, &edges, T).unwrap().shape;
    make_face_on(model, cylinder, &[wire], T).unwrap().shape
}

/// How far the mesh stands off the cylinder: the farthest any triangle's
/// middle lies from radius 2.
fn farthest_off(mesh: &ogeom::topo::Triangulation) -> f64 {
    mesh.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            let (x, y) = ((a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0);
            (x.hypot(y) - RADIUS).abs()
        })
        .fold(0.0, f64::max)
}

/// Asked for a chord of a ten-thousandth, the face's mesh stands off the
/// cylinder by most of the chords' sag along the cut: the rim follows the
/// chords, which the face's neighbours share, and the triangles fanned off
/// it span the turn the chords cut across. A mesh that says it met its
/// chord lies within it of the surface at every triangle's middle.
#[test]
fn a_mesh_off_its_surface_by_more_than_its_chord_says_it_missed() {
    let mut model = Model::new();
    let face = chorded_face(&mut model);
    let chord = 1e-4;
    let mesh =
        ogeom::mesh::triangulate(&model, &face, Deflection::with_chord(chord).unwrap(), T).unwrap();
    let off = farthest_off(&mesh);
    assert!(
        !mesh.deflection_met || off <= chord + T.confusion(),
        "the mesh says it met a chord of {chord} and stands {off:.3e} off the cylinder"
    );
}
