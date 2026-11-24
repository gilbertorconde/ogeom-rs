//! A point's distance to a free-form face whose surface's nearest point to
//! it lies outside the face, while a nearer point of the face than its
//! edges lies inside, at another local minimum.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_algo::{PointDistance, make_edge_between, make_face, make_vertex, make_wire};
use ogeom_core::Tolerances;
use ogeom_geom::{BSplineSurface, Curve, Line2d, PlanarCurve, Surface as _, SurfaceGeometry};
use ogeom_math::{Axis2, ControlGrid, Direction2, KnotVector, Point, Point2};
use ogeom_topo::{Location, Model, NodeData, Shape};

/// An edge, its exact chart line on the sheet, and its range.
type Charted = (Shape, PlanarCurve, (f64, f64));

const T: Tolerances = Tolerances::millimetres();

/// A bicubic sheet over `x` in [0, 4] and `y` in [0, 2] with two domes:
/// a taller one over (1, 1) and a lower one over (3, 1).
fn domes() -> BSplineSurface {
    let (nu, nv) = (9, 5);
    let mut points = Vec::with_capacity(nu * nv);
    for i in 0..nu {
        for j in 0..nv {
            #[allow(clippy::cast_precision_loss)]
            let (x, y) = (i as f64 * 0.5, j as f64 * 0.5);
            let z = match (i, j) {
                (2, 2) => 8.0,
                (6, 2) => 6.0,
                _ => 0.0,
            };
            points.push(Point::new(x, y, z));
        }
    }
    let grid = ControlGrid::new(points, nu, nv).unwrap();
    BSplineSurface::new(
        KnotVector::clamped_uniform(3, nu).unwrap(),
        KnotVector::clamped_uniform(3, nv).unwrap(),
        &grid,
        T,
    )
    .unwrap()
}

/// The ring of iso-curves of `surface` round the chart rectangle
/// `(u0, u1) x (v0, v1)`, counter-clockwise in the chart, and each edge's
/// chart line.
fn ring(
    model: &mut Model,
    surface: &BSplineSurface,
    (u0, u1): (f64, f64),
    (v0, v1): (f64, f64),
) -> (Shape, Vec<Charted>) {
    let corners = [(u0, v0), (u1, v0), (u1, v1), (u0, v1)];
    let vertices: Vec<Shape> = corners
        .iter()
        .map(|&(u, v)| make_vertex(model, surface.point_at(u, v, T).unwrap()).shape)
        .collect();
    let mut edges = Vec::new();
    let mut charted = Vec::new();
    for k in 0..4 {
        let (ka, kb) = (k, (k + 1) % 4);
        let (a, b) = (corners[ka], corners[kb]);
        // Along an iso-v line the curve's parameter is u, else v; the edge
        // runs from the corner at its range's start.
        let (curve, from, to, line) = if (a.1 - b.1).abs() < 1e-12 {
            let axis = Axis2::new(Point2::new(0.0, a.1), Direction2::X);
            (surface.iso_v_curve(a.1, T).unwrap(), a.0, b.0, axis)
        } else {
            let axis = Axis2::new(Point2::new(a.0, 0.0), Direction2::Y);
            (surface.iso_u_curve(a.0, T).unwrap(), a.1, b.1, axis)
        };
        let (start, end, range) = if from < to {
            (ka, kb, (from, to))
        } else {
            (kb, ka, (to, from))
        };
        let edge = make_edge_between(
            model,
            Curve::BSpline(curve),
            range,
            &vertices[start],
            &vertices[end],
            T,
        )
        .unwrap()
        .shape;
        let chart: PlanarCurve = Line2d::over(line, range.0, range.1).unwrap().into();
        charted.push((edge.clone(), chart, range));
        edges.push(if start == ka { edge } else { edge.reversed() });
    }
    (make_wire(model, &edges, T).unwrap().shape, charted)
}

/// The sheet with a hole round the taller dome: the surface's nearest
/// point to a point above the middle is on that dome, outside the face,
/// and the face's own nearest point is on the lower dome, inside it,
/// nearer than any point of its edges. The distance read is that point's,
/// against the least distance over a fine grid of the face.
#[test]
fn a_face_s_nearest_point_inside_its_trim_beats_its_edges() {
    let surface = domes();
    let mut model = Model::new();
    let (outer, mut charted) = ring(&mut model, &surface, (0.0, 1.0), (0.0, 1.0));
    let (hole, inner) = ring(&mut model, &surface, (0.1, 0.4), (0.25, 0.75));
    charted.extend(inner);
    let face = make_face(
        &mut model,
        SurfaceGeometry::BSpline(surface.clone()),
        &[outer, hole],
        T,
    )
    .unwrap()
    .shape;
    let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data().clone()) else {
        unreachable!()
    };
    for (edge, chart, range) in charted {
        ogeom_algo::attach_pcurve(
            &mut model,
            &edge,
            chart,
            data.surface,
            Location::identity(),
            range,
        )
        .unwrap();
    }
    let p = Point::new(2.0, 1.0, 4.0);
    let in_hole = |u: f64, v: f64| (0.1..=0.4).contains(&u) && (0.25..=0.75).contains(&v);
    // The least distance over the whole sheet, over the face, and over the
    // face's edges, read on a grid fine enough to place each.
    let n = 800;
    let (mut sheet, mut on_face, mut on_edges) = (f64::INFINITY, f64::INFINITY, f64::INFINITY);
    for i in 0..=n {
        for j in 0..=n {
            let (u, v) = (f64::from(i) / f64::from(n), f64::from(j) / f64::from(n));
            let d = surface.point_at(u, v, T).unwrap().distance(p);
            sheet = sheet.min(d);
            if !in_hole(u, v) {
                on_face = on_face.min(d);
            }
            let on_hole_ring = (0.1..=0.4).contains(&u)
                && ((v - 0.25).abs() < 1e-12 || (v - 0.75).abs() < 1e-12)
                || (0.25..=0.75).contains(&v)
                    && ((u - 0.1).abs() < 1e-12 || (u - 0.4).abs() < 1e-12);
            let on_outer = i == 0 || j == 0 || i == n || j == n;
            if on_hole_ring || on_outer {
                on_edges = on_edges.min(d);
            }
        }
    }
    // The setup: the sheet's nearest point is in the hole, and the face's
    // is clear of its edges.
    assert!(
        sheet < on_face - 0.05,
        "the sheet's {sheet} against the face's {on_face}"
    );
    assert!(
        on_face < on_edges - 0.05,
        "the face's {on_face} against its edges' {on_edges}"
    );
    let (got, at) = PointDistance::of(&model, &face, T)
        .unwrap()
        .nearest(p, T)
        .unwrap();
    eprintln!("nearest {got} at {at:?}; the face's grid {on_face}, its edges' {on_edges}");
    assert!(
        got <= on_face + 1e-9 && got >= on_face - 1e-4,
        "read {got} against the face's {on_face}"
    );
}
