//! A torus face cut by a plane near the tube's rim, where the cut runs
//! almost straight in space and bends in the torus's chart: its mesh keeps
//! to the face, on the face's side of the cutting plane.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    Spacing, attach_pcurve, interpolate, make_edge_between, make_face_on, make_vertex, make_wire,
};
use ogeom::core::Tolerances;
use ogeom::geom::{
    CircleCurve, Curve, Curve2d as _, Curve3d as _, Line2d, Surface as _, SurfaceGeometry,
    TorusSurface,
};
use ogeom::math::{Circle, Direction, Frame, Point, Point2, Torus, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Location, Model, Orientation, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

const MAJOR: f64 = 5.35;
const MINOR: f64 = 2.0;
/// The cutting plane `y = -CUT`, parallel to the torus's axis.
const CUT: f64 = 4.5;
/// The chord the face is drawn to: the cut is two spans at it.
const CHORD: f64 = 0.0666;

/// The quarter of the tube between its inner equator and its bottom, over
/// the turn about the axis between the two places the plane `y = -4.5`
/// meets the bottom, and on the side of the plane towards the axis. The cut
/// runs from the bottom at `x = -2.89` to the bottom at `x = 2.89` through
/// its nearest the axis at `x = 0`, and bows 0.005 off a straight line in
/// space either side of there; in the chart the straight chord of either
/// half lifts to the torus 0.2 beyond the plane.
fn cut_face(model: &mut Model) -> Shape {
    let tau = core::f64::consts::TAU;
    let pi = core::f64::consts::PI;
    let torus = TorusSurface::new(Torus::new(Frame::WORLD, MAJOR, MINOR, T).unwrap());
    let on = |u: f64, v: f64| torus.point_at(u, v, T).unwrap();
    let surface = model.geometry_mut().add_surface(torus.into());
    // The cut at `x`: the tube's bottom half where `y = -4.5`.
    let cut_at = |x: f64| {
        let rho = x.hypot(CUT);
        Point::new(x, -CUT, -(MINOR * MINOR - (rho - MAJOR).powi(2)).sqrt())
    };
    let x_low = -(MAJOR * MAJOR - CUT * CUT).sqrt();
    let (u_low, u_high) = ((-CUT).atan2(x_low) + tau, (-CUT).atan2(-x_low) + tau);
    // Counter-clockwise in the chart: along the inner equator, up the
    // meridian to one end of the cut, back along the cut and down the
    // other meridian.
    let corners = [
        (u_low, pi),
        (u_high, pi),
        (u_high, 1.5 * pi),
        (u_low, 1.5 * pi),
    ];
    let vertices: Vec<Shape> = corners
        .iter()
        .map(|&(u, v)| make_vertex(model, on(u, v)).shape)
        .collect();
    let line = |model: &mut Model, edge: &Shape, from: (f64, f64), to: (f64, f64), range| {
        let (a, b) = (Point2::new(from.0, from.1), Point2::new(to.0, to.1));
        let chart = Line2d::segment(a, b, T).unwrap();
        attach_pcurve(
            model,
            edge,
            chart.into(),
            surface,
            Location::identity(),
            range,
        )
        .unwrap();
    };
    // The tube's circle at longitude `u`, parameterized by the tube angle.
    let meridian = |u: f64| {
        let out = Vector::new(u.cos(), u.sin(), 0.0);
        let frame = Frame::new(
            Point::ORIGIN + out * MAJOR,
            Direction::new(out.cross(Vector::Z), T).unwrap(),
            Direction::new(out, T).unwrap(),
            T,
        )
        .unwrap();
        CircleCurve::new(Circle::new(frame, MINOR, T).unwrap())
    };
    let mut edges = Vec::new();
    // Along the inner equator, the circle of radius 3.35.
    let equator = CircleCurve::new(Circle::new(Frame::WORLD, MAJOR - MINOR, T).unwrap());
    let edge = make_edge_between(
        model,
        equator.into(),
        (u_low, u_high),
        &vertices[0],
        &vertices[1],
        T,
    )
    .unwrap()
    .shape;
    line(model, &edge, corners[0], corners[1], (0.0, u_high - u_low));
    edges.push(edge);
    // Up a meridian to the cut.
    let edge = make_edge_between(
        model,
        meridian(u_high).into(),
        (pi, 1.5 * pi),
        &vertices[1],
        &vertices[2],
        T,
    )
    .unwrap()
    .shape;
    line(model, &edge, corners[1], corners[2], (0.0, 0.5 * pi));
    edges.push(edge);
    // The cut, interpolated through points of it, and its pcurve
    // projected and checked to run between the corners.
    let points: Vec<Point> = (0..=64)
        .map(|i| cut_at(x_low * (1.0 - f64::from(i) / 32.0)))
        .collect();
    let spline: Curve = interpolate(&points, 3, Spacing::default(), T)
        .unwrap()
        .into();
    let range = spline.domain();
    let (pcurve, ..) = ogeom::algo::pcurve_fit::fit_projected_pcurve(
        &spline,
        range,
        &SurfaceGeometry::from(torus),
        T,
    )
    .unwrap();
    let ends = [
        pcurve.point_at(range.0, T).unwrap(),
        pcurve.point_at(range.1, T).unwrap(),
    ];
    for (end, corner) in ends.iter().zip([corners[3], corners[2]]) {
        assert!(
            end.distance(Point2::new(corner.0, corner.1)) < 1e-6,
            "the cut's pcurve ends at {end:?}, not at the corner {corner:?}"
        );
    }
    let edge = make_edge_between(model, spline, range, &vertices[3], &vertices[2], T)
        .unwrap()
        .shape;
    attach_pcurve(model, &edge, pcurve, surface, Location::identity(), range).unwrap();
    edges.push(edge.oriented(Orientation::Reversed));
    // Down the other meridian, from the bottom to the inner equator.
    let edge = make_edge_between(
        model,
        meridian(u_low).into(),
        (pi, 1.5 * pi),
        &vertices[0],
        &vertices[3],
        T,
    )
    .unwrap()
    .shape;
    line(model, &edge, corners[0], corners[3], (0.0, 0.5 * pi));
    edges.push(edge.oriented(Orientation::Reversed));
    let wire = make_wire(model, &edges, T).unwrap().shape;
    make_face_on(model, surface, &[wire], T).unwrap().shape
}

/// Drawn at a chord of 0.0666, the cut is two spans in space, and the
/// chart chord of each lifts 0.21 beyond the plane. The triangles along a
/// span may not split it, since the face across the cut shares it, and
/// were refined towards it, out to 0.215 beyond the plane. Every point of
/// a mesh keeping within its chord of the face stands within the chord of
/// the side of the plane the face lies on. The cut's points are split
/// where its chart chords sag, and are its own polyline's: stored, each
/// edge's polyline runs through the face's stored mesh, every point a
/// vertex of it and every step a side of a triangle.
#[test]
fn a_cut_bending_in_the_chart_keeps_the_mesh_on_the_face() {
    let mut model = Model::new();
    let face = cut_face(&mut model);
    ogeom::mesh::tessellate(&mut model, &face, Deflection::with_chord(CHORD).unwrap(), T).unwrap();
    let mesh = ogeom::mesh::triangulation_of(&model, &face).unwrap();
    let beyond = mesh
        .positions
        .iter()
        .map(|p| -CUT - p.y)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(
        beyond <= CHORD,
        "the mesh reaches {beyond:.3e} beyond the cutting plane, past its chord of {CHORD}"
    );
    let vertex = |p: Point| {
        mesh.positions
            .iter()
            .position(|q| q.distance(p) <= T.confusion())
            .unwrap_or_else(|| panic!("no vertex of the face's mesh at {p:?}"))
    };
    let sides: ogeom::core::FastSet<(usize, usize)> = mesh
        .triangles
        .iter()
        .flat_map(|t| {
            (0..3).map(move |k| {
                let (a, b) = (t[k] as usize, t[(k + 1) % 3] as usize);
                (a.min(b), a.max(b))
            })
        })
        .collect();
    for edge in explore_unique(&model, &face, ShapeType::Edge).unwrap() {
        let (points, _) = ogeom::mesh::polyline_of(&model, &edge).unwrap();
        for step in points.windows(2) {
            let (a, b) = (vertex(step[0]), vertex(step[1]));
            assert!(
                sides.contains(&(a.min(b), a.max(b))),
                "the polyline steps from {:?} to {:?} across the face's mesh",
                step[0],
                step[1]
            );
        }
    }
}
