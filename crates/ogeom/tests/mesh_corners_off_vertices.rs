//! A unit cube whose twelve edges are lines between its true corners, and
//! whose eight vertices stand a hundredth in from those corners along the
//! cube's diagonals, within the tolerance the vertices record: a file's
//! slop at its corners. Every face is bounded by straight edges alone, so
//! an edge drawn from vertex to vertex takes the whole face off its plane.
//! Meshed finer than the vertices stand off, the cube keeps its volume.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{make_edge_between, make_face_with_pcurves, make_shell, make_solid, make_vertex};
use ogeom::core::FastMap;
use ogeom::core::{Tolerance, Tolerances};
use ogeom::geom::{LineCurve, PlaneSurface};
use ogeom::math::{Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::{Deflection, triangulate};
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

/// How far each vertex stands in from its corner.
const OFF: f64 = 1e-2;

/// The corner with index `i + 2j + 4k` is at `(i, j, k)`.
fn corner(index: usize) -> Point {
    let bit = |b: usize| f64::from(u8::from(index & b != 0));
    Point::new(bit(1), bit(2), bit(4))
}

fn cube(model: &mut Model) -> Shape {
    let vertices: Vec<Shape> = (0..8)
        .map(|i| {
            let c = corner(i);
            let inward = Vector::new(0.5, 0.5, 0.5) - c.to_vector();
            let at = c + inward * (OFF / inward.magnitude());
            let vertex = make_vertex(model, at).shape;
            model
                .widen(&vertex, Tolerance::new(2.0 * OFF).unwrap())
                .unwrap();
            vertex
        })
        .collect();
    let mut edges: FastMap<(usize, usize), Shape> = FastMap::default();
    for a in 0..8 {
        for bit in [1, 2, 4] {
            if a & bit == 0 {
                let b = a | bit;
                let line = LineCurve::segment(corner(a), corner(b), T).unwrap();
                let edge = make_edge_between(
                    model,
                    line.into(),
                    (0.0, 1.0),
                    &vertices[a],
                    &vertices[b],
                    T,
                )
                .unwrap()
                .shape;
                edges.insert((a, b), edge);
            }
        }
    }
    // Each face's corners counter-clockwise seen from outside, its outward
    // normal and a direction in it.
    let sides: [([usize; 4], Vector, Vector); 6] = [
        ([0, 2, 3, 1], -Vector::Z, Vector::X),
        ([4, 5, 7, 6], Vector::Z, Vector::X),
        ([0, 1, 5, 4], -Vector::Y, Vector::X),
        ([2, 6, 7, 3], Vector::Y, Vector::X),
        ([0, 4, 6, 2], -Vector::X, Vector::Y),
        ([1, 3, 7, 5], Vector::X, Vector::Y),
    ];
    let faces: Vec<Shape> = sides
        .iter()
        .map(|(ring, normal, x)| {
            let edges: Vec<Shape> = (0..4)
                .map(|k| {
                    let (p, q) = (ring[k], ring[(k + 1) % 4]);
                    let edge = &edges[&(p.min(q), p.max(q))];
                    if p < q { edge.clone() } else { edge.reversed() }
                })
                .collect();
            let frame = Frame::new(
                corner(ring[0]),
                Direction::new(*normal, T).unwrap(),
                Direction::new(*x, T).unwrap(),
                T,
            )
            .unwrap();
            let plane = PlaneSurface::new(Plane::new(frame));
            make_face_with_pcurves(model, plane.into(), &[edges], T)
                .unwrap()
                .shape
        })
        .collect();
    let shell = make_shell(model, &faces).unwrap().shape;
    make_solid(model, &[shell]).unwrap().shape
}

/// Meshed at chords under how far the vertices stand off, each face stays
/// on its plane but for the step from a point on each edge, as far in
/// from the corner as the vertex stands off it, to the vertex. What the
/// mesh leaves out or takes in lies within a cube of that side at each
/// corner, and the rest is the unit cube.
#[test]
fn a_cube_whose_vertices_stand_off_its_corners_meshes_to_its_volume() {
    let mut model = Model::new();
    let solid = cube(&mut model);
    for chord in [1e-3, 1e-4] {
        let mesh = triangulate(&model, &solid, Deflection::with_chord(chord).unwrap(), T).unwrap();
        assert!(mesh.is_closed(), "the mesh at chord {chord} is open");
        let corners = 8.0 * OFF.powi(3);
        assert!(
            (mesh.volume() - 1.0).abs() <= corners,
            "volume {} against 1 at chord {chord}",
            mesh.volume()
        );
    }
}
