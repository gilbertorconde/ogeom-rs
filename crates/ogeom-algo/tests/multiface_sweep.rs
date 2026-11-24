//! Prisms and revolutions of a profile of several faces sharing edges,
//! against closed forms: each piece's volume, the wall a shared edge sweeps
//! held once by the solids on both sides, each profile edge's wall placed
//! where the edge's motion puts it, and profiles whose faces overlap
//! refused.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use ogeom_algo::sweep::roles::{SWEEP_SIDE, SWEEP_TOP};
use ogeom_algo::{
    check, make_compound, make_edge_between, make_face, make_prism, make_revolution, make_shell,
    make_wire, outer_solid, volume_properties,
};
use ogeom_core::{FastMap, Tolerances};
use ogeom_geom::{LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Axis, Direction, Frame, Plane, Point, Transform, Vector};
use ogeom_mesh::Deflection;
use ogeom_topo::{Filter, Model, Shape, ShapeType, explore, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A prism's volume over flat walls is exact up to rounding.
const EXACT: f64 = 1e-9;

/// A turned volume, relative to its closed form.
const CHORD: f64 = 1e-4;
const TURNED: f64 = 1e-9;

/// Points of one plane joined by shared edges: each pair of corners is one
/// edge however many faces walk it.
struct Sketch {
    frame: Frame,
    corners: Vec<Point>,
    vertices: Vec<Shape>,
    edges: FastMap<(usize, usize), Shape>,
}

impl Sketch {
    /// Corners `(u, v)` on `frame`'s xy plane.
    fn new(model: &mut Model, frame: Frame, corners: &[(f64, f64)]) -> Self {
        let corners: Vec<Point> = corners
            .iter()
            .map(|&(u, v)| frame.origin() + frame.x().vector() * u + frame.y().vector() * v)
            .collect();
        let vertices = corners
            .iter()
            .map(|&p| ogeom_algo::make_vertex(model, p).shape)
            .collect();
        Self {
            frame,
            corners,
            vertices,
            edges: FastMap::default(),
        }
    }

    /// The edge from corner `a` to corner `b`, walked that way.
    fn edge(&mut self, model: &mut Model, a: usize, b: usize) -> Shape {
        if let Some(edge) = self.edges.get(&(a, b)) {
            return edge.clone();
        }
        if let Some(edge) = self.edges.get(&(b, a)) {
            return edge.reversed();
        }
        let (from, to) = (self.corners[a], self.corners[b]);
        let length = (to - from).magnitude();
        let line = LineCurve::over(
            Axis {
                location: from,
                direction: Direction::new(to - from, T).unwrap(),
            },
            0.0,
            length,
        )
        .unwrap();
        let edge = make_edge_between(
            model,
            line.into(),
            (0.0, length),
            &self.vertices[a],
            &self.vertices[b],
            T,
        )
        .unwrap()
        .shape;
        self.edges.insert((a, b), edge.clone());
        edge
    }

    /// The face inside the polygon through `ring`'s corners.
    fn face(&mut self, model: &mut Model, ring: &[usize]) -> Shape {
        let edges: Vec<Shape> = (0..ring.len())
            .map(|i| self.edge(model, ring[i], ring[(i + 1) % ring.len()]))
            .collect();
        let wire = make_wire(model, &edges, T).unwrap().shape;
        let surface: SurfaceGeometry =
            PlaneSurface::over(Plane::new(self.frame), (-50.0, 50.0), (-50.0, 50.0))
                .unwrap()
                .into();
        make_face(model, surface, &[wire], T).unwrap().shape
    }

    fn shared(&self, a: usize, b: usize) -> Shape {
        self.edges
            .get(&(a, b))
            .or_else(|| self.edges.get(&(b, a)))
            .unwrap()
            .clone()
    }
}

fn point_of(model: &Model, vertex: &Shape) -> Point {
    let data = model.node(vertex).unwrap().data().as_vertex().unwrap();
    vertex.transform(model.datums()).unwrap().apply(data.point)
}

fn corners_of(model: &Model, shape: &Shape) -> Vec<Point> {
    explore_unique(model, shape, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| point_of(model, v))
        .collect()
}

/// Whether `points` and `expected` are the same set within `bound`.
fn same_points(points: &[Point], expected: &[Point], bound: f64) -> f64 {
    assert_eq!(points.len(), expected.len(), "{points:?} vs {expected:?}");
    let mut worst: f64 = 0.0;
    for p in expected {
        let nearest = points
            .iter()
            .map(|q| (*q - *p).magnitude())
            .fold(f64::INFINITY, f64::min);
        worst = worst.max(nearest);
    }
    assert!(worst <= bound, "{points:?} vs {expected:?}: {worst}");
    worst
}

fn solids_of(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore(model, shape, Filter::OfType(ShapeType::Solid)).unwrap()
}

fn faces_of(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore(model, shape, Filter::OfType(ShapeType::Face)).unwrap()
}

fn volume(model: &Model, solid: &Shape) -> f64 {
    let deflection = Deflection {
        chord: CHORD,
        ..Deflection::default()
    };
    volume_properties(model, solid, deflection, T).unwrap().mass
}

fn assert_valid(model: &Model, shape: &Shape) {
    let diagnosis = check(model, shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
}

/// The one wall `edge` generated, its role and source checked.
fn wall_of(model: &Model, history: &ogeom_algo::History, edge: &Shape) -> Shape {
    let walls: Vec<&Shape> = history
        .generated(edge)
        .iter()
        .filter(|s| model.kind_of(s).unwrap() == ShapeType::Face)
        .collect();
    assert_eq!(walls.len(), 1, "an edge generates one wall");
    let provenance = model.provenance_of(walls[0]).unwrap();
    assert_eq!(provenance.role(), Some(SWEEP_SIDE));
    walls[0].clone()
}

/// `wall` is one face node bounding both `a` and `b`, once each and from
/// opposite sides.
fn shared_by(model: &Model, wall: &Shape, a: &Shape, b: &Shape) {
    let find = |solid: &Shape| -> Vec<Shape> {
        faces_of(model, solid)
            .into_iter()
            .filter(|f| f.node() == wall.node())
            .collect()
    };
    let (in_a, in_b) = (find(a), find(b));
    assert_eq!((in_a.len(), in_b.len()), (1, 1));
    assert_ne!(in_a[0].orientation(), in_b[0].orientation());
}

/// Every face generates its own top cap and its solid.
fn caps_traced(model: &Model, history: &ogeom_algo::History, face: &Shape, solid: &Shape) {
    let made = history.generated(face);
    assert!(made.iter().any(|s| s.is_same(solid)));
    let tops: Vec<&Shape> = made
        .iter()
        .filter(|s| {
            model.kind_of(s).unwrap() == ShapeType::Face
                && model.provenance_of(s).and_then(|p| p.role()) == Some(SWEEP_TOP)
        })
        .collect();
    assert_eq!(tops.len(), 1);
    assert!(faces_of(model, solid).iter().any(|f| f.is_same(tops[0])));
}

/// Two unit squares sharing the edge x = 1, on the xy plane.
fn two_squares(model: &mut Model, frame: Frame, at: f64) -> (Sketch, Shape, Shape) {
    let mut sketch = Sketch::new(
        model,
        frame,
        &[
            (at, 0.0),
            (at + 1.0, 0.0),
            (at + 2.0, 0.0),
            (at, 1.0),
            (at + 1.0, 1.0),
            (at + 2.0, 1.0),
        ],
    );
    let left = sketch.face(model, &[0, 1, 4, 3]);
    let right = sketch.face(model, &[1, 2, 5, 4]);
    (sketch, left, right)
}

#[test]
fn two_squares_prism_into_two_solids_on_one_shared_wall() {
    let mut model = Model::new();
    let (sketch, left, right) = two_squares(&mut model, Frame::WORLD, 0.0);
    let profile = make_shell(&mut model, &[left.clone(), right.clone()])
        .unwrap()
        .shape;
    let travel = Vector::new(0.0, 0.0, 5.0);
    let built = make_prism(&mut model, &profile, travel, T).unwrap();

    assert_eq!(model.kind_of(&built.shape).unwrap(), ShapeType::CompSolid);
    assert_valid(&model, &built.shape);
    let solids = solids_of(&model, &built.shape);
    assert_eq!(solids.len(), 2);
    for solid in &solids {
        assert_valid(&model, solid);
        let v = volume(&model, solid);
        assert!((v - 5.0).abs() <= EXACT, "{v}");
    }
    // 6 + 6 faces, one of them held by both.
    assert_eq!(
        explore_unique(&model, &built.shape, ShapeType::Face)
            .unwrap()
            .len(),
        11
    );

    // Each profile edge sweeps its wall, the shared one included, over the
    // edge's ends and their translates.
    for ((a, b), edge) in &sketch.edges {
        let wall = wall_of(&model, &built.history, edge);
        let (p, q) = (sketch.corners[*a], sketch.corners[*b]);
        same_points(
            &corners_of(&model, &wall),
            &[p, q, p + travel, q + travel],
            EXACT,
        );
    }
    let shared = wall_of(&model, &built.history, &sketch.shared(1, 4));
    let (solid_left, solid_right) = (
        built
            .history
            .generated(&left)
            .iter()
            .find(|s| model.kind_of(s).unwrap() == ShapeType::Solid)
            .unwrap()
            .clone(),
        built
            .history
            .generated(&right)
            .iter()
            .find(|s| model.kind_of(s).unwrap() == ShapeType::Solid)
            .unwrap()
            .clone(),
    );
    shared_by(&model, &shared, &solid_left, &solid_right);
    caps_traced(&model, &built.history, &left, &solid_left);
    caps_traced(&model, &built.history, &right, &solid_right);
    assert!(
        built
            .history
            .generated(&profile)
            .iter()
            .any(|s| s.is_same(&built.shape))
    );

    // The two pieces as one: the outer walls and both caps of each, the
    // shared wall gone.
    let whole = outer_solid(&mut model, &built.shape).unwrap();
    assert_valid(&model, &whole.shape);
    let faces = faces_of(&model, &whole.shape);
    assert_eq!(faces.len(), 10);
    assert!(faces.iter().all(|f| f.node() != shared.node()));
    let v = volume(&model, &whole.shape);
    assert!((v - 10.0).abs() <= EXACT, "{v}");
}

#[test]
fn two_squares_revolve_into_two_solids_on_one_shared_wall() {
    // The squares on the xz plane, x from 1 to 3, about the z axis: the
    // shared edge x = 2 turns into a cylinder of radius 2.
    let frame = Frame::new(
        Point::ORIGIN,
        Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(),
        Direction::X,
        T,
    )
    .unwrap();
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Z,
    };
    for angle in [FRAC_PI_2, 2.0, TAU] {
        let mut model = Model::new();
        let (sketch, left, right) = two_squares(&mut model, frame, 1.0);
        let profile = make_shell(&mut model, &[left.clone(), right.clone()])
            .unwrap()
            .shape;
        let built = make_revolution(&mut model, &profile, axis, angle, T).unwrap();
        assert_eq!(model.kind_of(&built.shape).unwrap(), ShapeType::CompSolid);
        assert_valid(&model, &built.shape);

        let piece = |face: &Shape| {
            built
                .history
                .generated(face)
                .iter()
                .find(|s| model.kind_of(s).unwrap() == ShapeType::Solid)
                .unwrap()
                .clone()
        };
        let (inner, outer) = (piece(&left), piece(&right));
        // A turn of `angle` sweeps (angle / 2)(r1^2 - r0^2) h.
        for (solid, (r0, r1)) in [(&inner, (1.0_f64, 2.0_f64)), (&outer, (2.0, 3.0))] {
            assert_valid(&model, solid);
            let expected = angle / 2.0 * (r1 * r1 - r0 * r0);
            let v = volume(&model, solid);
            let off = (v - expected).abs() / expected;
            println!("revolve {angle:.4}: r {r0}..{r1}: {v} vs {expected}: {off:.2e}");
            assert!(off <= TURNED, "{v} vs {expected}");
        }
        let shared = wall_of(&model, &built.history, &sketch.shared(1, 4));
        shared_by(&model, &shared, &inner, &outer);
        if angle < TAU {
            caps_traced(&model, &built.history, &left, &inner);
            caps_traced(&model, &built.history, &right, &outer);
        }
        // Each edge's wall spans the edge's ends and their turned images.
        let turned = Transform::rotation(axis, angle);
        for ((a, b), edge) in &sketch.edges {
            let wall = wall_of(&model, &built.history, edge);
            let (p, q) = (sketch.corners[*a], sketch.corners[*b]);
            let expected: Vec<Point> = if angle < TAU {
                vec![p, q, turned.apply(p), turned.apply(q)]
            } else {
                vec![p, q]
            };
            same_points(&corners_of(&model, &wall), &expected, 1e-9);
        }
        let whole = outer_solid(&mut model, &built.shape).unwrap();
        assert_valid(&model, &whole.shape);
        let expected = angle / 2.0 * 8.0;
        let v = volume(&model, &whole.shape);
        assert!(
            (v - expected).abs() / expected <= TURNED,
            "{v} vs {expected}"
        );
        assert!(angle < TAU || (expected - 8.0 * PI).abs() < 1e-12);
    }
}

#[test]
fn three_faces_on_a_t_junction_share_three_walls() {
    // A 2 x 1 base whose top edge is split at (1, 1), and two unit squares
    // on it meeting along x = 1: the vertex (1, 1) is on all three faces.
    let mut model = Model::new();
    let mut sketch = Sketch::new(
        &mut model,
        Frame::WORLD,
        &[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (0.0, 2.0),
            (1.0, 2.0),
            (2.0, 2.0),
        ],
    );
    let base = sketch.face(&mut model, &[0, 1, 2, 3, 4]);
    let left = sketch.face(&mut model, &[4, 3, 6, 5]);
    let right = sketch.face(&mut model, &[3, 2, 7, 6]);
    let profile = make_compound(&mut model, &[base.clone(), left.clone(), right.clone()])
        .unwrap()
        .shape;
    let travel = Vector::new(0.0, 0.0, 3.0);
    let built = make_prism(&mut model, &profile, travel, T).unwrap();
    assert_eq!(model.kind_of(&built.shape).unwrap(), ShapeType::CompSolid);
    assert_valid(&model, &built.shape);

    let piece = |face: &Shape| {
        built
            .history
            .generated(face)
            .iter()
            .find(|s| model.kind_of(s).unwrap() == ShapeType::Solid)
            .unwrap()
            .clone()
    };
    let pieces = [
        (piece(&base), 6.0),
        (piece(&left), 3.0),
        (piece(&right), 3.0),
    ];
    for (solid, expected) in &pieces {
        assert_valid(&model, solid);
        let v = volume(&model, solid);
        assert!((v - expected).abs() <= EXACT, "{v} vs {expected}");
    }
    let [(b, _), (l, _), (r, _)] = &pieces;
    shared_by(
        &model,
        &wall_of(&model, &built.history, &sketch.shared(4, 3)),
        b,
        l,
    );
    shared_by(
        &model,
        &wall_of(&model, &built.history, &sketch.shared(3, 2)),
        b,
        r,
    );
    shared_by(
        &model,
        &wall_of(&model, &built.history, &sketch.shared(3, 6)),
        l,
        r,
    );
    // The rail from the junction is one edge bounding four walls: the base's
    // two halves of its top edge and both sides of x = 1.
    let junction = sketch.vertices[3].clone();
    let rail = explore_unique(&model, &built.shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (p, q) = ogeom_algo::edge_vertices(&model, e).unwrap().unwrap();
            let raised = corners_of(&model, e)
                .iter()
                .any(|p| (p.z - 3.0).abs() < 1e-9);
            (p.is_same(&junction) || q.is_same(&junction)) && raised
        });
    let rail = rail.expect("the junction sweeps a rail");
    let walls_on_rail = explore_unique(&model, &built.shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            explore(&model, f, Filter::OfType(ShapeType::Edge))
                .unwrap()
                .iter()
                .any(|e| e.is_same(&rail))
        })
        .count();
    assert_eq!(walls_on_rail, 3);
    for ((a, b), edge) in &sketch.edges {
        let wall = wall_of(&model, &built.history, edge);
        let (p, q) = (sketch.corners[*a], sketch.corners[*b]);
        same_points(
            &corners_of(&model, &wall),
            &[p, q, p + travel, q + travel],
            EXACT,
        );
    }
    let whole = outer_solid(&mut model, &built.shape).unwrap();
    assert_valid(&model, &whole.shape);
    let v = volume(&model, &whole.shape);
    assert!((v - 12.0).abs() <= EXACT, "{v}");
}

#[test]
fn faces_sharing_no_edge_sweep_apart() {
    let mut model = Model::new();
    let mut sketch = Sketch::new(
        &mut model,
        Frame::WORLD,
        &[
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (2.0, 0.0),
            (3.0, 0.0),
            (3.0, 1.0),
            (2.0, 1.0),
        ],
    );
    let a = sketch.face(&mut model, &[0, 1, 2, 3]);
    let b = sketch.face(&mut model, &[4, 5, 6, 7]);
    let profile = make_compound(&mut model, &[a, b]).unwrap().shape;
    let built = make_prism(&mut model, &profile, Vector::new(0.0, 0.0, 2.0), T).unwrap();
    assert_eq!(model.kind_of(&built.shape).unwrap(), ShapeType::Compound);
    let solids = solids_of(&model, &built.shape);
    assert_eq!(solids.len(), 2);
    for solid in &solids {
        assert_valid(&model, solid);
        assert!((volume(&model, solid) - 2.0).abs() <= EXACT);
    }
}

#[test]
fn overlapping_or_unjoined_profiles_are_refused() {
    let up = Vector::new(0.0, 0.0, 1.0);
    // Two squares overlapping by a quarter, their edges crossing.
    let mut model = Model::new();
    let mut sketch = Sketch::new(
        &mut model,
        Frame::WORLD,
        &[
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (0.5, 0.5),
            (1.5, 0.5),
            (1.5, 1.5),
            (0.5, 1.5),
        ],
    );
    let a = sketch.face(&mut model, &[0, 1, 2, 3]);
    let b = sketch.face(&mut model, &[4, 5, 6, 7]);
    let profile = make_compound(&mut model, &[a, b]).unwrap().shape;
    refused(make_prism(&mut model, &profile, up, T), "share no vertex");

    // A small square inside a large one.
    let mut model = Model::new();
    let mut sketch = Sketch::new(
        &mut model,
        Frame::WORLD,
        &[
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 4.0),
            (0.0, 4.0),
            (1.0, 1.0),
            (2.0, 1.0),
            (2.0, 2.0),
            (1.0, 2.0),
        ],
    );
    let a = sketch.face(&mut model, &[0, 1, 2, 3]);
    let b = sketch.face(&mut model, &[4, 5, 6, 7]);
    let profile = make_compound(&mut model, &[a, b]).unwrap().shape;
    refused(
        make_prism(&mut model, &profile, up, T),
        "holds part of another",
    );

    // A T junction whose base edge is not split where the stem ends.
    let mut model = Model::new();
    let mut sketch = Sketch::new(
        &mut model,
        Frame::WORLD,
        &[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (0.0, 1.0),
            (1.0, 1.0),
            (1.0, 2.0),
            (0.0, 2.0),
        ],
    );
    let base = sketch.face(&mut model, &[0, 1, 2, 3]);
    let stem = sketch.face(&mut model, &[3, 4, 5, 6]);
    let profile = make_compound(&mut model, &[base, stem]).unwrap().shape;
    refused(make_prism(&mut model, &profile, up, T), "running along");

    // Two squares sharing an edge, folded out of one plane.
    let mut model = Model::new();
    let (mut sketch, flat, _) = two_squares(&mut model, Frame::WORLD, 0.0);
    let raised = model_face_raised(&mut model, &mut sketch);
    let profile = make_compound(&mut model, &[flat, raised]).unwrap().shape;
    refused(
        make_prism(&mut model, &profile, Vector::new(1.0, 1.0, 1.0), T),
        "not in one plane",
    );

    // Squares on both sides of the axis they share an edge on.
    let mut model = Model::new();
    let frame = Frame::new(
        Point::ORIGIN,
        Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(),
        Direction::X,
        T,
    )
    .unwrap();
    let (_, left, right) = two_squares(&mut model, frame, -1.0);
    let profile = make_shell(&mut model, &[left, right]).unwrap().shape;
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Z,
    };
    refused(
        make_revolution(&mut model, &profile, axis, PI, T),
        "both sides of the axis",
    );
}

/// A unit square standing up from the sketch's shared edge x = 1, square
/// to the sketch's plane.
fn model_face_raised(model: &mut Model, sketch: &mut Sketch) -> Shape {
    let top = ogeom_algo::make_vertex(model, Point::new(1.0, 1.0, 1.0)).shape;
    let low = ogeom_algo::make_vertex(model, Point::new(1.0, 0.0, 1.0)).shape;
    sketch
        .corners
        .extend([Point::new(1.0, 1.0, 1.0), Point::new(1.0, 0.0, 1.0)]);
    sketch.vertices.extend([top, low]);
    let edges: Vec<Shape> = [(1, 4), (4, 6), (6, 7), (7, 1)]
        .into_iter()
        .map(|(a, b)| sketch.edge(model, a, b))
        .collect();
    let wire = make_wire(model, &edges, T).unwrap().shape;
    let frame = Frame::new(Point::new(1.0, 0.0, 0.0), Direction::X, Direction::Y, T).unwrap();
    let surface: SurfaceGeometry =
        PlaneSurface::over(Plane::new(frame), (-50.0, 50.0), (-50.0, 50.0))
            .unwrap()
            .into();
    make_face(model, surface, &[wire], T).unwrap().shape
}

/// `result` is an error saying `why`.
fn refused(result: ogeom_core::OgeomResult<ogeom_algo::Built>, why: &str) {
    match result {
        Ok(_) => panic!("built a profile that should be refused ({why})"),
        Err(e) => assert!(e.to_string().contains(why), "{e}"),
    }
}

#[test]
fn a_ring_and_the_disc_in_its_hole_share_the_circle() {
    // A 4 x 4 square with a hole of radius 1 at its centre, and the disc
    // filling the hole on the same circle edge: the circle is the ring's
    // hole and the disc's outline, and sweeps one cylinder held by both.
    let mut model = Model::new();
    let mut sketch = Sketch::new(
        &mut model,
        Frame::WORLD,
        &[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
    );
    let edges: Vec<Shape> = (0..4)
        .map(|i| sketch.edge(&mut model, i, (i + 1) % 4))
        .collect();
    let outline = make_wire(&mut model, &edges, T).unwrap().shape;
    let centre = Frame::new(Point::new(2.0, 2.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let curve: ogeom_geom::Curve =
        ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(centre, 1.0, T).unwrap()).into();
    let domain = ogeom_geom::Curve3d::domain(&curve);
    let circle = ogeom_algo::make_edge(&mut model, curve, domain, T)
        .unwrap()
        .shape;
    let ring_wire = make_wire(&mut model, std::slice::from_ref(&circle), T)
        .unwrap()
        .shape;
    let surface = || -> SurfaceGeometry {
        PlaneSurface::over(Plane::new(Frame::WORLD), (-50.0, 50.0), (-50.0, 50.0))
            .unwrap()
            .into()
    };
    let ring = make_face(&mut model, surface(), &[outline, ring_wire.clone()], T)
        .unwrap()
        .shape;
    let disc = make_face(&mut model, surface(), &[ring_wire], T)
        .unwrap()
        .shape;
    let profile = make_shell(&mut model, &[ring.clone(), disc.clone()])
        .unwrap()
        .shape;
    let built = make_prism(&mut model, &profile, Vector::new(0.0, 0.0, 2.0), T).unwrap();
    assert_eq!(model.kind_of(&built.shape).unwrap(), ShapeType::CompSolid);
    assert_valid(&model, &built.shape);
    let piece = |face: &Shape| {
        built
            .history
            .generated(face)
            .iter()
            .find(|s| model.kind_of(s).unwrap() == ShapeType::Solid)
            .unwrap()
            .clone()
    };
    let (outer, inner) = (piece(&ring), piece(&disc));
    for (solid, expected) in [(&outer, (16.0 - PI) * 2.0), (&inner, PI * 2.0)] {
        assert_valid(&model, solid);
        let v = volume(&model, solid);
        assert!(
            (v - expected).abs() / expected <= TURNED,
            "{v} vs {expected}"
        );
    }
    shared_by(
        &model,
        &wall_of(&model, &built.history, &circle),
        &outer,
        &inner,
    );
    let whole = outer_solid(&mut model, &built.shape).unwrap();
    assert_valid(&model, &whole.shape);
    assert!((volume(&model, &whole.shape) - 32.0).abs() <= 1e-9);
}

#[test]
fn stacked_squares_revolve_on_a_flat_shared_wall() {
    // Two squares one above the other on the xz plane, x from 1 to 2: the
    // shared edge z = 1 is square to the axis and turns into a flat ring.
    let frame = Frame::new(
        Point::ORIGIN,
        Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(),
        Direction::X,
        T,
    )
    .unwrap();
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Z,
    };
    for angle in [1.0, TAU] {
        let mut model = Model::new();
        let mut sketch = Sketch::new(
            &mut model,
            frame,
            &[
                (1.0, 0.0),
                (2.0, 0.0),
                (2.0, 1.0),
                (1.0, 1.0),
                (2.0, 2.0),
                (1.0, 2.0),
            ],
        );
        let low = sketch.face(&mut model, &[0, 1, 2, 3]);
        let high = sketch.face(&mut model, &[3, 2, 4, 5]);
        let profile = make_shell(&mut model, &[low.clone(), high.clone()])
            .unwrap()
            .shape;
        let built = make_revolution(&mut model, &profile, axis, angle, T).unwrap();
        assert_valid(&model, &built.shape);
        let piece = |face: &Shape| {
            built
                .history
                .generated(face)
                .iter()
                .find(|s| model.kind_of(s).unwrap() == ShapeType::Solid)
                .unwrap()
                .clone()
        };
        let (a, b) = (piece(&low), piece(&high));
        let expected = angle / 2.0 * 3.0;
        for solid in [&a, &b] {
            assert_valid(&model, solid);
            let v = volume(&model, solid);
            assert!(
                (v - expected).abs() / expected <= TURNED,
                "{v} vs {expected}"
            );
        }
        let wall = wall_of(&model, &built.history, &sketch.shared(2, 3));
        let data = model.node(&wall).unwrap().data().as_face().unwrap();
        assert_eq!(
            ogeom_geom::Surface::kind(model.geometry().surface(data.surface).unwrap()),
            ogeom_geom::SurfaceKind::Plane
        );
        shared_by(&model, &wall, &a, &b);
    }
}
