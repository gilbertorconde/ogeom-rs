//! The medial axis of general planar faces: holes, reflex corners, arcs
//! and free curves, each branch checked against the boundary it keeps its
//! distance from.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::Curve3d as _;
use ogeom_math::{Circle, Direction, Frame, Point};
use ogeom_topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

/// A boundary piece the test measures distance to in closed form (or, for
/// a free curve, by a dense sampling).
#[derive(Clone)]
enum Piece {
    Line(Point, Point),
    /// Centre, radius, start angle, end angle (counter-clockwise).
    Arc(Point, f64, f64, f64),
    Sampled(Vec<Point>),
}

impl Piece {
    fn distance(&self, p: Point) -> f64 {
        match self {
            Self::Line(a, b) => {
                let d = *b - *a;
                let s = ((p - *a).dot(d) / d.dot(d)).clamp(0.0, 1.0);
                p.distance(*a + d * s)
            }
            Self::Arc(c, r, a0, a1) => {
                let theta = (p.y - c.y).atan2(p.x - c.x);
                let rel = (theta - a0).rem_euclid(core::f64::consts::TAU);
                if rel <= a1 - a0 {
                    ((p - *c).magnitude() - r).abs()
                } else {
                    let at = |a: f64| Point::new(c.x + r * a.cos(), c.y + r * a.sin(), 0.0);
                    p.distance(at(*a0)).min(p.distance(at(*a1)))
                }
            }
            Self::Sampled(points) => points
                .windows(2)
                .map(|w| {
                    let d = w[1] - w[0];
                    let s = ((p - w[0]).dot(d) / d.dot(d)).clamp(0.0, 1.0);
                    p.distance(w[0] + d * s)
                })
                .fold(f64::INFINITY, f64::min),
        }
    }
}

/// A face on the xy plane from loops of pieces, each loop's pieces head
/// to tail. The first loop is the outer.
fn face_of(model: &mut Model, loops: &[Vec<Piece>]) -> Shape {
    use ogeom_geom::{PlaneSurface, SurfaceGeometry};
    let mut wires: Vec<Vec<Shape>> = Vec::new();
    for pieces in loops {
        let start = |p: &Piece| match p {
            Piece::Line(a, _) => *a,
            Piece::Arc(c, r, a0, _) => Point::new(c.x + r * a0.cos(), c.y + r * a0.sin(), 0.0),
            Piece::Sampled(points) => points[0],
        };
        let vertices: Vec<Shape> = pieces
            .iter()
            .map(|p| ogeom_algo::make_vertex(model, start(p)).shape)
            .collect();
        let mut edges = Vec::new();
        for (i, piece) in pieces.iter().enumerate() {
            let (from, to) = (&vertices[i], &vertices[(i + 1) % pieces.len()]);
            let (curve, range): (ogeom_geom::Curve, (f64, f64)) = match piece {
                Piece::Line(a, b) => {
                    let c: ogeom_geom::Curve =
                        ogeom_geom::LineCurve::segment(*a, *b, T).unwrap().into();
                    let d = c.domain();
                    (c, d)
                }
                Piece::Arc(c, r, a0, a1) => (
                    ogeom_geom::CircleCurve::new(
                        Circle::new(
                            Frame::new(*c, Direction::Z, Direction::X, T).unwrap(),
                            *r,
                            T,
                        )
                        .unwrap(),
                    )
                    .into(),
                    (*a0, *a1),
                ),
                Piece::Sampled(_) => unreachable!("free curves are built by their tests"),
            };
            edges.push(
                ogeom_algo::make_edge_between(model, curve, range, from, to, T)
                    .unwrap()
                    .shape,
            );
        }
        wires.push(edges);
    }
    ogeom_algo::make_face_with_pcurves(
        model,
        SurfaceGeometry::Plane(PlaneSurface::new(ogeom_math::Plane::new(Frame::WORLD))),
        &wires,
        T,
    )
    .unwrap()
    .shape
}

fn polygon(corners: &[(f64, f64)]) -> Vec<Piece> {
    (0..corners.len())
        .map(|i| {
            let (a, b) = (corners[i], corners[(i + 1) % corners.len()]);
            Piece::Line(Point::new(a.0, a.1, 0.0), Point::new(b.0, b.1, 0.0))
        })
        .collect()
}

/// The defining property: along every branch the clearance is the
/// distance to the nearest boundary piece, so no piece comes nearer than
/// the branch's own two. Returns the branches' total length.
fn check_axis(
    model: &Model,
    face: &Shape,
    boundary: &[Piece],
    within: f64,
) -> ogeom_algo::MedialGraph {
    let graph = ogeom_algo::medial_graph(model, face, 1e-4, T).unwrap();
    assert!(!graph.branches.is_empty(), "the axis has branches");
    assert!(graph.deviation <= 1e-4, "deviation {}", graph.deviation);
    let nearest = |p: Point| {
        boundary
            .iter()
            .map(|b| b.distance(p))
            .fold(f64::INFINITY, f64::min)
    };
    for (b, branch) in graph.branches.iter().enumerate() {
        for k in 0..=20 {
            let t = branch.range.0 + (branch.range.1 - branch.range.0) * f64::from(k) / 20.0;
            let p = branch.curve.point_at(t, T).unwrap();
            let clearance = graph.clearance_at(b, t, T).unwrap();
            let truth = nearest(p);
            assert!(
                (clearance - truth).abs() <= within,
                "branch {b} at {p:?}: clearance {clearance}, nearest boundary {truth}"
            );
        }
    }
    for v in &graph.vertices {
        let truth = nearest(v.point);
        assert!(
            (v.clearance - truth).abs() <= within,
            "vertex {:?}: clearance {}, nearest boundary {truth}",
            v.point,
            v.clearance
        );
    }
    graph
}

fn length(graph: &ogeom_algo::MedialGraph) -> f64 {
    graph
        .branches
        .iter()
        .map(|b| {
            let mut sum = 0.0;
            let mut last = b.curve.point_at(b.range.0, T).unwrap();
            for k in 1..=2000 {
                let t = b.range.0 + (b.range.1 - b.range.0) * f64::from(k) / 2000.0;
                let p = b.curve.point_at(t, T).unwrap();
                sum += last.distance(p);
                last = p;
            }
            sum
        })
        .sum()
}

#[test]
fn a_rectangles_graph_is_its_roof_line() {
    let mut model = Model::new();
    let pieces = polygon(&[(0.0, 0.0), (20.0, 0.0), (20.0, 10.0), (0.0, 10.0)]);
    let face = face_of(&mut model, std::slice::from_ref(&pieces));
    let graph = check_axis(&model, &face, &pieces, 1e-9);
    let want = 4.0 * 5.0 * core::f64::consts::SQRT_2 + 10.0;
    let got = length(&graph);
    assert!((got - want).abs() < 1e-6, "{got} against {want}");
}

#[test]
fn an_l_shapes_reflex_corner_bends_its_axis_along_a_parabola() {
    let mut model = Model::new();
    let pieces = polygon(&[
        (0.0, 0.0),
        (20.0, 0.0),
        (20.0, 10.0),
        (10.0, 10.0),
        (10.0, 20.0),
        (0.0, 20.0),
    ]);
    let face = face_of(&mut model, std::slice::from_ref(&pieces));
    let graph = check_axis(&model, &face, &pieces, 1e-9);
    assert!(
        graph
            .branches
            .iter()
            .any(|b| matches!(b.curve, ogeom_geom::Curve::Parabola(_))),
        "the reflex corner bisects its far walls along parabolas"
    );
}

#[test]
fn a_square_ring_holds_its_axis_round_the_hole() {
    let mut model = Model::new();
    let outer = polygon(&[(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0)]);
    let hole = polygon(&[(6.0, 6.0), (6.0, 14.0), (14.0, 14.0), (14.0, 6.0)]);
    let face = face_of(&mut model, &[outer.clone(), hole.clone()]);
    let mut all = outer;
    all.extend(hole);
    let graph = check_axis(&model, &face, &all, 1e-9);
    // The ring's midline at clearance 3 runs round the hole: some branch
    // stands at clearance 3 on each side.
    for probe in [
        Point::new(10.0, 3.0, 0.0),
        Point::new(17.0, 10.0, 0.0),
        Point::new(10.0, 17.0, 0.0),
        Point::new(3.0, 10.0, 0.0),
    ] {
        let near = graph.branches.iter().any(|b| {
            (0..=200).any(|k| {
                let t = b.range.0 + (b.range.1 - b.range.0) * f64::from(k) / 200.0;
                b.curve.point_at(t, T).unwrap().distance(probe) < 0.05
            })
        });
        assert!(near, "the axis passes {probe:?}");
    }
}

#[test]
fn a_stadiums_axis_joins_its_arc_centres() {
    let mut model = Model::new();
    let pi = core::f64::consts::PI;
    let pieces = vec![
        Piece::Line(Point::new(0.0, 0.0, 0.0), Point::new(20.0, 0.0, 0.0)),
        Piece::Arc(Point::new(20.0, 5.0, 0.0), 5.0, -pi / 2.0, pi / 2.0),
        Piece::Line(Point::new(20.0, 10.0, 0.0), Point::new(0.0, 10.0, 0.0)),
        Piece::Arc(Point::new(0.0, 5.0, 0.0), 5.0, pi / 2.0, 3.0 * pi / 2.0),
    ];
    let face = face_of(&mut model, std::slice::from_ref(&pieces));
    let graph = check_axis(&model, &face, &pieces, 1e-9);
    let got = length(&graph);
    assert!(
        (got - 20.0).abs() < 1e-6,
        "the axis is the 20 between the centres: {got}"
    );
    for v in &graph.vertices {
        assert!(
            (v.clearance - 5.0).abs() < 1e-9,
            "clearance {}",
            v.clearance
        );
    }
}

#[test]
fn a_plate_with_a_round_hole_bisects_along_parabolas() {
    let mut model = Model::new();
    let pi = core::f64::consts::PI;
    let outer = polygon(&[(0.0, 0.0), (30.0, 0.0), (30.0, 20.0), (0.0, 20.0)]);
    // The hole's circle in two halves, stored counter-clockwise like the
    // outer loop: the reader winds every loop itself.
    let hole = vec![
        Piece::Arc(Point::new(15.0, 10.0, 0.0), 3.0, 0.0, pi),
        Piece::Arc(Point::new(15.0, 10.0, 0.0), 3.0, pi, 2.0 * pi),
    ];
    let face = face_of(&mut model, &[outer.clone(), hole.clone()]);
    let mut all = outer;
    all.extend(hole);
    let graph = check_axis(&model, &face, &all, 1e-9);
    assert!(
        graph
            .branches
            .iter()
            .any(|b| matches!(b.curve, ogeom_geom::Curve::Parabola(_))),
        "the hole and a wall bisect along a parabola"
    );
}

#[test]
fn a_half_ellipse_fits_its_branches_to_the_tolerance() {
    let mut model = Model::new();
    let ellipse = ogeom_math::Ellipse::new(Frame::WORLD, 10.0, 5.0, T).unwrap();
    let curve: ogeom_geom::Curve = ogeom_geom::EllipseCurve::new(ellipse).into();
    let (a, b) = (Point::new(10.0, 0.0, 0.0), Point::new(-10.0, 0.0, 0.0));
    let va = ogeom_algo::make_vertex(&mut model, a).shape;
    let vb = ogeom_algo::make_vertex(&mut model, b).shape;
    let arc = ogeom_algo::make_edge_between(
        &mut model,
        curve.clone(),
        (0.0, core::f64::consts::PI),
        &va,
        &vb,
        T,
    )
    .unwrap()
    .shape;
    let line: ogeom_geom::Curve = ogeom_geom::LineCurve::segment(b, a, T).unwrap().into();
    let d = line.domain();
    let base = ogeom_algo::make_edge_between(&mut model, line, d, &vb, &va, T)
        .unwrap()
        .shape;
    let face = ogeom_algo::make_face_with_pcurves(
        &mut model,
        ogeom_geom::SurfaceGeometry::Plane(ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(
            Frame::WORLD,
        ))),
        &[vec![arc, base]],
        T,
    )
    .unwrap()
    .shape;
    let samples: Vec<Point> = (0..=40_000)
        .map(|i| {
            curve
                .point_at(core::f64::consts::PI * f64::from(i) / 40_000.0, T)
                .unwrap()
        })
        .collect();
    let pieces = vec![Piece::Sampled(samples), Piece::Line(b, a)];
    // The ellipse is measured through a polyline of chords under a micron
    // long, which sit within a hundredth of a micron of it.
    check_axis(&model, &face, &pieces, 2e-4);
}

/// A washer's axis is one closed branch between its two circles: the
/// ellipse with the circles' centres for foci, no branch point on it, its
/// clearance least where the hole comes nearest the rim.
#[test]
fn a_washers_axis_is_one_closed_branch_between_its_circles() {
    let tau = core::f64::consts::TAU;
    for (hole_centre, hole_radius, least) in [((0.0, 0.0), 7.0, 1.5), ((2.0, 0.0), 5.0, 1.5)] {
        let mut model = Model::new();
        let outer = vec![Piece::Arc(Point::ORIGIN, 10.0, 0.0, tau)];
        let hole = vec![Piece::Arc(
            Point::new(hole_centre.0, hole_centre.1, 0.0),
            hole_radius,
            0.0,
            tau,
        )];
        let face = face_of(&mut model, &[outer.clone(), hole.clone()]);
        let mut all = outer;
        all.extend(hole);
        let graph = check_axis(&model, &face, &all, 1e-9);
        assert_eq!(graph.branches.len(), 1);
        let branch = &graph.branches[0];
        assert_eq!(
            branch.ends[0], branch.ends[1],
            "the branch closes on itself"
        );
        assert!(((branch.range.1 - branch.range.0) - tau).abs() < 1e-12);
        let min = (0..=200)
            .map(|i| branch.range.0 + (branch.range.1 - branch.range.0) * f64::from(i) / 200.0)
            .map(|t| graph.clearance_at(0, t, T).unwrap())
            .fold(f64::INFINITY, f64::min);
        assert!((min - least).abs() < 1e-3, "least clearance {min}");
    }
}
