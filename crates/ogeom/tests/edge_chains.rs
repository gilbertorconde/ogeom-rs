//! Edges treated as one operation: chamfers mitred along a chain, and three
//! or more fillets closing their shared corner with the rolling ball's
//! patch.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{make_edge_between, make_face_with_pcurves, make_prism, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::{Curve, Curve3d, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Filter, Model, Shape, ShapeType, VertexData, explore, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A 20 × 20 × 10 box swept from a square on XY: its top vertices and
/// edges are its bottom ones placed by the sweep's travel, which is what
/// a caller's prism hands over.
fn prism_box(model: &mut Model) -> Shape {
    let corners = [(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0)];
    let points: Vec<Point> = corners
        .iter()
        .map(|(x, y)| Point::new(*x, *y, 0.0))
        .collect();
    let vertices: Vec<Shape> = points
        .iter()
        .map(|p| model.add_vertex(VertexData::new(*p)))
        .collect();
    let edges: Vec<Shape> = (0..4)
        .map(|i| {
            let j = (i + 1) % 4;
            let curve = LineCurve::segment(points[i], points[j], T).unwrap();
            let range = curve.domain();
            make_edge_between(
                model,
                Curve::Line(curve),
                range,
                &vertices[i],
                &vertices[j],
                T,
            )
            .unwrap()
            .shape
        })
        .collect();
    let plane = Plane::new(Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap());
    let surface = PlaneSurface::over(plane, (-1.0, 21.0), (-1.0, 21.0)).unwrap();
    let face = make_face_with_pcurves(model, SurfaceGeometry::Plane(surface), &[edges], T)
        .unwrap()
        .shape;
    make_prism(model, &face, Vector::new(0.0, 0.0, 10.0), T)
        .unwrap()
        .shape
}

/// The edge whose placed midpoint is nearest `at`.
fn edge_near(model: &Model, solid: &Shape, at: Point) -> Shape {
    let mid = |e: &Shape| -> Point {
        let data = model.node(e).unwrap().data().as_edge().unwrap().clone();
        let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            return Point::new(f64::INFINITY, 0.0, 0.0);
        };
        let local = model
            .geometry()
            .curve(*curve)
            .unwrap()
            .point_at(f64::midpoint(range.0, range.1), T)
            .unwrap();
        e.transform(model.datums()).unwrap().apply(local)
    };
    explore_unique(model, solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| mid(a).distance(at).total_cmp(&mid(b).distance(at)))
        .unwrap()
}

fn faces(model: &Model, shape: &Shape) -> usize {
    explore(model, shape, Filter::OfType(ShapeType::Face))
        .unwrap()
        .len()
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// The chamfer chain: the box's four top edges bevelled
/// by one millimetre as one operation mitre at every corner: ten faces,
/// the four prisms less their four corner overlaps.
#[test]
fn a_chamfer_chain_mitres_at_shared_vertices() {
    let mut model = Model::new();
    let solid = prism_box(&mut model);
    let top: Vec<Shape> = [
        Point::new(10.0, 0.0, 10.0),
        Point::new(20.0, 10.0, 10.0),
        Point::new(10.0, 20.0, 10.0),
        Point::new(0.0, 10.0, 10.0),
    ]
    .iter()
    .map(|at| edge_near(&model, &solid, *at))
    .collect();
    let bevelled = ogeom::fillet::chamfer_edges(&mut model, &solid, &top, 1.0, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &bevelled, T).unwrap().is_valid());
    assert_eq!(faces(&model, &bevelled), 10);
    let expected = 4000.0 - 4.0 * (0.5 * 20.0) + 4.0 / 3.0;
    let measured = volume(&model, &bevelled);
    assert!(
        (measured - expected).abs() < 1e-2,
        "volume {measured} against {expected}"
    );
}

/// Three bevels at a convex corner meet at one point, and every edge of
/// the box bevelled at once is the same inclusion and exclusion: the prisms,
/// less their pairwise overlaps at each corner, plus the triple overlap,
/// a quarter of the cube of the distance.
#[test]
fn chamfer_chains_meet_three_at_a_corner() {
    let mut model = Model::new();
    let solid = prism_box(&mut model);
    let corner: Vec<Shape> = [
        Point::new(10.0, 20.0, 10.0),
        Point::new(20.0, 10.0, 10.0),
        Point::new(20.0, 20.0, 5.0),
    ]
    .iter()
    .map(|at| edge_near(&model, &solid, *at))
    .collect();
    let bevelled = ogeom::fillet::chamfer_edges(&mut model, &solid, &corner, 1.0, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &bevelled, T).unwrap().is_valid());
    let expected = 4000.0 - (0.5 * (20.0 + 20.0 + 10.0) - 3.0 / 3.0 + 0.25);
    let measured = volume(&model, &bevelled);
    assert!(
        (measured - expected).abs() < 1e-2,
        "volume {measured} against {expected}"
    );

    let every = explore_unique(&model, &solid, ShapeType::Edge).unwrap();
    assert_eq!(every.len(), 12);
    let bevelled = ogeom::fillet::chamfer_edges(&mut model, &solid, &every, 1.0, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &bevelled, T).unwrap().is_valid());
    assert_eq!(faces(&model, &bevelled), 18, "six walls and twelve bevels");
    let expected = 4000.0 - (0.5 * (8.0 * 20.0 + 4.0 * 10.0) - 8.0 + 8.0 * 0.25);
    let measured = volume(&model, &bevelled);
    assert!(
        (measured - expected).abs() < 1e-2,
        "volume {measured} against {expected}"
    );
}

/// The fillet corner: three fillets at a box's top corner
/// asked together close it with the octant of the sphere a radius in from
/// all three faces: ten faces, and inside the corner cube nothing but that
/// sphere. The corner is the prism's far end, placed by the sweep, so the
/// corner tool must read where the vertex stands and not where its node
/// was built.
#[test]
fn three_fillets_at_a_corner_close_it_with_a_sphere() {
    let mut model = Model::new();
    let solid = prism_box(&mut model);
    let corner: Vec<Shape> = [
        Point::new(10.0, 20.0, 10.0),
        Point::new(20.0, 10.0, 10.0),
        Point::new(20.0, 20.0, 5.0),
    ]
    .iter()
    .map(|at| edge_near(&model, &solid, *at))
    .collect();
    let rounded = ogeom::fillet::fillet_edges(&mut model, &solid, &corner, 2.0, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &rounded, T).unwrap().is_valid());
    assert_eq!(faces(&model, &rounded), 10);

    let chord = 1e-2;
    let mesh = ogeom::mesh::triangulate(
        &model,
        &rounded,
        Deflection {
            chord,
            ..Deflection::default()
        },
        T,
    )
    .unwrap();
    let centre = Point::new(18.0, 18.0, 8.0);
    let inside =
        |p: Point| p.x > 18.0 && p.x < 20.0 && p.y > 18.0 && p.y < 20.0 && p.z > 8.0 && p.z < 10.0;
    let mut seen = 0;
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        let g = Point::new(
            (a.x + b.x + c.x) / 3.0,
            (a.y + b.y + c.y) / 3.0,
            (a.z + b.z + c.z) / 3.0,
        );
        if !inside(g) {
            continue;
        }
        seen += 1;
        for p in [a, b, c] {
            assert!(
                (p.distance(centre) - 2.0).abs() < 1e-6,
                "a vertex inside the corner off the sphere: {p:?}"
            );
        }
        assert!(
            (g.distance(centre) - 2.0).abs() < 2.0 * chord,
            "a triangle inside the corner off the sphere: {g:?}"
        );
    }
    assert!(seen > 0, "the corner is drawn");
}

/// At a vertex of more edges the corner still closes: a square pyramid's
/// four slopes filleted together round the apex with one sphere, and a
/// rectangular one's (no single ball touches its four slopes) with two
/// spheres and a cylinder between them.
#[test]
fn fillet_chains_close_an_apex() {
    for (half_y, spheres, cylinders) in [(10.0, 1, 4), (5.0, 2, 5)] {
        let mut model = Model::new();
        let base = [
            Point::new(-10.0, -half_y, 0.0),
            Point::new(10.0, -half_y, 0.0),
            Point::new(10.0, half_y, 0.0),
            Point::new(-10.0, half_y, 0.0),
        ];
        let apex = Point::new(0.0, 0.0, 15.0);
        let polygon = ogeom::algo::make_polygon(&mut model, &base, true, T)
            .unwrap()
            .shape;
        let tip = ogeom::algo::make_vertex(&mut model, apex).shape;
        let pyramid = ogeom::offset::make_loft(&mut model, &polygon, &tip, T)
            .unwrap()
            .shape;
        let slopes: Vec<Shape> = base
            .iter()
            .map(|b| edge_near(&model, &pyramid, *b + (apex - *b) * 0.5))
            .collect();
        let rounded = ogeom::fillet::fillet_edges(&mut model, &pyramid, &slopes, 1.5, T)
            .unwrap()
            .shape;
        assert!(ogeom::algo::check(&model, &rounded, T).unwrap().is_valid());
        let mut counts = (0, 0);
        for face in explore_unique(&model, &rounded, ShapeType::Face).unwrap() {
            let data = model.node(&face).unwrap().data().as_face().unwrap();
            match model.geometry().surface(data.surface) {
                Some(SurfaceGeometry::Sphere(_)) => counts.0 += 1,
                Some(SurfaceGeometry::Cylinder(_)) => counts.1 += 1,
                _ => {}
            }
        }
        assert_eq!(counts, (spheres, cylinders), "half width {half_y}");
    }
}

/// A 10 mm cube shaved by a drum of radius 6 standing on its centre: each
/// side face is flat between `5 - sqrt(11)` and `5 + sqrt(11)` and the drum
/// rounds the corners between.
fn shaved_cube(model: &mut Model) -> Shape {
    let block = ogeom::algo::make_box(model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let foot = Frame::new(Point::new(5.0, 5.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drum = ogeom::algo::make_cylinder(model, foot, 6.0, 12.0, T)
        .unwrap()
        .shape;
    ogeom::boolean::common(model, &block, &drum, T)
        .unwrap()
        .shape
}

/// The midpoint integral of `f` over `[a, b]` in `n` steps.
fn integral(a: f64, b: f64, n: u32, f: impl Fn(f64) -> f64) -> f64 {
    let h = (b - a) / f64::from(n);
    (0..n)
        .map(|i| f((f64::from(i) + 0.5).mul_add(h, a)))
        .sum::<f64>()
        * h
}

/// What rounding the shaved cube's front bottom edge and the upright edge
/// at its end takes off, as the union of the two blends each rounds alone.
/// The bottom edge's blend is a prism along x, an `r` square less a
/// quarter disc, capped flush where the flat ends at `x = 5 + sqrt(11)`.
/// The upright's is a prism along z the full height, its section the
/// sliver between the front plane, the drum and the ball touching both.
/// They overlap where that sliver, short of the cap, stands over the
/// bottom blend's section.
fn shaved_corner_removed(r: f64) -> f64 {
    let end = 5.0 + 11.0_f64.sqrt();
    let flat = 2.0 * 11.0_f64.sqrt();
    let bottom = flat * r * r * (1.0 - core::f64::consts::FRAC_PI_4);
    // The upright's ball: a radius off the front plane and in from the
    // drum, touching the drum where the line from its axis through the
    // ball's centre meets it.
    let centre_x = 5.0 + ((6.0 - r).powi(2) - (r - 5.0).powi(2)).sqrt();
    let top = (r - 5.0).mul_add(6.0 / (6.0 - r), 5.0);
    let ball = |y: f64| centre_x + r.mul_add(r, -(y - r).powi(2)).max(0.0).sqrt();
    let drum = |y: f64| 5.0 + (36.0 - (y - 5.0).powi(2)).sqrt();
    let sliver = integral(0.0, top, 200_000, |y| drum(y) - ball(y));
    let under = |y: f64| r - r.mul_add(r, -(y - r).powi(2)).max(0.0).sqrt();
    let overlap = integral(0.0, top.min(r), 200_000, |y| {
        under(y) * (drum(y).min(end) - ball(y)).max(0.0)
    });
    bottom + 10.0 * sliver - overlap
}

/// The front bottom edge of the shaved cube and the upright edge it ends
/// at, asked for in either order. Rounded first, the bottom blend caps
/// flush at its vertex, and that cap meets the drum along the upright's own
/// line, a corner of the cap's that is no piece of the upright: the upright
/// is found where it still lies between the front plane and the drum and
/// rounds the full height through the bottom blend. Rounded second, the
/// bottom blend runs on through the upright's band only as far as its own
/// vertex, since the drum leaves that vertex obliquely and the material
/// past it is no part of this edge's rounding. Both orders take off the
/// union of the two blends.
#[test]
fn an_edge_meeting_a_capped_blend_at_a_drum_is_not_its_cap() {
    let mut model = Model::new();
    let part = shaved_cube(&mut model);
    let before = volume(&model, &part);
    let end = 5.0 + 11.0_f64.sqrt();
    let bottom = edge_near(&model, &part, Point::new(5.0, 0.0, 0.0));
    let upright = edge_near(&model, &part, Point::new(end, 0.0, 5.0));
    for r in [0.18, 0.5, 1.0] {
        for (order, edges) in [
            ("bottom first", [bottom.clone(), upright.clone()]),
            ("upright first", [upright.clone(), bottom.clone()]),
        ] {
            let mut copy = model.clone();
            let rounded = ogeom::fillet::fillet_edges(&mut copy, &part, &edges, r, T)
                .unwrap()
                .shape;
            let diagnosis = ogeom::algo::check(&copy, &rounded, T).unwrap();
            assert!(
                diagnosis.is_valid(),
                "r {r}, {order}: {:?}",
                diagnosis.problems
            );
            let removed = before - volume(&copy, &rounded);
            let want = shaved_corner_removed(r);
            assert!(
                (removed - want).abs() < want * 1e-7,
                "r {r}, {order}: took off {removed}, the union is {want}"
            );
        }
    }
}
