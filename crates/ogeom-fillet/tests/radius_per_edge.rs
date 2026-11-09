//! A radius per edge in one call.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{fillet_edges, fillet_edges_with};
use ogeom_geom::{CircleCurve, Curve3d as _};
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn bx(m: &mut Model, s: [f64; 3]) -> Shape {
    ogeom_algo::make_box(m, at(0.0, 0.0, 0.0), (s[0], s[1], s[2]), T)
        .unwrap()
        .shape
}

fn vol(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(
        m,
        s,
        ogeom_mesh::Deflection {
            chord: 1e-3,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

/// The edge of `s` whose middle stands at `p`.
fn edge_at(m: &Model, s: &Shape, p: [f64; 3]) -> Shape {
    let q = Point::new(p[0], p[1], p[2]);
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            c.point_at(0.5 * (a + b), T).unwrap().distance(q) < 1e-9
        })
        .expect("an edge there")
}

/// The section a band of radius `r` takes off a right-angled crease.
fn section(r: f64) -> f64 {
    (1.0 - PI / 4.0) * r * r
}

#[test]
fn one_radius_on_every_edge_is_fillet_edges() {
    let mut m = Model::new();
    let b = bx(&mut m, [20.0, 20.0, 10.0]);
    let all = explore_unique(&m, &b, ShapeType::Edge).unwrap();
    let mut plain = m.clone();
    let a = fillet_edges(&mut plain, &b, &all, 2.0, T).unwrap();
    let with: Vec<(Shape, f64)> = all.iter().map(|e| (e.clone(), 2.0)).collect();
    let w = fillet_edges_with(&mut m, &b, &with, T).unwrap();
    let faces = |m: &Model, s: &Shape| explore_unique(m, s, ShapeType::Face).unwrap().len();
    assert_eq!(faces(&plain, &a.shape), faces(&m, &w.shape));
    assert_eq!(vol(&plain, &a.shape).to_bits(), vol(&m, &w.shape).to_bits());
}

/// Two creases of a box's top turning at a vertex, at radii 1 and 3: the
/// cut takes the union of the two bands, whichever order they are asked
/// in. Where the bands overlap at the vertex, the region under the
/// small band's spandrel within the large band's reach,
/// `∫ (1 − √(1 − (z − 9)²)) (3 − √(9 − (z − 7)²)) dz` over the top unit.
#[test]
fn two_radii_turning_at_a_vertex_take_the_union_of_their_bands() {
    let mut m = Model::new();
    let b = bx(&mut m, [20.0, 20.0, 10.0]);
    let front = edge_at(&m, &b, [10.0, 0.0, 10.0]);
    let left = edge_at(&m, &b, [0.0, 10.0, 10.0]);
    let n = 200_000;
    let h = 1.0 / f64::from(n);
    let overlap: f64 = (0..n)
        .map(|k| {
            let z = 9.0 + (f64::from(k) + 0.5) * h;
            let small = 1.0 - (1.0 - (z - 9.0).powi(2)).sqrt();
            let large = 3.0 - (9.0 - (z - 7.0).powi(2)).sqrt();
            small * large * h
        })
        .sum();
    let expected = 4000.0 - 20.0 * section(1.0) - 20.0 * section(3.0) + overlap;
    let mut volumes = Vec::new();
    for asked in [
        [(front.clone(), 1.0), (left.clone(), 3.0)],
        [(left.clone(), 3.0), (front.clone(), 1.0)],
    ] {
        let mut mm = m.clone();
        let built = fillet_edges_with(&mut mm, &b, &asked, T).unwrap();
        let diagnosis = ogeom_algo::check(&mm, &built.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        assert!(built.history.is_deleted(&front) && built.history.is_deleted(&left));
        volumes.push(vol(&mm, &built.shape));
    }
    for v in &volumes {
        assert!((v - expected).abs() < 1e-3, "{v} against {expected}");
    }
}

/// A box's three edges at one corner at the issue's three radii, 1, 3
/// and 2: no ball rolls round it, so each band stops where its section
/// clears the others and quads round a centre close the corner, tangent
/// to every band end and face, the same whatever the order asked.
#[test]
fn three_radii_at_a_corner_close_with_one_patch_in_any_order() {
    let mut m = Model::new();
    let b = bx(&mut m, [20.0, 20.0, 10.0]);
    let edges = [
        edge_at(&m, &b, [10.0, 0.0, 10.0]),
        edge_at(&m, &b, [0.0, 10.0, 10.0]),
        edge_at(&m, &b, [0.0, 0.0, 5.0]),
    ];
    let radii = [1.0, 3.0, 2.0];
    let asked: Vec<(Shape, f64)> = edges.iter().cloned().zip(radii).collect();
    let mut volumes = Vec::new();
    let mut faces = Vec::new();
    for order in [[0, 1, 2], [2, 0, 1], [1, 2, 0]] {
        let request: Vec<(Shape, f64)> = order.iter().map(|&i| asked[i].clone()).collect();
        let rounded = fillet_edges_with(&mut m, &b, &request, T).unwrap();
        let found = ogeom_algo::check(&m, &rounded.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        // The quads: the corner's B-spline faces, tangent across every edge.
        let quads: Vec<Shape> = explore_unique(&m, &rounded.shape, ShapeType::Face)
            .unwrap()
            .into_iter()
            .filter(|f| {
                m.node(f)
                    .and_then(|n| n.data().as_face())
                    .and_then(|d| m.geometry().surface(d.surface))
                    .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::BSpline(_)))
            })
            .collect();
        assert_eq!(quads.len(), 6, "a quad at each corner of the opening");
        let mut outside = 0;
        for quad in &quads {
            for c in ogeom_fillet::analyse_blend(&m, &rounded.shape, quad, 17, T).unwrap() {
                assert!(c.tangency_error < 2e-4, "a crease across the corner: {c:?}");
                if !quads.iter().any(|q| q.is_same(&c.neighbour)) {
                    outside += 1;
                }
            }
        }
        assert_eq!(
            outside, 12,
            "three band ends and three faces, each in halves"
        );
        volumes.push(vol(&m, &rounded.shape));
        faces.push(
            explore_unique(&m, &rounded.shape, ShapeType::Face)
                .unwrap()
                .len(),
        );
    }
    for (v, f) in volumes.iter().zip(&faces) {
        assert!((v - volumes[0]).abs() <= volumes[0] * 1e-9, "{volumes:?}");
        assert_eq!(*f, faces[0], "{faces:?}");
    }
    let mut at = |r: f64| {
        let request: Vec<(Shape, f64)> = edges.iter().cloned().map(|e| (e, r)).collect();
        let s = fillet_edges_with(&mut m, &b, &request, T).unwrap().shape;
        vol(&m, &s)
    };
    let (one, three) = (at(1.0), at(3.0));
    assert!(
        three < volumes[0] && volumes[0] < one,
        "{three} < {} < {one}",
        volumes[0]
    );
}

/// A stadium prism: a rectangle with semicircular ends, extruded up.
fn stadium(model: &mut Model, length: f64, r: f64, height: f64) -> Shape {
    let tl = ogeom_algo::make_vertex(model, Point::new(0.0, r, 0.0)).shape;
    let tr = ogeom_algo::make_vertex(model, Point::new(length, r, 0.0)).shape;
    let br = ogeom_algo::make_vertex(model, Point::new(length, -r, 0.0)).shape;
    let bl = ogeom_algo::make_vertex(model, Point::new(0.0, -r, 0.0)).shape;
    let arc = |model: &mut Model, centre: Point, x: Direction, from: &Shape, to: &Shape| {
        let frame = Frame::new(centre, Direction::Z, x, T).unwrap();
        let circle = Circle::new(frame, r, T).unwrap();
        let curve = ogeom_geom::Curve::Circle(CircleCurve::new(circle));
        ogeom_algo::make_edge_between(model, curve, (0.0, PI), from, to, T)
            .unwrap()
            .shape
    };
    let seg = |model: &mut Model, from: (&Shape, Point), to: (&Shape, Point)| {
        let line = ogeom_geom::LineCurve::segment(from.1, to.1, T).unwrap();
        let curve = ogeom_geom::Curve::Line(line);
        let domain = ogeom_geom::Curve3d::domain(&curve);
        ogeom_algo::make_edge_between(model, curve, domain, from.0, to.0, T)
            .unwrap()
            .shape
    };
    let top = seg(
        model,
        (&tl, Point::new(0.0, r, 0.0)),
        (&tr, Point::new(length, r, 0.0)),
    );
    let right = arc(
        model,
        Point::new(length, 0.0, 0.0),
        Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(),
        &br,
        &tr,
    );
    let bottom = seg(
        model,
        (&br, Point::new(length, -r, 0.0)),
        (&bl, Point::new(0.0, -r, 0.0)),
    );
    let left = arc(model, Point::new(0.0, 0.0, 0.0), Direction::Y, &tl, &bl);
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face_with_pcurves(
        model,
        plane.into(),
        &[vec![left, bottom.reversed(), right, top.reversed()]],
        T,
    )
    .unwrap()
    .shape;
    ogeom_algo::make_prism(model, &face, Vector::new(0.0, 0.0, height), T)
        .unwrap()
        .shape
}

/// A stadium's top rim, its straight runs at 1 and its ends at 2: the
/// radius steps at each tangent junction, each band ending in the plane
/// square to the rim there. The removed ring is each band's section swept
/// along its own run, the ends' by Pappus about the section's centroid.
#[test]
fn a_tangent_chain_steps_its_radius_at_the_junctions() {
    let (length, r, height) = (10.0, 5.0, 4.0);
    let mut m = Model::new();
    let s = stadium(&mut m, length, r, height);
    let before = vol(&m, &s);
    let asked = [
        (edge_at(&m, &s, [length / 2.0, r, height]), 1.0),
        (edge_at(&m, &s, [length + r, 0.0, height]), 2.0),
        (edge_at(&m, &s, [length / 2.0, -r, height]), 1.0),
        (edge_at(&m, &s, [-r, 0.0, height]), 2.0),
    ];
    let built = fillet_edges_with(&mut m, &s, &asked, T).unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let centroid = |b: f64| (10.0 - 3.0 * PI) / (12.0 - 3.0 * PI) * b;
    let expected =
        before - 2.0 * length * section(1.0) - 2.0 * PI * (r - centroid(2.0)) * section(2.0);
    let v = vol(&m, &built.shape);
    assert!((v - expected).abs() < 0.05, "{v} against {expected}");
}
