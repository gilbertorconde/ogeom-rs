//! The ways a part lies on a flat bed, against closed forms: a box, a
//! regular tetrahedron, a cylinder, a cone, a sphere, an L-bracket whose
//! centre of mass falls outside two of its faces, and a slot-shaped bar
//! whose flat sides end in rounds; each as a solid and, where the mesh is
//! the point, as a mesh built in code.
#![allow(clippy::unwrap_used, clippy::cast_precision_loss, reason = "test code")]

use core::f64::consts::PI;

use ogeom_algo::{
    make_box, make_cone, make_cylinder, make_face_with_pcurves, make_polygon, make_polyhedron,
    make_prism, make_sphere,
};
use ogeom_core::Tolerances;
use ogeom_fit::{Rest, resting_faces, resting_faces_of};
use ogeom_geom::{PlaneSurface, SurfaceGeometry};
use ogeom_math::{Frame, Plane, Point, Point2, Vector};
use ogeom_topo::{Model, Shape, Triangulation};

const T: Tolerances = Tolerances::millimetres();

fn down(rest: &Rest) -> Vector {
    rest.down.vector()
}

fn close(a: Vector, b: Vector) -> bool {
    (a - b).magnitude() < 1e-9
}

/// The shoelace area of a footprint.
fn polygon_area(points: &[Point2]) -> f64 {
    let n = points.len();
    (0..n)
        .map(|i| {
            let (a, b) = (points[i], points[(i + 1) % n]);
            a.x * b.y - b.x * a.y
        })
        .sum::<f64>()
        / 2.0
}

/// A closed mesh of a prism over a convex profile in `xy`, `height` tall,
/// the caps fanned from the profile's first point.
fn prism_mesh(profile: &[(f64, f64)], height: f64) -> Triangulation {
    let n = u32::try_from(profile.len()).unwrap();
    let mut mesh = Triangulation::new();
    for z in [0.0, height] {
        for &(x, y) in profile {
            mesh.positions.push(Point::new(x, y, z));
        }
    }
    for i in 1..n - 1 {
        mesh.triangles.push([0, i + 1, i]);
        mesh.triangles.push([n, n + i, n + i + 1]);
    }
    for i in 0..n {
        let j = (i + 1) % n;
        mesh.triangles.push([i, j, n + j]);
        mesh.triangles.push([i, n + j, n + i]);
    }
    mesh
}

/// A 30 × 20 × 10 box: six rests, all stable, the largest pair first and
/// the opposite faces of each pair equal, the one facing −z, −y or −x
/// first within its pair as the tie rule says.
#[test]
fn a_box_lies_six_ways() {
    let mut model = Model::new();
    let solid = make_box(&mut model, Frame::WORLD, (30.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let mesh =
        ogeom_mesh::triangulate(&model, &solid, ogeom_mesh::Deflection::default(), T).unwrap();
    for rests in [
        resting_faces_of(&model, &solid, T).unwrap(),
        resting_faces(&mesh, T).unwrap(),
    ] {
        assert_eq!(rests.len(), 6, "{rests:#?}");
        // (down, area, height, margin): the margin is the nearer half side.
        let expected = [
            (Vector::new(0.0, 0.0, -1.0), 600.0, 10.0, 10.0),
            (Vector::new(0.0, 0.0, 1.0), 600.0, 10.0, 10.0),
            (Vector::new(0.0, -1.0, 0.0), 300.0, 20.0, 5.0),
            (Vector::new(0.0, 1.0, 0.0), 300.0, 20.0, 5.0),
            (Vector::new(-1.0, 0.0, 0.0), 200.0, 30.0, 5.0),
            (Vector::new(1.0, 0.0, 0.0), 200.0, 30.0, 5.0),
        ];
        for (rest, (d, area, height, margin)) in rests.iter().zip(expected) {
            assert!(close(down(rest), d), "{:?} for {d:?}", rest.down);
            assert!(rest.stable);
            assert!((rest.area - area).abs() < 1e-9, "{}", rest.area);
            assert!((rest.height - height).abs() < 1e-9, "{}", rest.height);
            assert!((rest.centre_height - height / 2.0).abs() < 1e-9);
            assert!((rest.margin - margin).abs() < 1e-9, "{}", rest.margin);
            assert!((rest.tip - (margin / (height / 2.0)).atan()).abs() < 1e-12);
            assert!((polygon_area(&rest.footprint) - area).abs() < 1e-9);
        }
    }
}

/// A regular tetrahedron of edge 10, as a solid and as its four triangles:
/// four rests, each a face of area `√3/4 a²`, the part `a √(2/3)` tall,
/// its centre at a quarter of that and the margin the face's inradius.
#[test]
fn a_regular_tetrahedron_lies_on_each_face() {
    let a = 10.0;
    let corners = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(a, 0.0, 0.0),
        Point::new(a / 2.0, a * 3.0_f64.sqrt() / 2.0, 0.0),
        Point::new(
            a / 2.0,
            a / (2.0 * 3.0_f64.sqrt()),
            a * (2.0_f64 / 3.0).sqrt(),
        ),
    ];
    let mut model = Model::new();
    let rings = vec![vec![0, 1, 2], vec![0, 1, 3], vec![1, 2, 3], vec![0, 2, 3]];
    let solid = make_polyhedron(&mut model, &corners, &rings, T)
        .unwrap()
        .shape;
    let mut mesh = Triangulation::new();
    mesh.positions = corners.to_vec();
    mesh.triangles = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]];

    let area = 3.0_f64.sqrt() / 4.0 * a * a;
    let height = a * (2.0_f64 / 3.0).sqrt();
    let inradius = a / (2.0 * 3.0_f64.sqrt());
    for rests in [
        resting_faces_of(&model, &solid, T).unwrap(),
        resting_faces(&mesh, T).unwrap(),
    ] {
        assert_eq!(rests.len(), 4);
        for rest in &rests {
            assert!(rest.stable);
            assert!((rest.area - area).abs() < 1e-9, "{}", rest.area);
            assert!((rest.height - height).abs() < 1e-9, "{}", rest.height);
            assert!((rest.centre_height - height / 4.0).abs() < 1e-9);
            assert!((rest.margin - inradius).abs() < 1e-9, "{}", rest.margin);
            assert_eq!(rest.footprint.len(), 3);
        }
    }
}

/// A cylinder of radius 5 and height 20 stands on either cap and does not
/// lie on its side, where it rolls: two rests. The solid's contact area is
/// the disc's; the mesh's, of its 64-gon, within the chord's share of it,
/// and the two agree on the rest within the mesh's chord.
#[test]
fn a_cylinder_stands_on_its_caps_and_rolls_on_its_side() {
    let (r, h) = (5.0, 20.0);
    let mut model = Model::new();
    let solid = make_cylinder(&mut model, Frame::WORLD, r, h, T)
        .unwrap()
        .shape;
    let sides = 64;
    let profile: Vec<(f64, f64)> = (0..sides)
        .map(|k| {
            let t = 2.0 * PI * f64::from(k) / f64::from(sides);
            (r * t.cos(), r * t.sin())
        })
        .collect();
    let mesh = prism_mesh(&profile, h);
    // The solid's own mesh turns by at most 0.25 rad a facet, coarser
    // than the 64-gon: its sagitta bounds both.
    let chord = r * (1.0 - 0.125_f64.cos());

    let exact = resting_faces_of(&model, &solid, T).unwrap();
    let meshed = resting_faces(&mesh, T).unwrap();
    for rests in [&exact, &meshed] {
        assert_eq!(rests.len(), 2, "{rests:#?}");
        assert!(close(down(&rests[0]), Vector::new(0.0, 0.0, -1.0)));
        assert!(close(down(&rests[1]), Vector::new(0.0, 0.0, 1.0)));
        for rest in rests.iter() {
            assert!(rest.stable);
            assert!((rest.height - h).abs() < 1e-9);
            assert!((rest.centre_height - h / 2.0).abs() < 1e-6);
            assert!((rest.margin - r).abs() <= chord, "{}", rest.margin);
        }
    }
    for (e, m) in exact.iter().zip(&meshed) {
        assert!((e.area - PI * r * r).abs() < 1e-6, "{}", e.area);
        // The polygon's area falls short of the disc's by at most the
        // chord along the perimeter.
        assert!(
            m.area < e.area && e.area - m.area <= 2.0 * PI * r * chord,
            "{}",
            m.area
        );
        assert!((e.margin - m.margin).abs() <= chord);
    }
}

/// A cone of base radius 10 and height 20 stands on its base and nowhere
/// else: not on its tip, which is no plane of its hull, and not on its
/// side, where it rolls round its tip.
#[test]
fn a_cone_stands_on_its_base_only() {
    let mut model = Model::new();
    let solid = make_cone(&mut model, Frame::WORLD, 10.0, 0.0, 20.0, T)
        .unwrap()
        .shape;
    let rests = resting_faces_of(&model, &solid, T).unwrap();
    assert_eq!(rests.len(), 1, "{rests:#?}");
    let base = &rests[0];
    assert!(close(down(base), Vector::new(0.0, 0.0, -1.0)));
    assert!(base.stable);
    assert!((base.area - 100.0 * PI).abs() < 1e-6, "{}", base.area);
    assert!(
        (base.centre_height - 5.0).abs() < 1e-6,
        "{}",
        base.centre_height
    );
    assert!((base.height - 20.0).abs() < 1e-6, "{}", base.height);
}

/// A sphere has no flat pose: every facet of its mesh rolls.
#[test]
fn a_sphere_has_no_rest() {
    let mut model = Model::new();
    let solid = make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    assert!(resting_faces_of(&model, &solid, T).unwrap().is_empty());
    let mesh =
        ogeom_mesh::triangulate(&model, &solid, ogeom_mesh::Deflection::default(), T).unwrap();
    assert!(resting_faces(&mesh, T).unwrap().is_empty());
}

/// An L profile 20 wide and 50 tall, 2 thick, extruded 10: its centre of
/// mass `(62/17, 317/17, 5)` is past the 2 mm end of the foot and past the
/// 2 mm top of the upright, so those two faces tip; the inside faces of the
/// L are no rests at all, and the bed across the L's opening stands on the
/// two tips' edges with no contact area.
#[test]
fn an_l_brackets_end_faces_tip_and_its_inside_faces_are_no_rests() {
    let mut model = Model::new();
    let pts = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(20.0, 0.0, 0.0),
        Point::new(20.0, 2.0, 0.0),
        Point::new(2.0, 2.0, 0.0),
        Point::new(2.0, 50.0, 0.0),
        Point::new(0.0, 50.0, 0.0),
    ];
    let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
    let surface = PlaneSurface::over(Plane::XY, (-1.0, 21.0), (-1.0, 51.0)).unwrap();
    let edges = model.children_of(&wire).unwrap();
    let face = make_face_with_pcurves(&mut model, SurfaceGeometry::Plane(surface), &[edges], T)
        .unwrap()
        .shape;
    let solid: Shape = make_prism(&mut model, &face, Vector::new(0.0, 0.0, 10.0), T)
        .unwrap()
        .shape;

    let (cx, cy) = (62.0 / 17.0, 317.0 / 17.0);
    let across = Vector::new(48.0, 18.0, 0.0) / (48.0_f64.hypot(18.0));
    // The centre's foot along the opening's edge, from the foot's tip.
    let along = Vector::new(-18.0, 48.0, 0.0).dot(Vector::new(cx - 20.0, cy - 2.0, 0.0))
        / 48.0_f64.hypot(18.0);
    assert!(along > 5.0 && 48.0_f64.hypot(18.0) - along > 5.0);
    // (down, area, centre height, stable, margin), in the documented order.
    let expected = [
        (Vector::new(-1.0, 0.0, 0.0), 500.0, cx, true, 5.0),
        (Vector::new(0.0, -1.0, 0.0), 200.0, cy, true, cx),
        (Vector::new(0.0, 0.0, -1.0), 136.0, 5.0, true, cx),
        (Vector::new(0.0, 0.0, 1.0), 136.0, 5.0, true, cx),
        (
            across,
            0.0,
            across.dot(Vector::new(20.0 - cx, 2.0 - cy, 0.0)),
            true,
            5.0,
        ),
        (Vector::new(1.0, 0.0, 0.0), 20.0, 20.0 - cx, false, 2.0 - cy),
        (Vector::new(0.0, 1.0, 0.0), 20.0, 50.0 - cy, false, 2.0 - cx),
    ];
    let rests = resting_faces_of(&model, &solid, T).unwrap();
    assert_eq!(rests.len(), expected.len(), "{rests:#?}");
    for (rest, (d, area, centre, stable, margin)) in rests.iter().zip(expected) {
        assert!(close(down(rest), d), "{:?} for {d:?}", rest.down);
        assert_eq!(rest.stable, stable, "{d:?}");
        assert!((rest.area - area).abs() < 1e-9, "{d:?}: {}", rest.area);
        assert!(
            (rest.centre_height - centre).abs() < 1e-9,
            "{d:?}: {}",
            rest.centre_height
        );
        assert!(
            (rest.margin - margin).abs() < 1e-9,
            "{d:?}: {}",
            rest.margin
        );
    }
}

/// A slot-shaped bar, two flats 40 long joined by half rounds of radius 5
/// meshed in 32 facets each, 8 tall: it lies on its caps and on either
/// flat, whose ends run smoothly into the rounds but which it has to be
/// tipped 76° to leave; the rounds' facets roll.
#[test]
fn a_flat_ending_in_rounds_is_a_rest_and_the_rounds_roll() {
    let (length, r, h) = (40.0, 5.0, 8.0);
    let mut profile = Vec::new();
    for (cx, from) in [(length, -PI / 2.0), (0.0, PI / 2.0)] {
        for k in 0..=32 {
            let t = from + PI * f64::from(k) / 32.0;
            profile.push((cx + r * t.cos(), r * t.sin()));
        }
    }
    let rests = resting_faces(&prism_mesh(&profile, h), T).unwrap();
    assert_eq!(rests.len(), 4, "{rests:#?}");
    let flats: Vec<&Rest> = rests.iter().filter(|r| down(r).z.abs() < 0.5).collect();
    assert_eq!(flats.len(), 2);
    for flat in flats {
        assert!(flat.stable);
        assert!((flat.area - length * h).abs() < 1e-9, "{}", flat.area);
        assert!((flat.margin - h / 2.0).abs() < 1e-9);
        assert!((flat.centre_height - r).abs() < 1e-9);
    }
}
