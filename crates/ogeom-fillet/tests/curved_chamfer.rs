//! Chamfers of solids with curved faces through the blend engine: what the
//! history tells of each element, beside the M3 corpus's volumes.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn up_at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(
        Point::new(x, y, z),
        Direction::new(Vector::new(0.0, 0.0, 1.0), T).unwrap(),
        Direction::new(Vector::new(1.0, 0.0, 0.0), T).unwrap(),
        T,
    )
    .unwrap()
}

fn point_of(m: &Model, v: &Shape) -> Point {
    m.node(v)
        .and_then(|n| n.data().as_vertex().map(|d| d.point))
        .unwrap()
}

fn surface_kind(m: &Model, face: &Shape) -> ogeom_geom::SurfaceKind {
    use ogeom_geom::Surface as _;
    let data = m.node(face).unwrap().data().as_face().unwrap();
    m.geometry().surface(data.surface).unwrap().kind()
}

/// Every edge chamfered is deleted and generates one face on `bevel`, every
/// face of the solid is modified into one face, and the result is valid.
fn history_holds(
    m: &Model,
    solid: &Shape,
    built: &ogeom_algo::Built,
    edges: &[Shape],
    bevel: ogeom_geom::SurfaceKind,
) {
    assert!(ogeom_algo::check(m, &built.shape, T).unwrap().is_valid());
    for e in edges {
        assert!(built.history.is_deleted(e), "a chamfered edge left");
        let made = built.history.generated(e);
        assert_eq!(made.len(), 1, "a chamfered edge without its one bevel");
        assert_eq!(
            surface_kind(m, &made[0]),
            bevel,
            "a bevel on another surface"
        );
    }
    for f in explore_unique(m, solid, ShapeType::Face).unwrap() {
        assert_eq!(
            built.history.modified(&f).len(),
            1,
            "a face not carried to one face"
        );
    }
}

/// A tube of radius 3 bored to 2, 4 tall, its four rims chamfered 0.3:
/// each rim is deleted into one cone, each face carried to one face.
#[test]
fn a_tube_s_rims_chamfer_into_cones_their_faces_carried() {
    let mut m = Model::new();
    let outer = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 4.0, T)
        .unwrap()
        .shape;
    let inner = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, -1.0), 2.0, 6.0, T)
        .unwrap()
        .shape;
    let tube = ogeom_bool::cut(&mut m, &outer, &inner, T).unwrap().shape;
    let rims: Vec<Shape> = explore_unique(&m, &tube, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            matches!(
                ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap().0,
                ogeom_geom::Curve::Circle(_)
            )
        })
        .collect();
    assert!(rims.len() >= 4);
    let built = ogeom_fillet::chamfer_edges(&mut m, &tube, &rims, 0.3, T).unwrap();
    history_holds(&m, &tube, &built, &rims, ogeom_geom::SurfaceKind::Cone);
}

/// The D-shaped block (a drum of radius 3 cut in half along its axis, 10
/// long), its two straight edges chamfered 0.5: each edge deleted into one
/// plane, each face carried, and each cap's corners where the bevel ends
/// gone, the cap running round the bevel's end there.
#[test]
fn a_d_block_s_straight_edges_chamfer_into_planes_ending_on_the_caps() {
    let mut m = Model::new();
    let full = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 10.0, T)
        .unwrap()
        .shape;
    let half = ogeom_algo::make_box(&mut m, up_at(-4.0, -4.0, -1.0), (4.0, 8.0, 12.0), T)
        .unwrap()
        .shape;
    let d = ogeom_bool::cut(&mut m, &full, &half, T).unwrap().shape;
    let straight: Vec<Shape> = explore_unique(&m, &d, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
            let (p, q) = (point_of(&m, &a), point_of(&m, &b));
            p.x.abs() < 1e-9 && q.x.abs() < 1e-9 && (p.z - q.z).abs() > 1.0
        })
        .collect();
    assert_eq!(straight.len(), 2);
    let built = ogeom_fillet::chamfer_edges(&mut m, &d, &straight, 0.5, T).unwrap();
    history_holds(&m, &d, &built, &straight, ogeom_geom::SurfaceKind::Plane);
    for e in &straight {
        let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
        for v in [a, b] {
            assert!(
                built.history.is_deleted(&v),
                "a corner the bevel takes left"
            );
        }
    }
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(
        m,
        s,
        ogeom_mesh::Deflection {
            chord: 1e-5,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

/// The D-shaped block's straight edges chamfered 0.5 along the flat at 45
/// degrees from it: in the section the bevel leaves `(0, 2.5)` towards
/// the edge and into the material, meeting the wall of radius 3 where
/// `2w^2 + 5w - 2.75 = 0` for `w` its run along each axis. The half disc
/// less, twice along the 10, the triangle from the corner to the two
/// contacts and the wall's segment beyond the chord.
#[test]
fn a_d_block_s_edges_chamfer_by_distance_and_angle_on_the_flat() {
    let mut m = Model::new();
    let full = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 10.0, T)
        .unwrap()
        .shape;
    let half = ogeom_algo::make_box(&mut m, up_at(-4.0, -4.0, -1.0), (4.0, 8.0, 12.0), T)
        .unwrap()
        .shape;
    let d = ogeom_bool::cut(&mut m, &full, &half, T).unwrap().shape;
    let flat = explore_unique(&m, &d, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let (_, n) = ogeom_algo::face_normal(&m, f, T).unwrap();
            n.x < -0.5
        })
        .unwrap();
    let specs: Vec<(Shape, ogeom_fillet::Chamfer)> = explore_unique(&m, &d, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
            let (p, q) = (point_of(&m, &a), point_of(&m, &b));
            p.x.abs() < 1e-9 && q.x.abs() < 1e-9 && (p.z - q.z).abs() > 1.0
        })
        .map(|e| {
            (
                e,
                ogeom_fillet::Chamfer::Angle {
                    face: flat.clone(),
                    distance: 0.5,
                    angle: core::f64::consts::FRAC_PI_4,
                },
            )
        })
        .collect();
    assert_eq!(specs.len(), 2);
    let built = ogeom_fillet::chamfer_edges_with(&mut m, &d, &specs, T).unwrap();
    assert!(ogeom_algo::check(&m, &built.shape, T).unwrap().is_valid());
    for (e, _) in &specs {
        assert_eq!(
            built.history.generated(e).len(),
            1,
            "an edge without its bevel"
        );
    }
    let w = (-5.0 + (25.0_f64 + 22.0).sqrt()) / 4.0;
    let (corner, on_flat, on_wall) = ((0.0, 3.0), (0.0, 2.5), (w, 2.5 + w));
    let triangle = ((on_flat.0 - corner.0) * (on_wall.1 - corner.1)
        - (on_flat.1 - corner.1) * (on_wall.0 - corner.0))
        .abs()
        / 2.0;
    let phi = (w / 3.0).asin();
    let segment = 9.0 / 2.0 * (phi - phi.sin());
    let want = (core::f64::consts::PI * 9.0 / 2.0 - 2.0 * (triangle + segment)) * 10.0;
    let got = volume(&m, &built.shape);
    assert!((got - want).abs() < 1e-6, "volume {got} against {want}");
}

/// A drum's top rim chamfered 0.3 along the cap at 45 degrees from it:
/// on its square rim the same bevel as 0.3 along both faces, the same
/// solid's volume.
#[test]
fn a_rim_chamfers_by_distance_and_angle_as_its_square_corner_says() {
    let make = |m: &mut Model| {
        ogeom_algo::make_cylinder(m, up_at(0.0, 0.0, 0.0), 3.0, 4.0, T)
            .unwrap()
            .shape
    };
    let rim_and_cap = |m: &Model, s: &Shape| {
        let rim = explore_unique(m, s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let (a, _) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
                matches!(
                    ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap().0,
                    ogeom_geom::Curve::Circle(_)
                ) && point_of(m, &a).z > 3.0
            })
            .unwrap();
        let cap = explore_unique(m, s, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| ogeom_algo::face_normal(m, f, T).unwrap().1.z > 0.5)
            .unwrap();
        (rim, cap)
    };
    let mut m = Model::new();
    let drum = make(&mut m);
    let (rim, cap) = rim_and_cap(&m, &drum);
    let angled = ogeom_fillet::chamfer_edges_with(
        &mut m,
        &drum,
        &[(
            rim.clone(),
            ogeom_fillet::Chamfer::Angle {
                face: cap,
                distance: 0.3,
                angle: core::f64::consts::FRAC_PI_4,
            },
        )],
        T,
    )
    .unwrap();
    let square = ogeom_fillet::chamfer_edges(&mut m, &drum, &[rim], 0.3, T).unwrap();
    assert!(ogeom_algo::check(&m, &angled.shape, T).unwrap().is_valid());
    let (a, b) = (volume(&m, &angled.shape), volume(&m, &square.shape));
    assert!((a - b).abs() < 1e-9 * b, "{a} against {b}");
}

/// The D-shaped block (a drum of radius 3 halved along its axis, 10 long).
fn d_block(m: &mut Model) -> Shape {
    let full = ogeom_algo::make_cylinder(m, up_at(0.0, 0.0, 0.0), 3.0, 10.0, T)
        .unwrap()
        .shape;
    let half = ogeom_algo::make_box(m, up_at(-4.0, -4.0, -1.0), (4.0, 8.0, 12.0), T)
        .unwrap()
        .shape;
    ogeom_bool::cut(m, &full, &half, T).unwrap().shape
}

/// `s` with what lies past the plane through `at` of outward normal `n`
/// cut away, by a box standing on the plane.
fn cut_past(m: &mut Model, s: &Shape, at: Point, n: Vector) -> Shape {
    let z = n / n.magnitude();
    let x = if z.x.abs() < 0.9 {
        Vector::X
    } else {
        Vector::Y
    };
    let x = (x - z * x.dot(z)) / (x - z * x.dot(z)).magnitude();
    let y = z.cross(x);
    let corner = at - x * 30.0 - y * 30.0;
    let frame = Frame::new(
        corner,
        Direction::new(z, T).unwrap(),
        Direction::new(x, T).unwrap(),
        T,
    )
    .unwrap();
    let cutter = ogeom_algo::make_box(m, frame, (60.0, 60.0, 30.0), T)
        .unwrap()
        .shape;
    ogeom_bool::cut(m, s, &cutter, T).unwrap().shape
}

/// The D-shaped block with every edge chamfered 0.5: each corner of three
/// bevels (the straight edge, the cap's straight edge and the cap's arc)
/// closes with the facet through its faces' corners. The block is convex
/// and every bevel convex, so the result is the block cut by each bevel's
/// plane or, at a rim, its cone, and by each facet's plane: here cut by
/// booleans, which share nothing with the engine's assembly, to the same
/// volume. Each corner generates one triangular facet; where a cap's arc is
/// split at the wall's seam its bevels run on into each other, with none.
#[test]
fn a_d_block_chamfered_on_every_edge_closes_each_corner_with_its_facet() {
    let d = 0.5_f64;
    let mut m = Model::new();
    let block = d_block(&mut m);
    let edges = explore_unique(&m, &block, ShapeType::Edge).unwrap();
    let vertices = explore_unique(&m, &block, ShapeType::Vertex).unwrap();
    let built = ogeom_fillet::chamfer_edges(&mut m, &block, &edges, d, T).unwrap();
    assert!(ogeom_algo::check(&m, &built.shape, T).unwrap().is_valid());
    let back: Vec<Shape> = edges.iter().rev().cloned().collect();
    let again = ogeom_fillet::chamfer_edges(&mut m, &block, &back, d, T).unwrap();
    assert!(
        (volume(&m, &built.shape) - volume(&m, &again.shape)).abs() < 1e-9,
        "the order changes the solid"
    );
    // The four corners on the flat take a facet each; where a cap's arc is
    // split at the wall's seam, its two bevels run on into each other and
    // take none.
    let mut facets = 0;
    for v in &vertices {
        let made = built.history.generated(v);
        if point_of(&m, v).x.abs() < 1e-9 {
            assert_eq!(made.len(), 1, "a corner without its facet");
            let corners = explore_unique(&m, &made[0], ShapeType::Vertex).unwrap();
            assert_eq!(corners.len(), 3, "a facet that is not a triangle");
            facets += 1;
        } else {
            assert!(made.is_empty(), "a facet where two bevels run on");
        }
    }
    assert_eq!(facets, 4);

    // The oracle: the same planes and cones, by booleans.
    let mut o = Model::new();
    let s = d_block_cut_by_every_bevel(&mut o, d, true);
    let (got, want) = (volume(&m, &built.shape), volume(&o, &s));
    assert!(
        (got - want).abs() < 1e-6,
        "volume {got} against the booleans' {want}"
    );
}

/// The D-shaped block with every edge bevelled `d` by booleans: in common
/// with each rim's cone, cut past each straight edge's bevel plane and,
/// with `facets`, past each corner's facet plane.
fn d_block_cut_by_every_bevel(o: &mut Model, d: f64, facets: bool) -> Shape {
    let big = 3.0_f64;
    let mut s = d_block(o);
    let theta = d / big;
    let (wx, wy) = (big * theta.sin(), big * theta.cos());
    // The rims' cones: the material inside each.
    let top = ogeom_algo::make_cone(o, up_at(0.0, 0.0, -1.0), 13.5, 2.5, 11.0, T)
        .unwrap()
        .shape;
    s = ogeom_bool::common(o, &s, &top, T).unwrap().shape;
    let bottom = ogeom_algo::make_cone(o, up_at(0.0, 0.0, 0.0), 2.5, 13.5, 11.0, T)
        .unwrap()
        .shape;
    s = ogeom_bool::common(o, &s, &bottom, T).unwrap().shape;
    for side in [1.0_f64, -1.0] {
        // The straight edge's bevel, from the flat's contact to the wall's.
        let (p, q) = (
            Point::new(0.0, side * (big - d), 0.0),
            Point::new(wx, side * wy, 0.0),
        );
        let n = Vector::new(q.y - p.y, -(q.x - p.x), 0.0);
        let n = if n.dot(Vector::new(-1.0, side, 0.0)) < 0.0 {
            -n
        } else {
            n
        };
        s = cut_past(o, &s, p, n);
    }
    for (z, up) in [(10.0_f64, 1.0_f64), (0.0, -1.0)] {
        // The cap's straight edge's bevel.
        s = cut_past(o, &s, Point::new(d, 0.0, z), Vector::new(-1.0, 0.0, up));
        // The two corners' facets.
        for side in [1.0_f64, -1.0] {
            if !facets {
                break;
            }
            let on_cap = Point::new(d, side * (big - d).mul_add(big - d, -(d * d)).sqrt(), z);
            let on_flat = Point::new(0.0, side * (big - d), z - up * d);
            let on_wall = Point::new(wx, side * wy, z - up * d);
            let mut n = (on_flat - on_cap).cross(on_wall - on_cap);
            if n.dot(Vector::new(-1.0, side, up)) < 0.0 {
                n = -n;
            }
            s = cut_past(o, &s, on_cap, n);
        }
    }
    s
}

/// The D-shaped block with every edge chamfered 0.5 and mitred: at each of
/// its four corners on the flat the three bevels (two planes and a rim's
/// cone) run on to the point the three share and meet each other along
/// their sections from there, with no facet. The block is convex, so the
/// result is the block cut by each bevel's plane or cone alone: by
/// booleans, which share nothing with the engine, to the same volume, in
/// both orders. No corner generates a facet.
#[test]
fn a_d_block_chamfered_on_every_edge_mitres_each_corner_to_its_bevels_shared_point() {
    let d = 0.5_f64;
    let mut m = Model::new();
    let block = d_block(&mut m);
    let edges = explore_unique(&m, &block, ShapeType::Edge).unwrap();
    let vertices = explore_unique(&m, &block, ShapeType::Vertex).unwrap();
    let specs: Vec<(Shape, ogeom_fillet::Chamfer)> = edges
        .iter()
        .map(|e| (e.clone(), ogeom_fillet::Chamfer::Symmetric(d)))
        .collect();
    let mut o = Model::new();
    let s = d_block_cut_by_every_bevel(&mut o, d, false);
    let want = volume(&o, &s);
    let back: Vec<(Shape, ogeom_fillet::Chamfer)> = specs.iter().rev().cloned().collect();
    for asked in [&specs, &back] {
        let built = ogeom_fillet::chamfer_edges_cornered(
            &mut m,
            &block,
            asked,
            ogeom_fillet::ChamferCorner::Point,
            T,
        )
        .unwrap();
        let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        for v in &vertices {
            assert!(built.history.generated(v).is_empty(), "a facet");
        }
        let got = volume(&m, &built.shape);
        assert!(
            (got - want).abs() < 1e-6,
            "volume {got} against the booleans' {want}"
        );
    }
}

/// A disc of radius 5 and 2 thick, both rims chamfered 1.5: the bevels'
/// contacts cross on the wall, which goes; the two cones meet along their
/// ridge, the circle of radius 4.5 at mid height. The wall is deleted from
/// the history, each rim generates its one cone, and four faces are left:
/// the two caps and the two cones, the same solid to the bit in either
/// order. By Pappus the disc loses each rim's triangle of legs 1.5 round
/// its centroid, less their shared triangle `(4.5, 1), (5, 0.5), (5, 1.5)`
/// counted once.
#[test]
fn a_thin_disc_s_rim_bevels_consume_its_wall_and_meet_in_a_ridge() {
    use core::f64::consts::{PI, TAU};
    let mut m = Model::new();
    let disc = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 5.0, 2.0, T)
        .unwrap()
        .shape;
    let rims: Vec<Shape> = explore_unique(&m, &disc, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            matches!(
                ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap().0,
                ogeom_geom::Curve::Circle(_)
            )
        })
        .collect();
    assert_eq!(rims.len(), 2);
    let wall = explore_unique(&m, &disc, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| surface_kind(&m, f) == ogeom_geom::SurfaceKind::Cylinder)
        .unwrap();
    let built = ogeom_fillet::chamfer_edges(&mut m, &disc, &rims, 1.5, T).unwrap();
    let back: Vec<Shape> = rims.iter().rev().cloned().collect();
    let again = ogeom_fillet::chamfer_edges(&mut m, &disc, &back, 1.5, T).unwrap();
    assert!(ogeom_algo::check(&m, &built.shape, T).unwrap().is_valid());
    assert!(built.history.is_deleted(&wall), "the wall left");
    for e in &rims {
        let made = built.history.generated(e);
        assert_eq!(made.len(), 1, "a rim without its one bevel");
        assert_eq!(surface_kind(&m, &made[0]), ogeom_geom::SurfaceKind::Cone);
    }
    assert_eq!(
        explore_unique(&m, &built.shape, ShapeType::Face)
            .unwrap()
            .len(),
        4
    );
    let corners = |s: &Shape| {
        let mut out: Vec<[u64; 3]> = explore_unique(&m, s, ShapeType::Vertex)
            .unwrap()
            .iter()
            .map(|v| {
                let p = point_of(&m, v);
                [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
            })
            .collect();
        out.sort_unstable();
        out
    };
    assert!(
        corners(&built.shape) == corners(&again.shape),
        "the order changes the solid"
    );
    let lost = PI * 25.0 * 2.0 - volume(&m, &built.shape);
    let want = TAU * (2.0 * 1.125 * 4.5 - 0.25 * 29.0 / 6.0);
    assert!(
        (lost - want).abs() < 1e-6 * want,
        "lost {lost} against {want}"
    );
}

/// The D-shaped block with only its top cap's edges chamfered 0.5, the
/// straight one and the arc: at each end of the straight edge the plane
/// bevel and the arc's cone meet across the upright edge left sharp, along
/// the conic the plane cuts from the cone, with no facet. Valid, the same
/// solid in either order, each edge generating its bevel and the corners
/// on the flat none. The oracle: the block inside the arc's cone and under
/// the straight edge's bevel plane, by booleans.
#[test]
fn a_d_block_s_cap_chamfered_meets_its_bevels_across_the_sharp_uprights() {
    let d = 0.5_f64;
    let mut m = Model::new();
    let block = d_block(&mut m);
    let top: Vec<Shape> = explore_unique(&m, &block, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let b = ogeom_algo::shape_bounds(&m, e, T).unwrap();
            b.low().unwrap().z > 10.0 - 1e-6
        })
        .collect();
    assert!(top.len() >= 2, "the cap's straight edge and its arc");
    let built = ogeom_fillet::chamfer_edges(&mut m, &block, &top, d, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    let back: Vec<Shape> = top.iter().rev().cloned().collect();
    let again = ogeom_fillet::chamfer_edges(&mut m, &block, &back, d, T).unwrap();
    assert!(
        (volume(&m, &built.shape) - volume(&m, &again.shape)).abs() < 1e-9,
        "the order changes the solid"
    );
    for e in &top {
        assert_eq!(
            built.history.generated(e).len(),
            1,
            "an edge without its bevel"
        );
    }
    for v in explore_unique(&m, &block, ShapeType::Vertex).unwrap() {
        if point_of(&m, &v).x.abs() < 1e-9 && point_of(&m, &v).z > 10.0 - 1e-6 {
            assert!(
                built.history.generated(&v).is_empty(),
                "a facet at a corner of two"
            );
        }
    }
    let mut o = Model::new();
    let s = d_block_cap_bevelled(&mut o, d);
    let (got, want) = (volume(&m, &built.shape), volume(&o, &s));
    assert!(
        (got - want).abs() < 1e-6,
        "volume {got} against the booleans' {want}"
    );
}

/// The D-shaped block with its cap's arc and straight edge bevelled at
/// `d` by booleans: in common with the arc's 45-degree cone, and cut past
/// the straight edge's bevel plane.
fn d_block_cap_bevelled(o: &mut Model, d: f64) -> Shape {
    let s = d_block(o);
    let cone = ogeom_algo::make_cone(o, up_at(0.0, 0.0, -1.0), 14.0 - d, 3.0 - d, 11.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::common(o, &s, &cone, T).unwrap().shape;
    cut_past(o, &s, Point::new(d, 0.0, 10.0), Vector::new(-1.0, 0.0, 1.0))
}

/// The D-shaped block with its arc wall an extrusion of the half circle as
/// an exact rational B-spline, not a cylinder: the cap's arc is a crease
/// with no closed-form bevel, and the straight edge's bevel plane meets its
/// bevel across the sharp uprights.
fn spline_d_block(m: &mut Model) -> Shape {
    use ogeom_geom::{CircleCurve, Curve};
    use ogeom_math::{Axis, Circle, Plane};
    let half = core::f64::consts::FRAC_PI_2;
    let circle: Curve = CircleCurve::new(Circle::new(up_at(0.0, 0.0, 0.0), 3.0, T).unwrap()).into();
    let arc: Curve = circle.to_bspline_over((-half, half), T).unwrap().into();
    let (low, high) = (Point::new(0.0, -3.0, 0.0), Point::new(0.0, 3.0, 0.0));
    let (a, b) = (
        ogeom_algo::make_vertex(m, low).shape,
        ogeom_algo::make_vertex(m, high).shape,
    );
    let range = ogeom_geom::Curve3d::domain(&arc);
    let curved = ogeom_algo::make_edge_between(m, arc, range, &a, &b, T)
        .unwrap()
        .shape;
    let line: Curve =
        ogeom_geom::LineCurve::new(Axis::new(high, Direction::new(low - high, T).unwrap())).into();
    let straight = ogeom_algo::make_edge_between(m, line, (0.0, 6.0), &b, &a, T)
        .unwrap()
        .shape;
    let wire = ogeom_algo::make_wire(m, &[curved, straight], T)
        .unwrap()
        .shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(m, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    ogeom_algo::make_prism(m, &face, Vector::new(0.0, 0.0, 10.0), T)
        .unwrap()
        .shape
}

/// The spline D-block's cap chamfered 0.5 on its arc and straight edge:
/// the ruled bevel along the spline crease and the straight edge's plane
/// meet across the sharp uprights, on the curve the plane cuts from the
/// ruled bevel, marched. The wall is a circle's extrusion exactly, so the
/// solid is the cylinder's D-block bevelled by booleans, to within the
/// bevel's fit.
#[test]
fn a_spline_crease_s_bevel_meets_a_plane_bevel_across_the_sharp_uprights() {
    let d = 0.5_f64;
    let mut m = Model::new();
    let block = spline_d_block(&mut m);
    let found = ogeom_algo::check(&m, &block, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    let top: Vec<Shape> = explore_unique(&m, &block, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let b = ogeom_algo::shape_bounds(&m, e, T).unwrap();
            b.low().unwrap().z > 10.0 - 1e-6
        })
        .collect();
    assert_eq!(top.len(), 2, "the cap's arc and straight edge");
    let built = ogeom_fillet::chamfer_edges(&mut m, &block, &top, d, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    for e in &top {
        assert_eq!(
            built.history.generated(e).len(),
            1,
            "an edge without its bevel"
        );
    }
    let mut o = Model::new();
    let s = d_block_cap_bevelled(&mut o, d);
    let (got, want) = (volume(&m, &built.shape), volume(&o, &s));
    assert!(
        (got - want).abs() < 1e-7,
        "volume {got} against the booleans' {want}"
    );
}

/// The plane of a bevel along a convex straight edge through `at` along
/// the unit `u`, between faces of outward normals `na` and `nb`, set back
/// `da` and `db` square to the edge in each: an outward normal and its
/// offset.
fn bevel_plane(
    at: Point,
    u: Vector,
    (na, da): (Vector, f64),
    (nb, db): (Vector, f64),
) -> (Vector, f64) {
    let inward = |n: Vector, other: Vector| {
        let m = n.cross(u);
        let m = m / m.magnitude();
        if m.dot(other) > 0.0 { -m } else { m }
    };
    let (ca, cb) = (at + inward(na, nb) * da, at + inward(nb, na) * db);
    let mut n = u.cross(cb - ca);
    n = n / n.magnitude();
    if n.dot(na + nb) < 0.0 {
        n = -n;
    }
    (n, n.dot(ca.to_vector()))
}

/// The disc of radius 3 and 10 high cut to `x <= 1`, `y >= 0`, and by the
/// plane through the cap's corner `v = (1, √8, 10)` tilted to take off the
/// cap's edge on `x = 1` more the farther from `v`: four faces meet at `v`,
/// the cap, the wall, the flat `x = 1` and the tilted plane, along the cap's
/// rim, the wall's ruling and two lines. Returns the solid, `v` and the
/// tilted plane's outward normal.
fn tilted_corner(m: &mut Model) -> (Shape, Point, Vector) {
    let disc = ogeom_algo::make_cylinder(m, up_at(0.0, 0.0, 0.0), 3.0, 10.0, T)
        .unwrap()
        .shape;
    let s = cut_past(m, &disc, Point::new(1.0, 0.0, 0.0), Vector::X);
    let s = cut_past(m, &s, Point::ORIGIN, -Vector::Y);
    let v = Point::new(1.0, 8.0_f64.sqrt(), 10.0);
    let (a, b) = (
        Point::new(1.0, v.y - 4.0, 9.0),
        Point::new(0.0, v.y - 4.0, 10.0),
    );
    let mut n = (a - v).cross(b - v);
    if n.dot(Vector::new(1.0, 0.0, 1.0)) < 0.0 {
        n = -n;
    }
    let n = n / n.magnitude();
    (cut_past(m, &s, v, n), v, n)
}

/// [`tilted_corner`]'s four edges at `v` chamfered 0.3: the rim's bevel is
/// a cone, the other three planes, and their four corners at `v` lie off
/// one plane, so the corner closes on the patch through them, the cone's
/// side the curve it meets the plane through its two corners and the
/// hole's middle. Valid, `v` generating the patch's four sector faces, all
/// B-splines; every far end flush on the cut `y = 0` or the floor. The
/// patch takes off no less than the bevels' mitre, by booleans the solid
/// inside the rim's cone and under the three bevel planes, and no more
/// than that mitre's peak past a plane `3 d` in from `v`.
#[test]
fn a_corner_of_four_bevels_one_a_cone_closes_on_one_patch_off_a_plane() {
    let d = 0.3;
    let mut m = Model::new();
    let (s, v, tilt) = tilted_corner(&mut m);
    let at_v = |m: &Model, e: &Shape| {
        let (p, q) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
        point_of(m, &p).distance(v) < 1e-6 || point_of(m, &q).distance(v) < 1e-6
    };
    let edges: Vec<Shape> = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| at_v(&m, e))
        .collect();
    assert_eq!(edges.len(), 4, "four edges at the corner");
    let corner = explore_unique(&m, &s, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .filter(|x| point_of(&m, x).distance(v) < 1e-6)
        .collect::<Vec<_>>();
    let built = ogeom_fillet::chamfer_edges(&mut m, &s, &edges, d, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    let patch: Vec<Shape> = corner
        .iter()
        .flat_map(|x| built.history.generated(x).to_vec())
        .collect();
    assert_eq!(patch.len(), 4, "the patch's four sectors");
    for f in &patch {
        assert_eq!(surface_kind(&m, f), ogeom_geom::SurfaceKind::BSpline);
    }
    let got = volume(&m, &built.shape);

    // The mitre, by booleans.
    let mut o = Model::new();
    let (mut mitre, _, _) = tilted_corner(&mut o);
    let r = |z: f64| (3.0 - d) + (10.0 - z);
    let cone = ogeom_algo::make_cone(&mut o, up_at(0.0, 0.0, -1.0), r(-1.0), r(10.5), 11.5, T)
        .unwrap()
        .shape;
    mitre = ogeom_bool::common(&mut o, &mitre, &cone, T).unwrap().shape;
    let up = Vector::Z;
    let flat = Vector::X;
    let unit = |w: Vector| w / w.magnitude();
    // The ruling's bevel: upright, through the flat's contact `d` along it
    // and the wall's `d` round it.
    let phi = v.y.atan2(v.x) + d / 3.0;
    let (on_flat, on_wall) = (
        Point::new(1.0, v.y - d, 0.0),
        Point::new(3.0 * phi.cos(), 3.0 * phi.sin(), 0.0),
    );
    let mut n = (on_wall - on_flat).cross(up);
    if n.dot(Vector::new(1.0, 1.0, 0.0)) < 0.0 {
        n = -n;
    }
    mitre = cut_past(&mut o, &mitre, on_flat, n);
    for (u, na, nb) in [
        (unit(Vector::new(-1.0, -4.0, 0.0)), up, tilt),
        (unit(Vector::new(0.0, -4.0, -1.0)), flat, tilt),
    ] {
        let (n, off) = bevel_plane(v, u, (na, d), (nb, d));
        mitre = cut_past(&mut o, &mitre, Point::ORIGIN + n * off, n);
    }
    let most = volume(&o, &mitre);
    let out = unit(up + flat + tilt + Vector::new(v.x, v.y, 0.0) / 3.0);
    let peak = cut_past(&mut o, &mitre, v - out * (3.0 * d), -out);
    let least = most - volume(&o, &peak);
    assert!(
        least < got && got < most,
        "volume {got} outside ({least}, {most})"
    );
}

/// Two upright drums of radii 3 and 2, their axes 3.5 apart, 10 tall,
/// fused or, with `fused` false, their common lens; each drum cut by the
/// cone of its own cap rim's bevel at `d` first where `bevelled`.
fn twin_drums(m: &mut Model, fused: bool, bevelled: Option<f64>) -> Shape {
    let drum = |m: &mut Model, x: f64, r: f64| {
        let s = ogeom_algo::make_cylinder(m, up_at(x, 0.0, 0.0), r, 10.0, T)
            .unwrap()
            .shape;
        match bevelled {
            None => s,
            Some(d) => {
                // At 45 degrees: radius r - d at the cap, r at 10 - d.
                let low = r - d + 11.0;
                let cone = ogeom_algo::make_cone(m, up_at(x, 0.0, -1.0), low, r - d, 11.0, T)
                    .unwrap()
                    .shape;
                ogeom_bool::common(m, &s, &cone, T).unwrap().shape
            }
        }
    };
    let (a, b) = (drum(m, 0.0, 3.0), drum(m, 3.5, 2.0));
    let joined = if fused {
        ogeom_bool::fuse(m, &a, &b, T).unwrap().shape
    } else {
        ogeom_bool::common(m, &a, &b, T).unwrap().shape
    };
    ogeom_heal::unify_same_domain(m, &joined, T)
        .unwrap()
        .0
        .shape
}

/// The twin drums' cap rims chamfered 0.5, in both orders: valid, the same
/// solid either way, each arc's one bevel generated, and the volume of the
/// drums each cut by its own bevel's cone first, by booleans, which share
/// nothing with the engine.
fn twin_rims_chamfer_as_each_drum_alone(fused: bool) {
    let d = 0.5_f64;
    let mut m = Model::new();
    let block = twin_drums(&mut m, fused, None);
    let top: Vec<Shape> = explore_unique(&m, &block, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let b = ogeom_algo::shape_bounds(&m, e, T).unwrap();
            b.low().unwrap().z > 10.0 - 1e-6
        })
        .collect();
    assert!(top.len() >= 2, "the cap's two arcs");
    let built = ogeom_fillet::chamfer_edges(&mut m, &block, &top, d, T).unwrap();
    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    let back: Vec<Shape> = top.iter().rev().cloned().collect();
    let again = ogeom_fillet::chamfer_edges(&mut m, &block, &back, d, T).unwrap();
    assert!(
        (volume(&m, &built.shape) - volume(&m, &again.shape)).abs() < 1e-9,
        "the order changes the solid"
    );
    for e in &top {
        assert_eq!(
            built.history.generated(e).len(),
            1,
            "an edge without its bevel"
        );
    }
    let mut o = Model::new();
    let s = twin_drums(&mut o, fused, Some(d));
    let (got, want) = (volume(&m, &built.shape), volume(&o, &s));
    assert!(
        (got - want).abs() < 1e-6,
        "volume {got} against the booleans' {want}"
    );
}

/// Two fused drums of different radii, the cap's two arcs chamfered: the
/// two cone bevels meet across the concave upright where the walls cross,
/// and no plane mirrors one into the other, so they meet on the curve the
/// two cones share, marched. At a concave edge the material is the union
/// of what each side leaves: each drum cut by its own cone, fused.
#[test]
fn two_unequal_rim_bevels_meet_across_a_concave_upright_on_their_cones_curve() {
    twin_rims_chamfer_as_each_drum_alone(true);
}

/// The two drums' common lens, its cap's two arcs chamfered: the cone
/// bevels meet across the convex uprights on their shared curve. At a
/// convex edge the material is what both sides leave: each drum cut by its
/// own cone, in common.
#[test]
fn two_unequal_rim_bevels_meet_across_a_convex_upright_on_their_cones_curve() {
    twin_rims_chamfer_as_each_drum_alone(false);
}

/// A web of an arch: the annular sector between radii 10 and 11 about z,
/// from -15 to 15 degrees, 1 tall.
fn arch_web(m: &mut Model) -> Shape {
    use ogeom_geom::{CircleCurve, Curve, LineCurve};
    use ogeom_math::{Axis, Circle, Plane};
    let w = 15.0_f64.to_radians();
    let on = |r: f64, a: f64| Point::new(r * a.cos(), r * a.sin(), 0.0);
    let vertex = |m: &mut Model, p: Point| ogeom_algo::make_vertex(m, p).shape;
    let (oa, ob, ib, ia) = (
        vertex(m, on(11.0, -w)),
        vertex(m, on(11.0, w)),
        vertex(m, on(10.0, w)),
        vertex(m, on(10.0, -w)),
    );
    let arc = |m: &mut Model, r: f64, z: Direction, from: f64, a: &Shape, b: &Shape| {
        let x = Direction::new(Vector::new(from.cos(), from.sin(), 0.0), T).unwrap();
        let frame = Frame::new(Point::ORIGIN, z, x, T).unwrap();
        let curve: Curve = CircleCurve::new(Circle::new(frame, r, T).unwrap()).into();
        ogeom_algo::make_edge_between(m, curve, (0.0, 2.0 * w), a, b, T)
            .unwrap()
            .shape
    };
    let line = |m: &mut Model, p: Point, q: Point, a: &Shape, b: &Shape| {
        let curve: Curve = LineCurve::new(Axis::new(p, Direction::new(q - p, T).unwrap())).into();
        ogeom_algo::make_edge_between(m, curve, (0.0, p.distance(q)), a, b, T)
            .unwrap()
            .shape
    };
    let down = Direction::new(Vector::new(0.0, 0.0, -1.0), T).unwrap();
    let edges = [
        arc(m, 11.0, Direction::Z, -w, &oa, &ob),
        line(m, on(11.0, w), on(10.0, w), &ob, &ib),
        arc(m, 10.0, down, w, &ib, &ia),
        line(m, on(10.0, -w), on(11.0, -w), &ia, &oa),
    ];
    let wire = ogeom_algo::make_wire(m, &edges, T).unwrap().shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(m, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    ogeom_algo::make_prism(m, &face, Vector::new(0.0, 0.0, 1.0), T)
        .unwrap()
        .shape
}

/// The arch web's two outer uprights chamfered 0.5 on the end faces and
/// 3.2 round the outer wall: the two contacts on the wall cross (its arc is
/// about 2.88 a side), so the wall goes and the two bevel planes meet in
/// their ridge, an upright line midway. Valid, the wall deleted, the same
/// in both orders. Every face but the ends stands square to the ends, so
/// the volume is the section's area: the polygon from the inner arc's ends
/// to the bevels' feet on the ends and the ridge, less the inner arc's
/// segment.
#[test]
fn an_arch_web_s_outer_wall_consumed_by_its_uprights_bevels_meets_in_a_ridge() {
    let (on_end, on_wall) = (0.5_f64, 3.2_f64);
    let w = 15.0_f64.to_radians();
    let mut m = Model::new();
    let web = arch_web(&mut m);
    let found = ogeom_algo::check(&m, &web, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    let uprights: Vec<Shape> = explore_unique(&m, &web, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let b = ogeom_algo::shape_bounds(&m, e, T).unwrap();
            let (p, q) = (b.low().unwrap(), b.high().unwrap());
            q.z - p.z > 0.5 && q.x.hypot(q.y.abs().max(p.y.abs())) > 10.99
        })
        .collect();
    assert_eq!(uprights.len(), 2, "the outer uprights");
    let wall = explore_unique(&m, &web, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            surface_kind(&m, f) == ogeom_geom::SurfaceKind::Cylinder && {
                let b = ogeom_algo::shape_bounds(&m, f, T).unwrap();
                b.high().unwrap().x > 10.9
            }
        })
        .unwrap();
    let specs: Vec<(Shape, ogeom_fillet::Chamfer)> = uprights
        .iter()
        .map(|e| {
            let end = explore_unique(&m, &web, ShapeType::Face)
                .unwrap()
                .into_iter()
                .find(|f| {
                    surface_kind(&m, f) == ogeom_geom::SurfaceKind::Plane
                        && explore_unique(&m, f, ShapeType::Edge)
                            .unwrap()
                            .iter()
                            .any(|x| x.is_same(e))
                })
                .unwrap();
            (
                e.clone(),
                ogeom_fillet::Chamfer::Distances {
                    face: end,
                    on_face: on_end,
                    on_other: on_wall,
                },
            )
        })
        .collect();
    // The section: bevel A from its foot on the end at -15 degrees to its
    // contact on the wall, meeting its mirror at 0 degrees.
    let polar = |r: f64, a: f64| (r * a.cos(), r * a.sin());
    let foot = polar(11.0 - on_end, -w);
    let contact = polar(11.0, -w + on_wall / 11.0);
    let x = foot.0 - foot.1 * (contact.0 - foot.0) / (contact.1 - foot.1);
    assert!(x > 10.0 && x < 11.0, "the ridge stands in the web: {x}");
    let inner = polar(10.0, -w);
    let ring = [
        inner,
        foot,
        (x, 0.0),
        (foot.0, -foot.1),
        (inner.0, -inner.1),
    ];
    let mut twice = 0.0;
    for k in 0..ring.len() {
        let (p, q) = (ring[k], ring[(k + 1) % ring.len()]);
        twice += p.0 * q.1 - q.0 * p.1;
    }
    let segment = 50.0 * (2.0 * w - (2.0 * w).sin());
    let want = twice.abs() / 2.0 - segment;
    let back: Vec<(Shape, ogeom_fillet::Chamfer)> = specs.iter().rev().cloned().collect();
    for asked in [&specs, &back] {
        let built = ogeom_fillet::chamfer_edges_with(&mut m, &web, asked, T).unwrap();
        let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        assert!(built.history.is_deleted(&wall), "the wall stays");
        let got = volume(&m, &built.shape);
        assert!(
            (got - want).abs() < 1e-6,
            "volume {got} against the section's {want}"
        );
    }
}
