//! Asymmetric rounds: a different distance onto each face of a straight
//! edge between two planes, the section an ellipse or a conic.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2};

use ogeom_core::Tolerances;
use ogeom_fillet::{Blend, ChamferCorner, RoundProfile};
use ogeom_math::{Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn fine() -> ogeom_mesh::Deflection {
    ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    }
}

/// A vertex's point, placed.
fn point_of(m: &Model, v: &Shape) -> Point {
    let p = m
        .node(v)
        .and_then(|n| n.data().as_vertex().map(|d| d.point))
        .unwrap();
    v.transform(m.datums()).unwrap().apply(p)
}

/// The edge of `s` whose two ends satisfy `on`.
fn edge_where(m: &Model, s: &Shape, on: impl Fn(Point) -> bool) -> Shape {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(m, e)
                .unwrap()
                .is_some_and(|(a, b)| on(point_of(m, &a)) && on(point_of(m, &b)))
        })
        .expect("the edge")
}

/// The face of `s` holding `p` inside it.
fn face_at(m: &Model, s: &Shape, p: Point) -> Shape {
    explore_unique(m, s, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            ogeom_algo::classify_on_face(m, f, p, fine(), T)
                .is_ok_and(|c| c == ogeom_algo::Containment::In)
        })
        .expect("the face")
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, fine(), T).unwrap().mass
}

/// Every vertex of `s`, as bits, sorted: two solids built alike to the
/// bit list the same.
fn vertex_bits(m: &Model, s: &Shape) -> Vec<[u64; 3]> {
    let mut out: Vec<[u64; 3]> = explore_unique(m, s, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| {
            let p = point_of(m, v);
            [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
        })
        .collect();
    out.sort_unstable();
    out
}

/// The elliptical profile's weight on a dihedral of angle `alpha`.
fn elliptical(alpha: f64) -> f64 {
    (alpha * 0.5).sin()
}

/// The conic profile's weight for `rho`.
fn conic(rho: f64) -> f64 {
    rho / (1.0 - rho)
}

/// The area a conic section of weight `w` takes off (or adds to) the
/// corner between legs `d1` and `d2` meeting at `alpha`: the control
/// triangle less the segment between the arc and its chord. The segment's
/// share of the triangle is affine-invariant, read on the unit circle
/// (`w < 1`), the parabola (`w = 1`) or the unit hyperbola (`w > 1`).
fn corner_area(d1: f64, d2: f64, alpha: f64, w: f64) -> f64 {
    let triangle = 0.5 * d1 * d2 * alpha.sin();
    let share = if (w - 1.0).abs() < 1e-12 {
        2.0 / 3.0
    } else if w < 1.0 {
        let t = w.acos();
        (t - t.sin() * t.cos()) * t.cos() / t.sin().powi(3)
    } else {
        let p = w.acosh();
        (p.sinh() * p.cosh() - p) * p.cosh() / p.sinh().powi(3)
    };
    triangle * (1.0 - share)
}

/// The band's contact lines (its straight edges), each with the worst
/// angle between the band's and the face's normals along it.
fn contacts(m: &Model, s: &Shape) -> Vec<(Shape, f64)> {
    let mut out = Vec::new();
    for f in explore_unique(m, s, ShapeType::Face).unwrap() {
        let data = m.node(&f).unwrap().data().as_face().unwrap();
        if !matches!(
            m.geometry().surface(data.surface),
            Some(ogeom_geom::SurfaceGeometry::Extrusion(_))
        ) {
            continue;
        }
        for c in ogeom_fillet::analyse_blend(m, s, &f, 17, T).unwrap() {
            let (curve, _) = ogeom_algo::placed_edge_curve(m, &c.edge, T)
                .unwrap()
                .unwrap();
            if matches!(curve, ogeom_geom::Curve::Line(_)) {
                out.push((c.edge, c.tangency_error));
            }
        }
    }
    out
}

/// Checks the result is valid, holds `exact` within `within`, and meets
/// both faces tangent along its two contact lines, each through a point of
/// `lines`; its volume.
fn holds(m: &Model, s: &Shape, exact: f64, within: f64, lines: [(Point, Vector); 2]) -> f64 {
    let diagnosis = ogeom_algo::check(m, s, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let v = volume(m, s);
    assert!(
        (v - exact).abs() <= within,
        "volume {v} against the exact {exact}"
    );
    let found = contacts(m, s);
    assert_eq!(found.len(), 2, "the band meets two faces along lines");
    for (edge, tangency) in &found {
        assert!(
            tangency.abs() <= 1e-9,
            "the band meets a face at {tangency} rad"
        );
        let (a, b) = ogeom_algo::edge_vertices(m, edge).unwrap().unwrap();
        let on = |(through, along): (Point, Vector)| {
            [point_of(m, &a), point_of(m, &b)].iter().all(|p| {
                let off = *p - through;
                (off - along * off.dot(along)).magnitude() <= 1e-9
            })
        };
        assert!(
            lines.iter().any(|l| on(*l)),
            "a contact line stands where no distance puts it"
        );
    }
    v
}

/// A 20 x 20 x 10 box, its top front edge along x, and the top and front
/// faces.
fn block(m: &mut Model) -> (Shape, Shape, Shape, Shape) {
    let s = ogeom_algo::make_box(m, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edge = edge_where(m, &s, |p| p.y.abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9);
    let top = face_at(m, &s, Point::new(10.0, 10.0, 10.0));
    let front = face_at(m, &s, Point::new(10.0, 0.0, 5.0));
    (s, edge, top, front)
}

const X: Vector = Vector::new(1.0, 0.0, 0.0);
const Y: Vector = Vector::new(0.0, 1.0, 0.0);

#[test]
fn an_elliptical_round_meets_each_face_at_its_own_distance_whichever_is_named() {
    let mut m = Model::new();
    let (s, edge, top, front) = block(&mut m);
    let (d_top, d_front) = (3.0, 1.5);
    let profile = RoundProfile::Elliptical;
    let a =
        ogeom_fillet::fillet_edge_distances(&mut m, &s, &edge, &top, d_top, d_front, profile, T)
            .unwrap();
    let b =
        ogeom_fillet::fillet_edge_distances(&mut m, &s, &edge, &front, d_front, d_top, profile, T)
            .unwrap();
    let exact = 4000.0 - 20.0 * corner_area(d_top, d_front, FRAC_PI_2, elliptical(FRAC_PI_2));
    let lines = [
        (Point::new(0.0, d_top, 10.0), X),
        (Point::new(0.0, 0.0, 10.0 - d_front), X),
    ];
    let va = holds(&m, &a.shape, exact, 1e-8, lines);
    let vb = holds(&m, &b.shape, exact, 1e-8, lines);
    assert_eq!(va.to_bits(), vb.to_bits(), "naming the other face first");
    assert_eq!(vertex_bits(&m, &a.shape), vertex_bits(&m, &b.shape));
    assert!(a.history.is_deleted(&edge));
}

#[test]
fn a_conic_round_takes_off_its_conic_segment_for_every_rho() {
    for rho in [0.1, 0.3, 0.5, 0.7, 0.9] {
        let mut m = Model::new();
        let (s, edge, top, front) = block(&mut m);
        let (d_top, d_front) = (2.0, 4.0);
        let profile = RoundProfile::Conic(rho);
        let a = ogeom_fillet::fillet_edge_distances(
            &mut m, &s, &edge, &top, d_top, d_front, profile, T,
        )
        .unwrap();
        let b = ogeom_fillet::fillet_edge_distances(
            &mut m, &s, &edge, &front, d_front, d_top, profile, T,
        )
        .unwrap();
        let exact = 4000.0 - 20.0 * corner_area(d_top, d_front, FRAC_PI_2, conic(rho));
        let lines = [
            (Point::new(0.0, d_top, 10.0), X),
            (Point::new(0.0, 0.0, 10.0 - d_front), X),
        ];
        let va = holds(&m, &a.shape, exact, 1e-8, lines);
        let vb = holds(&m, &b.shape, exact, 1e-8, lines);
        assert_eq!(
            va.to_bits(),
            vb.to_bits(),
            "rho {rho}, naming the other face"
        );
    }
}

/// A prism along y whose top face meets a slanted face at a convex edge of
/// 150 degrees through (8, y, 10); the solid, the edge, the top and the
/// slanted face, the polygon's area.
fn slanted(m: &mut Model) -> (Shape, Shape, Shape, Shape, f64) {
    let rise = 8.0 * (30.0_f64).to_radians().tan();
    let points = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(0.0, 0.0, 10.0 - rise),
        Point::new(8.0, 0.0, 10.0),
        Point::new(20.0, 0.0, 10.0),
        Point::new(20.0, 0.0, 0.0),
    ];
    let wire = ogeom_algo::make_polygon(m, &points, true, T).unwrap().shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Y),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(m, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let s = ogeom_algo::make_prism(m, &face, Vector::new(0.0, 20.0, 0.0), T)
        .unwrap()
        .shape;
    let edge = edge_where(m, &s, |p| {
        (p.x - 8.0).abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9
    });
    let top = face_at(m, &s, Point::new(14.0, 10.0, 10.0));
    let slant = face_at(m, &s, Point::new(4.0, 10.0, 10.0 - rise * 0.5));
    let area = 20.0 * 10.0 - 0.5 * 8.0 * rise;
    (s, edge, top, slant, area)
}

#[test]
fn a_slanted_edge_takes_the_ellipse_of_its_own_dihedral() {
    let mut m = Model::new();
    let (s, edge, top, slant, area) = slanted(&mut m);
    let (d_top, d_slant) = (3.0, 2.0);
    let alpha = (150.0_f64).to_radians();
    let profile = RoundProfile::Elliptical;
    let a =
        ogeom_fillet::fillet_edge_distances(&mut m, &s, &edge, &top, d_top, d_slant, profile, T)
            .unwrap();
    let b =
        ogeom_fillet::fillet_edge_distances(&mut m, &s, &edge, &slant, d_slant, d_top, profile, T)
            .unwrap();
    let exact = 20.0 * (area - corner_area(d_top, d_slant, alpha, elliptical(alpha)));
    let down = Vector::new(
        -(30.0_f64).to_radians().cos(),
        0.0,
        -(30.0_f64).to_radians().sin(),
    );
    let lines = [
        (Point::new(8.0 + d_top, 0.0, 10.0), Y),
        (Point::new(8.0, 0.0, 10.0) + down * d_slant, Y),
    ];
    let va = holds(&m, &a.shape, exact, 1e-8, lines);
    let vb = holds(&m, &b.shape, exact, 1e-8, lines);
    assert_eq!(va.to_bits(), vb.to_bits(), "naming the other face first");
}

#[test]
fn equal_elliptical_distances_are_the_circular_round_to_the_bit() {
    // On the square edge the ball's radius is the distance; on the slanted
    // one it is the distance times tan(alpha / 2).
    let mut m = Model::new();
    let (s, edge, top, _) = block(&mut m);
    let asked = ogeom_fillet::fillet_edge_distances(
        &mut m,
        &s,
        &edge,
        &top,
        2.5,
        2.5,
        RoundProfile::Elliptical,
        T,
    )
    .unwrap();
    let round = ogeom_fillet::fillet_edge(&mut m, &s, &edge, 2.5, T).unwrap();
    assert_eq!(
        volume(&m, &asked.shape).to_bits(),
        volume(&m, &round.shape).to_bits()
    );
    assert_eq!(vertex_bits(&m, &asked.shape), vertex_bits(&m, &round.shape));

    let mut m = Model::new();
    let (s, edge, top, slant, _) = slanted(&mut m);
    let alpha = (150.0_f64).to_radians();
    let radius = 2.0 * alpha.sin() / (1.0 + alpha.cos());
    let asked = ogeom_fillet::fillet_edge_distances(
        &mut m,
        &s,
        &edge,
        &slant,
        2.0,
        2.0,
        RoundProfile::Elliptical,
        T,
    )
    .unwrap();
    let round = ogeom_fillet::fillet_edge(&mut m, &s, &edge, radius, T).unwrap();
    assert_eq!(
        volume(&m, &asked.shape).to_bits(),
        volume(&m, &round.shape).to_bits()
    );
    assert_eq!(vertex_bits(&m, &asked.shape), vertex_bits(&m, &round.shape));
    let _ = top;
}

#[test]
fn the_conic_band_of_the_circle_s_own_weight_is_the_circular_round() {
    // Not the same call: the conic band itself, at the weight the circle of
    // a square dihedral takes, against the cylinder.
    let rho = FRAC_1_SQRT_2 / (1.0 + FRAC_1_SQRT_2);
    let mut m = Model::new();
    let (s, edge, top, _) = block(&mut m);
    let conic = ogeom_fillet::fillet_edge_distances(
        &mut m,
        &s,
        &edge,
        &top,
        2.5,
        2.5,
        RoundProfile::Conic(rho),
        T,
    )
    .unwrap();
    let round = ogeom_fillet::fillet_edge(&mut m, &s, &edge, 2.5, T).unwrap();
    let (vc, vr) = (volume(&m, &conic.shape), volume(&m, &round.shape));
    assert!(
        (vc - vr).abs() <= 1e-9 * vr,
        "conic {vc} against the round {vr}"
    );
    // And through the round's arc middle, on the ball's diagonal.
    let r = 2.5;
    let middle = Point::new(10.0, r - r * FRAC_1_SQRT_2, 10.0 - r + r * FRAC_1_SQRT_2);
    let on_band = explore_unique(&m, &conic.shape, ShapeType::Face)
        .unwrap()
        .iter()
        .any(|f| {
            ogeom_algo::classify_on_face(&m, f, middle, fine(), T)
                .is_ok_and(|c| c != ogeom_algo::Containment::Out)
        });
    assert!(on_band, "the conic band passes through the circle's middle");
}

/// An L: a 20 x 20 x 20 block with the corner above z = 8 and past x = 8
/// taken out; its re-entrant edge along y through (8, y, 8), the floor and
/// the wall beside it.
fn l_bracket(m: &mut Model) -> (Shape, Shape, Shape, Shape) {
    let block = ogeom_algo::make_box(m, Frame::WORLD, (20.0, 20.0, 20.0), T)
        .unwrap()
        .shape;
    let corner = Frame::new(Point::new(8.0, -1.0, 8.0), Direction::Z, Direction::X, T).unwrap();
    let notch = ogeom_algo::make_box(m, corner, (20.0, 22.0, 20.0), T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(m, &block, &notch, T).unwrap().shape;
    let edge = edge_where(m, &s, |p| {
        (p.x - 8.0).abs() < 1e-9 && (p.z - 8.0).abs() < 1e-9
    });
    let floor = face_at(m, &s, Point::new(14.0, 10.0, 8.0));
    let wall = face_at(m, &s, Point::new(8.0, 10.0, 14.0));
    (s, edge, floor, wall)
}

#[test]
fn a_concave_edge_fills_with_its_conic_segment_and_nothing_past_the_faces() {
    for profile in [RoundProfile::Elliptical, RoundProfile::Conic(0.6)] {
        let mut m = Model::new();
        let (s, edge, floor, wall) = l_bracket(&mut m);
        let (d_floor, d_wall) = (3.0, 1.5);
        let a = ogeom_fillet::fillet_edge_distances(
            &mut m, &s, &edge, &floor, d_floor, d_wall, profile, T,
        )
        .unwrap();
        let b = ogeom_fillet::fillet_edge_distances(
            &mut m, &s, &edge, &wall, d_wall, d_floor, profile, T,
        )
        .unwrap();
        let w = match profile {
            RoundProfile::Elliptical => elliptical(FRAC_PI_2),
            RoundProfile::Conic(rho) => conic(rho),
        };
        let l_volume = 8000.0 - 12.0 * 12.0 * 20.0;
        let exact = l_volume + 20.0 * corner_area(d_floor, d_wall, FRAC_PI_2, w);
        let lines = [
            (Point::new(8.0 + d_floor, 0.0, 8.0), Y),
            (Point::new(8.0, 0.0, 8.0 + d_wall), Y),
        ];
        let va = holds(&m, &a.shape, exact, 1e-8, lines);
        let vb = holds(&m, &b.shape, exact, 1e-8, lines);
        assert_eq!(va.to_bits(), vb.to_bits(), "naming the other face first");
        // The fill stands within the corner's control triangle: inside the
        // solid's box, and no point of the band past either contact.
        let bounds = ogeom_algo::shape_bounds(&m, &a.shape, T).unwrap();
        let before = ogeom_algo::shape_bounds(&m, &s, T).unwrap();
        assert!(
            (bounds.low().unwrap() - before.low().unwrap()).magnitude() <= 1e-7
                && (bounds.high().unwrap() - before.high().unwrap()).magnitude() <= 1e-7,
            "the fill overhangs the bracket"
        );
        for v in explore_unique(&m, &a.shape, ShapeType::Vertex).unwrap() {
            let p = point_of(&m, &v);
            let (x, z) = (p.x - 8.0, p.z - 8.0);
            let fill = x > 1e-9 && z > 1e-9;
            assert!(
                !fill || x / d_floor + z / d_wall <= 1.0 + 1e-9,
                "a vertex of the fill at {p:?} stands past the chord"
            );
        }
    }
}

#[test]
fn what_is_not_built_is_refused_by_name() {
    let mut m = Model::new();
    let (s, edge, top, front) = block(&mut m);
    let e = ogeom_fillet::fillet_edge_distances;
    let refused = |r: ogeom_core::OgeomResult<ogeom_algo::Built>, words: &str| {
        let message = r.expect_err("refused").to_string();
        assert!(message.contains(words), "{message}");
    };
    let bottom = face_at(&m, &s, Point::new(10.0, 10.0, 0.0));
    refused(
        e(
            &mut m,
            &s,
            &edge,
            &bottom,
            2.0,
            1.0,
            RoundProfile::Elliptical,
            T,
        ),
        "not one of the edge's two faces",
    );
    for rho in [0.0, 1.0, -0.5, f64::NAN] {
        refused(
            e(
                &mut m,
                &s,
                &edge,
                &top,
                2.0,
                1.0,
                RoundProfile::Conic(rho),
                T,
            ),
            "strictly between 0 and 1",
        );
    }
    refused(
        e(
            &mut m,
            &s,
            &edge,
            &top,
            0.0,
            1.0,
            RoundProfile::Elliptical,
            T,
        ),
        "rounds nothing",
    );
    refused(
        e(
            &mut m,
            &s,
            &edge,
            &front,
            12.0,
            1.0,
            RoundProfile::Elliptical,
            T,
        ),
        "past that face's side",
    );

    // A circular rim.
    let drum = ogeom_algo::make_cylinder(&mut m, Frame::WORLD, 5.0, 4.0, T)
        .unwrap()
        .shape;
    let rim = explore_unique(&m, &drum, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (c, _) = ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap();
            matches!(c, ogeom_geom::Curve::Circle(_))
                && ogeom_algo::edge_vertices(&m, e)
                    .unwrap()
                    .is_some_and(|(a, _)| (point_of(&m, &a).z - 4.0).abs() < 1e-9)
        })
        .expect("the top rim");
    let cap = face_at(&m, &drum, Point::new(0.0, 0.0, 4.0));
    refused(
        e(
            &mut m,
            &drum,
            &rim,
            &cap,
            1.0,
            0.5,
            RoundProfile::Elliptical,
            T,
        ),
        "straight edge between two planes",
    );
}

#[test]
fn blend_edges_takes_an_asymmetric_round_beside_others_and_refuses_one_at_their_corner() {
    let profile = RoundProfile::Conic(0.4);
    let asked = |face: &Shape| Blend::RoundDistances {
        face: face.clone(),
        on_face: 3.0,
        on_other: 1.5,
        profile,
    };

    // Alone: the single-edge call.
    let mut m = Model::new();
    let (s, edge, top, _) = block(&mut m);
    let alone = ogeom_fillet::blend_edges(
        &mut m,
        &s,
        &[(edge.clone(), asked(&top))],
        ChamferCorner::Point,
        T,
    )
    .unwrap();
    let single =
        ogeom_fillet::fillet_edge_distances(&mut m, &s, &edge, &top, 3.0, 1.5, profile, T).unwrap();
    assert_eq!(
        volume(&m, &alone.shape).to_bits(),
        volume(&m, &single.shape).to_bits()
    );

    // Beside a round on the bottom back edge, which it does not meet.
    let mut m = Model::new();
    let (s, edge, top, _) = block(&mut m);
    let back = edge_where(&m, &s, |p| (p.y - 20.0).abs() < 1e-9 && p.z.abs() < 1e-9);
    let both = ogeom_fillet::blend_edges(
        &mut m,
        &s,
        &[
            (back.clone(), Blend::Round(2.0)),
            (edge.clone(), asked(&top)),
        ],
        ChamferCorner::Point,
        T,
    )
    .unwrap();
    let exact = 4000.0
        - 20.0 * (4.0 - core::f64::consts::PI)
        - 20.0 * corner_area(3.0, 1.5, FRAC_PI_2, conic(0.4));
    let diagnosis = ogeom_algo::check(&m, &both.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let v = volume(&m, &both.shape);
    assert!((v - exact).abs() <= 1e-8, "volume {v} against {exact}");

    // Meeting a round at the box's corner: refused by name.
    let mut m = Model::new();
    let (s, edge, top, _) = block(&mut m);
    let side = edge_where(&m, &s, |p| p.x.abs() < 1e-9 && p.y.abs() < 1e-9);
    let message = ogeom_fillet::blend_edges(
        &mut m,
        &s,
        &[(side, Blend::Round(1.0)), (edge, asked(&top))],
        ChamferCorner::Point,
        T,
    )
    .expect_err("refused")
    .to_string();
    assert!(
        message.contains("corner where they meet is not closed"),
        "{message}"
    );
}

#[test]
fn equal_elliptical_distances_in_blend_edges_round_a_corner_with_other_rounds() {
    // Three edges at the box's corner: two rounds and an elliptical round of
    // equal distances, which is the circular round and closes with them.
    let corner_edges = |m: &Model, s: &Shape| {
        [
            edge_where(m, s, |p| p.y.abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9),
            edge_where(m, s, |p| p.x.abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9),
            edge_where(m, s, |p| p.x.abs() < 1e-9 && p.y.abs() < 1e-9),
        ]
    };
    let mut m = Model::new();
    let (s, _, top, _) = block(&mut m);
    let [a, b, c] = corner_edges(&m, &s);
    let asked = ogeom_fillet::blend_edges(
        &mut m,
        &s,
        &[
            (
                a.clone(),
                Blend::RoundDistances {
                    face: top,
                    on_face: 2.0,
                    on_other: 2.0,
                    profile: RoundProfile::Elliptical,
                },
            ),
            (b.clone(), Blend::Round(2.0)),
            (c.clone(), Blend::Round(2.0)),
        ],
        ChamferCorner::Point,
        T,
    )
    .unwrap();
    let rounds = ogeom_fillet::blend_edges(
        &mut m,
        &s,
        &[
            (a, Blend::Round(2.0)),
            (b, Blend::Round(2.0)),
            (c, Blend::Round(2.0)),
        ],
        ChamferCorner::Point,
        T,
    )
    .unwrap();
    assert_eq!(
        volume(&m, &asked.shape).to_bits(),
        volume(&m, &rounds.shape).to_bits()
    );
    assert_eq!(
        vertex_bits(&m, &asked.shape),
        vertex_bits(&m, &rounds.shape)
    );
}

/// Where the segment between a conic arc of weight `w` and its chord has
/// its centroid: the share of the way from the chord's middle to the
/// crease. Affine-invariant, so read on the unit circle, the parabola
/// (two fifths of the way to the arc's middle, itself halfway) or the
/// unit hyperbola.
fn segment_centroid(w: f64) -> f64 {
    if (w - 1.0).abs() < 1e-12 {
        0.2
    } else if w < 1.0 {
        let t = w.acos();
        let x = (2.0 / 3.0) * t.sin().powi(3) / (t - t.sin() * t.cos());
        (x - t.cos()) / (1.0 / t.cos() - t.cos())
    } else {
        let p = w.acosh();
        let x = (2.0 / 3.0) * p.sinh().powi(3) / (p.sinh() * p.cosh() - p);
        (x - p.cosh()) / (1.0 / p.cosh() - p.cosh())
    }
}

/// The corner a conic section takes off a square edge, legs `d1` and `d2`:
/// its area, and its centroid's distance from each face, from the face of
/// the second leg and from the face of the first.
fn square_corner(d1: f64, d2: f64, w: f64) -> (f64, [f64; 2]) {
    let triangle = 0.5 * d1 * d2;
    let area = corner_area(d1, d2, FRAC_PI_2, w);
    let segment = triangle - area;
    let shrink = 1.0 - segment_centroid(w);
    let centroid = [d1, d2].map(|d| (triangle * d / 3.0 - segment * shrink * d * 0.5) / area);
    (area, centroid)
}

/// A prism standing on its profile in the plane x = 0, extruded 2 back
/// along -x, the profile's points given as (y, z).
fn profile_prism(m: &mut Model, profile: &[(f64, f64)]) -> Shape {
    let points: Vec<Point> = profile
        .iter()
        .map(|&(y, z)| Point::new(0.0, y, z))
        .collect();
    let wire = ogeom_algo::make_polygon(m, &points, true, T).unwrap().shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::X),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(m, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    ogeom_algo::make_prism(m, &face, Vector::new(-2.0, 0.0, 0.0), T)
        .unwrap()
        .shape
}

const PROFILES: [RoundProfile; 2] = [RoundProfile::Elliptical, RoundProfile::Conic(0.65)];

fn weight_of(profile: RoundProfile) -> f64 {
    match profile {
        RoundProfile::Elliptical => elliptical(FRAC_PI_2),
        RoundProfile::Conic(rho) => conic(rho),
    }
}

/// A floor edge whose floor gives way to a ramp of slope k at its end, the
/// wall carrying on: the band runs on along the wall until the ramp has
/// cut its whole section away. Past the end, a stretch s on, the ramp
/// takes what of the section stands below k s; summed, the band takes its
/// section A along the edge and A z / k past it, z its centroid's height.
#[test]
fn an_asymmetric_round_runs_out_along_its_wall_onto_a_ramp() {
    let (length, height) = (3.0, 2.0);
    let (d_wall, d_floor) = (0.8, 0.5);
    for slope in [3.0_f64.sqrt(), 1.0] {
        for profile in PROFILES {
            let mut m = Model::new();
            let prism = profile_prism(
                &mut m,
                &[
                    (0.0, 0.0),
                    (length, 0.0),
                    (length + height / slope, height),
                    (0.0, height),
                ],
            );
            let before = volume(&m, &prism);
            let edge = edge_where(&m, &prism, |p| {
                p.x.abs() < 1e-9 && p.z.abs() < 1e-9 && p.y <= length + 1e-9
            });
            let wall = face_at(&m, &prism, Point::new(0.0, 1.0, 1.0));
            let floor = face_at(&m, &prism, Point::new(-1.0, 1.0, 0.0));
            let (area, centroid) = square_corner(d_wall, d_floor, weight_of(profile));
            let want = area * length + area * centroid[0] / slope;
            let mut volumes = Vec::new();
            for (face, a, b) in [(&wall, d_wall, d_floor), (&floor, d_floor, d_wall)] {
                let r = ogeom_fillet::fillet_edge_distances(
                    &mut m, &prism, &edge, face, a, b, profile, T,
                )
                .unwrap();
                let diagnosis = ogeom_algo::check(&m, &r.shape, T).unwrap();
                assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
                let removed = before - volume(&m, &r.shape);
                assert!(
                    (removed - want).abs() < 1e-8,
                    "slope {slope}, {profile:?}: took {removed} against {want}"
                );
                volumes.push(removed.to_bits());
            }
            assert_eq!(volumes[0], volumes[1], "naming the other face first");
        }
    }
}

/// A rib's upright edge running up under a slanted roof that turns over
/// the rib's corner: the band stops on the roof, cut off by its plane. At
/// a point of the section b in from the rib's side the roof stands
/// h - k b, so the band takes A (h - k b), b at the section's centroid.
#[test]
fn an_asymmetric_round_running_up_under_a_roof_stops_on_it() {
    let height = 2.0;
    let (d_end, d_side) = (0.5, 0.3);
    for slope in [0.5, 2.0] {
        for profile in PROFILES {
            let mut m = Model::new();
            let prism = profile_prism(
                &mut m,
                &[
                    (-3.0, 0.0),
                    (0.0, 0.0),
                    (0.0, height),
                    (2.0, 2.0f64.mul_add(slope, height)),
                    (2.0, 2.0f64.mul_add(slope, height + 1.0)),
                    (-3.0, 2.0f64.mul_add(slope, height + 1.0)),
                ],
            );
            let before = volume(&m, &prism);
            let edge = edge_where(&m, &prism, |p| {
                p.x.abs() < 1e-9 && p.y.abs() < 1e-9 && p.z <= height + 1e-9
            });
            // The prism's end at x = 0, and the rib's side at y = 0.
            let end = face_at(&m, &prism, Point::new(0.0, -1.0, 1.0));
            let side = face_at(&m, &prism, Point::new(-1.0, 0.0, 1.0));
            // The leg on the end runs along -y, away from the side.
            let (area, centroid) = square_corner(d_end, d_side, weight_of(profile));
            let want = area * slope.mul_add(-centroid[0], height);
            let mut volumes = Vec::new();
            for (face, a, b) in [(&end, d_end, d_side), (&side, d_side, d_end)] {
                let r = ogeom_fillet::fillet_edge_distances(
                    &mut m, &prism, &edge, face, a, b, profile, T,
                )
                .unwrap();
                let diagnosis = ogeom_algo::check(&m, &r.shape, T).unwrap();
                assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
                let removed = before - volume(&m, &r.shape);
                assert!(
                    (removed - want).abs() < 1e-8,
                    "slope {slope}, {profile:?}: took {removed} against {want}"
                );
                volumes.push(removed.to_bits());
            }
            assert_eq!(volumes[0], volumes[1], "naming the other face first");
        }
    }
}
