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

    // The straight edge where a half drum's flat side meets its cylinder:
    // its sections run off the drum's ends.
    let drum = ogeom_algo::make_cylinder(&mut m, Frame::WORLD, 5.0, 4.0, T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(-6.0, -6.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let half = ogeom_algo::make_box(&mut m, at, (6.0, 12.0, 6.0), T)
        .unwrap()
        .shape;
    let d = ogeom_bool::cut(&mut m, &drum, &half, T).unwrap().shape;
    let side = edge_where(&m, &d, |p| p.x.abs() < 1e-9 && (p.y - 5.0).abs() < 1e-9);
    let flat = face_at(&m, &d, Point::new(0.0, 0.0, 2.0));
    refused(
        e(
            &mut m,
            &d,
            &side,
            &flat,
            1.0,
            0.5,
            RoundProfile::Elliptical,
            T,
        ),
        "sections do not seat all along the edge",
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
    corner(d1, d2, FRAC_PI_2, w)
}

/// The corner a conic section takes off an edge whose legs `d1` and `d2`
/// meet at `alpha`: its area, and its centroid as so far along each unit
/// leg from the crease. The triangle's centroid stands a third along each
/// leg and the segment's on the line from the chord's middle to the
/// crease; both are affine, so they hold in the legs' own frame.
fn corner(d1: f64, d2: f64, alpha: f64, w: f64) -> (f64, [f64; 2]) {
    let triangle = 0.5 * d1 * d2 * alpha.sin();
    let area = corner_area(d1, d2, alpha, w);
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

/// The circle edge of `s` of radius `r` whose centre stands at height `z`
/// on the z axis through (`cx`, `cy`).
fn ring_edge(m: &Model, s: &Shape, (cx, cy, z): (f64, f64, f64), r: f64) -> Shape {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (c, _) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
            let ogeom_geom::Curve::Circle(c) = c else {
                return false;
            };
            let o = c.circle().centre();
            (c.circle().radius() - r).abs() < 1e-9
                && (o.x - cx).abs() < 1e-9
                && (o.y - cy).abs() < 1e-9
                && (o.z - z).abs() < 1e-9
        })
        .expect("the rim")
}

/// Checks a rim's round: valid, its volume `exact` within 1e-8, and its
/// band tangent to the cap and the wall along two rings, each of the
/// radius and height one of `rings` names about the axis through
/// (`cx`, `cy`); its volume.
fn rim_holds(
    m: &Model,
    s: &Shape,
    exact: f64,
    (cx, cy): (f64, f64),
    rings: [(f64, f64); 2],
) -> f64 {
    let diagnosis = ogeom_algo::check(m, s, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let v = volume(m, s);
    assert!(
        (v - exact).abs() <= 1e-8,
        "volume {v} against the exact {exact}"
    );
    // The rings each contact runs along, as each ring's arcs name them.
    let mut found = [false; 2];
    for f in explore_unique(m, s, ShapeType::Face).unwrap() {
        let data = m.node(&f).unwrap().data().as_face().unwrap();
        if !matches!(
            m.geometry().surface(data.surface),
            Some(ogeom_geom::SurfaceGeometry::Revolution(_))
        ) {
            continue;
        }
        for c in ogeom_fillet::analyse_blend(m, s, &f, 17, T).unwrap() {
            let (curve, _) = ogeom_algo::placed_edge_curve(m, &c.edge, T)
                .unwrap()
                .unwrap();
            let ogeom_geom::Curve::Circle(circle) = curve else {
                continue;
            };
            assert!(
                c.tangency_error <= 1e-7,
                "the band meets a face at {} rad",
                c.tangency_error
            );
            let o = circle.circle().centre();
            let r = circle.circle().radius();
            assert!(
                (o.x - cx).abs() < 1e-9 && (o.y - cy).abs() < 1e-9,
                "a contact ring off the axis"
            );
            let ring = rings
                .iter()
                .position(|(rr, z)| (r - rr).abs() < 1e-9 && (o.z - z).abs() < 1e-9);
            let Some(ring) = ring else {
                panic!(
                    "a contact ring of radius {r} at {} where no distance puts one",
                    o.z
                );
            };
            found[ring] = true;
        }
    }
    assert_eq!(
        found, [true; 2],
        "the band meets the cap and the wall along rings"
    );
    v
}

/// Rounds `edge` of `s` with `profile`, naming `cap` with `d_cap` and then
/// `wall` with `d_wall`, and checks each as [`rim_holds`] does, and that
/// both are the same solid to the bit.
#[allow(clippy::too_many_arguments, reason = "test fixture")]
fn rim_both_ways(
    m: &mut Model,
    s: &Shape,
    edge: &Shape,
    (cap, d_cap): (&Shape, f64),
    (wall, d_wall): (&Shape, f64),
    profile: RoundProfile,
    exact: f64,
    axis: (f64, f64),
    rings: [(f64, f64); 2],
) {
    let a =
        ogeom_fillet::fillet_edge_distances(m, s, edge, cap, d_cap, d_wall, profile, T).unwrap();
    let b =
        ogeom_fillet::fillet_edge_distances(m, s, edge, wall, d_wall, d_cap, profile, T).unwrap();
    let va = rim_holds(m, &a.shape, exact, axis, rings);
    let vb = rim_holds(m, &b.shape, exact, axis, rings);
    assert_eq!(va.to_bits(), vb.to_bits(), "naming the other face first");
    assert_eq!(vertex_bits(m, &a.shape), vertex_bits(m, &b.shape));
    assert!(a.history.is_deleted(edge));
}

/// A drum's top rim, convex, cylinder and cap square: the band takes the
/// corner of area A whose centroid stands c in from the wall, and Pappus
/// sums it round the axis, 2 pi (R - c) A.
#[test]
fn a_drum_s_rim_takes_its_conic_corner_round_the_axis() {
    let (radius, height, d_cap, d_wall) = (5.0, 4.0, 1.5, 0.8);
    for profile in PROFILES {
        let mut m = Model::new();
        let drum = ogeom_algo::make_cylinder(&mut m, Frame::WORLD, radius, height, T)
            .unwrap()
            .shape;
        let rim = ring_edge(&m, &drum, (0.0, 0.0, height), radius);
        let cap = face_at(&m, &drum, Point::new(0.0, 0.0, height));
        let wall = face_at(
            &m,
            &drum,
            Point::new(radius * FRAC_1_SQRT_2, radius * FRAC_1_SQRT_2, height * 0.5),
        );
        let (area, c) = square_corner(d_cap, d_wall, weight_of(profile));
        let exact = core::f64::consts::PI * radius * radius * height
            - 2.0 * core::f64::consts::PI * (radius - c[0]) * area;
        rim_both_ways(
            &mut m,
            &drum,
            &rim,
            (&cap, d_cap),
            (&wall, d_wall),
            profile,
            exact,
            (0.0, 0.0),
            [(radius - d_cap, height), (radius, height - d_wall)],
        );
    }
}

/// A bore's rim through a block: convex, the material outside the wall, so
/// the corner's centroid stands c out from it, 2 pi (R + c) A.
#[test]
fn a_bore_s_rim_takes_its_conic_corner_round_the_axis() {
    let (radius, d_cap, d_wall) = (4.0, 1.2, 2.0);
    for profile in PROFILES {
        let mut m = Model::new();
        let block = ogeom_algo::make_box(&mut m, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let at = Frame::new(Point::new(10.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap();
        let tool = ogeom_algo::make_cylinder(&mut m, at, radius, 12.0, T)
            .unwrap()
            .shape;
        let s = ogeom_bool::cut(&mut m, &block, &tool, T).unwrap().shape;
        let rim = ring_edge(&m, &s, (10.0, 10.0, 10.0), radius);
        let cap = face_at(&m, &s, Point::new(1.0, 1.0, 10.0));
        let wall = face_at(&m, &s, Point::new(10.0 - radius, 10.0, 5.0));
        let (area, c) = square_corner(d_cap, d_wall, weight_of(profile));
        let exact = 4000.0
            - core::f64::consts::PI * radius * radius * 10.0
            - 2.0 * core::f64::consts::PI * (radius + c[0]) * area;
        rim_both_ways(
            &mut m,
            &s,
            &rim,
            (&cap, d_cap),
            (&wall, d_wall),
            profile,
            exact,
            (10.0, 10.0),
            [(radius + d_cap, 10.0), (radius, 10.0 - d_wall)],
        );
    }
}

/// A boss's foot on a plate: concave, the band fills the corner outside
/// the boss, 2 pi (R + c) A more, and nothing past its contacts.
#[test]
fn a_boss_s_foot_fills_with_its_conic_corner_round_the_axis() {
    let (radius, d_cap, d_wall) = (3.0, 1.5, 1.0);
    for profile in PROFILES {
        let mut m = Model::new();
        let plate = ogeom_algo::make_box(&mut m, Frame::WORLD, (20.0, 20.0, 4.0), T)
            .unwrap()
            .shape;
        let at = Frame::new(Point::new(10.0, 10.0, 4.0), Direction::Z, Direction::X, T).unwrap();
        let boss = ogeom_algo::make_cylinder(&mut m, at, radius, 6.0, T)
            .unwrap()
            .shape;
        let s = ogeom_bool::fuse(&mut m, &plate, &boss, T).unwrap().shape;
        let rim = ring_edge(&m, &s, (10.0, 10.0, 4.0), radius);
        let cap = face_at(&m, &s, Point::new(1.0, 1.0, 4.0));
        let wall = face_at(&m, &s, Point::new(10.0 - radius, 10.0, 7.0));
        let (area, c) = square_corner(d_cap, d_wall, weight_of(profile));
        let exact = 1600.0
            + core::f64::consts::PI * radius * radius * 6.0
            + 2.0 * core::f64::consts::PI * (radius + c[0]) * area;
        rim_both_ways(
            &mut m,
            &s,
            &rim,
            (&cap, d_cap),
            (&wall, d_wall),
            profile,
            exact,
            (10.0, 10.0),
            [(radius + d_cap, 4.0), (radius, 4.0 + d_wall)],
        );
    }
}

/// A frustum's top rim: the cone leans in, so the cap and the wall meet at
/// alpha past a right angle and the corner's centroid stands c1 in along
/// the cap and c2 down the wall, which leans out by its own share.
#[test]
fn a_frustum_s_rim_takes_the_conic_of_its_own_dihedral() {
    let (base, top, height, d_cap, d_wall) = (6.0, 4.0, 4.0, 1.0, 1.5);
    let down = Vector::new(base - top, 0.0, -height);
    let down = down / down.magnitude();
    let alpha = (-down.x).acos();
    for profile in [RoundProfile::Elliptical, RoundProfile::Conic(0.35)] {
        let w = match profile {
            RoundProfile::Elliptical => elliptical(alpha),
            RoundProfile::Conic(rho) => conic(rho),
        };
        let mut m = Model::new();
        let cone = ogeom_algo::make_cone(&mut m, Frame::WORLD, base, top, height, T)
            .unwrap()
            .shape;
        let rim = ring_edge(&m, &cone, (0.0, 0.0, height), top);
        let cap = face_at(&m, &cone, Point::new(0.0, 0.0, height));
        let wall = face_at(
            &m,
            &cone,
            Point::new(-(base + top) * 0.5, 0.0, height * 0.5),
        );
        let (area, c) = corner(d_cap, d_wall, alpha, w);
        let reach = top - c[0] + c[1] * down.x;
        let frustum = core::f64::consts::PI * height / 3.0 * (base * base + base * top + top * top);
        let exact = frustum - 2.0 * core::f64::consts::PI * reach * area;
        rim_both_ways(
            &mut m,
            &cone,
            &rim,
            (&cap, d_cap),
            (&wall, d_wall),
            profile,
            exact,
            (0.0, 0.0),
            [
                (top - d_cap, height),
                (top + d_wall * down.x, height + d_wall * down.z),
            ],
        );
    }
}

#[test]
fn equal_elliptical_distances_on_a_rim_are_its_torus_round() {
    // A drum's square rim: the radius is the distance.
    let mut m = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut m, Frame::WORLD, 5.0, 4.0, T)
        .unwrap()
        .shape;
    let rim = ring_edge(&m, &drum, (0.0, 0.0, 4.0), 5.0);
    let cap = face_at(&m, &drum, Point::new(0.0, 0.0, 4.0));
    let asked = ogeom_fillet::fillet_edge_distances(
        &mut m,
        &drum,
        &rim,
        &cap,
        1.2,
        1.2,
        RoundProfile::Elliptical,
        T,
    )
    .unwrap();
    let round = ogeom_fillet::fillet_edge(&mut m, &drum, &rim, 1.2, T).unwrap();
    assert_eq!(
        volume(&m, &asked.shape).to_bits(),
        volume(&m, &round.shape).to_bits()
    );
    assert_eq!(vertex_bits(&m, &asked.shape), vertex_bits(&m, &round.shape));

    // A frustum's rim: the radius is the distance times tan(alpha / 2).
    let mut m = Model::new();
    let cone = ogeom_algo::make_cone(&mut m, Frame::WORLD, 6.0, 4.0, 4.0, T)
        .unwrap()
        .shape;
    let rim = ring_edge(&m, &cone, (0.0, 0.0, 4.0), 4.0);
    let wall = face_at(&m, &cone, Point::new(-5.0, 0.0, 2.0));
    let down = Vector::new(2.0, 0.0, -4.0);
    let down = down / down.magnitude();
    let alpha = (-down.x).acos();
    let radius = 1.0 * alpha.sin() / (1.0 + alpha.cos());
    let asked = ogeom_fillet::fillet_edge_distances(
        &mut m,
        &cone,
        &rim,
        &wall,
        1.0,
        1.0,
        RoundProfile::Elliptical,
        T,
    )
    .unwrap();
    // That radius is tan(alpha / 2) rounded, and the call's own reading of
    // the angle may round its last bit otherwise: the same round within
    // rounding.
    let round = ogeom_fillet::fillet_edge(&mut m, &cone, &rim, radius, T).unwrap();
    let (va, vr) = (volume(&m, &asked.shape), volume(&m, &round.shape));
    assert!(
        (va - vr).abs() <= 1e-12 * vr,
        "{va} against the torus round's {vr}"
    );
    let points = |s: &Shape| {
        let mut p: Vec<Point> = explore_unique(&m, s, ShapeType::Vertex)
            .unwrap()
            .iter()
            .map(|v| point_of(&m, v))
            .collect();
        p.sort_by(|a, b| a.z.total_cmp(&b.z).then(a.x.total_cmp(&b.x)));
        p
    };
    let (pa, pr) = (points(&asked.shape), points(&round.shape));
    assert_eq!(pa.len(), pr.len());
    for (a, r) in pa.iter().zip(&pr) {
        assert!(a.distance(*r) <= 1e-12, "{a:?} against {r:?}");
    }
}

#[test]
fn the_revolved_conic_of_the_circle_s_own_weight_is_the_torus_round() {
    let rho = FRAC_1_SQRT_2 / (1.0 + FRAC_1_SQRT_2);
    let mut m = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut m, Frame::WORLD, 5.0, 4.0, T)
        .unwrap()
        .shape;
    let rim = ring_edge(&m, &drum, (0.0, 0.0, 4.0), 5.0);
    let cap = face_at(&m, &drum, Point::new(0.0, 0.0, 4.0));
    let conic = ogeom_fillet::fillet_edge_distances(
        &mut m,
        &drum,
        &rim,
        &cap,
        1.2,
        1.2,
        RoundProfile::Conic(rho),
        T,
    )
    .unwrap();
    let round = ogeom_fillet::fillet_edge(&mut m, &drum, &rim, 1.2, T).unwrap();
    let (vc, vr) = (volume(&m, &conic.shape), volume(&m, &round.shape));
    assert!(
        (vc - vr).abs() <= 1e-9 * vr,
        "conic {vc} against the torus {vr}"
    );
}

#[test]
fn what_a_rim_does_not_build_is_refused_by_name() {
    let mut m = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut m, Frame::WORLD, 5.0, 4.0, T)
        .unwrap()
        .shape;
    let rim = ring_edge(&m, &drum, (0.0, 0.0, 4.0), 5.0);
    let cap = face_at(&m, &drum, Point::new(0.0, 0.0, 4.0));
    let bottom = face_at(&m, &drum, Point::new(0.0, 0.0, 0.0));
    let e = ogeom_fillet::fillet_edge_distances;
    let refused = |r: ogeom_core::OgeomResult<ogeom_algo::Built>, words: &str| {
        let message = r.expect_err("refused").to_string();
        assert!(message.contains(words), "{message}");
    };
    refused(
        e(
            &mut m,
            &drum,
            &rim,
            &bottom,
            1.0,
            0.5,
            RoundProfile::Elliptical,
            T,
        ),
        "not one of the edge's two faces",
    );
    refused(
        e(
            &mut m,
            &drum,
            &rim,
            &cap,
            6.0,
            0.5,
            RoundProfile::Elliptical,
            T,
        ),
        "reaches the axis",
    );
    refused(
        e(
            &mut m,
            &drum,
            &rim,
            &cap,
            1.0,
            4.5,
            RoundProfile::Elliptical,
            T,
        ),
        "past the side of a face",
    );
}

/// A drum's top rim asked as an asymmetric round beside a round of its
/// bottom rim, which it does not meet: each takes its own corner, Pappus
/// for the conic and the ball's corner, 2 pi (R - c) (1 - pi / 4) r^2,
/// for the round.
#[test]
fn blend_edges_takes_an_asymmetric_rim_beside_a_round_rim() {
    use core::f64::consts::PI;
    let (radius, height, r) = (5.0, 4.0, 1.0);
    let (d_cap, d_wall) = (1.5, 0.8);
    let profile = RoundProfile::Conic(0.45);
    let mut m = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut m, Frame::WORLD, radius, height, T)
        .unwrap()
        .shape;
    let top = ring_edge(&m, &drum, (0.0, 0.0, height), radius);
    let bottom = ring_edge(&m, &drum, (0.0, 0.0, 0.0), radius);
    let cap = face_at(&m, &drum, Point::new(0.0, 0.0, height));
    let built = ogeom_fillet::blend_edges(
        &mut m,
        &drum,
        &[
            (bottom, Blend::Round(r)),
            (
                top,
                Blend::RoundDistances {
                    face: cap,
                    on_face: d_cap,
                    on_other: d_wall,
                    profile,
                },
            ),
        ],
        ChamferCorner::Point,
        T,
    )
    .unwrap();
    let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let (area, c) = square_corner(d_cap, d_wall, weight_of(profile));
    let ball = (1.0 - PI / 4.0) * r * r;
    let ball_reach = r * (10.0 - 3.0 * PI) / (3.0 * (4.0 - PI));
    let exact = PI * radius * radius * height
        - 2.0 * PI * (radius - c[0]) * area
        - 2.0 * PI * (radius - ball_reach) * ball;
    let v = volume(&m, &built.shape);
    assert!((v - exact).abs() <= 1e-8, "volume {v} against {exact}");
}

/// The band's contact curves with a varying law: each edge the band
/// shares with a face whose two ends stand `length` apart along `along`,
/// with the worst angle between the band's and the face's normals on it.
fn rails(m: &Model, s: &Shape, along: Vector, length: f64) -> Vec<(Shape, f64)> {
    let mut out = Vec::new();
    for f in explore_unique(m, s, ShapeType::Face).unwrap() {
        let data = m.node(&f).unwrap().data().as_face().unwrap();
        if !matches!(
            m.geometry().surface(data.surface),
            Some(ogeom_geom::SurfaceGeometry::BSpline(_))
        ) {
            continue;
        }
        for c in ogeom_fillet::analyse_blend(m, s, &f, 33, T).unwrap() {
            let (a, b) = ogeom_algo::edge_vertices(m, &c.edge).unwrap().unwrap();
            let run = (point_of(m, &b) - point_of(m, &a)).dot(along).abs();
            if (run - length).abs() <= 1e-9 {
                out.push((c.edge, c.tangency_error));
            }
        }
    }
    out
}

/// Where a contact curve crosses the plane `x = at`: its curve walked by
/// bisection, the rail running monotonically along x.
fn rail_at(m: &Model, rail: &Shape, at: f64) -> Point {
    use ogeom_geom::Curve3d as _;
    let (curve, (mut lo, mut hi)) = ogeom_algo::placed_edge_curve(m, rail, T).unwrap().unwrap();
    let x = |t: f64| curve.point_at(t, T).unwrap().x;
    if x(lo) > x(hi) {
        core::mem::swap(&mut lo, &mut hi);
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if x(mid) < at {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    curve.point_at(0.5 * (lo + hi), T).unwrap()
}

/// Linear laws on a box edge: at a fraction t the section is the conic
/// corner of legs d1(t), d2(t), its area d1 d2 times the square corner's
/// share, so the band takes that share of the length times the integral
/// of d1 d2: for d1 = a + b t and d2 = c + e t, a c + (a e + b c) / 2 +
/// b e / 3.
#[test]
fn distances_varying_linearly_take_the_integral_of_their_corners() {
    use ogeom_fillet::RadiusLaw;
    let (a, b, c, e) = (1.0, 2.0, 2.0, -1.0);
    let top_law = RadiusLaw::Linear {
        start: a,
        end: a + b,
    };
    let front_law = RadiusLaw::Linear {
        start: c,
        end: c + e,
    };
    for profile in PROFILES {
        let mut m = Model::new();
        let (s, edge, top, front) = block(&mut m);
        let x = ogeom_fillet::fillet_edge_distances_variable(
            &mut m, &s, &edge, &top, &top_law, &front_law, profile, T,
        )
        .unwrap();
        let y = ogeom_fillet::fillet_edge_distances_variable(
            &mut m, &s, &edge, &front, &front_law, &top_law, profile, T,
        )
        .unwrap();
        let unit = corner_area(1.0, 1.0, FRAC_PI_2, weight_of(profile));
        let integral = a * c + (a * e + b * c) / 2.0 + b * e / 3.0;
        let exact = 4000.0 - 20.0 * unit * integral;
        let mut volumes = Vec::new();
        for r in [&x, &y] {
            let diagnosis = ogeom_algo::check(&m, &r.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
            let v = volume(&m, &r.shape);
            assert!(
                (v - exact).abs() <= 1e-8,
                "{profile:?}: volume {v} against {exact}"
            );
            let found = rails(&m, &r.shape, X, 20.0);
            assert_eq!(
                found.len(),
                2,
                "the band meets the top and the front along curves"
            );
            for (_, tangency) in &found {
                assert!(*tangency <= 1e-7, "the band meets a face at {tangency} rad");
            }
            // Each contact stands where its law puts it, a third of the way
            // along.
            let at = |t: f64| 20.0 * t;
            let mut hit = [false; 2];
            for (rail, _) in &found {
                let p = rail_at(&m, rail, at(1.0 / 3.0));
                if (p.z - 10.0).abs() <= 1e-9 {
                    assert!((p.y - (a + b / 3.0)).abs() <= 1e-9, "top contact at {p:?}");
                    hit[0] = true;
                } else {
                    assert!(
                        p.y.abs() <= 1e-9 && (10.0 - p.z - (c + e / 3.0)).abs() <= 1e-9,
                        "front contact at {p:?}"
                    );
                    hit[1] = true;
                }
            }
            assert_eq!(hit, [true; 2]);
            volumes.push(v.to_bits());
        }
        assert_eq!(volumes[0], volumes[1], "naming the other face first");
        assert_eq!(vertex_bits(&m, &x.shape), vertex_bits(&m, &y.shape));
    }
}

/// The same laws on the L-bracket's concave edge fill that integral.
#[test]
fn distances_varying_linearly_fill_a_concave_edge_with_the_integral() {
    use ogeom_fillet::RadiusLaw;
    let floor_law = RadiusLaw::Linear {
        start: 1.0,
        end: 3.0,
    };
    let wall_law = RadiusLaw::Linear {
        start: 2.0,
        end: 1.0,
    };
    for profile in PROFILES {
        let mut m = Model::new();
        let (s, edge, floor, _) = l_bracket(&mut m);
        let r = ogeom_fillet::fillet_edge_distances_variable(
            &mut m, &s, &edge, &floor, &floor_law, &wall_law, profile, T,
        )
        .unwrap();
        let diagnosis = ogeom_algo::check(&m, &r.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        let unit = corner_area(1.0, 1.0, FRAC_PI_2, weight_of(profile));
        let integral = 2.0 + (-1.0 + 4.0) / 2.0 - 2.0 / 3.0;
        let exact = 8000.0 - 12.0 * 12.0 * 20.0 + 20.0 * unit * integral;
        let v = volume(&m, &r.shape);
        assert!(
            (v - exact).abs() <= 1e-8,
            "{profile:?}: volume {v} against {exact}"
        );
        let found = rails(&m, &r.shape, Y, 20.0);
        assert_eq!(found.len(), 2);
        for (_, tangency) in &found {
            assert!(*tangency <= 1e-7, "the band meets a face at {tangency} rad");
        }
        let bounds = ogeom_algo::shape_bounds(&m, &r.shape, T).unwrap();
        let before = ogeom_algo::shape_bounds(&m, &s, T).unwrap();
        assert!(
            (bounds.low().unwrap() - before.low().unwrap()).magnitude() <= 1e-7
                && (bounds.high().unwrap() - before.high().unwrap()).magnitude() <= 1e-7,
            "the fill overhangs the bracket"
        );
    }
}

/// A law of points passes through each of its distances: the contact on
/// the top stands at each point's distance at its fraction of the edge,
/// and the band takes no less than the least corner and no more than the
/// largest along the length.
#[test]
fn a_law_of_points_passes_its_contact_through_each_point() {
    use ogeom_fillet::RadiusLaw;
    let points = vec![(0.0, 1.0), (0.4, 2.5), (1.0, 1.5)];
    let top_law = RadiusLaw::Points(points.clone());
    let front_law = RadiusLaw::Linear {
        start: 1.2,
        end: 1.2,
    };
    let mut m = Model::new();
    let (s, edge, top, _) = block(&mut m);
    let r = ogeom_fillet::fillet_edge_distances_variable(
        &mut m,
        &s,
        &edge,
        &top,
        &top_law,
        &front_law,
        RoundProfile::Elliptical,
        T,
    )
    .unwrap();
    let diagnosis = ogeom_algo::check(&m, &r.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let found = rails(&m, &r.shape, X, 20.0);
    assert_eq!(found.len(), 2);
    let top_rail = found
        .iter()
        .find(|(rail, _)| (rail_at(&m, rail, 10.0).z - 10.0).abs() <= 1e-9)
        .expect("the top's contact");
    for (t, d) in points {
        let p = rail_at(&m, &top_rail.0, 20.0 * t);
        assert!(
            (p.y - d).abs() <= 1e-9,
            "at {t} the top's contact stands {} off, not {d}",
            p.y
        );
    }
    for (_, tangency) in &found {
        assert!(*tangency <= 1e-7, "the band meets a face at {tangency} rad");
    }
    let unit = corner_area(1.0, 1.0, FRAC_PI_2, elliptical(FRAC_PI_2));
    let removed = 4000.0 - volume(&m, &r.shape);
    assert!(removed >= 20.0 * unit * 1.0 * 1.2 && removed <= 20.0 * unit * 2.5 * 1.2);
}

#[test]
fn level_laws_are_the_constant_asymmetric_round_and_oblique_ends_are_refused() {
    use ogeom_fillet::RadiusLaw;
    let mut m = Model::new();
    let (s, edge, top, front) = block(&mut m);
    let level = ogeom_fillet::fillet_edge_distances_variable(
        &mut m,
        &s,
        &edge,
        &top,
        &RadiusLaw::Linear {
            start: 3.0,
            end: 3.0,
        },
        &RadiusLaw::Points(vec![(0.5, 1.5)]),
        RoundProfile::Conic(0.4),
        T,
    )
    .unwrap();
    let constant = ogeom_fillet::fillet_edge_distances(
        &mut m,
        &s,
        &edge,
        &top,
        3.0,
        1.5,
        RoundProfile::Conic(0.4),
        T,
    )
    .unwrap();
    assert_eq!(
        volume(&m, &level.shape).to_bits(),
        volume(&m, &constant.shape).to_bits()
    );

    let refused = |r: ogeom_core::OgeomResult<ogeom_algo::Built>, words: &str| {
        let message = r.expect_err("refused").to_string();
        assert!(message.contains(words), "{message}");
    };
    let up = RadiusLaw::Linear {
        start: 1.0,
        end: 2.0,
    };
    refused(
        ogeom_fillet::fillet_edge_distances_variable(
            &mut m,
            &s,
            &edge,
            &front,
            &RadiusLaw::Linear {
                start: 1.0,
                end: 12.0,
            },
            &up,
            RoundProfile::Elliptical,
            T,
        ),
        "past that face's side",
    );
    // A block whose end at x = 20 leans: its top front edge ends at a face
    // not square to it.
    let mut m = Model::new();
    let points = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(20.0, 0.0, 0.0),
        Point::new(25.0, 20.0, 0.0),
        Point::new(0.0, 20.0, 0.0),
    ];
    let wire = ogeom_algo::make_polygon(&mut m, &points, true, T)
        .unwrap()
        .shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(&mut m, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let leaning = ogeom_algo::make_prism(&mut m, &face, Vector::new(0.0, 0.0, 10.0), T)
        .unwrap()
        .shape;
    let edge = edge_where(&m, &leaning, |p| {
        p.y.abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9
    });
    let top = face_at(&m, &leaning, Point::new(10.0, 10.0, 10.0));
    refused(
        ogeom_fillet::fillet_edge_distances_variable(
            &mut m,
            &leaning,
            &edge,
            &top,
            &up,
            &up,
            RoundProfile::Elliptical,
            T,
        ),
        "capped square",
    );
}

/// Half a drum, its flat side through the axis: an arc of its top rim
/// ends square on the flat or flush where the next arc carries on, and the
/// band takes Pappus' corner over the arc's own turn.
#[test]
fn a_rim_arc_takes_its_conic_corner_over_its_own_turn() {
    use core::f64::consts::PI;
    let (radius, height, d_cap, d_wall) = (5.0, 4.0, 1.5, 0.8);
    for profile in PROFILES {
        let mut m = Model::new();
        let drum = ogeom_algo::make_cylinder(&mut m, Frame::WORLD, radius, height, T)
            .unwrap()
            .shape;
        let at = Frame::new(Point::new(-6.0, -6.0, -1.0), Direction::Z, Direction::X, T).unwrap();
        let half = ogeom_algo::make_box(&mut m, at, (6.0, 12.0, 6.0), T)
            .unwrap()
            .shape;
        let s = ogeom_bool::cut(&mut m, &drum, &half, T).unwrap().shape;
        let arc = ring_edge(&m, &s, (0.0, 0.0, height), radius);
        let cap = face_at(&m, &s, Point::new(2.0, 0.0, height));
        let wall = face_at(
            &m,
            &s,
            Point::new(radius * FRAC_1_SQRT_2, radius * FRAC_1_SQRT_2, 2.0),
        );
        // The arc's own turn: the drum's rim halves part the half turn.
        let (_, (t0, t1)) = ogeom_algo::placed_edge_curve(&m, &arc, T).unwrap().unwrap();
        let sweep = (t1 - t0).abs();
        assert!(
            sweep > 0.1 && sweep < PI + 1e-9,
            "the rim arc turns {sweep}"
        );
        let (area, c) = square_corner(d_cap, d_wall, weight_of(profile));
        let exact = PI * radius * radius * height / 2.0 - sweep * (radius - c[0]) * area;
        let a =
            ogeom_fillet::fillet_edge_distances(&mut m, &s, &arc, &cap, d_cap, d_wall, profile, T)
                .unwrap();
        let b =
            ogeom_fillet::fillet_edge_distances(&mut m, &s, &arc, &wall, d_wall, d_cap, profile, T)
                .unwrap();
        let va = rim_holds(
            &m,
            &a.shape,
            exact,
            (0.0, 0.0),
            [(radius - d_cap, height), (radius, height - d_wall)],
        );
        let vb = rim_holds(
            &m,
            &b.shape,
            exact,
            (0.0, 0.0),
            [(radius - d_cap, height), (radius, height - d_wall)],
        );
        assert_eq!(va.to_bits(), vb.to_bits(), "naming the other face first");
    }
}

/// A pad whose upright corners are rounded at R: its top loop is four
/// straight edges between planes and four quarter rims, tangent at every
/// vertex. Asked as one asymmetric round on each edge, naming the top,
/// the loop closes: along the straight edges the band takes A per unit
/// length, round each rim Pappus' quarter turn, so the whole takes
/// A (L + 2 pi (R - c)), L the straight edges' length.
#[test]
fn a_tangent_loop_of_asymmetric_rounds_closes_round_a_pad() {
    use core::f64::consts::PI;
    let (dx, dy, dz, corner) = (20.0, 12.0, 4.0, 3.0);
    let (d_top, d_side) = (1.0, 0.6);
    for profile in PROFILES {
        let mut m = Model::new();
        let block = ogeom_algo::make_box(&mut m, Frame::WORLD, (dx, dy, dz), T)
            .unwrap()
            .shape;
        let uprights: Vec<Shape> = explore_unique(&m, &block, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
                let (pa, pb) = (point_of(&m, &a), point_of(&m, &b));
                (pa.x - pb.x).abs() < 1e-9 && (pa.y - pb.y).abs() < 1e-9
            })
            .collect();
        assert_eq!(uprights.len(), 4);
        let pad = ogeom_fillet::fillet_edges(&mut m, &block, &uprights, corner, T)
            .unwrap()
            .shape;
        let before = volume(&m, &pad);
        let top = face_at(&m, &pad, Point::new(dx * 0.5, dy * 0.5, dz));
        let loop_edges: Vec<Shape> = explore_unique(&m, &top, ShapeType::Edge).unwrap();
        assert_eq!(
            loop_edges.len(),
            8,
            "the top loop of four sides and four arcs"
        );
        let request: Vec<(Shape, Blend)> = loop_edges
            .iter()
            .map(|e| {
                (
                    e.clone(),
                    Blend::RoundDistances {
                        face: top.clone(),
                        on_face: d_top,
                        on_other: d_side,
                        profile,
                    },
                )
            })
            .collect();
        let built =
            ogeom_fillet::blend_edges(&mut m, &pad, &request, ChamferCorner::Point, T).unwrap();
        let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "{profile:?}: {:?}",
            diagnosis.problems
        );
        let (area, c) = square_corner(d_top, d_side, weight_of(profile));
        let straight = 2.0 * ((dx - 2.0 * corner) + (dy - 2.0 * corner));
        let exact = before - area * (straight + 2.0 * PI * (corner - c[0]));
        let v = volume(&m, &built.shape);
        assert!(
            (v - exact).abs() <= 1e-8,
            "{profile:?}: volume {v} against {exact}"
        );
    }
}

/// A pocket whose upright corners are rounded at R: its floor loop is
/// four concave straight edges and four concave quarter rims, tangent at
/// every vertex. Asked as one asymmetric round on each, naming the floor,
/// the loop fills: A along the straight edges and Pappus' quarter turn
/// round each rim, the fill standing in toward the corner's axis.
#[test]
fn a_tangent_loop_of_asymmetric_rounds_fills_round_a_pocket_s_floor() {
    use core::f64::consts::PI;
    let (dx, dy, depth, corner) = (20.0, 12.0, 6.0, 3.0);
    let (d_floor, d_wall) = (1.0, 0.6);
    for profile in PROFILES {
        let mut m = Model::new();
        let block = ogeom_algo::make_box(&mut m, Frame::WORLD, (30.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let at = Frame::new(
            Point::new(5.0, 4.0, 10.0 - depth),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        let hole = ogeom_algo::make_box(&mut m, at, (dx, dy, depth + 1.0), T)
            .unwrap()
            .shape;
        let pocket = ogeom_bool::cut(&mut m, &block, &hole, T).unwrap().shape;
        let uprights: Vec<Shape> = explore_unique(&m, &pocket, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
                let (pa, pb) = (point_of(&m, &a), point_of(&m, &b));
                (pa.x - pb.x).abs() < 1e-9
                    && (pa.y - pb.y).abs() < 1e-9
                    && pa.x > 1.0
                    && pa.x < 29.0
                    && pa.y > 1.0
                    && pa.y < 19.0
            })
            .collect();
        assert_eq!(uprights.len(), 4, "the pocket's upright corners");
        let rounded = ogeom_fillet::fillet_edges(&mut m, &pocket, &uprights, corner, T)
            .unwrap()
            .shape;
        let before = volume(&m, &rounded);
        let floor = face_at(&m, &rounded, Point::new(15.0, 10.0, 10.0 - depth));
        let loop_edges: Vec<Shape> = explore_unique(&m, &floor, ShapeType::Edge).unwrap();
        assert_eq!(
            loop_edges.len(),
            8,
            "the floor loop of four sides and four arcs"
        );
        let request: Vec<(Shape, Blend)> = loop_edges
            .iter()
            .map(|e| {
                (
                    e.clone(),
                    Blend::RoundDistances {
                        face: floor.clone(),
                        on_face: d_floor,
                        on_other: d_wall,
                        profile,
                    },
                )
            })
            .collect();
        let built =
            ogeom_fillet::blend_edges(&mut m, &rounded, &request, ChamferCorner::Point, T).unwrap();
        let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "{profile:?}: {:?}",
            diagnosis.problems
        );
        let (area, c) = square_corner(d_floor, d_wall, weight_of(profile));
        let straight = 2.0 * ((dx - 2.0 * corner) + (dy - 2.0 * corner));
        let exact = before + area * (straight + 2.0 * PI * (corner - c[0]));
        let v = volume(&m, &built.shape);
        assert!(
            (v - exact).abs() <= 1e-8,
            "{profile:?}: volume {v} against {exact}"
        );
    }
}

#[test]
fn asymmetric_rounds_meeting_other_than_as_a_tangent_chain_are_refused() {
    let mut m = Model::new();
    let (s, edge, top, _) = block(&mut m);
    // The top's left edge meets the front edge square at the corner.
    let left = edge_where(&m, &s, |p| p.x.abs() < 1e-9 && (p.z - 10.0).abs() < 1e-9);
    let asked = |on_face: f64| Blend::RoundDistances {
        face: top.clone(),
        on_face,
        on_other: 1.0,
        profile: RoundProfile::Conic(0.4),
    };
    let message = ogeom_fillet::blend_edges(
        &mut m,
        &s,
        &[(edge.clone(), asked(2.0)), (left.clone(), asked(2.0))],
        ChamferCorner::Point,
        T,
    )
    .expect_err("refused")
    .to_string();
    assert!(
        message.contains("other than as one tangent chain"),
        "{message}"
    );

    // A pad's tangent loop asked with two distances on the top.
    let mut m = Model::new();
    let block = ogeom_algo::make_box(&mut m, Frame::WORLD, (20.0, 12.0, 4.0), T)
        .unwrap()
        .shape;
    let uprights: Vec<Shape> = explore_unique(&m, &block, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
            let (pa, pb) = (point_of(&m, &a), point_of(&m, &b));
            (pa.x - pb.x).abs() < 1e-9 && (pa.y - pb.y).abs() < 1e-9
        })
        .collect();
    let pad = ogeom_fillet::fillet_edges(&mut m, &block, &uprights, 3.0, T)
        .unwrap()
        .shape;
    let top = face_at(&m, &pad, Point::new(10.0, 6.0, 4.0));
    let request: Vec<(Shape, Blend)> = explore_unique(&m, &top, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .enumerate()
        .map(|(k, e)| {
            (
                e,
                Blend::RoundDistances {
                    face: top.clone(),
                    on_face: if k == 0 { 1.2 } else { 1.0 },
                    on_other: 0.6,
                    profile: RoundProfile::Elliptical,
                },
            )
        })
        .collect();
    let message = ogeom_fillet::blend_edges(&mut m, &pad, &request, ChamferCorner::Point, T)
        .expect_err("refused")
        .to_string();
    assert!(
        message.contains("other than as one tangent chain"),
        "{message}"
    );
}

/// The S-shaped crease's profile in `y = 0`: from the origin along `x`,
/// down by 1 and on along `x` to `x = 10`, its tangent along `x` at both
/// ends.
fn spline() -> ogeom_geom::Curve {
    let at = |x: f64, z: f64| Point::new(x, 0.0, z);
    ogeom_geom::BSplineCurve::new(
        ogeom_math::KnotVector::clamped_uniform(3, 4).unwrap(),
        vec![at(0.0, 0.0), at(3.0, 0.0), at(7.0, -1.0), at(10.0, -1.0)],
        T,
    )
    .unwrap()
    .into()
}

/// A slab 2 thick along `y`: its profile in `y = 0` under the spline,
/// down to `z = -5`, between end walls square to the spline's ends; the
/// slab, the crease along the spline at `y = 2`, the top face there and
/// the spline's swept wall.
fn slab(m: &mut Model) -> (Shape, Shape, Shape, Shape) {
    use ogeom_geom::LineCurve;
    let at = |x: f64, z: f64| Point::new(x, 0.0, z);
    let corners: Vec<Shape> = [at(0.0, 0.0), at(10.0, -1.0), at(10.0, -5.0), at(0.0, -5.0)]
        .into_iter()
        .map(|p| ogeom_algo::make_vertex(m, p).shape)
        .collect();
    let line = |m: &mut Model, from: (&Shape, Point), to: (&Shape, Point)| {
        let curve: ogeom_geom::Curve = LineCurve::new(ogeom_math::Axis::new(
            from.1,
            Direction::new(to.1 - from.1, T).unwrap(),
        ))
        .into();
        ogeom_algo::make_edge_between(m, curve, (0.0, from.1.distance(to.1)), from.0, to.0, T)
            .unwrap()
            .shape
    };
    let edges = vec![
        ogeom_algo::make_edge_between(m, spline(), (0.0, 1.0), &corners[0], &corners[1], T)
            .unwrap()
            .shape,
        line(
            m,
            (&corners[1], at(10.0, -1.0)),
            (&corners[2], at(10.0, -5.0)),
        ),
        line(
            m,
            (&corners[2], at(10.0, -5.0)),
            (&corners[3], at(0.0, -5.0)),
        ),
        line(m, (&corners[3], at(0.0, -5.0)), (&corners[0], at(0.0, 0.0))),
    ];
    let wire = ogeom_algo::make_wire(m, &edges, T).unwrap().shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Y),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(m, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let solid = ogeom_algo::make_prism(m, &face, Vector::new(0.0, 2.0, 0.0), T)
        .unwrap()
        .shape;
    let crease = explore_unique(m, &solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let b = ogeom_algo::shape_bounds(m, e, T).unwrap();
            let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
            (lo.y - 2.0).abs() < 1e-6
                && (hi.y - 2.0).abs() < 1e-6
                && lo.x.abs() < 1e-6
                && (hi.x - 10.0).abs() < 1e-6
                && hi.z - lo.z > 0.5
        })
        .expect("the crease along the spline's top");
    let beside: Vec<Shape> = explore_unique(m, &solid, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            explore_unique(m, f, ShapeType::Edge)
                .unwrap()
                .iter()
                .any(|x| x.is_same(&crease))
        })
        .collect();
    let is_plane = |f: &Shape| {
        matches!(
            m.node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| m.geometry().surface(d.surface)),
            Some(ogeom_geom::SurfaceGeometry::Plane(_))
        )
    };
    let top = beside
        .iter()
        .find(|f| is_plane(f))
        .expect("the top")
        .clone();
    let wall = beside
        .iter()
        .find(|f| !is_plane(f))
        .expect("the wall")
        .clone();
    (solid, crease, top, wall)
}

/// The spline's length, and its tangent's turn between its ends, by
/// Gauss-Legendre over 256 pieces.
fn length_and_turn() -> (f64, f64) {
    use ogeom_geom::Curve3d as _;
    let curve = spline();
    let nodes = [
        (-0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
        (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
    ];
    let mut length = 0.0;
    for k in 0..256 {
        let (a, b) = (f64::from(k) / 256.0, f64::from(k + 1) / 256.0);
        for (x, w) in nodes {
            let t = (b - a).mul_add(x, a + b) * 0.5;
            length += w * (b - a) * 0.5 * curve.d1_at(t, T).unwrap().magnitude();
        }
    }
    let heading = |t: f64| {
        let d = curve.d1_at(t, T).unwrap();
        d.z.atan2(d.x)
    };
    (length, heading(1.0) - heading(0.0))
}

/// The slab's spline crease, its top and swept wall square all along it:
/// the band's section is the square conic corner at every station, and the
/// crease turns through nothing between its ends, so the band takes the
/// corner's area times the crease's length, whichever face is named.
#[test]
fn an_open_spline_crease_takes_its_conic_corner_along_its_length() {
    let (length, turn) = length_and_turn();
    assert!(turn.abs() < 1e-12, "{turn}");
    let (d_top, d_wall) = (0.4, 0.7);
    for profile in PROFILES {
        let mut m = Model::new();
        let (solid, crease, top, wall) = slab(&mut m);
        let before = volume(&m, &solid);
        let (area, _) = square_corner(d_top, d_wall, weight_of(profile));
        let expected = before - area * length;
        let mut volumes = Vec::new();
        for (face, a, b) in [(&top, d_top, d_wall), (&wall, d_wall, d_top)] {
            let built = ogeom_fillet::fillet_edge_distances(
                &mut m, &solid, &crease, face, a, b, profile, T,
            )
            .unwrap();
            let diagnosis = ogeom_algo::check(&m, &built.shape, T).unwrap();
            assert!(
                diagnosis.is_valid(),
                "{profile:?}: {:?}",
                diagnosis.problems
            );
            assert!(built.history.is_deleted(&crease));
            let v = volume(&m, &built.shape);
            assert!(
                (v - expected).abs() <= 1e-5 * expected,
                "{profile:?}: {v} against {expected}"
            );
            volumes.push(v);
        }
        assert!(
            (volumes[0] - volumes[1]).abs() <= 1e-9 * volumes[0],
            "naming the other face first"
        );
    }
}

/// A drum of radius 10 and height 30 with a bore of radius 3 across it
/// along y at height 15: the bore's rim on the drum is a closed crease
/// only the march speaks.
fn holed_drum(m: &mut Model) -> (Shape, Shape, Shape, Shape) {
    let drum = ogeom_algo::make_cylinder(m, Frame::WORLD, 10.0, 30.0, T)
        .unwrap()
        .shape;
    let f = Frame::new(Point::new(0.0, -20.0, 15.0), Direction::Y, Direction::X, T).unwrap();
    let bore = ogeom_algo::make_cylinder(m, f, 3.0, 40.0, T).unwrap().shape;
    let s = ogeom_bool::cut(m, &drum, &bore, T).unwrap().shape;
    let q = Point::new(0.0, 9.5, 18.0);
    let rim = explore_unique(m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| {
            let d = |e: &Shape| {
                let (c, (lo, hi)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
                (0..=64)
                    .map(|i| {
                        use ogeom_geom::Curve3d as _;
                        c.point_at(lo + (hi - lo) * f64::from(i) / 64.0, T)
                            .unwrap()
                            .distance(q)
                    })
                    .fold(f64::INFINITY, f64::min)
            };
            d(a).total_cmp(&d(b))
        })
        .unwrap();
    let beside: Vec<Shape> = explore_unique(m, &s, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            explore_unique(m, f, ShapeType::Edge)
                .unwrap()
                .iter()
                .any(|x| x.is_same(&rim))
        })
        .collect();
    let is_drum = |f: &Shape| {
        matches!(
            m.node(f).and_then(|n| n.data().as_face()).and_then(|d| m.geometry().surface(d.surface)),
            Some(ogeom_geom::SurfaceGeometry::Cylinder(c)) if (c.cylinder().radius() - 10.0).abs() < 1e-9
        )
    };
    let drum_face = beside
        .iter()
        .find(|f| is_drum(f))
        .expect("the drum")
        .clone();
    let bore_face = beside
        .iter()
        .find(|f| !is_drum(f))
        .expect("the bore")
        .clone();
    (s, rim, drum_face, bore_face)
}

/// The bore's rim on a drum, a closed crease the march speaks, rounded
/// with a conic section: valid, the band's borders on the drum and the
/// bore within the band's fit, tangent to both within it, the same solid
/// whichever face is named, and the material it takes more than nothing
/// and less than the bevel of the same legs takes.
#[test]
fn a_bore_s_rim_on_a_drum_takes_a_marched_conic_band() {
    let (d_drum, d_bore) = (0.8, 0.4);
    for profile in PROFILES {
        let mut m = Model::new();
        let (s, rim, drum, bore) = holed_drum(&mut m);
        let before = volume(&m, &s);
        let a = ogeom_fillet::fillet_edge_distances(
            &mut m, &s, &rim, &drum, d_drum, d_bore, profile, T,
        )
        .unwrap();
        let b = ogeom_fillet::fillet_edge_distances(
            &mut m, &s, &rim, &bore, d_bore, d_drum, profile, T,
        )
        .unwrap();
        let diagnosis = ogeom_algo::check(&m, &a.shape, T).unwrap();
        assert!(
            diagnosis.is_valid(),
            "{profile:?}: {:?}",
            diagnosis.problems
        );
        assert!(a.history.is_deleted(&rim));
        let (va, vb) = (volume(&m, &a.shape), volume(&m, &b.shape));
        assert_eq!(va.to_bits(), vb.to_bits(), "naming the other face first");
        let spec = ogeom_fillet::Chamfer::Distances {
            face: drum.clone(),
            on_face: d_drum,
            on_other: d_bore,
        };
        let bevel =
            ogeom_fillet::chamfer_edges_with(&mut m, &s, &[(rim.clone(), spec)], T).unwrap();
        let vc = volume(&m, &bevel.shape);
        assert!(
            vc < va && va < before,
            "{profile:?}: bevel {vc}, round {va}, solid {before}"
        );
        let band = explore_unique(&m, &a.shape, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                let d = m.node(f).unwrap().data().as_face().unwrap();
                matches!(
                    m.geometry().surface(d.surface),
                    Some(ogeom_geom::SurfaceGeometry::BSpline(_))
                )
            })
            .expect("the band");
        let mut worst = 0.0_f64;
        for c in ogeom_fillet::analyse_blend(&m, &a.shape, &band, 33, T).unwrap() {
            assert!(c.gap <= 2e-4, "a rail stands {} off its face", c.gap);
            worst = worst.max(c.tangency_error);
        }
        assert!(
            worst <= 1e-3,
            "{profile:?}: the band meets a face at {worst} rad"
        );
    }
}
